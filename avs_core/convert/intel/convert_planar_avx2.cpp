// Avisynth v2.5.  Copyright 2002 Ben Rudiak-Gould et al.
// http://avisynth.nl

// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA, or visit
// http://www.gnu.org/copyleft/gpl.html .
//
// Linking Avisynth statically or dynamically with other modules is making a
// combined work based on Avisynth.  Thus, the terms and conditions of the GNU
// General Public License cover the whole combination.
//
// As a special exception, the copyright holders of Avisynth give you
// permission to link Avisynth with independent modules that communicate with
// Avisynth solely through the interfaces defined in avisynth.h, regardless of the license
// terms of these independent modules, and to copy and distribute the
// resulting combined work under terms of your choice, provided that
// every copy of the combined work is accompanied by a complete copy of
// the source code of Avisynth (the version of Avisynth used to produce the
// combined work), being distributed under the terms of the GNU General
// Public License plus this exception.  An independent module is a module
// which is not derived from or based on Avisynth, such as 3rd-party filters,
// import and export plugins, or graphical user interfaces.

// ConvertPlanar (c) 2005 by Klaus Post


#include <avs/alignment.h>
#ifdef _MSC_VER
    #include <intrin.h>
#else
    #include <x86intrin.h>
#endif
#include <immintrin.h>

#include "convert_planar_avx2.h"

#ifndef _mm256_set_m128i
#define _mm256_set_m128i(v0, v1) _mm256_insertf128_si256(_mm256_castsi128_si256(v1), (v0), 1)
#endif

DISABLE_WARNING_PUSH
DISABLE_WARNING_UNREFERENCED_LOCAL_VARIABLE

// packed rgb helper
static AVS_FORCEINLINE __m256i convert_yuv_to_rgb_avx2_core(const __m256i& px0189, const __m256i& px23AB, const __m256i& px45CD, const __m256i& px67EF, const __m256i& zero, const __m256i& matrix, const __m256i& round_mask_plus_rgb_offset) {
  //int b = (((int)m[0] * Y + (int)m[1] * U + (int)m[ 2] * V + 4096 + rgb_offset)>>13);

  //px01 - xx xx 00 V1 00 U1 00 Y1 xx xx 00 V0 00 U0 00 Y0

  auto low_lo = _mm256_madd_epi16(px0189, matrix); //xx*0 + v1*m2 | u1*m1 + y1*m0 | xx*0 + v0*m2 | u0*m1 + y0*m0
  auto low_hi = _mm256_madd_epi16(px23AB, matrix); //xx*0 + v3*m2 | u3*m1 + y3*m0 | xx*0 + v2*m2 | u2*m1 + y2*m0
  auto high_lo = _mm256_madd_epi16(px45CD, matrix);
  auto high_hi = _mm256_madd_epi16(px67EF, matrix);

  auto low_v = _mm256_castps_si256(_mm256_shuffle_ps(_mm256_castsi256_ps(low_lo), _mm256_castsi256_ps(low_hi), _MM_SHUFFLE(3, 1, 3, 1))); // v3*m2 | v2*m2 | v1*m2 | v0*m2
  auto high_v = _mm256_castps_si256(_mm256_shuffle_ps(_mm256_castsi256_ps(high_lo), _mm256_castsi256_ps(high_hi), _MM_SHUFFLE(3, 1, 3, 1)));

  auto low_yu = _mm256_castps_si256(_mm256_shuffle_ps(_mm256_castsi256_ps(low_lo), _mm256_castsi256_ps(low_hi), _MM_SHUFFLE(2, 0, 2, 0))); // u3*m1 + y3*m0 | u2*m1 + y2*m0 | u1*m1 + y1*m0 | u0*m1 + y0*m0
  auto high_yu = _mm256_castps_si256(_mm256_shuffle_ps(_mm256_castsi256_ps(high_lo), _mm256_castsi256_ps(high_hi), _MM_SHUFFLE(2, 0, 2, 0)));

  auto t_lo = _mm256_add_epi32(low_v, low_yu); // v3*m2 + u3*m1 + y3*m0...
  auto t_hi = _mm256_add_epi32(high_v, high_yu);

  // v3*m2 + u3*m1 + y3*m0 + 4096 + rgb_offset
  t_lo = _mm256_add_epi32(t_lo, round_mask_plus_rgb_offset);
  t_hi = _mm256_add_epi32(t_hi, round_mask_plus_rgb_offset);

  t_lo = _mm256_srai_epi32(t_lo, 13); // (v3*m2 + u3*m1 + y3*m0 + 4096) >> 13...
  t_hi = _mm256_srai_epi32(t_hi, 13);

  auto result = _mm256_packs_epi32(t_lo, t_hi);
  result = _mm256_packus_epi16(result, zero); //00 00 00 00 00 00 00 00 b15 b14 b14 b12 b11 b10 b9 b8 00 00 00 00 00 00 00 00 b7 b6 b5 b4 b3 b2 b1 b0
  return result;
}

template<int rgb_pixel_step, bool hasAlpha>
void convert_yv24_to_rgb_avx2(BYTE* dstp, const BYTE* srcY, const BYTE* srcU, const BYTE* srcV, const BYTE* srcA, size_t dst_pitch, size_t src_pitch_y, size_t src_pitch_uv, size_t src_pitch_a, size_t width, size_t height, const ConversionMatrix& matrix)
{
  dstp += dst_pitch * (height - 1);  // We start at last line

  size_t mod16_width = rgb_pixel_step == 3 ? width / 16 * 16 : width;
  // for rgb32 target we can process pixels beyond width, but we have 64bit alignment at target 16 pixels to 16*4=64 rgb pixels
  // if alignment would only be 32 bytes, we'd do the last cycle to process only 8 source pixels

  auto matrix_b = _mm256_set_epi16(0, matrix.v_b, matrix.u_b, matrix.y_b, 0, matrix.v_b, matrix.u_b, matrix.y_b, 0, matrix.v_b, matrix.u_b, matrix.y_b, 0, matrix.v_b, matrix.u_b, matrix.y_b);
  auto matrix_g = _mm256_set_epi16(0, matrix.v_g, matrix.u_g, matrix.y_g, 0, matrix.v_g, matrix.u_g, matrix.y_g, 0, matrix.v_g, matrix.u_g, matrix.y_g, 0, matrix.v_g, matrix.u_g, matrix.y_g);
  auto matrix_r = _mm256_set_epi16(0, matrix.v_r, matrix.u_r, matrix.y_r, 0, matrix.v_r, matrix.u_r, matrix.y_r, 0, matrix.v_r, matrix.u_r, matrix.y_r, 0, matrix.v_r, matrix.u_r, matrix.y_r);

  auto zero128 = _mm_setzero_si128();
  auto zero = _mm256_setzero_si256();

  // .13 bit frac integer arithmetic
  int round_mask_plus_rgb_offset_i = (1 << 12) + (matrix.offset_rgb << 13);
  auto round_mask_plus_rgb_offset = _mm256_set1_epi32(round_mask_plus_rgb_offset_i);
  auto offset = _mm256_set_epi16(0, -128, -128, matrix.offset_y, 0, -128, -128, matrix.offset_y, 0, -128, -128, matrix.offset_y, 0, -128, -128, matrix.offset_y);

  // 16 YUV(A) pixels --> 4x16 RGB quads = 64 bytes. Avisynth's alignment is 64 fortunately.

  for (size_t y = 0; y < height; ++y) {
    for (size_t x = 0; x < mod16_width; x += 16) {
      __m128i src_y = _mm_load_si128(reinterpret_cast<const __m128i*>(srcY + x)); //Y15 .. Y7 Y6 Y5 Y4 Y3 Y2 Y1 Y0
      __m128i src_u = _mm_load_si128(reinterpret_cast<const __m128i*>(srcU + x)); //U15 .. U7 U6 U5 U4 U3 U2 U1 U0
      __m128i src_v = _mm_load_si128(reinterpret_cast<const __m128i*>(srcV + x)); //V15 .. V7 V6 V5 V4 V3 V2 V1 V0
      [[maybe_unused]] __m128i src_a;
      if constexpr(hasAlpha)
        src_a = _mm_load_si128(reinterpret_cast<const __m128i*>(srcA + x)); //A15 .. A7 A6 A5 A4 A3 A2 A1 A0

      __m128i t1_lo = _mm_unpacklo_epi8(src_y, src_u); //U7 Y7 U6 Y6 U5 Y5 U4 Y4 U3 Y3 U2 Y2 U1 Y1 U0 Y0
      __m128i t1_hi = _mm_unpackhi_epi8(src_y, src_u); //U15 Y15 U14 Y14 U13 Y13 U12 Y12 U11 Y11 U10 Y10 U9 Y9 U8 Y8
      __m128i t2_lo = _mm_unpacklo_epi8(src_v, zero128);  //00 V7 00 V6 00 V5 00 V4 00 V3 00 V2 00 V1 00 V0
      __m128i t2_hi = _mm_unpackhi_epi8(src_v, zero128);  //00 V15 00 V14 00 V13 00 V12 00 V11 00 V10 00 V9 00 V8

      __m256i t1 = _mm256_set_m128i(t1_hi, t1_lo);
      __m256i t2 = _mm256_set_m128i(t2_hi, t2_lo);
      t1 = _mm256_permute4x64_epi64(t1, (0 << 0) | (2 << 2) | (1 << 4) | (3 << 6));
      t2 = _mm256_permute4x64_epi64(t2, (0 << 0) | (2 << 2) | (1 << 4) | (3 << 6));

      __m256i low = _mm256_unpacklo_epi16(t1, t2); //xx V11 U11 Y11 xx V10 U10 Y10 xx V9 U9 Y9 xx V8 U8 Y8    xx V3 U3 Y3 xx V2 U2 Y2 xx V1 U1 Y1 xx V0 U0 Y0
      __m256i high = _mm256_unpackhi_epi16(t1, t2); //xx V15 U15 Y15 xx V14 U14 Y14 xx V13 U13 Y13 xx V12 U12 Y12    xx V7 U7 Y7 xx V6 U6 Y6 xx V5 U5 Y5 xx V4 U4 Y4

      __m256i px0189 = _mm256_unpacklo_epi8(low, zero);  //xx xx 00 V1 00 U1 00 Y1 xx xx 00 V0 00 U0 00 Y0
      __m256i px23AB = _mm256_unpackhi_epi8(low, zero);  //xx xx 00 V3 00 U3 00 Y3 xx xx 00 V2 00 U2 00 Y2
      __m256i px45CD = _mm256_unpacklo_epi8(high, zero); //xx xx 00 V5 00 U5 00 Y5 xx xx 00 V4 00 U4 00 Y4
      __m256i px67EF = _mm256_unpackhi_epi8(high, zero); //xx xx 00 V7 00 U7 00 Y7 xx xx 00 V6 00 U6 00 Y6

      px0189 = _mm256_add_epi16(px0189, offset);
      px23AB = _mm256_add_epi16(px23AB, offset);
      px45CD = _mm256_add_epi16(px45CD, offset);
      px67EF = _mm256_add_epi16(px67EF, offset);

      __m256i result_b = convert_yuv_to_rgb_avx2_core(px0189, px23AB, px45CD, px67EF, zero, matrix_b, round_mask_plus_rgb_offset); // b15..0
      __m256i result_g = convert_yuv_to_rgb_avx2_core(px0189, px23AB, px45CD, px67EF, zero, matrix_g, round_mask_plus_rgb_offset); // g15..0
      __m256i result_r = convert_yuv_to_rgb_avx2_core(px0189, px23AB, px45CD, px67EF, zero, matrix_r, round_mask_plus_rgb_offset); // r15..0

      __m256i result_bg = _mm256_unpacklo_epi8(result_b, result_g); //g15 b15 g14 b14 g13 b13 g12 b12 g11 b11 g10 b10 g9 b9 g8 b8 | g7 b7 g6 b6 g5 b5 g4 b4 g3 b3 g2 b2 g1 b1 g0 b0
      __m256i alpha;
      if constexpr(hasAlpha) {
        __m128i a_lo = _mm_unpacklo_epi8(src_a, zero128);  //00 A7 00 A6 00 A5 00 A4 00 A3 00 A2 00 A1 00 A0
        __m128i a_hi = _mm_unpackhi_epi8(src_a, zero128);  //00 A15 00 A14 00 A13 00 A12 00 A11 00 A10 00 A9 00 A8
        alpha = _mm256_set_m128i(a_hi, a_lo); // a15 .. a0  at low of each m128i part
        alpha = _mm256_permute4x64_epi64(alpha, (0 << 0) | (2 << 2) | (1 << 4) | (3 << 6));
      }
      else
        alpha = _mm256_cmpeq_epi32(result_r, result_r); // FF FF FF FF ... default alpha transparent

      __m256i result_ra = _mm256_unpacklo_epi8(result_r, alpha);       //a7 r7 a6 r6 a5 r5 a4 r4 a3 r3 a2 r2 a1 r1 a0 r0

      __m256i result_lo = _mm256_unpacklo_epi16(result_bg, result_ra); // a11 r11 g11 b11 | a10 r10 g10 b10 | a9 r9 g9 b9 | a8 r8 g8 b8  |   a3 r3 g3 b3 | a2 r2 g2 b2 | a1 r1 g1 b1 | a0 r0 g0 b0
      __m256i result_hi = _mm256_unpackhi_epi16(result_bg, result_ra);
      // note: actual pixel indexes are different from the numbers in comments, they are kept to be follow more easily
      // Initial _mm256_permute4x64_epi64 ensures that the order here is properly argb7..argb0 and argb15..argb8

      if constexpr (rgb_pixel_step == 4) {
        //rgb32
        _mm256_store_si256(reinterpret_cast<__m256i*>(dstp + x * 4), result_lo);
        _mm256_store_si256(reinterpret_cast<__m256i*>(dstp + x * 4 + 32), result_hi);
      }
      else {
        // rgb24
        // 16*4 bytes to 16*3 bytes
        __m256i perm10 = _mm256_set_epi32(0, 0, 6, 5, 4, 2, 1, 0);
        __m256i shuffle_lo = _mm256_set_epi8(
            0, 0, 0, 0, 14, 13, 12, 10, 9, 8, 6, 5, 4, 2, 1, 0,
            0, 0, 0, 0, 14, 13, 12, 10, 9, 8, 6, 5, 4, 2, 1, 0
        );
        __m128i shuffle_hi_lo = _mm_set_epi8(9, 8, 6, 5,  4,  2,  1,  0, 0, 0, 0, 0, 0, 0, 0, 0);

        __m128i result_hi_lo = _mm256_extracti128_si256(result_hi, 0); // aBbBgBrB aAbAgArA a9b9g9r9 a8b8g8r8
        __m128i result_hi_hi = _mm256_extracti128_si256(result_hi, 1); // aFbFgFrF aEbEgErE aDbDgDrD aCbCgCrC
        
        __m256i result_lo_reorg = _mm256_shuffle_epi8(result_lo, shuffle_lo);
        // x  x  x  x  b7g7r7 b6g6r6 b5g5r5 b4g4r4 x  x  x  x  b3g3r3 b2g2r2 b1g1r1 b0g0r0
        result_lo_reorg = _mm256_permutevar8x32_epi32(result_lo_reorg, perm10);
        // x  x  x  x  x x x x b7g7r7 b6g6r6 b5g5r5 b4g4r4 b3g3r3 b2g2r2 b1g1r1 b0g0r0

        __m128i result_hi_lo_reorg = _mm_shuffle_epi8(result_hi_lo, shuffle_hi_lo);
        // gA rA b9 g9 r9 b8 g8 r8 x x x x x x x x

        __m256i dummy_y0 = _mm256_undefined_si256();
        auto result_hi_lo_reorg_2 = _mm256_inserti128_si256(dummy_y0, result_hi_lo_reorg, 1);
        //                                                                      
        // gA rA b9 g9|r9 b8 g8 r8|x x x  x|x x  x x|x  x x x |x x x  x|x x  x x|x  x x x   // result_hi_lo_reorg_2
        // x  x  x  x  x  x  x  x  b7g7r7 b6g6r6 b5g5r5 b4g4r4 b3g3r3 b2g2r2 b1g1r1 b0g0r0  // result_lo_reorg
        auto result_0_15 = _mm256_blend_epi32(result_lo_reorg, result_hi_lo_reorg_2, 0xC0); // 11000000
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(dstp + x * 3), result_0_15); // not necessarily 32 bytes aligned

        // => x  x  x  x  x  x  x  x  x  x  x  x  bB gB rB bA
        __m128i shuffle_hi_lo_2 = _mm_set_epi8((char)0x80, (char)0x80, (char)0x80, (char)0x80, (char)0x80, (char)0x80, (char)0x80, (char)0x80, (char)0x80, (char)0x80, (char)0x80, (char)0x80, 14, 13, 12, 10);
        __m128i xmm4 = _mm_shuffle_epi8(result_hi_lo, shuffle_hi_lo_2);

        // => bF gF rF bE gE rE bD gD rD bC gC rC x  x  x  x
        __m128i shuffle_hi_hi = _mm_set_epi8(14, 13, 12, 10, 9, 8, 6, 5, 4, 2, 1, 0, (char)0x80, (char)0x80, (char)0x80, (char)0x80);
        __m128i xmm5 = _mm_shuffle_epi8(result_hi_hi, shuffle_hi_hi);

        // => bF gF rF bE gE rE bD gD rD bC gC rC bB gB rB bA
        __m128i result_16_23 = _mm_or_si128(xmm4, xmm5);
        _mm_store_si128(reinterpret_cast<__m128i*>(dstp + x * 3 + 32), result_16_23);
#if 0
        // Intel compiler can cope with it 100% optimized. No store, just shuffles and blends as above.
        alignas(32) BYTE temp[64];
        _mm256_store_si256(reinterpret_cast<__m256i*>(temp), result_lo);
        _mm256_store_si256(reinterpret_cast<__m256i*>(temp + 32), result_hi);
        for (int i = 0; i < 16; ++i) {
          dstp[(x + i) * 3 + 0] = temp[i * 4 + 0];
          dstp[(x + i) * 3 + 1] = temp[i * 4 + 1];
          dstp[(x + i) * 3 + 2] = temp[i * 4 + 2];
        }
#endif
      }
    }

    if constexpr (rgb_pixel_step == 3) {
      // for rgb32 (pixel_step == 4) we processed full width and more, including padded bytes
      for (size_t x = mod16_width; x < width; ++x) {
        int Y = srcY[x] + matrix.offset_y;
        int U = srcU[x] - 128;
        int V = srcV[x] - 128;
        int b = (((int)matrix.y_b * Y + (int)matrix.u_b * U + (int)matrix.v_b * V + round_mask_plus_rgb_offset_i) >> 13);
        int g = (((int)matrix.y_g * Y + (int)matrix.u_g * U + (int)matrix.v_g * V + round_mask_plus_rgb_offset_i) >> 13);
        int r = (((int)matrix.y_r * Y + (int)matrix.u_r * U + (int)matrix.v_r * V + round_mask_plus_rgb_offset_i) >> 13);
        dstp[x * rgb_pixel_step + 0] = PixelClip(b);
        dstp[x * rgb_pixel_step + 1] = PixelClip(g);
        dstp[x * rgb_pixel_step + 2] = PixelClip(r);
        if constexpr (rgb_pixel_step == 4) { // n/a
          dstp[x * 4 + 3] = 255;
        }
      }
    }
    dstp -= dst_pitch;
    srcY += src_pitch_y;
    srcU += src_pitch_uv;
    srcV += src_pitch_uv;
    if constexpr(hasAlpha)
      srcA += src_pitch_a;
  }
}

//instantiate
//template<int rgb_pixel_step, bool hasAlpha>
template void convert_yv24_to_rgb_avx2<3, false>(BYTE* dstp, const BYTE* srcY, const BYTE* srcU, const BYTE* srcV, const BYTE* srcA, size_t dst_pitch, size_t src_pitch_y, size_t src_pitch_uv, size_t src_pitch_a, size_t width, size_t height, const ConversionMatrix& matrix);
template void convert_yv24_to_rgb_avx2<4, false>(BYTE* dstp, const BYTE* srcY, const BYTE* srcU, const BYTE* srcV, const BYTE* srcA, size_t dst_pitch, size_t src_pitch_y, size_t src_pitch_uv, size_t src_pitch_a, size_t width, size_t height, const ConversionMatrix& matrix);
template void convert_yv24_to_rgb_avx2<3, true>(BYTE* dstp, const BYTE* srcY, const BYTE* srcU, const BYTE* srcV, const BYTE* srcA, size_t dst_pitch, size_t src_pitch_y, size_t src_pitch_uv, size_t src_pitch_a, size_t width, size_t height, const ConversionMatrix& matrix);
template void convert_yv24_to_rgb_avx2<4, true>(BYTE* dstp, const BYTE* srcY, const BYTE* srcU, const BYTE* srcV, const BYTE* srcA, size_t dst_pitch, size_t src_pitch_y, size_t src_pitch_uv, size_t src_pitch_a, size_t width, size_t height, const ConversionMatrix& matrix);


template<int bits_per_pixel>
void convert_planarrgb_to_yuv_uint16_avx2(BYTE* (&dstp)[3], int(&dstPitch)[3], const BYTE* (&srcp)[3], const int(&srcPitch)[3], int width, int height, const ConversionMatrix& m)
{
  // generic for 10-16 bit uint16 
  // originally made only for 16 bits where unsigned 16 arithmetic makes things difficult

  __m256  half_f = _mm256_set1_ps((float)(1u << (bits_per_pixel - 1)));
  __m128i limit = _mm_set1_epi16((short)((1 << bits_per_pixel) - 1)); // 255
  __m256 offset_f = _mm256_set1_ps(m.offset_y_f);
  __m256i offset_rgb = _mm256_set1_epi32(m.offset_rgb);

  const bool has_offset_rgb = 0 != m.offset_rgb_f;

  //__m128i zero = _mm_setzero_si128();

  const int rowsize = width * sizeof(uint16_t);
  for (int yy = 0; yy < height; yy++) {
    for (int x = 0; x < rowsize; x += 8 * sizeof(uint16_t)) {
      __m256 g, b, r;
      // uint16_t: load 16 bytes: 8 pixels

      __m128i gi = _mm_load_si128(reinterpret_cast<const __m128i*>(srcp[0] + x));
      __m128i bi = _mm_load_si128(reinterpret_cast<const __m128i*>(srcp[1] + x));
      __m128i ri = _mm_load_si128(reinterpret_cast<const __m128i*>(srcp[2] + x));
      
      __m256i gi32 = _mm256_cvtepu16_epi32(gi);
      __m256i bi32 = _mm256_cvtepu16_epi32(bi);
      __m256i ri32 = _mm256_cvtepu16_epi32(ri);
      if (has_offset_rgb) {
        bi32 = _mm256_add_epi32(bi32, offset_rgb);
        gi32 = _mm256_add_epi32(gi32, offset_rgb);
        ri32 = _mm256_add_epi32(ri32, offset_rgb);
      }
      g = _mm256_cvtepi32_ps(gi32);
      b = _mm256_cvtepi32_ps(bi32);
      r = _mm256_cvtepi32_ps(ri32);
      /*
      int Y = m.offset_y + (int)(((sum_t)m.y_b * b + (sum_t)m.y_g * g + (sum_t)m.y_r * r + 16384)>>15);
      int U = half + (int)(((sum_t)m.u_b * b + (sum_t)m.u_g * g + (sum_t)m.u_r * r + 16384) >> 15);
      int V = half + (int)(((sum_t)m.v_b * b + (sum_t)m.v_g * g + (sum_t)m.v_r * r + 16384) >> 15);
      */
      // *Y*
      {
        auto mat_r = _mm256_set1_ps(m.y_r_f);
        auto mat_g = _mm256_set1_ps(m.y_g_f);
        auto mat_b = _mm256_set1_ps(m.y_b_f);
        __m256 y = _mm256_fmadd_ps(r, mat_r, _mm256_fmadd_ps(g, mat_g, _mm256_fmadd_ps(b, mat_b, offset_f)));
        __m256i yi = _mm256_cvtps_epi32(y); // no extra rounding, cvtps rounds to nearest
        yi = _mm256_packus_epi32(yi, _mm256_setzero_si256()); // 16x uint16_t
        __m128i res = _mm256_castsi256_si128(_mm256_permute4x64_epi64(yi, 0xD8));
        if constexpr (bits_per_pixel < 16) // albeit 10-14 bit have another function, make this general
          res = _mm_min_epi16(res, limit); // clamp 10,12,14 bit
        _mm_store_si128(reinterpret_cast<__m128i*>(dstp[0] + x), res);
      }
      // *U*
      {
        auto mat_r = _mm256_set1_ps(m.u_r_f);
        auto mat_g = _mm256_set1_ps(m.u_g_f);
        auto mat_b = _mm256_set1_ps(m.u_b_f);
        __m256 y = _mm256_fmadd_ps(r, mat_r, _mm256_fmadd_ps(g, mat_g, _mm256_fmadd_ps(b, mat_b, half_f)));
        __m256i yi = _mm256_cvtps_epi32(y); // no extra rounding, cvtps rounds to nearest
        yi = _mm256_packus_epi32(yi, _mm256_setzero_si256()); // 16x uint16_t
        __m128i res = _mm256_castsi256_si128(_mm256_permute4x64_epi64(yi, 0xD8));
        if constexpr (bits_per_pixel < 16) // albeit 10-14 bit have another function, make this general
          res = _mm_min_epi16(res, limit); // clamp 10,12,14 bit
        _mm_store_si128(reinterpret_cast<__m128i*>(dstp[1] + x), res);
      }
      // *V*
      {
        auto mat_r = _mm256_set1_ps(m.v_r_f);
        auto mat_g = _mm256_set1_ps(m.v_g_f);
        auto mat_b = _mm256_set1_ps(m.v_b_f);
        __m256 y = _mm256_fmadd_ps(r, mat_r, _mm256_fmadd_ps(g, mat_g, _mm256_fmadd_ps(b, mat_b, half_f)));
        __m256i yi = _mm256_cvtps_epi32(y); // no extra rounding, cvtps rounds to nearest
        yi = _mm256_packus_epi32(yi, _mm256_setzero_si256()); // 16x uint16_t
        __m128i res = _mm256_castsi256_si128(_mm256_permute4x64_epi64(yi, 0xD8));
        if constexpr (bits_per_pixel < 16) // albeit 10-14 bit have another function, make this general
          res = _mm_min_epi16(res, limit); // clamp 10,12,14 bit
        _mm_store_si128(reinterpret_cast<__m128i*>(dstp[2] + x), res);
      }
    }
    srcp[0] += srcPitch[0];
    srcp[1] += srcPitch[1];
    srcp[2] += srcPitch[2];
    dstp[0] += dstPitch[0];
    dstp[1] += dstPitch[1];
    dstp[2] += dstPitch[2];
  }
}

// Instantiate them
template void convert_planarrgb_to_yuv_uint16_avx2<10>(BYTE* (&dstp)[3], int(&dstPitch)[3], const BYTE* (&srcp)[3], const int(&srcPitch)[3], int width, int height, const ConversionMatrix& m);
template void convert_planarrgb_to_yuv_uint16_avx2<12>(BYTE* (&dstp)[3], int(&dstPitch)[3], const BYTE* (&srcp)[3], const int(&srcPitch)[3], int width, int height, const ConversionMatrix& m);
template void convert_planarrgb_to_yuv_uint16_avx2<14>(BYTE* (&dstp)[3], int(&dstPitch)[3], const BYTE* (&srcp)[3], const int(&srcPitch)[3], int width, int height, const ConversionMatrix& m);
template void convert_planarrgb_to_yuv_uint16_avx2<16>(BYTE* (&dstp)[3], int(&dstPitch)[3], const BYTE* (&srcp)[3], const int(&srcPitch)[3], int width, int height, const ConversionMatrix& m);

#define XP_LAMBDA_CAPTURE_FIX(x) (void)(x)

template<typename pixel_t, bool lessthan16bit, bool lessthan16bit_target, typename pixel_t_dst, YuvRgbConversionType conv_type>
static void convert_yuv_to_planarrgb_uintN_avx2_internal(BYTE* (&dstp)[3], int(&dstPitch)[3], const BYTE* (&srcp)[3], const int(&srcPitch)[3], int width, int height, const ConversionMatrix& m,
  const int bits_per_pixel, int bits_per_pixel_target)
{
  // 8 bit        uint8_t
  // 10,12,14 bit uint16_t (signed range)
  // 16 bit logic:
  /*
    1. pivot the pixel:
      Convert uint16 pixel Y to a signed int16 by flipping the MSB (this is mathematically identical to Y−32768).
        0 becomes −32768.
        65535 becomes 32767.
        All values now fit in int16 without saturation.
    2. reorganize the formula (offset is negative when exists, e.g. -4096):
       Original: (Y + offset_y) * Cy
       Substitute
       Y = (Ysigned​ + 32768)
       => (Ysigned​ + 32768 + offset_y) * Cy​
       => (Ysigned​ * Cy​) + (32768 + offset_y) * Cy​)
    3. correction after the madd section:
       Add ((32768 + offset_y) * Cy​) to the existing 32-bit rounding/offset constant.
    4. Output rgb offset is also added to the precalculated patch.
  */
  //const bool full_s = m.offset_y == 0;
  //const bool full_d = m.offset_rgb == 0;
  constexpr bool need_float_conversion = conv_type == YuvRgbConversionType::FLOAT_OUTPUT;
  constexpr bool need_int_conversion_narrow_range = conv_type == YuvRgbConversionType::BITCONV_INT_LIMITED;       // full_d is false
  constexpr bool need_int_conversion_full_range = conv_type == YuvRgbConversionType::BITCONV_INT_FULL; // full_d is true
  constexpr bool need_int_conversion = conv_type == YuvRgbConversionType::BITCONV_INT_FULL || conv_type == YuvRgbConversionType::BITCONV_INT_LIMITED;
  if constexpr (!need_int_conversion)
    bits_per_pixel_target = bits_per_pixel; // make it quasi constexpr for optimizer
  const int bit_diff = need_int_conversion ? bits_per_pixel_target - bits_per_pixel : 0;
  const int target_shift = need_int_conversion_narrow_range ? 13 - bit_diff : 13; // int->int narrow range: integrate the bit depth conversion into the scaling back

  const int ROUNDER = (need_float_conversion || need_int_conversion_full_range) ? 0 : (1 << (target_shift - 1)); // 0 when float internal calculation is involved

  const int round_mask_plus_rgb_offset_i = need_float_conversion ? 0 : ROUNDER + (m.offset_rgb << 13); // latter of same magnitude as coeffs, also in simd madd
  const float rgb_offset_f = m.offset_rgb_f_32; // for post-32-bit float conversion

  int half_pixel_offset;
  int max_pixel_value_source, max_pixel_value_target;
  if constexpr (sizeof(pixel_t) == 1) {
    // 8 bit quasi constexpr
    half_pixel_offset = 128;
    max_pixel_value_source = 255;
  }
  else {
    half_pixel_offset = 1 << (bits_per_pixel - 1);
    max_pixel_value_source = (1 << bits_per_pixel) - 1;
  }

  if constexpr (sizeof(pixel_t_dst) == 1) {
    max_pixel_value_target = 255;
  }
  else {
    max_pixel_value_target = (1 << bits_per_pixel_target) - 1;
  }

  constexpr int int_arithmetic_shift = 1 << 13;
  const float scale_f =
    need_float_conversion ? 1.0f / static_cast<float>(int_arithmetic_shift * m.target_span_f / m.target_span_f_32) : // X->X bits appear as float, before going to real 32-bit range
    need_int_conversion_full_range ? (float)max_pixel_value_target / max_pixel_value_source :
    1.0f; // n/a

  __m256i half = _mm256_set1_epi16((short)half_pixel_offset);  // 128
  __m256i limit = _mm256_set1_epi16((short)max_pixel_value_target); // 255
  __m256i offset;
  if constexpr (need_int_conversion_full_range)
    offset = _mm256_set1_epi32(m.offset_y);
  else
    offset = _mm256_set1_epi16((short)m.offset_y);

  // to be able to use it as signed 16 bit in madd; 4096+(16<<13) would not fit into i16
   // multiplier is 4096 instead of 1:
   // original   : 1      * 4096
   // needed     : 1      * (4096 + offset_rgb<<13)  (4096 + 131072 overflows i16)
   // changed to : 4096   * (1 + offset_rgb>>(13-1)
   // Except for exact 16 bit, where we do the output offset adjustment along with the 16 bit signed-unsigned pivot fix
  constexpr int ROUND_SCALE = 4096; // 1 << 12, for 13 bit integer arithmetic: "0.5"
  const __m256i m256i_round_scale = _mm256_set1_epi16(ROUND_SCALE);

  // Float conversion extra rules
  // - no rounding
  // - no integer scaling back 13 bits
  // - no clamping to bit depth limits
  // - the 13-bit scaling factor is integrated into the bits_per_pixel shift

  int round_mask_plus_rgb_offset_scaled_i;
  __m256i v_patch_G, v_patch_B, v_patch_R;
  __m256i sign_flip_mask = _mm256_set1_epi16((short)0x8000); // for 16 bit pivot

  // Full-range is full-float inside, so no need to do any of the above adjustments in integer arithmetic, just do the full float conversion and clamp at the end if needed
  if constexpr (!need_int_conversion_full_range) {
    if constexpr (lessthan16bit) {
      // 8-14 bit
      round_mask_plus_rgb_offset_scaled_i = round_mask_plus_rgb_offset_i / ROUND_SCALE;
      v_patch_G = v_patch_B = v_patch_R = _mm256_setzero_si256(); // No patch needed
    }
    else {
      // exact 16 bit
      // keep madd simple: only handle the rounding (ROUND_SCALE * 1)
      // rgb offset is handled later in the patch, after the madd, added to the same place as the pivot adjustment
      round_mask_plus_rgb_offset_scaled_i = ROUNDER / ROUND_SCALE; // effectively 1 or 0 (need_float_conversion)

      // move BOTH the pivot and the output RGB offset to the 32-bit patch
      // Since we have to do the patching anyway, we can combine both adjustments here
      const int luma_pivot = 32768 + m.offset_y;
      const int rgb_out_offset = need_float_conversion ? 0 : (m.offset_rgb << 13); // 32-bit post-conversion adds offset in float domain.

      // total patch = (pivot * coeff) + output RGB offset
      v_patch_G = _mm256_set1_epi32(luma_pivot * m.y_g + rgb_out_offset);
      v_patch_B = _mm256_set1_epi32(luma_pivot * m.y_b + rgb_out_offset);
      v_patch_R = _mm256_set1_epi32(luma_pivot * m.y_r + rgb_out_offset);
    }
  }

  const __m256 rgb_offset_f_avx2 = _mm256_set1_ps(rgb_offset_f);

  __m256i zero = _mm256_setzero_si256();
  const __m256 scale_f_avx2 = _mm256_set1_ps(scale_f);

  // not used in full-range conversion
  // all other int cases use madd: No bitdepth conversion, int-to-float, or int-to-int with narrow range scaling
  __m256i m_uy_G, m_vr_G, m_uy_B, m_vr_B, m_uy_R, m_vr_R; // integer arithmetic, including 32-bit float target
  __m256 m_y_g_f, m_y_b_f, m_y_r_f, m_u_g_f, m_u_b_f, m_u_r_f, m_v_g_f, m_v_b_f, m_v_r_f, m_offset_rgb_f; // full-range float-inside conversion

  if (!need_int_conversion_full_range) {
    // for 16 bit, the u/v coeffs are the same for G and R, and G and B respectively, so we can pack them together to save registers
    m_uy_G = _mm256_set1_epi32((static_cast<uint16_t>(m.y_g) << 16) | static_cast<uint16_t>(m.u_g));
    m_vr_G = _mm256_set1_epi32((static_cast<uint16_t>(round_mask_plus_rgb_offset_scaled_i) << 16) | static_cast<uint16_t>(m.v_g));
    m_uy_B = _mm256_set1_epi32((static_cast<uint16_t>(m.y_b) << 16) | static_cast<uint16_t>(m.u_b));
    m_vr_B = _mm256_set1_epi32((static_cast<uint16_t>(round_mask_plus_rgb_offset_scaled_i) << 16) | static_cast<uint16_t>(m.v_b));
    m_uy_R = _mm256_set1_epi32((static_cast<uint16_t>(m.y_r) << 16) | static_cast<uint16_t>(m.u_r));
    m_vr_R = _mm256_set1_epi32((static_cast<uint16_t>(round_mask_plus_rgb_offset_scaled_i) << 16) | static_cast<uint16_t>(m.v_r));
  }
  else {
    // and in full-range float-inside:
    m_y_g_f = _mm256_mul_ps(_mm256_set1_ps(m.y_g_f), scale_f_avx2);
    m_y_b_f = _mm256_mul_ps(_mm256_set1_ps(m.y_b_f), scale_f_avx2);
    m_y_r_f = _mm256_mul_ps(_mm256_set1_ps(m.y_r_f), scale_f_avx2);
    m_u_g_f = _mm256_mul_ps(_mm256_set1_ps(m.u_g_f), scale_f_avx2);
    m_u_b_f = _mm256_mul_ps(_mm256_set1_ps(m.u_b_f), scale_f_avx2);
    m_u_r_f = _mm256_mul_ps(_mm256_set1_ps(m.u_r_f), scale_f_avx2);
    m_v_g_f = _mm256_mul_ps(_mm256_set1_ps(m.v_g_f), scale_f_avx2);
    m_v_b_f = _mm256_mul_ps(_mm256_set1_ps(m.v_b_f), scale_f_avx2);
    m_v_r_f = _mm256_mul_ps(_mm256_set1_ps(m.v_r_f), scale_f_avx2);
    m_offset_rgb_f = _mm256_mul_ps(_mm256_set1_ps(m.offset_rgb_f), scale_f_avx2);
  }

  const int rowsize = width * sizeof(pixel_t);
  for (int yy = 0; yy < height; yy++) {
    // 32 pixels per loop: 32 * 8 bit = 32 bytes; 32 * 16 bit = 64 bytes

    // FIXME: when need_float_conversion, then processing 32 pixels at once breaks Avisynth+
    // FRAME_ALIGN == 64, so we have to separate the processing into 32 pixel then 16 pixel chunks.
    // See ConvertBits to float AVX2 version for reference.

    for (int x = 0; x < rowsize; x += 32 * sizeof(pixel_t)) {
      __m256i y1, u1, v1, y2, u2, v2;

      // Load and pivot for 16 bit if needed
      if constexpr (sizeof(pixel_t) == 1) {
        // Load 32 bytes (32 pixels), unpack to two 16-bit registers
        __m256i y_raw = _mm256_load_si256(reinterpret_cast<const __m256i*>(srcp[0] + x)); // 0..7 8..15 16..23 24..31
        __m256i u_raw = _mm256_load_si256(reinterpret_cast<const __m256i*>(srcp[1] + x));
        __m256i v_raw = _mm256_load_si256(reinterpret_cast<const __m256i*>(srcp[2] + x));
        if constexpr (sizeof(pixel_t_dst) == 2) {
          // 8->16 bit: pre-shuffle to avoid permutes at the end
          y_raw = _mm256_permute4x64_epi64(y_raw, 0xD8); // (0 << 0) | (2 << 2) | (1 << 4) | (3 << 6)
          // 0..7, 16..23, 8..15, 24..31
          u_raw = _mm256_permute4x64_epi64(u_raw, 0xD8);
          v_raw = _mm256_permute4x64_epi64(v_raw, 0xD8);
        }
        y1 = _mm256_unpacklo_epi8(y_raw, zero); y2 = _mm256_unpackhi_epi8(y_raw, zero);
        // y1: 0..7, 8..15, y2: 16..23, 24..31
        u1 = _mm256_unpacklo_epi8(u_raw, zero); u2 = _mm256_unpackhi_epi8(u_raw, zero);
        v1 = _mm256_unpacklo_epi8(v_raw, zero); v2 = _mm256_unpackhi_epi8(v_raw, zero);
      }
      else {
        // Load 64 bytes (32 pixels)
        y1 = _mm256_load_si256(reinterpret_cast<const __m256i*>(srcp[0] + x));
        u1 = _mm256_load_si256(reinterpret_cast<const __m256i*>(srcp[1] + x));
        v1 = _mm256_load_si256(reinterpret_cast<const __m256i*>(srcp[2] + x));
        y2 = _mm256_load_si256(reinterpret_cast<const __m256i*>(srcp[0] + x + 32));
        u2 = _mm256_load_si256(reinterpret_cast<const __m256i*>(srcp[1] + x + 32));
        v2 = _mm256_load_si256(reinterpret_cast<const __m256i*>(srcp[2] + x + 32));
      }

      // Avisynth FRAME_ALIGN == 64, so when need_float_conversion, we cannot process 32 pixels at once at the end of the line
      // Can I write all 32 pixels (blockA+B+C+D)?
      // yes: write all 4 blocks (blockA, B, C, D)
      // no : Write only first 2 blocks (blockA, B), 16 floats are guaranteed safe due to 64-byte alignment
      bool safe_last_16float_64bytes;
      if constexpr (need_float_conversion)
        safe_last_16float_64bytes = (x + 32 * sizeof(pixel_t)) <= rowsize;
      else
        safe_last_16float_64bytes = false;

      if constexpr (need_int_conversion_full_range) {
        // float workflow
        __m256i y1_32_lo = _mm256_add_epi32(_mm256_unpacklo_epi16(y1, zero), offset);
        __m256i y1_32_hi = _mm256_add_epi32(_mm256_unpackhi_epi16(y1, zero), offset);
        __m256i y2_32_lo = _mm256_add_epi32(_mm256_unpacklo_epi16(y2, zero), offset);
        __m256i y2_32_hi = _mm256_add_epi32(_mm256_unpackhi_epi16(y2, zero), offset);
        // this is SSE2, from SSE4.1 we could use _mm_cvtepi16_epi32 directly but it makes packus lane crossing correction difficult
        // and from avx512 _mm512_cvtepi16_ps (only signed 16 exists)
        // Create a register where each lane is 0xFFFF if negative, 0x0000 if positive
        u1 = _mm256_sub_epi16(u1, half); // no saturation!
        u2 = _mm256_sub_epi16(u2, half);
        __m256i u1_sign = _mm256_srai_epi16(u1, 15);
        __m256i u2_sign = _mm256_srai_epi16(u2, 15);
        v1 = _mm256_sub_epi16(v1, half);
        v2 = _mm256_sub_epi16(v2, half);
        __m256i v1_sign = _mm256_srai_epi16(v1, 15);
        __m256i v2_sign = _mm256_srai_epi16(v2, 15);
        // or cvtepi16_epi32 w/o sign ? does it keep lanes
        __m256i u1_32_lo = _mm256_unpacklo_epi16(u1, u1_sign);
        __m256i u1_32_hi = _mm256_unpackhi_epi16(u1, u1_sign);
        __m256i v1_32_lo = _mm256_unpacklo_epi16(v1, v1_sign);
        __m256i v1_32_hi = _mm256_unpackhi_epi16(v1, v1_sign);
        __m256i u2_32_lo = _mm256_unpacklo_epi16(u2, u2_sign);
        __m256i u2_32_hi = _mm256_unpackhi_epi16(u2, u2_sign);
        __m256i v2_32_lo = _mm256_unpacklo_epi16(v2, v2_sign);
        __m256i v2_32_hi = _mm256_unpackhi_epi16(v2, v2_sign);

        __m256 y1_f_lo = _mm256_cvtepi32_ps(y1_32_lo);
        __m256 y1_f_hi = _mm256_cvtepi32_ps(y1_32_hi);
        __m256 u1_f_lo = _mm256_cvtepi32_ps(u1_32_lo);
        __m256 u1_f_hi = _mm256_cvtepi32_ps(u1_32_hi);
        __m256 v1_f_lo = _mm256_cvtepi32_ps(v1_32_lo);
        __m256 v1_f_hi = _mm256_cvtepi32_ps(v1_32_hi);

        __m256 y2_f_lo = _mm256_cvtepi32_ps(y2_32_lo);
        __m256 y2_f_hi = _mm256_cvtepi32_ps(y2_32_hi);
        __m256 u2_f_lo = _mm256_cvtepi32_ps(u2_32_lo);
        __m256 u2_f_hi = _mm256_cvtepi32_ps(u2_32_hi);
        __m256 v2_f_lo = _mm256_cvtepi32_ps(v2_32_lo);
        __m256 v2_f_hi = _mm256_cvtepi32_ps(v2_32_hi);


        /*
        b_f = matrix.y_b_f * Y + matrix.u_b_f * U + matrix.v_b_f * V + matrix.offset_rgb_f;
        g_f = matrix.y_g_f * Y + matrix.u_g_f * U + matrix.v_g_f * V + matrix.offset_rgb_f;
        r_f = matrix.y_r_f * Y + matrix.u_r_f * U + matrix.v_r_f * V + matrix.offset_rgb_f;
        */
        // Blue 0-15
        __m256 b1_f_lo = _mm256_fmadd_ps(m_v_b_f, v1_f_lo, _mm256_fmadd_ps(m_u_b_f, u1_f_lo, _mm256_fmadd_ps(m_y_b_f, y1_f_lo, m_offset_rgb_f)));
        __m256 b1_f_hi = _mm256_fmadd_ps(m_v_b_f, v1_f_hi, _mm256_fmadd_ps(m_u_b_f, u1_f_hi, _mm256_fmadd_ps(m_y_b_f, y1_f_hi, m_offset_rgb_f)));
        // Green 0-15
        __m256 g1_f_lo = _mm256_fmadd_ps(m_v_g_f, v1_f_lo, _mm256_fmadd_ps(m_u_g_f, u1_f_lo, _mm256_fmadd_ps(m_y_g_f, y1_f_lo, m_offset_rgb_f)));
        __m256 g1_f_hi = _mm256_fmadd_ps(m_v_g_f, v1_f_hi, _mm256_fmadd_ps(m_u_g_f, u1_f_hi, _mm256_fmadd_ps(m_y_g_f, y1_f_hi, m_offset_rgb_f)));
        // Red 0-15
        __m256 r1_f_lo = _mm256_fmadd_ps(m_v_r_f, v1_f_lo, _mm256_fmadd_ps(m_u_r_f, u1_f_lo, _mm256_fmadd_ps(m_y_r_f, y1_f_lo, m_offset_rgb_f)));
        __m256 r1_f_hi = _mm256_fmadd_ps(m_v_r_f, v1_f_hi, _mm256_fmadd_ps(m_u_r_f, u1_f_hi, _mm256_fmadd_ps(m_y_r_f, y1_f_hi, m_offset_rgb_f)));
        // Blue 16-31
        __m256 b2_f_lo = _mm256_fmadd_ps(m_v_b_f, v2_f_lo, _mm256_fmadd_ps(m_u_b_f, u2_f_lo, _mm256_fmadd_ps(m_y_b_f, y2_f_lo, m_offset_rgb_f)));
        __m256 b2_f_hi = _mm256_fmadd_ps(m_v_b_f, v2_f_hi, _mm256_fmadd_ps(m_u_b_f, u2_f_hi, _mm256_fmadd_ps(m_y_b_f, y2_f_hi, m_offset_rgb_f)));
        // Green 16-31
        __m256 g2_f_lo = _mm256_fmadd_ps(m_v_g_f, v2_f_lo, _mm256_fmadd_ps(m_u_g_f, u2_f_lo, _mm256_fmadd_ps(m_y_g_f, y2_f_lo, m_offset_rgb_f)));
        __m256 g2_f_hi = _mm256_fmadd_ps(m_v_g_f, v2_f_hi, _mm256_fmadd_ps(m_u_g_f, u2_f_hi, _mm256_fmadd_ps(m_y_g_f, y2_f_hi, m_offset_rgb_f)));
        // Red 16-31
        __m256 r2_f_lo = _mm256_fmadd_ps(m_v_r_f, v2_f_lo, _mm256_fmadd_ps(m_u_r_f, u2_f_lo, _mm256_fmadd_ps(m_y_r_f, y2_f_lo, m_offset_rgb_f)));
        __m256 r2_f_hi = _mm256_fmadd_ps(m_v_r_f, v2_f_hi, _mm256_fmadd_ps(m_u_r_f, u2_f_hi, _mm256_fmadd_ps(m_y_r_f, y2_f_hi, m_offset_rgb_f)));

        /* already done in preparation
        // x bits integer arithmetic shift and the bit_depth correction is in scale_f.
        // In sse2 we pre-multiplied the matrix by scale_f_sse2, so no need to do it again here.
        // g_f = g_f * scale_f;
        // b_f = b_f * scale_f;
        // r_f = r_f * scale_f;
        b_f_lo = _mm_mul_ps(b_f_lo, scale_f_sse2);
        b_f_hi = _mm_mul_ps(b_f_hi, scale_f_sse2);
        g_f_lo = _mm_mul_ps(g_f_lo, scale_f_sse2);
        g_f_hi = _mm_mul_ps(g_f_hi, scale_f_sse2);
        r_f_lo = _mm_mul_ps(r_f_lo, scale_f_sse2);
        r_f_hi = _mm_mul_ps(r_f_hi, scale_f_sse2);
        */

        auto process_from_float_plane_avx2 = [&](BYTE* plane_ptr, __m256 lo_1, __m256 hi_1, __m256 lo_2, __m256 hi_2) {
          XP_LAMBDA_CAPTURE_FIX(zero);
          XP_LAMBDA_CAPTURE_FIX(limit);
          /*
          g = static_cast<int>(g_f + 0.5f);
          b = static_cast<int>(b_f + 0.5f);
          r = static_cast<int>(r_f + 0.5f);
          */
          __m256 float_rounder = _mm256_set1_ps(0.5f);
          __m256i res1_lo = _mm256_cvttps_epi32(_mm256_add_ps(lo_1, float_rounder));
          __m256i res1_hi = _mm256_cvttps_epi32(_mm256_add_ps(hi_1, float_rounder));
          __m256i res2_lo = _mm256_cvttps_epi32(_mm256_add_ps(lo_2, float_rounder));
          __m256i res2_hi = _mm256_cvttps_epi32(_mm256_add_ps(hi_2, float_rounder));
          const int pix_idx = x * sizeof(pixel_t_dst) / sizeof(pixel_t);

          // unlike SSE2, AVX2 has packus_epi32
          __m256i p1 = _mm256_packus_epi32(res1_lo, res1_hi);
          __m256i p2 = _mm256_packus_epi32(res2_lo, res2_hi);

          if constexpr (sizeof(pixel_t_dst) == 1) {
            // X->8 bits
            __m256i final8 = _mm256_packus_epi16(p1, p2);
            if constexpr (sizeof(pixel_t) == 2) {
              // 16->8 extra shuffle
              final8 = _mm256_permute4x64_epi64(final8, (0 << 0) | (2 << 2) | (1 << 4) | (3 << 6));
            }
            _mm256_store_si256(reinterpret_cast<__m256i*>(plane_ptr + pix_idx), final8);
          }
          else if constexpr (sizeof(pixel_t_dst) == 2) { 
            // x->16 bits
            if constexpr (lessthan16bit_target) {
              p1 = _mm256_min_epi16(p1, limit);
              p2 = _mm256_min_epi16(p2, limit);
            }
            if constexpr (sizeof(pixel_t) == 1) {
              // 8->16 bits:
              /* If we'd not pre-shuffled the input for 8->16 bit, we would need to do this shuffle at the end to heal the lanes:
              // p1: 0-7  16-23; p2: 8-15 24-31
              // two shuffles per R G B plane would be needed, so 6 total, Pre-shuffling Y U V is only 3
              __m256i temp_p1 = _mm256_permute2x128_si256(p1, p2, 0x20); // 0-7, 8-15
              __m256i temp_p2 = _mm256_permute2x128_si256(p1, p2, 0x31); // 16-23, 24-31
              p1 = temp_p1;
              p2 = temp_p2;
              */
            }

            _mm256_store_si256(reinterpret_cast<__m256i*>(plane_ptr + pix_idx), p1);
            _mm256_store_si256(reinterpret_cast<__m256i*>(plane_ptr + pix_idx + 32), p2);
          }
          };

        process_from_float_plane_avx2(dstp[0], g1_f_lo, g1_f_hi, g2_f_lo, g2_f_hi);
        process_from_float_plane_avx2(dstp[1], b1_f_lo, b1_f_hi, b2_f_lo, b2_f_hi);
        process_from_float_plane_avx2(dstp[2], r1_f_lo, r1_f_hi, r2_f_lo, r2_f_hi);

      }
      else {
        if constexpr (lessthan16bit) {
          y1 = _mm256_adds_epi16(y1, offset); // add, because offset is negative
          y2 = _mm256_adds_epi16(y2, offset); // add, because offset is negative
        }
        else {
          // make unsigned to signed by flipping MSB
          // pivot the luma, adjust the offset separately
          y1 = _mm256_xor_si256(y1, sign_flip_mask);
          y2 = _mm256_xor_si256(y2, sign_flip_mask);
        }
        u1 = _mm256_sub_epi16(u1, half);
        u2 = _mm256_sub_epi16(u2, half);
        v1 = _mm256_sub_epi16(v1, half);
        v2 = _mm256_sub_epi16(v2, half);

        // pre-unpack MADD pairs
        // These are common for all R, G, B planes
        // Pair 1: [U0 Y0 U1 Y1 ... | U8 Y8 U9 Y9 ...]
        __m256i uy1_lo = _mm256_unpacklo_epi16(u1, y1);
        __m256i uy1_hi = _mm256_unpackhi_epi16(u1, y1);
        __m256i uy2_lo = _mm256_unpacklo_epi16(u2, y2);
        __m256i uy2_hi = _mm256_unpackhi_epi16(u2, y2);

        // Pair 2: [V0 Rnd V1 Rnd ... | V8 Rnd V9 Rnd ...]
        __m256i vr1_lo = _mm256_unpacklo_epi16(v1, m256i_round_scale);
        __m256i vr1_hi = _mm256_unpackhi_epi16(v1, m256i_round_scale);
        __m256i vr2_lo = _mm256_unpacklo_epi16(v2, m256i_round_scale);
        __m256i vr2_hi = _mm256_unpackhi_epi16(v2, m256i_round_scale);


        // for 16 bit, the rgb_offset is merged into the post patch adjustment
        // 13 bit fixed point arithmetic rounder 0.5 is 4096.
        // Need1:  (m.y_b   m.u_b )     (m.y_b   m.u_b)     (m.y_b   m.u_b)     (m.y_b   m.u_b)   8x16 bit
        //         (  y3      u3  )     (  y2      u2 )     (  y1      u1 )     (   y0     u0 )   8x16 bit
        // res1=  (y_b*y3 + u_b*u3)   ...                                                         4x32 bit
        // Need2:  (m.v_b   round')     (m.y_b   round')     (m.y_b   round')     (m.y_b   round')
        //         (  v3     4096 )     (  v2     4096 )     (  v1     4096 )     (  v0     4096 )
        // res2=  (yv_b*v3 + round' )  ...  round' = round + rgb_offset

        // Processing lambda - checked and benchmarked to be inlined nicely -avoids code bloat

        // For v141_xp compatibility: forces the compiler to capture a const variable
        // that would otherwise be optimized out of nested lambda scopes.

        auto process_plane = [&](BYTE* plane_ptr, __m256i m_uy, __m256i m_vr, __m256i v_patch) {
          XP_LAMBDA_CAPTURE_FIX(limit);
          XP_LAMBDA_CAPTURE_FIX(safe_last_16float_64bytes);
          auto madd_scale = [&](__m256i uy, __m256i vr) {
            XP_LAMBDA_CAPTURE_FIX(v_patch);
            XP_LAMBDA_CAPTURE_FIX(target_shift);
            __m256i sum = _mm256_add_epi32(_mm256_madd_epi16(m_uy, uy), _mm256_madd_epi16(m_vr, vr));
            // 16-bit adjustment (signed patch, offset, output rgb offset)
            if constexpr (!lessthan16bit) sum = _mm256_add_epi32(sum, v_patch);
#ifdef XP_TLS
            if (!need_float_conversion)
#else
            if constexpr (!need_float_conversion)
#endif
              sum = _mm256_srai_epi32(sum, target_shift); // 13 bit fixed point shift
            return sum;
            };

          __m256i res1_lo = madd_scale(uy1_lo, vr1_lo); // Pixels 0-3, 8-11
          __m256i res1_hi = madd_scale(uy1_hi, vr1_hi); // Pixels 4-7, 12-15
          __m256i res2_lo = madd_scale(uy2_lo, vr2_lo); // Pixels 16-19, 24-27
          __m256i res2_hi = madd_scale(uy2_hi, vr2_hi); // Pixels 20-23, 28-31

#ifdef XP_TLS
          if (need_float_conversion) {
#else
          if constexpr (need_float_conversion) {
#endif
            // when float output is needed, convert after scaling, mimic a post-ConvertBits
            const int pix_idx = x / sizeof(pixel_t);
            float* f_dst = reinterpret_cast<float*>(plane_ptr) + pix_idx;

            // Define the blocks (8 pixels each) correctly by healing the lanes
            __m256i blockA = _mm256_permute2x128_si256(res1_lo, res1_hi, 0x20); // Pixels 0-7
            __m256i blockB, blockC, blockD;

            if constexpr (sizeof(pixel_t) == 1) {
              // 8-bit: The lanes were swapped across res1 and res2
              blockB = _mm256_permute2x128_si256(res2_lo, res2_hi, 0x20); // Pixels 8-15
              blockC = _mm256_permute2x128_si256(res1_lo, res1_hi, 0x31); // Pixels 16-23
              blockD = _mm256_permute2x128_si256(res2_lo, res2_hi, 0x31); // Pixels 24-31
            }
            else {
              // 16-bit: res1 is 0-15, res2 is 16-31
              blockB = _mm256_permute2x128_si256(res1_lo, res1_hi, 0x31); // Pixels 8-15
              blockC = _mm256_permute2x128_si256(res2_lo, res2_hi, 0x20); // Pixels 16-23
              blockD = _mm256_permute2x128_si256(res2_lo, res2_hi, 0x31); // Pixels 24-31
            }

            // Convert and Store
            auto store_block = [&](float* ptr, __m256i b) {
              _mm256_store_ps(ptr, _mm256_fmadd_ps(_mm256_cvtepi32_ps(b), scale_f_avx2, rgb_offset_f_avx2));
              };

            store_block(f_dst, blockA);
            store_block(f_dst + 8, blockB);

            // Safety check: only store the first half of the 32-pixel block if row width allows
            if (safe_last_16float_64bytes) {
              store_block(f_dst + 16, blockC);
              store_block(f_dst + 24, blockD);
            }
          }
          else {
            const int pix_idx = x * sizeof(pixel_t_dst) / sizeof(pixel_t);
            // unlike SSE2, AVX2 has packus_epi32
            __m256i p1 = _mm256_packus_epi32(res1_lo, res1_hi); // Auto-heals lanes, no permute needed
            __m256i p2 = _mm256_packus_epi32(res2_lo, res2_hi);

            if constexpr (sizeof(pixel_t_dst) == 1) {
              // X->8 bits
              __m256i final8 = _mm256_packus_epi16(p1, p2);
              if constexpr (sizeof(pixel_t) == 2) {
                // 16->8 extra shuffle
                final8 = _mm256_permute4x64_epi64(final8, (0 << 0) | (2 << 2) | (1 << 4) | (3 << 6));
              }
              _mm256_store_si256(reinterpret_cast<__m256i*>(plane_ptr + pix_idx), final8);
            }
            else { 
              // x->16 bits
              if constexpr (lessthan16bit_target) {
                p1 = _mm256_min_epi16(p1, limit);
                p2 = _mm256_min_epi16(p2, limit);
              }
              if constexpr (sizeof(pixel_t) == 1) {
                // 8->16 bits:
                /* If we'd not pre-shuffled the input for 8->16 bit, we would need to do this shuffle at the end to heal the lanes:
                // p1: 0-7  16-23; p2: 8-15 24-31
                // two shuffles per R G B plane would be needed, so 6 total, Pre-shuffling Y U V is only 3
                __m256i temp_p1 = _mm256_permute2x128_si256(p1, p2, 0x20); // 0-7, 8-15
                __m256i temp_p2 = _mm256_permute2x128_si256(p1, p2, 0x31); // 16-23, 24-31
                p1 = temp_p1;
                p2 = temp_p2;
                */
              }
              _mm256_store_si256(reinterpret_cast<__m256i*>(plane_ptr + pix_idx), p1);
              _mm256_store_si256(reinterpret_cast<__m256i*>(plane_ptr + pix_idx + 32), p2);
            }
          }
          };

        // Process planes, using pre-packed coefficient, and the 16 bit patch if needed
        process_plane(dstp[0], m_uy_G, m_vr_G, v_patch_G);
        process_plane(dstp[1], m_uy_B, m_vr_B, v_patch_B);
        process_plane(dstp[2], m_uy_R, m_vr_R, v_patch_R);
      }
    }
    srcp[0] += srcPitch[0];
    srcp[1] += srcPitch[1];
    srcp[2] += srcPitch[2];
    dstp[0] += dstPitch[0];
    dstp[1] += dstPitch[1];
    dstp[2] += dstPitch[2];
  }
}

#undef XP_LAMBDA_CAPTURE_FIX

// Further separating cases inside, dispatcher remains relatively simple
template<typename pixel_t_src, bool lessthan16bit>
void convert_yuv_to_planarrgb_uintN_avx2(BYTE* (&dstp)[3], int(&dstPitch)[3], const BYTE* (&srcp)[3], const int(&srcPitch)[3], int width, int height, const ConversionMatrix& m,
  int bits_per_pixel, int bits_per_pixel_target)
{
  const bool need_conversion = bits_per_pixel_target != bits_per_pixel;
  if (!need_conversion) {
    // no conversion, just YUV to RGB
    convert_yuv_to_planarrgb_uintN_avx2_internal<pixel_t_src, lessthan16bit, lessthan16bit, pixel_t_src, YuvRgbConversionType::NATIVE_INT>(dstp, dstPitch, srcp, srcPitch, width, height, m, bits_per_pixel, bits_per_pixel_target);
    return;
  }

  const bool full_d = m.offset_rgb == 0;

  if (bits_per_pixel_target >= 8 && bits_per_pixel <= 16 && bits_per_pixel_target <= 16) {
    // int->int conversion with range conversion (limited<->full), or upscale with no range conversion (full->full or limited->limited)
    if (bits_per_pixel_target == 8) {
      if (full_d)
        convert_yuv_to_planarrgb_uintN_avx2_internal<pixel_t_src, lessthan16bit, true, uint8_t, YuvRgbConversionType::BITCONV_INT_FULL>(dstp, dstPitch, srcp, srcPitch, width, height, m, bits_per_pixel, bits_per_pixel_target);
      else
        convert_yuv_to_planarrgb_uintN_avx2_internal<pixel_t_src, lessthan16bit, true, uint8_t, YuvRgbConversionType::BITCONV_INT_LIMITED>(dstp, dstPitch, srcp, srcPitch, width, height, m, bits_per_pixel, bits_per_pixel_target);
    }
    else if (bits_per_pixel_target < 16) {
      // lessthan16bit_target is false. Need signed pack and clamping
      if (full_d)
        convert_yuv_to_planarrgb_uintN_avx2_internal<pixel_t_src, lessthan16bit, true, uint16_t, YuvRgbConversionType::BITCONV_INT_FULL>(dstp, dstPitch, srcp, srcPitch, width, height, m, bits_per_pixel, bits_per_pixel_target);
      else
        convert_yuv_to_planarrgb_uintN_avx2_internal<pixel_t_src, lessthan16bit, true, uint16_t, YuvRgbConversionType::BITCONV_INT_LIMITED>(dstp, dstPitch, srcp, srcPitch, width, height, m, bits_per_pixel, bits_per_pixel_target);
    }
    else { // == 16
      if (full_d)
        convert_yuv_to_planarrgb_uintN_avx2_internal<pixel_t_src, lessthan16bit, false, uint16_t, YuvRgbConversionType::BITCONV_INT_FULL>(dstp, dstPitch, srcp, srcPitch, width, height, m, bits_per_pixel, bits_per_pixel_target);
      else
        convert_yuv_to_planarrgb_uintN_avx2_internal<pixel_t_src, lessthan16bit, false, uint16_t, YuvRgbConversionType::BITCONV_INT_LIMITED>(dstp, dstPitch, srcp, srcPitch, width, height, m, bits_per_pixel, bits_per_pixel_target);
    }
  }
  else {
    // int->float conversion
    // limited/full is handled automatically through scaling and offsets
    // lessthan16bit_target doesn't matter since float output has no clamping
    convert_yuv_to_planarrgb_uintN_avx2_internal<pixel_t_src, lessthan16bit, true, float, YuvRgbConversionType::FLOAT_OUTPUT>(dstp, dstPitch, srcp, srcPitch, width, height, m, bits_per_pixel, bits_per_pixel_target);
  }
}

//instantiate
//template<typename pixel_t, bool lessthan16bit>
template void convert_yuv_to_planarrgb_uintN_avx2<uint8_t, true>(BYTE* (&dstp)[3], int(&dstPitch)[3], const BYTE* (&srcp)[3], const int(&srcPitch)[3], int width, int height, const ConversionMatrix& m, int bits_per_pixel, int bits_per_pixel_target);
template void convert_yuv_to_planarrgb_uintN_avx2<uint16_t, true>(BYTE* (&dstp)[3], int(&dstPitch)[3], const BYTE* (&srcp)[3], const int(&srcPitch)[3], int width, int height, const ConversionMatrix& m, int bits_per_pixel, int bits_per_pixel_target);
template void convert_yuv_to_planarrgb_uintN_avx2<uint16_t, false>(BYTE* (&dstp)[3], int(&dstPitch)[3], const BYTE* (&srcp)[3], const int(&srcPitch)[3], int width, int height, const ConversionMatrix& m, int bits_per_pixel, int bits_per_pixel_target);

void convert_yuv_to_planarrgb_float_avx2(BYTE* (&dstp)[3], int(&dstPitch)[3], const BYTE* (&srcp)[3], const int(&srcPitch)[3], int width, int height, const ConversionMatrix& m)
{
  // 32 bit float
  // yes, we support "limited" float
  __m256  offset_f = _mm256_set1_ps(m.offset_y_f);
  __m256  offset_rgb = _mm256_set1_ps(m.offset_rgb_f);

  const __m256 mat_yG = _mm256_set1_ps(m.y_g_f);
  const __m256 mat_uG = _mm256_set1_ps(m.u_g_f);
  const __m256 mat_vG = _mm256_set1_ps(m.v_g_f);

  const __m256 mat_yB = _mm256_set1_ps(m.y_b_f);
  const __m256 mat_uB = _mm256_set1_ps(m.u_b_f);
  const __m256 mat_vB = _mm256_set1_ps(m.v_b_f);

  const __m256 mat_yR = _mm256_set1_ps(m.y_r_f);
  const __m256 mat_uR = _mm256_set1_ps(m.u_r_f);
  const __m256 mat_vR = _mm256_set1_ps(m.v_r_f);

  const int rowsize = width * sizeof(float);
  for (int yy = 0; yy < height; yy++) {
    for (int x = 0; x < rowsize; x += 8 * sizeof(float)) {
      __m256 resG, resB, resR;
      __m256 y, u, v;

      // float: load 32 bytes: 8 pixels
      y = _mm256_load_ps(reinterpret_cast<const float*>(srcp[0] + x));
      u = _mm256_load_ps(reinterpret_cast<const float*>(srcp[1] + x));
      v = _mm256_load_ps(reinterpret_cast<const float*>(srcp[2] + x));
      y = _mm256_add_ps(y, offset_f); // offset is negative*/

      // *G*
      resG = _mm256_fmadd_ps(y, mat_yG, offset_rgb); // mul_y
      resG = _mm256_fmadd_ps(u, mat_uG, resG); // mul_u
      resG = _mm256_fmadd_ps(v, mat_vG, resG); // mul_v
      // no clamp
      _mm256_store_ps(reinterpret_cast<float*>(dstp[0] + x), resG);
      // *B*
      resB = _mm256_fmadd_ps(y, mat_yB, offset_rgb); // mul_y
      resB = _mm256_fmadd_ps(u, mat_uB, resB); // mul_u
      resB = _mm256_fmadd_ps(v, mat_vB, resB); // mul_v
      // no clamp
      _mm256_store_ps(reinterpret_cast<float*>(dstp[1] + x), resB);
      // *R*
      resR = _mm256_fmadd_ps(y, mat_yR, offset_rgb); // mul_y
      resR = _mm256_fmadd_ps(u, mat_uR, resR); // mul_u
      resR = _mm256_fmadd_ps(v, mat_vR, resR); // mul_v
      // no clamp
      _mm256_store_ps(reinterpret_cast<float*>(dstp[2] + x), resR);
    }
    srcp[0] += srcPitch[0];
    srcp[1] += srcPitch[1];
    srcp[2] += srcPitch[2];
    dstp[0] += dstPitch[0];
    dstp[1] += dstPitch[1];
    dstp[2] += dstPitch[2];
  }
}

DISABLE_WARNING_POP
