// Avisynth+
// https://avs-plus.net
//
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

#include "avisynth.h"
#include "blend_common_avx2.h"
#include "../blend_common.h"

#include <stdint.h>

#ifdef _MSC_VER
#include <intrin.h>
#else
#include <x86intrin.h>
#endif
#include  <immintrin.h>




template<bool has_mask, typename pixel_t, bool lessthan16bits>
void overlay_blend_avx2_uint(BYTE* p1, const BYTE* p2, const BYTE* mask,
  const int p1_pitch, const int p2_pitch, const int mask_pitch,
  const int width, const int height, const int opacity, const float opacity_f, const int bits_per_pixel)
{
  const auto rounder = _mm256_set1_ps(0.5f);
  const int max_pixel_value = (1 << bits_per_pixel) - 1;
  const auto factor = has_mask ? opacity_f / max_pixel_value : opacity_f;
  const auto factor_v = _mm256_set1_ps(factor);

  constexpr int pixels_per_cycle = 16;
  constexpr int bytes_per_cycle = pixels_per_cycle * sizeof(pixel_t);
  const int wMod = (width * sizeof(pixel_t) / bytes_per_cycle) * bytes_per_cycle;

  // key formula: p1*(1-factor)+p2*factor -> p1 + (p2 - p1) * factor 
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < wMod; x += bytes_per_cycle) {
      if constexpr (lessthan16bits) {
        // 8-14 bits: int16 delta -> float factor -> int16 delta_scaled -->
        // int16 result = p1 + delta_scaled (all in int16 domain except factor!)
        __m256i p1_i16, p2_i16, delta_i16;
        __m256i mask_i16;

        if constexpr (sizeof(pixel_t) == 1) {
          // 8 bits
          __m128i p1_u8 = _mm_loadu_si128(reinterpret_cast<const __m128i*>(p1 + x));
          __m128i p2_u8 = _mm_loadu_si128(reinterpret_cast<const __m128i*>(p2 + x));

          p1_i16 = _mm256_cvtepu8_epi16(p1_u8);
          p2_i16 = _mm256_cvtepu8_epi16(p2_u8);

          if constexpr (has_mask) {
            __m128i mask_u8 = _mm_loadu_si128(reinterpret_cast<const __m128i*>(mask + x));
            mask_i16 = _mm256_cvtepu8_epi16(mask_u8);
          }
        }
        else { // 10/12/14 bits
          p1_i16 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(p1 + x));
          p2_i16 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(p2 + x));

          if constexpr (has_mask) {
            mask_i16 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(mask + x));
          }
        }

        delta_i16 = _mm256_sub_epi16(p2_i16, p1_i16);

        // split for 32-bit processing
        __m128i delta_i16_lo = _mm256_castsi256_si128(delta_i16);
        __m128i delta_i16_hi = _mm256_extracti128_si256(delta_i16, 1);

        // convert _only_ delta and mask to 32-bit
        __m256i delta_lo_i32 = _mm256_cvtepi16_epi32(delta_i16_lo);
        __m256i delta_hi_i32 = _mm256_cvtepi16_epi32(delta_i16_hi);

        // convert _only_ delta to float
        __m256 delta_lo_f = _mm256_cvtepi32_ps(delta_lo_i32);
        __m256 delta_hi_f = _mm256_cvtepi32_ps(delta_hi_i32);

        __m256 delta_scaled_lo, delta_scaled_hi;

        if constexpr (has_mask) {
          __m128i mask_i16_lo = _mm256_castsi256_si128(mask_i16);
          __m128i mask_i16_hi = _mm256_extracti128_si256(mask_i16, 1);

          __m256i mask_lo_i32 = _mm256_cvtepu16_epi32(mask_i16_lo);
          __m256i mask_hi_i32 = _mm256_cvtepu16_epi32(mask_i16_hi);

          __m256 mask_lo_f = _mm256_cvtepi32_ps(mask_lo_i32);
          __m256 mask_hi_f = _mm256_cvtepi32_ps(mask_hi_i32);

          __m256 new_factor_lo = _mm256_mul_ps(mask_lo_f, factor_v);
          __m256 new_factor_hi = _mm256_mul_ps(mask_hi_f, factor_v);

          // scale delta: delta * factor
          delta_scaled_lo = _mm256_mul_ps(delta_lo_f, new_factor_lo);
          delta_scaled_hi = _mm256_mul_ps(delta_hi_f, new_factor_hi);
        }
        else {
          delta_scaled_lo = _mm256_mul_ps(delta_lo_f, factor_v);
          delta_scaled_hi = _mm256_mul_ps(delta_hi_f, factor_v);
        }

        // round and truncate to int32
        delta_scaled_lo = _mm256_add_ps(delta_scaled_lo, rounder);
        delta_scaled_hi = _mm256_add_ps(delta_scaled_hi, rounder);
        __m256i delta_scaled_i32_lo = _mm256_cvttps_epi32(delta_scaled_lo);
        __m256i delta_scaled_i32_hi = _mm256_cvttps_epi32(delta_scaled_hi);

        // scaled delta back to int16
        __m256i delta_scaled_i16 = _mm256_packs_epi32(delta_scaled_i32_lo, delta_scaled_i32_hi);
        delta_scaled_i16 = _mm256_permute4x64_epi64(delta_scaled_i16, 0xD8);

        // +p1
        __m256i result_i16 = _mm256_add_epi16(p1_i16, delta_scaled_i16);

        if constexpr (sizeof(pixel_t) == 1) {
          __m128i i16_lo_128 = _mm256_castsi256_si128(result_i16);
          __m128i i16_hi_128 = _mm256_extracti128_si256(result_i16, 1);
          __m128i result_u8 = _mm_packus_epi16(i16_lo_128, i16_hi_128);

          _mm_storeu_si128(reinterpret_cast<__m128i*>(p1 + x), result_u8);
        }
        else {
          _mm256_storeu_si256(reinterpret_cast<__m256i*>(p1 + x), result_i16);
        }
      }
      else { // 16-bit overflow-safe, int32 and float
        static_assert(sizeof(pixel_t) == 2, "lessthan16bits=false only for uint16_t");

        // 16xuint16_t pixels
        __m256i p1_u16 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(p1 + x));
        __m256i p2_u16 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(p2 + x));

        __m128i p1_u16_lo = _mm256_castsi256_si128(p1_u16);
        __m128i p1_u16_hi = _mm256_extracti128_si256(p1_u16, 1);
        __m128i p2_u16_lo = _mm256_castsi256_si128(p2_u16);
        __m128i p2_u16_hi = _mm256_extracti128_si256(p2_u16, 1);

        // to int32 (8 pixels per __m256i)
        __m256i p1_i32_lo = _mm256_cvtepu16_epi32(p1_u16_lo);
        __m256i p1_i32_hi = _mm256_cvtepu16_epi32(p1_u16_hi);
        __m256i p2_i32_lo = _mm256_cvtepu16_epi32(p2_u16_lo);
        __m256i p2_i32_hi = _mm256_cvtepu16_epi32(p2_u16_hi);

        // p2-p1 in int32 and only then float
        __m256 delta_lo_32 = _mm256_sub_epi32(p2_i32_lo, p1_i32_lo);
        __m256 delta_hi_32 = _mm256_sub_epi32(p2_i32_hi, p1_i32_hi);
        __m256 delta_f_lo = _mm256_cvtepi32_ps(delta_lo_32);
        __m256 delta_f_hi = _mm256_cvtepi32_ps(delta_hi_32);

        __m256 delta_scaled_lo, delta_scaled_hi;
        if constexpr (has_mask) {
          __m256i mask_u16 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(mask + x));
          __m128i mask_u16_lo = _mm256_castsi256_si128(mask_u16);
          __m128i mask_u16_hi = _mm256_extracti128_si256(mask_u16, 1);

          __m256i mask_i32_lo = _mm256_cvtepu16_epi32(mask_u16_lo);
          __m256i mask_i32_hi = _mm256_cvtepu16_epi32(mask_u16_hi);
          __m256 mask_f_lo = _mm256_cvtepi32_ps(mask_i32_lo);
          __m256 mask_f_hi = _mm256_cvtepi32_ps(mask_i32_hi);

          __m256 new_factor_lo = _mm256_mul_ps(mask_f_lo, factor_v);
          __m256 new_factor_hi = _mm256_mul_ps(mask_f_hi, factor_v);

          delta_scaled_lo = _mm256_mul_ps(delta_f_lo, new_factor_lo);
          delta_scaled_hi = _mm256_mul_ps(delta_f_hi, new_factor_hi);
        }
        else {
          delta_scaled_lo = _mm256_mul_ps(delta_f_lo, factor_v);
          delta_scaled_hi = _mm256_mul_ps(delta_f_hi, factor_v);
        }

        // round and truncate to int32
        delta_scaled_lo = _mm256_add_ps(delta_scaled_lo, rounder);
        delta_scaled_hi = _mm256_add_ps(delta_scaled_hi, rounder);
        __m256i delta_scaled_i32_lo = _mm256_cvttps_epi32(delta_scaled_lo);
        __m256i delta_scaled_i32_hi = _mm256_cvttps_epi32(delta_scaled_hi);

        // p1 + delta_scaled still in int32
        __m256i result_i32_lo = _mm256_add_epi32(p1_i32_lo, delta_scaled_i32_lo);
        __m256i result_i32_hi = _mm256_add_epi32(p1_i32_hi, delta_scaled_i32_hi);

        // 8 int32 -> 8 uint16 in lower half
        __m256i i16_lo = _mm256_packus_epi32(result_i32_lo, result_i32_lo);
        __m256i i16_hi = _mm256_packus_epi32(result_i32_hi, result_i32_hi);

        // get correct layout
        i16_lo = _mm256_permute4x64_epi64(i16_lo, 0xD8);
        i16_hi = _mm256_permute4x64_epi64(i16_hi, 0xD8);

        // Extract lower 128 bits from each
        __m128i result_u16_lo = _mm256_castsi256_si128(i16_lo);
        __m128i result_u16_hi = _mm256_castsi256_si128(i16_hi);

        // Combine back into a single __m256i and store
        __m256i result_u16 = _mm256_setr_m128i(result_u16_lo, result_u16_hi);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(p1 + x), result_u16);
      }
    }

    // Scalar tail
    // Working with exact dimension, overlay is possible atop an existing frame at any position
    for (int x = wMod / sizeof(pixel_t); x < width; x++) {
      const float new_factor = has_mask ? static_cast<float>(reinterpret_cast<const pixel_t*>(mask)[x]) * factor : factor;
      const pixel_t pix1 = reinterpret_cast<pixel_t*>(p1)[x];
      const pixel_t pix2 = reinterpret_cast<const pixel_t*>(p2)[x];
      pixel_t result = pix1 + (int)((pix2 - pix1) * new_factor + 0.5f);
      reinterpret_cast<pixel_t*>(p1)[x] = result;
    }

    p1 += p1_pitch;
    p2 += p2_pitch;
    if (has_mask)
      mask += mask_pitch;
  }
}

// instantiate
// mask yes/no
template void overlay_blend_avx2_uint<true, uint8_t, true>(BYTE* p1, const BYTE* p2, const BYTE* mask,
  const int p1_pitch, const int p2_pitch, const int mask_pitch, const int width, const int height, const int opacity, const float opacity_f, const int bits_per_pixel);
template void overlay_blend_avx2_uint<true, uint16_t, true>(BYTE* p1, const BYTE* p2, const BYTE* mask,
  const int p1_pitch, const int p2_pitch, const int mask_pitch, const int width, const int height, const int opacity, const float opacity_f, const int bits_per_pixel);
template void overlay_blend_avx2_uint<true, uint16_t, false>(BYTE* p1, const BYTE* p2, const BYTE* mask,
  const int p1_pitch, const int p2_pitch, const int mask_pitch, const int width, const int height, const int opacity, const float opacity_f, const int bits_per_pixel);
//--
template void overlay_blend_avx2_uint<false, uint8_t, true>(BYTE* p1, const BYTE* p2, const BYTE* mask,
  const int p1_pitch, const int p2_pitch, const int mask_pitch, const int width, const int height, const int opacity, const float opacity_f, const int bits_per_pixel);
template void overlay_blend_avx2_uint<false, uint16_t, true>(BYTE* p1, const BYTE* p2, const BYTE* mask,
  const int p1_pitch, const int p2_pitch, const int mask_pitch, const int width, const int height, const int opacity, const float opacity_f, const int bits_per_pixel);
template void overlay_blend_avx2_uint<false, uint16_t, false>(BYTE* p1, const BYTE* p2, const BYTE* mask,
  const int p1_pitch, const int p2_pitch, const int mask_pitch, const int width, const int height, const int opacity, const float opacity_f, const int bits_per_pixel);

template<bool has_mask>
void overlay_blend_avx2_float(BYTE* p1, const BYTE* p2, const BYTE* mask,
  const int p1_pitch, const int p2_pitch, const int mask_pitch,
  const int width, const int height, const int opacity, const float opacity_f, const int bits_per_pixel)
{

  const int realwidth = width * sizeof(float);

  int wMod32 = (realwidth / 32) * 32;
  auto opacity_v = _mm256_set1_ps(opacity_f);

  for (int y = 0; y < height; y++) {
    for (int x = 0; x < wMod32; x += 32) {
      auto p1_f = _mm256_loadu_ps(reinterpret_cast<const float*>(p1 + x));
      auto p2_f = _mm256_loadu_ps(reinterpret_cast<const float*>(p2 + x));
      __m256 new_mask;
      if constexpr (has_mask) {
        new_mask = _mm256_loadu_ps(reinterpret_cast<const float*>(mask + x));
        new_mask = _mm256_mul_ps(new_mask, opacity_v);
      }
      else {
        new_mask = opacity_v;
      }
      auto result = _mm256_add_ps(p1_f, _mm256_mul_ps(_mm256_sub_ps(p2_f, p1_f), new_mask)); // p1*(1-mask) + p2*mask = p1+(p2-p1)*mask

      _mm256_storeu_ps(reinterpret_cast<float*>(p1 + x), result);
    }

    // Leftover value
    // Working with exact dimension, overlay is possible atop an existing frame at any position
    for (int x = wMod32 / sizeof(float); x < width; x++) {
      auto new_mask = has_mask ? reinterpret_cast<const float*>(mask)[x] * opacity_f : opacity_f;
      auto p1x = reinterpret_cast<float*>(p1)[x];
      auto p2x = reinterpret_cast<const float*>(p2)[x];
      auto result = p1x + (p2x - p1x) * new_mask; // p1x*(1-new_mask) + p2x*mask
      reinterpret_cast<float*>(p1)[x] = result;
    }


    p1 += p1_pitch;
    p2 += p2_pitch;
    if constexpr (has_mask)
      mask += mask_pitch;
  }
}

template void overlay_blend_avx2_float<false>(BYTE* p1, const BYTE* p2, const BYTE* mask,
  const int p1_pitch, const int p2_pitch, const int mask_pitch, const int width, const int height, const int opacity, const float opacity_f, const int bits_per_pixel);
template void overlay_blend_avx2_float<true>(BYTE* p1, const BYTE* p2, const BYTE* mask,
  const int p1_pitch, const int p2_pitch, const int mask_pitch, const int width, const int height, const int opacity, const float opacity_f, const int bits_per_pixel);

