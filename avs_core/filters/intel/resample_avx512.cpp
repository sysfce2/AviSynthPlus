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

#include <avs/config.h>
#include "../core/internal.h"

#include <avs/alignment.h>
#include <avs/minmax.h>

#include "check_avx512.h" // compiler avx512 directives check, basic f and bw is required
#include "resample_avx2.h"
#include "resample_avx512.h"

// Functions in this file are dispatched by AVS_CPUF_AVX512_FAST group feature flag.
// Assumes F, CD, BW, DQ, VL, VNNI, VBMI, VBMI2, BITALG, VPOPCNTDQ
// gcc/clang flags: -mavx512f -mavx512cd -mavx512bw -mavx512dq -mavx512vl -mavx512vnni -mavx512vbmi -mavx512vbmi2 -mavx512bitalg -mavx512vpopcntdq

//------- 512 bit float Horizontals

// Safe quad lane partial load with AVX512
// Read exactly N pixels (where N mod 4 is the template parameter), avoiding
// - reading beyond the end of the source buffer.
// - avoid NaN contamination by padding with zeros.
template <int Nmod4>
AVS_FORCEINLINE static __m512 _mm512_load_partial_safe_4_m128(const float* src_ptr_offsetted1, const float* src_ptr_offsetted2, const float* src_ptr_offsetted3, const float* src_ptr_offsetted4) {
  __m128 s1, s2, s3, s4;
  switch (Nmod4) {
  case 1:
    s1 = _mm_set_ps(0.0f, 0.0f, 0.0f, src_ptr_offsetted1[0]);
    s2 = _mm_set_ps(0.0f, 0.0f, 0.0f, src_ptr_offsetted2[0]);
    s3 = _mm_set_ps(0.0f, 0.0f, 0.0f, src_ptr_offsetted3[0]);
    s4 = _mm_set_ps(0.0f, 0.0f, 0.0f, src_ptr_offsetted4[0]);
    // ideally: movss
    break;
  case 2:
    s1 = _mm_set_ps(0.0f, 0.0f, src_ptr_offsetted1[1], src_ptr_offsetted1[0]);
    s2 = _mm_set_ps(0.0f, 0.0f, src_ptr_offsetted2[1], src_ptr_offsetted2[0]);
    s3 = _mm_set_ps(0.0f, 0.0f, src_ptr_offsetted3[1], src_ptr_offsetted3[0]);
    s4 = _mm_set_ps(0.0f, 0.0f, src_ptr_offsetted4[1], src_ptr_offsetted4[0]);
    // ideally: movsd
    break;
  case 3:
    s1 = _mm_set_ps(0.0f, src_ptr_offsetted1[2], src_ptr_offsetted1[1], src_ptr_offsetted1[0]);
    s2 = _mm_set_ps(0.0f, src_ptr_offsetted2[2], src_ptr_offsetted2[1], src_ptr_offsetted2[0]);
    s3 = _mm_set_ps(0.0f, src_ptr_offsetted3[2], src_ptr_offsetted3[1], src_ptr_offsetted3[0]);
    s4 = _mm_set_ps(0.0f, src_ptr_offsetted4[2], src_ptr_offsetted4[1], src_ptr_offsetted4[0]);
    // ideally: movss + movsd + shuffle or movsd + insert
    break;
  case 0:
    s1 = _mm_set_ps(src_ptr_offsetted1[3], src_ptr_offsetted1[2], src_ptr_offsetted1[1], src_ptr_offsetted1[0]);
    s2 = _mm_set_ps(src_ptr_offsetted2[3], src_ptr_offsetted2[2], src_ptr_offsetted2[1], src_ptr_offsetted2[0]);
    s3 = _mm_set_ps(src_ptr_offsetted3[3], src_ptr_offsetted3[2], src_ptr_offsetted3[1], src_ptr_offsetted3[0]);
    s4 = _mm_set_ps(src_ptr_offsetted4[3], src_ptr_offsetted4[2], src_ptr_offsetted4[1], src_ptr_offsetted4[0]);
    // ideally: movups
    break;
  default:
    s1 = _mm_setzero_ps(); // n/a cannot happen
    s2 = _mm_setzero_ps();
    s3 = _mm_setzero_ps();
    s4 = _mm_setzero_ps();
  }
  __m512 result = _mm512_castps128_ps512(s1); // Cast the first __m128 to __m512
  result = _mm512_insertf32x4(result, s2, 1); // Insert the second __m128 at position 1
  result = _mm512_insertf32x4(result, s3, 2); // Insert the third __m128 at position 2
  result = _mm512_insertf32x4(result, s4, 3); // Insert the fourth __m128 at position 3
  return result;
}

// same as resize_h_planar_float_avx2_gather_permutex_vstripe_ks4_pix16_check
// returns true if only transpose method is allowed, false if permutex method can be used
bool resize_h_planar_float_avx512_gather_permutex_vstripe_ks4_check(ResamplingProgram* program)
{
  // 'target_size_alignment' ensures we can safely access pixel_offset[] entries using offsets like
  // pixel_offset[x + 0] to pixel_offset[x + 15] per 16-pixel block processing
  assert(program->target_size_alignment >= 16);

  // Ensure that coefficient loading is safe for 4 float loads
  assert(program->filter_size_alignment >= 4);

  for (int x = 0; x < program->target_size; x += 16)
  {
    int start_off = program->pixel_offset[x + 0];
    // program->pixel_offset[x + 15] is still valid, since program->target_size_alignment >= 8
    // and pixel_offset[] values for x >= target_size are the same as for x=target_size-1.
    // This is ensured during the resampling program setup in resize_prepare_coeffs()
    int end_off = program->pixel_offset[x + 15] + 3;
    if ((end_off - start_off) >= 16) {
      return true; // only transpose is allowed
    }
  }
  return false; // permute is OK.
}

// Processes a horizontal resampling kernel of up to four coefficients for float pixel types.
// Supports BilinearResize, BicubicResize, or sinc with up to 2 taps (filter size <= 4).
// AVX512 optimization loads and processes four float coefficients and sixteen pixels simultaneously.
// The 'filtersizemod4' template parameter (0-3) helps optimize for different filter sizes modulo 4.
// This AVX512 requires only filter_size_alignment of 4.
template<int filtersizemod4>
void resize_h_planar_float_avx512_transpose_vstripe_ks4(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel) {
  assert(filtersizemod4 >= 0 && filtersizemod4 <= 3);

  const int filter_size = program->filter_size; // aligned, practically the coeff table stride

  src_pitch /= sizeof(float);
  dst_pitch /= sizeof(float);

  float* src = (float*)src8;
  float* dst = (float*)dst8;

  constexpr int PIXELS_AT_A_TIME = 16; // Process sixteen pixels in parallel using AVX512 (4x4 using m128 lanes)

  // 'source_overread_beyond_targetx' indicates if the filter kernel can read beyond the target width.
  // Even if the filter alignment allows larger reads, our safety boundary for unaligned loads starts at 4 pixels back
  // from the target width, as we load 4 floats at once conceptually with our safe load.
  const int width_safe_mod = (program->safelimit_4_pixels.overread_possible ? program->safelimit_4_pixels.source_overread_beyond_targetx : width) / PIXELS_AT_A_TIME * PIXELS_AT_A_TIME;

  // Preconditions:
  assert(program->filter_size_real <= 4);
  assert(program->target_size_alignment >= 16); // Must align for 16 pixel offsets
  assert(program->filter_size_alignment >= 4);
  assert(FRAME_ALIGN >= 64); // Adjusted for 16 pixels AviSynth+ default

  // Vertical stripe alignment
  constexpr int STRIPE_ALIGN = 16;

  const size_t cache_size_L2 = program->cache_size_L2;
  int max_scanlines = resampler_h_float_detect_optimal_scanline(program->source_size, program->target_size, cache_size_L2)
    / STRIPE_ALIGN * STRIPE_ALIGN;

  if (max_scanlines < STRIPE_ALIGN) max_scanlines = STRIPE_ALIGN;

  // --- outer loop: vertical stripes ---
  for (auto y_from = 0; y_from < height; y_from += max_scanlines) {
    size_t y_to = std::min(y_from + max_scanlines, height);

    // Reset current_coeff for the start of the stripe
    const float* AVS_RESTRICT current_coeff = program->pixel_coefficient_float;

    int x = 0;

    // Lambda for 512-bit Core
    auto do_h_float_core = [&](auto partial_load) {

      // load 4x4 sets of coefficients (16 pixels total)
      // at once before the height loop.
    // Coefficients for the source pixel offset (for src_ptr + begin1 [0..3], begin5 [0..3], begin9 [0..3], begin13 [0..3])
    __m512 coef_1_5_9_13 = _mm512_load_4_m128(current_coeff + filter_size * 0, current_coeff + filter_size * 4, current_coeff + filter_size * 8, current_coeff + filter_size * 12);
    __m512 coef_2_6_10_14 = _mm512_load_4_m128(current_coeff + filter_size * 1, current_coeff + filter_size * 5, current_coeff + filter_size * 9, current_coeff + filter_size * 13);
    __m512 coef_3_7_11_15 = _mm512_load_4_m128(current_coeff + filter_size * 2, current_coeff + filter_size * 6, current_coeff + filter_size * 10, current_coeff + filter_size * 14);
    __m512 coef_4_8_12_16 = _mm512_load_4_m128(current_coeff + filter_size * 3, current_coeff + filter_size * 7, current_coeff + filter_size * 11, current_coeff + filter_size * 15);

    _MM_TRANSPOSE16_LANE4_PS(coef_1_5_9_13, coef_2_6_10_14, coef_3_7_11_15, coef_4_8_12_16);

    // Pixel offsets for the current target x-positions.
    // Even for x >= width, these offsets are guaranteed to be within the allocated 'target_size_alignment'.
    const int begin1 = program->pixel_offset[x + 0];
    const int begin2 = program->pixel_offset[x + 1];
    const int begin3 = program->pixel_offset[x + 2];
    const int begin4 = program->pixel_offset[x + 3];
    const int begin5 = program->pixel_offset[x + 4];
    const int begin6 = program->pixel_offset[x + 5];
    const int begin7 = program->pixel_offset[x + 6];
    const int begin8 = program->pixel_offset[x + 7];
    const int begin9 = program->pixel_offset[x + 8];
    const int begin10 = program->pixel_offset[x + 9];
    const int begin11 = program->pixel_offset[x + 10];
    const int begin12 = program->pixel_offset[x + 11];
    const int begin13 = program->pixel_offset[x + 12];
    const int begin14 = program->pixel_offset[x + 13];
    const int begin15 = program->pixel_offset[x + 14];
    const int begin16 = program->pixel_offset[x + 15];

      int y = y_from;

      // Calculate pointers ONCE before the inner loop (Optimization from AVX2 version)
      float* AVS_RESTRICT dst_ptr = dst + y * dst_pitch + x;
      const float* src_ptr = src + y * src_pitch;

      // Inner loop: vertical processing. unroll 2 tested, no benefit
      for (; y < y_to; ++y) {

      __m512 data_1_5_9_13;
      __m512 data_2_6_10_14;
      __m512 data_3_7_11_15;
      __m512 data_4_8_12_16;

      if constexpr (partial_load) {
        // In the potentially unsafe zone (near the right edge of the image), we use a safe loading function
        // to prevent reading beyond the allocated source scanline.

        data_1_5_9_13 = _mm512_load_partial_safe_4_m128<filtersizemod4>(src_ptr + begin1, src_ptr + begin5, src_ptr + begin9, src_ptr + begin13);
        data_2_6_10_14 = _mm512_load_partial_safe_4_m128<filtersizemod4>(src_ptr + begin2, src_ptr + begin6, src_ptr + begin10, src_ptr + begin14);
        data_3_7_11_15 = _mm512_load_partial_safe_4_m128<filtersizemod4>(src_ptr + begin3, src_ptr + begin7, src_ptr + begin11, src_ptr + begin15);
        data_4_8_12_16 = _mm512_load_partial_safe_4_m128<filtersizemod4>(src_ptr + begin4, src_ptr + begin8, src_ptr + begin12, src_ptr + begin16);
      }
      else {
        // In the safe zone, we can directly load 4 pixels at a time for each of the four lanes.
        data_1_5_9_13 = _mm512_loadu_4_m128(src_ptr + begin1, src_ptr + begin5, src_ptr + begin9, src_ptr + begin13);
        data_2_6_10_14 = _mm512_loadu_4_m128(src_ptr + begin2, src_ptr + begin6, src_ptr + begin10, src_ptr + begin14);
        data_3_7_11_15 = _mm512_loadu_4_m128(src_ptr + begin3, src_ptr + begin7, src_ptr + begin11, src_ptr + begin15);
        data_4_8_12_16 = _mm512_loadu_4_m128(src_ptr + begin4, src_ptr + begin8, src_ptr + begin12, src_ptr + begin16);
      }

        // note: 256 bit simulation is slower
      _MM_TRANSPOSE16_LANE4_PS(data_1_5_9_13, data_2_6_10_14, data_3_7_11_15, data_4_8_12_16);

      __m512 result = _mm512_mul_ps(data_1_5_9_13, coef_1_5_9_13);
      result = _mm512_fmadd_ps(data_2_6_10_14, coef_2_6_10_14, result);
      result = _mm512_fmadd_ps(data_3_7_11_15, coef_3_7_11_15, result);
      result = _mm512_fmadd_ps(data_4_8_12_16, coef_4_8_12_16, result);

        _mm512_stream_ps(dst_ptr, result);

      dst_ptr += dst_pitch;
      src_ptr += src_pitch;
    } // y

      // Move to the next set of coefficients for the next 16 output pixels
      current_coeff += filter_size * 16;
      };

    // Process the 'safe zone' where direct full unaligned loads are acceptable.
    for (; x < width_safe_mod; x += PIXELS_AT_A_TIME)
    {
      do_h_float_core(std::false_type{});  // partial_load == false, use direct _mm512_loadu_ps
    }

    // Process the potentially 'unsafe zone' near the image edge, using safe loading.
    for (; x < width; x += PIXELS_AT_A_TIME)
    {
      do_h_float_core(std::true_type{}); // partial_load == true, use the safer _mm512_load_partial_safe_4_m128
    }
  }
}

//instatiate
template void resize_h_planar_float_avx512_transpose_vstripe_ks4<0>(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel);
template void resize_h_planar_float_avx512_transpose_vstripe_ks4<1>(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel);
template void resize_h_planar_float_avx512_transpose_vstripe_ks4<2>(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel);
template void resize_h_planar_float_avx512_transpose_vstripe_ks4<3>(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel);

// FIXME: make it safe + correct, like the avx2 counterpart
void resize_h_planar_float_avx512_permutex_vstripe_ks4(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel)
{
  const int filter_size = program->filter_size; // aligned, practically the coeff table stride

  src_pitch /= sizeof(float);
  dst_pitch /= sizeof(float);

  float* src = (float*)src8;
  float* dst = (float*)dst8;

  const float* AVS_RESTRICT current_coeff = (const float* AVS_RESTRICT)program->pixel_coefficient_float;

  constexpr int PIXELS_AT_A_TIME = 16; // Process sixteen pixels in parallel using AVX512 (4x4 using m128 lanes)

  // 'source_overread_beyond_targetx' indicates if the filter kernel can read beyond the target width.
  // Even if the filter alignment allows larger reads, our safety boundary for unaligned loads starts at 4 pixels back
  // from the target width, as we load 4 floats at once conceptually with our safe load.
  const int width_safe_mod = (program->safelimit_4_pixels.overread_possible ? program->safelimit_4_pixels.source_overread_beyond_targetx : width) / PIXELS_AT_A_TIME * PIXELS_AT_A_TIME;

  // Preconditions:
  assert(program->filter_size_real <= 4); // We preload all relevant coefficients (up to 4) before the height loop.

  // 'target_size_alignment' ensures we can safely access coefficients using offsets like
  // 'filter_size * 7' when processing 8 H pixels at a time or
  // 'filter_size * 15' when processing 16 H pixels at a time
  assert(program->target_size_alignment >= 16); // Adjusted for 16 pixels
  assert(FRAME_ALIGN >= 64); // Adjusted for 16 pixels AviSynth+ default

  // Ensure that coefficient loading beyond the valid target size is safe for 4x4 float loads.
  assert(program->filter_size_alignment >= 4);

  int x = 0;

  for (x = 0; x < width; x += 16)
  {
    // prepare coefs in transposed V-form
    __m512 coef_r0 = _mm512_load_4_m128(current_coeff + filter_size * 0, current_coeff + filter_size * 4, current_coeff + filter_size * 8, current_coeff + filter_size * 12);
    __m512 coef_r1 = _mm512_load_4_m128(current_coeff + filter_size * 1, current_coeff + filter_size * 5, current_coeff + filter_size * 9, current_coeff + filter_size * 13);
    __m512 coef_r2 = _mm512_load_4_m128(current_coeff + filter_size * 2, current_coeff + filter_size * 6, current_coeff + filter_size * 10, current_coeff + filter_size * 14);
    __m512 coef_r3 = _mm512_load_4_m128(current_coeff + filter_size * 3, current_coeff + filter_size * 7, current_coeff + filter_size * 11, current_coeff + filter_size * 15);

    _MM_TRANSPOSE16_LANE4_PS(coef_r0, coef_r1, coef_r2, coef_r3);

    // convert resampling program in H-form into permuting indexes for src transposition in V-form
    int iStart = program->pixel_offset[x + 0];

    __m512i perm_0 = _mm512_set_epi32(
      program->pixel_offset[x + 15] - iStart,
      program->pixel_offset[x + 14] - iStart,
      program->pixel_offset[x + 13] - iStart,
      program->pixel_offset[x + 12] - iStart,
      program->pixel_offset[x + 11] - iStart,
      program->pixel_offset[x + 10] - iStart,
      program->pixel_offset[x + 9] - iStart,
      program->pixel_offset[x + 8] - iStart,
      program->pixel_offset[x + 7] - iStart,
      program->pixel_offset[x + 6] - iStart,
      program->pixel_offset[x + 5] - iStart,
      program->pixel_offset[x + 4] - iStart,
      program->pixel_offset[x + 3] - iStart,
      program->pixel_offset[x + 2] - iStart,
      program->pixel_offset[x + 1] - iStart,
      0);

    __m512i one_epi32 = _mm512_set1_epi32(1);
    __m512i perm_1 = _mm512_add_epi32(perm_0, one_epi32);
    one_epi32 = _mm512_set1_epi32(program->pixel_offset[x + 2] - program->pixel_offset[x + 1]);
    __m512i perm_2 = _mm512_add_epi32(perm_1, one_epi32);
    one_epi32 = _mm512_set1_epi32(program->pixel_offset[x + 3] - program->pixel_offset[x + 2]);
    __m512i perm_3 = _mm512_add_epi32(perm_2, one_epi32);

    float* AVS_RESTRICT dst_ptr = dst + x;
    const float* src_ptr = src + program->pixel_offset[x + 0]; // all permute offsets relative to this start offset

    for (int y = 0; y < height; y++) // single row proc
    {
      __m512 data_src = _mm512_loadu_ps(src_ptr);
      __m512 data_src2 = _mm512_loadu_ps(src_ptr + 16); // not always needed for upscale also can cause end of buffer overread - need to add limitation (special end of buffer processing ?)

      __m512 data_0 = _mm512_permutex2var_ps(data_src, perm_0, data_src2);
      __m512 data_1 = _mm512_permutex2var_ps(data_src, perm_1, data_src2);
      __m512 data_2 = _mm512_permutex2var_ps(data_src, perm_2, data_src2);
      __m512 data_3 = _mm512_permutex2var_ps(data_src, perm_3, data_src2);

      __m512 result0 = _mm512_mul_ps(data_0, coef_r0);
      __m512 result1 = _mm512_mul_ps(data_2, coef_r2);

      result0 = _mm512_fmadd_ps(data_1, coef_r1, result0);
      result1 = _mm512_fmadd_ps(data_3, coef_r3, result1);

      _mm512_store_ps(dst_ptr, _mm512_add_ps(result0, result1));

      dst_ptr += dst_pitch;
      src_ptr += src_pitch;
    }

    current_coeff += filter_size * 16;
  }
}


//-------- 512 bit float Verticals

// base version, no horizontal unrolling
void resize_v_avx512_planar_float(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int target_height, int bits_per_pixel)
{
  AVS_UNUSED(bits_per_pixel);

  const int filter_size = program->filter_size;
  const float* AVS_RESTRICT current_coeff = program->pixel_coefficient_float;

  const float* src = (const float*)src8;
  float* AVS_RESTRICT dst = (float*)dst8;
  dst_pitch = dst_pitch / sizeof(float);
  src_pitch = src_pitch / sizeof(float);

  const int kernel_size = program->filter_size_real; // not the aligned
  const int kernel_size_mod2 = (kernel_size / 2) * 2; // Process pairs of rows for better efficiency
  const bool notMod2 = kernel_size_mod2 < kernel_size;

  for (int y = 0; y < target_height; y++) {
    int offset = program->pixel_offset[y];
    const float* src_ptr = src + offset * src_pitch;

    // 64 byte 16 floats (AVX512 register holds 16 floats)
    // no need for wmod8, alignment is safe 32 bytes at least - is it safe for 64 bytes ?
    for (int x = 0; x < width; x += 16) {
      __m512 result_single = _mm512_setzero_ps();
      __m512 result_single_2 = _mm512_setzero_ps();

      const float* AVS_RESTRICT src2_ptr = src_ptr + x; // __restrict here

      // Process pairs of rows for better efficiency (2 coeffs/cycle)
      // two result variables for potential parallel operation
      int i = 0;
      for (; i < kernel_size_mod2; i += 2) {
        __m512 coeff_even = _mm512_set1_ps(current_coeff[i]);
        __m512 coeff_odd = _mm512_set1_ps(current_coeff[i + 1]);

        __m512 src_even = _mm512_loadu_ps(src2_ptr);
        __m512 src_odd = _mm512_loadu_ps(src2_ptr + src_pitch);

        result_single = _mm512_fmadd_ps(src_even, coeff_even, result_single);
        result_single_2 = _mm512_fmadd_ps(src_odd, coeff_odd, result_single_2);

        src2_ptr += 2 * src_pitch;
      }

      result_single = _mm512_add_ps(result_single, result_single_2);

      // Process the last odd row if needed
      if (notMod2) {
        __m512 coeff = _mm512_set1_ps(current_coeff[i]);
        __m512 src_val = _mm512_loadu_ps(src2_ptr);
        result_single = _mm512_fmadd_ps(src_val, coeff, result_single);
      }

      _mm512_store_ps(dst + x, result_single);
    }

    dst += dst_pitch;
    current_coeff += filter_size;
  }
}

// memory-optimized version of resize_v_avx512_planar_float
void resize_v_avx512_planar_float_w_sr(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int target_height, int bits_per_pixel)
{
  AVS_UNUSED(bits_per_pixel);

  const int filter_size = program->filter_size;
  const float* AVS_RESTRICT current_coeff = program->pixel_coefficient_float;

  const float* src = (const float*)src8;
  float* AVS_RESTRICT dst = (float*)dst8;
  dst_pitch = dst_pitch / sizeof(float);
  src_pitch = src_pitch / sizeof(float);

  const int kernel_size = program->filter_size_real; // not the aligned
  const int kernel_size_mod2 = (kernel_size / 2) * 2; // Process pairs of rows for better efficiency
  const bool notMod2 = kernel_size_mod2 < kernel_size;

  for (int y = 0; y < target_height; y++) {
    int offset = program->pixel_offset[y];
    const float* src_ptr = src + offset * src_pitch;

    int x = 0;

    // 128-64-32 pixels prelude
    // we spare some coeff reload per unrolled loop
    // Perhaps for an i7-11700 which has only 1x512 as 2x256 fma unit, not ideal.

    // Process by 8x 512 (8 x 16 floats) to make memory read/write linear streams longer
    // 32x512 bit registers should be enough
    const int width_mod128 = (width / 128) * 128;
    for (; x < width_mod128; x += 128) {
      __m512 result_1 = _mm512_setzero_ps();
      __m512 result_2 = _mm512_setzero_ps();
      __m512 result_3 = _mm512_setzero_ps();
      __m512 result_4 = _mm512_setzero_ps();
      __m512 result_5 = _mm512_setzero_ps();
      __m512 result_6 = _mm512_setzero_ps();
      __m512 result_7 = _mm512_setzero_ps();
      __m512 result_8 = _mm512_setzero_ps();

      const float* AVS_RESTRICT src2_ptr = src_ptr + x;
      int i = 0;
      for (; i < kernel_size; i++) {
        __m512 coeff = _mm512_set1_ps(current_coeff[i]);

        __m512 src_1 = _mm512_load_ps(src2_ptr);
        __m512 src_2 = _mm512_load_ps(src2_ptr + 16);
        __m512 src_3 = _mm512_load_ps(src2_ptr + 32);
        __m512 src_4 = _mm512_load_ps(src2_ptr + 48);
        __m512 src_5 = _mm512_load_ps(src2_ptr + 64);
        __m512 src_6 = _mm512_load_ps(src2_ptr + 80);
        __m512 src_7 = _mm512_load_ps(src2_ptr + 96);
        __m512 src_8 = _mm512_load_ps(src2_ptr + 112);
        result_1 = _mm512_fmadd_ps(src_1, coeff, result_1);
        result_2 = _mm512_fmadd_ps(src_2, coeff, result_2);
        result_3 = _mm512_fmadd_ps(src_3, coeff, result_3);
        result_4 = _mm512_fmadd_ps(src_4, coeff, result_4);
        result_5 = _mm512_fmadd_ps(src_5, coeff, result_5);
        result_6 = _mm512_fmadd_ps(src_6, coeff, result_6);
        result_7 = _mm512_fmadd_ps(src_7, coeff, result_7);
        result_8 = _mm512_fmadd_ps(src_8, coeff, result_8);

        src2_ptr += src_pitch;
      }

      _mm512_store_ps(dst + x, result_1);
      _mm512_store_ps(dst + x + 16, result_2);
      _mm512_store_ps(dst + x + 32, result_3);
      _mm512_store_ps(dst + x + 48, result_4);
      _mm512_store_ps(dst + x + 64, result_5);
      _mm512_store_ps(dst + x + 80, result_6);
      _mm512_store_ps(dst + x + 96, result_7);
      _mm512_store_ps(dst + x + 112, result_8);
    }

    // Process by 4x512 (4 x 16 floats) to make memory read/write linear streams longer
    const int width_mod64 = (width / 64) * 64;
    for (; x < width_mod64; x += 64) {
      __m512 result_1 = _mm512_setzero_ps();
      __m512 result_2 = _mm512_setzero_ps();
      __m512 result_3 = _mm512_setzero_ps();
      __m512 result_4 = _mm512_setzero_ps();

      const float* AVS_RESTRICT src2_ptr = src_ptr + x;

      int i = 0;
      for (; i < kernel_size; i++) {
        __m512 coeff = _mm512_set1_ps(current_coeff[i]);

        __m512 src_1 = _mm512_load_ps(src2_ptr);
        __m512 src_2 = _mm512_load_ps(src2_ptr + 16);
        __m512 src_3 = _mm512_load_ps(src2_ptr + 32);
        __m512 src_4 = _mm512_load_ps(src2_ptr + 48);

        result_1 = _mm512_fmadd_ps(src_1, coeff, result_1);
        result_2 = _mm512_fmadd_ps(src_2, coeff, result_2);
        result_3 = _mm512_fmadd_ps(src_3, coeff, result_3);
        result_4 = _mm512_fmadd_ps(src_4, coeff, result_4);

        src2_ptr += src_pitch;
      }

      _mm512_store_ps(dst + x, result_1);
      _mm512_store_ps(dst + x + 16, result_2);
      _mm512_store_ps(dst + x + 32, result_3);
      _mm512_store_ps(dst + x + 48, result_4);
    }

    // Process by 2x512 (2 x 16 floats) to make memory read/write linear streams longer,
    const int width_mod32 = (width / 32) * 32;
    for (; x < width_mod32; x += 32) {
      __m512 result_1 = _mm512_setzero_ps();
      __m512 result_2 = _mm512_setzero_ps();

      const float* AVS_RESTRICT src2_ptr = src_ptr + x;

      int i = 0;
      for (; i < kernel_size; i++) {
        __m512 coeff = _mm512_set1_ps(current_coeff[i]);

        __m512 src_1 = _mm512_load_ps(src2_ptr);
        __m512 src_2 = _mm512_load_ps(src2_ptr + 16);

        result_1 = _mm512_fmadd_ps(src_1, coeff, result_1);
        result_2 = _mm512_fmadd_ps(src_2, coeff, result_2);

        src2_ptr += src_pitch;
      }

      _mm512_store_ps(dst + x, result_1);
      _mm512_store_ps(dst + x + 16, result_2);
    }

    // Process 1x512 dual
    // 64 byte 16 floats (AVX512 register holds 16 floats)
    // row alignment is 64 bytes - so it is safe to load mod16 of float32.
    for (; x < width; x += 16) {
      __m512 result_single = _mm512_setzero_ps();
      __m512 result_single_2 = _mm512_setzero_ps();

      const float* AVS_RESTRICT src2_ptr = src_ptr + x;

      // Process pairs of rows for better efficiency (2 coeffs/cycle)
      // two result variables for potential parallel operation
      int i = 0;
      for (; i < kernel_size_mod2; i += 2) {
        __m512 coeff_even = _mm512_set1_ps(current_coeff[i]);
        __m512 coeff_odd = _mm512_set1_ps(current_coeff[i + 1]);

        __m512 src_even = _mm512_load_ps(src2_ptr);
        __m512 src_odd = _mm512_load_ps(src2_ptr + src_pitch);

        result_single = _mm512_fmadd_ps(src_even, coeff_even, result_single);
        result_single_2 = _mm512_fmadd_ps(src_odd, coeff_odd, result_single_2);

        src2_ptr += 2 * src_pitch;
      }

      result_single = _mm512_add_ps(result_single, result_single_2);

      // Process the last odd row if needed
      if (notMod2) {
        __m512 coeff = _mm512_set1_ps(current_coeff[i]);
        __m512 src_val = _mm512_load_ps(src2_ptr);
        result_single = _mm512_fmadd_ps(src_val, coeff, result_single);
      }

      _mm512_store_ps(dst + x, result_single);
    }


    dst += dst_pitch;
    current_coeff += filter_size;
  }
}

//----------------------- generic horizontal avx512 float

// AVX512 Horizontal float

// Three helpers, each for processing 4 target pixels from 16, 8 and 4 source pixel/coeff pairs.

// Helper, _mm256_zextps128_ps256 exists only in AVX512 VL
// zero-extend 128-bit float vector to 256-bit float vector
AVS_FORCEINLINE static __m256 _mm256_zextps128_ps256_simul_avx(__m128 a)
{
  // Flags defines by MSVC at /AVX512 mode
  // other flags: __AVX512F__, __AVX512CD__, __AVX512VL__, __AVX512BW__, __AVX512DQ__
#ifdef __AVX512VL__
  return _mm256_zextps128_ps256(a);
#else
  __m256 zero_v = _mm256_setzero_ps();
  return _mm256_insertf128_ps(zero_v, a, 0);
#endif
}

// 4 target pixels, each from 16 source pixel/coeff pair
// Called only when accessing 16 source pixels and coefficients at a time is safe
AVS_FORCEINLINE static void process_pix4_coeff16_h_float_core_512(
  const float* src,
  int begin1, int begin2, int begin3, int begin4,
  const float* current_coeff,
  int filter_size,
  __m512& result1, __m512& result2, __m512& result3, __m512& result4)
{
  // 16 source floats for each of the four beginning source offsets
  __m512 data_1 = _mm512_loadu_ps(src + begin1);
  __m512 data_2 = _mm512_loadu_ps(src + begin2);
  __m512 data_3 = _mm512_loadu_ps(src + begin3);
  __m512 data_4 = _mm512_loadu_ps(src + begin4);

  // 16 coefficients for each of the four output pixels
  __m512 coeff_1 = _mm512_loadu_ps(current_coeff);               // 16 coeffs for pixel 1
  __m512 coeff_2 = _mm512_loadu_ps(current_coeff + 1 * filter_size); // 16 coeffs for pixel 2
  __m512 coeff_3 = _mm512_loadu_ps(current_coeff + 2 * filter_size); // 16 coeffs for pixel 3
  __m512 coeff_4 = _mm512_loadu_ps(current_coeff + 3 * filter_size); // 16 coeffs for pixel 4

  // multiply and accumulate
  result1 = _mm512_fmadd_ps(data_1, coeff_1, result1);
  result2 = _mm512_fmadd_ps(data_2, coeff_2, result2);
  result3 = _mm512_fmadd_ps(data_3, coeff_3, result3);
  result4 = _mm512_fmadd_ps(data_4, coeff_4, result4);
}

// 4 target pixels, each from 8 source pixel/coeff pair
// Called only when accessing 8 source pixels and coefficients at a time is safe
AVS_FORCEINLINE static void process_pix4_coeff8_h_float_core(
  const float* src,
  int begin1, int begin2, int begin3, int begin4,
  const float* current_coeff,
  int filter_size,
  __m256& result1, __m256& result2, __m256& result3, __m256& result4)
{
  // Load 8 source floats for each of the four beginning source offsets
  // Load 8 coefficients for each of the four output pixels
  __m256 data_1 = _mm256_loadu_ps(src + begin1);
  __m256 coeff_1 = _mm256_load_ps(current_coeff);                    // 8 coeffs for pixel 1
  result1 = _mm256_fmadd_ps(data_1, coeff_1, result1);

  __m256 data_2 = _mm256_loadu_ps(src + begin2);
  __m256 coeff_2 = _mm256_load_ps(current_coeff + 1 * filter_size); // 8 coeffs for pixel 2
  result2 = _mm256_fmadd_ps(data_2, coeff_2, result2);

  __m256 data_3 = _mm256_loadu_ps(src + begin3);
  __m256 coeff_3 = _mm256_load_ps(current_coeff + 2 * filter_size); // 8 coeffs for pixel 3
  result3 = _mm256_fmadd_ps(data_3, coeff_3, result3);

  __m256 data_4 = _mm256_loadu_ps(src + begin4);
  __m256 coeff_4 = _mm256_load_ps(current_coeff + 3 * filter_size); // 8 coeffs for pixel 4
  result4 = _mm256_fmadd_ps(data_4, coeff_4, result4);
}

// 4 target pixels, each from 4 source pixel/coeff pair.
// Called only for first iteration when results are not initialized.
// Otherwise same as process_pix4_coeff8_h_float_core.
AVS_FORCEINLINE static void process_pix4_coeff4_h_float_core_first(
  const float* src,
  int begin1, int begin2, int begin3, int begin4,
  const float* current_coeff,
  int filter_size,
  __m256& result1, __m256& result2, __m256& result3, __m256& result4)
{
  // Pixel 1: Load, Multiply, and Zero-Extend to __m256
  __m128 data_1 = _mm_loadu_ps(src + begin1);
  __m128 coeff_1 = _mm_load_ps(current_coeff);
  __m128 mul_result1 = _mm_mul_ps(data_1, coeff_1);
  result1 = _mm256_zextps128_ps256(mul_result1); // Sets upper 128 bits to zero

  // Pixel 2: Load, Multiply, and Zero-Extend to __m256
  __m128 data_2 = _mm_loadu_ps(src + begin2);
  __m128 coeff_2 = _mm_load_ps(current_coeff + 1 * filter_size);
  __m128 mul_result2 = _mm_mul_ps(data_2, coeff_2);
  result2 = _mm256_zextps128_ps256(mul_result2); // Sets upper 128 bits to zero

  // Pixel 3: Load, Multiply, and Zero-Extend to __m256
  __m128 data_3 = _mm_loadu_ps(src + begin3);
  __m128 coeff_3 = _mm_load_ps(current_coeff + 2 * filter_size);
  __m128 mul_result3 = _mm_mul_ps(data_3, coeff_3);
  result3 = _mm256_zextps128_ps256(mul_result3); // Sets upper 128 bits to zero

  // Pixel 4: Load, Multiply, and Zero-Extend to __m256
  __m128 data_4 = _mm_loadu_ps(src + begin4);
  __m128 coeff_4 = _mm_load_ps(current_coeff + 3 * filter_size);
  __m128 mul_result4 = _mm_mul_ps(data_4, coeff_4);
  result4 = _mm256_zextps128_ps256(mul_result4); // Sets upper 128 bits to zero
}

// filtersize_hint: special: 0..4 for 4,8,16,24,32. Generic: -1
// filter_size is an aligned value and always multiple of 8 (prerequisite)
// Processing rules:
// if filtersize_hint==0: filter size <=4, do one coeff4 step only
// if filtersize_hint>=2: do 1 or 2 coeff16 steps
// if filtersize_hint==1 or 3: do 1 coeff8 step (0*16 or 1*16 step done already)
// if filtersize_hint==-1: unknown filter size, do 16,8 steps as possible
template<bool safe_aligned_mode, int filtersize_hint>
AVS_FORCEINLINE static void process_four_pixels_h_float_pix4of16_ks_4_8_16(
  const float* src_ptr,
  int begin1, int begin2, int begin3, int begin4,
  float* current_coeff,
  int filter_size,
  __m256& result1, __m256& result2, __m256& result3, __m256& result4,
  int kernel_size)
{

  // very special case: filter size <= 4
  if constexpr (safe_aligned_mode) {
    if (filtersize_hint == 0) {
      // Process 4 target pixels and 4 source pixels/coefficients at a time
      // XMM-based loop internally, but returns __m256 with upper 128 cleared
      // Do not assume initialized zeros in result1..4, they will be set here.
      process_pix4_coeff4_h_float_core_first(
        src_ptr + 0, begin1, begin2, begin3, begin4,
        current_coeff + 0,
        filter_size,
        result1, result2, result3, result4);
      return;
    }
  }

  int i = 0;

  // do by 16 coeffs until possible
  if (filtersize_hint == -1 || filtersize_hint >= 2) {
    __m512 result1_512 = _mm512_setzero_ps();
    __m512 result2_512 = _mm512_setzero_ps();
    __m512 result3_512 = _mm512_setzero_ps();
    __m512 result4_512 = _mm512_setzero_ps();
    const int ksmod16 = safe_aligned_mode ? (filter_size / 16 * 16) : (kernel_size / 16 * 16);
    // Process 4 target pixels and 16 source pixels/coefficients at a time (ZMM-based loop)
    for (; i < ksmod16; i += 16) {
      process_pix4_coeff16_h_float_core_512(
        src_ptr + i, begin1, begin2, begin3, begin4,
        current_coeff + i,
        filter_size,
        result1_512, result2_512, result3_512, result4_512);
    }
    // Horizontal sum reduction from __m512 to __m256
    result1 = _mm256_add_ps(_mm512_castps512_ps256(result1_512), _mm512_extractf32x8_ps(result1_512, 1));
    result2 = _mm256_add_ps(_mm512_castps512_ps256(result2_512), _mm512_extractf32x8_ps(result2_512, 1));
    result3 = _mm256_add_ps(_mm512_castps512_ps256(result3_512), _mm512_extractf32x8_ps(result3_512, 1));
    result4 = _mm256_add_ps(_mm512_castps512_ps256(result4_512), _mm512_extractf32x8_ps(result4_512, 1));
  }

  // filter sizes 16 or 32 can return here
  if constexpr (safe_aligned_mode && (filtersize_hint == 2 || filtersize_hint == 4)) {
    return;
  }

  if constexpr (!safe_aligned_mode) {
    if (i == kernel_size) return; // kernel_size is not known compile time
  }

  // When to do the coeff8 step:
  // not safe-aligned mode: always. E.g. kernel_size == 28 -> 16 done, now 10 rest, do 8 next
  // filtersize_hint == -1: not-compile-time known filtersize (kernel_size / 16 * 16 done, rest follows)
  // filtersize_hint == 1 or 3: 0*16 or 1*16 done, now do 1*8
  if (!safe_aligned_mode || filtersize_hint == -1 || filtersize_hint == 1 || filtersize_hint == 3) {
    // 32 bytes contain 8 floats. We will use 256-bit registers (YMM).
    const int ksmod8 = safe_aligned_mode ? (filter_size / 8 * 8) : (kernel_size / 8 * 8);

    // Process 4 target pixels and 8 source pixels/coefficients at a time (YMM-based loop)
    for (; i < ksmod8; i += 8) {
      process_pix4_coeff8_h_float_core(
        src_ptr + i, begin1, begin2, begin3, begin4,
        current_coeff + i,
        filter_size,
        result1, result2, result3, result4);
    }
  }

  if constexpr (!safe_aligned_mode) {
    // Right edge case.
    // Coeffs are zero padded, reading them is no problem.
    // But if we read past the end of source then we can get possible NaN contamination.
    // Handle the remainder: 1 to 7 source/coefficient elements.
    // real_kernel_size is used here, it's guaranteed that reading real_kernel_size elements
    // from any pixel_offset[] is safe and ends within the source buffer.
    // Optional 4-2-1 processing loop.

    if (i == kernel_size) return;

    // --- Define Base Pointers for Source and Coefficients ---
    const float* src_ptr1 = src_ptr + begin1;
    const float* src_ptr2 = src_ptr + begin2;
    const float* src_ptr3 = src_ptr + begin3;
    const float* src_ptr4 = src_ptr + begin4;

    float* current_coeff2 = current_coeff + 1 * filter_size;
    float* current_coeff3 = current_coeff + 2 * filter_size;
    float* current_coeff4 = current_coeff + 3 * filter_size;

    const int ksmod4 = kernel_size / 4 * 4;

    // -------------------------------------------------------------------
    // Mod 4 Block (4 elements for four pixels using __m128)
    // -------------------------------------------------------------------
    if (i < ksmod4) {
      // Load 4 source floats and 4 coefficients for each of the four output pixels
      __m128 data_1 = _mm_loadu_ps(src_ptr1 + i);
      __m128 coeff_1 = _mm_loadu_ps(current_coeff + i);
      __m128 temp_result1 = _mm_mul_ps(data_1, coeff_1);

      __m128 data_2 = _mm_loadu_ps(src_ptr2 + i);
      __m128 coeff_2 = _mm_loadu_ps(current_coeff2 + i);
      __m128 temp_result2 = _mm_mul_ps(data_2, coeff_2);

      __m128 data_3 = _mm_loadu_ps(src_ptr3 + i);
      __m128 coeff_3 = _mm_loadu_ps(current_coeff3 + i);
      __m128 temp_result3 = _mm_mul_ps(data_3, coeff_3);

      __m128 data_4 = _mm_loadu_ps(src_ptr4 + i);
      __m128 coeff_4 = _mm_loadu_ps(current_coeff4 + i);
      __m128 temp_result4 = _mm_mul_ps(data_4, coeff_4);

      // --- Accumulate 128-bit results into 256-bit registers ---
      // Note: Since we are using __m256, we must zero the high 128-bits before insertion/addition.

      result1 = _mm256_add_ps(result1, _mm256_zextps128_ps256(temp_result1));
      result2 = _mm256_add_ps(result2, _mm256_zextps128_ps256(temp_result2));
      result3 = _mm256_add_ps(result3, _mm256_zextps128_ps256(temp_result3));
      result4 = _mm256_add_ps(result4, _mm256_zextps128_ps256(temp_result4));

      i += 4;
      if (i == kernel_size) return;
    }

    const int ksmod2 = kernel_size / 2 * 2;

    // -------------------------------------------------------------------
    // New Mod 2 Block (2 elements for four pixels using __m128)
    // -------------------------------------------------------------------
    if (i < ksmod2) {
      // We only need to load 2 elements (4 floats) for the __m128 load, 
      // but the low 2 elements of the __m128 register are used.
      // Since we use the scalar accumulation method, we load 4, but only the 
      // first 2 elements will hold non-zero data (or load 2, and rely on 
      // the two __m128 registers to contain the result).

      // Let's stick to using the low 2 elements of __m128 for 2 elements.

      // Load 2 source floats and 2 coefficients for each of the four output pixels
      __m128 data_1 = _mm_castpd_ps(_mm_load_sd(reinterpret_cast<const double*>(src_ptr1 + i))); // Load 2 floats (double)
      __m128 coeff_1 = _mm_castpd_ps(_mm_load_sd(reinterpret_cast<const double*>(current_coeff + i)));
      __m128 temp_result1 = _mm_mul_ps(data_1, coeff_1);

      __m128 data_2 = _mm_castpd_ps(_mm_load_sd(reinterpret_cast<const double*>(src_ptr2 + i)));
      __m128 coeff_2 = _mm_castpd_ps(_mm_load_sd(reinterpret_cast<const double*>(current_coeff2 + i)));
      __m128 temp_result2 = _mm_mul_ps(data_2, coeff_2);

      __m128 data_3 = _mm_castpd_ps(_mm_load_sd(reinterpret_cast<const double*>(src_ptr3 + i)));
      __m128 coeff_3 = _mm_castpd_ps(_mm_load_sd(reinterpret_cast<const double*>(current_coeff3 + i)));
      __m128 temp_result3 = _mm_mul_ps(data_3, coeff_3);

      __m128 data_4 = _mm_castpd_ps(_mm_load_sd(reinterpret_cast<const double*>(src_ptr4 + i)));
      __m128 coeff_4 = _mm_castpd_ps(_mm_load_sd(reinterpret_cast<const double*>(current_coeff4 + i)));
      __m128 temp_result4 = _mm_mul_ps(data_4, coeff_4);

      result1 = _mm256_add_ps(result1, _mm256_zextps128_ps256(temp_result1));
      result2 = _mm256_add_ps(result2, _mm256_zextps128_ps256(temp_result2));
      result3 = _mm256_add_ps(result3, _mm256_zextps128_ps256(temp_result3));
      result4 = _mm256_add_ps(result4, _mm256_zextps128_ps256(temp_result4));

      i += 2;
      if (i == kernel_size) return;
    }

    // -------------------------------------------------------------------
    // Fallback Scalar Operation (1 element remaining)
    // -------------------------------------------------------------------
    if (i < kernel_size) {

      // Optimized scalar loop for the single remaining element
      float final_scalar1 = src_ptr1[i] * current_coeff[i];
      float final_scalar2 = src_ptr2[i] * current_coeff2[i];
      float final_scalar3 = src_ptr3[i] * current_coeff3[i];
      float final_scalar4 = src_ptr4[i] * current_coeff4[i];

      __m128 s1_128 = _mm_set_ss(final_scalar1);
      __m128 s2_128 = _mm_set_ss(final_scalar2);
      __m128 s3_128 = _mm_set_ss(final_scalar3);
      __m128 s4_128 = _mm_set_ss(final_scalar4);

      result1 = _mm256_add_ps(result1, _mm256_zextps128_ps256(s1_128));
      result2 = _mm256_add_ps(result2, _mm256_zextps128_ps256(s2_128));
      result3 = _mm256_add_ps(result3, _mm256_zextps128_ps256(s3_128));
      result4 = _mm256_add_ps(result4, _mm256_zextps128_ps256(s4_128));

      // i is now equal to kernel_size (i++)
    }
  }
}


template<bool is_safe, int filtersize_hint>
AVS_FORCEINLINE static void process_sixteen_pixels_h_float_pix16_sub4_ks_4_8_16(
  const float* src, int x, float* current_coeff_base,
  int filter_size, // 8, 16, 24, 32 are quasi-constexpr here, others not compile-time known but still aligned to 8
  float* dst,
  ResamplingProgram* program)
{
  assert(program->filter_size_alignment == 8);

  float* current_coeff = current_coeff_base + x * filter_size;
  const int unaligned_kernel_size = program->filter_size_real;
  const __m256 zero256 = _mm256_setzero_ps();

  // --- Block 1: Pixels 0, 1, 2, 3 ---
  __m256 result0 = zero256;
  __m256 result1 = zero256;
  __m256 result2 = zero256;
  __m256 result3 = zero256;

  int begin0 = program->pixel_offset[x + 0];
  int begin1 = program->pixel_offset[x + 1];
  int begin2 = program->pixel_offset[x + 2];
  int begin3 = program->pixel_offset[x + 3];

  process_four_pixels_h_float_pix4of16_ks_4_8_16<is_safe, filtersize_hint>(
    src, begin0, begin1, begin2, begin3, current_coeff, filter_size,
    result0, result1, result2, result3, unaligned_kernel_size);
  current_coeff += 4 * filter_size;

  // --- Block 2: Pixels 4, 5, 6, 7 ---
  __m256 result4 = zero256;
  __m256 result5 = zero256;
  __m256 result6 = zero256;
  __m256 result7 = zero256;

  int begin4 = program->pixel_offset[x + 4];
  int begin5 = program->pixel_offset[x + 5];
  int begin6 = program->pixel_offset[x + 6];
  int begin7 = program->pixel_offset[x + 7];

  process_four_pixels_h_float_pix4of16_ks_4_8_16<is_safe, filtersize_hint>(
    src, begin4, begin5, begin6, begin7, current_coeff, filter_size,
    result4, result5, result6, result7, unaligned_kernel_size);
  current_coeff += 4 * filter_size;

  // ---------------------------------------------------------------------------
  // REDUCTION FOR PIXELS 0-7 (Result256_low)
  // ---------------------------------------------------------------------------

  // Round 1: Reduce pairs (8 vectors -> 4 vectors)
  __m256 sum01 = _mm256_hadd_ps(result0, result1);
  __m256 sum23 = _mm256_hadd_ps(result2, result3);
  __m256 sum45 = _mm256_hadd_ps(result4, result5);
  __m256 sum67 = _mm256_hadd_ps(result6, result7);

  // Round 2: Reduce quads (4 vectors -> 2 vectors)
  __m256 sum0123 = _mm256_hadd_ps(sum01, sum23);
  __m256 sum4567 = _mm256_hadd_ps(sum45, sum67);

  // Round 3: Final Merge (Add Lower 128-bit to Upper 128-bit)
  __m128 lo_0123 = _mm256_castps256_ps128(sum0123);
  __m128 lo_4567 = _mm256_castps256_ps128(sum4567);
  __m256 result_lo = _mm256_insertf128_ps(_mm256_castps128_ps256(lo_0123), lo_4567, 1);

  __m128 hi_0123 = _mm256_extractf128_ps(sum0123, 1);
  __m128 hi_4567 = _mm256_extractf128_ps(sum4567, 1);
  __m256 result_hi = _mm256_insertf128_ps(_mm256_castps128_ps256(hi_0123), hi_4567, 1);

  // Assemble the Low 256-bit result (Pixels 0-7)
  __m256 result256_low = _mm256_add_ps(result_lo, result_hi);


  // --- Block 3: Pixels 8, 9, 10, 11 ---
  __m256 result8 = zero256;
  __m256 result9 = zero256;
  __m256 result10 = zero256;
  __m256 result11 = zero256;

  int begin8 = program->pixel_offset[x + 8];
  int begin9 = program->pixel_offset[x + 9];
  int begin10 = program->pixel_offset[x + 10];
  int begin11 = program->pixel_offset[x + 11];

  process_four_pixels_h_float_pix4of16_ks_4_8_16<is_safe, filtersize_hint>(
    src, begin8, begin9, begin10, begin11, current_coeff, filter_size,
    result8, result9, result10, result11, unaligned_kernel_size);
  current_coeff += 4 * filter_size;

  // --- Block 4: Pixels 12, 13, 14, 15 ---
  __m256 result12 = zero256;
  __m256 result13 = zero256;
  __m256 result14 = zero256;
  __m256 result15 = zero256;

  int begin12 = program->pixel_offset[x + 12];
  int begin13 = program->pixel_offset[x + 13];
  int begin14 = program->pixel_offset[x + 14];
  int begin15 = program->pixel_offset[x + 15];

  process_four_pixels_h_float_pix4of16_ks_4_8_16<is_safe, filtersize_hint>(
    src, begin12, begin13, begin14, begin15, current_coeff, filter_size,
    result12, result13, result14, result15, unaligned_kernel_size);


  // ---------------------------------------------------------------------------
  // REDUCTION FOR PIXELS 8-15 (Result256_high)
  // ---------------------------------------------------------------------------

  // Round 1: Reduce pairs (8 vectors -> 4 vectors)
  __m256 sum89 = _mm256_hadd_ps(result8, result9);
  __m256 sum1011 = _mm256_hadd_ps(result10, result11);
  __m256 sum1213 = _mm256_hadd_ps(result12, result13);
  __m256 sum1415 = _mm256_hadd_ps(result14, result15);

  // Round 2: Reduce quads (4 vectors -> 2 vectors)
  __m256 sum8_11 = _mm256_hadd_ps(sum89, sum1011);
  __m256 sum12_15 = _mm256_hadd_ps(sum1213, sum1415);

  // Round 3: Final Merge (Add Lower 128-bit to Upper 128-bit)
  __m128 lo_8_11 = _mm256_castps256_ps128(sum8_11);
  __m128 lo_12_15 = _mm256_castps256_ps128(sum12_15);
  __m256 result_lo_high = _mm256_insertf128_ps(_mm256_castps128_ps256(lo_8_11), lo_12_15, 1);

  __m128 hi_8_11 = _mm256_extractf128_ps(sum8_11, 1);
  __m128 hi_12_15 = _mm256_extractf128_ps(sum12_15, 1);
  __m256 result_hi_high = _mm256_insertf128_ps(_mm256_castps128_ps256(hi_8_11), hi_12_15, 1);

  // Assemble the High 256-bit result (Pixels 8-15)
  __m256 result256_high = _mm256_add_ps(result_lo_high, result_hi_high);

  // ---------------------------------------------------------------------------
  // Stream the two 256-bit results
  // ---------------------------------------------------------------------------
  _mm256_stream_ps(reinterpret_cast<float*>(dst + x), result256_low);
  _mm256_stream_ps(reinterpret_cast<float*>(dst + x + 8), result256_high);
}

// filtersizealigned8: special: 0, 1..4, Generic : -1
template<int filtersize_hint>
static void internal_resizer_h_avx512_generic_float_pix16_sub4_ks_4_8_16(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel) {
  AVS_UNUSED(bits_per_pixel);
  // filter_size is aligned to 8 (prerequisite), contrary that we have a special case for filter size <=4

  // We note that when template is used, filter_size is quasi-constexpr if filtersize_hint != -1.
  // When filtersize_hint == -1, then program->filter_size is aligned to 8 anyway, but not known at compile time.
  const int filter_size =
    filtersize_hint == 0 ? 8 : // though we'll optimize for 4 internally, coeff buffer is still allocated for 8
    (filtersize_hint >= 1) ? filtersize_hint * 8 : program->filter_size; // this latter is always aligned to 8 as well

  const float* src = (float*)src8;
  float* dst = (float*)dst8;
  dst_pitch = dst_pitch / sizeof(float);
  src_pitch = src_pitch / sizeof(float);

  constexpr int PIXELS_AT_A_TIME = 16;
  // Align safe zone to 16 pixels
  const int w_safe_mod16 = (program->safelimit_16_pixels.overread_possible ? program->safelimit_16_pixels.source_overread_beyond_targetx : width) / PIXELS_AT_A_TIME * PIXELS_AT_A_TIME;

  for (int y = 0; y < height; y++) {
    float* current_coeff_base = program->pixel_coefficient_float;

    // Process safe aligned pixels
    for (int x = 0; x < w_safe_mod16; x += PIXELS_AT_A_TIME) {
      process_sixteen_pixels_h_float_pix16_sub4_ks_4_8_16<true, filtersize_hint>(src, x, current_coeff_base, filter_size, dst, program);
    }

    // Process up to the actual kernel size (unsafe zone)
    for (int x = w_safe_mod16; x < width; x += PIXELS_AT_A_TIME) {
      process_sixteen_pixels_h_float_pix16_sub4_ks_4_8_16<false, filtersize_hint>(src, x, current_coeff_base, filter_size, dst, program);
    }

    dst += dst_pitch;
    src += src_pitch;
  }
}

// Winner implementation: resizer_h_avx512_generic_float_pix16_sub4_ks_4_8_16;
// Other variants kept for reference, speed tested.
// Main test dimensions: pixels per cycle: 8,16,32 (pixX); sub-loops: 2,4,8 (subX); aligned filter sizes (ksX): 4, 8,16
// resizer_h_avx512_generic_float_pix8_sub8_ks16;
// resizer_h_avx512_generic_float_pix16_sub16_ks8;
// resizer_h_avx512_generic_float_pix32_sub8_ks8;
// resizer_h_avx2_generic_float_pix8_sub2_ks8; // like AVX2 version resizer_h_avx2_generic_float
// resizer_h_avx512_generic_float_pix8_sub2_ks8; // like AVX2 version with minor differences
// resizer_h_avx512_generic_float_pix8_sub4_ks8;
// resizer_h_avx512_generic_float_pix16_sub4_ks4;
// resizer_h_avx512_generic_float_pix16_sub4_ks8;

// Features of the chosen implementation:
// - 16 pixels per cycle
// - sub-loop 4 pixels per loop
// - filter size is aligned to 8 (prerequisite)
// - Special cases for aligned filter sizes 4,8,16,24,32
// - Depending on the filter size, calculates in chunks of 16, then 8, then 4 source pixels and coeffs at a time.
void resizer_h_avx512_generic_float_pix16_sub4_ks_4_8_16(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel) {
  const int filter_size = program->filter_size;
  // Expected alignment
  assert(program->filter_size_alignment == 8);

  // Dispatcher template now supports filter_size aligned to 8 (8, 16, 24, 32) and a special case for <=4
  // Larger filter sizes will use the generic method (-1) which still benefit from 16-8-4 coeff processing blocks.
  if (filter_size == 1 * 8)
    if (program->filter_size_real <= 4)
      internal_resizer_h_avx512_generic_float_pix16_sub4_ks_4_8_16<0>(dst8, src8, dst_pitch, src_pitch, program, width, height, bits_per_pixel); // Internally optimized for 4
    else
      internal_resizer_h_avx512_generic_float_pix16_sub4_ks_4_8_16<1>(dst8, src8, dst_pitch, src_pitch, program, width, height, bits_per_pixel); // Internally optimized for 8
  else if (filter_size == 2 * 8) // Internally optimized for 16
    internal_resizer_h_avx512_generic_float_pix16_sub4_ks_4_8_16<2>(dst8, src8, dst_pitch, src_pitch, program, width, height, bits_per_pixel);
  else if (filter_size == 3 * 8) // Internally optimized for 16+8
    internal_resizer_h_avx512_generic_float_pix16_sub4_ks_4_8_16<3>(dst8, src8, dst_pitch, src_pitch, program, width, height, bits_per_pixel);
  else if (filter_size == 4 * 8) // Internally optimized for 2*16
    internal_resizer_h_avx512_generic_float_pix16_sub4_ks_4_8_16<4>(dst8, src8, dst_pitch, src_pitch, program, width, height, bits_per_pixel);
  else // -1: basic method, use program->filter_size, internally optimized for calculating coeffs in N*16 + 8 + 4 + 2 + 1 blocks
    internal_resizer_h_avx512_generic_float_pix16_sub4_ks_4_8_16<-1>(dst8, src8, dst_pitch, src_pitch, program, width, height, bits_per_pixel);
}
