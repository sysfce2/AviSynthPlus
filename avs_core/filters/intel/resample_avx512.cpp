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
  int max_scanlines = detect_optimal_scanline(program->source_size, program->target_size, program->filter_size, cache_size_L2)
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

/* Universal function supporting 2 ways of processing depending on the max offset of the source samples to read in the resampling program :
1. For high upsampling ratios it uses low read (single 8 float source samples) and permute-transpose before V-fma
2. For downsample and no-resize convolution - use each input sequence gathering by direct addressing
*/
template<int filtersizemod4>
void resize_h_planar_float_avx512_gather_permutex_vstripe_ks4(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel)
{
  assert(filtersizemod4 >= 0 && filtersizemod4 <= 3);

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

  bool bDoGather = false;
  // Analyse input resampling program to select method of processing
  for (int x = 0; x < width - 16; x += 16) // -16 to save from vector overrread at program->pixel_offset[x + 15 + 3]; ?
  {
    int start_off = program->pixel_offset[x + 0];
    int end_off = program->pixel_offset[x + 15];
    if ((end_off - start_off) + (program->filter_size_real - 1) > 32) bDoGather = true;

    start_off = program->pixel_offset[x + 1];
    end_off = program->pixel_offset[x + 15 + 1];
    if ((end_off - start_off) + (program->filter_size_real - 1) > 32) bDoGather = true;

    start_off = program->pixel_offset[x + 2];
    end_off = program->pixel_offset[x + 15 + 2];
    if ((end_off - start_off) + (program->filter_size_real - 1) > 32) bDoGather = true;

    start_off = program->pixel_offset[x + 3];
    end_off = program->pixel_offset[x + 15 + 3];
    if ((end_off - start_off) + (program->filter_size_real - 1) > 32) bDoGather = true;
  }

  int x = 0;

  if (bDoGather)
  {
    // This 'auto' lambda construct replaces the need of templates
    auto do_h_float_core = [&](auto partial_load) {
      // Load up to 4x4 coefficients at once before the height loop.
      // Pre-loading and transposing coefficients keeps register usage efficient.
      // Assumes 'filter_size_aligned' is at least 4.

      // Coefficients for the source pixel offset (for src_ptr + begin1 [0..3], begin5 [0..3], begin9 [0..3], begin13 [0..3])
      __m512 coef_1_5_9_13 = _mm512_load_4_m128(current_coeff + filter_size * 0, current_coeff + filter_size * 4, current_coeff + filter_size * 8, current_coeff + filter_size * 12);
      __m512 coef_2_6_10_14 = _mm512_load_4_m128(current_coeff + filter_size * 1, current_coeff + filter_size * 5, current_coeff + filter_size * 9, current_coeff + filter_size * 13);
      __m512 coef_3_7_11_15 = _mm512_load_4_m128(current_coeff + filter_size * 2, current_coeff + filter_size * 6, current_coeff + filter_size * 10, current_coeff + filter_size * 14);
      __m512 coef_4_8_12_16 = _mm512_load_4_m128(current_coeff + filter_size * 3, current_coeff + filter_size * 7, current_coeff + filter_size * 11, current_coeff + filter_size * 15);

      _MM_TRANSPOSE16_LANE4_PS(coef_1_5_9_13, coef_2_6_10_14, coef_3_7_11_15, coef_4_8_12_16);

      float* AVS_RESTRICT dst_ptr = dst + x;
      const float* src_ptr = src;

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

      for (int y = 0; y < height; y++)
      {
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

        _MM_TRANSPOSE16_LANE4_PS(data_1_5_9_13, data_2_6_10_14, data_3_7_11_15, data_4_8_12_16);

        __m512 result = _mm512_mul_ps(data_1_5_9_13, coef_1_5_9_13);
        result = _mm512_fmadd_ps(data_2_6_10_14, coef_2_6_10_14, result);
        result = _mm512_fmadd_ps(data_3_7_11_15, coef_3_7_11_15, result);
        result = _mm512_fmadd_ps(data_4_8_12_16, coef_4_8_12_16, result);

        _mm512_store_ps(dst_ptr, result);

        dst_ptr += dst_pitch;
        src_ptr += src_pitch;
      } // y
      current_coeff += filter_size * 16; // Move to the next set of coefficients for the next 16 output pixels
      }; // end of lambda

    // Process the 'safe zone' where direct full unaligned loads are acceptable.
    for (; x < width_safe_mod; x += PIXELS_AT_A_TIME)
    {
      do_h_float_core(std::false_type{}); // partial_load == false, use direct _mm512_loadu_4_m128
    }

    // Process the potentially 'unsafe zone' near the image edge, using safe loading.
    for (; x < width; x += PIXELS_AT_A_TIME)
    {
      do_h_float_core(std::true_type{}); // partial_load == true, use the safer '_mm512_load_partial_safe_4_m128'
    }
  }
  else // if(bDoGather)
  {
    for (int x = 0; x < width; x += 16)
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
}

template void resize_h_planar_float_avx512_gather_permutex_vstripe_ks4<0>(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel);
template void resize_h_planar_float_avx512_gather_permutex_vstripe_ks4<1>(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel);
template void resize_h_planar_float_avx512_gather_permutex_vstripe_ks4<2>(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel);
template void resize_h_planar_float_avx512_gather_permutex_vstripe_ks4<3>(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel);
