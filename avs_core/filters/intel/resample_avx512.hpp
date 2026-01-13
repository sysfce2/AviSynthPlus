// AviSynth+.  Copyright 2026- AviSynth+ Project
// https://avs-plus.net
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

/*

This is a common source cpp include file (not header) for multi-arch AVX512 functions.
Functions here are static, they will be compiled into each translation unit including this file.

*/

// helper function for simulating _mm512_permutex2var_epi8 when VBMI is not available
// The MSB bit (128) zeroing effect is _not_ considered here, the indices must be all positive and within 0-127 range.
template<bool UseVBMI>
static AVS_FORCEINLINE __m512i _mm512_permutex2var_epi8_SIMUL(__m512i a, __m512i idx, __m512i b) {
  if constexpr (UseVBMI) {
    return _mm512_permutex2var_epi8(a, idx, b);
  }
  else {
    // Constants
    const __m512i v_one = _mm512_set1_epi16(1);

    // 1. Extract the byte indices for the first 32 and last 32 target pixels
    // We expand them to 16-bit so we can treat them as word-indices
    __m512i idx_lo = _mm512_cvtepu8_epi16(_mm512_castsi512_si256(idx));
    __m512i idx_hi = _mm512_cvtepu8_epi16(_mm512_extracti64x4_epi64(idx, 1));

    // Helper to process 32 bytes of the result at a time
    auto get_32_bytes = [&](__m512i target_idx) {
      // word_idx = byte_idx / 2
      __m512i word_idx = _mm512_srli_epi16(target_idx, 1);

      // VPERMT2W: Full 512-bit cross-lane word shuffle from 128-byte pool [a, b]
      __m512i words = _mm512_permutex2var_epi16(a, word_idx, b);

      // If the original byte index was odd, we need the High Byte of the word.
      // We shift those words right by 8 to put the High Byte into the Low Byte position.
      __mmask32 mask_odd = _mm512_test_epi16_mask(target_idx, v_one);
      words = _mm512_mask_srli_epi16(words, mask_odd, words, 8);

      // VPMOVWB: Truncates 32 words to 32 bytes LINEARLY (No lane scrambling)
      // Returns a __m256i
      return _mm512_cvtepi16_epi8(words);
      };

    // 2. Build the two 256-bit halves
    __m256i res_0_31 = get_32_bytes(idx_lo);
    __m256i res_32_63 = get_32_bytes(idx_hi);

    // 3. Combine into final __m512i
    return _mm512_inserti64x4(_mm512_castsi256_si512(res_0_31), res_32_63, 1);
  }
}


// H-Float-Resampler: 16 pixels, filter size 4, transpose 4x (4x_m128) to 4x_m512
// Transposes a 4x4 matrix of 4-float vectors (16x16 float matrix effectively).
// Input/Output: Four 512-bit vectors (16 floats each) passed by reference.
AVS_FORCEINLINE static void _MM_TRANSPOSE16_LANE4_PS(__m512& row0, __m512& row1, __m512& row2, __m512& row3) {
  // Stage 1: Interleave 32-bit (float) elements within 128-bit chunks (lanes)
  // t0 = (r0_lo, r1_lo) | t1 = (r0_hi, r1_hi)
  // t2 = (r2_lo, r3_lo) | t3 = (r2_hi, r3_hi)
  __m512 t0 = _mm512_unpacklo_ps(row0, row1);
  __m512 t1 = _mm512_unpackhi_ps(row0, row1);
  __m512 t2 = _mm512_unpacklo_ps(row2, row3);
  __m512 t3 = _mm512_unpackhi_ps(row2, row3);

  // Stage 2: Shuffle 128-bit chunks (lanes) to complete the transpose
  // We use _mm512_shuffle_ps which shuffles 64-bit blocks across the 512-bit register.
  // _MM_SHUFFLE(w, z, y, x) applies to the 64-bit pairs (4 floats) within each 128-bit lane.
  // Result: row0 = columns 0, 1, 2, 3
  row0 = _mm512_shuffle_ps(t0, t2, _MM_SHUFFLE(1, 0, 1, 0));
  // Result: row1 = columns 4, 5, 6, 7
  row1 = _mm512_shuffle_ps(t0, t2, _MM_SHUFFLE(3, 2, 3, 2));
  // Result: row2 = columns 8, 9, 10, 11
  row2 = _mm512_shuffle_ps(t1, t3, _MM_SHUFFLE(1, 0, 1, 0));
  // Result: row3 = columns 12, 13, 14, 15
  row3 = _mm512_shuffle_ps(t1, t3, _MM_SHUFFLE(3, 2, 3, 2));
}

// H-float-resampler: 16 pixels, filter size 8, transpose 8x (2x_m256) to 8x_m512
// Transposes an 8x8 matrix of 2-float vectors (16x16 float matrix).
// Input/Output: Eight 512-bit vectors (16 floats each) passed by reference.
AVS_FORCEINLINE static void _MM_TRANSPOSE8x16_PS(
  __m512& r0, __m512& r1, __m512& r2, __m512& r3,
  __m512& r4, __m512& r5, __m512& r6, __m512& r7)
{
  // --- Stage 1: Unpack 32-bit (Pairs of rows) ---
  __m512 t0 = _mm512_unpacklo_ps(r0, r1);
  __m512 t1 = _mm512_unpackhi_ps(r0, r1);
  __m512 t2 = _mm512_unpacklo_ps(r2, r3);
  __m512 t3 = _mm512_unpackhi_ps(r2, r3);
  __m512 t4 = _mm512_unpacklo_ps(r4, r5);
  __m512 t5 = _mm512_unpackhi_ps(r4, r5);
  __m512 t6 = _mm512_unpacklo_ps(r6, r7);
  __m512 t7 = _mm512_unpackhi_ps(r6, r7);

  // --- Stage 2: Unpack 64-bit (Quads of rows) ---
  // Uses _mm512_unpacklo/hi_pd for 64-bit (double) to interleave pairs of __m512 floats
  __m512 u0 = _mm512_castpd_ps(_mm512_unpacklo_pd(_mm512_castps_pd(t0), _mm512_castps_pd(t2)));
  __m512 u1 = _mm512_castpd_ps(_mm512_unpackhi_pd(_mm512_castps_pd(t0), _mm512_castps_pd(t2)));
  __m512 u2 = _mm512_castpd_ps(_mm512_unpacklo_pd(_mm512_castps_pd(t1), _mm512_castps_pd(t3)));
  __m512 u3 = _mm512_castpd_ps(_mm512_unpackhi_pd(_mm512_castps_pd(t1), _mm512_castps_pd(t3)));
  __m512 u4 = _mm512_castpd_ps(_mm512_unpacklo_pd(_mm512_castps_pd(t4), _mm512_castps_pd(t6)));
  __m512 u5 = _mm512_castpd_ps(_mm512_unpackhi_pd(_mm512_castps_pd(t4), _mm512_castps_pd(t6)));
  __m512 u6 = _mm512_castpd_ps(_mm512_unpacklo_pd(_mm512_castps_pd(t5), _mm512_castps_pd(t7)));
  __m512 u7 = _mm512_castpd_ps(_mm512_unpackhi_pd(_mm512_castps_pd(t5), _mm512_castps_pd(t7)));

  // --- Stage 3: Shuffle 128-bit lanes (Octets of rows) ---
  // _mm512_shuffle_f32x4 shuffles the 128-bit (f32x4) sub-vectors within and between two __m512 vectors.
  // 0x88 = (10001000)_2: selects lane 0 from first input and lane 0 from second input for lo/hi 256 bits.
  // 0xDD = (11011101)_2: selects lane 3 from first input and lane 3 from second input for lo/hi 256 bits.
  __m512 v0 = _mm512_shuffle_f32x4(u0, u4, 0x88); // Col 0, 4 (interleaved)
  __m512 v1 = _mm512_shuffle_f32x4(u0, u4, 0xDD); // Col 1, 5 (interleaved)
  __m512 v2 = _mm512_shuffle_f32x4(u1, u5, 0x88); // Col 2, 6 (interleaved)
  __m512 v3 = _mm512_shuffle_f32x4(u1, u5, 0xDD); // Col 3, 7 (interleaved)
  __m512 v4 = _mm512_shuffle_f32x4(u2, u6, 0x88); // Col 8, 12 (interleaved)
  __m512 v5 = _mm512_shuffle_f32x4(u2, u6, 0xDD); // Col 9, 13 (interleaved)
  __m512 v6 = _mm512_shuffle_f32x4(u3, u7, 0x88); // Col 10, 14 (interleaved)
  __m512 v7 = _mm512_shuffle_f32x4(u3, u7, 0xDD); // Col 11, 15 (interleaved)

  // --- Stage 4: Permute to Linearize Indices ---
  // Corrects the order of the 128-bit lanes to linearize the columns.
  // The columns are currently: (0, 4, 8, 12, 1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15)
  // The required order is: (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)
  __m512i idx = _mm512_setr_epi32(
    0, 4, 1, 5,    /* Lane 0: Rows 0, 1, 2, 3 */
    2, 6, 3, 7,    /* Lane 1: Rows 4, 5, 6, 7 */
    8, 12, 9, 13, /* Lane 2: Rows 8, 9, 10, 11 */
    10, 14, 11, 15 /* Lane 3: Rows 12, 13, 14, 15 */
  );

  // --- Final Assignment with Correct Mapping ---
  // Maps the permuted vector components back to the original row variables (now columns).
  r0 = _mm512_permutexvar_ps(idx, v0); // Col 0
  r1 = _mm512_permutexvar_ps(idx, v2); // Col 1
  r2 = _mm512_permutexvar_ps(idx, v4); // Col 2
  r3 = _mm512_permutexvar_ps(idx, v6); // Col 3
  r4 = _mm512_permutexvar_ps(idx, v1); // Col 4
  r5 = _mm512_permutexvar_ps(idx, v3); // Col 5
  r6 = _mm512_permutexvar_ps(idx, v5); // Col 6
  r7 = _mm512_permutexvar_ps(idx, v7); // Col 7
}

// Loads two 256-bit float vectors from registers (__m256) into a single 512-bit register.
// Equivalent to _mm512_insertf32x8(_mm512_castps256_ps512(lo), hi, 1)
AVS_FORCEINLINE static __m512 _mm512_insert_2_m256(__m256 lo, __m256 hi) {
  return _mm512_insertf32x8(_mm512_castps256_ps512(lo), hi, 1);
}

// Loads four 128-bit float vectors (unaligned) into a single 512-bit register.
AVS_FORCEINLINE static __m512 _mm512_loadu_4_m128(
  /* __m128 const* */ const float* addr1,
  /* __m128 const* */ const float* addr2,
  /* __m128 const* */ const float* addr3,
  /* __m128 const* */ const float* addr4)
{
  // The cast is needed for the first insertion to make the target a 512-bit register
  __m512 v = _mm512_castps128_ps512(_mm_loadu_ps(addr1));
  v = _mm512_insertf32x4(v, _mm_loadu_ps(addr2), 1);
  v = _mm512_insertf32x4(v, _mm_loadu_ps(addr3), 2);
  v = _mm512_insertf32x4(v, _mm_loadu_ps(addr4), 3);
  return v;
}

// Loads four 128-bit float vectors (aligned) into a single 512-bit register.
AVS_FORCEINLINE static __m512 _mm512_load_4_m128(
  /* __m128 const* */ const float* addr1,
  /* __m128 const* */ const float* addr2,
  /* __m128 const* */ const float* addr3,
  /* __m128 const* */ const float* addr4)
{
  // The cast is needed for the first insertion to make the target a 512-bit register
  __m512 v = _mm512_castps128_ps512(_mm_load_ps(addr1));
  v = _mm512_insertf32x4(v, _mm_load_ps(addr2), 1);
  v = _mm512_insertf32x4(v, _mm_load_ps(addr3), 2);
  v = _mm512_insertf32x4(v, _mm_load_ps(addr4), 3);
  return v;
}

// Loads two 256 - bit unaligned integer vectors from registers(__m256i) into a single 512i register.
AVS_FORCEINLINE static __m512i _mm512i_loadu_2_m256i(
  /* __m256i const* */ const __m256i* addr1,
  /* __m256i const* */ const __m256i* addr2)
{
  return _mm512_inserti64x4(_mm512_zextsi256_si512(_mm256_loadu_si256(addr1)), _mm256_loadu_si256(addr2), 1);
}

// Loads two 256 - bit aligned integer vectors from registers(__m256) into a single 512 - bit register.
AVS_FORCEINLINE static __m512i _mm512i_load_2_m256i(
  /* __m256i const* */ const __m256i* addr1,
  /* __m256i const* */ const __m256i* addr2)
{
  return _mm512_inserti64x4(_mm512_zextsi256_si512(_mm256_load_si256(addr1)), _mm256_load_si256(addr2), 1);
}

// Integers
// Loads four 128-bit integer vectors (unaligned) into a single 512-bit integer register.
AVS_FORCEINLINE static __m512i _mm512i_loadu_4_m128i(
  /* __m128i const* */ const __m128i* addr1,
  /* __m128i const* */ const __m128i* addr2,
  /* __m128i const* */ const __m128i* addr3,
  /* __m128i const* */ const __m128i* addr4)
{
  // The cast is needed for the first insertion to make the target a 512-bit register
  __m512i v = _mm512_zextsi128_si512(_mm_loadu_si128(addr1));
  v = _mm512_inserti32x4(v, _mm_loadu_si128(addr2), 1);
  v = _mm512_inserti32x4(v, _mm_loadu_si128(addr3), 2);
  v = _mm512_inserti32x4(v, _mm_loadu_si128(addr4), 3);
  return v;
}

