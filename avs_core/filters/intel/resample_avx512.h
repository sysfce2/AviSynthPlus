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

#ifndef __Resample_AVX512_H__
#define __Resample_AVX512_H__

#include <avisynth.h>
#include "../resample_functions.h"

bool resize_h_planar_float_avx512_gather_permutex_vstripe_ks4_check(ResamplingProgram* program);
void resizer_h_avx512_generic_float_pix16_sub4_ks_4_8_16(BYTE * dst8, const BYTE * src8, int dst_pitch, int src_pitch, ResamplingProgram * program, int width, int height, int bits_per_pixel);

// useful macros

// Full 512-bit version of transpose16
#define _MM_TRANSPOSE16_LANE4_PS(row0, row1, row2, row3) \
  do { \
    __m512 __t0, __t1, __t2, __t3; \
    __t0 = _mm512_unpacklo_ps(row0, row1); \
    __t1 = _mm512_unpackhi_ps(row0, row1); \
    __t2 = _mm512_unpacklo_ps(row2, row3); \
    __t3 = _mm512_unpackhi_ps(row2, row3); \
    row0 = _mm512_shuffle_ps(__t0, __t2, _MM_SHUFFLE(1, 0, 1, 0)); \
    row1 = _mm512_shuffle_ps(__t0, __t2, _MM_SHUFFLE(3, 2, 3, 2)); \
    row2 = _mm512_shuffle_ps(__t1, __t3, _MM_SHUFFLE(1, 0, 1, 0)); \
    row3 = _mm512_shuffle_ps(__t1, __t3, _MM_SHUFFLE(3, 2, 3, 2)); \
  } while (0)

#ifndef _mm512_loadu_4_m128
#define _mm512_loadu_4_m128(/* __m128 const* */ addr1, \
                            /* __m128 const* */ addr2, \
                            /* __m128 const* */ addr3, \
                            /* __m128 const* */ addr4) \
_mm512_insertf32x4(_mm512_insertf32x4(_mm512_insertf32x4(_mm512_castps128_ps512(_mm_loadu_ps(addr1)), _mm_loadu_ps(addr2), 1), _mm_loadu_ps(addr3), 2), _mm_loadu_ps(addr4), 3)
#endif

#ifndef _mm512_load_4_m128
#define _mm512_load_4_m128(/* __m128 const* */ addr1, \
                            /* __m128 const* */ addr2, \
                            /* __m128 const* */ addr3, \
                            /* __m128 const* */ addr4) \
_mm512_insertf32x4(_mm512_insertf32x4(_mm512_insertf32x4(_mm512_castps128_ps512(_mm_load_ps(addr1)), _mm_load_ps(addr2), 1), _mm_load_ps(addr3), 2), _mm_load_ps(addr4), 3)
#endif

// Transpose 4x4 blocks within each lane
#define _MM_TRANSPOSE8_LANE4_PS(row0, row1, row2, row3) \
  do { \
    __m256 __t0, __t1, __t2, __t3; \
    __t0 = _mm256_unpacklo_ps(row0, row1); \
    __t1 = _mm256_unpackhi_ps(row0, row1); \
    __t2 = _mm256_unpacklo_ps(row2, row3); \
    __t3 = _mm256_unpackhi_ps(row2, row3); \
    row0 = _mm256_shuffle_ps(__t0, __t2, _MM_SHUFFLE(1, 0, 1, 0)); \
    row1 = _mm256_shuffle_ps(__t0, __t2, _MM_SHUFFLE(3, 2, 3, 2)); \
    row2 = _mm256_shuffle_ps(__t1, __t3, _MM_SHUFFLE(1, 0, 1, 0)); \
    row3 = _mm256_shuffle_ps(__t1, __t3, _MM_SHUFFLE(3, 2, 3, 2)); \
  } while (0)

// a 256-bit version of transpose16 but not quicker than the full 512-bit version
#define _MM_TRANSPOSE16_LANE4_PS_256(row0, row1, row2, row3) \
  do { \
    /* Low Half (256-bit): Use cast for ZMM[0:255]. This is typically a zero-latency register rename on Intel. */ \
    __m256 row0_low = _mm512_castps512_ps256(row0); \
    __m256 row1_low = _mm512_castps512_ps256(row1); \
    __m256 row2_low = _mm512_castps512_ps256(row2); \
    __m256 row3_low = _mm512_castps512_ps256(row3); \
\
    /* 1. Transpose the LOW (left half) using the efficient AVX2 macro (parallel ports) */ \
    _MM_TRANSPOSE8_LANE4_PS(row0_low, row1_low, row2_low, row3_low); \
\
    /* High Half (256-bit): Assemble ZMM[256:511] by extracting the two high 128-bit chunks (index 2 and 3) */ \
    __m256 row0_high = _mm256_insertf32x4(_mm256_castps128_ps256(_mm512_extractf32x4_ps(row0, 2)), _mm512_extractf32x4_ps(row0, 3), 1); \
    __m256 row1_high = _mm256_insertf32x4(_mm256_castps128_ps256(_mm512_extractf32x4_ps(row1, 2)), _mm512_extractf32x4_ps(row1, 3), 1); \
    __m256 row2_high = _mm256_insertf32x4(_mm256_castps128_ps256(_mm512_extractf32x4_ps(row2, 2)), _mm512_extractf32x4_ps(row2, 3), 1); \
    __m256 row3_high = _mm256_insertf32x4(_mm256_castps128_ps256(_mm512_extractf32x4_ps(row3, 2)), _mm512_extractf32x4_ps(row3, 3), 1); \
\
    /* 2. Transpose the HIGH (right half) using the efficient AVX2 macro (parallel ports) */ \
    _MM_TRANSPOSE8_LANE4_PS(row0_high, row1_high, row2_high, row3_high); \
\
    /* 3. Re-assemble the results back into the 512-bit output vectors. */ \
    row0 = _mm512_insertf32x8(_mm512_castps256_ps512(row0_low), row0_high, 1); \
    row1 = _mm512_insertf32x8(_mm512_castps256_ps512(row1_low), row1_high, 1); \
    row2 = _mm512_insertf32x8(_mm512_castps256_ps512(row2_low), row2_high, 1); \
    row3 = _mm512_insertf32x8(_mm512_castps256_ps512(row3_low), row3_high, 1); \
\
  } while (0)

#endif // __Resample_AVX512_H__
