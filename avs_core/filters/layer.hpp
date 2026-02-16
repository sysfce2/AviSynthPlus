
// Basic utils in C++ language, to be included in both base and processor specific (e.g. avx2) source modules, where they can
// be optimized for the specific instruction set.
// Chroma placement to mask helpers, yuv add, subtract, mul, lighten, darken

DISABLE_WARNING_PUSH
DISABLE_WARNING_UNREFERENCED_LOCAL_VARIABLE

// Single function for all cases

// Helper function for calculation of effective mask value for a pixel, based on chroma placement
// integer 8-16 bits version, with in/out parameter for MPEG2 sliding window modes
template<MaskMode maskMode, typename pixel_t>
AVS_FORCEINLINE static int calculate_effective_mask(
  const pixel_t* ptr,
  int x,
  int pitch,
  int& right_value  // in/out parameter for MPEG2 sliding window modes
) {
  if constexpr (maskMode == MASK444) {
    // +------+
    // | 1.0  |
    // +------+
    return ptr[x];
  }
  else if constexpr (maskMode == MASK411) {
    // +------+------+------+------+
    // | 0.25 | 0.25 | 0.25 | 0.25 |
    // +------+------+------+------+
    return (ptr[x * 4] + ptr[x * 4 + 1] + ptr[x * 4 + 2] + ptr[x * 4 + 3] + 2) >> 2;
  }
  else if constexpr (maskMode == MASK420) {
    // +------+------+
    // | 0.25 | 0.25 |
    // |------+------|
    // | 0.25 | 0.25 |
    // +------+------+
    return (ptr[x * 2] + ptr[x * 2 + 1] + ptr[x * 2 + pitch] + ptr[x * 2 + 1 + pitch] + 2) >> 2;
  }
  else if constexpr (maskMode == MASK420_MPEG2) {
    // ------+------+-------+
    // 0.125 | 0.25 | 0.125 |
    // ------|------+-------|
    // 0.125 | 0.25 | 0.125 |
    // ------+------+-------+
    int left = right_value;
    const int mid = ptr[x * 2] + ptr[x * 2 + pitch];
    right_value = ptr[x * 2 + 1] + ptr[x * 2 + 1 + pitch];
    return (left + 2 * mid + right_value + 4) >> 3;
  }
  else if constexpr (maskMode == MASK422) {
    // +------+------+
    // | 0.5  | 0.5  |
    // +------+------+
    return (ptr[x * 2] + ptr[x * 2 + 1] + 1) >> 1;
  }
  else if constexpr (maskMode == MASK422_MPEG2) {
    // ------+------+-------+
    // 0.25  | 0.5  | 0.25  |
    // ------+------+-------+
    int left = right_value;
    const int mid = ptr[x * 2];
    right_value = ptr[x * 2 + 1];
    return (left + 2 * mid + right_value + 2) >> 2;
  }
}

// Helper function for calculation of effective mask value for a pixel, based on chroma placement
// Float version, with in/out parameter for MPEG2 sliding window modes
template<MaskMode maskMode>
AVS_FORCEINLINE static float calculate_effective_mask_f(
  const float* ptr,
  int x,
  int pitch,
  float& right_value  // in/out parameter for MPEG2 sliding window modes
) {
  if constexpr (maskMode == MASK444) {
    return ptr[x];
  }
  else if constexpr (maskMode == MASK411) {
    return (ptr[x * 4] + ptr[x * 4 + 1] + ptr[x * 4 + 2] + ptr[x * 4 + 3]) * 0.25f;
  }
  else if constexpr (maskMode == MASK420) {
    return (ptr[x * 2] + ptr[x * 2 + 1] + ptr[x * 2 + pitch] + ptr[x * 2 + 1 + pitch]) * 0.25f;
  }
  else if constexpr (maskMode == MASK420_MPEG2) {
    float left = right_value;
    const float mid = ptr[x * 2] + ptr[x * 2 + pitch];
    right_value = ptr[x * 2 + 1] + ptr[x * 2 + 1 + pitch];
    return (left + 2.0f * mid + right_value) * 0.125f;
  }
  else if constexpr (maskMode == MASK422) {
    return (ptr[x * 2] + ptr[x * 2 + 1]) * 0.5f;
  }
  else if constexpr (maskMode == MASK422_MPEG2) {
    float left = right_value;
    const float mid = ptr[x * 2];
    right_value = ptr[x * 2 + 1];
    return (left + 2.0f * mid + right_value) * 0.25f;
  }
}

// Helper function to prepare effective mask pointer for a complete row
// integer 8-16 bits version
template<MaskMode maskMode, typename pixel_t>
AVS_FORCEINLINE static const pixel_t* prepare_effective_mask_for_row(
  const pixel_t* maskp,
  int mask_pitch,
  int width,
  std::vector<pixel_t>& effective_mask_buffer
) {
  if constexpr (maskMode == MASK444) {
    // Direct access to original mask - no pre-calculation needed
    return maskp;
  }
  else {
    // Initialize sliding window state for MPEG2 modes
    int mask_right = 0; // initialized to silence warnings, actual value set below
    if constexpr (maskMode == MASK420_MPEG2) {
      mask_right = maskp[0] + maskp[0 + mask_pitch];
    }
    else if constexpr (maskMode == MASK422_MPEG2) {
      mask_right = maskp[0];
    }

    // Pre-calculate averaged mask values
    for (int x = 0; x < width; ++x) {
      effective_mask_buffer[x] = (pixel_t)calculate_effective_mask<maskMode>(maskp, x, mask_pitch, mask_right);
    }
    return effective_mask_buffer.data();
  }
}

// Helper function to prepare effective mask pointer for a complete row
// float version
template<MaskMode maskMode>
AVS_FORCEINLINE static const float* prepare_effective_mask_for_row_f(
  const float* maskp,
  int mask_pitch,
  int width,
  std::vector<float>& effective_mask_buffer
) {
  if constexpr (maskMode == MASK444) {
    // Direct access to original mask - no pre-calculation needed
    return maskp;
  }
  else {
    // Initialize sliding window state for MPEG2 modes
    float mask_right = 0.0f; // initialized to silence warnings, actual value set below
    if constexpr (maskMode == MASK420_MPEG2) {
      mask_right = maskp[0] + maskp[0 + mask_pitch];
    }
    else if constexpr (maskMode == MASK422_MPEG2) {
      mask_right = maskp[0];
    }

    // Pre-calculate averaged mask values
    for (int x = 0; x < width; ++x) {
      effective_mask_buffer[x] = calculate_effective_mask_f<maskMode>(maskp, x, mask_pitch, mask_right);
    }
    return effective_mask_buffer.data();
  }
}

// YUV(A) mul 8-16 bits
// when chroma is processed, one can use/not use source chroma,
// Only when use_alpha: maskMode defines mask generation for chroma planes
// When use_alpha == false maskMode ignored
template<MaskMode maskMode, typename pixel_t, bool lessthan16bits, bool is_chroma, bool use_chroma, bool has_alpha>
static void layer_yuv_mul_c(BYTE* dstp8, const BYTE* ovrp8, const BYTE* maskp8, int dst_pitch, int overlay_pitch, int mask_pitch, int width, int height, int level, int bits_per_pixel) {
  pixel_t* dstp = reinterpret_cast<pixel_t*>(dstp8);
  const pixel_t* ovrp = reinterpret_cast<const pixel_t*>(ovrp8);
  const pixel_t* maskp = reinterpret_cast<const pixel_t*>(maskp8);
  dst_pitch /= sizeof(pixel_t);
  overlay_pitch /= sizeof(pixel_t);
  mask_pitch /= sizeof(pixel_t);

  typedef typename std::conditional<lessthan16bits, int, int64_t>::type calc_t;

  if constexpr (sizeof(pixel_t) == 1)
    bits_per_pixel = 8; // make quasi constexpr
  else if constexpr (sizeof(pixel_t) == 2 && !lessthan16bits)
    bits_per_pixel = 16; // make quasi constexpr

  // precalculate mask buffer
  std::vector<pixel_t> effective_mask_buffer;
  if constexpr (has_alpha && maskMode != MASK444) {
    effective_mask_buffer.resize(width);
  }

  for (int y = 0; y < height; ++y) {
    const pixel_t* effective_mask_ptr = nullptr;

    // precalculate effective mask for this row
    if constexpr (has_alpha) {
      effective_mask_ptr = prepare_effective_mask_for_row<maskMode, pixel_t>(maskp, mask_pitch, width, effective_mask_buffer);
    }

    // Main blending loop - now simplified and vectorizable
    for (int x = 0; x < width; ++x) {
      int effective_mask = has_alpha ? effective_mask_ptr[x] : 0;
      int alpha_mask = has_alpha ? (int)(((calc_t)effective_mask * level + 1) >> bits_per_pixel) : level;

      // fixme: no rounding? (code from YUY2)
      // for mul: no.
      if constexpr (!is_chroma)
        dstp[x] = (pixel_t)(dstp[x] + ((((((calc_t)ovrp[x] * dstp[x]) >> bits_per_pixel) - dstp[x]) * alpha_mask) >> bits_per_pixel));
      else if constexpr (use_chroma) {
        // chroma mode + process chroma
        dstp[x] = (pixel_t)(dstp[x] + (((calc_t)(ovrp[x] - dstp[x]) * alpha_mask) >> bits_per_pixel));
        // U = U + ( ((Uovr - U)*level) >> 8 )
        // V = V + ( ((Vovr - V)*level) >> 8 )
      }
      else {
        // non-chroma mode + process chroma
        const int half = 1 << (bits_per_pixel - 1);
        dstp[x] = (pixel_t)(dstp[x] + (((calc_t)(half - dstp[x]) * (alpha_mask / 2)) >> bits_per_pixel));
        // U = U + ( ((128 - U)*(level/2)) >> 8 )
        // V = V + ( ((128 - V)*(level/2)) >> 8 )
      }
    }
    dstp += dst_pitch;
    ovrp += overlay_pitch;
    if constexpr (has_alpha) {
      if constexpr (maskMode == MASK420 || maskMode == MASK420_MPEG2)
        maskp += mask_pitch * 2;
      else
        maskp += mask_pitch;
    }
  }
}

// YUV(A) mul 32 bits
template<MaskMode maskMode, bool is_chroma, bool use_chroma, bool has_alpha>
static void layer_yuv_mul_f_c(BYTE* dstp8, const BYTE* ovrp8, const BYTE* maskp8, int dst_pitch, int overlay_pitch, int mask_pitch, int width, int height, float opacity) {
  float* dstp = reinterpret_cast<float*>(dstp8);
  const float* ovrp = reinterpret_cast<const float*>(ovrp8);
  const float* maskp = reinterpret_cast<const float*>(maskp8);
  dst_pitch /= sizeof(float);
  overlay_pitch /= sizeof(float);
  mask_pitch /= sizeof(float);

  // precalculate mask buffer
  std::vector<float> effective_mask_buffer;
  if constexpr (has_alpha && maskMode != MASK444) {
    effective_mask_buffer.resize(width);
  }

  for (int y = 0; y < height; ++y) {
    const float* effective_mask_ptr = nullptr;

    // precalculate effective mask for this row
    if constexpr (has_alpha) {
      effective_mask_ptr = prepare_effective_mask_for_row_f<maskMode>(maskp, mask_pitch, width, effective_mask_buffer);
    }

    // Main blending loop - now simplified and vectorizable
    for (int x = 0; x < width; ++x) {
      float effective_mask = has_alpha ? effective_mask_ptr[x] : 0.0f;
      float alpha_mask = has_alpha ? effective_mask * opacity : opacity;

      if constexpr (!is_chroma)
        dstp[x] = dstp[x] + (ovrp[x] * dstp[x] - dstp[x]) * alpha_mask;
      else if constexpr (use_chroma) {
        // chroma mode + process chroma
        dstp[x] = dstp[x] + (ovrp[x] - dstp[x]) * alpha_mask;
        // U = U + ( ((Uovr - U)*level) >> 8 )
        // V = V + ( ((Vovr - V)*level) >> 8 )
      }
      else {
        // non-chroma mode + process chroma
        constexpr float half = 0.0f;
        dstp[x] = dstp[x] + (half - dstp[x]) * (alpha_mask * 0.5f);
        // U = U + ( ((128 - U)*(level/2)) >> 8 )
        // V = V + ( ((128 - V)*(level/2)) >> 8 )
      }
    }
    dstp += dst_pitch;
    ovrp += overlay_pitch;
    if constexpr (has_alpha) {
      if constexpr (maskMode == MASK420 || maskMode == MASK420_MPEG2)
        maskp += mask_pitch * 2;
      else
        maskp += mask_pitch;
    }
  }
}

// YUV mul 8-16 bits
// when chroma is processed, one can use/not use source chroma
// Only when use_alpha: maskMode defines mask generation for chroma planes
// When use_alpha == false -> maskMode ignored


// Generated ASM analysis after doing chroma placement-dependent mask precalculation
// => vectorized blending now recognized!
//
// Refactoring separated mask calculation from blending, enabling compilers to auto-vectorize the main loop.
// Memory overhead: single row buffer (e.g., 1920 bytes for 1080p) fits comfortably in L1 cache.
// Estimated speedups were told by AI, but preparing the inputs for Layer takes time, i could not measure only the
// blending loop; but the speedup is significant (except 444 where no precalculation is needed)
//
// Mask calculation (e.g. MASK420_MPEG2 2x2 gather with sliding window):
//   - All compilers: Scalar with no or max. 2x unrolling (complex gather+dependency pattern blocks vectorization)
//   - Minimal overhead: ~1-2% of total runtime, one-time cost per row
//
// Main blending loop vectorization results:
//   MSVC 2022 (SSE4.1):  4-wide vectorization, ~80 instructions per 16 pixels
//     - Uses pmulld (SSE4.1) for 32-bit multiply, pmovzxbd for byte→dword extension
//     - Fallback path: 16→4→1 (main loop, then scalar cleanup)
//     - Estimated speedup: 3-4x vs scalar
//   
//   Intel C++ 2025 (SSE2):  16-wide vectorization, ~150 instructions per 16 pixels
//     - Workaround for missing pmulld: pmuludq + shuffle for odd/even dwords (complex!)
//     - Complex unpacking chain: 4× punpcklbw/punpckhbw + punpcklwd/punpckhwd for byte→dword
//     - Fallback path: 16→1 (main loop processes full 16, scalar cleanup for remainder)
//     - Estimated speedup: 6-8x vs scalar
//   
//   Intel C++ 2025 (AVX2):  16-wide vectorization (2×YMM), ~60 instructions per 16 pixels
//     - Native vpmulld on YMM, vpmovzxbd for clean extension, vpshufb+vpackusdw packing
//     - Fallback path: 16→4→1 (YMM main, XMM for 4-15 remaining, scalar for <4)
//     - Requires XMM6-XMM14 save/restore (ABI requirement), adds minimal function overhead
//     - Estimated speedup: 10-12x vs scalar
//   
//   Intel C++ 2025 (AVX-512):  16-wide vectorization (ZMM), ~40 instructions per 16 pixels
//     - Predicated execution via k-registers: vpcmpuq+kunpckbw for boundary checking
//     - Masked loads/stores (vmovdqu8 {k1}{z}) eliminate separate cleanup loops entirely!
//     - Native vpmovdb for efficient dword→byte packing, vpmovzxbd for extension
//     - Fallback path: 16→1 with masking (5+ remaining uses masked 16-wide, <5 uses scalar)
//     - No XMM register save overhead (ZMM registers are volatile), only vzeroupper at exit
//     - Estimated speedup: 12-15x vs scalar (requires Ice Lake+, Zen4+; may throttle on some CPUs)
//
// Sequential mask access pattern (vs. inline 2x2, 1x2, 2x3 gather) was critical for unlocking
// auto-vectorization. The 16-wide SIMD implementations process the same data in 1/10th the time
// despite the overhead of mask pre-calculation, proving the separation-of-concerns approach.
// 
// Instruction count <> performance; Intel SSE2 uses ~2x more instructions than 
// MSVC SSE4.1, but achieves better speedup due to:
//   1. Better instruction-level parallelism (ILP) - processes all 16 pixels in parallel
//   2. Lower loop overhead - single iteration vs 4 iterations for 16 pixels
//   3. Better memory bandwidth utilization - coalesced loads/stores
//   4. Aggressive register usage (all 16 XMM) reduces memory traffic
// The 16-wide approach's higher upfront cost is amortized by massive parallelism.

// Separated mask precalculation per row.
template<MaskMode maskMode, typename pixel_t, bool lessthan16bits, bool is_chroma, bool use_chroma, bool has_alpha, bool subtract>
static void layer_yuv_add_subtract_c(BYTE* dstp8, const BYTE* ovrp8, const BYTE* mask8, int dst_pitch, int overlay_pitch, int mask_pitch, int width, int height, int level, int bits_per_pixel) {
  pixel_t* dstp = reinterpret_cast<pixel_t*>(dstp8);
  const pixel_t* ovrp = reinterpret_cast<const pixel_t*>(ovrp8);
  const pixel_t* maskp = reinterpret_cast<const pixel_t*>(mask8);
  dst_pitch /= sizeof(pixel_t);
  overlay_pitch /= sizeof(pixel_t);
  mask_pitch /= sizeof(pixel_t);

  typedef typename std::conditional<lessthan16bits, int, int64_t>::type calc_t;

  if constexpr (sizeof(pixel_t) == 1)
    bits_per_pixel = 8; // make quasi constexpr
  else if constexpr (sizeof(pixel_t) == 2 && !lessthan16bits)
    bits_per_pixel = 16; // make quasi constexpr

  const int max_pixel_value = (1 << bits_per_pixel) - 1;
  const int rounder = 1 << (bits_per_pixel - 1);

  // Allocate temporary mask buffer
  std::vector<pixel_t> effective_mask_buffer;
  if constexpr (has_alpha && maskMode != MASK444) {
    effective_mask_buffer.resize(width);
  }

  for (int y = 0; y < height; ++y) {
    const pixel_t* effective_mask_ptr = nullptr;

    // precalculate effective mask for this row
    if constexpr (has_alpha) {
      effective_mask_ptr = prepare_effective_mask_for_row<maskMode, pixel_t>(maskp, mask_pitch, width, effective_mask_buffer);
    }

    // Main blending loop - now simplified
    for (int x = 0; x < width; ++x) {
      int alpha_mask;
      if constexpr (has_alpha) {
        int effective_mask = effective_mask_ptr[x];
        alpha_mask = (int)(((calc_t)effective_mask * level + 1) >> bits_per_pixel);
      }
      else {
        alpha_mask = level;
      }

      if constexpr (subtract) {
        if constexpr (!is_chroma || use_chroma) {
          if constexpr (!is_chroma)
            dstp[x] = (pixel_t)(dstp[x] + (((calc_t)(max_pixel_value - ovrp[x] - dstp[x]) * alpha_mask + rounder) >> bits_per_pixel));
          else
          {
            const int half = 1 << (bits_per_pixel - 1);
            dstp[x] = (pixel_t)(dstp[x] + (((calc_t)(2 * half - ovrp[x] - dstp[x]) * alpha_mask + rounder) >> bits_per_pixel));
          }
        }
        else {
          const int half = 1 << (bits_per_pixel - 1);
          dstp[x] = (pixel_t)(dstp[x] + (((calc_t)(half - dstp[x]) * alpha_mask + rounder) >> bits_per_pixel));
        }
      }
      else {
        if constexpr (!is_chroma || use_chroma)
          dstp[x] = (pixel_t)(dstp[x] + (((calc_t)(ovrp[x] - dstp[x]) * alpha_mask + rounder) >> bits_per_pixel));
        else {
          const int half = 1 << (bits_per_pixel - 1);
          dstp[x] = (pixel_t)(dstp[x] + (((calc_t)(half - dstp[x]) * alpha_mask + rounder) >> bits_per_pixel));
        }
      }
    }
    dstp += dst_pitch;
    ovrp += overlay_pitch;
    if constexpr (has_alpha) {
      if constexpr (maskMode == MASK420 || maskMode == MASK420_MPEG2)
        maskp += mask_pitch * 2;
      else
        maskp += mask_pitch;
    }
  }
}

// YUV(A) add/subtract 32 bits
template<MaskMode maskMode, bool is_chroma, bool use_chroma, bool has_alpha, bool subtract>
static void layer_yuv_add_subtract_f_c(BYTE* dstp8, const BYTE* ovrp8, const BYTE* maskp8, int dst_pitch, int overlay_pitch, int mask_pitch, int width, int height, float opacity) {
  float* dstp = reinterpret_cast<float*>(dstp8);
  const float* ovrp = reinterpret_cast<const float*>(ovrp8);
  const float* maskp = reinterpret_cast<const float*>(maskp8);
  dst_pitch /= sizeof(float);
  overlay_pitch /= sizeof(float);
  mask_pitch /= sizeof(float);

  // precalculate mask buffer
  std::vector<float> effective_mask_buffer;
  if constexpr (has_alpha && maskMode != MASK444) {
    effective_mask_buffer.resize(width);
  }

  for (int y = 0; y < height; ++y) {
    const float* effective_mask_ptr = nullptr;

    // precalculate effective mask for this row
    if constexpr (has_alpha) {
      if constexpr (maskMode == MASK444) {
        // Direct access to original mask - no pre-calculation needed
        effective_mask_ptr = maskp;
      }
      else {
        // Initialize sliding window state for MPEG2 modes
        float mask_right = 0.0f;
        if constexpr (maskMode == MASK420_MPEG2) {
          mask_right = maskp[0] + maskp[0 + mask_pitch];
        }
        else if constexpr (maskMode == MASK422_MPEG2) {
          mask_right = maskp[0];
        }

        // precalculate averaged mask values
        for (int x = 0; x < width; ++x) {
          effective_mask_buffer[x] = calculate_effective_mask_f<maskMode>(
            maskp, x, mask_pitch, mask_right
          );
        }
        effective_mask_ptr = effective_mask_buffer.data();
      }
    }

    // Main blending loop - now simplified and vectorizable
    for (int x = 0; x < width; ++x) {
      float effective_mask = has_alpha ? effective_mask_ptr[x] : 0.0f;
      float alpha_mask = has_alpha ? effective_mask * opacity : opacity;

      if constexpr (subtract) {
        constexpr float ref_pixel_value = is_chroma ? 0.0f : 1.0f;
        if constexpr (!is_chroma || use_chroma) {
          dstp[x] = dstp[x] + (ref_pixel_value - ovrp[x] - dstp[x]) * alpha_mask;
        }
        else {
          dstp[x] = dstp[x] + (ref_pixel_value - dstp[x]) * alpha_mask;
        }
      }
      else {
        if constexpr (!is_chroma || use_chroma) {
          dstp[x] = dstp[x] + (ovrp[x] - dstp[x]) * alpha_mask;
        }
        else {
          constexpr float half = 0.0f;
          dstp[x] = dstp[x] + (half - dstp[x]) * alpha_mask;
        }
      }
    }
    dstp += dst_pitch;
    ovrp += overlay_pitch;
    if constexpr (has_alpha) {
      if constexpr (maskMode == MASK420 || maskMode == MASK420_MPEG2)
        maskp += mask_pitch * 2;
      else
        maskp += mask_pitch;
    }
  }
}

// Unlike RGBA version, YUVA does not update destination alpha
template<int mode, MaskMode maskMode, typename pixel_t, bool lessthan16bits, bool lumaonly, bool has_alpha>
static void layer_yuv_lighten_darken_c(
  BYTE* dstp8, BYTE* dstp8_u, BYTE* dstp8_v,/* BYTE* dstp8_a,*/
  const BYTE* ovrp8, const BYTE* ovrp8_u, const BYTE* ovrp8_v, const BYTE* maskp8,
  int dst_pitch, int dst_pitchUV,
  int overlay_pitch, int overlay_pitchUV,
  int mask_pitch,
  int width, int height, int level, int thresh,
  int bits_per_pixel) {

  pixel_t* dstp = reinterpret_cast<pixel_t*>(dstp8);
  pixel_t* dstp_u = reinterpret_cast<pixel_t*>(dstp8_u);
  pixel_t* dstp_v = reinterpret_cast<pixel_t*>(dstp8_v);
  // pixel_t* dstp_a = reinterpret_cast<pixel_t *>(dstp8_a); // not destination alpha update

  const pixel_t* ovrp = reinterpret_cast<const pixel_t*>(ovrp8);
  const pixel_t* ovrp_u = reinterpret_cast<const pixel_t*>(ovrp8_u);
  const pixel_t* ovrp_v = reinterpret_cast<const pixel_t*>(ovrp8_v);
  const pixel_t* maskp = reinterpret_cast<const pixel_t*>(maskp8);

  dst_pitch /= sizeof(pixel_t);
  dst_pitchUV /= sizeof(pixel_t);
  overlay_pitch /= sizeof(pixel_t);
  overlay_pitchUV /= sizeof(pixel_t);
  mask_pitch /= sizeof(pixel_t);

  const int cwidth = (maskMode == MASK444) ? width : (maskMode == MASK411) ? width >> 2 : width >> 1; // 444:/1  420,422:/2  411:/4
  const int cheight = (maskMode == MASK444 || maskMode == MASK422 || maskMode == MASK422_MPEG2 || maskMode == MASK411) ? height : height >> 1; // 444,422,411:/1  420:/2

  // In lighten/darken we need 3 buffers:
  std::vector<pixel_t> ovr_buffer;
  std::vector<pixel_t> src_buffer;
  std::vector<pixel_t> mask_buffer;

  // precalculate mask buffer et al.
  if constexpr (maskMode != MASK444) {
    ovr_buffer.resize(cwidth);
    src_buffer.resize(cwidth);
    if constexpr (has_alpha)
      mask_buffer.resize(cwidth);
  }

  using calc_t = typename std::conditional < lessthan16bits, int, int64_t>::type; // for non-overflowing 16 bit alpha_mul
  if constexpr (sizeof(pixel_t) == 1)
    bits_per_pixel = 8; // make quasi constexpr
  else if constexpr (sizeof(pixel_t) == 2 && !lessthan16bits)
    bits_per_pixel = 16; // make quasi constexpr

  const int rounder = 1 << (bits_per_pixel - 1);

  // for subsampled color spaces first do chroma, luma is only used for decision
  // second pass will do luma only
  for (int y = 0; y < cheight; ++y) {

    // Prepare all three pointers using the helper
    const pixel_t* ovr_ptr = prepare_effective_mask_for_row<maskMode, pixel_t>(ovrp, overlay_pitch, cwidth, ovr_buffer);
    const pixel_t* src_ptr = prepare_effective_mask_for_row<maskMode, pixel_t>(dstp, dst_pitch, cwidth, src_buffer);
    const pixel_t* effective_mask_ptr = nullptr;
    if constexpr (has_alpha) {
      effective_mask_ptr = prepare_effective_mask_for_row<maskMode, pixel_t>(maskp, mask_pitch, cwidth, mask_buffer);
    }

    for (int x = 0; x < cwidth; ++x) {
      int ovr = ovr_ptr[x];
      int src = src_ptr[x];
      int effective_mask = has_alpha ? effective_mask_ptr[x] : 0;

      const int alpha = has_alpha ? (int)(((calc_t)effective_mask * level + 1) >> bits_per_pixel) : level;

      int alpha_mask;
      if constexpr (mode == LIGHTEN)
        alpha_mask = ovr > (src + thresh) ? alpha : 0; // YUY2 was wrong: alpha_mask = (thresh + ovr) > src ? level : 0;
      else // DARKEN
        alpha_mask = ovr < (src - thresh) ? alpha : 0; // YUY2 was wrong: alpha_mask = (thresh + src) > ovr ? level : 0;

      if constexpr (!lumaonly)
      {
        // chroma u,v
        dstp_u[x] = dstp_u[x] + (int)(((calc_t)(ovrp_u[x] - dstp_u[x]) * alpha_mask + rounder) >> bits_per_pixel);
        dstp_v[x] = dstp_v[x] + (int)(((calc_t)(ovrp_v[x] - dstp_v[x]) * alpha_mask + rounder) >> bits_per_pixel);
      }

      // for 444: update here, width/height is the same as for chroma
      if constexpr (maskMode == MASK444)
        dstp[x] = dstp[x] + (int)(((calc_t)(ovrp[x] - dstp[x]) * alpha_mask + rounder) >> bits_per_pixel);
    }
    if constexpr (maskMode == MASK420 || maskMode == MASK420_MPEG2) {
      dstp += dst_pitch * 2; // skip vertical subsampling
      ovrp += overlay_pitch * 2;
      if constexpr (has_alpha) {
        //dstp_a += dst_pitch * 2;
        maskp += mask_pitch * 2;
      }
    }
    else {
      dstp += dst_pitch;
      ovrp += overlay_pitch;
      if constexpr (has_alpha) {
        //dstp_a += dst_pitch;
        maskp += mask_pitch;
      }
    }

    if constexpr (!lumaonly) {
      dstp_u += dst_pitchUV;
      dstp_v += dst_pitchUV;
      ovrp_u += overlay_pitchUV;
      ovrp_v += overlay_pitchUV;
    }

  }

  dst_pitch *= sizeof(pixel_t);
  dst_pitchUV *= sizeof(pixel_t);
  overlay_pitch *= sizeof(pixel_t);
  overlay_pitchUV *= sizeof(pixel_t);
  mask_pitch *= sizeof(pixel_t);

  // make luma
  if constexpr (!lumaonly && maskMode != MASK444)
    layer_yuv_lighten_darken_c<mode, MASK444, pixel_t, lessthan16bits, true /* lumaonly*/, has_alpha>(
      dstp8, dstp8_u, dstp8_v, //dstp8_a,
      ovrp8, ovrp8_u, ovrp8_v, maskp8,
      dst_pitch, dst_pitchUV, overlay_pitch, overlay_pitchUV, mask_pitch,
      width, height, level, thresh, bits_per_pixel);
}

template<int mode, MaskMode maskMode, bool lumaonly, bool has_alpha>
static void layer_yuv_lighten_darken_f_c(
  BYTE* dstp8, BYTE* dstp8_u, BYTE* dstp8_v /*, BYTE* dstp8_a*/,
  const BYTE* ovrp8, const BYTE* ovrp8_u, const BYTE* ovrp8_v, const BYTE* maskp8,
  int dst_pitch, int dst_pitchUV,
  int overlay_pitch, int overlay_pitchUV,
  int mask_pitch,
  int width, int height, float opacity, float thresh) {

  float* dstp = reinterpret_cast<float*>(dstp8);
  float* dstp_u = reinterpret_cast<float*>(dstp8_u);
  float* dstp_v = reinterpret_cast<float*>(dstp8_v);
  //float* dstp_a = reinterpret_cast<float *>(dstp8_a);

  const float* ovrp = reinterpret_cast<const float*>(ovrp8);
  const float* ovrp_u = reinterpret_cast<const float*>(ovrp8_u);
  const float* ovrp_v = reinterpret_cast<const float*>(ovrp8_v);
  const float* maskp = reinterpret_cast<const float*>(maskp8);

  dst_pitch /= sizeof(float);
  dst_pitchUV /= sizeof(float);
  overlay_pitch /= sizeof(float);
  overlay_pitchUV /= sizeof(float);
  mask_pitch /= sizeof(float);

  const int cwidth = (maskMode == MASK444) ? width : (maskMode == MASK411) ? width >> 2 : width >> 1; // 444:/1  420,422:/2  411:/4
  const int cheight = (maskMode == MASK444 || maskMode == MASK422 || maskMode == MASK422_MPEG2 || maskMode == MASK411) ? height : height >> 1; // 444,422,411:/1  420:/2

  // In lighten/darken we need 3 buffers:
  std::vector<float> ovr_buffer;
  std::vector<float> src_buffer;
  std::vector<float> mask_buffer;

  // precalculate mask buffer et al.
  if constexpr (maskMode != MASK444) {
    ovr_buffer.resize(cwidth);
    src_buffer.resize(cwidth);
    if constexpr (has_alpha)
      mask_buffer.resize(cwidth);
  }

  // for subsampled color spaces first do chroma, because luma is used for decision
  // second pass will do luma only
  for (int y = 0; y < cheight; ++y) {
    const float* ovr_ptr = prepare_effective_mask_for_row_f<maskMode>(ovrp, overlay_pitch, cwidth, ovr_buffer);
    const float* src_ptr = prepare_effective_mask_for_row_f<maskMode>(dstp, dst_pitch, cwidth, src_buffer);
    const float* effective_mask_ptr = nullptr;
    if constexpr (has_alpha) {
      effective_mask_ptr = prepare_effective_mask_for_row_f<maskMode>(maskp, mask_pitch, cwidth, mask_buffer);
    }

    for (int x = 0; x < cwidth; ++x) {
      float ovr = ovr_ptr[x];
      float src = src_ptr[x];
      float effective_mask = has_alpha ? effective_mask_ptr[x] : 0;

      const float alpha = has_alpha ? effective_mask * opacity : opacity;

      float alpha_mask;
      if constexpr (mode == LIGHTEN)
        alpha_mask = ovr > (src + thresh) ? alpha : 0; // YUY2 was wrong: alpha_mask = (thresh + ovr) > src ? level : 0;
      else // DARKEN
        alpha_mask = ovr < (src - thresh) ? alpha : 0; // YUY2 was wrong: alpha_mask = (thresh + src) > ovr ? level : 0;

      if constexpr (!lumaonly)
      {
        // chroma u,v
        dstp_u[x] = dstp_u[x] + (ovrp_u[x] - dstp_u[x]) * alpha_mask;
        dstp_v[x] = dstp_v[x] + (ovrp_v[x] - dstp_v[x]) * alpha_mask;
        //dstp_a[x] = dstp_a[x] + (maskp[x] - dstp_a[x]) * alpha_mask;
      }

      // for 444: update here, width/height is the same as for chroma
      if constexpr (maskMode == MASK444)
        dstp[x] = dstp[x] + (ovrp[x] - dstp[x]) * alpha_mask;
    }
    if constexpr (maskMode == MASK420 || maskMode == MASK420_MPEG2) {
      dstp += dst_pitch * 2; // skip vertical subsampling
      ovrp += overlay_pitch * 2;
      if constexpr (has_alpha) {
        //dstp_a += dst_pitch * 2;
        maskp += mask_pitch * 2;
      }
    }
    else {
      dstp += dst_pitch;
      ovrp += overlay_pitch;
      if constexpr (has_alpha) {
        //dstp_a += dst_pitch;
        maskp += mask_pitch;
      }
    }

    if constexpr (!lumaonly) {
      dstp_u += dst_pitchUV;
      dstp_v += dst_pitchUV;
      ovrp_u += overlay_pitchUV;
      ovrp_v += overlay_pitchUV;
    }
  }

  dst_pitch *= sizeof(float);
  dst_pitchUV *= sizeof(float);
  overlay_pitch *= sizeof(float);
  overlay_pitchUV *= sizeof(float);
  mask_pitch *= sizeof(float);

  // make luma
  if constexpr (!lumaonly && maskMode != MASK444)
    layer_yuv_lighten_darken_f_c<mode, MASK444, true /* lumaonly*/, has_alpha>(
      dstp8, dstp8_u, dstp8_v, //dstp8_a,
      ovrp8, ovrp8_u, ovrp8_v, maskp8,
      dst_pitch, dst_pitchUV, overlay_pitch, overlay_pitchUV, mask_pitch,
      width, height, opacity, thresh);
}

DISABLE_WARNING_POP

// Dispatchers

static void get_layer_yuv_lighten_darken_functions(bool isLighten, int placement, VideoInfo& vi, int bits_per_pixel, /*out*/layer_yuv_lighten_darken_c_t** layer_fn, /*out*/layer_yuv_lighten_darken_f_c_t** layer_f_fn) {

#define YUV_LIGHTEN_DARKEN_DISPATCH(L_or_D, MaskType, lumaonly, has_alpha) \
      { if (bits_per_pixel == 8) \
        *layer_fn = layer_yuv_lighten_darken_c<L_or_D, MaskType, uint8_t, true /*lessthan16bits*/, lumaonly /*lumaonly*/, has_alpha /*has_alpha*/>; \
      else if (bits_per_pixel < 16) \
        *layer_fn = layer_yuv_lighten_darken_c<L_or_D, MaskType, uint16_t, true/*lessthan16bits*/, lumaonly /*lumaonly*/, has_alpha /*has_alpha*/>; \
      else if (bits_per_pixel == 16) \
        *layer_fn = layer_yuv_lighten_darken_c<L_or_D, MaskType, uint16_t, false/*lessthan16bits*/, lumaonly /*lumaonly*/, has_alpha /*has_alpha*/>; \
      else /* float */ \
        *layer_f_fn = layer_yuv_lighten_darken_f_c<L_or_D, MaskType, lumaonly /*lumaonly*/, has_alpha /*has_alpha*/>; \
}

  if (isLighten) {

    if (vi.IsYV411())
      *layer_fn = layer_yuv_lighten_darken_c<LIGHTEN, MASK411, uint8_t, true /*lessthan16bits*/, false /*lumaonly*/, false /*has_alpha*/>;
    else if (vi.Is420())
    {
      if (placement == PLACEMENT_MPEG1)
        YUV_LIGHTEN_DARKEN_DISPATCH(LIGHTEN, MASK420, false, false)
      else
        YUV_LIGHTEN_DARKEN_DISPATCH(LIGHTEN, MASK420_MPEG2, false, false)
        // PLACEMENT_MPEG2
    }
    else if (vi.Is422())
    {
      if (placement == PLACEMENT_MPEG1)
        YUV_LIGHTEN_DARKEN_DISPATCH(LIGHTEN, MASK422, false, false)
      else
        YUV_LIGHTEN_DARKEN_DISPATCH(LIGHTEN, MASK422_MPEG2, false, false)
        // PLACEMENT_MPEG2
    }
    else if (vi.Is444())
      YUV_LIGHTEN_DARKEN_DISPATCH(LIGHTEN, MASK444, false, false)
    else if (vi.IsY())
      YUV_LIGHTEN_DARKEN_DISPATCH(LIGHTEN, MASK444, true, false)
  }
  else {
    // darken
    if (vi.IsYV411())
      *layer_fn = layer_yuv_lighten_darken_c<DARKEN, MASK411, uint8_t, true /*lessthan16bits*/, false /*lumaonly*/, false /*has_alpha*/>;
    else if (vi.Is420())
    {
      if (placement == PLACEMENT_MPEG1)
        YUV_LIGHTEN_DARKEN_DISPATCH(DARKEN, MASK420, false, false)
      else // PLACEMENT_MPEG2
        YUV_LIGHTEN_DARKEN_DISPATCH(DARKEN, MASK420_MPEG2, false, false)
    }
    else if (vi.Is422())
    {
      if (placement == PLACEMENT_MPEG1)
        YUV_LIGHTEN_DARKEN_DISPATCH(DARKEN, MASK422, false, false)
      else // PLACEMENT_MPEG2
        YUV_LIGHTEN_DARKEN_DISPATCH(DARKEN, MASK422_MPEG2, false, false)
    }
    else if (vi.Is444())
      YUV_LIGHTEN_DARKEN_DISPATCH(DARKEN, MASK444, false, false)
    else if (vi.IsY())
      YUV_LIGHTEN_DARKEN_DISPATCH(DARKEN, MASK444, true, false)
  }
#undef YUV_LIGHTEN_DARKEN_DISPATCH
}


static void get_layer_yuv_mul_functions(
  bool is_chroma, bool use_chroma, bool hasAlpha,
  int placement, VideoInfo& vi, int bits_per_pixel,
  /*out*/layer_yuv_mul_c_t** layer_fn,
  /*out*/layer_yuv_mul_f_c_t** layer_f_fn)
{
#define YUV_MUL_DISPATCH(MaskType, is_chroma, use_chroma, has_alpha) \
  { if (bits_per_pixel == 8) \
    *layer_fn = layer_yuv_mul_c<MaskType, uint8_t, true /*lessthan16bits*/, is_chroma, use_chroma, has_alpha>; \
  else if (bits_per_pixel < 16) \
    *layer_fn = layer_yuv_mul_c<MaskType, uint16_t, true /*lessthan16bits*/, is_chroma, use_chroma, has_alpha>; \
  else if (bits_per_pixel == 16) \
    *layer_fn = layer_yuv_mul_c<MaskType, uint16_t, false /*lessthan16bits*/, is_chroma, use_chroma, has_alpha>; \
  else /* float */ \
    *layer_f_fn = layer_yuv_mul_f_c<MaskType, is_chroma, use_chroma, has_alpha>; \
  }

  if (is_chroma) // not luma channel
  {
    if (vi.IsYV411())
    {
      if (use_chroma)
        *layer_fn = layer_yuv_mul_c<MASK411, uint8_t, true /*lessthan16bits*/, true, true, false>;
      else
        *layer_fn = layer_yuv_mul_c<MASK411, uint8_t, true /*lessthan16bits*/, true, false, false>;
    }
    else if (vi.Is420())
    {
      if (placement == PLACEMENT_MPEG1) {
        if (hasAlpha) {
          if (use_chroma) YUV_MUL_DISPATCH(MASK420, true, true, true)
          else YUV_MUL_DISPATCH(MASK420, true, false, true)
        }
        else {
          if (use_chroma) YUV_MUL_DISPATCH(MASK420, true, true, false)
          else YUV_MUL_DISPATCH(MASK420, true, false, false)
        }
      }
      else {
        if (hasAlpha) {
          if (use_chroma) YUV_MUL_DISPATCH(MASK420_MPEG2, true, true, true)
          else YUV_MUL_DISPATCH(MASK420_MPEG2, true, false, true)
        }
        else {
          if (use_chroma) YUV_MUL_DISPATCH(MASK420_MPEG2, true, true, false)
          else YUV_MUL_DISPATCH(MASK420_MPEG2, true, false, false)
        }
      }
    }
    else if (vi.Is422())
    {
      if (placement == PLACEMENT_MPEG1) {
        if (hasAlpha) {
          if (use_chroma) YUV_MUL_DISPATCH(MASK422, true, true, true)
          else YUV_MUL_DISPATCH(MASK422, true, false, true)
        }
        else {
          if (use_chroma) YUV_MUL_DISPATCH(MASK422, true, true, false)
          else YUV_MUL_DISPATCH(MASK422, true, false, false)
        }
      }
      else {
        if (hasAlpha) {
          if (use_chroma) YUV_MUL_DISPATCH(MASK422_MPEG2, true, true, true)
          else YUV_MUL_DISPATCH(MASK422_MPEG2, true, false, true)
        }
        else {
          if (use_chroma) YUV_MUL_DISPATCH(MASK422_MPEG2, true, true, false)
          else YUV_MUL_DISPATCH(MASK422_MPEG2, true, false, false)
        }
      }
    }
    else if (vi.Is444())
    {
      if (hasAlpha) {
        if (use_chroma) YUV_MUL_DISPATCH(MASK444, true, true, true)
        else YUV_MUL_DISPATCH(MASK444, true, false, true)
      }
      else {
        if (use_chroma) YUV_MUL_DISPATCH(MASK444, true, true, false)
        else YUV_MUL_DISPATCH(MASK444, true, false, false)
      }
    }
  }
  else // luma channel
  {
    if (hasAlpha)
      YUV_MUL_DISPATCH(MASK444, false, false, true)
    else
      YUV_MUL_DISPATCH(MASK444, false, false, false)
  }
#undef YUV_MUL_DISPATCH
}

template<bool is_subtract>
static void get_layer_yuv_add_subtract_functions(
  bool is_chroma, bool use_chroma, bool hasAlpha,
  int placement, VideoInfo& vi, int bits_per_pixel,
  /*out*/layer_yuv_add_subtract_c_t** layer_fn,
  /*out*/layer_yuv_add_subtract_f_c_t** layer_f_fn)
{
#define YUV_ADD_SUBTRACT_DISPATCH(MaskType, is_chroma, use_chroma, has_alpha) \
  { if (bits_per_pixel == 8) \
    *layer_fn = layer_yuv_add_subtract_c<MaskType, uint8_t, true /*lessthan16bits*/, is_chroma, use_chroma, has_alpha, is_subtract>; \
  else if (bits_per_pixel < 16) \
    *layer_fn = layer_yuv_add_subtract_c<MaskType, uint16_t, true /*lessthan16bits*/, is_chroma, use_chroma, has_alpha, is_subtract>; \
  else if (bits_per_pixel == 16) \
    *layer_fn = layer_yuv_add_subtract_c<MaskType, uint16_t, false /*lessthan16bits*/, is_chroma, use_chroma, has_alpha, is_subtract>; \
  else /* float */ \
    *layer_f_fn = layer_yuv_add_subtract_f_c<MaskType, is_chroma, use_chroma, has_alpha, is_subtract>; \
  }

  if (is_chroma) // not luma channel
  {
    if (vi.IsYV411())
    {
      if (use_chroma)
        *layer_fn = layer_yuv_add_subtract_c<MASK411, uint8_t, true /*lessthan16bits*/, true, true, false, is_subtract>;
      else
        *layer_fn = layer_yuv_add_subtract_c<MASK411, uint8_t, true /*lessthan16bits*/, true, false, false, is_subtract>;
    }
    else if (vi.Is420())
    {
      if (placement == PLACEMENT_MPEG1) {
        if (hasAlpha) {
          if (use_chroma) YUV_ADD_SUBTRACT_DISPATCH(MASK420, true, true, true)
          else YUV_ADD_SUBTRACT_DISPATCH(MASK420, true, false, true)
        }
        else {
          if (use_chroma) YUV_ADD_SUBTRACT_DISPATCH(MASK420, true, true, false)
          else YUV_ADD_SUBTRACT_DISPATCH(MASK420, true, false, false)
        }
      }
      else {
        if (hasAlpha) {
          if (use_chroma) YUV_ADD_SUBTRACT_DISPATCH(MASK420_MPEG2, true, true, true)
          else YUV_ADD_SUBTRACT_DISPATCH(MASK420_MPEG2, true, false, true)
        }
        else {
          if (use_chroma) YUV_ADD_SUBTRACT_DISPATCH(MASK420_MPEG2, true, true, false)
          else YUV_ADD_SUBTRACT_DISPATCH(MASK420_MPEG2, true, false, false)
        }
      }
    }
    else if (vi.Is422())
    {
      if (placement == PLACEMENT_MPEG1) {
        if (hasAlpha) {
          if (use_chroma) YUV_ADD_SUBTRACT_DISPATCH(MASK422, true, true, true)
          else YUV_ADD_SUBTRACT_DISPATCH(MASK422, true, false, true)
        }
        else {
          if (use_chroma) YUV_ADD_SUBTRACT_DISPATCH(MASK422, true, true, false)
          else YUV_ADD_SUBTRACT_DISPATCH(MASK422, true, false, false)
        }
      }
      else {
        if (hasAlpha) {
          if (use_chroma) YUV_ADD_SUBTRACT_DISPATCH(MASK422_MPEG2, true, true, true)
          else YUV_ADD_SUBTRACT_DISPATCH(MASK422_MPEG2, true, false, true)
        }
        else {
          if (use_chroma) YUV_ADD_SUBTRACT_DISPATCH(MASK422_MPEG2, true, true, false)
          else YUV_ADD_SUBTRACT_DISPATCH(MASK422_MPEG2, true, false, false)
        }
      }
    }
    else if (vi.Is444())
    {
      if (hasAlpha) {
        if (use_chroma) YUV_ADD_SUBTRACT_DISPATCH(MASK444, true, true, true)
        else YUV_ADD_SUBTRACT_DISPATCH(MASK444, true, false, true)
      }
      else {
        if (use_chroma) YUV_ADD_SUBTRACT_DISPATCH(MASK444, true, true, false)
        else YUV_ADD_SUBTRACT_DISPATCH(MASK444, true, false, false)
      }
    }
  }
  else // luma channel
  {
    if (hasAlpha)
      YUV_ADD_SUBTRACT_DISPATCH(MASK444, false, false, true)
    else
      YUV_ADD_SUBTRACT_DISPATCH(MASK444, false, false, false)
  }
#undef YUV_ADD_SUBTRACT_DISPATCH
}

/* planar rgb */

template<int mode, typename pixel_t, bool lessthan16bits, bool has_alpha>
static void layer_planarrgb_lighten_darken_c(BYTE** dstp8, const BYTE** ovrp8, int dst_pitch, int overlay_pitch, int width, int height, int level, int thresh, int bits_per_pixel) {
  pixel_t* dstp_g = reinterpret_cast<pixel_t*>(dstp8[0]);
  pixel_t* dstp_b = reinterpret_cast<pixel_t*>(dstp8[1]);
  pixel_t* dstp_r = reinterpret_cast<pixel_t*>(dstp8[2]);
  pixel_t* dstp_a = reinterpret_cast<pixel_t*>(dstp8[3]);
  const pixel_t* ovrp_g = reinterpret_cast<const pixel_t*>(ovrp8[0]);
  const pixel_t* ovrp_b = reinterpret_cast<const pixel_t*>(ovrp8[1]);
  const pixel_t* ovrp_r = reinterpret_cast<const pixel_t*>(ovrp8[2]);
  const pixel_t* maskp = reinterpret_cast<const pixel_t*>(ovrp8[3]);

  dst_pitch /= sizeof(pixel_t);
  overlay_pitch /= sizeof(pixel_t);

  using calc_t = typename std::conditional < lessthan16bits, int, int64_t>::type; // for non-overflowing 16 bit alpha_mul
  if constexpr (sizeof(pixel_t) == 1)
    bits_per_pixel = 8; // make quasi constexpr
  else if constexpr (sizeof(pixel_t) == 2 && !lessthan16bits)
    bits_per_pixel = 16; // make quasi constexpr

  const int rounder = 1 << (bits_per_pixel - 1);

  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      calc_t alpha = has_alpha ? ((calc_t)maskp[x] * level + 1) >> bits_per_pixel : level;

      calc_t luma_ovr = (cyb * ovrp_b[x] + cyg * ovrp_g[x] + cyr * ovrp_r[x]) >> 15; // no rounding, not really needed here
      calc_t luma_src = (cyb * dstp_b[x] + cyg * dstp_g[x] + cyr * dstp_r[x]) >> 15;

      if constexpr (mode == LIGHTEN)
        alpha = luma_ovr > luma_src + thresh ? alpha : 0;
      else // DARKEN
        alpha = luma_ovr < luma_src - thresh ? alpha : 0;

      dstp_r[x] = (pixel_t)(dstp_r[x] + (((ovrp_r[x] - dstp_r[x]) * alpha + rounder) >> bits_per_pixel));
      dstp_g[x] = (pixel_t)(dstp_g[x] + (((ovrp_g[x] - dstp_g[x]) * alpha + rounder) >> bits_per_pixel));
      dstp_b[x] = (pixel_t)(dstp_b[x] + (((ovrp_b[x] - dstp_b[x]) * alpha + rounder) >> bits_per_pixel));
      if constexpr (has_alpha)
        dstp_a[x] = (pixel_t)(dstp_a[x] + (((maskp[x] - dstp_a[x]) * alpha + rounder) >> bits_per_pixel));
    }
    dstp_g += dst_pitch;
    dstp_b += dst_pitch;
    dstp_r += dst_pitch;
    if constexpr (has_alpha)
      dstp_a += dst_pitch;
    ovrp_g += overlay_pitch;
    ovrp_b += overlay_pitch;
    ovrp_r += overlay_pitch;
    if constexpr (has_alpha)
      maskp += overlay_pitch;
  }
}

template<int mode, bool has_alpha>
static void layer_planarrgb_lighten_darken_f_c(BYTE** dstp8, const BYTE** ovrp8, int dst_pitch, int overlay_pitch, int width, int height, float opacity, float thresh) {
  float* dstp_g = reinterpret_cast<float*>(dstp8[0]);
  float* dstp_b = reinterpret_cast<float*>(dstp8[1]);
  float* dstp_r = reinterpret_cast<float*>(dstp8[2]);
  float* dstp_a = reinterpret_cast<float*>(dstp8[3]);
  const float* ovrp_g = reinterpret_cast<const float*>(ovrp8[0]);
  const float* ovrp_b = reinterpret_cast<const float*>(ovrp8[1]);
  const float* ovrp_r = reinterpret_cast<const float*>(ovrp8[2]);
  const float* maskp = reinterpret_cast<const float*>(ovrp8[3]);

  dst_pitch /= sizeof(float);
  overlay_pitch /= sizeof(float);

  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      float alpha = has_alpha ? maskp[x] * opacity : opacity;

      float luma_ovr = cyb_f * ovrp_b[x] + cyg_f * ovrp_g[x] + cyr_f * ovrp_r[x];
      float luma_src = cyb_f * dstp_b[x] + cyg_f * dstp_g[x] + cyr_f * dstp_r[x];

      if constexpr (mode == LIGHTEN)
        alpha = luma_ovr > luma_src + thresh ? alpha : 0;
      else // DARKEN
        alpha = luma_ovr < luma_src - thresh ? alpha : 0;

      dstp_r[x] = dstp_r[x] + (ovrp_r[x] - dstp_r[x]) * alpha;
      dstp_g[x] = dstp_g[x] + (ovrp_g[x] - dstp_g[x]) * alpha;
      dstp_b[x] = dstp_b[x] + (ovrp_b[x] - dstp_b[x]) * alpha;
      if constexpr (has_alpha)
        dstp_a[x] = dstp_a[x] + (maskp[x] - dstp_a[x]) * alpha;
    }
    dstp_g += dst_pitch;
    dstp_b += dst_pitch;
    dstp_r += dst_pitch;
    if constexpr (has_alpha)
      dstp_a += dst_pitch;
    ovrp_g += overlay_pitch;
    ovrp_b += overlay_pitch;
    ovrp_r += overlay_pitch;
    if constexpr (has_alpha)
      maskp += overlay_pitch;
  }
}


static void get_layer_planarrgb_lighten_darken_functions(bool isLighten, bool hasAlpha, int bits_per_pixel, /*out*/layer_planarrgb_lighten_darken_c_t** layer_fn, /*out*/layer_planarrgb_lighten_darken_f_c_t** layer_f_fn) {

#define PLANARRGB_LD_DISPATCH(LorD, has_alpha) \
      { if (bits_per_pixel == 8) \
        *layer_fn = layer_planarrgb_lighten_darken_c<LorD, uint8_t, true /*lessthan16bits*/, has_alpha>; \
      else if (bits_per_pixel < 16) \
        *layer_fn = layer_planarrgb_lighten_darken_c<LorD, uint16_t, true /*lessthan16bits*/, has_alpha>; \
      else if (bits_per_pixel == 16) \
        *layer_fn = layer_planarrgb_lighten_darken_c<LorD, uint16_t, false /*lessthan16bits*/, has_alpha>; \
      else /* float */ \
        *layer_f_fn = layer_planarrgb_lighten_darken_f_c<LorD, has_alpha>; \
      }

  if (isLighten) {
    if (hasAlpha) {
      PLANARRGB_LD_DISPATCH(LIGHTEN, true)
    }
    else {
      PLANARRGB_LD_DISPATCH(LIGHTEN, false)
    }
  } // lighten end
  else {
    if (hasAlpha) {
      PLANARRGB_LD_DISPATCH(DARKEN, true)
    }
    else {
      PLANARRGB_LD_DISPATCH(DARKEN, false)
    }
  }
#undef PLANARRGB_LD_DISPATCH
}

template<typename pixel_t, bool lessthan16bits, bool chroma, bool has_alpha, bool subtract>
static void layer_planarrgb_add_subtract_c(BYTE** dstp8, const BYTE** ovrp8, int dst_pitch, int overlay_pitch, int width, int height, int level, int bits_per_pixel) {
  pixel_t* dstp_g = reinterpret_cast<pixel_t*>(dstp8[0]);
  pixel_t* dstp_b = reinterpret_cast<pixel_t*>(dstp8[1]);
  pixel_t* dstp_r = reinterpret_cast<pixel_t*>(dstp8[2]);
  pixel_t* dstp_a = reinterpret_cast<pixel_t*>(dstp8[3]);
  const pixel_t* ovrp_g = reinterpret_cast<const pixel_t*>(ovrp8[0]);
  const pixel_t* ovrp_b = reinterpret_cast<const pixel_t*>(ovrp8[1]);
  const pixel_t* ovrp_r = reinterpret_cast<const pixel_t*>(ovrp8[2]);
  const pixel_t* maskp = reinterpret_cast<const pixel_t*>(ovrp8[3]);

  dst_pitch /= sizeof(pixel_t);
  overlay_pitch /= sizeof(pixel_t);

  using calc_t = typename std::conditional < lessthan16bits, int, int64_t>::type; // for non-overflowing 16 bit alpha_mul
  if constexpr (sizeof(pixel_t) == 1)
    bits_per_pixel = 8; // make quasi constexpr
  else if constexpr (sizeof(pixel_t) == 2 && !lessthan16bits)
    bits_per_pixel = 16; // make quasi constexpr

  const int max_pixel_value = (1 << bits_per_pixel) - 1;

  const int rounder = 1 << (bits_per_pixel - 1);

  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      calc_t alpha = has_alpha ? ((calc_t)maskp[x] * level + 1) >> bits_per_pixel : level;
      if constexpr (subtract) {
        // subtract
        if constexpr (chroma) {
          dstp_r[x] = (pixel_t)(dstp_r[x] + (((max_pixel_value - ovrp_r[x] - dstp_r[x]) * alpha + rounder) >> bits_per_pixel));
          dstp_g[x] = (pixel_t)(dstp_g[x] + (((max_pixel_value - ovrp_g[x] - dstp_g[x]) * alpha + rounder) >> bits_per_pixel));
          dstp_b[x] = (pixel_t)(dstp_b[x] + (((max_pixel_value - ovrp_b[x] - dstp_b[x]) * alpha + rounder) >> bits_per_pixel));
          if constexpr (has_alpha) // fixme: to be decided. YUV does not update target alpha, rgb32 does
            dstp_a[x] = (pixel_t)(dstp_a[x] + (((max_pixel_value - maskp[x] - dstp_a[x]) * alpha + rounder) >> bits_per_pixel));
        }
        else { // use luma instead of overlay
          calc_t luma = (cyb * (max_pixel_value - ovrp_b[x]) + cyg * (max_pixel_value - ovrp_g[x]) + cyr * (max_pixel_value - ovrp_r[x])) >> 15; // no rounding not really needed here

          dstp_r[x] = (pixel_t)(dstp_r[x] + (((luma - dstp_r[x]) * alpha + rounder) >> bits_per_pixel));
          dstp_g[x] = (pixel_t)(dstp_g[x] + (((luma - dstp_g[x]) * alpha + rounder) >> bits_per_pixel));
          dstp_b[x] = (pixel_t)(dstp_b[x] + (((luma - dstp_b[x]) * alpha + rounder) >> bits_per_pixel));
          if constexpr (has_alpha)
            dstp_a[x] = (pixel_t)(dstp_a[x] + (((luma - dstp_a[x]) * alpha + rounder) >> bits_per_pixel));
        }
      }
      else {
        // add
        if constexpr (chroma) {
          dstp_r[x] = (pixel_t)(dstp_r[x] + (((ovrp_r[x] - dstp_r[x]) * alpha + rounder) >> bits_per_pixel));
          dstp_g[x] = (pixel_t)(dstp_g[x] + (((ovrp_g[x] - dstp_g[x]) * alpha + rounder) >> bits_per_pixel));
          dstp_b[x] = (pixel_t)(dstp_b[x] + (((ovrp_b[x] - dstp_b[x]) * alpha + rounder) >> bits_per_pixel));
          if constexpr (has_alpha)
            dstp_a[x] = (pixel_t)(dstp_a[x] + (((maskp[x] - dstp_a[x]) * alpha + rounder) >> bits_per_pixel));
        }
        else { // use luma instead of overlay
          calc_t luma = (cyb * ovrp_b[x] + cyg * ovrp_g[x] + cyr * ovrp_r[x]) >> 15; // no rounding not really needed here

          dstp_r[x] = (pixel_t)(dstp_r[x] + (((luma - dstp_r[x]) * alpha + rounder) >> bits_per_pixel));
          dstp_g[x] = (pixel_t)(dstp_g[x] + (((luma - dstp_g[x]) * alpha + rounder) >> bits_per_pixel));
          dstp_b[x] = (pixel_t)(dstp_b[x] + (((luma - dstp_b[x]) * alpha + rounder) >> bits_per_pixel));
          if constexpr (has_alpha)
            dstp_a[x] = (pixel_t)(dstp_a[x] + (((luma - dstp_a[x]) * alpha + rounder) >> bits_per_pixel));
        }
      }
    }
    dstp_g += dst_pitch;
    dstp_b += dst_pitch;
    dstp_r += dst_pitch;
    if constexpr (has_alpha)
      dstp_a += dst_pitch;
    ovrp_g += overlay_pitch;
    ovrp_b += overlay_pitch;
    ovrp_r += overlay_pitch;
    if constexpr (has_alpha)
      maskp += overlay_pitch;
  }
}

template<bool chroma, bool has_alpha, bool subtract>
static void layer_planarrgb_add_subtract_f_c(BYTE** dstp8, const BYTE** ovrp8, int dst_pitch, int overlay_pitch, int width, int height, float opacity) {
  float* dstp_g = reinterpret_cast<float*>(dstp8[0]);
  float* dstp_b = reinterpret_cast<float*>(dstp8[1]);
  float* dstp_r = reinterpret_cast<float*>(dstp8[2]);
  float* dstp_a = reinterpret_cast<float*>(dstp8[3]);
  const float* ovrp_g = reinterpret_cast<const float*>(ovrp8[0]);
  const float* ovrp_b = reinterpret_cast<const float*>(ovrp8[1]);
  const float* ovrp_r = reinterpret_cast<const float*>(ovrp8[2]);
  const float* maskp = reinterpret_cast<const float*>(ovrp8[3]);

  dst_pitch /= sizeof(float);
  overlay_pitch /= sizeof(float);

  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      float alpha = has_alpha ? maskp[x] * opacity : opacity;

      if constexpr (subtract) {
        // subtract
        if constexpr (chroma) {
          dstp_r[x] = dstp_r[x] + (1.0f - ovrp_r[x] - dstp_r[x]) * alpha;
          dstp_g[x] = dstp_g[x] + (1.0f - ovrp_g[x] - dstp_g[x]) * alpha;
          dstp_b[x] = dstp_b[x] + (1.0f - ovrp_b[x] - dstp_b[x]) * alpha;
          if constexpr (has_alpha)
            dstp_a[x] = dstp_a[x] + (1.0f - maskp[x] - dstp_a[x]) * alpha;
        }
        else { // use luma instead of overlay
          float luma = cyb_f * (1.0f - ovrp_b[x]) + cyg_f * (1.0f - ovrp_g[x]) + cyr_f * (1.0f - ovrp_r[x]);
          dstp_r[x] = dstp_r[x] + (luma - dstp_r[x]) * alpha;
          dstp_g[x] = dstp_g[x] + (luma - dstp_g[x]) * alpha;
          dstp_b[x] = dstp_b[x] + (luma - dstp_b[x]) * alpha;
          if constexpr (has_alpha)
            dstp_a[x] = dstp_a[x] + (luma * dstp_a[x] - dstp_a[x]) * alpha;
        }
      }
      else {
        // add
        if constexpr (chroma) {
          dstp_r[x] = dstp_r[x] + (ovrp_r[x] - dstp_r[x]) * alpha;
          dstp_g[x] = dstp_g[x] + (ovrp_g[x] - dstp_g[x]) * alpha;
          dstp_b[x] = dstp_b[x] + (ovrp_b[x] - dstp_b[x]) * alpha;
          if constexpr (has_alpha)
            dstp_a[x] = dstp_a[x] + (maskp[x] - dstp_a[x]) * alpha;
        }
        else { // use luma instead of overlay
          float luma = cyb_f * ovrp_b[x] + cyg_f * ovrp_g[x] + cyr_f * ovrp_r[x];
          dstp_r[x] = dstp_r[x] + (luma - dstp_r[x]) * alpha;
          dstp_g[x] = dstp_g[x] + (luma - dstp_g[x]) * alpha;
          dstp_b[x] = dstp_b[x] + (luma - dstp_b[x]) * alpha;
          if constexpr (has_alpha)
            dstp_a[x] = dstp_a[x] + (luma * dstp_a[x] - dstp_a[x]) * alpha;
        }
      }
    }
    dstp_g += dst_pitch;
    dstp_b += dst_pitch;
    dstp_r += dst_pitch;
    if constexpr (has_alpha)
      dstp_a += dst_pitch;
    ovrp_g += overlay_pitch;
    ovrp_b += overlay_pitch;
    ovrp_r += overlay_pitch;
    if constexpr (has_alpha)
      maskp += overlay_pitch;
  }
}

template<bool is_subtract>
static void get_layer_planarrgb_add_subtract_functions(
  bool chroma, bool hasAlpha, int bits_per_pixel,
  /*out*/layer_planarrgb_add_subtract_c_t** layer_fn,
  /*out*/layer_planarrgb_add_subtract_f_c_t** layer_f_fn)
{
#define PLANARRGB_ADD_SUBTRACT_DISPATCH(chroma, has_alpha) \
  { if (bits_per_pixel == 8) \
    *layer_fn = layer_planarrgb_add_subtract_c<uint8_t, true /*lessthan16bits*/, chroma, has_alpha, is_subtract>; \
  else if (bits_per_pixel < 16) \
    *layer_fn = layer_planarrgb_add_subtract_c<uint16_t, true /*lessthan16bits*/, chroma, has_alpha, is_subtract>; \
  else if (bits_per_pixel == 16) \
    *layer_fn = layer_planarrgb_add_subtract_c<uint16_t, false /*lessthan16bits*/, chroma, has_alpha, is_subtract>; \
  else /* float */ \
    *layer_f_fn = layer_planarrgb_add_subtract_f_c<chroma, has_alpha, is_subtract>; \
  }

  if (hasAlpha) {
    // planar RGBA
    if (chroma) {
      PLANARRGB_ADD_SUBTRACT_DISPATCH(true, true)
    }
    else {
      PLANARRGB_ADD_SUBTRACT_DISPATCH(false, true)
    }
  }
  else {
    // planar RGB
    if (chroma) {
      PLANARRGB_ADD_SUBTRACT_DISPATCH(true, false)
    }
    else {
      PLANARRGB_ADD_SUBTRACT_DISPATCH(false, false)
    }
  }
#undef PLANARRGB_ADD_SUBTRACT_DISPATCH
}

template<typename pixel_t, bool lessthan16bits, bool chroma, bool has_alpha>
static void layer_planarrgb_mul_c(BYTE** dstp8, const BYTE** ovrp8, int dst_pitch, int overlay_pitch, int width, int height, int level, int bits_per_pixel) {
  pixel_t* dstp_g = reinterpret_cast<pixel_t*>(dstp8[0]);
  pixel_t* dstp_b = reinterpret_cast<pixel_t*>(dstp8[1]);
  pixel_t* dstp_r = reinterpret_cast<pixel_t*>(dstp8[2]);
  pixel_t* dstp_a = reinterpret_cast<pixel_t*>(dstp8[3]);
  const pixel_t* ovrp_g = reinterpret_cast<const pixel_t*>(ovrp8[0]);
  const pixel_t* ovrp_b = reinterpret_cast<const pixel_t*>(ovrp8[1]);
  const pixel_t* ovrp_r = reinterpret_cast<const pixel_t*>(ovrp8[2]);
  const pixel_t* maskp = reinterpret_cast<const pixel_t*>(ovrp8[3]);

  dst_pitch /= sizeof(pixel_t);
  overlay_pitch /= sizeof(pixel_t);

  using calc_t = typename std::conditional < lessthan16bits, int, int64_t>::type;
  if constexpr (sizeof(pixel_t) == 1)
    bits_per_pixel = 8; // make quasi constexpr
  else if constexpr (sizeof(pixel_t) == 2 && !lessthan16bits)
    bits_per_pixel = 16; // make quasi constexpr

  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      calc_t alpha = has_alpha ? ((calc_t)maskp[x] * level + 1) >> bits_per_pixel : level;

      if constexpr (chroma) {
        dstp_r[x] = (pixel_t)(dstp_r[x] + ((((((calc_t)ovrp_r[x] * dstp_r[x]) >> bits_per_pixel) - dstp_r[x]) * alpha) >> bits_per_pixel));
        dstp_g[x] = (pixel_t)(dstp_g[x] + ((((((calc_t)ovrp_g[x] * dstp_g[x]) >> bits_per_pixel) - dstp_g[x]) * alpha) >> bits_per_pixel));
        dstp_b[x] = (pixel_t)(dstp_b[x] + ((((((calc_t)ovrp_b[x] * dstp_b[x]) >> bits_per_pixel) - dstp_b[x]) * alpha) >> bits_per_pixel));
        if constexpr (has_alpha)
          dstp_a[x] = (pixel_t)(dstp_a[x] + ((((((calc_t)maskp[x] * dstp_a[x]) >> bits_per_pixel) - dstp_a[x]) * alpha) >> bits_per_pixel));
      }
      else { // use luma instead of overlay
        calc_t luma = (cyb * ovrp_b[x] + cyg * ovrp_g[x] + cyr * ovrp_r[x]) >> 15; // no rounding not really needed here

        dstp_r[x] = (pixel_t)(dstp_r[x] + (((((luma * dstp_r[x]) >> bits_per_pixel) - dstp_r[x]) * alpha) >> bits_per_pixel));
        dstp_g[x] = (pixel_t)(dstp_g[x] + (((((luma * dstp_g[x]) >> bits_per_pixel) - dstp_g[x]) * alpha) >> bits_per_pixel));
        dstp_b[x] = (pixel_t)(dstp_b[x] + (((((luma * dstp_b[x]) >> bits_per_pixel) - dstp_b[x]) * alpha) >> bits_per_pixel));
        if constexpr (has_alpha)
          dstp_a[x] = (pixel_t)(dstp_a[x] + (((((luma * dstp_a[x]) >> bits_per_pixel) - dstp_a[x]) * alpha) >> bits_per_pixel));
      }
    }
    dstp_g += dst_pitch;
    dstp_b += dst_pitch;
    dstp_r += dst_pitch;
    if constexpr (has_alpha)
      dstp_a += dst_pitch;
    ovrp_g += overlay_pitch;
    ovrp_b += overlay_pitch;
    ovrp_r += overlay_pitch;
    if constexpr (has_alpha)
      maskp += overlay_pitch;
  }
}

template<bool chroma, bool has_alpha>
static void layer_planarrgb_mul_f_c(BYTE** dstp8, const BYTE** ovrp8, int dst_pitch, int overlay_pitch, int width, int height, float opacity) {
  float* dstp_g = reinterpret_cast<float*>(dstp8[0]);
  float* dstp_b = reinterpret_cast<float*>(dstp8[1]);
  float* dstp_r = reinterpret_cast<float*>(dstp8[2]);
  float* dstp_a = reinterpret_cast<float*>(dstp8[3]);
  const float* ovrp_g = reinterpret_cast<const float*>(ovrp8[0]);
  const float* ovrp_b = reinterpret_cast<const float*>(ovrp8[1]);
  const float* ovrp_r = reinterpret_cast<const float*>(ovrp8[2]);
  const float* maskp = reinterpret_cast<const float*>(ovrp8[3]);

  dst_pitch /= sizeof(float);
  overlay_pitch /= sizeof(float);

  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      float alpha = has_alpha ? maskp[x] * opacity : opacity;

      if constexpr (chroma) {
        dstp_r[x] = dstp_r[x] + (ovrp_r[x] * dstp_r[x] - dstp_r[x]) * alpha;
        dstp_g[x] = dstp_g[x] + (ovrp_g[x] * dstp_g[x] - dstp_g[x]) * alpha;
        dstp_b[x] = dstp_b[x] + (ovrp_b[x] * dstp_b[x] - dstp_b[x]) * alpha;
        if constexpr (has_alpha)
          dstp_a[x] = dstp_a[x] + (maskp[x] * dstp_a[x] - dstp_a[x]) * alpha;
      }
      else { // use luma instead of overlay
        float luma = cyb_f * ovrp_b[x] + cyg_f * ovrp_g[x] + cyr_f * ovrp_r[x];
        dstp_r[x] = dstp_r[x] + (luma * dstp_r[x] - dstp_r[x]) * alpha;
        dstp_g[x] = dstp_g[x] + (luma * dstp_g[x] - dstp_g[x]) * alpha;
        dstp_b[x] = dstp_b[x] + (luma * dstp_b[x] - dstp_b[x]) * alpha;
        if constexpr (has_alpha)
          dstp_a[x] = dstp_a[x] + (luma * dstp_a[x] - dstp_a[x]) * alpha;
      }
    }
    dstp_g += dst_pitch;
    dstp_b += dst_pitch;
    dstp_r += dst_pitch;
    if constexpr (has_alpha)
      dstp_a += dst_pitch;
    ovrp_g += overlay_pitch;
    ovrp_b += overlay_pitch;
    ovrp_r += overlay_pitch;
    if constexpr (has_alpha)
      maskp += overlay_pitch;
  }
}

static void get_layer_planarrgb_mul_functions(
  bool chroma, bool hasAlpha, int bits_per_pixel,
  /*out*/layer_planarrgb_mul_c_t** layer_fn,
  /*out*/layer_planarrgb_mul_f_c_t** layer_f_fn)
{
#define PLANARRGB_MUL_DISPATCH(chroma, has_alpha) \
  { if (bits_per_pixel == 8) \
    *layer_fn = layer_planarrgb_mul_c<uint8_t, true /*lessthan16bits*/, chroma, has_alpha>; \
  else if (bits_per_pixel < 16) \
    *layer_fn = layer_planarrgb_mul_c<uint16_t, true /*lessthan16bits*/, chroma, has_alpha>; \
  else if (bits_per_pixel == 16) \
    *layer_fn = layer_planarrgb_mul_c<uint16_t, false /*lessthan16bits*/, chroma, has_alpha>; \
  else /* float */ \
    *layer_f_fn = layer_planarrgb_mul_f_c<chroma, has_alpha>; \
  }

  if (hasAlpha) {
    // planar RGBA
    if (chroma) {
      PLANARRGB_MUL_DISPATCH(true, true)
    }
    else {
      PLANARRGB_MUL_DISPATCH(false, true)
    }
  }
  else {
    // planar RGB
    if (chroma) {
      PLANARRGB_MUL_DISPATCH(true, false)
    }
    else {
      PLANARRGB_MUL_DISPATCH(false, false)
    }
  }
#undef PLANARRGB_MUL_DISPATCH
}

