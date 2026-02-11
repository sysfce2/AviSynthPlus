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

#include <avisynth.h>
#include "../blend_common.h"
#include "blend_common_neon.h"
#include <arm_neon.h>
#include <cstdint>

/********************************
 ********* Blend Opaque *********
 ** Use for Lighten and Darken **
 ********************************/

// avs_core/filters/overlay/aarch64/blend_common_neon.cpp

template<typename pixel_t>
static AVS_FORCEINLINE void Eightpixels_to_Eightfloats_neon(const pixel_t* src, float32x4_t& src_lo, float32x4_t& src_hi) {
  if constexpr (sizeof(pixel_t) == 1) {
    // Load 8 bytes, widen to 8x16, then to 8x32, then to float
    uint8x8_t v = vld1_u8(reinterpret_cast<const uint8_t*>(src));
    uint16x8_t v16 = vmovl_u8(v);
    uint32x4_t v32_lo = vmovl_u16(vget_low_u16(v16));
    uint32x4_t v32_hi = vmovl_u16(vget_high_u16(v16));
    src_lo = vcvtq_f32_u32(v32_lo);
    src_hi = vcvtq_f32_u32(v32_hi);
  } else {
    // 16-bit: load 8x16, widen to 8x32, then to float
    uint16x8_t v16 = vld1q_u16(reinterpret_cast<const uint16_t*>(src));
    uint32x4_t v32_lo = vmovl_u16(vget_low_u16(v16));
    uint32x4_t v32_hi = vmovl_u16(vget_high_u16(v16));
    src_lo = vcvtq_f32_u32(v32_lo);
    src_hi = vcvtq_f32_u32(v32_hi);
  }
}

template<typename pixel_t>
static AVS_FORCEINLINE void Store_Eightpixels_neon(pixel_t* dst, float32x4_t what_lo, float32x4_t what_hi, const float32x4_t rounder) {
  // Add rounder and convert to int32
  what_lo = vaddq_f32(what_lo, rounder);
  what_hi = vaddq_f32(what_hi, rounder);
  int32x4_t si32_lo = vcvtq_s32_f32(what_lo);
  int32x4_t si32_hi = vcvtq_s32_f32(what_hi);

  if constexpr (sizeof(pixel_t) == 1) {
    // Clamp to [0,255], narrow to 16, then to 8, and store 8 bytes
    uint16x4_t u16_lo = vqmovun_s32(si32_lo);
    uint16x4_t u16_hi = vqmovun_s32(si32_hi);
    uint8x8_t u8 = vqmovn_u16(vcombine_u16(u16_lo, u16_hi));
    vst1_u8(reinterpret_cast<uint8_t*>(dst), u8);
  } else {
    // Clamp to [0,65535], narrow to 16, and store 8x16
    uint16x4_t u16_lo = vqmovun_s32(si32_lo);
    uint16x4_t u16_hi = vqmovun_s32(si32_hi);
    uint16x8_t u16 = vcombine_u16(u16_lo, u16_hi);
    vst1q_u16(reinterpret_cast<uint16_t*>(dst), u16);
  }
}

static AVS_FORCEINLINE float32x4_t overlay_blend_neon_core(const float32x4_t& p1_f, const float32x4_t& p2_f, const float32x4_t& factor) {
  // p1 + (p2 - p1) * factor
  return vmlaq_f32(p1_f, vsubq_f32(p2_f, p1_f), factor);
}

template<bool has_mask, typename pixel_t>
void overlay_blend_neon_uint_slow(BYTE* p1, const BYTE* p2, const BYTE* mask,
  const int p1_pitch, const int p2_pitch, const int mask_pitch,
  const int width, const int height, const int opacity, const float opacity_f, const int bits_per_pixel)
{
  const float32x4_t rounder = vdupq_n_f32(0.5f);
  const int max_pixel_value = (1 << bits_per_pixel) - 1;
  const float factor = has_mask ? (opacity_f / max_pixel_value) : opacity_f;
  const float32x4_t factor_v = vdupq_n_f32(factor);

  const int realwidth = width * sizeof(pixel_t);
  constexpr int bytes_per_cycle = 8 * sizeof(pixel_t);
  int wMod8 = (realwidth / bytes_per_cycle) * bytes_per_cycle;

  for (int y = 0; y < height; y++) {
    for (int x = 0; x < wMod8; x += bytes_per_cycle) {
      float32x4_t unpacked_p1_lo, unpacked_p1_hi;
      float32x4_t unpacked_p2_lo, unpacked_p2_hi;
      Eightpixels_to_Eightfloats_neon<pixel_t>((const pixel_t*)(p1 + x), unpacked_p1_lo, unpacked_p1_hi);
      Eightpixels_to_Eightfloats_neon<pixel_t>((const pixel_t*)(p2 + x), unpacked_p2_lo, unpacked_p2_hi);

      float32x4_t result_lo, result_hi;
      if constexpr (has_mask) {
        float32x4_t unpacked_mask_lo, unpacked_mask_hi;
        Eightpixels_to_Eightfloats_neon<pixel_t>((const pixel_t*)(mask + x), unpacked_mask_lo, unpacked_mask_hi);
        unpacked_mask_lo = vmulq_f32(unpacked_mask_lo, factor_v);
        unpacked_mask_hi = vmulq_f32(unpacked_mask_hi, factor_v);
        result_lo = overlay_blend_neon_core(unpacked_p1_lo, unpacked_p2_lo, unpacked_mask_lo);
        result_hi = overlay_blend_neon_core(unpacked_p1_hi, unpacked_p2_hi, unpacked_mask_hi);
      } else {
        result_lo = overlay_blend_neon_core(unpacked_p1_lo, unpacked_p2_lo, factor_v);
        result_hi = overlay_blend_neon_core(unpacked_p1_hi, unpacked_p2_hi, factor_v);
      }

      Store_Eightpixels_neon<pixel_t>((pixel_t*)(p1 + x), result_lo, result_hi, rounder);
    }

    // Leftover value
    for (int x = wMod8 / sizeof(pixel_t); x < width; x++) {
      const float new_factor = has_mask ? static_cast<float>(reinterpret_cast<const pixel_t*>(mask)[x]) * factor : factor;
      auto result = overlay_blend_c_core_simple(reinterpret_cast<pixel_t*>(p1)[x], reinterpret_cast<const pixel_t*>(p2)[x], new_factor);
      reinterpret_cast<pixel_t*>(p1)[x] = (pixel_t)(result + 0.5f);
    }

    p1 += p1_pitch;
    p2 += p2_pitch;
    if constexpr(has_mask)
      mask += mask_pitch;
  }
}

// Process_SixteenPixels_neon:
// To the greatest surprise, gcc 14 optimized the original (non-mask) C code so perfectly, that it took
// two days of refining to achieve their speed using these hand-crafted intrinsics :)
// Still, I leave this here, since it can work quickly in a debug build as well, and quicker than C from LLVM 19.

template<bool has_mask, typename pixel_t>
static AVS_FORCEINLINE void Process_SixteenPixels_neon(
  BYTE* p1_ptr, const BYTE* p2_ptr, const BYTE* mask_ptr, const float32x4_t factor_v, const float32x4_t rounder, const int offset)
{
  // The core blending function (p1 + (p2 - p1) * factor)
  auto overlay_blend_neon_core = [](const float32x4_t& p1_f, const float32x4_t& delta_f, const float32x4_t& factor) -> float32x4_t {
    // 1. Calculate the blend: p1_f + (delta_f * factor) (where delta_f = p2_f - p1_f)
    return vfmaq_f32(p1_f, delta_f, factor); // A + B * C
    };

  // --- CASE 1: 8-bit Pixels (uint8_t) ---
  if constexpr (sizeof(pixel_t) == 1) {
    // Pointers
    const uint8_t* src1 = reinterpret_cast<const uint8_t*>(p1_ptr + offset);
    const uint8_t* src2 = reinterpret_cast<const uint8_t*>(p2_ptr + offset);

    // 1. Load 16 bytes (128-bit)
    uint8x16_t q_p1 = vld1q_u8(src1);
    uint8x16_t q_p2 = vld1q_u8(src2);

    // 2. Calc Delta and Widen P1 (Keep everything in 128-bit Q regs)
    // Use _high intrinsics to generate 'usubl2' and 'uxtl2' instructions directly

    // Low 8 pixels
    // blending is p1 + (p2-p1)*f. We want p2-p1.
    int16x8_t delta_16_lo = vreinterpretq_s16_u16(vsubl_u8(vget_low_u8(q_p2), vget_low_u8(q_p1)));

    // High 8 pixels (Direct Hardware Instruction: usubl2)
    int16x8_t delta_16_hi = vreinterpretq_s16_u16(vsubl_high_u8(q_p2, q_p1));

    // Widen P1, cast to signed
    int16x8_t p1_16_lo = vreinterpretq_s16_u16(vmovl_u8(vget_low_u8(q_p1)));
    int16x8_t p1_16_hi = vreinterpretq_s16_u16(vmovl_high_u8(q_p1));

    float32x4_t res0, res1, res2, res3;

    if constexpr(has_mask) {
      // const float new_factor = has_mask ? static_cast<float>(reinterpret_cast<const pixel_t*>(mask)[x]) * factor : factor;
      // trouble: factor is fixed global here, mask is individual per pixel
      // Since mask byte is in [0,255], the original factor [0..1] was already scaled down by 255 (8 bit max_pixel_value)
      // So we don't have to take care here of max_pixel_value again.
      // Is it optimal to calculate the multiplied factor_v independently?
      // each 4 pixels must converted to 32 bit (uint8 -> signed 16 -> signed 32 -> float, like p1)

      // Next calculation needs to be changed accordingly, take care of minimal lo-hi loads,
      // stay in 128-bit pipeline as long as we can, multiply factor_v's with actual mask values

      const uint8_t* src_mask = reinterpret_cast<const uint8_t*>(mask_ptr + offset);

      // 1. Load 16 mask bytes
      uint8x16_t q_mask = vld1q_u8(src_mask);

      // 2. Widen 8-bit mask to 16-bit (low/high 8 pixels)
      // and cast to signed to spare instructions on 32 bit widening
      int16x8_t mask_16_lo = vreinterpretq_s16_u16(vmovl_u8(vget_low_u8(q_mask)));
      int16x8_t mask_16_hi = vreinterpretq_s16_u16(vmovl_high_u8(q_mask));

      // 3. Widen 16-bit mask to 32-bit (4 blocks of 4 pixels) and convert to float
      // Stay in Signed domain, uint8_t fits into int16_t without issues
      // Block 0 (Pixels 0-3) - new_factor0 = (float)mask[0-3] * factor_v
      float32x4_t new_factor0 = vcvtq_f32_s32(vmovl_s16(vget_low_s16(mask_16_lo)));

      // Block 1 (Pixels 4-7)
      float32x4_t new_factor1 = vcvtq_f32_s32(vmovl_high_s16(mask_16_lo));

      // Block 2 (Pixels 8-11)
      float32x4_t new_factor2 = vcvtq_f32_s32(vmovl_s16(vget_low_s16(mask_16_hi)));

      // Block 3 (Pixels 12-15)
      float32x4_t new_factor3 = vcvtq_f32_s32(vmovl_high_s16(mask_16_hi));

      // 4. Multiply mask float vectors by the constant global factor_v
      new_factor0 = vmulq_f32(new_factor0, factor_v);
      new_factor1 = vmulq_f32(new_factor1, factor_v);
      new_factor2 = vmulq_f32(new_factor2, factor_v);
      new_factor3 = vmulq_f32(new_factor3, factor_v);

      // 5. Blend using the new per-pixel factors
      // Block 0 (Pixels 0-3)
      res0 = overlay_blend_neon_core(
        vcvtq_f32_s32(vmovl_s16(vget_low_s16(p1_16_lo))),
        vcvtq_f32_s32(vmovl_s16(vget_low_s16(delta_16_lo))),
        new_factor0); // *** Use new_factor0 ***

      // Block 1 (Pixels 4-7)
      res1 = overlay_blend_neon_core(
        vcvtq_f32_s32(vmovl_high_s16(p1_16_lo)),
        vcvtq_f32_s32(vmovl_high_s16(delta_16_lo)),
        new_factor1); // *** Use new_factor1 ***

      // Block 2 (Pixels 8-11)
      res2 = overlay_blend_neon_core(
        vcvtq_f32_s32(vmovl_s16(vget_low_s16(p1_16_hi))),
        vcvtq_f32_s32(vmovl_s16(vget_low_s16(delta_16_hi))),
        new_factor2); // *** Use new_factor2 ***

      // Block 3 (Pixels 12-15)
      res3 = overlay_blend_neon_core(
        vcvtq_f32_s32(vmovl_high_s16(p1_16_hi)),
        vcvtq_f32_s32(vmovl_high_s16(delta_16_hi)),
        new_factor3); // *** Use new_factor3 ***

    }
    else {
      // Keep highly optimized non-mask path

      // 3. Convert to Float and Blend (4 blocks of 4 pixels)
      // We treat the 16-bit vectors as 32-bit vectors during conversion to save moves

      // Block 0 (Pixels 0-3)
      res0 = overlay_blend_neon_core(
        vcvtq_f32_s32(vmovl_s16(vget_low_s16(p1_16_lo))),
        vcvtq_f32_s32(vmovl_s16(vget_low_s16(delta_16_lo))),
        factor_v);

      // Block 1 (Pixels 4-7) - use _high to avoid 'dup'
      res1 = overlay_blend_neon_core(
        vcvtq_f32_s32(vmovl_high_s16(p1_16_lo)),
        vcvtq_f32_s32(vmovl_high_s16(delta_16_lo)),
        factor_v);

      // Block 2 (Pixels 8-11)
      res2 = overlay_blend_neon_core(
        vcvtq_f32_s32(vmovl_s16(vget_low_s16(p1_16_hi))),
        vcvtq_f32_s32(vmovl_s16(vget_low_s16(delta_16_hi))),
        factor_v);

      // Block 3 (Pixels 12-15) - use _high
      res3 = overlay_blend_neon_core(
        vcvtq_f32_s32(vmovl_high_s16(p1_16_hi)),
        vcvtq_f32_s32(vmovl_high_s16(delta_16_hi)),
        factor_v);
    }

    // 4. Rounding
    res0 = vaddq_f32(res0, rounder);
    res1 = vaddq_f32(res1, rounder);
    res2 = vaddq_f32(res2, rounder);
    res3 = vaddq_f32(res3, rounder);

    // 5. Convert to Int32
    int32x4_t s32_0 = vcvtq_s32_f32(res0);
    int32x4_t s32_1 = vcvtq_s32_f32(res1);
    int32x4_t s32_2 = vcvtq_s32_f32(res2);
    int32x4_t s32_3 = vcvtq_s32_f32(res3);

    // 6. FAST NARROWING (The "Compiler Trick")
    // Instead of vmovn (which goes 128->64 bit), we use vuzp1q (Unzip).
    // This packs the low 16-bits of every 32-bit element from two vectors into one 128-bit vector.
    // It stays in the 128-bit pipeline.

    // Narrow 32-bit -> 16-bit
    // Takes Pixels 0-3 (s32_0) and Pixels 4-7 (s32_1) -> Pixels 0-7 (u16_lo)
    int16x8_t u16_lo = vuzp1q_s16(vreinterpretq_s16_s32(s32_0), vreinterpretq_s16_s32(s32_1));

    // Takes Pixels 8-11 (s32_2) and Pixels 12-15 (s32_3) -> Pixels 8-15 (u16_hi)
    int16x8_t u16_hi = vuzp1q_s16(vreinterpretq_s16_s32(s32_2), vreinterpretq_s16_s32(s32_3));

    // Narrow 16-bit -> 8-bit
    // Takes Pixels 0-7 (u16_lo) and Pixels 8-15 (u16_hi) -> Pixels 0-15 (result)
    // Note: We cast to s8 to perform the unzip on bytes
    uint8x16_t result = vreinterpretq_u8_s8(
      vuzp1q_s8(vreinterpretq_s8_s16(u16_lo), vreinterpretq_s8_s16(u16_hi))
    );

    // 7. Store
    vst1q_u8(reinterpret_cast<uint8_t*>(p1_ptr + offset), result);
  }
  // --- CASE 2: 16-bit Pixels (uint16_t) ---
  else if constexpr (sizeof(pixel_t) == 2) {
    const uint16_t* src1 = reinterpret_cast<const uint16_t*>(p1_ptr + offset);
    const uint16_t* src2 = reinterpret_cast<const uint16_t*>(p2_ptr + offset);

    // 1. Bulk load + delta (like 8-bit success)
    uint16x8_t v_p1 = vld1q_u16(src1);
    uint16x8_t v_p2 = vld1q_u16(src2);

    // 2. Bulk signed deltas (p2-p1) - vsubl_u16 → sxtl
    int32x4_t delta_lo = vreinterpretq_s32_u32(vsubl_u16(vget_low_u16(v_p2), vget_low_u16(v_p1)));
    int32x4_t delta_hi = vreinterpretq_s32_u32(vsubl_high_u16(v_p2, v_p1));  // Direct USUBL2

    // 3. because we can have full 16 bit, we use unsigned 16 bit before 32 bit widening
    int32x4_t p1_s_lo = vreinterpretq_s32_u32(vmovl_u16(vget_low_u16(v_p1)));
    int32x4_t p1_s_hi = vreinterpretq_s32_u32(vmovl_high_u16(v_p1));

    // 4. Single scvtf path (matches 8-bit victory)
    float32x4_t p1_f_lo = vcvtq_f32_s32(p1_s_lo);
    float32x4_t delta_f_lo = vcvtq_f32_s32(delta_lo);
    float32x4_t p1_f_hi = vcvtq_f32_s32(p1_s_hi);
    float32x4_t delta_f_hi = vcvtq_f32_s32(delta_hi);

    float32x4_t res_lo, res_hi;

    if constexpr (has_mask) {
      const uint16_t* src_mask = reinterpret_cast<const uint16_t*>(mask_ptr + offset);
      // const float new_factor = has_mask ? static_cast<float>(reinterpret_cast<const pixel_t*>(mask)[x]) * factor : factor;
      // trouble: factor is fixed global here, mask is individual per pixel
      // Since mask byte is in [0,(max_pixel_value-1)], the original factor [0..1] was already scaled down by (2^N - 1) 1023,4095,16383,65535 (10-16 bit max_pixel_value)
      // So we don't have to take care here of max_pixel_value again.
      // Is it optimal to calculate the multiplied factor_v independently?
      // each 4 pixels must converted to 32 bit (uint16 -> signed 32 -> float, like p1)

      // Next calculation needs to be changed accordingly, take care of minimal lo-hi loads,
      // stay in 128-bit pipeline as long as we can, multiply factor_v's with actual mask values

      // 1. Load 16 mask bytes
      uint16x8_t v_mask = vld1q_u16(src_mask);

      // 2. because we can have full 16 bit, we use unsigned 16 bit before 32 bit widening
      int32x4_t mask_16_lo = vreinterpretq_s32_u32(vmovl_u16(vget_low_u16(v_mask))); // LSL, SXT
      int32x4_t mask_16_hi = vreinterpretq_s32_u32(vmovl_high_u16(v_mask));

      // 3. Widen 16-bit mask to 32-bit (4 blocks of 4 pixels) and convert to float
      // Unlike at 8 bit, we cannot stay in Signed domain

      // 4. Convert 32-bit mask to float
      float32x4_t mask_f_lo = vcvtq_f32_s32(mask_16_lo);
      float32x4_t mask_f_hi = vcvtq_f32_s32(mask_16_hi);

      // 4. Multiply mask float vectors by the constant global factor_v
      float32x4_t new_factor_lo = vmulq_f32(mask_f_lo, factor_v);
      float32x4_t new_factor_hi = vmulq_f32(mask_f_hi, factor_v);

      // 5. FMA + round
      // 5. FMA using the new per-pixel factors
      res_lo = vfmaq_f32(p1_f_lo, delta_f_lo, new_factor_lo); // *** Use new_factor_lo ***
      res_hi = vfmaq_f32(p1_f_hi, delta_f_hi, new_factor_hi); // *** Use new_factor_hi ***

    }
    else {
      // Keep highly optimized non-mask path

    // 5. FMA + round
      res_lo = vfmaq_f32(p1_f_lo, delta_f_lo, factor_v);
      res_hi = vfmaq_f32(p1_f_hi, delta_f_hi, factor_v);
    }

    // Rounding +0.5f
    res_lo = vaddq_f32(res_lo, rounder);
    res_hi = vaddq_f32(res_hi, rounder);

    // 6. FAST NON-SAT NARROW (like 8-bit uzp1 magic)
    // With factor in [0,1], and p1,p2 in [0,65535], rounder 0.5f, the truncated result is also in [0,65535].
    int32x4_t i32_lo = vcvtq_s32_f32(res_lo);
    int32x4_t i32_hi = vcvtq_s32_f32(res_hi);

    // Direct narrow (no vqmovun overhead)
    int16x4_t h16_lo = vmovn_s32(i32_lo);      // MOVN, 1 cycle
    int16x4_t h16_hi = vmovn_s32(i32_hi);
    uint16x8_t result = vreinterpretq_u16_s16(vcombine_s16(h16_lo, h16_hi));

    vst1q_u16(reinterpret_cast<uint16_t*>(p1_ptr + offset), result);

  }
}

/*
* Final new result asm:
* _Z23overlay_blend_neon_uintILb0EhEvPhPKhS2_iiiiiifi
  ldr	q7, [x0, x19]	//, MEM <__Uint8x16_t> [(unsigned char * {ref-all})p1_188 + ivtmp.288_938 * 1]
  ldr	q19, [x1, x19]	//, MEM <__Uint8x16_t> [(unsigned char * {ref-all})p2_928 + ivtmp.288_938 * 1]
  zip1	v16.16b, v7.16b, v31.16b	//,,
  usubl	v1.8h, v19.8b, v7.8b	//,,
  usubl2	v19.8h, v19.16b, v7.16b	//,,
  zip2	v7.16b, v7.16b, v31.16b	//,,
  sxtl	v18.4s, v16.4h	//,
  sxtl	v23.4s, v1.4h	//,
  sxtl	v20.4s, v19.4h	//,
  sxtl	v17.4s, v7.4h	//,
  sxtl2	v16.4s, v16.8h	//,
  sxtl2	v1.4s, v1.8h	//,
  sxtl2	v7.4s, v7.8h	//,
  sxtl2	v19.4s, v19.8h	//,
  scvtf	v18.4s, v18.4s	//,
  scvtf	v16.4s, v16.4s	//,
  scvtf	v23.4s, v23.4s	//,
  scvtf	v1.4s, v1.4s	//,
  scvtf	v17.4s, v17.4s	//,
  scvtf	v7.4s, v7.4s	//,
  scvtf	v20.4s, v20.4s	//,
  scvtf	v19.4s, v19.4s	//,
  fmla	v18.4s, v23.4s, v30.4s	//,,
  fmla	v16.4s, v1.4s, v30.4s	//,,
  fmla	v17.4s, v20.4s, v30.4s	//,,
  fmla	v7.4s, v19.4s, v30.4s	//,,
  fadd	v18.4s, v18.4s, v29.4s	//,,
  fadd	v16.4s, v16.4s, v29.4s	//,,
  fadd	v17.4s, v17.4s, v29.4s	//,,
  fadd	v7.4s, v7.4s, v29.4s	//,,
  fcvtzs	v18.4s, v18.4s	//,
  fcvtzs	v16.4s, v16.4s	//,
  fcvtzs	v17.4s, v17.4s	//,
  fcvtzs	v7.4s, v7.4s	//,
  uzp1	v16.8h, v18.8h, v16.8h	//,,
  uzp1	v7.8h, v17.8h, v7.8h	//,,
  uzp1	v7.16b, v16.16b, v7.16b	//,,
  str	q7, [x0, x19]	//, MEM <__Uint8x16_t> [(unsigned char * {ref-all})p1_188 + ivtmp.288_938 * 1]
pre SIMD:
  ldr	q3, [x0, x19]	//, MEM <__Uint8x16_t> [(unsigned char * {ref-all})p1_194 + ivtmp.289_958 * 1]
  ldr	q26, [x1, x19]	//, MEM <__Uint8x16_t> [(unsigned char * {ref-all})p2_948 + ivtmp.289_958 * 1]
  dup	d27, v3.d[1]	//,
  zip1	v18.16b, v3.16b, v31.16b	//,,
  dup	d1, v26.d[1]	//,
  usubl	v26.8h, v26.8b, v3.8b	//,,
  zip1	v17.16b, v27.16b, v31.16b	//,,
  zip1	v23.8h, v18.8h, v31.8h	//,,
  usubl	v1.8h, v1.8b, v27.8b	//,,
  sxtl	v22.4s, v26.4h	//,
  zip1	v19.8h, v17.8h, v31.8h	//,,
  zip2	v18.8h, v18.8h, v31.8h	//,,
  sxtl	v2.4s, v1.4h	//,
  zip2	v17.8h, v17.8h, v31.8h	//,,
  sxtl2	v26.4s, v26.8h	//,
  sxtl2	v1.4s, v1.8h	//,
  ucvtf	v19.4s, v19.4s	//,
  ucvtf	v17.4s, v17.4s	//,
  scvtf	v2.4s, v2.4s	//,
  scvtf	v1.4s, v1.4s	//,
  ucvtf	v23.4s, v23.4s	//,
  ucvtf	v18.4s, v18.4s	//,
  scvtf	v22.4s, v22.4s	//,
  scvtf	v26.4s, v26.4s	//,
  fmla	v19.4s, v2.4s, v30.4s	//,,
  fmla	v17.4s, v1.4s, v30.4s	//,,
  fmla	v23.4s, v22.4s, v30.4s	//,,
  fmla	v18.4s, v26.4s, v30.4s	//,,
  fadd	v19.4s, v19.4s, v29.4s	//,,
  fadd	v17.4s, v17.4s, v29.4s	//,,
  fadd	v23.4s, v23.4s, v29.4s	//,,
  fadd	v18.4s, v18.4s, v29.4s	//,,
  fcvtzs	v19.4s, v19.4s	//,
  fcvtzs	v17.4s, v17.4s	//,
  fcvtzs	v23.4s, v23.4s	//,
  fcvtzs	v18.4s, v18.4s	//,
  uzp1	v17.8h, v19.8h, v17.8h	//,,
  uzp1	v18.8h, v23.8h, v18.8h	//,,
  uzp1	v17.16b, v18.16b, v17.16b	//,,
  str	q17, [x0, x19]	//, MEM <__Uint8x16_t> [(unsigned char * {ref-all})p1_194 + ivtmp.289_958 * 1]

PrevSIMD:
  ldr	q27, [x0, x19]	//, MEM <__Uint8x16_t> [(unsigned char * {ref-all})p1_182 + ivtmp.289_958 * 1]
  ldr	q1, [x1, x19]	//, MEM <__Uint8x16_t> [(unsigned char * {ref-all})p2_948 + ivtmp.289_958 * 1]
  zip1	v18.16b, v27.16b, v31.16b	//,,
  usubl	v24.8h, v1.8b, v27.8b	//,,
  dup	d27, v27.d[1]	//,
  dup	d1, v1.d[1]	//,
  zip1	v25.8h, v18.8h, v31.8h	//,,
  sxtl	v23.4s, v24.4h	//,
  zip1	v17.16b, v27.16b, v31.16b	//,,
  usubl	v1.8h, v1.8b, v27.8b	//,,
  sxtl2	v24.4s, v24.8h	//,
  zip2	v18.8h, v18.8h, v31.8h	//,,
  zip1	v19.8h, v17.8h, v31.8h	//,,
  sxtl	v2.4s, v1.4h	//,
  zip2	v17.8h, v17.8h, v31.8h	//,,
  sxtl2	v1.4s, v1.8h	//,
  ucvtf	v25.4s, v25.4s	//,
  scvtf	v2.4s, v2.4s	//,
  ucvtf	v19.4s, v19.4s	//,
  scvtf	v1.4s, v1.4s	//,
  ucvtf	v17.4s, v17.4s	//,
  ucvtf	v18.4s, v18.4s	//,
  scvtf	v23.4s, v23.4s	//,
  scvtf	v24.4s, v24.4s	//,
  fmla	v19.4s, v2.4s, v30.4s	//,,
  fmla	v17.4s, v1.4s, v30.4s	//,,
  fmla	v25.4s, v23.4s, v30.4s	//,,
  fmla	v18.4s, v24.4s, v30.4s	//,,
  fadd	v19.4s, v19.4s, v29.4s	//,,
  fadd	v17.4s, v17.4s, v29.4s	//,,
  fadd	v25.4s, v25.4s, v29.4s	//,,
  fadd	v18.4s, v18.4s, v29.4s	//,,
  fcvtzs	v19.4s, v19.4s	//,
  fcvtzs	v17.4s, v17.4s	//,
  fcvtzs	v25.4s, v25.4s	//,
  fcvtzs	v18.4s, v18.4s	//,
  uzp1	v17.8h, v19.8h, v17.8h	//,,
  uzp1	v18.8h, v25.8h, v18.8h	//,,
  uzp1	v17.16b, v18.16b, v17.16b	//,,
  str	q17, [x0, x19]	//, MEM <__Uint8x16_t> [(unsigned char * {ref-all})p1_182 + ivtmp.289_958 * 1]

 And the quicker C version:

    ldr    q4, [x0, x16]    //, MEM <vector(16) unsigned char> [(BYTE *)p1_218 + ivtmp.60_46 * 1]
    ldr    q25, [x1, x16]    //, MEM <const vector(16) unsigned char> [(const BYTE *)p2_62 + ivtmp.60_46 * 1]
    zip1    v2.16b, v4.16b, v31.16b    //,,
    zip2    v26.16b, v4.16b, v31.16b    //,,
    usubl    v22.8h, v25.8b, v4.8b    //,,
    usubl2    v25.8h, v25.16b, v4.16b    //,,
    zip1    v3.8h, v2.8h, v31.8h    //,,
    zip1    v1.8h, v26.8h, v31.8h    //,,
    sxtl2    v21.4s, v22.8h    //,
    sxtl    v28.4s, v25.4h    //,
    zip2    v2.8h, v2.8h, v31.8h    //,,
    sxtl    v22.4s, v22.4h    //,
    zip2    v26.8h, v26.8h, v31.8h    //,,
    sxtl2    v25.4s, v25.8h    //,
    scvtf    v21.4s, v21.4s    //,
    scvtf    v22.4s, v22.4s    //,
    scvtf    v3.4s, v3.4s    //,
    scvtf    v2.4s, v2.4s    //,
    scvtf    v1.4s, v1.4s    //,
    scvtf    v26.4s, v26.4s    //,
    scvtf    v28.4s, v28.4s    //,
    scvtf    v25.4s, v25.4s    //,
    fmla    v3.4s, v22.4s, v30.4s    //,,
    fmla    v2.4s, v21.4s, v30.4s    //,,
    fmla    v1.4s, v28.4s, v30.4s    //,,
    fmla    v26.4s, v25.4s, v30.4s    //,,
    fadd    v3.4s, v3.4s, v29.4s    //,,
    fadd    v2.4s, v2.4s, v29.4s    //,,
    fadd    v1.4s, v1.4s, v29.4s    //,,
    fadd    v26.4s, v26.4s, v29.4s    //,,
    fcvtzs    v3.4s, v3.4s    //,
    fcvtzs    v2.4s, v2.4s    //,
    fcvtzs    v1.4s, v1.4s    //,
    fcvtzs    v26.4s, v26.4s    //,
    uzp1    v2.8h, v3.8h, v2.8h    //,,
    uzp1    v26.8h, v1.8h, v26.8h    //,,
    uzp1    v26.16b, v2.16b, v26.16b    //,,
    str    q26, [x0, x16]    //, MEM <vector(16) unsigned char> [(BYTE *)p1_218 + ivtmp.60_46 * 1]


*/


template<bool has_mask, typename pixel_t>
void overlay_blend_neon_uint(BYTE* p1, const BYTE* p2, const BYTE* mask,
  const int p1_pitch, const int p2_pitch, const int mask_pitch,
  const int width, const int height, const int opacity, const float opacity_f, const int bits_per_pixel)
{
  const float32x4_t rounder = vdupq_n_f32(0.5f);
  const int max_pixel_value = (1 << bits_per_pixel) - 1;
  const float factor = has_mask ? (opacity_f / max_pixel_value) : opacity_f;
  const float32x4_t factor_v = vdupq_n_f32(factor);

  // Calculate the vectorized limit (width rounded down to the nearest multiple of 16 bytes)
  constexpr int bytes_per_pixel = sizeof(pixel_t);
  constexpr int pixels_per_cycle = (bytes_per_pixel == 1) ? 16 : 8; // 16 bytes per cycle
  const int real_width_bytes = width * bytes_per_pixel;
  const int vector_limit_bytes = (real_width_bytes / (pixels_per_cycle * bytes_per_pixel)) * (pixels_per_cycle * bytes_per_pixel);

  // The vector processing size is always 16 bytes for either 8-bit (16 pixels) or 16-bit (8 pixels).
  constexpr int bytes_per_cycle = 16;

  for (int y = 0; y < height; y++) {
    // --- Vector Loop (16 bytes at a time) ---
    for (int x_bytes = 0; x_bytes < vector_limit_bytes; x_bytes += bytes_per_cycle) {
      // Call the core 16-byte processing function
      Process_SixteenPixels_neon<has_mask, pixel_t>(p1, p2, mask, factor_v, rounder, x_bytes);
    }

    // --- Scalar Tail (Pixel by pixel) ---
    const int start_pixel = vector_limit_bytes / bytes_per_pixel;
    for (int x = start_pixel; x < width; x++) {
      const float new_factor = has_mask ? static_cast<float>(reinterpret_cast<const pixel_t*>(mask)[x]) * factor : factor;
      auto result = overlay_blend_c_core_simple(reinterpret_cast<pixel_t*>(p1)[x], reinterpret_cast<const pixel_t*>(p2)[x], new_factor);
      reinterpret_cast<pixel_t*>(p1)[x] = (pixel_t)(result + 0.5f);
    }

    p1 += p1_pitch;
    p2 += p2_pitch;
    if constexpr (has_mask)
      mask += mask_pitch;
  }
}

// instantiate
template void overlay_blend_neon_uint<true, uint8_t>(BYTE* p1, const BYTE* p2, const BYTE* mask,
  const int p1_pitch, const int p2_pitch, const int mask_pitch, const int width, const int height, const int opacity, const float opacity_f, const int bits_per_pixel);
template void overlay_blend_neon_uint<true, uint16_t>(BYTE* p1, const BYTE* p2, const BYTE* mask,
  const int p1_pitch, const int p2_pitch, const int mask_pitch, const int width, const int height, const int opacity, const float opacity_f, const int bits_per_pixel);
template void overlay_blend_neon_uint<false, uint8_t>(BYTE* p1, const BYTE* p2, const BYTE* mask,
  const int p1_pitch, const int p2_pitch, const int mask_pitch, const int width, const int height, const int opacity, const float opacity_f, const int bits_per_pixel);
template void overlay_blend_neon_uint<false, uint16_t>(BYTE* p1, const BYTE* p2, const BYTE* mask,
  const int p1_pitch, const int p2_pitch, const int mask_pitch, const int width, const int height, const int opacity, const float opacity_f, const int bits_per_pixel);

// NEON equivalent of overlay_blend_opaque_sse2_core
AVS_FORCEINLINE uint8x16_t overlay_blend_opaque_neon_core(const uint8x16_t& p1, const uint8x16_t& p2, const uint8x16_t& mask) {
  // return (mask) ? p2 : p1;
  // mask: 0x00 = p1, 0xFF = p2 (same as SSE)
  // vbslq_u8(mask, p2, p1): for each bit, if mask=1, take from p2, else from p1
  return vbslq_u8(mask, p2, p1);
}

// For 16-bit pixels
AVS_FORCEINLINE uint16x8_t overlay_blend_opaque_neon_core_u16(const uint16x8_t& p1, const uint16x8_t& p2, const uint16x8_t& mask) {
  return vbslq_u16(mask, p2, p1);
}

// Compare functions for lighten and darken mode (8-bit)
AVS_FORCEINLINE static uint8x16_t overlay_darken_neon_cmp(const uint8x16_t& p1, const uint8x16_t& p2) {
  // mask = (p2 <= p1) ? 0xFF : 0x00
  return vorrq_u8(vcltq_u8(p2, p1), vceqq_u8(p2, p1)); // (p2 < p1) | (p2 == p1)
}

AVS_FORCEINLINE static uint8x16_t overlay_lighten_neon_cmp(const uint8x16_t& p1, const uint8x16_t& p2) {
  // mask = (p2 >= p1) ? 0xFF : 0x00
  return vorrq_u8(vcgtq_u8(p2, p1), vceqq_u8(p2, p1)); // (p2 > p1) | (p2 == p1)
}

// Compare functions for lighten and darken mode
AVS_FORCEINLINE static int overlay_darken_c_cmp(BYTE p1, BYTE p2) {
  return p2 <= p1;
}

AVS_FORCEINLINE static int overlay_lighten_c_cmp(BYTE p1, BYTE p2) {
  return p2 >= p1;
}

/***************************************
 ********* Mode: Lighten/Darken ********
 ***************************************/

using OverlayNeonCompare = uint8x16_t(const uint8x16_t& p1, const uint8x16_t& p2);
using OverlayCCompare = int(BYTE, BYTE);

// Main NEON version for 8-bit
template <OverlayNeonCompare compare, OverlayCCompare compare_c>
void overlay_darklighten_neon(BYTE *p1Y, BYTE *p1U, BYTE *p1V, const BYTE *p2Y, const BYTE *p2U, const BYTE *p2V,
                              int p1_pitch, int p2_pitch, int width, int height) {
  int wMod16 = (width / 16) * 16;
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < wMod16; x += 16) {
      uint8x16_t p1_y = vld1q_u8(p1Y + x);
      uint8x16_t p2_y = vld1q_u8(p2Y + x);
      uint8x16_t mask = compare(p1_y, p2_y);
      uint8x16_t result_y = overlay_blend_opaque_neon_core(p1_y, p2_y, mask);
      vst1q_u8(p1Y + x, result_y);

      uint8x16_t p1_u = vld1q_u8(p1U + x);
      uint8x16_t p2_u = vld1q_u8(p2U + x);
      uint8x16_t result_u = overlay_blend_opaque_neon_core(p1_u, p2_u, mask);
      vst1q_u8(p1U + x, result_u);

      uint8x16_t p1_v = vld1q_u8(p1V + x);
      uint8x16_t p2_v = vld1q_u8(p2V + x);
      uint8x16_t result_v = overlay_blend_opaque_neon_core(p1_v, p2_v, mask);
      vst1q_u8(p1V + x, result_v);
    }
    // Leftover value (fallback to C)
    for (int x = wMod16; x < width; x++) {
      int mask = compare_c(p1Y[x], p2Y[x]) ? 0xFF : 0x00;
      p1Y[x] = overlay_blend_opaque_c_core<uint8_t>(p1Y[x], p2Y[x], mask);
      p1U[x] = overlay_blend_opaque_c_core<uint8_t>(p1U[x], p2U[x], mask);
      p1V[x] = overlay_blend_opaque_c_core<uint8_t>(p1V[x], p2V[x], mask);
    }
    p1Y += p1_pitch; p1U += p1_pitch; p1V += p1_pitch;
    p2Y += p2_pitch; p2U += p2_pitch; p2V += p2_pitch;
  }
}


template<bool has_mask>
void overlay_blend_neon_float(BYTE* p1, const BYTE* p2, const BYTE* mask,
  const int p1_pitch, const int p2_pitch, const int mask_pitch,
  const int width, const int height, const int /*opacity*/, const float opacity_f, const int bits_per_pixel)
{
  const int realwidth = width * sizeof(float);
  int wMod16 = (realwidth / 16) * 16;
  float32x4_t opacity_v = vdupq_n_f32(opacity_f);

  for (int y = 0; y < height; y++) {
    for (int x = 0; x < wMod16; x += 16) {
      float32x4_t p1_f0 = vld1q_f32(reinterpret_cast<const float*>(p1 + x + 0));
      float32x4_t p1_f1 = vld1q_f32(reinterpret_cast<const float*>(p1 + x + 4));
      float32x4_t p1_f2 = vld1q_f32(reinterpret_cast<const float*>(p1 + x + 8));
      float32x4_t p1_f3 = vld1q_f32(reinterpret_cast<const float*>(p1 + x + 12));

      float32x4_t p2_f0 = vld1q_f32(reinterpret_cast<const float*>(p2 + x + 0));
      float32x4_t p2_f1 = vld1q_f32(reinterpret_cast<const float*>(p2 + x + 4));
      float32x4_t p2_f2 = vld1q_f32(reinterpret_cast<const float*>(p2 + x + 8));
      float32x4_t p2_f3 = vld1q_f32(reinterpret_cast<const float*>(p2 + x + 12));

      float32x4_t mask0, mask1, mask2, mask3;
      if constexpr (has_mask) {
        mask0 = vld1q_f32(reinterpret_cast<const float*>(mask + x + 0));
        mask1 = vld1q_f32(reinterpret_cast<const float*>(mask + x + 4));
        mask2 = vld1q_f32(reinterpret_cast<const float*>(mask + x + 8));
        mask3 = vld1q_f32(reinterpret_cast<const float*>(mask + x + 12));
        mask0 = vmulq_f32(mask0, opacity_v);
        mask1 = vmulq_f32(mask1, opacity_v);
        mask2 = vmulq_f32(mask2, opacity_v);
        mask3 = vmulq_f32(mask3, opacity_v);
      }
      else {
        mask0 = opacity_v;
        mask1 = opacity_v;
        mask2 = opacity_v;
        mask3 = opacity_v;
      }

      float32x4_t result0 = vaddq_f32(p1_f0, vmulq_f32(vsubq_f32(p2_f0, p1_f0), mask0));
      float32x4_t result1 = vaddq_f32(p1_f1, vmulq_f32(vsubq_f32(p2_f1, p1_f1), mask1));
      float32x4_t result2 = vaddq_f32(p1_f2, vmulq_f32(vsubq_f32(p2_f2, p1_f2), mask2));
      float32x4_t result3 = vaddq_f32(p1_f3, vmulq_f32(vsubq_f32(p2_f3, p1_f3), mask3));

      vst1q_f32(reinterpret_cast<float*>(p1 + x + 0), result0);
      vst1q_f32(reinterpret_cast<float*>(p1 + x + 4), result1);
      vst1q_f32(reinterpret_cast<float*>(p1 + x + 8), result2);
      vst1q_f32(reinterpret_cast<float*>(p1 + x + 12), result3);
    }

    // Leftover value
    for (int x = wMod16 / sizeof(float); x < width; x++) {
      auto new_mask = has_mask ? reinterpret_cast<const float*>(mask)[x] * opacity_f : opacity_f;
      auto p1x = reinterpret_cast<float*>(p1)[x];
      auto p2x = reinterpret_cast<const float*>(p2)[x];
      auto result = p1x + (p2x - p1x) * new_mask;
      reinterpret_cast<float*>(p1)[x] = result;
    }

    p1 += p1_pitch;
    p2 += p2_pitch;
    if constexpr (has_mask)
      mask += mask_pitch;
  }
}

// instantiate
template void overlay_blend_neon_float<false>(BYTE* p1, const BYTE* p2, const BYTE* mask,
  const int p1_pitch, const int p2_pitch, const int mask_pitch, const int width, const int height, const int /*opacity*/, const float opacity_f, const int bits_per_pixel);
template void overlay_blend_neon_float<true>(BYTE* p1, const BYTE* p2, const BYTE* mask,
  const int p1_pitch, const int p2_pitch, const int mask_pitch, const int width, const int height, const int /*opacity*/, const float opacity_f, const int bits_per_pixel);

void overlay_darken_neon(BYTE* p1Y, BYTE* p1U, BYTE* p1V, const BYTE* p2Y, const BYTE* p2U, const BYTE* p2V, int p1_pitch, int p2_pitch, int width, int height) {
  overlay_darklighten_neon<overlay_darken_neon_cmp, overlay_darken_c_cmp>(p1Y, p1U, p1V, p2Y, p2U, p2V, p1_pitch, p2_pitch, width, height);
}
void overlay_lighten_neon(BYTE* p1Y, BYTE* p1U, BYTE* p1V, const BYTE* p2Y, const BYTE* p2U, const BYTE* p2V, int p1_pitch, int p2_pitch, int width, int height) {
  overlay_darklighten_neon<overlay_lighten_neon_cmp, overlay_lighten_c_cmp>(p1Y, p1U, p1V, p2Y, p2U, p2V, p1_pitch, p2_pitch, width, height);
}
