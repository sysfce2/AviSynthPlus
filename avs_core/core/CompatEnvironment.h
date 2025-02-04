// Avisynth v2.5.  Copyright 2002-2009 Ben Rudiak-Gould et al.
// http://www.avisynth.org

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

#ifndef _AVS_COMPATENVIRONMENT_H_INCLUDED
#define _AVS_COMPATENVIRONMENT_H_INCLUDED

#include <avisynth.h>

// IScriptEnvironment_Avs25 is used internally.
// Entries 100% match with IScriptEnvironment to allow typecasting to IScriptEnvironment_Avs25.
// IScriptEnvironment_Avs25 ends at interface V6 changes.
// Differences involve three functions:
// 1. AddFunction25 : Marks the added function as originating from an Avs2.5 plugin, see also: arrays-deep-copy
// 2. Invoke25      : Argument arrays, which must not be deep-freed, baked AVSValue code incompatible with deep-array concept
// 3. ManageCache25 : Specially returns 1 for key MC_QueryAvs25 to check if called from AVS2.5 interface for cache GetFrame/GetAudio.

// See also: InternalEnvironment interface conversion methods comment.

// This compatibility camouflage class is *passed in place of the original IScriptEnvironment,
// providing the same function orders to match VMT table with the original.
// Differently named functions appear at specific indexes to perform different compatibility behaviors.
// *In real they all are static casted from InternalEnvironment, where all these pure virtual methods are implemented.

class IScriptEnvironment_Avs25 {
public:
  virtual ~IScriptEnvironment_Avs25() {}

  virtual /*static*/ int __stdcall GetCPUFlags() = 0;

  virtual char* __stdcall SaveString(const char* s, int length = -1) = 0;
  virtual char* Sprintf(const char* fmt, ...) = 0;
  virtual char* __stdcall VSprintf(const char* fmt, va_list val) = 0;

#ifdef AVS_WINDOWS
  __declspec(noreturn) virtual void ThrowError(const char* fmt, ...) = 0;
#else
  virtual void ThrowError(const char* fmt, ...) = 0;
#endif

  class NotFound /*exception*/ {};  // thrown by Invoke and GetVar

  typedef AVSValue(__cdecl* ApplyFunc)(AVSValue args, void* user_data, IScriptEnvironment* env);

  virtual void __stdcall AddFunction25(const char* name, const char* params, ApplyFunc apply, void* user_data) = 0;
  virtual bool __stdcall FunctionExists(const char* name) = 0;
  virtual AVSValue __stdcall Invoke25(const char* name, const AVSValue args, const char* const* arg_names = 0) = 0;

  virtual AVSValue __stdcall GetVar(const char* name) = 0;
  virtual bool __stdcall SetVar(const char* name, const AVSValue& val) = 0;
  virtual bool __stdcall SetGlobalVar(const char* name, const AVSValue& val) = 0;

  virtual void __stdcall PushContext(int level = 0) = 0;
  virtual void __stdcall PopContext() = 0;

  virtual PVideoFrame __stdcall NewVideoFrame(const VideoInfo& vi, int align = FRAME_ALIGN) = 0;

  virtual bool __stdcall MakeWritable(PVideoFrame* pvf) = 0;

  virtual void __stdcall BitBlt(BYTE* dstp, int dst_pitch, const BYTE* srcp, int src_pitch, int row_size, int height) = 0;

  typedef void(__cdecl* ShutdownFunc)(void* user_data, IScriptEnvironment* env);
  virtual void __stdcall AtExit(ShutdownFunc function, void* user_data) = 0;

  virtual void __stdcall CheckVersion(int version = AVISYNTH_CLASSIC_INTERFACE_VERSION_25) = 0;

  virtual PVideoFrame __stdcall Subframe(PVideoFrame src, int rel_offset, int new_pitch, int new_row_size, int new_height) = 0;

  virtual int __stdcall SetMemoryMax(int mem) = 0;

  virtual int __stdcall SetWorkingDir(const char* newdir) = 0;

  // specially returns 1 for key MC_QueryAvs25 to check if called from AVS2.5 interface
  virtual void* __stdcall ManageCache25(int key, void* data) = 0;

  enum PlanarChromaAlignmentMode {
    PlanarChromaAlignmentOff,
    PlanarChromaAlignmentOn,
    PlanarChromaAlignmentTest
  };

  virtual bool __stdcall PlanarChromaAlignment(IScriptEnvironment::PlanarChromaAlignmentMode key) = 0;

  virtual PVideoFrame __stdcall SubframePlanar(PVideoFrame src, int rel_offset, int new_pitch, int new_row_size,
    int new_height, int rel_offsetU, int rel_offsetV, int new_pitchUV) = 0;

  // Despite the name, we provide entries up to V6 in case someone requests
  // a V3 interface and still wants to use V5-V6 functions

  // **** AVISYNTH_INTERFACE_VERSION 5 **** defined since classic Avisynth 2.6 beta
  virtual void __stdcall DeleteScriptEnvironment() = 0;

  virtual void __stdcall ApplyMessage(PVideoFrame* frame, const VideoInfo& vi, const char* message, int size,
    int textcolor, int halocolor, int bgcolor) = 0;

  virtual const AVS_Linkage* __stdcall GetAVSLinkage() = 0;

  // **** AVISYNTH_INTERFACE_VERSION 6 **** defined since classic Avisynth 2.6
  // noThrow version of GetVar
  virtual AVSValue __stdcall GetVarDef(const char* name, const AVSValue& def = AVSValue()) = 0;

}; // end class IScriptEnvironment_Avs25. Order is important.


#endif // _AVS_COMPATENVIRONMENT_H_INCLUDED
