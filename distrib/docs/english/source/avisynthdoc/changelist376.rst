Changes from 3.7.5 to 3.7.6
---------------------------

Additions, changes
~~~~~~~~~~~~~~~~~~
- Bump year to 2026
- Overlay: "add" and "subtract" direct RGB mode
- Overlay: "add" and "subtract" supports 32-bit float.
- Added utf8 parameter to AddAutoLoadDir
- Added utf8 parameter to ListAutoLoadDirs
- Added utf8 parameter to LoadPlugin
- Added utf8 parameter to DumpFilterGraph
- "Info": Optimize AVX512 features display, group features, make a bit more compact
- "Info": add L2 cache size display
- "SetMaxCPU": add "avx512base" and "avx512fast" options to enable/disable AVX512 grouped features.
  see :ref:`SetMaxCPU <setmaxcpu>` . Users of base-only AVX512 CPUs can enable Avisynth-optimizations
  with SetMaxCPU("avx512base+").
- ARM64 (aarch64) area:

  * "Info": add ARMV8-A features display (NEON, DOTPROD, SVE2)

  * Add ArmV8-A cpu feature detection (NEON, DOTPROD, SVE2) on ARM64 Windows/Linux/macOS builds.
    On Windows, only up-to DOTPROD can be detected due to OS limitations.
  * New CPU flags in ``cpuid.h`` and ``avisynth_c.h``: CPUF_ARM_NEON, CPUF_ARM_DOTPROD, CPUF_ARM_SVE2
  * "SetMaxCPU": add "neon", "dotprod", "sve2" options to enable/disable ARM64 (aarch64) features.
- "ConvertToPlanarRGB": ``bits`` parameter: on-the-fly bit-depth conversions to YUV->RGB conversion.
  See :doc:`ConvertToPlanarRGB <./corefilters/convert>`.
- "ResetMask": add parameter float "opacity"
- "AddAlphaPlane": add parameter float "opacity"
- "Layer": YUY2 is handled as YV16 (lessen source code bloat)
- "Histogram" Color and Color2 mode additions and fixes:

  * Added ``matrix`` parameter (Color and Color2): specifies the YUV matrix for
    chroma interpretation (BT.601, BT.709, BT.2020, etc.), following the same
    convention as the ``ConvertToYUV`` family. If not set, the matrix is
    read from the clip's ``_Matrix`` and ``_ColorRange`` frame properties.
  * Added ``graticule`` string parameter (Color and Color2): controls the
    danger zone shading (Color) or valid chroma boundary square (Color2).
    ``"on"`` (default) always draws it, preserving pre-3.7.6 behavior;
    ``"off"`` never draws it; ``"auto"`` draws it only for limited-range
    clips and suppresses it for full-range clips.
  * Added ``circle`` parameter (Color and Color2, default false for Color,
    default true for Color2): draws the hue circle with 15° tick marks.
    In Color2 mode this was previously always drawn; it can now be disabled.
  * Added ``targets`` parameter (Color and Color2, default false): draws
    target boxes at the six 75%-amplitude ColorBars Cb/Cr positions
    (Yellow, Cyan, Green, Magenta, Red, Blue). Positions are computed from
    ground-truth linear RGB values through the active matrix, giving
    accurate coordinates at all bit depths.
  * Added ``axes`` parameter (Color and Color2, default false): draws
    horizontal and vertical crosshair lines through the vectorscope center.
  * Added ``iq`` parameter (Color and Color2, default false): draws target
    boxes for the NTSC −I, +I and +Q chroma-phase references, using the
    luma-corrected broadcast convention (Y raised until the most negative
    RGB component reaches Code 0, avoiding illegal RGB values).
  * Added ``iq_lines`` parameter (Color and Color2, default false): draws
    radial lines at the fixed NTSC subcarrier phase angles of 33°, 123°,
    213° and 303°, marking the I and Q chroma axes.
  * Fix: "Color" and "Color2" modes: copy alpha channel from source for
    alpha-carrying formats (YUVA, RGBPA, RGB32, RGB64); initialize alpha
    to zero in the histogram panel area.
  * Fix: accurate pixel positioning and scaling throughout the histogram
    panel, limited/full range aware. Marker positions use the same
    conversion path as the vectorscope signal dots, ensuring exact
    alignment at any bit depth and color range.
  * Vectorscope is now fully matrix-aware: all overlay marker positions are
    computed from ground-truth linear RGB values converted through the active
    YUV matrix, giving accurate Cb/Cr coordinates at all bit depths including
    32-bit float.


Build environment, Interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- introduce ``AVS_RESTRICT`` to ``avs/config.h`` (compiler invariant c++ ``__restrict``)
- AVX512: CMake to recognize ``*_avx512b.*`` and ``*_avx512.*`` file pattern, add compiler specific AVX512 
  compile flags accordingly (AVX512 Base and Ice Lake extensions)
- AVX512 support by conditional define.
  Define `INTEL_INTRINSICS_AVX512` if avx512 modules are enabled and compiler supports it.
  For MSVC,AVX512 support enabled only from MSVC 2019 16.2 (19.22) or newer.
- add ``.editorconfig``, update .gitignore to include the new .slnx format of Visual Studio 2026
- v12 interface: Global Lock support (https://github.com/AviSynth/AviSynthPlus/issues/444), 
  mainly for plugins using common fftw3 library:

  * ``env->AcquireGlobalLock``, ``env->ReleaseGlobalLock`` (C++),
  * ``avs_acquire_global_lock``, ``avs_release_global_lock`` (C)

  see :ref:`global lock support<cplusplus_acquiregloballock>`
- v12 interface: ApplyMessageEx supporting utf8 parameter.
  see :ref:`ApplyMessageEx<cplusplus_applymessage>`
- v12 interface: inform plugins about the effective thread count after Prefetch()
  via cache hints:

  * ``CachePolicyHint::CACHE_INFORM_NUM_THREADS`` (C++)
  * ``AVS_CACHE_INFORM_NUM_THREADS`` (C)

  See :ref:`SetCacheHints<cplusplus_setcachehints>` .
- the internal IScriptEnvironment2 methods AddAutoLoadDir and ListAutoLoadDirs explicitely
  work in UTF-8.
- New CPU flags: ``cpuid.h and ``avisynth_c.h``
  - added AVX512 group feature flags CPUF_AVX512_BASE and CPUF_AVX512_FAST (Ice Lake, usable AVX-512 since that point).
  - added many new AVX512 individual feature flags
  - added ARM64 feature flags CPUF_ARM_NEON, CPUF_ARM_DOTPROD, CPUF_ARM_SVE2
  - CPUF_xxxxx flags are now 64 bit, replace enum with constexpr.
- CMakeLists.txt: avx512 compile flag support for gcc/clang ("base" and "fast", latter is Ice Lake-like feature set).
- V12 interface: ``GetCPUFlagsEx`` returning 64 bit flags (too many AVX512 subfeatures to fit in 32 bit).
  C interface: ``avs_get_cpu_flags_ex``.
  see :ref:`GetCPUFlagsEx<cplusplus_getcpuflagsex>` and :ref:`GetCPUFlags<cplusplus_getcpuflags>`
- V12 interface: L2 cache size query support. New entry in ``AvsEnvProperty``: ``AEP_CACHESIZE_L2`` (C++), 
  ``AVS_AEP_CACHESIZE_L2`` (C) to query L2 cache size in bytes with ``IScriptEnvironment->GetEnvProperty()``. 
  x86/x64 architecture only for now. See :ref:`AvsEnvProperty<cplusplus_getenvproperty>` .

- Refactor CMakeLists.txt: 

  * Correct default of ``ENABLE_INTEL_SIMD`` for cross-compiling scenarios (e.g. ``ARM64`` target on ``x86_64`` host)
    Old logic relied on the host processor: ``${CMAKE_SYSTEM_PROCESSOR}``

  * Add back option to compile ``ARM64`` builds with Visual Studio on Windows. On VS2026 even clangcl (LLVM) is supported 
    out-of-box for ARM64 platform, in an easily cross-compilable way from an x64 machine.
  * VDubFilter: allow building on Windows only x86/x64 targets (and not for ARM64).
  * Fix LLVM/clangcl/Intel ICX compile warning: ``'WIN32' macro redefined as "#define WIN32 /D_WINDOWS /W3 /GR /EHsc 1 "``,
    when CMake injects a command-line macro wrongly and thus redefines WIN32 .
    The fix: converts global ``add_definitions("/D ...")`` and other option string magics into per-target ``target_compile_definitions()``
    and ``target_compile_options()``. Thus removing the accidental injection of ``${CMAKE_CXX_FLAGS}`` 
    into ``add_compile_options()``, and prevents the WIN32 macro redefinition.


Bugfixes
~~~~~~~~
- Fix: "Histogram" Color2 mode to copy alpha channel from source for alpha-carrying formats
  (YUVA, RGBPA, RGB32, RGB64); initialize alpha to zero in the histogram panel area.
  (Was: garbage)
- Fix: C-only vertical resampling code added more rounding than needed
  (regression since pre-3.7.5 20250427).
- Fix: "Invert": corrected chroma inversion to pivot around signed 0 instead
  of XORing with max_pixel_value.
- Fix: YUV->RGB limited range matrix accuracy for 10-16 bits, plus use a symmetric rounding in matrix 
  coefficient's integer approximation.
- Fix: inaccurate ColorBars 10+ bit values. Now they are derived from the 32-bit float 
  RGB definitions instead of upscaling a 8 bit precalculated YUV value. -I and +Q are still kept at
  legacy Avisynth values.
- Fix: inaccurate ColorBarsHD 10+ bit values. Now they are derived from the 32-bit float 
  RGB definitions instead of upscaling a 8 bit precalculated YUV value. Add 100% White after 
  Ramp section.
- Fix: GreyScale + SSE2 + RGB32 + matrix="RGB" overflow. 
  Rare usage; "RGB" matrix (Identity) uses a 1.0 coefficient which exceeds the signed 16-bit 
  SIMD limit of 32767 at 15-bit precision. Added bounds checking to fallback to C-code for any 
  coefficients >= 1.0 or < −1.0.
- Fix #448: Resolved an issue where MT_MULTI_INSTANCE filters using relative paths 
  (e.g. "video.mp4" or "../image.png") failed under Prefetch() when used in imported 
  scripts from different directories. The problem occurred because new thread instances did 
  not inherit the original working directory, causing path resolution to fail.
  Now, the current directory is captured at filter instantiation and passed to worker threads, 
  ensuring consistent path resolution.
- Fix #456: "Reverse" corrupts 24-bit audio (https://github.com/AviSynth/AviSynthPlus/issues/456)
- Fix BDF font rendering when it contains variable width characters like mixed Latin and CJK. 
  Preparing feature request #446 (https://github.com/AviSynth/AviSynthPlus/issues/446)
- Fix #462: Report: "AviSynth scripts don't work in a folder with a Unicode name."
  Plugin autoload folders are internally stored in UTF-8, regardless of which Windows ANSI codepage is set.

  * Folder names used in macros in AddAutoLoadDir (SCRIPTDIR, MAINSCRIPTDIR, PROGRAMDIR) and no longer 
    restricted to contain ANSI-only characters
  * Registry-backed macros that can contain plugin folder paths (USER_PLUS_PLUGINS, MACHINE_PLUS_PLUGINS, 
    USER_CLASSIC_PLUGINS, MACHINE_CLASSIC_PLUGINS) are read in Unicode friendly way as well. 
- Fix: Not existing registry entries won't appear as a macro string in auto-load path.
  E.g. Avisynth would automatically add ``USER_CLASSIC_PLUGINS`` at the beginning, but if no such entry
  exist, it kept being in the folder list as ``<current_directory>\USER_CLASSIC_PLUGINS\``. 
  Now this false entry is removed.
- Fix: Overlay give proper error message if 32-bit float is not supported in that mode.
- Change video-framebuffer over-allocation from 16 to 64 bytes. Allocate 64 bytes more than needed for 
  video frame buffer in order to be able to read 64 bytes safely with AVX512 without risking access violation 
  on the last pixels of the frame.
- Fix: The `Animate()` function now explicitly clamps interpolated values to ensure they remain 
  strictly between the start and end range. Due to the high precision of 64-bit `double` introduced 
  in v3.7.5, intermediate calculations could slightly exceed the boundary (e.g., 360.00000000000006 
  when interpolating from 0 to 360.0 in 564 steps), requiring this clamp to prevent out-of-range errors.


Optimizations
~~~~~~~~~~~~~
- "Layer" YUVA/YUV/RGBP "add", "subtract", "mul", "darken", "lighten": refactor
  chroma placement calculation and function dispatchers; main algorithm SIMD
  vectorization is now possible for non-444 chroma placements by precalculating
  a mask for the actual row. Add AVX2 path (LLVM/clangcl recommended).
- "Overlay" Blend: improved speed while keeping accuracy; use float arithmetic
  only where strictly needed.
- "Invert": planar formats no longer pre-copy all planes from source before
  conversion; only planes that are unchanged are copied, avoiding unnecessary work.
- TurnLeft, TurnRight: AVX2 support (1,5-3x speed on i7-11700 compared to SSE2 version)
- Turn180 AVX2 support (very slight speed gain)
- Resamplers: 

  * introduce a SIMD-like C header (avs_simd_c.h) for smart auto-vectorizing compilers.
  * restore vertical float performance (3.7.4 was slower than 3.7.3) + SSE2 special optimization
  * further optimize verticals, use ``AVS_RESTRICT``
  * (quicker RGB32/64 horizontal on AVX2 since TurnRigh/Left was optimized - packed RGB H-resize = TurnLeft-V-Resize-TurnRight)
  * optimize SSSE3 and AVX2 horizontal resampler for 32-bit float for small (<=4) kernel sizes
  * optimize 32-bit float vertical avx2
  * add AVX512 code path 
  
    - 32-bit float resamplers, verticals; horizontals up to kernel size 16.
    - 8-16-bit horizontal resamplers, for kernel size <= 16 and specific ratios; speed gain up to 300%+ (DTL2020)!!
    - 8-16-bit vertical resamplers
  * (Work In Progress) unify horizontal and vertical plane processing flow

- add NEON optimizations for ARM64 (aarch64) for TurnLeft/TurnRight/Turn180.
  (First aarch64 code in Avisynth)
- add NEON optimizatons for Overlay blend
- (filter graph) avoid MTGuard and CacheGuard creation if a filter returns one of its clip parameter unaltered.
- Add some avx2 stuff to Invert (no really gain, filter is too simple) and some Layer subfilter.

Documentation
~~~~~~~~~~~~~
- Build on Raspbian, Raspberry Pi 5 and llvm/gcc :doc:`This page <./contributing/posix>` 
  describes linux builds process.
- Add to :ref:`SIMD in ARM64 (aarch64) section <aarch64_simd_tiers>`
- Extend ``env->Allocate/Free`` see at :ref:`Allocate <cplusplus_allocate>`
- Interface V12 changes: see :ref:`API v12 changes<api_v12_whats_new>` for more details.
- Add folder macro description to AddAutoLoadPlugins
- Update :doc:`Overlay <./corefilters/overlay>`
- Update :ref:`SetCacheHints<cplusplus_setcachehints>` with ``CACHE_INFORM_NUM_THREADS``
- Update :ref:`GetCPUFlags<cplusplus_getcpuflags>`, add :ref:`GetCPUFlagsEx<cplusplus_getcpuflagsex>`
- Update :ref:`CPU Feature Flags<cplusplus_cpufeatureflags>` with AVX512 and ARM64 features
- Update :ref:`SetMaxCPU <setmaxcpu>` with AVX512 and ARM64 features
- Update :ref:`AvsEnvProperty<cplusplus_getenvproperty>` with L2 cache size entry
- Update Russian GPL notice in UTF-8 format
- Update :doc:`ConvertToPlanarRGB <./corefilters/convert>` with `"bits"` and "matrix" syntax `":same"`
- Update :doc:`ResetMask <./corefilters/mask>` with `"opacity"` and additional insights
- Update :doc:`AddAlphaPlane <./corefilters/mask>` with `"opacity"` and additional insights
- Update :doc:`Overlay <./corefilters/overlay>` with "add" and "subtract"
  direct RGB mode and 32-bit float support.
- Update :doc:`Layer <./corefilters/layer>` with "use_chroma" and opacity
  details, and YUY2/YV16 internal handling.
- Update :doc:`Histogram <./corefilters/histogram>` with new vectorscope parameters
- Update :doc:`ColorBars <./corefilters/colorbars>`
- Add another Ubuntu->Windows DLL cross-compilation guide:
  See :ref:`Ubuntu->Windows mingw crosscompilation<compiling_avsplus_crosscompiling2>`


Please report bugs at `github AviSynthPlus page`_ - or - `Doom9's AviSynth+
forum`_

$Date: 2026/02/14 23:42:00 $

.. _github AviSynthPlus page:
    https://github.com/AviSynth/AviSynthPlus
.. _Doom9's AviSynth+ forum:
    https://forum.doom9.org/showthread.php?t=181351
