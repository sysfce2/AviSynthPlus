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
  see :ref:`SetMaxCPU <setmaxcpu>` .
- ARM64 (aarch64) area:

  * "Info": add ARMV8-A features display (NEON, DOTPROD, SVE2)

  * Add ArmV8-A cpu feature detection (NEON, DOTPROD, SVE2) on ARM64 Windows/Linux/macOS builds.
    On Windows, only up-to DOTPROD can be detected due to OS limitations.
  * New CPU flags in ``cpuid.h`` and ``avisynth_c.h``: CPUF_ARM_NEON, CPUF_ARM_DOTPROD, CPUF_ARM_SVE2
  * "SetMaxCPU": add "neon", "dotprod", "sve2" options to enable/disable ARM64 (aarch64) features.


Build environment, Interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- introduce ``AVS_RESTRICT`` to ``avs/config.h`` (compiler invariant c++ ``__restrict``)
- AVX512: CMake to recognize ``*_avx512.*`` file pattern, add compiler specific AVX512 
  compile flags accordingly (AVX512 Base and Ice Lake extensions)
- AVX512 support by conditional define.
  Define `INTEL_INTRINSICS_AVX512` if avx512 modules are enabled and compiler supports it.
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
- CMakeLists.txt: avx512 compile flag support for gcc/clang ("fast" Ice Lake-like feature set).
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
- Resamplers: 

  * introduce a SIMD-like C header (avs_simd_c.h) for smart auto-vectorizing compilers.
  * restore vertical float performance (3.7.4 was slower than 3.7.3) + SSE2 special optimization
  * further optimize verticals, use ``AVS_RESTRICT``
  * optimize SSSE3 and AVX2 horizontal resampler for 32-bit float for small (<=4) kernel sizes
  * optimize 32-bit float vertical avx2
  * add AVX512 code path to 32-bit float resamplers
  * (Work In Progress) unify horizontal and vertical plane processing flow

- add NEON optimizations for ARM64 (aarch64) for TurnLeft/TurnRight/Turn180.
  (First aarch64 code in Avisynth)

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


Please report bugs at `github AviSynthPlus page`_ - or - `Doom9's AviSynth+
forum`_

$Date: 2026/01/01 10:42:00 $

.. _github AviSynthPlus page:
    https://github.com/AviSynth/AviSynthPlus
.. _Doom9's AviSynth+ forum:
    https://forum.doom9.org/showthread.php?t=181351
