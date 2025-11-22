Changes from 3.7.5 to 3.7.6
---------------------------

Additions, changes
~~~~~~~~~~~~~~~~~~
- Overlay: "add" and "subtract" direct RGB mode
- Overlay: "add" and "subtract" supports 32-bit float.
- Added utf8 parameter to AddAutoLoadDir
- Added utf8 parameter to ListAutoLoadDirs
- Added utf8 parameter to LoadPlugin
- Added utf8 parameter to DumpFilterGraph

Build environment, Interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- introduce ``AVS_RESTRICT`` to ``avs/config.h`` (compiler invariant c++ ``__restrict``)
- AVX512: CMake to recognize ``*_avx512.*`` file pattern, add compiler specific AVX512 
  compile flags accordingly (foundation and bw extension assumed)
- AVX512 support by conditional define.
  Define `INTEL_INTRINSICS_AVX512` if avx512 modules are enabled 
  (The conditional is undefined for non-intel arch., 32 bit, or for pre MSVC 2019 16.2 (19.22))
- add ``.editorconfig``, update .gitignore to include the new .slnx format of Visual Studio 2026
- v12 interface: Global Lock support (https://github.com/AviSynth/AviSynthPlus/issues/444), 
  mainly for plugins using common fftw3 library:

  * ``env->AcquireGlobalLock``, ``env->ReleaseGlobalLock`` (C++),
  * ``avs_acquire_global_lock``, ``avs_release_global_lock`` (C)

  see :ref:`global lock support<cplusplus_acquiregloballock>`
- v12 interface: ApplyMessageEx supporting utf8 parameter.
  see :ref:`ApplyMessageEx<cplusplus_applymessage>`
- the internal IScriptEnvironment2 methods AddAutoLoadDir and ListAutoLoadDirs explicitely
  work in UTF-8.

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
- Fix: Overlay specific modes give proper error message if 32-bit float is not supported

Optimizations
~~~~~~~~~~~~~
- Resamplers: 

  * introduce a SIMD-like C header (avs_simd_c.h) for smart auto-vectorizing compilers.
  * restore vertical float performance (3.7.4 was slower than 3.7.3) + SSE2 special optimization
  * further optimize verticals, use ``AVS_RESTRICT``
  * (Work In Progress) optimize horizontal resampler for small (<=4) kernel sizes.
  * (Work In Progress) add AVX512 code path
  * (Work In Progress) unifify horizontal and vertical place processing flow

Documentation
~~~~~~~~~~~~~
- Build on Raspbian, Raspberry Pi 5 and llvm/gcc :doc:`This page <./contributing/posix>` 
  describes linux builds process.
- Extend ``env->Allocate/Free`` see at :ref:`Allocate <cplusplus_allocate>`
- Interface V12 changes: see :ref:`api_v12_whats_new` for more details.
- Add folder macro description to AddAutoLoadPlugins
- Update :doc:`Overlay <./corefilters/overlay>` 


Please report bugs at `github AviSynthPlus page`_ - or - `Doom9's AviSynth+
forum`_

$Date: 2025/11/18 12:53:00 $

.. _github AviSynthPlus page:
    https://github.com/AviSynth/AviSynthPlus
.. _Doom9's AviSynth+ forum:
    https://forum.doom9.org/showthread.php?t=181351
