
Changes.
========


Changes from 3.7.3 to 3.7.4
---------------------------

Additions, changes
~~~~~~~~~~~~~~~~~~
- Bump year to 2025
- DirectShowSource new parameter ``utf8`` for utf8 filename support
- "propShow" ``font``, ``text_color``, ``halo_color``, ``bold`` new parameters for custom style
- "propShow" (#366): ``x``, ``y``, ``align`` new parameters for custom positioning
- "Info": ``cpu`` new parameter to disable showing CPU capabilities
- "Info" (#366): ``x``, ``y``, ``align`` new parameters for custom positioning
- Fix #368 Make proper vertical alignment for multiline text in Subtitle and Text 
  when vertical alignment is set to bottom or center.
- Studio RGB (narrow, limited) range will now be recognized (through _ColorRange=1)
  and utilized in conversions from RGB, such as in GreyScale, ConvertToY, ConvertToYUVxxx
- #392 "break" and "continue" in for-next and while loops
- Add "ArraySort" for sorting simple bool, numeric or string arrays
- SoundTouch: Increased the max number of channels to 32ch from the old limit of 16ch.
- Posix: Detect additional plugindirs from LD_LIBRARY_PATH environment variable
- #413 Add ListAutoloadDirs() script function returing a LF delimited string with directory list.
- SubTitle to accept real LF (``\r``) or CR LF (``\r\n``) control characters for line break.

Build environment, Intenface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- add LOONGARCH support
- Use system installs of DevIL and SoundTouch on all platforms, remove in-tree binaries/code
- avisynth.h: add ListAutoloadDirs() to internal interface declarations
- CMakeList.txt to accept Intel C++ Compiler 2025

Bugfixes
~~~~~~~~
- Fix "SetLogParams" defaults - mentioned in #391
- Fix corrupt Turn functions when a planar RGB turn would be followed by a YUV Turn.
  Regression since TurnXXXX supports planar RGB (2016.08.23; probably since r2081 commit dba954e2de0c9c6218d17fc5c4974f4c28b627c3)
- Fix #386: Interleave to call plugin destructor like StackXXXX (memory leak in case if script errors)
- Fix #384: swapped ShowGreen/ShowBlue for planar RGB sources
- Fix: allow use of "local" in ConditionalSelect string version (fixed wrong function signature)
- "Info" now can display a line which is only partially visible (instead of not showing it at all)
- "Text" use "lsp" parameter the same way as in SubTitle: in 1/8 pixel units, not in 1 pixels.
- "Text" vertical alignment position would be wrong for multiline strings containing even number of lines.
- Fix #365: Avisynth 2.5 plugins when NICE_FILTER would crash with "invalid response to CACHE_GETCHILD_AUDIO_MODE".
  Regression in 3.7.3 reintroduced audio cache.
- Fix #370: array size assert error in ConvertToYUY2 when internally ConvertToYUV422 is called
- Leave _ColorRange frame property as-is, when using matrix names "PC.709" or "PC.601", for example in ConvertToRGB32.
- Fix FadeIn, which behaved like FadeIn2


Optimizations
~~~~~~~~~~~~~
- avoid storing duplicated strings in internal string heap (related to issue #389)

Documentation
~~~~~~~~~~~~~
- Correct building DirectShowSource prerequisites (Release_MBCS)
- Update "DirectShowSource" with utf8 parameter
- Update "Info"
- Update rst docs with control structs if/else/for/while
- Update "ShowTime", "ShowSMPTE", "ShowFrameNumber" section with 3.7.3 changes
- Update most items at Syntax and internal functions sections, add arrays, function objects, 
  escaped string literals, multithreading, frame properties, debug functions
- Add if-else, do-while, for-next, break and continue
- Update Import (add utf8)
- update Conditional filters, Runtime functions
- update ShowAlpha/Red/...
- Update Intel 2024 compiler build process
- Add autoload helper functions to rst docs, update content and bring syntax_plugins.rst to new format
- update compile_avsplus.rst online docs with Intel C++ Compiler 2025

Please report bugs at `github AviSynthPlus page`_ - or - `Doom9's AviSynth+
forum`_

$Date: 2025/01/01 9:42:00 $

.. _github AviSynthPlus page:
    https://github.com/AviSynth/AviSynthPlus
.. _Doom9's AviSynth+ forum:
    https://forum.doom9.org/showthread.php?t=181351
