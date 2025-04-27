Changes from 3.7.5 to 3.7.6
---------------------------

Additions, changes
~~~~~~~~~~~~~~~~~~


Build environment, Interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- introduce `AVS_RESTRICT` to `avs/config.h` (compiler invariant c++ __restrict)

Bugfixes
~~~~~~~~

Optimizations
~~~~~~~~~~~~~
- Resizers: introduce a SIMD-like C header (avs_simd_c.h) for smart auto-vectorizing compilers.
- Resizers: add back vertical float performance (3.7.4 was slower than 3.7.3) + SSE2 special optimization
- Resizers: further optimize verticals, use AVS_RESTRICT

Documentation
~~~~~~~~~~~~~
- Build on Raspbian, Raspberry Pi 5 and llvm/gcc :doc:`This page <./contributing/posix>` 
  describes linux builds process.


Please report bugs at `github AviSynthPlus page`_ - or - `Doom9's AviSynth+
forum`_

$Date: 2025/04/27 12:50:00 $

.. _github AviSynthPlus page:
    https://github.com/AviSynth/AviSynthPlus
.. _Doom9's AviSynth+ forum:
    https://forum.doom9.org/showthread.php?t=181351
