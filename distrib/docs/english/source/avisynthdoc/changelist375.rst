Changes from 3.7.4 to 3.7.5
---------------------------

Additions, changes
~~~~~~~~~~~~~~~~~~

Build environment, Interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Bugfixes
~~~~~~~~
- Fix #434: YtoUV crash, regression in 3.7.4 (since 61d2c9a)
- Fix: crash in (Horizontal) resizers on non-Intel (e.g. aarch64) platforms

Optimizations
~~~~~~~~~~~~~
- (preliminary) convert source of C only horizontal and vertical resizers to 
  auto-vectorize compiler friendly (e.g. quicker on aarch64 until neon is supported)

Documentation
~~~~~~~~~~~~~



Please report bugs at `github AviSynthPlus page`_ - or - `Doom9's AviSynth+
forum`_

$Date: 2025/04/14 14:10:00 $

.. _github AviSynthPlus page:
    https://github.com/AviSynth/AviSynthPlus
.. _Doom9's AviSynth+ forum:
    https://forum.doom9.org/showthread.php?t=181351
