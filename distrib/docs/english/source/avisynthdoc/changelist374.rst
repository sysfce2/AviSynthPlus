
Changes.
========


Changes from 3.7.3 to 3.7.4
---------------------------

Additions, changes
~~~~~~~~~~~~~~~~~~
- DirectShowSource new parameter ``utf8`` for utf8 filename support
- "propShow" ``font``, ``text_color``, ``halo_color``, ``bold`` new parameters for custom style
- "propShow" (#366): ``x``, ``y``, ``align`` new parameters for custom positioning
- "Info": ``cpu`` new parameter to disable showing CPU capabilities
- "Info" (#366): ``x``, ``y``, ``align`` new parameters for custom positioning
- Fix #368 Make proper vertical alignment for multiline text in Subtitle and Text 
  when vertical alignment is set to bottom or center.
- Studio RGB (narrow, limited) range will now be recognized (through _ColorRange=1)
  and utilized in conversions from RGB, such as in GreyScale, ConvertToY, ConvertToYUVxxx

Build environment, Interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Bugfixes
~~~~~~~~
- "Info" now can display a line which is only partially visible (instead of not showing it at all)
- "Text" use "lsp" parameter the same way as in SubTitle: in 1/8 pixel units, not in 1 pixels.
- "Text" vertical alignment position would be wrong for multiline strings containing even number of lines.
- Fix #365: Avisynth 2.5 plugins when NICE_FILTER would crash with "invalid response to CACHE_GETCHILD_AUDIO_MODE".
  Regression in 3.7.3 reintroduced audio cache.
- Fix #370: array size assert error in ConvertToYUY2 when internally ConvertToYUV422 is called
- Leave _ColorRange frame property as-is, when using matrix names "PC.709" or "PC.601", for example in ConvertToRGB32.

Optimizations
~~~~~~~~~~~~~

Documentation
~~~~~~~~~~~~~
- Correct building DirectShowSource prerequisites (Release_MBCS)
- Update "DirectShowSource" with utf8 parameter
- Update "Info"
- Update rst docs with control structs if/else/for/while
- Update "ShowTime", "ShowSMPTE", "ShowFrameNumber" section with 3.7.3 changes

Please report bugs at `github AviSynthPlus page`_ - or - `Doom9's AviSynth+
forum`_

$Date: 2023/12/02 18:54:00 $

.. _github AviSynthPlus page:
    https://github.com/AviSynth/AviSynthPlus
.. _Doom9's AviSynth+ forum:
    https://forum.doom9.org/showthread.php?t=181351
