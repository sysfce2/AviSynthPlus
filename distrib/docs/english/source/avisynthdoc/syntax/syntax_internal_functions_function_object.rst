
AviSynth Syntax - Cast to function object
=========================================

Func
~~~~
::

    Func( script_function )

Converts an internal function to a valid function object.

For some functions (e.g. ScriptClip) it is convenient to use their "function object" parameter version.
With this function an internal function can be converted into it.

See also: :doc:`Function objects <syntax_function_objects>`

::

    # write frame properties with function object
    ScriptClip("""propSetInt("frameprop_from_str",func(YPlaneMax))""")
    # write frame properties with traditional script string
    ScriptClip(function[](clip c) { propSetInt("frameluma_sc_func",func(AverageLuma)) })

Changelog
~~~~~~~~~
+----------------+------------------------------------------------------------+
| Version        | Changes                                                    |
+================+============================================================+
| Avisynth 3.6.0 | Added function objects                                     |
+----------------+------------------------------------------------------------+

--------

Back to :doc:`Internal functions <syntax_internal_functions>`.


$Date: 2024/01/09 11:30:00 $

.. _planar: http://avisynth.org/mediawiki/Planar
.. _memory alignment used in the AVIFile output emulation (not yet written):
    http://avisynth.org/mediawiki/index.php?title=AVIFile_output_emulation
