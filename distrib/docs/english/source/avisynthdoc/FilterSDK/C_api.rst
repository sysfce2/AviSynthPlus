
C API
=====

.. toctree::
    :maxdepth: 3

.. contents:: Table of contents


The header, avisynth_c.h and some of its helpers in avs/ folder, declares all the classes, structures and
miscellaneous constants that you might need when writing a plugin or a client. All
external plugins should #include it:
::

    #include "avisynth_c.h"

or if proper paths are set to the installed package and SDK include files:

    #include <avisynth_c.h>

Useful source links
~~~~~~~~~~~~~~~~~~~

Note from 2025: until this part is updated properly, check these excellent
examples for using C API as a client or in a plugin:

- ffmpeg (client)

  https://github.com/FFmpeg/FFmpeg/blob/release/7.1/libavformat/avisynth.c

- x265mod (client)

  https://github.com/Patman86/x265-Mod-by-Patman/blob/master/source/input/avs.cpp

- x264 (client)

  https://code.videolan.org/videolan/x264/-/blob/master/input/avs.c?ref_type=heads
  
- avs2yuv (client)

  https://github.com/DJATOM/avs2yuv
  
- Many-many plugins from Asd-g:

  https://github.com/Asd-g/
  
- AvsInPaint plugin

  https://github.com/pinterf/AvsInpaint
  
- AssRender plugin

  https://github.com/pinterf/assrender


.. _c_avs_scriptenvironment:

AVS_ScriptEnvironment
~~~~~~~~~~~~~~~~~~~~~

.. _c_avs_add_function:

.. _c_avs_add_function_r:

avs_add_function
^^^^^^^^^^^^^^^^
avs_add_function_r
^^^^^^^^^^^^^^^^^^

::

    int avs_add_function(AVS_ScriptEnvironment *, 
                         const char * name, const char * params,
                         AVS_ApplyFunc apply, void * user_data);

    int avs_add_function_r(AVS_ScriptEnvironment *, 
                         const char * name, const char * params,
                         AVS_ApplyFuncR apply, void * user_data);


Both forms define the function name, parameter signature, and the callback function itself. The 
difference lies in the type of the callback function (apply).

``avs_add_function`` and ``avs_add_function_r`` are used to inform AviSynth of the existence of
our filter. These functions register a function with AviSynth's internal function table.

The base ``avs_add_function`` takes four arguments: the name of the new script function, the 
parameter-type string, the C function (callback) implementing the script function, and the 
user_data cookie.

The added function returns a type ``AVS_Value`` and can therefore return any AVS_Value type, Clip, string, 
integer, double, etc.. In this version the function returns AVS_Value directly.

The second form, ``avs_add_function_r``, is an alternative approach where the function result is 
provided by filling the ``AVS_Value`` result into a passed pointer.

This is particularly useful when interfacing with Python. In Python 3.13 (as of 2025), our callback
written in Python cannot return structs (like ``AVS_Value``) directly to the C caller via a function 
return value. However, it can accept and fill an ``AVS_Value`` C struct passed as a pointer.

The first form avs_add_function, the callback (``apply``) returns result as return value (``AVS_Value``).
This is the callback type used by ``avs_add_function``:
::

    typedef AVS_Value (AVSC_CC * AVS_ApplyFunc)(AVS_ScriptEnvironment *, 
                                                AVS_Value args, void * user_data);

Int the alternative form, the callback (``apply``) returns result in byref parameter (``AVS_Value *``)
This is the callback type used by ``avs_add_function_r``:
::

    typedef void(AVSC_CC* AVS_ApplyFuncR)(AVS_ScriptEnvironment*, 
                                          AVS_Value* ret, AVS_Value args, void* user_data);

Its main purpose in ``avisynth_c_plugin_init`` or ``avisynth_c_plugin_init2`` to add and define plugin filters.
But a client or a plugin can define its own functions or filters as well.

Example:
::

    static AVS_Value AVSC_CC Create_JincResize(AVS_ScriptEnvironment* env, AVS_Value args, void* param) {
    ...
    }

    avs_add_function(env, "JincResize", "cii[src_left]f[src_top]f[src_width]f[src_height]f[quant_x]i[quant_y]i"
                                        "[tap]i[blur]f[cplace]s[threads]i[opt]i", Create_JincResize, 0);

Example (alternative version): a simple (non-Clip oriented) example which would add X to the input.
::

    static void AVSC_CC Create_IncreaseBy(AVS_ScriptEnvironment* env, AVS_Value *retval, AVS_Value args, void* param) {
    ...
    }

    avs_add_function_r(env, "IncreaseBy", "[delta]i", Create_IncreaseBy, 0);


The added function is of type AVSValue and can therefore return any AVSValue. 

For more info and examples see also :ref:`cplusplus_addfunction` in C++ API.


Historical content
~~~~~~~~~~~~~~~~~~

* to be checked and updated *

source: http://forum.doom9.org/showthread.php?p=1464911#post1464911

compiling plugins: http://forum.doom9.org/showthread.php?p=1092380#post1092380

example: http://forum.doom9.org/showthread.php?p=1001260#post1001260

| For now the api is described here:
| http://www.kevina.org/avisynth_c/api.html.

| An example is given here:
| http://www.kevina.org/avisynth_c/example.html.

In v2.60 (AVISYNTH_INTERFACE_VERSION = 6) the following functions are
added to the C interface:
::

    avs_is_yv24
    avs_is_yv16
    avs_is_yv12
    avs_is_yv411
    avs_is_y8
    avs_is_color_space
    avs_get_plane_width_subsampling
    avs_get_plane_height_subsampling
    avs_bits_per_pixel
    avs_bytes_from_pixels
    avs_row_size
    avs_bmp_size
    avs_get_row_size_p
    avs_get_height_p

____

Back to :doc:`FilterSDK`

$Date: 2025/02/24 13:53:00 $
