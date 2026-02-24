
ConvertToXXXX function
======================

*RGB interleaved*
::

  ConvertToRGB(clip [, string matrix, bool interlaced,
         string ChromaInPlacement,
         string chromaresample,
         float param1, float param2, float param3] )

  ConvertToRGB24(clip [, string matrix, bool interlaced,
         string ChromaInPlacement,
         string chromaresample,
         float param1, float param2, float param3] )

  ConvertToRGB32(clip [, string matrix, bool interlaced,
         string ChromaInPlacement,
         string chromaresample,
         float param1, float param2, float param3] )

  ConvertToRGB48(clip, [ string matrix, bool interlaced,
         string ChromaInPlacement,
         string chromaresample,
         float param1, float param2, float param3 ] )
       
  ConvertToRGB64(clip, [ string matrix, bool interlaced,
         string ChromaInPlacement,
         string chromaresample,
         float param1, float param2, float param3 ] ) 


*RGB planar*
::

    ConvertToPlanarRGB(clip, [ string matrix, bool interlaced,
         string ChromaInPlacement,
         string chromaresample,
         float param1, float param2, float param3,
         int bits, bool quality] )
    ConvertToPlanarRGBA(clip, [ string matrix, bool interlaced,
         string ChromaInPlacement,
         string chromaresample,
         float param1, float param2, float param3,
         int bits, bool quality] )


*YUV444, YUVA444*
::

    ConvertToYV24(clip [, bool interlaced, string matrix,
         string ChromaInPlacement,
         string chromaresample,
         float param1, float param2, float param3] )
    ConvertToYUV444(clip, [ string matrix, bool interlaced,
         string ChromaInPlacement,
         string chromaresample,
         float param1, float param2, float param3] )
    ConvertToYUVA444(clip, [ string matrix, bool interlaced,
         string ChromaInPlacement,
         string chromaresample,
         float param1, float param2, float param3] )


*YUV422, YUVA422*
::

    ConvertToYV16(clip [, bool interlaced, string matrix,
         string ChromaInPlacement,
         string chromaresample,
         string ChromaOutPlacement
         float param1, float param2, float param3] )
    ConvertToYUV422(clip, [ string matrix, bool interlaced,
         string ChromaInPlacement,
         string chromaresample,
         string ChromaOutPlacement,
         float param1, float param2, float param3] )
    ConvertToYUVA422(clip, [ string matrix, bool interlaced,
         string ChromaInPlacement,
         string chromaresample,
         string ChromaOutPlacement,
         float param1, float param2, float param3] )

*YUY2*
::

    ConvertToYUY2(clip [, bool interlaced, string matrix,
         string ChromaInPlacement,
         string chromaresample,
         float param1, float param2, float param3] )
    ConvertBackToYUY2(clip [, string matrix ] )


*YUV420, YUVA420*
::

    ConvertToYV12(clip [, bool interlaced, string matrix,
         string ChromaInPlacement,
         string chromaresample,
         string ChromaOutPlacement,
         float param1, float param2, float param3] )
    ConvertToYUV420(clip, [ string matrix, bool interlaced,
         string ChromaInPlacement,
         string chromaresample,
         string ChromaOutPlacement,
         float param1, float param2, float param3] )
    ConvertToYUVA420(clip, [ string matrix, bool interlaced,
         string ChromaInPlacement,
         string chromaresample,
         string ChromaOutPlacement,
         float param1, float param2, float param3] )


*YUV411*
::

    ConvertToYV411(clip [, bool interlaced, string matrix,
         string ChromaInPlacement,
         string chromaresample,
         float param1, float param2, float param3] )
    ConvertToYUV411(clip [, bool interlaced, string matrix,
         string ChromaInPlacement,
         string chromaresample,
         float param1, float param2, float param3] )
         
(the 2nd one is just an alias)

*Y-only*
::

    ConvertToY8(clip [, string matrix] )
    ConvertToY(clip, [ string matrix ] )


Color formats
-------------

The following formats can be converted to and from.

Notes:

- Interleaved RGB formats (RGB24/32/48/64) are kept for compatibility, they come with fixed 8 and 16 bits.
- The successor RGB format is planar RGB/RGBA, which support any Avisynth+ high bit depth format.
- 8 bit YUV formats has their own old names, but can be written in the generic naming convention: 
  YV12=YUV420P8, YV16=YUV422P8, YV24=YUV444P8
- If possible, avoid using YUY2 which is kept for compatibility. Use YV16 instead.

+----------------+-----------+--------------+---------------------------------------------------------------+-------------+
| Color formats  | Bit depth | Sample ratio | Description                                                   | planar/     |
|                |           |              |                                                               | interleaved |
+================+===========+==============+===============================================================+=============+
| RGB24, RGB48   | 8, 16     | 4:4:4        | full chroma                                                   | interleaved |
+----------------+-----------+--------------+---------------------------------------------------------------+-------------+
| RGB32, RGB64   | 8, 16     | 4:4:4:4      | full chroma + alpha                                           | interleaved |
+----------------+-----------+--------------+---------------------------------------------------------------+-------------+
| RGBPxx         | 8-16, 32  | 4:4:4        | full chroma - known as planar RGB                             | planar      |
+----------------+-----------+--------------+---------------------------------------------------------------+-------------+
| RGBAPxx        | 8-16, 32  | 4:4:4:4      | full chroma + alpha - known as planar RGBA                    | planar      |
+----------------+-----------+--------------+---------------------------------------------------------------+-------------+
| YV24,YUV444Pxx | 8-16, 32  | 4:4:4        | full chroma                                                   | planar      |
+----------------+-----------+--------------+---------------------------------------------------------------+-------------+
| YUVA444Pxx     | 8-16, 32  | 4:4:4:4      | full chroma + alpha                                           | planar      |
+----------------+-----------+--------------+---------------------------------------------------------------+-------------+
| YV16,YUV422Pxx | 8-16, 32  | 4:2:2        | chroma shared between 2 pixels                                | planar      |
+----------------+-----------+--------------+---------------------------------------------------------------+-------------+
| YUVA422Pxx     | 8-16, 32  | 4:2:2:4      | chroma shared between 2 pixels + alpha                        | planar      |
+----------------+-----------+--------------+---------------------------------------------------------------+-------------+
| YV12,YUV420Pxx | 8-16, 32  | 4:2:0        | chroma shared between 2x2 pixels                              | planar      |
+----------------+-----------+--------------+---------------------------------------------------------------+-------------+
| YUVA420Pxx     | 8-16, 32  | 4:2:0:4      | chroma shared between 2x2 pixels + alpha                      | planar      |
+----------------+-----------+--------------+---------------------------------------------------------------+-------------+
| YV411,YUV411P8 | 8         | 4:1:1        | chroma shared between 4 pixels                                | planar      |
+----------------+-----------+--------------+---------------------------------------------------------------+-------------+
| Y8,Y10-16,Y32  | 8-16, 32  | 4:0:0        | no chroma                                                     | both        |
+----------------+-----------+--------------+---------------------------------------------------------------+-------------+
| YUY2           | 8         | 4:2:2        | chroma shared between 2 pixels. Deprecated format, use YV16   | interleaved |
+----------------+-----------+--------------+---------------------------------------------------------------+-------------+
| xx refers to the bit depth, where                                                                                       |
| 8-16 = 8,10,12,14,16 bit integer; S or 32 is the 32 bit float. Use "S" in format names, except for Y, where Y32 is used.|
+----------------+-------------+-----------+--------------+---------------------------------------------------------------+

When the target format is the same as the source format, the original clip will be returned unchanged, 
except for the cases where the ``ChromaInPlacement`` and ``ChromaOutPlacement`` parameters are different,
or the target placement is different from the source chroma placement read from ``_ChromaLocation`` frame property.

Such functions are ``ConvertToYV12``/``ConvertToYUV420``/``ConvertToYUVA420`` or 
``ConvertToYV16``/``ConvertToYUV422``/``ConvertToYUVA422``.

Note ConvertToRGB always converts to RGB32 – unless your source is RGB24, in which case no conversion is done. 
If you need 24-bit RGB for some reason, use ConvertToRGB24 explicitly.

Syntax and operation of ``ConvertToRGB24`` is identical to ``ConvertToRGB``,
except that the output format is 24-bit; if the source is RGB32, the alpha
channel will be stripped.

Use ConvertBackToYUY2 to convert back to YUY2 with minimal color-blurring when you have previously applied a YUY2 -> RGB conversion. 

There is no unique way of converting YUV to RGB or vice-versa. There are different conversion matrices in use.
The following should be correct in most cases, see 
`here <http://avisynth.nl/index.php/Colorimetry#How_can_I_see_if_the_correct_standard_is_used_upon_playback>`_ for more.

- Rec.601 should be used when your source is standard definition (usually defined as smaller than 720p): 
  ``ConvertToRGB(clip) (using default "Rec601")``

- Rec.709 should be used when your source is DVD or HDTV: 
  ``ConvertToRGB(clip, matrix="Rec709")``

- The special-purpose matrices PC.601 and PC.709 keep the range unchanged (!), 
  instead of converting between 0d-255d RGB and 16d-235d YUV, as is the normal practice. 
  Note that if you want to convert from Y and no ``_ColorRange`` frame property present then it will treat it as
  limited range.

- The special-purpose matrix AVERAGE is used for (rarely found) sources with unweighted luma, 
  where ``Y = (R + G + B) / 3``. 

How conversion detects whether full or limited (narrow range for RGB) conversion needed:

- Avisynth can use the ``_ColorRange`` frame property to detect whether the source clip is of full or limited (narrow) range.
- If frame property ``_ColorRange`` is not present, then
  
  - ``full`` for RGB sources
  - ``limited`` for YUV or Y (greyscale) sources
  
  is assumed.
- Along with frame properties, the matrix string can contain additional "hints", such as ``:f``, ``:l``, , ``:auto``, ``:same``.
- When no other hint is given, some old-style Avisynth matrix name can specify limited/full: e.g. ``"Rec.709"`` implies limited range.
  Note: unlike ``PC.709`` or ``PC.601``: these matrix names do not force the clip being full or limited.

See also: :doc:`ConvertBits <convertbits>` to convert between bit depths and/or between full-limited range. 

Syntax and parameters
---------------------

.. _matrix_parameter_syntax:


.. describe:: matrix

    ``string  matrix = "Rec601"``

    Controls the colour coefficients and scaling factors used in colour space conversions.

    **Old-style constants (still valid):**

    These names remain fully supported for backward compatibility.

    - ``"Rec601"``               : Rec.601 coefficients; limited-range YUV ↔ full-range RGB.
    - ``"Rec709"``               : Rec.709 (HD) coefficients; limited-range YUV ↔ full-range RGB.
    - ``"Rec2020"``              : Rec.2020 (UHD) coefficients; limited-range YUV ↔ full-range RGB.
    - ``"PC.601"``, ``"PC601"``  : Rec.601 coefficients; range unchanged (passes through as-is).
    - ``"PC.709"``, ``"PC709"``  : Rec.709 (HD) coefficients; range unchanged.
    - ``"PC.2020"``, ``"PC2020"``: Rec.2020 (UHD) coefficients; range unchanged.
    - ``"Average"``              : Averaged luma (``Y = (R+G+B)/3``); range unchanged.

    Using old-style matrix names implies the full/limited range behaviour described above.
    ``PC.601``, ``PC.709`` and ``PC.2020`` do not force a range conversion: the input range
    is passed through to the output unchanged, equivalent to the new-style ``":same"``
    modifier (see below).

    **New-style syntax:**

    Two forms are available::

        "matrixname[:rangehint]"
        "matrixname:rangehint=>matrixname:rangehint"

    The first (single) form covers all common RGB ↔ YUV conversions.
    The second (two-part) form is required when independent control of both sides is needed,
    in particular for YUV → YUV conversions or when the RGB side needs an explicit range
    different from its default.

    *matrixname* can be one of:

    +----------------+-------+--------------------------------------------------------------+
    | Name           | Int   | Description                                                  |
    +================+=======+==============================================================+
    | ``rgb``        | 0     | Identity / RGB (valid on either side of ``=>``)              |
    +----------------+-------+--------------------------------------------------------------+
    | ``709``        | 1     | BT.709 (HD)                                                  |
    +----------------+-------+--------------------------------------------------------------+
    | ``unspec``     | 2     | Unspecified                                                  |
    +----------------+-------+--------------------------------------------------------------+
    | ``470bg``      | 5     | BT.470 BG (same coefficients as Rec.601)                     |
    +----------------+-------+--------------------------------------------------------------+
    | ``601``        | 5     | Alias for ``470bg``                                          |
    +----------------+-------+--------------------------------------------------------------+
    | ``fcc``        | 4     | BT.470 M / FCC                                               |
    +----------------+-------+--------------------------------------------------------------+
    | ``bt470m``     | 4     | Alias for ``fcc``                                            |
    +----------------+-------+--------------------------------------------------------------+
    | ``170m``       | 6     | SMPTE 170M                                                   |
    +----------------+-------+--------------------------------------------------------------+
    | ``240m``       | 7     | SMPTE 240M                                                   |
    +----------------+-------+--------------------------------------------------------------+
    | ``2020ncl``    | 9     | BT.2020 Non-Constant Luminance                               |
    +----------------+-------+--------------------------------------------------------------+
    | ``2020``       | 9     | Alias for ``2020ncl``                                        |
    +----------------+-------+--------------------------------------------------------------+
    | ``2020cl``     | 10    | BT.2020 Constant Luminance (treated as NCL internally)       |
    +----------------+-------+--------------------------------------------------------------+
    | ``ycgco``      | 8     | YCgCo *(not supported)*                                      |
    +----------------+-------+--------------------------------------------------------------+
    | ``chromancl``  | 12    | Chromaticity Derived NCL *(not supported)*                   |
    +----------------+-------+--------------------------------------------------------------+
    | ``chromacl``   | 13    | Chromaticity Derived CL *(not supported)*                    |
    +----------------+-------+--------------------------------------------------------------+
    | ``ictcp``      | 14    | ICtCp *(not supported)*                                      |
    +----------------+-------+--------------------------------------------------------------+

    *rangehint* specifies the colour range for that side of the conversion:

    - ``full`` or ``f``    — force full range
    - ``limited`` or ``l`` — force limited range
    - ``auto``             — use ``_ColorRange`` frame property if present, otherwise use
                             the format default (full for RGB, limited for YUV/Y)
    - ``same``             — **output side only** in the two-part form; copies the resolved
                             input range to the output. Invalid on the input (left) side.

    When no rangehint is given at all, ``auto`` behaviour applies.

    **Single form — the range hint refers to the YUV side:**

    In the single ``"matrixname:rangehint"`` form the range hint always describes the
    **YUV side** of the conversion, whichever side that is:

    - **YUV → RGB**: the hint describes the *input* (YUV) range.
      The RGB output defaults to full range independently of the hint.
    - **RGB → YUV**: the hint describes the *output* (YUV) range.
      The RGB input defaults to full range independently of the hint,
      unless overridden by a ``_ColorRange`` frame property on the source.

    This means ``"709:l"`` is unambiguous regardless of direction: the YUV side is
    limited, using BT.709 coefficients. The RGB side resolves its range independently.

    The special value ``same`` in the single form propagates from the *input* to the
    output, regardless of direction. For YUV → RGB, the YUV input range (from frame
    property or default) is copied to the RGB output. For RGB → YUV, the RGB input range
    is copied to the YUV output. In both cases the net effect is range pass-through,
    equivalent to the legacy ``PC.709`` / ``PC.601`` behaviour.

    .. note::

       For **RGB → YUV** with an explicit range hint, the hint sets the YUV *output*
       range only. The RGB input range is always resolved independently from the
       ``_ColorRange`` frame property (defaulting to full if absent) and cannot be
       overridden by the hint. Use the two-part form if you need to declare both sides
       explicitly.

    **Two-part form** ``"matrixname:rangehint=>matrixname:rangehint"``:

    Both sides must specify a matrixname and a rangehint. No implicit defaults apply
    across the ``=>``. The ``same`` modifier is only valid on the *right (output)* side;
    using it on the left side is an error.

    This form is used when:

    - A **YUV → YUV** conversion is needed. Internally the two matrices are combined and
      applied as a single pass in the 32-bit float unclipped domain, so any out-of-range
      RGB intermediate values that would arise from a naïve chained conversion are handled
      safely without clipping artefacts.

      YUV→RGB→YUV conversion using the specified matrix and range for each leg::

          ConvertToYUV444(clip, matrix="601:l=>709:l")
          # SD limited input, HD limited output

          ConvertToYUV444(clip, matrix="709:auto=>2020:same")
          # detect input range from frame prop; output matches resolved input range

          ConvertToYUV444(clip, matrix="709:f=>709:l")
          # full-range YUV input, limited-range YUV output, same matrix

    - An **RGB side range** needs to be declared explicitly, overriding the default.
      The matrixname ``"rgb"`` is valid on either side::

          ConvertToPlanarRGB(clip, matrix="709:auto=>rgb:limited")
          # detect YUV input range from frame prop; RGB output is forced limited

          ConvertToYUV444(clip, matrix="rgb:limited=>709:l")
          # explicitly declare limited RGB input; limited YUV output

    .. note::

       In the two-part form, ``same`` on the right side resolves relative to the *left
       side's resolved range*, not relative to the source clip's format default.
       For example, ``"709:f=>2020:same"`` produces full-range output because the left
       side resolved to full, regardless of what the default for the output format would
       otherwise be.

    **Equivalence table** — old-style names and their new-style equivalents:

    +--------------------+-------------------------------+------------------------------------------+
    | Old name           | New equivalent                | Notes                                    |
    +====================+===============================+==========================================+
    | ``Rec601``         | ``470bg:l``                   | limited YUV ↔ full RGB                   |
    +--------------------+-------------------------------+------------------------------------------+
    | ``Rec709``         | ``709:l``                     | limited YUV ↔ full RGB                   |
    +--------------------+-------------------------------+------------------------------------------+
    | ``Rec2020``        | ``2020ncl:l``                 | limited YUV ↔ full RGB                   |
    +--------------------+-------------------------------+------------------------------------------+
    | ``PC.601``         | ``470bg:same``                | range preserved                          |
    +--------------------+-------------------------------+------------------------------------------+
    | ``PC.709``         | ``709:same``                  | range preserved                          |
    +--------------------+-------------------------------+------------------------------------------+
    | ``PC.2020``        | ``2020ncl:same``              | range preserved                          |
    +--------------------+-------------------------------+------------------------------------------+
    | ``Average``        | ``average``                   | no standard ``_Matrix`` equivalent       |
    +--------------------+-------------------------------+------------------------------------------+

    **Range detection defaults summary:**

    +--------------------+----------------------------+----------------------------------------------------------+
    | Side               | Default range              | Notes                                                    |
    +====================+============================+==========================================================+
    | YUV / Y input      | limited                    | overridden by ``_ColorRange`` frame prop                 |
    +--------------------+----------------------------+----------------------------------------------------------+
    | YUV / Y output     | limited                    | overridden by rangehint or ``same`` or two-part form     |
    +--------------------+----------------------------+----------------------------------------------------------+
    | RGB input          | full                       | overridden by ``_ColorRange`` frame prop                 |
    +--------------------+----------------------------+----------------------------------------------------------+
    | RGB output         | full                       | overridden by ``same`` or two-part form                  |
    +--------------------+----------------------------+----------------------------------------------------------+



.. describe:: interlaced

    bool  interlaced = false 
    
    If true, it is assumed that clip is interlaced; by default, it is assumed to be progressive. 
    This option is needed because for example, the following (assuming clip is interlaced YV12):
    ::

        SeparateFields(clip)
        ConvertToYV16
        Weave


    ...is upsampled incorrectly. Instead it is better to use: 
    ::

        ConvertToYV16(clip, interlaced=true)


    Note, interlaced=true has an effect only on YV12 <-> YV16/YUY2 or YV12 <-> RGB conversions.
    (and their high bit depth equivalents).
    More about that can be found here: 
    :doc:`Color conversions and interlaced / field-based video <../advancedtopics/interlaced_fieldbased>`.

.. describe:: ChromaInPlacement, ChromaOutPlacement

    string  ChromaInPlacement = "MPEG2"
    
    string  ChromaOutPlacement = "MPEG2"

    ChromaInPlacement determines the chroma placement in the clip when converting from YV12/YUV420 or YV16/YUV422.
    ChromaOutPlacement determines the chroma placement in the clip when converting to YV12/YUV420 or YV16/YUV422.
    
    The placement can be one of these strings: 

    - ``"MPEG2"`` (synonyms: ``"left"``)
      Subsampling used in MPEG-2 4:2:x and most other formats. Chroma samples are located on the left pixel column of the group (default).
    - ``"MPEG1"`` (synonyms: ``"jpeg"``, ``"center"``)
      Subsampling used in MPEG-1 4:2:0. Chroma samples are located on the center of each group of 4 pixels.
    - ``"DV"``
      Like MPEG-2, but U and V channels are co-sited vertically: V on the top row, and U on the bottom row. For 4:1:1, chroma is located on the leftmost column.
    - ``"top_left"``
      Subsampling used in UHD 4:2:0. Chroma samples are located on the top left pixel column of the group.
    - ``bottom_left`` 4:2:0 only
    - ``bottom``   4:2:0 only 

   See also the Frame properties section below.


.. describe:: chromaresample

    string  chromaresample = "bicubic"

    Determines which chroma resampler is used in the conversion. Only used when the chroma resolutions 
    of the source and target are different. All AviSynth :doc:`resizers <resize>` are allowed 
    ("point", "bilinear", "bicubic", "lanczos", "lanczos4", "blackman", "spline16", "spline36", "spline64", 
    "gauss" and "sinc", "sinpow", "sinclin2" and "userdefined2"). 
    
    Default is "bicubic". 

.. describe:: param1, param2, param3

    These 'float' type parameters can be the additional parameters for the chroma
    resamplers. Some resizer algorithms would need and can be fine tuned with up to 3 parameters.
    Their default values depend on the selected chromaresample resizer kernel,

.. describe:: bits

    Used by ConvertToPlanarRGB(A) to perform on-the-fly output bit-depth conversion.

    **Internal calculation methods of 8-16 bit sources** (when conversion is needed)

    ========================  ===================  ===============  ================================
    Target Range              Internal Math        Internal Math    Output Handling
                                                   (quality=true)
    ========================  ===================  ===============  ================================
    Full-range                32-bit float         32-bit float     Direct output
    Limited-range → integer   S18.13 fixed-point   32-bit float     Truncated to target bit depth
    Limited-range → float     S18.13 fixed-point   32-bit float     Converted to float (no truncation)
    ========================  ===================  ===============  ================================

    When ``quality=true`` (see below), the S18.13 fixed-point path is replaced by
    32-bit float processing regardless of target range or bit depth.

    Note: Limited-range to float conversion preserves the full precision of the
    S18.13 fixed-point calculation by converting directly to 32-bit float without
    the truncation that occurs with integer targets.

.. describe:: quality

    ``bool  quality = false``

    Available in ConvertToPlanarRGB(A) only.

    When ``false`` (default), the internal calculation method is chosen automatically
    based on the target range and bit depth, as described in the ``bits`` table above:
    full-range targets use 32-bit float, limited-range targets use the faster S18.13
    scaled integer path.

    When ``true``, all internal processing is forced to 32-bit float regardless of
    target range or bit depth. This avoids the rounding that occurs when the S18.13
    fixed-point path truncates to an integer target, at the cost of slightly more
    computation. Recommended when the output will undergo further processing and
    maximum precision is desired.


Frame properties
----------------

Since Avisynth v3.7.1 frame property (_ChromaLocation) support appears in selected filters (e.g. ConvertToYUV422). 
Property can be read and/or set. A frame property can replace default behaviour of location parameters and is set 
(or deleted) upon finishing conversion. Since a format without subsampling - such as 4:4:4 (YV24) - does not have 
chroma location, the property is deleted automatically when converting to YUV444 or RGB.

- "ChromaInPlacement" rules:

    * if source has _ChromaLocation frame property it will be used else the default is "mpeg2" ("left")
    * if parameter is "auto" or not given at all, ChromaInLocation will be set to the above mentioned default value
    * if parameter is explicitely given, it will be used 

- "ChromaOutPlacement" rules:

    * default is "mpeg2" ("left")
    * if parameter is "auto" or not given at all, ChromaOutLocation will be set to the above mentioned default value
    * if parameter is explicitely given, it will be used 

    Accepted values for "ChromaInPlacement" and "ChromaOutPlacement" (when source/target is a chroma subsampled format) 
    (full list):

    * "left" or "mpeg2"
    * "center" or "jpeg" or "mpeg1"
    * "top_left"
    * "dv"
    * "top"
    * "bottom_left"
    * "bottom"

- _ChromaLocation constants - just for info: as seen in propShow

    * AVS_CHROMA_LEFT = 0
    * AVS_CHROMA_CENTER = 1
    * AVS_CHROMA_TOP_LEFT = 2 (4:2:0 only)
    * AVS_CHROMA_TOP = 3 (4:2:0 only)
    * AVS_CHROMA_BOTTOM_LEFT = 4 (4:2:0 only)
    * AVS_CHROMA_BOTTOM = 5 (4:2:0 only)
    * AVS_CHROMA_DV = 6 Special to Avisynth 


Conversion paths
----------------

-   The *ChromaInPlacement*, *chromaresample* and *ChromaOutPlacement*
    options are only used in the 'planar conversion part' of the conversion
    path, and they process the chroma of the clip.

The following conversion paths occur

-   YUV planar -> RGB via YV24
-   YUV planar -> YUY2 via YV16
-   RGB -> YUV planar via YV24
-   YUY2 -> YUV planar via YV16

Suppose you have a YUY2 clip for example and you convert it to YV24.

The YUY2 will be converted to YV16 first without applying *ChromaInPlacement*,
*chromaresample* and *ChromaOutPlacement*. 

Then YV16 will be converted to YV24 while applying *chromaresample*. 
*ChromaOutPlacement* won't be used since our target is YV24.


Sampling
--------

:doc:`This part of the documentation <../advancedtopics/sampling>` covers the sampling methods and color formats in more detail.


Color conversions
-----------------

:doc:`This page <../advancedtopics/color_conversions>` covers the color conversions, "YUV <-> RGB", in more detail.

+----------+------------------------------------------------------------+
| Changes: |                                                            |
+==========+============================================================+
| v3.7.6   || Add "quality" parameter to ConvertToPlanarRGB(A)          |
|          || Add "bits" parameter to ConvertToPlanarRGB(A)             |
|          || Document ":same" in matrix specifier                      |
+----------+------------------------------------------------------------+
| v3.7.3   || Added "sinpow",  "sinclin2" and "userdefined2" to         |
|          |  chromaresampler options                                   |
|          || Add "param1", "param2" and "param3" to ConvertToXXXX where|
|          |  'chromaresample' exists (b,c,s,taps and p parameters can  |
|          |  be set, depending on the resizer.)                        |
|          || Add ConvertToYUVA420, ConvertToYUVA422, ConvertToYUVA444  |
+----------+------------------------------------------------------------+
| v3.7.1   || Added ChromaOutPlacement to 4:2:2 related functions       |
|          || Added new matrix constants, optional new syntax           |
|          || Added new chroma location constants                       |
|          || Added _ChromaLocation frame property                      |
+----------+------------------------------------------------------------+
| v2.60    || Added: ConvertToY8, ConvertToYV411,                       |
|          |  ConvertToYV16, ConvertToYV24,                             |
|          || Added ChromaInPlacement, ChromaOutPlacement and           |
|          |  chromaresample, matrix="AVERAGE"                          |
+----------+------------------------------------------------------------+
| v2.50    | ConvertToYV12                                              |
+----------+------------------------------------------------------------+

$Date: 2026/02/24 20:25:00 $
