
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
         int bits] )
    ConvertToPlanarRGBA(clip, [ string matrix, bool interlaced,
         string ChromaInPlacement,
         string chromaresample,
         float param1, float param2, float param3,
         int bits] )


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

.. describe:: matrix

    string  matrix = "Rec601"

    Default "Rec601". Controls the colour coefficients and scaling factors used in RGB - YUV conversions.

    Old-style constants:

    - "Rec601"  : Uses Rec.601 coefficients; scale full range [0d..255d] RGB ↔ TV range [16d..235d] YUV.
    - "Rec709"  : Uses Rec.709 (HD) coefficients; scale full range RGB ↔ TV range YUV.
    - "Rec2020" : Uses Rec.2020 (UHD) coefficients; scale full range RGB ↔ TV range YUV.
    - "PC.2020" : Uses Rec.2020 (UHD) coefficients; keep range unchanged.
    - "PC.601"  : Uses Rec.601 coefficients; keep range unchanged.
    - "PC.709"  : Uses Rec.709 (HD) coefficients; keep range unchanged.
    - "Average"  : Uses averaged coefficients (the luma becomes the average of the RGB channels); keep range unchanged. 

    Additional matrix constants:

    New syntax: more matrix string constants with separate full/limited range markers.
    ``"matrixname:full_or_limited_or_auto_or_same"`` where 
    ``"matrixname"`` can be set as (for developers, internal _Matrix integer constant are given in parenthesys)

    - "rgb" (0 - AVS_MATRIX_RGB)
    - "709" (1 - AVS_MATRIX_BT709)
    - "unspec" (2 - AVS_MATRIX_UNSPECIFIED)
    - "170m" (6 - AVS_MATRIX_ST170_M)
    - "240m" (7 - AVS_MATRIX_ST240_M)
    - "470bg" (5 - AVS_MATRIX_BT470_BG)
    - "601" (5 - AVS_MATRIX_BT470_BG)
    - "fcc" (4 - AVS_MATRIX_BT470_M)
    - "bt470m" (4 - AVS_MATRIX_BT470_M)
    - "ycgco" (8 - AVS_MATRIX_YCGCO not supported)
    - "2020ncl" (9 - AVS_MATRIX_BT2020_NCL)
    - "2020" (9 - AVS_MATRIX_BT2020_NCL)
    - "2020cl" (10 - AVS_MATRIX_BT2020_CL same as 2020ncl)
    - "chromacl" (13 - AVS_MATRIX_CHROMATICITY_DERIVED_CL not supported)
    - "chromancl" (12 - AVS_MATRIX_CHROMATICITY_DERIVED_NCL not supported)
    - "ictcp" (14 - AVS_MATRIX_ICTCP not supported) 

    The above "matrix" parameters can be followed by a 
    
    - ``"full"`` or ``"f"`` 
    - ``"limited"`` or ``"l"``
    - ``"auto"``
    - ``"same"``
    
    marker after a ``":"``.

    e.g. ``"709:l"`` means the same as the old "Rec709", since it defaults to limited to full conversion.

    When there is no limited-ness marker, or is set to "auto" then value of _ColorRange frame property is used.
    Using "same" will assume the input range for the output's range.

    Note: Avisynth+ defines a new matrix syntax, but old-style "matrix" parameter names are still valid.
    Using old-style matrix names imply the full/limited range, except ``"PC.601"`` and ``"PC.709"`` which 
    do not alter the input's range.

    For memo and the similar new string

    - "rec601" same as "470bg:l" (limited to full)
    - "rec709" "709:l" (limited to full)
    - "pc.601" and "pc601" same as "470bg:f" - but only if source has _ColorRange = 0 (full) 
    - "pc.709" and "pc709" same as "709:f" - but only if source has _ColorRange = 0 (full)
    - "pc.601" and "pc601" same as "470bg:same" (keep input range)
    - "pc.709" and "pc709" same as "709:same" (keep input range)
    - "average" - kept for compatibility, really it has no standard _Matrix equivalent
    - "rec2020" "2020cl:l"
    - "pc.2020" and "pc2020" "2020cl:f" - but only if source has _ColorRange = 0 (full)

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

    **Internal calculation methods:** (when conversion is needed)
        
        ========================  ===================  ================================
        Target Range              Internal Math        Output Handling
        ========================  ===================  ================================
        Full-range                32-bit float         Direct output
        Limited-range → integer   S18.13 fixed-point   Truncated to target bit depth
        Limited-range → float     S18.13 fixed-point   Converted to float (no truncation)
        ========================  ===================  ================================
        
        Note: Limited-range to float conversion preserves the full precision of the 
        S18.13 fixed-point calculation by converting directly to 32-bit float without 
        the truncation that occurs with integer targets.
    
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
| v3.7.6   || Add "bits" parameter to ConvertToPlanarRGB()              |
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

$Date: 2026/02/12 10:44:00 $
