
Using AviSynth+ on POSIX systems
================================

As of version 3.5, AviSynth+ can now be built and used natively
on Linux, macOS, and BSD.

.. toctree::
    :maxdepth: 3

.. contents:: Table of contents



AviSynth+ prerequisites
-----------------------

Depending on your OS or distribution, the commands to fetch
the necessary prerequisites for building AviSynth+ differ.

At a bare minimum:

* CMake 3.8 or higher.
* GCC 8 or higher, or similarly recent version of Clang or AppleClang.

.. note::
   The use of Ninja as the generator for CMake is a matter of personal preference.
   Feel free to use GNU Make if so compelled (i.e. just a plain 'cmake ..' invocation).

Linux
^^^^^

Ubuntu 19.10 or higher
~~~~~~~~~~~~~~~~~~~~~~

::

    sudo apt-get install build-essential cmake git ninja-build checkinstall


::

    git clone https://github.com/AviSynth/AviSynthPlus && \
    cd AviSynthPlus && \
    mkdir avisynth-build && \
    cd avisynth-build && \

    cmake ../ -G Ninja -DCMAKE_BUILD_TYPE=Release && \
    ninja && \
        sudo checkinstall --pkgname=avisynth --pkgversion="$(grep -r \
        Version avs_core/avisynth.pc | cut -f2 -d " ")-$(date --rfc-3339=date | \
        sed 's/-//g')-git" --backup=no --deldoc=yes --delspec=yes --deldesc=yes \
        --strip=yes --stripso=yes --addso=yes --fstrans=no --default ninja install


Ubuntu 18.04 LTS
~~~~~~~~~~~~~~~~

18.04 ships with GCC 7, which is not sufficient to build AviSynth+ without
the use of the `filesystem submodule`_.

::

    git clone --recursive https://github.com/AviSynth/AviSynthPlus && \
    cd AviSynthPlus && \
    mkdir avisynth-build && \
    cd avisynth-build && \

    cmake ../ -G Ninja -DCMAKE_BUILD_TYPE=Release && \
    ninja && \
        sudo checkinstall --pkgname=avisynth --pkgversion="$(grep -r \
        Version avs_core/avisynth.pc | cut -f2 -d " ")-$(date --rfc-3339=date | \
        sed 's/-//g')-git" --backup=no --deldoc=yes --delspec=yes --deldesc=yes \
        --strip=yes --stripso=yes --addso=yes --fstrans=no --default ninja install

Raspbian Raspberry Pi 5 + llvm + ninja
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Raspberry Pi 5 is an aarch64 architecture, presently (2025) comes with a gcc 12.2.
Unfortunately it is not able to optimize well, namely it does not recognize well
the vectorizable code, and even produces slower code with a so-called "vector attribute"
that without it. Probably the ArmV8-a support is not perfect.

Using llvm instead is the way to go. At the moment it comes with 14.0.6 version for this
distribution.

I'm still using here the -DCMAKE_BUILD_TYPE=Release flag, but after a May 2025 source
it is no longer needed to avoid a completely unoptimized build with Ninja.

First, get some basic stuff.

::

    sudo apt-get install build-essential cmake git ninja-build checkinstall

::

    git clone https://github.com/AviSynth/AviSynthPlus && \
    cd AviSynthPlus && \
    mkdir avisynth-build && \
    cd avisynth-build && \

From now on, instead of a copy-pastable content, we create script files with the following contents.

*config.sh*: grabs the missing components for compiling AviSynth with llvm.
The prefix path is very important.
::

    #!/bin/bash

    # Configuration variables for LLVM build
    LLVM_PACKAGE="llvm"
    LLVM_CMAKE_PREFIX_PATH="/usr/lib/llvm-$(ls /usr/lib/llvm-* 2>/dev/null | sed 's>

    echo "Using LLVM CMake prefix path: $LLVM_CMAKE_PREFIX_PATH"

    # Optional: Install LLVM, libc++-dev, and libc++abi-dev if not already present
    if ! command -v clang++ &> /dev/null || ! dpkg -s libc++-dev &> /dev/null || ! >
      echo "LLVM, libc++-dev, or libc++abi-dev not found. Installing..."
      sudo apt update
      sudo apt install -y "$LLVM_PACKAGE" libc++-dev libc++abi-dev
    fi


*build-llvm.sh* makes an actual clean build, assumes that avisynth was cloned
in the right folder in the previous step.

You can remove the general purge by remove the ``rm -rf *`` if everythings behaves as it should.

::

    #!/bin/bash

    # Source the configuration script
    source ./config.sh

    cd ~/AviSynthPlus/avisynth-build
    rm -rf *

    cmake ../ -G Ninja \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_C_COMPILER=clang \
        -DCMAKE_CXX_COMPILER=clang++ \
        -DCMAKE_PREFIX_PATH="$LLVM_CMAKE_PREFIX_PATH"

    ninja

    sudo checkinstall --pkgname=avisynth --pkgversion="$(grep -r \
        Version avs_core/avisynth.pc | cut -f2 -d " ")-$(date --rfc-3339=date | \
        sed 's/-//g')-git" --backup=no --deldoc=yes --delspec=yes --deldesc=yes \
        --strip=yes --stripso=yes --addso=yes --fstrans=no --default ninja install

    cd ~bin

    echo "AviSynth+ built and installed with LLVM!"


When you are happy with gcc, use this script:

build-gcc.sh:
::

    #!/bin/bash

    cd ~/AviSynthPlus/avisynth-build
    rm -rf *

    cmake ../ -G Ninja -DCMAKE_BUILD_TYPE=Release && \
    ninja && \
        sudo checkinstall --pkgname=avisynth --pkgversion="$(grep -r \
        Version avs_core/avisynth.pc | cut -f2 -d " ")-$(date --rfc-3339=date | \
        sed 's/-//g')-git" --backup=no --deldoc=yes --delspec=yes --deldesc=yes \
        --strip=yes --stripso=yes --addso=yes --fstrans=no --default ninja install

    cd ~bin





Distributions without checkinstall
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Not all Linux distributions have checkinstall in their repositories, either due to
a lack of checkinstall working with their package management system or simply due
to omission.  In these cases, the install step is a little different:

::

    sudo ninja install
    sudo ldconfig


macOS
^^^^^

| Requires Homebrew:
| `<https://brew.sh/>`_

::

    brew install cmake ninja gcc


GCC isn't strictly necessary for AviSynth+, but it can side-step
the need to use `an external implementation`_ on High Sierra and
Mojave.


10.13 High Sierra and 10.14 Mojave
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Apple's libc++ doesn't support the C++17 filesystem functionality
on either of these versions of macOS, so we have to resort to
using `an external implementation`_ as a submodule.

::

    git clone --recursive https://github.com/AviSynth/AviSynthPlus && \
    cd AviSynthPlus && \
    mkdir avisynth-build && \
    cd avisynth-build

    cmake ../ -G Ninja -DCMAKE_BUILD_TYPE=Release && \
    ninja && \
    sudo ninja install


10.15 Catalina and higher
~~~~~~~~~~~~~~~~~~~~~~~~~

C++17 filesystem support is available on Catalina, so it can
be built with the default Clang installation.

::

    cmake ../ -G Ninja -DCMAKE_BUILD_TYPE=Release && \
    ninja && \
    sudo ninja install


FreeBSD
^^^^^^^

Tested on FreeBSD 12.1.

::

    pkg install cmake git gmake ninja

    git clone https://github.com/AviSynth/AviSynthPlus && \
    cd AviSynthPlus && \
    mkdir avisynth-build && \
    cd avisynth-build


Building AviSynth+ (GNU Make)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    cmake ../ && \
    gmake -j$(nproc) && \
    gmake install


Building AviSynth+ (Ninja)
~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    cmake ../ -G Ninja && \
    ninja && \
    sudo ninja install


FFmpeg support
--------------

On all of these OSes, AviSynth+ can interface with FFmpeg.
This change was applied to the FFmpeg git master branch on
2020-04-05.

It is always useful to check their guide before the process:
https://trac.ffmpeg.org/wiki/CompilationGuide


To compile a basic build of FFmpeg that supports
AviSynth+, the following steps will suffice:

Prerequisites
^^^^^^^^^^^^^

Linux
~~~~~

Ubuntu
......

First, enable the Sources repository by either enabling it
using the Software Sources dialog or by uncommenting the
right lines in /etc/apt/sources.list.

::

    sudo apt-get build-dep ffmpeg
    sudo apt-get install nasm libsdl2-dev


macOS
~~~~~

Homebrew prerequisites:

::

    brew install xz sdl2 pkg-config nasm


FreeBSD
~~~~~~~

::

    pkg install nasm sdl2


Building FFmpeg
^^^^^^^^^^^^^^^

::

    git clone https://git.videolan.org/git/ffmpeg.git
    cd ffmpeg


Linux
~~~~~

Ubuntu
......

::

    ./configure --prefix=$HOME/ffmpeg_build --enable-gpl --enable-version3 \
    --disable-doc --disable-debug --enable-pic --enable-avisynth && \
    make -j$(nproc) && \
    make install


Installing FFmpeg to the system can be done by leaving out the `--prefix`
option and then using the following checkinstall command:

::

    sudo checkinstall --pkgname=ffmpeg --pkgversion="7:$(git rev-list \
    --count HEAD)-g$(git rev-parse --short HEAD)" --backup=no --deldoc=yes \
    --delspec=yes --deldesc=yes --strip=yes --stripso=yes --addso=yes \
    --fstrans=no --default


Raspbian
........

Raspberry Pi 5 is able to compile ffmpeg and its prerequisites within a reasonable time, 
which means 7-10 minutes (when you don't forget setting -j4 for multithreaded build process).

This is a general Ubuntu/Debian compilation guide by the ffmpeg project:
https://trac.ffmpeg.org/wiki/CompilationGuide/Ubuntu

You can also refer to this video (Raspberry Pi 5 - Compile FFMPEG on Raspberry PI OS):
https://www.youtube.com/live/if9UG4kJ9L4?si=Nq802tiUteRdlyHE

which helps you through the process.

Important: you'll need some additions, like extending the configuration with 
``--enable-avisynth`` and maybe ``--disable-doc``.



macOS
~~~~~

::

    ./configure --prefix=$HOME/ffmpeg_build --enable-gpl --enable-version3 --disable-doc \
    --disable-debug --enable-avisynth
    make -j$(nproc)
    make install

On Catalina, `--extra-cflags="-fno-stack-check"` is necessary when using AppleClang as the compiler.

FreeBSD
~~~~~~~

::

    ./configure --prefix=$HOME/ffmpeg_build --enable-gpl --enable-version3 --disable-doc \
    --disable-debug --enable-pic --enable-avisynth --cc=cc
    gmake -j$(nproc)
    gmake install


Testing the installation
------------------------

FFplay can be used to preview scripts in a pinch; if mpv or VLC is built against the patched
version of FFmpeg, those can be used to play back scripts in a more comfortable player
experience.

The easiest two scripts to test the installation are Version or Colorbars/ColorbarsHD.

::

    Version()


::

    Colorbars() # or ColorbarsHD()


And running this script in the test build of FFmpeg:

::

    cd ~/ffmpeg_build/bin


Create the script in this directory, for ease of testing.

To play the script::

    ./ffplay -i test.avs

To convert as usual::

    ./ffmpeg -i test.avs [encoding options]

Benchmark::

    ./ffmpeg -benchmark -i test.avs -f null -

Stream to a GUI, e.g. to VLC Media Player on a Raspberry Pi.
Choose "Open Network Stream" from menu, and enter ``udp://@127.0.0.1:10000`` 
for the network url.
::

    ./ffmpeg -i test.avs -pix_fmt yuv420p -f mpegts udp://127.0.0.1:10000


Troubleshooting:

- ffmpeg reports an unknown error. One possible reason is that libavisynth.so could
  not be loaded; the ``$LD_LIBRARY_PATH`` environment variable may be empty.
  Check with ``echo $LD_LIBRARY_PATH``.
  Fix it by (Raspberry Pi 5 Raspbian):

::

    echo 'export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH' >> ~/.bashrc

- ffmpeg does not recognize the AviSynth script. You probably typed ``ffmpeg`` 
  instead of ``./ffmpeg``. For example, Raspberry Pi 5 Raspbian comes with a 5.x ffmpeg 
  without AviSynth support, and your command finds that one.


Loading actual video sources will require a source filter.  FFMS2 doesn't require any porting
to these OSes, making it the most straightforward option at the moment.


Building FFMS2
--------------

FFMS2 doesn't require any additional prerequisites, so it can be
built straight away.

::

    git clone https://github.com/ffms/ffms2 && \
    cd ffms2


Linux
^^^^^

Ubuntu
~~~~~~

::

        PKG_CONFIG_PATH=$HOME/ffmpeg_build/lib/pkgconfig \
        CPPFLAGS="-I/usr/local/include/avisynth" \
        ./autogen.sh --enable-shared --enable-avisynth && \
    make -j$(nproc) && \
        sudo checkinstall --pkgname=ffms2 --pkgversion="1:$(./version.sh)-git" \
        --backup=no --deldoc=yes --delspec=yes --deldesc=yes --strip=yes --stripso=yes \
        --addso=yes --fstrans=no --default


macOS
^^^^^

::

    brew install autoconf automake libtool m4

        PKG_CONFIG_PATH=$HOME/ffmpeg_build/lib/pkgconfig \
        CPPFLAGS="-I/usr/local/include/avisynth" \
        ./autogen.sh --enable-shared --enable-avisynth && \
    make -j$(nproc) && \
    sudo make install


FreeBSD
^^^^^^^

::

    pkg install autoconf automake libtool m4

        PKG_CONFIG_PATH=$HOME/ffmpeg_build/lib/pkgconfig \
        CPPFLAGS="-I/usr/local/include/avisynth" \
        ./autogen.sh --enable-shared --enable-avisynth && \
    gmake -j$(nproc) && \
    gmake install


Plugin autoloading
------------------

AviSynth+ will use several directories for autoloading:
the `avisynth/` subdirectory where libavisynth.so was installed,
`$HOME/.avisynth`, and the directory given to the USER_AVS_PLUGINDIR_LOCATION
configuration option (defaults to `$HOME/.local/lib/avisynth`).
The latter of which can hold plugins (and symlinks to plugins)
or AVSI files without needing root permissions.

On FreeBSD, procfs needs to be mounted first in order for
autoloading to function.


Back to the :doc:`main page <../../index>`

$ Date: 2025-04-15 15:15:00 $

.. _an external implementation: https://github.com/gulrak/filesystem
.. _filesystem submodule: https://github.com/gulrak/filesystem
