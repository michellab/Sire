== Instructions for compiling and installing on windows ==

== MSYS2 BUILD ==

First, install Visual Studio 2015 (Community Edition is free for open source)

First, install msys2 by following the instructions on http://msys2.github.io
(follow the 64bit instructions as we should target 64bit windows)

Then ensure that you have fully updated, e.g. via

  # pacman -Syuu

Then install

  # pacman -S git
  # pacman -S cmake
  # pacman -S make

You need to ensure that you are using the mingw build tools. These are installed
into /mingw64, so you need to add /mingw64/bin to you PATH. We also need to 
install the complete mingw64 build toolchain (compilers, make, linkers etc.)

Do this by typing;

  # pacman -S mingw-w64-x86_64-gcc

Also install nano if you like this editor ;-)

  # pacman -S nano

Now install the dependencies of Sire - I am giving up getting these compiled
manually!

  # pacman -Sy mingw-w64-x86_64-intel-tbb
  # pacman -Sy mingw-w64-x86_64-qt5
  # pacman -Sy mingw-w64-x86_64-gsl
  # pacman -Sy mingw-w64-x86_64-boost
  # pacman -Sy mingw-w64-x86_64-python3
  # pacman -Sy mingw-w64-x86_64-pkg-config

If this has worked, then typing "which gcc" and "which g++" should return

  # /mingw64/bin/gcc
  # /mingw64/bin/g++

Then clone the Sire repository using

  # git clone https://github.com/michellab/Sire.git

If this fails with a "child_info_fork" error, then you need to fix
your msys2 installation. Do this by exiting from msys2 and then running
autorebase.bat which is in the C:\msys2 directory. Then go back into
msys2 and try again.

Run cmake using

cmake -G "MSYS Makefiles"

Then build using mingw32-make

This works, and builds corelib

For the wrappers, need

cmake -G "MSYS Makefiles" -DCMAKE_NEED_RESPONSE=1

Then, need to set the windows PATH. Can do this in a standard CMD shell, e.g.

PATH=C:\msys64\home\chzcjw\sire.app\bin;C:\msys64\home\chzcjw\sire.app\lib;C:\msys64\mingw64\lib;\msys64\mingw64\lib;%PATH%
set PYTHONHOME=C:\msys64\mingw64
set PYTHONPATH=C:\msys64\mingw64\lib\python3.5;C:\msys64\home\chzcjw\sire.app\lib\python\site-packages


