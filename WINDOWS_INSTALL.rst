== Instructions for compiling and installing on windows ==

== All of this is EXPERIMENTAL == 
This will be updated and replaced by instructions to compile against the
acaconda that is available as part of the WSL (Windows subsystem for Linux)
that is available as part of Windows 10
 
== MSYS2 BUILD ==

First, install msys2 by following the instructions on http://msys2.github.io
(follow the 64bit instructions as we should target 64bit windows)

Then ensure that you have fully updated, e.g. via

  # pacman -Syuu

Then install

  # pacman -S git
  # pacman -S mingw-w64-x86_64-cmake
  # pacman -S mingw-w64-x86_64-make
  # pacman -S mingw-w64-x86_64-toolchain

You need to ensure that you are using the mingw build tools. These are installed
into /mingw64, so you need to add /mingw64/bin to you PATH. We also need to 
install the complete mingw64 build toolchain (compilers, make, linkers etc.)

Do this by typing;

  # pacman -S mingw-w64-x86_64-gcc

Also install nano if you like this editor ;-)

  # pacman -S nano

Now install the dependencies of Sire - I gave up getting these compiled
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

Next switch to the c++15 branch that holds the development code needed
to support windows

  # git checkout c++15

Change into the corelib build directory using

   # mkdir build/corelib
   # cd build/corelib

Now we need to link mingw32-make to make

   # ln -s /mingw64/bin/mingw32-make.exe /mingw64/bin/make

Run cmake using

   # cmake -G "MSYS Makefiles" ../../corelib

Then build using 

   # mingw32-make install

 (note that you may have to copy sire.app/bundled/bin/libcpuid-10.dll to sire.app/bundled/lib/libcpuid.dll)


This works, and builds corelib

For the wrappers, change into the wrapper build directory

   # cd ../wrapper

Then run cmake again, using

   # cmake -G "MSYS Makefiles" -DCMAKE_NEED_RESPONSE=1 ../../wrapper

Then build using

   # mingw32-make install

While Sire compiles within the MSYS shell, I haven't got it to run from there
yet. You need to run from the standard windows cmd shell.

First, you need to set the windows PATH. Can do this in a standard CMD shell, e.g.

   # PATH=C:\msys64\home\chzcjw\sire.app\bin;C:\msys64\home\chzcjw\sire.app\lib;C:\msys64\mingw64\lib;C:\msys64\mingw64\bin;%PATH%
   # set PYTHONPATH=C:\msys64\mingw64\lib\python3.5;C:\msys64\home\chzcjw\sire.app\lib\python\site-packages

(note that you will need to update the above based on where you compiled Sire, and
 where msys64 is located on your machine)

Now you can run Sire. Try

    # sire_python

This should start up the Sire python shell. Into this, type

    # import Sire.Mol

If this runs, then Sire has successfully imported the Sire.Mol module,
meaning that all of the libraries have been found correctly.

To run waterswap or ligandswap, you need to run the associated python script
using the sire_python executable. For example

    # sire_python C:\msys64\home\chzcjw\sire.app\share\Sire\scripts\waterswap.py --help

** At the moment, this is failing because my MSYS python can't find the _struct
   module. This will need to be fixed... **

** FIX IS TO COPY SIRE_PYTHON AND ALL LIBRARIES INTO MSYS/BIN. NEED TO 
   INSTALL SIRE WITHIN MSYS!!! **

## Running tests - still not fully working

Install nose. (quite difficult as python easy_install breaks with error
VC 6.0 is not supported. I am afraid that I didn't record how I got nose installed!)

Then run within the 

cd \path\to\test_directory
\path\to\sire.app\bin\sire_python -m nose
