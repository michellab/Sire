== Instructions for compiling and installing on windows ==

We have updated our builds so that it is now possible to compile and link Sire on Windows.

=== Setting up the environment ===

To do this, you first need to install Visual Studio 2017 Build Tools (free as the community edition)

Download this from https://visualstudio.microsoft.com/vs/older-downloads

From here, download **Build Tools for Visual Studio 2017**. Run the installer, and 
select only **Visual Studio C++ build tools** (the top box in the left panel). Then
click the "Install" button on the bottom right. This will now download and install 
all of the Visual Studio build tools needed to compile Python modules that are compatible
from Python 3.7+

Next, download and install a conda environment. We recommend mambaforge, downloaded 
from here: https://github.com/conda-forge/miniforge

We use `Mambaforge-Windows-x86_64 <https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Windows-x86_64.exe>`__.

Install this using "Just for you" into a directory. We use `C:/Users/{username}/sire_app` (replace
`{username}` with your Windows user name)

Once it has installed, go to the start menu and look for **Visual Studio 2017 / Developer Command Prompt for VS2017**.
This will launch a command prompt.

In this prompt, change directory into your conda install, e.g.

.. code-block::

   cd C:\Users\{username}\sire_app


Now tell the CMD prompt about conda.

.. code-block::

   conda init cmd.exe


You only need to do this once, after you have installed your conda.

Now you can activate conda by typing 

.. code-block::

   conda activate


You should see your prompt change to include `(base)` at the start. Type 

.. code-block::

   python -V


and you should see that the version of Python you downloaded has been run, e.g. I see

.. code-block::

   Python 3.9.7


Next, you need to install git (if you don't have it already). You can do this either by installing git 
by `following the instructions here <https://git-scm.com/download/win>`__ or by installing through
conda via `conda install git`.

Now you can change back to your home directory and download Sire.

.. code-block::

   cd C:\Users\{username}

   git clone https://github.com/michellab/Sire

=== Compiling Sire ===

Change into the Sire directory.

.. code-block::

   cd Sire

Now run the `compile_sire.bat` script by typing

.. code-block::

   compile_sire

Sire will now compile and install itself into your conda directory. This will take a VERY long time (hours).

If everything works, you will see printed out at the end;

.. code-block::

   XXX

== Running Sire ==

You should now (hopefully) be able to import Sire into your script. To test, start 
`ipython` and then type in

>>> import sire as sr

If this imports ok then Sire is working. You are then able to go on to tutorial.

== Diagnosing problems ==

There are lots of things that can go wrong on Windows. If you see "_Qt.pyd cannot be found"
or similar errors then this is likely because Windows can't find the required libraries.

To debug, list the dependent libraries of `_Qt.pyd` via

.. code-block::

   dumpbin /dependents C:\Users\{username}\sire_app\lib\site-packages\sire\legacy\Qt\_Qt.pyd

You should see a list of libraries, e.g.

.. code-block::

   Microsoft (R) COFF/PE Dumper Version 14.16.27048.0
   Copyright (C) Microsoft Corporation.  All rights reserved.

   Dump of file C:\Users\{username}\sire_app\lib\site-packages\sire\legacy\Qt\_Qt.pyd

   File Type: DLL

   Image has the following dependencies:

    SirePython.dll
    boost_python39.dll
    python39.dll
    SireStream.dll
    Qt5Core_conda.dll
    VCRUNTIME140.dll
    api-ms-win-crt-heap-l1-1-0.dll
    api-ms-win-crt-runtime-l1-1-0.dll
    KERNEL32.dll

   Summary

       12000 .data
        5000 .pdata
       22000 .rdata
        1000 .reloc
        1000 .rsrc
       4C000 .text

If you see 

.. code-block::

   LINK : fatal error LNK1181: cannot open input file 'C:\Users\{username}\sire_app\lib\site-packages\sire\legacy\Qt\_Qt.pyd'

Then `_Qt.pyd` has not been installed by the Sire installation, which implies something big has gone wrong.
Please raise an issue about this on our repo.

Assuming you do see the shared libraries, then you should see the same libraries as we show above, e.g.
`SirePython.dll`, `boost_python39.dll` and `python39.dll` at a minimum. This shows that both Python and 
boost python have been dynamically linked to the `_Qy.pyd` module. This is very important, as Sire will not
work if boost python or Python are statically linked.

Next, check to see if the libaries are in their expected location. The Sire libraries (e.g. `SireStream.dll` and `SirePython.dll`)
should be in `C:\Users\{username}\sire_app\Library\bin`, along with the other DLLs that are imported by packages in conda.

If the libraries aren't there, then something has gone wrong with the install. Please get in touch with us and we can try to debug.

You can look to see what is missing (e.g. a missing depedent library) by downloading (the free) Dependency Walker and 
following the instructions in this excellent blog post - https://vxlabs.com/2017/12/06/how-to-debug-pyinstaller-dll-pyd-load-failed-issues-on-windows/

Dependency Walker can be downloaded from here: https://github.com/lucasg/Dependencies

Remember to drill down into every dependency, as they won't pop up to the top level. For example, we originally
saw that OpenMM.dll was missing from _Move.pyd. This is only visible in Dependency Walker by clicking on _Move.pyd and 
then drilling down through its dependencies (e.g. SireMove.dll).
