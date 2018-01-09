=========================
INSTALLATION INSTRUCTIONS
=========================

1. Install dependencies

   To compile and install Sire you need to have a UNIX or UNIX-like
   environment (e.g. working bash shell), a  working C++ compiler
   that fully supports C++-14 (e.g. gcc >= 5.0 or clang >= 3.7 are recommended), a working installation of
   cmake (version 3.0.0 minimum), and a working git client, so
   that you can download the source. You also need an internet connection
   to allow you to download Sire and for the Sire build to automatically
   download all of its dependencies.

   Note, on OS X you must make sure that you have installed XCode
   and the command line developer tools. Install the tools using
   "xcode-select --install" and following the instructions. Without
   the tools, you will find that some dependencies won't compile,
   with errors like "Cannot find stdio.h"

2. Download Sire using

   ``git clone git@github.com:michellab/Sire.git``

3. Change into the resulting Sire directory

   ``cd Sire``

4. Run the script to automatically compile and install Sire

   ``./compile_sire.sh``

   This script will ask you which directory you want to install
   Sire. By default, this is ${HOME}/sire.app. You can choose anywhere
   you want. In the documentation, $SIRE will refer to this
   installation directory. You should set this as an environment
   variable, e.g.

   ``export SIRE=$HOME/sire.app``

5. Running Sire
   
   To run a Sire script, e.g. script.py, simply using the Sire python 
   executable, e.g.

   ``$SIRE/bin/python script.py``

   Sire will automatically use all of the cores in a node to parallelise the job.

   Sire also comes with a set of installed scripts, that are linked to in the
   $SIRE/bin directory. These include the "waterswap" script. To get help
   with these scripts, use "--help", e.g.

   ``$SIRE/bin/waterswap --help``

6. Distributing your binaries

   To package your installation of Sire up into a self-extracting
   executable, type

   ``$SIRE/bin/package_sire``

   This will build a "sire.run" package that can be used to install Sire
   on any machine that uses the same operating system, C and C++ library
   as that on which you compiled the binary.

   To get further help, please get in touch with the authors
   via the Sire mailing lists, or via the email links on the
   Sire website, http://siremol.org
