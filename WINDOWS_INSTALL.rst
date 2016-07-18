== Instructions for compiling and installing on windows ==

First, install cygwin. With cygwin you must install

 * bash
 * cmake
 * gcc (core, c++)
 * git
 * curl
 * automake
 * make

Then, download Sire using git, e.g.

 # git clone https://github.com/michellib/Sire.git

Next, change into the Sire directory

 # cd Sire

Next, use the "compile_sire.sh" script to compile and install Sire into ~/sire.app

 # ./compile_sire.sh

This will download miniconda and will attempt to install it. Note that you will
be prompted with a graphical interface to install miniconda. Tell the installer
to install a copy of miniconda just for you, and choose the install directory
to be the one you have chosen to install Sire (typically $HOME/sire.app).

For me, this would be C:\cygwin\home\chris\sire.app

Then choose to not add this to the path, and to not register anaconda as the default python.



