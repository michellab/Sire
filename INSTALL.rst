=========================
INSTALLATION INSTRUCTIONS
=========================

Tested on Ubuntu 14.04
======================

Do you have OpenMM installed? If not make sure it is!

1. Check that zlib is installed
``sudo apt-get install zlib``

2. Check that openssl is installed ('the right package')
``sudo apt-get install libssl-dev``

3. Check your version of cmake
``cmake --version``

A minimum requirement is 2.8.12

4. Check your version of gcc
``gcc --version``

4.2<= gcc --version <=4.6


4. Create a corelib build directory and build/install
::
   mkdir build_core
   cd build_core
   cmake  -DSIRE_DISABLE_AVX=ON ../corelib
   make -j 4 #assuming you have 4 processor
   make install -j 4

5. To test that the corelib has installed correctly, check that the sire 
executable runs, e.g.

``:-> ~/sire.app/bin/sire``

(the application should run, print out some help, and then exit)

6. Create wrapper build directory and build/install
::
   mkdir build_wrapper
   cmake ../wrapper
   make -j 4
   make install -j 4

7.Note that cmake will look for sire.app in your home directory. If you
installed sire.app somewhere different, you will need to use

``:-> cmake ../../python -DSIRE_APP=/path/to/sire.app``

Check that your python3 installed properly, i.e. you didn't get an error
such as: 
::
   Python build finished, but the necessary bits to build these modules were not found:
   _ssl _xxx _xxx (some other modules here)
   To find the necessary bits, look in setup.py in detect_modules() for the module's name.
  
   Failed to build these modules:
   _ssl


If you checked step 1 and 2 there shouldn't be any problems. 

8. Now we can install setuptools and finally pip from the within sire python
::
   cd build_bundled/setuptools*/
   ~/sire.app/bundled/bin/python3 setup.py install
   cd ../../build_bundled/pip*/
   ~/sire.app/bundled/bin/python3 setup.py install

9. Now install all your favorite python packages via pip, e.g:
``~/sire.app/bundled/bin/pip install readline ipython`` 

10. To run a Sire script, e.g. script.py, simply using the Sire python 
executable, e.g.

``:-> ~/sire.app/bin/python script.py``

Sire will automatically use all of the cores in a node to parallelise the job.

Sire also comes with a set of installed scripts, that are linked to in the
sire.app/bin directory. These include the "waterswap" script. To get help
with these scripts, use "--help", e.g.

``:-> ~/sire.app/bin/waterswap --help``

11. Distributing your binaries
To package your installation of Sire up into a self-extracting
executable, type

``:-> ~/sire.app/bin/package_sire``

This will build a "sire.run" package that can be used to install Sire
on any machine that uses the same operating system, C and C++ library
as that on which you compiled the binary.

To get further help, please get in touch with the authors
via the Sire mailing lists, or via the email links on the
Sire website, http://siremol.org


