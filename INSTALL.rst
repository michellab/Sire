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

4. Create a corelib build directory and build/install
::
   mkdir build_core
   cd build_core
   cmake  -DSIRE_DISABLE_AVX=ON ../corelib
   make -j 4 #assuming you have 4 processor
   make install -j 4


5. Create wrapper build directory and build/install
::
   mkdir build_wrapper
   cmake ../wrapper
   make -j 4
   make install -j 4

Check that your python3 installed properly, i.e. you didn't get an error
such as: 
::
   Python build finished, but the necessary bits to build these modules were not found:
   _ssl _xxx _xxx (some other modules here)
   To find the necessary bits, look in setup.py in detect_modules() for the module's name.
  
   Failed to build these modules:
   _ssl


If you checked step 1 and 2 there shouldn't be any problems. 

6. Now we can install setuptools and finally pip from the within sire python
::
   cd build_bundled/setuptools*/
   ~/sire.app/bundled/bin/python3 setup.py install
   cd ../../build_bundled/pip*/
   ~/sire.app/bundled/bin/python3 setup.py install

7. Now install all your favorite python packages via pip, e.g:
``~/sire.app/bundled/bin/pip install ipython`` 


