****
`Sire <http://siremol.org>`__
****

.. image:: https://travis-ci.org/michellab/Sire.svg?branch=autobuild_feature
   :target: https://travis-ci.org/michellab/Sire


About
=====
`Sire <http://siremol.org>`__ is a free, open source, multiscale 
molecular simulation framework, written to allow computational 
modellers to quickly prototype and develop new algorithms for 
molecular simulation and molecular design. Sire is written 
as a collection of libraries, each of which contains self-contained 
and robust C++/Python building blocks. These building blocks are v
ectorised and thread-aware and can be streamed (saved/loaded) 
to and from a version-controlled and tagged binary format, 
thereby allowing them to be combined together easily to build 
custom multi-processor molecular simulation applications.

For more information about how to use Sire, and about application
built with Sire, please `visit the Sire website <http://siremol.org>`__.

Installation 
============

There are many `pre-built binary packages <http://siremol.org/Sire/Binaries.html>`__,
which are available for Linux and OS X, which are quick and easy to install.

However, as you are here, it is likely you want to download the latest,
greatest version of the code, which you will need to compile. To compile Sire,
you need a working C++ compiler (gcc or clang), `cmake <http://cmake.org>`__ 
(version 2.8.11.2 or above), a Git client to download the source,
and a working internet connection (needed by
the Sire compilation scripts to download additional dependencies).

The easy install option is::

    git clone git@github.com:michellab/Sire.git
    cd Sire
    ./compile_sire.sh

A small word of warning, the compilation can easily take over an hour!

The above will download and install a new miniconda python installation,
into which Sire will be compiled and deployed (together with its
dependencies). This is by far the easiest way to compile and install Sire,
and is the route we strongly recommend. If you have any problems with 
compiling and installing Sire, then please get in touch using the links below.

If you want to install Sire into an existing miniconda or Anaconda 
Python installation, please follow the instructions in `build/INSTALL_INTO_ANACONDA.rst`.

Support and Development
=======================

Bugs, Comments, Questions
--------------------------
For bug reports/sugguestions/complains please file an issue on 
`GitHub <http://github.com/michellab/Sire>`__.
or contact the developers via the google user group: `https://groups.google.com/forum/#!forum/sire-users`

Developers guide
-----------------
Please `visit the website <http://siremol.org>`__ for information on how to 
develop applications using Sire. 

Travis -- Autobuild feature
---------------------------
Since Sire is quite large, a build can take quite long and might not be neccessary 
if a commit is only fixing a couple of typos. Simply add the line `[ci skip]` 
to your commit message and Travis will not invoke an autobuild. 
