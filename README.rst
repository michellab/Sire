***********************************************
Sire is a versatile Molecular Mechanics library
***********************************************

.. image:: https://travis-ci.org/michellab/Sire.svg?branch=autobuild_feature
   :target: https://travis-ci.org/michellab/Sire


About
=====
Sire is a free, open source, multiscale molecular simulation framework, written to allow computational modelers to quickly prototype and develop new algorithms for molecular simulation and molecular design. Sire is written as a collection of libraries, each of which contains self-contained and robust C++/Python building blocks. These building blocks are vectorised and thread-aware and can be streamed (saved/loaded) to and from a version-controlled and tagged binary format, thereby allowing them to be combined together easily to build custom multi-processor molecular simulation applications.

Installation 
============

The easy install option is:
```
$: git clone git@github.com:michellab/Sire.git
$: cd Sire
$: ./compile_sire.sh
```
A small word of warning, the installation can easily take up to an hour!

Alternatively Sire can be installed into an existin Conda installation/Python installation
For these instructions plese refer to the respective INSTALL files. 


Support and Development
=======================

Bugs, Comments, Questions
--------------------------
For bug reports/sugguestions/complains please file an issue on 
`GitHub <http://github.com/michellab/Sire>`__.
or contact the developers via the google user group: `https://groups.google.com/forum/#!forum/sire-users`

Developers guide
-----------------
Please refer to the developers guide for active code development. 


Travis -- Autobuild feature
---------------------------

Since Sire is quite large, a build can take quite long and might not be neccessary if a commit is only fixing a couple of typos. Simply add the line `[ci skip]` to your commit message and Travis will not invoke an autobuild. 


External Libraries
------------------
* mdtraj (LGPLv3): https://mdtraj.org
* openMM (LGPLv3): https://openmm.org
