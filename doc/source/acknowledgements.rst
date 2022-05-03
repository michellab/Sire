================
Acknowledgements
================

We gratefully acknowledge funding from the EPSRC and BBSRC who contributed
funding to the development of Sire and related software.

We also thank the `Software Sustainability Institute <https://software.ac.uk>`__
for many useful discussions.

We thank the Universities of Bristol and Edinburgh for providing the
time to the members of staff who have contributed to Sire's development.

We thank `CCP-BioSim <https://ccpbiosim.ac.uk>`__ who have also provided
guidance and encouragement during the development of this software.

Third Party Software
====================

Sire depends on a lot of third party software, the details and licenses of
which can be found below. The software will be installed automatically
as part of the installing the Sire conda package, so you shouldn't
have to do anything yourself.

Sire is itself distributed under the terms of the GPL version 2.

Qt 5
----

Sire is built on top of Qt.Core from Qt 5. This is used under the terms
of the `LGPL 2 <http://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html>`__.
Sire does not modify Qt, so this license allows both
commercial and non-commercial use without fee. You can find out more
about Qt and its license `from here <https://www.qt.io/>`__.

Python
------

The Sire C++ library is wrapped up and made available for use within Python 3.
This is used under the terms of the `PSF license <https://docs.python.org/3/license.html>`__,
which is compatible with
the GPLv2. The license allows both commercial and non-commercial use
without fee. You can find out more about Python and its license
`from here <https://www.python.org/>`__.

boost
-----

Sire uses many of the components from the boost libraries, in particular
the boost::python module that is used to wrap up the C++ code.
This is used under the terms of the
`Boost Software License <http://www.boost.org/users/license.html>`__,
which allows both commercial and non-commercial
use without fee. This is compatible with the GPL. You can find out
more about boost and its license `from here <http://www.boost.org/>`__.

Py++
----

Sire uses Py++ to auto-generate all of the C++ python wrappers. Py++ uses
either GCCXML or CastXML, and, as it is used as a tool, its license does
not affect Sire. You can read more about Py++
`from here <http://pyplusplus.readthedocs.io/en/latest/>`__.

cmake
-----

Sire uses cmake as its build system. As it is used as a tool, its license
does not affect Sire. CMake is excellent. You can read more about it
`from here <https://cmake.org/>`__.

Anaconda
---------

Sire uses Anaconda Python (specifically miniconda and conda-forge) to
simplify the management and installation of Python and the various
modules on which Sire depends.

Anaconda (and miniconda) are distributed as
`open source projects <https://www.continuum.io/open-source-core-modern-software>`__.
As Sire does not explicitly link with them, the license is not an issue.
You can find out more about Anaconda `from here <https://www.continuum.io/>`__.

Threading Building Blocks (tbb)
-------------------------------

Sire uses the `Threading Building Blocks <https://www.threadingbuildingblocks.org/>`__
library for within-node
parallelisation. This is licensed under the open source
`Apache 2.0 license <https://www.threadingbuildingblocks.org/faq/10>`__.

Gnu Scientific Library (GSL)
----------------------------

Sire uses some of the routines from the Gnu Scientific Library.
This is used under the terms of the `GPL v3 <http://www.gnu.org/copyleft/gpl.html>`__
license. More information
about GSL and its license can be `found here <http://www.gnu.org/software/gsl/>`__.

Regress
-------

Sire uses the linear least squares regression library, `regress`, for
polynomial least squares fitting. This is used under the terms of
the GPLv3 license.

The source code for this module can be
`found here <https://github.com/michellab/Sire/blob/devel/corelib/src/libs/SireAnalysis/third_party/regress.cpp>`__.

eig3
----

Sire uses the eig3 library for eigenvector/eigenmatrix calculations by
Connelly Barnes. This is in the public domain, and is derived itself
from the Java matrix library JAMA (also public domain).

Information about this can be
`found here <http://barnesc.blogspot.co.uk/2007/02/eigenvectors-of-3x3-symmetric-matrix.html>`__,
with the license within Sire `found here <https://github.com/michellab/Sire/blob/devel/corelib/src/libs/SireMaths/third_party/eig3/readme.txt>`__.

Mersenne Twister
----------------

Sire uses the Mersenne Twister program by Richard Wagner for the generation
of random numbers. This is used under a BSD-style license, shown below.

::

 // Mersenne Twister random number generator -- a C++ class MTRand
 // Based on code by Makoto Matsumoto, Takuji Nishimura, and Shawn Cokus
 // Richard J. Wagner  v1.0  15 May 2003  rjwagner@writeme.com

 // The Mersenne Twister is an algorithm for generating random numbers.  It
 // was designed with consideration of the flaws in various other generators.
 // The period, 2^19937-1, and the order of equidistribution, 623 dimensions,
 // are far greater.  The generator is also fast; it avoids multiplication and
 // division, and it benefits from caches and pipelines.  For more information
 // see the inventors' web page at http://www.math.keio.ac.jp/~matumoto/emt.html

 // Reference
 // M. Matsumoto and T. Nishimura, "Mersenne Twister: A 623-Dimensionally
 // Equidistributed Uniform Pseudo-Random Number Generator", ACM Transactions on
 // Modeling and Computer Simulation, Vol. 8, No. 1, January 1998, pp 3-30.

 // Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
 // Copyright (C) 2000 - 2003, Richard J. Wagner
 // All rights reserved.
 //
 // Redistribution and use in source and binary forms, with or without
 // modification, are permitted provided that the following conditions
 // are met:
 //
 //   1. Redistributions of source code must retain the above copyright
 //      notice, this list of conditions and the following disclaimer.
 //
 //   2. Redistributions in binary form must reproduce the above copyright
 //      notice, this list of conditions and the following disclaimer in the
 //      documentation and/or other materials provided with the distribution.
 //
 //   3. The names of its contributors may not be used to endorse or promote
 //      products derived from this software without specific prior written
 //      permission.
 //
 // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 // "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 // LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 // A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 // CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 // EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 // PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 // PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 // LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 // NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 // SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 // The original code included the following notice:
 //
 //     When you use this, send an email to: matumoto@math.keio.ac.jp
 //     with an appropriate reference to your work.
 //
 // It would be nice to CC: rjwagner@writeme.com and Cokus@math.washington.edu
 // when you write.

I must remember to send them an email…

More information about Mersenne Twister can be
`found here <http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/VERSIONS/C-LANG/c-lang.html>`__.

sse_mathfun and neon_mathfun
----------------------------

Sire uses sse_mathfun and neon_mathfun for vectorising intrinsic maths
functions on processors that support SSE or Neon. These libraries were written
by Julien Pommier, and released under the BSD-style zlib license,
which is given here.

::

 /* Copyright (C) 2007  Julien Pommier
  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the authors be held liable for any damages
  arising from the use of this software.
  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:
  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.
  (this is the zlib license)
 */


avx_mathfun
-----------

This is an AVX library inspired by sse_mathfun, that extends support to
processors with AVX instructions. It was written by Giovanni Garberoglio,
and is also under a BSD-style zlib license.

::

  AVX implementation of sin, cos, sincos, exp and log
   Based on "sse_mathfun.h", by Julien Pommier
   http://gruntthepeon.free.fr/ssemath/
   Copyright (C) 2012 Giovanni Garberoglio
   Interdisciplinary Laboratory for Computational Science (LISC)
   Fondazione Bruno Kessler and University of Trento
   via Sommarive, 18
   I-38123 Trento (Italy)
  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the authors be held liable for any damages
  arising from the use of this software.
  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:
  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.
  (this is the zlib license)

BLAS
----

Sire bundles the reference version of the BLAS library. The license for which
is described below (it is ‘freely-available’).

::

 <http://www.netlib.org/blas/faq.html#2>

 The reference BLAS is a freely-available software package.
 It is available from netlib via anonymous ftp and the World Wide Web.
 Thus, it can be included in commercial software packages (and has been).
 We only ask that proper credit be given to the authors.

 Like all software, it is copyrighted. It is not trademarked, but we do ask the following:

 If you modify the source for these routines we ask that you change the name of the
 routine and comment the changes made to the original.

 We will gladly answer any questions regarding the software. If a modification is done,
 however, it is the responsibility of the person who modified the routine to provide support.

LAPACK
------

Sire also bundles a reference LAPACK, the license for which is given below;

::

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are
 met:

 - Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.

 - Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer listed
   in this license in the documentation and/or other materials
   provided with the distribution.

 - Neither the name of the copyright holders nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

LINPACK
-------

Sire also bundles some routines from LINPACK. Forum postings indicate
that LINPACK is BSD licensed.

::

 LINPACK is a collection of Fortran subroutines that analyze and
 solve linear equations and linear least-squares probles.  The
 package solves linear systems whose matrices are general, banded,
 symmetric indefinite, symmetric positive definite, triangular,
 and tridiagonal square.  In addition, the package computes
 the QR and singular value decompositions of rectangular matrices
 and applies them to least-squares problems.  LINPACK uses
 column-oriented algorithms to increase efficiency by preserving
 locality of reference.

 LINPACK was designed for supercomputers in use in the 1970s and
 early 1980s.  LINPACK has been largely superceded by LAPACK
 which has been designed to run efficiently on shared-memory, vector
 supercomputers.

 Developed by Jack Dongarra, Jim Bunch, Cleve Moler and Pete Stewart.
  1 Feb 84

 If you are interested in acquiring the entire LINPACK, it may
 make more sense to talk with NAG. NAG distribute the software
 on a mag tape for a nominal charge.
     NAG
     1400 Opus Place, Suite 200
     Downers Grove, IL  60515-5702
     708-971-2337, FAX 971-2706

See ThirdParty/LAPACK for BSD LAPACK license. Forum postings indicate
that LINPACK is also BSD licensed (e.g.
http://icl.cs.utk.edu/lapack-forum/archives/lapack/msg00301.html)

::

 [Lapack] Linpack license?
 From: Jakub Kurzak
 Date: Wed, 11 Jul 2007 09:50:02 -0400
 Basically it is BSD.
 Jakub

 On 7/10/07, Benjamin Collar <benjamin.collar@Domain.Removed> wrote:

 -----BEGIN PGP SIGNED MESSAGE-----
 Hash: SHA1

 Hello,

 I am trying to find out what license Linpack is released under. I
 checked the netlib/linpack site, but didn't find the answer there. Do
 you know?

 Thanks
 Benjamin
 -----BEGIN PGP SIGNATURE-----
 Version: GnuPG v1.4.7 (MingW32)
 Comment: http://firegpg.tuxfamily.org

 iD8DBQFGlE594CA25YTAYOMRAiG8AJ4uLe6DgAIszvyZl7fUdcmdv960MACdF1mA
 jfqWri2vpHcKrlsur6Uc1yg=
 =Xb40
 -----END PGP SIGNATURE-----
 _______________________________________________
 Lapack mailing list
 Lapack@Domain.Removed
 http://lists.cs.utk.edu/listinfo/lapack

 -------------- next part --------------
 An HTML attachment was scrubbed...
 URL:
 http://lists.cs.utk.edu/private/lapack/attachments/20070711/cbd659b6/attachment.htm

LAP (Linear Assignment Problem Solver)
--------------------------------------

Sire implements its own C++ version of the LAP library for solving the
linear assignment problem. This is `available here <https://github.com/michellab/Sire/blob/devel/corelib/src/libs/SireMaths/linearap.cpp>`__.

The original code is Freeware, with more information about it available
`from here <http://www.assignmentproblems.com/linearAP.htm>`__.

MD5
---

Sire uses the MD5 library written by L. Peter Deutsch.
It is used under a BSD-style license, given below.

::

  Copyright (C) 1999, 2002 Aladdin Enterprises.  All rights reserved.
  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the authors be held liable for any damages
  arising from the use of this software.
  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:
  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.
  L. Peter Deutsch
  ghost@aladdin.com

More information about MD5 libraries in general can be
`found here <http://userpages.umbc.edu/~mabzug1/cs/md5/md5.html>`__.

kabasch fitting
---------------

I have written a C++ implementation of the kabasch algorithm for alignment.
This was inspired by the calculate_rmsd python script written by
Jimmy Charnley Kromann and Lars Bratholm,
available https://github.com/charnley/rmsd, and under license;

::

        =====================
        Copyright (c) 2013, Jimmy Charnley Kromann <jimmy@charnley.dk> & Lars Bratholm
        All rights reserved.

        Redistribution and use in source and binary forms, with or without
        modification, are permitted provided that the following conditions are met:

        1. Redistributions of source code must retain the above copyright notice, this
           list of conditions and the following disclaimer.
        2. Redistributions in binary form must reproduce the above copyright notice,
           this list of conditions and the following disclaimer in the documentation
           and/or other materials provided with the distribution.

        THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
        ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
        WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
        DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
        ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
        (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
        LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
        ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
        (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
        SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
        =======================

ap (ascii plot)
---------------

Sire bundles the Python “ap” library for drawing ascii graphs.
This is available as “Sire.Tools.ap”

Version 0.9 written by M. Fouesneau is included, available
freely from `GitHub here <https://github.com/mfouesneau/asciiplot>`__.
The only change I've made is running this through Python's 2to3
program to make this code work with Python 3.

The header documentation reads;

::

 Package that allows you to plot simple graphs in ASCII, a la matplotlib.
 This package is a inspired from Imri Goldberg's ASCII-Plotter 1.0
 (https://pypi.python.org/pypi/ASCII-Plotter/1.0)
 At a time I was enoyed by security not giving me direct access to my computer,
 and thus to quickly make figures from python, I looked at how I could make
 quick and dirty ASCII figures. But if I were to develop something, I wanted
 something that can be used with just python and possible standard-ish packages
 (numpy, scipy).
 So I came up with this package after many iterations based of ASCII-plotter.
 I added the feature to show multiple curves on one plot with different markers.
 And I also made the usage, close to matplotlib, such that there is a plot,
 hist, hist2d and imshow functions.

 TODO:
     imshow does not plot axis yet.
     make a correct documentation
