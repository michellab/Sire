/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2006  Christopher Woods
  *
  *  This program is free software; you can redistribute it and/or modify
  *  it under the terms of the GNU General Public License as published by
  *  the Free Software Foundation; either version 2 of the License, or
  *  (at your option) any later version.
  *
  *  This program is distributed in the hope that it will be useful,
  *  but WITHOUT ANY WARRANTY; without even the implied warranty of
  *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  *  GNU General Public License for more details.
  *
  *  You should have received a copy of the GNU General Public License
  *  along with this program; if not, write to the Free Software
  *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
  *
  *  For full details of the license please see the COPYING file
  *  that should have come with this distribution.
  *
  *  You can contact the authors via the developer's mailing list
  *  at http://siremol.org
  *
\*********************************************/

#ifndef SIREMATHS_CONSTANTS_H
#define SIREMATHS_CONSTANTS_H

#include <gsl/gsl_math.h>

#include <limits>

#include "sireglobal.h"

SIRE_BEGIN_HEADER

/** expose GSL #define constants as const double values */
namespace SireMaths
{

const double e = M_E; /* e */

const double log2_e = M_LOG2E;   /* log_2 (e) */

const double log10_e = M_LOG10E;  /* log_10 (e) */

const double sqrt_two = M_SQRT2;  /* sqrt(2) */

const double sqrt_half = M_SQRT1_2;   /* sqrt(1/2) */

const double sqrt_three = M_SQRT3;   /* sqrt(3) */

const double pi = 3.141592653589793238462643383279;    /* pi */

const double two_pi = double(2) * pi;  /* 2 * pi */

const double pi_over_two = M_PI_2;   /* pi/2 */

const double pi_4 = M_PI_4;   /* pi/4 */

const double sqrtpi = M_SQRTPI;    /* sqrt(pi) */

const double two_sqrtpi = M_2_SQRTPI;  /* 2/sqrt(pi) */

const double one_over_pi = M_1_PI;   /* 1/pi */

const double two_over_pi = M_2_PI;   /* 2/pi */

const double ln_ten = M_LN10;    /* ln(10) */

const double ln_two = M_LN2;    /* ln(2) */

const double ln_pi = M_LNPI;    /* ln(pi) */

const double euler = M_EULER;   /* Euler constant */

//now expose some sizes...
/** A small number */
const double small = 0.001;

/** A tiny number */
const double tiny = 1.0e-8;

/** The smallest possible number */
const double smallest = std::numeric_limits<double>::min();

/** A large number */
const double large = 1000.0;

/** A huge number */
const double huge = 1.0e8;

/** The largest possible number */
const double largest = std::numeric_limits<double>::max();

/** The smallest possible int */
const int smallest_int = std::numeric_limits<int>::min();

/** The largest possible int */
const int largest_int = std::numeric_limits<int>::max();

}

SIRE_END_HEADER

#endif
