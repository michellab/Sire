/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2009  Christopher Woods
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

#ifndef SIREMATHS_GAMMA_H
#define SIREMATHS_GAMMA_H

#include "sireglobal.h"

SIRE_BEGIN_HEADER

namespace SireMaths
{

double Gamma(double alpha);
double gamma(double alpha);

double Gamma(double alpha, double x);
double gamma(double alpha, double x);

double incomplete_gamma_lower(double alpha, double x);
double incomplete_gamma_higher(double alpha, double x);

}

SIRE_EXPOSE_FUNCTION( SireMaths::Gamma )
SIRE_EXPOSE_FUNCTION( SireMaths::gamma )
SIRE_EXPOSE_FUNCTION( SireMaths::incomplete_gamma_lower )
SIRE_EXPOSE_FUNCTION( SireMaths::incomplete_gamma_higher )

SIRE_END_HEADER

#endif
