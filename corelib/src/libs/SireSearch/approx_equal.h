/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2022  Christopher Woods
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

#ifndef SIRESEARCH_APPROX_EQUAL_H
#define SIRESEARCH_APPROX_EQUAL_H

#include "sireglobal.h"

SIRE_BEGIN_HEADER

namespace SireSearch
{
    bool SIRESEARCH_EXPORT approx_equal(double val0, double val1);

    double SIRESEARCH_EXPORT get_approx_epsilon();
    void SIRESEARCH_EXPORT set_approx_epsilon(double eps);
}

SIRE_EXPOSE_FUNCTION( SireSearch::approx_equal )
SIRE_EXPOSE_FUNCTION( SireSearch::get_approx_epsilon )
SIRE_EXPOSE_FUNCTION( SireSearch::set_approx_epsilon )

SIRE_END_HEADER

#endif
