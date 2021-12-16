/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2018  Lester Hedges
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

#ifndef SIREIO_BIOSIMSPACE_H
#define SIREIO_BIOSIMSPACE_H

#include "sireglobal.h"

#include "SireBase/propertymap.h"

#include "SireMaths/vector.h"

#include "SireMol/select.h"

SIRE_BEGIN_HEADER

using namespace SireBase;
using namespace SireMaths;
using namespace SireMol;

namespace SireIO
{
    SIREIO_EXPORT bool isWater(const Molecule& molecule, const PropertyMap& map = PropertyMap());
    SIREIO_EXPORT bool isAmberWater(const Molecule& molecule, const PropertyMap& map = PropertyMap());
    SIREIO_EXPORT bool isGromacsWater(const Molecule& molecule, const PropertyMap& map = PropertyMap());

    SIREIO_EXPORT SelectResult setAmberWater(
            const SelectResult& molecules,
            const QString& model,
            const PropertyMap& map = PropertyMap());

    SIREIO_EXPORT SelectResult setGromacsWater(
            const SelectResult& molecules,
            const QString& model,
            const PropertyMap& map = PropertyMap());

    Vector cross(const Vector& v0, const Vector& v1);
}

SIRE_EXPOSE_FUNCTION( SireIO::isWater )
SIRE_EXPOSE_FUNCTION( SireIO::isAmberWater )
SIRE_EXPOSE_FUNCTION( SireIO::isGromacsWater )
SIRE_EXPOSE_FUNCTION( SireIO::setAmberWater )
SIRE_EXPOSE_FUNCTION( SireIO::setGromacsWater )

SIRE_END_HEADER

#endif
