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

#ifndef SIREMOL_ATOMFORCES_H
#define SIREMOL_ATOMFORCES_H

#include "atomproperty.hpp"

#include "SireMaths/vector3d.hpp"

#include "SireUnits/dimensions.h"

SIRE_BEGIN_HEADER

namespace SireMol
{

using SireUnits::Dimension::Force;

typedef SireMaths::Vector3D<Force> Force3D;

typedef AtomProperty<Force3D> AtomForces;

}

SIRE_EXPOSE_ALIAS( 
        (SireMaths::Vector3D<SireUnits::Dimension::PhysUnit<1, 1, -2, 0, 0, 0, 0> >),
        SireMol::Force3D )
                   
Q_DECLARE_METATYPE( SireMol::AtomForces );
Q_DECLARE_METATYPE( SireMol::Force3D );

SIRE_EXPOSE_ATOM_PROPERTY( SireMaths::Vector3D<SireUnits::Dimension::Force>,
                           SireMol::AtomForces )

#ifdef SIRE_INSTANTIATE_TEMPLATES
template class SireMaths::Vector3D<SireUnits::Dimension::Force>;
template class SireMol::AtomProperty<SireMol::Force3D>;
#endif

SIRE_END_HEADER

#endif
