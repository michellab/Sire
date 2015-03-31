/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2007  Christopher Woods
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

#ifndef SIREMM_INTERLJFF_H
#define SIREMM_INTERLJFF_H

#include "ljpotential.h"

#include "SireFF/inter2b3dff.hpp"
#include "SireFF/inter2b2g3dff.hpp"

SIRE_BEGIN_HEADER

namespace SireMM
{

using SireFF::Inter2BFF;
using SireFF::Inter2B3DFF;
using SireFF::Inter2B2GFF;
using SireFF::Inter2B2G3DFF;

typedef Inter2BFF< LJPotentialInterface<InterLJPotential> > InterLJFFBase;
typedef Inter2B3DFF< LJPotentialInterface<InterLJPotential> > InterLJFF;

typedef Inter2B2GFF< LJPotentialInterface<InterLJPotential> > InterGroupLJFFBase;
typedef Inter2B2G3DFF< LJPotentialInterface<InterLJPotential> > InterGroupLJFF;

}

Q_DECLARE_METATYPE(SireMM::InterLJFFBase);
Q_DECLARE_METATYPE(SireMM::InterLJFF);

Q_DECLARE_METATYPE(SireMM::InterGroupLJFFBase);
Q_DECLARE_METATYPE(SireMM::InterGroupLJFF);

SIRE_EXPOSE_ALIAS(SireMM::LJPotentialInterface<SireMM::InterLJPotential>,
                  SireMM::LJPotentialInterface_InterLJPotential_)

SIRE_EXPOSE_ALIAS( 
    SireFF::Inter2BFF<SireMM::LJPotentialInterface<SireMM::InterLJPotential> >,
    SireMM::InterLJFFBase )
    
SIRE_EXPOSE_ALIAS(
    SireFF::Inter2B3DFF<SireMM::LJPotentialInterface<SireMM::InterLJPotential> >,
    SireMM::InterLJFF )
    
SIRE_EXPOSE_ALIAS(
    SireFF::Inter2B2GFF<SireMM::LJPotentialInterface<SireMM::InterLJPotential> >,
    SireMM::InterGroupLJFFBase)
    
SIRE_EXPOSE_ALIAS(
    SireFF::Inter2B2G3DFF<SireMM::LJPotentialInterface<SireMM::InterLJPotential> >,
    SireMM::InterGroupLJFF )

#ifdef SIRE_INSTANTIATE_TEMPLATES
template class 
SireMM::LJPotentialInterface<SireMM::InterLJPotential>;

template class 
SireFF::Inter2BFF<SireMM::LJPotentialInterface<SireMM::InterLJPotential> >;

template class
SireFF::Inter2B3DFF<SireMM::LJPotentialInterface<SireMM::InterLJPotential> >;

template class
SireFF::Inter2B2GFF<SireMM::LJPotentialInterface<SireMM::InterLJPotential> >;

template class
SireFF::Inter2B2G3DFF<SireMM::LJPotentialInterface<SireMM::InterLJPotential> >;
#endif

SIRE_END_HEADER

#endif
