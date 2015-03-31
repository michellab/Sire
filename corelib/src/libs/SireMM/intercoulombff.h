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

#ifndef SIREMM_INTERCOULOMBFF_H
#define SIREMM_INTERCOULOMBFF_H

#include "coulombpotential.h"

#include "SireFF/inter2b3dff.hpp"
#include "SireFF/inter2b2g3dff.hpp"

SIRE_BEGIN_HEADER

namespace SireMM
{

using SireFF::Inter2BFF;
using SireFF::Inter2B3DFF;
using SireFF::Inter2B2GFF;
using SireFF::Inter2B2G3DFF;

typedef Inter2BFF< CoulombPotentialInterface<InterCoulombPotential> > 
InterCoulombFFBase;

typedef Inter2B3DFF< CoulombPotentialInterface<InterCoulombPotential> > 
InterCoulombFF;

typedef Inter2B2GFF< CoulombPotentialInterface<InterCoulombPotential> > 
InterGroupCoulombFFBase;

typedef Inter2B2G3DFF< CoulombPotentialInterface<InterCoulombPotential> > 
InterGroupCoulombFF;

}

Q_DECLARE_METATYPE(SireMM::InterCoulombFFBase);
Q_DECLARE_METATYPE(SireMM::InterCoulombFF);

Q_DECLARE_METATYPE(SireMM::InterGroupCoulombFFBase);
Q_DECLARE_METATYPE(SireMM::InterGroupCoulombFF);

SIRE_EXPOSE_ALIAS(SireMM::CoulombPotentialInterface<SireMM::InterCoulombPotential>,
                  SireMM::CoulombPotentialInterface_InterCoulombPotential_)

SIRE_EXPOSE_ALIAS( 
    SireFF::Inter2BFF<SireMM::CoulombPotentialInterface<SireMM::InterCoulombPotential> >,
    SireMM::InterCoulombFFBase )
    
SIRE_EXPOSE_ALIAS(
    SireFF::Inter2B3DFF<SireMM::CoulombPotentialInterface<SireMM::InterCoulombPotential> >,
    SireMM::InterCoulombFF )
    
SIRE_EXPOSE_ALIAS(
    SireFF::Inter2B2GFF<SireMM::CoulombPotentialInterface<SireMM::InterCoulombPotential> >,
    SireMM::InterGroupCoulombFFBase)
    
SIRE_EXPOSE_ALIAS(
    SireFF::Inter2B2G3DFF<SireMM::CoulombPotentialInterface<SireMM::InterCoulombPotential> >,
    SireMM::InterGroupCoulombFF )

#ifdef SIRE_INSTANTIATE_TEMPLATES
template class 
SireMM::CoulombPotentialInterface<SireMM::InterCoulombPotential>;

template class 
SireFF::Inter2BFF<SireMM::CoulombPotentialInterface<SireMM::InterCoulombPotential> >;

template class
SireFF::Inter2B3DFF<SireMM::CoulombPotentialInterface<SireMM::InterCoulombPotential> >;

template class
SireFF::Inter2B2GFF<SireMM::CoulombPotentialInterface<SireMM::InterCoulombPotential> >;

template class
SireFF::Inter2B2G3DFF<SireMM::CoulombPotentialInterface<SireMM::InterCoulombPotential> >;
#endif

SIRE_END_HEADER

#endif
