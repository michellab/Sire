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

#ifndef SIREMM_INTERCLJFF_H
#define SIREMM_INTERCLJFF_H

#include "cljpotential.h"

#include "SireFF/inter2b3dff.hpp"
#include "SireFF/inter2b2g3dff.hpp"

SIRE_BEGIN_HEADER

namespace SireMM
{

using SireFF::Inter2BFF;
using SireFF::Inter2B3DFF;
using SireFF::Inter2B2GFF;
using SireFF::Inter2B2G3DFF;

typedef Inter2BFF< CLJPotentialInterface<InterCLJPotential> > InterCLJFFBase;
typedef Inter2B3DFF< CLJPotentialInterface<InterCLJPotential> > InterCLJFF;

typedef Inter2B2GFF< CLJPotentialInterface<InterCLJPotential> > InterGroupCLJFFBase;
typedef Inter2B2G3DFF< CLJPotentialInterface<InterCLJPotential> > InterGroupCLJFF;

}

Q_DECLARE_METATYPE(SireMM::InterCLJFFBase);
Q_DECLARE_METATYPE(SireMM::InterCLJFF);

Q_DECLARE_METATYPE(SireMM::InterGroupCLJFFBase);
Q_DECLARE_METATYPE(SireMM::InterGroupCLJFF);

SIRE_EXPOSE_ALIAS(SireMM::CLJPotentialInterface<SireMM::InterCLJPotential>,
                  SireMM::CLJPotentialInterface_InterCLJPotential_)

SIRE_EXPOSE_ALIAS( 
    SireFF::Inter2BFF<SireMM::CLJPotentialInterface<SireMM::InterCLJPotential> >,
    SireMM::InterCLJFFBase )
    
SIRE_EXPOSE_ALIAS(
    SireFF::Inter2B3DFF<SireMM::CLJPotentialInterface<SireMM::InterCLJPotential> >,
    SireMM::InterCLJFF )
    
SIRE_EXPOSE_ALIAS(
    SireFF::Inter2B2GFF<SireMM::CLJPotentialInterface<SireMM::InterCLJPotential> >,
    SireMM::InterGroupCLJFFBase)
    
SIRE_EXPOSE_ALIAS(
    SireFF::Inter2B2G3DFF<SireMM::CLJPotentialInterface<SireMM::InterCLJPotential> >,
    SireMM::InterGroupCLJFF )

#ifdef SIRE_INSTANTIATE_TEMPLATES
template class 
SireMM::CLJPotentialInterface<SireMM::InterCLJPotential>;

template class 
SireFF::Inter2BFF<SireMM::CLJPotentialInterface<SireMM::InterCLJPotential> >;

template class
SireFF::Inter2B3DFF<SireMM::CLJPotentialInterface<SireMM::InterCLJPotential> >;

template class
SireFF::Inter2B2GFF<SireMM::CLJPotentialInterface<SireMM::InterCLJPotential> >;

template class
SireFF::Inter2B2G3DFF<SireMM::CLJPotentialInterface<SireMM::InterCLJPotential> >;
#endif

SIRE_END_HEADER

#endif
