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

#ifndef SIREMM_INTRACOULOMBFF_H
#define SIREMM_INTRACOULOMBFF_H

#include "coulombpotential.h"

#include "SireFF/intra2b3dff.hpp"
#include "SireFF/intra2b2g3dff.hpp"

SIRE_BEGIN_HEADER

namespace SireMM
{

using SireFF::Intra2BFF;
using SireFF::Intra2B3DFF;
using SireFF::Intra2B2GFF;
using SireFF::Intra2B2G3DFF;

typedef Intra2BFF< CoulombPotentialInterface<IntraCoulombPotential> > 
IntraCoulombFFBase;

typedef Intra2B3DFF< CoulombPotentialInterface<IntraCoulombPotential> > 
IntraCoulombFF;

typedef Intra2B2GFF< CoulombPotentialInterface<IntraCoulombPotential> > 
IntraGroupCoulombFFBase;

typedef Intra2B2G3DFF< CoulombPotentialInterface<IntraCoulombPotential> > 
IntraGroupCoulombFF;

}

Q_DECLARE_METATYPE(SireMM::IntraCoulombFFBase);
Q_DECLARE_METATYPE(SireMM::IntraCoulombFF);

Q_DECLARE_METATYPE(SireMM::IntraGroupCoulombFFBase);
Q_DECLARE_METATYPE(SireMM::IntraGroupCoulombFF);

SIRE_EXPOSE_ALIAS( SireMM::CoulombPotentialInterface<SireMM::IntraCoulombPotential>,
                   SireMM::CoulombPotentialInterface_IntraCoulombPotential_ )

SIRE_EXPOSE_ALIAS(
  SireFF::Intra2BFF<SireMM::CoulombPotentialInterface<SireMM::IntraCoulombPotential> >,
  SireMM::IntraCoulombFFBase )

SIRE_EXPOSE_ALIAS(
  SireFF::Intra2B3DFF<SireMM::CoulombPotentialInterface<SireMM::IntraCoulombPotential> >,
  SireMM::IntraCoulombFF )
  
SIRE_EXPOSE_ALIAS(
  SireFF::Intra2B2GFF<SireMM::CoulombPotentialInterface<SireMM::IntraCoulombPotential> >,
  SireMM::IntraGroupCoulombFFBase )

SIRE_EXPOSE_ALIAS(
  SireFF::Intra2B2G3DFF<SireMM::CoulombPotentialInterface<SireMM::IntraCoulombPotential> >,
  SireMM::IntraGroupCoulombFF )
  
#ifdef SIRE_INSTANTIATE_TEMPLATES

template class SireMM::CoulombPotentialInterface<SireMM::IntraCoulombPotential>;

template class
SireFF::Intra2BFF<SireMM::CoulombPotentialInterface<SireMM::IntraCoulombPotential> >;

template class
SireFF::Intra2B3DFF<SireMM::CoulombPotentialInterface<SireMM::IntraCoulombPotential> >;

template class
SireFF::Intra2B2GFF<SireMM::CoulombPotentialInterface<SireMM::IntraCoulombPotential> >;

template class
SireFF::Intra2B2G3DFF<SireMM::CoulombPotentialInterface<SireMM::IntraCoulombPotential> >;

#endif //SIRE_INSTANTIATE_TEMPLATES

SIRE_END_HEADER

#endif
