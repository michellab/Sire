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

#ifndef SIREMM_INTRACLJFF_H
#define SIREMM_INTRACLJFF_H

#include "cljpotential.h"

#include "SireFF/intra2b3dff.hpp"
#include "SireFF/intra2b2g3dff.hpp"

SIRE_BEGIN_HEADER

namespace SireMM
{

using SireFF::Intra2BFF;
using SireFF::Intra2B3DFF;
using SireFF::Intra2B2GFF;
using SireFF::Intra2B2G3DFF;

typedef Intra2BFF< CLJPotentialInterface<IntraCLJPotential> > IntraCLJFFBase;
typedef Intra2B3DFF< CLJPotentialInterface<IntraCLJPotential> > IntraCLJFF;

typedef Intra2B2GFF< CLJPotentialInterface<IntraCLJPotential> > IntraGroupCLJFFBase;
typedef Intra2B2G3DFF< CLJPotentialInterface<IntraCLJPotential> > IntraGroupCLJFF;

}

Q_DECLARE_METATYPE(SireMM::IntraCLJFFBase);
Q_DECLARE_METATYPE(SireMM::IntraCLJFF);

Q_DECLARE_METATYPE(SireMM::IntraGroupCLJFFBase);
Q_DECLARE_METATYPE(SireMM::IntraGroupCLJFF);

SIRE_EXPOSE_ALIAS( SireMM::CLJPotentialInterface<SireMM::IntraCLJPotential>,
                   SireMM::CLJPotentialInterface_IntraCLJPotential_ )

SIRE_EXPOSE_ALIAS(
  SireFF::Intra2BFF<SireMM::CLJPotentialInterface<SireMM::IntraCLJPotential> >,
  SireMM::IntraCLJFFBase )

SIRE_EXPOSE_ALIAS(
  SireFF::Intra2B3DFF<SireMM::CLJPotentialInterface<SireMM::IntraCLJPotential> >,
  SireMM::IntraCLJFF )
  
SIRE_EXPOSE_ALIAS(
  SireFF::Intra2B2GFF<SireMM::CLJPotentialInterface<SireMM::IntraCLJPotential> >,
  SireMM::IntraGroupCLJFFBase )

SIRE_EXPOSE_ALIAS(
  SireFF::Intra2B2G3DFF<SireMM::CLJPotentialInterface<SireMM::IntraCLJPotential> >,
  SireMM::IntraGroupCLJFF )
  
#ifdef SIRE_INSTANTIATE_TEMPLATES

template class SireMM::CLJPotentialInterface<SireMM::IntraCLJPotential>;

template class
SireFF::Intra2BFF<SireMM::CLJPotentialInterface<SireMM::IntraCLJPotential> >;

template class
SireFF::Intra2B3DFF<SireMM::CLJPotentialInterface<SireMM::IntraCLJPotential> >;

template class
SireFF::Intra2B2GFF<SireMM::CLJPotentialInterface<SireMM::IntraCLJPotential> >;

template class
SireFF::Intra2B2G3DFF<SireMM::CLJPotentialInterface<SireMM::IntraCLJPotential> >;

#endif //SIRE_INSTANTIATE_TEMPLATES

SIRE_END_HEADER

#endif
