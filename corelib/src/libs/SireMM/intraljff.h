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

#ifndef SIREMM_INTRALJFF_H
#define SIREMM_INTRALJFF_H

#include "ljpotential.h"

#include "SireFF/intra2b3dff.hpp"
#include "SireFF/intra2b2g3dff.hpp"

SIRE_BEGIN_HEADER

namespace SireMM
{

using SireFF::Intra2BFF;
using SireFF::Intra2B3DFF;
using SireFF::Intra2B2GFF;
using SireFF::Intra2B2G3DFF;

typedef Intra2BFF< LJPotentialInterface<IntraLJPotential> > IntraLJFFBase;
typedef Intra2B3DFF< LJPotentialInterface<IntraLJPotential> > IntraLJFF;

typedef Intra2B2GFF< LJPotentialInterface<IntraLJPotential> > IntraGroupLJFFBase;
typedef Intra2B2G3DFF< LJPotentialInterface<IntraLJPotential> > IntraGroupLJFF;

}

Q_DECLARE_METATYPE(SireMM::IntraLJFFBase);
Q_DECLARE_METATYPE(SireMM::IntraLJFF);

Q_DECLARE_METATYPE(SireMM::IntraGroupLJFFBase);
Q_DECLARE_METATYPE(SireMM::IntraGroupLJFF);

SIRE_EXPOSE_ALIAS( SireMM::LJPotentialInterface<SireMM::IntraLJPotential>,
                   SireMM::LJPotentialInterface_IntraLJPotential_ )

SIRE_EXPOSE_ALIAS(
  SireFF::Intra2BFF<SireMM::LJPotentialInterface<SireMM::IntraLJPotential> >,
  SireMM::IntraLJFFBase )

SIRE_EXPOSE_ALIAS(
  SireFF::Intra2B3DFF<SireMM::LJPotentialInterface<SireMM::IntraLJPotential> >,
  SireMM::IntraLJFF )
  
SIRE_EXPOSE_ALIAS(
  SireFF::Intra2B2GFF<SireMM::LJPotentialInterface<SireMM::IntraLJPotential> >,
  SireMM::IntraGroupLJFFBase )

SIRE_EXPOSE_ALIAS(
  SireFF::Intra2B2G3DFF<SireMM::LJPotentialInterface<SireMM::IntraLJPotential> >,
  SireMM::IntraGroupLJFF )
  
#ifdef SIRE_INSTANTIATE_TEMPLATES

template class SireMM::LJPotentialInterface<SireMM::IntraLJPotential>;

template class
SireFF::Intra2BFF<SireMM::LJPotentialInterface<SireMM::IntraLJPotential> >;

template class
SireFF::Intra2B3DFF<SireMM::LJPotentialInterface<SireMM::IntraLJPotential> >;

template class
SireFF::Intra2B2GFF<SireMM::LJPotentialInterface<SireMM::IntraLJPotential> >;

template class
SireFF::Intra2B2G3DFF<SireMM::LJPotentialInterface<SireMM::IntraLJPotential> >;

#endif //SIRE_INSTANTIATE_TEMPLATES

SIRE_END_HEADER

#endif
