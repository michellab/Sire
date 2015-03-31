/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2008  Christopher Woods
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

#ifndef SIREMM_INTRASOFTCLJFF_H
#define SIREMM_INTRASOFTCLJFF_H

#include "softcljpotential.h"

#include "SireFF/intra2b3dff.hpp"
#include "SireFF/intra2b2g3dff.hpp"

SIRE_BEGIN_HEADER

namespace SireMM
{

using SireFF::Intra2BFF;
using SireFF::Intra2B3DFF;
using SireFF::Intra2B2GFF;
using SireFF::Intra2B2G3DFF;

typedef Intra2BFF< SoftCLJPotentialInterface<IntraSoftCLJPotential> > IntraSoftCLJFFBase;
typedef Intra2B3DFF< SoftCLJPotentialInterface<IntraSoftCLJPotential> > IntraSoftCLJFF;

typedef Intra2B2GFF< SoftCLJPotentialInterface<IntraSoftCLJPotential> > IntraGroupSoftCLJFFBase;
                                                                
typedef Intra2B2G3DFF< SoftCLJPotentialInterface<IntraSoftCLJPotential> > IntraGroupSoftCLJFF;

}

Q_DECLARE_METATYPE(SireMM::IntraSoftCLJFFBase);
Q_DECLARE_METATYPE(SireMM::IntraSoftCLJFF);

Q_DECLARE_METATYPE(SireMM::IntraGroupSoftCLJFFBase);
Q_DECLARE_METATYPE(SireMM::IntraGroupSoftCLJFF);

SIRE_EXPOSE_ALIAS(SireMM::CLJPotentialInterface<SireMM::IntraSoftCLJPotential>,
                  SireMM::CLJPotentialInterface_IntraSoftCLJPotential_)

SIRE_EXPOSE_ALIAS(SireMM::SoftCLJPotentialInterface<SireMM::IntraSoftCLJPotential>,
                  SireMM::SoftCLJPotentialInterface_IntraSoftCLJPotential_)

SIRE_EXPOSE_ALIAS( 
    SireFF::Intra2BFF<SireMM::SoftCLJPotentialInterface<SireMM::IntraSoftCLJPotential> >,
    SireMM::IntraSoftCLJFFBase )
    
SIRE_EXPOSE_ALIAS(
    SireFF::Intra2B3DFF<SireMM::SoftCLJPotentialInterface<SireMM::IntraSoftCLJPotential> >,
    SireMM::IntraSoftCLJFF )
    
SIRE_EXPOSE_ALIAS(
 SireFF::Intra2B2GFF<SireMM::SoftCLJPotentialInterface<SireMM::IntraSoftCLJPotential> >,
 SireMM::IntraGroupSoftCLJFFBase)
    
SIRE_EXPOSE_ALIAS(
 SireFF::Intra2B2G3DFF<SireMM::SoftCLJPotentialInterface<SireMM::IntraSoftCLJPotential> >,
 SireMM::IntraGroupSoftCLJFF )

#ifdef SIRE_INSTANTIATE_TEMPLATES
template class 
SireMM::SoftCLJPotentialInterface<SireMM::IntraSoftCLJPotential>;

template class 
SireFF::Intra2BFF<SireMM::SoftCLJPotentialInterface<SireMM::IntraSoftCLJPotential> >;

template class
SireFF::Intra2B3DFF<SireMM::SoftCLJPotentialInterface<SireMM::IntraSoftCLJPotential> >;

template class
SireFF::Intra2B2GFF<SireMM::SoftCLJPotentialInterface<SireMM::IntraSoftCLJPotential> >;

template class
SireFF::Intra2B2G3DFF<SireMM::SoftCLJPotentialInterface<SireMM::IntraSoftCLJPotential> >;
#endif

SIRE_END_HEADER

#endif
