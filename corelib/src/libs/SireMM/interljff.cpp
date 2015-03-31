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

#include "interljff.h"

#include "SireMol/partialmolecule.h"

#include "SireMol/mover.hpp"

using namespace SireMM;
using namespace SireFF;

namespace SireMM
{
    template class LJPotentialInterface<InterLJPotential>;
}

namespace SireFF
{
    template class Inter2BFF< LJPotentialInterface<InterLJPotential> >;
    template class Inter2B3DFF< LJPotentialInterface<InterLJPotential> >;

    template class Inter2B2GFF< LJPotentialInterface<InterLJPotential> >;
    template class Inter2B2G3DFF< LJPotentialInterface<InterLJPotential> >;
}

static const RegisterMetaType<InterLJFF> r_interljff;
static const RegisterMetaType<InterGroupLJFF> r_intergroupljff;
