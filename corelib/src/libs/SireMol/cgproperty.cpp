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

#include "cgproperty.hpp"

using namespace SireMol;

namespace SireMol
{
    template class CGProperty<QString>;
    template class CGProperty<qint64>;
    template class CGProperty<double>;
    template class CGProperty<QVariant>;
    template class CGProperty<SireBase::PropertyPtr>;
}

static const RegisterMetaType<CGStringProperty> r_cgstring;
static const RegisterMetaType<CGIntProperty> r_cgint;
static const RegisterMetaType<CGFloatProperty> r_cgfloat;
static const RegisterMetaType<CGVariantProperty> r_cgvariant;
static const RegisterMetaType<CGPropertyProperty> r_cgproprop;
