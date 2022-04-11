/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2010  Christopher Woods
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

#ifndef SIREBASE_PACKEDARRAYS_H
#define SIREBASE_PACKEDARRAYS_H

#include "packedarray2d.hpp"
#include "propertylist.h"
#include "boost/throw_exception.hpp"

SIRE_BEGIN_HEADER

SIRE_EXPOSE_ALIAS(SireBase::PackedArray2D<double>,
                  SireBase::PackedArray2D_double_)

SIRE_EXPOSE_ALIAS(SireBase::detail::PackedArray2D_Array<double>,
                  SireBase::PackedArray2D_double_Array)

SIRE_EXPOSE_ALIAS(SireBase::PackedArray2D<int>,
                  SireBase::PackedArray2D_int_);

SIRE_EXPOSE_ALIAS(SireBase::detail::PackedArray2D_Array<int>,
                  SireBase::PackedArray2D_int_Array)

SIRE_EXPOSE_ALIAS(SireBase::PackedArray2D<QString>,
                  SireBase::PackedArray2D_QString_)

SIRE_EXPOSE_ALIAS(SireBase::detail::PackedArray2D_Array<QString>,
                  SireBase::PackedArray2D_QString_Array)

SIRE_EXPOSE_ALIAS(SireBase::PackedArray2D<QVariant>,
                  SireBase::PackedArray2D_QVariant_)

SIRE_EXPOSE_ALIAS(SireBase::detail::PackedArray2D_Array<QVariant>,
                  SireBase::PackedArray2D_QVariant_Array)

SIRE_EXPOSE_ALIAS(SireBase::PackedArray2D<SireBase::PropertyList>,
                  SireBase::PackedArray2D_PropertyList)

SIRE_EXPOSE_ALIAS(SireBase::detail::PackedArray2D_Array<SireBase::PropertyList>,
                  SireBase::PackedArray2D_PropertyList_Array)

SIRE_EXPOSE_ALIAS(SireBase::PackedArray2D<SireBase::DoubleArrayProperty>,
                  SireBase::PackedArray2D_DoubleArrayProperty)

SIRE_EXPOSE_ALIAS(SireBase::detail::PackedArray2D_Array<SireBase::DoubleArrayProperty>,
                  SireBase::PackedArray2D_DoubleArrayProperty_Array)

SIRE_EXPOSE_ALIAS(SireBase::PackedArray2D<SireBase::IntegerArrayProperty>,
                  SireBase::PackedArray2D_IntegerArrayProperty)

SIRE_EXPOSE_ALIAS(SireBase::detail::PackedArray2D_Array<SireBase::IntegerArrayProperty>,
                  SireBase::PackedArray2D_IntegerArrayProperty_Array)

SIRE_EXPOSE_ALIAS(SireBase::PackedArray2D<SireBase::StringArrayProperty>,
                  SireBase::PackedArray2D_StringArrayProperty)

SIRE_EXPOSE_ALIAS(SireBase::detail::PackedArray2D_Array<SireBase::StringArrayProperty>,
                  SireBase::PackedArray2D_StringArrayProperty_Array)


#ifdef SIRE_INSTANTIATE_TEMPLATES
template class
SireBase::PackedArray2D<double>;

template class
SireBase::detail::PackedArray2D_Array<double>;

template class
SireBase::PackedArray2D<int>;

template class
SireBase::detail::PackedArray2D_Array<int>;

template class
SireBase::PackedArray2D<QString>;

template class
SireBase::detail::PackedArray2D_Array<QString>;

template class
SireBase::PackedArray2D<QVariant>;

template class
SireBase::detail::PackedArray2D_Array<QVariant>;

template class
SireBase::PackedArray2D<SireBase::PropertyList>;

template class
SireBase::detail::PackedArray2D_Array<SireBase::PropertyList>;

template class
SireBase::PackedArray2D<SireBase::DoubleArrayProperty>;

template class
SireBase::detail::PackedArray2D_Array<SireBase::DoubleArrayProperty>;

template class
SireBase::PackedArray2D<SireBase::IntegerArrayProperty>;

template class
SireBase::detail::PackedArray2D_Array<SireBase::IntegerArrayProperty>;

template class
SireBase::PackedArray2D<SireBase::StringArrayProperty>;

template class
SireBase::detail::PackedArray2D_Array<SireBase::StringArrayProperty>;

#endif

SIRE_END_HEADER

#endif
