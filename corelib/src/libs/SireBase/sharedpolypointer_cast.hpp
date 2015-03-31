/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2006  Christopher Woods
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

#ifndef SIREBASE_SHAREDPOLYPOINTERCAST_HPP

#include "sharedpolypointer.hpp"

#include "SireError/errors.h"

SIRE_BEGIN_HEADER

namespace SireBase
{

/** Cast the pointer 'ptr' from type 'T' to an object of type 'S'

    \throw SireError::nullptr_error
    \throw SireError::invalid_cast
*/
template<class S, class T>
const S& sharedpolypointer_cast(const SharedPolyPointer<T> &ptr)
{
    if (not ptr)
        throw SireError::nullptr_error( QObject::tr("Cannot cast a null pointer!"),
                                        CODELOC );

    const S *sptr = dynamic_cast<const S*>(ptr.constData());

    if (not sptr)
        throw SireError::invalid_cast( QObject::tr(
                      "Cannot cast a SharedPolyPointer of type \"%1\" to "
                      "type \"%2\".")
                          .arg( SharedPolyPointerHelper<T>::what(*ptr),
                                SharedPolyPointerHelper<S>::typeName() ), CODELOC );

    return *sptr;
}

}

SIRE_END_HEADER

#endif
