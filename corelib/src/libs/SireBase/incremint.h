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

#ifndef SIREBASE_INCREMINT_H
#define SIREBASE_INCREMINT_H

#include "sireglobal.h"

#if QT_VERSION >= 0x040400
  #include <QAtomicInt>
#elif QT_VERSION >= 0x040100
  #include <QAtomic>
#else
  #error You need at least Qt Version 4.1
#endif

#include <boost/throw_exception.hpp>

SIRE_BEGIN_HEADER

namespace SireBase
{

/** This is a simple class that provides a thread-safe
    incrementable integer. This can be used, for example,
    to give a unique version of ID number to objects.
    
    @author Christopher Woods
*/
class Incremint
{
public:
    Incremint(int value=0);
    Incremint(const Incremint &other);
    
    ~Incremint();
    
    int increment();

private:
    /** The volatile integer that is being incremented */
    QAtomicInt atomic_int;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Constructor */
inline Incremint::Incremint(int value) : atomic_int(value)
{}

/** Copy constructor */
inline Incremint::Incremint(const Incremint &other) 
                 : atomic_int(other.atomic_int)
{}

/** Destructor */
inline Incremint::~Incremint()
{}

/** Increment the Incremint and return its new value */
inline int Incremint::increment()
{
    int value = atomic_int.fetchAndAddRelaxed(1);
    return value+1;
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

} // end of namespace SireBase

SIRE_EXPOSE_CLASS( SireBase::Incremint )

SIRE_END_HEADER

#endif
