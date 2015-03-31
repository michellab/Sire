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

#ifndef SIREID_NUMBER_H
#define SIREID_NUMBER_H

#include "sireglobal.h"

#include <limits>

SIRE_BEGIN_HEADER

namespace SireID
{
class Number;
}

QDataStream& operator<<(QDataStream&, const SireID::Number&);
QDataStream& operator>>(QDataStream&, SireID::Number&);

namespace SireID
{

/** This is the base class of all Number ID objects. A Number
    is used to provide an object with an identifying number.
    This could be the number of a residue in a molecule, a 
    user-supplied number of a CutGroup in a molecule, or
    perhaps the automatic unique ID numbers of molecules, 
    forcefields or molecule groups that are assigned by the
    program. The key point of a Number is to provide an ID
    that can be queried and compared rapidly, and that 
    does not change as the object is moved between different
    containers. Generally an object should keep its number
    throughout its lifetime.
    
    @author Christopher Woods
*/
class SIREID_EXPORT Number
{

friend QDataStream& ::operator<<(QDataStream&, const Number&);
friend QDataStream& ::operator>>(QDataStream&, Number&);

public:
    ~Number();
    
    operator qint32() const;
    
    static qint32 null();
    
    bool isNull() const;

    uint hash() const;
    
    qint32 value() const;
    
protected:
    Number(qint32 num=0);
    
    Number(const Number &other);
    
    /** The actual number */
    qint32 _num;
};

/** Return a hash of this Number */
inline uint qHash(const Number &number)
{
    return number.hash();
}

}

SIRE_EXPOSE_CLASS( SireID::Number )

SIRE_END_HEADER

#endif
