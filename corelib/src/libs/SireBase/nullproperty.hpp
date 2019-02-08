/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2016  Christopher Woods
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

#ifndef SIREBASE_NULLPROPERTY_HPP
#define SIREBASE_NULLPROPERTY_HPP

#include "property.h"
#include "SireStream/datastream.h"

SIRE_BEGIN_HEADER

namespace SireBase
{
template<class T>
class NullProp;
}

template<class T>
QDataStream& operator<<(QDataStream&, const SireBase::NullProp<T>&);
template<class T>
QDataStream& operator>>(QDataStream&, SireBase::NullProp<T>&);

namespace SireBase
{

/** This template simplifies the creation of a 'NullProperty' for each
    PropPtr<T> derived type
    
    @author Christopher Woods
*/
template<class T>
class SIREBASE_EXPORT NullProp : public ConcreteProperty< NullProp<T>, T >
{

friend SIREBASE_EXPORT QDataStream& ::operator<<<>(QDataStream&, const NullProp<T>&);
friend SIREBASE_EXPORT QDataStream& ::operator>><>(QDataStream&, NullProp<T>&);

public:
    NullProp() : ConcreteProperty< NullProp<T>, T >()
    {}
    
    NullProp(const NullProp<T> &other) : ConcreteProperty< NullProp<T>, T >(other)
    {}
    
    ~NullProp();
    
    NullProp<T>& operator=(const NullProp<T>&)
    {
        return *this;
    }
    
    bool operator==(const NullProp<T>&) const
    {
        return true;
    }
    
    bool operator!=(const NullProp<T>&) const
    {
        return false;
    }
    
    static const char* typeName()
    {
        return QMetaType::typeName( qMetaTypeId< NullProp<T> >() );
    }
    
    QString toString() const
    {
        return QObject::tr( "%1::null" ).arg(T::typeName());
    }

private:
    static const Sire::RegisterMetaType< NullProp<T> > r_null;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

template<class T>
SIRE_OUTOFLINE_TEMPLATE
NullProp<T>::~NullProp()
{}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
const Sire::RegisterMetaType< NullProp<T> > NullProp<T>::r_null;

#endif

}

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

template<class T>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator<<(QDataStream &ds, const SireBase::NullProp<T>&)
{
    SireStream::writeHeader(ds, SireBase::NullProp<T>::r_null, 1);
    return ds;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator>>(QDataStream &ds, SireBase::NullProp<T>&)
{
    SireStream::VersionID v = SireStream::readHeader(ds, SireBase::NullProp<T>::r_null);
    
    if (v != 1)
        throw SireStream::version_error(v, "1", SireBase::NullProp<T>::r_null, CODELOC);
    
    return ds;
}

#endif

SIRE_END_HEADER

#endif
