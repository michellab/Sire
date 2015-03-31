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

#ifndef SIREID_ID_H
#define SIREID_ID_H

#include "sireglobal.h"

#include <typeinfo>

#include <QString>
#include <QMetaType>

SIRE_BEGIN_HEADER

namespace SireID
{
class ID;
}

namespace SireID
{

/** This is the base class of all ID objects. An ID object
    is an object that is used to identify another object.
    
    @author Christopher Woods
*/
class SIREID_EXPORT ID
{
public:
    typedef ID ROOT;

    ID();
    ID(const ID&);
    
    virtual ~ID();
    
    static const char* typeName()
    {
        return "SireID::ID";
    }
    
    /** Return a clone of this ID object. You are responsible
        for managing the returned object. */
    virtual ID* clone() const=0;
    
    /** Return the type name of this ID object. */
    virtual const char* what() const=0;
    
    /** Return a hash for this ID object - this allows
        this object to be used as a key in a dictionary */
    virtual uint hash() const=0;
 
    /** Return a string representation of this ID */
    virtual QString toString() const=0;
 
    /** Return whether or not this ID is null */
    virtual bool isNull() const=0;
                            
    /** Comparison operator */
    virtual bool operator==(const ID &other) const=0;
    
    /** Comparison operator */
    virtual bool operator!=(const ID &other) const
    {
        return not this->operator==(other);
    }
    
    /** Use this function in your comparison operators, e.g.
    
        <code>
        bool Foo::operator==(const ID &other) const
        {
            return ID::compare(*this, other);
        }
        </code>
    */
    template<class T>
    static bool compare(const T &obj, const ID &other)
    {
        if ( typeid(obj).name() == typeid(other).name() )
        {
            const T *other_t = dynamic_cast<const T*>(&other);
            return (other_t != 0) and (obj == *other_t);
        }
        else
            return false;
    }
    
    template<class T>
    bool isA() const
    {
        return dynamic_cast<const T*>(this) != 0;
    }
    
    template<class T>
    const T& asA() const
    {
        return dynamic_cast<const T&>(*this);
    }

protected:
    ID& operator=(const ID&)
    {
        return *this;
    }
};

}

SIRE_EXPOSE_CLASS( SireID::ID )

SIRE_END_HEADER

#endif
