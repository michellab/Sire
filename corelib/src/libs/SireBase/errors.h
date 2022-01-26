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

#ifndef SIREBASE_ERRORS_H
#define SIREBASE_ERRORS_H

#include "SireError/exception.h"

SIRE_BEGIN_HEADER

namespace SireBase
{

/** This is the base class of all SireBase errors */
class SIREBASE_EXPORT sirebase_error : public SireError::exception
{
public:
    sirebase_error() : exception()
    {}

    sirebase_error(QString err, QString place = QString()) : exception(err,place)
    {}

    sirebase_error(const sirebase_error &other) : exception(other)
    {}

    ~sirebase_error() throw()
    {}

    static const char* typeName()
    {
        return "SireBase::sirebase_error";
    }
};

/** This exception is thrown when a request is made of a non-existant property

    @author Christopher Woods
*/
class SIREBASE_EXPORT missing_property : public sirebase_error
{
public:
    missing_property() : sirebase_error()
    {}

    missing_property(QString err, QString place = QString())
              : sirebase_error(err,place)
    {}

    missing_property(const missing_property &other) : sirebase_error(other)
    {}

    ~missing_property() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return missing_property::typeName();
    }

    void throwSelf() const
    {
        throw missing_property(*this);
    }
};


/** This exception is thrown when a request is made to duplicate a
    property when this would be inappropriate

    @author Christopher Woods
*/
class SIREBASE_EXPORT duplicate_property : public sirebase_error
{
public:
    duplicate_property() : sirebase_error()
    {}

    duplicate_property(QString err, QString place = QString())
              : sirebase_error(err,place)
    {}

    duplicate_property(const duplicate_property &other) : sirebase_error(other)
    {}

    ~duplicate_property() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return duplicate_property::typeName();
    }

    void throwSelf() const
    {
        throw duplicate_property(*this);
    }
};

}

Q_DECLARE_METATYPE(SireBase::missing_property)
Q_DECLARE_METATYPE(SireBase::duplicate_property)

SIRE_END_HEADER

#endif
