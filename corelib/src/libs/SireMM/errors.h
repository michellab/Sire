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

#ifndef SIREMM_ERRORS_H
#define SIREMM_ERRORS_H

#include "SireError/exception.h"

SIRE_BEGIN_HEADER

/**
This file contains the exceptions that can be thrown by the SireMM library.

@author Christopher Woods
*/

namespace SireMM
{

/** This is the base class of all SireMM errors */
class SIREMM_EXPORT siremm_error : public SireError::exception
{
public:
    siremm_error() : exception()
    {}

    siremm_error(QString err, QString place = QString()) : exception(err,place)
    {}

    siremm_error(const siremm_error &other) : exception(other)
    {}

    ~siremm_error() throw()
    {}

    static const char* typeName()
    {
        return "SireMM::siremm_error";
    }
};


/** This exception is thrown when a request is made of a non-existant bond

    @author Christopher Woods
*/
class SIREMM_EXPORT missing_bond : public siremm_error
{
public:
    missing_bond() : siremm_error()
    {}

    missing_bond(QString err, QString place = QString())
              : siremm_error(err,place)
    {}

    missing_bond(const missing_bond &other) : siremm_error(other)
    {}

    ~missing_bond() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return missing_bond::typeName();
    }

    void throwSelf() const
    {
        throw missing_bond(*this);
    }
};

/** This exception is thrown when a request is made of a non-existant angle

    @author Christopher Woods
*/
class SIREMM_EXPORT missing_angle : public siremm_error
{
public:
    missing_angle() : siremm_error()
    {}

    missing_angle(QString err, QString place = QString())
              : siremm_error(err,place)
    {}

    missing_angle(const missing_angle &other) : siremm_error(other)
    {}

    ~missing_angle() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return missing_angle::typeName();
    }

    void throwSelf() const
    {
        throw missing_angle(*this);
    }
};

/** This exception is thrown when a request is made of a non-existant dihedral

    @author Christopher Woods
*/
class SIREMM_EXPORT missing_dihedral : public siremm_error
{
public:
    missing_dihedral() : siremm_error()
    {}

    missing_dihedral(QString err, QString place = QString())
              : siremm_error(err,place)
    {}

    missing_dihedral(const missing_dihedral &other) : siremm_error(other)
    {}

    ~missing_dihedral() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return missing_dihedral::typeName();
    }

    void throwSelf() const
    {
        throw missing_dihedral(*this);
    }
};

}

Q_DECLARE_METATYPE(SireMM::missing_bond)
Q_DECLARE_METATYPE(SireMM::missing_angle)
Q_DECLARE_METATYPE(SireMM::missing_dihedral)

SIRE_END_HEADER

#endif
