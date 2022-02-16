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

#ifndef SIRESYSTEM_ERRORS_H
#define SIRESYSTEM_ERRORS_H

#include "SireError/exception.h"

SIRE_BEGIN_HEADER

/**
This file contains the exceptions that can be thrown by the SireSystem library.

@author Christopher Woods
*/

namespace SireSystem
{

/** This is the base class of all SireSystem errors */
class SIRESYSTEM_EXPORT siresystem_error : public SireError::exception
{
public:
    siresystem_error() : exception()
    {}

    siresystem_error(QString err, QString place = QString()) : exception(err,place)
    {}

    siresystem_error(const siresystem_error &other) : exception(other)
    {}

    ~siresystem_error() throw()
    {}

    static const char* typeName()
    {
        return "SireSystem::siresystem_error";
    }
};


/** This exception is thrown when a request is made of
    a non-existant monitor

    @author Christopher Woods
*/
class SIRESYSTEM_EXPORT missing_monitor : public siresystem_error
{
public:
    missing_monitor() : siresystem_error()
    {}

    missing_monitor(QString err, QString place = QString())
              : siresystem_error(err,place)
    {}

    missing_monitor(const missing_monitor &other) : siresystem_error(other)
    {}

    ~missing_monitor() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return missing_monitor::typeName();
    }

    void throwSelf() const
    {
        throw missing_monitor(*this);
    }
};

/** This exception is thrown when multiple monitors match an ID
    but only one monitor was requested

    @author Christopher Woods
*/
class SIRESYSTEM_EXPORT duplicate_monitor : public siresystem_error
{
public:
    duplicate_monitor() : siresystem_error()
    {}

    duplicate_monitor(QString err, QString place = QString())
              : siresystem_error(err,place)
    {}

    duplicate_monitor(const duplicate_monitor &other) : siresystem_error(other)
    {}

    ~duplicate_monitor() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return duplicate_monitor::typeName();
    }

    void throwSelf() const
    {
        throw duplicate_monitor(*this);
    }
};

/** This exception is thrown when a request is made of
    a non-existant system

    @author Christopher Woods
*/
class SIRESYSTEM_EXPORT missing_system : public siresystem_error
{
public:
    missing_system() : siresystem_error()
    {}

    missing_system(QString err, QString place = QString())
              : siresystem_error(err,place)
    {}

    missing_system(const missing_system &other) : siresystem_error(other)
    {}

    ~missing_system() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return missing_system::typeName();
    }

    void throwSelf() const
    {
        throw missing_system(*this);
    }
};

/** This exception is thrown when multiple systems match an ID
    but only one system was requested

    @author Christopher Woods
*/
class SIRESYSTEM_EXPORT duplicate_system : public siresystem_error
{
public:
    duplicate_system() : siresystem_error()
    {}

    duplicate_system(QString err, QString place = QString())
              : siresystem_error(err,place)
    {}

    duplicate_system(const duplicate_system &other) : siresystem_error(other)
    {}

    ~duplicate_system() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return duplicate_system::typeName();
    }

    void throwSelf() const
    {
        throw duplicate_system(*this);
    }
};

/** This exception is thrown when a violation of a constraint
    is detected

    @author Christopher Woods
*/
class SIRESYSTEM_EXPORT constraint_error : public siresystem_error
{
public:
    constraint_error() : siresystem_error()
    {}

    constraint_error(QString err, QString place = QString())
              : siresystem_error(err,place)
    {}

    constraint_error(const constraint_error &other) : siresystem_error(other)
    {}

    ~constraint_error() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return constraint_error::typeName();
    }

    void throwSelf() const
    {
        throw constraint_error(*this);
    }
};

}

Q_DECLARE_METATYPE(SireSystem::missing_monitor)
Q_DECLARE_METATYPE(SireSystem::duplicate_monitor)
Q_DECLARE_METATYPE(SireSystem::missing_system)
Q_DECLARE_METATYPE(SireSystem::duplicate_system)
Q_DECLARE_METATYPE(SireSystem::constraint_error)

SIRE_END_HEADER

#endif
