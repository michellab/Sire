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

#ifndef SIREFF_ERRORS_H
#define SIREFF_ERRORS_H

#include "SireError/exception.h"

SIRE_BEGIN_HEADER

/**
This file contains the exceptions that can be thrown by the SireFF library.

@author Christopher Woods
*/

namespace SireFF
{

/** This is the base class of all SireFF errors */
class SIREFF_EXPORT sireff_error : public SireError::exception
{
public:
    sireff_error() : exception()
    {}

    sireff_error(QString err, QString place = QString::null) : exception(err,place)
    {}

    sireff_error(const sireff_error &other) : exception(other)
    {}

    ~sireff_error() throw()
    {}

    static const char* typeName()
    {
        return "SireMol::sireff_error";
    }
};


/** This exception is thrown when a request is made of a non-existant
    forcefield component

    @author Christopher Woods
*/
class SIREFF_EXPORT missing_component : public sireff_error
{
public:
    missing_component() : sireff_error()
    {}

    missing_component(QString err, QString place = QString::null)
              : sireff_error(err,place)
    {}

    missing_component(const missing_component &other) : sireff_error(other)
    {}

    ~missing_component() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return missing_component::typeName();
    }
    
    void throwSelf() const
    {
        throw missing_component(*this);
    }
};

/** This exception is thrown when a request is made of a non-existant
    forcefield function

    @author Christopher Woods
*/
class SIREFF_EXPORT missing_function : public sireff_error
{
public:
    missing_function() : sireff_error()
    {}

    missing_function(QString err, QString place = QString::null)
              : sireff_error(err,place)
    {}

    missing_function(const missing_function &other) : sireff_error(other)
    {}

    ~missing_function() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return missing_function::typeName();
    }
    
    void throwSelf() const
    {
        throw missing_function(*this);
    }
};

/** This exception is thrown when a request is made of a non-existant
    forcefield

    @author Christopher Woods
*/
class SIREFF_EXPORT missing_forcefield : public sireff_error
{
public:
    missing_forcefield() : sireff_error()
    {}

    missing_forcefield(QString err, QString place = QString::null)
              : sireff_error(err,place)
    {}

    missing_forcefield(const missing_forcefield &other) : sireff_error(other)
    {}

    ~missing_forcefield() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return missing_forcefield::typeName();
    }
    
    void throwSelf() const
    {
        throw missing_forcefield(*this);
    }
};

/** This exception is thrown when a request is made of an unavailable
    forcefield parameter

    @author Christopher Woods
*/
class SIREFF_EXPORT missing_parameter : public sireff_error
{
public:
    missing_parameter() : sireff_error()
    {}

    missing_parameter(QString err, QString place = QString::null)
              : sireff_error(err,place)
    {}

    missing_parameter(const missing_parameter &other) : sireff_error(other)
    {}

    ~missing_parameter() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return missing_parameter::typeName();
    }
    
    void throwSelf() const
    {
        throw missing_parameter(*this);
    }
};

/** This exception is thrown when a request is made of an unavailable
    derivative of a function (either because the derivative has
    not been programmed, or because it is not mathematically available)

    @author Christopher Woods
*/
class SIREFF_EXPORT missing_derivative : public sireff_error
{
public:
    missing_derivative() : sireff_error()
    {}

    missing_derivative(QString err, QString place = QString::null)
              : sireff_error(err,place)
    {}

    missing_derivative(const missing_derivative &other) : sireff_error(other)
    {}

    ~missing_derivative() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return missing_derivative::typeName();
    }
    
    void throwSelf() const
    {
        throw missing_derivative(*this);
    }
};

/** This exception is thrown when a forcefield is added when
    one already exists!

    @author Christopher Woods
*/
class SIREFF_EXPORT duplicate_forcefield : public sireff_error
{
public:
    duplicate_forcefield() : sireff_error()
    {}

    duplicate_forcefield(QString err, QString place = QString::null)
              : sireff_error(err,place)
    {}

    duplicate_forcefield(const duplicate_forcefield &other) : sireff_error(other)
    {}

    ~duplicate_forcefield() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return duplicate_forcefield::typeName();
    }
    
    void throwSelf() const
    {
        throw duplicate_forcefield(*this);
    }
};

/** This exception is thrown when a forcefield function
    is added when one already exists!

    @author Christopher Woods
*/
class SIREFF_EXPORT duplicate_component : public sireff_error
{
public:
    duplicate_component() : sireff_error()
    {}

    duplicate_component(QString err, QString place = QString::null)
              : sireff_error(err,place)
    {}

    duplicate_component(const duplicate_component &other) : sireff_error(other)
    {}

    ~duplicate_component() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return duplicate_component::typeName();
    }
    
    void throwSelf() const
    {
        throw duplicate_component(*this);
    }
};

/** This exception is thrown when a component is added when
    one already exists!

    @author Christopher Woods
*/
class SIREFF_EXPORT duplicate_function : public sireff_error
{
public:
    duplicate_function() : sireff_error()
    {}

    duplicate_function(QString err, QString place = QString::null)
              : sireff_error(err,place)
    {}

    duplicate_function(const duplicate_function &other) : sireff_error(other)
    {}

    ~duplicate_function() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return duplicate_function::typeName();
    }
    
    void throwSelf() const
    {
        throw duplicate_function(*this);
    }
};

/** This exception is thrown when an invalid forcefield group
    is selected or used

    @author Christopher Woods
*/
class SIREFF_EXPORT invalid_group : public sireff_error
{
public:
    invalid_group() : sireff_error()
    {}

    invalid_group(QString err, QString place = QString::null)
              : sireff_error(err,place)
    {}

    invalid_group(const invalid_group &other) : sireff_error(other)
    {}

    ~invalid_group() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return invalid_group::typeName();
    }
    
    void throwSelf() const
    {
        throw invalid_group(*this);
    }
};

}

Q_DECLARE_METATYPE(SireFF::missing_component)
Q_DECLARE_METATYPE(SireFF::missing_function)
Q_DECLARE_METATYPE(SireFF::missing_forcefield)
Q_DECLARE_METATYPE(SireFF::missing_derivative)
Q_DECLARE_METATYPE(SireFF::missing_parameter)
Q_DECLARE_METATYPE(SireFF::duplicate_component)
Q_DECLARE_METATYPE(SireFF::duplicate_function)
Q_DECLARE_METATYPE(SireFF::duplicate_forcefield)
Q_DECLARE_METATYPE(SireFF::invalid_group)

SIRE_END_HEADER

#endif
