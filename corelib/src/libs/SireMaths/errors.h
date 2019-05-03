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

#ifndef SIREMATHS_ERRORS_H
#define SIREMATHS_ERRORS_H

#include "SireError/exception.h"

SIRE_BEGIN_HEADER

namespace SireMaths
{

/** This is the base class of all SireMaths errors */
class SIREMATHS_EXPORT siremaths_error : public SireError::exception
{
public:
    siremaths_error() : exception()
    {}
    
    siremaths_error(QString err, QString place = QString::null) 
                  : exception(err,place)
    {}
    
    siremaths_error(const siremaths_error &other) : exception(other)
    {}
    
    ~siremaths_error() throw()
    {}
    
    static const char* typeName()
    {
        return "SireMaths::siremaths_error";
    }
};

/** This class represents a general maths error */
class SIREMATHS_EXPORT math_error : public siremaths_error
{
public:
    math_error() : siremaths_error()
    {}
    
    math_error(QString err, QString place = QString::null) 
              : siremaths_error(err,place)
    {}
    
    math_error(const math_error &other) : siremaths_error(other)
    {}
    
    ~math_error() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return math_error::typeName();
    }
    
    void throwSelf() const
    {
        throw math_error(*this);
    }
};

/** This class represents a domain error */
class SIREMATHS_EXPORT domain_error : public siremaths_error
{
public:
    domain_error() : siremaths_error()
    {}
    
    domain_error(QString err, QString place = QString::null) 
              : siremaths_error(err,place)
    {}
    
    domain_error(const domain_error &other) : siremaths_error(other)
    {}
    
    ~domain_error() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return domain_error::typeName();
    }
    
    void throwSelf() const
    {
        throw domain_error(*this);
    }
};

}

Q_DECLARE_METATYPE(SireMaths::math_error)
Q_DECLARE_METATYPE(SireMaths::domain_error)

SIRE_END_HEADER

#endif
