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

#ifndef SIREIO_ERRORS_H
#define SIREIO_ERRORS_H

#include "SireError/exception.h"

SIRE_BEGIN_HEADER

/** This file contains the exceptions that can be thrown by the SireIO library.

    @author Christopher Woods
*/

namespace SireIO
{

/** This is the base class of all SireMM errors */
class SIREIO_EXPORT sireio_error : public SireError::exception
{
public:
    sireio_error() : exception()
    {}

    sireio_error(QString err, QString place = QString()) : exception(err,place)
    {}

    sireio_error(const sireio_error &other) : exception(other)
    {}

    ~sireio_error() throw()
    {}

    static const char* typeName()
    {
        return "SireIO::sireio_error";
    }
};

/** This exception is thrown when there is a non-recoverable error
    while parsing a file

    @author Christopher Woods
*/
class SIREIO_EXPORT parse_error : public sireio_error
{
public:
    parse_error() : sireio_error()
    {}

    parse_error(QString err, QString place = QString())
              : sireio_error(err,place)
    {}

    parse_error(const parse_error &other) : sireio_error(other)
    {}

    ~parse_error() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return parse_error::typeName();
    }

    void throwSelf() const
    {
        throw parse_error(*this);
    }
};

}

Q_DECLARE_METATYPE(SireIO::parse_error)

SIRE_END_HEADER

#endif
