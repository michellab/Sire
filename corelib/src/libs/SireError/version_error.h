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

#ifndef SIREERROR_VERSION_ERROR_H
#define SIREERROR_VERSION_ERROR_H

#include "SireError/exception.h"

SIRE_BEGIN_HEADER

namespace SireError
{

/** This exception is thrown whenever there is an error with a
    version number

    @author Christopher Woods
*/
class version_error : public exception
{
public:
    version_error() : exception()
    {}

    version_error(QString err, QString place = QString())
                  : exception(err,place)
    {}

    version_error(const version_error &other) : exception(other)
    {}

    ~version_error() throw()
    {}

    static const char* typeName();

    const char* what() const throw();

    void throwSelf() const
    {
        throw version_error(*this);
    }
};

}

Q_DECLARE_METATYPE(SireError::version_error)

SIRE_END_HEADER

#endif
