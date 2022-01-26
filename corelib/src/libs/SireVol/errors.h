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

#ifndef SIREVOL_ERRORS_H
#define SIREVOL_ERRORS_H

#include "SireError/exception.h"

SIRE_BEGIN_HEADER

namespace SireVol
{

/** This is the base class of all SireVol errors */
class SIREVOL_EXPORT sirevol_error : public SireError::exception
{
public:
    sirevol_error() : exception()
    {}

    sirevol_error(QString err, QString place = QString()) : exception(err,place)
    {}

    sirevol_error(const sirevol_error &other) : exception(other)
    {}

    ~sirevol_error() throw()
    {}

    static const char* typeName()
    {
        return "SireVol::sirevol_error";
    }
};

/** This exception is thrown when an attempt is made
    to interface or switch between incompatible spaces

    @author Christopher Woods
*/
class SIREVOL_EXPORT incompatible_space : public sirevol_error
{
public:
    incompatible_space() : sirevol_error()
    {}

    incompatible_space(QString err, QString place = QString())
              : sirevol_error(err,place)
    {}

    incompatible_space(const incompatible_space &other) : sirevol_error(other)
    {}

    ~incompatible_space() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return incompatible_space::typeName();
    }

    void throwSelf() const
    {
        throw incompatible_space(*this);
    }
};

}

Q_DECLARE_METATYPE(SireVol::incompatible_space)

SIRE_END_HEADER

#endif
