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

#ifndef SIREMOVE_ERRORS_H
#define SIREMOVE_ERRORS_H

#include "SireError/exception.h"

SIRE_BEGIN_HEADER

/**
This file contains the exceptions that can be thrown by the SireMove library.

@author Christopher Woods
*/

namespace SireMove
{

/** This is the base class of all SireMove errors */
class SIREMOVE_EXPORT siremove_error : public SireError::exception
{
public:
    siremove_error() : exception()
    {}

    siremove_error(QString err, QString place = QString::null) : exception(err,place)
    {}

    siremove_error(const siremove_error &other) : exception(other)
    {}

    ~siremove_error() throw()
    {}

    static const char* typeName()
    {
        return "SireMove::siremove_error";
    }
};

/** This exception is thrown when a z-matrix error is detected

    @author Christopher Woods
*/
class SIREMOVE_EXPORT zmatrix_error : public siremove_error
{
public:
    zmatrix_error() : siremove_error()
    {}

    zmatrix_error(QString err, QString place = QString::null)
              : siremove_error(err,place)
    {}

    zmatrix_error(const zmatrix_error &other) : siremove_error(other)
    {}

    ~zmatrix_error() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return zmatrix_error::typeName();
    }
    
    void throwSelf() const
    {
        throw zmatrix_error(*this);
    }
};

}

Q_DECLARE_METATYPE(SireMove::zmatrix_error)

SIRE_END_HEADER

#endif
