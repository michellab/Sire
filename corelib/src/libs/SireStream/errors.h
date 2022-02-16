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

#ifndef SIRESTREAM_ERRORS_H
#define SIRESTREAM_ERRORS_H

#include "SireError/exception.h"

SIRE_BEGIN_HEADER

namespace SireStream
{

/** This exception is thrown when corrupted data is detected

    @author Christopher Woods
*/
class SIRESTREAM_EXPORT corrupted_data : public SireError::exception
{
public:
    corrupted_data() : SireError::exception()
    {}

    corrupted_data(QString err, QString place = QString())
              : SireError::exception(err,place)
    {}

    corrupted_data(const corrupted_data &other) : SireError::exception(other)
    {}

    ~corrupted_data() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return corrupted_data::typeName();
    }

    void throwSelf() const
    {
        throw corrupted_data(*this);
    }
};

}

Q_DECLARE_METATYPE(SireStream::corrupted_data)

SIRE_END_HEADER

#endif
