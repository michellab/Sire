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

#ifndef SIREERROR_GETMAGIC_H
#define SIREERROR_GETMAGIC_H

#include <QObject>

#include "SireError/exception.h"

SIRE_BEGIN_HEADER

namespace SireStream
{

using SireError::exception;

/** This exception is thrown whenever there is an error with the magic number of
    the binary data streaming protocol.

    @author Christopher Woods
*/
class SIRESTREAM_EXPORT magic_error : public SireError::exception
{
public:
    magic_error() : exception()
    {}

    magic_error(QString err, QString place = QString::null)
                  : exception(err,place)
    {}

    magic_error(MagicID wrongid,
                const RegisterMetaTypeBase &info,
                QString place=QString::null)
            : exception(QObject::tr(
                    "Magic error for \"%1\". Got %2, but expected %3.")
                        .arg(info.typeName()).arg(wrongid).arg(info.magicID()), place)
    {}

    magic_error(MagicID wrongid, MagicID rightid, const char *type_name,
                QString place=QString::null)
            : exception(QObject::tr(
                    "Magic error for \"%1\". Got %2, but expected %3.")
                        .arg(type_name).arg(wrongid).arg(rightid), place)
    {}

    magic_error(const magic_error &other) : exception(other)
    {}

    ~magic_error() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return magic_error::typeName();
    }
    
    void throwSelf() const
    {
        throw magic_error(*this);
    }
};

}

Q_DECLARE_METATYPE(SireStream::magic_error);

SIRE_END_HEADER

#endif
