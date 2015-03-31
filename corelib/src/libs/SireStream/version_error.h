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

#ifndef SIRESTREAM_VERSION_ERROR_H
#define SIRESTREAM_VERSION_ERROR_H

#include "SireError/version_error.h"

#include "versionid.h"

SIRE_BEGIN_HEADER

namespace SireStream
{

/** This exception is thrown whenever there is an error with the version number of
    the binary data streaming protocol.

    @author Christopher Woods
*/
class SIRESTREAM_EXPORT version_error : public SireError::version_error
{
public:
    version_error() : SireError::version_error()
    {}

    version_error(QString err, QString place = QString::null)
                  : SireError::version_error(err,place)
    {}

    version_error(VersionID wrongid, QString supported_ids,
                  const RegisterMetaTypeBase &info,
                  QString place=QString::null)
            : SireError::version_error(QObject::tr(
                    "Incompatible version for \"%1\". Got %2, but can only support [ %3 ].")
                        .arg(info.typeName()).arg(wrongid).arg(supported_ids), place)
    {}

    version_error(VersionID wrongid, QString supported_ids,
                  const char *type_name, QString place=QString::null)
            : SireError::version_error(QObject::tr(
                    "Incompatible version for \"%1\". Got %2, but can only support [ %3 ].")
                        .arg(type_name).arg(wrongid).arg(supported_ids), place)
    {}

    version_error(const version_error &other) : SireError::version_error(other)
    {}

    ~version_error() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return SireStream::version_error::typeName();
    }
    
    void throwSelf() const
    {
        throw SireStream::version_error(*this);
    }
};

}

SIRE_END_HEADER

Q_DECLARE_METATYPE(SireStream::version_error)

#endif
