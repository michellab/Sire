/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2008  Christopher Woods
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

#include "majorminorversion.h"

#include "SireStream/datastream.h"

using namespace SireBase;
using namespace SireBase::detail;
using namespace SireStream;

static const RegisterMetaType<Version> r_version(NO_ROOT);

QDataStream SIREBASE_EXPORT &operator<<(QDataStream &ds, const Version &version)
{
    writeHeader(ds, r_version, 1);
    ds << version.maj << version.min;
    return ds;
}

QDataStream SIREBASE_EXPORT &operator>>(QDataStream &ds, Version &version)
{
    VersionID v = readHeader(ds, r_version);
    
    if (v == 1)
    {
        ds >> version.maj >> version.min;
    }
    else
        throw version_error(v, "1", r_version, CODELOC);
    
    return ds;
}

/** Constructor */
Version::Version(quint64 major, quint64 minor)
        : maj(major), min(minor)
{}

/** Copy constructor */
Version::Version(const Version &other)
        : maj(other.maj), min(other.min)
{}

/** Destructor */
Version::~Version()
{}

const char* Version::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Version>() );
}

/** Return a string representation of this version number */
QString Version::toString() const
{
    return QString("%1.%2").arg(maj).arg(min);
}

boost::shared_ptr<MajorMinorVersionData> 
MajorMinorVersion::shared_null( new MajorMinorVersionData() );

/** Null constructor */
MajorMinorVersion::MajorMinorVersion() 
                  : d(shared_null), v(0,0)
{}

/** Construct from a raw data object - this should only be called by 
    the registry function */
MajorMinorVersion::MajorMinorVersion(
                        const boost::shared_ptr<MajorMinorVersionData> &ptr)
                  : d(ptr)
{
    QMutexLocker lkr( &(d->version_mutex) );
    
    v = Version(d->last_major_version, d->last_minor_version);
}

/** Construct the object for a specific version */
MajorMinorVersion::MajorMinorVersion(quint64 vmaj, quint64 vmin)
                  : d( new MajorMinorVersionData(vmaj, vmin) ),
                    v(vmaj, vmin)
{}

/** Copy constructor */
MajorMinorVersion::MajorMinorVersion(const MajorMinorVersion &other)
                  : d(other.d), v(other.v)
{}

/** Destructor */
MajorMinorVersion::~MajorMinorVersion()
{}

const char* MajorMinorVersion::typeName()
{
    return QMetaType::typeName( qMetaTypeId<MajorMinorVersion>() );
}

/** Increment the major version number - this resets the 
    minor version number to 0 */
void MajorMinorVersion::incrementMajor()
{
    QMutexLocker lkr( &(d->version_mutex) );
    
    ++(d->last_major_version);
    d->last_minor_version = 0;
    
    v = Version(d->last_major_version, 0);
}

/** Increment the minor version number */
void MajorMinorVersion::incrementMinor()
{
    QMutexLocker lkr( &(d->version_mutex) );
    
    ++(d->last_minor_version);
    
    v = Version(v.majorVersion(), d->last_minor_version);
}
