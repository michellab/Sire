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

#ifndef SIRESTREAM_STREAMHELPER_HPP
#define SIRESTREAM_STREAMHELPER_HPP

#include <QByteArray>
#include <QLocale>
#include <QDateTime>
#include <QList>
#include <QStringList>

#include <boost/tuple/tuple.hpp>
#include <boost/shared_ptr.hpp>

#include "md5sum.h"

#include "sireglobal.h"

#include <QDebug>

SIRE_BEGIN_HEADER

namespace SireStream
{
class FileHeader;
}

QDataStream& operator<<(QDataStream&, const SireStream::FileHeader&);
QDataStream& operator>>(QDataStream&, SireStream::FileHeader&);

namespace SireStream
{

namespace detail
{
QByteArray streamDataSave( const void *object, const char *type_name );

QByteArray streamDataSave( 
                 const QList< boost::tuple<const void*,const char*> > &objects ); 

QByteArray streamDataSave(
                 const QList< boost::tuple<boost::shared_ptr<void>,QString> > &objects );

void streamDataSave( const void *object, const char *type_name, 
                     const QString &filename );

void streamDataSave(
                const QList< boost::tuple<const void*,const char*> > &objects,
                const QString &filename );

void streamDataSave(
                const QList< boost::tuple<boost::shared_ptr<void>,QString> > &objects,
                const QString &filename );

void throwStreamDataInvalidCast(const QString &load_type, 
                                const QString &cast_type);

struct void_deleter
{
public:
    void_deleter(int type_id) : id(type_id)
    {}
    
    ~void_deleter()
    {}
    
    void operator()(void const *ptr) const
    {
        QMetaType::destroy(id, const_cast<void*>(ptr));
    }

private:
    int id;
};

}

QList< boost::tuple<boost::shared_ptr<void>,QString> > load(const QByteArray &data);
QList< boost::tuple<boost::shared_ptr<void>,QString> > load(const QString &filename);

/** This class provides metadata about the binary representation
    of an object. This is to allow the owner of the data to identify
    it as belonging to themselves, to provide information about 
    what data is contained, when it was created and on what,
    and to provide some information about how much space the 
    data may require to load. The aim is to allow the provenance
    (well, at least the origin!) of a file to be determined.
    
    @author Christopher Woods
*/
class SIRESTREAM_EXPORT FileHeader
{

friend QDataStream& ::operator<<(QDataStream&, const FileHeader&);
friend QDataStream& ::operator>>(QDataStream&, FileHeader&);

friend QByteArray detail::streamDataSave( 
                        const QList< boost::tuple<const void*,const char*> >& );

friend QByteArray detail::streamDataSave( 
                        const QList< boost::tuple<boost::shared_ptr<void>,QString> >& );

public:
    FileHeader();
    FileHeader(const FileHeader &other);
    
    ~FileHeader();
    
    FileHeader& operator=(const FileHeader &other);
    
    QString toString() const;
    
    const QString& createdBy() const;
    const QDateTime& createdWhen() const;
    const QString& createdWhere() const;
    
    const QString& systemInfo() const;
    
    quint32 requiredMemory() const;
    
    double compressionRatio() const;
    
    QString dataType() const;
    const QStringList& dataTypes() const;
    
    QStringList requiredLibraries() const;

    bool requireLibrary(const QString &library) const;
    quint32 requiredVersion(const QString &library) const;

    const MD5Sum& digest() const;

    const QLocale& locale() const;

    QString repository() const;
    QString buildVersion() const;

    quint32 version() const;

    void assertCompatible() const;
    void assertNotCorrupted(const QByteArray &compressed_data) const;

private:
    FileHeader(const QString &type_name,
               const QByteArray &compressed_data,
               const QByteArray &raw_data);

    FileHeader(const QStringList &type_names,
               const QByteArray &compressed_data,
               const QByteArray &raw_data);

    /** The username of the person who created this data */
    QString created_by;
    
    /** When this data was created */
    QDateTime created_when;
    
    /** The hostname of the machine on which this data was created */
    QString created_where;

    /** A string giving information about the system and libraries */
    QString system_info;

    /** The typename(s) of the top-level object(s) stored in this data */
    QStringList type_names;
    
    /** The source repository from which this code was downloaded */
    QString build_repository;

    /** Version string giving the build version of
        the program used to create this file */
    QString build_version; 

    /** The names and versions of the libraries loaded when this
        data was written */
    QByteArray required_libraries;

    /** The locale of the system used to create this file */
    QLocale system_locale;

    /** The digest of the data - this is used to verify that
        the data is not corrupted */
    MD5Sum data_digest;

    /** The size of the compressed data */
    quint32 compressed_size;
    
    /** The size of the uncompressed data */
    quint32 uncompressed_size;
    
    /** The version number of the file */
    quint32 version_number;
};

FileHeader getDataHeader(const QByteArray &data);
FileHeader getDataHeader(const QString &filename);

void registerLibrary(const QString &library,
                     quint32 version, quint32 min_supported_version);

quint32 getLibraryVersion(const QString &library);
quint32 getMinimumSupportedVersion(const QString &library);

class SIRESTREAM_EXPORT RegisterLibrary
{
public:
    RegisterLibrary()
    {}

    RegisterLibrary(const QString &library_name,
                    quint32 version, quint32 min_supported_version)
    {
        registerLibrary(library_name, version, min_supported_version);
    }
    
    ~RegisterLibrary()
    {}
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** This loads an object of type T from the passed blob of binary
    data. Note that this data *must* have been created by the "save"
    function below, and T must match the type of the object saved 
    in this data. Also note that this type must have been registered
    with the metatype system (via RegisterMetaType)
    
    \throw SireError::invalid_cast
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
T loadType(const QByteArray &data)
{
    QList< boost::tuple<boost::shared_ptr<void>,QString> > new_objs
            = SireStream::load(data);

    if ( new_objs.isEmpty() )
    {
        detail::throwStreamDataInvalidCast("NULL", T::typeName());
    }
    else if ( QLatin1String(T::typeName()) != new_objs.at(0).get<1>() )
    {
        detail::throwStreamDataInvalidCast(new_objs.at(0).get<1>(), T::typeName());
    }

    return T( *(static_cast<const T*>(new_objs.at(0).get<0>().get())) );
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
T loadType(const QString &filename)
{
    QList< boost::tuple<boost::shared_ptr<void>,QString> > new_objs
            = SireStream::load(filename);

    if ( new_objs.isEmpty() )
    {
        detail::throwStreamDataInvalidCast("NULL", T::typeName());
    }
    else if ( QLatin1String(T::typeName()) != new_objs.at(0).get<1>() )
    {
        detail::throwStreamDataInvalidCast(new_objs.at(0).get<1>(), T::typeName());
    }

    return T( *(static_cast<const T*>(new_objs.at(0).get<0>().get())) );
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
QByteArray save(const T &old_obj)
{
    return detail::streamDataSave( &old_obj, old_obj.what() );
}

template<class T0, class T1>
SIRE_OUTOFLINE_TEMPLATE
QByteArray save(const T0 &obj0, const T1 &obj1)
{
    QList< boost::tuple<const void*,const char*> > objects;
    
    objects.append( boost::tuple<const void*,const char*>( &obj0, obj0.what() ) );
    objects.append( boost::tuple<const void*,const char*>( &obj1, obj1.what() ) );
    
    return detail::streamDataSave( objects );
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
void saveToFile(const T &old_obj, const QString &filename)
{
    detail::streamDataSave( &old_obj, old_obj.what(), filename );
}

template<class T0, class T1>
SIRE_OUTOFLINE_TEMPLATE
void saveToFile(const T0 &obj0, const T1 &obj1, const QString &filename)
{
    QList< boost::tuple<const void*,const char*> > objects;
    
    objects.append( boost::tuple<const void*,const char*>( &obj0, obj0.what() ) );
    objects.append( boost::tuple<const void*,const char*>( &obj1, obj1.what() ) );
    
    detail::streamDataSave( objects, filename );
}

#endif // SIRE_SKIP_INLINE_FUNCTIONS

}

SIRE_EXPOSE_FUNCTION( SireStream::getDataHeader )
SIRE_EXPOSE_FUNCTION( SireStream::getLibraryVersion )
SIRE_EXPOSE_FUNCTION( SireStream::getMinimumSupportedVersion )

SIRE_EXPOSE_CLASS( SireStream::FileHeader )

SIRE_END_HEADER

#endif
