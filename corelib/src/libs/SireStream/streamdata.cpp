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

#include <QByteArray>
#include <QFile>
#include <QDataStream>
#include <QList>
#include <QMutex>
#include <QSysInfo>
#include <QtGlobal>
#include <QProcess>

#include <cstdlib>

#include <memory>
#include <boost/config.hpp>

#ifdef Q_OS_UNIX
    #include <unistd.h>
    #include <sys/utsname.h>
#endif
    
#include "streamdata.hpp"

#include "tostring.h"

#include "shareddatastream.h"
#include "registeralternativename.h"

#include "sire_version.h"        // CONDITIONAL_INCLUDE

#include "SireStream/version_error.h"

#include "SireStream/errors.h"
#include "SireError/errors.h"

#include <QDebug>

using namespace SireStream;

using boost::tuple;
using boost::shared_ptr;

namespace SireStream
{
namespace detail
{

/** This internal class is used to record the names and versions of the 
    libraries that are linked in the executable. This information is attached
    to the binary data saved for an object as an aid to ensure compatibility
    when reading this data back.
    
    @author Christopher Woods
*/
class LibraryInfo
{
public:
    LibraryInfo();
    ~LibraryInfo();

    static void registerLibrary(const QString &library, 
                                quint32 version, quint32 minversion);
    
    static QString getSupportReport(const QList< tuple<QString,quint32> > &libraries);
                           
    static void assertSupported(const QString &library, quint32 version);

    static QByteArray getLibraryHeader();
    
    static QList< tuple<QString,quint32> > readLibraryHeader(const QByteArray &header);

    static void checkLibraryHeader(const QByteArray &header);

    static quint32 getLibraryVersion(const QString &library);
    static quint32 getMinimumSupportedVersion(const QString &library);

private:
    static LibraryInfo& libraryInfo();
    
    QHash< QString, tuple<quint32,quint32> > library_info;
    QByteArray library_header;
};

/////////
///////// Implementation of LibraryInfo
/////////

LibraryInfo::LibraryInfo()
{}

LibraryInfo::~LibraryInfo()
{}

Q_GLOBAL_STATIC( QMutex, dataMutex );
Q_GLOBAL_STATIC( LibraryInfo, globalLibrary );

LibraryInfo& LibraryInfo::libraryInfo()
{
    LibraryInfo *global_library = globalLibrary();
    
    if (not global_library)
        throw std::exception();

    if (global_library->library_info.isEmpty())
    {
        //while we are here, register the SireStream and SireError libraries.
        //This has to be done here to prevent crashes caused by registration
        //being attempted before the static data in this file is initialised.
        global_library->library_info.insert( "SireError", tuple<quint32,quint32>(1,1) );
        global_library->library_info.insert( "SireStream", tuple<quint32,quint32>(1,1) );
    }
        
    return *global_library;
}

/** This function is used to register libraries as they are loaded */
void LibraryInfo::registerLibrary(const QString &library, 
                                  quint32 version, quint32 minversion)
{
    QMutexLocker lkr( dataMutex() );
    
    if (minversion > version)
        minversion = version;
    
    libraryInfo().library_info.insert( library, 
                                       tuple<quint32,quint32>(version,minversion) );
                                       
    libraryInfo().library_header.clear();
}

/** This returns a report about whether or not the provided list of libraries
    are supported. This returns an empty string if they are all supported,
    or a string containing every detail of the lack of support */
QString LibraryInfo::getSupportReport(const QList< tuple<QString,quint32> > &libraries)
{
    QMutexLocker lkr( dataMutex() );

    QStringList problems;

    for (QList< tuple<QString,quint32> >::const_iterator it = libraries.constBegin();
         it != libraries.constEnd();
         ++it)
    {
        QString library = it->get<0>();
        quint32 version = it->get<1>();
    
        if (not libraryInfo().library_info.contains(library))
        {
            if (library == "SireDB")
                //this is one of the old, and now removed libraries
                continue;
        
            problems.append( QObject::tr(
                "  (%1) The required library \"%2\" is missing.")
                    .arg(problems.count() + 1).arg(library) );
        }
        else
        {
            quint32 max_version = libraryInfo().library_info[library].get<0>();
            quint32 min_version = libraryInfo().library_info[library].get<1>();
            
            if (version > max_version)
            {
                problems.append( QObject::tr(
                    "  (%1) We need a newer version of the library \"%2\" (%3) "
                    "than the one available (%4).")
                        .arg(problems.count() + 1)
                        .arg(library).arg(version).arg(max_version) );
            }
            else if (version < min_version)
            {
                problems.append( QObject::tr(
                    "   (%1) We need an older version of the library \"%2\" (%3) "
                    "than the one available (%4).")
                        .arg(problems.count() + 1)
                        .arg(library).arg(version).arg(min_version) );
            }
        }
    }
    
    if (problems.isEmpty())
        return QString();
    else
        return problems.join("\n");
}
                       
/** Assert that the library 'library' with version 'version' is supported 

    \throw SireStream::version_error
*/
void LibraryInfo::assertSupported(const QString &library, quint32 version)
{
    QMutexLocker lkr( dataMutex() );
    
    if (not libraryInfo().library_info.contains(library))
    {
        throw version_error( QObject::tr(
            "The required library (%1) is not available.")
                .arg(library), CODELOC );
    }
    else
    {
        quint32 max_version = libraryInfo().library_info[library].get<0>();
        quint32 min_version = libraryInfo().library_info[library].get<1>();
        
        if (version > max_version)
        {
            throw version_error( QObject::tr(
                    "We need a newer version of the library \"%1\" (%2) "
                    "than the one available (%3).")
                        .arg(library).arg(version).arg(max_version), CODELOC );
        }
        else if (version < min_version)
        {
            throw version_error( QObject::tr(
                    "We need an older version of the library \"%1\" (%2) "
                    "than the one available (%3).")
                        .arg(library).arg(version).arg(min_version), CODELOC );
        }
    }
}

/** Return the header data that describes the libraries linked with 
    this executable */
QByteArray LibraryInfo::getLibraryHeader()
{
    QMutexLocker lkr( dataMutex() );
    
    if (libraryInfo().library_header.isEmpty())
    {
        QByteArray header;
    
        QDataStream ds( &header, QIODevice::WriteOnly );
        
        //this is version 1
        ds << quint32(1);
        
        //hard code that version 1 uses the Qt 4.2 data format
        ds.setVersion(QDataStream::Qt_4_2);
        
        //write the number of libraries
        ds << quint32( libraryInfo().library_info.count() );
        
        for (QHash< QString,tuple<quint32,quint32> >::const_iterator 
                                         it = libraryInfo().library_info.constBegin();
             it != libraryInfo().library_info.constEnd();
             ++it)
        {
            ds << it.key() << it.value().get<0>();
        }
        
        libraryInfo().library_header = qCompress(header, 9);
    }
    
    return libraryInfo().library_header;
}

/** This reads the library header and returns the libraries and versions
    used when writing the binary data */
QList< tuple<QString,quint32> > LibraryInfo::readLibraryHeader(const QByteArray &header)
{
    if (header.isEmpty())
        throw SireError::program_bug( QObject::tr(
                "For some reason, the library header file is empty! This means that Sire "
                "is unable to read the saved file."), CODELOC );

    QByteArray unpacked_header = qUncompress(header);

    QDataStream ds(unpacked_header);
    
    quint32 version;
    
    ds >> version;
    
    QList< tuple<QString,quint32> > libraries;
    
    if (version == 1)
    {
        //version 1 uses the Qt 4.2 format
        ds.setVersion(QDataStream::Qt_4_2);

        quint32 nlibs;
        ds >> nlibs;
        
        for (quint32 i=0; i<nlibs; ++i)
        {
            QString library;
            quint32 lib_version;
            
            ds >> library >> lib_version;
            
            libraries.append( tuple<QString,quint32>(library, lib_version) );
        }
    }
    else
        throw version_error( QObject::tr(
            "Cannot even read the information about supported libraries, as not "
            "even that data is supported. This can only read version 1, but "
            "the data is version %1.").arg(version), CODELOC );
    
    return libraries;
}

/** This reads a library header and checks that all of the requirements are
    satisfied 
    
    \throw version_error
*/
void LibraryInfo::checkLibraryHeader(const QByteArray &header)
{
    QList< tuple<QString,quint32> > libraries = readLibraryHeader(header);
    
    QString report = getSupportReport(libraries);
    
    if (not report.isEmpty())
        throw version_error( QObject::tr(
            "There are incompatibilities between the libraries required "
            "to read this data and the libraries available to this program.\n%1")
                .arg(report), CODELOC );
}

quint32 LibraryInfo::getLibraryVersion(const QString &library)
{
    QMutexLocker lkr( dataMutex() );
    
    if (not libraryInfo().library_info.contains(library))
    {
        throw SireError::unsupported( QObject::tr(
            "The library %1 is not available to this program. Available libraries "
            "are %2.")
                .arg(library)
                .arg( Sire::toString(libraryInfo().library_info.keys()) ),
                    CODELOC );
    }
    
    return libraryInfo().library_info.value(library).get<0>();
}

quint32 LibraryInfo::getMinimumSupportedVersion(const QString &library)
{
    QMutexLocker lkr( dataMutex() );
    
    if (not libraryInfo().library_info.contains(library))
    {
        throw SireError::unsupported( QObject::tr(
            "The library %1 is not available to this program. Available libraries "
            "are %2.")
                .arg(library)
                .arg( Sire::toString(libraryInfo().library_info.keys()) ),
                    CODELOC );
    }
    
    return libraryInfo().library_info.value(library).get<1>();
}

} // end of namespace detail
} // end of namespace SireStream

/////////
///////// Implementation of FileHeader
/////////

/** Write the header to the data stream */
QDataStream SIRESTREAM_EXPORT &operator<<(QDataStream &ds,
                                          const FileHeader &header)
{
    ds << header.version();

    //versions 1 and 2 uses the Qt 4.2 data format
    ds.setVersion(QDataStream::Qt_4_2);

    QByteArray data;
    QDataStream ds2(&data, QIODevice::WriteOnly);
    
    ds2.setVersion(QDataStream::Qt_4_2);
    
    ds2 << header.created_by
        << header.created_when
        << header.created_where
        << header.system_info;
       
    if (header.version() == 2) 
        ds2 << header.type_names;
    
    else if (header.version() == 1)
    {
        if (header.type_names.count() != 1)
            throw SireError::version_error( QObject::tr(
                "Version 1 of the FileHeader format does not support multiple "
                "type names, yet the types ( %1 ) are being streamed!")
                    .arg(header.type_names.join(", ")), CODELOC );
                    
        ds2 << header.type_names.at(0);
    }
        
    ds2 << header.build_repository
        << header.build_version
        << header.required_libraries
        << header.system_locale
        << header.data_digest
        << header.compressed_size
        << header.uncompressed_size;

    ds << qCompress(data, 9);
    
    return ds;
}

/** Read the header from the data stream */
QDataStream SIRESTREAM_EXPORT &operator>>(QDataStream &ds,
                                          FileHeader &header)
{
    quint32 version;
    ds >> version;
    
    if (version == 2)
    {
        //Version 2 uses the Qt 4.2 data format
        ds.setVersion(QDataStream::Qt_4_2);
    
        QByteArray data;
        ds >> data;

        data = qUncompress(data);
    
        QDataStream ds2(data);
    
        ds2.setVersion(QDataStream::Qt_4_2);
    
        ds2 >> header.created_by
            >> header.created_when
            >> header.created_where
            >> header.system_info
            >> header.type_names
            >> header.build_repository
            >> header.build_version
            >> header.required_libraries
            >> header.system_locale
            >> header.data_digest
            >> header.compressed_size
            >> header.uncompressed_size;
            
        header.version_number = version;
    }
    else if (version == 1)
    {
        //Version 1 uses the Qt 4.2 data format
        ds.setVersion(QDataStream::Qt_4_2);
    
        QByteArray compressed_data;
        ds >> compressed_data;

        QByteArray data = qUncompress(compressed_data);
    
        QDataStream ds2(data);
    
        ds2.setVersion(QDataStream::Qt_4_2);
    
        QString type_name;
    
        ds2 >> header.created_by
            >> header.created_when
            >> header.created_where
            >> header.system_info
            >> type_name
            >> header.build_repository
            >> header.build_version
            >> header.required_libraries
            >> header.system_locale
            >> header.data_digest
            >> header.compressed_size
            >> header.uncompressed_size;
            
        header.type_names.clear();
        header.type_names.append(type_name);
            
        header.version_number = version;
    }
    else
        throw version_error( QObject::tr(
            "The header version (%1) is not recognised. Only header version "
            "1+2 are supported in this program.")
                .arg(version), CODELOC );
        
    return ds;
}

/** Null constructor */
FileHeader::FileHeader() : compressed_size(0), uncompressed_size(0), version_number(0)
{}

static QString *system_info(0);

static const QString& getSystemInfo()
{
    if (system_info)
        return *system_info;

    QStringList lines;
    
    #ifdef Q_OS_MAC
        lines.append( "Platform: Mac OS" );
    
        switch (QSysInfo::MacintoshVersion)
        {
            case QSysInfo::MV_9:
                lines.append( "Mac Version: OS 9" );
                break;
            case QSysInfo::MV_10_0:
                lines.append( "Mac Version: OS X 10.0 (Cheetah)" );
                break;
            case QSysInfo::MV_10_1:
                lines.append( "Mac Version: OS X 10.1 (Puma)" );
                break;
            case QSysInfo::MV_10_2:
                lines.append( "Mac Version: OS X 10.2 (Jaguar)" );
                break;
            case QSysInfo::MV_10_3:
                lines.append( "Mac Version: OS X 10.3 (Panther)" );
                break;
            case QSysInfo::MV_10_4:
                lines.append( "Mac Version: OS X 10.4 (Tiger)" );
                break;
            case QSysInfo::MV_10_5:
                lines.append( "Mac Version: OS X 10.5 (Leopard)" );
                break;
            default:
                lines.append( "Mac Version: Unknown" );
                break;
        }
    #else
    #ifdef Q_OS_LINUX
        lines.append( "Platform: Linux" );
    #else
    #ifdef Q_OS_UNIX
        lines.append( "Platform: UNIX" );

        #ifdef Q_OS_AIX
        lines.append( "UNIX flavour: AIX" );
        #endif
        
        #ifdef Q_OS_BSD4
        lines.append( "UNIX flavour: BSD 4.4" );
        #endif

        #ifdef Q_OS_BSDI
        lines.append( "UNIX flavour: BSD/OS" );
        #endif

        #ifdef Q_OS_DARWIN
        lines.append( "UNIX flavour: Darwin" );
        #endif

        #ifdef Q_OS_DGUX
        lines.append( "UNIX flavour: DG/UX" );
        #endif

        #ifdef Q_OS_DYNIX
        lines.append( "UNIX flavour: DYNIX/ptx" );
        #endif
        
        #ifdef Q_OS_FREEBSD
        lines.append( "UNIX flavour: FreeBSD" );
        #endif

        #ifdef Q_OS_HPUX
        lines.append( "UNIX flavour: HP-UX" );
        #endif

        #ifdef Q_OS_HURD
        lines.append( "UNIX flavour: GNU Hurd" );
        #endif

        #ifdef Q_OS_IRIX
        lines.append( "UNIX flavour: SGI Irix" );
        #endif

        #ifdef Q_OS_LYNX
        lines.append( "UNIX flavour: LynxOS" );
        #endif

        #ifdef Q_OS_NETBSD
        lines.append( "UNIX flavour: NetBSD" );
        #endif

        #ifdef Q_OS_OPENBSD
        lines.append( "UNIX flavour: OpenBSD" );
        #endif

        #ifdef Q_OS_OSF
        lines.append( "UNIX flavour: HP Tru64 UNIX" );
        #endif

        #ifdef Q_OS_QNX6
        lines.append( "UNIX flavour: QNX RTP 6.1" );
        #endif

        #ifdef Q_OS_QNX
        lines.append( "UNIX flavour: QNX" );
        #endif

        #ifdef Q_OS_RELIANT
        lines.append( "UNIX flavour: Reliant UNIX" );
        #endif
        
        #ifdef Q_OS_SCO
        lines.append( "UNIX flavour: SCO OpenServer 5" );
        #endif

        #ifdef Q_OS_SOLARIS
        lines.append( "UNIX flavour: Sun Solaris" );
        #endif

        #ifdef Q_OS_ULTRIX
        lines.append( "UNIX flavour: DEC Ultrix" );
        #endif

        #ifdef Q_OS_UNIXWARE
        lines.append( "UNIX flavour: UnixWare 7, Open UNIX 8" );
        #endif

    #else
    #ifdef Q_OS_WIN32
        lines.append( "Platform: Windows" );
        
        switch (QSysInfo::windowsVersion())
        {
            case QSysInfo::WV_32s:
                lines.append( "Windows Version: 3.1" );
                break;
            case QSysInfo::WV_95:
                lines.append( "Windows Version: 95" );
                break;
            case QSysInfo::WV_98:
                lines.append( "Windows Version: 98" );
                break;
            case QSysInfo::WV_Me:
                lines.append( "Windows Version: Me" );
                break;
            case QSysInfo::WV_NT:
                lines.append( "Windows Version: NT" );
                break;
            case QSysInfo::WV_2000:
                lines.append( "Windows Version: 2000" );
                break;
            case QSysInfo::WV_XP:
                lines.append( "Windows Version: XP" );
                break;
            case QSysInfo::WV_2003:
                lines.append( "Windows Version: Server 2003" );
                break;
            case QSysInfo::WV_VISTA:
                lines.append( "Windows Version: Vista" );
                break;
            case QSysInfo::WV_CE:
                lines.append( "Windows Version: CE" );
                break;
            case QSysInfo::WV_CENET:
                lines.append( "Windows Version: CE .NET" );
                break;
            case QSysInfo::WV_CE_5:
                lines.append( "Windows Version: CE 5" );
                break;
            case QSysInfo::WV_CE_6:
                lines.append( "Windows Version: CE 6" );
                break;
            default:
                lines.append( "Windows Version: Unknown" );
        }
    #else
    #ifdef Q_OS_OS2
        lines.append( "System: OS/2" );
    #endif  //OS/2
    #endif  //windows
    #endif  //unix
    #endif  //linux
    #endif  //mac

    #ifdef Q_OS_UNIX
        utsname sysname;
        
        if (uname(&sysname) != -1)
        {
            lines.append( QString("UNIX system: %1").arg( sysname.sysname ) );
            lines.append( QString("UNIX release: %1").arg( sysname.release ) );
            lines.append( QString("UNIX version: %1").arg( sysname.version ) );
            lines.append( QString("UNIX machine: %1").arg( sysname.machine ) );
        }
    #endif
        
    #ifdef Q_OS_CYGWIN
    lines.append( "UNIX flavour: Cygwin" );
    #endif
    
    #ifdef Q_CC_BOR
    lines.append( "Compiler: Borland/Turbo C++" );
    #endif

    #ifdef Q_CC_CDS
    lines.append( "Compiler: Reliant C++" );
    #endif

    #ifdef Q_CC_COMEAU
    lines.append( "Compiler: Comeau C++" );
    #endif

    #ifdef Q_CC_DEC
    lines.append( "Compiler: DEC C++" );
    #endif

    #ifdef Q_CC_EDG
    lines.append( "Compiler: Edison Design Group C++" );
    #endif

    #ifdef Q_CC_GHS
    lines.append( "Compiler: Green Hills Optimizing C++ Compiler" );
    #endif
    
    #ifdef Q_CC_GNU
    lines.append( QString("Compiler: GNU C++ (%1.%2.%3)")
                        .arg( __GNUC__)
                        .arg( __GNUC_MINOR__)
                        .arg( __GNUC_PATCHLEVEL__ ) );
    #endif

    #ifdef Q_CC_HIGHC
    lines.append( "Compiler: MetaWare High C/C++" );
    #endif

    #ifdef Q_CC_HPACC
    lines.append( "Compiler: HP aC++" );
    #endif
    
    #ifdef Q_CC_INTEL
    lines.append( "Compiler: Intel C++" );
    #endif

    #ifdef Q_CC_KAI
    lines.append( "Compiler: KAI C++" );
    #endif

    #ifdef Q_CC_MIPS
    lines.append( "Compiler: MIPSpro C++" );
    #endif

    #ifdef Q_CC_MSVC
    lines.append( "Compiler: Microsoft Visual C/C++" );
    #endif

    #ifdef Q_CC_MWERKS
    lines.append( "Compiler: Metrowerks CodeWarrior" );
    #endif

    #ifdef Q_CC_OC
    lines.append( "Compiler: CenterLine C++" );
    #endif

    #ifdef Q_CC_PGI
    lines.append( "Compiler: Portland Group C++" );
    #endif

    #ifdef Q_CC_SUN
    lines.append( "Compiler: Forte Developer or Sun Studio C++" );
    #endif

    #ifdef Q_CC_SYM
    lines.append( "Compiler: Digital Mars C/C++" );
    #endif

    #ifdef Q_CC_USLC
    lines.append( "Compiler: SCO OUDK and UDK" );
    #endif

    #ifdef Q_CC_WAT
    lines.append( "Compiler: Watcom C++" );
    #endif

    lines.append( QString("Wordsize: %1 bit").arg( QSysInfo::WordSize ) );
    
    switch (QSysInfo::ByteOrder)
    {
        case QSysInfo::BigEndian:
            lines.append( "ByteOrder: Big endian" );
            break;
        case QSysInfo::LittleEndian:
            lines.append( "ByteOrder: Little endian" );
            break;
        default:
            lines.append( "ByteOrder: Unknown" );
            break;
    }

    lines.append( QString("Qt runtime version: %1").arg( qVersion() ) );
    lines.append( QString("Qt compile version: %1").arg( QT_VERSION_STR ) );

    lines.append( QString("Sire compile version: %1.%2.%3")
                        .arg(SIRE_VERSION_MAJOR)
                        .arg(SIRE_VERSION_MINOR)
                        .arg(SIRE_VERSION_PATCH) );
                        
    lines.append( QString("Compile flags: %1").arg(SIRE_COMPILE_FLAGS) );
    lines.append( QString("Link flags: %1").arg(SIRE_LINK_FLAGS) );

    QString info = lines.join("\n");
    
    system_info = new QString(info);
    
    return *system_info;
}

/** Internal constructor used by streamDataSave */
FileHeader::FileHeader(const QStringList &typ_names,
                       const QByteArray &compressed_data,
                       const QByteArray &raw_data) : version_number(0)
{
    //these two may be UNIX only...
    created_by = std::getenv("USER");

    char buffer[128];
    gethostname(buffer, 128);
    created_where = buffer;

    created_when = QDateTime::currentDateTime();
    
    type_names = typ_names;
    
    build_repository = SIRE_REPOSITORY_URL;
    build_version = SIRE_REPOSITORY_VERSION;
    
    required_libraries = detail::LibraryInfo::getLibraryHeader();
    
    data_digest = MD5Sum(compressed_data);
    
    compressed_size = compressed_data.count();
    uncompressed_size = raw_data.count();

    system_info = getSystemInfo();
}

/** Copy constructor */
FileHeader::FileHeader(const FileHeader &other)
           : created_by(other.created_by), created_when(other.created_when),
             created_where(other.created_where), 
             system_info(other.system_info),
             type_names(other.type_names),
             build_repository(other.build_repository),
             build_version(other.build_version),
             required_libraries(other.required_libraries),
             system_locale(other.system_locale),
             data_digest(other.data_digest),
             compressed_size(other.compressed_size),
             uncompressed_size(other.uncompressed_size),
             version_number(other.version_number)
{}

/** Destructor */
FileHeader::~FileHeader()
{}

/** Copy assignment operator */
FileHeader& FileHeader::operator=(const FileHeader &other)
{
    if (this != &other)
    {
        created_by = other.created_by;
        created_when = other.created_when;
        created_where = other.created_where;
        system_info = other.system_info;
        type_names = other.type_names;
        build_repository = other.build_repository;
        build_version = other.build_version;
        system_locale = other.system_locale;
        required_libraries = other.required_libraries;
        data_digest = other.data_digest;
        compressed_size = other.compressed_size;
        uncompressed_size = other.uncompressed_size;
        version_number = other.version_number;
    }
    
    return *this;
}

/** Return a string representation of this header */
QString FileHeader::toString() const
{
    QStringList sysinfo = system_info.split("\n");

    return QObject::tr("*************************************\n"
                       "* Object type(s) : %1\n"
                       "* Created by     : %2\n"
                       "* Creation date  : %3\n"
                       "* Created on     : %4\n"
                       "* System         : %5\n"
                       "* Sire Version   : %6\n"
                       "* Repository     : %7\n"
                       "* Packed size    : %8 kB\n"
                       "* Unpacked size  : %9 kB\n"
                       "* Country        : %10\n"
                       "* Language       : %11\n"
                       "*************************************\n")
               .arg( type_names.join(" "), created_by, created_when.toString(),
                     created_where )
               .arg( sysinfo.join("\n*              : "), 
                     build_version, build_repository )
               .arg(double(compressed_size)/1024.0)
               .arg(double(uncompressed_size)/1024.0)
               .arg( QLocale::countryToString(system_locale.country()),
                     QLocale::languageToString(system_locale.language()) );
}

/** Return the username of whoever created this data */
const QString& FileHeader::createdBy() const
{
    return created_by;
}

/** Return the date and time when this was created */
const QDateTime& FileHeader::createdWhen() const
{
    return created_when;
}

/** Where this file was created (the name of the machine) */
const QString& FileHeader::createdWhere() const
{
    return created_where;
}

/** Return the minimum memory the will be necessary to read the file */
quint32 FileHeader::requiredMemory() const
{
    return 2*uncompressed_size + compressed_size;
}

/** Return the compression ratio of the file */
double FileHeader::compressionRatio() const
{
    if (uncompressed_size == 0)
        return 1;
    else
    {
        return double(compressed_size) / uncompressed_size;
    }
}

/** Return the name of the data type of the object in this data

    \throw SireError::invalid_state
*/
QString FileHeader::dataType() const
{
    if (type_names.count() != 1)
        throw SireError::invalid_state( QObject::tr(
            "To get the data type, there must be just one object in the data. "
            "The number of objects in this data is %1 ( %2 ).")
                .arg(type_names.count()).arg(type_names.join(", ")),
                    CODELOC );
                    
    return type_names.at(0);
}

/** Return the name(s) of the top-level data type in this data */
const QStringList& FileHeader::dataTypes() const
{
    return type_names;
}

/** Return information about the system on which this data was written */
const QString& FileHeader::systemInfo() const
{
    return system_info;
}

/** Return the locale in which this data was written. This is useful
    as it can help with the support of multiple languages
    (as the person who saved the file may not necessarily be using
     English ;-) */
const QLocale& FileHeader::locale() const
{
    return system_locale;
}

/** Return the list of libraries required to load this data */
QStringList FileHeader::requiredLibraries() const
{
    QList< tuple<QString,quint32> > libs =
                    detail::LibraryInfo::readLibraryHeader(required_libraries);
                    
    QStringList libraries;
    
    for (QList< tuple<QString,quint32> >::const_iterator it = libs.constBegin();
         it != libs.constEnd();
         ++it)
    {
        libraries.append( it->get<0>() );
    }
    
    return libraries;
}

/** Does this data require that the library 'library' be loaded? */
bool FileHeader::requireLibrary(const QString &library) const
{
    QStringList libraries = this->requiredLibraries();
    
    foreach (const QString &lib, libraries)
    {
        if (lib == library)
            return true;
    }
    
    return false;
}

/** Return the version number required of the library 'library'. This
    returns 0 if this library isn't required. */
quint32 FileHeader::requiredVersion(const QString &library) const
{
    QList< tuple<QString,quint32> > libs =
                    detail::LibraryInfo::readLibraryHeader(required_libraries);
                    
    for (QList< tuple<QString,quint32> >::const_iterator it = libs.constBegin();
         it != libs.constEnd();
         ++it)
    {
        if (it->get<0>() == library)
            return it->get<1>();
    }
    
    return 0;
}

/** Return the digest of the data - this is used to check 
    for any data corruption */
const MD5Sum& FileHeader::digest() const
{
    return data_digest;
}

/** Return the repository from which this source code was downloaded */
QString FileHeader::repository() const
{
    return build_repository;
}

/** Return the version of the source code from the repository */
QString FileHeader::buildVersion() const
{
    return build_version;
}

/** Assert that the libraries required are compatible with what has
    been loaded */
void FileHeader::assertCompatible() const
{
    detail::LibraryInfo::checkLibraryHeader(required_libraries);
}

/** Assert that the data in 'compressed_data' is not corrupt */
void FileHeader::assertNotCorrupted(const QByteArray &compressed_data) const
{
    if (quint32(compressed_data.size()) != compressed_size)
        throw SireStream::corrupted_data( QObject::tr(
            "The data for the object(s) [ %1 ] appears to be corrupt as it "
            "is the wrong size (%2 bytes vs. the expected %3 bytes)")
                .arg(type_names.join(", ")).arg(compressed_data.size())
                .arg(compressed_size), CODELOC );
                
    MD5Sum new_digest(compressed_data);
    
    if (data_digest != new_digest)
        throw SireStream::corrupted_data( QObject::tr(
            "The data for the object(s) [ %1 ] appears to be corrupt as the "
            "digests don't match (%2 vs. %3)")
                .arg(type_names.join(", "))
                .arg(new_digest.toString(), data_digest.toString()),
                    CODELOC );
}

/** Return the master version number for the file - this version number 
    is changed only when the file format is completely changed (e.g. we
    move away from using a compressed header, then the compressed object)
    
    Currently, we only use version 1, which has this format;
    
    SIRE_MAGIC_NUMBER  (quint32 = 251785387)
    VERSION_NUMBER     (quint32 = 2)
    QByteArray         (compressed array containing the file header)
    QByteArray         (compressed array containing the saved object)
    
    All of this is written using Qt datastream format for Qt 4.2
*/
quint32 FileHeader::version() const
{
    if (version_number == 0)
        //the version has not been set - so use the latest version
        //available - which is '2' in this case
        return 2;
    else
        return version_number;
}

/////////
///////// Implementation of free functions
/////////

namespace SireStream
{
namespace detail
{

void SIRESTREAM_EXPORT throwStreamDataInvalidCast(const QString &load_type,
                                                  const QString &cast_type)
{
    throw SireError::invalid_cast( QObject::tr(
        "Cannot load the binary data for the object of type %1 into "
        "an object of type %2.")
            .arg(load_type).arg(cast_type), CODELOC );
}

static quint32 SIRE_MAGIC_NUMBER(251785387);

static int RESERVE_SIZE = 48 * 1024 * 1024;

/** Save the object pointed to by 'object' with type 'type_name' to a binary
    array and return the array.
    
    Currently, we use version 1 of the format, which has;
    
    SIRE_MAGIC_NUMBER  (quint32 = 251785387)
    VERSION_NUMBER     (quint32 = 1)
    QByteArray         (compressed array containing the file header)
    QByteArray         (compressed array containing the saved object)
    
    All of this is written using Qt datastream format for Qt 4.2
*/
QByteArray SIRESTREAM_EXPORT streamDataSave( 
                               const QList< tuple<const void*,const char*> > &objects )
{
    FileHeader header;
    
    if (header.version() == 1 or header.version() == 2)
    {
        int nobjects = objects.count();
        
        if (nobjects == 0)
            return QByteArray();
            
        else if (nobjects > 1 and header.version() == 1)
            throw SireError::version_error( QObject::tr(
                "Version 1 of the global Sire format is incapable of saving "
                "multiple objects, while you are trying to save %1 objects.")
                    .arg(nobjects), CODELOC );

        //use the QMetaType streaming function to save this object
        QByteArray object_data;
    
        //reserve at least 48MB of space (most things shouldn't be this big)
        object_data.reserve( RESERVE_SIZE );

        if (object_data.capacity() != RESERVE_SIZE)
            qWarning() << "Possible memory allocation error!";
    
        QDataStream ds2(&object_data, QIODevice::WriteOnly);

        //version 1 of the format uses Qt 4.2 datastream format
        ds2.setVersion( QDataStream::Qt_4_2 );
        
        #ifdef BOOST_NO_CXX11_SMART_PTR
          std::auto_ptr<SharedDataStream> sds;
        #else
          std::unique_ptr<SharedDataStream> sds;
        #endif
        
        if (nobjects > 1)
            //create a shared data stream so that sub-objects in 
            //top-level objects can be shared
            sds.reset( new SharedDataStream(ds2) );
    
        QStringList type_names;
    
        for (int i=0; i<nobjects; ++i)
        {
            QString type_name = objects.at(i).get<1>();
        
            //get the ID number of this type
            int id = QMetaType::type( type_name.toLatin1().constData() );

            if ( id == 0 or not QMetaType::isRegistered(id) )
                throw SireError::unknown_type(QObject::tr(
                    "The object with type \"%1\" does not appear to have been "
                    "registered with QMetaType. It cannot be streamed! (%2, %3)")
                        .arg(type_name).arg(id).arg(QMetaType::isRegistered(id)), 
                            CODELOC);

            if (not QMetaType::save(ds2, id, objects.at(i).get<0>()))
                throw SireError::program_bug(QObject::tr(
                    "There was an error saving the object of type \"%1\". "
                    "Has the programmer remembered to add a RegisterMetaType "
                    "for this class?")
                        .arg(type_name), CODELOC);

            type_names.append(type_name);
        }

        //compress the object data (level 3 compression seems best, giving
        //about a ten-fold reduction for only a 30% increase in serialisation time)
        QByteArray compressed_object_data = qCompress(object_data, 3);

        //now write a header for the object
        header = FileHeader( type_names, compressed_object_data, object_data );
    
        //clear the uncompressed data to save space
        object_data = QByteArray();
    
        QByteArray data;
        data.reserve( RESERVE_SIZE );
        
        if (data.capacity() != RESERVE_SIZE)
            qWarning() << "Possible memory allocation error!";

        QDataStream ds(&data, QIODevice::WriteOnly);
    
        //write a magic number - then the header
        ds << SIRE_MAGIC_NUMBER;
        ds << header;

        //now write the object data
        ds << compressed_object_data;

        return data;
    }
    else
        throw version_error( QObject::tr(
            "Cannot write the object information, as it should be written using "
            "the global Sire format %1, while we can only write versions 1 and 2.")
                    .arg(header.version()), CODELOC );

    return QByteArray();
}

/** Save the object pointed to by 'object' with type 'type_name' to a binary
    array and return the array.
    
    Currently, we use version 1 of the format, which has;
    
    SIRE_MAGIC_NUMBER  (quint32 = 251785387)
    VERSION_NUMBER     (quint32 = 1)
    QByteArray         (compressed array containing the file header)
    QByteArray         (compressed array containing the saved object)
    
    All of this is written using Qt datastream format for Qt 4.2
*/
QByteArray SIRESTREAM_EXPORT streamDataSave( 
                               const QList< tuple<shared_ptr<void>,QString> > &objects )
{
    FileHeader header;
    
    if (header.version() == 1 or header.version() == 2)
    {
        int nobjects = objects.count();
        
        if (nobjects == 0)
            return QByteArray();
            
        else if (nobjects > 1 and header.version() == 1)
            throw SireError::version_error( QObject::tr(
                "Version 1 of the global Sire format is incapable of saving "
                "multiple objects, while you are trying to save %1 objects.")
                    .arg(nobjects), CODELOC );

        //use the QMetaType streaming function to save this object
        QByteArray object_data;
    
        //reserve at least 48MB of space (most things shouldn't be this big)
        object_data.reserve( RESERVE_SIZE );

        if (object_data.capacity() != RESERVE_SIZE)
            qWarning() << "Possible memory allocation error!";
    
        QDataStream ds2(&object_data, QIODevice::WriteOnly);

        //version 1 of the format uses Qt 4.2 datastream format
        ds2.setVersion( QDataStream::Qt_4_2 );

        #ifdef BOOST_NO_CXX11_SMART_PTR
          std::auto_ptr<SharedDataStream> sds;
        #else
          std::unique_ptr<SharedDataStream> sds;
        #endif
        
        if (nobjects > 1)
            //create a shared data stream so that sub-objects in 
            //top-level objects can be shared
            sds.reset( new SharedDataStream(ds2) );
    
        QStringList type_names;
    
        for (int i=0; i<nobjects; ++i)
        {
            QString type_name = objects.at(i).get<1>();
        
            //get the ID number of this type
            int id = QMetaType::type( type_name.toLatin1().constData() );

            if ( id == 0 or not QMetaType::isRegistered(id) )
                throw SireError::unknown_type(QObject::tr(
                    "The object with type \"%1\" does not appear to have been "
                    "registered with QMetaType. It cannot be streamed! (%2, %3)")
                        .arg(type_name).arg(id).arg(QMetaType::isRegistered(id)), 
                            CODELOC);

            if (not QMetaType::save(ds2, id, objects.at(i).get<0>().get()))
                throw SireError::program_bug(QObject::tr(
                    "There was an error saving the object of type \"%1\". "
                    "Has the programmer remembered to add a RegisterMetaType "
                    "for this class?")
                        .arg(type_name), CODELOC);

            type_names.append(type_name);
        }

        //compress the object data (level 3 compression seems best, giving
        //about a ten-fold reduction for only a 30% increase in serialisation time)
        QByteArray compressed_object_data = qCompress(object_data, 3);

        //now write a header for the object
        header = FileHeader( type_names, compressed_object_data, object_data );
    
        //clear the uncompressed data to save space
        object_data = QByteArray();
    
        QByteArray data;
        data.reserve( RESERVE_SIZE );
        
        if (data.capacity() != RESERVE_SIZE)
            qWarning() << "Possible memory allocation error!";

        QDataStream ds(&data, QIODevice::WriteOnly);
    
        //write a magic number - then the header
        ds << SIRE_MAGIC_NUMBER;
        ds << header;

        //now write the object data
        ds << compressed_object_data;

        return data;
    }
    else
        throw version_error( QObject::tr(
            "Cannot write the object information, as it should be written using "
            "the global Sire format %1, while we can only write versions 1 and 2.")
                    .arg(header.version()), CODELOC );

    return QByteArray();
}

/** Overloaded function that saves the object directly to a file, rather than
    to an array */
void SIRESTREAM_EXPORT streamDataSave( 
                            const QList< tuple<const void*,const char*> > &objects, 
                            const QString &filename )
{
    QFile f(filename);
    
    if (not f.open(QIODevice::WriteOnly))
        throw SireError::file_error(f, CODELOC);

    QByteArray data = streamDataSave(objects);

    if (f.write(data) == -1)
        throw SireError::file_error( QObject::tr(
            "There was an error writing to the file %1. Is there enough space to "
            "to write a file of %d bytes?").arg(filename).arg(data.count()),
                CODELOC );
}

/** Overloaded function that saves the object directly to a file, rather than
    to an array */
void SIRESTREAM_EXPORT streamDataSave( 
                            const QList< tuple<shared_ptr<void>,QString> > &objects, 
                            const QString &filename )
{
    QFile f(filename);
    
    if (not f.open(QIODevice::WriteOnly))
        throw SireError::file_error(f, CODELOC);

    QByteArray data = streamDataSave(objects);

    if (f.write(data) == -1)
        throw SireError::file_error( QObject::tr(
            "There was an error writing to the file %1. Is there enough space to "
            "to write a file of %d bytes?").arg(filename).arg(data.count()),
                CODELOC );
}

QByteArray SIRESTREAM_EXPORT streamDataSave( const void *object, const char *type_name )
{
    QList< tuple<const void*,const char *> > objects;
    
    objects.append( tuple<const void*,const char*>(object,type_name) );
    
    return streamDataSave(objects);
}

void SIRESTREAM_EXPORT streamDataSave( const void *object, const char *type_name,
                                       const QString &filename )
{
    QList< tuple<const void*,const char *> > objects;
    
    objects.append( tuple<const void*,const char*>(object,type_name) );
    
    streamDataSave(objects, filename);
}

} // end of namespace detail
} // end of namespace SireStream

namespace SireStream
{

/** This function is called by each Sire library when it is loaded
    to register the library with the streaming system. You should
    not call this function yourself! */
void SIRESTREAM_EXPORT registerLibrary(const QString &library,
                                       quint32 version, 
                                       quint32 min_supported_version)
{
    detail::LibraryInfo::registerLibrary(library, version, min_supported_version);
}

/** Return the version of the loaded library called 'library'

    \throw SireError::unsupported
*/
quint32 SIRESTREAM_EXPORT getLibraryVersion(const QString &library)
{
    return detail::LibraryInfo::getLibraryVersion(library);
}

/** Return the minimum data version that the library 'library' is capable
    of reading
    
    \throw SireError::unsupported
*/
quint32 SIRESTREAM_EXPORT getMinimumSupportedVersion(const QString &library)
{
    return detail::LibraryInfo::getMinimumSupportedVersion(library);
}

using namespace SireStream::detail;

/** This loads an object from the passed blob of binary data. This binary
    data *must* have been created by the "save" function below. */
QList< tuple<shared_ptr<void>,QString> > SIRESTREAM_EXPORT load(const QByteArray &data)
{
    QList< tuple<shared_ptr<void>,QString> > loaded_objects;

    QDataStream ds(data);
    
    //read the magic
    quint32 magic;
    ds >> magic;
    
    if (magic != SIRE_MAGIC_NUMBER)
        throw version_error( QObject::tr(
            "This data does not appear to have been written by the SireStream::save() "
            "function."), CODELOC );
    
    //read the version and the libraries and versions
    //(die if the required libraries aren't loaded!)
    FileHeader header;
    ds >> header;
    
    header.assertCompatible();
    
    if (header.dataTypes().isEmpty())
    {
        //this is a null pointer!
        loaded_objects.append( 
            tuple<shared_ptr<void>,QString>( shared_ptr<void>(), QString::null ) );
            
        return loaded_objects;
    }

    if (header.version() == 1 or header.version() == 2)
    {
        //read in the binary data containing all of the objects
        QByteArray compressed_data;
        compressed_data.reserve( RESERVE_SIZE );
    
        if (compressed_data.capacity() != RESERVE_SIZE)
            qWarning() << "Possible memory allocation error!";
    
        ds >> compressed_data;
    
        //validate that the data is correct
        header.assertNotCorrupted(compressed_data);
    
        //uncompress the data
        QByteArray object_data = qUncompress(compressed_data);
    
        QDataStream ds2(object_data);
    
        //version 1 uses Qt 4.2 datastream format
        ds2.setVersion( QDataStream::Qt_4_2 );

        int nobjects = header.dataTypes().count();
        
        #ifdef BOOST_NO_CXX11_SMART_PTR
          std::auto_ptr<SharedDataStream> sds;
        #else
          std::unique_ptr<SharedDataStream> sds;
        #endif        

        if (nobjects > 1)
            sds.reset( new SharedDataStream(ds2) );
        
        for (int i=0; i<nobjects; ++i)
        {
            QString datatype = header.dataTypes().at(i);
        
            //get the type that represents this name
            int id = QMetaType::type( datatype.toLatin1().constData() );

            if ( id == 0 or not QMetaType::isRegistered(id) )
            {
                // check for renamed classes
                QSet<QString> altnames = getAlternativeNames(datatype);
                
                foreach (QString altname, altnames)
                {
                    id = QMetaType::type(altname.toLatin1().constData());
                    
                    if (id != 0 and QMetaType::isRegistered(id))
                    {
                        datatype = altname;
                        break;
                    }
                }
            
                if (id == 0 or not QMetaType::isRegistered(id))
                    throw SireError::unknown_type( QObject::tr(
                        "Cannot deserialise an object of type \"%1\". "
                        "Ensure that the library or module containing "
                        "this type has been loaded and that it has been registered "
                        "with QMetaType.").arg(datatype), CODELOC );
            }
        
            //create a default-constructed object of this type
            shared_ptr<void> ptr( QMetaType::create(id,0), detail::void_deleter(id) );

            if (ptr.get() == 0)
                throw SireError::program_bug( QObject::tr(
                        "Could not create an object of type \"%1\" despite "
                        "this type having been registered with QMetaType. This is "
                        "a program bug!!!").arg(datatype), CODELOC );
    
            //load the object from the datastream
            if ( not QMetaType::load(ds2, id, ptr.get()) )
                throw SireError::program_bug(QObject::tr(
                    "There was an error loading the object of type \"%1\"")
                        .arg(datatype), CODELOC);

            loaded_objects.append( 
                    tuple<shared_ptr<void>,QString>( ptr, datatype ) );
        }
    }
    else
        throw version_error( QObject::tr(
            "Cannot read the object information, as it is written using "
            "the global Sire format %1, while we can only read versions 1 and 2.")
                    .arg(header.version()), CODELOC );


    return loaded_objects;
}

/** This loads an object from the specified file. This binary
    data *must* have been created by the "save" function below. */
QList< tuple<shared_ptr<void>,QString> > SIRESTREAM_EXPORT load(const QString &filename)
{
    QFile f(filename);
    
    if (not f.open( QIODevice::ReadOnly) )
        throw SireError::file_error(f, CODELOC);
        
    QByteArray data = f.readAll();
    
    if (data.isEmpty())
        throw SireError::file_error( QObject::tr(
            "There was an error reading data from the file %1. Either "
            "the file is empty, or some read error has occured.")
                .arg(filename), CODELOC );

    return load(data);
}

/** Return the header for the data */
FileHeader SIRESTREAM_EXPORT getDataHeader(const QByteArray &data)
{
    QDataStream ds(data);
    
    //read the magic
    quint32 magic;
    ds >> magic;
    
    if (magic != SIRE_MAGIC_NUMBER)
        throw version_error( QObject::tr(
            "This data does not appear to have been written by the SireStream::save() "
            "function."), CODELOC );
    
    //read the version and the libraries and versions
    //(die if the required libraries aren't loaded!)
    FileHeader header;
    ds >> header;

    return header;
}

/** Return the header of the file */
FileHeader SIRESTREAM_EXPORT getDataHeader(const QString &filename)
{
    QFile f(filename);
    
    if (not f.open( QIODevice::ReadOnly) )
        throw SireError::file_error(f, CODELOC);

    QDataStream ds(&f);
    
    //read the magic
    quint32 magic;
    ds >> magic;
    
    if (magic != SIRE_MAGIC_NUMBER)
        throw version_error( QObject::tr(
            "This data does not appear to have been written by the SireStream::save() "
            "function."), CODELOC );
    
    //read the version and the libraries and versions
    //(die if the required libraries aren't loaded!)
    FileHeader header;
    ds >> header;

    return header;
}

} //end of namespace SireStream
