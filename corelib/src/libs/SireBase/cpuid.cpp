/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2013  Christopher Woods
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

#include "cpuid.h"

#if defined(_WIN32)
    #define WIN32_LEAN_AND_MEAN
    #include <Windows.h>    // CONDITIONAL_INCLUDE
#elif defined(__APPLE__) || defined(__FreeBSD__)
    #include <sys/param.h>  // CONDITIONAL_INCLUDE
    #include <sys/sysctl.h> // CONDITIONAL_INCLUDE
#else
    #include <sys/sysinfo.h>    // CONDITIONAL_INCLUDE
#endif

#include "SireStream/shareddatastream.h"
#include "SireError/errors.h"

#include "tostring.h"

using namespace SireBase;
using namespace SireStream;
using namespace SireError;

static const RegisterMetaType<CPUID> r_cpuid;

QDataStream &operator<<(QDataStream &ds, const CPUID &cpuid)
{
    writeHeader(ds, r_cpuid, 1);

    SharedDataStream sds(ds);
    sds << cpuid.props;
    return ds;
}

QDataStream &operator>>(QDataStream &ds, CPUID &cpuid)
{
    VersionID v = readHeader(ds, r_cpuid);

    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> cpuid.props;
    }
    else
        throw SireStream::version_error(v, "1", r_cpuid, CODELOC);

    return ds;
}

static QHash<QString,QString> getCPUInfo()
{
    QHash<QString,QString> info;

    #if defined(__APPLE__) || defined(__FreeBSD__)
        char buffer[1024];
        size_t size=sizeof(buffer);
        int err = 0;

        err = sysctlbyname("machdep.cpu.brand_string", &buffer, &size, NULL, 0);
        if (err == 0)
        {
            info["CPU"] = QString(buffer);
        }

        uint64_t val = 0;
        size = sizeof(val);

        err = sysctlbyname("hw.cpufrequency", &val, &size, NULL, 0);
        if (err == 0)
        {
            info["CPU.frequency"] = QString("%1 GHz").arg(val / 1000000000);
        }

        err = sysctlbyname("hw.ncpu", &val, &size, NULL, 0);
        if (err == 0)
        {
            info["CPU.num_cores"] = QString("%1").arg(val);
        }

        err = sysctlbyname("hw.optional.arm64", &val, &size, NULL, 0);
        if (err == 0)
        {
            if (val == 1)
            {
                info["CPU.architecture"] = "ARM64";
            }
            else
            {
                info["CPU.architecture"] = "x86-64";
            }
        }

        err = sysctlbyname("hw.memsize", &val, &size, NULL, 0);
        if (err == 0)
        {
            info["memory"] = QString("%1 GB").arg(val / (1024*1024*1024));
        }
    #endif

    return info;
}

/** Return the list of all searchable supportable features */
QStringList CPUID::supportableFeatures() const
{
    return QStringList();
}

QHash<QString,QString>* CPUID::global_props = 0;

Q_GLOBAL_STATIC( QMutex, globalMutex )

QHash<QString,QString>* CPUID::getCPUID()
{
    if (not global_props)
    {
        QMutexLocker lkr( globalMutex() );

        global_props = new QHash<QString,QString>(getCPUInfo());
    }

    return global_props;
}

/** Constructor */
CPUID::CPUID() : ConcreteProperty<CPUID,Property>()
{
    props = *(getCPUID());
}

/** Copy constructor */
CPUID::CPUID(const CPUID &other) : ConcreteProperty<CPUID,Property>(), props(other.props)
{}

/** Destructor */
CPUID::~CPUID()
{}

/** Copy assignment operator */
CPUID& CPUID::operator=(const CPUID &other)
{
    props = other.props;
    return *this;
}

/** Comparison operator */
bool CPUID::operator==(const CPUID &other) const
{
    return props == other.props;
}

/** Comparison operator */
bool CPUID::operator!=(const CPUID &other) const
{
    return not CPUID::operator==(other);
}

CPUID* CPUID::clone() const
{
    return new CPUID(*this);
}

const char* CPUID::what() const
{
    return CPUID::typeName();
}

const char* CPUID::typeName()
{
    return QMetaType::typeName( qMetaTypeId<CPUID>() );
}

QString CPUID::toString() const
{
    return QObject::tr( "CPUID( %1 )" ).arg( Sire::toString(props) );
}

/** Returns whether or not the CPU supports the passed feature.
    Note that the passed feature must be one of the strings
    as returned by "supportableFeatures" */
bool CPUID::supports(const QString &feature) const
{
    return props.value(feature, "false") == "true";
}

/** Return the list of all features supported on this CPU */
QStringList CPUID::supportedFeatures() const
{
    QStringList supported;

    foreach( QString feature, supportableFeatures() )
    {
        if (supports(feature))
            supported.append(feature);
    }

    return supported;
}

/** Return the Vendor string for this CPU */
QString CPUID::vendor() const
{
    return props.value("vendor", "Unknown");
}

/** Return the Brand string for this CPU */
QString CPUID::brand() const
{
    return props.value("brand", "Unknown");
}

/** Return the clockspeed of this processor. A value of -1 is returned
    if this is not known */
int CPUID::clockSpeed() const
{
    return props.value("cpu_clock", "-1").toInt();
}

/** Return the number of cores of this processor. A value of 1 is returned
    if this is not known (as we must have at least 1 core!) */
int CPUID::numCores() const
{
    #if defined(_WIN32)
    SYSTEM_INFO system_info;
    GetSystemInfo(&system_info);
    return system_info.dwNumberOfProcessors;
    #elif defined(__APPLE__) || defined(__FreeBSD__)
    int mib[2] = { CTL_HW, HW_NCPU };
    int n_proc = 1;
    size_t len = sizeof(n_proc);
    sysctl(mib, 2, &n_proc, &len, NULL, 0);
    return n_proc;
    #else
    return get_nprocs();
    #endif
}

/** Return whether or not this processor supports SSE2 vector instructions */
bool CPUID::supportsSSE2() const
{
    return supports("sse2");
}

/** Return whether or not this processor supports AVX vector instructions */
bool CPUID::supportsAVX() const
{
    return supports("avx");
}
