/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2009  Christopher Woods
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

#include <QThread>
#include <QMutex>

#include <QFile>
#include <QTextStream>

#include "sire_config.h"

#include "meminfo.h"

#include <stdio.h>

#include <QDebug>

using namespace SireBase;

//////////////////////////////////////////////////
///
///  Implementation of SysInfoPvt
///    This is a platform dependent class
///    that is used to get information about
///    the amount of memory available to the system
///
//////////////////////////////////////////////////

#ifdef HAVE_SYSCTL

    ///////////
    /////////// Implementation for systems that have sysctl support
    /////////// for probing the status of the system
    ///////////
    namespace SireBase{ namespace detail {

        class SysInfoPvt
        {
        public:
            SysInfoPvt()
            {}

            ~SysInfoPvt()
            {}

            quint64 totalSystemMemory() const
            {
                return 0;
            }

            quint64 totalVirtualMemory() const
            {
                return 0;
            }

            quint64 usedSystemMemory() const
            {
                return 0;
            }

            quint64 usedVirtualMemory() const
            {
                return 0;
            }

            QString toString() const
            {
                return QString();
            }
        };

    }}

#else

    ///////////
    /////////// Implementation for systems with no memory info
    ///////////
    namespace SireBase{ namespace detail {

        class SysInfoPvt
        {
        public:
            SysInfoPvt()
            {}

            ~SysInfoPvt()
            {}

            quint64 totalSystemMemory() const
            {
                return 0;
            }

            quint64 totalVirtualMemory() const
            {
                return 0;
            }

            quint64 usedSystemMemory() const
            {
                return 0;
            }

            quint64 usedVirtualMemory() const
            {
                return 0;
            }

            QString toString() const
            {
                return QObject::tr("System memory: Unknown");
            }
        };

    }}

#endif

//////////////////////////////////////////////////
///
///  Implementation of MemInfoPvt
///    This is a platform dependent class
///    that is used to get information about
///    memory allocations via the malloc subsystem
///
//////////////////////////////////////////////////

#ifdef HAVE_MALLINFO

    ///////////
    /////////// Implementation for systems with mallinfo
    ///////////  (e.g. glibc systems)
    ///////////

    #ifdef HAVE_MALLOC_MALLOC_H
        #include <malloc/malloc.h>  // CONDITIONAL_INCLUDE
    #else
        #include <malloc.h>  // CONDITIONAL_INCLUDE
    #endif

    namespace SireBase{ namespace detail
    {
        /** Implementation of MemInfoPvt for mallinfo systems */
        class MemInfoPvt : public SysInfoPvt
        {
        public:
            MemInfoPvt() : SysInfoPvt()
            {
                minfo = mallinfo();
            }

            ~MemInfoPvt()
            {
            }

            struct ::mallinfo minfo;

            quint64 allocatedBytes() const
            {
                quint64 arena = minfo.arena;
                quint64 hblkhd = minfo.hblkhd;

                return arena + hblkhd;
            }

            quint64 mMappedBytes() const
            {
                return minfo.hblkhd;
            }

            quint64 usedBytes() const
            {
                qint64 uordblks = minfo.uordblks;
                qint64 hblkhd = minfo.hblkhd;

                return uordblks + hblkhd;
            }

            QString toString() const
            {
                float alloc = this->allocatedBytes() / (1024.0*1024.0);
                float used = this->usedBytes() / (1024.0*1024.0);

                QString this_string = QObject::tr("Memory usage: %1 MB allocated, "
                                                  "of which %2 MB are used (%3 \%)")
                                            .arg(alloc).arg(used)
                                            .arg(100.0*used/alloc);

                QString sys_string = SysInfoPvt::toString();

                if (sys_string.isEmpty())
                    return this_string;
                else
                    return QString("%1\n%2").arg(this_string, sys_string);
            }
        };

    }}

#elif HAVE_MSTATS

    ///////////
    /////////// Implementation for systems with mstats (e.g. OS X)
    ///////////

    #ifdef HAVE_MALLOC_MALLOC_H
        #include <malloc/malloc.h>  // CONDITIONAL_INCLUDE
    #else
        #include <malloc.h>  // CONDITIONAL_INCLUDE
    #endif

    namespace SireBase{ namespace detail
    {
        /** Implementation of MemInfoPvt for mallinfo systems */
        class MemInfoPvt : public SysInfoPvt
        {
        public:
            MemInfoPvt() : SysInfoPvt()
            {
                minfo = mstats();
            }

            ~MemInfoPvt()
            {
            }

            struct ::mstats minfo;

            quint64 allocatedBytes() const
            {
                return minfo.bytes_total;
            }

            quint64 mMappedBytes() const
            {
                return 0;
            }

            quint64 usedBytes() const
            {
                return minfo.bytes_used;
            }

            QString toString() const
            {
                float alloc = this->allocatedBytes() / (1024.0*1024.0);
                float used = this->usedBytes() / (1024.0*1024.0);

                QString this_string = QObject::tr("Memory usage: %1 MB allocated, "
                                                  "of which %2 MB are used (%3 \%)")
                                            .arg(alloc).arg(used)
                                            .arg(100.0*used/alloc);

                QString sys_string = SysInfoPvt::toString();

                if (sys_string.isEmpty())
                    return this_string;
                else
                    return QString("%1\n%2").arg(this_string, sys_string);
            }
        };

    }}

#else

    ///////////
    /////////// Implementation for systems with no support
    ///////////

    namespace SireBase{ namespace detail
    {

        /** Implementation of a null MemInfoPvt */
        class MemInfoPvt : public SysInfoPvt
        {
        public:
            MemInfoPvt() : SysInfoPvt()
            {}

            ~MemInfoPvt()
            {}

            QString toString() const
            {
                return QObject::tr("Memory usage: unknown\n%1")
                                    .arg( SysInfoPvt::toString() );
            }

            quint64 allocatedBytes() const
            {
                return 0;
            }

            quint64 mMappedBytes() const
            {
                return 0;
            }

            quint64 usedBytes() const
            {
                return 0;
            }
        };

    }}

#endif

///////////
/////////// Implementation of MemInfo
///////////

/** Null constructor */
MemInfo::MemInfo()
{}

/** Copy constructor */
MemInfo::MemInfo(const MemInfo &other) : d(other.d)
{}

/** Destructor */
MemInfo::~MemInfo()
{}

/** Copy assignment operator */
MemInfo& MemInfo::operator=(const MemInfo &other)
{
    d = other.d;
    return *this;
}

/** Return a string containing details of the memory usage. This is
    a platform dependent string, but can be useful to print to the user */
QString MemInfo::toString() const
{
    if (d.get() == 0)
        return QObject::tr("Memory usage not measured. Use MemInfo::takeMeasurement() "
                           "to take a reading.");

    else
        return d->toString();
}

/** Return the total number of bytes allocated to this process by
    memory subsystem (only on the heap - this ignores the stack and
    any memory allocated within libraries). Note that fragmentation
    may mean that not all of this memory is in use. */
quint64 MemInfo::allocatedBytes() const
{
    if (d.get() == 0)
        return 0;

    else
        return d->allocatedBytes();
}

/** Return the total number of bytes that have been allocated via
    mmap to this process by the memory subsystem. This total is included
    within "allocatedBytes", therefore "allocatedBytes() - mMappedBytes()
    will return the number of bytes that have not been allocated
    via mmap */
quint64 MemInfo::mMappedBytes() const
{
    if (d.get() == 0)
        return 0;
    else
        return d->mMappedBytes();
}

/** Return the total number of bytes in use by the program. This number
    will be less than or equal to the number of bytes returned by "allocatedBytes".
    This will be less if memory fragmentation means that pages of memory are
    not able to be returned to the operating system when they are freed
    by the program */
quint64 MemInfo::usedBytes() const
{
    if (d.get() == 0)
        return 0;
    else
        return d->usedBytes();
}

/** Take a measurement of the current memory usage of the process */
MemInfo MemInfo::takeMeasurement()
{
    MemInfo ret;

    ret.d.reset( new detail::MemInfoPvt() );

    return ret;
}

/** Return the total number of bytes of system memory
    (actual memory, not virtual memory). To ensure we don't page,
    we have to ensure that this processes memory usage is below
    the total system memory */
quint64 MemInfo::totalSystemMemory() const
{
    if (d.get() == 0)
        return 0;

    else
        return d->totalSystemMemory();
}

/** Return the total number of bytes of virtual memory available on this system
    (this includes system memory) */
quint64 MemInfo::totalVirtualMemory() const
{
    if (d.get() == 0)
        return 0;

    else
        return d->totalVirtualMemory();
}

/** Return the total amount of system memory in use.
    (actual memory, not virtual memory) */
quint64 MemInfo::usedSystemMemory() const
{
    if (d.get() == 0)
        return 0;

    else
        return d->usedSystemMemory();
}

/** Return the total amount of virtual memory in use.
    (includes system memory) */
quint64 MemInfo::usedVirtualMemory() const
{
    if (d.get() == 0)
        return 0;

    else
        return d->usedVirtualMemory();
}

class MemoryMonitor : public QThread
{
public:
    MemoryMonitor(int ms);
    MemoryMonitor(const QString &filename, int ms);

    ~MemoryMonitor();

    void stop();

protected:
    void run();

private:
    QString filename;
    int ms;
};

MemoryMonitor::MemoryMonitor(int time_ms)
              : QThread(), ms(time_ms)
{}

MemoryMonitor::MemoryMonitor(const QString &file, int time_ms)
              : QThread(), filename(file), ms(time_ms)
{}

MemoryMonitor::~MemoryMonitor()
{
    this->stop();
    this->wait();
}

void MemoryMonitor::stop()
{
    ms = -1;
}

static void openFile(const QString &filename, QFile &file)
{
    if (filename.isEmpty())
    {
        if (not file.open( stdout, QIODevice::WriteOnly | QIODevice::Unbuffered ))
            qWarning() << "Could not open STDOUT???";
    }
    else
    {
        if (not file.open( QIODevice::WriteOnly|QIODevice::Unbuffered|QIODevice::Append ))
        {
            openFile(QString(), file);
        }
    }
}

void MemoryMonitor::run()
{
    if (ms == -1)
        return;

    QFile file(filename);

    openFile(filename, file);

    QTextStream ts( &file );

    ts << " ** " << MemInfo::takeMeasurement().toString() << "\n";
    ts.flush();

    while (ms > 0)
    {
        int have_slept = 0;

        while (have_slept < ms)
        {
            QThread::msleep( qMin(100, ms - have_slept) );
            have_slept += 100;
        }

        if (ms > 0)
        {
            ts << " ** " << MemInfo::takeMeasurement().toString() << "\n";
            ts.flush();
        }
    }
}

Q_GLOBAL_STATIC( QMutex, memMonitorMutex );

static boost::shared_ptr<MemoryMonitor> monitor;

/** Start a monitor that prints the memory usage of the program
    out to the screen every 'ms' milliseconds */
void MemInfo::startMonitoring(int ms)
{
    QMutexLocker lkr( memMonitorMutex() );

    monitor.reset( new MemoryMonitor(ms) );

    monitor->start();
}

/** Start a monitor that prints the memory usage of the program
    out to the file 'filename' every 'ms' milliseconds */
void MemInfo::startMonitoring(const QString &filename, int ms)
{
    QMutexLocker lkr( memMonitorMutex() );

    monitor.reset( new MemoryMonitor(filename,ms) );

    monitor->start();
}

/** Stop monitoring the memory of this process */
void MemInfo::stopMonitoring()
{
    QMutexLocker lkr( memMonitorMutex() );

    monitor.reset();
}
