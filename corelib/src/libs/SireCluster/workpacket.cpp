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

#include <QTextStream>

#include "workpacket.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>

using namespace SireCluster;
using namespace SireStream;

///////////
/////////// Implementation of WorkPacketBase
///////////

static const RegisterMetaType<WorkPacketBase> r_workbase( MAGIC_ONLY,
                                                          WorkPacketBase::typeName() );

/** Serialise to a binary datastream */
QDataStream SIRECLUSTER_EXPORT &operator<<(QDataStream &ds, 
                                           const WorkPacketBase &workbase)
{
    writeHeader(ds, r_workbase, 1);
    
    SharedDataStream sds(ds);
    
    sds << workbase.current_progress;
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIRECLUSTER_EXPORT &operator>>(QDataStream &ds, WorkPacketBase &workbase)
{
    VersionID v = readHeader(ds, r_workbase);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> workbase.current_progress;
    }
    else
        throw version_error(v, "1", r_workbase, CODELOC);
        
    return ds;
}

/** Constructor */
WorkPacketBase::WorkPacketBase() 
               : QSharedData(), current_progress(0)
{}

/** Copy constructor */
WorkPacketBase::WorkPacketBase(const WorkPacketBase &other)
               : QSharedData(),
                 current_progress(other.current_progress)
{}

/** Destructor */
WorkPacketBase::~WorkPacketBase()
{}

/** Copy assignment operator */
WorkPacketBase& WorkPacketBase::operator=(const WorkPacketBase &other)
{
    if (this != &other)
    {
        current_progress = other.current_progress;
    }
    
    return *this;
}

/** Return whether or not this work packet should be stored 
    as a binary array - this is used by Promise to work out
    how to store the initial WorkPacket state. Only large
    packets should be binary packed (as they are then 
    compressed) */
bool WorkPacketBase::shouldPack() const
{
    return false;
}

/** Return the approximate maximum size (in bytes) of the WorkPacket. This
    doesn't have to exact (or indeed accurate) - it is used
    to help the WorkPacket::pack() function reserve enough
    space when serialising this packet to a binary array. 
    The only penalty of getting this wrong is that you'll
    either allocate too much space, or be reallocating while
    the packet is being written */
int WorkPacketBase::approximatePacketSize() const
{
    // adds 8 to give a little lee-way
    return sizeof(float) + sizeof(bool) + 8;
}

/** Return the current progress of the work (percentage) */
float WorkPacketBase::progress() const
{
    return current_progress;
}

/** Return whether or not this is an Error WorkPacket */
bool WorkPacketBase::isError() const
{
    return false;
}

/** Throw the error, if this is in an error state */
void WorkPacketBase::throwError() const
{}

/** Whether or not the job has been aborted */
bool WorkPacketBase::wasAborted() const
{
    return false;
}

/** Perform one chunk of the calculation - Any exceptions are
    caught in WorkPacket::runChunk, where that are converted
    into an ErrorPacket */
void WorkPacketBase::runChunk()
{
    if (this->wasAborted() or this->hasFinished())
    {
        return;
    }

    float new_progress = this->chunk();
    current_progress = qMin( float(0), qMax(new_progress,float(100)) );
}

///////////
/////////// Implementation of ErrorPacket
///////////

static const RegisterMetaType<ErrorPacket> r_errorpacket;

/** Serialise to a binary datastream */
QDataStream SIRECLUSTER_EXPORT &operator<<(QDataStream &ds, 
                                           const ErrorPacket &errorpacket)
{
    writeHeader(ds, r_errorpacket, 1);
    
    SharedDataStream sds(ds);

    sds << errorpacket.error_data
        << static_cast<const WorkPacketBase&>(errorpacket);
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIRECLUSTER_EXPORT &operator>>(QDataStream &ds, ErrorPacket &errorpacket)
{
    VersionID v = readHeader(ds, r_errorpacket);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> errorpacket.error_data
            >> static_cast<WorkPacketBase&>(errorpacket);
    }
    else
        throw version_error(v, "1", r_errorpacket, CODELOC);
        
    return ds;
}

/** Constructor */
ErrorPacket::ErrorPacket() : WorkPacketBase()
{}

/** Construct an ErrorPacket for the error 'e' */
ErrorPacket::ErrorPacket(const SireError::exception &e)
            : WorkPacketBase()
{
    error_data = e.pack();
}

/** Copy constructor */
ErrorPacket::ErrorPacket(const ErrorPacket &other)
            : WorkPacketBase(other), error_data(other.error_data)
{}

/** Destructor */
ErrorPacket::~ErrorPacket()
{}

/** Copy assignment operator */
ErrorPacket& ErrorPacket::operator=(const ErrorPacket &other)
{
    if (this != &other)
    {
        error_data = other.error_data;
        WorkPacketBase::operator=(other);
    }
    
    return *this;
}

/** Return the approximate maximum size (in bytes) of the WorkPacket. This
    doesn't have to exact (or indeed accurate) - it is used
    to help the WorkPacket::pack() function reserve enough
    space when serialising this packet to a binary array. 
    The only penalty of getting this wrong is that you'll
    either allocate too much space, or be reallocating while
    the packet is being written */
int ErrorPacket::approximatePacketSize() const
{
    return error_data.count() + WorkPacketBase::approximatePacketSize();
}

/** Return whether or not the work has finished */
bool ErrorPacket::hasFinished() const
{
    return true;
}

/** Return whether or not this is an error */
bool ErrorPacket::isError() const
{
    return not error_data.isEmpty();
}

/** Throw the error associated with this packet */
void ErrorPacket::throwError() const
{
    if (not error_data.isEmpty())
    {
        SireError::exception::unpackAndThrow(error_data);
    }
}

/** Perform one chunk of the work, returning the progress
    after the chunk */
float ErrorPacket::chunk()
{
    return 100;
}

///////////
/////////// Implementation of AbortPacket
///////////

static const RegisterMetaType<AbortPacket> r_abortpacket;

/** Serialise to a binary datastream */
QDataStream SIRECLUSTER_EXPORT &operator<<(QDataStream &ds, 
                                           const AbortPacket &abortpacket)
{
    writeHeader(ds, r_abortpacket, 1);
    
    SharedDataStream sds(ds);

    sds << static_cast<const WorkPacketBase&>(abortpacket);
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIRECLUSTER_EXPORT &operator>>(QDataStream &ds, AbortPacket &abortpacket)
{
    VersionID v = readHeader(ds, r_abortpacket);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> static_cast<WorkPacketBase&>(abortpacket);
    }
    else
        throw version_error(v, "1", r_abortpacket, CODELOC);
        
    return ds;
}

/** Constructor */
AbortPacket::AbortPacket() : WorkPacketBase()
{}

/** Copy constructor */
AbortPacket::AbortPacket(const AbortPacket &other) : WorkPacketBase(other)
{}

/** Destructor */
AbortPacket::~AbortPacket()
{}

/** Copy assignment operator */
AbortPacket& AbortPacket::operator=(const AbortPacket &other)
{
    WorkPacketBase::operator=(other);
    return *this;
}

/** Return whether or not the work has finished */
bool AbortPacket::hasFinished() const
{
    return true;
}

/** Return whether or not this was aborted (it obviously was!) */
bool AbortPacket::wasAborted() const
{
    return true;
}

/** Perform one chunk of the work, returning the progress
    after the chunk */
float AbortPacket::chunk()
{
    return 100;
}

///////////
/////////// Implementation of WorkPacket
///////////

static const RegisterMetaType<WorkPacket> r_workpacket(NO_ROOT);

/** Serialise to a binary datastream */
QDataStream SIRECLUSTER_EXPORT &operator<<(QDataStream &ds, 
                                           const WorkPacket &workpacket)
{
    writeHeader(ds, r_workpacket, 1);
    
    SharedDataStream sds(ds);
    
    sds << workpacket.d;
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIRECLUSTER_EXPORT &operator>>(QDataStream &ds, WorkPacket &workpacket)
{
    VersionID v = readHeader(ds, r_workpacket);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> workpacket.d;
    }
    else
        throw version_error(v, "1", r_workpacket, CODELOC);
        
    return ds;
}

/** Create a null work packet */
WorkPacket::WorkPacket()
{}

/** Construct from the passed work object */
WorkPacket::WorkPacket(const WorkPacketBase &work) : d(work)
{}

/** Copy constructor */
WorkPacket::WorkPacket(const WorkPacket &other) : d(other.d)
{}

/** Destructor */
WorkPacket::~WorkPacket()
{}

/** Copy assignment operator */
WorkPacket& WorkPacket::operator=(const WorkPacket &other)
{
    if (other.isNull())
    {
        d = 0;
    }
    else
        d = other.d;
        
    return *this;
}

/** Return whether or not this is the null (empty) work packet */
bool WorkPacket::isNull() const
{
    return d.constData() == 0;
}

/** Return whether or not we should pack this WorkPacket when
    we are storing it. */
bool WorkPacket::shouldPack() const
{
    if (this->isNull())
        return false;
        
    else
        return d->shouldPack();
}

/** Pack this WorkPacket into a binary array */
QByteArray WorkPacket::pack() const
{
    if (this->isNull())
        return QByteArray();

    QByteArray data;
    data.reserve( d->approximatePacketSize() );
    
    QDataStream ds(&data, QIODevice::WriteOnly);
    
    ds << *this;
    
    data = qCompress(data);
    
    return data;
}

/** Unpack a WorkPacket from the passed binary data. This binary
    data *MUST* have been created by WorkPacket::pack() */
WorkPacket WorkPacket::unpack(const QByteArray &data)
{
    if (data.isEmpty())
        return WorkPacket();
        
    WorkPacket workpacket;
    
    QByteArray uncompressed_data = qUncompress(data);
    
    QDataStream ds(uncompressed_data);
    
    ds >> workpacket;
    
    return workpacket;
}

/** Return whether or not this work is in an error state */
bool WorkPacket::isError() const
{
    if (not this->isNull())
    {
        return d->isError();
    }
    else
        return false;
}

/** Throw any error associated with this WorkPacket
    (this does nothing if there is no error) */
void WorkPacket::throwError() const
{
    if (not this->isNull())
    {
        d->throwError();
    }
}

/** Return whether or not the work was aborted */
bool WorkPacket::wasAborted() const
{
    if (not this->isNull())
    {
        return d->wasAborted();
    }
    else
        return false;
}

/** Return whether or not the work has finished (or is in an 
    error state, or was aborted) - essentially, is there any
    more of this work packet to run? */
bool WorkPacket::hasFinished() const
{
    if (not this->isNull())
    {
        return d->wasAborted() or d->hasFinished();
    }
    else
        return true;
}

/** Abort the work */
void WorkPacket::abort()
{
    if (not this->isNull())
    {
        d = AbortPacket();
    }
}

/** This sets the error state - this is done by replacing the 
    existing WorkPacket with an ErrorPacket that describes
    the error */
void WorkPacket::setError(const SireError::exception &e) throw()
{
    try
    {
        d = ErrorPacket(e);
    }
    catch(const SireError::exception &e2)
    {
        d = ErrorPacket(e2);
    }
    catch(...)
    {
        d = ErrorPacket( SireError::unknown_exception( QObject::tr(
                "An unknown error occured while creating an ErrorPacket."),
                    CODELOC ) );
    }
}

/** Run a chunk of work */
void WorkPacket::runChunk() throw()
{
    if (this->isNull())
        return;

    try
    {
        if (d->hasFinished() or d->wasAborted())
            return;
    
        d->runChunk();
    }
    catch(const SireError::exception &e)
    {
        this->setError(e);
    }
    catch(const std::exception &e)
    {
        this->setError( SireError::std_exception(e) );
    }
    catch(...)
    {
        this->setError( SireError::unknown_exception( QObject::tr(
                "There was an unknown exception thrown while running a chunk "
                "of the WorkPacket %1 (progress = %2 %%)")
                    .arg(d->what()).arg(d->progress()), CODELOC ) );
    }
}

/** Return the current progress of the calculation (percentage) */
float WorkPacket::progress() const
{
    if (not this->isNull())
    {
        return d->progress();
    }
    else
        return 100;
}

/** Return a reference to the underlying Worker object */
const WorkPacketBase& WorkPacket::base() const
{
    if (this->isNull())
        throw SireError::nullptr_error( QObject::tr(
            "The null WorkPacket has no base!"), CODELOC );
        
    return *d;
}

///////////
/////////// Implementation of WorkTest
///////////

static const RegisterMetaType<WorkTest> r_worktest;

/** Serialise to a binary datastream */
QDataStream SIRECLUSTER_EXPORT &operator<<(QDataStream &ds, const WorkTest &worktest)
{
    writeHeader(ds, r_worktest, 1);
    
    SharedDataStream sds(ds);

    sds << worktest.current << worktest.start 
        << worktest.end << worktest.step
        << static_cast<const WorkPacketBase&>(worktest);
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIRECLUSTER_EXPORT &operator>>(QDataStream &ds, WorkTest &worktest)
{
    VersionID v = readHeader(ds, r_worktest);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> worktest.current >> worktest.start >> worktest.end >> worktest.step
            >> static_cast<WorkPacketBase&>(worktest);
    }
    else
        throw version_error(v, "1", r_worktest, CODELOC);
        
    return ds;
}

/** Constructor */
WorkTest::WorkTest()
         : WorkPacketBase(), current(0), start(0), end(0), step(0)
{}

/** Construct a work test that counts from start to end in steps of 'step' */
WorkTest::WorkTest(int _start, int _end, int _step)
         : WorkPacketBase(), 
           current(_start), start(_start), end(_end), step(_step)
{}

/** Copy constructor */
WorkTest::WorkTest(const WorkTest &other)
         : WorkPacketBase(other), current(other.current),
           start(other.start), end(other.end), step(other.step)
{}

/** Destructor */
WorkTest::~WorkTest()
{}

/** Copy assignment operator */
WorkTest& WorkTest::operator=(const WorkTest &other)
{
    if (this != &other)
    {
        current = other.current;
        start = other.start;
        end = other.end;
        step = other.step;
        WorkPacketBase::operator=(other);
    }
    
    return *this;
}

/** Return the approximate maximum size (in bytes) of the WorkPacket. This
    doesn't have to exact (or indeed accurate) - it is used
    to help the WorkPacket::pack() function reserve enough
    space when serialising this packet to a binary array. 
    The only penalty of getting this wrong is that you'll
    either allocate too much space, or be reallocating while
    the packet is being written */
int WorkTest::approximatePacketSize() const
{
    return 4 * sizeof(qint32) + WorkPacketBase::approximatePacketSize();
}

/** Return whether or not the work has finished */
bool WorkTest::hasFinished() const
{
    return current == end;
}

/** Perform one chunk of the work, returning the progress
    after the chunk */
float WorkTest::chunk()
{
    if (step == 0)
        throw SireError::invalid_arg( QObject::tr(
                "You cannot use a step size of zero!"), CODELOC );

    if (start < end)
    {
        if (step < 0)
            throw SireError::invalid_arg( QObject::tr(
                "You cannot use a negative step size if start is less than end!"),
                    CODELOC );
    
        current = qMin( current+step, end );
        
        QTextStream ts(stdout);
        ts << "I've counted to " << current << "\n";
        sleep(1);
        
        return 100.0 - ( 100.0 * double(end - current) / double(end - start) );
    }
    else
    {
        if (step > 0)
            throw SireError::invalid_arg( QObject::tr(
                "You cannot use a positive step size if start is greater than end!"),
                    CODELOC );
    
        current = qMax( current+step, end );

        QTextStream ts(stdout);
        ts << "I've counted to " << current << "\n";
        sleep(1);
        
        return 100.0 - ( 100.0 * double(end - current) / double(end - start) );
    }
}

const char* WorkPacket::typeName()
{
    return QMetaType::typeName( qMetaTypeId<WorkPacket>() );
}

const char* ErrorPacket::typeName()
{
    return QMetaType::typeName( qMetaTypeId<ErrorPacket>() );
}

const char* AbortPacket::typeName()
{
    return QMetaType::typeName( qMetaTypeId<AbortPacket>() );
}

const char* WorkTest::typeName()
{
    return QMetaType::typeName( qMetaTypeId<WorkTest>() );
}

WorkTest* WorkTest::clone() const
{
    return new WorkTest(*this);
}


AbortPacket* AbortPacket::clone() const
{
    return new AbortPacket(*this);
}


ErrorPacket* ErrorPacket::clone() const
{
    return new ErrorPacket(*this);
}

