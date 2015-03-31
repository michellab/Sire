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

#include "frontend.h"
#include "backend.h"
#include "workpacket.h"

#include <QDebug>

using namespace SireCluster;

using boost::shared_ptr;

/////////
///////// Implementation of FrontendBase
/////////

/** Constructor */
FrontendBase::FrontendBase() : boost::noncopyable(), datamutex( QMutex::Recursive )
{}

/** Destructor */
FrontendBase::~FrontendBase()
{}

/////////
///////// Implementation of LocalFrontend
/////////

/** This private class implements the local Frontend that
    is used to communicate with a Backend that is in the same
    address space 
    
    @author Christopher Woods
*/
class LocalFrontend : public FrontendBase
{

public:
    LocalFrontend(const ActiveBackend &backend);
    ~LocalFrontend();
    
    bool isLocal() const;
    
    QUuid UID();
    
    void startJob(const WorkPacket &workpacket);
    
    void stopJob();
    void abortJob();
    
    void wait();
    bool wait(int timeout);
    
    float progress();
    WorkPacket interimResult();
    
    WorkPacket result();

private:
    /** The active backend */
    ActiveBackend backend;
};

/** Construct a Frontend for the local Backend 'backend' */
LocalFrontend::LocalFrontend(const ActiveBackend &_backend)
              : FrontendBase(), backend(_backend)
{}

/** Destructor */
LocalFrontend::~LocalFrontend()
{}

bool LocalFrontend::isLocal() const
{
    return true;
}

QUuid LocalFrontend::UID()
{
    return backend.UID();
}

void LocalFrontend::startJob(const WorkPacket &workpacket)
{
    backend.startJob(workpacket);
}

void LocalFrontend::stopJob()
{
    backend.stopJob();
}

void LocalFrontend::abortJob()
{
    backend.abortJob();
}

void LocalFrontend::wait()
{
    backend.wait();
}

bool LocalFrontend::wait(int timeout)
{
    return backend.wait(timeout);
}

float LocalFrontend::progress()
{
    return backend.progress();
}

WorkPacket LocalFrontend::interimResult()
{
    return backend.interimResult();
}

WorkPacket LocalFrontend::result()
{
    return backend.result();
}

/////////
///////// Implementation of Frontend
/////////

/** Null constructor */
Frontend::Frontend()
{}

/** Construct from the passed Frontend pointer */
Frontend::Frontend(const boost::shared_ptr<FrontendBase> &ptr)
         : d(ptr)
{}

/** Construct a local Frontend that talks to the local Backend 'backend'.
    This will block while the backend is busy talking to another frontend */
Frontend::Frontend(const Backend &backend)
{
    if (not backend.isNull())
    {
        d.reset( new LocalFrontend(backend.connect()) );
    }
}

/** Construct a local Frontend that talks to the local Backend. This only
    tries to make a connection - if the backend is busy then this gives
    up and a null Frontend is returned */
Frontend Frontend::tryAcquire(const Backend &backend)
{
    ActiveBackend active_backend = backend.tryConnect();
    
    Frontend frontend;
    
    if (not active_backend.isNull())
    {
        frontend.d.reset( new LocalFrontend(active_backend) );
    }
    
    return frontend;
}

/** Copy constructor */
Frontend::Frontend(const Frontend &other)
         : d(other.d)
{}

/** Destructor */
Frontend::~Frontend()
{}

/** Copy assignment operator */
Frontend& Frontend::operator=(const Frontend &other)
{
    d = other.d;
    return *this;
}

/** Comparison operator */
bool Frontend::operator==(const Frontend &other) const
{
    return d.get() == other.d.get();
}

/** Comparison operator */
bool Frontend::operator!=(const Frontend &other) const
{
    return d.get() != other.d.get();
}

/** Return whether or not the Backend is local (running
    in the same address space as the Frontend) */
bool Frontend::isLocal() const
{
    if (this->isNull())
        return false;
        
    else
        return d->isLocal();
}

/** Return whether or not this is a null frontend */
bool Frontend::isNull() const
{
    return d.get() == 0;
}

/** Return the UID of the backend */
QUuid Frontend::UID()
{
    if (not this->isNull())
    {
        QMutexLocker lkr( &(d->datamutex) );
        return d->UID();
    }
    else
        return QUuid();
}

/** Perform the work in 'workpacket' on the backend - this 
    blocks until the work has started */
void Frontend::startJob(const WorkPacket &workpacket)
{
    if (not this->isNull())
    {
        QMutexLocker lkr( &(d->datamutex) );
        d->startJob(workpacket);
    }
}

/** Stop the job running on the backend, and return the
    workpacket in the state it is now at now that the job
    has stopped */
void Frontend::stopJob()
{
    if (not this->isNull())
    {
        QMutexLocker lkr( &(d->datamutex) );
        d->stopJob();
    }
}

/** Abort the job running on the backend and return the
    state of the work once it has been aborted */
void Frontend::abortJob()
{
    if (not this->isNull())
    {
        QMutexLocker lkr( &(d->datamutex) );
        d->abortJob();
    }
}

/** Wait until the backend has finished processing the work */
void Frontend::wait()
{
    if (not this->isNull())
    {
        d->wait();
    }
}

/** Wait until the backend has finished processing the work, or
    until 'timeout' milliseconds have passed - this returns
    whether or not the job has finished */
bool Frontend::wait(int timeout)
{
    if (not this->isNull())
    {
        return d->wait(timeout);
    }
    else
        return true;
}

/** Return the current progress of the job */
float Frontend::progress()
{
    if (not this->isNull())
    {
        QMutexLocker lkr( &(d->datamutex) );
        return d->progress(); 
    }
    else
        return 0;
}

/** Return the work as it is at the moment */
WorkPacket Frontend::interimResult()
{
    if (not this->isNull())
    {
        QMutexLocker lkr( &(d->datamutex) );
        return d->interimResult();
    }
    else
        return WorkPacket();
}

/** Return the final result of the calculation */
WorkPacket Frontend::result()
{
    if (not this->isNull())
    {
        QMutexLocker lkr( &(d->datamutex) );
        return d->result();
    }
    else
        return WorkPacket();
}
