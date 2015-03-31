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

#ifdef SIRE_USE_MPI

#include <mpi.h>

#include <QTime>

#include "mpifrontend.h"
#include "mpicluster.h"

#include "SireCluster/workpacket.h"

#include "SireError/printerror.h"
#include "SireError/getbacktrace.h"

#include <QDebug>

using namespace SireCluster;
using namespace SireCluster::MPI;

/** Null constructor */
MPIFrontend::MPIFrontend() : FrontendBase()
{}

/** Construct a frontend that uses the passed  
    point-to-point communicator to communicate with 
    the backend */
MPIFrontend::MPIFrontend(const P2PComm &p2pcomm)
            : FrontendBase(), p2p(p2pcomm)
{}

/** Destructor */
MPIFrontend::~MPIFrontend()
{}

/** This is only local if the p2p communicator is local */
bool MPIFrontend::isLocal() const
{
    return const_cast<P2PComm*>(&p2p)->isLocal();
}

/** Return the UID of the backend */
QUuid MPIFrontend::UID()
{
    if (p2p.isNull())
        return QUuid();
    
    else if (cached_uid.isNull())
    {
        QMutexLocker lkr(&datamutex);

        p2p.sendMessage( P2PComm::GETUID );
        cached_uid = p2p.awaitResponse<QUuid>(true);
    }
    
    return cached_uid;
}

/** Start a job on the backend */
void MPIFrontend::startJob(const WorkPacket &workpacket)
{
    if (not p2p.isNull())
    {
        QMutexLocker lkr(&datamutex);

        p2p.sendMessage( P2PComm::START, workpacket );
        
        int result = p2p.awaitIntegerResponse(true);
        
        if (result != 0)
            qDebug() << SireError::getPIDString()
                     << "Starting a remote job got a weird response" << result;
    }
}

/** Stop the job on the backend */
void MPIFrontend::stopJob()
{
    if (not p2p.isNull())
    {
        QMutexLocker lkr(&datamutex);
        
        p2p.sendMessage( P2PComm::STOP );

        int result = p2p.awaitIntegerResponse();
        
        if (result != 0)
            qDebug() << SireError::getPIDString()
                     << "Stopping a remote job got a weird response" << result;
    }
}

/** Abort the job on the backend */
void MPIFrontend::abortJob()
{
    if (not p2p.isNull())
    {
        QMutexLocker lkr(&datamutex);
        
        p2p.sendMessage( P2PComm::ABORT );

        int result = p2p.awaitIntegerResponse();
        
        if (result != 0)
            qDebug() << SireError::getPIDString() 
                     << "Aborting a remote job got a weird response" << result;
    }
}

/** Wait for the job to finish */
void MPIFrontend::wait()
{
    if (p2p.isNull())
    {
        return;
    }
    
    QMutexLocker lkr( &datamutex );
    
    while (true)
    {
        p2p.sendMessage( P2PComm::IS_RUNNING );
        
        if (not p2p.awaitIntegerResponse())
            break;
            
        else
            //wait a second
            ::sleep(1);
    }
}

/** Wait for the job to finish, or until 'timeout'
    milliseconds have passed */
bool MPIFrontend::wait(int timeout)
{
    QTime t;
    t.start();

    if (timeout < 0)
    {
        this->wait();
        return true;
    }

    if (p2p.isNull())
    {
        return true;
    }
    
    QMutexLocker lkr( &datamutex );
    
    while (t.elapsed() < timeout)
    {
        p2p.sendMessage( P2PComm::IS_RUNNING );
        
        if (not p2p.awaitIntegerResponse())
            return true;
            
        //wait a second
        if (t.elapsed() + 1000 > timeout)
            break;
        
        ::sleep(1);
    }
    
    return false;
}

/** Return the progress of the work */
float MPIFrontend::progress()
{
    if (not p2p.isNull())
    {
        QMutexLocker lkr(&datamutex);
        
        p2p.sendMessage( P2PComm::PROGRESS );

        return p2p.awaitFloatResponse();
    }
    else
        return 0;
}

/** Return an interim result */
WorkPacket MPIFrontend::interimResult()
{
    if (not p2p.isNull())
    {
        QMutexLocker lkr(&datamutex);
        
        p2p.sendMessage( P2PComm::INTERIM );

        return p2p.awaitResponse<WorkPacket>();
    }
    else
        return WorkPacket();
}

/** Return the final result - this blocks until
    it is available */
WorkPacket MPIFrontend::result()
{
    if (not p2p.isNull())
    {
        QMutexLocker lkr(&datamutex);
        
        p2p.sendMessage( P2PComm::RESULT );

        return p2p.awaitResponse<WorkPacket>();
    }
    else
        return WorkPacket();
}

#endif // SIRE_USE_MPI
