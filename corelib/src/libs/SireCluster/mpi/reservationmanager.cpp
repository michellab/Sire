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

#include <mpi.h> // needed first by mpich

#include "SireCluster/cluster.h"
#include "reservationmanager.h"
#include "mpicluster.h"
#include "reply.h"
#include "p2pcomm.h"
#include "mpifrontend.h"

#include "SireCluster/frontend.h"

#include "SireMaths/rangenerator.h"

#include "SireError/errors.h"
#include "SireError/printerror.h"

#include <QDebug>

using namespace SireCluster;
using namespace SireCluster::MPI;
using namespace SireCluster::MPI::Messages;

using boost::tuple;
using boost::shared_ptr;

Q_GLOBAL_STATIC( ReservationManager, reservationManager );

/** Constructor */
ReservationManager::ReservationManager() : keep_going(false)
{}

/** Destructor */
ReservationManager::~ReservationManager()
{
    this->wait();
}

/** Start the reservation manager */
void ReservationManager::start()
{
    ReservationManager *d = reservationManager();

    if (MPICluster::isMaster())
    {
        QMutexLocker lkr( &(d->datamutex) );
    
        if (d->keep_going)
            //it is already running
            return;
        
        d->keep_going = true;
        
        static_cast<QThread*>(d)->start();
    }
}

/** Release all backends that have been reserved */
void ReservationManager::freeReservedBackends()
{
    reserved_backends.clear();
}

/** Shutdown the reservation manager */
void ReservationManager::shutdown()
{
    ReservationManager *d = reservationManager();

    if ( d->isRunning() )
    {
        QMutexLocker lkr( &(d->datamutex) );

        d->keep_going = false;
    
        d->waiter.wakeAll();

        lkr.unlock();

        d->wait();
    }
    
    //now free all of the backends that have been reserved
    d->freeReservedBackends();
}

/** Process the request 'request' - this sends out the request
    for availability and then processes the results */
void ReservationManager::processRequest(ReserveBackend request)
{
    QMutexLocker lkr( &reservation_mutex );

    //ask all of the processes to tell us what they have
    //available that matches this request
    RequestAvailability request_available(request);
    
    //create space for a reply to this broadcast
    Reply reply(request_available);
    
    //now broadcast this message to all nodes
    MPICluster::send(request_available);
    
    //this contains the list of available backends
    QList< tuple<int,QUuid> > backend_uids;

    //get all of the replies
    QHash<int,ReplyValue> replies = reply.replies();

    for (QHash<int,ReplyValue>::const_iterator it = replies.constBegin();
         it != replies.constEnd(); 
         ++it)
    {
        if (not it->isError())
        {
            QList<QUuid> available_backends = it->asA< QList<QUuid> >();
                                                        
            foreach (QUuid available_backend, available_backends)
            {
                backend_uids.append( tuple<int,QUuid>(it.key(),available_backend) );
            }
        }
    }
    
    //ok - we now have the list - work out if it is what we want
    Message message;
    
    if (request.requestedUID().isNull())
    {
        //we have requested 'n' nodes
        QList< tuple<int,QUuid> > reserve_uids;
        
        int navailable = qMin(backend_uids.count(), request.nBackends());
        
        for (int i=0; i<navailable; ++i)
        {
            reserve_uids.append( backend_uids[i] );
        }
        
        message = Reservation(reserve_uids, request);
    }
    else
    {
        if (backend_uids.count() == 1)
            message = Reservation(backend_uids, request);
        
        else    
            message = Reservation( QList< tuple<int,QUuid> >(), request );
    }
        
    Reply reservation_reply(message);

    //broadcast the reservation to the cluster - this allows nodes
    //that are not involved to remove themselves from the point-to-point
    //communicator, and to also tell them that they can free the nodes
    //that were temporarily reserved
    MPICluster::send(message);
        
    //wait for all of the nodes to have responded to the request
    reservation_reply.wait();

    //send a reply back to the original process that requested the
    //nodes, that contains the list of node UIDs that have connected
    replies = reservation_reply.replies();
    
    QList<QUuid> connected_backends;
    
    for (QHash<int,ReplyValue>::const_iterator it = replies.constBegin();
         it != replies.constEnd();
         ++it)
    {
        if (not it->isError())
        {
            connected_backends += it->asA< QList<QUuid> >();
        }
    }

    MPICluster::send( Result(request, connected_backends) );
}

/** The event loop that runs in a background thread - this 
    manages the reservations */
void ReservationManager::run()
{
    SireError::setThreadString( "ReservationManager" );
    SireMaths::seed_qrand();

    QMutexLocker lkr(&datamutex);
    
    while (keep_going)
    {
        if (not request_queue.isEmpty())
        {
            //we have a request to process
            ReserveBackend request = request_queue.dequeue();
            lkr.unlock();
            
            this->processRequest(request);
            
            lkr.relock();
        }
        else if (keep_going)
            waiter.wait( &datamutex, 50 );
    }
    
    //send a response back to any remaining requests saying that
    //no backends are available
    while (not request_queue.isEmpty())
    {
        ReserveBackend request = request_queue.dequeue();
        MPICluster::send( Result(request, QList<QUuid>()) );
    }
}

/** Process the request in 'request' */
void ReservationManager::reserveBackends( const ReserveBackend &request )
{
    if (not MPICluster::isMaster())
        throw SireError::program_bug( QObject::tr(
                "Only the master MPI process can read a ReserveBackend message."),
                    CODELOC );

    //is there a backend with the requested UID?
    if (not request.requestedUID().isNull())
    {
        if (not MPICluster::hasBackend(request.requestedUID()))
            throw SireError::unavailable_resource( QObject::tr(
                "There is no backend available with UID %1.")
                    .arg(request.requestedUID().toString()), CODELOC );
    }
      
    //queue this request
    ReservationManager *d = reservationManager();
    
    QMutexLocker lkr( &(d->datamutex) );
    d->request_queue.enqueue(request);

    d->waiter.wakeAll();
}

/** Reserve all available backends that match the request 'request'
    and return the UIDs of all matching reserved backends */                           
QList<QUuid> ReservationManager::findAvailable( const ReserveBackend &request )
{
    QList<Frontend> available_frontends = Cluster::localBackends();
    
    if (not request.requestedUID().isNull())
    {
        QMutableListIterator<Frontend> it(available_frontends);
        
        while (it.hasNext())
        {
            if (it.value().UID() != request.requestedUID())
                it.remove();
        }
    }

    if (available_frontends.isEmpty())
        return QList<QUuid>();
    else
    {
        ReservationManager *d = reservationManager();
    
        //we've found some frontends to connect to - save the 
        //reservation and return their UIDs
        QMutexLocker lkr( &(d->datamutex) );
        
        if ( d->reserved_backends.contains(request.subjectUID()) )
            qDebug() << SireError::getPIDString()
                     << "We are already holding a reservation for the UID"
                     << request.subjectUID().toString();

        QHash<QUuid,Frontend> reserved;
        
        foreach (Frontend frontend, available_frontends)
        {
            reserved.insert( frontend.UID(), frontend );
        }
                     
        d->reserved_backends[ request.subjectUID() ] = reserved;
        
        return reserved.keys();
    }
}

/** Collect on the request made with the message 'request', collecting
    the frontend to the backend with UID 'uid', and optionally cancelling
    the reservation of the remaining backends */
Frontend ReservationManager::collectReservation(const ReserveBackend &request,
                                                const QUuid &uid,
                                                bool dispose_of_rest)
{
    ReservationManager *d = reservationManager();
    
    QMutexLocker lkr( &(d->datamutex) );
    
    QHash<QUuid,Frontend> frontends = d->mpifrontends.take(request.subjectUID());
    
    Frontend frontend = frontends.take(uid);
    
    if (not dispose_of_rest)
    {
        d->mpifrontends.insert( request.subjectUID(), frontends );
    }
    else
    {
        frontends.clear();
    }

    return frontend;
}

/** Collect on the request made with the message 'request', collecting all
    the frontends to the backend whose UIDs are in 'uids', and optionally 
    cancelling the reservation of the remaining backends */
QList<Frontend> ReservationManager::collectReservation(const ReserveBackend &request,
                                                       const QList<QUuid> &uids,
                                                       bool dispose_of_rest)
{
    ReservationManager *d = reservationManager();
    
    QMutexLocker lkr( &(d->datamutex) );
    
    QHash<QUuid,Frontend> frontends = d->mpifrontends.take(request.subjectUID());
    
    QList<Frontend> my_frontends;
    
    foreach (QUuid uid, uids)
    {
        Frontend frontend = frontends.take(uid);
    
        if (not frontend.isNull())
            my_frontends.append(frontend);
    }
    
    if (not dispose_of_rest)
    {
        d->mpifrontends.insert( request.subjectUID(), frontends );
    }

    return my_frontends;
}

/** Establish the point-to-point connections between the MPI processes
    as detailed in 'details', in response to the request in 'request' */
void ReservationManager::establishConnections(const QList< tuple<int,QUuid> > &details,
                                              const ReserveBackend &request )
{
    int my_rank = MPICluster::getRank();

    ReservationManager *d = reservationManager();
    
    QMutexLocker lkr( &(d->datamutex) );

    for (int i=0; i<details.count(); ++i)
    {
        tuple<int,QUuid> detail = details[i];
        
        int rank = detail.get<0>();
        QUuid uid = detail.get<1>();
    
        if (uid.isNull())
            continue;
    
        //create a point-to-point communicator between the requested processes
        P2PComm comm = P2PComm( request.sender(), rank );
     
        if (comm.involves(my_rank))
        {
            if (comm.isMaster())
            {
                shared_ptr<FrontendBase> ptr( new MPIFrontend(comm) );
            
                d->mpifrontends[request.subjectUID()].insert(uid, Frontend(ptr));
            }

            if (comm.isSlave())
            {
                //give it the backend (handled via a local frontend)
                comm.setBackend( d->reserved_backends[request.subjectUID()]
                                                            .take(detail.get<1>()) );
                                                            
                d->mpibackends[request.subjectUID()].insert(uid, comm);
            }
        }
    }
    
    //ok, we've processed all of the requests - free up the remaining backends
    d->reserved_backends[ request.subjectUID() ].clear();
    d->reserved_backends.remove( request.subjectUID() );
    
    //also take the opportunity to remove any stale mpibackends
    QMutableHashIterator< QUuid, QHash<QUuid,P2PComm> > it( d->mpibackends );
    
    while (it.hasNext())
    {
        it.next();
        
        QMutableHashIterator<QUuid,P2PComm> it2( it.value() );
        
        while (it2.hasNext())
        {
            it2.next();
            
            if (it2.value().hasFinished())
                it2.remove();
        }
        
        if (it.value().isEmpty())
        {
            it.remove();
        }
    }
}

/** Return the UIDs of all of the backends that have connections
    established in response to the request 'request' */
QList<QUuid> ReservationManager::establishedConnections(const ReserveBackend &request)
{
    ReservationManager *d = reservationManager();
    
    QMutexLocker lkr( &(d->datamutex) );
    
    return d->mpifrontends.value(request.subjectUID()).keys() +
           d->mpibackends.value(request.subjectUID()).keys();
}

#endif // SIRE_USE_MPI
