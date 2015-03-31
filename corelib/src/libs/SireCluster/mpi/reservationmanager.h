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

#ifndef SIRECLUSTER_MPI_RESERVATIONMANAGER_H
#define SIRECLUSTER_MPI_RESERVATIONMANAGER_H

#ifdef SIRE_USE_MPI

#include <QUuid>
#include <QThread>
#include <QMutex>
#include <QWaitCondition>
#include <QQueue>
#include <QHash>

#include <boost/tuple/tuple.hpp>

#include "messages.h"

namespace SireCluster
{
namespace MPI
{

class Reply;
class P2PComm;

/** This class is used to handle the reservation of backends. This is
    used as the MPI code used to grab nodes is a three step process;
    
    (1) The process wanting a backend makes a request
    (2) This request is broadcast to all processes, and backends that
        match this request are then reserved
    (3) The process wanting a backend receives the list of matching
        backends, chooses which one(s) it wants, establishes connection(s),
        and then frees the rest.
        
    Because the list of what is available is separate from the actual
    connection, it is necessary to reserve the nodes that have been
    said to be available
    
    @author Christopher Woods
*/
class ReservationManager : private QThread
{
public:
    ReservationManager();
    ~ReservationManager();

    static void start();
    static void shutdown();

    static QList<Frontend> collectReservation( const Messages::ReserveBackend &request,
                                               const QList<QUuid> &uids,
                                               bool dispose_of_rest=true );

    static Frontend collectReservation( const Messages::ReserveBackend &request,
                                        const QUuid &uid,
                                        bool dispose_of_rest=true );

    static void reserveBackends(const Messages::ReserveBackend &request);

    static void awaitResponse( const Messages::ReserveBackend &request, 
                               const Reply &reply );
                               
    static QList<QUuid> findAvailable( const Messages::ReserveBackend &request );

    static void establishConnections(const QList< boost::tuple<int,QUuid> > &details,
                                     const Messages::ReserveBackend &request );
                                     
    static QList<QUuid> establishedConnections(const Messages::ReserveBackend &request);

protected:
    void run();

private:
    void processRequest( Messages::ReserveBackend request );

    void freeReservedBackends();

    /** Mutex to protect access to the data of this manager */
    QMutex datamutex;

    /** Mutex used to ensure that only one process makes a request
        at a time */
    QMutex reservation_mutex;

    /** Waiter used to wake up this manager when a new request
        is made */
    QWaitCondition waiter;

    /** The queue of requests to process */
    QQueue<Messages::ReserveBackend> request_queue;

    /** The collection of backends that have been reserved, indexed
        by the subject UID of the message used to request the 
        reservation (and then the UID of the backend) */
    QHash< QUuid, QHash<QUuid,Frontend> > reserved_backends;

    /** The set of MPI backends running on this process, indexed by
        the subject UID of the request */
    QHash< QUuid, QHash<QUuid,P2PComm> > mpibackends;
    
    /** The set of MPI frontends running on this process, indexed
        by the subject UID of the request */
    QHash< QUuid, QHash<QUuid,Frontend> > mpifrontends;

    /** Whether or not to keep going */
    bool keep_going;
};

} // end of namespace MPI
} // end of namespace SireCluster

#endif // SIRE_USE_MPI
#endif

