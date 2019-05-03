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

#include <mpi.h>   //mpich requires that mpi.h is included first

#include <QHash>
#include <QMutex>
#include <QUuid>
#include <QTextStream>

#include "SireCluster/cluster.h"
#include "SireCluster/frontend.h"
#include "SireCluster/backend.h"

#include "mpicluster.h"
#include "messages.h"
#include "sendqueue.h"
#include "receivequeue.h"
#include "reply.h"
#include "reservationmanager.h"
#include "p2pcomm.h"

#include "SireError/errors.h"
#include "SireError/printerror.h"

#include "sire_config.h"

#include <QDebug>

using namespace SireCluster;
using namespace SireCluster::MPI;

using boost::shared_ptr;

#ifndef HAVE_LSEEK64
    //////
    ////// add an lseek64 function stub to fill a function
    ////// that is missing - mpich needs lseek64 to be
    ////// defined, even if it is not available! Otherwise
    ////// dlopen errors as the symbol can't be found
    //////
    extern "C"
    {
        int lseek64(int fd, int offset, int whence)
        {
            throw SireError::program_bug( QObject::tr(
                "MPI implementation is calling lseek64 which is not supported "
                "on OS X (Leopard - 32bit)"), CODELOC );
            
            return 0;
        }
    }
#endif // HAVE_LSEEK64

/** Private implementation of MPICluster */
class MPIClusterPvt
{
public:
    MPIClusterPvt();
    
    ~MPIClusterPvt();
    
    /** Mutex to protect access to the data of this cluster */
    QMutex datamutex;
    
    /** Mutex used to protect access to the reply registry */
    QMutex replymutex;
    
    /** The UIDs of all backends, together with the rank of the 
        process that contains that backend */
    QHash<QUuid,int> backend_registry;

    /** All of the active replys on this process - this provides
        holders that will be filled with the replies to messages.
        This is indexed by subject UID */
    QHash<QUuid,ReplyPtr> reply_registry;

    /** The global (private) MPI communicator */
    MPI_Comm global_comm;

    /** The send message event loop */
    SendQueue *send_queue;
    
    /** The receive message event loops */
    ReceiveQueue *receive_queue;
    
    /** Whether or not we are being shutdown (or have shutdown) */
    bool already_shutting_down;
};

Q_GLOBAL_STATIC( QMutex, mpiGlobalMutex );

static MPIClusterPvt *global_cluster = 0;

static MPIClusterPvt* globalCluster()
{
    if (global_cluster == 0)
    {
        QMutexLocker lkr( mpiGlobalMutex() );
        
        if (global_cluster == 0)
            global_cluster = new MPIClusterPvt();
    }

    return global_cluster;
}

static void ensureMPIStarted()
{
    int initialized;
    MPI_Initialized(&initialized);

    if (not initialized)
    {
        int argc = 0;
        char **argv = 0;
        
        //Absolutely must use multi-threaded MPI
        int level;
        int mpierr = MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &level);
        
        if (mpierr != MPI_SUCCESS)
            throw SireError::unavailable_resource( QObject::tr(
                    "Unable to start multi-threaded MPI! MPI Error code %1.")
                        .arg(mpierr), CODELOC );
                        
        if (level != MPI_THREAD_MULTIPLE)
            throw SireError::unavailable_resource( QObject::tr(
                    "Multi-threaded MPI not available. MPI level %1.")
                        .arg(level), CODELOC );
    }
}

/** Construct the global cluster */
MPIClusterPvt::MPIClusterPvt() : send_queue(0), receive_queue(0),
                                 already_shutting_down(true)
{
    //create the private send and receive communicators
    ::ensureMPIStarted();

    //now make sure that we are using a version of MPI that
    //has multi-thread support
    int thread_level;
    MPI_Query_thread(&thread_level);
    
    if (thread_level != MPI_THREAD_MULTIPLE)
    {
        //we can't run without thread support, as we use multiple
        //communicators running in multiple threads
        if (MPICluster::isMaster())
        {
            QTextStream ts(stderr);
        
            ts << QObject::tr("Sire needs to use an MPI that supports threads.\n"
                    "This means that MPI must have been started using "
                    "MPI::Init_thread(MPI_THREAD_MULTIPLE), and that the MPI\n"
                    "library supports multithreaded MPI.\n\n"
                    "Check if your MPI library has been compiled with thread support.\n");
        }
        
        //stop MPI and exit
        MPI_Barrier( MPI_COMM_WORLD );
        MPI_Finalize();
        
        //kill the program
        std::exit(-1);
    }
        
    //get the global MPI communicator
    MPI_Group mpigroup;
    MPI_Comm_group(MPI_COMM_WORLD, &mpigroup);

    //create a new communicator that is just used by SireCluster
    MPI_Comm_create(MPI_COMM_WORLD, mpigroup, &global_comm);

    MPI_Comm send_comm, recv_comm;
        
    if ( MPICluster::isMaster() )
    {
        //create the send, then receive communicators
        MPI_Comm_create(global_comm, mpigroup, &send_comm);
        MPI_Comm_create(global_comm, mpigroup, &recv_comm);
    }
    else
    {
        //must be the other way around (as all other nodes
        //listen to the master)
        MPI_Comm_create(global_comm, mpigroup, &recv_comm);
        MPI_Comm_create(global_comm, mpigroup, &send_comm);
    }

    MPI_Group_free(&mpigroup);
    
    //create and start the reservation manager
    ReservationManager::start();
        
    //create and start the event loops
    send_queue = new SendQueue(send_comm);
    receive_queue = new ReceiveQueue(recv_comm);
        
    send_queue->start();
    receive_queue->start();

    //wait for everyone to get here
    MPI_Barrier(global_comm);
    
    already_shutting_down = false;
}

/** Destructor */
MPIClusterPvt::~MPIClusterPvt()
{
    //stop the reservation manager
    ReservationManager::shutdown();

    MPI_Comm_free(&global_comm);

    if (send_queue)
    {
        send_queue->stop();
        send_queue->wait();
        delete send_queue;
        send_queue = 0;
    }
    
    if (receive_queue)
    {
        receive_queue->stop();
        receive_queue->wait();
        delete receive_queue;
        receive_queue = 0;
    }
}

/** Start the MPI backend */
void MPICluster::start()
{
    ::ensureMPIStarted();
    globalCluster();
}

/** Synchronise the MPI processes - this can be used
    as a Barrier to ensure that all processes have reached
    the same point */
void MPICluster::sync()
{
    ::ensureMPIStarted();

    int is_finalized;
    MPI_Finalized(&is_finalized);

    if (not is_finalized)
        MPI_Barrier(globalCluster()->global_comm);
}

/** Create a new P2P communicator that allow for direct and
    private communicator between the processes with ranks
    'master_rank' and 'slave_rank' */
P2PComm MPICluster::createP2P(int master_rank, int slave_rank)
{
    int my_rank = MPICluster::getRank();
        
    bool is_master = (my_rank == master_rank);
    bool is_slave = (my_rank == slave_rank);
        
    if (master_rank == slave_rank)
    {   
        //this is an intra-process communicator - no 
        //need to do anything, unless it is us!
        if (is_master)
        {
            return P2PComm::createLocal();
        }
        else
            return P2PComm();
    }
    else
    {
        //we need to create a new communicator for this process
        // - this requires a collective operation
        QMutexLocker lkr( &(globalCluster()->datamutex) );
        
        int is_finalized;
        MPI_Finalized(&is_finalized);
        if (is_finalized)
            return P2PComm();
        
        int rank = 0;
        
        if (is_slave)
            rank = 1;
        
        MPI_Comm private_comm;
        MPI_Comm_split(globalCluster()->global_comm, 
                       (is_master or is_slave), is_slave, &private_comm);
        
        return P2PComm::create(private_comm, master_rank, slave_rank);
    }
}

/** Return the rank of this process in the MPI cluster */
int MPICluster::getRank()
{
    ::ensureMPIStarted();

    int finalized, rank;
    MPI_Finalized(&finalized);

    if (not finalized)
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        return rank;
    }
    else
        //MPI has finalized
        return -1;
}

/** Return the number of processes in the MPI cluster */
int MPICluster::getCount()
{
    ::ensureMPIStarted();
    
    int finalized, size;
    MPI_Finalized(&finalized);
    
    if (not finalized)
    {
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        return size;
    }
    else
        //MPI has shut down
        return 0;
}

/** Return the rank of the master process */
int MPICluster::master()
{
    return 0;
}

/** Return whether or not this is the master process */
bool MPICluster::isMaster()
{
    return MPICluster::getRank() == MPICluster::master();
}

/** Send the message 'message' (the destination is available
    in message.destination()) */
void MPICluster::send(const Message &message)
{
    if (message.isNull())
        return;
        
    globalCluster()->send_queue->send(message);
}

/** Receive the message 'message' */
void MPICluster::received(const Message &message)
{
    if (message.isNull())
        return;
        
    globalCluster()->receive_queue->received(message);
}

/** Return the reply object for the message 'message' - this
    will create a reply object for this message if one doesn't
    already exist */ 
Reply MPICluster::getReply(const Message &message)
{
    if (message.isNull())
        return Reply();

    QMutexLocker lkr( &(globalCluster()->replymutex) );
    
    Reply reply = globalCluster()->reply_registry.value(message.subjectUID());
    
    if (reply.isNull())
    {
        //we need to create the reply
        reply = Reply::create(message);
        
        globalCluster()->reply_registry.insert(message.subjectUID(), reply);
    }
    
    //take this opportunity to clear out any dead replies
    QMutableHashIterator<QUuid,ReplyPtr> it(globalCluster()->reply_registry);
    
    while (it.hasNext())
    {
        it.next();
        
        if (it.value().isNull())
        {
            it.remove();
        }
    }
    
    return reply;
}

/** Post the result contained in 'result_data' for the message with
    subject 'subject_uid', which has come from the process with 
    rank 'sender' */
void MPICluster::postResult(const QUuid &subject_uid, int sender,
                            const QByteArray &result_data)
{
    QMutexLocker lkr( &(globalCluster()->replymutex) );
    
    Reply reply = globalCluster()->reply_registry.value(subject_uid).lock();
    
    if (reply.isNull())
    {
        qDebug() << "There is no reply on process" << MPICluster::getRank()
                 << "awaiting the result on the subject"
                 << subject_uid.toString() << "sent by the process" << sender;
                 
        return;
    }
    
    if (not reply.isValidRank(sender))
    {
        qDebug() << "There is no space in the reply on process" << MPICluster::getRank()
                 << "for the result on the subject"
                 << subject_uid.toString() << "sent by the process" << sender;
                 
        return;
    }
    
    reply.setResultFrom(sender, result_data);
}
 
static void fatalError(const QByteArray &error_data,
                       const QByteArray &message_data, int sender)
{
    if (MPICluster::isMaster())
    {
        QTextStream ts(stdout);
        
        ts << QObject::tr(
                "\n*********************************************************\n"
                "There was a fatal error on the MPI process with rank %1.\n\n")
                    .arg(sender);
                    
        try
        {
            Message message = Message::unpack(message_data);
            
            ts << QObject::tr("The offending message is %1.\n"
                       "It has UID %2, and was sent by process %3 to "
                       "destination %4.\n\n")
                            .arg(message.toString(), message.UID().toString())
                            .arg(message.sender()).arg(message.destination());
        }
        catch(...)
        {
            ts << QObject::tr(
                    "No information about the offending message is available.\n\n");
        }
        
        ts << QObject::tr("_____ Here is the actual error _____\n");
        
        try
        {
            shared_ptr<SireError::exception> 
                        e = SireError::exception::unpack(error_data);
                        
            SireError::printError( *e );
        }
        catch(const SireError::exception &e2)
        {
            SireError::printError(e2);
        }
        catch(...)
        {
            SireError::printError( SireError::unknown_exception( QObject::tr(
                "Something went wrong when reading the error!"), CODELOC ) );
        }
        
        //now shutdown the cluster
        Cluster::shutdown();
    }
    else
    {
        MPICluster::send( Messages::Error(message_data, error_data) );
    }
}
 
/** Post the error contained in 'error_data' for the message with
    subject 'subject_uid' (contained in 'message_data'), which has 
    come from the process with rank 'sender' */
void MPICluster::postError(const QUuid &subject_uid, int sender,
                           const QByteArray &message_data,
                           const QByteArray &error_data)
{
    QMutexLocker lkr( &(globalCluster()->replymutex) );
    
    Reply reply = globalCluster()->reply_registry.value(subject_uid).lock();
    
    if (reply.isNull())
    {
        //there is no space for a reply - there is no way to report a problem,
        //so the safest thing is to send this to the master and shutdown the 
        //cluster
        ::fatalError(error_data, message_data, sender);
    }
    
    if (not reply.isValidRank(sender))
    {
        //there is no space for a reply from the sender - there is no way
        //to report a problem, so again it is best if we shut down
        ::fatalError(error_data, message_data, sender);
    }
    
    reply.setErrorFrom(sender, error_data);
}

/** Call this function on the master MPI process to register
    that the backend with UID 'uid' is on the MPI process
    with rank 'rank' */
void MPICluster::registerBackend(int rank, const QUuid &uid)
{
    if ( not MPICluster::isMaster() )
    {
        //why is this message here - it should have been
        //sent to the master!
        throw SireError::program_bug( QObject::tr(
            "A request to register the node with UID %1 on process %2 "
            "has ended up on process %3, while the master is process %4.")
                .arg(uid.toString())
                .arg(rank)
                .arg( MPICluster::getRank() )
                .arg( MPICluster::master() ),
                    CODELOC );
    }

    QMutexLocker lkr( &(globalCluster()->datamutex) );
    
    if (not globalCluster()->backend_registry.contains(uid))
    {
        globalCluster()->backend_registry.insert(uid, rank);
    }
}    

/** Register the local Backend 'backend' with the MPI cluster so that
    it can be connected to by any MPI-connected node */
void MPICluster::registerBackend(const Backend &backend)
{
    if ( MPICluster::isMaster() )
    {
        MPICluster::registerBackend( MPICluster::master(), backend.UID() );
    }
    else
    {
        MPICluster::send( Messages::RegisterBackend(backend.UID()) );
    }
}

/** Return the first available front end. This returns a null
    frontend if there are no backends available
*/
Frontend MPICluster::getFrontend()
{
    //the Cluster should already have looked for local nodes...

    //ask the master to reserve any backend for us
    //(the master will only reserve remote backends)
    Messages::ReserveBackend message(1);
    
    //create space to hold the reply (which contains the reservation)
    Reply reply(message);
    
    MPICluster::send(message);
    
    //wait for the reply
    reply.wait();
    
    //now get the UIDs of the nodes that have established a direct
    //connection to us
    QList<QUuid> uids = reply.from( MPICluster::master() )
                             .asA< QList<QUuid> >();

    if (uids.isEmpty())
    {
        return Frontend();
    }
    else
    {
        return ReservationManager::collectReservation(message, uids.first());
    }
}

/** Return the first 'n' available frontends. This returns an empty
    list if there are no backends available
*/
QList<Frontend> MPICluster::getFrontends(int n)
{
    if (n <= 0)
        return QList<Frontend>();

    //the Cluster should already have looked for local nodes...

    //ask the master to reserve any backend for us
    //(the master will only reserve remote backends)
    Messages::ReserveBackend message(n);
    
    //create space to hold the reply (which contains the reservation)
    Reply reply(message);
    
    MPICluster::send(message);
    
    //wait for the reply
    reply.wait();
    
    //now get the UIDs of the nodes that have established a direct
    //connection to us
    QList<QUuid> uids = reply.from( MPICluster::master() )
                             .asA< QList<QUuid> >();

    while (uids.count() > n)
    {
        uids.removeLast();
    }

    return ReservationManager::collectReservation(message, uids);
}

/** Return the frontend for backend with UID 'uid'.
    This returns a null frontend if this backend isn't
    currently available, and it raises an error if there
    is no backend associated with this UID
    
    \throw SireError::unavailable_resource
*/
Frontend MPICluster::getFrontend(const QUuid &uid)
{
    if (uid.isNull())
        throw SireError::unavailable_resource( QObject::tr(
            "There is no front end for the null backend!"), CODELOC );

    //this interface is only to get *REMOTE* backends
    BOOST_ASSERT( not Cluster::isLocal(uid) );

    //ask the master to arrange the connection to this backend
    Messages::ReserveBackend message(uid);
    
    //create space to hold the reply to this message
    Reply reply(message);
    
    MPICluster::send(message);
        
    //wait for the reply
    reply.wait();
    
    //now get the UIDs of the nodes that have established a direct
    //connection to us
    QList<QUuid> uids = reply.from( MPICluster::master() )
                             .asA< QList<QUuid> >();

    if (uids.isEmpty())
        return Frontend();
    else
        return ReservationManager::collectReservation(message, uids.first());
}

/** Return the list of all of the UIDs of all of the backends
    that are available via MPI */
QList<QUuid> MPICluster::UIDs()
{
    if (MPICluster::isMaster())
    {
        QMutexLocker lkr( &(globalCluster()->datamutex) );
        return globalCluster()->backend_registry.keys();
    }
    else
    {
        //ask the master for the list of UIDs
        Messages::GetUIDs message;
    
        //create space to hold the reply to this message
        Reply reply(message);
    
        MPICluster::send( message );

        //wait for all of the responses
        reply.wait();
    
        //return the result
        return reply.from( MPICluster::master() ).asA< QList<QUuid> >();
    }
}

/** Return whether or not this MPI cluster contains a backend
    with UID 'uid' */
bool MPICluster::hasBackend(const QUuid &uid)
{
    if (MPICluster::isMaster())
    {
        QMutexLocker lkr( &(globalCluster()->datamutex) );
        return globalCluster()->backend_registry.contains(uid);
    }
    else
    {
        return MPICluster::UIDs().contains(uid);
    }
}

/** Return whether or not the MPICluster is running */
bool MPICluster::isRunning()
{
    int initialized;
    MPI_Initialized(&initialized);

    if (not initialized)
        return false;
        
    else
    {
        return globalCluster()->send_queue->isRunning() or
               globalCluster()->receive_queue->isRunning();
    }
}

Q_GLOBAL_STATIC( QMutex, shutdownMutex );

/** This just shuts down this node */
void MPICluster::informedShutdown()
{
    //// Check we aren't already in the process of being shutdown
    {
        QMutexLocker lkr( shutdownMutex() );
        
        if (globalCluster()->already_shutting_down)
            return;
            
        globalCluster()->already_shutting_down;
    }

    if (not MPICluster::isRunning())
        return;

    globalCluster()->send_queue->stop();
    globalCluster()->receive_queue->stop();
    
    globalCluster()->send_queue->wait();
    globalCluster()->receive_queue->wait();

    ReservationManager::shutdown();

    ////////// make sure that no-one is still waiting for a reply
    {
        QMutexLocker lkr( &(globalCluster()->replymutex) );
        
        for (QHash<QUuid,ReplyPtr>::const_iterator 
                    it = globalCluster()->reply_registry.constBegin();
             it != globalCluster()->reply_registry.constEnd();
             ++it)
        {
            it.value().lock().shutdown();
        }
    }

    //////////
    {
        QMutexLocker lkr( &(globalCluster()->datamutex) );
        globalCluster()->backend_registry.clear();

        //finally, get rid of the global communicator
        MPI_Barrier( globalCluster()->global_comm );
        MPI_Comm_free( &(globalCluster()->global_comm) );
    }
    
    //shut down the cluster - this will call MPICluster::shutdown,
    //but this won't recurse, as MPICluster::isRunning() will be false
    Cluster::shutdown();
}

/** Shutdown the MPI cluster */
void MPICluster::shutdown()
{
    if (not MPICluster::isRunning())
        return;

    //this will shutdown the entire cluster!
    MPICluster::send( Messages::Shutdown() );
}

#endif // SIRE_USE_MPI
