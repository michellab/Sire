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

#include <mpi.h>  // must be first to satisy mpich

#include <QThread>
#include <QMutex>

#include "p2pcomm.h"
#include "mpicluster.h"

#include "SireCluster/frontend.h"
#include "SireCluster/workpacket.h"

#include "SireMaths/rangenerator.h"

#include "SireError/errors.h"
#include "SireError/printerror.h"

#include <QDebug>

using namespace SireCluster;
using namespace SireCluster::MPI;

namespace SireCluster
{
namespace MPI
{
namespace detail
{

/** Private implementation of P2PComm */
class P2PCommPvt : public QThread
{
public:
    P2PCommPvt() : keep_running(false)
    {
        //create a local P2P
        master_rank = MPICluster::getRank();
        slave_rank = master_rank;
        is_master = true;
        is_slave = true;
    }
    
    P2PCommPvt(MPI_Comm comm, int master, int slave)
            : private_comm(comm), master_rank(master), slave_rank(slave),
              keep_running(false)
    {
        int my_rank = MPICluster::getRank();

        if (private_comm == MPI_COMM_NULL)
        {
            throw SireError::program_bug( QObject::tr(
                "You cannot create a P2PComm with a null communicator!"),   
                    CODELOC );
        }
        
        //private_comm.Barrier();
        
        is_master = (my_rank == master_rank);
        is_slave = (my_rank == slave_rank);
    }
    
    ~P2PCommPvt()
    {
        datamutex.lock();
        keep_running = false;
        datamutex.unlock();
        
        //wait until the background thread has stopped
        this->wait();
    
        //free the communicator
        int finalized;
        MPI_Finalized(&finalized);
        
        if (not finalized)
        {
            if (private_comm != MPI_COMM_NULL)
            {
                if (is_master)
                {
                    //stop the backend
                    int envelope[2];
                    envelope[0] = P2PComm::EXIT;
                    envelope[1] = 0;
                    
                    MPI_Send(envelope, 2, MPI_INT, P2PComm::SLAVE, 1, private_comm);
                }
            
                //private_comm.Barrier();
                MPI_Comm_free(&private_comm);
            }
        }
    }
    
    void waitForResponse(int rank, int tag)
    {
        MPI_Status status;
        
        int received_message;
        MPI_Iprobe(rank, tag, private_comm, &received_message, &status);
        
        while (not received_message)
        {
            QThread::msleep(5);
            MPI_Iprobe(rank, tag, private_comm, &received_message, &status);
        }
    }
    
    /** Mutex to protect access to this communicator */
    QMutex datamutex;

    /** A local frontend to the backend that sits behind
        this communicator (this is only non-null on the slave) */
    Frontend local_backend;
    
    /** The private MPI communicator for this backend */
    MPI_Comm private_comm;
    
    /** The rank of the master process in this communicator */
    int master_rank;
    
    /** The rank of the slave process in this communicator */
    int slave_rank;
    
    /** Whether this is the master process */
    bool is_master;
    
    /** Whether this is the slave process */
    bool is_slave;
    
    /** Whether or not to keep looping */
    bool keep_running;

protected:
    void run();
    
    void sendIntegerResponse(int response);
    void sendFloatResponse(float response);
    
    template<class T>
    void sendResponse(const T &response);
};

} // end of namespace detail
} // end of namespace MPI
} // end of namespace SireCluster

using namespace SireCluster::MPI::detail;

/** Function used by the background thread of the slave to send
    answers back to the master */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void P2PCommPvt::sendResponse(const T &response)
{
    QByteArray data;
    
    QDataStream ds( &data, QIODevice::WriteOnly );
    
    ds << response;
    
    //send this to the master
    int size = data.count();
    
    MPI_Send(&size, 1, MPI_INT, P2PComm::MASTER, 1, private_comm);
    
    if (size > 0)
    {
        MPI_Send(data.data(), size, MPI_BYTE, P2PComm::MASTER, 1, private_comm);
    }
}

/** This just sends a quick integer response back to the master */
void P2PCommPvt::sendIntegerResponse(int response)
{
    MPI_Send(&response, 1, MPI_INT, P2PComm::MASTER, 1, private_comm);
}

/** This just sends a quick float response back to the master */
void P2PCommPvt::sendFloatResponse(float response)
{
    MPI_Send(&response, 1, MPI_FLOAT, P2PComm::MASTER, 1, private_comm);
}

template<class T>
static T decode(const QByteArray &data)
{
    T obj;
    QDataStream ds(data);
    ds >> obj;
    
    return obj;
}

/** Background thread containing the event loop for the slave P2PComm */
void P2PCommPvt::run()
{
    SireError::setThreadString( "P2PComm_Slave" );
    SireMaths::seed_qrand();

    QMutexLocker lkr( &datamutex );
    
    if (not is_slave)
        return;
    
    int envelope[2];
    int message;
    QByteArray message_data;
    
    MPI_Status status;
    
    while (keep_running)
    {
        lkr.unlock();

        //wait for a message from the master
        this->waitForResponse(P2PComm::MASTER, 1);
        
        MPI_Recv(envelope, 2, MPI_INT, P2PComm::MASTER, 1, private_comm, &status);

        if (envelope[1] > 0)
        {
            message_data.resize(envelope[1]+1);
            MPI_Recv(message_data.data(), envelope[1], MPI_BYTE,
                     P2PComm::MASTER, 1, private_comm, &status);
        }
        
        lkr.relock();
        
        message = envelope[0];
        
        if (message == P2PComm::EXIT)
                break;
    
        else if (message == P2PComm::GETUID)
        {
            this->sendResponse( local_backend.UID() );
        }
        else if (message == P2PComm::START)
        {
            local_backend.startJob( decode<WorkPacket>(message_data) );
            this->sendIntegerResponse(0);
        }
        else if (message == P2PComm::STOP)
        {
            local_backend.stopJob();
            this->sendIntegerResponse(0);
        }
        else if (message == P2PComm::ABORT)
        {
            local_backend.abortJob();
            this->sendIntegerResponse(0);
        }
        else if (message == P2PComm::IS_RUNNING)
        {
            if (local_backend.wait(100))
            {
                this->sendIntegerResponse(0);
            }
            else
            {
                this->sendIntegerResponse(1);
            }
        }
        else if (message == P2PComm::PROGRESS)
        {
            this->sendFloatResponse( local_backend.progress() );
        }
        else if (message == P2PComm::INTERIM)
        {
            this->sendResponse( local_backend.interimResult() );
        }
        else if (message == P2PComm::RESULT)
        {
            this->sendResponse( local_backend.result() );
        }
        else
        {
            qDebug() << SireError::getPIDString() 
                     << "Unrecognised message!" << message;
        }
    }

    //private_comm.Barrier();
    MPI_Comm_free(&private_comm);
    private_comm = MPI_COMM_NULL;
    
    //return the backend to the pool
    local_backend = Frontend();
}

/** Null constructor */
P2PComm::P2PComm()
{}

/** Construct the point-to-point communicator between the 
    processes with ranks 'master_rank' and 'slave_rank' in 
    the global MPICluster::communicator(). Note that this
    is a collective operation, so *must* be called on all
    processes in the MPI cluster at the same time. */
P2PComm::P2PComm(int master_rank, int slave_rank)
{
    this->operator=( MPICluster::createP2P(master_rank, slave_rank) );
}

/** Copy constructor */
P2PComm::P2PComm(const P2PComm &other) : d(other.d)
{}

/** Destructor */
P2PComm::~P2PComm()
{}

/** Copy assignment operator */
P2PComm& P2PComm::operator=(const P2PComm &other)
{
    d = other.d;
    return *this;
}

/** Return whether or not this is a null P2P communicator */
bool P2PComm::isNull() const
{
    return d.get() == 0;
}

/** Return whether or not this is the master side of the 
    communicator (the process holding the frontend)
    
    This can be both, if this is an intra-process communicator
*/
bool P2PComm::isMaster()
{
    if (this->isNull())
        return false;
        
    else
    {
        QMutexLocker lkr( &(d->datamutex) );
        return d->is_master;
    }
}

/** Return whether or not this is the slave side of the 
    communicator (the process holding the backend)
    
    This can be both, if this is an intra-process communicator
*/
bool P2PComm::isSlave()
{
    if (this->isNull())
        return false;
        
    else
    {
        QMutexLocker lkr( &(d->datamutex) );
        return d->is_slave;
    }
}

/** Return whether or not this communicator is local
    (is within the same MPI process) */
bool P2PComm::isLocal()
{
    if (this->isNull())
        return true;
        
    else
    {
        QMutexLocker lkr( &(d->datamutex) );
        return d->is_master and d->is_slave;
    }
}

/** Internal function called by MPICluster used to create a 
    point-to-point communicator that works for just this process */
P2PComm P2PComm::createLocal()
{
    P2PComm p2p;
    
    p2p.d.reset( new P2PCommPvt() );
    
    return p2p;
}

/** Internal function called by MPICluster used to create a 
    point-to-point communicator that works using the provided
    communicator, with this as the master process if 'is_master'
    is true */
P2PComm P2PComm::create(MPI_Comm private_comm,
                        int master_rank, int slave_rank)
{
    P2PComm p2p;
    
    int my_rank = MPICluster::getRank();
    
    if (my_rank == master_rank or my_rank == slave_rank)
    {
        p2p.d.reset( new P2PCommPvt(private_comm, master_rank, slave_rank) );
    }
    else
    {
        //we are not involved
        if (private_comm != MPI_COMM_NULL)
        {
            //private_comm.Barrier();
            MPI_Comm_free(&private_comm);
            private_comm = MPI_COMM_NULL;
        }
    }
    
    return p2p;
}

/** Return whether or not this communicator involves the process
    with MPI rank 'rank' */
bool P2PComm::involves(int rank)
{
    if (this->isNull())
        return false;
        
    else
    {
        int my_rank = MPICluster::getRank();
    
        QMutexLocker lkr( &(d->datamutex) );
        
        return (my_rank == d->master_rank) or 
               (my_rank == d->slave_rank);
    }
}

/** Send the message 'message' to the slave */
void P2PComm::sendMessage(int message)
{
    this->_pvt_sendMessage(message, QByteArray());
}

/** Internal function used to send a message to the slave that comes with
    some additional data */
void P2PComm::_pvt_sendMessage(int message, const QByteArray &data)
{
    int finalized;
    MPI_Finalized(&finalized);

    if (finalized)
        return;

    QMutexLocker lkr( &(d->datamutex) );
    
    if (not d->is_master)
        return;
    
    if (d->private_comm == MPI_COMM_NULL)
        return;
        
    int envelope[2];
    envelope[0] = message;
    envelope[1] = data.count();
    
    MPI_Send(envelope, 2, MPI_INT, P2PComm::SLAVE, 1, d->private_comm);

    if (not data.isEmpty())
    {
        MPI_Send(const_cast<char*>(data.data()), data.count(), 
                 MPI_BYTE, P2PComm::SLAVE, 1, d->private_comm);
    }
}

/** Await an integer response from the slave */
int P2PComm::awaitIntegerResponse(bool urgent)
{
    int finalized;
    MPI_Finalized(&finalized);

    if (finalized)
        return -1;

    QMutexLocker lkr( &(d->datamutex) );
    
    if (not d->is_master)
        return -1;
    
    if (d->private_comm == MPI_COMM_NULL)
        return -1;

    int message;
    
    if (not urgent)
        d->waitForResponse(P2PComm::SLAVE, 1);

    MPI_Status status;
    MPI_Recv(&message, 1, MPI_INT, P2PComm::SLAVE, 1, d->private_comm, &status);

    return message;
}

/** Await a float response from the slave */
float P2PComm::awaitFloatResponse(bool urgent)
{
    int finalized;
    MPI_Finalized(&finalized);

    if (finalized)
        return 0;

    QMutexLocker lkr( &(d->datamutex) );
    
    if (not d->is_master)
        return 0;
    
    if (d->private_comm == MPI_COMM_NULL)
        return 0;

    float message;
    
    if (not urgent)
        d->waitForResponse(P2PComm::SLAVE, 1);

    MPI_Status status;
    MPI_Recv(&message, 1, MPI_FLOAT, P2PComm::SLAVE, 1, d->private_comm, &status);

    return message;
}

/** Internal function used to wait for a response from the slave */
QByteArray P2PComm::_pvt_awaitResponse(bool urgent)
{
    int finalized;
    MPI_Finalized(&finalized);

    if (finalized)
        return QByteArray();

    QMutexLocker lkr( &(d->datamutex) );
    
    if (not d->is_master)
        return QByteArray();
    
    if (d->private_comm == MPI_COMM_NULL)
        return QByteArray();

    int size;
    
    if (not urgent)
        d->waitForResponse(P2PComm::SLAVE, 1);

    MPI_Status status;
    MPI_Recv(&size, 1, MPI_INT, P2PComm::SLAVE, 1, d->private_comm, &status);
    
    if (size <= 0)
        return QByteArray();
        
    QByteArray data( size+1, ' ' );

    MPI_Recv(data.data(), size, MPI_BYTE, P2PComm::SLAVE, 1, d->private_comm, &status);
    
    return data;
}

/** Set the backend that will be controlled by the MPI frontend
    (this is held using a local Frontend) */
void P2PComm::setBackend(const Frontend &backend)
{
    if (not this->isSlave())
        throw SireError::program_bug( QObject::tr(
            "Only the slave P2PComm holds the backends..."), CODELOC );

    QMutexLocker lkr( &(d->datamutex) );
    
    if (not d->local_backend.isNull())
        throw SireError::program_bug( QObject::tr(
            "You cannot give a second backend to a P2PComm!"), CODELOC );
            
    d->local_backend = backend;
    
    //start the background thread that listens for communication 
    //from the master
    d->keep_running = true;
    d->start();
}

/** Return whether or not the slave backend has finished */
bool P2PComm::hasFinished()
{
    QMutexLocker lkr( &(d->datamutex) );

    return d->local_backend.isNull();
}

#endif // SIRE_USE_MPI
