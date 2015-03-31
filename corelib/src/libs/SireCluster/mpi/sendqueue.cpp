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

#include <mpi.h>  // needs to be first for mpich

#include "sendqueue.h"
#include "mpicluster.h"

#include "SireMaths/rangenerator.h"

#include "SireError/errors.h"
#include "SireError/printerror.h"

#include <QDebug>

using namespace SireCluster::MPI;
using namespace SireCluster::MPI::Messages;

/** Construct, using the passed MPI communicator to send the 
    messages */
SendQueue::SendQueue(MPI_Comm comm)
          : QThread(), boost::noncopyable(),
            send_comm(comm), been_stopped(false)
{}

/** Destructor */
SendQueue::~SendQueue()
{
    datamutex.lock();
    message_queue.clear();
    been_stopped = true;
    datamutex.unlock();
    
    while (not QThread::wait(200))
    {
        waiter.wakeAll();
    }

    MPI_Comm_free( &send_comm );
}

/** Start the event loop in a background thread */
void SendQueue::start()
{
    QMutexLocker lkr(&datamutex);
    been_stopped = false;
    QThread::start();
}

/** Send the message 'message' */
void SendQueue::send(const Message &message)
{
    if (not message.isNull())
    {
        QMutexLocker lkr(&datamutex);
        
        if (not been_stopped)
            message_queue.enqueue(message);

        waiter.wakeAll();
    }
}

/** Wait until the queue has finished */
void SendQueue::wait()
{
    while (not QThread::wait(200))
    {
        waiter.wakeAll();
    }
}

/** Stop the queue - the clears all pending messages */
void SendQueue::stop()
{
    QMutexLocker lkr(&datamutex);
    message_queue.clear();
    been_stopped = true;
    waiter.wakeAll();
    lkr.unlock();
}

/** Return whether or not this queue is running */
bool SendQueue::isRunning()
{
    return QThread::isRunning();
}

/** This is the event loop */
void SendQueue::run()
{
    SireError::setThreadString("SendQueue");
    SireMaths::seed_qrand();
    
    QMutexLocker lkr(&datamutex);
    
    if (message_queue.isEmpty())
    {
        waiter.wait(&datamutex);
    }
    
    while (not message_queue.isEmpty())
    {
        if (been_stopped)
            //we've been stopped!
            break;
    
        Message message = message_queue.dequeue();
        lkr.unlock();

        //qDebug() << MPICluster::getRank() << "sending" << message.toString()
        //         << "to" << message.destination();

        try
        {
            if (message.destination() == MPICluster::getRank())
            {
                //the message is for us! - no need to send it
                MPICluster::received(message);
            }
            else if ( MPICluster::isMaster() )
            {
                //we can directly send the messages ourselves!
                if (message.destination() == -1 and 
                    not message.isA<Broadcast>())
                {
                    //this message is to be broadcast to everyone
                    message = Broadcast(message);
                }
                
                //the master sends the message to just the intended recipients
                QByteArray message_data = message.pack();
                
                int size = message_data.count();
                
                BOOST_ASSERT( size != 0 );
                
                if (message.isA<Broadcast>())
                {
                    QSet<int> recipients = message.asA<Broadcast>().recipients();
                    
                    foreach (int recipient, recipients)
                    {
                        if (recipient != MPICluster::master())
                        {
                            MPI_Request request;
                            MPI_Isend(&size, 1, MPI_INT, recipient, 1, 
                                      send_comm, &request);
                        }
                    }
                    
                    foreach (int recipient, recipients)
                    {
                        if (recipient != MPICluster::master())
                        {
                            MPI_Send(message_data.data(), message_data.count(),
                                     MPI_BYTE, recipient, 1, send_comm);
                        }
                    }

                    if (message.asA<Broadcast>().isRecipient(MPICluster::master()))
                    {
                        //if this is a shutdown, then process it here
                        //(to prevent a deadlock if the receive queue is still working)
                        if (message.isA<Shutdown>())
                        {
                            message.read();
                            been_stopped = true;
                        }
                        else
                            MPICluster::received(message);
                    }
                }
                else
                {
                    int recipient = message.destination();
                    
                    if (recipient != MPICluster::master())
                    {
                        MPI_Request request;
                        MPI_Isend(&size, 1, MPI_INT, recipient, 1, send_comm, &request);
                        MPI_Send(message_data.data(), message_data.count(),
                                 MPI_BYTE, recipient, 1, send_comm);
                    }
                }
            }
            else
            {
                QByteArray message_data;
            
                if (message.destination() == MPICluster::master())
                {
                    //send this message to the master directly
                    message_data = message.pack();
                }
                else
                {
                    //we need to send the message to the master, for 
                    //retransmission to the destination process
                    Message broadcast( Broadcast(message, message.destination()) );
                    message_data = broadcast.pack();
                }
                
                //maybe change to Isend so that we can kill the send if
                //'been_stopped' is true
                MPI_Send(message_data.data(), message_data.count(),
                         MPI_BYTE, MPICluster::master(), 1, send_comm);
            }
        }
        catch(const SireError::exception &e)
        {
            MPICluster::send( Messages::Error(message, e) );
        }
        catch(const std::exception &e)
        {
            MPICluster::send( Messages::Error(message, e) );
        }
        catch(...)
        {
            MPICluster::send( Messages::Error(message, CODELOC) );
        }
        
        lkr.relock();
        
        if (message_queue.isEmpty())
        {
            //wait until there are some more messages
            waiter.wait( &datamutex );
        }
    }
    
    //there are no more messages being broadcast - we need to tell
    //the backend nodes of this fact
    if ( MPICluster::isMaster() )
    {
        int quit = 0;
        //send_comm->Barrier();
        for (int i=0; i<MPICluster::getCount(); ++i)
        {
            if (i != MPICluster::master())
                MPI_Send(&quit, 1, MPI_INT, i, 1, send_comm);
        }
    }
    
    //we're not sending any more messages, so release the resources
    //held by the communicator
    MPI_Comm_free(&send_comm);
}

#endif // SIRE_USE_MPI
