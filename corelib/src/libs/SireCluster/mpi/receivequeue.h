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

#ifndef SIRECLUSTER_MPI_RECEIVEQUEUE_H
#define SIRECLUSTER_MPI_RECEIVEQUEUE_H

#ifdef SIRE_USE_MPI

#include <mpi.h> // mpich requires that this goes first

#include <QQueue>
#include <QThread>
#include <QMutex>
#include <QWaitCondition>

#include <boost/noncopyable.hpp>

#include "messages.h"

SIRE_BEGIN_HEADER

namespace SireCluster
{
namespace MPI
{

/** This is a simple queue that listens for messages from MPI
    processes and adds them onto a queue in which they are read
    (and processed). This ensures that only one message is 
    processed at a time from the global MPI recv communicator, and
    that only one message is received at a time (although the 
    next message may be being received while the current
    message is being processed - also, this is just the global
    communicator - point-to-point communicators can send and 
    receive whenever they want)
    
    @author Christopher Woods
*/
class ReceiveQueue : private QThread, public boost::noncopyable
{

protected:
class SecondThread;

friend class ReceiveQueue::SecondThread;

public:
    ReceiveQueue(MPI_Comm recv_comm);
    ~ReceiveQueue();
    
    void start();
    
    void received(const Message &message);
    
    void stop();

    void wait();
    
    bool isRunning();

protected:
    void run();
    void run2();
    
    /** We need to have two background threads running to overlay
        computation with communication - one thread reads messages
        and adds them to the queue, while another thread reads and
        processes messages from the queue */
    class SecondThread : public QThread
    {
    public:
        SecondThread(ReceiveQueue *parent) : QThread(), p(parent)
        {}
        
        ~SecondThread()
        {}
        
    protected:
        void run()
        {
            p->run2();
        }
    
    private:
        ReceiveQueue *p;
    };

private:
    Message unpackMessage(const QByteArray &message_data, int sender) const;

    SecondThread *secondthread;

    /** Mutex to protect access to the queue of received messages */
    QMutex datamutex;
    
    /** Wait condition used to sleep until there is a received message to read */
    QWaitCondition waiter;
    
    /** The communicator to use to receive messages */
    MPI_Comm recv_comm;
    
    /** The list of received messages */
    QQueue<Message> message_queue;

    /** Whether or not to keep listening for messages */
    bool keep_listening;
    
    /** Whether or not this has been stopped */
    bool been_stopped;
};

} // end of namespace MPI
} // end of namespace SireCluster

SIRE_END_HEADER

#endif // SIRE_USE_MPI
#endif
