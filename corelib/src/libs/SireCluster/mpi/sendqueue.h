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

#ifndef SIRECLUSTER_MPI_SENDQUEUE_H
#define SIRECLUSTER_MPI_SENDQUEUE_H

#ifdef SIRE_USE_MPI

#include <mpi.h>  // must be at the top, as that is what mpich needs

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

/** This is a simple queue that is used to send all of the MPI messages.
    This ensures that each node only sends one message at a time using
    the global Send communicator (other messages may be sent using
    point-to-point communicators).
    
    @author Christopher Woods
*/
class SendQueue : private QThread, public boost::noncopyable
{
public:
    SendQueue(MPI_Comm send_comm);
    ~SendQueue();
    
    void start();
    
    void send(const Message &message);
    
    void stop();
    
    void wait();
    
    bool isRunning();
    
protected:
    void run();
    
private:
    /** Mutex to protect access to the queue of messages to send */
    QMutex datamutex;
    
    /** Wait condition used to sleep until there is a message to send */
    QWaitCondition waiter;
    
    /** The communicator to use to send messages */
    MPI_Comm send_comm;
    
    /** The list of messages to send */
    QQueue<Message> message_queue;
    
    /** Whether or not the queue has been stopped */
    bool been_stopped;
};

} // end of namespace MPI

} // end of namespace SireCluster

SIRE_END_HEADER

#endif // SIRE_USE_MPI
#endif

