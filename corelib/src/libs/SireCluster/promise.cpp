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

#include <QThread>
#include <QMutex>
#include <QWaitCondition>

#include "promise.h"
#include "workpacket.h"
#include "node.h"

#include "SireMaths/rangenerator.h"

#include "SireError/errors.h"
#include "SireError/printerror.h"

using namespace SireCluster;

namespace SireCluster
{
namespace detail
{

/** Private implementation of Promise */
class PromisePvt : public QThread
{
public:
    PromisePvt() : QThread()
    {}
    
    ~PromisePvt()
    {
        this->wait();
    }
    
    /** Mutex to protect access to the data
        of this promise */
    QMutex datamutex;
    
    /** The node on which the calculation is running */
    Node node;
    
    /** Wait condition used to wait for a result of the calculation */
    QWaitCondition waiter;
    
    /** The initial state of the WorkPacket before the calculation.
        It will either be in this packet, or compressed into binary */
    WorkPacket initial_packet;
    QByteArray initial_data;
    
    /** The final state of the WorkPacket */
    WorkPacket result_packet;

protected:
    void run()
    {
        SireError::setThreadString("Promise");
        SireMaths::seed_qrand();
        
        QMutexLocker lkr(&datamutex);
        
        //get a local copy of the node
        Node my_node = node;
        lkr.unlock();
        
        //wait until the job has finished
        my_node.wait();
        
        //the node has finished - grab the result and
        //drop our copy of the node
        WorkPacket my_result = my_node.result();
        
        if (my_result.isNull())
        {
            //where did the result go???
            my_result = ErrorPacket( SireError::program_bug( QObject::tr(
                            "There was no result from the running calculation!!!"),
                                CODELOC ) );
        }
        
        //copy the result to the promise
        lkr.relock();
        result_packet = my_result;
        
        //drop the reference to the node
        node = Node();
        
        //wake anyone waiting for a result
        waiter.wakeAll();
    }
};

} // end of namespace detail
} // end of namespace SireCluster;

using namespace SireCluster::detail;

/** Construct a null promise */
Promise::Promise()
{}

/** Internal constructor called by Node that constructs a promise
    that is following the progress of the work in 'initial_packet' 
    as it is being processed by the node 'node' */
Promise::Promise(const Node &node, const WorkPacket &initial_workpacket)
{
    BOOST_ASSERT( not node.isNull() );
    
    d.reset( new PromisePvt() );
    
    d->node = node;
    
    if (initial_workpacket.shouldPack())
    {
        d->initial_data = initial_workpacket.pack();
    }
    else
    {
        d->initial_packet = initial_workpacket;
    }
    
    //now start a background thread that grabs the result
    //as soon as it is available
    d->start();
}

/** Copy constructor */
Promise::Promise(const Promise &other) : d(other.d)
{}

/** Destructor */
Promise::~Promise()
{}

/** Copy assignment operator */
Promise& Promise::operator=(const Promise &other)
{
    d = other.d;
    return *this;
}

/** Comparison operator */
bool Promise::operator==(const Promise &other) const
{
    return d.get() == other.d.get();
}

/** Comparison operator */
bool Promise::operator!=(const Promise &other) const
{
    return d.get() != other.d.get();
}

/** Return whether or not this promise is null */
bool Promise::isNull() const
{
    return d.get() == 0;
}

/** Abort this job */
void Promise::abort()
{
    if (d.get() == 0)
        return;
        
    d->datamutex.lock();
    Node my_node = d->node;
    d->datamutex.unlock();
    
    my_node.abortJob();
}

/** Stop this job */
void Promise::stop()
{
    if (d.get() == 0)
        return;
        
    d->datamutex.lock();
    Node my_node = d->node;
    d->datamutex.unlock();
    
    my_node.stopJob();
}

/** Wait for the job to have completed */
void Promise::wait()
{
    if (d.get() == 0)
        return;
        
    QMutexLocker lkr( &(d->datamutex) );
    
    if (d->result_packet.isNull())
    {
        //we still don't have the result
        while (not d->waiter.wait( &(d->datamutex), 2500 ))
        {
            if (not d->result_packet.isNull())
                //we've got the result!
                return;
        }
    }
}

/** Wait until the job has completed, or until 'timeout' milliseconds
    has passed. This returns whether or not the job has finished */
bool Promise::wait(int timeout)
{
    if (d.get() == 0)
        return true;
        
    QMutexLocker lkr( &(d->datamutex) );
    
    if (d->result_packet.isNull())
    {
        //we still don't have the result
        d->waiter.wait( &(d->datamutex), timeout );
        
        return not d->result_packet.isNull();
    }
    else
        return true;
}

/** Return whether or not the calculation is still in progress */
bool Promise::isRunning()
{
    if (d.get() == 0)
        return false;
        
    QMutexLocker lkr( &(d->datamutex) );
    
    //if we are running, then we still have a handle on the node
    return not d->node.isNull();
}

/** Return whether or not the result is an error.
    This blocks until the result is available */
bool Promise::isError()
{
    if (d.get() == 0)
        return false;

    this->wait();
    
    QMutexLocker lkr( &(d->datamutex) );
    
    return d->result_packet.isError();
}

/** Throw any errors associated with this promise - does
    nothing if there is no error */
void Promise::throwError()
{
    if (d.get() != 0)
    {
        this->wait();
        
        QMutexLocker lkr( &(d->datamutex) );
        
        d->result_packet.throwError();
    }
}

/** Return whether or not the job was stopped.
    This blocks until a result is available */
bool Promise::wasStopped()
{
    if (d.get() == 0)
        return false;

    this->wait();
    
    QMutexLocker lkr( &(d->datamutex) );
    
    return not d->result_packet.hasFinished();
}

/** Return whether or not the job was aborted.
    If it was, then you can rerun the job using
    the initial state of the WorkPacket stored in
    this Promise. This blocks until a result is available */
bool Promise::wasAborted()
{
    if (d.get() == 0)
        return false;

    this->wait();
    
    QMutexLocker lkr( &(d->datamutex) );
    
    return d->result_packet.wasAborted();
}

/** Return the progress of the calculation */
float Promise::progress()
{
    if (d.get() == 0)
        return 1;

    QMutexLocker lkr( &(d->datamutex) );
    
    if (d->result_packet.isNull())
    {
        Node my_node = d->node;
        lkr.unlock();
        
        float current_progress = my_node.progress();
        
        lkr.relock();
        
        if (not d->result_packet.isNull())
            //the result came in while we were getting the progress
            return d->result_packet.progress();
        else
            return current_progress;
    }
    else
        return d->result_packet.progress();
}

/** Return the WorkPacket in the state it was in at the 
    start of the job. You can use this to restart jobs that
    failed in error, or were aborted, or if you just want
    to try to run the job again */
WorkPacket Promise::input()
{
    if (d.get() == 0)
        return WorkPacket();
        
    QMutexLocker lkr( &(d->datamutex) );
    
    if (d->initial_data.isEmpty())
    {
        return d->initial_packet;
    }
    else
    {
        return WorkPacket::unpack( d->initial_data );
    }
}

/** Return an interim result of the calculation */
WorkPacket Promise::interimResult()
{
    if (d.get() == 0)
        return WorkPacket();
        
    QMutexLocker lkr( &(d->datamutex) );
    
    if (not d->result_packet.isNull())
        //we already have the final result!
        return d->result_packet;
        
    else
    {
        Node my_node = d->node;
        
        lkr.unlock();
        
        WorkPacket interim_result = my_node.interimResult();
        
        lkr.relock();
        
        if (not d->result_packet.isNull())
            //we got the final result while waiting for the interim result
            return d->result_packet;
        else
            return interim_result;
    }
}

/** Return the result of the work - this blocks until  
    the work has completed */
WorkPacket Promise::result()
{
    if (d.get() == 0)
        return WorkPacket();
        
    this->wait();
    
    QMutexLocker lkr( &(d->datamutex) );
    
    return d->result_packet;
}
