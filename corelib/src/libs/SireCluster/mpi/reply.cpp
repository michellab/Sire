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

#include <QTime>
#include <QMutex>
#include <QWaitCondition>
#include <QHash>
#include <QSet>

#include "reply.h"
#include "messages.h"
#include "mpicluster.h"

#include "SireError/errors.h"

#include "tostring.h"

#include <QDebug>

using namespace SireCluster::MPI;

using boost::shared_ptr;
using boost::weak_ptr;

namespace SireCluster
{
namespace MPI
{
namespace detail
{

/** Private implementation of Reply */
class ReplyPvt
{
public:
    ReplyPvt()
    {}
    
    ~ReplyPvt()
    {}
    
    /** Mutex to protect access to this response */
    QMutex datamutex;
    
    /** Wait conditions so that we can wait for the response
        from each process */
    QHash< int, shared_ptr<QWaitCondition> > waiters;
    
    /** The responses received so far... */
    QHash< int, ReplyValue > responses;
};

} // end of namespace detail
} // end of namespace MPI
} // end of namespace SireCluster

using namespace SireCluster::MPI::detail;

///////////
/////////// Implementation of ReplyValue
///////////

/** Null constructor */
ReplyValue::ReplyValue() : is_error(false)
{}

/** Copy constructor */
ReplyValue::ReplyValue(const ReplyValue &other)
           : result_data(other.result_data), is_error(other.is_error)
{}
     
/** Destructor */
ReplyValue::~ReplyValue()
{}

/** Copy assignment operator */
ReplyValue& ReplyValue::operator=(const ReplyValue &other)
{
    if (this != &other)
    {
        result_data = other.result_data;
        is_error = other.is_error;
    }
    
    return *this;
}

/** Comparison operator */
bool ReplyValue::operator==(const ReplyValue &other) const
{
    return is_error == other.is_error and
           result_data == other.result_data;
}

/** Comparison operator */
bool ReplyValue::operator!=(const ReplyValue &other) const
{
    return not this->operator==(other);
}

/** Create a ReplyValue that represents the result contained
    in 'result_data' */
ReplyValue ReplyValue::result(const QByteArray &result_data)
{
    ReplyValue val;
    val.is_error = false;
    val.result_data = result_data;
    
    return val;
}

/** Create a ReplyValue that represents the error contained
    in 'error_data' */
ReplyValue ReplyValue::error(const QByteArray &error_data)
{
    ReplyValue error;
    error.is_error = true;
    error.result_data = error_data;
    
    return error;
}

/** Return whether or not this result is an error */
bool ReplyValue::isError() const
{
    return is_error;
}

///////////
/////////// Implementation of Reply
///////////

/** Null constructor */
Reply::Reply()
{}

/** Copy assignment operator */
Reply& Reply::operator=(const Reply &other)
{
    d = other.d;
    return *this;
}

/** Construct the reply for the message 'message' */
Reply::Reply(const Message &message)
{
    this->operator=( MPICluster::getReply(message) );
}

/** Construct from the passed pointer */
Reply::Reply(const ReplyPtr &ptr)
{
    this->operator=( ptr.lock() );
}

/** Internal constructor */
Reply::Reply(const shared_ptr<ReplyPvt> &ptr) : d(ptr)
{}

/** Copy constructor */
Reply::Reply(const Reply &other) : d(other.d)
{}

/** Destructor */
Reply::~Reply()
{}

/** Comparison operator */
bool Reply::operator==(const Reply &other) const
{
    return d.get() == other.d.get();
}

/** Comparison operator */
bool Reply::operator!=(const Reply &other) const
{
    return d.get() != other.d.get();
}

/** Return whether or not this reply is null */
bool Reply::isNull() const
{
    return d.get() == 0;
}

/** Internal function called by MPICluster used to create 
    a new reply for the message 'message' */
Reply Reply::create(const Message &message)
{
    Reply reply;
    
    reply.d.reset( new ReplyPvt() );
    
    QMutexLocker lkr( &(reply.d->datamutex) );
    
    QSet<int> recipients = message.recipients();
    
    foreach (int recipient, recipients)
    {
        reply.d->waiters.insert(recipient, 
                                shared_ptr<QWaitCondition>(new QWaitCondition()) );
    }
    
    return reply;
}

/** Inform the reply that an error occured while trying to get the
    reply from the process with rank 'rank' */
void Reply::setErrorFrom(int rank, const QByteArray &error_data)
{
    if (this->isNull())
        return;

    QMutexLocker lkr( &(d->datamutex) );
    
    if (d->waiters.contains(rank))
    {
        d->responses.insert(rank, ReplyValue::error(error_data));
        d->waiters.value(rank)->wakeAll();
    }
}

/** Inform the reply that a reply has been received from the 
     process with rank 'rank' */
void Reply::setResultFrom(int rank, const QByteArray &result_data)
{
    if (this->isNull())
        return;

    QMutexLocker lkr( &(d->datamutex) );
    
    if (d->waiters.contains(rank))
    {
        d->responses.insert(rank, ReplyValue::result(result_data));
        d->waiters.value(rank)->wakeAll();
    }
}

/** Inform this reply that the process with rank 'rank' has
    gone down (and so we should stop listening for replies) */
void Reply::setProcessDown(int rank)
{
    if (this->isNull())
        return;

    QMutexLocker lkr( &(d->datamutex) );
    
    if (d->waiters.contains(rank))
    {
        if (not d->responses.contains(rank))
        {
            SireError::unavailable_resource e( QObject::tr(
                    "The node with rank %1 has gone down while waiting "
                    "for a reply.").arg(rank), CODELOC );
        
            d->responses.insert(rank, ReplyValue::error( e.pack() ));
        }
        
        //wake everyone waiting for this process
        d->waiters.value(rank)->wakeAll();
    }
}

/** Wake all replies as we are about to be shut down */
void Reply::shutdown()
{
    if (this->isNull())
        return;
        
    QMutexLocker lkr( &(d->datamutex) );
    
    QList<int> ranks = d->waiters.keys();
    
    if (d->responses.count() != ranks.count())
    {
        //we need to fill in some missing responses
        SireError::unavailable_resource e( QObject::tr(
                "There are no more replies available as we are being shutdown."),
                    CODELOC );
                    
        QByteArray error_data = e.pack();
        
        foreach (int rank, ranks)
        {
            if (not d->responses.contains(rank))
            {
                d->responses.insert(rank, ReplyValue::error(error_data));
            }
        }
    }
    
    //make sure everyone is awake
    foreach (int rank, ranks)
    {
        d->waiters.value(rank)->wakeAll();
    }
}

/** Return whether or not the rank 'rank' is valid for this reply */
bool Reply::isValidRank(int rank)
{
    if (this->isNull())
        return false;

    QMutexLocker lkr( &(d->datamutex) );
    return d->waiters.contains(rank);
}

/** Wait for a reply from the process with rank 'rank' */
void Reply::waitFrom(int rank)
{
    if (this->isNull())
        return;

    QMutexLocker lkr( &(d->datamutex) );
    
    if (d->waiters.contains(rank))
    {
        if (not d->responses.contains(rank))
        {
            d->waiters.value(rank)->wait( &(d->datamutex) );
        }
    }
}

/** Wait for a reply from the process with rank 'rank', or until
    'timeout' milliseconds have passed. This returns whether or not
    there is a reply. */
bool Reply::waitFrom(int rank, int timeout)
{
    QTime t;
    t.start();

    if (timeout < 0)
    {
        this->waitFrom(rank);
        return true;
    }

    if (this->isNull())
        return true;

    QMutexLocker lkr( &(d->datamutex) );
    
    if (d->waiters.contains(rank))
    {
        if (not d->responses.contains(rank))
        {
            timeout -= t.elapsed();
            if (timeout < 0)
                return false;
        
            return d->waiters.value(rank)->wait( &(d->datamutex), timeout );
        }
        else
            return true;
    }
    else
        return true;
}

/** Wait for all of the necessary replies to have been returned */
void Reply::wait()
{
    if (this->isNull())
        return;
    
    QMutexLocker lkr( &(d->datamutex) );
    
    QList<int> ranks = d->waiters.keys();
    
    foreach (int rank, ranks)
    {
        if (not d->responses.contains(rank))
        {
            d->waiters.value(rank)->wait( &(d->datamutex) );
        }
    }
}

/** Wait for all of the necessary replies to have been returned,
    or until 'timeout' milliseconds have passed. This returns whether or not
    all of the replies have been received. */
bool Reply::wait(int timeout)
{
    QTime t;
    t.start();

    if (timeout < 0)
    {
        this->wait();
        return true;
    }

    if (this->isNull())
        return true;

    QMutexLocker lkr( &(d->datamutex) );
    
    QList<int> ranks = d->waiters.keys();
    
    foreach (int rank, ranks)
    {
        if (not d->responses.contains(rank))
        {
            timeout -= t.elapsed();
            t.start();
            
            if (timeout < 0)
                return false;
                
            if (not d->waiters.value(rank)->wait( &(d->datamutex), timeout ))
            {
                return false;
            }
        }
    }
    
    return true;
}

/** Return the value of the reply from the process with 
    rank 'rank'. This blocks until the reply is available 
    
    \throw SireError::invalid_arg
*/
ReplyValue Reply::from(int rank)
{
    if (this->isNull())
        throw SireError::invalid_arg( QObject::tr(
            "There is no reply from the process with rank %1 in a null Reply!")
                .arg(rank), CODELOC );

    QMutexLocker lkr( &(d->datamutex) );
    
    if (not d->waiters.contains(rank))
        throw SireError::invalid_arg( QObject::tr(
            "There is no reply from the process with rank %1 in this reply. "
            "We are only waiting for replies from the processes with ranks %2.")
                .arg(rank).arg(Sire::toString(d->waiters.keys())), CODELOC );
                
    while (not d->responses.contains(rank))
    {
        //we need to wait for a response
        d->waiters.value(rank)->wait( &(d->datamutex) );
    }
    
    return d->responses.value(rank);
}

/** Return the list of all replies indexed by the rank of the 
    process that supplied the reply - this blocks until they are
    all available */
QHash<int,ReplyValue> Reply::replies()
{
    if (this->isNull())
        return QHash<int,ReplyValue>();
    
    this->wait();
    
    QMutexLocker lkr( &(d->datamutex) );

    return d->responses;
}

///////////
/////////// Implementation of ReplyPtr
///////////

/** Construct a null pointer */
ReplyPtr::ReplyPtr()
{}

/** Construct to point to the reply 'reply' */
ReplyPtr::ReplyPtr(const Reply &reply) : d(reply.d)
{}

/** Copy constructor */
ReplyPtr::ReplyPtr(const ReplyPtr &other) : d(other.d)
{}

/** Destructor */
ReplyPtr::~ReplyPtr()
{}

/** Copy assignment operator */
ReplyPtr& ReplyPtr::operator=(const ReplyPtr &other)
{
    d = other.d;
    return *this;
}

/** Return whether or not this is a null pointer */
bool ReplyPtr::isNull() const
{
    return d.expired();
}

/** Return the reply - this will return a null reply
    if this is a null pointer */
Reply ReplyPtr::operator*() const
{
    return Reply( d.lock() );
}

/** Return the reply - this will return a null reply
    if this is a null pointer */
Reply ReplyPtr::lock() const
{
    return Reply( d.lock() );
}

#endif // SIRE_USE_MPI
