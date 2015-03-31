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

#ifndef SIRECLUSTER_MPI_REPLY_H
#define SIRECLUSTER_MPI_REPLY_H

#ifdef SIRE_USE_MPI

#include "SireError/exception.h"

#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>

#include <QByteArray>
#include <QDataStream>

SIRE_BEGIN_HEADER

namespace SireCluster
{
namespace MPI
{

class MPICluster;
class Message;
class ReplyPtr;

namespace detail
{
class ReplyPvt;
}

/** This class provides the holder for the value of a reply */
class ReplyValue
{
public:
    ReplyValue();
    ReplyValue(const ReplyValue &other);
    
    ~ReplyValue();
    
    ReplyValue& operator=(const ReplyValue &other);
    
    bool operator==(const ReplyValue &other) const;
    bool operator!=(const ReplyValue &other) const;

    static ReplyValue result(const QByteArray &result_data);
    static ReplyValue error(const QByteArray &error_data);
    
    bool isError() const;
    
    template<class T>
    bool isA() const
    {
        if (this->is_error)
            return false;
            
        try
        {
            QDataStream ds(this->result_data);
            T value;
            ds >> value;
            
            return true;
        }
        catch(...)
        {
            return false;
        }
    }
    
    template<class T>
    T asA() const
    {
        if (this->is_error)
            SireError::exception::unpackAndThrow(result_data);
            
        QDataStream ds(this->result_data);
        T value;
        ds >> value;
        
        return value;
    }
    
private:
    /** The binary representation of the value */
    QByteArray result_data;
    
    /** Is this reply actually an error? */
    bool is_error;
};

/** This class hold the reply (or replies) to a message. This class
    provides a space in which the replies will be received, and can
    be used as a block to wait for a reply
    
    Replies are arranged according to the rank of the MPI process
    that sent the reply
    
    @author Christopher Woods
*/
class Reply
{

friend class MPICluster;
friend class ReplyPtr;

public:
    Reply();
    Reply(const Message &message);
    Reply(const ReplyPtr &ptr);

    Reply(const Reply &other);
    
    ~Reply();
    
    Reply& operator=(const Reply &other);
    
    bool operator==(const Reply &other) const;
    bool operator!=(const Reply &other) const;
    
    bool isNull() const;
    
    bool isValidRank(int rank);
    
    void waitFrom(int rank);
    bool waitFrom(int rank, int timeout);
    
    void wait();
    bool wait(int timeout);

    ReplyValue from(int rank);

    QHash<int,ReplyValue> replies();

protected:
    static Reply create(const Message &message);

    void setErrorFrom(int rank, const QByteArray &error_data);
    void setResultFrom(int rank, const QByteArray &result_data);

    void setProcessDown(int rank);
    void shutdown();

private:
    Reply(const boost::shared_ptr<detail::ReplyPvt> &ptr);

    /** Private implementation of this object */
    boost::shared_ptr<detail::ReplyPvt> d;
};

/** This is a weak pointer to a Reply */
class ReplyPtr
{
public:
    ReplyPtr();
    ReplyPtr(const Reply &reply);
    
    ReplyPtr(const ReplyPtr &other);
    
    ~ReplyPtr();
    
    ReplyPtr& operator=(const ReplyPtr &other);
    
    bool isNull() const;
    
    Reply operator*() const;

    Reply lock() const;

private:
    /** The weak pointer itself */
    boost::weak_ptr<detail::ReplyPvt> d;
};

} // end of namespace MPI
} // end of namespace SireCluster

SIRE_END_HEADER

#endif // SIRE_USE_MPI
#endif
