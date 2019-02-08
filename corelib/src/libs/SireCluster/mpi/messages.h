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

#ifndef SIRECLUSTER_MPI_MESSAGES_H
#define SIRECLUSTER_MPI_MESSAGES_H

#ifdef SIRE_USE_MPI

#include "SireBase/sharedpolypointer.hpp"

#include <QSharedData>
#include <QUuid>
#include <QSet>
#include <QList>

#include <boost/tuple/tuple.hpp>

SIRE_BEGIN_HEADER

namespace SireCluster
{

class Frontend;

namespace MPI
{

class MessageBase;
class Message;

namespace Messages
{
class RegisterBackend;
class Shutdown;
class Broadcast;
class GetUIDs;
class ReserveBackend;
class RequestAvailability;
class Reservation;

class Error;
class Result;
}

}
}

QDataStream& operator<<(QDataStream&, const SireCluster::MPI::Message&);
QDataStream& operator>>(QDataStream&, SireCluster::MPI::Message&);

QDataStream& operator<<(QDataStream&, const SireCluster::MPI::MessageBase&);
QDataStream& operator>>(QDataStream&, SireCluster::MPI::MessageBase&);

QDataStream& operator<<(QDataStream&, 
                        const SireCluster::MPI::Messages::RegisterBackend&);
QDataStream& operator>>(QDataStream&, 
                        SireCluster::MPI::Messages::RegisterBackend&);

QDataStream& operator<<(QDataStream&, 
                        const SireCluster::MPI::Messages::Shutdown&);
QDataStream& operator>>(QDataStream&, 
                        SireCluster::MPI::Messages::Shutdown&);

QDataStream& operator<<(QDataStream&, 
                       const SireCluster::MPI::Messages::Broadcast&);
QDataStream& operator>>(QDataStream&, 
                        SireCluster::MPI::Messages::Broadcast&);

QDataStream& operator<<(QDataStream&, 
                       const SireCluster::MPI::Messages::ReserveBackend&);
QDataStream& operator>>(QDataStream&, 
                        SireCluster::MPI::Messages::ReserveBackend&);

QDataStream& operator<<(QDataStream&, 
                       const SireCluster::MPI::Messages::RequestAvailability&);
QDataStream& operator>>(QDataStream&, 
                        SireCluster::MPI::Messages::RequestAvailability&);

QDataStream& operator<<(QDataStream&, 
                       const SireCluster::MPI::Messages::Reservation&);
QDataStream& operator>>(QDataStream&, 
                        SireCluster::MPI::Messages::Reservation&);

QDataStream& operator<<(QDataStream&, 
                       const SireCluster::MPI::Messages::GetUIDs&);
QDataStream& operator>>(QDataStream&, 
                        SireCluster::MPI::Messages::GetUIDs&);

QDataStream& operator<<(QDataStream&, 
                       const SireCluster::MPI::Messages::Error&);
QDataStream& operator>>(QDataStream&, 
                        SireCluster::MPI::Messages::Error&);

QDataStream& operator<<(QDataStream&, 
                       const SireCluster::MPI::Messages::Result&);
QDataStream& operator>>(QDataStream&, 
                        SireCluster::MPI::Messages::Result&);

namespace SireError
{
class exception;
}

namespace SireCluster
{
namespace MPI
{

/** This is the virtual base class of all MPI messages sent between processes
    in the MPI cluster 
    
    @author Christopher Woods
*/
class MessageBase : public QSharedData
{

friend SIRECLUSTER_EXPORT QDataStream& ::operator<<(QDataStream&, const MessageBase&);
friend SIRECLUSTER_EXPORT QDataStream& ::operator>>(QDataStream&, MessageBase&);

public:
    MessageBase(const MessageBase &other);
    
    virtual ~MessageBase();
    
    static const char* typeName()
    {
        return "SireCluster::MPI::MessageBase";
    }
    
    virtual const char* what() const=0;
    
    virtual MessageBase* clone() const=0;
    
    virtual void read()=0;
    
    virtual QString toString() const;
    
    virtual bool hasReply() const;
    
    virtual Message reply() const;
    
    const QUuid& UID() const;
    const QUuid& subjectUID() const;
    
    int sender() const;
    int destination() const;
    
    virtual bool isRecipient(int rank) const;
    
    virtual QSet<int> recipients() const;
    
    template<class T>
    bool isA() const;
    
    template<class T>
    T asA() const;
    
protected:
    MessageBase(int destination);
    MessageBase(int destination, const QUuid &subject_uid);
    MessageBase();

    MessageBase& operator=(const MessageBase &other);

private:
    /** The unique ID of this message (this allows messages and
        responses to be tallied */
    QUuid uid;

    /** The ID of the subject that this message relates to
         - this allows several related messages to be exchanged
           and tracked */
    QUuid subject_uid;

    /** The MPI rank of the process that sent this message */
    qint32 sent_by;
    
    /** The MPI rank of the process that should receive this message */
    qint32 dest;
};

/** The polymorphic holder of all of the message classes

    @author Christopher Woods
*/
class Message
{

friend SIRECLUSTER_EXPORT QDataStream& ::operator<<(QDataStream&, const Message&);
friend SIRECLUSTER_EXPORT QDataStream& ::operator>>(QDataStream&, Message&);

public:
    Message();
    Message(const MessageBase &message);
    
    Message(const Message &other);
    
    ~Message();
    
    Message& operator=(const Message &other);
    
    static const char* typeName();
    
    const char* what() const
    {
        if (this->isNull())
            return Message::typeName();
        else
            return d->what();
    }
    
    QString toString() const;
    
    bool isNull() const;

    QByteArray pack() const;
    static Message unpack(const QByteArray &data);

    void read();
    
    bool hasReply() const;
    
    Message reply() const;
    
    const QUuid& UID() const;
    const QUuid& subjectUID() const;
    
    int sender() const;
    int destination() const;
    
    bool isRecipient(int rank) const;
    
    QSet<int> recipients() const;
    
    template<class T>
    bool isA() const
    {
        if (this->isNull())
            return false;
        else
            return base().isA<T>();
    }
    
    template<class T>
    T asA() const
    {
        return base().asA<T>();
    }

private:
    const MessageBase& base() const;

    /** Polymorphic, implicitly shared pointer to the message */
    SireBase::SharedPolyPointer<MessageBase> d;
};

namespace Messages
{

/** This is the message sent by the master node when it re-broadcasts
    messages sent by the slave nodes (as the slave nodes can't talk
    to each other directly - not until they open up a point-to-point
    communicator)
    
    @author Christopher Woods
*/
class Broadcast : public MessageBase
{

friend SIRECLUSTER_EXPORT QDataStream& ::operator<<(QDataStream&, const Broadcast&);
friend SIRECLUSTER_EXPORT QDataStream& ::operator>>(QDataStream&, Broadcast&);

public:
    Broadcast();
    Broadcast(const Message &message);
    Broadcast(const Message &message, int rank);
    Broadcast(const Message &message, const QSet<qint32> &ranks);
    
    Broadcast(const Broadcast &other);
    
    ~Broadcast();
    
    Broadcast& operator=(const Broadcast &other);
    
    static const char* typeName();
    
    const char* what() const
    {
        return Broadcast::typeName();
    }
    
    Broadcast* clone() const;
    
    QString toString() const;
    
    void read();
    
    bool hasReply() const;
    
    Message reply() const;
    
    bool isRecipient(int rank) const;
    
    QSet<int> recipients() const;
    
    const QString& messageType() const;
    
    template<class T>
    bool messageIsA() const
    {
        return message_type == T::typeName();
    }
    
    template<class T>
    T messageAsA() const
    {
        return Message::unpack(message_data).asA<T>();
    }
    
private:
    /** The actual destination of this message - this is empty
        if it is meant to go to everybody */
    QSet<qint32> destinations;
    
    /** The type of the message being broadcast */
    QString message_type;
    
    /** The data containing the message being broadcast */
    QByteArray message_data;
    
    /** Any reply to this message */
    Message replymsg;
};

/** This is the message sent by a slave node to tell the master
    node that a new Backend has been created on the slave which
    the master should register in the global registry 
    
    @author Christopher Woods
*/
class RegisterBackend : public MessageBase
{

friend SIRECLUSTER_EXPORT QDataStream& ::operator<<(QDataStream&, const RegisterBackend&);
friend SIRECLUSTER_EXPORT QDataStream& ::operator>>(QDataStream&, RegisterBackend&);

public:
    RegisterBackend();
    RegisterBackend(const QUuid &node_uid);

    RegisterBackend(const RegisterBackend &other);
    
    ~RegisterBackend();
    
    RegisterBackend& operator=(const RegisterBackend &other);
    
    static const char* typeName();
    
    const char* what() const
    {
        return RegisterBackend::typeName();
    }
    
    RegisterBackend* clone() const;
    
    QString toString() const;
    
    void read();

private:
    /** The UID of the node to be registered */
    QUuid node_uid;
    
    /** The rank of the MPI process that contains this node */
    qint32 node_rank;
};

/** This message is sent by a node to the master to request the complete
    list of UIDs of all of the backends. This message expects a response,
    and will block waiting for that response
    
    @author Christopher Woods
*/
class GetUIDs : public MessageBase
{

friend SIRECLUSTER_EXPORT QDataStream& ::operator<<(QDataStream&, const GetUIDs&);
friend SIRECLUSTER_EXPORT QDataStream& ::operator>>(QDataStream&, GetUIDs&);

public:
    GetUIDs();
    
    GetUIDs(const GetUIDs &other);
    
    ~GetUIDs();
    
    GetUIDs& operator=(const GetUIDs &other);
    
    static const char* typeName();
    
    const char* what() const
    {
        return GetUIDs::typeName();
    }
    
    GetUIDs* clone() const;
    
    void read();
    
    bool hasReply() const;
    Message reply() const;
};

/** This message is sent to the master to request the reservation of
    backends from remote nodes 
    
    @author Christopher Woods
*/
class ReserveBackend : public MessageBase
{

friend SIRECLUSTER_EXPORT QDataStream& ::operator<<(QDataStream&, const ReserveBackend&);
friend SIRECLUSTER_EXPORT QDataStream& ::operator>>(QDataStream&, ReserveBackend&);

public:
    ReserveBackend();
    ReserveBackend(int nbackends);
    ReserveBackend(const QUuid &backend_uid);

    ReserveBackend(const ReserveBackend &other);
    
    ~ReserveBackend();
    
    ReserveBackend& operator=(const ReserveBackend &other);
    
    static const char* typeName();
    
    const char* what() const
    {
        return ReserveBackend::typeName();
    }
    
    ReserveBackend* clone() const;
    
    const QUuid& requestedUID() const;
    
    int nBackends() const;
    
    QString toString() const;
    
    void read();

private:
    /** The UID of the backend to be reserved - null if we
        don't care which backend we get */
    QUuid backend_uid;
    
    /** The number of backends to reserve */
    qint32 nbackends;
};

/** This message is broadcast by the master to all processes to
    request that they return their availability to meet the 
    
    
    @author Christopher Woods
*/
class RequestAvailability : public MessageBase
{

friend SIRECLUSTER_EXPORT QDataStream& ::operator<<(QDataStream&, const RequestAvailability&);
friend SIRECLUSTER_EXPORT QDataStream& ::operator>>(QDataStream&, RequestAvailability&);

public:
    RequestAvailability();
    RequestAvailability(const ReserveBackend &reservation_request);

    RequestAvailability(const RequestAvailability &other);
    
    ~RequestAvailability();
    
    RequestAvailability& operator=(const RequestAvailability &other);
    
    static const char* typeName();
    
    const char* what() const
    {
        return RequestAvailability::typeName();
    }
    
    RequestAvailability* clone() const;
    
    QString toString() const;
    
    void read();

    bool hasReply() const;
    
    Message reply() const;

private:
    /** The reservation request */
    ReserveBackend request;
    
    /** The available backends on this process */
    QList<QUuid> available_backends;
};

/** This message is sent back from the master, containing the
    reservation details saying which backends are available
    for the connection (and, indeed, must be connected to)
            
    @author Christopher Woods
*/
class Reservation : public MessageBase
{

friend SIRECLUSTER_EXPORT QDataStream& ::operator<<(QDataStream&, const Reservation&);
friend SIRECLUSTER_EXPORT QDataStream& ::operator>>(QDataStream&, Reservation&);

public:
    Reservation();
    
    Reservation(const QList< boost::tuple<int,QUuid> > &reserved_backends,
                const ReserveBackend &request);

    Reservation(const Reservation &other);
    
    ~Reservation();
    
    Reservation& operator=(const Reservation &other);
    
    static const char* typeName();
    
    const char* what() const
    {
        return Reservation::typeName();
    }
    
    Reservation* clone() const;
    
    QString toString() const;
    
    void read();

    bool hasReply() const;
    
    Message reply() const;

private:
    /** The initial request that generated this reservation */
    ReserveBackend request;

    /** The reservation details for each acquired backend */
    QList< boost::tuple<int,QUuid> > details;
};

/** This message is sent to return a result
    
    @author Christopher Woods
*/
class Result : public MessageBase
{

friend SIRECLUSTER_EXPORT QDataStream& ::operator<<(QDataStream&, const Result&);
friend SIRECLUSTER_EXPORT QDataStream& ::operator>>(QDataStream&, Result&);

public:
    Result();
    
    template<class T>
    Result(const Message &message, const T &value)
          : MessageBase(message.sender(), message.subjectUID())
    {
        QDataStream ds(&result_data, QIODevice::WriteOnly);
        ds << value;
    }
    
    Result(const Result &other);
    
    ~Result();
    
    Result& operator=(const Result &other);
    
    static const char* typeName();
    
    const char* what() const
    {
        return Result::typeName();
    }
    
    Result* clone() const;
    
    void read();

private:
    /** A binary representation of the result */
    QByteArray result_data;
};

/** This message is sent to return an error!
    
    @author Christopher Woods
*/
class Error : public MessageBase
{

friend SIRECLUSTER_EXPORT QDataStream& ::operator<<(QDataStream&, const Error&);
friend SIRECLUSTER_EXPORT QDataStream& ::operator>>(QDataStream&, Error&);

public:
    Error();
    Error(const QString &code_loc);
    Error(const SireError::exception &e);
    Error(const std::exception &e);
    
    Error(const Message &message);
    Error(const Message &message, const QString &code_loc);
    Error(const Message &message, const SireError::exception &e);
    Error(const Message &message, const std::exception &e);
    
    Error(int sender);
    Error(int sender, const QString &code_loc);
    Error(int sender, const SireError::exception &e);
    Error(int sender, const std::exception &e);
    
    Error(const QByteArray &message_data, const QByteArray &error_data);
    
    Error(const Error &other);
    
    ~Error();
    
    Error& operator=(const Error &other);
    
    static const char* typeName();
    
    const char* what() const
    {
        return Error::typeName();
    }
    
    Error* clone() const;
    
    void read();

private:
    void setError(const SireError::exception &e);
    void setError(const std::exception &e);
    void setUnknownError(const QString &code_loc);

    /** A binary representation of the message that caused the error */
    QByteArray message_data;
    
    /** A binary representation of the error */
    QByteArray error_data;
};

/** This is the message sent to a node to tell it to shutdown. This
    is used by the master node to tell all of the slaves
    to shutdown when the program is exiting 
    
    @author Christopher Woods
*/
class Shutdown : public MessageBase
{

friend SIRECLUSTER_EXPORT QDataStream& ::operator<<(QDataStream&, const Shutdown&);
friend SIRECLUSTER_EXPORT QDataStream& ::operator>>(QDataStream&, Shutdown&);

public:
    Shutdown();
    Shutdown(const Shutdown &other);
    
    ~Shutdown();
    
    Shutdown& operator=(const Shutdown &other);
    
    static const char* typeName();
    
    const char* what() const
    {
        return Shutdown::typeName();
    }
    
    Shutdown* clone() const;

    void read();
};

}  // end of namespace Messages

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

template<>
inline bool MessageBase::isA<Messages::Broadcast>() const
{
    return dynamic_cast<const Messages::Broadcast*>(this) != 0;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool MessageBase::isA() const
{
    if (dynamic_cast<const T*>(this) == 0)
    {
        if (this->isA<Messages::Broadcast>())
        {
            return this->asA<Messages::Broadcast>().messageIsA<T>();
        }
        else
            return false;
    }
    else
        return true;
}

template<>
inline Messages::Broadcast MessageBase::asA<Messages::Broadcast>() const
{
    return *(dynamic_cast<const Messages::Broadcast*>(this));
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
T MessageBase::asA() const
{
    const T *ptr = dynamic_cast<const T*>(this);
    
    if (ptr == 0)
    {
        if (this->isA<Messages::Broadcast>())
        {
            return this->asA<Messages::Broadcast>().messageAsA<T>();
        }
        else
            return *ptr;
    }
    else
        return *ptr;
}

#endif // SIRE_SKIP_INLINE_FUNCTIONS

}  // end of namespace MPI
}  // end of namespace SireCluster

Q_DECLARE_METATYPE( SireCluster::MPI::Messages::Broadcast )
Q_DECLARE_METATYPE( SireCluster::MPI::Messages::Result )
Q_DECLARE_METATYPE( SireCluster::MPI::Messages::Error )

Q_DECLARE_METATYPE( SireCluster::MPI::Messages::RegisterBackend )
Q_DECLARE_METATYPE( SireCluster::MPI::Messages::ReserveBackend )
Q_DECLARE_METATYPE( SireCluster::MPI::Messages::RequestAvailability )
Q_DECLARE_METATYPE( SireCluster::MPI::Messages::Reservation )
Q_DECLARE_METATYPE( SireCluster::MPI::Messages::Shutdown )
Q_DECLARE_METATYPE( SireCluster::MPI::Messages::GetUIDs )

Q_DECLARE_METATYPE( SireCluster::MPI::Message )

SIRE_END_HEADER

#endif // SIRE_USE_MPI

#endif
