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

#include <QThread>
#include <QMutex>
#include <QWaitCondition>

#include "messages.h"
#include "mpicluster.h"
#include "reply.h"
#include "reservationmanager.h"

#include "SireMaths/rangenerator.h"

#include "SireCluster/cluster.h"

#include "SireError/exception.h"
#include "SireError/printerror.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>

using namespace SireCluster;
using namespace SireCluster::MPI;
using namespace SireCluster::MPI::Messages;

using namespace SireBase;
using namespace SireStream;

using boost::tuple;

/////////
///////// Implementation of MessageBase
/////////

static const RegisterMetaType<MessageBase> r_msgbase( MAGIC_ONLY, 
                                                      MessageBase::typeName() );
                                                      
/** Serialise to a binary datastream */
QDataStream& operator<<(QDataStream &ds, const MessageBase &msgbase)
{
    writeHeader(ds, r_msgbase, 1);
    
    ds << msgbase.uid << msgbase.subject_uid
       << msgbase.sent_by << msgbase.dest;
       
    return ds;
}

/** Extract from a binary datastream */
QDataStream& operator>>(QDataStream &ds, MessageBase &msgbase)
{
    VersionID v = readHeader(ds, r_msgbase);
    
    if (v == 1)
    {
        ds >> msgbase.uid >> msgbase.subject_uid
           >> msgbase.sent_by >> msgbase.dest;
    }
    else
        throw version_error( v, "1", r_msgbase, CODELOC );
        
    return ds;
}

/** Null constructor */
MessageBase::MessageBase()
            : QSharedData(), sent_by(-1), dest(-1)
{}

static void assertValidDestination(int rank)
{
    if (rank != -1)
    {
        if (rank < 0 or rank >= MPICluster::getCount())
            throw SireError::unavailable_resource( QObject::tr(
                "There is no MPI process with rank %1. Available ranks are "
                "from 0 to %2.")
                    .arg(rank).arg(MPICluster::getCount()-1), CODELOC );
    }
}

/** Construct a message that is intended to go to the MPI process
    with rank 'destination' */
MessageBase::MessageBase(int destination) 
            : QSharedData(), 
              uid( QUuid::createUuid() ),
              subject_uid( QUuid::createUuid() ),
              sent_by( MPICluster::getRank() ),
              dest(destination)
{
    assertValidDestination(dest);
}

/** Construct a message that is intended to go to the MPI process
    with rank 'destination', and that is a response to the subject
    with UID 'uid' */
MessageBase::MessageBase(int destination, const QUuid &subuid)
            : QSharedData(),
              uid( QUuid::createUuid() ),
              subject_uid(subuid),
              sent_by( MPICluster::getRank() ),
              dest(destination)
{
    assertValidDestination(dest);
}

/** Copy constructor */
MessageBase::MessageBase(const MessageBase &other)
            : QSharedData(),
              uid(other.uid),
              subject_uid(other.subject_uid),
              sent_by(other.sent_by), dest(other.dest)
{}

/** Destructor */
MessageBase::~MessageBase()
{}

/** Copy assignment operator */
MessageBase& MessageBase::operator=(const MessageBase &other)
{
    if (this != &other)
    {
        uid = other.uid;
        subject_uid = other.subject_uid;
        sent_by = other.sent_by;
        dest = other.dest;
    }
    
    return *this;
}

/** Return a string representation of this message */
QString MessageBase::toString() const
{
    return this->what();
}

/** Return whether or not this message has a reply */
bool MessageBase::hasReply() const
{
    return false;
}

/** Return the reply to this message (which will be null, if hasReply() is false) */
Message MessageBase::reply() const
{
    return Message();
}

/** Return the unique ID of this message */
const QUuid& MessageBase::UID() const
{
    return uid;
}

/** Return the unique ID of the subject to which this 
    message relates */
const QUuid& MessageBase::subjectUID() const
{
    return subject_uid;
}

/** Return the rank of the MPI process that sent this message */
int MessageBase::sender() const
{
    return sent_by;
}

/** Return the rank of the MPI process to which this message
    should be sent */
int MessageBase::destination() const
{
    return dest;
}

/** Return whether this is one of the recipients of this message
    (and is thus expected to read the message!) */
bool MessageBase::isRecipient(int rank) const
{
    return (dest == -1) or (rank == dest);
}

/** Return the set of recipients for this message. This is not the same
    as the destination - the destination is where the message should be
    sent (which can either be a specific process, or broadcast to all processes),
    while the recipients are the actual processes that are expected to read
    the message, rather than just receive it. */
QSet<int> MessageBase::recipients() const
{
    QSet<int> ranks;
    
    if (dest == -1)
    {
        int nprocs = MPICluster::getCount();
        ranks.reserve(nprocs);
        
        for (int i=0; i<nprocs; ++i)
        {
            ranks.insert(i);
        }
    }
    else
        ranks.insert(dest);
        
    return ranks;
}

/////////
///////// Implementation of Message
/////////

static const RegisterMetaType<Message> r_message;
                                                      
/** Serialise to a binary datastream */
QDataStream& operator<<(QDataStream &ds, const Message &message)
{
    writeHeader(ds, r_message, 1);
    
    SharedDataStream sds(ds);
    sds << message.d;
       
    return ds;
}

/** Extract from a binary datastream */
QDataStream& operator>>(QDataStream &ds, Message &message)
{
    VersionID v = readHeader(ds, r_message);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> message.d;
    }
    else
        throw version_error( v, "1", r_message, CODELOC );
        
    return ds;
}

/** Construct a null message */
Message::Message()
{}

/** Construct from the message 'message' */
Message::Message(const MessageBase &message)
        : d(message)
{}

/** Copy constructor */
Message::Message(const Message &other)
        : d(other.d)
{}

/** Destructor */
Message::~Message()
{}

/** Copy assignment operator */
Message& Message::operator=(const Message &other)
{
    d = other.d;
    return *this;
}

/** Return whether or not this is a null message */
bool Message::isNull() const
{
    return d.constData() == 0;
}

/** Return a string representation of this message */
QString Message::toString() const
{
    if (this->isNull())
        return "Message::null";
    else
        return d->toString();
}

/** Internal function used to return the base object of the message */
const MessageBase& Message::base() const
{
    return *d;
}

/** Pack this message into a QByteArray */
QByteArray Message::pack() const
{
    QByteArray data;
    QDataStream ds(&data, QIODevice::WriteOnly);
    ds << *this;
    
    return data;
}

/** Load and return a message from the passed data */
Message Message::unpack(const QByteArray &data)
{
    if (data.isEmpty())
    {
        return Message();
    }
    else
    {
        QDataStream ds(data);
    
        Message message;
        ds >> message;
    
        return message;
    }
}

/** Read the message - this performs any actions contained therein */
void Message::read()
{
    if (not this->isNull())
    {
        d->read();
    }
}

/** Return whether or not this message has a reply */
bool Message::hasReply() const
{
    if (not this->isNull())
    {
        return d->hasReply();
    }
    else
        return false;
}

/** Return any reply associated with this message
    (a null Message will be returned if there is no reply) */
Message Message::reply() const
{
    if (not this->isNull())
    {
        return d->reply();
    }
    else
        return Message();
}

static QUuid null_uid;    

/** Return the unique ID of this message */
const QUuid& Message::UID() const
{
    if (not this->isNull())
    {
        return d->UID();
    }
    else
        return null_uid;
}

/** Return the unique ID of the subject of this message */
const QUuid& Message::subjectUID() const
{
    if (not this->isNull())
    {
        return d->subjectUID();
    }
    else
        return null_uid;
}

/** Return the rank of the MPI process that originally
    sent this message (-1 if this is not known) */
int Message::sender() const
{
    if (not this->isNull())
    {
        return d->sender();
    }
    else
        return -1;
}

/** Return the rank of the MPI process to which this
    message should be sent (-1 if it is to go to everyone) */
int Message::destination() const
{
    if (not this->isNull())
    {
        return d->destination();
    }
    else
        return -1;
}

/** Return whether this is one of the recipients of this message
    (and is thus expected to read the message!) */
bool Message::isRecipient(int rank) const
{
    if (not this->isNull())
    {
        return d->isRecipient(rank);
    }
    else
        return false;
}

/** Return the set of recipients for this message. This is not the same
    as the destination - the destination is where the message should be
    sent (which can either be a specific process, or broadcast to all processes),
    while the recipients are the actual processes that will really read
    the message, rather than just receive it. */
QSet<int> Message::recipients() const
{
    if (not this->isNull())
    {
        return d->recipients();
    }
    else
        return QSet<int>();
}

/////////
///////// Implementation of Broadcast
/////////

static const RegisterMetaType<Broadcast> r_broadcast;
                                                      
/** Serialise to a binary datastream */
QDataStream& operator<<(QDataStream &ds, const Broadcast &broadcast)
{
    writeHeader(ds, r_broadcast, 1);
    
    SharedDataStream sds(ds);
    sds << broadcast.destinations
        << broadcast.message_type
        << broadcast.message_data
        << broadcast.replymsg
        << static_cast<const MessageBase&>(broadcast);
       
    return ds;
}

/** Extract from a binary datastream */
QDataStream& operator>>(QDataStream &ds, Broadcast &broadcast)
{
    VersionID v = readHeader(ds, r_broadcast);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> broadcast.destinations
            >> broadcast.message_type
            >> broadcast.message_data
            >> broadcast.replymsg
            >> static_cast<MessageBase&>(broadcast);
    }
    else
        throw version_error( v, "1", r_broadcast, CODELOC );
        
    return ds;
}

/** Constructor */
Broadcast::Broadcast() : MessageBase()
{}

/** Construct to broadcast the message 'message' to
    every single MPI process */
Broadcast::Broadcast(const Message &message)
          : MessageBase(-1)  // -1 means broadcast to all
{
    if (message.isA<Broadcast>())
    {
        this->operator=( message.asA<Broadcast>() );
    }
    else
    {
        message_type = message.what();
        message_data = message.pack();
    }
}       

/** Construct to broadcast this message to just the MPI 
    process with rank 'rank' */
Broadcast::Broadcast(const Message &message, int rank)
          : MessageBase(rank),
            message_data( message.pack() )
{
    assertValidDestination(rank);

    if (rank != -1) // -1 means broadcast to everyone
    {
        destinations.insert(rank);
    }
    
    if (message.isA<Broadcast>())
    {
        message_data = message.asA<Broadcast>().message_data;
    }
    else
    {
        message_type = message.what();
        message_data = message.pack();
    }
}

/** Construct to broadcast this message to the MPI processes
    with ranks that are in 'ranks' */
Broadcast::Broadcast(const Message &message, const QSet<qint32> &ranks)
          : MessageBase(-1),
            destinations(ranks),
            message_data( message.pack() )
{
    foreach (qint32 rank, ranks)
    {
        assertValidDestination(rank);
    }

    if (destinations.contains(-1))
    {
        //this should be sent to everyone
        destinations = QSet<qint32>();
    }

    if (message.isA<Broadcast>())
    {
        message_data = message.asA<Broadcast>().message_data;
    }
    else
    {
        message_type = message.what();
        message_data = message.pack();
    }
}

/** Copy constructor */
Broadcast::Broadcast(const Broadcast &other)
          : MessageBase(other),
            destinations(other.destinations),
            message_type(other.message_type),
            message_data(other.message_data), 
            replymsg(other.replymsg)
{}

/** Destructor */
Broadcast::~Broadcast()
{}

/** Copy assignment operator */
Broadcast& Broadcast::operator=(const Broadcast &other)
{
    if (this != &other)
    {
        destinations = other.destinations;
        message_type = other.message_type;
        message_data = other.message_data;
        replymsg = other.replymsg;
        MessageBase::operator=(other);
    }
    
    return *this;
}

/** Return whether or not the process with MPI rank 'rank' 
    is an intended recipient of this broadcast */
bool Broadcast::isRecipient(int rank) const
{
    if (destinations.isEmpty())
    {
        return MessageBase::isRecipient(rank);
    }
    else
        return destinations.contains(rank);
}

/** Return the intended recipients of this broadcast */
QSet<int> Broadcast::recipients() const
{
    if (destinations.isEmpty())
    {
        return MessageBase::recipients();
    }
    else
    {
        QSet<int> procs;
        procs.reserve(destinations.count());
        
        foreach (qint32 dest, destinations)
        {
            procs.insert(dest);
        }
        
        return procs;
    }
}

/** Return a string representation of this message */
QString Broadcast::toString() const
{
    return QString( "Broadcast( %1 )" )
                .arg( Message::unpack(message_data).toString() );
}

/** Return the type of the message being broadcast */
const QString& Broadcast::messageType() const
{
    return message_type;
}

/** Read this message - this will only read this message if
    it is really intended for this node */
void Broadcast::read()
{
    if ( this->isRecipient(MPICluster::getRank()) )
    {
        //yes - we are an intended recipient
        Message message = Message::unpack(message_data);
        
        //read the real message
        message.read();
        
        //save any reply
        if (message.hasReply())
        {
            replymsg = message.reply();
        }
    }
}

/** Return whether or not the message resulted in a reply */
bool Broadcast::hasReply() const
{
    return not replymsg.isNull();
}

/** Return the reply that resulted from this message */
Message Broadcast::reply() const
{
    return replymsg;
}

/////////
///////// Implementation of RegisterBackend
/////////

static const RegisterMetaType<RegisterBackend> r_regbackend;
                                                      
/** Serialise to a binary datastream */
QDataStream& operator<<(QDataStream &ds, const RegisterBackend &regbackend)
{
    writeHeader(ds, r_regbackend, 1);
    
    ds << regbackend.node_uid << regbackend.node_rank
       << static_cast<const MessageBase&>(regbackend);
       
    return ds;
}

/** Extract from a binary datastream */
QDataStream& operator>>(QDataStream &ds, RegisterBackend &regbackend)
{
    VersionID v = readHeader(ds, r_regbackend);
    
    if (v == 1)
    {
        ds >> regbackend.node_uid >> regbackend.node_rank
           >> static_cast<MessageBase&>(regbackend);
    }
    else
        throw version_error( v, "1", r_regbackend, CODELOC );
        
    return ds;
}

/** Constructor */
RegisterBackend::RegisterBackend() : node_rank(-1)
{}

/** Construct the message to register the node with unique
    ID 'node_uid' that resides on this MPI process */
RegisterBackend::RegisterBackend(const QUuid &uid)
                : MessageBase( MPICluster::master() ),  // must be sent to the master
                  node_uid(uid),
                  node_rank( MPICluster::getRank() )
{}

/** Copy constructor */
RegisterBackend::RegisterBackend(const RegisterBackend &other)
                : MessageBase(other),
                  node_uid(other.node_uid),
                  node_rank(other.node_rank)
{}

/** Destructor */
RegisterBackend::~RegisterBackend()
{}

/** Copy assignment operator */
RegisterBackend& RegisterBackend::operator=(const RegisterBackend &other)
{
    if (this != &other)
    {
        node_uid = other.node_uid;
        node_rank = other.node_rank;
        MessageBase::operator=(other);
    }
    
    return *this;
}

/** Return a string representation of this message */
QString RegisterBackend::toString() const
{
    return QString("RegisterBackend( uid=%1, process=%2 )")
                .arg(node_uid.toString())
                .arg(node_rank);
}

/** Read this message - this must only occur on the master process! */
void RegisterBackend::read()
{
    if (not MPICluster::isMaster())
        throw SireError::program_bug( QObject::tr(
            "Only the master MPI process is allowed to read a RegisterBackend message."),
                CODELOC );

    //register the node
    MPICluster::registerBackend(node_rank, node_uid);
}

/////////
///////// Implementation of ReserveBackend
/////////

static const RegisterMetaType<ReserveBackend> r_reservebackend;
                                                      
/** Serialise to a binary datastream */
QDataStream& operator<<(QDataStream &ds, const ReserveBackend &reservebackend)
{
    writeHeader(ds, r_reservebackend, 1);
    
    ds << reservebackend.backend_uid << reservebackend.nbackends
       << static_cast<const MessageBase&>(reservebackend);
       
    return ds;
}

/** Extract from a binary datastream */
QDataStream& operator>>(QDataStream &ds, ReserveBackend &reservebackend)
{
    VersionID v = readHeader(ds, r_reservebackend);
    
    if (v == 1)
    {
        ds >> reservebackend.backend_uid >> reservebackend.nbackends
           >> static_cast<MessageBase&>(reservebackend);
    }
    else
        throw version_error( v, "1", r_reservebackend, CODELOC );
        
    return ds;
}

/** Constructor */
ReserveBackend::ReserveBackend() 
               : MessageBase( MPICluster::master() ), // must be sent to the master
                 nbackends(0)
{}

/** Construct the message to request that 'n' backends are reserved
    for use by this process */
ReserveBackend::ReserveBackend(int n) 
               : MessageBase( MPICluster::master() ),  // must be sent to the master
                 nbackends(n)
{}

/** Construct the message to request that the backend with ID 'uid'
    is reserved for use by this process */
ReserveBackend::ReserveBackend(const QUuid &uid)
               : MessageBase( MPICluster::master() ), // must be sent to the master
                 backend_uid(uid),
                 nbackends(1)
{
    if (uid.isNull())
        nbackends = 0;
}

/** Copy constructor */
ReserveBackend::ReserveBackend(const ReserveBackend &other)
               : MessageBase(other),
                 backend_uid(other.backend_uid),
                 nbackends(other.nbackends)
{}

/** Destructor */
ReserveBackend::~ReserveBackend()
{}

/** Copy assignment operator */
ReserveBackend& ReserveBackend::operator=(const ReserveBackend &other)
{
    if (this != &other)
    {
        backend_uid = other.backend_uid;
        nbackends = other.nbackends;
        MessageBase::operator=(other);
    }
    
    return *this;
}

/** Return a string representation of this message */
QString ReserveBackend::toString() const
{
    if (backend_uid.isNull())
        return QObject::tr("ReserveBackend( nBackends()=%1 )").arg(nbackends);
    else
        return QObject::tr("ReserveBackend( UID()=%1 )").arg( backend_uid.toString() );
}

/** Return the requested UID - this will be null if we don't
    care which backend we get */
const QUuid& ReserveBackend::requestedUID() const
{
    return backend_uid;
}

/** Return the requested number of nodes */
int ReserveBackend::nBackends() const
{
    return nbackends;
}

/** Read this message - this must only occur on the master process! */
void ReserveBackend::read()
{
    //tell the reservation manager to process this request
    ReservationManager::reserveBackends(*this);
}

/////////
///////// Implementation of RequestAvailability
/////////

static const RegisterMetaType<RequestAvailability> r_requestavail;
                                                      
/** Serialise to a binary datastream */
QDataStream& operator<<(QDataStream &ds, const RequestAvailability &requestavail)
{
    writeHeader(ds, r_requestavail, 1);
    
    ds << requestavail.request
       << static_cast<const MessageBase&>(requestavail);
       
    //no need to save available_backends
       
    return ds;
}

/** Extract from a binary datastream */
QDataStream& operator>>(QDataStream &ds, RequestAvailability &requestavail)
{
    VersionID v = readHeader(ds, r_requestavail);
    
    if (v == 1)
    {
        ds >> requestavail.request
           >> static_cast<MessageBase&>(requestavail);
    }
    else
        throw version_error( v, "1", r_requestavail, CODELOC );
        
    return ds;
}

/** Constructor */
RequestAvailability::RequestAvailability() 
                    : MessageBase(-1) // will be broadcast to everyone
{}

/** Construct a request to find out what is available that matches the
    request 'request' */
RequestAvailability::RequestAvailability(const ReserveBackend &reservation_request) 
                    : MessageBase(-1),  // will be broadcast to everyone
                      request(reservation_request)
{}

/** Copy constructor */
RequestAvailability::RequestAvailability(const RequestAvailability &other)
                    : MessageBase(other),
                      request(other.request),
                      available_backends(other.available_backends)
{}

/** Destructor */
RequestAvailability::~RequestAvailability()
{}

/** Copy assignment operator */
RequestAvailability& RequestAvailability::operator=(const RequestAvailability &other)
{
    if (this != &other)
    {
        request = other.request;
        available_backends = other.available_backends;
        MessageBase::operator=(other);
    }
    
    return *this;
}

/** Return a string representation of this message */
QString RequestAvailability::toString() const
{
    return QString( "RequestAvailability( %1 )" )
                    .arg(request.toString());
}

/** This message always has a reply */
bool RequestAvailability::hasReply() const
{
    return true;
}

/** Return the message containing the reply to this message */
Message RequestAvailability::reply() const
{
    return Result( *this, available_backends );
}

/** Read this message - finds all matching backends on this process
    and reserves them in case they are selected */
void RequestAvailability::read()
{
    available_backends = ReservationManager::findAvailable(request);
}

/////////
///////// Implementation of Reservation
/////////

static const RegisterMetaType<Reservation> r_reservation;
                                                      
/** Serialise to a binary datastream */
QDataStream& operator<<(QDataStream &ds, const Reservation &reservation)
{
    writeHeader(ds, r_reservation, 1);
    
    ds << reservation.request;
    
    quint32 nbackends = reservation.details.count();
    
    ds << nbackends;
    
    for (quint32 i=0; i<nbackends; ++i)
    {
        ds << qint32(reservation.details[i].get<0>())
           << reservation.details[i].get<1>();
    }
    
    ds << static_cast<const MessageBase&>(reservation);
       
    return ds;
}

/** Extract from a binary datastream */
QDataStream& operator>>(QDataStream &ds, Reservation &reservation)
{
    VersionID v = readHeader(ds, r_reservation);
    
    if (v == 1)
    {
        ds >> reservation.request;

        quint32 nbackends;
        
        ds >> nbackends;
        
        reservation.details.clear();
        
        for (quint32 i=0; i<nbackends; ++i)
        {
            qint32 rank;
            QUuid uid;
            
            ds >> rank >> uid;
            
            reservation.details.append( tuple<int,QUuid>(rank,uid) );
        }

        ds >> static_cast<MessageBase&>(reservation);
    }
    else
        throw version_error( v, "1", r_reservation, CODELOC );
        
    return ds;
}

/** Constructor */
Reservation::Reservation() 
            : MessageBase(-1) // will be broadcast to everyone
{}

/** Construct the message to be broadcast that instructs that the passed
    backends must now make the connections as requested in 'request' */
Reservation::Reservation(const QList< tuple<int,QUuid> > &reserved_backends,
                         const ReserveBackend &reservation_request)
            : MessageBase(-1),  // will be broadcast to everyone
              request(reservation_request),
              details(reserved_backends)
{}

/** Copy constructor */
Reservation::Reservation(const Reservation &other)
                    : MessageBase(other),
                      request(other.request),
                      details(other.details)
{}

/** Destructor */
Reservation::~Reservation()
{}

/** Copy assignment operator */
Reservation& Reservation::operator=(const Reservation &other)
{
    if (this != &other)
    {
        request = other.request;
        details = other.details;
        MessageBase::operator=(other);
    }
    
    return *this;
}

/** Return a string representation of this message */
QString Reservation::toString() const
{
    return QString( "Reservation( %1 : nMatches()=%2 )" )
                    .arg(request.toString()).arg(details.count());
}

/** This message always has a reply */
bool Reservation::hasReply() const
{
    return true;
}

/** Return the message containing the reply to this message */
Message Reservation::reply() const
{
    return Result( *this, ReservationManager::establishedConnections(request) );
}

/** Read this message - finds all matching backends on this process
    and reserves them in case they are selected */
void Reservation::read()
{
    ReservationManager::establishConnections( details, request );
}

/////////
///////// Implementation of Shutdown
/////////

static const RegisterMetaType<Shutdown> r_shutdown;
                                                      
/** Serialise to a binary datastream */
QDataStream& operator<<(QDataStream &ds, const Shutdown &shutdown)
{
    writeHeader(ds, r_shutdown, 1);
    
    ds << static_cast<const MessageBase&>(shutdown);
       
    return ds;
}

/** Extract from a binary datastream */
QDataStream& operator>>(QDataStream &ds, Shutdown &shutdown)
{
    VersionID v = readHeader(ds, r_shutdown);
    
    if (v == 1)
    {
        ds >> static_cast<MessageBase&>(shutdown);
    }
    else
        throw version_error( v, "1", r_shutdown, CODELOC );
        
    return ds;
}

/** Construct a message to tell everyone to shut down */
Shutdown::Shutdown() : MessageBase(-1)  // goes to everyone
{}

/** Copy constructor */
Shutdown::Shutdown(const Shutdown &other) : MessageBase(other)
{}

/** Destructor */
Shutdown::~Shutdown()
{}

/** Copy assignment operator */
Shutdown& Shutdown::operator=(const Shutdown &other)
{
    MessageBase::operator=(other);
    return *this;
}

/** This is a simple class that provide one-shot
    running of functions in a background thread - this
    is necessary when shutting down the cluster, as otherwise
    we'll block waiting for the call to shut down the cluster
    to finish... */
class BG : private QThread
{
public:
    static void call( void (*func)() )
    {
        QMutexLocker lkr(&runmutex);
        QMutexLocker lkr1(&datamutex);
        runner.func = func;
        runner.start();
        runwaiter.wait(&datamutex);
        lkr1.unlock();
        lkr.unlock();
    }
    
private:
    BG()
    {}
    
    ~BG()
    {}

    void run()
    {
        SireError::setThreadString("BG");
        SireMaths::seed_qrand();
    
        datamutex.lock();
        void (*f)() = func;
        runwaiter.wakeAll();
        datamutex.unlock();
        
        (*f)();
    }

    void (*func)();
    
    static QMutex runmutex;
    static QMutex datamutex;
    static QWaitCondition runwaiter;
    static BG runner;
};

QMutex BG::runmutex;
QMutex BG::datamutex;
QWaitCondition BG::runwaiter;
BG BG::runner;

/** Read this message */
void Shutdown::read()
{
    BG::call( &MPICluster::informedShutdown );
}

/////////
///////// Implementation Error
/////////

static const RegisterMetaType<Error> r_error;

/** Serialise to a binary datastream */
QDataStream& operator<<(QDataStream &ds, const Error &error)
{
    writeHeader(ds, r_error, 1);
    
    ds << error.message_data << error.error_data
       << static_cast<const MessageBase&>(error);
       
    return ds;
}

/** Extract from a binary datastream */
QDataStream& operator>>(QDataStream &ds, Error &error)
{
    VersionID v = readHeader(ds, r_error);
    
    if (v == 1)
    {
        ds >> error.message_data >> error.error_data
           >> static_cast<MessageBase&>(error);
    }
    else
        throw version_error(v, "1", r_error, CODELOC);
        
    return ds;
}

void Error::setError(const SireError::exception &e)
{
    error_data = e.pack();
}

void Error::setUnknownError(const QString &code_loc)
{
    this->setError( SireError::unknown_exception( QObject::tr(
                            "An unknown error occured!"), code_loc ) );
}

void Error::setError(const std::exception &e)
{
    this->setError( SireError::std_exception(e) );
}

/** Construct a message containing an unknown error 
    in response to an unknown sender */
Error::Error() : MessageBase( MPICluster::master() )
{
    this->setUnknownError(CODELOC);
}

/** Construct a message containing an unknown error that occured 
    at the point 'code_loc' in response to an unknown sender */
Error::Error(const QString &code_loc) : MessageBase( MPICluster::master() )
{
    this->setUnknownError(code_loc);
}

/** Construct a message containing the error 'e'
    in response to an unknown sender */
Error::Error(const SireError::exception &e) : MessageBase( MPICluster::master() )
{
    this->setError(e);
}

/** Construct a message containing the error 'e'
    in response to an unknown sender */
Error::Error(const std::exception &e) : MessageBase( MPICluster::master() )
{
    this->setError(e);
}

/** Construct a message containing an unknown error
    in response to the sender that has rank 'sender' */
Error::Error(int sender) 
      : MessageBase( sender )
{
    this->setUnknownError(CODELOC);
}

/** Construct a message containing an unknown error that occured
    at the point 'code_loc' in response to the sender that has rank 'sender' */
Error::Error(int sender, const QString &code_loc) 
      : MessageBase( sender )
{
    this->setUnknownError(code_loc);
}

/** Construct a message containing the error 'e'
    in response to the sender that has rank 'sender' */
Error::Error(int sender, const SireError::exception &e) 
      : MessageBase( sender )
{
    this->setError(e);
}

/** Construct a message containing the error 'e'
    in response to the sender that has rank 'sender' */
Error::Error(int sender, const std::exception &e) 
      : MessageBase( sender )
{
    this->setError(e);
}

/** Construct a message containing an unknown error
    in response to the message 'message' */
Error::Error(const Message &message) 
      : MessageBase( message.sender(), message.subjectUID() )
{
    message_data = message.pack();
    this->setUnknownError(CODELOC);
}

/** Construct a message containing an unknown error that occured at
    the point 'code_loc' in response to the message 'message' */
Error::Error(const Message &message, const QString &code_loc) 
      : MessageBase( message.sender(), message.subjectUID() )
{
    message_data = message.pack();
    this->setUnknownError(code_loc);
}

/** Construct a message containing the error 'e'
    in response to the message 'message' */
Error::Error(const Message &message, const SireError::exception &e) 
      : MessageBase( message.sender(), message.subjectUID() )
{
    message_data = message.pack();
    this->setError(e);
}

/** Construct a message containing the error 'e'
    in response to the message 'message' */
Error::Error(const Message &message, const std::exception &e) 
      : MessageBase( message.sender(), message.subjectUID() )
{
    message_data = message.pack();
    this->setError(e);
}

/** Construct an error that will just pass the message and error
    data on to the master - this is used when an error that can't 
    be handled is received */
Error::Error(const QByteArray &message, const QByteArray &error)
      : MessageBase( MPICluster::master() ),
        message_data(message), error_data(error)
{}

/** Copy constructor */
Error::Error(const Error &other) 
             : MessageBase(other),
               message_data(other.message_data), error_data(other.error_data)
{}

/** Destructor */
Error::~Error()
{}

/** Copy assignment operator */
Error& Error::operator=(const Error &other)
{
    if (this != &other)
    {
        message_data = other.message_data;
        error_data = other.error_data;
        MessageBase::operator=(other);
    }
    
    return *this;
}

/** Read this message */
void Error::read()
{
    MPICluster::postError( this->subjectUID(), this->sender(),
                           message_data, error_data );
}

/////////
///////// Implementation Result
/////////

static const RegisterMetaType<Result> r_result;

/** Serialise to a binary datastream */
QDataStream& operator<<(QDataStream &ds, const Result &result)
{
    writeHeader(ds, r_result, 1);
    
    ds << result.result_data
       << static_cast<const MessageBase&>(result);
       
    return ds;
}

/** Extract from a binary datastream */
QDataStream& operator>>(QDataStream &ds, Result &result)
{
    VersionID v = readHeader(ds, r_result);
    
    if (v == 1)
    {
        ds >> result.result_data
           >> static_cast<MessageBase&>(result);
    }
    else
        throw version_error(v, "1", r_result, CODELOC);
        
    return ds;
}

/** Construct a message containing an unknown result 
    in response to an unknown sender */
Result::Result() : MessageBase( MPICluster::master() )
{}

/** Copy constructor */
Result::Result(const Result &other) 
             : MessageBase(other),
               result_data(other.result_data)
{}

/** Destructor */
Result::~Result()
{}

/** Copy assignment operator */
Result& Result::operator=(const Result &other)
{
    if (this != &other)
    {
        result_data = other.result_data;
        MessageBase::operator=(other);
    }
    
    return *this;
}

/** Read this message */
void Result::read()
{
    MPICluster::postResult( this->subjectUID(), this->sender(),
                            result_data );
}

/////////
///////// Implementation GetUIDs
/////////

static const RegisterMetaType<GetUIDs> r_getuids;

/** Serialise to a binary datastream */
QDataStream& operator<<(QDataStream &ds, const GetUIDs &getuids)
{
    writeHeader(ds, r_getuids, 1);
    
    ds << static_cast<const MessageBase&>(getuids);
       
    return ds;
}

/** Extract from a binary datastream */
QDataStream& operator>>(QDataStream &ds, GetUIDs &getuids)
{
    VersionID v = readHeader(ds, r_getuids);
    
    if (v == 1)
    {
        ds >> static_cast<MessageBase&>(getuids);
    }
    else
        throw version_error(v, "1", r_getuids, CODELOC);
        
    return ds;
}

/** Construct a message that asks the master to return the UIDs
    of all of the backends */
GetUIDs::GetUIDs() : MessageBase( MPICluster::master() )
{}

/** Copy constructor */
GetUIDs::GetUIDs(const GetUIDs &other) : MessageBase(other)
{}

/** Destructor */
GetUIDs::~GetUIDs()
{}

/** Copy assignment operator */
GetUIDs& GetUIDs::operator=(const GetUIDs &other)
{
    MessageBase::operator=(other);
    return *this;
}

/** This message does have a reply */
bool GetUIDs::hasReply() const
{
    return true;
}

/** Get the reply to this message */
Message GetUIDs::reply() const
{
    QList<QUuid> uids = MPICluster::UIDs();

    return Result(*this, uids);
}

/** Read this message */
void GetUIDs::read()
{
    //ensure that this is the master process
    if (not MPICluster::isMaster())
    {
        throw SireError::program_bug( QObject::tr(
            "A GetUIDs message can only be read by the master process!"),
                CODELOC );
    }
}

const char* Broadcast::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Broadcast>() );
}

const char* Result::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Result>() );
}

const char* Error::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Error>() );
}

const char* RegisterBackend::typeName()
{
    return QMetaType::typeName( qMetaTypeId<RegisterBackend>() );
}

const char* ReserveBackend::typeName()
{
    return QMetaType::typeName( qMetaTypeId<ReserveBackend>() );
}

const char* RequestAvailability::typeName()
{
    return QMetaType::typeName( qMetaTypeId<RequestAvailability>() );
}

const char* Reservation::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Reservation>() );
}

const char* Shutdown::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Shutdown>() );
}

const char* GetUIDs::typeName()
{
    return QMetaType::typeName( qMetaTypeId<GetUIDs>() );
}

const char* Message::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Message>() );
}

Reservation* Reservation::clone() const
{
    return new Reservation(*this);
}

RequestAvailability* RequestAvailability::clone() const
{
    return new RequestAvailability(*this);
}

Result* Result::clone() const
{
    return new Result(*this);
}

GetUIDs* GetUIDs::clone() const
{
    return new GetUIDs(*this);
}

Shutdown* Shutdown::clone() const
{
    return new Shutdown(*this);
}

Broadcast* Broadcast::clone() const
{
    return new Broadcast(*this);
}

Error* Error::clone() const
{
    return new Error(*this);
}

RegisterBackend* RegisterBackend::clone() const
{
    return new RegisterBackend(*this);
}

ReserveBackend* ReserveBackend::clone() const
{
    return new ReserveBackend(*this);
}

#endif // SIRE_USE_MPI
