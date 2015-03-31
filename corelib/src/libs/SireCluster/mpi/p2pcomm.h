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

#ifndef SIRECLUSTER_MPI_P2PCOMM_H
#define SIRECLUSTER_MPI_P2PCOMM_H

#ifdef SIRE_USE_MPI

#include <mpi.h> // must be first to satisfy mpich

#include <QUuid>
#include <QByteArray>
#include <QDataStream>

#include "sireglobal.h"

#include <boost/shared_ptr.hpp>

SIRE_BEGIN_HEADER

namespace SireCluster
{

class Frontend;

namespace MPI
{

class MPICluster;

namespace detail
{
class P2PCommPvt;
}

/** This class holds a point-to-point communicator, which is
    used to allow direct, private communication between two
    processes in the MPI cluster
    
    @author Christopher Woods
*/
class P2PComm
{

friend class MPICluster;

public:
    enum { MASTER            = 0,  // the rank of the master
           SLAVE             = 1   // the rank of the slave
         };

    enum { EXIT              = 0,  // exit and drop the backend
           GETUID            = 1,  // get the UID of the backend
           START             = 2,  // start a new job
           STOP              = 3,  // stop a running job
           ABORT             = 4,  // abort a running job
           PROGRESS          = 5,  // get a progress update
           INTERIM           = 6,  // get an interim result
           RESULT            = 7,  // get the result
           OK                = 8,  // everything is OK
           IS_RUNNING        = 9
         };

    P2PComm();
    P2PComm(int master_rank, int slave_rank);
    
    P2PComm(const P2PComm &other);
    
    ~P2PComm();
    
    P2PComm& operator=(const P2PComm &other);
    
    bool involves(int rank);
    
    bool isMaster();
    bool isSlave();
    
    bool isLocal();
    
    bool isNull() const;
    
    void setBackend(const Frontend &backend);
    
    bool hasFinished();
    
    void sendMessage(int message);
    
    template<class T>
    void sendMessage(int message, const T &data);
    
    template<class T>
    T awaitResponse(bool urgent=false);
    
    int awaitIntegerResponse(bool urgent=false);
    float awaitFloatResponse(bool urgent=false);
    
protected:
    static P2PComm createLocal();
    static P2PComm create(MPI_Comm private_comm, 
                          int master_rank, int slave_rank);
    
private:
    void _pvt_sendMessage(int message, const QByteArray &data);
    
    QByteArray _pvt_awaitResponse(bool urgent);
    
    /** Shared pointer to the implementation */
    boost::shared_ptr<detail::P2PCommPvt> d;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Send the message 'message', which comes with the associated
    data 'data' */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void P2PComm::sendMessage(int message, const T &data)
{
    QByteArray bindata;
    QDataStream ds( &bindata, QIODevice::WriteOnly );
    ds << data;
    
    this->_pvt_sendMessage(message, bindata);
}
 
/** Wait for a response of type 'T', and return that 
    response */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
T P2PComm::awaitResponse(bool urgent)
{
    QByteArray bindata = this->_pvt_awaitResponse(urgent);
    
    if (bindata.isEmpty())
        return T();
    
    QDataStream ds(bindata);
    
    T response;
    ds >> response;
    
    return response;
}

#endif // SIRE_SKIP_INLINE_FUNCTIONS

} // end of namespace MPI

} // end of namespace SireCluster

SIRE_END_HEADER

#endif // SIRE_USE_MPI
#endif
