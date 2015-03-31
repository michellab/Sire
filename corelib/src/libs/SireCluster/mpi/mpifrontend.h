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

#ifndef SIRECLUSTER_MPI_MPIFRONTEND_H
#define SIRECLUSTER_MPI_MPIFRONTEND_H

#ifdef SIRE_USE_MPI

#include <QMutex>

#include "p2pcomm.h"

#include "SireCluster/frontend.h"

SIRE_BEGIN_HEADER

namespace SireCluster
{
namespace MPI
{

/** This is a Frontend that is specialised to communicate with 
    a backend over an MPI connection
    
    @author Christopher Woods
*/
class MPIFrontend : public FrontendBase
{
public:
    MPIFrontend();
    MPIFrontend(const P2PComm &p2pcomm);
    
    ~MPIFrontend();
    
    bool isLocal() const;
    
    QUuid UID();
    
    void startJob(const WorkPacket &workpacket);
    
    void stopJob();
    void abortJob();
    
    void wait();
    bool wait(int timeout);
    
    float progress();
    WorkPacket interimResult();
    
    WorkPacket result();

private:
    /** A mutex used to protect access to the communicator */
    QMutex datamutex;

    /** The cached QUuid */
    QUuid cached_uid;

    /** The point-to-point communicator used to communicate
        with the remote backend */
    P2PComm p2p;
};

} // end of namespace MPI
} // end of namespace SireCluster

SIRE_END_HEADER

#endif // SIRE_USE_MPI
#endif
