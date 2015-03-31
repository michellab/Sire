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

#ifndef SIRECLUSTER_CLUSTER_H
#define SIRECLUSTER_CLUSTER_H

#include "sireglobal.h"

#include <QUuid>
#include <QList>

#include <boost/shared_ptr.hpp>

SIRE_BEGIN_HEADER

namespace SireCluster
{

class Nodes;

class Backend;
class Frontend;

namespace MPI
{
class ReservationManager;
}

/** This class provides the global registry for all nodes in the cluster.
    A node is defined as a resource that can run a WorkPacket. A node
    consists of a Backend (the object in which the WorkPacket is 
    run) and a Frontend (the object that allows the node to communicate
    with the Backend)
    
    @author Christopher Woods
*/
class SIRECLUSTER_EXPORT Cluster
{

friend class Backend;
friend class SireCluster::MPI::ReservationManager;

public:
    static QList<QUuid> localUIDs();
    
    static QList<QUuid> UIDs();

    static bool isLocal(const QUuid &uid);

    static int getRank();
    static int getCount();
    
    static void start(int ppn=1);

    static void shutdown();

    static void wait();
    
    static bool isRunning();

    static bool supportsMPI();

    static Nodes getNode();
    static Nodes getNode(int timeout);

    static Nodes getNode(const QUuid &uid);
    static Nodes getNode(const QUuid &uid, int timeout);
    
    static Nodes getNodes(int nnodes);
    static Nodes getNodes(int nnodes, int timeout);
    
    static Nodes getNodes(const QList<QUuid> &uids);
    static Nodes getNodes(const QList<QUuid> &uids, int timeout);
    
    static Nodes getAllNodes();
    static Nodes getAllNodes(int timeout);

protected:    
    static void registerBackend(const Backend &backend);  // called by Backend
    
    static QList<Frontend> localBackends(); // call by ReservationManager

private:
    static Frontend _pvt_getFrontend();
    static Frontend _pvt_getFrontend(const QUuid &uid);

    static QList<Frontend> _pvt_getFrontends(int n);
    static QList<Frontend> _pvt_getFrontends(int n, const QUuid &uid);

    static Frontend getFrontend();
    static Frontend getFrontend(int timeout);

    static QList<Frontend> getFrontends(int n);
    static QList<Frontend> getFrontends(int n, int timeout);
    
    static Frontend getFrontend(const QUuid &uid);
    static Frontend getFrontend(const QUuid &uid, int timeout);
};

}

SIRE_EXPOSE_CLASS( SireCluster::Cluster )

SIRE_END_HEADER

#endif
