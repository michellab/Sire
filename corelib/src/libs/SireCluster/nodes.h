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

#ifndef SIRECLUSTER_NODES_H
#define SIRECLUSTER_NODES_H

#include "sireglobal.h"

#include <QUuid>

#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>

SIRE_BEGIN_HEADER

namespace SireCluster
{

class Frontend;
class Node;
class NodesPtr;

class ThisThread;

namespace detail
{
class NodePvt;
class NodesPvt;
class ThisThreadPvt;
}

/** This class holds, and schedules, a collection of Node objects.
    
    @author Christopher Woods
*/
class SIRECLUSTER_EXPORT Nodes
{

friend class NodesPtr;
friend class Cluster;

friend class detail::ThisThreadPvt;

public:
    Nodes();
    
    Nodes(const Nodes &other);
    
    ~Nodes();
    
    Nodes& operator=(const Nodes &other);
    
    bool operator==(const Nodes &other) const;
    bool operator!=(const Nodes &other) const;
    
    bool isEmpty();
    
    QString toString() const;
    
    Node getNode();
    QList<Node> getNodes(int n);
    QList<Node> getAllNodes();
    
    Node getNode(int timeout);
    QList<Node> getNodes(int n, int timeout);
    QList<Node> getAllNodes(int timeout);
    
    void waitUntilAllFree();
    bool waitUntilAllFree(int timeout);
    
    int nFree();
    int nBusy();
    int nNodes();
    int count();

    void add(Node node);
    void remove(Node node);

    void add(Nodes &nodes);

    void addNode();
    void addNode(int timeout);
    
    void addNodes(int n);
    void addNodes(int n, int timeout);

    void removeAll();

    ThisThread borrowThisThread();

protected:
    Nodes(const boost::shared_ptr<detail::NodesPvt> &ptr); // called by NodesPvt

    Nodes(Frontend frontend);                // called by Cluster
    Nodes(const QList<Frontend> &frontends); // called by Cluster

    QUuid createThisThread();                  // called by ThisThreadPvt
    void reclaimThisThread(const QUuid &uid);  // called by ThisThreadPvt

private:
    Node _pvt_getNode();

    /** Private implementation of the Nodes */
    boost::shared_ptr<detail::NodesPvt> d;
};

/** This class holds a weak pointer to the Nodes */
class NodesPtr
{

friend class detail::NodePvt;

public:
    NodesPtr();
    NodesPtr(const Nodes &nodes);
    
    NodesPtr(const NodesPtr &other);
    
    ~NodesPtr();
    
    NodesPtr& operator=(const NodesPtr &other);
    
    Nodes lock() const;
    Nodes operator*() const;

    bool expired() const;

    void reset();
    
protected:
    void returnFrontend(Frontend frontend);  // called by NodePvt

private:
    /** Weak pointer to the nodes */
    boost::weak_ptr<detail::NodesPvt> d;
};

/** This simple class is used to allow the current thread
    to be made available to a Nodes object
    
    @author Christopher Woods
*/
class SIRECLUSTER_EXPORT ThisThread
{
public:
    ThisThread();
    ThisThread(const ThisThread &other);
    
    ThisThread(const Nodes &nodes);
    
    ThisThread& operator=(const ThisThread &other);
    
    ~ThisThread();
    
    void reclaim();

private:
    /** Pointer to the implementation */
    boost::shared_ptr<detail::ThisThreadPvt> d;
};

}

SIRE_EXPOSE_CLASS( SireCluster::Nodes )
SIRE_EXPOSE_CLASS( SireCluster::ThisThread )

SIRE_END_HEADER

#endif
