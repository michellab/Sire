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

#include <QMutex>
#include <QSemaphore>
#include <QWaitCondition>
#include <QUuid>
#include <QThreadStorage>

#include <QHash>
#include <QSet>
#include <QList>

#include <QTime>

#include "nodes.h"
#include "node.h"
#include "frontend.h"
#include "backend.h"
#include "cluster.h"

#include "SireError/errors.h"
#include "SireError/printerror.h"

#include <QDebug>

using namespace SireCluster;

using boost::shared_ptr;
using boost::weak_ptr;

namespace SireCluster
{
namespace detail
{

/** Private implementation of Nodes */
class NodesPvt
{
public:
    NodesPvt()
    {}
    
    ~NodesPvt()
    {}
    
    /** Mutex to protect access to the main data
        (the busy and free queues) */
    QMutex datamutex;
    
    /** Pointer to a semaphore that is used to control
        the reservation and allocation of nodes */
    shared_ptr<QSemaphore> nodesem;
    
    /** WaitCondition used to wait until all of the nodes are free */
    QWaitCondition waiter;
    
    /** The collection of all non-null Frontends in this
        nodes object, indexed by the UID of the backend they
        are connected to */
    QHash<QUuid,Frontend> frontends;
    
    /** The set of UIDs of the busy frontends */
    QSet<QUuid> busy_frontends;
    
    /** The list of UIDs of the free frontends - local front
        ends will tend to be at the beginning of this list */
    QList<QUuid> free_frontends;
};

} // end of namespace detail
} // end of namespace SireCluster

using namespace SireCluster::detail;

//////////////
////////////// Implementation of Nodes
//////////////

/** Construct an empty set of nodes */
Nodes::Nodes()
{}

/** Construct to hold the node that connects to the backend
    using 'frontend' */
Nodes::Nodes(Frontend frontend)
{
    QUuid uid = frontend.UID();
    
    if (uid.isNull())
        return;
        
    d.reset( new NodesPvt() );
    
    QMutexLocker lkr( &(d->datamutex) );
    
    d->frontends.insert( uid, frontend );
    d->free_frontends.append(uid);
    
    d->nodesem.reset( new QSemaphore(d->frontends.count()) );
}

/** Construct to hold the nodes that connect to the backends
    using the frontends in 'frontends' */
Nodes::Nodes(const QList<Frontend> &frontends) : d( new NodesPvt() )
{
    ///// Add all of the frontends
    {
        QMutexLocker lkr( &(d->datamutex) );

        foreach (Frontend frontend, frontends)
        {
            QUuid uid = frontend.UID();
            
            if (uid.isNull())
                continue;
                
            d->frontends.insert( uid, frontend );
            
            if (frontend.isLocal())
            {
                d->free_frontends.prepend(uid);
            }
            else
            {
                d->free_frontends.append(uid);
            }
        }
        
        if (not d->frontends.isEmpty())
        {
            d->nodesem.reset( new QSemaphore(d->frontends.count()) );
        }
    }
    
    if (d->frontends.isEmpty())
    {
        //this is an empty set of Nodes
        d.reset();
    }
}

/** Construct from the passed pointer */
Nodes::Nodes(const shared_ptr<NodesPvt> &ptr) : d(ptr)
{}

/** Copy constructor - Nodes are explicitly shared */
Nodes::Nodes(const Nodes &other) : d(other.d)
{}

/** Destructor */
Nodes::~Nodes()
{
    if (d.unique())
    {
        this->removeAll();
    }
}

/** Copy assignment operator - Nodes are explicitly shared */
Nodes& Nodes::operator=(const Nodes &other)
{
    d = other.d;
    return *this;
}

/** Comparison operator */
bool Nodes::operator==(const Nodes &other) const
{
    return d.get() == other.d.get();
}

/** Comparison operator */
bool Nodes::operator!=(const Nodes &other) const
{
    return d.get() != other.d.get();
}

/** Return whether or not this is empty (contains no nodes) */
bool Nodes::isEmpty()
{
    if (d.get() == 0)
        return true;
        
    else
    {
        QMutexLocker lkr( &(d->datamutex) );
        return d->frontends.isEmpty();
    }
}

/** Return a string representation of these nodes */
QString Nodes::toString() const
{
    if (d.get() == 0)
    {
        return QObject::tr("Nodes( nBusy() == 0, nFree() == 0 )");
    }
    else
    {
        Nodes *nonconst_this = const_cast<Nodes*>(this);
    
        QStringList lines;
        
        QMutexLocker lkr( &(nonconst_this->d->datamutex) );
    
        lines.append( QObject::tr("Nodes( nBusy() == %1, nFree() == %2 )")
                            .arg(nonconst_this->d->busy_frontends.count())
                            .arg(nonconst_this->d->free_frontends.count()) );
                            
        if (nonconst_this->d->busy_frontends.count() > 0)
        {
            lines.append( QObject::tr("****** Busy nodes ******") );
            
            foreach (QUuid uid, nonconst_this->d->busy_frontends)
            {
                lines.append( QObject::tr("* Node %1").arg(uid.toString()) );
            }
        }

        if (nonconst_this->d->free_frontends.count() > 0)
        {
            lines.append( QObject::tr("****** Free nodes ******") );
            
            foreach (QUuid uid, nonconst_this->d->free_frontends)
            {
                lines.append( QObject::tr("* Node %1").arg(uid.toString()) );
            }
        }
        
        if (lines.count() > 1)
            lines.append( QObject::tr("************************") );
            
        return lines.join("\n");
    }
}

/** Internal function used to return a Node - you must be
    holding the datamutex to call this function */
Node Nodes::_pvt_getNode()
{
    QUuid uid = d->free_frontends.takeFirst();
    
    if (uid.isNull())
        throw SireError::program_bug( QObject::tr(
            "How have we managed to reserve a node, but there are no "
            "free nodes available???"), CODELOC );
        
    d->busy_frontends.insert(uid);
    
    Frontend frontend = d->frontends.value(uid);

    return Node::create(*this, frontend);
}

/** Return a free node - this blocks until a free node
    is available. In certain circumstances this will fail
    and return a null Node (e.g. at program shutdown, or
    if all nodes are removed from this set)
*/
Node Nodes::getNode()
{
    if (this->isEmpty())
        return Node();
    
    QMutexLocker lkr( &(d->datamutex) );
    shared_ptr<QSemaphore> nodesem = d->nodesem;

    bool reserved_node = false;
    
    //keep trying to reserve a node
    while (not reserved_node)
    {
        lkr.unlock();

        if (nodesem.get() == 0)
            //all of the nodes have been removed
            return Node();

        //reserve a node
        #if QT_VERSION >= 0x040300
        while (not nodesem->tryAcquire(1, 1000))
        #else
        while (not nodesem->tryAcquire(1))
        #endif
        {
            lkr.relock();
            nodesem = d->nodesem;
            lkr.unlock();
            
            if (nodesem.get() == 0)
                //all of the nodes have been removed
                return Node();

            #if QT_VERSION < 0x40300
            sleep(1);
            #endif
        }

        lkr.relock();
        
        //check that the semaphore hasn't changed while we were waiting
        if (d->nodesem.get() != nodesem.get())
        {
            //yes - it has changed - we didn't get the node
            nodesem->release();
            nodesem = d->nodesem;
        }
        else
            //we got a valid reservation
            reserved_node = true;
    }

    //now collect the reservation
    return this->_pvt_getNode();
}

/** Return 'n' free nodes - this blocks until all of the 
    nodes are available. In some circumstances this will fail,
    e.g. if there aren't enough nodes to fulfill this request
*/
QList<Node> Nodes::getNodes(int n)
{
    QList<Node> nodes;

    if (n <= 0 or this->isEmpty())
    {
        return nodes;
    }
    else if (n == 1)
    {
        nodes.append( this->getNode() );
    }
    else
    {
        QMutexLocker lkr( &(d->datamutex) );
        
        shared_ptr<QSemaphore> nodesem = d->nodesem;
        
        bool reserved_nodes = false;
        
        while (not reserved_nodes)
        {
            lkr.unlock();
            
            if (nodesem.get() == 0)
                //all of the nodes have been removed!
                return nodes;
                
            //reserve n nodes
            #if QT_VERSION >= 0x040300
            while (not nodesem->tryAcquire(n, 2000))
            #else
            while (not nodesem->tryAcquire(n))
            #endif
            {
                lkr.relock();
                nodesem = d->nodesem;
                lkr.unlock();
                
                if (nodesem.get() == 0)
                    //all of the nodes have been removed!
                    return nodes;
                    
                #if QT_VERSION < 0x40300
                sleep(2);
                #endif
            }
            
            lkr.relock();
            
            if (nodesem.get() != d->nodesem.get())
            {
                //the nodes available changed while we were making
                //this reservation - we've got to try again!
                nodesem->release(n);
                nodesem = d->nodesem;
            }
            else
                //we got the reservation!
                reserved_nodes = true;
        }

        //now grab all n nodes
        for (int i=0; i<n; ++i)
        {
            nodes.append( this-> _pvt_getNode() );
        }
    }
    
    return nodes;
}

/** Return all of the nodes - this blocks until all of the 
    nodes are available. Remember that this will return an
    empty list if there are no nodes.
*/
QList<Node> Nodes::getAllNodes()
{
    if (this->isEmpty())
        return QList<Node>();

    QMutexLocker lkr( &(d->datamutex) );
    
    shared_ptr<QSemaphore> nodesem = d->nodesem;
    
    bool reserved_nodes = false;
    
    while (not reserved_nodes)
    {
        lkr.unlock();
        
        if (nodesem.get() == 0)
            //all of the nodes have been removed!
            return QList<Node>();
            
        //reserve all nodes
        #if QT_VERSION >= 0x040300
        while (not nodesem->tryAcquire( nodesem->available(), 2000 ))
        #else
        while (not nodesem->tryAcquire( nodesem->available() ))
        #endif
        {
            lkr.relock();
            nodesem = d->nodesem;
            lkr.unlock();
            
            if (nodesem.get() == 0)
                //all of the nodes have been rmeoved
                return QList<Node>();
                
            #if QT_VERSION < 0x40300
            sleep(2);
            #endif
        }
        
        lkr.relock();
        
        if (nodesem.get() != d->nodesem.get())
        {
            //the nodes available changed while we were making
            //this reservation - we've got to try again!
            nodesem->release( nodesem->available() );
            nodesem = d->nodesem;
        }
        else
            //we got the reservation!
            reserved_nodes = true;
    }

    //now grab all n nodes
    QList<Node> nodes;

    int n = nodesem->available();

    for (int i=0; i<n; ++i)
    {
        nodes.append( this-> _pvt_getNode() );
    }
    
    return nodes;
}

/** Try to get a free node - giving only 'timeout' milliseconds
    to get that node. This returns a null node if the call
    is unsuccessful */
Node Nodes::getNode(int timeout)
{
    QTime t;
    t.start();
    
    if (timeout <= 0)
        return this->getNode();

    if (this->isEmpty())
        return Node();
    
    QMutexLocker lkr( &(d->datamutex) );
    
    shared_ptr<QSemaphore> nodesem = d->nodesem;
    
    bool reserved_node = false;
    
    while (not reserved_node)
    {
        lkr.unlock();
        
        if (nodesem.get() == 0)
            //all of the nodes have gone!
            return Node();
            
        int new_timeout = timeout - t.elapsed();
        if (new_timeout <= 0)
            //we've run out of time
            return Node();
        
        #if QT_VERSION >= 0x040300
            //try to reserve a node
            if (not nodesem->tryAcquire(1, new_timeout))
                //we ran out of time
                return Node();
        #else
            if (not nodesem->tryAcquire(1))
            {
                bool acquired_node = false;
            
                while (new_timeout > 0)
                {
                    if (nodesem->tryAcquire(1))
                    {
                        acquired_node = true;
                        break;
                    }
                    
                    --new_timeout;
                    sleep(1);
                }
                
                if (not acquired_node)
                    //we ran out of time
                    return Node();
            }
        #endif
            
        lkr.relock();
        
        if (t.elapsed() >= timeout)
        {
            //we've run out of time
            nodesem->release();
            return Node();
        }
        
        if (nodesem.get() != d->nodesem.get())
        {
            //the nodes changed while we were waiting 
            //our reservation isn't valid
            nodesem->release();
            nodesem = d->nodesem;
        }
        else
            //we've got a valid reservation!
            reserved_node = true;
    }

    //collect the reservation
    return this->_pvt_getNode();
}

/** Try to get 'n' free nodes, within the time 'timeout'.
    If this fails, then *no* nodes are returned */
QList<Node> Nodes::getNodes(int n, int timeout)
{
    QTime t;
    t.start();

    if (timeout <= 0)
        return this->getNodes(n);

    QList<Node> nodes;
    
    if (n <= 0)
    {
        return nodes;
    }
    else if (n == 1)
    {
        nodes.append( this->getNode(timeout) );
    }
    else
    {
        QMutexLocker lkr( &(d->datamutex) );

        shared_ptr<QSemaphore> nodesem = d->nodesem;
        
        bool reserved_nodes = false;
        
        while (not reserved_nodes)
        {
            lkr.unlock();
            
            if (nodesem.get() == 0)
                //there are no nodes
                return nodes;
                
            int new_timeout = timeout - t.elapsed();
            if (new_timeout <= 0)
                //we've run out of time
                return nodes;
                
            //reserve the nodes
            if (not nodesem->tryAcquire(n))
                //we ran out of time
                return nodes;
            
            lkr.relock();
                
            if (t.elapsed() >= timeout)
            {
                //we've run out of time
                nodesem->release(n);
                return nodes;
            }
            
            if (nodesem.get() != d->nodesem.get())
            {
                //the nodes changed, so our reservation isn't valid
                nodesem->release(n);
                nodesem = d->nodesem;
            }
            else
                //we've got a reservation!
                reserved_nodes = true;
        }

        //grab onto our reserved nodes
        for (int i=0; i<n; ++i)
        {
            nodes.append( _pvt_getNode() );
        }
    }
    
    return nodes;
}

/** Try to get all the nodes, within the time 'timeout'.
    If this fails, then *no* nodes are returned */
QList<Node> Nodes::getAllNodes(int timeout)
{
    QTime t;
    t.start();

    if (timeout <= 0)
        return this->getAllNodes();

    if (this->isEmpty())
        return QList<Node>();

    QMutexLocker lkr( &(d->datamutex) );

    shared_ptr<QSemaphore> nodesem = d->nodesem;
    
    bool reserved_nodes = false;
    
    while (not reserved_nodes)
    {
        lkr.unlock();
        
        if (nodesem.get() == 0)
            //there are no nodes
            return QList<Node>();
            
        int new_timeout = timeout - t.elapsed();
        if (new_timeout <= 0)
            //we've run out of time
            return QList<Node>();
            
        //reserve the nodes
        if (not nodesem->tryAcquire( nodesem->available()) )
            //we ran out of time
            return QList<Node>();
        
        lkr.relock();
            
        if (t.elapsed() >= timeout)
        {
            //we've run out of time
            nodesem->release( nodesem->available() );
            return QList<Node>();
        }
        
        if (nodesem.get() != d->nodesem.get())
        {
            //the nodes changed, so our reservation isn't valid
            nodesem->release( nodesem->available() );
            nodesem = d->nodesem;
        }
        else
            //we've got a reservation!
            reserved_nodes = true;
    }

    //grab onto our reserved nodes
    QList<Node> nodes;

    int n = nodesem->available();

    for (int i=0; i<n; ++i)
    {
        nodes.append( _pvt_getNode() );
    }
    
    return nodes;
}

/** Wait until all of the nodes are free */
void Nodes::waitUntilAllFree()
{
    if (d.get() == 0)
        return;

    QMutexLocker lkr( &(d->datamutex) );

    while (not d->busy_frontends.isEmpty())
    {
        d->waiter.wait( &(d->datamutex) );
    }
}

/** Wait until all of the nodes are free, or until
    timeout milliseconds have passed - this returns
    whether or not all of the nodes are free */
bool Nodes::waitUntilAllFree(int timeout)
{
    QTime t;
    t.start();

    if (timeout < 0)
    {
        this->waitUntilAllFree();
        return true;
    }

    if (d.get() == 0)
        return true;
        
    QMutexLocker lkr( &(d->datamutex) );
    
    while (not d->busy_frontends.isEmpty())
    {
        int new_timeout = timeout - t.elapsed();
        if (new_timeout <= 0)
            return false;
    
        if (not d->waiter.wait( &(d->datamutex), new_timeout) )
            return false;
    }
    
    return true;
}

/** Return the number of free nodes */
int Nodes::nFree()
{
    QMutexLocker lkr( &(d->datamutex) );
    return d->free_frontends.count();
}

/** Return the number of busy nodes */
int Nodes::nBusy()
{
    QMutexLocker lkr( &(d->datamutex) );
    return d->busy_frontends.count();
}

/** Return the total number of nodes available */
int Nodes::nNodes()
{
    if (d.get() == 0)
    {
        return 0;
    }
    else
    {
        QMutexLocker lkr( &(d->datamutex) );
        return d->frontends.count();
    }
}

/** Return the total number of nodes available */
int Nodes::count()
{
    return this->nNodes();
}

/** Add the node 'node' to this set. This will be added as a busy
    Node (as obviously someone is holding onto it!) This does
    nothing if this node is already part of this set */
void Nodes::add(Node node)
{
    if (node.isNull())
        return;
        
    if (d.get() == 0)
    {
        d.reset( new NodesPvt() );
    }
    
    QUuid uid = node.UID();

    if (uid.isNull())
        //something is dodgy with this node...
        return;
    
    QMutexLocker lkr( &(d->datamutex) );
    
    Nodes old_nodes = node.nodes();
    
    //are we 'old_nodes'?
    if (old_nodes.d.get() == d.get())
    {
        //yes we are! - but are we sure?
        if (not d->frontends.contains(uid))
            throw SireError::program_bug( QObject::tr(
                "The node %1 thinks it is part of this Nodes object, but "
                "this Nodes object isn't so sure!").arg(uid.toString()),
                    CODELOC );
    
        //nothing to do
        return;
    }
    
    //remove the node from its old set
    old_nodes.remove(node);
    
    //add it to this set
    d->frontends.insert(uid, node.frontend());
    d->busy_frontends.insert(uid);
    
    //increase the semaphore count
    shared_ptr<QSemaphore> old_nodesem = d->nodesem;

    d->nodesem.reset( new QSemaphore(d->frontends.count()) );
    
    if (old_nodesem.get() != 0)
    {
        old_nodesem->release( old_nodesem->available() );
    }
    
    //tell the node that it is now living here
    node.rehome(*this);
}

/** Add all of the nodes in 'nodes' to this set */
void Nodes::add(Nodes &nodes)
{
    while (nodes.nNodes() > 0)
    {
        Node node = nodes.getNode(250);
        
        if (not node.isNull())
            this->add(node);
    }
}

/** Try to add one more node to this set by taking a node from the 
    pool - this only looks for immediately available nodes, and may
    not work! */
void Nodes::addNode()
{
    Nodes newnode = Cluster::getNode();
    this->add(newnode);
}

/** Try to add one more node to this set by taking a node from the
    pool - this only tries to find an available node for 
    'timeout' milliseconds, and so it may fail */
void Nodes::addNode(int timeout)
{
    Nodes newnode = Cluster::getNode(timeout);
    this->add(newnode);
}

/** Try to add up to 'n' nodes to this set, by taking the nodes
    from the pool - this only looks for immediately available
    nodes, so you may get less than 'n' (you may even get zero!) */
void Nodes::addNodes(int n)
{
    Nodes newnodes = Cluster::getNodes(n);
    this->add(newnodes);
}

/** Try to add up to 'n' nodes to this set, by taking the nodes
    from the pool - this only looks for available
    nodes for 'timeout' milliseconds, so you may get less 
    than 'n' (you may even get zero!) */
void Nodes::addNodes(int n, int timeout)
{
    Nodes newnodes = Cluster::getNodes(n, timeout);
    this->add(newnodes);
}

/** Remove the node 'node' from this set. This doesn't abort
    the job running on the node, but the node will be automatically
    returned to the Cluster pool once it is destroyed.
      
    This does nothing if this node isn't in this set.
*/
void Nodes::remove(Node node)
{
    if (d.get() == 0)
        return;
        
    QUuid uid = node.UID();
        
    QMutexLocker lkr( &(d->datamutex) );
    
    if (d->frontends.contains(uid))
    {
        //remove this frontend
        d->frontends.remove(uid);
        d->busy_frontends.remove(uid);
        d->free_frontends.removeAll(uid);
        
        //reset the semaphore
        if (d->frontends.isEmpty())
        {
            d->nodesem.reset();
        }
        else
            d->nodesem.reset( new QSemaphore(d->frontends.count()) );
    
        //tell the node!
        node.evict();
    }
}

/** Remove all nodes from this set - this aborts any running
    jobs, then disconnects from all of the backends. 
    Note that this does not block  */
void Nodes::removeAll()
{
    if (d.get() == 0)
        return;

    QMutexLocker lkr( &(d->datamutex) );

    //remove the current semaphore - this stop anyone asking for more nodes
    shared_ptr<QSemaphore> old_nodesem = d->nodesem;
    d->nodesem.reset();

    QList<Frontend> frontends = d->frontends.values();
    
    d->frontends.clear();
    d->busy_frontends.clear();
    d->free_frontends.clear();
    
    //wake anyone waiting for all the nodes to become free
    d->waiter.wakeAll();
    
    lkr.unlock();

    //wake up everyone waiting for a node, as they are now
    //out of luck
    if (old_nodesem.get() != 0)
        old_nodesem->release( old_nodesem->available() );
    
    //finally, tell each frontend to abort it current job,
    //so that the node will exit and the backend will be
    //returned to the pool
    foreach (Frontend frontend, frontends)
    {
        frontend.abortJob();
    }
}



//////////////
////////////// Implementation of NodesPtr
//////////////

/** Construct a null pointer */
NodesPtr::NodesPtr()
{}

/** Construct to point to 'nodes' */
NodesPtr::NodesPtr(const Nodes &nodes) : d( nodes.d )
{}

/** Copy constructor */
NodesPtr::NodesPtr(const NodesPtr &other) : d( other.d )
{}

/** Destructor */
NodesPtr::~NodesPtr()
{}

/** Copy assignment operator */
NodesPtr& NodesPtr::operator=(const NodesPtr &other)
{
    d = other.d;
    return *this;
}

/** Return the Nodes object - an empty set of Nodes will
    be returned if this is a null pointer */
Nodes NodesPtr::operator*() const
{
    return Nodes( d.lock() );
}

/** Return the Nodes object - an empty set of Nodes will
    be returned if this is a null pointer */
Nodes NodesPtr::lock() const
{
    return Nodes( d.lock() );
}

/** Return whether or not this pointer is null */
bool NodesPtr::expired() const
{
    return d.expired();
}

/** Reset this pointer to null */
void NodesPtr::reset()
{
    d.reset();
}

/** Function called by the destructor of NodePvt to return
    a front end back to the Nodes from whence it came */
void NodesPtr::returnFrontend(Frontend frontend)
{
    shared_ptr<NodesPvt> nodes = d.lock();
    
    if (nodes.get() == 0)
        //we don't have a parent Nodes object any more
        return;
    
    QUuid uid = frontend.UID();
                
    QMutexLocker lkr( &(nodes->datamutex) );
    
    if (nodes->busy_frontends.contains(uid))
    {
        //the parent Nodes object still contains this Node
    
        //put the local front ends at the beginning of 
        //the list so that they are used first

        if (frontend.isLocal())
            nodes->free_frontends.prepend(uid);
        else
            nodes->free_frontends.append(uid);
    
        nodes->busy_frontends.remove(uid);
        
        if (nodes->busy_frontends.isEmpty())
        {
            //wake up anyone waiting for all of the nodes to become free
            nodes->waiter.wakeAll();
        }
        
        if (nodes->nodesem.get() != 0)
            nodes->nodesem->release();
    }
}

static QThreadStorage<QUuid*> this_thread_uids;

/** Let this Nodes scheduler borrow this thread to run WorkPackets.
    Technically, this doesn't use the current thread, but instead 
    creates a duplicate, so you are still able to use your thread.
    However, you should avoid doing anything compute intensive in
    your thread while the Node has borrowed it, as otherwise you 
    risk using more CPU than is available.
    
    This returns a ThisThread holder that is used to ask the Nodes
    object to return this thread when the ThisThread object is
    deleted, or when ThisThread::reclaim() is called.
*/
ThisThread Nodes::borrowThisThread()
{
    if (this_thread_uids.hasLocalData())
        //this thread is already being borrowed
        return ThisThread();
        
    else
    {
        if (d.get() == 0)
            d.reset( new NodesPvt() );

        return ThisThread(*this);
    }
}

/** Create a temporary backend for this thread, and return
    the UID of that backend */
QUuid Nodes::createThisThread()
{
    BOOST_ASSERT( d.get() != 0 );

    Frontend frontend( Backend::createLocalOnly() );

    QMutexLocker lkr( &(d->datamutex) );
    
    //unlock anyone waiting for a node
    if (d->nodesem.get() != 0)
        d->nodesem->release( d->nodesem->available() );
    
    //add the new frontend to the set of available frontends
    QUuid uid = frontend.UID();
    
    d->frontends.insert( uid, frontend );
    d->free_frontends.prepend(uid);
    
    d->nodesem.reset( new QSemaphore(d->frontends.count()) );
    
    return uid;
}

/** Delete the temporary backend for this thread, which
    was created with the UID 'uid' */
void Nodes::reclaimThisThread(const QUuid &uid)
{
    BOOST_ASSERT( d.get() != 0);

    //ensure that the local thread has UID 'uid'
    if (this_thread_uids.hasLocalData())
    {
        if (uid != *(this_thread_uids.localData()))
        {
            throw SireError::program_bug( QObject::tr(
                "We cannot reclaim the local thread (UID == %1) as "
                "the UID of the local thread should be %2.") 
                    .arg(uid.toString(), this_thread_uids.localData()->toString()),
                        CODELOC );
        }
        
        //yes - it is safe to remove this thread
        QMutexLocker lkr( &(d->datamutex) );
        
        if (d->nodesem.get() != 0)
            d->nodesem->release( d->nodesem->available() );
        
        d->frontends.remove(uid);
        d->free_frontends.removeAll(uid);
        d->busy_frontends.remove(uid);

        if (d->frontends.isEmpty())
            d->nodesem.reset();
            
        else
            d->nodesem.reset( new QSemaphore(d->frontends.count()) );
    }
}

//////////////
////////////// Implementation of ThisThread
//////////////

namespace SireCluster
{
namespace detail
{

class ThisThreadPvt
{
public:
    ThisThreadPvt(Nodes nodes)
    {
        if (not this_thread_uids.hasLocalData())
        { 
            QUuid uid = nodes.createThisThread();
            this_thread_uids.setLocalData( new QUuid(uid) );
            nodesptr = nodes;
        }
    }
    
    ~ThisThreadPvt()
    {
        Nodes nodes = nodesptr.lock();

        if (not nodes.isEmpty())
        {
            QUuid *uid = this_thread_uids.localData();
    
            //we need to remove this thread from the nodes object
            if (uid)
                nodes.reclaimThisThread(*uid);
        
            this_thread_uids.setLocalData(0);
        }
    }
    
    NodesPtr nodesptr;
};

} // end of namespace detail
} // end of namespace SireCluster

/** Construct a null locker */
ThisThread::ThisThread()
{}

/** Construct a lock for this thread being added to 
    the nodes object 'nodes' */
ThisThread::ThisThread(const Nodes &nodes) : d( new ThisThreadPvt(nodes) )
{}

/** Copy constructor */
ThisThread::ThisThread(const ThisThread &other) : d(other.d)
{}

/** Destructor */
ThisThread::~ThisThread()
{}

/** Copy assignment operator */
ThisThread& ThisThread::operator=(const ThisThread &other)
{
    d = other.d;
    return *this;
}

/** Reclaim this thread */
void ThisThread::reclaim()
{
    d.reset();
}
