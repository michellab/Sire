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
#include <mpi.h>             // CONDITIONAL_INCLUDE
#include "mpi/mpicluster.h"  // CONDITIONAL_INCLUDE
#endif

#include <QHash>
#include <QMutex>
#include <QWaitCondition>

#include <QTime>

#include "cluster.h"
#include "backend.h"
#include "frontend.h"
#include "nodes.h"

#include "SireError/errors.h"
#include "SireError/printerror.h"

#ifdef Q_OS_WIN
   #include <windows.h>    // CONDITIONAL_INCLUDE
#else
   #ifdef Q_OS_UNIX
       #include <unistd.h>  // CONDITIONAL_INCLUDE
   #endif
#endif


#include <QDebug>

using namespace SireCluster;

using boost::shared_ptr;

/** Private implementation of Cluster */
class ClusterPvt
{
public:
    ClusterPvt();
    
    ~ClusterPvt();
    
    /** Mutex to protect access to the registry */
    QMutex datamutex;
    
    /** Waitcondition used to wait until the cluster has been
        shut down - this allows the 'exec' function to work */
    QWaitCondition execwaiter;
    
    /** The registry of all backends that are local to this 
        address space */
    QHash<QUuid,Backend> local_backends;
};

using namespace SireCluster::detail;

/** Return whether or not there is any point using MPI
    (only worth it if we have more than one MPI process!) */
static bool usingMPI()
{
    #ifdef SIRE_USE_MPI
        if (not ::MPI::Is_initialized())
        {
            int argc = 0;
            char **argv = 0;
            ::MPI::Init(argc,argv);
        }
        
        return ::MPI::COMM_WORLD.Get_size() > 1;
    #else
        return false;
    #endif
}

Q_GLOBAL_STATIC( QMutex, globalMutex );
Q_GLOBAL_STATIC( QMutex, execMutex );

static ClusterPvt *global_cluster(0);
static bool global_cluster_is_running = false;

ClusterPvt* globalCluster()
{
    QMutexLocker lkr( globalMutex() );
    
    if (global_cluster == 0)
    {
        global_cluster = new ClusterPvt();
        
        lkr.unlock();
        
        execMutex()->lock();
        global_cluster_is_running = true;
        execMutex()->unlock();
        
        #ifdef SIRE_USE_MPI
            if (::usingMPI())
                SireCluster::MPI::MPICluster::start();
        #endif
        
        //create a new backend
        Backend::create();
        
        #ifdef SIRE_USE_MPI
            //now wait until we have all got here
            if (::usingMPI())
                SireCluster::MPI::MPICluster::sync();
        #endif
    }
    
    return global_cluster;
}

/** Constructor for the global cluster */
ClusterPvt::ClusterPvt()
{}

/** Destructor */
ClusterPvt::~ClusterPvt()
{
    execwaiter.wakeAll();
}

/** Start the cluster - this is like exec, but it doesn't
    block until the cluster has been shutdown */
void Cluster::start(int ppn)
{
    globalCluster();
    
    if (ppn > 1)
    {
        //add extra threads to this process
        for (int i=1; i<ppn; ++i)
        {
            Backend::create();
        }
    }
}

/** Return whether or not the cluster is running */
bool Cluster::isRunning()
{
    QMutexLocker lkr( execMutex() );
    return global_cluster_is_running;
}

/** Wait for the global cluster to stop running */
void Cluster::wait()
{
    QMutexLocker lkr( execMutex() );
    
    if (global_cluster_is_running)
        globalCluster()->execwaiter.wait( execMutex() );
}

/** Return whether or not this cluster supports MPI
    (this is true if MPI is available, and there is more than
     one MPI process) */
bool Cluster::supportsMPI()
{
    return ::usingMPI();
}

/** Return the rank of the process - this is either the 
    rank of this process in the MPI group, or it is 0 */
int Cluster::getRank()
{
    #ifdef SIRE_USE_MPI
        if (::usingMPI())
            return SireCluster::MPI::MPICluster::getRank();
        else
            return 0;
    #else
        return 0;
    #endif
}

/** Return the number of processes - this is either the 
    size of the MPI group, or it is 1 */
int Cluster::getCount()
{
    #ifdef SIRE_USE_MPI
        if (::usingMPI())
            return SireCluster::MPI::MPICluster::getCount();
        else
            return 1;
    #else
        return 1;
    #endif
}

/** Register the backend 'backend' with the cluster */
void Cluster::registerBackend(const Backend &backend)
{
    if (backend.isNull())
        return;

    QMutexLocker lkr( &(globalCluster()->datamutex) );

    if (not globalCluster()->local_backends.contains(backend.UID()))
    {
        globalCluster()->local_backends.insert( backend.UID(), backend );
        
        #ifdef SIRE_USE_MPI
            //now inform the MPI connected nodes that a backend with this
            //UID is available on this node
            if (::usingMPI())
            {
                SireCluster::MPI::MPICluster::registerBackend(backend);
            }
        #endif
    }
}

/** Return all of the available local backends, connected
    to local frontends - this is used by ReservationManager to
    reserve local frontends */
QList<Frontend> Cluster::localBackends()
{
    QList<QUuid> local_uids = Cluster::localUIDs();
    
    QList<Frontend> frontends;
    
    foreach (QUuid uid, local_uids)
    {
        Frontend frontend = Cluster::_pvt_getFrontend(uid);
        
        if (not frontend.isNull())
            frontends.append(frontend);
    }
    
    return frontends;
}

/** Return a Frontend that will allow us to communicate with
    an available backend (this gets the first available backend,
    though does prefer to find a local backend if possible) 
    
    A null frontend is returned no resources are currently available.
*/
Frontend Cluster::_pvt_getFrontend()
{
    QMutexLocker lkr( &(globalCluster()->datamutex ) );
    
    //first loop through the local backends to see if any
    //of them are available
    for (QHash<QUuid,Backend>::const_iterator 
                    it = globalCluster()->local_backends.constBegin();
         it != globalCluster()->local_backends.constEnd();
         ++it)
    {
        Frontend frontend = Frontend::tryAcquire( it.value() );

        if (not frontend.isNull())
            //we successfully grabbed this backend
            return frontend;
    }
    
    //ok, there were no locally available frontends
    lkr.unlock();

    #ifdef SIRE_USE_MPI
         if (::usingMPI())
            //see if there are any available remote backends
            return SireCluster::MPI::MPICluster::getFrontend();
    #endif

    //nothing could be found - return a null frontend
    return Frontend();
}

/** Return up to 'n' Frontends that will allow us to communicate with
    up to 'n' available backends (this gets the first available backends,
    though does prefer to find a local backend if possible) 
    
    An empty list is returned no resources are currently available.
*/
QList<Frontend> Cluster::_pvt_getFrontends(int n)
{
    QMutexLocker lkr( &(globalCluster()->datamutex ) );
    
    QList<Frontend> frontends;
    
    //first loop through the local backends to see if any
    //of them are available
    for (QHash<QUuid,Backend>::const_iterator 
                    it = globalCluster()->local_backends.constBegin();
         it != globalCluster()->local_backends.constEnd();
         ++it)
    {
        Frontend frontend = Frontend::tryAcquire( it.value() );
        
        if (not frontend.isNull())
        {
            //we successfully grabbed this backend
            frontends.append(frontend);
            --n;
            if (n == 0)
                return frontends;
        }
    }
    
    //ok, there were no locally available frontends
    lkr.unlock();

    #ifdef SIRE_USE_MPI
         if (::usingMPI())
         {
            //see if there are any available remote backends
            QList<Frontend> mpifrontends = SireCluster::MPI::MPICluster::getFrontends(n);
            
            if (not mpifrontends.isEmpty())
                frontends += mpifrontends;
         }
    #endif

    //return all the frontends that have been found
    return frontends;
}

/** Return a Frontend that will allow us to communicate with the 
    Backend with UID 'uid'. A null frontend will be returned
    if this backend is busy. An error will be raised if there
    is no backend associated with this UID
    
    \throw SireError::unavailable_resource
*/
Frontend Cluster::_pvt_getFrontend(const QUuid &uid)
{
    if (uid.isNull())
        throw SireError::unavailable_resource( QObject::tr(
            "There is no front end for the null backend!"), CODELOC );

    QMutexLocker lkr( &(globalCluster()->datamutex) );
    
    if (globalCluster()->local_backends.contains(uid))
    {
        //return a local frontend for this local backend
        return Frontend::tryAcquire( globalCluster()->local_backends.value(uid) );
    }
    else
    {
        lkr.unlock();
        
        #ifdef SIRE_USE_MPI
            if (::usingMPI())
                //see if this node exists on any of the MPI nodes...
                return SireCluster::MPI::MPICluster::getFrontend(uid);
            else
                return Frontend();
        #else
            return Frontend();
        #endif
    }
}

/** Return a Frontend that will allow us to communicate with
    an available backend (this gets the first available backend,
    though does prefer to find a local backend if possible) 
    
    This will keep trying to get a frontend for up to 'timeout'
    milliseconds. Use a negative timeout to wait forever
    (well, for-ages - which is until it is available or until
    the cluster is shut down) 
    
    A null frontend is returned no resources are currently available
    within the specified timeout
*/
Frontend Cluster::getFrontend(int timeout)
{
    //ensure that the global cluster has been created
    globalCluster();

    if (Cluster::UIDs().isEmpty())
        //there are no frontends in the cluster!
        return Frontend();

    if (timeout < 0)
    {
        while (true)
        {
            if (not Cluster::isRunning())
                //oh dear - the cluster has stopped
                return Frontend();
        
            Frontend frontend = Cluster::_pvt_getFrontend();
            
            if (not frontend.isNull())
                //we've found a frontend
                return frontend;

            //wait a second
            #ifdef Q_OS_WIN
                Sleep(1);
            #else
                sleep(1);
            #endif
        }
    }
    else
    {
        QTime t;
        t.start();
        
        while (t.elapsed() < timeout)
        {
            if (not Cluster::isRunning())
                //oh dear - the cluster has stopped
                return Frontend();

            Frontend frontend = Cluster::_pvt_getFrontend();
            
            if (not frontend.isNull())
                //we've found a frontend
                return frontend;
        
            //only try once a second
            if (t.elapsed() + 1000 > timeout)
                return Frontend();;
                
            #ifdef Q_OS_WIN
                Sleep(1);
            #else
                sleep(1);
            #endif
        }
    }
    
    return Frontend();
}

/** Return up to 'n' Frontends that will allow us to communicate with
    up to 'n' available backends (this gets the first available backend,
    though does prefer to find a local backend if possible) 
    
    This will keep trying to get the frontends for up to 'timeout'
    milliseconds. Use a negative timeout to wait forever
    (well, for-ages - which is until it is available or until
    the cluster is shut down) 
    
    An empty list is returned if no resources are currently available
    within the specified timeout
*/
QList<Frontend> Cluster::getFrontends(int n, int timeout)
{
    //ensure that the global cluster has been created
    globalCluster();

    QList<Frontend> frontends;

    n = qMin( n, Cluster::UIDs().count() );

    if (timeout < 0)
    {
        while (frontends.count() < n)
        {
            if (not Cluster::isRunning())
                //oh dear - the cluster has stopped
                return frontends;
        
            int nremaining = n - frontends.count();
        
            QList<Frontend> active_frontends = Cluster::_pvt_getFrontends(nremaining);
            
            if (not active_frontends.isEmpty())
                //we've found a frontend
                frontends += active_frontends;

            if (frontends.count() < n)
            {
                //wait a second
                #ifdef Q_OS_WIN
                    Sleep(1);
                #else
                    sleep(1);
                #endif
            }
        }
    }
    else
    {
        QTime t;
        t.start();
        
        while (frontends.count() < n and t.elapsed() < timeout)
        {
            if (not Cluster::isRunning())
                //oh dear - the cluster has stopped
                return frontends;

            int nremaining = n - frontends.count();

            QList<Frontend> active_frontends = Cluster::_pvt_getFrontends(nremaining);
            
            if (not active_frontends.isEmpty())
                //we've found a frontend
                return frontends += active_frontends;
        
            //only try once a second
            if (t.elapsed() + 1000 > timeout)
                return frontends;
                
            if (frontends.count() < n)
            {
                #ifdef Q_OS_WIN
                    Sleep(1);
                #else
                    sleep(1);
                #endif
            }
        }
    }
    
    return frontends;
}

/** Return a Frontend that will allow us to communicate with the 
    Backend with UID 'uid'. A null frontend will be returned
    if this backend is busy. An error will be raised if there
    is no backend associated with this UID
    
    This will keep trying to get the frontend for up to 'timeout'
    milliseconds. Use a negative timeout to wait forever
    (well, for-ages - which is until it is available or until
    the cluster is shut down) 
    
    \throw SireError::unavailable_resource
*/
Frontend Cluster::getFrontend(const QUuid &uid, int timeout)
{
    //ensure that the global cluster has been created
    globalCluster();

    if (timeout < 0)
    {
        while (true)
        {
            if (not Cluster::isRunning())
                //oh dear - the cluster has stopped
                return Frontend();
        
            Frontend frontend = Cluster::_pvt_getFrontend(uid);
            
            if (not frontend.isNull())
                //we've found a frontend
                return frontend;
                
            //wait a second
            #ifdef Q_OS_WIN
                Sleep(1);
            #else
                sleep(1);
            #endif
        }
    }
    else
    {
        QTime t;
        t.start();
        
        while (t.elapsed() < timeout)
        {
            if (not Cluster::isRunning())
                //oh dear - the cluster has stopped
                return Frontend();

            Frontend frontend = Cluster::_pvt_getFrontend(uid);
            
            if (not frontend.isNull())
                //we've found a frontend
                return frontend;
        
            //only try once a second
            if (t.elapsed() + 1000 > timeout)
                return Frontend();;
                
            #ifdef Q_OS_WIN
                Sleep(1);
            #else
                sleep(1);
            #endif
        }
    }
    
    return Frontend();
}

/** Return a Frontend that will allow us to communicate with
    an available backend (this gets the first available backend,
    though does prefer to find a local backend if possible) 
    
    A null frontend is returned no resources are currently available
*/
Frontend Cluster::getFrontend()
{
    //ensure that the global cluster has been created
    globalCluster();

    if (not Cluster::isRunning())
        //oh dear - the cluster has stopped
        return Frontend();
        
    return Cluster::_pvt_getFrontend();
}

/** Return up to 'n' Frontends that will allow us to communicate with
    up to 'n' available backends (this gets the first available backend,
    though does prefer to find a local backend if possible) 
    
    An empty list is returned no resources are currently available
*/
QList<Frontend> Cluster::getFrontends(int n)
{
    //ensure that the global cluster has been created
    globalCluster();

    if (not Cluster::isRunning())
        //oh dear - the cluster has stopped
        return QList<Frontend>();
        
    return Cluster::_pvt_getFrontends(n);
}

/** Return a Frontend that will allow us to communicate with the 
    Backend with UID 'uid'. A null frontend will be returned
    if this backend is busy. An error will be raised if there
    is no backend associated with this UID
    
    \throw SireError::unavailable_resource
*/
Frontend Cluster::getFrontend(const QUuid &uid)
{
    //ensure that the global cluster has been created
    globalCluster();

    if (not Cluster::isRunning())
        //oh dear - the cluster has stopped
        return Frontend();
        
    return Cluster::_pvt_getFrontend(uid);
}

/** Return a Nodes object that contains just a single node.
    This doesn't block - it just grabs the first available 
    node, and if there are none available then it returns
    immediately, returning an empty Nodes object */
Nodes Cluster::getNode()
{
    //get the first available frontend
    Frontend frontend = Cluster::getFrontend();
    
    if (frontend.isNull())
        return Nodes();
        
    else
        return Nodes(frontend);
}

/** Return a Nodes object that contains just a single node. 
    This blocks until a Node is available, or until 'timeout'
    milliseconds has passed (use negative timeout to wait forever).
    
    There are some cases where a Node is just not available, in which
    case an empty Nodes object will be returned
*/
Nodes Cluster::getNode(int timeout)
{
    //get the first available frontend
    Frontend frontend = Cluster::getFrontend(timeout);
    
    if (frontend.isNull())
        return Nodes();
        
    else
        return Nodes(frontend);
}

/** Return a Nodes object that contains the node with UID 'uid'.
    
    There are some cases where a Node is just not available, in which
    case an empty Nodes object will be returned.
    
    \throw SireError::unavailable_resource
*/
Nodes Cluster::getNode(const QUuid &uid)
{
    Frontend frontend = Cluster::getFrontend(uid);
    
    if (frontend.isNull())
        return Nodes();
        
    else
        return Nodes(frontend);
}

/** Return a Nodes object that contains the node with UID 'uid'.
    This blocks until the Node is available, or until 'timeout'
    milliseconds has passed (use negative timeout to wait forever).
    
    There are some cases where a Node is just not available, in which
    case an empty Nodes object will be returned.
    
    \throw SireError::unavailable_resource
*/
Nodes Cluster::getNode(const QUuid &uid, int timeout)
{
    Frontend frontend = Cluster::getFrontend(uid, timeout);
    
    if (frontend.isNull())
        return Nodes();
        
    else
        return Nodes(frontend);
}

/** Return a Nodes object containing up to 'nnodes' nodes. This function
    will do its best, but you may end with less than you asked for
    (or even none at all!). */
Nodes Cluster::getNodes(int nnodes)
{
    QList<Frontend> frontends = Cluster::getFrontends(nnodes);
    
    if (frontends.isEmpty())
        return Nodes();
        
    else
        return Nodes(frontends);
}

/** Return a Nodes object containing up to 'nnodes' nodes. This function
    will do its best, but you may end with less than you asked for
    (or even none at all!). It is a bad idea to ask for more nodes
    than there are backends using an infinite timeout... */
Nodes Cluster::getNodes(int nnodes, int timeout)
{
    QList<Frontend> frontends;
    
    if (timeout < 0)
    {
        //keep going for-ages
        while (frontends.count() < nnodes)
        {
            int nremaining = nnodes - frontends.count();
        
            QList<Frontend> active_frontends = Cluster::getFrontends(nremaining, timeout);

            if (not active_frontends.isEmpty())
                frontends += active_frontends;
        }
    }
    else
    {
        QTime t;
        t.start();
        
        int remaining_time = timeout - t.elapsed();
        
        while (frontends.count() < nnodes and remaining_time > 0)
        {
            int nremaining = nnodes - frontends.count();
        
            QList<Frontend> active_frontends = Cluster::getFrontends(nremaining, 
                                                                     remaining_time);
            
            if (not active_frontends.isEmpty())
                frontends += active_frontends;
            
            remaining_time = timeout - t.elapsed();
        }
    }
    
    if (frontends.isEmpty())
        return Nodes();
        
    else
        return Nodes(frontends);
}

/** Return a Nodes object that contains as many of the nodes with 
    UIDs from 'uids' as possible, within the time allowed. Note that
    this may not give you all of the nodes (it may give you none!).
        
    \throw SireError::unavailable_resource
*/
Nodes Cluster::getNodes(const QList<QUuid> &uids, int timeout)
{
    QHash<QUuid,Frontend> frontends;
    
    if (timeout < 0)
    {
        foreach( QUuid uid, uids )
        {
            if (frontends.contains(uid))
                continue;
                
            Frontend frontend = Cluster::getFrontend(uid, timeout);
            
            if (not frontend.isNull())
                frontends.insert(uid, frontend);
        }
    }
    else
    {
        QTime t;
        t.start();
    
        bool got_them_all = false;
        
        int remaining_time = timeout - t.elapsed();
        
        while (remaining_time > 0 and (not got_them_all))
        {
            foreach ( QUuid uid, uids )
            {
                if (remaining_time <= 0)
                    break;

                if (frontends.contains(uid))
                    continue;
                    
                //only try to get the node for 1/8th of the time
                Frontend frontend = Cluster::getFrontend(uid, 1 + (remaining_time/8));
                
                if (not frontend.isNull())
                    frontends.insert(uid, frontend);
                    
                remaining_time = timeout - t.elapsed();
            }
        }
    }
    
    if (frontends.isEmpty())
        return Nodes();
        
    else
        return Nodes( frontends.values() );
}

/** Return a Nodes object that contains as many of the nodes with 
    UIDs from 'uids' as possible, within the time allowed. Note that
    this may not give you all of the nodes (it may give you none!).
        
    \throw SireError::unavailable_resource
*/
Nodes Cluster::getNodes(const QList<QUuid> &uids)
{
    QList<Frontend> frontends;

    foreach (QUuid uid, uids)
    {
        Frontend frontend = Cluster::getFrontend(uid);
        
        if (not frontend.isNull())
            frontends.append(frontend);
    }
    
    if (frontends.isEmpty())
        return Nodes();
        
    else
        return Nodes( frontends );
}

/** Try to get hold of all of the nodes that are available on this
    cluster */
Nodes Cluster::getAllNodes()
{
    //how many nodes are available?
    int nnodes = Cluster::UIDs().count();
    
    return Cluster::getNodes(nnodes);
}

/** Try to get hold of all of the nodes that are available on this
    cluster (within the specified timeout) */
Nodes Cluster::getAllNodes(int timeout)
{
    //how many nodes are available?
    int nnodes = Cluster::UIDs().count();
    
    return Cluster::getNodes(nnodes, timeout);
}

/** Return the list of all of the UIDs of the local nodes
    (the nodes that exist in this address space) */
QList<QUuid> Cluster::localUIDs()
{
    QMutexLocker lkr( &(globalCluster()->datamutex) );

    return globalCluster()->local_backends.keys();
}

/** Return whether or not the backend with unique ID 'uid' 
    is local to this process */
bool Cluster::isLocal(const QUuid &uid)
{
    QMutexLocker lkr( &(globalCluster()->datamutex) );
    
    return globalCluster()->local_backends.contains(uid);
}

/** Return the list of all of the UIDs of all of the nodes 
    in this entire cluster */
QList<QUuid> Cluster::UIDs()
{
    #ifdef SIRE_USE_MPI
        if (::usingMPI())
            return SireCluster::MPI::MPICluster::UIDs();
        else
            return Cluster::localUIDs();
    #else
        return Cluster::localUIDs();
    #endif
}

/** Shutdown this cluster */
void Cluster::shutdown()
{
    /////check to see if we are still running
    {
        QMutexLocker lkr( execMutex() );
        
        if (not global_cluster_is_running)
            return;
            
        global_cluster_is_running = false;
    }

    #ifdef SIRE_USE_MPI
        if (::usingMPI())
            //shutdown MPI - this stops the backends from 
            //receiving any more work from remote nodes
            SireCluster::MPI::MPICluster::shutdown();
    #endif

    QMutexLocker lkr( &(globalCluster()->datamutex) );
    
    //shutdown all of the backends
    for (QHash<QUuid,Backend>::iterator it = globalCluster()->local_backends.begin();
         it != globalCluster()->local_backends.end();
         ++it)
    {
        it.value().shutdown();
    }

    //wake all threads waiting for the cluster to be shutdown
    globalCluster()->execwaiter.wakeAll();
}
