/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2016  Christopher Woods
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

#include "repexmove2.h"

#include "replica.h"
#include "replicas.h"

#include "suprasystem.h"

#include "suprasubsim.h"

#include "SireUnits/units.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireError/errors.h"

#include "tbb/task.h"

#include <QDebug>

using namespace SireMove;
using namespace SireMaths;
using namespace SireBase;
using namespace SireVol;
using namespace SireUnits;
using namespace SireUnits::Dimension;
using namespace SireStream;

////////////
//////////// Implementation of RepExMove
////////////

static const RegisterMetaType<RepExMove2> r_repexmove2;

/** Serialise to a binary datastream */
QDataStream SIREMOVE_EXPORT &operator<<(QDataStream &ds, const RepExMove2 &repexmove2)
{
    writeHeader(ds, r_repexmove2, 1);

    SharedDataStream sds(ds);
    
    sds << repexmove2.rangenerator
        << repexmove2.naccept
        << repexmove2.nreject
        << repexmove2.swap_monitors
        << repexmove2.disable_swaps
        << static_cast<const SupraMove&>(repexmove2);
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOVE_EXPORT &operator>>(QDataStream &ds, RepExMove2 &repexmove2)
{
    VersionID v = readHeader(ds, r_repexmove2);

    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> repexmove2.rangenerator
            >> repexmove2.naccept
            >> repexmove2.nreject
            >> repexmove2.swap_monitors
            >> repexmove2.disable_swaps
            >> static_cast<SupraMove&>(repexmove2);
    }
    else
        throw version_error(v, "1", r_repexmove2, CODELOC);
        
    return ds;
}

/** Constructor */
RepExMove2::RepExMove2()
           : ConcreteProperty<RepExMove2,SupraMove>(),
             naccept(0), nreject(0), swap_monitors(false), disable_swaps(false)
{}

/** Copy constructor */
RepExMove2::RepExMove2(const RepExMove2 &other)
           : ConcreteProperty<RepExMove2,SupraMove>(other),
             naccept(other.naccept), nreject(other.nreject),
             swap_monitors(other.swap_monitors),
             disable_swaps(other.disable_swaps)
{}

/** Destructor */
RepExMove2::~RepExMove2()
{}

/** Copy assignment operator */
RepExMove2& RepExMove2::operator=(const RepExMove2 &other)
{
    if (this != &other)
    {
        SupraMove::operator=(other);
        
        naccept = other.naccept;
        nreject = other.nreject;
        swap_monitors = other.swap_monitors;
        disable_swaps = other.disable_swaps;
    }
    
    return *this;
}

/** Comparison operator */
bool RepExMove2::operator==(const RepExMove2 &other) const
{
    return (this == &other) or
           (naccept == other.naccept and nreject == other.nreject and
            swap_monitors == other.swap_monitors and
            disable_swaps == other.disable_swaps and SupraMove::operator==(other));
}

/** Comparison operator */
bool RepExMove2::operator!=(const RepExMove2 &other) const
{
    return not this->operator==(other);
}

/** Return the total number of accepted replica exchange tests */
int RepExMove2::nAttempted() const
{
    return naccept + nreject;
}

/** Return the total number of accepted replica exchange tests */
int RepExMove2::nAccepted() const
{
    return naccept;
}

/** Return the total number of rejected replica exchange tests */
int RepExMove2::nRejected() const
{
    return nreject;
}

/** Return whether or not swap moves are disabled */
bool RepExMove2::swapMovesDisabled() const
{
    return disable_swaps;
}

/** Set disabling of swap moves */
void RepExMove2::setDisableSwaps(bool disable)
{
    disable_swaps = disable;
}

/** Return the average acceptance ratio of the replica exchange
    tests over all replicas */
double RepExMove2::acceptanceRatio() const
{
    if (this->nAttempted() > 0)
    {
        return double(this->nAccepted()) / double(this->nAttempted());
    }
    else
        return 0;
}

/** Return a string representation of this move */
QString RepExMove2::toString() const
{
    return QObject::tr("RepExMove2( %1 accepted, %2 rejected : %3 %% )")
                .arg(this->nAccepted())
                .arg(this->nRejected())
                .arg(100 * this->acceptanceRatio());
}

/** Clear the move statistics */
void RepExMove2::clearStatistics()
{
    naccept = 0;
    nreject = 0;
    SupraMove::clearStatistics();
}

/** Set the random number generator used for the replica exchange tests */
void RepExMove2::setGenerator(const RanGenerator &generator)
{
    rangenerator = generator;
}

/** Return the random number generator used for the replica exchange tests */
const RanGenerator& RepExMove2::generator() const
{
    return rangenerator;
}

/** Set whether or not to swap the system monitors when we swap the systems */
void RepExMove2::setSwapMonitors(bool swap)
{
    swap_monitors = swap;
}

/** Simple task that runs all of the moves of a single replica */
class RunReplicaTask : public tbb::task
{
    SupraSubSystemPtr *replica;
    bool record_stats;
    
public:
    RunReplicaTask(SupraSubSystemPtr *replica_, bool record_stats_)
        : replica(replica_), record_stats(record_stats_)
    {}
    
    tbb::task* execute()
    {
        if (replica == 0)
            return;

        //actually perform the replica moves, recording statistics
        //if record_stats is true
        result->edit().subMove(record_stats);
    }
};

/** Simple task that wraps up the running of all replicas
    in parallel as a single parent task. This spawns 
    repliacs.nReplicas() child tasks that run the
    simulations */
class RunReplicasTask : public tbb::task
{
    SupraSubSystemPtr *replicas;
    const int nreplicas;
    bool record_stats;

public:
    RunReplicasTask(QVector<SupraSubSystemPtr> &replicas_,
                    bool record_stats_)
     : replicas( replicas_.data() ), nreplicas(replicas_.count(), record_stats( record_stats_ )
    {}
    
    task* execute()
    {
        if (nreplicas <= 0 or replicas == 0)
            return 0;
        
        //set the reference count to nreplicas + 1 (add one as we will wait
        //for all replicas to finish)
        this->set_ref_count( nreplicas + 1 );
        
        //create one child task for each replica to be run
        tbb::task_list tasks;
        
        for (int i=0; i<nreplicas; ++i)
        {
            RunReplicaTask &task = *new( this->allocate_child() )
                                     RunReplicaTask(replicas[i], record_stats);
        
            tasks.push_back(task);
        }
        
        //spawn all child tasks and wait for them to complete
        this->spawn_and_wait_for_all(tasks);
    }
};

/** More complex task that runs replicas, with pairs of replicas run as children
    that perform the replica exchange test and mark whether or not the replicas
    should be swapped */
class RunReplicasWithSwapTask : public tbb::task
{
    SupraSubSystemPtr *replicas;
    const int nreplicas;
    bool record_stats;
    bool even_pairs;
    QVector< std::pair<int,int> > *to_swap;
    int *naccept;
    int *nreject;
    QMutex *swap_mutex;
    
public:
    RunReplicasWithSwapTask(QVector<SupraSubSystemPtr> &replicas_,
                            bool record_stats_, bool even_pairs_,
                            QVector< std::pair<int,int> > *to_swap_,
                            int *naccept_, int *nreject_, QMutex *swap_mutex_)
     : replicas( replicas_.data() ), nreplicas(replicas_.count(), record_stats( record_stats_ ),
       even_pairs(even_pairs_), to_swap(to_swap_), naccept(naccept_), nreject(nreject_),
       swap_mutex(swap_mutex_)
    {}
    
    tbb::task* execute()
    {
        if (nreplicas <= 0 or replicas == 0)
            return 0;
        
        //create a list of tasks to be run
        tbb::task_list tasks;
        
        int start = 0;
        
        if (not even_pairs)
        {
            //create a child task to run replica 0 on its own
            RunReplicaTask &task = *new( this->allocate_child() )
                                     RunReplicaTask(replicas[0], record_stats);
            
            tasks.push_back(task);
            
            start = 1;
        }
        
        for (int i=start; i<nreplicas-1; i+=2)
        {
            //create a child task that runs a pair of replicas, and then
            //performs a replica exchange test between the pair
            RunSwapTask &task = *new( this->allocate_child() )
                                    RunSwapTask(replicas[i], replicas[i+1],
                                                to_swap, naccept, nreject, swap_mutex);
            
            tasks.push_back(task);
        }
        
        //make sure that any dangling end replica is also run
        if ( (even_pairs and nreplicas % 2 == 1) or
             ((not even_pairs) and nreplicas % 2 == 0) )
        {
            RunReplicaTask &task = *new( this->allocate_child() )
                                      RunReplicaTask(replicas[nreplicas-1], record_stats);
            
            tasks.push_back(task);
        }

        //set the reference count to nreplicas + 1 (add one as we will wait
        //for all replicas to finish)
        this->set_ref_count( tasks.count() + 1 );
        
        //spawn all child tasks and wait for them to complete
        this->spawn_and_wait_for_all(tasks);
    }
};

/** Internal function that performs a single block of sampling on all
    replicas (recording statistics if 'record_stats' is true), using
    threaded building blocks to parallelise the tasks, 
    and then performing replica exchange moves between
    pairs */
void RepExMove2::performMove(Replicas &replicas, bool record_stats)
{
    if (replicas.nReplicas() == 0)
        return;

    //will we swap even pairs or odd pairs?
    bool even_pairs = true;
    
    if (replicas.nReplicas() > 2)
        even_pairs = rangenerator.randBool();

    //create space to hold all of the post-move results
    QVector<SupraSubSystemPtr> results( replicas.nReplicas() );
    
    //create space to hold all of the pairs that should be swapped
    QVector< std::pair<int,int> > to_swap;
    
    //create a mutex to serialise access to 'to_swap', 'naccept' and 'nreject'
    QMutex swap_mutex;
    
    for (int i=0; i<replicas.nReplicas(); ++i)
    {
        results[i] = replicas[i];
    }

    if (replicas.nReplicas() == 1 or disable_swaps)
    {
        //just run each replica in parallel without any swaps
        RunReplicasTask &parent_task = *new( tbb::task::allocate_root() )
                                            RunReplicasTask(results, record_stats);
        
        tbb::task::spawn_root_and_wait(parent_task);
    }
    else
    {
        //we will perform replica exchange moves between pairs
        RunReplicasWithSwapsTask &parent_task = *new( tbb::task::allocate_root() )
                                RunReplicasWithSwapsTask(results, record_stats,
                                                         generator().randBool(),
                                                         &to_swap, &naccept, &nreject,
                                                         &swap_mutex);

        tbb::task::spawn_root_and_wait(parent_task);
    }
    
    //copy the results back into the replicas
    int nreplicas = replicas.count();
    
    for (int i=0; i<nreplicas; ++i)
    {
        submoves[i] = subsims[i].result().subMoves().asA<SameSupraSubMoves>()[0]
                                         .asA<RepExSubMove>();
        
        replicas.setReplica( i, results[i].read().subSystem() );
    }

    //clear up the results to save memory
    results.clear();
    
    //now perform any swaps
    foreach( const std::pair<int,int> &swap, to_swap )
    {
        replicas.swapSystems(swap.first, swap.second, swap_monitors);
    }

    to_swap.clear();
    
    //now collect any necessary statistics
    if (record_stats)
        replicas.collectSupraStats();
}

/** Perform 'nmoves' replica exchange moves (block of sampling for all
    replicas, then replica exchange test between all pairs),
    of the system 'system' (which must be a Replicas object), optionally
    recording statistics if 'record_stats' is true 
    
    \throw SireError::invalid_cast
*/
void RepExMove2::move(SupraSystem &system, int nmoves, bool record_stats)
{
    Replicas &replicas = system.asA<Replicas>();
    
    if (replicas.nReplicas() == 0 or nmoves <= 0)
        return;
    
    SupraSystemPtr old_replicas = replicas.clone();
    SupraMovePtr old_state = this->clone();
    
    try
    {
        for (int i=0; i<nmoves; ++i)
        {
            this->performMove(replicas, record_stats);
        }

        SupraMove::incrementNMoves(nmoves);
    }
    catch(...)
    {
        replicas.copy(*old_replicas);
        this->copy(*old_state);
        throw;
    }
}

const char* RepExMove2::typeName()
{
    return QMetaType::typeName( qMetaTypeId<RepExMove2>() );
}


