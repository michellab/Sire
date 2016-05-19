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

#include <QVector>
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
    RunReplicaTask(SupraSubSystemPtr &replica_, bool record_stats_)
        : replica(&replica_), record_stats(record_stats_)
    {}
    
    tbb::task* execute()
    {
        if (replica == 0)
            return 0;

        //actually perform the replica moves, recording statistics
        //if record_stats is true
        replica->edit().subMove(record_stats);
        
        return 0;
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
     : replicas( replicas_.data() ), nreplicas(replicas_.count()), record_stats( record_stats_ )
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
        
        return 0;
    }
};

/** Internal function used to test the passed pair of replicas - 
    this returns whether or not the test has passed */
bool replicaTest(Replica &replica_a, Replica &replica_b,
                 const RanGenerator &rangenerator)
{
    //get the ensembles of the two replicas
    const Ensemble &ensemble_a = replica_a.ensemble();
    const Ensemble &ensemble_b = replica_b.ensemble();
    
    if ( (ensemble_a.isNVT() and ensemble_a.isNVT()) or 
         (ensemble_b.isNPT() and ensemble_b.isNPT()) )
    {
        bool need_pv = (ensemble_a.isNPT() and ensemble_b.isNPT());
    
        //get the values of the thermodynamic parameters
        double beta_a = 1.0 / (k_boltz * ensemble_a.temperature()).value();
        double beta_b = 1.0 / (k_boltz * ensemble_b.temperature()).value();
        
        Pressure p_a(0);
        Pressure p_b(0);
        
        if (need_pv)
        {
            p_a = ensemble_a.pressure();
            p_b = ensemble_b.pressure();
        }
        
        //now get the values of the system properties at their current state,
        //and at their swapped states
        double H_a_i = replica_a.energy().value();
        double H_b_i = replica_b.energy().value();

        double H_a_j = H_a_i;
        double H_b_j = H_b_i;

        if (replica_a.lambdaValue() != replica_b.lambdaValue() or
            replica_a.energyComponent() != replica_b.energyComponent())
        {
            //there will be a change in energy associated with the swap
            System swapped = replica_a.subSystem();

            //evaluate the energy of replica A swapped into the B state
            swapped.setComponent(replica_a.lambdaComponent(), replica_b.lambdaValue());
            H_a_j = swapped.energy( replica_b.energyComponent() ).value();
            
            //evaluate the energy of replica B swapped into the A state
            swapped = replica_b.subSystem();
            swapped.setComponent(replica_b.lambdaComponent(), replica_a.lambdaValue());
            H_b_j = swapped.energy( replica_a.energyComponent() ).value();
        }
        
        double V_a_i(0);
        double V_a_j(0);

        double V_b_i(0);
        double V_b_j(0);
        
        if (need_pv)
        {
            V_a_i = replica_a.volume().value();
            V_a_j = V_a_i;
            
            V_b_i = replica_b.volume().value();
            V_b_j = V_b_i;

            if (replica_a.spaceProperty() != replica_b.spaceProperty())
            {
                //the space property changes, so volume could change
                V_a_j = replica_a.subSystem().property(replica_b.spaceProperty())
                                 .asA<Space>().volume().value();

                V_b_j = replica_b.subSystem().property(replica_a.spaceProperty())
                                 .asA<Space>().volume().value();
            }
        }
        
        //now calculate delta needed for the Monte Carlo test
        //
        //  For derivation see Appendix C of Christopher Woods' thesis
        //   (or original replica exchange literature of course!)
        //
        //  delta = beta_b * [ H_b_i - H_b_j + P_b (V_b_i - V_b_j) ] + 
        //          beta_a * [ H_a_i - H_a_j + P_a (V_a_i - V_a_j) ]
        
        double delta = beta_b * ( H_b_i - H_b_j + p_b*(V_b_i - V_b_j) ) +
                       beta_a * ( H_a_i - H_a_j + p_a*(V_a_i - V_a_j) );
        
        bool move_passed = ( delta > 0 or (std::exp(delta) >= rangenerator.rand()) );
        
        return move_passed;
    }
    else
    {
        throw SireError::incompatible_error( QObject::tr(
            "There is no available replica exchange test that allows tests between "
            "replicas with ensembles %1 and %2.")
                .arg(ensemble_a.toString(), ensemble_b.toString()), CODELOC );
                
    }

    return false;
}

/** This is a task that runs two replicas, and then tests whether or not
    they should be swapped */
class RunSwapTask : public tbb::task
{
    SupraSubSystemPtr *replica_a;
    SupraSubSystemPtr *replica_b;
    bool record_stats;
    const RanGenerator *rangenerator;
    QVector< std::pair<int,int> > *to_swap;
    int *naccept;
    int *nreject;
    QMutex *swap_mutex;
    int idx_a;
    int idx_b;

public:
    RunSwapTask(SupraSubSystemPtr &replica_a_, SupraSubSystemPtr &replica_b_,
                bool record_stats_, const RanGenerator *rangenerator_,
                QVector< std::pair<int,int> > *to_swap_,
                int *naccept_, int *nreject_, QMutex *swap_mutex_,
                int idx_a_, int idx_b_)
        : replica_a(&replica_a_), replica_b(&replica_b_), record_stats(record_stats_),
          rangenerator(rangenerator_),
          to_swap(to_swap_), naccept(naccept_), nreject(nreject_), swap_mutex(swap_mutex_),
          idx_a(idx_a_), idx_b(idx_b_)
    {}
    
    tbb::task* execute()
    {
        if (replica_a == 0 or replica_b == 0 or replica_a == replica_b)
            return 0;
        
        //set the reference count to 2 + 1 (two child tasks, plus add one as we will wait
        //for both replicas to finish)
        this->set_ref_count( 2 + 1 );
        
        //spawn tasks for both replicas
        RunReplicaTask &task_a = *new( this->allocate_child() )
                                       RunReplicaTask(*replica_a, record_stats);
        
        RunReplicaTask &task_b = *new( this->allocate_child() )
                                       RunReplicaTask(*replica_b, record_stats);
        
        this->spawn(task_a);
        this->spawn_and_wait_for_all(task_b);
        
        //now apply the replica exchange test
        bool test_passed = replicaTest(replica_a->edit().asA<Replica>(),
                                       replica_b->edit().asA<Replica>(),
                                       *rangenerator);
        
        //serialise access to to_swap, naccept and nreject
        QMutexLocker lkr(swap_mutex);
        
        if (test_passed)
        {
            to_swap->append( std::pair<int,int>(idx_a,idx_b) );
            *naccept += 1;
        }
        else
        {
            *nreject += 1;
        }
        
        return 0;
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
    const RanGenerator *rangenerator;
    QVector< std::pair<int,int> > *to_swap;
    int *naccept;
    int *nreject;
    QMutex *swap_mutex;
    
public:
    RunReplicasWithSwapTask(QVector<SupraSubSystemPtr> &replicas_,
                            bool record_stats_,
                            const RanGenerator &rangenerator_,
                            QVector< std::pair<int,int> > *to_swap_,
                            int *naccept_, int *nreject_, QMutex *swap_mutex_)
     : replicas( replicas_.data() ), nreplicas(replicas_.count()), record_stats( record_stats_ ),
       even_pairs(rangenerator_.randBool()), rangenerator(&rangenerator_),
       to_swap(to_swap_), naccept(naccept_), nreject(nreject_),
       swap_mutex(swap_mutex_)
    {}
    
    tbb::task* execute()
    {
        if (nreplicas <= 0 or replicas == 0)
            return 0;
        
        //create a list of tasks to be run
        tbb::task_list tasks;
        int ntasks = 0;
        
        int start = 0;
        
        if (not even_pairs)
        {
            //create a child task to run replica 0 on its own
            RunReplicaTask &task = *new( this->allocate_child() )
                                     RunReplicaTask(replicas[0], record_stats);
            
            tasks.push_back(task);
            ntasks += 1;
            
            start = 1;
        }
        
        for (int i=start; i<nreplicas-1; i+=2)
        {
            //create a child task that runs a pair of replicas, and then
            //performs a replica exchange test between the pair
            RunSwapTask &task = *new( this->allocate_child() )
                                    RunSwapTask(replicas[i], replicas[i+1], record_stats,
                                                rangenerator, to_swap, naccept, nreject,
                                                swap_mutex, i, i+1);
            
            tasks.push_back(task);
            ntasks += 1;
        }
        
        //make sure that any dangling end replica is also run
        if ( (even_pairs and nreplicas % 2 == 1) or
             ((not even_pairs) and nreplicas % 2 == 0) )
        {
            RunReplicaTask &task = *new( this->allocate_child() )
                                      RunReplicaTask(replicas[nreplicas-1], record_stats);
            
            tasks.push_back(task);
            ntasks += 1;
        }

        //set the reference count to nreplicas + 1 (add one as we will wait
        //for all replicas to finish)
        this->set_ref_count( ntasks + 1 );
        
        //spawn all child tasks and wait for them to complete
        this->spawn_and_wait_for_all(tasks);
        
        return 0;
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
        int naccepted = 0;
        int nrejected = 0;
    
        //we will perform replica exchange moves between pairs
        RunReplicasWithSwapTask &parent_task = *new( tbb::task::allocate_root() )
                                RunReplicasWithSwapTask(results, record_stats,
                                                        generator(),
                                                        &to_swap, &naccepted, &nrejected,
                                                        &swap_mutex);

        tbb::task::spawn_root_and_wait(parent_task);
        
        naccept += naccepted;
        nreject += nrejected;
    }
    
    //copy the results back into the replicas
    int nreplicas = replicas.count();
    
    for (int i=0; i<nreplicas; ++i)
    {
        replicas.setReplica( i, results[i].read() );
    }

    //clear up the results to save memory
    results.clear();
    
    //now perform any swaps
    for( const std::pair<int,int> &swap : to_swap )
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


