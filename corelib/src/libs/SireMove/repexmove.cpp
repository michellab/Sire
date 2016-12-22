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

#include <QPair>

#include "repexmove.h"

#include "replica.h"
#include "replicas.h"

#include "suprasystem.h"

#include "suprasubsim.h"

#include "SireCluster/cluster.h"
#include "SireCluster/node.h"
#include "SireCluster/nodes.h"

#include "SireUnits/units.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireError/errors.h"

#include <QDebug>

using namespace SireMove;
using namespace SireMaths;
using namespace SireCluster;
using namespace SireBase;
using namespace SireVol;
using namespace SireUnits;
using namespace SireUnits::Dimension;
using namespace SireStream;

////////////
//////////// Implementation of RepExSubMove
////////////

static const RegisterMetaType<RepExSubMove> r_repexsubmove;

/** Serialise to a binary datastream */
QDataStream SIREMOVE_EXPORT &operator<<(QDataStream &ds,
                                        const RepExSubMove &repexsubmove)
{
    writeHeader(ds, r_repexsubmove, 2);
    
    SharedDataStream sds(ds);
    
    sds << repexsubmove.partner_properties
        << repexsubmove.have_new_vals;
        
    if (repexsubmove.have_new_vals)
    {
        sds << repexsubmove.new_volume_i.to(angstrom3)
            << repexsubmove.new_energy_i.to(kcal_per_mol)
            << repexsubmove.new_volume_j.to(angstrom3) 
            << repexsubmove.new_energy_j.to(kcal_per_mol);
    }

    sds << repexsubmove.need_volume;
    
    sds << static_cast<const SupraSubMove&>(repexsubmove);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOVE_EXPORT &operator>>(QDataStream &ds, RepExSubMove &repexsubmove)
{
    VersionID v = readHeader(ds, r_repexsubmove);
    
    if (v == 1 or v == 2)
    {
        RepExSubMove new_submove;
    
        SharedDataStream sds(ds);
        
        sds >> new_submove.partner_properties
            >> new_submove.have_new_vals;
            
        if (new_submove.have_new_vals)
        {
            double new_volume_i, new_energy_i, new_volume_j, new_energy_j;
            
            sds >> new_volume_i >> new_energy_i
                >> new_volume_j >> new_energy_j;

            new_submove.new_volume_i = new_volume_i * angstrom3;
            new_submove.new_energy_i = new_energy_i * kcal_per_mol;
            
            new_submove.new_volume_j = new_volume_j * angstrom3;
            new_submove.new_energy_j = new_energy_j * kcal_per_mol;
        }
        
        if (v == 2)
            sds >> new_submove.need_volume;
        else
            new_submove.need_volume = true;
        
        sds >> static_cast<SupraSubMove&>(new_submove);
        
        //check that all of the partner properties are valid...
        for (QList< QPair<quint32,QVariant> >::const_iterator 
                                it = new_submove.partner_properties.constBegin();
             it != new_submove.partner_properties.constEnd();
             ++it)
        {
            switch (it->first)
            {
                case RepExSubMove::LAMBDA_VALUE:
                case RepExSubMove::NRG_COMPONENT:
                case RepExSubMove::SPACE_PROPERTY:
                    break;
                    
                default:
                    throw version_error( QObject::tr(
                        "Version 1 of SireMove::RepExSubMove does not support "
                        "the partner property with ID %2.")
                            .arg(it->first), CODELOC );
            }
        }
        
        repexsubmove = new_submove;
    }
    else
        throw version_error(v, "1,2", r_repexsubmove, CODELOC);
        
    return ds;
}

/** Constructor */
RepExSubMove::RepExSubMove()
             : ConcreteProperty<RepExSubMove,SupraSubMove>(),
               new_volume_i(0), new_energy_i(0),
               new_volume_j(0), new_energy_j(0),
               have_new_vals(false), need_volume(false)
{}

/** Internal function used to add a property of our partner replica
    that differs from the value in this replica */
template<class T>
void RepExSubMove::addPartnerProperty(quint32 property, const T &value)
{
    partner_properties.append( QPair<quint32,QVariant>(property,
                                                       QVariant::fromValue<T>(value)) );
}

/** Construct the sub-move that will perform a move on 'replica_a', after
    which it will then calculate the values necessary to test the 
    swap from 'replica_a' to 'replica_b' */
RepExSubMove::RepExSubMove(const Replica &replica_a, const Replica &replica_b)
             : ConcreteProperty<RepExSubMove,SupraSubMove>(),
               new_volume_i(0), new_energy_i(0),
               new_volume_j(0), new_energy_j(0),
               have_new_vals(false)
{
    need_volume = replica_a.ensemble().isConstantPressure() and
                  replica_b.ensemble().isConstantPressure();

    if (replica_b.lambdaValue() != replica_a.lambdaValue())
    {
        addPartnerProperty( LAMBDA_VALUE, replica_b.lambdaValue() );
    }

    if (replica_b.energyComponent() != replica_a.energyComponent())
    {
        addPartnerProperty( NRG_COMPONENT, replica_b.energyComponent() );
    }
    
    if (need_volume)
    {
        if (replica_b.spaceProperty() != replica_a.spaceProperty())
        {
            addPartnerProperty( SPACE_PROPERTY, replica_b.spaceProperty() );
        }
    }
}

/** Copy constructor */
RepExSubMove::RepExSubMove(const RepExSubMove &other)
             : ConcreteProperty<RepExSubMove,SupraSubMove>(other),
               new_volume_i(other.new_volume_i), new_energy_i(other.new_energy_i),
               new_volume_j(other.new_volume_j), new_energy_j(other.new_energy_j),
               partner_properties(other.partner_properties),
               have_new_vals(other.have_new_vals), need_volume(other.need_volume)
{}

/** Destructor */
RepExSubMove::~RepExSubMove()
{}

/** Copy assignment operator */
RepExSubMove& RepExSubMove::operator=(const RepExSubMove &other)
{
    if (this != &other)
    {
        partner_properties = other.partner_properties;

        new_volume_i = other.new_volume_i;
        new_energy_i = other.new_energy_i;
        
        new_volume_j = other.new_volume_j;
        new_energy_j = other.new_energy_j;
        
        have_new_vals = other.have_new_vals;
        need_volume = other.need_volume;
        
        SupraSubMove::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool RepExSubMove::operator==(const RepExSubMove &other) const
{
    return (this == &other) or
           ( have_new_vals == other.have_new_vals and
             need_volume == other.need_volume and
             new_volume_i == other.new_volume_i and
             new_energy_i == other.new_energy_i and
             new_volume_j == other.new_volume_j and
             new_energy_j == other.new_energy_j and
             SupraSubMove::operator==(other) );
}

/** Comparison operator */
bool RepExSubMove::operator!=(const RepExSubMove &other) const
{
    return not this->operator==(other);
}

/** Return a string representation of this move */
QString RepExSubMove::toString() const
{
    if (have_new_vals)
    {
        if (need_volume)
        {
            return QObject::tr( "RepExSubMove( E_i = %1 kcal mol-1, V_i = %2 A^3 : "
                                "E_j = %3 kcal mol-1, V_j = %4 A^3 )")
                        .arg( new_energy_i.to(kcal_per_mol) )
                        .arg( new_volume_i.to(angstrom3) )
                        .arg( new_energy_j.to(kcal_per_mol) )
                        .arg( new_volume_j.to(angstrom3) );
        }
        else
        {
            return QObject::tr( "RepExSubMove( E_i = %1 kcal mol-1 : "
                                "E_j = %2 kcal mol-1 )")
                        .arg( new_energy_i.to(kcal_per_mol) )
                        .arg( new_energy_j.to(kcal_per_mol) );
        }
    }
    else
        return QObject::tr( "RepExSubMove()" );
}

static void throwNoValues(const char *value, const QString &codeloc)
{
    throw SireError::invalid_state( QObject::tr(
        "Cannot get the value of %1 as new values have not been evaluated!")
            .arg(value), codeloc );
}

static void throwNoVolume(const char *value, const QString &codeloc)
{
    throw SireError::invalid_state( QObject::tr(
        "Cannot get the volume %1 as the move indicated that volumes weren't necessary!")
            .arg(value), codeloc );
}

/** Return the energy of the replica in its normal state at the 
    end of the block of moves */
MolarEnergy RepExSubMove::energy_i() const
{
    if (not have_new_vals)
        ::throwNoValues( "E_i", CODELOC );
        
    return new_energy_i;
}

/** Return the volume of the replica in its normal state at
    the end of the block of moves */
Volume RepExSubMove::volume_i() const
{
    if (not have_new_vals)
        ::throwNoValues( "V_i", CODELOC );
    
    if (not need_volume)
        ::throwNoVolume( "V_i", CODELOC );
                
    return new_volume_i;
}

/** Return the energy of the replica in its partner state 
    at the end of the block of moves */
MolarEnergy RepExSubMove::energy_j() const
{
    if (not have_new_vals)
        ::throwNoValues( "E_j", CODELOC );
        
    return new_energy_j;
}

/** Return the volume of the replica in its partner state
    at the end of the block of moves */
Volume RepExSubMove::volume_j() const
{
    if (not have_new_vals)
        ::throwNoValues( "V_j", CODELOC );
    
    if (not need_volume)
        ::throwNoVolume( "V_j", CODELOC );
        
    return new_volume_j;
}

/** Internal function used to extract a value of type 'T' from
    the passed QVariant, throwing an exception if this is not possible
    
    \throw SireError::invalid_cast
*/
template<class T>
static T convert(const QVariant &value)
{
    if (not value.canConvert<T>())
        throw SireError::invalid_cast( QObject::tr(
            "Cannot apply a deferred command as the argument of type %1 "
            "cannot be cast to a value of type %2.")
                .arg( typeid(T).name() )
                .arg( QVariant::typeToName(value.type()) ), CODELOC );
                
    return value.value<T>();
}

/** Evaluate the energy and volume of this replica after
    it has been swapped into its partner state */
void RepExSubMove::evaluateSwappedState(const Replica &replica)
{
    if (have_new_vals)
        return;

    Replica state_i = replica;

    if (state_i.isPacked())
        state_i.unpack();

    //get the energy and volume at this state
    new_energy_i = state_i.energy();

    if (need_volume)
        new_volume_i = state_i.volume();

    if (partner_properties.isEmpty())
    {
        new_volume_j = new_volume_i;
        new_energy_j = new_energy_i;
    }
    else
    {
        System state_j = state_i.subSystem();
        
        Symbol nrg_component = state_i.energyComponent();
        PropertyName space_property = state_i.spaceProperty();
        
        for (QList< QPair<quint32,QVariant> >::const_iterator 
                                                it = partner_properties.constBegin();
             it != partner_properties.constEnd();
             ++it)
        {
            switch (it->first)
            {
                case LAMBDA_VALUE:
                {
                    if (state_i.lambdaComponent().isNull())
                        throw SireError::incompatible_error( QObject::tr(
                            "Cannot set the lambda value for a replica that doesn't "
                            "have a lambda component!"), CODELOC );
                
                    //set a new lambda value
                    state_j.setComponent( state_i.lambdaComponent(),
                                          ::convert<double>(it->second) );
                    break;
                }   
                case NRG_COMPONENT:
                    //set a new Hamiltonian (represented by the component)
                    nrg_component = ::convert<Symbol>(it->second);
                    break;
                    
                case SPACE_PROPERTY:
                    //set a new space property
                    space_property = ::convert<PropertyName>(it->second);
                    break;
                    
                default:
                    throw SireError::unsupported( QObject::tr(
                        "A request was made of an unsuppoted action in RepExSubMove. "
                        "The action with ID %1 was requested, but this is not "
                        "supported with this version of RepExSubMove.")
                            .arg(it->first), CODELOC );
            }
        }

        new_energy_j = state_j.energy(nrg_component);

        if (need_volume)
            new_volume_j = state_j.property(space_property).asA<Space>().volume();
    }

    have_new_vals = true;
}

/** Perform the sub-moves on the passed sub-system */
void RepExSubMove::move(SupraSubSystem &system, int n_supra_moves,
                        int n_supra_moves_per_block, bool record_stats)
{
    //replica exchange moves work only with Replica objects
    Replica &replica = system.asA<Replica>();

    SupraSubSystemPtr old_replica = replica.clone();
    SupraSubMovePtr old_state = this->clone();

    try
    {
        have_new_vals = false;

        if (n_supra_moves <= 0)
        {
            if (n_supra_moves_per_block <= 0)
                this->evaluateSwappedState(replica);
        
            return;
        }

        //unpack the system, if necessary
        bool replica_was_packed = replica.isPacked();
        replica.unpack();
    
        //perform the moves
        for (int i=0; i<n_supra_moves; ++i)
        {
            replica.subMove(record_stats);
        }

        //if we have finished a block of sub-moves, then collect
        //the information necessary to perform the replica exchange test
        if (n_supra_moves >= n_supra_moves_per_block)
            this->evaluateSwappedState(replica);
        
        //repack the system, if necessary
        if (replica_was_packed)
            replica.pack();
    }
    catch(...)
    {
        replica.copy(*old_replica);
        this->copy(*old_state);
        
        throw;
    }
}

const char* RepExSubMove::typeName()
{
    return QMetaType::typeName( qMetaTypeId<RepExSubMove>() );
}

////////////
//////////// Implementation of RepExMove
////////////

static const RegisterMetaType<RepExMove> r_repexmove;

/** Serialise to a binary datastream */
QDataStream SIREMOVE_EXPORT &operator<<(QDataStream &ds, const RepExMove &repexmove)
{
    writeHeader(ds, r_repexmove, 3);

    SharedDataStream sds(ds);
    
    sds << repexmove.rangenerator
        << repexmove.naccept
        << repexmove.nreject
        << repexmove.swap_monitors
        << repexmove.disable_swaps
        << static_cast<const SupraMove&>(repexmove);
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOVE_EXPORT &operator>>(QDataStream &ds, RepExMove &repexmove)
{
    VersionID v = readHeader(ds, r_repexmove);

    repexmove.disable_swaps = false;

    if (v == 3)
    {
        SharedDataStream sds(ds);
        
        sds >> repexmove.rangenerator
            >> repexmove.naccept
            >> repexmove.nreject
            >> repexmove.swap_monitors
            >> repexmove.disable_swaps
            >> static_cast<SupraMove&>(repexmove);
    }
    else if (v == 2)
    {
        SharedDataStream sds(ds);
        
        sds >> repexmove.rangenerator
            >> repexmove.naccept
            >> repexmove.nreject
            >> repexmove.swap_monitors
            >> static_cast<SupraMove&>(repexmove);
    }
    else if (v == 1)
    {
        SharedDataStream sds(ds);

        quint32 nmoves;
        
        sds >> repexmove.rangenerator
            >> nmoves
            >> repexmove.naccept
            >> repexmove.nreject;
            
        repexmove.swap_monitors = false;
    }
    else
        throw version_error(v, "1-3", r_repexmove, CODELOC);
        
    return ds;
}

/** Constructor */
RepExMove::RepExMove()
          : ConcreteProperty<RepExMove,SupraMove>(),
            naccept(0), nreject(0), swap_monitors(false), disable_swaps(false)
{}

/** Copy constructor */
RepExMove::RepExMove(const RepExMove &other)
          : ConcreteProperty<RepExMove,SupraMove>(other),
            naccept(other.naccept), nreject(other.nreject),
            swap_monitors(other.swap_monitors),
            disable_swaps(other.disable_swaps)
{}

/** Destructor */
RepExMove::~RepExMove()
{}

/** Copy assignment operator */
RepExMove& RepExMove::operator=(const RepExMove &other)
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
bool RepExMove::operator==(const RepExMove &other) const
{
    return (this == &other) or
           (naccept == other.naccept and nreject == other.nreject and
            swap_monitors == other.swap_monitors and
            disable_swaps == other.disable_swaps and SupraMove::operator==(other));
}

/** Comparison operator */
bool RepExMove::operator!=(const RepExMove &other) const
{
    return not this->operator==(other);
}

/** Return the total number of accepted replica exchange tests */
int RepExMove::nAttempted() const
{
    return naccept + nreject;
}

/** Return the total number of accepted replica exchange tests */
int RepExMove::nAccepted() const
{
    return naccept;
}

/** Return the total number of rejected replica exchange tests */
int RepExMove::nRejected() const
{
    return nreject;
}

/** Return whether or not swap moves are disabled */
bool RepExMove::swapMovesDisabled() const
{
    return disable_swaps;
}

/** Set disabling of swap moves */
void RepExMove::setDisableSwaps(bool disable)
{
    disable_swaps = disable;
}

/** Return the average acceptance ratio of the replica exchange
    tests over all replicas */
double RepExMove::acceptanceRatio() const
{
    if (this->nAttempted() > 0)
    {
        return double(this->nAccepted()) / double(this->nAttempted());
    }
    else
        return 0;
}

/** Return a string representation of this move */
QString RepExMove::toString() const
{
    return QObject::tr("RepExMove( %1 accepted, %2 rejected : %3 %% )")
                .arg(this->nAccepted())
                .arg(this->nRejected())
                .arg(100 * this->acceptanceRatio());
}

/** Clear the move statistics */
void RepExMove::clearStatistics()
{
    naccept = 0;
    nreject = 0;
    SupraMove::clearStatistics();
}

/** Set the random number generator used for the replica exchange tests */
void RepExMove::setGenerator(const RanGenerator &generator)
{
    rangenerator = generator;
}

/** Return the random number generator used for the replica exchange tests */
const RanGenerator& RepExMove::generator() const
{
    return rangenerator;
}

/** Set whether or not to swap the system monitors when we swap the systems */
void RepExMove::setSwapMonitors(bool swap)
{
    swap_monitors = swap;
}

/** Internal function used to submit the simulation in 'replica' */
static SupraSubSim submitSimulation(Nodes &nodes, const Replica &replica,
                                    bool record_stats)
{
    Node node = nodes.getNode();

    return SupraSubSim::run( node, replica, RepExSubMove(), 1, record_stats );
}

/** Internal function used to submit the two simulations for the replicas
    'replica_a' and 'replica_b', telling each simulation to follow the 
    simulation by evaluating the values necessary to perform a swap test
    between these two replicas */
static QPair<SupraSubSim,SupraSubSim> submitSimulation(
                                    Nodes &nodes, const Replica &replica_a,
                                    const Replica &replica_b, bool record_stats)
{
    QPair<SupraSubSim,SupraSubSim> sims;

    //submit each simulation in a separate scope so that
    //the node is released once the simulation has finished
    {
        Node node_a = nodes.getNode();
        
        sims.first = SupraSubSim::run( node_a, replica_a,
                                       RepExSubMove(replica_a, replica_b),
                                       1, record_stats );
    }
     
    Node node_b = nodes.getNode();

    sims.second = SupraSubSim::run( node_b, replica_b,
                                    RepExSubMove(replica_b, replica_a),
                                                 1, record_stats );

    return sims;
}

/** Internal function used to submit all of the replica simulations in 'replicas'
    to the nodes 'nodes', returning an array of running simulations. This
    will set up the moves to swap even pairs if 'even_pairs' is true, or
    odd pairs if 'even_pairs' is false. Statistics will be recorded during
    the simulation only if 'record_stats' is true */
static QVector<SupraSubSim> submitSimulations(Nodes &nodes, Replicas &replicas,
                                              bool even_pairs, bool record_stats)
{
    int nreplicas = replicas.nReplicas();

    QVector<SupraSubSim> subsims(nreplicas);
    
    //pair the replicas up...
    //do we swap the even pairs or the odd pairs?
    int start = 1;

    if (even_pairs)
        start = 0;

    if (start == 1)
    {
        //the first replica hasn't got a partner - start it on its own
        subsims[0] = ::submitSimulation(nodes, replicas[0], record_stats);
    }

    if ( (nreplicas-start) % 2 == 1 )
    {
        //the last replica hasn't got a partner - start it on its own too
        subsims[nreplicas-1] = ::submitSimulation(nodes, replicas[nreplicas-1],
                                                  record_stats );
    }

    for (int i=start; i<nreplicas-1; i+=2)
    {
        //submit the simulations for replica i and replica i+1
        QPair<SupraSubSim,SupraSubSim> sims = ::submitSimulation( nodes, replicas[i], 
                                                                  replicas[i+1],
                                                                  record_stats );
        
        subsims[i] = sims.first;
        subsims[i+1] = sims.second;
    }

    return subsims;
}

/** Internal function used to wait until all of the simulations
    in 'subsims' have finished, and to optionally restart broken
    simulations up to 'max_tries' times using the nodes in 'nodes' */
static void waitUntilFinished(Nodes &nodes, QVector<SupraSubSim> &subsims, int max_tries)
{
    bool all_finished = false;
    int ntries = 0;
    
    int nreplicas = subsims.count();
    
    while (not all_finished)
    {
        ++ntries;
        
        if (ntries > max_tries)
            return;
            
        all_finished = true;
        
        for (int i=0; i<nreplicas; ++i)
        {
            SupraSubSim &subsim = subsims[i];
        
            subsim.wait();
            
            if (subsim.isError() or subsim.wasAborted())
            {
                //resubmit this calculation
                Node node = nodes.getNode();
                subsim = SupraSubSim::run(node, subsim.input());
            
                all_finished = false;
            }
            else if (not subsim.hasFinished())
            {
                //continue the calculation from where it finished
                SupraSubSimPacket simpacket = subsim.result();
                
                Node node = nodes.getNode();
                subsim = SupraSubSim::run( node, simpacket );
                
                all_finished = false;
            }
        }
        
        if (not all_finished)
        {
            //wait for the resubmitted calculations to finish
            for (int i=0; i<nreplicas; ++i)
            {
                subsims[i].wait();
            }
        }
    }
}

/** Internal function used to test the passed pair of replicas - 
    this returns whether or not the test has passed */
bool RepExMove::testPair(const Replica &replica_a, const RepExSubMove &move_a, 
                         const Replica &replica_b, const RepExSubMove &move_b) const
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
        MolarEnergy H_a_i = move_a.energy_i();
        MolarEnergy H_a_j = move_a.energy_j();
        
        MolarEnergy H_b_i = move_b.energy_i();
        MolarEnergy H_b_j = move_b.energy_j();
        
        Volume V_a_i(0);
        Volume V_a_j(0);

        Volume V_b_i(0);
        Volume V_b_j(0);
        
        if (need_pv)
        {
            V_a_i = move_a.volume_i();
            V_a_j = move_a.volume_j();
            
            V_b_i = move_b.volume_i();
            V_b_j = move_b.volume_j();
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

/** Internal function used to test and swap all pairs of replicas */
void RepExMove::testAndSwap(Replicas &replicas, const QVector<RepExSubMove> &submoves,
                            bool even_pairs, bool record_stats)
{
    int nreplicas = replicas.nReplicas();

    if (nreplicas > 1)
    {
        int start = 1;
        
        if (even_pairs)
            start = 0;
            
        //loop over all pairs
        for (int i=start; i<nreplicas-1; i+=2)
        {
            qDebug() << "Test replicas" << i << (i+1);
        
            if (this->testPair(replicas[i], submoves.at(i),
                               replicas[i+1], submoves.at(i+1) ))
            {
                //swap the replicas
                replicas.swapSystems(i, i+1, swap_monitors);
                ++naccept;
            }
            else
                ++nreject;
        }
    }
}

/** Internal function that performs a single block of sampling on all
    replicas (recording statistics if 'record_stats' is true), using the
    nodes in 'nodes', and then performing replica exchange moves between
    pairs */
void RepExMove::performMove(Nodes &nodes, Replicas &replicas, bool record_stats)
{
    //will we swap even pairs or odd pairs?
    bool even_pairs = true;
    
    if (replicas.nReplicas() > 2)
        even_pairs = rangenerator.randBool();

    //submit all of the simulations
    QVector<SupraSubSim> subsims = ::submitSimulations(nodes, replicas,
                                                       even_pairs, record_stats);
        
    //wait for all of the simulations to finish (retrying broken simulations
    //just five times)
    ::waitUntilFinished(nodes, subsims, 5);
    
    //copy the results back into the replicas
    int nreplicas = replicas.count();
    
    QVector<RepExSubMove> submoves(nreplicas);
    submoves.squeeze();
    
    for (int i=0; i<nreplicas; ++i)
    {
        submoves[i] = subsims[i].result().subMoves().asA<SameSupraSubMoves>()[0]
                                         .asA<RepExSubMove>();
        
        replicas.setReplica( i, subsims[i].result().subSystem() );
    }
    
    //get rid of the simulation handles, as they are no longer needed
    //(and we need to save memory if we can)
    subsims = QVector<SupraSubSim>();
    
    //now perform all of the replica exchange tests
    if (not disable_swaps)
        this->testAndSwap(replicas, submoves, even_pairs, record_stats);
    
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
void RepExMove::move(SupraSystem &system, int nmoves, bool record_stats)
{
    Replicas &replicas = system.asA<Replicas>();
    
    if (replicas.nReplicas() == 0 or nmoves <= 0)
        return;
    
    SupraSystemPtr old_replicas = replicas.clone();
    SupraMovePtr old_state = this->clone();
    
    try
    {
        //try to get as many nodes as possible to run the moves
        Nodes nodes = Cluster::getNodes( replicas.nReplicas() - 1, 5000 );
        
        ///hold this_thread in a local scope to ensure it is deleted before 'nodes'
        {
            ThisThread this_thread = nodes.borrowThisThread();
        
            for (int i=0; i<nmoves; ++i)
            {
                this->performMove(nodes, replicas, record_stats);
            }

            SupraMove::incrementNMoves(nmoves);
        }
    }
    catch(...)
    {
        replicas.copy(*old_replicas);
        this->copy(*old_state);
        throw;
    }
}

const char* RepExMove::typeName()
{
    return QMetaType::typeName( qMetaTypeId<RepExMove>() );
}
