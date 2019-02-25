/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2007  Christopher Woods
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

#include "weightedmoves.h"
#include "move.h"

#include "SireUnits/units.h"
#include "SireUnits/temperature.h"

#include "SireSystem/system.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireError/errors.h"

#include <QElapsedTimer>

using namespace SireMove;
using namespace SireMaths;
using namespace SireSystem;
using namespace SireBase;
using namespace SireUnits;
using namespace SireUnits::Dimension;
using namespace SireStream;

using boost::tuples::tuple;

static const RegisterMetaType<WeightedMoves> r_weightedmoves;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, 
                                        const WeightedMoves &weightedmoves)
{
    writeHeader(ds, r_weightedmoves, 3);
    
    SharedDataStream sds(ds);
    
    qint32 nmoves = weightedmoves.mvs.count();
    const tuple<MovePtr,double> *mvs_array = weightedmoves.mvs.constData();
    
    sds << nmoves;
    
    for (int i=0; i<nmoves; ++i)
    {
        sds << mvs_array[i].get<0>() << mvs_array[i].get<1>();
    }
    
    sds << weightedmoves.avgtimes << weightedmoves.rangenerator << weightedmoves.combined_space
        << static_cast<const Moves&>(weightedmoves);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, WeightedMoves &weightedmoves)
{
    VersionID v = readHeader(ds, r_weightedmoves);
    
    if (v == 3)
    {
        SharedDataStream sds(ds);
        
        qint32 nmoves;
        
        sds >> nmoves;
        
        QVector< tuple<MovePtr,double> > mvs(nmoves);
        tuple<MovePtr,double> *mvs_array = mvs.data();
        
        for (int i=0; i<nmoves; ++i)
        {
            MovePtr mv; double weight;
            
            sds >> mv >> weight;
            mvs_array[i] = tuple<MovePtr,double>(mv, weight);
        }
        
        weightedmoves.mvs = mvs;
        
        sds >> weightedmoves.avgtimes
            >> weightedmoves.rangenerator
            >> weightedmoves.combined_space
            >> static_cast<Moves&>(weightedmoves);
        
        weightedmoves.recalculateWeights();
    }
    else if (v == 1 or v == 2)
    {
        SharedDataStream sds(ds);
        
        qint32 nmoves;
        
        sds >> nmoves;
        
        QVector< tuple<MovePtr,double> > mvs(nmoves);
        tuple<MovePtr,double> *mvs_array = mvs.data();
        
        for (int i=0; i<nmoves; ++i)
        {
            MovePtr mv; double weight;
            
            sds >> mv >> weight;
            mvs_array[i] = tuple<MovePtr,double>(mv, weight);
        }
        
        weightedmoves.mvs = mvs;
        
        sds >> weightedmoves.rangenerator;
        
        if (v == 2)
            sds >> weightedmoves.combined_space;
        else
            weightedmoves.combined_space = PropertyName();
        
        sds >> static_cast<Moves&>(weightedmoves);
        
        weightedmoves.recalculateWeights();
    }
    else
        throw version_error(v, "1", r_weightedmoves, CODELOC);

    return ds;
}

/** Constructor */
WeightedMoves::WeightedMoves()
              : ConcreteProperty<WeightedMoves,Moves>(), maxweight(0)
{}

/** Copy constructor */
WeightedMoves::WeightedMoves(const WeightedMoves &other)
              : ConcreteProperty<WeightedMoves,Moves>(other),
                mvs(other.mvs), avgtimes(other.avgtimes),
                rangenerator(other.rangenerator),
                combined_space(other.combined_space),
                maxweight(other.maxweight)
{}

/** Destructor */
WeightedMoves::~WeightedMoves()
{}

/** Copy assignment operator */
WeightedMoves& WeightedMoves::operator=(const WeightedMoves &other)
{
    mvs = other.mvs;
    avgtimes = other.avgtimes;
    rangenerator = other.rangenerator;
    combined_space = other.combined_space;
    maxweight = other.maxweight;
    
    Moves::operator=(other);
    
    return *this;
}

inline bool compare( const QVector< boost::tuple<MovePtr,double> > &o1,
                     const QVector< boost::tuple<MovePtr,double> > &o2 )
{
    if (o1.constData() == o2.constData())
        return true;
    
    else if (o1.count() != o2.count())
        return false;
    
    else
    {
        for (int i=0; i<o1.count(); ++i)
        {
            if (o1.constData()[i].get<0>() != o2.constData()[i].get<0>() or
                o1.constData()[i].get<1>() != o2.constData()[i].get<1>())
                return false;
        }
        
        return true;
    }
}

/** Comparison operator */
bool WeightedMoves::operator==(const WeightedMoves &other) const
{
    return compare(mvs,other.mvs) and avgtimes == other.avgtimes and
           combined_space == other.combined_space and
           Moves::operator==(other);
}

/** Comparison operator */
bool WeightedMoves::operator!=(const WeightedMoves &other) const
{
    return not this->operator==(other);
}

/** Return a string representation */
QString WeightedMoves::toString() const
{
    QStringList moves;
    
    for (int i=0; i<mvs.count(); ++i)
    {
        moves.append( QObject::tr("  %1 : weight == %2, timing = %3 ms\n"
                                  "       %4").arg(i+1)
                          .arg(mvs.at(i).get<1>())
                          .arg((avgtimes.at(i).average() * nanosecond).to(millisecond))
                          .arg(mvs.at(i).get<0>()->toString()) );
    }
    
    return QObject::tr("WeightedMoves{\n%1\n}")
                .arg(moves.join("\n"));
}

/** Recalculate all of the weights */
void WeightedMoves::recalculateWeights()
{
    maxweight = 0;

    int nmoves = mvs.count();

    if (nmoves == 0)
        return;

    const tuple<MovePtr,double> *mvs_array = mvs.constData();
    
    for (int i=0; i<nmoves; ++i)
    {
        maxweight = qMax( maxweight, mvs_array[i].get<1>() );
    }
}

/** Set the random number generator used to pick moves, and also
    used by the moves themselves during the simulation */
void WeightedMoves::setGenerator(const RanGenerator &rangen)
{
    rangenerator = rangen;
    
    for (int i=0; i<mvs.count(); ++i)
    {
        mvs[i].get<0>().edit().setGenerator(rangen);
    } 
}

/** Return the random number generator used to pick moves. This
    may not be the same as the generator used by the moves themselves.
    To ensure it is the same, run;
    
    weightedmoves.setGenerator( weightedmoves.generator() )
*/
const RanGenerator& WeightedMoves::generator() const
{
    return rangenerator;
}

/** Add the move 'move' to the list of moves, with the weight 'weight' */
void WeightedMoves::add(const Move &move, double weight)
{
    if (weight <= 0)
        return;

    int nmoves = mvs.count();
    
    tuple<MovePtr,double> *mvs_array = mvs.data();
    
    for (int i=0; i<nmoves; ++i)
    {
        if (mvs_array[i].get<0>()->equals(move))
        {
            mvs_array[i].get<1>() += weight;
            this->recalculateWeights();
            return;
        }
    }
    
    mvs.append( tuple<MovePtr,double>(move, weight) );
    avgtimes.append( Average() );
    
    this->recalculateWeights();
}

/** Completely clear all of the move statistics */
void WeightedMoves::clearStatistics()
{
    tuple<MovePtr,double> *mvs_array = mvs.data();
    int nmoves = mvs.count();
    
    for (int i=0; i<nmoves; ++i)
    {
        mvs_array[i].get<0>().edit().clearStatistics();
    }
}

/** Perform 'nmoves' moves on the system 'system' and return the result */
System WeightedMoves::move(const System &system, int nmoves, bool record_stats)
{
    if (mvs.isEmpty())
        return system;

    WeightedMoves old_state(*this);
    System run_system(system);

    try
    {
        QElapsedTimer t;

        Moves::preCheck(run_system);
    
        int n = mvs.count();
        tuple<MovePtr,double> *mvs_array = mvs.data();

        if (n == 1)
        {
            t.start();
            mvs_array[0].get<0>().edit().move(run_system, nmoves, record_stats);
            
            qint64 ns = t.nsecsElapsed();
            
            for (int i=0; i<nmoves; ++i)
            {
                avgtimes[0].accumulate( (1.0*ns)/nmoves );
            }
            
            Moves::postCheck(run_system);
            
            return run_system;
        }

        for (int i=0; i<nmoves; ++i)
        {
            //use the von Neumann rejection method to choose a random move
            //using the supplied weights
            //  rejection method - choose random
            //  move, then choose random number from 0 to maxweight. If
            //  probability of the move <= the random number, then accept
            //  this move, else go back to the beginning and try again...
            while (true)
            {
                quint32 idx = generator().randInt(n-1);
    
                tuple<MovePtr,double> &move = mvs_array[idx];
    
                if ( generator().rand(maxweight) <= move.get<1>() )
                {
                    //use this move
                    t.start();
                    move.get<0>().edit().move(run_system, 1, record_stats);
                    qint64 ns = t.nsecsElapsed();
                    
                    avgtimes[idx].accumulate( 1.0*ns );
                    
                    break;
                }
            }
        }
        
        Moves::postCheck(run_system);
    }
    catch(...)
    {
        this->operator=(old_state);
        throw;
    }

    return run_system;
}

/** Return the energy component used by these moves. An exception
    will be raised if the component moves use different energy
    components to one another
    
    \throw SireError::incompatible_error
*/
const Symbol& WeightedMoves::energyComponent() const
{
    int nmoves = mvs.count();
    
    const tuple<MovePtr,double> *mvs_array = mvs.constData();
    
    if (nmoves == 0)
        return Move::null().energyComponent();
    
    else
    {
        const Symbol &nrg = mvs_array[0].get<0>().read().energyComponent();
        
        QSet<Symbol> symbols;
        
        for (int i=1; i<nmoves; ++i)
        {
            if (nrg != mvs_array[i].get<0>().read().energyComponent())
            {
                symbols.insert(nrg);
                symbols.insert(mvs_array[i].get<0>().read().energyComponent());
                
                for (int j=i+1; j<nmoves; ++j)
                {
                    symbols.insert( mvs_array[j].get<0>().read().energyComponent() );
                }
                
                throw SireError::incompatible_error( QObject::tr(
                        "Different moves are sampling from different energy "
                        "components (%1) so a single energy component cannot "
                        "be found for %2.")
                            .arg( Sire::toString(symbols), this->toString() ),
                                CODELOC );
            }
        }
        
        return nrg;
    }
}

/** Return the space property used by these moves. An exception
    will be raised if the component moves use different space
    properties to one another
    
    \throw SireError::incompatible_error
*/
const PropertyName& WeightedMoves::spaceProperty() const
{
    if (not combined_space.isNull())
        //we are using a combined space property
        return combined_space;

    int nmoves = mvs.count();
    
    const tuple<MovePtr,double> *mvs_array = mvs.constData();
    
    if (nmoves == 0)
        return Move::null().spaceProperty();
    
    else
    {
        const PropertyName &space = mvs_array[0].get<0>().read().spaceProperty();
        
        QList<PropertyName> spaces;
        
        for (int i=1; i<nmoves; ++i)
        {
            if (space != mvs_array[i].get<0>().read().spaceProperty())
            {
                spaces.append(space);
                spaces.append(mvs_array[i].get<0>().read().spaceProperty());
                
                for (int j=i+1; j<nmoves; ++j)
                {
                    const PropertyName &jspace = mvs_array[j].get<0>()
                                                             .read().spaceProperty();
                
                    if (not spaces.contains(jspace))
                        spaces.append(jspace);
                }
                
                throw SireError::incompatible_error( QObject::tr(
                        "Different moves are sampling from different space "
                        "properties (%1) so a single space property cannot "
                        "be found for %2.")
                            .arg( Sire::toString(spaces), this->toString() ),
                                CODELOC );
            }
        }
        
        return space;
    }
}

/** Set the energy component of all of the moves to 'component' */
void WeightedMoves::setEnergyComponent(const Symbol &component)
{
    int nmoves = mvs.count();
    
    if (nmoves > 0)
    {
        tuple<MovePtr,double> *mvs_array = mvs.data();
        
        for (int i=0; i<nmoves; ++i)
        {
            if (mvs_array[i].get<0>().read().energyComponent() != component)
            {
                mvs_array[i].get<0>().edit().setEnergyComponent(component);
            }
        }
    }
}

/** Set the name of the property that all of the moves will use to 
    find the simulation space (simulation box) to 'spaceproperty' */
void WeightedMoves::setSpaceProperty(const PropertyName &spaceproperty)
{
    if (spaceproperty == combined_space)
        //nothing needs to be done
        return;

    int nmoves = mvs.count();
    
    if (nmoves > 0)
    {
        tuple<MovePtr,double> *mvs_array = mvs.data();
        
        for (int i=0; i<nmoves; ++i)
        {
            if (mvs_array[i].get<0>().read().spaceProperty() != spaceproperty)
            {
                mvs_array[i].get<0>().edit().setSpaceProperty(spaceproperty);
            }
        }
    }
    
    combined_space = PropertyName();
}

/** Set the combined space property - this tells this moves object
    to return a different space that represents the combined space
    of all of the sub-moves. Note that this does not change the
    space used in the sub-moves */
void WeightedMoves::setCombinedSpaceProperty(const PropertyName &space)
{
    combined_space = space;
}

/** Return whether or not these moves use a combined space to 
    calculate the volume */
bool WeightedMoves::hasCombinedSpaceProperty() const
{
    return not combined_space.isNull();
}

/** Set the temperature for all moves that have a constant temperature
    to 'temperature'. It has already been checked that these moves
    between them sample at constant temperature */
void WeightedMoves::_pvt_setTemperature(const Temperature &temperature)
{
    int nmoves = mvs.count();
    
    if (nmoves > 0)
    {
        tuple<MovePtr,double> *mvs_array = mvs.data();
        
        for (int i=0; i<nmoves; ++i)
        {
            if (mvs_array[i].get<0>().read().isConstantTemperature())
            {
                if (mvs_array[i].get<0>().read().temperature() != temperature)
                {
                    mvs_array[i].get<0>().edit().setTemperature(temperature);
                }
            }
        }
    }
}

/** Set the pressure for all moves that have a constant pressure
    to 'pressure'. It has already been checked that these moves
    between them sample at constant pressure */
void WeightedMoves::_pvt_setPressure(const Pressure &pressure)
{
    int nmoves = mvs.count();
    
    if (nmoves > 0)
    {
        tuple<MovePtr,double> *mvs_array = mvs.data();
        
        for (int i=0; i<nmoves; ++i)
        {
            if (mvs_array[i].get<0>().read().isConstantPressure())
            {
                if (mvs_array[i].get<0>().read().pressure() != pressure)
                {
                    mvs_array[i].get<0>().edit().setPressure(pressure);
                }
            }
        }
    }
}

/** Set the fugacity for all moves that have a constant fugacity
    to 'fugacity'. It has already been checked that these moves
    between them sample at constant fugacity */
void WeightedMoves::_pvt_setFugacity(const Pressure &fugacity)
{
    int nmoves = mvs.count();
    
    if (nmoves > 0)
    {
        tuple<MovePtr,double> *mvs_array = mvs.data();
        
        for (int i=0; i<nmoves; ++i)
        {
            if (mvs_array[i].get<0>().read().isConstantFugacity())
            {
                if (mvs_array[i].get<0>().read().fugacity() != fugacity)
                {
                    mvs_array[i].get<0>().edit().setFugacity(fugacity);
                }
            }
        }
    }
}

/** Return the moves available in this set */
QList<MovePtr> WeightedMoves::moves() const
{
    int nmoves = mvs.count();
    const tuple<MovePtr,double> *mvs_array = mvs.constData();
    
    QList<MovePtr> moves;
    
    for (int i=0; i<nmoves; ++i)
    {
        moves.append( mvs_array[i].get<0>() );
    }
    
    return moves;
}

/** Return the average time to perform each move */
QList<SireUnits::Dimension::Time> WeightedMoves::timing() const
{
    QList<SireUnits::Dimension::Time> times;
    
    for (int i=0; i<avgtimes.count(); ++i)
    {
        times.append( avgtimes[i].average() * nanosecond );
    }
    
    return times;
}

/** Clear all of the timing information */
void WeightedMoves::clearTiming()
{
    for (int i=0; i<avgtimes.count(); ++i)
    {
        avgtimes[i] = Average();
    }
}

const char* WeightedMoves::typeName()
{
    return QMetaType::typeName( qMetaTypeId<WeightedMoves>() );
}
