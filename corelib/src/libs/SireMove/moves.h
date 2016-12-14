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

#ifndef SIREMOVE_MOVES_H
#define SIREMOVE_MOVES_H

#include <QList>

#include "move.h"

#include "SireUnits/dimensions.h"
#include "SireMaths/accumulator.h"

SIRE_BEGIN_HEADER

namespace SireMove
{
class Moves;
class SameMoves;
}

QDataStream &operator<<(QDataStream&, const SireMove::Moves&);
QDataStream &operator>>(QDataStream&, SireMove::Moves&);

QDataStream &operator<<(QDataStream&, const SireMove::SameMoves&);
QDataStream &operator>>(QDataStream&, SireMove::SameMoves&);

namespace SireMove
{

using SireCAS::Symbol;

/** This is the base class of all Moves objects. These are objects
    that contain a collection of moves that are applied to a system
    
    @author Christopher Woods
*/
class SIREMOVE_EXPORT Moves : public SireBase::Property
{

friend QDataStream& ::operator<<(QDataStream&, const Moves&);
friend QDataStream& ::operator>>(QDataStream&, Moves&);

public:
    Moves();
    
    Moves(const Moves &other);
    
    virtual ~Moves();
    
    virtual Moves* clone() const=0;
    
    static const char* typeName()
    {
        return "SireMove::Moves";
    }
    
    const Move& operator[](int i) const;
    
    int count() const;
    int size() const;
    int nMoveTypes() const;
    
    virtual int nMoves() const;
    
    virtual QString toString() const=0;
    
    virtual const Symbol& energyComponent() const=0;
    virtual const PropertyName& spaceProperty() const=0;
    
    virtual void setEnergyComponent(const Symbol &nrg_component);
    virtual void setSpaceProperty(const PropertyName &space_property);
    
    SireUnits::Dimension::MolarEnergy energy(System &system) const;
    SireUnits::Dimension::Volume volume(const System &system) const;
    
    Ensemble ensemble() const;
    
    bool isConstantEnergy() const;
    bool isConstantTemperature() const;
    bool isConstantVolume() const;
    bool isConstantPressure() const;
    bool isConstantChemicalPotential() const;
    bool isConstantFugacity() const;

    bool isConstantLambda(const Symbol &lam) const;
    
    SireUnits::Dimension::Temperature temperature() const;
    SireUnits::Dimension::Pressure pressure() const;
    SireUnits::Dimension::Pressure fugacity() const;
    SireUnits::Dimension::MolarEnergy chemicalPotential() const;

    virtual void setGenerator(const RanGenerator &rangenerator)=0;

    void setTemperature(const SireUnits::Dimension::Temperature &temperature);
    void setPressure(const SireUnits::Dimension::Pressure &pressure);
    void setChemicalPotential(
                    const SireUnits::Dimension::MolarEnergy &chemical_potential);
    void setFugacity(const SireUnits::Dimension::Pressure &fugacity);
        
    virtual System move(const System &system, 
                        int nmoves=1, bool record_stats=false)=0;
    
    virtual void clearStatistics()=0;

    void setCheckRunningTotal(bool on);
    void setAcceptableDelta(SireUnits::Dimension::MolarEnergy delta);
    
    SireUnits::Dimension::MolarEnergy acceptableDelta() const;
    
    bool checkingRunningTotal() const;
    
    virtual QList<MovePtr> moves() const=0;

    SireUnits::Dimension::Time timing(int i) const;

    virtual QList<SireUnits::Dimension::Time> timing() const=0;
    virtual void clearTiming()=0;

    static const SameMoves& null();

protected:
    void preCheck(System &system) const;
    void postCheck(System &system) const;

    Moves& operator=(const Moves &other);
    
    bool operator==(const Moves &other) const;

    /** Set the temperature for all moves that have a constant temperature
        to 'temperature'. It has already been checked that these moves
        between them sample at constant temperature */
    virtual void _pvt_setTemperature(
                            const SireUnits::Dimension::Temperature &temperature)=0;
                            
    /** Set the pressure for all moves that have a constant pressure
        to 'pressure'. It has already been checked that these moves
        between them sample at constant pressure */
    virtual void _pvt_setPressure(
                            const SireUnits::Dimension::Pressure &pressure)=0;
                            
    /** Set the fugacity for all moves that have a constant fugacity
        to 'fugacity'. It has already been checked that these moves
        between them sample at constant fugacity */
    virtual void _pvt_setFugacity(
                            const SireUnits::Dimension::Pressure &fugacity)=0;

private:
    /** The acceptable delta between the running total and recalcualted total energy */
    double acceptable_delta;

    /** Whether or not to check the running total against the recalculated
        total energy in the post-check of the moves */
    bool check_running_total;
};

/** This is a Moves class that just applies the same move over and
    over again
    
    @author Christopher Woods
*/
class SIREMOVE_EXPORT SameMoves 
         : public SireBase::ConcreteProperty<SameMoves,Moves>
{

friend QDataStream& ::operator<<(QDataStream&, const SameMoves&);
friend QDataStream& ::operator>>(QDataStream&, SameMoves&);

public:
    SameMoves();
    SameMoves(const Move &move);
    
    SameMoves(const SameMoves &other);
    
    ~SameMoves();
    
    static const char* typeName();
    
    SameMoves& operator=(const SameMoves &other);
    
    bool operator==(const SameMoves &other) const;
    bool operator!=(const SameMoves &other) const;
    
    SameMoves* clone() const;
    
    QString toString() const;
    
    using Moves::move;
    
    void setEnergyComponent(const Symbol &component);
    void setSpaceProperty(const PropertyName &spaceproperty); 

    const Symbol& energyComponent() const;
    const PropertyName& spaceProperty() const;
    
    void setGenerator(const RanGenerator &rangenerator);
    
    System move(const System &system, int nmoves=1, bool record_stats=false);
    
    void clearStatistics();
    
    QList<MovePtr> moves() const;

    QList<SireUnits::Dimension::Time> timing() const;
    void clearTiming();

private:
    void _pvt_setTemperature(const SireUnits::Dimension::Temperature &temperature);
    void _pvt_setPressure(const SireUnits::Dimension::Pressure &pressure);
    void _pvt_setFugacity(const SireUnits::Dimension::Pressure &fugacity);

    bool sampleSameSpaceProperty() const;

    /** The move that will be repeatedly applied */
    MovePtr mv;
    
    /** The average time to run each move, in nanoseconds */
    SireMaths::Average avgtime;
};

typedef SireBase::PropPtr<Moves> MovesPtr;

}

Q_DECLARE_METATYPE( SireMove::SameMoves )

SIRE_EXPOSE_CLASS( SireMove::Moves )
SIRE_EXPOSE_CLASS( SireMove::SameMoves )

SIRE_EXPOSE_PROPERTY( SireMove::MovesPtr, SireMove::Moves )

SIRE_END_HEADER

#endif
