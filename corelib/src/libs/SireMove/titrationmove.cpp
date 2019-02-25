/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2013  Christopher Woods
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

#include "titrationmove.h"
#include "titrator.h"

#include "SireUnits/units.h"
#include "SireUnits/temperature.h"

#include "SireSystem/system.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>

using namespace SireMove;
using namespace SireSystem;
using namespace SireBase;
using namespace SireUnits;
using namespace SireStream;

static RegisterMetaType<TitrationMove> r_titrationmove;

QDataStream &operator<<(QDataStream &ds, const TitrationMove &move)
{
    writeHeader(ds, r_titrationmove, 1);
    
    ds << static_cast<const MonteCarlo&>(move);
    
    return ds;
}

QDataStream &operator>>(QDataStream &ds, TitrationMove &move)
{
    VersionID v = readHeader(ds, r_titrationmove);
    
    if (v == 1)
    {
        ds >> static_cast<MonteCarlo&>(move);
    }
    else
        throw version_error(v, "1", r_titrationmove, CODELOC);
    
    return ds;
}

/** Null constructor */
TitrationMove::TitrationMove() : ConcreteProperty<TitrationMove,MonteCarlo>()
{
    MonteCarlo::setEnsemble( Ensemble::NVT(25*celsius) );
}

/** Copy constructor */
TitrationMove::TitrationMove(const TitrationMove &other)
              : ConcreteProperty<TitrationMove,MonteCarlo>(other)
{}

/** Destructor */
TitrationMove::~TitrationMove()
{}

void TitrationMove::_pvt_setTemperature(const Temperature &temperature)
{
    MonteCarlo::setEnsemble( Ensemble::NVT(temperature) );
}

const char* TitrationMove::typeName()
{
    return QMetaType::typeName( qMetaTypeId<TitrationMove>() );
}

const char* TitrationMove::what() const
{
    return TitrationMove::typeName();
}

/** Copy assignment operator */
TitrationMove& TitrationMove::operator=(const TitrationMove &other)
{
    if (this != &other)
    {
        MonteCarlo::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool TitrationMove::operator==(const TitrationMove &other) const
{
    return MonteCarlo::operator==(other);
}

/** Comparison operator */
bool TitrationMove::operator!=(const TitrationMove &other) const
{
    return not TitrationMove::operator==(other);
}

/** Return a string representation of the move */
QString TitrationMove::toString() const
{
    return QObject::tr("TitrationMove( nAccept() == %1, nReject() == %2 )")
                .arg(this->nAccepted()).arg(this->nRejected());
}

/** Actually perform the move */
void TitrationMove::move(System &system, int nmoves, bool record_stats)
{
    if (nmoves <= 0)
        return;
 
    //save our, and the system's, current state
    TitrationMove old_state(*this);
    System old_system_state(system);
    
    try
    {
        const PropertyMap &map = Move::propertyMap();

        //first, get hold of the titration property of the system. This provides
        //all of the information about the titratable molecules, and the current
        //state of all of the titratable molecules
        Titrator titrator = system.property( map["titrator"] ).asA<Titrator>();
        
        if (titrator.nIons() == 0 and titrator.nNeutrals() == 0)
            //there is nothing available to move
            return;
        
        for (int i=0; i<nmoves; ++i)
        {
            //get the old total energy of the system
            double old_nrg = system.energy( this->energyComponent() );

            //save the old system
            System old_system(system);

            double old_bias = 1;
            double new_bias = 1;

            //now choose two random groups from the set of titratable groups
            int ion_index = titrator.getIonIndex( generator().randInt(0, titrator.nIons()-1) );
            int neutral_index = titrator.getNeutralIndex(
                                            generator().randInt(0, titrator.nNeutrals()-1) );

            //swap the charges of these two groups
            titrator.swapCharge(neutral_index, ion_index);
    
            //apply the titrator to the system. This will update
            //the charges of the titratable groups and will update
            //the copy of titrator stored in the system. The return value
            //is the change in the 'zero' energy caused by this change. This
            //is used to account for changes in internal state caused by
            //a titration, e.g. the intraresidue change from an aspartate
            //to an aspartic acid is not accounted for only by the change
            //in intramolecular energy from an MM forcefield
            double delta_zero = titrator.applyTo(system);
    
            //calculate the energy of the system
            double new_nrg = system.energy( this->energyComponent() );

            //accept or reject the move based on the change of energy
            //and the biasing factors
            if (not this->test(new_nrg + delta_zero, old_nrg, new_bias, old_bias))
            {
                //the move has been rejected - reset the state
                system = old_system;
                titrator = system.property( map["titrator"] ).asA<Titrator>();
            }

            if (record_stats)
            {
                system.collectStats();
            }
        }
    }
    catch(...)
    {
        this->operator=(old_state);
        system = old_system_state;
        throw;
    }
}
