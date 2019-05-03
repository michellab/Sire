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

#include "volumemove.h"

#include "SireSystem/system.h"
#include "SireVol/space.h"

#include "SireMol/molecule.h"
#include "SireMol/mover.hpp"
#include "SireMol/moleculegroup.h"
#include "SireMol/mgname.h"

#include "SireUnits/dimensions.h"
#include "SireUnits/units.h"
#include "SireUnits/temperature.h"

#include "SireBase/savestate.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireError/errors.h"

using namespace SireMove;
using namespace SireMol;
using namespace SireVol;
using namespace SireSystem;
using namespace SireBase;
using namespace SireUnits;
using namespace SireUnits::Dimension;
using namespace SireStream;

static const RegisterMetaType<VolumeMove> r_volmove;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const VolumeMove &volmove)
{
    writeHeader(ds, r_volmove, 2);
    
    ds << volmove.volchanger
       << double(volmove.maxchange.to(angstrom3)) 
       << static_cast<const MonteCarlo&>(volmove);
       
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, VolumeMove &volmove)
{
    VersionID v = readHeader(ds, r_volmove);
    
    if (v == 2)
    {
        double maxchange;
        ds >> volmove.volchanger >> maxchange
           >> static_cast<MonteCarlo&>(volmove);
           
        volmove.maxchange = maxchange * angstrom3;
    }
    else if (v == 1)
    {
        double maxchange;
        ds >> maxchange >> static_cast<MonteCarlo&>(volmove);
        
        volmove.maxchange = maxchange * angstrom3;
        volmove.volchanger = ScaleVolumeFromCenter( MGIdentifier() );
    }
    else
        throw version_error( v, "1", r_volmove, CODELOC );
        
    return ds;
}

/** Null constructor */
VolumeMove::VolumeMove(const PropertyMap &map)
           : ConcreteProperty<VolumeMove,MonteCarlo>(map), maxchange(0)
{
    MonteCarlo::setEnsemble( Ensemble::NPT( 25 * celsius, 1 * atm ) );
}

/** Construct a volume move that can be used to generate the ensemble
    for a temperature of 25 C, pressure of 1 atm, and with a maximum 
    change of 100 A^3 by moving the molecules in the
    molecule groups that match the ID 'mgid' 
    using a ScaleVolumeFromCenter centered on the origin */
VolumeMove::VolumeMove(const MGID &mgid, const PropertyMap &map)
           : ConcreteProperty<VolumeMove,MonteCarlo>(map),
             volchanger( ScaleVolumeFromCenter(mgid) ),
             maxchange(100*angstrom3)
{
    MonteCarlo::setEnsemble( Ensemble::NPT( 25 * celsius, 1 * atm ) );
}

/** Construct a volume move that can be used to generate the ensemble
    for a temperature of 25 C, pressure of 1 atm, and with a maximum 
    change of 100 A^3 by moving the molecules in 'molgroup' 
    using a ScaleVolumeFromCenter centered on the origin */
VolumeMove::VolumeMove(const MoleculeGroup &molgroup,
                       const PropertyMap &map)
           : ConcreteProperty<VolumeMove,MonteCarlo>(map),
             volchanger( ScaleVolumeFromCenter(molgroup) ),
             maxchange(100*angstrom3)
{
    MonteCarlo::setEnsemble( Ensemble::NPT( 25 * celsius, 1 * atm ) );
}

/** Construct a volume move that can be used to generate the ensemble
    for a temperature of 25 C, pressure of 1 atm, and with a maximum 
    change of 100 A^3 using the passed volume changer */
VolumeMove::VolumeMove(const VolumeChanger &volumechanger,
                       const PropertyMap &map)
           : ConcreteProperty<VolumeMove,MonteCarlo>(map),
             volchanger(volumechanger),
             maxchange(100*angstrom3)
{
    MonteCarlo::setEnsemble( Ensemble::NPT( 25 * celsius, 1 * atm ) );
}

/** Copy constructor */
VolumeMove::VolumeMove(const VolumeMove &other)
           : ConcreteProperty<VolumeMove,MonteCarlo>(other),
             volchanger(other.volchanger),
             maxchange(other.maxchange)
{}

/** Destructor */
VolumeMove::~VolumeMove()
{}

/** Copy assignment operator */
VolumeMove& VolumeMove::operator=(const VolumeMove &other)
{
    volchanger = other.volchanger;
    maxchange = other.maxchange;
    MonteCarlo::operator=(other);
    
    return *this;
}

/** Comparison operator */
bool VolumeMove::operator==(const VolumeMove &other) const
{
    return volchanger == other.volchanger and
           maxchange == other.maxchange and
           MonteCarlo::operator==(other);
}

/** Comparison operator */
bool VolumeMove::operator!=(const VolumeMove &other) const
{
    return volchanger != other.volchanger or
           maxchange != other.maxchange or
           MonteCarlo::operator!=(other);
}

/** Return a string representation of this move */
QString VolumeMove::toString() const
{
    return QObject::tr("VolumeMove( maximumVolumeChange() = %1 A^3 "
                       "nAccepted() = %2 nRejected() = %3 )")
                  .arg(this->maximumVolumeChange().to(angstrom3))
                  .arg(this->nAccepted())
                  .arg(this->nRejected());
}

/** Internal function called to set the temperature */
void VolumeMove::_pvt_setTemperature(const Temperature &temperature)
{
    MonteCarlo::setEnsemble( Ensemble::NPT( temperature, this->pressure() ) );
}

/** Internal function called to set the pressure */
void VolumeMove::_pvt_setPressure(const Pressure &pressure)
{
    MonteCarlo::setEnsemble( Ensemble::NPT( this->temperature(), pressure ) );
}

/** Set the maximum change in volume */
void VolumeMove::setMaximumVolumeChange(const Volume &delta)
{
    maxchange = delta;
}

/** Return the maximum change of volume attempted by a move */
const Volume& VolumeMove::maximumVolumeChange() const
{
    return maxchange;
}

/** Set the volume changer used to change the volume to 'volchanger' */
void VolumeMove::setVolumeChanger(const VolumeChanger &new_volchanger)
{
    volchanger = new_volchanger;
}

/** Set the volume changer used to change the volume to a 
    ScaleVolumeFromCenter that scales the molecules in 'molgroup'
    from the center of a box centered at (0,0,0) */
void VolumeMove::setVolumeChanger(const MoleculeGroup &molgroup)
{
    volchanger = ScaleVolumeFromCenter(molgroup, Vector(0));
}

/** Return the volume changer used to change the volume */
const VolumeChanger& VolumeMove::volumeChanger() const
{
    return volchanger.read();
}

/** Return the ID that matches the molecule groups that
    will be affected by this move */
const MGID& VolumeMove::groupID() const
{
    return volchanger.read().groupID();
}

/** Set the random number generator used by this move */
void VolumeMove::setGenerator(const RanGenerator &rangenerator)
{
    MonteCarlo::setGenerator(rangenerator);
    volchanger.edit().setGenerator(this->generator());
}

/** Perform 'nmoves' volume moves on the passed system, optionally
    recording simulation statistics if 'record_stats' is true */
void VolumeMove::move(System &system, int nmoves, bool record_stats)
{
    if (nmoves <= 0)
        return;

    SaveState old_system_state = SaveState::save(system);

    VolumeMove old_state(*this);
    
    try
    {
        const PropertyMap &map = this->propertyMap();

        for (int i=0; i<nmoves; ++i)
        {
            System old_system(system);
        
            //calculate the old energy and volume
            double old_nrg = this->energy(system);
            Volume old_vol = this->volume(system);
            
            double old_bias = 1;
            double new_bias = 1;
            
            int nmols = this->volumeChanger()
                            .randomChangeVolume(system, maxchange, 
                                                new_bias, old_bias, map);
            
            //calculate the new energy and volume
            double new_nrg = this->energy(system);
            Volume new_vol = this->volume(system);
            
            if (not this->test(new_nrg, old_nrg, nmols, 
                               new_vol, old_vol,
                               new_bias, old_bias))
            {
                //move failed - go back to the last step
                system = old_system;
            }
            
            if (record_stats)
            {
                system.collectStats();
            }
        }
    }
    catch(...)
    {
        old_system_state.restore(system);
        this->operator=(old_state);
        throw;
    }
    
}

const char* VolumeMove::typeName()
{
    return QMetaType::typeName( qMetaTypeId<VolumeMove>() );
}
