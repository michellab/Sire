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

#ifndef SIREMOVE_VOLUMEMOVE_H
#define SIREMOVE_VOLUMEMOVE_H

#include "montecarlo.h"
#include "volumechanger.h"

SIRE_BEGIN_HEADER

namespace SireMove
{
class VolumeMove;
}

SIREMOVE_EXPORT QDataStream& operator<<(QDataStream&, const SireMove::VolumeMove&);
SIREMOVE_EXPORT QDataStream& operator>>(QDataStream&, SireMove::VolumeMove&);

namespace SireMove
{

using SireMol::MGID;
using SireMol::MoleculeGroup;

/** This is a Monte Carlo volume move. This is used to allow
    the pressure to be kept constant
    
    @author Christopher Woods
*/
class SIREMOVE_EXPORT VolumeMove
            : public SireBase::ConcreteProperty<VolumeMove,MonteCarlo>
{

friend SIREMOVE_EXPORT QDataStream& ::operator<<(QDataStream&, const VolumeMove&);
friend SIREMOVE_EXPORT QDataStream& ::operator>>(QDataStream&, VolumeMove&);

public:
    VolumeMove(const PropertyMap &map = PropertyMap());

    VolumeMove(const MGID &mgid, const PropertyMap &map = PropertyMap());
    VolumeMove(const MoleculeGroup &molgroup,
               const PropertyMap &map = PropertyMap());
    
    VolumeMove(const VolumeChanger &volchanger,
               const PropertyMap &map = PropertyMap());
    
    VolumeMove(const VolumeMove &other);
    
    ~VolumeMove();
    
    VolumeMove& operator=(const VolumeMove &other);
    
    static const char* typeName();
        
    bool operator==(const VolumeMove &other) const;
    bool operator!=(const VolumeMove &other) const;
    
    QString toString() const;

    void setVolumeChanger(const VolumeChanger &volchanger);
    void setVolumeChanger(const MoleculeGroup &molgroup);
    
    const VolumeChanger& volumeChanger() const;

    const MGID& groupID() const;
    
    void setGenerator(const RanGenerator &rangenerator);
    
    void setMaximumVolumeChange(const SireUnits::Dimension::Volume &delta);
    const SireUnits::Dimension::Volume& maximumVolumeChange() const;
    
    void move(System &system, int nmoves, bool record_stats=true);

protected:
    void _pvt_setTemperature(const SireUnits::Dimension::Temperature &temperature);
    void _pvt_setPressure(const SireUnits::Dimension::Pressure &pressure);
    
private:
    /** The volume changing function used to change the volume of  
        the system */
    VolumeChangerPtr volchanger;

    #ifndef SKIP_BROKEN_GCCXML_PARTS
    /** The maximum volume change */
    SireUnits::Dimension::Volume maxchange;
    #endif
};

}

Q_DECLARE_METATYPE( SireMove::VolumeMove )

SIRE_EXPOSE_CLASS( SireMove::VolumeMove )

SIRE_END_HEADER

#endif
