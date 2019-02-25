/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2014  Christopher Woods
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

#ifndef SIRESYSTEM_VOLMAPMONITOR_H
#define SIRESYSTEM_VOLMAPMONITOR_H

#include "systemmonitor.h"

#include "SireUnits/dimensions.h"
#include "SireBase/propertymap.h"
#include "SireMol/volumemap.h"

SIRE_BEGIN_HEADER

namespace SireSystem
{
class VolMapMonitor;
}

SIRESYSTEM_EXPORT QDataStream& operator<<(QDataStream&, const SireSystem::VolMapMonitor&);
SIRESYSTEM_EXPORT QDataStream& operator>>(QDataStream&, SireSystem::VolMapMonitor&);

namespace SireSystem
{

using SireMol::MoleculeGroup;
using SireMol::VolumeMap;
using SireBase::PropertyMap;

/** This class create an occupation volume map showing the 
    regions of space that are occupied by the monitored atoms
    on average during the simulation
    
    @author Christopher Woods
*/
class SIRESYSTEM_EXPORT VolMapMonitor
        : public SireBase::ConcreteProperty<VolMapMonitor,SystemMonitor>
{

friend SIRESYSTEM_EXPORT QDataStream& ::operator<<(QDataStream&, const VolMapMonitor&);
friend SIRESYSTEM_EXPORT QDataStream& ::operator>>(QDataStream&, VolMapMonitor&);

public:
    VolMapMonitor();
    
    VolMapMonitor(const SireUnits::Dimension::Length &grid_spacing);
    
    VolMapMonitor(const MoleculeGroup &group,
                  const PropertyMap &map = PropertyMap());
    
    VolMapMonitor(const MoleculeGroup &group,
                  const SireUnits::Dimension::Length &grid_spacing,
                  const PropertyMap &map = PropertyMap());
    
    VolMapMonitor(const MoleculeGroup &group,
                  bool skip_light_atoms,
                  const PropertyMap &map = PropertyMap());

    VolMapMonitor(const MoleculeGroup &group,
                  const SireUnits::Dimension::Length &grid_spacing,
                  bool skip_light_atoms,
                  const PropertyMap &map = PropertyMap());
    
    VolMapMonitor(const VolMapMonitor &other);
    
    ~VolMapMonitor();
    
    VolMapMonitor& operator=(const VolMapMonitor &other);
    
    bool operator==(const VolMapMonitor &other) const;
    bool operator!=(const VolMapMonitor &other) const;
    
    static const char* typeName();
    
    const char* what() const;
    
    QString toString() const;
    
    const MoleculeGroup &group() const;
    
    bool skippingLightAtoms() const;
    
    void setSkipLightAtoms(bool on);
    
    void setGridSpacing(const SireUnits::Dimension::Length &grid_spacing);
    
    void setGroup(const MoleculeGroup &group, const PropertyMap &map=PropertyMap());
    
    SireVol::GridInfo gridInfo() const;
    
    SireUnits::Dimension::Length gridSpacing() const;
    
    PropertyMap map() const;
    
    void setPropertyMap(const PropertyMap &map);
    
    QVector<float> averageOccupancy() const;
    
    qint64 nSamples() const;
    
    void clearStatistics();
    
    VolumeMap volumeMap() const;
    
    void monitor(System &system);

private:
    /** The molecule group being monitored */
    SireMol::MolGroupPtr molgroup;
    
    /** The PropertyMap used to find the right properties of the atoms */
    SireBase::PropertyMap propmap;
    
    /** The volume map under construction */
    VolumeMap volmap;
};

}

Q_DECLARE_METATYPE( SireSystem::VolMapMonitor )

SIRE_EXPOSE_CLASS( SireSystem::VolMapMonitor )

SIRE_END_HEADER

#endif
