/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2009  Christopher Woods
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

#ifndef SIRESYSTEM_CLOSEMOLS_H
#define SIRESYSTEM_CLOSEMOLS_H

#include <QSet>

#include "SireFF/point.h"

#include "SireBase/propertymap.h"

#include "SireMol/moleculegroup.h"
#include "SireVol/space.h"

SIRE_BEGIN_HEADER

namespace SireSystem
{
class CloseMols;
}

QDataStream& operator<<(QDataStream&, const SireSystem::CloseMols&);
QDataStream& operator>>(QDataStream&, SireSystem::CloseMols&);

namespace SireSystem
{

class System;

using SireFF::Point;

using SireMol::MoleculeGroup;
using SireMol::MolNum;
using SireMol::Molecules;

using SireVol::Space;

using SireBase::PropertyMap;

/** This class is used to maintain a list of the closest molecules
    to a specified point
    
    @author Christopher Woods
*/
class SIRESYSTEM_EXPORT CloseMols
{

friend QDataStream& ::operator<<(QDataStream&, const CloseMols&);
friend QDataStream& ::operator>>(QDataStream&, CloseMols&);

public:
    CloseMols();

    CloseMols(const SireFF::PointRef &point, const MoleculeGroup &molgroup,
              int nclosest=1, const PropertyMap &map = PropertyMap());

    CloseMols(const SireFF::PointRef &point, const MoleculeGroup &molgroup,
              const Space &space, int nclosest=1, const PropertyMap &map = PropertyMap());
    
    CloseMols(const CloseMols &other);
    
    ~CloseMols();
    
    CloseMols& operator=(const CloseMols &other);
    
    bool operator==(const CloseMols &other) const;
    bool operator!=(const CloseMols &other) const;
    
    static const char* typeName();
    
    const char* what() const;
    
    CloseMols* clone() const;
    
    const Point& point() const;
    
    const MoleculeGroup& moleculeGroup() const;
    
    const Space& space() const;
    
    const PropertyMap& propertyMap() const;
    
    int nClosest() const;
    
    const QHash<MolNum,double>& closeMolecules() const;
    
    bool isClose(MolNum molnum) const;
    
    bool update(const System &system);
    bool update(const System &system, MolNum changed_mol);
    bool update(const System &system, const Molecules &molecules);

private:
    bool recalculate();
    bool recalculate(MolNum changed_mol);
    bool recalculate(const Molecules &molecules);

    void getNewFurthestMolNum();

    bool differentMolecules(const QHash<MolNum,double> &mols) const;
    
    void updateData(const System &system,
                    bool &point_changed, bool &space_changed,
                    bool &molgroup_major_changed,
                    bool &molgroup_minor_changed);
    
    /** The point in space from which to find the closest molecules */
    SireFF::PointPtr p;
    
    /** The molecule group containing the molecules */
    SireMol::MolGroupPtr molgroup;
    
    /** The space used to calculate the distances between the 
        molecules and the point */
    SireVol::SpacePtr spce;
    
    /** The number of molecules to record */
    quint32 nclosest;
    
    /** The property map containing the location of the coordinates 
        and space properties */
    PropertyMap map;
    
    /** The numbers of the 'nclosest' closest molecules, 
        together with the distance from the molecule to the point */
    QHash<MolNum,double> close_mols;
    
    /** The molecule number of the furthest recorded molecule */
    MolNum furthest_molnum;
    
    /** The distance^2 from the point to the furthest recorded molecule */
    double cutoff_dist2;
};

/** Return whether or not the molecule with number 'molnum' is
    one of the close molecules */
inline bool CloseMols::isClose(MolNum molnum) const
{
    return close_mols.contains(molnum);
}

}

Q_DECLARE_METATYPE( SireSystem::CloseMols )

SIRE_EXPOSE_CLASS( SireSystem::CloseMols )

SIRE_END_HEADER

#endif
