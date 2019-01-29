/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2012  Christopher Woods
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

#ifndef SIRESYSTEM_IDASSIGNER_H
#define SIRESYSTEM_IDASSIGNER_H

#include <QVector>

#include "closemols.h"

#include "SireMol/moleculegroup.h"
#include "SireMol/partialmolecule.h"
#include "SireMaths/nvector.h"
#include "SireVol/space.h"
#include "SireBase/propertymap.h"

SIRE_BEGIN_HEADER

namespace SireSystem
{
class IDAssigner;
}

SIRESYSTEM_EXPORT QDataStream& operator<<(QDataStream&, const SireSystem::IDAssigner&);
SIRESYSTEM_EXPORT QDataStream& operator>>(QDataStream&, SireSystem::IDAssigner&);

namespace SireSystem
{

using SireBase::PropertyMap;
using SireMol::MoleculeGroup;
using SireMol::PartialMolecule;
using SireFF::PointPtr;
using SireFF::PointRef;

/** This class uses the machinery of the identity point to 
    pick out molecules that are associated with identity points.
    This is useful if you want to monitor a property or energy,
    but don't actually want to change the coordinates of atoms
    
    @author Christopher Woods
*/
class SIRESYSTEM_EXPORT IDAssigner 
         : public SireBase::ConcreteProperty<IDAssigner,SireBase::Property>
{

friend QDataStream& ::operator<<(QDataStream&, const IDAssigner&);
friend QDataStream& ::operator>>(QDataStream&, IDAssigner&);

public:
    IDAssigner();
    
    IDAssigner(const PointRef &point,
               const MoleculeGroup &molgroup,
               const PropertyMap &map = PropertyMap());
    
    IDAssigner(const QVector<SireFF::PointPtr> &points,
               const MoleculeGroup &molgroup,
               const PropertyMap &map = PropertyMap());
    
    IDAssigner(const PointRef &point,
               const MoleculeGroup &molgroup,
               const SireVol::Space &space,
               const PropertyMap &map = PropertyMap());
    
    IDAssigner(const QVector<SireFF::PointPtr> &points,
               const MoleculeGroup &molgroup,
               const SireVol::Space &space,
               const PropertyMap &map = PropertyMap());
    
    IDAssigner(const IDAssigner &other);

    ~IDAssigner();
    
    static const char* typeName();
    
    const char* what() const;
    
    IDAssigner* clone() const;
    
    IDAssigner& operator=(const IDAssigner &other);
    
    bool operator==(const IDAssigner &other) const;
    bool operator!=(const IDAssigner &other) const;
    
    QString toString() const;
    
    const MoleculeGroup& moleculeGroup() const;
    
    QVector<SireFF::PointPtr> points() const;

    int nPoints() const;

    const PropertyMap& propertyMap() const;

    QVector<PartialMolecule> identifiedMolecules() const;

    void update(const System &system);

    const Space& space() const;

private:
    void validatePoints(const QVector<SireFF::PointPtr> &points) const;
    void validateGroup(const MoleculeGroup &new_group) const;

    bool updatePoints(const System &system);
    bool updateGroup(const System &system);
    bool updateSpace(const System &system);

    void rebuildMolToMolNum();

    void recalculateDistances();

    void assignMoleculesToPoints();

    /** The molecule group that will be scanned for molecules */
    SireMol::MolGroupPtr molgroup;
    
    /** The set of identity points */
    QVector<PointPtr> identity_points;
 
    /** The space being used to calculate the distances */
    SireVol::SpacePtr spce;

    /** The property map used to find the properties used
        by this constraint */
    PropertyMap map;

    /** The collection of the 'npoints' closest molecules to each
        of the identity points */
    QVector<CloseMols> points_with_mols;
    
    /** The mapping of index to molecule number - this is necessary
        as the mapping above is just over the close molecules - we
        need to record the numbers of the close molecules */
    QVector<SireMol::MolNum> mol_to_molnum;
    
    /** The distances between all molecules in 'mol_to_molnum' and 
        all of the identity points */
    QHash<SireMol::MolNum,SireMaths::NVector> point_distances;
    
    /** The current mapping of points to molecules. This is the index
        of the molecule number of the matching molecule for each point
        in the "mol_to_molnum" array, e.g. point_to_mol[3] will return
        the molecule matched to the fourth point, with point_to_mol[3]
        being the index in mol_to_molnum of the molecule number of that
        molecule */
    QVector<int> point_to_mol;
    
    /** Whether or not the distances have changed since the last 
        time this constraint was updated - if they have, then the
        mapping of molecules to points needs to be recalculated */
    bool distances_changed;
};

}

Q_DECLARE_METATYPE( SireSystem::IDAssigner )

SIRE_EXPOSE_CLASS( SireSystem::IDAssigner )

SIRE_END_HEADER

#endif
