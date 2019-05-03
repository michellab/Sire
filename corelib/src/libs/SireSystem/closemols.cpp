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

#include "closemols.h"
#include "system.h"

#include "SireMol/molecule.h"
#include "SireMol/moleculedata.h"
#include "SireMol/atomcoords.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>

using namespace SireSystem;
using namespace SireFF;
using namespace SireMol;
using namespace SireVol;
using namespace SireStream;

static const RegisterMetaType<CloseMols> r_closemols(NO_ROOT);

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const CloseMols &closemols)
{
    writeHeader(ds, r_closemols, 1);
    
    SharedDataStream sds(ds);
    
    sds << closemols.p << closemols.molgroup << closemols.spce
        << closemols.nclosest << closemols.map;
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, CloseMols &closemols)
{
    VersionID v = readHeader(ds, r_closemols);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        PointPtr p;
        MolGroupPtr molgroup;
        SpacePtr spce;
        quint32 nclosest;
        PropertyMap map;
        
        sds >> p >> molgroup >> spce >> nclosest >> map;
            
        closemols = CloseMols(p, molgroup, spce, nclosest, map);
    }
    else
        throw version_error(v, "1", r_closemols, CODELOC);
        
    return ds;
}

/** Internal function used to work out which of the close molecules
    is now furthest from the point */
void CloseMols::getNewFurthestMolNum()
{
    double cutoff_dist = 0;
    furthest_molnum = MolNum();
    
    for (QHash<MolNum,double>::const_iterator it = close_mols.constBegin();
         it != close_mols.constEnd();
         ++it)
    {
        if (it.value() > cutoff_dist)
        {
            furthest_molnum = it.key();
            cutoff_dist = it.value();
        }
        else if (it.value() == cutoff_dist and it.key() > furthest_molnum)
        {
            furthest_molnum = it.key();
        }
    }
    
    cutoff_dist2 = cutoff_dist * cutoff_dist;
}

/** Internal function used to return whether or not the passed collection
    of molecules all belong to the close molecules of this CloseMols object */
bool CloseMols::differentMolecules(const QHash<MolNum,double> &molecules) const
{
    if (molecules.count() == close_mols.count())
    {
        for (QHash<MolNum,double>::const_iterator it = close_mols.constBegin();
             it != close_mols.constEnd();
             ++it)
        {
            if (not molecules.contains(it.key()))
            {
                return true;
            }
        }
        
        return false;
    }
    else
        return true;
}

/** Internal function that rescans through the molecules to find 
    the closest ones to the point - this returns whether or not this
    changes the identity of the close molecules */
bool CloseMols::recalculate()
{
    QHash<MolNum,double> old_close_mols = close_mols;

    //clear the current list
    close_mols = QHash<MolNum,double>();
    furthest_molnum = MolNum();
    cutoff_dist2 = 0;

    if (nclosest <= 0)
        return false;

    const Molecules &molecules = molgroup.read().molecules();
    
    const quint32 nmols = molecules.nMolecules();

    const Vector &point = p.read().point();
    
    const Space &space = spce.read();
    
    const PropertyName &coords_property = map["coordinates"];
    
    if (nmols <= nclosest)
    {
        //we are selecting all of the molecules
        close_mols.reserve(nmols);
        
        for (Molecules::const_iterator it = molecules.constBegin();
             it != molecules.constEnd();
             ++it)
        {
            //just get the center of the whole molecule
            const Vector &center = it->data().property(coords_property)
                                             .asA<AtomCoords>()
                                             .array().aaBox().center();
                                             
            const double dist = space.calcDist(point, center);
            close_mols.insert(it.key(), dist);
        }
        
        cutoff_dist2 = std::numeric_limits<double>::max();
        
        return this->differentMolecules(old_close_mols);
    }

    close_mols.reserve(nclosest);

    quint32 nfound = 0;
    
    for (Molecules::const_iterator it = molecules.constBegin();
         it != molecules.constEnd();
         ++it)
    {
        //just get the center of the whole molecule
        const Vector &center = it->data().property(coords_property)
                                         .asA<AtomCoords>()
                                         .array().aaBox().center();
        
        const double dist2 = space.calcDist2(point, center);
        
        if (nfound < nclosest)
        {
            //there is still space for this molecule
            ++nfound;
            
            close_mols.insert(it.key(), std::sqrt(dist2));
            
            if (dist2 > cutoff_dist2)
            {
                cutoff_dist2 = dist2;
                furthest_molnum = it.key();
            }
            else if (dist2 == cutoff_dist2 and it.key() > furthest_molnum)
            {
                furthest_molnum = it.key();
            }
        }
        else
        {
            //there isn't space - is this molecule closer than one
            //of the existing molecules?
            if (dist2 < cutoff_dist2 or
                (dist2 == cutoff_dist2 and it.key() < furthest_molnum))
            {
                //it is - replace the current furthest molecule with this molecule
                close_mols.remove(furthest_molnum);
                close_mols.insert(it.key(), std::sqrt(dist2));
                
                this->getNewFurthestMolNum();
            }
        }
    }
    
    return this->differentMolecules(old_close_mols);
}

/** Internal function that recalculates the distance of the molecule
    with number 'molnum' and updates the list of close molecules. This
    returns whether or not this changes the identity of the close
    molecules */
bool CloseMols::recalculate(MolNum changed_mol)
{
    const Molecules &molecules = molgroup.read().molecules();

    Molecules::const_iterator it = molecules.constFind(changed_mol);
    
    if (it == molecules.constEnd())
        //this molecule is not contained
        return false;

    const Vector &point = p.read().point();
    const Space &space = spce.read();

    const PropertyName &coords_property = map["coordinates"];

    //calculate the distance from the new molecule to the point
    const Vector &center = it->data().property(coords_property)
                                     .asA<AtomCoords>()
                                     .array().aaBox().center();
    
    const double dist2 = space.calcDist2(point, center);
    
    //is this already a close molecule?
    if (close_mols.contains(changed_mol))
    {
        if (dist2 > cutoff_dist2 or
            (dist2 == cutoff_dist2 and changed_mol > furthest_molnum))
            //the molecule has moved beyond the cutoff, so other
            //molecules in the group may now be closer - we now need
            //to rescan all molecules
            return this->recalculate();
        else
        {
            //the molecule has moved, but this cannot affect the order
            //(as it hasn't gone past the cutoff)
            close_mols[changed_mol] = std::sqrt(dist2);
            
            if (dist2 == cutoff_dist2 or changed_mol == furthest_molnum)
            {
                //the identity of the assigned molecule furthest from 
                //the point may have changed
                this->getNewFurthestMolNum();
            }
            
            return false;
        }
    }
    else
    {
        //has this molecule replaced any of the existing close molecules?
        if ( dist2 < cutoff_dist2 or
             (dist2 == cutoff_dist2 and changed_mol < furthest_molnum) )
        {
            //yes it has!
            close_mols.remove(furthest_molnum);
            close_mols.insert(changed_mol, std::sqrt(dist2));
            this->getNewFurthestMolNum();
            return true;
        }
        else
            //no, it hasn't!
            return false;
    }
}

/** Internal function that recalculates the distances of the molecules
    in 'changed_mols' and updates the list of close molecules. This
    returns whether or not this changes the identity of the close
    molecules */
bool CloseMols::recalculate(const Molecules &changed_mols)
{
    if (changed_mols.isEmpty())
        return this->recalculate();
    
    else if (changed_mols.nMolecules() == 1)
        return this->recalculate( changed_mols.constBegin().key() );

    const Vector &point = p.read().point();
    const Space &space = spce.read();

    const PropertyName &coords_property = map["coordinates"];

    const Molecules &molecules = molgroup.read().molecules();

    bool changed_order = false;

    QHash<MolNum,double> old_close_mols = close_mols;

    for (Molecules::const_iterator it = changed_mols.constBegin();
         it != changed_mols.constEnd();
         ++it)
    {
        const MolNum changed_mol = it.key();
        
        if (not molecules.contains(changed_mol))
            continue;
            
        //calculate the distance from the new molecule to the point
        const Vector &center = molecules.constFind(changed_mol)
                                        ->data().property(coords_property)
                                                .asA<AtomCoords>()
                                                .array().aaBox().center();

        const double dist2 = space.calcDist2(point, center);
    
        //is this already a close molecule?
        if (close_mols.contains(changed_mol))
        {
            if (dist2 > cutoff_dist2 or
                (dist2 == cutoff_dist2 and changed_mol > furthest_molnum))
                //the molecule has moved beyond the cutoff, so other
                //molecules in the group may now be closer - we now need
                //to rescan all molecules
                return this->recalculate();
            else
            {
                //the molecule has moved, but this cannot affect the order
                //(as it hasn't gone past the cutoff)
                close_mols[changed_mol] = std::sqrt(dist2);
                
                if (dist2 == cutoff_dist2 or changed_mol == furthest_molnum)
                {
                    //the identity of the assigned molecule furthest from 
                    //the point may have changed
                    this->getNewFurthestMolNum();
                }
            }
        }
        else
        {
            //has this molecule replaced any of the existing close molecules?
            if ( dist2 < cutoff_dist2 or
                 (dist2 == cutoff_dist2 and changed_mol < furthest_molnum) )
            {
                //yes it has!
                close_mols.remove(furthest_molnum);
                close_mols.insert(changed_mol, std::sqrt(dist2));
                this->getNewFurthestMolNum();
                changed_order = true;
            }
        }
    }
    
    if (changed_order)
        return this->differentMolecules(old_close_mols);
    else
        return false;
}

/** Constructor */
CloseMols::CloseMols() : nclosest(0), cutoff_dist2(0)
{}

/** Construct to find the 'nclosest' molecules from the molecule group
    'molgroup' to the point 'point' */
CloseMols::CloseMols(const PointRef &point, const MoleculeGroup &mgroup,
                     int nclose, const PropertyMap &propmap)
          : p(point), molgroup(mgroup), nclosest(nclose), 
            map(propmap), cutoff_dist2(0)
{
    if (nclose < 0)
        nclosest = 0;
    
    if (p.read().usesMoleculesIn(mgroup))
        p.edit().update(mgroup);
    
    if (nclosest > 0)
        this->recalculate();
}

/** Construct to find the 'nclosest' molecules from the molecule group
    'molgroup' to the point 'point', using the space 'space' to calculate
    the distances between the molecules and the point */
CloseMols::CloseMols(const PointRef &point, const MoleculeGroup &mgroup,
                     const Space &space, int nclose, const PropertyMap &propmap)
          : p(point), molgroup(mgroup), spce(space), nclosest(nclose), 
            map(propmap), cutoff_dist2(0)
{
    if (nclose < 0)
        nclosest = 0;
    
    if (p.read().usesMoleculesIn(mgroup))
        p.edit().update(mgroup);
    
    if (nclosest > 0)
        this->recalculate();
}

/** Copy constructor */
CloseMols::CloseMols(const CloseMols &other)
          : p(other.p), molgroup(other.molgroup), spce(other.spce),
            nclosest(other.nclosest), map(other.map), 
            close_mols(other.close_mols),
            furthest_molnum(other.furthest_molnum),
            cutoff_dist2(other.cutoff_dist2)
{}

/** Destructor */
CloseMols::~CloseMols()
{}

/** Copy assignment operator */
CloseMols& CloseMols::operator=(const CloseMols &other)
{
    if (this != &other)
    {
        p = other.p;
        molgroup = other.molgroup;
        spce = other.spce;
        nclosest = other.nclosest;
        map = other.map;
        close_mols = other.close_mols;
        furthest_molnum = other.furthest_molnum;
        cutoff_dist2 = other.cutoff_dist2;
    }
    
    return *this;
}

/** Comparison operator */
bool CloseMols::operator==(const CloseMols &other) const
{
    return (this == &other) or
           (p == other.p and molgroup == other.molgroup and
            spce == other.spce and nclosest == other.nclosest and
            map == other.map);
}

/** Comparison operator */
bool CloseMols::operator!=(const CloseMols &other) const
{
    return not this->operator==(other);
}

const char* CloseMols::typeName()
{
    return QMetaType::typeName( qMetaTypeId<CloseMols>() );
}

const char* CloseMols::what() const
{
    return CloseMols::typeName();
}

CloseMols* CloseMols::clone() const
{
    return new CloseMols(*this);
}

/** Return the point used to find the closest molecules */
const Point& CloseMols::point() const
{
    return p.read();
}

/** Return the molecule group that contains the molecules */
const MoleculeGroup& CloseMols::moleculeGroup() const
{
    return molgroup.read();
}

/** Return the space used to calculate the distances between the 
    molecules and the point */
const Space& CloseMols::space() const
{
    return spce.read();
}

/** Return the property map used to find the coordinates and
    space properties */
const PropertyMap& CloseMols::propertyMap() const
{
    return map;
}

/** Return the number of molecules to record */
int CloseMols::nClosest() const
{
    return nclosest;
}

/** Return the set of close molecules, together with the 
    distances from the molecule to the point */
const QHash<MolNum,double>& CloseMols::closeMolecules() const
{
    return close_mols;
}

/** Internal function used to update the data for this object from
    the passed system - this returns (via arguments) whether or not
    the location of the point has changed (point_changed), whether
    or not the space has changed (space_changed), and whether the
    major and/or minor version numbers of the molecule group
    has changed (molgroup_major_change and molgroup_minor_change) */
void CloseMols::updateData(const System &system,
                           bool &point_changed, bool &space_changed,
                           bool &molgroup_major_changed,
                           bool &molgroup_minor_changed)
{
    point_changed = false;
    space_changed = false;
    molgroup_major_changed = false;
    molgroup_minor_changed = false;
    
    if (p.read().usesMoleculesIn(system))
        point_changed = p.edit().update(system);

    const Space &new_space = system.property( map["space"] )
                                   .asA<Space>();
                                   
    if (not spce.read().equals(new_space))
    {
        spce = new_space;
        space_changed = true;
    }
    
    if (system.contains(molgroup.read().number()))
    {
        const MoleculeGroup &newgroup = system[molgroup.read().number()];
        
        molgroup_major_changed = newgroup.version().majorVersion() != 
                                 molgroup.read().version().majorVersion();
                                 
        if (molgroup_major_changed)
            molgroup_minor_changed = true;
        else
        {
            molgroup_minor_changed = newgroup.version().minorVersion() !=
                                     molgroup.read().version().minorVersion();
        }
        
        if (molgroup_major_changed or molgroup_minor_changed)
        {
            molgroup = newgroup;
        }
    }
    else
    {
        //ensure that the molecule group has the same version of molecules
        //as in the passed system (so it is compatible with the molecules
        //in the point)
        QList<Molecule> changed_mols = molgroup.edit().update(system.molecules());
        
        molgroup_minor_changed = (not changed_mols.isEmpty());
    }
}

/** Update from the passed system - this returns whether or not this
    update changes the identity of the close molecules */
bool CloseMols::update(const System &system)
{
    CloseMols old_state(*this);
    
    try
    {
        bool point_changed(false), space_changed(false);
        bool molgroup_major_changed(false), molgroup_minor_changed(false);
        
        this->updateData(system, point_changed, space_changed,
                         molgroup_major_changed, molgroup_minor_changed);
    
        if (point_changed or space_changed or 
            molgroup_major_changed or molgroup_minor_changed)
        {
            return this->recalculate();
        }
    }
    catch(...)
    {
        this->operator=(old_state);
        throw;
    }
    
    return false;
}

/** Update from the passed system, providing the hint that only the
    molecule with number 'molnum' has changed since the last update.
    It is a bad idea to provide an incorrect hint! This returns whether
    or not this changes the identity of the close molecules */
bool CloseMols::update(const System &system, MolNum changed_mol)
{
    CloseMols old_state(*this);
    
    try
    {
        bool point_changed(false), space_changed(false);
        bool molgroup_major_changed(false), molgroup_minor_changed(false);
        
        this->updateData(system, point_changed, space_changed,
                         molgroup_major_changed, molgroup_minor_changed);
    
        if (point_changed or space_changed or molgroup_major_changed)
        {
            //more than just the molecule has changed!
            return this->recalculate();
        }
        else if (molgroup_minor_changed)
        {
            //we will trust the user that only the molecule with 
            //number 'changed_mol' has changed
            return this->recalculate(changed_mol);
        }
    }
    catch(...)
    {
        this->operator=(old_state);
        throw;
    }
    
    return false;
}

/** Update from the passed system, providing the hint that only the
    molecules in 'changed_mols' have changed since the last update.
    It is a bad idea to provide an incorrect hint! This returns whether
    or not this changes the identity of the close molecules */
bool CloseMols::update(const System &system, const Molecules &changed_mols)
{
    CloseMols old_state(*this);
    
    try
    {
        bool point_changed(false), space_changed(false);
        bool molgroup_major_changed(false), molgroup_minor_changed(false);
        
        this->updateData(system, point_changed, space_changed,
                         molgroup_major_changed, molgroup_minor_changed);
    
        if (point_changed or space_changed or molgroup_major_changed)
        {
            //more than just the molecules have changed!
            return this->recalculate();
        }
        else if (molgroup_minor_changed)
        {
            //we will trust the user that only the molecules in 
            //'changed_mols' have changed
            return this->recalculate(changed_mols);
        }
    }
    catch(...)
    {
        this->operator=(old_state);
        throw;
    }
    
    return false;
}
