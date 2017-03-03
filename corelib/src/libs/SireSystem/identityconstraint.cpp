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

#include <QVarLengthArray>

#include "identityconstraint.h"
#include "system.h"
#include "closemols.h"
#include "delta.h"

#include "SireSystem/errors.h"

#include "SireMaths/linearap.h"
#include "SireMaths/nmatrix.h"
#include "SireMaths/nvector.h"

#include "SireMol/molecules.h"
#include "SireMol/molecule.h"
#include "SireMol/moleculegroup.h"
#include "SireMol/moleditor.h"
#include "SireMol/viewsofmol.h"

#include "SireBase/refcountdata.h"
#include "SireBase/shareddatapointer.hpp"

#include "SireVol/space.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QTime>
#include <QDebug>

using namespace SireSystem;
using namespace SireFF;
using namespace SireMol;
using namespace SireBase;
using namespace SireMaths;
using namespace SireVol;
using namespace SireStream;

namespace SireSystem
{
namespace detail
{

class NullIdentityConstraintHelper;

/** This is the virtual base class of the IdentityConstraint helper
    classes - these allow the IdentityConstraint to use different
    strategies to solve the constraint depending on the ratio
    of the number of points to the number of molecules */
class IdentityConstraintPvt : public RefCountData
{
public:
    IdentityConstraintPvt();
    
    IdentityConstraintPvt(const MoleculeGroup &molgroup,
                          const PointRef &point,
                          const PropertyMap &map);
    
    IdentityConstraintPvt(const MoleculeGroup &molgroup,
                          const QVector<PointPtr> &points,
                          const PropertyMap &map);
    
    IdentityConstraintPvt(const IdentityConstraintPvt &other);
    
    virtual ~IdentityConstraintPvt();
    
    static const char* typeName()
    {
        return "SireSystem::detail::IdentityConstraintPvt";
    }
    
    virtual const char* what() const=0;
    
    virtual IdentityConstraintPvt* clone() const=0;
    
    virtual Molecules update(const System &system, bool new_system)=0;
    virtual Molecules update(const System &system, MolNum changed_mol,
                             bool new_system)=0;
    virtual Molecules update(const System &system, const Molecules &molecules,
                             bool new_system)=0;
    Molecules update(const System &system, const QList<MolNum> &changed_mols,
                     bool new_system);

    virtual Molecules applyConstraint() const=0;

    const MolGroupPtr& molGroupPtr() const;

    const MoleculeGroup& moleculeGroup() const;
    
    const QVector<PointPtr>& points() const;
    
    const PropertyMap& propertyMap() const;

    const Space& space() const;

    static const NullIdentityConstraintHelper& null();

    virtual bool isNull() const
    {
        return false;
    }

protected:
    bool updatePoints(const System &system);
    bool updateGroup(const System &system);
    bool updateSpace(const System &system);

    void validateGroup(const MoleculeGroup &new_group) const;

    void validatePoints(const QVector<PointPtr> &points) const;

    /** The molecule group containing the molecules whose identities
        are being constrained */
    MolGroupPtr molgroup;
    
    /** The set of points that provide the locations that
        identify the molecules. These 'n' points constrain
        the identity of the first 'n' molecule in 'molgroup'
        (with point 'i' constraining the identity of 
        molecule 'i') */
    QVector<PointPtr> identity_points;

    /** The space being used to calculate the distances */
    SpacePtr spce;

    /** The property map used to find the properties used
        by this constraint */
    PropertyMap map;
};

/** This is the null identity constraint helper used when
    there are no points against which the constraint is applied
    
    @author Christopher Woods
*/
class NullIdentityConstraintHelper : public IdentityConstraintPvt
{
public:
    NullIdentityConstraintHelper();
    
    NullIdentityConstraintHelper(const NullIdentityConstraintHelper &other);
    
    ~NullIdentityConstraintHelper();
    
    static const char* typeName()
    {
        return "SireSystem::detail::NullIdentityConstraintHelper";
    }
    
    const char* what() const;
    
    NullIdentityConstraintHelper* clone() const;

    Molecules update(const System &system, bool new_system);
    Molecules update(const System &system, MolNum changed_mol, bool new_system);
    Molecules update(const System &system, const Molecules &molecules, bool new_system);
    Molecules update(const System &system, const QList<MolNum> &changed_mols, 
                     bool new_system);

    Molecules applyConstraint() const;
    
    bool isNull() const
    {
        return true;
    }
};

/** This is the IdentityConstraint helper class that is used
    when the number of points is comparable to the number
    of molecules (or if the points are widely dispersed)
    
    In these cases it is best just to calculate all of the 
    distances to all of the points all of the time
    
    @author Christopher Woods
*/
class ManyPointsHelper : public IdentityConstraintPvt
{
public:
    ManyPointsHelper();
    
    ManyPointsHelper(const MoleculeGroup &molgroup,
                     const PropertyMap &map);
    
    ManyPointsHelper(const MoleculeGroup &molgroup,
                     const QVector<PointPtr> &points,
                     const PropertyMap &map);
                               
    ManyPointsHelper(const ManyPointsHelper &other);
    
    ~ManyPointsHelper();
    
    static const char* typeName()
    {
        return "SireSystem::detail::ManyPointsHelper";
    }
    
    const char* what() const;
    
    ManyPointsHelper* clone() const;

    Molecules update(const System &system, bool new_system);
    Molecules update(const System &system, MolNum changed_mol, bool new_system);
    Molecules update(const System &system, const Molecules &molecules, bool new_system);
    Molecules update(const System &system, const QList<MolNum> &changed_mols, 
                     bool new_system);

    Molecules applyConstraint() const;

private:
    void recalculateDistances();
    void recalculateDistances(MolNum molnum);
    void recalculateDistances(const Molecules &molecules);

    void assignMoleculesToPoints();

    /** The (squared) distances between all molecules and all points */
    QHash<MolNum,NVector> point_distances;
    
    /** The current mapping of molecules to points - if this is
        empty than the molecules are correctly mapped so that the
        nth molecule is assigned to the nth point */
    QVector<int> mol_to_point;
    
    /** Whether or not the set of distances have changed since the
        last time the constraint was applied - if they have, then
        the order of molecules needs to be recalculated */
    bool distances_changed;
};

/** This is the identity constraint helper used when the number
    of points is significantly less than the number of molecules,
    and the points are close together. In this case, the 
    constraint can be simplified by only calculating distances
    between molecules which are close to the points */
class FewPointsHelper : public IdentityConstraintPvt
{
public:
    FewPointsHelper();
    FewPointsHelper(const MoleculeGroup &molgroup,
                    const QVector<PointPtr> &points,
                    const PropertyMap &map);
                               
    FewPointsHelper(const FewPointsHelper &other);
    
    ~FewPointsHelper();
    
    static const char* typeName()
    {
        return "SireSystem::detail::FewPointsHelper";
    }
    
    const char* what() const;
    
    FewPointsHelper* clone() const;

    Molecules update(const System &system, bool new_system);
    Molecules update(const System &system, MolNum changed_mol, bool new_system);
    Molecules update(const System &system, const Molecules &molecules, bool new_system);
    Molecules update(const System &system, const QList<MolNum> &changed_mols, 
                     bool new_system);
    
    Molecules applyConstraint() const;

private:
    void rebuildMolToMolNum();

    void recalculateDistances();
    void recalculateDistances(MolNum changed_mol);
    void recalculateDistances(const Molecules &changed_mols);

    void assignMoleculesToPoints();

    /** The collection of the 'npoints' closest molecules to each
        of the identity points */
    QVector<CloseMols> points_with_mols;
    
    /** The mapping of index to molecule number - this is necessary
        as the mapping above is just over the close molecules - we
        need to record the numbers of the close molecules */
    QVector<MolNum> mol_to_molnum;
    
    /** The distances between all molecules in 'mol_to_molnum' and 
        all of the identity points */
    QHash<MolNum,NVector> point_distances;
    
    /** The current mapping of molecules to points - if this is empty
        then the molecules are correctly mapped (i.e. the nth molecule
        is mapped to the nth point) */
    QVector<int> mol_to_point;
    
    /** Whether or not the distances have changed since the last 
        time this constraint was updated - if they have, then the
        mapping of molecules to points needs to be recalculated */
    bool distances_changed;
};

/** This is the identity constraint helper class used when there
    is only a single point against which the molecule is constrained.
    This is the easiest case, as the constraint just has to look
    for the closest molecule */
class SinglePointHelper : public IdentityConstraintPvt
{
public:
    SinglePointHelper();
    
    SinglePointHelper(const MoleculeGroup &molgroup,
                      const PointRef &point,
                      const PropertyMap &map);
    
    SinglePointHelper(const MoleculeGroup &molgroup,
                      const QVector<PointPtr> &points,
                      const PropertyMap &map);
                               
    SinglePointHelper(const SinglePointHelper &other);
    
    ~SinglePointHelper();
    
    static const char* typeName()
    {
        return "SireSystem::detail::SinglePointHelper";
    }
    
    const char* what() const;
    
    SinglePointHelper* clone() const;

    Molecules update(const System &system, bool new_system);
    Molecules update(const System &system, MolNum changed_mol, bool new_system);
    Molecules update(const System &system, const Molecules &molecules, bool new_system);
    Molecules update(const System &system, const QList<MolNum> &changed_mols, 
                     bool new_system);
    
    Molecules applyConstraint() const;

private:
    void recalculateDistances();
    void recalculateDistances(MolNum changed_mol);
    void recalculateDistances(const Molecules &molecules);

    /** The number of the current molecule closest to the point */
    MolNum closest_molnum;
    
    /** The distance^2 between the closest molecule and the point */
    double closest_distance2;
};


} // end of namespace detail
} // end of namespace SireSystem

using namespace SireSystem::detail;

/////////
///////// Implementation of IdentityConstraintPvt
/////////

void IdentityConstraintPvt::validatePoints(const QVector<PointPtr> &points) const
{
    for (QVector<PointPtr>::const_iterator it = points.constBegin();
         it != points.constEnd();
         ++it)
    {
        if (it->isNull())
            throw SireError::incompatible_error( QObject::tr(
                    "You cannot create an identity constraint with "
                    "null identity points!"), CODELOC );
                    
        else if ( not (it->read().isIntraMoleculePoint() or
                       it->read().isExtraMoleculePoint()) )
        {
            throw SireError::incompatible_error( QObject::tr(
                    "The identity constraint can only be used with "
                    "intramolecular or extramolecular identity points."),
                        CODELOC );
        }
    }
}

void IdentityConstraintPvt::validateGroup(const MoleculeGroup &new_group) const
{
    //run through each molecule and ensure that the coordinate properties
    //of them all are compatible
    
    if (new_group.nMolecules() < 2)
        return;
    
    const PropertyName &coords_property = map["coordinates"];
    
    const ViewsOfMol &first_mol = new_group.moleculeAt(0);
    
    const AtomCoords &coords = first_mol.data().property(coords_property)
                                               .asA<AtomCoords>();
    
    for ( Molecules::const_iterator it = new_group.constBegin(); 
          it != new_group.constEnd(); 
          ++it )
    {
        if (not coords.isCompatibleWith(it->data().info()))
            throw SireError::incompatible_error( QObject::tr(
                    "The coordinates property (%1) of molecule number %2 "
                    "is not compatible with that of molecule number %3. "
                    "Only molecules with completely compatible coordinates "
                    "properties can be used in an IdentityConstraint.")
                        .arg(coords_property.toString())
                        .arg(first_mol.number())
                        .arg(it.key()), CODELOC );
    }
}

/** Constructor */
IdentityConstraintPvt::IdentityConstraintPvt() : RefCountData()
{}

/** Constructor */
IdentityConstraintPvt::IdentityConstraintPvt(const MoleculeGroup &group,
                                             const QVector<PointPtr> &points,
                                             const PropertyMap &property_map)
                      : RefCountData(),
                        molgroup(group), identity_points(points), map(property_map)
{
    this->validateGroup(group);
    this->validatePoints(points);
    identity_points.squeeze();
}

/** Constructor */
IdentityConstraintPvt::IdentityConstraintPvt(const MoleculeGroup &group,
                                             const PointRef &point,
                                             const PropertyMap &property_map)
                      : RefCountData(),
                        molgroup(group), map(property_map)
{
    this->validateGroup(group);
    identity_points.append( PointPtr(point) );
    this->validatePoints(identity_points);
    identity_points.squeeze();
}

/** Copy constructor */
IdentityConstraintPvt::IdentityConstraintPvt(const IdentityConstraintPvt &other)
                      : RefCountData(),
                        molgroup(other.molgroup), 
                        identity_points(other.identity_points), 
                        spce(other.spce), map(other.map)
{}

/** Destructor */
IdentityConstraintPvt::~IdentityConstraintPvt()
{}

/** Return the shared pointer to the molecule group */
const MolGroupPtr& IdentityConstraintPvt::molGroupPtr() const
{
    return molgroup;
}

Q_GLOBAL_STATIC( MoleculeGroup, nullMoleculeGroup );

/** Return the molecule group operated on by this constraint */
const MoleculeGroup& IdentityConstraintPvt::moleculeGroup() const
{
    if (molgroup.isNull())
        return *(nullMoleculeGroup());
    else
        return molgroup.read();
}

/** Return the points used to identify the molecules */
const QVector<PointPtr>& IdentityConstraintPvt::points() const
{
    return identity_points;
}

/** Return the space used to calculate distances between the molecules
    and the identity points */
const Space& IdentityConstraintPvt::space() const
{
    return spce.read();
}

/** Return the property map used to find the properties required
    by this constraint */
const PropertyMap& IdentityConstraintPvt::propertyMap() const
{
    return map;
}

/** Update the space used to calculate the distances between
    the points and the molecules - this returns whether or 
    not the space has changed */
bool IdentityConstraintPvt::updateSpace(const System &system)
{
    const Space &new_space = system.property(map["space"]).asA<Space>();
    
    if (spce != new_space)
    {
        spce = new_space;
        return true;
    }
    else
        return false;
}

/** Update the points in this constraint from the passed system - 
    this returns whether or not this changes any points */
bool IdentityConstraintPvt::updatePoints(const System &system)
{
    int npoints = identity_points.count();
    
    const PointPtr *const_points_array = identity_points.constData();
    
    bool need_update = false;
    
    for (int i=0; i<npoints; ++i)
    {
        if (const_points_array[i].read().usesMoleculesIn(system))
        {
            need_update = true;
            break;
        }
    }
    
    if (need_update)
    {
        bool changed = false;
        
        PointPtr *points_array = identity_points.data();
        
        for (int i=0; i<npoints; ++i)
        {
            bool this_changed = points_array[i].edit().update(system);
            
            changed = changed or this_changed;
        }
        
        return changed;
    }
    else
        return false;
}

/** Update the molecule group whose molecules are constrained to equal
    the version from the passed system. Note that all of the molecules
    in the molecule group *must* have exactly the same atomic
    layout, as the coordinates are just straight swapped */
bool IdentityConstraintPvt::updateGroup(const System &system)
{
    const MoleculeGroup &old_group = molgroup.read();

    if (not system.contains(old_group.number()))
    {
        molgroup.edit().update(system.molecules());
        return true;
    }
        
    const MoleculeGroup &new_group = system[old_group.number()];
    
    if (new_group.version() == old_group.version())
        return false;
        
    else if (new_group.version().majorVersion() != old_group.version().majorVersion())
    {
        //the group's contents have changed - ensure that all
        //of the molecules are still compatible
        this->validateGroup(new_group);
    }
    
    molgroup = new_group;
    return true;
}

Molecules IdentityConstraintPvt::update(const System &system, 
                                        const QList<MolNum> &changed_mols,
                                        bool new_system)
{
    if (new_system)
        return this->update(system, true);

    else if (changed_mols.isEmpty())
        return Molecules();
    
    else if (changed_mols.count() == 1)
        return this->update(system, changed_mols.at(0), false);
    
    else
    {
        Molecules mols;
        
        foreach (MolNum molnum, changed_mols)
        {
            mols.add( system[molnum].molecule() );
        }
        
        return this->update(system, mols, false);
    }
}

/////////
///////// Implementation of NullIdentityConstraintHelper
/////////

/** Constructor */
NullIdentityConstraintHelper::NullIdentityConstraintHelper()
                             : IdentityConstraintPvt()
{}

/** Copy constructor */
NullIdentityConstraintHelper::NullIdentityConstraintHelper(
                                        const NullIdentityConstraintHelper&)
                             : IdentityConstraintPvt()
{}

/** Destructor */
NullIdentityConstraintHelper::~NullIdentityConstraintHelper()
{}

const char* NullIdentityConstraintHelper::what() const
{
    return NullIdentityConstraintHelper::typeName();
}

NullIdentityConstraintHelper* NullIdentityConstraintHelper::clone() const
{
    return new NullIdentityConstraintHelper(*this);
}

Molecules NullIdentityConstraintHelper::update(const System&, bool)
{
    return Molecules();
}

Molecules NullIdentityConstraintHelper::update(const System&, MolNum, bool)
{
    return Molecules();
}

Molecules NullIdentityConstraintHelper::update(const System&, const Molecules&, bool)
{
    return Molecules();
}

Molecules NullIdentityConstraintHelper::update(const System&, const QList<MolNum>&, bool)
{
    return Molecules();
}

Molecules NullIdentityConstraintHelper::applyConstraint() const
{
    return Molecules();
}

Q_GLOBAL_STATIC( NullIdentityConstraintHelper, nullHelper );

/** Return the global null helper */
const NullIdentityConstraintHelper& IdentityConstraintPvt::null()
{
    return *(nullHelper());
}

/////////
///////// Implementation of SinglePointHelper
/////////

/** Constructor */
SinglePointHelper::SinglePointHelper() 
                  : IdentityConstraintPvt(),
                    closest_distance2( std::numeric_limits<double>::max() )
{}

/** Constructor */
SinglePointHelper::SinglePointHelper(const MoleculeGroup &molgroup,
                                     const PointRef &point,
                                     const PropertyMap &map)
                  : IdentityConstraintPvt(molgroup, point, map),
                    closest_distance2( std::numeric_limits<double>::max() )
{}

/** Constructor */
SinglePointHelper::SinglePointHelper(const MoleculeGroup &molgroup,
                                     const QVector<PointPtr> &points,
                                     const PropertyMap &map)
                  : IdentityConstraintPvt(molgroup, points, map),
                    closest_distance2( std::numeric_limits<double>::max() )
{
    if (points.count() > 1)
        throw SireError::program_bug( QObject::tr(
                "It is wrong to try and create a SinglePointHelper using "
                "%1 points!!!").arg(points.count()), CODELOC );
}

/** Copy constructor */
SinglePointHelper::SinglePointHelper(const SinglePointHelper &other)
                  : IdentityConstraintPvt(other),
                    closest_molnum(other.closest_molnum),
                    closest_distance2(other.closest_distance2)
{}

/** Destructor */
SinglePointHelper::~SinglePointHelper()
{}

const char* SinglePointHelper::what() const
{
    return SinglePointHelper::typeName();
}

SinglePointHelper* SinglePointHelper::clone() const
{
    return new SinglePointHelper(*this);
}

/** Recalulate all of the distances from scratch to find the closest
    molecule to the identity constraint point */
void SinglePointHelper::recalculateDistances()
{
    if (identity_points.isEmpty())
        return;

    //loop through every molecule in the group and calculate
    //its distance to the constraint point
    closest_distance2 = std::numeric_limits<double>::max();
    closest_molnum = MolNum();

    const MoleculeGroup &molgroup = this->moleculeGroup();
    
    const Vector &point = identity_points.at(0).read().point();
    
    const PropertyName &coords_property = map["coordinates"];
    
    for (Molecules::const_iterator it = molgroup.constBegin();
         it != molgroup.constEnd();
         ++it)
    {
        //get the coordinates property of this molecule
        const AtomCoords &coords = it->data().property(coords_property)
                                             .asA<AtomCoords>();
    
        double new_dist2 = space().calcDist2( point, coords.array().aaBox().center() );
        
        if (new_dist2 < closest_distance2)
        {
            closest_distance2 = new_dist2;
            closest_molnum = it.key();
        }
        else if (new_dist2 == closest_distance2)
        {
            //choose the molecule with the lowest number
            if (it.key() < closest_molnum)
            {
                closest_molnum = it.key();
            }
        }
    }
}

/** Recalculate the distances to find the identity constraint point,
    using the hint that only the molecule with number 'changed_mol' 
    has moved */
void SinglePointHelper::recalculateDistances(MolNum changed_mol)
{
    if (closest_molnum.isNull())
    {
        this->recalculateDistances();
        return;
    }

    if (not this->moleculeGroup().contains(changed_mol))
    {
        this->recalculateDistances();
        return;
    }

    //get the coordinates property of this molecule
    const PropertyName &coords_property = map["coordinates"];
    
    const ViewsOfMol &molecule = this->moleculeGroup()[changed_mol];
    
    const AtomCoords &coords = molecule.data().property(coords_property)
                                              .asA<AtomCoords>();

    const Vector &point = identity_points.at(0).read().point();

    double new_dist2 = space().calcDist2( point, coords.array().aaBox().center() );
    
    if (closest_molnum == changed_mol)
    {
        if (new_dist2 > closest_distance2)
            //we need to recalculate all of the distances
            this->recalculateDistances();
        else
            closest_distance2 = new_dist2;
    }
    else if (new_dist2 < closest_distance2)
    {
        closest_distance2 = new_dist2;
        closest_molnum = changed_mol;
    }
    else if (new_dist2 == closest_distance2)
    {
        //choose the molecule with the lowest molecule number
        if (changed_mol < closest_molnum)
            closest_molnum = changed_mol;
    }
}

/** Recalculate the distances to find the identity constraint point,
    using the hint that only the molecules in 'molecules' have changed */
void SinglePointHelper::recalculateDistances(const Molecules &molecules)
{
    if (closest_molnum.isNull())
    {
        this->recalculateDistances();
        return;
    }

    const PropertyName &coords_property = map["coordinates"];
    const Vector &point = identity_points.at(0).read().point();

    if (molecules.contains(closest_molnum))
    {
        //check the current closest molecule first, as if this
        //has moved away, then we need to recalculate from scratch

        if (not this->moleculeGroup().contains(closest_molnum))
        {
            this->recalculateDistances();
            return;
        }

        const ViewsOfMol &molecule = this->moleculeGroup()[closest_molnum];
    
        const AtomCoords &coords = molecule.data().property(coords_property)
                                                  .asA<AtomCoords>();

        double new_dist2 = space().calcDist2( point, coords.array().aaBox().center() );

        if (new_dist2 > closest_distance2)
        {
            this->recalculateDistances();
            return;
        }
        else
            closest_distance2 = new_dist2;
    }   

    MolNum old_closest_molnum = closest_molnum;

    for (Molecules::const_iterator it = molecules.constBegin();
         it != molecules.constEnd();
         ++it)
    {
        if (it.key() == old_closest_molnum)
            continue;

        if (not this->moleculeGroup().contains(it.key()))
            continue;
    
        const ViewsOfMol &molecule = this->moleculeGroup()[it.key()];
    
        //get the coordinates property of this molecule
        const AtomCoords &coords = molecule.data().property(coords_property)
                                           .asA<AtomCoords>();

        double new_dist2 = space().calcDist2( point, coords.array().aaBox().center() );
    
        if (new_dist2 < closest_distance2)
        {
            closest_distance2 = new_dist2;
            closest_molnum = it.key();
        }
        else if (new_dist2 == closest_distance2)
        {
            //choose the molecule with the lowest molecule number
            if (it.key() < closest_molnum)
                closest_molnum = it.key();
        }
    }
}

/** Actually apply the constraint - if the first molecule is not
    the closest to the point, then this swaps the coordinates of 
    the closest molecule with the first molecule */
Molecules SinglePointHelper::applyConstraint() const
{
    if (this->moleculeGroup().isEmpty())
        return Molecules();
        
    else if (this->moleculeGroup().molNumAt(0) == closest_molnum)
        //this constraint is already satisfied
        return Molecules();
        
    else
    {
        //we need to swap the coordinates of the first molecule
        //in the group with the closest molecule to the constraint
        //point
        Molecules changed_mols;
        
        Molecule first_mol = this->moleculeGroup().moleculeAt(0).molecule();
        Molecule closest_mol = this->moleculeGroup()[closest_molnum].molecule();
        
        const PropertyName &coords_property = map["coordinates"];
        
        AtomCoords first_mol_coords = first_mol.property(coords_property)
                                               .asA<AtomCoords>();
                                               
        AtomCoords closest_mol_coords = closest_mol.property(coords_property)
                                                   .asA<AtomCoords>();
                                                   
        first_mol = first_mol.edit().setProperty(coords_property, closest_mol_coords);
        closest_mol = closest_mol.edit().setProperty(coords_property, first_mol_coords);
        
        changed_mols.add(first_mol);
        changed_mols.add(closest_mol);
        
        return changed_mols;
    }
}

/** Update this constraint, returning what needs to be changed */
Molecules SinglePointHelper::update(const System &system, bool new_system)
{
    if (identity_points.isEmpty())
        return Molecules();
        
    bool new_group = this->updateGroup(system);
    bool new_points = this->updatePoints(system);
    bool new_space = false;
    
    if (new_system)
        new_space = this->updateSpace(system);

    if (new_points or new_group or new_space)
        //recalculate all distances to find the closest molecule
        this->recalculateDistances();

    return this->applyConstraint();
}

/** Update this constraint, returning what needs to be changed */
Molecules SinglePointHelper::update(const System &system, MolNum changed_mol,
                                    bool new_system)
{
    if (identity_points.isEmpty())
        return Molecules();

    bool new_group = this->updateGroup(system);
    bool new_points = this->updatePoints(system);

    bool new_space = false;
    
    if (new_system)
        new_space = this->updateSpace(system);

    if (new_points or new_space)
    {
        //points have changed - need to recalculate all distances
        this->recalculateDistances();
    }
    else if (new_group)
    {
        this->recalculateDistances(changed_mol);
    }

    return this->applyConstraint();
}

/** Update this constraint, returning what needs to be changed */
Molecules SinglePointHelper::update(const System &system, const Molecules &molecules,
                                    bool new_system)
{
    if (identity_points.isEmpty())
        return Molecules();

    bool new_group = this->updateGroup(system);
    bool new_points = this->updatePoints(system);

    bool new_space = false;
    
    if (new_system)
        new_space = this->updateSpace(system);

    if (new_points or new_space)
    {
        //points have changed - need to recalculate all distances
        this->recalculateDistances();
    }
    else if (new_group)
    {
        this->recalculateDistances(molecules);
    }

    return this->applyConstraint();
}

/////////
///////// Implementation of FewPointsHelper
/////////

/** Constructor */
FewPointsHelper::FewPointsHelper() : IdentityConstraintPvt()
{}

/** Construct to constrain the passed molecules to the passed groups - 
    this constrains the ith molecule to the ith point */
FewPointsHelper::FewPointsHelper(const MoleculeGroup &molgroup,
                                 const QVector<PointPtr> &points,
                                 const PropertyMap &map)
                : IdentityConstraintPvt(molgroup, points, map)
{
    if (not identity_points.isEmpty())
    {
        const int npoints = identity_points.count();
        
        //construct a CloseMols to record the closest 'npoints' molecules
        //to each of the identity points
        points_with_mols.reserve(npoints);
        
        for (int i=0; i<npoints; ++i)
        {
            points_with_mols.append( CloseMols(identity_points.at(i), 
                                               molgroup, npoints, map) );
        }
        
        points_with_mols.squeeze();
        
        this->rebuildMolToMolNum();
        this->recalculateDistances();
    }
}
  
/** Copy constructor */                         
FewPointsHelper::FewPointsHelper(const FewPointsHelper &other)
                : IdentityConstraintPvt(other),
                  points_with_mols(other.points_with_mols),
                  mol_to_molnum(other.mol_to_molnum),
                  point_distances(other.point_distances),
                  mol_to_point(other.mol_to_point),
                  distances_changed(other.distances_changed)
{}

/** Destructor */
FewPointsHelper::~FewPointsHelper()
{}

const char* FewPointsHelper::what() const
{
    return FewPointsHelper::typeName();
}

FewPointsHelper* FewPointsHelper::clone() const
{
    return new FewPointsHelper(*this);
}

/** Internal function used to recalculate the distances from all
    of the points to the molecule with number 'molnum' - this does
    nothing if the molecule is not in 'mol_to_molnum' (which is
    implied by point_distances not already containing distances
    for this molecule) */
void FewPointsHelper::recalculateDistances(MolNum molnum)
{
    if (not point_distances.contains(molnum))
        return;
        
    const int npoints = identity_points.count();
    const PointPtr *const_points_array = identity_points.constData();

    const Molecules &molecules = this->moleculeGroup().molecules();

    const PropertyName &coords_property = map["coordinates"];

    Molecules::const_iterator it = molecules.constFind(molnum);
    BOOST_ASSERT( it != molecules.constEnd() );

    const AtomCoords &coords = it->data().property(coords_property)
                                         .asA<AtomCoords>();
                                         
    Vector center = coords.array().aaBox().center();
    
    NVector &distances = point_distances[molnum];

    BOOST_ASSERT(distances.count() == npoints);
    
    double *distances_array = distances.data();
    
    bool any_changes = false;
    
    for (int i=0; i<npoints; ++i)
    {
        double dist2 = space().calcDist2( center, const_points_array[i].read().point() );
        
        if (dist2 != distances_array[i])
        {
            distances_array[i] = dist2;
            any_changes = true;
        }
    }
    
    if (any_changes)
        distances_changed = true;
}

/** Internal function used to recalculate the distances from all
    of the points to all of the molecules in 'changed_mols' - this does
    nothing for the molecules that are not in 'mol_to_molnum' (which is
    implied by point_distances not already containing distances
    for the molecule) */
void FewPointsHelper::recalculateDistances(const Molecules &changed_mols)
{
    if (changed_mols.isEmpty())
        return;

    else if (changed_mols.count() == 1)
    {
        this->recalculateDistances( changed_mols.constBegin().key() );
        return;
    }

    const Molecules &molecules = this->moleculeGroup().molecules();
    
    const int npoints = identity_points.count();
    const PointPtr *const_points_array = identity_points.constData();
    
    QVector< Vector > points(npoints);
    Vector *points_array = points.data();
    
    for (int i=0; i<npoints; ++i)
    {
        points_array[i] = const_points_array[i].read().point();
    }

    const PropertyName &coords_property = map["coordinates"];

    bool any_changes = false;

    for (Molecules::const_iterator it = changed_mols.constBegin();
         it != changed_mols.constEnd();
         ++it)
    {
        const MolNum molnum = it.key();
        
        if (not point_distances.contains(molnum))
            continue;
            
        Molecules::const_iterator it2 = molecules.constFind(molnum);
        BOOST_ASSERT( it2 != molecules.constEnd() );

        const AtomCoords &coords = it2->data().property(coords_property)
                                              .asA<AtomCoords>();
                                         
        Vector center = coords.array().aaBox().center();

        NVector &distances = point_distances[molnum];
    
        BOOST_ASSERT(distances.count() == npoints);
    
        double *distances_array = distances.data();
    
        for (int i=0; i<npoints; ++i)
        {
            double dist2 = space().calcDist2( center, points_array[i] );
            
            if (dist2 != distances_array[i])
            {
                distances_array[i] = dist2;
                any_changes = true;
            }
        }
    }

    if (any_changes)
        distances_changed = true;
}

/** Internal function used to recalculate all of the distances^2 
    between the molecules in 'mol_to_molnum' and all of the identity points */
void FewPointsHelper::recalculateDistances()
{
    point_distances = QHash<MolNum,NVector>();
    
    const Molecules &molecules = this->moleculeGroup().molecules();
    
    const int nmols = mol_to_molnum.count();
    const MolNum *mol_to_molnum_array = mol_to_molnum.constData();
    
    point_distances.reserve(nmols);
    
    const int npoints = identity_points.count();
    const PointPtr *const_points_array = identity_points.constData();
    
    QVector< Vector > points(npoints);
    Vector *points_array = points.data();
    
    for (int i=0; i<npoints; ++i)
    {
        points_array[i] = const_points_array[i].read().point();
    }
    
    const PropertyName &coords_property = map["coordinates"];
    
    for (int i=0; i<nmols; ++i)
    {
        const MolNum &molnum = mol_to_molnum_array[i];
    
        Molecules::const_iterator it = molecules.constFind(molnum);
        BOOST_ASSERT( it != molecules.constEnd() );
    
        const AtomCoords &coords = it->data().property(coords_property)
                                             .asA<AtomCoords>();
                                             
        Vector center = coords.array().aaBox().center();
        
        NVector distances(npoints);
        double *distances_array = distances.data();
        
        for (int j=0; j<npoints; ++j)
        {
            double dist2 = space().calcDist2(center, points_array[j]);
            distances_array[j] = dist2;
        }
        
        point_distances.insert( molnum, distances );
    }
    
    distances_changed = true;
}

/** Internal function used to get the identities of the candidate molecules
    for the points (the first 'npoints' molecules in the molecule group, plus
    the closest 'npoints' molecules to each identity point) and the mapping
    of candidate molecule index to molecule number (with molecule number
    matching the order in the molecule group) */
void FewPointsHelper::rebuildMolToMolNum()
{
    mol_to_molnum = QVector<MolNum>();
    
    const int npoints = points_with_mols.count();
    const QVector<MolNum> &molnums = molgroup.read().molNums();
    const int nmols = molnums.count();
    
    //reserve space - worst case is we have to record all npoints molecules
    //from all npoints points, together with the npoints assigned molecules
    mol_to_molnum.reserve( qMin(nmols, npoints*(npoints+1)) );
    
    //loop through the molecules in the order they appear in the 
    //molecule group and see if they are a candidate - the first 'npoints'
    //molecules are automatically candidates (as they are assigned)
    for (int i=0; i<nmols; ++i)
    {
        const MolNum &molnum = molnums.at(i);
    
        if (i < npoints)
        {
            mol_to_molnum.append(molnum);
        }
        else
        {
            //is this molecule one of the close molecules to any 
            //of the points?
            for (int j=0; j<npoints; ++j)
            {
                if (points_with_mols.at(j).isClose(molnum))
                {
                    mol_to_molnum.append(molnum);
                    break;
                }
            }
        }
    }

    mol_to_molnum.squeeze();
}

/** Internal function used to invert a mapping of points to molecules
    into a mapping of molecules to points */
static QVector<int> invert(const QVector<int> &point_to_mol)
{
    QVector<int> mol_to_point = point_to_mol;
    
    const int *point_to_mol_array = point_to_mol.constData();
    int *mol_to_point_array = mol_to_point.data();
    
    const int npoints = point_to_mol.count();
    
    for (int i=0; i<npoints; ++i)
    {
        mol_to_point_array[ point_to_mol_array[i] ] = i;
    }
    
    return mol_to_point;
}

/** This function uses the distances between all points and molecules
    stored in 'point_distances' to work out the optimum assignment
    of molecules to points such that the total distance between
    each molecule and each point is minimised */
void FewPointsHelper::assignMoleculesToPoints()
{
    if (not distances_changed)
        return;
    
    //the order of molecules may have changed - recalculate 
    //the correct order

    //use the 'mol_to_molnum' array as this holds the numbers
    //of all of the candidate molecules in the same order as
    //they appear in the molecule group
    const int nmols = mol_to_molnum.count();
    const MolNum *mol_to_molnum_array = mol_to_molnum.constData();
    
    const int npoints = identity_points.count();
    
    //construct the matrix that contains the distances between every
    //candidate molecule and every point - one molecule per row, one point
    //per column - this has to be a square matrix, so missing rows/columns
    //are given a value of 0
    NMatrix distmatrix;
    
    if (nmols == npoints)
    {
        distmatrix = NMatrix(nmols, nmols);
        distmatrix = distmatrix.transpose(); // change to row-major memory order

        for (int i=0; i<nmols; ++i)
        {
            const MolNum &molnum = mol_to_molnum_array[i];
            
            QHash<MolNum,NVector>::const_iterator 
                                       it = point_distances.constFind(molnum);
                                                
            BOOST_ASSERT( it != point_distances.constEnd() );
            
            distmatrix.setRow(i, it.value());
        }
    }
    else if (nmols > npoints)
    {
        //there are more molecules than points - we create some extra
        //points which have zero distance to all molecules
        distmatrix = NMatrix(nmols, nmols);
        distmatrix = distmatrix.transpose(); // change to row-major memory order
        
        const int nzeroes = nmols - npoints;
        NVector new_row(npoints + nzeroes);
        
        for (int i=0; i<nmols; ++i)
        {
            const MolNum &molnum = mol_to_molnum_array[i];
            
            QHash<MolNum,NVector>::const_iterator 
                                   it = point_distances.constFind(molnum);
                                                
            BOOST_ASSERT( it != point_distances.constEnd() );
            
            const NVector &distances = it.value();
            
            memcpy( new_row.data(), distances.constData(), npoints*sizeof(double) );
            
            distmatrix.setRow( i, new_row );
        }
    }
    else
    {
        //there are more points than molecules - we create some extra
        //molecules that are all equally a very long way from all of the points
        distmatrix = NMatrix(npoints, npoints, std::numeric_limits<double>::max());
        distmatrix = distmatrix.transpose(); // change to row-major memory order

        //copy the distances to the real molecules
        for (int i=0; i<nmols; ++i)
        {
            const MolNum &molnum = mol_to_molnum_array[i];
            
            QHash<MolNum,NVector>::const_iterator 
                                        it = point_distances.constFind(molnum);
                                                
            BOOST_ASSERT( it != point_distances.constEnd() );
            
            distmatrix.setRow(i, it.value());
        }
    }

    //find the smallest and second smallest absolute difference between distances
    // - this is used to see if there are any degeneracies
    //if (false)
    {
        double delta0 = std::numeric_limits<double>::max();
        double delta1 = delta0;
    
        const int nrows = distmatrix.nRows();
        const int ncolumns = qMin( distmatrix.nColumns(), npoints );
        
        const double *distmatrix_array = distmatrix.constData();
        
        for (int i=0; i<nrows; ++i)
        {
            const double *row = &(distmatrix_array[ distmatrix.offset(i,0) ]);
            
            for (int j=0; j<ncolumns-1; ++j)
            {
                for (int k=j+1; k<ncolumns; ++k)
                {
                    const double delta = std::abs(row[k] - row[j]);
                    
                    if (delta <= delta0)
                    {
                        if (delta < delta0)
                        {
                            delta1 = delta0;
                            delta0 = delta;
                        }
                    }
                    else if (delta < delta1)
                    {
                        delta1 = delta;
                    }
                }
            }
        }
    
        if (delta0 == 0)
        {
            //there are degeneracies - some molecules are an identical
            //distance from some points
            if (delta1 == std::numeric_limits<double>::max())
                //really degenerate!
                delta1 = 0.1;
    
            //add a small penalty function that make molecule 'i' prefer
            //point 'i', or the point as close to point 'i' as possible. This
            //is used to help remove degeneracies caused by using points with
            //the same coordinates. If multiple points have the same coordinates,
            //then the molecule with index closest to the index of the point
            //will be preferred - this penalty function has to be kept to 
            //about 10% of the smallest non-zero difference between distances,
            //so that it doesn't affect the assignment of points with no degeneracies

            const double scl = 0.1 * delta1 / (ncolumns * nrows);

            for (int i=0; i<nrows; ++i)
            {
                for (int j=0; j<ncolumns; ++j)
                {
                    distmatrix(i,j) += scl * (i-j) * (i-j);
                }
            }
        }
    }
    
    //now calculate optimum assignment of molecules to points that
    //minimises the total distance between each molecule and its
    //assigned point
    QVector<int> point_to_mol = solve_linear_assignment(distmatrix, true);

    //point_to_mol maps points to molecules - we need to invert this
    //so that we map molecules to points
    mol_to_point = ::invert(point_to_mol);

    //if there are more molecules than points, then the last set of 
    //molecules are not associated with points. To ensure deterministic
    //mapping the extra molecules must be ordered so that they are
    //assigned to the null points in order
    if (nmols > npoints)
    {
        qSort( mol_to_point.data() + npoints, mol_to_point.data() + nmols );
    }

    //do we have the correct arrangement? (the nth point maps
    //to the nth molecule)
    const int *mol_to_point_array = mol_to_point.constData();
    
    bool correct_order = true;
    
    int n_to_match = qMin(npoints, nmols);
    
    for (int i=0; i<n_to_match; ++i)
    {
        if (mol_to_point_array[i] != i)
        {
            correct_order = false;
            break;
        }
    }

    if (correct_order)
        //the nth point matches up with the nth molecule
        // - we can indicate this by clearing the mol_to_point array
        mol_to_point = QVector<int>();

    distances_changed = false;
}

/** Copy the coordinates of the molecule with number 'mol_with_coords' to
    the molecule with number 'molnum' */
static Molecule swapCoordinatesTo(const Molecules &molecules,
                                  MolNum molnum,
                                  MolNum mol_with_coords,
                                  const PropertyName &coords_property)
{
    Molecule molecule = molecules[molnum].molecule();
    
    if (molnum != mol_with_coords)
    {
        return molecule.edit()
                       .setProperty(coords_property,
                                    molecules[mol_with_coords].data()
                                                              .property(coords_property)
                                   ).commit();
    }
    else
        return molecule;
}

/** Internal function that uses the molecule-point distances calculated
    and stored in 'point_distances' to work out which are the best molecules
    to maintain the constraint. This then returns which molecules must
    change to maintain the constraint */
Molecules FewPointsHelper::applyConstraint() const
{
    if (mol_to_point.isEmpty())
        //nothing needs to be changed as the constraint is satisfied
        return Molecules();

    //qDebug() << "ORDER IS INCORRECT" << mol_to_point;
        
    //the order is incorrect - we need to swap the molecules around
    const Molecules molecules = this->moleculeGroup().molecules();
    
    int n_to_match = qMin( identity_points.count(), molecules.count() );
    
    const int *mol_to_point_array = mol_to_point.constData();
    
    //the match uses mol_to_molnum, as only a subset of molecules
    //are candidate molecules
    const MolNum *mol_to_molnum_array = mol_to_molnum.constData();
    
    Molecules changed_mols;
    
    const PropertyName &coords_property = map["coordinates"];
    
    QVarLengthArray<int,4> lost_coords;
    QVarLengthArray<int,4> lost_mols;
    
    for (int i=0; i<n_to_match; ++i)
    {
        const int new_i = mol_to_point_array[i];
    
        if (new_i != i)
        {
            //we need to swap coordinates so that the ith candidate molecule is 
            //associated with the ith point
            //qDebug() << "Giving coordinates of" << new_i << "to molecule" << i;
            changed_mols.add( ::swapCoordinatesTo(molecules,
                                                  mol_to_molnum_array[i],
                                                  mol_to_molnum_array[new_i],
                                                  coords_property) );

            if (new_i >= n_to_match)
                //molecule at 'new_i' will not get new coordinates
                //unless further action is taken...
                lost_mols.append(new_i);
        
            //we now need to find the new index of the ith molecule, to
            //see if we need to swap it now (as it has moved to not be associated
            //with a point)
            int new_j = mol_to_point.indexOf(i, n_to_match);
            
            if (new_j != -1)
                //the coordinates of molecule i will be lost unless
                //further action is taken (they need to be matched up
                //with one of the molecules in 'lost_mols'
                lost_coords.append(i);
        }
    }

    BOOST_ASSERT( lost_mols.count() == lost_coords.count() );
    
    int n = lost_mols.count();
    
    for (int j=0; j<n; ++j)
    {
        const int i = lost_coords[j];
        const int new_i = lost_mols[j];
    
        //qDebug() << "...and coordinates of" << i << "to molecule" << new_i;
            
        changed_mols.add( ::swapCoordinatesTo(molecules,
                                              mol_to_molnum_array[new_i],
                                              mol_to_molnum_array[i],
                                              coords_property) );
    }
    
    //qDebug() << "IdentityConstraint" << changed_mols.molNums();
    
    return changed_mols;
}

/** Update this constraint from the passed system and return the molecules
    that need to change to maintain this constraint */
Molecules FewPointsHelper::update(const System &system, bool new_system)
{
    if (identity_points.isEmpty())
        return Molecules();

    bool new_group = this->updateGroup(system);
    bool new_points = this->updatePoints(system);

    bool new_space = false;

    if (new_system)
        new_space = this->updateSpace(system);

    if (new_group or new_points or new_space)
    {
        //update all of the points
        const int npoints = points_with_mols.count();

        bool closemols_changed = false;
    
        for (int i=0; i<npoints; ++i)
        {
            bool this_changed = points_with_mols[i].update(system);
            closemols_changed = closemols_changed or this_changed;
        }

        if (closemols_changed)
            this->rebuildMolToMolNum();

        this->recalculateDistances();
    }

    this->assignMoleculesToPoints();

    return this->applyConstraint();
}

/** Update this constraint from the passed system and return the molecules
    that need to change to maintain this constraint */
Molecules FewPointsHelper::update(const System &system, MolNum changed_mol, 
                                  bool new_system)
{
    if (identity_points.isEmpty())
        return Molecules();

    else if (point_distances.isEmpty())
        return this->update(system, new_system);

    bool new_group = this->updateGroup(system);
    bool new_points = this->updatePoints(system);

    bool new_space = false;
    
    if (new_system)
        new_space = this->updateSpace(system);

    if (new_points or new_space or new_group)
    {
        const int npoints = points_with_mols.count();
        bool closemols_changed = false;

        if (new_points or new_space)
        {
            //update all of the points
            for (int i=0; i<npoints; ++i)
            {
                bool this_changed = points_with_mols[i].update(system);
                closemols_changed = closemols_changed or this_changed;
            }
        }
        else
        {
            //guided update of all of the points
            for (int i=0; i<npoints; ++i)
            {
                bool this_changed = points_with_mols[i].update(system, changed_mol);
                closemols_changed = closemols_changed or this_changed;
            }
        }

        if (closemols_changed)
            this->rebuildMolToMolNum();

        if (closemols_changed or new_points or new_space)
            this->recalculateDistances();
        else
            this->recalculateDistances(changed_mol);
    }

    this->assignMoleculesToPoints();
    return this->applyConstraint();
}

/** Update this constraint from the passed system and return the molecules
    that need to change to maintain this constraint */
Molecules FewPointsHelper::update(const System &system, const Molecules &molecules, 
                                  bool new_system)
{
    if (molecules.isEmpty())
        return this->update(system, new_system);

    else if (identity_points.isEmpty())
        return Molecules();

    else if (point_distances.isEmpty())
        return this->update(system, new_system);

    bool new_group = this->updateGroup(system);
    bool new_points = this->updatePoints(system);

    bool new_space = false;
    
    if (new_system)
        new_space = this->updateSpace(system);

    if (new_points or new_space or new_group)
    {
        const int npoints = points_with_mols.count();
        bool closemols_changed = false;

        if (new_points or new_space)
        {
            //update all of the points
            for (int i=0; i<npoints; ++i)
            {
                bool this_changed = points_with_mols[i].update(system);
                closemols_changed = closemols_changed or this_changed;
            }
        }
        else
        {
            //guided update of all of the points
            for (int i=0; i<npoints; ++i)
            {
                bool this_changed = points_with_mols[i].update(system, molecules);
                closemols_changed = closemols_changed or this_changed;
            }
        }

        if (closemols_changed)
            this->rebuildMolToMolNum();

        if (closemols_changed or new_points or new_space)
            this->recalculateDistances();
        else
            this->recalculateDistances(molecules);
    }

    this->assignMoleculesToPoints();
    return this->applyConstraint();
}

/////////
///////// Implementation of ManyPointsHelper
/////////

/** Constructor */
ManyPointsHelper::ManyPointsHelper() : IdentityConstraintPvt(), distances_changed(false)
{}

/** Construct to constrain the identity of all of the molecules of the passed group
     - this automatically generates the points from the current configuration */
ManyPointsHelper::ManyPointsHelper(const MoleculeGroup &molgroup,
                                   const PropertyMap &map)
                 : IdentityConstraintPvt(), distances_changed(false)
{
    if (molgroup.isEmpty())
        return;

    int nmols = molgroup.nMolecules();

    QVector<PointPtr> points( nmols );
    points.squeeze();
    
    PointPtr *points_array = points.data();
    
    const PropertyName &coords_property = map["coordinates"];
    
    for (int i=0; i<nmols; ++i)
    {
        const ViewsOfMol &mol = molgroup.moleculeAt(i);
        
        const AtomCoords &coords = mol.data().property(coords_property)
                                             .asA<AtomCoords>();
                                             
        points_array[i] = new VectorPoint(coords.array().aaBox().center());
    }
    
    this->validateGroup(molgroup);
    
    points.squeeze();
    identity_points = points;
}

/** Construct to constrain the molecules in 'molgroup' to keep their
    identities using the points in 'points' */
ManyPointsHelper::ManyPointsHelper(const MoleculeGroup &molgroup,
                                   const QVector<PointPtr> &points,
                                   const PropertyMap &map)
                 : IdentityConstraintPvt(molgroup, points, map), distances_changed(false)
{}

/** Copy constructor */                          
ManyPointsHelper::ManyPointsHelper(const ManyPointsHelper &other)
                 : IdentityConstraintPvt(other),
                   point_distances(other.point_distances),
                   mol_to_point(other.mol_to_point),
                   distances_changed(other.distances_changed)
{}

/** Destructor */
ManyPointsHelper::~ManyPointsHelper()
{}

const char* ManyPointsHelper::what() const
{
    return ManyPointsHelper::typeName();
}

ManyPointsHelper* ManyPointsHelper::clone() const
{
    return new ManyPointsHelper(*this);
}

/** Recalculate the distances between the molecule with number
    'molnum' and all of the points (if the molecule is in the
    group affected by this constraint) */
void ManyPointsHelper::recalculateDistances(MolNum molnum)
{
    const Molecules &molecules = this->moleculeGroup().molecules();
    
    if (not molecules.contains(molnum))
        return;

    const PropertyName &coords_property = map["coordinates"];

    const AtomCoords &coords = molecules[molnum].data().property(coords_property)
                                                .asA<AtomCoords>();
                                             
    Vector center = coords.array().aaBox().center();
        
    NVector distances = point_distances.value(molnum);
    
    const double *const_distances_array = distances.constData();
    
    int npoints = identity_points.count();
    const PointPtr *points_array = identity_points.constData();
    
    bool these_distances_changed = false;
    
    for (int i=0; i<npoints; ++i)
    {
        double dist2 = space().calcDist2(center, points_array[i].read().point());
        
        if (dist2 != const_distances_array[i])
        {
            distances.data()[i] = dist2;
            these_distances_changed = true;
        }
    }
    
    if (these_distances_changed)
    {
        point_distances.insert(molnum, distances);
        distances_changed = true;
    }
}

/** Recalculate all of the distances of the molecules from 'molecules'
    that are in the molecule group affected by this constraint */
void ManyPointsHelper::recalculateDistances(const Molecules &new_molecules)
{
    if (new_molecules.isEmpty())
        return;
        
    else if (new_molecules.count() == 1)
    {
        this->recalculateDistances( new_molecules.constBegin()->number() );
    }

    int npoints = identity_points.count();
    const PointPtr *const_points_array = identity_points.constData();
    
    QVector< Vector > points(npoints);
    Vector *points_array = points.data();
    
    for (int i=0; i<npoints; ++i)
    {
        points_array[i] = const_points_array[i].read().point();
    }
    
    const PropertyName &coords_property = map["coordinates"];

    const Molecules &molecules = this->moleculeGroup().molecules();

    for (Molecules::const_iterator it = new_molecules.constBegin();
         it != new_molecules.constEnd();
         ++it)
    {
        MolNum molnum = it.key();

        if (not molecules.contains(molnum))
            continue;
    
        const AtomCoords &coords = molecules[molnum].data().property(coords_property)
                                                    .asA<AtomCoords>();
        
        Vector center = coords.array().aaBox().center();
            
        NVector distances = point_distances.value(molnum);
        
        const double *const_distances_array = distances.constData();
        
        bool these_distances_changed = false;
        
        for (int i=0; i<npoints; ++i)
        {
            double dist2 = space().calcDist2(center, points_array[i]);
            
            if (dist2 != const_distances_array[i])
            {
                distances.data()[i] = dist2;
                these_distances_changed = true;
            }
        }
        
        if (these_distances_changed)
        {
            point_distances.insert(molnum, distances);
            distances_changed = true;
        }
    }
}

/** Recalculate all of the distances between all of the molecules
    and all of the points */
void ManyPointsHelper::recalculateDistances()
{
    point_distances = QHash<MolNum,NVector>();
    
    const Molecules &molecules = this->moleculeGroup().molecules();
    
    point_distances.reserve(molecules.nMolecules());
    
    int npoints = identity_points.count();
    const PointPtr *const_points_array = identity_points.constData();
    
    QVector< Vector > points(npoints);
    Vector *points_array = points.data();
    
    for (int i=0; i<npoints; ++i)
    {
        points_array[i] = const_points_array[i].read().point();
    }
    
    const PropertyName &coords_property = map["coordinates"];
    
    for (Molecules::const_iterator it = molecules.constBegin();
         it != molecules.constEnd();
         ++it)
    {
        const AtomCoords &coords = it->data().property(coords_property)
                                             .asA<AtomCoords>();
                                             
        Vector center = coords.array().aaBox().center();
        
        NVector distances(npoints);
        double *distances_array = distances.data();
        
        for (int i=0; i<npoints; ++i)
        {
            double dist2 = space().calcDist2(center, points_array[i]);
            distances_array[i] = dist2;
        }
        
        point_distances.insert( it.key(), distances );
    }
    
    distances_changed = true;
}

/** This function uses the distances between all points and molecules
    stored in 'point_distances' to work out the optimum assignment
    of molecules to points such that the total distance between
    each molecule and each point is minimised */
void ManyPointsHelper::assignMoleculesToPoints()
{
    if (not distances_changed)
        return;
    
    //the order of molecules may have changed - recalculate 
    //the correct order

    //get the current order of the molecules in the group - this 
    //is used to ensure that the nth molecule in the group is
    //allocated to the nth point
    const QVector<MolNum> &molnums = this->moleculeGroup().molNums();
    const MolNum *molnums_array = molnums.constData();
    
    const int nmols = molnums.count();

    const int npoints = identity_points.count();
    
    //construct the matrix that contains the distances between every
    //molecule and every point - one molecule per row, one point
    //per column - this has to be a square matrix, so missing rows/columns
    //are given a value of 0
    NMatrix distmatrix;
    
    if (nmols == npoints)
    {
        distmatrix = NMatrix(nmols, nmols);
        distmatrix.transpose(); // change to row-major memory order
        
        for (int i=0; i<nmols; ++i)
        {
            const MolNum &molnum = molnums_array[i];
            
            QHash<MolNum,NVector>::const_iterator 
                                       it = point_distances.constFind(molnum);
                                                
            BOOST_ASSERT( it != point_distances.constEnd() );
            
            distmatrix.setRow(i, it.value());
        }
    }
    else if (nmols > npoints)
    {
        //there are more molecules than points - we create some extra
        //points which have zero distance to all molecules
        distmatrix = NMatrix(nmols, nmols);
        distmatrix.transpose(); // change to row-major memory order
        
        const int nzeroes = nmols - npoints;

        NVector new_row(npoints + nzeroes, 0.0);
        
        for (int i=0; i<nmols; ++i)
        {
            const MolNum &molnum = molnums_array[i];
            
            QHash<MolNum,NVector>::const_iterator 
                                     it = point_distances.constFind(molnum);
                                                
            BOOST_ASSERT( it != point_distances.constEnd() );
            
            const NVector &distances = it.value();
            
            BOOST_ASSERT( distances.count() == npoints );
            
            memcpy( new_row.data(), distances.constData(), npoints*sizeof(double) );

            distmatrix.setRow(i, new_row);
        }
    }
    else
    {
        //there are more points than molecules - we create some extra
        //molecules that are all equally a very long way from all of the points
        distmatrix = NMatrix(npoints, npoints, std::numeric_limits<double>::max());
        distmatrix.transpose(); // change to row-major memory order

        //copy the distances to the real molecules
        for (int i=0; i<nmols; ++i)
        {
            const MolNum &molnum = molnums_array[i];
            
            QHash<MolNum,NVector>::const_iterator 
                                        it = point_distances.constFind(molnum);
                                                
            BOOST_ASSERT( it != point_distances.constEnd() );
            
            distmatrix.setRow(i, it.value());
        }
    }

    //now calculate optimum assignment of molecules to points that
    //minimises the total distance between each molecule and its
    //assigned point
    QVector<int> point_to_mol = solve_linear_assignment(distmatrix);
    
    //point_to_mol maps points to molecules - we need to swap this
    //so that mol_to_point maps molecules to points
    mol_to_point = ::invert(point_to_mol);

    //do we have the correct arrangement? (the nth point maps
    //to the nth molecule)
    const int *mol_to_point_array = mol_to_point.constData();
    
    bool correct_order = true;
    
    int n_to_match = qMin(npoints, nmols);
    
    for (int i=0; i<n_to_match; ++i)
    {
        if (mol_to_point_array[i] != i)
        {
            correct_order = false;
            break;
        }
    }

    if (correct_order)
        //the nth point matches up with the nth molecule
        // - we can indicate this by clearing the mol_to_point array
        mol_to_point = QVector<int>();

    distances_changed = false;
}

/** Internal function that uses the molecule-point distances calculated
    and stored in 'point_distances' to work out which are the best molecules
    to maintain the constraint. This then returns which molecules must
    change to maintain the constraint */
Molecules ManyPointsHelper::applyConstraint() const
{
    if (mol_to_point.isEmpty())
        //nothing needs to be changed as the constraint is satisfied
        return Molecules();
        
    //the order is incorrect - we need to swap the molecules around
    const Molecules molecules = this->moleculeGroup().molecules();
    
    int n_to_match = qMin( identity_points.count(), molecules.count() );
    
    const int *mol_to_point_array = mol_to_point.constData();
    
    const QVector<MolNum> &molnums = this->moleculeGroup().molNums();
    const MolNum *molnums_array = molnums.constData();
    
    Molecules changed_mols;
    
    const PropertyName &coords_property = map["coordinates"];
    
    QVarLengthArray<int, 4> lost_mols;
    QVarLengthArray<int, 4> lost_coords;
    
    for (int i=0; i<n_to_match; ++i)
    {
        const int new_i = mol_to_point_array[i];
    
        if (new_i != i)
        {
            //we need to swap coordinates so that the ith molecule is 
            //associated with the ith point
            changed_mols.add( ::swapCoordinatesTo(molecules,
                                                  molnums_array[i],
                                                  molnums_array[new_i],
                                                  coords_property) );

            if (new_i >= n_to_match)
                //molecule at 'new_i' will not get new coordinates
                //unless further action is taken...
                lost_mols.append(new_i);
        
            //we now need to find the new index of the ith molecule, to
            //see if we need to swap it now (as it has moved to not be associated
            //with a point)
            int new_j = mol_to_point.indexOf(i, n_to_match);
            
            if (new_j != -1)
                //the coordinates of molecule i will be lost unless
                //further action is taken (they need to be matched up
                //with one of the molecules in 'lost_mols'
                lost_coords.append(i);
        }
    }

    BOOST_ASSERT( lost_mols.count() == lost_coords.count() );
    
    int n = lost_mols.count();
    
    for (int j=0; j<n; ++j)
    {
        const int i = lost_coords[j];
        const int new_i = lost_mols[j];
            
        changed_mols.add( ::swapCoordinatesTo(molecules,
                                              molnums_array[new_i],
                                              molnums_array[i],
                                              coords_property) );
    }
    
    return changed_mols;
}

/** Update this constraint, returning what needs to change */
Molecules ManyPointsHelper::update(const System &system, bool new_system)
{
    if (identity_points.isEmpty())
        return Molecules();

    bool new_group = this->updateGroup(system);
    bool new_points = this->updatePoints(system);

    bool new_space = false;
    
    if (new_system)
        new_space = this->updateSpace(system);

    if (new_group or new_points or new_space)
    {
        this->recalculateDistances();
    }

    this->assignMoleculesToPoints();
    return this->applyConstraint();
}

/** Update this constraint, returning what needs to change */
Molecules ManyPointsHelper::update(const System &system, MolNum changed_mol,
                                   bool new_system)
{
    if (identity_points.isEmpty())
        return Molecules();

    else if (point_distances.isEmpty())
        return this->update(system, true);

    bool new_group = this->updateGroup(system);
    bool new_points = this->updatePoints(system);

    bool new_space = false;
    
    if (new_system)
        new_space = this->updateSpace(system);

    if (new_points or new_space)
    {
        this->recalculateDistances();
    }
    else if (new_group)
    {
        this->recalculateDistances(changed_mol);
    }

    this->assignMoleculesToPoints();
    return this->applyConstraint();
}

/** Update this constraint, returning what needs to change */
Molecules ManyPointsHelper::update(const System &system, const Molecules &molecules,
                                   bool new_system)
{
    if (identity_points.isEmpty())
        return Molecules();

    else if (point_distances.isEmpty())
        return this->update(system, true);

    bool new_group = this->updateGroup(system);
    bool new_points = this->updatePoints(system);

    bool new_space = false;
    
    if (new_system)
        new_space = this->updateSpace(system);

    if (new_points or new_space)
    {
        this->recalculateDistances();
    }
    else if (new_group)
    {
        this->recalculateDistances(molecules);
    }

    this->assignMoleculesToPoints();
    return this->applyConstraint();
}

/////////
///////// Implementation of IdentityConstraint
/////////

static const RegisterMetaType<IdentityConstraint> r_identityconstraint;

/** Serialise to a binary datastream */
QDataStream SIRESYSTEM_EXPORT &operator<<(QDataStream &ds,
                                          const IdentityConstraint &identityconstraint)
{
    writeHeader(ds, r_identityconstraint, 1);
    
    SharedDataStream sds(ds);
    
    sds << identityconstraint.d->molGroupPtr()
        << identityconstraint.d->points()
        << identityconstraint.d->propertyMap()
        << static_cast<const MoleculeConstraint&>(identityconstraint);
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIRESYSTEM_EXPORT &operator>>(QDataStream &ds,
                                          IdentityConstraint &identityconstraint)
{
    VersionID v = readHeader(ds, r_identityconstraint);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        MolGroupPtr molgroup;
        QVector<PointPtr> points;
        PropertyMap map;
        
        sds >> molgroup >> points >> map;

        IdentityConstraint new_constraint(points, molgroup, map);
        
        sds >> static_cast<MoleculeConstraint&>(new_constraint);
        
        identityconstraint = new_constraint;
    }
    else
        throw version_error( v, "1", r_identityconstraint, CODELOC );
        
    return ds;
}

/** Constructor */
IdentityConstraint::IdentityConstraint()
                   : ConcreteProperty<IdentityConstraint,MoleculeConstraint>(),
                     d( IdentityConstraintPvt::null() ),
                     space_property("space")
{}

/** Construct the constraint that constrains the identities of all
    of the molecules in the passed molecule group. This uses the current
    locations of the molecules to apply the constraint. The (optionally
    supplied) property map is used to find the properties required
    of this constraint */
IdentityConstraint::IdentityConstraint(const MoleculeGroup &molgroup,
                                       const PropertyMap &map)
                   : ConcreteProperty<IdentityConstraint,MoleculeConstraint>(),
                     d( static_cast<IdentityConstraintPvt*>(
                                new ManyPointsHelper(molgroup, map) ) )
{
    space_property = map["space"];
}

/** Construct the constraint that constrains the identity of a single
    molecule in the passed molecule group - this sets the identity 
    of the first molecule to be that of the one closest to the
    passed point. The (optionally supplied) property map is used to 
    find the properties required of this constraint */
IdentityConstraint::IdentityConstraint(const PointRef &point,
                                       const MoleculeGroup &molgroup,
                                       const PropertyMap &map)
                   : ConcreteProperty<IdentityConstraint,MoleculeConstraint>(),
                     d( static_cast<IdentityConstraintPvt*>(
                                new SinglePointHelper(molgroup, point, map) ) )
{
    space_property = map["space"];
}

/** Construct the constraint that constrains the identities of the 
    points.count() molecules from the passed molecule group so that
    the first molecule is identified by the first point, the second
    molecule is identified by the second point, and the nth molecule
    is identified by the nth point. The (optionally supplied) property 
    map is used to find the properties required of this constraint */
IdentityConstraint::IdentityConstraint(const QVector<PointPtr> &points,
                                       const MoleculeGroup &molgroup,
                                       const PropertyMap &map)
                   : ConcreteProperty<IdentityConstraint,MoleculeConstraint>()
{
    if (points.isEmpty())
        d = IdentityConstraintPvt::null();
    
    else if (points.count() == 1)
        d = static_cast<IdentityConstraintPvt*>(
                    new SinglePointHelper(molgroup, points.first(), map) );
    
    else if (points.count() > 0.5 * molgroup.nMolecules())
        d = static_cast<IdentityConstraintPvt*>(
                    new ManyPointsHelper(molgroup, points, map) );
        
    else
        d = static_cast<IdentityConstraintPvt*>(
                    new FewPointsHelper(molgroup, points, map) );

    space_property = map["space"];
}                  

/** Function used for debugging that switches this object over
    to using the many points algorithm to apply the constraint */
void IdentityConstraint::useManyPointsAlgorithm()
{
    d = static_cast<IdentityConstraintPvt*>(
                new ManyPointsHelper( this->moleculeGroup(),
                                      this->points(),
                                      this->propertyMap() ) );

    this->setSystem( System() );
}

/** Function used for debugging that switches this object over
    to using the few points algorithm to apply the constraint */
void IdentityConstraint::useFewPointsAlgorithm()
{
    d = static_cast<IdentityConstraintPvt*>(
                new FewPointsHelper( this->moleculeGroup(),
                                     this->points(),
                                     this->propertyMap() ) );

    this->setSystem( System() );
}

/** Function used for debugging that switches this object over
    to using the single point algorithm to apply the constraint 
    
    \throw SireError::invalid_state
*/
void IdentityConstraint::useSinglePointAlgorithm()
{
    if (this->points().count() != 1)
        throw SireError::invalid_state( QObject::tr(    
                "The single point algorithm can only be used when there "
                "is just a single point - not when there are %1 points.")
                    .arg(this->points().count()), CODELOC );
                    
    d = static_cast<IdentityConstraintPvt*>(
                        new SinglePointHelper( this->moleculeGroup(),
                                               this->points().first(),
                                               this->propertyMap() ) );

    this->setSystem( System() );
}
  
/** Construct the constraint that constrains the identities of the 
    points.count() molecules from the passed molecule group so that
    the first molecule is identified by the first point, the second
    molecule is identified by the second point, and the nth molecule
    is identified by the nth point. The (optionally supplied) property 
    map is used to find the properties required of this constraint */
IdentityConstraint::IdentityConstraint(const QList<PointPtr> &points,
                                       const MoleculeGroup &molgroup,
                                       const PropertyMap &map)
                   : ConcreteProperty<IdentityConstraint,MoleculeConstraint>()
{
    if (points.isEmpty())
        d = IdentityConstraintPvt::null();
    
    else if (points.count() == 1)
        d = static_cast<IdentityConstraintPvt*>(
                        new SinglePointHelper(molgroup, points.first(), map) );
    
    else
    {
        IdentityConstraint::operator=( IdentityConstraint(points.toVector(),
                                                          molgroup, map) );
    }
    
    space_property = map["space"];
}

/** Copy constructor */
IdentityConstraint::IdentityConstraint(const IdentityConstraint &other)
                   : ConcreteProperty<IdentityConstraint,MoleculeConstraint>(other),
                     d(other.d), 
                     space_property(other.space_property),
                     changed_mols(other.changed_mols)
{}

/** Destructor */
IdentityConstraint::~IdentityConstraint()
{}

/** Copy assignment operator */
IdentityConstraint& IdentityConstraint::operator=(const IdentityConstraint &other)
{
    if (this != &other)
    {
        MoleculeConstraint::operator=(other);
        d = other.d;
        space_property = other.space_property;
        changed_mols = other.changed_mols;
    }
    
    return *this;
}

/** Comparison operator */
bool IdentityConstraint::operator==(const IdentityConstraint &other) const
{
    return this == &other or
           (MoleculeConstraint::operator==(other) and
            QLatin1String(d->what()) == QLatin1String(other.d->what()) and
            this->moleculeGroup() == other.moleculeGroup() and
            this->points() == other.points() and
            this->propertyMap() == other.propertyMap());
}

/** Comparison operator */
bool IdentityConstraint::operator!=(const IdentityConstraint &other) const
{
    return not IdentityConstraint::operator==(other);
}

const char* IdentityConstraint::typeName()
{
    return QMetaType::typeName( qMetaTypeId<IdentityConstraint>() );
}

/** Return a string representation of this constraint */
QString IdentityConstraint::toString() const
{
    return QString("IdentityConstraint( nPoints() == %1, nMolecules() == %2 )")
                .arg(this->points().count())
                .arg(this->moleculeGroup().nMolecules());
}

/** Return the molecule group acted on by this constraint */
const MoleculeGroup& IdentityConstraint::moleculeGroup() const
{
    return d->moleculeGroup();
}

/** Return the points used to identify the molecules */
QVector<SireFF::PointPtr> IdentityConstraint::points() const
{
    return d->points();
}

/** Return the property map used to find the properties used
    by this constraint */
const PropertyMap& IdentityConstraint::propertyMap() const
{
    return d->propertyMap();
}

/** Update this constraint so that it is applied to the system 'system' */
void IdentityConstraint::setSystem(const System &system)
{
    if ( Constraint::wasLastSystem(system) and Constraint::wasLastSubVersion(system) )
        return;

    else if (d.constData() == 0)
    {
        Constraint::setSatisfied(system, true);
        return;
    }

    changed_mols = d->update(system, true);
    
    Constraint::setSatisfied(system, changed_mols.isEmpty());
}

static bool pointsChanged(const QVector<PointPtr> &points, 
                          const Delta &delta, quint32 last_subversion)
{
    for (QVector<PointPtr>::const_iterator it = points.constBegin();
         it != points.constEnd();
         ++it)
    {
        if (delta.sinceChanged(it->read(), last_subversion))
            return true;
    }

    return false;
}

bool IdentityConstraint::mayChange(const Delta &delta, quint32 last_subversion) const
{
    if ( d.constData() == 0 or d->isNull() )
        return false;

    else if (not changed_mols.isEmpty())
        return true;
        
    else
        return delta.sinceChanged(space_property, last_subversion) or
               delta.sinceChanged(moleculeGroup().molecules(), last_subversion) or
               ::pointsChanged(d->points(), delta, last_subversion);
}

bool IdentityConstraint::fullApply(Delta &delta)
{
    if (d.constData() == 0)
        return false;

    this->setSystem(delta.deltaSystem());

    bool changed = false;
    
    int i = 0;
    while (not changed_mols.isEmpty())
    {
        bool this_changed = delta.update(changed_mols);
        changed = changed or this_changed;
        changed_mols = d->update(delta.deltaSystem(), changed_mols, false);

        ++i;
        
        if (i > 10)
            throw SireSystem::constraint_error( QObject::tr(
                    "The identity constraint could not be solved self-consistently!"),
                        CODELOC );
    }

    return changed;
}

bool IdentityConstraint::deltaApply(Delta &delta, quint32 last_subversion)
{
    if (d.constData() == 0)
        return false;

    else if ( (not changed_mols.isEmpty() ) or
               delta.sinceChanged(space_property, last_subversion) or
               ::pointsChanged(d->points(), delta, last_subversion) )
    {
        return this->fullApply(delta);
    }
    else
    {
        QTime t;
        t.start();

        QList<MolNum> changed_molnums = delta.changedMoleculesSince(
                                                moleculeGroup().molecules(),
                                                last_subversion);

        if (changed_molnums.isEmpty())
            return false;
    
        else if (changed_molnums.count() == 1)
        {
            changed_mols = d->update(delta.deltaSystem(), changed_molnums.at(0),
                                     false);
        }
        else
        {
            changed_mols = d->update(delta.deltaSystem(), changed_molnums, false);
        }
        
        bool changed = false;
        int i = 0;
        
        while (not changed_mols.isEmpty())
        {
            bool this_changed = delta.update(changed_mols);
            changed = changed or this_changed;
            changed_mols = d->update(delta.deltaSystem(), changed_mols, false);
        
            ++i;
        
            if (i > 10)
                throw SireSystem::constraint_error( QObject::tr(
                        "The identity constraint could not be solved self-consistently!"),
                            CODELOC );
        }
        
        return changed;
    }
}

static MolGroupPtr constrain(const MoleculeGroup &molgroup, 
                             IdentityConstraint &constraint,
                             const PropertyMap &map)
{
    System tmp_system;
    tmp_system.add(molgroup);
    
    if (not map["space"].hasValue())
        tmp_system.setProperty("space", Cartesian());
    
    tmp_system.add(constraint);
    tmp_system.applyConstraints();
    
    if (not tmp_system.constraintsSatisfied())
        qDebug() << "WARNING - constraints not satisfied!";
    
    MolGroupPtr new_molgroup( molgroup );
    new_molgroup.edit().update(tmp_system.molecules());
    
    return new_molgroup;
}

/** Static function used to constrain the identities of the molecules
    in 'molgroup' against the point 'point'. This makes the first molecule
    in the group have the identity that matches this point */
MolGroupPtr IdentityConstraint::constrain(const MoleculeGroup &molgroup,
                                          const PointRef &point,
                                          const PropertyMap &map)
{
    IdentityConstraint constraint(point, molgroup, map);
    return ::constrain(molgroup, constraint, map);
}

/** Static function used to constrain the identities of the molecules
    in 'molgroup' against the identity points in 'points' - the
    first npoints molecules in the group are constrained in order
    against the points */
MolGroupPtr IdentityConstraint::constrain(const MoleculeGroup &molgroup,
                                          const QVector<PointPtr> &points,
                                          const PropertyMap &map)
{
    IdentityConstraint constraint(points, molgroup, map);
    return ::constrain(molgroup, constraint, map);
}

/** Static function used to constrain the identities of the molecules
    in 'molgroup' against the identity points in 'points' - the
    first npoints molecules in the group are constrained in order
    against the points */
MolGroupPtr IdentityConstraint::constrain(const MoleculeGroup &molgroup,
                                          const QList<PointPtr> &points,
                                          const PropertyMap &map)
{
    IdentityConstraint constraint(points, molgroup, map);
    return ::constrain(molgroup, constraint, map);
}
