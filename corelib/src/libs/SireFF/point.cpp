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

#include "point.h"

#include "forcetable.h"

#include "SireMol/evaluator.h"
#include "SireMol/molecule.h"
#include "SireMol/moleculegroup.h"
#include "SireMol/moleculegroups.h"
#include "SireMol/mgidx.h"

#include "SireVol/aabox.h"

#include <boost/tuple/tuple.hpp>

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireVol/errors.h"

using namespace SireFF;
using namespace SireMol;
using namespace SireMaths;
using namespace SireVol;
using namespace SireBase;
using namespace SireStream;

using boost::tuples::tuple;

//////////////
////////////// Implementation of Point
//////////////

static const RegisterMetaType<Point> r_point( MAGIC_ONLY, Point::typeName() );

/** Serialise to a binary datastream */
QDataStream SIREFF_EXPORT &operator<<(QDataStream &ds, const Point &point)
{
    writeHeader(ds, r_point, 1);
    
    SharedDataStream sds(ds);
    
    sds << point.p << point.spce
        << static_cast<const Property&>(point);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREFF_EXPORT &operator>>(QDataStream &ds, Point &point)
{
    VersionID v = readHeader(ds, r_point);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
    
        sds >> point.p >> point.spce
            >> static_cast<Property&>(point);
    }
    else
        throw version_error(v, "1", r_point, CODELOC);
        
    return ds;
}

/** Constructor */
Point::Point() : Property()
{}

/** Construct with the passed initial point */
Point::Point(const Vector &point) : Property(), p(point)
{}

/** Copy constructor */
Point::Point(const Point &other) 
      : Property(other), p(other.p), spce(other.spce)
{}

/** Destructor */
Point::~Point()
{}

/** Copy assignment operator */
Point& Point::operator=(const Point &other)
{
    p = other.p;
    spce = other.spce;
    Property::operator=(other);
    
    return *this;
}

/** Comparison operator */
bool Point::operator==(const Point &other) const
{
    return p == other.p and spce.read().equals(other.spce) and 
           Property::operator==(other);
}

/** Comparison operator */
bool Point::operator!=(const Point &other) const
{
    return not Point::operator==(other);
}

/** Return the point in 3D space */
const Vector& Point::operator()() const
{
    return p;
}

/** Update the point to 'point', returning whether or not this
    changes the point */
bool Point::updatePoint(const Vector &point)
{
    if (p != point)
    {
        p = point;
        return true;
    }
    else
        return false;
}

/** Return the point in 3D space */
const Vector& Point::point() const
{
    return Point::operator()();
}

/** Return the 3D space in which this point is calculated
    (although note that this 3D point, like the molecules,
    exists in the infinite cartesian space) */
const Space& Point::space() const
{
    return spce.read();
}

/** Set the 3D space in which this point is calculated 
    (although note that this 3D point, like the molecules,
    exists in the infinite cartesian space) */
void Point::setSpace(const Space &space)
{
    spce = space;
}

/** Return whether or not the points 'point0' and 'point1' are
    both within the same molecule (so together are intra-molecule points) */
bool Point::areIntraMoleculePoints(const Point &point0, const Point &point1)
{
    if (point0.isIntraMoleculePoint() and point1.isIntraMoleculePoint())
    {
        return point0.molecules().constBegin()->number() ==
               point1.molecules().constBegin()->number();
    }
    else
        return false;
}

/** Return whether or not the points 'point0', 'point1' and 'point2' are
    within the same molecule (so together are intra-molecule points) */
bool Point::areIntraMoleculePoints(const Point &point0, const Point &point1,
                                   const Point &point2)
{
    if (point0.isIntraMoleculePoint() and point1.isIntraMoleculePoint() and
        point2.isIntraMoleculePoint())
    {
        MolNum mol0 = point0.molecules().constBegin()->number();
    
        return point1.molecules().constBegin()->number() == mol0 and
               point2.molecules().constBegin()->number() == mol0;
    }
    else
        return false;
}

/** Return whether or not the points 'point0', 'point1', 'point2' and 'point3' are
    within the same molecule (so together are intra-molecule points) */
bool Point::areIntraMoleculePoints(const Point &point0, const Point &point1,
                                   const Point &point2, const Point &point3)
{
    if (point0.isIntraMoleculePoint() and point1.isIntraMoleculePoint() and
        point2.isIntraMoleculePoint() and point3.isIntraMoleculePoint())
    {
        MolNum mol0 = point0.molecules().constBegin()->number();
    
        return point1.molecules().constBegin()->number() == mol0 and
               point2.molecules().constBegin()->number() == mol0 and
               point3.molecules().constBegin()->number() == mol0;
    }
    else
        return false;
}

/** Return whether or not the points in 'points' are all
    within the same molecule (so together are intra-molecule points) */
bool Point::areIntraMoleculePoints(const QVector<PointPtr> &points)
{
    if (points.isEmpty())
        return false;

    QVector<PointPtr>::const_iterator it = points.constBegin();
    
    MolNum mol0;
    
    for ( ; it != points.constEnd(); ++it )
    {
        if (not it->isNull())
        {
            if (it->read().isIntraMoleculePoint())
                mol0 = it->read().molecules().constBegin()->number();
            else
                return false;
                
            break;
        }
    }

    for ( ; it != points.constEnd(); ++it )
    {
        if (not it->isNull())
        {
            if (it->read().isIntraMoleculePoint())
            {
                if (it->read().molecules().constBegin()->number() != mol0)
                    return false;
            }
            else
                return false;
        }
    }
    
    return true;
}

Q_GLOBAL_STATIC( VectorPoint, vectorPoint )

const VectorPoint& Point::null()
{
    return *(vectorPoint());
}

//////////////
////////////// Implementation of AtomPoint
//////////////

static const RegisterMetaType<AtomPoint> r_atompoint;

/** Serialise to a binary datastream */
QDataStream SIREFF_EXPORT &operator<<(QDataStream &ds, const AtomPoint &atompoint)
{
    writeHeader(ds, r_atompoint, 1);
    
    SharedDataStream sds(ds);
    
    sds << atompoint.atm << atompoint.coords_property
        << static_cast<const Point&>(atompoint);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREFF_EXPORT &operator>>(QDataStream &ds, AtomPoint &atompoint)
{
    VersionID v = readHeader(ds, r_atompoint);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> atompoint.atm >> atompoint.coords_property
            >> static_cast<Point&>(atompoint);
    }
    else
        throw version_error(v, "1", r_atompoint, CODELOC);
        
    return ds;
}

/** Constructor */
AtomPoint::AtomPoint() : ConcreteProperty<AtomPoint,Point>()
{}

/** Construct for the passed atom */
AtomPoint::AtomPoint(const Atom &atom, const PropertyMap &map)
          : ConcreteProperty<AtomPoint,Point>(), atm(atom)
{
    coords_property = map["coordinates"];
    AtomPoint::update(atom.data());
}

/** Copy constructor */
AtomPoint::AtomPoint(const AtomPoint &other)
          : ConcreteProperty<AtomPoint,Point>(other), 
            atm(other.atm), coords_property(other.coords_property)
{}

/** Destructor */
AtomPoint::~AtomPoint()
{}

/** Copy assignment operator */
AtomPoint& AtomPoint::operator=(const AtomPoint &other)
{
    Point::operator=(other);
    atm = other.atm;
    coords_property = other.coords_property;
    
    return *this;
}

/** Comparison operator */
bool AtomPoint::operator==(const AtomPoint &other) const
{
    return atm == other.atm and coords_property == other.coords_property and
           Point::operator==(other);
}

/** Comparison operator */
bool AtomPoint::operator!=(const AtomPoint &other) const
{
    return atm != other.atm or coords_property != other.coords_property or
           Point::operator!=(other);
}

const char* AtomPoint::typeName()
{
    return QMetaType::typeName( qMetaTypeId<AtomPoint>() );
}

/** Return a string representation */
QString AtomPoint::toString() const
{
    return QString("AtomPoint{ %1 : %2 }").arg(atm.toString(),
                                               this->point().toString());
}

/** Update this point, returning whether or not this changes
    the location of this point 

    \throw SireBase::missing_property
    \throw SireError::incompatible_error
    \throw SireError::invalid_cast
*/
bool AtomPoint::update(const MoleculeData &moldata)
{
    if (atm.data().number() == moldata.number())
    {
        atm.update(moldata);
        return this->updatePoint(atm.property<Vector>(coords_property));
    }
    else
        return false;
}
                    
/** Update this point, returning whether or not this changes
    the location of this point 

    \throw SireBase::missing_property
    \throw SireError::incompatible_error
    \throw SireError::invalid_cast
*/
bool AtomPoint::update(const Molecules &molecules)
{
    if (molecules.contains(atm.data().number()))
    {
        return this->update( molecules[atm.data().number()].data() );
    }
    else
        return false;
}
                    
/** Update this point, returning whether or not this changes
    the location of this point 

    \throw SireBase::missing_property
    \throw SireError::incompatible_error
    \throw SireError::invalid_cast
*/
bool AtomPoint::update(const MoleculeGroup &molgroup)
{
    if (molgroup.contains(atm.data().number()))
    {
        return this->update( molgroup[atm.data().number()].data() );
    }
    else
        return false;
}
                    
/** Update this point, returning whether or not this changes
    the location of this point 

    \throw SireBase::missing_property
    \throw SireError::incompatible_error
    \throw SireError::invalid_cast
*/
bool AtomPoint::update(const MolGroupsBase &molgroups)
{
    if (molgroups.contains(atm.data().number()))
    {
        const QList<MGNum> &mgnums = molgroups.groupsContaining(atm.data().number());
        
        if (mgnums.isEmpty())
            return false;
        
        return this->update( molgroups[mgnums.first()][atm.data().number()].data() );
    }
    else
        return false;
}

/** Return whether or not this changes
    the location of this point 

    \throw SireBase::missing_property
    \throw SireError::incompatible_error
    \throw SireError::invalid_cast
*/
bool AtomPoint::wouldUpdate(const MoleculeData &moldata) const
{
    if (atm.data().number() == moldata.number())
    {
        return moldata.property(coords_property)
                      .asA<AtomCoords>()[atm.cgAtomIdx()] != this->point();
    }
    else
        return false;
}
                    
/** Return whether or not this changes
    the location of this point 

    \throw SireBase::missing_property
    \throw SireError::incompatible_error
    \throw SireError::invalid_cast
*/
bool AtomPoint::wouldUpdate(const Molecules &molecules) const
{
    if (molecules.contains(atm.data().number()))
    {
        return this->wouldUpdate( molecules[atm.data().number()].data() );
    }
    else
        return false;
}
                    
/** Return whether or not this changes
    the location of this point 

    \throw SireBase::missing_property
    \throw SireError::incompatible_error
    \throw SireError::invalid_cast
*/
bool AtomPoint::wouldUpdate(const MoleculeGroup &molgroup) const
{
    if (molgroup.contains(atm.data().number()))
    {
        return this->wouldUpdate( molgroup[atm.data().number()].data() );
    }
    else
        return false;
}
                    
/** Return whether or not this changes
    the location of this point 

    \throw SireBase::missing_property
    \throw SireError::incompatible_error
    \throw SireError::invalid_cast
*/
bool AtomPoint::wouldUpdate(const MolGroupsBase &molgroups) const
{
    if (molgroups.contains(atm.data().number()))
    {
        const QList<MGNum> &mgnums = molgroups.groupsContaining(atm.data().number());
        
        if (mgnums.isEmpty())
            return false;
        
        return this->wouldUpdate( molgroups[mgnums.first()][atm.data().number()].data() );
    }
    else
        return false;
}

/** Return the molecules needed to get this point */
Molecules AtomPoint::molecules() const
{
    Molecules mols;

    if (not atm.isEmpty())
        mols += atm;
    
    return mols;
}

/** Return the number of molecules needed to get this point */
int AtomPoint::nMolecules() const
{
    if (not atm.isEmpty())
        return 1;
    else
        return 0;
}

/** Does this point require information from the molecule with number 'molnum' */
bool AtomPoint::contains(MolNum molnum) const
{
    return atm.data().number() == molnum;
}

/** Does this point require information from the molecule with ID 'molid' */
bool AtomPoint::contains(const MolID &molid) const
{
    try
    {
        return not molid.map( this->molecules() ).isEmpty();
    }
    catch(...)
    {
        return false;
    }
}

/** Return whether or not this point uses data from any of the 
    molecules in the passed forcetable */
bool AtomPoint::usesMoleculesIn(const ForceTable &forcetable) const
{
    if (not atm.isEmpty())
        return forcetable.containsTable(atm.data().number());
    else
        return false;
}

/** Return whether or not this point uses data from any of the 
    molecules in 'molecules' */
bool AtomPoint::usesMoleculesIn(const Molecules &molecules) const
{
    if (not atm.isEmpty())
        return molecules.contains(atm.data().number());
    else
        return false;
}

/** Return whether or not this point uses data from any of the 
    molecules in the group 'molgroup' */
bool AtomPoint::usesMoleculesIn(const MoleculeGroup &molgroup) const
{
    if (not atm.isEmpty())
        return molgroup.contains(atm.data().number());
    else
        return false;
}

/** Return whether or not this point uses data from any of the 
    molecules in the groups in 'molgroups' */
bool AtomPoint::usesMoleculesIn(const MolGroupsBase &molgroups) const
{
    if (not atm.isEmpty())
        return molgroups.contains(atm.data().number());
    else
        return false;
}

/** Return the actual atom */
const Atom& AtomPoint::atom() const
{
    return atm;
}

/** Add the force acting on this atom to the passed table (if it is 
    the table for the molecule containing the atom */
bool AtomPoint::addForce(MolForceTable &molforces, const Vector &force) const
{
    if (molforces.molNum() == atm.data().number())
    {
        return molforces.add( atm.cgAtomIdx(), force );
    }
    else
        return false;
} 

/** Add the force acting on this atom to the passed table (if it contains 
    the table for the molecule containing the atom */
bool AtomPoint::addForce(ForceTable &forces, const Vector &force) const
{
    if (forces.containsTable(atm.data().number()))
    {
        return forces.getTable(atm.data().number()).add( atm.cgAtomIdx(), force );
    }
    else
        return false;
}

/** Return whether this is an intramolecular point (it depends on coordinates
    of atoms in just one molecule) */
bool AtomPoint::isIntraMoleculePoint() const
{
    return not atm.isEmpty();
}

/** Return whether or not this is an intermolecular point (it depends on
    coordinates of atoms from than one molecule) */
bool AtomPoint::isInterMoleculePoint() const
{
    return false;
}

/** Return whether or not this is an extramolecular point (it is independent
    of the coordinates of atoms in any molecule, i.e. it is just a point in space) */
bool AtomPoint::isExtraMoleculePoint() const
{
    return atm.isEmpty();
}

//////////////
////////////// Implementation of VectorPoint
//////////////

static const RegisterMetaType<VectorPoint> r_vectorpoint;

/** Serialise to a binary datastream */
QDataStream SIREFF_EXPORT &operator<<(QDataStream &ds, const VectorPoint &vectorpoint)
{
    writeHeader(ds, r_vectorpoint, 1);
    
    ds << static_cast<const Point&>(vectorpoint);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREFF_EXPORT &operator>>(QDataStream &ds, VectorPoint &vectorpoint)
{
    VersionID v = readHeader(ds, r_vectorpoint);
    
    if (v == 1)
    {
        ds >> static_cast<Point&>(vectorpoint);
    }
    else
        throw version_error( v, "1", r_vectorpoint, CODELOC );
        
    return ds;
}

/** Constructor */
VectorPoint::VectorPoint() : ConcreteProperty<VectorPoint,Point>()
{}

/** Constructor for the specified point */
VectorPoint::VectorPoint(const Vector &point)
            : ConcreteProperty<VectorPoint,Point>(point)
{}

/** Copy constructor */
VectorPoint::VectorPoint(const VectorPoint &other)
            : ConcreteProperty<VectorPoint,Point>(other)
{}

/** Destructor */
VectorPoint::~VectorPoint()
{}

/** Copy assignment operator */
VectorPoint& VectorPoint::operator=(const VectorPoint &other)
{
    Point::operator=(other);
    return *this;
}

/** Comparison operator */
bool VectorPoint::operator==(const VectorPoint &other) const
{
    return Point::operator==(other);
}

/** Comparison operator */
bool VectorPoint::operator!=(const VectorPoint &other) const
{
    return Point::operator!=(other);
}

const char* VectorPoint::typeName()
{
    return QMetaType::typeName( qMetaTypeId<VectorPoint>() );
}

/** Return a string representation */
QString VectorPoint::toString() const
{
    return QString("VectorPoint{ %1 }").arg(this->point().toString());
}

/** A VectorPoint is not updatable */
bool VectorPoint::update(const MoleculeData&)
{
    return false;
}
                    
/** A VectorPoint is not updatable */
bool VectorPoint::update(const Molecules&)
{
    return false;
}

/** A VectorPoint is not updatable */
bool VectorPoint::update(const MoleculeGroup&)
{
    return false;
}
                    
/** A VectorPoint is not updatable */
bool VectorPoint::update(const MolGroupsBase&)
{
    return false;
}

/** A VectorPoint is not updatable */
bool VectorPoint::wouldUpdate(const MoleculeData&) const
{
    return false;
}
                    
/** A VectorPoint is not updatable */
bool VectorPoint::wouldUpdate(const Molecules&) const
{
    return false;
}

/** A VectorPoint is not updatable */
bool VectorPoint::wouldUpdate(const MoleculeGroup&) const
{
    return false;
}
                    
/** A VectorPoint is not updatable */
bool VectorPoint::wouldUpdate(const MolGroupsBase&) const
{
    return false;
}

/** No molecules are needed to create this point */
Molecules VectorPoint::molecules() const
{
    return Molecules();
}

/** No molecules are needed to create this point */
int VectorPoint::nMolecules() const
{
    return 0;
}

/** No molecules are needed to create this point */
bool VectorPoint::contains(MolNum) const
{
    return false;
}

/** No molecules are needed to create this point */
bool VectorPoint::contains(const MolID&) const
{
    return false;
}

/** No molecules are needed to create this point */
bool VectorPoint::usesMoleculesIn(const ForceTable&) const
{
    return false;
}

/** No molecules are needed to create this point */
bool VectorPoint::usesMoleculesIn(const Molecules&) const
{
    return false;
}

/** No molecules are needed to create this point */
bool VectorPoint::usesMoleculesIn(const MoleculeGroup&) const
{
    return false;
}

/** No molecules are needed to create this point */
bool VectorPoint::usesMoleculesIn(const MolGroupsBase&) const
{
    return false;
}

/** No forces on a point */
bool VectorPoint::addForce(MolForceTable&, const Vector&) const
{
    return false;
} 

/** No forces on a point */
bool VectorPoint::addForce(ForceTable&, const Vector&) const
{
    return false;
}

/** Return whether this is an intramolecular point (it depends on coordinates
    of atoms in just one molecule) */
bool VectorPoint::isIntraMoleculePoint() const
{
    return false;
}

/** Return whether or not this is an intermolecular point (it depends on
    coordinates of atoms from than one molecule) */
bool VectorPoint::isInterMoleculePoint() const
{
    return false;
}

/** Return whether or not this is an extramolecular point (it is independent
    of the coordinates of atoms in any molecule, i.e. it is just a point in space) */
bool VectorPoint::isExtraMoleculePoint() const
{
    return true;
}

//////////////
////////////// Implementation of Center
//////////////

static const RegisterMetaType<Center> r_center;

/** Serialise to a binary datastream */
QDataStream SIREFF_EXPORT &operator<<(QDataStream &ds, const Center &center)
{
    writeHeader(ds, r_center, 1);
    
    SharedDataStream sds(ds);
    
    sds << center.mols << center.property_map
        << static_cast<const Point&>(center);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREFF_EXPORT &operator>>(QDataStream &ds, Center &center)
{
    VersionID v = readHeader(ds, r_center);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> center.mols >> center.property_map
            >> static_cast<Point&>(center);
    }
    else
        throw version_error( v, "1", r_center, CODELOC );
        
    return ds;
}

/** Constructor */
Center::Center() : ConcreteProperty<Center,Point>()
{}

/** Construct to get the center of the molecule view 'molview' using the
    passed property map to find the required properties 
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
Center::Center(const MoleculeView &molview, const PropertyMap &map)
       : ConcreteProperty<Center,Point>(), mols(molview), property_map(map)
{
    this->recalculatePoint();
}

/** Construct to get the center of the molecules in 'molecules', using the
    passed property map to find the required properties 
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
Center::Center(const Molecules &molecules, const PropertyMap &map)
       : ConcreteProperty<Center,Point>(), mols(molecules), property_map(map)
{
    this->recalculatePoint();
}

/** Copy constructor */
Center::Center(const Center &other)
       : ConcreteProperty<Center,Point>(other), mols(other.mols),
         property_map(other.property_map)
{}

/** Destructor */
Center::~Center()
{}

/** Copy assignment operator */
Center& Center::operator=(const Center &other)
{
    mols = other.mols;
    property_map = other.property_map;
    Point::operator=(other);
    
    return *this;
}

/** Comparison operator */
bool Center::operator==(const Center &other) const
{
    return mols == other.mols and property_map == other.property_map and
           Point::operator==(other);
}

/** Comparison operator */
bool Center::operator!=(const Center &other) const
{
    return mols != other.mols or property_map != other.property_map or
           Point::operator!=(other);
}

const char* Center::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Center>() );
}

/** Return a string representation */
QString Center::toString() const
{
    return QString("Center{ %1 : %2 }").arg(mols.toString())
                                       .arg(this->point().toString());
}

/** Set the space - if there is more than one molecule, then this
    point can only be used with a cartesian, non-periodic space */
void Center::setSpace(const Space &space)
{
    if (mols.nMolecules() > 1)
    {
        if (space.isPeriodic())
            throw SireVol::incompatible_space( QObject::tr(
                "A Center point with more than one molecule cannot be used "
                "with the periodic space \"%1\". Either use a Center point "
                "built from a single molecule, or switch to a Cartesian space.")
                    .arg(space.what()), CODELOC );
    
        if (not space.isCartesian())
            throw SireVol::incompatible_space( QObject::tr(
                "A Center point with more than one molecule cannot be used "
                "with the non-cartesian space \"%1\". Either use a Center point "
                "built from a single molecule, or switch to a Cartesian space.")
                    .arg(space.what()), CODELOC );
    }
    
    Point::setSpace(space);
}

/** Recalculate the point using the current space */
bool Center::recalculatePoint()
{
    if (mols.isEmpty())
        return false;

    Molecules::const_iterator it = mols.constBegin();
    
    AABox aabox = it->evaluate().center(property_map);
    
    for (++it; it != mols.constEnd(); ++it)
    {
        aabox += it->evaluate().center(property_map);
    }
    
    return this->updatePoint( aabox.center() );
}

/** Update the molecules used to create this point */
bool Center::update(const MoleculeData &moldata)
{
    if (mols.contains(moldata.number()))
    {
        if (mols.update(moldata))
        {
            return this->recalculatePoint();
        }
    }
    
    return false;
}
                    
/** Update the molecules used to create this point */
bool Center::update(const Molecules &molecules)
{
    if (not mols.update(molecules).isEmpty())
    {
        return this->recalculatePoint();
    }
    
    return false;
}
                    
/** Update the molecules used to create this point */
bool Center::update(const MoleculeGroup &molgroup)
{
    return this->update(molgroup.molecules());
}

/** Update the molecules used to create this point */
bool Center::update(const MolGroupsBase &molgroups)
{
    return this->update(molgroups.molecules());
}

/** Return whether or not the passed molecule would change this point */
bool Center::wouldUpdate(const MoleculeData &moldata) const
{
    if (mols.contains(moldata.number()))
    {
        Center new_cent(*this);
        return new_cent.update(moldata);
    }
    
    return false;
}
                    
/** Return whether or not the passed molecules would change this point */
bool Center::wouldUpdate(const Molecules &molecules) const
{
    if (mols.intersects(molecules))
    {
        Center new_cent(*this);
        return new_cent.update(molecules);
    }
    
    return false;
}
                    
/** Return whether or not the passed molecules would change this point */
bool Center::wouldUpdate(const MoleculeGroup &molgroup) const
{
    return this->wouldUpdate(molgroup.molecules());
}

/** Return whether or not the passed molecules would change this point */
bool Center::wouldUpdate(const MolGroupsBase &molgroups) const
{
    return this->wouldUpdate(molgroups.molecules());
}

/** Return all of the molecules used to generate this point */
Molecules Center::molecules() const
{
    return mols;
}

/** Return the number of molecules needed to generate this point */
int Center::nMolecules() const
{
    return mols.nMolecules();
}

/** Return whether or not the molecule with number 'molnum' is
    needed to generate this point */
bool Center::contains(MolNum molnum) const
{
    return mols.contains(molnum);
}

/** Return whether or not this molecule with ID 'molid' is 
    needed to generate this point */
bool Center::contains(const MolID &molid) const
{
    try
    {
        return not molid.map(mols).isEmpty();
    }
    catch(...)
    {
        return false;
    }
}

/** Return whether or not this point uses data from any of the 
    molecules in the passed forcetable */
bool Center::usesMoleculesIn(const ForceTable &forcetable) const
{
    for (Molecules::const_iterator it = mols.constBegin();
         it != mols.constEnd();
         ++it)
    {
        if (forcetable.containsTable(it->number()))
            return true;
    }
    
    return false;
}

/** Return whether or not this point uses data from any of the 
    molecules in 'molecules' */
bool Center::usesMoleculesIn(const Molecules &molecules) const
{
    for (Molecules::const_iterator it = mols.constBegin();
         it != mols.constEnd();
         ++it)
    {
        if (molecules.contains(it->number()))
            return true;
    }
    
    return false;
}

/** Return whether or not this point uses data from any of the
    molecules in the group 'molgroup' */
bool Center::usesMoleculesIn(const MoleculeGroup &molgroup) const
{
    return this->usesMoleculesIn(molgroup.molecules());
}

/** Return whether or not this point uses data from any of the
    molecules in the groups in 'molgroups' */
bool Center::usesMoleculesIn(const MolGroupsBase &molgroups) const
{
    for (MGIdx i(0); i<molgroups.nGroups(); ++i)
    {
        if (this->usesMoleculesIn(molgroups[i]))
            return true;
    }
    
    return false;
}

/** Decompose the force 'force' acting on this point from the
    molecule whose forces are in 'molforces' and add the
    force onto the table */
bool Center::addForce(MolForceTable &molforces, const Vector &force) const
{
    throw SireError::incomplete_code( QObject::tr(
            "Need to work out how to decompose forces..."), CODELOC );
            
    return false;
} 

/** Decompose the force 'force' into the forces acting on 
    the molecules that contribute to this point and add those
    forces onto the table 'forces' */
bool Center::addForce(ForceTable &forces, const Vector &force) const
{
    throw SireError::incomplete_code( QObject::tr(
            "Need to work out how to decompose forces..."), CODELOC );
            
    return false;
}

/** Return whether this is an intramolecular point (it depends on coordinates
    of atoms in just one molecule) */
bool Center::isIntraMoleculePoint() const
{
    return mols.nMolecules() == 1;
}

/** Return whether or not this is an intermolecular point (it depends on
    coordinates of atoms from than one molecule) */
bool Center::isInterMoleculePoint() const
{
    return mols.nMolecules() == 2;
}

/** Return whether or not this is an extramolecular point (it is independent
    of the coordinates of atoms in any molecule, i.e. it is just a point in space) */
bool Center::isExtraMoleculePoint() const
{
    return mols.isEmpty();
}

//////////////
////////////// Implementation of CenterOfGeometry
//////////////

static const RegisterMetaType<CenterOfGeometry> r_cog;

/** Serialise to a binary datastream */
QDataStream SIREFF_EXPORT &operator<<(QDataStream &ds, const CenterOfGeometry &cog)
{
    writeHeader(ds, r_cog, 1);
    
    SharedDataStream sds(ds);
    
    sds << cog.mols << cog.property_map
        << static_cast<const Point&>(cog);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREFF_EXPORT &operator>>(QDataStream &ds, CenterOfGeometry &cog)
{
    VersionID v = readHeader(ds, r_cog);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> cog.mols >> cog.property_map
            >> static_cast<Point&>(cog);
    }
    else
        throw version_error( v, "1", r_cog, CODELOC );
        
    return ds;
}

/** Constructor */
CenterOfGeometry::CenterOfGeometry() : ConcreteProperty<CenterOfGeometry,Point>()
{}

/** Construct to get the center of the molecule view 'molview' using the
    passed property map to find the required properties 
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
CenterOfGeometry::CenterOfGeometry(const MoleculeView &molview, const PropertyMap &map)
                 : ConcreteProperty<CenterOfGeometry,Point>(), 
                   mols(molview), property_map(map)
{
    this->recalculatePoint();
}

/** Construct to get the center of the molecules in 'molecules', using the
    passed property map to find the required properties 
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
CenterOfGeometry::CenterOfGeometry(const Molecules &molecules, const PropertyMap &map)
                 : ConcreteProperty<CenterOfGeometry,Point>(), 
                   mols(molecules), property_map(map)
{
    this->recalculatePoint();
}

/** Copy constructor */
CenterOfGeometry::CenterOfGeometry(const CenterOfGeometry &other)
                 : ConcreteProperty<CenterOfGeometry,Point>(other), 
                   mols(other.mols), property_map(other.property_map)
{}

/** Destructor */
CenterOfGeometry::~CenterOfGeometry()
{}

/** Copy assignment operator */
CenterOfGeometry& CenterOfGeometry::operator=(const CenterOfGeometry &other)
{
    mols = other.mols;
    property_map = other.property_map;
    Point::operator=(other);
    
    return *this;
}

/** Comparison operator */
bool CenterOfGeometry::operator==(const CenterOfGeometry &other) const
{
    return mols == other.mols and property_map == other.property_map and
           Point::operator==(other);
}

/** Comparison operator */
bool CenterOfGeometry::operator!=(const CenterOfGeometry &other) const
{
    return mols != other.mols or property_map != other.property_map or
           Point::operator!=(other);
}

const char* CenterOfGeometry::typeName()
{
    return QMetaType::typeName( qMetaTypeId<CenterOfGeometry>() );
}

/** Return a string representation */
QString CenterOfGeometry::toString() const
{
    return QString("CenterOfGeometry{ %1 : %2 }").arg(mols.toString())
                                                 .arg(this->point().toString());
}

static Vector averagePoints(const QList< tuple<Vector,double> > &points)
{
    double sum_of_weights = 0;
    
    for (QList< tuple<Vector,double> >::const_iterator it = points.constBegin();
         it != points.constEnd();
         ++it)
    {
        sum_of_weights += it->get<1>();
    }
    
    sum_of_weights = 1.0 / sum_of_weights;
    
    Vector point;
    
    for (QList< tuple<Vector,double> >::const_iterator it = points.constBegin();
         it != points.constEnd();
         ++it)
    {
        point += (sum_of_weights * it->get<1>()) * it->get<0>();
    }
    
    return point;
}

/** Set the space used by this point - a CenterOfGeometry cannot
    be calculated for periodic or non-cartesian spaces if there
    is more than one molecule */
void CenterOfGeometry::setSpace(const Space &space)
{
    if (mols.nMolecules() > 1)
    {
        if (space.isPeriodic())
            throw SireVol::incompatible_space( QObject::tr(
                "A CenterOfGeometry point with more than one molecule cannot be used "
                "with the periodic space \"%1\". Either use a CenterOfGeometry point "
                "built from a single molecule, or switch to a Cartesian space.")
                    .arg(space.what()), CODELOC );
    
        if (not space.isCartesian())
            throw SireVol::incompatible_space( QObject::tr(
               "A CenterOfGeometry point with more than one molecule cannot be used "
               "with the non-cartesian space \"%1\". Either use a CenterOfGeometry point "
               "built from a single molecule, or switch to a Cartesian space.")
                   .arg(space.what()), CODELOC );
    }
    
    Point::setSpace(space);
}

/** Internal function used to recalculate the center */
bool CenterOfGeometry::recalculatePoint()
{
    if (mols.isEmpty())
        return false;

    else if (mols.count() == 1)
    {
        return Point::updatePoint( mols.constBegin()->evaluate()
                                       .centerOfGeometry(property_map) );
    }

    QList< tuple<Vector,double> > points;

    for (Molecules::const_iterator it = mols.constBegin();
         it != mols.constEnd();
         ++it)
    {
        points.append( tuple<Vector,double>(it->evaluate().centerOfGeometry(property_map),
                                            it->selection().nAtoms()) );
    }
    
    return Point::updatePoint( ::averagePoints(points) );
}

/** Update the molecules used to create this point */
bool CenterOfGeometry::update(const MoleculeData &moldata)
{
    if (mols.contains(moldata.number()))
    {
        if (mols.update(moldata))
            return this->recalculatePoint();
    }
    
    return false;
}
                    
/** Update the molecules used to create this point */
bool CenterOfGeometry::update(const Molecules &molecules)
{
    if (not mols.update(molecules).isEmpty())
    {
        return this->recalculatePoint();
    }
    
    return false;
}
                    
/** Update the molecules used to create this point */
bool CenterOfGeometry::update(const MoleculeGroup &molgroup)
{
    return this->update(molgroup.molecules());
}

/** Update the molecules used to create this point */
bool CenterOfGeometry::update(const MolGroupsBase &molgroups)
{
    return this->update(molgroups.molecules());
}

/** Return whether or not the passed molecule would change this point */
bool CenterOfGeometry::wouldUpdate(const MoleculeData &moldata) const
{
    if (mols.contains(moldata.number()))
    {
        CenterOfGeometry new_cent(*this);
        return new_cent.update(moldata);
    }
    
    return false;
}
                    
/** Return whether or not the passed molecules would change this point */
bool CenterOfGeometry::wouldUpdate(const Molecules &molecules) const
{
    if (mols.intersects(molecules))
    {
        CenterOfGeometry new_cent(*this);
        return new_cent.update(molecules);
    }
    
    return false;
}
                    
/** Return whether or not the passed molecules would change this point */
bool CenterOfGeometry::wouldUpdate(const MoleculeGroup &molgroup) const
{
    return this->wouldUpdate(molgroup.molecules());
}

/** Return whether or not the passed molecules would change this point */
bool CenterOfGeometry::wouldUpdate(const MolGroupsBase &molgroups) const
{
    return this->wouldUpdate(molgroups.molecules());
}

/** Return all of the molecules used to generate this point */
Molecules CenterOfGeometry::molecules() const
{
    return mols;
}

/** Return the number of molecules needed to generate this point */
int CenterOfGeometry::nMolecules() const
{
    return mols.nMolecules();
}

/** Return whether or not the molecule with number 'molnum' is
    needed to generate this point */
bool CenterOfGeometry::contains(MolNum molnum) const
{
    return mols.contains(molnum);
}

/** Return whether or not this molecule with ID 'molid' is 
    needed to generate this point */
bool CenterOfGeometry::contains(const MolID &molid) const
{
    try
    {
        return not molid.map(mols).isEmpty();
    }
    catch(...)
    {
        return false;
    }
}

/** Return whether or not this point uses data from any of the 
    molecules in the passed forcetable */
bool CenterOfGeometry::usesMoleculesIn(const ForceTable &forcetable) const
{
    for (Molecules::const_iterator it = mols.constBegin();
         it != mols.constEnd();
         ++it)
    {
        if (forcetable.containsTable(it->number()))
            return true;
    }
    
    return false;
}

/** Return whether or not this point uses data from any of the 
    molecules in 'molecules' */
bool CenterOfGeometry::usesMoleculesIn(const Molecules &molecules) const
{
    for (Molecules::const_iterator it = mols.constBegin();
         it != mols.constEnd();
         ++it)
    {
        if (molecules.contains(it->number()))
            return true;
    }
    
    return false;
}

/** Return whether or not this point uses data from any of the
    molecules in the group 'molgroup' */
bool CenterOfGeometry::usesMoleculesIn(const MoleculeGroup &molgroup) const
{
    return this->usesMoleculesIn(molgroup.molecules());
}

/** Return whether or not this point uses data from any of the
    molecules in the groups in 'molgroups' */
bool CenterOfGeometry::usesMoleculesIn(const MolGroupsBase &molgroups) const
{
    for (MGIdx i(0); i<molgroups.nGroups(); ++i)
    {
        if (this->usesMoleculesIn(molgroups[i]))
            return true;
    }
    
    return false;
}

/** Decompose the force 'force' acting on this point from the
    molecule whose forces are in 'molforces' and add the
    force onto the table */
bool CenterOfGeometry::addForce(MolForceTable &molforces, const Vector &force) const
{
    throw SireError::incomplete_code( QObject::tr(
            "Need to work out how to decompose forces..."), CODELOC );
            
    return false;
} 

/** Decompose the force 'force' into the forces acting on 
    the molecules that contribute to this point and add those
    forces onto the table 'forces' */
bool CenterOfGeometry::addForce(ForceTable &forces, const Vector &force) const
{
    throw SireError::incomplete_code( QObject::tr(
            "Need to work out how to decompose forces..."), CODELOC );
            
    return false;
}

/** Return whether this is an intramolecular point (it depends on coordinates
    of atoms in just one molecule) */
bool CenterOfGeometry::isIntraMoleculePoint() const
{
    return mols.nMolecules() == 1;
}

/** Return whether or not this is an intermolecular point (it depends on
    coordinates of atoms from than one molecule) */
bool CenterOfGeometry::isInterMoleculePoint() const
{
    return mols.nMolecules() == 2;
}

/** Return whether or not this is an extramolecular point (it is independent
    of the coordinates of atoms in any molecule, i.e. it is just a point in space) */
bool CenterOfGeometry::isExtraMoleculePoint() const
{
    return mols.isEmpty();
}

//////////////
////////////// Implementation of CenterOfMass
//////////////

static const RegisterMetaType<CenterOfMass> r_com;

/** Serialise to a binary datastream */
QDataStream SIREFF_EXPORT &operator<<(QDataStream &ds, const CenterOfMass &com)
{
    writeHeader(ds, r_com, 1);
    
    SharedDataStream sds(ds);
    
    sds << com.mols << com.property_map
        << static_cast<const Point&>(com);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREFF_EXPORT &operator>>(QDataStream &ds, CenterOfMass &com)
{
    VersionID v = readHeader(ds, r_com);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> com.mols >> com.property_map
            >> static_cast<Point&>(com);
    }
    else
        throw version_error( v, "1", r_com, CODELOC );
        
    return ds;
}

/** Constructor */
CenterOfMass::CenterOfMass() : ConcreteProperty<CenterOfMass,Point>()
{}

/** Construct to get the center of the molecule view 'molview' using the
    passed property map to find the required properties 
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
CenterOfMass::CenterOfMass(const MoleculeView &molview, const PropertyMap &map)
                 : ConcreteProperty<CenterOfMass,Point>(), 
                   mols(molview), property_map(map)
{
    this->recalculatePoint();
}

/** Construct to get the center of the molecules in 'molecules', using the
    passed property map to find the required properties 
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
CenterOfMass::CenterOfMass(const Molecules &molecules, const PropertyMap &map)
                 : ConcreteProperty<CenterOfMass,Point>(), 
                   mols(molecules), property_map(map)
{
    this->recalculatePoint();
}

/** Copy constructor */
CenterOfMass::CenterOfMass(const CenterOfMass &other)
                 : ConcreteProperty<CenterOfMass,Point>(other), 
                   mols(other.mols), property_map(other.property_map)
{}

/** Destructor */
CenterOfMass::~CenterOfMass()
{}

/** Copy assignment operator */
CenterOfMass& CenterOfMass::operator=(const CenterOfMass &other)
{
    mols = other.mols;
    property_map = other.property_map;
    Point::operator=(other);
    
    return *this;
}

/** Comparison operator */
bool CenterOfMass::operator==(const CenterOfMass &other) const
{
    return mols == other.mols and property_map == other.property_map and
           Point::operator==(other);
}

/** Comparison operator */
bool CenterOfMass::operator!=(const CenterOfMass &other) const
{
    return mols != other.mols or property_map != other.property_map or
           Point::operator!=(other);
}

const char* CenterOfMass::typeName()
{
    return QMetaType::typeName( qMetaTypeId<CenterOfMass>() );
}

/** Return a string representation */
QString CenterOfMass::toString() const
{
    return QString("CenterOfMass{ %1 : %2 }").arg(mols.toString())
                                             .arg(this->point().toString());
}

/** Set the space used by this point - a CenterOfGeometry cannot
    be calculated for periodic or non-cartesian spaces if there
    is more than one molecule */
void CenterOfMass::setSpace(const Space &space)
{
    if (mols.nMolecules() > 1)
    {
        if (space.isPeriodic())
            throw SireVol::incompatible_space( QObject::tr(
                "A CenterOfMass point with more than one molecule cannot be used "
                "with the periodic space \"%1\". Either use a CenterOfMass point "
                "built from a single molecule, or switch to a Cartesian space.")
                    .arg(space.what()), CODELOC );
    
        if (not space.isCartesian())
            throw SireVol::incompatible_space( QObject::tr(
                "A CenterOfMass point with more than one molecule cannot be used "
                "with the non-cartesian space \"%1\". Either use a CenterOfMass point "
                "built from a single molecule, or switch to a Cartesian space.")
                    .arg(space.what()), CODELOC );
    }
    
    Point::setSpace(space);
}

/** Internal function used to recalculate the center */
bool CenterOfMass::recalculatePoint()
{
    if (mols.isEmpty())
        return false;

    else if (mols.count() == 1)
    {
        return Point::updatePoint( mols.constBegin()->evaluate()
                                       .centerOfMass(property_map) );
    }

    QList< tuple<Vector,double> > points;

    for (Molecules::const_iterator it = mols.constBegin();
         it != mols.constEnd();
         ++it)
    {
        Evaluator evaluator(*it);
    
        points.append( tuple<Vector,double>( evaluator.centerOfMass(property_map),
                                             evaluator.mass() ) );
    }
    
    return Point::updatePoint( ::averagePoints(points) );
}

/** Update the molecules used to create this point */
bool CenterOfMass::update(const MoleculeData &moldata)
{
    if (mols.contains(moldata.number()))
    {
        if (mols.update(moldata))
            return this->recalculatePoint();
    }
    
    return false;
}
                    
/** Update the molecules used to create this point */
bool CenterOfMass::update(const Molecules &molecules)
{
    if (not mols.update(molecules).isEmpty())
    {
        return this->recalculatePoint();
    }
    
    return false;
}
                    
/** Update the molecules used to create this point */
bool CenterOfMass::update(const MoleculeGroup &molgroup)
{
    return this->update(molgroup.molecules());
}

/** Update the molecules used to create this point */
bool CenterOfMass::update(const MolGroupsBase &molgroups)
{
    return this->update(molgroups.molecules());
}

/** Return whether or not the passed molecule would change this point */
bool CenterOfMass::wouldUpdate(const MoleculeData &moldata) const
{
    if (mols.contains(moldata.number()))
    {
        CenterOfMass new_cent(*this);
        return new_cent.update(moldata);
    }
    
    return false;
}
                    
/** Return whether or not the passed molecules would change this point */
bool CenterOfMass::wouldUpdate(const Molecules &molecules) const
{
    if (mols.intersects(molecules))
    {
        CenterOfMass new_cent(*this);
        return new_cent.update(molecules);
    }
    
    return false;
}
                    
/** Return whether or not the passed molecules would change this point */
bool CenterOfMass::wouldUpdate(const MoleculeGroup &molgroup) const
{
    return this->wouldUpdate(molgroup.molecules());
}

/** Return whether or not the passed molecules would change this point */
bool CenterOfMass::wouldUpdate(const MolGroupsBase &molgroups) const
{
    return this->wouldUpdate(molgroups.molecules());
}

/** Return all of the molecules used to generate this point */
Molecules CenterOfMass::molecules() const
{
    return mols;
}

/** Return the number of molecules needed to generate this point */
int CenterOfMass::nMolecules() const
{
    return mols.nMolecules();
}

/** Return whether or not the molecule with number 'molnum' is
    needed to generate this point */
bool CenterOfMass::contains(MolNum molnum) const
{
    return mols.contains(molnum);
}

/** Return whether or not this molecule with ID 'molid' is 
    needed to generate this point */
bool CenterOfMass::contains(const MolID &molid) const
{
    try
    {
        return not molid.map(mols).isEmpty();
    }
    catch(...)
    {
        return false;
    }
}

/** Return whether or not this point uses data from any of the 
    molecules in the passed forcetable */
bool CenterOfMass::usesMoleculesIn(const ForceTable &forcetable) const
{
    for (Molecules::const_iterator it = mols.constBegin();
         it != mols.constEnd();
         ++it)
    {
        if (forcetable.containsTable(it->number()))
            return true;
    }
    
    return false;
}

/** Return whether or not this point uses data from any of the 
    molecules in 'molecules' */
bool CenterOfMass::usesMoleculesIn(const Molecules &molecules) const
{
    for (Molecules::const_iterator it = mols.constBegin();
         it != mols.constEnd();
         ++it)
    {
        if (molecules.contains(it->number()))
            return true;
    }
    
    return false;
}

/** Return whether or not this point uses data from any of the
    molecules in the group 'molgroup' */
bool CenterOfMass::usesMoleculesIn(const MoleculeGroup &molgroup) const
{
    return this->usesMoleculesIn(molgroup.molecules());
}

/** Return whether or not this point uses data from any of the
    molecules in the groups in 'molgroups' */
bool CenterOfMass::usesMoleculesIn(const MolGroupsBase &molgroups) const
{
    for (MGIdx i(0); i<molgroups.nGroups(); ++i)
    {
        if (this->usesMoleculesIn(molgroups[i]))
            return true;
    }
    
    return false;
}

/** Decompose the force 'force' acting on this point from the
    molecule whose forces are in 'molforces' and add the
    force onto the table */
bool CenterOfMass::addForce(MolForceTable &molforces, const Vector &force) const
{
    throw SireError::incomplete_code( QObject::tr(
            "Need to work out how to decompose forces..."), CODELOC );
            
    return false;
} 

/** Decompose the force 'force' into the forces acting on 
    the molecules that contribute to this point and add those
    forces onto the table 'forces' */
bool CenterOfMass::addForce(ForceTable &forces, const Vector &force) const
{
    throw SireError::incomplete_code( QObject::tr(
            "Need to work out how to decompose forces..."), CODELOC );
            
    return false;
}

/** Return whether this is an intramolecular point (it depends on coordinates
    of atoms in just one molecule) */
bool CenterOfMass::isIntraMoleculePoint() const
{
    return mols.nMolecules() == 1;
}

/** Return whether or not this is an intermolecular point (it depends on
    coordinates of atoms from than one molecule) */
bool CenterOfMass::isInterMoleculePoint() const
{
    return mols.nMolecules() == 2;
}

/** Return whether or not this is an extramolecular point (it is independent
    of the coordinates of atoms in any molecule, i.e. it is just a point in space) */
bool CenterOfMass::isExtraMoleculePoint() const
{
    return mols.isEmpty();
}

//////////////
////////////// Implementation of PointRef
//////////////

/** Construct from the passed atom */
PointRef::PointRef(const Atom &atom) : ptr( AtomPoint(atom) )
{}

/** Construct from the passed point */
PointRef::PointRef(const Vector &point) : ptr( VectorPoint(point) )
{}

/** Construct from the passed point */
PointRef::PointRef(const Point &point) : ptr(point)
{}

/** Construct from the passed point */
PointRef::PointRef(const PointPtr &point) : ptr(point)
{}

/** Destructor */
PointRef::~PointRef()
{}

/** Return the point in 3D space */
const Vector& PointRef::operator()() const
{
    return ptr.read().operator()();
}

/** Return the point in 3D space */
const Vector& PointRef::point() const
{
    return ptr.read().point();
}

/** Allow automatic casting to a Point */
PointRef::operator const Point&() const
{
    return ptr.read();
}

/** Decompose the force 'force' acting on this point from the
    molecule whose forces are in 'molforces' and add the
    force onto the table */
bool PointRef::addForce(MolForceTable &molforces, const Vector &force) const
{
    return ptr.read().addForce(molforces, force);
} 

/** Decompose the force 'force' into the forces acting on 
    the molecules that contribute to this point and add those
    forces onto the table 'forces' */
bool PointRef::addForce(ForceTable &forces, const Vector &force) const
{
    return ptr.read().addForce(forces, force);
}
