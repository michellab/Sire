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

#ifndef SIREFF_POINT_H
#define SIREFF_POINT_H

#include "SireBase/property.h"

#include "SireMol/atom.h"
#include "SireMol/molecules.h"
#include "SireMol/selector.hpp"
#include "SireMol/mover.hpp"

#include "SireVol/space.h"

SIRE_BEGIN_HEADER

namespace SireFF
{

class Point;
class PointBase;

class AtomPoint;
class VectorPoint;
class Center;
class CenterOfMass;
class CenterOfGeometry;

}

SIREFF_EXPORT QDataStream& operator<<(QDataStream&, const SireFF::Point&);
SIREFF_EXPORT QDataStream& operator>>(QDataStream&, SireFF::Point&);

SIREFF_EXPORT QDataStream& operator<<(QDataStream&, const SireFF::AtomPoint&);
SIREFF_EXPORT QDataStream& operator>>(QDataStream&, SireFF::AtomPoint&);

SIREFF_EXPORT QDataStream& operator<<(QDataStream&, const SireFF::VectorPoint&);
SIREFF_EXPORT QDataStream& operator>>(QDataStream&, SireFF::VectorPoint&);

SIREFF_EXPORT QDataStream& operator<<(QDataStream&, const SireFF::Center&);
SIREFF_EXPORT QDataStream& operator>>(QDataStream&, SireFF::Center&);

SIREFF_EXPORT QDataStream& operator<<(QDataStream&, const SireFF::CenterOfGeometry&);
SIREFF_EXPORT QDataStream& operator>>(QDataStream&, SireFF::CenterOfGeometry&);

SIREFF_EXPORT QDataStream& operator<<(QDataStream&, const SireFF::CenterOfMass&);
SIREFF_EXPORT QDataStream& operator>>(QDataStream&, SireFF::CenterOfMass&);

namespace SireFF
{
class MolForceTable;
class ForceTable;
}

namespace SireFF
{

using SireMol::Atom;
using SireMol::MoleculeData;
using SireMol::MoleculeView;
using SireMol::Molecules;
using SireMol::MoleculeGroup;
using SireMol::MolGroupsBase;
using SireMol::MolNum;
using SireMol::MolID;

using SireBase::PropertyMap;

using SireVol::Space;

using SireMaths::Vector;

class Point;
typedef SireBase::PropPtr<Point> PointPtr;

/** This is the base class of all Points. A Point is a class that 
    allows a view of a molecule (or molecules) to be turned
    into a 3D point. This 3D point can then be used within a restraint
    class (thereby seperating out the code that implements the 
    restraint from the code that selects the parts of the molecule(s)
    being restrained)
    
    @author Christopher Woods
*/
class SIREFF_EXPORT Point : public SireBase::Property
{

friend SIREFF_EXPORT QDataStream& ::operator<<(QDataStream&, const Point&);
friend SIREFF_EXPORT QDataStream& ::operator>>(QDataStream&, Point&);

public:
    Point();
    
    virtual ~Point();

    virtual Point* clone() const=0;
    
    const Vector& operator()() const;
    
    static const char* typeName()
    {
        return "SireFF::Point";
    }
    
    const Vector& point() const;
    
    const Space& space() const;
    
    virtual void setSpace(const Space &space);
    
    virtual bool update(const MoleculeData &moldata)=0;
    virtual bool update(const Molecules &molecules)=0;
    virtual bool update(const MoleculeGroup &molgroup)=0;
    virtual bool update(const MolGroupsBase &molgroups)=0;
    
    virtual bool wouldUpdate(const MoleculeData &moldata) const=0;
    virtual bool wouldUpdate(const Molecules &molecules) const=0;
    virtual bool wouldUpdate(const MoleculeGroup &molgroup) const=0;
    virtual bool wouldUpdate(const MolGroupsBase &molgroups) const=0;
    
    virtual Molecules molecules() const=0;
    
    virtual int nMolecules() const=0;
    
    virtual bool contains(MolNum molnum) const=0;
    virtual bool contains(const MolID &molid) const=0;
    
    virtual bool usesMoleculesIn(const ForceTable &forcetable) const=0;
    virtual bool usesMoleculesIn(const Molecules &molecules) const=0;
    virtual bool usesMoleculesIn(const MoleculeGroup &molgroup) const=0;
    virtual bool usesMoleculesIn(const MolGroupsBase &molgroups) const=0;
    
    virtual bool addForce(MolForceTable &molforces, const Vector &force) const=0;
    virtual bool addForce(ForceTable &forces, const Vector &force) const=0;
    
    static const VectorPoint& null();

    virtual bool isExtraMoleculePoint() const=0;
    virtual bool isIntraMoleculePoint() const=0;
    virtual bool isInterMoleculePoint() const=0;

    static bool areIntraMoleculePoints(const Point &point0, const Point &point1);
    static bool areIntraMoleculePoints(const Point &point0, const Point &point1,
                                       const Point &point2);
    static bool areIntraMoleculePoints(const Point &point0, const Point &point1,
                                       const Point &point2, const Point &point3);
                                       
    static bool areIntraMoleculePoints(const QVector<PointPtr> &points);
    
protected:
    Point(const Vector &point);
    Point(const Point &other);
    
    Point& operator=(const Point &other);
    
    bool operator==(const Point &other) const;
    bool operator!=(const Point &other) const;

    bool updatePoint(const Vector &point);

private:
    /** The actual point in space */
    Vector p;
    
    /** The 3D space in which this point is calculated and exists */
    SireVol::SpacePtr spce;
};

/** This is a small class used to help convert different molecule
    views into a Point class in function signatures, e.g.
    
    DistanceRestraint(const PointRef &point0, const PointRef &point1);
    
    would allow distance restraints to be created between pairs
    of atoms, or an atom and a vector etc.
    
    @author Christopher Woods
*/
class SIREFF_EXPORT PointRef
{
public:
    PointRef(const Atom &atom);
    PointRef(const Vector &point);
    PointRef(const Point &point);
    PointRef(const PointPtr &point);
    
    ~PointRef();
    
    const Vector& operator()() const;
    
    const Vector& point() const;
    
    operator const Point&() const;
    
    bool addForce(MolForceTable &molforces, const Vector &force) const;
    bool addForce(ForceTable &forces, const Vector &force) const;
    
private:
    /** Pointer to the implementation of this point */
    PointPtr ptr;
};

/** This point returns the location of an atom */
class SIREFF_EXPORT AtomPoint : public SireBase::ConcreteProperty<AtomPoint,Point>
{

friend SIREFF_EXPORT QDataStream& ::operator<<(QDataStream&, const AtomPoint&);
friend SIREFF_EXPORT QDataStream& ::operator>>(QDataStream&, AtomPoint&);

public:
    AtomPoint();
    AtomPoint(const Atom &atom, const PropertyMap &map=PropertyMap());
    
    AtomPoint(const AtomPoint &other);
    
    ~AtomPoint();
    
    AtomPoint& operator=(const AtomPoint &other);
    
    bool operator==(const AtomPoint &other) const;
    bool operator!=(const AtomPoint &other) const;
    
    static const char* typeName();
    
    QString toString() const;
    
    bool update(const MoleculeData &moldata);
    bool update(const Molecules &molecules);
    bool update(const MoleculeGroup &molgroup);
    bool update(const MolGroupsBase &molgroups);
    
    bool wouldUpdate(const MoleculeData &moldata) const;
    bool wouldUpdate(const Molecules &molecules) const;
    bool wouldUpdate(const MoleculeGroup &molgroup) const;
    bool wouldUpdate(const MolGroupsBase &molgroups) const;
    
    Molecules molecules() const;
    
    int nMolecules() const;
    
    bool contains(MolNum molnum) const;
    bool contains(const MolID &molid) const;

    bool usesMoleculesIn(const ForceTable &forcetable) const;
    bool usesMoleculesIn(const Molecules &molecules) const;
    bool usesMoleculesIn(const MoleculeGroup &molgroup) const;
    bool usesMoleculesIn(const MolGroupsBase &molgroups) const;

    const Atom& atom() const;
    
    bool addForce(MolForceTable &molforces, const Vector &force) const;
    bool addForce(ForceTable &forces, const Vector &force) const;

    bool isExtraMoleculePoint() const;
    bool isIntraMoleculePoint() const;
    bool isInterMoleculePoint() const;

private:
    /** The actual atom! */
    Atom atm;
    
    /** The coordinates property */
    SireBase::PropertyName coords_property;
};

/** This is a simple wrapper for a point in space */
class SIREFF_EXPORT VectorPoint : public SireBase::ConcreteProperty<VectorPoint,Point>
{

friend SIREFF_EXPORT QDataStream& ::operator<<(QDataStream&, const VectorPoint&);
friend SIREFF_EXPORT QDataStream& ::operator>>(QDataStream&, VectorPoint&);

public:
    VectorPoint();
    VectorPoint(const Vector &point);
    
    VectorPoint(const VectorPoint &other);
    
    ~VectorPoint();
    
    VectorPoint& operator=(const VectorPoint &other);
    
    bool operator==(const VectorPoint &other) const;
    bool operator!=(const VectorPoint &other) const;
    
    static const char* typeName();
    
    QString toString() const;
    
    bool update(const MoleculeData &moldata);
    bool update(const Molecules &molecules);
    bool update(const MoleculeGroup &molgroup);
    bool update(const MolGroupsBase &molgroups);
    
    bool wouldUpdate(const MoleculeData &moldata) const;
    bool wouldUpdate(const Molecules &molecules) const;
    bool wouldUpdate(const MoleculeGroup &molgroup) const;
    bool wouldUpdate(const MolGroupsBase &molgroups) const;
    
    Molecules molecules() const;
    
    int nMolecules() const;
    
    bool contains(MolNum molnum) const;
    bool contains(const MolID &molid) const;

    bool usesMoleculesIn(const ForceTable &forcetable) const;
    bool usesMoleculesIn(const Molecules &molecules) const;
    bool usesMoleculesIn(const MoleculeGroup &molgroup) const;
    bool usesMoleculesIn(const MolGroupsBase &molgroups) const;
    
    bool addForce(MolForceTable &molforces, const Vector &force) const;
    bool addForce(ForceTable &forces, const Vector &force) const;

    bool isExtraMoleculePoint() const;
    bool isIntraMoleculePoint() const;
    bool isInterMoleculePoint() const;
};

/** This point returns the center of a view of a molecule, or group
    of molecules */
class SIREFF_EXPORT Center : public SireBase::ConcreteProperty<Center,Point>
{

friend SIREFF_EXPORT QDataStream& ::operator<<(QDataStream&, const Center&);
friend SIREFF_EXPORT QDataStream& ::operator>>(QDataStream&, Center&);

public:
    Center();
    Center(const MoleculeView &molview, const PropertyMap &map = PropertyMap());
    Center(const Molecules &molecules, const PropertyMap &map = PropertyMap());
    
    Center(const Center &other);
    
    ~Center();
    
    Center& operator=(const Center &other);
    
    bool operator==(const Center &other) const;
    bool operator!=(const Center &other) const;
    
    static const char* typeName();

    QString toString() const;

    void setSpace(const Space &space);
    
    bool update(const MoleculeData &moldata);
    bool update(const Molecules &molecules);
    bool update(const MoleculeGroup &molgroup);
    bool update(const MolGroupsBase &molgroups);
    
    bool wouldUpdate(const MoleculeData &moldata) const;
    bool wouldUpdate(const Molecules &molecules) const;
    bool wouldUpdate(const MoleculeGroup &molgroup) const;
    bool wouldUpdate(const MolGroupsBase &molgroups) const;
    
    Molecules molecules() const;
    
    int nMolecules() const;
    
    bool contains(MolNum molnum) const;
    bool contains(const MolID &molid) const;

    bool usesMoleculesIn(const ForceTable &forcetable) const;
    bool usesMoleculesIn(const Molecules &molecules) const;
    bool usesMoleculesIn(const MoleculeGroup &molgroup) const;
    bool usesMoleculesIn(const MolGroupsBase &molgroups) const;
    
    bool addForce(MolForceTable &molforces, const Vector &force) const;
    bool addForce(ForceTable &forces, const Vector &force) const;

    bool isExtraMoleculePoint() const;
    bool isIntraMoleculePoint() const;
    bool isInterMoleculePoint() const;

private:
    bool recalculatePoint();

    /** The molecules whose center is to be determined */
    Molecules mols;

    /** The map used to find the properties */
    PropertyMap property_map;
};

/** This point returns the center of geometry of a view of a molecule,
    or group of molecules */
class SIREFF_EXPORT CenterOfGeometry 
            : public SireBase::ConcreteProperty<CenterOfGeometry,Point>
{

friend SIREFF_EXPORT QDataStream& ::operator<<(QDataStream&, const CenterOfGeometry&);
friend SIREFF_EXPORT QDataStream& ::operator>>(QDataStream&, CenterOfGeometry&);

public:
    CenterOfGeometry();
    CenterOfGeometry(const MoleculeView &molview,
                     const PropertyMap &map = PropertyMap());
    CenterOfGeometry(const Molecules &molecules,
                     const PropertyMap &map = PropertyMap());
                     
    CenterOfGeometry(const CenterOfGeometry &other);
    
    ~CenterOfGeometry();
    
    CenterOfGeometry& operator=(const CenterOfGeometry &other);

    bool operator==(const CenterOfGeometry &other) const;
    bool operator!=(const CenterOfGeometry &other) const;

    static const char* typeName();

    QString toString() const;

    void setSpace(const Space &space);
    
    bool update(const MoleculeData &moldata);
    bool update(const Molecules &molecules);
    bool update(const MoleculeGroup &molgroup);
    bool update(const MolGroupsBase &molgroups);
    
    bool wouldUpdate(const MoleculeData &moldata) const;
    bool wouldUpdate(const Molecules &molecules) const;
    bool wouldUpdate(const MoleculeGroup &molgroup) const;
    bool wouldUpdate(const MolGroupsBase &molgroups) const;
    
    Molecules molecules() const;
    
    int nMolecules() const;
    
    bool contains(MolNum molnum) const;
    bool contains(const MolID &molid) const;

    bool usesMoleculesIn(const ForceTable &forcetable) const;
    bool usesMoleculesIn(const Molecules &molecules) const;
    bool usesMoleculesIn(const MoleculeGroup &molgroup) const;
    bool usesMoleculesIn(const MolGroupsBase &molgroups) const;
    
    bool addForce(MolForceTable &molforces, const Vector &force) const;
    bool addForce(ForceTable &forces, const Vector &force) const;

    bool isExtraMoleculePoint() const;
    bool isIntraMoleculePoint() const;
    bool isInterMoleculePoint() const;

private:
    bool recalculatePoint();

    /** The molecules whose center of geometry is to be determined */
    Molecules mols;

    /** The map used to find the properties */
    PropertyMap property_map;
};

/** This point returns the center of mass of a view of a molecule,
    or group of molecules */
class SIREFF_EXPORT CenterOfMass 
            : public SireBase::ConcreteProperty<CenterOfMass,Point>
{

friend SIREFF_EXPORT QDataStream& ::operator<<(QDataStream&, const CenterOfMass&);
friend SIREFF_EXPORT QDataStream& ::operator>>(QDataStream&, CenterOfMass&);

public:
    CenterOfMass();
    CenterOfMass(const MoleculeView &molview,
                 const PropertyMap &map = PropertyMap());
    CenterOfMass(const Molecules &molecules,
                 const PropertyMap &map = PropertyMap());
                     
    CenterOfMass(const CenterOfMass &other);
    
    ~CenterOfMass();
    
    CenterOfMass& operator=(const CenterOfMass &other);

    bool operator==(const CenterOfMass &other) const;
    bool operator!=(const CenterOfMass &other) const;

    static const char* typeName();

    QString toString() const;

    void setSpace(const Space &space);
    
    bool update(const MoleculeData &moldata);
    bool update(const Molecules &molecules);
    bool update(const MoleculeGroup &molgroup);
    bool update(const MolGroupsBase &molgroups);
    
    bool wouldUpdate(const MoleculeData &moldata) const;
    bool wouldUpdate(const Molecules &molecules) const;
    bool wouldUpdate(const MoleculeGroup &molgroup) const;
    bool wouldUpdate(const MolGroupsBase &molgroups) const;
    
    Molecules molecules() const;
    
    int nMolecules() const;
    
    bool contains(MolNum molnum) const;
    bool contains(const MolID &molid) const;

    bool usesMoleculesIn(const ForceTable &forcetable) const;
    bool usesMoleculesIn(const Molecules &molecules) const;
    bool usesMoleculesIn(const MoleculeGroup &molgroup) const;
    bool usesMoleculesIn(const MolGroupsBase &molgroups) const;
    
    bool addForce(MolForceTable &molforces, const Vector &force) const;
    bool addForce(ForceTable &forces, const Vector &force) const;

    bool isExtraMoleculePoint() const;
    bool isIntraMoleculePoint() const;
    bool isInterMoleculePoint() const;

private:
    bool recalculatePoint();

    /** The molecules whose center of geometry is to be determined */
    Molecules mols;

    /** The map used to find the properties */
    PropertyMap property_map;
};

}

Q_DECLARE_METATYPE( SireFF::AtomPoint )
Q_DECLARE_METATYPE( SireFF::VectorPoint )
Q_DECLARE_METATYPE( SireFF::Center )
Q_DECLARE_METATYPE( SireFF::CenterOfGeometry )
Q_DECLARE_METATYPE( SireFF::CenterOfMass )

SIRE_EXPOSE_CLASS( SireFF::Point )
SIRE_EXPOSE_CLASS( SireFF::PointRef )
SIRE_EXPOSE_CLASS( SireFF::AtomPoint )
SIRE_EXPOSE_CLASS( SireFF::VectorPoint )
SIRE_EXPOSE_CLASS( SireFF::Center )
SIRE_EXPOSE_CLASS( SireFF::CenterOfGeometry )
SIRE_EXPOSE_CLASS( SireFF::CenterOfMass )

SIRE_EXPOSE_PROPERTY( SireFF::PointPtr, SireFF::Point )

SIRE_END_HEADER

#endif
