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

#include "evaluator.h"
#include "atomcoords.h"
#include "atommasses.h"
#include "atomcharges.h"
#include "atomelements.h"
#include "atommatcher.h"
#include "atommatchers.h"
#include "bondid.h"
#include "angleid.h"
#include "dihedralid.h"
#include "connectivity.h"
#include "molecule.h"
#include "mover.hpp"
#include "editor.hpp"

#include "SireVol/coordgroup.h"

#include "SireMaths/sphere.h"
#include "SireMaths/axisset.h"
#include "SireMaths/line.h"
#include "SireMaths/triangle.h"
#include "SireMaths/torsion.h"
#include "SireMaths/accumulator.h"

#include "SireBase/errors.h"
#include "SireMol/errors.h"

#include "SireUnits/dimensions.h"
#include "SireUnits/units.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"

#include <QDebug>
#include <QElapsedTimer>

using namespace SireMol;
using namespace SireMaths;
using namespace SireVol;
using namespace SireUnits;
using namespace SireUnits::Dimension;
using namespace SireStream;

static const RegisterMetaType<Evaluator> r_eval;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                       const Evaluator &eval)
{
    writeHeader(ds, r_eval, 1);

    SharedDataStream sds(ds);

    sds << eval.selected_atoms << static_cast<const MoleculeView&>(eval);

    return ds;
}

/** Deserialise from a binary datastream */
QDataStream &operator>>(QDataStream &ds, Evaluator &eval)
{
    VersionID v = readHeader(ds, r_eval);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> eval.selected_atoms >> static_cast<MoleculeView&>(eval);
    }
    else
        throw version_error(v, "1", r_eval, CODELOC);

    return ds;
}

/** Null constructor */
Evaluator::Evaluator() : ConcreteProperty<Evaluator,MoleculeView>()
{}

/** Construct from the passed molecule view */
Evaluator::Evaluator(const MoleculeView &molecule)
          : ConcreteProperty<Evaluator,MoleculeView>(molecule),
            selected_atoms(molecule.selection())
{}

/** Construct to evaluate for the entire molecule in 'moldata' */
Evaluator::Evaluator(const MoleculeData &moldata)
          : ConcreteProperty<Evaluator,MoleculeView>(moldata),
            selected_atoms(moldata)
{}

/** Construct to evaluate properties of the passed selected atoms
    of the molecule viewed in 'molecule' */
Evaluator::Evaluator(const MoleculeView &molecule,
                     const AtomSelection &atoms)
          : ConcreteProperty<Evaluator,MoleculeView>(molecule), selected_atoms(atoms)
{
    selected_atoms.assertCompatibleWith(this->data());
}

/** Construct to evaluate properties of the selected atoms of the
    passed molecule */
Evaluator::Evaluator(const MoleculeData &moldata,
                     const AtomSelection &atoms)
          : ConcreteProperty<Evaluator,MoleculeView>(moldata), selected_atoms(atoms)
{
    selected_atoms.assertCompatibleWith(this->data());
}

/** Copy constructor */
Evaluator::Evaluator(const Evaluator &other)
          : ConcreteProperty<Evaluator,MoleculeView>(other),
            selected_atoms(other.selected_atoms)
{}

/** Destructor */
Evaluator::~Evaluator()
{}

/** Copy assignment from another evaluator */
Evaluator& Evaluator::operator=(const Evaluator &other)
{
    if (this != &other)
    {
        MoleculeView::operator=(other);
        selected_atoms = other.selected_atoms;
    }

    return *this;
}

/** Copy assignment from another molecule */
Evaluator& Evaluator::operator=(const MoleculeView &other)
{
    MoleculeView::operator=(other);
    selected_atoms = other.selection();

    return *this;
}

/** Return a string representation of this evaluator */
QString Evaluator::toString() const
{
    return QObject::tr( "Evaluator( nAtoms() == %1 )" )
            .arg( selected_atoms.nSelected() );
}

/** Return whether or not this is empty */
bool Evaluator::isEmpty() const
{
    return selected_atoms.selectedNone();
}

/** Return whether or not this contains the whole molecule */
bool Evaluator::selectedAll() const
{
    return selected_atoms.selectedAll();
}

/** Return the selected atoms over which the properties
    will be evaluated */
AtomSelection Evaluator::selection() const
{
    return selected_atoms;
}

static void getMinMax(const CoordGroup &cgroup, Vector &min, Vector &max)
{
    //we can cheat by using the CoordGroup's aabox!
    min.setMin(cgroup.aaBox().minCoords());
    max.setMax(cgroup.aaBox().maxCoords());
}

static void getMinMax(const CoordGroup &cgroup, const QSet<Index> &idxs,
                      Vector &min, Vector &max)
{
    const Vector *cgroup_array = cgroup.constData();

    foreach (Index i, idxs)
    {
        const Vector &coords = cgroup_array[i];

        min.setMin(coords);
        max.setMax(coords);
    }
}

/** Return the axis-aligned box that just contains all of the
    atoms in this view

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
AABox Evaluator::aaBox(const PropertyMap &map) const
{
    if (selected_atoms.selectedNone())
        return AABox();

    //get the coordinates of the atoms
    const Property &prop = d->property(map["coordinates"]);
    const AtomCoords &coords = prop.asA<AtomCoords>();

    const CoordGroup *coords_array = coords.constData();
    int ncg = coords.count();

    //now get the minimum and maximum coordinates...
    Vector mincoords( std::numeric_limits<double>::max() );
    Vector maxcoords( -std::numeric_limits<double>::max() );

    if (selected_atoms.selectedAll())
    {
        for (int i=0; i<ncg; ++i)
        {
            getMinMax(coords_array[i], mincoords, maxcoords);
        }
    }
    else if (selected_atoms.selectedAllCutGroups())
    {
        for (CGIdx i(0); i<ncg; ++i)
        {
            if (selected_atoms.selectedAll(i))
            {
                getMinMax(coords_array[i], mincoords, maxcoords);
            }
            else
            {
                getMinMax(coords_array[i], selected_atoms.selectedAtoms(i),
                          mincoords, maxcoords);
            }
        }
    }
    else
    {
        foreach (CGIdx cgidx, selected_atoms.selectedCutGroups())
        {
            if (selected_atoms.selectedAll(cgidx))
            {
                getMinMax(coords_array[cgidx], mincoords, maxcoords);
            }
            else
            {
                getMinMax(coords_array[cgidx], selected_atoms.selectedAtoms(cgidx),
                          mincoords, maxcoords);
            }
        }
    }

    return AABox::from(mincoords, maxcoords);
}

/** Return the center of the selected atoms,
    using the passed property map to find the coordinates
    property of the molecule (the center is the point
    that is exactly in the middle of the atoms - i.e.
    halfway between the maximum and minimum coordinates

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
Vector Evaluator::center(const PropertyMap &map) const
{
    return this->aaBox(map).center();
}

/** Return the sphere that just encloses all of the atoms in this view

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
Sphere Evaluator::boundingSphere(const PropertyMap &map) const
{
    return this->aaBox(map).boundingSphere();
}

inline double getMass(const MolarMass &mass)
{
    return mass.value();
}

inline double getMass(const Element &element)
{
    return element.mass().value();
}

inline double getCharge(const Charge &charge)
{
    return charge.value();
}

/** Internal function used to calculate the center of mass of the selected atoms */
template<class T>
static Vector getCOM(const AtomCoords &coords, const AtomProperty<T> &masses,
                     const AtomSelection &selected_atoms)
{
    if (selected_atoms.selectedNone())
        return Vector(0);

    //calculate the center of mass
    Vector com(0);
    double mass(0);

    if (selected_atoms.selectedAll())
    {
        const Vector *coords_array = coords.array().constCoordsData();
        const T *masses_array = masses.array().constValueData();

        const int nats = coords.nAtoms();

        for (int i=0; i<nats; ++i)
        {
            const double m = ::getMass(masses_array[i]);
            com += m*coords_array[i];
            mass += m;
        }
    }
    else if (selected_atoms.selectedAllCutGroups())
    {
        for (CGIdx i(0); i<coords.nCutGroups(); ++i)
        {
            const Vector *coords_array = coords.constData(i);
            const T *masses_array = masses.constData(i);

            if (selected_atoms.selectedAll(i))
            {
                const int nats = coords.nAtoms(i);

                for (int j=0; j<nats; ++j)
                {
                    const double m = ::getMass(masses_array[j]);
                    com += m*coords_array[j];
                    mass += m;
                }
            }
            else
            {
                foreach (Index j, selected_atoms.selectedAtoms(i))
                {
                    const double m = ::getMass(masses_array[j]);
                    com += m*coords_array[j];
                    mass += m;
                }
            }
        }
    }
    else
    {
        foreach (CGIdx i, selected_atoms.selectedCutGroups())
        {
            const Vector *coords_array = coords.constData(i);
            const T *masses_array = masses.constData(i);

            if (selected_atoms.selectedAll(i))
            {
                const int nats = coords.nAtoms(i);

                for (int j=0; j<nats; ++j)
                {
                    const double m = ::getMass(masses_array[j]);
                    com += m*coords_array[j];
                    mass += m;
                }
            }
            else
            {
                foreach (Index j, selected_atoms.selectedAtoms(i))
                {
                    const double m = ::getMass(masses_array[j]);
                    com += m*coords_array[j];
                    mass += m;
                }
            }
        }
    }

    return com / mass;
}

static void addToInertia(const Vector &d, double m, Matrix &inertia)
{
    double *inertia_array = inertia.data();

    inertia_array[0] += m * (d.y()*d.y() + d.z()*d.z());
    inertia_array[4] += m * (d.x()*d.x() + d.z()*d.z());
    inertia_array[8] += m * (d.x()*d.x() + d.y()*d.y());

    inertia_array[1] -= m * d.x() * d.y();
    inertia_array[2] -= m * d.x() * d.z();
    inertia_array[5] -= m * d.y() * d.z();
}

// JM 08/14 buggy
//static void getPrincipalAxes(Matrix &inertia, Vector principal_moments)
//{
//    double *inertia_array = inertia.data();
//
//    //remove near-zero elements
//    for (int i=0; i<9; ++i)
//    {
//        if (inertia_array[i] < 1e-6 and inertia_array[i] > -1e-6)
//            inertia_array[i] = 0;
//    }
//
//    //symmetric matrix
//    //
//    //   0 1 2
//    //   3 4 5
//    //   6 7 8
//    //
//    inertia_array[3] = inertia_array[1];
//    inertia_array[6] = inertia_array[2];
//   inertia_array[7] = inertia_array[5];
//
//    std::pair<Vector,Matrix> eigs = inertia.diagonalise();
//
//    principal_moments = eigs.first;
//    inertia = eigs.second;
//
//    //if one or more of the eigenvalues is zero then we may have a problem
//    //because the wrong eigenvector direction may be chosen - in this case,
//    //we will build this eigenvector using a cross product to ensure that
//    //the right-hand-rule definition of our axes is maintained
//    //
//    // Also, even if we have three eigenvalues, we still need to make sure
//    // that a right-hand-rule set is chosen, rather than the left-hand set
//    bool zero_x = std::abs(principal_moments[0]) < 1e-6;
//    bool zero_y = std::abs(principal_moments[1]) < 1e-6;
//    bool zero_z = std::abs(principal_moments[2]) < 1e-6;
//
//    if (zero_x){ principal_moments.setX(0); }
//    if (zero_y){ principal_moments.setY(0); }
//    if (zero_z){ principal_moments.setZ(0); }
//
//    int n_zeroes = int(zero_x) + int(zero_y) + int(zero_z);
//
//    if (n_zeroes == 3)
//    {
//        //no axes!
//        inertia = Matrix(1);
//    }
//    else if (n_zeroes == 2)
//    {
//        //just one well-defined axis - I don't know how to handle this...
//        throw SireError::incomplete_code( QObject::tr(
//                "The code to get principal axes for molecules with only a single "
//                "eigenvalue has yet to be written... (%1 and %2)")
//                    .arg(principal_moments.toString(),
//                         inertia.toString()), CODELOC );
//    }
//    else if (n_zeroes == 1)
//    {
//        Vector r0 = inertia.row0();
//        Vector r1 = inertia.row1();
//        Vector r2 = inertia.row2();
//
//        if (zero_x)
//            r0 = Vector::cross(r1,r2);
//        else if (zero_y)
//            r1 = Vector::cross(r2,r0);
//        else if (zero_z)
//            r2 = Vector::cross(r0,r1);
//
//        inertia = Matrix(r0, r1, r2);
//    }
//    else
//    {
//        Vector r0 = inertia.row0();
//        Vector r1 = inertia.row1();
//
//        inertia = Matrix( r0, r1, Vector::cross(r0,r1) );
//    }
//}

/** Internal function used to get the principal axes of the selected atoms */
//template<class T>
//static AxisSet getPrincipalAxes(const AtomCoords &coords,
//                                const AtomProperty<T> &masses,
//                                const AtomSelection &selected_atoms,
//                                Vector &principal_moments)
//{
//    if (selected_atoms.selectedNone())
//        return AxisSet();
//
//    Vector com = ::getCOM(coords, masses, selected_atoms);
//
//   Matrix inertia(0);
//
//    if (selected_atoms.selectedAll())
//    {
//        const Vector *coords_array = coords.array().constCoordsData();
//        const T *masses_array = masses.array().constValueData();
//
//        const int nats = coords.nAtoms();
//
//        for (int i=0; i<nats; ++i)
//        {
//            ::addToInertia(coords_array[i]-com, ::getMass(masses_array[i]), inertia);
//        }
//    }
//    else if (selected_atoms.selectedAllCutGroups())
//    {
//        for (CGIdx i(0); i<coords.nCutGroups(); ++i)
//        {
//            const Vector *coords_array = coords.constData(i);
//            const T *masses_array = masses.constData(i);
//
//            if (selected_atoms.selectedAll(i))
//            {
//                const int nats = coords.nAtoms(i);
//
//                for (int j=0; j<nats; ++j)
//                {
//                    ::addToInertia(coords_array[j]-com, ::getMass(masses_array[j]),
//                                   inertia);
//                }
//            }
//            else
//            {
//                foreach (Index j, selected_atoms.selectedAtoms(i))
//                {
//                    ::addToInertia(coords_array[j]-com, ::getMass(masses_array[j]),
//                                   inertia);
//                }
//            }
//        }
//    }
//    else
//    {
//        foreach (CGIdx i, selected_atoms.selectedCutGroups())
//        {
//            const Vector *coords_array = coords.constData(i);
//           const T *masses_array = masses.constData(i);
//
//            if (selected_atoms.selectedAll(i))
//            {
//                const int nats = coords.nAtoms(i);
//
//                for (int j=0; j<nats; ++j)
//                {
//                    ::addToInertia(coords_array[j]-com, ::getMass(masses_array[j]),
//                                   inertia);
//                }
//            }
//            else
//            {
//                foreach (Index j, selected_atoms.selectedAtoms(i))
//                {
//                    ::addToInertia(coords_array[j]-com, ::getMass(masses_array[j]),
//                                   inertia);
//                }
//            }
//        }
//   }
//
//    ::getPrincipalAxes(inertia, principal_moments);
//
//    return AxisSet(inertia, com);
//}

/** Internal function used to calculate the total mass of the selected atoms */
template<class T>
static MolarMass getMass(const AtomProperty<T> &masses,
                         const AtomSelection &selected_atoms)
{
    if (selected_atoms.selectedNone())
        return MolarMass(0);

    double mass(0);

    if (selected_atoms.selectedAll())
    {
        const T *masses_array = masses.array().constValueData();

        const int nats = masses.nAtoms();

        for (int i=0; i<nats; ++i)
        {
            mass += ::getMass(masses_array[i]);
        }
    }
    else if (selected_atoms.selectedAllCutGroups())
    {
        for (CGIdx i(0); i<masses.nCutGroups(); ++i)
        {
            const T *masses_array = masses.constData(i);

            if (selected_atoms.selectedAll(i))
            {
                const int nats = masses.nAtoms(i);

                for (int j=0; j<nats; ++j)
                {
                    mass += ::getMass(masses_array[j]);
                }
            }
            else
            {
                foreach (Index j, selected_atoms.selectedAtoms(i))
                {
                    mass += ::getMass(masses_array[j]);
                }
            }
        }
    }
    else
    {
        foreach (CGIdx i, selected_atoms.selectedCutGroups())
        {
            const T *masses_array = masses.constData(i);

            if (selected_atoms.selectedAll(i))
            {
                const int nats = masses.nAtoms(i);

                for (int j=0; j<nats; ++j)
                {
                    mass += ::getMass(masses_array[j]);
                }
            }
            else
            {
                foreach (Index j, selected_atoms.selectedAtoms(i))
                {
                    mass += ::getMass(masses_array[j]);
                }
            }
        }
    }

    return MolarMass(mass);
}

/** Internal function used to calculate the total charge of the selected atoms */
template<class T>
static Charge getCharge(const AtomProperty<T> &charges,
                        const AtomSelection &selected_atoms)
{
    if (selected_atoms.selectedNone())
        return Charge(0);

    double charge(0);

    if (selected_atoms.selectedAll())
    {
        const T *charges_array = charges.array().constValueData();

        const int nats = charges.nAtoms();

        for (int i=0; i<nats; ++i)
        {
            charge += ::getCharge(charges_array[i]);
        }
    }
    else if (selected_atoms.selectedAllCutGroups())
    {
        for (CGIdx i(0); i<charges.nCutGroups(); ++i)
        {
            const T *charges_array = charges.constData(i);

            if (selected_atoms.selectedAll(i))
            {
                const int nats = charges.nAtoms(i);

                for (int j=0; j<nats; ++j)
                {
                    charge += ::getCharge(charges_array[j]);
                }
            }
            else
            {
                foreach (Index j, selected_atoms.selectedAtoms(i))
                {
                    charge += ::getCharge(charges_array[j]);
                }
            }
        }
    }
    else
    {
        foreach (CGIdx i, selected_atoms.selectedCutGroups())
        {
            const T *charges_array = charges.constData(i);

            if (selected_atoms.selectedAll(i))
            {
                const int nats = charges.nAtoms(i);

                for (int j=0; j<nats; ++j)
                {
                    charge += ::getCharge(charges_array[j]);
                }
            }
            else
            {
                foreach (Index j, selected_atoms.selectedAtoms(i))
                {
                    charge += ::getCharge(charges_array[j]);
                }
            }
        }
    }

    return Charge(charge);
}

/** Return the mass of the selected part of this molecule, using
    the supplied map to find either the mass property, or if that
    does not exist, using the element property

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
MolarMass Evaluator::mass(const PropertyMap &map) const
{
    const PropertyName mass_property = map["mass"];

    if (d->hasProperty(mass_property))
    {
        const AtomMasses &masses = d->property(mass_property).asA<AtomMasses>();

        return ::getMass(masses, selected_atoms);
    }
    else
    {
        const AtomElements &elements = d->property(map["element"]).asA<AtomElements>();

        return ::getMass(elements, selected_atoms);
    }
}

/** Return the total charge of the selected part of the molecule, using
    the supplied map to find the "charge" property

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
Charge Evaluator::charge(const PropertyMap &map) const
{
    const AtomCharges &charges = d->property(map["charge"]).asA<AtomCharges>();

    return ::getCharge(charges, selected_atoms);
}

/** Return the centroid of these atoms - this is the average
    of the coordinates

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
Vector Evaluator::centroid(const PropertyMap &map) const
{
    const AtomCoords &coords = d->property(map["coordinates"]).asA<AtomCoords>();

    if (selected_atoms.selectedNone())
        return Vector(0);

    if (selected_atoms.selectedNone())
        return Vector(0);

    Vector cent(0);
    int natoms(0);

    if (selected_atoms.selectedAll())
    {
        const Vector *coords_array = coords.array().constCoordsData();

        for (int i=0; i<coords.nAtoms(); ++i)
        {
            cent += coords_array[i];
            ++natoms;
        }
    }
    else if (selected_atoms.selectedAllCutGroups())
    {
        for (CGIdx i(0); i<coords.nCutGroups(); ++i)
        {
            const Vector *coords_array = coords.constData(i);

            if (selected_atoms.selectedAll(i))
            {
                for (int j=0; j<coords.nAtoms(i); ++j)
                {
                    cent += coords_array[j];
                    ++natoms;
                }
            }
            else
            {
                foreach (Index j, selected_atoms.selectedAtoms(i))
                {
                    cent += coords_array[j];
                    ++natoms;
                }
            }
        }
    }
    else
    {
        foreach (CGIdx i, selected_atoms.selectedCutGroups())
        {
            const Vector *coords_array = coords.constData(i);

            if (selected_atoms.selectedAll(i))
            {
                for (int j=0; j<coords.nAtoms(i); ++j)
                {
                    cent += coords_array[j];
                    ++natoms;
                }
            }
            else
            {
                foreach (Index j, selected_atoms.selectedAtoms(i))
                {
                    cent += coords_array[j];
                    ++natoms;
                }
            }
        }
    }

    return cent / natoms;
}

/** Return the center of geometry of this part of the molecule.
    This is the mid-point between the maximum coordinates and
    minimum coordinates

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
Vector Evaluator::centerOfGeometry(const PropertyMap &map) const
{
    const AtomCoords &coords = d->property(map["coordinates"]).asA<AtomCoords>();

    if (selected_atoms.selectedNone())
        return Vector(0);

    if (selected_atoms.selectedNone())
        return Vector(0);

    Vector mincoords( std::numeric_limits<double>::max() );
    Vector maxcoords( -std::numeric_limits<double>::max() );

    if (selected_atoms.selectedAll())
    {
        const Vector *coords_array = coords.array().constCoordsData();

        for (int i=0; i<coords.nAtoms(); ++i)
        {
            mincoords.setMin(coords_array[i]);
            maxcoords.setMax(coords_array[i]);
        }
    }
    else if (selected_atoms.selectedAllCutGroups())
    {
        for (CGIdx i(0); i<coords.nCutGroups(); ++i)
        {
            const Vector *coords_array = coords.constData(i);

            if (selected_atoms.selectedAll(i))
            {
                for (int j=0; j<coords.nAtoms(i); ++j)
                {
                    mincoords.setMin(coords_array[j]);
                    maxcoords.setMax(coords_array[j]);
                }
            }
            else
            {
                foreach (Index j, selected_atoms.selectedAtoms(i))
                {
                    mincoords.setMin(coords_array[j]);
                    maxcoords.setMax(coords_array[j]);
                }
            }
        }
    }
    else
    {
        foreach (CGIdx i, selected_atoms.selectedCutGroups())
        {
            const Vector *coords_array = coords.constData(i);

            if (selected_atoms.selectedAll(i))
            {
                for (int j=0; j<coords.nAtoms(i); ++j)
                {
                    mincoords.setMin(coords_array[j]);
                    maxcoords.setMax(coords_array[j]);
                }
            }
            else
            {
                foreach (Index j, selected_atoms.selectedAtoms(i))
                {
                    mincoords.setMin(coords_array[j]);
                    maxcoords.setMax(coords_array[j]);
                }
            }
        }
    }

    return mincoords + 0.5*(maxcoords-mincoords);
}

/** Return the center of mass of this part of the molecule

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
Vector Evaluator::centerOfMass(const PropertyMap &map) const
{
    const AtomCoords &coords = d->property(map["coordinates"]).asA<AtomCoords>();

    const PropertyName mass_property = map["mass"];

    if (d->hasProperty(mass_property))
    {
        const AtomMasses &masses = d->property(mass_property).asA<AtomMasses>();

        return ::getCOM(coords, masses, selected_atoms);
    }
    else
    {
        const AtomElements &elements = d->property(map["element"]).asA<AtomElements>();

        return ::getCOM(coords, elements, selected_atoms);
    }
}

/** Return the principal axes of this view - this uses
    the "coordinates", and "mass" or "element" properties
    to find the moment of inertia tensor for this view, and
    then diagonalises that to obtain the principal axes. These
    axes are constructed to follow the right-hand-rule.
    This returns the principal moments of inertia in
    'principal_moments' */
AxisSet Evaluator::principalAxes(Vector &principal_moments,
                                 const PropertyMap &map) const
{
    const PropertyName coords_property = map["coordinates"];
    const PropertyName mass_property = map["mass"];

    const AtomCoords &coords = d->property(map["coordinates"])
                                    .asA<AtomCoords>();

    if (d->hasProperty(mass_property))
    {
        const AtomMasses &masses = d->property(mass_property).asA<AtomMasses>();

        // JM 08/14 buggy in current code
#ifdef _MSC_VER
        #pragma WARNING(BUGGY IN CURRENT CODE)
#else
        #warning BUGGY IN CURRENT CODE
#endif
        throw SireError::program_bug( QObject::tr("CODE IS BROKEN"), CODELOC );
        //return ::getPrincipalAxes(coords, masses, selected_atoms,
        //                          principal_moments);
    }
    else
    {
        const AtomElements &elements = d->property(map["element"]).asA<AtomElements>();

        // JM 08/14 buggy in current code
#ifdef _MSC_VER
        #pragma WARNING(BUGGY IN CURRENT CODE)
#else
        #warning BUGGY IN CURRENT CODE
#endif
        throw SireError::program_bug( QObject::tr("CODE IS BROKEN"), CODELOC );
        //return ::getPrincipalAxes(coords, elements, selected_atoms,
        //                          principal_moments);
    }

}

/** Return the principal axes of this view - this uses
    the "coordinates", and "mass" or "element" properties
    to find the moment of inertia tensor for this view, and
    then diagonalises that to obtain the principal axes. These
    axes are constructed to follow the right-hand-rule.
*/
AxisSet Evaluator::principalAxes(const PropertyMap &map) const
{
    Vector principal_moments;
    return this->principalAxes(principal_moments,map);
}

AxisSet Evaluator::alignmentAxes(const MoleculeView &other,
                                 const AtomMatcher &matcher,
                                 const PropertyMap &map) const
{
    throw SireError::incomplete_code(CODELOC);
    return AxisSet();
}

AxisSet Evaluator::alignmentAxes(const MoleculeView &other,
                                 const AtomMatcher &matcher,
                                 const PropertyMap &map0,
                                 const PropertyMap &map1) const
{
    throw SireError::incomplete_code(CODELOC);
    return AxisSet();
}

static CGAtomIdx selectOnly(const AtomID &atom, const MoleculeInfoData &molinfo,
                            const AtomSelection &selected_atoms)
{
    QList<AtomIdx> atomidxs = molinfo.map(atom);

    QMutableListIterator<AtomIdx> it(atomidxs);

    while (it.hasNext())
    {
        if (not selected_atoms.selected(it.next()))
            it.remove();
    }

    if (atomidxs.isEmpty())
        throw SireMol::missing_atom( QObject::tr(
                "None of the selected atoms match the ID \"%1\".")
                    .arg(atom.toString()), CODELOC );

    else if (atomidxs.count() > 1)
        throw SireMol::duplicate_atom( QObject::tr(
                "More than one selected atom matches the ID \"%1\".")
                    .arg(atom.toString()), CODELOC );

    return molinfo.cgAtomIdx(atomidxs.at(0));
}

/** Measure the distance between the atoms 'atom0' and 'atom1'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
Length Evaluator::measure(const AtomID &atom0, const AtomID &atom1,
                          const PropertyMap &map) const
{
    const AtomCoords &coords = data().property(map["coordinates"])
                                     .asA<AtomCoords>();

    return Length( Line( coords[::selectOnly(atom0, data().info(), selected_atoms)],
                         coords[::selectOnly(atom1, data().info(), selected_atoms)] )
                            .length()
                 );
}

/** Measure the angle between the atoms 'atom0', 'atom1' and 'atom2'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
Angle Evaluator::measure(const AtomID &atom0, const AtomID &atom1,
                         const AtomID &atom2, const PropertyMap &map) const
{
    const AtomCoords &coords = data().property(map["coordinates"])
                                     .asA<AtomCoords>();

    return Triangle( coords[::selectOnly(atom0, data().info(), selected_atoms)],
                     coords[::selectOnly(atom1, data().info(), selected_atoms)],
                     coords[::selectOnly(atom2, data().info(), selected_atoms)]
                   ).angle();
}

/** Measure the dihedral between the atoms 'atom0', 'atom1', 'atom2' and 'atom3'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
Angle Evaluator::measure(const AtomID &atom0, const AtomID &atom1,
                         const AtomID &atom2, const AtomID &atom3,
                         const PropertyMap &map) const
{
    const AtomCoords &coords = data().property(map["coordinates"])
                                     .asA<AtomCoords>();

    return Torsion( coords[::selectOnly(atom0, data().info(), selected_atoms)],
                    coords[::selectOnly(atom1, data().info(), selected_atoms)],
                    coords[::selectOnly(atom2, data().info(), selected_atoms)],
                    coords[::selectOnly(atom3, data().info(), selected_atoms)]
                  ).angle();
}

/** Measure the length of the bond 'bond'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
Length Evaluator::measure(const BondID &bond, const PropertyMap &map) const
{
    return measure(bond.atom0(), bond.atom1(), map);
}

/** Measure the size of the angle 'angle'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
Angle Evaluator::measure(const AngleID &angle, const PropertyMap &map) const
{
    return measure(angle.atom0(), angle.atom1(), angle.atom2(), map);
}

/** Measure the size of the dihedral 'dihedral'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
Angle Evaluator::measure(const DihedralID &dihedral, const PropertyMap &map) const
{
    return measure(dihedral.atom0(), dihedral.atom1(),
                   dihedral.atom2(), dihedral.atom3(), map);
}

/** Find the maximum common substructure of this molecule view with 'other'. This
    returns the mapping from this structure to 'other' for the matching parts,
    using the optionally supplied propertymap to find the elements, masses,
    connectivity and coordinates of the two molecules, with the passed 'atommatcher'
    used to pre-match atoms before the common substructure search (useful to speed
    up the search and to enforce matching sub-parts) */
QHash<AtomIdx,AtomIdx> Evaluator::findMCS(const MoleculeView &other,
                                          const AtomMatcher &matcher,
                                          const PropertyMap &map0,
                                          const PropertyMap &map1,
                                          bool verbose) const
{
    return this->findMCS(other, matcher, 5*second, map0, map1, verbose);
}

/** Find the maximum common substructure of this molecule view with 'other'. This
    returns the mapping from this structure to 'other' for the matching parts,
    using the optionally supplied propertymap to find the elements, masses,
    connectivity and coordinates of the two molecules */
QHash<AtomIdx,AtomIdx> Evaluator::findMCS(const MoleculeView &other,
                                          const PropertyMap &map,
                                          bool verbose) const
{
    return this->findMCS(other, AtomMultiMatcher(), map, map, verbose);
}

/** Find the maximum common substructure of this molecule view with 'other'. This
    returns the mapping from this structure to 'other' for the matching parts,
    using map0 and map1 to find the elements, masses,
    connectivity and coordinates of the two molecules respectively */
QHash<AtomIdx,AtomIdx> Evaluator::findMCS(const MoleculeView &other,
                                          const PropertyMap &map0,
                                          const PropertyMap &map1,
                                          bool verbose) const
{
    return this->findMCS(other, AtomMultiMatcher(), map0, map1, verbose);
}

/** Find the maximum common substructure of this molecule view with 'other'. This
    returns the mapping from this structure to 'other' for the matching parts,
    using the optionally supplied propertymap to find the elements, masses,
    connectivity and coordinates of the two molecules, with the passed 'atommatcher'
    used to pre-match atoms before the common substructure search (useful to speed
    up the search and to enforce matching sub-parts) */
QHash<AtomIdx,AtomIdx> Evaluator::findMCS(const MoleculeView &other,
                                          const AtomMatcher &atommatcher,
                                          const PropertyMap &map,
                                          bool verbose) const
{
    return this->findMCS(other, atommatcher, map, map, verbose);
}

/** Find the maximum common substructure of this molecule view with 'other'. This
    returns the mapping from this structure to 'other' for the matching parts,
    using the optionally supplied propertymap to find the elements, masses,
    connectivity and coordinates of the two molecules. Terminate the calculation
    returning the best match found within 'timeout'. */
QHash<AtomIdx,AtomIdx> Evaluator::findMCS(const MoleculeView &other,
                                          const Time &timeout,
                                          const PropertyMap &map,
                                          bool verbose) const
{
    return this->findMCS(other, AtomMultiMatcher(), timeout, map, map, verbose);
}

/** Find the maximum common substructure of this molecule view with 'other'. This
    returns the mapping from this structure to 'other' for the matching parts,
    using map0 and map1 to find the elements, masses,
    connectivity and coordinates of the two molecules respectively. Terminate the calculation
    returning the best match found within 'timeout'. */
QHash<AtomIdx,AtomIdx> Evaluator::findMCS(const MoleculeView &other,
                                          const Time &timeout,
                                          const PropertyMap &map0,
                                          const PropertyMap &map1,
                                          bool verbose) const
{
    return this->findMCS(other, AtomMultiMatcher(), timeout, map0, map1, verbose);
}

/** Find the maximum common substructure of this molecule view with 'other'. This
    returns the mapping from this structure to 'other' for the matching parts,
    using the optionally supplied propertymap to find the elements, masses,
    connectivity and coordinates of the two molecules, with the passed 'atommatcher'
    used to pre-match atoms before the common substructure search (useful to speed
    up the search and to enforce matching sub-parts). Terminate the calculation
    returning the best match found within 'timeout'. */
QHash<AtomIdx,AtomIdx> Evaluator::findMCS(const MoleculeView &other,
                                          const AtomMatcher &atommatcher,
                                          const Time &timeout,
                                          const PropertyMap &map,
                                          bool verbose) const
{
    return this->findMCS(other, atommatcher, timeout, map, map, verbose);
}

/** Find the maximum common substructure of this molecule view with 'other'. This
    returns the mapping from this structure to 'other' for the matching parts,
    using the optionally supplied propertymap to find the elements, masses,
    connectivity and coordinates of the two molecules, with the passed 'atommatcher'
    used to pre-match atoms before the common substructure search (useful to speed
    up the search and to enforce matching sub-parts).

    If 'match_light_atoms' is true, then include light atoms (e.g. hydrogen)
    in the match. This may make things slower...
*/
QHash<AtomIdx,AtomIdx> Evaluator::findMCS(const MoleculeView &other,
                                          const AtomMatcher &matcher,
                                          bool match_light_atoms,
                                          const PropertyMap &map0,
                                          const PropertyMap &map1,
                                          int min_heavy_protons,
                                          bool verbose) const
{
    return this->findMCS(other, matcher, 5*second, match_light_atoms, map0, map1, min_heavy_protons, verbose);
}

/** Find the maximum common substructure of this molecule view with 'other'. This
    returns the mapping from this structure to 'other' for the matching parts,
    using the optionally supplied propertymap to find the elements, masses,
    connectivity and coordinates of the two molecules

    If 'match_light_atoms' is true, then include light atoms (e.g. hydrogen)
    in the match. This may make things slower...
*/
QHash<AtomIdx,AtomIdx> Evaluator::findMCS(const MoleculeView &other,
                                          bool match_light_atoms,
                                          const PropertyMap &map,
                                          int min_heavy_protons,
                                          bool verbose) const
{
    return this->findMCS(other, AtomMultiMatcher(), match_light_atoms, map, map, min_heavy_protons, verbose);
}

/** Find the maximum common substructure of this molecule view with 'other'. This
    returns the mapping from this structure to 'other' for the matching parts,
    using map0 and map1 to find the elements, masses,
    connectivity and coordinates of the two molecules respectively

    If 'match_light_atoms' is true, then include light atoms (e.g. hydrogen)
    in the match. This may make things slower...

*/
QHash<AtomIdx,AtomIdx> Evaluator::findMCS(const MoleculeView &other,
                                          bool match_light_atoms,
                                          const PropertyMap &map0,
                                          const PropertyMap &map1,
                                          int min_heavy_protons,
                                          bool verbose) const
{
    return this->findMCS(other, AtomMultiMatcher(), match_light_atoms, map0, map1, min_heavy_protons, verbose);
}

/** Find the maximum common substructure of this molecule view with 'other'. This
    returns the mapping from this structure to 'other' for the matching parts,
    using the optionally supplied propertymap to find the elements, masses,
    connectivity and coordinates of the two molecules, with the passed 'atommatcher'
    used to pre-match atoms before the common substructure search (useful to speed
    up the search and to enforce matching sub-parts)

    If 'match_light_atoms' is true, then include light atoms (e.g. hydrogen)
    in the match. This may make things slower...
*/
QHash<AtomIdx,AtomIdx> Evaluator::findMCS(const MoleculeView &other,
                                          const AtomMatcher &atommatcher,
                                          bool match_light_atoms,
                                          const PropertyMap &map,
                                          int min_heavy_protons,
                                          bool verbose) const
{
    return this->findMCS(other, atommatcher, match_light_atoms, map, map, min_heavy_protons, verbose);
}

/** Find the maximum common substructure of this molecule view with 'other'. This
    returns the mapping from this structure to 'other' for the matching parts,
    using the optionally supplied propertymap to find the elements, masses,
    connectivity and coordinates of the two molecules. Terminate the calculation
    returning the best match found within 'timeout'.

    If 'match_light_atoms' is true, then include light atoms (e.g. hydrogen)
    in the match. This may make things slower...
*/
QHash<AtomIdx,AtomIdx> Evaluator::findMCS(const MoleculeView &other,
                                          const Time &timeout,
                                          bool match_light_atoms,
                                          const PropertyMap &map,
                                          int min_heavy_protons,
                                          bool verbose) const
{
    return this->findMCS(other, AtomMultiMatcher(), timeout, match_light_atoms, map, map, min_heavy_protons, verbose);
}

/** Find the maximum common substructure of this molecule view with 'other'. This
    returns the mapping from this structure to 'other' for the matching parts,
    using map0 and map1 to find the elements, masses,
    connectivity and coordinates of the two molecules respectively. Terminate the calculation
    returning the best match found within 'timeout'.

    If 'match_light_atoms' is true, then include light atoms (e.g. hydrogen)
    in the match. This may make things slower...
*/
QHash<AtomIdx,AtomIdx> Evaluator::findMCS(const MoleculeView &other,
                                          const Time &timeout,
                                          bool match_light_atoms,
                                          const PropertyMap &map0,
                                          const PropertyMap &map1,
                                          int min_heavy_protons,
                                          bool verbose) const
{
    return this->findMCS(other, AtomMultiMatcher(), timeout, match_light_atoms, map0, map1, min_heavy_protons, verbose);
}

/** Find the maximum common substructure of this molecule view with 'other'. This
    returns the mapping from this structure to 'other' for the matching parts,
    using the optionally supplied propertymap to find the elements, masses,
    connectivity and coordinates of the two molecules, with the passed 'atommatcher'
    used to pre-match atoms before the common substructure search (useful to speed
    up the search and to enforce matching sub-parts). Terminate the calculation
    returning the best match found within 'timeout'.

    If 'match_light_atoms' is true, then include light atoms (e.g. hydrogen)
    in the match. This may make things slower...
*/
QHash<AtomIdx,AtomIdx> Evaluator::findMCS(const MoleculeView &other,
                                          const AtomMatcher &atommatcher,
                                          const Time &timeout,
                                          bool match_light_atoms,
                                          const PropertyMap &map,
                                          int min_heavy_protons,
                                          bool verbose) const
{
    return this->findMCS(other, atommatcher, timeout, match_light_atoms, map, map, verbose);
}

/** Find the maximum common substructure of this molecule view with 'other'. This
    returns the mapping from this structure to 'other' for the matching parts,
    using the optionally supplied propertymap to find the elements, masses,
    connectivity and coordinates of the two molecules, with the passed 'atommatcher'
    used to pre-match atoms before the common substructure search (useful to speed
    up the search and to enforce matching sub-parts) */
QHash<AtomIdx,AtomIdx> Evaluator::findMCS(const MoleculeView &other,
                                          const AtomMatcher &matcher,
                                          const Time &timeout,
                                          const PropertyMap &map0,
                                          const PropertyMap &map1,
                                          bool verbose) const
{
    return this->findMCS(other, matcher, timeout, false, map0, map1, verbose);
}

/** Find the maximum common substructure of this molecule view with 'other'. This
    returns all mappings from this structure to 'other' for the matching parts,
    using the optionally supplied propertymap to find the elements, masses,
    connectivity and coordinates of the two molecules, with the passed 'atommatcher'
    used to pre-match atoms before the common substructure search (useful to speed
    up the search and to enforce matching sub-parts) */
QVector<QHash<AtomIdx,AtomIdx> > Evaluator::findMCSmatches(const MoleculeView &other,
                                                           const AtomMatcher &matcher,
                                                           const PropertyMap &map0,
                                                           const PropertyMap &map1,
                                                           bool verbose) const
{
    return this->findMCSmatches(other, matcher, 5*second, map0, map1, verbose);
}

/** Find the maximum common substructure of this molecule view with 'other'. This
    returns all mappings from this structure to 'other' for the matching parts,
    using the optionally supplied propertymap to find the elements, masses,
    connectivity and coordinates of the two molecules */
QVector<QHash<AtomIdx,AtomIdx> > Evaluator::findMCSmatches(const MoleculeView &other,
                                                           const PropertyMap &map,
                                                           bool verbose) const
{
    return this->findMCSmatches(other, AtomMultiMatcher(), map, map, verbose);
}

/** Find the maximum common substructure of this molecule view with 'other'. This
    returns all mappings from this structure to 'other' for the matching parts,
    using map0 and map1 to find the elements, masses,
    connectivity and coordinates of the two molecules respectively */
QVector<QHash<AtomIdx,AtomIdx> > Evaluator::findMCSmatches(const MoleculeView &other,
                                                           const PropertyMap &map0,
                                                           const PropertyMap &map1,
                                                           bool verbose) const
{
    return this->findMCSmatches(other, AtomMultiMatcher(), map0, map1, verbose);
}

/** Find the maximum common substructure of this molecule view with 'other'. This
    returns all mappings from this structure to 'other' for the matching parts,
    using the optionally supplied propertymap to find the elements, masses,
    connectivity and coordinates of the two molecules, with the passed 'atommatcher'
    used to pre-match atoms before the common substructure search (useful to speed
    up the search and to enforce matching sub-parts) */
QVector<QHash<AtomIdx,AtomIdx> > Evaluator::findMCSmatches(const MoleculeView &other,
                                                           const AtomMatcher &atommatcher,
                                                           const PropertyMap &map,
                                                           bool verbose) const
{
    return this->findMCSmatches(other, atommatcher, map, map, verbose);
}

/** Find the maximum common substructure of this molecule view with 'other'. This
    returns all mappings from this structure to 'other' for the matching parts,
    using the optionally supplied propertymap to find the elements, masses,
    connectivity and coordinates of the two molecules. Terminate the calculation
    returning the best match found within 'timeout'. */
QVector<QHash<AtomIdx,AtomIdx> > Evaluator::findMCSmatches(const MoleculeView &other,
                                                           const Time &timeout,
                                                           const PropertyMap &map,
                                                           bool verbose) const
{
    return this->findMCSmatches(other, AtomMultiMatcher(), timeout, map, map, verbose);
}

/** Find the maximum common substructure of this molecule view with 'other'. This
    returns all mappings from this structure to 'other' for the matching parts,
    using map0 and map1 to find the elements, masses,
    connectivity and coordinates of the two molecules respectively. Terminate the calculation
    returning the best match found within 'timeout'. */
QVector<QHash<AtomIdx,AtomIdx> > Evaluator::findMCSmatches(const MoleculeView &other,
                                                           const Time &timeout,
                                                           const PropertyMap &map0,
                                                           const PropertyMap &map1,
                                                           bool verbose) const
{
    return this->findMCSmatches(other, AtomMultiMatcher(), timeout, map0, map1, verbose);
}

/** Find the maximum common substructure of this molecule view with 'other'. This
    returns all mappings from this structure to 'other' for the matching parts,
    using the optionally supplied propertymap to find the elements, masses,
    connectivity and coordinates of the two molecules, with the passed 'atommatcher'
    used to pre-match atoms before the common substructure search (useful to speed
    up the search and to enforce matching sub-parts). Terminate the calculation
    returning the best match found within 'timeout'. */
QVector<QHash<AtomIdx,AtomIdx> > Evaluator::findMCSmatches(const MoleculeView &other,
                                                           const AtomMatcher &atommatcher,
                                                           const Time &timeout,
                                                           const PropertyMap &map,
                                                           bool verbose) const
{
    return this->findMCSmatches(other, atommatcher, timeout, map, map, verbose);
}

/** Find the maximum common substructure of this molecule view with 'other'. This
    returns all mappings from this structure to 'other' for the matching parts,
    using the optionally supplied propertymap to find the elements, masses,
    connectivity and coordinates of the two molecules, with the passed 'atommatcher'
    used to pre-match atoms before the common substructure search (useful to speed
    up the search and to enforce matching sub-parts).

    If 'match_light_atoms' is true, then include light atoms (e.g. hydrogen)
    in the match. This may make things slower...
*/
QVector<QHash<AtomIdx,AtomIdx> > Evaluator::findMCSmatches(const MoleculeView &other,
                                                           const AtomMatcher &matcher,
                                                           bool match_light_atoms,
                                                           const PropertyMap &map0,
                                                           const PropertyMap &map1,
                                                           int min_heavy_protons,
                                                           bool verbose) const
{
    return this->findMCSmatches(other, matcher, 5*second, match_light_atoms, map0, map1, min_heavy_protons, verbose);
}

/** Find the maximum common substructure of this molecule view with 'other'. This
    returns all mappings from this structure to 'other' for the matching parts,
    using the optionally supplied propertymap to find the elements, masses,
    connectivity and coordinates of the two molecules

    If 'match_light_atoms' is true, then include light atoms (e.g. hydrogen)
    in the match. This may make things slower...
*/
QVector<QHash<AtomIdx,AtomIdx> > Evaluator::findMCSmatches(const MoleculeView &other,
                                                           bool match_light_atoms,
                                                           const PropertyMap &map,
                                                           int min_heavy_protons,
                                                           bool verbose) const
{
    return this->findMCSmatches(other, AtomMultiMatcher(), match_light_atoms, map, map, min_heavy_protons, verbose);
}

/** Find the maximum common substructure of this molecule view with 'other'. This
    returns all mappings from this structure to 'other' for the matching parts,
    using map0 and map1 to find the elements, masses,
    connectivity and coordinates of the two molecules respectively

    If 'match_light_atoms' is true, then include light atoms (e.g. hydrogen)
    in the match. This may make things slower...

*/
QVector<QHash<AtomIdx,AtomIdx> > Evaluator::findMCSmatches(const MoleculeView &other,
                                                           bool match_light_atoms,
                                                           const PropertyMap &map0,
                                                           const PropertyMap &map1,
                                                           int min_heavy_protons,
                                                           bool verbose) const
{
    return this->findMCSmatches(other, AtomMultiMatcher(), match_light_atoms, map0, map1, min_heavy_protons, verbose);
}

/** Find the maximum common substructure of this molecule view with 'other'. This
    returns all mappings from this structure to 'other' for the matching parts,
    using the optionally supplied propertymap to find the elements, masses,
    connectivity and coordinates of the two molecules, with the passed 'atommatcher'
    used to pre-match atoms before the common substructure search (useful to speed
    up the search and to enforce matching sub-parts)

    If 'match_light_atoms' is true, then include light atoms (e.g. hydrogen)
    in the match. This may make things slower...
*/
QVector<QHash<AtomIdx,AtomIdx> > Evaluator::findMCSmatches(const MoleculeView &other,
                                                           const AtomMatcher &atommatcher,
                                                           bool match_light_atoms,
                                                           const PropertyMap &map,
                                                           int min_heavy_protons,
                                                           bool verbose) const
{
    return this->findMCSmatches(other, atommatcher, match_light_atoms, map, map, min_heavy_protons, verbose);
}

/** Find the maximum common substructure of this molecule view with 'other'. This
    returns all mappings from this structure to 'other' for the matching parts,
    using the optionally supplied propertymap to find the elements, masses,
    connectivity and coordinates of the two molecules. Terminate the calculation
    returning the best match found within 'timeout'.

    If 'match_light_atoms' is true, then include light atoms (e.g. hydrogen)
    in the match. This may make things slower...
*/
QVector<QHash<AtomIdx,AtomIdx> > Evaluator::findMCSmatches(const MoleculeView &other,
                                                           const Time &timeout,
                                                           bool match_light_atoms,
                                                           const PropertyMap &map,
                                                           int min_heavy_protons,
                                                           bool verbose) const
{
    return this->findMCSmatches(other, AtomMultiMatcher(), timeout, match_light_atoms, map, map, min_heavy_protons, verbose);
}

/** Find the maximum common substructure of this molecule view with 'other'. This
    returns all mappings from this structure to 'other' for the matching parts,
    using map0 and map1 to find the elements, masses,
    connectivity and coordinates of the two molecules respectively. Terminate the calculation
    returning the best match found within 'timeout'.

    If 'match_light_atoms' is true, then include light atoms (e.g. hydrogen)
    in the match. This may make things slower...
*/
QVector<QHash<AtomIdx,AtomIdx> > Evaluator::findMCSmatches(const MoleculeView &other,
                                                           const Time &timeout,
                                                           bool match_light_atoms,
                                                           const PropertyMap &map0,
                                                           const PropertyMap &map1,
                                                           int min_heavy_protons,
                                                           bool verbose) const
{
    return this->findMCSmatches(other, AtomMultiMatcher(), timeout, match_light_atoms, map0, map1, min_heavy_protons, verbose);
}

/** Find the maximum common substructure of this molecule view with 'other'. This
    returns all mappings from this structure to 'other' for the matching parts,
    using the optionally supplied propertymap to find the elements, masses,
    connectivity and coordinates of the two molecules, with the passed 'atommatcher'
    used to pre-match atoms before the common substructure search (useful to speed
    up the search and to enforce matching sub-parts). Terminate the calculation
    returning the best match found within 'timeout'.

    If 'match_light_atoms' is true, then include light atoms (e.g. hydrogen)
    in the match. This may make things slower...
*/
QVector<QHash<AtomIdx,AtomIdx> > Evaluator::findMCSmatches(const MoleculeView &other,
                                                           const AtomMatcher &atommatcher,
                                                           const Time &timeout,
                                                           bool match_light_atoms,
                                                           const PropertyMap &map,
                                                           int min_heavy_protons,
                                                           bool verbose) const
{
    return this->findMCSmatches(other, atommatcher, timeout, match_light_atoms, map, map, min_heavy_protons, verbose);
}

/** Find the maximum common substructure of this molecule view with 'other'. This
    returns all mappings from this structure to 'other' for the matching parts,
    using the optionally supplied propertymap to find the elements, masses,
    connectivity and coordinates of the two molecules, with the passed 'atommatcher'
    used to pre-match atoms before the common substructure search (useful to speed
    up the search and to enforce matching sub-parts) */
QVector<QHash<AtomIdx,AtomIdx> > Evaluator::findMCSmatches(const MoleculeView &other,
                                                           const AtomMatcher &matcher,
                                                           const Time &timeout,
                                                           const PropertyMap &map0,
                                                           const PropertyMap &map1,
                                                           bool verbose) const
{
    return this->findMCSmatches(other, matcher, timeout, false, map0, map1, verbose);
}

/** Return the root mean square deviation (RMSD) of the atoms in this view against
    the atoms in 'other', using the passed AtomMatcher to match atoms in this
    view against 'other', and using the passed property maps to find the required
    properties */
SireUnits::Dimension::Length Evaluator::rmsd(const MoleculeView &other,
                                             const AtomMatcher &atommatcher,
                                             const PropertyMap &map0,
                                             const PropertyMap &map1) const
{
    const AtomCoords &c0 = this->data().property( map0["coordinates"] ).asA<AtomCoords>();
    const AtomCoords &c1 = other.data().property( map1["coordinates"] ).asA<AtomCoords>();

    QHash<AtomIdx,AtomIdx> match = atommatcher.match(*this, map0, other, map1);

    const AtomSelection &sel0 = this->selection();
    const AtomSelection &sel1 = other.selection();

    Average msd;

    for (QHash<AtomIdx,AtomIdx>::const_iterator it = match.constBegin();
         it != match.constEnd();
         ++it)
    {
        const AtomIdx atm0 = it.key();
        const AtomIdx atm1 = it.value();

        if (sel0.selected(atm0) and sel1.selected(atm1))
        {
            Vector v0 = c0.get( this->data().info().cgAtomIdx(atm0) );
            Vector v1 = c1.get( this->data().info().cgAtomIdx(atm1) );

            msd.accumulate( Vector::distance2(v0,v1) );
        }
    }

    return Length( std::sqrt( msd.average() ) );
}

/** Return the root mean square deviation (RMSD) of the atoms in this view against
    the atoms in 'other', using the passed property map to find the required
    properties */
SireUnits::Dimension::Length Evaluator::rmsd(const MoleculeView &other,
                                             const PropertyMap &map) const
{
    return this->rmsd(other, AtomIdxMatcher(), map, map);
}

/** Return the root mean square deviation (RMSD) of the atoms in this view against
    the atoms in 'other', using the passed property maps to find the required
    properties */
SireUnits::Dimension::Length Evaluator::rmsd(const MoleculeView &other,
                                             const PropertyMap &map0,
                                             const PropertyMap &map1) const
{
    return this->rmsd(other, AtomIdxMatcher(), map0, map1);
}

/** Return the root mean square deviation (RMSD) of the atoms in this view against
    the atoms in 'other', using the passed AtomMatcher to match atoms in this
    view against 'other', and using the passed property map to find the required
    properties */
SireUnits::Dimension::Length Evaluator::rmsd(const MoleculeView &other,
                                             const AtomMatcher &atommatcher,
                                             const PropertyMap &map) const
{
    return this->rmsd(other, atommatcher, map, map);
}

const char* Evaluator::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Evaluator>() );
}

Evaluator* Evaluator::clone() const
{
    return new Evaluator(*this);
}
