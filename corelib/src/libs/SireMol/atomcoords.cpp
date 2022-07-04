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

#include "atomcoords.h"

#include "SireMaths/align.h"

#include "SireVol/space.h"

#include "SireBase/quickcopy.hpp"
#include "SireMaths/vectorproperty.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireStream;
using namespace SireVol;
using namespace SireMaths;
using namespace SireMol;
using namespace SireBase;

static const RegisterMetaType< AtomProperty<Vector> > r_atomcoords;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                       const AtomProperty<Vector> &atomcoords)
{
    writeHeader(ds, r_atomcoords, 1);

    SharedDataStream sds(ds);

    sds << atomcoords.coords
        << static_cast<const MolViewProperty&>(atomcoords);

    return ds;
}

/** Deserialise from a binary datastream */
QDataStream &operator>>(QDataStream &ds,
                                       AtomProperty<Vector> &atomcoords)
{
    VersionID v = readHeader(ds, r_atomcoords);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> atomcoords.coords
            >> static_cast<MolViewProperty&>(atomcoords);
    }
    else
        throw version_error(v, "1", r_atomcoords, CODELOC);

    return ds;
}

/** Null constructor */
AtomProperty<Vector>::AtomProperty()
                     : ConcreteProperty<AtomProperty<Vector>,AtomProp>()
{}

/** Create space for the coordinates for all of the atoms described in
    'molinfo' - the coordinates are all initially (0,0,0) */
AtomProperty<Vector>::AtomProperty(const MoleculeInfoData &molinfo)
                     : ConcreteProperty<AtomProperty<Vector>,AtomProp>()
{
    int ncg = molinfo.nCutGroups();

    if (ncg > 0)
    {
        QVector<CoordGroup> tmp_coords(ncg);
        CoordGroup *tmp_coords_array = tmp_coords.data();

        for (CGIdx i(0); i<ncg; ++i)
        {
            tmp_coords_array[i] = CoordGroup(molinfo.nAtoms(i));
        }

        coords = CoordGroupArray(tmp_coords);
    }
}

/** Construct from the passed set of coordinates (arranged into
    a single CutGroup */
AtomProperty<Vector>::AtomProperty(const CoordGroup &cgroup)
                     : ConcreteProperty<AtomProperty<Vector>,AtomProp>(),
                       coords(cgroup)
{}

/** Construct from the passed set of coordinates (arranged into
    CutGroups) */
AtomProperty<Vector>::AtomProperty(const CoordGroupArray &cgroups)
                     : ConcreteProperty<AtomProperty<Vector>,AtomProp>(),
                       coords(cgroups)
{}

/** Copy constructor */
AtomProperty<Vector>::AtomProperty(const AtomProperty<Vector> &other)
                     : ConcreteProperty<AtomProperty<Vector>,AtomProp>(other),
                       coords(other.coords)
{}

/** Destructor */
AtomProperty<Vector>::~AtomProperty()
{}

/** Copy assignment operator */
AtomProperty<Vector>& AtomProperty<Vector>::operator=(const AtomProperty<Vector> &other)
{
    MolViewProperty::operator=(other);
    coords = other.coords;
    return *this;
}

/** Comparison operator */
bool AtomProperty<Vector>::operator==(const AtomProperty<Vector> &other) const
{
    return coords == other.coords;
}

/** Comparison operator */
bool AtomProperty<Vector>::operator!=(const AtomProperty<Vector> &other) const
{
    return coords != other.coords;
}

/** Return the CoordGroup for the CutGroup at index 'cgidx' */
const CoordGroup& AtomProperty<Vector>::operator[](CGIdx cgidx) const
{
    return coords.constData()[cgidx.map(coords.count())];
}

/** Merge all of the atomic properties into a single array, with
    the properties arranged in AtomIdx order */
PropertyPtr AtomProperty<Vector>::merge(const MoleculeInfoData &moldata) const
{
    this->assertCompatibleWith(moldata);

    QVector<Vector> vals( moldata.nAtoms() );

    Vector *vals_array = vals.data();

    for (AtomIdx i(0); i<moldata.nAtoms(); ++i)
    {
        vals_array[i] = this->at( moldata.cgAtomIdx(i) );
    }

    return AtomProperty<Vector>(vals);
}

/** Divide the AtomProperty into beads according to the passed atom selections,
    and returning the properties in AtomIdx order within each bead

    \throw SireError::incompatible_error
*/
PropertyPtr AtomProperty<Vector>::divide(const QVector<AtomSelection> &beads) const
{
    if (beads.isEmpty())
        return PropertyPtr();

    const int nbeads = beads.count();
    const AtomSelection *beads_array = beads.constData();

    QVector< QVector<Vector> > bead_vals(nbeads);
    QVector<Vector> *bead_vals_array = bead_vals.data();

    for (int i=0; i<nbeads; ++i)
    {
        const AtomSelection &bead = beads_array[i];

        bead.assertCompatibleWith<Vector>(*this);

        QVector<Vector> vals( bead.nSelected() );
        Vector *vals_array = vals.data();

        if (bead.selectedAll())
        {
            for (AtomIdx j(0); j<bead.nSelected(); ++j)
            {
                vals_array[j] = this->at( bead.info().cgAtomIdx(j) );
            }
        }
        else
        {
            foreach (const AtomIdx &j, bead.selectedAtoms())
            {
                *vals_array = this->at( bead.info().cgAtomIdx(j) );
                ++vals_array;
            }
        }

        bead_vals_array[i] = vals;
    }

    return AtomProperty<Vector>(bead_vals);
}

/** Divide the properties into residues. This returns the values in
    Residue/Index order

    \throw SireError::incompatible_error
*/
PropertyPtr AtomProperty<Vector>::divideByResidue(const MoleculeInfoData &molinfo) const
{
    this->assertCompatibleWith(molinfo);

    QVector< QVector<Vector> > res_vals( molinfo.nResidues() );
    QVector<Vector> *res_vals_array = res_vals.data();

    for (ResIdx i(0); i<molinfo.nResidues(); ++i)
    {
        const int nats = molinfo.nAtoms(i);

        QVector<Vector> vals(nats);
        Vector *vals_array = vals.data();

        for (int j=0; j<nats; ++j)
        {
            vals_array[j] = this->at( molinfo.cgAtomIdx(molinfo.getAtom(i,j)) );
        }

        res_vals_array[i] = vals;
    }

    return AtomProperty<Vector>(res_vals);
}

/** Convert the coordinates to an array of array of vectors - each
    array holds the coordinates of one CutGroup */
AtomProperty<QVariant> AtomProperty<Vector>::toVariant() const
{
    if (coords.count() == 0)
        return PackedArray2D<QVariant>();

    int ngroups = coords.count();

    QVector< QVector<QVariant> > converted_coords(ngroups);
    converted_coords.squeeze();

    const CoordGroup *coords_array = coords.constData();
    QVector<QVariant> *converted_coords_array = converted_coords.data();

    for (int i=0; i<ngroups; ++i)
    {
        const CoordGroup &cgroup = coords_array[i];
        int nats = cgroup.count();

        if (nats == 0)
            converted_coords_array[i] = QVector<QVariant>();
        else
        {
            const Vector *cgroup_array = cgroup.constData();

            QVector<QVariant> converted_vals(nats);
            converted_vals.squeeze();
            QVariant *converted_vals_array = converted_vals.data();

            for (int j=0; j<nats; ++j)
            {
                converted_vals_array[j].setValue<Vector>(cgroup_array[j]);
            }

            converted_coords_array[i] = converted_vals;
        }
    }

    return AtomProperty<QVariant>( PackedArray2D<QVariant>(converted_coords) );
}

static void assertCanConvert(const QVariant &value)
{
    if (not value.canConvert<Vector>())
    {
        throw SireError::invalid_cast( QObject::tr(
            "Cannot convert an object of type %1 to a SireMaths::Vector, "
            "as is required to allow this object to be part of an "
            "AtomCoords property.")
                .arg(value.typeName()), CODELOC );
    }
}

static CoordGroup makeCoordGroup(const PackedArray2D<QVariant>::Array &values)
{
    if (values.isEmpty())
        return CoordGroup();

    int nvals = values.count();
    const QVariant *values_array = values.constData();

    QVector<Vector> tmp_coords(nvals);

    for (int i=0; i<nvals; ++i)
    {
        const QVariant &value = values_array[i];
        assertCanConvert(value);

        tmp_coords[i] = value.value<Vector>();
    }

    return CoordGroup(tmp_coords);
}

/** Return a set of AtomCoords that have been converted from an array of QVariants */
AtomProperty<Vector> AtomProperty<Vector>::fromVariant(
                                        const AtomProperty<QVariant> &variantprop)
{
    const PackedArray2D<QVariant> &variant = variantprop.array();

    if (variant.isEmpty())
        return AtomProperty<Vector>();

    int ngroups = variant.count();

    QVector<CoordGroup> coords(ngroups);
    const PackedArray2D<QVariant>::Array *variant_array = variant.constData();

    for (int i=0; i<ngroups; ++i)
    {
        coords[i] = ::makeCoordGroup( variant_array[i] );
    }

    return AtomProperty<Vector>( CoordGroupArray(coords) );
}

void AtomProperty<Vector>::assignFrom(const AtomProperty<QVariant> &variant)
{
    this->operator=( AtomProperty<Vector>::fromVariant(variant) );
}

/** Return whether or not this set of coordinates is compatible with
    the molecule layout in 'molinfo' */
bool AtomProperty<Vector>::isCompatibleWith(const MoleculeInfoData &molinfo) const
{
    int ncg = molinfo.nCutGroups();

    if (ncg != coords.count())
        return false;

    const CoordGroup *coords_array = coords.constData();

    for (CGIdx i(0); i<ncg; ++i)
    {
        if (molinfo.nAtoms(i) != coords_array[i].count())
            return false;
    }

    return true;
}

/** Return whether or not the passed QVariant is compatible with
    this property */
bool AtomProperty<Vector>::canConvert(const QVariant &value) const
{
    return value.canConvert<Vector>();
}

/** Assert that the passed QVariant is compatible with this property */
void AtomProperty<Vector>::assertCanConvert(const QVariant &value) const
{
    return ::assertCanConvert(value);
}

/** Return the CoordGroup for the CutGroup at index 'cgidx' */
const CoordGroup& AtomProperty<Vector>::at(CGIdx cgidx) const
{
    return this->operator[](cgidx);
}

/** Return the CoordGroup for the CutGroup at index 'cgidx' */
const CoordGroup& AtomProperty<Vector>::get(CGIdx cgidx) const
{
    return this->operator[](cgidx);
}

/** Return the coordinates of the atom at index 'cgatomidx' */
const Vector& AtomProperty<Vector>::operator[](const CGAtomIdx &cgatomidx) const
{
    return this->operator[](cgatomidx.cutGroup()).at(cgatomidx.atom());
}

/** Return the coordinates of the atom at index 'cgatomidx' */
const Vector& AtomProperty<Vector>::at(const CGAtomIdx &cgatomidx) const
{
    return this->operator[](cgatomidx);
}

/** Return the coordinates of the atom at index 'cgatomidx' */
const Vector& AtomProperty<Vector>::get(const CGAtomIdx &cgatomidx) const
{
    return this->operator[](cgatomidx);
}

/** Return the coordinates as a QVariant */
QVariant AtomProperty<Vector>::getAsVariant(const CGAtomIdx &cgatomidx) const
{
    const Vector &c = this->get(cgatomidx);
    return QVariant::fromValue(c);
}

/** Retunr the coordinates as a Property */
PropertyPtr AtomProperty<Vector>::getAsProperty(const CGAtomIdx &cgatomidx) const
{
    return PropertyPtr(VectorProperty(this->get(cgatomidx)));
}

/** Set the coordinates of the atom at index 'cgatomidx' to the value 'value' */
AtomProperty<Vector>& AtomProperty<Vector>::set(const CGAtomIdx &cgatomidx,
                                                const Vector &value)
{
    quint32 cgidx = cgatomidx.cutGroup().map(coords.count());

    const CoordGroup &cgroup = coords.at(cgidx);
    quint32 atomidx = cgatomidx.atom().map(cgroup.count());

    coords.update( cgidx, cgroup.edit().setCoordinates(atomidx, value) );

    return *this;
}

/** Set the coordinates of the atoms in the CutGroup at index 'cgidx'
    to have the value 'values' - note that the number of values
    must equal the number of atoms!!! */
AtomProperty<Vector>& AtomProperty<Vector>::set(CGIdx cgidx,
                                                const QVector<Vector> &values)
{
    quint32 i = cgidx.map(coords.count());

    const CoordGroup &group = coords.constData()[i];

    if (group.count() != values.count())
        throw SireError::incompatible_error( QObject::tr(
            "Cannot set the coordinates of the atoms in the CutGroup as "
            "the number of coordinates (%1) is not the same as the number "
            "of atoms! (%2)")
                .arg(values.count()).arg(group.count()), CODELOC );

    coords.update(i, values);

    return *this;
}

/** Set the coordinates of atoms in the CutGroup at index 'cgidx' from
    the passed CoordGroup - note that the number of values must
    equal the number of atoms!!! */
AtomProperty<Vector>& AtomProperty<Vector>::set(CGIdx cgidx, const CoordGroup &cgroup)
{
    quint32 i = cgidx.map(coords.count());

    const CoordGroup &group = coords.data()[i];

    if (group.count() != cgroup.count())
        throw SireError::incompatible_error( QObject::tr(
            "Cannot set the coordinates of the atoms in the CutGroup as "
            "the number of coordinates (%1) is not the same as the number "
            "of atoms! (%2)")
                .arg(cgroup.count()).arg(group.count()), CODELOC );

    coords.update(i, cgroup);

    return *this;
}

/** Translate all of the atoms in this container by 'delta' */
void AtomProperty<Vector>::translate(const Vector &delta)
{
    coords.translate(delta);
}

/** Translate all of the atoms in the CutGroup at index 'cgidx' by 'delta'

    \throw SireError::invalid_index
*/
void AtomProperty<Vector>::translate(CGIdx cgidx, const Vector &delta)
{
    coords.translate( cgidx.map(coords.count()), delta );
}

/** Rotate all of the atoms in this container using the quaternion 'quat'
    about the point 'point' */
void AtomProperty<Vector>::rotate(const Quaternion &quat, const Vector &point)
{
    coords.rotate(quat, point);
}

/** Transform all of the atoms in this container using the transformation 't' */
void AtomProperty<Vector>::transform(const Transform &t)
{
    coords.transform(t);
}

/** Rotate all of the atoms in this container using the matrix 'rotmat'
    about the point 'point' */
void AtomProperty<Vector>::rotate(const Matrix &rotmat, const Vector &point)
{
    coords.rotate(rotmat, point);
}

/** Rotate all of the atoms in the CutGroup at index 'cgidx' using
    the quaternion 'quat' about the point 'point'

    \throw SireError::invalid_index
*/
void AtomProperty<Vector>::rotate(CGIdx cgidx, const Quaternion &quat,
                                  const Vector &point)
{
    coords.rotate( cgidx.map(coords.count()), quat, point );
}

/** Transform all of the atoms in the CutGroup at index 'cgidx' using
    the transformation 't'

    \throw SireError::invalid_index
*/
void AtomProperty<Vector>::transform(CGIdx cgidx, const Transform &t)
{
    coords.transform( cgidx.map(coords.count()), t );
}

/** Rotate all of the atoms in the CutGroup at index 'cgidx' using
    the matrix 'rotmat' about the point 'point'

    \throw SireError::invalid_index
*/
void AtomProperty<Vector>::rotate(CGIdx cgidx, const Matrix &rotmat,
                                  const Vector &point)
{
    coords.rotate( cgidx.map(coords.count()), rotmat, point );
}

/** Map all of the atoms in this container into the coordinate
    space represented by 'axes' */
void AtomProperty<Vector>::mapInto(const AxisSet &axes)
{
    coords.mapInto(axes);
}

/** Map all of the atoms in the CutGroup at index 'cgidx' into
    the coordinate space represented by 'axes'

    \throw SireError::invalid_index
*/
void AtomProperty<Vector>::mapInto(CGIdx cgidx, const AxisSet &axes)
{
    coords.mapInto( cgidx.map(coords.count()), axes );
}

/** Change all of the atoms in this container from residing in the
    coordinate frame 'from_frame' and move them into the coordinate
    frame 'to_frame' */
void AtomProperty<Vector>::changeFrame(const AxisSet &from_frame,
                                       const AxisSet &to_frame)
{
    coords.changeFrame(from_frame, to_frame);
}

/** Change all of the atoms in the CutGroup at index 'cgidx' so that
    they move from residing in the coordinate frame 'from_frame' into
    the coordinate frame 'to_frame'

    \throw SireError::invalid_index
*/
void AtomProperty<Vector>::changeFrame(CGIdx cgidx, const AxisSet &from_frame,
                                       const AxisSet &to_frame)
{
    coords.changeFrame( cgidx.map(coords.count()), from_frame, to_frame );
}

/** Return a raw pointer to the array of CoordGroups */
const CoordGroup* AtomProperty<Vector>::data() const
{
    return coords.constData();
}

/** Return a raw pointer to the array of CoordGroups */
const CoordGroup* AtomProperty<Vector>::constData() const
{
    return coords.constData();
}

/** Return a raw pointer to the array of coordinates for the
    CoordGroup at index 'cgidx' */
const Vector* AtomProperty<Vector>::data(CGIdx cgidx) const
{
    return this->operator[](cgidx).constData();
}

/** Return a raw pointer to the array of coordinates for the
    CoordGroup at index 'cgidx' */
const Vector* AtomProperty<Vector>::constData(CGIdx cgidx) const
{
    return this->data(cgidx);
}

/** Return the number of CoordGroups in this set */
int AtomProperty<Vector>::size() const
{
    return coords.count();
}

/** Return the number of CoordGroups in this set */
int AtomProperty<Vector>::count() const
{
    return coords.count();
}

/** Return the number of CoordGroups in this set */
int AtomProperty<Vector>::nCutGroups() const
{
    return coords.count();
}

/** Return the number of atoms in this set */
int AtomProperty<Vector>::nAtoms() const
{
    int nats = 0;

    for (int i=0; i<coords.count(); ++i)
    {
        nats += coords.constData()[i].count();
    }

    return nats;
}

/** Return the number of atoms in the CutGroup at index 'cgidx' */
int AtomProperty<Vector>::nAtoms(CGIdx cgidx) const
{
    return this->operator[](cgidx).count();
}

/** Convert this atom property to an array of values. The values
    are written in CGAtomIdx order */
QVector<Vector> AtomProperty<Vector>::toVector() const
{
    int nats = this->nAtoms();

    if (nats == 0)
        return QVector<Vector>();

    QVector<Vector> vals(nats);

    quickCopy<Vector>( vals.data(), coords.constCoordsData(), nats );

    return vals;
}

/** Convert the properties of the atoms selected in 'selection' to an
    array of values. The values are written in CGAtomIdx order

    \throw SireError::incompatible_error
*/
QVector<Vector> AtomProperty<Vector>::toVector(const AtomSelection &selected_atoms) const
{
    selected_atoms.assertCompatibleWith(*this);

    if (selected_atoms.selectedAll())
        return this->toVector();

    else if (selected_atoms.selectedAllCutGroups())
    {
        int nselected = selected_atoms.nSelected();
        QVector<Vector> vals(nselected);

        Vector *value = vals.data();

        const int ncg = selected_atoms.nCutGroups();

        const CoordGroup *cgroup_array = coords.constData();

        for (CGIdx i(0); i<ncg; ++i)
        {
            const Vector *group_coords = cgroup_array[i].constData();

            if (selected_atoms.selectedAll(i))
            {
                const int nats = cgroup_array[i].count();

                quickCopy(value, group_coords, nats);
                value += nats;
            }
            else
            {
                QList<Index> idxs = selected_atoms.selectedAtoms(i).values();
                std::sort(idxs.begin(), idxs.end());

                foreach (Index idx, idxs)
                {
                    *value = group_coords[idx];
                    ++value;
                }
            }
        }

        return vals;
    }
    else
    {
        int nselected = selected_atoms.nSelected();
        QVector<Vector> vals(nselected);

        Vector *value = vals.data();

        const CoordGroup *cgroup_array = coords.constData();

        QList<CGIdx> cgidxs = selected_atoms.selectedCutGroups();
        std::sort(cgidxs.begin(), cgidxs.end());

        foreach (CGIdx i, cgidxs)
        {
            const Vector *group_coords = cgroup_array[i].constData();

            if (selected_atoms.selectedAll(i))
            {
                const int nats = cgroup_array[i].count();

                quickCopy(value, group_coords, nats);
                value += nats;
            }
            else
            {
                QList<Index> idxs = selected_atoms.selectedAtoms(i).values();
                std::sort(idxs.begin(), idxs.end());

                foreach (Index idx, idxs)
                {
                    *value = group_coords[idx];
                    ++value;
                }
            }
        }

        return vals;
    }
}

/** Copy into this atom property set the values from 'values'. The values
    are copied in CGAtomIdx order, and there must be as many values
    as there are atoms

    \throw SireError::incompatible_error
*/
void AtomProperty<Vector>::copyFrom(const QVector<Vector> &values)
{
    if (values.count() != this->nAtoms())
        this->throwIncorrectNumberOfAtoms(values.count(), this->nAtoms());

    const int ncg = coords.nCoordGroups();

    const Vector *value = values.constData();

    for (int i=0; i<ncg; ++i)
    {
        const int nats = coords.at(i).count();

        coords.update(i, value, nats);
        value += nats;
    }
}

/** Copy into this atom property set the values from 'values', but only
    for the atoms selected in 'selection'. This copies the properties
    in in CGAtomIdx order, and there must be the same number of values
    as there are selected atoms

    \throw SireError::incompatible_error
*/
void AtomProperty<Vector>::copyFrom(const QVector<Vector> &values,
                                    const AtomSelection &selected_atoms)
{
    selected_atoms.assertCompatibleWith(*this);

    if (selected_atoms.selectedAll())
    {
        this->copyFrom(values);
        return;
    }

    if (values.count() != selected_atoms.nSelected())
        this->throwIncorrectNumberOfSelectedAtoms(values.count(),
                                                  selected_atoms.nSelected());

    if (selected_atoms.selectedAllCutGroups())
    {
        const Vector *value = values.constData();

        const int ncg = coords.nCoordGroups();

        for (CGIdx i(0); i<ncg; ++i)
        {
            if (selected_atoms.selectedAll(i))
            {
                const int nats = coords.at(i).count();

                coords.update(i, value, nats);
                value += nats;
            }
            else
            {
                QVector<Vector> new_coords = coords.at(i).toVector();

                QList<Index> idxs = selected_atoms.selectedAtoms(i).values();
                std::sort(idxs.begin(), idxs.end());

                foreach (Index idx, idxs)
                {
                    new_coords[idx] = *value;
                    ++value;
                }

                coords.update(i, new_coords);
            }
        }
    }
    else
    {
        QList<CGIdx> cgidxs = selected_atoms.selectedCutGroups();
        std::sort(cgidxs.begin(), cgidxs.end());

        const Vector *value = values.constData();

        foreach (CGIdx i, cgidxs)
        {
            if (selected_atoms.selectedAll(i))
            {
                const int nats = coords.at(i).count();

                coords.update(i, value, nats);
                value += nats;
            }
            else
            {
                QVector<Vector> new_coords = coords.at(i).toVector();

                QList<Index> idxs = selected_atoms.selectedAtoms(i).values();
                std::sort(idxs.begin(), idxs.end());

                foreach (Index idx, idxs)
                {
                    new_coords[idx] = *value;
                    ++value;
                }

                coords.update(i, new_coords);
            }
        }
    }
}

const char* AtomCoords::typeName()
{
    return QMetaType::typeName( qMetaTypeId<AtomCoords>() );
}

AtomProperty<Vector>* AtomCoords::clone() const
{
    return new AtomProperty<Vector>(*this);
}
