/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2010  Christopher Woods
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

#include "fieldtable.h"

#include "SireMol/moleculeview.h"
#include "SireMol/moleculegroup.h"
#include "SireMol/atomselection.h"

#include "SireMol/mover.hpp"

#include "SireID/index.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireMol/errors.h"
#include "SireError/errors.h"

using namespace SireFF;
using namespace SireMol;
using namespace SireMaths;
using namespace SireBase;
using namespace SireID;
using namespace SireStream;

/////////
///////// Implementation of MolFieldTable
/////////

static const RegisterMetaType<MolFieldTable> r_moltable(NO_ROOT);

QDataStream &operator<<(QDataStream &ds,
                                      const MolFieldTable &moltable)
{
    writeHeader(ds, r_moltable, 1);

    SharedDataStream sds(ds);

    sds << moltable.molnum << moltable.moluid << moltable.ncgroups
        << moltable.cgidx_to_idx << static_cast<const PackedArray2D<Vector>&>(moltable);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, MolFieldTable &moltable)
{
    VersionID v = readHeader(ds, r_moltable);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> moltable.molnum >> moltable.moluid >> moltable.ncgroups
            >> moltable.cgidx_to_idx >> static_cast<PackedArray2D<Vector>&>(moltable);
    }
    else
        throw version_error( v, "1", r_moltable, CODELOC );

    return ds;
}

/** Null constructor */
MolFieldTable::MolFieldTable() : PackedArray2D<Vector>(), molnum(0), ncgroups(0)
{}

/** Construct to hold the field acting at the points of all of the atoms
    of all of the cutgroups viewed in 'molview' */
MolFieldTable::MolFieldTable(const MoleculeView &molview)
              : PackedArray2D<Vector>(),
                molnum(molview.data().number()),
                moluid(molview.data().info().UID()),
                ncgroups(molview.data().info().nCutGroups())
{
    //build arrays for each selected CutGroup
    AtomSelection selected_atoms = molview.selection();

    if (selected_atoms.selectedAllCutGroups())
    {
        QVector< QVector<Vector> > fields(ncgroups);
        QVector<Vector> *fields_array = fields.data();

        for (CGIdx i(0); i<ncgroups; ++i)
        {
            fields_array[i] = QVector<Vector>(molview.data().info().nAtoms(i),
                                              Vector(0));
        }

        PackedArray2D<Vector>::operator=(fields);
    }
    else
    {
        QVector< QVector<Vector> > fields(selected_atoms.nSelectedCutGroups());
        cgidx_to_idx.reserve(selected_atoms.nSelectedCutGroups());

        QVector<Vector> *fields_array = fields.data();
        qint32 idx = 0;

        foreach (CGIdx i, selected_atoms.selectedCutGroups())
        {
            fields_array[i] = QVector<Vector>(molview.data().info().nAtoms(i),
                                              Vector(0));

            cgidx_to_idx.insert(i, idx);
            ++idx;
        }

        PackedArray2D<Vector>::operator=(fields);
    }
}

/** Copy constructor */
MolFieldTable::MolFieldTable(const MolFieldTable &other)
              : PackedArray2D<Vector>(other),
                molnum(other.molnum), moluid(other.moluid),
                ncgroups(other.ncgroups), cgidx_to_idx(other.cgidx_to_idx)
{}

/** Destructor */
MolFieldTable::~MolFieldTable()
{}

/** Copy assignment operator */
MolFieldTable& MolFieldTable::operator=(const MolFieldTable &other)
{
    if (this != &other)
    {
        PackedArray2D<Vector>::operator=(other);
        molnum = other.molnum;
        moluid = other.moluid;
        ncgroups = other.ncgroups;
        cgidx_to_idx = other.cgidx_to_idx;
    }

    return *this;
}

/** Set the field at all points equal to 'field' */
MolFieldTable& MolFieldTable::operator=(const Vector &field)
{
    this->setAll(field);
    return *this;
}

/** Comparison operator */
bool MolFieldTable::operator==(const MolFieldTable &other) const
{
    return this == &other or
           (molnum == other.molnum and moluid == other.moluid and
            ncgroups == other.ncgroups and cgidx_to_idx == other.cgidx_to_idx and
            PackedArray2D<Vector>::operator==(other));
}

/** Comparison operator */
bool MolFieldTable::operator!=(const MolFieldTable &other) const
{
    return not MolFieldTable::operator==(other);
}

/** Add the fields in the passed table onto this table */
MolFieldTable& MolFieldTable::operator+=(const MolFieldTable &other)
{
    this->add(other);
    return *this;
}

/** Subtract the fields in the passed table from this table */
MolFieldTable& MolFieldTable::operator-=(const MolFieldTable &other)
{
    this->subtract(other);
    return *this;
}

/** Return the sum of the fields of this table with 'other' */
MolFieldTable MolFieldTable::operator+(const MolFieldTable &other) const
{
    MolFieldTable ret(*this);
    ret += other;
    return ret;
}

/** Return the difference of the fields of this table with 'other' */
MolFieldTable MolFieldTable::operator-(const MolFieldTable &other) const
{
    MolFieldTable ret(*this);
    ret -= other;
    return *this;
}

/** Add 'field' to all of the points in this table */
MolFieldTable& MolFieldTable::operator+=(const Vector &field)
{
    this->add(field);
    return *this;
}

/** Subtract 'field' from all of the points in this table */
MolFieldTable& MolFieldTable::operator-=(const Vector &field)
{
    this->subtract(field);
    return *this;
}

/** Return a table that has 'field' added to all of the points */
MolFieldTable MolFieldTable::operator+(const Vector &field) const
{
    MolFieldTable ret(*this);
    ret += field;
    return ret;
}

MolFieldTable operator+(const Vector &field, const MolFieldTable &table)
{
    return table + field;
}

/** Return a table that has 'field' subtracted from all of the points */
MolFieldTable MolFieldTable::operator-(const Vector &field) const
{
    MolFieldTable ret(*this);
    ret -= field;
    return ret;
}

/** Multiply the field at all points by the scalar 'value' */
MolFieldTable& MolFieldTable::operator*=(double value)
{
    this->multiply(value);
    return *this;
}

/** Divide the field at all points by the scalar 'value' */
MolFieldTable& MolFieldTable::operator/=(double value)
{
    this->divide(value);
    return *this;
}

/** Return the field that has been multiplied by 'value' */
MolFieldTable MolFieldTable::operator*(double value) const
{
    MolFieldTable ret(*this);
    ret *= value;
    return ret;
}

/** Return the field that has been multiplied by 'value' */
MolFieldTable operator*(double value, const MolFieldTable &table)
{
    return table * value;
}

/** Return the field that has been divided by 'value' */
MolFieldTable MolFieldTable::operator/(double value) const
{
    MolFieldTable ret(*this);
    ret /= value;
    return ret;
}

/** Return the negative of this table */
MolFieldTable MolFieldTable::operator-() const
{
    MolFieldTable ret(*this);
    ret *= -1;
    return ret;
}

const char* MolFieldTable::typeName()
{
    return QMetaType::typeName( qMetaTypeId<MolFieldTable>() );
}

/** Initialise this table - this clears all of the fields, resetting them to zero */
void MolFieldTable::initialise()
{
    int nvals = PackedArray2D<Vector>::nValues();

    if (nvals > 0)
    {
        Vector *vals = PackedArray2D<Vector>::valueData();

        for (int i=0; i<nvals; ++i)
        {
            vals[i] = Vector(0);
        }
    }
}

/** Return all of the fields in this table in a single array */
QVector<Vector> MolFieldTable::toVector() const
{
    return PackedArray2D<Vector>::toQVector();
}


/** Return an array of all of the fields at the location of
    the atoms selected in 'selection'

    \throw SireError::incompatible_error
*/
QVector<Vector> MolFieldTable::toVector(const AtomSelection &selection) const
{
    this->assertCompatibleWith(selection);

    if (selection.selectedAll())
    {
        if (not this->selectedAll())
            throw SireMol::missing_atom( QObject::tr(
                "Cannot return the forces on all atoms as not all of the atoms "
                "are selected in this forcetable."), CODELOC );

        return this->toVector();
    }

    QVector<Vector> vals( selection.nSelected() );
    Vector *value = vals.data();

    if (this->selectedAll())
    {
        if (selection.selectedAllCutGroups())
        {
            const int ncg = selection.nCutGroups();

            for (CGIdx i(0); i<ncg; ++i)
            {
                const Vector *groupfields = PackedArray2D<Vector>::constData(i);

                if (selection.selectedAll(i))
                {
                    const int nats = PackedArray2D<Vector>::nValues(i);

                    quickCopy<Vector>(value, groupfields, nats);
                    value += nats;
                }
                else
                {
                    QList<Index> idxs = selection.selectedAtoms(i).values();
                    std::sort(idxs.begin(), idxs.end());

                    foreach (Index idx, idxs)
                    {
                        *value = groupfields[idx];
                        ++value;
                    }
                }
            }
        }
        else
        {
            QList<CGIdx> cgidxs = selection.selectedCutGroups();
            std::sort(cgidxs.begin(), cgidxs.end());

            foreach (CGIdx i, cgidxs)
            {
                const Vector *groupfields = PackedArray2D<Vector>::constData(i);

                if (selection.selectedAll(i))
                {
                    const int nats = PackedArray2D<Vector>::nValues(i);

                    quickCopy<Vector>(value, groupfields, nats);
                    value += nats;
                }
                else
                {
                    QList<Index> idxs = selection.selectedAtoms(i).values();
                    std::sort(idxs.begin(), idxs.end());

                    foreach (Index idx, idxs)
                    {
                        *value = groupfields[idx];
                        ++value;
                    }
                }
            }
        }
    }
    else
    {
        if (selection.selectedAllCutGroups())
            throw SireMol::missing_atom( QObject::tr(
                "Cannot return the forces as while all CutGroups are selected, "
                "not all CutGroups are present in the forcetable."), CODELOC );

        QList<CGIdx> cgidxs = selection.selectedCutGroups();
        std::sort(cgidxs.begin(), cgidxs.end());

        foreach (CGIdx cgidx, cgidxs)
        {
            int i = cgidx_to_idx.value(cgidx, -1);

            if (i == -1)
                throw SireMol::missing_atom( QObject::tr(
                    "Cannot return the forces as while atoms in CutGroup %1 "
                    "are selected, this CutGroup is not present in the forcetable.")
                        .arg(cgidx), CODELOC );

            const Vector *groupfields = PackedArray2D<Vector>::constData(i);

            if (selection.selectedAll(cgidx))
            {
                const int nats = PackedArray2D<Vector>::nValues(i);

                quickCopy<Vector>(value, groupfields, nats);
                value += nats;
            }
            else
            {
                QList<Index> idxs = selection.selectedAtoms(cgidx).values();
                std::sort(idxs.begin(), idxs.end());

                foreach (Index idx, idxs)
                {
                    *value = groupfields[idx];
                    ++value;
                }
            }
        }
    }

    return vals;
}

/** Add the field 'field' onto this table - this returns whether or not the
    atom is in this table

    \throw SireError::invalid_index
*/
bool MolFieldTable::add(const CGAtomIdx &cgatomidx, const Vector &field)
{
    CGIdx cgidx( cgatomidx.cutGroup().map(this->nCutGroups()) );

    int i = -1;

    if (this->selectedAll())
    {
        i = cgidx;
    }
    else if (cgidx_to_idx.contains(cgidx))
    {
        i = cgidx_to_idx.value(cgidx);
    }
    else
    {
        return false;
    }

    int j = cgatomidx.atom().map( this->nValues(i) );

    this->operator()(i, j) += field;

    return true;
}

/** Subtract the field 'field' from this table - this returns whether or not the
    atom is in this table

    \throw SireError::invalid_index
*/
bool MolFieldTable::subtract(const CGAtomIdx &cgatomidx, const Vector &field)
{
    return this->add( cgatomidx, -field );
}

static void addField(const Vector &field, Vector *fields, const int nats)
{
    for (int i=0; i<nats; ++i)
    {
        fields[i] += field;
    }
}

/** Add the field 'field' onto this table for all of the atoms
    in 'selected_atoms - this returns whether
    or not any selected atoms are in this table

    \throw SireError::incompatible_error
*/
bool MolFieldTable::add(const AtomSelection &selected_atoms, const Vector &field)
{
    this->assertCompatibleWith(selected_atoms);

    if (selected_atoms.selectedNone() or this->isEmpty())
        return false;

    bool changed_atoms = false;

    if (selected_atoms.selectedAll())
    {
        //this is easy - all atoms are selected for updating,
        //so just update all of the forces in this table
        ::addField(field, this->valueData(), this->nValues());

        changed_atoms = true;
    }
    else if (this->selectedAll())
    {
        //easy(ish) case - all atoms are in this forcetable,
        //so we only need to update the forces of the selected atoms

        if (selected_atoms.selectedAllCutGroups())
        {
            for (CGIdx i(0); i<ncgroups; ++i)
            {
                if (selected_atoms.selectedAll(i))
                {
                    ::addField(field, this->data(i), this->nValues(i));
                    changed_atoms = true;
                }
                else
                {
                    QSet<Index> idxs = selected_atoms.selectedAtoms(i);

                    Vector *atomfields = this->data(i);
                    const int nats = this->nValues(i);

                    foreach (Index idx, idxs)
                    {
                        BOOST_ASSERT( idx >= 0 and idx < nats );
                        atomfields[idx] += field;
                    }

                    changed_atoms = true;
                }
            }
        }
        else
        {
            QList<CGIdx> cgidxs = selected_atoms.selectedCutGroups();

            foreach (CGIdx i, cgidxs)
            {
                if (selected_atoms.selectedAll(i))
                {
                    ::addField(field, this->data(i), this->nValues(i));
                    changed_atoms = true;
                }
                else
                {
                    QSet<Index> idxs = selected_atoms.selectedAtoms(i);

                    Vector *atomfields = this->data(i);
                    const int nats = this->nValues(i);

                    foreach (Index idx, idxs)
                    {
                        BOOST_ASSERT( idx >= 0 and idx < nats );
                        atomfields[idx] += field;
                    }

                    changed_atoms = true;
                }
            }
        }
    }
    else
    {
        //harder case - not all atoms are in this fieldtable
        //and not all atoms are selected for updating

        if (selected_atoms.selectedAllCutGroups())
        {
            for (QHash<CGIdx,qint32>::const_iterator it = cgidx_to_idx.constBegin();
                 it != cgidx_to_idx.constEnd();
                 ++it)
            {
                const CGIdx cgidx = it.key();
                const int i = it.value();

                if (selected_atoms.selectedAll(cgidx))
                {
                    ::addField(field, this->data(i), this->nValues(i));
                    changed_atoms = true;
                }
                else
                {
                    QSet<Index> idxs = selected_atoms.selectedAtoms(cgidx);

                    Vector *atomfields = this->data(i);
                    const int nats = this->nValues(i);

                    foreach (Index idx, idxs)
                    {
                        BOOST_ASSERT( idx >= 0 and idx < nats );
                        atomfields[idx] += field;
                    }
                }

                changed_atoms = true;
            }
        }
        else
        {
            for (QHash<CGIdx,qint32>::const_iterator it = cgidx_to_idx.constBegin();
                 it != cgidx_to_idx.constEnd();
                 ++it)
            {
                const CGIdx cgidx = it.key();
                const int i = it.value();

                if (selected_atoms.selectedAll(cgidx))
                {
                    ::addField(field, this->data(i), this->nValues(i));
                    changed_atoms = true;
                }
                else if (selected_atoms.selected(cgidx))
                {
                    QSet<Index> idxs = selected_atoms.selectedAtoms(cgidx);

                    Vector *atomfields = this->data(i);
                    const int nats = this->nValues(i);

                    foreach (Index idx, idxs)
                    {
                        BOOST_ASSERT( idx >= 0 and idx < nats );
                        atomfields[idx] += field;
                    }

                    changed_atoms = true;
                }
            }
        }
    }

    return changed_atoms;
}

/** Subtract the field 'field' from this table for all of the atoms
    in 'selected_atoms' - this returns whether
    or not any selected atoms are in this table

    \throw SireError::incompatible_error
*/
bool MolFieldTable::subtract(const AtomSelection &selected_atoms, const Vector &field)
{
    return MolFieldTable::add( selected_atoms, -field );
}

/** Add the field 'field' onto all of the atom points in this table */
void MolFieldTable::add(const Vector &field)
{
    int nvals = PackedArray2D<Vector>::nValues();

    if (nvals > 0)
    {
        Vector *vals = PackedArray2D<Vector>::valueData();

        for (int i=0; i<nvals; ++i)
        {
            vals[i] += field;
        }
    }
}

/** Subtract the field 'field' from all of the atom points in this table */
void MolFieldTable::subtract(const Vector &field)
{
    this->add( -field );
}

/** Multiply the field at all atom points by 'value' */
void MolFieldTable::multiply(double value)
{
    int nvals = PackedArray2D<Vector>::nValues();

    if (nvals > 0)
    {
        Vector *vals = PackedArray2D<Vector>::valueData();

        for (int i=0; i<nvals; ++i)
        {
            vals[i] *= value;
        }
    }
}

/** Divide the field at all atom points by 'value' */
void MolFieldTable::divide(double value)
{
    this->multiply( 1.0 / value );
}

/** Set all of the fields at the atom points equal to 'field' */
void MolFieldTable::setAll(const Vector &field)
{
    int nvals = PackedArray2D<Vector>::nValues();

    if (nvals > 0)
    {
        Vector *vals = PackedArray2D<Vector>::valueData();

        for (int i=0; i<nvals; ++i)
        {
            vals[i] = field;
        }
    }
}

void MolFieldTable::assertCompatibleWith(const AtomSelection &selection) const
{
    if (not selection.selectedAll())
    {
        AtomSelection new_selection(selection);
        new_selection.selectAll();
        this->assertCompatibleWith(new_selection);
        return;
    }

    bool compatible = true;

    if (selection.nCutGroups() != ncgroups)
    {
        compatible = false;
    }
    else if (this->selectedAll())
    {
        for (CGIdx i(0); i<ncgroups; ++i)
        {
            if (selection.nSelected(i) != PackedArray2D<Vector>::nValues(i))
            {
                compatible = false;
                break;
            }
        }
    }
    else
    {
        for (QHash<CGIdx,qint32>::const_iterator it = cgidx_to_idx.constBegin();
             it != cgidx_to_idx.constEnd();
             ++it)
        {
            if (selection.nSelected(it.key()) != PackedArray2D<Vector>::nValues(it.key()))
            {
                compatible = false;
                break;
            }
        }
    }

    if (not compatible)
        throw SireError::incompatible_error( QObject::tr(
            "This MolForceTable is incompatible with the passed atom selection."),
                CODELOC );
}

/** Add the fields contained in 'other' onto this field table. This will only
    add the fields for CutGroups that are in both tables */
void MolFieldTable::add(const MolFieldTable &other)
{
    if (this == &other)
    {
        //just double everything
        this->operator*=(2);
        return;
    }

    if (molnum != other.molnum)
        throw SireError::incompatible_error( QObject::tr(
                "You cannot combine the field table for molecule %1 with the "
                "field table for molecule %2. The molecules must be the same.")
                    .arg(molnum).arg(other.molnum), CODELOC );

    if (moluid != other.moluid)
        throw SireError::incompatible_error( QObject::tr(
                "You cannot combine together the tables for molecule %1 as the "
                "layout UIDs are different (%2 vs. %3). They must be the same.")
                    .arg(molnum).arg(moluid.toString(), other.moluid.toString()),
                        CODELOC );

    if (this->selectedAll() and other.selectedAll())
    {
        int nvals = PackedArray2D<Vector>::nValues();

        BOOST_ASSERT( nvals == other.nValues() );

        if (nvals > 0)
        {
            Vector *vals = PackedArray2D<Vector>::valueData();
            const Vector *other_vals = other.constValueData();

            for (int i=0; i<nvals; ++i)
            {
                vals[i] += other_vals[i];
            }
        }
    }
    else if (this->selectedAll())
    {
        for (CGIdx i(0); i<ncgroups; ++i)
        {
            int idx = other.map(i);

            if (idx != -1)
            {
                int nvals = this->nValues(i);
                BOOST_ASSERT( nvals == other.nValues(idx) );

                Vector *vals = PackedArray2D<Vector>::data(i);
                const Vector *other_vals = other.constData(idx);

                for (int j=0; j<nvals; ++j)
                {
                    vals[j] += other_vals[j];
                }
            }
        }
    }
    else
    {
        for (QHash<CGIdx,qint32>::const_iterator it = cgidx_to_idx.constBegin();
             it != cgidx_to_idx.constEnd();
             ++it)
        {
            int idx = other.map(it.key());

            if (idx != -1)
            {
                int nvals = this->nValues(it.value());
                BOOST_ASSERT( nvals == other.nValues(idx) );

                Vector *vals = PackedArray2D<Vector>::data(it.key());
                const Vector *other_vals = other.constData(idx);

                for (int j=0; j<nvals; ++j)
                {
                    vals[j] += other_vals[j];
                }
            }
        }
    }
}

/** Subtract the fields contained in 'other' from this field table. This will only
    subtract the fields for CutGroups that are in both tables */
void MolFieldTable::subtract(const MolFieldTable &other)
{
    if (this == &other)
    {
        this->setAll( Vector(0) );
        return;
    }

    this->add( -other );
}

/////////
///////// Implementation of GridFieldTable
/////////

static const RegisterMetaType<GridFieldTable> r_gridtable(NO_ROOT);

QDataStream &operator<<(QDataStream &ds, const GridFieldTable &gridtable)
{
    writeHeader(ds, r_gridtable, 1);

    SharedDataStream sds(ds);

    sds << gridtable.grd << gridtable.fieldvals;

    return ds;
}

QDataStream &operator>>(QDataStream &ds, GridFieldTable &gridtable)
{
    VersionID v = readHeader(ds, r_gridtable);

    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> gridtable.grd >> gridtable.fieldvals;
    }
    else
        throw version_error(v, "1", r_gridtable, CODELOC);

    return ds;
}

/** Null constructor */
GridFieldTable::GridFieldTable()
{}

/** Construct to hold the field at each of the points of the passed grid */
GridFieldTable::GridFieldTable(const Grid &grid) : grd(grid)
{
    if (grd.read().nPoints() > 0)
    {
        fieldvals = QVector<Vector>(grd.read().nPoints(), Vector(0));
        fieldvals.squeeze();
    }
}

/** Copy constructor */
GridFieldTable::GridFieldTable(const GridFieldTable &other)
               : grd(other.grd), fieldvals(other.fieldvals)
{}

/** Destructor */
GridFieldTable::~GridFieldTable()
{}

/** Copy assignment operator */
GridFieldTable& GridFieldTable::operator=(const GridFieldTable &other)
{
    if (this != &other)
    {
        grd = other.grd;
        fieldvals = other.fieldvals;
    }

    return *this;
}

/** Set the field at each of the grid points equal to 'field' */
GridFieldTable& GridFieldTable::operator=(const Vector &field)
{
    this->setAll(field);
    return *this;
}

/** Comparison operator */
bool GridFieldTable::operator==(const GridFieldTable &other) const
{
    return this == &other or
           (grd == other.grd and fieldvals == other.fieldvals);
}

/** Comparison operator */
bool GridFieldTable::operator!=(const GridFieldTable &other) const
{
    return not this->operator==(other);
}

/** Add the fields on the passed table onto this table. Note that this
    only adds the fields if both tables use the same grid points */
GridFieldTable& GridFieldTable::operator+=(const GridFieldTable &other)
{
    this->add(other);
    return *this;
}

/** Subtract the fields of 'other' from this table. Note that this
    only subtracts the fields if both tables use the same grid points */
GridFieldTable& GridFieldTable::operator-=(const GridFieldTable &other)
{
    this->subtract(other);
    return *this;
}

/** Return the sum of the fields of this table and 'other' - note that
    this returns the first table if the two tables use different grids */
GridFieldTable GridFieldTable::operator+(const GridFieldTable &other) const
{
    GridFieldTable ret(*this);
    ret += other;
    return ret;
}

/** Return the difference of the fields of this table and 'other' - note that
    this returns the first table if the two tables use different grids */
GridFieldTable GridFieldTable::operator-(const GridFieldTable &other) const
{
    GridFieldTable ret(*this);
    ret -= other;
    return ret;
}

/** Add the field 'field' to all of the grid points */
GridFieldTable& GridFieldTable::operator+=(const Vector &field)
{
    this->add(field);
    return *this;
}

/** Subtract the field 'field' from all of the grid points */
GridFieldTable& GridFieldTable::operator-=(const Vector &field)
{
    this->subtract(field);
    return *this;
}

/** Return the table where 'field' has been added to all of the grid points */
GridFieldTable GridFieldTable::operator+(const Vector &field) const
{
    GridFieldTable ret(*this);
    ret += field;
    return ret;
}

/** Return the table where 'field' has been subtracted from all of the grid points */
GridFieldTable GridFieldTable::operator-(const Vector &field) const
{
    GridFieldTable ret(*this);
    ret -= field;
    return ret;
}

/** Multiply the field at all grid points by 'value' */
GridFieldTable& GridFieldTable::operator*=(double value)
{
    this->multiply(value);
    return *this;
}

/** Divide the field at all grid points by 'value' */
GridFieldTable& GridFieldTable::operator/=(double value)
{
    this->divide(value);
    return *this;
}

/** Return the table where the field at all grid points has been
    multiplied by 'value' */
GridFieldTable GridFieldTable::operator*(double value) const
{
    GridFieldTable ret(*this);
    ret *= value;
    return ret;
}

/** Return the table where the field at all grid points has been
    divided by 'value' */
GridFieldTable GridFieldTable::operator/(double value) const
{
    GridFieldTable ret(*this);
    ret /= value;
    return ret;
}

/** Return the table where the field at all grid points is negated */
GridFieldTable GridFieldTable::operator-() const
{
    GridFieldTable ret(*this);

    for (QVector<Vector>::iterator it = ret.fieldvals.begin();
         it != ret.fieldvals.end();
         ++it)
    {
        *it = -(*it);
    }

    return ret;
}

/** Return a modifiable reference to the ith grid point's field value

    \throw SireError::invalid_index
*/
Vector& GridFieldTable::operator[](int i)
{
    return fieldvals[ Index(i).map(fieldvals.count()) ];
}

/** Return the field value of the ith grid point

    \throw SireError::invalid_index
*/
const Vector& GridFieldTable::operator[](int i) const
{
    return fieldvals.at( Index(i).map(fieldvals.count()) );
}

/** Return the field value of the ith grid point

    \throw SireError::invalid_index
*/
const Vector& GridFieldTable::at(int i) const
{
    return this->operator[](i);
}

const char* GridFieldTable::typeName()
{
    return QMetaType::typeName( qMetaTypeId<GridFieldTable>() );
}

/** Initialise the field at each grid point to equal 0 */
void GridFieldTable::initialise()
{
    for (QVector<Vector>::iterator it = fieldvals.begin();
         it != fieldvals.end();
         ++it)
    {
        *it = Vector(0);
    }
}

/** Return the number of grid points (and thus field values) */
int GridFieldTable::nPoints() const
{
    return fieldvals.count();
}

/** Return the number of grid points (and thus field values) */
int GridFieldTable::count() const
{
    return fieldvals.count();
}

/** Return the grid that contains the points at which the field is
    evaluated - the order of points in the grid is the same as the order
    of field values in this table */
const Grid& GridFieldTable::grid() const
{
    if (grd.constData() == 0)
        return Grid::null();
    else
        return grd.read();
}

/** Return the array of field values - the order is the same
    as the order of points in the grid */
QVector<Vector> GridFieldTable::toVector() const
{
    return fieldvals;
}

/** Add the field 'field' onto the field for the ipoint'th grid point

    \throw SireError::invalid_index
*/
void GridFieldTable::add(int ipoint, const Vector &field)
{
    fieldvals[ Index(ipoint).map(fieldvals.count()) ] += field;
}

/** Subtract the field 'field' from the field for the ipoint'th grid point

    \throw SireError::invalid_index
*/
void GridFieldTable::subtract(int ipoint, const Vector &field)
{
    fieldvals[ Index(ipoint).map(fieldvals.count()) ] -= field;
}

/** Add the field in 'other' onto that for this table - this only
    adds the field if the two grids are identical */
void GridFieldTable::add(const GridFieldTable &other)
{
    if (grd == other.grd)
    {
        int nvals = fieldvals.count();
        BOOST_ASSERT( nvals == other.fieldvals.count() );

        if (nvals > 0)
        {
            Vector *data = fieldvals.data();
            const Vector *other_data = other.fieldvals.constData();

            for (int i=0; i<nvals; ++i)
            {
                data[i] += other_data[i];
            }
        }
    }
}

/** Subtract the field in 'other' from that for this table - this only
    subtracts the field if the two grids are identical */
void GridFieldTable::subtract(const GridFieldTable &other)
{
    if (grd == other.grd)
    {
        int nvals = fieldvals.count();
        BOOST_ASSERT( nvals == other.fieldvals.count() );

        if (nvals > 0)
        {
            Vector *data = fieldvals.data();
            const Vector *other_data = other.fieldvals.constData();

            for (int i=0; i<nvals; ++i)
            {
                data[i] -= other_data[i];
            }
        }
    }
}

/** Add the field 'field' to all of the points in this table */
void GridFieldTable::add(const Vector &field)
{
    for (QVector<Vector>::iterator it = fieldvals.begin();
         it != fieldvals.end();
         ++it)
    {
        *it += field;
    }
}

/** Subtract the field 'field' from all of the points in this table */
void GridFieldTable::subtract(const Vector &field)
{
    for (QVector<Vector>::iterator it = fieldvals.begin();
         it != fieldvals.end();
         ++it)
    {
        *it -= field;
    }
}

/** Set the field at all of the points in this table equal to 'field' */
void GridFieldTable::setAll(const Vector &field)
{
    for (QVector<Vector>::iterator it = fieldvals.begin();
         it != fieldvals.end();
         ++it)
    {
        *it = field;
    }
}

/** Multiply the field at all of the points in this table by 'value' */
void GridFieldTable::multiply(double value)
{
    for (QVector<Vector>::iterator it = fieldvals.begin();
         it != fieldvals.end();
         ++it)
    {
        *it *= value;
    }
}

/** Divide the field at all of the points in this table by 'value' */
void GridFieldTable::divide(double value)
{
    this->multiply( 1 / value );
}

/** Return a raw pointer to the array of field values */
Vector* GridFieldTable::data()
{
    return fieldvals.data();
}

/** Return a raw pointer to the array of field values */
const Vector* GridFieldTable::data() const
{
    return fieldvals.constData();
}

/** Return a raw pointer to the array of field values */
const Vector* GridFieldTable::constData() const
{
    return fieldvals.constData();
}

GridFieldTable::iterator GridFieldTable::begin()
{
    return fieldvals.begin();
}

GridFieldTable::iterator GridFieldTable::end()
{
    return fieldvals.end();
}

GridFieldTable::const_iterator GridFieldTable::begin() const
{
    return fieldvals.constBegin();
}

GridFieldTable::const_iterator GridFieldTable::end() const
{
    return fieldvals.constEnd();
}

GridFieldTable::const_iterator GridFieldTable::constBegin() const
{
    return fieldvals.constBegin();
}

GridFieldTable::const_iterator GridFieldTable::constEnd() const
{
    return fieldvals.constEnd();
}

/////////
///////// Implementation of FieldTable
/////////

static const RegisterMetaType<FieldTable> r_fieldtable(NO_ROOT);

QDataStream &operator<<(QDataStream &ds, const FieldTable &fieldtable)
{
    writeHeader(ds, r_fieldtable, 1);

    SharedDataStream sds(ds);

    sds << fieldtable.moltables_by_idx << fieldtable.gridtables;

    return ds;
}

QDataStream &operator>>(QDataStream &ds, FieldTable &fieldtable)
{
    VersionID v = readHeader(ds, r_fieldtable);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> fieldtable.moltables_by_idx >> fieldtable.gridtables;

        QHash<MolNum,qint32> molnum_to_idx;
        molnum_to_idx.reserve(fieldtable.moltables_by_idx.count());

        quint32 i = 0;

        for (QVector<MolFieldTable>::const_iterator
                                    it = fieldtable.moltables_by_idx.constBegin();
             it != fieldtable.moltables_by_idx.constEnd();
             ++it)
        {
            molnum_to_idx.insert( it->molNum(), i );
            i += 1;
        }

        fieldtable.molnum_to_idx = molnum_to_idx;
    }
    else
        throw version_error( v, "1", r_fieldtable, CODELOC );

    return ds;
}

/** Null constructor */
FieldTable::FieldTable()
{}

void FieldTable::setGroup(const MoleculeGroup &molgroup)
{
    if (molgroup.isEmpty())
        return;

    int nmols = molgroup.nMolecules();

    moltables_by_idx = QVector<MolFieldTable>(nmols);
    moltables_by_idx.squeeze();

    molnum_to_idx = QHash<MolNum,qint32>();
    molnum_to_idx.reserve(nmols);

    MolFieldTable *moltables_by_idx_array = moltables_by_idx.data();

    quint32 i = 0;

    for (MoleculeGroup::const_iterator it = molgroup.constBegin();
         it != molgroup.constEnd();
         ++it)
    {
        moltables_by_idx_array[i] = MolFieldTable(*it);
        molnum_to_idx.insert(it->data().number(), i);
        ++i;
    }
}

/** Construct the table to hold the fields at the points of all
    of the atoms in the CutGroups that are viewed in the molecules
    in 'molgroup' */
FieldTable::FieldTable(const MoleculeGroup &molgroup)
{
    this->setGroup(molgroup);
}

/** Construct the table to hold the fields at all of the points
    in the passed grid */
FieldTable::FieldTable(const Grid &grid)
{
    gridtables.append( GridFieldTable(grid) );
    gridtables.squeeze();
}

/** Construct the table to hold the fields at all of the points
    of all of the passed grids */
FieldTable::FieldTable(const QVector<GridPtr> &grids)
{
    for (QVector<GridPtr>::const_iterator it = grids.constBegin();
         it != grids.constEnd();
         ++it)
    {
        if (it->constData() != 0)
        {
            if (not this->contains(it->read()))
                gridtables.append( GridFieldTable(it->read()) );
        }
    }

    gridtables.squeeze();
}

/** Construct the table to hold the fields at the points of all
    of the atoms in the CutGroups that are viewed in the molecules
    in 'molgroup', and all of the grid points in the passed grid */
FieldTable::FieldTable(const MoleculeGroup &molgroup, const Grid &grid)
{
    this->setGroup(molgroup);
    gridtables.append( GridFieldTable(grid) );
    gridtables.squeeze();
}

/** Construct the table to hold the fields at the points of all
    of the atoms in the CutGroups that are viewed in the molecules
    in 'molgroup', and all of the grid points in the passed grids */
FieldTable::FieldTable(const MoleculeGroup &molgroup, const QVector<GridPtr> &grids)
{
    this->setGroup(molgroup);

    for (QVector<GridPtr>::const_iterator it = grids.constBegin();
         it != grids.constEnd();
         ++it)
    {
        if (it->constData() != 0)
        {
            if (not this->contains(it->read()))
                gridtables.append( GridFieldTable(it->read()) );
        }
    }

    gridtables.squeeze();
}

/** Copy constructor */
FieldTable::FieldTable(const FieldTable &other)
           : moltables_by_idx(other.moltables_by_idx),
             gridtables(other.gridtables),
             molnum_to_idx(other.molnum_to_idx)
{}

/** Destructor */
FieldTable::~FieldTable()
{}

const char* FieldTable::typeName()
{
    return QMetaType::typeName( qMetaTypeId<FieldTable>() );
}

/** Copy assignment operator */
FieldTable& FieldTable::operator=(const FieldTable &other)
{
    if (this != &other)
    {
        moltables_by_idx = other.moltables_by_idx;
        gridtables = other.gridtables;
        molnum_to_idx = other.molnum_to_idx;
    }

    return *this;
}

/** Set the field at all points in this table equal to 'field' */
FieldTable& FieldTable::operator=(const Vector &field)
{
    this->setAll(field);
    return *this;
}

/** Comparison operator */
bool FieldTable::operator==(const FieldTable &other) const
{
    return this == &other or
           (moltables_by_idx == other.moltables_by_idx and
            gridtables == other.gridtables);
}

/** Comparison operator */
bool FieldTable::operator!=(const FieldTable &other) const
{
    return not this->operator==(other);
}

/** Add the fields from 'other' onto this table. This only adds the fields
    for molecules / grids that are in both tables */
FieldTable& FieldTable::operator+=(const FieldTable &other)
{
    this->add(other);
    return *this;
}

/** Subtract the fields from 'other' from this table. This only subtracts
    the fields for molecules / grids that are in both tables */
FieldTable& FieldTable::operator-=(const FieldTable &other)
{
    this->subtract(other);
    return *this;
}

/** Return the sum of this table with 'other' - this only adds the
    fields from 'other' to this table for molecules / grids that are
    in both tables */
FieldTable FieldTable::operator+(const FieldTable &other) const
{
    FieldTable ret(*this);
    ret += other;
    return ret;
}

/** Return the difference of this table with 'other' - this only subtracts the
    fields from 'other' to this table for molecules / grids that are
    in both tables */
FieldTable FieldTable::operator-(const FieldTable &other) const
{
    FieldTable ret(*this);
    ret -= other;
    return ret;
}

/** Add the field 'field' to all of the atom and grid points in this table */
FieldTable& FieldTable::operator+=(const Vector &field)
{
    this->add(field);
    return *this;
}

/** Substract the field 'field' from all of the atom and grid
    points in this table */
FieldTable& FieldTable::operator-=(const Vector &field)
{
    this->subtract(field);
    return *this;
}

/** Return the result of adding 'field' onto all of the atom
    and grid points in this table */
FieldTable FieldTable::operator+(const Vector &field) const
{
    FieldTable ret(*this);
    ret += field;
    return ret;
}

/** Return the result of subtracting 'field' from all of the atom
    and grid points in this table */
FieldTable FieldTable::operator-(const Vector &field) const
{
    FieldTable ret(*this);
    ret -= field;
    return ret;
}

/** Multiply the fields at all points in this table by 'value' */
FieldTable& FieldTable::operator*=(double value)
{
    this->multiply(value);
    return *this;
}

/** Divide the fields at all points in this table by 'value' */
FieldTable& FieldTable::operator/=(double value)
{
    this->divide(value);
    return *this;
}

/** Return the result of multiplying the fields at all points by 'value' */
FieldTable FieldTable::operator*(double value) const
{
    FieldTable ret(*this);
    ret *= value;
    return ret;
}

/** Return the result of dividing the fields at all points by 'value' */
FieldTable FieldTable::operator/(double value) const
{
    FieldTable ret(*this);
    ret /= value;
    return ret;
}

/** Return the result of negating the field at all points */
FieldTable FieldTable::operator-() const
{
    FieldTable ret(*this);

    for (QVector<MolFieldTable>::iterator it = ret.moltables_by_idx.begin();
         it != ret.moltables_by_idx.end();
         ++it)
    {
        *it = -(*it);
    }

    for (QVector<GridFieldTable>::iterator it = ret.gridtables.begin();
         it != ret.gridtables.end();
         ++it)
    {
        *it = -(*it);
    }

    return ret;
}

/** Return whether or not this contains a table for the passed grid */
bool FieldTable::contains(const Grid &grid) const
{
    for (QVector<GridFieldTable>::const_iterator it = gridtables.constBegin();
         it != gridtables.constEnd();
         ++it)
    {
        if (it->grid().equals(grid))
            return true;
    }

    return false;
}

/** Initialise all of the tables to have a zero field */
void FieldTable::initialiseTables()
{
    for (QVector<MolFieldTable>::iterator it = moltables_by_idx.begin();
         it != moltables_by_idx.end();
         ++it)
    {
        it->initialise();
    }

    for (QVector<GridFieldTable>::iterator it = gridtables.begin();
         it != gridtables.end();
         ++it)
    {
        it->initialise();
    }
}

/** Initialise the table for the molecule with number 'molnum' */
void FieldTable::initialiseTable(MolNum molnum)
{
    int idx = molnum_to_idx.value(molnum, -1);

    if (idx == -1)
        assertContainsTableFor(molnum);

    moltables_by_idx[idx].initialise();
}

/** Initialise the table for the grid 'grid' */
void FieldTable::initialiseTable(const Grid &grid)
{
    assertContainsTableFor(grid);

    for (QVector<GridFieldTable>::iterator it = gridtables.begin();
         it != gridtables.end();
         ++it)
    {
        if (it->grid().equals(grid))
        {
            it->initialise();
            return;
        }
    }
}

/** Return the field table for the passed grid

    \throw SireError::unavailable_resource
*/
GridFieldTable& FieldTable::getTable(const Grid &grid)
{
    this->assertContainsTableFor(grid);

    for (QVector<GridFieldTable>::iterator it = gridtables.begin();
         it != gridtables.end();
         ++it)
    {
        if (it->grid().equals(grid))
            return *it;
    }

    BOOST_ASSERT( false );

    //this line needed to remove warning about lack of return value
    return *( (GridFieldTable*) 0 );
}

/** Return the field table for the passed grid

    \throw SireError::unavailable_resource
*/
const GridFieldTable& FieldTable::getTable(const Grid &grid) const
{
    this->assertContainsTableFor(grid);

    for (QVector<GridFieldTable>::const_iterator it = gridtables.constBegin();
         it != gridtables.constEnd();
         ++it)
    {
        if (it->grid().equals(grid))
            return *it;
    }

    BOOST_ASSERT( false );

    //this line needed to remove warning about lack of return value
    return *( (GridFieldTable*) 0 );
}

/** Return the field table for the passed grid

    \throw SireError::unavailable_resource
*/
const GridFieldTable& FieldTable::constGetTable(const Grid &grid) const
{
    return this->getTable(grid);
}

/** Return whether or not this table is empty */
bool FieldTable::isEmpty() const
{
    return moltables_by_idx.isEmpty() and gridtables.isEmpty();
}

/** Assert that this contains a table for the molecule with number 'molnum'

    \throw SireError::unavailable_resource
*/
void FieldTable::assertContainsTableFor(MolNum molnum) const
{
    if (not molnum_to_idx.contains(molnum))
        throw SireError::unavailable_resource( QObject::tr(
                "This FieldTable does not contain an entry for the molecule "
                "with number %1.")
                    .arg(molnum), CODELOC );
}

/** Assert that this contains a table for the passed grid

    \throw SireError::unavailable_resource
*/
void FieldTable::assertContainsTableFor(const Grid &grid) const
{
    for (QVector<GridFieldTable>::const_iterator it = gridtables.constBegin();
         it != gridtables.constEnd();
         ++it)
    {
        if (it->grid().equals(grid))
            return;
    }

    throw SireError::unavailable_resource( QObject::tr(
            "This field table does not contain an entry for the grid %1.")
                .arg(grid.toString()), CODELOC );
}

/** Add the contents of the table 'other' onto this table. This will only
    add the fields for the molecules / grids that are in both tables */
void FieldTable::add(const FieldTable &other)
{
    for (QHash<MolNum,qint32>::const_iterator it = other.molnum_to_idx.constBegin();
         it != other.molnum_to_idx.constEnd();
         ++it)
    {
        int idx = molnum_to_idx.value(it.key(), -1);

        if (idx != -1)
            moltables_by_idx[idx] += other.moltables_by_idx[it.value()];
    }

    for (QVector<GridFieldTable>::const_iterator it = other.gridtables.constBegin();
         it != other.gridtables.constEnd();
         ++it)
    {
        for (int i=0; i<gridtables.count(); ++i)
        {
            if (gridtables.at(i).grid().equals(it->grid()))
            {
                gridtables[i] += *it;
                continue;
            }
        }
    }
}

/** Subtract the contents of the table 'other' from this table. This will only
    subtract the fields for the molecules / grids that are in both tables */
void FieldTable::subtract(const FieldTable &other)
{
    for (QHash<MolNum,qint32>::const_iterator it = other.molnum_to_idx.constBegin();
         it != other.molnum_to_idx.constEnd();
         ++it)
    {
        int idx = molnum_to_idx.value(it.key(), -1);

        if (idx != -1)
            moltables_by_idx[idx] -= other.moltables_by_idx[it.value()];
    }

    for (QVector<GridFieldTable>::const_iterator it = other.gridtables.constBegin();
         it != other.gridtables.constEnd();
         ++it)
    {
        for (int i=0; i<gridtables.count(); ++i)
        {
            if (gridtables.at(i).grid().equals(it->grid()))
            {
                gridtables[i] -= *it;
                continue;
            }
        }
    }
}

/** Add the field 'field' onto all of the atom / grid points in this table */
void FieldTable::add(const Vector &field)
{
    for (QVector<MolFieldTable>::iterator it = moltables_by_idx.begin();
         it != moltables_by_idx.end();
         ++it)
    {
        it->add(field);
    }

    for (QVector<GridFieldTable>::iterator it = gridtables.begin();
         it != gridtables.end();
         ++it)
    {
        it->add(field);
    }
}

/** Subtract the field 'field' from all of the atom / grid points in this table */
void FieldTable::subtract(const Vector &field)
{
    for (QVector<MolFieldTable>::iterator it = moltables_by_idx.begin();
         it != moltables_by_idx.end();
         ++it)
    {
        it->subtract(field);
    }

    for (QVector<GridFieldTable>::iterator it = gridtables.begin();
         it != gridtables.end();
         ++it)
    {
        it->subtract(field);
    }
}

/** Set the field at all atom and grid points equal to 'field' */
void FieldTable::setAll(const Vector &field)
{
    for (QVector<MolFieldTable>::iterator it = moltables_by_idx.begin();
         it != moltables_by_idx.end();
         ++it)
    {
        it->setAll(field);
    }

    for (QVector<GridFieldTable>::iterator it = gridtables.begin();
         it != gridtables.end();
         ++it)
    {
        it->setAll(field);
    }
}

/** Multiply the field at all atom and grid points by 'value' */
void FieldTable::multiply(double value)
{
    for (QVector<MolFieldTable>::iterator it = moltables_by_idx.begin();
         it != moltables_by_idx.end();
         ++it)
    {
        it->multiply(value);
    }

    for (QVector<GridFieldTable>::iterator it = gridtables.begin();
         it != gridtables.end();
         ++it)
    {
        it->multiply(value);
    }
}

/** Divide the field at all atom and grid points by 'value' */
void FieldTable::divide(double value)
{
    for (QVector<MolFieldTable>::iterator it = moltables_by_idx.begin();
         it != moltables_by_idx.end();
         ++it)
    {
        it->divide(value);
    }

    for (QVector<GridFieldTable>::iterator it = gridtables.begin();
         it != gridtables.end();
         ++it)
    {
        it->divide(value);
    }
}

/** Return the index of the molecule with number 'molnum' in this table

    \throw SireMol::missing_molecule
*/
int FieldTable::indexOf(MolNum molnum) const
{
    QHash<MolNum,qint32>::const_iterator it = molnum_to_idx.constFind(molnum);

    if (it == molnum_to_idx.constEnd())
        throw SireMol::missing_molecule( QObject::tr(
            "There is no molecule with number %1 in this field table.")
                .arg(molnum), CODELOC );

    return it.value();
}
