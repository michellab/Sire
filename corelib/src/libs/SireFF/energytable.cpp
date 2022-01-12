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

#include "energytable.h"

#include "SireMol/moleculeview.h"
#include "SireMol/moleculegroup.h"
#include "SireMol/atomselection.h"

#include "SireMol/mover.hpp"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireMol/errors.h"
#include "SireError/errors.h"

using namespace SireFF;
using namespace SireMol;
using namespace SireMaths;
using namespace SireBase;
using namespace SireStream;

/////////
///////// Implementation of MolEnergyTable
/////////

static const RegisterMetaType<MolEnergyTable> r_molenergytable(NO_ROOT);

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                      const MolEnergyTable &molenergytable)
{
    writeHeader(ds, r_molenergytable, 1);

    SharedDataStream sds(ds);

    sds << molenergytable.molnum
        << molenergytable.moluid
        << molenergytable.ncgroups
        << molenergytable.cgidx_to_idx
        << static_cast<const PackedArray2D<Vector>&>(molenergytable);

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds,
                                      MolEnergyTable &molenergytable)
{
    VersionID v = readHeader(ds, r_molenergytable);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> molenergytable.molnum
            >> molenergytable.moluid
            >> molenergytable.ncgroups
            >> molenergytable.cgidx_to_idx
            >> static_cast<PackedArray2D<Vector>&>(molenergytable);
    }
    else
        throw version_error(v, "1", r_molenergytable, CODELOC);

    return ds;
}

/** Constructor */
MolEnergyTable::MolEnergyTable() : PackedArray2D<Vector>(),
                                 molnum(0), ncgroups(0)
{}

/** Construct a table to hold the forces on all of the CutGroups that
    are selected in 'molview'. The forces are initialised to zero */
MolEnergyTable::MolEnergyTable(const MoleculeView &molview)
              : PackedArray2D<Vector>(),
                molnum(molview.data().number()),
                moluid(molview.data().info().UID()),
                ncgroups(molview.data().info().nCutGroups())
{
    //build arrays for each selected CutGroup
    AtomSelection selected_atoms = molview.selection();

    if (selected_atoms.selectedAllCutGroups())
    {
        QVector< QVector<Vector> > forces(ncgroups);
        QVector<Vector> *forces_array = forces.data();

        for (CGIdx i(0); i<ncgroups; ++i)
        {
            forces_array[i] = QVector<Vector>(molview.data().info().nAtoms(i),
                                              Vector(0));
        }

        PackedArray2D<Vector>::operator=(forces);
    }
    else
    {
        QVector< QVector<Vector> > forces(selected_atoms.nSelectedCutGroups());
        cgidx_to_idx.reserve(selected_atoms.nSelectedCutGroups());

        QVector<Vector> *forces_array = forces.data();
        qint32 idx = 0;

        foreach (CGIdx i, selected_atoms.selectedCutGroups())
        {
            forces_array[i] = QVector<Vector>(molview.data().info().nAtoms(i),
                                              Vector(0));

            cgidx_to_idx.insert(i, idx);
            ++idx;
        }

        PackedArray2D<Vector>::operator=(forces);
    }
}

/** Copy constructor */
MolEnergyTable::MolEnergyTable(const MolEnergyTable &other)
              : PackedArray2D<Vector>(other),
                molnum(other.molnum),
                moluid(other.moluid),
                ncgroups(other.ncgroups),
                cgidx_to_idx(other.cgidx_to_idx)
{}

/** Destructor */
MolEnergyTable::~MolEnergyTable()
{}

/** Copy assignment operator */
MolEnergyTable& MolEnergyTable::operator=(const MolEnergyTable &other)
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

/** Comparison operator */
bool MolEnergyTable::operator==(const MolEnergyTable &other) const
{
    return molnum == other.molnum and moluid == other.moluid
             and cgidx_to_idx == other.cgidx_to_idx
             and PackedArray2D<Vector>::operator==(other);
}

/** Comparison operator */
bool MolEnergyTable::operator!=(const MolEnergyTable &other) const
{
    return not this->operator==(other);
}

/** Set the force at all points in this table equal to 'force' */
MolEnergyTable& MolEnergyTable::operator=(const Vector &force)
{
    this->setAll(force);
    return *this;
}

/** Add the forces in 'other' onto this table - this only adds forces
    to atoms that are in this table

    \throw SireError::incompatible_error
*/
MolEnergyTable& MolEnergyTable::operator+=(const MolEnergyTable &other)
{
    this->add(other);
    return *this;
}

/** Subtract the forces in 'other' from this table - this only subtracts forces
    from atoms that are in this table

    \throw SireError::incompatible_error
*/
MolEnergyTable& MolEnergyTable::operator-=(const MolEnergyTable &other)
{
    this->subtract(other);
    return *this;
}

/** Return the table where the forces in 'other' have been added to the
    forces on the atoms in this table

    \throw SireError::incompatible_error
*/
MolEnergyTable MolEnergyTable::operator+(const MolEnergyTable &other) const
{
    MolEnergyTable ret(*this);
    ret += other;
    return ret;
}

/** Return the table where the forces in 'other' have been subtracted from the
    forces on the atoms in this table

    \throw SireError::incompatible_error
*/
MolEnergyTable MolEnergyTable::operator-(const MolEnergyTable &other) const
{
    MolEnergyTable ret(*this);
    ret -= other;
    return ret;
}

/** Add 'force' to all of the points in this table */
MolEnergyTable& MolEnergyTable::operator+=(const Vector &force)
{
    this->add(force);
    return *this;
}

/** Subtract 'force' from all of the points in this table */
MolEnergyTable& MolEnergyTable::operator-=(const Vector &force)
{
    this->subtract(force);
    return *this;
}

/** Return the table where 'force' has been added to all of the
    points in this table */
MolEnergyTable MolEnergyTable::operator+(const Vector &force) const
{
    MolEnergyTable ret(*this);
    ret += force;
    return ret;
}

/** Return the table where 'force' has been subtracted from all of the
    points in this table */
MolEnergyTable MolEnergyTable::operator-(const Vector &force) const
{
    MolEnergyTable ret(*this);
    ret -= force;
    return ret;
}

/** Multiply the force at each point in this table by 'value' */
MolEnergyTable& MolEnergyTable::operator*=(double value)
{
    this->multiply(value);
    return *this;
}

/** Divide the force at each point in this table by 'value' */
MolEnergyTable& MolEnergyTable::operator/=(double value)
{
    this->divide(value);
    return *this;
}

/** Return the table where the force at each point has been
    mulitiplied by 'value' */
MolEnergyTable MolEnergyTable::operator*(double value) const
{
    MolEnergyTable ret(*this);
    ret *= value;
    return ret;
}

/** Return the table where the force at each point has been
    divided by 'value' */
MolEnergyTable MolEnergyTable::operator/(double value) const
{
    MolEnergyTable ret(*this);
    ret /= value;
    return ret;
}

/** Return the negative of this table */
MolEnergyTable MolEnergyTable::operator-() const
{
    MolEnergyTable ret(*this);
    ret *= -1;
    return ret;
}

/** Initialise this table - this resets all of the forces back to zero */
void MolEnergyTable::initialise()
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

/** Return an array of all of the forces on the atoms, in CGAtomIdx order */
QVector<Vector> MolEnergyTable::toVector() const
{
    return PackedArray2D<Vector>::toQVector();
}

/** Add the force 'force' onto this table. This ignores
    forces calculated for atoms that are in CutGroups that are
    not in this table - this returns whether or not the
    atom is in this table

    \throw SireError::invalid_index
*/
bool MolEnergyTable::add(const CGAtomIdx &cgatomidx, const Vector &force)
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

    this->operator()(i, j) += force;

    return true;
}

/** Subtract the force 'force' from this table. This ignores
    forces calculated for atoms that are in CutGroups that are
    not in this table - this returns whether or not the
    atom is in this table

    \throw SireError::invalid_index
*/
bool MolEnergyTable::subtract(const CGAtomIdx &cgatomidx, const Vector &force)
{
    return this->add( cgatomidx, -force );
}

static void addEnergy(const Vector &force, Vector *forces, const int nats)
{
    for (int i=0; i<nats; ++i)
    {
        forces[i] += force;
    }
}

/** Add the force 'force' onto this table for all of the atoms
    in 'selected_atoms'. This ignores forces calculated for atoms
    that are in CutGroups that are not in this table - this returns whether
    or not any selected atoms are in this table

    \throw SireError::incompatible_error
*/
bool MolEnergyTable::add(const AtomSelection &selected_atoms, const Vector &force)
{
    this->assertCompatibleWith(selected_atoms);

    if (selected_atoms.selectedNone() or this->isEmpty())
        return false;

    bool changed_atoms = false;

    if (selected_atoms.selectedAll())
    {
        //this is easy - all atoms are selected for updating,
        //so just update all of the forces in this table
        ::addEnergy(force, this->valueData(), this->nValues());

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
                    ::addEnergy(force, this->data(i), this->nValues(i));
                    changed_atoms = true;
                }
                else
                {
                    QSet<Index> idxs = selected_atoms.selectedAtoms(i);

                    Vector *atomforces = this->data(i);
                    const int nats = this->nValues(i);

                    foreach (Index idx, idxs)
                    {
                        BOOST_ASSERT( idx >= 0 and idx < nats );
                        atomforces[idx] += force;
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
                    ::addEnergy(force, this->data(i), this->nValues(i));
                    changed_atoms = true;
                }
                else
                {
                    QSet<Index> idxs = selected_atoms.selectedAtoms(i);

                    Vector *atomforces = this->data(i);
                    const int nats = this->nValues(i);

                    foreach (Index idx, idxs)
                    {
                        BOOST_ASSERT( idx >= 0 and idx < nats );
                        atomforces[idx] += force;
                    }

                    changed_atoms = true;
                }
            }
        }
    }
    else
    {
        //harder case - not all atoms are in this forcetable
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
                    ::addEnergy(force, this->data(i), this->nValues(i));
                    changed_atoms = true;
                }
                else
                {
                    QSet<Index> idxs = selected_atoms.selectedAtoms(cgidx);

                    Vector *atomforces = this->data(i);
                    const int nats = this->nValues(i);

                    foreach (Index idx, idxs)
                    {
                        BOOST_ASSERT( idx >= 0 and idx < nats );
                        atomforces[idx] += force;
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
                    ::addEnergy(force, this->data(i), this->nValues(i));
                    changed_atoms = true;
                }
                else if (selected_atoms.selected(cgidx))
                {
                    QSet<Index> idxs = selected_atoms.selectedAtoms(cgidx);

                    Vector *atomforces = this->data(i);
                    const int nats = this->nValues(i);

                    foreach (Index idx, idxs)
                    {
                        BOOST_ASSERT( idx >= 0 and idx < nats );
                        atomforces[idx] += force;
                    }

                    changed_atoms = true;
                }
            }
        }
    }

    return changed_atoms;
}

/** Subtract the force 'force' from this table for all of the atoms
    in 'selected_atoms'. This ignores forces calculated for atoms
    that are in CutGroups that are not in this table - this returns whether
    or not any selected atoms are in this table

    \throw SireError::incompatible_error
*/
bool MolEnergyTable::subtract(const AtomSelection &selected_atoms, const Vector &force)
{
    return MolEnergyTable::add( selected_atoms, -force );
}

/** Add the forces contained in 'other' onto this force table. This will only
    add the forces for CutGroups that are in both tables */
void MolEnergyTable::add(const MolEnergyTable &other)
{
    if (this == &other)
    {
        //just double everything
        this->operator*=(2);
        return;
    }

    if (molnum != other.molnum)
        throw SireError::incompatible_error( QObject::tr(
                "You cannot combine the force table for molecule %1 with the "
                "force table for molecule %2. The molecules must be the same.")
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

/** Subtract the forces contained in 'other' from this force table. This will only
    subtract the forces for CutGroups that are in both tables */
void MolEnergyTable::subtract(const MolEnergyTable &other)
{
    if (this == &other)
    {
        this->setAll( Vector(0) );
        return;
    }

    this->add( -other );
}

/** Add the force 'force' onto all of the atom points in this table */
void MolEnergyTable::add(const Vector &force)
{
    int nvals = PackedArray2D<Vector>::nValues();

    if (nvals > 0)
    {
        Vector *vals = PackedArray2D<Vector>::valueData();

        for (int i=0; i<nvals; ++i)
        {
            vals[i] += force;
        }
    }
}

/** Subtract the force 'force' from all of the atom points in this table */
void MolEnergyTable::subtract(const Vector &force)
{
    this->add( -force );
}

/** Set all of the forces at the atom points equal to 'force' */
void MolEnergyTable::setAll(const Vector &force)
{
    int nvals = PackedArray2D<Vector>::nValues();

    if (nvals > 0)
    {
        Vector *vals = PackedArray2D<Vector>::valueData();

        for (int i=0; i<nvals; ++i)
        {
            vals[i] = force;
        }
    }
}

/** Multiply the force at all atom points by 'value' */
void MolEnergyTable::multiply(double value)
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

/** Divide the force at all atom points by 'value' */
void MolEnergyTable::divide(double value)
{
    this->multiply( 1.0 / value );
}

void MolEnergyTable::assertCompatibleWith(const AtomSelection &selection) const
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
            "This MolEnergyTable is incompatible with the passed atom selection."),
                CODELOC );
}

/** Return an array of all of the forces on the atoms selected in 'selection'

    \throw SireError::incompatible_error
*/
QVector<Vector> MolEnergyTable::toVector(const AtomSelection &selection) const
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
                const Vector *groupforces = PackedArray2D<Vector>::constData(i);

                if (selection.selectedAll(i))
                {
                    const int nats = PackedArray2D<Vector>::nValues(i);

                    quickCopy<Vector>(value, groupforces, nats);
                    value += nats;
                }
                else
                {
                    QList<Index> idxs = selection.selectedAtoms(i).toList();
                    std::sort(idxs.begin(), idxs.end());

                    foreach (Index idx, idxs)
                    {
                        *value = groupforces[idx];
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
                const Vector *groupforces = PackedArray2D<Vector>::constData(i);

                if (selection.selectedAll(i))
                {
                    const int nats = PackedArray2D<Vector>::nValues(i);

                    quickCopy<Vector>(value, groupforces, nats);
                    value += nats;
                }
                else
                {
                    QList<Index> idxs = selection.selectedAtoms(i).toList();
                    std::sort(idxs.begin(), idxs.end());

                    foreach (Index idx, idxs)
                    {
                        *value = groupforces[idx];
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

            const Vector *groupforces = PackedArray2D<Vector>::constData(i);

            if (selection.selectedAll(cgidx))
            {
                const int nats = PackedArray2D<Vector>::nValues(i);

                quickCopy<Vector>(value, groupforces, nats);
                value += nats;
            }
            else
            {
                QList<Index> idxs = selection.selectedAtoms(cgidx).toList();
                std::sort(idxs.begin(), idxs.end());

                foreach (Index idx, idxs)
                {
                    *value = groupforces[idx];
                    ++value;
                }
            }
        }
    }

    return vals;
}

const char* MolEnergyTable::typeName()
{
    return QMetaType::typeName( qMetaTypeId<MolEnergyTable>() );
}

////////
//////// Implementation of EnergyTable
////////

static const RegisterMetaType<EnergyTable> r_forcetable(NO_ROOT);

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                      const EnergyTable &forcetable)
{
    writeHeader(ds, r_forcetable, 1);

    SharedDataStream sds(ds);

    sds << forcetable.tables_by_idx << forcetable.molnum_to_idx;

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds,
                                      EnergyTable &forcetable)
{
    VersionID v = readHeader(ds, r_forcetable);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> forcetable.tables_by_idx >> forcetable.molnum_to_idx;
    }
    else
        throw version_error(v, "1", r_forcetable, CODELOC);

    return ds;
}

/** Constructor */
EnergyTable::EnergyTable()
{}

/** Construct a table that holds all of the forces on all of the atoms
    for all of the CutGroups viewed in all of the molecules in the passed
    molecule group */
EnergyTable::EnergyTable(const MoleculeGroup &molgroup)
{
    if (molgroup.isEmpty())
        return;

    int nmols = molgroup.nMolecules();

    tables_by_idx = QVector<MolEnergyTable>(nmols);
    molnum_to_idx.reserve(nmols);

    MolEnergyTable *tables_by_idx_array = tables_by_idx.data();

    qint32 i = 0;

    for (MoleculeGroup::const_iterator it = molgroup.constBegin();
         it != molgroup.constEnd();
         ++it)
    {
        tables_by_idx_array[i] = MolEnergyTable(*it);
        molnum_to_idx.insert(it->data().number(), i);
        ++i;
    }
}

/** Copy constructor */
EnergyTable::EnergyTable(const EnergyTable &other)
           : tables_by_idx(other.tables_by_idx),
             molnum_to_idx(other.molnum_to_idx)
{}

/** Destructor */
EnergyTable::~EnergyTable()
{}

/** Copy assignment operator */
EnergyTable& EnergyTable::operator=(const EnergyTable &other)
{
    tables_by_idx = other.tables_by_idx;
    molnum_to_idx = other.molnum_to_idx;

    return *this;
}

/** Set all of the forces at all of the points in this table equal to 'force' */
EnergyTable& EnergyTable::operator=(const Vector &force)
{
    this->setAll(force);
    return *this;
}

/** Comparison operator */
bool EnergyTable::operator==(const EnergyTable &other) const
{
    return tables_by_idx == other.tables_by_idx;
}

/** Comparison operator */
bool EnergyTable::operator!=(const EnergyTable &other) const
{
    return tables_by_idx != other.tables_by_idx;
}

/** Add the forces from 'other' onto this table. This only adds the forces
    for molecules / grids that are in both tables */
EnergyTable& EnergyTable::operator+=(const EnergyTable &other)
{
    this->add(other);
    return *this;
}

/** Subtract the forces from 'other' from this table. This only subtracts
    the forces for molecules / grids that are in both tables */
EnergyTable& EnergyTable::operator-=(const EnergyTable &other)
{
    this->subtract(other);
    return *this;
}

/** Return the sum of this table with 'other' - this only adds the
    forces from 'other' to this table for molecules / grids that are
    in both tables */
EnergyTable EnergyTable::operator+(const EnergyTable &other) const
{
    EnergyTable ret(*this);
    ret += other;
    return ret;
}

/** Return the difference of this table with 'other' - this only subtracts the
    forces from 'other' to this table for molecules / grids that are
    in both tables */
EnergyTable EnergyTable::operator-(const EnergyTable &other) const
{
    EnergyTable ret(*this);
    ret -= other;
    return ret;
}

/** Add the force 'force' to all of the atom and grid points in this table */
EnergyTable& EnergyTable::operator+=(const Vector &force)
{
    this->add(force);
    return *this;
}

/** Substract the force 'force' from all of the atom and grid
    points in this table */
EnergyTable& EnergyTable::operator-=(const Vector &force)
{
    this->subtract(force);
    return *this;
}

/** Return the result of adding 'force' onto all of the atom
    and grid points in this table */
EnergyTable EnergyTable::operator+(const Vector &force) const
{
    EnergyTable ret(*this);
    ret += force;
    return ret;
}

/** Return the result of subtracting 'force' from all of the atom
    and grid points in this table */
EnergyTable EnergyTable::operator-(const Vector &force) const
{
    EnergyTable ret(*this);
    ret -= force;
    return ret;
}

/** Multiply the forces at all points in this table by 'value' */
EnergyTable& EnergyTable::operator*=(double value)
{
    this->multiply(value);
    return *this;
}

/** Divide the forces at all points in this table by 'value' */
EnergyTable& EnergyTable::operator/=(double value)
{
    this->divide(value);
    return *this;
}

/** Return the result of multiplying the forces at all points by 'value' */
EnergyTable EnergyTable::operator*(double value) const
{
    EnergyTable ret(*this);
    ret *= value;
    return ret;
}

/** Return the result of dividing the forces at all points by 'value' */
EnergyTable EnergyTable::operator/(double value) const
{
    EnergyTable ret(*this);
    ret /= value;
    return ret;
}

/** Return the result of negating the force at all points */
EnergyTable EnergyTable::operator-() const
{
    EnergyTable ret(*this);

    for (QVector<MolEnergyTable>::iterator it = ret.tables_by_idx.begin();
         it != ret.tables_by_idx.end();
         ++it)
    {
        *it = -(*it);
    }

    return ret;
}

/** Add the contents of the table 'other' onto this table. This will only
    add the forces for the molecules / grids that are in both tables */
void EnergyTable::add(const EnergyTable &other)
{
    for (QHash<MolNum,qint32>::const_iterator it = other.molnum_to_idx.constBegin();
         it != other.molnum_to_idx.constEnd();
         ++it)
    {
        int idx = molnum_to_idx.value(it.key(), -1);

        if (idx != -1)
            tables_by_idx[idx] += other.tables_by_idx[it.value()];
    }
}

/** Subtract the contents of the table 'other' from this table. This will only
    subtract the forces for the molecules / grids that are in both tables */
void EnergyTable::subtract(const EnergyTable &other)
{
    for (QHash<MolNum,qint32>::const_iterator it = other.molnum_to_idx.constBegin();
         it != other.molnum_to_idx.constEnd();
         ++it)
    {
        int idx = molnum_to_idx.value(it.key(), -1);

        if (idx != -1)
            tables_by_idx[idx] -= other.tables_by_idx[it.value()];
    }
}

/** Add the force 'force' onto all of the atom / grid points in this table */
void EnergyTable::add(const Vector &force)
{
    for (QVector<MolEnergyTable>::iterator it = tables_by_idx.begin();
         it != tables_by_idx.end();
         ++it)
    {
        it->add(force);
    }
}

/** Subtract the force 'force' from all of the atom / grid points in this table */
void EnergyTable::subtract(const Vector &force)
{
    for (QVector<MolEnergyTable>::iterator it = tables_by_idx.begin();
         it != tables_by_idx.end();
         ++it)
    {
        it->subtract(force);
    }
}

/** Set the force at all atom and grid points equal to 'force' */
void EnergyTable::setAll(const Vector &force)
{
    for (QVector<MolEnergyTable>::iterator it = tables_by_idx.begin();
         it != tables_by_idx.end();
         ++it)
    {
        it->setAll(force);
    }
}

/** Multiply the force at all atom and grid points by 'value' */
void EnergyTable::multiply(double value)
{
    for (QVector<MolEnergyTable>::iterator it = tables_by_idx.begin();
         it != tables_by_idx.end();
         ++it)
    {
        it->multiply(value);
    }
}

/** Divide the force at all atom and grid points by 'value' */
void EnergyTable::divide(double value)
{
    for (QVector<MolEnergyTable>::iterator it = tables_by_idx.begin();
         it != tables_by_idx.end();
         ++it)
    {
        it->divide(value);
    }
}

/** Initialise all of the tables - this resets all of the forces
    back to zero */
void EnergyTable::initialiseTables()
{
    int nmols = tables_by_idx.count();

    //    qDebug() << " In initialiseTables nmols is " << nmols;

    if (nmols > 0)
    {
        MolEnergyTable *tables_by_idx_array = tables_by_idx.data();

        for (int i=0; i<nmols; ++i)
        {
            tables_by_idx_array[i].initialise();
        }
    }
}

/** Initialise all of the forces for the table for the molecule
    with number 'molnum'

    \throw SireMol::missing_molecule
*/
void EnergyTable::initialiseTable(MolNum molnum)
{
    this->getTable(molnum).initialise();
}

/** Return the index of the molecule with number 'molnum' in this table

    \throw SireMol::missing_molecule
*/
int EnergyTable::indexOf(MolNum molnum) const
{
    QHash<MolNum,qint32>::const_iterator it = molnum_to_idx.constFind(molnum);

    if (it == molnum_to_idx.constEnd())
        throw SireMol::missing_molecule( QObject::tr(
            "There is no molecule with number %1 in this forcetable.")
                .arg(molnum), CODELOC );

    return it.value();
}

/** Assert that this forcetable contains a table for the
    forces for the molecule at number 'molnum'

    \throw SireMol::missing_molecule
*/
void EnergyTable::assertContainsTableFor(MolNum molnum) const
{
    QList<MolNum> molnums = molnum_to_idx.keys();

    qDebug() << " THE MOLNUMS ARE ";

    for (int i=0; i<molnums.length() ; i++)
      qDebug() << " molnum " << molnums[i].toString() ;


    if (not this->containsTable(molnum))
        throw SireMol::missing_molecule( QObject::tr(
            "This energy table does not contain a table for the "
            "molecule with number %1.")
                .arg(molnum), CODELOC );
}

const char* EnergyTable::typeName()
{
    return QMetaType::typeName( qMetaTypeId<EnergyTable>() );
}
