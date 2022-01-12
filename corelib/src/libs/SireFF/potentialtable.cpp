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

#include "potentialtable.h"

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
///////// Implementation of MolPotentialTable
/////////

static const RegisterMetaType<MolPotentialTable> r_moltable(NO_ROOT);

QDataStream &operator<<(QDataStream &ds,
                                      const MolPotentialTable &moltable)
{
    writeHeader(ds, r_moltable, 1);

    SharedDataStream sds(ds);

    sds << moltable.molnum << moltable.moluid << moltable.ncgroups
        << moltable.cgidx_to_idx
        << static_cast<const PackedArray2D<MolarEnergy>&>(moltable);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, MolPotentialTable &moltable)
{
    VersionID v = readHeader(ds, r_moltable);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> moltable.molnum >> moltable.moluid >> moltable.ncgroups
            >> moltable.cgidx_to_idx
            >> static_cast<PackedArray2D<MolarEnergy>&>(moltable);
    }
    else
        throw version_error( v, "1", r_moltable, CODELOC );

    return ds;
}

/** Null constructor */
MolPotentialTable::MolPotentialTable()
                  : PackedArray2D<MolarEnergy>(), molnum(0), ncgroups(0)
{}

/** Construct to hold the potential acting at the points of all of the atoms
    of all of the cutgroups viewed in 'molview' */
MolPotentialTable::MolPotentialTable(const MoleculeView &molview)
              : PackedArray2D<MolarEnergy>(),
                molnum(molview.data().number()),
                moluid(molview.data().info().UID()),
                ncgroups(molview.data().info().nCutGroups())
{
    //build arrays for each selected CutGroup
    AtomSelection selected_atoms = molview.selection();

    if (selected_atoms.selectedAllCutGroups())
    {
        QVector< QVector<MolarEnergy> > potentials(ncgroups);
        QVector<MolarEnergy> *potentials_array = potentials.data();

        for (CGIdx i(0); i<ncgroups; ++i)
        {
            potentials_array[i] = QVector<MolarEnergy>(molview.data().info().nAtoms(i),
                                                       MolarEnergy(0));
        }

        PackedArray2D<MolarEnergy>::operator=(potentials);
    }
    else
    {
        QVector< QVector<MolarEnergy> > potentials(selected_atoms.nSelectedCutGroups());
        cgidx_to_idx.reserve(selected_atoms.nSelectedCutGroups());

        QVector<MolarEnergy> *potentials_array = potentials.data();
        qint32 idx = 0;

        foreach (CGIdx i, selected_atoms.selectedCutGroups())
        {
            potentials_array[i] = QVector<MolarEnergy>(molview.data().info().nAtoms(i),
                                                       MolarEnergy(0));

            cgidx_to_idx.insert(i, idx);
            ++idx;
        }

        PackedArray2D<MolarEnergy>::operator=(potentials);
    }
}

/** Copy constructor */
MolPotentialTable::MolPotentialTable(const MolPotentialTable &other)
                  : PackedArray2D<MolarEnergy>(other),
                    molnum(other.molnum), moluid(other.moluid),
                    ncgroups(other.ncgroups), cgidx_to_idx(other.cgidx_to_idx)
{}

/** Destructor */
MolPotentialTable::~MolPotentialTable()
{}

/** Copy assignment operator */
MolPotentialTable& MolPotentialTable::operator=(const MolPotentialTable &other)
{
    if (this != &other)
    {
        PackedArray2D<MolarEnergy>::operator=(other);
        molnum = other.molnum;
        moluid = other.moluid;
        ncgroups = other.ncgroups;
        cgidx_to_idx = other.cgidx_to_idx;
    }

    return *this;
}

/** Set the potential at all points equal to 'potential' */
MolPotentialTable& MolPotentialTable::operator=(const MolarEnergy &potential)
{
    this->setAll(potential);
    return *this;
}

/** Comparison operator */
bool MolPotentialTable::operator==(const MolPotentialTable &other) const
{
    return this == &other or
           (molnum == other.molnum and moluid == other.moluid and
            ncgroups == other.ncgroups and cgidx_to_idx == other.cgidx_to_idx and
            PackedArray2D<MolarEnergy>::operator==(other));
}

/** Comparison operator */
bool MolPotentialTable::operator!=(const MolPotentialTable &other) const
{
    return not MolPotentialTable::operator==(other);
}

/** Add the potentials in the passed table onto this table */
MolPotentialTable& MolPotentialTable::operator+=(const MolPotentialTable &other)
{
    this->add(other);
    return *this;
}

/** Subtract the potentials in the passed table from this table */
MolPotentialTable& MolPotentialTable::operator-=(const MolPotentialTable &other)
{
    this->subtract(other);
    return *this;
}

/** Return the sum of the potentials of this table with 'other' */
MolPotentialTable MolPotentialTable::operator+(const MolPotentialTable &other) const
{
    MolPotentialTable ret(*this);
    ret += other;
    return ret;
}

/** Return the difference of the potentials of this table with 'other' */
MolPotentialTable MolPotentialTable::operator-(const MolPotentialTable &other) const
{
    MolPotentialTable ret(*this);
    ret -= other;
    return *this;
}

/** Add 'potential' to all of the points in this table */
MolPotentialTable& MolPotentialTable::operator+=(const MolarEnergy &potential)
{
    this->add(potential);
    return *this;
}

/** Subtract 'potential' from all of the points in this table */
MolPotentialTable& MolPotentialTable::operator-=(const MolarEnergy &potential)
{
    this->subtract(potential);
    return *this;
}

/** Return a table that has 'potential' added to all of the points */
MolPotentialTable MolPotentialTable::operator+(const MolarEnergy &potential) const
{
    MolPotentialTable ret(*this);
    ret += potential;
    return ret;
}

MolPotentialTable operator+(const MolarEnergy &potential,
                                          const MolPotentialTable &table)
{
    return table + potential;
}

/** Return a table that has 'potential' subtracted from all of the points */
MolPotentialTable MolPotentialTable::operator-(const MolarEnergy &potential) const
{
    MolPotentialTable ret(*this);
    ret -= potential;
    return ret;
}

/** Multiply the potential at all points by the scalar 'value' */
MolPotentialTable& MolPotentialTable::operator*=(double value)
{
    this->multiply(value);
    return *this;
}

/** Divide the potential at all points by the scalar 'value' */
MolPotentialTable& MolPotentialTable::operator/=(double value)
{
    this->divide(value);
    return *this;
}

/** Return the potential that has been multiplied by 'value' */
MolPotentialTable MolPotentialTable::operator*(double value) const
{
    MolPotentialTable ret(*this);
    ret *= value;
    return ret;
}

/** Return the potential that has been multiplied by 'value' */
MolPotentialTable operator*(double value, const MolPotentialTable &table)
{
    return table * value;
}

/** Return the potential that has been divided by 'value' */
MolPotentialTable MolPotentialTable::operator/(double value) const
{
    MolPotentialTable ret(*this);
    ret /= value;
    return ret;
}

/** Return the negative of this table */
MolPotentialTable MolPotentialTable::operator-() const
{
    MolPotentialTable ret(*this);
    ret *= -1;
    return ret;
}

const char* MolPotentialTable::typeName()
{
    return QMetaType::typeName( qMetaTypeId<MolPotentialTable>() );
}

/** Initialise this table - this clears all of the potentials, resetting them to zero */
void MolPotentialTable::initialise()
{
    int nvals = PackedArray2D<MolarEnergy>::nValues();

    if (nvals > 0)
    {
        MolarEnergy *vals = PackedArray2D<MolarEnergy>::valueData();

        for (int i=0; i<nvals; ++i)
        {
            vals[i] = MolarEnergy(0);
        }
    }
}

/** Return all of the potentials in this table in a single array */
QVector<MolarEnergy> MolPotentialTable::toVector() const
{
    return PackedArray2D<MolarEnergy>::toQVector();
}


/** Return an array of all of the potentials at the location of
    the atoms selected in 'selection'

    \throw SireError::incompatible_error
*/
QVector<MolarEnergy> MolPotentialTable::toVector(const AtomSelection &selection) const
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

    QVector<MolarEnergy> vals( selection.nSelected() );
    MolarEnergy *value = vals.data();

    if (this->selectedAll())
    {
        if (selection.selectedAllCutGroups())
        {
            const int ncg = selection.nCutGroups();

            for (CGIdx i(0); i<ncg; ++i)
            {
                const MolarEnergy *grouppotentials
                                = PackedArray2D<MolarEnergy>::constData(i);

                if (selection.selectedAll(i))
                {
                    const int nats = PackedArray2D<MolarEnergy>::nValues(i);

                    quickCopy<MolarEnergy>(value, grouppotentials, nats);
                    value += nats;
                }
                else
                {
                    QList<Index> idxs = selection.selectedAtoms(i).toList();
                    std::sort(idxs.begin(), idxs.end());

                    foreach (Index idx, idxs)
                    {
                        *value = grouppotentials[idx];
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
                const MolarEnergy *grouppotentials = PackedArray2D<MolarEnergy>::constData(i);

                if (selection.selectedAll(i))
                {
                    const int nats = PackedArray2D<MolarEnergy>::nValues(i);

                    quickCopy<MolarEnergy>(value, grouppotentials, nats);
                    value += nats;
                }
                else
                {
                    QList<Index> idxs = selection.selectedAtoms(i).toList();
                    std::sort(idxs.begin(), idxs.end());

                    foreach (Index idx, idxs)
                    {
                        *value = grouppotentials[idx];
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
                "Cannot return the potentials as while all CutGroups are selected, "
                "not all CutGroups are present in the forcetable."), CODELOC );

        QList<CGIdx> cgidxs = selection.selectedCutGroups();
        std::sort(cgidxs.begin(), cgidxs.end());

        foreach (CGIdx cgidx, cgidxs)
        {
            int i = cgidx_to_idx.value(cgidx, -1);

            if (i == -1)
                throw SireMol::missing_atom( QObject::tr(
                    "Cannot return the potentials as while atoms in CutGroup %1 "
                    "are selected, this CutGroup is not present in the forcetable.")
                        .arg(cgidx), CODELOC );

            const MolarEnergy *grouppotentials = PackedArray2D<MolarEnergy>::constData(i);

            if (selection.selectedAll(cgidx))
            {
                const int nats = PackedArray2D<MolarEnergy>::nValues(i);

                quickCopy<MolarEnergy>(value, grouppotentials, nats);
                value += nats;
            }
            else
            {
                QList<Index> idxs = selection.selectedAtoms(cgidx).toList();
                std::sort(idxs.begin(), idxs.end());

                foreach (Index idx, idxs)
                {
                    *value = grouppotentials[idx];
                    ++value;
                }
            }
        }
    }

    return vals;
}

/** Add the potential 'potential' onto this table - this returns whether or not the
    atom is in this table

    \throw SireError::invalid_index
*/
bool MolPotentialTable::add(const CGAtomIdx &cgatomidx, const MolarEnergy &potential)
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

    this->operator()(i, j) += potential;

    return true;
}

/** Subtract the potential 'potential' from this table - this returns whether or not the
    atom is in this table

    \throw SireError::invalid_index
*/
bool MolPotentialTable::subtract(const CGAtomIdx &cgatomidx, const MolarEnergy &potential)
{
    return this->add( cgatomidx, -potential );
}

static void addField(const MolarEnergy &potential, MolarEnergy *potentials,
                     const int nats)
{
    for (int i=0; i<nats; ++i)
    {
        potentials[i] += potential;
    }
}

/** Add the potential 'potential' onto this table for all of the atoms
    in 'selected_atoms - this returns whether
    or not any selected atoms are in this table

    \throw SireError::incompatible_error
*/
bool MolPotentialTable::add(const AtomSelection &selected_atoms,
                            const MolarEnergy &potential)
{
    this->assertCompatibleWith(selected_atoms);

    if (selected_atoms.selectedNone() or this->isEmpty())
        return false;

    bool changed_atoms = false;

    if (selected_atoms.selectedAll())
    {
        //this is easy - all atoms are selected for updating,
        //so just update all of the forces in this table
        ::addField(potential, this->valueData(), this->nValues());

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
                    ::addField(potential, this->data(i), this->nValues(i));
                    changed_atoms = true;
                }
                else
                {
                    QSet<Index> idxs = selected_atoms.selectedAtoms(i);

                    MolarEnergy *atompotentials = this->data(i);
                    const int nats = this->nValues(i);

                    foreach (Index idx, idxs)
                    {
                        BOOST_ASSERT( idx >= 0 and idx < nats );
                        atompotentials[idx] += potential;
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
                    ::addField(potential, this->data(i), this->nValues(i));
                    changed_atoms = true;
                }
                else
                {
                    QSet<Index> idxs = selected_atoms.selectedAtoms(i);

                    MolarEnergy *atompotentials = this->data(i);
                    const int nats = this->nValues(i);

                    foreach (Index idx, idxs)
                    {
                        BOOST_ASSERT( idx >= 0 and idx < nats );
                        atompotentials[idx] += potential;
                    }

                    changed_atoms = true;
                }
            }
        }
    }
    else
    {
        //harder case - not all atoms are in this potentialtable
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
                    ::addField(potential, this->data(i), this->nValues(i));
                    changed_atoms = true;
                }
                else
                {
                    QSet<Index> idxs = selected_atoms.selectedAtoms(cgidx);

                    MolarEnergy *atompotentials = this->data(i);
                    const int nats = this->nValues(i);

                    foreach (Index idx, idxs)
                    {
                        BOOST_ASSERT( idx >= 0 and idx < nats );
                        atompotentials[idx] += potential;
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
                    ::addField(potential, this->data(i), this->nValues(i));
                    changed_atoms = true;
                }
                else if (selected_atoms.selected(cgidx))
                {
                    QSet<Index> idxs = selected_atoms.selectedAtoms(cgidx);

                    MolarEnergy *atompotentials = this->data(i);
                    const int nats = this->nValues(i);

                    foreach (Index idx, idxs)
                    {
                        BOOST_ASSERT( idx >= 0 and idx < nats );
                        atompotentials[idx] += potential;
                    }

                    changed_atoms = true;
                }
            }
        }
    }

    return changed_atoms;
}

/** Subtract the potential 'potential' from this table for all of the atoms
    in 'selected_atoms' - this returns whether
    or not any selected atoms are in this table

    \throw SireError::incompatible_error
*/
bool MolPotentialTable::subtract(const AtomSelection &selected_atoms,
                                 const MolarEnergy &potential)
{
    return MolPotentialTable::add( selected_atoms, -potential );
}

/** Add the potential 'potential' onto all of the atom points in this table */
void MolPotentialTable::add(const MolarEnergy &potential)
{
    int nvals = PackedArray2D<MolarEnergy>::nValues();

    if (nvals > 0)
    {
        MolarEnergy *vals = PackedArray2D<MolarEnergy>::valueData();

        for (int i=0; i<nvals; ++i)
        {
            vals[i] += potential;
        }
    }
}

/** Subtract the potential 'potential' from all of the atom points in this table */
void MolPotentialTable::subtract(const MolarEnergy &potential)
{
    this->add( -potential );
}

/** Multiply the potential at all atom points by 'value' */
void MolPotentialTable::multiply(double value)
{
    int nvals = PackedArray2D<MolarEnergy>::nValues();

    if (nvals > 0)
    {
        MolarEnergy *vals = PackedArray2D<MolarEnergy>::valueData();

        for (int i=0; i<nvals; ++i)
        {
            vals[i] *= value;
        }
    }
}

/** Divide the potential at all atom points by 'value' */
void MolPotentialTable::divide(double value)
{
    this->multiply( 1.0 / value );
}

/** Set all of the potentials at the atom points equal to 'potential' */
void MolPotentialTable::setAll(const MolarEnergy &potential)
{
    int nvals = PackedArray2D<MolarEnergy>::nValues();

    if (nvals > 0)
    {
        MolarEnergy *vals = PackedArray2D<MolarEnergy>::valueData();

        for (int i=0; i<nvals; ++i)
        {
            vals[i] = potential;
        }
    }
}

void MolPotentialTable::assertCompatibleWith(const AtomSelection &selection) const
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
            if (selection.nSelected(i) != PackedArray2D<MolarEnergy>::nValues(i))
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
            if (selection.nSelected(it.key()) != PackedArray2D<MolarEnergy>::nValues(it.key()))
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

/** Add the potentials contained in 'other' onto this potential table. This will only
    add the potentials for CutGroups that are in both tables */
void MolPotentialTable::add(const MolPotentialTable &other)
{
    if (this == &other)
    {
        //just double everything
        this->operator*=(2);
        return;
    }

    if (molnum != other.molnum)
        throw SireError::incompatible_error( QObject::tr(
                "You cannot combine the potential table for molecule %1 with the "
                "potential table for molecule %2. The molecules must be the same.")
                    .arg(molnum).arg(other.molnum), CODELOC );

    if (moluid != other.moluid)
        throw SireError::incompatible_error( QObject::tr(
                "You cannot combine together the tables for molecule %1 as the "
                "layout UIDs are different (%2 vs. %3). They must be the same.")
                    .arg(molnum).arg(moluid.toString(), other.moluid.toString()),
                        CODELOC );

    if (this->selectedAll() and other.selectedAll())
    {
        int nvals = PackedArray2D<MolarEnergy>::nValues();

        BOOST_ASSERT( nvals == other.nValues() );

        if (nvals > 0)
        {
            MolarEnergy *vals = PackedArray2D<MolarEnergy>::valueData();
            const MolarEnergy *other_vals = other.constValueData();

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

                MolarEnergy *vals = PackedArray2D<MolarEnergy>::data(i);
                const MolarEnergy *other_vals = other.constData(idx);

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

                MolarEnergy *vals = PackedArray2D<MolarEnergy>::data(it.key());
                const MolarEnergy *other_vals = other.constData(idx);

                for (int j=0; j<nvals; ++j)
                {
                    vals[j] += other_vals[j];
                }
            }
        }
    }
}

/** Subtract the potentials contained in 'other' from this potential table. This will only
    subtract the potentials for CutGroups that are in both tables */
void MolPotentialTable::subtract(const MolPotentialTable &other)
{
    if (this == &other)
    {
        this->setAll( MolarEnergy(0) );
        return;
    }

    return this->add( -other );
}

/////////
///////// Implementation of GridPotentialTable
/////////

static const RegisterMetaType<GridPotentialTable> r_gridtable(NO_ROOT);

QDataStream &operator<<(QDataStream &ds, const GridPotentialTable &gridtable)
{
    writeHeader(ds, r_gridtable, 1);

    SharedDataStream sds(ds);

    sds << gridtable.grd << gridtable.potentialvals;

    return ds;
}

QDataStream &operator>>(QDataStream &ds, GridPotentialTable &gridtable)
{
    VersionID v = readHeader(ds, r_gridtable);

    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> gridtable.grd >> gridtable.potentialvals;
    }
    else
        throw version_error(v, "1", r_gridtable, CODELOC);

    return ds;
}

/** Null constructor */
GridPotentialTable::GridPotentialTable()
{}

/** Construct to hold the potential at each of the points of the passed grid */
GridPotentialTable::GridPotentialTable(const Grid &grid) : grd(grid)
{
    if (grd.read().nPoints() > 0)
    {
        potentialvals = QVector<MolarEnergy>(grd.read().nPoints(), MolarEnergy(0));
        potentialvals.squeeze();
    }
}

/** Copy constructor */
GridPotentialTable::GridPotentialTable(const GridPotentialTable &other)
               : grd(other.grd), potentialvals(other.potentialvals)
{}

/** Destructor */
GridPotentialTable::~GridPotentialTable()
{}

/** Copy assignment operator */
GridPotentialTable& GridPotentialTable::operator=(const GridPotentialTable &other)
{
    if (this != &other)
    {
        grd = other.grd;
        potentialvals = other.potentialvals;
    }

    return *this;
}

/** Set the potential at each of the grid points equal to 'potential' */
GridPotentialTable& GridPotentialTable::operator=(const MolarEnergy &potential)
{
    this->setAll(potential);
    return *this;
}

/** Comparison operator */
bool GridPotentialTable::operator==(const GridPotentialTable &other) const
{
    return this == &other or
           (grd == other.grd and potentialvals == other.potentialvals);
}

/** Comparison operator */
bool GridPotentialTable::operator!=(const GridPotentialTable &other) const
{
    return not this->operator==(other);
}

/** Add the potentials on the passed table onto this table. Note that this
    only adds the potentials if both tables use the same grid points */
GridPotentialTable& GridPotentialTable::operator+=(const GridPotentialTable &other)
{
    this->add(other);
    return *this;
}

/** Subtract the potentials of 'other' from this table. Note that this
    only subtracts the potentials if both tables use the same grid points */
GridPotentialTable& GridPotentialTable::operator-=(const GridPotentialTable &other)
{
    this->subtract(other);
    return *this;
}

/** Return the sum of the potentials of this table and 'other' - note that
    this returns the first table if the two tables use different grids */
GridPotentialTable GridPotentialTable::operator+(const GridPotentialTable &other) const
{
    GridPotentialTable ret(*this);
    ret += other;
    return ret;
}

/** Return the difference of the potentials of this table and 'other' - note that
    this returns the first table if the two tables use different grids */
GridPotentialTable GridPotentialTable::operator-(const GridPotentialTable &other) const
{
    GridPotentialTable ret(*this);
    ret -= other;
    return ret;
}

/** Add the potential 'potential' to all of the grid points */
GridPotentialTable& GridPotentialTable::operator+=(const MolarEnergy &potential)
{
    this->add(potential);
    return *this;
}

/** Subtract the potential 'potential' from all of the grid points */
GridPotentialTable& GridPotentialTable::operator-=(const MolarEnergy &potential)
{
    this->subtract(potential);
    return *this;
}

/** Return the table where 'potential' has been added to all of the grid points */
GridPotentialTable GridPotentialTable::operator+(const MolarEnergy &potential) const
{
    GridPotentialTable ret(*this);
    ret += potential;
    return ret;
}

/** Return the table where 'potential' has been subtracted from all of the grid points */
GridPotentialTable GridPotentialTable::operator-(const MolarEnergy &potential) const
{
    GridPotentialTable ret(*this);
    ret -= potential;
    return ret;
}

/** Multiply the potential at all grid points by 'value' */
GridPotentialTable& GridPotentialTable::operator*=(double value)
{
    this->multiply(value);
    return *this;
}

/** Divide the potential at all grid points by 'value' */
GridPotentialTable& GridPotentialTable::operator/=(double value)
{
    this->divide(value);
    return *this;
}

/** Return the table where the potential at all grid points has been
    multiplied by 'value' */
GridPotentialTable GridPotentialTable::operator*(double value) const
{
    GridPotentialTable ret(*this);
    ret *= value;
    return ret;
}

/** Return the table where the potential at all grid points has been
    divided by 'value' */
GridPotentialTable GridPotentialTable::operator/(double value) const
{
    GridPotentialTable ret(*this);
    ret /= value;
    return ret;
}

/** Return the table where the potential at all grid points is negated */
GridPotentialTable GridPotentialTable::operator-() const
{
    GridPotentialTable ret(*this);

    for (QVector<MolarEnergy>::iterator it = ret.potentialvals.begin();
         it != ret.potentialvals.end();
         ++it)
    {
        *it = -(*it);
    }

    return ret;
}

/** Return a modifiable reference to the ith grid point's potential value

    \throw SireError::invalid_index
*/
MolarEnergy& GridPotentialTable::operator[](int i)
{
    return potentialvals[ Index(i).map(potentialvals.count()) ];
}

/** Return the potential value of the ith grid point

    \throw SireError::invalid_index
*/
const MolarEnergy& GridPotentialTable::operator[](int i) const
{
    return potentialvals.at( Index(i).map(potentialvals.count()) );
}

/** Return the potential value of the ith grid point

    \throw SireError::invalid_index
*/
const MolarEnergy& GridPotentialTable::at(int i) const
{
    return this->operator[](i);
}

const char* GridPotentialTable::typeName()
{
    return QMetaType::typeName( qMetaTypeId<GridPotentialTable>() );
}

/** Initialise the potential at each grid point to equal 0 */
void GridPotentialTable::initialise()
{
    for (QVector<MolarEnergy>::iterator it = potentialvals.begin();
         it != potentialvals.end();
         ++it)
    {
        *it = MolarEnergy(0);
    }
}

/** Return the number of grid points (and thus potential values) */
int GridPotentialTable::nPoints() const
{
    return potentialvals.count();
}

/** Return the number of grid points (and thus potential values) */
int GridPotentialTable::count() const
{
    return potentialvals.count();
}

/** Return the grid that contains the points at which the potential is
    evaluated - the order of points in the grid is the same as the order
    of potential values in this table */
const Grid& GridPotentialTable::grid() const
{
    if (grd.constData() == 0)
        return Grid::null();
    else
        return grd.read();
}

/** Return the array of potential values - the order is the same
    as the order of points in the grid */
QVector<MolarEnergy> GridPotentialTable::toVector() const
{
    return potentialvals;
}

/** Add the potential 'potential' onto the potential for the ipoint'th grid point

    \throw SireError::invalid_index
*/
void GridPotentialTable::add(int ipoint, const MolarEnergy &potential)
{
    potentialvals[ Index(ipoint).map(potentialvals.count()) ] += potential;
}

/** Subtract the potential 'potential' from the potential for the ipoint'th grid point

    \throw SireError::invalid_index
*/
void GridPotentialTable::subtract(int ipoint, const MolarEnergy &potential)
{
    potentialvals[ Index(ipoint).map(potentialvals.count()) ] -= potential;
}

/** Add the potential in 'other' onto that for this table - this only
    adds the potential if the two grids are identical */
void GridPotentialTable::add(const GridPotentialTable &other)
{
    if (grd == other.grd)
    {
        int nvals = potentialvals.count();
        BOOST_ASSERT( nvals == other.potentialvals.count() );

        if (nvals > 0)
        {
            MolarEnergy *data = potentialvals.data();
            const MolarEnergy *other_data = other.potentialvals.constData();

            for (int i=0; i<nvals; ++i)
            {
                data[i] += other_data[i];
            }
        }
    }
}

/** Subtract the potential in 'other' from that for this table - this only
    subtracts the potential if the two grids are identical */
void GridPotentialTable::subtract(const GridPotentialTable &other)
{
    if (grd == other.grd)
    {
        int nvals = potentialvals.count();
        BOOST_ASSERT( nvals == other.potentialvals.count() );

        if (nvals > 0)
        {
            MolarEnergy *data = potentialvals.data();
            const MolarEnergy *other_data = other.potentialvals.constData();

            for (int i=0; i<nvals; ++i)
            {
                data[i] -= other_data[i];
            }
        }
    }
}

/** Add the potential 'potential' to all of the points in this table */
void GridPotentialTable::add(const MolarEnergy &potential)
{
    for (QVector<MolarEnergy>::iterator it = potentialvals.begin();
         it != potentialvals.end();
         ++it)
    {
        *it += potential;
    }
}

/** Subtract the potential 'potential' from all of the points in this table */
void GridPotentialTable::subtract(const MolarEnergy &potential)
{
    for (QVector<MolarEnergy>::iterator it = potentialvals.begin();
         it != potentialvals.end();
         ++it)
    {
        *it -= potential;
    }
}

/** Set the potential at all of the points in this table equal to 'potential' */
void GridPotentialTable::setAll(const MolarEnergy &potential)
{
    for (QVector<MolarEnergy>::iterator it = potentialvals.begin();
         it != potentialvals.end();
         ++it)
    {
        *it = potential;
    }
}

/** Multiply the potential at all of the points in this table by 'value' */
void GridPotentialTable::multiply(double value)
{
    for (QVector<MolarEnergy>::iterator it = potentialvals.begin();
         it != potentialvals.end();
         ++it)
    {
        *it *= value;
    }
}

/** Divide the potential at all of the points in this table by 'value' */
void GridPotentialTable::divide(double value)
{
    this->multiply( 1 / value );
}

/** Return a raw pointer to the array of potential values */
MolarEnergy* GridPotentialTable::data()
{
    return potentialvals.data();
}

/** Return a raw pointer to the array of potential values */
const MolarEnergy* GridPotentialTable::data() const
{
    return potentialvals.constData();
}

/** Return a raw pointer to the array of potential values */
const MolarEnergy* GridPotentialTable::constData() const
{
    return potentialvals.constData();
}

GridPotentialTable::iterator GridPotentialTable::begin()
{
    return potentialvals.begin();
}

GridPotentialTable::iterator GridPotentialTable::end()
{
    return potentialvals.end();
}

GridPotentialTable::const_iterator GridPotentialTable::begin() const
{
    return potentialvals.constBegin();
}

GridPotentialTable::const_iterator GridPotentialTable::end() const
{
    return potentialvals.constEnd();
}

GridPotentialTable::const_iterator GridPotentialTable::constBegin() const
{
    return potentialvals.constBegin();
}

GridPotentialTable::const_iterator GridPotentialTable::constEnd() const
{
    return potentialvals.constEnd();
}

/////////
///////// Implementation of PotentialTable
/////////

static const RegisterMetaType<PotentialTable> r_potentialtable(NO_ROOT);

QDataStream &operator<<(QDataStream &ds, const PotentialTable &potentialtable)
{
    writeHeader(ds, r_potentialtable, 1);

    SharedDataStream sds(ds);

    sds << potentialtable.moltables_by_idx << potentialtable.gridtables;

    return ds;
}

QDataStream &operator>>(QDataStream &ds, PotentialTable &potentialtable)
{
    VersionID v = readHeader(ds, r_potentialtable);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> potentialtable.moltables_by_idx >> potentialtable.gridtables;

        QHash<MolNum,qint32> molnum_to_idx;
        molnum_to_idx.reserve(potentialtable.moltables_by_idx.count());

        quint32 i = 0;

        for (QVector<MolPotentialTable>::const_iterator
                                    it = potentialtable.moltables_by_idx.constBegin();
             it != potentialtable.moltables_by_idx.constEnd();
             ++it)
        {
            molnum_to_idx.insert( it->molNum(), i );
            i += 1;
        }

        potentialtable.molnum_to_idx = molnum_to_idx;
    }
    else
        throw version_error( v, "1", r_potentialtable, CODELOC );

    return ds;
}

/** Null constructor */
PotentialTable::PotentialTable()
{}

void PotentialTable::setGroup(const MoleculeGroup &molgroup)
{
    if (molgroup.isEmpty())
        return;

    int nmols = molgroup.nMolecules();

    moltables_by_idx = QVector<MolPotentialTable>(nmols);
    moltables_by_idx.squeeze();

    molnum_to_idx = QHash<MolNum,qint32>();
    molnum_to_idx.reserve(nmols);

    MolPotentialTable *moltables_by_idx_array = moltables_by_idx.data();

    quint32 i = 0;

    for (MoleculeGroup::const_iterator it = molgroup.constBegin();
         it != molgroup.constEnd();
         ++it)
    {
        moltables_by_idx_array[i] = MolPotentialTable(*it);
        molnum_to_idx.insert(it->data().number(), i);
        ++i;
    }
}

/** Construct the table to hold the potentials at the points of all
    of the atoms in the CutGroups that are viewed in the molecules
    in 'molgroup' */
PotentialTable::PotentialTable(const MoleculeGroup &molgroup)
{
    this->setGroup(molgroup);
}

/** Construct the table to hold the potentials at all of the points
    in the passed grid */
PotentialTable::PotentialTable(const Grid &grid)
{
    gridtables.append( GridPotentialTable(grid) );
    gridtables.squeeze();
}

/** Construct the table to hold the potentials at all of the points
    of all of the passed grids */
PotentialTable::PotentialTable(const QVector<GridPtr> &grids)
{
    for (QVector<GridPtr>::const_iterator it = grids.constBegin();
         it != grids.constEnd();
         ++it)
    {
        if (it->constData() != 0)
        {
            if (not this->contains(it->read()))
                gridtables.append( GridPotentialTable(it->read()) );
        }
    }

    gridtables.squeeze();
}

/** Construct the table to hold the potentials at the points of all
    of the atoms in the CutGroups that are viewed in the molecules
    in 'molgroup', and all of the grid points in the passed grid */
PotentialTable::PotentialTable(const MoleculeGroup &molgroup, const Grid &grid)
{
    this->setGroup(molgroup);
    gridtables.append( GridPotentialTable(grid) );
    gridtables.squeeze();
}

/** Construct the table to hold the potentials at the points of all
    of the atoms in the CutGroups that are viewed in the molecules
    in 'molgroup', and all of the grid points in the passed grids */
PotentialTable::PotentialTable(const MoleculeGroup &molgroup, const QVector<GridPtr> &grids)
{
    this->setGroup(molgroup);

    for (QVector<GridPtr>::const_iterator it = grids.constBegin();
         it != grids.constEnd();
         ++it)
    {
        if (it->constData() != 0)
        {
            if (not this->contains(it->read()))
                gridtables.append( GridPotentialTable(it->read()) );
        }
    }

    gridtables.squeeze();
}

/** Copy constructor */
PotentialTable::PotentialTable(const PotentialTable &other)
           : moltables_by_idx(other.moltables_by_idx),
             gridtables(other.gridtables),
             molnum_to_idx(other.molnum_to_idx)
{}

/** Destructor */
PotentialTable::~PotentialTable()
{}

const char* PotentialTable::typeName()
{
    return QMetaType::typeName( qMetaTypeId<PotentialTable>() );
}

/** Copy assignment operator */
PotentialTable& PotentialTable::operator=(const PotentialTable &other)
{
    if (this != &other)
    {
        moltables_by_idx = other.moltables_by_idx;
        gridtables = other.gridtables;
        molnum_to_idx = other.molnum_to_idx;
    }

    return *this;
}

/** Set the potential at all points in this table equal to 'potential' */
PotentialTable& PotentialTable::operator=(const MolarEnergy &potential)
{
    this->setAll(potential);
    return *this;
}

/** Comparison operator */
bool PotentialTable::operator==(const PotentialTable &other) const
{
    return this == &other or
           (moltables_by_idx == other.moltables_by_idx and
            gridtables == other.gridtables);
}

/** Comparison operator */
bool PotentialTable::operator!=(const PotentialTable &other) const
{
    return not this->operator==(other);
}

/** Add the potentials from 'other' onto this table. This only adds the potentials
    for molecules / grids that are in both tables */
PotentialTable& PotentialTable::operator+=(const PotentialTable &other)
{
    this->add(other);
    return *this;
}

/** Subtract the potentials from 'other' from this table. This only subtracts
    the potentials for molecules / grids that are in both tables */
PotentialTable& PotentialTable::operator-=(const PotentialTable &other)
{
    this->subtract(other);
    return *this;
}

/** Return the sum of this table with 'other' - this only adds the
    potentials from 'other' to this table for molecules / grids that are
    in both tables */
PotentialTable PotentialTable::operator+(const PotentialTable &other) const
{
    PotentialTable ret(*this);
    ret += other;
    return ret;
}

/** Return the difference of this table with 'other' - this only subtracts the
    potentials from 'other' to this table for molecules / grids that are
    in both tables */
PotentialTable PotentialTable::operator-(const PotentialTable &other) const
{
    PotentialTable ret(*this);
    ret -= other;
    return ret;
}

/** Add the potential 'potential' to all of the atom and grid points in this table */
PotentialTable& PotentialTable::operator+=(const MolarEnergy &potential)
{
    this->add(potential);
    return *this;
}

/** Substract the potential 'potential' from all of the atom and grid
    points in this table */
PotentialTable& PotentialTable::operator-=(const MolarEnergy &potential)
{
    this->subtract(potential);
    return *this;
}

/** Return the result of adding 'potential' onto all of the atom
    and grid points in this table */
PotentialTable PotentialTable::operator+(const MolarEnergy &potential) const
{
    PotentialTable ret(*this);
    ret += potential;
    return ret;
}

/** Return the result of subtracting 'potential' from all of the atom
    and grid points in this table */
PotentialTable PotentialTable::operator-(const MolarEnergy &potential) const
{
    PotentialTable ret(*this);
    ret -= potential;
    return ret;
}

/** Multiply the potentials at all points in this table by 'value' */
PotentialTable& PotentialTable::operator*=(double value)
{
    this->multiply(value);
    return *this;
}

/** Divide the potentials at all points in this table by 'value' */
PotentialTable& PotentialTable::operator/=(double value)
{
    this->divide(value);
    return *this;
}

/** Return the result of multiplying the potentials at all points by 'value' */
PotentialTable PotentialTable::operator*(double value) const
{
    PotentialTable ret(*this);
    ret *= value;
    return ret;
}

/** Return the result of dividing the potentials at all points by 'value' */
PotentialTable PotentialTable::operator/(double value) const
{
    PotentialTable ret(*this);
    ret /= value;
    return ret;
}

/** Return the result of negating the potential at all points */
PotentialTable PotentialTable::operator-() const
{
    PotentialTable ret(*this);

    for (QVector<MolPotentialTable>::iterator it = ret.moltables_by_idx.begin();
         it != ret.moltables_by_idx.end();
         ++it)
    {
        *it = -(*it);
    }

    for (QVector<GridPotentialTable>::iterator it = ret.gridtables.begin();
         it != ret.gridtables.end();
         ++it)
    {
        *it = -(*it);
    }

    return ret;
}

/** Return whether or not this contains a table for the passed grid */
bool PotentialTable::contains(const Grid &grid) const
{
    for (QVector<GridPotentialTable>::const_iterator it = gridtables.constBegin();
         it != gridtables.constEnd();
         ++it)
    {
        if (it->grid().equals(grid))
            return true;
    }

    return false;
}

/** Initialise all of the tables to have a zero potential */
void PotentialTable::initialiseTables()
{
    for (QVector<MolPotentialTable>::iterator it = moltables_by_idx.begin();
         it != moltables_by_idx.end();
         ++it)
    {
        it->initialise();
    }

    for (QVector<GridPotentialTable>::iterator it = gridtables.begin();
         it != gridtables.end();
         ++it)
    {
        it->initialise();
    }
}

/** Initialise the table for the molecule with number 'molnum' */
void PotentialTable::initialiseTable(MolNum molnum)
{
    int idx = molnum_to_idx.value(molnum, -1);

    if (idx == -1)
        assertContainsTableFor(molnum);

    moltables_by_idx[idx].initialise();
}

/** Initialise the table for the grid 'grid' */
void PotentialTable::initialiseTable(const Grid &grid)
{
    assertContainsTableFor(grid);

    for (QVector<GridPotentialTable>::iterator it = gridtables.begin();
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

/** Return the potential table for the passed grid

    \throw SireError::unavailable_resource
*/
GridPotentialTable& PotentialTable::getTable(const Grid &grid)
{
    this->assertContainsTableFor(grid);

    for (QVector<GridPotentialTable>::iterator it = gridtables.begin();
         it != gridtables.end();
         ++it)
    {
        if (it->grid().equals(grid))
            return *it;
    }

    BOOST_ASSERT( false );

    //this line needed to remove warning about lack of return value
    return *( (GridPotentialTable*) 0 );
}

/** Return the potential table for the passed grid

    \throw SireError::unavailable_resource
*/
const GridPotentialTable& PotentialTable::getTable(const Grid &grid) const
{
    this->assertContainsTableFor(grid);

    for (QVector<GridPotentialTable>::const_iterator it = gridtables.constBegin();
         it != gridtables.constEnd();
         ++it)
    {
        if (it->grid().equals(grid))
            return *it;
    }

    BOOST_ASSERT( false );

    //this line needed to remove warning about lack of return value
    return *( (GridPotentialTable*) 0 );
}

/** Return the potential table for the passed grid

    \throw SireError::unavailable_resource
*/
const GridPotentialTable& PotentialTable::constGetTable(const Grid &grid) const
{
    return this->getTable(grid);
}

/** Return whether or not this table is empty */
bool PotentialTable::isEmpty() const
{
    return moltables_by_idx.isEmpty() and gridtables.isEmpty();
}

/** Assert that this contains a table for the molecule with number 'molnum'

    \throw SireError::unavailable_resource
*/
void PotentialTable::assertContainsTableFor(MolNum molnum) const
{
    if (not molnum_to_idx.contains(molnum))
        throw SireError::unavailable_resource( QObject::tr(
                "This PotentialTable does not contain an entry for the molecule "
                "with number %1.")
                    .arg(molnum), CODELOC );
}

/** Assert that this contains a table for the passed grid

    \throw SireError::unavailable_resource
*/
void PotentialTable::assertContainsTableFor(const Grid &grid) const
{
    for (QVector<GridPotentialTable>::const_iterator it = gridtables.constBegin();
         it != gridtables.constEnd();
         ++it)
    {
        if (it->grid().equals(grid))
            return;
    }

    throw SireError::unavailable_resource( QObject::tr(
            "This potential table does not contain an entry for the grid %1.")
                .arg(grid.toString()), CODELOC );
}

/** Add the contents of the table 'other' onto this table. This will only
    add the potentials for the molecules / grids that are in both tables */
void PotentialTable::add(const PotentialTable &other)
{
    for (QHash<MolNum,qint32>::const_iterator it = other.molnum_to_idx.constBegin();
         it != other.molnum_to_idx.constEnd();
         ++it)
    {
        int idx = molnum_to_idx.value(it.key(), -1);

        if (idx != -1)
            moltables_by_idx[idx] += other.moltables_by_idx[it.value()];
    }

    for (QVector<GridPotentialTable>::const_iterator it = other.gridtables.constBegin();
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
    subtract the potentials for the molecules / grids that are in both tables */
void PotentialTable::subtract(const PotentialTable &other)
{
    for (QHash<MolNum,qint32>::const_iterator it = other.molnum_to_idx.constBegin();
         it != other.molnum_to_idx.constEnd();
         ++it)
    {
        int idx = molnum_to_idx.value(it.key(), -1);

        if (idx != -1)
            moltables_by_idx[idx] -= other.moltables_by_idx[it.value()];
    }

    for (QVector<GridPotentialTable>::const_iterator it = other.gridtables.constBegin();
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

/** Add the potential 'potential' onto all of the atom / grid points in this table */
void PotentialTable::add(const MolarEnergy &potential)
{
    for (QVector<MolPotentialTable>::iterator it = moltables_by_idx.begin();
         it != moltables_by_idx.end();
         ++it)
    {
        it->add(potential);
    }

    for (QVector<GridPotentialTable>::iterator it = gridtables.begin();
         it != gridtables.end();
         ++it)
    {
        it->add(potential);
    }
}

/** Subtract the potential 'potential' from all of the atom / grid points in this table */
void PotentialTable::subtract(const MolarEnergy &potential)
{
    for (QVector<MolPotentialTable>::iterator it = moltables_by_idx.begin();
         it != moltables_by_idx.end();
         ++it)
    {
        it->subtract(potential);
    }

    for (QVector<GridPotentialTable>::iterator it = gridtables.begin();
         it != gridtables.end();
         ++it)
    {
        it->subtract(potential);
    }
}

/** Set the potential at all atom and grid points equal to 'potential' */
void PotentialTable::setAll(const MolarEnergy &potential)
{
    for (QVector<MolPotentialTable>::iterator it = moltables_by_idx.begin();
         it != moltables_by_idx.end();
         ++it)
    {
        it->setAll(potential);
    }

    for (QVector<GridPotentialTable>::iterator it = gridtables.begin();
         it != gridtables.end();
         ++it)
    {
        it->setAll(potential);
    }
}

/** Multiply the potential at all atom and grid points by 'value' */
void PotentialTable::multiply(double value)
{
    for (QVector<MolPotentialTable>::iterator it = moltables_by_idx.begin();
         it != moltables_by_idx.end();
         ++it)
    {
        it->multiply(value);
    }

    for (QVector<GridPotentialTable>::iterator it = gridtables.begin();
         it != gridtables.end();
         ++it)
    {
        it->multiply(value);
    }
}

/** Divide the potential at all atom and grid points by 'value' */
void PotentialTable::divide(double value)
{
    for (QVector<MolPotentialTable>::iterator it = moltables_by_idx.begin();
         it != moltables_by_idx.end();
         ++it)
    {
        it->divide(value);
    }

    for (QVector<GridPotentialTable>::iterator it = gridtables.begin();
         it != gridtables.end();
         ++it)
    {
        it->divide(value);
    }
}

/** Return the index of the molecule with number 'molnum' in this table

    \throw SireMol::missing_molecule
*/
int PotentialTable::indexOf(MolNum molnum) const
{
    QHash<MolNum,qint32>::const_iterator it = molnum_to_idx.constFind(molnum);

    if (it == molnum_to_idx.constEnd())
        throw SireMol::missing_molecule( QObject::tr(
            "There is no molecule with number %1 in this potential table.")
                .arg(molnum), CODELOC );

    return it.value();
}
