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

#ifndef SIREFF_FIELDTABLE_H
#define SIREFF_FIELDTABLE_H

#include <QHash>
#include <QVector>
#include <QUuid>

#include "SireBase/packedarray2d.hpp"
#include "SireMaths/vector.h"
#include "SireVol/grid.h"

#include "SireMol/molnum.h"
#include "SireMol/cgidx.h"

#include "SireUnits/dimensions.h"

SIRE_BEGIN_HEADER

namespace SireFF
{
class MolFieldTable;
class GridFieldTable;
class FieldTable;
}

SIREFF_EXPORT QDataStream& operator<<(QDataStream&, const SireFF::MolFieldTable&);
SIREFF_EXPORT QDataStream& operator>>(QDataStream&, SireFF::MolFieldTable&);

SIREFF_EXPORT QDataStream& operator<<(QDataStream&, const SireFF::GridFieldTable&);
SIREFF_EXPORT QDataStream& operator>>(QDataStream&, SireFF::GridFieldTable&);

SIREFF_EXPORT QDataStream& operator<<(QDataStream&, const SireFF::FieldTable&);
SIREFF_EXPORT QDataStream& operator>>(QDataStream&, SireFF::FieldTable&);

namespace SireMol
{
class MoleculeView;
class MoleculeGroup;
}

namespace SireFF
{

using SireMol::MoleculeView;
using SireMol::CGIdx;
using SireMol::CGAtomIdx;
using SireMol::MolNum;
using SireMol::MoleculeGroup;
using SireMol::AtomSelection;

using SireMaths::Vector;

using SireVol::Grid;
using SireVol::GridPtr;

/** This class holds the field at the points of all of the atoms of 
    selected CutGroups in a molecule. The MolFieldTable is used
    to accumulate all of the fields acting on these atoms during
    a field evaluation, and also to control which fields are
    evaluated (as only the fields on atoms in selected CutGroups
    are evaluated). This allows you to provide some control over
    the calculation, e.g. only placing a few protein residues into
    the field table, thereby preventing the fields on all atoms
    in a protein from being evaluated if they aren't actually 
    necessary.
    
    @author Christopher Woods
*/
class SIREFF_EXPORT MolFieldTable : public SireBase::PackedArray2D<SireMaths::Vector>
{

friend QDataStream& ::operator<<(QDataStream&, const MolFieldTable&);
friend QDataStream& ::operator>>(QDataStream&, MolFieldTable&);

public:
    typedef SireBase::PackedArray2D<SireMaths::Vector>::Array Array;

    MolFieldTable();
    
    MolFieldTable(const MoleculeView &molview);
    
    MolFieldTable(const MolFieldTable &other);
    
    ~MolFieldTable();
    
    MolFieldTable& operator=(const MolFieldTable &other);
    MolFieldTable& operator=(const Vector &field);

    bool operator==(const MolFieldTable &other) const;
    bool operator!=(const MolFieldTable &other) const;

    MolFieldTable& operator+=(const MolFieldTable &other);
    MolFieldTable& operator-=(const MolFieldTable &other);
    
    MolFieldTable operator+(const MolFieldTable &other) const;
    MolFieldTable operator-(const MolFieldTable &other) const;

    MolFieldTable& operator+=(const Vector &field);
    MolFieldTable& operator-=(const Vector &field);
    
    MolFieldTable operator+(const Vector &field) const;
    MolFieldTable operator-(const Vector &field) const;

    MolFieldTable& operator*=(double value);
    MolFieldTable& operator/=(double value);

    MolFieldTable operator*(double value) const;
    MolFieldTable operator/(double value) const;

    MolFieldTable operator-() const;

    static const char* typeName();
    
    const char* what() const
    {
        return MolFieldTable::typeName();
    }

    void initialise();

    int nCutGroups() const;
    int nSelectedCutGroups() const;

    bool selectedAll() const;

    bool selected(CGIdx cgidx) const;

    MolNum molNum() const;

    const QUuid& molUID() const;
    
    int map(CGIdx cgidx) const;

    QVector<Vector> toVector() const;
    QVector<Vector> toVector(const AtomSelection &selection) const;
    
    bool add(const CGAtomIdx &cgatomidx, const Vector &field);
    bool subtract(const CGAtomIdx &cgatomidx, const Vector &field);

    bool add(const AtomSelection &selected_atoms, const Vector &field);
    bool subtract(const AtomSelection &selected_atoms, const Vector &field);

    void add(const MolFieldTable &other);
    void subtract(const MolFieldTable &other);

    void add(const Vector &field);
    void subtract(const Vector &field);
    
    void setAll(const Vector &value);
    
    void multiply(double value);
    void divide(double value);

private:
    void assertCompatibleWith(const AtomSelection &selection) const;

    /** The number of this molecule */
    MolNum molnum;
    
    /** The UID of the molecular layout that the molecule possessed
        when this table was last constructed */
    QUuid moluid;
    
    /** The total number of CutGroups in this molecule */
    qint32 ncgroups;
    
    /** Index mapping CGIdx to index in forcetable
        of that CutGroup. If this is empty then there
        is a one-to-one mapping (all CutGroups are present) */
    QHash<CGIdx,qint32> cgidx_to_idx;
};

/** A GridFieldTable contains the fields at point specified by a grid */
class SIREFF_EXPORT GridFieldTable
{

friend QDataStream& ::operator<<(QDataStream&, const GridFieldTable&);
friend QDataStream& ::operator>>(QDataStream&, GridFieldTable&);

public:
    typedef QVector<Vector>::const_iterator const_iterator;
    typedef QVector<Vector>::iterator iterator;

    GridFieldTable();
    GridFieldTable(const Grid &grid);
    GridFieldTable(const GridFieldTable &other);
    
    ~GridFieldTable();
    
    GridFieldTable& operator=(const GridFieldTable &other);
    GridFieldTable& operator=(const Vector &field);
    
    bool operator==(const GridFieldTable &other) const;
    bool operator!=(const GridFieldTable &other) const;
    
    GridFieldTable& operator+=(const GridFieldTable &other);
    GridFieldTable& operator-=(const GridFieldTable &other);
    
    GridFieldTable operator+(const GridFieldTable &other) const;
    GridFieldTable operator-(const GridFieldTable &other) const;
    
    GridFieldTable& operator+=(const Vector &field);
    GridFieldTable& operator-=(const Vector &field);
    
    GridFieldTable operator+(const Vector &field) const;
    GridFieldTable operator-(const Vector &field) const;
    
    GridFieldTable& operator*=(double value);
    GridFieldTable& operator/=(double value);
    
    GridFieldTable operator*(double value) const;
    GridFieldTable operator/(double value) const;
    
    GridFieldTable operator-() const;
    
    Vector& operator[](int i);
    const Vector& operator[](int i) const;
    
    const Vector& at(int i) const;
    
    static const char* typeName();
    
    void initialise();
    
    int nPoints() const;
    int count() const;
    
    const Grid& grid() const;
    
    QVector<Vector> toVector() const;
    
    void add(int ipoint, const Vector &field);
    void subtract(int ipoint, const Vector &field);
    
    void add(const GridFieldTable &other);
    void subtract(const GridFieldTable &other);
    
    void add(const Vector &field);
    void subtract(const Vector &field);
    
    void setAll(const Vector &field);
    
    void multiply(double value);
    void divide(double value);

    Vector* data();
    const Vector* data() const;

    const Vector* constData() const;

    iterator begin();
    iterator end();
    
    const_iterator begin() const;
    const_iterator end() const;
    
    const_iterator constBegin() const;
    const_iterator constEnd() const;

private:
    /** The grid that contains the points at which the field is evaluated */
    GridPtr grd;
    
    /** The potential at each of the grid points */
    QVector<Vector> fieldvals;
};

/** A FieldTable is a workspace within which all of the fields acting 
    at the points of atoms in molecules, or the points on a grid
    may be stored. A FieldTable is used as storing the fields 
    may require lots of memory, and continually 
    creating a deleting such large amouts of memory would be inefficient. 
    Also, using a FieldTable allows for fields to be accumalated directly, 
    rather than requiring intermediate storage space for the 
    individual components.

    You create a fieldtable to hold all of the fields on all of 
    the atoms of all of the molecules in a specified MoleculeGroup,
    or at all of the points of a passed Grid.
    
    The fields are held in an array that holds the fields for 
    the molecules in the same order as the molecules appear
    in the molecule group, or in an array that holds the fields
    in the same order as they appear in the grid. 
    The fieldtable also comes with an index so you can quickly 
    look up the field for a specific molecule.

    @author Christopher Woods
*/
class SIREFF_EXPORT FieldTable
{

friend QDataStream& ::operator<<(QDataStream&, const FieldTable&);
friend QDataStream& ::operator>>(QDataStream&, FieldTable&);

public:
    FieldTable();
    FieldTable(const MoleculeGroup &molgroup);

    FieldTable(const Grid &grid);
    FieldTable(const QVector<GridPtr> &grids);
    
    FieldTable(const MoleculeGroup &molgroup, const Grid &grid);
    FieldTable(const MoleculeGroup &molgroup, const QVector<GridPtr> &grids);
    
    FieldTable(const FieldTable &other);
    
    ~FieldTable();

    static const char* typeName();

    const char* what() const
    {
        return FieldTable::typeName();
    }

    FieldTable& operator=(const FieldTable &other);
    FieldTable& operator=(const Vector &field);

    bool operator==(const FieldTable &other) const;
    bool operator!=(const FieldTable &other) const;

    FieldTable& operator+=(const FieldTable &other);
    FieldTable& operator-=(const FieldTable &other);
    
    FieldTable operator+(const FieldTable &other) const;
    FieldTable operator-(const FieldTable &other) const;

    FieldTable& operator+=(const Vector &field);
    FieldTable& operator-=(const Vector &field);
    
    FieldTable operator+(const Vector &field) const;
    FieldTable operator-(const Vector &field) const;
    
    FieldTable& operator*=(double value);
    FieldTable& operator/=(double value);
    
    FieldTable operator*(double value) const;
    FieldTable operator/(double value) const;

    FieldTable operator-() const;

    bool contains(MolNum molnum) const;
    bool contains(const Grid &grid) const;

    void initialiseTables();

    void initialiseTable(MolNum molnum);
    void initialiseTable(const Grid &grid);

    MolFieldTable& getTable(MolNum molnum);
    GridFieldTable& getTable(const Grid &grid);

    const MolFieldTable& getTable(MolNum molnum) const;
    const GridFieldTable& getTable(const Grid &grid) const;

    const MolFieldTable& constGetTable(MolNum molnum) const;
    const GridFieldTable& constGetTable(const Grid &grid) const;

    bool isEmpty() const;

    const QHash<MolNum,qint32>& index() const;
    
    int indexOf(MolNum molnum) const;
    
    QList<MolNum> molNums() const;

    int nMolecules() const;

    MolFieldTable* moleculeData();
    const MolFieldTable* moleculeData() const;
    const MolFieldTable* constMoleculeData() const;

    int nGrids() const;

    GridFieldTable* gridData();
    const GridFieldTable* gridData() const;
    const GridFieldTable* constGridData() const;

    void assertContainsTableFor(MolNum molnum) const;
    void assertContainsTableFor(const Grid &grid) const;
    
    void add(const FieldTable &other);
    void subtract(const FieldTable &other);
    
    void add(const Vector &field);
    void subtract(const Vector &field);
    
    void setAll(const Vector &field);
    
    void multiply(double value);
    void divide(double value);

private:
    void setGroup(const MoleculeGroup &molgroup);

    /** All of the molecule tables */
    QVector<MolFieldTable> moltables_by_idx;

    /** All of the grid tables */
    QVector<GridFieldTable> gridtables;

    /** Index mapping from the number of the Molecule to 
        the index of its force table in the above array */
    QHash<MolNum,qint32> molnum_to_idx;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

//////
////// Inline functions for MolFieldTable
//////

/** Return the molecule number of the molecule whose fields are
    contained in this table */
inline MolNum MolFieldTable::molNum() const
{
    return molnum;
}

/** Return the UID of the molecular layout of the molecule
    whose fields are contained in this table */
inline const QUuid& MolFieldTable::molUID() const
{
    return moluid;
}

/** Return the total number of CutGroups in the molecule whose
    fields are contained in this table */
inline int MolFieldTable::nCutGroups() const
{
    return ncgroups;
}

/** Return the number of selected CutGroups in this table 
    (the number of CutGroups for which fields are held) */
inline int MolFieldTable::nSelectedCutGroups() const
{
    return SireBase::PackedArray2D<SireMaths::Vector>::count();
}

/** Return whether or not this table contains fields for all
    of the CutGroups in the molecule */
inline bool MolFieldTable::selectedAll() const
{
    return this->nCutGroups() == this->nSelectedCutGroups();
}

/** Return whether or not the CutGroup at index 'cgidx' has been
    selected 
    
    \throw SireError::invalid_index
*/
inline bool MolFieldTable::selected(CGIdx cgidx) const
{
    cgidx = CGIdx(cgidx.map(this->nCutGroups()));
    
    return this->selectedAll() or cgidx_to_idx.contains(cgidx);
}

/** Return the index of the fieldtable for the CutGroup at index
    'cgidx'. This returns -1 if there is no fieldtable for the
    specified CutGroup
    
    \throw SireError::invalid_index
*/
inline int MolFieldTable::map(CGIdx cgidx) const
{
    if (this->selectedAll())
        return cgidx.map( this->nCutGroups() );
    else
        return cgidx_to_idx.value( CGIdx(cgidx.map(this->nCutGroups())), -1 );
}

//////
////// Inline functions for FieldTable
//////

/** Return whether or not this contains a table for the 
    molecule with number 'molnum' */
inline bool FieldTable::contains(MolNum molnum) const
{
    return molnum_to_idx.contains(molnum);
}

/** Return the table for the molecule with number 'molnum'

    \throw SireMol::missing_molecule
*/
inline MolFieldTable& FieldTable::getTable(MolNum molnum)
{
    QHash<MolNum,qint32>::const_iterator it = molnum_to_idx.constFind(molnum);
    
    if (it == molnum_to_idx.constEnd())
        this->assertContainsTableFor(molnum);
        
    return moltables_by_idx.data()[ it.value() ];
}

/** Return the table for the molecule with number 'molnum'

    \throw SireMol::missing_molecule
*/
inline const MolFieldTable& FieldTable::getTable(MolNum molnum) const
{
    QHash<MolNum,qint32>::const_iterator it = molnum_to_idx.constFind(molnum);
    
    if (it == molnum_to_idx.constEnd())
        this->assertContainsTableFor(molnum);
        
    return moltables_by_idx.constData()[ it.value() ];
}

/** Return the table for the molecule with number 'molnum'

    \throw SireMol::missing_molecule
*/
inline const MolFieldTable& FieldTable::constGetTable(MolNum molnum) const
{
    return this->getTable(molnum);
}

/** Return the number of molecule tables in this object */
inline int FieldTable::nMolecules() const
{
    return moltables_by_idx.count();
}

/** Return the number of grid tables in this object */
inline int FieldTable::nGrids() const
{
    return gridtables.count();
}

/** Return the index used to find the index into the field tables array 
    for the field table for the molecule with a specified number */
inline const QHash<MolNum,qint32>& FieldTable::index() const
{
    return molnum_to_idx;
}

/** Return the numbers of molecules that have field tables in this
    table */
inline QList<MolNum> FieldTable::molNums() const
{
    return molnum_to_idx.keys();
}

/** Return a raw point to the array of field tables for each molecule */
inline MolFieldTable* FieldTable::moleculeData()
{
    return moltables_by_idx.data();
}

/** Return a raw point to the array of field tables for each molecule */
inline const MolFieldTable* FieldTable::moleculeData() const
{
    return moltables_by_idx.constData();
}

/** Return a raw point to the array of field tables for each molecule */
inline const MolFieldTable* FieldTable::constMoleculeData() const
{
    return moltables_by_idx.constData();
}

/** Return a raw pointer to the array of field tables for each grid */
inline GridFieldTable* FieldTable::gridData()
{
    return gridtables.data();
}

/** Return a raw pointer to the array of field tables for each grid */
inline const GridFieldTable* FieldTable::gridData() const
{
    return gridtables.constData();
}

/** Return a raw pointer to the array of field tables for each grid */
inline const GridFieldTable* FieldTable::constGridData() const
{
    return gridtables.constData();
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

SireFF::MolFieldTable operator+(const SireMaths::Vector &force,
                                const SireFF::MolFieldTable &table);


SireFF::MolFieldTable operator*(double value, const SireFF::MolFieldTable &table);

SireFF::FieldTable operator+(const SireMaths::Vector &force,
                             const SireFF::FieldTable &table);


SireFF::FieldTable operator*(double value, const SireFF::FieldTable &table);

Q_DECLARE_METATYPE(SireFF::FieldTable);
Q_DECLARE_METATYPE(SireFF::GridFieldTable);
Q_DECLARE_METATYPE(SireFF::MolFieldTable);

SIRE_EXPOSE_CLASS( SireFF::FieldTable )
SIRE_EXPOSE_CLASS( SireFF::GridFieldTable )
SIRE_EXPOSE_CLASS( SireFF::MolFieldTable )

SIRE_END_HEADER

#endif
