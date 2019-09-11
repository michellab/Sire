/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2006  Christopher Woods
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

#ifndef SIREFF_FORCETABLE_H
#define SIREFF_FORCETABLE_H

#include <QHash>
#include <QVector>
#include <QUuid>

#include "SireBase/packedarray2d.hpp"
#include "SireMaths/vector.h"

#include "SireMol/molnum.h"
#include "SireMol/cgidx.h"

SIRE_BEGIN_HEADER

namespace SireFF
{
class MolForceTable;
class ForceTable;
}

SIREFF_EXPORT QDataStream& operator<<(QDataStream&, const SireFF::MolForceTable&);
SIREFF_EXPORT QDataStream& operator>>(QDataStream&, SireFF::MolForceTable&);

SIREFF_EXPORT QDataStream& operator<<(QDataStream&, const SireFF::ForceTable&);
SIREFF_EXPORT QDataStream& operator>>(QDataStream&, SireFF::ForceTable&);

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

/** This class holds the forces acting on all of the atoms of 
    selected CutGroups in a molecule. The MolForceTable is used
    to accumulate all of the forces acting on these atoms during
    a force evaluation, and also to control which forces are
    evaluated (as only the forces on atoms in selected CutGroups
    are evaluated). This allows you to provide some control over
    the calculation, e.g. only placing a few protein residues into
    the force table, thereby preventing the forces on all atoms
    in a protein from being evaluated if they aren't actually 
    necessary.
    
    @author Christopher Woods
*/
class SIREFF_EXPORT MolForceTable : public SireBase::PackedArray2D<SireMaths::Vector>
{

friend SIREFF_EXPORT QDataStream& ::operator<<(QDataStream&, const MolForceTable&);
friend SIREFF_EXPORT QDataStream& ::operator>>(QDataStream&, MolForceTable&);

public:
    typedef SireBase::PackedArray2D<SireMaths::Vector>::Array Array;

    MolForceTable();
    
    MolForceTable(const MoleculeView &molview);
    
    MolForceTable(const MolForceTable &other);
    
    ~MolForceTable();
    
    MolForceTable& operator=(const MolForceTable &other);
    MolForceTable& operator=(const Vector &force);

    bool operator==(const MolForceTable &other) const;
    bool operator!=(const MolForceTable &other) const;

    MolForceTable& operator+=(const MolForceTable &other);
    MolForceTable& operator-=(const MolForceTable &other);
    
    MolForceTable operator+(const MolForceTable &other) const;
    MolForceTable operator-(const MolForceTable &other) const;

    MolForceTable& operator+=(const Vector &force);
    MolForceTable& operator-=(const Vector &force);
    
    MolForceTable operator+(const Vector &force) const;
    MolForceTable operator-(const Vector &force) const;
    
    MolForceTable& operator*=(double value);
    MolForceTable& operator/=(double value);
    
    MolForceTable operator*(double value) const;
    MolForceTable operator/(double value) const;

    MolForceTable operator-() const;

    static const char* typeName();
    
    const char* what() const
    {
        return MolForceTable::typeName();
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
    
    bool add(const CGAtomIdx &cgatomidx, const Vector &force);
    bool subtract(const CGAtomIdx &cgatomidx, const Vector &force);

    bool add(const AtomSelection &selected_atoms, const Vector &force);
    bool subtract(const AtomSelection &selected_atoms, const Vector &force);

    void add(const MolForceTable &other);
    void subtract(const MolForceTable &other);
    
    void add(const Vector &force);
    void subtract(const Vector &force);
    
    void setAll(const Vector &force);
    
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

/** A ForceTable is a workspace within which all of the forces acting 
    on the atoms of several molecules may be stored. A ForceTable is 
    used as storing the forces requires lots of memory, and continually 
    creating a deleting such large amouts of memory would be inefficient. 
    Also, using a ForceTable allows for forces to be accumalated directly, 
    rather than requiring intermediate storage space for the 
    individual components.

    You create a forcetable to hold all of the forces on all of 
    the atoms of all of the molecules in a specified MoleculeGroup.
    The forces are held in an array that holds the forces for 
    the molecules in the same order as the molecules appear
    in the molecule group. The forcetable also comes with
    an index so you can quickly look up the forces for
    a specific molecule.

    @author Christopher Woods
*/
class SIREFF_EXPORT ForceTable
{

friend SIREFF_EXPORT QDataStream& ::operator<<(QDataStream&, const ForceTable&);
friend SIREFF_EXPORT QDataStream& ::operator>>(QDataStream&, ForceTable&);

public:
    ForceTable();
    ForceTable(const MoleculeGroup &molgroup);
    
    ForceTable(const ForceTable &other);
    
    ~ForceTable();

    static const char* typeName();

    const char* what() const
    {
        return ForceTable::typeName();
    }

    ForceTable& operator=(const ForceTable &other);
    ForceTable& operator=(const Vector &force);

    bool operator==(const ForceTable &other) const;
    bool operator!=(const ForceTable &other) const;

    ForceTable& operator+=(const ForceTable &other);
    ForceTable& operator-=(const ForceTable &other);
    
    ForceTable operator+(const ForceTable &other) const;
    ForceTable operator-(const ForceTable &other) const;

    ForceTable& operator+=(const Vector &force);
    ForceTable& operator-=(const Vector &force);
    
    ForceTable operator+(const Vector &force) const;
    ForceTable operator-(const Vector &force) const;
    
    ForceTable& operator*=(double value);
    ForceTable& operator/=(double value);
    
    ForceTable operator*(double value) const;
    ForceTable operator/(double value) const;
    
    ForceTable operator-() const;
    
    bool containsTable(MolNum molnum) const;

    void initialiseTables();

    void initialiseTable(MolNum molnum);

    MolForceTable& getTable(MolNum molnum);

    void setTable(MolNum molnum, MolForceTable& table);

    const MolForceTable& getTable(MolNum molnum) const;

    const MolForceTable& constGetTable(MolNum molnum) const;

    int count() const;

    const QHash<MolNum,qint32>& index() const;
    
    int indexOf(MolNum molnum) const;
    
    QList<MolNum> molNums() const;

    MolForceTable* data();
    const MolForceTable* data() const;
    const MolForceTable* constData() const;

    void assertContainsTableFor(MolNum molnum) const;

    void add(const ForceTable &other);
    void subtract(const ForceTable &other);
    
    void add(const Vector &force);
    void subtract(const Vector &force);
    
    void setAll(const Vector &force);
    
    void multiply(double value);
    void divide(double value);

private:
    /** All of the tables */
    QVector<MolForceTable> tables_by_idx;

    /** Index mapping from the number of the Molecule to 
        the index of its force table in the above array */
    QHash<MolNum,qint32> molnum_to_idx;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

//////
////// Inline functions for MolForceTable
//////

/** Return the molecule number of the molecule whose forces are
    contained in this table */
SIRE_ALWAYS_INLINE MolNum MolForceTable::molNum() const
{
    return molnum;
}

/** Return the UID of the molecular layout of the molecule
    whose forces are contained in this table */
SIRE_ALWAYS_INLINE const QUuid& MolForceTable::molUID() const
{
    return moluid;
}

/** Return the total number of CutGroups in the molecule whose
    forces are contained in this table */
SIRE_ALWAYS_INLINE int MolForceTable::nCutGroups() const
{
    return ncgroups;
}

/** Return the number of selected CutGroups in this table 
    (the number of CutGroups for which forces are held) */
SIRE_ALWAYS_INLINE int MolForceTable::nSelectedCutGroups() const
{
    return SireBase::PackedArray2D<SireMaths::Vector>::count();
}

/** Return whether or not this table contains forces for all
    of the CutGroups in the molecule */
SIRE_ALWAYS_INLINE bool MolForceTable::selectedAll() const
{
    return this->nCutGroups() == this->nSelectedCutGroups();
}

/** Return whether or not the CutGroup at index 'cgidx' has been
    selected 
    
    \throw SireError::invalid_index
*/
SIRE_ALWAYS_INLINE bool MolForceTable::selected(CGIdx cgidx) const
{
    cgidx = CGIdx(cgidx.map(this->nCutGroups()));
    
    return this->selectedAll() or cgidx_to_idx.contains(cgidx);
}

/** Return the index of the forcetable for the CutGroup at index
    'cgidx'. This returns -1 if there is no forcetable for the
    specified CutGroup
    
    \throw SireError::invalid_index
*/
SIRE_ALWAYS_INLINE int MolForceTable::map(CGIdx cgidx) const
{
    if (this->selectedAll())
        return cgidx.map( this->nCutGroups() );
    else
        return cgidx_to_idx.value( CGIdx(cgidx.map(this->nCutGroups())), -1 );
}

//////
////// Inline functions for ForceTable
//////

/** Return whether or not this contains a table for the 
    molecule with number 'molnum' */
SIRE_ALWAYS_INLINE bool ForceTable::containsTable(MolNum molnum) const
{
    return molnum_to_idx.contains(molnum);
}

/** Return the table for the molecule with number 'molnum'

    \throw SireMol::missing_molecule
*/
SIRE_ALWAYS_INLINE MolForceTable& ForceTable::getTable(MolNum molnum)
{
    QHash<MolNum,qint32>::const_iterator it = molnum_to_idx.constFind(molnum);
    
    if (it == molnum_to_idx.constEnd())
        this->assertContainsTableFor(molnum);
        
    return tables_by_idx.data()[ it.value() ];
}


/** Return the table for the molecule with number 'molnum'

    \throw SireMol::missing_molecule
*/
SIRE_ALWAYS_INLINE const MolForceTable& ForceTable::getTable(MolNum molnum) const
{
    QHash<MolNum,qint32>::const_iterator it = molnum_to_idx.constFind(molnum);
    
    if (it == molnum_to_idx.constEnd())
        this->assertContainsTableFor(molnum);
        
    return tables_by_idx.constData()[ it.value() ];
}

/** Return the table for the molecule with number 'molnum'

    \throw SireMol::missing_molecule
*/
SIRE_ALWAYS_INLINE const MolForceTable& ForceTable::constGetTable(MolNum molnum) const
{
    return this->getTable(molnum);
}

/** Return the number of tables (molecules) in this object */
SIRE_ALWAYS_INLINE int ForceTable::count() const
{
    return tables_by_idx.count();
}

/** Return the index used to find the index into the force tables array 
    for the force table for the molecule with a specified number */
SIRE_ALWAYS_INLINE const QHash<MolNum,qint32>& ForceTable::index() const
{
    return molnum_to_idx;
}

/** Return the numbers of molecules that have force tables in this
    table */
SIRE_ALWAYS_INLINE QList<MolNum> ForceTable::molNums() const
{
    return molnum_to_idx.keys();
}

/** Return a raw point to the array of force tables for each molecule */
SIRE_ALWAYS_INLINE MolForceTable* ForceTable::data()
{
    return tables_by_idx.data();
}

/** Return a raw point to the array of force tables for each molecule */
SIRE_ALWAYS_INLINE const MolForceTable* ForceTable::data() const
{
    return tables_by_idx.constData();
}

/** Return a raw point to the array of force tables for each molecule */
SIRE_ALWAYS_INLINE const MolForceTable* ForceTable::constData() const
{
    return tables_by_idx.constData();
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

SireFF::MolForceTable operator+(const SireMaths::Vector &force,
                                const SireFF::MolForceTable &table);


SireFF::MolForceTable operator*(double value, const SireFF::MolForceTable &table);

SireFF::ForceTable operator+(const SireMaths::Vector &force,
                             const SireFF::ForceTable &table);


SireFF::ForceTable operator*(double value, const SireFF::ForceTable &table);

Q_DECLARE_METATYPE(SireFF::ForceTable);
Q_DECLARE_METATYPE(SireFF::MolForceTable);

SIRE_EXPOSE_CLASS( SireFF::ForceTable )
SIRE_EXPOSE_CLASS( SireFF::MolForceTable )

SIRE_END_HEADER

#endif
