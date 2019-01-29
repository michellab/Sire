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

#ifndef SIREFF_ENERGYTABLE_H
#define SIREFF_ENERGYTABLE_H

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
class MolEnergyTable;
class EnergyTable;
}

SIREFF_EXPORT QDataStream& operator<<(QDataStream&, const SireFF::MolEnergyTable&);
SIREFF_EXPORT QDataStream& operator>>(QDataStream&, SireFF::MolEnergyTable&);

SIREFF_EXPORT QDataStream& operator<<(QDataStream&, const SireFF::EnergyTable&);
SIREFF_EXPORT QDataStream& operator>>(QDataStream&, SireFF::EnergyTable&);

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

/** This class holds the energies of all of the atoms of 
    selected CutGroups in a molecule. The MolEnergyTable is used
    to accumulate all of the energies of these atoms during
    an energy evaluation, and also to control which energies are
    evaluated (as only the energies on atoms in selected CutGroups
    are evaluated). This allows you to provide some control over
    the calculation, e.g. only placing a few protein residues into
    the energy table, thereby preventing the energy of all atoms
    in a protein from being evaluated if they aren't actually 
    necessary.
    
    @author Christopher Woods
*/
class SIREFF_EXPORT MolEnergyTable : public SireBase::PackedArray2D<SireMaths::Vector>
{

friend QDataStream& ::operator<<(QDataStream&, const MolEnergyTable&);
friend QDataStream& ::operator>>(QDataStream&, MolEnergyTable&);

public:
    typedef SireBase::PackedArray2D<SireMaths::Vector>::Array Array;

    MolEnergyTable();
    
    MolEnergyTable(const MoleculeView &molview);
    
    MolEnergyTable(const MolEnergyTable &other);
    
    ~MolEnergyTable();
    
    MolEnergyTable& operator=(const MolEnergyTable &other);
    MolEnergyTable& operator=(const Vector &force);

    bool operator==(const MolEnergyTable &other) const;
    bool operator!=(const MolEnergyTable &other) const;

    MolEnergyTable& operator+=(const MolEnergyTable &other);
    MolEnergyTable& operator-=(const MolEnergyTable &other);
    
    MolEnergyTable operator+(const MolEnergyTable &other) const;
    MolEnergyTable operator-(const MolEnergyTable &other) const;

    MolEnergyTable& operator+=(const Vector &force);
    MolEnergyTable& operator-=(const Vector &force);
    
    MolEnergyTable operator+(const Vector &force) const;
    MolEnergyTable operator-(const Vector &force) const;
    
    MolEnergyTable& operator*=(double value);
    MolEnergyTable& operator/=(double value);
    
    MolEnergyTable operator*(double value) const;
    MolEnergyTable operator/(double value) const;

    MolEnergyTable operator-() const;

    static const char* typeName();
    
    const char* what() const
    {
        return MolEnergyTable::typeName();
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

    void add(const MolEnergyTable &other);
    void subtract(const MolEnergyTable &other);
    
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

/** A EnergyTable is a workspace within which all of the energies acting 
    on the atoms of several molecules may be stored. A EnergyTable is 
    used as storing the energies requires lots of memory, and continually 
    creating a deleting such large amouts of memory would be inefficient. 
    Also, using a EnergyTable allows for energies to be accumulated directly, 
    rather than requiring intermediate storage space for the 
    individual components.

    You create an energy table to hold all of the energies of all of 
    the atoms of all of the molecules in a specified MoleculeGroup.
    The energies are held in an array that holds the energies for 
    the molecules in the same order as the molecules appear
    in the molecule group. The energytable also comes with
    an index so you can quickly look up the energies for
    a specific molecule.

    @author Christopher Woods
*/
class SIREFF_EXPORT EnergyTable
{

friend QDataStream& ::operator<<(QDataStream&, const EnergyTable&);
friend QDataStream& ::operator>>(QDataStream&, EnergyTable&);

public:
    EnergyTable();
    EnergyTable(const MoleculeGroup &molgroup);
    
    EnergyTable(const EnergyTable &other);
    
    ~EnergyTable();

    static const char* typeName();

    const char* what() const
    {
        return EnergyTable::typeName();
    }

    EnergyTable& operator=(const EnergyTable &other);
    EnergyTable& operator=(const Vector &force);

    bool operator==(const EnergyTable &other) const;
    bool operator!=(const EnergyTable &other) const;

    EnergyTable& operator+=(const EnergyTable &other);
    EnergyTable& operator-=(const EnergyTable &other);
    
    EnergyTable operator+(const EnergyTable &other) const;
    EnergyTable operator-(const EnergyTable &other) const;

    EnergyTable& operator+=(const Vector &force);
    EnergyTable& operator-=(const Vector &force);
    
    EnergyTable operator+(const Vector &force) const;
    EnergyTable operator-(const Vector &force) const;
    
    EnergyTable& operator*=(double value);
    EnergyTable& operator/=(double value);
    
    EnergyTable operator*(double value) const;
    EnergyTable operator/(double value) const;
    
    EnergyTable operator-() const;
    
    bool containsTable(MolNum molnum) const;

    void initialiseTables();

    void initialiseTable(MolNum molnum);

    MolEnergyTable& getTable(MolNum molnum);

    const MolEnergyTable& getTable(MolNum molnum) const;

    const MolEnergyTable& constGetTable(MolNum molnum) const;

    int count() const;

    const QHash<MolNum,qint32>& index() const;
    
    int indexOf(MolNum molnum) const;
    
    QList<MolNum> molNums() const;

    MolEnergyTable* data();
    const MolEnergyTable* data() const;
    const MolEnergyTable* constData() const;

    void assertContainsTableFor(MolNum molnum) const;

    void add(const EnergyTable &other);
    void subtract(const EnergyTable &other);
    
    void add(const Vector &force);
    void subtract(const Vector &force);
    
    void setAll(const Vector &force);
    
    void multiply(double value);
    void divide(double value);

private:
    /** All of the tables */
    QVector<MolEnergyTable> tables_by_idx;

    /** Index mapping from the number of the Molecule to 
        the index of its force table in the above array */
    QHash<MolNum,qint32> molnum_to_idx;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

//////
////// Inline functions for MolEnergyTable
//////

/** Return the molecule number of the molecule whose forces are
    contained in this table */
inline MolNum MolEnergyTable::molNum() const
{
    return molnum;
}

/** Return the UID of the molecular layout of the molecule
    whose forces are contained in this table */
inline const QUuid& MolEnergyTable::molUID() const
{
    return moluid;
}

/** Return the total number of CutGroups in the molecule whose
    forces are contained in this table */
inline int MolEnergyTable::nCutGroups() const
{
    return ncgroups;
}

/** Return the number of selected CutGroups in this table 
    (the number of CutGroups for which forces are held) */
inline int MolEnergyTable::nSelectedCutGroups() const
{
    return SireBase::PackedArray2D<SireMaths::Vector>::count();
}

/** Return whether or not this table contains forces for all
    of the CutGroups in the molecule */
inline bool MolEnergyTable::selectedAll() const
{
    return this->nCutGroups() == this->nSelectedCutGroups();
}

/** Return whether or not the CutGroup at index 'cgidx' has been
    selected 
    
    \throw SireError::invalid_index
*/
inline bool MolEnergyTable::selected(CGIdx cgidx) const
{
    cgidx = CGIdx(cgidx.map(this->nCutGroups()));
    
    return this->selectedAll() or cgidx_to_idx.contains(cgidx);
}

/** Return the index of the forcetable for the CutGroup at index
    'cgidx'. This returns -1 if there is no forcetable for the
    specified CutGroup
    
    \throw SireError::invalid_index
*/
inline int MolEnergyTable::map(CGIdx cgidx) const
{
    if (this->selectedAll())
        return cgidx.map( this->nCutGroups() );
    else
        return cgidx_to_idx.value( CGIdx(cgidx.map(this->nCutGroups())), -1 );
}

//////
////// Inline functions for EnergyTable
//////

/** Return whether or not this contains a table for the 
    molecule with number 'molnum' */
inline bool EnergyTable::containsTable(MolNum molnum) const
{
    return molnum_to_idx.contains(molnum);
}

/** Return the table for the molecule with number 'molnum'

    \throw SireMol::missing_molecule
*/
inline MolEnergyTable& EnergyTable::getTable(MolNum molnum)
{
    QHash<MolNum,qint32>::const_iterator it = molnum_to_idx.constFind(molnum);
    
    if (it == molnum_to_idx.constEnd())
        this->assertContainsTableFor(molnum);
        
    return tables_by_idx.data()[ it.value() ];
}

/** Return the table for the molecule with number 'molnum'

    \throw SireMol::missing_molecule
*/
inline const MolEnergyTable& EnergyTable::getTable(MolNum molnum) const
{
    QHash<MolNum,qint32>::const_iterator it = molnum_to_idx.constFind(molnum);
    
    if (it == molnum_to_idx.constEnd())
        this->assertContainsTableFor(molnum);
        
    return tables_by_idx.constData()[ it.value() ];
}

/** Return the table for the molecule with number 'molnum'

    \throw SireMol::missing_molecule
*/
inline const MolEnergyTable& EnergyTable::constGetTable(MolNum molnum) const
{
    return this->getTable(molnum);
}

/** Return the number of tables (molecules) in this object */
inline int EnergyTable::count() const
{
    return tables_by_idx.count();
}

/** Return the index used to find the index into the force tables array 
    for the force table for the molecule with a specified number */
inline const QHash<MolNum,qint32>& EnergyTable::index() const
{
    return molnum_to_idx;
}

/** Return the numbers of molecules that have force tables in this
    table */
inline QList<MolNum> EnergyTable::molNums() const
{
    return molnum_to_idx.keys();
}

/** Return a raw point to the array of force tables for each molecule */
inline MolEnergyTable* EnergyTable::data()
{
    return tables_by_idx.data();
}

/** Return a raw point to the array of force tables for each molecule */
inline const MolEnergyTable* EnergyTable::data() const
{
    return tables_by_idx.constData();
}

/** Return a raw point to the array of force tables for each molecule */
inline const MolEnergyTable* EnergyTable::constData() const
{
    return tables_by_idx.constData();
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

SireFF::MolEnergyTable operator+(const SireMaths::Vector &force,
                                const SireFF::MolEnergyTable &table);


SireFF::MolEnergyTable operator*(double value, const SireFF::MolEnergyTable &table);

SireFF::EnergyTable operator+(const SireMaths::Vector &force,
                             const SireFF::EnergyTable &table);


SireFF::EnergyTable operator*(double value, const SireFF::EnergyTable &table);

Q_DECLARE_METATYPE(SireFF::EnergyTable);
Q_DECLARE_METATYPE(SireFF::MolEnergyTable);

SIRE_EXPOSE_CLASS( SireFF::EnergyTable )
SIRE_EXPOSE_CLASS( SireFF::MolEnergyTable )

SIRE_END_HEADER

#endif
