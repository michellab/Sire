/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2008  Christopher Woods
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

#include "ffmolecules.h"

#include "SireMol/atomselection.h"
#include "SireMol/molnum.h"
#include "SireMol/moleculedata.h"
#include "SireMol/moleculeinfodata.h"
#include "SireMol/molecule.h"

#include "SireMol/mover.hpp"

#include "SireError/errors.h"
#include "SireMol/errors.h"

#include "tostring.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>

using namespace SireFF;
using namespace SireFF::detail;

using namespace SireStream;


void 
SireFF::detail::throwFFMoleculesIncompatibleParameterNames(const PropertyMap &old_map, 
                                                           const PropertyMap &new_map)
{
    throw SireError::incompatible_error( QObject::tr(
        "You cannot change the parameter source names when you add "
        "atoms to a molecule. The existing molecule uses the parameter "
        "source names %1, while you are trying to add atoms using "
        "the parameter sources names %2.")
            .arg(old_map.toString(), new_map.toString()), CODELOC );
}

void 
SireFF::detail::assertCompatible(const FFMoleculeBase &mol0,
                                 const FFMoleculeBase &mol1)
{
    if (mol0.number() != mol1.number())
        throw SireError::incompatible_error( QObject::tr(
            "You cannot construct a ChangedMolecule from two different "
            "molecules! (%1 vs. %2)")
                .arg(mol0.number()).arg(mol1.number()),
                    CODELOC );
}

/////////
///////// Implementation of FFMoleculeBase
/////////

static const RegisterMetaType<FFMoleculeBase> r_ffmolbase(MAGIC_ONLY, NO_ROOT,
                                                  "SireFF::detail::FFMoleculeBase");

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                      const FFMoleculeBase &ffmol)
{
    writeHeader(ds, r_ffmolbase, 1);
    
    SharedDataStream sds(ds);

    sds << ffmol.mol << ffmol.idx_to_cgidx;
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds,
                                      FFMoleculeBase &ffmol)
{
    VersionID v = readHeader(ds, r_ffmolbase);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> ffmol.mol >> ffmol.idx_to_cgidx;
    }
    else
        throw version_error(v, "1", r_ffmolbase, CODELOC);
        
    return ds;
}

/** Null constructor */
FFMoleculeBase::FFMoleculeBase()
{}

/** Construct a specialised version of the view of the molecule
    in 'molview' */
FFMoleculeBase::FFMoleculeBase(const PartialMolecule &molview)
           : mol(molview)
{
    if (not molview.selection().selectedAllCutGroups())
        idx_to_cgidx = molview.selection().selectedCutGroups();
}

/** Destructor */
FFMoleculeBase::~FFMoleculeBase()
{}

/** Copy assignment operator */
FFMoleculeBase& FFMoleculeBase::operator=(const FFMoleculeBase &other)
{
    mol = other.mol;
    idx_to_cgidx = other.idx_to_cgidx;
    return *this;
}

/** Comparison operator */
bool FFMoleculeBase::operator==(const FFMoleculeBase &other) const
{
    return mol == other.mol;
}

/** Comparison operator */
bool FFMoleculeBase::operator!=(const FFMoleculeBase &other) const
{
    return mol != other.mol;
}

/** Return a reference to the non-specialised molecule */
const PartialMolecule& FFMoleculeBase::molecule() const
{
    return mol;
}

/** Return the number of the molecule */
MolNum FFMoleculeBase::number() const
{
    return mol.number();
}

/** Return whether or not this molecule is empty (contains no atoms) */
bool FFMoleculeBase::isEmpty() const
{
    return mol.nAtoms() == 0;
}

/** Return the CGIdx of the ith CutGroup that is selected as
    part of this molecule
    
    \throw SireError::invalid_index
*/
CGIdx FFMoleculeBase::cgIdx(quint32 i) const
{
    if (this->isEmpty())
        throw SireError::invalid_index( QObject::tr(
            "Invalid index (%1) as there are no CutGroups selected "
            "as part of this molecule!").arg(i), CODELOC );

    else if (idx_to_cgidx.isEmpty())
    {
        //all of the CutGroups are selected
        if (i >= quint32(mol.data().info().nCutGroups()))
            throw SireError::invalid_index( QObject::tr(
                "Invalid index (%1) as the number of CutGroups in "
                "this molecule is only %2!")
                    .arg(i).arg(mol.data().info().nCutGroups()), CODELOC );
                    
        return CGIdx(i);
    }
    else
    {
        if (i >= quint32(idx_to_cgidx.count()))
            throw SireError::invalid_index( QObject::tr(
                "Invalid index (%1) as the number of CutGroups "
                "selected as part of this molecule is only %2!")
                    .arg(i).arg(idx_to_cgidx.count()), CODELOC );
                    
        return idx_to_cgidx.at(i);
    }
}

/** Return the number of CutGroups in this view of the molecule
    (this is the number of CutGroups that contain at least
    one selected atom) */
int FFMoleculeBase::nCutGroups() const
{
    return mol.nCutGroups();
}

/** Change this view to use the same molecule data version
    as the passed molecule. 
    
    Note that you can only change the molecule layout ID of the
    molecule if the entire molecule is included in this
    group. If only part of the molecule is included then
    it is not possible, as the code is not capable of
    updating an atom selection to a new molecule layout ID.

    \throw SireError::incompatible_error
*/
bool FFMoleculeBase::change(const SireMol::Molecule &molecule)
{
    if (molecule.number() != mol.number())
        throw SireError::incompatible_error( QObject::tr(
            "You cannot change the molecule %1 using the molecule %2 "
            "as they are different molecules! (their numbers are "
            "different)")
                .arg(mol.number()).arg(molecule.number()), CODELOC );

    if (mol.data().version() != molecule.data().version())
    {
        if (mol.selection().selectedAll())
            mol = molecule;
            
        else if (mol.data().info().UID() != molecule.data().info().UID())
        {
            throw SireError::incompatible_error( QObject::tr(
                "Cannot change the molecule %1 as its molecule layout "
                "ID has changed (%2 to %3), and only part of the molecule "
                "exists in the forcefield.")
                    .arg(mol.number()).arg(mol.data().info().UID().toString())
                    .arg(molecule.data().info().UID().toString()), CODELOC );
        }
    
        mol = PartialMolecule(molecule.data(), mol.selection());
        return true;
    }
    else
        return false;
}

/** Add the selected atoms 'added_atoms' to this view. This does
    nothing if these atoms are already part of the view.
    
    Returns whether this changes the view.
    
    \throw SireError::incompatible_error
*/
bool FFMoleculeBase::add(const AtomSelection &added_atoms)
{
    if (mol.selection().contains(added_atoms))
        return false;
        
    mol = PartialMolecule(mol.data(), mol.selection() + added_atoms);

    if (mol.selection().selectedAllCutGroups())
        idx_to_cgidx = QList<CGIdx>();
    else
        idx_to_cgidx = mol.selection().selectedCutGroups();

    return true;
}

/** Remove the selected atoms 'removed_atoms' from this view.
    This does nothing if none of these atoms are part of this view.
    
    \throw SireError::incompatible_error
*/
bool FFMoleculeBase::remove(const AtomSelection &removed_atoms)
{
    if (mol.selection().selectedNone())
        return false;

    else if (not mol.selection().intersects(removed_atoms))
        return false;
    
    else if (removed_atoms.contains(mol.selection()))
    {
        //we are removing the entire molecule
        mol = PartialMolecule();
        idx_to_cgidx = QList<CGIdx>();
        return true;
    }
    else
    {
        mol = PartialMolecule(mol.data(), mol.selection() - removed_atoms);

        if (mol.selection().selectedAllCutGroups())
            idx_to_cgidx = QList<CGIdx>();
        else
            idx_to_cgidx = mol.selection().selectedCutGroups();

        return true;
    }
}

/** Private function used to restore the old molecule in case of 
    an exception being thrown (so to restore the original state) */
void FFMoleculeBase::restore(const PartialMolecule &oldmol)
{
    mol = oldmol;
}

////////
//////// Implementation of FFMoleculesBase
////////

static const RegisterMetaType<FFMoleculesBase> r_ffmolsbase( MAGIC_ONLY, NO_ROOT,
                                                "SireFF::detail::FFMoleculesBase" );
                                                
/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                      const FFMoleculesBase &ffmols)
{
    writeHeader(ds, r_ffmolsbase, 1);
    
    SharedDataStream sds(ds);
    sds << ffmols.molnums_by_idx
        << ffmols.parameter_names;
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds,
                                      FFMoleculesBase &ffmols)
{
    VersionID v = readHeader(ds, r_ffmolsbase);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> ffmols.molnums_by_idx
            >> ffmols.parameter_names;
        
        ffmols._pvt_reindex();
    }
    else
        throw version_error(v, "1", r_ffmolsbase, CODELOC);
        
    return ds;
}

/** Null constructor */
FFMoleculesBase::FFMoleculesBase()
{}

/** Construct space for all of the contained molecules */
FFMoleculesBase::FFMoleculesBase(const MoleculeGroup &molgroup, const PropertyMap &map)
{
    int nmols = molgroup.nMolecules();

    parameter_names = QVector<PropertyMap>(nmols, map);

    molnums_by_idx = QVector<MolNum>(nmols);
    idxs_by_molnum.reserve(nmols);
    
    quint32 i = 0;
    MolNum *molnums_by_idx_array = molnums_by_idx.data();
    
    for (MoleculeGroup::const_iterator it = molgroup.constBegin();
         it != molgroup.constEnd();
         ++it)
    {
        molnums_by_idx_array[i] = it->number();
        idxs_by_molnum.insert( it->number(), i );
        ++i;
    }
}

/** Copy constructor */
FFMoleculesBase::FFMoleculesBase(const FFMoleculesBase &other)
                : parameter_names(other.parameter_names),
                  molnums_by_idx(other.molnums_by_idx),
                  idxs_by_molnum(other.idxs_by_molnum)
{}

/** Destructor */
FFMoleculesBase::~FFMoleculesBase()
{}

/** Copy assignment operator */
FFMoleculesBase& FFMoleculesBase::operator=(const FFMoleculesBase &other)
{
    parameter_names = other.parameter_names;
    molnums_by_idx = other.molnums_by_idx;
    idxs_by_molnum = other.idxs_by_molnum;
    return *this;
}

/** Comparison operator */
bool FFMoleculesBase::operator==(const FFMoleculesBase &other) const
{
    return molnums_by_idx == other.molnums_by_idx;
}

/** Comparison operator */
bool FFMoleculesBase::operator!=(const FFMoleculesBase &other) const
{
    return molnums_by_idx != other.molnums_by_idx;
}

/** Return the number of molecules in this group */
int FFMoleculesBase::count() const
{
    return molnums_by_idx.count();
}

/** Return whether or not this collection of molecules is empty */
bool FFMoleculesBase::isEmpty() const
{
    return molnums_by_idx.isEmpty();
}

/** Return the array that maps the index in the group to the 
    molecule number */
const QVector<MolNum> FFMoleculesBase::molNumsByIndex() const
{
    return molnums_by_idx;
}

/** Return the hash that maps molecule numbers to their index 
    in the array */
const QHash<MolNum,quint32>& FFMoleculesBase::indexesByMolNum() const
{
    return idxs_by_molnum;
}

/** Return the index of the molecule with number 'molnum'

    \throw SireMol::missing_molecule
*/
quint32 FFMoleculesBase::indexOf(MolNum molnum) const
{
    QHash<MolNum,quint32>::const_iterator it = idxs_by_molnum.constFind(molnum);
    
    if (it == idxs_by_molnum.constEnd())
        throw SireMol::missing_molecule( QObject::tr(
            "There is no molecule with number %1 in this group. "
            "The available molecules are %2.")
                .arg(molnum).arg( Sire::toString(idxs_by_molnum.keys()) ),
                    CODELOC );
                    
    return *it;
}

/** Completely remove all molecules */
void FFMoleculesBase::clear()
{
    molnums_by_idx.clear();
    idxs_by_molnum.clear();
}

/** Return whether or not this group contains the molecule with 
    number 'molnum' */
bool FFMoleculesBase::contains(MolNum molnum) const
{
    return idxs_by_molnum.contains(molnum);
}

/** Assert that this group contains the molecule with number 'molnum'

    \throw SireMol::missing_molecule
*/
void FFMoleculesBase::assertContains(MolNum molnum) const
{
    if (not idxs_by_molnum.contains(molnum))
        throw SireMol::missing_molecule( QObject::tr(
            "There is no molecule with number %1 in this group. "
            "The available molecules are %2.")
                .arg(molnum).arg( Sire::toString(idxs_by_molnum.keys()) ),
                    CODELOC );
}

/** Add the molecule 'molnum' to the array, returning its index */
quint32 FFMoleculesBase::_pvt_add(MolNum molnum)
{
    QHash<MolNum,quint32>::const_iterator it = idxs_by_molnum.constFind(molnum);
    
    if (it == idxs_by_molnum.constEnd())
    {
        //this is a new molecule - add it onto the end of the array
        molnums_by_idx.append(molnum);
        idxs_by_molnum.insert(molnum, molnums_by_idx.count() - 1);
        
        return molnums_by_idx.count() - 1;
    }
    else
        return *it;
}
    
/** Return the array of the names of the properties used to get the parameters
    for each molecule, arranged in the order that the molecules appear
    in this group */
const QVector<PropertyMap>& FFMoleculesBase::parameterNamesByIndex() const
{
    return parameter_names;
}


/** Reindex this group */
void FFMoleculesBase::_pvt_reindex()
{
    idxs_by_molnum.clear();
    
    if (molnums_by_idx.isEmpty())
        return;
        
    quint32 nmols = molnums_by_idx.count();
    const MolNum *molnums_array = molnums_by_idx.constData();
    
    idxs_by_molnum.reserve(nmols);
    
    for (quint32 i=0; i<nmols; ++i)
    {
        idxs_by_molnum.insert(molnums_array[i], i);
    }
}

/** Completely remove the molecule with number 'molnum' from this index.
    This does nothing if this molecule is not in this group */
void FFMoleculesBase::_pvt_remove(MolNum molnum)
{
    QHash<MolNum,quint32>::const_iterator it = idxs_by_molnum.constFind(molnum);
    
    if (it != idxs_by_molnum.constEnd())
    {
        quint32 idx = *it;
        
        idxs_by_molnum.remove(molnum);
        
        molnums_by_idx.remove(idx);
        
        if (molnums_by_idx.isEmpty())
        {
            molnums_by_idx.clear();
            idxs_by_molnum.clear();
        }
        else if (idx != quint32(molnums_by_idx.count()))
        {
            //we need to reindex the arrays...
            this->_pvt_reindex();
        }
    }
}

/** Remove all of the molecules whose numbers are in 'molnums' */
void FFMoleculesBase::_pvt_remove(const QList<MolNum> &molnums)
{
    if (molnums.isEmpty())
        return;
    else if (molnums.count() == 1)
    {
        this->_pvt_remove( *(molnums.begin()) );
        return;
    }
    
    QHash<MolNum,quint32>::const_iterator it;
    bool removed_some = false;
    
    foreach (MolNum molnum, molnums)
    {
        it = idxs_by_molnum.constFind(molnum);
        
        if (it != idxs_by_molnum.constEnd())
        {
            molnums_by_idx.remove(*it);
            idxs_by_molnum.remove(molnum);
            removed_some = true;
        }
    }
    
    if (removed_some)
        this->_pvt_reindex();
}
