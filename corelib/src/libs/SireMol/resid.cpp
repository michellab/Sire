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

#include "resid.h"
#include "residentifier.h"

#include "chainresid.h"
#include "groupatomids.h"
#include "groupgroupids.h"

#include "atom.h"
#include "selector.hpp"

#include "molinfo.h"

#include "mover.hpp"
#include "editor.hpp"

#include "partialmolecule.h"
#include "segment.h"
#include "chain.h"
#include "residue.h"
#include "cutgroup.h"
#include "atom.h"

#include "molecules.h"
#include "moleculegroup.h"
#include "moleculegroups.h"

#include "withres.h"

#include "SireMol/errors.h"

#include "SireStream/datastream.h"

#include "tostring.h"

using namespace SireMol;
using namespace SireID;

/** Constructor */
ResID::ResID() : ID()
{}

/** Copy constructor */
ResID::ResID(const ResID &other) : ID(other)
{}

/** Destructor */
ResID::~ResID()
{}
  
/** Return a specific residue that matches this ID */
Specify<ResID> ResID::operator[](int i) const
{
    return Specify<ResID>(*this, i);
}

/** Return a specific residue that matches this ID */
Specify<ResID> ResID::operator()(int i) const
{
    return this->operator[](i);
}

/** Return a range of residues that match this ID */
Specify<ResID> ResID::operator()(int i, int j) const
{
    return Specify<ResID>(*this, i, j);
}

/** Combine with another ID type */
IDAndSet<ResID> ResID::operator+(const ResID &other) const
{
    return IDAndSet<ResID>(*this, other);
}

/** Syntactic sugar for operator+ */
IDAndSet<ResID> ResID::operator&&(const ResID &other) const
{
    return this->operator+(other);
}

/** Syntactic sugar for operator+ */
IDAndSet<ResID> ResID::operator&(const ResID &other) const
{
    return this->operator+(other);
}

/** Combine with another ID type */
ChainResID ResID::operator+(const ChainID &other) const
{
    return ChainResID(other, *this);
}

/** Syntactic sugar for operator+ */
ChainResID ResID::operator&&(const ChainID &other) const
{
    return this->operator+(other);
}

/** Syntactic sugar for operator+ */
ChainResID ResID::operator&(const ChainID &other) const
{
    return this->operator+(other);
}

/** Combine with another ID type */
GroupAtomID<ResID,AtomID> ResID::operator+(const AtomID &other) const
{
    return GroupAtomID<ResID,AtomID>(*this, other);
}

/** Syntactic sugar for operator+ */
GroupAtomID<ResID,AtomID> ResID::operator&&(const AtomID &other) const
{
    return this->operator+(other);
}

/** Syntactic sugar for operator+ */
GroupAtomID<ResID,AtomID> ResID::operator&(const AtomID &other) const
{
    return this->operator+(other);
}

/** Combine with another ID type */
GroupGroupID<SegID,ResID> ResID::operator+(const SegID &other) const
{
    return GroupGroupID<SegID,ResID>(other, *this);
}

/** Syntactic sugar for operator+ */
GroupGroupID<SegID,ResID> ResID::operator&&(const SegID &other) const
{
    return this->operator+(other);
}

/** Syntactic sugar for operator+ */
GroupGroupID<SegID,ResID> ResID::operator&(const SegID &other) const
{
    return this->operator+(other);
}

/** Combine with another ID type */
GroupGroupID<CGID,ResID> ResID::operator+(const CGID &other) const
{
    return GroupGroupID<CGID,ResID>(other, *this);
}

/** Syntactic sugar for operator+ */
GroupGroupID<CGID,ResID> ResID::operator&&(const CGID &other) const
{
    return this->operator+(other);
}

/** Syntactic sugar for operator+ */
GroupGroupID<CGID,ResID> ResID::operator&(const CGID &other) const
{
    return this->operator+(other);
}

/** Return the match for this ID or 'other' */
IDOrSet<ResID> ResID::operator*(const ResID &other) const
{
    return IDOrSet<ResID>(*this, other);
}

/** Syntactic sugar for operator* */
IDOrSet<ResID> ResID::operator||(const ResID &other) const
{
    return this->operator*(other);
}

/** Syntactic sugar for operator* */
IDOrSet<ResID> ResID::operator|(const ResID &other) const
{
    return this->operator*(other);
}

/** Return the match for this ID or 'other' */
IDOrSet<AtomID> ResID::operator*(const AtomID &other) const
{
    return other * *this;
}

/** Syntactic sugar for operator* */
IDOrSet<AtomID> ResID::operator||(const AtomID &other) const
{
    return this->operator*(other);
}

/** Syntactic sugar for operator* */
IDOrSet<AtomID> ResID::operator|(const AtomID &other) const
{
    return this->operator*(other);
}

/** Return the match for this ID or 'other' */
IDOrSet<ResID> ResID::operator*(const ChainID &other) const
{
    return IDOrSet<ResID>(*this, other+ResID::any());
}

/** Syntactic sugar for operator* */
IDOrSet<ResID> ResID::operator||(const ChainID &other) const
{
    return this->operator*(other);
}

/** Syntactic sugar for operator* */
IDOrSet<ResID> ResID::operator|(const ChainID &other) const
{
    return this->operator*(other);
}

/** Invert this match */
InvertMatch<ResID> ResID::invert() const
{
    return InvertMatch<ResID>(*this);
}

/** Invert this match */
InvertMatch<ResID> ResID::inverse() const
{
    return this->invert();
}

/** Invert this match */
InvertMatch<ResID> ResID::operator!() const
{
    return this->invert();
}

/** Return a match for all residues */
MatchAll<ResID> ResID::any()
{
    return MatchAll<ResID>();
}

/** Match everything */
QList<ResIdx> ResID::matchAll(const MolInfo &molinfo)
{
    return molinfo.getResidues();
}

/** Return the atoms in the matching residues */
AtomsIn<ResID> ResID::atoms() const
{
    return AtomsIn<ResID>(*this);
}

/** Return a specific atom in the matching residues */
AtomsIn<ResID> ResID::atom(int i) const
{
    return AtomsIn<ResID>(*this, i);
}

/** Return a range of atoms in the matching residues */
AtomsIn<ResID> ResID::atoms(int i, int j) const
{
    return AtomsIn<ResID>(*this, i, j);
}

/** Return a Chain ID that matches chains that contain residues
    that match this Residue ID */
ChainsWithRes ResID::chains() const
{
    return ChainsWithRes(*this);
}

void ResID::processMatches(QList<ResIdx> &matches, const MolInfo&) const
{
    if (matches.isEmpty())
        throw SireMol::missing_residue( QObject::tr(
            "There are no residues that match the ID \"%1\"")
                .arg(this->toString()), CODELOC );
                
    qSort(matches);
}

/** Map this ResID to the atoms in the passed molecule view

    \throw SireMol::missing_residue
    \throw SireError::invalid_index
*/
QList<ResIdx> ResID::map(const MoleculeView &molview, const PropertyMap&) const
{
    QList<ResIdx> residxs = this->map( molview.data().info() );
    
    if (molview.selectedAll())
        return residxs;
    else 
    {
        QMutableListIterator<ResIdx> it(residxs);
        
        const AtomSelection selected_atoms = molview.selection();
        
        while (it.hasNext())
        {
            it.next();
            
            if (not selected_atoms.selected(it.value()))
                it.remove();
        }
        
        if (residxs.isEmpty())
            throw SireMol::missing_residue( QObject::tr(
                    "No atoms matching %1 can be found in the passed molecule.")
                        .arg(this->toString()), CODELOC );
                        
        return residxs;
    }
}

/** Select the atom from the passed view that matches this ID

    \throw SireMol::missing_residue
    \throw SireError::invalid_index
    \throw SireMol::duplicate_residue
*/
Residue ResID::selectFrom(const MoleculeView &molview, const PropertyMap &map) const
{
    QList<ResIdx> residxs = this->map(molview, map);
    
    if (residxs.count() > 1)
        throw SireMol::duplicate_residue( QObject::tr(
                "More than one atom matches the ID %1 (atoms %2).")
                    .arg(this->toString()).arg(Sire::toString(residxs)),
                        CODELOC );
                        
    return Residue(molview.data(), residxs.at(0)); 
}

/** Select all the atoms from the passed view that match this ID

    \throw SireMol::missing_residue
    \throw SireError::invalid_index
    \throw SireMol::duplicate_residue
*/
Selector<Residue> ResID::selectAllFrom(const MoleculeView &molview, 
                                     const PropertyMap &map) const
{
    QList<ResIdx> residxs = this->map(molview, map);
                        
    return Selector<Residue>(molview.data(), residxs); 
}

/** Return all of the atoms from the 'molecules' that match
    this ID
    
    \throw SireMol::missing_residue
*/
QHash< MolNum,Selector<Residue> >
ResID::selectAllFrom(const Molecules &molecules, const PropertyMap &map) const
{
    QHash< MolNum,Selector<Residue> > selected_atoms;
    
    //loop over all molecules...
    for (Molecules::const_iterator it = molecules.constBegin();
         it != molecules.constEnd();
         ++it)
    {
        try
        {
            //try to find this atom in this molecule
            selected_atoms.insert( it.key(), this->selectAllFrom(*it,map) );
        }
        catch(...)
        {}
    }
    
    if (selected_atoms.isEmpty())
        throw SireMol::missing_residue( QObject::tr(
            "There was no atom matching the ID \"%1\" in "
            "the set of molecules.")
                .arg(this->toString()), CODELOC );

    return selected_atoms;
}

/** Return the atom from the molecules 'molecules' that matches
    this ID
    
    \throw SireMol::missing_residue
    \throw SireMol::duplicate_residue
*/
Residue ResID::selectFrom(const Molecules &molecules, const PropertyMap &map) const
{
    QHash< MolNum,Selector<Residue> > mols = this->selectAllFrom(molecules, map);
    
    if (mols.count() > 1)
        throw SireMol::duplicate_residue( QObject::tr(
            "More than one molecule contains an atom that "
            "matches this ID (%1). These molecules have numbers %2.")
                .arg(this->toString()).arg(Sire::toString(mols.keys())),
                    CODELOC );
                    
    const Selector<Residue> &atoms = *(mols.constBegin());
    
    if (atoms.count() > 1)
        throw SireMol::duplicate_residue( QObject::tr(
            "While only one molecule (MolNum == %1) "
            "contains an atom that matches this ID (%2), it contains "
            "more than one atom that matches.")
                .arg(atoms.data().number()).arg(this->toString()),
                    CODELOC );
                    
    return atoms(0);
}

/** Return the atom from the molecule group 'molgroup' that matches
    this ID
    
    \throw SireMol::missing_residue
    \throw SireMol::duplicate_residue
*/
Residue ResID::selectFrom(const MoleculeGroup &molgroup, 
                        const PropertyMap &map) const
{
    return this->selectFrom(molgroup.molecules(), map);
}

/** Return the atoms from the molecule group 'molgroup' that match
    this ID
    
    \throw SireMol::missing_residue
*/
QHash< MolNum,Selector<Residue> >
ResID::selectAllFrom(const MoleculeGroup &molgroup,
                      const PropertyMap &map) const
{
    return this->selectAllFrom(molgroup.molecules(), map);
}

/** Return the atom from the molecule groups 'molgroups' that matches 
    this ID
    
    \throw SireMol::missing_residue
    \throw SireMol::duplicate_residue
*/
Residue ResID::selectFrom(const MolGroupsBase &molgroups,
                        const PropertyMap &map) const
{
    return this->selectFrom(molgroups.molecules(), map);
}

/** Return the set of atoms that match this ID in the molecule groups
    set 'molgroups' 
    
    \throw SireMol::missing_residue
*/
QHash< MolNum,Selector<Residue> >
ResID::selectAllFrom(const MolGroupsBase &molgroups,
                      const PropertyMap &map) const
{
    return this->selectAllFrom(molgroups.molecules(), map);
}

//fully instantiate Specify<ResID> and AtomsIn<ResID>
namespace SireID
{
    template class Specify<ResID>;
    template class IDAndSet<ResID>;
    template class IDOrSet<ResID>;
}

namespace SireMol
{
    template class AtomsIn<ResID>;
}

static const RegisterMetaType< Specify<ResID> > r_specify_resid;
static const RegisterMetaType< AtomsIn<ResID> > r_atomsin_resid;
static const RegisterMetaType< IDAndSet<ResID> > r_idandset_resid;
static const RegisterMetaType< IDOrSet<ResID> > r_idorset_resid;
