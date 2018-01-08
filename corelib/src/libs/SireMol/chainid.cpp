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

#include "chainid.h"
#include "chainidentifier.h"

#include "chainresid.h"
#include "groupatomids.h"
#include "groupgroupids.h"

#include "atom.h"
#include "selector.hpp"

#include "mover.hpp"
#include "editor.hpp"

#include "molinfo.h"

#include "partialmolecule.h"
#include "segment.h"
#include "chain.h"
#include "residue.h"
#include "cutgroup.h"
#include "atom.h"

#include "molecules.h"
#include "moleculegroup.h"
#include "moleculegroups.h"

#include "SireMol/errors.h"

#include "tostring.h"

using namespace SireMol;
using namespace SireID;

/** Constructor */
ChainID::ChainID() : ID()
{}

/** Copy constructor */
ChainID::ChainID(const ChainID &other) : ID(other)
{}

/** Destructor */
ChainID::~ChainID()
{}

/** Combine with another ID object */
IDAndSet<ChainID> ChainID::operator+(const ChainID &other) const
{
    return IDAndSet<ChainID>(*this, other);
}

/** Syntactic sugar for operator+ */
IDAndSet<ChainID> ChainID::operator&&(const ChainID &other) const
{
    return this->operator+(other);
}

/** Syntactic sugar for operator+ */
IDAndSet<ChainID> ChainID::operator&(const ChainID &other) const
{
    return this->operator+(other);
}

/** Combine with another ID object */
ChainResID ChainID::operator+(const ResID &other) const
{
    return ChainResID(*this, other);
}

/** Syntactic sugar for operator+ */
ChainResID ChainID::operator&&(const ResID &other) const
{
    return this->operator+(other);
}

/** Syntactic sugar for operator+ */
ChainResID ChainID::operator&(const ResID &other) const
{
    return this->operator+(other);
}

/** Combine with another ID object */
GroupAtomID<ChainID,AtomID> ChainID::operator+(const AtomID &other) const
{
    return GroupAtomID<ChainID,AtomID>(*this, other);
}

/** Syntactic sugar for operator+ */
GroupAtomID<ChainID,AtomID> ChainID::operator&&(const AtomID &other) const
{
    return this->operator+(other);
}

/** Syntactic sugar for operator+ */
GroupAtomID<ChainID,AtomID> ChainID::operator&(const AtomID &other) const
{
    return this->operator+(other);
}

/** Combine with another ID object */
GroupGroupID<SegID,ChainID> ChainID::operator+(const SegID &other) const
{
    return GroupGroupID<SegID,ChainID>(other, *this);
}

/** Syntactic sugar for operator+ */
GroupGroupID<SegID,ChainID> ChainID::operator&&(const SegID &other) const
{
    return this->operator+(other);
}

/** Syntactic sugar for operator+ */
GroupGroupID<SegID,ChainID> ChainID::operator&(const SegID &other) const
{
    return this->operator+(other);
}

/** Combine with another ID object */
GroupGroupID<CGID,ChainID> ChainID::operator+(const CGID &other) const
{
    return GroupGroupID<CGID,ChainID>(other, *this);
}

/** Syntactic sugar for operator+ */
GroupGroupID<CGID,ChainID> ChainID::operator&&(const CGID &other) const
{
    return this->operator+(other);
}

/** Syntactic sugar for operator+ */
GroupGroupID<CGID,ChainID> ChainID::operator&(const CGID &other) const
{
    return this->operator+(other);
}

/** Return the combination of this ID or other */
IDOrSet<ChainID> ChainID::operator*(const ChainID &other) const
{
    return IDOrSet<ChainID>(*this, other);
}

/** Syntactic sugar for operator* */
IDOrSet<ChainID> ChainID::operator||(const ChainID &other) const
{
    return this->operator*(other);
}

/** Syntactic sugar for operator* */
IDOrSet<ChainID> ChainID::operator|(const ChainID &other) const
{
    return this->operator*(other);
}

/** Return the combination of this ID or other */
IDOrSet<ResID> ChainID::operator*(const ResID &other) const
{
    return other * *this;
}

/** Syntactic sugar for operator* */
IDOrSet<ResID> ChainID::operator||(const ResID &other) const
{
    return this->operator*(other);
}

/** Syntactic sugar for operator* */
IDOrSet<ResID> ChainID::operator|(const ResID &other) const
{
    return this->operator*(other);
}

/** Return the combination of this ID or other */
IDOrSet<AtomID> ChainID::operator*(const AtomID &other) const
{
    return other * *this;
}

/** Syntactic sugar for operator* */
IDOrSet<AtomID> ChainID::operator||(const AtomID &other) const
{
    return this->operator*(other);
}

/** Syntactic sugar for operator* */
IDOrSet<AtomID> ChainID::operator|(const AtomID &other) const
{
    return this->operator*(other);
}

/** Inverse this match */
InvertMatch<ChainID> ChainID::invert() const
{
    return InvertMatch<ChainID>(*this);
}

/** Inverse this match */
InvertMatch<ChainID> ChainID::inverse() const
{
    return this->invert();
}

/** Inverse this match */
InvertMatch<ChainID> ChainID::operator!() const
{
    return this->invert();
}

/** Return this and not other */
IDAndSet<ChainID> ChainID::operator-(const ChainID &other) const
{
    return this->operator+(other.inverse());
}

/** Return this and not other */
ChainResID ChainID::operator-(const ResID &other) const
{
    return this->operator+(other.inverse());
}

/** Return this and not other */
GroupAtomID<ChainID,AtomID> ChainID::operator-(const AtomID &other) const
{
    return this->operator+(other.inverse());
}

/** Return this and not other */
GroupGroupID<SegID,ChainID> ChainID::operator-(const SegID &other) const
{
    return this->operator+(other.inverse());
}

/** Return this and not other */
GroupGroupID<CGID,ChainID> ChainID::operator-(const CGID &other) const
{
    return this->operator+(other.inverse());
}

/** Return not this */
SireID::InvertMatch<ChainID> ChainID::operator-() const
{
    return InvertMatch<ChainID>(*this);
}

/** Return a match for any chains */
MatchAll<ChainID> ChainID::any()
{
    return MatchAll<ChainID>();
}

/** Match everything */
QList<ChainIdx> ChainID::matchAll(const MolInfo &molinfo)
{
    return molinfo.getChains();
}

/** Return the atoms in the matching residues */
AtomsIn<ChainID> ChainID::atoms() const
{
    return AtomsIn<ChainID>(*this);
}

/** Return a specific atom in the matching residues */
AtomsIn<ChainID> ChainID::atom(int i) const
{
    return AtomsIn<ChainID>(*this, i);
}

/** Return a range of atoms in the matching residues */
AtomsIn<ChainID> ChainID::atoms(int i, int j) const
{
    return AtomsIn<ChainID>(*this, i, j);
}

/** Return the atoms in the matching residues */
ResIn<ChainID> ChainID::residues() const
{
    return ResIn<ChainID>(*this);
}

/** Return a specific atom in the matching residues */
ResIn<ChainID> ChainID::residue(int i) const
{
    return ResIn<ChainID>(*this, i);
}

/** Return a range of atoms in the matching residues */
ResIn<ChainID> ChainID::residues(int i, int j) const
{
    return ResIn<ChainID>(*this, i, j);
}

/** Return a specific object that matches this ID */
Specify<ChainID> ChainID::operator[](qint64 i) const
{
    return Specify<ChainID>(*this, i);
}

/** Return a range of objects that match this ID */
Specify<ChainID> ChainID::operator[](const SireBase::Range &range) const
{
    return Specify<ChainID>(*this, range);
}

/** Return a range of objects that match this ID */
Specify<ChainID> ChainID::operator()(const SireBase::Range &range) const
{
    return Specify<ChainID>(*this, range);
}

/** Return a specific object that matches this ID */
Specify<ChainID> ChainID::operator()(qint64 i) const
{
    return this->operator[](i);
}

/** Return a range of objects that match this ID */
Specify<ChainID> ChainID::operator()(qint64 start, qint64 end) const
{
    return Specify<ChainID>(*this, start, end);
}

/** Return a range of objects that match this ID */
Specify<ChainID> ChainID::operator()(qint64 start, qint64 end, qint64 increment) const
{
    return Specify<ChainID>(*this, start, end, increment);
}

void ChainID::processMatches(QList<ChainIdx> &matches, const MolInfo&) const
{
    if (matches.isEmpty())
        throw SireMol::missing_chain( QObject::tr(
            "There are no chains that match the ID \"%1\" in the passed molecule.")
                .arg(this->toString()), CODELOC );

    qSort(matches);
}

/** Map this ChainID to the chains in the passed molecule view

    \throw SireMol::missing_chain
    \throw SireError::invalid_index
*/
QList<ChainIdx> ChainID::map(const MoleculeView &molview, const PropertyMap&) const
{
    QList<ChainIdx> chainidxs = this->map( molview.data().info() );
    
    if (molview.selectedAll())
        return chainidxs;
    else 
    {
        QMutableListIterator<ChainIdx> it(chainidxs);
        
        const AtomSelection selected_chains = molview.selection();
        
        while (it.hasNext())
        {
            it.next();
            
            if (not selected_chains.selected(it.value()))
                it.remove();
        }
        
        if (chainidxs.isEmpty())
            throw SireMol::missing_chain( QObject::tr(
                    "No chains matching %1 can be found in the passed molecule.")
                        .arg(this->toString()), CODELOC );
                        
        return chainidxs;
    }
}

/** Select the chain from the passed view that matches this ID

    \throw SireMol::missing_chain
    \throw SireError::invalid_index
    \throw SireMol::duplicate_chain
*/
Chain ChainID::selectFrom(const MoleculeView &molview, const PropertyMap &map) const
{
    QList<ChainIdx> chainidxs = this->map(molview, map);
    
    if (chainidxs.count() > 1)
        throw SireMol::duplicate_chain( QObject::tr(
                "More than one chain matches the ID %1 (chains %2).")
                    .arg(this->toString()).arg(Sire::toString(chainidxs)),
                        CODELOC );
                        
    return Chain(molview.data(), chainidxs.at(0)); 
}

/** Select all the chains from the passed view that match this ID

    \throw SireMol::missing_chain
    \throw SireError::invalid_index
    \throw SireMol::duplicate_chain
*/
Selector<Chain> ChainID::selectAllFrom(const MoleculeView &molview, 
                                     const PropertyMap &map) const
{
    QList<ChainIdx> chainidxs = this->map(molview, map);
                        
    return Selector<Chain>(molview.data(), chainidxs); 
}

/** Return all of the chains from the 'molecules' that match
    this ID
    
    \throw SireMol::missing_chain
*/
QHash< MolNum,Selector<Chain> >
ChainID::selectAllFrom(const Molecules &molecules, const PropertyMap &map) const
{
    QHash< MolNum,Selector<Chain> > selected_chains;
    
    //loop over all molecules...
    for (Molecules::const_iterator it = molecules.constBegin();
         it != molecules.constEnd();
         ++it)
    {
        try
        {
            //try to find this chain in this molecule
            selected_chains.insert( it.key(), this->selectAllFrom(*it,map) );
        }
        catch(...)
        {}
    }
    
    if (selected_chains.isEmpty())
        throw SireMol::missing_chain( QObject::tr(
            "There was no chain matching the ID \"%1\" in "
            "the set of molecules.")
                .arg(this->toString()), CODELOC );

    return selected_chains;
}

/** Return the chain from the molecules 'molecules' that matches
    this ID
    
    \throw SireMol::missing_chain
    \throw SireMol::duplicate_chain
*/
Chain ChainID::selectFrom(const Molecules &molecules, const PropertyMap &map) const
{
    QHash< MolNum,Selector<Chain> > mols = this->selectAllFrom(molecules, map);
    
    if (mols.count() > 1)
        throw SireMol::duplicate_chain( QObject::tr(
            "More than one molecule contains an chain that "
            "matches this ID (%1). These molecules have numbers %2.")
                .arg(this->toString()).arg(Sire::toString(mols.keys())),
                    CODELOC );
                    
    const Selector<Chain> &chains = *(mols.constBegin());
    
    if (chains.count() > 1)
        throw SireMol::duplicate_chain( QObject::tr(
            "While only one molecule (MolNum == %1) "
            "contains an chain that matches this ID (%2), it contains "
            "more than one chain that matches.")
                .arg(chains.data().number()).arg(this->toString()),
                    CODELOC );
                    
    return chains(0);
}

/** Return the chain from the molecule group 'molgroup' that matches
    this ID
    
    \throw SireMol::missing_chain
    \throw SireMol::duplicate_chain
*/
Chain ChainID::selectFrom(const MoleculeGroup &molgroup, 
                        const PropertyMap &map) const
{
    return this->selectFrom(molgroup.molecules(), map);
}

/** Return the chains from the molecule group 'molgroup' that match
    this ID
    
    \throw SireMol::missing_chain
*/
QHash< MolNum,Selector<Chain> >
ChainID::selectAllFrom(const MoleculeGroup &molgroup,
                      const PropertyMap &map) const
{
    return this->selectAllFrom(molgroup.molecules(), map);
}

/** Return the chain from the molecule groups 'molgroups' that matches 
    this ID
    
    \throw SireMol::missing_chain
    \throw SireMol::duplicate_chain
*/
Chain ChainID::selectFrom(const MolGroupsBase &molgroups,
                        const PropertyMap &map) const
{
    return this->selectFrom(molgroups.molecules(), map);
}

/** Return the set of chains that match this ID in the molecule groups
    set 'molgroups' 
    
    \throw SireMol::missing_chain
*/
QHash< MolNum,Selector<Chain> >
ChainID::selectAllFrom(const MolGroupsBase &molgroups,
                      const PropertyMap &map) const
{
    return this->selectAllFrom(molgroups.molecules(), map);
}

//fully instantiate template classes
namespace SireID
{
    template class Specify<ChainID>;
    template class IDAndSet<ChainID>;
    template class IDOrSet<ChainID>;
}

namespace SireMol
{
    template class AtomsIn<ChainID>;
    template class ResIn<ChainID>;
}

static const RegisterMetaType< Specify<ChainID> > r_specify_chainid;
static const RegisterMetaType< AtomsIn<ChainID> > r_atomsin_chainid;
static const RegisterMetaType< ResIn<ChainID> > r_resin_chainid;
static const RegisterMetaType< IDAndSet<ChainID> > r_idandset_chainid;
static const RegisterMetaType< IDOrSet<ChainID> > r_idorset_chainid;
