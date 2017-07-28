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

#include "segid.h"
#include "segidentifier.h"

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

#include "SireMol/errors.h"

#include "SireStream/datastream.h"

#include "tostring.h"

using namespace SireMol;
using namespace SireID;

/** Constructor */
SegID::SegID() : ID()
{}

/** Copy constructor */
SegID::SegID(const SegID &other) : ID(other)
{}

/** Destructor */
SegID::~SegID()
{}
  
/** Return a specific atom that matches this ID */
Specify<SegID> SegID::operator[](int i) const
{
    return Specify<SegID>(*this, i);
}

/** Return a specific atom that matches this ID */
Specify<SegID> SegID::operator()(int i) const
{
    return this->operator[](i);
}

/** Return a range of atoms that match this ID */
Specify<SegID> SegID::operator()(int i, int j) const
{
    return Specify<SegID>(*this, i, j);
}

/** Combine two ID types */
IDAndSet<SegID> SegID::operator+(const SegID &other) const
{
    return IDAndSet<SegID>(*this, other);
}

/** Syntactic sugar for operator+ */
IDAndSet<SegID> SegID::operator&&(const SegID &other) const
{
    return this->operator+(other);
}

/** Syntactic sugar for operator+ */
IDAndSet<SegID> SegID::operator&(const SegID &other) const
{
    return this->operator+(other);
}

/** Combine two ID types */
GroupAtomID<SegID,AtomID> SegID::operator+(const AtomID &other) const
{
    return GroupAtomID<SegID,AtomID>(*this, other);
}

/** Syntactic sugar for operator+ */
GroupAtomID<SegID,AtomID> SegID::operator&&(const AtomID &other) const
{
    return this->operator+(other);
}

/** Syntactic sugar for operator+ */
GroupAtomID<SegID,AtomID> SegID::operator&(const AtomID &other) const
{
    return this->operator+(other);
}

/** Combine two ID types */
GroupGroupID<SegID,CGID> SegID::operator+(const CGID &other) const
{
    return GroupGroupID<SegID,CGID>(*this, other);
}

/** Syntactic sugar for operator+ */
GroupGroupID<SegID,CGID> SegID::operator&&(const CGID &other) const
{
    return this->operator+(other);
}

/** Syntactic sugar for operator+ */
GroupGroupID<SegID,CGID> SegID::operator&(const CGID &other) const
{
    return this->operator+(other);
}

/** Combine two ID types */
GroupGroupID<SegID,ResID> SegID::operator+(const ResID &other) const
{
    return GroupGroupID<SegID,ResID>(*this, other);
}

/** Syntactic sugar for operator+ */
GroupGroupID<SegID,ResID> SegID::operator&&(const ResID &other) const
{
    return this->operator+(other);
}

/** Syntactic sugar for operator+ */
GroupGroupID<SegID,ResID> SegID::operator&(const ResID &other) const
{
    return this->operator+(other);
}

/** Combine two ID types */
GroupGroupID<SegID,ChainID> SegID::operator+(const ChainID &other) const
{
    return GroupGroupID<SegID,ChainID>(*this, other);
}

/** Syntactic sugar for operator+ */
GroupGroupID<SegID,ChainID> SegID::operator&&(const ChainID &other) const
{
    return this->operator+(other);
}

/** Syntactic sugar for operator+ */
GroupGroupID<SegID,ChainID> SegID::operator&(const ChainID &other) const
{
    return this->operator+(other);
}

/** Return the object to search for the match of this or 'other' */
IDOrSet<SegID> SegID::operator*(const SegID &other) const
{
    return IDOrSet<SegID>(*this, other);
}

/** Syntactic sugar for operator* */
IDOrSet<SegID> SegID::operator||(const SegID &other) const
{
    return this->operator*(other);
}

/** Syntactic sugar for operator* */
IDOrSet<SegID> SegID::operator|(const SegID &other) const
{
    return this->operator*(other);
}

/** Return the atoms in the matching residues */
AtomsIn<SegID> SegID::atoms() const
{
    return AtomsIn<SegID>(*this);
}

/** Return a specific atom in the matching residues */
AtomsIn<SegID> SegID::atom(int i) const
{
    return AtomsIn<SegID>(*this, i);
}

/** Return a range of atoms in the matching residues */
AtomsIn<SegID> SegID::atoms(int i, int j) const
{
    return AtomsIn<SegID>(*this, i, j);
}

void SegID::processMatches(QList<SegIdx> &matches, const MolInfo &molinfo) const
{
    if (matches.isEmpty())
        throw SireMol::missing_segment( QObject::tr(
            "There are no segments that match the ID \"%1\".")
                .arg(this->toString()), CODELOC );
                
    qSort(matches);
}

/** Map this SegID to the atoms in the passed molecule view

    \throw SireMol::missing_segment
    \throw SireError::invalid_index
*/
QList<SegIdx> SegID::map(const MoleculeView &molview, const PropertyMap&) const
{
    QList<SegIdx> segidxs = this->map( molview.data().info() );
    
    if (molview.selectedAll())
        return segidxs;
    else 
    {
        QMutableListIterator<SegIdx> it(segidxs);
        
        const AtomSelection selected_atoms = molview.selection();
        
        while (it.hasNext())
        {
            it.next();
            
            if (not selected_atoms.selected(it.value()))
                it.remove();
        }
        
        if (segidxs.isEmpty())
            throw SireMol::missing_segment( QObject::tr(
                    "No atoms matching %1 can be found in the passed molecule.")
                        .arg(this->toString()), CODELOC );
                        
        return segidxs;
    }
}

/** Select the atom from the passed view that matches this ID

    \throw SireMol::missing_segment
    \throw SireError::invalid_index
    \throw SireMol::duplicate_segment
*/
Segment SegID::selectFrom(const MoleculeView &molview, const PropertyMap &map) const
{
    QList<SegIdx> segidxs = this->map(molview, map);
    
    if (segidxs.count() > 1)
        throw SireMol::duplicate_segment( QObject::tr(
                "More than one atom matches the ID %1 (atoms %2).")
                    .arg(this->toString()).arg(Sire::toString(segidxs)),
                        CODELOC );
                        
    return Segment(molview.data(), segidxs.at(0)); 
}

/** Select all the atoms from the passed view that match this ID

    \throw SireMol::missing_segment
    \throw SireError::invalid_index
    \throw SireMol::duplicate_segment
*/
Selector<Segment> SegID::selectAllFrom(const MoleculeView &molview, 
                                     const PropertyMap &map) const
{
    QList<SegIdx> segidxs = this->map(molview, map);
                        
    return Selector<Segment>(molview.data(), segidxs); 
}

/** Return all of the atoms from the 'molecules' that match
    this ID
    
    \throw SireMol::missing_segment
*/
QHash< MolNum,Selector<Segment> >
SegID::selectAllFrom(const Molecules &molecules, const PropertyMap &map) const
{
    QHash< MolNum,Selector<Segment> > selected_atoms;
    
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
        throw SireMol::missing_segment( QObject::tr(
            "There was no atom matching the ID \"%1\" in "
            "the set of molecules.")
                .arg(this->toString()), CODELOC );

    return selected_atoms;
}

/** Return the atom from the molecules 'molecules' that matches
    this ID
    
    \throw SireMol::missing_segment
    \throw SireMol::duplicate_segment
*/
Segment SegID::selectFrom(const Molecules &molecules, const PropertyMap &map) const
{
    QHash< MolNum,Selector<Segment> > mols = this->selectAllFrom(molecules, map);
    
    if (mols.count() > 1)
        throw SireMol::duplicate_segment( QObject::tr(
            "More than one molecule contains an atom that "
            "matches this ID (%1). These molecules have numbers %2.")
                .arg(this->toString()).arg(Sire::toString(mols.keys())),
                    CODELOC );
                    
    const Selector<Segment> &atoms = *(mols.constBegin());
    
    if (atoms.count() > 1)
        throw SireMol::duplicate_segment( QObject::tr(
            "While only one molecule (MolNum == %1) "
            "contains an atom that matches this ID (%2), it contains "
            "more than one atom that matches.")
                .arg(atoms.data().number()).arg(this->toString()),
                    CODELOC );
                    
    return atoms(0);
}

/** Return the atom from the molecule group 'molgroup' that matches
    this ID
    
    \throw SireMol::missing_segment
    \throw SireMol::duplicate_segment
*/
Segment SegID::selectFrom(const MoleculeGroup &molgroup, 
                        const PropertyMap &map) const
{
    return this->selectFrom(molgroup.molecules(), map);
}

/** Return the atoms from the molecule group 'molgroup' that match
    this ID
    
    \throw SireMol::missing_segment
*/
QHash< MolNum,Selector<Segment> >
SegID::selectAllFrom(const MoleculeGroup &molgroup,
                      const PropertyMap &map) const
{
    return this->selectAllFrom(molgroup.molecules(), map);
}

/** Return the atom from the molecule groups 'molgroups' that matches 
    this ID
    
    \throw SireMol::missing_segment
    \throw SireMol::duplicate_segment
*/
Segment SegID::selectFrom(const MolGroupsBase &molgroups,
                        const PropertyMap &map) const
{
    return this->selectFrom(molgroups.molecules(), map);
}

/** Return the set of atoms that match this ID in the molecule groups
    set 'molgroups' 
    
    \throw SireMol::missing_segment
*/
QHash< MolNum,Selector<Segment> >
SegID::selectAllFrom(const MolGroupsBase &molgroups,
                      const PropertyMap &map) const
{
    return this->selectAllFrom(molgroups.molecules(), map);
}

//fully instantiate template classes
namespace SireID
{
    template class Specify<SegID>;
    template class IDAndSet<SegID>;
    template class IDOrSet<SegID>;
}

namespace SireMol
{
    template class AtomsIn<SegID>;
}

static const RegisterMetaType< Specify<SegID> > r_specify_segid;
static const RegisterMetaType< AtomsIn<SegID> > r_atomsin_segid;
static const RegisterMetaType< IDAndSet<SegID> > r_idandset_segid;
static const RegisterMetaType< IDOrSet<SegID> > r_idorset_segid;
