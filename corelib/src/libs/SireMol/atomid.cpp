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

#include "atomid.h"
#include "atomidentifier.h"

#include "groupatomids.h"
#include "molatomid.h"

#include "atom.h"
#include "selector.hpp"

#include "mover.hpp"
#include "editor.hpp"

#include "segment.h"
#include "chain.h"
#include "residue.h"
#include "cutgroup.h"
#include "atom.h"

#include "molinfo.h"
#include "molecules.h"
#include "moleculegroup.h"
#include "moleculegroups.h"

#include "withatoms.h"

#include "SireMol/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "tostring.h"

#include <QDebug>

using namespace SireMol;
using namespace SireID;
using namespace SireStream;

/////////
///////// Implementation of AtomID
/////////

/** Constructor */
AtomID::AtomID() : ID()
{}

/** Copy constructor */
AtomID::AtomID(const AtomID &other) : ID(other)
{}

/** Destructor */
AtomID::~AtomID()
{}
  
/** Return a specific atom that matches this ID */
Specify<AtomID> AtomID::operator[](int i) const
{
    return Specify<AtomID>(*this, i);
}

/** Return a specific atom that matches this ID */
Specify<AtomID> AtomID::operator()(int i) const
{
    return this->operator[](i);
}

/** Return a range of atoms that match this ID */
Specify<AtomID> AtomID::operator()(int i, int j) const
{
    return Specify<AtomID>(*this, i, j);
}

/** Combine with other ID types */
IDAndSet<AtomID> AtomID::operator+(const AtomID &other) const
{
    return IDAndSet<AtomID>(*this, other);
}

/** Syntactic sugar for *this + other */
IDAndSet<AtomID> AtomID::operator&&(const AtomID &other) const
{
    return this->operator+(other);
}

/** Syntactic sugar for *this + other */
IDAndSet<AtomID> AtomID::operator&(const AtomID &other) const
{
    return this->operator+(other);
}

/** Combine with other ID types */
CGAtomID AtomID::operator+(const CGID &other) const
{
    return CGAtomID(other, *this);
}

/** Syntactic sugar for *this + other */
CGAtomID AtomID::operator&&(const CGID &other) const
{
    return this->operator+(other);
}

/** Syntactic sugar for *this + other */
CGAtomID AtomID::operator&(const CGID &other) const
{
    return this->operator+(other);
}

/** Combine with other ID types */
ResAtomID AtomID::operator+(const ResID &other) const
{
    return ResAtomID(other, *this);
}

/** Syntactic sugar for *this + other */
ResAtomID AtomID::operator&&(const ResID &other) const
{
    return this->operator+(other);
}

/** Syntactic sugar for *this + other */
ResAtomID AtomID::operator&(const ResID &other) const
{
    return this->operator+(other);
}

/** Combine with other ID types */
ChainAtomID AtomID::operator+(const ChainID &other) const
{
    return ChainAtomID(other, *this);
}

/** Syntactic sugar for *this + other */
ChainAtomID AtomID::operator&&(const ChainID &other) const
{
    return this->operator+(other);
}

/** Syntactic sugar for *this + other */
ChainAtomID AtomID::operator&(const ChainID &other) const
{
    return this->operator+(other);
}

/** Combine with other ID types */
SegAtomID AtomID::operator+(const SegID &other) const
{
    return SegAtomID(other, *this);
}

/** Syntactic sugar for *this + other */
SegAtomID AtomID::operator&&(const SegID &other) const
{
    return this->operator+(other);
}

/** Syntactic sugar for *this + other */
SegAtomID AtomID::operator&(const SegID &other) const
{
    return this->operator+(other);
}

/** Combine with other ID types */
MolAtomID AtomID::operator+(const MolID &other) const
{
    return MolAtomID(other, *this);
}

/** Syntactic sugar for *this + other */
MolAtomID AtomID::operator&&(const MolID &other) const
{
    return this->operator+(other);
}

/** Syntactic sugar for *this + other */
MolAtomID AtomID::operator&(const MolID &other) const
{
    return this->operator+(other);
}

/** Return the selection that matches this atom or 'other' */
IDOrSet<AtomID> AtomID::operator*(const AtomID &other) const
{
    return IDOrSet<AtomID>(*this, other);
}

/** Syntactic sugar for *this * other */
IDOrSet<AtomID> AtomID::operator||(const AtomID &other) const
{
    return this->operator*(other);
}

/** Syntactic sugar for *this * other */
IDOrSet<AtomID> AtomID::operator|(const AtomID &other) const
{
    return this->operator*(other);
}

/** Return the selection that matches this atom or all atoms in 'other' */
IDOrSet<AtomID> AtomID::operator*(const CGID &other) const
{
    return IDOrSet<AtomID>(*this, other+AtomID::any());
}

/** Syntactic sugar for *this * other */
IDOrSet<AtomID> AtomID::operator||(const CGID &other) const
{
    return this->operator*(other);
}

/** Syntactic sugar for *this * other */
IDOrSet<AtomID> AtomID::operator|(const CGID &other) const
{
    return this->operator*(other);
}

/** Return the selection that matches this atom or all atoms in 'other' */
IDOrSet<AtomID> AtomID::operator*(const ResID &other) const
{
    return IDOrSet<AtomID>(*this, other+AtomID::any());
}

/** Syntactic sugar for *this * other */
IDOrSet<AtomID> AtomID::operator||(const ResID &other) const
{
    return this->operator*(other);
}

/** Syntactic sugar for *this * other */
IDOrSet<AtomID> AtomID::operator|(const ResID &other) const
{
    return this->operator*(other);
}

/** Return the selection that matches this atom or all atoms in 'other' */
IDOrSet<AtomID> AtomID::operator*(const ChainID &other) const
{
    return IDOrSet<AtomID>(*this, other+AtomID::any());
}

/** Syntactic sugar for *this * other */
IDOrSet<AtomID> AtomID::operator||(const ChainID &other) const
{
    return this->operator*(other);
}

/** Syntactic sugar for *this * other */
IDOrSet<AtomID> AtomID::operator|(const ChainID &other) const
{
    return this->operator*(other);
}

/** Return the selection that matches this atom or all atoms in 'other' */
IDOrSet<AtomID> AtomID::operator*(const SegID &other) const
{
    return IDOrSet<AtomID>(*this, other+AtomID::any());
}

/** Syntactic sugar for *this * other */
IDOrSet<AtomID> AtomID::operator||(const SegID &other) const
{
    return this->operator*(other);
}

/** Syntactic sugar for *this * other */
IDOrSet<AtomID> AtomID::operator|(const SegID &other) const
{
    return this->operator*(other);
}

/** Return the selection that matches this atom or all atoms in 'other' */
IDOrSet<AtomID> AtomID::operator*(const MolID &other) const
{
    return IDOrSet<AtomID>(*this, other+AtomID::any());
}

/** Syntactic sugar for *this * other */
IDOrSet<AtomID> AtomID::operator||(const MolID &other) const
{
    return this->operator*(other);
}

/** Syntactic sugar for *this * other */
IDOrSet<AtomID> AtomID::operator|(const MolID &other) const
{
    return this->operator*(other);
}

/** Return the invert (not) of this match */
InvertMatch<AtomID> AtomID::operator!() const
{
    return InvertMatch<AtomID>(*this);
}

/** Return the invert (not) of this match */
InvertMatch<AtomID> AtomID::invert() const
{
    return InvertMatch<AtomID>(*this);
}

/** Syntactic sugar for AtomID::invert() */
InvertMatch<AtomID> AtomID::inverse() const
{
    return this->invert();
}

/** Return an AtomID that matches everything */
MatchAll<AtomID> AtomID::any()
{
    return MatchAll<AtomID>();
}

void AtomID::processMatches(QList<AtomIdx> &matches, const MolInfo&) const
{
    if (matches.isEmpty())
        throw SireMol::missing_atom( QObject::tr(
            "There are no atoms that match the ID \"%1\" in the passed molecule.")
                .arg(this->toString()), CODELOC );

    qSort(matches);
}

/** Function used by MatchAll<AtomID> to return everything */
QList<AtomIdx> AtomID::matchAll(const MolInfo &molinfo)
{
    QList<AtomIdx> idxs = molinfo.getAtoms();
    return idxs;
}

/** Map this AtomID to the atoms in the passed molecule view

    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
QList<AtomIdx> AtomID::map(const MoleculeView &molview, const PropertyMap&) const
{
    QList<AtomIdx> atomidxs = this->map( molview.data().info() );
    
    if (molview.selectedAll())
        return atomidxs;
    else 
    {
        QMutableListIterator<AtomIdx> it(atomidxs);
        
        const AtomSelection selected_atoms = molview.selection();
        
        while (it.hasNext())
        {
            it.next();
            
            if (not selected_atoms.selected(it.value()))
                it.remove();
        }
        
        if (atomidxs.isEmpty())
            throw SireMol::missing_atom( QObject::tr(
                    "No atoms matching %1 can be found in the passed molecule.")
                        .arg(this->toString()), CODELOC );
                        
        return atomidxs;
    }
}

/** Select the atom from the passed view that matches this ID

    \throw SireMol::missing_atom
    \throw SireError::invalid_index
    \throw SireMol::duplicate_atom
*/
Atom AtomID::selectFrom(const MoleculeView &molview, const PropertyMap &map) const
{
    QList<AtomIdx> atomidxs = this->map(molview, map);
    
    if (atomidxs.count() > 1)
        throw SireMol::duplicate_atom( QObject::tr(
                "More than one atom matches the ID %1 (atoms %2).")
                    .arg(this->toString()).arg(Sire::toString(atomidxs)),
                        CODELOC );
                        
    return Atom(molview.data(), atomidxs.at(0)); 
}

/** Select all the atoms from the passed view that match this ID

    \throw SireMol::missing_atom
    \throw SireError::invalid_index
    \throw SireMol::duplicate_atom
*/
Selector<Atom> AtomID::selectAllFrom(const MoleculeView &molview, 
                                     const PropertyMap &map) const
{
    QList<AtomIdx> atomidxs = this->map(molview, map);
                        
    return Selector<Atom>(molview.data(), atomidxs); 
}

/** Return all of the atoms from the 'molecules' that match
    this ID
    
    \throw SireMol::missing_atom
*/
QHash< MolNum,Selector<Atom> >
AtomID::selectAllFrom(const Molecules &molecules, const PropertyMap &map) const
{
    QHash< MolNum,Selector<Atom> > selected_atoms;
    
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
        throw SireMol::missing_atom( QObject::tr(
            "There was no atom matching the ID \"%1\" in "
            "the set of molecules.")
                .arg(this->toString()), CODELOC );

    return selected_atoms;
}

/** Return the atom from the molecules 'molecules' that matches
    this ID
    
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
*/
Atom AtomID::selectFrom(const Molecules &molecules, const PropertyMap &map) const
{
    QHash< MolNum,Selector<Atom> > mols = this->selectAllFrom(molecules, map);
    
    if (mols.count() > 1)
        throw SireMol::duplicate_atom( QObject::tr(
            "More than one molecule contains an atom that "
            "matches this ID (%1). These molecules have numbers %2.")
                .arg(this->toString()).arg(Sire::toString(mols.keys())),
                    CODELOC );
                    
    const Selector<Atom> &atoms = *(mols.constBegin());
    
    if (atoms.count() > 1)
        throw SireMol::duplicate_atom( QObject::tr(
            "While only one molecule (MolNum == %1) "
            "contains an atom that matches this ID (%2), it contains "
            "more than one atom that matches.")
                .arg(atoms.data().number()).arg(this->toString()),
                    CODELOC );
                    
    return atoms(0);
}

/** Return the atom from the molecule group 'molgroup' that matches
    this ID
    
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
*/
Atom AtomID::selectFrom(const MoleculeGroup &molgroup, 
                        const PropertyMap &map) const
{
    return this->selectFrom(molgroup.molecules(), map);
}

/** Return the atoms from the molecule group 'molgroup' that match
    this ID
    
    \throw SireMol::missing_atom
*/
QHash< MolNum,Selector<Atom> >
AtomID::selectAllFrom(const MoleculeGroup &molgroup,
                      const PropertyMap &map) const
{
    return this->selectAllFrom(molgroup.molecules(), map);
}

/** Return the atom from the molecule groups 'molgroups' that matches 
    this ID
    
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
*/
Atom AtomID::selectFrom(const MolGroupsBase &molgroups,
                        const PropertyMap &map) const
{
    return this->selectFrom(molgroups.molecules(), map);
}

/** Return the set of atoms that match this ID in the molecule groups
    set 'molgroups' 
    
    \throw SireMol::missing_atom
*/
QHash< MolNum,Selector<Atom> >
AtomID::selectAllFrom(const MolGroupsBase &molgroups,
                      const PropertyMap &map) const
{
    return this->selectAllFrom(molgroups.molecules(), map);
}

/** Return a Residue ID that matches residues that contain atoms
    that match this Atom ID */
ResWithAtoms AtomID::residues() const
{
    return ResWithAtoms(*this);
}

/** Return a CutGroup ID that matches CutGroups that contain atoms
    that match this Atom ID */
CGsWithAtoms AtomID::cutGroups() const
{
    return CGsWithAtoms(*this);
}

/** Return a Chain ID that matches chains that contain atoms
    that match this Atom ID */
ChainsWithAtoms AtomID::chains() const
{
    return ChainsWithAtoms(*this);
}

/** Return a Segment ID that matches segments that contain atoms
    that match this Atom ID */
SegsWithAtoms AtomID::segments() const
{
    return SegsWithAtoms(*this);
}

/////////
///////// Implementation of AtomIdx
/////////

static const RegisterMetaType<AtomIdx> r_atomidx;

/** Serialise to a binary datastream */
QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds, const AtomIdx &atomidx)
{
    writeHeader(ds, r_atomidx, 1);
    
    ds << static_cast<const SireID::Index_T_<AtomIdx>&>(atomidx);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds, AtomIdx &atomidx)
{
    VersionID v = readHeader(ds, r_atomidx);
    
    if (v == 1)
    {
        ds >> static_cast<SireID::Index_T_<AtomIdx>&>(atomidx);
    }
    else
        throw version_error( v, "1", r_atomidx, CODELOC );
        
    return ds;
}

QList<AtomIdx> AtomIdx::map(const MolInfo &molinfo) const
{
    return molinfo.map(*this);
}

AtomIdx::AtomIdx() : SireID::Index_T_<AtomIdx>(), AtomID()
{}

AtomIdx::AtomIdx(qint32 idx) : SireID::Index_T_<AtomIdx>(idx), AtomID()
{}

AtomIdx::AtomIdx(const AtomIdx &other) : SireID::Index_T_<AtomIdx>(other), AtomID(other)
{}

AtomIdx::~AtomIdx()
{}

AtomIdx AtomIdx::null()
{
    return AtomIdx();
}

bool AtomIdx::isNull() const
{
    return SireID::Index_T_<AtomIdx>::isNull();
}

uint AtomIdx::hash() const
{
    return SireID::Index_T_<AtomIdx>::hash();
}

QString AtomIdx::toString() const
{
    return QString("AtomIdx(%1)").arg(_idx);
}

AtomIdx& AtomIdx::operator=(const AtomIdx &other)
{
    SireID::IndexBase::operator=(other);
    AtomID::operator=(other);
    return *this;
}

bool AtomIdx::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare<AtomIdx>(*this, other);
}

const char* AtomIdx::typeName()
{
    return QMetaType::typeName( qMetaTypeId<AtomIdx>() );
}

/////////
///////// Implementation of AtomNum
/////////

static const RegisterMetaType<AtomNum> r_atomnum;

/** Serialise to a binary datastream */
QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds, const AtomNum &atomnum)
{
    writeHeader(ds, r_atomnum, 1);
    
    ds << static_cast<const SireID::Number&>(atomnum);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds, AtomNum &atomnum)
{
    VersionID v = readHeader(ds, r_atomnum);
    
    if (v == 1)
    {
        ds >> static_cast<SireID::Number&>(atomnum);
    }
    else
        throw version_error( v, "1", r_atomnum, CODELOC );
        
    return ds;
}

AtomNum::AtomNum() : SireID::Number(), AtomID()
{}

AtomNum::AtomNum(quint32 num) : SireID::Number(num), AtomID()
{}

AtomNum::AtomNum(const AtomNum &other) : SireID::Number(other), AtomID(other)
{}

AtomNum::~AtomNum()
{}

bool AtomNum::isNull() const
{
    return Number::isNull();
}

uint AtomNum::hash() const
{
    return ::qHash( static_cast<const Number&>(*this) );
}

QString AtomNum::toString() const
{
    return QString("AtomNum(%1)").arg(_num);
}

AtomNum& AtomNum::operator=(const AtomNum &other)
{
    SireID::Number::operator=(other);
    AtomID::operator=(other);
    return *this;
}

bool AtomNum::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare<AtomNum>(*this, other);
}

bool AtomNum::operator==(const AtomNum &other) const
{
    return _num == other._num;
}

bool AtomNum::operator!=(const AtomNum &other) const
{
    return _num != other._num;
}

bool AtomNum::operator<(const AtomNum &other) const
{
    return _num < other._num;
}

bool AtomNum::operator<=(const AtomNum &other) const
{
    return _num <= other._num;
}

bool AtomNum::operator>(const AtomNum &other) const
{
    return _num > other._num;
}

bool AtomNum::operator>=(const AtomNum &other) const
{
    return _num >= other._num;
}

QList<AtomIdx> AtomNum::map(const MolInfo &molinfo) const
{
    return molinfo.map(*this);
}

const char* AtomNum::typeName()
{
    return QMetaType::typeName( qMetaTypeId<AtomNum>() );
}

/////////
///////// Implementation of AtomName
/////////

static const RegisterMetaType<AtomName> r_atomname;

/** Serialise to a binary datastream */
QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds, const AtomName &atomname)
{
    writeHeader(ds, r_atomname, 1);
    
    ds << static_cast<const SireID::Name&>(atomname);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds, AtomName &atomname)
{
    VersionID v = readHeader(ds, r_atomname);
    
    if (v == 1)
    {
        ds >> static_cast<SireID::Name&>(atomname);
    }
    else
        throw version_error( v, "1", r_atomname, CODELOC );
        
    return ds;
}

AtomName::AtomName() : SireID::Name(), AtomID()
{}

AtomName::AtomName(const QString &name) : SireID::Name(name), AtomID()
{}

AtomName::AtomName(const QString &name, SireID::CaseSensitivity case_sensitivity)
         : SireID::Name(name, case_sensitivity), AtomID()
{}

AtomName::AtomName(const AtomName &other) : SireID::Name(other), AtomID(other)
{}

AtomName::~AtomName()
{}

bool AtomName::isNull() const
{
    return SireID::Name::isNull();
}

uint AtomName::hash() const
{
    return ::qHash(_name);
}

QString AtomName::toString() const
{
    if (case_sensitive)
        return QString("AtomName('%1')").arg(_name);
    else
        return QString("AtomName('%1', isCaseSensitive=False)").arg(_name);
}

AtomName& AtomName::operator=(const AtomName &other)
{
    SireID::Name::operator=(other);
    AtomID::operator=(other);
    return *this;
}

bool AtomName::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare<AtomName>(*this, other);
}

bool AtomName::operator==(const AtomName &other) const
{
    return SireID::Name::operator==(other);
}

bool AtomName::operator!=(const AtomName &other) const
{
    return SireID::Name::operator!=(other);
}

QList<AtomIdx> AtomName::map(const MolInfo &molinfo) const
{
    return molinfo.map(*this);
}

const char* AtomName::typeName()
{
    return QMetaType::typeName( qMetaTypeId<AtomName>() );
}

//fully instantiate Specify<AtomID>
namespace SireID
{
    template class Specify<AtomID>;
    template class IDAndSet<AtomID>;
    template class IDOrSet<AtomID>;
}

static const RegisterMetaType< Specify<AtomID> > r_specify_atomid;
static const RegisterMetaType< IDAndSet<AtomID> > r_idandset_atomid;
static const RegisterMetaType< IDOrSet<AtomID> > r_idorset_atomid;

AtomNum* AtomNum::clone() const
{
    return new AtomNum(*this);
}


AtomName* AtomName::clone() const
{
    return new AtomName(*this);
}


AtomIdx* AtomIdx::clone() const
{
    return new AtomIdx(*this);
}

