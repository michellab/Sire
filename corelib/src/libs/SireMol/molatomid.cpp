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

#include "molatomid.h"

#include "selector.hpp"
#include "atom.h"

#include "moleculegroup.h"
#include "moleculegroups.h"

#include "SireMol/errors.h"
#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMol;
using namespace SireID;
using namespace SireStream;

using boost::tuples::tuple;

#include <QDebug>

static const RegisterMetaType<MolAtomID> r_molatomid;

/** Serialise to a binary datastream */
QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds, const MolAtomID &molatomid)
{
    writeHeader(ds, r_molatomid, 1);
    
    SharedDataStream sds(ds);
    
    sds << molatomid.molid << molatomid.atomid;
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds, MolAtomID &molatomid)
{
    VersionID v = readHeader(ds, r_molatomid);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> molatomid.molid >> molatomid.atomid;
    }
    else
        throw version_error( v, "1", r_molatomid, CODELOC );
        
    return ds;
}

/** Collapse any nested IDs together */
void MolAtomID::collapse()
{
    const MolAtomID *molatomid = dynamic_cast<const MolAtomID*>( &(atomid.base()) );

    if ( molatomid != 0 )
    {
        if (not molatomid->molID().isNull())
        {
            if (molid.isNull())
                molid = molatomid->molID();
                
            else
                molid = IDAndSet<MolID>( molid, molatomid->molID() );
        }
        
        atomid = molatomid->atomID();
    }
}

/** Construct a MolAtomID that matches everything */
MolAtomID::MolAtomID() : AtomID()
{}

/** Construct a MolAtomID that matches the atoms identified by 'atomid'
    in the molecules identified by 'molid' */
MolAtomID::MolAtomID(const MolID &mol_id, const AtomID &atom_id)
          : AtomID(), molid(mol_id), atomid(atom_id)
{}

/** Construct a MolAtomID that matches the atoms identified by 'atomid'
    in the molecules identified by 'molid' */
MolAtomID::MolAtomID(const AtomID &atom_id, const MolID &mol_id)
          : AtomID(), molid(mol_id), atomid(atom_id)
{}

/** Construct a MolAtomID that matches the specified atoms in the specified
    molecules */
MolAtomID::MolAtomID(const tuple<MolIdentifier,AtomIdentifier> &molatomid)
          : AtomID(), molid(molatomid.get<0>()), atomid(molatomid.get<1>())
{}

/** Construct a MolAtomID that matches the specified atoms in the specified
    molecules */
MolAtomID::MolAtomID(const boost::tuple<AtomIdentifier,MolIdentifier> &molatomid)
          : AtomID(), molid(molatomid.get<1>()), atomid(molatomid.get<0>())
{}

/** Copy constructor */
MolAtomID::MolAtomID(const MolAtomID &other)
          : AtomID(other), molid(other.molid), atomid(other.atomid)
{}

/** Destructor */
MolAtomID::~MolAtomID()
{}

/** Return whether or not this is null */
bool MolAtomID::isNull() const
{
    return molid.isNull() and atomid.isNull();
}

/** Return a hash of this ID */
uint MolAtomID::hash() const
{
    return (molid.hash() << 16) | (atomid.hash() & 0x0000FFFF);
}
            
/** Return a string representation of this ID */
QString MolAtomID::toString() const
{
    if (atomid.isNull())
        return QObject::tr("Atoms in %1").arg(molid.toString());
    
    else if (molid.isNull())
        return atomid.toString();
        
    else
        return QObject::tr("%1 and %2")
                            .arg(molid.toString(), atomid.toString());
}

/** Return the AtomID part of this match */
const AtomID& MolAtomID::atomID() const
{
    return atomid.base();
}

/** Return the MolID part of this match */
const MolID& MolAtomID::molID() const
{
    return molid.base();
}

/** Copy assignment operator */
MolAtomID& MolAtomID::operator=(const MolAtomID &other)
{
    molid = other.molid;
    atomid = other.atomid;
    return *this;
}

/** Comparison operator */
bool MolAtomID::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare<MolAtomID>(*this, other);
}

/** Comparison operator */
bool MolAtomID::operator==(const MolAtomID &other) const
{
    return molid == other.molid and atomid == other.atomid;
}

/** Comparison operator */
bool MolAtomID::operator!=(const MolAtomID &other) const
{
    return molid != other.molid or atomid != other.atomid;
}

/** Map this ID to the indicies of the matching atoms 

    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
QList<AtomIdx> MolAtomID::map(const MolInfo &molinfo) const
{
    //there isn't enough information to check that the molecule
    //matches the MolID part of this ID!!!
    return atomid.map(molinfo);
}

QHash< MolNum,Selector<Atom> > 
MolAtomID::selectAllFrom(const Molecules &molecules, const PropertyMap &map) const
{
    QHash< MolNum,Selector<Atom> > selected_atoms;
    
    QList<MolNum> molnums = molid.map(molecules);
    
    foreach (MolNum molnum, molnums)
    {
        const ViewsOfMol &mol = molecules[molnum];

        try
        {
            //try to find this atom in this molecule
            selected_atoms.insert( molnum,
                                   mol.selectAll(*this,map) );
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

QHash< MolNum,Selector<Atom> > 
MolAtomID::selectAllFrom(const MoleculeGroup &molgroup, const PropertyMap &map) const
{
    QHash< MolNum,Selector<Atom> > selected_atoms;
    
    QList<MolNum> molnums = molid.map(molgroup);
    
    foreach (MolNum molnum, molnums)
    {
        const ViewsOfMol &mol = molgroup[molnum];

        try
        {
            //try to find this atom in this molecule
            selected_atoms.insert( molnum,
                                   mol.selectAll(*this,map) );
        }
        catch(...)
        {}
    }

    if (selected_atoms.isEmpty())
        throw SireMol::missing_atom( QObject::tr(
            "There was no atom matching the ID \"%1\" in "
            "the MoleculeGroup %2.")
                .arg(this->toString(), molgroup.toString()), CODELOC );

    return selected_atoms;
}

QHash< MolNum,Selector<Atom> > 
MolAtomID::selectAllFrom(const MolGroupsBase &molgroups, const PropertyMap &map) const
{
    QHash< MolNum,Selector<Atom> > selected_atoms;
    
    QList<MolNum> molnums = molid.map(molgroups);
    
    foreach (MolNum molnum, molnums)
    {
        ViewsOfMol mol = molgroups[molnum];

        try
        {
            //try to find this atom in this molecule
            selected_atoms.insert( molnum,
                                   mol.selectAll(*this,map) );
        }
        catch(...)
        {}
    }

    if (selected_atoms.isEmpty())
        throw SireMol::missing_atom( QObject::tr(
            "There was no atom matching the ID \"%1\" in "
            "the passed MoleculeGroups.")
                .arg(this->toString()), CODELOC );

    return selected_atoms;
}

const char* MolAtomID::typeName()
{
    return QMetaType::typeName( qMetaTypeId<MolAtomID>() );
}

MolAtomID* MolAtomID::clone() const
{
    return new MolAtomID(*this);
}

