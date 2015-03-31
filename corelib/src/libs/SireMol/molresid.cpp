/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2014  Christopher Woods
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

#include "molresid.h"

#include "selector.hpp"
#include "residue.h"

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

static const RegisterMetaType<MolResID> r_molresid;

/** Serialise to a binary datastream */
QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds, const MolResID &molresid)
{
    writeHeader(ds, r_molresid, 1);
    
    SharedDataStream sds(ds);
    
    sds << molresid.molid << molresid.resid;
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds, MolResID &molresid)
{
    VersionID v = readHeader(ds, r_molresid);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> molresid.molid >> molresid.resid;
    }
    else
        throw version_error( v, "1", r_molresid, CODELOC );
        
    return ds;
}

/** Collapse any nested IDs together */
void MolResID::collapse()
{
    const MolResID *molresid = dynamic_cast<const MolResID*>( &(resid.base()) );

    if ( molresid != 0 )
    {
        if (not molresid->molID().isNull())
        {
            if (molid.isNull())
                molid = molresid->molID();
                
            else
                molid = IDAndSet<MolID>( molid, molresid->molID() );
        }
        
        resid = molresid->resID();
    }
}

/** Construct a MolResID that matches everything */
MolResID::MolResID() : ResID()
{}

/** Construct a MolResID that matches the residues identified by 'resid'
    in the molecules identified by 'molid' */
MolResID::MolResID(const MolID &mol_id, const ResID &res_id)
          : ResID(), molid(mol_id), resid(res_id)
{}

/** Construct a MolResID that matches the residues identified by 'resid'
    in the molecules identified by 'molid' */
MolResID::MolResID(const ResID &res_id, const MolID &mol_id)
          : ResID(), molid(mol_id), resid(res_id)
{}

/** Construct a MolResID that matches the specified atoms in the specified
    molecules */
MolResID::MolResID(const tuple<MolIdentifier,ResIdentifier> &id)
          : ResID(), molid(id.get<0>()), resid(id.get<1>())
{}

/** Construct a MolResID that matches the specified atoms in the specified
    molecules */
MolResID::MolResID(const boost::tuple<ResIdentifier,MolIdentifier> &id)
          : ResID(), molid(id.get<1>()), resid(id.get<0>())
{}

/** Copy constructor */
MolResID::MolResID(const MolResID &other)
          : ResID(other), molid(other.molid), resid(other.resid)
{}

/** Destructor */
MolResID::~MolResID()
{}

/** Return whether or not this is null */
bool MolResID::isNull() const
{
    return molid.isNull() and resid.isNull();
}

/** Return a hash of this ID */
uint MolResID::hash() const
{
    return (molid.hash() << 16) | (resid.hash() & 0x0000FFFF);
}
            
/** Return a string representation of this ID */
QString MolResID::toString() const
{
    if (resid.isNull())
        return QObject::tr("Residues in %1").arg(molid.toString());
    
    else if (molid.isNull())
        return resid.toString();
        
    else
        return QObject::tr("%1 and %2")
                            .arg(molid.toString(), resid.toString());
}

/** Return the ResID part of this match */
const ResID& MolResID::resID() const
{
    return resid.base();
}

/** Return the MolID part of this match */
const MolID& MolResID::molID() const
{
    return molid.base();
}

/** Copy assignment operator */
MolResID& MolResID::operator=(const MolResID &other)
{
    molid = other.molid;
    resid = other.resid;
    return *this;
}

/** Comparison operator */
bool MolResID::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare<MolResID>(*this, other);
}

/** Comparison operator */
bool MolResID::operator==(const MolResID &other) const
{
    return molid == other.molid and resid == other.resid;
}

/** Comparison operator */
bool MolResID::operator!=(const MolResID &other) const
{
    return molid != other.molid or resid != other.resid;
}

/** Map this ID to the indicies of the matching residues

    \throw SireMol::missing_residue
    \throw SireError::invalid_index
*/
QList<ResIdx> MolResID::map(const MolInfo &molinfo) const
{
    //there isn't enough information to check that the molecule
    //matches the MolID part of this ID!!!
    return resid.map(molinfo);
}

QHash< MolNum,Selector<Residue> >
MolResID::selectAllFrom(const Molecules &molecules, const PropertyMap &map) const
{
    QHash< MolNum,Selector<Residue> > selected_res;
    
    QList<MolNum> molnums = molid.map(molecules);
    
    foreach (MolNum molnum, molnums)
    {
        const ViewsOfMol &mol = molecules[molnum];

        try
        {
            //try to find this atom in this molecule
            selected_res.insert( molnum,
                                 mol.selectAll(*this,map) );
        }
        catch(...)
        {}
    }

    if (selected_res.isEmpty())
        throw SireMol::missing_residue( QObject::tr(
            "There was no residue matching the ID \"%1\" in "
            "the set of molecules.")
                .arg(this->toString()), CODELOC );

    return selected_res;
}

QHash< MolNum,Selector<Residue> >
MolResID::selectAllFrom(const MoleculeGroup &molgroup, const PropertyMap &map) const
{
    QHash< MolNum,Selector<Residue> > selected_res;
    
    QList<MolNum> molnums = molid.map(molgroup);
    
    foreach (MolNum molnum, molnums)
    {
        const ViewsOfMol &mol = molgroup[molnum];

        try
        {
            //try to find this atom in this molecule
            selected_res.insert( molnum,
                                 mol.selectAll(*this,map) );
        }
        catch(...)
        {}
    }

    if (selected_res.isEmpty())
        throw SireMol::missing_residue( QObject::tr(
            "There was no residue matching the ID \"%1\" in "
            "the MoleculeGroup %2.")
                .arg(this->toString(), molgroup.toString()), CODELOC );

    return selected_res;
}

QHash< MolNum,Selector<Residue> >
MolResID::selectAllFrom(const MolGroupsBase &molgroups, const PropertyMap &map) const
{
    QHash< MolNum,Selector<Residue> > selected_res;
    
    QList<MolNum> molnums = molid.map(molgroups);
    
    foreach (MolNum molnum, molnums)
    {
        ViewsOfMol mol = molgroups[molnum];

        try
        {
            //try to find this atom in this molecule
            selected_res.insert( molnum,
                                 mol.selectAll(*this,map) );
        }
        catch(...)
        {}
    }

    if (selected_res.isEmpty())
        throw SireMol::missing_residue( QObject::tr(
            "There was no residue matching the ID \"%1\" in "
            "the passed MoleculeGroups.")
                .arg(this->toString()), CODELOC );

    return selected_res;
}

const char* MolResID::typeName()
{
    return QMetaType::typeName( qMetaTypeId<MolResID>() );
}

MolResID* MolResID::clone() const
{
    return new MolResID(*this);
}

/////////////
///////////// Implementation of MolResNum
/////////////

static const RegisterMetaType<MolResNum> r_molresnum;

/** Serialise to a binary datastream */
QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds, const MolResNum &molresnum)
{
    writeHeader(ds, r_molresnum, 1);
    
    ds << molresnum.molnum << molresnum.resnum;
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds, MolResNum &molresnum)
{
    VersionID v = readHeader(ds, r_molresnum);
    
    if (v == 1)
    {
        ds >> molresnum.molnum >> molresnum.resnum;
    }
    else
        throw version_error( v, "1", r_molresnum, CODELOC );
        
    return ds;
}

/** Construct a MolResNum that matches everything */
MolResNum::MolResNum() : ResID()
{}

/** Construct a MolResNum that matches the residue 'resnum' in the molecule 'molnum' */
MolResNum::MolResNum(const MolNum &mol_num, const ResNum &res_num)
          : ResID(), molnum(mol_num), resnum(res_num)
{}

/** Construct a MolResNum that matches the residue 'resnum' in the molecule 'molnum' */
MolResNum::MolResNum(const ResNum &res_num, const MolNum &mol_num)
          : ResID(), molnum(mol_num), resnum(res_num)
{}

/** Copy constructor */
MolResNum::MolResNum(const MolResNum &other)
          : ResID(other), molnum(other.molnum), resnum(other.resnum)
{}

MolResID::MolResID(const MolResNum &molresnum)
         : ResID(), molid(molresnum.molNum()), resid(molresnum.resNum())
{}

/** Destructor */
MolResNum::~MolResNum()
{}

/** Return whether or not this is null */
bool MolResNum::isNull() const
{
    return molnum.isNull() and resnum.isNull();
}

/** Return a hash of this ID */
uint MolResNum::hash() const
{
    return (molnum.hash() << 16) | (resnum.hash() & 0x0000FFFF);
}
            
/** Return a string representation of this ID */
QString MolResNum::toString() const
{
    if (resnum.isNull())
        return QObject::tr("Residues in %1").arg(molnum.toString());
    
    else if (molnum.isNull())
        return resnum.toString();
        
    else
        return QObject::tr("%1 and %2")
                            .arg(molnum.toString(), resnum.toString());
}

/** Return the ResNum part of this match */
const ResNum& MolResNum::resNum() const
{
    return resnum;
}

/** Return the MolNum part of this match */
const MolNum& MolResNum::molNum() const
{
    return molnum;
}

/** Copy assignment operator */
MolResNum& MolResNum::operator=(const MolResNum &other)
{
    molnum = other.molnum;
    resnum = other.resnum;
    return *this;
}

/** Comparison operator */
bool MolResNum::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare<MolResNum>(*this, other);
}

/** Comparison operator */
bool MolResNum::operator==(const MolResNum &other) const
{
    return molnum == other.molnum and resnum == other.resnum;
}

/** Comparison operator */
bool MolResNum::operator!=(const MolResNum &other) const
{
    return not MolResNum::operator==(other);
}

/** Map this ID to the indicies of the matching residues

    \throw SireMol::missing_residue
    \throw SireError::invalid_index
*/
QList<ResIdx> MolResNum::map(const MolInfo &molinfo) const
{
    //there isn't enough information to check that the molecule
    //matches the MolID part of this ID!!!
    return resnum.map(molinfo);
}

QHash< MolNum,Selector<Residue> >
MolResNum::selectAllFrom(const Molecules &molecules, const PropertyMap &map) const
{
    return MolResID(*this).selectAllFrom(molecules, map);
}

QHash< MolNum,Selector<Residue> >
MolResNum::selectAllFrom(const MoleculeGroup &molgroup, const PropertyMap &map) const
{
    return MolResID(*this).selectAllFrom(molgroup, map);
}

QHash< MolNum,Selector<Residue> >
MolResNum::selectAllFrom(const MolGroupsBase &molgroups, const PropertyMap &map) const
{
    return MolResID(*this).selectAllFrom(molgroups, map);
}

const char* MolResNum::typeName()
{
    return QMetaType::typeName( qMetaTypeId<MolResNum>() );
}

MolResNum* MolResNum::clone() const
{
    return new MolResNum(*this);
}
