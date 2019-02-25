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

#include "withatoms.h"

#include "SireMol/errors.h"

#include "SireStream/datastream.h"

using namespace SireMol;
using namespace SireID;
using namespace SireStream;

/////////
///////// Implementation of ResWithAtoms
/////////

static const RegisterMetaType<ResWithAtoms> r_resid;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const ResWithAtoms &resid)
{
    writeHeader(ds, r_resid, 1);
    
    ds << resid.atomid;
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, ResWithAtoms &resid)
{
    VersionID v = readHeader(ds, r_resid);
    
    if (v == 1)
    {
        ds >> resid.atomid;
    }
    else
        throw version_error( v, "1", r_resid, CODELOC );
        
    return ds;
}

/** Null constructor */
ResWithAtoms::ResWithAtoms() : ResID()
{}

/** Construct from the passed AtomID */
ResWithAtoms::ResWithAtoms(const AtomID &id)
              : ResID(), atomid(id)
{}

/** Copy constructor */
ResWithAtoms::ResWithAtoms(const ResWithAtoms &other)
              : ResID(other), atomid(other.atomid)
{}

/** Destructor */
ResWithAtoms::~ResWithAtoms()
{}

/** Is this selection null? */
bool ResWithAtoms::isNull() const
{
    return atomid.isNull();
}

/** Return a hash of this identifier */
uint ResWithAtoms::hash() const
{
    return atomid.hash();
}
            
/** Return a string representatio of this ID */
QString ResWithAtoms::toString() const
{
    return QObject::tr("ResWithAtoms( %1 )").arg(atomid.toString());
}

/** Return the atom ID */
const AtomID& ResWithAtoms::atomID() const
{
    return atomid.base();
}

/** Copy assignment operator */
ResWithAtoms& ResWithAtoms::operator=(const ResWithAtoms &other)
{
    atomid = other.atomid;
    return *this;
}

/** Comparison operator */
bool ResWithAtoms::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare<ResWithAtoms>(*this, other);
}

/** Comparison operator */
bool ResWithAtoms::operator==(const ResWithAtoms &other) const
{
    return atomid == other.atomid;
}

/** Comparison operator */
bool ResWithAtoms::operator!=(const ResWithAtoms &other) const
{
    return not ResWithAtoms::operator==(other);
}

/** Map this ID to the list of indicies of residues that match this ID

    \throw SireMol::missing_atom
    \throw SireMol::missing_residue
    \throw SireError::invalid_index
*/
QList<ResIdx> ResWithAtoms::map(const MolInfo &molinfo) const
{
    QSet<ResIdx> resids;
    QList<ResIdx> residxs;
    
    foreach (AtomIdx atomidx, molinfo.map(atomid))
    {
        try
        {
            ResIdx residx = molinfo.parentResidue(atomidx);
        
            if (not resids.contains(residx))
            {
                residxs.append(residx);
                resids.insert(residx);
            }
        }
        catch(...)
        {}
    }
    
    if (residxs.isEmpty())
        throw SireMol::missing_residue( QObject::tr(
                "There are no residues that contain atoms that match "
                "the Atom ID \"%1\".")
                    .arg(atomid.toString()), CODELOC );
    
    qSort(residxs);
    
    return residxs;
}

const char* ResWithAtoms::typeName()
{
    return QMetaType::typeName( qMetaTypeId<ResWithAtoms>() );
}

ResWithAtoms* ResWithAtoms::clone() const
{
    return new ResWithAtoms(*this);
}

/////////
///////// Implementation of CGsWithAtoms
/////////

static const RegisterMetaType<CGsWithAtoms> r_cgid;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const CGsWithAtoms &cgid)
{
    writeHeader(ds, r_cgid, 1);
    
    ds << cgid.atomid;
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, CGsWithAtoms &cgid)
{
    VersionID v = readHeader(ds, r_cgid);
    
    if (v == 1)
    {
        ds >> cgid.atomid;
    }
    else
        throw version_error( v, "1", r_cgid, CODELOC );
        
    return ds;
}

/** Null constructor */
CGsWithAtoms::CGsWithAtoms() : CGID()
{}

/** Construct from the passed AtomID */
CGsWithAtoms::CGsWithAtoms(const AtomID &id)
             : CGID(), atomid(id)
{}

/** Copy constructor */
CGsWithAtoms::CGsWithAtoms(const CGsWithAtoms &other)
             : CGID(other), atomid(other.atomid)
{}

/** Destructor */
CGsWithAtoms::~CGsWithAtoms()
{}

/** Is this selection null? */
bool CGsWithAtoms::isNull() const
{
    return atomid.isNull();
}

/** Return a hash of this identifier */
uint CGsWithAtoms::hash() const
{
    return atomid.hash();
}
            
/** Return a string representatio of this ID */
QString CGsWithAtoms::toString() const
{
    return QObject::tr( "CGsWithAtoms( %1 )" ).arg(atomid.toString());
}

/** Return the atom ID */
const AtomID& CGsWithAtoms::atomID() const
{
    return atomid.base();
}

/** Copy assignment operator */
CGsWithAtoms& CGsWithAtoms::operator=(const CGsWithAtoms &other)
{
    atomid = other.atomid;
    return *this;
}

/** Comparison operator */
bool CGsWithAtoms::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare<CGsWithAtoms>(*this, other);
}

/** Comparison operator */
bool CGsWithAtoms::operator==(const CGsWithAtoms &other) const
{
    return atomid == other.atomid;
}

/** Comparison operator */
bool CGsWithAtoms::operator!=(const CGsWithAtoms &other) const
{
    return not CGsWithAtoms::operator==(other);
}

/** Map this ID to the list of indicies of CutGroups that match this ID

    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
*/
QList<CGIdx> CGsWithAtoms::map(const MolInfo &molinfo) const
{
    QSet<CGIdx> cgids;
    QList<CGIdx> cgidxs;
    
    foreach (AtomIdx atomidx, molinfo.map(atomid))
    {
        try
        {
            CGIdx cgidx = molinfo.parentCutGroup(atomidx);
        
            if (not cgids.contains(cgidx))
            {
                cgidxs.append(cgidx);
                cgids.insert(cgidx);
            }
        }
        catch(...)
        {}
    }
    
    if (cgidxs.isEmpty())
        throw SireMol::missing_cutgroup( QObject::tr(
                "There are no CutGroups that contain atoms that match "
                "the Atom ID \"%1\".")
                    .arg(atomid.toString()), CODELOC );
    
    qSort(cgidxs);
    
    return cgidxs;
}

const char* CGsWithAtoms::typeName()
{
    return QMetaType::typeName( qMetaTypeId<CGsWithAtoms>() );
}

CGsWithAtoms* CGsWithAtoms::clone() const
{
    return new CGsWithAtoms(*this);
}

/////////
///////// Implementation of ChainsWithAtoms
/////////

static const RegisterMetaType<ChainsWithAtoms> r_chainid;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const ChainsWithAtoms &chainid)
{
    writeHeader(ds, r_chainid, 1);
    
    ds << chainid.atomid;
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, ChainsWithAtoms &chainid)
{
    VersionID v = readHeader(ds, r_chainid);
    
    if (v == 1)
    {
        ds >> chainid.atomid;
    }
    else
        throw version_error( v, "1", r_chainid, CODELOC );
        
    return ds;
}

/** Null constructor */
ChainsWithAtoms::ChainsWithAtoms() : ChainID()
{}

/** Construct from the passed atom ID */
ChainsWithAtoms::ChainsWithAtoms(const AtomID &id)
               : ChainID(), atomid(id)
{}

/** Copy constructor */
ChainsWithAtoms::ChainsWithAtoms(const ChainsWithAtoms &other)
               : ChainID(other), atomid(other.atomid)
{}

/** Destructor */
ChainsWithAtoms::~ChainsWithAtoms()
{}

/** Is this selection null? */
bool ChainsWithAtoms::isNull() const
{
    return atomid.isNull();
}

/** Return a hash of this identifier */
uint ChainsWithAtoms::hash() const
{
    return atomid.hash();
}
            
/** Return a string representatio of this ID */
QString ChainsWithAtoms::toString() const
{
    return QObject::tr("ChainsWithAtoms( %1 )").arg(atomid.toString());
}

/** Return the atom ID */
const AtomID& ChainsWithAtoms::atomID() const
{
    return atomid.base();
}

/** Copy assignment operator */
ChainsWithAtoms& ChainsWithAtoms::operator=(const ChainsWithAtoms &other)
{
    atomid = other.atomid;
    return *this;
}

/** Comparison operator */
bool ChainsWithAtoms::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare<ChainsWithAtoms>(*this, other);
}

/** Comparison operator */
bool ChainsWithAtoms::operator==(const ChainsWithAtoms &other) const
{
    return atomid == other.atomid;
}

/** Comparison operator */
bool ChainsWithAtoms::operator!=(const ChainsWithAtoms &other) const
{
    return not ChainsWithAtoms::operator==(other);
}

/** Map this ID to the list of indicies of chains that match this ID

    \throw SireMol::missing_chain
    \throw SireError::invalid_index
*/
QList<ChainIdx> ChainsWithAtoms::map(const MolInfo &molinfo) const
{
    QSet<ChainIdx> chainids;
    QList<ChainIdx> chainidxs;
    
    foreach (AtomIdx atomidx, molinfo.map(atomid))
    {
        try
        {
            ChainIdx chainidx = molinfo.parentChain(atomidx);
        
            if (not chainids.contains(chainidx))
            {
                chainidxs.append(chainidx);
                chainids.insert(chainidx);
            }
        }
        catch(...)
        {}
    }
    
    if (chainidxs.isEmpty())
        throw SireMol::missing_chain( QObject::tr(
                "There are no chains that contain atoms that match "
                "the Atom ID \"%1\".")
                    .arg(atomid.toString()), CODELOC );
    
    qSort(chainidxs);
    
    return chainidxs;
}

const char* ChainsWithAtoms::typeName()
{
    return QMetaType::typeName( qMetaTypeId<ChainsWithAtoms>() );
}

ChainsWithAtoms* ChainsWithAtoms::clone() const
{
    return new ChainsWithAtoms(*this);
}

/////////
///////// Implementation of SegsWithAtoms
/////////

static const RegisterMetaType<SegsWithAtoms> r_segid;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const SegsWithAtoms &segid)
{
    writeHeader(ds, r_segid, 1);
    
    ds << segid.atomid;
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, SegsWithAtoms &segid)
{
    VersionID v = readHeader(ds, r_segid);
    
    if (v == 1)
    {
        ds >> segid.atomid;
    }
    else
        throw version_error( v, "1", r_segid, CODELOC );
        
    return ds;
}

/** Null constructor */
SegsWithAtoms::SegsWithAtoms() : SegID()
{}

/** Construct from the passed atom ID */
SegsWithAtoms::SegsWithAtoms(const AtomID &id)
              : SegID(), atomid(id)
{}

/** Copy constructor */
SegsWithAtoms::SegsWithAtoms(const SegsWithAtoms &other)
              : SegID(other), atomid(other.atomid)
{}

/** Destructor */
SegsWithAtoms::~SegsWithAtoms()
{}

/** Is this selection null? */
bool SegsWithAtoms::isNull() const
{
    return atomid.isNull();
}

/** Return a hash of this identifier */
uint SegsWithAtoms::hash() const
{
    return atomid.hash();
}
            
/** Return a string representatio of this ID */
QString SegsWithAtoms::toString() const
{
    return QObject::tr("SegsWithAtoms( %1 )").arg(atomid.toString());
}

/** Return the atom ID */
const AtomID& SegsWithAtoms::atomID() const
{
    return atomid.base();
}

/** Copy assignment operator */
SegsWithAtoms& SegsWithAtoms::operator=(const SegsWithAtoms &other)
{
    atomid = other.atomid;
    return *this;
}

/** Comparison operator */
bool SegsWithAtoms::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare<SegsWithAtoms>(*this, other);
}

/** Comparison operator */
bool SegsWithAtoms::operator==(const SegsWithAtoms &other) const
{
    return atomid == other.atomid;
}

/** Comparison operator */
bool SegsWithAtoms::operator!=(const SegsWithAtoms &other) const
{
    return not SegsWithAtoms::operator==(other);
}

/** Map this ID to the list of indicies of segments that match this ID

    \throw SireMol::missing_segment
    \throw SireError::invalid_index
*/
QList<SegIdx> SegsWithAtoms::map(const MolInfo &molinfo) const
{
    QSet<SegIdx> segids;
    QList<SegIdx> segidxs;
    
    foreach (AtomIdx atomidx, molinfo.map(atomid))
    {
        try
        {
            SegIdx segidx = molinfo.parentSegment(atomidx);
        
            if (not segids.contains(segidx))
            {
                segidxs.append(segidx);
                segids.insert(segidx);
            }
        }
        catch(...)
        {}
    }
    
    if (segidxs.isEmpty())
        throw SireMol::missing_segment( QObject::tr(
                "There are no segments that contain atoms that match "
                "the Atom ID \"%1\".")
                    .arg(atomid.toString()), CODELOC );
    
    qSort(segidxs);
    
    return segidxs;
}

const char* SegsWithAtoms::typeName()
{
    return QMetaType::typeName( qMetaTypeId<SegsWithAtoms>() );
}

SegsWithAtoms* SegsWithAtoms::clone() const
{
    return new SegsWithAtoms(*this);
}
