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

#include "withres.h"

#include "SireMol/errors.h"

#include "SireStream/datastream.h"

using namespace SireMol;
using namespace SireID;
using namespace SireStream;

static const RegisterMetaType<ChainsWithRes> r_chainid;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const ChainsWithRes &chainid)
{
    writeHeader(ds, r_chainid, 1);
    
    ds << chainid.resid;
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, ChainsWithRes &chainid)
{
    VersionID v = readHeader(ds, r_chainid);
    
    if (v == 1)
    {
        ds >> chainid.resid;
    }
    else
        throw version_error( v, "1", r_chainid, CODELOC );
        
    return ds;
}

/** Null constructor */
ChainsWithRes::ChainsWithRes() : ChainID()
{}

/** Construct from the passed residue ID */
ChainsWithRes::ChainsWithRes(const ResID &id)
               : ChainID(), resid(id)
{}

/** Copy constructor */
ChainsWithRes::ChainsWithRes(const ChainsWithRes &other)
               : ChainID(other), resid(other.resid)
{}

/** Destructor */
ChainsWithRes::~ChainsWithRes()
{}

/** Is this selection null? */
bool ChainsWithRes::isNull() const
{
    return resid.isNull();
}

/** Return a hash of this identifier */
uint ChainsWithRes::hash() const
{
    return resid.hash();
}
            
/** Return a string representatio of this ID */
QString ChainsWithRes::toString() const
{
    return QObject::tr("ChainsWithRes( %1 )").arg(resid.toString());
}

/** Return the residue ID */
const ResID& ChainsWithRes::resID() const
{
    return resid.base();
}

/** Copy assignment operator */
ChainsWithRes& ChainsWithRes::operator=(const ChainsWithRes &other)
{
    resid = other.resid;
    return *this;
}

/** Comparison operator */
bool ChainsWithRes::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare<ChainsWithRes>(*this, other);
}

/** Comparison operator */
bool ChainsWithRes::operator==(const ChainsWithRes &other) const
{
    return resid == other.resid;
}

/** Comparison operator */
bool ChainsWithRes::operator!=(const ChainsWithRes &other) const
{
    return not ChainsWithRes::operator==(other);
}

/** Map this ID to the list of indicies of chains that match this ID

    \throw SireMol::missing_chain
    \throw SireError::invalid_index
*/
QList<ChainIdx> ChainsWithRes::map(const MolInfo &molinfo) const
{
    QSet<ChainIdx> chainids;
    QList<ChainIdx> chainidxs;
    
    foreach (ResIdx residx, molinfo.map(resid))
    {
        try
        {
            ChainIdx chainidx = molinfo.parentChain(residx);
        
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
                "There are no chains that contain residues that match "
                "the residue ID \"%1\".")
                    .arg(resid.toString()), CODELOC );
    
    qSort(chainidxs);
    
    return chainidxs;
}

const char* ChainsWithRes::typeName()
{
    return QMetaType::typeName( qMetaTypeId<ChainsWithRes>() );
}

ChainsWithRes* ChainsWithRes::clone() const
{
    return new ChainsWithRes(*this);
}
