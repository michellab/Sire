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

#ifndef SIREMOL_GROUPGROUPIDS_H
#define SIREMOL_GROUPGROUPIDS_H

#include "groupatomids.h"

#include "atomid.h"
#include "cgid.h"
#include "resid.h"
#include "chainid.h"
#include "segid.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
template<class G0, class G1>
class GroupGroupID;
}

template<class G0, class G1>
QDataStream& operator<<(QDataStream&, const SireMol::GroupGroupID<G0,G1>&);
template<class G0, class G1>
QDataStream& operator>>(QDataStream&, SireMol::GroupGroupID<G0,G1>&);

namespace SireMol
{

/** This class holds the combination of two different ID types
    for two different groups */
template<class G0, class G1>
class SIREMOL_EXPORT GroupGroupID : public GroupAtomIDBase
{

friend SIREMOL_EXPORT QDataStream& ::operator<<<>(QDataStream&, const GroupGroupID<G0,G1>&);
friend SIREMOL_EXPORT QDataStream& ::operator>><>(QDataStream&, GroupGroupID<G0,G1>&);

public:
    GroupGroupID();
    
    GroupGroupID(const G0 &group0, const G1 &group1);
    
    GroupGroupID(const GroupGroupID<G0,G1> &other);
    
    ~GroupGroupID();
    
    static const char* typeName();
    
    const char* what() const;
       
    GroupGroupID<G0,G1>* clone() const;
    
    using SireID::ID::operator!=;
    
    bool operator==(const SireID::ID &other) const;
    
    bool operator==(const GroupGroupID<G0,G1> &other) const;
    
    bool operator!=(const GroupGroupID<G0,G1> &other) const;
    
    QString toString() const;
    
    uint hash() const;
    
    bool isNull() const;
             
    QList<AtomIdx> map(const MolInfo &molinfo) const;

private:
    typename G0::Identifier g0;
    typename G1::Identifier g1;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Constructor */
template<class G0, class G1>
SIRE_OUTOFLINE_TEMPLATE
GroupGroupID<G0,G1>::GroupGroupID() : GroupAtomIDBase()
{}

/** Construct the combination of the two groups */
template<class G0, class G1>
SIRE_OUTOFLINE_TEMPLATE
GroupGroupID<G0,G1>::GroupGroupID(const G0 &group0, const G1 &group1)
                    : GroupAtomIDBase(), g0(group0), g1(group1)
{}

/** Copy constructor */
template<class G0, class G1>
SIRE_OUTOFLINE_TEMPLATE
GroupGroupID<G0,G1>::GroupGroupID(const GroupGroupID<G0,G1> &other)
                    : GroupAtomIDBase(),
                      g0(other.g0), g1(other.g1)
{}

/** Destructor */
template<class G0, class G1>
SIRE_OUTOFLINE_TEMPLATE
GroupGroupID<G0,G1>::~GroupGroupID()
{}
    
template<class G0, class G1>
SIRE_OUTOFLINE_TEMPLATE
const char* GroupGroupID<G0,G1>::typeName()
{
    return QMetaType::typeName( qMetaTypeId< GroupGroupID<G0,G1> >() );
}

template<class G0, class G1>
SIRE_OUTOFLINE_TEMPLATE
const char* GroupGroupID<G0,G1>::what() const
{
    return GroupGroupID<G0,G1>::typeName();
}
   
template<class G0, class G1>
SIRE_OUTOFLINE_TEMPLATE
GroupGroupID<G0,G1>* GroupGroupID<G0,G1>::clone() const
{
    return new GroupGroupID<G0,G1>(*this);
}

template<class G0, class G1>
SIRE_OUTOFLINE_TEMPLATE
bool GroupGroupID<G0,G1>::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare< GroupGroupID<G0,G1> >(*this, other);
}

template<class G0, class G1>
SIRE_OUTOFLINE_TEMPLATE
bool GroupGroupID<G0,G1>::operator==(const GroupGroupID<G0,G1> &other) const
{
    return g0 == other.g0 and g1 == other.g1;
}

template<class G0, class G1>
SIRE_OUTOFLINE_TEMPLATE
bool GroupGroupID<G0,G1>::operator!=(const GroupGroupID<G0,G1> &other) const
{
    return g0 != other.g0 or g1 != other.g1;
}

template<class G0, class G1>
SIRE_OUTOFLINE_TEMPLATE
uint GroupGroupID<G0,G1>::hash() const
{
    return (g0.hash() << 16) | (g1.hash() & 0x0000FFFF);
}

template<class G0, class G1>
SIRE_OUTOFLINE_TEMPLATE
bool GroupGroupID<G0,G1>::isNull() const
{
    return g0.isNull() and g1.isNull();
}

/** Return a string representation of this match */
template<class G0, class G1>
SIRE_OUTOFLINE_TEMPLATE
QString GroupGroupID<G0,G1>::toString() const
{
    if (this->isNull())
        return "*";
    else if (g0.isNull())
        return g1.toString();
    else if (g1.isNull())
        return g0.toString();
    else
        return QString("%1 and %2").arg(g0.toString(), g1.toString());
}

/** Map this pair of group IDs to the list of indicies of atoms that are
    in both of these groups.
  
    \throw ???::missing_G0
    \throw ???::missing_G1
    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
template<class G0, class G1>
SIRE_OUTOFLINE_TEMPLATE
QList<AtomIdx> GroupGroupID<G0,G1>::map(const MolInfo &molinfo) const
{
    if (this->isNull())
        return molinfo.getAtoms();
    else if (g0.isNull())
        return molinfo.getAtomsIn(g1);
    else if (g1.isNull())
        return molinfo.getAtomsIn(g0);
        
    QList<AtomIdx> atomidxs = 
                     MolInfo::intersection(molinfo.getAtomsIn(g0),
                                           molinfo.getAtomsIn(g1) );
        
    if (atomidxs.isEmpty())
        this->throwMissingAtom(molinfo);
            
    qSort(atomidxs);
    return atomidxs;
}

#endif

typedef GroupGroupID<SegID,ResID> SegResID;
typedef GroupGroupID<SegID,ChainID> SegChainID;
typedef GroupGroupID<SegID,CGID> SegCGID;

typedef GroupGroupID<CGID,ResID> CGResID;
typedef GroupGroupID<CGID,ChainID> CGChainID;

}

#ifndef SIRE_SKIP_INLINE_FUNCTIONS
/** Serialise to a binary datastream */
template<class G0, class G1>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator<<(QDataStream &ds, const SireMol::GroupGroupID<G0,G1> &groupid)
{
    ds << groupid.g0 << groupid.g1;
    return ds;
}

/** Extract from a binary datastream */
template<class G0, class G1>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator>>(QDataStream &ds, SireMol::GroupGroupID<G0,G1> &groupid)
{
    ds >> groupid.g0 >> groupid.g1;
    return ds;
}
#endif // SIRE_SKIP_INLINE_FUNCTIONS

Q_DECLARE_METATYPE(SireMol::SegResID);
Q_DECLARE_METATYPE(SireMol::SegChainID);
Q_DECLARE_METATYPE(SireMol::SegCGID);
Q_DECLARE_METATYPE(SireMol::CGResID);
Q_DECLARE_METATYPE(SireMol::CGChainID);

SIRE_EXPOSE_ALIAS( (SireMol::GroupGroupID<SireMol::SegID, SireMol::ResID>), 
                    SireMol::SegResID )
SIRE_EXPOSE_ALIAS( (SireMol::GroupGroupID<SireMol::SegID, SireMol::ChainID>), 
                    SireMol::SegChainID )
SIRE_EXPOSE_ALIAS( (SireMol::GroupGroupID<SireMol::SegID, SireMol::CGID>), 
                    SireMol::SegCGID )
SIRE_EXPOSE_ALIAS( (SireMol::GroupGroupID<SireMol::CGID, SireMol::ResID>), 
                    SireMol::CGResID )
SIRE_EXPOSE_ALIAS( (SireMol::GroupGroupID<SireMol::CGID, SireMol::ChainID>), 
                    SireMol::CGChainID )

#ifdef SIRE_INSTANTIATE_TEMPLATES
template class SireMol::GroupGroupID<SireMol::SegID,SireMol::ResID>;
template class SireMol::GroupGroupID<SireMol::SegID,SireMol::ChainID>;
template class SireMol::GroupGroupID<SireMol::SegID,SireMol::CGID>;
template class SireMol::GroupGroupID<SireMol::CGID,SireMol::ResID>;
template class SireMol::GroupGroupID<SireMol::CGID,SireMol::ChainID>;
#endif

SIRE_END_HEADER

#endif
