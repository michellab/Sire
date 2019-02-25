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

#ifndef SIREMOL_GROUPATOMIDS_H
#define SIREMOL_GROUPATOMIDS_H

#include <QSet>

#include "molinfo.h"
#include "atomidx.h"

#include "atomid.h"
#include "cgid.h"
#include "resid.h"
#include "chainid.h"
#include "segid.h"

#include "atomidentifier.h"
#include "chainidentifier.h"
#include "cgidentifier.h"
#include "residentifier.h"
#include "segidentifier.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
template<class GROUP,class ATOM>
class GroupAtomID;
}

template<class GROUP,class ATOM>
QDataStream& operator<<(QDataStream&, const SireMol::GroupAtomID<GROUP,ATOM>&);
template<class GROUP,class ATOM>
QDataStream& operator>>(QDataStream&, SireMol::GroupAtomID<GROUP,ATOM>&);

namespace SireMol
{

class ResID;
class ChainID;
class SegID;
class CGID;

/** This is the base class of GroupAtomID, used to abstract   
    template-independent parts away from the template code.
    
    @author Christopher Woods
*/
class SIREMOL_EXPORT GroupAtomIDBase : public AtomID
{
public:
    GroupAtomIDBase();
    
    GroupAtomIDBase(const GroupAtomIDBase &other);
    
    ~GroupAtomIDBase();
    
protected:
    void throwMissingAtom(const MolInfo &molinfo) const;
};

/** This class represents an Atom ID that is comprised of both 
    an AtomID part and an ID of a group in the molecule (e.g.
    a Residue, Segment, Chain or CutGroup)
    
    @author Christopher Woods
*/
template<class GROUP, class ATOM>
class SIREMOL_EXPORT GroupAtomID : public GroupAtomIDBase
{

friend SIREMOL_EXPORT QDataStream& ::operator<<<>(QDataStream&, const GroupAtomID<GROUP,ATOM>&);
friend SIREMOL_EXPORT QDataStream& ::operator>><>(QDataStream&, GroupAtomID<GROUP,ATOM>&);

public:
    GroupAtomID();
    
    GroupAtomID(const GROUP &groupid, const ATOM &atomid);
    
    GroupAtomID(const GroupAtomID &other);
    
    ~GroupAtomID();
    
    static const char* typeName();
    
    const char* what() const;
    
    GroupAtomID<GROUP,ATOM>* clone() const;
    
    bool operator==(const GroupAtomID<GROUP,ATOM> &other) const;
    
    bool operator!=(const GroupAtomID<GROUP,ATOM> &other) const;
    
    bool operator==(const SireID::ID &other) const;
    
    uint hash() const;
    
    bool isNull() const;
    
    QString toString() const;
    
    QList<AtomIdx> map(const MolInfo &molinfo) const;
 
private:
    typename GROUP::Identifier groupid;
    typename ATOM::Identifier atomid;
};

/** Null constructor */
template<class GROUP, class ATOM>
SIRE_OUTOFLINE_TEMPLATE
GroupAtomID<GROUP,ATOM>::GroupAtomID()
                        : GroupAtomIDBase()
{}

/** Construct using the passed group and atom IDs */
template<class GROUP, class ATOM>
SIRE_OUTOFLINE_TEMPLATE
GroupAtomID<GROUP,ATOM>::GroupAtomID(const GROUP &group,
                                     const ATOM &atom)
                        : GroupAtomIDBase(),
                          groupid(group), atomid(atom)
{}

/** Copy constructor */
template<class GROUP, class ATOM>
SIRE_OUTOFLINE_TEMPLATE
GroupAtomID<GROUP,ATOM>::GroupAtomID(const GroupAtomID<GROUP,ATOM> &other)
                        : GroupAtomIDBase(other),
                          groupid(other.groupid), atomid(other.atomid)
{}

/** Destructor */
template<class GROUP, class ATOM>
SIRE_OUTOFLINE_TEMPLATE
GroupAtomID<GROUP,ATOM>::~GroupAtomID()
{}
    
template<class GROUP, class ATOM>
SIRE_OUTOFLINE_TEMPLATE
const char* GroupAtomID<GROUP,ATOM>::typeName()
{
    return QMetaType::typeName( qMetaTypeId< GroupAtomID<GROUP,ATOM> >() );
}

template<class GROUP, class ATOM>
SIRE_OUTOFLINE_TEMPLATE
const char* GroupAtomID<GROUP,ATOM>::what() const
{
    return GroupAtomID<GROUP,ATOM>::typeName();
}

template<class GROUP, class ATOM>
SIRE_OUTOFLINE_TEMPLATE
GroupAtomID<GROUP,ATOM>* GroupAtomID<GROUP,ATOM>::clone() const
{
    return new GroupAtomID<GROUP,ATOM>(*this);
}

template<class GROUP, class ATOM>
SIRE_OUTOFLINE_TEMPLATE
bool GroupAtomID<GROUP,ATOM>::operator==(const GroupAtomID<GROUP,ATOM> &other) const
{
    return groupid == other.groupid and
           atomid == other.atomid;
}

template<class GROUP, class ATOM>
SIRE_OUTOFLINE_TEMPLATE
bool GroupAtomID<GROUP,ATOM>::operator!=(const GroupAtomID<GROUP,ATOM> &other) const
{
    return not this->operator==(other);
}

template<class GROUP, class ATOM>
SIRE_OUTOFLINE_TEMPLATE
bool GroupAtomID<GROUP,ATOM>::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare< GroupAtomID<GROUP,ATOM> >(*this, other);
}

template<class GROUP, class ATOM>
SIRE_OUTOFLINE_TEMPLATE
uint GroupAtomID<GROUP,ATOM>::hash() const
{
    return (groupid.hash() << 16) | (atomid.hash() & 0x0000FFFF);
}

template<class GROUP, class ATOM>
SIRE_OUTOFLINE_TEMPLATE
bool GroupAtomID<GROUP,ATOM>::isNull() const
{
    return groupid.isNull() and atomid.isNull();
}

/** Return a string representation of this ID */
template<class GROUP, class ATOM>
SIRE_OUTOFLINE_TEMPLATE
QString GroupAtomID<GROUP,ATOM>::toString() const
{
    return QString("%1 and %2").arg(groupid.toString(), atomid.toString());
}

/** Map this combined ID back to the indicies of the atoms that match this ID 

    \throw ???::missing_GROUP
    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
template<class GROUP, class ATOM>
SIRE_OUTOFLINE_TEMPLATE
QList<AtomIdx> GroupAtomID<GROUP,ATOM>::map(const MolInfo &molinfo) const
{
    if (this->isNull())
        return molinfo.getAtoms();
    else if (atomid.isNull())
        return molinfo.getAtomsIn(groupid);
    else if (groupid.isNull())
        return atomid.map(molinfo);
    
    QList<AtomIdx> atomidxs = 
                MolInfo::intersection(atomid.map(molinfo),
                                      molinfo.getAtomsIn(groupid) );
                                             
    if (atomidxs.isEmpty())
        this->throwMissingAtom(molinfo);
            
    return atomidxs;
}

typedef GroupAtomID<ResID,AtomID> ResAtomID;
typedef GroupAtomID<ChainID,AtomID> ChainAtomID;
typedef GroupAtomID<SegID,AtomID> SegAtomID;
typedef GroupAtomID<CGID,AtomID> CGAtomID;

}

#ifndef SIRE_SKIP_INLINE_FUNCTIONS
/** Serialise to a binary datastream */
template<class GROUP, class ATOM>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator<<(QDataStream &ds, 
                        const SireMol::GroupAtomID<GROUP,ATOM> &groupatomid)
{
    ds << groupatomid.groupid << groupatomid.atomid;
    return ds;
}

/** Extract from a binary datastream */
template<class GROUP, class ATOM>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator>>(QDataStream &ds, 
                        SireMol::GroupAtomID<GROUP,ATOM> &groupatomid)
{
    ds >> groupatomid.groupid >> groupatomid.atomid;
    return ds;
}
#endif // SIRE_SKIP_INLINE_FUNCTIONS

Q_DECLARE_METATYPE(SireMol::ResAtomID);
Q_DECLARE_METATYPE(SireMol::ChainAtomID);
Q_DECLARE_METATYPE(SireMol::SegAtomID);
Q_DECLARE_METATYPE(SireMol::CGAtomID);

SIRE_EXPOSE_CLASS( SireMol::GroupAtomIDBase )
SIRE_EXPOSE_ALIAS( (SireMol::GroupAtomID<SireMol::ResID, SireMol::AtomID>), 
                    SireMol::ResAtomID )
SIRE_EXPOSE_ALIAS( (SireMol::GroupAtomID<SireMol::ChainID, SireMol::AtomID>), 
                    SireMol::ChainAtomID )
SIRE_EXPOSE_ALIAS( (SireMol::GroupAtomID<SireMol::SegID, SireMol::AtomID>), 
                    SireMol::SegAtomID )
SIRE_EXPOSE_ALIAS( (SireMol::GroupAtomID<SireMol::CGID, SireMol::AtomID>), 
                    SireMol::CGAtomID )

#ifdef SIRE_INSTANTIATE_TEMPLATES
template class SireMol::GroupAtomID<SireMol::ResID,SireMol::AtomID>;
template class SireMol::GroupAtomID<SireMol::ChainID,SireMol::AtomID>;
template class SireMol::GroupAtomID<SireMol::SegID,SireMol::AtomID>;
template class SireMol::GroupAtomID<SireMol::CGID,SireMol::AtomID>;
#endif

SIRE_END_HEADER

#endif
