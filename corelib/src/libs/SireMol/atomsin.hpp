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

#ifndef SIREMOL_ATOMSIN_HPP
#define SIREMOL_ATOMSIN_HPP

#include "SireID/index.h"

#include "atomidx.h"
#include "atomidentifier.h"
#include "molinfo.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
template<class GROUP>
class AtomsIn;
}

template<class GROUP>
QDataStream& operator<<(QDataStream&, const SireMol::AtomsIn<GROUP>&);
template<class GROUP>
QDataStream& operator>>(QDataStream&, SireMol::AtomsIn<GROUP>&);

namespace SireMol
{

class MoleculeInfoData;

/** This helper class is used to provide the '.atoms()' functionality
    of the group ID classes. This allows the class to an atom, or
    range of atoms by index from the group that has been identified.
    
    @author Christopher Woods
*/
template<class GROUP>
class SIREMOL_EXPORT AtomsIn : public AtomID
{

friend SIREMOL_EXPORT QDataStream& ::operator<<<>(QDataStream&, const AtomsIn<GROUP>&);
friend SIREMOL_EXPORT QDataStream& ::operator>><>(QDataStream&, AtomsIn<GROUP>&);

public:
    AtomsIn();
    
    AtomsIn(const GROUP &id);
    AtomsIn(const GROUP &id, qint32 i);
    AtomsIn(const GROUP &id, qint32 i, qint32 j);
    
    AtomsIn(const AtomsIn<GROUP> &other);
    
    ~AtomsIn();
    
    static const char* typeName();
    
    const char* what() const;
    
    AtomsIn<GROUP>* clone() const;
    
    AtomsIn<GROUP>& operator=(const AtomsIn<GROUP> &other);

    bool operator==(const AtomsIn<GROUP> &other) const;
    bool operator==(const SireID::ID &other) const;
    
    bool operator!=(const AtomsIn<GROUP> &other) const;
    bool operator!=(const SireID::ID &other) const;
    
    bool isNull() const;
    
    uint hash() const;
    
    QString toString() const;
    
    QList<AtomIdx> map(const MolInfo &molinfo) const;

private:
    /** The ID of the group that contains the atoms */
    typename GROUP::Identifier groupid;

    /** The indicies of the range of atoms that this specifies */
    SireID::Index strt,end;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Null constructor */
template<class GROUP>
SIRE_OUTOFLINE_TEMPLATE
AtomsIn<GROUP>::AtomsIn()
               : AtomID(), strt(0), end(-1)
{}

/** Construct to get all of the atoms in the passed group */
template<class GROUP>
SIRE_OUTOFLINE_TEMPLATE
AtomsIn<GROUP>::AtomsIn(const GROUP &id)
               : AtomID(), groupid(id), strt(0), end(-1)
{}

/** Construct for a specified atom */
template<class GROUP>
SIRE_OUTOFLINE_TEMPLATE
AtomsIn<GROUP>::AtomsIn(const GROUP &id, qint32 i)
            : AtomID(), groupid(id), strt(i), end(i)
{}

/** Construct for a range of atoms */
template<class GROUP>
SIRE_OUTOFLINE_TEMPLATE
AtomsIn<GROUP>::AtomsIn(const GROUP &id, qint32 i, qint32 j)
            : AtomID(), groupid(id), strt(i), end(j)
{}

/** Copy constructor */
template<class GROUP>
SIRE_OUTOFLINE_TEMPLATE
AtomsIn<GROUP>::AtomsIn(const AtomsIn<GROUP> &other)
            : AtomID(other), groupid(other.groupid), 
              strt(other.strt), end(other.end)
{}

/** Destructor */
template<class GROUP>
SIRE_OUTOFLINE_TEMPLATE
AtomsIn<GROUP>::~AtomsIn()
{}

/** Return a string representation of this ID */
template<class GROUP>
SIRE_OUTOFLINE_TEMPLATE
QString AtomsIn<GROUP>::toString() const
{
    if (strt == 0 and end == -1)
        return QString("(%1).atoms()").arg(groupid.toString());
    else if (strt == end)
        return QString("(%1).atom(%2)").arg(groupid.toString()).arg(strt);
    else
        return QString("(%1).atoms(%2,%3)")
                    .arg(groupid.toString()).arg(strt).arg(end);
}

/** Copy assignment operator */
template<class GROUP>
SIRE_OUTOFLINE_TEMPLATE
AtomsIn<GROUP>& AtomsIn<GROUP>::operator=(const AtomsIn<GROUP> &other)
{
    if (this != &other)
    {
        AtomID::operator=(other);
        groupid = other.groupid;
        strt = other.strt;
        end = other.end;
    }

    return *this;
}

/** Comparison operator */
template<class GROUP>
SIRE_OUTOFLINE_TEMPLATE
bool AtomsIn<GROUP>::operator==(const AtomsIn<GROUP> &other) const
{
    return strt == other.strt and end == other.end and
           groupid == other.groupid;
}

/** Comparison operator */
template<class GROUP>
SIRE_OUTOFLINE_TEMPLATE
bool AtomsIn<GROUP>::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare< AtomsIn<GROUP> >(*this, other);
}
    
template<class GROUP>
SIRE_OUTOFLINE_TEMPLATE
const char* AtomsIn<GROUP>::typeName()
{
    return QMetaType::typeName( qMetaTypeId< AtomsIn<GROUP> >() );
}

template<class GROUP>
SIRE_OUTOFLINE_TEMPLATE
const char* AtomsIn<GROUP>::what() const
{
    return AtomsIn<GROUP>::typeName();
}

template<class GROUP>
SIRE_OUTOFLINE_TEMPLATE
AtomsIn<GROUP>* AtomsIn<GROUP>::clone() const
{
    return new AtomsIn<GROUP>(*this);
}

template<class GROUP>
SIRE_OUTOFLINE_TEMPLATE
bool AtomsIn<GROUP>::operator!=(const AtomsIn<GROUP> &other) const
{
    return not this->operator==(other);
}

template<class GROUP>
SIRE_OUTOFLINE_TEMPLATE
bool AtomsIn<GROUP>::operator!=(const SireID::ID &other) const
{
    return not this->operator==(other);
}

template<class GROUP>
SIRE_OUTOFLINE_TEMPLATE
bool AtomsIn<GROUP>::isNull() const
{
    return groupid.isNull() and strt == 0 and end == -1;
}

template<class GROUP>
SIRE_OUTOFLINE_TEMPLATE
uint AtomsIn<GROUP>::hash() const
{
    return groupid.hash() + strt + end;
}

/** Map this ID to the indicies of matching atoms 

    \throw ???::missing_ID
    \throw SireError::invalid_index
*/
template<class GROUP>
SIRE_OUTOFLINE_TEMPLATE
QList<AtomIdx> AtomsIn<GROUP>::map(const MolInfo &molinfo) const
{
    //first get the list of the indicies of the matching groups
    QList<typename GROUP::Index> idxs = groupid.map(molinfo);
    
    //now get a list of the indicies of all of the atoms in these groups
    QList<AtomIdx> atomidxs;
    
    foreach (typename GROUP::Index idx, idxs)
    {
        atomidxs += molinfo.getAtomsIn(idx);
    }
    
    //now map _i and _j to the indicies...
    int nats = atomidxs.count();
    
    int sane_strt = strt.map(nats);
    int sane_end = end.map(nats);
    
    if (sane_strt > sane_end)
        qSwap(sane_strt, sane_end);
    
    //now extract only the desired atom indicies
    if (sane_end - sane_strt == nats)
    {
        return atomidxs;
    }
    else
    {
        QList<AtomIdx> specified_atomidxs;
    
        for (int i=sane_strt; i<=sane_end; ++i)
        {
            specified_atomidxs.append( atomidxs[i] );
        }
    
        return specified_atomidxs;
    }
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

/** Serialise to a binary datastream */
template<class GROUP>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator<<(QDataStream &ds, const SireMol::AtomsIn<GROUP> &atoms)
{
    ds << atoms.groupid << atoms.strt << atoms.end;
    return ds;
}

/** Extract from a binary datastream */
template<class GROUP>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator>>(QDataStream &ds, SireMol::AtomsIn<GROUP> &atoms)
{
    ds >> atoms.groupid >> atoms.strt >> atoms.end;
    return ds;
}

SIRE_END_HEADER

#endif
