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

#ifndef SIREMOL_RESIN_HPP
#define SIREMOL_RESIN_HPP

#include "SireID/index.h"

#include "residx.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
template<class GROUP>
class ResIn;
}

template<class GROUP>
QDataStream& operator<<(QDataStream&, const SireMol::ResIn<GROUP>&);
template<class GROUP>
QDataStream& operator>>(QDataStream&, SireMol::ResIn<GROUP>&);

namespace SireMol
{

class MolInfo;

/** This helper class is used to provide the '.residues()' functionality
    of the group ID classes. This allows the class to a residue, or
    range of residues by index from the group that has been identified.
    
    @author Christopher Woods
*/
template<class GROUP>
class SIREMOL_EXPORT ResIn : public ResID
{

friend SIREMOL_EXPORT QDataStream& ::operator<<<>(QDataStream&, const ResIn<GROUP>&);
friend SIREMOL_EXPORT QDataStream& ::operator>><>(QDataStream&, ResIn<GROUP>&);

public:
    ResIn();
    
    ResIn(const GROUP &id);
    ResIn(const GROUP &id, qint32 i);
    ResIn(const GROUP &id, qint32 i, qint32 j);
    
    ResIn(const ResIn<GROUP> &other);
    
    ~ResIn();
    
    static const char* typeName();
    
    const char* what() const;
    
    ResIn<GROUP>* clone() const;
    
    ResIn<GROUP>& operator=(const ResIn<GROUP> &other);

    bool operator==(const ResIn<GROUP> &other) const;
    bool operator==(const SireID::ID &other) const;
    
    bool operator!=(const ResIn<GROUP> &other) const;
    
    bool operator!=(const SireID::ID &other) const;
    
    bool isNull() const;
    
    uint hash() const;
    
    QString toString() const;
    
    QList<ResIdx> map(const MolInfo &molinfo) const;

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
ResIn<GROUP>::ResIn()
             : ResID(), strt(0), end(-1)
{}

/** Construct to get all of the residues in the passed group */
template<class GROUP>
SIRE_OUTOFLINE_TEMPLATE
ResIn<GROUP>::ResIn(const GROUP &id)
             : ResID(), groupid(id), strt(0), end(-1)
{}

/** Construct for a specified residue */
template<class GROUP>
SIRE_OUTOFLINE_TEMPLATE
ResIn<GROUP>::ResIn(const GROUP &id, qint32 i)
             : ResID(), groupid(id), strt(i), end(i)
{}

/** Construct for a range of residues */
template<class GROUP>
SIRE_OUTOFLINE_TEMPLATE
ResIn<GROUP>::ResIn(const GROUP &id, qint32 i, qint32 j)
             : ResID(), groupid(id), strt(i), end(j)
{}

/** Copy constructor */
template<class GROUP>
SIRE_OUTOFLINE_TEMPLATE
ResIn<GROUP>::ResIn(const ResIn<GROUP> &other)
            : ResID(other), groupid(other.groupid), 
              strt(other.strt), end(other.end)
{}

/** Destructor */
template<class GROUP>
SIRE_OUTOFLINE_TEMPLATE
ResIn<GROUP>::~ResIn()
{}
    
template<class GROUP>
SIRE_OUTOFLINE_TEMPLATE
const char* ResIn<GROUP>::typeName()
{
    return QMetaType::typeName( qMetaTypeId< ResIn<GROUP> >() );
}

template<class GROUP>
SIRE_OUTOFLINE_TEMPLATE
const char* ResIn<GROUP>::what() const
{
    return ResIn<GROUP>::typeName();
}

template<class GROUP>
SIRE_OUTOFLINE_TEMPLATE
ResIn<GROUP>* ResIn<GROUP>::clone() const
{
    return new ResIn<GROUP>(*this);
}

template<class GROUP>
SIRE_OUTOFLINE_TEMPLATE
bool ResIn<GROUP>::operator!=(const ResIn<GROUP> &other) const
{
    return not this->operator==(other);
}

template<class GROUP>
SIRE_OUTOFLINE_TEMPLATE
bool ResIn<GROUP>::operator!=(const SireID::ID &other) const
{
    return not this->operator==(other);
}

template<class GROUP>
SIRE_OUTOFLINE_TEMPLATE
bool ResIn<GROUP>::isNull() const
{
    return groupid.isNull() and strt == 0 and end == -1;
}

template<class GROUP>
SIRE_OUTOFLINE_TEMPLATE
uint ResIn<GROUP>::hash() const
{
    return groupid.hash() + strt + end;
}

/** Return a string representation of this ID */
template<class GROUP>
SIRE_OUTOFLINE_TEMPLATE
QString ResIn<GROUP>::toString() const
{
    if (strt == 0 and end == -1)
        return QString("(%1).residues()").arg(groupid.toString());
    else if (strt == end)
        return QString("(%1).residue(%2)").arg(groupid.toString()).arg(strt);
    else
        return QString("(%1).residues(%2,%3)")
                    .arg(groupid.toString()).arg(strt).arg(end);
}

/** Copy assignment operator */
template<class GROUP>
SIRE_OUTOFLINE_TEMPLATE
ResIn<GROUP>& ResIn<GROUP>::operator=(const ResIn<GROUP> &other)
{
    if (this != &other)
    {
        ResID::operator=(other);
        groupid = other.groupid;
        strt = other.strt;
        end = other.end;
    }

    return *this;
}

/** Comparison operator */
template<class GROUP>
SIRE_OUTOFLINE_TEMPLATE
bool ResIn<GROUP>::operator==(const ResIn<GROUP> &other) const
{
    return strt == other.strt and end == other.end and
           groupid == other.groupid;
}

/** Comparison operator */
template<class GROUP>
SIRE_OUTOFLINE_TEMPLATE
bool ResIn<GROUP>::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare< ResIn<GROUP> >(*this, other);
}

/** Map this ID to the indicies of matching residues
    
    \throw ???::missing_ID
    \throw SireError::invalid_index
*/
template<class GROUP>
SIRE_OUTOFLINE_TEMPLATE
QList<ResIdx> ResIn<GROUP>::map(const MolInfo &molinfo) const
{
    //first get the list of the indicies of the matching groups
    QList<typename GROUP::Index> idxs = groupid.map(molinfo);
    
    //now get a list of the indicies of all of the residues in these groups
    QList<ResIdx> residxs;
    
    foreach (typename GROUP::Index idx, idxs)
    {
        residxs += molinfo.getResiduesIn(idx);
    }
    
    //now map _i and _j to the indicies...
    int nres = residxs.count();
    
    int sane_strt = strt.map(nres);
    int sane_end = end.map(nres);
    
    if (sane_strt > sane_end)
        qSwap(sane_strt, sane_end);
    
    //now extract only the desired atom indicies
    if (sane_end - sane_strt == nres)
    {
        return residxs;
    }
    else
    {
        QList<ResIdx> specified_residxs;
    
        for (int i=sane_strt; i<=sane_end; ++i)
        {
            specified_residxs.append( residxs[i] );
        }
    
        return specified_residxs;
    }
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

#ifndef SIRE_SKIP_INLINE_FUNCTIONS
template<class GROUP>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator<<(QDataStream &ds, const SireMol::ResIn<GROUP> &res)
{
    ds << res.groupid << res.strt << res.end;
    return ds;
}

template<class GROUP>
SIRE_OUTOFLINE_TEMPLATE
QDataStream& operator>>(QDataStream &ds, SireMol::ResIn<GROUP> &res)
{
    ds >> res.groupid >> res.strt >> res.end;
    return ds;
}
#endif // SIRE_SKIP_INLINE_FUNCTIONS

SIRE_END_HEADER

#endif
