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

#ifndef SIREMOL_MOLINFO_H
#define SIREMOL_MOLINFO_H

#include "sireglobal.h"

#include <QList>

SIRE_BEGIN_HEADER

namespace SireMol
{

class AtomIdx;
class AtomName;
class AtomNum;
class AtomID;

class ResIdx;
class ResName;
class ResNum;
class ResID;

class CGIdx;
class CGName;
class CGID;

class ChainIdx;
class ChainName;
class ChainID;

class SegIdx;
class SegName;
class SegID;

class AtomSelection;

/** This is the pure virtual interface to the metadata object
    that can be used to track part of a molecule back to its
    source based on an ID, e.g. being able to find the index
    of an atom based on its name and the name of residue
    that contains it.

    @author Christopher Woods
*/
class SIREMOL_EXPORT MolInfo
{
public:
    MolInfo();
    virtual ~MolInfo();

    virtual QList<AtomIdx> map(const AtomName &name) const=0;
    virtual QList<AtomIdx> map(AtomNum num) const=0;
    virtual QList<AtomIdx> map(AtomIdx idx) const=0;
    virtual QList<AtomIdx> map(const AtomID &atomid) const=0;

    virtual QList<ResIdx> map(const ResName &name) const=0;
    virtual QList<ResIdx> map(ResNum num) const=0;
    virtual QList<ResIdx> map(ResIdx idx) const=0;
    virtual QList<ResIdx> map(const ResID &resid) const=0;

    virtual QList<CGIdx> map(const CGName &name) const=0;
    virtual QList<CGIdx> map(CGIdx idx) const=0;
    virtual QList<CGIdx> map(const CGID &cgid) const=0;

    virtual QList<ChainIdx> map(const ChainName &name) const=0;
    virtual QList<ChainIdx> map(ChainIdx idx) const=0;
    virtual QList<ChainIdx> map(const ChainID &chainid) const=0;

    virtual QList<SegIdx> map(const SegName &name) const=0;
    virtual QList<SegIdx> map(SegIdx idx) const=0;
    virtual QList<SegIdx> map(const SegID &segid) const=0;

    virtual QList<AtomIdx> getAtoms() const=0;

    virtual AtomIdx getAtom(CGIdx cgidx, int i) const=0;
    virtual AtomIdx getAtom(ResIdx residx, int i) const=0;
    virtual AtomIdx getAtom(ChainIdx chainidx, int i) const=0;
    virtual AtomIdx getAtom(SegIdx segidx, int i) const=0;

    virtual QList<AtomIdx> getAtomsIn(const ResID &resid) const=0;
    virtual QList<AtomIdx> getAtomsIn(const CGID &cgid) const=0;
    virtual QList<AtomIdx> getAtomsIn(const ChainID &chainid) const=0;
    virtual QList<AtomIdx> getAtomsIn(const SegID &segid) const=0;

    virtual QList<ResIdx> getResidues() const=0;

    virtual ResIdx getResidue(ChainIdx chainidx, int i) const=0;

    virtual QList<ResIdx> getResiduesIn(const ChainID &chainid) const=0;

    virtual QList<CGIdx> getCutGroups() const=0;
    virtual QList<ChainIdx> getChains() const=0;
    virtual QList<SegIdx> getSegments() const=0;

    virtual ChainIdx parentChain(ResIdx residx) const=0;
    virtual ChainIdx parentChain(const ResID &resid) const=0;

    virtual ChainIdx parentChain(AtomIdx atomidx) const=0;
    virtual ChainIdx parentChain(const AtomID &atomid) const=0;

    virtual ResIdx parentResidue(AtomIdx atomidx) const=0;
    virtual ResIdx parentResidue(const AtomID &atomid) const=0;

    virtual SegIdx parentSegment(AtomIdx atomidx) const=0;
    virtual SegIdx parentSegment(const AtomID &atomid) const=0;

    virtual CGIdx parentCutGroup(AtomIdx atomidx) const=0;
    virtual CGIdx parentCutGroup(const AtomID &atomid) const=0;

    virtual AtomIdx atomIdx(const AtomID &atomid) const=0;
    virtual CGIdx cgIdx(const CGID &cgid) const=0;
    virtual ResIdx resIdx(const ResID &resid) const=0;
    virtual ChainIdx chainIdx(const ChainID &chainid) const=0;
    virtual SegIdx segIdx(const SegID &segid) const=0;

    virtual void assertCompatibleWith(const AtomSelection &selected_atoms) const=0;

    void assertSingleAtom(const QList<AtomIdx> &atomidxs) const;
    void assertSingleResidue(const QList<ResIdx> &residxs) const;
    void assertSingleChain(const QList<ChainIdx> &chainidxs) const;
    void assertSingleCutGroup(const QList<CGIdx> &cgidxs) const;
    void assertSingleSegment(const QList<SegIdx> &segidxs) const;

    template<class T>
    static QList<T> intersection(const QList<T> &list0, const QList<T> &list1);

};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Return the intersection of the objects in 'list0' and 'list1'. This returns
    a list that contains only those items that are in both lists. The order
    of the list is the same as the input, namely if both list0 and list1 are
    sorted, then the returned list will be sorted as well. */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QList<T> MolInfo::intersection(const QList<T> &list0, const QList<T> &list1)
{
    QList<T> intersection_list;

    if (list0.count() <= list1.count())
    {
        QSet<T> set1(list1.constBegin(), list1.constEnd());

        for (typename QList<T>::const_iterator it = list0.constBegin();
             it != list0.constEnd();
             ++it)
        {
            if (set1.contains(*it))
                intersection_list.append(*it);
        }
    }
    else
    {
        QSet<T> set0(list0.constBegin(), list0.constEnd());

        for (typename QList<T>::const_iterator it = list1.constBegin();
             it != list1.constEnd();
             ++it)
        {
            if (set0.contains(*it))
                intersection_list.append(*it);
        }
    }

    return intersection_list;
}

namespace detail
{

template<class T>
QList<typename T::Index> getAll(const MolInfo &molinfo);

template<class T>
QList<typename T::Index> getAll(const MolInfo &molinfo,
                                const AtomSelection &selected_atoms);

} //end of namespace detail

#endif //SIRE_SKIP_INLINE_FUNCTIONS

SIRE_EXPOSE_CLASS( SireMol::MolInfo )

}

SIRE_END_HEADER

#endif
