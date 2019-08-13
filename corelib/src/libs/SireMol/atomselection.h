/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2006  Christopher Woods
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

#ifndef SIREMOL_ATOMSELECTION_H
#define SIREMOL_ATOMSELECTION_H

#include <QSet>
#include <QHash>

#include "molviewproperty.h"

#include "SireBase/shareddatapointer.hpp"

#include "SireID/index.h"

#include "cgidx.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class AtomSelection;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::AtomSelection&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::AtomSelection&);

namespace SireMol
{

class Molecule;
class MoleculeView;
class MoleculeInfoData;
class MoleculeData;

class AtomIdx;
class CGIdx;
class ResIdx;
class ChainIdx;
class SegIdx;

class AtomID;
class CGID;
class ResID;
class ChainID;
class SegID;

class CGAtomIdx;

template<class T>
class AtomProperty;

using SireID::Index;

using SireBase::ConcreteProperty;
using SireBase::SharedDataPointer;

/** This class holds information about a selection of atoms in a Molecule.
    The selection is held in the most memory-efficient manner possible,
    and takes advantage of the CutGroup-based layout of Molecule objects.

    This is a const-class, which returns new AtomSelections that
    represent any change.

    @author Christopher Woods
*/
class SIREMOL_EXPORT AtomSelection
           : public ConcreteProperty<AtomSelection,MoleculeProperty>
{

friend SIREMOL_EXPORT QDataStream& ::operator<<(QDataStream&, const AtomSelection&);
friend SIREMOL_EXPORT QDataStream& ::operator>>(QDataStream&, AtomSelection&);

friend class SelectionFromMol; //so can modify a single AtomSelection!

public:
    AtomSelection();

    AtomSelection(const MoleculeView &molecule);
    AtomSelection(const MoleculeData &moldata);
    AtomSelection(const MoleculeInfoData &molinfo);

    AtomSelection(const AtomSelection &other);

    ~AtomSelection();

    static const char* typeName();

    AtomSelection& operator=(const AtomSelection &other);

    bool operator==(const AtomSelection &other) const;
    bool operator!=(const AtomSelection &other) const;

    bool isEmpty() const;
    bool isNull() const;

    const MoleculeInfoData& info() const;

    int nSelected() const;

    int nSelected(CGIdx cgidx) const;
    int nSelected(AtomIdx atomidx) const;
    int nSelected(ResIdx residx) const;
    int nSelected(ChainIdx chainidx) const;
    int nSelected(SegIdx segidx) const;
    
    int nSelected(const CGID &cgid) const;
    int nSelected(const AtomID &atomid) const;
    int nSelected(const ResID &resid) const;
    int nSelected(const ChainID &chainid) const;
    int nSelected(const SegID &segid) const;
    
    int nSelected(const AtomSelection &selection) const;

    int nSelectedAtoms() const;
    int nSelectedCutGroups() const;
    int nSelectedResidues() const;
    int nSelectedChains() const;
    int nSelectedSegments() const;

    int nAtoms() const;
    int nCutGroups() const;
    int nResidues() const;
    int nChains() const;
    int nSegments() const;

    bool selectedAllAtoms() const;
    bool selectedAllCutGroups() const;
    bool selectedAllResidues() const;
    bool selectedAllChains() const;
    bool selectedAllSegments() const;

    bool selected(const CGAtomIdx &cgatomidx) const;
    bool selected(AtomIdx atomidx) const;
    bool selected(const AtomID &atomid) const;

    bool selected(CGIdx cgidx) const;
    bool selected(ResIdx residx) const;
    bool selected(ChainIdx chainidx) const;
    bool selected(SegIdx segidx) const;
    
    bool selected(const CGID &cgid) const;
    bool selected(const ResID &resid) const;
    bool selected(const ChainID &chainid) const;
    bool selected(const SegID &segid) const;

    bool selected(const AtomSelection &selection) const;

    bool selectedAll() const;

    bool selectedAll(AtomIdx atomidx) const;
    bool selectedAll(CGIdx cgidx) const;
    bool selectedAll(ResIdx residx) const;
    bool selectedAll(ChainIdx chainidx) const;
    bool selectedAll(SegIdx segidx) const;
    
    bool selectedAll(const AtomID &atomid) const;
    bool selectedAll(const CGID &cgid) const;
    bool selectedAll(const ResID &resid) const;
    bool selectedAll(const ChainID &chainid) const;
    bool selectedAll(const SegID &segid) const;

    bool selectedAll(const AtomSelection &selection) const;

    bool selectedNone() const;

    bool selectedNone(AtomIdx atomidx) const;
    bool selectedNone(CGIdx cgidx) const;
    bool selectedNone(ResIdx residx) const;
    bool selectedNone(ChainIdx chainidx) const;
    bool selectedNone(SegIdx segidx) const;
    
    bool selectedNone(const AtomID &atomid) const;
    bool selectedNone(const CGID &cgid) const;
    bool selectedNone(const ResID &resid) const;
    bool selectedNone(const ChainID &chainid) const;
    bool selectedNone(const SegID &segid) const;

    bool selectedNone(const AtomSelection &selection) const;

    AtomSelection& selectAll();
    AtomSelection& deselectAll();
    AtomSelection& selectNone();

    AtomSelection& select(AtomIdx atomidx);
    AtomSelection& deselect(AtomIdx atomidx);
    AtomSelection& selectOnly(AtomIdx atomidx);

    AtomSelection& select(CGIdx cgidx);
    AtomSelection& deselect(CGIdx cgidx);
    AtomSelection& selectOnly(CGIdx cgidx);

    AtomSelection& select(ResIdx residx);
    AtomSelection& deselect(ResIdx residx);
    AtomSelection& selectOnly(ResIdx residx);

    AtomSelection& select(ChainIdx chainidx);
    AtomSelection& deselect(ChainIdx chainidx);
    AtomSelection& selectOnly(ChainIdx chainidx);
    
    AtomSelection& select(SegIdx segidx);
    AtomSelection& deselect(SegIdx segidx);
    AtomSelection& selectOnly(SegIdx segidx);

    AtomSelection& select(const QSet<AtomIdx> &atomidxs);
    AtomSelection& deselect(const QSet<AtomIdx> &atomidxs);
    AtomSelection& selectOnly(const QSet<AtomIdx> &atomidxs);

    AtomSelection& select(const QSet<CGIdx> &cgidxs);
    AtomSelection& deselect(const QSet<CGIdx> &cgidxs);
    AtomSelection& selectOnly(const QSet<CGIdx> &cgidxs);

    AtomSelection& select(const QSet<ResIdx> &residxs);
    AtomSelection& deselect(const QSet<ResIdx> &residxs);
    AtomSelection& selectOnly(const QSet<ResIdx> &residxs);

    AtomSelection& select(const QSet<ChainIdx> &chainidxs);
    AtomSelection& deselect(const QSet<ChainIdx> &chainidxs);
    AtomSelection& selectOnly(const QSet<ChainIdx> &chainidxs);

    AtomSelection& select(const QSet<SegIdx> &segidxs);
    AtomSelection& deselect(const QSet<SegIdx> &segidxs);
    AtomSelection& selectOnly(const QSet<SegIdx> &segidxs);

    AtomSelection& select(const QList<AtomIdx> &atomidxs);
    AtomSelection& deselect(const QList<AtomIdx> &atomidxs);
    AtomSelection& selectOnly(const QList<AtomIdx> &atomidxs);

    AtomSelection& select(const QList<CGIdx> &cgidxs);
    AtomSelection& deselect(const QList<CGIdx> &cgidxs);
    AtomSelection& selectOnly(const QList<CGIdx> &cgidxs);

    AtomSelection& select(const QList<ResIdx> &residxs);
    AtomSelection& deselect(const QList<ResIdx> &residxs);
    AtomSelection& selectOnly(const QList<ResIdx> &residxs);

    AtomSelection& select(const QList<ChainIdx> &chainidxs);
    AtomSelection& deselect(const QList<ChainIdx> &chainidxs);
    AtomSelection& selectOnly(const QList<ChainIdx> &chainidxs);

    AtomSelection& select(const QList<SegIdx> &segidxs);
    AtomSelection& deselect(const QList<SegIdx> &segidxs);
    AtomSelection& selectOnly(const QList<SegIdx> &segidxs);

    AtomSelection& select(const AtomID &atomid);
    AtomSelection& deselect(const AtomID &atomid);
    AtomSelection& selectOnly(const AtomID &atomid);

    AtomSelection& select(const CGID &cgid);
    AtomSelection& deselect(const CGID &cgid);
    AtomSelection& selectOnly(const CGID &cgid);

    AtomSelection& select(const ResID &resid);
    AtomSelection& deselect(const ResID &resid);
    AtomSelection& selectOnly(const ResID &resid);

    AtomSelection& select(const ChainID &chainid);
    AtomSelection& deselect(const ChainID &chainid);
    AtomSelection& selectOnly(const ChainID &chainid);
    
    AtomSelection& select(const SegID &segid);
    AtomSelection& deselect(const SegID &segid);
    AtomSelection& selectOnly(const SegID &segid);

    AtomSelection& select(const AtomSelection &selection);
    AtomSelection& deselect(const AtomSelection &selection);
    AtomSelection& selectOnly(const AtomSelection &selection);

    AtomSelection& invert();
 
    bool intersects(AtomIdx atomidx) const;
    bool intersects(CGIdx cgidx) const;
    bool intersects(ResIdx residx) const;
    bool intersects(ChainIdx chainidx) const;
    bool intersects(SegIdx segidx) const;
    
    bool intersects(const AtomID &atomid) const;
    bool intersects(const CGID &cgid) const;
    bool intersects(const ResID &resid) const;
    bool intersects(const ChainID &chainid) const;
    bool intersects(const SegID &segid) const;

    bool intersects(const AtomSelection &selection) const;
    
    bool contains(AtomIdx atomidx) const;
    bool contains(CGIdx cgidx) const;
    bool contains(ResIdx residx) const;
    bool contains(ChainIdx chainidx) const;
    bool contains(SegIdx segidx) const;
    
    bool contains(const AtomID &atomid) const;
    bool contains(const CGID &cgid) const;
    bool contains(const ResID &resid) const;
    bool contains(const ChainID &chainid) const;
    bool contains(const SegID &segid) const;

    bool contains(const AtomSelection &selection) const;

    AtomSelection& intersect(AtomIdx atomidx);
    AtomSelection& intersect(CGIdx cgidx);
    AtomSelection& intersect(ResIdx residx);
    AtomSelection& intersect(ChainIdx chainidx);
    AtomSelection& intersect(SegIdx segidx);

    AtomSelection& intersect(const QSet<AtomIdx> &atomidx);
    AtomSelection& intersect(const QSet<CGIdx> &cgidx);
    AtomSelection& intersect(const QSet<ResIdx> &residx);
    AtomSelection& intersect(const QSet<ChainIdx> &chainidx);
    AtomSelection& intersect(const QSet<SegIdx> &segidx);

    AtomSelection& intersect(const QList<AtomIdx> &atomidx);
    AtomSelection& intersect(const QList<CGIdx> &cgidx);
    AtomSelection& intersect(const QList<ResIdx> &residx);
    AtomSelection& intersect(const QList<ChainIdx> &chainidx);
    AtomSelection& intersect(const QList<SegIdx> &segidx);
    
    AtomSelection& intersect(const AtomID &atomid);
    AtomSelection& intersect(const CGID &cgid);
    AtomSelection& intersect(const ResID &resid);
    AtomSelection& intersect(const ChainID &chainid);
    AtomSelection& intersect(const SegID &segid);

    AtomSelection& intersect(const AtomSelection &selection);

    AtomSelection& unite(AtomIdx atomidx);
    AtomSelection& unite(CGIdx cgidx);
    AtomSelection& unite(ResIdx residx);
    AtomSelection& unite(ChainIdx chainidx);
    AtomSelection& unite(SegIdx segidx);

    AtomSelection& unite(const QSet<AtomIdx> &atomidx);
    AtomSelection& unite(const QSet<CGIdx> &cgidx);
    AtomSelection& unite(const QSet<ResIdx> &residx);
    AtomSelection& unite(const QSet<ChainIdx> &chainidx);
    AtomSelection& unite(const QSet<SegIdx> &segidx);

    AtomSelection& unite(const QList<AtomIdx> &atomidx);
    AtomSelection& unite(const QList<CGIdx> &cgidx);
    AtomSelection& unite(const QList<ResIdx> &residx);
    AtomSelection& unite(const QList<ChainIdx> &chainidx);
    AtomSelection& unite(const QList<SegIdx> &segidx);
    
    AtomSelection& unite(const AtomID &atomid);
    AtomSelection& unite(const CGID &cgid);
    AtomSelection& unite(const ResID &resid);
    AtomSelection& unite(const ChainID &chainid);
    AtomSelection& unite(const SegID &segid);

    AtomSelection& unite(const AtomSelection &selection);

    AtomSelection& unite(const QList<AtomSelection> &selections);

    AtomSelection& subtract(AtomIdx atomidx);
    AtomSelection& subtract(CGIdx cgidx);
    AtomSelection& subtract(ResIdx residx);
    AtomSelection& subtract(ChainIdx chainidx);
    AtomSelection& subtract(SegIdx segidx);

    AtomSelection& subtract(const QSet<AtomIdx> &atomidx);
    AtomSelection& subtract(const QSet<CGIdx> &cgidx);
    AtomSelection& subtract(const QSet<ResIdx> &residx);
    AtomSelection& subtract(const QSet<ChainIdx> &chainidx);
    AtomSelection& subtract(const QSet<SegIdx> &segidx);

    AtomSelection& subtract(const QList<AtomIdx> &atomidx);
    AtomSelection& subtract(const QList<CGIdx> &cgidx);
    AtomSelection& subtract(const QList<ResIdx> &residx);
    AtomSelection& subtract(const QList<ChainIdx> &chainidx);
    AtomSelection& subtract(const QList<SegIdx> &segidx);
    
    AtomSelection& subtract(const AtomID &atomid);
    AtomSelection& subtract(const CGID &cgid);
    AtomSelection& subtract(const ResID &resid);
    AtomSelection& subtract(const ChainID &chainid);
    AtomSelection& subtract(const SegID &segid);

    AtomSelection& subtract(const AtomSelection &selection);

    AtomSelection& mask(AtomIdx atomidx);
    AtomSelection& mask(CGIdx cgidx);
    AtomSelection& mask(ResIdx residx);
    AtomSelection& mask(ChainIdx chainidx);
    AtomSelection& mask(SegIdx segidx);

    AtomSelection& mask(const QSet<AtomIdx> &atomidx);
    AtomSelection& mask(const QSet<CGIdx> &cgidx);
    AtomSelection& mask(const QSet<ResIdx> &residx);
    AtomSelection& mask(const QSet<ChainIdx> &chainidx);
    AtomSelection& mask(const QSet<SegIdx> &segidx);

    AtomSelection& mask(const QList<AtomIdx> &atomidx);
    AtomSelection& mask(const QList<CGIdx> &cgidx);
    AtomSelection& mask(const QList<ResIdx> &residx);
    AtomSelection& mask(const QList<ChainIdx> &chainidx);
    AtomSelection& mask(const QList<SegIdx> &segidx);
    
    AtomSelection& mask(const AtomID &atomid);
    AtomSelection& mask(const CGID &cgid);
    AtomSelection& mask(const ResID &resid);
    AtomSelection& mask(const ChainID &chainid);
    AtomSelection& mask(const SegID &segid);

    AtomSelection& mask(const AtomSelection &selection);

    QVector<AtomIdx> selectedAtoms() const;

    QSet<Index> selectedAtoms(CGIdx cgid) const;
    
    QList<CGIdx> selectedCutGroups() const;
    QList<ResIdx> selectedResidues() const;
    QList<ChainIdx> selectedChains() const;
    QList<SegIdx> selectedSegments() const;

    void assertSelected(AtomIdx atomidx) const;
    void assertSelected(const AtomID &atomid) const;
    
    bool isCompatibleWith(const MoleculeInfoData &molinfo) const;

    void assertCompatibleWith(const MoleculeData &moldata) const;
    void assertCompatibleWith(const MoleculeView &molview) const;
    void assertCompatibleWith(const MoleculeInfoData &molinfo) const;
    void assertCompatibleWith(const AtomSelection &other) const;

    template<class T>
    void assertCompatibleWith(const AtomProperty<T> &prop) const;

private:
    bool _pvt_selected(const CGAtomIdx &cgatomidx) const;
    bool _pvt_selected(AtomIdx atomidx) const;
    
    bool _pvt_selectedAll(CGIdx cgidx) const;
    bool _pvt_selectedAll(const QVector<CGAtomIdx> &atomidxs) const;

    int _pvt_nSelected(ResIdx residx) const;

    void _pvt_select(AtomIdx atomidx);
    void _pvt_deselect(AtomIdx atomidx);
    
    void _pvt_select(CGIdx cgidx);
    void _pvt_deselect(CGIdx cgidx);
    
    void _pvt_select(const CGAtomIdx &cgatomidx);
    void _pvt_deselect(const CGAtomIdx &cgatomidx);
    
    void _pvt_select(CGIdx cgidx, const QSet<Index> &atoms);
    
    void _pvt_select(const QVector<CGAtomIdx> &cgatomidxs);
    void _pvt_deselect(const QVector<CGAtomIdx> &cgatomidxs);

    void _pvt_select(const AtomSelection &selection);

    template<class IDXS>
    void _pvt_selectAtoms(const IDXS &atoms);
    
    template<class IDXS>
    void _pvt_deselectAtoms(const IDXS &atoms);

    /** The indicies of selected atoms, arranged by CGIdx */
    QHash< CGIdx, QSet<Index> > selected_atoms;

    /** The MoleculeInfo describing the molecule whose parts
        are being selected by this object */
    SharedDataPointer<MoleculeInfoData> d;

    /** The total number of selected atoms */
    qint32 nselected;
};

}

SIRE_ALWAYS_INLINE SireMol::AtomSelection operator+(const SireMol::AtomSelection &a,
                                        const SireMol::AtomSelection &b)
{
    SireMol::AtomSelection ret(a);
    return ret.unite(b);
}

SIRE_ALWAYS_INLINE SireMol::AtomSelection operator-(const SireMol::AtomSelection &a,
                                        const SireMol::AtomSelection &b)
{
    SireMol::AtomSelection ret(a);
    return ret.subtract(b);
}

Q_DECLARE_METATYPE(SireMol::AtomSelection);

SIRE_EXPOSE_CLASS( SireMol::AtomSelection )

SIRE_END_HEADER

#endif
