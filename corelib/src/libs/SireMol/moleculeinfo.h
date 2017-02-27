/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2017  Christopher Woods
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

#ifndef SIREMOL_MOLECULEINFO_H
#define SIREMOL_MOLECULEINFO_H

#include "moleculeinfodata.h"
#include "atomselection.h"

#include "SireBase/shareddatapointer.hpp"

namespace SireMol
{
class MoleculeInfo;
}

QDataStream& operator<<(QDataStream&, const SireMol::MoleculeInfo&);
QDataStream& operator>>(QDataStream&, SireMol::MoleculeInfo&);

namespace SireMol
{

/** This is the class that is used to provide information
    about the arrangement of atoms in a molecule, specifically
    how one method of indexing the atoms can be mapped to another method.
    In so doing, a MoleculeInfo object contains all of the information
    about which atoms are in a molecule, what they are called, how they
    are arranged into residues, chains and segments, and how all
    of the different parts of the molecule are identified. However,
    importantly, this object does not contain any information about
    how the atoms in the molecule are connected together, where
    the atoms are in space, or any additional properties that are
    associated with the molecule (or indeed the name or number of the molecule!)
    
    Each layout is given a unique ID (UID) number, which is unique within a
    single invocation of Sire, and uniquely identifies a MoleculeInfo
    layout within the program (thus allowing for a quick and simple
    test to ensure that molecules have the same layout of data).
    
    @author Christopher Woods
*/
class SIREMOL_EXPORT MoleculeInfo
        : public SireBase::ConcreteProperty<MoleculeInfo,SireBase::Property>
{

friend QDataStream& ::operator<<(QDataStream&, const MoleculeInfo&);
friend QDataStream& ::operator>>(QDataStream&, MoleculeInfo&);

public:
    MoleculeInfo();
    
    MoleculeInfo(const MoleculeView &molecule);
    
    MoleculeInfo(const SireBase::SharedDataPointer<MoleculeInfoData> &ptr);
    
    MoleculeInfo(const MoleculeInfo &other);
    
    ~MoleculeInfo();
    
    static const char* typeName();
    
    const char* what() const;
    
    MoleculeInfo& operator=(const MoleculeInfo &other);
    
    bool operator==(const MoleculeInfo &other) const;
    bool operator!=(const MoleculeInfo &other) const;
    
    const QUuid& UID() const;
    
    const ChainName& name(const ChainID &chainid) const;
    const ChainName& name(ChainIdx chainidx) const;
    
    const SegName& name(const SegID &segid) const;
    const SegName& name(SegIdx segidx) const;
    
    const ResName& name(const ResID &resid) const;
    const ResName& name(ResIdx residx) const;
    
    const CGName& name(const CGID &cgid) const;
    const CGName& name(CGIdx cgidx) const;
    
    const AtomName& name(const AtomID &atomid) const;
    const AtomName& name(AtomIdx atomidx) const;
    
    ResNum number(const ResID &resid) const;
    ResNum number(ResIdx residx) const;
    
    AtomNum number(const AtomID &atomid) const;
    AtomNum number(AtomIdx atomidx) const;
    
    MoleculeInfo rename(AtomIdx atomidx, const AtomName &newname) const;
    MoleculeInfo renumber(AtomIdx atomidx, const AtomNum &newnum) const;
    
    MoleculeInfo rename(ResIdx residx, const ResName &newname) const;
    MoleculeInfo renumber(ResIdx residx, const ResNum &newnum) const;
    
    MoleculeInfo rename(CGIdx cgidx, const CGName &newname) const;
    MoleculeInfo rename(ChainIdx chainidx, const ChainName &newname) const;
    MoleculeInfo rename(SegIdx segidx, const SegName &newname) const;
    
    const CGAtomIdx& cgAtomIdx(AtomIdx atomidx) const;
    const CGAtomIdx& cgAtomIdx(const AtomID &atomid) const;
    
    QVector<CGAtomIdx> cgAtomIdxs(AtomIdx atomidx) const;
    QVector<CGAtomIdx> cgAtomIdxs(CGIdx cgidx) const;
    QVector<CGAtomIdx> cgAtomIdxs(ResIdx residx) const;
    QVector<CGAtomIdx> cgAtomIdxs(ChainIdx chainidx) const;
    QVector<CGAtomIdx> cgAtomIdxs(SegIdx segidx) const;
    
    QVector<CGAtomIdx> cgAtomIdxs(const AtomID &atomid) const;
    QVector<CGAtomIdx> cgAtomIdxs(const CGID &cgid) const;
    QVector<CGAtomIdx> cgAtomIdxs(const ResID &resid) const;
    QVector<CGAtomIdx> cgAtomIdxs(const ChainID &chainid) const;
    QVector<CGAtomIdx> cgAtomIdxs(const SegID &segid) const;

    AtomIdx atomIdx(const AtomID &atomid) const;
    AtomIdx atomIdx(const CGAtomIdx &cgatomidx) const;

    ResIdx resIdx(const ResID &resid) const;
    ChainIdx chainIdx(const ChainID &chainid) const;
    SegIdx segIdx(const SegID &segid) const;
    CGIdx cgIdx(const CGID &cgid) const;

    QList<SegIdx> getSegments() const;
    QList<CGIdx> getCutGroups() const;
    QList<ChainIdx> getChains() const;
    QList<ResIdx> getResidues() const;
    
    const QList<ResIdx>& getResiduesIn(ChainIdx chainidx) const;
    QList<ResIdx> getResiduesIn(const ChainID &chainid) const;

    QList<AtomIdx> getAtoms() const;

    AtomIdx getAtom(CGIdx cgidx, int i) const;
    AtomIdx getAtom(ResIdx residx, int i) const;
    AtomIdx getAtom(ChainIdx chainidx, int i) const;
    AtomIdx getAtom(SegIdx segidx, int i) const;
    
    ResIdx getResidue(ChainIdx chainidx, int i) const;

    const QList<AtomIdx>& getAtomsIn(ResIdx residx) const;
    QList<AtomIdx> getAtomsIn(const ResID &resid) const;
    QList<AtomIdx> getAtomsIn(ResIdx residx, const AtomName &name) const;
    QList<AtomIdx> getAtomsIn(const ResID &resid,
                              const AtomName &atomname) const;

    QList<AtomIdx> getAtomsIn(ChainIdx chainidx) const;
    QList<AtomIdx> getAtomsIn(const ChainID &chainid) const;
    QList<AtomIdx> getAtomsIn(ChainIdx chainidx, 
                              const AtomName &atomname) const;
    QList<AtomIdx> getAtomsIn(const ChainID &chainid,
                              const AtomName &atomname) const;
                              
    const QList<AtomIdx>& getAtomsIn(CGIdx cgidx) const;
    QList<AtomIdx> getAtomsIn(const CGID &cgid) const;
    
    const QList<AtomIdx>& getAtomsIn(SegIdx segidx) const;
    QList<AtomIdx> getAtomsIn(const SegID &segid) const;

    bool isWithinResidue(AtomIdx atomidx) const;
    bool isWithinResidue(const AtomID &atomid) const;
    
    bool isWithinChain(AtomIdx atomidx) const;
    bool isWithinChain(const AtomID &atomid) const;
    
    bool isWithinSegment(AtomIdx atomidx) const;
    bool isWithinSegment(const AtomID &atomid) const;
    
    bool isWithinChain(ResIdx residx) const;
    bool isWithinChain(const ResID &resid) const;

    ChainIdx parentChain(ResIdx residx) const;
    ChainIdx parentChain(const ResID &resid) const;
    
    ChainIdx parentChain(AtomIdx atomidx) const;
    ChainIdx parentChain(const AtomID &atomid) const;
    
    ResIdx parentResidue(AtomIdx atomidx) const;
    ResIdx parentResidue(const AtomID &atomid) const;
    
    SegIdx parentSegment(AtomIdx atomidx) const;
    SegIdx parentSegment(const AtomID &atomid) const;
    
    CGIdx parentCutGroup(AtomIdx atomidx) const;
    CGIdx parentCutGroup(const AtomID &atomid) const;
    
    bool contains(ResIdx residx, AtomIdx atomidx) const;
    bool contains(ChainIdx chainidx, AtomIdx atomidx) const;
    bool contains(SegIdx segidx, AtomIdx atomidx) const;
    bool contains(CGIdx cgidx, AtomIdx atomidx) const;
    bool contains(ChainIdx chainidx, ResIdx residx) const;
    
    bool contains(ResIdx residx, const AtomID &atomid) const;
    bool contains(ChainIdx chainidx, const AtomID &atomid) const;
    bool contains(SegIdx segidx, const AtomID &atomid) const;
    bool contains(CGIdx cgidx, const AtomID &atomid) const;
    bool contains(ChainIdx chainidx, const ResID &resid) const;

    bool intersects(ResIdx residx, const AtomID &atomid) const;
    bool intersects(ChainIdx chainidx, const AtomID &atomid) const;
    bool intersects(SegIdx segidx, const AtomID &atomid) const;
    bool intersects(CGIdx cgidx, const AtomID &atomid) const;
    bool intersects(ChainIdx chainidx, const ResID &resid) const;
    
    int nAtoms() const;
    
    int nAtoms(const ChainID &chainid) const;
    int nAtoms(ChainIdx chainidx) const;
    
    int nAtoms(const ResID &resid) const;
    int nAtoms(ResIdx residx) const;
    
    int nAtoms(const SegID &segid) const;
    int nAtoms(SegIdx segidx) const;
    
    int nAtoms(const CGID &cgid) const;
    int nAtoms(CGIdx cgidx) const;
    
    int nResidues() const;
    int nResidues(const ChainID &chainid) const;
    int nResidues(ChainIdx chainidx) const;
    
    int nChains() const;
    int nCutGroups() const;
    int nSegments() const;
    
    QList<ResIdx> map(const ResName &name) const;
    QList<ResIdx> map(ResNum num) const;
    QList<ResIdx> map(ResIdx idx) const;
    QList<ResIdx> map(const ResID &resid) const;
    
    QList<ChainIdx> map(const ChainName &name) const;
    QList<ChainIdx> map(ChainIdx idx) const;
    QList<ChainIdx> map(const ChainID &chainid) const;

    QList<SegIdx> map(const SegName &name) const;
    QList<SegIdx> map(SegIdx idx) const;
    QList<SegIdx> map(const SegID &segid) const;
   
    QList<CGIdx> map(const CGName &name) const;
    QList<CGIdx> map(CGIdx idx) const;
    QList<CGIdx> map(const CGID &cgid) const;
    
    QList<AtomIdx> map(const AtomName &name) const;
    QList<AtomIdx> map(AtomNum num) const;
    QList<AtomIdx> map(AtomIdx idx) const;
    QList<AtomIdx> map(const AtomID &atomid) const;

    void squeeze(const MoleculeInfo &other) const;

    void assertCompatibleWith(const AtomSelection &selected_atoms) const;
    void assertCompatibleWith(const MoleculeView &molecule) const;

    void assertContains(AtomIdx atomidx) const;
    void assertContains(CGIdx cgidx) const;
    void assertContains(ResIdx residx) const;
    void assertContains(ChainIdx chainidx) const;
    void assertContains(SegIdx segidx) const;

    void assertEqualTo(const MoleculeInfo &other) const;
    void assertEqualTo(const MoleculeInfoData &other) const;

    operator const MoleculeInfoData&() const;

    const MoleculeInfoData& data() const;

private:
    SireBase::SharedDataPointer<MoleculeInfoData> d;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Allow automatic casting to a MoleculeInfoData object */
inline MoleculeInfo::operator const MoleculeInfoData&() const
{
    return *d;
}

/** Retrieve a reference to the underlying data object - this is used for 
    compaibility with old code. It should not be used in new code */
inline const MoleculeInfoData& MoleculeInfo::data() const
{
    return *d;
}

#endif // SIRE_SKIP_INLINE_FUNCTIONS

}

Q_DECLARE_METATYPE( SireMol::MoleculeInfo )

SIRE_EXPOSE_CLASS( SireMol::MoleculeInfo )

#endif
