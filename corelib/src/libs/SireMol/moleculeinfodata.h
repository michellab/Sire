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

#ifndef SIREMOL_MOLECULEINFODATA_H
#define SIREMOL_MOLECULEINFODATA_H

#include <QVector>
#include <QMultiHash>
#include <QSet>
#include <QUuid>

#include "molinfo.h"

#include "cgatomidx.h"
#include "molname.h"

#include "segname.h"
#include "segidx.h"
#include "segidentifier.h"

#include "chainname.h"
#include "chainidx.h"
#include "chainidentifier.h"

#include "cgname.h"
#include "cgidx.h"
#include "cgidentifier.h"

#include "resname.h"
#include "resnum.h"
#include "residx.h"
#include "residentifier.h"

#include "atomname.h"
#include "atomnum.h"
#include "atomidx.h"
#include "atomidentifier.h"

#include "atomsin.hpp"
#include "resin.hpp"

#include "SireBase/shareddatapointer.hpp"
#include "SireBase/refcountdata.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class MoleculeInfoData;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::MoleculeInfoData&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::MoleculeInfoData&);

namespace SireMol
{

class Atom;
class CutGroup;
class Residue;
class Chain;
class Segment;

class MoleculeData;
class AtomSelection;

class StructureEditor;

namespace detail
{
class MolInfoRegistry;

class AtomInfo;
class ResInfo;
class CGInfo;
class ChainInfo;
class SegInfo;
}

/** This is the implicitly shared class that is used to provide information
    about the arrangement of atoms in a molecule, specifically
    how one method of indexing the atoms can be mapped to another method.
    In so doing, a MoleculeInfoData object contains all of the information
    about which atoms are in a molecule, what they are called, how they
    are arranged into residues, cutgroups and segments, and how all
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
class SIREMOL_EXPORT MoleculeInfoData : public MolInfo, public SireBase::RefCountData
{

friend SIREMOL_EXPORT QDataStream& ::operator<<(QDataStream&, const MoleculeInfoData&);
friend SIREMOL_EXPORT QDataStream& ::operator>>(QDataStream&, MoleculeInfoData&);

public:
    MoleculeInfoData();

    MoleculeInfoData(const StructureEditor &editor);

    MoleculeInfoData(const MoleculeInfoData &other);

    ~MoleculeInfoData();

    static const char* typeName();

    const char* what() const
    {
        return MoleculeInfoData::typeName();
    }

    MoleculeInfoData& operator=(const MoleculeInfoData &other);

    bool operator==(const MoleculeInfoData &other) const;
    bool operator!=(const MoleculeInfoData &other) const;

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

    SegIdx number(const SegID &segid) const;
    SegIdx number(SegIdx segidx) const;

    CGIdx number(const CGID &cgid) const;
    CGIdx number(CGIdx cgidx) const;

    ChainIdx number(const ChainID &chainid) const;
    ChainIdx number(ChainIdx chainidx) const;

    ResNum number(const ResID &resid) const;
    ResNum number(ResIdx residx) const;

    AtomNum number(const AtomID &atomid) const;
    AtomNum number(AtomIdx atomidx) const;

    MoleculeInfoData rename(AtomIdx atomidx, const AtomName &newname) const;
    MoleculeInfoData renumber(AtomIdx atomidx, const AtomNum &newnum) const;

    MoleculeInfoData rename(ResIdx residx, const ResName &newname) const;
    MoleculeInfoData renumber(ResIdx residx, const ResNum &newnum) const;

    MoleculeInfoData rename(CGIdx cgidx, const CGName &newname) const;
    MoleculeInfoData rename(ChainIdx chainidx, const ChainName &newname) const;
    MoleculeInfoData rename(SegIdx segidx, const SegName &newname) const;

    MoleculeInfoData renumber( const QHash<AtomNum,AtomNum> &atomnums ) const;
    MoleculeInfoData renumber( const QHash<ResNum,ResNum> &resnums) const;
    MoleculeInfoData renumber( const QHash<AtomNum,AtomNum> &atomnums,
                               const QHash<ResNum,ResNum> &resnums ) const;

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
    CGIdx cgIdx(const ResIdx &residx) const;
    CGIdx cgIdx(const ResID &resid) const;

    bool isAtomCutting() const;
    bool isResidueCutting() const;
    bool isMoleculeCutting() const;

    bool isResidueCutting(const ResIdx &residx) const;
    bool isResidueCutting(const ResID &resid) const;

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

    bool isEmpty() const;

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

    QList<AtomIdx> mapNoThrow(const AtomName &name) const;
    QList<AtomIdx> mapNoThrow(const AtomNum &num) const;
    QList<AtomIdx> mapNoThrow(const AtomIdx &idx) const;
    QList<AtomIdx> mapNoThrow(const AtomID &atomid) const;

    QList<ResIdx> mapNoThrow(const ResName &name) const;
    QList<ResIdx> mapNoThrow(const ResNum &num) const;
    QList<ResIdx> mapNoThrow(const ResIdx &idx) const;
    QList<ResIdx> mapNoThrow(const ResID &resid) const;

    QList<ChainIdx> mapNoThrow(const ChainName &name) const;
    QList<ChainIdx> mapNoThrow(const ChainIdx &idx) const;
    QList<ChainIdx> mapNoThrow(const ChainID &chainid) const;

    QList<SegIdx> mapNoThrow(const SegName &name) const;
    QList<SegIdx> mapNoThrow(const SegIdx &idx) const;
    QList<SegIdx> mapNoThrow(const SegID &segid) const;

    QList<CGIdx> mapNoThrow(const CGName &name) const;
    QList<CGIdx> mapNoThrow(const CGIdx &idx) const;
    QList<CGIdx> mapNoThrow(const CGID &cgid) const;

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

    void squeeze(const MoleculeInfoData &other) const;

    void assertCompatibleWith(const AtomSelection &selected_atoms) const;

    void assertContains(AtomIdx atomidx) const;
    void assertContains(CGIdx cgidx) const;
    void assertContains(ResIdx residx) const;
    void assertContains(ChainIdx chainidx) const;
    void assertContains(SegIdx segidx) const;

    void assertEqualTo(const MoleculeInfoData &other) const;

    static const MoleculeInfoData& null();

private:

    void rebuildNameAndNumberIndexes();

    bool _pvt_hasSameFingerprint(const MoleculeInfoData &other);

    QList<AtomIdx> _pvt_getAtomsIn(const QList<ResIdx> &residxs) const;
    QList<AtomIdx> _pvt_getAtomsIn(const QList<ResIdx> &residxs,
                                   const AtomName &name) const;

    int _pvt_nAtoms(const QVector<ResIdx> &residxs) const;
    int _pvt_nAtoms(const QList<ResIdx> &residxs) const;
    int _pvt_nAtoms(ChainIdx chainidx) const;

    QVector<CGAtomIdx> _pvt_cgAtomIdxs(const QList<AtomIdx> &atomidxs) const;

    /** The unique ID that identifies this particular
        molecule layout. */
    QUuid uid;

    /** All of the atoms in the molecule, in the order they were
        added to the molecule */
    QVector<detail::AtomInfo> atoms_by_index;

    /** Hash mapping atom names to atom indicies */
    QMultiHash<QString,AtomIdx> atoms_by_name;

    /** Hash mapping atom numbers to atom indicies */
    QMultiHash<AtomNum,AtomIdx> atoms_by_num;

    /** All of the residues in this molecule, arranged in the
        order that they appear in this molecule */
    QVector<detail::ResInfo> res_by_index;

    /** Hash mapping residue names to residue indicies */
    QMultiHash<QString,ResIdx> res_by_name;

    /** Hash mapping residue numbers to residue indicies */
    QMultiHash<ResNum,ResIdx> res_by_num;

    /** All of the chains in this molecule, arranged in the
        order that they appear in this molecule */
    QVector<detail::ChainInfo> chains_by_index;

    /** Hash mapping chain names to chain indicies */
    QMultiHash<QString,ChainIdx> chains_by_name;

    /** All of the segments in this molecule, arranged in the
        order that they appear in this molecule */
    QVector<detail::SegInfo> seg_by_index;

    /** Hash mapping segment names to segment indicies */
    QMultiHash<QString,SegIdx> seg_by_name;

    /** All of the CutGroups in this molecule, arranged in the
        order that they appear in this molecule */
    QVector<detail::CGInfo> cg_by_index;

    /** Hash mapping CutGroup name to CutGroup indicies */
    QMultiHash<QString,CGIdx> cg_by_name;

    /** The cutting scheme for the molecule. This is either
        0 = unknown, 1 = atom, 2 = residue or 3 = molecule */
    qint32 cutting_scheme;
};

} //end of namespace SireMol


Q_DECLARE_METATYPE(SireMol::MoleculeInfoData);

SIRE_END_HEADER

#endif
