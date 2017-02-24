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

#include "moleculeinfo.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMol;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<MoleculeInfo> r_molinfo;

QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds, const MoleculeInfo &molinfo)
{
    writeHeader(ds, r_molinfo, 1);
    
    SharedDataStream sds(ds);
    sds << molinfo.d;
    return ds;
}

QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds, MoleculeInfo &molinfo)
{
    VerisonID v = readHeader(ds, r_molinfo);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> molinfo.d;
    }
    else
        throw version_error(v, "1", r_molinfo, CODELOC);
    
    return ds;
}

/** Null constructor */
MoleculeInfo::MoleculeInfo()
{}

/** Construct to get the MoleculeInfo for the passed molecule */
MoleculeInfo::MoleculeInfo(const MoleculeView &molecule) : d(molecule.data().info())
{}

/** Construct from the passed shared pointer */
MoleculeInfo::MoleculeInfo(const SharedDataPointer<MoleculeInfoData> &ptr) : d(ptr)
{}

/** Copy constructor */
MoleculeInfo::MoleculeInfo(const MoleculeInfo &other) : d(other.d)
{}

/** Destructor */
MoleculeInfo::~MoleculeInfo()
{}

const char* MoleculeInfo::typeName()
{
    return QMetaType::typeName( qMetaTypeId<MoleculeInfo>() );
}

const char* MoleculeInfo::what() const
{
    return MoleculeInfo::typeName();
}

MoleculeInfo& MoleculeInfo::operator=(const MoleculeInfo &other)
{
    d = other.d;
    return *this;
}

bool MoleculeInfo::operator==(const MoleculeInfo &other) const
{
    return this->UID() == other.UID();
}

bool MoleculeInfo::operator!=(const MoleculeInfo &other) const
{
    return not operator==(other);
}

/** Return the unique ID of this layout. Each unique layout has its
    own unique ID. You can use this to quickly check if two molecules
    have the same layout */
const QUuid& MoleculeInfo::UID() const
{
    return d->UID();
}

/** Return the name of the identified chain */
const ChainName& MoleculeInfo::name(const ChainID &chainid) const
{
    return d->name(chainid);
}

/** Return the name of the identified chain */
const ChainName& MoleculeInfo::name(ChainIdx chainidx) const
{
    return d->name(chainidx);
}

/** Return the name of the identified segment */
const SegName& MoleculeInfo::name(const SegID &segid) const
{
    return d->name(segid);
}

/** Return the name of the identified segment */
const SegName& MoleculeInfo::name(SegIdx segidx) const
{
    return d->name(segidx);
}

/** Return the name of the identified residue */
const ResName& MoleculeInfo::name(const ResID &resid) const
{
    return d->name(resid);
}

/** Return the name of the identified residue */
const ResName& MoleculeInfo::name(ResIdx residx) const
{
    return d->name(residx);
}

/** Return the name of the identified CutGroup */
const CGName& MoleculeInfo::name(const CGID &cgid) const
{
    return d->name(cgid);
}

/** Return the name of the identified CutGroup */
const CGName& MoleculeInfo::name(CGIdx cgidx) const
{
    return d->name(cgidx);
}

/** Return the name of the identified atom */
const AtomName& MoleculeInfo::name(const AtomID &atomid) const
{
    return d->name(atomid);
}

/** Return the name of the identified atom */
const AtomName& MoleculeInfo::name(AtomIdx atomidx) const
{
    return d->name(atomidx);
}

/** Return the number of the identified residue */
ResNum MoleculeInfo::number(const ResID &resid) const
{
    return d->number(resid);
}

/** Return the number of the identified residue */
ResNum MoleculeInfo::number(ResIdx residx) const
{
    return d->number(residx);
}

/** Return the number of the identified atom */
AtomNum MoleculeInfo::number(const AtomID &atomid) const
{
    return d->number(atomid);
}

/** Return the number of the identified atom */
AtomNum MoleculeInfo::number(AtomIdx atomidx) const
{
    return d->number(atomidx);
}

/** Return a copy of this MoleculeInfo where the identified atom has been
    renamed to 'newname' */
MoleculeInfo MoleculeInfo::rename(AtomIdx atomidx, const AtomName &newname) const
{
    return MoleculeInfo( d->rename(atomidx,newname) );
}

/** Return a copy of this MoleculeInfo where the identified atom has been
    renumbered to 'newnum' */
MoleculeInfo MoleculeInfo::renumber(AtomIdx atomidx, const AtomNum &newnum) const
{
    return MoleculeInfo( d->renumber(atomidx,newnum) );
}

/** Return a copy of this MoleculeInfo where the identified residue has been
    renamed to 'newname' */
MoleculeInfo MoleculeInfo::rename(ResIdx residx, const ResName &newname) const
{
    return MoleculeInfo( d->rename(residx,newname) );
}

/** Return a copy of this MoleculeInfo where the identified residue has been
    renumbered to 'newnum' */
MoleculeInfo MoleculeInfo::renumber(ResIdx residx, const ResNum &newnum) const
{
    return MoleculeInfo( d->renumber(residx,newnum) );
}

/** Return a copy of this MoleculeInfo where the identified CutGroup has been
    renamed to 'newname' */
MoleculeInfo MoleculeInfo::rename(CGIdx cgidx, const CGName &newname) const
{
    return MoleculeInfo( d->rename(cgidx,newname) );
}

/** Return a copy of this MoleculeInfo where the identified chain has been
    renamed to 'newname' */
MoleculeInfo MoleculeInfo::rename(ChainIdx chainidx, const ChainName &newname) const
{
    return MoleculeInfo( d->rename(chainidx,newname) );
}

/** Return a copy of this MoleculeInfo where the identified segment has been
    renamed to 'newname' */
MoleculeInfo MoleculeInfo::rename(SegIdx segidx, const SegName &newname) const
{
    return MoleculeInfo( d->rename(segidx,newname) );
}

/** Return the combined CutGroup / AtomIndex of the identified atom */
const CGAtomIdx& MoleculeInfo::cgAtomIdx(AtomIdx atomidx) const
{
    return d->cgAtomIdx(atomidx);
}

/** Return the combined CutGroup / AtomIndex of the identified atom */
const CGAtomIdx& MoleculeInfo::cgAtomIdx(const AtomID &atomid) const
{
    return d->cgAtomIdx(atomid);
}

/** Return the CGAtomIdxs of all of the identified atoms */
QVector<CGAtomIdx> MoleculeInfo::cgAtomIdxs(AtomIdx atomidx) const
{
    return d->cgAtomIdxs(atomidx);
}

/** Return the CGAtomIdxs of the atoms in all of the identified CutGroups */
QVector<CGAtomIdx> MoleculeInfo::cgAtomIdxs(CGIdx cgidx) const
{
    return d->cgAtomIdxs(cgidx);
}

/** Return the CGAtomIdxs of the atoms in all of the identified residues */
QVector<CGAtomIdx> MoleculeInfo::cgAtomIdxs(ResIdx residx) const
{
    return d->cgAtomIdxs(residx);
}

/** Return the CGAtomIdxs of the atoms in all of the identified chains */
QVector<CGAtomIdx> MoleculeInfo::cgAtomIdxs(ChainIdx chainidx) const
{
    return d->cgAtomIdxs(chainidx);
}

/** Return the CGAtomIdxs of the atoms in all of the identified segments */
QVector<CGAtomIdx> MoleculeInfo::cgAtomIdxs(SegIdx segidx) const
{
    return d->cgAtomIdxs(segidx);
}

/** Return the CGAtomIdxs of all of the identified atoms */
QVector<CGAtomIdx> MoleculeInfo::cgAtomIdxs(const AtomID &atomid) const
{
    return d->cgAtomIdxs(atomid);
}

/** Return the CGAtomIdxs of the atoms in all of the identified CutGroups */
QVector<CGAtomIdx> MoleculeInfo::cgAtomIdxs(const CGID &cgid) const
{
    return d->cgAtomIdxs(cgid);
}

/** Return the CGAtomIdxs of the atoms in all of the identified residues */
QVector<CGAtomIdx> MoleculeInfo::cgAtomIdxs(const ResID &resid) const
{
    return d->cgAtomIdxs(resid);
}

/** Return the CGAtomIdxs of the atoms in all of the identified chains */
QVector<CGAtomIdx> MoleculeInfo::cgAtomIdxs(const ChainID &chainid) const
{
    return d->cgAtomIdxs(chainid);
}

/** Return the CGAtomIdxs of the atoms in all of the identified segments */
QVector<CGAtomIdx> MoleculeInfo::cgAtomIdxs(const SegID &segid) const
{
    return d->cgAtomIdxs(segid);
}

/** Return the index of the identified atom */
AtomIdx MoleculeInfo::atomIdx(const AtomID &atomid) const
{
    return d->atomIdx(atomid);
}

/** Return the index of the identified atom */
AtomIdx MoleculeInfo::atomIdx(const CGAtomIdx &cgatomidx) const
{
    return d->atomIdx(cgatomidx);
}

/** Return the index of the identified residue */
ResIdx MoleculeInfo::resIdx(const ResID &resid) const
{
    return d->resIdx(resid);
}

/** Return the index of the identified chain */
ChainIdx MoleculeInfo::chainIdx(const ChainID &chainid) const
{
    return d->chainIdx(chainid);
}

/** Return the index of the identified segment */
SegIdx MoleculeInfo::segIdx(const SegID &segid) const
{
    return d->segIdx(segid);
}

/** Return the index of the identified CutGroup */
CGIdx MoleculeInfo::cgIdx(const CGID &cgid) const
{
    return d->cgIdx(cgid);
}

/** Return a list of the indicies of all segments */
QList<SegIdx> MoleculeInfo::getSegments() const
{
    return d->getSegments();
}

/** Return a list of the indicies of all CutGroups */
QList<CGIdx> MoleculeInfo::getCutGroups() const
{
    return d->getCutGroups();
}

/** Return a list of the indicies of all chains */
QList<ChainIdx> MoleculeInfo::getChains() const
{
    return d->getChains();
}

/** Return a list of the indicies of all residues */
QList<ResIdx> MoleculeInfo::getResidues() const
{
    return d->getResidues();
}

/** Return the indicies of residues in the identified chain(s) */
const QList<ResIdx>& MoleculeInfo::getResiduesIn(ChainIdx chainidx) const
{
    return d->getResiduesIn(chainidx);
}

/** Return the indicies of residues in the identified chain(s) */
QList<ResIdx> MoleculeInfo::getResiduesIn(const ChainID &chainid) const
{
    return d->getResiduesIn(chainid);
}

/** Return the indicies of atoms */
QList<AtomIdx> MoleculeInfo::getAtoms() const
{
    return d->getAtoms();
}

/** Return the index of the ith atom in the specified CutGroup */
AtomIdx MoleculeInfo::getAtom(CGIdx cgidx, int i) const
{
    return d->getAtom(cgidx,i);
}

/** Return the index of the ith atom in the specified residue */
AtomIdx MoleculeInfo::getAtom(ResIdx residx, int i) const
{
    return d->getAtom(residx,i);
}

/** Return the index of the ith atom in the specified chain */
AtomIdx MoleculeInfo::getAtom(ChainIdx chainidx, int i) const
{
    return d->getAtom(chainidx,i);
}

/** Return the index of the ith atom in the specified segment */
AtomIdx MoleculeInfo::getAtom(SegIdx segidx, int i) const
{
    return d->getAtom(segidx,i);
}

/** Return the index of the ith residue in the specified chain */
ResIdx MoleculeInfo::getResidue(ChainIdx chainidx, int i) const
{
    return d->getResidue(chainidx,i);
}

/** Return the indicies of all atoms in the specified residue(s) */
const QList<AtomIdx>& MoleculeInfo::getAtomsIn(ResIdx residx) const
{
    return d->getAtomsIn(residx);
}

/** Return the indicies of all atoms in the specified residue(s) */
QList<AtomIdx> MoleculeInfo::getAtomsIn(const ResID &resid) const
{
    return d->getAtomsIn(resid);
}

/** Return the indicies of all atoms called 'name' in the specified residue(s) */
QList<AtomIdx> MoleculeInfo::getAtomsIn(ResIdx residx, const AtomName &name) const
{
    return d->getAtomsIn(residx,name);
}

/** Return the indicies of all atoms called 'name' in the specified residue(s) */
QList<AtomIdx> MoleculeInfo::getAtomsIn(const ResID &resid,
                                        const AtomName &name) const
{
    return d->getAtomsIn(resid,name);
}

/** Return the indicies of all atoms in the specified chain(s) */
QList<AtomIdx> MoleculeInfo::getAtomsIn(ChainIdx chainidx) const
{
    return d->getAtomsIn(chainidx);
}

/** Return the indicies of all atoms in the specified chain(s) */
QList<AtomIdx> MoleculeInfo::getAtomsIn(const ChainID &chainid) const
{
    return d->getAtomsIn(chainid);
}

/** Return the indicies of all atoms called 'name' in the specified residue(s) */
QList<AtomIdx> MoleculeInfo::getAtomsIn(ChainIdx chainidx,
                                        const AtomName &name) const
{
    return d->getAtomsIn(chainidx,name);
}

/** Return the indicies of all atoms called 'name' in the specified residue(s) */
QList<AtomIdx> MoleculeInfo::getAtomsIn(const ChainID &chainid,
                                        const AtomName &name) const
{
    return d->getAtomsIn(chainid,name);
}
                          
/** Return the indicies of all atoms in the specified chain(s) */
const QList<AtomIdx>& MoleculeInfo::getAtomsIn(CGIdx cgidx) const
{
    return d->getAtomsIn(cgidx);
}

/** Return the indicies of all atoms in the specified CutGroup(s) */
QList<AtomIdx> MoleculeInfo::getAtomsIn(const CGID &cgid) const
{
    return d->getAtomsIn(cgid);
}

/** Return the indicies of all atoms in the specified segment(s) */
const QList<AtomIdx>& MoleculeInfo::getAtomsIn(SegIdx segidx) const
{
    return d->getAtomsIn(segidx);
}

/** Return the indicies of all atoms in the specified segment(s) */
QList<AtomIdx> MoleculeInfo::getAtomsIn(const SegID &segid) const
{
    return d->getAtomsIn(segid);
}

bool MoleculeInfo::isWithinResidue(AtomIdx atomidx) const;
bool MoleculeInfo::isWithinResidue(const AtomID &atomid) const;

bool MoleculeInfo::isWithinChain(AtomIdx atomidx) const;
bool MoleculeInfo::isWithinChain(const AtomID &atomid) const;

bool MoleculeInfo::isWithinSegment(AtomIdx atomidx) const;
bool MoleculeInfo::isWithinSegment(const AtomID &atomid) const;

bool MoleculeInfo::isWithinChain(ResIdx residx) const;
bool MoleculeInfo::isWithinChain(const ResID &resid) const;

ChainIdx MoleculeInfo::parentChain(ResIdx residx) const;
ChainIdx MoleculeInfo::parentChain(const ResID &resid) const;

ChainIdx MoleculeInfo::parentChain(AtomIdx atomidx) const;
ChainIdx MoleculeInfo::parentChain(const AtomID &atomid) const;

ResIdx MoleculeInfo::parentResidue(AtomIdx atomidx) const;
ResIdx MoleculeInfo::parentResidue(const AtomID &atomid) const;

SegIdx MoleculeInfo::parentSegment(AtomIdx atomidx) const;
SegIdx MoleculeInfo::parentSegment(const AtomID &atomid) const;

CGIdx MoleculeInfo::parentCutGroup(AtomIdx atomidx) const;
CGIdx MoleculeInfo::parentCutGroup(const AtomID &atomid) const;

bool MoleculeInfo::contains(ResIdx residx, AtomIdx atomidx) const;
bool MoleculeInfo::contains(ChainIdx chainidx, AtomIdx atomidx) const;
bool MoleculeInfo::contains(SegIdx segidx, AtomIdx atomidx) const;
bool MoleculeInfo::contains(CGIdx cgidx, AtomIdx atomidx) const;
bool MoleculeInfo::contains(ChainIdx chainidx, ResIdx residx) const;

bool MoleculeInfo::contains(ResIdx residx, const AtomID &atomid) const;
bool MoleculeInfo::contains(ChainIdx chainidx, const AtomID &atomid) const;
bool MoleculeInfo::contains(SegIdx segidx, const AtomID &atomid) const;
bool MoleculeInfo::contains(CGIdx cgidx, const AtomID &atomid) const;
bool MoleculeInfo::contains(ChainIdx chainidx, const ResID &resid) const;

bool MoleculeInfo::intersects(ResIdx residx, const AtomID &atomid) const;
bool MoleculeInfo::intersects(ChainIdx chainidx, const AtomID &atomid) const;
bool MoleculeInfo::intersects(SegIdx segidx, const AtomID &atomid) const;
bool MoleculeInfo::intersects(CGIdx cgidx, const AtomID &atomid) const;
bool MoleculeInfo::intersects(ChainIdx chainidx, const ResID &resid) const;

int MoleculeInfo::nAtoms() const;

int MoleculeInfo::nAtoms(const ChainID &chainid) const;
int MoleculeInfo::nAtoms(ChainIdx chainidx) const;

int MoleculeInfo::nAtoms(const ResID &resid) const;
int MoleculeInfo::nAtoms(ResIdx residx) const;

int MoleculeInfo::nAtoms(const SegID &segid) const;
int MoleculeInfo::nAtoms(SegIdx segidx) const;

int MoleculeInfo::nAtoms(const CGID &cgid) const;
int MoleculeInfo::nAtoms(CGIdx cgidx) const;

int MoleculeInfo::nResidues() const;
int MoleculeInfo::nResidues(const ChainID &chainid) const;
int MoleculeInfo::nResidues(ChainIdx chainidx) const;

int MoleculeInfo::nChains() const;
int MoleculeInfo::nCutGroups() const;
int MoleculeInfo::nSegments() const;

QList<ResIdx> MoleculeInfo::map(const ResName &name) const;
QList<ResIdx> MoleculeInfo::map(ResNum num) const;
QList<ResIdx> MoleculeInfo::map(ResIdx idx) const;
QList<ResIdx> MoleculeInfo::map(const ResID &resid) const;

QList<ChainIdx> MoleculeInfo::map(const ChainName &name) const;
QList<ChainIdx> MoleculeInfo::map(ChainIdx idx) const;
QList<ChainIdx> MoleculeInfo::map(const ChainID &chainid) const;

QList<SegIdx> MoleculeInfo::map(const SegName &name) const;
QList<SegIdx> MoleculeInfo::map(SegIdx idx) const;
QList<SegIdx> MoleculeInfo::map(const SegID &segid) const;

QList<CGIdx> MoleculeInfo::map(const CGName &name) const;
QList<CGIdx> MoleculeInfo::map(CGIdx idx) const;
QList<CGIdx> MoleculeInfo::map(const CGID &cgid) const;

QList<AtomIdx> MoleculeInfo::map(const AtomName &name) const;
QList<AtomIdx> MoleculeInfo::map(AtomNum num) const;
QList<AtomIdx> MoleculeInfo::map(AtomIdx idx) const;
QList<AtomIdx> MoleculeInfo::map(const AtomID &atomid) const;

void MoleculeInfo::squeeze(const MoleculeInfo &other) const;

void MoleculeInfo::assertCompatibleWith(const AtomSelection &selected_atoms) const;
void MoleculeInfo::assertCompatibleWith(const MoleculeView &molecule) const;

void MoleculeInfo::assertContains(AtomIdx atomidx) const;
void MoleculeInfo::assertContains(CGIdx cgidx) const;
void MoleculeInfo::assertContains(ResIdx residx) const;
void MoleculeInfo::assertContains(ChainIdx chainidx) const;
void MoleculeInfo::assertContains(SegIdx segidx) const;

void MoleculeInfo::assertEqualTo(const MoleculeInfo &other) const;



