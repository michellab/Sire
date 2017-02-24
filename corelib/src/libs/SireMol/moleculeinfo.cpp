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

const CGName& MoleculeInfo::name(const CGID &cgid) const;
const CGName& MoleculeInfo::name(CGIdx cgidx) const;

const AtomName& MoleculeInfo::name(const AtomID &atomid) const;
const AtomName& MoleculeInfo::name(AtomIdx atomidx) const;

ResNum MoleculeInfo::number(const ResID &resid) const;
ResNum MoleculeInfo::number(ResIdx residx) const;

AtomNum MoleculeInfo::number(const AtomID &atomid) const;
AtomNum MoleculeInfo::number(AtomIdx atomidx) const;

MoleculeInfo MoleculeInfo::rename(AtomIdx atomidx, const AtomName &newname) const;
MoleculeInfo MoleculeInfo::renumber(AtomIdx atomidx, const AtomNum &newnum) const;

MoleculeInfo MoleculeInfo::rename(ResIdx residx, const ResName &newname) const;
MoleculeInfo MoleculeInfo::renumber(ResIdx residx, const ResNum &newnum) const;

MoleculeInfo MoleculeInfo::rename(CGIdx cgidx, const CGName &newname) const;
MoleculeInfo MoleculeInfo::rename(ChainIdx chainidx, const ChainName &newname) const;
MoleculeInfo MoleculeInfo::rename(SegIdx segidx, const SegName &newname) const;

const CGAtomIdx& MoleculeInfo::cgAtomIdx(AtomIdx atomidx) const;
const CGAtomIdx& MoleculeInfo::cgAtomIdx(const AtomID &atomid) const;

QVector<CGAtomIdx> MoleculeInfo::cgAtomIdxs(AtomIdx atomidx) const;
QVector<CGAtomIdx> MoleculeInfo::cgAtomIdxs(CGIdx cgidx) const;
QVector<CGAtomIdx> MoleculeInfo::cgAtomIdxs(ResIdx residx) const;
QVector<CGAtomIdx> MoleculeInfo::cgAtomIdxs(ChainIdx chainidx) const;
QVector<CGAtomIdx> MoleculeInfo::cgAtomIdxs(SegIdx segidx) const;

QVector<CGAtomIdx> MoleculeInfo::cgAtomIdxs(const AtomID &atomid) const;
QVector<CGAtomIdx> MoleculeInfo::cgAtomIdxs(const CGID &cgid) const;
QVector<CGAtomIdx> MoleculeInfo::cgAtomIdxs(const ResID &resid) const;
QVector<CGAtomIdx> MoleculeInfo::cgAtomIdxs(const ChainID &chainid) const;
QVector<CGAtomIdx> MoleculeInfo::cgAtomIdxs(const SegID &segid) const;

AtomIdx MoleculeInfo::atomIdx(const AtomID &atomid) const;
AtomIdx MoleculeInfo::atomIdx(const CGAtomIdx &cgatomidx) const;

ResIdx MoleculeInfo::resIdx(const ResID &resid) const;
ChainIdx MoleculeInfo::chainIdx(const ChainID &chainid) const;
SegIdx MoleculeInfo::segIdx(const SegID &segid) const;
CGIdx MoleculeInfo::cgIdx(const CGID &cgid) const;

QList<SegIdx> MoleculeInfo::getSegments() const;
QList<CGIdx> MoleculeInfo::getCutGroups() const;
QList<ChainIdx> MoleculeInfo::getChains() const;
QList<ResIdx> MoleculeInfo::getResidues() const;

const QList<ResIdx>& MoleculeInfo::getResiduesIn(ChainIdx chainidx) const;
QList<ResIdx> MoleculeInfo::getResiduesIn(const ChainID &chainid) const;

QList<AtomIdx> MoleculeInfo::getAtoms() const;

AtomIdx MoleculeInfo::getAtom(CGIdx cgidx, int i) const;
AtomIdx MoleculeInfo::getAtom(ResIdx residx, int i) const;
AtomIdx MoleculeInfo::getAtom(ChainIdx chainidx, int i) const;
AtomIdx MoleculeInfo::getAtom(SegIdx segidx, int i) const;

ResIdx MoleculeInfo::getResidue(ChainIdx chainidx, int i) const;

const QList<AtomIdx>& MoleculeInfo::getAtomsIn(ResIdx residx) const;
QList<AtomIdx> MoleculeInfo::getAtomsIn(const ResID &resid) const;
QList<AtomIdx> MoleculeInfo::getAtomsIn(ResIdx residx, const AtomName &name) const;
QList<AtomIdx> MoleculeInfo::getAtomsIn(const ResID &resid,
                          const AtomName &atomname) const;

QList<AtomIdx> MoleculeInfo::getAtomsIn(ChainIdx chainidx) const;
QList<AtomIdx> MoleculeInfo::getAtomsIn(const ChainID &chainid) const;
QList<AtomIdx> MoleculeInfo::getAtomsIn(ChainIdx chainidx,
                          const AtomName &atomname) const;
QList<AtomIdx> MoleculeInfo::getAtomsIn(const ChainID &chainid,
                          const AtomName &atomname) const;
                          
const QList<AtomIdx>& MoleculeInfo::getAtomsIn(CGIdx cgidx) const;
QList<AtomIdx> MoleculeInfo::getAtomsIn(const CGID &cgid) const;

const QList<AtomIdx>& MoleculeInfo::getAtomsIn(SegIdx segidx) const;
QList<AtomIdx> MoleculeInfo::getAtomsIn(const SegID &segid) const;

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



