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

#include "moleculeview.h"
#include "atomselection.h"
#include "atom.h"
#include "cutgroup.h"
#include "residue.h"
#include "chain.h"
#include "segment.h"
#include "selector.hpp"
#include "molecule.h"

#include "SireBase/errors.h"
#include "SireError/errors.h"
#include "SireMol/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>

using namespace SireStream;
using namespace SireBase;
using namespace SireMol;

RegisterMetaType<MoleculeView> r_molview( MAGIC_ONLY,
                                          "SireMol::MoleculeView" );

/** Serialise to a binary datastream */
QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds, 
                                       const MoleculeView &molview)
{
    writeHeader(ds, r_molview, 2);

    SharedDataStream sds(ds);
    
    sds << molview.d
        << static_cast<const Property&>(molview);

    return ds;
}

/** Deserialise from a binary datastream */
QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds,
                                       MoleculeView &molview)
{
    VersionID v = readHeader(ds, r_molview);

    if (v == 2)
    {
        SharedDataStream sds(ds);
        
        sds >> molview.d >> static_cast<Property&>(molview);
    }
    else if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> molview.d;
    }
    else
        throw version_error(v, "1", r_molview, CODELOC);

    return ds;
}

/** Null constructor */
MoleculeView::MoleculeView() : Property(), d( MoleculeData::null() )
{}

/** Construct a view of the molecule whose data is in 'moldata' */
MoleculeView::MoleculeView(const MoleculeData &moldata)
             : Property(), d(moldata)
{}

/** Copy constructor */
MoleculeView::MoleculeView(const MoleculeView &other)
             : Property(other), d(other.d)
{}

/** Destructor */
MoleculeView::~MoleculeView()
{}

/** Copy assignment operator */
MoleculeView& MoleculeView::operator=(const MoleculeView &other)
{
    Property::operator=(other);
    d = other.d;
    return *this;
}

/** Comparison operator */
bool MoleculeView::operator==(const MoleculeView &other) const
{
    return d == other.d or *d == *(other.d);
}

/** Comparison operator */
bool MoleculeView::operator!=(const MoleculeView &other) const
{
    return d != other.d and *d != *(other.d);
}

/** Return whether or not this molecule view is null */
bool MoleculeView::isNull() const
{
    return *d == *MoleculeData::null();
}

/** Return whether or not this view is of the same molecule as 'other'
    (albeit perhaps a different version of the molecule) */
bool MoleculeView::isSameMolecule(const MoleculeData &other) const
{
    return d->number() == other.number();
}

/** Return whether or not this view is of the same molecule as 'other'
    (albeit perhaps a different version of the molecule) */
bool MoleculeView::isSameMolecule(const MoleculeView &other) const
{
    return this->isSameMolecule(other.data());
}

/** Assert that this view is looking at the molecule whose data is 
    in 'other' (albeit perhaps a different version of that molecule)
    
    \throw SireError::incompatible_error
*/
void MoleculeView::assertSameMolecule(const MoleculeData &other) const
{
    if (d->number() != other.number())
        //these are different molecules!
        throw SireError::incompatible_error( QObject::tr(
            "The molecules \"%1\", number %2, and \"%3\", number %3, "
            "are different, and therefore incompatible.")
                .arg(d->name()).arg(d->number())
                .arg(other.name()).arg(other.number()),
                    CODELOC );
}

/** Assert that this is a view of the same molecule as 'other'
    (albeit at a different version)
    
    \throw SireError::incompatible_error
*/
void MoleculeView::assertSameMolecule(const MoleculeView &other) const
{
    this->assertSameMolecule(other.data());
}

/** Update this view with a new version of the molecule. You 
    can only update the molecule if it has the same layout UID
    (so same atoms, residues, cutgroups etc.)
    
    \throw SireError::incompatible_error
*/
void MoleculeView::update(const MoleculeData &moldata)
{
    this->assertSameMolecule(moldata);
    d->info().assertEqualTo(moldata.info());
    
    d = moldata;
}

/** Return the type of the property at key 'key'

    \throw SireBase::missing_property 
*/
const char* MoleculeView::propertyType(const PropertyName &key) const
{
    return d->property(key).what();
}

/** Return the type of the metadata at metakey 'metakey' 

    \throw SireBase::missing_property
*/
const char* MoleculeView::metadataType(const PropertyName &metakey) const
{
    return d->metadata(metakey).what();
}

/** Return the type of the metadata at metakey 'metakey' 
    for the property at key 'key'
    
    \throw SireBase::missing_property
*/
const char* MoleculeView::metadataType(const PropertyName &key,
                                       const PropertyName &metakey) const
{
    return d->metadata(key, metakey).what();
}

/** Assert that this view contains the atom at index 'atomidx'

    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
void MoleculeView::assertContains(AtomIdx atomidx) const
{
    if (not this->selection().selected(atomidx))
        throw SireMol::missing_atom( QObject::tr(
            "This view of the molecule \"%1\" (%2) does not "
            "contain the atom at index %3.")
                .arg(d->name()).arg(d->number())
                .arg(atomidx), CODELOC );
}

/** Assert that this contains a property at key 'key' 

    \throw SireBase::missing_property
*/
void MoleculeView::assertHasProperty(const PropertyName &key) const
{
    if (not this->hasProperty(key))
        throw SireBase::missing_property( QObject::tr(
            "This view of the molecule \"%1\" (view type %2) "
            "does not have a valid property at key \"%3\".")
                .arg(d->name())
                .arg(this->what())
                .arg(key.toString()), CODELOC );
}

/** Assert that this contains some metadata at metakey 'metakey' 

    \throw SireBase::missing_property
*/
void MoleculeView::assertHasMetadata(const PropertyName &metakey) const
{
    if (not this->hasMetadata(metakey))
        throw SireBase::missing_property( QObject::tr(
            "This view of the molecule \"%1\" (view type %2) "
            "does not have some valid metadata at metakey \"%3\".")
                .arg(d->name())
                .arg(this->what())
                .arg(metakey.toString()), CODELOC );
}

/** Assert that this contains some metadata at metakey 'metakey' 
    for the property at key 'key'

    \throw SireBase::missing_property
*/
void MoleculeView::assertHasMetadata(const PropertyName &key,
                                     const PropertyName &metakey) const
{
    if (not this->hasMetadata(key,metakey))
        throw SireBase::missing_property( QObject::tr(
            "This view of the molecule \"%1\" (view type %2) "
            "does not have some valid metadata at metakey \"%3\" "
            "for the property at key \"%4\".")
                .arg(d->name())
                .arg(this->what())
                .arg(metakey.toString())
                .arg(key.toString()), CODELOC );
}

/** Return the atom in this view that matches the ID 'atomid'

    \throw SireMol::missing_atom
    \throw SireError::invalid_index
    \throw SireMol::duplicate_atom
*/
Atom MoleculeView::atom(const AtomID &atomid, const PropertyMap &map) const
{
    return atomid.selectFrom(*this, map);
}

/** Return the atoms from this view that match the ID 'atomid'

    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
Selector<Atom> MoleculeView::atoms(const AtomID &atomid,
                                   const PropertyMap &map) const
{
    return atomid.selectAllFrom(*this, map);
}

/** Return this view as a Atom - this will only work if
    this view contains only a single atom
    
    \throw SireMol::duplicate_atom
*/
Atom MoleculeView::atom() const
{
    QVector<AtomIdx> selected_atoms = this->selection().selectedAtoms();
    
    if (selected_atoms.isEmpty())
        throw SireMol::missing_atom( QObject::tr(
                "This view does not contain any atoms."), CODELOC );
    
    else if (selected_atoms.count() != 1)
        throw SireMol::duplicate_atom( QObject::tr(
                "Cannot convert this view (%1) into an Atom as "
                "we can only do this is just one atom is selected. "
                "These atoms are selected ; %2.")
                    .arg(this->toString(), Sire::toString(selected_atoms)), 
                        CODELOC );
    
    return Atom(this->data(), selected_atoms.at(0));
}

/** Return all of the atoms in this view 

    \throw SireMol::missing_atom
*/
Selector<Atom> MoleculeView::atoms() const
{
    AtomSelection selected_atoms = this->selection();
    
    if (selected_atoms.selectedNone())
        throw SireMol::missing_atom( QObject::tr(
                "No atoms are available in this view (%1).")
                    .arg(this->toString()), CODELOC );

    return Selector<Atom>(this->data(), this->selection());
}

/** Return the CutGroup whose atoms are in this view that matches
    the ID in 'cgid'
    
    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
    \throw SireMol::duplicate_cutgroup
*/
CutGroup MoleculeView::cutGroup(const CGID &cgid, const PropertyMap &map) const
{
    return cgid.selectFrom(*this, map);
}

/** Return the CutGroups whose atoms are in this view that match
    the ID in 'cgid'
    
    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
*/
Selector<CutGroup> MoleculeView::cutGroups(const CGID &cgid,
                                           const PropertyMap &map) const
{
    return cgid.selectAllFrom(*this, map);
}

/** Return the CutGroup that contains the atom(s) in this view

    \throw SireMol::missing_cutgroup
    \throw SireMol::duplicate_cutgroup
*/
CutGroup MoleculeView::cutGroup() const
{
    QList<CGIdx> selected_cgs = this->selection().selectedCutGroups();
    
    if (selected_cgs.isEmpty())
        throw SireMol::missing_cutgroup( QObject::tr(
                "This view does not contain any CutGroups."), CODELOC );
    
    else if (selected_cgs.count() != 1)
        throw SireMol::duplicate_cutgroup( QObject::tr(
                "Cannot convert this view (%1) into a CutGroup as "
                "we can only do this is just one CutGroup is selected. "
                "These CutGroups are selected ; %2.")
                    .arg(this->toString(), Sire::toString(selected_cgs)), 
                        CODELOC );
    
    return CutGroup(this->data(), selected_cgs.at(0));
}

/** Return all of the CutGroups that are involved in this view

    \throw SireMol::missing_cutgroup
*/
Selector<CutGroup> MoleculeView::cutGroups() const
{
    QList<CGIdx> selected_cgs = this->selection().selectedCutGroups();
    
    if (selected_cgs.isEmpty())
        throw SireMol::missing_cutgroup( QObject::tr(
                "This view does not contain any CutGroups."), CODELOC );

    return Selector<CutGroup>(this->data(), selected_cgs);
}

/** Return the residue from this view that matches the ID 'resid'

    \throw SireMol::missing_residue
    \throw SireError::invalid_index
    \throw SireMol::duplicate_residue
*/
Residue MoleculeView::residue(const ResID &resid, const PropertyMap &map) const
{
    return resid.selectFrom(*this, map);
}

/** Return the residues from this view that match the ID 'resid'

    \throw SireMol::missing_residue
    \throw SireError::invalid_index
*/
Selector<Residue> MoleculeView::residues(const ResID &resid,
                                         const PropertyMap &map) const
{
    return resid.selectAllFrom(*this, map);
}

/** Return the residue that is part of this view

    \throw SireMol::missing_residue
    \throw SireMol::duplicate_residue
*/
Residue MoleculeView::residue() const
{
    QList<ResIdx> selected_res = this->selection().selectedResidues();
    
    if (selected_res.isEmpty())
        throw SireMol::missing_residue( QObject::tr(
                "This view does not contain any residues."), CODELOC );
    
    else if (selected_res.count() != 1)
        throw SireMol::duplicate_residue( QObject::tr(
                "Cannot convert this view (%1) into a Residue as "
                "we can only do this is just one Residue is selected. "
                "These Residues are selected ; %2.")
                    .arg(this->toString(), Sire::toString(selected_res)), 
                        CODELOC );
    
    return Residue(this->data(), selected_res.at(0));
}

/** Return all of the residues that are involved with this view

    \throw SireMol::missing_residue
*/
Selector<Residue> MoleculeView::residues() const
{
    QList<ResIdx> selected_res = this->selection().selectedResidues();
    
    if (selected_res.isEmpty())
        throw SireMol::missing_residue( QObject::tr(
                "This view does not contain any residues."), CODELOC );
    
    return Selector<Residue>(this->data(), selected_res);
}

/** Return the chain that is involved with this view that matches
    the ID 'chainid'

    \throw SireMol::missing_chain
    \throw SireError::invalid_index
    \throw SireMol::duplicate_chain
*/
Chain MoleculeView::chain(const ChainID &chainid, const PropertyMap &map) const
{
    return chainid.selectFrom(*this, map);
}

/** Return the chains that are involved with this view that match
    the ID 'chainid'

    \throw SireMol::missing_chain
    \throw SireError::invalid_index
    \throw SireMol::duplicate_chain
*/
Selector<Chain> MoleculeView::chains(const ChainID &chainid,
                                     const PropertyMap &map) const
{
    return chainid.selectAllFrom(*this, map);
}

/** Return the chain that is involved with this view

    \throw SireMol::missing_chain
    \throw SireMol::duplicate_chain
*/
Chain MoleculeView::chain() const
{
    QList<ChainIdx> selected_chn = this->selection().selectedChains();
    
    if (selected_chn.isEmpty())
        throw SireMol::missing_chain( QObject::tr(
                "This view does not contain any chains."), CODELOC );
    
    else if (selected_chn.count() != 1)
        throw SireMol::duplicate_chain( QObject::tr(
                "Cannot convert this view (%1) into a Chain as "
                "we can only do this is just one Chain is selected. "
                "These Chains are selected ; %2.")
                    .arg(this->toString(), Sire::toString(selected_chn)), 
                        CODELOC );
    
    return Chain(this->data(), selected_chn.at(0));
}

/** Return the chains that are involved with this view

    \throw SireMol::missing_chain
*/
Selector<Chain> MoleculeView::chains() const
{
    QList<ChainIdx> selected_chn = this->selection().selectedChains();
    
    if (selected_chn.isEmpty())
        throw SireMol::missing_chain( QObject::tr(
                "This view does not contain any chains."), CODELOC );
    
    return Selector<Chain>(this->data(), selected_chn);
}

/** Return the segment that is involved with this view that matches
    the ID 'segid'

    \throw SireMol::missing_segment
    \throw SireError::invalid_index
    \throw SireMol::duplicate_segment
*/
Segment MoleculeView::segment(const SegID &segid, const PropertyMap &map) const
{
    return segid.selectFrom(*this, map);
}

/** Return the segments that are involved with this view that match
    the ID 'segid'

    \throw SireMol::missing_segment
    \throw SireError::invalid_index
*/
Selector<Segment> MoleculeView::segments(const SegID &segid, const PropertyMap &map) const
{
    return segid.selectAllFrom(*this, map);
}

/** Return the segment that is involved with this view

    \throw SireMol::missing_segment
    \throw SireMol::duplicate_segment
*/
Segment MoleculeView::segment() const
{
    QList<SegIdx> selected_seg = this->selection().selectedSegments();
    
    if (selected_seg.isEmpty())
        throw SireMol::missing_segment( QObject::tr(
                "This view does not contain any segments."), CODELOC );
    
    else if (selected_seg.count() != 1)
        throw SireMol::duplicate_segment( QObject::tr(
                "Cannot convert this view (%1) into a Segment as "
                "we can only do this is just one Segment is selected. "
                "These Segments are selected ; %2.")
                    .arg(this->toString(), Sire::toString(selected_seg)), 
                        CODELOC );
    
    return Segment(this->data(), selected_seg.at(0));
}

/** Return the segments that are involved with this view

    \throw SireMol::missing_segment
*/
Selector<Segment> MoleculeView::segments() const
{
    QList<SegIdx> selected_seg = this->selection().selectedSegments();
    
    if (selected_seg.isEmpty())
        throw SireMol::missing_segment( QObject::tr(
                "This view does not contain any segments."), CODELOC );
    
    return Selector<Segment>(this->data(), selected_seg);
}

/** Return the molecule involved with this view */
Molecule MoleculeView::molecule() const
{
    return Molecule(this->data());
}

/** Return the CutGroup whose atoms are in this view that matches
    the ID in 'cgid'
    
    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
    \throw SireMol::duplicate_cutgroup
*/
CutGroup MoleculeView::select(const CGID &cgid, const PropertyMap &map) const
{
    return this->cutGroup(cgid, map);
}

/** Return the residue from this view that matches the ID 'resid'

    \throw SireMol::missing_residue
    \throw SireError::invalid_index
    \throw SireMol::duplicate_residue
*/
Residue MoleculeView::select(const ResID &resid, const PropertyMap &map) const
{
    return this->residue(resid, map);
}

/** Return the chain that is involved with this view that matches
    the ID 'chainid'

    \throw SireMol::missing_chain
    \throw SireError::invalid_index
    \throw SireMol::duplicate_chain
*/
Chain MoleculeView::select(const ChainID &chainid, const PropertyMap &map) const
{
    return this->chain(chainid, map);
}

/** Return the segment that is involved with this view that matches
    the ID 'segid'

    \throw SireMol::missing_segment
    \throw SireError::invalid_index
    \throw SireMol::duplicate_segment
*/
Segment MoleculeView::select(const SegID &segid, const PropertyMap &map) const
{
    return this->segment(segid, map);
}

/** Return the atom in this view that matches the ID 'atomid'

    \throw SireMol::missing_atom
    \throw SireError::invalid_index
    \throw SireMol::duplicate_atom
*/
Atom MoleculeView::select(const AtomID &atomid, const PropertyMap &map) const
{
    return this->atom(atomid, map);
}

/** Return the atoms from this view that match the ID 'atomid'

    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
Selector<Atom> MoleculeView::selectAll(const AtomID &atomid,
                                       const PropertyMap &map) const
{
    return this->atoms(atomid, map);
}

/** Return all of the atoms in this view 

    \throw SireMol::missing_atom
*/
Selector<Atom> MoleculeView::selectAll() const
{
    return this->atoms();
}

/** Return all of the atoms in this view 

    \throw SireMol::missing_atom
*/
Selector<Atom> MoleculeView::selectAllAtoms() const
{
    return this->atoms();
}

/** Return the CutGroups whose atoms are in this view that match
    the ID in 'cgid'
    
    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
*/
Selector<CutGroup> MoleculeView::selectAll(const CGID &cgid,
                                           const PropertyMap &map) const
{
    return this->cutGroups(cgid, map);
}

/** Return all of the CutGroups that are involved in this view

    \throw SireMol::missing_cutgroup
*/
Selector<CutGroup> MoleculeView::selectAllCutGroups() const
{
    return this->cutGroups();
}

/** Return the residues from this view that match the ID 'resid'

    \throw SireMol::missing_residue
    \throw SireError::invalid_index
*/
Selector<Residue> MoleculeView::selectAll(const ResID &resid,
                                          const PropertyMap &map) const
{
    return this->residues(resid, map);
}

/** Return all of the residues that are involved with this view

    \throw SireMol::missing_residue
*/
Selector<Residue> MoleculeView::selectAllResidues() const
{
    return this->residues();
}

/** Return the chains that are involved with this view that match
    the ID 'chainid'

    \throw SireMol::missing_chain
    \throw SireError::invalid_index
    \throw SireMol::duplicate_chain
*/
Selector<Chain> MoleculeView::selectAll(const ChainID &chainid,
                                        const PropertyMap &map) const
{
    return this->chains(chainid, map);
}

/** Return the chains that are involved with this view

    \throw SireMol::missing_chain
*/
Selector<Chain> MoleculeView::selectAllChains() const
{
    return this->chains();
}

/** Return the segments that are involved with this view that match
    the ID 'segid'

    \throw SireMol::missing_segment
    \throw SireError::invalid_index
*/
Selector<Segment> MoleculeView::selectAll(const SegID &segid,
                                          const PropertyMap &map) const
{
    return this->segments(segid, map);
}

/** Return the segments that are involved with this view

    \throw SireMol::missing_segment
*/
Selector<Segment> MoleculeView::selectAllSegments() const
{
    return this->segments();
}
