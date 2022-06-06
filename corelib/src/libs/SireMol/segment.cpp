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

#include "segment.h"

#include "atom.h"
#include "molecule.h"

#include "mover.hpp"
#include "selector.hpp"
#include "evaluator.h"
#include "segeditor.h"

#include "groupatomids.h"

#include "SireBase/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "mover_metaid.h"

using namespace SireMol;
using namespace SireStream;

///////
/////// Implementation of SegProp
///////

static const RegisterMetaType<SegProp> r_segprop(MAGIC_ONLY,
                                                 "SireMol::SegProp");

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const SegProp &segprop)
{
    writeHeader(ds, r_segprop, 1)
         << static_cast<const MolViewProperty&>(segprop);

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, SegProp &segprop)
{
    VersionID v = readHeader(ds, r_segprop);

    if (v == 1)
    {
        ds >> static_cast<MolViewProperty&>(segprop);
    }
    else
        throw version_error(v, "1", r_segprop, CODELOC);

    return ds;
}

SegProp::SegProp() : MolViewProperty()
{}

SegProp::SegProp(const SegProp &other) : MolViewProperty(other)
{}

SegProp::~SegProp()
{}

///////
/////// Implementation of Segment
///////

RegisterMetaType<Segment> r_seg;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const Segment &seg)
{
    writeHeader(ds, r_seg, 1);

    SharedDataStream sds(ds);

    sds << seg.segidx << static_cast<const MoleculeView&>(seg);

    return ds;
}

/** Deserialise from a binary datastream */
QDataStream &operator>>(QDataStream &ds, Segment &seg)
{
    VersionID v = readHeader(ds, r_seg);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> seg.segidx >> static_cast<MoleculeView&>(seg);

        seg.selected_atoms = AtomSelection(seg.data());
        seg.selected_atoms.selectOnly(seg.segidx);
    }
    else
        throw version_error(v, "1", r_seg, CODELOC);

    return ds;
}

/** Null constructor */
Segment::Segment() : ConcreteProperty<Segment,MoleculeView>(), segidx( SegIdx::null() )
{}

/** Construct the Segment at ID 'cgid' in the molecule whose data
    is in 'moldata'

    \throw SireMol::missing_Segment
    \throw SireMol::duplicate_Segment
    \throw SireError::invalid_index
*/
Segment::Segment(const MoleculeData &moldata, const SegID &segid)
      : ConcreteProperty<Segment,MoleculeView>(moldata),
        segidx( moldata.info().segIdx(segid) )
{
    selected_atoms = AtomSelection(moldata);
    selected_atoms.selectOnly(segidx);
}

/** Copy constructor */
Segment::Segment(const Segment &other)
        : ConcreteProperty<Segment,MoleculeView>(other), segidx(other.segidx),
          selected_atoms(other.selected_atoms)
{}

/** Destructor */
Segment::~Segment()
{}

/** Copy assignment operator */
Segment& Segment::operator=(const Segment &other)
{
    MoleculeView::operator=(other);
    segidx = other.segidx;
    selected_atoms = other.selected_atoms;
    return *this;
}

/** Comparison operator */
bool Segment::operator==(const Segment &other) const
{
    return segidx == other.segidx and
           MoleculeView::operator==(other);
}

/** Comparison operator */
bool Segment::operator!=(const Segment &other) const
{
    return segidx != other.segidx or
           MoleculeView::operator!=(other);
}

/** Return a string representation of this segment */
QString Segment::toString() const
{
    return QObject::tr( "Segment( %1 num_atoms=%2 )" )
                .arg( this->name() ).arg(this->nAtoms());
}

/** Return whether or not this segment is empty */
bool Segment::isEmpty() const
{
    return selected_atoms.selectedNone();
}

/** Return whether or not this segment contains the entire molecule */
bool Segment::selectedAll() const
{
    return selected_atoms.selectedAll();
}

MolViewPtr Segment::toSelector() const
{
    return MolViewPtr( Selector<Segment>(*this) );
}

/** Return the atoms that are in this Segment */
AtomSelection Segment::selection() const
{
    return selected_atoms;
}

/** Update this segment with the passed molecule data.

    \throw SireError::incompatible_error
*/
void Segment::update(const MoleculeData &moldata)
{
    //check that the new data is compatible (has same molecule
    //number and info ID number)
    if (d->number() != moldata.number() or
        d->info().UID() != moldata.info().UID())
    {
        throw SireError::incompatible_error( QObject::tr(
            "You can only update a segment with the molecule data "
            "for the same molecule (same molecule number) and that "
            "has a .info() object that has the same UID. You are "
            "trying to update segment %1 in molecule %2 with UID %3 "
            "with molecule %4 with UID %5.")
                .arg(segidx).arg(d->number()).arg(d->info().UID().toString())
                .arg(moldata.number()).arg(moldata.info().UID().toString()),
                    CODELOC );
    }

    d = moldata;
}

/** Return the name of this Segment */
const SegName& Segment::name() const
{
    return d->info().name(segidx);
}

/** Return the index of this Segment in the molecule */
SegIdx Segment::index() const
{
    return segidx;
}

/** Return the number of this segment (same as its index) */
SegIdx Segment::number() const
{
    return segidx;
}

/** Return an object that can move a copy of this Segment */
Mover<Segment> Segment::move() const
{
    return Mover<Segment>(*this);
}

/** Return an evaluator that can evaluate properties
    of this Segment */
Evaluator Segment::evaluate() const
{
    return Evaluator(*this);
}

/** Return an editor that can edit this Segment */
SegEditor Segment::edit() const
{
    return SegEditor(*this);
}

/** Return a selector that can be used to change the selection
    of segments from the molecule */
Selector<Segment> Segment::selector() const
{
    return Selector<Segment>(*this);
}

/** Return the number of atoms in this Segment */
int Segment::nAtoms() const
{
    return d->info().nAtoms(segidx);
}

/** Return the indicies of the atoms in this segment, in the
    order that they appear in this segment */
const QList<AtomIdx>& Segment::atomIdxs() const
{
    return d->info().getAtomsIn(segidx);
}

/** Return whether or not this segment contains the atom
    at index 'atomidx' */
bool Segment::contains(AtomIdx atomidx) const
{
    return d->info().contains(segidx, atomidx);
}

/** Return whether or not this segment contains all of
    the atoms identified by the ID 'atomid' */
bool Segment::contains(const AtomID &atomid) const
{
    return d->info().contains(segidx, atomid);
}

/** Return whether or not this segment contains some of
    the atoms identified by the ID 'atomid' */
bool Segment::intersects(const AtomID &atomid) const
{
    return d->info().intersects(segidx, atomid);
}

/** Return the specified property as a QVariant */
QVariant Segment::propertyAsVariant(const PropertyName &key) const
{
    const Property &property = d->property(key);
    return property.asA<SegProp>().getAsVariant(segidx);
}

/** Return the specified property as a PropertyPtr */
PropertyPtr Segment::propertyAsProperty(const PropertyName &key) const
{
    const Property &property = d->property(key);
    return property.asA<SegProp>().getAsProperty(segidx);
}

/** Return whether or not there is a SegProperty at key 'key' */
bool Segment::hasProperty(const PropertyName &key) const
{
    return d->hasPropertyOfType<SegProp>(key);
}

/** Return whether or not there is a SegProperty at metakey 'metakey' */
bool Segment::hasMetadata(const PropertyName &metakey) const
{
    return d->hasMetadataOfType<SegProp>(metakey);
}

/** Return whether the metadata at metakey 'metakey' for the property
    at key 'key' is a SegProperty

    \throw SireBase::missing_property
*/
bool Segment::hasMetadata(const PropertyName &key,
                       const PropertyName &metakey) const
{
    return d->hasMetadataOfType<SegProp>(key, metakey);
}

/** Return the keys of all SegProperty properties */
QStringList Segment::propertyKeys() const
{
    return d->properties().propertyKeysOfType<SegProp>();
}

/** Return the metakeys of all SegProperty metadata */
QStringList Segment::metadataKeys() const
{
    return d->properties().metadataKeysOfType<SegProp>();
}

/** Return the metakeys of all SegProperty metadata for
    the property at key 'key'

    \throw SireBase::missing_property
*/
QStringList Segment::metadataKeys(const PropertyName &key) const
{
    return d->properties().metadataKeysOfType<SegProp>(key);
}

/** Assert that this segment has an SegProperty at key 'key'

    \throw SireBase::missing_property
*/
void Segment::assertContainsProperty(const PropertyName &key) const
{
    if (not this->hasProperty(key))
        throw SireBase::missing_property( QObject::tr(
            "There is no SegProperty at key '%1' for this segment.")
                .arg(key.toString()), CODELOC );
}

/** Assert that this segment has an SegProperty piece of metadata
    at metakey 'metakey'

    \throw SireBase::missing_property
*/
void Segment::assertContainsMetadata(const PropertyName &metakey) const
{
    if (not this->hasMetadata(metakey))
        throw SireBase::missing_property( QObject::tr(
            "There is no SegProperty metadata at metakey '%1' for "
            "this segment.")
                .arg(metakey.toString()), CODELOC );
}

/** Assert that the property at key 'key' has an SegProperty
    piece of metadata at metakey 'metakey'

    \throw SireBase::missing_property
*/
void Segment::assertContainsMetadata(const PropertyName &key,
                                     const PropertyName &metakey) const
{
    if (not this->hasMetadata(key, metakey))
        throw SireBase::missing_property( QObject::tr(
            "There is no SegProperty metadata at metakey '%1' "
            "for the property at key '%2' for this segment.")
                .arg(metakey.toString(), key.toString()), CODELOC );
}

const char* Segment::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Segment>() );
}

bool SireMol::detail::has_property(const Segment*, const MoleculeData &moldata,
                                                  const PropertyName &key)
{
    return moldata.hasPropertyOfType<SegProp>(key);
}

bool SireMol::detail::has_metadata(const Segment*, const MoleculeData &moldata,
                                                  const PropertyName &metakey)
{
    return moldata.hasMetadataOfType<SegProp>(metakey);
}

bool SireMol::detail::has_metadata(const Segment*, const MoleculeData &moldata,
                                                  const PropertyName &key, const PropertyName &metakey)
{
    return moldata.hasMetadataOfType<SegProp>(key, metakey);
}

namespace SireMol
{
    /////// explicitly instantiate the templates
    template class Mover<Segment>;
    template class Selector<Segment>;

    template class Mover< Selector<Segment> >;
}
