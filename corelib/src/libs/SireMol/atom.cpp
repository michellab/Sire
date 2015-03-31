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

#include "atom.h"
#include "atomeditor.h"
#include "atomproperty.hpp"

#include "atomcharges.h"

#include "mover.hpp"
#include "evaluator.h"

#include "cutgroup.h"
#include "residue.h"
#include "chain.h"
#include "segment.h"
#include "molecule.h"
#include "selector.hpp"

#include "SireBase/errors.h"
#include "SireMol/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "mover_metaid.h"

#include <QDebug>

using namespace SireMol;
using namespace SireStream;

///////
/////// Implementation of Atom
///////

static const RegisterMetaType<Atom> r_atom;

/** Serialise to a binary datastream */
QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds, const Atom &atom)
{
    writeHeader(ds, r_atom, 1);

    ds << atom.atomidx << static_cast<const MoleculeView&>(atom);

    return ds;
}

/** Deserialise from a binary datastream */
QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds, Atom &atom)
{
    VersionID v = readHeader(ds, r_atom);

    if (v == 1)
    {
        ds >> atom.atomidx >> static_cast<MoleculeView&>(atom);
    }
    else
        throw version_error(v, "1", r_atom, CODELOC);

    return ds;
}

void SireMol::detail::assertSameSize(Atom*, int nats, int nprops)
{
    if (nats != nprops)
        throw SireError::incompatible_error( QObject::tr(
            "The number of supplied properties (%1) is not the same "
            "as the number of atoms (%2).")
                .arg(nprops).arg(nats), CODELOC );
}

/** Null constructor */
Atom::Atom() : ConcreteProperty<Atom,MoleculeView>(), atomidx( Index::null() )
{}

/** Construct the atom that that is identified by ID 'atomid'
    in the view 'molview' - this atom must be within this view
    
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
Atom::Atom(const MoleculeView &molview, const AtomID &atomid)
     : ConcreteProperty<Atom,MoleculeView>(molview)
{
    atomidx = d->info().atomIdx(atomid);
    molview.assertContains(atomidx);
}

/** Construct the atom that is identified by ID 'atomid'
    in the molecule whose data is in 'moldata'
    
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
Atom::Atom(const MoleculeData &moldata, const AtomID &atomid)
     : ConcreteProperty<Atom,MoleculeView>(moldata), 
       atomidx( moldata.info().atomIdx(atomid) )
{}

/** Copy constructor */
Atom::Atom(const Atom &other)
     : ConcreteProperty<Atom,MoleculeView>(other), atomidx(other.atomidx)
{}

/** Destructor */
Atom::~Atom()
{}

/** Copy assignment operator */
Atom& Atom::operator=(const Atom &other)
{
    MoleculeView::operator=(other);
    atomidx = other.atomidx;
    
    return *this;
}

/** Comparison operator */
bool Atom::operator==(const Atom &other) const
{
    return atomidx == other.atomidx and MoleculeView::operator==(other);
}

/** Comparison operator */
bool Atom::operator!=(const Atom &other) const
{
    return atomidx != other.atomidx or MoleculeView::operator!=(other);
}

/** Return a string representation of this atom */
QString Atom::toString() const
{
    return QObject::tr( "Atom( %1 : %2 )" ).arg(this->name())
                                           .arg(this->number());
}

/** Is this atom empty? */
bool Atom::isEmpty() const
{
    return atomidx.isNull();
}

/** Is this atom a view of the whole (1 atom) molecule? */
bool Atom::selectedAll() const
{
    return (not atomidx.isNull()) and d->info().nAtoms() == 1;
}

/** Return the selected atom! */
AtomSelection Atom::selection() const
{
    AtomSelection selected_atoms(data());
    selected_atoms.selectOnly(atomidx);
    
    return selected_atoms;
}

/** Return the name of the atom */
AtomName Atom::name() const
{
    return d->info().name(atomidx);
}

/** Return the number of the atom */
AtomNum Atom::number() const
{
    return d->info().number(atomidx);
}

/** Return the index number of this atom in the molecule */
AtomIdx Atom::index() const
{
    return atomidx;
}

/** Return the CGAtomIdx of this atom */
const CGAtomIdx& Atom::cgAtomIdx() const
{
    return d->info().cgAtomIdx(atomidx);
}

/** Return a Mover that can be used to move this atom */
Mover<Atom> Atom::move() const
{
    return Mover<Atom>(*this);
}

/** Return an evaluator that can be used to evaluate properties
    of this atom */
Evaluator Atom::evaluate() const
{
    return Evaluator(*this);
}

/** Return an editor that can be used to edit this atom */
AtomEditor Atom::edit() const
{
    return AtomEditor(*this);
}

/** Return a selector that can change the atom selection */
Selector<Atom> Atom::selector() const
{
    return Selector<Atom>(*this);
}

/** Return whether or not this atom is part of a residue */
bool Atom::isWithinResidue() const
{
    return d->info().isWithinResidue(atomidx);
}

/** Return whether or not this atom is part of a chain */
bool Atom::isWithinChain() const
{
    return d->info().isWithinChain(atomidx);
}

/** Return whether or not this atom is part of a segment */
bool Atom::isWithinSegment() const
{
    return d->info().isWithinSegment(atomidx);
}

/** Return the residue that this atom is in 

    \throw SireMol::missing_residue
*/
Residue Atom::residue() const
{
    return Residue(*d, d->info().parentResidue(atomidx));
}

/** Return the chain this atom is in 

    \throw SireMol::missing_chain
*/
Chain Atom::chain() const
{
    return Chain(*d, d->info().parentChain(atomidx));
}

/** Return the segment this atom is in 

    \throw SireMol::missing_segment
*/
Segment Atom::segment() const
{
    return Segment(*d, d->info().parentSegment(atomidx));
}

/** Return the CutGroup this atom is in */
CutGroup Atom::cutGroup() const
{
    return CutGroup(*d, d->info().parentCutGroup(atomidx));
}

/** Return the molecule that contains this atom */
Molecule Atom::molecule() const
{
    return Molecule(*d);
}

/** Update this atom with the passed molecule data.

    \throw SireError::incompatible_error
*/
void Atom::update(const MoleculeData &moldata)
{
    //check that the new data is compatible (has same molecule
    //number and info ID number)
    if (d->number() != moldata.number() or
        d->info().UID() != moldata.info().UID())
    {
        throw SireError::incompatible_error( QObject::tr(
            "You can only update an atom with the molecule data "
            "for the same molecule (same molecule number) and that "
            "has a .info() object that has the same UID. You are "
            "trying to update atom %1 in molecule %2 with UID %3 "
            "with molecule %4 with UID %5.")
                .arg(atomidx).arg(d->number()).arg(d->info().UID().toString())
                .arg(moldata.number()).arg(moldata.info().UID().toString()),
                    CODELOC );
    }
    
    d = moldata;
}

/** Return whether or not there is an AtomProperty at key 'key' */
bool Atom::hasProperty(const PropertyName &key) const
{
    return d->hasPropertyOfType<AtomProp>(key);
}

/** Return whether or not there is an AtomProperty at metakey 'metakey' */
bool Atom::hasMetadata(const PropertyName &metakey) const
{
    return d->hasMetadataOfType<AtomProp>(metakey);
}

/** Return whether the metadata at metakey 'metakey' for the property
    at key 'key' is an AtomProperty
    
    \throw SireBase::missing_property
*/
bool Atom::hasMetadata(const PropertyName &key,
                       const PropertyName &metakey) const
{
    return d->hasMetadataOfType<AtomProp>(key, metakey);
}

/** Return the keys of all AtomProperty properties */
QStringList Atom::propertyKeys() const
{
    return d->properties().propertyKeysOfType<AtomProp>();
}

/** Return the metakeys of all AtomProperty metadata */
QStringList Atom::metadataKeys() const
{
    return d->properties().metadataKeysOfType<AtomProp>();
}

/** Return the metakeys of all AtomProperty metadata for 
    the property at key 'key'
    
    \throw SireBase::missing_property
*/
QStringList Atom::metadataKeys(const PropertyName &key) const
{
    return d->properties().metadataKeysOfType<AtomProp>(key);
}

/** Assert that this atom is the atom at index 'atomidx'

    \throw SireMol::missing_atom
*/
void Atom::assertContains(AtomIdx atom) const
{
    if (atomidx != atom.map(d->info().nAtoms()))
        throw SireMol::missing_atom( QObject::tr(
            "This atom (index %1) is not the atom at index %2.")
                .arg(atomidx).arg(atom), CODELOC );
}

/** Assert that this atom has an AtomProperty at key 'key'

    \throw SireBase::missing_property
*/
void Atom::assertContainsProperty(const PropertyName &key) const
{
    if (not this->hasProperty(key))
        throw SireBase::missing_property( QObject::tr(
            "There is no AtomProperty at key '%1' for this atom.")
                .arg(key.toString()), CODELOC );
}

/** Assert that this atom has an AtomProperty piece of metadata
    at metakey 'metakey'
    
    \throw SireBase::missing_property
*/
void Atom::assertContainsMetadata(const PropertyName &metakey) const
{
    if (not this->hasMetadata(metakey))
        throw SireBase::missing_property( QObject::tr(
            "There is no AtomProperty metadata at metakey '%1' for "
            "this atom.")
                .arg(metakey.toString()), CODELOC );
}

/** Assert that the property at key 'key' has an AtomProperty
    piece of metadata at metakey 'metakey'
    
    \throw SireBase::missing_property
*/
void Atom::assertContainsMetadata(const PropertyName &key,
                                  const PropertyName &metakey) const
{
    if (not this->hasMetadata(key, metakey))
        throw SireBase::missing_property( QObject::tr(
            "There is no AtomProperty metadata at metakey '%1' "
            "for the property at key '%2' for this atom.")
                .arg(metakey.toString(), key.toString()), CODELOC );
}

const char* Atom::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Atom>() );
}

namespace SireMol
{
namespace detail
{

bool SIREMOL_EXPORT has_property(const Atom*, const MoleculeData &moldata,
                                 const PropertyName &key)
{
    return moldata.hasPropertyOfType<AtomProp>(key);
}

bool SIREMOL_EXPORT has_metadata(const Atom*, const MoleculeData &moldata,
                                 const PropertyName &metakey)
{
    return moldata.hasMetadataOfType<AtomProp>(metakey);
}

bool SIREMOL_EXPORT has_metadata(const Atom*, const MoleculeData &moldata,
                                 const PropertyName &key, const PropertyName &metakey)
{
    return moldata.hasMetadataOfType<AtomProp>(key, metakey);
}

} // end of namespace detail
} // end of namespace SireMol

///////
/////// Explicitly instantiate the Atom manipulator classes
///////

namespace SireMol
{
    template class Selector<Atom>;
    template class Mover<Atom>;

    template class Mover< Selector<Atom> >;
}

Atom* Atom::clone() const
{
    return new Atom(*this);
}
