/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2010  Christopher Woods
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

#include "bead.h"
#include "beadeditor.h"
#include "beads.h"
#include "mover.hpp"
#include "selector.hpp"
#include "atom.h"
#include "residue.h"
#include "chain.h"
#include "cutgroup.h"
#include "segment.h"
#include "partialmolecule.h"

#include "SireBase/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "tostring.h"

using namespace SireMol;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<Bead> r_bead;

QDataStream &operator<<(QDataStream &ds, const Bead &bead)
{
    writeHeader(ds, r_bead, 1);

    SharedDataStream sds(ds);

    sds << bead.beadidx << bead.bdng
        << bead.beading_property << bead.selected_atoms
        << static_cast<const MoleculeView&>(bead);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, Bead &bead)
{
    VersionID v = readHeader(ds, r_bead);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> bead.beadidx >> bead.bdng
            >> bead.beading_property >> bead.selected_atoms
            >> static_cast<MoleculeView&>(bead);
    }
    else
        throw version_error(v, "1", r_bead, CODELOC);

    return ds;
}

/** Null constructor */
Bead::Bead() : ConcreteProperty<Bead,MoleculeView>(), beading_property("beading")
{}

/** Construct a bead view of the passed molecule, of the bead with index 'bead',
    using the passed property map to find the required beading property */
Bead::Bead(const MoleculeData &moldata, const BeadIdx &bead,
           const PropertyMap &map)
     : ConcreteProperty<Bead,MoleculeView>(moldata), beadidx(bead)
{
    beading_property = map["beading"];

    if (not moldata.hasProperty(beading_property))
    {
        bdng = ResidueBeading();
    }
    else
    {
        bdng = moldata.property(beading_property).asA<Beading>();
    }

    beadidx = BeadIdx( beadidx.map(bdng.read().nBeads(moldata.info())) );
    selected_atoms = bdng.read().selection(moldata.info(), beadidx);
}

/** Internal constructor */
Bead::Bead(const MoleculeData &moldata, BeadIdx idx,
           const Beading &beading, const PropertyName &prop)
     : ConcreteProperty<Bead,MoleculeView>(moldata),
       bdng(beading), beading_property(prop)
{
    beadidx = BeadIdx( idx.map(beading.nBeads(moldata.info())) );
    selected_atoms = beading.selection(moldata.info(), beadidx);
}

/** Copy constructor */
Bead::Bead(const Bead &other)
     : ConcreteProperty<Bead,MoleculeView>(other),
       beadidx(other.beadidx), bdng(other.bdng),
       beading_property(other.beading_property),
       selected_atoms(other.selected_atoms)
{}

/** Destructor */
Bead::~Bead()
{}

/** Copy assignment operator */
Bead& Bead::operator=(const Bead &other)
{
    if (this != &other)
    {
        beadidx = other.beadidx;
        bdng = other.bdng;
        beading_property = other.beading_property;
        selected_atoms = other.selected_atoms;
        MoleculeView::operator=(other);
    }

    return *this;
}

/** Comparison operator */
bool Bead::operator==(const Bead &other) const
{
    return beadidx == other.beadidx and MoleculeView::operator==(other) and
           beading_property == other.beading_property;
}

/** Comparison operator */
bool Bead::operator!=(const Bead &other) const
{
    return not Bead::operator==(other);
}

const char* Bead::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Bead>() );
}

Bead* Bead::clone() const
{
    return new Bead(*this);
}

MolViewPtr Bead::toSelector() const
{
    return PartialMolecule(*this).toSelector();
}

/** Return a string representation of this bead */
QString Bead::toString() const
{
    if (beadidx.isNull())
        return QObject::tr("Bead::null");
    else
        return QObject::tr("Bead( %1 : %2 : %3 )")
                    .arg(beadidx).arg(this->data().name())
                    .arg(bdng.read().toString());
}

/** Return if this is an empty bead */
bool Bead::isEmpty() const
{
    return selected_atoms.selectedNone();
}

/** Return whether or not this bead includes all of the atoms in the molecule */
bool Bead::selectedAll() const
{
    return selected_atoms.selectedAll();
}

/** Return the selection of atoms that are part of this bead */
AtomSelection Bead::selection() const
{
    return selected_atoms;
}

/** Update this bead using the passed molecule data */
void Bead::update(const MoleculeData &moldata)
{
    BeadingPtr new_beading;

    if (moldata.hasProperty(beading_property))
    {
        new_beading = moldata.property(beading_property).asA<Beading>();
    }
    else
    {
        new_beading = ResidueBeading();
    }

    new_beading.read().assertValidIndex(beadidx, moldata.info());

    MoleculeView::update(moldata);

    if (not new_beading.read().equals(bdng.read()))
    {
        bdng = new_beading;
        selected_atoms = bdng.read().selection(moldata.info());
    }
}

/** Return the index of this bead */
BeadIdx Bead::index() const
{
    return beadidx;
}

/** Return the set of all beads */
Beads Bead::beads() const
{
    return Beads(*this);
}

/** Return whether or not this bead has a property called 'key' */
bool Bead::hasProperty(const PropertyName &key) const
{
    if (this->data().hasProperty(key))
    {
        const Property &prop = this->data().property(key);

        if (prop.isA<BeadProp>())
        {
            if (prop.asA<BeadProp>().beading().equals(bdng.read()))
                return true;
        }
    }

    return false;
}

/** Return whether or not this bead had some metadata at 'metakey' */
bool Bead::hasMetadata(const PropertyName &metakey) const
{
    if (this->data().hasMetadata(metakey))
    {
        const Property &prop = this->data().metadata(metakey);

        if (prop.isA<BeadProp>())
        {
            if (prop.asA<BeadProp>().beading().equals(bdng.read()))
                return true;
        }
    }

    return false;
}

/** Return whether or not this bead has some metadata at 'key':'metakey' */
bool Bead::hasMetadata(const PropertyName &key,
                       const PropertyName &metakey) const
{
    if (this->data().hasMetadata(key,metakey))
    {
        const Property &prop = this->data().metadata(key,metakey);

        if (prop.isA<BeadProp>())
        {
            if (prop.asA<BeadProp>().beading().equals(bdng.read()))
                return true;
        }
    }

    return false;
}

/** Return a list of all of the properties associated with this bead */
QStringList Bead::propertyKeys() const
{
    QStringList beadprops = d->properties().propertyKeysOfType<BeadProp>();

    QMutableStringListIterator it(beadprops);

    while (it.hasNext())
    {
        const QString &key = it.next();

        if (not this->data().property(key).asA<BeadProp>().beading()
                            .equals(bdng.read()))
        {
            it.remove();
        }
    }

    return beadprops;
}

/** Return a list of all of the metadata properties associated with this bead */
QStringList Bead::metadataKeys() const
{
    QStringList beadprops = d->properties().metadataKeysOfType<BeadProp>();

    QMutableStringListIterator it(beadprops);

    while (it.hasNext())
    {
        const QString &metakey = it.next();

        if (not this->data().metadata(metakey).asA<BeadProp>().beading()
                            .equals(bdng.read()))
        {
            it.remove();
        }
    }

    return beadprops;
}

/** Return a list of all of the metadata properties of the key 'key' that
    are associated with this bead */
QStringList Bead::metadataKeys(const PropertyName &key) const
{
    QStringList beadprops = d->properties().metadataKeysOfType<BeadProp>(key);

    QMutableStringListIterator it(beadprops);

    while (it.hasNext())
    {
        const QString &metakey = it.next();

        if (not this->data().metadata(key, metakey).asA<BeadProp>().beading()
                            .equals(bdng.read()))
        {
            it.remove();
        }
    }

    return beadprops;
}

/** Return the mover for this bead */
Mover<Bead> Bead::move() const
{
    return Mover<Bead>(*this);
}

/** Return the evaluator for this bead */
Evaluator Bead::evaluate() const
{
    return Evaluator(*this);
}

/** Return the editor for this bead */
BeadEditor Bead::edit() const
{
    return BeadEditor(*this);
}

/** Return the number of atoms in this bead */
int Bead::nAtoms() const
{
    return selected_atoms.nSelected();
}

/** Return the ith atom in this bead

    \throw SireError::invalid_index
*/
MolViewPtr Bead::operator[](int i) const
{
    return Atom(data(), bdng.read().atomIdx(data().info(), beadidx, i));
}

int Bead::nViews() const
{
    return this->nAtoms();
}

/** Return the ith atom in this bead

    \throw SireError::invalid_index
*/
Atom Bead::atom(int i) const
{
    return Atom(data(), bdng.read().atomIdx(data().info(), beadidx, i));
}

/** Return the beading function used to bead up the molecule */
const Beading& Bead::beading() const
{
    return bdng.read();
}

/** Return the list of atom indexes of the atoms in this bead */
QList<AtomIdx> Bead::atomIdxs() const
{
    return bdng.read().atomIdxs(this->data().info());
}

/** Return whether or not this bead contains the atom with index 'atomidx' */
bool Bead::contains(AtomIdx atomidx) const
{
    return selected_atoms.selected(atomidx);
}

/** Return whether or not this bead contains the atom with ID 'atomid' */
bool Bead::contains(const AtomID &atomid) const
{
    return selected_atoms.selected(atomid);
}

/** Return whether or not this bead contains the atom with ID 'atomid' */
bool Bead::intersects(const AtomID &atomid) const
{
    return selected_atoms.selected(atomid);
}

/** Assert that this bead contains the property with key 'key'

    \throw SireBase::missing_property
*/
void Bead::assertContainsProperty(const PropertyName &key) const
{
    if (not this->hasProperty(key))
        throw SireBase::missing_property( QObject::tr(
                "The bead %1 does not contain the property %2. Available "
                "properties are %3.")
                    .arg(this->toString(), key.toString(),
                         Sire::toString(this->propertyKeys())), CODELOC );
}

/** Assert that this bead contains the metadata with key 'metakey'

    \throw SireBase::missing_property
*/
void Bead::assertContainsMetadata(const PropertyName &metakey) const
{
    if (not this->hasMetadata(metakey))
        throw SireBase::missing_property( QObject::tr(
                "The bead %1 does not contain the metadata %2. Available "
                "metadata values are %3.")
                    .arg(this->toString(), metakey.toString(),
                         Sire::toString(this->metadataKeys())), CODELOC );
}

/** Assert that this bead contains the metadata property with key 'key':'metakey'

    \throw SireBase::missing_property
*/
void Bead::assertContainsMetadata(const PropertyName &key,
                                  const PropertyName &metakey) const
{
    if (not this->hasMetadata(key,metakey))
        throw SireBase::missing_property( QObject::tr(
                "The bead %1 does not contain the metadata %2:%3. Available "
                "properties are %4.")
                    .arg(this->toString(), key.toString(), metakey.toString(),
                         Sire::toString(this->metadataKeys(key))), CODELOC );
}
