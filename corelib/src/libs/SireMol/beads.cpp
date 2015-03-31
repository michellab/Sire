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

#include "beads.h"
#include "mover.hpp"

#include "SireStream/datastream.h"
#include "SireStream/datastream.h"

using namespace SireMol;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<Beads> r_beads;

QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds, const Beads &beads)
{
    writeHeader(ds, r_beads, 1);
    
    SharedDataStream sds(ds);
    
    sds << beads.bdng << beads.beading_property << beads.selected_atoms
        << static_cast<const MoleculeView&>(beads);
        
    return ds;
}

QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds, Beads &beads)
{
    VersionID v = readHeader(ds, r_beads);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> beads.bdng >> beads.beading_property >> beads.selected_atoms
            >> static_cast<MoleculeView&>(beads);
    }
    else
        throw version_error(v, "1", r_beads, CODELOC);
        
    return ds;
}

/** Null constructor */
Beads::Beads() : ConcreteProperty<Beads,MoleculeView>()
{}

/** Construct the beads view of the passed molecule, using
    the passed (optional) property map to find the "beading" property */
Beads::Beads(const MoleculeData &moldata, const PropertyMap &map)
      : ConcreteProperty<Beads,MoleculeView>(moldata)
{
    beading_property = map["beading"];
    
    if (moldata.hasProperty(beading_property))
    {
        bdng = moldata.property(beading_property).asA<Beading>();
    }
    else
    {
        bdng = ResidueBeading();
    }
    
    selected_atoms = bdng.read().selection(moldata.info());
}
    
/** Copy constructor */
Beads::Beads(const Beads &other)
      : ConcreteProperty<Beads,MoleculeView>(other),
        bdng(other.bdng), beading_property(other.beading_property),
        selected_atoms(other.selected_atoms)
{}

/** Destructor */
Beads::~Beads()
{}

/** Copy assignment operator */
Beads& Beads::operator=(const Beads &other)
{
    if (this != &other)
    {
        MoleculeView::operator=(other);
        bdng = other.bdng;
        beading_property = other.beading_property;
        selected_atoms = other.selected_atoms;
    }
    
    return *this;
}

/** Comparison operator */
bool Beads::operator==(const Beads &other) const
{
    return MoleculeView::operator==(other) and 
           beading_property == other.beading_property;
}

/** Comparison operator */
bool Beads::operator!=(const Beads &other) const
{
    return not Beads::operator==(other);
}

const char* Beads::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Beads>() );
}

Beads* Beads::clone() const
{
    return new Beads(*this);
}

/** Return the bead at index 'beadidx' 

    \throw SireError::invalid_index
*/
Bead Beads::operator[](BeadIdx beadidx) const
{
    return Bead(this->data(), beadidx, bdng.read(), beading_property); 
}

/** Return the bead at index 'beadidx' 

    \throw SireError::invalid_index
*/
Bead Beads::at(BeadIdx beadidx) const
{
    return this->operator[](beadidx);
}

/** Return the bead at index 'beadidx' 

    \throw SireError::invalid_index
*/
Bead Beads::bead(BeadIdx beadidx) const
{
    return this->operator[](beadidx);
}

/** Return a string representation of these beads */
QString Beads::toString() const
{
    if (selected_atoms.isEmpty())
        return QObject::tr("Beads::null");
    else
        return QObject::tr("Beads( %1 : nBeads() == %2, beading() == %3 )")
                    .arg(this->data().name())
                    .arg(this->nBeads())
                    .arg(this->beading().toString());
}

/** Return whether or not this is empty */
bool Beads::isEmpty() const
{
    return selected_atoms.isEmpty();
}

/** Return whether or not these beads contain all atoms */
bool Beads::selectedAll() const
{
    return selected_atoms.selectedAll();
}

/** Return the selection of atoms that are part of the beads */
AtomSelection Beads::selection() const
{
    return selected_atoms;
}

/** Return the number of beads */
int Beads::nBeads() const
{
    if (isEmpty())
        return 0;
    else
        return bdng.read().nBeads( this->data().info() );
}

/** Return the number of beads */
int Beads::count() const
{
    return nBeads();
}

/** Return the number of beads */
int Beads::size() const
{
    return nBeads();
}

/** Return the number of atoms in the beads */
int Beads::nAtoms() const
{
    return selected_atoms.nSelected();
}

/** Update these beads with the new molecule data */
void Beads::update(const MoleculeData &moldata)
{
    BeadingPtr new_beading = bdng;

    if (moldata.hasProperty(beading_property))
    {
        new_beading = moldata.property(beading_property).asA<Beading>();
    }
    else if (data().hasProperty(beading_property))
    {
        new_beading = ResidueBeading();
    }
    
    MoleculeView::update(moldata);
    
    if (not new_beading.read().equals(bdng.read()))
    {
        bdng = new_beading;
        selected_atoms = bdng.read().selection(moldata.info());
    }
}

/** Return a mover that acts of all of the beads */
Mover<Beads> Beads::move() const
{
    return Mover<Beads>(*this);
}

/** Return an evaluator for all of the beads */
Evaluator Beads::evaluate() const
{
    return Evaluator(*this);
}

/** Return the indicies of all of the atoms in all of the beads */
QList<AtomIdx> Beads::atomIdxs() const
{
    return bdng.read().atomIdxs(this->data().info());
}

/** Return the beading function used to create the beads */
const Beading& Beads::beading() const
{
    return bdng.read();
}

/** Return whether any of the beads contains the atom with index 'atomidx' */
bool Beads::contains(AtomIdx atomidx) const
{
    return selected_atoms.selected(atomidx);
}

/** Return whether or not any of the beads contains the atom with ID 'atomid' */
bool Beads::contains(const AtomID &atomid) const
{
    return selected_atoms.selected(atomid);
}

/** Return whether or not any of the beads contains the atom with ID 'atomid' */
bool Beads::intersects(const AtomID &atomid) const
{
    return selected_atoms.selected(atomid);
}

/** Return the atom properties for all of the atoms in the beads, in 
    BeadIdx/Index order */
PropertyPtr Beads::atomProperty(const PropertyName &key) const
{
    return bdng.read().atomProperty(this->data(), key);
}

/** At the moment, the "Beads" object has no properties or metadata */
bool Beads::hasProperty(const PropertyName &key) const
{
    return false;
}

/** At the moment, the "Beads" object has no properties or metadata */
bool Beads::hasMetadata(const PropertyName &metakey) const
{
    return false;
}

/** At the moment, the "Beads" object has no properties or metadata */
bool Beads::hasMetadata(const PropertyName &key,
                        const PropertyName &metakey) const
{
    return false;
}
                 
/** At the moment, the "Beads" object has no properties or metadata */
QStringList Beads::propertyKeys() const
{
    return QStringList();
}

/** At the moment, the "Beads" object has no properties or metadata */
QStringList Beads::metadataKeys() const
{
    return QStringList();
}

/** At the moment, the "Beads" object has no properties or metadata */
QStringList Beads::metadataKeys(const PropertyName &key) const
{
    return QStringList();
}
