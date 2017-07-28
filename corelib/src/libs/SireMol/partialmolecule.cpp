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

#include "partialmolecule.h"
#include "molecule.h"
#include "atom.h"
#include "cutgroup.h"
#include "residue.h"
#include "chain.h"
#include "segment.h"

#include "evaluator.h"

#include "mover.hpp"
#include "selector.hpp"
#include "editor.hpp"

#include "SireError/errors.h"
#include "SireMol/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "tostring.h"

#include <QDebug>

using namespace SireMol;
using namespace SireBase;
using namespace SireStream;

RegisterMetaType<PartialMolecule> r_partialmol;

/** Serialise to a binary datastream */
QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds, 
                                       const PartialMolecule &partialmol)
{
    writeHeader(ds, r_partialmol, 1);

    SharedDataStream sds(ds);
    
    sds << partialmol.selected_atoms
        << static_cast<const MoleculeView&>(partialmol);

    return ds;
}

/** Deserialise from a binary datastream */
QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds, 
                                       PartialMolecule &partialmol)
{
    VersionID v = readHeader(ds, r_partialmol);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> partialmol.selected_atoms
            >> static_cast<MoleculeView&>(partialmol);
    }
    else
        throw version_error(v, "1", r_partialmol, CODELOC);

    return ds;
}

/** Null constructor */
PartialMolecule::PartialMolecule() : ConcreteProperty<PartialMolecule,MoleculeView>()
{}

/** Construct from the passed view */
PartialMolecule::PartialMolecule(const MoleculeView &molview)
                : ConcreteProperty<PartialMolecule,MoleculeView>(molview),
                  selected_atoms(molview.selection())
{}

/** Construct from the selected atoms of the passed molecule whose
    data is in 'moldata' */
PartialMolecule::PartialMolecule(const MoleculeData &moldata,
                                 const AtomSelection &atoms)
                : ConcreteProperty<PartialMolecule,MoleculeView>(moldata),
                  selected_atoms(atoms)
{
    atoms.assertCompatibleWith(moldata);
}

/** Copy constructor */
PartialMolecule::PartialMolecule(const PartialMolecule &other)
                : ConcreteProperty<PartialMolecule,MoleculeView>(other),
                  selected_atoms(other.selected_atoms)
{}           

/** Destructor */
PartialMolecule::~PartialMolecule()
{}

/** Copy assignment operator */
PartialMolecule& PartialMolecule::operator=(const MoleculeView &other)
{
    MoleculeView::operator=(other);
    selected_atoms = other.selection();
    return *this;
}

/** Copy assignment operator */
PartialMolecule& PartialMolecule::operator=(const PartialMolecule &other)
{
    MoleculeView::operator=(other);
    selected_atoms = other.selected_atoms;
    return *this;
}

/** Comparison operator */
bool PartialMolecule::operator==(const PartialMolecule &other) const
{
    return MoleculeView::operator==(other) and
           selected_atoms == other.selected_atoms;
}

/** Comparison operator */
bool PartialMolecule::operator!=(const PartialMolecule &other) const
{
    return MoleculeView::operator!=(other) or
           selected_atoms != other.selected_atoms;
}

/** Return a string representation of this molecule */
QString PartialMolecule::toString() const
{
    return QObject::tr( "PartialMolecule( %1 : %2 : nAtoms() == %3 )" )
                .arg(this->name())
                .arg(this->number())
                .arg(selected_atoms.nSelected());
}

/** Return whether or not this is empty */
bool PartialMolecule::isEmpty() const
{
    return selected_atoms.selectedNone();
}

/** Return whether or not this contains the entire molecule */
bool PartialMolecule::selectedAll() const
{
    return selected_atoms.selectedAll();
}

/** Return the name of this molecule */
const MolName& PartialMolecule::name() const
{
    return d->name();
}

/** Return the identifying number of this molecule */
MolNum PartialMolecule::number() const
{
    return d->number();
}
 
/** Return the version number of this molecule - all molecules
    with the same ID number and version number must be identical */
quint64 PartialMolecule::version() const
{
    return d->version();
}
 
/** Return the version number of the property at key 'key'.
    All molecules with the same ID number and same property version
    number must have the same value of this property
    (although this says nothing about any metadata associated
    with this property)
    
    \throw SireBase::missing_property 
*/
quint64 PartialMolecule::version(const PropertyName &key) const
{
    return d->version(key);
}

/** Return the number of atoms in this view */
int PartialMolecule::nAtoms() const
{
    return selected_atoms.nSelected();
}

/** Return the number of CutGroups in this view */
int PartialMolecule::nCutGroups() const
{
    return selected_atoms.nSelectedCutGroups();
}

/** Return the number of residues in this view */
int PartialMolecule::nResidues() const
{
    return selected_atoms.nSelectedResidues();
}

/** Return the number of chains in this view */
int PartialMolecule::nChains() const
{
    return selected_atoms.nSelectedChains();
}

/** Return the number of segments in this view */
int PartialMolecule::nSegments() const
{
    return selected_atoms.nSelectedSegments();
}

/** Return a mover that can move all of the atoms in this view */
Mover<PartialMolecule> PartialMolecule::move() const
{
    return Mover<PartialMolecule>(*this);
}

/** Return an evaluator that can evaluate properties 
    over all of the atoms in this view */
Evaluator PartialMolecule::evaluate() const
{
    return Evaluator(*this);
}

/** Return the atoms that are part of this view */
AtomSelection PartialMolecule::selection() const
{
    return selected_atoms;
}

/** Return whether or not this molecule has a property at key 'key' */
bool PartialMolecule::hasProperty(const PropertyName &key) const
{
    return d->hasProperty(key);
}

/** Return whether or not this molecule has some metadata
    at the metakey 'metakey' */
bool PartialMolecule::hasMetadata(const PropertyName &metakey) const
{
    return d->hasMetadata(metakey);
}

/** Return whether or not the property at key 'key' has
    some metadata at metakey 'metakey'
    
    \throw SireBase::missing_property
*/
bool PartialMolecule::hasMetadata(const PropertyName &key, 
                                  const PropertyName &metakey) const
{
    return d->hasMetadata(key, metakey);
}

/** Return the keys of all of the properties contained in this molecule */
QStringList PartialMolecule::propertyKeys() const
{
    return d->propertyKeys();
}

/** Extract a copy of this PartialMolecule which contains only the currently
    selected atoms. This allows the used to pull out parts of a larger molecule,
    e.g. if they want to have only selected residues in a protein and do not
    want to have to store or manipulate the larger protein molecule */
PartialMolecule PartialMolecule::extract() const
{
    if (this->isEmpty())
        return PartialMolecule();

    else if (this->selectedAll())
        return *this;
    
    else
    {
        PartialMolecule ret;
    
        ret.d = d->extract(selected_atoms);
        ret.selected_atoms = AtomSelection( *(ret.d) );
        
        return ret;
    }
}

/** Return the keys of all of the metadata contained directly by 
    this molecule */
QStringList PartialMolecule::metadataKeys() const
{
    return d->metadataKeys();
}

/** Return the keys of all metadata for the property at key 'key' 

    \throw SireBase::missing_property
*/
QStringList PartialMolecule::metadataKeys(const PropertyName &key) const
{
    return d->metadataKeys(key);
}

/** Return the property at key 'key'. Note that this returns the 
    property for the molecule - no attempt is made to mask this
    property to match the current selection
    
    \throw SireBase::missing_property
*/
const Property& PartialMolecule::property(const PropertyName &key) const
{
    return d->property(key);
}

/** Return the metadata at metakey 'metakey'. Note that this returns the 
    metadata for the molecule - no attempt is made to mask this
    metadata to match the current selection
    
    \throw SireBase::missing_property
*/
const Property& PartialMolecule::metadata(const PropertyName &metakey) const
{
    return d->metadata(metakey);
}

/** Return the metadata at the metakey 'metakey' for the property
    at key 'key'. Note that this returns the 
    metadata for the molecule - no attempt is made to mask this
    metadata to match the current selection
    
    \throw SireBase::missing_property
*/
const Property& PartialMolecule::metadata(const PropertyName &key,
                                          const PropertyName &metakey) const
{
    return d->metadata(key, metakey);
}

const char* PartialMolecule::typeName()
{
    return QMetaType::typeName( qMetaTypeId<PartialMolecule>() );
}

/** Return a copy of this PartialMolecule that has been reduced to its unit
    type, i.e. if this is a single Atom, this returns the Atom, if this is a single
    residue, this returns the Residue etc. */
MolViewPtr PartialMolecule::toUnit() const
{
    if (this->selectedAll())
    {
        return this->molecule();
    }
    else if (this->nAtoms() == 1)
    {
        return this->atom();
    }

    if (this->nResidues() == 1)
    {
        const auto res = this->residue();

        if (this->selection().selectedAll(res.number()))
        {
            return res;
        }
    }

    if (this->nCutGroups() == 1)
    {
        const auto cg = this->cutGroup();

        if (this->selection().selectedAll(cg.name()))
        {
            return cg;
        }
    }

    if (this->nChains() == 1)
    {
        const auto chain = this->chain();
        
        if (this->selection().selectedAll(chain.name()))
        {
            return chain;
        }
    }

    if (this->nSegments() == 1)
    {
        const auto segment = this->segment();
        
        if (this->selection().selectedAll(segment.name()))
        {
            return segment;
        }
    }

    return *this;
}

namespace SireMol
{
    /////////
    ///////// Explicit instantiation of the template classes
    /////////

    template class Mover<PartialMolecule>;
}

PartialMolecule* PartialMolecule::clone() const
{
    return new PartialMolecule(*this);
}
