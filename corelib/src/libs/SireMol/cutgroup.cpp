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

#include "cutgroup.h"

#include "atom.h"
#include "molecule.h"
#include "cgeditor.h"

#include "mover.hpp"
#include "selector.hpp"
#include "evaluator.h"

#include "cgatomidx.h"
#include "groupatomids.h"

#include "SireBase/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>

#include "mover_metaid.h"

using namespace SireMol;
using namespace SireStream;

///////
/////// Implementation of CGProp
///////

static const RegisterMetaType<CGProp> r_cgprop(MAGIC_ONLY,
                                               "SireMol::CGProp");
                                                   
/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const CGProp &cgprop)
{
    writeHeader(ds, r_cgprop, 1)
         << static_cast<const MolViewProperty&>(cgprop);
         
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, CGProp &cgprop)
{
    VersionID v = readHeader(ds, r_cgprop);
    
    if (v == 1)
    {
        ds >> static_cast<MolViewProperty&>(cgprop);
    }
    else
        throw version_error(v, "1", r_cgprop, CODELOC);
        
    return ds;
}

CGProp::CGProp() : MolViewProperty()
{}

CGProp::CGProp(const CGProp &other) : MolViewProperty(other)
{}

CGProp::~CGProp()
{}

///////
/////// Implementation of CutGroup
///////

RegisterMetaType<CutGroup> r_cg;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const CutGroup &cg)
{
    writeHeader(ds, r_cg, 1);

    SharedDataStream sds(ds);
    
    sds << cg.cgidx << static_cast<const MoleculeView&>(cg);

    return ds;
}

/** Deserialise from a binary datastream */
QDataStream &operator>>(QDataStream &ds, CutGroup &cg)
{
    VersionID v = readHeader(ds, r_cg);

    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> cg.cgidx >> static_cast<MoleculeView&>(cg);
    }
    else
        throw version_error(v, "1", r_cg, CODELOC);

    return ds;
}

/** Null constructor */
CutGroup::CutGroup() : ConcreteProperty<CutGroup,MoleculeView>(), cgidx( CGIdx::null() )
{}

/** Construct the CutGroup at ID 'cgid' in the molecule whose data
    is in 'moldata'
    
    \throw SireMol::missing_CutGroup
    \throw SireMol::duplicate_CutGroup
    \throw SireError::invalid_index
*/
CutGroup::CutGroup(const MoleculeData &moldata, const CGID &cgid)
         : ConcreteProperty<CutGroup,MoleculeView>(moldata), 
           cgidx( moldata.info().cgIdx(cgid) )
{}

/** Copy constructor */
CutGroup::CutGroup(const CutGroup &other)
         : ConcreteProperty<CutGroup,MoleculeView>(other), cgidx(other.cgidx)
{}

/** Destructor */
CutGroup::~CutGroup()
{}

/** Copy assignment operator */
CutGroup& CutGroup::operator=(const CutGroup &other)
{
    MoleculeView::operator=(other);
    cgidx = other.cgidx;
    return *this;
}

/** Comparison operator */
bool CutGroup::operator==(const CutGroup &other) const
{
    return cgidx == other.cgidx and 
           MoleculeView::operator==(other);
}

/** Comparison operator */
bool CutGroup::operator!=(const CutGroup &other) const
{
    return cgidx != other.cgidx or
           MoleculeView::operator!=(other);
}

/** Return a string representation of this CutGroup */
QString CutGroup::toString() const
{
    return QObject::tr( "CutGroup( %1 )" ).arg(this->name());
}

/** Is this CutGroup empty? */
bool CutGroup::isEmpty() const
{
    return cgidx.isNull();
}

/** Is this CutGroup the whole molecule? */
bool CutGroup::selectedAll() const
{
    return MoleculeView::data().info().nCutGroups() == 1;
}

/** Return the atoms that are in this CutGroup */
AtomSelection CutGroup::selection() const
{
    AtomSelection selected_atoms(this->data());
    selected_atoms.selectOnly(cgidx);
    
    return selected_atoms;
}

/** Update this CutGroup with the passed molecule data.

    \throw SireError::incompatible_error
*/
void CutGroup::update(const MoleculeData &moldata)
{
    //check that the new data is compatible (has same molecule
    //number and info ID number)
    if (d->number() != moldata.number() or
        d->info().UID() != moldata.info().UID())
    {
        throw SireError::incompatible_error( QObject::tr(
            "You can only update a CutGroup with the molecule data "
            "for the same molecule (same molecule number) and that "
            "has a .info() object that has the same UID. You are "
            "trying to update CutGroup %1 in molecule %2 with UID %3 "
            "with molecule %4 with UID %5.")
                .arg(cgidx).arg(d->number()).arg(d->info().UID().toString())
                .arg(moldata.number()).arg(moldata.info().UID().toString()),
                    CODELOC );
    }
    
    d = moldata;
}

/** Return the name of this CutGroup */
const CGName& CutGroup::name() const
{
    return d->info().name(cgidx);
}

/** Return the index of this CutGroup in the molecule */
CGIdx CutGroup::index() const
{
    return cgidx;
}

/** Return an object that can move a copy of this CutGroup */
Mover<CutGroup> CutGroup::move() const
{
    return Mover<CutGroup>(*this);
}

/** Return an evaluator that can evaluate properties 
    of this CutGroup */
Evaluator CutGroup::evaluate() const
{
    return Evaluator(*this);
}

/** Return an editor that can edit this CutGroup */
CGEditor CutGroup::edit() const
{
    return CGEditor(*this);
}

/** Return a selector that can change the selection of CutGroups */
Selector<CutGroup> CutGroup::selector() const
{
    return Selector<CutGroup>(*this);
}

/** Return the number of atoms in this CutGroup */
int CutGroup::nAtoms() const
{
    return d->info().nAtoms(cgidx);
}

/** Return the indicies of the atoms in this CutGroup, in the
    order they appear in this CutGroup */
const QList<AtomIdx>& CutGroup::atomIdxs() const
{
    return d->info().getAtomsIn(cgidx);
}

/** Return whether or not this CutGroup contains the atom 
    at index 'atomidx' in the molecule */
bool CutGroup::contains(AtomIdx atomidx) const
{
    return d->info().contains(cgidx, atomidx);
}

/** Return whether or not this CutGroup contains all of 
    the atoms that match the ID 'atomid' */
bool CutGroup::contains(const AtomID &atomid) const
{
    return d->info().contains(cgidx, atomid);
}

/** Return whether or not this CutGroup contains some of
    the atoms that match the ID 'atomid' */
bool CutGroup::intersects(const AtomID &atomid) const
{
    return d->info().intersects(cgidx, atomid);
}

/** Return whether or not there is a CGProperty at key 'key' */
bool CutGroup::hasProperty(const PropertyName &key) const
{
    return d->hasPropertyOfType<CGProp>(key);
}

/** Return whether or not there is a CGProperty at metakey 'metakey' */
bool CutGroup::hasMetadata(const PropertyName &metakey) const
{
    return d->hasMetadataOfType<CGProp>(metakey);
}

/** Return the keys of all CGProperty properties */
QStringList CutGroup::propertyKeys() const
{
    return d->properties().propertyKeysOfType<CGProp>();
}

/** Return the metakeys of all CGProperty metadata */
QStringList CutGroup::metadataKeys() const
{
    return d->properties().metadataKeysOfType<CGProp>();
}

/** Return the metakeys of all CGProperty metadata for 
    the property at key 'key'
    
    \throw SireBase::missing_property
*/
QStringList CutGroup::metadataKeys(const PropertyName &key) const
{
    return d->properties().metadataKeysOfType<CGProp>(key);
}

/** Return whether the metadata at metakey 'metakey' for the property
    at key 'key' is a CGProperty
    
    \throw SireBase::missing_property
*/
bool CutGroup::hasMetadata(const PropertyName &key,
                       const PropertyName &metakey) const
{
    return d->hasMetadataOfType<CGProp>(key, metakey);
}

/** Assert that this CutGroup has an CGProperty at key 'key'

    \throw SireBase::missing_property
*/
void CutGroup::assertContainsProperty(const PropertyName &key) const
{
    if (not this->hasProperty(key))
        throw SireBase::missing_property( QObject::tr(
            "There is no CGProperty at key '%1' for this CutGroup.")
                .arg(key.toString()), CODELOC );
}

/** Assert that this CutGroup has an CGProperty piece of metadata
    at metakey 'metakey'
    
    \throw SireBase::missing_property
*/
void CutGroup::assertContainsMetadata(const PropertyName &metakey) const
{
    if (not this->hasMetadata(metakey))
        throw SireBase::missing_property( QObject::tr(
            "There is no CGProperty metadata at metakey '%1' for "
            "this CutGroup.")
                .arg(metakey.toString()), CODELOC );
}

/** Assert that the property at key 'key' has an CGProperty
    piece of metadata at metakey 'metakey'
    
    \throw SireBase::missing_property
*/
void CutGroup::assertContainsMetadata(const PropertyName &key,
                                  const PropertyName &metakey) const
{
    if (not this->hasMetadata(key, metakey))
        throw SireBase::missing_property( QObject::tr(
            "There is no CGProperty metadata at metakey '%1' "
            "for the property at key '%2' for this CutGroup.")
                .arg(metakey.toString(), key.toString()), CODELOC );
}

const char* CutGroup::typeName()
{
    return QMetaType::typeName( qMetaTypeId<CutGroup>() );
}

namespace SireMol
{
namespace detail
{

bool has_property(const CutGroup*, const MoleculeData &moldata,
                                 const PropertyName &key)
{
    return moldata.hasPropertyOfType<CGProp>(key);
}

bool has_metadata(const CutGroup*, const MoleculeData &moldata,
                                 const PropertyName &metakey)
{
    return moldata.hasMetadataOfType<CGProp>(metakey);
}

bool has_metadata(const CutGroup*, const MoleculeData &moldata,
                                 const PropertyName &key, const PropertyName &metakey)
{
    return moldata.hasMetadataOfType<CGProp>(key, metakey);
}

} // end of namespace detail
} // end of namespace SireMol

namespace SireMol
{
    ////// explicitly instantiate the CutGroup templates
    template class Mover<CutGroup>;
    template class Selector<CutGroup>;

    template class Mover< Selector<CutGroup> >;
}

CutGroup* CutGroup::clone() const
{
    return new CutGroup(*this);
}
