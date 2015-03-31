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

#include <QDataStream>
#include <QSharedData>

#include "moleculedata.h"
#include "moleculeinfodata.h"
#include "structureeditor.h"
#include "moleculeview.h"
#include "atomselection.h"
#include "moleditor.h"

#include "SireBase/incremint.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>

using namespace SireStream;
using namespace SireMol;
using namespace SireBase;

static const RegisterMetaType<MoleculeData> r_moldata(NO_ROOT);

/////////
///////// Implementation of MoleculeData::PropVersions
/////////

/** Increment the global version number of the molecule,
    and return it */
quint64 MoleculeData::PropVersions::increment()
{
    QMutexLocker lkr(&mutex);
    
    ++version;
    
    return version;
}

/** Increment the version number for the property with key 'key' */
quint64 MoleculeData::PropVersions::increment(const QString &key)
{
    QHash<QString,quint64>::iterator it = property_version.find(key);
    
    if (it == property_version.end())
    {
        property_version.insert(key, 1);
        return 1;
    }
    else
    {
        it.value() = it.value() + 1;
        return it.value();
    }
}

/** Increment the version number of the property with key 'key',
    also updating and returning the global version number of 
    the molecule */
quint64 MoleculeData::PropVersions::increment(const QString &key,
                                              quint64 &molversion)
{
    QMutexLocker lkr(&mutex);
    
    //now increment the global version
    ++version;
    
    //return the values
    molversion = version;
    
    return this->increment(key);
}

/** Increment all of the IDs for all of the properties contained
    in the passed MoleculeData */
void MoleculeData::PropVersions::incrementAll(MoleculeData &moldata)
{
    QMutexLocker lkr(&mutex);
    
    //increment the version for all of the keys in the molecule
    moldata.prop_vrsns.clear();
    
    foreach (QString key, moldata.props.propertyKeys())
    {
        moldata.prop_vrsns.insert(key, this->increment(key));
    }
    
    //now increment the global version
    ++version;
    
    moldata.vrsn = version;
}

/////////
///////// Implementation of MoleculeData
/////////

/** Serialise to a binary data stream */
QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds, const MoleculeData &moldata)
{
    writeHeader(ds, r_moldata, 1);
    
    SharedDataStream sds(ds);
    
    sds << moldata.molinfo << moldata.props 
        << moldata.molname << moldata.molnum;
    
    return ds;
}

/** Deserialise from a binary data stream */
QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds, MoleculeData &moldata)
{
    VersionID v = readHeader(ds, r_moldata);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> moldata.molinfo;

        sds >> moldata.props;

        sds >> moldata.molname;

        sds >> moldata.molnum;

        //get a versions object for this molecule
        moldata.vrsns = MoleculeData::registerMolecule(moldata.molnum);

        //now get the version numbers for all parts of this molecule
        moldata.vrsns->incrementAll(moldata);
    }
    else
        throw version_error(v, "1", r_moldata, CODELOC);

    return ds;
}

QHash< MolNum, boost::weak_ptr<MoleculeData::PropVersions> >
MoleculeData::version_registry;

QMutex MoleculeData::version_registry_mutex;

/** Return the version object that is used to get version numbers
    for the molecule with ID number 'molnum' */
boost::shared_ptr<MoleculeData::PropVersions>
MoleculeData::registerMolecule(MolNum molnum)
{
    QMutexLocker lkr(&version_registry_mutex);

    boost::shared_ptr<PropVersions> vrsns 
                    = version_registry[molnum].lock();

    if (not vrsns)
    {
        vrsns.reset( new PropVersions() );
        version_registry[molnum] = vrsns;
        
        if (version_registry.capacity() - version_registry.count() < 10)
        {
            //ok, its time to try and clean out - remove all expired
            //molecules from the registry - this prevents us filling
            //up memory, and also may prevent the next malloc of the
            //hash (as the test is based on the remaining capacity
            //of the hash)
            QMutableHashIterator< MolNum, boost::weak_ptr<PropVersions > >
                                                it( version_registry );
            
            while( it.hasNext() )
            {
                it.next();
            
                if (it.value().expired())
                    it.remove();
            }
        }
    }
    
    return vrsns;
}

static Incremint last_registered_molnum(0);

MolNum MoleculeData::createUniqueMolNum()
{
    QMutexLocker lkr(&version_registry_mutex);
    
    MolNum molnum( last_registered_molnum.increment() );
    
    while (version_registry.contains(molnum))
    {
        molnum = MolNum( last_registered_molnum.increment() );
    }
    
    return molnum;
}

MolNum MolNum::getUniqueNumber()
{
    return MoleculeData::createUniqueMolNum();
}

/** Null constructor */
MoleculeData::MoleculeData()
             : QSharedData(),
               molinfo( MoleculeInfoData::null() ),
               vrsn(0), 
               molnum(0),
               vrsns( MoleculeData::registerMolecule(molnum) )
{}

/** Get the underlying data viewed by the passed molecule view */
MoleculeData::MoleculeData(const MoleculeView &molview)
             : QSharedData()
{
    this->operator=( molview.data() );
}

static SharedDataPointer<MoleculeData> shared_null( new MoleculeData() );

SharedDataPointer<MoleculeData> MoleculeData::null()
{
    return shared_null;
}

/** Give this molecule a brand new unique ID number! */
void MoleculeData::renumber()
{
    //get the new ID number...
    molnum = MolNum::getUniqueNumber();
    
    vrsns = MoleculeData::registerMolecule(molnum);
    vrsns->incrementAll(*this);
}

/** Renumber this molecule to have the number 'molnum' */
void MoleculeData::renumber(MolNum newnum)
{
    if (newnum == molnum)
        //nothing to do
        return;
        
    molnum = newnum;
    vrsns = MoleculeData::registerMolecule(molnum);
    vrsns->incrementAll(*this);
}

/** Construct from the passed StructureEditor */
MoleculeData::MoleculeData(const StructureEditor &editor)
             : QSharedData(), vrsn(0), molnum(0)
{
    //create the info object from this editor
    molinfo = SharedDataPointer<MoleculeInfoData>( new MoleculeInfoData(editor) );
    
    //now copy across the properties...
    props = editor.properties();
    
    //copy across the name
    molname = editor.molName();
    
    //finally, sort out the molecule number - this also
    //sets up all of the version numbers and performs
    //the registration of the molecule
    if (editor.molNum().isNull())
        this->renumber();
    else
        this->renumber(editor.molNum());
}

/** Copy constructor */
MoleculeData::MoleculeData(const MoleculeData &other)
                : QSharedData(),
                  molinfo(other.molinfo),
                  props(other.props),
                  vrsn(other.vrsn),
                  molname(other.molname),
                  molnum(other.molnum),
                  prop_vrsns(other.prop_vrsns),
                  vrsns(other.vrsns)
{}

/** Destructor */
MoleculeData::~MoleculeData()
{}

/** Assignment operator */
MoleculeData& MoleculeData::operator=(const MoleculeData &other)
{
    if (this != &other)
    {
        molinfo = other.molinfo;
        props = other.props;
        vrsn = other.vrsn;
        molname = other.molname;
        molnum = other.molnum;
        prop_vrsns = other.prop_vrsns;
        vrsns = other.vrsns;
    }
    
    return *this;
}

/** Assign equal to a copy of the molecule being edited in 'editor' */
MoleculeData& MoleculeData::operator=(const StructureEditor &editor)
{
    return this->operator=( MoleculeData(editor) );
}

/** Comparison operator - two molecules are the same if they have
    the same ID and version numbers. */
bool MoleculeData::operator==(const MoleculeData &other) const
{
    return molnum == other.molnum and 
           vrsn == other.vrsn;
}

/** Comparison operator - two molecules are the same if they have
    the same ID and version numbers. */
bool MoleculeData::operator!=(const MoleculeData &other) const
{
    return molnum != other.molnum or vrsn != other.vrsn;
}

/** Return a new MoleculeData that contains only the passed selected
    atoms. This allows parts of the molecule to be pulled out and used independently */
MoleculeData MoleculeData::extract(const AtomSelection &selected_atoms) const
{
    selected_atoms.assertCompatibleWith(*this);
    
    if (selected_atoms.selectedAll())
        return *this;
    
    else if (selected_atoms.selectedNone())
        return MoleculeData();
    
    // edit a copy of this molecule
    MolStructureEditor editor = MolStructureEditor( Molecule(*this) );
    
    // delete all of the atoms that don't exist
    for (int i=molinfo->nAtoms()-1; i>=0; --i)
    {
        AtomIdx idx(i);
        
        if (not selected_atoms.selected(idx))
        {
            editor = editor.remove(idx);
        }
    }
    
    // delete all empty cutgroups
    for (int i=molinfo->nCutGroups()-1; i>=0; --i)
    {
        CGIdx idx(i);
        
        if (not selected_atoms.selected(idx))
        {
            editor = editor.remove(idx);
        }
    }
    
    //delete all empty segments
    for (int i=molinfo->nSegments()-1; i>=0; --i)
    {
        SegIdx idx(i);
        
        if (not selected_atoms.selected(idx))
        {
            editor = editor.remove(idx);
        }
    }
    
    //delete all empty chains
    for (int i=molinfo->nChains()-1; i>=0; --i)
    {
        ChainIdx idx(i);
        
        if (not selected_atoms.selected(idx))
        {
            editor = editor.remove(idx);
        }
    }
    
    //delete all empty residues
    for (int i=molinfo->nResidues()-1; i>=0; --i)
    {
        ResIdx idx(i);
        
        if (not selected_atoms.selected(idx))
        {
            editor = editor.remove(idx);
        }
    }
    
    return editor.commit().data();
}

/** Return the version number of the property at key 'key'.
    If there is no such key in this molecule, or 
    the value is supplied by the key itself, then
    a version number of 0 is returned */
quint64 MoleculeData::version(const PropertyName &key) const
{
    if (key.hasValue())
        return 0;
    else
        return prop_vrsns.value(key.source(), 0);
}

/** Return whether this molecule contains a property at key 'key' */
bool MoleculeData::hasProperty(const PropertyName &key) const
{
    return props.hasProperty(key);
}

/** Return whether this molecule contains metadata at 
    metakey 'metakey' */
bool MoleculeData::hasMetadata(const PropertyName &metakey) const
{
    return props.hasMetadata(metakey);
}

/** Return whether this molecule has metadata at metakey 'metakey'
    for the property at key 'key' 
    
    \throw SireBase::missing_property
*/
bool MoleculeData::hasMetadata(const PropertyName &key,
                               const PropertyName &metakey) const
{
    return props.hasMetadata(key, metakey);
}
                 
/** Return the type name of the property at key 'key'. 

    \throw SireBase::missing_property
*/
const char* MoleculeData::propertyType(const PropertyName &key) const
{
    return props.propertyType(key);
}

/** Return all of the keys of the properties of this molecule */
QStringList MoleculeData::propertyKeys() const
{
    return props.propertyKeys();
}

/** Return all of the metadata keys of the metadata of this molecule */
QStringList MoleculeData::metadataKeys() const
{
    return props.metadataKeys();
}

/** Return all of the metadata keys for the property at key 'key'

    \throw SireBase::missing_property
*/
QStringList MoleculeData::metadataKeys(const PropertyName &key) const
{
    return props.metadataKeys(key);
}

/** Return the type name of the metadata at metakey 'metakey'

    \throw SireBase::missing_property
*/
const char* MoleculeData::metadataType(const PropertyName &metakey) const
{
    return props.metadataType(metakey);
}

/** Return the type name of the metadata at metakey 'metakey'
    for the property at key 'key'
    
    \throw SireBase::missing_property
*/
const char* MoleculeData::metadataType(const PropertyName &key,
                                       const PropertyName &metakey) const
{
    return props.metadataType(key, metakey);
}

/** Return the property at key 'key' 

    \throw SireBase::missing_property
*/
const Property& MoleculeData::property(const PropertyName &key) const
{
    return props.property(key);
}

/** Return the property at key 'key', or 'default_value' if there
    is no such property */ 
const Property& MoleculeData::property(const PropertyName &key,
                                       const Property &default_value) const
{
    return props.property(key, default_value);
}

/** Return the metadata at metakey 'metakey'

    \throw SireBase::missing_property
*/
const Property& MoleculeData::metadata(const PropertyName &metakey) const
{
    return props.metadata(metakey);
}

/** Return the metadata at metakey 'metakey' for the property
    at key 'key'
    
    \throw SireBase::missing_property
*/
const Property& MoleculeData::metadata(const PropertyName &key,
                                       const PropertyName &metakey) const
{
    return props.metadata(key, metakey);
}

/** Return the metadata at metakey 'metakey', or 'default_value'
    if there is no such value */
const Property& MoleculeData::metadata(const PropertyName &metakey,
                                       const Property &default_value) const
{
    return props.metadata(metakey, default_value);
}
                         
/** Return the metadata at metakey 'metakey' for the property
    at key 'key', or 'default_value' if there is no such value
    
    \throw SireBase::missing_property
*/
const Property& MoleculeData::metadata(const PropertyName &key,
                                       const PropertyName &metakey,
                                       const Property &default_value) const
{
    return props.metadata(key, metakey, default_value);
}

/** Rename this molecule to 'newname'. This changes the info().UID()
    number, and the version number, but doesn't change this->number() */
void MoleculeData::rename(const MolName &newname)
{
    if (newname != molname)
    {
        molname = newname;
        vrsn = vrsns->increment();
    }
}

/** Rename the atom at index 'atomidx' to 'newname'

    \throw SireError::invalid_index
*/
void MoleculeData::rename(AtomIdx atomidx, const AtomName &newname)
{
    MoleculeInfoData newinfo = molinfo->rename(atomidx, newname);
    
    if (newinfo.UID() != molinfo.constData()->UID())
    {
        molinfo = newinfo;
        vrsn = vrsns->increment();
    }
}

/** Rename all atoms identified by 'atomid' to 'newname'

    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
void MoleculeData::rename(const AtomID &atomid, const AtomName &newname)
{
    MoleculeInfoData newinfo( *(molinfo.constData()) );
    
    foreach (AtomIdx atomidx, atomid.map(info()))
    {
        newinfo = newinfo.rename(atomidx, newname);
    }
    
    if (newinfo.UID() != molinfo.constData()->UID())
    {
        molinfo = newinfo;
        vrsn = vrsns->increment();
    }
}

/** Rename the CutGroup at index 'cgidx' to 'newname'

    \throw SireError::invalid_index
*/
void MoleculeData::rename(CGIdx cgidx, const CGName &newname)
{
    MoleculeInfoData newinfo = molinfo->rename(cgidx, newname);
    
    if (newinfo.UID() != molinfo.constData()->UID())
    {
        molinfo = newinfo;
        vrsn = vrsns->increment();
    }
}

/** Rename all CutGroups identified by 'cgid' to 'newname'

    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
*/
void MoleculeData::rename(const CGID &cgid, const CGName &newname)
{
    MoleculeInfoData newinfo( *(molinfo.constData()) );
    
    foreach (CGIdx cgidx, cgid.map(info()))
    {
        newinfo = newinfo.rename(cgidx, newname);
    }
    
    if (newinfo.UID() != molinfo.constData()->UID())
    {
        molinfo = newinfo;
        vrsn = vrsns->increment();
    }
}

/** Rename the residue at index 'residx' to 'newname'

    \throw SireError::invalid_index
*/
void MoleculeData::rename(ResIdx residx, const ResName &newname)
{
    MoleculeInfoData newinfo = molinfo->rename(residx, newname);
    
    if (newinfo.UID() != molinfo.constData()->UID())
    {
        molinfo = newinfo;
        vrsn = vrsns->increment();
    }
}

/** Rename all residues identified by 'resid' to 'newname'

    \throw SireMol::missing_residue
    \throw SireError::invalid_index
*/
void MoleculeData::rename(const ResID &resid, const ResName &newname)
{
    MoleculeInfoData newinfo( *(molinfo.constData()) );
    
    foreach (ResIdx residx, resid.map(info()))
    {
        newinfo = newinfo.rename(residx, newname);
    }
    
    if (newinfo.UID() != molinfo.constData()->UID())
    {
        molinfo = newinfo;
        vrsn = vrsns->increment();
    }
}

/** Rename the chain at index 'chainidx' to 'newname'

    \throw SireError::invalid_index
*/
void MoleculeData::rename(ChainIdx chainidx, const ChainName &newname)
{
    MoleculeInfoData newinfo = molinfo->rename(chainidx, newname);
    
    if (newinfo.UID() != molinfo.constData()->UID())
    {
        molinfo = newinfo;
        vrsn = vrsns->increment();
    }
}

/** Rename all chains identified by 'chainid' to 'newname'

    \throw SireMol::missing_chain
    \throw SireError::invalid_index
*/
void MoleculeData::rename(const ChainID &chainid, const ChainName &newname)
{
    MoleculeInfoData newinfo( *(molinfo.constData()) );
    
    foreach (ChainIdx chainidx, chainid.map(info()))
    {
        newinfo = newinfo.rename(chainidx, newname);
    }
    
    if (newinfo.UID() != molinfo.constData()->UID())
    {
        molinfo = newinfo;
        vrsn = vrsns->increment();
    }
}

/** Rename the segment at index 'segidx' to 'newname'

    \throw SireError::invalid_index
*/
void MoleculeData::rename(SegIdx segidx, const SegName &newname)
{
    MoleculeInfoData newinfo = molinfo->rename(segidx, newname);
    
    if (newinfo.UID() != molinfo.constData()->UID())
    {
        molinfo = newinfo;
        vrsn = vrsns->increment();
    }
}

/** Rename all segments identified by 'atomid' to 'newname'

    \throw SireMol::missing_segment
    \throw SireError::invalid_index
*/
void MoleculeData::rename(const SegID &segid, const SegName &newname)
{
    MoleculeInfoData newinfo( *(molinfo.constData()) );
    
    foreach (SegIdx segidx, segid.map(info()))
    {
        newinfo = newinfo.rename(segidx, newname);
    }
    
    if (newinfo.UID() != molinfo.constData()->UID())
    {
        molinfo = newinfo;
        vrsn = vrsns->increment();
    }
}

/** Renumber the atom at index 'atomidx' to 'newnum'

    \throw SireError::invalid_index
*/
void MoleculeData::renumber(AtomIdx atomidx, AtomNum newnum)
{
    MoleculeInfoData newinfo = molinfo->renumber(atomidx, newnum);
    
    if (newinfo.UID() != molinfo.constData()->UID())
    {
        molinfo = newinfo;
        vrsn = vrsns->increment();
    }
}

/** Renumber all of the atoms that match the ID 'atomid' to 'newnum'

    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
void MoleculeData::renumber(const AtomID &atomid, AtomNum newnum)
{
    MoleculeInfoData newinfo( *(molinfo.constData()) );
    
    foreach (AtomIdx atomidx, atomid.map(info()))
    {
        newinfo = newinfo.renumber(atomidx, newnum);
    }
    
    if (newinfo.UID() != molinfo.constData()->UID())
    {
        molinfo = newinfo;
        vrsn = vrsns->increment();
    }
}

/** Renumber the residue at index 'residx' to 'newnum'

    \throw SireError::invalid_index
*/
void MoleculeData::renumber(ResIdx residx, ResNum newnum)
{
    MoleculeInfoData newinfo = molinfo->renumber(residx, newnum);
    
    if (newinfo.UID() != molinfo.constData()->UID())
    {
        molinfo = newinfo;
        vrsn = vrsns->increment();
    }
}

/** Renumber all of the residues that match the ID 'resid' to 'newnum'

    \throw SireMol::missing_residue
    \throw SireError::invalid_index
*/
void MoleculeData::renumber(const ResID &resid, ResNum newnum)
{
    MoleculeInfoData newinfo( *(molinfo.constData()) );
    
    foreach (ResIdx residx, resid.map(info()))
    {
        newinfo = newinfo.renumber(residx, newnum);
    }
    
    if (newinfo.UID() != molinfo.constData()->UID())
    {
        molinfo = newinfo;
        vrsn = vrsns->increment();
    }
}

/** Set the property with the key 'key' to the value 'value'.
    This replaces any current property with that key. This
    also checks to ensure that this property is compatible
    with this molecule (e.g. it ensures that if the property
    is an atomic property, then it has the right number of
    values for the atoms)
    
    \throw SireError::incompatible_error
*/
void MoleculeData::setProperty(const QString &key, 
                               const Property &value,
                               bool clear_metadata)
{
    if (key.isEmpty())
        throw SireError::invalid_arg( QObject::tr(
            "You cannot set a property with an empty key!"), CODELOC );
        
    //now the property version number
    prop_vrsns.insert( key, vrsns->increment(key, vrsn) );
    
    //now save the property itself
    props.setProperty(key, value, clear_metadata);
}

/** Completely remove the property at key 'key', together
    with all of its metadata. This does nothing if there is
    no property with this key */
void MoleculeData::removeProperty(const QString &key)
{
    if (props.hasProperty(key))
    {
        props.removeProperty(key);
        prop_vrsns.remove(key);
    
        //do not remove from the shared version numbers, in
        //case the user re-adds a property with this key - 
        //we have to still ensure that the version number is
        //unique :-)
        
        //increment the global version number
        vrsn = vrsns->increment();
    }
}

/** Remove the property at key 'key', returning the value
    of that property
    
    \throw SireBase::missing_property
*/
PropertyPtr MoleculeData::takeProperty(const QString &key)
{
    PropertyPtr value = this->property(key);
    this->removeProperty(key);
    
    return value;
}

/** Remove the metadata at metakey 'metakey', returning the value
    of the metadata
    
    \throw SireBase::missing_property
*/
PropertyPtr MoleculeData::takeMetadata(const QString &metakey)
{
    PropertyPtr value = this->metadata(metakey);
    this->removeMetadata(metakey);
    
    return value;
}

/** Remove the metadata at metakey 'metakey' from the property
    at key 'key', returning the value of the metadata
    
    \throw SireBase::missing_property
*/
PropertyPtr MoleculeData::takeMetadata(const QString &key, const QString &metakey)
{
    PropertyPtr value = this->metadata(key, metakey);
    this->removeMetadata(key,metakey);
    
    return value;
}

/** Set the value of the metadata at metakey 'metakey' to 
    the value 'value' */
void MoleculeData::setMetadata(const QString &metakey, const Property &value)
{
    if (metakey.isEmpty())
        throw SireError::invalid_arg( QObject::tr(
            "You cannot set some metadata with an empty metakey!"), CODELOC );

    props.setMetadata(metakey, value);

    //increment the global version number
    vrsn = vrsns->increment();
}

/** Set the value of the metadata at metakey 'metakey' of the 
    property at key 'key' to the value 'value'
    
    \throw SireBase::missing_property
*/ 
void MoleculeData::setMetadata(const QString &key, const QString &metakey, 
                               const Property &value)
{
    if (key.isNull())
        throw SireError::invalid_arg( QObject::tr(
            "You cannot set some metadata with an empty key!"), CODELOC );
    else if (metakey.isEmpty())
        throw SireError::invalid_arg( QObject::tr(
            "You cannot set some metadata with an empty metakey!"), CODELOC );

    props.setMetadata(key, metakey, value);
    
    //increment the global version number
    vrsn = vrsns->increment();
}

/** Remove the metadata at metakey 'metakey' */
void MoleculeData::removeMetadata(const QString &metakey)
{
    if (props.hasMetadata(metakey))
    {
        props.removeMetadata(metakey);
        vrsn = vrsns->increment();
    }
}

/** Remove the metadata at metakey 'metakey' from the 
    property at key 'key'
    
    \throw SireBase::missing_property
*/
void MoleculeData::removeMetadata(const QString &key, const QString &metakey)
{
    if (props.hasMetadata(key, metakey))
    {
        props.removeMetadata(key, metakey);
        vrsn = vrsns->increment();
    }
}

const char* MoleculeData::typeName()
{
    return QMetaType::typeName( qMetaTypeId<MoleculeData>() );
}
