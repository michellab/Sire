/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2006  Christopher Woods
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

#ifndef SIREMOL_MOLECULEDATA_H
#define SIREMOL_MOLECULEDATA_H

#include <QMutex>

#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>

#include "SireBase/properties.h"
#include "SireBase/propertymap.h"
#include "SireBase/refcountdata.h"

#include "SireBase/shareddatapointer.hpp"

#include "molname.h"
#include "molnum.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class MoleculeData;
}

QDataStream& operator<<(QDataStream&, const SireMol::MoleculeData&);
QDataStream& operator>>(QDataStream&, SireMol::MoleculeData&);

namespace SireMol
{

class MoleculeInfoData;

class AtomMatcher;

class AtomNum;
class AtomIdx;
class AtomName;
class AtomID;

class CGIdx;
class CGName;
class CGID;

class MoleculeView;

class ResNum;
class ResIdx;
class ResName;
class ResID;

class ChainName;
class ChainIdx;
class ChainID;

class SegName;
class SegIdx;
class SegID;

class StructureEditor;

using SireBase::PropertyName;
using SireBase::Property;
using SireBase::PropertyPtr;
using SireBase::Properties;

/**
This class holds the shared molecule data for the Molecule classes
(which are all just views on this MolculeData class).

This class is not part of the standard API of the program and should not
be used in your own code.

@author Christopher Woods
*/
class SIREMOL_EXPORT MoleculeData : public SireBase::RefCountData
{

friend QDataStream& ::operator<<(QDataStream&, const MoleculeData&);
friend QDataStream& ::operator>>(QDataStream&, MoleculeData&);

public:
    MoleculeData();
    MoleculeData(const StructureEditor &editor);

    MoleculeData(const MoleculeView &molview);

    MoleculeData(const MoleculeData &other);

    ~MoleculeData();

    static const char* typeName();
    
    const char* what() const
    {
        return MoleculeData::typeName();
    }

    MoleculeData& operator=(const MoleculeData &other);
    MoleculeData& operator=(const StructureEditor &editor);

    bool operator==(const MoleculeData &other) const;
    bool operator!=(const MoleculeData &other) const;

    /** Return the name of the molecule */
    const MolName& name() const
    {
        return molname;
    }

    /** The ID number of this molecule - two molecules with
        the same ID number are the same (though potentially
        at different versions) */
    MolNum number() const
    {
        return molnum;
    }

    /** The version number of this molecule - two molecules
        with the same version number and ID number are identical */
    quint64 version() const
    {
        return vrsn;
    }

    quint64 version(const PropertyName &key) const;

    /** Return the info object that contains all of the
        metainformation about the atoms, residues, chains,
        cutgroups and segments that make up this molecule,
        and also gives the name and number of this molecule */
    const MoleculeInfoData& info() const
    {
        return *molinfo;
    }

    bool hasProperty(const PropertyName &key) const;
    bool hasMetadata(const PropertyName &metakey) const;
    bool hasMetadata(const PropertyName &key,
                     const PropertyName &metakey) const;
                     
    template<class T>
    bool hasPropertyOfType(const PropertyName &key) const;
                     
    template<class T>
    bool hasMetadataOfType(const PropertyName &metakey) const;
    
    template<class T>
    bool hasMetadataOfType(const PropertyName &key,
                           const PropertyName &metakey) const;
                     
    const char* propertyType(const PropertyName &key) const;
    const char* metadataType(const PropertyName &metakey) const;
    const char* metadataType(const PropertyName &key,
                             const PropertyName &metakey) const;

    /** Return all of the properties of this molecule - this
        includes the coordinates, connectivity and any
        forcefield parameters that have been assigned to
        this molecule */
    const Properties& properties() const
    {
        return props;
    }

    MoleculeData extract(const AtomSelection &selected_atoms) const;

    QStringList propertyKeys() const;
    
    QStringList metadataKeys() const;
    QStringList metadataKeys(const PropertyName &key) const;

    const Property& property(const PropertyName &key) const;
    
    const Property& property(const PropertyName &key,
                             const Property &default_value) const;
    
    const Property& metadata(const PropertyName &metakey) const;

    const Property& metadata(const PropertyName &key,
                             const PropertyName &metakey) const;

    const Property& metadata(const PropertyName &metakey,
                             const Property &default_value) const;
                             
    const Property& metadata(const PropertyName &key,
                             const PropertyName &metakey,
                             const Property &default_value) const;

    void rename(const MolName &newname);

    void renumber();
    void renumber(MolNum newnum);

    void rename(AtomIdx atomidx, const AtomName &newname);
    void rename(const AtomID &atomid, const AtomName &newname);

    void rename(CGIdx cgidx, const CGName &newname);
    void rename(const CGID &cgid, const CGName &newname);
    
    void rename(ResIdx residx, const ResName &newname);
    void rename(const ResID &resid, const ResName &newname);
    
    void rename(ChainIdx chainidx, const ChainName &newname);
    void rename(const ChainID &chainid, const ChainName &newname);
    
    void rename(SegIdx segix, const SegName &newname);
    void rename(const SegID &segid, const SegName &newname);
    
    void renumber(AtomIdx atomidx, AtomNum newnum);
    void renumber(const AtomID &atomid, AtomNum newnum);
    
    void renumber(ResIdx residx, ResNum newnum);
    void renumber(const ResID &resid, ResNum newnum);

    void renumber(const QHash<AtomNum,AtomNum> &atomnums);
    void renumber(const QHash<ResNum,ResNum> &resnums);
    void renumber(const QHash<AtomNum,AtomNum> &atomnums,
                  const QHash<ResNum,ResNum> &resnums);

    void setProperty(const QString &key, 
                     const Property &value, bool clear_metadata=false);

    void removeProperty(const QString &key);

    void setMetadata(const QString &metakey, const Property &value);
    void setMetadata(const QString &key, const QString &metakey, 
                     const Property &value);

    void removeMetadata(const QString &metakey);
    void removeMetadata(const QString &key, const QString &metakey);

    PropertyPtr takeProperty(const QString &key);
    PropertyPtr takeMetadata(const QString &metakey);
    PropertyPtr takeMetadata(const QString &key, const QString &metakey);

    void assertHasProperty(const PropertyName &key) const;
    void assertHasMetadata(const PropertyName &metakey) const;
    void assertHasMetadata(const PropertyName &key,
                           const PropertyName &metakey) const;

    void updatePropertyMolInfo();
    void updatePropertyMolInfo(const AtomMatcher &matcher);

    /** Return the shared null MoleculeData */
    static SireBase::SharedDataPointer<MoleculeData> null();

private:
    /** The metainfo about the molecule - this contains the names of the molecule,
        residue and all atoms, and additional metainfo about all of the residues
        and atoms. This object may also be used to map from atom or residue IDs
        to CGAtomIDs (which are used to lookup the coordinates) */
    SireBase::SharedDataPointer<MoleculeInfoData> molinfo;

    /** All of the properties of this molecule - this includes
        the coordinates of the atoms, their connectivity and
        all of the molecule's forcefield parameters (for all forcefields that
        it is to be added to) and all additional properties that may need to
        be added to the molecule (together with their versions held
        as a metadata property "version" with each property) */
    Properties props;

    /** The version number of this molecule - this changes 
        whenever the molecule is changed in any way. If two molecules
        have the same molecule number and version then they must be the same */
    quint64 vrsn;

    /** The name of this molecule. This can be anything you want! */
    MolName molname;

    /** The ID number of this molecule. This is used to identify a molecule
        in a group */
    MolNum molnum;

    friend class MolNum;
    static MolNum createUniqueMolNum();

    class PropVersions
    {
    public:
        PropVersions() : version(0)
        {}
        
        ~PropVersions()
        {}
        
        quint64 increment();
        quint64 increment(const QString &key, quint64 &mol);
        
        void incrementAll(MoleculeData &moldata);
        
        quint64 reset(QHash<QString,quint64> &prop_vrsns);
        
    private:
        quint64 increment(const QString &key);
    
        /** Mutex used to serialise access to the last version
            number */
        QMutex mutex;
        
        /** The last version number assigned to 
            this molecule */ 
        quint64 version;
        
        /** The last version number assigned to
            each property of the molecule */
        QHash<QString,quint64> property_version;
    };

    static QHash< MolNum, boost::weak_ptr<PropVersions> > version_registry;
    static QMutex version_registry_mutex;

    static boost::shared_ptr<PropVersions> registerMolecule(MolNum molnum);

    /** The version number of each of the properties in 
        this molecule */
    QHash<QString,quint64> prop_vrsns;

    /** The incremints that are used to update the version numbers
        of the different parts of this molecule */
    boost::shared_ptr<PropVersions> vrsns;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Return whether or not this molecule has a property called 'key'
    that is of type 'T' */
template<class T>
SIRE_INLINE_TEMPLATE
bool MoleculeData::hasPropertyOfType(const PropertyName &key) const
{
    return props.hasPropertyOfType<T>(key);
}
                 
/** Return whether or not this molecule has some metadata at metakey
    'metakey' that is of type 'T' */
template<class T>
SIRE_INLINE_TEMPLATE
bool MoleculeData::hasMetadataOfType(const PropertyName &metakey) const
{
    return props.hasMetadataOfType<T>(metakey);
}

/** Return whether or not the property at key 'key' has some metadata
    at metakey 'metakey' that is of type 'T'
    
    \throw SireBase::missing_property
*/
template<class T>
SIRE_INLINE_TEMPLATE
bool MoleculeData::hasMetadataOfType(const PropertyName &key,
                                     const PropertyName &metakey) const
{
    return props.hasMetadataOfType<T>(key, metakey);
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

Q_DECLARE_METATYPE(SireMol::MoleculeData)

SIRE_END_HEADER

#endif
