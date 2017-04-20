/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2008  Christopher Woods
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

#include <QHash>
#include <QMutex>

#include "structureeditor.h"

#include "atom.h"
#include "cutgroup.h"
#include "residue.h"
#include "chain.h"
#include "segment.h"
#include "molecule.h"
#include "mover.hpp"
#include "selector.hpp"

#include "atomeditor.h"
#include "cgeditor.h"
#include "reseditor.h"
#include "chaineditor.h"
#include "segeditor.h"
#include "moleditor.h"

#include "atomname.h"
#include "atomidx.h"
#include "atomnum.h"

#include "cgname.h"
#include "cgidx.h"

#include "resname.h"
#include "resnum.h"
#include "residx.h"

#include "chainname.h"
#include "chainidx.h"

#include "segname.h"
#include "segidx.h"

#include "molname.h"
#include "molnum.h"

#include "atommatcher.h"

#include "atomproperty.hpp"
#include "cgproperty.hpp"
#include "resproperty.hpp"
#include "chainproperty.hpp"
#include "segproperty.hpp"

#include "SireBase/properties.h"

#include "tostring.h"

#include <QDebug>
#include <QReadWriteLock>

#include "SireMol/errors.h"
#include "SireBase/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMol;
using namespace SireBase;
using namespace SireID;
using namespace SireStream;

using boost::tuple;

/////////
///////// Implementation of name_cache
/////////

typedef QHash<QString,QString> NameCacheType;
Q_GLOBAL_STATIC( NameCacheType, getNameCache )
Q_GLOBAL_STATIC( QReadWriteLock, getNameCacheLock )

/** This function is used to cache all name strings of molecules.
    This is useful, as most names used in molecules are repeated many
    times (think about how many "O" atoms there are in a box of water!).
    This allows each equal string to share the same data, rather than
    have multiple copies */
QString SIREMOL_EXPORT SireMol::cacheName(const QString &name)
{
    if (name.isEmpty())
        return name;

    NameCacheType *name_cache = getNameCache();
    QReadWriteLock *lock = getNameCacheLock();

    QReadLocker r_lkr(lock);

    if (name_cache->contains(name))
        return name_cache->value(name);
 
    r_lkr.unlock();
    
    QWriteLocker w_lkr(lock);
    
    if (name_cache->contains(name))
        return name_cache->value(name);
    else
    {
        name_cache->insert(name, name);
        return name;
    }
}

/////////
///////// Implementation of detail::EditMolData
/////////
namespace SireMol
{
namespace detail
{

/** This class holds the editable data of an Atom */
class EditAtomData
{
public:
    EditAtomData();
    EditAtomData(const MoleculeInfoData &molinfo, AtomIdx i);
    
    EditAtomData(const EditAtomData &other);
    
    ~EditAtomData();
    
    AtomName name;
    AtomNum number;
    
    quint32 cg_parent;
    quint32 res_parent;
    quint32 seg_parent;
    
    QHash<QString,QVariant> properties;

    QHash<QString,QVariant> molecule_metadata;
    QHash< QString,QHash<QString,QVariant> > property_metadata;
};

/** This class holds the editable data of a CutGroup */
class EditCGData
{
public:
    EditCGData();
    EditCGData(const MoleculeInfoData &molinfo, CGIdx i, 
               const EditMolData &editmol);

    EditCGData(const EditCGData &other);
    
    ~EditCGData();
    
    CGName name;
    
    QList<quint32> atoms;
    
    QHash<QString,QVariant> properties;

    QHash<QString,QVariant> molecule_metadata;
    QHash< QString,QHash<QString,QVariant> > property_metadata;
};

/** This class holds the editable data of a Residue */
class EditResData
{
public:
    EditResData();
    EditResData(const MoleculeInfoData &molinfo, ResIdx residx, 
                const EditMolData &editmol);
    
    EditResData(const EditResData &other);
    
    ~EditResData();
    
    ResName name;
    ResNum number;
    
    quint32 chain_parent;
    
    QList<quint32> atoms;
    
    QHash<QString,QVariant> properties;

    QHash<QString,QVariant> molecule_metadata;
    QHash< QString,QHash<QString,QVariant> > property_metadata;
};

/** This class holds the editable data of a Chain */
class EditChainData
{
public:
    EditChainData();
    EditChainData(const MoleculeInfoData &molinfo, ChainIdx chainidx,
                  const EditMolData &editmol);
    
    EditChainData(const EditChainData &other);
    
    ~EditChainData();
    
    ChainName name;
    
    QList<quint32> residues;
    
    QHash<QString,QVariant> properties;

    QHash<QString,QVariant> molecule_metadata;
    QHash< QString,QHash<QString,QVariant> > property_metadata;
};

/** This class holds the editable data of a Segment */
class EditSegData
{
public:
    EditSegData();
    EditSegData(const MoleculeInfoData &molinfo, SegIdx segidx,
                const EditMolData &editmol);
    
    EditSegData(const EditSegData &other);
    
    ~EditSegData();
    
    SegName name;
    
    QList<quint32> atoms;
    
    QHash<QString,QVariant> properties;

    QHash<QString,QVariant> molecule_metadata;
    QHash< QString,QHash<QString,QVariant> > property_metadata;
};

/** This private class is used to hold the explicitly shared
    data of the StructureEditor. */
class EditMolData
{
public:
    EditMolData();
    
    EditMolData(const MoleculeData &moldata);
    
    EditMolData(const EditMolData &other);
    
    ~EditMolData();
    
    const EditAtomData& atom(quint32 uid) const;
    const EditResData& residue(quint32 uid) const;
    const EditCGData& cutGroup(quint32 uid) const;
    const EditChainData& chain(quint32 uid) const;
    const EditSegData& segment(quint32 uid) const;

    EditAtomData& atom(quint32 uid);
    EditResData& residue(quint32 uid);
    EditCGData& cutGroup(quint32 uid);
    EditChainData& chain(quint32 uid);
    EditSegData& segment(quint32 uid);
    
    CGAtomIdx cgAtomIdx(quint32 atomuid, const EditAtomData &atom) const;
    CGIdx cgIdx(const EditAtomData &atom) const;
    ResIdx resIdx(const EditAtomData &atom) const;
    SegIdx segIdx(const EditAtomData &atom) const;
    
    ChainIdx chainIdx(const EditResData &residue) const;
    
    QList<AtomIdx> atomIdxsFromUIDs(const QList<quint32> &uids) const;
    QList<ResIdx> resIdxsFromUIDs(const QList<quint32> &uids) const;
    
    quint32 getNewUID();
    
    AtomVariantProperty atomProperty(const QString &key) const;
    AtomVariantProperty atomMetadata(const QString &key) const;
    AtomVariantProperty atomMetadata(const QString &key,
                                     const QString &metakey) const;
                                              
    CGVariantProperty cgProperty(const QString &key) const;
    CGVariantProperty cgMetadata(const QString &metakey) const;
    CGVariantProperty cgMetadata(const QString &key, 
                                 const QString &metakey) const;
                                              
    ResVariantProperty resProperty(const QString &key) const;
    ResVariantProperty resMetadata(const QString &metakey) const;
    ResVariantProperty resMetadata(const QString &key, 
                                   const QString &metakey) const;
                                              
    ChainVariantProperty chainProperty(const QString &key) const;
    ChainVariantProperty chainMetadata(const QString &metakey) const;
    ChainVariantProperty chainMetadata(const QString &key, 
                                       const QString &metakey) const;
                                              
    SegVariantProperty segProperty(const QString &key) const;
    SegVariantProperty segMetadata(const QString &metakey) const;
    SegVariantProperty segMetadata(const QString &key, 
                                   const QString &metakey) const;
    
    void assertHasAtomProperty(const QString &key) const;
    void assertHasAtomMetadata(const QString &metakey) const;
    void assertHasAtomMetadata(const QString &key, const QString &metakey) const;

    void assertHasCGProperty(const QString &key) const;
    void assertHasCGMetadata(const QString &metakey) const;
    void assertHasCGMetadata(const QString &key, const QString &metakey) const;

    void assertHasResProperty(const QString &key) const;
    void assertHasResMetadata(const QString &metakey) const;
    void assertHasResMetadata(const QString &key, const QString &metakey) const;

    void assertHasChainProperty(const QString &key) const;
    void assertHasChainMetadata(const QString &metakey) const;
    void assertHasChainMetadata(const QString &key, const QString &metakey) const;

    void assertHasSegProperty(const QString &key) const;
    void assertHasSegMetadata(const QString &metakey) const;
    void assertHasSegMetadata(const QString &key, const QString &metakey) const;
    
    QHash<AtomIdx,AtomIdx> getOldToNewAtomMapping() const;
    
    MolName molname;
    MolNum molnum;
    
    QHash<quint32,EditAtomData> atoms;
    QHash<quint32,EditCGData> cutgroups;
    QHash<quint32,EditResData> residues;
    QHash<quint32,EditChainData> chains;
    QHash<quint32,EditSegData> segments;
    
    QList<quint32> atoms_by_index;
    QList<quint32> cg_by_index;
    QList<quint32> res_by_index;
    QList<quint32> chains_by_index;
    QList<quint32> seg_by_index; 
    
    /** This stores the AtomIdx of each atom when the below
        properties object was constructed. This allows properties
        to be rebuilt by mapping from the old AtomIdx to the new AtomIdx */
    QHash<quint32, AtomIdx> old_atomidxs;
    
    Properties properties;

    SharedDataPointer<MoleculeInfoData> cached_molinfo;

    quint32 last_uid;

private:
    void extractProperties(const Properties &properties);

    void extractProperty(const QString &key, const AtomProp &atom_property);
    void extractProperty(const QString &key, const ResProp &res_property);
    void extractProperty(const QString &key, const CGProp &cg_property);
    void extractProperty(const QString &key, const ChainProp &chain_property);
    void extractProperty(const QString &key, const SegProp &seg_property);
    void extractProperty(const QString &key, const AtomSelection &selected_atoms);

    void extractMetadata(const QString &key, const AtomProp &atom_property);
    void extractMetadata(const QString &key, const ResProp &res_property);
    void extractMetadata(const QString &key, const CGProp &cg_property);
    void extractMetadata(const QString &key, const ChainProp &chain_property);
    void extractMetadata(const QString &key, const SegProp &seg_property);
    void extractMetadata(const QString &key, const AtomSelection &selected_atoms);

    void extractMetadata(const QString &key, const QString &metakey,
                         const AtomProp &atom_property);
    void extractMetadata(const QString &key, const QString &metakey,
                         const ResProp &res_property);
    void extractMetadata(const QString &key, const QString &metakey,
                         const CGProp &cg_property);
    void extractMetadata(const QString &key, const QString &metakey,
                         const ChainProp &chain_property);
    void extractMetadata(const QString &key, const QString &metakey,
                         const SegProp &seg_property);
    void extractMetadata(const QString &key, const QString &metakey,
                         const AtomSelection &selected_atoms);
};

} // end of namespace detail
} // end of namespace SireMol

using namespace SireMol::detail;

/////////
///////// Implementation of EditAtomData
/////////

QDataStream& operator<<(QDataStream &ds, const EditAtomData &editatom)
{
    ds << editatom.name << editatom.number
       << editatom.cg_parent << editatom.res_parent
       << editatom.seg_parent << editatom.properties
       << editatom.molecule_metadata << editatom.property_metadata;
       
    return ds;
}

QDataStream& operator>>(QDataStream &ds, EditAtomData &editatom)
{
    ds >> editatom.name >> editatom.number
       >> editatom.cg_parent >> editatom.res_parent
       >> editatom.seg_parent >> editatom.properties
       >> editatom.molecule_metadata >> editatom.property_metadata;
       
    return ds;
}

EditAtomData::EditAtomData()
             : name( QString::null ), number( AtomNum::null() ),
               cg_parent(0), res_parent(0), seg_parent(0)
{}
   
EditAtomData::EditAtomData(const MoleculeInfoData &molinfo, AtomIdx i)
             : name(molinfo.name(i)), number(molinfo.number(i)),
               cg_parent(0), res_parent(0), seg_parent(0)
{}
   
EditAtomData::EditAtomData(const EditAtomData &other)
             : name(other.name), number(other.number),
               cg_parent(other.cg_parent), res_parent(other.res_parent),
               seg_parent(other.seg_parent), properties(other.properties),
               molecule_metadata(other.molecule_metadata),
               property_metadata(other.property_metadata)
{}
    
EditAtomData::~EditAtomData()
{}

/////////
///////// Implementation of EditCGData
/////////

QDataStream& operator<<(QDataStream &ds, const EditCGData &editcg)
{
    ds << editcg.name << editcg.atoms << editcg.properties
       << editcg.molecule_metadata << editcg.property_metadata;
       
    return ds;
}

QDataStream& operator>>(QDataStream &ds, EditCGData &editcg)
{
    ds >> editcg.name >> editcg.atoms >> editcg.properties
       >> editcg.molecule_metadata >> editcg.property_metadata;
       
    return ds;
}

EditCGData::EditCGData() : name( QString::null )
{}

EditCGData::EditCGData(const MoleculeInfoData &molinfo, CGIdx i,
                       const EditMolData &editmol)
           : name(molinfo.name(i))
{
    foreach (AtomIdx atomidx, molinfo.getAtomsIn(i))
    {
        quint32 uid = editmol.atoms_by_index.at(atomidx);
        atoms.append(uid);
    }
}

EditCGData::EditCGData(const EditCGData &other)
           : name(other.name), atoms(other.atoms), properties(other.properties),
             molecule_metadata(other.molecule_metadata),
             property_metadata(other.property_metadata)
{}
    
EditCGData::~EditCGData()
{}

/////////
///////// Implementation of EditResData
/////////

QDataStream& operator<<(QDataStream &ds, const EditResData &editres)
{
    ds << editres.name << editres.number << editres.chain_parent
       << editres.atoms << editres.properties
       << editres.molecule_metadata << editres.property_metadata;
       
    return ds;
}

QDataStream& operator>>(QDataStream &ds, EditResData &editres)
{
    ds >> editres.name >> editres.number >> editres.chain_parent
       >> editres.atoms >> editres.properties
       >> editres.molecule_metadata >> editres.property_metadata;
       
    return ds;
}

EditResData::EditResData()
            : name( QString::null ), number( ResNum::null() ),
              chain_parent(0)
{}

EditResData::EditResData(const MoleculeInfoData &molinfo, ResIdx i,
                         const EditMolData &editmol)
            : name(molinfo.name(i)), number(molinfo.number(i)), chain_parent(0)
{
    foreach (AtomIdx atomidx, molinfo.getAtomsIn(i))
    {
        quint32 uid = editmol.atoms_by_index.at(atomidx);
        atoms.append(uid);
    }
}

EditResData::EditResData(const EditResData &other)
            : name(other.name), number(other.number), 
              chain_parent(other.chain_parent), atoms(other.atoms),
              properties(other.properties),
              molecule_metadata(other.molecule_metadata),
              property_metadata(other.property_metadata)
{}
    
EditResData::~EditResData()
{}

/////////
///////// Implementation of EditChainData
/////////

QDataStream& operator<<(QDataStream &ds, const EditChainData &editchain)
{
    ds << editchain.name << editchain.residues << editchain.properties
       << editchain.molecule_metadata << editchain.property_metadata;
    
    return ds;
}

QDataStream& operator>>(QDataStream &ds, EditChainData &editchain)
{
    ds >> editchain.name >> editchain.residues >> editchain.properties
       >> editchain.molecule_metadata >> editchain.property_metadata;
    
    return ds;
}

EditChainData::EditChainData() : name( QString::null )
{}

EditChainData::EditChainData(const MoleculeInfoData &molinfo, ChainIdx i,
                             const EditMolData &editmol)
              : name(molinfo.name(i))
{
    foreach (ResIdx residx, molinfo.getResiduesIn(i))
    {
        quint32 uid = editmol.res_by_index.at(residx);
        residues.append(uid);
    }
}

EditChainData::EditChainData(const EditChainData &other)
              : name(other.name), residues(other.residues),
                properties(other.properties),
                molecule_metadata(other.molecule_metadata),
                property_metadata(other.property_metadata)
{}
    
EditChainData::~EditChainData()
{}
    
/////////
///////// Implementation of EditSegData
/////////

QDataStream& operator<<(QDataStream &ds, const EditSegData &editseg)
{
    ds << editseg.name << editseg.atoms << editseg.properties
       << editseg.molecule_metadata << editseg.property_metadata;
    
    return ds;
}

QDataStream& operator>>(QDataStream &ds, EditSegData &editseg)
{
    ds >> editseg.name >> editseg.atoms >> editseg.properties
       >> editseg.molecule_metadata >> editseg.property_metadata;
    
    return ds;
}

EditSegData::EditSegData() : name( QString::null )
{}

EditSegData::EditSegData(const MoleculeInfoData &molinfo, SegIdx i,
                         const EditMolData &editmol)
            : name(molinfo.name(i))
{
    foreach (AtomIdx atomidx, molinfo.getAtomsIn(i))
    {
        quint32 uid = editmol.atoms_by_index.at(atomidx);
        atoms.append(uid);
    }
}

EditSegData::EditSegData(const EditSegData &other)
            : name(other.name), atoms(other.atoms),
              properties(other.properties),
              molecule_metadata(other.molecule_metadata),
              property_metadata(other.property_metadata)
{}
    
EditSegData::~EditSegData()
{}

/////////
///////// Implementation of EditMolData
/////////

/** Serialise to a binary datastream */
QDataStream& operator<<(QDataStream &ds, const EditMolData &editmol)
{
    SharedDataStream sds(ds);

    sds << editmol.molname << editmol.molnum
        << editmol.atoms << editmol.cutgroups << editmol.residues
        << editmol.chains << editmol.segments
        << editmol.atoms_by_index << editmol.cg_by_index
        << editmol.res_by_index << editmol.chains_by_index
        << editmol.seg_by_index
        << editmol.old_atomidxs
        << editmol.properties << editmol.cached_molinfo
        << editmol.last_uid;
       
    return ds;
}

/** Extract from a binary datastream */
QDataStream& operator>>(QDataStream &ds, EditMolData &editmol)
{
    SharedDataStream sds(ds);

    sds >> editmol.molname >> editmol.molnum
        >> editmol.atoms >> editmol.cutgroups >> editmol.residues
        >> editmol.chains >> editmol.segments
        >> editmol.atoms_by_index >> editmol.cg_by_index
        >> editmol.res_by_index >> editmol.chains_by_index
        >> editmol.seg_by_index
        >> editmol.old_atomidxs
        >> editmol.properties >> editmol.cached_molinfo
        >> editmol.last_uid;
       
    return ds;
}

/** Constructor */
EditMolData::EditMolData() : molname( QString::null ), 
                             molnum( MolNum::null() ), last_uid(0)
{}

/** Copy the properties for each atom to the individual EditAtomData objects */
void EditMolData::extractProperty(const QString &key, 
                                  const AtomProp &atom_property)
{
    //convert the properties for each atom into an array of array
    //of QVariants. These are arranged in CutGroups, in CGAtomIdx order...
    PackedArray2D<QVariant> values = atom_property.toVariant().array();
    
    int ngroups = values.count();
    BOOST_ASSERT( ngroups == cg_by_index.count() );
    
    const PackedArray2D<QVariant>::Array *values_array = values.constData();
    
    for (int i=0; i<ngroups; ++i)
    {
        const PackedArray2D<QVariant>::Array &group_values = values_array[i];
    
        //get the CutGroup at index i
        const EditCGData &cgroup = this->cutGroup(cg_by_index.at(i));
        
        //now loop over the atoms...
        int nats = group_values.count();
        BOOST_ASSERT( nats == cgroup.atoms.count() );
        const QVariant *group_values_array = group_values.constData();
        
        for (int j=0; j<nats; ++j)
        {
            EditAtomData &atom = this->atom(cgroup.atoms.at(j));
            
            atom.properties.insert(key, group_values_array[j]);
        }
    }
}

/** Copy the properties for each residue to the individual EditResData objects */
void EditMolData::extractProperty(const QString &key, 
                                  const ResProp &res_property)
{
    //convert each property to a QVariant, in ResIdx order
    ResVariantProperty values = res_property.toVariant();
    
    int nres = values.count();
    BOOST_ASSERT( nres == res_by_index.count() );
    
    const QVariant *values_array = values.constData();

    for (int i=0; i<nres; ++i)
    {
        EditResData &residue = this->residue(res_by_index.at(i));
        residue.properties.insert(key, values_array[i]);
    }
}

/** Copy the properties for each atom to the individual EditCGData objects */
void EditMolData::extractProperty(const QString &key, 
                                  const CGProp &cg_property)
{
    //convert each property to a QVariant, in CGIdx order
    CGVariantProperty values = cg_property.toVariant();
    
    int ncg = values.count();
    BOOST_ASSERT( ncg == cg_by_index.count() );
    
    const QVariant *values_array = values.constData();

    for (int i=0; i<ncg; ++i)
    {
        EditCGData &cgroup = this->cutGroup(cg_by_index.at(i));
        cgroup.properties.insert(key, values_array[i]);
    }
}

/** Copy the properties for each atom to the individual EditChainData objects */
void EditMolData::extractProperty(const QString &key, 
                                  const ChainProp &chain_property)
{
    //convert each property to a QVariant, in ChainIdx order
    ChainVariantProperty values = chain_property.toVariant();
    
    int nchains = values.count();
    BOOST_ASSERT( nchains == chains_by_index.count() );
    
    const QVariant *values_array = values.constData();

    for (int i=0; i<nchains; ++i)
    {
        EditChainData &chain = this->chain(chains_by_index.at(i));
        chain.properties.insert(key, values_array[i]);
    }
}

/** Copy the properties for each atom to the individual EditSegData objects */
void EditMolData::extractProperty(const QString &key, 
                                  const SegProp &seg_property)
{
    //convert each property to a QVariant, in ResIdx order
    SegVariantProperty values = seg_property.toVariant();
    
    int nseg = values.count();
    BOOST_ASSERT( nseg == seg_by_index.count() );
    
    const QVariant *values_array = values.constData();

    for (int i=0; i<nseg; ++i)
    {
        EditSegData &segment = this->segment(seg_by_index.at(i));
        segment.properties.insert(key, values_array[i]);
    }
}

/** Copy the properties for each atom to the individual EditAtomData objects */
void EditMolData::extractProperty(const QString &key, 
                                  const AtomSelection &selected_atoms)
{
    if (selected_atoms.selectedAll())
    {
        //select every atom
        for (QHash<quint32,EditAtomData>::iterator it = atoms.begin();
             it != atoms.end();
             ++it)
        {
            it->properties.insert(key, QVariant(true));
        }
    }
    else if (selected_atoms.selectedNone())
    {
        //deselect every atom
        for (QHash<quint32,EditAtomData>::iterator it = atoms.begin();
             it != atoms.end();
             ++it)
        {
            it->properties.insert(key, QVariant(false));
        }
    }
    else
    {
        //we need to go through each atom in turn...
        int nats = atoms_by_index.count();
        
        for (AtomIdx i(0); i<nats; ++i)
        {
            EditAtomData &atom = this->atom( atoms_by_index.at(i) );
            
            atom.properties.insert(key, selected_atoms.selected(i));
        }
    }
}

/** Copy the metadata for each atom to the individual EditAtomData objects */
void EditMolData::extractMetadata(const QString &key, 
                                  const AtomProp &atom_property)
{
    //convert the properties for each atom into an array of array
    //of QVariants. These are arranged in CutGroups, in CGAtomIdx order...
    PackedArray2D<QVariant> values = atom_property.toVariant().array();
    
    int ngroups = values.count();
    BOOST_ASSERT( ngroups == cg_by_index.count() );
    
    const PackedArray2D<QVariant>::Array *values_array = values.constData();
    
    for (int i=0; i<ngroups; ++i)
    {
        const PackedArray2D<QVariant>::Array &group_values = values_array[i];
    
        //get the CutGroup at index i
        const EditCGData &cgroup = this->cutGroup(cg_by_index.at(i));
        
        //now loop over the atoms...
        int nats = group_values.count();
        BOOST_ASSERT( nats == cgroup.atoms.count() );
        const QVariant *group_values_array = group_values.constData();
        
        for (int j=0; j<nats; ++j)
        {
            EditAtomData &atom = this->atom(cgroup.atoms.at(j));
            
            atom.molecule_metadata.insert(key, group_values_array[j]);
        }
    }
}

/** Copy the metadata for each residue to the individual EditResData objects */
void EditMolData::extractMetadata(const QString &key, 
                                  const ResProp &res_property)
{
    //convert each property to a QVariant, in ResIdx order
    ResVariantProperty values = res_property.toVariant();
    
    int nres = values.count();
    BOOST_ASSERT( nres == res_by_index.count() );
    
    const QVariant *values_array = values.constData();

    for (int i=0; i<nres; ++i)
    {
        EditResData &residue = this->residue(res_by_index.at(i));
        residue.molecule_metadata.insert(key, values_array[i]);
    }
}

/** Copy the metadata for each atom to the individual EditCGData objects */
void EditMolData::extractMetadata(const QString &key, 
                                  const CGProp &cg_property)
{
    //convert each property to a QVariant, in CGIdx order
    CGVariantProperty values = cg_property.toVariant();
    
    int ncg = values.count();
    BOOST_ASSERT( ncg == cg_by_index.count() );
    
    const QVariant *values_array = values.constData();

    for (int i=0; i<ncg; ++i)
    {
        EditCGData &cgroup = this->cutGroup(cg_by_index.at(i));
        cgroup.molecule_metadata.insert(key, values_array[i]);
    }
}

/** Copy the metadata for each atom to the individual EditChainData objects */
void EditMolData::extractMetadata(const QString &key, 
                                  const ChainProp &chain_property)
{
    //convert each property to a QVariant, in ChainIdx order
    ChainVariantProperty values = chain_property.toVariant();
    
    int nchains = values.count();
    BOOST_ASSERT( nchains == chains_by_index.count() );
    
    const QVariant *values_array = values.constData();

    for (int i=0; i<nchains; ++i)
    {
        EditChainData &chain = this->chain(chains_by_index.at(i));
        chain.molecule_metadata.insert(key, values_array[i]);
    }
}

/** Copy the metadata for each atom to the individual EditSegData objects */
void EditMolData::extractMetadata(const QString &key, 
                                  const SegProp &seg_property)
{
    //convert each property to a QVariant, in ResIdx order
    SegVariantProperty values = seg_property.toVariant();
    
    int nseg = values.count();
    BOOST_ASSERT( nseg == seg_by_index.count() );
    
    const QVariant *values_array = values.constData();

    for (int i=0; i<nseg; ++i)
    {
        EditSegData &segment = this->segment(seg_by_index.at(i));
        segment.molecule_metadata.insert(key, values_array[i]);
    }
}

/** Copy the properties for each atom to the individual EditAtomData objects */
void EditMolData::extractMetadata(const QString &key, 
                                  const AtomSelection &selected_atoms)
{
    if (selected_atoms.selectedAll())
    {
        //select every atom
        for (QHash<quint32,EditAtomData>::iterator it = atoms.begin();
             it != atoms.end();
             ++it)
        {
            it->molecule_metadata.insert(key, QVariant(true));
        }
    }
    else if (selected_atoms.selectedNone())
    {
        //deselect every atom
        for (QHash<quint32,EditAtomData>::iterator it = atoms.begin();
             it != atoms.end();
             ++it)
        {
            it->molecule_metadata.insert(key, QVariant(false));
        }
    }
    else
    {
        //we need to go through each atom in turn...
        int nats = atoms_by_index.count();
        
        for (AtomIdx i(0); i<nats; ++i)
        {
            EditAtomData &atom = this->atom( atoms_by_index.at(i) );
            
            atom.molecule_metadata.insert(key, selected_atoms.selected(i));
        }
    }
}

/** Copy the metadata for each atom to the individual EditAtomData objects */
void EditMolData::extractMetadata(const QString &key, const QString &metakey,
                                  const AtomProp &atom_property)
{
    //convert the properties for each atom into an array of array
    //of QVariants. These are arranged in CutGroups, in CGAtomIdx order...
    PackedArray2D<QVariant> values = atom_property.toVariant().array();
    
    int ngroups = values.count();
    BOOST_ASSERT( ngroups == cg_by_index.count() );
    
    const PackedArray2D<QVariant>::Array *values_array = values.constData();
    
    for (int i=0; i<ngroups; ++i)
    {
        const PackedArray2D<QVariant>::Array &group_values = values_array[i];
    
        //get the CutGroup at index i
        const EditCGData &cgroup = this->cutGroup(cg_by_index.at(i));
        
        //now loop over the atoms...
        int nats = group_values.count();
        BOOST_ASSERT( nats == cgroup.atoms.count() );
        const QVariant *group_values_array = group_values.constData();
        
        for (int j=0; j<nats; ++j)
        {
            EditAtomData &atom = this->atom(cgroup.atoms.at(j));
            
            atom.property_metadata[key].insert(metakey, group_values_array[j]);
        }
    }
}

/** Copy the metadata for each residue to the individual EditResData objects */
void EditMolData::extractMetadata(const QString &key, const QString &metakey,
                                  const ResProp &res_property)
{
    //convert each property to a QVariant, in ResIdx order
    ResVariantProperty values = res_property.toVariant();
    
    int nres = values.count();
    BOOST_ASSERT( nres == res_by_index.count() );
    
    const QVariant *values_array = values.constData();

    for (int i=0; i<nres; ++i)
    {
        EditResData &residue = this->residue(res_by_index.at(i));
        residue.property_metadata[key].insert(metakey, values_array[i]);
    }
}

/** Copy the metadata for each atom to the individual EditCGData objects */
void EditMolData::extractMetadata(const QString &key, const QString &metakey,
                                  const CGProp &cg_property)
{
    //convert each property to a QVariant, in CGIdx order
    CGVariantProperty values = cg_property.toVariant();
    
    int ncg = values.count();
    BOOST_ASSERT( ncg == cg_by_index.count() );
    
    const QVariant *values_array = values.constData();

    for (int i=0; i<ncg; ++i)
    {
        EditCGData &cgroup = this->cutGroup(cg_by_index.at(i));
        cgroup.property_metadata[key].insert(metakey, values_array[i]);
    }
}

/** Copy the metadata for each atom to the individual EditChainData objects */
void EditMolData::extractMetadata(const QString &key, const QString &metakey,
                                  const ChainProp &chain_property)
{
    //convert each property to a QVariant, in ChainIdx order
    ChainVariantProperty values = chain_property.toVariant();
    
    int nchains = values.count();
    BOOST_ASSERT( nchains == chains_by_index.count() );
    
    const QVariant *values_array = values.constData();

    for (int i=0; i<nchains; ++i)
    {
        EditChainData &chain = this->chain(chains_by_index.at(i));
        chain.property_metadata[key].insert(metakey, values_array[i]);
    }
}

/** Copy the metadata for each atom to the individual EditSegData objects */
void EditMolData::extractMetadata(const QString &key, const QString &metakey,
                                  const SegProp &seg_property)
{
    //convert each property to a QVariant, in ResIdx order
    SegVariantProperty values = seg_property.toVariant();
    
    int nseg = values.count();
    BOOST_ASSERT( nseg == seg_by_index.count() );
    
    const QVariant *values_array = values.constData();

    for (int i=0; i<nseg; ++i)
    {
        EditSegData &segment = this->segment(seg_by_index.at(i));
        segment.property_metadata[key].insert(metakey, values_array[i]);
    }
}

/** Copy the properties for each atom to the individual EditAtomData objects */
void EditMolData::extractMetadata(const QString &key, const QString &metakey,
                                  const AtomSelection &selected_atoms)
{
    if (selected_atoms.selectedAll())
    {
        //select every atom
        for (QHash<quint32,EditAtomData>::iterator it = atoms.begin();
             it != atoms.end();
             ++it)
        {
            it->property_metadata[key].insert(metakey, QVariant(true));
        }
    }
    else if (selected_atoms.selectedNone())
    {
        //deselect every atom
        for (QHash<quint32,EditAtomData>::iterator it = atoms.begin();
             it != atoms.end();
             ++it)
        {
            it->property_metadata[key].insert(metakey, QVariant(false));
        }
    }
    else
    {
        //we need to go through each atom in turn...
        int nats = atoms_by_index.count();
        
        for (AtomIdx i(0); i<nats; ++i)
        {
            EditAtomData &atom = this->atom( atoms_by_index.at(i) );
            
            atom.property_metadata[key].insert(metakey, 
                                               selected_atoms.selected(i));
        }
    }
}

/** Extract all of the properties and attach them to the correct
    parts of the molecule */
void EditMolData::extractProperties(const Properties &props)
{
    //first copy the properties to our value
    properties = props;

    for (Properties::const_iterator it = properties.constBegin();
         it != properties.constEnd();
         ++it)
    {
        if (it.value()->isA<AtomProp>())
            this->extractProperty(it.key(), it.value()->asA<AtomProp>());
            
        else if (it.value()->isA<CGProp>())
            this->extractProperty(it.key(), it.value()->asA<CGProp>());
            
        else if (it.value()->isA<ResProp>())
            this->extractProperty(it.key(), it.value()->asA<ResProp>());
            
        else if (it.value()->isA<ChainProp>())
            this->extractProperty(it.key(), it.value()->asA<ChainProp>());
            
        else if (it.value()->isA<SegProp>())
            this->extractProperty(it.key(), it.value()->asA<SegProp>());
            
        else if (it.value()->isA<AtomSelection>())
            this->extractProperty(it.key(), it.value()->asA<AtomSelection>());
                    
        //now do the same for all of the metadata attached to this
        //property
        const Properties &metadata = props.allMetadata(it.key());
        
        for (Properties::const_iterator it2 = metadata.constBegin();
             it2 != metadata.constEnd();
             ++it2)
        {
            if (it2.value()->isA<AtomProp>())
                this->extractMetadata(it.key(), it2.key(),
                                      it2.value()->asA<AtomProp>());
            
            else if (it2.value()->isA<CGProp>())
                this->extractMetadata(it.key(), it2.key(),
                                      it2.value()->asA<CGProp>());
            
            else if (it2.value()->isA<ResProp>())
                this->extractMetadata(it.key(), it2.key(),
                                      it2.value()->asA<ResProp>());
            
            else if (it2.value()->isA<ChainProp>())
                this->extractMetadata(it.key(), it2.key(),
                                      it2.value()->asA<ChainProp>());
            
            else if (it2.value()->isA<SegProp>())
                this->extractMetadata(it.key(), it2.key(),
                                      it2.value()->asA<SegProp>());
            
            else if (it2.value()->isA<AtomSelection>())
                this->extractMetadata(it.key(), it2.key(),
                                      it2.value()->asA<AtomSelection>());
        }
    }
    
    //now we've done the properties and their metadata, it is now time 
    //to extract the molecule's metadata as well (as this may also 
    //be attached to various molecular subgroups)
    const Properties &metadata = properties.allMetadata();
    
    for (Properties::const_iterator it = metadata.constBegin();
         it != metadata.constEnd();
         ++it)
    {
        if (it.value()->isA<AtomProp>())
            this->extractMetadata(it.key(), it.value()->asA<AtomProp>());
            
        else if (it.value()->isA<CGProp>())
            this->extractMetadata(it.key(), it.value()->asA<CGProp>());
            
        else if (it.value()->isA<ResProp>())
            this->extractMetadata(it.key(), it.value()->asA<ResProp>());
            
        else if (it.value()->isA<ChainProp>())
            this->extractMetadata(it.key(), it.value()->asA<ChainProp>());
            
        else if (it.value()->isA<SegProp>())
            this->extractMetadata(it.key(), it.value()->asA<SegProp>());
            
        else if (it.value()->isA<AtomSelection>())
            this->extractMetadata(it.key(), it.value()->asA<AtomSelection>());
    }
}

/** Construct from a MoleculeData */
EditMolData::EditMolData(const MoleculeData &moldata)
            : molname(moldata.name()), molnum(moldata.number()),
              cached_molinfo(moldata.info()), last_uid(0)
{
    const MoleculeInfoData &molinfo = moldata.info();
    
    int nats = molinfo.nAtoms();
    
    for (AtomIdx i(0); i<nats; ++i)
    {
        EditAtomData atom(molinfo, i);
    
        quint32 uid = getNewUID();
        atoms.insert(uid, atom);
        atoms_by_index.append(uid);
        
        old_atomidxs.insert(uid, i);
    }
    
    int nres = molinfo.nResidues();
    
    for (ResIdx i(0); i<nres; ++i)
    {
        EditResData residue(molinfo, i, *this);
        
        quint32 uid = getNewUID();
        residues.insert(uid, residue);
        res_by_index.append(uid);
        
        foreach (quint32 atomuid, residue.atoms)
        {
            this->atom(atomuid).res_parent = uid;
        }
    }
    
    int ncg = molinfo.nCutGroups();
    
    for (CGIdx i(0); i<ncg; ++i)
    {
        EditCGData cgroup(molinfo, i, *this);
        
        quint32 uid = getNewUID();
        cutgroups.insert(uid, cgroup);
        cg_by_index.append(uid);
        
        foreach (quint32 atomuid, cgroup.atoms)
        {
            this->atom(atomuid).cg_parent = uid;
        }
    }
    
    int nchains = molinfo.nChains();
    
    for (ChainIdx i(0); i<nchains; ++i)
    {
        EditChainData chain(molinfo, i, *this);
        
        quint32 uid = getNewUID();
        chains.insert(uid, chain);
        chains_by_index.append(uid);
        
        foreach (quint32 resuid, chain.residues)
        {
            this->residue(resuid).chain_parent = uid;
        }
    }
    
    int nseg = molinfo.nSegments();
    
    for (SegIdx i(0); i<nseg; ++i)
    {
        EditSegData segment(molinfo, i, *this);
        
        quint32 uid = getNewUID();
        segments.insert(uid, segment);
        seg_by_index.append(uid);
        
        foreach (quint32 atomuid, segment.atoms)
        {
            this->atom(atomuid).seg_parent = uid;
        }
    }
    
    //need to convert them...
    this->extractProperties(moldata.properties());
}

/** Copy constructor */
EditMolData::EditMolData(const EditMolData &other)
            : molname(other.molname), molnum(other.molnum),
              atoms(other.atoms), cutgroups(other.cutgroups),
              residues(other.residues), chains(other.chains),
              segments(other.segments),
              atoms_by_index(other.atoms_by_index),
              cg_by_index(other.cg_by_index),
              res_by_index(other.res_by_index),
              chains_by_index(other.chains_by_index),
              seg_by_index(other.seg_by_index),
              old_atomidxs(other.old_atomidxs),
              properties(other.properties),
              cached_molinfo(other.cached_molinfo),
              last_uid(other.last_uid)
{}

/** Destructor */
EditMolData::~EditMolData()
{}

/** Get a new unique ID number */
quint32 EditMolData::getNewUID()
{
    if (last_uid == std::numeric_limits<quint32>::max())
        throw SireError::program_bug( QObject::tr(
            "An EditMolData can only identify %1 unique objects!")
                .arg(std::numeric_limits<quint32>::max()), CODELOC );
                
    ++last_uid;
    
    return last_uid;
}
    
const EditAtomData& EditMolData::atom(quint32 uid) const
{
    QHash<quint32,EditAtomData>::const_iterator it = atoms.find(uid);
    
    if (it == atoms.end())
        throw SireMol::missing_atom( QObject::tr(
            "There is no atom in this molecule that is identified by "
            "the UID %1.").arg(uid), CODELOC );
            
    return it.value();
}

const EditResData& EditMolData::residue(quint32 uid) const
{
    QHash<quint32,EditResData>::const_iterator it = residues.find(uid);
    
    if (it == residues.end())
        throw SireMol::missing_residue( QObject::tr(
            "There is no residue in this molecule that is identified by "
            "the UID %1.").arg(uid), CODELOC );
            
    return it.value();
}

const EditCGData& EditMolData::cutGroup(quint32 uid) const
{
    QHash<quint32,EditCGData>::const_iterator it = cutgroups.find(uid);
    
    if (it == cutgroups.end())
        throw SireMol::missing_cutgroup( QObject::tr(
            "There is no CutGroup in this molecule that is identified by "
            "the UID %1.").arg(uid), CODELOC );
            
    return it.value();
}

const EditChainData& EditMolData::chain(quint32 uid) const
{
    QHash<quint32,EditChainData>::const_iterator it = chains.find(uid);
    
    if (it == chains.end())
        throw SireMol::missing_chain( QObject::tr(
            "There is no chain in this molecule that is identified by "
            "the UID %1.").arg(uid), CODELOC );
            
    return it.value();
}

const EditSegData& EditMolData::segment(quint32 uid) const
{
    QHash<quint32,EditSegData>::const_iterator it = segments.find(uid);
    
    if (it == segments.end())
        throw SireMol::missing_segment( QObject::tr(
            "There is no segment in this molecule that is identified by "
            "the UID %1.").arg(uid), CODELOC );
            
    return it.value();
}

EditAtomData& EditMolData::atom(quint32 uid)
{
    QHash<quint32,EditAtomData>::iterator it = atoms.find(uid);
    
    if (it == atoms.end())
        throw SireMol::missing_atom( QObject::tr(
            "There is no atom in this molecule that is identified by "
            "the UID %1.").arg(uid), CODELOC );
            
    return it.value();
}

EditResData& EditMolData::residue(quint32 uid)
{
    QHash<quint32,EditResData>::iterator it = residues.find(uid);
    
    if (it == residues.end())
        throw SireMol::missing_residue( QObject::tr(
            "There is no residue in this molecule that is identified by "
            "the UID %1.").arg(uid), CODELOC );
            
    return it.value();
}

EditCGData& EditMolData::cutGroup(quint32 uid)
{
    QHash<quint32,EditCGData>::iterator it = cutgroups.find(uid);
    
    if (it == cutgroups.end())
        throw SireMol::missing_cutgroup( QObject::tr(
            "There is no CutGroup in this molecule that is identified by "
            "the UID %1.").arg(uid), CODELOC );
            
    return it.value();
}

EditChainData& EditMolData::chain(quint32 uid)
{
    QHash<quint32,EditChainData>::iterator it = chains.find(uid);
    
    if (it == chains.end())
        throw SireMol::missing_chain( QObject::tr(
            "There is no chain in this molecule that is identified by "
            "the UID %1.").arg(uid), CODELOC );
            
    return it.value();
}

EditSegData& EditMolData::segment(quint32 uid)
{
    QHash<quint32,EditSegData>::iterator it = segments.find(uid);
    
    if (it == segments.end())
        throw SireMol::missing_segment( QObject::tr(
            "There is no segment in this molecule that is identified by "
            "the UID %1.").arg(uid), CODELOC );
            
    return it.value();
}

QList<AtomIdx> EditMolData::atomIdxsFromUIDs(const QList<quint32> &uids) const
{
    QList<AtomIdx> atomidxs;
    QList<quint32> invalid_uids;
    
    foreach (quint32 uid, uids)
    {
        int idx = atoms_by_index.indexOf(uid);
        
        if (idx >= 0)
            atomidxs.append( AtomIdx(idx) );
        else
            invalid_uids.append(uid);
    }
    
    if (not invalid_uids.isEmpty())
        throw SireError::invalid_index( QObject::tr( 
            "There were some invalid Atom UID numbers! %1")
                .arg( Sire::toString(invalid_uids) ), CODELOC );
                
    return atomidxs;
}

QList<ResIdx> EditMolData::resIdxsFromUIDs(const QList<quint32> &uids) const
{
    QList<ResIdx> residxs;
    QList<quint32> invalid_uids;
    
    foreach (quint32 uid, uids)
    {
        int idx = res_by_index.indexOf(uid);
        
        if (idx >= 0)
            residxs.append( ResIdx(idx) );
        else
            invalid_uids.append(uid);
    }
    
    if (not invalid_uids.isEmpty())
        throw SireError::invalid_index( QObject::tr( 
            "There were some invalid residue UID numbers! %1")
                .arg( Sire::toString(invalid_uids) ), CODELOC );
                
    return residxs;
}
 
void EditMolData::assertHasAtomProperty(const QString &key) const
{
    const Property &property = properties.property(key);
    
    if (not property.isA<AtomProp>())
        throw SireError::invalid_cast( QObject::tr(
            "The property at key %1 of type %2 is not an AtomProperty!")
                .arg(key, property.what()), CODELOC );
}

void EditMolData::assertHasAtomMetadata(const QString &metakey) const
{
    const Property &property = properties.metadata(metakey);
    
    if (not property.isA<AtomProp>())
        throw SireError::invalid_cast( QObject::tr(
            "The metadata at metakey %1 of type %2 is not an AtomProperty!")
                .arg(metakey, property.what()), CODELOC );
}

void EditMolData::assertHasAtomMetadata(const QString &key,
                                        const QString &metakey) const
{
    const Property &property = properties.metadata(key, metakey);
    
    if (not property.isA<AtomProp>())
        throw SireError::invalid_cast( QObject::tr(
            "The metadata at key %1, metakey %2 of type %3 is not an AtomProperty!")
                .arg(key, metakey, property.what()), CODELOC );
}

void EditMolData::assertHasCGProperty(const QString &key) const
{
    const Property &property = properties.property(key);
    
    if (not property.isA<CGProp>())
        throw SireError::invalid_cast( QObject::tr(
            "The property at key %1 of type %2 is not a CGProperty!")
                .arg(key, property.what()), CODELOC );
}

void EditMolData::assertHasCGMetadata(const QString &metakey) const
{
    const Property &property = properties.metadata(metakey);
    
    if (not property.isA<CGProp>())
        throw SireError::invalid_cast( QObject::tr(
            "The metadata at metakey %1 of type %2 is not a CGProperty!")
                .arg(metakey, property.what()), CODELOC );
}

void EditMolData::assertHasCGMetadata(const QString &key,
                                      const QString &metakey) const
{
    const Property &property = properties.metadata(key, metakey);
    
    if (not property.isA<CGProp>())
        throw SireError::invalid_cast( QObject::tr(
            "The metadata at key %1, metakey %2 of type %3 is not a CGProperty!")
                .arg(key, metakey, property.what()), CODELOC );
}

void EditMolData::assertHasResProperty(const QString &key) const
{
    const Property &property = properties.property(key);
    
    if (not property.isA<ResProp>())
        throw SireError::invalid_cast( QObject::tr(
            "The property at key %1 of type %2 is not a ResProperty!")
                .arg(key, property.what()), CODELOC );
}

void EditMolData::assertHasResMetadata(const QString &metakey) const
{
    const Property &property = properties.metadata(metakey);
    
    if (not property.isA<ResProp>())
        throw SireError::invalid_cast( QObject::tr(
            "The metadata at metakey %1 of type %2 is not a ResProperty!")
                .arg(metakey, property.what()), CODELOC );
}

void EditMolData::assertHasResMetadata(const QString &key,
                                       const QString &metakey) const
{
    const Property &property = properties.metadata(key, metakey);
    
    if (not property.isA<ResProp>())
        throw SireError::invalid_cast( QObject::tr(
            "The metadata at key %1, metakey %2 of type %3 is not a ResProperty!")
                .arg(key, metakey, property.what()), CODELOC );
}

void EditMolData::assertHasChainProperty(const QString &key) const
{
    const Property &property = properties.property(key);
    
    if (not property.isA<ChainProp>())
        throw SireError::invalid_cast( QObject::tr(
            "The property at key %1 of type %2 is not a ChainProperty!")
                .arg(key, property.what()), CODELOC );
}

void EditMolData::assertHasChainMetadata(const QString &metakey) const
{
    const Property &property = properties.metadata(metakey);
    
    if (not property.isA<ChainProp>())
        throw SireError::invalid_cast( QObject::tr(
            "The metadata at metakey %1 of type %2 is not a ChainProperty!")
                .arg(metakey, property.what()), CODELOC );
}

void EditMolData::assertHasChainMetadata(const QString &key,
                                         const QString &metakey) const
{
    const Property &property = properties.metadata(key, metakey);
    
    if (not property.isA<ChainProp>())
        throw SireError::invalid_cast( QObject::tr(
            "The metadata at key %1, metakey %2 of type %3 is not an ChainProperty!")
                .arg(key, metakey, property.what()), CODELOC );
}

void EditMolData::assertHasSegProperty(const QString &key) const
{
    const Property &property = properties.property(key);
    
    if (not property.isA<SegProp>())
        throw SireError::invalid_cast( QObject::tr(
            "The property at key %1 of type %2 is not a SegProperty!")
                .arg(key, property.what()), CODELOC );
}

void EditMolData::assertHasSegMetadata(const QString &metakey) const
{
    const Property &property = properties.metadata(metakey);
    
    if (not property.isA<SegProp>())
        throw SireError::invalid_cast( QObject::tr(
            "The metadata at metakey %1 of type %2 is not a SegProperty!")
                .arg(metakey, property.what()), CODELOC );
}

void EditMolData::assertHasSegMetadata(const QString &key,
                                       const QString &metakey) const
{
    const Property &property = properties.metadata(key, metakey);
    
    if (not property.isA<SegProp>())
        throw SireError::invalid_cast( QObject::tr(
            "The metadata at key %1, metakey %2 of type %3 is not a SegProperty!")
                .arg(key, metakey, property.what()), CODELOC );
}

/** Return the values of the atom property at key 'key' 

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
AtomVariantProperty EditMolData::atomProperty(const QString &key) const
{
    //ensure that this atom property exists...
    this->assertHasAtomProperty(key);
    
    int ngroups = cg_by_index.count();
    
    QVector< QVector<QVariant> > values(ngroups);
    
    QVector<QVariant> *values_array = values.data();
    
    for (int i=0; i<ngroups; ++i)
    {
        const EditCGData &cgroup = this->cutGroup( cg_by_index.at(i) );
        
        int nats = cgroup.atoms.count();
        
        if (nats == 0)
            values_array[i] = QVector<QVariant>();
        else
        {
            QVector<QVariant> group_vals(nats);
            QVariant *group_vals_array = group_vals.data();
            
            for (int j=0; j<nats; ++j)
            {
                group_vals_array[j] = this->atom(cgroup.atoms.at(j))
                                             .properties.value(key);
            }
         
            values_array[i] = group_vals;
        }
    }
    
    return AtomVariantProperty( PackedArray2D<QVariant>(values) );
}

/** Return the values of the atom metadata at metakey 'metakey'

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
AtomVariantProperty EditMolData::atomMetadata(const QString &metakey) const
{
    //ensure that this atom metadata exists...
    this->assertHasAtomMetadata(metakey);
    
    int ngroups = cg_by_index.count();
    
    QVector< QVector<QVariant> > values(ngroups);
    
    QVector<QVariant> *values_array = values.data();
    
    for (int i=0; i<ngroups; ++i)
    {
        const EditCGData &cgroup = this->cutGroup( cg_by_index.at(i) );
        
        int nats = cgroup.atoms.count();
        
        if (nats == 0)
            values_array[i] = QVector<QVariant>();
        else
        {
            QVector<QVariant> group_vals(nats);
            QVariant *group_vals_array = group_vals.data();
            
            for (int j=0; j<nats; ++j)
            {
                group_vals_array[j] = this->atom(cgroup.atoms.at(j))
                                             .molecule_metadata.value(metakey);
            }
         
            values_array[i] = group_vals;
        }
    }
    
    return AtomVariantProperty( PackedArray2D<QVariant>(values) );
}

/** Return the values of the atom metadata for the property at
    key 'key', and metadata key 'metakey'

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
AtomVariantProperty EditMolData::atomMetadata(const QString &key,
                                              const QString &metakey) const
{
    //ensure that this atom metadata exists...
    this->assertHasAtomMetadata(key, metakey);
    
    int ngroups = cg_by_index.count();
    
    QVector< QVector<QVariant> > values(ngroups);
    
    QVector<QVariant> *values_array = values.data();
    
    for (int i=0; i<ngroups; ++i)
    {
        const EditCGData &cgroup = this->cutGroup( cg_by_index.at(i) );
        
        int nats = cgroup.atoms.count();
        
        if (nats == 0)
            values_array[i] = QVector<QVariant>();
        else
        {
            QVector<QVariant> group_vals(nats);
            QVariant *group_vals_array = group_vals.data();
            
            for (int j=0; j<nats; ++j)
            {
                group_vals_array[j] = this->atom(cgroup.atoms.at(j))
                                             .property_metadata[key].value(metakey);
            }
         
            values_array[i] = group_vals;
        }
    }
    
    return AtomVariantProperty( PackedArray2D<QVariant>(values) );
}                                                     
            
/** Return the values of the CutGroup property at key 'key'

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
CGVariantProperty EditMolData::cgProperty(const QString &key) const
{
    this->assertHasCGProperty(key);
    
    int ncg = cg_by_index.count();
    
    QVector<QVariant> values(ncg);
    values.squeeze();
    
    QVariant *values_array = values.data();
    
    for (int i=0; i<ncg; ++i)
    {
        values_array[i] = this->cutGroup(cg_by_index.at(i))
                                  .properties.value(key);
    }
    
    return values;
}

/** Return the values of the CutGroup metadata at metakey 'metakey'

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
CGVariantProperty EditMolData::cgMetadata(const QString &metakey) const
{
    this->assertHasCGMetadata(metakey);
    
    int ncg = cg_by_index.count();
    
    QVector<QVariant> values(ncg);
    values.squeeze();
    
    QVariant *values_array = values.data();
    
    for (int i=0; i<ncg; ++i)
    {
        values_array[i] = this->cutGroup(cg_by_index.at(i))
                                  .molecule_metadata.value(metakey);
    }
    
    return values;
}

/** Return the values of the CutGroup metadata at metakey 'metakey'
    for the property at key 'key'

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
CGVariantProperty EditMolData::cgMetadata(const QString &key, 
                                          const QString &metakey) const
{
    this->assertHasCGMetadata(key, metakey);
    
    int ncg = cg_by_index.count();
    
    QVector<QVariant> values(ncg);
    values.squeeze();
    
    QVariant *values_array = values.data();
    
    for (int i=0; i<ncg; ++i)
    {
        values_array[i] = this->cutGroup(cg_by_index.at(i))
                                  .property_metadata[key].value(metakey);
    }
    
    return values;
}
                                          
/** Return the values of the residue property at key 'key'

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
ResVariantProperty EditMolData::resProperty(const QString &key) const
{
    this->assertHasResProperty(key);
    
    int nres = res_by_index.count();
    
    QVector<QVariant> values(nres);
    values.squeeze();
    
    QVariant *values_array = values.data();
    
    for (int i=0; i<nres; ++i)
    {
        values_array[i] = this->residue(res_by_index.at(i))
                                  .properties.value(key);
    }
    
    return values;
}

/** Return the values of the residue metadata at metakey 'metakey'

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
ResVariantProperty EditMolData::resMetadata(const QString &metakey) const
{
    this->assertHasResMetadata(metakey);
    
    int nres = res_by_index.count();
    
    QVector<QVariant> values(nres);
    values.squeeze();
    
    QVariant *values_array = values.data();
    
    for (int i=0; i<nres; ++i)
    {
        values_array[i] = this->residue(res_by_index.at(i))
                                  .molecule_metadata.value(metakey);
    }
    
    return values;
}

/** Return the values of the residue metadata at metakey 'metakey'
    for the property at key 'key'

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
ResVariantProperty EditMolData::resMetadata(const QString &key, 
                                           const QString &metakey) const
{
    this->assertHasResMetadata(key, metakey);
    
    int nres = res_by_index.count();
    
    QVector<QVariant> values(nres);
    values.squeeze();
    
    QVariant *values_array = values.data();
    
    for (int i=0; i<nres; ++i)
    {
        values_array[i] = this->residue(res_by_index.at(i))
                                  .property_metadata[key].value(metakey);
    }
    
    return values;
}
                                          
/** Return the values of the chain property at key 'key'

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
ChainVariantProperty EditMolData::chainProperty(const QString &key) const
{
    this->assertHasChainProperty(key);
    
    int nchains = chains_by_index.count();
    
    QVector<QVariant> values(nchains);
    values.squeeze();
    
    QVariant *values_array = values.data();
    
    for (int i=0; i<nchains; ++i)
    {
        values_array[i] = this->chain(chains_by_index.at(i))
                                  .properties.value(key);
    }
    
    return values;
}

/** Return the values of the chain metadata at metakey 'metakey'

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
ChainVariantProperty EditMolData::chainMetadata(const QString &metakey) const
{
    this->assertHasChainMetadata(metakey);
    
    int nchains = chains_by_index.count();
    
    QVector<QVariant> values(nchains);
    values.squeeze();
    
    QVariant *values_array = values.data();
    
    for (int i=0; i<nchains; ++i)
    {
        values_array[i] = this->chain(chains_by_index.at(i))
                                  .molecule_metadata.value(metakey);
    }
    
    return values;
}

/** Return the values of the chain metadata at metakey 'metakey'
    for the property at key 'key'

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
ChainVariantProperty EditMolData::chainMetadata(const QString &key, 
                                             const QString &metakey) const
{
    this->assertHasChainMetadata(key,metakey);
    
    int nchains = chains_by_index.count();
    
    QVector<QVariant> values(nchains);
    values.squeeze();
    
    QVariant *values_array = values.data();
    
    for (int i=0; i<nchains; ++i)
    {
        values_array[i] = this->chain(chains_by_index.at(i))
                                  .property_metadata[key].value(metakey);
    }
    
    return values;
}
                                          
/** Return the values of the segment property at key 'key'

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
SegVariantProperty EditMolData::segProperty(const QString &key) const
{
    this->assertHasSegProperty(key);
    
    int nseg = seg_by_index.count();
    
    QVector<QVariant> values(nseg);
    values.squeeze();
    
    QVariant *values_array = values.data();
    
    for (int i=0; i<nseg; ++i)
    {
        values_array[i] = this->segment(seg_by_index.at(i))
                                  .properties.value(key);
    }
    
    return values;
}

/** Return the values of the segment metadata at metakey 'metakey'

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
SegVariantProperty EditMolData::segMetadata(const QString &metakey) const
{
    this->assertHasSegMetadata(metakey);
    
    int nseg = seg_by_index.count();
    
    QVector<QVariant> values(nseg);
    values.squeeze();
    
    QVariant *values_array = values.data();
    
    for (int i=0; i<nseg; ++i)
    {
        values_array[i] = this->segment(seg_by_index.at(i))
                                  .molecule_metadata.value(metakey);
    }
    
    return values;
}

/** Return the values of the segment metadata at metakey 'metakey'
    for the property at key 'key'

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
SegVariantProperty EditMolData::segMetadata(const QString &key, 
                                           const QString &metakey) const
{
    this->assertHasSegMetadata(key, metakey);
    
    int nseg = seg_by_index.count();
    
    QVector<QVariant> values(nseg);
    values.squeeze();
    
    QVariant *values_array = values.data();
    
    for (int i=0; i<nseg; ++i)
    {
        values_array[i] = this->segment(seg_by_index.at(i))
                                  .molecule_metadata.value(key, metakey);
    }
    
    return values;
}

/** Return the mapping from the AtomIdx of the atom before editing to 
    the AtomIdx in the current (new version) of the molecule. Note that 
    this will only contain entries for atoms that exist both in the old
    and new versions, and show how to map from the old AtomIdx to the new AtomIdx */
QHash<AtomIdx,AtomIdx> EditMolData::getOldToNewAtomMapping() const
{
    QHash<AtomIdx,AtomIdx> mapping;
    mapping.reserve(atoms_by_index.count());

    for (int i=0; i<atoms_by_index.count(); ++i)
    {
        quint32 uid = atoms_by_index[i];
        
        if (old_atomidxs.contains(uid))
        {
            mapping.insert( old_atomidxs[uid], AtomIdx(i) );
        }
    }
    
    return mapping;
}

/////////
///////// Implementation of EditMolInfo
/////////

/** Null constructor */
EditMolInfo::EditMolInfo() : StructureEditor(), MolInfo()
{}

/** Construct from the passed editor */
EditMolInfo::EditMolInfo(const StructureEditor &editor)
            : StructureEditor(editor), MolInfo()
{}

/** Copy constructor */
EditMolInfo::EditMolInfo(const EditMolInfo &other)
            : StructureEditor(other), MolInfo()
{}

/** Destructor */
EditMolInfo::~EditMolInfo()
{}

/** Assign from the passed editor */
EditMolInfo& EditMolInfo::operator=(const StructureEditor &editor)
{
    StructureEditor::operator=(editor);
    return *this;
}

/** Copy assignment operator */
EditMolInfo& EditMolInfo::operator=(const EditMolInfo &other)
{
    StructureEditor::operator=(other);
    return *this;
}

/** Return a string representation */
QString EditMolInfo::toString() const
{
    return QObject::tr( "EditMolInfo( %1 : %2 )" )
                .arg( this->molName() )
                .arg( this->molNum() );
}

/** Clone this info */
EditMolInfo* EditMolInfo::clone() const
{
    return new EditMolInfo(*this);
}

/** Is this the entire molecule */
bool EditMolInfo::selectedAll() const
{
    return not StructureEditor::isEmpty();
}

/** Return the indicies of all of the atoms that have the name 'name'

    \throw SireMol::missing_atom
*/
QList<AtomIdx> EditMolInfo::map(const AtomName &name) const
{
    int nats = d->atoms_by_index.count();
    
    QList<AtomIdx> atomidxs;
    
    if (name.isCaseSensitive())
    {
        for (AtomIdx i(0); i<nats; ++i)
        {
            if (d->atom(d->atoms_by_index.at(i)).name == name)
                atomidxs.append(i);
        }
    }
    else
    {
        QString lower_name = QString(name).toLower();
    
        for (AtomIdx i(0); i<nats; ++i)
        {
            if ( QString(d->atom(d->atoms_by_index.at(i)).name).toLower() == lower_name )
                atomidxs.append(i);
        }
    }
    
    if (atomidxs.isEmpty())
        throw SireMol::missing_atom( QObject::tr(
            "There is no atom called \"%1\" in the molecule called \"%2\".")
                .arg(name, d->molname), CODELOC );
                
    return atomidxs;
}

/** Return the indicies of all atoms that have the number 'num'

    \throw SireMol::missing_atom
*/
QList<AtomIdx> EditMolInfo::map(AtomNum num) const
{
    int nats = d->atoms_by_index.count();
 
    QList<AtomIdx> atomidxs;
    
    for (AtomIdx i(0); i<nats; ++i)
    {
        if (d->atom(d->atoms_by_index.at(i)).number == num)
            atomidxs.append(i);
    }
    
    if (atomidxs.isEmpty())
        throw SireMol::missing_atom( QObject::tr(
            "There is no atom with number %1 in the molecule called \"%2\".")
                .arg(num).arg(d->molname), CODELOC );
                
    return atomidxs;
}

/** Simple function used to return just the index, if it is valid!

    \throw SireError::invalid_index
*/
QList<AtomIdx> EditMolInfo::map(AtomIdx idx) const
{
    QList<AtomIdx> atomidxs;
    atomidxs.append( AtomIdx(idx.map(d->atoms_by_index.count())) );
    
    return atomidxs;
}

/** Return the indicies of all of the atoms in this molecule
    that match the ID 'atomid'
    
    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
QList<AtomIdx> EditMolInfo::map(const AtomID &atomid) const
{
    return atomid.map(*this);
}

/** Return the indicies of all of the residues in this molecule
    that have the name 'resname'
    
    \throw SireMol::missing_residue
*/
QList<ResIdx> EditMolInfo::map(const ResName &name) const
{
    int nres = d->res_by_index.count();
    
    QList<ResIdx> residxs;
    
    if (name.isCaseSensitive())
    {
        for (ResIdx i(0); i<nres; ++i)
        {
            if (d->residue(d->res_by_index.at(i)).name == name)
                residxs.append(i);
        }
    }
    else
    {
        QString lower_name = QString(name).toLower();
    
        for (ResIdx i(0); i<nres; ++i)
        {
            if ( QString(d->residue(d->res_by_index.at(i)).name).toLower() 
                                                                      == lower_name )
                residxs.append(i);
        }
    }
    
    if (residxs.isEmpty())
        throw SireMol::missing_residue( QObject::tr(
            "There are no residues called \"%1\" in the molecule called \"%2\".")
                .arg(name, d->molname), CODELOC );
                
    return residxs;
}

/** Return the indicies of all of the residues that have number 'num'

    \throw SireMol::missing_residue
*/
QList<ResIdx> EditMolInfo::map(ResNum num) const
{
    int nres = d->res_by_index.count();
    QList<ResIdx> residxs;
    
    for (ResIdx i(0); i<nres; ++i)
    {
        if (d->residue(d->res_by_index.at(i)).number == num)
            residxs.append(i);
    }
    
    if (residxs.isEmpty())
        throw SireMol::missing_residue( QObject::tr(
            "There are no residues with number %1 in the molecule called \"%2\".")
                .arg(num).arg(d->molname), CODELOC );
                
    return residxs;
}

/** Simple overload function

    \throw SireError::invalid_index
*/
QList<ResIdx> EditMolInfo::map(ResIdx idx) const
{
    QList<ResIdx> residxs;
    residxs.append( ResIdx(idx.map(d->res_by_index.count())) );
    
    return residxs;
}

/** Return the indicies of all of the residues that match the ID 'resid'

    \throw SireMol::missing_residue
    \throw SireError::invalid_index
*/
QList<ResIdx> EditMolInfo::map(const ResID &resid) const
{
    return resid.map(*this);
}

/** Return the indicies of all of the CutGroups in this molecule
    that have the name 'name'
    
    \throw SireMol::missing_cutgroup
*/
QList<CGIdx> EditMolInfo::map(const CGName &name) const
{
    int ncg = d->cg_by_index.count();
    
    QList<CGIdx> cgidxs;
    
    if (name.isCaseSensitive())
    {
        for (CGIdx i(0); i<ncg; ++i)
        {
            if (d->cutGroup(d->cg_by_index.at(i)).name == name)
                cgidxs.append(i);
        }
    }
    else
    {
        QString lower_name = QString(name).toLower();
    
        for (CGIdx i(0); i<ncg; ++i)
        {
            if ( QString(d->cutGroup(d->cg_by_index.at(i)).name).toLower() 
                                                                      == lower_name )
                cgidxs.append(i);
        }
    }

    if (cgidxs.isEmpty())
        throw SireMol::missing_cutgroup( QObject::tr(
            "There are no CutGroups called \"%1\" in the molecule called \"%2\".")
                .arg(name, d->molname), CODELOC );
                
    return cgidxs;
}

/** Simple overload function

    \throw SireError::invalid_index
*/
QList<CGIdx> EditMolInfo::map(CGIdx idx) const
{
    QList<CGIdx> cgidxs;
    cgidxs.append( CGIdx(idx.map(d->cg_by_index.count())) );
    
    return cgidxs;
}

/** Return the indicies of all of the CutGroups that match the ID 'cgid'

    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
*/
QList<CGIdx> EditMolInfo::map(const CGID &cgid) const
{
    return cgid.map(*this);
}

/** Return the indicies of all of the Chains in this molecule
    that have the name 'name'
    
    \throw SireMol::missing_chain
*/
QList<ChainIdx> EditMolInfo::map(const ChainName &name) const
{
    int nchains = d->chains_by_index.count();
    
    QList<ChainIdx> chainidxs;
    
    if (name.isCaseSensitive())
    {
        for (ChainIdx i(0); i<nchains; ++i)
        {
            if (d->chain(d->chains_by_index.at(i)).name == name)
                chainidxs.append(i);
        }
    }
    else
    {
        QString lower_name = QString(name).toLower();
    
        for (ChainIdx i(0); i<nchains; ++i)
        {
            if ( QString(d->residue(d->chains_by_index.at(i)).name).toLower() 
                                                                      == lower_name )
                chainidxs.append(i);
        }
    }

    if (chainidxs.isEmpty())
        throw SireMol::missing_chain( QObject::tr(
            "There are no chains called \"%1\" in the molecule called \"%2\".")
                .arg(name, d->molname), CODELOC );
                
    return chainidxs;
}

/** Simple overload function

    \throw SireError::invalid_index
*/
QList<ChainIdx> EditMolInfo::map(ChainIdx idx) const
{
    QList<ChainIdx> chainidxs;
    chainidxs.append( ChainIdx(idx.map(d->chains_by_index.count())) );
    
    return chainidxs;
}

/** Return the indicies of all of the Chains that match the ID 'chainid'

    \throw SireMol::missing_chain
    \throw SireError::invalid_index
*/
QList<ChainIdx> EditMolInfo::map(const ChainID &chainid) const
{
    return chainid.map(*this);
}

/** Return the indicies of all of the Segments in this molecule
    that have the name 'name'
    
    \throw SireMol::missing_segment
*/
QList<SegIdx> EditMolInfo::map(const SegName &name) const
{
    int nsegs = d->seg_by_index.count();
    
    QList<SegIdx> segidxs;
    
    if (name.isCaseSensitive())
    {
        for (SegIdx i(0); i<nsegs; ++i)
        {
            if (d->segment(d->seg_by_index.at(i)).name == name)
                segidxs.append(i);
        }
    }
    else
    {
        QString lower_name = QString(name).toLower();
    
        for (SegIdx i(0); i<nsegs; ++i)
        {
            if ( QString(d->segment(d->seg_by_index.at(i)).name).toLower() 
                                                                      == lower_name )
                segidxs.append(i);
        }
    }

    if (segidxs.isEmpty())
        throw SireMol::missing_segment( QObject::tr(
            "There are no segments called \"%1\" in the molecule called \"%2\".")
                .arg(name, d->molname), CODELOC );
                
    return segidxs;
}

/** Simple overload function

    \throw SireError::invalid_index
*/
QList<SegIdx> EditMolInfo::map(SegIdx idx) const
{
    QList<SegIdx> segidxs;
    segidxs.append( SegIdx(idx.map(d->seg_by_index.count())) );
    
    return segidxs;
}

/** Return the indicies of all of the Segments that match the ID 'segid'

    \throw SireMol::missing_segment
    \throw SireError::invalid_index
*/
QList<SegIdx> EditMolInfo::map(const SegID &segid) const
{
    return segid.map(*this);
}

/** Return the indicies of all of the atoms in this molecule */
QList<AtomIdx> EditMolInfo::getAtoms() const
{
    int nats = d->atoms_by_index.count();

    if (nats == 0)
        throw SireMol::missing_atom( QObject::tr(
            "There are no atoms in this molecule!"), CODELOC );

    QList<AtomIdx> atomidxs;
    
    for (AtomIdx i(0); i<nats; ++i)
    {
        atomidxs.append(i);
    }
    
    return atomidxs;
}

/** Return the index of the ith atom in the CutGroup at index 'cgidx'

    \throw SireError::invalid_index
*/
AtomIdx EditMolInfo::getAtom(CGIdx cgidx, int i) const
{
    const EditCGData &cgroup = d->cutGroup( getUID(cgidx) );
    
    quint32 atomuid = cgroup.atoms.at( Index(i).map(cgroup.atoms.count()) );
    
    return StructureEditor::atomIdx(atomuid);
}

/** Return the index of the ith atom in the residue at index 'residx'

    \throw SireError::invalid_index
*/
AtomIdx EditMolInfo::getAtom(ResIdx residx, int i) const
{
    const EditResData &residue = d->residue( getUID(residx) );
    
    quint32 atomuid = residue.atoms.at( Index(i).map(residue.atoms.count()) );
    
    return StructureEditor::atomIdx(atomuid);
}

/** Return the index of the ith atom in the chain at index 'chainidx'

    \throw SireError::invalid_index
*/
AtomIdx EditMolInfo::getAtom(ChainIdx chainidx, int i) const
{
    QList<AtomIdx> atomidxs = this->getAtomsIn(chainidx);
    
    return atomidxs.at( Index(i).map(atomidxs.count()) );
}

/** Return the index of the ith atom in the segment at index 'segidx'

    \throw SireError::invalid_index
*/
AtomIdx EditMolInfo::getAtom(SegIdx segidx, int i) const
{
    const EditSegData &segment = d->segment( getUID(segidx) );
    
    quint32 atomuid = segment.atoms.at( Index(i).map(segment.atoms.count()) );
    
    return StructureEditor::atomIdx(atomuid);
}

/** Return the indicies of the atoms in the residues identified by 'resid' 

    \throw SireMol::missing_residue
    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
QList<AtomIdx> EditMolInfo::getAtomsIn(const ResID &resid) const
{
    QList<ResIdx> residxs = resid.map(*this);
    
    QList<AtomIdx> atomidxs;
    
    foreach (ResIdx residx, residxs)
    {
        atomidxs += d->atomIdxsFromUIDs( d->residue(getUID(residx)).atoms );
    }
    
    if (atomidxs.isEmpty())
        throw SireMol::missing_atom( QObject::tr(
            "There are no atoms in the residues identified by %1.")
                .arg(resid.toString()), CODELOC );
    
    return atomidxs;
}

/** Return the indicies of the atoms in the CutGroups identified by 'cgid' 

    \throw SireMol::missing_cutgroup
    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
QList<AtomIdx> EditMolInfo::getAtomsIn(const CGID &cgid) const
{
    QList<CGIdx> cgidxs = cgid.map(*this);
    
    QList<AtomIdx> atomidxs;
    
    foreach (CGIdx cgidx, cgidxs)
    {
        atomidxs += d->atomIdxsFromUIDs( d->cutGroup(getUID(cgidx)).atoms );
    }
    
    if (atomidxs.isEmpty())
        throw SireMol::missing_atom( QObject::tr(
            "There are no atoms in the CutGroups identified by %1.")
                .arg(cgid.toString()), CODELOC );
    
    return atomidxs;
}

/** Return the list of atoms in the chains identified by 'chainidx'

    \throw SireMol::missing_chain
    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
QList<AtomIdx> EditMolInfo::getAtomsIn(const ChainID &chainid) const
{
    QList<ChainIdx> chainidxs = chainid.map(*this);
    
    QList<AtomIdx> atomidxs;
    
    foreach (ChainIdx chainidx, chainidxs)
    {
        foreach (quint32 resuid, d->chain(getUID(chainidx)).residues)
        {
            atomidxs += d->atomIdxsFromUIDs( d->residue(resuid).atoms );
        }
    }
    
    if (atomidxs.isEmpty())
        throw SireMol::missing_atom( QObject::tr(
            "There are no atoms in the chains identified by %1.")
                .arg(chainid.toString()), CODELOC );
    
    return atomidxs;
}

/** Return the indicies of all of the atoms in the segment(s) identified
    by 'segid'
    
    \throw SireMol::missing_segment
    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
QList<AtomIdx> EditMolInfo::getAtomsIn(const SegID &segid) const
{
    QList<SegIdx> segidxs = segid.map(*this);
    
    QList<AtomIdx> atomidxs;
    
    foreach (SegIdx segidx, segidxs)
    {
        atomidxs += d->atomIdxsFromUIDs( d->segment(getUID(segidx)).atoms );
    }
    
    if (atomidxs.isEmpty())
        throw SireMol::missing_atom( QObject::tr(
            "There are no atoms in the segments identified by %1.")
                .arg(segid.toString()), CODELOC );
    
    return atomidxs;
}

/** Return the index of the chain that contains the passed residue

    \throw SireError::invalid_index
    \throw SireMol::missing_chain
*/
ChainIdx EditMolInfo::parentChain(ResIdx residx) const
{
    const EditResData &residue = d->residue( getUID(residx) );

    ChainIdx chainidx = d->chainIdx(residue);
    
    if (chainidx.isNull())
        throw SireMol::missing_residue( QObject::tr(
                "The residue %1:%2 is not part of any chain.")
                    .arg(residue.name).arg(residue.number), CODELOC );
                    
    return chainidx;
}

/** Return the index of the chain that contains the passed residue

    \throw SireError::invalid_index
    \throw SireMol::missing_residue
    \throw SireMol::duplicate_residue
    \throw SireMol::missing_chain
*/
ChainIdx EditMolInfo::parentChain(const ResID &resid) const
{
    return this->parentChain( this->resIdx(resid) );
}

/** Return the index of the residue that contains the passed atom

    \throw SireError::invalid_index
    \throw SireMol::missing_residue
*/
ResIdx EditMolInfo::parentResidue(AtomIdx atomidx) const
{
    const EditAtomData &atom = d->atom( getUID(atomidx) );

    ResIdx residx = d->resIdx(atom);
    
    if (residx.isNull())
        throw SireMol::missing_residue( QObject::tr(
                "The atom %1:%2 is not part of any residue.")
                    .arg(atom.name).arg(atom.number), CODELOC );
                    
    return residx;
}

/** Return the index of the chain that contains the passed atom

    \throw SireError::invalid_index
    \throw SireMol::missing_chain
*/
ChainIdx EditMolInfo::parentChain(AtomIdx atomidx) const
{
    return this->parentChain( this->parentResidue(atomidx) );
}

/** Return the index of the CutGroup that contains the passed atom

    \throw SireError::invalid_index
    \throw SireMol::missing_cutgroup
*/
CGIdx EditMolInfo::parentCutGroup(AtomIdx atomidx) const
{
    const EditAtomData &atom = d->atom( getUID(atomidx) );

    CGIdx cgidx = d->cgIdx(atom);
    
    if (cgidx.isNull())
        throw SireMol::missing_cutgroup( QObject::tr(
                "The atom %1:%2 is not part of any CutGroup.")
                    .arg(atom.name).arg(atom.number), CODELOC );
                    
    return cgidx;
}

/** Return the index of the segment that contains the passed atom

    \throw SireError::invalid_index
    \throw SireMol::missing_segment
*/
SegIdx EditMolInfo::parentSegment(AtomIdx atomidx) const
{
    const EditAtomData &atom = d->atom( getUID(atomidx) );

    SegIdx segidx = d->segIdx(atom);
    
    if (segidx.isNull())
        throw SireMol::missing_segment( QObject::tr(
                "The atom %1:%2 is not part of any segment.")
                    .arg(atom.name).arg(atom.number), CODELOC );
                    
    return segidx;
}

/** Return the index of the chain that contains the passed atom

    \throw SireError::invalid_index
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireMol::missing_chain
*/
ChainIdx EditMolInfo::parentChain(const AtomID &atomid) const
{
    return this->parentChain( this->atomIdx(atomid) );
}

/** Return the index of the residue that contains the passed atom

    \throw SireError::invalid_index
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireMol::missing_residue
*/
ResIdx EditMolInfo::parentResidue(const AtomID &atomid) const
{
    return this->parentResidue( this->atomIdx(atomid) );
}

/** Return the index of the CutGroup that contains the passed atom

    \throw SireError::invalid_index
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireMol::missing_cutgroup
*/
CGIdx EditMolInfo::parentCutGroup(const AtomID &atomid) const
{
    return this->parentCutGroup( this->atomIdx(atomid) );
}

/** Return the index of the segment that contains the passed atom

    \throw SireError::invalid_index
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireMol::missing_segment
*/
SegIdx EditMolInfo::parentSegment(const AtomID &atomid) const
{
    return this->parentSegment( this->atomIdx(atomid) );
}

/** Return the indicies of all of the residues in this molecule */
QList<ResIdx> EditMolInfo::getResidues() const
{
    int nres = d->res_by_index.count();
    
    if (nres == 0)
        throw SireMol::missing_residue( QObject::tr(
            "There are no residues in this molecule."), CODELOC );
            
    QList<ResIdx> residxs;
    
    for (ResIdx i(0); i<nres; ++i)
    {
        residxs.append(i);
    }
    
    return residxs;
}

/** Return the ith residue in the chain at index 'chainidx'

    \throw SireError::invalid_index
*/
ResIdx EditMolInfo::getResidue(ChainIdx chainidx, int i) const
{
    const EditChainData &chain = d->chain( getUID(chainidx) );
    
    quint32 resuid = chain.residues.at( Index(i).map(chain.residues.count()) );
    
    return StructureEditor::resIdx(resuid);
}

/** Return the indicies of the residues in the chains identified by 'chainid'

    \throw SireMol::missing_chain
    \throw SireMol::missing_residue
    \throw SireError::invalid_index
*/
QList<ResIdx> EditMolInfo::getResiduesIn(const ChainID &chainid) const
{
    QList<ChainIdx> chainidxs = chainid.map(*this);
    
    QList<ResIdx> residxs;
    
    foreach (ChainIdx chainidx, chainidxs)
    {
        residxs += d->resIdxsFromUIDs( d->chain(getUID(chainidx)).residues );
    }
    
    if (residxs.isEmpty())
        throw SireMol::missing_residue( QObject::tr(
            "There are no residues in the chains identified by %1.")
                .arg(chainid.toString()), CODELOC );
    
    return residxs;
}

/** Return a list of all of the indicies of the CutGroups in this molecule */
QList<CGIdx> EditMolInfo::getCutGroups() const
{
    int ncg = d->cg_by_index.count();
    
    if (ncg == 0)
        throw SireMol::missing_cutgroup( QObject::tr(
            "There are no CutGroups in this molecule."), CODELOC );
            
    QList<CGIdx> cgidxs;
    
    for (CGIdx i(0); i<ncg; ++i)
    {
        cgidxs.append(i);
    }
    
    return cgidxs;
}

/** Return a list of all of the indicies of the chains in this molecule */
QList<ChainIdx> EditMolInfo::getChains() const
{
    int nchains = d->chains_by_index.count();
    
    if (nchains == 0)
        throw SireMol::missing_chain( QObject::tr(
            "There are no chains in this molecule."), CODELOC );
            
    QList<ChainIdx> chainidxs;
    
    for (ChainIdx i(0); i<nchains; ++i)
    {
        chainidxs.append(i);
    }
    
    return chainidxs;
}

/** Return a list of all of the indicies of the segments in this molecule */
QList<SegIdx> EditMolInfo::getSegments() const
{
    int nseg = d->seg_by_index.count();
    
    if (nseg == 0)
        throw SireMol::missing_segment( QObject::tr(
            "There are no segments in this molecule."), CODELOC );
            
    QList<SegIdx> segidxs;
    
    for (SegIdx i(0); i<nseg; ++i)
    {
        segidxs.append(i);
    }
    
    return segidxs;
}

AtomIdx EditMolInfo::atomIdx(const AtomID &atomid) const
{
    return StructureEditor::atomIdx(atomid);
}

CGIdx EditMolInfo::cgIdx(const CGID &cgid) const
{
    return StructureEditor::cgIdx(cgid);
}

ResIdx EditMolInfo::resIdx(const ResID &resid) const
{
    return StructureEditor::resIdx(resid);
}

ChainIdx EditMolInfo::chainIdx(const ChainID &chainid) const
{
    return StructureEditor::chainIdx(chainid);
}

SegIdx EditMolInfo::segIdx(const SegID &segid) const
{
    return StructureEditor::segIdx(segid);
}

void EditMolInfo::assertCompatibleWith(const AtomSelection &selected_atoms) const
{
    throw SireError::incomplete_code(CODELOC);
}

/////////
///////// Implementation of StructureEditor
/////////

static const RegisterMetaType<StructureEditor> r_editor(MAGIC_ONLY,
                                                        "SireMol::StructureEditor");
                                                  
/** Serialise to a binary datastream */
QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds,
                                       const StructureEditor &editor)
{
    writeHeader(ds, r_editor, 1);
    
    SharedDataStream sds(ds);
    
    sds << editor.d;
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds,
                                       StructureEditor &editor)
{
    VersionID v = readHeader(ds, r_editor);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> editor.d;
    }
    else
        throw version_error( v, "1", r_editor, CODELOC );
        
    return ds;
}

/** Null constructor */
StructureEditor::StructureEditor() : d( new EditMolData() )
{}

/** Assign so that this will edit a copy of 'moldata' */
StructureEditor& StructureEditor::operator=(const MoleculeData &moldata)
{
    d.reset( new EditMolData(moldata) );
    
    return *this;
}

/** Construct an editor that edits a copy of 'moldata' */
StructureEditor::StructureEditor(const MoleculeData &moldata)
{
    this->operator=(moldata);
}

/** Copy constructor - this is fast as this class is explicitly
    shared */
StructureEditor::StructureEditor(const StructureEditor &other)
                : d(other.d)
{}

/** Destructor */
StructureEditor::~StructureEditor()
{}

/** Assert that we have a valid shared pointer */
void StructureEditor::assertSane() const
{
    if (d.get() == 0)
        throw SireError::program_bug( QObject::tr(
            "PROGRAM BUG! Null shared_ptr<detail::EditMolData>!"), CODELOC );
}

/** Copy assignment operator - this is fast as this class is 
    explicitly shared */
StructureEditor& StructureEditor::operator=(const StructureEditor &other)
{
    d = other.d;
    return *this;
}

Properties& StructureEditor::_pvt_properties()
{
    return d->properties;
}

CGAtomIdx EditMolData::cgAtomIdx(quint32 atomuid, const EditAtomData &atom) const
{
    if (atom.cg_parent == 0)
        return CGAtomIdx::null();
        
    else
        return CGAtomIdx( CGIdx(cg_by_index.indexOf(atom.cg_parent)), 
                          Index(cutGroup(atom.cg_parent).atoms.indexOf(atomuid) ) );
}

CGIdx EditMolData::cgIdx(const EditAtomData &atom) const
{
    if (atom.cg_parent == 0)
        return CGIdx::null();
    else
        return CGIdx( cg_by_index.indexOf(atom.cg_parent) );
}

ResIdx EditMolData::resIdx(const EditAtomData &atom) const
{
    if (atom.res_parent == 0)
        return ResIdx::null();
    else
        return ResIdx( res_by_index.indexOf(atom.res_parent) );
}

SegIdx EditMolData::segIdx(const EditAtomData &atom) const
{
    if (atom.seg_parent == 0)
        return SegIdx::null();
    else
        return SegIdx( seg_by_index.indexOf(atom.seg_parent) );
}

/** Return all of the metadata about the atom at index 'atomidx'

    \throw SireError::invalid_index
*/
tuple<AtomName,AtomNum,CGAtomIdx,ResIdx,SegIdx> 
StructureEditor::getAtomData(AtomIdx atomidx) const
{
    this->assertSane();

    quint32 atomuid = getUID(atomidx);
    const EditAtomData &atom = d->atom(atomuid);
    
    return tuple<AtomName,AtomNum,CGAtomIdx,ResIdx,SegIdx>(
                 atom.name, atom.number, d->cgAtomIdx(atomuid, atom),
                 d->resIdx(atom), d->segIdx(atom) );
}

/** Return all of the metadata about the CutGroup at index 'cgidx' 

    \throw SireError::invalid_index
*/
boost::tuple< CGName,QList<AtomIdx> >
StructureEditor::getCGData(CGIdx cgidx) const
{
    this->assertSane();

    const EditCGData &cgroup = d->cutGroup( getUID(cgidx) );
    
    if (cgroup.atoms.isEmpty())
        return tuple< CGName,QList<AtomIdx> >(cgroup.name, QList<AtomIdx>());
    else
    {
        QList<AtomIdx> atomidxs;
        
        foreach (quint32 atom, cgroup.atoms)
        {
            atomidxs.append( AtomIdx(d->atoms_by_index.indexOf(atom)) );
        }
        
        return tuple< CGName,QList<AtomIdx> >(cgroup.name, atomidxs);
    }
}

ChainIdx EditMolData::chainIdx(const EditResData &residue) const
{
    if (residue.chain_parent == 0)
        return ChainIdx::null();
    else
        return ChainIdx( chains_by_index.indexOf(residue.chain_parent) );
}

/** Return the metadata for the residue at index 'residx' 

    \throw SireError::invalid_index
*/
boost::tuple< ResName,ResNum,ChainIdx,QList<AtomIdx> >
StructureEditor::getResData(ResIdx residx) const
{
    this->assertSane();

    const EditResData &residue = d->residue( getUID(residx) );
    
    QList<AtomIdx> atomidxs;
    
    foreach (quint32 atom, residue.atoms)
    {
        atomidxs.append( AtomIdx(d->atoms_by_index.indexOf(atom) ) );
    }
    
    return tuple< ResName,ResNum,ChainIdx,QList<AtomIdx> >(
                    residue.name, residue.number, d->chainIdx(residue), atomidxs );
}

/** Return the metadata for the chain at index 'chainidx'

    \throw SireError::invalid_index
*/
boost::tuple< ChainName,QList<ResIdx> >
StructureEditor::getChainData(ChainIdx chainidx) const
{
    this->assertSane();

    const EditChainData &chain = d->chain( getUID(chainidx) );
    
    QList<ResIdx> residxs;
    
    foreach (quint32 residue, chain.residues)
    {
        residxs.append( ResIdx(d->res_by_index.indexOf(residue)) );
    }
    
    return tuple< ChainName,QList<ResIdx> >(chain.name, residxs);
}

/** Return the metadata for the segment at index 'segidx' 

    \throw SireError::invalid_index
*/
boost::tuple< SegName,QList<AtomIdx> >
StructureEditor::getSegData(SegIdx segidx) const
{
    this->assertSane();

    const EditSegData &segment = d->segment( getUID(segidx) );
    
    QList<AtomIdx> atomidxs;
    
    foreach (quint32 atom, segment.atoms)
    {
        atomidxs.append( AtomIdx(d->atoms_by_index.indexOf(atom)) );
    }
    
    return tuple< SegName,QList<AtomIdx> >(segment.name, atomidxs);
}

/** Return whether or not the MoleculeInfoData info molecule layout object
    needs to be rebuilt */
bool StructureEditor::needsInfoRebuild() const
{
    this->assertSane();

    return d->cached_molinfo.data() == 0;
}

/** Return the cached MoleculeInfoData object

    \throw SireError::program_bug
*/
const MoleculeInfoData& StructureEditor::info() const
{
    this->assertSane();

    if (d->cached_molinfo.constData() == 0)
        throw SireError::program_bug( QObject::tr(
            "We must never access the info object when it is invalid!"),
                CODELOC );
                
    return *(d->cached_molinfo);
}

/** Commit the MoleculeInfo object - this creates the molecule 
    layout from the current system */
const MoleculeInfoData& StructureEditor::commitInfo()
{
    this->assertSane();

    if (this->needsInfoRebuild())
        d->cached_molinfo = new MoleculeInfoData(*this);
        
    return *(d->cached_molinfo);
}

/** Commit the changes - this creates a MoleculeData object that contains
    all of the data in this editor */
MoleculeData StructureEditor::commitChanges() const
{
    this->assertSane();
    return MoleculeData(*this);
}

/** Return the number of atoms in the residue identified by 'uid'

    \throw SireMol::missing_residue
*/
int StructureEditor::nAtomsInResidue(quint32 uid) const
{
    this->assertSane();

    return d->residue(uid).atoms.count();
}

/** Return the number of atoms in the CutGroup identified by 'uid'

    \throw SireMol::missing_cutgroup
*/
int StructureEditor::nAtomsInCutGroup(quint32 uid) const
{
    this->assertSane();

    return d->cutGroup(uid).atoms.count();
}

/** Return the number of residues in the chain identified by 'uid'

    \throw SireMol::missing_chain
*/
int StructureEditor::nResiduesInChain(quint32 uid) const
{
    this->assertSane();

    return d->chain(uid).residues.count();
}

/** Return the number of atoms in the chain identified by 'uid'

    \throw SireMol::missing_chain
*/
int StructureEditor::nAtomsInChain(quint32 uid) const
{
    this->assertSane();

    int nats = 0;
    
    foreach (quint32 resuid, d->chain(uid).residues)
    {
        nats += d->residue(resuid).atoms.count();
    }
    
    return nats;
}

/** Return the number of atoms in the segment identified by 'uid'

    \throw SireMol::missing_segment
*/
int StructureEditor::nAtomsInSegment(quint32 uid) const
{
    this->assertSane();

    return d->segment(uid).atoms.count();
}

/** Return the number of atoms this molecule */
int StructureEditor::nAtomsInMolecule() const
{
    this->assertSane();

    return d->atoms.count();
}

/** Return whether or not this molecule is empty */
bool StructureEditor::isEmpty() const
{
    return this->nAtomsInMolecule() == 0;
}

/** Return the number of residues in this molecule */
int StructureEditor::nResiduesInMolecule() const
{
    this->assertSane();

    return d->residues.count();
}

/** Return the number of chains in this molecule */
int StructureEditor::nChainsInMolecule() const
{
    this->assertSane();

    return d->chains.count();
}

/** Return the number of segments in this molecule */
int StructureEditor::nSegmentsInMolecule() const
{
    this->assertSane();

    return d->segments.count();
}

/** Return the number of CutGroups in this molecule */
int StructureEditor::nCutGroupsInMolecule() const
{
    this->assertSane();

    return d->cutgroups.count();
}

/** Assert that the UID 'uid' refers to a valid atom 

    \throw SireMol::missing_atom
*/
void StructureEditor::assertValidAtom(quint32 uid) const
{
    this->assertSane();

    d->atom(uid);
}

/** Assert that the UID 'uid' refers to a valid CutGroup

    \throw SireMol::missing_cutgroup
*/
void StructureEditor::assertValidCutGroup(quint32 uid) const
{
    this->assertSane();

    d->cutGroup(uid);
}

/** Assert that the UID 'uid' refers to a valid residue

    \throw SireMol::missing_residue
*/
void StructureEditor::assertValidResidue(quint32 uid) const
{
    this->assertSane();

    d->residue(uid);
}

/** Assert that the UID 'uid' refers to a valid chain

    \throw SireMol::missing_chain
*/
void StructureEditor::assertValidChain(quint32 uid) const
{
    this->assertSane();

    d->chain(uid);
}

/** Assert that the UID 'uid' refers to a valid segment

    \throw SireMol::missing_segment
*/
void StructureEditor::assertValidSegment(quint32 uid) const
{
    this->assertSane();

    d->segment(uid);
}

/** Return the UID of the atom at index 'atomidx'

    \throw SireError::invalid_index
*/
quint32 StructureEditor::getUID(AtomIdx atomidx) const
{
    this->assertSane();

    return d->atoms_by_index.at( atomidx.map(nAtomsInMolecule()) );
}

/** Return the UID of the CutGroup at index 'cgidx'

    \throw SireError::invalid_index
*/
quint32 StructureEditor::getUID(CGIdx cgidx) const
{
    this->assertSane();

    return d->cg_by_index.at( cgidx.map(nCutGroupsInMolecule()) );
}

/** Return the UID of the residue at index 'residx'

    \throw SireError::invalid_index
*/
quint32 StructureEditor::getUID(ResIdx residx) const
{
    this->assertSane();

    return d->res_by_index.at( residx.map(nResiduesInMolecule()) );
}

/** Return the UID of the chain at index 'chainidx'

    \throw SireError::invalid_index
*/
quint32 StructureEditor::getUID(ChainIdx chainidx) const
{
    this->assertSane();

    return d->chains_by_index.at( chainidx.map(nChainsInMolecule()) );
}

/** Return the UID of the segment at index 'segidx'

    \throw SireError::invalid_index
*/
quint32 StructureEditor::getUID(SegIdx segidx) const
{
    this->assertSane();

    return d->seg_by_index.at( segidx.map(nSegmentsInMolecule()) );
}

/** Return the UID of the atom that matches the ID 'atomid'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
quint32 StructureEditor::getUID(const AtomID &atomid) const
{
    this->assertSane();

    return this->getUID( atomIdx(atomid) );
}

/** Return the UID of the CutGroup that matches the ID 'cgid'

    \throw SireMol::missing_cutgroup
    \throw SireMol::duplicate_cutgroup
    \throw SireError::invalid_index
*/
quint32 StructureEditor::getUID(const CGID &cgid) const
{
    this->assertSane();

    return this->getUID( cgIdx(cgid) );
}

/** Return the UID of the residue that matches the ID 'resid'

    \throw SireMol::missing_residue
    \throw SireMol::duplicate_residue
    \throw SireError::invalid_index
*/
quint32 StructureEditor::getUID(const ResID &resid) const
{
    this->assertSane();

    return this->getUID( resIdx(resid) );
}

/** Return the UID of the chain that matches the ID 'chainid'

    \throw SireMol::missing_chain
    \throw SireMol::duplicate_chain
    \throw SireError::invalid_index
*/
quint32 StructureEditor::getUID(const ChainID &chainid) const
{
    this->assertSane();

    return this->getUID( chainIdx(chainid) );
}

/** Return the UID of the segment that matches the ID 'segid'

    \throw SireMol::missing_segment
    \throw SireMol::duplicate_segment
    \throw SireError::invalid_index
*/
quint32 StructureEditor::getUID(const SegID &segid) const
{
    this->assertSane();

    return this->getUID( segIdx(segid) );
}

/** Return the UID of the ith atom in the CutGroup identified by 'uid'

    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
*/
quint32 StructureEditor::atomInCutGroup(quint32 uid, int i) const
{
    this->assertSane();

    const EditCGData &cgroup = d->cutGroup(uid);
    
    return cgroup.atoms.at( Index(i).map(cgroup.atoms.count()) );
}

/** Return the UID of the ith atom in the residue identified by 'uid'

    \throw SireMol::missing_residue
    \throw SireError::invalid_index
*/
quint32 StructureEditor::atomInResidue(quint32 uid, int i) const
{
    this->assertSane();

    const EditResData &residue = d->residue(uid);
    
    return residue.atoms.at( Index(i).map(residue.atoms.count()) );
}

/** Return the UID of the ith atom in the segment identified
    by 'uid'
    
    \throw SireMol::missing_segment
    \throw SireError::invalid_index
*/
quint32 StructureEditor::atomInSegment(quint32 uid, int i) const
{
    this->assertSane();

    const EditSegData &segment = d->segment(uid);
    
    return segment.atoms.at( Index(i).map(segment.atoms.count()) );
}

/** Return the UID of the ith residue in the chain identified
    by 'uid'
    
    \throw SireMol::missing_chain
    \throw SireError::invalid_index
*/
quint32 StructureEditor::residueInChain(quint32 uid, int i) const
{
    this->assertSane();

    const EditChainData &chain = d->chain(uid);
    
    return chain.residues.at( Index(i).map(chain.residues.count()) );
}

/** Return the UID of the CutGroup that contains the atom identified
    by 'uid'
    
    \throw SireMol::missing_atom
    \throw SireMol::missing_cutgroup
*/
quint32 StructureEditor::cutGroupParentOfAtom(quint32 uid) const
{
    this->assertSane();

    const EditAtomData &atom = d->atom(uid);
    
    if (atom.cg_parent == 0)
        throw SireMol::missing_cutgroup( QObject::tr(
            "The atom at index %1 is not part of a CutGroup.")
                .arg(this->atomIdx(uid)), CODELOC );
                
    return atom.cg_parent;
}

/** Return the UID of the residue that contains the atom identified
    by 'uid'
    
    \throw SireMol::missing_atom
    \throw SireMol::missing_residue
*/
quint32 StructureEditor::residueParentOfAtom(quint32 uid) const
{
    this->assertSane();

    const EditAtomData &atom = d->atom(uid);
    
    if (atom.res_parent == 0)
        throw SireMol::missing_residue( QObject::tr(
            "The atom at index %1 is not part of a residue.")
                .arg(this->atomIdx(uid)), CODELOC );
                
    return atom.res_parent;
}

/** Return the UID of the chain that contains the residue
    identified by the ID 'uid'
    
    \throw SireMol::missing_residue
    \throw SireMol::missing_chain
*/
quint32 StructureEditor::chainParentOfResidue(quint32 uid) const
{
    this->assertSane();

    const EditResData &residue = d->residue(uid);
    
    if (residue.chain_parent == 0)
        throw SireMol::missing_chain( QObject::tr(
            "The residue at index %1 is not part of a chain.")
                .arg(this->resIdx(uid)), CODELOC );
                
    return residue.chain_parent;
}

/** Return the UID of the chain that contains the atom identified
    by the UID 'uid'
    
    \throw SireMol::missing_atom
    \throw SireMol::missing_residue
    \throw SireMol::missing_chain
*/
quint32 StructureEditor::chainParentOfAtom(quint32 uid) const
{
    return this->chainParentOfResidue( this->residueParentOfAtom(uid) );
}

/** Return the UID of the segment that contains the atom identified
    by the UID 'uid'
    
    \throw SireMol::missing_atom
    \throw SireMol::missing_segment
*/
quint32 StructureEditor::segmentParentOfAtom(quint32 uid) const
{
    this->assertSane();

    const EditAtomData &atom = d->atom(uid);
    
    if (atom.seg_parent == 0)
        throw SireMol::missing_segment( QObject::tr(
            "The atom at index %1 is not part of a segment.")
                .arg(this->atomIdx(uid)), CODELOC );
                
    return atom.seg_parent;
}

/** Return the name of this molecule */
const MolName& StructureEditor::molName() const
{
    this->assertSane();

    return d->molname;
}

/** Return the number of this molecule */
MolNum StructureEditor::molNum() const
{
    this->assertSane();

    return d->molnum;
}

/** Return the name of the atom identified by 'uid'

    \throw SireMol::missing_atom
*/
const AtomName& StructureEditor::atomName(quint32 uid) const
{
    this->assertSane();

    return d->atom(uid).name;
}

/** Return the number of the atom identified by 'uid' 
    
    \throw SireMol::missing_atom
*/
AtomNum StructureEditor::atomNum(quint32 uid) const
{
    this->assertSane();

    return d->atom(uid).number;
}

/** Return the index of the atom identified by 'uid'

    \throw SireMol::missing_atom
*/
AtomIdx StructureEditor::atomIdx(quint32 uid) const
{
    this->assertSane();

    int i = d->atoms_by_index.indexOf(uid);
    
    if (i == -1)
        throw SireMol::missing_atom( QObject::tr(
            "There is no atom identified by the UID %1 in this molecule.")
                .arg(uid), CODELOC );
                
    return AtomIdx(i);
}

/** Return the name of the CutGroup identified by 'uid'

    \throw SireMol::missing_cutgroup
*/
const CGName& StructureEditor::cgName(quint32 uid) const
{
    this->assertSane();

    return d->cutGroup(uid).name;
}

/** Return the index of the CutGroup identified by 'uid'

    \throw SireMol::missing_cutgroup
*/
CGIdx StructureEditor::cgIdx(quint32 uid) const
{
    this->assertSane();

    int i = d->cg_by_index.indexOf(uid);
    
    if (i == -1)
        throw SireMol::missing_cutgroup( QObject::tr(
            "There is no CutGroup identified by the UID %1 in this molecule.")
                .arg(uid), CODELOC );
                
    return CGIdx(i);
}

/** Return the name of the residue identified by 'uid'

    \throw SireMol::missing_residue
*/
const ResName& StructureEditor::resName(quint32 uid) const
{
    this->assertSane();

    return d->residue(uid).name;
}

/** Return the number of the residue identified by 'uid'

    \throw SireMol::missing_residue
*/
ResNum StructureEditor::resNum(quint32 uid) const
{
    this->assertSane();

    return d->residue(uid).number;
}

/** Return the index of the residue identified by 'uid'

    \throw SireMol::missing_residue
*/
ResIdx StructureEditor::resIdx(quint32 uid) const
{
    this->assertSane();

    int i = d->res_by_index.indexOf(uid);
    
    if (i == -1)
        throw SireMol::missing_residue( QObject::tr(
            "There is no residue identified by the UID %1 in this molecule.")
                .arg(uid), CODELOC );
                
    return ResIdx(i);
}

/** Return the name of the chain identified by 'uid'

    \throw SireMol::missing_chain
*/
const ChainName& StructureEditor::chainName(quint32 uid) const
{
    this->assertSane();

    return d->chain(uid).name;
}

/** Return the index of the chain identified by 'uid'

    \throw SireMol::missing_chain
*/
ChainIdx StructureEditor::chainIdx(quint32 uid) const
{
    this->assertSane();

    int i = d->chains_by_index.indexOf(uid);
    
    if (i == -1)
        throw SireMol::missing_chain( QObject::tr(
            "There is no chain identified by the UID %1 in this molecule.")
                .arg(uid), CODELOC );
                
    return ChainIdx(i);
}

/** Return the name of the segment identified by 'uid'

    \throw SireMol::missing_segment
*/
const SegName& StructureEditor::segName(quint32 uid) const
{
    this->assertSane();

    return d->segment(uid).name;
}

/** Return the index of the segment identified by 'uid'

    \throw SireMol::missing_segment
*/
SegIdx StructureEditor::segIdx(quint32 uid) const
{
    this->assertSane();

    int i = d->seg_by_index.indexOf(uid);
    
    if (i == -1)
        throw SireMol::missing_segment( QObject::tr(
            "There is no segment identified by the UID %1 in this molecule.")
                .arg(uid), CODELOC );
                
    return SegIdx(i);
}

/** Return the index of the atom that matches the ID 'atomid'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
AtomIdx StructureEditor::atomIdx(const AtomID &atomid) const
{
    this->assertSane();

    QList<AtomIdx> atomidxs = atomid.map( EditMolInfo(*this) );
    
    if (atomidxs.count() > 1)
        throw SireMol::duplicate_atom( QObject::tr(
            "More than one atom in the molecule matches the ID \"%1\" %2")
                .arg(atomid.toString(), Sire::toString(atomidxs)), CODELOC );
                
    return atomidxs.first();
}

/** Return the index of the CutGroup that matches the ID 'cgid'

    \throw SireMol::missing_cutgroup
    \throw SireMol::duplicate_cutgroup
    \throw SireError::invalid_index
*/
CGIdx StructureEditor::cgIdx(const CGID &cgid) const
{
    this->assertSane();

    QList<CGIdx> cgidxs = cgid.map( EditMolInfo(*this) );
    
    if (cgidxs.count() > 1)
        throw SireMol::duplicate_cutgroup( QObject::tr(
            "More than one CutGroup in the molecule matches the ID \"%1\" %2")
                .arg(cgid.toString(), Sire::toString(cgidxs)), CODELOC );
                
    return cgidxs.first();
}

/** Return the index of the residue that matches the ID 'resid'

    \throw SireMol::missing_residue
    \throw SireMol::duplicate_residue
    \throw SireError::invalid_index
*/
ResIdx StructureEditor::resIdx(const ResID &resid) const
{
    this->assertSane();

    QList<ResIdx> residxs = resid.map( EditMolInfo(*this) );
    
    if (residxs.count() > 1)
        throw SireMol::duplicate_residue( QObject::tr(
            "More than one residue in the molecule matches the ID \"%1\" %2")
                .arg(resid.toString(), Sire::toString(residxs)), CODELOC );
                
    return residxs.first();
}

/** Return the index of the chain that matches the ID 'chainid'

    \throw SireMol::missing_chain
    \throw SireMol::duplicate_chain
    \throw SireError::invalid_index
*/
ChainIdx StructureEditor::chainIdx(const ChainID &chainid) const
{
    this->assertSane();

    QList<ChainIdx> chainidxs = chainid.map( EditMolInfo(*this) );
    
    if (chainidxs.count() > 1)
        throw SireMol::duplicate_chain( QObject::tr(
            "More than one chain in the molecule matches the ID \"%1\" %2")
                .arg(chainid.toString(), Sire::toString(chainidxs)), CODELOC );
                
    return chainidxs.first();
}

/** Return the index of the segment that matches the ID 'segid'

    \throw SireMol::missing_segment
    \throw SireMol::duplicate_segment
    \throw SireError::invalid_index
*/
SegIdx StructureEditor::segIdx(const SegID &segid) const
{
    this->assertSane();

    QList<SegIdx> segidxs = segid.map( EditMolInfo(*this) );
    
    if (segidxs.count() > 1)
        throw SireMol::duplicate_segment( QObject::tr(
            "More than one segment in the molecule matches the ID \"%1\" %2")
                .arg(segid.toString(), Sire::toString(segidxs)), CODELOC );
                
    return segidxs.first();
}

/** Rename this molecule to 'newname' */
void StructureEditor::renameMolecule(const MolName &newname)
{
    this->assertSane();

    d->molname = MolName( cacheName(newname) );
}

/** Give this molecule a new, unique, number */
void StructureEditor::renumberMolecule()
{
    this->assertSane();

    d->molnum = MolNum::getUniqueNumber();
}

/** Renumber this molecule to 'newnum' */
void StructureEditor::renumberMolecule(MolNum newnum)
{
    this->assertSane();

    d->molnum = newnum;
}

/** Rename the atom identified by 'uid' to 'newname'

    \throw SireMol::missing_atom
*/
void StructureEditor::renameAtom(quint32 uid, const AtomName &newname)
{
    this->assertSane();

    EditAtomData &atom = d->atom(uid);
    
    if (atom.name != newname)
    {
        atom.name = AtomName( cacheName(newname) );
        d->cached_molinfo = 0;
    }
}

/** Renumber the atom identified by 'uid' to 'newnum'

    \throw SireMol::missing_atom
*/
void StructureEditor::renumberAtom(quint32 uid, AtomNum newnum)
{
    this->assertSane();

    EditAtomData &atom = d->atom(uid);
    
    if (atom.number != newnum)
    {
        atom.number = newnum;
        d->cached_molinfo = 0;
    }
}

static void changeIndex(QList<quint32> &uids, quint32 uid, int newidx)
{
    if (uids.removeAll(uid) == 0)
        throw SireError::program_bug( QObject::tr(
            "Couldn't remove UID %1 from UIDs %2")
                .arg(uid).arg( Sire::toString(uids) ), CODELOC );
                
    if (newidx < 0)
    {
        newidx = uids.count() + newidx;
        
        if (newidx < 0)
            newidx = 0;
    }
    else if (newidx > uids.count())
    {
        newidx = uids.count();
    }
    
    uids.insert( newidx, uid );
}

/** Change the index of the atom identified by 'uid' to 'newidx'

    \throw SireMol::missing_atom
*/
void StructureEditor::reindexAtom(quint32 uid, AtomIdx newidx)
{
    this->assertSane();

    AtomIdx old_idx = this->atomIdx(uid);

    changeIndex( d->atoms_by_index, uid, newidx );

    if (this->atomIdx(uid) != old_idx)
        d->cached_molinfo = 0;
}

/** Rename the CutGroup identified by 'uid' to 'newname'

    \throw SireMol::missing_cutgroup
*/
void StructureEditor::renameCutGroup(quint32 uid, const CGName &newname)
{
    this->assertSane();

    if ( this->cgName(uid) != newname )
    {
        d->cutGroup(uid).name = CGName( cacheName(newname) );
        d->cached_molinfo = 0;
    }
}

/** Change the index of the atom identified by 'uid' to 'newidx'

    \throw SireMol::missing_cutgroup
*/
void StructureEditor::reindexCutGroup(quint32 uid, CGIdx newidx)
{
    this->assertSane();

    CGIdx old_idx = this->cgIdx(uid);

    changeIndex( d->cg_by_index, uid, newidx );

    if (this->cgIdx(uid) != old_idx)
        d->cached_molinfo = 0;
}

/** Rename the residue identified by 'uid' to 'newname'

    \throw SireMol::missing_residue
*/
void StructureEditor::renameResidue(quint32 uid, const ResName &newname)
{
    this->assertSane();

    if (this->resName(uid) != newname)
    {
        d->residue(uid).name = ResName( cacheName(newname) );
        d->cached_molinfo = 0;
    }
}

/** Renumber the residue identified by 'uid' to 'newnum'

    \throw SireMol::missing_residue
*/
void StructureEditor::renumberResidue(quint32 uid, ResNum newnum)
{
    this->assertSane();

    if (this->resNum(uid) != newnum)
    {
        d->residue(uid).number = newnum;
        d->cached_molinfo = 0;
    }
}

/** Change the index of the residue identified by 'uid' to 'newidx'

    \throw SireMol::missing_residue
*/
void StructureEditor::reindexResidue(quint32 uid, ResIdx newidx)
{
    this->assertSane();

    ResIdx old_idx = this->resIdx(uid);
    
    changeIndex( d->res_by_index, uid, newidx );

    if (this->resIdx(uid) != old_idx)
        d->cached_molinfo = 0;
}

/** Rename the chain identified by 'uid' to 'newname'

    \throw SireMol::missing_chain
*/
void StructureEditor::renameChain(quint32 uid, const ChainName &newname)
{
    this->assertSane();

    if (this->chainName(uid) != newname)
    {
        d->chain(uid).name = ChainName( cacheName(newname) );
        d->cached_molinfo = 0;
    }
}

/** Change the index of the chain identified by 'uid' to 'newidx'

    \throw SireMol::missing_chain
*/
void StructureEditor::reindexChain(quint32 uid, ChainIdx newidx)
{
    this->assertSane();

    ChainIdx old_idx = this->chainIdx(uid);
    
    changeIndex( d->chains_by_index, uid, newidx );

    if (this->chainIdx(uid) != old_idx)
        d->cached_molinfo = 0;
}

/** Rename the segment identified by 'uid' to 'newname'

    \throw SireMol::missing_segment
*/
void StructureEditor::renameSegment(quint32 uid, const SegName &newname)
{
    this->assertSane();

    if (this->segName(uid) != newname)
    {
        d->segment(uid).name = SegName( cacheName(newname) );
        d->cached_molinfo = 0;
    }
}

/** Change the index of the segment identified by 'uid' to 'newidx' */
void StructureEditor::reindexSegment(quint32 uid, SegIdx newidx)
{
    this->assertSane();

    SegIdx old_idx = this->segIdx(uid);
    
    changeIndex( d->seg_by_index, uid, newidx );
    
    if (this->segIdx(uid) != old_idx)
        d->cached_molinfo = 0;
}

/** Remove the atom identified by 'uid'

    \throw SireMol::missing_atom
*/
void StructureEditor::removeAtom(quint32 uid)
{
    this->assertSane();

    const EditAtomData &atom = d->atom(uid);
    
    //remove the atom from its parent groups...
    if (atom.cg_parent != 0)
        d->cutGroup(atom.cg_parent).atoms.removeAll(uid);
        
    if (atom.res_parent != 0)
        d->residue(atom.res_parent).atoms.removeAll(uid);
        
    if (atom.seg_parent != 0)
        d->segment(atom.seg_parent).atoms.removeAll(uid);
    
    //now remove the atom itself
    d->atoms.remove(uid);
    d->atoms_by_index.removeAll(uid);
    
    d->cached_molinfo = 0;
}

/** Remove the CutGroup identified by 'uid'. This only
    removes the CutGroup - it doesn't remove the atoms
    that are contained in this CutGroup
    
    \throw SireMol::missing_cutgroup
*/
void StructureEditor::removeCutGroup(quint32 uid)
{
    this->assertSane();

    const EditCGData &cgroup = d->cutGroup(uid);
    
    //tell the atoms that they don't now have a CG parent
    foreach (quint32 atom, cgroup.atoms)
    {
        d->atom(atom).cg_parent = 0;
    }
    
    d->cutgroups.remove(uid);
    d->cg_by_index.removeAll(uid);

    d->cached_molinfo = 0;
}

/** Remove the residue identified by 'uid'

    \throw SireMol:missing_residue
*/
void StructureEditor::removeResidue(quint32 uid)
{
    this->assertSane();

    const EditResData &residue = d->residue(uid);
    
    //tell the parent Chain that this residue is leaving...
    if (residue.chain_parent != 0)
        d->chain(residue.chain_parent).residues.removeAll(uid);
        
    //now tell the atoms that this residue is leaving
    foreach (quint32 atom, residue.atoms)
    {
        d->atom(atom).res_parent = 0;
    }
    
    d->residues.remove(uid);
    d->res_by_index.removeAll(uid);

    d->cached_molinfo = 0;
}

/** Remove the chain identified by 'uid'

    \throw SireMol::missing_chain
*/
void StructureEditor::removeChain(quint32 uid)
{
    this->assertSane();

    const EditChainData &chain = d->chain(uid);
    
    //tell each residue that this chain is leaving
    foreach (quint32 residue, chain.residues)
    {
        d->residue(residue).chain_parent = 0;
    }
    
    d->chains.remove(uid);
    d->chains_by_index.removeAll(uid);

    d->cached_molinfo = 0;
}

/** Remove the segment identified by 'uid'

    \throw SireMol::missing_segment
*/
void StructureEditor::removeSegment(quint32 uid)
{
    this->assertSane();

    const EditSegData &segment = d->segment(uid);
    
    //tell each atom that this segment is leaving
    foreach (quint32 atom, segment.atoms)
    {
        d->atom(atom).seg_parent = 0;
    }
    
    d->segments.remove(uid);
    d->seg_by_index.removeAll(uid);

    d->cached_molinfo = 0;
}

/** Remove all atoms identified by 'atomid'

    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
void StructureEditor::removeAtoms(const AtomID &atomid)
{
    this->assertSane();

    QList<AtomIdx> atomidxs = atomid.map( EditMolInfo(*this) );
    
    if (atomidxs.count() == 1)
    {
        this->removeAtom( getUID(atomidxs.first()) );
        return;
    }
    
    //convert the atomidx to uid
    QList<quint32> uids;
    
    foreach (AtomIdx atomidx, atomidxs)
    {
        uids.append( getUID(atomidx) );
    }
    
    //now remove all of the atoms
    foreach (quint32 uid, uids)
    {
        this->removeAtom(uid);
    }

    d->cached_molinfo = 0;
}

/** Remove all CutGroups identified by 'cgid'

    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
*/
void StructureEditor::removeCutGroups(const CGID &cgid)
{
    this->assertSane();

    QList<CGIdx> cgidxs = cgid.map( EditMolInfo(*this) );
    
    if (cgidxs.count() == 1)
    {
        this->removeCutGroup( getUID(cgidxs.first()) );
        return;
    }
    
    //convert the cgidx to uid
    QList<quint32> uids;
    
    foreach (CGIdx cgidx, cgidxs)
    {
        uids.append( getUID(cgidx) );
    }
    
    //now remove all of the CutGroups
    foreach (quint32 uid, uids)
    {
        this->removeCutGroup(uid);
    }

    d->cached_molinfo = 0;
}

/** Remove all residues identified by 'resid'

    \throw SireMol::missing_residue
    \throw SireError::invalid_index
*/
void StructureEditor::removeResidues(const ResID &resid)
{
    this->assertSane();

    QList<ResIdx> residxs = resid.map( EditMolInfo(*this) );
    
    if (residxs.count() == 1)
    {
        this->removeResidue( getUID(residxs.first()) );
        return;
    }
    
    //convert the residx to uid
    QList<quint32> uids;
    
    foreach (ResIdx residx, residxs)
    {
        uids.append( getUID(residx) );
    }
    
    //now remove all of the residues
    foreach (quint32 uid, uids)
    {
        this->removeResidue(uid);
    }

    d->cached_molinfo = 0;
}

/** Remove all chains identified by 'chainid'

    \throw SireMol::missing_chain
    \throw SireError::invalid_index
*/
void StructureEditor::removeChains(const ChainID &chainid)
{
    this->assertSane();

    QList<ChainIdx> chainidxs = chainid.map( EditMolInfo(*this) );
    
    if (chainidxs.count() == 1)
    {
        this->removeChain( getUID(chainidxs.first()) );
        return;
    }
    
    //convert the chainidx to uid
    QList<quint32> uids;
    
    foreach (ChainIdx chainidx, chainidxs)
    {
        uids.append( getUID(chainidx) );
    }
    
    //now remove all of the chains
    foreach (quint32 uid, uids)
    {
        this->removeChain(uid);
    }

    d->cached_molinfo = 0;
}

/** Remove all segments identified by 'segid'

    \throw SireMol::missing_segment
    \throw SireError::invalid_index
*/
void StructureEditor::removeSegments(const SegID &segid)
{
    this->assertSane();

    QList<SegIdx> segidxs = segid.map( EditMolInfo(*this) );
    
    if (segidxs.count() == 1)
    {
        this->removeSegment( getUID(segidxs.first()) );
        return;
    }
    
    //convert the segidx to uid
    QList<quint32> uids;
    
    foreach (SegIdx segidx, segidxs)
    {
        uids.append( getUID(segidx) );
    }
    
    //now remove all of the segments
    foreach (quint32 uid, uids)
    {
        this->removeSegment(uid);
    }

    d->cached_molinfo = 0;
}

/** Remove all atoms from this molecule */
void StructureEditor::removeAllAtoms()
{
    this->assertSane();

    QList<quint32> atoms_by_index = d->atoms_by_index;
    
    if (not atoms_by_index.isEmpty())
    {
        foreach (quint32 atom, atoms_by_index)
        {
            this->removeAtom(atom);
        }
        
        d->cached_molinfo = 0;
    }
}

/** Remove all CutGroups from this molecule */
void StructureEditor::removeAllCutGroups()
{
    this->assertSane();

    QList<quint32> cg_by_index = d->cg_by_index;
    
    if (not cg_by_index.isEmpty())
    {
        foreach (quint32 cg, cg_by_index)
        {
            this->removeCutGroup(cg);
        }
        
        d->cached_molinfo = 0;
    }
}

/** Remove all residues from this molecule */
void StructureEditor::removeAllResidues()
{
    this->assertSane();

    QList<quint32> res_by_index = d->res_by_index;
    
    if (not res_by_index.isEmpty())
    {
        foreach (quint32 res, res_by_index)
        {
            this->removeResidue(res);
        }
        
        d->cached_molinfo = 0;
    }
}

/** Remove all chains from this molecule */
void StructureEditor::removeAllChains()
{
    this->assertSane();

    QList<quint32> chains_by_index = d->chains_by_index;
    
    if (not chains_by_index.isEmpty())
    {
        foreach (quint32 chain, chains_by_index)
        {
            this->removeChain(chain);
        }
        
        d->cached_molinfo = 0;
    }
}

/** Remove all segments from this molecule */
void StructureEditor::removeAllSegments()
{
    this->assertSane();

    QList<quint32> seg_by_index = d->seg_by_index;
    
    if (not seg_by_index.isEmpty())
    {
        foreach (quint32 seg, seg_by_index)
        {
            this->removeSegment(seg);
        }
        
        d->cached_molinfo = 0;
    }
}

/** Move the atom identified by 'uid' into the CutGroup at index 'cgidx'

    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
void StructureEditor::reparentAtom(quint32 uid, CGIdx cgidx)
{
    this->assertSane();

    EditAtomData &atom = d->atom(uid);
    
    quint32 cg_uid = this->getUID(cgidx);
    
    if (atom.cg_parent == cg_uid)
        return;
        
    //remove the atom from the old CutGroup
    if (atom.cg_parent != 0)    
        d->cutGroup(atom.cg_parent).atoms.removeAll(uid);
        
    //add the atom to the new CutGroup
    if (cg_uid != 0)
        d->cutGroup(cg_uid).atoms.append(uid);
        
    atom.cg_parent = cg_uid;

    d->cached_molinfo = 0;
}

/** Move the atom identified by 'uid' into the residue at index 'residx'

    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
void StructureEditor::reparentAtom(quint32 uid, ResIdx residx)
{
    this->assertSane();

    EditAtomData &atom = d->atom(uid);
    
    quint32 res_uid = this->getUID(residx);
    
    if (atom.res_parent == res_uid)
        return;
        
    //remove the atom from the old residue
    if (atom.res_parent != 0)    
        d->residue(atom.res_parent).atoms.removeAll(uid);
        
    //add the atom to the new residue
    if (res_uid != 0)
        d->residue(res_uid).atoms.append(uid);
        
    atom.res_parent = res_uid;

    d->cached_molinfo = 0;
}

/** Move the atom identified by 'uid' into the segment at index 'segidx'

    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
void StructureEditor::reparentAtom(quint32 uid, SegIdx segidx)
{
    this->assertSane();

    EditAtomData &atom = d->atom(uid);
    
    quint32 seg_uid = this->getUID(segidx);
    
    if (atom.seg_parent == seg_uid)
        return;
        
    //remove the atom from the old segment
    if (atom.seg_parent != 0)    
        d->segment(atom.seg_parent).atoms.removeAll(uid);
        
    //add the atom to the new segment
    if (seg_uid != 0)
        d->segment(seg_uid).atoms.append(uid);
        
    atom.seg_parent = seg_uid;

    d->cached_molinfo = 0;
}

/** Move the residue identified by 'uid' into the chain at index 'chainidx'

    \throw SireMol::missing_residue
    \throw SireError::invalid_index
*/
void StructureEditor::reparentResidue(quint32 uid, ChainIdx chainidx)
{
    this->assertSane();

    EditResData &residue = d->residue(uid);
    
    quint32 chain_uid = this->getUID(chainidx);
    
    if (residue.chain_parent == chain_uid)
        return;
        
    //remove the residue from the old chain
    if (residue.chain_parent != 0)    
        d->chain(residue.chain_parent).residues.removeAll(uid);
        
    //add the residue to the new chain
    if (chain_uid != 0)
        d->chain(chain_uid).residues.append(uid);
        
    residue.chain_parent = chain_uid;

    d->cached_molinfo = 0;
}

/** Add a new atom to the molecule, returning an editor for that atom */
AtomStructureEditor StructureEditor::addAtom()
{
    this->assertSane();

    quint32 uid = d->getNewUID();
    d->atoms.insert( uid, EditAtomData() );
    d->atoms_by_index.append(uid);
    d->cached_molinfo = 0;
    
    return AtomStructureEditor( *this, AtomIdx(nAtomsInMolecule()-1) );
}

/** Add a new CutGroup to the molecule, returning an editor for that CutGroup */
CGStructureEditor StructureEditor::addCutGroup()
{
    this->assertSane();

    quint32 uid = d->getNewUID();
    d->cutgroups.insert( uid, EditCGData() );
    d->cg_by_index.append(uid);
    d->cached_molinfo = 0;
    
    return CGStructureEditor( *this, CGIdx(nCutGroupsInMolecule()-1) );
}

/** Add a new residue to the molecule, returning an editor for that residue */
ResStructureEditor StructureEditor::addResidue()
{
    this->assertSane();

    quint32 uid = d->getNewUID();
    d->residues.insert( uid, EditResData() );
    d->res_by_index.append(uid);
    d->cached_molinfo = 0;
    
    return ResStructureEditor( *this, ResIdx(nResiduesInMolecule()-1) );
}

/** Add a new chain to the molecule, returning an editor for that chain */
ChainStructureEditor StructureEditor::addChain()
{
    this->assertSane();

    quint32 uid = d->getNewUID();
    d->chains.insert( uid, EditChainData() );
    d->chains_by_index.append(uid);
    d->cached_molinfo = 0;
    
    return ChainStructureEditor( *this, ChainIdx(nChainsInMolecule()-1) );
}

/** Add a new segment to the molecule, returning an editor for that segment */
SegStructureEditor StructureEditor::addSegment()
{
    this->assertSane();

    quint32 uid = d->getNewUID();
    d->segments.insert( uid, EditSegData() );
    d->seg_by_index.append(uid);
    d->cached_molinfo = 0;
    
    return SegStructureEditor( *this, SegIdx(nSegmentsInMolecule()-1) );
}

/** Return an AtomSelection based on the bool AtomProperty 'values' */
AtomSelection StructureEditor::extractAtomSelection(
                                    const AtomVariantProperty &values) const
{
    this->assertSane();

    //create an AtomSelection using the current MoleculeInfo object
    if (this->needsInfoRebuild())
        const_cast<StructureEditor*>(this)->commitInfo();

    AtomSelection selected_atoms( this->info() );
    selected_atoms.selectNone();
    
    int ncg = values.count();
    const PackedArray2D<QVariant>::Array *values_array = values.array().constData();
    
    for (CGIdx i(0); i<ncg; ++i)
    {
        const PackedArray2D<QVariant>::Array &group_values = values_array[i];
        int nats = group_values.count();
        const QVariant *group_values_array = group_values.constData();
        
        for (Index j(0); j<nats; ++j)
        {
            const QVariant &value = group_values_array[j];
            
            if (not value.canConvert<bool>())
                throw SireError::invalid_cast( QObject::tr(
                    "Cannot convert the object of type %1 to a bool, "
                    "as required for extraction into an AtomSelection!")
                        .arg(value.typeName()), CODELOC );
                        
            if (value.value<bool>())
            {
                selected_atoms.select( CGAtomIdx(i,j) );
            }
        }
    }
    
    return selected_atoms;
}

/** This is a custom AtomMatcher that is used exclusively by
    StructureEditor to migrate properties between different
    versions of the same molecule
    
    @author Christopher Woods
*/
class StructureEditorAtomMatcher
         : public SireBase::ConcreteProperty<StructureEditorAtomMatcher,AtomMatcher>
{

public:
    StructureEditorAtomMatcher()
          : SireBase::ConcreteProperty<StructureEditorAtomMatcher,AtomMatcher>()
    {}
    
    StructureEditorAtomMatcher(const QHash<AtomIdx,AtomIdx> &map)
          : SireBase::ConcreteProperty<StructureEditorAtomMatcher,AtomMatcher>(),
            mapping(map)
    {}
    
    StructureEditorAtomMatcher(const StructureEditorAtomMatcher &other)
          : SireBase::ConcreteProperty<StructureEditorAtomMatcher,AtomMatcher>(other),
            mapping(other.mapping)
    {}
    
    ~StructureEditorAtomMatcher()
    {}
    
    static const char* typeName()
    {
        return "SireMol::StructureEditorAtomMatcher";
    }
    
    const char* what() const
    {
        return StructureEditorAtomMatcher::typeName();
    }

    StructureEditorAtomMatcher& operator=(const StructureEditorAtomMatcher &other)
    {
        mapping = other.mapping;
        return *this;
    }
    
    bool operator==(const StructureEditorAtomMatcher &other) const
    {
        return mapping == other.mapping;
    }
    
    bool operator!=(const StructureEditorAtomMatcher &other) const
    {
        return mapping != other.mapping;
    }
    
protected:
    bool pvt_changesAtomOrder(const MoleculeInfoData &molinfo0,
                              const MoleculeInfoData &molinfo1) const
    {
        return true;
    }
    
    bool pvt_changesAtomOrder(const MoleculeView &mol0,
                              const PropertyMap &map0,
                              const MoleculeView &mol1,
                              const PropertyMap &map1) const
    {
        return true;
    }
    
    QHash<AtomIdx,AtomIdx> pvt_match(const MoleculeInfoData &molinfo0,
                                     const MoleculeInfoData &molinfo1) const
    {
        return mapping;
    }

    QHash<AtomIdx,AtomIdx> pvt_match(const MoleculeView &mol0,
                                     const PropertyMap &map0,
                                     const MoleculeView &mol1,
                                     const PropertyMap &map1) const
    {
        return mapping;
    }

private:
    QHash<AtomIdx,AtomIdx> mapping;
};

/** Return the properties of this molecule - this is a slow function
    as it has to convert all of the properties from an editable to
    a fixed format. Also note that changing the molecule in any
    way will mean that the properties will have to be recalculated */
Properties StructureEditor::properties() const
{
    this->assertSane();

    //make sure that the MoleculeInfo object is up to date
    if (this->needsInfoRebuild())
        const_cast<StructureEditor*>(this)->commitInfo();

    //now get the mapping from the original AtomIdx indicies used when
    //the old properties were constructed to the new AtomIdx indicies of the
    //atoms in the edited molecule
    StructureEditorAtomMatcher mapping( d->getOldToNewAtomMapping() );

    //go through each property in turn and extract it based
    //on its current type
    Properties updated_properties;
    
    for (Properties::const_iterator it = d->properties.constBegin();
         it != d->properties.constEnd();
         ++it)
    {
        PropertyPtr updated_property = it.value();
        const QString &key = it.key();
    
        if (updated_property->isA<AtomProp>())
        {
            updated_property.edit()
                .asA<AtomProp>().assignFrom( d->atomProperty(key) );
        }
        else if (updated_property->isA<CGProp>())
        {
            updated_property.edit()
                .asA<CGProp>().assignFrom( d->cgProperty(key) );
        }
        else if (updated_property->isA<ResProp>())
        {
            updated_property.edit()
                .asA<ResProp>().assignFrom( d->resProperty(key) );
        }
        else if (updated_property->isA<ChainProp>())
        {
            updated_property.edit()
                .asA<ChainProp>().assignFrom( d->chainProperty(key) );
        }
        else if (updated_property->isA<SegProp>())
        {
            updated_property.edit()
                .asA<SegProp>().assignFrom( d->segProperty(key) );
        }
        else if (updated_property->isA<AtomSelection>())
        {
            updated_property = this->extractAtomSelection(
                                            d->atomProperty(key) );
        }
        else if (updated_property->isA<MolViewProperty>())
        {
            try
            {
                updated_property = updated_property->asA<MolViewProperty>()
                                            .makeCompatibleWith(this->info(), mapping);
            }
            catch(...)
            {
                updated_property = NullProperty();
            }
        }
        
        if (not updated_property.isNull())
            updated_properties.setProperty(key, updated_property, false);
        
        //now update all of the metadata for this particular property
        const Properties &metadata = d->properties.allMetadata(it.key());
        
        for (Properties::const_iterator it2 = metadata.constBegin();
             it2 != metadata.constEnd();
             ++it2)
        {
            updated_property = it2.value();
            const QString &metakey = it2.key();
            
            if (updated_property->isA<AtomProp>())
            {
                updated_property.edit()
                    .asA<AtomProp>().assignFrom( d->atomMetadata(key, metakey) );
            }
            else if (updated_property->isA<CGProp>())
            {
                updated_property.edit()
                    .asA<CGProp>().assignFrom( d->cgMetadata(key, metakey) );
            }
            else if (updated_property->isA<ResProp>())
            {
                updated_property.edit()
                    .asA<ResProp>().assignFrom( d->resMetadata(key, metakey) );
            }
            else if (updated_property->isA<ChainProp>())
            {
                updated_property.edit()
                    .asA<ChainProp>().assignFrom( d->chainMetadata(key, metakey) );
            }
            else if (updated_property->isA<SegProp>())
            {
                updated_property.edit()
                    .asA<SegProp>().assignFrom( d->segMetadata(key, metakey) );
            }
            else if (updated_property->isA<AtomSelection>())
            {
                updated_property = this->extractAtomSelection(
                                                d->atomMetadata(key, metakey) );
            }
            else if (updated_property->isA<MolViewProperty>())
            {
                try
                {
                    updated_property = updated_property->asA<MolViewProperty>()
                                                .makeCompatibleWith(this->info(), mapping);
                }
                catch(...)
                {
                    qDebug() << QObject::tr("WARNING: Could not convert the old metadata "
                                            "at key %1:%3 with type %2 to match the new, edited "
                                            "molecule. This metadata has been deleted.")
                            .arg(key, updated_property->what(), metakey);
                
                    updated_property = 0;
                }
            }
            else
            {
                qDebug() << QObject::tr("WARNING: The old metadata at key %1.%3 with type %2 "
                                        "is not a molecule property (MolViewProperty) so has "
                                        "not been edited by the change in molecule. This may mean "
                                        "that this metadata is incompatible with the new molecule.")
                                .arg(key, updated_property->what(), metakey);
            }
            
            if (not updated_property.isNull())
                updated_properties.setMetadata(key, metakey, updated_property);
        }
    }
    
    //the last step is converting all of the molecule metadata
    const Properties &metadata = d->properties.allMetadata();
    
    for (Properties::const_iterator it = metadata.constBegin();
         it != metadata.constEnd();
         ++it)
    {
        PropertyPtr updated_property = it.value();
        const QString &metakey = it.key();
    
        if (updated_property->isA<AtomProp>())
        {
            updated_property.edit()
                .asA<AtomProp>().assignFrom( d->atomMetadata(metakey) );
        }
        else if (updated_property->isA<CGProp>())
        {
            updated_property.edit()
                .asA<CGProp>().assignFrom( d->cgMetadata(metakey) );
        }
        else if (updated_property->isA<ResProp>())
        {
            updated_property.edit()
                .asA<ResProp>().assignFrom( d->resMetadata(metakey) );
        }
        else if (updated_property->isA<ChainProp>())
        {
            updated_property.edit()
                .asA<ChainProp>().assignFrom( d->chainMetadata(metakey) );
        }
        else if (updated_property->isA<SegProp>())
        {
            updated_property.edit()
                .asA<SegProp>().assignFrom( d->segMetadata(metakey) );
        }
        else if (updated_property->isA<AtomSelection>())
        {
            updated_property = this->extractAtomSelection(
                                            d->atomMetadata(metakey) );
        }
        else if (updated_property->isA<MolViewProperty>())
        {
            try
            {
                updated_property = updated_property->asA<MolViewProperty>()
                                            .makeCompatibleWith(this->info(), mapping);
            }
            catch(...)
            {
                qDebug() << QObject::tr("WARNING: Could not convert the old metadata "
                                        "at key %1 with type %2 to match the new, edited "
                                        "molecule. This metadata has been deleted.")
                        .arg(metakey, updated_property->what());
            
                updated_property = 0;
            }
        }
        else
        {
            qDebug() << QObject::tr("WARNING: The old metadata at key %1 with type %2 "
                                    "is not a molecule property (MolViewProperty) so has "
                                    "not been edited by the change in molecule. This may mean "
                                    "that this metadata is incompatible with the new molecule.")
                            .arg(metakey, updated_property->what());
        }
        
        if (not updated_property.isNull())
            updated_properties.setMetadata(metakey, updated_property);
    }
    
    return updated_properties;
}

/** Return the value of the property at key 'key' for the atom
    identified by the UID 'uid'
    
    \throw SireMol::missing_atom
    \throw SireBase::missing_property
*/
const QVariant& StructureEditor::getAtomProperty(quint32 uid, 
                                                 const QString &key) const
{
    this->assertSane();

    const EditAtomData &atom = d->atom(uid);

    QHash<QString,QVariant>::const_iterator it = atom.properties.constFind(key);
                                
    if (it == atom.properties.constEnd())
        throw SireBase::missing_property( QObject::tr(
            "There is no AtomProperty identified by the key '%1'. "
            "Available AtomProperties are [%2].")
                .arg(key, QStringList(
                            atom.properties.keys()).join(", ")), CODELOC );
                
    return it.value();
}

/** Return the value of the property at key 'key' for the residue
    identified by the UID 'uid'
    
    \throw SireMol::missing_residue
    \throw SireBase::missing_property
*/
const QVariant& StructureEditor::getResProperty(quint32 uid, 
                                                const QString &key) const
{
    this->assertSane();

    const EditResData &residue = d->residue(uid);

    QHash<QString,QVariant>::const_iterator it = residue.properties.constFind(key);
                                
    if (it == residue.properties.constEnd())
        throw SireBase::missing_property( QObject::tr(
            "There is no ResProperty identified by the key '%1'. "
            "Available ResProperties are [%2].")
                .arg(key, QStringList(
                            residue.properties.keys()).join(", ")), CODELOC );
                
    return it.value();
}

/** Return the value of the property at key 'key' for the CutGroup
    identified by the UID 'uid'
    
    \throw SireMol::missing_cutgroup
    \throw SireBase::missing_property
*/
const QVariant& StructureEditor::getCGProperty(quint32 uid,
                                               const QString &key) const
{
    this->assertSane();

    const EditCGData &cgroup = d->cutGroup(uid);

    QHash<QString,QVariant>::const_iterator it = cgroup.properties.constFind(key);
                                
    if (it == cgroup.properties.constEnd())
        throw SireBase::missing_property( QObject::tr(
            "There is no CGProperty identified by the key '%1'. "
            "Available CGProperties are [%2].")
                .arg(key, QStringList(
                            cgroup.properties.keys()).join(", ")), CODELOC );
                
    return it.value();
}

/** Return the value of the property at key 'key' for the chain
    identified by the UID 'uid'
    
    \throw SireMol::missing_chain
    \throw SireBase::missing_property
*/
const QVariant& StructureEditor::getChainProperty(quint32 uid, 
                                                  const QString &key) const
{
    this->assertSane();

    const EditChainData &chain = d->chain(uid);

    QHash<QString,QVariant>::const_iterator it = chain.properties.constFind(key);
                                
    if (it == chain.properties.constEnd())
        throw SireBase::missing_property( QObject::tr(
            "There is no ChainProperty identified by the key '%1'. "
            "Available ChainProperties are [%2].")
                .arg(key, QStringList(
                            chain.properties.keys()).join(", ")), CODELOC );
                
    return it.value();
}

/** Return the value of the property at key 'key' for the segment
    identified by the UID 'uid'
    
    \throw SireMol::missing_segment
    \throw SireBase::missing_property
*/
const QVariant& StructureEditor::getSegProperty(quint32 uid, 
                                                const QString &key) const
{
    this->assertSane();

    const EditSegData &segment = d->segment(uid);

    QHash<QString,QVariant>::const_iterator it = segment.properties.constFind(key);
                                
    if (it == segment.properties.constEnd())
        throw SireBase::missing_property( QObject::tr(
            "There is no SegProperty identified by the key '%1'. "
            "Available SegProperties are [%2].")
                .arg(key, QStringList(
                            segment.properties.keys()).join(", ")), CODELOC );
                
    return it.value();
}

/** Return the value of the metadata at metakey 'metakey' for the atom
    identified by the UID 'uid'
    
    \throw SireMol::missing_atom
    \throw SireBase::missing_property
*/
const QVariant& StructureEditor::getAtomMetadata(quint32 uid, 
                                                 const QString &metakey) const
{
    this->assertSane();

    const EditAtomData &atom = d->atom(uid);

    QHash<QString,QVariant>::const_iterator 
                        it = atom.molecule_metadata.constFind(metakey);
                                
    if (it == atom.molecule_metadata.constEnd())
        throw SireBase::missing_property( QObject::tr(
            "There is no AtomProperty metadata identified by the metakey '%1'. "
            "Available AtomProperties are [%2].")
                .arg(metakey, QStringList(
                            atom.molecule_metadata.keys()).join(", ")), CODELOC );
                
    return it.value();
}

/** Return the value of the metadata at metakey 'metakey' for the residue
    identified by the UID 'uid'
    
    \throw SireMol::missing_residue
    \throw SireBase::missing_property
*/
const QVariant& StructureEditor::getResMetadata(quint32 uid, 
                                                const QString &metakey) const
{
    this->assertSane();

    const EditResData &residue = d->residue(uid);

    QHash<QString,QVariant>::const_iterator 
                           it = residue.molecule_metadata.constFind(metakey);
                                
    if (it == residue.molecule_metadata.constEnd())
        throw SireBase::missing_property( QObject::tr(
            "There is no ResProperty metadata identified by the metakey '%1'. "
            "Available ResProperties are [%2].")
                .arg(metakey, QStringList(
                            residue.molecule_metadata.keys()).join(", ")), CODELOC );
                
    return it.value();
}

/** Return the value of the metadata at metakey 'metakey' for the CutGroup
    identified by the UID 'uid'
    
    \throw SireMol::missing_cutgroup
    \throw SireBase::missing_property
*/
const QVariant& StructureEditor::getCGMetadata(quint32 uid,
                                               const QString &metakey) const
{
    this->assertSane();

    const EditCGData &cgroup = d->cutGroup(uid);

    QHash<QString,QVariant>::const_iterator 
                            it = cgroup.molecule_metadata.constFind(metakey);
                                
    if (it == cgroup.molecule_metadata.constEnd())
        throw SireBase::missing_property( QObject::tr(
            "There is no CGProperty metadata identified by the metakey '%1'. "
            "Available CGProperties are [%2].")
                .arg(metakey, QStringList(
                            cgroup.molecule_metadata.keys()).join(", ")), CODELOC );
                
    return it.value();
}

/** Return the value of the metadata at metakey 'metakey' for the chain
    identified by the UID 'uid'
    
    \throw SireMol::missing_chain
    \throw SireBase::missing_property
*/
const QVariant& StructureEditor::getChainMetadata(quint32 uid, 
                                                  const QString &metakey) const
{
    this->assertSane();

    const EditChainData &chain = d->chain(uid);

    QHash<QString,QVariant>::const_iterator 
                        it = chain.molecule_metadata.constFind(metakey);
                                
    if (it == chain.molecule_metadata.constEnd())
        throw SireBase::missing_property( QObject::tr(
            "There is no ChainProperty metadata identified by the metakey '%1'. "
            "Available ChainProperties are [%2].")
                .arg(metakey, QStringList(
                            chain.molecule_metadata.keys()).join(", ")), CODELOC );
                
    return it.value();
}

/** Return the value of the metadata at metakey 'metakey' for the segment
    identified by the UID 'uid'
    
    \throw SireMol::missing_segment
    \throw SireBase::missing_property
*/
const QVariant& StructureEditor::getSegMetadata(quint32 uid, 
                                                const QString &metakey) const
{
    this->assertSane();

    const EditSegData &segment = d->segment(uid);

    QHash<QString,QVariant>::const_iterator 
                        it = segment.molecule_metadata.constFind(metakey);
                                
    if (it == segment.molecule_metadata.constEnd())
        throw SireBase::missing_property( QObject::tr(
            "There is no SegProperty metadata identified by the metakey '%1'. "
            "Available SegProperties are [%2].")
                .arg(metakey, QStringList(
                            segment.molecule_metadata.keys()).join(", ")), CODELOC );
                
    return it.value();
}

/** Return the value of the metadata at metakey 'metakey' for the 
    property at key 'key' for the atom identified by the UID 'uid'
    
    \throw SireMol::missing_atom
    \throw SireBase::missing_property
*/
const QVariant& StructureEditor::getAtomMetadata(quint32 uid, const QString &key,
                                                 const QString &metakey) const
{
    this->assertSane();

    const EditAtomData &atom = d->atom(uid);

    QHash< QString,QHash<QString,QVariant> >::const_iterator
                        it2 = atom.property_metadata.constFind(key);
                        
    if (it2 == atom.property_metadata.constEnd())
        throw SireBase::missing_property( QObject::tr(
            "There is no metadata for the property at key '%1'. Available "
            "atom metadata for this property is [%2].")
                .arg(key, QStringList(
                            atom.property_metadata.keys()).join(", ")), CODELOC );

    QHash<QString,QVariant>::const_iterator 
                        it = it2->constFind(metakey);
                                
    if (it == it2->constEnd())
        throw SireBase::missing_property( QObject::tr(
            "There is no AtomProperty metadata for the property '%1' "
            "identified by the metakey '%2'. "
            "Available AtomProperties are [%3].")
                .arg(key, metakey, 
                     QStringList(it2->keys()).join(", ")), CODELOC );
                
    return it.value();
}

/** Return the value of the metadata at metakey 'metakey' for the 
    property at key 'key' for the residue identified by the UID 'uid'
    
    \throw SireMol::missing_residue
    \throw SireBase::missing_property
*/
const QVariant& StructureEditor::getResMetadata(quint32 uid, const QString &key,
                                                const QString &metakey) const
{
    this->assertSane();

    const EditResData &residue = d->residue(uid);

    QHash< QString,QHash<QString,QVariant> >::const_iterator
                        it2 = residue.property_metadata.constFind(key);
                        
    if (it2 == residue.property_metadata.constEnd())
        throw SireBase::missing_property( QObject::tr(
            "There is no metadata for the property at key '%1'. Available "
            "residue metadata for this property is [%2].")
                .arg(key, QStringList(
                            residue.property_metadata.keys()).join(", ")), CODELOC );

    QHash<QString,QVariant>::const_iterator 
                        it = it2->constFind(metakey);
                                
    if (it == it2->constEnd())
        throw SireBase::missing_property( QObject::tr(
            "There is no ResProperty metadata for the property '%1' "
            "identified by the metakey '%2'. "
            "Available ResProperties are [%3].")
                .arg(key, metakey, 
                     QStringList(it2->keys()).join(", ")), CODELOC );
                
    return it.value();
}

/** Return the value of the metadata at metakey 'metakey' for the 
    property at key 'key' for the CutGroup identified by the UID 'uid'
    
    \throw SireMol::missing_cutgroup
    \throw SireBase::missing_property
*/
const QVariant& StructureEditor::getCGMetadata(quint32 uid, const QString &key,
                                               const QString &metakey) const
{
    this->assertSane();

    const EditCGData &cgroup = d->cutGroup(uid);

    QHash< QString,QHash<QString,QVariant> >::const_iterator
                        it2 = cgroup.property_metadata.constFind(key);
                        
    if (it2 == cgroup.property_metadata.constEnd())
        throw SireBase::missing_property( QObject::tr(
            "There is no metadata for the property at key '%1'. Available "
            "CutGroup metadata for this property is [%2].")
                .arg(key, QStringList(
                            cgroup.property_metadata.keys()).join(", ")), CODELOC );

    QHash<QString,QVariant>::const_iterator 
                        it = it2->constFind(metakey);
                                
    if (it == it2->constEnd())
        throw SireBase::missing_property( QObject::tr(
            "There is no CGProperty metadata for the property '%1' "
            "identified by the metakey '%2'. "
            "Available CGProperties are [%3].")
                .arg(key, metakey, 
                     QStringList(it2->keys()).join(", ")), CODELOC );
                
    return it.value();
}

/** Return the value of the metadata at metakey 'metakey' for the 
    property at key 'key' for the chain identified by the UID 'uid'
    
    \throw SireMol::missing_chain
    \throw SireBase::missing_property
*/
const QVariant& StructureEditor::getChainMetadata(quint32 uid, const QString &key,
                                                  const QString &metakey) const
{
    this->assertSane();

    const EditChainData &chain = d->chain(uid);

    QHash< QString,QHash<QString,QVariant> >::const_iterator
                        it2 = chain.property_metadata.constFind(key);
                        
    if (it2 == chain.property_metadata.constEnd())
        throw SireBase::missing_property( QObject::tr(
            "There is no metadata for the property at key '%1'. Available "
            "chain metadata for this property is [%2].")
                .arg(key, QStringList(
                            chain.property_metadata.keys()).join(", ")), CODELOC );

    QHash<QString,QVariant>::const_iterator 
                        it = it2->constFind(metakey);
                                
    if (it == it2->constEnd())
        throw SireBase::missing_property( QObject::tr(
            "There is no ChainProperty metadata for the property '%1' "
            "identified by the metakey '%2'. "
            "Available ChainProperties are [%3].")
                .arg(key, metakey, 
                     QStringList(it2->keys()).join(", ")), CODELOC );
                
    return it.value();
}

/** Return the value of the metadata at metakey 'metakey' for the 
    property at key 'key' for the segment identified by the UID 'uid'
    
    \throw SireMol::missing_segment
    \throw SireBase::missing_property
*/
const QVariant& StructureEditor::getSegMetadata(quint32 uid, const QString &key,
                                                const QString &metakey) const
{
    this->assertSane();

    const EditSegData &segment = d->segment(uid);

    QHash< QString,QHash<QString,QVariant> >::const_iterator
                        it2 = segment.property_metadata.constFind(key);
                        
    if (it2 == segment.property_metadata.constEnd())
        throw SireBase::missing_property( QObject::tr(
            "There is no metadata for the property at key '%1'. Available "
            "segment metadata for this property is [%2].")
                .arg(key, QStringList(
                            segment.property_metadata.keys()).join(", ")), CODELOC );

    QHash<QString,QVariant>::const_iterator 
                        it = it2->constFind(metakey);
                                
    if (it == it2->constEnd())
        throw SireBase::missing_property( QObject::tr(
            "There is no SegProperty metadata for the property '%1' "
            "identified by the metakey '%2'. "
            "Available SegProperties are [%3].")
                .arg(key, metakey, 
                     QStringList(it2->keys()).join(", ")), CODELOC );
                
    return it.value();
}

/** Throw an invalid cast exception - used to move this out of template
    functions */
void StructureEditor::_pvt_invalidPropertyCast(const QString &key,
                                               const QString &existing_type,
                                               const QString &new_type)
{
    throw SireError::invalid_cast( QObject::tr(
        "Cannot cast the property at key %1 into a %2, as it "
        "is currently a %3.")
            .arg(key, new_type, existing_type), CODELOC );
}

/** Protected function used to set the property at key 'key' 
    of the atom identified by 'uid' to the value 'value'. You
    *must* have already checked that this a valid thing to do,
    as this function does *no* error checking! */
void StructureEditor::_pvt_setAtomProperty(quint32 uid, const QString &key, 
                                           const QVariant &value)
{
    this->assertSane();

    d->atom(uid).properties.insert(key, value);
}

/** Protected function used to set the metadata at metakey 'metakey' 
    of the atom identified by 'uid' to the value 'value'. You
    *must* have already checked that this a valid thing to do,
    as this function does *no* error checking! */
void StructureEditor::_pvt_setAtomMetadata(quint32 uid, const QString &metakey,
                                           const QVariant &value)
{
    this->assertSane();

    d->atom(uid).molecule_metadata.insert(metakey, value);
}
                          
/** Protected function used to set the metadata at metakey
    'metakey' of the property at key 'key' 
    of the atom identified by 'uid' to the value 'value'. You
    *must* have already checked that this a valid thing to do,
    as this function does *no* error checking! */
void StructureEditor::_pvt_setAtomMetadata(quint32 uid, const QString &key,
                                           const QString &metakey,
                                           const QVariant &value)
{
    this->assertSane();

    d->atom(uid).property_metadata[key].insert(metakey, value);
}

/** Protected function used to set the property at key 'key' 
    of the residue identified by 'uid' to the value 'value'. You
    *must* have already checked that this a valid thing to do,
    as this function does *no* error checking! */
void StructureEditor::_pvt_setResProperty(quint32 uid, const QString &key, 
                                          const QVariant &value)
{
    this->assertSane();

    d->residue(uid).properties.insert(key, value);
}

/** Protected function used to set the metadata at metakey 'metakey' 
    of the residue identified by 'uid' to the value 'value'. You
    *must* have already checked that this a valid thing to do,
    as this function does *no* error checking! */
void StructureEditor::_pvt_setResMetadata(quint32 uid, const QString &metakey,
                                          const QVariant &value)
{
    this->assertSane();

    d->residue(uid).molecule_metadata.insert(metakey, value);
}
                          
/** Protected function used to set the metadata at metakey
    'metakey' of the property at key 'key' 
    of the residue identified by 'uid' to the value 'value'. You
    *must* have already checked that this a valid thing to do,
    as this function does *no* error checking! */
void StructureEditor::_pvt_setResMetadata(quint32 uid, const QString &key,
                                          const QString &metakey,
                                          const QVariant &value)
{
    this->assertSane();

    d->residue(uid).property_metadata[key].insert(metakey, value);
}

/** Protected function used to set the property at key 'key' 
    of the CutGroup identified by 'uid' to the value 'value'. You
    *must* have already checked that this a valid thing to do,
    as this function does *no* error checking! */
void StructureEditor::_pvt_setCGProperty(quint32 uid, const QString &key, 
                                         const QVariant &value)
{
    this->assertSane();

    d->cutGroup(uid).properties.insert(key, value);
}

/** Protected function used to set the metadata at metakey 'metakey' 
    of the CutGroup identified by 'uid' to the value 'value'. You
    *must* have already checked that this a valid thing to do,
    as this function does *no* error checking! */
void StructureEditor::_pvt_setCGMetadata(quint32 uid, const QString &metakey,
                                         const QVariant &value)
{
    this->assertSane();

    d->cutGroup(uid).molecule_metadata.insert(metakey, value);
}
                          
/** Protected function used to set the metadata at metakey
    'metakey' of the property at key 'key' 
    of the CutGroup identified by 'uid' to the value 'value'. You
    *must* have already checked that this a valid thing to do,
    as this function does *no* error checking! */
void StructureEditor::_pvt_setCGMetadata(quint32 uid, const QString &key,
                                         const QString &metakey,
                                         const QVariant &value)
{
    this->assertSane();

    d->cutGroup(uid).property_metadata[key].insert(metakey, value);
}

/** Protected function used to set the property at key 'key' 
    of the chain identified by 'uid' to the value 'value'. You
    *must* have already checked that this a valid thing to do,
    as this function does *no* error checking! */
void StructureEditor::_pvt_setChainProperty(quint32 uid, const QString &key, 
                                            const QVariant &value)
{
    this->assertSane();

    d->chain(uid).properties.insert(key, value);
}

/** Protected function used to set the metadata at metakey 'metakey' 
    of the chain identified by 'uid' to the value 'value'. You
    *must* have already checked that this a valid thing to do,
    as this function does *no* error checking! */
void StructureEditor::_pvt_setChainMetadata(quint32 uid, const QString &metakey,
                                            const QVariant &value)
{
    this->assertSane();

    d->chain(uid).molecule_metadata.insert(metakey, value);
}
                          
/** Protected function used to set the metadata at metakey
    'metakey' of the property at key 'key' 
    of the chain identified by 'uid' to the value 'value'. You
    *must* have already checked that this a valid thing to do,
    as this function does *no* error checking! */
void StructureEditor::_pvt_setChainMetadata(quint32 uid, const QString &key,
                                            const QString &metakey,
                                            const QVariant &value)
{
    this->assertSane();

    d->chain(uid).property_metadata[key].insert(metakey, value);
}

/** Protected function used to set the property at key 'key' 
    of the segment identified by 'uid' to the value 'value'. You
    *must* have already checked that this a valid thing to do,
    as this function does *no* error checking! */
void StructureEditor::_pvt_setSegProperty(quint32 uid, const QString &key, 
                                          const QVariant &value)
{
    this->assertSane();

    d->segment(uid).properties.insert(key, value);
}

/** Protected function used to set the metadata at metakey 'metakey' 
    of the segment identified by 'uid' to the value 'value'. You
    *must* have already checked that this a valid thing to do,
    as this function does *no* error checking! */
void StructureEditor::_pvt_setSegMetadata(quint32 uid, const QString &metakey,
                                          const QVariant &value)
{
    this->assertSane();

    d->segment(uid).molecule_metadata.insert(metakey, value);
}
                          
/** Protected function used to set the metadata at metakey
    'metakey' of the property at key 'key' 
    of the segment identified by 'uid' to the value 'value'. You
    *must* have already checked that this a valid thing to do,
    as this function does *no* error checking! */
void StructureEditor::_pvt_setSegMetadata(quint32 uid, const QString &key,
                                           const QString &metakey,
                                           const QVariant &value)
{
    this->assertSane();

    d->segment(uid).property_metadata[key].insert(metakey, value);
}
