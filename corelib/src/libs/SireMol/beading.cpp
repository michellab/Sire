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

#include "beading.h"
#include "beadidx.h"
#include "atomidx.h"
#include "atombeads.h"
#include "atomselection.h"
#include "moleculedata.h"
#include "moleculeinfodata.h"

#include <boost/noncopyable.hpp>

#include "SireError/errors.h"
#include "SireBase/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMol;
using namespace SireBase;
using namespace SireStream;

//////////
////////// Implementation of Beading
//////////

static const RegisterMetaType<Beading> r_beading( MAGIC_ONLY,
                                                  Beading::typeName() );
                                                  
QDataStream &operator<<(QDataStream &ds, const Beading &beading)
{
    writeHeader(ds, r_beading, 1);
    
    ds << static_cast<const MolViewProperty&>(beading);
    
    return ds; 
}

QDataStream &operator>>(QDataStream &ds, Beading &beading)
{
    VersionID v = readHeader(ds, r_beading);
    
    if (v == 1)
    {
        ds >> static_cast<MolViewProperty&>(beading);
    }
    else
        throw version_error(v, "1", r_beading, CODELOC);
        
    return ds;
}

/** Constructor */
Beading::Beading() : MolViewProperty()
{}

/** Copy constructor */
Beading::Beading(const Beading &other) : MolViewProperty(other)
{}

/** Destructor */
Beading::~Beading()
{}

void Beading::assertValidIndex(BeadIdx idx, const MoleculeInfoData &molinfo) const
{
    idx.map( this->nBeads(molinfo) );
}

bool Beading::isCompatibleWith(const MoleculeInfoData&) const
{
    return true;
}

/** Copy assignment operator */
Beading& Beading::operator=(const Beading &other)
{
    MolViewProperty::operator=(other);
    return *this;
}

/** Comparison operator */
bool Beading::operator==(const Beading &other) const
{
    return MolViewProperty::operator==(other);
}

/** Comparison operator */
bool Beading::operator!=(const Beading &other) const
{
    return MolViewProperty::operator!=(other);
}

/** By default the bead number is the index + 1 (e.g. if
    the BeadIdx is 0, the BeadNum is 1) */
BeadNum Beading::beadNum(const MoleculeInfoData &moldata,
                         const BeadIdx &bead) const
{
    return BeadNum( 1 + bead.map( this->nBeads(moldata) ) );
}

const char* Beading::typeName()
{
    return "SireMol::Beading";
}

NullBeading Beading::null()
{
    return NullBeading();
}

//////////
////////// Implementation of MoleculeBeading
//////////

static const RegisterMetaType<MoleculeBeading> r_molbeading;

QDataStream &operator<<(QDataStream &ds, const MoleculeBeading &molbeading)
{
    writeHeader(ds, r_molbeading, 1);
    
    ds << static_cast<const Beading&>(molbeading);
    
    return ds;
}

QDataStream &operator>>(QDataStream &ds, MoleculeBeading &molbeading)
{
    VersionID v = readHeader(ds, r_molbeading);
    
    if (v == 1)
    {
        ds >> static_cast<Beading&>(molbeading);
    }
    else
        throw version_error(v, "1", r_molbeading, CODELOC);
        
    return ds;
}

/** Constructor */
MoleculeBeading::MoleculeBeading() : ConcreteProperty<MoleculeBeading,Beading>()
{}

/** Copy constructor */
MoleculeBeading::MoleculeBeading(const MoleculeBeading &other)
                : ConcreteProperty<MoleculeBeading,Beading>(other)
{}

/** Destructor */
MoleculeBeading::~MoleculeBeading()
{}

/** Copy assignment operator */
MoleculeBeading& MoleculeBeading::operator=(const MoleculeBeading &other)
{
    Beading::operator=(other);
    return *this;
}

/** Comparison operator */
bool MoleculeBeading::operator==(const MoleculeBeading &other) const
{
    return Beading::operator==(other);
}

/** Comparison operator */
bool MoleculeBeading::operator!=(const MoleculeBeading &other) const
{
    return Beading::operator!=(other);
}

const char* MoleculeBeading::typeName()
{
    return QMetaType::typeName( qMetaTypeId<MoleculeBeading>() );
}

/** Return the number of beads */
int MoleculeBeading::nBeads(const MoleculeInfoData &moldata) const
{
    if (moldata.nAtoms() == 0)
        return 0;
    else
        return 1;
}

/** Return the index of the ith atom in the ith bead */
AtomIdx MoleculeBeading::atomIdx(const MoleculeInfoData &moldata,
                                 const BeadIdx &bead, int i) const
{
    bead.map( MoleculeBeading::nBeads(moldata) );
    
    return AtomIdx( AtomIdx(i).map(moldata.nAtoms()) );
}

/** Return the values of the atom property 'key' for all beads, 
    arranged in bead:index order */
PropertyPtr MoleculeBeading::atomProperty(const MoleculeData &moldata,
                                          const SireBase::PropertyName &key) const
{
    return moldata.property(key).asA<AtomProp>().merge(moldata.info());
}
                         
/** Return all of the atoms that are selected in the beads */
AtomSelection MoleculeBeading::selection(const MoleculeInfoData &moldata) const
{
    return AtomSelection(moldata);
}

/** Return the atoms that are part of the bead with index 'bead' 

    \throw SireError::invalid_index
*/
AtomSelection MoleculeBeading::selection(const MoleculeInfoData &moldata,
                                         const BeadIdx &bead) const
{
    bead.map( MoleculeBeading::nBeads(moldata) );
    return AtomSelection(moldata);
}

/** Return the indicies of all of the atoms in the beads */
QList<AtomIdx> MoleculeBeading::atomIdxs(const MoleculeInfoData &moldata) const
{
    QList<AtomIdx> atoms;
    
    for (int i=0; i<moldata.nAtoms(); ++i)
    {
        atoms.append( AtomIdx(i) );
    }
    
    return atoms;
}

/** Return the indicies of all of the atoms in the bead with index 'bead'

    \throw SireError::invalid_index
*/
QList<AtomIdx> MoleculeBeading::atomIdxs(const MoleculeInfoData &moldata,
                                         const BeadIdx &bead) const
{
    bead.map( MoleculeBeading::nBeads(moldata) );
    return MoleculeBeading::atomIdxs(moldata);
}

//////////
////////// Implementation of ResidueBeading
//////////

static const RegisterMetaType<ResidueBeading> r_resbeading;

QDataStream &operator<<(QDataStream &ds, const ResidueBeading &resbeading)
{
    writeHeader(ds, r_resbeading, 1);
    
    ds << static_cast<const Beading&>(resbeading);
    
    return ds;
}

QDataStream &operator>>(QDataStream &ds, ResidueBeading &resbeading)
{
    VersionID v = readHeader(ds, r_resbeading);
    
    if (v == 1)
    {
        ds >> static_cast<Beading&>(resbeading);
    }
    else
        throw version_error(v, "1", r_resbeading, CODELOC);
        
    return ds;
}

/** Constructor */
ResidueBeading::ResidueBeading() : ConcreteProperty<ResidueBeading,Beading>()
{}

/** Copy constructor */
ResidueBeading::ResidueBeading(const ResidueBeading &other)
               : ConcreteProperty<ResidueBeading,Beading>(other)
{}

/** Destructor */
ResidueBeading::~ResidueBeading()
{}

/** Copy assignment operator */
ResidueBeading& ResidueBeading::operator=(const ResidueBeading &other)
{
    Beading::operator=(other);
    return *this;
}

/** Comparison operator */
bool ResidueBeading::operator==(const ResidueBeading &other) const
{
    return Beading::operator==(other);
}

/** Comparison operator */
bool ResidueBeading::operator!=(const ResidueBeading &other) const
{
    return Beading::operator!=(other);
}

const char* ResidueBeading::typeName()
{
    return QMetaType::typeName( qMetaTypeId<ResidueBeading>() );
}

/** Return the bead number - for residue beads, the bead number
    is the same as the residue number */
BeadNum ResidueBeading::beadNum(const MoleculeInfoData &moldata,
                                const BeadIdx &bead) const
{
    return BeadNum( moldata.number( ResIdx(bead.value()) ).value() );
}

/** Return the number of beads */
int ResidueBeading::nBeads(const MoleculeInfoData &moldata) const
{
    return moldata.nResidues();
}

/** Return the index of the ith atom in the bead with index 'i' */
AtomIdx ResidueBeading::atomIdx(const MoleculeInfoData &moldata,
                                const BeadIdx &bead, int i) const
{
    return moldata.getAtom( ResIdx(bead.value()), i );
}

/** Return the atom properties in bead/atom order */
PropertyPtr ResidueBeading::atomProperty(const MoleculeData &moldata,
                                         const SireBase::PropertyName &key) const
{
    return moldata.property(key).asA<AtomProp>().divideByResidue(moldata.info());
}
                         
/** Return the atoms that are in all of the beads */
AtomSelection ResidueBeading::selection(const MoleculeInfoData &moldata) const
{
    return AtomSelection(moldata);
}

/** Return the atoms that are in the bead with index 'bead' 

    \throw SireError::invalid_index
*/                         
AtomSelection ResidueBeading::selection(const MoleculeInfoData &moldata,
                                        const BeadIdx &bead) const
{
    AtomSelection selected_atoms(moldata);
    return selected_atoms.selectNone().select( ResIdx(bead.value()) );
}

/** Return the indexes of all of the atoms in the beads */
QList<AtomIdx> ResidueBeading::atomIdxs(const MoleculeInfoData &moldata) const
{
    QList<AtomIdx> atoms;
    
    for (int i=0; i<moldata.nAtoms(); ++i)
    {
        atoms.append( AtomIdx(i) );
    }
    
    return atoms;
}

/** Return the indexes of the atoms in the bead with index 'bead'

    \throw SireError::invalid_index
*/
QList<AtomIdx> ResidueBeading::atomIdxs(const MoleculeInfoData &moldata,
                                        const BeadIdx &bead) const
{
    return moldata.getAtomsIn( ResIdx(bead.value()) );
}

//////////
////////// Implementation of UserBeading
//////////

namespace SireMol{ namespace detail {

    /** This is an internal helper class used to cache the information
        required for UserBeading */
    class UserBeadingInfo : public boost::noncopyable
    {
    public:
        UserBeadingInfo();
    
        UserBeadingInfo(const AtomBeads &beads,
                        const MoleculeInfoData &molinfo);
                        
        ~UserBeadingInfo();

        int nBeads() const;

        BeadNum beadNum(const BeadIdx &bead) const;

        AtomIdx atomIdx(const BeadIdx &bead, int i) const;

        PropertyPtr atomProperty(const MoleculeData &moldata,
                                 const SireBase::PropertyName &key) const;
                                 
        AtomSelection selection() const;
                                 
        AtomSelection selection(const BeadIdx &bead) const;

        QList<AtomIdx> atomIdxs() const;
        QList<AtomIdx> atomIdxs(const BeadIdx &bead) const;

    private:
        /** The atoms selected in each bead */
        QVector<AtomSelection> bead_atoms;
        
        /** All of the atoms in the beads */
        AtomSelection all_atoms;

        /** The map of bead index to number */
        QVector<BeadNum> beadidx_to_num;
        
        /** The indicies of the atoms in each bead */
        QVector< QVector<AtomIdx> > bead_atomidxs;
    };
    
    /** This internal class holds the cache of UserBeadInfos
        for each MoleculeInfoData object */
    class UserBeadingInfoRegistry : public boost::noncopyable
    {
    public:
        UserBeadingInfoRegistry(const AtomBeads &beads);
        ~UserBeadingInfoRegistry();
        
        const AtomBeads& atomBeads() const;
        
        const UserBeadingInfo& getInfoFor(const MoleculeInfoData &moldata);
        
    private:
        /** The AtomBeads to use to bead up the molecule */
        AtomBeads atom_beads;
    
        /** The actual database of beading infos */
        QHash< QUuid,boost::shared_ptr<UserBeadingInfo> > reg;
        
        /** Lock to protect access to the database */
        QMutex datamutex;
    };

} } // end of namespaces SireMol::detail and SireMol

using namespace SireMol::detail;

UserBeadingInfoRegistry::UserBeadingInfoRegistry(const AtomBeads &beads)
                        : boost::noncopyable(), atom_beads(beads)
{}

UserBeadingInfoRegistry::~UserBeadingInfoRegistry()
{
    QMutexLocker lkr(&datamutex);
    reg.clear();
}

const AtomBeads& UserBeadingInfoRegistry::atomBeads() const
{
    return atom_beads;
}
        
const UserBeadingInfo& UserBeadingInfoRegistry::getInfoFor(
                                        const MoleculeInfoData &moldata)
{
    if (atom_beads.isEmpty())
    {
        static const UserBeadingInfo default_info;
        return default_info;
    }

    QMutexLocker lkr(&datamutex);
    
    if (not reg.contains(moldata.UID()))
    {
        lkr.unlock();
        
        boost::shared_ptr<UserBeadingInfo> info(
                                new UserBeadingInfo(atom_beads, moldata) );
                                
        lkr.relock();
        
        if (not reg.contains(moldata.UID()))
            reg.insert(moldata.UID(), info);
    }
    
    return *(reg.value(moldata.UID()));
}

/** Null constructor */
UserBeadingInfo::UserBeadingInfo() : boost::noncopyable()
{}

/** Construct the bead information for the passed atom beads applied
    to the passed molecule info */
UserBeadingInfo::UserBeadingInfo(const AtomBeads &beads,
                                 const MoleculeInfoData &molinfo)
                : boost::noncopyable()
{
    throw SireError::incomplete_code( "TODO", CODELOC );
}

/** Destructor */                
UserBeadingInfo::~UserBeadingInfo()
{}

/** Return the number of beads */
int UserBeadingInfo::nBeads() const
{
    return beadidx_to_num.count();
}

/** Return the number of the bead with index 'bead' 

    \throw SireError::invalid_index
*/
BeadNum UserBeadingInfo::beadNum(const BeadIdx &bead) const
{
    return beadidx_to_num.constData()[ bead.map(beadidx_to_num.count()) ];
}

/** Return the index of the ith atom in the bead with index 'bead'

    \throw SireError::invalid_index
*/
AtomIdx UserBeadingInfo::atomIdx(const BeadIdx &bead, int i) const
{
    const QVector<AtomIdx> &atoms = bead_atomidxs.constData()
                                            [ bead.map(bead_atomidxs.count()) ];
                                            
    return atoms.constData()[ Index(i).map(atoms.count()) ];
}

/** Return the property with key 'key' of all of the atoms in the beads
    in bead/index order */
PropertyPtr UserBeadingInfo::atomProperty(const MoleculeData &moldata,
                                          const SireBase::PropertyName &key) const
{
    return moldata.property(key).asA<AtomProp>().divide(bead_atoms);
}
                         
/** Return the selection of all of the atoms that are
    in all of the beads */
AtomSelection UserBeadingInfo::selection() const
{
    return all_atoms;
}
                     
/** Return the selection of the atoms that are in the bead with
    index 'bead'
    
    \throw SireError::invalid_index
*/
AtomSelection UserBeadingInfo::selection(const BeadIdx &bead) const
{
    return bead_atoms.constData()[ bead.map(bead_atoms.count()) ];
}

/** Return the indices of all of the atoms in all of the beads */
QList<AtomIdx> UserBeadingInfo::atomIdxs() const
{
    QList<AtomIdx> atoms;
    
    for (int i=0; i<bead_atomidxs.count(); ++i)
    {
        atoms += bead_atomidxs.constData()[i].toList();
    }
    
    return atoms;
}

/** Return the indicies of the atoms that are in the bead 
    with index 'bead'
    
    \throw SireError::invalid_index
*/
QList<AtomIdx> UserBeadingInfo::atomIdxs(const BeadIdx &bead) const
{
    return bead_atomidxs.constData()[ bead.map(bead_atomidxs.count()) ].toList();
}

static const RegisterMetaType<UserBeading> r_userbeading;

QDataStream &operator<<(QDataStream &ds, const UserBeading &userbeading)
{
    writeHeader(ds, r_userbeading, 1);
    
    SharedDataStream sds(ds);
    
    if (userbeading.registry.get() == 0)
    {
        sds << AtomBeads();
    }
    else
    {
        sds << userbeading.registry->atomBeads();
    }
    
    sds << static_cast<const Beading&>(userbeading);
    
    return ds;
}

QDataStream &operator>>(QDataStream &ds, UserBeading &userbeading)
{
    VersionID v = readHeader(ds, r_userbeading);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        AtomBeads atom_beads;
        
        sds >> atom_beads;
        
        userbeading = UserBeading(atom_beads);
        
        sds >> static_cast<Beading&>(userbeading);
    }
    else
        throw version_error(v, "1", r_userbeading, CODELOC);
        
    return ds;
}

/** Constructor - this looks for the bead property in "bead" */
UserBeading::UserBeading()
            : ConcreteProperty<UserBeading,Beading>()
{}

/** Constructor used to specify the beads for each atom */
UserBeading::UserBeading(const AtomBeads &beads)
            : ConcreteProperty<UserBeading,Beading>()
{
    if (not beads.isEmpty())
    {
        registry.reset( new UserBeadingInfoRegistry(beads) );
    }
}

/** Copy constructor */
UserBeading::UserBeading(const UserBeading &other)
            : ConcreteProperty<UserBeading,Beading>(other),
              registry(other.registry)
{}

/** Destructor */
UserBeading::~UserBeading()
{}

/** Copy assignment operator */
UserBeading& UserBeading::operator=(const UserBeading &other)
{
    if (this != &other)
    {
        registry = other.registry;
        Beading::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool UserBeading::operator==(const UserBeading &other) const
{
    return this->atomBeads() == other.atomBeads() and
           Beading::operator==(other);
}

/** Comparison operator */
bool UserBeading::operator!=(const UserBeading &other) const
{
    return not UserBeading::operator==(other);
}

const char* UserBeading::typeName()
{
    return QMetaType::typeName( qMetaTypeId<UserBeading>() );
}

/** Return the bead specification for each atom */
const AtomBeads& UserBeading::atomBeads() const
{
    if (registry.get() == 0)
    {
        static const AtomBeads atom_beads;
        return atom_beads;
    }
    else
    {
        return registry->atomBeads();
    }
}

/** Return whether or not this beading is compatible with the passed molecule info */
bool UserBeading::isCompatibleWith(const MoleculeInfoData &molinfo) const
{
    return this->atomBeads().isCompatibleWith(molinfo);
}

const UserBeadingInfo& UserBeading::getUserBeadingInfo(
                                const MoleculeInfoData &moldata) const
{
    if (registry.get() == 0)
    {
        static const UserBeadingInfo default_info;
        return default_info;
    }
    else
    {
        return registry->getInfoFor(moldata);
    }
}

/** Return the number of the bead with index 'bead'

    \throw SireError::invalid_index
*/
BeadNum UserBeading::beadNum(const MoleculeInfoData &moldata, const BeadIdx &bead) const
{
    return getUserBeadingInfo(moldata).beadNum(bead);
}

/** Return the number of beads in the passed molecule */
int UserBeading::nBeads(const MoleculeInfoData &moldata) const
{
    return getUserBeadingInfo(moldata).nBeads();
}

/** Return the index of the ith atom in the bead with index 'bead'

    \throw SireError::invalid_index
*/
AtomIdx UserBeading::atomIdx(const MoleculeInfoData &moldata,
                             const BeadIdx &bead, int i) const
{
    return getUserBeadingInfo(moldata).atomIdx(bead, i);
}

/** Return the properties of all of the atoms of all of the beads,
    in bead index, atom in bead index order */
PropertyPtr UserBeading::atomProperty(const MoleculeData &moldata,
                                      const SireBase::PropertyName &key) const
{
    return getUserBeadingInfo(moldata.info()).atomProperty(moldata, key);
}
                         
/** Return the selection of all of the atoms in all of the beads */
AtomSelection UserBeading::selection(const MoleculeInfoData &moldata) const
{
    return getUserBeadingInfo(moldata).selection();
}
                         
/** Return the selection of atoms in the bead with index 'bead'

    \throw SireError::invalid_index
*/
AtomSelection UserBeading::selection(const MoleculeInfoData &moldata,
                                     const BeadIdx &bead) const
{
    return getUserBeadingInfo(moldata).selection(bead);
}

/** Return the indicies of all of the atoms in all of the beads */
QList<AtomIdx> UserBeading::atomIdxs(const MoleculeInfoData &moldata) const
{
    return getUserBeadingInfo(moldata).atomIdxs();
}

/** Return the indicies of the atoms in the bead with index 'bead'

    \throw SireError::invalid_index
*/
QList<AtomIdx> UserBeading::atomIdxs(const MoleculeInfoData &moldata,
                                     const BeadIdx &bead) const
{
    return getUserBeadingInfo(moldata).atomIdxs(bead);
}

//////////
////////// Implementation of NullBeading
//////////

static const RegisterMetaType<NullBeading> r_nullbeading;

QDataStream &operator<<(QDataStream &ds, const NullBeading &nullbeading)
{
    writeHeader(ds, r_nullbeading, 1);
    
    ds << static_cast<const Beading&>(nullbeading);
    
    return ds;
}

QDataStream &operator>>(QDataStream &ds, NullBeading &nullbeading)
{
    VersionID v = readHeader(ds, r_nullbeading);
    
    if (v == 1)
    {
        ds >> static_cast<Beading&>(nullbeading);
    }
    else
        throw version_error(v, "1", r_nullbeading, CODELOC);
        
    return ds;
}

/** Constructor */
NullBeading::NullBeading() : ConcreteProperty<NullBeading,Beading>()
{}

/** Copy constructor */
NullBeading::NullBeading(const NullBeading &other)
            : ConcreteProperty<NullBeading,Beading>(other)
{}

/** Destructor */
NullBeading::~NullBeading()
{}

/** Copy assignment operator */
NullBeading& NullBeading::operator=(const NullBeading &other)
{
    Beading::operator=(other);
    return *this;
}

/** Comparison operator */
bool NullBeading::operator==(const NullBeading &other) const
{
    return Beading::operator==(other);
}

/** Comparison operator */
bool NullBeading::operator!=(const NullBeading &other) const
{
    return Beading::operator!=(other);
}

const char* NullBeading::typeName()
{
    return QMetaType::typeName( qMetaTypeId<NullBeading>() );
}

/** Return the number of beads in the passed molecule */
int NullBeading::nBeads(const MoleculeInfoData &moldata) const
{
    return 0;
}

/** Return the AtomIdx of the ith atom in the bead with index 'bead' */
AtomIdx NullBeading::atomIdx(const MoleculeInfoData &moldata,
                             const BeadIdx &bead, int i) const
{
    throw SireError::invalid_index( QObject::tr(
            "NullBeading creates zero beads, so there are no atoms in any bead!"),
                CODELOC );
                
    return AtomIdx();
}

/** Return the atom properties of all of the beads in bead, then atom order */
PropertyPtr NullBeading::atomProperty(const MoleculeData &moldata,
                                      const PropertyName &key) const
{
    throw SireBase::missing_property( QObject::tr(
            "There are no Beads using a NullBeading, so no property with "
            "key %1.").arg(key.toString()), CODELOC );
            
    return PropertyPtr();
}
                         
/** Return the atoms from the molecule that are part of the beads */
AtomSelection NullBeading::selection(const MoleculeInfoData &moldata) const
{
    return AtomSelection();
}
                         
/** Return the atoms from the molecule that in the bead with index 'bead'

    \throw SireError::invalid_index
*/
AtomSelection NullBeading::selection(const MoleculeInfoData &moldata,
                                     const BeadIdx &bead) const
{
    throw SireError::invalid_index( QObject::tr(
            "The NullBeading beader creates no beads, so there is no bead with "
            "index %1.").arg(bead.value()), CODELOC );
            
    return AtomSelection();
}

/** Return the atom indexes of the atoms that are part of the beads */
QList<AtomIdx> NullBeading::atomIdxs(const MoleculeInfoData &moldata) const
{
    return QList<AtomIdx>();
}

/** Return the atom indexes that are part of the bead with index 'bead' */
QList<AtomIdx> NullBeading::atomIdxs(const MoleculeInfoData &moldata,
                                     const BeadIdx &bead) const
{
    throw SireError::invalid_index( QObject::tr(
            "The NullBeading beader creates no beads, so there is no bead with "
            "index %1.").arg(bead), CODELOC );

    return QList<AtomIdx>();
}

