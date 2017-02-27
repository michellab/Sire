/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2009  Christopher Woods
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

#include "moldeleter.h"
#include "uniformsampler.h"

#include "SireMol/moleculegroup.h"
#include "SireMol/molecule.h"

#include "SireSystem/system.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMove;
using namespace SireSystem;
using namespace SireMol;
using namespace SireMaths;
using namespace SireBase;
using namespace SireStream;

/////////////
///////////// Implementation of MolDeleter
/////////////

static const RegisterMetaType<MolDeleter> r_moldeleter( MAGIC_ONLY, 
                                                        MolDeleter::typeName() );

/** Serialise to a binary datastream */
QDataStream SIREMOVE_EXPORT &operator<<(QDataStream &ds,
                                        const MolDeleter &moldeleter)
{
    writeHeader(ds, r_moldeleter, 1);
    
    ds << static_cast<const Property&>(moldeleter);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOVE_EXPORT &operator>>(QDataStream &ds, MolDeleter &moldeleter)
{
    VersionID v = readHeader(ds, r_moldeleter);
    
    if (v == 1)
    {
        ds >> static_cast<Property&>(moldeleter);
    }
    else
        throw version_error(v, "1", r_moldeleter, CODELOC);
        
    return ds;
}

/** Constructor */
MolDeleter::MolDeleter() : Property()
{}

/** Copy constructor */
MolDeleter::MolDeleter(const MolDeleter &other) : Property(other)
{}

/** Destructor */
MolDeleter::~MolDeleter()
{}

/** Copy assignment operator */
MolDeleter& MolDeleter::operator=(const MolDeleter &other)
{
    Property::operator=(other);
    return *this;
}

/** Comparison operator */
bool MolDeleter::operator==(const MolDeleter &other) const
{
    return true;
}

/** Comparison operator */
bool MolDeleter::operator!=(const MolDeleter &other) const
{
    return false;
}

const char* MolDeleter::typeName()
{
    return "SireMove::MolDeleter";
}

/** Return the global null MolDeleter */
const NullDeleter& MolDeleter::null()
{
    return *(create_shared_null<NullDeleter>());
}

/////////////
///////////// Implementation of NullDeleter
/////////////

static const RegisterMetaType<NullDeleter> r_nulldeleter;

/** Serialise to a binary datastream */
QDataStream SIREMOVE_EXPORT &operator<<(QDataStream &ds, const NullDeleter &nulldeleter)
{
    writeHeader(ds, r_nulldeleter, 1);
    
    ds << static_cast<const MolDeleter&>(nulldeleter);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOVE_EXPORT &operator>>(QDataStream &ds, NullDeleter &nulldeleter)
{
    VersionID v = readHeader(ds, r_nulldeleter);
    
    if (v == 1)
    {
        ds >> static_cast<MolDeleter&>(nulldeleter);
    }
    else
        throw version_error(v, "1", r_nulldeleter, CODELOC);

    return ds;
}

/** Constructor */
NullDeleter::NullDeleter() : ConcreteProperty<NullDeleter,MolDeleter>()
{}

/** Copy constructor */
NullDeleter::NullDeleter(const NullDeleter &other)
            : ConcreteProperty<NullDeleter,MolDeleter>(other)
{}

/** Destructor */
NullDeleter::~NullDeleter()
{}

const char* NullDeleter::typeName()
{
    return QMetaType::typeName( qMetaTypeId<NullDeleter>() );
}

/** Copy operator */
NullDeleter& NullDeleter::operator=(const NullDeleter &other)
{
    MolDeleter::operator=(other);
    return *this;
}

/** Comparison operator */
bool NullDeleter::operator==(const NullDeleter &other) const
{
    return MolDeleter::operator==(other);
}

/** Comparison operator */
bool NullDeleter::operator!=(const NullDeleter &other) const
{
    return MolDeleter::operator!=(other);
}

/** Set the generator used to select random molecules 
    (this does nothing!) */
void NullDeleter::setGenerator(const RanGenerator&)
{}

/** Return the generator used to select random molecules
    (this just returns the global generator) */
const RanGenerator& NullDeleter::generator() const
{
    return RanGenerator::global();
}

/** Delete a molecule from the system - well this does nothing too! */
tuple<Molecule,double> NullDeleter::deleteFrom(System&)
{
    return tuple<Molecule,double>( Molecule(), 0.0 );
}

/////////////
///////////// Implementation of SystemWideDeleter
/////////////

static const RegisterMetaType<SystemWideDeleter> r_syswide;

/** Serialise to a binary datastream */
QDataStream SIREMOVE_EXPORT &operator<<(QDataStream &ds,
                                        const SystemWideDeleter &syswide)
{
    writeHeader(ds, r_syswide, 1);
    
    SharedDataStream sds(ds);
    
    sds << syswide.smplr << static_cast<const MolDeleter&>(syswide);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOVE_EXPORT &operator>>(QDataStream &ds,
                                        SystemWideDeleter &syswide)
{
    VersionID v = readHeader(ds, r_syswide);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> syswide.smplr >> static_cast<MolDeleter&>(syswide);
    }
    else
        throw version_error(v, "1", r_syswide, CODELOC);
        
    return ds;
}

/** Constructor */
SystemWideDeleter::SystemWideDeleter()
                  : ConcreteProperty<SystemWideDeleter,MolDeleter>()
{}

/** Construct a deleter that deletes molecules uniformly chosen from
    the passed molecule group from the entire system */
SystemWideDeleter::SystemWideDeleter(const MoleculeGroup &molgroup)
                  : ConcreteProperty<SystemWideDeleter,MolDeleter>(),
                    smplr( UniformSampler(molgroup) )
{}

/** Construct a deleter that deletes molecules chosen using the 
    passed sampler from the entire system */
SystemWideDeleter::SystemWideDeleter(const Sampler &sampler)
                  : ConcreteProperty<SystemWideDeleter,MolDeleter>(),
                    smplr(sampler)
{}

/** Copy constructor */
SystemWideDeleter::SystemWideDeleter(const SystemWideDeleter &other)
                  : ConcreteProperty<SystemWideDeleter,MolDeleter>(other),
                    smplr(other.smplr)
{}

/** Destructor */
SystemWideDeleter::~SystemWideDeleter()
{}

const char* SystemWideDeleter::typeName()
{
    return QMetaType::typeName( qMetaTypeId<SystemWideDeleter>() );
}

/** Copy assignment operator */
SystemWideDeleter& SystemWideDeleter::operator=(const SystemWideDeleter &other)
{
    if (this != &other)
    {
        smplr = other.smplr;
        MolDeleter::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool SystemWideDeleter::operator==(const SystemWideDeleter &other) const
{
    return smplr == other.smplr and MolDeleter::operator==(other);
}

/** Comparison operator */
bool SystemWideDeleter::operator!=(const SystemWideDeleter &other) const
{
    return not SystemWideDeleter::operator==(other);
}

/** Set the random number generator used by the sampler */
void SystemWideDeleter::setGenerator(const RanGenerator &generator)
{
    smplr.edit().setGenerator(generator);
}

/** Return the random number generator used by the sampler */
const RanGenerator& SystemWideDeleter::generator() const
{
    return smplr.read().generator();
}

/** Return the sampler used to pick molecules to be deleted */
const Sampler& SystemWideDeleter::sampler() const
{
    return smplr.read();
}

/** Return the molecule group from which molecules to be 
    deleted are chosen */
const MoleculeGroup& SystemWideDeleter::group() const
{
    return smplr.read().group();
}

/** Set the sampler used to pick molecules to be deleted  */
void SystemWideDeleter::setSampler(const Sampler &sampler)
{
    smplr = sampler;
}

/** Set the sampler to be the one that selects molecules uniformly
    from the passed molecule group */
void SystemWideDeleter::setSampler(const MoleculeGroup &molgroup)
{
    smplr = UniformSampler(molgroup);
}

/** Set the molecule group that will be sampled by the sampler */
void SystemWideDeleter::setGroup(const MoleculeGroup &molgroup)
{
    smplr.edit().setGroup(molgroup);
}

/** Delete a molecule from the system. This returns the molecule that
    was deleted, and the probability with which it was sampled
    (normalised so that a probability of 1 is returned if the molecule
     was picked purely randomly). This deleter deletes the molecule
     from the entire system. This returns an empty molecule if
     the molecule was not contained in the system and nothing 
     was deleted */
tuple<Molecule,double> SystemWideDeleter::deleteFrom(System &system)
{
    smplr.edit().updateFrom(system);
    
    tuple<Molecule,double> picked_mol = smplr.read().sampleMolecule();
    
    if (picked_mol.get<1>() == 0)
        //nothing was picked
        return picked_mol;
        
    else if ( system.remove(picked_mol.get<0>().number()) )
    {
        //now normalise the picking probability
        picked_mol.get<1>() /= group().nMolecules();
        
        //update the sampler again, as it now has fewer molecules
        smplr.edit().updateFrom(system);
        
        return picked_mol;
    }
    else 
        //the picked molecule isn't in the system - this is an edge
        //case where the molecule group in the sampler is not actually
        //in the system either!
        return tuple<Molecule,double>( Molecule(), 0.0 );
}

/////////////
///////////// Implementation of SpecifiedGroupsDeleter
/////////////

static const RegisterMetaType<SpecifiedGroupsDeleter> r_specgroups;

/** Serialise to a binary datastream */
QDataStream SIREMOVE_EXPORT &operator<<(QDataStream &ds,
                                        const SpecifiedGroupsDeleter &specgroups)
{
    writeHeader(ds, r_specgroups, 1);
    
    SharedDataStream sds(ds);
    
    sds << specgroups.smplr << specgroups.mgid
        << static_cast<const MolDeleter&>(specgroups);
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOVE_EXPORT &operator>>(QDataStream &ds,
                                        SpecifiedGroupsDeleter &specgroups)
{
    VersionID v = readHeader(ds, r_specgroups);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> specgroups.smplr >> specgroups.mgid
            >> static_cast<MolDeleter&>(specgroups);
    }
    else
        throw version_error(v, "1", r_specgroups, CODELOC);
        
    return ds;
}

/** Constructor */
SpecifiedGroupsDeleter::SpecifiedGroupsDeleter()
                       : ConcreteProperty<SpecifiedGroupsDeleter,MolDeleter>()
{}

/** Construct a deleter that deletes molecules selected at random from the 
    molecule group 'molgroup' from the molecule groups that are identified
    by the ID in 'mgid' */
SpecifiedGroupsDeleter::SpecifiedGroupsDeleter(const MoleculeGroup &molgroup, 
                                               const MGID &mgroup_id)
                       : ConcreteProperty<SpecifiedGroupsDeleter,MolDeleter>(),
                         smplr( UniformSampler(molgroup) ),
                         mgid(mgroup_id)
{}

/** Construct a deleter that deletes molecules selected at random using the 
    passed sampler from the molecule groups that are identified by the 
    ID in 'mgid' */
SpecifiedGroupsDeleter::SpecifiedGroupsDeleter(const Sampler &sampler, 
                                               const MGID &mgroup_id)
                       : ConcreteProperty<SpecifiedGroupsDeleter,MolDeleter>(),
                         smplr(sampler), mgid(mgroup_id)
{}

/** Copy constructor */
SpecifiedGroupsDeleter::SpecifiedGroupsDeleter(const SpecifiedGroupsDeleter &other)
                       : ConcreteProperty<SpecifiedGroupsDeleter,MolDeleter>(other),
                         smplr(other.smplr), mgid(other.mgid)
{}

/** Destructor */
SpecifiedGroupsDeleter::~SpecifiedGroupsDeleter()
{}

const char* SpecifiedGroupsDeleter::typeName()
{
    return QMetaType::typeName( qMetaTypeId<SpecifiedGroupsDeleter>() );
}

/** Copy assignment operator */
SpecifiedGroupsDeleter& SpecifiedGroupsDeleter::operator=(
                                            const SpecifiedGroupsDeleter &other)
{
    if (this != &other)
    {
        smplr = other.smplr;
        mgid = other.mgid;
        MolDeleter::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool SpecifiedGroupsDeleter::operator==(const SpecifiedGroupsDeleter &other) const
{
    return smplr == other.smplr and mgid == other.mgid and
           MolDeleter::operator==(other);
}

/** Comparison operator */
bool SpecifiedGroupsDeleter::operator!=(const SpecifiedGroupsDeleter &other) const
{
    return not SpecifiedGroupsDeleter::operator==(other);
}

/** Set the random number generator used by the sampler to randomly
    select molecules to be deleted */
void SpecifiedGroupsDeleter::setGenerator(const RanGenerator &generator)
{
    smplr.edit().setGenerator(generator);
}

/** Return the random number generator used by the sampler to randomly
    select molecules to be deleted */
const RanGenerator& SpecifiedGroupsDeleter::generator() const
{
    return smplr.read().generator();
}

/** Return the sampler used to randomly select molecules to be deleted */
const Sampler& SpecifiedGroupsDeleter::sampler() const
{
    return smplr.read();
}

/** Return the molecule group from which molecules are randomly
    selected to be deleted */
const MoleculeGroup& SpecifiedGroupsDeleter::group() const
{
    return smplr.read().group();
}

/** Set the sampler used to randomly select molecules to be deleted */
void SpecifiedGroupsDeleter::setSampler(const Sampler &sampler)
{
    smplr = sampler;
}

/** Set the sampler to one that selects molecules uniformly from the
    passed molecule group */
void SpecifiedGroupsDeleter::setSampler(const MoleculeGroup &molgroup)
{
    smplr = UniformSampler(molgroup);
}

/** Set the molecule group from which molecules are selected randomly */
void SpecifiedGroupsDeleter::setGroup(const MoleculeGroup &molgroup)
{
    smplr.edit().setGroup(molgroup);
}

/** Return the ID of the groups from which molecules are actually deleted */
const MGID& SpecifiedGroupsDeleter::specifiedGroups() const
{
    return mgid.base();
}

/** Set the ID of the groups from which molecules are actually deleted */
void SpecifiedGroupsDeleter::setSpecifiedGroups(const MGID &mgroup_id)
{
    mgid = mgroup_id;
}

/** Delete a molecule from the system. This returns the molecule that
    was deleted, and the probability with which it was sampled
    (normalised so that a probability of 1 is returned if the molecule
     was picked purely randomly). This deleter deletes the molecule
     from the entire system. This returns an empty molecule if
     the molecule was not contained in the system and nothing 
     was deleted */
tuple<Molecule,double> SpecifiedGroupsDeleter::deleteFrom(System &system)
{
    smplr.edit().updateFrom(system);
    
    tuple<Molecule,double> picked_mol = smplr.read().sampleMolecule();
    
    if (picked_mol.get<1>() == 0)
        //nothing was picked
        return picked_mol;
        
    else if ( system.remove(picked_mol.get<0>().number(), mgid) )
    {
        //now normalise the picking probability
        picked_mol.get<1>() /= group().nMolecules();
        
        //update the sampler again, as it now has fewer molecules
        smplr.edit().updateFrom(system);
        
        return picked_mol;
    }
    else 
        //the picked molecule isn't in the system - this is an edge
        //case where the molecule group in the sampler is not actually
        //in the system either!
        return tuple<Molecule,double>( Molecule(), 0.0 );
}
