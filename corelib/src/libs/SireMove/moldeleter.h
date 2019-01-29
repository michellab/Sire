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

#ifndef SIREMOVE_MOLDELETER_H
#define SIREMOVE_MOLDELETER_H

#include "SireBase/property.h"

#include "SireMaths/rangenerator.h"

#include "SireMol/mgidentifier.h"

#include "sampler.h"

SIRE_BEGIN_HEADER

namespace SireMove
{
class MolDeleter;
class NullDeleter;

class SystemWideDeleter;
class SpecifiedGroupsDeleter;
}

SIREMOVE_EXPORT QDataStream& operator<<(QDataStream&, const SireMove::MolDeleter&);
SIREMOVE_EXPORT QDataStream& operator>>(QDataStream&, SireMove::MolDeleter&);

SIREMOVE_EXPORT QDataStream& operator<<(QDataStream&, const SireMove::NullDeleter&);
SIREMOVE_EXPORT QDataStream& operator>>(QDataStream&, SireMove::NullDeleter&);

SIREMOVE_EXPORT QDataStream& operator<<(QDataStream&, const SireMove::SpecifiedGroupsDeleter&);
SIREMOVE_EXPORT QDataStream& operator>>(QDataStream&, SireMove::SpecifiedGroupsDeleter&);

SIREMOVE_EXPORT QDataStream& operator<<(QDataStream&, const SireMove::SystemWideDeleter&);
SIREMOVE_EXPORT QDataStream& operator>>(QDataStream&, SireMove::SystemWideDeleter&);

namespace SireSystem
{
class System;
}

namespace SireMove
{

using SireMol::MGID;
using SireMol::MGIdentifier;

using SireMaths::RanGenerator;

using SireSystem::System;

/** This is the base class of all molecule deleters. A molecule deleter
    is a manipulator that deletes (removes) molecules from Systems
    (or just from molecule groups within a System)
    
    @author Christopher Woods
*/
class SIREMOVE_EXPORT MolDeleter : public SireBase::Property
{

friend QDataStream& ::operator<<(QDataStream&, const MolDeleter&);
friend QDataStream& ::operator>>(QDataStream&, MolDeleter&);

public:
    MolDeleter();
    MolDeleter(const MolDeleter &other);
    
    virtual ~MolDeleter();
    
    static const char* typeName();
    
    virtual MolDeleter* clone() const=0;

    virtual void setGenerator(const RanGenerator &generator)=0;
    virtual const RanGenerator& generator() const=0;
    
    virtual tuple<Molecule,double> deleteFrom(System &system)=0;
    
    static const NullDeleter& null();
    
protected:
    MolDeleter& operator=(const MolDeleter &other);
    
    bool operator==(const MolDeleter &other) const;
    bool operator!=(const MolDeleter &other) const;

private:
    /** The random number generator used by the sampler */
    RanGenerator rangen;
};

/** This is a null deleter - this deletes nothing! */
class SIREMOVE_EXPORT NullDeleter 
            : public SireBase::ConcreteProperty<NullDeleter,MolDeleter>
{

friend QDataStream& ::operator<<(QDataStream&, const NullDeleter&);
friend QDataStream& ::operator>>(QDataStream&, NullDeleter&);

public:
    NullDeleter();
    NullDeleter(const NullDeleter &other);
    
    ~NullDeleter();
    
    static const char* typeName();
    
    NullDeleter& operator=(const NullDeleter &other);
    
    bool operator==(const NullDeleter &other) const;
    bool operator!=(const NullDeleter &other) const;

    void setGenerator(const RanGenerator &generator);
    const RanGenerator& generator() const;
    
    tuple<Molecule,double> deleteFrom(System &system);
};

/** This deleter selects a molecule at random from an identified
    molecule group (using the passed sampler) and then deletes
    that molecule completely from the system
    
    @author Christopher Woods
*/
class SIREMOVE_EXPORT SystemWideDeleter
        : public SireBase::ConcreteProperty<SystemWideDeleter,MolDeleter>
{

friend QDataStream& ::operator<<(QDataStream&, const SystemWideDeleter&);
friend QDataStream& ::operator>>(QDataStream&, SystemWideDeleter&);

public:
    SystemWideDeleter();

    SystemWideDeleter(const MoleculeGroup &molgroup);
    SystemWideDeleter(const Sampler &sampler);
    
    SystemWideDeleter(const SystemWideDeleter &other);
    
    ~SystemWideDeleter();
    
    static const char* typeName();
    
    SystemWideDeleter& operator=(const SystemWideDeleter &other);
    
    bool operator==(const SystemWideDeleter &other) const;
    bool operator!=(const SystemWideDeleter &other) const;

    void setGenerator(const RanGenerator &generator);
    const RanGenerator& generator() const;
    
    const Sampler& sampler() const;
    const MoleculeGroup& group() const;
    
    void setSampler(const Sampler &sampler);
    void setSampler(const MoleculeGroup &molgroup);
    
    void setGroup(const MoleculeGroup &molgroup);
    
    tuple<Molecule,double> deleteFrom(System &system);

private:
    /** The sampler that is used to select a molecule to be deleted */
    SamplerPtr smplr;
};

/** This is a molecule deleter that selects a molecule at random
    from a specified molecule group (using the contained sampler)
    and then deletes that molecule from a specific set of 
    molecule groups in the system
    
    @author Christopher Woods
*/
class SIREMOVE_EXPORT SpecifiedGroupsDeleter
            : public SireBase::ConcreteProperty<SpecifiedGroupsDeleter,MolDeleter>
{

friend QDataStream& ::operator<<(QDataStream&, const SpecifiedGroupsDeleter&);
friend QDataStream& ::operator>>(QDataStream&, SpecifiedGroupsDeleter&);

public:
    SpecifiedGroupsDeleter();
    
    SpecifiedGroupsDeleter(const MoleculeGroup &molgroup, const MGID &mgid);
    SpecifiedGroupsDeleter(const Sampler &sampler, const MGID &mgid);
    
    SpecifiedGroupsDeleter(const SpecifiedGroupsDeleter &other);
    
    ~SpecifiedGroupsDeleter();

    static const char* typeName();
    
    SpecifiedGroupsDeleter& operator=(const SpecifiedGroupsDeleter &other);
    
    bool operator==(const SpecifiedGroupsDeleter &other) const;
    bool operator!=(const SpecifiedGroupsDeleter &other) const;

    void setGenerator(const RanGenerator &generator);
    const RanGenerator& generator() const;
    
    const Sampler& sampler() const;
    const MoleculeGroup& group() const;
    
    void setSampler(const Sampler &sampler);
    void setSampler(const MoleculeGroup &molgroup);
    
    void setGroup(const MoleculeGroup &molgroup);
    
    const MGID& specifiedGroups() const;
    
    void setSpecifiedGroups(const MGID &mgid);

    tuple<Molecule,double> deleteFrom(System &system);

private:
    /** The sampler used to select molecules to be deleted */
    SamplerPtr smplr;
    
    /** The ID used to find the groups from which the molecules
        will be deleted */
    MGIdentifier mgid;
};

typedef SireBase::PropPtr<MolDeleter> MolDeleterPtr;

}

Q_DECLARE_METATYPE( SireMove::NullDeleter )
Q_DECLARE_METATYPE( SireMove::SpecifiedGroupsDeleter )
Q_DECLARE_METATYPE( SireMove::SystemWideDeleter )

SIRE_EXPOSE_CLASS( SireMove::MolDeleter )
SIRE_EXPOSE_CLASS( SireMove::NullDeleter )
SIRE_EXPOSE_CLASS( SireMove::SpecifiedGroupsDeleter )
SIRE_EXPOSE_CLASS( SireMove::SystemWideDeleter )

SIRE_EXPOSE_PROPERTY( SireMove::MolDeleterPtr, SireMove::MolDeleter )

SIRE_END_HEADER

#endif
