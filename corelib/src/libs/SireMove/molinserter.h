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

#ifndef SIREMOVE_MOLINSERTER_H
#define SIREMOVE_MOLINSERTER_H

#include "SireBase/property.h"

#include "SireMaths/rangenerator.h"

#include "SireMol/mgidsandmaps.h"

#include <QStringList>

SIRE_BEGIN_HEADER

namespace SireMove
{
class MolInserter;
class NullInserter;

class UniformInserter;
}

SIREMOVE_EXPORT QDataStream& operator<<(QDataStream&, const SireMove::MolInserter&);
SIREMOVE_EXPORT QDataStream& operator>>(QDataStream&, SireMove::MolInserter&);

SIREMOVE_EXPORT QDataStream& operator<<(QDataStream&, const SireMove::NullInserter&);
SIREMOVE_EXPORT QDataStream& operator>>(QDataStream&, SireMove::NullInserter&);

SIREMOVE_EXPORT QDataStream& operator<<(QDataStream&, const SireMove::UniformInserter&);
SIREMOVE_EXPORT QDataStream& operator>>(QDataStream&, SireMove::UniformInserter&);

namespace SireMol
{
class Molecule;
class PartialMolecule;
}

namespace SireVol
{
class Space;
}

namespace SireSystem
{
class System;
}

namespace SireMove
{

using SireMaths::RanGenerator;
using SireMol::MGIDsAndMaps;

using SireMol::Molecule;
using SireMol::PartialMolecule;
using SireVol::Space;

using SireSystem::System;

/** This is the base class of all molecule inserters. These are
    manipulator classes that are used to insert (add) molecules
    to a system or molecule group(s) during a running simulation.
    e.g. This is useful for Grand Canonial or Gibbs Ensemble simulations.
    
    @author Christopher Woods
*/
class SIREMOVE_EXPORT MolInserter : public SireBase::Property
{

friend QDataStream& ::operator<<(QDataStream&, const MolInserter&);
friend QDataStream& ::operator>>(QDataStream&, MolInserter&);

public:
    MolInserter();
    MolInserter(const MolInserter &other);
    
    virtual ~MolInserter();
    
    static const char* typeName();
    
    virtual MolInserter* clone() const=0;
    
    virtual void setGenerator(const RanGenerator &generator);
    virtual void setGroups(const MGIDsAndMaps &mgids);
    
    const RanGenerator& generator() const;
    const MGIDsAndMaps& groups() const;
    
    /** Insert the entire molecule 'molecule' into the system 'system'
        using the space 'space', using the information about which
        groups to add the molecule to (and which properties to use)
        which are contained in this inserter). This returns the
        probability of where the molecule was inserted (normalised
        so that a probability of 1 is returned if the molecule
        was added at a uniformly random orientation and position) */
    virtual double insert(const Molecule &molecule,
                          System &system, const Space &space)=0;
    
    /** Insert the partial molecule 'molecule' into the system 'system'
        using the space 'space', using the information about which
        groups to add the molecule to (and which properties to use)
        which are contained in this inserter). This returns the
        probability of where the molecule was inserted (normalised
        so that a probability of 1 is returned if the molecule
        was added at a uniformly random orientation and position) */
    virtual double insert(const PartialMolecule &molecule,
                          System &system, const Space &space)=0;
    
    static const NullInserter& null();
    
protected:
    MolInserter& operator=(const MolInserter &other);
    
    bool operator==(const MolInserter &other) const;
    bool operator!=(const MolInserter &other) const;
    
    void rebuildCoordsProperties();
    
    const QStringList& coordsProperties() const;
    
private:
    /** The random number generator used to randomly position
        (and/or orientate) the molecule when it is inserted */
    RanGenerator rangen;
    
    /** The identities of which groups the molecule should be added
        to, together with the property maps to use for each group
        to find the correct properties */
    MGIDsAndMaps mgids;
    
    /** The list of all coordinates properties */
    QStringList coords_properties;
};

/** This is the null inserter - this doesn't insert anything 
    into anything
    
    @author Christopher Woods
*/
class SIREMOVE_EXPORT NullInserter 
            : public SireBase::ConcreteProperty<NullInserter,MolInserter>
{

friend QDataStream& ::operator<<(QDataStream&, const NullInserter&);
friend QDataStream& ::operator>>(QDataStream&, NullInserter&);

public:
    NullInserter();
    NullInserter(const NullInserter &other);
    
    ~NullInserter();
    
    static const char* typeName();
    
    NullInserter& operator=(const NullInserter &other);
    
    bool operator==(const NullInserter &other) const;
    bool operator!=(const NullInserter &other) const;
    
    double insert(const Molecule &molecule, System &system,
                  const Space &space);
                  
    double insert(const PartialMolecule &molecule, System &system,
                  const Space &space);
};

/** This inserter inserts a molecule to a random point in space, using
    a random orientation, chosen uniformly at random
    
    @author Christopher Woods
*/
class SIREMOVE_EXPORT UniformInserter
            : public SireBase::ConcreteProperty<UniformInserter,MolInserter>
{

friend QDataStream& ::operator<<(QDataStream&, const UniformInserter&);
friend QDataStream& ::operator>>(QDataStream&, UniformInserter&);

public:
    UniformInserter();
    
    UniformInserter(const MGIDsAndMaps &mgids);
    
    UniformInserter(const UniformInserter &other);
    
    ~UniformInserter();

    static const char* typeName();
    
    UniformInserter& operator=(const UniformInserter &other);
    
    bool operator==(const UniformInserter &other) const;
    bool operator!=(const UniformInserter &other) const;
    
    double insert(const Molecule &molecule, System &system,
                  const Space &space);
                  
    double insert(const PartialMolecule &molecule, System &system,
                  const Space &space);

private:
    template<class T>
    void uniform_insert(const T &molecule, System &system,
                        const Space &space) const;
};

typedef SireBase::PropPtr<MolInserter> MolInserterPtr;

}

Q_DECLARE_METATYPE( SireMove::NullInserter )
Q_DECLARE_METATYPE( SireMove::UniformInserter )

SIRE_EXPOSE_CLASS( SireMove::MolInserter )
SIRE_EXPOSE_CLASS( SireMove::NullInserter )
SIRE_EXPOSE_CLASS( SireMove::UniformInserter )

SIRE_EXPOSE_PROPERTY( SireMove::MolInserterPtr, SireMove::MolInserter )

SIRE_END_HEADER

#endif
