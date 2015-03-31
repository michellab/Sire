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

#include "sampler.h"
#include "uniformsampler.h"

#include "SireMol/molecule.h"
#include "SireMol/partialmolecule.h"

#include "SireSystem/system.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMove;
using namespace SireMol;
using namespace SireBase;
using namespace SireStream;

//////////
////////// Implementation of Sampler
//////////

static const RegisterMetaType<Sampler> r_sampler(MAGIC_ONLY,
                                                 "SireMove::Sampler");

/** Serialise to a binary datastream */
QDataStream SIREMOVE_EXPORT &operator<<(QDataStream &ds, 
                                        const Sampler &sampler)
{
    writeHeader(ds, r_sampler, 1);

    SharedDataStream sds(ds);
    sds << sampler.molgroup << sampler.rangen
        << static_cast<const Property&>(sampler);

    return ds;
}

/** Deserialise from a binary datastream */
QDataStream SIREMOVE_EXPORT &operator>>(QDataStream &ds, 
                                        Sampler &sampler)
{
    VersionID v = readHeader(ds, r_sampler);

    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> sampler.molgroup >> sampler.rangen
            >> static_cast<Property&>(sampler);
    }
    else
        throw version_error(v, "1", r_sampler, CODELOC);

    return ds;
}

/** Empty constructor */
Sampler::Sampler() : Property()
{}

/** Construct a sampler that picks molecules at random from the 
    molecule group 'molgroup' */
Sampler::Sampler(const MoleculeGroup &moleculegroup)
        : Property(), molgroup(moleculegroup)
{}

/** Copy constructor */
Sampler::Sampler(const Sampler &other) 
        : Property(other), molgroup(other.molgroup), rangen(other.rangen)
{}

/** Destructor */
Sampler::~Sampler()
{}

/** Copy assignment */
Sampler& Sampler::operator=(const Sampler &other)
{
    Property::operator=(other);

    molgroup = other.molgroup;
    rangen = other.rangen;

    return *this;
}

/** Comparison operator */
bool Sampler::operator==(const Sampler &other) const
{
    return molgroup == other.molgroup;
}

/** Comparison operator */
bool Sampler::operator!=(const Sampler &other) const
{
    return molgroup != other.molgroup;
}

/** Set the random number generator used by this sampler */
void Sampler::setGenerator(const RanGenerator &generator)
{
    rangen = generator;
}

/** Set the molecule group from which random molecules will be sampled */    
void Sampler::setGroup(const MoleculeGroup &moleculegroup)
{
    molgroup = moleculegroup;
}

/** Update this sampler so that it matches the state of the molecules 
    in the System 'system' */
void Sampler::updateFrom(const System &system)
{
    if (system.contains( molgroup.read().number() ))
    {
        const MoleculeGroup &new_group = system[molgroup.read().number()];
        
        if (new_group.version() != molgroup.read().version())
            //the molecule group has changed - update it
            this->setGroup( system[molgroup.read().number()] );
    }
}
