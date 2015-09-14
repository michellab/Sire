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

#ifndef SIREMOVE_SAMPLER_H
#define SIREMOVE_SAMPLER_H

#include <boost/tuple/tuple.hpp>

#include "SireMol/moleculegroup.h"

#include "SireBase/property.h"

#include "SireMaths/rangenerator.h"

SIRE_BEGIN_HEADER

namespace SireMove
{
class Sampler;
}

QDataStream& operator<<(QDataStream&, const SireMove::Sampler&);
QDataStream& operator>>(QDataStream&, SireMove::Sampler&);

namespace SireMol
{
class Molecule;
class PartialMolecule;
}

namespace SireSystem
{
class System;
}

namespace SireMove
{

using boost::tuple;

using SireBase::Property;

using SireMaths::RanGenerator;

using SireMol::Molecule;
using SireMol::PartialMolecule;
using SireMol::MolGroupPtr;
using SireMol::MoleculeGroup;

using SireSystem::System;

class UniformSampler;

/** This is the base class of all Samplers. A Sampler is used
    to pick a random molecule from a MoleculeGroup

    @author Christopher Woods
*/
class SIREMOVE_EXPORT Sampler : public Property
{

friend QDataStream& ::operator<<(QDataStream&, const Sampler&);
friend QDataStream& ::operator>>(QDataStream&, Sampler&);

public:
    Sampler();
    Sampler(const MoleculeGroup &molgroup);

    Sampler(const Sampler &other);

    virtual ~Sampler();

    virtual Sampler* clone() const=0;

    static const char* typeName()
    {
        return "SireMove::Sampler";
    }

    void setGenerator(const RanGenerator &generator);
    const RanGenerator& generator() const;

    const MoleculeGroup& group() const;
    
    virtual void setGroup(const MoleculeGroup &molgroup);

    virtual void updateFrom(const System &system);

    virtual tuple<PartialMolecule,double> sample() const=0;
    virtual tuple<Molecule,double> sampleMolecule() const=0;

    virtual double probabilityOf(const PartialMolecule &molecule) const=0;
    virtual double probabilityOfMolecule(const Molecule &molecule) const=0;

    virtual bool isBiased() const;

    static const UniformSampler& null();

protected:
    Sampler& operator=(const Sampler &other);
    
    bool operator==(const Sampler &other) const;
    bool operator!=(const Sampler &other) const;

private:
    /** The molecule group from which molecules are sampled */
    MolGroupPtr molgroup;

    /** The random number generator used by the sampler */
    RanGenerator rangen;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Return the molecule group from which a molecule will be sampled */
inline const MoleculeGroup& Sampler::group() const
{
    return molgroup;
}

/** Internal function used to return a reference to the random
    number generator used by this sampler */
inline const RanGenerator& Sampler::generator() const
{
    return rangen;
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

typedef SireBase::PropPtr<Sampler> SamplerPtr;

}

SIRE_EXPOSE_CLASS( SireMove::Sampler )

SIRE_EXPOSE_PROPERTY( SireMove::SamplerPtr, SireMove::Sampler )

SIRE_END_HEADER

#include "uniformsampler.h"

#endif
