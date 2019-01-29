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

#ifndef SIREMOVE_INTERNALMOVE_H
#define SIREMOVE_INTERNALMOVE_H

#include "montecarlo.h"
#include "sampler.h"
#include "getpoint.h"

SIRE_BEGIN_HEADER

namespace SireMove
{
class InternalMove;
class DofID;
}

SIREMOVE_EXPORT QDataStream& operator<<(QDataStream&, const SireMove::InternalMove&);
SIREMOVE_EXPORT QDataStream& operator>>(QDataStream&, SireMove::InternalMove&);

namespace SireMol
{
class MoleculeGroup;
class PartialMolecule;
class AtomIdx;
class Connectivity;
class AtomID;
class BondID;
class AngleID;
class DihedralID;
}

namespace SireMove
{

class Sampler;
class ZMatrixCoords;

using SireMol::MoleculeGroup;
using SireMol::PartialMolecule;
using SireMol::Connectivity;
using SireMol::AtomIdx;
using SireMol::AtomID;
using SireMol::BondID;
using SireMol::AngleID;
using SireMol::DihedralID;

/** This class implements an intramolecular Monte Carlo move that uses 
    the move() method to perturb intramolecular degrees of freedom and 
    that may be applied to a random molecule (or part of a molecule)
    within a MoleculeGroup. It is based on the ZMatMove class.
    
    @author Julien Michel
*/
class SIREMOVE_EXPORT InternalMove 
            : public SireBase::ConcreteProperty<InternalMove,MonteCarlo>
{

friend QDataStream& ::operator<<(QDataStream&, const InternalMove&);
friend QDataStream& ::operator>>(QDataStream&, InternalMove&);

public:
    InternalMove(const PropertyMap &map = PropertyMap() );
    
    InternalMove(const MoleculeGroup &molgroup,
                 const PropertyMap &map = PropertyMap() );

    InternalMove(const Sampler &sampler, 
                 const PropertyMap &map = PropertyMap() );
    
    InternalMove(const InternalMove &other);
    
    ~InternalMove();
    
    InternalMove& operator=(const InternalMove &other);
    
    static const char* typeName();

    bool operator==(const InternalMove &other) const;
    bool operator!=(const InternalMove &other) const;

    QString toString() const;

    void setSampler(const Sampler &sampler);
    void setSampler(const MoleculeGroup &molgroup);

    const Sampler& sampler() const;
    const MoleculeGroup& moleculeGroup() const;
    
    const PropertyName& flexibilityProperty() const;

    void setFlexibilityProperty(const PropertyName &property);

    void setGenerator(const RanGenerator &rangenerator);

    void setCenterOfMolecule(const GetPoint &center_function);
    const GetPoint& centerOfMolecule() const;

    void move(System &system, int nmoves, bool record_stats=true);

protected:
    void _pvt_setTemperature(const SireUnits::Dimension::Temperature &temperature);

private:
    /** The sampler used to select random molecules for the move */
    SamplerPtr smplr;

    /** The function used to get the center of each molecule */
    GetPointPtr center_function;

    /** The name of the property that contains the flexibility*/
    PropertyName flexibility_property;
};

}

Q_DECLARE_METATYPE( SireMove::InternalMove )

SIRE_EXPOSE_CLASS( SireMove::InternalMove )

SIRE_END_HEADER

#endif
