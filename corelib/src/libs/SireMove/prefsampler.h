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

#ifndef SIREMOVE_PREFSAMPLER_H
#define SIREMOVE_PREFSAMPLER_H

#include "sampler.h"

#include "SireMol/partialmolecule.h"

#include "SireBase/propertymap.h"

#include "SireVol/space.h"

#include "SireCAS/expression.h"

#include "SireMaths/vector.h"

#include "SireUnits/dimensions.h"

SIRE_BEGIN_HEADER

namespace SireMove
{
class PrefSampler;
}

SIREMOVE_EXPORT QDataStream& operator<<(QDataStream&, const SireMove::PrefSampler&);
SIREMOVE_EXPORT QDataStream& operator>>(QDataStream&, SireMove::PrefSampler&);

namespace SireCAS
{
class Symbol;
}

namespace SireMove
{

using SireMol::MoleculeView;
using SireMol::PartialMolecule;

using SireVol::SpacePtr;

using SireMaths::Vector;

using SireBase::PropertyName;

/** This is a sampler that is used to select molecules
    using the distance-based preferential sampling algorithm.
    
    This sampler can be used to bias the choice of molecules
    such that those closest to a target molecule are chosen
    with a higher probability than those that are further
    away
    
    @author Christopher Woods
*/
class SIREMOVE_EXPORT PrefSampler
          : public SireBase::ConcreteProperty<PrefSampler,Sampler>
{

friend QDataStream& ::operator<<(QDataStream&, const PrefSampler&);
friend QDataStream& ::operator>>(QDataStream&, PrefSampler&);

public:
    PrefSampler();
    PrefSampler(SireUnits::Dimension::Area k);
    PrefSampler(const SireCAS::Expression &f);
    PrefSampler(const SireCAS::Expression &f, SireUnits::Dimension::Area k);
    
    PrefSampler(const Vector &point);
    PrefSampler(const Vector &point, SireUnits::Dimension::Area k);
    PrefSampler(const Vector &point, const SireCAS::Expression &f);
    PrefSampler(const Vector &point, const SireCAS::Expression &f,
                SireUnits::Dimension::Area k);

    PrefSampler(const Vector &point, const MoleculeGroup &molgroup);
    PrefSampler(const Vector &point, const MoleculeGroup &molgroup,
                SireUnits::Dimension::Area k);
    PrefSampler(const Vector &point, const MoleculeGroup &molgroup,
                const SireCAS::Expression &f);
    PrefSampler(const Vector &point, const MoleculeGroup &molgroup,
                const SireCAS::Expression &f, SireUnits::Dimension::Area k);
    
    PrefSampler(const MoleculeView &molview);
    PrefSampler(const MoleculeView &molview, SireUnits::Dimension::Area k);
    PrefSampler(const MoleculeView &molview, const SireCAS::Expression &f);
    PrefSampler(const MoleculeView &molview, const SireCAS::Expression &f,
                SireUnits::Dimension::Area k);
    
    PrefSampler(const MoleculeView &molview, const MoleculeGroup &molgroup);
    PrefSampler(const MoleculeView &molview, const MoleculeGroup &molgroup,
                SireUnits::Dimension::Area k);
    PrefSampler(const MoleculeView &molview, const MoleculeGroup &molgroup,
                const SireCAS::Expression &f);
    PrefSampler(const MoleculeView &molview, const MoleculeGroup &molgroup,
                const SireCAS::Expression &f, SireUnits::Dimension::Area k);
    
    PrefSampler(const PrefSampler &other);
    
    ~PrefSampler();
    
    PrefSampler& operator=(const PrefSampler &other);
    
    static const char* typeName();
    
    bool operator==(const PrefSampler &other) const;
    bool operator!=(const PrefSampler &other) const;
    
    static SireCAS::Symbol r();
    static SireCAS::Symbol k();
    
    void setGroup(const MoleculeGroup &molgroup);

    void updateFrom(const System &system);
    
    tuple<PartialMolecule,double> sample() const;
    tuple<Molecule,double> sampleMolecule() const;

    double probabilityOf(const PartialMolecule &molecule) const;
    double probabilityOfMolecule(const Molecule &molecule) const;

    void setFocalMolecule(const MoleculeView &molview);
    void setFocalPoint(const Vector &point);
    
    void setCoordinatesProperty(const PropertyName &coords_property);
    void setSpaceProperty(const PropertyName &space_property);
    
    void setSamplingConstant(SireUnits::Dimension::Area k);
    void setBiasingFunction(const SireCAS::Expression &f);

    SireCAS::Expression biasingFunction() const;
    SireUnits::Dimension::Area samplingConstant() const;
    
    bool usingFocalMolecule() const;
    bool usingFocalPoint() const;
    
    const Vector& focalPoint() const;
    const PartialMolecule& focalMolecule() const;

    const PropertyName& coordinatesProperty() const;
    const PropertyName& spaceProperty() const;

    bool isBiased() const;

private:
    void updateWeights(const MoleculeGroup &new_group);
    void recalculateWeights();

    /** The view of the molecule, the center of which is used
        as the focal point for the preferential sampling algorithm */
    PartialMolecule focal_molecule;
    
    /** The actual focal point of the preferential sampling algorithm */
    Vector focal_point;
    
    /** The coordinates property used to find the coordinates
        of the molecules */
    PropertyName coords_property;
    
    /** The property used to find the system space */
    PropertyName space_property;
    
    /** The preferential sampling expression */
    SireCAS::Expression sampling_expression;
    
    /** The preferential sampling constant */
    double sampling_constant;

    /** The sum of all of the weights */
    double sum_of_weights;

    /** The maximum weight */
    double max_weight;

    /** The current weights for all of the molecules - the
        index matches the viewAt() index of MoleculeGroup */
    QVector<double> molweights;

    /** The current space that is used to calculate distances */
    SpacePtr current_space;
    
    /** Whether or not this is dirty (requires a complete recalculation
        of the weights) */
    bool is_dirty;
};

}

Q_DECLARE_METATYPE( SireMove::PrefSampler )

SIRE_EXPOSE_CLASS( SireMove::PrefSampler )

SIRE_END_HEADER

#endif
