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

#ifndef SQUIRE_QMPOTENTIAL_H
#define SQUIRE_QMPOTENTIAL_H

#include "SireMol/element.h"

#include "SireBase/property.h"
#include "SireBase/properties.h"

#include "SireVol/space.h"

#include "SireFF/ffcomponent.h"

#include "SireFF/detail/ffmolecules3d.h"
#include "SireFF/detail/atomicparameters3d.hpp"

#include "SireUnits/dimensions.h"

SIRE_BEGIN_HEADER

namespace Squire
{
class QMPotential;
class QMProgram;
class QMComponent;
}

QDataStream& operator<<(QDataStream&, const Squire::QMPotential&);
QDataStream& operator>>(QDataStream&, Squire::QMPotential&);

QDataStream& operator<<(QDataStream&, const Squire::QMProgram&);
QDataStream& operator>>(QDataStream&, Squire::QMProgram&);

QDataStream& operator<<(QDataStream&, const Squire::QMComponent&);
QDataStream& operator>>(QDataStream&, Squire::QMComponent&);

namespace SireFF
{
class ForceTable;
class FieldTable;
class PotentialTable;
class Probe;
}

namespace SireMM
{
class CoulombProbe;
}

namespace SireCAS
{
class Symbol;
}

namespace Squire
{

using SireBase::PropertyMap;
using SireBase::Property;
using SireBase::Property;
using SireBase::Properties;

using SireMol::PartialMolecule;
using SireMol::MoleculeGroup;

using SireCAS::Symbol;

using SireVol::Space;
using SireVol::SpacePtr;

using SireFF::ForceTable;
using SireFF::FieldTable;
using SireFF::PotentialTable;
using SireFF::Probe;

class QMComponent;
class QMProgram;
class LatticeCharges;

typedef SireFF::ComponentEnergy<QMComponent> QMEnergy;

/** This class represents a QM energy */
class SQUIRE_EXPORT QMComponent : public SireFF::FFComponent
{
public:
    QMComponent(const SireFF::FFName &ffname = SireFF::FFName());
    QMComponent(const SireCAS::Symbol &symbol);
    
    QMComponent(const QMComponent &other);
    
    ~QMComponent();
    
    static const char* typeName()
    {
        return "Squire::QMComponent";
    }
    
    const char* what() const
    {
        return QMComponent::typeName();
    }
    
    QMComponent* clone() const
    {
        return new QMComponent(*this);
    }
    
    const QMComponent& total() const
    {
        return *this;
    }

    void setEnergy(SireFF::FF &ff, const QMEnergy &qmnrg) const;
    void changeEnergy(SireFF::FF &ff, const QMEnergy &qmnrg) const;
    
    SireCAS::Symbols symbols() const
    {
        return *this;
    }
};

/** This class provides the default name of the 
    property that contains the element parameters */
class SQUIRE_EXPORT ElementParameterName
{
public:
    ElementParameterName()
    {}
    
    ~ElementParameterName()
    {}
    
    const QString& element() const
    {
        return element_param;
    }
    
private:
    static QString element_param;
};

/** This class provides the default names of the element
    and coordinates properties that contain these parameters */
class SQUIRE_EXPORT ElementParameterName3D
          : public ElementParameterName, public SireFF::detail::Coords3DParameterName
{
public:
    ElementParameterName3D() : ElementParameterName(),
                               SireFF::detail::Coords3DParameterName()
    {}
    
    ~ElementParameterName3D()
    {}
};

typedef SireBase::PropPtr<QMProgram> QMProgPtr;

/** This is a QM potential. This provides the classes and code necessary
    to calculate QM energies and forces
    
    @author Christopher Woods
*/
class SQUIRE_EXPORT QMPotential
{

friend QDataStream& ::operator<<(QDataStream&, const QMPotential&);
friend QDataStream& ::operator>>(QDataStream&, QMPotential&);

public:
    
    typedef QMEnergy Energy;
    typedef Energy::Components Components;
    typedef ElementParameterName3D ParameterNames;
    
    typedef SireMM::CoulombProbe Probe;
    
    typedef SireMol::Element Parameter;
    typedef SireFF::detail::AtomicParameters3D<Parameter> Parameters;
    
    typedef SireFF::detail::FFMolecule3D<QMPotential> Molecule;
    typedef SireFF::detail::FFMolecules3D<QMPotential> Molecules;
    typedef SireFF::detail::ChangedMolecule<Molecule> ChangedMolecule;

    QMPotential();
    
    QMPotential(const QMPotential &other);
    
    virtual ~QMPotential();

    static ParameterNames parameters()
    {
        return ParameterNames();
    }

    static const char* typeName()
    {
        return "Squire::QMPotential";
    }
    
    virtual const char* what() const
    {
        return QMPotential::typeName();
    }

    QMPotential& operator=(const QMPotential &other);

    bool setProperty(const QString &name, const Property &value);
    const Property& property(const QString &name) const;
    bool containsProperty(const QString &name) const;
    const Properties& properties() const;

    bool setSpace(const Space &space);
    bool setQuantumProgram(const QMProgram &program);
    bool setZeroEnergy(SireUnits::Dimension::MolarEnergy zero_energy);
    
    const Space& space() const;
    const QMProgram& quantumProgram() const;
    SireUnits::Dimension::MolarEnergy zeroEnergy() const;

    QMPotential::Parameters 
    getParameters(const PartialMolecule &mol,
                  const PropertyMap &map = PropertyMap());
                  
    QMPotential::Parameters
    updateParameters(const QMPotential::Parameters &old_params,
                     const PartialMolecule &old_mol,
                     const PartialMolecule &new_mol,
                     const PropertyMap &map = PropertyMap());
                     
    QMPotential::Parameters
    updateParameters(const QMPotential::Parameters &old_params,
                     const PartialMolecule &old_molecule,
                     const PartialMolecule &new_molecule,
                     const PropertyMap &old_map, const PropertyMap &new_map);
                     
    QMPotential::Molecule
    parameterise(const PartialMolecule &molecule,
                 const PropertyMap &map = PropertyMap());
    
    QMPotential::Molecules 
    parameterise(const MoleculeGroup &molecules,
                 const PropertyMap &map = PropertyMap());

    void calculateForce(const Molecules &molecules, 
                        ForceTable &forcetable, 
                        double scale_force=1) const;
                        
    void calculateForce(const Molecules &molecules, 
                        ForceTable &forcetable,
                        const Symbol &symbol, 
                        const Components &components,
                        double scale_force=1) const;

    void calculatePotential(const Molecules &molecules,
                            PotentialTable &pottable,
                            const SireFF::Probe &probe,
                            double scale_potential=1) const;

    void calculatePotential(const Molecules &molecules,
                            PotentialTable &pottable,
                            const SireFF::Probe &probe,
                            const Symbol &symbol,
                            const Components &components,
                            double scale_potential=1) const;

    void calculateField(const Molecules &molecules,
                        FieldTable &fieldtable,
                        const SireFF::Probe &probe,
                        double scale_field=1) const;

    void calculateField(const Molecules &molecules,
                        FieldTable &fieldtable,
                        const SireFF::Probe &probe,
                        const Symbol &symbol,
                        const Components &components,
                        double scale_field=1) const;
    
    void calculateEnergy(const Molecules &molecules, Energy &nrg,
                         double scale_energy=1) const;

    QString energyCommandFile(const Molecules &molecules) const;

    QString forceCommandFile(const Molecules &molecules,
                             const ForceTable &forcetable) const;
    
    QString fieldCommandFile(const Molecules &molecules,
                             const FieldTable &fieldtable,
                             const SireFF::Probe &probe) const;

    QString potentialCommandFile(const Molecules &molecules,
                                 const PotentialTable &pottable,
                                 const SireFF::Probe &probe) const;

protected:
    virtual void changedPotential()=0;

    Molecules mapIntoSpace(const Molecules &molecules) const;
    
private:
    /** The properties that define this potential */
    Properties props;
    
    /** The space in which the molecules exist */
    SpacePtr spce;
    
    /** The QM program that is used to calculate the energy */
    QMProgPtr qmprog;
    
    /** The absolute value of the energy that corresponds
        to 'zero' - this is necessary to shift the QM energy
        into the same range as the MM energy. This is a constant
        that is added on to every QM energy evaluated */
    double zero_energy;
};

} // end of namespace Squire

Q_DECLARE_METATYPE( Squire::QMComponent )

SIRE_EXPOSE_CLASS( Squire::QMComponent )

SIRE_EXPOSE_CLASS( Squire::ElementParameterName )
SIRE_EXPOSE_CLASS( Squire::ElementParameterName3D )

SIRE_END_HEADER

#endif
