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

#ifndef SIREMM_COULOMBPOTENTIAL_H
#define SIREMM_COULOMBPOTENTIAL_H

#include "SireBase/properties.h"
#include "SireBase/propertymap.h"
#include "SireBase/pairmatrix.hpp"
#include "SireBase/packedarray2d.hpp"

#include "SireVol/space.h"

#include "SireMol/atomproperty.hpp"
#include "SireMol/atomcharges.h"

#include "SireUnits/dimensions.h"

#include "cljcomponent.h"
#include "cljnbpairs.h"
#include "cljprobe.h"

#include "detail/intrascaledatomicparameters.hpp"

#include "switchingfunction.h"

#include "SireFF/energytable.h"
#include "SireFF/forcetable.h"
#include "SireFF/fieldtable.h"
#include "SireFF/potentialtable.h"
#include "SireFF/detail/ffmolecules3d.h"

SIRE_BEGIN_HEADER

namespace SireMM
{
class InterCoulombPotential;
class IntraCoulombPotential;
class CoulombPotential;

template<class CoulPot>
class CoulombPotentialInterface;

namespace detail{ class CoulombParameter; }
}

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::InterCoulombPotential&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::InterCoulombPotential&);

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::IntraCoulombPotential&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::IntraCoulombPotential&);

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::CoulombPotential&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::CoulombPotential&);

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::detail::CoulombParameter&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::detail::CoulombParameter&);

template<class CoulPot>
QDataStream& operator<<(QDataStream&, const SireMM::CoulombPotentialInterface<CoulPot>&);
template<class CoulPot>
QDataStream& operator>>(QDataStream&, SireMM::CoulombPotentialInterface<CoulPot>&);

namespace SireMol
{
class PartialMolecule;
class MoleculeGroup;
}

namespace SireMM
{

using SireBase::Properties;
using SireBase::Property;
using SireBase::PropertyMap;
using SireBase::PropertyName;
using SireBase::Property;

using SireCAS::Symbol;

using SireMaths::Vector;

using SireVol::Space;
using SireVol::SpacePtr;
using SireVol::CoordGroup;

using SireMol::PartialMolecule;
using SireMol::MoleculeGroup;

using SireFF::MolEnergyTable;
using SireFF::MolForceTable;
using SireFF::MolFieldTable;
using SireFF::MolPotentialTable;
using SireFF::GridFieldTable;
using SireFF::GridPotentialTable;

using SireMol::AtomCharges;

/** This class provides the default name of the 
    property that contains the charge parameters */
class SIREMM_EXPORT ChargeParameterName
{
public:
    ChargeParameterName()
    {}
    
    ~ChargeParameterName()
    {}
    
    const QString& charge() const
    {
        return chg_param;
    }
    
private:
    static QString chg_param;
};

/** This class provides the default name of the properties
    that contain the charge, LJ and 3D coordinates properties */
class SIREMM_EXPORT ChargeParameterName3D : public ChargeParameterName,
                                            public SireFF::detail::Coords3DParameterName
{
public:
    ChargeParameterName3D() : ChargeParameterName(),
                              SireFF::detail::Coords3DParameterName()
    {}
    
    ~ChargeParameterName3D()
    {}
};

/** This class provides the default name of the properties 
    that contain the charge and intramolecular NB scale parameters and
    3D coordinates properties */
class SIREMM_EXPORT ScaledChargeParameterNames3D : public ChargeParameterName3D,
                                                   public detail::IntraScaleParameterName
{
public:
    ScaledChargeParameterNames3D() : ChargeParameterName3D(),
                                     detail::IntraScaleParameterName()
    {}
    
    ~ScaledChargeParameterNames3D()
    {}
};

namespace detail
{

/** This class holds the charge parameter used by both the Inter- and Intra-
    CoulombPotentials. It is just the charge of the atom in internal units
    (charge divided by sqrt(4 pi epsilon_0))
    
    @author Christopher Woods
*/
class SIREMM_EXPORT ChargeParameter
{
public:
    ChargeParameter(double charge=0) : reduced_charge(charge)
    {}
    
    ChargeParameter(const ChargeParameter &other) : reduced_charge(other.reduced_charge)
    {}
    
    ~ChargeParameter()
    {}
    
    bool operator==(const ChargeParameter &other) const
    {
        return reduced_charge == other.reduced_charge;
    }
    
    bool operator!=(const ChargeParameter &other) const
    {
        return reduced_charge != other.reduced_charge;
    }
    
    double reduced_charge;
};

} // end of namespace detail

/** This is the common base class of InterCoulombPotential and IntraCoulombPotential

    @author Christopher Woods
*/
class SIREMM_EXPORT CoulombPotential
{

friend QDataStream& ::operator<<(QDataStream&, const CoulombPotential&);
friend QDataStream& ::operator>>(QDataStream&, CoulombPotential&);

public:
    virtual ~CoulombPotential();

    const Properties& properties() const;
    const Property& property(const QString &name) const;
    bool containsProperty(const QString &name) const;
    
    bool setProperty(const QString &name, const Property &value);

    bool setSpace(const Space &new_space);
    bool setSwitchingFunction(const SwitchingFunction &new_switchfunc);
    bool setShiftElectrostatics(bool switchelectro);
    
    const Space& space() const;
    const SwitchingFunction& switchingFunction() const;
    bool shiftElectrostatics() const;

protected:
    CoulombPotential();
    CoulombPotential(const CoulombPotential &other);
    
    CoulombPotential& operator=(const CoulombPotential &other);

    void startEvaluation();
    void finishedEvaluation();
    
    virtual void changedPotential()=0;

    /** The current values of the properties of this functional */
    Properties props;
    
    /** The space in which this functional operates */
    SpacePtr spce;
    
    /** The nonbonded switching function */
    SwitchFuncPtr switchfunc;
    
    /** Whether or not electrostatic potential shifting is used
        (this shifts the entire electrostatic potential so that it
        is zero at the cutoff distance) */
    bool use_electrostatic_shifting;
};

/** This class provides all of the functions and containers  
    necessary to provide an interface to calculate the
    intermolecular interatomic potentials using a Coulomb potential 
    
    This is a 3D potential class, namely it requires that
    the atoms possess 3D coordinates, thereby allowing this
    potential to also be used to calculate 3D forces on the atoms.

    This potential has the following properties (parameters);
    
    (1) space               : This is the 3D space in which the molecules exist
    (2) switchingFunction   : This is the switching function used to scale the 
                              energies / forces to zero at the cutoff
    (3) shiftElectrostatics : This is a boolean - if it is true then the 
                              group-group electrostatic interactions are scaled 
                              so that they are zero at the cutoff
                              
    The molecules used with this potential must contain the following 
    properties (defined in parameters(), e.g. parameters().coordinates())
    
    (1) .coordinates()      : These are the 3D coordinates of each atom, must be
                              a property of type AtomCoords
    (2) .charge()           : These are the partial charges of each atom, must
                              be a property of type AtomCharges
    
    @author Christopher Woods
*/
class SIREMM_EXPORT InterCoulombPotential : public CoulombPotential
{

friend QDataStream& ::operator<<(QDataStream&, const InterCoulombPotential&);
friend QDataStream& ::operator>>(QDataStream&, InterCoulombPotential&);

public:

    typedef CoulombEnergy Energy;
    typedef Energy::Components Components;

    typedef ChargeParameterName3D ParameterNames;

    typedef detail::ChargeParameter Parameter;
    typedef SireFF::detail::AtomicParameters3D<Parameter> Parameters;
    
    typedef SireBase::PairMatrix<double> EnergyWorkspace;
    typedef SireBase::PairMatrix<SireMaths::DistVector> ForceWorkspace;
    typedef SireBase::PairMatrix<SireMaths::DistVector> FieldWorkspace;
    typedef SireBase::PairMatrix<double> PotentialWorkspace;

    typedef CoulombProbe Probe;

    typedef SireFF::detail::FFMolecule3D<InterCoulombPotential> Molecule;
    typedef SireFF::detail::FFMolecules3D<InterCoulombPotential> Molecules;

    typedef SireFF::detail::ChangedMolecule<Molecule> ChangedMolecule;

    InterCoulombPotential();
    
    InterCoulombPotential(const InterCoulombPotential &other);
    
    ~InterCoulombPotential();
    
    InterCoulombPotential& operator=(const InterCoulombPotential &other);

    static const char* typeName()
    {
        return "SireMM::InterCoulombPotential";
    }

    const char* what() const
    {
        return InterCoulombPotential::typeName();
    }
    
    static ParameterNames parameters()
    {
        return ParameterNames();
    }
    
    InterCoulombPotential::Parameters 
    getParameters(const PartialMolecule &molecule,
                  const PropertyMap &map = PropertyMap());
    
    InterCoulombPotential::Parameters
    updateParameters(const InterCoulombPotential::Parameters &old_params,
                     const PartialMolecule &old_molecule,
                     const PartialMolecule &new_molecule,
                     const PropertyMap &map = PropertyMap());
                     
    InterCoulombPotential::Parameters
    updateParameters(const InterCoulombPotential::Parameters &old_params,
                     const PartialMolecule &old_molecule,
                     const PartialMolecule &new_molecule,
                     const PropertyMap &old_map, const PropertyMap &new_map);
    
    InterCoulombPotential::Molecule
    parameterise(const PartialMolecule &molecule,
                 const PropertyMap &map = PropertyMap());
    
    InterCoulombPotential::Molecules 
    parameterise(const MoleculeGroup &molecules,
                 const PropertyMap &map = PropertyMap());

    void calculateEnergy(const InterCoulombPotential::Molecule &mol0, 
                         const InterCoulombPotential::Molecule &mol1,
                         InterCoulombPotential::Energy &energy, 
                         InterCoulombPotential::EnergyWorkspace &workspace,
                         double scale_energy=1) const;

    void calculateEnergy(const InterCoulombPotential::Molecule &mol0, 
                         const InterCoulombPotential::Molecule &mol1,
                         MolEnergyTable &energies0, 
                         InterCoulombPotential::EnergyWorkspace &workspace,
                         double scale_energy=1) const;

    void calculateForce(const InterCoulombPotential::Molecule &mol0, 
                        const InterCoulombPotential::Molecule &mol1,
                        MolForceTable &forces0,
                        InterCoulombPotential::ForceWorkspace &workspace,
                        double scale_force=1) const;

    void calculateForce(const InterCoulombPotential::Molecule &mol0,
                        const InterCoulombPotential::Molecule &mol1,
                        MolForceTable &forces0,
                        const Symbol &symbol,
                        const Components &components,
                        InterCoulombPotential::ForceWorkspace &workspace,
                        double scale_force=1) const;

    void calculateCoulombForce(const InterCoulombPotential::Molecule &mol0, 
                               const InterCoulombPotential::Molecule &mol1,
                               MolForceTable &forces0,
                               InterCoulombPotential::ForceWorkspace &workspace,
                               double scale_force=1) const;

    void calculateField(const InterCoulombPotential::Molecule &mol0, 
                        const InterCoulombPotential::Molecule &mol1,
                        const CoulombProbe &probe,
                        MolFieldTable &forces0,
                        InterCoulombPotential::FieldWorkspace &workspace,
                        double scale_field=1) const;

    void calculateField(const InterCoulombPotential::Molecule &mol0,
                        const InterCoulombPotential::Molecule &mol1,
                        const CoulombProbe &probe,
                        MolFieldTable &forces0,
                        const Symbol &symbol,
                        const Components &components,
                        InterCoulombPotential::FieldWorkspace &workspace,
                        double scale_field=1) const;

    void calculateField(const InterCoulombPotential::Molecule &mol0,
                        const CoulombProbe &probe,
                        GridFieldTable &fields,
                        InterCoulombPotential::FieldWorkspace &workspace,
                        double scale_field=1) const;

    void calculateField(const InterCoulombPotential::Molecule &mol0,
                        const CoulombProbe &probe,
                        GridFieldTable &fields,
                        const Symbol &symbol,
                        const Components &components,
                        InterCoulombPotential::FieldWorkspace &workspace,
                        double scale_field=1) const;

    void calculateCoulombField(const InterCoulombPotential::Molecule &mol0, 
                               const InterCoulombPotential::Molecule &mol1,
                               const CoulombProbe &probe,
                               MolFieldTable &fields0,
                               InterCoulombPotential::FieldWorkspace &workspace,
                               double scale_field=1) const;

    void calculateCoulombField(const InterCoulombPotential::Molecule &mol0, 
                               const CoulombProbe &probe,
                               GridFieldTable &fields,
                               InterCoulombPotential::FieldWorkspace &workspace,
                               double scale_field=1) const;

    void calculatePotential(const InterCoulombPotential::Molecule &mol0, 
                            const InterCoulombPotential::Molecule &mol1,
                            const CoulombProbe &probe,
                            MolPotentialTable &pots0,
                            InterCoulombPotential::PotentialWorkspace &workspace,
                            double scale_potential=1) const;

    void calculatePotential(const InterCoulombPotential::Molecule &mol0,
                            const InterCoulombPotential::Molecule &mol1,
                            const CoulombProbe &probe,
                            MolPotentialTable &pots0,
                            const Symbol &symbol,
                            const Components &components,
                            InterCoulombPotential::PotentialWorkspace &workspace,
                            double scale_potential=1) const;

    void calculatePotential(const InterCoulombPotential::Molecule &mol0,
                            const CoulombProbe &probe,
                            GridPotentialTable &pots,
                            InterCoulombPotential::PotentialWorkspace &workspace,
                            double scale_potential=1) const;

    void calculatePotential(const InterCoulombPotential::Molecule &mol0,
                            const CoulombProbe &probe,
                            GridPotentialTable &pots,
                            const Symbol &symbol,
                            const Components &components,
                            InterCoulombPotential::PotentialWorkspace &workspace,
                            double scale_potential=1) const;

    void calculateCoulombPotential(const InterCoulombPotential::Molecule &mol0, 
                                   const InterCoulombPotential::Molecule &mol1,
                                   const CoulombProbe &probe,
                                   MolPotentialTable &pots0,
                                   InterCoulombPotential::PotentialWorkspace &workspace,
                                   double scale_potential=1) const;

    void calculateCoulombPotential(const InterCoulombPotential::Molecule &mol0, 
                                   const CoulombProbe &probe,
                                   GridPotentialTable &pots,
                                   InterCoulombPotential::PotentialWorkspace &workspace,
                                   double scale_potential=1) const;

private:
    double totalCharge(const InterCoulombPotential::Parameters::Array &params) const;

    void throwMissingForceComponent(const Symbol &symbol,
                                    const Components &components) const;
    void throwMissingFieldComponent(const Symbol &symbol,
                                    const Components &components) const;
    void throwMissingPotentialComponent(const Symbol &symbol,
                                        const Components &components) const;

    void _pvt_calculateEnergy(const InterCoulombPotential::Molecule &mol0, 
                              const InterCoulombPotential::Molecule &mol1,
                              InterCoulombPotential::Energy &energy, 
                              InterCoulombPotential::EnergyWorkspace &workspace,
                              double scale_energy) const;

    void _pvt_calculateCoulombForce(const InterCoulombPotential::Molecule &mol0, 
                                    const InterCoulombPotential::Molecule &mol1,
                                    MolForceTable &forces0, 
                                    InterCoulombPotential::ForceWorkspace &workspace,
                                    double scale_force) const;

    void _pvt_calculateCoulombField(const InterCoulombPotential::Molecule &mol0, 
                                    const InterCoulombPotential::Molecule &mol1,
                                    const CoulombProbe &probe,
                                    MolFieldTable &forces0, 
                                    InterCoulombPotential::ForceWorkspace &workspace,
                                    double scale_force) const;

    void _pvt_calculateCoulombField(const InterCoulombPotential::Molecule &mol, 
                                    const CoulombProbe &probe,
                                    GridFieldTable &forces0, 
                                    InterCoulombPotential::ForceWorkspace &workspace,
                                    double scale_force) const;

    void _pvt_calculateCoulombPotential(
                                    const InterCoulombPotential::Molecule &mol0, 
                                    const InterCoulombPotential::Molecule &mol1,
                                    const CoulombProbe &probe,
                                    MolPotentialTable &pots0, 
                                    InterCoulombPotential::PotentialWorkspace &workspace,
                                    double scale_potential) const;

    void _pvt_calculateCoulombPotential(
                                    const InterCoulombPotential::Molecule &mol, 
                                    const CoulombProbe &probe,
                                    GridPotentialTable &pots, 
                                    InterCoulombPotential::PotentialWorkspace &workspace,
                                    double scale_potential) const;
};

/** This class provides all of the functions and containers  
    necessary to provide an interface to calculate the
    intramolecular interatomic potentials using a Coulomb potential
    
    This is a 3D potential class, namely it requires that
    the atoms possess 3D coordinates, thereby allowing this
    potential to also be used to calculate 3D forces on the atoms.

    This potential has the following properties (parameters);
    
    (1) space               : This is the 3D space in which the molecules exist
    (2) switchingFunction   : This is the switching function used to scale the 
                              energies / forces to zero at the cutoff
    (3) shiftElectrostatics : This is a boolean - if it is true then the 
                              group-group electrostatic interactions are scaled 
                              so that they are zero at the cutoff
                              
    The molecules used with this potential must contain the following 
    properties (defined in parameters(), e.g. parameters().coordinates())
    
    (1) .coordinates()       : These are the 3D coordinates of each atom, must be
                               a property of type AtomCoords
    (2) .charge()            : These are the partial charges of each atom, must
                               be a property of type AtomCharges
    (3) .intraScaleFactors() : These are the intramolecular atomic scaling factors,
                               used to scale the intramolecular coulomb
                               energies, must be a property of type CoulombNBPairs 
                               or CLJNBPairs
    
    @author Christopher Woods
*/
class SIREMM_EXPORT IntraCoulombPotential : public CoulombPotential
{

friend QDataStream& ::operator<<(QDataStream&, const IntraCoulombPotential&);
friend QDataStream& ::operator>>(QDataStream&, IntraCoulombPotential&);

public:
    typedef CoulombEnergy Energy;
    typedef Energy::Components Components;
    
    typedef ScaledChargeParameterNames3D ParameterNames;

    typedef detail::ChargeParameter Parameter;
    
    typedef detail::IntraScaledAtomicParameters<
                  SireFF::detail::AtomicParameters3D<Parameter>,
                  detail::IntraScaledParameters<CoulombNBPairs> > Parameters;
        
    typedef SireBase::PairMatrix<double> EnergyWorkspace;
    typedef SireBase::PairMatrix<SireMaths::DistVector> ForceWorkspace;
    typedef SireBase::PairMatrix<SireMaths::DistVector> FieldWorkspace;
    typedef SireBase::PairMatrix<double> PotentialWorkspace;

    typedef CoulombProbe Probe;

    typedef SireFF::detail::FFMolecule3D<IntraCoulombPotential> Molecule;
    typedef SireFF::detail::FFMolecules3D<IntraCoulombPotential> Molecules;

    typedef SireFF::detail::ChangedMolecule<Molecule> ChangedMolecule;

    IntraCoulombPotential();
    
    IntraCoulombPotential(const IntraCoulombPotential &other);
    
    ~IntraCoulombPotential();
    
    IntraCoulombPotential& operator=(const IntraCoulombPotential &other);
    
    static const char* typeName()
    {
        return "SireMM::IntraCoulombPotential";
    }
    
    const char* what() const
    {
        return IntraCoulombPotential::typeName();
    }
    
    static ParameterNames parameters()
    {
        return ParameterNames();
    }
    
    IntraCoulombPotential::Parameters 
    getParameters(const PartialMolecule &molecule,
                  const PropertyMap &map = PropertyMap());
    
    IntraCoulombPotential::Parameters
    updateParameters(const IntraCoulombPotential::Parameters &old_params,
                     const PartialMolecule &old_molecule,
                     const PartialMolecule &new_molecule,
                     const PropertyMap &map = PropertyMap());
                     
    IntraCoulombPotential::Parameters
    updateParameters(const IntraCoulombPotential::Parameters &old_params,
                     const PartialMolecule &old_molecule,
                     const PartialMolecule &new_molecule,
                     const PropertyMap &old_map, const PropertyMap &new_map);
    
    IntraCoulombPotential::Molecule
    parameterise(const PartialMolecule &molecule,
                 const PropertyMap &map = PropertyMap());
    
    IntraCoulombPotential::Molecules 
    parameterise(const MoleculeGroup &molecules,
                 const PropertyMap &map = PropertyMap());

    void calculateEnergy(const IntraCoulombPotential::Molecule &mol, 
                         IntraCoulombPotential::Energy &energy,
                         IntraCoulombPotential::EnergyWorkspace &workspace,
                         double scale_energy=1) const;

    void calculateEnergy(const IntraCoulombPotential::Molecule &mol,
                         const IntraCoulombPotential::Molecule &rest_of_mol,
                         IntraCoulombPotential::Energy &energy,
                         IntraCoulombPotential::EnergyWorkspace &workspace,
                         double scale_energy=1) const;

    void calculateForce(const IntraCoulombPotential::Molecule &mol, 
                        MolForceTable &forces,
                        IntraCoulombPotential::ForceWorkspace &workspace,
                        double scale_force=1) const;

    void calculateForce(const IntraCoulombPotential::Molecule &mol,
                        const IntraCoulombPotential::Molecule &rest_of_mol,
                        MolForceTable &forces,
                        IntraCoulombPotential::ForceWorkspace &workspace,
                        double scale_force=1) const;

    void calculateForce(const IntraCoulombPotential::Molecule &mol, 
                        MolForceTable &forces,
                        const Symbol &symbol,
                        const Components &components,
                        IntraCoulombPotential::ForceWorkspace &workspace,
                        double scale_force=1) const;

    void calculateForce(const IntraCoulombPotential::Molecule &mol,
                        const IntraCoulombPotential::Molecule &rest_of_mol,
                        MolForceTable &forces,
                        const Symbol &symbol,
                        const Components &components,
                        IntraCoulombPotential::ForceWorkspace &workspace,
                        double scale_force=1) const;

    void calculateCoulombForce(const IntraCoulombPotential::Molecule &mol,
                               const IntraCoulombPotential::Molecule &rest_of_mol,
                               MolForceTable &forces,
                               IntraCoulombPotential::ForceWorkspace &workspace,
                               double scale_force=1) const;
                               
    void calculateCoulombForce(const IntraCoulombPotential::Molecule &mol,
                               MolForceTable &forces,
                               IntraCoulombPotential::ForceWorkspace &workspace,
                               double scale_force=1) const;

    void calculateField(const IntraCoulombPotential::Molecule &mol, 
                        const CoulombProbe &probe,
                        MolFieldTable &fields,
                        IntraCoulombPotential::FieldWorkspace &workspace,
                        double scale_field=1) const;

    void calculateField(const IntraCoulombPotential::Molecule &mol,
                        const IntraCoulombPotential::Molecule &rest_of_mol,
                        const CoulombProbe &probe,
                        MolFieldTable &fields,
                        IntraCoulombPotential::FieldWorkspace &workspace,
                        double scale_field=1) const;

    void calculateField(const IntraCoulombPotential::Molecule &mol, 
                        const CoulombProbe &probe,
                        MolFieldTable &fields,
                        const Symbol &symbol,
                        const Components &components,
                        IntraCoulombPotential::FieldWorkspace &workspace,
                        double scale_field=1) const;

    void calculateField(const IntraCoulombPotential::Molecule &mol,
                        const IntraCoulombPotential::Molecule &rest_of_mol,
                        const CoulombProbe &probe,
                        MolFieldTable &fields,
                        const Symbol &symbol,
                        const Components &components,
                        IntraCoulombPotential::FieldWorkspace &workspace,
                        double scale_field=1) const;

    void calculateField(const IntraCoulombPotential::Molecule &mol, 
                        const CoulombProbe &probe,
                        GridFieldTable &fields,
                        IntraCoulombPotential::FieldWorkspace &workspace,
                        double scale_field=1) const;

    void calculateField(const IntraCoulombPotential::Molecule &mol, 
                        const CoulombProbe &probe,
                        GridFieldTable &fields,
                        const Symbol &symbol,
                        const Components &components,
                        IntraCoulombPotential::FieldWorkspace &workspace,
                        double scale_field=1) const;

    void calculatePotential(const IntraCoulombPotential::Molecule &mol, 
                            const CoulombProbe &probe,
                            MolPotentialTable &potentials,
                            IntraCoulombPotential::PotentialWorkspace &workspace,
                            double scale_potential=1) const;

    void calculatePotential(const IntraCoulombPotential::Molecule &mol,
                            const IntraCoulombPotential::Molecule &rest_of_mol,
                            const CoulombProbe &probe,
                            MolPotentialTable &potentials,
                            IntraCoulombPotential::PotentialWorkspace &workspace,
                            double scale_potential=1) const;

    void calculatePotential(const IntraCoulombPotential::Molecule &mol, 
                            const CoulombProbe &probe,
                            MolPotentialTable &potentials,
                            const Symbol &symbol,
                            const Components &components,
                            IntraCoulombPotential::PotentialWorkspace &workspace,
                            double scale_potential=1) const;

    void calculatePotential(const IntraCoulombPotential::Molecule &mol,
                            const IntraCoulombPotential::Molecule &rest_of_mol,
                            const CoulombProbe &probe,
                            MolPotentialTable &potentials,
                            const Symbol &symbol,
                            const Components &components,
                            IntraCoulombPotential::PotentialWorkspace &workspace,
                            double scale_potential=1) const;

    void calculatePotential(const IntraCoulombPotential::Molecule &mol, 
                            const CoulombProbe &probe,
                            GridPotentialTable &potentials,
                            IntraCoulombPotential::PotentialWorkspace &workspace,
                            double scale_potential=1) const;

    void calculatePotential(const IntraCoulombPotential::Molecule &mol, 
                            const CoulombProbe &probe,
                            GridPotentialTable &potentials,
                            const Symbol &symbol,
                            const Components &components,
                            IntraCoulombPotential::PotentialWorkspace &workspace,
                            double scale_potential=1) const;

    void calculateCoulombField(const IntraCoulombPotential::Molecule &mol,
                               const CoulombProbe &probe,
                               MolFieldTable &fields,
                               IntraCoulombPotential::FieldWorkspace &workspace,
                               double scale_field=1) const;

    void calculateCoulombField(const IntraCoulombPotential::Molecule &mol,
                               const CoulombProbe &probe,
                               GridFieldTable &fields,
                               IntraCoulombPotential::FieldWorkspace &workspace,
                               double scale_field=1) const;

    void calculateCoulombField(const IntraCoulombPotential::Molecule &mol,
                               const IntraCoulombPotential::Molecule &rest_of_mol,
                               const CoulombProbe &probe,
                               MolFieldTable &fields,
                               IntraCoulombPotential::FieldWorkspace &workspace,
                               double scale_field=1) const;

    void calculateCoulombPotential(const IntraCoulombPotential::Molecule &mol,
                                   const CoulombProbe &probe,
                                   MolPotentialTable &potentials,
                                   IntraCoulombPotential::PotentialWorkspace &workspace,
                                   double scale_potential=1) const;

    void calculateCoulombPotential(const IntraCoulombPotential::Molecule &mol,
                                   const CoulombProbe &probe,
                                   GridPotentialTable &fields,
                                   IntraCoulombPotential::PotentialWorkspace &workspace,
                                   double scale_potential=1) const;

    void calculateCoulombPotential(const IntraCoulombPotential::Molecule &mol,
                                   const IntraCoulombPotential::Molecule &rest_of_mol,
                                   const CoulombProbe &probe,
                                   MolPotentialTable &potentials,
                                   IntraCoulombPotential::PotentialWorkspace &workspace,
                                   double scale_potential=1) const;
                               
private:
    double totalCharge(const IntraCoulombPotential::Parameters::Array &params) const;

    void assertCompatible(const IntraCoulombPotential::Molecule &mol,
                          const IntraCoulombPotential::Molecule &rest_of_mol) const;

    void calculateEnergy(const CoulombNBPairs::CGPairs &group_pairs,
                         IntraCoulombPotential::EnergyWorkspace &workspace,
                         const Parameter *params0_array,
                         const Parameter *params1_array,
                         const quint32 nats0, const quint32 nats1,
                         double &icnrg) const;

    void calculateEnergy(const CoulombNBPairs::CGPairs &group_pairs,
                         const QSet<SireID::Index> &atoms0, 
                         const QSet<SireID::Index> &atoms1,
                         IntraCoulombPotential::EnergyWorkspace &workspace,
                         const Parameter *params0_array,
                         const Parameter *params1_array,
                         const quint32 nats0, const quint32 nats1,
                         double &icnrg) const;

    void calculateForce(const CoulombNBPairs::CGPairs &group_pairs,
                        const CoordGroup &group0, const CoordGroup &group1,
                        const double mindist,
                        IntraCoulombPotential::ForceWorkspace &workspace,
                        const IntraCoulombPotential::Parameter *params0_array,
                        const IntraCoulombPotential::Parameter *params1_array,
                        const quint32 nats0, const quint32 nats1,
                        const double shift_coul,
                        Vector *group_forces0_array,
                        const double scale_force) const;

    void calculateForce(const CoulombNBPairs::CGPairs &group_pairs,
                        const QSet<SireID::Index> &atoms0,
                        const QSet<SireID::Index> &atoms1,
                        const CoordGroup &group0, const CoordGroup &group1,
                        const double mindist,
                        IntraCoulombPotential::ForceWorkspace &workspace,
                        const IntraCoulombPotential::Parameter *params0_array,
                        const IntraCoulombPotential::Parameter *params1_array,
                        const double shift_coul,
                        Vector *group_forces0_array,
                        const double scale_force) const;

    void calculateCoulombForce(const CoulombNBPairs::CGPairs &group_pairs,
                               const CoordGroup &group0, const CoordGroup &group1,
                               const double mindist,
                               IntraCoulombPotential::ForceWorkspace &workspace,
                               const IntraCoulombPotential::Parameter *params0_array,
                               const IntraCoulombPotential::Parameter *params1_array,
                               const quint32 nats0, const quint32 nats1,
                               const double shift_coul,
                               Vector *group_forces0_array,
                               const double scale_force) const;

    void calculateCoulombForce(const CoulombNBPairs::CGPairs &group_pairs,
                               const QSet<SireID::Index> &atoms0,
                               const QSet<SireID::Index> &atoms1,
                               const CoordGroup &group0, const CoordGroup &group1,
                               const double mindist,
                               IntraCoulombPotential::ForceWorkspace &workspace,
                               const IntraCoulombPotential::Parameter *params0_array,
                               const IntraCoulombPotential::Parameter *params1_array,
                               const double shift_coul,
                               Vector *group_forces0_array,
                               const double scale_force) const;

    void throwMissingForceComponent(const Symbol &symbol,
                                    const Components &components) const;
};

/** This small class is used to hide most of the public interfaces of the 
    CoulombPotential derived class, so that only the property-related functions
    are publically available. This provides a better interface when deriving
    a full forcefield class from a CLJ potential.
    
    @author Christopher Woods
*/
template<class CoulPot>
class CoulombPotentialInterface : protected CoulPot
{

friend QDataStream& ::operator<<<>(QDataStream&, 
                                   const CoulombPotentialInterface<CoulPot>&);
friend QDataStream& ::operator>><>(QDataStream&, CoulombPotentialInterface<CoulPot>&);

public:
    CoulombPotentialInterface() : CoulPot()
    {}
    
    CoulombPotentialInterface(const CoulombPotentialInterface &other) : CoulPot(other)
    {}
    
    ~CoulombPotentialInterface()
    {}
    
    static typename CoulPot::ParameterNames parameters()
    {
        return CoulPot::parameters();
    }
    
    const Properties& properties() const
    {
        return CoulPot::properties();
    }
    
    const Property& property(const QString &name) const
    {
        return CoulPot::property(name);
    }
    
    bool containsProperty(const QString &name) const
    {
        return CoulPot::containsProperty(name);
    }
    
    bool setProperty(const QString &name, const Property &value)
    {
        return CoulPot::setProperty(name, value);
    }

    bool setSpace(const Space &new_space)
    {
        return CoulPot::setSpace(new_space);
    }
    
    bool setSwitchingFunction(const SwitchingFunction &new_switchfunc)
    {
        return CoulPot::setSwitchingFunction(new_switchfunc);
    }
    
    bool setShiftElectrostatics(bool switchelectro)
    {
        return CoulPot::setShiftElectrostatics(switchelectro);
    }
    
    const Space& space() const
    {
        return CoulPot::space();
    }
    
    const SwitchingFunction& switchingFunction() const
    {
        return CoulPot::switchingFunction();
    }

    bool shiftElectrostatics() const
    {
        return CoulPot::shiftElectrostatics();
    }
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

//////
////// Inline functions of InterCoulombPotential
//////

/** Calculate the coulomb energy between the passed pair
    of molecules and add this onto 'energy'. This uses
    the passed workspace to perform the calculation */
inline void 
InterCoulombPotential::calculateEnergy(const InterCoulombPotential::Molecule &mol0,
                                       const InterCoulombPotential::Molecule &mol1,
                                       InterCoulombPotential::Energy &energy,
                                       InterCoulombPotential::EnergyWorkspace &workspace,
                                       double scale_energy) const
{
    if (scale_energy != 0 and 
        not (mol0.isEmpty() or mol1.isEmpty()))
    {
        this->_pvt_calculateEnergy(mol0, mol1, energy, workspace, scale_energy);
    }
}

/** Calculate the coulomb and LJ energy between the passed pair
    of molecules and add these energies on mol0 onto energies. This uses
    the passed workspace to perform the calculation */
inline void 
InterCoulombPotential::calculateEnergy(const InterCoulombPotential::Molecule &mol0,
				       const InterCoulombPotential::Molecule &mol1,
				       MolEnergyTable &energies0,
				       InterCoulombPotential::EnergyWorkspace &workspace,
				       double scale_energy) const
{
    throw SireError::incomplete_code( QObject::tr(
            "InterCoulombPotential does not yet support this energy calculations!"), CODELOC );
}

/** Calculate the coulomb forces on the atoms between the passed pair
    of molecules and add the forces on 'mol0' onto 'forces'. This uses
    the passed workspace to perform the calculation. The forces
    are scaled by the optional 'scaled_forces' */
inline void 
InterCoulombPotential::calculateCoulombForce(const InterCoulombPotential::Molecule &mol0, 
                                             const InterCoulombPotential::Molecule &mol1,
                                             MolForceTable &forces0, 
                                             InterCoulombPotential::ForceWorkspace &workspace,
                                             double scale_force) const
{
    if ( scale_force != 0 and 
         not (mol0.isEmpty() or mol1.isEmpty()) and
         not spce->beyond(switchfunc->cutoffDistance(),
                          mol0.aaBox(), mol1.aaBox()) )
    {
        this->_pvt_calculateCoulombForce(mol0, mol1, forces0,
                                         workspace, scale_force);
    }
}

/** Calculate the coulomb forces on the atoms between the passed pair
    of molecules and add the forces on 'mol0' onto 'forces'. This uses
    the passed workspace to perform the calculation. The forces
    are scaled by the optional 'scaled_forces' */
inline void 
InterCoulombPotential::calculateForce(const InterCoulombPotential::Molecule &mol0, 
                                      const InterCoulombPotential::Molecule &mol1,
                                      MolForceTable &forces0, 
                                      InterCoulombPotential::ForceWorkspace &workspace,
                                      double scale_force) const
{
    this->calculateCoulombForce(mol0, mol1, forces0, workspace, scale_force);
}

/** Calculate the component of the force represented by 'symbol' between the 
    passed pair of molecules, and add the forces on 'mol0' onto 'forces0'.
    This uses the passed workspace to perform the calculation. The forces
    are scaled by the optional 'scaled_forces' */
inline void 
InterCoulombPotential::calculateForce(const InterCoulombPotential::Molecule &mol0,
                                      const InterCoulombPotential::Molecule &mol1,
                                      MolForceTable &forces0,
                                      const Symbol &symbol,
                                      const InterCoulombPotential::Components &components,
                                      InterCoulombPotential::ForceWorkspace &workspace,
                                      double scale_force) const
{
    if (symbol == components.total())
        this->calculateForce(mol0, mol1, forces0, workspace, scale_force);
        
    else
        throwMissingForceComponent(symbol, components);
}

inline void 
InterCoulombPotential::calculateCoulombField(const InterCoulombPotential::Molecule &mol0, 
                                      const InterCoulombPotential::Molecule &mol1,
                                      const CoulombProbe &probe,
                                      MolFieldTable &fields0, 
                                      InterCoulombPotential::FieldWorkspace &workspace,
                                      double scale_field) const
{
    if ( scale_field != 0 and 
         not (mol0.isEmpty() or mol1.isEmpty()) and
         not spce->beyond(switchfunc->cutoffDistance(),
                          mol0.aaBox(), mol1.aaBox()) )
    {
        this->_pvt_calculateCoulombField(mol0, mol1, probe, fields0,
                                         workspace, scale_field);
    }
}

inline void 
InterCoulombPotential::calculateField(const InterCoulombPotential::Molecule &mol0, 
                                      const InterCoulombPotential::Molecule &mol1,
                                      const CoulombProbe &probe,
                                      MolFieldTable &fields0, 
                                      InterCoulombPotential::FieldWorkspace &workspace,
                                      double scale_field) const
{
    this->calculateCoulombField(mol0, mol1, probe, fields0, workspace, scale_field);
}

inline void 
InterCoulombPotential::calculateField(const InterCoulombPotential::Molecule &mol0,
                                      const InterCoulombPotential::Molecule &mol1,
                                      const CoulombProbe &probe,
                                      MolFieldTable &fields0,
                                      const Symbol &symbol,
                                      const InterCoulombPotential::Components &components,
                                      InterCoulombPotential::FieldWorkspace &workspace,
                                      double scale_field) const
{
    if (symbol == components.total())
        this->calculateField(mol0, mol1, probe, fields0, workspace, scale_field);
        
    else
        throwMissingFieldComponent(symbol, components);
}

inline void 
InterCoulombPotential::calculateCoulombField(
                                      const InterCoulombPotential::Molecule &mol, 
                                      const CoulombProbe &probe,
                                      GridFieldTable &fields, 
                                      InterCoulombPotential::FieldWorkspace &workspace,
                                      double scale_field) const
{
    if ( scale_field != 0 and 
         not (mol.isEmpty()) and
         not spce->beyond(switchfunc->cutoffDistance(),
                          mol.aaBox(), fields.grid().aaBox()) )
    {
        this->_pvt_calculateCoulombField(mol, probe, fields,
                                         workspace, scale_field);
    }
}

inline void 
InterCoulombPotential::calculateField(const InterCoulombPotential::Molecule &mol, 
                                      const CoulombProbe &probe,
                                      GridFieldTable &fields0, 
                                      InterCoulombPotential::FieldWorkspace &workspace,
                                      double scale_field) const
{
    this->calculateCoulombField(mol, probe, fields0, workspace, scale_field);
}

inline void 
InterCoulombPotential::calculateField(const InterCoulombPotential::Molecule &mol,
                                      const CoulombProbe &probe,
                                      GridFieldTable &fields0,
                                      const Symbol &symbol,
                                      const InterCoulombPotential::Components &components,
                                      InterCoulombPotential::FieldWorkspace &workspace,
                                      double scale_field) const
{
    if (symbol == components.total())
        this->calculateField(mol, probe, fields0, workspace, scale_field);
        
    else
        throwMissingFieldComponent(symbol, components);
}

inline void 
InterCoulombPotential::calculateCoulombPotential(
                                    const InterCoulombPotential::Molecule &mol0, 
                                    const InterCoulombPotential::Molecule &mol1,
                                    const CoulombProbe &probe,
                                    MolPotentialTable &pots0, 
                                    InterCoulombPotential::PotentialWorkspace &workspace,
                                    double scale_potential) const
{
    if ( scale_potential != 0 and 
         not (mol0.isEmpty() or mol1.isEmpty()) and
         not spce->beyond(switchfunc->cutoffDistance(),
                          mol0.aaBox(), mol1.aaBox()) )
    {
        this->_pvt_calculateCoulombPotential(mol0, mol1, probe, pots0,
                                             workspace, scale_potential);
    }
}

inline void 
InterCoulombPotential::calculatePotential(
                                   const InterCoulombPotential::Molecule &mol0, 
                                   const InterCoulombPotential::Molecule &mol1,
                                   const CoulombProbe &probe,
                                   MolPotentialTable &pots0, 
                                   InterCoulombPotential::PotentialWorkspace &workspace,
                                   double scale_potential) const
{
    this->calculateCoulombPotential(mol0, mol1, probe, pots0, workspace, scale_potential);
}

inline void 
InterCoulombPotential::calculatePotential(
                                const InterCoulombPotential::Molecule &mol0,
                                const InterCoulombPotential::Molecule &mol1,
                                const CoulombProbe &probe,
                                MolPotentialTable &pots0,
                                const Symbol &symbol,
                                const InterCoulombPotential::Components &components,
                                InterCoulombPotential::PotentialWorkspace &workspace,
                                double scale_potential) const
{
    if (symbol == components.total())
        this->calculatePotential(mol0, mol1, probe, pots0, workspace, scale_potential);
        
    else
        throwMissingFieldComponent(symbol, components);
}

inline void 
InterCoulombPotential::calculateCoulombPotential(
                                   const InterCoulombPotential::Molecule &mol, 
                                   const CoulombProbe &probe,
                                   GridPotentialTable &pots, 
                                   InterCoulombPotential::PotentialWorkspace &workspace,
                                   double scale_potential) const
{
    if ( scale_potential != 0 and 
         not (mol.isEmpty()) and
         not spce->beyond(switchfunc->cutoffDistance(),
                          mol.aaBox(), pots.grid().aaBox()) )
    {
        this->_pvt_calculateCoulombPotential(mol, probe, pots,
                                             workspace, scale_potential);
    }
}

inline void 
InterCoulombPotential::calculatePotential(
                                const InterCoulombPotential::Molecule &mol, 
                                const CoulombProbe &probe,
                                GridPotentialTable &pots, 
                                InterCoulombPotential::PotentialWorkspace &workspace,
                                double scale_potential) const
{
    this->calculateCoulombPotential(mol, probe, pots, workspace, scale_potential);
}

inline void 
InterCoulombPotential::calculatePotential(
                                const InterCoulombPotential::Molecule &mol,
                                const CoulombProbe &probe,
                                GridPotentialTable &pots,
                                const Symbol &symbol,
                                const InterCoulombPotential::Components &components,
                                InterCoulombPotential::PotentialWorkspace &workspace,
                                double scale_potential) const
{
    if (symbol == components.total())
        this->calculatePotential(mol, probe, pots, workspace, scale_potential);
        
    else
        throwMissingFieldComponent(symbol, components);
}

//////
////// Inline functions of IntraCLJPotential
//////

/** Calculate the forces represented by the symbol 'symbol' between the 
    atoms in the molecule 'mol' and add these forces onto 'forces'. This uses
    the passed workspace to perform the calculation */
inline void 
IntraCoulombPotential::calculateForce(const IntraCoulombPotential::Molecule &mol, 
                                      MolForceTable &forces,
                                      const Symbol &symbol,
                                      const IntraCoulombPotential::Components &components,
                                      IntraCoulombPotential::ForceWorkspace &workspace,
                                      double scale_force) const
{
    if (symbol == components.total())
        this->calculateForce(mol, forces, workspace, scale_force);
        
    else
        throwMissingForceComponent(symbol, components);
}

/** Calculate the forces represented by the symbol 'symbol' acting on the 
    atoms in 'mol1' caused by the atoms in the rest of the same molecule 
    in 'rest_of_mol', and add these forces onto 'forces'. This uses the 
    passed workspace to perform the calculation
    
    \throw SireError::incompatible_error
*/
inline void
IntraCoulombPotential::calculateForce(const IntraCoulombPotential::Molecule &mol,
                                      const IntraCoulombPotential::Molecule &rest_of_mol,
                                      MolForceTable &forces,
                                      const Symbol &symbol,
                                      const IntraCoulombPotential::Components &components,
                                      IntraCoulombPotential::ForceWorkspace &workspace,
                                      double scale_force) const
{
    if (symbol == components.total())
        this->calculateForce(mol, rest_of_mol, forces, workspace, scale_force);
        
    else
        throwMissingForceComponent(symbol, components);
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::detail::ChargeParameter&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::detail::ChargeParameter&);

template<class CoulPot>
QDataStream& operator<<(QDataStream &ds,
                        const SireMM::CoulombPotentialInterface<CoulPot> &coulpot)
{
    ds << static_cast<const CoulPot&>(coulpot);
    return ds;
}

template<class CoulPot>
QDataStream& operator>>(QDataStream &ds,
                        SireMM::CoulombPotentialInterface<CoulPot> &coulpot)
{
    ds >> static_cast<CoulPot&>(coulpot);
    return ds;
}

Q_DECLARE_TYPEINFO( SireMM::detail::ChargeParameter, Q_MOVABLE_TYPE );

SIRE_EXPOSE_CLASS( SireMM::ChargeParameterName )
SIRE_EXPOSE_CLASS( SireMM::ChargeParameterName3D )
SIRE_EXPOSE_CLASS( SireMM::ScaledChargeParameterNames3D )

SIRE_END_HEADER

#endif
