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

#ifndef SIREMM_LJPOTENTIAL_H
#define SIREMM_LJPOTENTIAL_H

#include "SireBase/properties.h"
#include "SireBase/propertymap.h"
#include "SireBase/pairmatrix.hpp"
#include "SireBase/packedarray2d.hpp"

#include "SireVol/space.h"

#include "SireMol/atomproperty.hpp"

#include "SireUnits/dimensions.h"

#include "cljcomponent.h"
#include "cljnbpairs.h"
#include "cljprobe.h"

#include "ljparameter.h"
#include "atomljs.h"

#include "detail/intrascaledatomicparameters.hpp"

#include "ljparameterdb.h"
#include "switchingfunction.h"

#include "SireFF/energytable.h"
#include "SireFF/forcetable.h"
#include "SireFF/fieldtable.h"
#include "SireFF/potentialtable.h"
#include "SireFF/detail/ffmolecules3d.h"

SIRE_BEGIN_HEADER

namespace SireMM
{
class InterLJPotential;
class IntraLJPotential;
class LJPotential;

template<class LJPot>
class LJPotentialInterface;

namespace detail{ class LJParamID; }
}

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::InterLJPotential&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::InterLJPotential&);

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::IntraLJPotential&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::IntraLJPotential&);

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::LJPotential&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::LJPotential&);

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::detail::LJParamID&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::detail::LJParamID&);

template<class LJPot>
QDataStream& operator<<(QDataStream&, const SireMM::LJPotentialInterface<LJPot>&);
template<class LJPot>
QDataStream& operator>>(QDataStream&, SireMM::LJPotentialInterface<LJPot>&);

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

/** This class provides the default name of the 
    property that contains the LJ parameters */
class SIREMM_EXPORT LJParameterName
{
public:
    LJParameterName()
    {}
    
    ~LJParameterName()
    {}
    
    const QString& lj() const
    {
        return lj_param;
    }
    
private:
    static QString lj_param;
};

/** This class provides the default name of the properties
    that contain the LJ and 3D coordinates properties */
class SIREMM_EXPORT LJParameterName3D : public LJParameterName,
                                        public SireFF::detail::Coords3DParameterName
{
public:
    LJParameterName3D() : LJParameterName(),
                          SireFF::detail::Coords3DParameterName()
    {}
    
    ~LJParameterName3D()
    {}
};

/** This class provides the default name of the properties 
    that contain the LJ, intramolecular NB scale parameters and
    3D coordinates properties */
class SIREMM_EXPORT ScaledLJParameterNames3D : public LJParameterName3D,
                                               public detail::IntraScaleParameterName
{
public:
    ScaledLJParameterNames3D() : LJParameterName3D(),
                                 detail::IntraScaleParameterName()
    {}
    
    ~ScaledLJParameterNames3D()
    {}
};

namespace detail
{

/** This class holds the LJ parameter used by both the Inter- and Intra-
    CLJPotentials. It is just the ID of the LJ
    parameter in the singleton LJParameterDB database
    
    @author Christopher Woods
*/
class SIREMM_EXPORT LJParamID
{
public:
    LJParamID(quint32 lj_id=0) : ljid(lj_id)
    {}
    
    LJParamID(const LJParamID &other) : ljid(other.ljid)
    {}
    
    ~LJParamID()
    {}
    
    bool operator==(const LJParamID &other) const
    {
        return ljid == other.ljid;
    }
    
    bool operator!=(const LJParamID &other) const
    {
        return ljid != other.ljid;
    }
    
    quint32 ljid;
};

} // end of namespace detail

/** This is the common base class of InterLJPotential and IntraLJPotential

    @author Christopher Woods
*/
class SIREMM_EXPORT LJPotential
{

friend SIREMM_EXPORT QDataStream& ::operator<<(QDataStream&, const LJPotential&);
friend SIREMM_EXPORT QDataStream& ::operator>>(QDataStream&, LJPotential&);

public:
    virtual ~LJPotential();

    const Properties& properties() const;
    const Property& property(const QString &name) const;
    bool containsProperty(const QString &name) const;
    
    bool setProperty(const QString &name, const Property &value);

    bool setSpace(const Space &new_space);
    bool setSwitchingFunction(const SwitchingFunction &new_switchfunc);
    bool setCombiningRules(const QString &combiningrules);
    
    const Space& space() const;
    const SwitchingFunction& switchingFunction() const;
    const QString& combiningRules() const;

protected:
    LJPotential();
    LJPotential(const LJPotential &other);
    
    LJPotential& operator=(const LJPotential &other);

    void startEvaluation();
    void finishedEvaluation();
    
    virtual void changedPotential()=0;
    
    /** All possible LJ parameter pair combinations, arranged
        in a symmetric 2D array */
    LJPairMatrix ljpairs;

    /** The current values of the properties of this functional */
    Properties props;
    
    /** The space in which this functional operates */
    SpacePtr spce;
    
    /** The nonbonded switching function */
    SwitchFuncPtr switchfunc;

    /** The combining rules to use to get mixed LJ parameters */
    LJParameterDB::CombiningRules combining_rules;
    
    /** Whether or not the LJ pair matrix needs to be rebuilt */
    bool need_update_ljpairs;
};

/** This class provides all of the functions and containers  
    necessary to provide an interface to calculate the
    intermolecular interatomic potentials using
    Lennard Jones functions.
    
    This is a 3D potential class, namely it requires that
    the atoms possess 3D coordinates, thereby allowing this
    potential to also be used to calculate 3D forces on the atoms.

    This potential has the following properties (parameters);
    
    (1) space               : This is the 3D space in which the molecules exist
    (2) switchingFunction   : This is the switching function used to scale the 
                              energies / forces to zero at the cutoff
    (3) combiningRules      : This is a string specifying the LJ combining rules
                              (currently "arithmetic" or "geometric")
                              
    The molecules used with this potential must contain the following 
    properties (defined in parameters(), e.g. parameters().coordinates())
    
    (1) .coordinates()      : These are the 3D coordinates of each atom, must be
                              a property of type AtomCoords
    (2) .lj()               : These are the LJ parameters for each atom, must
                              be a property of type AtomLJs
    
    @author Christopher Woods
*/
class SIREMM_EXPORT InterLJPotential : public LJPotential
{

friend SIREMM_EXPORT QDataStream& ::operator<<(QDataStream&, const InterLJPotential&);
friend SIREMM_EXPORT QDataStream& ::operator>>(QDataStream&, InterLJPotential&);

public:

    typedef LJEnergy Energy;
    typedef Energy::Components Components;

    typedef LJParameterName3D ParameterNames;

    typedef detail::LJParamID Parameter;
    typedef SireFF::detail::AtomicParameters3D<Parameter> Parameters;
    
    typedef SireBase::PairMatrix<double> EnergyWorkspace;
    typedef SireBase::PairMatrix<SireMaths::DistVector> ForceWorkspace;
    typedef SireBase::PairMatrix<SireMaths::DistVector> FieldWorkspace;
    typedef SireBase::PairMatrix<double> PotentialWorkspace;

    typedef LJProbe Probe;

    typedef SireFF::detail::FFMolecule3D<InterLJPotential> Molecule;
    typedef SireFF::detail::FFMolecules3D<InterLJPotential> Molecules;

    typedef SireFF::detail::ChangedMolecule<Molecule> ChangedMolecule;

    InterLJPotential();
    
    InterLJPotential(const InterLJPotential &other);
    
    ~InterLJPotential();
    
    InterLJPotential& operator=(const InterLJPotential &other);

    static const char* typeName()
    {
        return "SireMM::InterLJPotential";
    }

    const char* what() const
    {
        return InterLJPotential::typeName();
    }
    
    static ParameterNames parameters()
    {
        return ParameterNames();
    }
    
    InterLJPotential::Parameters 
    getParameters(const PartialMolecule &molecule,
                  const PropertyMap &map = PropertyMap());
    
    InterLJPotential::Parameters
    updateParameters(const InterLJPotential::Parameters &old_params,
                     const PartialMolecule &old_molecule,
                     const PartialMolecule &new_molecule,
                     const PropertyMap &map = PropertyMap());
                     
    InterLJPotential::Parameters
    updateParameters(const InterLJPotential::Parameters &old_params,
                     const PartialMolecule &old_molecule,
                     const PartialMolecule &new_molecule,
                     const PropertyMap &old_map, const PropertyMap &new_map);
    
    InterLJPotential::Molecule
    parameterise(const PartialMolecule &molecule,
                 const PropertyMap &map = PropertyMap());
    
    InterLJPotential::Molecules 
    parameterise(const MoleculeGroup &molecules,
                 const PropertyMap &map = PropertyMap());

    void calculateEnergy(const InterLJPotential::Molecule &mol0, 
                         const InterLJPotential::Molecule &mol1,
                         InterLJPotential::Energy &energy, 
                         InterLJPotential::EnergyWorkspace &workspace,
                         double scale_energy=1) const;

    void calculateEnergy(const InterLJPotential::Molecule &mol0, 
                         const InterLJPotential::Molecule &mol1,
                         MolEnergyTable &energies0, 
                         InterLJPotential::EnergyWorkspace &workspace,
                         double scale_energy=1) const;

    void calculateForce(const InterLJPotential::Molecule &mol0, 
                        const InterLJPotential::Molecule &mol1,
                        MolForceTable &forces0,
                        InterLJPotential::ForceWorkspace &workspace,
                        double scale_force=1) const;

    void calculateForce(const InterLJPotential::Molecule &mol0,
                        const InterLJPotential::Molecule &mol1,
                        MolForceTable &forces0,
                        const Symbol &symbol,
                        const Components &components,
                        InterLJPotential::ForceWorkspace &workspace,
                        double scale_force=1) const;

    void calculateLJForce(const InterLJPotential::Molecule &mol0, 
                          const InterLJPotential::Molecule &mol1,
                          MolForceTable &forces0,
                          InterLJPotential::ForceWorkspace &workspace,
                          double scale_force=1) const;

    void calculateField(const InterLJPotential::Molecule &mol0, 
                        const InterLJPotential::Molecule &mol1,
                        const LJProbe &probe,
                        MolFieldTable &forces0,
                        InterLJPotential::FieldWorkspace &workspace,
                        double scale_field=1) const;

    void calculateField(const InterLJPotential::Molecule &mol0,
                        const InterLJPotential::Molecule &mol1,
                        const LJProbe &probe,
                        MolFieldTable &forces0,
                        const Symbol &symbol,
                        const Components &components,
                        InterLJPotential::FieldWorkspace &workspace,
                        double scale_field=1) const;

    void calculateField(const InterLJPotential::Molecule &mol0,
                        const LJProbe &probe,
                        GridFieldTable &fields,
                        InterLJPotential::FieldWorkspace &workspace,
                        double scale_field=1) const;

    void calculateField(const InterLJPotential::Molecule &mol0,
                        const LJProbe &probe,
                        GridFieldTable &fields,
                        const Symbol &symbol,
                        const Components &components,
                        InterLJPotential::FieldWorkspace &workspace,
                        double scale_field=1) const;

    void calculateLJField(const InterLJPotential::Molecule &mol0, 
                          const InterLJPotential::Molecule &mol1,
                          const LJProbe &probe,
                          MolFieldTable &fields0,
                          InterLJPotential::ForceWorkspace &workspace,
                          double scale_field=1) const;

    void calculateLJField(const InterLJPotential::Molecule &mol0, 
                          const LJProbe &probe,
                          GridFieldTable &fields,
                          InterLJPotential::ForceWorkspace &workspace,
                          double scale_field=1) const;

    void calculatePotential(const InterLJPotential::Molecule &mol0, 
                            const InterLJPotential::Molecule &mol1,
                            const LJProbe &probe,
                            MolPotentialTable &pots0,
                            InterLJPotential::PotentialWorkspace &workspace,
                            double scale_potential=1) const;

    void calculatePotential(const InterLJPotential::Molecule &mol0,
                            const InterLJPotential::Molecule &mol1,
                            const LJProbe &probe,
                            MolPotentialTable &pots0,
                            const Symbol &symbol,
                            const Components &components,
                            InterLJPotential::PotentialWorkspace &workspace,
                            double scale_potential=1) const;

    void calculatePotential(const InterLJPotential::Molecule &mol0,
                            const LJProbe &probe,
                            GridPotentialTable &pots,
                            InterLJPotential::PotentialWorkspace &workspace,
                            double scale_potential=1) const;

    void calculatePotential(const InterLJPotential::Molecule &mol0,
                            const LJProbe &probe,
                            GridPotentialTable &pots,
                            const Symbol &symbol,
                            const Components &components,
                            InterLJPotential::PotentialWorkspace &workspace,
                            double scale_potential=1) const;

    void calculateLJPotential(const InterLJPotential::Molecule &mol0, 
                              const InterLJPotential::Molecule &mol1,
                              const LJProbe &probe,
                              MolPotentialTable &pots0,
                              InterLJPotential::PotentialWorkspace &workspace,
                              double scale_potential=1) const;

    void calculateLJPotential(const InterLJPotential::Molecule &mol0, 
                              const LJProbe &probe,
                              GridPotentialTable &pots,
                              InterLJPotential::PotentialWorkspace &workspace,
                              double scale_potential=1) const;

private:
    void throwMissingForceComponent(const Symbol &symbol,
                                    const Components &components) const;
    void throwMissingFieldComponent(const Symbol &symbol,
                                    const Components &components) const;
    void throwMissingPotentialComponent(const Symbol &symbol,
                                        const Components &components) const;

    void _pvt_calculateEnergy(const InterLJPotential::Molecule &mol0, 
                              const InterLJPotential::Molecule &mol1,
                              InterLJPotential::Energy &energy, 
                              InterLJPotential::EnergyWorkspace &workspace,
                              double scale_energy) const;

    void _pvt_calculateLJForce(const InterLJPotential::Molecule &mol0, 
                               const InterLJPotential::Molecule &mol1,
                               MolForceTable &forces0, 
                               InterLJPotential::ForceWorkspace &workspace,
                               double scale_force) const;

    void _pvt_calculateLJField(const InterLJPotential::Molecule &mol0, 
                               const InterLJPotential::Molecule &mol1,
                               const LJProbe &probe,
                               MolFieldTable &forces0, 
                               InterLJPotential::ForceWorkspace &workspace,
                               double scale_force) const;

    void _pvt_calculateLJField(const InterLJPotential::Molecule &mol, 
                               const LJProbe &probe,
                               GridFieldTable &forces0, 
                               InterLJPotential::ForceWorkspace &workspace,
                               double scale_force) const;

    void _pvt_calculateLJPotential(const InterLJPotential::Molecule &mol0, 
                                   const InterLJPotential::Molecule &mol1,
                                   const LJProbe &probe,
                                   MolPotentialTable &pots0, 
                                   InterLJPotential::PotentialWorkspace &workspace,
                                   double scale_potential) const;

    void _pvt_calculateLJPotential(const InterLJPotential::Molecule &mol, 
                                   const LJProbe &probe,
                                   GridPotentialTable &pots, 
                                   InterLJPotential::PotentialWorkspace &workspace,
                                   double scale_potential) const;
};

/** This class provides all of the functions and containers  
    necessary to provide an interface to calculate the
    intramolecular interatomic potentials using
    Lennard Jones functions.
    
    This is a 3D potential class, namely it requires that
    the atoms possess 3D coordinates, thereby allowing this
    potential to also be used to calculate 3D forces on the atoms.

    This potential has the following properties (parameters);
    
    (1) space               : This is the 3D space in which the molecules exist
    (2) switchingFunction   : This is the switching function used to scale the 
                              energies / forces to zero at the cutoff
    (3) combiningRules      : This is a string specifying the LJ combining rules
                              (currently "arithmetic" or "geometric")
                              
    The molecules used with this potential must contain the following 
    properties (defined in parameters(), e.g. parameters().coordinates())
    
    (1) .coordinates()       : These are the 3D coordinates of each atom, must be
                               a property of type AtomCoords
    (2) .lj()                : These are the LJ parameters for each atom, must
                               be a property of type AtomLJs
    (3) .intraScaleFactors() : These are the intramolecular atomic scaling factors,
                               used to scale the intramolecular LJ
                               energies, must be a property of type LJNBPairs or 
                               CLJNBPairs
    
    @author Christopher Woods
*/
class SIREMM_EXPORT IntraLJPotential : public LJPotential
{

friend SIREMM_EXPORT QDataStream& ::operator<<(QDataStream&, const IntraLJPotential&);
friend SIREMM_EXPORT QDataStream& ::operator>>(QDataStream&, IntraLJPotential&);

public:
    typedef LJEnergy Energy;
    typedef Energy::Components Components;
    
    typedef ScaledLJParameterNames3D ParameterNames;

    typedef detail::LJParamID Parameter;
    
    typedef detail::IntraScaledAtomicParameters<
                  SireFF::detail::AtomicParameters3D<Parameter>,
                  detail::IntraScaledParameters<LJNBPairs> > Parameters;
        
    typedef SireBase::PairMatrix<double> EnergyWorkspace;
    typedef SireBase::PairMatrix<SireMaths::DistVector> ForceWorkspace;
    typedef SireBase::PairMatrix<SireMaths::DistVector> FieldWorkspace;
    typedef SireBase::PairMatrix<double> PotentialWorkspace;

    typedef LJProbe Probe;

    typedef SireFF::detail::FFMolecule3D<IntraLJPotential> Molecule;
    typedef SireFF::detail::FFMolecules3D<IntraLJPotential> Molecules;

    typedef SireFF::detail::ChangedMolecule<Molecule> ChangedMolecule;

    IntraLJPotential();
    
    IntraLJPotential(const IntraLJPotential &other);
    
    ~IntraLJPotential();
    
    IntraLJPotential& operator=(const IntraLJPotential &other);
    
    static const char* typeName()
    {
        return "SireMM::IntraLJPotential";
    }
    
    const char* what() const
    {
        return IntraLJPotential::typeName();
    }
    
    static ParameterNames parameters()
    {
        return ParameterNames();
    }
    
    IntraLJPotential::Parameters 
    getParameters(const PartialMolecule &molecule,
                  const PropertyMap &map = PropertyMap());
    
    IntraLJPotential::Parameters
    updateParameters(const IntraLJPotential::Parameters &old_params,
                     const PartialMolecule &old_molecule,
                     const PartialMolecule &new_molecule,
                     const PropertyMap &map = PropertyMap());
                     
    IntraLJPotential::Parameters
    updateParameters(const IntraLJPotential::Parameters &old_params,
                     const PartialMolecule &old_molecule,
                     const PartialMolecule &new_molecule,
                     const PropertyMap &old_map, const PropertyMap &new_map);
    
    IntraLJPotential::Molecule
    parameterise(const PartialMolecule &molecule,
                 const PropertyMap &map = PropertyMap());
    
    IntraLJPotential::Molecules 
    parameterise(const MoleculeGroup &molecules,
                 const PropertyMap &map = PropertyMap());

    void calculateEnergy(const IntraLJPotential::Molecule &mol, 
                         IntraLJPotential::Energy &energy,
                         IntraLJPotential::EnergyWorkspace &workspace,
                         double scale_energy=1) const;

    void calculateEnergy(const IntraLJPotential::Molecule &mol,
                         const IntraLJPotential::Molecule &rest_of_mol,
                         IntraLJPotential::Energy &energy,
                         IntraLJPotential::EnergyWorkspace &workspace,
                         double scale_energy=1) const;

    void calculateForce(const IntraLJPotential::Molecule &mol, 
                        MolForceTable &forces,
                        IntraLJPotential::ForceWorkspace &workspace,
                        double scale_force=1) const;

    void calculateForce(const IntraLJPotential::Molecule &mol,
                        const IntraLJPotential::Molecule &rest_of_mol,
                        MolForceTable &forces,
                        IntraLJPotential::ForceWorkspace &workspace,
                        double scale_force=1) const;

    void calculateForce(const IntraLJPotential::Molecule &mol, 
                        MolForceTable &forces,
                        const Symbol &symbol,
                        const Components &components,
                        IntraLJPotential::ForceWorkspace &workspace,
                        double scale_force=1) const;

    void calculateForce(const IntraLJPotential::Molecule &mol,
                        const IntraLJPotential::Molecule &rest_of_mol,
                        MolForceTable &forces,
                        const Symbol &symbol,
                        const Components &components,
                        IntraLJPotential::ForceWorkspace &workspace,
                        double scale_force=1) const;
                               
    void calculateLJForce(const IntraLJPotential::Molecule &mol,
                          const IntraLJPotential::Molecule &rest_of_mol,
                          MolForceTable &forces,
                          IntraLJPotential::ForceWorkspace &workspace,
                          double scale_force=1) const;
                               
    void calculateLJForce(const IntraLJPotential::Molecule &mol,
                          MolForceTable &forces,
                          IntraLJPotential::ForceWorkspace &workspace,
                          double scale_force=1) const;

    void calculateField(const IntraLJPotential::Molecule &mol, 
                        const LJProbe &probe,
                        MolFieldTable &fields,
                        IntraLJPotential::FieldWorkspace &workspace,
                        double scale_field=1) const;

    void calculateField(const IntraLJPotential::Molecule &mol, 
                        const LJProbe &probe,
                        MolFieldTable &fields,
                        const Symbol &symbol,
                        const Components &components,
                        IntraLJPotential::FieldWorkspace &workspace,
                        double scale_field=1) const;

    void calculateField(const IntraLJPotential::Molecule &mol0, 
                        const IntraLJPotential::Molecule &mol1,
                        const LJProbe &probe,
                        MolFieldTable &forces0,
                        IntraLJPotential::FieldWorkspace &workspace,
                        double scale_field=1) const;

    void calculateField(const IntraLJPotential::Molecule &mol0,
                        const IntraLJPotential::Molecule &mol1,
                        const LJProbe &probe,
                        MolFieldTable &forces0,
                        const Symbol &symbol,
                        const Components &components,
                        IntraLJPotential::FieldWorkspace &workspace,
                        double scale_field=1) const;

    void calculateField(const IntraLJPotential::Molecule &mol0,
                        const LJProbe &probe,
                        GridFieldTable &fields,
                        IntraLJPotential::FieldWorkspace &workspace,
                        double scale_field=1) const;

    void calculateField(const IntraLJPotential::Molecule &mol0,
                        const LJProbe &probe,
                        GridFieldTable &fields,
                        const Symbol &symbol,
                        const Components &components,
                        IntraLJPotential::FieldWorkspace &workspace,
                        double scale_field=1) const;

    void calculateLJField(const IntraLJPotential::Molecule &mol0, 
                          const IntraLJPotential::Molecule &mol1,
                          const LJProbe &probe,
                          MolFieldTable &fields0,
                          IntraLJPotential::ForceWorkspace &workspace,
                          double scale_field=1) const;

    void calculateLJField(const IntraLJPotential::Molecule &mol0, 
                          const LJProbe &probe,
                          GridFieldTable &fields,
                          IntraLJPotential::ForceWorkspace &workspace,
                          double scale_field=1) const;

    void calculatePotential(const IntraLJPotential::Molecule &mol, 
                            const LJProbe &probe,
                            MolPotentialTable &potentials,
                            IntraLJPotential::PotentialWorkspace &workspace,
                            double scale_potential=1) const;

    void calculatePotential(const IntraLJPotential::Molecule &mol, 
                            const LJProbe &probe,
                            MolPotentialTable &potentials,
                            const Symbol &symbol,
                            const Components &components,
                            IntraLJPotential::PotentialWorkspace &workspace,
                            double scale_potential=1) const;

    void calculatePotential(const IntraLJPotential::Molecule &mol0, 
                            const IntraLJPotential::Molecule &mol1,
                            const LJProbe &probe,
                            MolPotentialTable &pots0,
                            IntraLJPotential::PotentialWorkspace &workspace,
                            double scale_potential=1) const;

    void calculatePotential(const IntraLJPotential::Molecule &mol0,
                            const IntraLJPotential::Molecule &mol1,
                            const LJProbe &probe,
                            MolPotentialTable &pots0,
                            const Symbol &symbol,
                            const Components &components,
                            IntraLJPotential::PotentialWorkspace &workspace,
                            double scale_potential=1) const;

    void calculatePotential(const IntraLJPotential::Molecule &mol0,
                            const LJProbe &probe,
                            GridPotentialTable &pots,
                            IntraLJPotential::PotentialWorkspace &workspace,
                            double scale_potential=1) const;

    void calculatePotential(const IntraLJPotential::Molecule &mol0,
                            const LJProbe &probe,
                            GridPotentialTable &pots,
                            const Symbol &symbol,
                            const Components &components,
                            IntraLJPotential::PotentialWorkspace &workspace,
                            double scale_potential=1) const;

    void calculateLJPotential(const IntraLJPotential::Molecule &mol0, 
                              const IntraLJPotential::Molecule &mol1,
                              const LJProbe &probe,
                              MolPotentialTable &pots0,
                              IntraLJPotential::PotentialWorkspace &workspace,
                              double scale_potential=1) const;

    void calculateLJPotential(const IntraLJPotential::Molecule &mol0, 
                              const LJProbe &probe,
                              GridPotentialTable &pots,
                              IntraLJPotential::PotentialWorkspace &workspace,
                              double scale_potential=1) const;
                               
private:
    void assertCompatible(const IntraLJPotential::Molecule &mol,
                          const IntraLJPotential::Molecule &rest_of_mol) const;

    void calculateEnergy(const LJNBPairs::CGPairs &group_pairs,
                         IntraLJPotential::EnergyWorkspace &workspace,
                         const Parameter *params0_array,
                         const Parameter *params1_array,
                         const quint32 nats0, const quint32 nats1,
                         double &iljnrg) const;

    void calculateEnergy(const LJNBPairs::CGPairs &group_pairs,
                         const QSet<SireID::Index> &atoms0, 
                         const QSet<SireID::Index> &atoms1,
                         IntraLJPotential::EnergyWorkspace &workspace,
                         const Parameter *params0_array,
                         const Parameter *params1_array,
                         const quint32 nats0, const quint32 nats1,
                         double &iljnrg) const;

    void calculateForce(const LJNBPairs::CGPairs &group_pairs,
                        const CoordGroup &group0, const CoordGroup &group1,
                        const double mindist,
                        IntraLJPotential::ForceWorkspace &workspace,
                        const IntraLJPotential::Parameter *params0_array,
                        const IntraLJPotential::Parameter *params1_array,
                        const quint32 nats0, const quint32 nats1,
                        Vector *group_forces0_array,
                        const double scale_force) const;

    void calculateForce(const LJNBPairs::CGPairs &group_pairs,
                        const QSet<SireID::Index> &atoms0,
                        const QSet<SireID::Index> &atoms1,
                        const CoordGroup &group0, const CoordGroup &group1,
                        const double mindist,
                        IntraLJPotential::ForceWorkspace &workspace,
                        const IntraLJPotential::Parameter *params0_array,
                        const IntraLJPotential::Parameter *params1_array,
                        Vector *group_forces0_array,
                        const double scale_force) const;

    void calculateLJForce(const LJNBPairs::CGPairs &group_pairs,
                          const CoordGroup &group0, const CoordGroup &group1,
                          const double mindist,
                          IntraLJPotential::ForceWorkspace &workspace,
                          const IntraLJPotential::Parameter *params0_array,
                          const IntraLJPotential::Parameter *params1_array,
                          const quint32 nats0, const quint32 nats1,
                          Vector *group_forces0_array,
                          const double scale_force) const;

    void calculateLJForce(const LJNBPairs::CGPairs &group_pairs,
                          const QSet<SireID::Index> &atoms0,
                          const QSet<SireID::Index> &atoms1,
                          const CoordGroup &group0, const CoordGroup &group1,
                          const double mindist,
                          IntraLJPotential::ForceWorkspace &workspace,
                          const IntraLJPotential::Parameter *params0_array,
                          const IntraLJPotential::Parameter *params1_array,
                          Vector *group_forces0_array,
                          const double scale_force) const;

    void throwMissingForceComponent(const Symbol &symbol,
                                    const Components &components) const;
};

/** This small class is used to hide most of the public interfaces of the 
    CLJPotential derived class, so that only the property-related functions
    are publically available. This provides a better interface when deriving
    a full forcefield class from a CLJ potential.
    
    @author Christopher Woods
*/
template<class LJPot>
class LJPotentialInterface : protected LJPot
{

friend SIREMM_EXPORT QDataStream& ::operator<<<>(QDataStream&, const LJPotentialInterface<LJPot>&);
friend SIREMM_EXPORT QDataStream& ::operator>><>(QDataStream&, LJPotentialInterface<LJPot>&);

public:
    LJPotentialInterface() : LJPot()
    {}
    
    LJPotentialInterface(const LJPotentialInterface &other) : LJPot(other)
    {}
    
    ~LJPotentialInterface()
    {}
    
    static typename LJPot::ParameterNames parameters()
    {
        return LJPot::parameters();
    }
    
    const Properties& properties() const
    {
        return LJPot::properties();
    }
    
    const Property& property(const QString &name) const
    {
        return LJPot::property(name);
    }
    
    bool containsProperty(const QString &name) const
    {
        return LJPot::containsProperty(name);
    }
    
    bool setProperty(const QString &name, const Property &value)
    {
        return LJPot::setProperty(name, value);
    }

    bool setSpace(const Space &new_space)
    {
        return LJPot::setSpace(new_space);
    }
    
    bool setSwitchingFunction(const SwitchingFunction &new_switchfunc)
    {
        return LJPot::setSwitchingFunction(new_switchfunc);
    }
    
    bool setCombiningRules(const QString &combiningrules)
    {
        return LJPot::setCombiningRules(combiningrules);
    }
    
    const Space& space() const
    {
        return LJPot::space();
    }
    
    const SwitchingFunction& switchingFunction() const
    {
        return LJPot::switchingFunction();
    }
    
    const QString& combiningRules() const
    {
        return LJPot::combiningRules();
    }
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

//////
////// Inline functions of InterLJPotential
//////

/** Calculate the LJ energy between the passed pair
    of molecules and add these energies onto 'energy'. This uses
    the passed workspace to perform the calculation */
SIRE_ALWAYS_INLINE void 
InterLJPotential::calculateEnergy(const InterLJPotential::Molecule &mol0,
                                  const InterLJPotential::Molecule &mol1,
                                  InterLJPotential::Energy &energy,
                                  InterLJPotential::EnergyWorkspace &workspace,
                                  double scale_energy) const
{
    if (scale_energy != 0 and 
        not (mol0.isEmpty() or mol1.isEmpty()))
    {
        this->_pvt_calculateEnergy(mol0, mol1, energy, workspace, scale_energy);
    }
}

SIRE_ALWAYS_INLINE void 
InterLJPotential::calculateEnergy(const InterLJPotential::Molecule &mol0,
                                  const InterLJPotential::Molecule &mol1,
                                  MolEnergyTable &energies0,
                                  InterLJPotential::EnergyWorkspace &workspace,
                                  double scale_energy) const
{
    throw SireError::incomplete_code( QObject::tr(
            "InterLJPotential does not yet support this energy calculations!"), CODELOC );
}

/** Calculate the LJ forces on the atoms between the passed pair
    of molecules and add the forces on 'mol0' onto 'forces'. This uses
    the passed workspace to perform the calculation. The forces
    are scaled by the optional 'scaled_forces' */
SIRE_ALWAYS_INLINE void 
InterLJPotential::calculateLJForce(const InterLJPotential::Molecule &mol0, 
                                   const InterLJPotential::Molecule &mol1,
                                   MolForceTable &forces0, 
                                   InterLJPotential::ForceWorkspace &workspace,
                                   double scale_force) const
{
    if ( scale_force != 0 and 
         not (mol0.isEmpty() or mol1.isEmpty()) and
         not spce->beyond(switchfunc->cutoffDistance(),
                          mol0.aaBox(), mol1.aaBox()) )
    {
        this->_pvt_calculateLJForce(mol0, mol1, forces0,
                                    workspace, scale_force);
    }
}

/** Calculate the LJ forces on the atoms between the passed pair
    of molecules and add the forces on 'mol0' onto 'forces'. This uses
    the passed workspace to perform the calculation. The forces
    are scaled by the optional 'scaled_forces' */
SIRE_ALWAYS_INLINE void 
InterLJPotential::calculateForce(const InterLJPotential::Molecule &mol0, 
                                 const InterLJPotential::Molecule &mol1,
                                 MolForceTable &forces0, 
                                 InterLJPotential::ForceWorkspace &workspace,
                                 double scale_force) const
{
    this->calculateLJForce(mol0, mol1, forces0, workspace, scale_force);
}

/** Calculate the component of the force represented by 'symbol' between the 
    passed pair of molecules, and add the forces on 'mol0' onto 'forces0'.
    This uses the passed workspace to perform the calculation. The forces
    are scaled by the optional 'scaled_forces' */
SIRE_ALWAYS_INLINE void 
InterLJPotential::calculateForce(const InterLJPotential::Molecule &mol0,
                                 const InterLJPotential::Molecule &mol1,
                                 MolForceTable &forces0,
                                 const Symbol &symbol,
                                 const InterLJPotential::Components &components,
                                 InterLJPotential::ForceWorkspace &workspace,
                                 double scale_force) const
{
    if (symbol == components.total())
        this->calculateForce(mol0, mol1, forces0, workspace, scale_force);
        
    else
        throwMissingForceComponent(symbol, components);
}

SIRE_ALWAYS_INLINE void 
InterLJPotential::calculateLJField(const InterLJPotential::Molecule &mol0, 
                                   const InterLJPotential::Molecule &mol1,
                                   const LJProbe &probe,
                                   MolFieldTable &fields0, 
                                   InterLJPotential::FieldWorkspace &workspace,
                                   double scale_field) const
{
    if ( scale_field != 0 and 
         not (mol0.isEmpty() or mol1.isEmpty()) and
         not spce->beyond(switchfunc->cutoffDistance(),
                          mol0.aaBox(), mol1.aaBox()) )
    {
        this->_pvt_calculateLJField(mol0, mol1, probe, fields0,
                                         workspace, scale_field);
    }
}

SIRE_ALWAYS_INLINE void 
InterLJPotential::calculateField(const InterLJPotential::Molecule &mol0, 
                                 const InterLJPotential::Molecule &mol1,
                                 const LJProbe &probe,
                                 MolFieldTable &fields0, 
                                 InterLJPotential::FieldWorkspace &workspace,
                                 double scale_field) const
{
    this->calculateLJField(mol0, mol1, probe, fields0, workspace, scale_field);
}

SIRE_ALWAYS_INLINE void 
InterLJPotential::calculateField(const InterLJPotential::Molecule &mol0,
                                 const InterLJPotential::Molecule &mol1,
                                 const LJProbe &probe,
                                 MolFieldTable &fields0,
                                 const Symbol &symbol,
                                 const InterLJPotential::Components &components,
                                 InterLJPotential::FieldWorkspace &workspace,
                                 double scale_field) const
{
    if (symbol == components.total())
        this->calculateField(mol0, mol1, probe, fields0, workspace, scale_field);
        
    else
        throwMissingFieldComponent(symbol, components);
}

SIRE_ALWAYS_INLINE void 
InterLJPotential::calculateLJField(const InterLJPotential::Molecule &mol, 
                                   const LJProbe &probe,
                                   GridFieldTable &fields, 
                                   InterLJPotential::FieldWorkspace &workspace,
                                   double scale_field) const
{
    if ( scale_field != 0 and 
         not (mol.isEmpty()) and
         not spce->beyond(switchfunc->cutoffDistance(),
                          mol.aaBox(), fields.grid().aaBox()) )
    {
        this->_pvt_calculateLJField(mol, probe, fields,
                                         workspace, scale_field);
    }
}

SIRE_ALWAYS_INLINE void 
InterLJPotential::calculateField(const InterLJPotential::Molecule &mol, 
                                 const LJProbe &probe,
                                 GridFieldTable &fields0, 
                                 InterLJPotential::FieldWorkspace &workspace,
                                 double scale_field) const
{
    this->calculateLJField(mol, probe, fields0, workspace, scale_field);
}

SIRE_ALWAYS_INLINE void 
InterLJPotential::calculateField(const InterLJPotential::Molecule &mol,
                                 const LJProbe &probe,
                                 GridFieldTable &fields0,
                                 const Symbol &symbol,
                                 const InterLJPotential::Components &components,
                                 InterLJPotential::FieldWorkspace &workspace,
                                 double scale_field) const
{
    if (symbol == components.total())
        this->calculateField(mol, probe, fields0, workspace, scale_field);
        
    else
        throwMissingFieldComponent(symbol, components);
}

SIRE_ALWAYS_INLINE void 
InterLJPotential::calculateLJPotential(
                                    const InterLJPotential::Molecule &mol0, 
                                    const InterLJPotential::Molecule &mol1,
                                    const LJProbe &probe,
                                    MolPotentialTable &pots0, 
                                    InterLJPotential::PotentialWorkspace &workspace,
                                    double scale_potential) const
{
    if ( scale_potential != 0 and 
         not (mol0.isEmpty() or mol1.isEmpty()) and
         not spce->beyond(switchfunc->cutoffDistance(),
                          mol0.aaBox(), mol1.aaBox()) )
    {
        this->_pvt_calculateLJPotential(mol0, mol1, probe, pots0,
                                             workspace, scale_potential);
    }
}

SIRE_ALWAYS_INLINE void 
InterLJPotential::calculatePotential(
                                   const InterLJPotential::Molecule &mol0, 
                                   const InterLJPotential::Molecule &mol1,
                                   const LJProbe &probe,
                                   MolPotentialTable &pots0, 
                                   InterLJPotential::PotentialWorkspace &workspace,
                                   double scale_potential) const
{
    this->calculateLJPotential(mol0, mol1, probe, pots0, workspace, scale_potential);
}

SIRE_ALWAYS_INLINE void 
InterLJPotential::calculatePotential(
                                const InterLJPotential::Molecule &mol0,
                                const InterLJPotential::Molecule &mol1,
                                const LJProbe &probe,
                                MolPotentialTable &pots0,
                                const Symbol &symbol,
                                const InterLJPotential::Components &components,
                                InterLJPotential::PotentialWorkspace &workspace,
                                double scale_potential) const
{
    if (symbol == components.total())
        this->calculatePotential(mol0, mol1, probe, pots0, workspace, scale_potential);
        
    else
        throwMissingFieldComponent(symbol, components);
}

SIRE_ALWAYS_INLINE void 
InterLJPotential::calculateLJPotential(
                                   const InterLJPotential::Molecule &mol, 
                                   const LJProbe &probe,
                                   GridPotentialTable &pots, 
                                   InterLJPotential::PotentialWorkspace &workspace,
                                   double scale_potential) const
{
    if ( scale_potential != 0 and 
         not (mol.isEmpty()) and
         not spce->beyond(switchfunc->cutoffDistance(),
                          mol.aaBox(), pots.grid().aaBox()) )
    {
        this->_pvt_calculateLJPotential(mol, probe, pots,
                                             workspace, scale_potential);
    }
}

SIRE_ALWAYS_INLINE void 
InterLJPotential::calculatePotential(
                                const InterLJPotential::Molecule &mol, 
                                const LJProbe &probe,
                                GridPotentialTable &pots, 
                                InterLJPotential::PotentialWorkspace &workspace,
                                double scale_potential) const
{
    this->calculateLJPotential(mol, probe, pots, workspace, scale_potential);
}

SIRE_ALWAYS_INLINE void 
InterLJPotential::calculatePotential(
                                const InterLJPotential::Molecule &mol,
                                const LJProbe &probe,
                                GridPotentialTable &pots,
                                const Symbol &symbol,
                                const InterLJPotential::Components &components,
                                InterLJPotential::PotentialWorkspace &workspace,
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
SIRE_ALWAYS_INLINE void 
IntraLJPotential::calculateForce(const IntraLJPotential::Molecule &mol, 
                                 MolForceTable &forces,
                                 const Symbol &symbol,
                                 const IntraLJPotential::Components &components,
                                 IntraLJPotential::ForceWorkspace &workspace,
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
SIRE_ALWAYS_INLINE void
IntraLJPotential::calculateForce(const IntraLJPotential::Molecule &mol,
                                 const IntraLJPotential::Molecule &rest_of_mol,
                                 MolForceTable &forces,
                                 const Symbol &symbol,
                                 const IntraLJPotential::Components &components,
                                 IntraLJPotential::ForceWorkspace &workspace,
                                 double scale_force) const
{
    if (symbol == components.total())
        this->calculateForce(mol, rest_of_mol, forces, workspace, scale_force);
        
    else
        throwMissingForceComponent(symbol, components);
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::detail::LJParamID&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::detail::LJParamID&);

template<class LJPot>
QDataStream& operator<<(QDataStream &ds,
                        const SireMM::LJPotentialInterface<LJPot> &ljpot)
{
    ds << static_cast<const LJPot&>(ljpot);
    return ds;
}

template<class LJPot>
QDataStream& operator>>(QDataStream &ds,
                        SireMM::LJPotentialInterface<LJPot> &ljpot)
{
    ds >> static_cast<LJPot&>(ljpot);
    return ds;
}

Q_DECLARE_TYPEINFO( SireMM::detail::LJParamID, Q_MOVABLE_TYPE );

SIRE_EXPOSE_CLASS( SireMM::LJParameterName )
SIRE_EXPOSE_CLASS( SireMM::LJParameterName3D )
SIRE_EXPOSE_CLASS( SireMM::ScaledLJParameterNames3D )

SIRE_END_HEADER

#endif
