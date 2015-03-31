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

#include "SireBase/sparsematrix.hpp"

#include "cljpotential.h"
#include "ljparameter.h"
#include "switchingfunction.h"

#include "SireMol/mover.hpp"
#include "SireMol/atomcoords.h"

#include "SireVol/cartesian.h"

#include "SireBase/countflops.h"

#include "SireMaths/maths.h"

#include "SireUnits/units.h"

#include "SireBase/errors.h"
#include "SireError/errors.h"
#include "SireFF/errors.h"

#include "SireStream/datastream.h"

//#undef SIRE_USE_SSE // for some reason, this makes things slower...

#ifdef SIRE_USE_SSE
    #ifdef __SSE__
        #include <emmintrin.h>   // CONDITIONAL_INCLUDE
    #else
        #undef SIRE_USE_SSE
    #endif
#endif

#include <QDebug>

using namespace SireMM;
using namespace SireMM::detail;

using namespace SireFF;
using namespace SireFF::detail;

using namespace SireUnits;

using namespace SireMol;
using namespace SireVol;

using namespace SireMaths;

using namespace SireBase;

using namespace SireStream;

///////
/////// Static data for IntraScaleParameterName
///////

QString IntraScaleParameterName::nbscl_param( "intrascale" );

///////
/////// Completely instantiate the CLJPotential ancillary classes
///////

namespace SireMM
{
    namespace detail
    {
        template
        class IntraScaledParameters<CLJNBPairs>;

        template
        class IntraScaledAtomicParameters< AtomicParameters3D<CLJParameter>,
                                           IntraScaledParameters<CLJNBPairs> >;
    }
}

namespace SireFF
{
    namespace detail
    {
        template
        class AtomicParameters3D<CLJParameter>;

        template
        class FFMolecule3D<InterCLJPotential>;

        template
        class FFMolecules3D<InterCLJPotential>;

        template
        class ChangedMolecule<InterCLJPotential::Molecule>;

        template
        class FFMolecule3D<IntraCLJPotential>;

        template
        class FFMolecules3D<IntraCLJPotential>;

        template
        class ChangedMolecule<IntraCLJPotential::Molecule>;
    }
}

/** Streaming functions for CLJParameter */
QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, 
                                      const SireMM::detail::CLJParameter &cljparam)
{
    ds << static_cast<const SireMM::detail::ChargeParameter&>(cljparam)
       << static_cast<const SireMM::detail::LJParamID&>(cljparam);
    
    return ds;
}

QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, 
                                      SireMM::detail::CLJParameter &cljparam)
{
    ds >> static_cast<SireMM::detail::ChargeParameter&>(cljparam)
       >> static_cast<SireMM::detail::LJParamID&>(cljparam);
    
    return ds;
}

/** Internal function used to get the charge and LJ parameters from a molecule
    and convert them into a PackedArray of reduced charges and LJ parameter IDs */
PackedArray2D<CLJParameter> 
CLJPotential::getCLJParameters(const PartialMolecule &molecule,
                               const PropertyName &charge_property,
                               const PropertyName &lj_property)
{
    const AtomCharges &chgs = molecule.property(charge_property).asA<AtomCharges>();

    const AtomLJs &ljs = molecule.property(lj_property).asA<AtomLJs>();
    
    const AtomSelection &selected_atoms = molecule.selection();
    
    if (selected_atoms.selectedNone())
        return PackedArray2D<CLJParameter>();
    
    //create space for the parameters - only need space for CutGroups
    //that contain at least one selected atom
    QVector< QVector<CLJParameter> > cljparams( selected_atoms.nSelectedCutGroups() );
    QVector<CLJParameter>* cljparams_array = cljparams.data();

    const double sqrt_4pieps0 = std::sqrt(SireUnits::one_over_four_pi_eps0);

    try
    {

    LJParameterDB::lock();

    if (selected_atoms.selectedAllCutGroups())
    {
        const int ncg = molecule.data().info().nCutGroups();
    
        for (CGIdx i(0); i<ncg; ++i)
        {
            const int nats = molecule.data().info().nAtoms(i);
            
            QVector<CLJParameter> group_cljs(nats);
            CLJParameter *group_cljs_array = group_cljs.data();
            
            //get the arrays containing the charge and LJ parameters
            //for this CutGroup
            const SireUnits::Dimension::Charge *group_chgs = chgs.constData(i);
            const LJParameter *group_ljs = ljs.constData(i);
            
            if (selected_atoms.selectedAll(i))
            {
                for (Index j(0); j<nats; ++j)
                {
                    group_cljs_array[j].reduced_charge = group_chgs[j] * sqrt_4pieps0;
                    group_cljs_array[j].ljid = 
                            LJParameterDB::_locked_addLJParameter(group_ljs[j]);
                }
            }
            else
            {
                foreach (Index j, selected_atoms.selectedAtoms(i))
                {
                    group_cljs_array[j].reduced_charge = group_chgs[j] * sqrt_4pieps0;
                    group_cljs_array[j].ljid =
                            LJParameterDB::_locked_addLJParameter(group_ljs[j]);
                }
            }
            
            cljparams_array[i] = group_cljs;
        }
    }
    else
    {
        int igroup = 0;
    
        foreach (CGIdx i, selected_atoms.selectedCutGroups())
        {
            const int nats = molecule.data().info().nAtoms(i);
            
            QVector<CLJParameter> group_cljs(nats);
            CLJParameter *group_cljs_array = group_cljs.data();
            
            //get the arrays containing the charge and LJ parameters
            //for this CutGroup
            const SireUnits::Dimension::Charge *group_chgs = chgs.constData(i);
            const LJParameter *group_ljs = ljs.constData(i);
            
            if (selected_atoms.selectedAll(i))
            {
                for (Index j(0); j<nats; ++j)
                {
                    group_cljs_array[j].reduced_charge = group_chgs[j] * sqrt_4pieps0;
                    group_cljs_array[j].ljid = 
                            LJParameterDB::_locked_addLJParameter(group_ljs[j]);
                }
            }
            else
            {
                foreach (Index j, selected_atoms.selectedAtoms(i))
                {
                    group_cljs_array[j].reduced_charge = group_chgs[j] * sqrt_4pieps0;
                    group_cljs_array[j].ljid =
                            LJParameterDB::_locked_addLJParameter(group_ljs[j]);
                }
            }
            
            cljparams_array[igroup] = group_cljs;
            igroup += 1;
        }
    }
    
    LJParameterDB::unlock();
    
    }
    catch(...)
    {
        LJParameterDB::unlock();
        throw;
    }
    
    return PackedArray2D<CLJParameter>( cljparams );
}

/////////////
///////////// Implementation of CLJPotential
/////////////

static const RegisterMetaType<CLJPotential> r_cljpot( MAGIC_ONLY, NO_ROOT,
                                                      "SireMM::CLJPotential" );

/** Serialise to a binary datastream */
QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds,
                                      const CLJPotential &cljpot)
{
    writeHeader(ds, r_cljpot, 1);
    
    ds << cljpot.props;

    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds,
                                      CLJPotential &cljpot)
{
    VersionID v = readHeader(ds, r_cljpot);
    
    if (v == 1)
    {
        ds >> cljpot.props;
    
        //extract all of the properties
        cljpot.spce = cljpot.props.property("space").asA<Space>();
        cljpot.switchfunc = cljpot.props.property("switchingFunction")
                                        .asA<SwitchingFunction>();
    
        cljpot.rf_dielectric_constant = cljpot.props.property("reactionFieldDielectric")
                                            .asA<VariantProperty>().convertTo<double>();
                                            
        cljpot.use_reaction_field = cljpot.props.property("useReactionField")
                                        .asA<VariantProperty>().convertTo<bool>();
    
        cljpot.combining_rules = LJParameterDB::interpret(
                                    cljpot.props.property("combiningRules")
                                        .asA<VariantProperty>().convertTo<QString>() );

        cljpot.use_electrostatic_shifting = cljpot.props.property("shiftElectrostatics")
                                        .asA<VariantProperty>().convertTo<bool>();
        
        cljpot.use_atomistic_cutoff = cljpot.props.property("useAtomisticCutoff")
                                        .asA<VariantProperty>().convertTo<bool>();
        
        cljpot.need_update_ljpairs = true;
    }
    else 
        throw version_error(v, "1", r_cljpot, CODELOC);
    
    return ds;
}

/** Constructor */
CLJPotential::CLJPotential()
             : spce( Space::null() ), switchfunc( SwitchingFunction::null() ),
               combining_rules( LJParameterDB::interpret("arithmetic") ),
               rf_dielectric_constant(1), use_reaction_field(false),
               need_update_ljpairs(true), use_electrostatic_shifting(false),
               use_atomistic_cutoff(false)
{
    //record the defaults
    props.setProperty( "space", spce );
    props.setProperty( "switchingFunction", switchfunc );
    props.setProperty( "combiningRules", 
                       VariantProperty( LJParameterDB::toString(combining_rules) ) );
    props.setProperty( "shiftElectrostatics",
                       VariantProperty(use_electrostatic_shifting) );
    props.setProperty("useReactionField", VariantProperty(use_reaction_field));
    props.setProperty("reactionFieldDielectric", VariantProperty(rf_dielectric_constant));
    props.setProperty("useAtomisticCutoff", VariantProperty(use_atomistic_cutoff));
    props.setProperty("useGroupCutoff", VariantProperty(true));
}

/** Copy constructor */
CLJPotential::CLJPotential(const CLJPotential &other)
             : ljpairs(other.ljpairs), props(other.props),
               spce(other.spce), switchfunc(other.switchfunc),
               combining_rules(other.combining_rules),
               rf_dielectric_constant(other.rf_dielectric_constant),
               use_reaction_field(other.use_reaction_field),
               need_update_ljpairs(other.need_update_ljpairs),
               use_electrostatic_shifting(other.use_electrostatic_shifting),
               use_atomistic_cutoff(other.use_atomistic_cutoff)
{}

/** Destructor */
CLJPotential::~CLJPotential()
{}

/** Copy assignment operator */
CLJPotential& CLJPotential::operator=(const CLJPotential &other)
{
    if (this != &other)
    {
        ljpairs = other.ljpairs;
        props = other.props;
        spce = other.spce;
        switchfunc = other.switchfunc;
        combining_rules = other.combining_rules;
        rf_dielectric_constant = other.rf_dielectric_constant;
        use_reaction_field = other.use_reaction_field;
        need_update_ljpairs = other.need_update_ljpairs;
        use_electrostatic_shifting = other.use_electrostatic_shifting;
        use_atomistic_cutoff = other.use_atomistic_cutoff;
    }
    
    return *this;
}

/** You need to call this function before you start a block of 
    energy of force evaluation using this forcefield. You should
    also call 'finishedEvaluation()' once you have finished. */
void CLJPotential::startEvaluation()
{
    if (need_update_ljpairs)
    {
        //get the LJPairs array from the LJParameterDB
        ljpairs = LJParameterDB::getLJPairs(combining_rules);
        need_update_ljpairs = false;
    }
}

/** You should call this function once you have finished a block of
    force or energy evaluation using this potential */
void CLJPotential::finishedEvaluation()
{}

/** Return all of the properties set in this forcefield */
const Properties& CLJPotential::properties() const
{
    return props;
}

/** Return whether or not this potential has a property called 'name' */
bool CLJPotential::containsProperty(const QString &name) const
{
    return props.hasProperty(name);
}

/** Return the property with name 'name'

    \throw SireBase::missing_property
*/
const Property& CLJPotential::property(const QString &name) const
{
    return props.property(name);
}

/** Set the 3D space in which the molecules in this potential are evaluated */
bool CLJPotential::setSpace(const Space &new_space)
{
    if ( not spce->equals(new_space) )
    {
        spce = new_space;
        props.setProperty( "space", new_space );
        this->changedPotential();
        return true;
    }
    else
        return false;
}

/** Set the switching function used to scale the interactions between
    CutGroups to zero at the cutoff */
bool CLJPotential::setSwitchingFunction(const SwitchingFunction &new_switchfunc)
{
    if ( not switchfunc->equals(new_switchfunc) )
    {
        switchfunc = new_switchfunc;
        props.setProperty( "switchingFunction", new_switchfunc );
        
        this->changedPotential();
        return true;
    }
    else
        return false;
}

/** Return whether or not an atomistic cutoff is being used */
bool CLJPotential::useAtomisticCutoff() const
{
    return use_atomistic_cutoff;
}

/** Return whether or not a group-based cutoff is being used. This
    is the default cutoff method of the CLJ potentials */
bool CLJPotential::useGroupCutoff() const
{
    return not (use_atomistic_cutoff or use_reaction_field or use_electrostatic_shifting);
}

/** Switch the use of the group-based cutoff. If turning on the cutoff, this will
    disable the atomistic, electrostatic shifting and reaction field cutoffs. If turning
    off this cutoff, this will switch over to an atomistic cutoff */
bool CLJPotential::setUseGroupCutoff(bool switchgroup)
{
    if (useGroupCutoff() != switchgroup)
    {
        if (switchgroup)
        {
            use_atomistic_cutoff = false;
            use_reaction_field = false;
            rf_dielectric_constant = 1;
            use_electrostatic_shifting = false;
            
            props.setProperty( "shiftElectrostatics",
                               VariantProperty(use_electrostatic_shifting) );
            props.setProperty("useReactionField", VariantProperty(use_reaction_field));
            props.setProperty("reactionFieldDielectric", VariantProperty(rf_dielectric_constant));
            props.setProperty("useAtomisticCutoff", VariantProperty(use_atomistic_cutoff));
            props.setProperty("useGroupCutoff", VariantProperty(true));
        }
        else
            setUseAtomisticCutoff(true);
        
        this->changedPotential();
        return true;
    }
    else
        return false;
}

/** Switch the use of an atomistic cutoff. If turning on the cutoff, this will disable
    the reaction field, electrostatic shifting and group based cutoff methods.
    If turning off the cutoff, this will switch over to the group based cutoff */
bool CLJPotential::setUseAtomisticCutoff(bool switchatomistic)
{
    if (use_atomistic_cutoff != switchatomistic)
    {
        if (switchatomistic)
        {
            use_atomistic_cutoff = true;
            use_reaction_field = false;
            rf_dielectric_constant = 1;
            use_electrostatic_shifting = false;
            
            props.setProperty( "shiftElectrostatics",
                               VariantProperty(use_electrostatic_shifting) );
            props.setProperty("useReactionField", VariantProperty(use_reaction_field));
            props.setProperty("reactionFieldDielectric", VariantProperty(rf_dielectric_constant));
            props.setProperty("useAtomisticCutoff", VariantProperty(use_atomistic_cutoff));
            props.setProperty("useGroupCutoff", VariantProperty(false));
        }
        else
        {
            //switch back to the group-based cutoff
            setUseGroupCutoff(true);
        }
        
        this->changedPotential();
        
        return true;
    }
    else
        return false;
}

/** Set whether or not to shift the electrostatics between CutGroups so that
    the group-group net electrostatic interaction energy between CutGroups
    is zero at the cutoff */
bool CLJPotential::setShiftElectrostatics(bool switchelectro)
{
    if (use_electrostatic_shifting != switchelectro)
    {
        if (switchelectro)
        {
            use_atomistic_cutoff = false;
            use_reaction_field = false;
            rf_dielectric_constant = 1;
            use_electrostatic_shifting = true;
            
            props.setProperty( "shiftElectrostatics",
                               VariantProperty(use_electrostatic_shifting) );
            props.setProperty("useReactionField", VariantProperty(use_reaction_field));
            props.setProperty("reactionFieldDielectric", VariantProperty(rf_dielectric_constant));
            props.setProperty("useAtomisticCutoff", VariantProperty(use_atomistic_cutoff));
            props.setProperty("useGroupCutoff", VariantProperty(false));
        }
        else
        {
            //switch back to the group-based cutoff
            setUseGroupCutoff(true);
        }
        
        this->changedPotential();
        
        return true;
    }
    else
        return false;
}

/** Set the combining rules to use to obtain mixed LJ parameters */
bool CLJPotential::setCombiningRules(const QString &combiningrules)
{
    LJParameterDB::CombiningRules new_rules = LJParameterDB::interpret(combiningrules);
    
    if (new_rules != combining_rules)
    {
        combining_rules = new_rules;
        need_update_ljpairs = true;
        props.setProperty( "combiningRules", VariantProperty(combiningrules) );
        this->changedPotential();
        return true;
    }
    else
        return false;
}

/** Set the property 'name' to the value 'value'. Returns whether or not
    this changes this forcefield.

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
bool CLJPotential::setProperty(const QString &name, const Property &value)
{
    if (name == QLatin1String("space"))
    {
        return this->setSpace( value.asA<Space>() );
    }
    else if (name == QLatin1String("switchingFunction"))
    {
        return this->setSwitchingFunction( value.asA<SwitchingFunction>() );
    }
    else if (name == QLatin1String("shiftElectrostatics"))
    {
        return this->setShiftElectrostatics( value.asA<VariantProperty>()
                                                     .convertTo<bool>() );
    }
    else if (name == QLatin1String("useReactionField"))
    {
        return this->setUseReactionField( value.asA<VariantProperty>()
                                                    .convertTo<bool>() );
    }
    else if (name == QLatin1String("reactionFieldDielectric"))
    {
        return this->setReactionFieldDielectric( value.asA<VariantProperty>()
                                                        .convertTo<double>() );
    }
    else if (name == QLatin1String("combiningRules"))
    {
        return this->setCombiningRules( value.asA<VariantProperty>()
                                                     .convertTo<QString>() );
    }
    else
        throw SireBase::missing_property( QObject::tr(
            "The CLJ potentials do not have a property called \"%1\" that "
            "can be changed. Available properties are [ %2 ].")
                .arg(name, QStringList(props.propertyKeys()).join(", ")), CODELOC );
                
    return false;
}

/** Return the 3D space in which this potential is evaluated */
const Space& CLJPotential::space() const
{
    return *spce;
}

/** Return the switching function used to scale the group-group
    interactions to zero at the cutoff */
const SwitchingFunction& CLJPotential::switchingFunction() const
{
    return *switchfunc;
}

/** Return whether or not the net group-group electrostatic interaction
    energy is shifted so that it is zero at the cutoff */
bool CLJPotential::shiftElectrostatics() const
{
    return use_electrostatic_shifting;
}

/** Return the string identifying the combining rules used to 
    obtain the mixed LJ parameters */
const QString& CLJPotential::combiningRules() const
{
    return LJParameterDB::toString(combining_rules);
}

/** Turn on the use of the reaction field method for handling the cutoff */
bool CLJPotential::setUseReactionField(bool switchrf)
{
    if (use_reaction_field != switchrf)
    {
        if (switchrf)
        {
            use_atomistic_cutoff = false;
            use_reaction_field = true;
            use_electrostatic_shifting = false;
            
            props.setProperty( "shiftElectrostatics",
                               VariantProperty(use_electrostatic_shifting) );
            props.setProperty("useReactionField", VariantProperty(use_reaction_field));
            props.setProperty("reactionFieldDielectric", VariantProperty(rf_dielectric_constant));
            props.setProperty("useAtomisticCutoff", VariantProperty(use_atomistic_cutoff));
            props.setProperty("useGroupCutoff", VariantProperty(false));
        }
        else
        {
            //switch back to the group-based cutoff
            setUseGroupCutoff(true);
        }
        
        this->changedPotential();
        
        return true;
    }
    else
        return false;
}

/** Return whether or not the reaction field method is being used */
bool CLJPotential::useReactionField() const
{
    return use_reaction_field;
}

/** Set the dielectric constant to use for the reaction field. Note that this
    only has an effect if the reaction field cutoff method is being used */
bool CLJPotential::setReactionFieldDielectric(double dielectric)
{
    if (dielectric != rf_dielectric_constant)
    {
        rf_dielectric_constant = dielectric;
        props.setProperty("reactionFieldDielectric", VariantProperty(dielectric) );
        
        if (use_reaction_field)
            this->changedPotential();
            
        return true;
    }
    else
        return false;
}

/** Return the dielectric constant used for the reaction field cutoff. This 
    only has an effect if the reaction field cutoff method is being used */
double CLJPotential::reactionFieldDielectric() const
{
    return rf_dielectric_constant;
}

/////////////
///////////// Implementation of InterCLJPotential
/////////////

static const RegisterMetaType<InterCLJPotential> r_interclj( MAGIC_ONLY, NO_ROOT,
                                            InterCLJPotential::typeName() );

/** Serialise to a binary datastream */
QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds,
                                      const InterCLJPotential &interclj)
{
    writeHeader(ds, r_interclj, 1);
    
    ds << static_cast<const CLJPotential&>(interclj);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds,
                                      InterCLJPotential &interclj)
{
    VersionID v = readHeader(ds, r_interclj);
    
    if (v == 1)
    {
        ds >> static_cast<CLJPotential&>(interclj);
    }
    else
        throw version_error(v, "1", r_interclj, CODELOC);
        
    return ds;
}

/** Constructor */
InterCLJPotential::InterCLJPotential() : CLJPotential()
{}

/** Copy constructor */
InterCLJPotential::InterCLJPotential(const InterCLJPotential &other)
                  : CLJPotential(other)
{}

/** Destructor */
InterCLJPotential::~InterCLJPotential()
{}

/** Copy assignment operator */
InterCLJPotential& InterCLJPotential::operator=(const InterCLJPotential &other)
{
    CLJPotential::operator=(other);
    return *this;
}

void InterCLJPotential::throwMissingForceComponent(const Symbol &symbol,
                              const InterCLJPotential::Components &components) const
{
    throw SireFF::missing_component( QObject::tr(
        "There is no force component in potential %1 - available "
        "components are %2, %3 and %4.")
            .arg(this->what())
            .arg(components.total().toString(), components.coulomb().toString(),
                 components.lj().toString()), CODELOC );
}

void InterCLJPotential::throwMissingFieldComponent(const Symbol &symbol,
                              const InterCLJPotential::Components &components) const
{
    throw SireFF::missing_component( QObject::tr(
        "There is no field component in potential %1 - available "
        "components are %2, %3 and %4.")
            .arg(this->what())
            .arg(components.total().toString(), components.coulomb().toString(),
                 components.lj().toString()), CODELOC );
}

void InterCLJPotential::throwMissingPotentialComponent(const Symbol &symbol,
                              const InterCLJPotential::Components &components) const
{
    throw SireFF::missing_component( QObject::tr(
        "There is no potential component in potential %1 - available "
        "components are %2, %3 and %4.")
            .arg(this->what())
            .arg(components.total().toString(), components.coulomb().toString(),
                 components.lj().toString()), CODELOC );
}

/** Return all of the parameters needed by this potential for 
    the molecule 'molecule', using the supplied property map to
    find the properties that contain those parameters
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
InterCLJPotential::Parameters 
InterCLJPotential::getParameters(const PartialMolecule &molecule,
                                 const PropertyMap &map)
{
    need_update_ljpairs = true;

    return Parameters( molecule, map[parameters().coordinates()],
                       getCLJParameters(molecule, map[parameters().charge()],
                                        map[parameters().lj()]) );
}

/** Update the parameters for the molecule going from 'old_molecule' to 
    'new_molecule', with the parameters found using the property map 'map'
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
InterCLJPotential::Parameters
InterCLJPotential::updateParameters(const InterCLJPotential::Parameters &old_params,
                                    const PartialMolecule &old_molecule,
                                    const PartialMolecule &new_molecule,
                                    const PropertyMap &map)
{
    if (old_molecule.selection() != new_molecule.selection())
        //the selection has changed - just get completely new parameters
        return this->getParameters(new_molecule, map);

    Parameters new_params = old_params;

    //get the property names
    const PropertyName &coords_property = map[parameters().coordinates()];
    const PropertyName &chg_property = map[parameters().charge()];
    const PropertyName &lj_property = map[parameters().lj()];
    
    //get what has changed
    bool new_coords = old_molecule.version(coords_property) !=
                         new_molecule.version(coords_property);
                             
    bool new_clj = ( old_molecule.version(chg_property) !=
                         new_molecule.version(chg_property) ) or
                   ( old_molecule.version(lj_property) !=
                         new_molecule.version(lj_property) );

    if (new_coords)
    {
        new_params.setAtomicCoordinates( AtomicCoords3D(new_molecule, 
                                                        coords_property) );
    }

    if (new_clj)
    {
        new_params.setAtomicParameters( getCLJParameters(new_molecule,
                                                         chg_property, lj_property) );
        need_update_ljpairs = true;
    }

    return new_params;
}
                 
/** Update the parameters for the molecule going from 'old_molecule' to 
    'new_molecule', also while the parameters of 'old_molecule'
    where found in 'old_map', now get the parameters using 'new_map'
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
InterCLJPotential::Parameters
InterCLJPotential::updateParameters(const InterCLJPotential::Parameters &old_params,
                                    const PartialMolecule &old_molecule,
                                    const PartialMolecule &new_molecule,
                                    const PropertyMap &old_map, 
                                    const PropertyMap &new_map)
{
    if (old_molecule.selection() != new_molecule.selection())
        //the selection has changed - just get completely new parameters
        return this->getParameters(new_molecule, new_map);

    Parameters new_params = old_params;

    //get the property names
    const PropertyName &old_coords = old_map[parameters().coordinates()];
    const PropertyName &old_chg = old_map[parameters().charge()];
    const PropertyName &old_lj = old_map[parameters().lj()];
    
    const PropertyName &new_coords = new_map[parameters().coordinates()];
    const PropertyName &new_chg = new_map[parameters().charge()];
    const PropertyName &new_lj = new_map[parameters().lj()];
    
    //get what has changed
    bool changed_coords = (new_coords != old_coords) or
                           old_molecule.version(old_coords) !=
                           new_molecule.version(old_coords);
                             
    bool changed_clj = (new_chg != old_chg or new_lj != old_lj) or
                       ( old_molecule.version(old_chg) !=
                         new_molecule.version(old_chg) ) or
                       ( old_molecule.version(old_lj) !=
                         new_molecule.version(old_lj) );

    if (changed_coords)
        new_params.setAtomicCoordinates( AtomicCoords3D(new_molecule, 
                                                        new_coords) );

    if (changed_clj)
    {
        new_params.setAtomicParameters( getCLJParameters(new_molecule,
                                                         new_chg, new_lj) );
        
        need_update_ljpairs = true;
    }

    return new_params;
}

/** Return the InterCLJPotential::Molecule representation of 'molecule',
    using the supplied PropertyMap to find the properties that contain
    the necessary forcefield parameters
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
InterCLJPotential::Molecule
InterCLJPotential::parameterise(const PartialMolecule &molecule,
                                const PropertyMap &map)
{
    return InterCLJPotential::Molecule(molecule, *this, map);
}

/** Convert the passed group of molecules into InterCLJPotential::Molecules,
    using the supplied PropertyMap to find the properties that contain
    the necessary forcefield parameters in each molecule
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
InterCLJPotential::Molecules 
InterCLJPotential::parameterise(const MoleculeGroup &molecules,
                                const PropertyMap &map)
{
    return InterCLJPotential::Molecules(molecules, *this, map);
}

/** Return the total charge of the parameters for the group in 'params' */
double InterCLJPotential::totalCharge(
                        const InterCLJPotential::Parameters::Array &params) const
{
    int n = params.count();
    const Parameter *params_array = params.constData();
    
    double chg = 0;
    
    for (int i=0; i<n; ++i)
    {
        chg += params_array[i].reduced_charge;
    }
    
    return chg;
}

/** Calculate the coulomb and LJ energy between the passed pair
    of molecules and add these energies onto 'energy'. This uses
    the passed workspace to perform the calculation */
void InterCLJPotential::_pvt_calculateEnergy(const InterCLJPotential::Molecule &mol0,
                                             const InterCLJPotential::Molecule &mol1,
                                             InterCLJPotential::Energy &energy,
                                             InterCLJPotential::EnergyWorkspace &distmat,
                                             double scale_energy) const
{
    const quint32 ngroups0 = mol0.coordinates().count();
    const quint32 ngroups1 = mol1.coordinates().count();
    
    if (ngroups0 < ngroups1)
    {
        //this is faster if the larger molecule is on the outer loop
        this->_pvt_calculateEnergy(mol1, mol0, energy, distmat, scale_energy);
        return;
    }
    
    const CoordGroup *groups0_array = mol0.coordinates().constData();
    const CoordGroup *groups1_array = mol1.coordinates().constData();
    
    BOOST_ASSERT( mol0.parameters().atomicParameters().count() == int(ngroups0) );
    BOOST_ASSERT( mol1.parameters().atomicParameters().count() == int(ngroups1) );
    
    const Parameters::Array *molparams0_array 
                                = mol0.parameters().atomicParameters().constData();
    const Parameters::Array *molparams1_array 
                                = mol1.parameters().atomicParameters().constData();
    
    double cnrg = 0;
    double ljnrg = 0;
    
    #ifdef SIRE_TIME_ROUTINES
    int nflops = 0;
    #endif

    const double Rcoul = qMax(1e-5,qMin(1e9,
                            switchfunc->electrostaticCutoffDistance().to(angstrom)));
    const double Rlj = qMax(1e-5,qMin(1e9,switchfunc->vdwCutoffDistance().to(angstrom)) );
    const double Rc = qMax(Rcoul,Rlj);

    //loop over all pairs of CutGroups in the two molecules
    if (use_electrostatic_shifting)
    {
        //we use the force shifted coulomb energy described
        //in Fennell and Gezelter, J. Chem. Phys., 124, 234104, 2006
        //We use alpha=0, as I have seen that a 25 A cutoff gives stable results
        //with alpha=0, and this way we avoid changing the hamiltonian significantly
        //by having an erfc function
        const double one_over_Rcoul = double(1) / Rcoul;
        const double one_over_Rcoul2 = double(1) / (Rcoul*Rcoul);
        
        for (quint32 igroup=0; igroup<ngroups0; ++igroup)
        {
            const Parameters::Array &params0 = molparams0_array[igroup];
            const CoordGroup &group0 = groups0_array[igroup];
            const AABox &aabox0 = group0.aaBox();
            const quint32 nats0 = group0.count();
            const Parameter *params0_array = params0.constData();
            
            for (quint32 jgroup=0; jgroup<ngroups1; ++jgroup)
            {
                const CoordGroup &group1 = groups1_array[jgroup];
                const Parameters::Array &params1 = molparams1_array[jgroup];
                
                //check first that these two CoordGroups could be within cutoff
                //(if there is only one CutGroup in both molecules then this
                //test has already been performed and passed)
                const bool within_cutoff = (ngroups0 == 1 and ngroups1 == 1) or not
                                            spce->beyond(Rc, aabox0, group1.aaBox());
                
                if (not within_cutoff)
                    //this CutGroup is beyond the cutoff distance
                    continue;
                
                //calculate all of the interatomic distances
                const double mindist = spce->calcDist(group0, group1, distmat);
                
                if (mindist > Rc)
                {
                    //all of the atoms are definitely beyond cutoff
                    continue;
                }
                   
                double icnrg = 0;
                double iljnrg = 0;
                
                //loop over all interatomic pairs and calculate the energies
                const quint32 nats1 = group1.count();
                const Parameter *params1_array = params1.constData();

                #ifdef SIRE_USE_SSE
                {
                    const int remainder = nats1 % 2;
                    
                    __m128d sse_cnrg = { 0, 0 };
                    __m128d sse_ljnrg = { 0, 0 };

                    const __m128d sse_one = { 1.0, 1.0 };

                    const __m128d sse_Rcoul = _mm_set1_pd(Rcoul);
                    const __m128d sse_one_over_Rcoul = _mm_set1_pd(one_over_Rcoul);
                    const __m128d sse_one_over_Rcoul2 = _mm_set1_pd(one_over_Rcoul2);

                    const __m128d sse_Rlj = _mm_set1_pd(Rlj);
                    
                    for (quint32 i=0; i<nats0; ++i)
                    {
                        distmat.setOuterIndex(i);
                        const Parameter &param0 = params0_array[i];
                        
                        __m128d sse_chg0 = _mm_set_pd( param0.reduced_charge, 
                                                       param0.reduced_charge );
                                             
                        //process atoms in pairs (so can then use SSE)
                        for (quint32 j=0; j<nats1-1; j += 2)
                        {
                            const Parameter &param10 = params1_array[j];
                            const Parameter &param11 = params1_array[j+1];
                            
                            const __m128d sse_r = _mm_set_pd( distmat[j], distmat[j+1] );
                            const __m128d sse_one_over_r = _mm_div_pd(sse_one, sse_r);

                            //coulomb calculation
                            {
                                const __m128d sse_in_cutoff = _mm_cmplt_pd(sse_r, sse_Rcoul);

                                __m128d nrg = _mm_sub_pd(sse_r, sse_Rcoul);
                                nrg = _mm_mul_pd(nrg, sse_one_over_Rcoul2);
                                nrg = _mm_add_pd(nrg, sse_one_over_r);
                                nrg = _mm_sub_pd(nrg, sse_one_over_Rcoul);

                                __m128d sse_chg = _mm_set_pd( param10.reduced_charge,
                                                              param11.reduced_charge );
                        
                                sse_chg = _mm_mul_pd(sse_chg, sse_chg0);
                        
                                nrg = _mm_mul_pd(sse_chg, nrg);
                                nrg = _mm_and_pd(nrg, sse_in_cutoff);
                                
                                sse_cnrg = _mm_add_pd(sse_cnrg, nrg);
                            }
                            
                            //LJ calculation
                            {
                                const __m128d sse_in_cutoff = _mm_cmplt_pd(sse_r, sse_Rlj);

                                const LJPair &ljpair0 = ljpairs.constData()[
                                                            ljpairs.map(param0.ljid,
                                                                        param10.ljid)];
                        
                                const LJPair &ljpair1 = ljpairs.constData()[
                                                            ljpairs.map(param0.ljid,
                                                                        param11.ljid)];
                        
                                __m128d sse_sig = _mm_set_pd( ljpair0.sigma(), ljpair1.sigma() );
                                __m128d sse_eps = _mm_set_pd( ljpair0.epsilon(),
                                                              ljpair1.epsilon() );
                                                        
                                //calculate (sigma/r)^6 and (sigma/r)^12
                                __m128d sse_sig_over_dist2 = _mm_mul_pd(sse_sig, sse_one_over_r);
                                sse_sig_over_dist2 = _mm_mul_pd( sse_sig_over_dist2,
                                                                 sse_sig_over_dist2 );
                                                         
                                __m128d sse_sig_over_dist6 = _mm_mul_pd(sse_sig_over_dist2,
                                                                        sse_sig_over_dist2);
                                                            
                                sse_sig_over_dist6 = _mm_mul_pd(sse_sig_over_dist6,
                                                                sse_sig_over_dist2);
                                                         
                                __m128d sse_sig_over_dist12 = _mm_mul_pd(sse_sig_over_dist6,
                                                                         sse_sig_over_dist6);
                                                  
                                __m128d nrg = _mm_sub_pd(sse_sig_over_dist12,
                                                         sse_sig_over_dist6);
                                                     
                                nrg = _mm_mul_pd(nrg, sse_eps);
                                nrg = _mm_and_pd(nrg, sse_in_cutoff);
                                sse_ljnrg = _mm_add_pd(sse_ljnrg, nrg);
                            }
                        }
                              
                        if (remainder == 1)
                        {
                            const Parameter &param1 = params1_array[nats1-1];

                            const double r = distmat[nats1-1];

                            if (r < Rc)
                            {
                                const double one_over_r = double(1) / r;

                                if (r < Rcoul)
                                {
                                    icnrg += param0.reduced_charge * param1.reduced_charge *
                                         (one_over_r - one_over_Rcoul + one_over_Rcoul2*(r-Rcoul));
                                }

                                if (r < Rlj)
                                {
                                    const LJPair &ljpair = ljpairs.constData()[
                                                            ljpairs.map(param0.ljid,
                                                                        param1.ljid)];

                                    double sig_over_dist6 = pow_6(ljpair.sigma()*one_over_r);
                                    double sig_over_dist12 = pow_2(sig_over_dist6);
            
                                    iljnrg += ljpair.epsilon() * (sig_over_dist12 -
                                                                  sig_over_dist6);
                                }
                            }
                        }
                    }
                    
                    icnrg += *((const double*)&sse_cnrg) +
                             *( ((const double*)&sse_cnrg) + 1 );
                             
                    iljnrg += *((const double*)&sse_ljnrg) +
                              *( ((const double*)&sse_ljnrg) + 1 );
                }
                #else
                {
                    for (quint32 i=0; i<nats0; ++i)
                    {
                        distmat.setOuterIndex(i);
                        const Parameter &param0 = params0_array[i];
                    
                        for (quint32 j=0; j<nats1; ++j)
                        {
                            const Parameter &param1 = params1_array[j];

                            const double r = distmat[j];
                            
                            if (r < Rc)
                            {
                                const double one_over_r = double(1) / r;
                            
                                if (r < Rcoul)
                                {
                                    icnrg += param0.reduced_charge * param1.reduced_charge *
                                        (one_over_r - one_over_Rcoul + one_over_Rcoul2*(r-Rcoul));
                                }

                                if (r < Rlj)
                                {
                                    const LJPair &ljpair = ljpairs.constData()[
                                                            ljpairs.map(param0.ljid,
                                                                        param1.ljid)];

                                    double sig_over_dist6 = pow_6(ljpair.sigma()*one_over_r);
                                    double sig_over_dist12 = pow_2(sig_over_dist6);
            
                                    iljnrg += ljpair.epsilon() * (sig_over_dist12 -
                                                                  sig_over_dist6);
                                }
                            }
                        }
                    }
                }
                #endif
                
                cnrg += icnrg;
                ljnrg += iljnrg;
            }
        }
    }
    else if (use_atomistic_cutoff)
    {
        for (quint32 igroup=0; igroup<ngroups0; ++igroup)
        {
            const Parameters::Array &params0 = molparams0_array[igroup];
            const CoordGroup &group0 = groups0_array[igroup];
            const AABox &aabox0 = group0.aaBox();
            const quint32 nats0 = group0.count();
            const Parameter *params0_array = params0.constData();
            
            for (quint32 jgroup=0; jgroup<ngroups1; ++jgroup)
            {
                const CoordGroup &group1 = groups1_array[jgroup];
                const Parameters::Array &params1 = molparams1_array[jgroup];
                
                //check first that these two CoordGroups could be within cutoff
                //(if there is only one CutGroup in both molecules then this
                //test has already been performed and passed)
                const bool within_cutoff = (ngroups0 == 1 and ngroups1 == 1) or not
                                            spce->beyond(Rc, aabox0, group1.aaBox());
                
                if (not within_cutoff)
                    //this CutGroup is beyond the cutoff distance
                    continue;
                
                //calculate all of the interatomic distances
                const double mindist = spce->calcDist(group0, group1, distmat);
                
                if (mindist > Rc)
                {
                    //all of the atoms are definitely beyond cutoff
                    continue;
                }
                   
                double icnrg = 0;
                double iljnrg = 0;
                
                //loop over all interatomic pairs and calculate the energies
                const quint32 nats1 = group1.count();
                const Parameter *params1_array = params1.constData();

                #ifdef SIRE_USE_SSE
                {
                    const int remainder = nats1 % 2;
                    
                    __m128d sse_cnrg = { 0, 0 };
                    __m128d sse_ljnrg = { 0, 0 };

                    const __m128d sse_one = { 1.0, 1.0 };
                    const __m128d sse_Rcoul = _mm_set1_pd(Rcoul);
                    const __m128d sse_Rlj = _mm_set1_pd(Rlj);
                    
                    for (quint32 i=0; i<nats0; ++i)
                    {
                        distmat.setOuterIndex(i);
                        const Parameter &param0 = params0_array[i];
                        
                        __m128d sse_chg0 = _mm_set_pd( param0.reduced_charge, 
                                                       param0.reduced_charge );
                                             
                        //process atoms in pairs (so can then use SSE)
                        for (quint32 j=0; j<nats1-1; j += 2)
                        {
                            const Parameter &param10 = params1_array[j];
                            const Parameter &param11 = params1_array[j+1];
                            
                            const __m128d sse_r = _mm_set_pd( distmat[j], distmat[j+1] );
                            const __m128d sse_one_over_r = _mm_div_pd(sse_one, sse_r);

                            // coulomb calculation
                            {
                                const __m128d sse_in_cutoff = _mm_cmplt_pd(sse_r, sse_Rcoul);

                                __m128d nrg = _mm_set_pd( param10.reduced_charge,
                                                          param11.reduced_charge );
                        
                                nrg = _mm_mul_pd(nrg, sse_chg0);
                        
                                nrg = _mm_mul_pd(nrg, sse_one_over_r);
                                nrg = _mm_and_pd(nrg, sse_in_cutoff);
                                
                                sse_cnrg = _mm_add_pd(sse_cnrg, nrg);
                            }
                            
                            // lj calculation
                            {
                                const __m128d sse_in_cutoff = _mm_cmplt_pd(sse_r, sse_Rlj);

                                const LJPair &ljpair0 = ljpairs.constData()[
                                                         ljpairs.map(param0.ljid,
                                                                     param10.ljid)];
                        
                                const LJPair &ljpair1 = ljpairs.constData()[
                                                         ljpairs.map(param0.ljid,
                                                                     param11.ljid)];
                        
                                __m128d sse_sig = _mm_set_pd( ljpair0.sigma(), ljpair1.sigma() );
                                __m128d sse_eps = _mm_set_pd( ljpair0.epsilon(),
                                                              ljpair1.epsilon() );
                                                        
                                //calculate (sigma/r)^6 and (sigma/r)^12
                                __m128d sse_sig_over_dist2 = _mm_mul_pd(sse_sig, sse_one_over_r);
                                sse_sig_over_dist2 = _mm_mul_pd( sse_sig_over_dist2,
                                                                 sse_sig_over_dist2 );
                                                         
                                __m128d sse_sig_over_dist6 = _mm_mul_pd(sse_sig_over_dist2,
                                                                        sse_sig_over_dist2);
                                                            
                                sse_sig_over_dist6 = _mm_mul_pd(sse_sig_over_dist6,
                                                                sse_sig_over_dist2);
                                                         
                                __m128d sse_sig_over_dist12 = _mm_mul_pd(sse_sig_over_dist6,
                                                                         sse_sig_over_dist6);
                                                  
                                __m128d nrg = _mm_sub_pd(sse_sig_over_dist12,
                                                         sse_sig_over_dist6);
                                                     
                                nrg = _mm_mul_pd(nrg, sse_eps);
                                nrg = _mm_and_pd(nrg, sse_in_cutoff);
                                sse_ljnrg = _mm_add_pd(sse_ljnrg, nrg);
                            }
                        }
                              
                        if (remainder == 1)
                        {
                            const Parameter &param1 = params1_array[nats1-1];

                            const double r = distmat[nats1-1];
                            
                            if (r < Rc)
                            {
                                const double one_over_r = double(1) / r;
                            
                                if (r < Rcoul)
                                {
                                    icnrg += param0.reduced_charge * param1.reduced_charge *
                                                one_over_r;
                                }

                                if (r < Rlj)
                                {
                                    const LJPair &ljpair = ljpairs.constData()[
                                                            ljpairs.map(param0.ljid,
                                                                        param1.ljid)];

                                    double sig_over_dist6 = pow_6(ljpair.sigma()*one_over_r);
                                    double sig_over_dist12 = pow_2(sig_over_dist6);
            
                                    iljnrg += ljpair.epsilon() * (sig_over_dist12 -
                                                                  sig_over_dist6);
                                }
                            }
                        }
                    }
                    
                    icnrg += *((const double*)&sse_cnrg) +
                             *( ((const double*)&sse_cnrg) + 1 );
                             
                    iljnrg += *((const double*)&sse_ljnrg) +
                              *( ((const double*)&sse_ljnrg) + 1 );
                }
                #else
                {
                    for (quint32 i=0; i<nats0; ++i)
                    {
                        distmat.setOuterIndex(i);
                        const Parameter &param0 = params0_array[i];
                    
                        for (quint32 j=0; j<nats1; ++j)
                        {
                            const Parameter &param1 = params1_array[j];

                            const double r = distmat[j];
                            
                            if (r < Rc)
                            {
                                const double one_over_r = double(1) / r;
                            
                                if (r < Rcoul)
                                {
                                    icnrg += param0.reduced_charge * param1.reduced_charge *
                                                one_over_r;
                                }

                                if (r < Rlj)
                                {
                                    const LJPair &ljpair = ljpairs.constData()[
                                                            ljpairs.map(param0.ljid,
                                                                        param1.ljid)];

                                    double sig_over_dist6 = pow_6(ljpair.sigma()*one_over_r);
                                    double sig_over_dist12 = pow_2(sig_over_dist6);
            
                                    iljnrg += ljpair.epsilon() * (sig_over_dist12 -
                                                                  sig_over_dist6);
                                }
                            }
                        }
                    }
                }
                #endif
                
                cnrg += icnrg;
                ljnrg += iljnrg;
            }
        }
    }
    else if (use_reaction_field)
    {
        //use the reaction field potential
        // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
        // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
        // c = (1/r_c) * (3 eps)/(2 eps + 1)
        const double k_rf = (1.0 / pow_3(Rcoul)) * ( (rf_dielectric_constant-1) /
                                                     (2*rf_dielectric_constant + 1) );
        const double c_rf = (1.0 / Rcoul) * ( (3*rf_dielectric_constant) /
                                              (2*rf_dielectric_constant + 1) );
        
        for (quint32 igroup=0; igroup<ngroups0; ++igroup)
        {
            const Parameters::Array &params0 = molparams0_array[igroup];
            const CoordGroup &group0 = groups0_array[igroup];
            const AABox &aabox0 = group0.aaBox();
            const quint32 nats0 = group0.count();
            const Parameter *params0_array = params0.constData();
            
            for (quint32 jgroup=0; jgroup<ngroups1; ++jgroup)
            {
                const CoordGroup &group1 = groups1_array[jgroup];
                const Parameters::Array &params1 = molparams1_array[jgroup];
                
                //check first that these two CoordGroups could be within cutoff
                //(if there is only one CutGroup in both molecules then this
                //test has already been performed and passed)
                const bool within_cutoff = (ngroups0 == 1 and ngroups1 == 1) or not
                                            spce->beyond(Rc, aabox0, group1.aaBox());
                
                if (not within_cutoff)
                    //this CutGroup is beyond the cutoff distance
                    continue;
                
                //calculate all of the interatomic distances
                const double mindist = spce->calcDist(group0, group1, distmat);
                
                if (mindist > Rc)
                {
                    //all of the atoms are definitely beyond cutoff
                    continue;
                }
                   
                double icnrg = 0;
                double iljnrg = 0;
                
                //loop over all interatomic pairs and calculate the energies
                const quint32 nats1 = group1.count();
                const Parameter *params1_array = params1.constData();

                #ifdef SIRE_USE_SSE
                {
                    const int remainder = nats1 % 2;
                    
                    __m128d sse_cnrg = { 0, 0 };
                    __m128d sse_ljnrg = { 0, 0 };

                    const __m128d sse_one = { 1.0, 1.0 };
                    
                    const __m128d sse_Rcoul = _mm_set1_pd(Rcoul);
                    const __m128d sse_Rlj = _mm_set1_pd(Rlj);
                    
                    const __m128d sse_k_rf = _mm_set1_pd(k_rf);
                    const __m128d sse_c_rf = _mm_set1_pd(c_rf);
                    
                    for (quint32 i=0; i<nats0; ++i)
                    {
                        distmat.setOuterIndex(i);
                        const Parameter &param0 = params0_array[i];
                        
                        __m128d sse_chg0 = _mm_set_pd( param0.reduced_charge, 
                                                       param0.reduced_charge );
                                             
                        //process atoms in pairs (so can then use SSE)
                        for (quint32 j=0; j<nats1-1; j += 2)
                        {
                            const Parameter &param10 = params1_array[j];
                            const Parameter &param11 = params1_array[j+1];
                            
                            const __m128d sse_r = _mm_set_pd( distmat[j], distmat[j+1] );
                            const __m128d sse_one_over_r = _mm_div_pd(sse_one, sse_r);

                            //coulomb calculation
                            {
                                const __m128d sse_in_cutoff = _mm_cmplt_pd(sse_r, sse_Rcoul);
                                
                                __m128d nrg = _mm_mul_pd(sse_r, sse_r);
                                nrg = _mm_mul_pd(nrg, sse_k_rf);
                                nrg = _mm_sub_pd(nrg, sse_c_rf);
                                nrg = _mm_add_pd(nrg, sse_one_over_r);

                                __m128d sse_chg = _mm_set_pd( param10.reduced_charge,
                                                              param11.reduced_charge );
                        
                                sse_chg = _mm_mul_pd(sse_chg, sse_chg0);
                        
                                nrg = _mm_mul_pd(sse_chg, nrg);

                                nrg = _mm_and_pd(nrg, sse_in_cutoff);
                                
                                sse_cnrg = _mm_add_pd(sse_cnrg, nrg);
                            }
                            
                            //lj calculation
                            {
                                const __m128d sse_in_cutoff = _mm_cmplt_pd(sse_r, sse_Rlj);

                                const LJPair &ljpair0 = ljpairs.constData()[
                                                            ljpairs.map(param0.ljid,
                                                                        param10.ljid)];
                        
                                const LJPair &ljpair1 = ljpairs.constData()[
                                                            ljpairs.map(param0.ljid,
                                                                        param11.ljid)];
                        
                                __m128d sse_sig = _mm_set_pd( ljpair0.sigma(), ljpair1.sigma() );
                                __m128d sse_eps = _mm_set_pd( ljpair0.epsilon(),
                                                              ljpair1.epsilon() );
                                                        
                                //calculate (sigma/r)^6 and (sigma/r)^12
                                __m128d sse_sig_over_dist2 = _mm_mul_pd(sse_sig, sse_one_over_r);
                                sse_sig_over_dist2 = _mm_mul_pd( sse_sig_over_dist2,
                                                                 sse_sig_over_dist2 );
                                                         
                                __m128d sse_sig_over_dist6 = _mm_mul_pd(sse_sig_over_dist2,
                                                                        sse_sig_over_dist2);
                                                            
                                sse_sig_over_dist6 = _mm_mul_pd(sse_sig_over_dist6,
                                                                sse_sig_over_dist2);
                                                         
                                __m128d sse_sig_over_dist12 = _mm_mul_pd(sse_sig_over_dist6,
                                                                         sse_sig_over_dist6);
                                                  
                                __m128d nrg = _mm_sub_pd(sse_sig_over_dist12,
                                                         sse_sig_over_dist6);
                                                     
                                nrg = _mm_mul_pd(nrg, sse_eps);
                                nrg = _mm_and_pd(nrg, sse_in_cutoff);
                                sse_ljnrg = _mm_add_pd(sse_ljnrg, nrg);
                            }
                        }
                              
                        if (remainder == 1)
                        {
                            const Parameter &param1 = params1_array[nats1-1];

                            const double r = distmat[nats1-1];
                            
                            if (r < Rc)
                            {
                                const double one_over_r = double(1) / r;
                            
                                if (r < Rcoul)
                                {
                                    icnrg += param0.reduced_charge * param1.reduced_charge *
                                                (one_over_r + k_rf*r*r - c_rf);
                                }

                                if (r < Rlj)
                                {
                                    const LJPair &ljpair = ljpairs.constData()[
                                                            ljpairs.map(param0.ljid,
                                                                        param1.ljid)];

                                    double sig_over_dist6 = pow_6(ljpair.sigma()*one_over_r);
                                    double sig_over_dist12 = pow_2(sig_over_dist6);
            
                                    iljnrg += ljpair.epsilon() * (sig_over_dist12 -
                                                                  sig_over_dist6);
                                }
                            }
                        }
                    }
                    
                    icnrg += *((const double*)&sse_cnrg) +
                             *( ((const double*)&sse_cnrg) + 1 );
                             
                    iljnrg += *((const double*)&sse_ljnrg) +
                              *( ((const double*)&sse_ljnrg) + 1 );
                }
                #else
                {
                    for (quint32 i=0; i<nats0; ++i)
                    {
                        distmat.setOuterIndex(i);
                        const Parameter &param0 = params0_array[i];
                    
                        for (quint32 j=0; j<nats1; ++j)
                        {
                            const Parameter &param1 = params1_array[j];

                            const double r = distmat[j];
                            
                            if (r < Rc)
                            {
                                const double one_over_r = double(1) / r;
                            
                                if (r < Rcoul)
                                {
                                    icnrg += param0.reduced_charge * param1.reduced_charge *
                                                (one_over_r + k_rf*r*r - c_rf);
                                }

                                if (r < Rlj)
                                {
                                    const LJPair &ljpair = ljpairs.constData()[
                                                            ljpairs.map(param0.ljid,
                                                                        param1.ljid)];

                                    double sig_over_dist6 = pow_6(ljpair.sigma()*one_over_r);
                                    double sig_over_dist12 = pow_2(sig_over_dist6);
            
                                    iljnrg += ljpair.epsilon() * (sig_over_dist12 -
                                                                  sig_over_dist6);
                                }
                            }
                        }
                    }
                }
                #endif
                
                cnrg += icnrg;
                ljnrg += iljnrg;
            }
        }
    }
    else
    {
        //using the group-based cutoff with feather function
        for (quint32 igroup=0; igroup<ngroups0; ++igroup)
        {
            const Parameters::Array &params0 = molparams0_array[igroup];

            const CoordGroup &group0 = groups0_array[igroup];
            const AABox &aabox0 = group0.aaBox();
            const quint32 nats0 = group0.count();
            const Parameter *params0_array = params0.constData();
        
            for (quint32 jgroup=0; jgroup<ngroups1; ++jgroup)
            {
                const CoordGroup &group1 = groups1_array[jgroup];
                const Parameters::Array &params1 = molparams1_array[jgroup];

                //check first that these two CoordGroups could be within cutoff
                //(if there is only one CutGroup in both molecules then this
                //test has already been performed and passed)
                const bool within_cutoff = (ngroups0 == 1 and ngroups1 == 1) or not
                                            spce->beyond(Rc,aabox0, group1.aaBox());
                
                if (not within_cutoff)
                    //this CutGroup is either the cutoff distance
                    continue;
                
                //calculate all of the interatomic distances
                const double mindist = spce->calcDist(group0, group1, distmat);
                
                if (mindist > Rc)
                {
                    //all of the atoms are definitely beyond cutoff
                    continue;
                }
                   
                double icnrg = 0;
                double iljnrg = 0;
                
                //loop over all interatomic pairs and calculate the energies
                const quint32 nats1 = group1.count();
                const Parameter *params1_array = params1.constData();

                #ifdef SIRE_USE_SSE
                {
                    const int remainder = nats1 % 2;
                    
                    __m128d sse_cnrg = { 0, 0 };
                    __m128d sse_ljnrg = { 0, 0 };

                    const __m128d sse_one = { 1.0, 1.0 };
                    
                    for (quint32 i=0; i<nats0; ++i)
                    {
                        distmat.setOuterIndex(i);
                        const Parameter &param0 = params0_array[i];
                        
                        __m128d sse_chg0 = _mm_set_pd( param0.reduced_charge, 
                                                       param0.reduced_charge );
                                             
                        //process atoms in pairs (so can then use SSE)
                        for (quint32 j=0; j<nats1-1; j += 2)
                        {
                            const Parameter &param10 = params1_array[j];
                            const Parameter &param11 = params1_array[j+1];
                            
                            __m128d sse_dist = _mm_set_pd( distmat[j], distmat[j+1] );
                            __m128d sse_chg1 = _mm_set_pd( param10.reduced_charge,
                                                           param11.reduced_charge );
                                               
                            const LJPair &ljpair0 = ljpairs.constData()[
                                                    ljpairs.map(param0.ljid,
                                                                param10.ljid)];
                        
                            const LJPair &ljpair1 = ljpairs.constData()[
                                                    ljpairs.map(param0.ljid,
                                                                param11.ljid)];
                        
                            __m128d sse_sig = _mm_set_pd( ljpair0.sigma(), ljpair1.sigma() );
                            __m128d sse_eps = _mm_set_pd( ljpair0.epsilon(), 
                                                          ljpair1.epsilon() );
                            
                            sse_dist = _mm_div_pd(sse_one, sse_dist);
                            
                            //calculate the coulomb energy
                            __m128d tmp = _mm_mul_pd(sse_chg0, sse_chg1);
                            tmp = _mm_mul_pd(tmp, sse_dist);
                            sse_cnrg = _mm_add_pd(sse_cnrg, tmp);
                            
                            //calculate (sigma/r)^6 and (sigma/r)^12
                            __m128d sse_sig_over_dist2 = _mm_mul_pd(sse_sig, sse_dist);
                            sse_sig_over_dist2 = _mm_mul_pd( sse_sig_over_dist2,  
                                                             sse_sig_over_dist2 );
                                                         
                            __m128d sse_sig_over_dist6 = _mm_mul_pd(sse_sig_over_dist2,
                                                                    sse_sig_over_dist2);
                                                            
                            sse_sig_over_dist6 = _mm_mul_pd(sse_sig_over_dist6,
                                                            sse_sig_over_dist2);
                                                         
                            __m128d sse_sig_over_dist12 = _mm_mul_pd(sse_sig_over_dist6,
                                                                     sse_sig_over_dist6);
                                                  
                            //calculate LJ energy (the factor of 4 is added later)
                            tmp = _mm_sub_pd(sse_sig_over_dist12, 
                                             sse_sig_over_dist6);
                                                     
                            tmp = _mm_mul_pd(tmp, sse_eps);
                            sse_ljnrg = _mm_add_pd(sse_ljnrg, tmp);
                        }
                              
                        if (remainder == 1)
                        {
                            const Parameter &param1 = params1_array[nats1-1];

                            const double invdist = double(1) / distmat[nats1-1];
                            
                            icnrg += param0.reduced_charge * param1.reduced_charge 
                                        * invdist;

                            const LJPair &ljpair = ljpairs.constData()[
                                                    ljpairs.map(param0.ljid,
                                                                param1.ljid)];
                            
                            double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                            double sig_over_dist12 = pow_2(sig_over_dist6);
        
                            iljnrg += ljpair.epsilon() * (sig_over_dist12 - 
                                                          sig_over_dist6);
                        }
                    }
                    
                    icnrg += *((const double*)&sse_cnrg) +
                             *( ((const double*)&sse_cnrg) + 1 );
                             
                    iljnrg += *((const double*)&sse_ljnrg) +
                              *( ((const double*)&sse_ljnrg) + 1 );
                }
                #else
                {
                    for (quint32 i=0; i<nats0; ++i)
                    {
                        distmat.setOuterIndex(i);
                        const Parameter &param0 = params0_array[i];
                    
                        for (quint32 j=0; j<nats1; ++j)
                        {
                            const Parameter &param1 = params1_array[j];

                            const double invdist = double(1) / distmat[j];
                            
                            icnrg += param0.reduced_charge * param1.reduced_charge 
                                        * invdist;

                            const LJPair &ljpair = ljpairs.constData()[
                                                    ljpairs.map(param0.ljid,
                                                                param1.ljid)];

                            double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                            double sig_over_dist12 = pow_2(sig_over_dist6);
        
                            iljnrg += ljpair.epsilon() * (sig_over_dist12 - 
                                                          sig_over_dist6);
                        }
                    }
                }
                #endif
                
                //now add these energies onto the total for the molecule,
                //scaled by any non-bonded feather factor
                if (mindist > switchfunc->electrostaticFeatherDistance())
                {
                    cnrg += switchfunc->electrostaticScaleFactor( Length(mindist) ) * icnrg;
                }
                else
                {
                    cnrg += icnrg;
                }
                
                if (mindist > switchfunc->vdwFeatherDistance())
                {
                    ljnrg += switchfunc->vdwScaleFactor( Length(mindist) ) * iljnrg;
                }
                else
                {
                    ljnrg += iljnrg;
                }
            }
        }
    }
    
    //add this molecule pair's energy onto the total
    //(also multiply LJ by 4 as it is 4 * epsilon ((sig/r)^12 - (sig/r)^6))
    Energy cljnrg(scale_energy * cnrg, 4 * scale_energy * ljnrg);

    //energy += Energy(scale_energy * cnrg, 4 * scale_energy * ljnrg);
    energy = energy + cljnrg;
}

/** Add to the energies in 'energies0' the energy on 'mol0' caused
    by 'mol1' */
void InterCLJPotential::_pvt_calculateEnergy(const InterCLJPotential::Molecule &mol0, 
                                            const InterCLJPotential::Molecule &mol1,
                                            MolEnergyTable &energies0, 
                                            InterCLJPotential::EnergyWorkspace &distmat,
                                            double scale_energy) const
{
    BOOST_ASSERT( mol0.molecule().data().info().nCutGroups() == energies0.nCutGroups() );
    BOOST_ASSERT( mol0.molecule().data().number() == energies0.molNum() );

    const quint32 ngroups0 = mol0.nCutGroups();
    const quint32 ngroups1 = mol1.nCutGroups();
    
    const CoordGroup *groups0_array = mol0.coordinates().constData();
    const CoordGroup *groups1_array = mol1.coordinates().constData();
    
    BOOST_ASSERT(mol0.parameters().atomicParameters().count() == int(ngroups0));
    BOOST_ASSERT(mol1.parameters().atomicParameters().count() == int(ngroups1));
    
    const Parameters::Array *molparams0_array 
                                    = mol0.parameters().atomicParameters().constData();
    const Parameters::Array *molparams1_array
                                    = mol1.parameters().atomicParameters().constData();
    
    const MolEnergyTable::Array *energies0_array = energies0.constData();
    
    if (use_reaction_field)
    {
      //loop over all pairs of CutGroups in the two molecules
      for (quint32 igroup=0; igroup<ngroups0; ++igroup)
      {
	//get the CGIdx of this group
	CGIdx cgidx_igroup = mol0.cgIdx(igroup);
	
        //get the index of this CutGroup in the forces array
        int energy0_idx = energies0.map(cgidx_igroup);
        
        if (energy0_idx == -1)
            //there is no space for the energies on this CutGroup in 
            //the energytable - were are therefore not interested in
            //this CutGroup
            continue;

        const Parameters::Array &params0 = molparams0_array[igroup];

        const CoordGroup &group0 = groups0_array[igroup];
        const AABox &aabox0 = group0.aaBox();
        const quint32 nats0 = group0.count();
        const Parameter *params0_array = params0.constData();
    
        //get the table that holds the energies of all the
        //atoms of this CutGroup (tables are indexed by CGIdx)
        BOOST_ASSERT(energies0_array[energy0_idx].count() == int(nats0));
    
        Vector *group_energies0_array = energies0.data(energy0_idx);

        //ok, we are interested in the energies of this CutGroup
        // - calculate all of the energies of this group due to interactions
        //   with all of the CutGroups in mol1 
        for (quint32 jgroup=0; jgroup<ngroups1; ++jgroup)
        {
            const CoordGroup &group1 = groups1_array[jgroup];
            const Parameters::Array &params1 = molparams1_array[jgroup];

            //check first that these two CoordGroups could be within cutoff
            //(if there is only one CutGroup in both molecules then this
            //test has already been performed and passed)

	    // JM dec 12. It doesn't look like the within cutoff test has been performed and 
	    // passed. 
	    // Expression below might speed things up a bit. 
            //const bool within_cutoff = not spce->beyond(switchfunc->cutoffDistance(), 
	    //					aabox0, group1.aaBox());
            const bool within_cutoff = (ngroups0 == 1 and ngroups1 == 1) or not
                                        spce->beyond(switchfunc->cutoffDistance(), 
                                                     aabox0, group1.aaBox());
            
            if (not within_cutoff)
                //this CutGroup is beyond the cutoff distance
                continue;
            
            //calculate all of the interatomic distances
            const double mindist = spce->calcDist(group0, group1, distmat);
            
            if (mindist > switchfunc->cutoffDistance())
                //all of the atoms are definitely beyond cutoff
                continue;

	    // JM dec 12. We are going to need the reaction field parameters. 
	    // This could be initialised somewhere else...
	    double Rc = switchfunc->electrostaticCutoffDistance();
        
	    if (Rc != switchfunc->vdwCutoffDistance())
	      throw SireError::unsupported( QObject::tr(
							"This code does not support having a reaction field together "
							"with different coulomb and vdw cutoffs..."), CODELOC );
	    if (Rc > 1e9)
	      {
		Rc = 1e9;
	      }

	    const double k_rf = (1.0 / pow_3(Rc)) * ( (rf_dielectric_constant-1) /
						(2*rf_dielectric_constant + 1) );
	    const double c_rf = (1.0 / Rc) * ( (3*rf_dielectric_constant) /
					       (2*rf_dielectric_constant + 1) );
	    
            const quint32 nats1 = group1.count();
            
            //loop over all interatomic pairs and calculate the energies
            const Parameter *params1_array = params1.constData();
	    
	    for (quint32 i=0; i<nats0; ++i)
	    {
		distmat.setOuterIndex(i);
		const Parameter &param0 = params0_array[i];
                
		double icnrg = 0;
		double iljnrg = 0;

		for (quint32 j=0; j<nats1; ++j)
		{
		  //do both coulomb and LJ

		  const Parameter &param1 = params1_array[j];
		      
		  const double r = distmat[j];
                            
		  if (r < Rc)
		    {
		      const double one_over_r = double(1) / r;
                      
		      icnrg += param0.reduced_charge * param1.reduced_charge *
			(one_over_r + k_rf*r*r - c_rf);
		      
		      const LJPair &ljpair = ljpairs.constData()[
					      ljpairs.map(param0.ljid,
						     param1.ljid)];
		      
		      double sig_over_dist6 = pow_6(ljpair.sigma()*one_over_r);
		      double sig_over_dist12 = pow_2(sig_over_dist6);
			      
		      iljnrg += ljpair.epsilon() * (sig_over_dist12 - sig_over_dist6);
		    }
		}

		iljnrg = scale_energy * 0.5 * 4 * iljnrg ;
		icnrg = scale_energy * 0.5 * icnrg;

		Vector nrg = Vector( icnrg + iljnrg, icnrg, iljnrg);

		group_energies0_array[i] += nrg;

	    }

        } // end of loop over jgroup CutGroups
	
      } // end of loop over igroup CutGroups
    }
    else
    {
      //loop over all pairs of CutGroups in the two molecules
      for (quint32 igroup=0; igroup<ngroups0; ++igroup)
      {
	//get the CGIdx of this group
	CGIdx cgidx_igroup = mol0.cgIdx(igroup);
	
        //get the index of this CutGroup in the forces array
        int energy0_idx = energies0.map(cgidx_igroup);
        
        if (energy0_idx == -1)
            //there is no space for the energies on this CutGroup in 
            //the energytable - were are therefore not interested in
            //this CutGroup
            continue;

        const Parameters::Array &params0 = molparams0_array[igroup];

        const CoordGroup &group0 = groups0_array[igroup];
        const AABox &aabox0 = group0.aaBox();
        const quint32 nats0 = group0.count();
        const Parameter *params0_array = params0.constData();
    
        //get the table that holds the energies of all the
        //atoms of this CutGroup (tables are indexed by CGIdx)
        BOOST_ASSERT(energies0_array[energy0_idx].count() == int(nats0));
    
        Vector *group_energies0_array = energies0.data(energy0_idx);

        //ok, we are interested in the energies of this CutGroup
        // - calculate all of the energies of this group due to interactions
        //   with all of the CutGroups in mol1 
        for (quint32 jgroup=0; jgroup<ngroups1; ++jgroup)
        {
            const CoordGroup &group1 = groups1_array[jgroup];
            const Parameters::Array &params1 = molparams1_array[jgroup];

            //check first that these two CoordGroups could be within cutoff
            //(if there is only one CutGroup in both molecules then this
            //test has already been performed and passed)
            const bool within_cutoff = (ngroups0 == 1 and ngroups1 == 1) or not
                                        spce->beyond(switchfunc->cutoffDistance(), 
                                                     aabox0, group1.aaBox());
            
            if (not within_cutoff)
                //this CutGroup is beyond the cutoff distance
                continue;
            
            //calculate all of the interatomic distances
            const double mindist = spce->calcDist(group0, group1, distmat);
            
            if (mindist > switchfunc->cutoffDistance())
                //all of the atoms are definitely beyond cutoff
                continue;

            const quint32 nats1 = group1.count();
            
            //loop over all interatomic pairs and calculate the energies
            const Parameter *params1_array = params1.constData();
	    
	    for (quint32 i=0; i<nats0; ++i)
	    {
		distmat.setOuterIndex(i);
		const Parameter &param0 = params0_array[i];
                
		double icnrg = 0;
		double iljnrg = 0;

		for (quint32 j=0; j<nats1; ++j)
		{
		    const Parameter &param1 = params1_array[j];
		    
		    const double invdist = double(1) / distmat[j];
			
		    icnrg += param0.reduced_charge * param1.reduced_charge 
                                    * invdist;
                    
		    const LJPair &ljpair = ljpairs.constData()[
						   ljpairs.map(param0.ljid,
							       param1.ljid)];

		    double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
		    double sig_over_dist12 = pow_2(sig_over_dist6);
    
		    iljnrg += ljpair.epsilon() * (sig_over_dist12 - 
                                                      sig_over_dist6);
		}

		//are we shifting the electrostatic potential?
		if (use_electrostatic_shifting)
		  {
		    icnrg -= this->totalCharge(params0) * this->totalCharge(params1)
		      / switchfunc->electrostaticCutoffDistance();
		    
		  }

		//now add these energies onto the total for the molecule,
		//scaled by any non-bonded feather factor
		if (mindist > switchfunc->featherDistance())
		  {
		    icnrg = switchfunc->electrostaticScaleFactor( Length(mindist) ) * icnrg;
		    iljnrg = switchfunc->vdwScaleFactor( Length(mindist) ) * iljnrg;
		  }
		
		iljnrg = scale_energy * 0.5 * 4 * iljnrg ;
		icnrg = scale_energy * 0.5 * icnrg;

		Vector nrg = Vector( icnrg + iljnrg, icnrg, iljnrg);

		group_energies0_array[i] += nrg;
	    }

        } // end of loop over jgroup CutGroups
	
      } // end of loop over igroup CutGroups
    }
}

/** Add to the forces in 'forces0' the forces acting on 'mol0' caused
    by 'mol1' */
void InterCLJPotential::_pvt_calculateForce(const InterCLJPotential::Molecule &mol0, 
                                            const InterCLJPotential::Molecule &mol1,
                                            MolForceTable &forces0, 
                                            InterCLJPotential::ForceWorkspace &distmat,
                                            double scale_force) const
{
    BOOST_ASSERT( mol0.molecule().data().info().nCutGroups() == forces0.nCutGroups() );
    BOOST_ASSERT( mol0.molecule().data().number() == forces0.molNum() );

    const quint32 ngroups0 = mol0.nCutGroups();
    const quint32 ngroups1 = mol1.nCutGroups();
    
    const CoordGroup *groups0_array = mol0.coordinates().constData();
    const CoordGroup *groups1_array = mol1.coordinates().constData();
    
    BOOST_ASSERT(mol0.parameters().atomicParameters().count() == int(ngroups0));
    BOOST_ASSERT(mol1.parameters().atomicParameters().count() == int(ngroups1));
    
    const Parameters::Array *molparams0_array 
                                    = mol0.parameters().atomicParameters().constData();
    const Parameters::Array *molparams1_array
                                    = mol1.parameters().atomicParameters().constData();
    
    const MolForceTable::Array *forces0_array = forces0.constData();
    

    if (use_reaction_field)
    {
       //qDebug() << " Using a reaction field ";

      //loop over all pairs of CutGroups in the two molecules
      for (quint32 igroup=0; igroup<ngroups0; ++igroup)
      {
	//get the CGIdx of this group
        CGIdx cgidx_igroup = mol0.cgIdx(igroup);

        //get the index of this CutGroup in the forces array
        int force0_idx = forces0.map(cgidx_igroup);
        
        if (force0_idx == -1)
            //there is no space for the forces on this CutGroup in 
            //the forcetable - were are therefore not interested in
            //this CutGroup
            continue;

        const Parameters::Array &params0 = molparams0_array[igroup];

        const CoordGroup &group0 = groups0_array[igroup];
        const AABox &aabox0 = group0.aaBox();
        const quint32 nats0 = group0.count();
        const Parameter *params0_array = params0.constData();
    
        //get the table that holds the forces acting on all of the
        //atoms of this CutGroup (tables are indexed by CGIdx)
        BOOST_ASSERT(forces0_array[force0_idx].count() == int(nats0));
    
        Vector *group_forces0_array = forces0.data(force0_idx);

        //ok, we are interested in the forces acting on this CutGroup
        // - calculate all of the forces on this group interacting
        //   with all of the CutGroups in mol1 
        for (quint32 jgroup=0; jgroup<ngroups1; ++jgroup)
        {
            const CoordGroup &group1 = groups1_array[jgroup];
            const Parameters::Array &params1 = molparams1_array[jgroup];

            //check first that these two CoordGroups could be within cutoff
            //(if there is only one CutGroup in both molecules then this
            //test has already been performed and passed)
            const bool within_cutoff = (ngroups0 == 1 and ngroups1 == 1) or not
                                        spce->beyond(switchfunc->cutoffDistance(), 
                                                     aabox0, group1.aaBox());
            
            if (not within_cutoff)
                //this CutGroup is beyond the cutoff distance
                continue;
            
            //calculate all of the interatomic distances
            const double mindist = spce->calcDistVectors(group0, group1, distmat);
            
            if (mindist > switchfunc->cutoffDistance())
                //all of the atoms are definitely beyond cutoff
                continue;

	    // JM Dec 12 Define reaction field parameters
	    //use the reaction field potential
	    // dE/dr = (q1 q2 / 4 pi eps_0) * ( -1/r^2 + 2 k r )
	    // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)

	    double Rc = switchfunc->electrostaticCutoffDistance();
        
	    if (Rc != switchfunc->vdwCutoffDistance())
	      throw SireError::unsupported( QObject::tr(
			     "This code does not support having electrostatic shifting together "
				"with different coulomb and vdw cutoffs..."), CODELOC );
        
	    if (Rc > 1e9)
	      {
		Rc = 1e9;
	      }

	    const double k_rf = (1.0 / pow_3(Rc)) * ( (rf_dielectric_constant-1) /
						      (2*rf_dielectric_constant + 1) );
	    
            const quint32 nats1 = group1.count();
            
            //loop over all interatomic pairs and calculate the energies
            const Parameter *params1_array = params1.constData();
                
	    //no feather region with reaction field, can calculate directly
	    for (quint32 i=0; i<nats0; ++i)
	      {
		distmat.setOuterIndex(i);
		const Parameter &param0 = params0_array[i];

		Vector total_force;
                
		//qDebug() << " atom " << i << " total force initialised ";

		if (param0.ljid == 0)
		  {
		    //null LJ parameter - only add on the coulomb energy
		    for (quint32 j=0; j<nats1; ++j)
		      {
			const double q2 = param0.reduced_charge * 
			  params1_array[j].reduced_charge;
                        
			//calculate the coulomb force
			if (q2 != 0)
			  {
			    //Vector cforce = -(q2 / distmat[j].length2()) *
			    //  distmat[j].direction();
			    double r2 = distmat[j].length2();
			    double r = distmat[j].length();

			    Vector cforce = ( q2 *  ( (-1/r2) + 2*k_rf*r) ) *
			      distmat[j].direction();

			    total_force += cforce;
			    //qDebug() << " i no LJ...after  " << j << " force " << total_force.toString();
			  }
		      }
		  }
		else
		  {
		    for (quint32 j=0; j<nats1; ++j)
		      {
			//do both coulomb and LJ
			const Parameter &param1 = params1_array[j];
                        
			const double dist = distmat[j].length();
			const double invdist = double(1) / dist;
			const double invdist2 = pow_2(invdist);
			
			//calculate the force
			//Vector force = -(param0.reduced_charge * 
			//		 param1.reduced_charge * invdist2) 
			//  * distmat[j].direction();
			
			Vector force = ( (param0.reduced_charge *param1.reduced_charge) * ( (-invdist2) + 2*k_rf*dist) ) 
			  *  distmat[j].direction();

			//qDebug() <<  " ...coul only" << j << " force " << force.toString();
			
			if (param1.ljid != 0)
			  {
			    const LJPair &ljpair = ljpairs.constData()[
								       ljpairs.map(param0.ljid,
										   param1.ljid)];
			    
			    double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
			    double sig_over_dist12 = pow_2(sig_over_dist6);
			    
			    // dU/dr requires an extra power of r
			    sig_over_dist6 *= invdist;
			    sig_over_dist12 *= invdist;
			    
			    force += (4 * ljpair.epsilon() * (6.0*sig_over_dist6 - 
                                                              12.0*sig_over_dist12))
			      * distmat[j].direction();
			    
			    //qDebug() << "ljpair.sigma " << ljpair.sigma();
			    //qDebug() << "ljpair.epsilon" << ljpair.epsilon();
			    //qDebug() << "distmat[j].direction() " << distmat[j].direction().toString();
			    //qDebug() << " invdist " << invdist;
			    //qDebug() << " param0.ljid " << param0.ljid;
			    //qDebug() << " param1.ljid " << param1.ljid;
			    //qDebug() <<  " ...after LJ " << j << " force " << force.toString();
			    //force = Vector(0);
			    
			  }
                        
			total_force += force;
			//qDebug() << " ...after  " << j << " force " << total_force.toString();
		      }
		  }
		
		group_forces0_array[i] += scale_force * total_force;
		//qDebug() << " atom " << i << " force " << total_force.toString();
		//exit(0);
	      } // end of loop over i atoms
	}
      }
    }
    else
    {
      //loop over all pairs of CutGroups in the two molecules
      for (quint32 igroup=0; igroup<ngroups0; ++igroup)
	{
	  //get the CGIdx of this group
	  CGIdx cgidx_igroup = mol0.cgIdx(igroup);

	  //get the index of this CutGroup in the forces array
	  int force0_idx = forces0.map(cgidx_igroup);
        
	  if (force0_idx == -1)
            //there is no space for the forces on this CutGroup in 
            //the forcetable - were are therefore not interested in
            //this CutGroup
            continue;

	  const Parameters::Array &params0 = molparams0_array[igroup];

	  const CoordGroup &group0 = groups0_array[igroup];
	  const AABox &aabox0 = group0.aaBox();
	  const quint32 nats0 = group0.count();
	  const Parameter *params0_array = params0.constData();
    
	  //get the table that holds the forces acting on all of the
	  //atoms of this CutGroup (tables are indexed by CGIdx)
	  BOOST_ASSERT(forces0_array[force0_idx].count() == int(nats0));
    
	  Vector *group_forces0_array = forces0.data(force0_idx);

	  //ok, we are interested in the forces acting on this CutGroup
	  // - calculate all of the forces on this group interacting
	  //   with all of the CutGroups in mol1 
	  for (quint32 jgroup=0; jgroup<ngroups1; ++jgroup)
	    {
	      const CoordGroup &group1 = groups1_array[jgroup];
	      const Parameters::Array &params1 = molparams1_array[jgroup];

	      //check first that these two CoordGroups could be within cutoff
	      //(if there is only one CutGroup in both molecules then this
	      //test has already been performed and passed)
	      const bool within_cutoff = (ngroups0 == 1 and ngroups1 == 1) or not
                                        spce->beyond(switchfunc->cutoffDistance(), 
                                                     aabox0, group1.aaBox());
            
	      if (not within_cutoff)
                //this CutGroup is beyond the cutoff distance
                continue;
            
	      //calculate all of the interatomic distances
	      const double mindist = spce->calcDistVectors(group0, group1, distmat);
            
	      if (mindist > switchfunc->cutoffDistance())
                //all of the atoms are definitely beyond cutoff
                continue;

	      const quint32 nats1 = group1.count();
            
	      //loop over all interatomic pairs and calculate the energies
	      const Parameter *params1_array = params1.constData();

	      if (mindist > switchfunc->featherDistance())
		{
		  //we need to calculate the forces taking into account
		  //the derivative of the switching function!
                
		  //calculate the switching scale factors and their 
		  //derivatives
		  const double scl_coul = switchfunc->electrostaticScaleFactor( 
                                                                        Length(mindist) );
		  const double scl_lj = switchfunc->vdwScaleFactor( Length(mindist) );
                
		  Vector group_sep = (group1.aaBox().center() - aabox0.center())
                                        .normalise();
                
		  Vector dscl_coul = switchfunc->dElectrostaticScaleFactor( 
                                                                    Length(mindist) ) 
                                     * group_sep;
                                     
		  Vector dscl_lj = switchfunc->dVDWScaleFactor( Length(mindist) )
                                     * group_sep;
                
		  double shift_coul = 0;

		  if (use_electrostatic_shifting)
                    shift_coul = this->totalCharge(params0) * this->totalCharge(params1)
                                    / switchfunc->electrostaticCutoffDistance();
                
		  for (quint32 i=0; i<nats0; ++i)
		    {
		      distmat.setOuterIndex(i);
		      const Parameter &param0 = params0_array[i];
                
		      Vector total_force;
                
		      if (param0.ljid == 0)
			{
			  //null LJ parameter - only add on the coulomb energy
			  for (quint32 j=0; j<nats1; ++j)
			    {
			      const double q2 = param0.reduced_charge *
                                              params1_array[j].reduced_charge;
                            
			      if (q2 != 0)
				{
				  //calculate the coulomb energy
				  const double cnrg = q2 / distmat[j].length();
                                               
				  //calculate the coulomb force
				  Vector cforce = (scl_coul * -cnrg / distmat[j].length() *
                                                 distmat[j].direction()) +
                                             
                                                ((cnrg-shift_coul) * dscl_coul);
                        
				  total_force += cforce;
				}
			    }
			}
		      else
			{
			  for (quint32 j=0; j<nats1; ++j)
			    {
			      //do both coulomb and LJ
			      const Parameter &param1 = params1_array[j];
                        
			      const double invdist = double(1) / distmat[j].length();
                        
			      Vector force;
                            
			      const double q2 = param0.reduced_charge *
				param1.reduced_charge;
                        
			      if (q2 != 0)
				{
				  //calculate the energy
				  const double cnrg = q2 * invdist;
                        
				  //calculate the force
				  force = (scl_coul * -cnrg / distmat[j].length() *
                                         distmat[j].direction()) +
                                             
                                         ((cnrg-shift_coul) * dscl_coul);
				}
                              
			      if (param1.ljid != 0)
				{
				  const LJPair &ljpair = ljpairs.constData()[
                                                        ljpairs.map(param0.ljid,
                                                                    param1.ljid)];
                        
				  double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
				  double sig_over_dist12 = pow_2(sig_over_dist6);

				  //calculate the energy
				  const double ljnrg = 4 * ljpair.epsilon() *
                                                      (sig_over_dist12 - sig_over_dist6);

				  // dU/dr requires an extra power of r
				  sig_over_dist6 *= invdist;
				  sig_over_dist12 *= invdist;

				  force += ((4 * scl_lj * ljpair.epsilon() * 
                                            (6.0*sig_over_dist6 - 12.0*sig_over_dist12))
                                            * distmat[j].direction())
                                            
				    + (ljnrg * dscl_lj);
				}

			      total_force += force;
			    }
			}
                    
		      //update the forces array
		      group_forces0_array[i] += scale_force * total_force;
		    }
		}
	      else
		{
		  //not in the feather region, so can calculate the forces
		  //directly (also, no need to calculate shift, as 
		  //the shifting function is constant, so does not
		  //affect the gradient)
		  for (quint32 i=0; i<nats0; ++i)
		    {
		      distmat.setOuterIndex(i);
		      const Parameter &param0 = params0_array[i];

		      Vector total_force;
                
		      if (param0.ljid == 0)
			{
			  //null LJ parameter - only add on the coulomb energy
			  for (quint32 j=0; j<nats1; ++j)
			    {
			      const double q2 = param0.reduced_charge * 
                                              params1_array[j].reduced_charge;
                        
			      //calculate the coulomb force
			      if (q2 != 0)
				{
				  Vector cforce = -(q2 / distmat[j].length2()) *
                                                    distmat[j].direction();
                        
				  total_force += cforce;
				}
			    }
			}
		      else
			{
			  for (quint32 j=0; j<nats1; ++j)
			    {
			      //do both coulomb and LJ
			      const Parameter &param1 = params1_array[j];
                        
			      const double invdist = double(1) / distmat[j].length();
			      const double invdist2 = pow_2(invdist);
                        
			      //calculate the force
			      Vector force = -(param0.reduced_charge * 
					       param1.reduced_charge * invdist2) 
                                            
                                            * distmat[j].direction();
                              
			      if (param1.ljid != 0)
				{
				  const LJPair &ljpair = ljpairs.constData()[
                                                        ljpairs.map(param0.ljid,
                                                                    param1.ljid)];
                        
				  double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
				  double sig_over_dist12 = pow_2(sig_over_dist6);

				  // dU/dr requires an extra power of r
				  sig_over_dist6 *= invdist;
				  sig_over_dist12 *= invdist;

				  force += (4 * ljpair.epsilon() * (6.0*sig_over_dist6 - 
                                                              12.0*sig_over_dist12))
                                        * distmat[j].direction();
				}
			      
			      total_force += force;
			    }
			}
		      
		      group_forces0_array[i] += scale_force * total_force;

		    } // end of loop over i atoms

		} // end of if within feather

	    } // end of loop over jgroup CutGroups

	} // end of loop over igroup CutGroups
    }//end of else
}

/** Add to the forces in 'forces0' the forces acting on 'mol0' caused
    by 'mol1' */
void InterCLJPotential::_pvt_calculateCoulombForce(
                                            const InterCLJPotential::Molecule &mol0, 
                                            const InterCLJPotential::Molecule &mol1,
                                            MolForceTable &forces0, 
                                            InterCLJPotential::ForceWorkspace &distmat,
                                            double scale_force) const
{
    BOOST_ASSERT( mol0.molecule().data().info().nCutGroups() == forces0.nCutGroups() );
    BOOST_ASSERT( mol0.molecule().data().number() == forces0.molNum() );

    const quint32 ngroups0 = mol0.nCutGroups();
    const quint32 ngroups1 = mol1.nCutGroups();
    
    const CoordGroup *groups0_array = mol0.coordinates().constData();
    const CoordGroup *groups1_array = mol1.coordinates().constData();
    
    BOOST_ASSERT(mol0.parameters().atomicParameters().count() == int(ngroups0));
    BOOST_ASSERT(mol1.parameters().atomicParameters().count() == int(ngroups1));
    
    const Parameters::Array *molparams0_array 
                                    = mol0.parameters().atomicParameters().constData();
    const Parameters::Array *molparams1_array
                                    = mol1.parameters().atomicParameters().constData();
    
    const MolForceTable::Array *forces0_array = forces0.constData();
    
    //loop over all pairs of CutGroups in the two molecules
    for (quint32 igroup=0; igroup<ngroups0; ++igroup)
    {
        //get the CGIdx of this group
        CGIdx cgidx_igroup = mol0.cgIdx(igroup);

        //get the index of this CutGroup in the forces array
        int force0_idx = forces0.map(cgidx_igroup);
        
        if (force0_idx == -1)
            //there is no space for the forces on this CutGroup in 
            //the forcetable - were are therefore not interested in
            //this CutGroup
            continue;

        const Parameters::Array &params0 = molparams0_array[igroup];

        const CoordGroup &group0 = groups0_array[igroup];
        const AABox &aabox0 = group0.aaBox();
        const quint32 nats0 = group0.count();
        const Parameter *params0_array = params0.constData();
    
        //get the table that holds the forces acting on all of the
        //atoms of this CutGroup (tables are indexed by CGIdx)
        BOOST_ASSERT(forces0_array[force0_idx].count() == int(nats0));
    
        Vector *group_forces0_array = forces0.data(force0_idx);

        //ok, we are interested in the forces acting on this CutGroup
        // - calculate all of the forces on this group interacting
        //   with all of the CutGroups in mol1 
        for (quint32 jgroup=0; jgroup<ngroups1; ++jgroup)
        {
            const CoordGroup &group1 = groups1_array[jgroup];
            const Parameters::Array &params1 = molparams1_array[jgroup];

            //check first that these two CoordGroups could be within cutoff
            //(if there is only one CutGroup in both molecules then this
            //test has already been performed and passed)
            const bool within_cutoff = (ngroups0 == 1 and ngroups1 == 1) or not
                                        spce->beyond(switchfunc->cutoffDistance(), 
                                                     aabox0, group1.aaBox());
            
            if (not within_cutoff)
                //this CutGroup is beyond the cutoff distance
                continue;
            
            //calculate all of the interatomic distances
            const double mindist = spce->calcDistVectors(group0, group1, distmat);
            
            if (mindist > switchfunc->cutoffDistance())
                //all of the atoms are definitely beyond cutoff
                continue;

            const quint32 nats1 = group1.count();
            
            //loop over all interatomic pairs and calculate the energies
            const Parameter *params1_array = params1.constData();

            if (mindist > switchfunc->featherDistance())
            {
                //we need to calculate the forces taking into account
                //the derivative of the switching function!
                
                //calculate the switching scale factors and their 
                //derivatives
                const double scl_coul = switchfunc->electrostaticScaleFactor(
                                                                    Length(mindist) );
                
                Vector group_sep = (group1.aaBox().center() - aabox0.center())
                                        .normalise();
                
                Vector dscl_coul = switchfunc->dElectrostaticScaleFactor(
                                                                    Length(mindist) ) 
                                     * group_sep;
                
                double shift_coul = 0;

                if (use_electrostatic_shifting)
                    shift_coul = this->totalCharge(params0) * this->totalCharge(params1)
                                    / switchfunc->electrostaticCutoffDistance();
                
                for (quint32 i=0; i<nats0; ++i)
                {
                    distmat.setOuterIndex(i);
                    const Parameter &param0 = params0_array[i];
                
                    Vector total_force;
                
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        const double q2 = param0.reduced_charge *
                                          params1_array[j].reduced_charge;
                          
                        if (q2 != 0)
                        {
                            //calculate the coulomb energy
                            const double cnrg = q2 / distmat[j].length();
                                               
                            //calculate the coulomb force
                            Vector cforce = (scl_coul * -cnrg / distmat[j].length() *
                                             distmat[j].direction()) +
                                             
                                            ((cnrg-shift_coul) * dscl_coul);
                        
                            total_force += cforce;
                        }
                    }
                    
                    //update the forces array
                    group_forces0_array[i] += scale_force * total_force;
                }
            }
            else
            {
                //not in the feather region, so can calculate the forces
                //directly (also, no need to calculate shift, as 
                //the shifting function is constant, so does not
                //affect the gradient)
                for (quint32 i=0; i<nats0; ++i)
                {
                    distmat.setOuterIndex(i);
                    const Parameter &param0 = params0_array[i];

                    Vector total_force;
                
                    //null LJ parameter - only add on the coulomb energy
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        const double q2 = param0.reduced_charge * 
                                          params1_array[j].reduced_charge;
                        
                        //calculate the coulomb force
                        if (q2 != 0)
                        {
                            Vector cforce = -(q2 / distmat[j].length2()) *
                                                distmat[j].direction();
                        
                            total_force += cforce;
                        }
                    }
                    
                    group_forces0_array[i] += scale_force * total_force;

                } // end of loop over i atoms

            } // end of if within feather

        } // end of loop over jgroup CutGroups

    } // end of loop over igroup CutGroups
}

/** Add to the forces in 'forces0' the forces acting on 'mol0' caused
    by 'mol1' */
void InterCLJPotential::_pvt_calculateLJForce(
                                            const InterCLJPotential::Molecule &mol0, 
                                            const InterCLJPotential::Molecule &mol1,
                                            MolForceTable &forces0, 
                                            InterCLJPotential::ForceWorkspace &distmat,
                                            double scale_force) const
{
    BOOST_ASSERT( mol0.molecule().data().info().nCutGroups() == forces0.nCutGroups() );
    BOOST_ASSERT( mol0.molecule().data().number() == forces0.molNum() );

    const quint32 ngroups0 = mol0.nCutGroups();
    const quint32 ngroups1 = mol1.nCutGroups();
    
    const CoordGroup *groups0_array = mol0.coordinates().constData();
    const CoordGroup *groups1_array = mol1.coordinates().constData();
    
    BOOST_ASSERT(mol0.parameters().atomicParameters().count() == int(ngroups0));
    BOOST_ASSERT(mol1.parameters().atomicParameters().count() == int(ngroups1));
    
    const Parameters::Array *molparams0_array 
                                    = mol0.parameters().atomicParameters().constData();
    const Parameters::Array *molparams1_array
                                    = mol1.parameters().atomicParameters().constData();
    
    const MolForceTable::Array *forces0_array = forces0.constData();
    
    //loop over all pairs of CutGroups in the two molecules
    for (quint32 igroup=0; igroup<ngroups0; ++igroup)
    {
        //get the CGIdx of this group
        CGIdx cgidx_igroup = mol0.cgIdx(igroup);

        //get the index of this CutGroup in the forces array
        int force0_idx = forces0.map(cgidx_igroup);
        
        if (force0_idx == -1)
            //there is no space for the forces on this CutGroup in 
            //the forcetable - were are therefore not interested in
            //this CutGroup
            continue;

        const Parameters::Array &params0 = molparams0_array[igroup];

        const CoordGroup &group0 = groups0_array[igroup];
        const AABox &aabox0 = group0.aaBox();
        const quint32 nats0 = group0.count();
        const Parameter *params0_array = params0.constData();
    
        //get the table that holds the forces acting on all of the
        //atoms of this CutGroup (tables are indexed by CGIdx)
        BOOST_ASSERT(forces0_array[force0_idx].count() == int(nats0));
    
        Vector *group_forces0_array = forces0.data(force0_idx);

        //ok, we are interested in the forces acting on this CutGroup
        // - calculate all of the forces on this group interacting
        //   with all of the CutGroups in mol1 
        for (quint32 jgroup=0; jgroup<ngroups1; ++jgroup)
        {
            const CoordGroup &group1 = groups1_array[jgroup];
            const Parameters::Array &params1 = molparams1_array[jgroup];

            //check first that these two CoordGroups could be within cutoff
            //(if there is only one CutGroup in both molecules then this
            //test has already been performed and passed)
            const bool within_cutoff = (ngroups0 == 1 and ngroups1 == 1) or not
                                        spce->beyond(switchfunc->cutoffDistance(), 
                                                     aabox0, group1.aaBox());
            
            if (not within_cutoff)
                //this CutGroup is beyond the cutoff distance
                continue;
            
            //calculate all of the interatomic distances
            const double mindist = spce->calcDistVectors(group0, group1, distmat);
            
            if (mindist > switchfunc->cutoffDistance())
                //all of the atoms are definitely beyond cutoff
                continue;

            const quint32 nats1 = group1.count();
            
            //loop over all interatomic pairs and calculate the energies
            const Parameter *params1_array = params1.constData();

            if (mindist > switchfunc->featherDistance())
            {
                //we need to calculate the forces taking into account
                //the derivative of the switching function!
                
                //calculate the switching scale factors and their 
                //derivatives
                const double scl_lj = switchfunc->vdwScaleFactor( Length(mindist) );
                
                Vector group_sep = (group1.aaBox().center() - aabox0.center())
                                        .normalise();
                
                Vector dscl_lj = switchfunc->dVDWScaleFactor( Length(mindist) )
                                     * group_sep;
                
                for (quint32 i=0; i<nats0; ++i)
                {
                    distmat.setOuterIndex(i);
                    const Parameter &param0 = params0_array[i];
                
                    Vector total_force;
                
                    if (param0.ljid != 0)
                    {
                        for (quint32 j=0; j<nats1; ++j)
                        {
                            const Parameter &param1 = params1_array[j];

                            if (param1.ljid != 0)
                            {
                                const double invdist = double(1) / distmat[j].length();

                                const LJPair &ljpair = ljpairs.constData()[
                                                        ljpairs.map(param0.ljid,
                                                                    param1.ljid)];
                        
                                double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                                double sig_over_dist12 = pow_2(sig_over_dist6);

                                //calculate the energy
                                const double ljnrg = 4 * ljpair.epsilon() *
                                                      (sig_over_dist12 - sig_over_dist6);

                                // dU/dr requires an extra power of r
                                sig_over_dist6 *= invdist;
                                sig_over_dist12 *= invdist;

                                Vector force = ((scl_lj * 4 * ljpair.epsilon() * 
                                            (6.0*sig_over_dist6 - 12.0*sig_over_dist12))
                                               * distmat[j].direction())
                                            
                                             + (ljnrg * dscl_lj);

                                total_force += force;
                            }
                        }
                    }
                    
                    //update the forces array
                    group_forces0_array[i] += scale_force * total_force;
                }
            }
            else
            {
                //not in the feather region, so can calculate the forces
                //directly
                for (quint32 i=0; i<nats0; ++i)
                {
                    distmat.setOuterIndex(i);
                    const Parameter &param0 = params0_array[i];

                    Vector total_force;
                
                    if (param0.ljid != 0)
                    {
                        for (quint32 j=0; j<nats1; ++j)
                        {
                            //do both coulomb and LJ
                            const Parameter &param1 = params1_array[j];
                              
                            if (param1.ljid != 0)
                            {
                                const double invdist = double(1) / distmat[j].length();

                                const LJPair &ljpair = ljpairs.constData()[
                                                        ljpairs.map(param0.ljid,
                                                                    param1.ljid)];
                        
                                double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                                double sig_over_dist12 = pow_2(sig_over_dist6);

                                // dU/dr requires an extra power of r
                                sig_over_dist6 *= invdist;
                                sig_over_dist12 *= invdist;

                                Vector force = (4 * ljpair.epsilon() * 
                                                               (6.0*sig_over_dist6 - 
                                                                12.0*sig_over_dist12))
                                                * distmat[j].direction();

                                total_force += force;
                            }
                        }
                    }
                    
                    group_forces0_array[i] += scale_force * total_force;

                } // end of loop over i atoms

            } // end of if within feather

        } // end of loop over jgroup CutGroups

    } // end of loop over igroup CutGroups
}

/** Add to the potentials in 'pots0' the LJ potential acting on 'mol0' caused
    by 'mol1' */
void InterCLJPotential::_pvt_calculateLJPotential(
                                        const InterCLJPotential::Molecule &mol0, 
                                        const InterCLJPotential::Molecule &mol1,
                                        const CLJProbe &probe,
                                        MolPotentialTable &pots0, 
                                        InterCLJPotential::PotentialWorkspace &distmat,
                                        double scale_potential) const
{
    if (probe.lj().isDummy())
        return;

    BOOST_ASSERT( mol0.molecule().data().info().nCutGroups() == pots0.nCutGroups() );
    BOOST_ASSERT( mol0.molecule().data().number() == pots0.molNum() );

    const quint32 ngroups0 = mol0.nCutGroups();
    const quint32 ngroups1 = mol1.nCutGroups();
    
    const CoordGroup *groups0_array = mol0.coordinates().constData();
    const CoordGroup *groups1_array = mol1.coordinates().constData();
    
    const Parameters::Array *molparams1_array
                                    = mol1.parameters().atomicParameters().constData();
    
    //loop over all pairs of CutGroups in the two molecules
    for (quint32 igroup=0; igroup<ngroups0; ++igroup)
    {
        //get the CGIdx of this group
        CGIdx cgidx_igroup = mol0.cgIdx(igroup);

        //get the index of this CutGroup in the fields array
        int pot0_idx = pots0.map(cgidx_igroup);
        
        if (pot0_idx == -1)
            //there is no space for the potentials on this CutGroup in 
            //the potentialtable - were are therefore not interested in
            //this CutGroup
            continue;

        const CoordGroup &group0 = groups0_array[igroup];
        const AABox &aabox0 = group0.aaBox();
        const quint32 nats0 = group0.count();
    
        //get the table that holds the potentials acting at all of the
        //atoms of this CutGroup (tables are indexed by CGIdx)
        MolarEnergy *group_pots0_array = pots0.data(pot0_idx);

        //ok, we are interested in the potentials acting on this CutGroup
        // - calculate all of the potentials on this group interacting
        //   with all of the CutGroups in mol1 
        for (quint32 jgroup=0; jgroup<ngroups1; ++jgroup)
        {
            const CoordGroup &group1 = groups1_array[jgroup];
            const Parameters::Array &params1 = molparams1_array[jgroup];

            //check first that these two CoordGroups could be within cutoff
            //(if there is only one CutGroup in both molecules then this
            //test has already been performed and passed)
            const bool within_cutoff = (ngroups0 == 1 and ngroups1 == 1) or not
                                        spce->beyond(switchfunc->cutoffDistance(), 
                                                     aabox0, group1.aaBox());
            
            if (not within_cutoff)
                //this CutGroup is beyond the cutoff distance
                continue;
            
            //calculate all of the interatomic distances
            const double mindist = spce->calcDist2(group0, group1, distmat);
            
            if (mindist > switchfunc->cutoffDistance())
                //all of the atoms are definitely beyond cutoff
                continue;

            const quint32 nats1 = group1.count();
            
            //loop over all interatomic pairs and calculate the energies
            const Parameter *params1_array = params1.constData();

            if (mindist > switchfunc->featherDistance())
            {
                //calculate the switching scale factors
                const double scl_lj = switchfunc->vdwScaleFactor( Length(mindist) );

                for (quint32 i=0; i<nats0; ++i)
                {
                    distmat.setOuterIndex(i);
                
                    double total_potential = 0;
                
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        const Parameter &param1 = params1_array[j];
                        
                        if (param1.ljid != 0)
                        {
                            LJPair ljpair( ljpairs.constData()[
                                                ljpairs.map(param1.ljid,
                                                            param1.ljid)],
                                           probe.lj(),
                                           combining_rules );
                        
                            double sig_over_dist6 = pow_3(ljpair.sigma()*ljpair.sigma()*
                                                            distmat[j]);
                            double sig_over_dist12 = pow_2(sig_over_dist6);

                            //calculate the energy
                            const double ljnrg = 4 * ljpair.epsilon() *
                                                  (sig_over_dist12 - sig_over_dist6);

                            total_potential += (ljnrg * scl_lj);
                        }
                    }

                    //update the fields array
                    group_pots0_array[i] += MolarEnergy(scale_potential * 
                                                        total_potential);
                }
            }
            else
            {
                //not in the feather region, so can calculate the potential
                //directly
                for (quint32 i=0; i<nats0; ++i)
                {
                    double total_potential = 0;
                        
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        const Parameter &param1 = params1_array[j];
                        
                        if (param1.ljid != 0)
                        {
                            LJPair ljpair( ljpairs.constData()[
                                                    ljpairs.map(param1.ljid,
                                                                param1.ljid)],
                                           probe.lj(),
                                           combining_rules );
                       
                            double sig_over_dist6 = pow_6(ljpair.sigma()*
                                                          ljpair.sigma() / distmat[j]);
                            double sig_over_dist12 = pow_2(sig_over_dist6);

                            total_potential += (4 * ljpair.epsilon() * (
                                                          12.0*sig_over_dist12 -
                                                           6.0*sig_over_dist6) ); 
                        }
                    }

                    group_pots0_array[i] += MolarEnergy(scale_potential * 
                                                        total_potential);
                } // end of loop over i atoms

            } // end of if within feather

        } // end of loop over jgroup CutGroups

    } // end of loop over igroup CutGroups
}

/** Add to the potentials in 'pots0' the coulomb potential acting on 'mol0' caused
    by 'mol1' */
void InterCLJPotential::_pvt_calculateCoulombPotential(
                                        const InterCLJPotential::Molecule &mol0, 
                                        const InterCLJPotential::Molecule &mol1,
                                        const CLJProbe &probe,
                                        MolPotentialTable &pots0, 
                                        InterCLJPotential::PotentialWorkspace &distmat,
                                        double scale_potential) const
{
    if (probe.reducedCharge() == 0)
        return;

    BOOST_ASSERT( mol0.molecule().data().info().nCutGroups() == pots0.nCutGroups() );
    BOOST_ASSERT( mol0.molecule().data().number() == pots0.molNum() );

    const quint32 ngroups0 = mol0.nCutGroups();
    const quint32 ngroups1 = mol1.nCutGroups();
    
    const CoordGroup *groups0_array = mol0.coordinates().constData();
    const CoordGroup *groups1_array = mol1.coordinates().constData();
    
    const Parameters::Array *molparams1_array
                                    = mol1.parameters().atomicParameters().constData();
    
    //loop over all pairs of CutGroups in the two molecules
    for (quint32 igroup=0; igroup<ngroups0; ++igroup)
    {
        //get the CGIdx of this group
        CGIdx cgidx_igroup = mol0.cgIdx(igroup);

        //get the index of this CutGroup in the fields array
        int pot0_idx = pots0.map(cgidx_igroup);
        
        if (pot0_idx == -1)
            //there is no space for the potentials on this CutGroup in 
            //the potentialtable - were are therefore not interested in
            //this CutGroup
            continue;

        const CoordGroup &group0 = groups0_array[igroup];
        const AABox &aabox0 = group0.aaBox();
        const quint32 nats0 = group0.count();
    
        //get the table that holds the potentials acting at all of the
        //atoms of this CutGroup (tables are indexed by CGIdx)
        MolarEnergy *group_pots0_array = pots0.data(pot0_idx);

        //ok, we are interested in the potentials acting on this CutGroup
        // - calculate all of the potentials on this group interacting
        //   with all of the CutGroups in mol1 
        for (quint32 jgroup=0; jgroup<ngroups1; ++jgroup)
        {
            const CoordGroup &group1 = groups1_array[jgroup];
            const Parameters::Array &params1 = molparams1_array[jgroup];

            //check first that these two CoordGroups could be within cutoff
            //(if there is only one CutGroup in both molecules then this
            //test has already been performed and passed)
            const bool within_cutoff = (ngroups0 == 1 and ngroups1 == 1) or not
                                        spce->beyond(switchfunc->cutoffDistance(), 
                                                     aabox0, group1.aaBox());
            
            if (not within_cutoff)
                //this CutGroup is beyond the cutoff distance
                continue;
            
            //calculate all of the interatomic distances
            const double mindist = spce->calcDist(group0, group1, distmat);
            
            if (mindist > switchfunc->cutoffDistance())
                //all of the atoms are definitely beyond cutoff
                continue;

            const quint32 nats1 = group1.count();
            
            //loop over all interatomic pairs and calculate the energies
            const Parameter *params1_array = params1.constData();

            double shift_coul = 0;

            if (use_electrostatic_shifting)
                shift_coul = probe.charge() * this->totalCharge(params1)
                                / switchfunc->electrostaticCutoffDistance();

            if (mindist > switchfunc->featherDistance())
            {
                //calculate the switching scale factors
                const double scl_coul = switchfunc->electrostaticScaleFactor( 
                                                                        Length(mindist) );
                for (quint32 i=0; i<nats0; ++i)
                {
                    distmat.setOuterIndex(i);
                        
                    double total_potential = 0;
                        
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        const double q2 = probe.reducedCharge() *
                                          params1_array[j].reduced_charge;
                            
                        if (q2 != 0)
                        {
                            //calculate the coulomb energy
                            const double cnrg = q2 / distmat[j];
                                               
                            total_potential += scl_coul * (cnrg - shift_coul);
                        }
                    }

                    //update the fields array
                    group_pots0_array[i] += MolarEnergy(scale_potential * 
                                                        total_potential);
                }
            }
            else
            {
                //not in the feather region, so can calculate the potentials
                //directly
                for (quint32 i=0; i<nats0; ++i)
                {
                    double total_potential = 0;
                      
                    distmat.setOuterIndex(i);

                    for (quint32 j=0; j<nats1; ++j)
                    {
                        const double q2 = probe.reducedCharge() * 
                                          params1_array[j].reduced_charge;
                        
                        //calculate the coulomb potential
                        if (q2 != 0)
                            total_potential += (q2 / distmat[j]) - shift_coul;
                    }

                    group_pots0_array[i] += MolarEnergy(scale_potential * 
                                                            total_potential);
                }

            } // end of if within feather

        } // end of loop over jgroup CutGroups

    } // end of loop over igroup CutGroups
}

/** Add to the potentials in 'pots0' the fields acting on 'mol0' caused
    by 'mol1' */
void InterCLJPotential::_pvt_calculatePotential(
                                        const InterCLJPotential::Molecule &mol0, 
                                        const InterCLJPotential::Molecule &mol1,
                                        const CLJProbe &probe,
                                        MolPotentialTable &pots0, 
                                        InterCLJPotential::PotentialWorkspace &distmat,
                                        double scale_potential) const
{
    if (probe.lj().isDummy())
    {
        this->_pvt_calculateCoulombPotential(mol0, mol1, probe, pots0,
                                             distmat, scale_potential);
        return;
    }
    else if (probe.reducedCharge() == 0)
    {
        this->_pvt_calculateLJPotential(mol0, mol1, probe, pots0,
                                        distmat, scale_potential);
        return;
    }


    BOOST_ASSERT( mol0.molecule().data().info().nCutGroups() == pots0.nCutGroups() );
    BOOST_ASSERT( mol0.molecule().data().number() == pots0.molNum() );

    const quint32 ngroups0 = mol0.nCutGroups();
    const quint32 ngroups1 = mol1.nCutGroups();
    
    const CoordGroup *groups0_array = mol0.coordinates().constData();
    const CoordGroup *groups1_array = mol1.coordinates().constData();
    
    const Parameters::Array *molparams1_array
                                    = mol1.parameters().atomicParameters().constData();
    
    //loop over all pairs of CutGroups in the two molecules
    for (quint32 igroup=0; igroup<ngroups0; ++igroup)
    {
        //get the CGIdx of this group
        CGIdx cgidx_igroup = mol0.cgIdx(igroup);

        //get the index of this CutGroup in the fields array
        int pot0_idx = pots0.map(cgidx_igroup);
        
        if (pot0_idx == -1)
            //there is no space for the potentials on this CutGroup in 
            //the potentialtable - were are therefore not interested in
            //this CutGroup
            continue;

        const CoordGroup &group0 = groups0_array[igroup];
        const AABox &aabox0 = group0.aaBox();
        const quint32 nats0 = group0.count();
    
        //get the table that holds the potentials acting at all of the
        //atoms of this CutGroup (tables are indexed by CGIdx)
        MolarEnergy *group_pots0_array = pots0.data(pot0_idx);

        //ok, we are interested in the potentials acting on this CutGroup
        // - calculate all of the potentials on this group interacting
        //   with all of the CutGroups in mol1 
        for (quint32 jgroup=0; jgroup<ngroups1; ++jgroup)
        {
            const CoordGroup &group1 = groups1_array[jgroup];
            const Parameters::Array &params1 = molparams1_array[jgroup];

            //check first that these two CoordGroups could be within cutoff
            //(if there is only one CutGroup in both molecules then this
            //test has already been performed and passed)
            const bool within_cutoff = (ngroups0 == 1 and ngroups1 == 1) or not
                                        spce->beyond(switchfunc->cutoffDistance(), 
                                                     aabox0, group1.aaBox());
            
            if (not within_cutoff)
                //this CutGroup is beyond the cutoff distance
                continue;
            
            //calculate all of the interatomic distances
            const double mindist = spce->calcDist(group0, group1, distmat);
            
            if (mindist > switchfunc->cutoffDistance())
                //all of the atoms are definitely beyond cutoff
                continue;

            const quint32 nats1 = group1.count();
            
            //loop over all interatomic pairs and calculate the energies
            const Parameter *params1_array = params1.constData();

            double shift_coul = 0;

            if (use_electrostatic_shifting)
                shift_coul = probe.charge() * this->totalCharge(params1)
                                / switchfunc->electrostaticCutoffDistance();

            if (mindist > switchfunc->featherDistance())
            {
                //calculate the switching scale factors
                const double scl_coul = switchfunc->electrostaticScaleFactor( 
                                                                        Length(mindist) );
                const double scl_lj = switchfunc->vdwScaleFactor( Length(mindist) );

                if (probe.lj().isDummy())
                {
                    for (quint32 i=0; i<nats0; ++i)
                    {
                        distmat.setOuterIndex(i);
                        
                        double total_potential = 0;
                        
                        //null LJ probe - only add on the coulomb field
                        for (quint32 j=0; j<nats1; ++j)
                        {
                            const double q2 = probe.reducedCharge() *
                                              params1_array[j].reduced_charge;
                            
                            if (q2 != 0)
                            {
                                //calculate the coulomb energy
                                const double cnrg = q2 / distmat[j];
                                               
                                total_potential += scl_coul * (cnrg - shift_coul);
                            }
                        }

                        //update the fields array
                        group_pots0_array[i] += MolarEnergy(scale_potential * 
                                                            total_potential);
                    }
                }
                else
                {
                    for (quint32 i=0; i<nats0; ++i)
                    {
                        distmat.setOuterIndex(i);
                
                        double total_potential = 0;
                
                        for (quint32 j=0; j<nats1; ++j)
                        {
                            //do both coulomb and LJ
                            const Parameter &param1 = params1_array[j];
                        
                            const double invdist = double(1) / distmat[j];
                        
                            double potential = 0;
                            
                            const double q2 = probe.reducedCharge() * 
                                              param1.reduced_charge;
                                              
                            if (q2 != 0)
                                potential += scl_coul * (q2 * invdist - shift_coul);
                              
                            if (param1.ljid != 0)
                            {
                                LJPair ljpair( ljpairs.constData()[
                                                    ljpairs.map(param1.ljid,
                                                                param1.ljid)],
                                               probe.lj(),
                                               combining_rules );
                        
                                double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                                double sig_over_dist12 = pow_2(sig_over_dist6);

                                //calculate the energy
                                const double ljnrg = 4 * ljpair.epsilon() *
                                                      (sig_over_dist12 - sig_over_dist6);

                                potential += (ljnrg * scl_lj);
                            }

                            total_potential += potential;
                        }
                    
                        //update the fields array
                        group_pots0_array[i] += MolarEnergy(scale_potential * 
                                                            total_potential);
                    }
                }
            }
            else
            {
                //not in the feather region, so can calculate the potentials
                //directly

                if (probe.lj().isDummy())
                {
                    for (quint32 i=0; i<nats0; ++i)
                    {
                        double total_potential = 0;
                        
                        distmat.setOuterIndex(i);

                        //null LJ probe - only add on the coulomb field
                        for (quint32 j=0; j<nats1; ++j)
                        {
                            const double q2 = probe.reducedCharge() * 
                                              params1_array[j].reduced_charge;
                        
                            //calculate the coulomb potential
                            if (q2 != 0)
                                total_potential += (q2 / distmat[j]) - shift_coul;
                        }

                        group_pots0_array[i] += MolarEnergy(scale_potential * 
                                                            total_potential);
                    }
                }
                else
                {
                    for (quint32 i=0; i<nats0; ++i)
                    {
                        double total_potential = 0;
                        
                        for (quint32 j=0; j<nats1; ++j)
                        {
                            //do both coulomb and LJ
                            const Parameter &param1 = params1_array[j];
                        
                            const double invdist = double(1) / distmat[j];
                        
                            //calculate the potential
                            double potential = 0;
                            
                            const double q2 = probe.reducedCharge() * 
                                              param1.reduced_charge;
                                              
                            if (q2 != 0)
                                potential += (q2 * invdist) - shift_coul;
                              
                            if (param1.ljid != 0)
                            {
                                LJPair ljpair( ljpairs.constData()[
                                                        ljpairs.map(param1.ljid,
                                                                    param1.ljid)],
                                               probe.lj(),
                                               combining_rules );
                        
                                double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                                double sig_over_dist12 = pow_2(sig_over_dist6);

                                potential += (4 * ljpair.epsilon() * (
                                                              12.0*sig_over_dist12 -
                                                               6.0*sig_over_dist6) ); 
                            }
                        
                            total_potential += potential;
                        }

                        group_pots0_array[i] += MolarEnergy(scale_potential * 
                                                            total_potential);
                    }

                } // end of loop over i atoms

            } // end of if within feather

        } // end of loop over jgroup CutGroups

    } // end of loop over igroup CutGroups
}

/** Add to the potentials in 'pots0' the potential on the passed  
    grid caused by 'mol' */
void InterCLJPotential::_pvt_calculateLJPotential(
                                        const InterCLJPotential::Molecule &mol, 
                                        const CLJProbe &probe,
                                        GridPotentialTable &pots, 
                                        InterCLJPotential::PotentialWorkspace &distmat,
                                        double scale_potential) const
{
    if (probe.lj().isDummy())
        return;

    const int ngroups = mol.nCutGroups();
    const CoordGroup *groups_array = mol.coordinates().constData();
    const Parameters::Array *molparams_array
                                    = mol.parameters().atomicParameters().constData();

    const Grid &grid = pots.grid();
    const int npoints = grid.nPoints();
    const Vector *gridpoints_array = grid.constData();
    
    if (npoints == 0 or ngroups == 0)
        return;

    MolarEnergy *grid_pot_array = pots.data();
    
    for (int i=0; i<ngroups; ++i)
    {
        const CoordGroup &group = groups_array[i];
        const AABox &aabox = group.aaBox();
        const int nats = group.count();
        
        //check first that these two CoordGroups could be within cutoff
        //(if there is only one CutGroup in both molecules then this
        //test has already been performed and passed)
        const bool within_cutoff = (ngroups == 1) or not
                                        spce->beyond(switchfunc->cutoffDistance(), 
                                                     aabox, grid.aaBox());
            
        if (not within_cutoff)
            //this CutGroup is beyond the cutoff distance
            continue;

        const Parameters::Array &params = molparams_array[i];
        const Parameter *params_array = params.constData();
        
        for (int j=0; j<npoints; ++j)
        {
            const Vector &gridpoint = gridpoints_array[j];
            
            const double mindist = spce->calcDist2(group, gridpoint, distmat);

            double total_potential = 0;
            
            if (mindist > switchfunc->cutoffDistance())
                continue;
                
            else if (mindist > switchfunc->featherDistance())
            {
                //we need to calculate the field taking into account
                //the derivative of the switching function!
            
                //calculate the switching scale factors and their 
                //derivatives
                const double scl_lj = switchfunc->vdwScaleFactor( Length(mindist) );
            
                distmat.setOuterIndex(0);
        
                for (int k=0; k<nats; ++k)
                {
                    //do both coulomb and LJ
                    const Parameter &param = params_array[k];
                
                    if (param.ljid != 0)
                    {
                        LJPair ljpair( ljpairs.constData()[
                                            ljpairs.map(param.ljid,
                                                        param.ljid)],
                                       probe.lj(),
                                       combining_rules );
                
                        double sig_over_dist6 = pow_3(ljpair.sigma()*
                                                      ljpair.sigma() / 
                                                      distmat[j]);
                                                      
                        double sig_over_dist12 = pow_2(sig_over_dist6);

                        //calculate the energy
                        total_potential +=  scl_lj * 4 * ljpair.epsilon() *
                                                    (sig_over_dist12 - sig_over_dist6);
                    }
                
                } // end of loop over atoms
            }
            else
            {
                //no need to worry about the switching function :-)
                distmat.setOuterIndex(0);
            
                for (int k=0; k<nats; ++k)
                {
                    const Parameter &param = params_array[k];
                
                    if (param.ljid != 0)
                    {
                        LJPair ljpair( ljpairs.constData()[
                                                ljpairs.map(param.ljid,
                                                            param.ljid)],
                                       probe.lj(),
                                       combining_rules );
                
                        double sig_over_dist6 = pow_3(ljpair.sigma()*
                                                      ljpair.sigma()/
                                                      distmat[j]);
                        double sig_over_dist12 = pow_2(sig_over_dist6);

                        total_potential += 4 * ljpair.epsilon() * (
                                                   12.0*sig_over_dist12 -
                                                    6.0*sig_over_dist6 );
                    }

                } // end of loop over atoms
            }
            
            grid_pot_array[j] += MolarEnergy(scale_potential * total_potential);

        } // end of loop over grid points
    } // end of loop over CutGroups
}

/** Add to the potentials in 'pots0' the coulomb potential on the passed  
    grid caused by 'mol' */
void InterCLJPotential::_pvt_calculateCoulombPotential(
                                        const InterCLJPotential::Molecule &mol, 
                                        const CLJProbe &probe,
                                        GridPotentialTable &pots, 
                                        InterCLJPotential::PotentialWorkspace &distmat,
                                        double scale_potential) const
{
    if (probe.reducedCharge() == 0)
        return;

    const int ngroups = mol.nCutGroups();
    const CoordGroup *groups_array = mol.coordinates().constData();
    const Parameters::Array *molparams_array
                                    = mol.parameters().atomicParameters().constData();

    const Grid &grid = pots.grid();
    const int npoints = grid.nPoints();
    const Vector *gridpoints_array = grid.constData();
    
    if (npoints == 0 or ngroups == 0)
        return;

    MolarEnergy *grid_pot_array = pots.data();
    
    for (int i=0; i<ngroups; ++i)
    {
        const CoordGroup &group = groups_array[i];
        const AABox &aabox = group.aaBox();
        const int nats = group.count();
        
        //check first that these two CoordGroups could be within cutoff
        //(if there is only one CutGroup in both molecules then this
        //test has already been performed and passed)
        const bool within_cutoff = (ngroups == 1) or not
                                        spce->beyond(switchfunc->cutoffDistance(), 
                                                     aabox, grid.aaBox());
            
        if (not within_cutoff)
            //this CutGroup is beyond the cutoff distance
            continue;

        const Parameters::Array &params = molparams_array[i];
        const Parameter *params_array = params.constData();
        
        for (int j=0; j<npoints; ++j)
        {
            const Vector &gridpoint = gridpoints_array[j];
            
            const double mindist = spce->calcDist(group, gridpoint, distmat);

            double total_potential = 0;
            
            if (mindist > switchfunc->cutoffDistance())
                continue;
                
            else if (mindist > switchfunc->featherDistance())
            {
                //we need to calculate the field taking into account
                //the derivative of the switching function!
            
                //calculate the switching scale factors
                const double scl_coul = switchfunc->electrostaticScaleFactor( 
                                                                    Length(mindist) );
            
                double shift_coul = 0;

                if (use_electrostatic_shifting)
                    shift_coul = probe.charge() * this->totalCharge(params)
                                    / switchfunc->electrostaticCutoffDistance();

                distmat.setOuterIndex(0);
        
                for (int k=0; k<nats; ++k)
                {
                    //do both coulomb and LJ
                    const Parameter &param = params_array[k];
                
                    const double invdist = double(1) / distmat[k];
                
                    const double q2 = probe.reducedCharge() *
                                      param.reduced_charge;
                
                    if (q2 != 0)
                        total_potential += scl_coul * (q2 * invdist - shift_coul);
                      
                } // end of loop over atoms
            }
            else
            {
                //no need to worry about the switching function :-)
                distmat.setOuterIndex(0);
            
                double shift_coul = 0;

                if (use_electrostatic_shifting)
                    shift_coul = probe.charge() * this->totalCharge(params)
                                    / switchfunc->electrostaticCutoffDistance();
        
                for (int k=0; k<nats; ++k)
                {
                    const Parameter &param = params_array[k];
                
                    const double invdist = double(1) / distmat[k];
                
                    const double q2 = probe.reducedCharge() * param.reduced_charge;
                    
                    if (q2 != 0)
                        total_potential += q2 * invdist - shift_coul;

                } // end of loop over atoms
            }
            
            grid_pot_array[j] += MolarEnergy(scale_potential * total_potential);

        } // end of loop over grid points
    } // end of loop over CutGroups
}

/** Add to the potentials in 'pots0' the potential on the passed  
    grid caused by 'mol' */
void InterCLJPotential::_pvt_calculatePotential(
                                        const InterCLJPotential::Molecule &mol, 
                                        const CLJProbe &probe,
                                        GridPotentialTable &pots, 
                                        InterCLJPotential::PotentialWorkspace &distmat,
                                        double scale_potential) const
{
    if (probe.lj().isDummy())
    {
        this->_pvt_calculateCoulombPotential(mol, probe, pots, distmat, scale_potential);
        return;
    }
    else if (probe.reducedCharge() == 0)
    {
        this->_pvt_calculateLJPotential(mol, probe, pots, distmat, scale_potential);
        return;
    }

    const int ngroups = mol.nCutGroups();
    const CoordGroup *groups_array = mol.coordinates().constData();
    const Parameters::Array *molparams_array
                                    = mol.parameters().atomicParameters().constData();

    const Grid &grid = pots.grid();
    const int npoints = grid.nPoints();
    const Vector *gridpoints_array = grid.constData();
    
    if (npoints == 0 or ngroups == 0)
        return;

    MolarEnergy *grid_pot_array = pots.data();
    
    for (int i=0; i<ngroups; ++i)
    {
        const CoordGroup &group = groups_array[i];
        const AABox &aabox = group.aaBox();
        const int nats = group.count();
        
        //check first that these two CoordGroups could be within cutoff
        //(if there is only one CutGroup in both molecules then this
        //test has already been performed and passed)
        const bool within_cutoff = (ngroups == 1) or not
                                        spce->beyond(switchfunc->cutoffDistance(), 
                                                     aabox, grid.aaBox());
            
        if (not within_cutoff)
            //this CutGroup is beyond the cutoff distance
            continue;

        const Parameters::Array &params = molparams_array[i];
        const Parameter *params_array = params.constData();
        
        for (int j=0; j<npoints; ++j)
        {
            const Vector &gridpoint = gridpoints_array[j];
            
            const double mindist = spce->calcDist(group, gridpoint, distmat);

            double total_potential = 0;
            
            if (mindist > switchfunc->cutoffDistance())
                continue;
                
            else if (mindist > switchfunc->featherDistance())
            {
                //we need to calculate the field taking into account
                //the derivative of the switching function!
            
                //calculate the switching scale factors and their 
                //derivatives
                const double scl_coul = switchfunc->electrostaticScaleFactor( 
                                                                    Length(mindist) );
                const double scl_lj = switchfunc->vdwScaleFactor( Length(mindist) );
            
                double shift_coul = 0;

                if (use_electrostatic_shifting)
                    shift_coul = probe.charge() * this->totalCharge(params)
                                    / switchfunc->electrostaticCutoffDistance();

                distmat.setOuterIndex(0);
        
                for (int k=0; k<nats; ++k)
                {
                    //do both coulomb and LJ
                    const Parameter &param = params_array[k];
                
                    const double invdist = double(1) / distmat[k];
                
                    double potential = 0;
                    
                    const double q2 = probe.reducedCharge() *
                                      param.reduced_charge;
                
                    if (q2 != 0)
                        potential += scl_coul * (q2 * invdist - shift_coul);
                      
                    if (param.ljid != 0)
                    {
                        LJPair ljpair( ljpairs.constData()[
                                            ljpairs.map(param.ljid,
                                                        param.ljid)],
                                       probe.lj(),
                                       combining_rules );
                
                        double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                        double sig_over_dist12 = pow_2(sig_over_dist6);

                        //calculate the energy
                        potential +=  scl_lj * 4 * ljpair.epsilon() *
                                              (sig_over_dist12 - sig_over_dist6);
                    }

                    total_potential += potential;
                
                } // end of loop over atoms
            }
            else
            {
                //no need to worry about the switching function :-)
                distmat.setOuterIndex(0);
            
                double shift_coul = 0;

                if (use_electrostatic_shifting)
                    shift_coul = probe.charge() * this->totalCharge(params)
                                    / switchfunc->electrostaticCutoffDistance();
        
                for (int k=0; k<nats; ++k)
                {
                    const Parameter &param = params_array[k];
                
                    const double invdist = double(1) / distmat[k];
                
                    //calculate the potential
                    double potential = 0;
                    
                    const double q2 = probe.reducedCharge() * param.reduced_charge;
                    
                    if (q2 != 0)
                        potential += q2 * invdist - shift_coul;
                      
                    if (param.ljid != 0)
                    {
                        LJPair ljpair( ljpairs.constData()[
                                                ljpairs.map(param.ljid,
                                                            param.ljid)],
                                       probe.lj(),
                                       combining_rules );
                
                        double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                        double sig_over_dist12 = pow_2(sig_over_dist6);

                        potential += 4 * ljpair.epsilon() * (
                                                12.0*sig_over_dist12 -
                                                 6.0*sig_over_dist6 );
                    }

                    total_potential += potential;

                } // end of loop over atoms
            }
            
            grid_pot_array[j] += MolarEnergy(scale_potential * total_potential);

        } // end of loop over grid points
    } // end of loop over CutGroups
}

/** Add to the fields in 'fields0' the fields acting on 'mol0' caused
    by 'mol1' */
void InterCLJPotential::_pvt_calculateField(const InterCLJPotential::Molecule &mol0, 
                                            const InterCLJPotential::Molecule &mol1,
                                            const CLJProbe &probe,
                                            MolFieldTable &fields0, 
                                            InterCLJPotential::FieldWorkspace &distmat,
                                            double scale_field) const
{
    BOOST_ASSERT( mol0.molecule().data().info().nCutGroups() == fields0.nCutGroups() );
    BOOST_ASSERT( mol0.molecule().data().number() == fields0.molNum() );

    const quint32 ngroups0 = mol0.nCutGroups();
    const quint32 ngroups1 = mol1.nCutGroups();
    
    const CoordGroup *groups0_array = mol0.coordinates().constData();
    const CoordGroup *groups1_array = mol1.coordinates().constData();
    
    const Parameters::Array *molparams1_array
                                    = mol1.parameters().atomicParameters().constData();
    
    const MolFieldTable::Array *fields0_array = fields0.constData();
    
    //loop over all pairs of CutGroups in the two molecules
    for (quint32 igroup=0; igroup<ngroups0; ++igroup)
    {
        //get the CGIdx of this group
        CGIdx cgidx_igroup = mol0.cgIdx(igroup);

        //get the index of this CutGroup in the fields array
        int field0_idx = fields0.map(cgidx_igroup);
        
        if (field0_idx == -1)
            //there is no space for the fields on this CutGroup in 
            //the fieldtable - were are therefore not interested in
            //this CutGroup
            continue;

        const CoordGroup &group0 = groups0_array[igroup];
        const AABox &aabox0 = group0.aaBox();
        const quint32 nats0 = group0.count();
    
        //get the table that holds the fields acting at all of the
        //atoms of this CutGroup (tables are indexed by CGIdx)
        BOOST_ASSERT(fields0_array[field0_idx].count() == int(nats0));
    
        Vector *group_fields0_array = fields0.data(field0_idx);

        //ok, we are interested in the fields acting on this CutGroup
        // - calculate all of the fields on this group interacting
        //   with all of the CutGroups in mol1 
        for (quint32 jgroup=0; jgroup<ngroups1; ++jgroup)
        {
            const CoordGroup &group1 = groups1_array[jgroup];
            const Parameters::Array &params1 = molparams1_array[jgroup];

            //check first that these two CoordGroups could be within cutoff
            //(if there is only one CutGroup in both molecules then this
            //test has already been performed and passed)
            const bool within_cutoff = (ngroups0 == 1 and ngroups1 == 1) or not
                                        spce->beyond(switchfunc->cutoffDistance(), 
                                                     aabox0, group1.aaBox());
            
            if (not within_cutoff)
                //this CutGroup is beyond the cutoff distance
                continue;
            
            //calculate all of the interatomic distances
            const double mindist = spce->calcDistVectors(group0, group1, distmat);
            
            if (mindist > switchfunc->cutoffDistance())
                //all of the atoms are definitely beyond cutoff
                continue;

            const quint32 nats1 = group1.count();
            
            //loop over all interatomic pairs and calculate the energies
            const Parameter *params1_array = params1.constData();

            //note that we are ignoring the switching function when calculating
            //fields...

            if (mindist > switchfunc->featherDistance())
            {
                //we need to calculate the forces taking into account
                //the derivative of the switching function!
                
                //calculate the switching scale factors and their 
                //derivatives
                const double scl_coul = switchfunc->electrostaticScaleFactor( 
                                                                        Length(mindist) );
                const double scl_lj = switchfunc->vdwScaleFactor( Length(mindist) );
                
                Vector group_sep = (group1.aaBox().center() - aabox0.center())
                                        .normalise();
                
                Vector dscl_coul = switchfunc->dElectrostaticScaleFactor( 
                                                                    Length(mindist) ) 
                                     * group_sep;
                                     
                Vector dscl_lj = switchfunc->dVDWScaleFactor( Length(mindist) )
                                     * group_sep;
                
                double shift_coul = 0;

                if (use_electrostatic_shifting)
                    shift_coul = probe.charge() * this->totalCharge(params1)
                                    / switchfunc->electrostaticCutoffDistance();

                if (probe.lj().isDummy())
                {
                    for (quint32 i=0; i<nats0; ++i)
                    {
                        distmat.setOuterIndex(i);
                        
                        Vector field;
                        
                        //null LJ probe - only add on the coulomb field
                        for (quint32 j=0; j<nats1; ++j)
                        {
                            const double q2 = probe.reducedCharge() *
                                              params1_array[j].reduced_charge;
                            
                            if (q2 != 0)
                            {
                                //calculate the coulomb energy
                                const double cnrg = q2 / distmat[j].length();
                                               
                                //calculate the coulomb field
                                Vector cfield = (scl_coul * -cnrg / distmat[j].length() *
                                                 distmat[j].direction()) +
                                             
                                                ((cnrg-shift_coul) * dscl_coul);
                        
                                field += cfield;
                            }
                        }
                    }
                }
                else
                {
                    for (quint32 i=0; i<nats0; ++i)
                    {
                        distmat.setOuterIndex(i);
                
                        Vector total_field;
                
                        for (quint32 j=0; j<nats1; ++j)
                        {
                            //do both coulomb and LJ
                            const Parameter &param1 = params1_array[j];
                        
                            const double invdist = double(1) / distmat[j].length();
                        
                            Vector field;
                            
                            const double q2 = probe.reducedCharge() *
                                              param1.reduced_charge;
                        
                            if (q2 != 0)
                            {
                                //calculate the energy
                                const double cnrg = q2 * invdist;
                        
                                //calculate the field
                                field = (scl_coul * -cnrg / distmat[j].length() *
                                         distmat[j].direction()) +
                                             
                                         ((cnrg-shift_coul) * dscl_coul);
                            }
                              
                            if (param1.ljid != 0)
                            {
                                LJPair ljpair( ljpairs.constData()[
                                                    ljpairs.map(param1.ljid,
                                                                param1.ljid)],
                                               probe.lj(),
                                               combining_rules );
                        
                                double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                                double sig_over_dist12 = pow_2(sig_over_dist6);

                                //calculate the energy
                                const double ljnrg = 4 * ljpair.epsilon() *
                                                      (sig_over_dist12 - sig_over_dist6);

                                // dU/dr requires an extra power of r
                                sig_over_dist6 *= invdist;
                                sig_over_dist12 *= invdist;

                                field += ((scl_lj * ljpair.epsilon() * 
                                            (6.0*sig_over_dist6 - 12.0*sig_over_dist12))
                                            * distmat[j].direction())
                                            
                                          + (ljnrg * dscl_lj);
                            }

                            total_field += field;
                        }
                    
                        //update the fields array
                        group_fields0_array[i] += scale_field * total_field;
                    }
                }
            }
            else
            {
                //not in the feather region, so can calculate the fields
                //directly (also, no need to calculate shift, as 
                //the shifting function is constant, so does not
                //affect the gradient)
                Vector total_field;

                if (probe.lj().isDummy())
                {
                    for (quint32 i=0; i<nats0; ++i)
                    {
                        distmat.setOuterIndex(i);

                        //null LJ probe - only add on the coulomb field
                        for (quint32 j=0; j<nats1; ++j)
                        {
                            const double q2 = probe.reducedCharge() * 
                                              params1_array[j].reduced_charge;
                        
                            //calculate the coulomb field
                            if (q2 != 0)
                            {
                                Vector cfield = -(q2 / distmat[j].length2()) *
                                                    distmat[j].direction();
                        
                                total_field += cfield;
                            }
                        }
                    }
                }
                else
                {
                    for (quint32 i=0; i<nats0; ++i)
                    {
                        for (quint32 j=0; j<nats1; ++j)
                        {
                            //do both coulomb and LJ
                            const Parameter &param1 = params1_array[j];
                        
                            const double invdist = double(1) / distmat[j].length();
                            const double invdist2 = pow_2(invdist);
                        
                            //calculate the field
                            Vector field = -(probe.reducedCharge() * 
                                             param1.reduced_charge * invdist2) 
                                            
                                            * distmat[j].direction();
                              
                            if (param1.ljid != 0)
                            {
                                LJPair ljpair( ljpairs.constData()[
                                                        ljpairs.map(param1.ljid,
                                                                    param1.ljid)],
                                               probe.lj(),
                                               combining_rules );
                        
                                double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                                double sig_over_dist12 = pow_2(sig_over_dist6);

                                // dU/dr requires an extra power of r
                                sig_over_dist6 *= invdist;
                                sig_over_dist12 *= invdist;

                                field += (4 * ljpair.epsilon() * (6.0*sig_over_dist6 - 
                                                              12.0*sig_over_dist12))
                                        * distmat[j].direction();
                            }
                        
                            total_field += field;
                        }
                    
                        group_fields0_array[i] += scale_field * total_field;
                    }

                } // end of loop over i atoms

            } // end of if within feather

        } // end of loop over jgroup CutGroups

    } // end of loop over igroup CutGroups
}

/** Calculate the field caused by the molecule 'mol' on the grid points in 
    'fields' */
void InterCLJPotential::_pvt_calculateLJField(const InterCLJPotential::Molecule &mol,
                                              const CLJProbe &probe,
                                              GridFieldTable &fields,
                                              InterCLJPotential::FieldWorkspace &distmat,
                                              double scale_field) const
{
    if (probe.lj().isDummy())
        return;

    const int ngroups = mol.nCutGroups();
    const CoordGroup *groups_array = mol.coordinates().constData();
    const Parameters::Array *molparams_array
                                    = mol.parameters().atomicParameters().constData();

    const Grid &grid = fields.grid();
    const int npoints = grid.nPoints();
    const Vector *gridpoints_array = grid.constData();
    Vector *grid_field_array = fields.data();
    
    if (npoints == 0 or ngroups == 0)
        return;
    
    for (int i=0; i<ngroups; ++i)
    {
        const CoordGroup &group = groups_array[i];
        const AABox &aabox = group.aaBox();
        const int nats = group.count();
        
        //check first that these two CoordGroups could be within cutoff
        //(if there is only one CutGroup in both molecules then this
        //test has already been performed and passed)
        const bool within_cutoff = (ngroups == 1) or not
                                        spce->beyond(switchfunc->cutoffDistance(), 
                                                     aabox, grid.aaBox());
            
        if (not within_cutoff)
            //this CutGroup is beyond the cutoff distance
            continue;

        const Parameters::Array &params = molparams_array[i];
        const Parameter *params_array = params.constData();
        
        for (int j=0; j<npoints; ++j)
        {
            const Vector &gridpoint = gridpoints_array[j];
            
            const double mindist = spce->calcDistVectors(group, gridpoint, distmat);

            Vector total_field;
            
            if (mindist > switchfunc->cutoffDistance())
                continue;
                
            else if (mindist > switchfunc->featherDistance())
            {
                //we need to calculate the field taking into account
                //the derivative of the switching function!
            
                //calculate the switching scale factors and their 
                //derivatives
                const double scl_lj = switchfunc->vdwScaleFactor( Length(mindist) );
            
                Vector group_sep = (group.aaBox().center() - gridpoint).normalise();

                Vector dscl_lj = switchfunc->dVDWScaleFactor( Length(mindist) )
                                                    * group_sep;

                distmat.setOuterIndex(0);
        
                for (int k=0; k<nats; ++k)
                {
                    //do both coulomb and LJ
                    const Parameter &param = params_array[k];
                
                    const double invdist = double(1) / distmat[k].length();
                
                    if (param.ljid != 0)
                    {
                        LJPair ljpair( ljpairs.constData()[
                                            ljpairs.map(param.ljid,
                                                        param.ljid)],
                                       probe.lj(),
                                       combining_rules );
                
                        double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                        double sig_over_dist12 = pow_2(sig_over_dist6);

                        //calculate the energy
                        const double ljnrg = 4 * ljpair.epsilon() *
                                              (sig_over_dist12 - sig_over_dist6);

                        // dU/dr requires an extra power of r
                        sig_over_dist6 *= invdist;
                        sig_over_dist12 *= invdist;

                        Vector field = ((scl_lj * 4 * ljpair.epsilon() * 
                                        (6.0*sig_over_dist6 - 12.0*sig_over_dist12))
                                        * distmat[j].direction())
                                     
                                        + (ljnrg * dscl_lj);

                        total_field += field;
                    }
                
                } // end of loop over atoms
            }
            else
            {
                //no need to worry about the switching function :-)
                distmat.setOuterIndex(0);
        
                for (int k=0; k<nats; ++k)
                {
                    const Parameter &param = params_array[k];
                
                    const double invdist = double(1) / distmat[k].length();
                
                    if (param.ljid != 0)
                    {
                        LJPair ljpair( ljpairs.constData()[
                                                ljpairs.map(param.ljid,
                                                            param.ljid)],
                                       probe.lj(),
                                       combining_rules );
                
                        double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                        double sig_over_dist12 = pow_2(sig_over_dist6);

                        // dU/dr requires an extra power of r
                        sig_over_dist6 *= invdist;
                        sig_over_dist12 *= invdist;

                        Vector field = (4 * ljpair.epsilon() * (6.0*sig_over_dist6 - 
                                                               12.0*sig_over_dist12))
                                          * distmat[j].direction();
                        
                        total_field += field;
                    }

                } // end of loop over atoms
            }
            
            grid_field_array[j] += scale_field * total_field;

        } // end of loop over grid points
    } // end of loop over CutGroups
}

/** Calculate the coulomb field caused by the molecule 'mol' on the grid points in 
    'fields' */
void InterCLJPotential::_pvt_calculateCoulombField(
                                            const InterCLJPotential::Molecule &mol,
                                            const CLJProbe &probe,
                                            GridFieldTable &fields,
                                            InterCLJPotential::FieldWorkspace &distmat,
                                            double scale_field) const
{
    if (probe.reducedCharge() == 0)
        return;

    const int ngroups = mol.nCutGroups();
    const CoordGroup *groups_array = mol.coordinates().constData();
    const Parameters::Array *molparams_array
                                    = mol.parameters().atomicParameters().constData();

    const Grid &grid = fields.grid();
    const int npoints = grid.nPoints();
    const Vector *gridpoints_array = grid.constData();
    Vector *grid_field_array = fields.data();
    
    if (npoints == 0 or ngroups == 0)
        return;
    
    for (int i=0; i<ngroups; ++i)
    {
        const CoordGroup &group = groups_array[i];
        const AABox &aabox = group.aaBox();
        const int nats = group.count();
        
        //check first that these two CoordGroups could be within cutoff
        //(if there is only one CutGroup in both molecules then this
        //test has already been performed and passed)
        const bool within_cutoff = (ngroups == 1) or not
                                        spce->beyond(switchfunc->cutoffDistance(), 
                                                     aabox, grid.aaBox());
            
        if (not within_cutoff)
            //this CutGroup is beyond the cutoff distance
            continue;

        const Parameters::Array &params = molparams_array[i];
        const Parameter *params_array = params.constData();
        
        for (int j=0; j<npoints; ++j)
        {
            const Vector &gridpoint = gridpoints_array[j];
            
            const double mindist = spce->calcDistVectors(group, gridpoint, distmat);

            Vector total_field;
            
            if (mindist > switchfunc->cutoffDistance())
                continue;
                
            else if (mindist > switchfunc->featherDistance())
            {
                //we need to calculate the field taking into account
                //the derivative of the switching function!
            
                //calculate the switching scale factors and their 
                //derivatives
                const double scl_coul = switchfunc->electrostaticScaleFactor( 
                                                                    Length(mindist) );
            
                Vector group_sep = (group.aaBox().center() - gridpoint).normalise();

                Vector dscl_coul = switchfunc->dElectrostaticScaleFactor( 
                                                                    Length(mindist) ) 
                                                * group_sep;
                                 
                double shift_coul = 0;

                if (use_electrostatic_shifting)
                    shift_coul = probe.charge() * this->totalCharge(params)
                                    / switchfunc->electrostaticCutoffDistance();

                distmat.setOuterIndex(0);
        
                for (int k=0; k<nats; ++k)
                {
                    //do both coulomb and LJ
                    const Parameter &param = params_array[k];
                
                    const double invdist = double(1) / distmat[k].length();
                
                    const double q2 = probe.reducedCharge() *
                                      param.reduced_charge;
                
                    if (q2 != 0)
                    {
                        //calculate the energy
                        const double cnrg = q2 * invdist;
                
                        //calculate the field
                        Vector field = (scl_coul * -cnrg / distmat[j].length() *
                                        distmat[j].direction()) +
                                     
                                       ((cnrg-shift_coul) * dscl_coul);

                        total_field += field;
                    }
                
                } // end of loop over atoms
            }
            else
            {
                //no need to worry about the switching function :-)
                distmat.setOuterIndex(0);
        
                for (int k=0; k<nats; ++k)
                {
                    const Parameter &param = params_array[k];
                
                    const double invdist = double(1) / distmat[k].length();
                    const double invdist2 = pow_2(invdist);
                
                    //calculate the field
                    Vector field = -(probe.reducedCharge() * 
                                     param.reduced_charge * invdist2) 
                                    
                                    * distmat[j].direction();

                    total_field += field;

                } // end of loop over atoms
            }
            
            grid_field_array[j] += scale_field * total_field;

        } // end of loop over grid points
    } // end of loop over CutGroups
}

/** Calculate the field caused by the molecule 'mol' on the grid points in 
    'fields' */
void InterCLJPotential::_pvt_calculateField(const InterCLJPotential::Molecule &mol,
                                            const CLJProbe &probe,
                                            GridFieldTable &fields,
                                            InterCLJPotential::FieldWorkspace &distmat,
                                            double scale_field) const
{
    if (probe.lj().isDummy())
    {
        this->_pvt_calculateCoulombField(mol, probe, fields, distmat, scale_field);
        return;
    }
    else if (probe.reducedCharge() == 0)
    {
        this->_pvt_calculateLJField(mol, probe, fields, distmat, scale_field);
        return;
    }

    const int ngroups = mol.nCutGroups();
    const CoordGroup *groups_array = mol.coordinates().constData();
    const Parameters::Array *molparams_array
                                    = mol.parameters().atomicParameters().constData();

    const Grid &grid = fields.grid();
    const int npoints = grid.nPoints();
    const Vector *gridpoints_array = grid.constData();
    Vector *grid_field_array = fields.data();
    
    if (npoints == 0 or ngroups == 0)
        return;
    
    for (int i=0; i<ngroups; ++i)
    {
        const CoordGroup &group = groups_array[i];
        const AABox &aabox = group.aaBox();
        const int nats = group.count();
        
        //check first that these two CoordGroups could be within cutoff
        //(if there is only one CutGroup in both molecules then this
        //test has already been performed and passed)
        const bool within_cutoff = (ngroups == 1) or not
                                        spce->beyond(switchfunc->cutoffDistance(), 
                                                     aabox, grid.aaBox());
            
        if (not within_cutoff)
            //this CutGroup is beyond the cutoff distance
            continue;

        const Parameters::Array &params = molparams_array[i];
        const Parameter *params_array = params.constData();
        
        for (int j=0; j<npoints; ++j)
        {
            const Vector &gridpoint = gridpoints_array[j];
            
            const double mindist = spce->calcDistVectors(group, gridpoint, distmat);

            Vector total_field;
            
            if (mindist > switchfunc->cutoffDistance())
                continue;
                
            else if (mindist > switchfunc->featherDistance())
            {
                //we need to calculate the field taking into account
                //the derivative of the switching function!
            
                //calculate the switching scale factors and their 
                //derivatives
                const double scl_coul = switchfunc->electrostaticScaleFactor( 
                                                                    Length(mindist) );
                const double scl_lj = switchfunc->vdwScaleFactor( Length(mindist) );
            
                Vector group_sep = (group.aaBox().center() - gridpoint).normalise();

                Vector dscl_coul = switchfunc->dElectrostaticScaleFactor( 
                                                                    Length(mindist) ) 
                                                * group_sep;
                                 
                Vector dscl_lj = switchfunc->dVDWScaleFactor( Length(mindist) )
                                                    * group_sep;
            
                double shift_coul = 0;

                if (use_electrostatic_shifting)
                    shift_coul = probe.charge() * this->totalCharge(params)
                                    / switchfunc->electrostaticCutoffDistance();

                distmat.setOuterIndex(0);
        
                for (int k=0; k<nats; ++k)
                {
                    //do both coulomb and LJ
                    const Parameter &param = params_array[k];
                
                    const double invdist = double(1) / distmat[k].length();
                
                    Vector field;
                    
                    const double q2 = probe.reducedCharge() *
                                      param.reduced_charge;
                
                    if (q2 != 0)
                    {
                        //calculate the energy
                        const double cnrg = q2 * invdist;
                
                        //calculate the field
                        field = (scl_coul * -cnrg / distmat[j].length() *
                                 distmat[j].direction()) +
                                     
                                 ((cnrg-shift_coul) * dscl_coul);
                    }
                      
                    if (param.ljid != 0)
                    {
                        LJPair ljpair( ljpairs.constData()[
                                            ljpairs.map(param.ljid,
                                                        param.ljid)],
                                       probe.lj(),
                                       combining_rules );
                
                        double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                        double sig_over_dist12 = pow_2(sig_over_dist6);

                        //calculate the energy
                        const double ljnrg = 4 * ljpair.epsilon() *
                                              (sig_over_dist12 - sig_over_dist6);

                        // dU/dr requires an extra power of r
                        sig_over_dist6 *= invdist;
                        sig_over_dist12 *= invdist;

                        field += ((scl_lj * 4 * ljpair.epsilon() * 
                                    (6.0*sig_over_dist6 - 12.0*sig_over_dist12))
                                    * distmat[j].direction())
                                    
                                  + (ljnrg * dscl_lj);
                    }

                    total_field += field;
                
                } // end of loop over atoms
            }
            else
            {
                //no need to worry about the switching function :-)
                distmat.setOuterIndex(0);
        
                for (int k=0; k<nats; ++k)
                {
                    const Parameter &param = params_array[k];
                
                    const double invdist = double(1) / distmat[k].length();
                    const double invdist2 = pow_2(invdist);
                
                    //calculate the field
                    Vector field = -(probe.reducedCharge() * 
                                     param.reduced_charge * invdist2) 
                                    
                                    * distmat[j].direction();
                      
                    if (param.ljid != 0)
                    {
                        LJPair ljpair( ljpairs.constData()[
                                                ljpairs.map(param.ljid,
                                                            param.ljid)],
                                       probe.lj(),
                                       combining_rules );
                
                        double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                        double sig_over_dist12 = pow_2(sig_over_dist6);

                        // dU/dr requires an extra power of r
                        sig_over_dist6 *= invdist;
                        sig_over_dist12 *= invdist;

                        field += (4 * ljpair.epsilon() * (6.0*sig_over_dist6 - 
                                                      12.0*sig_over_dist12))
                                * distmat[j].direction();
                    }

                    total_field += field;

                } // end of loop over atoms
            }
            
            grid_field_array[j] += scale_field * total_field;

        } // end of loop over grid points
    } // end of loop over CutGroups
}

/** Add to the fields in 'fields0' the fields acting on the passed probe
    at the atom points in 'mol0' caused by 'mol1' */
void InterCLJPotential::_pvt_calculateCoulombField(
                                            const InterCLJPotential::Molecule &mol0, 
                                            const InterCLJPotential::Molecule &mol1,
                                            const CLJProbe &probe,
                                            MolFieldTable &fields0, 
                                            InterCLJPotential::FieldWorkspace &distmat,
                                            double scale_field) const
{
    BOOST_ASSERT( mol0.molecule().data().info().nCutGroups() == fields0.nCutGroups() );
    BOOST_ASSERT( mol0.molecule().data().number() == fields0.molNum() );

    if (probe.reducedCharge() == 0)
        return;

    const quint32 ngroups0 = mol0.nCutGroups();
    const quint32 ngroups1 = mol1.nCutGroups();
    
    const CoordGroup *groups0_array = mol0.coordinates().constData();
    const CoordGroup *groups1_array = mol1.coordinates().constData();
    
    BOOST_ASSERT(mol1.parameters().atomicParameters().count() == int(ngroups1));
    
    const Parameters::Array *molparams1_array
                                    = mol1.parameters().atomicParameters().constData();
    
    const MolFieldTable::Array *fields0_array = fields0.constData();
    
    //loop over all pairs of CutGroups in the two molecules
    for (quint32 igroup=0; igroup<ngroups0; ++igroup)
    {
        //get the CGIdx of this group
        CGIdx cgidx_igroup = mol0.cgIdx(igroup);

        //get the index of this CutGroup in the fields array
        int field0_idx = fields0.map(cgidx_igroup);
        
        if (field0_idx == -1)
            //there is no space for the fields on this CutGroup in 
            //the fieldtable - were are therefore not interested in
            //this CutGroup
            continue;

        const CoordGroup &group0 = groups0_array[igroup];
        const AABox &aabox0 = group0.aaBox();
        const quint32 nats0 = group0.count();
    
        //get the table that holds the fields acting on all of the
        //atoms of this CutGroup (tables are indexed by CGIdx)
        BOOST_ASSERT(fields0_array[field0_idx].count() == int(nats0));
    
        Vector *group_fields0_array = fields0.data(field0_idx);

        //ok, we are interested in the fields acting on this CutGroup
        // - calculate all of the fieldss on this group interacting
        //   with all of the CutGroups in mol1 
        for (quint32 jgroup=0; jgroup<ngroups1; ++jgroup)
        {
            const CoordGroup &group1 = groups1_array[jgroup];
            const Parameters::Array &params1 = molparams1_array[jgroup];

            //check first that these two CoordGroups could be within cutoff
            //(if there is only one CutGroup in both molecules then this
            //test has already been performed and passed)
            const bool within_cutoff = (ngroups0 == 1 and ngroups1 == 1) or not
                                        spce->beyond(switchfunc->cutoffDistance(), 
                                                     aabox0, group1.aaBox());
            
            if (not within_cutoff)
                //this CutGroup is beyond the cutoff distance
                continue;
            
            //calculate all of the interatomic distances
            const double mindist = spce->calcDistVectors(group0, group1, distmat);
            
            if (mindist > switchfunc->cutoffDistance())
                //all of the atoms are definitely beyond cutoff
                continue;

            const quint32 nats1 = group1.count();
            
            //loop over all interatomic pairs and calculate the energies
            const Parameter *params1_array = params1.constData();

            if (mindist > switchfunc->featherDistance())
            {
                //we need to calculate the fields taking into account
                //the derivative of the switching function!
                
                //calculate the switching scale factors and their 
                //derivatives
                const double scl_coul = switchfunc->electrostaticScaleFactor(
                                                                    Length(mindist) );
                
                Vector group_sep = (group1.aaBox().center() - aabox0.center())
                                        .normalise();
                
                Vector dscl_coul = switchfunc->dElectrostaticScaleFactor(
                                                                    Length(mindist) ) 
                                     * group_sep;
                
                double shift_coul = 0;

                if (use_electrostatic_shifting)
                    shift_coul = probe.charge() * this->totalCharge(params1)
                                    / switchfunc->electrostaticCutoffDistance();
                
                for (quint32 i=0; i<nats0; ++i)
                {
                    distmat.setOuterIndex(i);
                
                    Vector total_field;
                
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        const double q2 = probe.reducedCharge() *
                                          params1_array[j].reduced_charge;
                          
                        if (q2 != 0)
                        {
                            //calculate the coulomb energy
                            const double cnrg = q2 / distmat[j].length();
                                               
                            //calculate the coulomb force
                            Vector cfield = (scl_coul * -cnrg / distmat[j].length() *
                                             distmat[j].direction()) +
                                             
                                            ((cnrg-shift_coul) * dscl_coul);
                        
                            total_field += cfield;
                        }
                    }
                    
                    //update the fields array
                    group_fields0_array[i] += scale_field * total_field;
                }
            }
            else
            {
                //not in the feather region, so can calculate the fields
                //directly (also, no need to calculate shift, as 
                //the shifting function is constant, so does not
                //affect the gradient)
                for (quint32 i=0; i<nats0; ++i)
                {
                    distmat.setOuterIndex(i);

                    Vector total_field;
                
                    //null LJ parameter - only add on the coulomb energy
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        const double q2 = params1_array[j].reduced_charge;
                        
                        //calculate the coulomb force
                        if (q2 != 0)
                        {
                            total_field += (-q2 / distmat[j].length2()) *
                                                  distmat[j].direction();
                        }
                    }
                    
                    group_fields0_array[i] += (scale_field * probe.reducedCharge()) *
                                                            total_field;

                } // end of loop over i atoms

            } // end of if within feather

        } // end of loop over jgroup CutGroups

    } // end of loop over igroup CutGroups
}

/** Add to the fields in 'fields0' the fields acting on the probe
    at all of the atoms sites of 'mol0' caused by 'mol1' */
void InterCLJPotential::_pvt_calculateLJField(
                                            const InterCLJPotential::Molecule &mol0, 
                                            const InterCLJPotential::Molecule &mol1,
                                            const CLJProbe &probe,
                                            MolFieldTable &fields0, 
                                            InterCLJPotential::FieldWorkspace &distmat,
                                            double scale_field) const
{
    BOOST_ASSERT( mol0.molecule().data().info().nCutGroups() == fields0.nCutGroups() );
    BOOST_ASSERT( mol0.molecule().data().number() == fields0.molNum() );

    if (probe.lj().isDummy())
        return;

    const quint32 ngroups0 = mol0.nCutGroups();
    const quint32 ngroups1 = mol1.nCutGroups();
    
    const CoordGroup *groups0_array = mol0.coordinates().constData();
    const CoordGroup *groups1_array = mol1.coordinates().constData();
    
    BOOST_ASSERT(mol1.parameters().atomicParameters().count() == int(ngroups1));
    
    const Parameters::Array *molparams1_array
                                    = mol1.parameters().atomicParameters().constData();
    
    const MolFieldTable::Array *fields0_array = fields0.constData();
    
    //loop over all pairs of CutGroups in the two molecules
    for (quint32 igroup=0; igroup<ngroups0; ++igroup)
    {
        //get the CGIdx of this group
        CGIdx cgidx_igroup = mol0.cgIdx(igroup);

        //get the index of this CutGroup in the forces array
        int field0_idx = fields0.map(cgidx_igroup);
        
        if (field0_idx == -1)
            //there is no space for the fields on this CutGroup in 
            //the fieldtable - were are therefore not interested in
            //this CutGroup
            continue;

        const CoordGroup &group0 = groups0_array[igroup];
        const AABox &aabox0 = group0.aaBox();
        const quint32 nats0 = group0.count();
    
        //get the table that holds the fields acting on all of the
        //atoms of this CutGroup (tables are indexed by CGIdx)
        BOOST_ASSERT(fields0_array[field0_idx].count() == int(nats0));
    
        Vector *group_fields0_array = fields0.data(field0_idx);

        //ok, we are interested in the fields acting on this CutGroup
        // - calculate all of the fields on this group interacting
        //   with all of the CutGroups in mol1 
        for (quint32 jgroup=0; jgroup<ngroups1; ++jgroup)
        {
            const CoordGroup &group1 = groups1_array[jgroup];
            const Parameters::Array &params1 = molparams1_array[jgroup];

            //check first that these two CoordGroups could be within cutoff
            //(if there is only one CutGroup in both molecules then this
            //test has already been performed and passed)
            const bool within_cutoff = (ngroups0 == 1 and ngroups1 == 1) or not
                                        spce->beyond(switchfunc->cutoffDistance(), 
                                                     aabox0, group1.aaBox());
            
            if (not within_cutoff)
                //this CutGroup is beyond the cutoff distance
                continue;
            
            //calculate all of the interatomic distances
            const double mindist = spce->calcDistVectors(group0, group1, distmat);
            
            if (mindist > switchfunc->cutoffDistance())
                //all of the atoms are definitely beyond cutoff
                continue;

            const quint32 nats1 = group1.count();
            
            //loop over all interatomic pairs and calculate the energies
            const Parameter *params1_array = params1.constData();

            if (mindist > switchfunc->featherDistance())
            {
                //we need to calculate the forces taking into account
                //the derivative of the switching function!
                
                //calculate the switching scale factors and their 
                //derivatives
                const double scl_lj = switchfunc->vdwScaleFactor( Length(mindist) );
                
                Vector group_sep = (group1.aaBox().center() - aabox0.center())
                                        .normalise();
                
                Vector dscl_lj = switchfunc->dVDWScaleFactor( Length(mindist) )
                                     * group_sep;
                
                for (quint32 i=0; i<nats0; ++i)
                {
                    distmat.setOuterIndex(i);

                    Vector total_field;
                
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        const Parameter &param1 = params1_array[j];

                        if (param1.ljid != 0)
                        {
                            const double invdist = double(1) / distmat[j].length();

                            LJPair ljpair( ljpairs.constData()[
                                                   ljpairs.map(param1.ljid,
                                                               param1.ljid)],
                                           probe.lj(), combining_rules );
                        
                            double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                            double sig_over_dist12 = pow_2(sig_over_dist6);

                            //calculate the energy
                            const double ljnrg = 4 * ljpair.epsilon() *
                                                  (sig_over_dist12 - sig_over_dist6);

                            // dU/dr requires an extra power of r
                            sig_over_dist6 *= invdist;
                            sig_over_dist12 *= invdist;

                            Vector field = ((scl_lj * 4 * ljpair.epsilon() * 
                                            (6.0*sig_over_dist6 - 12.0*sig_over_dist12))
                                               * distmat[j].direction())
                                            
                                             + (ljnrg * dscl_lj);

                            total_field += field;
                        }
                    }
                    
                    //update the fields array
                    group_fields0_array[i] += scale_field * total_field;
                }
            }
            else
            {
                //not in the feather region, so can calculate the fields
                //directly
                for (quint32 i=0; i<nats0; ++i)
                {
                    distmat.setOuterIndex(i);

                    Vector total_field;
                
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        const Parameter &param1 = params1_array[j];
                              
                        if (param1.ljid != 0)
                        {
                            const double invdist = double(1) / distmat[j].length();

                            LJPair ljpair( ljpairs.constData()[
                                                     ljpairs.map(param1.ljid,
                                                                 param1.ljid)],
                                           probe.lj(), combining_rules );
                        
                            double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                            double sig_over_dist12 = pow_2(sig_over_dist6);

                            // dU/dr requires an extra power of r
                            sig_over_dist6 *= invdist;
                            sig_over_dist12 *= invdist;

                            Vector field = (4 * ljpair.epsilon() * 
                                                       (6.0*sig_over_dist6 - 
                                                       12.0*sig_over_dist12))
                                        * distmat[j].direction();

                            total_field += field;
                        }
                    }
                    
                    group_fields0_array[i] += scale_field * total_field;

                } // end of loop over i atoms

            } // end of if within feather

        } // end of loop over jgroup CutGroups

    } // end of loop over igroup CutGroups
}

/////////////
///////////// Implementation of IntraCLJPotential
/////////////

static const RegisterMetaType<IntraCLJPotential> r_intraclj( MAGIC_ONLY, NO_ROOT,
                                            IntraCLJPotential::typeName() );

/** Serialise to a binary datastream */
QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds,
                                      const IntraCLJPotential &intraclj)
{
    writeHeader(ds, r_intraclj, 1);
    
    ds << static_cast<const CLJPotential&>(intraclj);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds,
                                      IntraCLJPotential &intraclj)
{
    VersionID v = readHeader(ds, r_intraclj);
    
    if (v == 1)
    {
        ds >> static_cast<CLJPotential&>(intraclj);
    }
    else
        throw version_error(v, "1", r_intraclj, CODELOC);
        
    return ds;
}

/** Constructor */
IntraCLJPotential::IntraCLJPotential() : CLJPotential()
{}

/** Copy constructor */
IntraCLJPotential::IntraCLJPotential(const IntraCLJPotential &other)
                  : CLJPotential(other)
{}

/** Destructor */
IntraCLJPotential::~IntraCLJPotential()
{}

/** Copy assignment operator */
IntraCLJPotential& IntraCLJPotential::operator=(const IntraCLJPotential &other)
{
    CLJPotential::operator=(other);
    return *this;
}

void IntraCLJPotential::throwMissingForceComponent(const Symbol &symbol,
                              const IntraCLJPotential::Components &components) const
{
    throw SireFF::missing_component( QObject::tr(
        "There is no force component in potential %1 - available "
        "components are %1, %2 and %3.")
            .arg(this->what())
            .arg(components.total().toString(), components.coulomb().toString(),
                 components.lj().toString()), CODELOC );
}

/** Assert that 'rest_of_mol' is compatible with 'mol'. They are only 
    compatible if they are both part of the same molecule (not necessarily
    the same version) with the same layout UID.
    
    \throw SireError::incompatible_error
*/
void IntraCLJPotential::assertCompatible(const IntraCLJPotential::Molecule &mol,
                                const IntraCLJPotential::Molecule &rest_of_mol) const
{
    if (mol.number() != rest_of_mol.number() or
        mol.molecule().data().info().UID() 
                        != rest_of_mol.molecule().data().info().UID())
    {
        throw SireError::incompatible_error( QObject::tr(
            "Problem adding the molecule %1 (%2). It doesn't appear to be "
            "in the same molecule as %3 (%4), or the molecule layout "
            "IDs may be different!")
                .arg(mol.molecule().name()).arg(mol.number())
                .arg(rest_of_mol.molecule().name()).arg(rest_of_mol.number()),
                    CODELOC );
    }
    else if (mol.parameters().intraScaleFactors() != 
                rest_of_mol.parameters().intraScaleFactors())
    {
        throw SireError::incompatible_error( QObject::tr(
            "The non-bonded scaling factors for the molecule %1 (%2) "
            "are different to the ones for the rest of the molecule. "
            "This is probably caused by a program bug...")
                .arg(mol.molecule().name()).arg(mol.number()),
                    CODELOC );
    }
}

/** Return all of the parameters needed by this potential for 
    the molecule 'molecule', using the supplied property map to
    find the properties that contain those parameters
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
IntraCLJPotential::Parameters 
IntraCLJPotential::getParameters(const PartialMolecule &molecule,
                                 const PropertyMap &map)
{
    need_update_ljpairs = true;

    return Parameters( AtomicParameters3D<CLJParameter>(
                               molecule, map[parameters().coordinates()],
                               getCLJParameters(molecule, 
                                                map[parameters().charge()],
                                                map[parameters().lj()]) ),
              IntraScaledParameters<CLJNBPairs>(
                               molecule, map[parameters().intraScaleFactors()] )
                     );
}

/** Update the parameters for the molecule going from 'old_molecule' to 
    'new_molecule', with the parameters found using the property map 'map'
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
IntraCLJPotential::Parameters
IntraCLJPotential::updateParameters(const IntraCLJPotential::Parameters &old_params,
                                    const PartialMolecule &old_molecule,
                                    const PartialMolecule &new_molecule,
                                    const PropertyMap &map)
{
    if (old_molecule.selection() != new_molecule.selection())
        //the selection has changed - just get completely new parameters
        return this->getParameters(new_molecule, map);

    Parameters new_params = old_params;

    //get the property names
    const PropertyName &coords_property = map[parameters().coordinates()];
    const PropertyName &chg_property = map[parameters().charge()];
    const PropertyName &lj_property = map[parameters().lj()];
    const PropertyName &scl_property = map[parameters().intraScaleFactors()];
    
    //get what has changed
    bool new_coords = old_molecule.version(coords_property) !=
                         new_molecule.version(coords_property);
                             
    bool new_clj = ( old_molecule.version(chg_property) !=
                         new_molecule.version(chg_property) ) or
                   ( old_molecule.version(lj_property) !=
                         new_molecule.version(lj_property) );

    bool new_scl = ( old_molecule.version(scl_property) !=
                         new_molecule.version(scl_property) );

    if (new_coords)
        new_params.setAtomicCoordinates( AtomicCoords3D(new_molecule, 
                                                        coords_property) );

    if (new_clj)
    {
        new_params.setAtomicParameters( getCLJParameters(new_molecule,
                                            chg_property, lj_property) );
        
        need_update_ljpairs = true;
    }

    if (new_scl)
        new_params.setIntraScaleFactors( 
                IntraScaledParameters<CLJNBPairs>(new_molecule, scl_property) );

    return new_params;
}
                 
/** Update the parameters for the molecule going from 'old_molecule' to 
    'new_molecule', also while the parameters of 'old_molecule'
    where found in 'old_map', now get the parameters using 'new_map'
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
IntraCLJPotential::Parameters
IntraCLJPotential::updateParameters(const IntraCLJPotential::Parameters &old_params,
                                    const PartialMolecule &old_molecule,
                                    const PartialMolecule &new_molecule,
                                    const PropertyMap &old_map, 
                                    const PropertyMap &new_map)
{
    if (old_molecule.selection() != new_molecule.selection())
        //the selection has changed - just get completely new parameters
        return this->getParameters(new_molecule, new_map);

    Parameters new_params = old_params;

    //get the property names
    const PropertyName &old_coords = old_map[parameters().coordinates()];
    const PropertyName &old_chg = old_map[parameters().charge()];
    const PropertyName &old_lj = old_map[parameters().lj()];
    const PropertyName &old_scl = old_map[parameters().intraScaleFactors()];
    
    const PropertyName &new_coords = new_map[parameters().coordinates()];
    const PropertyName &new_chg = new_map[parameters().charge()];
    const PropertyName &new_lj = new_map[parameters().lj()];
    const PropertyName &new_scl = new_map[parameters().intraScaleFactors()];
    
    //get what has changed
    bool changed_coords = (new_coords != old_coords) or
                           old_molecule.version(old_coords) !=
                           new_molecule.version(old_coords);
                             
    bool changed_clj = (new_chg != old_chg or new_lj != old_lj) or
                       ( old_molecule.version(old_chg) !=
                         new_molecule.version(old_chg) ) or
                       ( old_molecule.version(old_lj) !=
                         new_molecule.version(old_lj) );

    bool changed_scl = (new_scl != old_scl) or
                        old_molecule.version(old_scl) !=
                        new_molecule.version(old_scl);

    if (changed_coords)
        new_params.setAtomicCoordinates( AtomicCoords3D(new_molecule, new_coords) );

    if (changed_clj)
    {
        new_params.setAtomicParameters( getCLJParameters(new_molecule,
                                                         new_chg, new_lj) );
                                                         
        need_update_ljpairs = true;
    }

    if (changed_scl)
        new_params.setIntraScaleFactors( 
                        IntraScaledParameters<CLJNBPairs>(new_molecule, new_scl) );

    return new_params;
}

/** Return the IntraCLJPotential::Molecule representation of 'molecule',
    using the supplied PropertyMap to find the properties that contain
    the necessary forcefield parameters
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
IntraCLJPotential::Molecule
IntraCLJPotential::parameterise(const PartialMolecule &molecule,
                                const PropertyMap &map)
{
    return IntraCLJPotential::Molecule(molecule, *this, map);
}

/** Concert the passed group of molecules into IntraCLJPotential::Molecules,
    using the supplied PropertyMap to find the properties that contain
    the necessary forcefield parameters in each molecule
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
IntraCLJPotential::Molecules 
IntraCLJPotential::parameterise(const MoleculeGroup &molecules,
                                const PropertyMap &map)
{
    return IntraCLJPotential::Molecules(molecules, *this, map);
}

/** Return the total charge of the parameters for the group in 'params' */
double IntraCLJPotential::totalCharge(
                        const IntraCLJPotential::Parameters::Array &params) const
{
    int n = params.count();
    const Parameter *params_array = params.constData();
    
    double chg = 0;
    
    for (int i=0; i<n; ++i)
    {
        chg += params_array[i].reduced_charge;
    }
    
    return chg;
}

void IntraCLJPotential::calculateEnergy(const CLJNBPairs::CGPairs &group_pairs, 
                            IntraCLJPotential::EnergyWorkspace &distmat,
                            const IntraCLJPotential::Parameter *params0_array, 
                            const IntraCLJPotential::Parameter *params1_array,
                            const quint32 nats0, const quint32 nats1, 
                            double &cnrg, double &ljnrg) const
{
    double icnrg = 0;
    double iljnrg = 0;

    const double Rcoul = qMax(1e-5,qMin(1e9,
                            switchfunc->electrostaticCutoffDistance().to(angstrom)));
    const double Rlj = qMax(1e-5,qMin(1e9, switchfunc->vdwCutoffDistance().to(angstrom)) );
    const double Rc = qMax(Rcoul,Rlj);
        
    if (group_pairs.isEmpty())
    {
        //there is a constant scale factor between groups
        CLJScaleFactor cljscl = group_pairs.defaultValue();

        if (cljscl.coulomb() == 0 and cljscl.lj() == 0)
            return;

        if (use_electrostatic_shifting)
        {
            const double one_over_Rcoul = double(1) / Rcoul;
            const double one_over_Rcoul2 = double(1) / (Rcoul*Rcoul);
        
            for (quint32 i=0; i<nats0; ++i)
            {
                distmat.setOuterIndex(i);
                const Parameter &param0 = params0_array[i];
                    
                if (param0.ljid == 0)
                {
                    //null LJ parameter - only add on the coulomb energy
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        const Parameter &param1 = params1_array[j];

                        const double r = distmat[j];
                        
                        if (r < Rcoul)
                        {
                            const double one_over_r = double(1) / r;
                        
                            //if (cljscl.coulomb() != 1)
                            //    icnrg += cljscl.coulomb() *
                            //                param0.reduced_charge * param1.reduced_charge *
                            //                    one_over_r;
                            //else
                                icnrg += cljscl.coulomb() *
                                         param0.reduced_charge * param1.reduced_charge *
                                        (one_over_r - one_over_Rcoul + one_over_Rcoul2*(r-Rcoul));
                        }
                    }
                }
                else
                {
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        //do both coulomb and LJ
                        const Parameter &param1 = params1_array[j];
                            
                        const double r = distmat[j];
                        
                        if (r < Rc)
                        {
                            const double one_over_r = double(1) / r;

                            if (r < Rcoul)
                            {
                                //if (cljscl.coulomb() != 1)
                                //    icnrg += cljscl.coulomb() *
                                //            param0.reduced_charge * param1.reduced_charge *
                                //                one_over_r;
                                //else
                                    icnrg += cljscl.coulomb() *
                                             param0.reduced_charge * param1.reduced_charge *
                                        (one_over_r - one_over_Rcoul + one_over_Rcoul2*(r-Rcoul));
                            }

                            if (param1.ljid != 0 and r < Rlj)
                            {
                                const LJPair &ljpair = ljpairs.constData()[
                                                          ljpairs.map(param0.ljid,
                                                                      param1.ljid)];
                                
                                double sig_over_dist6 = pow_6(ljpair.sigma()*one_over_r);
                                double sig_over_dist12 = pow_2(sig_over_dist6);

                                iljnrg += cljscl.lj() *
                                              ljpair.epsilon() * (sig_over_dist12 - 
                                                                  sig_over_dist6);
                            }
                        }
                    }
                }
            }
        }
        else if (use_reaction_field)
        {
            const double k_rf = (1.0 / pow_3(Rcoul)) * ( (rf_dielectric_constant-1) /
                                                         (2*rf_dielectric_constant + 1) );
            const double c_rf = (1.0 / Rcoul) * ( (3*rf_dielectric_constant) /
                                                  (2*rf_dielectric_constant + 1) );
        
            for (quint32 i=0; i<nats0; ++i)
            {
                distmat.setOuterIndex(i);
                const Parameter &param0 = params0_array[i];
                    
                if (param0.ljid == 0)
                {
                    //null LJ parameter - only add on the coulomb energy
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        const Parameter &param1 = params1_array[j];

                        const double r = distmat[j];
                        
                        if (r < Rcoul)
                        {
                            const double one_over_r = double(1) / r;
                        
                            //if (cljscl.coulomb() != 1)
                            //    icnrg += cljscl.coulomb() *
                            //                param0.reduced_charge * param1.reduced_charge *
                            //                    one_over_r;
                            //else
                                icnrg += cljscl.coulomb() *
                                         param0.reduced_charge * param1.reduced_charge *
                                            (one_over_r + k_rf*r*r - c_rf);
                        }
                    }
                }
                else
                {
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        //do both coulomb and LJ
                        const Parameter &param1 = params1_array[j];
                            
                        const double r = distmat[j];
                        
                        if (r < Rc)
                        {
                            const double one_over_r = double(1) / r;
                            
                            if (r < Rcoul)
                            {
                                //if (cljscl.coulomb() != 1)
                                //    icnrg += cljscl.coulomb() *
                                //                param0.reduced_charge * param1.reduced_charge *
                                //                    one_over_r;
                                //else
                                    icnrg += cljscl.coulomb() *
                                              param0.reduced_charge * param1.reduced_charge *
                                                    (one_over_r + k_rf*r*r - c_rf);
                            }
                                  
                            if (param1.ljid != 0 and r < Rlj)
                            {
                                const LJPair &ljpair = ljpairs.constData()[
                                                          ljpairs.map(param0.ljid,
                                                                      param1.ljid)];
                                
                                double sig_over_dist6 = pow_6(ljpair.sigma()*one_over_r);
                                double sig_over_dist12 = pow_2(sig_over_dist6);

                                iljnrg += cljscl.lj() *
                                              ljpair.epsilon() * (sig_over_dist12 - 
                                                                  sig_over_dist6);
                            }
                        }
                    }
                }
            }
        }
        else if (use_atomistic_cutoff)
        {
            for (quint32 i=0; i<nats0; ++i)
            {
                distmat.setOuterIndex(i);
                const Parameter &param0 = params0_array[i];
                    
                if (param0.ljid == 0)
                {
                    //null LJ parameter - only add on the coulomb energy
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        const Parameter &param1 = params1_array[j];

                        const double r = distmat[j];
                        
                        if (r < Rcoul)
                        {
                            const double one_over_r = double(1) / r;
                        
                            icnrg += cljscl.coulomb() *
                                        param0.reduced_charge * param1.reduced_charge *
                                            one_over_r;
                        }
                    }
                }
                else
                {
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        //do both coulomb and LJ
                        const Parameter &param1 = params1_array[j];
                            
                        const double r = distmat[j];
                        
                        if (r < Rc)
                        {
                            const double one_over_r = double(1) / r;

                            if (r < Rcoul)
                            {
                                icnrg += cljscl.coulomb() *
                                            param0.reduced_charge * param1.reduced_charge *
                                                one_over_r;
                            }
                            
                            if (param1.ljid != 0 and r < Rlj)
                            {
                                const LJPair &ljpair = ljpairs.constData()[
                                                          ljpairs.map(param0.ljid,
                                                                      param1.ljid)];
                                
                                double sig_over_dist6 = pow_6(ljpair.sigma()*one_over_r);
                                double sig_over_dist12 = pow_2(sig_over_dist6);

                                iljnrg += cljscl.lj() *
                                              ljpair.epsilon() * (sig_over_dist12 - 
                                                                  sig_over_dist6);
                            }
                        }
                    }
                }
            }
        }
        else
        {
            //group-based feathered cutoff
            for (quint32 i=0; i<nats0; ++i)
            {
                distmat.setOuterIndex(i);
                const Parameter &param0 = params0_array[i];
                    
                if (param0.ljid == 0)
                {
                    //null LJ parameter - only add on the coulomb energy
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        icnrg += cljscl.coulomb() *
                                 param0.reduced_charge * 
                                 params1_array[j].reduced_charge / distmat[j];
                    }
                }
                else
                {
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        //do both coulomb and LJ
                        const Parameter &param1 = params1_array[j];
                            
                        const double invdist = double(1) / distmat[j];
                            
                        icnrg += cljscl.coulomb() *
                                 param0.reduced_charge * param1.reduced_charge
                                                       * invdist;
                                  
                        if (param1.ljid != 0)
                        {
                            const LJPair &ljpair = ljpairs.constData()[
                                                      ljpairs.map(param0.ljid,
                                                                  param1.ljid)];
                            
                            double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                            double sig_over_dist12 = pow_2(sig_over_dist6);

                            iljnrg += cljscl.lj() *
                                          ljpair.epsilon() * (sig_over_dist12 - 
                                                              sig_over_dist6);
                        }
                    }
                }
            }
        } // end of if use_electrostatic_shifting
    }
    else
    {
        //there are different nb scale factors between
        //the atoms. We need to calculate the energies using
        //them...
        if (use_electrostatic_shifting)
        {
            const double one_over_Rcoul = double(1) / Rcoul;
            const double one_over_Rcoul2 = double(1) / (Rcoul*Rcoul);
        
            for (quint32 i=0; i<nats0; ++i)
            {
                distmat.setOuterIndex(i);
                const Parameter &param0 = params0_array[i];
                    
                if (param0.ljid == 0)
                {
                    //null LJ parameter - only add on the coulomb energy
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        const CLJScaleFactor &cljscl = group_pairs(i,j);
                                
                        if (cljscl.coulomb() != 0)
                        {
                            const Parameter &param1 = params1_array[j];

                            const double r = distmat[j];
                            
                            if (r < Rcoul)
                            {
                                const double one_over_r = double(1) / r;
                            
                                if (cljscl.coulomb() != 1)
                                    icnrg += cljscl.coulomb() *
                                                param0.reduced_charge * param1.reduced_charge *
                                                    one_over_r;
                                else
                                    icnrg += param0.reduced_charge * param1.reduced_charge *
                                        (one_over_r - one_over_Rcoul + one_over_Rcoul2*(r-Rcoul));
                            }
                        }
                    }
                }
                else
                {
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        //do both coulomb and LJ
                        const CLJScaleFactor &cljscl = group_pairs(i,j);

                        if (cljscl.coulomb() != 0 or cljscl.lj() != 0)
                        {
                            const Parameter &param1 = params1_array[j];

                            const double r = distmat[j];
                            
                            if (r < Rc)
                            {
                                const double one_over_r = double(1) / r;
                            
                                //coulomb
                                if (r < Rcoul)
                                {
                                    if (cljscl.coulomb() != 1)
                                        icnrg += cljscl.coulomb() *
                                                    param0.reduced_charge * param1.reduced_charge *
                                                        one_over_r;
                                    else
                                        icnrg += param0.reduced_charge * param1.reduced_charge *
                                          (one_over_r - one_over_Rcoul + one_over_Rcoul2*(r-Rcoul));
                                }
                            
                                //lj
                                if (cljscl.lj() != 0 and param1.ljid != 0 and r < Rlj)
                                {
                                    const LJPair &ljpair = ljpairs.constData()[
                                                                ljpairs.map(param0.ljid,
                                                                            param1.ljid)];
                                
                                    double sig_over_dist6 = pow_6(ljpair.sigma()*one_over_r);
                                    double sig_over_dist12 = pow_2(sig_over_dist6);

                                    iljnrg += cljscl.lj() * ljpair.epsilon() *
                                            (sig_over_dist12 - sig_over_dist6);
                                }
                            }
                        }
                    }
                }
            }
        }
        else if (use_reaction_field)
        {
            const double k_rf = (1.0 / pow_3(Rcoul)) * ( (rf_dielectric_constant-1) /
                                                         (2*rf_dielectric_constant + 1) );
            const double c_rf = (1.0 / Rcoul) * ( (3*rf_dielectric_constant) /
                                                  (2*rf_dielectric_constant + 1) );
        
            for (quint32 i=0; i<nats0; ++i)
            {
                distmat.setOuterIndex(i);
                const Parameter &param0 = params0_array[i];
                    
                if (param0.ljid == 0)
                {
                    //null LJ parameter - only add on the coulomb energy
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        const CLJScaleFactor &cljscl = group_pairs(i,j);
                                
                        if (cljscl.coulomb() != 0)
                        {
                            const Parameter &param1 = params1_array[j];

                            const double r = distmat[j];
                            
                            if (r < Rcoul)
                            {
                                const double one_over_r = double(1) / r;
                            
                                if (cljscl.coulomb() != 1)
                                    icnrg += cljscl.coulomb() *
                                                param0.reduced_charge * param1.reduced_charge *
                                                    one_over_r;
                                else
                                    icnrg += param0.reduced_charge * param1.reduced_charge *
                                                (one_over_r + k_rf*r*r - c_rf);
                            }
                        }
                    }
                }
                else
                {
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        //do both coulomb and LJ
                        const CLJScaleFactor &cljscl = group_pairs(i,j);

                        if (cljscl.coulomb() != 0 or cljscl.lj() != 0)
                        {
                            const Parameter &param1 = params1_array[j];

                            const double r = distmat[j];

                            if (r < Rc)
                            {
                                const double one_over_r = double(1) / r;
                            
                                if (r < Rcoul)
                                {
                                    if (cljscl.coulomb() != 1)
                                        icnrg += cljscl.coulomb() *
                                                    param0.reduced_charge * param1.reduced_charge *
                                                        one_over_r;
                                    else
                                        icnrg += param0.reduced_charge * param1.reduced_charge *
                                                    (one_over_r + k_rf*r*r - c_rf);
                                }

                                if (cljscl.lj() != 0 and param1.ljid != 0 and r < Rlj)
                                {
                                    const LJPair &ljpair = ljpairs.constData()[
                                                            ljpairs.map(param0.ljid,
                                                                        param1.ljid)];
                                
                                    double sig_over_dist6 = pow_6(ljpair.sigma()*one_over_r);
                                    double sig_over_dist12 = pow_2(sig_over_dist6);

                                    iljnrg += cljscl.lj() * ljpair.epsilon() *
                                            (sig_over_dist12 - sig_over_dist6);
                                }
                            }
                        }
                    }
                }
            }
        }
        else if (use_atomistic_cutoff)
        {
            for (quint32 i=0; i<nats0; ++i)
            {
                distmat.setOuterIndex(i);
                const Parameter &param0 = params0_array[i];
                    
                if (param0.ljid == 0)
                {
                    //null LJ parameter - only add on the coulomb energy
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        const CLJScaleFactor &cljscl = group_pairs(i,j);
                                
                        if (cljscl.coulomb() != 0)
                        {
                            const Parameter &param1 = params1_array[j];

                            const double r = distmat[j];
                            
                            if (r < Rcoul)
                            {
                                const double one_over_r = double(1) / r;
                            
                                icnrg += cljscl.coulomb() *
                                            param0.reduced_charge * param1.reduced_charge *
                                                one_over_r;
                            }
                        }
                    }
                }
                else
                {
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        //do both coulomb and LJ
                        const CLJScaleFactor &cljscl = group_pairs(i,j);

                        if (cljscl.coulomb() != 0 or cljscl.lj() != 0)
                        {
                            const Parameter &param1 = params1_array[j];

                            const double r = distmat[j];
                            
                            if (r < Rc)
                            {
                                const double one_over_r = double(1) / r;
                            
                                if (r < Rcoul)
                                    icnrg += cljscl.coulomb() *
                                                param0.reduced_charge * param1.reduced_charge *
                                                    one_over_r;

                                if (cljscl.lj() != 0 and param1.ljid != 0 and r < Rlj)
                                {
                                    const LJPair &ljpair = ljpairs.constData()[
                                                             ljpairs.map(param0.ljid,
                                                                         param1.ljid)];
                                
                                    double sig_over_dist6 = pow_6(ljpair.sigma()*one_over_r);
                                    double sig_over_dist12 = pow_2(sig_over_dist6);

                                    iljnrg += cljscl.lj() * ljpair.epsilon() * 
                                               (sig_over_dist12 - sig_over_dist6);
                                }
                            }
                        }
                    }
                }
            }
        }
        else // if use_electrostatic_shifting
        {
            //using the group-based feathered cutoff
        
            for (quint32 i=0; i<nats0; ++i)
            {
                distmat.setOuterIndex(i);
                const Parameter &param0 = params0_array[i];
                    
                if (param0.ljid == 0)
                {
                    //null LJ parameter - only add on the coulomb energy
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        const CLJScaleFactor &cljscl = group_pairs(i,j);
                                
                        if (cljscl.coulomb() != 0)
                        {
                            icnrg += cljscl.coulomb() * 
                                        param0.reduced_charge * 
                                        params1_array[j].reduced_charge / distmat[j];
                        }
                    }
                }
                else
                {
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        //do both coulomb and LJ
                        const CLJScaleFactor &cljscl = group_pairs(i,j);

                        if (cljscl.coulomb() != 0 or cljscl.lj() != 0)
                        {
                            const Parameter &param1 = params1_array[j];
                            
                            const double invdist = double(1) / distmat[j];
                            
                            icnrg += cljscl.coulomb() *  
                                     param0.reduced_charge * 
                                     param1.reduced_charge * invdist;

                            if (cljscl.lj() != 0 and param1.ljid != 0)
                            {
                                const LJPair &ljpair = ljpairs.constData()[
                                                         ljpairs.map(param0.ljid,
                                                                     param1.ljid)];
                            
                                double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                                double sig_over_dist12 = pow_2(sig_over_dist6);

                                iljnrg += cljscl.lj() * ljpair.epsilon() * 
                                           (sig_over_dist12 - sig_over_dist6);
                            }
                        }
                    }
                }
            }
        } // end of if use_electrostatic_shifting
    }
    
    cnrg += icnrg;
    ljnrg += 4*iljnrg;
}

void IntraCLJPotential::calculateEnergy(const CLJNBPairs::CGPairs &group_pairs, 
                            const QSet<Index> &atoms0, const QSet<Index> &atoms1,
                            IntraCLJPotential::EnergyWorkspace &distmat,
                            const IntraCLJPotential::Parameter *params0_array, 
                            const IntraCLJPotential::Parameter *params1_array,
                            const quint32 nats0, const quint32 nats1, 
                            double &cnrg, double &ljnrg) const
{
    if (atoms0.isEmpty() or atoms1.isEmpty())
        return;

    double icnrg = 0;
    double iljnrg = 0;

    const double Rcoul = qMax(1e-5,qMin(1e9,
                            switchfunc->electrostaticCutoffDistance().to(angstrom)));
    const double Rlj = qMax(1e-5,qMin(1e9, switchfunc->vdwCutoffDistance().to(angstrom)) );
    const double Rc = qMax(Rcoul,Rlj);

    if (group_pairs.isEmpty())
    {
        //there is a constant scale factor between groups
        CLJScaleFactor cljscl = group_pairs.defaultValue();

        if (cljscl.coulomb() == 0 and cljscl.lj() == 0)
            return;

        if (use_electrostatic_shifting)
        {
            const double one_over_Rcoul = double(1) / Rcoul;
            const double one_over_Rcoul2 = double(1) / (Rcoul*Rcoul);
        
            foreach (Index i, atoms0)
            {
                distmat.setOuterIndex(i);
                const Parameter &param0 = params0_array[i];
                    
                if (param0.ljid == 0)
                {
                    //null LJ parameter - only add on the coulomb energy
                    foreach (Index j, atoms1)
                    {
                        const Parameter &param1 = params1_array[j];

                        const double r = distmat[j];
                        
                        if (r < Rcoul)
                        {
                            const double one_over_r = double(1) / r;
                        
                            if (cljscl.coulomb() != 1)
                                icnrg += cljscl.coulomb() *
                                            param0.reduced_charge * param1.reduced_charge *
                                                one_over_r;
                            else
                                icnrg += param0.reduced_charge * param1.reduced_charge *
                                    (one_over_r - one_over_Rcoul + one_over_Rcoul2*(r-Rcoul));
                        }
                    }
                }
                else
                {
                    foreach (Index j, atoms1)
                    {
                        //do both coulomb and LJ
                        const Parameter &param1 = params1_array[j];

                        const double r = distmat[j];
                        
                        if (r < Rc)
                        {
                            const double one_over_r = double(1) / r;
                        
                            if (r < Rcoul)
                            {
                                if (cljscl.coulomb() != 1)
                                    icnrg += cljscl.coulomb() *
                                                param0.reduced_charge * param1.reduced_charge *
                                                    one_over_r;
                                else
                                    icnrg += param0.reduced_charge * param1.reduced_charge *
                                        (one_over_r - one_over_Rcoul + one_over_Rcoul2*(r-Rcoul));
                            }
                            
                            if (param1.ljid != 0 and r < Rlj)
                            {
                                const LJPair &ljpair = ljpairs.constData()[
                                                          ljpairs.map(param0.ljid,
                                                                      param1.ljid)];
                                
                                double sig_over_dist6 = pow_6(ljpair.sigma()*one_over_r);
                                double sig_over_dist12 = pow_2(sig_over_dist6);

                                iljnrg += cljscl.lj() *
                                              ljpair.epsilon() * (sig_over_dist12 - 
                                                                  sig_over_dist6);
                            }
                        }
                    }
                }
            }
        }
        else if (use_reaction_field)
        {
            const double k_rf = (1.0 / pow_3(Rcoul)) * ( (rf_dielectric_constant-1) /
                                                         (2*rf_dielectric_constant + 1) );
            const double c_rf = (1.0 / Rcoul) * ( (3*rf_dielectric_constant) /
                                                  (2*rf_dielectric_constant + 1) );
        
            foreach (Index i, atoms0)
            {
                distmat.setOuterIndex(i);
                const Parameter &param0 = params0_array[i];
                    
                if (param0.ljid == 0)
                {
                    //null LJ parameter - only add on the coulomb energy
                    foreach (Index j, atoms1)
                    {
                        const Parameter &param1 = params1_array[j];

                        const double r = distmat[j];
                        
                        if (r < Rcoul)
                        {
                            const double one_over_r = double(1) / r;
                        
                            if (cljscl.coulomb() != 1)
                                icnrg += cljscl.coulomb() *
                                            param0.reduced_charge * param1.reduced_charge *
                                                  (one_over_r);
                            else
                                icnrg += param0.reduced_charge * param1.reduced_charge *
                                                  (one_over_r + k_rf*r*r - c_rf);
                        }
                    }
                }
                else
                {
                    foreach (Index j, atoms1)
                    {
                        //do both coulomb and LJ
                        const Parameter &param1 = params1_array[j];

                        const double r = distmat[j];
                        
                        if (r < Rc)
                        {
                            const double one_over_r = double(1) / r;
                        
                            if (r < Rcoul)
                            {
                                if (cljscl.coulomb() != 1)
                                    icnrg += cljscl.coulomb() *
                                                param0.reduced_charge * param1.reduced_charge *
                                                    one_over_r;
                                else
                                    icnrg += param0.reduced_charge * param1.reduced_charge *
                                                (one_over_r + k_rf*r*r - c_rf);
                            }
                            
                            if (param1.ljid != 0 and r < Rlj)
                            {
                                const LJPair &ljpair = ljpairs.constData()[
                                                          ljpairs.map(param0.ljid,
                                                                      param1.ljid)];
                                
                                double sig_over_dist6 = pow_6(ljpair.sigma()*one_over_r);
                                double sig_over_dist12 = pow_2(sig_over_dist6);

                                iljnrg += cljscl.lj() *
                                              ljpair.epsilon() * (sig_over_dist12 - 
                                                                  sig_over_dist6);
                            }
                        }
                    }
                }
            }
        }
        else if (use_atomistic_cutoff)
        {
            foreach (Index i, atoms0)
            {
                distmat.setOuterIndex(i);
                const Parameter &param0 = params0_array[i];
                    
                if (param0.ljid == 0)
                {
                    //null LJ parameter - only add on the coulomb energy
                    foreach (Index j, atoms1)
                    {
                        const Parameter &param1 = params1_array[j];

                        const double r = distmat[j];
                        
                        if (r < Rcoul)
                        {
                            const double one_over_r = double(1) / r;
                        
                            icnrg += cljscl.coulomb() *
                                        param0.reduced_charge * param1.reduced_charge *
                                              (one_over_r);
                        }
                    }
                }
                else
                {
                    foreach (Index j, atoms1)
                    {
                        //do both coulomb and LJ
                        const Parameter &param1 = params1_array[j];

                        const double r = distmat[j];
                        
                        if (r < Rc)
                        {
                            const double one_over_r = double(1) / r;
                        
                            if (r < Rcoul)
                            {
                                icnrg += cljscl.coulomb() *
                                            param0.reduced_charge * param1.reduced_charge *
                                                one_over_r;
                            }
                                  
                            if (param1.ljid != 0 and r < Rlj)
                            {
                                const LJPair &ljpair = ljpairs.constData()[
                                                          ljpairs.map(param0.ljid,
                                                                      param1.ljid)];
                                
                                double sig_over_dist6 = pow_6(ljpair.sigma()*one_over_r);
                                double sig_over_dist12 = pow_2(sig_over_dist6);

                                iljnrg += cljscl.lj() *
                                              ljpair.epsilon() * (sig_over_dist12 - 
                                                                  sig_over_dist6);
                            }
                        }
                    }
                }
            }
        }
        else
        {
            //group-based feathered cutoff
        
            foreach (Index i, atoms0)
            {
                distmat.setOuterIndex(i);
                const Parameter &param0 = params0_array[i];
                    
                if (param0.ljid == 0)
                {
                    //null LJ parameter - only add on the coulomb energy
                    foreach (Index j, atoms1)
                    {
                        icnrg += cljscl.coulomb() * 
                                 param0.reduced_charge * 
                                 params1_array[j].reduced_charge / distmat[j];
                    }
                }
                else
                {
                    foreach (Index j, atoms1)
                    {
                        //do both coulomb and LJ
                        const Parameter &param1 = params1_array[j];
                            
                        const double invdist = double(1) / distmat[j];
                            
                        icnrg += cljscl.coulomb() *
                                 param0.reduced_charge * param1.reduced_charge
                                                       * invdist;
                                  
                        if (param1.ljid != 0)
                        {
                            const LJPair &ljpair = ljpairs.constData()[
                                                      ljpairs.map(param0.ljid,
                                                                  param1.ljid)];
                            
                            double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                            double sig_over_dist12 = pow_2(sig_over_dist6);

                            iljnrg += cljscl.lj() *
                                          ljpair.epsilon() * (sig_over_dist12 - 
                                                              sig_over_dist6);
                        }
                    }
                }
            }
        } // end of if use_electrostatic_shifting
    }
    else
    {
        //there are different nb scale factors between
        //the atoms. We need to calculate the energies using
        //them...
        if (use_electrostatic_shifting)
        {
            const double one_over_Rcoul = double(1) / Rcoul;
            const double one_over_Rcoul2 = double(1) / (Rcoul*Rcoul);

            foreach (Index i, atoms0)
            {
                distmat.setOuterIndex(i);
                const Parameter &param0 = params0_array[i];
                    
                if (param0.ljid == 0)
                {
                    //null LJ parameter - only add on the coulomb energy
                    foreach (Index j, atoms1)
                    {
                        const CLJScaleFactor &cljscl = group_pairs(i,j);
                                
                        if (cljscl.coulomb() != 0)
                        {
                            const Parameter &param1 = params1_array[j];

                            const double r = distmat[j];
                            
                            if (r < Rcoul)
                            {
                                const double one_over_r = double(1) / r;
                            
                                if (cljscl.coulomb() != 1)
                                    icnrg += cljscl.coulomb() *
                                                param0.reduced_charge * param1.reduced_charge *
                                                    one_over_r;
                                else
                                    icnrg += param0.reduced_charge * param1.reduced_charge *
                                         (one_over_r - one_over_Rcoul + one_over_Rcoul2*(r-Rcoul));
                            }
                        }
                    }
                }
                else
                {
                    foreach (Index j, atoms1)
                    {
                        //do both coulomb and LJ
                        const CLJScaleFactor &cljscl = group_pairs(i,j);

                        if (cljscl.coulomb() != 0 or cljscl.lj() != 0)
                        {
                            const Parameter &param1 = params1_array[j];

                            const double r = distmat[j];
                            
                            if (r < Rc)
                            {
                                const double one_over_r = double(1) / r;
                            
                                if (cljscl.coulomb() != 0 and r < Rcoul)
                                {
                                    if (cljscl.coulomb() != 1)
                                        icnrg += cljscl.coulomb() *
                                                    param0.reduced_charge * param1.reduced_charge *
                                                        one_over_r;
                                    else
                                        icnrg += param0.reduced_charge * param1.reduced_charge *
                                         (one_over_r - one_over_Rcoul + one_over_Rcoul2*(r-Rcoul));
                                }

                                if (cljscl.lj() != 0 and param1.ljid != 0 and r < Rlj)
                                {
                                    const LJPair &ljpair = ljpairs.constData()[
                                                             ljpairs.map(param0.ljid,
                                                                         param1.ljid)];
                                
                                    double sig_over_dist6 = pow_6(ljpair.sigma()*one_over_r);
                                    double sig_over_dist12 = pow_2(sig_over_dist6);

                                    iljnrg += cljscl.lj() * ljpair.epsilon() * 
                                               (sig_over_dist12 - sig_over_dist6);
                                }
                            }
                        }
                    }
                }
            }
        }
        else if (use_reaction_field)
        {
            const double k_rf = (1.0 / pow_3(Rcoul)) * ( (rf_dielectric_constant-1) /
                                                         (2*rf_dielectric_constant + 1) );
            const double c_rf = (1.0 / Rcoul) * ( (3*rf_dielectric_constant) /
                                                  (2*rf_dielectric_constant + 1) );

            foreach (Index i, atoms0)
            {
                distmat.setOuterIndex(i);
                const Parameter &param0 = params0_array[i];
                    
                if (param0.ljid == 0)
                {
                    //null LJ parameter - only add on the coulomb energy
                    foreach (Index j, atoms1)
                    {
                        const CLJScaleFactor &cljscl = group_pairs(i,j);
                                
                        if (cljscl.coulomb() != 0)
                        {
                            const Parameter &param1 = params1_array[j];

                            const double r = distmat[j];
                            
                            if (r < Rc)
                            {
                                const double one_over_r = double(1) / r;
                            
                                if (cljscl.coulomb() != 1)
                                    icnrg += cljscl.coulomb() *
                                                param0.reduced_charge * param1.reduced_charge *
                                                    one_over_r;
                                else
                                    icnrg += param0.reduced_charge * param1.reduced_charge *
                                                (one_over_r + k_rf*r*r - c_rf);
                            }
                        }
                    }
                }
                else
                {
                    foreach (Index j, atoms1)
                    {
                        //do both coulomb and LJ
                        const CLJScaleFactor &cljscl = group_pairs(i,j);

                        if (cljscl.coulomb() != 0 or cljscl.lj() != 0)
                        {
                            const Parameter &param1 = params1_array[j];

                            const double r = distmat[j];
                            
                            if (r < Rc)
                            {
                                const double one_over_r = double(1) / r;
                            
                                if (cljscl.coulomb() != 0 and r < Rcoul)
                                {
                                    if (cljscl.coulomb() != 1)
                                        icnrg += cljscl.coulomb() *
                                                    param0.reduced_charge * param1.reduced_charge *
                                                        one_over_r;
                                    else
                                        icnrg += param0.reduced_charge * param1.reduced_charge *
                                                    (one_over_r + k_rf*r*r - c_rf);
                                }
                                
                                if (cljscl.lj() != 0 and param1.ljid != 0 and r < Rlj)
                                {
                                    const LJPair &ljpair = ljpairs.constData()[
                                                             ljpairs.map(param0.ljid,
                                                                         param1.ljid)];
                                
                                    double sig_over_dist6 = pow_6(ljpair.sigma()*one_over_r);
                                    double sig_over_dist12 = pow_2(sig_over_dist6);

                                    iljnrg += cljscl.lj() * ljpair.epsilon() * 
                                               (sig_over_dist12 - sig_over_dist6);
                                }
                            }
                        }
                    }
                }
            }
        }
        else if (use_atomistic_cutoff)
        {
            foreach (Index i, atoms0)
            {
                distmat.setOuterIndex(i);
                const Parameter &param0 = params0_array[i];
                    
                if (param0.ljid == 0)
                {
                    //null LJ parameter - only add on the coulomb energy
                    foreach (Index j, atoms1)
                    {
                        const CLJScaleFactor &cljscl = group_pairs(i,j);
                                
                        if (cljscl.coulomb() != 0)
                        {
                            const Parameter &param1 = params1_array[j];

                            const double r = distmat[j];
                            
                            if (r < Rcoul)
                            {
                                const double one_over_r = double(1) / r;
                            
                                icnrg += cljscl.coulomb() *
                                            param0.reduced_charge * param1.reduced_charge *
                                                one_over_r;
                            }
                        }
                    }
                }
                else
                {
                    foreach (Index j, atoms1)
                    {
                        //do both coulomb and LJ
                        const CLJScaleFactor &cljscl = group_pairs(i,j);

                        if (cljscl.coulomb() != 0 or cljscl.lj() != 0)
                        {
                            const Parameter &param1 = params1_array[j];

                            const double r = distmat[j];
                            
                            if (r < Rc)
                            {
                                const double one_over_r = double(1) / r;
                            
                                if (cljscl.coulomb() != 0 and r < Rcoul)
                                {
                                    icnrg += cljscl.coulomb() *
                                                param0.reduced_charge * param1.reduced_charge *
                                                    one_over_r;
                                }
                                      
                                if (cljscl.lj() != 0 and param1.ljid != 0 and r < Rlj)
                                {
                                    const LJPair &ljpair = ljpairs.constData()[
                                                             ljpairs.map(param0.ljid,
                                                                         param1.ljid)];
                                
                                    double sig_over_dist6 = pow_6(ljpair.sigma()*one_over_r);
                                    double sig_over_dist12 = pow_2(sig_over_dist6);

                                    iljnrg += cljscl.lj() * ljpair.epsilon() * 
                                               (sig_over_dist12 - sig_over_dist6);
                                }
                            }
                        }
                    }
                }
            }
        }
        else
        {
            //group-based feathered cutoff
        
            foreach (Index i, atoms0)
            {
                distmat.setOuterIndex(i);
                const Parameter &param0 = params0_array[i];
                    
                if (param0.ljid == 0)
                {
                    //null LJ parameter - only add on the coulomb energy
                    foreach (Index j, atoms1)
                    {
                        const CLJScaleFactor &cljscl = group_pairs(i,j);
                                
                        if (cljscl.coulomb() != 0)
                               icnrg += cljscl.coulomb() * 
                                        param0.reduced_charge * 
                                        params1_array[j].reduced_charge / distmat[j];
                    }
                }
                else
                {
                    foreach (Index j, atoms1)
                    {
                        //do both coulomb and LJ
                        const CLJScaleFactor &cljscl = group_pairs(i,j);

                        if (cljscl.coulomb() != 0 or cljscl.lj() != 0)
                        {
                            const Parameter &param1 = params1_array[j];
                            
                            const double invdist = double(1) / distmat[j];
                            
                            icnrg += cljscl.coulomb() *  
                                     param0.reduced_charge * 
                                     param1.reduced_charge * invdist;
                                  
                            if (cljscl.lj() != 0 and param1.ljid != 0)
                            {
                                const LJPair &ljpair = ljpairs.constData()[
                                                         ljpairs.map(param0.ljid,
                                                                     param1.ljid)];
                            
                                double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                                double sig_over_dist12 = pow_2(sig_over_dist6);

                                iljnrg += cljscl.lj() * ljpair.epsilon() * 
                                           (sig_over_dist12 - sig_over_dist6);
                            }
                        }
                    }
                }
            }
        } // end of if use_electrostatic_shifting
    }
    
    cnrg += icnrg;
    ljnrg += 4*iljnrg;
}

/** Calculate the intramolecular CLJ energy of the passed molecule, and
    add this onto 'energy'. This uses the passed workspace when
    performing the calculation */
void IntraCLJPotential::calculateEnergy(const IntraCLJPotential::Molecule &mol,
                                        IntraCLJPotential::Energy &energy,
                                        IntraCLJPotential::EnergyWorkspace &distmat,
                                        double scale_energy) const
{
    if (scale_energy == 0 or mol.isEmpty())
        return;

    const quint32 ngroups = mol.nCutGroups();
    
    const CoordGroup *groups_array = mol.coordinates().constData();
    
    BOOST_ASSERT( mol.parameters().atomicParameters().count() == int(ngroups) );
    const Parameters::Array *molparams_array 
                            = mol.parameters().atomicParameters().constData();

    const CLJNBPairs &nbpairs = mol.parameters().intraScaleFactors();
    
    double cnrg = 0;
    double ljnrg = 0;
      
    //loop over all pairs of CutGroups in the molecule
    for (quint32 igroup=0; igroup<ngroups; ++igroup)
    {
        const Parameters::Array &params0 = molparams_array[igroup];

        const CoordGroup &group0 = groups_array[igroup];
        const AABox &aabox0 = group0.aaBox();
        const quint32 nats0 = group0.count();
        const Parameter *params0_array = params0.constData();
    
        CGIdx cgidx_igroup = mol.cgIdx(igroup);
    
        for (quint32 jgroup=igroup; jgroup<ngroups; ++jgroup)
        {
            const CoordGroup &group1 = groups_array[jgroup];
            const Parameters::Array &params1 = molparams_array[jgroup];

            //check first that these two CoordGroups could be within cutoff
            //(don't test igroup==jgroup as this is the same CutGroup
            // and definitely within cutoff!)
            const bool within_cutoff = (igroup == jgroup) or not
                                        spce->beyond(switchfunc->cutoffDistance(), 
                                                     aabox0, group1.aaBox());
            
            if (not within_cutoff)
                //this CutGroup is beyond the cutoff distance
                continue;
            
            //calculate all of the interatomic distances
            const double mindist = spce->calcDist(group0, group1, distmat);
            
            if (mindist > switchfunc->cutoffDistance())
                //all of the atoms are definitely beyond cutoff
                continue;
                
            CGIdx cgidx_jgroup = mol.cgIdx(jgroup);
                
            //get the non-bonded scale factors for all pairs of atoms
            //between these groups (or within this group, if igroup == jgroup)
            const CLJNBPairs::CGPairs &group_pairs = nbpairs(cgidx_igroup,
                                                             cgidx_jgroup);

            double icnrg = 0;
            double iljnrg = 0;
            
            //loop over all intraatomic pairs and calculate the energies
            const quint32 nats1 = group1.count();
            const Parameter *params1_array = params1.constData();
            
            calculateEnergy(group_pairs, distmat, params0_array, params1_array,
                            nats0, nats1, icnrg, iljnrg);
            
            //if this is the same group then half the energies to 
            //correct for double-counting
            if (igroup == jgroup)
            {
                icnrg *= 0.5;
                iljnrg *= 0.5;
            }

            //now add these energies onto the total for the molecule,
            //scaled by any non-bonded feather factor if using the group-based cutoff
            if (not (use_electrostatic_shifting or use_reaction_field or use_atomistic_cutoff))
            {
                if (mindist > switchfunc->electrostaticFeatherDistance())
                {
                    cnrg += switchfunc->electrostaticScaleFactor( Length(mindist) ) * icnrg;
                }
                else
                {
                    cnrg += icnrg;
                }
                
                if (mindist > switchfunc->vdwFeatherDistance())
                {
                    ljnrg += switchfunc->vdwScaleFactor( Length(mindist) ) * iljnrg;
                }
                else
                {
                    ljnrg += iljnrg;
                }
            }
            else
            {
                cnrg += icnrg;
                ljnrg += iljnrg;
            }
        }
    }
    
    //add this molecule pair's energy onto the total
    energy += Energy(scale_energy * cnrg, scale_energy * ljnrg);
}

/** Calculate the intramolecular CLJ energy of the passed molecule
    interacting with the rest of the same molecule in 'rest_of_mol', and
    add this onto 'energy'. This uses the passed workspace when
    performing the calculation. Note that mol and rest_of_mol should
    not contain any overlapping atoms, and that they should both be
    part of the same molecule (albeit potentially at different versions,
    but with the same layout UID)
    
    \throw SireError::incompatible_error
*/
void IntraCLJPotential::calculateEnergy(const IntraCLJPotential::Molecule &mol,
                                        const IntraCLJPotential::Molecule &rest_of_mol,
                                        IntraCLJPotential::Energy &energy,
                                        IntraCLJPotential::EnergyWorkspace &distmat,
                                        double scale_energy) const
{
    if (scale_energy == 0 or mol.isEmpty() or rest_of_mol.isEmpty())
        return;

    //ensure that this is the same molecule, with the same layout UID
    this->assertCompatible(mol, rest_of_mol);

    const quint32 ngroups0 = mol.nCutGroups();
    const CoordGroup *groups0_array = mol.coordinates().constData();
    
    BOOST_ASSERT( mol.parameters().atomicParameters().count() == int(ngroups0) );
    const Parameters::Array *molparams0_array 
                            = mol.parameters().atomicParameters().constData();

    const quint32 ngroups1 = rest_of_mol.nCutGroups();
    const CoordGroup *groups1_array = rest_of_mol.coordinates().constData();
    
    BOOST_ASSERT( rest_of_mol.parameters().atomicParameters().count() == int(ngroups1) );
    const Parameters::Array *molparams1_array
                            = rest_of_mol.parameters().atomicParameters().constData();

    //the CLJNBPairs must be the same in both molecules - this is checked
    //as part of assertCompatible(..)
    const CLJNBPairs &nbpairs = mol.parameters().intraScaleFactors();
    
    double cnrg = 0;
    double ljnrg = 0;
    
    //calculate the energy of all of the atoms in 'mol' interacting with
    //all of the atoms in 'rest_of_mol' that aren't in 'mol'
    for (quint32 igroup=0; igroup<ngroups0; ++igroup)
    {
        const Parameters::Array &params0 = molparams0_array[igroup];

        const CoordGroup &group0 = groups0_array[igroup];
        const AABox &aabox0 = group0.aaBox();
        const quint32 nats0 = group0.count();
        const Parameter *params0_array = params0.constData();
    
        //get the CGIdx of this CutGroup
        CGIdx cgidx_igroup = mol.cgIdx(igroup);
    
        for (quint32 jgroup=0; jgroup<ngroups1; ++jgroup)
        {
            const CoordGroup &group1 = groups1_array[jgroup];
            const Parameters::Array &params1 = molparams1_array[jgroup];

            CGIdx cgidx_jgroup = rest_of_mol.cgIdx(jgroup);

            //skip this CutGroup if it is in 'mol'
            if (mol.molecule().selection().selectedAll(cgidx_jgroup))
                continue;

            //check first that these two CoordGroups could be within cutoff
            //(don't test igroup==jgroup as this is the same CutGroup
            // and definitely within cutoff!)
            const bool within_cutoff = (cgidx_igroup == cgidx_jgroup) or not
                                        spce->beyond(switchfunc->cutoffDistance(), 
                                                     aabox0, group1.aaBox());
            
            if (not within_cutoff)
                //this CutGroup is beyond the cutoff distance
                continue;
            
            //calculate all of the interatomic distances
            const double mindist = spce->calcDist(group0, group1, distmat);
            
            if (mindist > switchfunc->cutoffDistance())
                //all of the atoms are definitely beyond cutoff
                continue;
                
            //get the non-bonded scale factors for all pairs of atoms
            //between these groups (or within this group, if igroup == jgroup)
            const CLJNBPairs::CGPairs &group_pairs = nbpairs(cgidx_igroup,
                                                             cgidx_jgroup);

            double icnrg = 0;
            double iljnrg = 0;
            
            //loop over all intraatomic pairs and calculate the energies
            const quint32 nats1 = group1.count();
            const Parameter *params1_array = params1.constData();
            
            if (cgidx_igroup == cgidx_jgroup
                or mol.molecule().selection().selected(cgidx_jgroup))
            {
                //this is the same CutGroup, or some of CutGroup jgroup
                //is present in 'mol' - we must be careful not to 
                //double-count the energy between atoms in both 'mol'
                //and 'rest_of_mol'
                
                //get the atoms from this CutGroup that are contained in 
                //each part of the molecule
                QSet<Index> atoms0 = mol.molecule()
                                        .selection().selectedAtoms(cgidx_igroup);
                                      
                QSet<Index> mol_atoms1 = atoms0;
                
                if (cgidx_igroup != cgidx_jgroup)
                    mol_atoms1 = mol.molecule().selection().selectedAtoms(cgidx_jgroup);
                                            
                QSet<Index> atoms1 = rest_of_mol.molecule()
                                        .selection().selectedAtoms(cgidx_jgroup);
                                        
                //remove from atoms1 atoms that are part of 'mol'
                atoms1 -= mol_atoms1;
                
                calculateEnergy(group_pairs, atoms0, atoms1, distmat,
                                params0_array, params1_array, 
                                nats0, nats1, icnrg, iljnrg);
            }
            else
            {
                calculateEnergy(group_pairs, distmat,
                                params0_array, params1_array,
                                nats0, nats1, icnrg, iljnrg);
            }

            //now add these energies onto the total for the molecule,
            //scaled by any non-bonded feather factor if using the group-based cutoff
            if (not (use_electrostatic_shifting or use_reaction_field or use_atomistic_cutoff))
            {
                if (mindist > switchfunc->electrostaticFeatherDistance())
                {
                    cnrg += switchfunc->electrostaticScaleFactor( Length(mindist) ) * icnrg;
                }
                else
                {
                    cnrg += icnrg;
                }
                
                if (mindist > switchfunc->vdwFeatherDistance())
                {
                    ljnrg += switchfunc->vdwScaleFactor( Length(mindist) ) * iljnrg;
                }
                else
                {
                    ljnrg += iljnrg;
                }
            }
            else
            {
                cnrg += icnrg;
                ljnrg += iljnrg;
            }
        }
    }
    
    //add the molecule's energy onto the total
    energy += Energy(scale_energy * cnrg, scale_energy * ljnrg);
}

void IntraCLJPotential::calculateForce(const CLJNBPairs::CGPairs &group_pairs,
                             const CoordGroup &group0, const CoordGroup &group1,
                             const double mindist,
                             IntraCLJPotential::ForceWorkspace &distmat,
                             const IntraCLJPotential::Parameter *params0_array,
                             const IntraCLJPotential::Parameter *params1_array,
                             const quint32 nats0, const quint32 nats1,
                             const double shift_coul,
                             Vector *group_forces0_array,
                             const double scale_force) const
{
    if (mindist > switchfunc->featherDistance())
    {
        //we need to calculate the forces taking into account
        //the derivative of the switching function!
                
        //calculate the switching scale factors and their 
        //derivatives
        const double scl_coul = switchfunc->electrostaticScaleFactor( Length(mindist) );
        const double scl_lj = switchfunc->vdwScaleFactor( Length(mindist) );
                
        Vector group_sep = (group1.aaBox().center() -
                            group0.aaBox().center()).normalise(); 
                
        Vector dscl_coul = switchfunc->dElectrostaticScaleFactor( Length(mindist) ) 
                                     * group_sep;
                                     
        Vector dscl_lj = switchfunc->dVDWScaleFactor( Length(mindist) )
                                     * group_sep;

        if (group_pairs.isEmpty())
        {
            //there are no scale factors between atoms in these groups
            for (quint32 i=0; i<nats0; ++i)
            {
                distmat.setOuterIndex(i);
                const Parameter &param0 = params0_array[i];

                Vector total_force;
                
                if (param0.ljid == 0)
                {
                    //null LJ parameter - only add on the coulomb energy
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        const double q2 = param0.reduced_charge *
                                          params1_array[j].reduced_charge;
                                
                        if (q2 != 0)
                        {
                            //calculate the coulomb energy
                            const double cnrg = scl_coul * q2 /
                                                distmat[j].length();

                            //calculate the coulomb force
                            Vector cforce = (-cnrg / distmat[j].length() *
                                                     distmat[j].direction()) +
                                             
                                                    ((cnrg-shift_coul) * dscl_coul);

                            total_force += cforce;
                        }
                    }
                }
                else
                {
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        //do both coulomb and LJ
                        const Parameter &param1 = params1_array[j];
                        
                        const double invdist = double(1) / distmat[j].length();
                        
                        Vector force;

                        const double q2 = param0.reduced_charge *
                                          param1.reduced_charge;
                                
                        if (q2 != 0)
                        {
                            //calculate the energy
                            const double cnrg = scl_coul * q2 * invdist;

                            //calculate the force
                            force = (-cnrg / distmat[j].length() *
                                             distmat[j].direction()) +
                                             
                                             ((cnrg-shift_coul) * dscl_coul);
                        }
                              
                        if (param1.ljid != 0)
                        {
                            const LJPair &ljpair = ljpairs.constData()[
                                                    ljpairs.map(param0.ljid,
                                                                param1.ljid)];
                        
                            double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                            double sig_over_dist12 = pow_2(sig_over_dist6);

                            //calculate the energy
                            const double ljnrg = 4 * ljpair.epsilon() *
                                                   (sig_over_dist12 - sig_over_dist6);

                            // dU/dr requires an extra power of r
                            sig_over_dist6 *= invdist;
                            sig_over_dist12 *= invdist;

                            force += ((scl_lj * 4 * ljpair.epsilon() * 
                                       (6.0*sig_over_dist6 - 12.0*sig_over_dist12))
                                           * distmat[j].direction())
                                            
                                     + (ljnrg * dscl_lj);
                        }
                        
                        total_force += force;
                    }
                }
                        
                group_forces0_array[i] += scale_force * total_force;
            }
        }
        else
        {
            //there are different nb scale factors between
            //the atoms. We need to calculate the forces
            //using them...
            for (quint32 i=0; i<nats0; ++i)
            {
                distmat.setOuterIndex(i);
                const Parameter &param0 = params0_array[i];
                
                Vector total_force;
                
                if (param0.ljid == 0)
                {
                    //null LJ parameter - only add on the coulomb energy
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        const CLJScaleFactor &cljscl = group_pairs(i,j);
                            
                        if (cljscl.coulomb() != 0)
                        {
                            const double q2 = param0.reduced_charge * 
                                              params1_array[j].reduced_charge;
                                                      
                            if (q2 != 0)
                            {
                                //calculate the coulomb energy
                                const double cnrg = cljscl.coulomb() *
                                                      scl_coul * q2
                                                      / distmat[j].length();
                                               
                                //calculate the coulomb force
                                Vector cforce = (-cnrg 
                                                  / distmat[j].length() *
                                                    distmat[j].direction()) +
                                             
                                                  ((cnrg-shift_coul) * dscl_coul);
                            
                                total_force += cforce;
                            }
                        }
                    }
                }
                else
                {
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        //do both coulomb and LJ
                        const CLJScaleFactor &cljscl = group_pairs(i,j);
                                
                        if (cljscl.coulomb() != 0 or cljscl.lj() != 0)
                        {
                            const Parameter &param1 = params1_array[j];
                        
                            const double invdist = double(1) 
                                                     / distmat[j].length();
                        
                            Vector force;
                                    
                            const double q2 = param0.reduced_charge *
                                              param1.reduced_charge;

                            if (q2 != 0)
                            {
                                //calculate the energy
                                const double cnrg = cljscl.coulomb() *
                                                            scl_coul * q2 * invdist;
                        
                                //calculate the force
                                force = (scl_coul * -cnrg 
                                           / distmat[j].length() *
                                             distmat[j].direction()) +
                                             
                                        ((cnrg-shift_coul) * dscl_coul);
                            }
                            
                            if (cljscl.lj() != 0 and param1.ljid != 0)
                            {
                                const LJPair &ljpair = ljpairs.constData()[
                                                        ljpairs.map(param0.ljid,
                                                                    param1.ljid)];
                            
                                double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                                double sig_over_dist12 = pow_2(sig_over_dist6);

                                //calculate the energy
                                const double ljnrg = cljscl.lj() *
                                                      4 * ljpair.epsilon() *
                                                   (sig_over_dist12 - sig_over_dist6);

                                // dU/dr requires an extra power of r
                                sig_over_dist6 *= invdist;
                                sig_over_dist12 *= invdist;

                                force += ((scl_lj * 4 * ljpair.epsilon() * 
                                         (6.0*sig_over_dist6 - 12.0*sig_over_dist12))
                                          * distmat[j].direction())
                                            
                                          + (ljnrg * dscl_lj);
                            }
                        
                            total_force += force;
                        }
                    }
                }
                        
                group_forces0_array[i] += scale_force * total_force;
            }
        }
    }
    else
    {
        //not in the feather region, so can calculate the forces
        //directly (also, no need to calculate shift, as 
        //the shifting function is constant, so does not
        //affect the gradient)
                
        if (group_pairs.isEmpty())
        {
            //no nb scale factors to worry about
            for (quint32 i=0; i<nats0; ++i)
            {
                distmat.setOuterIndex(i);
                const Parameter &param0 = params0_array[i];
                
                Vector total_force;
                
                if (param0.ljid == 0)
                {
                    //null LJ parameter - only add on the coulomb energy
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        //calculate the coulomb force
                        Vector cforce = -(param0.reduced_charge *
                                          params1_array[j].reduced_charge / 
                                          distmat[j].length2()) *
                                             
                                          distmat[j].direction();

                        total_force += cforce;
                    }
                }
                else
                {
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        //do both coulomb and LJ
                        const Parameter &param1 = params1_array[j];
                            
                        const double invdist = double(1) / distmat[j].length();
                        const double invdist2 = pow_2(invdist);
                        
                        //calculate the force
                        Vector force = -(param0.reduced_charge * 
                                         param1.reduced_charge * invdist2) 
                                            
                                         * distmat[j].direction();
                              
                        if (param1.ljid != 0)
                        {
                            const LJPair &ljpair = ljpairs.constData()[
                                                    ljpairs.map(param0.ljid,
                                                                param1.ljid)];
                        
                            double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                            double sig_over_dist12 = pow_2(sig_over_dist6);

                            // dU/dr requires an extra power of r
                            sig_over_dist6 *= invdist;
                            sig_over_dist12 *= invdist;

                            force += (4 * ljpair.epsilon() * (6.0*sig_over_dist6 - 
                                                             12.0*sig_over_dist12))
                                      * distmat[j].direction();
                        }

                        total_force += force;
                    }
                }
                        
                group_forces0_array[i] += scale_force * total_force;
            }
        }
        else
        {
            //there are different nb scale factors between
            //different atoms...
            for (quint32 i=0; i<nats0; ++i)
            {
                distmat.setOuterIndex(i);
                const Parameter &param0 = params0_array[i];
                
                Vector total_force;
                
                if (param0.ljid == 0)
                {
                    //null LJ parameter - only add on the coulomb energy
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        const CLJScaleFactor &cljscl = group_pairs(i,j);
                            
                        if (cljscl.coulomb() != 0)
                        {
                            //calculate the coulomb force
                            Vector cforce = -(cljscl.coulomb() *
                                              param0.reduced_charge *
                                              params1_array[j].reduced_charge / 
                                              distmat[j].length2()) *
                                             
                                              distmat[j].direction();

                            total_force += cforce;
                        }
                    }
                }
                else
                {
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        //do both coulomb and LJ
                        const CLJScaleFactor &cljscl = group_pairs(i,j);
                                
                        if (cljscl.coulomb() != 0 or cljscl.lj() != 0)
                        {
                            const Parameter &param1 = params1_array[j];
                            
                            const double invdist = double(1)    
                                                    / distmat[j].length();
                            const double invdist2 = pow_2(invdist);
                        
                            //calculate the force
                            Vector force = -(cljscl.coulomb() *
                                             param0.reduced_charge * 
                                             param1.reduced_charge * invdist2) 
                                            
                                             * distmat[j].direction();
                              
                            if (cljscl.lj() != 0 and param1.ljid != 0)
                            {
                                const LJPair &ljpair = ljpairs.constData()[
                                                         ljpairs.map(param0.ljid,
                                                                     param1.ljid)];
                        
                                double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                                double sig_over_dist12 = pow_2(sig_over_dist6);

                                // dU/dr requires an extra power of r
                                sig_over_dist6 *= invdist;
                                sig_over_dist12 *= invdist;

                                force += (cljscl.lj() *
                                          4 * ljpair.epsilon() * (6.0*sig_over_dist6 - 
                                                                 12.0*sig_over_dist12))
                                           * distmat[j].direction();
                            }

                            total_force += force;
                        }
                    }
                }
                        
                group_forces0_array[i] += scale_force * total_force;

            } // end of loop over i atoms

        } // end of whether there are intra scale factors

    } // end of whether within feather region
}

void IntraCLJPotential::calculateForce(const CLJNBPairs::CGPairs &group_pairs,
                             const QSet<Index> &atoms0, const QSet<Index> &atoms1,
                             const CoordGroup &group0, const CoordGroup &group1,
                             const double mindist,
                             IntraCLJPotential::ForceWorkspace &distmat,
                             const IntraCLJPotential::Parameter *params0_array,
                             const IntraCLJPotential::Parameter *params1_array,
                             const double shift_coul,
                             Vector *group_forces0_array,
                             const double scale_force) const
{
    if (atoms0.isEmpty() or atoms1.isEmpty())
        return;

    if (mindist > switchfunc->featherDistance())
    {
        //we need to calculate the forces taking into account
        //the derivative of the switching function!
                
        //calculate the switching scale factors and their 
        //derivatives
        const double scl_coul = switchfunc->electrostaticScaleFactor( Length(mindist) );
        const double scl_lj = switchfunc->vdwScaleFactor( Length(mindist) );
                
        Vector group_sep = (group1.aaBox().center() -
                            group0.aaBox().center()).normalise(); 
                
        Vector dscl_coul = switchfunc->dElectrostaticScaleFactor( Length(mindist) ) 
                                     * group_sep;
                                     
        Vector dscl_lj = switchfunc->dVDWScaleFactor( Length(mindist) )
                                     * group_sep;

        if (group_pairs.isEmpty())
        {
            //there are no scale factors between atoms in these groups
            foreach (Index i, atoms0)
            {
                distmat.setOuterIndex(i);
                const Parameter &param0 = params0_array[i];

                Vector total_force;
                
                if (param0.ljid == 0)
                {
                    //null LJ parameter - only add on the coulomb energy
                    foreach (Index j, atoms1)
                    {
                        const double q2 = param0.reduced_charge *
                                          params1_array[j].reduced_charge;
                                
                        if (q2 != 0)
                        {
                            //calculate the coulomb energy
                            const double cnrg = scl_coul * q2 /
                                                distmat[j].length();

                            //calculate the coulomb force
                            Vector cforce = (-cnrg / distmat[j].length() *
                                                     distmat[j].direction()) +
                                             
                                                    ((cnrg-shift_coul) * dscl_coul);

                            total_force += cforce;
                        }
                    }
                }
                else
                {
                    foreach (Index j, atoms1)
                    {
                        //do both coulomb and LJ
                        const Parameter &param1 = params1_array[j];
                        
                        const double invdist = double(1) / distmat[j].length();
                        
                        Vector force;

                        const double q2 = param0.reduced_charge *
                                          param1.reduced_charge;
                                
                        if (q2 != 0)
                        {
                            //calculate the energy
                            const double cnrg = scl_coul * q2 * invdist;

                            //calculate the force
                            force = (-cnrg / distmat[j].length() *
                                             distmat[j].direction()) +
                                             
                                             ((cnrg-shift_coul) * dscl_coul);
                        }
                              
                        if (param1.ljid != 0)
                        {
                            const LJPair &ljpair = ljpairs.constData()[
                                                    ljpairs.map(param0.ljid,
                                                                param1.ljid)];
                        
                            double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                            double sig_over_dist12 = pow_2(sig_over_dist6);

                            //calculate the energy
                            const double ljnrg = 4 * ljpair.epsilon() *
                                                   (sig_over_dist12 - sig_over_dist6);

                            // dU/dr requires an extra power of r
                            sig_over_dist6 *= invdist;
                            sig_over_dist12 *= invdist;

                            force += ((scl_lj * 4 * ljpair.epsilon() * 
                                       (6.0*sig_over_dist6 - 12.0*sig_over_dist12))
                                           * distmat[j].direction())
                                            
                                     + (ljnrg * dscl_lj);
                        }
                        
                        total_force += force;
                    }
                }
                        
                group_forces0_array[i] += scale_force * total_force;
            }
        }
        else
        {
            //there are different nb scale factors between
            //the atoms. We need to calculate the forces
            //using them...
            foreach (Index i, atoms0)
            {
                distmat.setOuterIndex(i);
                const Parameter &param0 = params0_array[i];
                
                Vector total_force;
                
                if (param0.ljid == 0)
                {
                    //null LJ parameter - only add on the coulomb energy
                    foreach (Index j, atoms1)
                    {
                        const CLJScaleFactor &cljscl = group_pairs(i,j);
                            
                        if (cljscl.coulomb() != 0)
                        {
                            const double q2 = param0.reduced_charge * 
                                              params1_array[j].reduced_charge;
                                                      
                            if (q2 != 0)
                            {
                                //calculate the coulomb energy
                                const double cnrg = cljscl.coulomb() *
                                                      scl_coul * q2
                                                      / distmat[j].length();
                                               
                                //calculate the coulomb force
                                Vector cforce = (-cnrg 
                                                  / distmat[j].length() *
                                                    distmat[j].direction()) +
                                             
                                                  ((cnrg-shift_coul) * dscl_coul);
                            
                                total_force += cforce;
                            }
                        }
                    }
                }
                else
                {
                    foreach (Index j, atoms1)
                    {
                        //do both coulomb and LJ
                        const CLJScaleFactor &cljscl = group_pairs(i,j);
                                
                        if (cljscl.coulomb() != 0 or cljscl.lj() != 0)
                        {
                            const Parameter &param1 = params1_array[j];
                        
                            const double invdist = double(1) 
                                                     / distmat[j].length();
                        
                            Vector force;
                                    
                            const double q2 = param0.reduced_charge *
                                              param1.reduced_charge;

                            if (q2 != 0)
                            {
                                //calculate the energy
                                const double cnrg = cljscl.coulomb() *
                                                            scl_coul * q2 * invdist;
                        
                                //calculate the force
                                force = (scl_coul * -cnrg 
                                           / distmat[j].length() *
                                             distmat[j].direction()) +
                                             
                                        ((cnrg-shift_coul) * dscl_coul);
                            }
                            
                            if (cljscl.lj() != 0 and param1.ljid != 0)
                            {
                                const LJPair &ljpair = ljpairs.constData()[
                                                        ljpairs.map(param0.ljid,
                                                                    param1.ljid)];
                            
                                double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                                double sig_over_dist12 = pow_2(sig_over_dist6);

                                //calculate the energy
                                const double ljnrg = cljscl.lj() *
                                                      4 * ljpair.epsilon() *
                                                   (sig_over_dist12 - sig_over_dist6);

                                // dU/dr requires an extra power of r
                                sig_over_dist6 *= invdist;
                                sig_over_dist12 *= invdist;

                                force += ((scl_lj * 4 * ljpair.epsilon() * 
                                         (6.0*sig_over_dist6 - 12.0*sig_over_dist12))
                                          * distmat[j].direction())
                                            
                                          + (ljnrg * dscl_lj);
                            }
                        
                            total_force += force;
                        }
                    }
                }
                        
                group_forces0_array[i] += scale_force * total_force;
            }
        }
    }
    else
    {
        //not in the feather region, so can calculate the forces
        //directly (also, no need to calculate shift, as 
        //the shifting function is constant, so does not
        //affect the gradient)
                
        if (group_pairs.isEmpty())
        {
            //no nb scale factors to worry about
            foreach (Index i, atoms0)
            {
                distmat.setOuterIndex(i);
                const Parameter &param0 = params0_array[i];
                
                Vector total_force;
                
                if (param0.ljid == 0)
                {
                    //null LJ parameter - only add on the coulomb energy
                    foreach (Index j, atoms1)
                    {
                        //calculate the coulomb force
                        Vector cforce = -(param0.reduced_charge *
                                          params1_array[j].reduced_charge / 
                                          distmat[j].length2()) *
                                             
                                          distmat[j].direction();

                        total_force += cforce;
                    }
                }
                else
                {
                    foreach (Index j, atoms1)
                    {
                        //do both coulomb and LJ
                        const Parameter &param1 = params1_array[j];
                            
                        const double invdist = double(1) / distmat[j].length();
                        const double invdist2 = pow_2(invdist);
                        
                        //calculate the force
                        Vector force = -(param0.reduced_charge * 
                                         param1.reduced_charge * invdist2) 
                                            
                                         * distmat[j].direction();
                              
                        if (param1.ljid != 0)
                        {
                            const LJPair &ljpair = ljpairs.constData()[
                                                    ljpairs.map(param0.ljid,
                                                                param1.ljid)];
                        
                            double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                            double sig_over_dist12 = pow_2(sig_over_dist6);

                            // dU/dr requires an extra power of r
                            sig_over_dist6 *= invdist;
                            sig_over_dist12 *= invdist;

                            force += (4 * ljpair.epsilon() * (6.0*sig_over_dist6 - 
                                                             12.0*sig_over_dist12))
                                      * distmat[j].direction();
                        }

                        total_force += force;
                    }
                }
                        
                group_forces0_array[i] += scale_force * total_force;
            }
        }
        else
        {
            //there are different nb scale factors between
            //different atoms...
            foreach (Index i, atoms0)
            {
                distmat.setOuterIndex(i);
                const Parameter &param0 = params0_array[i];
                
                Vector total_force;
                
                if (param0.ljid == 0)
                {
                    //null LJ parameter - only add on the coulomb energy
                    foreach (Index j, atoms1)
                    {
                        const CLJScaleFactor &cljscl = group_pairs(i,j);
                            
                        if (cljscl.coulomb() != 0)
                        {
                            //calculate the coulomb force
                            Vector cforce = -(cljscl.coulomb() *
                                              param0.reduced_charge *
                                              params1_array[j].reduced_charge / 
                                              distmat[j].length2()) *
                                             
                                              distmat[j].direction();

                            total_force += cforce;
                        }
                    }
                }
                else
                {
                    foreach (Index j, atoms1)
                    {
                        //do both coulomb and LJ
                        const CLJScaleFactor &cljscl = group_pairs(i,j);
                                
                        if (cljscl.coulomb() != 0 or cljscl.lj() != 0)
                        {
                            const Parameter &param1 = params1_array[j];
                            
                            const double invdist = double(1)    
                                                    / distmat[j].length();
                            const double invdist2 = pow_2(invdist);
                        
                            //calculate the force
                            Vector force = -(cljscl.coulomb() *
                                             param0.reduced_charge * 
                                             param1.reduced_charge * invdist2) 
                                            
                                             * distmat[j].direction();
                              
                            if (cljscl.lj() != 0 and param1.ljid != 0)
                            {
                                const LJPair &ljpair = ljpairs.constData()[
                                                         ljpairs.map(param0.ljid,
                                                                     param1.ljid)];
                        
                                double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                                double sig_over_dist12 = pow_2(sig_over_dist6);

                                // dU/dr requires an extra power of r
                                sig_over_dist6 *= invdist;
                                sig_over_dist12 *= invdist;

                                force += (cljscl.lj() *
                                          4 * ljpair.epsilon() * (6.0*sig_over_dist6 - 
                                                                 12.0*sig_over_dist12))
                                           * distmat[j].direction();
                            }

                            total_force += force;
                        }
                    }
                }
                        
                group_forces0_array[i] += scale_force * total_force;

            } // end of loop over i atoms

        } // end of whether there are intra scale factors

    } // end of whether within feather region
}

/** Calculate the coulomb and LJ forces between the atoms in the molecule 'mol'
    and add these forces onto 'forces'. This uses
    the passed workspace to perform the calculation */
void IntraCLJPotential::calculateForce(const IntraCLJPotential::Molecule &mol,
                                       MolForceTable &forces, 
                                       IntraCLJPotential::ForceWorkspace &distmat,
                                       double scale_force) const
{
    if (scale_force == 0 or mol.isEmpty())
        return;
    
    const quint32 ngroups = mol.nCutGroups();
    
    const CoordGroup *groups_array = mol.coordinates().constData();
    
    const Parameters::Array *molparams_array 
                            = mol.parameters().atomicParameters().constData();
    
    BOOST_ASSERT(forces.nCutGroups() == mol.molecule().data().info().nCutGroups());
    BOOST_ASSERT(forces.molNum() == mol.molecule().data().number());
    
    const MolForceTable::Array *forces_array = forces.constData();

    const CLJNBPairs &nbpairs = mol.parameters().intraScaleFactors();
    
    //loop over all pairs of CutGroups in the molecule
    for (quint32 igroup=0; igroup<ngroups; ++igroup)
    {
        const Parameters::Array &params0 = molparams_array[igroup];

        const CoordGroup &group0 = groups_array[igroup];
        const AABox &aabox0 = group0.aaBox();
        const quint32 nats0 = group0.count();
        const Parameter *params0_array = params0.constData();
    
        //get the CGIdx of this CutGroup
        CGIdx cgidx_igroup = mol.cgIdx(igroup);
        
        //now get the index of the force table for this CutGroup
        int force_idx = forces.map(cgidx_igroup);
        
        if (force_idx == -1)
            //we are not interested in the forces on this CutGroup
            continue;
        
        //get the table that holds the forces acting on all of the
        //atoms of this CutGroup
        BOOST_ASSERT(forces_array[force_idx].count() == int(nats0));
        Vector *group_forces0_array = forces.data(force_idx);
    
        //get the forces acting on this group from all other groups
        //(yes, we are doing a full n2 loop, and not taking advantage
        // of equal and opposite forces)
        for (quint32 jgroup=0; jgroup<ngroups; ++jgroup)
        {
            const CoordGroup &group1 = groups_array[jgroup];
            const Parameters::Array &params1 = molparams_array[jgroup];

            //check first that these two CoordGroups could be within cutoff
            const bool outside_cutoff = igroup != jgroup and 
                                        spce->beyond(switchfunc->cutoffDistance(), 
                                                     aabox0, group1.aaBox());
            
            if (outside_cutoff)
                //this CutGroup is beyond the cutoff distance
                continue;
            
            //calculate all of the interatomic distances
            const double mindist = spce->calcDistVectors(group0, group1, distmat);
            
            if (mindist > switchfunc->cutoffDistance())
                //all of the atoms are definitely beyond cutoff
                continue;

            const quint32 nats1 = group1.count();
            
            //loop over all interatomic pairs and calculate the energies
            const Parameter *params1_array = params1.constData();

            CGIdx cgidx_jgroup = mol.cgIdx(jgroup);

            //get the non-bonded scale factors for all pairs of atoms
            //between these two groups (or within this group, if igroup == jgroup)
            const CLJNBPairs::CGPairs &group_pairs = nbpairs(cgidx_igroup,
                                                             cgidx_jgroup);

            double shift_coul = 0;

            if (use_electrostatic_shifting and igroup != jgroup)
                shift_coul = this->totalCharge(params0) * this->totalCharge(params1)
                              / switchfunc->electrostaticCutoffDistance();

            //now calculate the forces acting on group0 caused by group1
            calculateForce(group_pairs, group0, group1,
                           mindist, distmat, params0_array, params1_array,
                           nats0, nats1, shift_coul, group_forces0_array,
                           scale_force);
             
        } // end of loop over CutGroups (jgroup)

    } // end of loop over CutGroups (igroup)
}

/** Calculate the total forces acting on the atoms in 'mol' caused by the 
    other atoms in the same molecule contained in 'rest_of_mol'. This calculates
    the forces and adds them onto 'forces' (which are for 'mol'). Note that they must
    use the same layout UID and same intra-nonbonded scaling factors
    
    \throw SireError::incompatible_error
*/
void IntraCLJPotential::calculateForce(const IntraCLJPotential::Molecule &mol,
                                       const IntraCLJPotential::Molecule &rest_of_mol,
                                       MolForceTable &forces,
                                       IntraCLJPotential::ForceWorkspace &distmat,
                                       double scale_force) const
{
    if (scale_force == 0 or mol.isEmpty() or rest_of_mol.isEmpty())
        return;
        
    const quint32 ngroups0 = mol.nCutGroups();
    const CoordGroup *groups0_array = mol.coordinates().constData();
    
    const Parameters::Array *molparams0_array
                            = mol.parameters().atomicParameters().constData();
                            
    BOOST_ASSERT(forces.nCutGroups() == mol.molecule().data().info().nCutGroups());
    BOOST_ASSERT(forces.molNum() == mol.molecule().data().number());
    
    const MolForceTable::Array *forces_array = forces.constData();
    
    this->assertCompatible(mol, rest_of_mol);
    
    const quint32 ngroups1 = rest_of_mol.nCutGroups();
    const CoordGroup *groups1_array = rest_of_mol.coordinates().constData();
    
    const Parameters::Array *molparams1_array
                            = rest_of_mol.parameters().atomicParameters().constData();

    //this is now guaranteed to be the same in both passed parts of 
    //the molecule
    const CLJNBPairs &nbpairs = mol.parameters().intraScaleFactors();

    //calculate the forces acting on the atoms in 'mol' caused by 
    //the atoms in 'rest_of_mol'
    for (quint32 igroup=0; igroup<ngroups0; ++igroup)
    {
        const Parameters::Array &params0 = molparams0_array[igroup];

        const CoordGroup &group0 = groups0_array[igroup];
        const AABox &aabox0 = group0.aaBox();
        const quint32 nats0 = group0.count();
        const Parameter *params0_array = params0.constData();
    
        //get the CGIdx of this CutGroup
        CGIdx cgidx_igroup = mol.cgIdx(igroup);
        
        //now get the index of the force table for this CutGroup
        int force_idx = forces.map(cgidx_igroup);
        
        if (force_idx == -1)
            //we are not interested in the forces on this CutGroup
            continue;
        
        //get the table that holds the forces acting on all of the
        //atoms of this CutGroup
        BOOST_ASSERT(forces_array[force_idx].count() == int(nats0));
        Vector *group_forces0_array = forces.data(force_idx);
        
        for (quint32 jgroup=0; jgroup<ngroups1; ++jgroup)
        {
            CGIdx cgidx_jgroup = rest_of_mol.cgIdx(jgroup);

            if (mol.molecule().selection().selectedAll(cgidx_jgroup))
                //all of this CutGroup is in 'mol', so don't evaluate
                //the force from it
                continue;

            const CoordGroup &group1 = groups1_array[jgroup];
            const Parameters::Array &params1 = molparams1_array[jgroup];

            //check first that these two CoordGroups could be within cutoff
            const bool outside_cutoff = cgidx_igroup != cgidx_jgroup and 
                                        spce->beyond(switchfunc->cutoffDistance(), 
                                                     aabox0, group1.aaBox());
            
            if (outside_cutoff)
                //this CutGroup is beyond the cutoff distance
                continue;
            
            //calculate all of the interatomic distances
            const double mindist = spce->calcDistVectors(group0, group1, distmat);
            
            if (mindist > switchfunc->cutoffDistance())
                //all of the atoms are definitely beyond cutoff
                continue;

            const quint32 nats1 = group1.count();
            
            //loop over all interatomic pairs and calculate the energies
            const Parameter *params1_array = params1.constData();

            //get the non-bonded scale factors for all pairs of atoms
            //between these two groups (or within this group, if igroup == jgroup)
            const CLJNBPairs::CGPairs &group_pairs = nbpairs(cgidx_igroup,
                                                             cgidx_jgroup);
            
            if (cgidx_igroup == cgidx_jgroup or
                mol.molecule().selection().selected(cgidx_jgroup))
            {
                //some of the atoms in jgroup are selected as part of 'mol'.
                //We must be careful not to double-count the interactions
                //with these atoms
                
                QSet<Index> atoms0 = mol.molecule().selection()
                                               .selectedAtoms(cgidx_igroup);
                                  
                QSet<Index> mol_atoms1 = atoms0;
                
                if (cgidx_igroup != cgidx_jgroup)
                    mol_atoms1 = mol.molecule().selection().selectedAtoms(cgidx_jgroup);
                                                                         
                QSet<Index> atoms1 = rest_of_mol.molecule().selection()
                                               .selectedAtoms(cgidx_jgroup);
                                         
                //remove the atoms from 'rest_of_mol' that are part of 'mol'
                atoms1 -= mol_atoms1;
                                                           
                double shift_coul = 0;
            
                if (use_electrostatic_shifting and cgidx_igroup != cgidx_jgroup)
                    shift_coul = this->totalCharge(params0) * this->totalCharge(params1)
                                / switchfunc->electrostaticCutoffDistance();
                
                calculateForce(group_pairs, atoms0, atoms1,
                               group0, group1, mindist, distmat,
                               params0_array, params1_array,
                               shift_coul,
                               group_forces0_array, scale_force);
            }
            else
            {
                double shift_coul = 0;
            
                if (use_electrostatic_shifting)
                    shift_coul = this->totalCharge(params0) * this->totalCharge(params1)
                                / switchfunc->electrostaticCutoffDistance();

                calculateForce(group_pairs, group0, group1,
                               mindist, distmat,
                               params0_array, params1_array,
                               nats0, nats1, shift_coul,
                               group_forces0_array, scale_force);
            }
        }
    }
}

void IntraCLJPotential::calculateCoulombForce(const CLJNBPairs::CGPairs &group_pairs,
                             const CoordGroup &group0, const CoordGroup &group1,
                             const double mindist,
                             IntraCLJPotential::ForceWorkspace &distmat,
                             const IntraCLJPotential::Parameter *params0_array,
                             const IntraCLJPotential::Parameter *params1_array,
                             const quint32 nats0, const quint32 nats1,
                             const double shift_coul,
                             Vector *group_forces0_array,
                             const double scale_force) const
{
    if (mindist > switchfunc->featherDistance())
    {
        //we need to calculate the forces taking into account
        //the derivative of the switching function!
                
        //calculate the switching scale factors and their 
        //derivatives
        const double scl_coul = switchfunc->electrostaticScaleFactor( Length(mindist) );
                
        Vector group_sep = (group1.aaBox().center() -
                            group0.aaBox().center()).normalise(); 
                
        Vector dscl_coul = switchfunc->dElectrostaticScaleFactor( Length(mindist) ) 
                              * group_sep;

        if (group_pairs.isEmpty())
        {
            //there are no scale factors between atoms in these groups
            for (quint32 i=0; i<nats0; ++i)
            {
                const Parameter &param0 = params0_array[i];

                if (param0.reduced_charge == 0)
                    return;

                distmat.setOuterIndex(i);

                Vector total_force;
                
                for (quint32 j=0; j<nats1; ++j)
                {
                    const double q2 = param0.reduced_charge *
                                      params1_array[j].reduced_charge;
                                
                    if (q2 != 0)
                    {
                        //calculate the coulomb energy
                        const double cnrg = scl_coul * q2 /
                                            distmat[j].length();

                        //calculate the coulomb force
                        Vector cforce = (-cnrg / distmat[j].length() *
                                                 distmat[j].direction()) +
                                             
                                                ((cnrg-shift_coul) * dscl_coul);

                        total_force += cforce;
                    }
                }
                        
                group_forces0_array[i] += scale_force * total_force;
            }
        }
        else
        {
            //there are different nb scale factors between
            //the atoms. We need to calculate the forces
            //using them...
            for (quint32 i=0; i<nats0; ++i)
            {
                const Parameter &param0 = params0_array[i];

                if (param0.reduced_charge == 0)
                    continue;
                
                distmat.setOuterIndex(i);

                Vector total_force;
                
                for (quint32 j=0; j<nats1; ++j)
                {
                    const CLJScaleFactor &cljscl = group_pairs(i,j);
                            
                    if (cljscl.coulomb() != 0)
                    {
                        const double q2 = param0.reduced_charge * 
                                          params1_array[j].reduced_charge;
                                                      
                        if (q2 != 0)
                        {
                            //calculate the coulomb energy
                            const double cnrg = cljscl.coulomb() *
                                                  scl_coul * q2
                                                  / distmat[j].length();
                                               
                            //calculate the coulomb force
                            Vector cforce = (-cnrg 
                                             / distmat[j].length() *
                                               distmat[j].direction()) +
                                             
                                              ((cnrg-shift_coul) * dscl_coul);
                            
                            total_force += cforce;
                        }
                    }
                }
                        
                group_forces0_array[i] += scale_force * total_force;
            }
        }
    }
    else
    {
        //not in the feather region, so can calculate the forces
        //directly (also, no need to calculate shift, as 
        //the shifting function is constant, so does not
        //affect the gradient)
                
        if (group_pairs.isEmpty())
        {
            //no nb scale factors to worry about
            for (quint32 i=0; i<nats0; ++i)
            {
                const Parameter &param0 = params0_array[i];

                if (param0.reduced_charge == 0)
                    continue;
                    
                distmat.setOuterIndex(i);
                
                Vector total_force;
                
                for (quint32 j=0; j<nats1; ++j)
                {
                    if (params1_array[j].reduced_charge != 0)
                    {
                        //calculate the coulomb force
                        Vector cforce = -(param0.reduced_charge *
                                          params1_array[j].reduced_charge / 
                                          distmat[j].length2()) *
                                             
                                          distmat[j].direction();

                        total_force += cforce;
                    }
                }
                        
                group_forces0_array[i] += scale_force * total_force;
            }
        }
        else
        {
            //there are different nb scale factors between
            //different atoms...
            for (quint32 i=0; i<nats0; ++i)
            {
                const Parameter &param0 = params0_array[i];
                
                if (param0.reduced_charge == 0)    
                    continue;
                
                distmat.setOuterIndex(i);

                Vector total_force;
                
                for (quint32 j=0; j<nats1; ++j)
                {
                    const CLJScaleFactor &cljscl = group_pairs(i,j);
                          
                    double cscale = cljscl.coulomb() 
                                      * params1_array[j].reduced_charge;
                           
                    if (cscale != 0)
                    {
                        //calculate the coulomb force
                        Vector cforce = -(cscale *
                                          param0.reduced_charge /
                                          distmat[j].length2()) *
                                             
                                          distmat[j].direction();

                        total_force += cforce;
                    }
                }
                        
                group_forces0_array[i] += scale_force * total_force;

            } // end of loop over i atoms
        } // end of whether there are intra scale factors
    } // end of whether within feather region
}

void IntraCLJPotential::calculateCoulombForce(const CLJNBPairs::CGPairs &group_pairs,
                             const QSet<Index> &atoms0, const QSet<Index> &atoms1,
                             const CoordGroup &group0, const CoordGroup &group1,
                             const double mindist,
                             IntraCLJPotential::ForceWorkspace &distmat,
                             const IntraCLJPotential::Parameter *params0_array,
                             const IntraCLJPotential::Parameter *params1_array,
                             const double shift_coul,
                             Vector *group_forces0_array,
                             const double scale_force) const
{
    if (atoms0.isEmpty() or atoms1.isEmpty())
        return;

    if (mindist > switchfunc->featherDistance())
    {
        //we need to calculate the forces taking into account
        //the derivative of the switching function!
                
        //calculate the switching scale factors and their 
        //derivatives
        const double scl_coul = switchfunc->electrostaticScaleFactor( Length(mindist) );
                
        Vector group_sep = (group1.aaBox().center() -
                            group0.aaBox().center()).normalise(); 
                
        Vector dscl_coul = switchfunc->dElectrostaticScaleFactor( Length(mindist) ) 
                              * group_sep;

        if (group_pairs.isEmpty())
        {
            //there are no scale factors between atoms in these groups
            foreach (Index i, atoms0)
            {
                const Parameter &param0 = params0_array[i];

                if (param0.reduced_charge == 0)
                    return;

                distmat.setOuterIndex(i);

                Vector total_force;
                
                foreach (Index j, atoms1)
                {
                    const double q2 = param0.reduced_charge *
                                      params1_array[j].reduced_charge;
                                
                    if (q2 != 0)
                    {
                        //calculate the coulomb energy
                        const double cnrg = scl_coul * q2 /
                                            distmat[j].length();

                        //calculate the coulomb force
                        Vector cforce = (-cnrg / distmat[j].length() *
                                                 distmat[j].direction()) +
                                             
                                                ((cnrg-shift_coul) * dscl_coul);

                        total_force += cforce;
                    }
                }
                        
                group_forces0_array[i] += scale_force * total_force;
            }
        }
        else
        {
            //there are different nb scale factors between
            //the atoms. We need to calculate the forces
            //using them...
            foreach (Index i, atoms0)
            {
                const Parameter &param0 = params0_array[i];

                if (param0.reduced_charge == 0)
                    continue;
                
                distmat.setOuterIndex(i);

                Vector total_force;
                
                foreach (Index j, atoms1)
                {
                    const CLJScaleFactor &cljscl = group_pairs(i,j);
                            
                    if (cljscl.coulomb() != 0)
                    {
                        const double q2 = param0.reduced_charge * 
                                          params1_array[j].reduced_charge;
                                                      
                        if (q2 != 0)
                        {
                            //calculate the coulomb energy
                            const double cnrg = cljscl.coulomb() *
                                                  scl_coul * q2
                                                  / distmat[j].length();
                                               
                            //calculate the coulomb force
                            Vector cforce = (-cnrg 
                                             / distmat[j].length() *
                                               distmat[j].direction()) +
                                             
                                              ((cnrg-shift_coul) * dscl_coul);
                            
                            total_force += cforce;
                        }
                    }
                }
                        
                group_forces0_array[i] += scale_force * total_force;
            }
        }
    }
    else
    {
        //not in the feather region, so can calculate the forces
        //directly (also, no need to calculate shift, as 
        //the shifting function is constant, so does not
        //affect the gradient)
                
        if (group_pairs.isEmpty())
        {
            //no nb scale factors to worry about
            foreach (Index i, atoms0)
            {
                const Parameter &param0 = params0_array[i];

                if (param0.reduced_charge == 0)
                    continue;
                    
                distmat.setOuterIndex(i);
                
                Vector total_force;
                
                foreach (Index j, atoms1)
                {
                    if (params1_array[j].reduced_charge != 0)
                    {
                        //calculate the coulomb force
                        Vector cforce = -(param0.reduced_charge *
                                          params1_array[j].reduced_charge / 
                                          distmat[j].length2()) *
                                             
                                          distmat[j].direction();

                        total_force += cforce;
                    }
                }
                        
                group_forces0_array[i] += scale_force * total_force;
            }
        }
        else
        {
            //there are different nb scale factors between
            //different atoms...
            foreach (Index i, atoms0)
            {
                const Parameter &param0 = params0_array[i];
                
                if (param0.reduced_charge == 0)    
                    continue;
                
                distmat.setOuterIndex(i);

                Vector total_force;
                
                foreach (Index j, atoms1)
                {
                    const CLJScaleFactor &cljscl = group_pairs(i,j);
                          
                    double cscale = cljscl.coulomb() 
                                      * params1_array[j].reduced_charge;
                           
                    if (cscale != 0)
                    {
                        //calculate the coulomb force
                        Vector cforce = -(cscale *
                                          param0.reduced_charge /
                                          distmat[j].length2()) *
                                             
                                          distmat[j].direction();

                        total_force += cforce;
                    }
                }
                        
                group_forces0_array[i] += scale_force * total_force;

            } // end of loop over i atoms
        } // end of whether there are intra scale factors
    } // end of whether within feather region
}

/** Calculate the coulomb forces between the atoms in the molecule 'mol'
    and add these forces onto 'forces'. This uses
    the passed workspace to perform the calculation */
void IntraCLJPotential::calculateCoulombForce(
                                       const IntraCLJPotential::Molecule &mol,
                                       MolForceTable &forces, 
                                       IntraCLJPotential::ForceWorkspace &distmat,
                                       double scale_force) const
{
    if (scale_force == 0 or mol.isEmpty())
        return;
    
    const quint32 ngroups = mol.nCutGroups();
    
    const CoordGroup *groups_array = mol.coordinates().constData();
    
    const Parameters::Array *molparams_array 
                            = mol.parameters().atomicParameters().constData();
    
    BOOST_ASSERT(forces.nCutGroups() == mol.molecule().data().info().nCutGroups());
    BOOST_ASSERT(forces.molNum() == mol.molecule().data().number());
    
    const MolForceTable::Array *forces_array = forces.constData();

    const CLJNBPairs &nbpairs = mol.parameters().intraScaleFactors();
    
    //loop over all pairs of CutGroups in the molecule
    for (quint32 igroup=0; igroup<ngroups; ++igroup)
    {
        const Parameters::Array &params0 = molparams_array[igroup];

        const CoordGroup &group0 = groups_array[igroup];
        const AABox &aabox0 = group0.aaBox();
        const quint32 nats0 = group0.count();
        const Parameter *params0_array = params0.constData();
    
        //get the CGIdx of this CutGroup
        CGIdx cgidx_igroup = mol.cgIdx(igroup);
        
        //now get the index of the force table for this CutGroup
        int force_idx = forces.map(cgidx_igroup);
        
        if (force_idx == -1)
            //we are not interested in the forces on this CutGroup
            continue;
        
        //get the table that holds the forces acting on all of the
        //atoms of this CutGroup
        BOOST_ASSERT(forces_array[force_idx].count() == int(nats0));
        Vector *group_forces0_array = forces.data(force_idx);
    
        //get the forces acting on this group from all other groups
        //(yes, we are doing a full n2 loop, and not taking advantages
        // of equal and opposite forces)
        for (quint32 jgroup=0; jgroup<ngroups; ++jgroup)
        {
            const CoordGroup &group1 = groups_array[jgroup];
            const Parameters::Array &params1 = molparams_array[jgroup];

            //check first that these two CoordGroups could be within cutoff
            const bool outside_cutoff = igroup != jgroup and 
                                        spce->beyond(switchfunc->cutoffDistance(), 
                                                     aabox0, group1.aaBox());
            
            if (outside_cutoff)
                //this CutGroup is beyond the cutoff distance
                continue;
            
            //calculate all of the interatomic distances
            const double mindist = spce->calcDistVectors(group0, group1, distmat);
            
            if (mindist > switchfunc->cutoffDistance())
                //all of the atoms are definitely beyond cutoff
                continue;

            const quint32 nats1 = group1.count();
            
            //loop over all interatomic pairs and calculate the energies
            const Parameter *params1_array = params1.constData();

            CGIdx cgidx_jgroup = mol.cgIdx(jgroup);

            //get the non-bonded scale factors for all pairs of atoms
            //between these two groups (or within this group, if igroup == jgroup)
            const CLJNBPairs::CGPairs &group_pairs = nbpairs(cgidx_igroup,
                                                             cgidx_jgroup);

            double shift_coul = 0;

            if (use_electrostatic_shifting and igroup != jgroup)
               shift_coul = this->totalCharge(params0) * this->totalCharge(params1)
                           / switchfunc->electrostaticCutoffDistance();

            //calculate the forces acting on group0 caused by group1
            calculateCoulombForce(group_pairs, group0, group1,
                                  mindist, distmat, 
                                  params0_array, params1_array,
                                  nats0, nats1, shift_coul, 
                                  group_forces0_array, scale_force);

        } // end of loop over CutGroups (jgroup)

    } // end of loop over CutGroups (igroup)
}

/** Calculate the coulomb force acting on the part of the molecule
    in 'mol' caused by the rest of the molecule in 'rest_of_mol'. Note
    that these must both be of the same molecule, with the same
    layout UID and same nonbonded scale factors
    
    \throw SireError::incompatible_error
*/
void IntraCLJPotential::calculateCoulombForce(
                                        const IntraCLJPotential::Molecule &mol,
                                        const IntraCLJPotential::Molecule &rest_of_mol,
                                        MolForceTable &forces,
                                        IntraCLJPotential::ForceWorkspace &distmat,
                                        double scale_force) const
{
    if (scale_force == 0 or mol.isEmpty() or rest_of_mol.isEmpty())
        return;
        
    const quint32 ngroups0 = mol.nCutGroups();
    const CoordGroup *groups0_array = mol.coordinates().constData();
    
    const Parameters::Array *molparams0_array
                            = mol.parameters().atomicParameters().constData();
                            
    BOOST_ASSERT(forces.nCutGroups() == mol.molecule().data().info().nCutGroups());
    BOOST_ASSERT(forces.molNum() == mol.molecule().data().number());
    
    const MolForceTable::Array *forces_array = forces.constData();
    
    this->assertCompatible(mol, rest_of_mol);
    
    const quint32 ngroups1 = rest_of_mol.nCutGroups();
    const CoordGroup *groups1_array = rest_of_mol.coordinates().constData();
    
    const Parameters::Array *molparams1_array
                            = rest_of_mol.parameters().atomicParameters().constData();

    //this is now guaranteed to be the same in both passed parts of 
    //the molecule
    const CLJNBPairs &nbpairs = mol.parameters().intraScaleFactors();

    //calculate the forces acting on the atoms in 'mol' caused by 
    //the atoms in 'rest_of_mol'
    for (quint32 igroup=0; igroup<ngroups0; ++igroup)
    {
        const Parameters::Array &params0 = molparams0_array[igroup];

        const CoordGroup &group0 = groups0_array[igroup];
        const AABox &aabox0 = group0.aaBox();
        const quint32 nats0 = group0.count();
        const Parameter *params0_array = params0.constData();
    
        //get the CGIdx of this CutGroup
        CGIdx cgidx_igroup = mol.cgIdx(igroup);
        
        //now get the index of the force table for this CutGroup
        int force_idx = forces.map(cgidx_igroup);
        
        if (force_idx == -1)
            //we are not interested in the forces on this CutGroup
            continue;
        
        //get the table that holds the forces acting on all of the
        //atoms of this CutGroup
        BOOST_ASSERT(forces_array[force_idx].count() == int(nats0));
        Vector *group_forces0_array = forces.data(force_idx);
        
        for (quint32 jgroup=0; jgroup<ngroups1; ++jgroup)
        {
            CGIdx cgidx_jgroup = rest_of_mol.cgIdx(jgroup);

            if (mol.molecule().selection().selectedAll(cgidx_jgroup))
                //all of this CutGroup is part of 'mol', so ignore it
                continue;

            const CoordGroup &group1 = groups1_array[jgroup];
            const Parameters::Array &params1 = molparams1_array[jgroup];

            //check first that these two CoordGroups could be within cutoff
            const bool outside_cutoff = cgidx_igroup != cgidx_jgroup and 
                                        spce->beyond(switchfunc->cutoffDistance(), 
                                                     aabox0, group1.aaBox());
            
            if (outside_cutoff)
                //this CutGroup is beyond the cutoff distance
                continue;
            
            //calculate all of the interatomic distances
            const double mindist = spce->calcDistVectors(group0, group1, distmat);
            
            if (mindist > switchfunc->cutoffDistance())
                //all of the atoms are definitely beyond cutoff
                continue;

            const quint32 nats1 = group1.count();
            
            //loop over all interatomic pairs and calculate the energies
            const Parameter *params1_array = params1.constData();

            //get the non-bonded scale factors for all pairs of atoms
            //between these two groups (or within this group, if igroup == jgroup)
            const CLJNBPairs::CGPairs &group_pairs = nbpairs(cgidx_igroup,
                                                             cgidx_jgroup);
            
            if (cgidx_igroup == cgidx_jgroup or
                mol.molecule().selection().selected(cgidx_jgroup))
            {
                //some of the atoms in jgroup are selected as part of 'mol'.
                
                QSet<Index> atoms0 = mol.molecule().selection()
                                               .selectedAtoms(cgidx_igroup);
                                               
                QSet<Index> mol_atoms1 = atoms0;
                
                if (cgidx_igroup != cgidx_jgroup)
                    mol_atoms1 = mol.molecule().selection().selectedAtoms(cgidx_jgroup);
                                               
                QSet<Index> atoms1 = rest_of_mol.molecule().selection()
                                               .selectedAtoms(cgidx_jgroup);
                                         
                //remove the atoms in jgroup that are part of 'mol'
                atoms1 -= mol_atoms1;

                double shift_coul = 0;
            
                if (use_electrostatic_shifting and cgidx_igroup != cgidx_jgroup)
                    shift_coul = this->totalCharge(params0) * this->totalCharge(params1)
                                / switchfunc->electrostaticCutoffDistance();

                calculateCoulombForce(group_pairs, atoms0, atoms1,
                                      group0, group1,
                                      mindist, distmat,
                                      params0_array, params1_array,
                                      shift_coul,
                                      group_forces0_array, scale_force);
            }
            else
            {
                double shift_coul = 0;
            
                if (use_electrostatic_shifting)
                    shift_coul = this->totalCharge(params0) * this->totalCharge(params1)
                                / switchfunc->electrostaticCutoffDistance();

                calculateCoulombForce(group_pairs, group0, group1,
                                      mindist, distmat,
                                      params0_array, params1_array,
                                      nats0, nats1, shift_coul,
                                      group_forces0_array, scale_force);
            }
        }
    }
}

void IntraCLJPotential::calculateLJForce(const CLJNBPairs::CGPairs &group_pairs,
                             const CoordGroup &group0, const CoordGroup &group1,
                             const double mindist,
                             IntraCLJPotential::ForceWorkspace &distmat,
                             const IntraCLJPotential::Parameter *params0_array,
                             const IntraCLJPotential::Parameter *params1_array,
                             const quint32 nats0, const quint32 nats1,
                             Vector *group_forces0_array,
                             const double scale_force) const
{
    if (mindist > switchfunc->featherDistance())
    {
        //we need to calculate the forces taking into account
        //the derivative of the switching function!
        
        //calculate the switching scale factors and their 
        //derivatives
        const double scl_lj = switchfunc->vdwScaleFactor( Length(mindist) );
        
        Vector group_sep = (group1.aaBox().center() -
                            group0.aaBox().center()).normalise(); 
                             
        Vector dscl_lj = switchfunc->dVDWScaleFactor( Length(mindist) )
                             * group_sep;

        if (group_pairs.isEmpty())
        {
            //there are no scale factors between atoms in these groups
            for (quint32 i=0; i<nats0; ++i)
            {
                const Parameter &param0 = params0_array[i];

                Vector total_force;
        
                if (param0.ljid != 0)
                {
                    distmat.setOuterIndex(i);

                    for (quint32 j=0; j<nats1; ++j)
                    {
                        const Parameter &param1 = params1_array[j];
                      
                        if (param1.ljid != 0)
                        {
                            const double invdist = double(1) / distmat[j].length();
                
                            const LJPair &ljpair = ljpairs.constData()[
                                                   ljpairs.map(param0.ljid,
                                                               param1.ljid)];
                
                            double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                            double sig_over_dist12 = pow_2(sig_over_dist6);

                            //calculate the energy
                            const double ljnrg = 4 * ljpair.epsilon() *
                                              (sig_over_dist12 - sig_over_dist6);

                            // dU/dr requires an extra power of r
                            sig_over_dist6 *= invdist;
                            sig_over_dist12 *= invdist;

                            const Vector force = 
                                    ((scl_lj * 4 * ljpair.epsilon() * 
                                    (6.0*sig_over_dist6 - 12.0*sig_over_dist12))
                                    * distmat[j].direction())
                                    
                                    + (ljnrg * dscl_lj);

                            total_force += force;
                        }
                    }
                }
                
                group_forces0_array[i] += scale_force * total_force;
            }
        }
        else
        {
            //there are different nb scale factors between
            //the atoms. We need to calculate the forces
            //using them...
            for (quint32 i=0; i<nats0; ++i)
            {
                const Parameter &param0 = params0_array[i];
        
                Vector total_force;
        
                if (param0.ljid != 0)
                {
                    distmat.setOuterIndex(i);

                    for (quint32 j=0; j<nats1; ++j)
                    {
                        const Parameter &param1 = params1_array[j];
                        const CLJScaleFactor &cljscl = group_pairs(i,j);
                        
                        if (cljscl.lj() != 0 and param1.ljid != 0)
                        {
                            const double invdist = double(1) / distmat[j].length();
                            
                            const LJPair &ljpair = ljpairs.constData()[
                                                   ljpairs.map(param0.ljid,
                                                               param1.ljid)];
                
                            double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                            double sig_over_dist12 = pow_2(sig_over_dist6);

                            //calculate the energy
                            const double ljnrg = cljscl.lj() *
                                                 4 * ljpair.epsilon() *
                                          (sig_over_dist12 - sig_over_dist6);

                            // dU/dr requires an extra power of r
                            sig_over_dist6 *= invdist;
                            sig_over_dist12 *= invdist;

                            const Vector force = 
                                   ((scl_lj * 4 * ljpair.epsilon() * 
                                    (6.0*sig_over_dist6 - 12.0*sig_over_dist12))
                                   * distmat[j].direction())
                                   
                                   + (ljnrg * dscl_lj);

                            total_force += force;
                        }
                    }
                }
                
                group_forces0_array[i] += scale_force * total_force;
            }
        }
    }
    else
    {
        //not in the feather region, so can calculate the forces
        //directly
        if (group_pairs.isEmpty())
        {
            //no nb scale factors to worry about
            for (quint32 i=0; i<nats0; ++i)
            {
                distmat.setOuterIndex(i);
                const Parameter &param0 = params0_array[i];
        
                Vector total_force;
        
                if (param0.ljid != 0)
                {
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        const Parameter &param1 = params1_array[j];
                
                        if (param1.ljid != 0)
                        {
                            const double invdist = double(1) / distmat[j].length();
                            
                            const LJPair &ljpair = ljpairs.constData()[
                                                   ljpairs.map(param0.ljid,
                                                               param1.ljid)];
                
                            double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                            double sig_over_dist12 = pow_2(sig_over_dist6);

                            // dU/dr requires an extra power of r
                            sig_over_dist6 *= invdist;
                            sig_over_dist12 *= invdist;

                            const Vector force =
                                (4 * ljpair.epsilon() * (6.0*sig_over_dist6 - 
                                                        12.0*sig_over_dist12))
                                 * distmat[j].direction();

                            total_force += force;
                        }
                    }
                }
                
                group_forces0_array[i] += scale_force * total_force;
            }
        }
        else
        {
            //there are different nb scale factors between
            //different atoms...
            for (quint32 i=0; i<nats0; ++i)
            {
                const Parameter &param0 = params0_array[i];
        
                Vector total_force;
        
                if (param0.ljid != 0)
                {
                    distmat.setOuterIndex(i);

                    for (quint32 j=0; j<nats1; ++j)
                    {
                        const Parameter &param1 = params1_array[j];
                        const CLJScaleFactor &cljscl = group_pairs(i,j);
                        
                        if (cljscl.lj() != 0 and param1.ljid != 0)
                        {
                            const double invdist = double(1) / distmat[j].length();

                            const LJPair &ljpair = ljpairs.constData()[
                                                   ljpairs.map(param0.ljid,
                                                               param1.ljid)];
                
                            double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                            double sig_over_dist12 = pow_2(sig_over_dist6);

                            // dU/dr requires an extra power of r
                            sig_over_dist6 *= invdist;
                            sig_over_dist12 *= invdist;

                            const Vector force = 
                                  (cljscl.lj() *
                                   4 * ljpair.epsilon() * (6.0*sig_over_dist6 - 
                                                          12.0*sig_over_dist12))
                                   * distmat[j].direction();

                            total_force += force;
                        }
                    }
                }
                
                group_forces0_array[i] += scale_force * total_force;

            } // end of loop over i atoms

        } // end of whether there are intra scale factors

    } // end of whether within feather region
}

void IntraCLJPotential::calculateLJForce(const CLJNBPairs::CGPairs &group_pairs,
                             const QSet<Index> &atoms0, const QSet<Index> &atoms1,
                             const CoordGroup &group0, const CoordGroup &group1,
                             const double mindist,
                             IntraCLJPotential::ForceWorkspace &distmat,
                             const IntraCLJPotential::Parameter *params0_array,
                             const IntraCLJPotential::Parameter *params1_array,
                             Vector *group_forces0_array,
                             const double scale_force) const
{
    if (atoms0.isEmpty() or atoms1.isEmpty())
        return;

    if (mindist > switchfunc->featherDistance())
    {
        //we need to calculate the forces taking into account
        //the derivative of the switching function!
        
        //calculate the switching scale factors and their 
        //derivatives
        const double scl_lj = switchfunc->vdwScaleFactor( Length(mindist) );
        
        Vector group_sep = (group1.aaBox().center() -
                            group0.aaBox().center()).normalise(); 
                             
        Vector dscl_lj = switchfunc->dVDWScaleFactor( Length(mindist) )
                             * group_sep;

        if (group_pairs.isEmpty())
        {
            //there are no scale factors between atoms in these groups
            foreach (Index i, atoms0)
            {
                const Parameter &param0 = params0_array[i];

                Vector total_force;
        
                if (param0.ljid != 0)
                {
                    distmat.setOuterIndex(i);

                    foreach (Index j, atoms1)
                    {
                        const Parameter &param1 = params1_array[j];
                      
                        if (param1.ljid != 0)
                        {
                            const double invdist = double(1) / distmat[j].length();
                
                            const LJPair &ljpair = ljpairs.constData()[
                                                   ljpairs.map(param0.ljid,
                                                               param1.ljid)];
                
                            double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                            double sig_over_dist12 = pow_2(sig_over_dist6);

                            //calculate the energy
                            const double ljnrg = 4 * ljpair.epsilon() *
                                              (sig_over_dist12 - sig_over_dist6);

                            // dU/dr requires an extra power of r
                            sig_over_dist6 *= invdist;
                            sig_over_dist12 *= invdist;

                            const Vector force = 
                                    ((scl_lj * 4 * ljpair.epsilon() * 
                                    (6.0*sig_over_dist6 - 12.0*sig_over_dist12))
                                    * distmat[j].direction())
                                    
                                    + (ljnrg * dscl_lj);

                            total_force += force;
                        }
                    }
                }
                
                group_forces0_array[i] += scale_force * total_force;
            }
        }
        else
        {
            //there are different nb scale factors between
            //the atoms. We need to calculate the forces
            //using them...
            foreach (Index i, atoms0)
            {
                const Parameter &param0 = params0_array[i];
        
                Vector total_force;
        
                if (param0.ljid != 0)
                {
                    distmat.setOuterIndex(i);

                    foreach (Index j, atoms1)
                    {
                        const Parameter &param1 = params1_array[j];
                        const CLJScaleFactor &cljscl = group_pairs(i,j);
                        
                        if (cljscl.lj() != 0 and param1.ljid != 0)
                        {
                            const double invdist = double(1) / distmat[j].length();
                            
                            const LJPair &ljpair = ljpairs.constData()[
                                                   ljpairs.map(param0.ljid,
                                                               param1.ljid)];
                
                            double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                            double sig_over_dist12 = pow_2(sig_over_dist6);

                            //calculate the energy
                            const double ljnrg = cljscl.lj() *
                                                 4 * ljpair.epsilon() *
                                          (sig_over_dist12 - sig_over_dist6);

                            // dU/dr requires an extra power of r
                            sig_over_dist6 *= invdist;
                            sig_over_dist12 *= invdist;

                            const Vector force = 
                                   ((scl_lj * 4 * ljpair.epsilon() * 
                                    (6.0*sig_over_dist6 - 12.0*sig_over_dist12))
                                   * distmat[j].direction())
                                   
                                   + (ljnrg * dscl_lj);

                            total_force += force;
                        }
                    }
                }
                
                group_forces0_array[i] += scale_force * total_force;
            }
        }
    }
    else
    {
        //not in the feather region, so can calculate the forces
        //directly
        if (group_pairs.isEmpty())
        {
            //no nb scale factors to worry about
            foreach (Index i, atoms0)
            {
                distmat.setOuterIndex(i);
                const Parameter &param0 = params0_array[i];
        
                Vector total_force;
        
                if (param0.ljid != 0)
                {
                    foreach (Index j, atoms1)
                    {
                        const Parameter &param1 = params1_array[j];
                
                        if (param1.ljid != 0)
                        {
                            const double invdist = double(1) / distmat[j].length();
                            
                            const LJPair &ljpair = ljpairs.constData()[
                                                   ljpairs.map(param0.ljid,
                                                               param1.ljid)];
                
                            double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                            double sig_over_dist12 = pow_2(sig_over_dist6);

                            // dU/dr requires an extra power of r
                            sig_over_dist6 *= invdist;
                            sig_over_dist12 *= invdist;

                            const Vector force =
                                (4 * ljpair.epsilon() * (6.0*sig_over_dist6 - 
                                                        12.0*sig_over_dist12))
                                 * distmat[j].direction();

                            total_force += force;
                        }
                    }
                }
                
                group_forces0_array[i] += scale_force * total_force;
            }
        }
        else
        {
            //there are different nb scale factors between
            //different atoms...
            foreach (Index i, atoms0)
            {
                const Parameter &param0 = params0_array[i];
        
                Vector total_force;
        
                if (param0.ljid != 0)
                {
                    distmat.setOuterIndex(i);

                    foreach (Index j, atoms1)
                    {
                        const Parameter &param1 = params1_array[j];
                        const CLJScaleFactor &cljscl = group_pairs(i,j);
                        
                        if (cljscl.lj() != 0 and param1.ljid != 0)
                        {
                            const double invdist = double(1) / distmat[j].length();

                            const LJPair &ljpair = ljpairs.constData()[
                                                   ljpairs.map(param0.ljid,
                                                               param1.ljid)];
                
                            double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                            double sig_over_dist12 = pow_2(sig_over_dist6);

                            // dU/dr requires an extra power of r
                            sig_over_dist6 *= invdist;
                            sig_over_dist12 *= invdist;

                            const Vector force = 
                                  (cljscl.lj() *
                                   4 * ljpair.epsilon() * (6.0*sig_over_dist6 - 
                                                          12.0*sig_over_dist12))
                                   * distmat[j].direction();

                            total_force += force;
                        }
                    }
                }
                
                group_forces0_array[i] += scale_force * total_force;

            } // end of loop over i atoms

        } // end of whether there are intra scale factors

    } // end of whether within feather region
}

/** Calculate the LJ forces between the atoms in the molecule 'mol'
    and add these forces onto 'forces'. This uses
    the passed workspace to perform the calculation */
void IntraCLJPotential::calculateLJForce(const IntraCLJPotential::Molecule &mol,
                                         MolForceTable &forces, 
                                         IntraCLJPotential::ForceWorkspace &distmat,
                                         double scale_force) const
{
    if (scale_force == 0 or mol.isEmpty())
        return;
    
    const quint32 ngroups = mol.nCutGroups();
    
    const CoordGroup *groups_array = mol.coordinates().constData();
    
    const Parameters::Array *molparams_array 
                            = mol.parameters().atomicParameters().constData();
    
    BOOST_ASSERT(forces.nCutGroups() == mol.molecule().data().info().nCutGroups());
    BOOST_ASSERT(forces.molNum() == mol.molecule().data().number());
    
    const MolForceTable::Array *forces_array = forces.constData();

    const CLJNBPairs &nbpairs = mol.parameters().intraScaleFactors();
    
    //loop over all pairs of CutGroups in the molecule
    for (quint32 igroup=0; igroup<ngroups; ++igroup)
    {
        const Parameters::Array &params0 = molparams_array[igroup];

        const CoordGroup &group0 = groups_array[igroup];
        const AABox &aabox0 = group0.aaBox();
        const quint32 nats0 = group0.count();
        const Parameter *params0_array = params0.constData();
    
        //get the CGIdx of this CutGroup
        CGIdx cgidx_igroup = mol.cgIdx(igroup);
        
        //now get the index of the force table for this CutGroup
        int force_idx = forces.map(cgidx_igroup);
        
        if (force_idx == -1)
            //we are not interested in the forces on this CutGroup
            continue;
        
        //get the table that holds the forces acting on all of the
        //atoms of this CutGroup
        BOOST_ASSERT(forces_array[force_idx].count() == int(nats0));
        Vector *group_forces0_array = forces.data(force_idx);
    
        //get the forces acting on this group from all other groups
        //(yes, we are doing a full n2 loop, and not taking advantages
        // of equal and opposite forces)
        for (quint32 jgroup=0; jgroup<ngroups; ++jgroup)
        {
            const CoordGroup &group1 = groups_array[jgroup];
            const Parameters::Array &params1 = molparams_array[jgroup];

            //check first that these two CoordGroups could be within cutoff
            const bool outside_cutoff = igroup != jgroup and 
                                        spce->beyond(switchfunc->cutoffDistance(), 
                                                     aabox0, group1.aaBox());
            
            if (outside_cutoff)
                //this CutGroup is beyond the cutoff distance
                continue;
            
            //calculate all of the interatomic distances
            const double mindist = spce->calcDistVectors(group0, group1, distmat);
            
            if (mindist > switchfunc->cutoffDistance())
                //all of the atoms are definitely beyond cutoff
                continue;

            const quint32 nats1 = group1.count();
            
            //loop over all interatomic pairs and calculate the energies
            const Parameter *params1_array = params1.constData();

            CGIdx cgidx_jgroup = mol.cgIdx(jgroup);

            //get the non-bonded scale factors for all pairs of atoms
            //between these two groups (or within this group, if igroup == jgroup)
            const CLJNBPairs::CGPairs &group_pairs = nbpairs(cgidx_igroup,
                                                             cgidx_jgroup);

            //calculate the forces acting on group0 caused by group1
            calculateLJForce(group_pairs, group0, group1,
                             mindist, distmat, 
                             params0_array, params1_array,
                             nats0, nats1, group_forces0_array, scale_force);
            
        } // end of loop over CutGroups (jgroup)

    } // end of loop over CutGroups (igroup)
}

/** Calculate the LJ force acting on the part of the molecule
    in 'mol' caused by the rest of the molecule in 'rest_of_mol'. Note
    that these must both be of the same molecule, with the same
    layout UID and same nonbonded scale factors
    
    \throw SireError::incompatible_error
*/
void IntraCLJPotential::calculateLJForce(
                                        const IntraCLJPotential::Molecule &mol,
                                        const IntraCLJPotential::Molecule &rest_of_mol,
                                        MolForceTable &forces,
                                        IntraCLJPotential::ForceWorkspace &distmat,
                                        double scale_force) const
{
    if (scale_force == 0 or mol.isEmpty() or rest_of_mol.isEmpty())
        return;
        
    const quint32 ngroups0 = mol.nCutGroups();
    const CoordGroup *groups0_array = mol.coordinates().constData();
    
    const Parameters::Array *molparams0_array
                            = mol.parameters().atomicParameters().constData();
                            
    BOOST_ASSERT(forces.nCutGroups() == mol.molecule().data().info().nCutGroups());
    BOOST_ASSERT(forces.molNum() == mol.molecule().data().number());
    
    const MolForceTable::Array *forces_array = forces.constData();
    
    this->assertCompatible(mol, rest_of_mol);
    
    const quint32 ngroups1 = rest_of_mol.nCutGroups();
    const CoordGroup *groups1_array = rest_of_mol.coordinates().constData();
    
    const Parameters::Array *molparams1_array
                            = rest_of_mol.parameters().atomicParameters().constData();

    //this is now guaranteed to be the same in both passed parts of 
    //the molecule
    const CLJNBPairs &nbpairs = mol.parameters().intraScaleFactors();

    //calculate the forces acting on the atoms in 'mol' caused by 
    //the atoms in 'rest_of_mol'
    for (quint32 igroup=0; igroup<ngroups0; ++igroup)
    {
        const Parameters::Array &params0 = molparams0_array[igroup];

        const CoordGroup &group0 = groups0_array[igroup];
        const AABox &aabox0 = group0.aaBox();
        const quint32 nats0 = group0.count();
        const Parameter *params0_array = params0.constData();
    
        //get the CGIdx of this CutGroup
        CGIdx cgidx_igroup = mol.cgIdx(igroup);
        
        //now get the index of the force table for this CutGroup
        int force_idx = forces.map(cgidx_igroup);
        
        if (force_idx == -1)
            //we are not interested in the forces on this CutGroup
            continue;
        
        //get the table that holds the forces acting on all of the
        //atoms of this CutGroup
        BOOST_ASSERT(forces_array[force_idx].count() == int(nats0));
        Vector *group_forces0_array = forces.data(force_idx);
        
        for (quint32 jgroup=0; jgroup<ngroups1; ++jgroup)
        {
            CGIdx cgidx_jgroup = rest_of_mol.cgIdx(jgroup);

            if (mol.molecule().selection().selectedAll(cgidx_jgroup))
                continue;

            const CoordGroup &group1 = groups1_array[jgroup];
            const Parameters::Array &params1 = molparams1_array[jgroup];

            //check first that these two CoordGroups could be within cutoff
            const bool outside_cutoff = cgidx_igroup != cgidx_jgroup and 
                                        spce->beyond(switchfunc->cutoffDistance(), 
                                                     aabox0, group1.aaBox());
            
            if (outside_cutoff)
                //this CutGroup is beyond the cutoff distance
                continue;
            
            //calculate all of the interatomic distances
            const double mindist = spce->calcDistVectors(group0, group1, distmat);
            
            if (mindist > switchfunc->cutoffDistance())
                //all of the atoms are definitely beyond cutoff
                continue;

            const quint32 nats1 = group1.count();
            
            //loop over all interatomic pairs and calculate the energies
            const Parameter *params1_array = params1.constData();

            //get the non-bonded scale factors for all pairs of atoms
            //between these two groups (or within this group, if igroup == jgroup)
            const CLJNBPairs::CGPairs &group_pairs = nbpairs(cgidx_igroup,
                                                             cgidx_jgroup);
            
            if (cgidx_igroup == cgidx_jgroup or
                mol.molecule().selection().selected(cgidx_jgroup))
            {
                //some of the atoms in jgroup are also selected in 'mol'
                
                QSet<Index> atoms0 = mol.molecule().selection()
                                               .selectedAtoms(cgidx_igroup);
                                               
                QSet<Index> mol_atoms1 = atoms0;
                
                if (cgidx_jgroup != cgidx_igroup)
                    mol_atoms1 = mol.molecule().selection().selectedAtoms(cgidx_jgroup);
                                               
                QSet<Index> atoms1 = rest_of_mol.molecule().selection()
                                               .selectedAtoms(cgidx_jgroup);
                                         
                //remove the atoms in 'rest_of_mol' that are part of 'mol'
                atoms1 -= mol_atoms1;
                                                           
                calculateLJForce(group_pairs, atoms0, atoms1,
                                 group0, group1,
                                 mindist, distmat,
                                 params0_array, params1_array,
                                 group_forces0_array, scale_force);
            }
            else
            {
                calculateLJForce(group_pairs, group0, group1,
                                 mindist, distmat,
                                 params0_array, params1_array,
                                 nats0, nats1, group_forces0_array, scale_force);
            }
        }
    }
}

void IntraCLJPotential::calculateField(const IntraCLJPotential::Molecule &mol, 
                    const CLJProbe &probe,
                    MolFieldTable &fields,
                    IntraCLJPotential::FieldWorkspace &workspace,
                    double scale_field) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code necessary to calculate intramolecular coulomb "
                "and LJ fields has not yet been written..."), CODELOC );
}

void IntraCLJPotential::calculateField(const IntraCLJPotential::Molecule &mol,
                    const IntraCLJPotential::Molecule &rest_of_mol,
                    const CLJProbe &probe,
                    MolFieldTable &fields,
                    IntraCLJPotential::FieldWorkspace &workspace,
                    double scale_field) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code necessary to calculate intramolecular coulomb "
                "and LJ fields has not yet been written..."), CODELOC );
}

void IntraCLJPotential::calculateField(const IntraCLJPotential::Molecule &mol, 
                    const CLJProbe &probe,
                    MolFieldTable &fields,
                    const Symbol &symbol,
                    const Components &components,
                    IntraCLJPotential::FieldWorkspace &workspace,
                    double scale_field) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code necessary to calculate intramolecular coulomb "
                "and LJ fields has not yet been written..."), CODELOC );
}

void IntraCLJPotential::calculateField(const IntraCLJPotential::Molecule &mol,
                    const IntraCLJPotential::Molecule &rest_of_mol,
                    const CLJProbe &probe,
                    MolFieldTable &fields,
                    const Symbol &symbol,
                    const Components &components,
                    IntraCLJPotential::FieldWorkspace &workspace,
                    double scale_field) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code necessary to calculate intramolecular coulomb "
                "and LJ fields has not yet been written..."), CODELOC );
}

void IntraCLJPotential::calculateField(const IntraCLJPotential::Molecule &mol, 
                    const CLJProbe &probe,
                    GridFieldTable &fields,
                    IntraCLJPotential::FieldWorkspace &workspace,
                    double scale_field) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code necessary to calculate intramolecular coulomb "
                "and LJ fields has not yet been written..."), CODELOC );
}

void IntraCLJPotential::calculateField(const IntraCLJPotential::Molecule &mol, 
                    const CLJProbe &probe,
                    GridFieldTable &fields,
                    const Symbol &symbol,
                    const Components &components,
                    IntraCLJPotential::FieldWorkspace &workspace,
                    double scale_field) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code necessary to calculate intramolecular coulomb "
                "and LJ fields has not yet been written..."), CODELOC );
}

void IntraCLJPotential::calculatePotential(const IntraCLJPotential::Molecule &mol, 
                        const CLJProbe &probe,
                        MolPotentialTable &potentials,
                        IntraCLJPotential::PotentialWorkspace &workspace,
                        double scale_potential) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code necessary to calculate intramolecular coulomb "
                "and LJ potentials has not yet been written..."), CODELOC );
}

void IntraCLJPotential::calculatePotential(const IntraCLJPotential::Molecule &mol,
                        const IntraCLJPotential::Molecule &rest_of_mol,
                        const CLJProbe &probe,
                        MolPotentialTable &potentials,
                        IntraCLJPotential::PotentialWorkspace &workspace,
                        double scale_potential) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code necessary to calculate intramolecular coulomb "
                "and LJ fields has not yet been written..."), CODELOC );
}

void IntraCLJPotential::calculatePotential(const IntraCLJPotential::Molecule &mol, 
                        const CLJProbe &probe,
                        MolPotentialTable &potentials,
                        const Symbol &symbol,
                        const Components &components,
                        IntraCLJPotential::PotentialWorkspace &workspace,
                        double scale_potential) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code necessary to calculate intramolecular coulomb "
                "and LJ potentials has not yet been written..."), CODELOC );
}

void IntraCLJPotential::calculatePotential(const IntraCLJPotential::Molecule &mol,
                        const IntraCLJPotential::Molecule &rest_of_mol,
                        const CLJProbe &probe,
                        MolPotentialTable &potentials,
                        const Symbol &symbol,
                        const Components &components,
                        IntraCLJPotential::PotentialWorkspace &workspace,
                        double scale_potential) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code necessary to calculate intramolecular coulomb "
                "and LJ potentials has not yet been written..."), CODELOC );
}

void IntraCLJPotential::calculatePotential(const IntraCLJPotential::Molecule &mol, 
                        const CLJProbe &probe,
                        GridPotentialTable &potentials,
                        IntraCLJPotential::PotentialWorkspace &workspace,
                        double scale_potential) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code necessary to calculate intramolecular coulomb "
                "and LJ potentials has not yet been written..."), CODELOC );
}

void IntraCLJPotential::calculatePotential(const IntraCLJPotential::Molecule &mol, 
                        const CLJProbe &probe,
                        GridPotentialTable &potentials,
                        const Symbol &symbol,
                        const Components &components,
                        IntraCLJPotential::PotentialWorkspace &workspace,
                        double scale_potential) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code necessary to calculate intramolecular coulomb "
                "and LJ potentials has not yet been written..."), CODELOC );
}

void IntraCLJPotential::calculateCoulombField(const IntraCLJPotential::Molecule &mol,
                           const CLJProbe &probe,
                           MolFieldTable &fields,
                           IntraCLJPotential::FieldWorkspace &workspace,
                           double scale_field) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code necessary to calculate intramolecular coulomb "
                "and LJ fields has not yet been written..."), CODELOC );
}

void IntraCLJPotential::calculateCoulombField(const IntraCLJPotential::Molecule &mol,
                           const CLJProbe &probe,
                           GridFieldTable &fields,
                           IntraCLJPotential::FieldWorkspace &workspace,
                           double scale_field) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code necessary to calculate intramolecular coulomb "
                "and LJ fields has not yet been written..."), CODELOC );
}

void IntraCLJPotential::calculateCoulombField(const IntraCLJPotential::Molecule &mol,
                           const IntraCLJPotential::Molecule &rest_of_mol,
                           const CLJProbe &probe,
                           MolFieldTable &fields,
                           IntraCLJPotential::FieldWorkspace &workspace,
                           double scale_field) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code necessary to calculate intramolecular coulomb "
                "and LJ fields has not yet been written..."), CODELOC );
}

void IntraCLJPotential::calculateLJField(const IntraCLJPotential::Molecule &mol,
                      const CLJProbe &probe,
                      MolFieldTable &fields,
                      IntraCLJPotential::FieldWorkspace &workspace,
                      double scale_field) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code necessary to calculate intramolecular coulomb "
                "and LJ fields has not yet been written..."), CODELOC );
}

void IntraCLJPotential::calculateLJField(const IntraCLJPotential::Molecule &mol,
                      const CLJProbe &probe,
                      GridFieldTable &fields,
                      IntraCLJPotential::FieldWorkspace &workspace,
                      double scale_field) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code necessary to calculate intramolecular coulomb "
                "and LJ fields has not yet been written..."), CODELOC );
}

void IntraCLJPotential::calculateLJField(const IntraCLJPotential::Molecule &mol,
                      const IntraCLJPotential::Molecule &rest_of_mol,
                      const CLJProbe &probe,
                      MolFieldTable &fields,
                      IntraCLJPotential::FieldWorkspace &workspace,
                      double scale_field) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code necessary to calculate intramolecular coulomb "
                "and LJ fields has not yet been written..."), CODELOC );
}

void IntraCLJPotential::calculateCoulombPotential(const IntraCLJPotential::Molecule &mol,
                               const CLJProbe &probe,
                               MolPotentialTable &potentials,
                               IntraCLJPotential::PotentialWorkspace &workspace,
                               double scale_potential) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code necessary to calculate intramolecular coulomb "
                "and LJ potentials has not yet been written..."), CODELOC );
}

void IntraCLJPotential::calculateCoulombPotential(const IntraCLJPotential::Molecule &mol,
                               const CLJProbe &probe,
                               GridPotentialTable &fields,
                               IntraCLJPotential::PotentialWorkspace &workspace,
                               double scale_potential) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code necessary to calculate intramolecular coulomb "
                "and LJ potentials has not yet been written..."), CODELOC );
}

void IntraCLJPotential::calculateCoulombPotential(const IntraCLJPotential::Molecule &mol,
                               const IntraCLJPotential::Molecule &rest_of_mol,
                               const CLJProbe &probe,
                               MolPotentialTable &potentials,
                               IntraCLJPotential::PotentialWorkspace &workspace,
                               double scale_potential) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code necessary to calculate intramolecular coulomb "
                "and LJ potentials has not yet been written..."), CODELOC );
}

void IntraCLJPotential::calculateLJPotential(const IntraCLJPotential::Molecule &mol,
                          const CLJProbe &probe,
                          MolPotentialTable &potentials,
                          IntraCLJPotential::PotentialWorkspace &workspace,
                          double scale_potential) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code necessary to calculate intramolecular coulomb "
                "and LJ potentials has not yet been written..."), CODELOC );
}

void IntraCLJPotential::calculateLJPotential(const IntraCLJPotential::Molecule &mol,
                          const CLJProbe &probe,
                          GridPotentialTable &potentials,
                          IntraCLJPotential::PotentialWorkspace &workspace,
                          double scale_potential) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code necessary to calculate intramolecular coulomb "
                "and LJ potentials has not yet been written..."), CODELOC );
}

void IntraCLJPotential::calculateLJPotential(const IntraCLJPotential::Molecule &mol,
                          const IntraCLJPotential::Molecule &rest_of_mol,
                          const CLJProbe &probe,
                          MolPotentialTable &potentials,
                          IntraCLJPotential::PotentialWorkspace &workspace,
                          double scale_potential) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code necessary to calculate intramolecular coulomb "
                "and LJ potentials has not yet been written..."), CODELOC );
}
