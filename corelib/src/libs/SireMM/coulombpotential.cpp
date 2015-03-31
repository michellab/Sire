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

#include "coulombpotential.h"
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

using namespace SireMol;
using namespace SireVol;

using namespace SireMaths;

using namespace SireBase;

using namespace SireStream;

///////
/////// Static data for ChargeParameterName
///////

QString ChargeParameterName::chg_param( "charge" );

///////
/////// Completely instantiate the CoulombPotential ancillary classes
///////

namespace SireMM
{
    namespace detail
    {
        template
        class IntraScaledParameters<CoulombNBPairs>;

        template
        class IntraScaledAtomicParameters< AtomicParameters3D<ChargeParameter>,
                                           IntraScaledParameters<CoulombNBPairs> >;
    }
}

namespace SireFF
{
    namespace detail
    {
        template
        class AtomicParameters3D<ChargeParameter>;

        template
        class FFMolecule3D<InterCoulombPotential>;

        template
        class FFMolecules3D<InterCoulombPotential>;

        template
        class ChangedMolecule<InterCoulombPotential::Molecule>;

        template
        class FFMolecule3D<IntraCoulombPotential>;

        template
        class FFMolecules3D<IntraCoulombPotential>;

        template
        class ChangedMolecule<IntraCoulombPotential::Molecule>;
    }
}

/** Streaming functions for ChargeParameter */
QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, 
                                      const SireMM::detail::ChargeParameter &chgparam)
{
    ds << chgparam.reduced_charge;

    return ds;
}

QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, 
                                      SireMM::detail::ChargeParameter &chgparam)
{
    ds >> chgparam.reduced_charge;
    
    return ds;
}

/** Internal function used to get the charge parameters from a molecule
    and convert them into a PackedArray of reduced charge */
static PackedArray2D<ChargeParameter> 
getChargeParameters(const PartialMolecule &molecule,
                    const PropertyName &charge_property)
{
    const AtomCharges &chgs = molecule.property(charge_property).asA<AtomCharges>();

    const AtomSelection &selected_atoms = molecule.selection();
    
    if (selected_atoms.selectedNone())
        return PackedArray2D<ChargeParameter>();
    
    //create space for the parameters - only need space for CutGroups
    //that contain at least one selected atom
    QVector< QVector<ChargeParameter> > chgparams( selected_atoms.nSelectedCutGroups() );
    QVector<ChargeParameter>* chgparams_array = chgparams.data();

    const double sqrt_4pieps0 = std::sqrt(SireUnits::one_over_four_pi_eps0);

    if (selected_atoms.selectedAllCutGroups())
    {
        const int ncg = molecule.data().info().nCutGroups();
    
        for (CGIdx i(0); i<ncg; ++i)
        {
            const int nats = molecule.data().info().nAtoms(i);
            
            QVector<ChargeParameter> group_charges(nats);
            ChargeParameter *group_charges_array = group_charges.data();
            
            //get the arrays containing the charge parameters
            //for this CutGroup
            const SireUnits::Dimension::Charge *group_chgs = chgs.constData(i);
            
            if (selected_atoms.selectedAll(i))
            {
                for (Index j(0); j<nats; ++j)
                {
                    group_charges_array[j].reduced_charge = group_chgs[j] * sqrt_4pieps0;
                }
            }
            else
            {
                foreach (Index j, selected_atoms.selectedAtoms(i))
                {
                    group_charges_array[j].reduced_charge = group_chgs[j] * sqrt_4pieps0;
                }
            }
            
            chgparams_array[i] = group_charges;
        }
    }
    else
    {
        foreach (CGIdx i, selected_atoms.selectedCutGroups())
        {
            const int nats = molecule.data().info().nAtoms(i);
            
            QVector<ChargeParameter> group_charges(nats);
            ChargeParameter *group_charges_array = group_charges.data();
            
            //get the arrays containing the charge parameters
            //for this CutGroup
            const SireUnits::Dimension::Charge *group_chgs = chgs.constData(i);
            
            if (selected_atoms.selectedAll(i))
            {
                for (Index j(0); j<nats; ++j)
                {
                    group_charges_array[j].reduced_charge = group_chgs[j] * sqrt_4pieps0;
                }
            }
            else
            {
                foreach (Index j, selected_atoms.selectedAtoms(i))
                {
                    group_charges_array[j].reduced_charge = group_chgs[j] * sqrt_4pieps0;
                }
            }
            
            chgparams_array[i] = group_charges;
        }
    }
    
    return PackedArray2D<ChargeParameter>( chgparams );
}

/////////////
///////////// Implementation of CoulombPotential
/////////////

static const RegisterMetaType<CoulombPotential> r_coulpot( MAGIC_ONLY, NO_ROOT,
                                                           "SireMM::CoulombPotential" );

/** Serialise to a binary datastream */
QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds,
                                      const CoulombPotential &coulpot)
{
    writeHeader(ds, r_coulpot, 1);
    
    ds << coulpot.props;

    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds,
                                      CoulombPotential &coulpot)
{
    VersionID v = readHeader(ds, r_coulpot);
    
    if (v == 1)
    {
        ds >> coulpot.props;
    
        //extract all of the properties
        coulpot.spce = coulpot.props.property("space").asA<Space>();
        coulpot.switchfunc = coulpot.props.property("switchingFunction")
                                          .asA<SwitchingFunction>();
    
        coulpot.use_electrostatic_shifting = coulpot.props.property("shiftElectrostatics")
                                        .asA<VariantProperty>().convertTo<bool>();
    }
    else 
        throw version_error(v, "1", r_coulpot, CODELOC);
    
    return ds;
}

/** Constructor */
CoulombPotential::CoulombPotential()
                 : spce( Space::null() ), switchfunc( SwitchingFunction::null() ),
                   use_electrostatic_shifting(false)
{
    //record the defaults
    props.setProperty( "space", spce );
    props.setProperty( "switchingFunction", switchfunc );
    props.setProperty( "shiftElectrostatics",
                       VariantProperty(use_electrostatic_shifting) );
}

/** Copy constructor */
CoulombPotential::CoulombPotential(const CoulombPotential &other)
                 : props(other.props),
                   spce(other.spce), switchfunc(other.switchfunc),
                   use_electrostatic_shifting(other.use_electrostatic_shifting)
{}

/** Destructor */
CoulombPotential::~CoulombPotential()
{}

/** Copy assignment operator */
CoulombPotential& CoulombPotential::operator=(const CoulombPotential &other)
{
    if (this != &other)
    {
        props = other.props;
        spce = other.spce;
        switchfunc = other.switchfunc;
        use_electrostatic_shifting = other.use_electrostatic_shifting;
    }
    
    return *this;
}

/** You need to call this function before you start a block of 
    energy of force evaluation using this forcefield. You should
    also call 'finishedEvaluation()' once you have finished. */
void CoulombPotential::startEvaluation()
{}

/** You should call this function once you have finished a block of
    force or energy evaluation using this potential */
void CoulombPotential::finishedEvaluation()
{}

/** Return all of the properties set in this forcefield */
const Properties& CoulombPotential::properties() const
{
    return props;
}

/** Return whether or not this potential has a property called 'name' */
bool CoulombPotential::containsProperty(const QString &name) const
{
    return props.hasProperty(name);
}

/** Return the property with name 'name'

    \throw SireBase::missing_property
*/
const Property& CoulombPotential::property(const QString &name) const
{
    return props.property(name);
}

/** Set the 3D space in which the molecules in this potential are evaluated */
bool CoulombPotential::setSpace(const Space &new_space)
{
    if (not spce->equals(new_space))
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
bool CoulombPotential::setSwitchingFunction(const SwitchingFunction &new_switchfunc)
{
    if (not switchfunc->equals(new_switchfunc))
    {
        switchfunc = new_switchfunc;
        props.setProperty( "switchingFunction", new_switchfunc );
        
        this->changedPotential();
        return true;
    }
    else
        return false;
}

/** Set whether or not to shift the electrostatics between CutGroups so that
    the group-group net electrostatic interaction energy between CutGroups
    is zero at the cutoff */
bool CoulombPotential::setShiftElectrostatics(bool switchelectro)
{
    if (switchelectro != use_electrostatic_shifting)
    {
        use_electrostatic_shifting = switchelectro;
        props.setProperty( "shiftElectrostatics", VariantProperty(switchelectro) );
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
bool CoulombPotential::setProperty(const QString &name, const Property &value)
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
    else
        throw SireBase::missing_property( QObject::tr(
            "The CLJ potentials do not have a property called \"%1\" that "
            "can be changed. Available properties are [ %2 ].")
                .arg(name, QStringList(props.propertyKeys()).join(", ")), CODELOC );
                
    return false;
}

/** Return the 3D space in which this potential is evaluated */
const Space& CoulombPotential::space() const
{
    return *spce;
}

/** Return the switching function used to scale the group-group
    interactions to zero at the cutoff */
const SwitchingFunction& CoulombPotential::switchingFunction() const
{
    return *switchfunc;
}

/** Return whether or not the net group-group electrostatic interaction
    energy is shifted so that it is zero at the cutoff */
bool CoulombPotential::shiftElectrostatics() const
{
    return use_electrostatic_shifting;
}

/////////////
///////////// Implementation of InterCoulombPotential
/////////////

static const RegisterMetaType<InterCoulombPotential> r_intercoul( MAGIC_ONLY, NO_ROOT,
                                            InterCoulombPotential::typeName() );

/** Serialise to a binary datastream */
QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds,
                                      const InterCoulombPotential &intercoul)
{
    writeHeader(ds, r_intercoul, 1);
    
    ds << static_cast<const CoulombPotential&>(intercoul);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds,
                                      InterCoulombPotential &intercoul)
{
    VersionID v = readHeader(ds, r_intercoul);
    
    if (v == 1)
    {
        ds >> static_cast<CoulombPotential&>(intercoul);
    }
    else
        throw version_error(v, "1", r_intercoul, CODELOC);
        
    return ds;
}

/** Constructor */
InterCoulombPotential::InterCoulombPotential() : CoulombPotential()
{}

/** Copy constructor */
InterCoulombPotential::InterCoulombPotential(const InterCoulombPotential &other)
                      : CoulombPotential(other)
{}

/** Destructor */
InterCoulombPotential::~InterCoulombPotential()
{}

/** Copy assignment operator */
InterCoulombPotential& 
InterCoulombPotential::operator=(const InterCoulombPotential &other)
{
    CoulombPotential::operator=(other);
    return *this;
}

void InterCoulombPotential::throwMissingForceComponent(const Symbol &symbol,
                              const IntraCoulombPotential::Components &components) const
{
    throw SireFF::missing_component( QObject::tr(
        "There is no force component in potential %1 - available "
        "components are %2.")
            .arg(this->what())
            .arg(components.total().toString()), CODELOC );
}

void InterCoulombPotential::throwMissingFieldComponent(const Symbol &symbol,
                              const IntraCoulombPotential::Components &components) const
{
    throw SireFF::missing_component( QObject::tr(
        "There is no field component in potential %1 - available "
        "components are %2.")
            .arg(this->what())
            .arg(components.total().toString()), CODELOC );
}

void InterCoulombPotential::throwMissingPotentialComponent(const Symbol &symbol,
                              const IntraCoulombPotential::Components &components) const
{
    throw SireFF::missing_component( QObject::tr(
        "There is no potential component in potential %1 - available "
        "components are %2.")
            .arg(this->what())
            .arg(components.total().toString()), CODELOC );
}

/** Return all of the parameters needed by this potential for 
    the molecule 'molecule', using the supplied property map to
    find the properties that contain those parameters
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
InterCoulombPotential::Parameters 
InterCoulombPotential::getParameters(const PartialMolecule &molecule,
                                     const PropertyMap &map)
{
    return Parameters( molecule, map[parameters().coordinates()],
                       getChargeParameters(molecule, map[parameters().charge()]) );
}

/** Update the parameters for the molecule going from 'old_molecule' to 
    'new_molecule', with the parameters found using the property map 'map'
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
InterCoulombPotential::Parameters
InterCoulombPotential::updateParameters(const InterCoulombPotential::Parameters &old_params,
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
    
    //get what has changed
    bool new_coords = old_molecule.version(coords_property) !=
                         new_molecule.version(coords_property);
                             
    bool new_chg = ( old_molecule.version(chg_property) !=
                         new_molecule.version(chg_property) );

    if (new_coords)
    {
        new_params.setAtomicCoordinates( AtomicCoords3D(new_molecule, 
                                                        coords_property) );
    }

    if (new_chg)
    {
        new_params.setAtomicParameters( getChargeParameters(new_molecule,chg_property) );
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
InterCoulombPotential::Parameters
InterCoulombPotential::updateParameters(const InterCoulombPotential::Parameters &old_params,
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
    
    const PropertyName &new_coords = new_map[parameters().coordinates()];
    const PropertyName &new_chg = new_map[parameters().charge()];
    
    //get what has changed
    bool changed_coords = (new_coords != old_coords) or
                           old_molecule.version(old_coords) !=
                           new_molecule.version(old_coords);
                             
    bool changed_chg = (new_chg != old_chg) or
                       ( old_molecule.version(old_chg) !=
                         new_molecule.version(old_chg) );

    if (changed_coords)
        new_params.setAtomicCoordinates( AtomicCoords3D(new_molecule, 
                                                        new_coords) );

    if (changed_chg)
        new_params.setAtomicParameters( getChargeParameters(new_molecule,
                                                            new_chg) );

    return new_params;
}

/** Return the InterCoulombPotential::Molecule representation of 'molecule',
    using the supplied PropertyMap to find the properties that contain
    the necessary forcefield parameters
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
InterCoulombPotential::Molecule
InterCoulombPotential::parameterise(const PartialMolecule &molecule,
                                    const PropertyMap &map)
{
    return InterCoulombPotential::Molecule(molecule, *this, map);
}

/** Convert the passed group of molecules into InterCoulombPotential::Molecules,
    using the supplied PropertyMap to find the properties that contain
    the necessary forcefield parameters in each molecule
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
InterCoulombPotential::Molecules 
InterCoulombPotential::parameterise(const MoleculeGroup &molecules,
                                    const PropertyMap &map)
{
    return InterCoulombPotential::Molecules(molecules, *this, map);
}

/** Return the total charge of the parameters for the group in 'params' */
double InterCoulombPotential::totalCharge(
                        const InterCoulombPotential::Parameters::Array &params) const
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
void InterCoulombPotential::_pvt_calculateEnergy(
                                        const InterCoulombPotential::Molecule &mol0,
                                        const InterCoulombPotential::Molecule &mol1,
                                        InterCoulombPotential::Energy &energy,
                                        InterCoulombPotential::EnergyWorkspace &distmat,
                                        double scale_energy) const
{
    const quint32 ngroups0 = mol0.coordinates().count();
    const quint32 ngroups1 = mol1.coordinates().count();
    
    if (ngroups0 < ngroups1)
    {
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
    
    #ifdef SIRE_TIME_ROUTINES
    int nflops = 0;
    #endif

    //loop over all pairs of CutGroups in the two molecules
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
                                        spce->beyond(switchfunc->cutoffDistance(), 
                                                     aabox0, group1.aaBox());
            
            if (not within_cutoff)
                //this CutGroup is either the cutoff distance
                continue;
            
            //calculate all of the interatomic distances
            const double mindist = spce->calcDist(group0, group1, distmat);
            
            if (mindist > switchfunc->cutoffDistance())
            {
                //all of the atoms are definitely beyond cutoff
                continue;
            }
               
            double icnrg = 0;
            
            //loop over all interatomic pairs and calculate the energies
            const quint32 nats1 = group1.count();
            const Parameter *params1_array = params1.constData();

            #ifdef SIRE_USE_SSE
            {
                const int remainder = nats1 % 2;
                
                __m128d sse_cnrg = { 0, 0 };

                const __m128d sse_one = { 1.0, 1.0 };
                
                for (quint32 i=0; i<nats0; ++i)
                {
                    distmat.setOuterIndex(i);
                    const Parameter &param0 = params0_array[i];
                    
                    __m128d sse_chg0 = _mm_set_pd(param0.reduced_charge, 
                                                  param0.reduced_charge);
                                         
                    //process atoms in pairs (so can then use SSE)
                    for (quint32 j=0; j<nats1-1; j += 2)
                    {
                        const Parameter &param10 = params1_array[j];
                        const Parameter &param11 = params1_array[j+1];
                        
                        __m128d sse_dist = _mm_set_pd( distmat[j], distmat[j+1] );
                        __m128d sse_chg1 = _mm_set_pd( param10.reduced_charge,
                                                       param11.reduced_charge );
                                           
                        sse_dist = _mm_div_pd(sse_one, sse_dist);
                        
                        //calculate the coulomb energy
                        __m128d tmp = _mm_mul_pd(sse_chg0, sse_chg1);
                        tmp = _mm_mul_pd(tmp, sse_dist);
                        sse_cnrg = _mm_add_pd(sse_cnrg, tmp);
                        
                        #ifdef SIRE_TIME_ROUTINES
                        nflops += 8;
                        #endif
                    }
                          
                    if (remainder == 1)
                    {
                        const Parameter &param1 = params1_array[nats1-1];

                        const double invdist = double(1) / distmat[nats1-1];
                        
                        icnrg += param0.reduced_charge * param1.reduced_charge 
                                    * invdist;
                    
                        #ifdef SIRE_TIME_ROUTINES
                        nflops += 4;
                        #endif
                    }
                }
                
                icnrg += *((const double*)&sse_cnrg) +
                         *( ((const double*)&sse_cnrg) + 1 );
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
                    
                        #ifdef SIRE_TIME_ROUTINES
                        nflops += 4;
                        #endif
                    }
                }
            }
            #endif
            
            //are we shifting the electrostatic potential?
            if (use_electrostatic_shifting)
            {
                icnrg -= this->totalCharge(params0) * this->totalCharge(params1)
                              / switchfunc->electrostaticCutoffDistance();
                        
                #ifdef SIRE_TIME_ROUTINES      
                nflops += 3;
                #endif
            }
            
            //now add these energies onto the total for the molecule,
            //scaled by any non-bonded feather factor
            if (mindist > switchfunc->featherDistance())
            {
                cnrg += switchfunc->electrostaticScaleFactor( Length(mindist) ) * icnrg;
                
                #ifdef SIRE_TIME_ROUTINES
                nflops += 2;
                #endif
            }
            else
            {
                cnrg += icnrg;
                
                #ifdef SIRE_TIME_ROUTINES
                nflops += 1;
                #endif
            }
        }
    }
    
    //add this molecule pair's energy onto the total
    energy += Energy(scale_energy * cnrg);
    
    #ifdef SIRE_TIME_ROUTINES
    nflops += 2;
    ADD_FLOPS(nflops);
    #endif
}

/** Add to the forces in 'forces0' the forces acting on 'mol0' caused
    by 'mol1' */
void InterCoulombPotential::_pvt_calculateCoulombForce(
                                        const InterCoulombPotential::Molecule &mol0, 
                                        const InterCoulombPotential::Molecule &mol1,
                                        MolForceTable &forces0, 
                                        InterCoulombPotential::ForceWorkspace &distmat,
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

/** Add to the potentials in 'pots0' the coulomb potential acting on 'mol0' caused
    by 'mol1' */
void InterCoulombPotential::_pvt_calculateCoulombPotential(
                                   const InterCoulombPotential::Molecule &mol0, 
                                   const InterCoulombPotential::Molecule &mol1,
                                   const CoulombProbe &probe,
                                   MolPotentialTable &pots0, 
                                   InterCoulombPotential::PotentialWorkspace &distmat,
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

/** Add to the potentials in 'pots0' the coulomb potential on the passed  
    grid caused by 'mol' */
void InterCoulombPotential::_pvt_calculateCoulombPotential(
                                    const InterCoulombPotential::Molecule &mol, 
                                    const CoulombProbe &probe,
                                    GridPotentialTable &pots, 
                                    InterCoulombPotential::PotentialWorkspace &distmat,
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

/** Calculate the coulomb field caused by the molecule 'mol' on the grid points in 
    'fields' */
void InterCoulombPotential::_pvt_calculateCoulombField(
                                        const InterCoulombPotential::Molecule &mol,
                                        const CoulombProbe &probe,
                                        GridFieldTable &fields,
                                        InterCoulombPotential::FieldWorkspace &distmat,
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

/** Add to the fields in 'fields0' the fields acting on the passed probe
    at the atom points in 'mol0' caused by 'mol1' */
void InterCoulombPotential::_pvt_calculateCoulombField(
                                       const InterCoulombPotential::Molecule &mol0, 
                                       const InterCoulombPotential::Molecule &mol1,
                                       const CoulombProbe &probe,
                                       MolFieldTable &fields0, 
                                       InterCoulombPotential::FieldWorkspace &distmat,
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
                            total_field -= (q2 / distmat[j].length2()) *
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

/////////////
///////////// Implementation of IntraCoulombPotential
/////////////

static const RegisterMetaType<IntraCoulombPotential> r_intracoul( MAGIC_ONLY, NO_ROOT,
                                            IntraCoulombPotential::typeName() );

/** Serialise to a binary datastream */
QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds,
                                      const IntraCoulombPotential &intracoul)
{
    writeHeader(ds, r_intracoul, 1);
    
    ds << static_cast<const CoulombPotential&>(intracoul);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds,
                                      IntraCoulombPotential &intracoul)
{
    VersionID v = readHeader(ds, r_intracoul);
    
    if (v == 1)
    {
        ds >> static_cast<CoulombPotential&>(intracoul);
    }
    else
        throw version_error(v, "1", r_intracoul, CODELOC);
        
    return ds;
}

/** Constructor */
IntraCoulombPotential::IntraCoulombPotential() : CoulombPotential()
{}

/** Copy constructor */
IntraCoulombPotential::IntraCoulombPotential(const IntraCoulombPotential &other)
                      : CoulombPotential(other)
{}

/** Destructor */
IntraCoulombPotential::~IntraCoulombPotential()
{}

/** Copy assignment operator */
IntraCoulombPotential& IntraCoulombPotential::operator=(const IntraCoulombPotential &other)
{
    CoulombPotential::operator=(other);
    return *this;
}

void IntraCoulombPotential::throwMissingForceComponent(const Symbol &symbol,
                              const IntraCoulombPotential::Components &components) const
{
    throw SireFF::missing_component( QObject::tr(
        "There is no force component in potential %1 - available "
        "components are %2.")
            .arg(this->what())
            .arg(components.total().toString()), CODELOC );
}

/** Assert that 'rest_of_mol' is compatible with 'mol'. They are only 
    compatible if they are both part of the same molecule (not necessarily
    the same version) with the same layout UID.
    
    \throw SireError::incompatible_error
*/
void IntraCoulombPotential::assertCompatible(const IntraCoulombPotential::Molecule &mol,
                                const IntraCoulombPotential::Molecule &rest_of_mol) const
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
IntraCoulombPotential::Parameters 
IntraCoulombPotential::getParameters(const PartialMolecule &molecule,
                                 const PropertyMap &map)
{
    return Parameters( AtomicParameters3D<ChargeParameter>(
                               molecule, map[parameters().coordinates()],
                               getChargeParameters(molecule,map[parameters().charge()]) ),
                       IntraScaledParameters<CoulombNBPairs>(
                               molecule, map[parameters().intraScaleFactors()] )
                     );
}

/** Update the parameters for the molecule going from 'old_molecule' to 
    'new_molecule', with the parameters found using the property map 'map'
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
IntraCoulombPotential::Parameters
IntraCoulombPotential::updateParameters(
                                    const IntraCoulombPotential::Parameters &old_params,
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
    const PropertyName &scl_property = map[parameters().intraScaleFactors()];
    
    //get what has changed
    bool new_coords = old_molecule.version(coords_property) !=
                         new_molecule.version(coords_property);
                             
    bool new_chg = ( old_molecule.version(chg_property) !=
                         new_molecule.version(chg_property) );

    bool new_scl = ( old_molecule.version(scl_property) !=
                         new_molecule.version(scl_property) );

    if (new_coords)
        new_params.setAtomicCoordinates( AtomicCoords3D(new_molecule, 
                                                        coords_property) );

    if (new_chg)
        new_params.setAtomicParameters( getChargeParameters(new_molecule,chg_property) );

    if (new_scl)
        new_params.setIntraScaleFactors( 
                IntraScaledParameters<CoulombNBPairs>(new_molecule, scl_property) );

    return new_params;
}
                 
/** Update the parameters for the molecule going from 'old_molecule' to 
    'new_molecule', also while the parameters of 'old_molecule'
    where found in 'old_map', now get the parameters using 'new_map'
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
IntraCoulombPotential::Parameters
IntraCoulombPotential::updateParameters(
                                    const IntraCoulombPotential::Parameters &old_params,
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
    const PropertyName &old_scl = old_map[parameters().intraScaleFactors()];
    
    const PropertyName &new_coords = new_map[parameters().coordinates()];
    const PropertyName &new_chg = new_map[parameters().charge()];
    const PropertyName &new_scl = new_map[parameters().intraScaleFactors()];
    
    //get what has changed
    bool changed_coords = (new_coords != old_coords) or
                           old_molecule.version(old_coords) !=
                           new_molecule.version(old_coords);
                             
    bool changed_chg = (new_chg != old_chg) or
                       ( old_molecule.version(old_chg) !=
                         new_molecule.version(old_chg) );

    bool changed_scl = (new_scl != old_scl) or
                        old_molecule.version(old_scl) !=
                        new_molecule.version(old_scl);

    if (changed_coords)
        new_params.setAtomicCoordinates( AtomicCoords3D(new_molecule, new_coords) );

    if (changed_chg)
        new_params.setAtomicParameters( getChargeParameters(new_molecule,
                                                            new_chg) );

    if (changed_scl)
        new_params.setIntraScaleFactors( 
                        IntraScaledParameters<CoulombNBPairs>(new_molecule, new_scl) );

    return new_params;
}

/** Return the IntraCoulombPotential::Molecule representation of 'molecule',
    using the supplied PropertyMap to find the properties that contain
    the necessary forcefield parameters
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
IntraCoulombPotential::Molecule
IntraCoulombPotential::parameterise(const PartialMolecule &molecule,
                                    const PropertyMap &map)
{
    return IntraCoulombPotential::Molecule(molecule, *this, map);
}

/** Concert the passed group of molecules into IntraCoulombPotential::Molecules,
    using the supplied PropertyMap to find the properties that contain
    the necessary forcefield parameters in each molecule
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
IntraCoulombPotential::Molecules 
IntraCoulombPotential::parameterise(const MoleculeGroup &molecules,
                                    const PropertyMap &map)
{
    return IntraCoulombPotential::Molecules(molecules, *this, map);
}

/** Return the total charge of the parameters for the group in 'params' */
double IntraCoulombPotential::totalCharge(
                        const IntraCoulombPotential::Parameters::Array &params) const
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

void IntraCoulombPotential::calculateEnergy(const CoulombNBPairs::CGPairs &group_pairs, 
                            IntraCoulombPotential::EnergyWorkspace &distmat,
                            const IntraCoulombPotential::Parameter *params0_array, 
                            const IntraCoulombPotential::Parameter *params1_array,
                            const quint32 nats0, const quint32 nats1, 
                            double &icnrg) const
{
    if (group_pairs.isEmpty())
    {
        //there are no scale factors between groups...
        for (quint32 i=0; i<nats0; ++i)
        {
            distmat.setOuterIndex(i);
            const Parameter &param0 = params0_array[i];
                
            for (quint32 j=0; j<nats1; ++j)
            {
                icnrg += param0.reduced_charge * 
                         params1_array[j].reduced_charge / distmat[j];
            }
        }
    }
    else
    {
        //there are different nb scale factors between
        //the atoms. We need to calculate the energies using
        //them...
        for (quint32 i=0; i<nats0; ++i)
        {
            distmat.setOuterIndex(i);
            const Parameter &param0 = params0_array[i];
                
            for (quint32 j=0; j<nats1; ++j)
            {
                const CoulombScaleFactor &coulscl = group_pairs(i,j);
                            
                if (coulscl.coulomb() != 0)
                       icnrg += coulscl.coulomb() * 
                                param0.reduced_charge * 
                                params1_array[j].reduced_charge / distmat[j];
            }
        }
    }
}

void IntraCoulombPotential::calculateEnergy(const CoulombNBPairs::CGPairs &group_pairs, 
                            const QSet<Index> &atoms0, const QSet<Index> &atoms1,
                            IntraCoulombPotential::EnergyWorkspace &distmat,
                            const IntraCoulombPotential::Parameter *params0_array, 
                            const IntraCoulombPotential::Parameter *params1_array,
                            const quint32 nats0, const quint32 nats1, 
                            double &icnrg) const
{
    if (atoms0.isEmpty() or atoms1.isEmpty())
        return;

    if (group_pairs.isEmpty())
    {
        //there are no scale factors between groups...
        foreach (Index i, atoms0)
        {
            distmat.setOuterIndex(i);
            const Parameter &param0 = params0_array[i];
                
            foreach (Index j, atoms1)
            {
                icnrg += param0.reduced_charge * 
                         params1_array[j].reduced_charge / distmat[j];
            }
        } 
    }
    else
    {
        //there are different nb scale factors between
        //the atoms. We need to calculate the energies using
        //them...
        foreach (Index i, atoms0)
        {
            distmat.setOuterIndex(i);
            const Parameter &param0 = params0_array[i];
                
            foreach (Index j, atoms1)
            {
                const CoulombScaleFactor &coulscl = group_pairs(i,j);
                            
                if (coulscl.coulomb() != 0)
                    icnrg += coulscl.coulomb() * 
                             param0.reduced_charge * 
                             params1_array[j].reduced_charge / distmat[j];
            }
        }
    }
}

/** Calculate the intramolecular coulomb energy of the passed molecule, and
    add this onto 'energy'. This uses the passed workspace when
    performing the calculation */
void IntraCoulombPotential::calculateEnergy(const IntraCoulombPotential::Molecule &mol,
                                        IntraCoulombPotential::Energy &energy,
                                        IntraCoulombPotential::EnergyWorkspace &distmat,
                                        double scale_energy) const
{
    if (scale_energy == 0 or mol.isEmpty())
        return;

    const quint32 ngroups = mol.nCutGroups();
    
    const CoordGroup *groups_array = mol.coordinates().constData();
    
    BOOST_ASSERT( mol.parameters().atomicParameters().count() == int(ngroups) );
    const Parameters::Array *molparams_array 
                            = mol.parameters().atomicParameters().constData();

    const CoulombNBPairs &nbpairs = mol.parameters().intraScaleFactors();
    
    double cnrg = 0;
      
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
            const CoulombNBPairs::CGPairs &group_pairs = nbpairs(cgidx_igroup,
                                                                 cgidx_jgroup);

            double icnrg = 0;
            
            //loop over all intraatomic pairs and calculate the energies
            const quint32 nats1 = group1.count();
            const Parameter *params1_array = params1.constData();
            
            calculateEnergy(group_pairs, distmat, params0_array, params1_array,
                            nats0, nats1, icnrg);
            
            //if this is the same group then half the energies to 
            //correct for double-counting
            if (igroup == jgroup)
            {
                icnrg *= 0.5;
            }

            //are we shifting the electrostatic potential?
            if (use_electrostatic_shifting and igroup != jgroup)
                icnrg -= this->totalCharge(params0) * this->totalCharge(params1)
                              / switchfunc->electrostaticCutoffDistance();

            //now add these energies onto the total for the molecule,
            //scaled by any non-bonded feather factor
            if (mindist > switchfunc->featherDistance())
            {
                cnrg += switchfunc->electrostaticScaleFactor( Length(mindist) ) * icnrg;
            }
            else
            {
                cnrg += icnrg;
            }
        }
    }
    
    //add this molecule pair's energy onto the total
    energy += Energy(scale_energy * cnrg);
}

/** Calculate the intramolecular coulomb energy of the passed molecule
    interacting with the rest of the same molecule in 'rest_of_mol', and
    add this onto 'energy'. This uses the passed workspace when
    performing the calculation. Note that mol and rest_of_mol should
    not contain any overlapping atoms, and that they should both be
    part of the same molecule (albeit potentially at different versions,
    but with the same layout UID)
    
    \throw SireError::incompatible_error
*/
void IntraCoulombPotential::calculateEnergy(const IntraCoulombPotential::Molecule &mol,
                                        const IntraCoulombPotential::Molecule &rest_of_mol,
                                        IntraCoulombPotential::Energy &energy,
                                        IntraCoulombPotential::EnergyWorkspace &distmat,
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

    //the CoulombNBPairs must be the same in both molecules - this is checked
    //as part of assertCompatible(..)
    const CoulombNBPairs &nbpairs = mol.parameters().intraScaleFactors();
    
    double cnrg = 0;
    
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
            const CoulombNBPairs::CGPairs &group_pairs = nbpairs(cgidx_igroup,
                                                                 cgidx_jgroup);

            double icnrg = 0;
            
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
                                nats0, nats1, icnrg);
            }
            else
            {
                calculateEnergy(group_pairs, distmat,
                                params0_array, params1_array,
                                nats0, nats1, icnrg);
            }

            //are we shifting the electrostatic potential?
            if (use_electrostatic_shifting and cgidx_igroup != cgidx_jgroup)
                icnrg -= this->totalCharge(params0) * this->totalCharge(params1)
                              / switchfunc->electrostaticCutoffDistance();

            //now add these energies onto the total for the molecule,
            //scaled by any non-bonded feather factor
            if (mindist > switchfunc->featherDistance())
            {
                cnrg += switchfunc->electrostaticScaleFactor( Length(mindist) ) * icnrg;
            }
            else
            {
                cnrg += icnrg;
            }
        }
    }

    //add the molecule's energy onto the total
    energy += Energy(scale_energy * cnrg);
}

void IntraCoulombPotential::calculateCoulombForce(
                             const CoulombNBPairs::CGPairs &group_pairs,
                             const CoordGroup &group0, const CoordGroup &group1,
                             const double mindist,
                             IntraCoulombPotential::ForceWorkspace &distmat,
                             const IntraCoulombPotential::Parameter *params0_array,
                             const IntraCoulombPotential::Parameter *params1_array,
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
                    const CoulombScaleFactor &coulscl = group_pairs(i,j);
                            
                    if (coulscl.coulomb() != 0)
                    {
                        const double q2 = param0.reduced_charge * 
                                          params1_array[j].reduced_charge;
                                                      
                        if (q2 != 0)
                        {
                            //calculate the coulomb energy
                            const double cnrg = coulscl.coulomb() *
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
                    const CoulombScaleFactor &coulscl = group_pairs(i,j);
                          
                    double cscale = coulscl.coulomb() 
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

void IntraCoulombPotential::calculateCoulombForce(
                             const CoulombNBPairs::CGPairs &group_pairs,
                             const QSet<Index> &atoms0, const QSet<Index> &atoms1,
                             const CoordGroup &group0, const CoordGroup &group1,
                             const double mindist,
                             IntraCoulombPotential::ForceWorkspace &distmat,
                             const IntraCoulombPotential::Parameter *params0_array,
                             const IntraCoulombPotential::Parameter *params1_array,
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
                    const CoulombScaleFactor &coulscl = group_pairs(i,j);
                            
                    if (coulscl.coulomb() != 0)
                    {
                        const double q2 = param0.reduced_charge * 
                                          params1_array[j].reduced_charge;
                                                      
                        if (q2 != 0)
                        {
                            //calculate the coulomb energy
                            const double cnrg = coulscl.coulomb() *
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
                    const CoulombScaleFactor &coulscl = group_pairs(i,j);
                          
                    double cscale = coulscl.coulomb() 
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
void IntraCoulombPotential::calculateCoulombForce(
                                       const IntraCoulombPotential::Molecule &mol,
                                       MolForceTable &forces, 
                                       IntraCoulombPotential::ForceWorkspace &distmat,
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

    const CoulombNBPairs &nbpairs = mol.parameters().intraScaleFactors();
    
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
            const CoulombNBPairs::CGPairs &group_pairs = nbpairs(cgidx_igroup,
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

/** Calculate the coulomb forces between the atoms in the molecule 'mol'
    and add these forces onto 'forces'. This uses
    the passed workspace to perform the calculation */
void IntraCoulombPotential::calculateForce(
                                       const IntraCoulombPotential::Molecule &mol,
                                       MolForceTable &forces, 
                                       IntraCoulombPotential::ForceWorkspace &distmat,
                                       double scale_force) const
{
    this->calculateCoulombForce(mol, forces, distmat, scale_force);
}

/** Calculate the coulomb force acting on the part of the molecule
    in 'mol' caused by the rest of the molecule in 'rest_of_mol'. Note
    that these must both be of the same molecule, with the same
    layout UID and same nonbonded scale factors
    
    \throw SireError::incompatible_error
*/
void IntraCoulombPotential::calculateCoulombForce(
                                        const IntraCoulombPotential::Molecule &mol,
                                        const IntraCoulombPotential::Molecule &rest_of_mol,
                                        MolForceTable &forces,
                                        IntraCoulombPotential::ForceWorkspace &distmat,
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
    const CoulombNBPairs &nbpairs = mol.parameters().intraScaleFactors();

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
            const CoulombNBPairs::CGPairs &group_pairs = nbpairs(cgidx_igroup,
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

/** Calculate the coulomb force acting on the part of the molecule
    in 'mol' caused by the rest of the molecule in 'rest_of_mol'. Note
    that these must both be of the same molecule, with the same
    layout UID and same nonbonded scale factors
    
    \throw SireError::incompatible_error
*/
void IntraCoulombPotential::calculateForce(
                                        const IntraCoulombPotential::Molecule &mol,
                                        const IntraCoulombPotential::Molecule &rest_of_mol,
                                        MolForceTable &forces,
                                        IntraCoulombPotential::ForceWorkspace &distmat,
                                        double scale_force) const
{
    this->calculateCoulombForce(mol, rest_of_mol, forces, distmat, scale_force);
}

void IntraCoulombPotential::calculateField(const IntraCoulombPotential::Molecule &mol, 
                                      const CoulombProbe &probe,
                                      MolFieldTable &fields,
                                      IntraCoulombPotential::FieldWorkspace &workspace,
                                      double scale_field) const
{
    throw SireError::incomplete_code( QObject::tr(
                "Calculating the intramolecular coulomb field has not "
                "yet been implemented."), CODELOC );
}

void IntraCoulombPotential::calculateField(const IntraCoulombPotential::Molecule &mol,
                    const IntraCoulombPotential::Molecule &rest_of_mol,
                    const CoulombProbe &probe,
                    MolFieldTable &fields,
                    IntraCoulombPotential::FieldWorkspace &workspace,
                    double scale_field) const
{
    throw SireError::incomplete_code( QObject::tr(
                "Calculating the intramolecular coulomb field has not "
                "yet been implemented."), CODELOC );
}

void IntraCoulombPotential::calculateField(const IntraCoulombPotential::Molecule &mol, 
                    const CoulombProbe &probe,
                    MolFieldTable &fields,
                    const Symbol &symbol,
                    const Components &components,
                    IntraCoulombPotential::FieldWorkspace &workspace,
                    double scale_field) const
{
    throw SireError::incomplete_code( QObject::tr(
                "Calculating the intramolecular coulomb field has not "
                "yet been implemented."), CODELOC );
}

void IntraCoulombPotential::calculateField(const IntraCoulombPotential::Molecule &mol,
                    const IntraCoulombPotential::Molecule &rest_of_mol,
                    const CoulombProbe &probe,
                    MolFieldTable &fields,
                    const Symbol &symbol,
                    const Components &components,
                    IntraCoulombPotential::FieldWorkspace &workspace,
                    double scale_field) const
{
    throw SireError::incomplete_code( QObject::tr(
                "Calculating the intramolecular coulomb field has not "
                "yet been implemented."), CODELOC );
}

void IntraCoulombPotential::calculateField(const IntraCoulombPotential::Molecule &mol, 
                    const CoulombProbe &probe,
                    GridFieldTable &fields,
                    IntraCoulombPotential::FieldWorkspace &workspace,
                    double scale_field) const
{
    throw SireError::incomplete_code( QObject::tr(
                "Calculating the intramolecular coulomb field has not "
                "yet been implemented."), CODELOC );
}

void IntraCoulombPotential::calculateField(const IntraCoulombPotential::Molecule &mol, 
                    const CoulombProbe &probe,
                    GridFieldTable &fields,
                    const Symbol &symbol,
                    const Components &components,
                    IntraCoulombPotential::FieldWorkspace &workspace,
                    double scale_field) const
{
    throw SireError::incomplete_code( QObject::tr(
                "Calculating the intramolecular coulomb field has not "
                "yet been implemented."), CODELOC );
}

void IntraCoulombPotential::calculatePotential(const IntraCoulombPotential::Molecule &mol, 
                        const CoulombProbe &probe,
                        MolPotentialTable &potentials,
                        IntraCoulombPotential::PotentialWorkspace &workspace,
                        double scale_potential) const
{
    throw SireError::incomplete_code( QObject::tr(
                "Calculating the intramolecular coulomb potential has not "
                "yet been implemented."), CODELOC );
}

void IntraCoulombPotential::calculatePotential(const IntraCoulombPotential::Molecule &mol,
                        const IntraCoulombPotential::Molecule &rest_of_mol,
                        const CoulombProbe &probe,
                        MolPotentialTable &potentials,
                        IntraCoulombPotential::PotentialWorkspace &workspace,
                        double scale_potential) const
{
    throw SireError::incomplete_code( QObject::tr(
                "Calculating the intramolecular coulomb potential has not "
                "yet been implemented."), CODELOC );
}

void IntraCoulombPotential::calculatePotential(const IntraCoulombPotential::Molecule &mol, 
                        const CoulombProbe &probe,
                        MolPotentialTable &potentials,
                        const Symbol &symbol,
                        const Components &components,
                        IntraCoulombPotential::PotentialWorkspace &workspace,
                        double scale_potential) const
{
    throw SireError::incomplete_code( QObject::tr(
                "Calculating the intramolecular coulomb potential has not "
                "yet been implemented."), CODELOC );
}

void IntraCoulombPotential::calculatePotential(const IntraCoulombPotential::Molecule &mol,
                        const IntraCoulombPotential::Molecule &rest_of_mol,
                        const CoulombProbe &probe,
                        MolPotentialTable &potentials,
                        const Symbol &symbol,
                        const Components &components,
                        IntraCoulombPotential::PotentialWorkspace &workspace,
                        double scale_potential) const
{
    throw SireError::incomplete_code( QObject::tr(
                "Calculating the intramolecular coulomb potential has not "
                "yet been implemented."), CODELOC );
}

void IntraCoulombPotential::calculatePotential(const IntraCoulombPotential::Molecule &mol, 
                        const CoulombProbe &probe,
                        GridPotentialTable &potentials,
                        IntraCoulombPotential::PotentialWorkspace &workspace,
                        double scale_potential) const
{
    throw SireError::incomplete_code( QObject::tr(
                "Calculating the intramolecular coulomb potential has not "
                "yet been implemented."), CODELOC );
}

void IntraCoulombPotential::calculatePotential(const IntraCoulombPotential::Molecule &mol, 
                        const CoulombProbe &probe,
                        GridPotentialTable &potentials,
                        const Symbol &symbol,
                        const Components &components,
                        IntraCoulombPotential::PotentialWorkspace &workspace,
                        double scale_potential) const
{
    throw SireError::incomplete_code( QObject::tr(
                "Calculating the intramolecular coulomb potential has not "
                "yet been implemented."), CODELOC );
}

void IntraCoulombPotential::calculateCoulombField(const IntraCoulombPotential::Molecule &mol,
                           const CoulombProbe &probe,
                           MolFieldTable &fields,
                           IntraCoulombPotential::FieldWorkspace &workspace,
                           double scale_field) const
{
    throw SireError::incomplete_code( QObject::tr(
                "Calculating the intramolecular coulomb field has not "
                "yet been implemented."), CODELOC );
}

void IntraCoulombPotential::calculateCoulombField(const IntraCoulombPotential::Molecule &mol,
                           const CoulombProbe &probe,
                           GridFieldTable &fields,
                           IntraCoulombPotential::FieldWorkspace &workspace,
                           double scale_field) const
{
    throw SireError::incomplete_code( QObject::tr(
                "Calculating the intramolecular coulomb field has not "
                "yet been implemented."), CODELOC );
}

void IntraCoulombPotential::calculateCoulombField(const IntraCoulombPotential::Molecule &mol,
                           const IntraCoulombPotential::Molecule &rest_of_mol,
                           const CoulombProbe &probe,
                           MolFieldTable &fields,
                           IntraCoulombPotential::FieldWorkspace &workspace,
                           double scale_field) const
{
    throw SireError::incomplete_code( QObject::tr(
                "Calculating the intramolecular coulomb field has not "
                "yet been implemented."), CODELOC );
}

void IntraCoulombPotential::calculateCoulombPotential(const IntraCoulombPotential::Molecule &mol,
                               const CoulombProbe &probe,
                               MolPotentialTable &potentials,
                               IntraCoulombPotential::PotentialWorkspace &workspace,
                               double scale_potential) const
{
    throw SireError::incomplete_code( QObject::tr(
                "Calculating the intramolecular coulomb potential has not "
                "yet been implemented."), CODELOC );
}

void IntraCoulombPotential::calculateCoulombPotential(const IntraCoulombPotential::Molecule &mol,
                               const CoulombProbe &probe,
                               GridPotentialTable &fields,
                               IntraCoulombPotential::PotentialWorkspace &workspace,
                               double scale_potential) const
{
    throw SireError::incomplete_code( QObject::tr(
                "Calculating the intramolecular coulomb potential has not "
                "yet been implemented."), CODELOC );
}

void IntraCoulombPotential::calculateCoulombPotential(const IntraCoulombPotential::Molecule &mol,
                               const IntraCoulombPotential::Molecule &rest_of_mol,
                               const CoulombProbe &probe,
                               MolPotentialTable &potentials,
                               IntraCoulombPotential::PotentialWorkspace &workspace,
                               double scale_potential) const
{
    throw SireError::incomplete_code( QObject::tr(
                "Calculating the intramolecular coulomb potential has not "
                "yet been implemented."), CODELOC );
}
