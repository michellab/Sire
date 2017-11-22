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

#include "ljpotential.h"
#include "ljparameter.h"
#include "switchingfunction.h"

#include "ljfunction.h"

#include "SireMol/mover.hpp"
#include "SireMol/atomcoords.h"

#include "SireVol/cartesian.h"

#include "SireBase/countflops.h"
#include "SireBase/propertylist.h"

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

#undef SIRE_USE_SSE

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
/////// LJParameterName
///////

QString LJParameterName::lj_param( "LJ" );

///////
/////// Completely instantiate the LJPotential ancillary classes
///////

namespace SireFF
{
    namespace detail
    {
        template class AtomicParameters3D<LJParamID>;

        template class FFMolecule3D<InterLJPotential>;

        template class FFMolecules3D<InterLJPotential>;

        template class ChangedMolecule<InterLJPotential::Molecule>;

        template class FFMolecule3D<IntraLJPotential>;

        template class FFMolecules3D<IntraLJPotential>;

        template class ChangedMolecule<IntraLJPotential::Molecule>;
    }
}

namespace SireMM
{
    namespace detail
    {
        template class IntraScaledParameters<LJNBPairs>;

        template class IntraScaledAtomicParameters< AtomicParameters3D<LJParamID>,
                                                    IntraScaledParameters<LJNBPairs> >;
    }
}

/** Streaming functions for LJParamID - these must convert the 
    LJID number to and from an actual LJParameter (as the LJParameterDB
    will have different ID numbers of different processors) */
QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, 
                                      const SireMM::detail::LJParamID &ljparam)
{
    ds << LJParameterDB::getLJParameter(ljparam.ljid);
    return ds;
}

QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, 
                                      SireMM::detail::LJParamID &ljparam)
{
    LJParameter lj;

    ds >> lj;
    
    ljparam.ljid = LJParameterDB::addLJParameter(lj);
    
    return ds;
}

/** Internal function used to get the LJ parameters from a molecule
    and convert them into a PackedArray of LJ parameter IDs */
static PackedArray2D<LJParamID> getLJParamIDs(const PartialMolecule &molecule,
                                              const PropertyName &lj_property)
{
    const AtomLJs &ljs = molecule.property(lj_property).asA<AtomLJs>();
    
    const AtomSelection &selected_atoms = molecule.selection();
    
    if (selected_atoms.selectedNone())
        return PackedArray2D<LJParamID>();
    
    //create space for the parameters - only need space for CutGroups
    //that contain at least one selected atom
    QVector< QVector<LJParamID> > ljparams( selected_atoms.nSelectedCutGroups() );
    QVector<LJParamID>* ljparams_array = ljparams.data();

    try
    {

    LJParameterDB::lock();

    if (selected_atoms.selectedAllCutGroups())
    {
        const int ncg = molecule.data().info().nCutGroups();
    
        for (CGIdx i(0); i<ncg; ++i)
        {
            const int nats = molecule.data().info().nAtoms(i);
            
            QVector<LJParamID> group_ljones(nats);
            LJParamID *group_ljones_array = group_ljones.data();
            
            //get the arrays containing the LJ parameters
            //for this CutGroup
            const LJParameter *group_ljs = ljs.constData(i);
            
            if (selected_atoms.selectedAll(i))
            {
                for (Index j(0); j<nats; ++j)
                {
                    group_ljones_array[j].ljid = 
                            LJParameterDB::_locked_addLJParameter(group_ljs[j]);
                }
            }
            else
            {
                foreach (Index j, selected_atoms.selectedAtoms(i))
                {
                    group_ljones_array[j].ljid =
                            LJParameterDB::_locked_addLJParameter(group_ljs[j]);
                }
            }
            
            ljparams_array[i] = group_ljones;
        }
    }
    else
    {
        foreach (CGIdx i, selected_atoms.selectedCutGroups())
        {
            const int nats = molecule.data().info().nAtoms(i);
            
            QVector<LJParamID> group_cljones(nats);
            LJParamID *group_cljones_array = group_cljones.data();
            
            //get the arrays containing the LJ parameters
            //for this CutGroup
            const LJParameter *group_ljs = ljs.constData(i);
            
            if (selected_atoms.selectedAll(i))
            {
                for (Index j(0); j<nats; ++j)
                {
                    group_cljones_array[j].ljid = 
                            LJParameterDB::_locked_addLJParameter(group_ljs[j]);
                }
            }
            else
            {
                foreach (Index j, selected_atoms.selectedAtoms(i))
                {
                    group_cljones_array[j].ljid =
                            LJParameterDB::_locked_addLJParameter(group_ljs[j]);
                }
            }
            
            ljparams_array[i] = group_cljones;
        }
    }
    
    LJParameterDB::unlock();
    
    }
    catch(...)
    {
        LJParameterDB::unlock();
        throw;
    }
    
    return PackedArray2D<LJParamID>( ljparams );
}

/////////////
///////////// Implementation of LJPotential
/////////////

static const RegisterMetaType<LJPotential> r_ljpot( MAGIC_ONLY, NO_ROOT,
                                                    "SireMM::LJPotential" );

/** Serialise to a binary datastream */
QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds,
                                     const LJPotential &ljpot)
{
    writeHeader(ds, r_ljpot, 1);
    
    ds << ljpot.props;

    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds,
                                      LJPotential &ljpot)
{
    VersionID v = readHeader(ds, r_ljpot);
    
    if (v == 1)
    {
        ds >> ljpot.props;
    
        //extract all of the properties
        ljpot.spce = ljpot.props.property("space").asA<Space>();
        ljpot.switchfunc = ljpot.props.property("switchingFunction")
                                      .asA<SwitchingFunction>();
    
        ljpot.combining_rules = LJParameterDB::interpret(
                                    ljpot.props.property("combiningRules")
                                               .asAString() );
                                        
        ljpot.need_update_ljpairs = true;
    }
    else 
        throw version_error(v, "1", r_ljpot, CODELOC);
    
    return ds;
}

/** Constructor */
LJPotential::LJPotential()
             : spce( Space::null() ), switchfunc( SwitchingFunction::null() ),
               combining_rules( LJParameterDB::interpret("arithmetic") ),
               need_update_ljpairs(true)
{
    //record the defaults
    props.setProperty( "space", spce );
    props.setProperty( "switchingFunction", switchfunc );
    props.setProperty( "combiningRules", wrap( LJParameterDB::toString(combining_rules) ) );
}

/** Copy constructor */
LJPotential::LJPotential(const LJPotential &other)
             : ljpairs(other.ljpairs), props(other.props),
               spce(other.spce), switchfunc(other.switchfunc),
               combining_rules(other.combining_rules),
               need_update_ljpairs(other.need_update_ljpairs)
{}

/** Destructor */
LJPotential::~LJPotential()
{}

/** Copy assignment operator */
LJPotential& LJPotential::operator=(const LJPotential &other)
{
    if (this != &other)
    {
        ljpairs = other.ljpairs;
        props = other.props;
        spce = other.spce;
        switchfunc = other.switchfunc;
        combining_rules = other.combining_rules;
        need_update_ljpairs = other.need_update_ljpairs;
    }
    
    return *this;
}

/** You need to call this function before you start a block of 
    energy of force evaluation using this forcefield. You should
    also call 'finishedEvaluation()' once you have finished. */
void LJPotential::startEvaluation()
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
void LJPotential::finishedEvaluation()
{}

/** Return all of the properties set in this forcefield */
const Properties& LJPotential::properties() const
{
    return props;
}

/** Return whether or not this potential has a property called 'name' */
bool LJPotential::containsProperty(const QString &name) const
{
    return props.hasProperty(name);
}

/** Return the property with name 'name'

    \throw SireBase::missing_property
*/
const Property& LJPotential::property(const QString &name) const
{
    return props.property(name);
}

/** Set the 3D space in which the molecules in this potential are evaluated */
bool LJPotential::setSpace(const Space &new_space)
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
bool LJPotential::setSwitchingFunction(const SwitchingFunction &new_switchfunc)
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

/** Set the combining rules to use to obtain mixed LJ parameters */
bool LJPotential::setCombiningRules(const QString &combiningrules)
{
    LJParameterDB::CombiningRules new_rules = LJParameterDB::interpret(combiningrules);
    
    if (new_rules != combining_rules)
    {
        combining_rules = new_rules;
        need_update_ljpairs = true;
        props.setProperty( "combiningRules", wrap(combiningrules) );
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
bool LJPotential::setProperty(const QString &name, const Property &value)
{
    if (name == QLatin1String("space"))
    {
        return this->setSpace( value.asA<Space>() );
    }
    else if (name == QLatin1String("switchingFunction"))
    {
        return this->setSwitchingFunction( value.asA<SwitchingFunction>() );
    }
    else if (name == QLatin1String("combiningRules"))
    {
        return this->setCombiningRules( value.asAString() );
    }
    else
        throw SireBase::missing_property( QObject::tr(
            "The CLJ potentials do not have a property called \"%1\" that "
            "can be changed. Available properties are [ %2 ].")
                .arg(name, QStringList(props.propertyKeys()).join(", ")), CODELOC );
                
    return false;
}

/** Return the 3D space in which this potential is evaluated */
const Space& LJPotential::space() const
{
    return *spce;
}

/** Return the switching function used to scale the group-group
    interactions to zero at the cutoff */
const SwitchingFunction& LJPotential::switchingFunction() const
{
    return *switchfunc;
}

/** Return the string identifying the combining rules used to 
    obtain the mixed LJ parameters */
const QString& LJPotential::combiningRules() const
{
    return LJParameterDB::toString(combining_rules);
}

/////////////
///////////// Implementation of InterLJPotential
/////////////

static const RegisterMetaType<InterLJPotential> r_interlj( MAGIC_ONLY, NO_ROOT,
                                            InterLJPotential::typeName() );

/** Serialise to a binary datastream */
QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds,
                                      const InterLJPotential &interlj)
{
    writeHeader(ds, r_interlj, 1);
    
    ds << static_cast<const LJPotential&>(interlj);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds,
                                      InterLJPotential &interlj)
{
    VersionID v = readHeader(ds, r_interlj);
    
    if (v == 1)
    {
        ds >> static_cast<LJPotential&>(interlj);
    }
    else
        throw version_error(v, "1", r_interlj, CODELOC);
        
    return ds;
}

/** Constructor */
InterLJPotential::InterLJPotential() : LJPotential()
{}

/** Copy constructor */
InterLJPotential::InterLJPotential(const InterLJPotential &other)
                 : LJPotential(other)
{}

/** Destructor */
InterLJPotential::~InterLJPotential()
{}

/** Copy assignment operator */
InterLJPotential& InterLJPotential::operator=(const InterLJPotential &other)
{
    LJPotential::operator=(other);
    return *this;
}

void InterLJPotential::throwMissingForceComponent(const Symbol &symbol,
                              const IntraLJPotential::Components &components) const
{
    throw SireFF::missing_component( QObject::tr(
        "There is no force component in potential %1 - available "
        "components are %2.")
            .arg(this->what())
            .arg(components.total().toString()), CODELOC );
}

void InterLJPotential::throwMissingFieldComponent(const Symbol &symbol,
                              const IntraLJPotential::Components &components) const
{
    throw SireFF::missing_component( QObject::tr(
        "There is no field component in potential %1 - available "
        "components are %2.")
            .arg(this->what())
            .arg(components.total().toString()), CODELOC );
}

void InterLJPotential::throwMissingPotentialComponent(const Symbol &symbol,
                              const IntraLJPotential::Components &components) const
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
InterLJPotential::Parameters 
InterLJPotential::getParameters(const PartialMolecule &molecule,
                                const PropertyMap &map)
{
    need_update_ljpairs = true;
    
    return Parameters( molecule, map[parameters().coordinates()],
                       getLJParamIDs(molecule, map[parameters().lj()]) );
}

/** Update the parameters for the molecule going from 'old_molecule' to 
    'new_molecule', with the parameters found using the property map 'map'
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
InterLJPotential::Parameters
InterLJPotential::updateParameters(const InterLJPotential::Parameters &old_params,
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
    const PropertyName &lj_property = map[parameters().lj()];
    
    //get what has changed
    bool new_coords = old_molecule.version(coords_property) !=
                         new_molecule.version(coords_property);
                             
    bool new_lj = ( old_molecule.version(lj_property) !=
                    new_molecule.version(lj_property) );

    if (new_coords)
    {
        new_params.setAtomicCoordinates( AtomicCoords3D(new_molecule, 
                                                        coords_property) );
    }

    if (new_lj)
    {
        new_params.setAtomicParameters( getLJParamIDs(new_molecule, lj_property) );
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
InterLJPotential::Parameters
InterLJPotential::updateParameters(const InterLJPotential::Parameters &old_params,
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
    const PropertyName &old_lj = old_map[parameters().lj()];
    
    const PropertyName &new_coords = new_map[parameters().coordinates()];
    const PropertyName &new_lj = new_map[parameters().lj()];
    
    //get what has changed
    bool changed_coords = (new_coords != old_coords) or
                           old_molecule.version(old_coords) !=
                           new_molecule.version(old_coords);
                             
    bool changed_lj = (new_lj != old_lj) or
                      ( old_molecule.version(old_lj) !=
                        new_molecule.version(old_lj) );

    if (changed_coords)
        new_params.setAtomicCoordinates( AtomicCoords3D(new_molecule, 
                                                        new_coords) );

    if (changed_lj)
    {
        new_params.setAtomicParameters( getLJParamIDs(new_molecule, new_lj) );
        need_update_ljpairs = true;
    }

    return new_params;
}

/** Return the InterLJPotential::Molecule representation of 'molecule',
    using the supplied PropertyMap to find the properties that contain
    the necessary forcefield parameters
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
InterLJPotential::Molecule
InterLJPotential::parameterise(const PartialMolecule &molecule,
                               const PropertyMap &map)
{
    return InterLJPotential::Molecule(molecule, *this, map);
}

/** Convert the passed group of molecules into InterLJPotential::Molecules,
    using the supplied PropertyMap to find the properties that contain
    the necessary forcefield parameters in each molecule
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
InterLJPotential::Molecules 
InterLJPotential::parameterise(const MoleculeGroup &molecules,
                               const PropertyMap &map)
{
    return InterLJPotential::Molecules(molecules, *this, map);
}

/** Calculate the coulomb and LJ energy between the passed pair
    of molecules and add these energies onto 'energy'. This uses
    the passed workspace to perform the calculation */
void InterLJPotential::_pvt_calculateEnergy(const InterLJPotential::Molecule &mol0,
                                            const InterLJPotential::Molecule &mol1,
                                            InterLJPotential::Energy &energy,
                                            InterLJPotential::EnergyWorkspace &distmat,
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
    
    double ljnrg = 0;
    
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
            const double mindist = spce->calcDist2(group0, group1, distmat);
            
            if (mindist > switchfunc->cutoffDistance())
            {
                //all of the atoms are definitely beyond cutoff
                continue;
            }
               
            double iljnrg = 0;
            
            //loop over all interatomic pairs and calculate the energies
            const quint32 nats1 = group1.count();
            const Parameter *params1_array = params1.constData();

            #ifdef SIRE_USE_SSE
            {
                const int remainder = nats1 % 2;
                
                __m128d sse_ljnrg = { 0, 0 };
                
                for (quint32 i=0; i<nats0; ++i)
                {
                    distmat.setOuterIndex(i);
                    const Parameter &param0 = params0_array[i];

                    if (param0.ljid == 0)
                        //skip dummy atoms
                        continue;
                    
                    //process atoms in pairs (so can then use SSE)
                    for (quint32 j=0; j<nats1-1; j += 2)
                    {
                        const Parameter &param10 = params1_array[j];
                        const Parameter &param11 = params1_array[j+1];

                        const LJPair &ljpair0 = ljpairs.constData()[
                                                ljpairs.map(param0.ljid,
                                                            param10.ljid)];
                    
                        const LJPair &ljpair1 = ljpairs.constData()[
                                                ljpairs.map(param0.ljid,
                                                            param11.ljid)];

                        __m128d sse_nrg = calcLJEnergy(distmat[j], distmat[j+1],
                                                       ljpair0, ljpair1);
                        
                        sse_ljnrg = _mm_add_pd( sse_ljnrg, sse_nrg );
                    }
                          
                    if (remainder == 1)
                    {
                        const Parameter &param1 = params1_array[nats1-1];

                        const LJPair &ljpair = ljpairs.constData()[
                                                ljpairs.map(param0.ljid,
                                                            param1.ljid)];
                    
                        iljnrg += calcLJEnergy(distmat[nats1-1], ljpair);
                    }
                }
                         
                iljnrg += *((const double*)&sse_ljnrg) +
                          *( ((const double*)&sse_ljnrg) + 1 );
            }
            #else
            {
                for (quint32 i=0; i<nats0; ++i)
                {
                    distmat.setOuterIndex(i);
                    const Parameter &param0 = params0_array[i];
                    
                    if (param0.ljid == 0)
                        continue;
                
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        const Parameter &param1 = params1_array[j];
                        
                        if (param1.ljid == 0)
                            continue;

                        iljnrg += calcLJEnergy(distmat[j], ljpairs.constData()[
                                                                ljpairs.map(param0.ljid,
                                                                        param1.ljid)]);
                    }
                }
            }
            #endif
            
            //now add these energies onto the total for the molecule,
            //scaled by any non-bonded feather factor
            if (mindist > switchfunc->featherDistance())
            {
                ljnrg += switchfunc->vdwScaleFactor( Length(mindist) ) * iljnrg;
            }
            else
            {
                ljnrg += iljnrg;
            }
        }
    }
    
    //add this molecule pair's energy onto the total
    energy += Energy(scale_energy * ljnrg);
}

/** Add to the forces in 'forces0' the forces acting on 'mol0' caused
    by 'mol1' */
void InterLJPotential::_pvt_calculateLJForce(
                                       const InterLJPotential::Molecule &mol0, 
                                       const InterLJPotential::Molecule &mol1,
                                       MolForceTable &forces0, 
                                       InterLJPotential::ForceWorkspace &distmat,
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
void InterLJPotential::_pvt_calculateLJPotential(
                                        const InterLJPotential::Molecule &mol0, 
                                        const InterLJPotential::Molecule &mol1,
                                        const LJProbe &probe,
                                        MolPotentialTable &pots0, 
                                        InterLJPotential::PotentialWorkspace &distmat,
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

/** Add to the potentials in 'pots0' the potential on the passed  
    grid caused by 'mol' */
void InterLJPotential::_pvt_calculateLJPotential(
                                        const InterLJPotential::Molecule &mol, 
                                        const LJProbe &probe,
                                        GridPotentialTable &pots, 
                                        InterLJPotential::PotentialWorkspace &distmat,
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

/** Calculate the field caused by the molecule 'mol' on the grid points in 
    'fields' */
void InterLJPotential::_pvt_calculateLJField(const InterLJPotential::Molecule &mol,
                                             const LJProbe &probe,
                                             GridFieldTable &fields,
                                             InterLJPotential::FieldWorkspace &distmat,
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

/** Add to the fields in 'fields0' the fields acting on the probe
    at all of the atoms sites of 'mol0' caused by 'mol1' */
void InterLJPotential::_pvt_calculateLJField(
                                            const InterLJPotential::Molecule &mol0, 
                                            const InterLJPotential::Molecule &mol1,
                                            const LJProbe &probe,
                                            MolFieldTable &fields0, 
                                            InterLJPotential::FieldWorkspace &distmat,
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
///////////// Implementation of IntraLJPotential
/////////////

static const RegisterMetaType<IntraLJPotential> r_intralj( MAGIC_ONLY, NO_ROOT,
                                            IntraLJPotential::typeName() );

/** Serialise to a binary datastream */
QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds,
                                      const IntraLJPotential &intralj)
{
    writeHeader(ds, r_intralj, 1);
    
    ds << static_cast<const LJPotential&>(intralj);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds,
                                      IntraLJPotential &intralj)
{
    VersionID v = readHeader(ds, r_intralj);
    
    if (v == 1)
    {
        ds >> static_cast<LJPotential&>(intralj);
    }
    else
        throw version_error(v, "1", r_intralj, CODELOC);
        
    return ds;
}

/** Constructor */
IntraLJPotential::IntraLJPotential() : LJPotential()
{}

/** Copy constructor */
IntraLJPotential::IntraLJPotential(const IntraLJPotential &other)
                  : LJPotential(other)
{}

/** Destructor */
IntraLJPotential::~IntraLJPotential()
{}

/** Copy assignment operator */
IntraLJPotential& IntraLJPotential::operator=(const IntraLJPotential &other)
{
    LJPotential::operator=(other);
    return *this;
}

void IntraLJPotential::throwMissingForceComponent(const Symbol &symbol,
                              const IntraLJPotential::Components &components) const
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
void IntraLJPotential::assertCompatible(const IntraLJPotential::Molecule &mol,
                                const IntraLJPotential::Molecule &rest_of_mol) const
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
IntraLJPotential::Parameters 
IntraLJPotential::getParameters(const PartialMolecule &molecule,
                                const PropertyMap &map)
{
    need_update_ljpairs = true;

    return Parameters( AtomicParameters3D<LJParamID>(
                               molecule, map[parameters().coordinates()],
                               getLJParamIDs(molecule, map[parameters().lj()]) ),
              IntraScaledParameters<LJNBPairs>(
                               molecule, map[parameters().intraScaleFactors()] )
                     );
}

/** Update the parameters for the molecule going from 'old_molecule' to 
    'new_molecule', with the parameters found using the property map 'map'
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
IntraLJPotential::Parameters
IntraLJPotential::updateParameters(const IntraLJPotential::Parameters &old_params,
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
    const PropertyName &lj_property = map[parameters().lj()];
    const PropertyName &scl_property = map[parameters().intraScaleFactors()];
    
    //get what has changed
    bool new_coords = old_molecule.version(coords_property) !=
                         new_molecule.version(coords_property);
                             
    bool new_lj = ( old_molecule.version(lj_property) !=
                    new_molecule.version(lj_property) );

    bool new_scl = ( old_molecule.version(scl_property) !=
                         new_molecule.version(scl_property) );

    if (new_coords)
        new_params.setAtomicCoordinates( AtomicCoords3D(new_molecule, 
                                                        coords_property) );

    if (new_lj)
    {
        new_params.setAtomicParameters( getLJParamIDs(new_molecule, lj_property) );
        need_update_ljpairs = true;
    }

    if (new_scl)
        new_params.setIntraScaleFactors( 
                IntraScaledParameters<LJNBPairs>(new_molecule, scl_property) );

    return new_params;
}
                 
/** Update the parameters for the molecule going from 'old_molecule' to 
    'new_molecule', also while the parameters of 'old_molecule'
    where found in 'old_map', now get the parameters using 'new_map'
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
IntraLJPotential::Parameters
IntraLJPotential::updateParameters(const IntraLJPotential::Parameters &old_params,
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
    const PropertyName &old_lj = old_map[parameters().lj()];
    const PropertyName &old_scl = old_map[parameters().intraScaleFactors()];
    
    const PropertyName &new_coords = new_map[parameters().coordinates()];
    const PropertyName &new_lj = new_map[parameters().lj()];
    const PropertyName &new_scl = new_map[parameters().intraScaleFactors()];
    
    //get what has changed
    bool changed_coords = (new_coords != old_coords) or
                           old_molecule.version(old_coords) !=
                           new_molecule.version(old_coords);
                             
    bool changed_lj = (new_lj != old_lj) or
                      ( old_molecule.version(old_lj) !=
                        new_molecule.version(old_lj) );

    bool changed_scl = (new_scl != old_scl) or
                        old_molecule.version(old_scl) !=
                        new_molecule.version(old_scl);

    if (changed_coords)
        new_params.setAtomicCoordinates( AtomicCoords3D(new_molecule, new_coords) );

    if (changed_lj)
    {
        new_params.setAtomicParameters( getLJParamIDs(new_molecule, new_lj) );
        need_update_ljpairs = true;
    }

    if (changed_scl)
        new_params.setIntraScaleFactors( 
                        IntraScaledParameters<LJNBPairs>(new_molecule, new_scl) );

    return new_params;
}

/** Return the IntraLJPotential::Molecule representation of 'molecule',
    using the supplied PropertyMap to find the properties that contain
    the necessary forcefield parameters
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
IntraLJPotential::Molecule
IntraLJPotential::parameterise(const PartialMolecule &molecule,
                               const PropertyMap &map)
{
    return IntraLJPotential::Molecule(molecule, *this, map);
}

/** Concert the passed group of molecules into IntraLJPotential::Molecules,
    using the supplied PropertyMap to find the properties that contain
    the necessary forcefield parameters in each molecule
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
IntraLJPotential::Molecules 
IntraLJPotential::parameterise(const MoleculeGroup &molecules,
                               const PropertyMap &map)
{
    return IntraLJPotential::Molecules(molecules, *this, map);
}

void IntraLJPotential::calculateEnergy(const LJNBPairs::CGPairs &group_pairs, 
                            IntraLJPotential::EnergyWorkspace &distmat,
                            const IntraLJPotential::Parameter *params0_array, 
                            const IntraLJPotential::Parameter *params1_array,
                            const quint32 nats0, const quint32 nats1, 
                            double &iljnrg) const
{
    if (group_pairs.isEmpty())
    {
        //there are no scale factors between groups...
        for (quint32 i=0; i<nats0; ++i)
        {
            distmat.setOuterIndex(i);
            const Parameter &param0 = params0_array[i];
                
            if (param0.ljid == 0)
                continue;
                
            for (quint32 j=0; j<nats1; ++j)
            {
                const Parameter &param1 = params1_array[j];
                    
                if (param1.ljid != 0)
                {
                    iljnrg += calcLJEnergy(distmat[j], ljpairs.constData()[
                                                            ljpairs.map(param0.ljid,
                                                                    param1.ljid)]);
                }
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
                
            if (param0.ljid == 0)
                continue;
                
            for (quint32 j=0; j<nats1; ++j)
            {
                const LJScaleFactor &ljscl = group_pairs(i,j);

                const Parameter &param1 = params1_array[j];
                
                if (ljscl.lj() != 0 and param1.ljid != 0)
                {
                    iljnrg += ljscl.lj() * 
                                calcLJEnergy(distmat[j], ljpairs.constData()[
                                                            ljpairs.map(param0.ljid,
                                                                    param1.ljid)]);
                }
            }
        }
    }
}

void IntraLJPotential::calculateEnergy(const LJNBPairs::CGPairs &group_pairs, 
                            const QSet<Index> &atoms0, const QSet<Index> &atoms1,
                            IntraLJPotential::EnergyWorkspace &distmat,
                            const IntraLJPotential::Parameter *params0_array, 
                            const IntraLJPotential::Parameter *params1_array,
                            const quint32 nats0, const quint32 nats1, 
                            double &iljnrg) const
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

            if (param0.ljid == 0)
                continue;
                
            foreach (Index j, atoms1)
            {
                //do both coulomb and LJ
                const Parameter &param1 = params1_array[j];
                    
                if (param1.ljid != 0)
                {
                    iljnrg += calcLJEnergy(distmat[j], ljpairs.constData()[
                                                            ljpairs.map(param0.ljid,
                                                                    param1.ljid)]);
                }
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
                
            if (param0.ljid == 0)
                continue;
                
            foreach (Index j, atoms1)
            {
                //do both coulomb and LJ
                const LJScaleFactor &ljscl = group_pairs(i,j);

                const Parameter &param1 = params1_array[j];
                
                if (ljscl.lj() != 0 and param1.ljid != 0)
                {
                    iljnrg += ljscl.lj() * 
                                calcLJEnergy(distmat[j], ljpairs.constData()[
                                                            ljpairs.map(param0.ljid,
                                                                    param1.ljid)]);
                }
            }
        }
    }
}

/** Calculate the intramolecular LJ energy of the passed molecule, and
    add this onto 'energy'. This uses the passed workspace when
    performing the calculation */
void IntraLJPotential::calculateEnergy(const IntraLJPotential::Molecule &mol,
                                       IntraLJPotential::Energy &energy,
                                       IntraLJPotential::EnergyWorkspace &distmat,
                                       double scale_energy) const
{
    if (scale_energy == 0 or mol.isEmpty())
        return;

    const quint32 ngroups = mol.nCutGroups();
    
    const CoordGroup *groups_array = mol.coordinates().constData();
    
    BOOST_ASSERT( mol.parameters().atomicParameters().count() == int(ngroups) );
    const Parameters::Array *molparams_array 
                            = mol.parameters().atomicParameters().constData();

    const LJNBPairs &nbpairs = mol.parameters().intraScaleFactors();
    
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
            const double mindist = spce->calcDist2(group0, group1, distmat);
            
            if (mindist > switchfunc->cutoffDistance())
                //all of the atoms are definitely beyond cutoff
                continue;
                
            CGIdx cgidx_jgroup = mol.cgIdx(jgroup);
                
            //get the non-bonded scale factors for all pairs of atoms
            //between these groups (or within this group, if igroup == jgroup)
            const LJNBPairs::CGPairs &group_pairs = nbpairs(cgidx_igroup,
                                                            cgidx_jgroup);

            double iljnrg = 0;
            
            //loop over all intraatomic pairs and calculate the energies
            const quint32 nats1 = group1.count();
            const Parameter *params1_array = params1.constData();
            
            calculateEnergy(group_pairs, distmat, params0_array, params1_array,
                            nats0, nats1, iljnrg);
            
            //if this is the same group then half the energies to 
            //correct for double-counting
            if (igroup == jgroup)
            {
                iljnrg *= 0.5;
            }

            //now add these energies onto the total for the molecule,
            //scaled by any non-bonded feather factor
            if (mindist > switchfunc->featherDistance())
            {
                ljnrg += switchfunc->vdwScaleFactor( Length(mindist) ) * iljnrg;
            }
            else
            {
                ljnrg += iljnrg;
            }
        }
    }
    
    //add this molecule pair's energy onto the total
    energy += Energy(scale_energy * ljnrg);
}

/** Calculate the intramolecular LJ energy of the passed molecule
    interacting with the rest of the same molecule in 'rest_of_mol', and
    add this onto 'energy'. This uses the passed workspace when
    performing the calculation. Note that mol and rest_of_mol should
    not contain any overlapping atoms, and that they should both be
    part of the same molecule (albeit potentially at different versions,
    but with the same layout UID)
    
    \throw SireError::incompatible_error
*/
void IntraLJPotential::calculateEnergy(const IntraLJPotential::Molecule &mol,
                                       const IntraLJPotential::Molecule &rest_of_mol,
                                       IntraLJPotential::Energy &energy,
                                       IntraLJPotential::EnergyWorkspace &distmat,
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

    //the LJNBPairs must be the same in both molecules - this is checked
    //as part of assertCompatible(..)
    const LJNBPairs &nbpairs = mol.parameters().intraScaleFactors();
    
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
            const double mindist = spce->calcDist2(group0, group1, distmat);
            
            if (mindist > switchfunc->cutoffDistance())
                //all of the atoms are definitely beyond cutoff
                continue;
                
            //get the non-bonded scale factors for all pairs of atoms
            //between these groups (or within this group, if igroup == jgroup)
            const LJNBPairs::CGPairs &group_pairs = nbpairs(cgidx_igroup,
                                                            cgidx_jgroup);

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
                                nats0, nats1, iljnrg);
            }
            else
            {
                calculateEnergy(group_pairs, distmat,
                                params0_array, params1_array,
                                nats0, nats1, iljnrg);
            }

            //now add these energies onto the total for the molecule,
            //scaled by any non-bonded feather factor
            if (mindist > switchfunc->featherDistance())
            {
                ljnrg += switchfunc->vdwScaleFactor( Length(mindist) ) * iljnrg;
            }
            else
            {
                ljnrg += iljnrg;
            }
        }
    }

    //add the molecule's energy onto the total
    energy += Energy(scale_energy * ljnrg);
}

void IntraLJPotential::calculateLJForce(const LJNBPairs::CGPairs &group_pairs,
                             const CoordGroup &group0, const CoordGroup &group1,
                             const double mindist,
                             IntraLJPotential::ForceWorkspace &distmat,
                             const IntraLJPotential::Parameter *params0_array,
                             const IntraLJPotential::Parameter *params1_array,
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
                        const LJScaleFactor &ljscl = group_pairs(i,j);
                        
                        if (ljscl.lj() != 0 and param1.ljid != 0)
                        {
                            const double invdist = double(1) / distmat[j].length();
                            
                            const LJPair &ljpair = ljpairs.constData()[
                                                   ljpairs.map(param0.ljid,
                                                               param1.ljid)];
                
                            double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                            double sig_over_dist12 = pow_2(sig_over_dist6);

                            //calculate the energy
                            const double ljnrg = ljscl.lj() *
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
                        const LJScaleFactor &ljscl = group_pairs(i,j);
                        
                        if (ljscl.lj() != 0 and param1.ljid != 0)
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
                                  (ljscl.lj() *
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

void IntraLJPotential::calculateLJForce(const LJNBPairs::CGPairs &group_pairs,
                             const QSet<Index> &atoms0, const QSet<Index> &atoms1,
                             const CoordGroup &group0, const CoordGroup &group1,
                             const double mindist,
                             IntraLJPotential::ForceWorkspace &distmat,
                             const IntraLJPotential::Parameter *params0_array,
                             const IntraLJPotential::Parameter *params1_array,
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
                        const LJScaleFactor &ljscl = group_pairs(i,j);
                        
                        if (ljscl.lj() != 0 and param1.ljid != 0)
                        {
                            const double invdist = double(1) / distmat[j].length();
                            
                            const LJPair &ljpair = ljpairs.constData()[
                                                   ljpairs.map(param0.ljid,
                                                               param1.ljid)];
                
                            double sig_over_dist6 = pow_6(ljpair.sigma()*invdist);
                            double sig_over_dist12 = pow_2(sig_over_dist6);

                            //calculate the energy
                            const double ljnrg = ljscl.lj() *
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
                        const LJScaleFactor &ljscl = group_pairs(i,j);
                        
                        if (ljscl.lj() != 0 and param1.ljid != 0)
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
                                  (ljscl.lj() *
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
void IntraLJPotential::calculateLJForce(const IntraLJPotential::Molecule &mol,
                                         MolForceTable &forces, 
                                         IntraLJPotential::ForceWorkspace &distmat,
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

    const LJNBPairs &nbpairs = mol.parameters().intraScaleFactors();
    
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
            const LJNBPairs::CGPairs &group_pairs = nbpairs(cgidx_igroup,
                                                             cgidx_jgroup);

            //calculate the forces acting on group0 caused by group1
            calculateLJForce(group_pairs, group0, group1,
                             mindist, distmat, 
                             params0_array, params1_array,
                             nats0, nats1, group_forces0_array, scale_force);
            
        } // end of loop over CutGroups (jgroup)

    } // end of loop over CutGroups (igroup)
}

/** Calculate the LJ forces between the atoms in the molecule 'mol'
    and add these forces onto 'forces'. This uses
    the passed workspace to perform the calculation */
void IntraLJPotential::calculateForce(const IntraLJPotential::Molecule &mol,
                                      MolForceTable &forces, 
                                      IntraLJPotential::ForceWorkspace &distmat,
                                      double scale_force) const
{
    this->calculateLJForce(mol, forces, distmat, scale_force);
}

/** Calculate the LJ force acting on the part of the molecule
    in 'mol' caused by the rest of the molecule in 'rest_of_mol'. Note
    that these must both be of the same molecule, with the same
    layout UID and same nonbonded scale factors
    
    \throw SireError::incompatible_error
*/
void IntraLJPotential::calculateLJForce(const IntraLJPotential::Molecule &mol,
                                        const IntraLJPotential::Molecule &rest_of_mol,
                                        MolForceTable &forces,
                                        IntraLJPotential::ForceWorkspace &distmat,
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
    const LJNBPairs &nbpairs = mol.parameters().intraScaleFactors();

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
            const LJNBPairs::CGPairs &group_pairs = nbpairs(cgidx_igroup,
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

/** Calculate the LJ force acting on the part of the molecule
    in 'mol' caused by the rest of the molecule in 'rest_of_mol'. Note
    that these must both be of the same molecule, with the same
    layout UID and same nonbonded scale factors
    
    \throw SireError::incompatible_error
*/
void IntraLJPotential::calculateForce(const IntraLJPotential::Molecule &mol,
                                      const IntraLJPotential::Molecule &rest_of_mol,
                                      MolForceTable &forces,
                                      IntraLJPotential::ForceWorkspace &distmat,
                                      double scale_force) const
{
    this->calculateLJForce(mol, rest_of_mol, forces, distmat, scale_force);
}

void IntraLJPotential::calculateField(const IntraLJPotential::Molecule &mol, 
                    const LJProbe &probe,
                    MolFieldTable &fields,
                    IntraLJPotential::FieldWorkspace &workspace,
                    double scale_field) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code to calculate intramolecular LJ fields "
                "has yet to be written..."), CODELOC );
}

void IntraLJPotential::calculateField(const IntraLJPotential::Molecule &mol, 
                    const LJProbe &probe,
                    MolFieldTable &fields,
                    const Symbol &symbol,
                    const Components &components,
                    IntraLJPotential::FieldWorkspace &workspace,
                    double scale_field) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code to calculate intramolecular LJ fields "
                "has yet to be written..."), CODELOC );
}

void IntraLJPotential::calculateField(const IntraLJPotential::Molecule &mol0, 
                    const IntraLJPotential::Molecule &mol1,
                    const LJProbe &probe,
                    MolFieldTable &forces0,
                    IntraLJPotential::FieldWorkspace &workspace,
                    double scale_field) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code to calculate intramolecular LJ fields "
                "has yet to be written..."), CODELOC );
}

void IntraLJPotential::calculateField(const IntraLJPotential::Molecule &mol0,
                    const IntraLJPotential::Molecule &mol1,
                    const LJProbe &probe,
                    MolFieldTable &forces0,
                    const Symbol &symbol,
                    const Components &components,
                    IntraLJPotential::FieldWorkspace &workspace,
                    double scale_field) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code to calculate intramolecular LJ fields "
                "has yet to be written..."), CODELOC );
}

void IntraLJPotential::calculateField(const IntraLJPotential::Molecule &mol0,
                    const LJProbe &probe,
                    GridFieldTable &fields,
                    IntraLJPotential::FieldWorkspace &workspace,
                    double scale_field) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code to calculate intramolecular LJ fields "
                "has yet to be written..."), CODELOC );
}

void IntraLJPotential::calculateField(const IntraLJPotential::Molecule &mol0,
                    const LJProbe &probe,
                    GridFieldTable &fields,
                    const Symbol &symbol,
                    const Components &components,
                    IntraLJPotential::FieldWorkspace &workspace,
                    double scale_field) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code to calculate intramolecular LJ fields "
                "has yet to be written..."), CODELOC );
}

void IntraLJPotential::calculateLJField(const IntraLJPotential::Molecule &mol0, 
                      const IntraLJPotential::Molecule &mol1,
                      const LJProbe &probe,
                      MolFieldTable &fields0,
                      IntraLJPotential::ForceWorkspace &workspace,
                      double scale_field) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code to calculate intramolecular LJ fields "
                "has yet to be written..."), CODELOC );
}

void IntraLJPotential::calculateLJField(const IntraLJPotential::Molecule &mol0, 
                      const LJProbe &probe,
                      GridFieldTable &fields,
                      IntraLJPotential::ForceWorkspace &workspace,
                      double scale_field) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code to calculate intramolecular LJ fields "
                "has yet to be written..."), CODELOC );
}

void IntraLJPotential::calculatePotential(const IntraLJPotential::Molecule &mol, 
                        const LJProbe &probe,
                        MolPotentialTable &potentials,
                        IntraLJPotential::PotentialWorkspace &workspace,
                        double scale_potential) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code to calculate intramolecular LJ potentials "
                "has yet to be written..."), CODELOC );
}

void IntraLJPotential::calculatePotential(const IntraLJPotential::Molecule &mol, 
                        const LJProbe &probe,
                        MolPotentialTable &potentials,
                        const Symbol &symbol,
                        const Components &components,
                        IntraLJPotential::PotentialWorkspace &workspace,
                        double scale_potential) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code to calculate intramolecular LJ potentials "
                "has yet to be written..."), CODELOC );
}

void IntraLJPotential::calculatePotential(const IntraLJPotential::Molecule &mol0, 
                        const IntraLJPotential::Molecule &mol1,
                        const LJProbe &probe,
                        MolPotentialTable &pots0,
                        IntraLJPotential::PotentialWorkspace &workspace,
                        double scale_potential) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code to calculate intramolecular LJ potentials "
                "has yet to be written..."), CODELOC );
}

void IntraLJPotential::calculatePotential(const IntraLJPotential::Molecule &mol0,
                        const IntraLJPotential::Molecule &mol1,
                        const LJProbe &probe,
                        MolPotentialTable &pots0,
                        const Symbol &symbol,
                        const Components &components,
                        IntraLJPotential::PotentialWorkspace &workspace,
                        double scale_potential) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code to calculate intramolecular LJ potentials "
                "has yet to be written..."), CODELOC );
}

void IntraLJPotential::calculatePotential(const IntraLJPotential::Molecule &mol0,
                        const LJProbe &probe,
                        GridPotentialTable &pots,
                        IntraLJPotential::PotentialWorkspace &workspace,
                        double scale_potential) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code to calculate intramolecular LJ potentials "
                "has yet to be written..."), CODELOC );
}

void IntraLJPotential::calculatePotential(const IntraLJPotential::Molecule &mol0,
                        const LJProbe &probe,
                        GridPotentialTable &pots,
                        const Symbol &symbol,
                        const Components &components,
                        IntraLJPotential::PotentialWorkspace &workspace,
                        double scale_potential) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code to calculate intramolecular LJ potentials "
                "has yet to be written..."), CODELOC );
}

void IntraLJPotential::calculateLJPotential(const IntraLJPotential::Molecule &mol0, 
                          const IntraLJPotential::Molecule &mol1,
                          const LJProbe &probe,
                          MolPotentialTable &pots0,
                          IntraLJPotential::PotentialWorkspace &workspace,
                          double scale_potential) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code to calculate intramolecular LJ potentials "
                "has yet to be written..."), CODELOC );
}

void IntraLJPotential::calculateLJPotential(const IntraLJPotential::Molecule &mol0, 
                          const LJProbe &probe,
                          GridPotentialTable &pots,
                          IntraLJPotential::PotentialWorkspace &workspace,
                          double scale_potential) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code to calculate intramolecular LJ potentials "
                "has yet to be written..."), CODELOC );
}
