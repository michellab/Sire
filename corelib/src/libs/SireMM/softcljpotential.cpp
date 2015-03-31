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

#include "softcljpotential.h"
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

#ifdef SIRE_USE_SSE
    #ifdef __SSE__
        #include <emmintrin.h>   // CONDITIONAL_INCLUDE
    #else
        #undef SIRE_USE_SSE
    #endif
#endif

//#undef SIRE_USE_SSE

#include <QDebug>

using namespace SireMM;
using namespace SireMM::detail;

using namespace SireFF;
using namespace SireFF::detail;

using namespace SireMol;
using namespace SireVol;

using namespace SireMaths;

using namespace SireBase;

using namespace SireUnits;

using namespace SireStream;

///////
/////// Completely instantiate the SoftCLJPotential ancillary classes
///////

namespace SireFF
{
    namespace detail
    {
        template
        class FFMolecule3D<InterSoftCLJPotential>;

        template
        class FFMolecules3D<InterSoftCLJPotential>;

        template
        class ChangedMolecule<InterSoftCLJPotential::Molecule>;
    }
}

/////////////
///////////// Implementation of SoftCLJPotential
/////////////

static const RegisterMetaType<SoftCLJPotential> r_cljpot( MAGIC_ONLY, NO_ROOT,
                                                          "SireMM::SoftCLJPotential" );

/** Serialise to a binary datastream */
QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds,
                                      const SoftCLJPotential &cljpot)
{
    writeHeader(ds, r_cljpot, 2);
    
    ds << static_cast<const CLJPotential&>(cljpot);

    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds,
                                      SoftCLJPotential &cljpot)
{
    VersionID v = readHeader(ds, r_cljpot);
    
    if (v == 1 or v == 2)
    {
        ds >> static_cast<CLJPotential&>(cljpot);
    
        //extract all of the properties
        QVector<double> alpha_values;
        QVector<qint32> alpha_index(MAX_ALPHA_VALUES, -1);
        
        for (int i=0; i<MAX_ALPHA_VALUES; ++i)
        {
            QString propname = QString( "alpha%1" ).arg(i);
            
            if (cljpot.props.hasProperty(propname))
            {
                const Property &prop = cljpot.props.property(propname);
                
                if (prop.isA<NullProperty>())
                    continue;
            
                double alpha = prop.asA<VariantProperty>().convertTo<double>();
                                     
                int idx = alpha_values.indexOf(alpha);
                
                if (idx == -1)
                {
                    idx = alpha_values.count();
                    alpha_values.append(alpha);
                }
                
                alpha_index[i] = idx;
            }
        }
        
        if (cljpot.props.hasProperty("alpha"))
        {
            const Property &prop = cljpot.props.property("alpha");
            
            if (not prop.isA<NullProperty>())
            {
                double alpha = prop.asA<VariantProperty>().convertTo<double>();
                                 
                if (alpha_values.count() == 0)
                {
                    alpha_values.append(alpha);
                    alpha_index[0] = 0;
                }
                else if (alpha_values.count() > 1)
                {
                    throw SireError::program_bug( QObject::tr(
                        "How can we have a single value of alpha (%1) "
                        "when multiple values (%2) have been read??")
                            .arg(alpha)
                            .arg( Sire::toString(alpha_values) ),
                                CODELOC );
                }
                else if (alpha_values.at(0) != alpha)
                {
                    throw SireError::program_bug( QObject::tr(
                        "How can the value of alpha (%1) be different to "
                        "the only alpha component value that has been read? (%2)")
                            .arg(alpha).arg(alpha_values.at(0)), CODELOC );
                }
            }
        }
    
        cljpot.alpha_values = alpha_values;
        cljpot.alpha_index = alpha_index;
        
        cljpot.clearOrphanedAlpha();
        cljpot.rebuildAlphaProperties();
    
        cljpot.shift_delta = cljpot.props.property("shiftDelta")
                                  .asA<VariantProperty>().convertTo<double>();
                
        cljpot.coul_power = cljpot.props.property("coulombPower")
                                  .asA<VariantProperty>().convertTo<quint32>();

        cljpot.lj_power = cljpot.props.property("ljPower")
                                  .asA<VariantProperty>().convertTo<quint32>();
    }
    else 
        throw version_error(v, "1,2", r_cljpot, CODELOC);
    
    return ds;
}

/** Constructor */
SoftCLJPotential::SoftCLJPotential()
                 : CLJPotential(), shift_delta(1.0 * angstrom),
                                   coul_power(1), lj_power(1)
{
    //construct the default alpha values
    alpha_values = QVector<double>(1, 0.0);
    alpha_index = QVector<qint32>(MAX_ALPHA_VALUES, -1);
    alpha_index[0] = 0;
    
    this->rebuildAlphaProperties();

    //record the defaults
    props.setProperty( "shiftDelta", VariantProperty(shift_delta) );
    props.setProperty( "coulombPower", VariantProperty(coul_power) );
    props.setProperty( "ljPower", VariantProperty(lj_power) );
}

/** Copy constructor */
SoftCLJPotential::SoftCLJPotential(const SoftCLJPotential &other)
                 : CLJPotential(other),
                   alpha_values(other.alpha_values), 
                   alpha_index(other.alpha_index),
                   shift_delta(other.shift_delta),
                   coul_power(other.coul_power), lj_power(other.lj_power)
{}

/** Destructor */
SoftCLJPotential::~SoftCLJPotential()
{}

/** Copy assignment operator */
SoftCLJPotential& SoftCLJPotential::operator=(const SoftCLJPotential &other)
{
    if (this != &other)
    {
        CLJPotential::operator=(other);
        alpha_values = other.alpha_values;
        alpha_index = other.alpha_index;
        shift_delta = other.shift_delta;
        coul_power = other.coul_power;
        lj_power = other.lj_power;
    }
    
    return *this;
}

/** Return the value of alpha - this will raise an exception if this
    potential calculates multiple alpha values at the same time
    
    \throw SireError::invalid_state
*/
double SoftCLJPotential::alpha() const
{
    if (alpha_values.count() != 1)
        throw SireError::invalid_state( QObject::tr(
                "This potential has multiple alpha values (%1), so it "
                "is not possible to return just one.")
                    .arg(Sire::toString(alpha_values)), CODELOC );
                    
    return alpha_values.at(0);
}

/** Return whether or not the ith alpha component has 
    an alpha value
    
    \throw SireError::invalid_index
*/
bool SoftCLJPotential::hasAlphaValue(int i) const
{
    return alpha_index.at( Index(i).map(alpha_index.count()) ) != -1;
}

/** Return the value of alpha for the ith alpha component

    \throw SireError::invalid_index
    \throw SireError::invalid_state
*/
double SoftCLJPotential::alpha(int i) const
{
    if (not this->hasAlphaValue(i))
        throw SireError::invalid_state( QObject::tr(
            "There is no alpha value associated with the alpha component "
            "at index %1.").arg(i), CODELOC );
            
    return alpha_values.at( alpha_index.at(Index(i).map(alpha_index.count())) );
}

/** Return the number of active alpha components (the number of 
    alpha components that have a value of alpha - those that don't
    will not contribute to this potential) */
int SoftCLJPotential::nActiveAlphaComponents() const
{
    int n = 0;
    
    for (int i=0; i<alpha_index.count(); ++i)
    {
        if (alpha_index.at(i) != -1)
            ++n;
    }
    
    return n;
}

/** Internal function used to remove orphaned alpha values */
void SoftCLJPotential::clearOrphanedAlpha()
{
    for (qint32 i=0; i<alpha_values.count(); ++i)
    {
        if (not alpha_index.contains(i))
        {
            //this is an orphaned alpha value - remove it
            alpha_values.remove(i);
            
            //update the alpha index
            for (int j=0; j<alpha_index.count(); ++j)
            {
                if (alpha_index.at(j) > i)
                    alpha_index[j] -= 1;
            }
            
            //use recursion to find other missing alpha values
            this->clearOrphanedAlpha();
            return;
        } 
    }
}

/** Internal function used to rebuild the properties representing the 
    alpha components */
void SoftCLJPotential::rebuildAlphaProperties()
{
    if (alpha_values.count() == 1)
    {
        props.setProperty("alpha", VariantProperty(alpha_values.at(0)));
    }
    else
    {
        props.setProperty("alpha", NullProperty());
    }
    
    for (int i=0; i<MAX_ALPHA_VALUES; ++i)
    {
        int idx = alpha_index.at(i);
        
        if (idx == -1)
        {
            props.setProperty( QString("alpha%1").arg(i), NullProperty() );
        }
        else
        {
            props.setProperty( QString("alpha%1").arg(i), 
                               VariantProperty(alpha_values.at(idx)) );
        }
    }
}

/** Set the value of alpha to 'alpha' - this clears all of the alpha
    values of all of the components, and sets only the first component
    to use this value of alpha. This returns whether or not this
    changes the forcefield */
bool SoftCLJPotential::setAlpha(double alpha)
{
    if (alpha_values.count() == 1)
    {
        if (alpha_values.at(0) == alpha)
            return false;
    }
    
    alpha_values = QVector<double>(1, alpha);
    
    alpha_index = QVector<qint32>(MAX_ALPHA_VALUES, -1);
    alpha_index[0] = 0;
    
    this->rebuildAlphaProperties();
    this->changedPotential();
    
    return true;
}

/** Set the value of alpha for the ith alpha component to 'alpha'

    \throw SireError::invalid_index
*/
bool SoftCLJPotential::setAlpha(int i, double alpha)
{
    i = Index(i).map( MAX_ALPHA_VALUES );
    
    int idx = alpha_values.indexOf(alpha);
    
    if (idx == -1)
    {
        idx = alpha_values.count();
        alpha_values.append(alpha);
    }
    
    if (alpha_index[i] != idx)
    {
        alpha_index[i] = idx;
        this->clearOrphanedAlpha();
        this->rebuildAlphaProperties();
        this->changedPotential();
        
        return true;
    }
    else
        return false;
}

/** Remove the value of alpha for the ith alpha component

    \throw SireError::invalid_index
*/
bool SoftCLJPotential::removeAlpha(int i)
{
    i = Index(i).map(MAX_ALPHA_VALUES);
    
    if (alpha_index.at(i) != -1)
    {
        alpha_index[i] = -1;
        this->clearOrphanedAlpha();
        this->rebuildAlphaProperties();
        this->changedPotential();
        
        return true;
    }
    else
        return false;
}

/** Clear all of the alpha values - this removes the values of alpha
    for all of the alpha components, and then sets the alpha value
    of the first component to 0 (as we need to have at least one 
    alpha value) */
void SoftCLJPotential::clearAlphas()
{
    QVector<double> new_alpha_values(1, 0.0);
    QVector<qint32> new_alpha_index(MAX_ALPHA_VALUES, -1);
    
    new_alpha_index[0] = 0;
        
    if (alpha_values != new_alpha_values or 
        alpha_index != new_alpha_index)
    {
        alpha_values = new_alpha_values;
        alpha_index = new_alpha_index;
        
        this->rebuildAlphaProperties();
        this->changedPotential();
    }
}

/** Set the value of the delta value used in the LJ shift function */
bool SoftCLJPotential::setShiftDelta(double new_delta)
{
    if (shift_delta != new_delta)
    {
        shift_delta = new_delta;
        props.setProperty("shiftDelta", VariantProperty(shift_delta));
        this->changedPotential();
        return true;
    }
    else
        return false;
}

/** Set the coulomb power, which is used to control how strongly
    the coulomb terms are softened */
bool SoftCLJPotential::setCoulombPower(int power)
{
    if (power < 0)
        throw SireError::invalid_arg( QObject::tr(
            "You cannot use a negative Coulomb power (%1)").arg(power),
                CODELOC );
                
    quint32 new_power(power);
    
    if (coul_power != new_power)
    {
        coul_power = new_power;
        props.setProperty("coulombPower", VariantProperty(coul_power));
        this->changedPotential();
        return true;
    }
    else
        return false;
}

/** Set the LJ power, which is used to control how strongly
    the LJ terms are softened */
bool SoftCLJPotential::setLJPower(int power)
{
    if (power < 0)
        throw SireError::invalid_arg( QObject::tr(
            "You cannot use a negative LJ power (%1)").arg(power),
                CODELOC );
                
    quint32 new_power(power);
    
    if (lj_power != new_power)
    {
        lj_power = new_power;
        props.setProperty("ljPower", VariantProperty(lj_power));
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
bool SoftCLJPotential::setProperty(const QString &name, const Property &value)
{
    if (name == QLatin1String("alpha"))
    {
        return this->setAlpha( value.asA<VariantProperty>().convertTo<double>() );
    }
    else if (name == QLatin1String("shiftDelta"))
    {
        return this->setShiftDelta( value.asA<VariantProperty>().convertTo<double>() );
    }
    else if (name == QLatin1String("coulombPower"))
    {
        return this->setCoulombPower( value.asA<VariantProperty>().convertTo<int>() );
    }
    else if (name == QLatin1String("ljPower"))
    {
        return this->setLJPower( value.asA<VariantProperty>().convertTo<int>() );
    }
    else
    {
        //is this one of the alpha properties?
        for (int i=0; i<alpha_index.count(); ++i)
        {
            QString propname = QString("alpha%1").arg(i);
            
            if (name == propname)
            {
                return this->setAlpha(i, value.asA<VariantProperty>()
                                              .convertTo<double>() );
            }
        }
    
        //no it's not - see if the CLJPotential recognises this property
        return CLJPotential::setProperty(name, value);
    }
}

/** Return the delta value used in the LJ shifting function */
double SoftCLJPotential::shiftDelta() const
{
    return shift_delta;
}

/** Return the coulomb power */
int SoftCLJPotential::coulombPower() const
{
    return coul_power;
}

/** Return the LJ power */
int SoftCLJPotential::ljPower() const
{
    return lj_power;
}

/////////////
///////////// Implementation of InterSoftCLJPotential
/////////////

static const RegisterMetaType<InterSoftCLJPotential> r_interclj( MAGIC_ONLY, NO_ROOT,
                                                   InterSoftCLJPotential::typeName() );

/** Serialise to a binary datastream */
QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds,
                                      const InterSoftCLJPotential &interclj)
{
    writeHeader(ds, r_interclj, 1);
    
    ds << static_cast<const SoftCLJPotential&>(interclj);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds,
                                      InterSoftCLJPotential &interclj)
{
    VersionID v = readHeader(ds, r_interclj);
    
    if (v == 1)
    {
        ds >> static_cast<SoftCLJPotential&>(interclj);
    }
    else
        throw version_error(v, "1", r_interclj, CODELOC);
        
    return ds;
}

/** Constructor */
InterSoftCLJPotential::InterSoftCLJPotential() : SoftCLJPotential()
{}

/** Copy constructor */
InterSoftCLJPotential::InterSoftCLJPotential(const InterSoftCLJPotential &other)
                      : SoftCLJPotential(other)
{}

/** Destructor */
InterSoftCLJPotential::~InterSoftCLJPotential()
{}

/** Copy assignment operator */
InterSoftCLJPotential& 
InterSoftCLJPotential::operator=(const InterSoftCLJPotential &other)
{
    SoftCLJPotential::operator=(other);
    return *this;
}

void InterSoftCLJPotential::throwMissingForceComponent(const Symbol &symbol,
                              const InterSoftCLJPotential::Components &components) const
{
    throw SireFF::missing_component( QObject::tr(
        "There is no force component in potential %1 - available "
        "components are %2, %3 and %4.")
            .arg(this->what())
            .arg(components.total().toString(), components.coulomb().toString(),
                 components.lj().toString()), CODELOC );
}

void InterSoftCLJPotential::throwMissingFieldComponent(const Symbol &symbol,
                              const InterSoftCLJPotential::Components &components) const
{
    throw SireFF::missing_component( QObject::tr(
        "There is no field component in potential %1 - available "
        "components are %2, %3 and %4.")
            .arg(this->what())
            .arg(components.total().toString(), components.coulomb().toString(),
                 components.lj().toString()), CODELOC );
}

void InterSoftCLJPotential::throwMissingPotentialComponent(const Symbol &symbol,
                              const InterSoftCLJPotential::Components &components) const
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
InterSoftCLJPotential::Parameters 
InterSoftCLJPotential::getParameters(const PartialMolecule &molecule,
                                     const PropertyMap &map)
{
    return Parameters( molecule, map[parameters().coordinates()],
                       CLJPotential::getCLJParameters(molecule, 
                                                      map[parameters().charge()],
                                                      map[parameters().lj()]) );
}

/** Update the parameters for the molecule going from 'old_molecule' to 
    'new_molecule', with the parameters found using the property map 'map'
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
InterSoftCLJPotential::Parameters
InterSoftCLJPotential::updateParameters(
                                    const InterSoftCLJPotential::Parameters &old_params,
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
        new_params.setAtomicParameters( CLJPotential::getCLJParameters(new_molecule,
                                                         chg_property, lj_property) );
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
InterSoftCLJPotential::Parameters
InterSoftCLJPotential::updateParameters(
                                    const InterSoftCLJPotential::Parameters &old_params,
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
        new_params.setAtomicParameters( CLJPotential::getCLJParameters(new_molecule,
                                                                       new_chg, new_lj) );

    return new_params;
}

/** Return the InterSoftCLJPotential::Molecule representation of 'molecule',
    using the supplied PropertyMap to find the properties that contain
    the necessary forcefield parameters
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
InterSoftCLJPotential::Molecule
InterSoftCLJPotential::parameterise(const PartialMolecule &molecule,
                                    const PropertyMap &map)
{
    return InterSoftCLJPotential::Molecule(molecule, *this, map);
}

/** Convert the passed group of molecules into InterSoftCLJPotential::Molecules,
    using the supplied PropertyMap to find the properties that contain
    the necessary forcefield parameters in each molecule
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
InterSoftCLJPotential::Molecules 
InterSoftCLJPotential::parameterise(const MoleculeGroup &molecules,
                                    const PropertyMap &map)
{
    return InterSoftCLJPotential::Molecules(molecules, *this, map);
}

/** Return the total charge of the parameters for the group in 'params' */
double InterSoftCLJPotential::totalCharge(
                            const InterSoftCLJPotential::Parameters::Array &params) const
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
void InterSoftCLJPotential::_pvt_calculateEnergy(
                                          const InterSoftCLJPotential::Molecule &mol0,
                                          const InterSoftCLJPotential::Molecule &mol1,
                                          InterSoftCLJPotential::Energy &energy,
                                          InterSoftCLJPotential::EnergyWorkspace &distmat,
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

    //the alpha_values array contains all of the unique alpha values
    const double *alfa = alpha_values.constData();
    const int nalpha = alpha_values.count();
    
    if (nalpha <= 0)
        return;
    
    double cnrg[nalpha];
    double ljnrg[nalpha];
    
    for (int i=0; i<nalpha; ++i)
    {
        cnrg[i] = 0;
        ljnrg[i] = 0;
    }

    const double Rcoul = qMax(1e-5,qMin(1e9,
                            switchfunc->electrostaticCutoffDistance().to(angstrom)));
    const double Rlj = qMax(1e-5,qMin(1e9, switchfunc->vdwCutoffDistance().to(angstrom)) );
    const double Rc = qMax(Rcoul,Rlj);
    const double Rlj2 = Rlj*Rlj;
    
    //this uses the following potentials
    //           Zacharias and McCammon, J. Chem. Phys., 1994, and also,
    //           Michel et al., JCTC, 2007
    //
    //  V_{LJ}(r) = 4 epsilon [ ( sigma^12 / (delta*sigma + r^2)^6 ) - 
    //                          ( sigma^6  / (delta*sigma + r^2)^3 ) ]
    //
    //  delta = shift_delta * alpha
    //
    //  V_{coul}(r) = (1-alpha)^n q_i q_j / 4 pi eps_0 (alpha+r^2)^(1/2)
    //
    // If force shifting is used, then I have modified this to fit the
    // force-shifted coulomb equation by substituting
    //
    // sr = (alpha+r^2)^(1/2)
    // sRc = (alpha+Rc^2)^(1/2)
    //
    // into
    //
    // V_{coul}{sr} = (1-alpha)^n q_i q_j [ 1/sr - 1/sRc + 1/sRc^2 (sr - sRc) ]
    //
    // This contrasts to Rich T's LJ softcore function, which was;
    //
    //  V_{LJ}(r) = 4 epsilon [ (sigma^12 / (alpha^m sigma^6 + r^6)^2) - 
    //                          (sigma^6  / (alpha^m sigma^6 + r^6) ) ]

    double one_minus_alfa_to_n[nalpha];
    double delta[nalpha];
    
    for (int i=0; i<nalpha; ++i)
    {
        one_minus_alfa_to_n[i] = SireMaths::pow(1 - alfa[i], int(coul_power));
        delta[i] = shift_delta * alfa[i];
    }
    
    #ifdef SIRE_TIME_ROUTINES
    int nflops = 0;
    #endif

    if (use_electrostatic_shifting)
    {
        double sRcoul[nalpha];
        double one_over_sRcoul[nalpha];
        double one_over_sRcoul2[nalpha];
        
        for (int i=0; i<nalpha; ++i)
        {
            sRcoul[i] = std::sqrt(alfa[i] + Rcoul*Rcoul);
            one_over_sRcoul[i] = double(1) / sRcoul[i];
            one_over_sRcoul2[i] = double(1) / (sRcoul[i]*sRcoul[i]);
        }
    
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
                                            spce->beyond(Rc,aabox0, group1.aaBox());
                
                if (not within_cutoff)
                    //this CutGroup is either the cutoff distance
                    continue;
                
                //calculate all of the interatomic distances^2
                const double mindist = spce->calcDist2(group0, group1, distmat);
                
                if (mindist > Rc)
                {
                    //all of the atoms are definitely beyond cutoff
                    continue;
                }
                   
                double icnrg[nalpha];
                double iljnrg[nalpha];
                
                for (int i=0; i<nalpha; ++i)
                {
                    icnrg[i] = 0;
                    iljnrg[i] = 0;
                }
                
                //loop over all interatomic pairs and calculate the energies
                const quint32 nats1 = group1.count();
                const Parameter *params1_array = params1.constData();

                #ifdef SIRE_USE_SSE
                {
                    const int remainder = nats1 % 2;
                    
                    const __m128d sse_one = { 1.0, 1.0 };
                    
                    __m128d sse_cnrg[nalpha];
                    __m128d sse_ljnrg[nalpha];
                    __m128d sse_alpha[nalpha];
                    __m128d sse_delta[nalpha];
                    
                    __m128d sse_Rlj2 = _mm_set1_pd(Rlj2);
                    __m128d sse_sRcoul[nalpha];
                    __m128d sse_one_over_sRcoul[nalpha];
                    __m128d sse_one_over_sRcoul2[nalpha];
                    
                    for (int i=0; i<nalpha; ++i)
                    {
                        sse_cnrg[i] = _mm_set1_pd(0);
                        sse_ljnrg[i] = _mm_set1_pd(0);
                        
                        sse_alpha[i] = _mm_set1_pd(alfa[i]);
                        sse_delta[i] = _mm_set1_pd(delta[i]);
                        
                        sse_sRcoul[i] = _mm_set1_pd(sRcoul[i]);
                        sse_one_over_sRcoul[i] = _mm_set1_pd(one_over_sRcoul[i]);
                        sse_one_over_sRcoul2[i] = _mm_set1_pd(one_over_sRcoul2[i]);
                    }
                    
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
                            
                            const __m128d sse_r2 = _mm_set_pd(distmat[j], distmat[j+1]);
                                                       
                            const LJPair &ljpair0 = ljpairs.constData()[
                                                    ljpairs.map(param0.ljid,
                                                                param10.ljid)];
                        
                            const LJPair &ljpair1 = ljpairs.constData()[
                                                    ljpairs.map(param0.ljid,
                                                                param11.ljid)];
                        
                            __m128d sse_sig = _mm_set_pd(ljpair0.sigma(), 
                                                         ljpair1.sigma());
                            __m128d sse_eps = _mm_set_pd(ljpair0.epsilon(), 
                                                         ljpair1.epsilon());

                            const __m128d sse_sig2 = _mm_mul_pd(sse_sig, sse_sig);
                            const __m128d sse_sig3 = _mm_mul_pd(sse_sig2, sse_sig);
                            const __m128d sse_sig6 = _mm_mul_pd(sse_sig3, sse_sig3);
                            
                            for (int k=0; k<nalpha; ++k)
                            {
                                //coulomb energy
                                {
                                    __m128d sse_sr = _mm_sqrt_pd( _mm_add_pd(sse_r2,sse_alpha[k]) );
                                    __m128d sse_one_over_sr = _mm_div_pd(sse_one, sse_sr);
                            
                                    __m128d nrg = _mm_sub_pd(sse_sr, sse_sRcoul[k]);
                                    nrg = _mm_mul_pd(nrg, sse_one_over_sRcoul2[k]);
                                    nrg = _mm_add_pd(nrg, sse_one_over_sr);
                                    nrg = _mm_sub_pd(nrg, sse_one_over_sRcoul[k]);

                                    __m128d sse_chg = _mm_set_pd( param10.reduced_charge,
                                                                  param11.reduced_charge );
                        
                                    sse_chg = _mm_mul_pd(sse_chg, sse_chg0);
                        
                                    nrg = _mm_mul_pd(sse_chg, nrg);

                                    const __m128d sse_in_cutoff = _mm_cmplt_pd(sse_sr,
                                                                               sse_sRcoul[k]);
                                    nrg = _mm_and_pd(nrg, sse_in_cutoff);
                                
                                    sse_cnrg[k] = _mm_add_pd(sse_cnrg[k], nrg);
                                }
                                
                                //lj energy
                                {
                                    const __m128d sse_in_cutoff = _mm_cmplt_pd(sse_r2, sse_Rlj2);
                                
                                    //calculate shift = alpha * sigma * shift_delta
                                    const __m128d sse_shift = _mm_mul_pd(sse_sig, sse_delta[k]);

                                    __m128d lj_denom = _mm_add_pd(sse_r2, sse_shift);
                                    __m128d lj_denom2 = _mm_mul_pd(lj_denom, lj_denom);
                                    lj_denom = _mm_mul_pd(lj_denom, lj_denom2);
                            
                                    const __m128d sig6_over_denom = _mm_div_pd(sse_sig6,
                                                                               lj_denom);
                                                                           
                                    const __m128d sig12_over_denom2 = _mm_mul_pd(sig6_over_denom,
                                                                                 sig6_over_denom);
                                                  
                                    //calculate LJ energy (the factor of 4 is added later)
                                    __m128d nrg = _mm_sub_pd(sig12_over_denom2, sig6_over_denom);
                                                         
                                    nrg = _mm_mul_pd(sse_eps, nrg);
                                    nrg = _mm_and_pd(nrg, sse_in_cutoff);
                                    sse_ljnrg[k] = _mm_add_pd(sse_ljnrg[k], nrg);
                                }
                            }
                        }
                              
                        if (remainder == 1)
                        {
                            const Parameter &param1 = params1_array[nats1-1];

                            const double r2 = distmat[nats1-1];
                            
                            const double q2 = param0.reduced_charge * param1.reduced_charge;

                            const LJPair &ljpair = ljpairs.constData()[
                                                    ljpairs.map(param0.ljid,
                                                                param1.ljid)];
                            
                            const double sig2 = ljpair.sigma() * ljpair.sigma();
                            const double sig6 = sig2 * sig2 * sig2;
                            
                            for (int k=0; k<nalpha; ++k)
                            {
                                //coulomb energy
                                {
                                    const double sr = std::sqrt(alfa[k] + r2);

                                    if (sr < sRcoul[k])
                                    {
                                        const double one_over_sr = double(1) / sr;
                                
                                        icnrg[k] += q2 *
                                                (one_over_sr - one_over_sRcoul[k] +
                                                 one_over_sRcoul2[k]*(sr-sRcoul[k]));
                                    }
                                }
                                
                                //lj energy
                                if (r2 < Rlj2)
                                {
                                    const double shift = ljpair.sigma() * delta[k];
                                
                                    double lj_denom = r2 + shift;
                                    lj_denom = lj_denom * lj_denom * lj_denom;
                                
                                    const double sig6_over_denom = sig6 / lj_denom;
                                    const double sig12_over_denom2 = sig6_over_denom * 
                                                                     sig6_over_denom;
                                
                                    iljnrg[k] += ljpair.epsilon() * (sig12_over_denom2 - 
                                                                     sig6_over_denom);
                                }
                            }
                        }
                    }
                    
                    for (int k=0; k<nalpha; ++k)
                    {
                        icnrg[k] += *((const double*)&(sse_cnrg[k])) +
                                    *( ((const double*)&(sse_cnrg[k])) + 1 );
                             
                        iljnrg[k] += *((const double*)&(sse_ljnrg[k])) +
                                     *( ((const double*)&(sse_ljnrg[k])) + 1 );
                    }
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

                            const double r2 = distmat[j];
                            
                            const double q2 = param0.reduced_charge * param1.reduced_charge;

                            const LJPair &ljpair = ljpairs.constData()[
                                                    ljpairs.map(param0.ljid,
                                                                param1.ljid)];
                            
                            const double sig2 = ljpair.sigma() * ljpair.sigma();
                            const double sig6 = sig2 * sig2 * sig2;
                            
                            for (int k=0; k<nalpha; ++k)
                            {
                                //coulomb energy
                                {
                                    const double sr = std::sqrt(alfa[k] + r2);

                                    if (sr < sRcoul[k])
                                    {
                                        const double one_over_sr = double(1) / sr;

                                        icnrg[k] += q2 * (one_over_sr - one_over_sRcoul[k] +
                                                          one_over_sRcoul2[k]*(sr-sRcoul[k]));
                                    }
                                }
                                
                                //lj energy
                                if (r2 < Rlj2)
                                {
                                    const double shift = ljpair.sigma() * delta[k];
                                
                                    double lj_denom = r2 + shift;
                                    lj_denom = lj_denom * lj_denom * lj_denom;
                                
                                    const double sig6_over_denom = sig6 / lj_denom;
                                    const double sig12_over_denom2 = sig6_over_denom *
                                                                     sig6_over_denom;
            
                                    iljnrg[k] += ljpair.epsilon() * (sig12_over_denom2 - 
                                                                     sig6_over_denom);
                                }
                            }
                        }
                    }
                }
                #endif
                
                //now add these energies onto the total for the molecule,
                //scaled by any non-bonded feather factor
                for (int i=0; i<nalpha; ++i)
                {
                    cnrg[i] += icnrg[i];
                    ljnrg[i] += iljnrg[i];
                }
            }
        }
    }
    else if (use_reaction_field)
    {
        //use the reaction field potential
        // E = (q1 q2 / 4 pi eps_0) * ( 1/r + k r^2 - c )
        // where k = (1 / r_c^3) * (eps - 1)/(2 eps + 1)
        // c = (1/r_c) * (3 eps)/(2 eps + 1)
        double sRcoul[nalpha];
        double k_rf[nalpha];
        double c_rf[nalpha];
        
        for (int i=0; i<nalpha; ++i)
        {
            sRcoul[i] = std::sqrt(alfa[i] + Rcoul*Rcoul);
            k_rf[i] = (1.0 / pow_3(sRcoul[i])) * ( (rf_dielectric_constant-1) /
                                                   (2*rf_dielectric_constant + 1) );
            c_rf[i] = (1.0 / sRcoul[i]) * ( (3*rf_dielectric_constant) /
                                            (2*rf_dielectric_constant + 1) );
        }
    
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
                                            spce->beyond(Rc,aabox0, group1.aaBox());
                
                if (not within_cutoff)
                    //this CutGroup is either the cutoff distance
                    continue;
                
                //calculate all of the interatomic distances^2
                const double mindist = spce->calcDist2(group0, group1, distmat);
                
                if (mindist > Rc)
                {
                    //all of the atoms are definitely beyond cutoff
                    continue;
                }
                   
                double icnrg[nalpha];
                double iljnrg[nalpha];
                
                for (int i=0; i<nalpha; ++i)
                {
                    icnrg[i] = 0;
                    iljnrg[i] = 0;
                }
                
                //loop over all interatomic pairs and calculate the energies
                const quint32 nats1 = group1.count();
                const Parameter *params1_array = params1.constData();

                #ifdef SIRE_USE_SSE
                {
                    const int remainder = nats1 % 2;
                    
                    const __m128d sse_one = { 1.0, 1.0 };
                    
                    __m128d sse_cnrg[nalpha];
                    __m128d sse_ljnrg[nalpha];
                    __m128d sse_alpha[nalpha];
                    __m128d sse_delta[nalpha];
                    
                    const __m128d sse_Rlj2 = _mm_set1_pd(Rlj2);
                    __m128d sse_sRcoul[nalpha];
                    __m128d sse_k_rf[nalpha];
                    __m128d sse_c_rf[nalpha];
                    
                    for (int i=0; i<nalpha; ++i)
                    {
                        sse_cnrg[i] = _mm_set_pd(0, 0);
                        sse_ljnrg[i] = _mm_set_pd(0, 0);
                        
                        sse_alpha[i] = _mm_set_pd(alfa[i], alfa[i]);
                        sse_delta[i] = _mm_set_pd(delta[i], delta[i]);
                        
                        sse_sRcoul[i] = _mm_set1_pd(sRcoul[i]);
                        sse_k_rf[i] = _mm_set1_pd(k_rf[i]);
                        sse_c_rf[i] = _mm_set1_pd(c_rf[i]);
                    }
                    
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
                            
                            const __m128d sse_r2 = _mm_set_pd(distmat[j], distmat[j+1]);

                                                       
                            const LJPair &ljpair0 = ljpairs.constData()[
                                                    ljpairs.map(param0.ljid,
                                                                param10.ljid)];
                        
                            const LJPair &ljpair1 = ljpairs.constData()[
                                                    ljpairs.map(param0.ljid,
                                                                param11.ljid)];
                        
                            __m128d sse_sig = _mm_set_pd(ljpair0.sigma(), 
                                                         ljpair1.sigma());
                            __m128d sse_eps = _mm_set_pd(ljpair0.epsilon(), 
                                                         ljpair1.epsilon());

                            const __m128d sse_sig2 = _mm_mul_pd(sse_sig, sse_sig);
                            const __m128d sse_sig3 = _mm_mul_pd(sse_sig2, sse_sig);
                            const __m128d sse_sig6 = _mm_mul_pd(sse_sig3, sse_sig3);
                            
                            for (int k=0; k<nalpha; ++k)
                            {
                                //coulomb energy
                                {
                                    __m128d sse_sr = _mm_sqrt_pd( _mm_add_pd(sse_r2,sse_alpha[k]) );
                                    const __m128d sse_in_cutoff = _mm_cmplt_pd(sse_sr,
                                                                               sse_sRcoul[k]);

                                    __m128d sse_one_over_sr = _mm_div_pd(sse_one, sse_sr);
                            
                                    __m128d nrg = _mm_mul_pd(sse_sr, sse_sr);
                                    nrg = _mm_mul_pd(nrg, sse_k_rf[k]);
                                    nrg = _mm_sub_pd(nrg, sse_c_rf[k]);
                                    nrg = _mm_add_pd(nrg, sse_one_over_sr);

                                    __m128d sse_chg = _mm_set_pd( param10.reduced_charge,
                                                                  param11.reduced_charge );
                        
                                    sse_chg = _mm_mul_pd(sse_chg, sse_chg0);
                        
                                    nrg = _mm_mul_pd(sse_chg, nrg);

                                    nrg = _mm_and_pd(nrg, sse_in_cutoff);
                                
                                    sse_cnrg[k] = _mm_add_pd(sse_cnrg[k], nrg);
                                }
                                
                                //lj energy
                                {
                                    const __m128d sse_in_cutoff = _mm_cmplt_pd(sse_r2, sse_Rlj2);

                                    //calculate shift = alpha * sigma * shift_delta
                                    const __m128d sse_shift = _mm_mul_pd(sse_sig, sse_delta[k]);

                                    __m128d lj_denom = _mm_add_pd(sse_r2, sse_shift);
                                    __m128d lj_denom2 = _mm_mul_pd(lj_denom, lj_denom);
                                    lj_denom = _mm_mul_pd(lj_denom, lj_denom2);
                            
                                    const __m128d sig6_over_denom = _mm_div_pd(sse_sig6,
                                                                               lj_denom);
                                                                           
                                    const __m128d sig12_over_denom2 = _mm_mul_pd(sig6_over_denom,
                                                                                 sig6_over_denom);
                                                  
                                    //calculate LJ energy (the factor of 4 is added later)
                                    __m128d tmp = _mm_sub_pd(sig12_over_denom2,
                                                             sig6_over_denom);
                                                         
                                    tmp = _mm_mul_pd(sse_eps, tmp);
                                    tmp = _mm_and_pd(tmp, sse_in_cutoff);
                                    sse_ljnrg[k] = _mm_add_pd(sse_ljnrg[k], tmp);
                                }
                            }
                        }
                              
                        if (remainder == 1)
                        {
                            const Parameter &param1 = params1_array[nats1-1];

                            const double r2 = distmat[nats1-1];
                            
                            const double q2 = param0.reduced_charge * param1.reduced_charge;

                            const LJPair &ljpair = ljpairs.constData()[
                                                    ljpairs.map(param0.ljid,
                                                                param1.ljid)];
                            
                            const double sig2 = ljpair.sigma() * ljpair.sigma();
                            const double sig6 = sig2 * sig2 * sig2;
                            
                            for (int k=0; k<nalpha; ++k)
                            {
                                //coulomb energy
                                {
                                    const double sr = std::sqrt(alfa[k] + r2);

                                    if (sr < sRcoul[k])
                                    {
                                        const double one_over_sr = double(1) / sr;
                                    
                                        icnrg[k] += q2 * (one_over_sr + sr*sr*k_rf[k] - c_rf[k]);
                                    }
                                }
                                
                                //lj energy
                                if (r2 < Rlj2)
                                {
                                    const double shift = ljpair.sigma() * delta[k];
                                
                                    double lj_denom = r2 + shift;
                                    lj_denom = lj_denom * lj_denom * lj_denom;
                                
                                    const double sig6_over_denom = sig6 / lj_denom;
                                    const double sig12_over_denom2 = sig6_over_denom *
                                                                     sig6_over_denom;
                                
                                    iljnrg[k] += ljpair.epsilon() * (sig12_over_denom2 -
                                                                     sig6_over_denom);
                                }
                            }
                        }
                    }
                    
                    for (int k=0; k<nalpha; ++k)
                    {
                        icnrg[k] += *((const double*)&(sse_cnrg[k])) +
                                    *( ((const double*)&(sse_cnrg[k])) + 1 );
                             
                        iljnrg[k] += *((const double*)&(sse_ljnrg[k])) +
                                     *( ((const double*)&(sse_ljnrg[k])) + 1 );
                    }
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

                            const double r2 = distmat[j];
                            
                            const double q2 = param0.reduced_charge * param1.reduced_charge;

                            const LJPair &ljpair = ljpairs.constData()[
                                                    ljpairs.map(param0.ljid,
                                                                param1.ljid)];
                            
                            const double sig2 = ljpair.sigma() * ljpair.sigma();
                            const double sig6 = sig2 * sig2 * sig2;
                            
                            for (int k=0; k<nalpha; ++k)
                            {
                                //coulomb energy
                                {
                                    const double sr = std::sqrt(alfa[k] + r2);
                                    
                                    if (sr < sRcoul[k])
                                    {
                                        const double one_over_sr = double(1) / sr;
                                
                                        icnrg[k] += q2 * (one_over_sr + sr*sr*k_rf[k] - c_rf[k]);
                                    }
                                }
                                
                                //lj energy
                                if (r2 < Rlj2)
                                {
                                    const double shift = ljpair.sigma() * delta[k];
                                
                                    double lj_denom = r2 + shift;
                                    lj_denom = lj_denom * lj_denom * lj_denom;
                                
                                    const double sig6_over_denom = sig6 / lj_denom;
                                    const double sig12_over_denom2 = sig6_over_denom *
                                                                     sig6_over_denom;
            
                                    iljnrg[k] += ljpair.epsilon() * (sig12_over_denom2 - 
                                                                     sig6_over_denom);
                                }
                            }
                        }
                    }
                }
                #endif
                
                //now add these energies onto the total for the molecule,
                //scaled by any non-bonded feather factor
                for (int i=0; i<nalpha; ++i)
                {
                    cnrg[i] += icnrg[i];
                    ljnrg[i] += iljnrg[i];
                }
            }
        }
    }
    else if (use_atomistic_cutoff)
    {
        //use a straight atomistic cutoff
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
                                            spce->beyond(Rc, aabox0, group1.aaBox());
                
                if (not within_cutoff)
                    //this CutGroup is either the cutoff distance
                    continue;
                
                //calculate all of the interatomic distances^2
                const double mindist = spce->calcDist2(group0, group1, distmat);
                
                if (mindist > Rc)
                {
                    //all of the atoms are definitely beyond cutoff
                    continue;
                }
                   
                double icnrg[nalpha];
                double iljnrg[nalpha];
                
                for (int i=0; i<nalpha; ++i)
                {
                    icnrg[i] = 0;
                    iljnrg[i] = 0;
                }
                
                //loop over all interatomic pairs and calculate the energies
                const quint32 nats1 = group1.count();
                const Parameter *params1_array = params1.constData();

                #ifdef SIRE_USE_SSE
                {
                    const int remainder = nats1 % 2;
                    
                    const __m128d sse_one = { 1.0, 1.0 };
                    
                    __m128d sse_cnrg[nalpha];
                    __m128d sse_ljnrg[nalpha];
                    __m128d sse_alpha[nalpha];
                    __m128d sse_delta[nalpha];
                    
                    const __m128d sse_Rcoul = _mm_set1_pd(Rcoul);
                    const __m128d sse_Rlj2 = _mm_set1_pd(Rlj2);
                    
                    for (int i=0; i<nalpha; ++i)
                    {
                        sse_cnrg[i] = _mm_set_pd(0, 0);
                        sse_ljnrg[i] = _mm_set_pd(0, 0);
                        
                        sse_alpha[i] = _mm_set_pd(alfa[i], alfa[i]);
                        sse_delta[i] = _mm_set_pd(delta[i], delta[i]);
                    }
                    
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
                            
                            const __m128d sse_r2 = _mm_set_pd(distmat[j], distmat[j+1]);

                                                       
                            const LJPair &ljpair0 = ljpairs.constData()[
                                                    ljpairs.map(param0.ljid,
                                                                param10.ljid)];
                        
                            const LJPair &ljpair1 = ljpairs.constData()[
                                                    ljpairs.map(param0.ljid,
                                                                param11.ljid)];
                        
                            __m128d sse_sig = _mm_set_pd(ljpair0.sigma(), 
                                                         ljpair1.sigma());
                            __m128d sse_eps = _mm_set_pd(ljpair0.epsilon(), 
                                                         ljpair1.epsilon());

                            const __m128d sse_sig2 = _mm_mul_pd(sse_sig, sse_sig);
                            const __m128d sse_sig3 = _mm_mul_pd(sse_sig2, sse_sig);
                            const __m128d sse_sig6 = _mm_mul_pd(sse_sig3, sse_sig3);
                            
                            for (int k=0; k<nalpha; ++k)
                            {
                                //coulomb energy
                                {
                                    __m128d sse_sr = _mm_sqrt_pd( _mm_add_pd(sse_r2,sse_alpha[k]) );
                                    const __m128d sse_in_cutoff = _mm_cmplt_pd(sse_sr, sse_Rcoul);
                            
                                    __m128d nrg = _mm_div_pd(sse_one, sse_sr);

                                    __m128d sse_chg = _mm_set_pd( param10.reduced_charge,
                                                                  param11.reduced_charge );
                        
                                    sse_chg = _mm_mul_pd(sse_chg, sse_chg0);
                        
                                    nrg = _mm_mul_pd(sse_chg, nrg);

                                    nrg = _mm_and_pd(nrg, sse_in_cutoff);
                                
                                    sse_cnrg[k] = _mm_add_pd(sse_cnrg[k], nrg);
                                }
                                
                                //lj energy
                                {
                                    const __m128d sse_in_cutoff = _mm_cmplt_pd(sse_r2, sse_Rlj2);

                                    //calculate shift = alpha * sigma * shift_delta
                                    const __m128d sse_shift = _mm_mul_pd(sse_sig, sse_delta[k]);

                                    __m128d lj_denom = _mm_add_pd(sse_r2, sse_shift);
                                    __m128d lj_denom2 = _mm_mul_pd(lj_denom, lj_denom);
                                    lj_denom = _mm_mul_pd(lj_denom, lj_denom2);
                            
                                    const __m128d sig6_over_denom = _mm_div_pd(sse_sig6,
                                                                               lj_denom);
                                                                           
                                    const __m128d sig12_over_denom2 = _mm_mul_pd(sig6_over_denom,
                                                                                 sig6_over_denom);
                                                  
                                    //calculate LJ energy (the factor of 4 is added later)
                                    __m128d tmp = _mm_sub_pd(sig12_over_denom2,
                                                             sig6_over_denom);
                                                         
                                    tmp = _mm_mul_pd(sse_eps, tmp);
                                    tmp = _mm_and_pd(tmp, sse_in_cutoff);
                                    sse_ljnrg[k] = _mm_add_pd(sse_ljnrg[k], tmp);
                                }
                            }
                        }
                              
                        if (remainder == 1)
                        {
                            const Parameter &param1 = params1_array[nats1-1];

                            const double r2 = distmat[nats1-1];
                            
                            const double q2 = param0.reduced_charge * param1.reduced_charge;


                            const LJPair &ljpair = ljpairs.constData()[
                                                    ljpairs.map(param0.ljid,
                                                                param1.ljid)];
                            
                            const double sig2 = ljpair.sigma() * ljpair.sigma();
                            const double sig6 = sig2 * sig2 * sig2;
                            
                            for (int k=0; k<nalpha; ++k)
                            {
                                //coulomb energy
                                {
                                    const double sr = std::sqrt(alfa[k] + r2);
                                    
                                    if (sr < Rcoul)
                                    {
                                        const double one_over_sr = double(1) / sr;
                                    
                                        icnrg[k] += q2 * one_over_sr;
                                    }
                                }
                                
                                //lj energy
                                if (r2 < Rlj2)
                                {
                                    const double shift = ljpair.sigma() * delta[k];
                                
                                    double lj_denom = r2 + shift;
                                    lj_denom = lj_denom * lj_denom * lj_denom;
                                
                                    const double sig6_over_denom = sig6 / lj_denom;
                                    const double sig12_over_denom2 = sig6_over_denom * 
                                                                     sig6_over_denom;
                                
                                    iljnrg[k] += ljpair.epsilon() * (sig12_over_denom2 - 
                                                                     sig6_over_denom);
                                }
                            }
                        }
                    }
                    
                    for (int k=0; k<nalpha; ++k)
                    {
                        icnrg[k] += *((const double*)&(sse_cnrg[k])) +
                                    *( ((const double*)&(sse_cnrg[k])) + 1 );
                             
                        iljnrg[k] += *((const double*)&(sse_ljnrg[k])) +
                                     *( ((const double*)&(sse_ljnrg[k])) + 1 );
                    }
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

                            const double r2 = distmat[j];
                            
                            const double q2 = param0.reduced_charge * param1.reduced_charge;

                            const LJPair &ljpair = ljpairs.constData()[
                                                    ljpairs.map(param0.ljid,
                                                                param1.ljid)];
                            
                            const double sig2 = ljpair.sigma() * ljpair.sigma();
                            const double sig6 = sig2 * sig2 * sig2;
                            
                            for (int k=0; k<nalpha; ++k)
                            {
                                //coulomb energy
                                {
                                    const double sr = std::sqrt(alfa[k] + r2);
                                    
                                    if (sr < Rcoul)
                                    {
                                        const double one_over_sr = double(1) / sr;

                                        icnrg[k] += q2 * one_over_sr;
                                    }
                                }
                                
                                //lj energy
                                if (r2 < Rlj2)
                                {
                                    const double shift = ljpair.sigma() * delta[k];
                                
                                    double lj_denom = r2 + shift;
                                    lj_denom = lj_denom * lj_denom * lj_denom;
                                
                                    const double sig6_over_denom = sig6 / lj_denom;
                                    const double sig12_over_denom2 = sig6_over_denom *
                                                                     sig6_over_denom;
            
                                    iljnrg[k] += ljpair.epsilon() * (sig12_over_denom2 - 
                                                                     sig6_over_denom);
                                }
                            }
                        }
                    }
                }
                #endif
                
                //now add these energies onto the total for the molecule,
                //scaled by any non-bonded feather factor
                for (int i=0; i<nalpha; ++i)
                {
                    cnrg[i] += icnrg[i];
                    ljnrg[i] += iljnrg[i];
                }
            }
        }
    }
    else // group-based feathered cutoff
    {
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
                
                //calculate all of the interatomic distances^2
                const double mindist = spce->calcDist2(group0, group1, distmat);
                
                if (mindist > switchfunc->cutoffDistance())
                {
                    //all of the atoms are definitely beyond cutoff
                    continue;
                }
                   
                double icnrg[nalpha];
                double iljnrg[nalpha];
                
                for (int i=0; i<nalpha; ++i)
                {
                    icnrg[i] = 0;
                    iljnrg[i] = 0;
                }
                
                //loop over all interatomic pairs and calculate the energies
                const quint32 nats1 = group1.count();
                const Parameter *params1_array = params1.constData();

                #ifdef SIRE_USE_SSE
                {
                    const int remainder = nats1 % 2;
                    
                    __m128d sse_cnrg[nalpha];
                    __m128d sse_ljnrg[nalpha];
                    __m128d sse_alpha[nalpha];
                    __m128d sse_delta[nalpha];
                                    
                    for (int i=0; i<nalpha; ++i)
                    {
                        sse_cnrg[i] = _mm_set1_pd(0);
                        sse_ljnrg[i] = _mm_set1_pd(0);
                        
                        sse_alpha[i] = _mm_set1_pd(alfa[i]);
                        sse_delta[i] = _mm_set1_pd(delta[i]);
                    }
                    
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
                            
                            __m128d sse_dist2 = _mm_set_pd(distmat[j], distmat[j+1]);
                            __m128d sse_chg1 = _mm_set_pd(param10.reduced_charge,
                                                          param11.reduced_charge);
                                               
                            const LJPair &ljpair0 = ljpairs.constData()[
                                                    ljpairs.map(param0.ljid,
                                                                param10.ljid)];
                        
                            const LJPair &ljpair1 = ljpairs.constData()[
                                                    ljpairs.map(param0.ljid,
                                                                param11.ljid)];
                        
                            __m128d sse_sig = _mm_set_pd(ljpair0.sigma(), 
                                                         ljpair1.sigma());
                            __m128d sse_eps = _mm_set_pd(ljpair0.epsilon(), 
                                                         ljpair1.epsilon());
                            
                            //calculate the coulomb energy
                            __m128d sse_q2 = _mm_mul_pd(sse_chg0, sse_chg1);
                            
                            for (int k=0; k<nalpha; ++k)
                            {
                                __m128d tmp = _mm_add_pd( sse_dist2, sse_alpha[k] );
                                __m128d coul_denom = _mm_sqrt_pd( tmp );
                            
                                coul_denom = _mm_div_pd(sse_q2, coul_denom);
                            
                                sse_cnrg[k] = _mm_add_pd(sse_cnrg[k], coul_denom);
                            }

                            const __m128d sse_sig2 = _mm_mul_pd(sse_sig, sse_sig);
                            const __m128d sse_sig3 = _mm_mul_pd(sse_sig2, sse_sig);
                            const __m128d sse_sig6 = _mm_mul_pd(sse_sig3, sse_sig3);
                            
                            for (int k=0; k<nalpha; ++k)
                            {
                                //calculate shift = alpha * sigma * shift_delta
                                const __m128d sse_shift = _mm_mul_pd(sse_sig, sse_delta[k]);

                                __m128d lj_denom = _mm_add_pd(sse_dist2, sse_shift);
                                __m128d lj_denom2 = _mm_mul_pd(lj_denom, lj_denom);
                                lj_denom = _mm_mul_pd(lj_denom, lj_denom2);
                            
                                const __m128d sig6_over_denom = _mm_div_pd(sse_sig6, 
                                                                           lj_denom);
                                                                           
                                const __m128d sig12_over_denom2 = _mm_mul_pd(sig6_over_denom, 
                                                                             sig6_over_denom);
                                                  
                                //calculate LJ energy (the factor of 4 is added later)
                                __m128d tmp = _mm_sub_pd(sig12_over_denom2,
                                                         sig6_over_denom);
                                                         
                                tmp = _mm_mul_pd(sse_eps, tmp);
                                sse_ljnrg[k] = _mm_add_pd(sse_ljnrg[k], tmp);
                            }
                        }
                              
                        if (remainder == 1)
                        {
                            const Parameter &param1 = params1_array[nats1-1];

                            const double dist2 = distmat[nats1-1];
                            
                            const double q2 = param0.reduced_charge * param1.reduced_charge;
                            
                            for (int k=0; k<nalpha; ++k)
                            {
                                icnrg[k] += q2 / std::sqrt( alfa[k] + dist2 );
                            }

                            const LJPair &ljpair = ljpairs.constData()[
                                                    ljpairs.map(param0.ljid,
                                                                param1.ljid)];
                            
                            const double sig2 = ljpair.sigma() * ljpair.sigma();
                            const double sig6 = sig2 * sig2 * sig2;

                            for (int k=0; k<nalpha; ++k)
                            {
                                const double shift = ljpair.sigma() * delta[k];
                            
                                double lj_denom = dist2 + shift; 
                                lj_denom = lj_denom * lj_denom * lj_denom;
                            
                                const double sig6_over_denom = sig6 / lj_denom;
                                const double sig12_over_denom2 = sig6_over_denom * 
                                                                 sig6_over_denom;
                            
                                iljnrg[k] += ljpair.epsilon() * (sig12_over_denom2 - 
                                                                 sig6_over_denom);
                            }
                        }
                    }
                    
                    for (int k=0; k<nalpha; ++k)
                    {
                        icnrg[k] += *((const double*)&(sse_cnrg[k])) +
                                    *( ((const double*)&(sse_cnrg[k])) + 1 );
                             
                        iljnrg[k] += *((const double*)&(sse_ljnrg[k])) +
                                     *( ((const double*)&(sse_ljnrg[k])) + 1 );
                    }
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

                            const double dist2 = distmat[j];
                            
                            const double q2 = param0.reduced_charge * param1.reduced_charge;
                            
                            for (int k=0; k<nalpha; ++k)
                            {
                                icnrg[k] += q2 / std::sqrt(alfa[k] + dist2);
                            }

                            const LJPair &ljpair = ljpairs.constData()[
                                                    ljpairs.map(param0.ljid,
                                                                param1.ljid)];
                            
                            const double sig2 = ljpair.sigma() * ljpair.sigma();
                            const double sig6 = sig2 * sig2 * sig2;
                            
                            for (int k=0; k<nalpha; ++k)
                            {
                                const double shift = ljpair.sigma() * delta[k];
                            
                                double lj_denom = dist2 + shift;
                                lj_denom = lj_denom * lj_denom * lj_denom;
                            
                                const double sig6_over_denom = sig6 / lj_denom;
                                const double sig12_over_denom2 = sig6_over_denom *
                                                                 sig6_over_denom;
        
                                iljnrg[k] += ljpair.epsilon() * (sig12_over_denom2 - 
                                                                 sig6_over_denom);
                            }
                        }
                    }
                }
                #endif
                
                //now add these energies onto the total for the molecule,
                //scaled by any non-bonded feather factor
                if (mindist > switchfunc->electrostaticFeatherDistance())
                {
                    const double cscl = switchfunc->electrostaticScaleFactor(Length(mindist));
                    
                    for (int i=0; i<nalpha; ++i)
                    {
                        cnrg[i] += cscl * icnrg[i];
                    }
                }
                else
                {
                    for (int i=0; i<nalpha; ++i)
                    {
                        cnrg[i] += icnrg[i];
                    }
                }
                
                if (mindist > switchfunc->vdwFeatherDistance())
                {
                    const double ljscl = switchfunc->vdwScaleFactor(Length(mindist));
                    
                    for (int i=0; i<nalpha; ++i)
                    {
                        ljnrg[i] += ljscl * iljnrg[i];
                    }
                }
                else
                {
                    for (int i=0; i<nalpha; ++i)
                    {
                        ljnrg[i] += iljnrg[i];
                    }
                }
            }
        } // end of if use_electrostatic_shifting
    }
    
    //add this molecule pair's energy onto the total
    //(also multiply LJ by 4 as it is 4 * epsilon ((sig/r)^12 - (sig/r)^6))
    //(also multiply COUL by (1-alpha)^n)
    for (int i=0; i<nalpha; ++i)
    {
        cnrg[i] *= scale_energy * one_minus_alfa_to_n[i];
        ljnrg[i] *= 4 * scale_energy;
    }

    //now copy the calculated energies back to the Energy object
    Energy soft_energy;
    
    for (int i=0; i<alpha_index.count(); ++i)
    {
        int idx = alpha_index.at(i);
        
        if (idx >= 0)
            soft_energy.setEnergy(i, cnrg[idx], ljnrg[idx]);
    }
    
    energy += soft_energy;
}

/** Add to the forces in 'forces0' the forces acting on 'mol0' caused
    by 'mol1' */
void InterSoftCLJPotential::_pvt_calculateForce(
                                         const InterSoftCLJPotential::Molecule &mol0, 
                                         const InterSoftCLJPotential::Molecule &mol1,
                                         MolForceTable &forces0, 
                                         InterSoftCLJPotential::ForceWorkspace &distmat,
                                         double scale_force) const
{
    throw SireError::incomplete_code( QObject::tr(
            "Soft-core forces are not yet implemented!"), CODELOC );

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
                        
                                total_force -= cforce;
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

                                force += ((scl_lj * ljpair.epsilon() * 
                                            (6.0*sig_over_dist6 - 12.0*sig_over_dist12))
                                            * distmat[j].direction())
                                            
                                          + (ljnrg * dscl_lj);
                            }

                            total_force -= force;
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
                        
                                total_force -= cforce;
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
                        
                            total_force -= force;
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
void InterSoftCLJPotential::_pvt_calculateCoulombForce(
                                         const InterSoftCLJPotential::Molecule &mol0, 
                                         const InterSoftCLJPotential::Molecule &mol1,
                                         MolForceTable &forces0, 
                                         InterSoftCLJPotential::ForceWorkspace &distmat,
                                         double scale_force) const
{
    throw SireError::incomplete_code( QObject::tr(
            "Soft-core forces are not yet implemented!"), CODELOC );

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
                        
                            total_force -= cforce;
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
                        
                            total_force -= cforce;
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
void InterSoftCLJPotential::_pvt_calculateLJForce(
                                         const InterSoftCLJPotential::Molecule &mol0, 
                                         const InterSoftCLJPotential::Molecule &mol1,
                                         MolForceTable &forces0, 
                                         InterSoftCLJPotential::ForceWorkspace &distmat,
                                         double scale_force) const
{
    throw SireError::incomplete_code( QObject::tr(
            "Soft-core forces are not yet implemented!"), CODELOC );

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

                                total_force -= force;
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

                                total_force -= force;
                            }
                        }
                    }
                    
                    group_forces0_array[i] += scale_force * total_force;

                } // end of loop over i atoms

            } // end of if within feather

        } // end of loop over jgroup CutGroups

    } // end of loop over igroup CutGroups
}

void InterSoftCLJPotential::_pvt_calculateField(
                         const InterSoftCLJPotential::Molecule &mol0, 
                         const InterSoftCLJPotential::Molecule &mol1,
                         const CLJProbe &probe,
                         MolFieldTable &fields0, 
                         InterSoftCLJPotential::FieldWorkspace &workspace,
                         double scale_field) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code to calculate soft coulomb and LJ fields has not "
                "yet been written..."), CODELOC );
}

void InterSoftCLJPotential::_pvt_calculateField(const InterSoftCLJPotential::Molecule &mol,
                         const CLJProbe &probe,
                         GridFieldTable &fields,
                         InterSoftCLJPotential::FieldWorkspace &workspace,
                         double scale_field) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code to calculate soft coulomb and LJ fields has not "
                "yet been written..."), CODELOC );
}

void InterSoftCLJPotential::_pvt_calculateCoulombField(const InterSoftCLJPotential::Molecule &mol0, 
                                const InterSoftCLJPotential::Molecule &mol1,
                                const CLJProbe &probe,
                                MolFieldTable &fields0, 
                                InterSoftCLJPotential::FieldWorkspace &workspace,
                                double scale_field) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code to calculate soft coulomb and LJ fields has not "
                "yet been written..."), CODELOC );
}

void InterSoftCLJPotential::_pvt_calculateCoulombField(const InterSoftCLJPotential::Molecule &mol,
                                const CLJProbe &probe,
                                GridFieldTable &fields,
                                InterSoftCLJPotential::FieldWorkspace &workspace,
                                double scale_field) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code to calculate soft coulomb and LJ fields has not "
                "yet been written..."), CODELOC );
}

void InterSoftCLJPotential::_pvt_calculateLJField(const InterSoftCLJPotential::Molecule &mol0, 
                           const InterSoftCLJPotential::Molecule &mol1,
                           const CLJProbe &probe,
                           MolFieldTable &fields0, 
                           InterSoftCLJPotential::FieldWorkspace &workspace,
                           double scale_field) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code to calculate soft coulomb and LJ fields has not "
                "yet been written..."), CODELOC );
}

void InterSoftCLJPotential::_pvt_calculateLJField(const InterSoftCLJPotential::Molecule &mol,
                           const CLJProbe &probe,
                           GridFieldTable &fields,
                           InterSoftCLJPotential::FieldWorkspace &workspace,
                           double scale_field) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code to calculate soft coulomb and LJ fields has not "
                "yet been written..."), CODELOC );
}

void InterSoftCLJPotential::_pvt_calculatePotential(const InterSoftCLJPotential::Molecule &mol0, 
                             const InterSoftCLJPotential::Molecule &mol1,
                             const CLJProbe &probe,
                             MolPotentialTable &pots0, 
                             InterSoftCLJPotential::PotentialWorkspace &workspace,
                             double scale_potential) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code to calculate soft coulomb and LJ fields has not "
                "yet been written..."), CODELOC );
}

void InterSoftCLJPotential::_pvt_calculatePotential(const InterSoftCLJPotential::Molecule &mol,
                             const CLJProbe &probe,
                             GridPotentialTable &fields,
                             InterSoftCLJPotential::PotentialWorkspace &workspace,
                             double scale_potential) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code to calculate soft coulomb and LJ fields has not "
                "yet been written..."), CODELOC );
}

void InterSoftCLJPotential::_pvt_calculateCoulombPotential(const InterSoftCLJPotential::Molecule &mol0, 
                                    const InterSoftCLJPotential::Molecule &mol1,
                                    const CLJProbe &probe,
                                    MolPotentialTable &pots0, 
                                    InterSoftCLJPotential::PotentialWorkspace &workspace,
                                    double scale_potential) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code to calculate soft coulomb and LJ fields has not "
                "yet been written..."), CODELOC );
}

void InterSoftCLJPotential::_pvt_calculateCoulombPotential(const InterSoftCLJPotential::Molecule &mol,
                                    const CLJProbe &probe,
                                    GridPotentialTable &fields,
                                    InterSoftCLJPotential::PotentialWorkspace &workspace,
                                    double scale_potential) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code to calculate soft coulomb and LJ fields has not "
                "yet been written..."), CODELOC );
}

void InterSoftCLJPotential::_pvt_calculateLJPotential(const InterSoftCLJPotential::Molecule &mol0, 
                               const InterSoftCLJPotential::Molecule &mol1,
                               const CLJProbe &probe,
                               MolPotentialTable &pots0, 
                               InterSoftCLJPotential::PotentialWorkspace &workspace,
                               double scale_potential) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code to calculate soft coulomb and LJ fields has not "
                "yet been written..."), CODELOC );
}

void InterSoftCLJPotential::_pvt_calculateLJPotential(const InterSoftCLJPotential::Molecule &mol,
                               const CLJProbe &probe,
                               GridPotentialTable &fields,
                               InterSoftCLJPotential::PotentialWorkspace &workspace,
                               double scale_potential) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code to calculate soft coulomb and LJ fields has not "
                "yet been written..."), CODELOC );
}

//////////////
////////////// Implementation of IntraSoftCLJPotential
/////////////

static const RegisterMetaType<IntraSoftCLJPotential> r_intraclj( MAGIC_ONLY, NO_ROOT,
                                                   IntraSoftCLJPotential::typeName() );

/** Serialise to a binary datastream */
QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds,
                                      const IntraSoftCLJPotential &intraclj)
{
    writeHeader(ds, r_intraclj, 1);
    
    ds << static_cast<const SoftCLJPotential&>(intraclj);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds,
                                      IntraSoftCLJPotential &intraclj)
{
    VersionID v = readHeader(ds, r_intraclj);
    
    if (v == 1)
    {
        ds >> static_cast<SoftCLJPotential&>(intraclj);
    }
    else
        throw version_error(v, "1", r_intraclj, CODELOC);
        
    return ds;
}

/** Constructor */
IntraSoftCLJPotential::IntraSoftCLJPotential() : SoftCLJPotential()
{}

/** Copy constructor */
IntraSoftCLJPotential::IntraSoftCLJPotential(const IntraSoftCLJPotential &other)
                  : SoftCLJPotential(other)
{}

/** Destructor */
IntraSoftCLJPotential::~IntraSoftCLJPotential()
{}

/** Copy assignment operator */
IntraSoftCLJPotential& 
IntraSoftCLJPotential::operator=(const IntraSoftCLJPotential &other)
{
    SoftCLJPotential::operator=(other);
    return *this;
}

/** Return all of the parameters needed by this potential for 
    the molecule 'molecule', using the supplied property map to
    find the properties that contain those parameters
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
IntraSoftCLJPotential::Parameters 
IntraSoftCLJPotential::getParameters(const PartialMolecule &molecule,
				     const PropertyMap &map)
{
    need_update_ljpairs = true;

    return Parameters( AtomicParameters3D<CLJParameter>(
                               molecule, map[parameters().coordinates()],
                               CLJPotential::getCLJParameters(molecule, 
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
IntraSoftCLJPotential::Parameters
IntraSoftCLJPotential::updateParameters(const IntraSoftCLJPotential::Parameters &old_params,
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
        new_params.setAtomicParameters( CLJPotential::getCLJParameters(new_molecule,
                                            chg_property, lj_property) );

	//need_update_ljpairs = true;
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
IntraSoftCLJPotential::Parameters
IntraSoftCLJPotential::updateParameters(const IntraSoftCLJPotential::Parameters &old_params,
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
      new_params.setAtomicParameters( CLJPotential::getCLJParameters(new_molecule,
                                                         new_chg, new_lj) );
                                                         
      //       need_update_ljpairs = true;
    }

    if (changed_scl)
        new_params.setIntraScaleFactors( 
                        IntraScaledParameters<CLJNBPairs>(new_molecule, new_scl) );

    return new_params;
}

/** Return the IntraSoftCLJPotential::Molecule representation of 'molecule',
    using the supplied PropertyMap to find the properties that contain
    the necessary forcefield parameters
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
IntraSoftCLJPotential::Molecule
IntraSoftCLJPotential::parameterise(const PartialMolecule &molecule,
                                    const PropertyMap &map)
{
    return IntraSoftCLJPotential::Molecule(molecule, *this, map);
}

/** Convert the passed group of molecules into IntraSoftCLJPotential::Molecules,
    using the supplied PropertyMap to find the properties that contain
    the necessary forcefield parameters in each molecule
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
IntraSoftCLJPotential::Molecules 
IntraSoftCLJPotential::parameterise(const MoleculeGroup &molecules,
                                    const PropertyMap &map)
{
    return IntraSoftCLJPotential::Molecules(molecules, *this, map);
}


/** Assert that 'rest_of_mol' is compatible with 'mol'. They are only 
    compatible if they are both part of the same molecule (not necessarily
    the same version) with the same layout UID.
    
    \throw SireError::incompatible_error
*/
void IntraSoftCLJPotential::assertCompatible(const IntraSoftCLJPotential::Molecule &mol,
					     const IntraSoftCLJPotential::Molecule &rest_of_mol) const
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

/** Return the total charge of the parameters for the group in 'params' */
double IntraSoftCLJPotential::totalCharge(
                            const IntraSoftCLJPotential::Parameters::Array &params) const
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

void IntraSoftCLJPotential::_pvt_calculateEnergy(const CLJNBPairs::CGPairs &group_pairs, 
						 IntraSoftCLJPotential::EnergyWorkspace &distmat,
						 const IntraSoftCLJPotential::Parameter *params0_array, 
						 const IntraSoftCLJPotential::Parameter *params1_array,
						 const quint32 nats0, const quint32 nats1, 
						 double icnrg[], double iljnrg[],
						 const double alfa[], double delta[], const int nalpha) const
{
    const double Rcoul = qMax(1e-5,qMin(1e9,
                            switchfunc->electrostaticCutoffDistance().to(angstrom)));
    const double Rlj = qMax(1e-5,qMin(1e9, switchfunc->vdwCutoffDistance().to(angstrom)) );
    const double Rlj2 = Rlj*Rlj;
    const double Rcoul2 = Rcoul*Rcoul;

    if (group_pairs.isEmpty())
    {
        //there is a constant scale factor between groups
        CLJScaleFactor cljscl = group_pairs.defaultValue();

        if (cljscl.coulomb() == 0 and cljscl.lj() == 0)
            return;

        if (use_electrostatic_shifting)
        {
            double sRcoul[nalpha];
            double one_over_sRcoul[nalpha];
            double one_over_sRcoul2[nalpha];
        
            for (int i=0; i<nalpha; ++i)
            {
                sRcoul[i] = std::sqrt(alfa[i] + Rcoul*Rcoul);
                one_over_sRcoul[i] = double(1) / sRcoul[i];
                one_over_sRcoul2[i] = double(1) / (sRcoul[i]*sRcoul[i]);
            }

            for (quint32 i=0; i<nats0; ++i)
            {
                distmat.setOuterIndex(i);
                const Parameter &param0 = params0_array[i];
	    
                if (param0.ljid == 0)
                {
                    //null LJ parameter - only add on the coulomb energy
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        const double dist2 = distmat[j];
                        const double q2 = cljscl.coulomb() *
                            param0.reduced_charge * params1_array[j].reduced_charge;

                        for (int k=0; k<nalpha; ++k)
                        {
                            const double sr = std::sqrt(alfa[k] + dist2);
                            const double one_over_sr = double(1) / sr;
 
                            // JM Jan 13. No reaction field on 1,4 atoms?
                            if (sr < sRcoul[k])
                            {
                                if (cljscl.coulomb() != 1)
                                    icnrg[k] += q2 * (one_over_sr);
                                else
                                    icnrg[k] += q2 * (one_over_sr - one_over_sRcoul[k] +
                                                      one_over_sRcoul2[k]*(sr-sRcoul[k]));
                            }
                        }
                    }
                }
                else  // do both
                {
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        //do both coulomb and LJ
                        const Parameter &param1 = params1_array[j];

                        const double dist2 = distmat[j];
                        const double q2 = cljscl.coulomb() *
                                param0.reduced_charge * param1.reduced_charge;
                        
                        for (int k=0; k<nalpha; ++k)
                        {
                            const double sr = std::sqrt(alfa[k] + dist2);
                            const double one_over_sr = double(1) / sr;
 
                            // JM Jan 13. No reaction field on 1,4 atoms?
                            if (sr < sRcoul[k])
                            {
                                if (cljscl.coulomb() != 1)
                                    icnrg[k] += q2 * (one_over_sr);
                                else
                                    icnrg[k] += q2 * (one_over_sr - one_over_sRcoul[k] +
                                                      one_over_sRcoul2[k]*(sr-sRcoul[k]));
                            }
                        }
		    
                        if (param1.ljid != 0 and dist2 < Rlj2)
                        {
                            const LJPair &ljpair = ljpairs.constData()[
                                                        ljpairs.map(param0.ljid,
                                                                    param1.ljid)];
                            const double sig2 = ljpair.sigma() * ljpair.sigma();
                            const double sig6 = sig2 * sig2 * sig2;
                            
                            for (int k=0; k<nalpha; ++k)
                            {
                                const double shift = ljpair.sigma() * delta[k];
			  
                                double lj_denom = dist2 + shift;
                                lj_denom = lj_denom * lj_denom * lj_denom;
			  
                                const double sig6_over_denom = sig6 / lj_denom;
                                const double sig12_over_denom2 = sig6_over_denom *
                                                                 sig6_over_denom;
			  
                                iljnrg[k] += cljscl.lj() * ljpair.epsilon() * (sig12_over_denom2 -
                                                sig6_over_denom);
                            }
                        }
                    }// quint j
                }//do both
            }//quint i
        }
        else if (use_reaction_field)
        {
            double sRcoul[nalpha];
            double k_rf[nalpha];
            double c_rf[nalpha];
        
            for (int i=0; i<nalpha; ++i)
            {
                sRcoul[i] = std::sqrt(alfa[i] + Rcoul*Rcoul);
                k_rf[i] = (1.0 / pow_3(sRcoul[i])) * ( (rf_dielectric_constant-1) /
                                                       (2*rf_dielectric_constant + 1) );
                c_rf[i] = (1.0 / sRcoul[i]) * ( (3*rf_dielectric_constant) /
                                                (2*rf_dielectric_constant + 1) );
            }

            for (quint32 i=0; i<nats0; ++i)
            {
                distmat.setOuterIndex(i);
                const Parameter &param0 = params0_array[i];
	    
                if (param0.ljid == 0)
                {
                    //null LJ parameter - only add on the coulomb energy
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        const double dist2 = distmat[j];
                        const double q2 = cljscl.coulomb() *
                            param0.reduced_charge * params1_array[j].reduced_charge;

                        for (int k=0; k<nalpha; ++k)
                        {
                            const double sr = std::sqrt(alfa[k] + dist2);
                            const double one_over_sr = double(1) / sr;
 
                            // JM Jan 13. No reaction field on 1,4 atoms?
                            if (sr < sRcoul[k])
                            {
                                if (cljscl.coulomb() != 1)
                                    icnrg[k] += q2 * (one_over_sr);
                                else
                                    icnrg[k] += q2 * (one_over_sr + sr*sr*k_rf[k] - c_rf[k]);
                            }
                        }
                    }
                }
                else  // do both
                {
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        //do both coulomb and LJ
                        const Parameter &param1 = params1_array[j];

                        const double dist2 = distmat[j];
                        const double q2 = cljscl.coulomb() *
                                param0.reduced_charge * param1.reduced_charge;
                        
                        for (int k=0; k<nalpha; ++k)
                        {
                            const double sr = std::sqrt(alfa[k] + dist2);
                            const double one_over_sr = double(1) / sr;
 
                            // JM Jan 13. No reaction field on 1,4 atoms?
                            if (sr < sRcoul[k])
                            {
                                if (cljscl.coulomb() != 1)
                                    icnrg[k] += q2 * (one_over_sr);
                                else
                                    icnrg[k] += q2 * (one_over_sr + sr*sr*k_rf[k] - c_rf[k]);
                            }
                        }
		    
                        if (param1.ljid != 0 and dist2 < Rlj2)
                        {
                            const LJPair &ljpair = ljpairs.constData()[
                                                        ljpairs.map(param0.ljid,
                                                                    param1.ljid)];
                            const double sig2 = ljpair.sigma() * ljpair.sigma();
                            const double sig6 = sig2 * sig2 * sig2;
                            
                            for (int k=0; k<nalpha; ++k)
                            {
                                const double shift = ljpair.sigma() * delta[k];
			  
                                double lj_denom = dist2 + shift;
                                lj_denom = lj_denom * lj_denom * lj_denom;
			  
                                const double sig6_over_denom = sig6 / lj_denom;
                                const double sig12_over_denom2 = sig6_over_denom *
                                                                 sig6_over_denom;
			  
                                iljnrg[k] += cljscl.lj() * ljpair.epsilon() * (sig12_over_denom2 -
                                                sig6_over_denom);
                            }
                        }
                    }// quint j
                }//do both
            }//quint i
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
                        const double dist2 = distmat[j];
                        const double q2 = cljscl.coulomb() *
                            param0.reduced_charge * params1_array[j].reduced_charge;

                        for (int k=0; k<nalpha; ++k)
                        {
                            const double sr = std::sqrt(alfa[k] + dist2);
                            const double one_over_sr = double(1) / sr;
 
                            // JM Jan 13. No reaction field on 1,4 atoms?
                            if (dist2 < Rcoul2)
                            {
                                icnrg[k] += q2 * (one_over_sr);
                            }
                        }
                    }
                }
                else  // do both
                {
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        //do both coulomb and LJ
                        const Parameter &param1 = params1_array[j];

                        const double dist2 = distmat[j];
                        const double q2 = cljscl.coulomb() *
                                param0.reduced_charge * param1.reduced_charge;
                        
                        for (int k=0; k<nalpha; ++k)
                        {
                            const double sr = std::sqrt(alfa[k] + dist2);
                            const double one_over_sr = double(1) / sr;
 
                            // JM Jan 13. No reaction field on 1,4 atoms?
                            if (dist2 < Rcoul2)
                            {
                                icnrg[k] += q2 * (one_over_sr);
                            }
                        }
		    
                        if (param1.ljid != 0 and dist2 < Rlj2)
                        {
                            const LJPair &ljpair = ljpairs.constData()[
                                                        ljpairs.map(param0.ljid,
                                                                    param1.ljid)];
                            const double sig2 = ljpair.sigma() * ljpair.sigma();
                            const double sig6 = sig2 * sig2 * sig2;
                            
                            for (int k=0; k<nalpha; ++k)
                            {
                                const double shift = ljpair.sigma() * delta[k];
			  
                                double lj_denom = dist2 + shift;
                                lj_denom = lj_denom * lj_denom * lj_denom;
			  
                                const double sig6_over_denom = sig6 / lj_denom;
                                const double sig12_over_denom2 = sig6_over_denom *
                                                                 sig6_over_denom;
			  
                                iljnrg[k] += cljscl.lj() * ljpair.epsilon() * (sig12_over_denom2 -
                                                sig6_over_denom);
                            }
                        }
                    }// quint j
                }//do both
            }//quint i
        }
        else
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
                        const double dist2 = distmat[j];
                        const double q2 = cljscl.coulomb() *
                            param0.reduced_charge * params1_array[j].reduced_charge;
		  
                        for (int k=0; k<nalpha; ++k)
                        {
                            icnrg[k] += q2 / std::sqrt(alfa[k] + dist2);
                        }
                    }
                }
                else
                {
                    for (quint32 j=0; j<nats1; ++j)
                    {
                        //do both coulomb and LJ
                        const Parameter &param1 = params1_array[j];

                        const double dist2 = distmat[j];
                        const double q2 = cljscl.coulomb() *
                                param0.reduced_charge * param1.reduced_charge;
                        
                        for (int k=0; k<nalpha; ++k)
                        {
                            icnrg[k] += q2 / std::sqrt(alfa[k] + dist2);
                        }
		    
                        if (param1.ljid != 0)
                        {
                            const LJPair &ljpair = ljpairs.constData()[
                                                    ljpairs.map(param0.ljid,
                                                                param1.ljid)];

                            const double sig2 = ljpair.sigma() * ljpair.sigma();
                            const double sig6 = sig2 * sig2 * sig2;

                            for (int k=0; k<nalpha; ++k)
                            {
                                const double shift = ljpair.sigma() * delta[k];
			  
                                double lj_denom = dist2 + shift;
                                lj_denom = lj_denom * lj_denom * lj_denom;
			  
                                const double sig6_over_denom = sig6 / lj_denom;
                                const double sig12_over_denom2 = sig6_over_denom *
                                                                 sig6_over_denom;
			  
                                iljnrg[k] += cljscl.lj() * ljpair.epsilon() * (sig12_over_denom2 -
                                                                               sig6_over_denom);
                            }
                        }
                    }
                }
            }
        }
	}
    else
    {
        //there are different nb scale factors between
        //the atoms. We need to calculate the energies using
        //them...
        if (use_electrostatic_shifting)
        {
            double sRcoul[nalpha];
            double one_over_sRcoul[nalpha];
            double one_over_sRcoul2[nalpha];
        
            for (int i=0; i<nalpha; ++i)
            {
                sRcoul[i] = std::sqrt(alfa[i] + Rcoul*Rcoul);
                one_over_sRcoul[i] = double(1) / sRcoul[i];
                one_over_sRcoul2[i] = double(1) / (sRcoul[i]*sRcoul[i]);
            }
	
            for (quint32 i=0; i<nats0; ++i)
            {
                distmat.setOuterIndex(i);
                const Parameter &param0 = params0_array[i];
                
                for (quint32 j=0; j<nats1; ++j)
                {
                    const CLJScaleFactor &cljscl = group_pairs(i,j);

                    if (cljscl.coulomb() != 0 or cljscl.lj() != 0)
                    {
                        const Parameter &param1 = params1_array[j];
		    
                        const double dist2 = distmat[j];
  		        
                        const double q2 = cljscl.coulomb() *
                            param0.reduced_charge * params1_array[j].reduced_charge;

                        // ADAPT FOR RF HERE
                        for (int k=0; k<nalpha; ++k)
                        {
                            const double sr = std::sqrt(alfa[k] + dist2);
                            const double one_over_sr = double(1) / sr;
 
                            // JM Jan 13. No reaction field on 1,4 atoms?
                            if (sr < sRcoul[k])
                            {
                                if (cljscl.coulomb() != 1)
                                    icnrg[k] += q2 * (one_over_sr);
                                else
                                    icnrg[k] += q2 * (one_over_sr - one_over_sRcoul[k] +
                                                      one_over_sRcoul2[k]*(sr-sRcoul[k]));
                            }
                        }
		    
                        if (cljscl.lj() != 0 and param1.ljid != 0 and dist2 < Rlj2)
                        {
                            const LJPair &ljpair = ljpairs.constData()[
                                                        ljpairs.map(param0.ljid,
                                                                    param1.ljid)];
                            
                            const double sig2 = ljpair.sigma() * ljpair.sigma();
                            const double sig6 = sig2 * sig2 * sig2;
			  
                            for (int k=0; k<nalpha; ++k)
                            {
                                const double shift = ljpair.sigma() * delta[k];
			        
                                double lj_denom = dist2 + shift;
                                lj_denom = lj_denom * lj_denom * lj_denom;
			        
                                const double sig6_over_denom = sig6 / lj_denom;
                                const double sig12_over_denom2 = sig6_over_denom *
                                                                 sig6_over_denom;
			        
                                iljnrg[k] += cljscl.lj() * ljpair.epsilon() * (sig12_over_denom2 -
                                                                               sig6_over_denom);
                            }
                        }
                    }
                }//for (quint32 j
            }// for quint 32 i
        }
        else if (use_reaction_field)
        {
            double sRcoul[nalpha];
            double k_rf[nalpha];
            double c_rf[nalpha];
        
            for (int i=0; i<nalpha; ++i)
            {
                sRcoul[i] = std::sqrt(alfa[i] + Rcoul*Rcoul);
                k_rf[i] = (1.0 / pow_3(sRcoul[i])) * ( (rf_dielectric_constant-1) /
                                                       (2*rf_dielectric_constant + 1) );
                c_rf[i] = (1.0 / sRcoul[i]) * ( (3*rf_dielectric_constant) /
                                                (2*rf_dielectric_constant + 1) );
            }
	
            for (quint32 i=0; i<nats0; ++i)
            {
                distmat.setOuterIndex(i);
                const Parameter &param0 = params0_array[i];
                
                for (quint32 j=0; j<nats1; ++j)
                {
                    const CLJScaleFactor &cljscl = group_pairs(i,j);

                    if (cljscl.coulomb() != 0 or cljscl.lj() != 0)
                    {
                        const Parameter &param1 = params1_array[j];
		    
                        const double dist2 = distmat[j];
  		        
                        const double q2 = cljscl.coulomb() *
                            param0.reduced_charge * params1_array[j].reduced_charge;

                        // ADAPT FOR RF HERE
                        for (int k=0; k<nalpha; ++k)
                        {
                            const double sr = std::sqrt(alfa[k] + dist2);
                            const double one_over_sr = double(1) / sr;
 
                            // JM Jan 13. No reaction field on 1,4 atoms?
                            if (sr < sRcoul[k])
                            {
                                if (cljscl.coulomb() != 1)
                                    icnrg[k] += q2 * (one_over_sr);
                                else
                                    icnrg[k] += q2 * (one_over_sr + sr*sr*k_rf[k] - c_rf[k]);
                            }
                        }
		    
                        if (cljscl.lj() != 0 and param1.ljid != 0 and dist2 < Rlj2)
                        {
                            const LJPair &ljpair = ljpairs.constData()[
                                                        ljpairs.map(param0.ljid,
                                                                    param1.ljid)];
                            
                            const double sig2 = ljpair.sigma() * ljpair.sigma();
                            const double sig6 = sig2 * sig2 * sig2;
			  
                            for (int k=0; k<nalpha; ++k)
                            {
                                const double shift = ljpair.sigma() * delta[k];
			        
                                double lj_denom = dist2 + shift;
                                lj_denom = lj_denom * lj_denom * lj_denom;
			        
                                const double sig6_over_denom = sig6 / lj_denom;
                                const double sig12_over_denom2 = sig6_over_denom *
                                                                 sig6_over_denom;
			        
                                iljnrg[k] += cljscl.lj() * ljpair.epsilon() * (sig12_over_denom2 -
                                                                               sig6_over_denom);
                            }
                        }
                    }
                }//for (quint32 j
            }// for quint 32 i
        }
        else if (use_atomistic_cutoff)
        {
            for (quint32 i=0; i<nats0; ++i)
            {
                distmat.setOuterIndex(i);
                const Parameter &param0 = params0_array[i];
                
                for (quint32 j=0; j<nats1; ++j)
                {
                    const CLJScaleFactor &cljscl = group_pairs(i,j);

                    if (cljscl.coulomb() != 0 or cljscl.lj() != 0)
                    {
                        const Parameter &param1 = params1_array[j];
		    
                        const double dist2 = distmat[j];
  		        
                        const double q2 = cljscl.coulomb() *
                            param0.reduced_charge * params1_array[j].reduced_charge;

                        // ADAPT FOR RF HERE
                        for (int k=0; k<nalpha; ++k)
                        {
                            const double sr = std::sqrt(alfa[k] + dist2);
                            const double one_over_sr = double(1) / sr;
 
                            // JM Jan 13. No reaction field on 1,4 atoms?
                            if (dist2 < Rcoul2)
                            {
                                icnrg[k] += q2 * (one_over_sr);
                            }
                        }
		    
                        if (cljscl.lj() != 0 and param1.ljid != 0 and dist2 < Rlj2)
                        {
                            const LJPair &ljpair = ljpairs.constData()[
                                                        ljpairs.map(param0.ljid,
                                                                    param1.ljid)];
                            
                            const double sig2 = ljpair.sigma() * ljpair.sigma();
                            const double sig6 = sig2 * sig2 * sig2;
			  
                            for (int k=0; k<nalpha; ++k)
                            {
                                const double shift = ljpair.sigma() * delta[k];
			        
                                double lj_denom = dist2 + shift;
                                lj_denom = lj_denom * lj_denom * lj_denom;
			        
                                const double sig6_over_denom = sig6 / lj_denom;
                                const double sig12_over_denom2 = sig6_over_denom *
                                                                 sig6_over_denom;
			        
                                iljnrg[k] += cljscl.lj() * ljpair.epsilon() * (sig12_over_denom2 -
                                                                               sig6_over_denom);
                            }
                        }
                    }
                }//for (quint32 j
            }// for quint 32 i
        }
        else //JM normal group based cutoff
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
                            const double dist2 = distmat[j];
		      
                            const double q2 = cljscl.coulomb() *
                                    param0.reduced_charge * params1_array[j].reduced_charge;
		      
                            for (int k=0; k<nalpha; ++k)
                            {
                                icnrg[k] += q2 / std::sqrt(alfa[k] + dist2);
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
                            const double dist2 = distmat[j];
  		        
                            const double q2 = cljscl.coulomb() *
                                param0.reduced_charge * params1_array[j].reduced_charge;
			
                            for (int k=0; k<nalpha; ++k)
                            {
                                icnrg[k] += q2 / std::sqrt(alfa[k] + dist2);
                            }

                            if (cljscl.lj() != 0 and param1.ljid != 0)
                            {
                                const LJPair &ljpair = ljpairs.constData()[
                                                        ljpairs.map(param0.ljid,
                                                                    param1.ljid)];

                                const double sig2 = ljpair.sigma() * ljpair.sigma();
                                const double sig6 = sig2 * sig2 * sig2;
              
                                for (int k=0; k<nalpha; ++k)
                                {
                                    const double shift = ljpair.sigma() * delta[k];
			        
                                    double lj_denom = dist2 + shift;
                                    lj_denom = lj_denom * lj_denom * lj_denom;
			        
                                    const double sig6_over_denom = sig6 / lj_denom;
                                    const double sig12_over_denom2 = sig6_over_denom *
                                                                     sig6_over_denom;
			        
                                    iljnrg[k] += cljscl.lj() * ljpair.epsilon() *
                                            (sig12_over_denom2 - sig6_over_denom);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void IntraSoftCLJPotential::_pvt_calculateEnergy(const CLJNBPairs::CGPairs &group_pairs, 
						 const QSet<Index> &atoms0, const QSet<Index> &atoms1,
						 IntraSoftCLJPotential::EnergyWorkspace &distmat,
						 const IntraSoftCLJPotential::Parameter *params0_array, 
						 const IntraSoftCLJPotential::Parameter *params1_array,
						 const quint32 nats0, const quint32 nats1, 
						 double icnrg[], double iljnrg[],
						 const double alfa[], double delta[], const int nalpha) const
{
    if (atoms0.isEmpty() or atoms1.isEmpty())
        return;

    const double Rcoul = qMax(1e-5,qMin(1e9,
                            switchfunc->electrostaticCutoffDistance().to(angstrom)));
    const double Rlj = qMax(1e-5,qMin(1e9, switchfunc->vdwCutoffDistance().to(angstrom)) );
    const double Rlj2 = Rlj*Rlj;
    const double Rcoul2 = Rcoul*Rcoul;

    if (group_pairs.isEmpty())
    {
        //there is a constant scale factor between groups
        CLJScaleFactor cljscl = group_pairs.defaultValue();

        if (cljscl.coulomb() == 0 and cljscl.lj() == 0)
            return;
        
        if (use_electrostatic_shifting)
        {
            double sRcoul[nalpha];
            double one_over_sRcoul[nalpha];
            double one_over_sRcoul2[nalpha];
        
            for (int i=0; i<nalpha; ++i)
            {
                sRcoul[i] = std::sqrt(alfa[i] + Rcoul*Rcoul);
                one_over_sRcoul[i] = double(1) / sRcoul[i];
                one_over_sRcoul2[i] = double(1) / (sRcoul[i]*sRcoul[i]);
            }
	  
            foreach (Index i, atoms0)
            {
                distmat.setOuterIndex(i);
                const Parameter &param0 = params0_array[i];
	    
                foreach (Index j, atoms1)
                {
                    //do both coulomb and LJ
                    const Parameter &param1 = params1_array[j];
                    const double dist2 = distmat[j];
                        
                    const double q2 = cljscl.coulomb() *
                    param0.reduced_charge * param1.reduced_charge;
                        
                    for (int k=0; k<nalpha; ++k)
                    {
                        //coulomb calculation
                        {
                            const double sr = std::sqrt(alfa[k] + dist2);
                            const double one_over_sr = double(1) / sr;
 
                            // JM Jan 13. No reaction field on 1,4 atoms?
                            if (sr < sRcoul[k])
                            {
                                if (cljscl.coulomb() != 1)
                                    icnrg[k] += q2 * (one_over_sr);
                                else
                                    icnrg[k] += q2 * (one_over_sr - one_over_sRcoul[k] +
                                                      one_over_sRcoul2[k]*(sr-sRcoul[k]));
                            }
                        }
		  
                        //lj calculation
                        if (param1.ljid != 0 and dist2 < Rlj2)
                        {
                            const LJPair &ljpair = ljpairs.constData()[
                                        ljpairs.map(param0.ljid,
                                                    param1.ljid)];
                            const double sig2 = ljpair.sigma() * ljpair.sigma();
                            const double sig6 = sig2 * sig2 * sig2;
			                        
                            for (int k=0; k<nalpha; ++k)
                            {
                                const double shift = ljpair.sigma() * delta[k];
			    
                                double lj_denom = dist2 + shift;
                                lj_denom = lj_denom * lj_denom * lj_denom;
			    
                                const double sig6_over_denom = sig6 / lj_denom;
                                const double sig12_over_denom2 = sig6_over_denom *
                                                                 sig6_over_denom;
			    
                                iljnrg[k] +=  cljscl.lj() * ljpair.epsilon() * (sig12_over_denom2 -
                                                                                sig6_over_denom);
                            }
                        }
                    }
                }// for each j
            }//for each i
        }
        else if (use_reaction_field)
        {
            double sRcoul[nalpha];
            double k_rf[nalpha];
            double c_rf[nalpha];
        
            for (int i=0; i<nalpha; ++i)
            {
                sRcoul[i] = std::sqrt(alfa[i] + Rcoul*Rcoul);
                k_rf[i] = (1.0 / pow_3(sRcoul[i])) * ( (rf_dielectric_constant-1) /
                            (2*rf_dielectric_constant + 1) );
                c_rf[i] = (1.0 / sRcoul[i]) * ( (3*rf_dielectric_constant) /
                            (2*rf_dielectric_constant + 1) );
            }
	  
            foreach (Index i, atoms0)
            {
                distmat.setOuterIndex(i);
                const Parameter &param0 = params0_array[i];
	    
                foreach (Index j, atoms1)
                {
                    //do both coulomb and LJ
                    const Parameter &param1 = params1_array[j];
                    const double dist2 = distmat[j];
                        
                    const double q2 = cljscl.coulomb() *
                    param0.reduced_charge * param1.reduced_charge;
                        
                    for (int k=0; k<nalpha; ++k)
                    {
                        //coulomb calculation
                        {
                            const double sr = std::sqrt(alfa[k] + dist2);
                            const double one_over_sr = double(1) / sr;
 
                            // JM Jan 13. No reaction field on 1,4 atoms?
                            if (sr < sRcoul[k])
                            {
                                if (cljscl.coulomb() != 1)
                                    icnrg[k] += q2 * (one_over_sr);
                                else
                                    icnrg[k] += q2 * (one_over_sr + sr*sr*k_rf[k] - c_rf[k]);
                            }
                        }
		  
                        //lj calculation
                        if (param1.ljid != 0 and dist2 < Rlj2)
                        {
                            const LJPair &ljpair = ljpairs.constData()[
                                        ljpairs.map(param0.ljid,
                                                    param1.ljid)];
                            const double sig2 = ljpair.sigma() * ljpair.sigma();
                            const double sig6 = sig2 * sig2 * sig2;
			                        
                            for (int k=0; k<nalpha; ++k)
                            {
                                const double shift = ljpair.sigma() * delta[k];
			    
                                double lj_denom = dist2 + shift;
                                lj_denom = lj_denom * lj_denom * lj_denom;
			    
                                const double sig6_over_denom = sig6 / lj_denom;
                                const double sig12_over_denom2 = sig6_over_denom *
                                                                 sig6_over_denom;
			    
                                iljnrg[k] +=  cljscl.lj() * ljpair.epsilon() * (sig12_over_denom2 -
                                                                                sig6_over_denom);
                            }
                        }
                    }
                }// for each j
            }//for each i
        }
        else if (use_atomistic_cutoff)
        {
            foreach (Index i, atoms0)
            {
                distmat.setOuterIndex(i);
                const Parameter &param0 = params0_array[i];
	    
                foreach (Index j, atoms1)
                {
                    //do both coulomb and LJ
                    const Parameter &param1 = params1_array[j];
                    const double dist2 = distmat[j];
                        
                    const double q2 = cljscl.coulomb() *
                    param0.reduced_charge * param1.reduced_charge;
                        
                    for (int k=0; k<nalpha; ++k)
                    {
                        //coulomb calculation
                        {
                            const double sr = std::sqrt(alfa[k] + dist2);
                            const double one_over_sr = double(1) / sr;
 
                            // JM Jan 13. No reaction field on 1,4 atoms?
                            if (dist2 < Rcoul2)
                            {
                                icnrg[k] += q2 * (one_over_sr);
                            }
                        }
		  
                        //lj calculation
                        if (param1.ljid != 0 and dist2 < Rlj2)
                        {
                            const LJPair &ljpair = ljpairs.constData()[
                                        ljpairs.map(param0.ljid,
                                                    param1.ljid)];
                            const double sig2 = ljpair.sigma() * ljpair.sigma();
                            const double sig6 = sig2 * sig2 * sig2;
			                        
                            for (int k=0; k<nalpha; ++k)
                            {
                                const double shift = ljpair.sigma() * delta[k];
			    
                                double lj_denom = dist2 + shift;
                                lj_denom = lj_denom * lj_denom * lj_denom;
			    
                                const double sig6_over_denom = sig6 / lj_denom;
                                const double sig12_over_denom2 = sig6_over_denom *
                                                                 sig6_over_denom;
			    
                                iljnrg[k] +=  cljscl.lj() * ljpair.epsilon() * (sig12_over_denom2 -
                                                                                sig6_over_denom);
                            }
                        }
                    }
                }// for each j
            }//for each i
        }
        else //group-based feathered cutoff
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
                        const double dist2 = distmat[j];
		    
                        const double q2 = cljscl.coulomb() *
                            param0.reduced_charge * params1_array[j].reduced_charge;
		    
                        for (int k=0; k<nalpha; ++k)
                        {
                            icnrg[k] += q2 / std::sqrt(alfa[k] + dist2);
                        }
                    }
                }
                else
                {
                    foreach (Index j, atoms1)
                    {
                        //do both coulomb and LJ
                        const Parameter &param1 = params1_array[j];
                        const double dist2 = distmat[j];
                        
                        const double q2 = cljscl.coulomb() *
                            param0.reduced_charge * param1.reduced_charge;
                        
                        for (int k=0; k<nalpha; ++k)
                        {
                            icnrg[k] += q2 / std::sqrt(alfa[k] + dist2);
                        }

                        if (param1.ljid != 0)
                        {
                            const LJPair &ljpair = ljpairs.constData()[
                                                    ljpairs.map(param0.ljid,
                                                                param1.ljid)];
                        
                            const double sig2 = ljpair.sigma() * ljpair.sigma();
                            const double sig6 = sig2 * sig2 * sig2;
			                        
                            for (int k=0; k<nalpha; ++k)
                            {
                                const double shift = ljpair.sigma() * delta[k];
			    
                                double lj_denom = dist2 + shift;
                                lj_denom = lj_denom * lj_denom * lj_denom;
			    
                                const double sig6_over_denom = sig6 / lj_denom;
                                const double sig12_over_denom2 = sig6_over_denom *
                                                                 sig6_over_denom;
			    
                                iljnrg[k] +=  cljscl.lj() * ljpair.epsilon() * (sig12_over_denom2 -
                                                                                sig6_over_denom);
                            }
                        }
                    }
                }
            }
        } 
    }
    else
    {
        //there are different nb scale factors between
        //the atoms. We need to calculate the energies using
        //them...
        if (use_electrostatic_shifting)
        {
            double sRcoul[nalpha];
            double one_over_sRcoul[nalpha];
            double one_over_sRcoul2[nalpha];
        
            for (int i=0; i<nalpha; ++i)
            {
                sRcoul[i] = std::sqrt(alfa[i] + Rcoul*Rcoul);
                one_over_sRcoul[i] = double(1) / sRcoul[i];
                one_over_sRcoul2[i] = double(1) / (sRcoul[i]*sRcoul[i]);
            }
	
            foreach (Index i, atoms0)
            {
                distmat.setOuterIndex(i);
                const Parameter &param0 = params0_array[i];

                if (param0.ljid == 0)
                {
                    foreach (Index j, atoms1)
                    {
                        const CLJScaleFactor &cljscl = group_pairs(i,j);
                            
                        if (cljscl.coulomb() != 0)
                        {
                            const double dist2 = distmat[j];
			
                            const double q2 = cljscl.coulomb() *
                                param0.reduced_charge * params1_array[j].reduced_charge;
			
                            for (int k=0; k<nalpha; ++k)
                            {
                                const double sr = std::sqrt(alfa[k] + dist2);
                                const double one_over_sr = double(1) / sr;
 
                                // JM Jan 13. No reaction field on 1,4 atoms?
                                if (sr < sRcoul[k])
                                {
                                    if (cljscl.coulomb() != 1)
                                        icnrg[k] += q2 * (one_over_sr);
                                    else
                                        icnrg[k] += q2 * (one_over_sr - one_over_sRcoul[k] +
                                                          one_over_sRcoul2[k]*(sr-sRcoul[k]));
                                }
                            }
                        }
                    }//foreach j
                }
                else //do both
                {
                    //do both coulomb and LJ
                    foreach (Index j, atoms1)
                    {
                        const CLJScaleFactor &cljscl = group_pairs(i,j);

                        if (cljscl.coulomb() != 0 or cljscl.lj() != 0)
                        {
                            const Parameter &param1 = params1_array[j];
		    
                            const double dist2 = distmat[j];
			
                            const double q2 = cljscl.coulomb() *
                                param0.reduced_charge * params1_array[j].reduced_charge;
			
                            for (int k=0; k<nalpha; ++k)
                            {
                                const double sr = std::sqrt(alfa[k] + dist2);
                                const double one_over_sr = double(1) / sr;
 
                                // JM Jan 13. No reaction field on 1,4 atoms?
                                if (sr < sRcoul[k])
                                {
                                    if (cljscl.coulomb() != 1)
                                        icnrg[k] += q2 * (one_over_sr);
                                    else
                                        icnrg[k] += q2 * (one_over_sr - one_over_sRcoul[k] +
                                                          one_over_sRcoul2[k]*(sr-sRcoul[k]));
                                }
                            }

                            if (cljscl.lj() != 0 and param1.ljid != 0 and dist2 < Rlj2)
                            {
                                const LJPair &ljpair = ljpairs.constData()[
                                                        ljpairs.map(param0.ljid,
                                                                    param1.ljid)];
                                const double sig2 = ljpair.sigma() * ljpair.sigma();
                                const double sig6 = sig2 * sig2 * sig2;
			                            
                                for (int k=0; k<nalpha; ++k)
                                {
                                    const double shift = ljpair.sigma() * delta[k];
			    
                                    double lj_denom = dist2 + shift;
                                    lj_denom = lj_denom * lj_denom * lj_denom;
			    
                                    const double sig6_over_denom = sig6 / lj_denom;
                                    const double sig12_over_denom2 = sig6_over_denom *
                                                                     sig6_over_denom;
			    
                                    iljnrg[k] +=  cljscl.lj() * ljpair.epsilon() *
                                                    (sig12_over_denom2 - sig6_over_denom);
                                }
                            }
                        }
                    } //for each j
                }
            }//for each i
        }
        else if (use_reaction_field)
        {
            double sRcoul[nalpha];
            double k_rf[nalpha];
            double c_rf[nalpha];
        
            for (int i=0; i<nalpha; ++i)
            {
                sRcoul[i] = std::sqrt(alfa[i] + Rcoul*Rcoul);
                k_rf[i] = (1.0 / pow_3(sRcoul[i])) * ( (rf_dielectric_constant-1) /
                                                       (2*rf_dielectric_constant + 1) );
                c_rf[i] = (1.0 / sRcoul[i]) * ( (3*rf_dielectric_constant) /
                                                (2*rf_dielectric_constant + 1) );
            }
	
            foreach (Index i, atoms0)
            {
                distmat.setOuterIndex(i);
                const Parameter &param0 = params0_array[i];

                if (param0.ljid == 0)
                {
                    foreach (Index j, atoms1)
                    {
                        const CLJScaleFactor &cljscl = group_pairs(i,j);
                            
                        if (cljscl.coulomb() != 0)
                        {
                            const double dist2 = distmat[j];
			
                            const double q2 = cljscl.coulomb() *
                                param0.reduced_charge * params1_array[j].reduced_charge;
			
                            for (int k=0; k<nalpha; ++k)
                            {
                                const double sr = std::sqrt(alfa[k] + dist2);
                                const double one_over_sr = double(1) / sr;
 
                                // JM Jan 13. No reaction field on 1,4 atoms?
                                if (sr < sRcoul[k])
                                {
                                    if (cljscl.coulomb() != 1)
                                        icnrg[k] += q2 * (one_over_sr);
                                    else
                                        icnrg[k] += q2 * (one_over_sr + sr*sr*k_rf[k] - c_rf[k]);
                                }
                            }
                        }
                    }//foreach j
                }
                else //do both
                {
                    //do both coulomb and LJ
                    foreach (Index j, atoms1)
                    {
                        const CLJScaleFactor &cljscl = group_pairs(i,j);

                        if (cljscl.coulomb() != 0 or cljscl.lj() != 0)
                        {
                            const Parameter &param1 = params1_array[j];
		    
                            const double dist2 = distmat[j];
			
                            const double q2 = cljscl.coulomb() *
                                param0.reduced_charge * params1_array[j].reduced_charge;
			
                            for (int k=0; k<nalpha; ++k)
                            {
                                const double sr = std::sqrt(alfa[k] + dist2);
                                const double one_over_sr = double(1) / sr;
 
                                // JM Jan 13. No reaction field on 1,4 atoms?
                                if (sr < sRcoul[k])
                                {
                                    if (cljscl.coulomb() != 1)
                                        icnrg[k] += q2 * (one_over_sr);
                                    else
                                        icnrg[k] += q2 * (one_over_sr + sr*sr*k_rf[k] - c_rf[k]);
                                }
                            }

                            if (cljscl.lj() != 0 and param1.ljid != 0 and dist2 < Rlj2)
                            {
                                const LJPair &ljpair = ljpairs.constData()[
                                                        ljpairs.map(param0.ljid,
                                                                    param1.ljid)];
                                const double sig2 = ljpair.sigma() * ljpair.sigma();
                                const double sig6 = sig2 * sig2 * sig2;
			                            
                                for (int k=0; k<nalpha; ++k)
                                {
                                    const double shift = ljpair.sigma() * delta[k];
			    
                                    double lj_denom = dist2 + shift;
                                    lj_denom = lj_denom * lj_denom * lj_denom;
			    
                                    const double sig6_over_denom = sig6 / lj_denom;
                                    const double sig12_over_denom2 = sig6_over_denom *
                                                                     sig6_over_denom;
			    
                                    iljnrg[k] +=  cljscl.lj() * ljpair.epsilon() *
                                                    (sig12_over_denom2 - sig6_over_denom);
                                }
                            }
                        }
                    } //for each j
                }
            }//for each i
        }
        else if (use_atomistic_cutoff)
        {
            foreach (Index i, atoms0)
            {
                distmat.setOuterIndex(i);
                const Parameter &param0 = params0_array[i];

                if (param0.ljid == 0)
                {
                    foreach (Index j, atoms1)
                    {
                        const CLJScaleFactor &cljscl = group_pairs(i,j);
                            
                        if (cljscl.coulomb() != 0)
                        {
                            const double dist2 = distmat[j];
			
                            const double q2 = cljscl.coulomb() *
                                param0.reduced_charge * params1_array[j].reduced_charge;
			
                            for (int k=0; k<nalpha; ++k)
                            {
                                const double sr = std::sqrt(alfa[k] + dist2);
                                const double one_over_sr = double(1) / sr;
 
                                // JM Jan 13. No reaction field on 1,4 atoms?
                                if (dist2 < Rcoul2)
                                {
                                    icnrg[k] += q2 * (one_over_sr);
                                }
                            }
                        }
                    }//foreach j
                }
                else //do both
                {
                    //do both coulomb and LJ
                    foreach (Index j, atoms1)
                    {
                        const CLJScaleFactor &cljscl = group_pairs(i,j);

                        if (cljscl.coulomb() != 0 or cljscl.lj() != 0)
                        {
                            const Parameter &param1 = params1_array[j];
		    
                            const double dist2 = distmat[j];
			
                            const double q2 = cljscl.coulomb() *
                                param0.reduced_charge * params1_array[j].reduced_charge;
			
                            for (int k=0; k<nalpha; ++k)
                            {
                                const double sr = std::sqrt(alfa[k] + dist2);
                                const double one_over_sr = double(1) / sr;
 
                                // JM Jan 13. No reaction field on 1,4 atoms?
                                if (dist2 < Rcoul2)
                                {
                                    icnrg[k] += q2 * (one_over_sr);
                                }
                            }

                            if (cljscl.lj() != 0 and param1.ljid != 0 and dist2 < Rlj2)
                            {
                                const LJPair &ljpair = ljpairs.constData()[
                                                        ljpairs.map(param0.ljid,
                                                                    param1.ljid)];
                                const double sig2 = ljpair.sigma() * ljpair.sigma();
                                const double sig6 = sig2 * sig2 * sig2;
			                            
                                for (int k=0; k<nalpha; ++k)
                                {
                                    const double shift = ljpair.sigma() * delta[k];
			    
                                    double lj_denom = dist2 + shift;
                                    lj_denom = lj_denom * lj_denom * lj_denom;
			    
                                    const double sig6_over_denom = sig6 / lj_denom;
                                    const double sig12_over_denom2 = sig6_over_denom *
                                                                     sig6_over_denom;
			    
                                    iljnrg[k] +=  cljscl.lj() * ljpair.epsilon() *
                                                    (sig12_over_denom2 - sig6_over_denom);
                                }
                            }
                        }
                    } //for each j
                }
            }//for each i
        }
        else // use the group-based feathered cutoff
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
                            const double dist2 = distmat[j];
			
                            const double q2 = cljscl.coulomb() *
                                param0.reduced_charge * params1_array[j].reduced_charge;
			
                            for (int k=0; k<nalpha; ++k)
                            {
                                icnrg[k] += q2 / std::sqrt(alfa[k] + dist2);
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

                            const double dist2 = distmat[j];
			
                            const double q2 = cljscl.coulomb() *
                                param0.reduced_charge * params1_array[j].reduced_charge;
			
                            for (int k=0; k<nalpha; ++k)
                            {
                                icnrg[k] += q2 / std::sqrt(alfa[k] + dist2);
                            }

                            if (cljscl.lj() != 0 and param1.ljid != 0)
                            {
                                const LJPair &ljpair = ljpairs.constData()[
                                                        ljpairs.map(param0.ljid,
                                                                    param1.ljid)];
                                const double sig2 = ljpair.sigma() * ljpair.sigma();
                                const double sig6 = sig2 * sig2 * sig2;
			                            
                                for (int k=0; k<nalpha; ++k)
                                {
                                    const double shift = ljpair.sigma() * delta[k];
			        
                                    double lj_denom = dist2 + shift;
                                    lj_denom = lj_denom * lj_denom * lj_denom;
			        
                                    const double sig6_over_denom = sig6 / lj_denom;
                                    const double sig12_over_denom2 = sig6_over_denom *
                                                                     sig6_over_denom;
				
                                    iljnrg[k] +=  cljscl.lj() * ljpair.epsilon() *
                                                        (sig12_over_denom2 - sig6_over_denom);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

/** Calculate the intramolecular CLJ energy of the passed molecule, and
    add this onto 'energy'. This uses the passed workspace when
    performing the calculation */
void IntraSoftCLJPotential::calculateEnergy(const IntraSoftCLJPotential::Molecule &mol,
					    IntraSoftCLJPotential::Energy &energy,
					    IntraSoftCLJPotential::EnergyWorkspace &distmat,
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

    //the alpha_values array contains all of the unique alpha values
    const double *alfa = alpha_values.constData();
    const int nalpha = alpha_values.count();
        
    if (nalpha <= 0)
        return;
        
    double cnrg[nalpha];
    double ljnrg[nalpha];
        
    for (int i=0; i<nalpha; ++i)
    {
        cnrg[i] = 0;
        ljnrg[i] = 0;
    }
        
    //this uses the following potentials
    //           Zacharias and McCammon, J. Chem. Phys., 1994, and also,
    //           Michel et al., JCTC, 2007
    //
    //  V_{LJ}(r) = 4 epsilon [ ( sigma^12 / (delta*sigma + r^2)^6 ) - 
    //                          ( sigma^6  / (delta*sigma + r^2)^3 ) ]
    //
    //  delta = shift_delta * alpha
    //
    //  V_{coul}(r) = (1-alpha)^n q_i q_j / 4 pi eps_0 (alpha+r^2)^(1/2)
    //
    // This contrasts to Rich T's LJ softcore function, which was;
    //
    //  V_{LJ}(r) = 4 epsilon [ (sigma^12 / (alpha^m sigma^6 + r^6)^2) - 
    //                          (sigma^6  / (alpha^m sigma^6 + r^6) ) ]
    
    double one_minus_alfa_to_n[nalpha];
    double delta[nalpha];
        
    for (int i=0; i<nalpha; ++i)
    {
        one_minus_alfa_to_n[i] = SireMaths::pow(1 - alfa[i], int(coul_power));
        delta[i] = shift_delta * alfa[i];
    }
      
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
            
            //calculate all of the interatomic distances^2
            const double mindist = spce->calcDist2(group0, group1, distmat);
            
            if (mindist > switchfunc->cutoffDistance())
                //all of the atoms are definitely beyond cutoff
                continue;
                
            CGIdx cgidx_jgroup = mol.cgIdx(jgroup);
                
            //get the non-bonded scale factors for all pairs of atoms
            //between these groups (or within this group, if igroup == jgroup)
            const CLJNBPairs::CGPairs &group_pairs = nbpairs(cgidx_igroup,
                                                             cgidx_jgroup);

            double icnrg[nalpha];
            double iljnrg[nalpha];

            for (int i=0; i<nalpha; ++i)
            {
                icnrg[i] = 0;
                iljnrg[i] = 0;
            }
            
            //loop over all intraatomic pairs and calculate the energies
            const quint32 nats1 = group1.count();
            const Parameter *params1_array = params1.constData();
            
            _pvt_calculateEnergy(group_pairs, distmat, params0_array, params1_array,
                                 nats0, nats1, icnrg, iljnrg, alfa, delta, nalpha);
            
            //if this is the same group then half the energies to 
            //correct for double-counting
            if (igroup == jgroup)
            {
                for (int i=0; i<nalpha; ++i)
                {
                    icnrg[i] *= 0.5;
                    iljnrg[i] *= 0.5;
                }
            }

            //now add these energies onto the total for the molecule,
            //scaled by any non-bonded feather factor
            if (not (use_electrostatic_shifting or use_reaction_field or use_atomistic_cutoff))
            {
                if (mindist > switchfunc->electrostaticFeatherDistance())
                {
                    const double cscl = switchfunc->electrostaticScaleFactor(Length(mindist));
                    
                    for (int i=0; i<nalpha; ++i)
                    {
                        icnrg[i] *= cscl;
                    }
                }
                
                if (mindist > switchfunc->vdwFeatherDistance())
                {
                    const double ljscl = switchfunc->vdwScaleFactor(Length(mindist));
                    
                    for (int i=0; i<nalpha; ++i)
                    {
                        iljnrg[i] *= ljscl;
                    }
                }
            }

            for (int i=0; i<nalpha; ++i)
            {
                cnrg[i] += icnrg[i];
                ljnrg[i] += iljnrg[i];
            }
        }
    }
    
    //add this molecule pair's energy onto the total
    //energy += Energy(scale_energy * cnrg, scale_energy * ljnrg);
    for (int i=0; i<nalpha; ++i)
    {
        cnrg[i] *= scale_energy * one_minus_alfa_to_n[i];
        ljnrg[i] *= 4 * scale_energy;
    }

    //now copy the calculated energies back to the Energy object
    Energy soft_energy;
    
    for (int i=0; i<alpha_index.count(); ++i)
    {
        int idx = alpha_index.at(i);
      
        if (idx >= 0)
            soft_energy.setEnergy(i, cnrg[idx], ljnrg[idx]);
    }
    
    energy += soft_energy;
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
void IntraSoftCLJPotential::calculateEnergy(const IntraSoftCLJPotential::Molecule &mol,
					    const IntraSoftCLJPotential::Molecule &rest_of_mol,
					    IntraSoftCLJPotential::Energy &energy,
					    IntraSoftCLJPotential::EnergyWorkspace &distmat,
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
    
    //the alpha_values array contains all of the unique alpha values
    const double *alfa = alpha_values.constData();
    const int nalpha = alpha_values.count();
    
    if (nalpha <= 0)
        return;
    
    double cnrg[nalpha];
    double ljnrg[nalpha];
    
    for (int i=0; i<nalpha; ++i)
    {
        cnrg[i] = 0;
        ljnrg[i] = 0;
    }

    //this uses the following potentials
    //           Zacharias and McCammon, J. Chem. Phys., 1994, and also,
    //           Michel et al., JCTC, 2007
    //
    //  V_{LJ}(r) = 4 epsilon [ ( sigma^12 / (delta*sigma + r^2)^6 ) - 
    //                          ( sigma^6  / (delta*sigma + r^2)^3 ) ]
    //
    //  delta = shift_delta * alpha
    //
    //  V_{coul}(r) = (1-alpha)^n q_i q_j / 4 pi eps_0 (alpha+r^2)^(1/2)
    //
    // This contrasts to Rich T's LJ softcore function, which was;
    //
    //  V_{LJ}(r) = 4 epsilon [ (sigma^12 / (alpha^m sigma^6 + r^6)^2) - 
    //                          (sigma^6  / (alpha^m sigma^6 + r^6) ) ]
    //

    double one_minus_alfa_to_n[nalpha];
    double delta[nalpha];
        
    for (int i=0; i<nalpha; ++i)
    {
        one_minus_alfa_to_n[i] = SireMaths::pow(1 - alfa[i], int(coul_power));
        delta[i] = shift_delta * alfa[i];
    }

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
            
            //calculate all of the interatomic distances^2
            const double mindist = spce->calcDist2(group0, group1, distmat);
            
            if (mindist > switchfunc->cutoffDistance())
                //all of the atoms are definitely beyond cutoff
                continue;
                
            //get the non-bonded scale factors for all pairs of atoms
            //between these groups (or within this group, if igroup == jgroup)
            const CLJNBPairs::CGPairs &group_pairs = nbpairs(cgidx_igroup,
                                                             cgidx_jgroup);

            double icnrg[nalpha];
            double iljnrg[nalpha];
	                
            for (int i=0; i<nalpha; ++i)
            {
                icnrg[i] = 0;
                iljnrg[i] = 0;
            }
            
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
                
                _pvt_calculateEnergy(group_pairs, atoms0, atoms1, distmat,
				     params0_array, params1_array, 
				     nats0, nats1, icnrg, iljnrg, alfa, delta, nalpha);
            }
            else
            {
                _pvt_calculateEnergy(group_pairs, distmat,
                    params0_array, params1_array,
                    nats0, nats1, icnrg, iljnrg, alfa, delta, nalpha);
            }

            //now add these energies onto the total for the molecule,
            //scaled by any non-bonded feather factor
            //now add these energies onto the total for the molecule,
            //scaled by any non-bonded feather factor
            if (not (use_electrostatic_shifting or use_reaction_field or use_atomistic_cutoff))
            {
                if (mindist > switchfunc->electrostaticFeatherDistance())
                {
                    const double cscl = switchfunc->electrostaticScaleFactor(Length(mindist));
                    
                    for (int i=0; i<nalpha; ++i)
                    {
                        icnrg[i] *= cscl;
                    }
                }
                
                if (mindist > switchfunc->vdwFeatherDistance())
                {
                    const double ljscl = switchfunc->vdwScaleFactor(Length(mindist));
                    
                    for (int i=0; i<nalpha; ++i)
                    {
                        iljnrg[i] *= ljscl;
                    }
                }
            }

            for (int i=0; i<nalpha; ++i)
            {
                cnrg[i] += icnrg[i];
                ljnrg[i] += iljnrg[i];
            }
        }
    }
    
    //add the molecule's energy onto the total
    //energy += Energy(scale_energy * cnrg, scale_energy * ljnrg);
    for (int i=0; i<nalpha; ++i)
    {
        cnrg[i] *= scale_energy * one_minus_alfa_to_n[i];
        ljnrg[i] *= 4 * scale_energy;
    }
    
    //now copy the calculated energies back to the Energy object
    Energy soft_energy;
        
    for (int i=0; i<alpha_index.count(); ++i)
    {
        int idx = alpha_index.at(i);
        
        if (idx >= 0)
            soft_energy.setEnergy(i, cnrg[idx], ljnrg[idx]);
    }
    
    energy += soft_energy;
}

/** Calculate the coulomb and LJ forces between the atoms in the molecule 'mol'
    and add these forces onto 'forces'. This uses
    the passed workspace to perform the calculation */
void IntraSoftCLJPotential::calculateForce(const IntraSoftCLJPotential::Molecule &mol,
					   MolForceTable &forces, 
					   IntraSoftCLJPotential::ForceWorkspace &distmat,
					   double scale_force) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code necessary to calculate intramolecular soft coulomb "
                "and LJ forces has not yet been written..."), CODELOC );
}

/** Calculate the total forces acting on the atoms in 'mol' caused by the 
    other atoms in the same molecule contained in 'rest_of_mol'. This calculates
    the forces and adds them onto 'forces' (which are for 'mol'). Note that they must
    use the same layout UID and same intra-nonbonded scaling factors
    
    \throw SireError::incompatible_error
*/
void IntraSoftCLJPotential::calculateForce(const IntraSoftCLJPotential::Molecule &mol,
					   const IntraSoftCLJPotential::Molecule &rest_of_mol,
					   MolForceTable &forces,
					   IntraSoftCLJPotential::ForceWorkspace &distmat,
					   double scale_force) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code necessary to calculate intramolecular soft coulomb "
                "and LJ forces has not yet been written..."), CODELOC );
}

void IntraSoftCLJPotential::calculateField(const IntraSoftCLJPotential::Molecule &mol, 
					   const CLJProbe &probe,
					   MolFieldTable &fields,
					   IntraSoftCLJPotential::FieldWorkspace &workspace,
					   double scale_field) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code necessary to calculate intramolecular soft coulomb "
                "and LJ fields has not yet been written..."), CODELOC );
}

void IntraSoftCLJPotential::calculateField(const IntraSoftCLJPotential::Molecule &mol,
					   const IntraSoftCLJPotential::Molecule &rest_of_mol,
					   const CLJProbe &probe,
					   MolFieldTable &fields,
					   IntraSoftCLJPotential::FieldWorkspace &workspace,
					   double scale_field) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code necessary to calculate intramolecular soft coulomb "
                "and LJ fields has not yet been written..."), CODELOC );
}

void IntraSoftCLJPotential::calculateField(const IntraSoftCLJPotential::Molecule &mol, 
					   const CLJProbe &probe,
					   MolFieldTable &fields,
					   const Symbol &symbol,
					   const Components &components,
					   IntraSoftCLJPotential::FieldWorkspace &workspace,
					   double scale_field) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code necessary to calculate intramolecular soft coulomb "
                "and LJ fields has not yet been written..."), CODELOC );
}

void IntraSoftCLJPotential::calculateField(const IntraSoftCLJPotential::Molecule &mol,
					   const IntraSoftCLJPotential::Molecule &rest_of_mol,
					   const CLJProbe &probe,
					   MolFieldTable &fields,
					   const Symbol &symbol,
					   const Components &components,
					   IntraCLJPotential::FieldWorkspace &workspace,
					   double scale_field) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code necessary to calculate intramolecular soft coulomb "
                "and LJ fields has not yet been written..."), CODELOC );
}

void IntraSoftCLJPotential::calculateField(const IntraSoftCLJPotential::Molecule &mol, 
					   const CLJProbe &probe,
					   GridFieldTable &fields,
					   IntraSoftCLJPotential::FieldWorkspace &workspace,
					   double scale_field) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code necessary to calculate intramolecular soft coulomb "
                "and LJ fields has not yet been written..."), CODELOC );
}

void IntraSoftCLJPotential::calculateField(const IntraSoftCLJPotential::Molecule &mol, 
					   const CLJProbe &probe,
					   GridFieldTable &fields,
					   const Symbol &symbol,
					   const Components &components,
					   IntraSoftCLJPotential::FieldWorkspace &workspace,
                    double scale_field) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code necessary to calculate intramolecular soft coulomb "
                "and LJ fields has not yet been written..."), CODELOC );
}

void IntraSoftCLJPotential::calculatePotential(const IntraSoftCLJPotential::Molecule &mol, 
					       const CLJProbe &probe,
					       MolPotentialTable &potentials,
					       IntraSoftCLJPotential::PotentialWorkspace &workspace,
					       double scale_potential) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code necessary to calculate intramolecular soft coulomb "
                "and LJ potentials has not yet been written..."), CODELOC );
}

void IntraSoftCLJPotential::calculatePotential(const IntraSoftCLJPotential::Molecule &mol,
					       const IntraSoftCLJPotential::Molecule &rest_of_mol,
					       const CLJProbe &probe,
					       MolPotentialTable &potentials,
					       IntraSoftCLJPotential::PotentialWorkspace &workspace,
					       double scale_potential) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code necessary to calculate intramolecular soft coulomb "
                "and LJ fields has not yet been written..."), CODELOC );
}

void IntraSoftCLJPotential::calculatePotential(const IntraSoftCLJPotential::Molecule &mol, 
					       const CLJProbe &probe,
					       MolPotentialTable &potentials,
					       const Symbol &symbol,
					       const Components &components,
					       IntraSoftCLJPotential::PotentialWorkspace &workspace,
					       double scale_potential) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code necessary to calculate intramolecular soft coulomb "
                "and LJ potentials has not yet been written..."), CODELOC );
}

void IntraSoftCLJPotential::calculatePotential(const IntraSoftCLJPotential::Molecule &mol,
                        const IntraSoftCLJPotential::Molecule &rest_of_mol,
                        const CLJProbe &probe,
                        MolPotentialTable &potentials,
                        const Symbol &symbol,
                        const Components &components,
                        IntraSoftCLJPotential::PotentialWorkspace &workspace,
                        double scale_potential) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code necessary to calculate intramolecular soft coulomb "
                "and LJ potentials has not yet been written..."), CODELOC );
}

void IntraSoftCLJPotential::calculatePotential(const IntraSoftCLJPotential::Molecule &mol, 
                        const CLJProbe &probe,
                        GridPotentialTable &potentials,
                        IntraSoftCLJPotential::PotentialWorkspace &workspace,
                        double scale_potential) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code necessary to calculate intramolecular soft coulomb "
                "and LJ potentials has not yet been written..."), CODELOC );
}

void IntraSoftCLJPotential::calculatePotential(const IntraSoftCLJPotential::Molecule &mol,
                        const CLJProbe &probe,
                        GridPotentialTable &potentials,
                        const Symbol &symbol,
                        const Components &components,
                        IntraSoftCLJPotential::PotentialWorkspace &workspace,
                        double scale_potential) const
{
    throw SireError::incomplete_code( QObject::tr(
                "The code necessary to calculate intramolecular soft coulomb "
                "and LJ potentials has not yet been written..."), CODELOC );
}
