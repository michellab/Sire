/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2009  Christopher Woods
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

#include "ljperturbation.h"
#include "atomljs.h"

#include "SireMol/molecule.h"
#include "SireMol/moleditor.h"
#include "SireMol/mover.hpp"

#include "SireCAS/values.h"
#include "SireCAS/identities.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMM;
using namespace SireMol;
using namespace SireCAS;
using namespace SireBase;
using namespace SireUnits::Dimension;
using namespace SireStream;

static const RegisterMetaType<LJPerturbation> r_ljpert;

QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds,
                                       const LJPerturbation &ljpert)
{
    writeHeader(ds, r_ljpert, 1);
    
    SharedDataStream sds(ds);
    
    sds << ljpert.sigma_mapfunc << quint32(ljpert.maptype)
        << static_cast<const Perturbation&>(ljpert);
    
    return ds;
}

QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds,
                                       LJPerturbation &ljpert)
{
    VersionID v = readHeader(ds, r_ljpert);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        quint32 maptype;
        sds >> ljpert.sigma_mapfunc >> maptype >> static_cast<Perturbation&>(ljpert);
        
        switch (maptype)
        {
            case LJPerturbation::MAP_SIGMA_AND_EPSILON:
                ljpert.maptype = LJPerturbation::MAP_SIGMA_AND_EPSILON;
                break;
            case LJPerturbation::MAP_RMIN_AND_EPSILON:
                ljpert.maptype = LJPerturbation::MAP_RMIN_AND_EPSILON;
                break;
            case LJPerturbation::MAP_A_AND_B:
                ljpert.maptype = LJPerturbation::MAP_A_AND_B;
                break;
            default:
                ljpert.maptype = LJPerturbation::MAP_SIGMA_AND_EPSILON;
        }
    }
    else
        throw version_error(v, "1", r_ljpert, CODELOC);
        
    return ds;
}

/** Constructor - this creates a LJ perturbation that 
    perturbs from LJs in "initial_LJ" to LJs in
    "final_LJ", placing the current LJs in "LJ",
    and using Perturbation::defaultEquation() to map the
    sigma and epsilon values of the LJ.  */
LJPerturbation::LJPerturbation(const PropertyMap &map)
               : ConcreteProperty<LJPerturbation,Perturbation>(map),
                 sigma_mapfunc( Perturbation::defaultFunction() ),
                 maptype( MAP_SIGMA_AND_EPSILON )
{}

/** Construct, using the passed map to find the properties used
    by this perturbation */
LJPerturbation::LJPerturbation(MapType typ, const PropertyMap &map)
               : ConcreteProperty<LJPerturbation,Perturbation>(map),
                 sigma_mapfunc( Perturbation::defaultFunction() ),
                 maptype(typ)
{}

/** Construct, using the passed map to find the properties used
    by this perturbation and the passed mapping function to map
    the LJs between the states */
LJPerturbation::LJPerturbation(const Expression &mapping_function,
                               const PropertyMap &map)
               : ConcreteProperty<LJPerturbation,Perturbation>(mapping_function, map),
                 sigma_mapfunc(mapping_function),
                 maptype( MAP_SIGMA_AND_EPSILON )
{}

/** Construct, using the passed map to find the properties used
    by this perturbation and the passed mapping function to map
    the LJs between the states */
LJPerturbation::LJPerturbation(const Expression &mapping_function, MapType typ,
                               const PropertyMap &map)
               : ConcreteProperty<LJPerturbation,Perturbation>(mapping_function, map),
                 sigma_mapfunc(mapping_function),
                 maptype(typ)
{}

/** Construct, using the passed map to find the properties used
    by this perturbation and the passed mapping function to map
    the LJs between the states */
LJPerturbation::LJPerturbation(const Expression &sigma_function,
                               const Expression &epsilon_function,
                               const PropertyMap &map)
               : ConcreteProperty<LJPerturbation,Perturbation>(epsilon_function, map),
                 sigma_mapfunc(sigma_function),
                 maptype( MAP_SIGMA_AND_EPSILON )
{}

/** Construct, using the passed map to find the properties used
    by this perturbation and the passed mapping function to map
    the LJs between the states */
LJPerturbation::LJPerturbation(const Expression &sigma_function, 
                               const Expression &epsilon_function,
                               MapType typ,
                               const PropertyMap &map)
               : ConcreteProperty<LJPerturbation,Perturbation>(epsilon_function, map),
                 sigma_mapfunc(sigma_function), maptype(typ)
{}

/** Copy constructor */
LJPerturbation::LJPerturbation(const LJPerturbation &other)
               : ConcreteProperty<LJPerturbation,Perturbation>(other),
                 sigma_mapfunc(other.sigma_mapfunc), maptype(other.maptype)
{}

/** Destructor */
LJPerturbation::~LJPerturbation()
{}

const char* LJPerturbation::typeName()
{
    return QMetaType::typeName( qMetaTypeId<LJPerturbation>() );
}

/** Copy assignment operator */
LJPerturbation& LJPerturbation::operator=(const LJPerturbation &other)
{
    if (this != &other)
    {
        sigma_mapfunc = other.sigma_mapfunc;
        maptype = other.maptype;
        Perturbation::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool LJPerturbation::operator==(const LJPerturbation &other) const
{
    return sigma_mapfunc == other.sigma_mapfunc and maptype == other.maptype and
           Perturbation::operator==(other);
}

/** Comparison operator */
bool LJPerturbation::operator!=(const LJPerturbation &other) const
{
    return not LJPerturbation::operator==(other);
}

PerturbationPtr LJPerturbation::recreate(const Expression &mapping_function) const
{
    PerturbationPtr ret = Perturbation::recreate(mapping_function);
    ret.edit().asA<LJPerturbation>().sigma_mapfunc = mapping_function;
    
    return ret;
}

PerturbationPtr LJPerturbation::recreate(const Expression &mapping_function,
                                         const PropertyMap &map) const
{
    PerturbationPtr ret = Perturbation::recreate(mapping_function, map);
    ret.edit().asA<LJPerturbation>().sigma_mapfunc = mapping_function;
    
    return ret;
}

/** Substitute the identities in 'identities' in all of the mapping functions 
    used by this perturbation. This is useful if, for example, you want to 
    switch from using 'lambda' to control the perturbation to using 'alpha', e.g.
    
    alpha_perturbations = lambda_perturbations.substitute( lam == Expression(alpha) );
*/
PerturbationPtr LJPerturbation::substitute(const Identities &identities) const
{
    PerturbationPtr ret = Perturbation::substitute(identities);
    
    ret.edit().asA<LJPerturbation>().sigma_mapfunc = sigma_mapfunc.substitute(identities);

    return ret;
}

/** Return whether or not this maps sigma and epsilon */
bool LJPerturbation::mapSigmaEpsilon() const
{
    return maptype == MAP_SIGMA_AND_EPSILON;
}

/** Return whether or not this maps r_min and epsilon */
bool LJPerturbation::mapRMinEpsilon() const
{
    return maptype == MAP_RMIN_AND_EPSILON;
}

/** Return whether or not this maps A and B */
bool LJPerturbation::mapAB() const
{
    return maptype == MAP_A_AND_B;
}

/** Return a string representation of this perturbation */
QString LJPerturbation::toString() const
{
    if (sigma_mapfunc == Perturbation::mappingFunction())
    {
        switch (maptype)
        {
            case MAP_SIGMA_AND_EPSILON:
                return QObject::tr("LJPerturbation( sigma+epsilon => %1 )")
                            .arg(sigma_mapfunc.toString());
            case MAP_RMIN_AND_EPSILON:
                return QObject::tr("LJPerturbation( r_min+epsilon => %1 )")
                            .arg(sigma_mapfunc.toString());
            case MAP_A_AND_B:
                return QObject::tr("LJPerturbation( A+B => %1 )")
                            .arg(sigma_mapfunc.toString());
        }
    }
    else
    {
        switch (maptype)
        {
            case MAP_SIGMA_AND_EPSILON:
                return QObject::tr("LJPerturbation( sigma => %1, epsilon => %2 )")
                            .arg(sigma_mapfunc.toString(), 
                                 Perturbation::mappingFunction().toString());
            case MAP_RMIN_AND_EPSILON:
                return QObject::tr("LJPerturbation( r_min => %1, epsilon => %2 )")
                            .arg(sigma_mapfunc.toString(), 
                                 Perturbation::mappingFunction().toString());
            case MAP_A_AND_B:
                return QObject::tr("LJPerturbation( A => %1, B => %2 )")
                            .arg(sigma_mapfunc.toString(), 
                                 Perturbation::mappingFunction().toString());
        }
    }
    
    return QObject::tr( "LJParameter( ??? )" );
}

/** Return the mapping function

    \throw SireError::invalid_state
*/
const Expression& LJPerturbation::mappingFunction() const
{
    if ( Perturbation::mappingFunction() != sigma_mapfunc )
        throw SireError::invalid_state( QObject::tr(
            "This LJ perturbation uses different functions to map the "
            "two parts of the LJ parameter (%1 and %2). It is thus not "
            "possible to return a single mapping function.")
                .arg(sigma_mapfunc.toString(),
                     Perturbation::mappingFunction().toString()), CODELOC );

    return sigma_mapfunc;
}

/** Return the function used to map r_min

    \throw SireError::invalid_state
*/
const Expression& LJPerturbation::rMinMappingFunction() const
{
    if (not mapRMinEpsilon())
        throw SireError::invalid_state( QObject::tr(
            "This LJ perturbation (%1) does not map r_min.")
                .arg(toString()), CODELOC );

    return sigma_mapfunc;
}

/** Return the function used to map sigma 

    \throw SireError::invalid_state
*/
const Expression& LJPerturbation::sigmaMappingFunction() const
{
    if (not mapSigmaEpsilon())
        throw SireError::invalid_state( QObject::tr(
            "This LJ perturbation (%1) does not map sigma.")
                .arg(toString()), CODELOC );

    return sigma_mapfunc;
}

/** Return the function used to map epsilon

    \throw SireError::invalid_state
*/
const Expression& LJPerturbation::epsilonMappingFunction() const
{
    if ( not (mapRMinEpsilon() or mapSigmaEpsilon()) )
        throw SireError::invalid_state( QObject::tr(
            "This LJ perturbation (%1) does not map epsilon.")
                .arg(toString()), CODELOC );

    return Perturbation::mappingFunction();
}

/** Return the function used to map 'A'

    \throw SireError::invalid_state
*/
const Expression& LJPerturbation::A_MappingFunction() const
{
    if (not mapAB())
        throw SireError::invalid_state( QObject::tr(
            "This LJ perturbation (%1) does not map A.")
                .arg(toString()), CODELOC );

    return sigma_mapfunc;
}

/** Return the function used to map 'B'

    \throw SireError::invalid_state
*/
const Expression& LJPerturbation::B_MappingFunction() const
{
    if (not mapAB())
        throw SireError::invalid_state( QObject::tr(
            "This LJ perturbation (%1) does not map B.")
                .arg(toString()), CODELOC );

    return Perturbation::mappingFunction();
}

/** Return the properties required or changed by this perturbation */
QSet<QString> LJPerturbation::requiredProperties() const
{
    QSet<QString> props;
    
    PropertyName prop = propertyMap()["LJ"];
    
    if (prop.hasSource())
        props.insert( prop.source() );
        
    prop = propertyMap()["initial_LJ"];
    
    if (prop.hasSource())
        props.insert( prop.source() );
        
    prop = propertyMap()["final_LJ"];
    
    if (prop.hasSource())
        props.insert( prop.source() );
        
    return props;
}

/** Return whether or not this perturbation with the passed values would
    change the molecule 'molecule' */
bool LJPerturbation::wouldChange(const Molecule &molecule, 
                                 const Values &values) const
{
    try
    {
        const AtomLJs &initial_ljs = molecule.property( 
                                            propertyMap()["initial_LJ"] )
                                                .asA<AtomLJs>();
                                           
        const AtomLJs &final_ljs = molecule.property( 
                                            propertyMap()["final_LJ"] )
                                                .asA<AtomLJs>();

        const AtomLJs &ljs = molecule.property( 
                                            propertyMap()["LJ"] )
                                                .asA<AtomLJs>();
                                                
        const Expression &f0 = sigma_mapfunc;
        const Expression &f1 = Perturbation::mappingFunction();
        const Symbol &initial = this->symbols().initial();
        const Symbol &final = this->symbols().final();
    
        for (CGIdx i(0); i<initial_ljs.nCutGroups(); ++i)
        {
            for (Index j(0); j<initial_ljs.nAtoms(i); ++j)
            {
                CGAtomIdx atomidx(i,j);

                const LJParameter &initial_lj = initial_ljs[atomidx];
                const LJParameter &final_lj = final_ljs[atomidx];
                const LJParameter &lj = ljs[atomidx];

                if (initial_lj != final_lj)
                {
                    if ( mapSigmaEpsilon() )
                    {
                        Values atom_values = values + 
                                            (initial == initial_lj.sigma().value()) +
                                            (final == final_lj.sigma().value());
            
                        double new_sigma = f0(atom_values);
                
                        atom_values = values +
                                    (initial == initial_lj.epsilon().value()) +
                                    (final == final_lj.epsilon().value());
                                
                        double new_epsilon = f1(atom_values);
            
                        if (lj != LJParameter(Length(new_sigma),MolarEnergy(new_epsilon)))
                            return true;
                    }
                    else if ( mapAB() )
                    {
                        Values atom_values = values + 
                                            (initial == initial_lj.A()) +
                                            (final == final_lj.A());
            
                        double new_A = f0(atom_values);
                
                        atom_values = values +
                                    (initial == initial_lj.B()) +
                                    (final == final_lj.B());
                                
                        double new_B = f1(atom_values);
            
                        if (lj != LJParameter::fromAAndB(new_A, new_B))
                            return true;
                    }
                    else
                    {
                        Values atom_values = values + 
                                            (initial == initial_lj.rmin().value()) +
                                            (final == final_lj.rmin().value());
            
                        double new_rmin = f0(atom_values);
                
                        atom_values = values +
                                    (initial == initial_lj.epsilon().value()) +
                                    (final == final_lj.epsilon().value());
                                
                        double new_epsilon = f1(atom_values);
            
                        if (lj != LJParameter::fromRMinAndEpsilon(Length(new_rmin), 
                                                           MolarEnergy(new_epsilon)))
                        {
                            return true;
                        }
                    }
                }
                else if (initial_lj != lj)
                {
                    return true;
                }
            }
        }
        
        return false;
    }
    catch(...)
    {
        return false;
    }
}

/** Perturb the LJs in the passed molecule using the reaction
    coordinate(s) in 'values' 
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void LJPerturbation::perturbMolecule(MolEditor &molecule, const Values &values) const
{
    const AtomLJs &initial_ljs = molecule.property( propertyMap()["initial_LJ"] )
                                         .asA<AtomLJs>();
                                           
    const AtomLJs &final_ljs = molecule.property( propertyMap()["final_LJ"] )
                                       .asA<AtomLJs>();
                                            
    AtomLJs ljs(initial_ljs);
    
    const Expression &f0 = sigma_mapfunc;
    const Expression &f1 = Perturbation::mappingFunction();
    const Symbol &initial = this->symbols().initial();
    const Symbol &final = this->symbols().final();
    
    for (CGIdx i(0); i<initial_ljs.nCutGroups(); ++i)
    {
        for (Index j(0); j<initial_ljs.nAtoms(i); ++j)
        {
            CGAtomIdx atomidx(i,j);

            const LJParameter &initial_lj = initial_ljs[atomidx];
            const LJParameter &final_lj = final_ljs[atomidx];

            if (initial_lj != final_lj)
            {
                if ( mapSigmaEpsilon() )
                {
                    Values atom_values = values + 
                                        (initial == initial_lj.sigma().value()) +
                                        (final == final_lj.sigma().value());
        
                    double new_sigma = f0(atom_values);
            
                    atom_values = values +
                                (initial == initial_lj.epsilon().value()) +
                                (final == final_lj.epsilon().value());
                            
                    double new_epsilon = f1(atom_values);

                    ljs.set( atomidx, LJParameter::fromSigmaAndEpsilon(Length(new_sigma), 
                                                            MolarEnergy(new_epsilon)) );
                }
                else if ( mapAB() )
                {
                    Values atom_values = values + 
                                        (initial == initial_lj.A()) +
                                        (final == final_lj.A());
        
                    double new_A = f0(atom_values);
            
                    atom_values = values +
                                (initial == initial_lj.B()) +
                                (final == final_lj.B());
                            
                    double new_B = f1(atom_values);
        
                    ljs.set( atomidx, LJParameter::fromAAndB(new_A, new_B) );
                }
                else
                {
                    Values atom_values = values + 
                                        (initial == initial_lj.rmin().value()) +
                                        (final == final_lj.rmin().value());
        
                    double new_rmin = f0(atom_values);
            
                    atom_values = values +
                                (initial == initial_lj.epsilon().value()) +
                                (final == final_lj.epsilon().value());
                            
                    double new_epsilon = f1(atom_values);
        
                    ljs.set( atomidx, LJParameter::fromRMinAndEpsilon(Length(new_rmin), 
                                                            MolarEnergy(new_epsilon)) );
                }
            }
        }
    }
    
    molecule.setProperty( propertyMap()["LJ"].source(), ljs );
}
