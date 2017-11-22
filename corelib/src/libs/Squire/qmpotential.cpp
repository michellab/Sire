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

#include "qmpotential.h"
#include "qmprogram.h"

#include "SireFF/ff.h"

#include "SireMol/atomelements.h"

#include "SireBase/propertylist.h"

#include "SireError/errors.h"
#include "SireFF/errors.h"
#include "SireBase/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace Squire;
using namespace SireFF;
using namespace SireFF::detail;
using namespace SireMol;
using namespace SireBase;
using namespace SireUnits::Dimension;
using namespace SireStream;

///////
/////// Completely instantiate the QMPotential ancillary classes
///////

namespace SireFF
{
    namespace detail
    {
        template class AtomicParameters3D<Element>;
        template class FFMolecule3D<QMPotential>;
        template class FFMolecules3D<QMPotential>;
        template class ChangedMolecule<QMPotential::Molecule>;
    }
}

/////////
///////// Implementation of QMComponent
/////////

static const RegisterMetaType<QMComponent> r_qmcomp(NO_ROOT);

/** Serialise to a binary datastream */
QDataStream SQUIRE_EXPORT &operator<<(QDataStream &ds, const QMComponent &qmcomp)
{
    writeHeader(ds, r_qmcomp, 1);
    ds << static_cast<const FFComponent&>(qmcomp);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SQUIRE_EXPORT &operator>>(QDataStream &ds, QMComponent &qmcomp)
{
    VersionID v = readHeader(ds, r_qmcomp);
    
    if (v == 1)
    {
        ds >> static_cast<FFComponent&>(qmcomp);
    }
    else
        throw version_error(v, "1", r_qmcomp, CODELOC);
        
    return ds;
}

/** Constructor for the forcefield called 'ffname' */
QMComponent::QMComponent(const FFName &ffname) 
            : FFComponent(ffname, QLatin1String("QM"))
{}

/** Construct from a symbol

    \throw SireError::incompatible_error
*/
QMComponent::QMComponent(const SireCAS::Symbol &symbol)
            : FFComponent(symbol, QLatin1String("QM"))
{}

/** Copy constructor */
QMComponent::QMComponent(const QMComponent &other)
            : FFComponent(other)
{}

/** Destructor */
QMComponent::~QMComponent()
{}

/** Set the QM component of the energy in the forcefield 'ff'
    to equal to the passed QMEnergy */
void QMComponent::setEnergy(FF &ff, const QMEnergy &qmnrg) const
{
    FFComponent::setEnergy(ff, this->total(), qmnrg);
}

/** Change the QM component of the energy in the forcefield 'ff'
    by 'delta' */
void QMComponent::changeEnergy(FF &ff, const QMEnergy &delta) const
{
    FFComponent::changeEnergy(ff, this->total(), delta);
}

/////////
///////// Implementation of ElementParameterName
/////////

/** This is the default name of the property containing
    the element type of an atom */
QString ElementParameterName::element_param("element");

/////////
///////// Implementation of QMPotential
/////////

static const RegisterMetaType<QMPotential> r_qmpot( MAGIC_ONLY, NO_ROOT,
                                                    "Squire::QMPotential" );
                                                    
/** Serialise to a binary datastream */
QDataStream SQUIRE_EXPORT &operator<<(QDataStream &ds, const QMPotential &qmpot)
{
    writeHeader(ds, r_qmpot, 1);
    
    SharedDataStream sds(ds);
    
    sds << qmpot.props;
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SQUIRE_EXPORT &operator>>(QDataStream &ds, QMPotential &qmpot)
{
    VersionID v = readHeader(ds, r_qmpot);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> qmpot.props;
        
        qmpot.spce = qmpot.props.property("space");
        qmpot.qmprog = qmpot.props.property("quantum program");
        qmpot.zero_energy = qmpot.props.property("zero energy")
                                       .asADouble();
    }
    else
        throw version_error(v, "1", r_qmpot, CODELOC);
        
    return ds;
}

/** Constructor */
QMPotential::QMPotential() : zero_energy(0)
{
    props.setProperty( "space", spce );
    props.setProperty( "quantum program", qmprog );
    props.setProperty( "zero energy", wrap(zero_energy) );
}

/** Copy constructor */
QMPotential::QMPotential(const QMPotential &other)
            : props(other.props), spce(other.spce), 
              qmprog(other.qmprog), zero_energy(other.zero_energy)
{}

/** Destructor */
QMPotential::~QMPotential()
{}

/** Copy assignment operator */
QMPotential& QMPotential::operator=(const QMPotential &other)
{
    if (this != &other)
    {
        props = other.props;
        spce = other.spce;
        qmprog = other.qmprog;
        zero_energy = other.zero_energy;
    }

    return *this;
}

/** Return the packed array that contains the elements for just the 
    selected atoms of 'molecule' - the elements are contained in the
    property called 'element_property' 
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
static PackedArray2D<Element> getQMParameters(const PartialMolecule &molecule,
                                              const PropertyName &element_property)
{
    const AtomElements &elements = molecule.property(element_property)
                                           .asA<AtomElements>();
    
    const AtomSelection &selected_atoms = molecule.selection();
    
    if (selected_atoms.selectedNone())
        return PackedArray2D<Element>();
    
    //create space for the parameters - only need space for CutGroups
    //that contain at least one selected atom
    QVector< QVector<Element> > elementparams( selected_atoms.nSelectedCutGroups() );
    QVector<Element>* elementparams_array = elementparams.data();

    if (selected_atoms.selectedAllCutGroups())
    {
        const int ncg = molecule.data().info().nCutGroups();
    
        for (CGIdx i(0); i<ncg; ++i)
        {
            const int nats = molecule.data().info().nAtoms(i);
            
            QVector<Element> group_elements(nats);
            Element *group_elements_array = group_elements.data();
            
            //get the arrays containing the element parameters
            //for this CutGroup
            const Element *group_e = elements.constData(i);
            
            if (selected_atoms.selectedAll(i))
            {
                for (Index j(0); j<nats; ++j)
                {
                    group_elements_array[j] = group_e[j];
                }
            }
            else
            {
                foreach (Index j, selected_atoms.selectedAtoms(i))
                {
                    group_elements_array[j] = group_e[j];
                }
            }
            
            elementparams_array[i] = group_elements;
        }
    }
    else
    {
        foreach (CGIdx i, selected_atoms.selectedCutGroups())
        {
            const int nats = molecule.data().info().nAtoms(i);
            
            QVector<Element> group_elements(nats);
            Element *group_elements_array = group_elements.data();
            
            //get the arrays containing the element parameters
            //for this CutGroup
            const Element *group_e = elements.constData(i);
            
            if (selected_atoms.selectedAll(i))
            {
                for (Index j(0); j<nats; ++j)
                {
                    group_elements_array[j] = group_e[j];
                }
            }
            else
            {
                foreach (Index j, selected_atoms.selectedAtoms(i))
                {
                    group_elements_array[j] = group_e[j];
                }
            }
            
            elementparams_array[i] = group_elements;
        }
    }
    
    return PackedArray2D<Element>( elementparams );
}

/** Return the parameters for the molecule 'mol', using the property names
    contained in 'map' to find the properties that locate the parameters 
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
QMPotential::Parameters QMPotential::getParameters(const PartialMolecule &molecule,
                                                   const PropertyMap &map)
{
    return Parameters( molecule, map[parameters().coordinates()],
                       getQMParameters(molecule, map[parameters().element()]) );
}

/** Update the parameters for the molecule going from 'old_molecule' to 
    'new_molecule', with the parameters found using the property map 'map'
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
QMPotential::Parameters 
QMPotential::updateParameters(const QMPotential::Parameters &old_params,
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
    const PropertyName &element_property = map[parameters().element()];
    
    //get what has changed
    bool new_coords = old_molecule.version(coords_property) !=
                         new_molecule.version(coords_property);
                             
    bool new_elements = old_molecule.version(element_property) !=
                          new_molecule.version(element_property);

    if (new_coords)
    {
        new_params.setAtomicCoordinates( AtomicCoords3D(new_molecule, 
                                                        coords_property) );
    }

    if (new_elements)
    {
        new_params.setAtomicParameters( getQMParameters(new_molecule,
                                                        element_property) );
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
QMPotential::Parameters
QMPotential::updateParameters(const QMPotential::Parameters &old_params,
                              const PartialMolecule &old_molecule,
                              const PartialMolecule &new_molecule,
                              const PropertyMap &old_map, const PropertyMap &new_map)
{
    if (old_molecule.selection() != new_molecule.selection())
        //the selection has changed - just get completely new parameters
        return this->getParameters(new_molecule, new_map);

    Parameters new_params = old_params;

    //get the property names
    const PropertyName &old_coords = old_map[parameters().coordinates()];
    const PropertyName &old_elements = old_map[parameters().element()];
    
    const PropertyName &new_coords = new_map[parameters().coordinates()];
    const PropertyName &new_elements = new_map[parameters().element()];
    
    //get what has changed
    bool changed_coords = (new_coords != old_coords) or
                           old_molecule.version(old_coords) !=
                           new_molecule.version(old_coords);
                             
    bool changed_elements = (new_elements != old_elements) or
                             old_molecule.version(old_elements) !=
                             new_molecule.version(old_elements);

    if (changed_coords)
        new_params.setAtomicCoordinates( AtomicCoords3D(new_molecule, 
                                                        new_coords) );

    if (changed_elements)
        new_params.setAtomicParameters( getQMParameters(new_molecule,
                                                        new_elements) );

    return new_params;
}

/** Return the QMPotential::Molecule representation of 'molecule',
    using the supplied PropertyMap to find the properties that contain
    the necessary forcefield parameters
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
QMPotential::Molecule QMPotential::parameterise(const PartialMolecule &molecule,
                                                const PropertyMap &map)
{
    return QMPotential::Molecule(molecule, *this, map);
}

/** Convert the passed group of molecules into QMPotential::Molecules,
    using the supplied PropertyMap to find the properties that contain
    the necessary forcefield parameters in each molecule
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
QMPotential::Molecules QMPotential::parameterise(const MoleculeGroup &molecules,
                                                 const PropertyMap &map)
{
    return QMPotential::Molecules(molecules, *this, map);
}

/** Return the space within which the molecules in this potential exist */
const Space& QMPotential::space() const
{
    return spce.read();
}

/** Return the handle to the quantum chemical program that is used 
    by this potential to calculate the QM energies and forces */
const QMProgram& QMPotential::quantumProgram() const
{
    return qmprog.read();
}

/** Return the absolute value of the energy which is considered
    as zero (on the relative energy scale used by this potential).
    A relative scale is used so that the QM energy can be shifted
    so that it is comparable to an MM energy */
MolarEnergy QMPotential::zeroEnergy() const
{
    return MolarEnergy(zero_energy);
}

/** Set the space within which all of the molecules in this potential
    will exist. This returns whether or not this changes the
    potential. */
bool QMPotential::setSpace(const Space &space)
{
    if ( not spce->equals(space) )
    {
        spce = space;
        props.setProperty("space", spce);
        this->changedPotential();
        return true;
    }
    else
        return false;
}

/** Set the handle to the quantum chemical program that will be
    used by this potential to calculate the QM energies and forces.
    This returns whether or not this changes this potential */
bool QMPotential::setQuantumProgram(const QMProgram &program)
{
    if ( not qmprog->equals(program) )
    {
        qmprog = program;
        props.setProperty("quantum program", qmprog);
        this->changedPotential();
        return true;
    }
    else
        return false;
}

/** Set the absolute value of the energy which is considered
    as zero (on the relative energy scale used by this potential).
    A relative scale is used so that the QM energy can be shifted
    so that it is comparable to an MM energy */
bool QMPotential::setZeroEnergy(MolarEnergy zeroenergy)
{
    if ( zero_energy != zeroenergy )
    {
        zero_energy = zeroenergy;
        props.setProperty("zero energy", wrap(zero_energy));
        this->changedPotential();
        return true;
    }
    else
        return false;
}

/** Set the property 'name' to the value 'value'. This returns whether
    or not this changes the potential

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
bool QMPotential::setProperty(const QString &name, const Property &value)
{
    if (name == QLatin1String("space"))
    {
        return this->setSpace( value.asA<Space>() );
    }
    else if (name == QLatin1String("quantum program"))
    {
        return this->setQuantumProgram( value.asA<QMProgram>() );
    }
    else if (name == QLatin1String("zero energy"))
    {
        return this->setZeroEnergy( MolarEnergy(value.asADouble()) );
    }
    else
        throw SireBase::missing_property( QObject::tr(
            "The QM potential does not have a property called \"%1\" that "
            "can be changed. Available properties are [ %2 ].")
                .arg(name, QStringList(props.propertyKeys()).join(", ")), CODELOC );
                
    return false;
}

/** Return the property of this potential called 'name'

    \throw SireBase::missing_property
*/
const Property& QMPotential::property(const QString &name) const
{
    return props.property(name);
}

/** Return whether or not this potential contains a property called 'name' */
bool QMPotential::containsProperty(const QString &name) const
{
    return props.hasProperty(name);
}

/** Return all of the properties of this potential */
const Properties& QMPotential::properties() const
{
    return props;
}

/** Map the molecules in 'molecules' into the space for this potential */
QMPotential::Molecules 
QMPotential::mapIntoSpace(const QMPotential::Molecules &mols) const
{
    #warning Need to write QMPotential::mapIntoSpace()
    return mols;
}

/** Calculate the QM forces on the molecules in 'molecules' and add them 
    onto the forces in the passed force table (scaled by the optional
    scaling constant 'scale_force') */
void QMPotential::calculateForce(const QMPotential::Molecules &molecules, 
                                 ForceTable &forcetable, 
                                 double scale_force) const
{
    if (scale_force == 0)
        return;

    //map all of the molecules so that they are in this space
    QMPotential::Molecules mapped_qmmols = QMPotential::mapIntoSpace(molecules);
    
    quantumProgram().calculateForce(mapped_qmmols, forcetable, scale_force);
}
                    
/** Calculate the QM forces on the molecules in 'molecules' and add them 
    onto the forces in the passed force table (scaled by the optional
    scaling constant 'scale_force') */
void QMPotential::calculateForce(const QMPotential::Molecules &molecules, 
                                 ForceTable &forcetable,
                                 const Symbol &symbol,
                                 const QMPotential::Components &components,
                                 double scale_force) const
{
    if (symbol == components.total())
        this->calculateForce(molecules, forcetable, scale_force);
        
    else
        throw SireFF::missing_component( QObject::tr(
            "There is no force component in potential %1 - available "
            "components are %2.")
                .arg(this->what())
                .arg(components.total().toString()), CODELOC );
}

/** Calculate the QM energy of the collection of molecules in 'molecules'
    and add it onto 'nrg' (scaled by the optional scaling constant 'scale_energy') */
void QMPotential::calculateEnergy(const QMPotential::Molecules &molecules, 
                                  QMPotential::Energy &nrg, double scale_energy) const
{
    if (scale_energy == 0)
        return;

    //map all of the molecules so that they are in this space
    Molecules mapped_mols = this->mapIntoSpace(molecules);

    //now tell the qm_program to calculate the energy that we need
    double qmnrg = qmprog->calculateEnergy(mapped_mols);
    
    nrg += Energy( scale_energy * (qmnrg - zero_energy) );
}

/** Calculate the fields on the atoms */
void QMPotential::calculateField(const QMPotential::Molecules &molecules,
                                 FieldTable &fieldtable, const SireFF::Probe &probe,
                                 double scale_field) const
{
    if (scale_field == 0)
        return;

    //map all of the molecules so that they are in this space
    QMPotential::Molecules mapped_qmmols = QMPotential::mapIntoSpace(molecules);
    
    quantumProgram().calculateField(mapped_qmmols, fieldtable, probe, scale_field);
}

/** Calculate the fields on the atoms */
void QMPotential::calculateField(const QMPotential::Molecules &molecules,
                                 FieldTable &fieldtable, const SireFF::Probe &probe,
                                 const Symbol &symbol, 
                                 const QMPotential::Components &components,
                                 double scale_field) const
{
    if (symbol == components.total())
        this->calculateField(molecules, fieldtable, probe, scale_field);
        
    else
        throw SireFF::missing_component( QObject::tr(
            "There is no field component in potential %1 - available "
            "components are %2.")
                .arg(this->what())
                .arg(components.total().toString()), CODELOC );
}

/** Calculate the potentials on the atoms */
void QMPotential::calculatePotential(const QMPotential::Molecules &molecules,
                                     PotentialTable &pottable, 
                                     const SireFF::Probe &probe,
                                     double scale_potential) const
{
    if (scale_potential == 0)
        return;

    //map all of the molecules so that they are in this space
    QMPotential::Molecules mapped_qmmols = QMPotential::mapIntoSpace(molecules);
    
    quantumProgram().calculatePotential(mapped_qmmols, pottable, probe, scale_potential);
}

/** Calculate the potentials on the atoms */
void QMPotential::calculatePotential(const QMPotential::Molecules &molecules,
                                     PotentialTable &pottable, 
                                     const SireFF::Probe &probe,
                                     const Symbol &symbol, 
                                     const QMPotential::Components &components,
                                     double scale_potential) const
{
    if (symbol == components.total())
        this->calculatePotential(molecules, pottable, probe, scale_potential);
        
    else
        throw SireFF::missing_component( QObject::tr(
            "There is no potential component in potential %1 - available "
            "components are %2.")
                .arg(this->what())
                .arg(components.total().toString()), CODELOC );
}

/** Return the command file that would be used to calculate the forces on
    the molecules in 'molecules' */
QString QMPotential::forceCommandFile(const QMPotential::Molecules &molecules,
                                      const ForceTable &forcetable) const
{
    //map all of the molecules so that they are in this space
    Molecules mapped_mols = this->mapIntoSpace(molecules);

    return qmprog->forceCommandFile(mapped_mols, forcetable);
}

/** Return the command file that would be used to calculate the energy of
    the molecules in 'molecules' */
QString QMPotential::energyCommandFile(const QMPotential::Molecules &molecules) const
{
    //map all of the molecules so that they are in this space
    Molecules mapped_mols = this->mapIntoSpace(molecules);

    return qmprog->energyCommandFile(mapped_mols);
}
    
/** Return the command file that would be used to calculate the field
    for the passed molecules, if they are in the passed fieldtable */
QString QMPotential::fieldCommandFile(const Molecules &molecules,
                                      const FieldTable &fieldtable,
                                      const SireFF::Probe &probe) const
{
    Molecules mapped_mols = this->mapIntoSpace(molecules);
    
    return qmprog->fieldCommandFile(mapped_mols, fieldtable, probe);
}

/** Return the command file that would be used to calculate the potential
    for the passed molecules, if they are in the passed potential table */
QString QMPotential::potentialCommandFile(const Molecules &molecules,
                                          const PotentialTable &pottable,
                                          const SireFF::Probe &probe) const
{
    Molecules mapped_mols = this->mapIntoSpace(molecules);
    
    return qmprog->potentialCommandFile(mapped_mols, pottable, probe);
}
