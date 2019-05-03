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

#include "integratorworkspace.h"
#include "velocitygenerator.h"

#include "SireSystem/system.h"

#include "SireMol/moleculeview.h"
#include "SireMol/molecule.h"
#include "SireMol/partialmolecule.h"
#include "SireMol/moleditor.h"
#include "SireMol/atomcoords.h"
#include "SireMol/atomelements.h"
#include "SireMol/atomvelocities.h"
#include "SireMol/atommasses.h"
#include "SireMol/atomcoords.h"
#include "SireMol/molidx.h"

#include "SireBase/quickcopy.hpp"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireMol/errors.h"

#include <QDebug>

using namespace SireMove;
using namespace SireSystem;
using namespace SireFF;
using namespace SireMol;
using namespace SireCAS;
using namespace SireBase;
using namespace SireStream;
using namespace SireUnits::Dimension;

//////////
////////// Implementation of IntegratorWorkspace
//////////

static const RegisterMetaType<IntegratorWorkspace> r_intws(MAGIC_ONLY,
                                                    "SireMove::IntegratorWorkspace");

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                        const IntegratorWorkspace &intws)
{
    writeHeader(ds, r_intws, 1);
    
    SharedDataStream sds(ds);
    
    if (intws.need_new_forces)
    {
        sds << true << intws.sys << intws.molgroup << intws.map
            << static_cast<const Property&>(intws);
    }
    else
    {
        sds << false << intws.sys << intws.molgroup << intws.map
            << intws.molforces << intws.last_nrg_component
            << static_cast<const Property&>(intws);
    }
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, IntegratorWorkspace &intws)
{
    VersionID v = readHeader(ds, r_intws);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);

        bool need_new_forces;
        
        sds >> need_new_forces;
        
        if (need_new_forces)
        {
            intws.need_new_forces = true;
        
            sds >> intws.sys >> intws.molgroup >> intws.map
                >> static_cast<Property&>(intws);
        }
        else
        {
            intws.need_new_forces = false;
            
            sds >> intws.sys >> intws.molgroup >> intws.map
                >> intws.molforces >> intws.last_nrg_component
                >> static_cast<Property&>(intws);
        }
    }
    else
        throw version_error( v, "1", r_intws, CODELOC );

    return ds;
}

/** Constructor */
IntegratorWorkspace::IntegratorWorkspace(const PropertyMap &m) 
                    : Property(), map(m), need_new_forces(true)
{}

/** Construct to hold the variables used to integrate the molecules in 'molgroup' */
IntegratorWorkspace::IntegratorWorkspace(const MoleculeGroup &molecule_group,
                                         const PropertyMap &m)
                    : Property(), molgroup(molecule_group), 
                      molforces(molecule_group), map(m), need_new_forces(true)
{}

/** Copy constructor */
IntegratorWorkspace::IntegratorWorkspace(const IntegratorWorkspace &other)
                    : Property(other), 
                      sys(other.sys),
                      molgroup(other.molgroup),
                      molforces(other.molforces),
                      last_nrg_component(other.last_nrg_component),
                      map(other.map),
                      need_new_forces(other.need_new_forces)
{}

/** Destructor */
IntegratorWorkspace::~IntegratorWorkspace()
{}

/** Copy assignment operator */
IntegratorWorkspace& IntegratorWorkspace::operator=(const IntegratorWorkspace &other)
{
    if (this != &other)
    {
        Property::operator=(other);
        sys = other.sys;
        molgroup = other.molgroup;
        molforces = other.molforces;
        last_nrg_component = other.last_nrg_component;
        map = other.map;
        need_new_forces = other.need_new_forces;
    }
    
    return *this;
}

/** Comparison operator */
bool IntegratorWorkspace::operator==(const IntegratorWorkspace &other) const
{
    return (this == &other) or
           (sys == other.sys and
            molgroup == other.molgroup and 
            last_nrg_component == other.last_nrg_component and
            need_new_forces == other.need_new_forces and 
            map == other.map and
            molforces == other.molforces);
}

/** Comparison operator */
bool IntegratorWorkspace::operator!=(const IntegratorWorkspace &other) const
{
    return not this->operator==(other);
}

/** Return the molecule group whose molecules will be moved by this
    integrator */
const MoleculeGroup& IntegratorWorkspace::moleculeGroup() const
{
    return molgroup.read();
}

/** Return the force table containing the forces on the molecule */
const ForceTable& IntegratorWorkspace::forceTable() const
{
    return molforces;
}

/** Set the system to be integrated - this updates the molecules in 
    the passed molecule group - this returns whether or not the
    system has changed */
bool IntegratorWorkspace::setSystem(const System &system)
{
    if (system.subVersion() != 0)
        throw SireError::incompatible_error( QObject::tr(
                    "You can not give an integrator workspace (%1) a System "
                    "which is in a subversion state (%2).")
                        .arg(this->toString()).arg(system.toString()), CODELOC );

    if (sys.UID() == system.UID() and sys.version() == system.version())
        //nothing needs to change
        return false;
    
    if (system.contains(molgroup.read().number()))
    {
        MolGroupPtr new_molgroup = system[molgroup.read().number()];
        
        if (new_molgroup.read().version().majorVersion() 
                    != molgroup.read().version().majorVersion())
        {
            //we need to reallocate new space for the forces
            molforces = ForceTable(new_molgroup.read());
        }
        
        molgroup = new_molgroup;
    }
    else
    {
        molgroup.edit().update(system.molecules());
    }
                
    need_new_forces = true;
    last_nrg_component = Symbol();
    sys = system;
    
    return true;
}

/** Return the system being integrated */
const System& IntegratorWorkspace::system() const
{
    return sys;
}

/** Return the system being integrated */
System& IntegratorWorkspace::nonConstsystem() 
{
    return sys;
}

/** Return the property map used to find the properties that are
    required for integration */
const PropertyMap& IntegratorWorkspace::propertyMap() const
{
    return map;
}

/** Set the property map that is used to find the properties
    that are used for integration */
void IntegratorWorkspace::setPropertyMap(const PropertyMap &m)
{
    map = m;
}

/** Set the random number generator that is used during integration */
void IntegratorWorkspace::setGenerator(const RanGenerator&)
{}

/** Function called when a property is changed */
void IntegratorWorkspace::changedProperty(const QString&)
{}

/** Set the property used to find the coordinates of the molecules */
void IntegratorWorkspace::setCoordinatesProperty(const PropertyName &source)
{
    if (map["coordinates"] != source)
    {
        map.set("coordinates", source);
        this->changedProperty("coordinates");
    }
}

/** Set the property used to find the system space */
void IntegratorWorkspace::setSpaceProperty(const PropertyName &source)
{
    if (map["space"] != source)
    {
        map.set("space", source);
        this->changedProperty("space");
    }
}

/** Set the property used to find the velocities of the molecules */
void IntegratorWorkspace::setVelocitiesProperty(const PropertyName &source)
{
    if (map["velocity"] != source)
    {
        map.set("velocity", source);
        this->changedProperty("velocity");
    }
}

/** Set the property used to find the masses of the molecules */
void IntegratorWorkspace::setMassesProperty(const PropertyName &source)
{
    if (map["mass"] != source)
    {
        map.set("mass", source);
        this->changedProperty("mass");
    }
}

/** Set the property used to find the elements of the atoms in the molecule */
void IntegratorWorkspace::setElementsProperty(const PropertyName &source)
{
    if (map["element"] != source)
    {
        map.set("element", source);
        this->changedProperty("element");
    }
}

/** Set the property used to generate new velocities */
void IntegratorWorkspace::setVelocityGeneratorProperty(const PropertyName &source)
{
    if (map["velocity generator"] != source)
    {
        map.set("velocity generator", source);
        this->changedProperty("velocity generator");
    }
}

/** Return the property that contains the molecule coordinates */
PropertyName IntegratorWorkspace::coordinatesProperty() const
{
    return map["coordinates"];
}

/** Return the property that contains the system space */
PropertyName IntegratorWorkspace::spaceProperty() const
{
    return map["space"];
}

/** Return the property that contains the molecule velocities */
PropertyName IntegratorWorkspace::velocitiesProperty() const
{
    return map["velocity"];
}

/** Return the property that contains the molecule masses */
PropertyName IntegratorWorkspace::massesProperty() const
{
    return map["mass"];
}

/** Return the property that contains the molecule elements */
PropertyName IntegratorWorkspace::elementsProperty() const
{
    return map["element"];
}

/** Return the property used to generate missing velocities */
PropertyName IntegratorWorkspace::velocityGeneratorProperty() const
{
    return map["velocity generator"];
}

/** Calculate the current forces on the molecules in the molecule
    group using the energy component 'nrg_component' */
bool IntegratorWorkspace::calculateForces(const Symbol &nrg_component)
{
    if (need_new_forces or last_nrg_component != nrg_component)
    {
        molforces.initialiseTables();
        sys.force(molforces, nrg_component);
        last_nrg_component = nrg_component;
        need_new_forces = false;
        
        return true;
    }
    else
        return false;
}

/** Return whether or not the forces need calculating for the energy
    component 'nrg_component' */
bool IntegratorWorkspace::forcesNeedCalculating(const Symbol &nrg_component) const
{
    return need_new_forces or nrg_component != last_nrg_component;
}

/** Set it so that the forces must now be recalculated from scratch */
void IntegratorWorkspace::mustNowRecalculateFromScratch()
{
    need_new_forces = true;
    last_nrg_component = Symbol();
}

/** Internal function used to update the system and molecule group with
    changed molecules */
void IntegratorWorkspace::pvt_update(const Molecules &changed_mols)
{
    sys.update(changed_mols);
    
    if (sys.contains(molgroup.read().number()))
    {
        molgroup = sys[molgroup.read().number()];
    }
    else
    {
        molgroup.edit().update(changed_mols);
    }
    
    need_new_forces = true;
    last_nrg_component = Symbol();
}

/** Tell the contained system to collect statistics */
void IntegratorWorkspace::collectStatistics()
{
    sys.collectStats();
}

Q_GLOBAL_STATIC( NullIntegratorWorkspace, nullIntegratorWorkspace );

/** Return the global null workspace */
const NullIntegratorWorkspace& IntegratorWorkspace::null()
{
    return *(nullIntegratorWorkspace());
}

//////////
////////// Implementation of NullIntegratorWorkspace
//////////

static const RegisterMetaType<NullIntegratorWorkspace> r_nullintws;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, 
                                        const NullIntegratorWorkspace &nullintws)
{
    writeHeader(ds, r_nullintws, 1);
    
    ds << static_cast<const IntegratorWorkspace&>(nullintws);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, 
                                        NullIntegratorWorkspace &nullintws)
{
    VersionID v = readHeader(ds, r_nullintws);
    
    if (v == 1)
    {
        ds >> static_cast<IntegratorWorkspace&>(nullintws);
    }
    else
        throw version_error( v, "1", r_nullintws, CODELOC );
        
    return ds;
}

/** Constructor */
NullIntegratorWorkspace::NullIntegratorWorkspace()
                        : ConcreteProperty<NullIntegratorWorkspace,IntegratorWorkspace>()
{}

/** Copy constructor */
NullIntegratorWorkspace::NullIntegratorWorkspace(const NullIntegratorWorkspace &other)
          : ConcreteProperty<NullIntegratorWorkspace,IntegratorWorkspace>(other)
{}

/** Destructor */
NullIntegratorWorkspace::~NullIntegratorWorkspace()
{}

/** Assignment operator */
NullIntegratorWorkspace& NullIntegratorWorkspace::operator=(
                                            const NullIntegratorWorkspace &other)
{
    IntegratorWorkspace::operator=(other);
    return *this;
}

/** Comparison operator */
bool NullIntegratorWorkspace::operator==(const NullIntegratorWorkspace&) const
{
    return true;
}

/** Comparison operator */
bool NullIntegratorWorkspace::operator!=(const NullIntegratorWorkspace &) const
{
    return false;
}

/** Regenerate the velocities using the passed generator */
void NullIntegratorWorkspace::regenerateVelocities(const VelocityGenerator&)
{}

/** Zero kinetic energy */
MolarEnergy NullIntegratorWorkspace::kineticEnergy() const
{
    return MolarEnergy(0);
}

/** This contains no molecules 

    \throw SireMol::missing_molecule
*/
MolarEnergy NullIntegratorWorkspace::kineticEnergy(MolNum molnum) const
{
    throw SireMol::missing_molecule( QObject::tr(
        "The null integrator workspace does not contain any molecules, so it "
        "definitely does not contain %1.")
            .arg(molnum.toString()), CODELOC );
            
    return MolarEnergy(0);
}

/** This contains no molecules 

    \throw SireMol::missing_molecule
*/
MolarEnergy NullIntegratorWorkspace::kineticEnergy(const MoleculeView &molview) const
{
    throw SireMol::missing_molecule( QObject::tr(
        "The null integrator workspace does not contain any molecules, so it "
        "definitely does not contain %1.")
            .arg(molview.toString()), CODELOC );
            
    return MolarEnergy(0);
}

const char* NullIntegratorWorkspace::typeName()
{
    return QMetaType::typeName( qMetaTypeId<NullIntegratorWorkspace>() );
}

//////////
////////// Implementation of AtomicVelocityWorkspace
//////////

static const RegisterMetaType<AtomicVelocityWorkspace> r_atvelws;

QDataStream &operator<<(QDataStream &ds, 
                                        const AtomicVelocityWorkspace &atvelws)
{
    writeHeader(ds, r_atvelws, 1);
    
    SharedDataStream sds(ds);
    
    sds << static_cast<const IntegratorWorkspace&>(atvelws);
    
    return ds;
}

QDataStream &operator>>(QDataStream &ds,
                                        AtomicVelocityWorkspace &atvelws)
{
    VersionID v = readHeader(ds, r_atvelws);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        AtomicVelocityWorkspace ws;
        
        sds >> static_cast<IntegratorWorkspace&>(ws);
        
        ws.rebuildFromScratch();
        
        atvelws = ws;
    }
    else
        throw version_error(v, "1", r_atvelws, CODELOC);
        
    return ds;
}

static QVector<double> getMasses(const QVector<MolarMass> &masses)
{
    int sz = masses.count();
    QVector<double> atom_masses(sz);
    atom_masses.squeeze();
    
    const MolarMass *masses_array = masses.constData();
    double *atom_masses_array = atom_masses.data();
    
    for (int i=0; i<sz; ++i)
    {
        atom_masses_array[i] = masses_array[i].value();
    }
    
    return atom_masses;
}

static QVector<double> getMasses(const QVector<Element> &elements)
{
    int sz = elements.count();
    QVector<double> atom_masses(sz);
    atom_masses.squeeze();
    
    const Element *elements_array = elements.constData();
    double *atom_masses_array = atom_masses.data();
    
    for (int i=0; i<sz; ++i)
    {
        atom_masses_array[i] = elements_array[i].mass().value();
    }
    
    return atom_masses;
}

static QVector<Vector> getMomenta(const QVector<Velocity3D> &velocities,
                                  const QVector<double> &masses)
{
    int nats = velocities.count();
    
    QVector<Vector> momenta(nats);
    momenta.squeeze();
    
    const Velocity3D *vel_array = velocities.constData();
    const double *mass_array = masses.constData();
    
    Vector *mom_array = momenta.data();
    
    for (int i=0; i<nats; ++i)
    {
        mom_array[i] = vel_array[i].value() * mass_array[i];
    }
    
    return momenta;
}

static QVector<Velocity3D> getVelocities(const QVector<Vector> &momenta,
                                         const QVector<double> &masses)
{
    int nats = momenta.count();
    
    QVector<Velocity3D> velocities(nats);
    velocities.squeeze();
    
    const Vector *mom_array = momenta.constData();
    const double *mass_array = masses.constData();
    
    Velocity3D *vel_array = velocities.data();
    
    for (int i=0; i<nats; ++i)
    {
        if (mass_array[i] != 0)
            vel_array[i] = Velocity3D(mom_array[i] / mass_array[i]);
        else
            vel_array[i] = Velocity3D(0);
    }
    
    return velocities;
}

/** Internal function used to rebuild all of the arrays from the
    current properties and current system */
void AtomicVelocityWorkspace::rebuildFromScratch()
{
    const System &sys = this->system();
    
    PropertyName coords_property = this->coordinatesProperty();
    PropertyName mass_property = this->massesProperty();
    PropertyName element_property = this->elementsProperty();
    PropertyName velocity_property = this->velocitiesProperty();
    PropertyName velgen_property = this->velocityGeneratorProperty();
    
    const MoleculeGroup &molgroup = this->moleculeGroup();
    const ForceTable &forcetable = this->forceTable();
    
    int nmols = molgroup.nMolecules();
    
    atom_coords = QVector< QVector<Vector> >(nmols);
    atom_momenta = QVector< QVector<Vector> >(nmols);
    atom_masses = QVector< QVector<double> >(nmols);
    
    atom_coords.squeeze();
    atom_momenta.squeeze();
    atom_masses.squeeze();
    
    QVector<Vector> *atom_coords_array = atom_coords.data();
    QVector<Vector> *atom_mom_array = atom_momenta.data();
    QVector<double> *atom_masses_array = atom_masses.data();

    atom_forces = QVector< QVector<Vector> >();
    QVector<Vector> *atom_forces_array = 0;
    
    if (sys.containsProperty(velgen_property))
        vel_generator = sys.property(velgen_property).asA<VelocityGenerator>();
    else
        vel_generator = NullVelocityGenerator();
    
    for (int i=0; i<nmols; ++i)
    {
        MolNum molnum = molgroup.molNumAt(i);
        const ViewsOfMol &mol = molgroup[molnum].data();
        
        const MoleculeData &moldata = mol.data();
        
        if (mol.selectedAll())
        {
            atom_coords_array[i] = moldata.property(coords_property)
                                          .asA<AtomCoords>().toVector();

            if (moldata.hasProperty(mass_property))
            {
                atom_masses_array[i] = ::getMasses(
                                            moldata.property(mass_property)
                                                   .asA<AtomMasses>().toVector());
            }
            else
            {
                atom_masses_array[i] = ::getMasses(
                                            moldata.property(element_property)
                                                   .asA<AtomElements>().toVector());
            }
                                          
            if (moldata.hasProperty(velocity_property))
            {
                atom_mom_array[i] = ::getMomenta(
                                        moldata.property(velocity_property)
                                               .asA<AtomVelocities>().toVector(),
                                        atom_masses_array[i]);
            }
            else
            {
                atom_mom_array[i] = ::getMomenta(
                                        vel_generator.read().generate(mol, propertyMap())
                                                     .toVector(),
                                        atom_masses_array[i]);
            }
        }
        else
        {
            AtomSelection selected_atoms = mol.selection();
            
            atom_coords_array[i] = moldata.property(coords_property)
                                          .asA<AtomCoords>().toVector(selected_atoms);
            
            if (moldata.hasProperty(mass_property))
            {
                atom_masses_array[i] = ::getMasses(
                                            moldata.property(mass_property)
                                                   .asA<AtomMasses>()
                                                   .toVector(selected_atoms));
            }
            else
            {
                atom_masses_array[i] = ::getMasses(
                                            moldata.property(element_property)
                                                   .asA<AtomElements>()
                                                   .toVector(selected_atoms));
            }
                                          
            if (moldata.hasProperty(velocity_property))
            {
                atom_mom_array[i] = ::getMomenta(
                                        moldata.property(velocity_property)
                                         .asA<AtomVelocities>().toVector(selected_atoms),
                                            atom_masses_array[i]);
            }
            else
            {
                atom_mom_array[i] = ::getMomenta(
                                        vel_generator.read().generate(mol, propertyMap())
                                                         .toVector(selected_atoms),
                                            atom_masses_array[i]);
                                                         
            }
            
            if (atom_forces_array == 0)
            {
                atom_forces = QVector< QVector<Vector> >(nmols);
                atom_forces.squeeze();
                atom_forces_array = atom_forces.data();
            }
            
            atom_forces_array[i] = forcetable.getTable(molnum).toVector(selected_atoms);
        }
    }
}

/** Calculate the forces caused by the passed energy component */
bool AtomicVelocityWorkspace::calculateForces(const Symbol &nrg_component)
{
    if (not IntegratorWorkspace::calculateForces(nrg_component))
        return false;
        
    else if (atom_forces.isEmpty())
    {
        //there is nothing to do, as we have no partial molecules,
        //so all of the forces can be obtained direct from the forcetable
        return true;
    }
    else
    {
        int nmols = atom_forces.count();
        QVector<Vector> *atom_forces_array = atom_forces.data();
        
        const MoleculeGroup &molgroup = moleculeGroup();
        const ForceTable &forcetable = forceTable();
        
        for (int i=0; i<nmols; ++i)
        {
            MolNum molnum = molgroup.molNumAt(i);
            const ViewsOfMol &mol = molgroup[molnum].data();
        
            if (mol.selectedAll())
                atom_forces_array[i] = QVector<Vector>();

            else
                atom_forces_array[i] = forcetable.getTable(molnum)
                                                 .toVector(mol.selection());
        }
        
        return true;
    }
}

/** Regenerate the velocities using the passed generator */
void AtomicVelocityWorkspace::regenerateVelocities(const VelocityGenerator &generator)
{
    const MoleculeGroup &molgroup = moleculeGroup();
    
    int nmols = molgroup.nMolecules();

    QVector<Vector> *atom_mom_array = atom_momenta.data();
    
    const QVector<double> *atom_masses_array = atom_masses.constData();
    
    for (int i=0; i<nmols; ++i)
    {
        MolNum molnum = molgroup.molNumAt(i);
        const ViewsOfMol &mol = molgroup[molnum].data();
        
        if (mol.selectedAll())
        {
            atom_mom_array[i] = ::getMomenta(
                                        generator.generate(mol, propertyMap()).toVector(),
                                            atom_masses_array[i]);
        }
        else
        {
            AtomSelection selected_atoms = mol.selection();

            atom_mom_array[i] = ::getMomenta(
                                        generator.generate(mol, propertyMap())
                                                     .toVector(mol.selection()),
                                            atom_masses_array[i]);
        }
    }
}

/** Construct an empty workspace */
AtomicVelocityWorkspace::AtomicVelocityWorkspace(const PropertyMap &map)
       : ConcreteProperty<AtomicVelocityWorkspace,IntegratorWorkspace>(map)
{}

/** Construct a workspace to operate on the passed molecule group */
AtomicVelocityWorkspace::AtomicVelocityWorkspace(const MoleculeGroup &molgroup,
                                                 const PropertyMap &map)
       : ConcreteProperty<AtomicVelocityWorkspace,IntegratorWorkspace>(molgroup, map)
{
    this->rebuildFromScratch();
}

/** Copy constructor */
AtomicVelocityWorkspace::AtomicVelocityWorkspace(const AtomicVelocityWorkspace &other)
       : ConcreteProperty<AtomicVelocityWorkspace,IntegratorWorkspace>(other),
         atom_coords(other.atom_coords), atom_momenta(other.atom_momenta),
         atom_forces(other.atom_forces), atom_masses(other.atom_masses),
         vel_generator(other.vel_generator)
{}

/** Destructor */
AtomicVelocityWorkspace::~AtomicVelocityWorkspace()
{}

/** Copy assignment operator */
AtomicVelocityWorkspace& 
AtomicVelocityWorkspace::operator=(const AtomicVelocityWorkspace &other)
{
    if (this != &other)
    {
        atom_coords = other.atom_coords;
        atom_momenta = other.atom_momenta;
        atom_forces = other.atom_forces;
        atom_masses = other.atom_masses;
        vel_generator = other.vel_generator;
        IntegratorWorkspace::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool AtomicVelocityWorkspace::operator==(const AtomicVelocityWorkspace &other) const
{
    return IntegratorWorkspace::operator==(other);
}

/** Comparison operator */
bool AtomicVelocityWorkspace::operator!=(const AtomicVelocityWorkspace &other) const
{
    return not AtomicVelocityWorkspace::operator==(other);
}

const char* AtomicVelocityWorkspace::typeName()
{
    return QMetaType::typeName( qMetaTypeId<AtomicVelocityWorkspace>() );
}

/** This is called whenever a property is changed */
void AtomicVelocityWorkspace::changedProperty(const QString&)
{
    this->rebuildFromScratch();
}

static double getKineticEnergy(const QVector<double> &masses,
                               const QVector<Vector> &momenta)
{
    int nats = masses.count();
    BOOST_ASSERT( momenta.count() == nats );
    
    const double *masses_array = masses.constData();
    const Vector *mom_array = momenta.constData();
    
    double nrg = 0;
    
    // Kinetic energy is p**2 / 2m
    //
    // Internal units are self-consistent 
    // (Angstrom / AKMA, and g mol-1, which give
    // energies in kcal mol-1
    
    for (int i=0; i<nats; ++i)
    {
        //qDebug() << " i " << i << " m[i] " << masses_array[i] << "mom_array[i] " << mom_array[i].toString();
        if (masses_array[i] != 0)
            nrg += mom_array[i].length2() / masses_array[i];
    }
    
    return 0.5 * nrg;
}

/** Return the total kinetic energy of the molecules in the molecule group */
MolarEnergy AtomicVelocityWorkspace::kineticEnergy() const
{
    int nmols = atom_momenta.count();
    BOOST_ASSERT( atom_masses.count() == nmols );
    
    const QVector<double> *masses_array = atom_masses.constData();
    const QVector<Vector> *mom_array = atom_momenta.constData();
    
    double nrg = 0;
    
    for (int i=0; i<nmols; ++i)
    {
        nrg += ::getKineticEnergy(masses_array[i], mom_array[i]);
    }
    
    return MolarEnergy(nrg);
}

/** Return the total kinetic energy of the molecule with number 'molnum'

    \throw SireMol::missing_molecule
*/
MolarEnergy AtomicVelocityWorkspace::kineticEnergy(MolNum molnum) const
{
    int i = this->moleculeGroup().indexOf(molnum);
    
    return MolarEnergy( ::getKineticEnergy(atom_masses[i], atom_momenta[i]) );
}

/** Return the total kinetic energy of the atoms in the the molecule viewed
    in 'molview'
    
    \throw SireMol::missing_molecule
*/
MolarEnergy AtomicVelocityWorkspace::kineticEnergy(const MoleculeView &molview) const
{
    throw SireError::incomplete_code( QObject::tr("This needs writing!"), CODELOC );
    return MolarEnergy(0);
}

/** Return the number of molecules that are being integrated */
int AtomicVelocityWorkspace::nMolecules() const
{
    return atom_masses.count();
}

/** Return the number of atoms of the ith molecule. This does not
    check that 'i' is a valid index - use of an invalid index
    will lead to undefined results (e.g. crash or worse) */
int AtomicVelocityWorkspace::nAtoms(int i) const
{
    return atom_masses.constData()[i].count();
}

/** Return the array of the coordinates of the ith molecule. This does not
    check that 'i' is a valid index - use of an invalid index
    will lead to undefined results (e.g. crash or worse) */
Vector* AtomicVelocityWorkspace::coordsArray(int i)
{
    return atom_coords.data()[i].data();
}

/** Return the array of momenta of the ith molecule. This does not
    check that 'i' is a valid index - use of an invalid index
    will lead to undefined results (e.g. crash or worse) */
Vector* AtomicVelocityWorkspace::momentaArray(int i)
{
    return atom_momenta.data()[i].data();
}

/** Return the array of coordinates of the ith molecule. This does not
    check that 'i' is a valid index - use of an invalid index
    will lead to undefined results (e.g. crash or worse) */
const Vector* AtomicVelocityWorkspace::coordsArray(int i) const
{
    return atom_coords.constData()[i].constData();
}

/** Return the array of forces of the ith molecule. This does not
    check that 'i' is a valid index - use of an invalid index
    will lead to undefined results (e.g. crash or worse) */
const Vector* AtomicVelocityWorkspace::forceArray(int i) const
{
    if (atom_forces.isEmpty())
        return forceTable().getTable( moleculeGroup().molNumAt(i) ).constValueData();

    const QVector<Vector> &forces = atom_forces.constData()[i];
    
    if (forces.isEmpty())
        //we can get the forces straight from the forcetable
        return forceTable().getTable( moleculeGroup().molNumAt(i) ).constValueData();

    else
        return forces.constData();
}

/** Return the array of momenta of the ith molecule. This does not
    check that 'i' is a valid index - use of an invalid index
    will lead to undefined results (e.g. crash or worse) */
const Vector* AtomicVelocityWorkspace::momentaArray(int i) const
{
    return atom_momenta.constData()[i].constData();
}

/** Return the array of masses of the ith molecule. This does not
    check that 'i' is a valid index - use of an invalid index
    will lead to undefined results (e.g. crash or worse) */
const double* AtomicVelocityWorkspace::massArray(int i) const
{
    return atom_masses.constData()[i].constData();
}

/** Return the array of coordinates of the ith molecule. This does not
    check that 'i' is a valid index - use of an invalid index
    will lead to undefined results (e.g. crash or worse) */
const Vector* AtomicVelocityWorkspace::constCoordsArray(int i) const
{
    return AtomicVelocityWorkspace::coordsArray(i);
}

/** Return the array of forces of the ith molecule. This does not
    check that 'i' is a valid index - use of an invalid index
    will lead to undefined results (e.g. crash or worse) */
const Vector* AtomicVelocityWorkspace::constForceArray(int i) const
{
    return AtomicVelocityWorkspace::forceArray(i);
}

/** Return the array of momenta of the ith molecule. This does not
    check that 'i' is a valid index - use of an invalid index
    will lead to undefined results (e.g. crash or worse) */
const Vector* AtomicVelocityWorkspace::constMomentaArray(int i) const
{
    return AtomicVelocityWorkspace::momentaArray(i);
}

/** Return the array of masses of the ith molecule. This does not
    check that 'i' is a valid index - use of an invalid index
    will lead to undefined results (e.g. crash or worse) */
const double* AtomicVelocityWorkspace::constMassArray(int i) const
{
    return AtomicVelocityWorkspace::massArray(i);
}

/** Set the system that is being integrated */
bool AtomicVelocityWorkspace::setSystem(const System &new_system)
{
    if (IntegratorWorkspace::setSystem(new_system))
    {
        this->rebuildFromScratch();
        return true;
    }
    else
        return false;
}

/** Save the coordinates back to the system */
void AtomicVelocityWorkspace::commitCoordinates()
{
    int nmols = atom_coords.count();
    
    const MoleculeGroup &molgroup = moleculeGroup();
    const Molecules &molecules = molgroup.molecules();
    
    BOOST_ASSERT( molgroup.nMolecules() == nmols );
    
    const QVector<Vector> *coords_array = atom_coords.constData();
    
    PropertyName coords_property = coordinatesProperty();
    
    Molecules changed_mols;
    changed_mols.reserve(nmols);
    
    for (int i=0; i<nmols; ++i)
    {
        MolNum molnum = molgroup.molNumAt(i);
        
        const ViewsOfMol &mol = molecules[molnum];
        
        AtomCoords coords = mol.data().property(coords_property)
                                      .asA<AtomCoords>();
                                          
        if (mol.selectedAll())
            coords.copyFrom(coords_array[i]);
        else
            coords.copyFrom(coords_array[i], mol.selection());

        changed_mols.add( mol.molecule().edit()
                             .setProperty(coords_property, coords)
                             .commit() );
    }
    
    IntegratorWorkspace::pvt_update(changed_mols);
}

/** Save the velocities back to the system */
void AtomicVelocityWorkspace::commitVelocities()
{
    int nmols = atom_coords.count();
    
    const MoleculeGroup &molgroup = moleculeGroup();
    const Molecules &molecules = molgroup.molecules();
    
    BOOST_ASSERT( molgroup.nMolecules() == nmols );
    
    const QVector<Vector> *mom_array = atom_momenta.constData();
    const QVector<double> *mass_array = atom_masses.constData();
    
    PropertyName vels_property = velocitiesProperty();
    
    Molecules changed_mols;
    changed_mols.reserve(nmols);
    
    for (int i=0; i<nmols; ++i)
    {
        MolNum molnum = molgroup.molNumAt(i);
        
        const ViewsOfMol &mol = molecules[molnum];
        
        AtomVelocities vels;
        
        if (mol.data().hasProperty(vels_property))
        {
            vels = mol.data().property(vels_property)
                             .asA<AtomVelocities>();
        }
        else
        {
            vels = AtomVelocities(mol.data().info());
        }
                                          
        if (mol.selectedAll())
            vels.copyFrom( ::getVelocities(mom_array[i], mass_array[i]) );
        else
            vels.copyFrom( ::getVelocities(mom_array[i], mass_array[i]), 
                           mol.selection() );

        changed_mols.add( mol.molecule().edit()
                             .setProperty(vels_property, vels)
                             .commit() );
    }
    
    IntegratorWorkspace::pvt_update(changed_mols);
}

/** Save both the coordinates and velocities back to the system */
void AtomicVelocityWorkspace::commitCoordinatesAndVelocities()
{
    int nmols = atom_coords.count();
    
    const MoleculeGroup &molgroup = moleculeGroup();
    const Molecules &molecules = molgroup.molecules();
    
    BOOST_ASSERT( molgroup.nMolecules() == nmols );
    
    const QVector<Vector> *coords_array = atom_coords.constData();
    const QVector<Vector> *mom_array = atom_momenta.constData();
    const QVector<double> *mass_array = atom_masses.constData();
    
    PropertyName coords_property = coordinatesProperty();
    PropertyName vels_property = velocitiesProperty();
    
    Molecules changed_mols;
    changed_mols.reserve(nmols);
    
    for (int i=0; i<nmols; ++i)
    {
        MolNum molnum = molgroup.molNumAt(i);
        
        const ViewsOfMol &mol = molecules[molnum];
        
        AtomCoords coords = mol.data().property(coords_property)
                                      .asA<AtomCoords>();

        AtomVelocities vels;
        
        if (mol.data().hasProperty(vels_property))
        {
            vels = mol.data().property(vels_property)
                             .asA<AtomVelocities>();
        }
        else
        {
            vels = AtomVelocities(mol.data().info());
        }
                                          
        if (mol.selectedAll())
        {
            coords.copyFrom(coords_array[i]);
            vels.copyFrom( ::getVelocities(mom_array[i],mass_array[i]) );
        }
        else
        {
            coords.copyFrom(coords_array[i], mol.selection());
            vels.copyFrom( ::getVelocities(mom_array[i],mass_array[i]), 
                           mol.selection() );
        }

        changed_mols.add( mol.molecule().edit()
                             .setProperty(coords_property, coords)
                             .setProperty(vels_property, vels)
                             .commit() );
    }
    
    IntegratorWorkspace::pvt_update(changed_mols);
}

/** Save both the coordinates and velocities back to the system */
void AtomicVelocityWorkspace::commitBufferedCoordinatesAndVelocities(  QVector < QVector< QVector< Vector > > > &buffered_coords )
{
    int nmols = atom_coords.count();
    
    const MoleculeGroup &molgroup = moleculeGroup();
    const Molecules &molecules = molgroup.molecules();
    
    BOOST_ASSERT( molgroup.nMolecules() == nmols );
    
    const QVector<Vector> *coords_array = atom_coords.constData();
    const QVector<Vector> *mom_array = atom_momenta.constData();
    const QVector<double> *mass_array = atom_masses.constData();
    
    PropertyName coords_property = coordinatesProperty();
    PropertyName vels_property = velocitiesProperty();
    
    Molecules changed_mols;
    changed_mols.reserve(nmols);
    
    //qDebug() << " buffered_coords has " << buffered_coords.size() << " elements ";// number of mols
    //qDebug() << " buffered_coords[0] has " << buffered_coords[0].size() << " elements ";// number of atoms in mol 1
    //qDebug() << " buffered_coords[0][0] has " << buffered_coords[0][0].size() << " elements "; // number of coords for atom

    for (int i=0; i<nmols; ++i)
    {
      //qDebug() << " Doing mol " << i ;
        MolNum molnum = molgroup.molNumAt(i);
        
        const ViewsOfMol &mol = molecules[molnum];
        
        AtomCoords coords = mol.data().property(coords_property)
                                      .asA<AtomCoords>();

	//qDebug() << " buffered_coords[0][ " << i << " ] has " << buffered_coords[0][i].size() << " elements ";
	
	QVector< AtomCoords > buffered_molcoords;
	
        AtomVelocities vels;
        
        if (mol.data().hasProperty(vels_property))
        {
            vels = mol.data().property(vels_property)
                             .asA<AtomVelocities>();
        }
        else
        {
            vels = AtomVelocities(mol.data().info());
        }
                                          
        if (mol.selectedAll())
        {
            coords.copyFrom(coords_array[i]);
            vels.copyFrom( ::getVelocities(mom_array[i],mass_array[i]) );

	    for (int k=0; k < buffered_coords.size() ; k++)
	      {
	    	AtomCoords framecoords = AtomCoords( mol.data().info() );
	    	framecoords.copyFrom( buffered_coords[k][i] );
	    	buffered_molcoords.append( framecoords );
	      }
        }
        else
        {
            coords.copyFrom(coords_array[i], mol.selection());
            vels.copyFrom( ::getVelocities(mom_array[i],mass_array[i]), 
                           mol.selection() );
        }

	MolEditor editmol = mol.molecule().edit();
	editmol.setProperty(coords_property, coords);
	editmol.setProperty(vels_property, vels);

	for (int k=0; k < buffered_molcoords.size() ; k++)
	  {
	    PropertyName buffered_property = PropertyName( "buffered_coord_" + QString::number(k) );
	    editmol.setProperty(buffered_property, buffered_molcoords[k] );
	  }
	changed_mols.add( editmol.commit() );

        //changed_mols.add( mol.molecule().edit()
	//                             .setProperty(coords_property, coords)
	//                             .setProperty(vels_property, vels)
	//                             .commit() );

	buffered_molcoords.clear();
	//qDebug() << " Looping ";
    }
    
    IntegratorWorkspace::pvt_update(changed_mols);
}
