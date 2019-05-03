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

#include "qmmmff.h"
#include "qmprogram.h"

#include "SireMM/cljprobe.h"

#include "SireFF/energytable.h"
#include "SireFF/forcetable.h"
#include "SireFF/fieldtable.h"
#include "SireFF/potentialtable.h"

#include "SireBase/booleanproperty.h"

#include "SireFF/errors.h"
#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>

using namespace Squire;
using namespace SireMM;
using namespace SireFF;
using namespace SireBase;
using namespace SireUnits::Dimension;
using namespace SireStream;

static const RegisterMetaType<QMMMFF> r_qmmmff;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const QMMMFF &qmmmff)
{
    writeHeader(ds, r_qmmmff, 2);
    
    SharedDataStream sds(ds);
    
    sds << static_cast<const G2FF&>(qmmmff)
        << static_cast<const QMMMElecEmbedPotential&>(qmmmff)
        << qmmmff.qmmols
        << qmmmff.mmmols
        << qmmmff.intermolecular_only;
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, QMMMFF &qmmmff)
{
    VersionID v = readHeader(ds, r_qmmmff);
    
    if (v <= 2)
    {
        SharedDataStream sds(ds);
        
        sds >> static_cast<G2FF&>(qmmmff)
            >> static_cast<QMMMElecEmbedPotential&>(qmmmff)
            >> qmmmff.qmmols
            >> qmmmff.mmmols;
        
        if (v == 2)
            sds >> qmmmff.intermolecular_only;
        else
            qmmmff.intermolecular_only = false;
        
        qmmmff._pvt_updateName();
    }
    else
        throw version_error(v, "1,2", r_qmmmff, CODELOC);
        
    return ds;
}

/** Construct an empty, unnamed QM/MM forcefield */
QMMMFF::QMMMFF() : ConcreteProperty<QMMMFF,G2FF>(), FF3D(), 
                   QMMMElecEmbedPotential(), intermolecular_only(false)
{
    this->_pvt_updateName();
}

/** Construct an empty, named QM/MM forcefield */
QMMMFF::QMMMFF(const QString &name) : ConcreteProperty<QMMMFF,G2FF>(), FF3D(), 
                                      QMMMElecEmbedPotential(), intermolecular_only(false)
{
    G2FF::setName(name);
}

/** Copy constructor */
QMMMFF::QMMMFF(const QMMMFF &other)
       : ConcreteProperty<QMMMFF,G2FF>(other), FF3D(),
         QMMMElecEmbedPotential(other),
         ffcomponents(other.ffcomponents),
         qmmols(other.qmmols), mmmols(other.mmmols),
         intermolecular_only(other.intermolecular_only)
{}

/** Destructor */
QMMMFF::~QMMMFF()
{}

/** Copy assignment operator */    
QMMMFF& QMMMFF::operator=(const QMMMFF &other)
{
    if (this != &other)
    {
        G2FF::operator=(other);
        QMMMElecEmbedPotential::operator=(other);
        
        qmmols = other.qmmols;
        mmmols = other.mmmols;
        
        intermolecular_only = other.intermolecular_only;
    }
    
    return *this;
}

/** Comparison operator */
bool QMMMFF::operator==(const QMMMFF &other) const
{
    return FF::operator==(other) and intermolecular_only == other.intermolecular_only;
}

/** Comparison operator */
bool QMMMFF::operator!=(const QMMMFF &other) const
{
    return FF::operator!=(other) or intermolecular_only != other.intermolecular_only;
}

void QMMMFF::throwInvalidGroup(int group_id) const
{
    throw SireError::program_bug( QObject::tr(
        "Program bug! The only groups in QM/MM forcefield are groups 0 and 1. "
        "There is no group %1.").arg(group_id), CODELOC );
}

/** Return the symbols representing the energy components of this forcefield */
const QMMMFF::Components& QMMMFF::components() const
{
    return ffcomponents;
}

/** Return the space within which the QM molecules exist */
const Space& QMMMFF::space() const
{
    return QMMMElecEmbedPotential::space();
}

/** Return the switching function used to provide the nonbonded cutoff
    between the QM and MM regions */
const SwitchingFunction& QMMMFF::switchingFunction() const
{
    return QMMMElecEmbedPotential::switchingFunction();
}

/** Return the QM program that will be used to calculate the 
    energies and forces on the molecules */
const QMProgram& QMMMFF::quantumProgram() const
{
    return QMMMElecEmbedPotential::quantumProgram();
}

/** Return the absolute value of the energy which is considered
    as zero (on the relative energy scale used by this potential).
    A relative scale is used so that the QM energy can be shifted
    so that it is comparable to an MM energy */
MolarEnergy QMMMFF::zeroEnergy() const
{
    return QMMMElecEmbedPotential::zeroEnergy();
}

/** Return the amount by which the MM charges are scaled in the QM/MM interaction */
double QMMMFF::chargeScalingFactor() const
{
    return QMMMElecEmbedPotential::chargeScalingFactor();
}

/** Set the space within which the QM molecules exist */
bool QMMMFF::setSpace(const Space &space)
{
    return QMMMElecEmbedPotential::setSpace(space);
}

/** Set the switching function used to provide the 
    cutoff between the QM and MM regions */
bool QMMMFF::setSwitchingFunction(const SwitchingFunction &switchfunc)
{
    return QMMMElecEmbedPotential::setSwitchingFunction(switchfunc);
}

/** Set the QM program that will be used to calculate the 
    energies and forces */
bool QMMMFF::setQuantumProgram(const QMProgram &qmprog)
{
    return QMMMElecEmbedPotential::setQuantumProgram(qmprog);
}

/** Set the absolute value of the energy which is considered
    as zero (on the relative energy scale used by this potential).
    A relative scale is used so that the QM energy can be shifted
    so that it is comparable to an MM energy */
bool QMMMFF::setZeroEnergy(MolarEnergy zero_energy)
{
    if (intermolecular_only.value())
        //do not set zero energy if we are only calculating intermolecular energies
        return false;
    
    return QMMMElecEmbedPotential::setZeroEnergy(zero_energy);
}

/** Set the scaling factor for the MM charges in the QM/MM interaction */
bool QMMMFF::setChargeScalingFactor(double scale_factor)
{
    return QMMMElecEmbedPotential::setChargeScalingFactor(scale_factor);
}

/** Set whether or not we only calculate the intermolecular energy
    (energy between the QM and MM atoms) */
bool QMMMFF::setIntermolecularOnly(bool on)
{
    if (intermolecular_only.value() != on)
    {
        if (on)
            //there is no zero energy if we are only calculating intermolecular energies
            setZeroEnergy( MolarEnergy(0) );

        intermolecular_only = on;
        mustNowRecalculateFromScratch();
        return true;
    }
    else
        return false;
}

/** Return whether or not we only calculate the intermolecular energy
    (energy between the QM and MM atoms) */
bool QMMMFF::isIntermolecularOnly() const
{
    return intermolecular_only.value();
}

/** Set the property 'name' to the value 'value'

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
bool QMMMFF::setProperty(const QString &name, const Property &value)
{
    if (name == "intermolecularOnly")
        return this->setIntermolecularOnly(value.asA<BooleanProperty>().value());
    else
        return QMMMElecEmbedPotential::setProperty(name, value);
}

/** Return the value of the property with name 'name'

    \throw SireBase::missing_property
*/
const Property& QMMMFF::property(const QString &name) const
{
    if (name == "intermolecularOnly")
        return intermolecular_only;
    else
        return QMMMElecEmbedPotential::property(name);
}

/** Return whether or not this forcefield contains a property
    called 'name' */
bool QMMMFF::containsProperty(const QString &name) const
{
    return name == "intermolecularOnly" or QMMMElecEmbedPotential::containsProperty(name);
}

static Properties global_props;

/** Return the properties available in this forcefield (and their values) */
const Properties& QMMMFF::properties() const
{
    //dangerous, but will be fixed by me getting rid of const Properties& from the API!
    global_props = QMMMElecEmbedPotential::properties();
    global_props.setProperty("intermolecularOnly", intermolecular_only);
    return global_props;
}

/** Signal that this forcefield must recalculate the energy from scratch */
void QMMMFF::mustNowRecalculateFromScratch()
{
    G2FF::setDirty();
}

void QMMMFF::energy(EnergyTable &energytable, double scale_energy)
{
    throw SireError::incomplete_code( QObject::tr(
						  "QMMMFF does not yet support energy calculations!"), CODELOC );
}

void QMMMFF::energy(EnergyTable &energytable, const Symbol &symbol,
		    double scale_energy)
{
    throw SireError::incomplete_code( QObject::tr(
            "QMMMFF does not yet support energy calculations!"), CODELOC );
}

/** Calculate the QM/MM forces on the molecules in this forcefield
    and add the results to the forces for the molecules contained
    in the table 'forcetable' - this scales the forces by
    the optional 'scale_force' */
void QMMMFF::force(ForceTable &forcetable, double scale_force)
{
    if (intermolecular_only.value())
        throw SireError::unsupported( QObject::tr("Support for QM/MM intermolecular "
                 "only forces is not yet available."), CODELOC );
    else
        QMMMElecEmbedPotential::calculateForce(qmmols, mmmols, forcetable, scale_force);
}

/** Calculate the QM/MM forces on the molecules in this forcefield
    and add the results to the forces for the molecules contained
    in the table 'forcetable' - this scales the forces by
    the optional 'scale_force' */
void QMMMFF::force(ForceTable &forcetable, const Symbol &symbol,
                   double scale_force)
{
    if (intermolecular_only.value())
        throw SireError::unsupported( QObject::tr("Support for QM/MM intermolecular "
                 "only forces is not yet available."), CODELOC );
    else
        QMMMElecEmbedPotential::calculateForce(qmmols, mmmols,
                                               forcetable, symbol,
                                               this->components(),
                                               scale_force);
}

/** This function is called whenever the underlying potential is changed */
void QMMMFF::changedPotential()
{
    G2FF::incrementVersion();
    this->mustNowRecalculateFromScratch();
}

////
//// Virtual functions from SireFF::FF
////

/** Recalculate the QM energy */
void QMMMFF::recalculateEnergy()
{
    //QM/MM energies are always recalculated from scratch
    QMEnergy nrg(0);

    QMMMElecEmbedPotential::calculateEnergy(qmmols, mmmols, nrg);
    
    if (intermolecular_only.value())
    {
        //also calculate only the QM energy, and subtract this from the QM/MM energy
        QMEnergy qmnrg(0);
        QMMMElecEmbedPotential::calculateEnergy(qmmols, MMMolecules(), qmnrg);
        
        nrg -= qmnrg;
    }
    
    this->components().setEnergy(*this, nrg);
    this->setClean();
}


/** Function used to update the symbols representing the forcefield
    components whenever the name of the forcefield changes */
void QMMMFF::_pvt_updateName()
{
    ffcomponents = Components(this->name());
    G2FF::_pvt_updateName();
}

////
//// Virtual functions from SireFF::G2FF
////

/** Record the fact that the molecule 'mol' has been added to this forcefield 

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void QMMMFF::_pvt_added(quint32 group_id, 
                        const SireMol::PartialMolecule &molecule, 
                        const SireBase::PropertyMap &map)
{
    //add the molecule (don't record changes as everything
    //is recalculated from scratch)
    if (group_id == 0)
    {
        qmmols.add(molecule, map, *this, false);
        G2FF::setDirty();
    }
    else if (group_id == 1)
    {
        mmmols.add(molecule, map, *this, false);
        G2FF::setDirty();
    }
    else
        throwInvalidGroup(group_id);
}
                
/** Record the fact that the molecule 'mol' has been removed from this forcefield */
void QMMMFF::_pvt_removed(quint32 group_id, 
                          const SireMol::PartialMolecule &molecule)
{
    if (group_id == 0)
    {
        qmmols.remove(molecule, *this, false);
        G2FF::setDirty();
    }
    else if (group_id == 1)
    {
        mmmols.remove(molecule, *this, false);
        G2FF::setDirty();
    }
    else
        throwInvalidGroup(group_id);
}
                  
/** Record that fact that the molecule 'molecule' has been changed

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void QMMMFF::_pvt_changed(quint32 group_id, const SireMol::Molecule &molecule, bool auto_update)
{
    if (group_id == 0)
    {
        qmmols.change(molecule, *this, false);
        G2FF::setDirty();
    }   
    else if (group_id == 1)
    {
        mmmols.change(molecule, *this, false);
        G2FF::setDirty();
    }
    else
        throwInvalidGroup(group_id);
}

/** Record that the provided list of molecules have changed 

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void QMMMFF::_pvt_changed(quint32 group_id, const QList<SireMol::Molecule> &molecules,
                          bool auto_update)
{
    if (group_id == 0)
    {
        QMMolecules old_mols;
    
        try
        {
            for (QList<SireMol::Molecule>::const_iterator it = molecules.constBegin();
                 it != molecules.constEnd();
                 ++it)
            {
                qmmols.change(*it, *this, false);
            }
            
            G2FF::setDirty();
        }
        catch(...)
        {
            //restore the state
            qmmols = old_mols;
            throw;
        }
    }
    else if (group_id == 1)
    {
        MMMolecules old_mols;
    
        try
        {
            for (QList<SireMol::Molecule>::const_iterator it = molecules.constBegin();
                 it != molecules.constEnd();
                 ++it)
            {
                mmmols.change(*it, *this, false);
            }
            
            G2FF::setDirty();
        }
        catch(...)
        {
            //restore the state
            mmmols = old_mols;
            throw;
        }
    }   
    else
        throwInvalidGroup(group_id);
}

/** Record that all of the molecules have been removed */
void QMMMFF::_pvt_removedAll(quint32 group_id)
{
    if (group_id == 0)
    {
        qmmols.clear();
        G2FF::setDirty();
    }
    else if (group_id == 1)
    {
        mmmols.clear();
        G2FF::setDirty();
    }
    else
        throwInvalidGroup(group_id);
}

/** Return whether or not the supplied property map contains different
    properties for the molecule with number 'molnum' */       
bool QMMMFF::_pvt_wouldChangeProperties(quint32 group_id, 
                                        SireMol::MolNum molnum, 
                                        const SireBase::PropertyMap &map) const
{
    if (group_id == 0)
        return qmmols.wouldChangeProperties(molnum, map);
        
    else if (group_id == 1)
        return mmmols.wouldChangeProperties(molnum, map);
        
    else
    {
        throwInvalidGroup(group_id);
        return false;
    }
}

/** Return the command file that would be used to calculate the energy
    of the molecules in this forcefield */
QString QMMMFF::energyCommandFile() const
{
    return QMMMElecEmbedPotential::energyCommandFile(qmmols, mmmols);
}

/** Return the command file that would be used to calculate the forces
    on the molecules in this forcefield */
QString QMMMFF::forceCommandFile(const ForceTable &forcetable) const
{
    return QMMMElecEmbedPotential::forceCommandFile(qmmols, mmmols, forcetable);
}

/** Return the command file that would be used to calculate the potential
    of the molecules in this forcefield */
QString QMMMFF::potentialCommandFile(const PotentialTable &potentialtable,
                                     const SireFF::Probe &probe) const
{
    return QMMMElecEmbedPotential::potentialCommandFile(qmmols, mmmols, 
                                                        potentialtable, probe);
}

/** Return the command file that would be used to calculate the potential
    of the molecules in this forcefield */
QString QMMMFF::potentialCommandFile(const PotentialTable &potentialtable) const
{
    return QMMMFF::potentialCommandFile(potentialtable, QMMMElecEmbedPotential::Probe());
}

/** Return the command file that would be used to calculate the fields
    of the molecules in this forcefield */
QString QMMMFF::fieldCommandFile(const FieldTable &fieldtable,
                                 const SireFF::Probe &probe) const
{
    return QMMMElecEmbedPotential::fieldCommandFile(qmmols, mmmols, fieldtable, probe);
}

/** Return the command file that would be used to calculate the fields
    of the molecules in this forcefield */
QString QMMMFF::fieldCommandFile(const FieldTable &fieldtable) const
{
    return QMMMFF::fieldCommandFile(fieldtable, QMMMElecEmbedPotential::Probe());
}

/** Calculate the field from this forcefield in the passed fieldtable */
void QMMMFF::field(FieldTable &fieldtable, const SireFF::Probe &probe, double scale_field)
{
    if (scale_field != 0)
        QMMMElecEmbedPotential::calculateField(qmmols, mmmols, fieldtable, 
                                               probe, scale_field);
}

/** Calculate the field from this forcefield in the passed fieldtable */
void QMMMFF::field(FieldTable &fieldtable, const Symbol &component,
                   const SireFF::Probe &probe, double scale_field)
{
    if (scale_field != 0)
        QMMMElecEmbedPotential::calculateField(qmmols, mmmols, fieldtable, probe,
                                               component, this->components(),
                                               scale_field);
}

/** Calculate the potential from this forcefield in the passed potentialtable */
void QMMMFF::potential(PotentialTable &potentialtable, const SireFF::Probe &probe,
                       double scale_potential)
{
    if (scale_potential != 0)
        QMMMElecEmbedPotential::calculatePotential(qmmols, mmmols, potentialtable, 
                                                   probe, scale_potential);
}

/** Calculate the potential from this forcefield in the passed potentialtable */
void QMMMFF::potential(PotentialTable &potentialtable, const Symbol &component,
                       const SireFF::Probe &probe, double scale_potential)
{
    if (scale_potential != 0)
        QMMMElecEmbedPotential::calculatePotential(qmmols, mmmols, potentialtable, 
                                                   probe, component, 
                                                   this->components(), scale_potential);
}

/** Calculate the field from this forcefield in the passed fieldtable */
void QMMMFF::field(FieldTable &fieldtable, double scale_field)
{
    QMMMFF::field(fieldtable, QMMMElecEmbedPotential::Probe(), scale_field);
}

/** Calculate the field from this forcefield in the passed fieldtable */
void QMMMFF::field(FieldTable &fieldtable, const Symbol &component,
                   double scale_field)
{
    QMMMFF::field(fieldtable, component, QMMMElecEmbedPotential::Probe(), 
                  scale_field);
}
           
/** Calculate the potential from this forcefield in the passed potentialtable */
void QMMMFF::potential(PotentialTable &potentialtable, double scale_potential)
{
    QMMMFF::potential(potentialtable, QMMMElecEmbedPotential::Probe(), scale_potential);
}

/** Calculate the potential from this forcefield in the passed potentialtable */
void QMMMFF::potential(PotentialTable &potentialtable, const Symbol &component,
                       double scale_potential)
{
    QMMMFF::potential(potentialtable, component, QMMMElecEmbedPotential::Probe(), 
                      scale_potential);
}

const char* QMMMFF::typeName()
{
    return QMetaType::typeName( qMetaTypeId<QMMMFF>() );
}
