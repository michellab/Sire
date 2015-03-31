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

#include "qmff.h"
#include "qmprogram.h"

#include "SireMol/partialmolecule.h"
#include "SireMol/moleculegroup.h"

#include "SireFF/energytable.h"
#include "SireFF/forcetable.h"
#include "SireFF/fieldtable.h"
#include "SireFF/potentialtable.h"

#include "SireMM/cljprobe.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>

using namespace Squire;
using namespace SireFF;
using namespace SireMol;
using namespace SireBase;
using namespace SireUnits::Dimension;
using namespace SireStream;

static const RegisterMetaType<QMFF> r_qmff;

/** Serialise to a binary datastream */
QDataStream SQUIRE_EXPORT &operator<<(QDataStream &ds, const QMFF &qmff)
{
    writeHeader(ds, r_qmff, 1);
    
    SharedDataStream sds(ds);
    
    sds << static_cast<const G1FF&>(qmff)
        << static_cast<const QMPotential&>(qmff)
        << qmff.qmmols;
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream SQUIRE_EXPORT &operator>>(QDataStream &ds, QMFF &qmff)
{
    VersionID v = readHeader(ds, r_qmff);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> static_cast<G1FF&>(qmff)
            >> static_cast<QMPotential&>(qmff)
            >> qmff.qmmols;

        qmff._pvt_updateName();
    }
    else
        throw version_error(v, "1", r_qmff, CODELOC);
        
    return ds;
}

/** Construct an empty, unnamed QM forcefield */
QMFF::QMFF() : ConcreteProperty<QMFF,G1FF>(), FF3D(), QMPotential()
{
    this->_pvt_updateName();
}

/** Construct an empty QM forcefield called 'name' */
QMFF::QMFF(const QString &name)
     : ConcreteProperty<QMFF,G1FF>(), FF3D(), QMPotential()
{
    FF::setName(name);
}

/** Copy constructor */
QMFF::QMFF(const QMFF &other)
     : ConcreteProperty<QMFF,G1FF>(other), FF3D(other), QMPotential(other),
       ffcomponents(other.ffcomponents),
       qmmols(other.qmmols)
{}

/** Destructor */
QMFF::~QMFF()
{}

/** Copy assignment operator */
QMFF& QMFF::operator=(const QMFF &other)
{
    if (this != &other)
    {
        G1FF::operator=(other);
        FF3D::operator=(other);
        QMPotential::operator=(other);

        ffcomponents = other.ffcomponents;
        qmmols = other.qmmols;
    }
    
    return *this;
}

/** Comparison operator */
bool QMFF::operator==(const QMFF &other) const
{
    return FF::operator==(other);
}

/** Comparison operator */
bool QMFF::operator!=(const QMFF &other) const
{
    return FF::operator!=(other);
}

/** Return the space within which the QM molecules exist */
const Space& QMFF::space() const
{
    return QMPotential::space();
}

/** Return the QM program that will be used to calculate the 
    energies and forces on the molecules */
const QMProgram& QMFF::quantumProgram() const
{
    return QMPotential::quantumProgram();
}

/** Return the absolute value of the energy which is considered
    as zero (on the relative energy scale used by this potential).
    A relative scale is used so that the QM energy can be shifted
    so that it is comparable to an MM energy */
MolarEnergy QMFF::zeroEnergy() const
{
    return QMPotential::zeroEnergy();
}

/** Set the space within which the QM molecules exist */
bool QMFF::setSpace(const Space &space)
{
    return QMPotential::setSpace(space);
}

/** Set the QM program that will be used to calculate the 
    energies and forces */
bool QMFF::setQuantumProgram(const QMProgram &qmprog)
{
    return QMPotential::setQuantumProgram(qmprog);
}

/** Set the absolute value of the energy which is considered
    as zero (on the relative energy scale used by this potential).
    A relative scale is used so that the QM energy can be shifted
    so that it is comparable to an MM energy */
bool QMFF::setZeroEnergy(MolarEnergy zero_energy)
{
    return QMPotential::setZeroEnergy(zero_energy);
}

/** Set the property 'name' to the value 'value'

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
bool QMFF::setProperty(const QString &name, const Property &value)
{
    return QMPotential::setProperty(name, value);
}

/** Return the value of the property with name 'name'

    \throw SireBase::missing_property
*/
const Property& QMFF::property(const QString &name) const
{
    return QMPotential::property(name);
}

/** Return whether or not this forcefield contains a property
    called 'name' */
bool QMFF::containsProperty(const QString &name) const
{
    return QMPotential::containsProperty(name);
}

/** Return the properties available in this forcefield (and their values) */
const Properties& QMFF::properties() const
{
    return QMPotential::properties();
}

/** Return the energy components available to this forcefield */
const QMFF::Components& QMFF::components() const
{
    return ffcomponents;
}

/** Trigger a complete recalculation of the QM energy */
void QMFF::mustNowRecalculateFromScratch()
{
    //QM energies are always recalculated from scratch
    this->setDirty();
}

/** Function used to tell the forcefield that the underlying QM
    potential has changed */
void QMFF::changedPotential()
{
    G1FF::incrementVersion();
    this->mustNowRecalculateFromScratch();
}

void QMFF::energy(EnergyTable &energytable, double scale_energy)
{
    throw SireError::incomplete_code( QObject::tr(
            "QMFF does not yet support energy calculations!"), CODELOC );
}

void QMFF::energy(EnergyTable &energytable, const Symbol &symbol,
		  double scale_energy)
{
    throw SireError::incomplete_code( QObject::tr(
            "QMFF does not yet support energy calculations!"), CODELOC );
}

/** Calculate the QM forces on the molecules in this forcefield
    and add the results to the forces for the molecules contained
    in the table 'forcetable' - this scales the forces by
    the optional 'scale_force' */
void QMFF::force(ForceTable &forcetable, double scale_force)
{
    QMPotential::calculateForce(qmmols, forcetable, scale_force);
}

/** Calculate the QM forces on the molecules in this forcefield
    and add the results to the forces for the molecules contained
    in the table 'forcetable' - this scales the forces by
    the optional 'scale_force' */
void QMFF::force(ForceTable &forcetable, const Symbol &symbol,
                 double scale_force)
{
    QMPotential::calculateForce(qmmols, forcetable,
                                symbol, this->components(),
                                scale_force);
}

////
//// Virtual functions from SireFF::FF
////

/** Recalculate the QM energy */
void QMFF::recalculateEnergy()
{
    //QM energies are always recalculated from scratch
    QMEnergy nrg(0);
    QMPotential::calculateEnergy(qmmols, nrg);
    
    this->components().setEnergy(*this, nrg);
    
    this->setClean();
}

/** Function used to update the symbols representing the forcefield
    components whenever the name of the forcefield changes */
void QMFF::_pvt_updateName()
{
    ffcomponents = Components(this->name());
    G1FF::_pvt_updateName();
}

////
//// Virtual functions from SireFF::G1FF
////

/** Record the fact that the molecule 'mol' has been added to this forcefield 

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void QMFF::_pvt_added(const SireMol::PartialMolecule &molecule, 
                      const SireBase::PropertyMap &map)
{
    //add the molecule (don't record changes as everything
    //is recalculated from scratch)
    qmmols.add(molecule, map, *this, false);
    G1FF::setDirty();
}

/** Record the fact that the molecule 'mol' has been removed from this forcefield */
void QMFF::_pvt_removed(const SireMol::PartialMolecule &molecule)
{
    //remove the molecule, again without recording changes
    qmmols.remove(molecule, *this, false);
    G1FF::setDirty();
}

/** Record that fact that the molecule 'molecule' has been changed

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void QMFF::_pvt_changed(const SireMol::Molecule &molecule, bool auto_update)
{
    //change the molecule, again without recording the change
    qmmols.change(molecule, *this, false);
    G1FF::setDirty();
}

/** Record that the provided list of molecules have changed 

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void QMFF::_pvt_changed(const QList<SireMol::Molecule> &mols, bool auto_update)
{
    QMPotential::Molecules old_mols;
    
    try
    {
        for (QList<SireMol::Molecule>::const_iterator it = mols.constBegin();
             it != mols.constEnd();
             ++it)
        {
            qmmols.change(*it, *this, false);
        }
        
        G1FF::setDirty();
    }
    catch(...)
    {
        //restore the state
        qmmols = old_mols;
        throw;
    }
}

/** Record that all of the molecules have been removed */
void QMFF::_pvt_removedAll()
{
    qmmols.clear();
    G1FF::setDirty();
}

/** Return whether or not the supplied property map contains different
    properties for the molecule with number 'molnum' */       
bool QMFF::_pvt_wouldChangeProperties(SireMol::MolNum molnum, 
                                      const SireBase::PropertyMap &map) const
{
    return qmmols.wouldChangeProperties(molnum, map);
}

/** Return the command file that would be used to calculate the energy
    of the molecules in this forcefield */
QString QMFF::energyCommandFile() const
{
    return QMPotential::energyCommandFile(qmmols);
}

/** Return the command file that would be used to calculate the forces
    on the molecules in this forcefield */
QString QMFF::forceCommandFile(const ForceTable &forcetable) const
{
    return QMPotential::forceCommandFile(qmmols, forcetable);
}

/** Return the command file that would be used to calculate the potential
    of the molecules in this forcefield */
QString QMFF::potentialCommandFile(const PotentialTable &potentialtable,
                                   const SireFF::Probe &probe) const
{
    return QMPotential::potentialCommandFile(qmmols, potentialtable, probe);
}

/** Return the command file that would be used to calculate the potential
    of the molecules in this forcefield */
QString QMFF::potentialCommandFile(const PotentialTable &potentialtable) const
{
    return QMFF::potentialCommandFile(potentialtable, QMPotential::Probe());
}

/** Return the command file that would be used to calculate the fields
    of the molecules in this forcefield */
QString QMFF::fieldCommandFile(const FieldTable &fieldtable, 
                               const SireFF::Probe &probe) const
{
    return QMPotential::fieldCommandFile(qmmols, fieldtable, probe);
}

/** Return the command file that would be used to calculate the fields
    of the molecules in this forcefield */
QString QMFF::fieldCommandFile(const FieldTable &fieldtable) const
{
    return QMFF::fieldCommandFile(fieldtable, QMPotential::Probe());
}

/** Calculate the field from this forcefield in the passed fieldtable */
void QMFF::field(FieldTable &fieldtable, const SireFF::Probe &probe, double scale_field)
{
    if (scale_field != 0)
        QMPotential::calculateField(qmmols, fieldtable, probe, scale_field);
}

/** Calculate the field from this forcefield in the passed fieldtable */
void QMFF::field(FieldTable &fieldtable, const Symbol &component,
                 const SireFF::Probe &probe, double scale_field)
{
    if (scale_field != 0)
        QMPotential::calculateField(qmmols, fieldtable, probe, component, 
                                    this->components(), scale_field);
}

/** Calculate the potential from this forcefield in the passed potentialtable */
void QMFF::potential(PotentialTable &potentialtable, const SireFF::Probe &probe,
                     double scale_potential)
{
    if (scale_potential != 0)
        QMPotential::calculatePotential(qmmols, potentialtable, probe, scale_potential);
}

/** Calculate the potential from this forcefield in the passed potentialtable */
void QMFF::potential(PotentialTable &potentialtable, const Symbol &component,
                     const SireFF::Probe &probe, double scale_potential)
{
    if (scale_potential != 0)
        QMPotential::calculatePotential(qmmols, potentialtable, probe, component,
                                        this->components(), scale_potential);
}

/** Calculate the field from this forcefield in the passed fieldtable */
void QMFF::field(FieldTable &fieldtable, double scale_field)
{
    QMFF::field(fieldtable, QMPotential::Probe(), scale_field);
}

/** Calculate the field from this forcefield in the passed fieldtable */
void QMFF::field(FieldTable &fieldtable, const Symbol &component,
                 double scale_field)
{
    QMFF::field(fieldtable, component, QMPotential::Probe(), scale_field);
}
           
/** Calculate the potential from this forcefield in the passed potentialtable */
void QMFF::potential(PotentialTable &potentialtable, double scale_potential)
{
    QMFF::potential(potentialtable, QMPotential::Probe(), scale_potential);
}

/** Calculate the potential from this forcefield in the passed potentialtable */
void QMFF::potential(PotentialTable &potentialtable, const Symbol &component,
                     double scale_potential)
{
    QMFF::potential(potentialtable, component, QMPotential::Probe(), scale_potential);
}

const char* QMFF::typeName()
{
    return QMetaType::typeName( qMetaTypeId<QMFF>() );
}
