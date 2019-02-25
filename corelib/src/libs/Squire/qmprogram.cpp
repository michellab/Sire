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

#include <QMutex>

#include "qmprogram.h"
#include "latticecharges.h"

#include "SireMol/molecule.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace Squire;
using namespace SireBase;
using namespace SireUnits::Dimension;
using namespace SireStream;

///////
/////// Implementation of QMProgram
///////

static const RegisterMetaType<QMProgram> r_qmprog( MAGIC_ONLY,
                                                   "Squire::QMProgram" );

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const QMProgram &qmprog)
{
    writeHeader(ds, r_qmprog, 1);
    
    ds << static_cast<const Property&>(qmprog);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, QMProgram &qmprog)
{
    VersionID v = readHeader(ds, r_qmprog);
    
    if (v == 1)
    {
        ds >> static_cast<Property&>(qmprog);
    }
    else
        throw version_error(v, "1", r_qmprog, CODELOC);
        
    return ds;
}

/** Constructor */
QMProgram::QMProgram() : Property()
{}

/** Copy constructor */
QMProgram::QMProgram(const QMProgram &other) : Property(other)
{}

/** Destructor */
QMProgram::~QMProgram()
{}

/** Return the maximum number of MM atoms supported by this QM program. This
    returns -1 if there is no limit */
int QMProgram::numberOfMMAtomsLimit() const
{
    return -1;
}

/** Return the maximum number of QM atoms supported by this QM program, 
    given the supplied number of MM atoms. This returns -1 if there is
    no limit */
int QMProgram::numberOfMMAtomsLimit(int num_qm_atoms) const
{
    return -1;
}

/** Calculate the charges on the molecule 'molecule' using the properties
    specified in the passed property map */
AtomCharges QMProgram::calculateCharges(const Molecule &molecule,
                                        const PropertyMap &map) const
{
    throw SireError::unsupported( QObject::tr(
            "Calculating the charges via this interface (%1) is not "
            "yet supported.")
                .arg(this->what()), CODELOC );
}

/** Calculate the charges on the molecule 'molecule' using the default
    property locations */
AtomCharges QMProgram::calculateCharges(const Molecule &molecule) const
{
    return this->calculateCharges( molecule, PropertyMap() );
}

/** Return the QM energy of the molecules 'molecules' surrounded by the 
    field of point charges 'lattice_charges' */
double QMProgram::calculateEnergy(const QMPotential::Molecules &molecules,
                                  const LatticeCharges &lattice_charges,
                                  int) const
{
    throw SireError::unsupported( QObject::tr(
        "This QM program (%1) does not support the use of point lattice charges.")
            .arg(this->what()), CODELOC );
}

/** Calculate the forces on the passed molecules, and place them into the passed
    force table, optionally scaled by 'scale_force' */
void QMProgram::calculateForce(const QMPotential::Molecules&, ForceTable &,
                               double, int) const
{
    throw SireError::unsupported( QObject::tr(
        "This QM program (%1) does not support the calculation of forces.")
            .arg(this->what()), CODELOC );
}

/** Calculate the forces on the passed molecules, and place them into the passed
    force table, optionally scaled by 'scale_force', and return the accompanying
    forces on the passed lattice points, also scaled by 'scale_force' */
QVector<Vector> QMProgram::calculateForce(const QMPotential::Molecules&,
                                          const LatticeCharges&, ForceTable&,
                                          double, int) const
{
    throw SireError::unsupported( QObject::tr(
        "This QM program (%1) does not support the calculation of forces with "
        "associated lattice charges.")
            .arg(this->what()), CODELOC );
    
    return QVector<Vector>();
}

/** Calculate the field around the passed molecules, and place them into the passed
    field table, optionally scaled by 'scale_field' */
void QMProgram::calculateField(const QMPotential::Molecules&, FieldTable&,
                               const SireFF::Probe&, double, int) const
{
    throw SireError::unsupported( QObject::tr(
        "This QM program (%1) does not support the calculation of fields.")
            .arg(this->what()), CODELOC );
}

/** Calculate the fields around the passed molecules, and place them into the passed
    field table, optionally scaled by 'scale_field', and return the accompanying
    fields on the passed lattice points, also scaled by 'scale_field' */
QVector<Vector> QMProgram::calculateField(const QMPotential::Molecules&,
                                          const LatticeCharges&, FieldTable&,
                                          const SireFF::Probe&, double, int) const
{
    throw SireError::unsupported( QObject::tr(
        "This QM program (%1) does not support the calculation of fields with "
        "associated lattice charges.")
            .arg(this->what()), CODELOC );
    
    return QVector<Vector>();
}

/** Calculate the potential around the passed molecules, and place them into the passed
    potential table, optionally scaled by 'scale_potential' */
void QMProgram::calculatePotential(const QMPotential::Molecules&, PotentialTable&,
                                   const SireFF::Probe&, double, int) const
{
    throw SireError::unsupported( QObject::tr(
        "This QM program (%1) does not support the calculation of potentials.")
            .arg(this->what()), CODELOC );
}

/** Calculate the potentials around the passed molecules, and place them into the passed
    potential table, optionally scaled by 'scale_potential', and return the accompanying
    potentials on the passed lattice points, also scaled by 'scale_potential' */
QVector<MolarEnergy> QMProgram::calculatePotential(const QMPotential::Molecules&,
                                                   const LatticeCharges&, 
                                                   PotentialTable&,
                                                   const SireFF::Probe&, 
                                                   double, int) const
{
    throw SireError::unsupported( QObject::tr(
        "This QM program (%1) does not support the calculation of potentials with "
        "associated lattice charges.")
            .arg(this->what()), CODELOC );
    
    return QVector<MolarEnergy>();
}

/** Return the command file that would be used to calculate the atomic
    partial charges of the passed molecule */
QString QMProgram::chargeCommandFile(const Molecule &molecule,
                                     const PropertyMap &map) const
{
    throw SireError::unsupported( QObject::tr(
            "Calculating the charges via this interface (%1) is not "
            "yet supported.")
                .arg(this->what()), CODELOC );
}

/** Return the command file that would be used to calculate the atomic
    partial charges of the passed molecule */
QString QMProgram::chargeCommandFile(const Molecule &molecule) const
{
    return this->chargeCommandFile(molecule, PropertyMap());
}

/** Return the command file that would be used to calculate the energy
    of the molecules in 'molecules' in the field of point charges in
    'lattice_charges' */
QString QMProgram::energyCommandFile(const QMPotential::Molecules &molecules,
                                    const LatticeCharges &lattice_charges) const
{
    throw SireError::unsupported( QObject::tr(
        "This QM program (%1) does not support the use of point lattice charges.")
            .arg(this->what()), CODELOC );
}

/** Return the command file that would be used to calculate the forces
    of the molecules in 'molecules' in the field of point charges in
    'lattice_charges' (and the forces on the charges themselves) */
QString QMProgram::forceCommandFile(const QMPotential::Molecules&,
                                    const LatticeCharges&,
                                    const ForceTable&) const
{
    throw SireError::unsupported( QObject::tr(
        "This QM program (%1) does not support the use of point lattice charges.")
            .arg(this->what()), CODELOC );
}

/** Return the command file that would be used to calculate the fields
    of the molecules in 'molecules' in the field of point charges in
    'lattice_charges' */
QString QMProgram::fieldCommandFile(const QMPotential::Molecules&,
                                    const LatticeCharges&,
                                    const FieldTable&,
                                    const SireFF::Probe&) const
{
    throw SireError::unsupported( QObject::tr(
        "This QM program (%1) does not support the use of point lattice charges.")
            .arg(this->what()), CODELOC );
}

/** Return the command file that would be used to calculate the fields
    of the molecules in 'molecules' in the field of point charges in
    'lattice_charges' */
QString QMProgram::potentialCommandFile(const QMPotential::Molecules&,
                                        const LatticeCharges&,
                                        const PotentialTable&,
                                        const SireFF::Probe&) const
{
    throw SireError::unsupported( QObject::tr(
        "This QM program (%1) does not support the use of point lattice charges.")
            .arg(this->what()), CODELOC );
}

///////
/////// Implementation of NullQM
///////

static const RegisterMetaType<NullQM> r_nullqm;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const NullQM &nullqm)
{
    writeHeader(ds, r_nullqm, 1);
    ds << static_cast<const QMProgram&>(nullqm);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, NullQM &nullqm)
{
    VersionID v = readHeader(ds, r_nullqm);
    
    if (v == 1)
    {
        ds >> static_cast<QMProgram&>(nullqm);
    }
    else
        throw version_error(v, "1", r_nullqm, CODELOC);

    return ds;
}

/** Constructor */
NullQM::NullQM() : ConcreteProperty<NullQM,QMProgram>()
{}

static SharedPolyPointer<NullQM> shared_null;

const NullQM& QMProgram::null()
{
    return *(create_shared_null<NullQM>());
}

/** Copy constructor */
NullQM::NullQM(const NullQM &other)
         : ConcreteProperty<NullQM,QMProgram>(other)
{}

/** Destructor */
NullQM::~NullQM()
{}

/** Copy assignment operator */
NullQM& NullQM::operator=(const NullQM &other)
{
    return *this;
}

/** Comparison operator */
bool NullQM::operator==(const NullQM &other) const
{
    return true;
}

/** Comparison operator */
bool NullQM::operator!=(const NullQM &other) const
{
    return false;
}

/** Calculate and return the QM energy of all of the molecules
    in 'molecules' */
double NullQM::calculateEnergy(const QMPotential::Molecules&, int) const
{
    return 0;
}

/** Return the QM energy of the molecules 'molecules' surrounded by the 
    field of point charges 'lattice_charges' */
double NullQM::calculateEnergy(const QMPotential::Molecules&,
                               const LatticeCharges&,
                               int) const
{
    return 0;
}

/** Return the command file that would be used to calculate the energy of
    the molecules in 'molecules' */
QString NullQM::energyCommandFile(const QMPotential::Molecules&) const
{
    return QString::null;
}

/** Return the command file that would be used to calculate the forces on
    the molecules in 'molecules' */
QString NullQM::forceCommandFile(const QMPotential::Molecules&,
                                 const ForceTable&) const
{
    return QString::null;
}

/** Return the command file that would be used to calculate the fields on
    the molecules in 'molecules' */
QString NullQM::fieldCommandFile(const QMPotential::Molecules&,
                                 const FieldTable&, const SireFF::Probe&) const
{
    return QString::null;
}

/** Return the command file that would be used to calculate the fields on
    the molecules in 'molecules' */
QString NullQM::potentialCommandFile(const QMPotential::Molecules&,
                                     const PotentialTable&, const SireFF::Probe&) const
{
    return QString::null;
}

/** Return the command file that would be used to calculate the energy
    of the molecules in 'molecules' in the field of point charges in
    'lattice_charges' */
QString NullQM::energyCommandFile(const QMPotential::Molecules&,
                                  const LatticeCharges &lattice_charges) const
{
    return QString::null;
}

/** Return the command file that would be used to calculate the forces
    of the molecules in 'molecules' in the field of point charges in
    'lattice_charges' (and the forces on the charges themselves) */
QString NullQM::forceCommandFile(const QMPotential::Molecules&,
                                 const LatticeCharges&,
                                 const ForceTable&) const
{
    return QString::null;
}

/** Return the command file that would be used to calculate the fields
    of the molecules in 'molecules' in the field of point charges in
    'lattice_charges' (and the fields on the charges themselves) */
QString NullQM::fieldCommandFile(const QMPotential::Molecules&,
                                 const LatticeCharges&,
                                 const FieldTable&, const SireFF::Probe&) const
{
    return QString::null;
}

/** Return the command file that would be used to calculate the potentials
    of the molecules in 'molecules' in the field of point charges in
    'lattice_charges' (and the potentials on the charges themselves) */
QString NullQM::potentialCommandFile(const QMPotential::Molecules&,
                                     const LatticeCharges&,
                                     const PotentialTable&, const SireFF::Probe&) const
{
    return QString::null;
}

const char* NullQM::typeName()
{
    return QMetaType::typeName( qMetaTypeId<NullQM>() );
}

NullQM* NullQM::clone() const
{
    return new NullQM(*this);
}

