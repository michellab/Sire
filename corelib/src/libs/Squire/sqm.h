/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2014  Christopher Woods
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

#ifndef SQUIRE_SQM_H
#define SQUIRE_SQM_H

#include <QHash>
#include <QString>
#include <QProcess>

#include "qmprogram.h"
#include "latticecharges.h"

SIRE_BEGIN_HEADER

namespace Squire 
{
class SQM;
}

SQUIRE_EXPORT QDataStream& operator<<(QDataStream&, const Squire::SQM&);
SQUIRE_EXPORT QDataStream& operator>>(QDataStream&, Squire::SQM&);

class QFile;

namespace SireBase
{
class TempDir;
}

namespace Squire
{

/** This is a wrapper that allows SQM to be used to calculate
    QM and QM/MM energies (SQM is the semiempirical QM program
    that comes free with AmberTools)
    
    @author Christopher Woods
*/
class SQUIRE_EXPORT SQM : public SireBase::ConcreteProperty<SQM,QMProgram>
{

friend SQUIRE_EXPORT QDataStream& ::operator<<(QDataStream&, const SQM&);
friend SQUIRE_EXPORT QDataStream& ::operator>>(QDataStream&, SQM&);

public:
    SQM();
    SQM(const QString &SQM);
    
    SQM(const SQM &other);
    
    ~SQM();
    
    static const char* typeName();
    
    SQM& operator=(const SQM &other);
    
    bool operator==(const SQM &other) const;
    bool operator!=(const SQM &other) const;
    
    QString toString() const;
    
    void setExecutable(const QString &SQM_exe);
    void setEnvironment(const QString &variable, const QString &value);
    
    void setMaximumRunTime(int max_runtime);
    int maximumRunTime() const;
    
    void setMaximumNumberOfSQMInputLines(int numlines);
    void setExpectedNumberOfQMAtoms(int natoms);
    
    int maximumNumberOfSQMInputLines() const;
    int expectedNumberOfQMAtoms() const;
    
    int numberOfMMAtomsLimit() const;
    int numberOfMMAtomsLimit(int num_qm_atoms) const;
    
    QString executable() const;
    
    const QHash<QString,QString>& environment() const;
    
    QString environment(const QString &variable) const;
    
    void setMethod(const QString &method);
    
    const QString& method() const;
    
    void setTotalCharge(int charge);
    int totalCharge() const;
    
    void setEnergyTemplate(const QString &energy_template);
    
    const QString& energyTemplate() const;
    
    void setForceTemplate(const QString &force_template);
    
    const QString& forceTemplate() const;

    bool supportsLatticeCharges() const
    {
        return true;
    }

protected:
    double calculateEnergy(const QMPotential::Molecules &molecules,
                           int ntries = 5) const;
    double calculateEnergy(const QMPotential::Molecules &molecules,
                           const LatticeCharges &lattice_charges,
                           int ntries = 5) const;

    QHash<QString,double> calculatePotential(const QString &cmdfile,
                                            int ntries=5) const;

    void calculatePotential(const QMPotential::Molecules &molecules,
                            PotentialTable &pottable,
                            const SireFF::Probe &probe,
                            double scale_potential=1, 
                            int ntries=5) const;

    QVector<SireUnits::Dimension::MolarEnergy> 
                    calculatePotential(const QMPotential::Molecules &molecules,
                                       const LatticeCharges &lattice_charges,
                                       PotentialTable &pottable,
                                       const SireFF::Probe &probe,
                                       double scale_potential=1, 
                                       int ntries=5) const;

    QString energyCommandFile(const QMPotential::Molecules &molecules) const;
    QString energyCommandFile(const QMPotential::Molecules &molecules,
                              const LatticeCharges &lattice_charges) const;
    
    QString forceCommandFile(const QMPotential::Molecules &molecules,
                             const ForceTable &forcetable) const;
    QString forceCommandFile(const QMPotential::Molecules &molecules,
                             const LatticeCharges &lattice_charges,
                             const ForceTable &forcetable) const;
    
    QString fieldCommandFile(const QMPotential::Molecules &molecules,
                             const FieldTable &fieldtable,
                             const SireFF::Probe &probe) const;
    QString fieldCommandFile(const QMPotential::Molecules &molecules,
                             const LatticeCharges &lattice_charges,
                             const FieldTable &fieldtable,
                             const SireFF::Probe &probe) const;
    
    QString potentialCommandFile(const QMPotential::Molecules &molecules,
                                 const PotentialTable &pottable,
                                 const SireFF::Probe &probe) const;
    QString potentialCommandFile(const QMPotential::Molecules &molecules,
                                 const LatticeCharges &lattice_charges,
                                 const PotentialTable &pottable,
                                 const SireFF::Probe &probe) const;

private:
    QString createCommandFile(QString cmd_template,
                        const QMPotential::Molecules &molecules,
                        const LatticeCharges &lattice_charges = LatticeCharges()) const;

    QString writeShellFile(const SireBase::TempDir &tempdir) const;

    double extractEnergy(QFile &SQM_output) const;

    double calculateEnergy(const QString &cmd_file, int ntries) const;

    /** The environmental variables to hold when running SQM */
    QHash<QString,QString> env_variables;
    
    /** The full path to the SQM executable to run (including
        any necessary command line arguments) */
    QString sqm_exe;
    
    /** The QM method to use to calculate the energy */
    QString qm_method;
    
    /** The template command file used for the energy calculations.
        The basis set, QM method and atom coordinates are substituted
        into this template */
    QString energy_template;
    
    /** The template command file used for the force calculations.
        The basis set, QM method and atom coordinates are substituted
        into this template */
    QString force_template;
    
    /** The total charge of the system */
    qint32 total_charge;
    
    /** The maximum amount of time to wait for a SQM
        job to complete (15 minutes) in milliseconds */
    quint32 max_sqm_runtime;
    
    /** The maximum number of lines that can be parsed in an SQM command file */
    qint32 max_sqm_lines;
    
    /** The expected number of QM atoms. This is used, together with the maximum
        number of SQM command file lines, to work out the maximum number of MM
        atoms that can be supported */
    qint32 expected_n_qm;
};

}

Q_DECLARE_METATYPE( Squire::SQM )

SIRE_EXPOSE_CLASS( Squire::SQM )

SIRE_END_HEADER

#endif
