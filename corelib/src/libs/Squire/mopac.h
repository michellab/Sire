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

#ifndef SQUIRE_MOPACK_H
#define SQUIRE_MOPACK_H

#include <QHash>
#include <QString>
#include <QStringList>
#include <QProcess>
#include <QList>
#include <QPair>

#include "qmprogram.h"
#include "latticecharges.h"

SIRE_BEGIN_HEADER

namespace Squire 
{
class Mopac;
}

SQUIRE_EXPORT QDataStream& operator<<(QDataStream&, const Squire::Mopac&);
SQUIRE_EXPORT QDataStream& operator>>(QDataStream&, Squire::Mopac&);

class QFile;

namespace SireBase
{
class TempDir;
}

namespace Squire
{

/** This is a wrapper that allows Mopac to be used to calculate
    semiempirical QM energies
    
    @author Christopher Woods
*/
class SQUIRE_EXPORT Mopac : public SireBase::ConcreteProperty<Mopac,QMProgram>
{

friend QDataStream& ::operator<<(QDataStream&, const Mopac&);
friend QDataStream& ::operator>>(QDataStream&, Mopac&);

public:
    Mopac();
    Mopac(const QString &mopac);
    
    Mopac(const Mopac &other);
    
    ~Mopac();
    
    static const char* typeName();
    
    Mopac& operator=(const Mopac &other);
    
    bool operator==(const Mopac &other) const;
    bool operator!=(const Mopac &other) const;
    
    void setExecutable(const QString &mopac_exe);
    void setEnvironment(const QString &variable, const QString &value);
    
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

    void setChargeTemplate(const QString &charge_template);
    
    const QString &chargeTemplate() const;

    void setMopacInputFilename(const QString &filename);
    void setMopacOutputFilename(const QString &filename);
    
    const QString& mopacInputFilename() const;
    const QString& mopacOutputFilename() const;

    AtomCharges calculateCharges(const Molecule &molecule,
                                 const PropertyMap &map) const;

protected:
    double calculateEnergy(const QMPotential::Molecules &molecules,
                           int ntries = 5) const;

    QString energyCommandFile(const QMPotential::Molecules &molecules) const;

    QString forceCommandFile(const QMPotential::Molecules &molecules,
                             const ForceTable &forcetable) const;

    QString fieldCommandFile(const QMPotential::Molecules &molecules,
                             const FieldTable &fieldtable,
                             const SireFF::Probe &probe) const;

    QString potentialCommandFile(const QMPotential::Molecules &molecules,
                                 const PotentialTable &pottable,
                                 const SireFF::Probe &probe) const;

    QString chargeCommandFile(const Molecule &molecule,
                              const PropertyMap &map) const;

private:
    QString createCommandFile(QString cmd_template,
                     const QList< QPair<Vector,SireMol::Element> > &atoms) const;

    QString createCommandFile(QString cmd_template, 
                              const QMPotential::Molecules &molecules) const;  

    QString writeShellFile(const SireBase::TempDir &tempdir) const;

    QStringList runMopac(const QString &cmdfile) const;

    QString _pvt_chargeCommandFile(const Molecule &molecule,
                                   const PropertyMap &map,
                                   QHash<int,int> *atom_idx) const;

    double extractEnergy(const QStringList &mopac_output) const;

    double calculateEnergy(const QString &cmd_file, int ntries) const;

    /** The environmental variables to hold when running Mopac */
    QHash<QString,QString> env_variables;
    
    /** The full path to the mopac executable to run (including
        any necessary command line arguments) */
    QString mopac_exe;
    
    /** The QM method to use to calculate the energy */
    QString qm_method;
    
    /** The template command file used for the energy calculations.
        The QM method and atom coordinates are substituted
        into this template */
    QString energy_template;
    
    /** The template command file used for the force calculations.
        The QM method and atom coordinates are substituted
        into this template */
    QString force_template;
    
    /** The template command file used for the atomic charge
        calculations. The QM method and atom coordinates
        are substituted into this template */
    QString charge_template;
    
    /** The name mopac uses for the input file (e.g. FOR005) */
    QString mopac_input_filename;
    
    /** The name mopac uses for the output file (e.g. FOR006) */
    QString mopac_output_filename;
    
    /** The total charge of the system */
    qint32 total_charge;
};

}

Q_DECLARE_METATYPE( Squire::Mopac )

SIRE_EXPOSE_CLASS( Squire::Mopac )

SIRE_END_HEADER

#endif
