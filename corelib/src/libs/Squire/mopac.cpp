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

#include "mopac.h"
#include "qmpotential.h"

#include "SireMol/element.h"
#include "SireMol/atomcharges.h"
#include "SireMol/atomelements.h"
#include "SireMol/atomcoords.h"

#include "SireBase/tempdir.h"
#include "SireBase/findexe.h"
#include "SireBase/process.h"

#include "SireUnits/units.h"
#include "SireUnits/dimensions.h"

#include "tostring.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>

using namespace Squire;
using namespace SireFF;
using namespace SireMol;
using namespace SireVol;
using namespace SireBase;
using namespace SireUnits;
using namespace SireUnits::Dimension;
using namespace SireStream;
using namespace SireStream;

static const RegisterMetaType<Mopac> r_mopac;

/** Serialise to a binary datastream */
QDataStream SQUIRE_EXPORT &operator<<(QDataStream &ds, const Mopac &mopac)
{
    writeHeader(ds, r_mopac, 1);
    
    SharedDataStream sds(ds);
    
    sds << mopac.env_variables
        << mopac.mopac_exe << mopac.qm_method
        << mopac.energy_template << mopac.force_template
        << mopac.charge_template
        << mopac.mopac_input_filename << mopac.mopac_output_filename
        << mopac.total_charge;
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream SQUIRE_EXPORT &operator>>(QDataStream &ds, Mopac &mopac)
{
    VersionID v = readHeader(ds, r_mopac);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> mopac.env_variables
            >> mopac.mopac_exe >> mopac.qm_method
            >> mopac.energy_template >> mopac.force_template
            >> mopac.charge_template
            >> mopac.mopac_input_filename >> mopac.mopac_output_filename
            >> mopac.total_charge;
    }
    else
        throw version_error(v, "1", r_mopac, CODELOC);
        
    return ds;
}

static const QString default_energy_template =
       "@QM_METHOD@ 1SCF CHARGE=@QM_CHARGE@ GEO-OK MMOK PRECISE\n"
       "SCF single point calculation. Created by Sire for Mopac 6\n\n"
       "@QM_COORDS@\n";

static const QString default_force_template = "! NEEDS TO BE WRITTEN";

static const QString default_charge_template =
       "@QM_METHOD@ 1SCF CHARGE=@QM_CHARGE@ GEO-OK MMOK PRECISE MULLIK\n"
       "SCF Mulliken point charge calculation. Created by Sire for Mopac 6\n\n"
       "@QM_COORDS@\n";

/** Constructor */
Mopac::Mopac() 
      : ConcreteProperty<Mopac,QMProgram>(),
        qm_method("AM1"),
        energy_template(default_energy_template),
        force_template(default_force_template),
        charge_template(default_charge_template),
        mopac_input_filename("FOR005"),
        mopac_output_filename("FOR006"),
        total_charge(0)
{}

/** Construct, passing in the location of the Mopac executable */
Mopac::Mopac(const QString &mopac)
      : ConcreteProperty<Mopac,QMProgram>(),
        qm_method("AM1"),
        energy_template(default_energy_template),
        force_template(default_force_template),
        charge_template(default_charge_template),
        mopac_input_filename("FOR005"),
        mopac_output_filename("FOR006"),
        total_charge(0)
{
    this->setExecutable(mopac);
}

/** Copy constructor */
Mopac::Mopac(const Mopac &other)
      : ConcreteProperty<Mopac,QMProgram>(other),
        env_variables(other.env_variables), mopac_exe(other.mopac_exe),
        qm_method(other.qm_method),
        energy_template(other.energy_template),
        force_template(other.force_template),
        charge_template(other.charge_template),
        mopac_input_filename(other.mopac_input_filename),
        mopac_output_filename(other.mopac_output_filename),
        total_charge(other.total_charge)
{}

/** Destructor */
Mopac::~Mopac()
{}

/** Copy assignment operator */
Mopac& Mopac::operator=(const Mopac &other)
{
    if (this != &other)
    {
        env_variables = other.env_variables;
        mopac_exe = other.mopac_exe;
        qm_method = other.qm_method;
        energy_template = other.energy_template;
        force_template = other.force_template;
        charge_template = other.charge_template;
        mopac_input_filename = other.mopac_input_filename;
        mopac_output_filename = other.mopac_output_filename;
        total_charge = other.total_charge;
    }
    
    return *this;
}

/** Comparison operator */
bool Mopac::operator==(const Mopac &other) const
{
    return this == &other or 
           (env_variables == other.env_variables and
            mopac_exe == other.mopac_exe and
            qm_method == other.qm_method and
            energy_template == other.energy_template and
            force_template == other.force_template and
            charge_template == other.charge_template and
            mopac_input_filename == other.mopac_input_filename and
            mopac_output_filename == other.mopac_output_filename and
            total_charge == other.total_charge);
}

/** Comparison operator */
bool Mopac::operator!=(const Mopac &other) const
{
    return not this->operator==(other);
}

/** Set the Mopac executable (full path and also arguments) to be used */
void Mopac::setExecutable(const QString &mopac_executable)
{
    mopac_exe = mopac_executable;
}

/** Return the executable (full path and also arguments) to be used. This
    is null if the executable is searched for in the path */
QString Mopac::executable() const
{
    return mopac_exe;
}

/** Set the environmental variable 'variable' to have the value 'value'
    when the Mopac executable is run. This replaces any existing
    value of this environmental variable */
void Mopac::setEnvironment(const QString &variable, const QString &value)
{
    env_variables[variable] = value;
}

/** Return all of the environmental variables that are to be set explicitly
    when Mopac is run. This does not include any environmental variables
    that have not been explicitly set, but do have values */
const QHash<QString,QString>& Mopac::environment() const
{
    return env_variables;
}

/** Return the value of the explicitly set environmental variable 'variable'.
    A null string is returned if this variable has not been set 
    explicitly (this does not mean the variable doesn't exist - merely
    that a specific value has not been set) */
QString Mopac::environment(const QString &variable) const
{
    return env_variables.value(variable);
}

/** Set the QM method to be used to calculate the energy or
    force (e.g. AM1, PM3). This will substitute for 
    @QM_METHOD@ in the energy and force command file templates */
void Mopac::setMethod(const QString &method)
{
    qm_method = method;
}

/** Return the QM method to be used to calculate the energy or
    force (e.g. AM1, PM3). This will substitute for 
    @QM_METHOD@ in the energy and force command file templates */
const QString& Mopac::method() const
{
    return qm_method;
}

/** Set the total charge of the system (in unit charges) */
void Mopac::setTotalCharge(int charge)
{
    total_charge = charge;
}

/** Tell this interface which name Mopac uses for the input file. This
    is hardcoded into Mopac and depends on the system on which this
    program is running, e.g. 'FOR005' is the name for the Mopac
    I am using on OS X */
void Mopac::setMopacInputFilename(const QString &f)
{
    mopac_input_filename = f;
}

/** Tell this interface which name Mopac uses for the output file. This
    is hardcoded into Mopac and depends on the system on which this
    program is running, e.g. 'FOR006' is the name for the Mopac
    I am using on OS X */
void Mopac::setMopacOutputFilename(const QString &f)
{
    mopac_output_filename = f;
}

/** Return the name of the file mopac is using for input. This
    is hardcoded into Mopac and depends on the system on which this
    program is running, e.g. 'FOR005' is the name for the Mopac
    I am using on OS X */
const QString& Mopac::mopacInputFilename() const
{
    return mopac_input_filename;
}

/** Return the name of the file mopac is using for output. This
    is hardcoded into Mopac and depends on the system on which this
    program is running, e.g. 'FOR006' is the name for the Mopac
    I am using on OS X */
const QString& Mopac::mopacOutputFilename() const
{
    return mopac_output_filename;
}

/** Return the total charge of the system */
int Mopac::totalCharge() const
{
    return total_charge;
}

/** Set the template for the command file to be used to get
    Mopac to calculate an energy. The following tags will
    be substituted in the template;
    
    @QM_METHOD@          - the desired QM method (e.g. AM1 or PM3)
    @QM_COORDS@          - the list of elements and coordinates of QM atoms
    @QM_CHARGE@          - the total charge of the system
*/
void Mopac::setEnergyTemplate(const QString &cmd_template)
{
    energy_template = cmd_template;
}

/** Return the template for the command file to be used to get Mopac
    to calculate the energy. */
const QString& Mopac::energyTemplate() const
{
    return energy_template;
}

/** Set the template for the command file to be used to get
    Mopac to calculate the forces. The following tags will
    be substituted in the template;
    
    @QM_METHOD@          - the desired QM method (e.g. AM1 or PM3)
    @QM_COORDS@          - the list of elements and coordinates of QM atoms
    @QM_CHARGE@          - the total charge of the system
*/
void Mopac::setForceTemplate(const QString &cmd_template)
{
    force_template = cmd_template;
}

/** Return the template for the command file to be used to get Mopac
    to calculate the forces. */
const QString& Mopac::forceTemplate() const
{
    return force_template;
}

/** Set the template for the command file to be used to get
    Mopac to calculate atomic partial charges. The following tags will
    be substituted in the template;
    
    @QM_METHOD@          - the desired QM method (e.g. AM1 or PM3)
    @QM_COORDS@          - the list of elements and coordinates of QM atoms
    @QM_CHARGE@          - the total charge of the system
*/
void Mopac::setChargeTemplate(const QString &cmd_template)
{
    charge_template = cmd_template;
}

/** Return the template for the command file to be used to get Mopac
    to calculate atomic partial charges. */
const QString& Mopac::chargeTemplate() const
{
    return force_template;
}

static void writeMopacLine(const Element &element,
                           double x, double y, double z,
                           QTextStream &ts)
{
    ts.setFieldWidth(5);
    
    ts << element.symbol();

    ts.setFieldWidth(12);
    ts << x;

    ts.setFieldWidth(0);
    ts << "  0  ";

    ts.setFieldWidth(12);
    ts << y;

    ts.setFieldWidth(0);
    ts << "  0  ";

    ts.setFieldWidth(12);
    ts << z;

    ts.setFieldWidth(0);
    ts << "  0\n";
}

/** Function used to substitute in the atom coordinates
    into the provided mopac command template */
QString Mopac::createCommandFile(QString cmd_template,
                                 const QList< QPair<Vector,Element> > &atoms) const
{
    if (atoms.isEmpty())
        return QString::null;

    cmd_template.replace( QLatin1String("@QM_METHOD@"),
                          qm_method, Qt::CaseInsensitive );

    cmd_template.replace( QLatin1String("@QM_CHARGE@"),
                          QString::number(total_charge), Qt::CaseInsensitive );
                  
    QString coords_string;
    QTextStream ts(&coords_string);

    ts.setFieldAlignment(QTextStream::AlignRight);
    ts.setRealNumberNotation(QTextStream::FixedNotation);
    ts.setRealNumberPrecision(4);
    
    //if there are less than 4 atoms, then we need to write a z-matrix.
    // Otherwise, we write out cartesian coordinates
    if (atoms.count() < 4)
    {
        //first atom is at the origin
        ::writeMopacLine( atoms.at(0).second, 0, 0, 0, ts );
        
        if (atoms.count() > 1)
        {
            //the second atom is along the x-axis, bonded to the first
            ::writeMopacLine( atoms.at(1).second, 
            
                              Length(Vector::distance(atoms.at(1).first, 
                                                      atoms.at(0).first)).to(angstrom),
                                                      
                              0, 0, ts );

            if (atoms.count() > 2)
            {
                //the third atom is bonded to the second, angled to the first
                ::writeMopacLine( atoms.at(2).second,
                
                        Length( Vector::distance(atoms.at(2).first, 
                                                 atoms.at(1).first) ).to(angstrom),
                                                 
                        Vector::angle(atoms.at(2).first, atoms.at(1).first,
                                      atoms.at(0).first).to(degrees),
                        
                        0, ts );
            }
        }
    }
    else
    {
        for (QList< QPair<Vector,Element> >::const_iterator it = atoms.constBegin();
             it != atoms.constEnd();
             ++it)
        {
            const Vector &v = it->first;
        
            ::writeMopacLine( it->second,
                              Length(v.x()).to(angstrom),
                              Length(v.y()).to(angstrom),
                              Length(v.z()).to(angstrom), ts );
        }
    }

    cmd_template.replace( QLatin1String("@QM_COORDS@"),
                          coords_string, Qt::CaseInsensitive );
                                       
    return cmd_template;
}

/** Function used to substitute in the atom coordinates
    into the provided mopac command template */
QString Mopac::createCommandFile(QString cmd_template,
                                 const QMPotential::Molecules &molecules) const
{
    int nmols = molecules.count();
    const ChunkedVector<QMPotential::Molecule> &molecules_array 
                                                    = molecules.moleculesByIndex();

    //first, collect together all of the coordinates and elements...
    QList< QPair<Vector,Element> > atoms;
    
    for (int i=0; i<nmols; ++i)
    {
        const QMPotential::Molecule &molecule = molecules_array[i];

        //loop through the atoms...
        const CoordGroupArray &coords = molecule.coordinates();
        const PackedArray2D<Element> &elements = molecule.parameters().atomicParameters();
        
        int natoms = coords.nCoords();
        
        BOOST_ASSERT( natoms == elements.nValues() );
        
        const Vector *coords_array = coords.constCoordsData();
        const Element *elements_array = elements.constValueData();
        
        for (int j=0; j<natoms; ++j)
        {
            const Element &element = elements_array[j];
        
            if (element.nProtons() > 0)
            {
                //this is not a dummy atom!
                const Vector &c = coords_array[j];
                
                atoms.append( QPair<Vector,Element>(c, element) );
            }
        }
    }
    
    return this->createCommandFile(cmd_template, atoms);
}

/** Return the command file that will be used to calculate the energy of the 
    molecules in 'molecules' */
QString Mopac::energyCommandFile(const QMPotential::Molecules &molecules) const
{
    return createCommandFile(energy_template, molecules);
}

/** Return the command files that will be used to calculate the forces on the  
    atoms of the molecules in 'molecules' */
QString Mopac::forceCommandFile(const QMPotential::Molecules &molecules,
                                const ForceTable&) const
{
    return createCommandFile(force_template, molecules);
}

/** Return the command files that will be used to calculate the fields on the  
    atoms of the molecules in 'molecules' */
QString Mopac::fieldCommandFile(const QMPotential::Molecules &molecules,
                                const FieldTable &fieldtable,
                                const SireFF::Probe &probe) const
{
    throw SireError::unsupported( QObject::tr(
            "Calculating QM fields using the Mopac() interface is not yet "
            "supported."), CODELOC );
            
    return QString::null;
}

/** Return the command files that will be used to calculate the potentials on the  
    atoms of the molecules in 'molecules' */
QString Mopac::potentialCommandFile(const QMPotential::Molecules &molecules,
                                    const PotentialTable &pottable,
                                    const SireFF::Probe &probe) const
{
    throw SireError::unsupported( QObject::tr(
            "Calculating QM potentials using the Mopac() interface is not yet "
            "supported."), CODELOC );
            
    return QString::null;
}

/** Extract the energy from the mopac output in 'mopac_output' */
double Mopac::extractEnergy(const QStringList &lines) const
{
    // Looking for the line;
    // FINAL HEAT OF FORMATION =        -59.22261 KCAL
    QRegExp regexp("FINAL HEAT OF FORMATION\\s*=\\s*([-\\d\\.]+)\\s* KCAL");

    foreach (const QString &line, lines)
    {
        if (regexp.indexIn(line) != -1)
        {
            //we've found the FINAL HEAT OF FORMATION line
            QString num = regexp.cap(1);
        
            bool ok;
        
            double nrg = num.toDouble(&ok);
        
            if (not ok)
                throw SireError::process_error( QObject::tr(
                    "The energy obtained from Mopac is garbled (%1) - %2.")
                        .arg(regexp.cap(1), regexp.cap(0)), CODELOC );
        
            //the energy is already in kcal per mol
            return nrg * kcal_per_mol;
        }
    }

    //the FINAL HEAT OF FORMATION line was not found!
    throw SireError::process_error( QObject::tr(
            "Could not find the total energy in the mopac output!\n"
                "%1").arg(lines.join("\n")), CODELOC );
}

/** Internal function used to write the shell script that is used to
    run the mopac job and collect the output
*/
QString Mopac::writeShellFile(const TempDir &tempdir) const
{
    QString cmdfile = QString("%1/run_mopac.cmd").arg(tempdir.path());
    
    QFile f(cmdfile);
    
    if (not f.open(QIODevice::WriteOnly))
        throw SireError::file_error(f, CODELOC);
    
    QTextStream ts(&f);

    //set the environmental variables of the job
    for (QHash<QString,QString>::const_iterator it = env_variables.constBegin();
         it != env_variables.constEnd();
         ++it)
    {
        ts << "export " << it.key() << "=\"" << it.value() << "\"\n";
    }

    //set the script to change into the run directory of the job
    ts << QString("\ncd %1").arg(tempdir.path()) << "\n\n";

    //write the line used to run mopac
    if (mopac_exe.isEmpty())
    {
        //the user hasn't specified a mopac executable - try to find one
        QString found_mopac = SireBase::findExe("mopac").absoluteFilePath();
        ts << found_mopac << " > " << mopac_output_filename << "\n";
    }
    else
        ts << mopac_exe << " > " << mopac_output_filename << "\n";

    //this 'sync' costs 100ms to perform, which is quite long compared
    //to the mopac calculation, but may be necessary on some filesystems...
    //ts << "sync\n";
    
    f.close();
    
    return cmdfile;
}

static QByteArray readAll(const QString &file)
{
    QFile f(file);
    
    if (f.open(QIODevice::ReadOnly))
    {
        return f.readAll();
    }
    else
        return QByteArray();
}

QStringList Mopac::runMopac(const QString &cmdfile) const
{
    //create a temporary directory in which to run Mopac
    QString tmppath = env_variables.value("TMPDIR");
    
    if (tmppath.isEmpty())
        tmppath = QDir::temp().absolutePath();

    TempDir tmpdir(tmppath);
    //tmpdir.doNotDelete();

    //write the file processed by the shell used to run the job
    QString shellfile = this->writeShellFile(tmpdir);

    {
        QFile f( QString("%1/%2").arg(tmpdir.path(), mopac_input_filename) );
        
        if (not f.open( QIODevice::WriteOnly ))
            throw SireError::file_error(f, CODELOC);
   
        //write the command file
        f.write( cmdfile.toUtf8() );
        f.close();
    }

    //run the shell file...
    Process p = Process::run( "sh", shellfile );

    //wait until the job has finished
    p.wait();
    
    if (p.wasKilled())
    {
        throw SireError::process_error( QObject::tr(
            "The Mopac job was killed."), CODELOC );
    }
    
    if (p.isError())
    {
        QByteArray shellcontents = ::readAll(shellfile);
        QByteArray cmdcontents = ::readAll(QString("%1/%2")
                                              .arg(tmpdir.path(), mopac_input_filename));

        QByteArray outputcontents = ::readAll(QString("%1/%2")
                                              .arg(tmpdir.path(), mopac_output_filename));

        throw SireError::process_error( QObject::tr(
            "There was an error running the Mopac job.\n"
            "The shell script used to run the job was;\n"
            "*****************************************\n"
            "%1\n"
            "*****************************************\n"
            "The mopac input used to run the job was;\n"
            "*****************************************\n"
            "%2\n"
            "*****************************************\n"
            "The mopac output was;\n"
            "*****************************************\n"
            "%3\n"
            "*****************************************\n"
            )
                .arg( QLatin1String(shellcontents),
                      QLatin1String(cmdcontents),
                      QLatin1String(outputcontents) ), CODELOC );
    }

    //read all of the output
    QFile f( QString("%1/%2").arg(tmpdir.path(), mopac_output_filename) );
    
    if ( not (f.exists() and f.open(QIODevice::ReadOnly)) )
    {
        QByteArray shellcontents = ::readAll(shellfile);
        QByteArray cmdcontents = ::readAll(QString("%1/%2")
                                            .arg(tmpdir.path(), mopac_input_filename));
    
        throw SireError::process_error( QObject::tr(
            "There was an error running the Mopac job - no output was created.\n"
            "The shell script used to run the job was;\n"
            "*****************************************\n"
            "%1\n"
            "*****************************************\n"
            "The mopac input used to run the job was;\n"
            "*****************************************\n"
            "%2\n"
            "*****************************************\n")
                .arg( QLatin1String(shellcontents),
                      QLatin1String(cmdcontents) ), CODELOC );
    }

    QStringList lines;
    
    QTextStream ts(&f);
    
    while (not ts.atEnd())
    {
        lines.append( ts.readLine() );
    }
    
    return lines;
}

/** Return the energy calculate according to the Mopac command
    file 'cmd_file' (this is the contents of the file, not
    the path to the file) */
double Mopac::calculateEnergy(const QString &cmdfile, int ntries) const
{
    QStringList lines = this->runMopac(cmdfile);

    try
    {
        //parse the output to get the energy
        return this->extractEnergy(lines);
    }
    catch(...)
    {
        qDebug() << "Mopac process error. Number of remaining attempts = " << ntries;

        //print out the last twenty lines of output
        const int nlines_to_print = 20;

        qDebug() << "Printing out the last" << nlines_to_print << "lines of output...";

        int start = qMax( 0, lines.count() - nlines_to_print );

        for (int i=start; i < lines.count(); ++i)
        {
            qDebug() << lines.at(i);
        }

        if (ntries <= 0)
            //don't bother trying again - it's not going to work!
            throw;
            
        //give it one more go - you never know, it may work
        return this->calculateEnergy(cmdfile, ntries-1);
    }
}

/** Run Mopac and use it to calculate the energy of the molecules in 
    'molecules'. This blocks until Mopac has completed */
double Mopac::calculateEnergy(const QMPotential::Molecules &molecules,
                               int ntries) const
{
    //create the command file to be used by Mopac
    QString cmdfile = this->energyCommandFile(molecules);
    
    return this->calculateEnergy(cmdfile, ntries);
}

QString Mopac::_pvt_chargeCommandFile(const Molecule &molecule,
                                      const PropertyMap &map,
                                      QHash<int,int> *atom_idxs) const
{
    const AtomElements &elements = molecule.property( map["element"] )
                                           .asA<AtomElements>();
                                           
    const AtomCoords &coords = molecule.property( map["coordinates"] )
                                       .asA<AtomCoords>();
                                       
    QList< QPair<Vector,Element> > atoms;
    
    int nats = coords.nAtoms();
    
    const Vector *coords_array = coords.array().constCoordsData();
    const Element *elements_array = elements.array().constValueData();
    
    if (atom_idxs)
    {
        atom_idxs->clear();
        atom_idxs->reserve(nats);
    
        int j = 1;
    
        for (int i=0; i<nats; ++i)
        {
            if (elements_array[i].nProtons() > 0)
            {
                atoms.append( QPair<Vector,Element>(coords_array[i], elements_array[i]) );
                atom_idxs->insert( j, i );
                ++j;
            }
        }
    }
    else
    {
        for (int i=0; i<nats; ++i)
        {
            if (elements_array[i].nProtons() > 0)
                atoms.append( QPair<Vector,Element>(coords_array[i], elements_array[i]) );
        }
    }
    
    return createCommandFile(charge_template, atoms);
}

/** Return the command file used to calculate partial atomic charges */
QString Mopac::chargeCommandFile(const Molecule &molecule, const PropertyMap &map) const
{
    return this->_pvt_chargeCommandFile(molecule, map, 0);
}

/** Use Mopac to calculate the partial charges for the passed molecule
    (using the total charge set in this program) */    
AtomCharges Mopac::calculateCharges(const Molecule &molecule,
                                    const PropertyMap &map) const
{
    QHash<int,int> atom_idxs;
    
    QString cmdfile = this->_pvt_chargeCommandFile(molecule, map, &atom_idxs);

    const int ntries = 5;
    
    int itry = 0;
    
    while (true)
    {
        try
        {
            const QStringList lines = this->runMopac(cmdfile);

            //we are looking for the line
            //          MULLIKEN POPULATIONS AND CHARGES
            int idx = lines.indexOf( 
                                QRegExp("\\s*MULLIKEN POPULATIONS AND CHARGES\\s*") );
    
            if (idx == -1)
                throw SireError::process_error( QObject::tr(
                        "Could not find the partial charges in the mopac output!\n"
                        "%1").arg(lines.join("\n")), CODELOC );

            const MoleculeInfoData &molinfo = molecule.data().info();

            AtomCharges chgs(molinfo);
    
            for (int i=idx+1; i<lines.count(); ++i)
            {
                const QStringList words = lines[i].split(" ",QString::SkipEmptyParts);
            
                if (words.count() < 3)
                    //we've finished
                    break;
                    
                bool ok;
                    
                int iatm = words[0].toInt(&ok);

                if (not ok)
                    throw SireError::process_error( QObject::tr(
                            "Could not find the partial charge from the line "
                            "%1 (full output is;%2)")
                                .arg(lines[i], lines.join("\n")), CODELOC );
                    
                chgs.set( molinfo.cgAtomIdx( AtomIdx(atom_idxs[iatm]) ),
                          words[2].toDouble(&ok) * mod_electron );
                          
                if (not ok)
                    throw SireError::process_error( QObject::tr(
                            "Could not find the partial charge from the line "
                            "%1 (full output is;%2)")
                                .arg(lines[i], lines.join("\n")), CODELOC );
            }

            return chgs;
        }
        catch(...)
        {
            ++itry;
        
            qDebug() << "PROBLEM RUNNING MOPAC - TRY" << itry << "of" << ntries;
        
            if (itry >= ntries)
                throw;
        }
    }
    
    return AtomCharges();
}
                                    
const char* Mopac::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Mopac>() );
}
