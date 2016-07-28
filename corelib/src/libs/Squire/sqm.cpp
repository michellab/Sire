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

#include "SireMaths/vector.h"

#include "sqm.h"
#include "qmpotential.h"
#include "latticecharges.h"

#include "SireMol/element.h"

#include "SireMM/cljprobe.h"

#include "SireFF/potentialtable.h"

#include "SireVol/grid.h"

#include "SireBase/tempdir.h"
#include "SireBase/findexe.h"
#include "SireBase/sire_process.h"

#include "SireUnits/units.h"

#include "tostring.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QTime>
#include <QDebug>

using namespace Squire;
using namespace SireMM;
using namespace SireFF;
using namespace SireMol;
using namespace SireVol;
using namespace SireBase;
using namespace SireUnits;
using namespace SireUnits::Dimension;
using namespace SireStream;
using namespace SireStream;

static const RegisterMetaType<SQM> r_sqm;

/** Serialise to a binary datastream */
QDataStream SQUIRE_EXPORT &operator<<(QDataStream &ds, const SQM &sqm)
{
    writeHeader(ds, r_sqm, 1);
    
    SharedDataStream sds(ds);
    
    sds << sqm.env_variables
        << sqm.sqm_exe << sqm.qm_method
        << sqm.energy_template << sqm.force_template
        << sqm.total_charge
        << sqm.max_sqm_runtime
        << sqm.max_sqm_lines
        << sqm.expected_n_qm;
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream SQUIRE_EXPORT &operator>>(QDataStream &ds, SQM &sqm)
{
    VersionID v = readHeader(ds, r_sqm);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> sqm.env_variables
            >> sqm.sqm_exe >> sqm.qm_method
            >> sqm.energy_template >> sqm.force_template
            >> sqm.total_charge >> sqm.max_sqm_runtime
            >> sqm.max_sqm_lines >> sqm.expected_n_qm;
    }
    else
        throw version_error(v, "1", r_sqm, CODELOC);
        
    return ds;
}

//here is the template for the energy command file. We try to save as many lines
//as possible as SQM has a 1000 line limit for the command file!!!
static const QString default_energy_template =
"single point energy calculation created for QM or QM/MM calcs using Sire\n"
"&qmmm\n"
" qm_theory = '@QM_METHOD@', qmcharge = @QM_CHARGE@, maxcyc = 0, qmmm_int = @USE_LATTICE_POINTS@,\n"
" /\n"  // need space before backslash or it is not recognised
"@QM_COORDS@\n"
"@LATTICE_POINTS@\n";

//the number of lines in the energy template command file, excluding the lines needed
//for the QM atoms or MM atoms
static const int n_energy_template_lines = 4;

static const QString default_force_template = "! NEEDS TO BE WRITTEN";

/** Constructor */
SQM::SQM() 
       : ConcreteProperty<SQM,QMProgram>(),
         qm_method("AM1"),
         energy_template(default_energy_template),
         force_template(default_force_template),
         total_charge(0),
         max_sqm_runtime( 5 * 60 * 1000 ),
         max_sqm_lines(999),
         expected_n_qm(50)
{}

/** Construct, passing in the location of the SQM executable */
SQM::SQM(const QString &sqm_exe)
       : ConcreteProperty<SQM,QMProgram>(),
         qm_method("AM1"),
         energy_template(default_energy_template),
         force_template(default_force_template),
         total_charge(0),
         max_sqm_runtime( 5 * 60 * 1000 ),
         max_sqm_lines(999),
         expected_n_qm(50)
{
    this->setExecutable(sqm_exe);
}

/** Copy constructor */
SQM::SQM(const SQM &other)
       : ConcreteProperty<SQM,QMProgram>(other),
         env_variables(other.env_variables), sqm_exe(other.sqm_exe),
         qm_method(other.qm_method),
         energy_template(other.energy_template),
         force_template(other.force_template),
         total_charge(other.total_charge),
         max_sqm_runtime(other.max_sqm_runtime),
         max_sqm_lines(other.max_sqm_lines),
         expected_n_qm(other.expected_n_qm)
{}

/** Destructor */
SQM::~SQM()
{}

/** Copy assignment operator */
SQM& SQM::operator=(const SQM &other)
{
    if (this != &other)
    {
        env_variables = other.env_variables;
        sqm_exe = other.sqm_exe;
        qm_method = other.qm_method;
        energy_template = other.energy_template;
        force_template = other.force_template;
        total_charge = other.total_charge;
        max_sqm_runtime = other.max_sqm_runtime;
        max_sqm_lines = other.max_sqm_lines;
        expected_n_qm = other.expected_n_qm;
    }
    
    return *this;
}

/** Comparison operator */
bool SQM::operator==(const SQM &other) const
{
    return this == &other or 
           (env_variables == other.env_variables and
            sqm_exe == other.sqm_exe and
            qm_method == other.qm_method and
            energy_template == other.energy_template and
            force_template == other.force_template and
            total_charge == other.total_charge and
            max_sqm_runtime == other.max_sqm_runtime and
            max_sqm_lines == other.max_sqm_lines and
            expected_n_qm == other.expected_n_qm);
}

/** Comparison operator */
bool SQM::operator!=(const SQM &other) const
{
    return not this->operator==(other);
}

/** Set the SQM executable (full path and also arguments) to be used */
void SQM::setExecutable(const QString &sqm_executable)
{
    sqm_exe = sqm_executable;
}

/** Return the executable (full path and also arguments) to be used. This
    is null if the executable is searched for in the path */
QString SQM::executable() const
{
    return sqm_exe;
}

/** Set the environmental variable 'variable' to have the value 'value'
    when the SQM executable is run. This replaces any existing
    value of this environmental variable */
void SQM::setEnvironment(const QString &variable, const QString &value)
{
    env_variables[variable] = value;
}

/** Return all of the environmental variables that are to be set explicitly
    when SQM is run. This does not include any environmental variables
    that have not been explicitly set, but do have values */
const QHash<QString,QString>& SQM::environment() const
{
    return env_variables;
}

/** Return the value of the explicitly set environmental variable 'variable'.
    A null string is returned if this variable has not been set 
    explicitly (this does not mean the variable doesn't exist - merely
    that a specific value has not been set) */
QString SQM::environment(const QString &variable) const
{
    return env_variables.value(variable);
}

/** Set the QM method to be used to calculate the energy or
    force (e.g. AM1, PM3, AM1/d etc. See the AmberTools documentation
    for SQM to find the supported methods and the string used to 
    specify that method). This will substitute for
    @QM_METHOD@ in the command file templates, and should be the same
    string used in SQM as specified in the SQM documentation */
void SQM::setMethod(const QString &method)
{
    qm_method = method;
}

/** Return the QM method to be used to calculate the energy or
    force (e.g. AM1, PM3, AM1/d etc.). This will substitute for
    @QM_METHOD@ in the command file templates */
const QString& SQM::method() const
{
    return qm_method;
}

QString SQM::toString() const
{
    return QObject::tr("SQM( method = %1 )").arg(method());
}

/** Set the total charge of the system (in unit charges) */
void SQM::setTotalCharge(int charge)
{
    total_charge = charge;
}

/** Return the total charge of the system */
int SQM::totalCharge() const
{
    return total_charge;
}

/** Set the maximum allowed runtime for the SQM job - this is used
    to detect hangs - if the SQM job takes longer than this 
    time then it is killed and an exception raised. The maximum
    runtime is measured in milliseconds */
void SQM::setMaximumRunTime(int ms)
{
    if (ms > 0)
        max_sqm_runtime = ms;
    else
        max_sqm_runtime = 0;
}

/** Return the maximum runtime allowed for a SQM job, in milliseconds */
int SQM::maximumRunTime() const
{
    return max_sqm_runtime;
}

/** Set the maximum number of lines that can be parsed from an SQM input file.
    Currently, SQM has a hard-coded limit of 1000 lines! */
void SQM::setMaximumNumberOfSQMInputLines(int numlines)
{
    if (numlines < 0)
        max_sqm_lines = -1;
    else
        max_sqm_lines = numlines;
}

/** Set the maximum number of expected QM atoms. This is used, together with
    the maximum number of lines in a SQM input file, to work out the maximum
    number of supported MM atoms */
void SQM::setExpectedNumberOfQMAtoms(int natoms)
{
    if (natoms < 0)
        expected_n_qm = -1;
    else
        expected_n_qm = natoms;
}

/** Return the maximum number of supported SQM input lines. This returns
    -1 if SQM doesn't have a file size limit */
int SQM::maximumNumberOfSQMInputLines() const
{
    return max_sqm_lines;
}

/** Return the maximum number of expected QM atoms. This returns -1 if
    we don't expect any QM atoms */
int SQM::expectedNumberOfQMAtoms() const
{
    return expected_n_qm;
}

/** Return the maximum number of MM atoms supported by SQM. This returns
    -1 if there is no limit on the number of atoms */
int SQM::numberOfMMAtomsLimit() const
{
    if (max_sqm_lines < 0 or expected_n_qm < 0)
        return -1;
    else
        return max_sqm_lines - n_energy_template_lines - expected_n_qm - 2;
        //need -2 as have two extra lines if we have MM atoms, #EXCHARGES and #END
}

/** Return the maximum number of MM atoms supported by SQM if there
    are 'num_qm_atoms' QM atoms. This returns
    -1 if there is no limit on the number of atoms */
int SQM::numberOfMMAtomsLimit(int num_qm_atoms) const
{
    if (max_sqm_lines < 0 or num_qm_atoms < 0)
        return -1;
    else
        return max_sqm_lines - n_energy_template_lines - num_qm_atoms - 2;
        //need -2 as have two extra lines if we have MM atoms, #EXCHARGES and #END
}

/** Set the template for the command file to be used to get
    SQM to calculate an energy. The following tags will
    be substituted in the template;
    
    @QM_METHOD@          - the desired QM method (e.g. AM1)
    @QM_COORDS@          - the list of elements and coordinates of QM atoms
    @QM_CHARGE@          - the total charge of the system
    @LATTICE_POINTS@     - the coordinates and charges of the lattice points (MM atoms)
    @USE_LATTICE_POINTS@ - whether or not lattice points are used (0 or 1)
*/
void SQM::setEnergyTemplate(const QString &cmd_template)
{
    energy_template = cmd_template;
}

/** Return the template for the command file to be used to get SQM
    to calculate the energy. */
const QString& SQM::energyTemplate() const
{
    return energy_template;
}

/** Set the template for the command file to be used to get
    SQM to calculate the forces. The following tags will
    be substituted in the template;
    
    @QM_METHOD@          - the desired QM method (e.g. AM1)
    @QM_COORDS@          - the list of elements and coordinates of QM atoms
    @QM_CHARGE@          - the total charge of the system
    @LATTICE_POINTS@     - the coordinates and charges of the lattice points (MM atoms)
    @USE_LATTICE_POINTS@ - whether or not lattice points are used (0 or 1)
*/
void SQM::setForceTemplate(const QString &cmd_template)
{
    force_template = cmd_template;
}

/** Return the template for the command file to be used to get SQM
    to calculate the forces. */
const QString& SQM::forceTemplate() const
{
    return force_template;
}

/** Function used to substitute in the atom and lattice coordinates
    into the provided SQM command template */
QString SQM::createCommandFile(QString cmd_template,
                               const QMPotential::Molecules &molecules,
                               const LatticeCharges &lattice_charges) const
{
    //replace the easy things...
    cmd_template.replace( QLatin1String("@QM_METHOD@"),
                          qm_method, Qt::CaseInsensitive );

    cmd_template.replace( QLatin1String("@QM_CHARGE@"),
                          QString::number(total_charge), Qt::CaseInsensitive );

    if (lattice_charges.isEmpty())
    {
        //there are no lattice charges
        cmd_template.replace( QLatin1String("@USE_LATTICE_POINTS@"),
                              QString::number(0), Qt::CaseInsensitive );
                          
        cmd_template.replace( QLatin1String("@LATTICE_POINTS@"),
                              QLatin1String(""), Qt::CaseInsensitive );
    }
                  
    //now build the list of all of the atoms
    QStringList atom_coords;
    
    int nmols = molecules.count();
    const ChunkedVector<QMPotential::Molecule> &molecules_array 
                                                    = molecules.moleculesByIndex();
    
    // format for atom lines;
    // atomic_number  atom_name  x_coords  y_coords  z_coords  (free format)
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
    
                atom_coords.append( QString("%1  %2  %3  %4  %5")
                                        .arg(QString::number(element.nProtons()),
                                             element.symbol(),
                                             QString::number(c.x(), 'f', 8),
                                             QString::number(c.y(), 'f', 8),
                                             QString::number(c.z(), 'f', 8) ) );
            }
        }
    }
    
    if (max_sqm_lines > 0 and
        atom_coords.count() + n_energy_template_lines > max_sqm_lines)
    {
        throw SireError::unsupported( QObject::tr(
                "SQM has a %1-line limit for its input command file, so it is not possible "
                "to have more than %2 QM atoms. You currently have %3 atoms...")
                    .arg(max_sqm_lines)
                    .arg(max_sqm_lines - n_energy_template_lines)
                    .arg(atom_coords.count()), CODELOC );
    }
    
    cmd_template.replace( QLatin1String("@QM_COORDS@"),
                          atom_coords.join("\n"), Qt::CaseInsensitive );
    
    //put the lattice charges in now (as they can make the command
    //file *very* large)
    if (not lattice_charges.isEmpty())
    {
        QStringList charges;
        
        int ncharges = 0;
        const LatticeCharge *charges_array = lattice_charges.constData();
        
        for (int i=0; i<lattice_charges.count(); ++i)
        {
            const LatticeCharge &charge = charges_array[i];
        
            if (charge.charge() != 0)
            {
                ncharges += 1;
            
                Element element = charge.element();
                
                charges.append( QString("%1 %2 %3 %4 %5 %6")
                                    .arg(QString::number(element.nProtons()),
                                         element.symbol(),
                                         QString::number(charge.x(), 'f', 8),
                                         QString::number(charge.y(), 'f', 8),
                                         QString::number(charge.z(), 'f', 8),
                                         QString::number(charge.charge(), 'f', 8) ) );
            }
        }

        if (max_sqm_lines > 0 and
            atom_coords.count()+n_energy_template_lines+ncharges+2 > max_sqm_lines)
        {
            throw SireError::unsupported( QObject::tr(
                    "SQM has a %1-line limit for its input command file, so it is not possible "
                    "to have more than %2 MM atoms when you have %3 QM atoms. You currently have "
                    "%4 MM atoms. Try to reduce the cutoff so that you have fewer MM atoms.")
                        .arg(max_sqm_lines)
                        .arg(max_sqm_lines - n_energy_template_lines - atom_coords.count() - 2)
                        .arg(atom_coords.count())
                        .arg(ncharges), CODELOC );
        }

        cmd_template.replace( QLatin1String("@USE_LATTICE_POINTS@"),
                              QString::number(1), Qt::CaseInsensitive );
            
        cmd_template.replace( QLatin1String("@LATTICE_POINTS@"),
                              QString("#EXCHARGES\n%1\n#END").arg(charges.join("\n")),
                              Qt::CaseInsensitive );
    }
                                       
    return cmd_template;
}

/** Return the command file that will be used to calculate the energy of the 
    molecules in 'molecules' */
QString SQM::energyCommandFile(const QMPotential::Molecules &molecules) const
{
    return this->createCommandFile(energy_template, molecules);
}

/** Return the command file that will be used to calculate the energy of the 
    molecules in 'molecules' in the field of point charges in 'lattice_charges' */
QString SQM::energyCommandFile(const QMPotential::Molecules &molecules,
                               const LatticeCharges &lattice_charges) const
{
    return this->createCommandFile(energy_template, molecules, lattice_charges);
}

/** Return the command files that will be used to calculate the forces on the  
    atoms of the molecules in 'molecules' */
QString SQM::forceCommandFile(const QMPotential::Molecules &molecules,
                              const ForceTable&) const
{
    throw SireError::unsupported( QObject::tr(
            "Calculating QM forces using the SQM() interface is "
            "currently not supported."), CODELOC );
            
    return QString::null;
}

/** Return the command files that will be used to calculate the forces on the  
    atoms of the molecules in 'molecules' in the field of point charges
    in 'lattice_charges' - this also calculates the forces on those
    point charges */
QString SQM::forceCommandFile(const QMPotential::Molecules &molecules,
                              const LatticeCharges &lattice_charges,
                              const ForceTable&) const
{
    throw SireError::unsupported( QObject::tr(
            "Calculating QM forces using the SQM() interface is "
            "currently not supported."), CODELOC );
            
    return QString::null;
}

/** Return the command files that will be used to calculate the fields on the  
    atoms of the molecules in 'molecules' */
QString SQM::fieldCommandFile(const QMPotential::Molecules &molecules,
                              const FieldTable &fieldtable,
                              const SireFF::Probe &probe) const
{
    throw SireError::unsupported( QObject::tr(
            "Calculating QM fields using the SQM() interface is "
            "currently not supported."), CODELOC );
            
    return QString::null;
}

/** Return the command files that will be used to calculate the fields on the  
    atoms of the molecules in 'molecules' */
QString SQM::fieldCommandFile(const QMPotential::Molecules &molecules,
                              const LatticeCharges &lattice_charges,
                              const FieldTable &fieldtable,
                              const SireFF::Probe &probe) const
{
    throw SireError::unsupported( QObject::tr(
            "Calculating QM fields using the SQM() interface is "
            "currently not supported."), CODELOC );
            
    return QString::null;
}

/** Return the command files that will be used to calculate the potentials on the  
    atoms of the molecules in 'molecules' */
QString SQM::potentialCommandFile(const QMPotential::Molecules &molecules,
                                  const PotentialTable &pottable,
                                  const SireFF::Probe &probe) const
{
    throw SireError::unsupported( QObject::tr(
            "Calculating QM potentials using the SQM() interface is "
            "currently not supported."), CODELOC );
            
    return QString::null;
}

/** Return the command files that will be used to calculate the potentials on the  
    atoms of the molecules in 'molecules' */
QString SQM::potentialCommandFile(const QMPotential::Molecules &molecules,
                                     const LatticeCharges &lattice_charges,
                                     const PotentialTable &pottable,
                                     const SireFF::Probe &probe) const
{
    throw SireError::unsupported( QObject::tr(
            "Calculating QM potentials using the SQM() interface is "
            "currently not supported."), CODELOC );
            
    return QString::null;
}

/** Extract the energy from the SQM output in 'sqm_output' */
double SQM::extractEnergy(QFile &sqm_output) const
{
    QTextStream ts(&sqm_output);

    //the energy is written on the line that looks like this;
    //Total SCF energy    =      -8037.75735154 kcal/mol  (     -348.54331345 eV)
    QRegExp regexp("Total\\s+SCF\\s+energy\\s*=\\s*([-\\d\\.]+)", Qt::CaseInsensitive);

    QStringList lines;
    
    while (not ts.atEnd())
    {
        QString line = ts.readLine();
        lines.append(line);

        if (regexp.indexIn(line) != -1)
        {
            //we've found the "Total SCF energy" line
            QString num = regexp.cap(1);
        
            bool ok;
        
            double nrg = num.toDouble(&ok);
        
            if (not ok)
                throw SireError::process_error( QObject::tr(
                    "The energy obtained from SQM is garbled (%1) - %2.")
                        .arg(regexp.cap(1), regexp.cap(0)), CODELOC );
        
            //the energy is already in kcal mol-1
            return nrg;
        }
    }

    //ok - something went wrong - we couldn't find the output
    // Lets try reading the file again (in case it wasn't fully written
    // to the disc)
    QTextStream ts2(&sqm_output);

    lines = QStringList();
    
    while (not ts2.atEnd())
    {
        QString line = ts2.readLine();
        lines.append(line);

        if (regexp.indexIn(line) != -1)
        {
            //we've found the "SCF energy" line
            QString num = regexp.cap(1);
        
            bool ok;
        
            double nrg = num.toDouble(&ok);
        
            if (not ok)
                throw SireError::process_error( QObject::tr(
                    "The energy obtained from SQM is garbled (%1) - %2.")
                        .arg(regexp.cap(1), regexp.cap(0)), CODELOC );
        
            //the energy is already in kcal mol-1
            return nrg;
        }
    }

    //the "SCF energy" line was not found!
    throw SireError::process_error( QObject::tr(
            "Could not find the total energy in the SQM output!\n"
                "%1").arg(lines.join("\n")), CODELOC );
}

/** Internal function used to write the shell script that is used to
    run the SQM job and collect the output
*/
QString SQM::writeShellFile(const TempDir &tempdir) const
{
    QString cmdfile = QString("%1/run_sqm.cmd").arg(tempdir.path());
    
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

    //write the line used to run SQM
    if (sqm_exe.isEmpty())
    {
        //the user hasn't specified a SQM executable - try to find one
        QString found_sqm = SireBase::findExe("sqm").absoluteFilePath();
        ts << QString("%1 -O -i sqm_input -o sqm_output\n")
                    .arg(found_sqm);
    }
    else
        ts << QString("%1 -O -i sqm_input -o sqm_output\n")
                        .arg(sqm_exe);

    ts << "sync\n";
    
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

/** Return the energy calculate according to the SQM command
    file 'cmd_file' (this is the contents of the file, not
    the path to the file) */
double SQM::calculateEnergy(const QString &cmdfile, int ntries) const
{
    //create a temporary directory in which to run SQM
    QString tmppath = env_variables.value("TMPDIR");
    
    if (tmppath.isEmpty())
        tmppath = QDir::temp().absolutePath();

    TempDir tmpdir(tmppath);
    //tmpdir.doNotDelete();

    //write the file processed by the shell used to run the job
    QString shellfile = this->writeShellFile(tmpdir);

    {
        QFile f( QString("%1/sqm_input").arg(tmpdir.path()) );
        
        if (not f.open( QIODevice::WriteOnly ))
            throw SireError::file_error(f, CODELOC);
   
        //write the command file
        f.write( cmdfile.toUtf8() );
        f.close();
    }

    //run the shell file...
    //QTime t;
    //qDebug() << "Running SQM...";
    //t.start();
    Process p = Process::run( "sh", shellfile );

    //wait until the job has finished
    if (not p.wait(max_sqm_runtime))
    {
        qDebug() << "Maximum SQM runtime was exceeded - has it hung?";
        p.kill();
        
        if (ntries > 0)
            return this->calculateEnergy(cmdfile, ntries-1);
    }
    
    //int ms = t.elapsed();
    //qDebug() << "SQM finised. Took" << ms << "ms";
    
    if (p.wasKilled())
    {
        throw SireError::process_error( QObject::tr(
            "The SQM job was killed."), CODELOC );
    }
    
    if (p.isError())
    {
        QByteArray shellcontents = ::readAll(shellfile);
        QByteArray cmdcontents = ::readAll(QString("%1/sqm_input").arg(tmpdir.path()));
        QByteArray outputcontents = ::readAll(QString("%1/sqm_output")
                                                    .arg(tmpdir.path()));

        throw SireError::process_error( QObject::tr(
            "There was an error running the SQM job.\n"
            "The shell script used to run the job was;\n"
            "*****************************************\n"
            "%1\n"
            "*****************************************\n"
            "The SQM input used to run the job was;\n"
            "*****************************************\n"
            "%2\n"
            "*****************************************\n"
            "The SQM output was;\n"
            "*****************************************\n"
            "%3\n"
            "*****************************************\n"
            )
                .arg( QLatin1String(shellcontents),
                      QLatin1String(cmdcontents),
                      QLatin1String(outputcontents) ), CODELOC );
    }

    //read all of the output
    QFile f( QString("%1/sqm_output").arg(tmpdir.path()) );
    
    if ( not (f.exists() and f.open(QIODevice::ReadOnly)) )
    {
        QByteArray shellcontents = ::readAll(shellfile);
        QByteArray cmdcontents = ::readAll(QString("%1/sqm_input").arg(tmpdir.path()));
    
        throw SireError::process_error( QObject::tr(
            "There was an error running the SQM job - no output was created.\n"
            "The shell script used to run the job was;\n"
            "*****************************************\n"
            "%1\n"
            "*****************************************\n"
            "The SQM input used to run the job was;\n"
            "*****************************************\n"
            "%2\n"
            "*****************************************\n")
                .arg( QLatin1String(shellcontents),
                      QLatin1String(cmdcontents) ), CODELOC );
    }


    try
    {
        //parse the output to get the energy
        return this->extractEnergy(f);
    }
    catch(...)
    {
        qDebug() << "SQM process error. Number of remaining attempts = " << ntries;

        //print out the last 150 lines of output
        const int nlines_to_print = 150;

        qDebug() << "Printing out the last" << nlines_to_print << "lines of output...";

        QFile f( QString("%1/sqm_output").arg(tmpdir.path()) );

        if ( not (f.exists() and f.open(QIODevice::ReadOnly)) )    
            qDebug() << "Could not read the file" << tmpdir.path() << "/sqm_output";
        else
        {
            QStringList lines;
 
            QTextStream ts(&f);

            while (not ts.atEnd())
            {
                lines.append( ts.readLine() );

                if (lines.count() > nlines_to_print)
                    lines.removeFirst();
            }

            foreach (QString line, lines){ qDebug() << qPrintable(line); }
        }

        if (ntries <= 0)
            //don't bother trying again - it's not going to work!
            throw;
            
        //give it one more go - you never know, it may work
        return this->calculateEnergy(cmdfile, ntries-1);
    }
}

/** Run SQM and use it to calculate the energy of the molecules in 
    'molecules'. This blocks until SQM has completed */
double SQM::calculateEnergy(const QMPotential::Molecules &molecules,
                            int ntries) const
{
    if (molecules.count() == 0)
        return 0;

    //create the command file to be used by SQM
    QString cmdfile = this->energyCommandFile(molecules);
    
    return this->calculateEnergy(cmdfile, ntries);
}

/** Calculate the SQM QM energy of the molecules in 'molecules'
    in the field of point charges in 'lattice_charges' */
double SQM::calculateEnergy(const QMPotential::Molecules &molecules,
                            const LatticeCharges &lattice_charges,
                            int ntries) const
{
    if (molecules.count() == 0)
        return 0;

    //create the command file to be used by SQM
    QString cmdfile = this->energyCommandFile(molecules, lattice_charges);
    
    return this->calculateEnergy(cmdfile, ntries);
}

/** Internal function to calculate the potentials specified in the 
    passed command file, and to return then together with the points
    at which they were evaluated */
QHash<QString,double> SQM::calculatePotential(const QString &cmdfile,
                                              int ntries) const
{
    throw SireError::unsupported( QObject::tr(
            "Calculating QM potentials using the SQM() interface is "
            "currently not supported."), CODELOC );
            
    return QHash<QString,double>();
}

/** Calculate the potential around the passed molecules, and place them into the passed
    potential table, optionally scaled by 'scale_potential' */
void SQM::calculatePotential(const QMPotential::Molecules &molecules,
                             PotentialTable &pottable,
                             const SireFF::Probe &probe,
                             double scale_potential, int ntries) const
{
    throw SireError::unsupported( QObject::tr(
            "Calculating QM potentials using the SQM() interface is "
            "currently not supported."), CODELOC );
}

/** Calculate the potentials around the passed molecules, and place them into the passed
    potential table, optionally scaled by 'scale_potential', and return the accompanying
    potentials on the passed lattice points, also scaled by 'scale_potential' */
QVector<MolarEnergy> SQM::calculatePotential(const QMPotential::Molecules &molecules,
                                             const LatticeCharges &lattice_charges,
                                             PotentialTable &pottable,
                                             const SireFF::Probe &probe,
                                             double scale_potential, int ntries) const
{
    throw SireError::unsupported( QObject::tr(
            "Calculating QM potentials using the SQM() interface is "
            "currently not supported."), CODELOC );
            
    return QVector<MolarEnergy>();
}

const char* SQM::typeName()
{
    return QMetaType::typeName( qMetaTypeId<SQM>() );
}
