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

#include <cstdlib>

#include "am1bcc.h"

#include "SireMol/moleculedata.h"
#include "SireMol/moleculeinfodata.h"
#include "SireMol/atomname.h"
#include "SireMol/resname.h"
#include "SireMol/groupatomids.h"
#include "SireMol/atomcoords.h"

#include "SireUnits/units.h"
#include "SireUnits/dimensions.h"

#include "SireBase/tempdir.h"
#include "SireBase/sire_process.h"

#include "SireIO/pdb.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace Squire;
using namespace SireMol;
using namespace SireIO;
using namespace SireUnits;
using namespace SireUnits::Dimension;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<AM1BCC> r_am1bcc;

/** Serialise to a binary datastream */
QDataStream SQUIRE_EXPORT &operator<<(QDataStream &ds, const AM1BCC &am1bcc)
{
    writeHeader(ds, r_am1bcc, 1);
    
    SharedDataStream sds(ds);
    
    sds << am1bcc.mopac
        << static_cast<const QMChargeCalculator&>(am1bcc);
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream SQUIRE_EXPORT &operator>>(QDataStream &ds, AM1BCC &am1bcc)
{
    VersionID v = readHeader(ds, r_am1bcc);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> am1bcc.mopac
            >> static_cast<QMChargeCalculator&>(am1bcc);
    }
    else
        throw version_error(v, "1", r_am1bcc, CODELOC);
        
    return ds;
}

/** Constructor */
AM1BCC::AM1BCC() : ConcreteProperty<AM1BCC,QMChargeCalculator>()
{}

/** Copy constructor */
AM1BCC::AM1BCC(const AM1BCC &other)
       : ConcreteProperty<AM1BCC,QMChargeCalculator>(other),
         mopac(other.mopac)
{}

/** Destructor */
AM1BCC::~AM1BCC()
{}

const char* AM1BCC::typeName()
{
    return QMetaType::typeName( qMetaTypeId<AM1BCC>() );
}

/** Comparison operator */
AM1BCC& AM1BCC::operator=(const AM1BCC &other)
{
    if (this != &other)
    {
        mopac = other.mopac;
    
        QMChargeCalculator::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool AM1BCC::operator==(const AM1BCC &other) const
{
    return this == &other or
           (mopac == other.mopac and QMChargeCalculator::operator==(other));
}

/** Comparison operator */
bool AM1BCC::operator!=(const AM1BCC &other) const
{
    return not AM1BCC::operator==(other);
}

/** Set the environmental variable 'variable' to 'value'. It is 
    important that "AMBERHOME" is set, as this is needed by antechamber
    to generate the charges */
void AM1BCC::setEnvironment(const QString &variable, const QString &value)
{
    mopac.setEnvironment(variable, value);
}

/** Return the environmental variables that will be overridden when
    the charges are calculated */
const QHash<QString,QString>& AM1BCC::environment() const
{
    return mopac.environment();
}

/** Return the value of the environmental variable 'variable' */
QString AM1BCC::environment(const QString &variable) const
{
    return mopac.environment(variable);
}
    
/** Set the total charge of the molecule whose charges are being generated */
void AM1BCC::setTotalCharge(int charge)
{
    mopac.setTotalCharge(charge);
    
    if (charge != 0)
        //we cannot scale non-zero total charge molecules
        QMChargeCalculator::setScaleFactor(1.0);
}

/** Return the total charge of the molecule whose charges are being generated */
int AM1BCC::totalCharge() const
{
    return mopac.totalCharge();
}

/** Set the scale factor for the charges */
void AM1BCC::setScaleFactor(double sclfac)
{
    if (sclfac != 1 and this->totalCharge() != 0)
        throw SireError::incompatible_error( QObject::tr(
            "It is not possible to use a scaling factor (%1) with a molecule "
            "with a non-zero total charge (charge == %2)")
                .arg(sclfac).arg(this->totalCharge()), CODELOC );
                
    QMChargeCalculator::setScaleFactor(sclfac);
}

static QString getOutput(const QString &filename)
{
    QFile f(filename);
    
    if (not f.open(QIODevice::ReadOnly))
    {
        return QObject::tr("Unable to open file %1.").arg(filename);
    }
    else
    {
        QTextStream ts(&f);
        
        return ts.readAll();
    }
}

static void runProcess(const QString &cmd, const QString &path,
               const QHash<QString,QString> &environment = QHash<QString,QString>())
{
    QString shellfile = QString("%1/shell_cmd_file").arg(path);

    //write the shell file to run the command
    {
        QFile f( shellfile );
    
        if (not f.open(QIODevice::WriteOnly))
            throw SireError::file_error(f, CODELOC);
    
        QTextStream ts(&f);
    
        //set the environmental variables of the job
        for (QHash<QString,QString>::const_iterator it = environment.constBegin();
             it != environment.constEnd();
             ++it)
        {
            ts << "export " << it.key() << "=\"" << it.value() << "\"\n";
        }

        //set the script to change into the run directory of the job
        ts << QString("\ncd %1").arg(path) << "\n\n";
    
        //set the command to run
        ts << cmd << " > runcommand.log";
    
        f.close();
    }
    
    //now run the command
    {
        Process p = Process::run("sh", shellfile);
        
        //wait until the job has finished
        p.wait();
        
        if (p.wasKilled())
        {
            throw SireError::process_error( QObject::tr(
                "The command {%1} was killed.").arg(cmd), CODELOC );
        }
        
        if (p.isError())
        {
            throw SireError::process_error( QObject::tr(
                "There was an error running the command {%1}. Here is the output "
                "from running the command:\n%2")
                    .arg(cmd, getOutput(QString("%1/runcommand.log").arg(path))), 
                        CODELOC );
        }
    }
}

static void addChargesToAC(const AtomCharges &mulliken_charges,
                           const QString &infilename, const QString &outfilename,
                           const Molecule &molecule)
{
    QFile infile( infilename );
    
    if (not infile.open(QIODevice::ReadOnly))
        throw SireError::file_error(infile, CODELOC);
    
    QFile outfile( outfilename );
    
    if (not outfile.open(QIODevice::WriteOnly))
        throw SireError::file_error(outfile, CODELOC);
        
    QTextStream instream( &infile );
    QTextStream outstream( &outfile );

    outstream.setFieldAlignment( QTextStream::AlignRight );
    outstream.setRealNumberNotation( QTextStream::FixedNotation );
    outstream.setRealNumberPrecision(6);
    
    const MoleculeInfoData &molinfo = molecule.data().info();
    
    while (not instream.atEnd())
    {
        QString line = instream.readLine();
        
        if (line.startsWith("ATOM", Qt::CaseInsensitive))
        {
            // found an antechamber AC ATOM line
            //ATOM      1  C01 MEO     1       0.000   0.000   0.000  0.000000        c3
        
            QStringList words = line.split(" ",QString::SkipEmptyParts);
            
            if (words.count() < 10)
                throw SireError::file_error( QObject::tr(   
                    "The AC file line \"%1\" does not look like a valid antechamber "
                    "format ATOM line.").arg(line), CODELOC );
                    
            AtomName atmnam( words[2] );
            ResName resnam( words[3] );
            
            //get the CGAtomIdx of the atom with this atom and residue name
            CGAtomIdx cgatomidx = molinfo.cgAtomIdx( atmnam + resnam );
            
            //get the mulliken charge of this atom
            Charge atmchg = mulliken_charges[cgatomidx];
            
            //substitute that into the AC file
            outstream << line.left(54);
            
            outstream.setFieldWidth(10);
            outstream << atmchg.to( mod_electron );
            
            outstream.setFieldWidth(0);
            outstream << line.mid(62,-1) << "\n";
        }
        else
        {
            outstream << line << "\n";
        }
    }
}

static AtomCharges extractAM1BCC(const QString &file, const Molecule &molecule,
                                 const double sclfac)
{
    QFile f(file);

    if (not f.open(QIODevice::ReadOnly))
        throw SireError::file_error(f, CODELOC);

    const MoleculeInfoData &molinfo = molecule.data().info();
        
    AtomCharges am1bcc_chgs( molinfo );
        
    QTextStream ts( &f );
    
    while (not ts.atEnd())
    {
        QString line = ts.readLine();

        if (line.startsWith("ATOM"))
        {
            QStringList words = line.split(" ", QString::SkipEmptyParts);
            
            if (words.count() < 10)
            {
                f.close();
                
                throw SireError::file_error( QObject::tr(
                        "The AC Atom line \"%1\" does not look like a valid "
                        "antechamber format ATOM line. Here's the complete "
                        "file;\n%2")
                            .arg(line).arg(::getOutput(file)), CODELOC );
            }
                    
            AtomName atmnam( words[2] );
            ResName resnam( words[3] );
            
            //get the CGAtomIdx of the atom with this atom and residue name
            CGAtomIdx cgatomidx = molinfo.cgAtomIdx( atmnam + resnam );
            
            bool ok;
            
            Charge chg = sclfac * words[8].toDouble(&ok) * mod_electron;

            if (not ok)
            {
                f.close();
                
                throw SireError::file_error( QObject::tr(
                        "The AC Atom line \"%1\" does not look like a valid "
                        "antechamber format ATOM line. Here's the complete "
                        "file;\n%2")
                            .arg(line).arg(::getOutput(file)), CODELOC );
            }
            
            am1bcc_chgs.set( cgatomidx, chg );
        }
    }
    
    return am1bcc_chgs;
}

/** Internal function used to convert AM1 mulliken charges to AM1-BCC charges */
AtomCharges AM1BCC::convertAM1MullikenToAM1BCC(const AtomCharges &mulliken_charges,
                                               const Molecule &molecule,
                                               const PropertyMap &map,
                                               const QString &amberhome) const
{
    //we will use antechamber to create an AC file for this
    //molecule (with zero charges). We will then copy the 
    //AM1 mulliken charges into this file, run "am1bcc" to
    //get the AM1-BCC charges, and will then parse the output
    //file to extract those charges. 

    //antechamber needs the AMBERHOME environmental variable to be set
    QHash<QString,QString> env = this->environment();
    env.insert( "AMBERHOME", amberhome );
    
    //all of this will be performed in a temporary directory!
    TempDir tmpdir;
    //tmpdir.doNotDelete();
    
    const QString pdbfile = QString("%1/molecule.pdb").arg(tmpdir.path());
    const QString acfile = QString("%1/molecule.AC").arg(tmpdir.path());
    const QString am1file = QString("%1/am1mulliken.AC").arg(tmpdir.path());
    const QString am1bccfile = QString("%1/am1bcc.AC").arg(tmpdir.path());
    
    //first we need a PDB of the molecule to input to antechamber
    PDB().write(molecule, pdbfile, map);

    //now we run antechamber on this PDB to create the AC file (needed by am1bcc)
    ::runProcess( QString("%1/bin/antechamber -i molecule.pdb -fi pdb "
                          "-o molecule.AC -j 4 -fo ac -nc %2")
                  .arg(amberhome).arg(mopac.totalCharge()),
            tmpdir.path(), env );
    
    //now we edit the resulting AC file to insert the AM1 mulliken charges...
    ::addChargesToAC( mulliken_charges, acfile, am1file, molecule );
    
    //use "am1bcc" to convert the charges to AM1-BCC charges
    ::runProcess(
            QString("%1/bin/am1bcc -i am1mulliken.AC -o am1bcc.AC -f ac "
                    "-p %1/dat/antechamber/BCCPARM.DAT -j 4")
                    .arg(amberhome),
            tmpdir.path(), env);

    //finally(!) read the output AC file and extract all of the AM1-BCC charges
    return ::extractAM1BCC( am1bccfile, molecule, this->scaleFactor() );
}

/** Return the amber directory (AMBERHOME) */
QString AM1BCC::getAmberHome() const
{
    if (this->environment().contains("AMBERHOME"))
        return this->environment().value("AMBERHOME");
    else
    {
        const char *amberhome = ::getenv("AMBERHOME");
        
        if (not amberhome)
            throw SireError::process_error( QObject::tr(
                    "It is not possible to run the mopac and antechamber "
                    "programs required to generate charges unless the "
                    "AMBERHOME environmental variable is set (this must point "
                    "to the location of the amber directory, e.g. "
                    "$HOME/ambertools/amber10"), CODELOC );
                    
        return QString(amberhome);
    }
}

/** This calculates and returns the AM1BCC charges for the atoms in 
    molecule 'molecule'. This only calculates the charges for the selected
    atoms in the molecule. Either default (0) charges, or the original
    charges are use for atoms that are not selected */
AtomCharges AM1BCC::operator()(const PartialMolecule &molecule,
                               const PropertyMap &map) const
{
    //we use the tools provided by antechamber - all in AMBERHOME
    QString amberhome = this->getAmberHome();

    //we use mopac to calculate all of the AM1 mulliken charges for the entire molecule
    Molecule whole_mol( molecule );
    
    Mopac my_mopac( mopac );
    my_mopac.setExecutable( QString("%1/bin/mopac").arg(amberhome) );
    
    AtomCharges mopac_chgs = my_mopac.calculateCharges(whole_mol, map);

    //we now need to convert these AM1 charges to AM1-BCC charges
    //using the am1bcc program from antechamber
    AtomCharges am1bcc_chgs = convertAM1MullikenToAM1BCC(mopac_chgs, whole_mol, map,
                                                         amberhome);

    if (not molecule.selectedAll())
    {
        //merge any original charges for unselected atoms
        //with the new charges for selected atoms - note that this may
        //(probably will!) make the total charge on the molecule
        //non-integer
        
        const PropertyName &charge_property = map["charge"];
        
        AtomCharges old_charges;
        bool got_charges = false;
        
        if (whole_mol.hasProperty(charge_property))
        {
            const Property &prop = whole_mol.property(charge_property);
            
            if (prop.isA<AtomCharges>())
            {
                old_charges = prop.asA<AtomCharges>();
                got_charges = true;
            }
        }
        
        const MoleculeInfoData &molinfo = whole_mol.data().info();
        
        if (not got_charges)
        {
            old_charges = AtomCharges(molinfo);
        }
        
        if (molecule.selection().nSelected() > molecule.nAtoms()/2)
        {
            for (CGIdx i(0); i < molinfo.nCutGroups(); ++i)
            { 
                for (Index j(0); j < molinfo.nAtoms(i); ++j)
                {
                    CGAtomIdx cgatomidx(i,j);
            
                    if (not molecule.selection().selected(cgatomidx))
                    {
                        am1bcc_chgs.set(cgatomidx, old_charges[cgatomidx]);
                    }
                }
            }
        }
        else
        {
            for (CGIdx i(0); i<molinfo.nCutGroups(); ++i)
            {
                for (Index j(0); j<molinfo.nAtoms(i); ++j)
                {
                    CGAtomIdx cgatomidx(i,j);
            
                    if (molecule.selection().selected(cgatomidx))
                    {
                        old_charges.set(cgatomidx, am1bcc_chgs[cgatomidx]);
                    }
                }
            }
            
            am1bcc_chgs = old_charges;
        }
    }
    
    return am1bcc_chgs;
}

/** Simple internal function to return whether 'val0' and 'val1'
    are comparable */
static bool comparable(const double val0, const double val1)
{
    const double diff = val0 - val1;
    const double tol = 1e-5;

    return ( diff > -tol and diff < tol);
}

/** This returns whether or not the charges will change when going
    from 'oldmol' to 'newmol' - note that this assumes that the 
    charges in 'oldmol' are already AM1BCC charges! If they are
    not, then this will give the wrong answer! */
bool AM1BCC::mayChangeCharges(const PartialMolecule &oldmol,
                              const PartialMolecule &newmol,
                              const PropertyMap &map) const
{
    const MoleculeData &olddata = oldmol.data();
    const MoleculeData &newdata = newmol.data();

    //the charges won't change if the molecule hasn't changed!
    if (olddata.version() == newdata.version())
    {
        return false;
    }
        
    //if atoms have been added or removed, then the charges will change
    if (olddata.info() != newdata.info())
    {
        return true;
    }

    //the charges won't change if the elements or coordinates
    //properties haven't changed
    const PropertyName &coords_property = map["coordinates"];
    const PropertyName &element_property = map["element"];
    
    if (olddata.hasProperty(element_property))
    {
        if (not newdata.hasProperty(element_property))
            return true;
            
        if (olddata.version(element_property) != newdata.version(element_property))
            return true;
    }
    else if (newdata.hasProperty(element_property))
    {
        return true;
    }
    
    if (olddata.hasProperty(coords_property))
    {
        if (not newdata.hasProperty(coords_property))
            return true;
            
        if (olddata.version(coords_property) != newdata.version(coords_property))
        {
            //the coordinates have changed - this will only affect the charges
            //if the conformation of the molecule has changed - we can test
            //this by building rudimentary z-matricies of the two versions 
            //of the molecules and making comparisons that way
            const QVector<Vector> old_coords = olddata.property(coords_property)
                                                      .asA<AtomCoords>().toVector();
                                                  
            const QVector<Vector> new_coords = newdata.property(coords_property)
                                                      .asA<AtomCoords>().toVector();
                                                  
            const int nats = old_coords.count();
            BOOST_ASSERT( new_coords.count() == nats );
            
            if (nats <= 1)
                return false;
            
            const Vector *old_coords_array = old_coords.constData();
            const Vector *new_coords_array = new_coords.constData();
            
            //do atom 2
            if ( not ::comparable(
                        Vector::distance2(old_coords_array[1], old_coords_array[0]),
                        Vector::distance2(new_coords_array[1], new_coords_array[0])) )
            {
                return true;
            }
            else if (nats < 3)
            {
                return false;
            }
            
            //do atom 3
            if ( not (::comparable(
                        Vector::distance2(old_coords_array[2], old_coords_array[1]),
                        Vector::distance2(new_coords_array[2], new_coords_array[1])) and
                 
                      ::comparable(
                        Vector::angle(old_coords_array[2], old_coords_array[1],
                                      old_coords_array[0]),
                        Vector::angle(new_coords_array[2], new_coords_array[1],
                                      new_coords_array[0])) ) )
            {
                return true;
            }
            else if (nats < 4)
            {
                return false;
            }
            
            //now do the remaining atoms
            for (int i=3; i<nats; ++i)
            {
                const Vector &old_v0 = old_coords_array[i];
                const Vector &old_v1 = old_coords_array[i-1];
                const Vector &old_v2 = old_coords_array[i-2];
                const Vector &old_v3 = old_coords_array[i-3];
                
                const Vector &new_v0 = new_coords_array[i];
                const Vector &new_v1 = new_coords_array[i-1];
                const Vector &new_v2 = new_coords_array[i-2];
                const Vector &new_v3 = new_coords_array[i-3];
                
                const double old_bond = Vector::distance2(old_v0, old_v1);
                const double new_bond = Vector::distance2(new_v0, new_v1);
                
                const Angle old_ang = Vector::angle(old_v0, old_v1, old_v2);
                const Angle new_ang = Vector::angle(new_v0, new_v1, new_v2);
                
                const Angle old_dih = Vector::dihedral(old_v0, old_v1, old_v2, old_v3);
                const Angle new_dih = Vector::dihedral(new_v0, new_v1, new_v2, new_v3);
                
                if ( not (::comparable(old_bond, new_bond) and
                          ::comparable(old_ang, new_ang) and
                          ::comparable(old_dih, new_dih)) )
                {
                    return true;
                }
            }
            
            return false;
        }
    }
    else if (newdata.hasProperty(coords_property))
    {
        return true;
    }

    return false;
}
