/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2007  Christopher Woods
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

#include "SireMM/cljnbpairs.h"

#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QTextStream>
#include <QByteArray>

#include "protoms.h"

#include "pdb.h"

#include "SireMol/molecule.h"
#include "SireMol/molecules.h"
#include "SireMol/moleditor.h"
#include "SireMol/reseditor.h"
#include "SireMol/atomeditor.h"
#include "SireMol/atomname.h"
#include "SireMol/resname.h"
#include "SireMol/resnum.h"
#include "SireMol/groupatomids.h"
#include "SireMol/geometryperturbation.h"
#include "SireMol/chargeperturbation.h"
#include "SireMol/connectivity.h"
#include "SireMol/mover.hpp"
#include "SireMol/selector.hpp"
#include "SireMol/cgeditor.h"

#include "SireMol/atomcharges.h"

#include "SireMove/zmatrix.h"

#include "SireCAS/trigfuncs.h"

#include "SireMM/ljparameter.h"
#include "SireMM/atomljs.h"
#include "SireMM/internalff.h"
#include "SireMM/ljperturbation.h"
#include "SireMM/internalperturbation.h"

#include "SireUnits/units.h"

#include "SireBase/tempdir.h"
#include "SireBase/findexe.h"
#include "SireBase/sire_process.h"

#include "SireMove/errors.h"
#include "SireMol/errors.h"
#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireIO;
using namespace SireMol;
using namespace SireMM;
using namespace SireMove;
using namespace SireID;
using namespace SireCAS;
using namespace SireUnits;
using namespace SireBase;
using namespace SireStream;

///////////
/////////// Implementation of ProtoMSParameters
///////////

ProtoMSParameters::ProtoMSParameters()
{}

ProtoMSParameters::~ProtoMSParameters()
{}

PropertyName ProtoMSParameters::charge_property( "charge" );
PropertyName ProtoMSParameters::lj_property( "LJ" );

PropertyName ProtoMSParameters::initial_charge_property( "initial_charge" );
PropertyName ProtoMSParameters::initial_lj_property( "initial_LJ" );

PropertyName ProtoMSParameters::final_charge_property( "final_charge" );
PropertyName ProtoMSParameters::final_lj_property( "final_LJ" );

PropertyName ProtoMSParameters::connectivity_property( "connectivity" );
PropertyName ProtoMSParameters::coords_property( "coordinates" );

PropertyName ProtoMSParameters::bond_property( "bond" );
PropertyName ProtoMSParameters::angle_property( "angle" );
PropertyName ProtoMSParameters::dihedral_property( "dihedral" );
PropertyName ProtoMSParameters::ub_property( "Urey-Bradley" );

PropertyName ProtoMSParameters::zmatrix_property( "z-matrix" );
PropertyName ProtoMSParameters::nb_property( "intrascale" );

PropertyName ProtoMSParameters::perts_property( "perturbations" );

///////////
/////////// Implementation of ProtoMS
///////////

static const RegisterMetaType<ProtoMS> r_protoms(NO_ROOT);

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const ProtoMS &protoms)
{
    writeHeader(ds, r_protoms, 1);

    SharedDataStream sds(ds);

    sds << protoms.paramfiles << protoms.protoms_exe;

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, ProtoMS &protoms)
{
    VersionID v = readHeader(ds, r_protoms);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> protoms.paramfiles >> protoms.protoms_exe;
    }
    else
        throw version_error( v, "1", r_protoms, CODELOC );

    return ds;
}

ProtoMSParameters ProtoMS::protoms_parameters;

/** Constructor */
ProtoMS::ProtoMS()
{}

/** Constructor, specifying the location of ProtoMS */
ProtoMS::ProtoMS(const QString &protoms)
{
    this->setExecutable(protoms);
}

/** Destructor */
ProtoMS::~ProtoMS()
{}

/** Add a parameter file to the list of files which will be used
    to parameterise the molecules */
void ProtoMS::addParameterFile(const QString &paramfile)
{
    //get the absolute path to this file
    QFile f(paramfile);

    if (not f.open( QIODevice::ReadOnly))
    {
        throw SireError::file_error(f, CODELOC);
    }

    QString absfilepath = QFileInfo(f).absoluteFilePath();

    f.close();

    if (not paramfiles.contains(absfilepath))
        paramfiles.append(absfilepath);
}

/** Return the list of parameter files which will be used to
    parameterise the molecules, in the order that they will
    be read */
QStringList ProtoMS::parameterFiles() const
{
    return paramfiles;
}

/** Set the path to the ProtoMS executable that will be used
    to parameterise the molecules */
void ProtoMS::setExecutable(const QString &protoms)
{
    protoms_exe = protoms;
}

/** Return the command file used to run ProtoMS on the passed molecule as the passed type */
QString ProtoMS::parameterisationCommandFile(const Molecule &molecule,
                                             int type) const
{
    QString contents;

    QTextStream ts(&contents, QIODevice::WriteOnly);

    QString name = molecule.name();

    if (name.isEmpty())
        name = "molecule";

    for (int i=0; i<paramfiles.count(); ++i)
    {
        ts << QString("parfile%1 %2\n").arg(i+1).arg(paramfiles[i]);
    }

    switch (type)
    {
    case PROTEIN:
        ts << "protein1 " << name << "\n";
        break;
    case SOLUTE:
        ts << "solute1 " << name << "\n";
        break;
    case SOLVENT:
        ts << "solvent1 " << name << "\n";
        break;

    default:
        throw SireError::program_bug( QObject::tr(
                "Unrecognised ProtoMS molecule type (%1)").arg(type), CODELOC );
    }

    ts << "streamheader off\n"
          "streaminfo stdout\n"
          "streamwarning stdout\n"
          "streamfatal stdout\n"
          "streamparameters stdout\n\n"
          "chunk1 printparameters\n";

    return contents;
}

/** Internal function used to write the ProtoMS command file */
QString ProtoMS::writeCommandFile(const TempDir &tmpdir,
                                  const Molecule &molecule, int type) const
{
    //write a PDB of the molecule to the TMPDIR for ProtoMS to read
    // ProtoMS solutes can be named by setting the name of the
    // PDB file
    QString name = molecule.name();

    if (name.isEmpty())
        name = "molecule";

    {
        QFile f( QString("%1/%2").arg(tmpdir.path(),name) );
        f.open(QIODevice::WriteOnly);

        f.write( QString("header %1\n").arg(name).toUtf8().constData() );

        PDB().write( molecule, f );

        f.close();
    }

    QString cmdfile = QString("%1/protoms_input").arg(tmpdir.path());

    QFile f(cmdfile);
    f.open(QIODevice::WriteOnly);

    QTextStream ts(&f);

    ts << this->parameterisationCommandFile(molecule, type);

    return cmdfile;
}

/** Internal function used to write a shell file used to run ProtoMS */
QString ProtoMS::writeShellFile(const TempDir &tempdir,
                                const QString &cmdfile) const
{
    QString shellfile = QString("%1/run_protoms.cmd").arg(tempdir.path());

    QFile f(shellfile);
    f.open(QIODevice::WriteOnly);

    QTextStream ts(&f);

    //set the script to change into the run directory of the job
    ts << QString("\ncd %1").arg(tempdir.path()) << "\n\n";

    //write the line used to run ProtoMS
    if (protoms_exe.isEmpty())
    {
        //the user hasn't specified a ProtoMS executable - try to find one
        QString found_protoms = SireBase::findExe("protoms2").absoluteFilePath();
        ts << QString("%1 %2 > protoms_output\n")
                    .arg(found_protoms, cmdfile);
    }
    else
        ts << QString("%1 %2 > protoms_output\n")
                        .arg(protoms_exe, cmdfile);

    f.close();

    return shellfile;
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

namespace SireIO
{
namespace detail
{

class ProtoMSWorkspace
{
public:
    ProtoMSWorkspace()
    {}

    ~ProtoMSWorkspace()
    {}

    QHash< QString, QHash<QString,AtomIdx> > solute_atom_cache;
    QHash< ResNum, QHash<QString,AtomIdx> > protein_atom_cache;
};

} // end of namespace detail
} // end of namespace SireIO

using namespace SireIO::detail;

static AtomIdx getSoluteAtom(const Molecule &molecule,
                             const QString &atomname, const QString &resname,
                             ProtoMSWorkspace &workspace)
{
    if (workspace.solute_atom_cache[resname].contains(atomname))
    {
        return workspace.solute_atom_cache[resname][atomname];
    }

    const MoleculeInfoData &molinfo = molecule.data().info();

    QList<AtomIdx> atomidxs = molinfo.getAtomsIn( ResName(resname,CaseInsensitive) );

    QString lower_atomname = atomname.toLower();

    foreach (AtomIdx atomidx, atomidxs)
    {
        if (molinfo.name(atomidx).value().toLower() == lower_atomname)
        {
            workspace.solute_atom_cache[resname].insert(atomname,atomidx);
            return atomidx;
        }
    }

    throw SireMol::missing_atom( QObject::tr(
            "There is no atom in residue %1 with name %2.")
                .arg(resname).arg(atomname), CODELOC );

    return AtomIdx();
}

static AtomIdx getProteinAtom(const Molecule &molecule,
                              const QString &atomname, const QString &resnum,
                              ProtoMSWorkspace &workspace)
{
    bool ok;

    ResNum rnum( resnum.toInt(&ok) );

    if (not ok)
        throw SireError::program_bug( QObject::tr(
            "The string '%1' should contain a residue number!").arg(resnum),
                CODELOC );

    if (workspace.protein_atom_cache[rnum].contains(atomname))
    {
        return workspace.protein_atom_cache[rnum][atomname];
    }

    const MoleculeInfoData &molinfo = molecule.data().info();

    QList<AtomIdx> atomidxs = molinfo.getAtomsIn(rnum);

    QString lower_atomname = atomname.toLower();

    foreach (AtomIdx atomidx, atomidxs)
    {
        if (molinfo.name(atomidx).value().toLower() == lower_atomname)
        {
            workspace.protein_atom_cache[rnum].insert(atomname,atomidx);
            return atomidx;
        }
    }

    throw SireMol::missing_atom( QObject::tr(
            "There is no atom in residue %1 with name %2.")
                .arg(resnum).arg(atomname), CODELOC );

    return AtomIdx();
}

static AtomEditor getSolventAtom(MolEditor &molecule, const QString &atomname,
                                 ProtoMSWorkspace &workspace)
{
    return molecule.atom( AtomName(atomname, CaseInsensitive) );
}

static AtomEditor getSoluteAtom(MolEditor &molecule,
                                const QString &atomname, const QString &resname,
                                 ProtoMSWorkspace &workspace)
{
    AtomIdx atomidx = ::getSoluteAtom( static_cast<const Molecule&>(molecule),
                                       atomname, resname, workspace );

    return molecule.atom(atomidx);
}

static AtomEditor getProteinAtom(MolEditor &molecule,
                                 const QString &atomname, const QString &resnum,
                                 ProtoMSWorkspace &workspace)
{
    bool ok;

    ResNum rnum( resnum.toInt(&ok) );

    if (not ok)
        throw SireError::program_bug( QObject::tr(
            "The string '%1' should contain a residue number!").arg(resnum),
                CODELOC );

    AtomIdx atomidx = ::getProteinAtom( static_cast<const Molecule&>(molecule),
                                        atomname, resnum, workspace );

    return molecule.atom(atomidx);
}

/** This processes the output line that contains the z-matrix */
void ProtoMS::processZMatrixLine(const QStringList &words, const Molecule &mol,
                                 int type, ZMatrix &zmatrix,
                                 ProtoMSWorkspace &workspace) const
{
    AtomIdx atom, bond, angle, dihedral;

    if (type == SOLVENT)
        return;

    else if (type == SOLUTE)
    {
        atom = getSoluteAtom(mol, words[2], words[4], workspace);
        bond = getSoluteAtom(mol, words[7], words[9], workspace);
        angle = getSoluteAtom(mol, words[12], words[14], workspace);
        dihedral = getSoluteAtom(mol, words[17], words[19], workspace);
    }
    else if (type == PROTEIN)
    {
        atom = getProteinAtom(mol, words[2], words[4], workspace);
        bond = getProteinAtom(mol, words[7], words[9], workspace);
        angle = getProteinAtom(mol, words[12], words[14], workspace);
        dihedral = getProteinAtom(mol, words[17], words[19], workspace);
    }

    zmatrix.add( atom, bond, angle, dihedral );
}

void ProtoMS::processZMatrixPertLine(const QStringList &words,
                                     const Molecule &mol, int type,
                                     QList<SireMol::GeomPertPtr> &geom_perturbations,
                                     const ZMatrix &zmatrix,
                                     const PropertyMap &pert_map,
                                     ProtoMSWorkspace &workspace) const
{
    AtomIdx atom;

    if (type == SOLVENT)
        return;

    else if (type == SOLUTE)
    {
        atom = getSoluteAtom(mol, words[2], words[4], workspace);
    }
    else if (type == PROTEIN)
    {
        atom = getProteinAtom(mol, words[2], words[4], workspace);
    }

    const ZMatrixLine &line = zmatrix[atom];

    bool ok1, ok2;
    double initial = words[7].toDouble(&ok1);
    double final = words[8].toDouble(&ok2);

    if (not (ok1 and ok2))
        throw SireError::io_error( QObject::tr(
            "Could not interpret the range of a variable perturbation "
            "from the line \"%1\".")
                .arg(words.join(" ")), CODELOC );

    if ( words[6] == "BOND" )
    {
        geom_perturbations.append( BondPerturbation( line.atom(), line.bond(),
                                                     initial*angstrom,
                                                     final*angstrom,
                                                     pert_map ) );
    }
    else if ( words[6] == "ANGLE" )
    {
        geom_perturbations.append( AnglePerturbation( line.atom(), line.bond(),
                                                      line.angle(),
                                                      initial*degrees,
                                                      final*degrees,
                                                      pert_map ) );
    }
    else if ( words[6] == "DIHEDRAL" )
    {
        geom_perturbations.append( DihedralPerturbation( line.atom(), line.bond(),
                                                         line.angle(), line.dihedral(),
                                                         initial*degrees,
                                                         final*degrees,
                                                         pert_map ) );
    }
    else
    {
        throw SireError::io_error( QObject::tr(
            "Could not interpret the variable type from the line \"%1\".")
                .arg(words.join(" ")), CODELOC );
    }
}

void ProtoMS::processAtomPertLine(const QStringList &words,
                                  MolEditor &mol, int type,
                                  const QString &initial_charge_property,
                                  const QString &final_charge_property,
                                  const QString &initial_lj_property,
                                  const QString &final_lj_property,
                                  ProtoMSWorkspace &workspace) const
{
    AtomEditor atom;

    if (type == SOLUTE)
        atom = getSoluteAtom(mol, words[2], words[4], workspace);

    else if (type == SOLVENT)
        atom = getSolventAtom(mol, words[2], workspace);

    else if (type == PROTEIN)
        atom = getProteinAtom(mol, words[2], words[4], workspace);

    atom.setProperty( initial_charge_property,
                      words[7].toDouble() * mod_electron );

    atom.setProperty( final_charge_property,
                      words[8].toDouble() * mod_electron );

    atom.setProperty( initial_lj_property,
                      LJParameter( words[10].toDouble() * angstrom,
                                   words[13].toDouble() * kcal_per_mol ) );

    atom.setProperty( final_lj_property,
                      LJParameter( words[11].toDouble() * angstrom,
                                   words[14].toDouble() * kcal_per_mol ) );

    mol = atom.molecule();
}

void ProtoMS::processBondDeltaLine(const QStringList &words, const Molecule &mol,
                                   int type, ZMatrix &zmatrix,
                                   ProtoMSWorkspace &workspace) const
{
    AtomIdx atom, bond;

    if (type == SOLVENT)
        return;

    else if (type == SOLUTE)
    {
        atom = getSoluteAtom(mol, words[2], words[4], workspace);
        bond = getSoluteAtom(mol, words[7], words[9], workspace);
    }
    else if (type == PROTEIN)
    {
        atom = getProteinAtom(mol, words[2], words[4], workspace);
        bond = getProteinAtom(mol, words[7], words[9], workspace);
    }

    try
    {
        zmatrix.setBondDelta( atom, bond,
                              words[12].toDouble() * angstrom );
    }
    catch(const SireMove::zmatrix_error&)
    {
        //skip z-matrix errors...
        return;
    }
}

void ProtoMS::processAngleDeltaLine(const QStringList &words, const Molecule &mol,
                                    int type, ZMatrix &zmatrix,
                                    ProtoMSWorkspace &workspace) const
{
    AtomIdx atom, bond, angle;

    if (type == SOLVENT)
        return;

    else if (type == SOLUTE)
    {
        atom = getSoluteAtom(mol, words[2], words[4], workspace);
        bond = getSoluteAtom(mol, words[7], words[9], workspace);
        angle = getSoluteAtom(mol, words[12], words[14], workspace);
    }
    else if (type == PROTEIN)
    {
        atom = getProteinAtom(mol, words[2], words[4], workspace);
        bond = getProteinAtom(mol, words[7], words[9], workspace);
        angle = getProteinAtom(mol, words[12], words[14], workspace);
    }

    try
    {
        zmatrix.setAngleDelta( atom, bond, angle,
                               words[17].toDouble() * degrees );
    }
    catch(const SireMove::zmatrix_error&)
    {
        //skip zmatrix errors
        return;
    }
}

void ProtoMS::processDihedralDeltaLine(const QStringList &words, const Molecule &mol,
                                       int type, ZMatrix &zmatrix,
                                       ProtoMSWorkspace &workspace) const
{
    AtomIdx atom, bond, angle, dihedral;

    if (type == SOLVENT)
        return;

    else if (type == SOLUTE)
    {
        atom = getSoluteAtom(mol, words[2], words[4], workspace);
        bond = getSoluteAtom(mol, words[7], words[9], workspace);
        angle = getSoluteAtom(mol, words[12], words[14], workspace);
        dihedral = getSoluteAtom(mol, words[17], words[19], workspace);
    }
    else if (type == PROTEIN)
    {
        atom = getProteinAtom(mol, words[2], words[4], workspace);
        bond = getProteinAtom(mol, words[7], words[9], workspace);
        angle = getProteinAtom(mol, words[12], words[14], workspace);
        dihedral = getProteinAtom(mol, words[17], words[19], workspace);
    }

    try
    {
        zmatrix.setDihedralDelta( atom, bond,
                                  angle, dihedral,
                                  words[22].toDouble() * degrees );
    }
    catch(const SireMove::zmatrix_error &e)
    {
        qDebug() << e.toString();

        //ignore zmatrix errors
        return;
    }
}

/** This processes the output line that contains the atom parameters */
void ProtoMS::processAtomLine(const QStringList &words, MolEditor &editmol,
                              int type, const QString &charge_property,
                              const QString &lj_property,
                              ProtoMSWorkspace &workspace) const
{
    AtomEditor atom;

    if (type == SOLUTE)
        atom = getSoluteAtom(editmol, words[2], words[4], workspace);

    else if (type == SOLVENT)
        atom = getSolventAtom(editmol, words[2], workspace);

    else if (type == PROTEIN)
        atom = getProteinAtom(editmol, words[2], words[4], workspace);

    SireUnits::Dimension::Charge chg = words[7].toDouble() * mod_electron;
    LJParameter lj( words[9].toDouble() * angstrom,
                    words[11].toDouble() * kcal_per_mol );

    atom.setProperty( charge_property, chg );
    atom.setProperty( lj_property, lj );

    editmol = atom.molecule();
}

/** This processes the output line that contains the bond parameters */
void ProtoMS::processBondLine(const QStringList &words, const Molecule &mol,
                              int type,
                              TwoAtomFunctions &bondfuncs,
                              ProtoMSWorkspace &workspace) const
{
    if (type == SOLVENT)
        return;

    if (words.count() < 15)
        throw SireError::io_error( QObject::tr(
            "Cannot understand the ProtoMS bond line\n%1")
                .arg(words.join(" ")), CODELOC );

    Symbol r = InternalPotential::symbols().bond().r();

    double k = words[12].toDouble();
    double r0 = words[14].toDouble();

    Expression bondfunc = k * SireMaths::pow_2( r - r0 );

    AtomIdx atom0, atom1;

    if (type == SOLUTE)
    {
        atom0 = getSoluteAtom(mol, words[2], words[4], workspace);
        atom1 = getSoluteAtom(mol, words[7], words[9], workspace);
    }
    else if (type == PROTEIN)
    {
        atom0 = getProteinAtom(mol, words[2], words[4], workspace);
        atom1 = getProteinAtom(mol, words[7], words[9], workspace);
    }

    bondfuncs.set( atom0, atom1, bondfunc );
}

/** This processes the output line that contains the bond perturbation parameters */
PerturbationPtr ProtoMS::processBondPertLine(const QStringList &words,
                                             const Molecule &mol,
                                             int type,
                                             const PropertyName &bond_property,
                                             ProtoMSWorkspace &workspace) const
{
    if (type == SOLVENT)
        return PerturbationPtr();

    if (words.count() < 19)
        throw SireError::io_error( QObject::tr(
            "Cannot understand the ProtoMS bond perturbation line\n%1")
                .arg(words.join(" ")), CODELOC );

    double kb = words[12].toDouble();
    double r0b = words[14].toDouble();

    double kf = words[16].toDouble();
    double r0f = words[18].toDouble();

    AtomIdx atom0, atom1;

    if (type == SOLUTE)
    {
        atom0 = getSoluteAtom(mol, words[2], words[4], workspace);
        atom1 = getSoluteAtom(mol, words[7], words[9], workspace);
    }
    else if (type == PROTEIN)
    {
        atom0 = getProteinAtom(mol, words[2], words[4], workspace);
        atom1 = getProteinAtom(mol, words[7], words[9], workspace);
    }

    if (kb != kf or r0b != r0f)
    {
        Symbol r = InternalPotential::symbols().bond().r();
        Symbol k("k");
        Symbol r0("r0");

        Expression base = k * SireMaths::pow_2(r0 - r);
        Identities initial_forms, final_forms;

        initial_forms.set(k, kb);
        initial_forms.set(r0, r0b);

        final_forms.set(k, kf);
        final_forms.set(r0, r0f);

        return TwoAtomPerturbation( atom0, atom1,
                                    base, initial_forms, final_forms,
                                    PropertyMap("parameters",bond_property) );
    }
    else
        return PerturbationPtr();
}

void ProtoMS::processConnectLine(const QStringList &words, const Molecule &mol,
                                 int type,
                                 ConnectivityEditor &connectivity,
                                 ProtoMSWorkspace &workspace) const
{
    if (type == SOLVENT)
        return;

    if (words.count() < 11)
        throw SireError::io_error( QObject::tr(
            "Cannot understand the ProtoMS bond line\n%1")
                .arg(words.join(" ")), CODELOC );

    AtomIdx atom0, atom1;

    if (type == SOLUTE)
    {
        atom0 = getSoluteAtom(mol, words[2], words[4], workspace);
        atom1 = getSoluteAtom(mol, words[7], words[9], workspace);
    }
    else if (type == PROTEIN)
    {
        atom0 = getProteinAtom(mol, words[2], words[4], workspace);
        atom1 = getProteinAtom(mol, words[7], words[9], workspace);
    }

    connectivity.connect(atom0, atom1);
}

/** This processes the output line that contains the angle parameters */
void ProtoMS::processAngleLine(const QStringList &words, const Molecule &mol,
                               int type, ThreeAtomFunctions &anglefuncs,
                               ProtoMSWorkspace &workspace) const
{
    if (type == SOLVENT)
        return;

    Symbol theta = InternalPotential::symbols().angle().theta();

    Expression anglefunc = words[17].toDouble()
                * SireMaths::pow_2( theta - (words[19].toDouble()*degrees).to(radians) );

    AtomIdx atom0, atom1, atom2;

    if (type == SOLUTE)
    {
        atom0 = getSoluteAtom(mol, words[2], words[4], workspace);
        atom1 = getSoluteAtom(mol, words[7], words[9], workspace);
        atom2 = getSoluteAtom(mol, words[12], words[14], workspace);
    }
    else if (type == PROTEIN)
    {
        atom0 = getProteinAtom(mol, words[2], words[4], workspace);
        atom1 = getProteinAtom(mol, words[7], words[9], workspace);
        atom2 = getProteinAtom(mol, words[12], words[14],workspace);
    }

    anglefuncs.set( atom0, atom1, atom2, anglefunc );
}

/** This processes the output line that contains the angle perturbation parameters */
PerturbationPtr ProtoMS::processAnglePertLine(const QStringList &words,
                                              const Molecule &mol, int type,
                                              const PropertyName &angle_property,
                                              ProtoMSWorkspace &workspace) const
{
    if (type == SOLVENT)
        return PerturbationPtr();

    double kb = words[17].toDouble();
    double t0b = (words[19].toDouble() * degrees).to(radians);
    double kf = words[21].toDouble();
    double t0f = (words[23].toDouble() * degrees).to(radians);

    AtomIdx atom0, atom1, atom2;

    if (type == SOLUTE)
    {
        atom0 = getSoluteAtom(mol, words[2], words[4], workspace);
        atom1 = getSoluteAtom(mol, words[7], words[9], workspace);
        atom2 = getSoluteAtom(mol, words[12], words[14], workspace);
    }
    else if (type == PROTEIN)
    {
        atom0 = getProteinAtom(mol, words[2], words[4], workspace);
        atom1 = getProteinAtom(mol, words[7], words[9], workspace);
        atom2 = getProteinAtom(mol, words[12], words[14],workspace);
    }

    if (kb != kf or t0b != t0f)
    {
        Symbol theta = InternalPotential::symbols().angle().theta();
        Symbol k("k");
        Symbol theta0("theta0");

        Expression base = k * SireMaths::pow_2(theta - theta0);

        Identities initial_forms, final_forms;

        initial_forms.set(k, kb);
        initial_forms.set(theta0, t0b);
        final_forms.set(k, kf);
        final_forms.set(theta0, t0f);

        return ThreeAtomPerturbation( atom0, atom1, atom2,
                                      base, initial_forms, final_forms,
                                      PropertyMap("parameters", angle_property) );
    }
    else
        return PerturbationPtr();
}

/** This processes the output line that contains the UB parameters */
void ProtoMS::processUBLine(const QStringList &words, const Molecule &mol,
                            int type, TwoAtomFunctions &ubfuncs,
                            ProtoMSWorkspace &workspace) const
{
    if (type == SOLVENT)
        return;

    Symbol r = InternalPotential::symbols().ureyBradley().r();

    Expression ubfunc = words[12].toDouble()
                                * SireMaths::pow_2( r - words[14].toDouble() );

    AtomIdx atom0, atom1;

    if (type == SOLUTE)
    {
        atom0 = getSoluteAtom(mol, words[2], words[4], workspace);
        atom1 = getSoluteAtom(mol, words[7], words[9], workspace);
    }
    else if (type == PROTEIN)
    {
        atom0 = getProteinAtom(mol, words[2], words[4], workspace);
        atom1 = getProteinAtom(mol, words[7], words[9], workspace);
    }

    ubfuncs.set( atom0, atom1, ubfunc );
}

/** This processes the output line that contains the UB perturbation parameters */
PerturbationPtr ProtoMS::processUBPertLine(const QStringList &words, const Molecule &mol,
                                           int type, const PropertyName &ub_property,
                                           ProtoMSWorkspace &workspace) const
{
    if (type == SOLVENT)
        return PerturbationPtr();

    double kb = words[12].toDouble();
    double r0b = words[14].toDouble();

    double kf = words[16].toDouble();
    double r0f = words[18].toDouble();

    AtomIdx atom0, atom1;

    if (type == SOLUTE)
    {
        atom0 = getSoluteAtom(mol, words[2], words[4], workspace);
        atom1 = getSoluteAtom(mol, words[7], words[9], workspace);
    }
    else if (type == PROTEIN)
    {
        atom0 = getProteinAtom(mol, words[2], words[4], workspace);
        atom1 = getProteinAtom(mol, words[7], words[9], workspace);
    }

    if (kb != kf or r0b != r0f)
    {
        Symbol r = InternalPotential::symbols().ureyBradley().r();
        Symbol k("k");
        Symbol r0("r0");

        Expression base = k * SireMaths::pow_2(r - r0);

        Identities initial_forms, final_forms;

        initial_forms.set(k, kb);
        initial_forms.set(r0, r0b);
        final_forms.set(k, kf);
        final_forms.set(r0, r0b);

        return TwoAtomPerturbation(atom0, atom1,
                                   base, initial_forms, final_forms,
                                   PropertyMap("parameters", ub_property));
    }
    else
        return PerturbationPtr();
}

/** This processes the lines that contain the dihedral parameters */
QString ProtoMS::processDihedralLine(QTextStream &ts, const QStringList &words,
                                     const Molecule &mol, int type,
                                     FourAtomFunctions &dihedralfuncs,
                                     const PropertyName &dihedral_property,
                                     QList<PerturbationPtr> &perturbations,
                                     ProtoMSWorkspace &workspace) const
{
    //read in the parameter - each cosine term is on a separate line
    QString line = ts.readLine();

    Expression dihedralfunc, dihedralfunc_f, dihedralfunc_b;

    Symbol phi = InternalPotential::symbols().dihedral().phi();

    enum TYPE { REFERENCE, BACKWARDS, FORWARDS };

    bool have_perturbation = false;

    while (not line.isNull())
    {

        QStringList words = line.split(" ", Qt::SkipEmptyParts);

        if (words[0] != "PARAMS")
        {
            line = ts.readLine();
            continue;
        }

        TYPE type;

        if (words[1] == "DihedralParameter")
        {
            type = REFERENCE;
        }
        else if (words[1] == "DihedralParameterPertF")
        {
            type = FORWARDS;
            have_perturbation = true;
        }
        else if (words[1] == "DihedralParameterPertB")
        {
            type = BACKWARDS;
            have_perturbation = true;
        }
        else
            break;

        //there are four parameters, k0, k1, k2, k3
        // The cosine function is;
        //  k0 { 1 + k1 [ cos(k2*phi + k3) ] }
        double k0 = words[3].toDouble();
        double k1 = words[5].toDouble();
        double k2 = (words[7].toDouble() * radians).to(radians);
        double k3 = words[9].toDouble();

        Expression term = k0 * ( 1 + k1*( Cos(k2*phi + k3) ) );

        switch (type)
        {
            case REFERENCE:
                dihedralfunc += term;
                break;
            case FORWARDS:
                dihedralfunc_f += term;
                break;
            case BACKWARDS:
                dihedralfunc_b += term;
                break;
        }

        line = ts.readLine();
    }

    AtomIdx atom0, atom1, atom2, atom3;

    if (type == SOLVENT)
        return line;

    else if (type == SOLUTE)
    {
        atom0 = getSoluteAtom(mol, words[2], words[4], workspace);
        atom1 = getSoluteAtom(mol, words[7], words[9], workspace);
        atom2 = getSoluteAtom(mol, words[12], words[14], workspace);
        atom3 = getSoluteAtom(mol, words[17], words[19], workspace);
    }
    else if (type == PROTEIN)
    {
        atom0 = getProteinAtom(mol, words[2], words[4], workspace);
        atom1 = getProteinAtom(mol, words[7], words[9], workspace);
        atom2 = getProteinAtom(mol, words[12], words[14], workspace);
        atom3 = getProteinAtom(mol, words[17], words[19], workspace);
    }

    dihedralfuncs.set( atom0, atom1, atom2, atom3, dihedralfunc );

    if (have_perturbation and dihedralfunc_f != dihedralfunc_b)
    {
        perturbations.append( FourAtomPerturbation(atom0, atom1, atom2, atom3,
                                   dihedralfunc_b, dihedralfunc_f,
                                   PropertyMap("parameters",dihedral_property)) );
    }

    return line;
}

/** Process the line that contains information about the
    intramolecular non-bonded pairs */
void ProtoMS::processNBLine(const QStringList &words, const Molecule &mol, int type,
                            CLJNBPairs &cljpairs,
                            ProtoMSWorkspace &workspace) const
{
    AtomIdx atom0, atom1;

    if (type == SOLVENT)
        return;

    else if (type == SOLUTE)
    {
        atom0 = getSoluteAtom(mol, words[2], words[4], workspace);
        atom1 = getSoluteAtom(mol, words[7], words[9], workspace);
    }
    else if (type == PROTEIN)
    {
        atom0 = getProteinAtom(mol, words[2], words[4], workspace);
        atom1 = getProteinAtom(mol, words[7], words[9], workspace);
    }

    double cscl = words[12].toDouble();
    double ljscl = words[14].toDouble();

    cljpairs.set( mol.data().info().cgAtomIdx(atom0),
                  mol.data().info().cgAtomIdx(atom1),
                  CLJScaleFactor(cscl, ljscl) );
}

/** Internal function used to run ProtoMS to get it to
    parameterise a molecule */
Molecule ProtoMS::runProtoMS(const Molecule &molecule, int type,
                             const PropertyMap &map) const
{
    //get the names of the properties that we need
    QString charge_property = map[ parameters().charge() ].source();
    QString lj_property = map[ parameters().lj() ].source();

    QString initial_charge_property = map[ parameters().initialCharge() ].source();
    QString final_charge_property = map[ parameters().finalCharge() ].source();
    QString initial_lj_property = map[ parameters().initialLJ() ].source();
    QString final_lj_property = map[ parameters().finalLJ() ].source();

    QString bond_property = map[ parameters().bond() ].source();
    QString angle_property = map[ parameters().angle() ].source();
    QString dihedral_property = map[ parameters().dihedral() ].source();
    QString ub_property = map[ parameters().ureyBradley() ].source();
    QString nb_property = map[ parameters().nonBonded() ].source();

    QString zmatrix_property = map[ parameters().zmatrix() ].source();

    QString coords_property = map[ parameters().coordinates() ].source();
    QString connectivity_property = map[ parameters().connectivity() ].source();

    QString perts_property = map[ parameters().perturbations() ].source();

    //create a temporary directory in which to run ProtoMS
    QString tmppath = QDir::temp().absolutePath();

    if (tmppath.isEmpty())
        tmppath = QDir::temp().absolutePath();

    TempDir tmpdir(tmppath);

    //write the ProtoMS command file
    QString cmdfile = this->writeCommandFile(tmpdir, molecule, type);

    //write the file processed by the shell used to run the job
    QString shellfile = this->writeShellFile(tmpdir, cmdfile);

    Process p = Process::run("sh", shellfile);

    p.wait();

    if (p.wasKilled())
        throw SireError::process_error( QObject::tr(
                "The ProtoMS job was killed!"), CODELOC );

    if (p.isError())
    {
        QByteArray shellcontents = ::readAll(shellfile);
        QByteArray cmdcontents = ::readAll(QString("%1/protoms_input")
                                                .arg(tmpdir.path()));
        QByteArray outputcontents = ::readAll(QString("%1/protoms_output")
                                                    .arg(tmpdir.path()));

        throw SireError::process_error( QObject::tr(
            "There was an error running the ProtoMS - no output was created.\n"
            "The shell script used to run the job was;\n"
            "*****************************************\n"
            "%1\n"
            "*****************************************\n"
            "The ProtoMS input used to run the job was;\n"
            "*****************************************\n"
            "%2\n"
            "*****************************************\n"
            "The ProtoMS output was;\n"
            "*****************************************\n"
            "%3\n"
            "*****************************************\n"
            )
                .arg( QLatin1String(shellcontents),
                      QLatin1String(cmdcontents),
                      QLatin1String(outputcontents) ), CODELOC );
    }

    QFile f( QString("%1/protoms_output").arg(tmpdir.path()) );

    if ( not (f.exists() and f.open(QIODevice::ReadOnly)) )
    {
        QByteArray shellcontents = ::readAll(shellfile);
        QByteArray cmdcontents = ::readAll(QString("%1/protoms_input").arg(tmpdir.path()));

        throw SireError::process_error( QObject::tr(
            "There was an error running the ProtoMS - no output was created.\n"
            "The shell script used to run the job was;\n"
            "*****************************************\n"
            "%1\n"
            "*****************************************\n"
            "The ProtoMS input used to run the job was;\n"
            "*****************************************\n"
            "%2\n"
            "*****************************************\n")
                .arg( QLatin1String(shellcontents),
                      QLatin1String(cmdcontents) ), CODELOC );
    }

    //now read the output to parameterise the molecule
    QTextStream ts(&f);

    QString line = ts.readLine();

    MolEditor editmol = molecule.edit();

    //if we are parameterising a solute, then we need to add dummy atoms
    if (type == SOLUTE)
    {
        MolStructureEditor edit_structure_mol( editmol );

        ResStructureEditor dummies = edit_structure_mol.add( ResName("DUM") );
        dummies.renumber( ResNum(0) );

        AtomStructureEditor dm1 = dummies.add( AtomName("DM1") );
        AtomStructureEditor dm2 = dummies.add( AtomName("DM2") );
        AtomStructureEditor dm3 = dummies.add( AtomName("DM3") );

        CGStructureEditor dummy_cg = edit_structure_mol.add( CGName("DUM") );

        dm1.reparent( CGName("DUM") );
        dm2.reparent( CGName("DUM") );
        dm3.reparent( CGName("DUM") );

        //get the center of the molecule
        Vector center = molecule.evaluate().center();

        dm1.setProperty(coords_property, center);
        dm2.setProperty(coords_property, center + Vector(1,0,0));
        dm3.setProperty(coords_property, center + Vector(0,1,0));

        dm1.setProperty(charge_property, 0.0*mod_electron);
        dm2.setProperty(charge_property, 0.0*mod_electron);
        dm3.setProperty(charge_property, 0.0*mod_electron);

        dm1.setProperty(lj_property, LJParameter::dummy());
        dm2.setProperty(lj_property, LJParameter::dummy());
        dm3.setProperty(lj_property, LJParameter::dummy());

        editmol = edit_structure_mol;
    }

    ZMatrix zmatrix( editmol );

    TwoAtomFunctions bondfuncs(editmol);
    ThreeAtomFunctions anglefuncs(editmol);
    FourAtomFunctions dihedralfuncs(editmol);
    TwoAtomFunctions ubfuncs(editmol);

    ConnectivityEditor connectivity = Connectivity(editmol.data()).edit();

    CLJNBPairs nbpairs;

    if (type == SOLUTE)
    {
        //by default, say that all atom pairs are bonded
        nbpairs = CLJNBPairs( editmol.data().info(), CLJScaleFactor(0,0) );
    }
    else if (type == PROTEIN)
    {
        //by default, say that all atom pairs are non-bonded...
        nbpairs = CLJNBPairs( editmol.data().info(), CLJScaleFactor(1,1) );

        //...except for intra-residue pairs
        int nres = editmol.nResidues();

        for (ResIdx i(0); i<nres; ++i)
        {
            Residue residue = editmol.residue(i);

            Selector<Atom> atoms = residue.atoms();

            int nats = atoms.count();

            for (int j=0; j<nats; ++j)
            {
                Atom atom0 = atoms(j);

                nbpairs.set( atom0.cgAtomIdx(), atom0.cgAtomIdx(),
                             CLJScaleFactor(0,0) );

                for (int k=j+1; k<nats; ++k)
                {
                    Atom atom1 = atoms(k);

                    nbpairs.set( atom0.cgAtomIdx(), atom1.cgAtomIdx(),
                                 CLJScaleFactor(0,0) );
                    nbpairs.set( atom1.cgAtomIdx(), atom1.cgAtomIdx(),
                                 CLJScaleFactor(0,0) );
                }
            }
        }
    }

    ProtoMSWorkspace workspace;

    QStringList fatal_errors;

    QList<GeomPertPtr> geom_perturbations;
    QList<PerturbationPtr> perturbations;

    QList<QStringList> atom_pert_lines;

    PropertyMap pert_map;
    pert_map.set("coordinates", coords_property);
    pert_map.set("connectivity", connectivity_property);

    AtomSelection anchors(editmol);
    anchors.selectOnly( AtomIdx(0) );
    pert_map.set("anchors", anchors);

    while (not line.isNull())
    {
        if (line.startsWith("PARAMS "))
        {
            QStringList words = line.split(" ", Qt::SkipEmptyParts);

            if (words[1] == "ZMATRIX")
                this->processZMatrixLine(words, editmol, type, zmatrix, workspace);

            else if (words[1] == "ZMATRIXPERT")
                this->processZMatrixPertLine(words, editmol, type, geom_perturbations,
                                             zmatrix, pert_map, workspace);

            else if (words[1] == "Atom")
                this->processAtomLine(words, editmol, type,
                                      charge_property, lj_property, workspace);

            else if (words[1] == "AtomPert")
                atom_pert_lines.append(words);

            else if (words[1] == "Bond")
                this->processBondLine(words, editmol, type,
                                      bondfuncs, workspace);

            else if (words[1] == "BondPert")
            {
                PerturbationPtr pert = this->processBondPertLine(words, editmol, type,
                                                                 bond_property,
                                                                 workspace);

                if (pert.constData() != 0)
                    perturbations.append(pert);
            }
            else if (words[1] == "Connect")
                this->processConnectLine(words, editmol, type,
                                         connectivity, workspace);

            else if (words[1] == "BondDelta")
                this->processBondDeltaLine(words, editmol, type, zmatrix, workspace);

            else if (words[1] == "Angle")
                this->processAngleLine(words, editmol, type, anglefuncs, workspace);

            else if (words[1] == "AnglePert")
            {
                PerturbationPtr pert = this->processAnglePertLine(words, editmol, type,
                                                                  angle_property,
                                                                  workspace);

                if (pert.constData() != 0)
                    perturbations.append(pert);
            }
            else if (words[1] == "AngleDelta")
                this->processAngleDeltaLine(words, editmol, type, zmatrix, workspace);

            else if (words[1] == "UreyBradley")
                this->processUBLine(words, editmol, type, ubfuncs, workspace);

            else if (words[1] == "UreyBradleyPert")
            {
                PerturbationPtr pert = this->processUBPertLine(words, editmol, type,
                                                               ub_property, workspace);

                if (pert.constData() != 0)
                    perturbations.append(pert);
            }
            else if (words[1] == "Dihedral")
            {
                line = this->processDihedralLine(ts, words, editmol,
                                                 type, dihedralfuncs,
                                                 dihedral_property,
                                                 perturbations,
                                                 workspace);
                continue;
            }

            else if (words[1] == "DihedralDelta")
                this->processDihedralDeltaLine(words, editmol, type,
                                               zmatrix, workspace);

            else if (words[1] == "NB")
                this->processNBLine(words, editmol, type, nbpairs, workspace);
        }
        else if (line.startsWith("FATAL"))
        {
            fatal_errors.append(line);
        }

        line = ts.readLine();
    }

    if (not fatal_errors.isEmpty())
    {
        //something went wrong in ProtoMS
        throw SireError::process_error( QObject::tr(
            "Something went wrong in ProtoMS when parameterising the molecule. "
            "Here is the error output from ProtoMS.\n%1")
                .arg(fatal_errors.join("\n")), CODELOC );
    }

    if ( not (editmol.hasProperty(charge_property) and
              editmol.hasProperty(lj_property)) )
    {
        QStringList errors;
        errors.append( QObject::tr("The molecule is missing either the %1 (charge) "
                 "or the %2 (LJ) properties. The available properties are %3.\n"
                 "Here's the output from ProtoMS. Can you see what's gone wrong?\n")
                    .arg(charge_property, lj_property)
                    .arg(Sire::toString(editmol.propertyKeys())) );

        QFileInfo finfo( QString("%1/protoms_output").arg(tmpdir.path()) );

        QFile f2( finfo.absoluteFilePath() );

        if (not f2.open(QIODevice::ReadOnly))
        {

            errors.append( QString("Cannot open the file %1.")
                                .arg(finfo.absoluteFilePath()) );

            if (finfo.exists())
                errors.append( QString("This is weird, as the file exists "
                                       "(size == %1 bytes)").arg(finfo.size()) );
            else
                errors.append( QString(
                                    "This is because the output file doesn't exist.") );
        }
        else
        {
            QTextStream ts2(&f2);

            bool read_line = false;

            while (not ts2.atEnd())
            {
                QString line2 = ts2.readLine();
                errors.append(line2);
                read_line = true;
            }

            if (not read_line)
                errors.append(
                   QString("The output file %1 appears to be empty "
                           "(size == %d bytes)")
                                  .arg(finfo.absoluteFilePath())
                                  .arg(finfo.size()) );
        }

        throw SireError::process_error( errors.join("\n"), CODELOC );
    }

    if (not atom_pert_lines.isEmpty())
    {
        //these are the charge and LJ perturbations - we need to
        //copy the current charges and then apply the differences
        editmol.setProperty( initial_charge_property, editmol.property(charge_property) );
        editmol.setProperty( final_charge_property, editmol.property(charge_property) );

        editmol.setProperty( initial_lj_property, editmol.property(lj_property) );
        editmol.setProperty( final_lj_property, editmol.property(lj_property) );

        foreach (QStringList words, atom_pert_lines)
        {
            this->processAtomPertLine(words, editmol, type,
                                      initial_charge_property,
                                      final_charge_property,
                                      initial_lj_property,
                                      final_lj_property,
                                      workspace);
        }

        PropertyMap charge_map;
        charge_map.set("initial_charge", initial_charge_property);
        charge_map.set("charge", charge_property);
        charge_map.set("final_charge", final_charge_property);

        perturbations.append( ChargePerturbation(charge_map) );

        PropertyMap lj_map;
        lj_map.set("initial_LJ", initial_lj_property);
        lj_map.set("LJ", lj_property);
        lj_map.set("final_LJ", final_lj_property);

        perturbations.append( LJPerturbation(lj_map) );
    }

    if (type == SOLUTE or type == PROTEIN)
    {
        editmol.setProperty( zmatrix_property, zmatrix );
        editmol.setProperty( connectivity_property, connectivity.commit() );
        editmol.setProperty( bond_property, bondfuncs );
        editmol.setProperty( angle_property, anglefuncs );
        editmol.setProperty( dihedral_property, dihedralfuncs );
        editmol.setProperty( ub_property, ubfuncs );
        editmol.setProperty( nb_property, nbpairs );

        if (not geom_perturbations.isEmpty())
            perturbations.append( GeometryPerturbations(geom_perturbations) );

        if (not perturbations.isEmpty())
            editmol.setProperty( perts_property, Perturbations(perturbations) );
    }

    return editmol.commit();
}

/** Parameterise the molecule 'molecule' as a 'type' type of
    molecule (PROTEIN, SOLUTE or SOLVENT) */
Molecule ProtoMS::parameterise(const Molecule &molecule, int type,
                               const PropertyMap &map)
{
    if (type != PROTEIN and type != SOLUTE and type != SOLVENT)
        throw SireError::invalid_arg( QObject::tr(
            "Unrecognised ProtoMS molecule type (%1). Only "
            "ProtoMS::PROTEIN, ProtoMS::SOLUTE and ProtoMS::SOLVENT "
            "are supported.").arg(type), CODELOC );

    return this->runProtoMS(molecule, type, map);
}

/** Parameterise the molecules 'molecules' as 'type' type of
    molecules (PROTEIN, SOLUTE or SOLVENT) */
Molecules ProtoMS::parameterise(const Molecules &molecules, int type,
                                const PropertyMap &map)
{
    Molecules new_molecules = molecules;

    for (Molecules::const_iterator it = molecules.constBegin();
         it != molecules.constEnd();
         ++it)
    {
        new_molecules.update( this->parameterise(it->molecule(), type, map) );
    }

    return new_molecules;
}

const char* ProtoMS::typeName()
{
    return QMetaType::typeName( qMetaTypeId<ProtoMS>() );
}
