/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2017  Christopher Woods
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

#include "SireIO/charmmpsf.h"

#include "SireSystem/system.h"

#include "SireBase/parallel.h"
#include "SireBase/stringproperty.h"

#include "SireCAS/trigfuncs.h"
#include "SireCAS/sum.h"

#include "SireError/errors.h"
#include "SireIO/errors.h"

#include "SireMM/internalff.h"
#include "SireMM/twoatomfunctions.h"
#include "SireMM/threeatomfunctions.h"
#include "SireMM/fouratomfunctions.h"

#include "SireMol/atomcharges.h"
#include "SireMol/atomelements.h"
#include "SireMol/atommasses.h"
#include "SireMol/connectivity.h"
#include "SireMol/molecule.h"
#include "SireMol/moleditor.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireUnits/units.h"

#include "SireVol/periodicbox.h"

#include <QDateTime>
#include <QFileInfo>
#include <QtMath>

using namespace SireBase;
using namespace SireCAS;
using namespace SireIO;
using namespace SireMM;
using namespace SireMol;
using namespace SireStream;
using namespace SireSystem;
using namespace SireVol;

const RegisterParser<CharmmPSF> register_psf;
static const RegisterMetaType<CharmmPSF> r_psf;
static const RegisterMetaType<CharmmParam> r_charmmparam(NO_ROOT);
static const RegisterMetaType<PSFAtom> r_psfatom(NO_ROOT);

QDataStream SIREIO_EXPORT &operator<<(QDataStream &ds, const PSFAtom &psfatom)
{
    writeHeader(ds, r_psfatom, 1);

    SharedDataStream sds(ds);

    sds << psfatom.index << psfatom.mol_idx << psfatom.number << psfatom.segment
        << psfatom.res_num << psfatom.res_name << psfatom.name << psfatom.type
        << psfatom.charge << psfatom.mass;

    return ds;
}

QDataStream SIREIO_EXPORT &operator>>(QDataStream &ds, PSFAtom &psfatom)
{
    VersionID v = readHeader(ds, r_psfatom);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> psfatom.index >> psfatom.mol_idx >> psfatom.number >> psfatom.segment
            >> psfatom.res_num >> psfatom.res_name >> psfatom.name >> psfatom.type
            >> psfatom.charge >> psfatom.mass;
    }
    else
        throw version_error(v, "1", r_psfatom, CODELOC);

    return ds;
}

QDataStream SIREIO_EXPORT &operator<<(QDataStream &ds, const CharmmParam &charmmparam)
{
    writeHeader(ds, r_charmmparam, 1);

    SharedDataStream sds(ds);

    sds << charmmparam.atoms << charmmparam.params << charmmparam.type;

    return ds;
}

QDataStream SIREIO_EXPORT &operator>>(QDataStream &ds, CharmmParam &charmmparam)
{
    VersionID v = readHeader(ds, r_charmmparam);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> charmmparam.atoms >> charmmparam.params >> charmmparam.type;
    }
    else
        throw version_error(v, "1", r_charmmparam, CODELOC);

    return ds;
}

QDataStream SIREIO_EXPORT &operator<<(QDataStream &ds, const CharmmPSF &psf)
{
    writeHeader(ds, r_psf, 1);

    ds << psf.atoms << psf.bonds << psf.mol_bonds << psf.angles << psf.mol_angles
       << psf.dihedrals << psf.mol_dihedrals << psf.impropers << psf.mol_impropers
       << psf.cross_terms << psf.mol_cross_terms << psf.num_to_idx << psf.molecules
       << psf.charmm_params << psf.box << psf.has_box << static_cast<const MoleculeParser&>(psf);

    return ds;
}

QDataStream SIREIO_EXPORT &operator>>(QDataStream &ds, CharmmPSF &psf)
{
    VersionID v = readHeader(ds, r_psf);

    if (v == 1)
    {
        ds >> psf.atoms >> psf.bonds >> psf.mol_bonds >> psf.angles >> psf.mol_angles
           >> psf.dihedrals >> psf.mol_dihedrals >> psf.impropers >> psf.mol_impropers
           >> psf.cross_terms >> psf.mol_cross_terms >> psf.num_to_idx >> psf.molecules
           >> psf.charmm_params >> psf.box >> psf.has_box >> static_cast<MoleculeParser&>(psf);
    }
    else
        throw version_error(v, "1", r_psf, CODELOC);

    return ds;
}

/** Default constructor. */
PSFAtom::PSFAtom() :
    index(0),
    mol_idx(0),
    number(0),
    res_num(0),
    type("X"),
    charge(0),
    mass(0)
{
}

/** Constructor. */
PSFAtom::PSFAtom(const QString &line, int index, QStringList &errors) :
    index(index),
    mol_idx(0),
    number(0),
    res_num(0),
    type("X"),
    charge(0),
    mass(0)
{
    // Tokenize the line, splitting using a single whitespace character.
    QStringList data = line.simplified().split(QRegExp("\\s"));

    // There must be at least 8 data records.
    if (data.count() < 8)
    {
        errors.append(QObject::tr("This doesn't look like a PSF atom "
            "record! There should be at least 8 records, found %1: %2")
            .arg(data.count()).arg(line));

        return;
    }

    // Extract the atom number.
    bool ok;
    number = data[0].toInt(&ok);

    if (not ok)
    {
        errors.append(QObject::tr("Cannot extract the ID number for the atom "
            "from part (%1) of line '%2'").arg(data[0]).arg(line));

        return;
    }

    // Extract the segment name.
    segment = data[1];

    // Extract the residue number.
    res_num = data[2].toInt(&ok);

    if (not ok)
    {
        errors.append(QObject::tr("Cannot extract the residue ID number for the atom "
            "from part (%1) of line '%2'").arg(data[2]).arg(line));

        return;
    }

    // Extract the residue name.
    res_name = data[3];

    // Extract the atom name.
    name = data[4];

    // Extract atom type.
    type = data[5];

    // Extract atom charge.
    charge = data[6].toDouble(&ok);

    if (not ok)
    {
        errors.append(QObject::tr("Cannot extract the atom charge "
            "from part (%1) of line '%2'").arg(data[6]).arg(line));

        return;
    }

    // Extract atom mass.
    mass = data[7].toDouble(&ok);

    if (not ok)
    {
        errors.append(QObject::tr("Cannot extract the atom mass "
            "from part (%1) of line '%2'").arg(data[7]).arg(line));

        return;
    }
}

/** Constructor. */
PSFAtom::PSFAtom(const SireMol::Atom &atom, const QString &segment,
    QStringList &errors, const PropertyMap &map) :
    index(0),
    mol_idx(0),
    number(atom.number().value()),
    segment(segment),
    res_num(0),
    name(atom.name().value().toUpper()),
    type("X"),
    charge(0),
    mass(0)
{
    // The atom is within a residue.
    if (atom.isWithinResidue())
    {
        res_name = atom.residue().name().value().left(4);
        res_num  = atom.residue().number().value();
    }

    // Extract the atom type.
    if (atom.hasProperty(map["atomtype"]))
    {
        type = atom.property<QString>(map["atomtype"]);
    }

    // Extract the atomic charge.
    if (atom.hasProperty("charge"))
    {
        charge = atom.property<SireUnits::Dimension::Charge>(map["charge"]).value();
    }

    // Extract the mass.
    if (atom.hasProperty(map["mass"]))
    {
        mass = atom.property<SireUnits::Dimension::MolarMass>(map["mass"]).value();
    }
}

/** Generate a PSF record from the atom data. */
QString PSFAtom::toPSFRecord() const
{
    QString line;

    line.append(QString("%1").arg(number, 8));
    line.append(QString(" %1").arg(segment, -4));
    line.append(QString(" %1").arg(res_num, -4));
    line.append(QString(" %1").arg(res_name, -4));
    line.append(QString(" %1").arg(name, -4));
    line.append(QString(" %1").arg(type, -4));
    line.append(QString(" %1").arg(charge, 10, 'f', 6));
    line.append(QString(" %1").arg(mass, 13, 'f', 4));
    line.append(QString("%1").arg("0", 12));

    return line;
}

const char* PSFAtom::typeName()
{
    return QMetaType::typeName( qMetaTypeId<PSFAtom>() );
}

/** Get the atom index. */
int PSFAtom::getIndex() const
{
    return index;
}

/** Get the molecule index. */
int PSFAtom::getMolIndex() const
{
    return mol_idx;
}

/** Set the molecule index. */
void PSFAtom::setMolIndex(int index)
{
    mol_idx = index;
}

/** Get the atom number. */
int PSFAtom::getNumber() const
{
    return number;
}

/** Get the segment name. */
QString PSFAtom::getSegment() const
{
    return segment;
}

/** Get the residue number. */
qint64 PSFAtom::getResNum() const
{
    return res_num;
}

/** Get the residue name. */
QString PSFAtom::getResName() const
{
    return res_name;
}

/** Get the atom name. */
QString PSFAtom::getName() const
{
    return name;
}

/** Get the atom type. */
QString PSFAtom::getType() const
{
    return type;
}

/** Get the atom charge. */
double PSFAtom::getCharge() const
{
    return charge;
}

/** Get the atom mass. */
double PSFAtom::getMass() const
{
    return mass;
}

/** Default constructor. */
CharmmParam::CharmmParam() : type(-1)
{
}

/** Constructor.

    @param line
        A record line from a parameter file.

    @param type
        The type of parameter record: 0 = bond
                                      1 = angle
                                      2 = dihedral
                                      3 = improper
                                      4 = non-bonded

    @param errors
        A list of parse errors.

    @param is_xplor
        Whether this is an xplor format record.
 */
CharmmParam::CharmmParam(const QString& line, int type, QStringList &errors, bool is_xplor)
    : type(type)
{
    // Tokenize the line.
    // First split on the comment identifier '!', take the first record,
    // i.e. the data, then split this by a single whitespace character.
    QStringList data = line.simplified().split("!")[0].split(QRegExp("\\s"));

    // Bond parameters.
    if (type == 0)
    {
        // Check record count.
        if (data.count() < 4)
        {
            errors.append(QObject::tr("This doesn't look like a CHARMM bond "
                "parameter record! There should be at least 4 entries, found %1: %2")
                .arg(data.count()).arg(line));

            return;
        }

        // Atom data.
        atoms.append(data[0]);
        atoms.append(data[1]);

        // Parameter data.
        bool ok1, ok2;

        // Attempt to read the parameter values.
        double p1 = data[2].toDouble(&ok1);
        double p2 = data[3].toDouble(&ok2);

        if (not ok1 or not ok2)
        {
            errors.append(QObject::tr("Could not read CHARMM bond parameter record! %1")
                .arg(line));

            return;
        }

        // Append the parameters.
        params.append(p1);
        params.append(p2);
    }

    // Angle parameters.
    else if (type == 1)
    {
        // Check record count.
        if (data.count() < 5)
        {
            errors.append(QObject::tr("This doesn't look like a CHARMM angle "
                "parameter record! There should be at least 5 entries, found %1: %2")
                .arg(data.count()).arg(line));

            return;
        }

        // Atom data.
        atoms.append(data[0]);
        atoms.append(data[1]);
        atoms.append(data[2]);

        // Parameter data.
        bool ok1, ok2;

        // Attempt to read the parameter values.
        double p1 = data[3].toDouble(&ok1);
        double p2 = data[4].toDouble(&ok2);

        if (not ok1 or not ok2)
        {
            errors.append(QObject::tr("Could not read CHARMM angle parameter record! %1")
                .arg(line));

            return;
        }

        // Append the parameters.
        params.append(p1);
        params.append(p2);

        // Check whether there are Urey-Bradley parameters.
        if ((data.count() > 6) and (data[5].at(0) != '!'))
        {
            // Attempt to read the parameter values.

            if (data[5] == "UB")
            {
                p1 = data[6].toDouble(&ok1);
                p2 = data[7].toDouble(&ok2);
            }
            else
            {
                p1 = data[5].toDouble(&ok1);
                p2 = data[6].toDouble(&ok2);
            }

            if (not ok1 or not ok2)
            {
                errors.append(QObject::tr("Could not read CHARMM Urey-Bradley angle "
                    "parameters '%1' and '%2' from line '%3'")
                    .arg(data[5]).arg(data[6]).arg(line));

                return;
            }

            // Append the parameters.
            params.append(p1);
            params.append(p2);
        }
    }

    // Dihedral parameters.
    else if (type == 2)
    {
        // Check record count.
        if (data.count() < 7)
        {
            errors.append(QObject::tr("This doesn't look like a CHARMM dihedral "
                "parameter record! There should be at least 7 entries, found %1: %2")
                .arg(data.count()).arg(line));

            return;
        }

        // Atom data.
        atoms.append(data[0]);
        atoms.append(data[1]);
        atoms.append(data[2]);
        atoms.append(data[3]);

        // Parameter data.
        bool ok1, ok2, ok3;

        // Attempt to read the parameter values.
        double p1 = data[4].toDouble(&ok1);
        double p2 = data[5].toDouble(&ok2);
        double p3 = data[6].toDouble(&ok3);

        if (not ok1 or not ok2 or not ok3)
        {
            errors.append(QObject::tr("Could not read CHARMM dihedral parameter record! %1")
                .arg(line));

            return;
        }

        // Append the parameters.
        params.append(p1);
        params.append(p2);
        params.append(p3);
    }

    // Improper parameters.
    else if (type == 3)
    {
        // Check record count.
        if (data.count() < 7)
        {
            errors.append(QObject::tr("This doesn't look like a CHARMM improper "
                "parameter record! There should be at least 6 entries, found %1: %2")
                .arg(data.count()).arg(line));

            return;
        }

        // Atom data.
        atoms.append(data[0]);
        atoms.append(data[1]);
        atoms.append(data[2]);
        atoms.append(data[3]);

        // Parameter data.
        bool ok1, ok2;

        // Attempt to read the parameter values.
        double p1 = data[4].toDouble(&ok1);
        double p2 = data[6].toDouble(&ok2);

        if (not ok1 or not ok2)
        {
            errors.append(QObject::tr("Could not read CHARMM improper parameter record! %1")
                .arg(line));

            return;
        }

        // Append the parameters.
        params.append(p1);
        params.append(p2);
    }

    // Non-bonded parameters.
    else if (type == 4)
    {
        int min_records = 3;
        if (not is_xplor) min_records = 4;

        // Check record count.
        if (data.count() < 3)
        {
            errors.append(QObject::tr("This doesn't look like a CHARMM non-bonded "
                "parameter record! There should be at least %1 entries, found %2: %3")
                .arg(min_records).arg(data.count()).arg(line));

            return;
        }

        // Atom data.
        atoms.append(data[0]);

        // Parameter data.
        bool ok1, ok2;
        double p1, p2;

        // Attempt to read the parameter values.
        if (is_xplor)
        {
            p1 = data[1].toDouble(&ok1);
            p2 = data[2].toDouble(&ok2);
        }
        else
        {
            p1 = data[2].toDouble(&ok1);
            p2 = data[3].toDouble(&ok2);
        }

        if (not ok1 or not ok2)
        {
            errors.append(QObject::tr("Could not read CHARMM non-bonded parameter record! %1")
                .arg(line));

            return;
        }

        // Append the parameters.
        params.append(p1);
        params.append(p2);

        // Now check for 1-4 non-bonded parameters.
        bool is_nb14 = false;

        if (is_xplor)
        {
            if ((data.count() > 4) and (data[3].at(0) != '!'))
                is_nb14 = true;
        }
        else
        {
            if ((data.count() > 6) and (data[4].at(0) != '!'))
                is_nb14 = true;
        }

        if (is_nb14)
        {
            if (is_xplor)
            {
                p1 = data[3].toDouble(&ok1);
                p2 = data[4].toDouble(&ok2);
            }
            else
            {
                p1 = data[5].toDouble(&ok1);
                p2 = data[6].toDouble(&ok2);
            }

            if (not ok1 or not ok2)
            {
                errors.append(QObject::tr("Could not read CHARMM non-bonded 1-4 parameter record! %1")
                    .arg(line));

                return;
            }

            // Append the parameters.
            params.append(p1);
            params.append(p2);
        }
    }

    // Uknown parameter type.
    else
    {
        throw SireError::program_bug(QObject::tr("Unknown parameter type (%1). "
            "Valid types are 0, 1, 2, 3, 4").arg(type), CODELOC);
    }
}

/** Return the atoms to which the parameters apply. */
const QVector<QString>& CharmmParam::getAtoms() const
{
    return atoms;
}

/** Return the parameter terms. */
const QVector<double>& CharmmParam::getParams() const
{
    return params;
}

/** Return the parameter type. */
qint64 CharmmParam::getType() const
{
    return type;
}

const char* CharmmParam::typeName()
{
    return QMetaType::typeName( qMetaTypeId<CharmmParam>() );
}

/** Constructor */
CharmmPSF::CharmmPSF() :
    ConcreteProperty<CharmmPSF,MoleculeParser>(), has_box(false)
{}

/** Construct to read in the data from the file called 'filename'. The
passed property map can be used to pass extra parameters to control
the parsing */
CharmmPSF::CharmmPSF(const QString &filename, const PropertyMap &map) :
    ConcreteProperty<CharmmPSF,MoleculeParser>(filename, map), has_box(false)
{
    //the file has been read into memory and is available via
    //the MoleculeParser::lines() function.

    //a parameter has also been read in MoleculeParser to say whether
    //we are allowed to use multiple cores to parse the file, e.g.
    //MoleculeParser::usesParallel() will be true

    //parse the data in the parse function
    this->parseLines(map);

    //now make sure that everything is correct with this object
    this->assertSane();
}

/** Construct to read in the data from the passed text lines. The
    passed property map can be used to pass extra parameters to control
    the parsing */
CharmmPSF::CharmmPSF(const QStringList &lines, const PropertyMap &map) :
    ConcreteProperty<CharmmPSF,MoleculeParser>(lines, map), has_box(false)
{
    //the file has been read into memory and is available via
    //the MoleculeParser::lines() function.

    //a parameter has also been read in MoleculeParser to say whether
    //we are allowed to use multiple cores to parse the file, e.g.
    //MoleculeParser::usesParallel() will be true

    //parse the data in the parse function
    this->parseLines(map);

    //now make sure that everything is correct with this object
    this->assertSane();
}

/** Construct this parser by extracting all necessary information from the
    passed SireSystem::System, looking for the properties that are specified
    in the passed property map */
CharmmPSF::CharmmPSF(const SireSystem::System &system, const PropertyMap &map) :
    ConcreteProperty<CharmmPSF,MoleculeParser>(map), has_box(false)
{
    // Get the MolNums of each molecule in the System - this returns the
    // numbers in MolIdx order.
    const QVector<MolNum> molnums = system.getMoleculeNumbers().toVector();

    // Store the number of molecules.
    const int nmols = molnums.count();

    // No molecules in the system.
    if (nmols == 0)
    {
        this->operator=(CharmmPSF());
        return;
    }

    // The unique CHARMM parameter strings.
    QSet<QString> bond_params;
    QSet<QString> angle_params;
    QSet<QString> dihedral_params;
    QSet<QString> improper_params;

    // Reconstruct the record data for each molecule.
    for (int i=0; i<nmols; ++i)
    {
        // Create local data objects.
        QVector<PSFAtom> local_atoms;
        QVector<QVector<qint64> > local_bonds;
        QVector<QVector<qint64> > local_angles;
        QVector<QVector<qint64> > local_dihedrals;
        QVector<QVector<qint64> > local_impropers;

        // Parse the molecular system.
        parseMolecule(i, system[molnums[i]].molecule(), local_atoms, local_bonds,
            local_angles, local_dihedrals, local_impropers, bond_params, angle_params,
            dihedral_params, improper_params, parse_warnings, map);

        atoms     += local_atoms;
        bonds     += local_bonds;
        angles    += local_angles;
        dihedrals += local_dihedrals;
        impropers += local_impropers;
    }

    // Generate the PSF record lines.
    QVector<QString> lines = toLines();

    // Convert the vector of records to a list.
    QStringList lines_list;
	for (const auto &line : lines)
		lines_list << line;

    // Extract the existing CHARMM parameters.
    try
    {
        // Extract the existing parameters from the system.
        QString params = system.property(map["charmm-params"]).asA<StringProperty>().value();

        // Convert into a list.
        charmm_params = params.split("\n");
    }
    catch (...)
    {
        // Generate a new list of CHARMM parameter records.
        // The records strings are sorted before appending to the list.

        // Bonds.
        if (bond_params.count() > 0)
        {
            charmm_params.append("BONDS");
            charmm_params.append("!");
            charmm_params.append("!V(bond) = Kb(b - b0)**2");
            charmm_params.append("!");
            charmm_params.append("!Kb: kcal/mole/A**2");
            charmm_params.append("!");
            charmm_params.append("!atom type   Kb          b0");
            charmm_params.append("!");
            QStringList params = bond_params.toList();
            params.sort();
            charmm_params.append(params);
            charmm_params.append("");
        }

        // Angles.
        if (angle_params.count() > 0)
        {
            charmm_params.append("ANGLES");
            charmm_params.append("!");
            charmm_params.append("!V(angle) = Ktheta(Theta - Theta0)**2");
            charmm_params.append("!");
            charmm_params.append("!V(Urey-Bradley) = Kub(S - S0)**2");
            charmm_params.append("!");
            charmm_params.append("!Ktheta: kcal/mole/rad**2");
            charmm_params.append("!Theta0: degrees");
            charmm_params.append("!Kub: kcal/mole/A**2 (Urey-Bradley)");
            charmm_params.append("!S0: A");
            charmm_params.append("!");
            charmm_params.append("!atom types        Ktheta   Theta0       Kub        S0");
            charmm_params.append("!");
            QStringList params = angle_params.toList();
            params.sort();
            charmm_params.append(params);
            charmm_params.append("");
        }

        // Dihedrals.
        if (dihedral_params.count() > 0)
        {
            charmm_params.append("DIHEDRALS");
            charmm_params.append("!");
            charmm_params.append("!V(dihedral) = Kchi(1 + cos(n(chi) - delta))");
            charmm_params.append("!");
            charmm_params.append("!Kchi: kcal/mole");
            charmm_params.append("!n: multiplicity");
            charmm_params.append("!delta: degrees");
            charmm_params.append("!");
            charmm_params.append("!atom types             Kchi     n     delta");
            charmm_params.append("!");
            QStringList params = dihedral_params.toList();
            params.sort();
            charmm_params.append(params);
            charmm_params.append("");
        }

        // Impropers.
        if (improper_params.count() > 0)
        {
            charmm_params.append("IMPROPER");
            charmm_params.append("!");
            charmm_params.append("!V(improper) = Kpsi(psi - psi0)**2");
            charmm_params.append("!");
            charmm_params.append("!Kpsi: kcal/mole/rad**2");
            charmm_params.append("!psi0: degrees");
            charmm_params.append("!note that the second column of numbers (0) is ignored");
            charmm_params.append("!");
            charmm_params.append("!atom types             Kpsi                  psi0");
            charmm_params.append("!");
            QStringList params = improper_params.toList();
            params.sort();
            charmm_params.append(params);
            charmm_params.append("");
        }

        if (charmm_params.count() > 0)
            charmm_params.append("END");
    }

    // Flag whether the system has a periodic simulation box.
    has_box = false;

    try
    {
        box = system.property(map["space"]).asA<PeriodicBox>();
        has_box = true;
    }
    catch(...)
    {}

    // Reparse the lines as a self-consistency check.
    CharmmPSF parsed(lines_list, map);

    // Set the CHARMM parameters.
    parsed.charmm_params = charmm_params;

    // Set the box.
    if (has_box)
    {
        parsed.box = box;
        parsed.has_box = true;
    }

    this->operator=(parsed);
}

/** Copy constructor */
CharmmPSF::CharmmPSF(const CharmmPSF &other) :
    ConcreteProperty<CharmmPSF,MoleculeParser>(other),
    atoms(other.atoms),
    bonds(other.bonds),
    mol_bonds(other.mol_bonds),
    angles(other.angles),
    mol_angles(other.mol_angles),
    dihedrals(other.dihedrals),
    mol_dihedrals(other.mol_dihedrals),
    impropers(other.impropers),
    mol_impropers(other.mol_impropers),
    cross_terms(other.cross_terms),
    mol_cross_terms(other.mol_cross_terms),
    num_to_idx(other.num_to_idx),
    molecules(other.molecules),
    charmm_params(other.charmm_params),
    box(other.box),
    has_box(other.has_box),
    parse_warnings(other.parse_warnings)
{}

/** Destructor */
CharmmPSF::~CharmmPSF()
{}

/** Copy assignment operator */
CharmmPSF& CharmmPSF::operator=(const CharmmPSF &other)
{
    if (this != &other)
    {
        molecules = other.molecules;
        bonds = other.bonds;
        mol_bonds = other.mol_bonds;
        angles = other.angles;
        mol_angles = other.mol_angles;
        dihedrals = other.dihedrals;
        mol_dihedrals = other.mol_dihedrals;
        impropers = other.impropers;
        mol_impropers = other.mol_impropers;
        cross_terms = other.cross_terms;
        mol_cross_terms = other.mol_cross_terms;
        num_to_idx = other.num_to_idx;
        molecules = other.molecules;
        charmm_params = other.charmm_params;
        box = other.box;
        has_box = other.has_box;
        parse_warnings = other.parse_warnings;

        MoleculeParser::operator=(other);
    }

    return *this;
}

/** Comparison operator */
bool CharmmPSF::operator==(const CharmmPSF &other) const
{
    return MoleculeParser::operator==(other);
}

/** Comparison operator */
bool CharmmPSF::operator!=(const CharmmPSF &other) const
{
    return not operator==(other);
}

/** Return the C++ name for this class */
const char* CharmmPSF::typeName()
{
    return QMetaType::typeName( qMetaTypeId<CharmmPSF>() );
}

/** Return the C++ name for this class */
const char* CharmmPSF::what() const
{
    return CharmmPSF::typeName();
}

/** Return the parser that has been constructed by reading in the passed
    file using the passed properties */
MoleculeParserPtr CharmmPSF::construct(const QString &filename,
                                  const PropertyMap &map) const
{
    return CharmmPSF(filename,map);
}

/** Return the parser that has been constructed by reading in the passed
    text lines using the passed properties */
MoleculeParserPtr CharmmPSF::construct(const QStringList &lines,
                                  const PropertyMap &map) const
{
    return CharmmPSF(lines,map);
}

/** Return the parser that has been constructed by extract all necessary
    data from the passed SireSystem::System using the specified properties */
MoleculeParserPtr CharmmPSF::construct(const SireSystem::System &system,
                                  const PropertyMap &map) const
{
    return CharmmPSF(system,map);
}

/** Return a string representation of this parser */
QString CharmmPSF::toString() const
{
    if (lines().isEmpty())
        return QObject::tr("CharmmPSF::null");
    else
    {
        return QObject::tr("CharmmPSF( nMolecules() = %1, nAtoms() = %2, "
            "nBonds() = %3, nAngles() = %4, nDihedrals() = %5, "
            "nImpropers() = %6, nCrossTerms() = %7 )")
            .arg(nMolecules()).arg(nAtoms()).arg(nBonds()).arg(nAngles())
            .arg(nDihedrals()).arg(nImpropers()).arg(nCrossTerms());
    }
}

/** Convert the parsed data to a collection of PSF record lines. */
QVector<QString> CharmmPSF::toLines() const
{
    // A good description of the PSF format is given here:
    // http://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-win-html/node24.html

    // No records.
    if ((nAtoms() + nBonds() + nAngles() + nDihedrals()
        + nImpropers() + nCrossTerms()) == 0)
    {
        return QVector<QString>();
    }

    // Helper function to reconstruct record lines.
    auto generate_line = [](const QVector<QVector<qint64> > &data, int iline,
        int num_lines, int num_records, int record_width)
    {
        // The record line.
        QString line;

        // The default number of columns per line.
        int num_cols = 8;

        // Triples have three entries per line, so 9 columns.
        if (record_width == 3) num_cols = 9;

        // Work out the number of records per line.
        int records_per_line = num_cols / record_width;

        // Work out the row index in the data vector for the first record.
        int start = records_per_line * iline;

        // Adjust the records per line if we're on the last record line and
        // the total number of records is not a multiple of the records per line.
        if (iline == (num_lines - 1))
            records_per_line = num_records - records_per_line*(num_lines - 1);

        // Add each record to the line.
        for (int i=0; i<records_per_line; ++i)
        {
            for (int j=0; j<record_width; j++)
            {
                line.append(QString("%1").arg(data[start+i][j], 8));
            }
        }

        return line;
    };

    // The vector of record lines.
    QVector<QString> lines;

    // Add header information.
    // TODO: Add Sire version info.
    lines.append("PSF");
    lines.append("");
    lines.append("       1 !NTITLE");
    lines.append(QString(" REMARKS DATE:%1    created by Sire")
         .arg(QDateTime::currentDateTime().toString("dd-MMM-yy  hh:mm:ss")));
    lines.append("");

    // Add any atom records.
    if (nAtoms() > 0)
    {
        const int num_records = nAtoms();

        QVector<QString> record_lines(num_records + 2);
        record_lines[0] = QString("%1 !NATOM").arg(num_records, 8);

        if (usesParallel())
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, num_records),
                            [&](const tbb::blocked_range<int> &r)
            {
                for (int i=r.begin(); i<r.end(); ++i)
                {
                    record_lines[i+1] = atoms[i].toPSFRecord();
                }
            });
        }
        else
        {
            for (int i=0; i<num_records; ++i)
            {
                record_lines[i+1] = atoms[i].toPSFRecord();
            }
        }

        // Append the atom record lines.
        lines += record_lines;
    }

    // Add any bond records.
    if (nBonds() > 0)
    {
        const int num_records = nBonds();

        // There are 4 bond records per line.
        const int num_lines = qCeil(num_records/4.0);

        QVector<QString> record_lines(num_lines + 2);
        record_lines[0] = QString("%1 !NBOND: bonds").arg(num_records, 8);

        if (usesParallel())
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, num_lines),
                            [&](const tbb::blocked_range<int> &r)
            {
                for (int i=r.begin(); i<r.end(); ++i)
                {
                    record_lines[i+1] = generate_line(bonds, i, num_lines, num_records, 2);
                }
            });
        }
        else
        {
            for (int i=0; i<num_lines; ++i)
            {
                record_lines[i+1] = generate_line(bonds, i, num_lines, num_records, 2);
            }
        }

        // Append the bond record lines.
        lines += record_lines;
    }

    // Add any angle records.
    if (nAngles() > 0)
    {
        const int num_records = nAngles();

        // There are 3 angle records per line.
        const int num_lines = qCeil(num_records/3.0);

        QVector<QString> record_lines(num_lines + 2);
        record_lines[0] = QString("%1 !NTHETA: angles").arg(num_records, 8);

        if (usesParallel())
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, num_lines),
                            [&](const tbb::blocked_range<int> &r)
            {
                for (int i=r.begin(); i<r.end(); ++i)
                {
                    record_lines[i+1] = generate_line(angles, i, num_lines, num_records, 3);
                }
            });
        }
        else
        {
            for (int i=0; i<num_lines; ++i)
            {
                record_lines[i+1] = generate_line(angles, i, num_lines, num_records, 3);
            }
        }

        // Append the angle record lines.
        lines += record_lines;
    }

    // Add any dihedral records.
    if (nDihedrals() > 0)
    {
        const int num_records = nDihedrals();

        // There are 2 dihedral records per line.
        const int num_lines = qCeil(num_records/2.0);

        QVector<QString> record_lines(num_lines + 2);
        record_lines[0] = QString("%1 !NPHI: dihedrals").arg(num_records, 8);

        if (usesParallel())
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, num_lines),
                            [&](const tbb::blocked_range<int> &r)
            {
                for (int i=r.begin(); i<r.end(); ++i)
                {
                    record_lines[i+1] = generate_line(dihedrals, i, num_lines, num_records, 4);
                }
            });
        }
        else
        {
            for (int i=0; i<num_lines; ++i)
            {
                record_lines[i+1] = generate_line(dihedrals, i, num_lines, num_records, 4);
            }
        }

        // Append the dihedral record lines.
        lines += record_lines;
    }

    // Add any improper records.
    if (nImpropers() > 0)
    {
        const int num_records = nImpropers();

        // There are 2 improper records per line.
        const int num_lines = qCeil(num_records/2.0);

        QVector<QString> record_lines(num_lines + 2);
        record_lines[0] = QString("%1 !NIMPHI: impropers").arg(num_records, 8);

        if (usesParallel())
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, num_lines),
                            [&](const tbb::blocked_range<int> &r)
            {
                for (int i=r.begin(); i<r.end(); ++i)
                {
                    record_lines[i+1] = generate_line(impropers, i, num_lines, num_records, 4);
                }
            });
        }
        else
        {
            for (int i=0; i<num_lines; ++i)
            {
                record_lines[i+1] = generate_line(impropers, i, num_lines, num_records, 4);
            }
        }

        // Append the improper record lines.
        lines += record_lines;
    }

    // Add any cross-term records.
    if (nCrossTerms() > 0)
    {
        const int num_records = nCrossTerms();

        // There are 2 improper records per line.
        const int num_lines = qCeil(num_records/2.0);

        QVector<QString> record_lines(num_lines + 2);
        record_lines[0] = QString("%1 !NCRTERM: cross-terms").arg(num_records, 8);

        if (usesParallel())
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, num_lines),
                            [&](const tbb::blocked_range<int> &r)
            {
                for (int i=r.begin(); i<r.end(); ++i)
                {
                    record_lines[i+1] = generate_line(cross_terms, i, num_lines, num_records, 4);
                }
            });
        }
        else
        {
            for (int i=0; i<num_lines; ++i)
            {
                record_lines[i+1] = generate_line(cross_terms, i, num_lines, num_records, 4);
            }
        }

        // Append the cross-term record lines.
        lines += record_lines;
    }

    return lines;
}

/** Return the format name that is used to identify this file format within Sire */
QString CharmmPSF::formatName() const
{
    return "PSF";
}

/** Return a description of the file format */
QString CharmmPSF::formatDescription() const
{
    return QObject::tr("Charmm PSF format files.");
}

/** Return the suffixes that these files are normally associated with */
QStringList CharmmPSF::formatSuffix() const
{
    static const QStringList suffixes = { "psf" };
    return suffixes;
}

/** Return whether or not this is a lead parser. The lead parser is responsible
    for starting the process of turning the parsed file into the System. There
    must be one and one-only lead parser in a set of parsers creating a System */
bool CharmmPSF::isLead() const
{
    return true;
}

/** Return whether or not this parser can follow another lead parser, and add
    data to an existing molecular system. The CharmmPSF parser cannot follow. */
bool CharmmPSF::canFollow() const
{
    return false;
}

/** Return the number of molecules. */
int CharmmPSF::nMolecules() const
{
    return molecules.count();
}

/** Return the number of atom records. */
int CharmmPSF::nAtoms() const
{
    return atoms.count();
}

/** Return the number of atoms in molecule i. */
int CharmmPSF::nAtoms(int i) const
{
    return molecules[i].count();
}

/** Return the number of bond records. */
int CharmmPSF::nBonds() const
{
    return bonds.count();
}

/** Return the number of bonds in molecule 'i'. */
int CharmmPSF::nBonds(int i) const
{
    return mol_bonds[i].count();
}

/** Return the number of angle records. */
int CharmmPSF::nAngles() const
{
    return angles.count();
}

/** Return the number of angles in molecule 'i'. */
int CharmmPSF::nAngles(int i) const
{
    return mol_angles[i].count();
}

/** Return the number of dihedral records. */
int CharmmPSF::nDihedrals() const
{
    return dihedrals.count();
}

/** Return the number of dihedrals in molecule 'i'. */
int CharmmPSF::nDihedrals(int i) const
{
    return mol_dihedrals[i].count();
}

/** Return the number of improper records. */
int CharmmPSF::nImpropers() const
{
    return impropers.count();
}

/** Return the number of impropers in molecule 'i'. */
int CharmmPSF::nImpropers(int i) const
{
    return mol_impropers[i].count();
}

/** Return the number of cross-term records. */
int CharmmPSF::nCrossTerms() const
{
    return cross_terms.count();
}

/** Return the number of cross-terms in molecule 'i'. */
int CharmmPSF::nCrossTerms(int i) const
{
    return mol_cross_terms[i].count();
}

/** Function that is called to assert that this object is sane. This
    should raise an exception if the parser is in an invalid state */
void CharmmPSF::assertSane() const
{
    //check state, raise SireError::program_bug if we are in an invalid state
}

/** Internal function that is used to actually parse the data contained
    in the lines of the file */
void CharmmPSF::parseLines(const PropertyMap &map)
{
    int num_atoms = 0;

    // Helper function to parse record lines.
    auto parse_line = [](const QString& line, QVector<QVector<qint64> > &data, int iline,
        int num_lines, int num_records, int record_width, QStringList &errors)
    {
        // The default number of columns per line.
        int num_cols = 8;

        // Triples have three entries per line, so 9 columns.
        if (record_width == 3) num_cols = 9;

        // Work out the number of records per line.
        int records_per_line = num_cols / record_width;

        // Work out the row index in the data vector for the first record.
        int start = records_per_line * iline;

        // Adjust the records per line if we're on the last record line and
        // the total number of records is not a multiple of the records per line.
        if (iline == (num_lines - 1))
            records_per_line = num_records - records_per_line*(num_lines - 1);

        // Tokenize the line, splitting using a single whitespace character.
        QStringList records = line.simplified().split(QRegExp("\\s"));

        // Check that the line has the right number of columns.
        if (records.count() != (record_width*records_per_line))
        {
            errors.append(QObject::tr("The number of records on the line is "
                "incorrect. Expected %1, found %2: %3")
                .arg(num_cols * record_width).arg(records.count()).arg(line));

            return;
        }

        bool ok;

        // Validate and store each record.
        for (int i=0; i<records_per_line; ++i)
        {
            for (int j=0; j<record_width; j++)
            {
                data[start+i][j] = records[i*record_width + j].toInt(&ok);

                if (not ok)
                {
                    errors.append(QObject::tr("Could not extract data record %1 "
                        "from line '%2'").arg(records[i]).arg(line));

                    return;
                }
            }
        }
    };

    // Loop through all lines in the file.
    for (int iline=0; iline<lines().count(); ++iline)
    {
        // Store the line.
        const QString line = lines()[iline];

        // Tokenize the line, splitting using a single whitespace character.
        QStringList data = lines()[iline].simplified().split(QRegExp("\\s"));

        // Atom records.
        if (line.contains("!NATOM"))
        {
            // Extract the number of atoms (should be the first record).
            bool ok;
            num_atoms = data.first().toInt(&ok);

            if (not ok)
            {
                parse_warnings.append(QObject::tr("Cannot extract number of atoms "
                    "from part (%1) from line '%2'").arg(data.first()).arg(lines()[iline]));

                return;
            }

            if (num_atoms == 0)
            {
                parse_warnings.append(QObject::tr("The molecule contains no atoms!: %1")
                    .arg(lines()[iline]));

                return;
            }

            // Resize the atoms vector.
            atoms.resize(num_atoms);

            if (usesParallel())
            {
                ++iline;

                QMutex mutex;

                tbb::parallel_for(tbb::blocked_range<int>(0, num_atoms),
                                [&](const tbb::blocked_range<int> &r)
                {
                    // Create local data objects.
                    QStringList local_errors;

                    for (int i=r.begin(); i<r.end(); ++i)
                    {
                        // Parse the data from the atom record.
                        atoms[i] = PSFAtom(lines()[iline+i], i, local_errors);
                    }

                    if (not local_errors.isEmpty())
                    {
                        // Acquire a lock.
                        QMutexLocker lkr(&mutex);

                        // Update the warning messages.
                        parse_warnings += local_errors;
                    }
                });

                // Fast-forward the line index.
                iline += (num_atoms - 1);
            }
            else
            {
                for (int i=0; i<num_atoms; ++i)
                {
                    // Parse the data from the atom record.
                    atoms[i] = PSFAtom(lines()[++iline], i, parse_warnings);
                }
            }
        }

        // Bond records.
        else if (line.contains("!NBOND"))
        {
            // Extract the number of bonds;
            bool ok;
            int num_bonds = data.first().toInt(&ok);

            if (not ok)
            {
                parse_warnings.append(QObject::tr("Cannot extract number of bonds "
                    "from part (%1) from line '%2'").arg(data.first()).arg(lines()[iline]));

                return;
            }

            // Resize data structures.
            bonds.resize(num_bonds);
            for (int i=0; i<num_bonds; ++i)
                bonds[i].resize(2);

            // Work out the number of record lines.
            // There are 4 bond records per line.
            const int num_lines = qCeil(num_bonds/4.0);

            ++iline;

            if (usesParallel())
            {
                QMutex mutex;

                tbb::parallel_for(tbb::blocked_range<int>(0, num_lines),
                                [&](const tbb::blocked_range<int> &r)
                {
                    // Create local data objects.
                    QStringList local_errors;

                    for (int i=r.begin(); i<r.end(); ++i)
                    {
                        // Parse a bond record line.
                        parse_line(lines()[iline+i], bonds, i,
                            num_lines, num_bonds, 2, local_errors);
                    }

                    if (not local_errors.isEmpty())
                    {
                        // Acquire a lock.
                        QMutexLocker lkr(&mutex);

                        // Update the warning messages.
                        parse_warnings += local_errors;
                    }
                });

                // Fast-forward the line index.
                iline += (num_lines - 1);
            }
            else
            {
                for (int i=0; i<num_lines; ++i)
                {
                    // Parse a bond record line.
                    parse_line(lines()[iline+i], bonds, i,
                        num_lines, num_bonds, 2, parse_warnings);
                }
            }
        }

        // Angle records.
        else if (line.contains("!NTHETA"))
        {
            // Extract the number of bonds;
            bool ok;
            int num_angles = data.first().toInt(&ok);

            if (not ok)
            {
                parse_warnings.append(QObject::tr("Cannot extract number of angles "
                    "from part (%1) from line '%2'").arg(data.first()).arg(lines()[iline]));

                return;
            }

            // Resize data structures.
            angles.resize(num_angles);
            for (int i=0; i<num_angles; ++i)
                angles[i].resize(3);

            // Work out the number of record lines.
            // There are 3 angle records per line.
            const int num_lines = qCeil(num_angles/3.0);

            ++iline;

            if (usesParallel())
            {
                QMutex mutex;

                tbb::parallel_for(tbb::blocked_range<int>(0, num_lines),
                                [&](const tbb::blocked_range<int> &r)
                {
                    // Create local data objects.
                    QStringList local_errors;

                    for (int i=r.begin(); i<r.end(); ++i)
                    {
                        // Parse an angle record line.
                        parse_line(lines()[iline+i], angles, i,
                            num_lines, num_angles, 3, local_errors);
                    }

                    if (not local_errors.isEmpty())
                    {
                        // Acquire a lock.
                        QMutexLocker lkr(&mutex);

                        // Update the warning messages.
                        parse_warnings += local_errors;
                    }
                });

                // Fast-forward the line index.
                iline += (num_lines - 1);
            }
            else
            {
                for (int i=0; i<num_lines; ++i)
                {
                    // Parse an angle record line.
                    parse_line(lines()[iline+i], angles, i,
                        num_lines, num_angles, 3, parse_warnings);
                }
            }
        }

        // Dihedral records.
        else if (line.contains("!NPHI"))
        {
            // Extract the number of dihedrals;
            bool ok;
            int num_dihedrals = data.first().toInt(&ok);

            if (not ok)
            {
                parse_warnings.append(QObject::tr("Cannot extract number of dihedrals "
                    "from part (%1) from line '%2'").arg(data.first()).arg(lines()[iline]));

                return;
            }

            // Resize data structures.
            dihedrals.resize(num_dihedrals);
            for (int i=0; i<num_dihedrals; ++i)
                dihedrals[i].resize(4);

            // Work out the number of record lines.
            // There are 2 dihedral records per line.
            const int num_lines = qCeil(num_dihedrals/2.0);

            ++iline;

            if (usesParallel())
            {
                QMutex mutex;

                tbb::parallel_for(tbb::blocked_range<int>(0, num_lines),
                                [&](const tbb::blocked_range<int> &r)
                {
                    // Create local data objects.
                    QStringList local_errors;

                    for (int i=r.begin(); i<r.end(); ++i)
                    {
                        // Parse a dihedral record line.
                        parse_line(lines()[iline+i], dihedrals, i,
                            num_lines, num_dihedrals, 4, local_errors);
                    }

                    if (not local_errors.isEmpty())
                    {
                        // Acquire a lock.
                        QMutexLocker lkr(&mutex);

                        // Update the warning messages.
                        parse_warnings += local_errors;
                    }
                });

                // Fast-forward the line index.
                iline += (num_lines - 1);
            }
            else
            {
                for (int i=0; i<num_lines; ++i)
                {
                    // Parse a dihedral record line.
                    parse_line(lines()[iline+i], dihedrals, i,
                        num_lines, num_dihedrals, 4, parse_warnings);
                }
            }
        }

        // Improper records.
        else if (line.contains("!NIMPHI"))
        {
            // Extract the number of impropers;
            bool ok;
            int num_impropers = data.first().toInt(&ok);

            if (not ok)
            {
                parse_warnings.append(QObject::tr("Cannot extract number of impropers "
                    "from part (%1) from line '%2'").arg(data.first()).arg(lines()[iline]));

                return;
            }

            // Resize data structures.
            impropers.resize(num_impropers);
            for (int i=0; i<num_impropers; ++i)
                impropers[i].resize(4);

            // Work out the number of record lines.
            // There are 2 improper records per line.
            const int num_lines = qCeil(num_impropers/2.0);

            ++iline;

            if (usesParallel())
            {
                QMutex mutex;

                tbb::parallel_for(tbb::blocked_range<int>(0, num_lines),
                                [&](const tbb::blocked_range<int> &r)
                {
                    // Create local data objects.
                    QStringList local_errors;

                    for (int i=r.begin(); i<r.end(); ++i)
                    {
                        // Parse an improper record line.
                        parse_line(lines()[iline+i], impropers, i,
                            num_lines, num_impropers, 4, local_errors);
                    }

                    if (not local_errors.isEmpty())
                    {
                        // Acquire a lock.
                        QMutexLocker lkr(&mutex);

                        // Update the warning messages.
                        parse_warnings += local_errors;
                    }
                });

                // Fast-forward the line index.
                iline += (num_lines - 1);
            }
            else
            {
                for (int i=0; i<num_lines; ++i)
                {
                    // Parse an improper record line.
                    parse_line(lines()[iline+i], impropers, i,
                        num_lines, num_impropers, 4, parse_warnings);
                }
            }
        }

        // Cross-term records.
        else if (line.contains("!NCRTERM"))
        {
            // Extract the number of cross terms;
            bool ok;
            int num_cross = data.first().toInt(&ok);

            if (not ok)
            {
                parse_warnings.append(QObject::tr("Cannot extract number of cross terms "
                    "from part (%1) from line '%2'").arg(data.first()).arg(lines()[iline]));

                return;
            }

            // Resize data structures.
            cross_terms.resize(num_cross);
            for (int i=0; i<num_cross; ++i)
                cross_terms[i].resize(4);

            // Work out the number of record lines.
            // There are 2 cross term records per line.
            const int num_lines = qCeil(num_cross/2.0);

            ++iline;

            if (usesParallel())
            {
                QMutex mutex;

                tbb::parallel_for(tbb::blocked_range<int>(0, num_lines),
                                [&](const tbb::blocked_range<int> &r)
                {
                    // Create local data objects.
                    QStringList local_errors;

                    for (int i=r.begin(); i<r.end(); ++i)
                    {
                        // Parse a cross-term record line.
                        parse_line(lines()[iline+i], cross_terms, i,
                            num_lines, num_cross, 4, local_errors);
                    }

                    if (not local_errors.isEmpty())
                    {
                        // Acquire a lock.
                        QMutexLocker lkr(&mutex);

                        // Update the warning messages.
                        parse_warnings += local_errors;
                    }
                });

                // Fast-forward the line index.
                iline += (num_lines - 1);
            }
            else
            {
                for (int i=0; i<num_lines; ++i)
                {
                    // Parse a cross-term record line.
                    parse_line(lines()[iline+i], cross_terms, i,
                        num_lines, num_cross, 4, parse_warnings);
                }
            }
        }
    }

    // Use bonding information to break the parsed data into molecules.
    findMolecules();

    this->setScore(num_atoms);
}

/** Internal function that is used to parse a CHARMM parameter file.
        @param lines
            The lines from the parameter file.

        @param bond_params
            A multi-hash of the parsed bond parameters.

        @param angle_params
            A multi-hash of the parsed angle parameters.

        @param dihedral_params
            A multi-hash of the parsed dihedral parameters.

        @param improper_params
            A multi-hash of the parsed improper parameters.

        @param nonbonded_params
            A multi-hash of the parsed non-bonded parameters.

        @param box
            A periodic box object.

        @param has_box_params
            Whether the parameters contained box record information.
 */
void CharmmPSF::parseParameters(
    const QVector<QString> &parameter_lines,
    QMultiHash<QString, CharmmParam> &bond_params,
    QMultiHash<QString, CharmmParam> &angle_params,
    QMultiHash<QString, CharmmParam> &dihedral_params,
    QMultiHash<QString, CharmmParam> &improper_params,
    QMultiHash<QString, CharmmParam> &nonbonded_params,
    PeriodicBox &box, bool &has_box_params) const
{
    /* CHARMM parameter files are split in sections for different record types,
       i.e. bonds, angles, etc., separated by blank lines.

       The specifics of the formatting is dependent on the type of CHARMM
       forcefiled. Currently, we support the following types:

       1) As in CHARMM 22:
            - Record sections start with an upper-case record type identifier,
              e.g. BONDS.

       2) As in CHARMM 19:
            - Record sections have no starting identifier, but each record
              line in a section begins with a lower-case record indentifier,
              e.g. bond

       We assume that parameter files are formatted correctly, i.e. there
       is no overlap between record sections, and all sections must be
       separated by at least one blank line. (There are no blank lines within
       a record section.)

       We also parse NAMD XSC box data, which is added to the Sire System
       as a PeriodicBox Space property.
     */

    QStringList errors;

    for (int i=0; i<parameter_lines.count(); ++i)
    {
        // Extract the first word in the line.
        QStringList data = parameter_lines[i].simplified().split(QRegExp("\\s"));
        QString start = data[0];

        // Bond parameters.
        if (start == "BONDS")
        {
            bool is_blank_line = false;

            // Advance to the next line.
            i++;

            // Step-forward until we end the section.
            while (not is_blank_line)
            {
                QString line = parameter_lines[i].simplified();

                // Blank line.
                if (line.count() == 0)
                {
                    is_blank_line = true;
                    break;
                }

                // Not a comment.
                if (not (line[0] == '!'))
                {
                    CharmmParam param(line, 0, errors);
                    bond_params.insert(generateKey(param.getAtoms()), param);
                }

                i++;
            }
        }
        else if ((start == "bond") or (start == "BOND"))
        {
            data.removeFirst();
            CharmmParam param(data.join(" "), 0, errors);
            bond_params.insert(generateKey(param.getAtoms()), param);
        }

        // Angle parameters.
        else if (start == "ANGLES")
        {
            bool is_blank_line = false;

            // Advance to the next line.
            i++;

            // Step-forward until we end the section.
            while (not is_blank_line)
            {
                QString line = parameter_lines[i].simplified();

                // Blank line.
                if (line.count() == 0)
                {
                    is_blank_line = true;
                    break;
                }

                // Not a comment.
                if (not (line[0] == '!'))
                {
                    CharmmParam param(line, 1, errors);
                    angle_params.insert(generateKey(param.getAtoms()), param);
                }

                i++;
            }
        }
        else if ((start == "angle") or (start == "ANGLE"))
        {
            data.removeFirst();
            CharmmParam param(data.join(" "), 1, errors);
            angle_params.insert(generateKey(param.getAtoms()), param);
        }

        // Dihedral parameters.
        else if (start == "DIHEDRALS")
        {
            bool is_blank_line = false;

            // Advance to the next line.
            i++;

            // Step-forward until we end the section.
            while (not is_blank_line)
            {
                QString line = parameter_lines[i].simplified();

                // Blank line.
                if (line.count() == 0)
                {
                    is_blank_line = true;
                    break;
                }

                // Not a comment.
                if (not (line[0] == '!'))
                {
                    CharmmParam param(line, 2, errors);
                    dihedral_params.insert(generateKey(param.getAtoms()), param);
                }

                i++;
            }
        }
        else if ((start == "dihe") or (start == "DIHEDRAL"))
        {
            // Special handling for parameter records spanning multiple lines.
            if (data[5].toUpper() == "MULTIPLE=")
            {
                // Store the atom names.
                QStringList atoms = QStringList() << data[1] << data[2] << data[3] << data[4];

                // Process the first record.
                QStringList tmp = QStringList() << atoms << data[7] << data[8] << data[9];
                CharmmParam param(tmp.join(" "), 2, errors);
                dihedral_params.insert(generateKey(param.getAtoms()), param);

                // Now find all of the accompanying records.
                // We loop through the file until we ecounter the next dihedral record.
                while (true)
                {
                    // Move to the next line in the parameter file.
                    i++;

                    // Extract the first word in the line.
                    QStringList data = parameter_lines[i].simplified().split(QRegExp("\\s"));
                    start = data[0];

                    if (start == "DIHEDRAL")
                    {
                        // Go back a line and exit the loop.
                        i--;
                        break;
                    }

                    // Avoid blank lines.
                    if (data.count() > 0)
                    {
                        // This isn't a comment line.
                        if (data.first() != "!")
                        {
                            tmp = QStringList() << atoms << data[0] << data[1] << data[2];
                            CharmmParam param(tmp.join(" "), 2, errors);
                            dihedral_params.insert(generateKey(param.getAtoms()), param);
                        }
                    }
                }
            }

            // Standard, CHARMM-style dihedral records.
            else
            {
                data.removeFirst();
                CharmmParam param(data.join(" "), 2, errors);
                dihedral_params.insert(generateKey(param.getAtoms()), param);
            }
        }

        // Improper parameters.
        else if ((start == "IMPROPER") and (data.count() == 1))
        {
            bool is_blank_line = false;

            // Advance to the next line.
            i++;

            // Step-forward until we end the section.
            while (not is_blank_line)
            {
                QString line = parameter_lines[i].simplified();

                // Blank line.
                if (line.count() == 0)
                {
                    is_blank_line = true;
                    break;
                }

                // Not a comment.
                if (not (line[0] == '!'))
                {
                    CharmmParam param(line, 3, errors);
                    improper_params.insert(generateKey(param.getAtoms()), param);
                }

                i++;
            }
        }
        else if ((start == "impr") or (start == "IMPROPER"))
        {
            data.removeFirst();
            CharmmParam param(data.join(" "), 3, errors);
            improper_params.insert(generateKey(param.getAtoms()), param);
        }

        // Non-bonded parameters.
        else if (start.toUpper() == "NONBONDED")
        {
            // There must be more than one record on the line.
            if (data.count() > 1)
            {
                // This is a CHARMM format parameter file.
                if (data[1] == "nbxmod")
                {
                    bool is_blank_line = false;

                    // Skip two lines.
                    i += 2;

                    // Step-forward until we end the section.
                    while (not is_blank_line)
                    {
                        QString line = parameter_lines[i].simplified();

                        // Blank line.
                        if (line.count() == 0)
                        {
                            is_blank_line = true;
                            break;
                        }

                        // Not a comment.
                        if (not (line[0] == '!'))
                        {
                            CharmmParam param(line, 4, errors);
                            nonbonded_params.insert(generateKey(param.getAtoms()), param);
                        }

                        i++;
                    }
                }
                else
                {
                    data.removeFirst();
                    CharmmParam param(data.join(" "), 4, errors, true);
                    nonbonded_params.insert(generateKey(param.getAtoms()), param);
                }
            }
        }

        // Periodic box record data.
        // TODO:
        //  1) Figure out how to make this work if the file doesn't contain
        //     a $LABEL comment preceeding the cell record data.
        //
        //  2) How to handle boxes where the centre isn't at (0, 0, 0)?
        else if (start == "#$LABELS")
        {
            // Get the next line.
            i++;
            QStringList data = parameter_lines[i].simplified().split(QRegExp("\\s"));

            // Now try to read the cell parameters.
            // Here we are assuming a periodic, rectilinear simulation box.
            // TODO: Add support for other box types.
            bool ok_x, ok_y, ok_z;

            double x = data[1].toDouble(&ok_x);
            double y = data[5].toDouble(&ok_y);
            double z = data[9].toDouble(&ok_z);

            if (not ok_x or not ok_y or not ok_z)
            {
                errors.append(QObject::tr("Could not read NAMD XSC record! %1")
                    .arg(parameter_lines[i]));

                return;
            }

            else
            {
                // Create a periodic box object.
                box = PeriodicBox(Vector(x, y, z));
                has_box_params = true;
            }
        }
    }

    if (not errors.isEmpty())
    {
        throw SireIO::parse_error(QObject::tr("There were errors reading the CHARMM "
          "parameter file:\n%1").arg(errors.join("\n\n")), CODELOC);
    }
}

/** Internal function that is used to parameterise an existing molecule using
    the data from a CHARMM parameter file.
        @param sire_mol
            The molecule to parameterise.

        @param bond_params
            A multi-hash of the parsed bond parameters.

        @param angle_params
            A multi-hash of the parsed angle parameters.

        @param dihedral_params
            A multi-hash of the parsed dihedral parameters.

        @param improper_params
            A multi-hash of the parsed improper parameters.

        @param nonbonded_params
            A multi-hash of the non-bonded parameters.

        @return
            The parameterised molecule.
 */
SireMol::Molecule CharmmPSF::parameteriseMolecule(
    int imol,
    const SireMol::Molecule &sire_mol,
    const QMultiHash<QString, CharmmParam> &bond_params,
    const QMultiHash<QString, CharmmParam> &angle_params,
    const QMultiHash<QString, CharmmParam> &dihedral_params,
    const QMultiHash<QString, CharmmParam> &improper_params,
    const QMultiHash<QString, CharmmParam> &nonbonded_params,
    const PropertyMap &map) const
{
    // Get an editable version of the molecule.
    MolEditor edit_mol = sire_mol.edit();

    // Get the info object that can map between AtomNum to AtomIdx etc.
    const auto molinfo = sire_mol.info();

    // Check for "atomtype" property.
    if (not sire_mol.hasProperty(map["atomtype"]))
    {
        throw SireError::program_bug(QObject::tr("The molecule is missing "
            "property: \"atomtype\""), CODELOC);
    }

    // Get the atom type property object.
    const auto &atom_types = sire_mol.property(map["atomtype"]).asA<AtomStringProperty>();

    // Symbols for the various molecular potentials.
    const auto R     = InternalPotential::symbols().bond().r();
    const auto R_UB  = InternalPotential::symbols().ureyBradley().r();
    const auto Theta = InternalPotential::symbols().angle().theta();
    const auto Phi   = InternalPotential::symbols().dihedral().phi();

    // The molecule has bond connectivity information. Parameterise the bonds.
    if (sire_mol.hasProperty(map["connectivity"]))
    {
        // Get the atom connectivity object.
        const auto &connectivity = sire_mol.property(map["connectivity"])
                                           .asA<ConnectivityEditor>();

        // Initialise the bond parameter object.
        TwoAtomFunctions bond_funcs(edit_mol);

        // Loop over all of the bonds.
        for (const auto &bond : connectivity.getBonds())
        {
            // Extract the cgAtomIdx for the two atoms.
            const auto idx0 = molinfo.cgAtomIdx(bond.atom0());
            const auto idx1 = molinfo.cgAtomIdx(bond.atom1());

            // Now get the type for each atom.
            const auto atom0 = atom_types[idx0];
            const auto atom1 = atom_types[idx1];

            // Create a vector of the atom types.
            QVector<QString> bond_atoms({atom0, atom1});

            // Find the parameters for the bond.
            auto matches = findParameters(bond_atoms, bond_params, 0);

            // No matches!
            if (matches.count() == 0)
            {
                throw SireError::incompatible_error(QObject::tr("Missing bond parameters "
                    "for atom types \'%1\' and \'%2\'").arg(atom0).arg(atom1), CODELOC);
            }

            // Get the bond parameters.
            auto params = matches[0].getParams();

            // Create the expression for the bond function.
            Expression func = params[0] * SireMaths::pow_2( R - params[1] );

            // Set the bond function parameter.
            bond_funcs.set(idx0, idx1, func);
        }

        // Set the bond property.
        edit_mol.setProperty(map["bond"], bond_funcs);
    }

    // Parameterise the angles.
    if (nAngles(imol) > 0)
    {
        // Initialise the angle parameter object.
        ThreeAtomFunctions angle_funcs(edit_mol);

        // Intialise Urey-Bradley parameter object.
        TwoAtomFunctions ub_funcs(edit_mol);

        // Whether Urey-Bradley terms are present.
        bool has_ub = false;

        // Loop over all of the angles.
        for (int i=0; i<nAngles(imol); ++i)
        {
            // Get the angle index.
            const int idx = mol_angles[imol][i];

            // Get the indices of the three atoms.
            const int idx0 = num_to_idx[angles[idx][0]];
            const int idx1 = num_to_idx[angles[idx][1]];
            const int idx2 = num_to_idx[angles[idx][2]];

            // Create a vector of the atom types.
            QVector<QString> angle_atoms({atoms[idx0].getType(),
                                          atoms[idx1].getType(),
                                          atoms[idx2].getType()});

            // Find the parameters for the angle.
            auto matches = findParameters(angle_atoms, angle_params, 1);

            // No matches!
            if (matches.count() == 0)
            {
                throw SireError::incompatible_error(QObject::tr("Missing angle parameters "
                    "for atom types \'%1\', \'%2\', and \'%3\'")
                    .arg(angle_atoms[0]).arg(angle_atoms[1]).arg(angle_atoms[2]), CODELOC);
            }

            // Get the angle parameters.
            auto params = matches[0].getParams();

            // Create the expression for the angle function, converting the angle to radians.
            Expression func = params[0] * SireMaths::pow_2( Theta - qDegreesToRadians(params[1]) );

            // Set the angle function parameter.
            angle_funcs.set(AtomNum(atoms[idx0].getNumber()),
                            AtomNum(atoms[idx1].getNumber()),
                            AtomNum(atoms[idx2].getNumber()),
                            func);

            // Add a Urey-Bradley bond function.
            if (params.count() > 2)
            {
                func = params[2] * SireMaths::pow_2( R_UB - params[3] );

                // Set the Urey-Bradley function parameter.
                ub_funcs.set(AtomNum(atoms[idx0].getNumber()),
                             AtomNum(atoms[idx2].getNumber()),
                             func);

                // Flag that we've found Urey-Bradley parameters.
                has_ub = true;
            }
        }

        // Set the angle property...
        edit_mol.setProperty(map["angle"], angle_funcs);

        // ...and Urey-Bradley, if present.
        if (has_ub)
        {
            edit_mol.setProperty(map["urey-bradley"], ub_funcs);
        }
    }

    // Parameterise the dihedrals.
    if (nDihedrals(imol) > 0)
    {
        // Initialise the dihedral parameter object.
        FourAtomFunctions dihedral_funcs(edit_mol);

        // Loop over all of the dihedrals.
        for (int i=0; i<nDihedrals(imol); ++i)
        {
            // Get the dihedral index.
            const int idx = mol_dihedrals[imol][i];

            // Get the indices of the four atoms.
            const int idx0 = num_to_idx[dihedrals[idx][0]];
            const int idx1 = num_to_idx[dihedrals[idx][1]];
            const int idx2 = num_to_idx[dihedrals[idx][2]];
            const int idx3 = num_to_idx[dihedrals[idx][3]];

            // Create a vector of the atom types.
            QVector<QString> dihedral_atoms({atoms[idx0].getType(),
                                             atoms[idx1].getType(),
                                             atoms[idx2].getType(),
                                             atoms[idx3].getType()});

            // Find the parameters for the dihedral.
            auto matches = findParameters(dihedral_atoms, dihedral_params, 2);

            // No matches!
            if (matches.count() == 0)
            {
                throw SireError::incompatible_error(QObject::tr("Missing dihedral parameters "
                    "for atom types \'%1\', \'%2\', \'%3\', and \'%4\'")
                    .arg(dihedral_atoms[0]).arg(dihedral_atoms[1]).arg(dihedral_atoms[2])
                    .arg(dihedral_atoms[3]), CODELOC);
            }

            // Intialise the function object.
            Expression func;

            // Loop over all matches.
            // Potentially mutliple multiplicity values.
            for (const auto &match : matches)
            {
                // Get the dihedral parameters.
                auto params = match.getParams();

                // Update the function, converting the phase shift to radians.
                func += params[0] * (1 + Cos(( params[1] * Phi ) - qDegreesToRadians(params[2]) ));
            }

            // Set the dihedral function parameter.
            dihedral_funcs.set(AtomNum(atoms[idx0].getNumber()),
                               AtomNum(atoms[idx1].getNumber()),
                               AtomNum(atoms[idx2].getNumber()),
                               AtomNum(atoms[idx3].getNumber()),
                               func);
        }

        // Set the dihedral property.
        edit_mol.setProperty(map["dihedral"], dihedral_funcs);
    }

    // Parameterise the impropers.
    if (nImpropers(imol) > 0)
    {
        // Initialise the improper parameter object.
        FourAtomFunctions improper_funcs(edit_mol);

        // Loop over all of the impropers.
        for (int i=0; i<nImpropers(imol); ++i)
        {
            // Get the improper index.
            const int idx = mol_impropers[imol][i];

            // Get the indices of the four atoms.
            const int idx0 = num_to_idx[impropers[idx][0]];
            const int idx1 = num_to_idx[impropers[idx][1]];
            const int idx2 = num_to_idx[impropers[idx][2]];
            const int idx3 = num_to_idx[impropers[idx][3]];

            // Create a vector of the atom types.
            QVector<QString> improper_atoms({atoms[idx0].getType(),
                                             atoms[idx1].getType(),
                                             atoms[idx2].getType(),
                                             atoms[idx3].getType()});

            // Find the parameters for the improper.
            auto matches = findParameters(improper_atoms, improper_params, 3);

            // No matches!
            if (matches.count() == 0)
            {
                throw SireError::incompatible_error(QObject::tr("Missing improper parameters "
                    "for atom types \'%1\', \'%2\', \'%3\', and \'%4\'")
                    .arg(improper_atoms[0]).arg(improper_atoms[1]).arg(improper_atoms[2])
                    .arg(improper_atoms[3]), CODELOC);
            }

            // Get the improper parameters.
            auto params = matches[0].getParams();

            // Intialise the function object, converting the out of plane angle to radians.
            Expression func = params[0] * SireMaths::pow_2( Phi - qDegreesToRadians(params[1]) );

            // Set the improper function parameter.
            improper_funcs.set(AtomNum(atoms[idx0].getNumber()),
                               AtomNum(atoms[idx1].getNumber()),
                               AtomNum(atoms[idx2].getNumber()),
                               AtomNum(atoms[idx3].getNumber()),
                               func);

        }

        // Set the improper property.
        edit_mol.setProperty(map["improper"], improper_funcs);
    }

    // Parameterise the cross-terms.
    if (nCrossTerms(imol) > 0)
    {
        // TODO: Work out what to do here...
    }

    return edit_mol.commit();
}

/** Internal function to generate CHARMM PSF record data from a Sire molecule

    @param imol
        The molecule index.

    @param sire_mol
        A reference to the Sire molecule.

    @param local_atoms
        The PSFAtom vector local to the molecule.

    @param local_bonds
        The bond records vector local to the molecule.

    @param local_angles
        The angle records vector local to the molecule.

    @param local_dihedrals
        The dihedral records vector local to the molecule.

    @param local_impropers
        The improper records vector local to the molecule.

    @param bond_params
        The set of CHARMM bond parameter strings.

    @param angle_params
        The set of CHARMM angle parameter strings.

    @param dihedral_params
        The set of CHARMM dihedral parameter strings.

    @param improper_params
        The set of CHARMM improper parameter strings.

    @param local_errors
        A list of error messages local to the molecule.

    @param map
        The map of user properties.
 */
void CharmmPSF::parseMolecule(
        int imol,
        const SireMol::Molecule &sire_mol,
        QVector<PSFAtom> &local_atoms,
        QVector<QVector<qint64> > &local_bonds,
        QVector<QVector<qint64> > &local_angles,
        QVector<QVector<qint64> > &local_dihedrals,
        QVector<QVector<qint64> > &local_impropers,
        QSet<QString> &bond_params,
        QSet<QString> &angle_params,
        QSet<QString> &dihedral_params,
        QSet<QString> &improper_params,
        QStringList &local_errors,
        const PropertyMap &map)
{
    // Store the number of atoms in the molecule.
    int num_atoms = sire_mol.nAtoms();

    // Early exit.
    if (num_atoms == 0) return;

    // Extract the molecule data.
    const auto moldata = sire_mol.data();

    // Resize the atoms vector.
    local_atoms.resize(num_atoms);

    // Set a default molecule name (used for a "segment" ID in the PSF file).
    QString segment = QString("M%1").arg(imol+1);

    // Extract the existing molecule name.
    if (sire_mol.hasProperty(map["mol-name"]))
    {
        segment = sire_mol.property(map["mol-name"]).toString().simplified().toUpper();
    }

    if (usesParallel())
    {
        QMutex mutex;

        tbb::parallel_for(tbb::blocked_range<int>(0, num_atoms),
                           [&](const tbb::blocked_range<int> r)
        {
            // Create local data objects.
            QStringList local_errors;

            for (int i=r.begin(); i<r.end(); ++i)
            {
                // Parse atom data.
                local_atoms[i] = PSFAtom(sire_mol.atom(AtomIdx(i)), segment, local_errors, map);
            }

            if (not local_errors.isEmpty())
            {
                // Acquire a lock.
                QMutexLocker lkr(&mutex);

                // Update global error data.
                parse_warnings += local_errors;
            }
        });
    }
    else
    {
        for (int i=0; i<num_atoms; ++i)
        {
            // Parse atom data.
            local_atoms[i] = PSFAtom(sire_mol.atom(AtomIdx(i)), segment, local_errors, map);
        }
    }

    bool has_bonds, has_angles, has_ubs, has_dihedrals, has_impropers;

    // Get the required properties.
    const auto bond_funcs = getProperty<TwoAtomFunctions>(map["bond"], moldata, &has_bonds);
    const auto angle_funcs = getProperty<ThreeAtomFunctions>(map["angle"], moldata, &has_angles);
    const auto ub_funcs = getProperty<TwoAtomFunctions>(map["urey-bradley"], moldata, &has_ubs);
    const auto dihedral_funcs = getProperty<FourAtomFunctions>(map["dihedral"], moldata, &has_dihedrals);
    const auto improper_funcs = getProperty<FourAtomFunctions>(map["improper"], moldata, &has_impropers);

    // Now construct the PSF record data...

    // Bonds.
    if (has_bonds)
        getBondsFrom(bond_funcs, sire_mol, local_bonds, bond_params, map);

    // Angles.
    if (has_angles)
        getAnglesFrom(angle_funcs, ub_funcs, sire_mol, local_angles, angle_params, map);

    // Dihedrals.
    if (has_dihedrals)
        getFourAtomFrom(dihedral_funcs, sire_mol, local_dihedrals, dihedral_params, map);

    // Impropers.
    if (has_impropers)
        getFourAtomFrom(improper_funcs, sire_mol, local_impropers, improper_params, map, true);
}

/** Find the index of the parameters associated with the 'search_atoms'.
        @param search_atoms
            The list of atoms to parameterise.

        @param params
            A multi-hash of the parameterisation records.

        @param type
            The type of parameterisation record.

        @return
            A vector of the indices of the matching parameter sets.
  */
QList<CharmmParam> CharmmPSF::findParameters(const QVector<QString> &search_atoms,
    const QMultiHash<QString, CharmmParam> &params, int type) const
{
    // Bond params.
    if (type == 0)
    {
        // Make sure that there are two atoms to search for.
        if (search_atoms.count() != 2)
        {
            throw SireError::program_bug(QObject::tr("Cannot search for CHARMM "
                "bond parameters, incorrect number of atoms passed. Given "
                "%1, expected 2").arg(search_atoms.count()), CODELOC);
        }

        // Generate the key for this set of atoms.
        QString key = generateKey(search_atoms);

        // Return the matches.
        return params.values(key);
    }

    // Angle params.
    else if (type == 1)
    {
        // Make sure that there are three atoms to search for.
        if (search_atoms.count() != 3)
        {
            throw SireError::program_bug(QObject::tr("Cannot search for CHARMM "
                "angle parameters, incorrect number of atoms passed. Given "
                "%1, expected 3").arg(search_atoms.count()), CODELOC);
        }

        // Generate the key for this set of atoms.
        QString key = generateKey(search_atoms);

        // Return the matches.
        return params.values(key);
    }

    // Dihedral params.
    else if ((type == 2) or (type == 3))
    {
        // Make sure that there are four atoms to search for.
        if (search_atoms.count() != 4)
        {
            throw SireError::program_bug(QObject::tr("Cannot search for CHARMM "
                "%1 parameters, incorrect number of atoms passed. Given "
                "%2, expected 4")
                .arg((type == 2) ? "dihedral" : "improper")
                .arg(search_atoms.count()), CODELOC);
        }

        // Generate the key for this set of atoms.
        QString key = generateKey(search_atoms);

        // Check for matches.
        auto matches = params.values(key);

        // No matches, try wildcards.
        if (matches.count() == 0)
        {
            // Backup the original vector of atom types.
            auto copy_atoms = search_atoms;

            // Single wildcard.
            for (int i=0; i<4; ++i)
            {
                // Replace the atom type with a wildcard.
                copy_atoms[i] = "X";

                // Generate the key for this set of atoms.
                QString key = generateKey(copy_atoms);

                // Check for matches.
                matches = params.values(key);

                // Return the matches.
                if (matches.count() > 0)
                    return matches;

                // Reset the atom type.
                copy_atoms[i] = search_atoms[i];
            }

            // Double wildcard.
            for (int i=0; i<4; ++i)
            {
                // Replace atom type 'i' with a wildcard.
                copy_atoms[i] = "X";

                for (int j=i+i; j<4; ++j)
                {
                    // Replace atom type 'j' with a wildcard.
                    copy_atoms[j] = "X";

                    // Generate the key for this set of atoms.
                    QString key = generateKey(copy_atoms);

                    // Check for matches.
                    matches = params.values(key);

                    // Return the matches.
                    if (matches.count() > 0)
                        return matches;

                    // Reset the type for atom 'j'.
                    copy_atoms[j] = search_atoms[j];
                }

                // Reset the type for atom 'i'.
                copy_atoms[i] = search_atoms[i];
            }

            // Triple wildcard.
            for (int i=0; i<4; ++i)
            {
                // Replace atom type 'i' with a wildcard.
                copy_atoms[i] = "X";

                for (int j=i+i; j<4; ++j)
                {
                    // Replace atom type 'j' with a wildcard.
                    copy_atoms[j] = "X";

                    for (int k=j+1; k<4; ++k)
                    {
                        // Replace atom type 'k' with a wildcard.
                        copy_atoms[k] = "X";

                        // Generate the key for this set of atoms.
                        QString key = generateKey(copy_atoms);

                        // Check for matches.
                        matches = params.values(key);

                        // Return the matches.
                        if (matches.count() > 0)
                            return matches;

                        // Reset the type for atom 'k'.
                        copy_atoms[k] = search_atoms[k];
                    }

                    // Reset the type for atom 'j'.
                    copy_atoms[j] = search_atoms[j];
                }

                // Reset the type for atom 'i'.
                copy_atoms[i] = search_atoms[i];
            }

            // Quadruple wildcard.
            // Unlikley, but what the heck!

            // Check for matches.
            matches = params.values("X;X;X;X");

            // If we've got this far, then there are no matches.
            // Return the empty list.
            return matches;
        }
        else return matches;
    }

    else if (type == 4)
    {
        // TODO: Not sure what these parameter records look like. Find an example
        //       In the CHARMM forcefield files.

        return QList<CharmmParam>();
    }

    // Uknown parameter type.
    else
    {
        throw SireError::program_bug(QObject::tr("Unknown parameter type (%1). "
            "Valid types are 0, 1, 2, 3, 4").arg(type), CODELOC);
    }
}

/** Generate a key from the vector of words.
        @param words
            The list of words.

        @return
            The key.
 */
QString CharmmPSF::generateKey(QVector<QString> words) const
{
    if (words.count() == 0)
        return QString();

    // Sort the words.
    qSort(words);

    // Initialise the key string.
    QString key(words[0]);

    // Append each word to the key, semicolon separated.
    for (int i=1; i<words.count(); ++i)
        key += QString(";%1").arg(words[i]);

    return key;
}

/** Use the data contained in this parser to create a new System of molecules,
    assigning properties based on the mapping in 'map' */
System CharmmPSF::startSystem(const PropertyMap &map) const
{
    // This assumes that we've already run 'findMolecules()' to
    // break the parsed data into molecules using the bond records.

    const int nmols = nMolecules();

    if (nmols == 0)
        return System();

    QVector<Molecule> mols(nmols);
    Molecule *mols_array = mols.data();

    if (usesParallel())
    {
        tbb::parallel_for(tbb::blocked_range<int>(0, nmols),
                           [&](tbb::blocked_range<int> r)
        {
            for (int i=r.begin(); i<r.end(); ++i)
            {
                mols_array[i] = this->getMolecule(i, map);
            }
        });
    }
    else
    {
        for (int i=0; i<nmols; ++i)
        {
            mols_array[i] = this->getMolecule(i, map);
        }
    }

    MoleculeGroup molgroup("all");

    for (auto mol : mols)
    {
        molgroup.add(mol);
    }

    System system;
    system.add(molgroup);
    system.setProperty(map["fileformat"].source(), StringProperty(this->formatName()));

    return system;
}

/** Use the data contained in this parser to create a new System of molecules,
    using the supplementary data records from 'lines' and assigning properties
    based on the mapping in 'map' */
System CharmmPSF::startSystem(const QVector<QString> &param_lines, const PropertyMap &map) const
{
    // This assumes that we've already run 'findMolecules()' to
    // break the parsed data into molecules using the bond records.

    // Potential parameter file data is contained in the 'param_lines' vector.
    // This is parsed separately by 'parseParameters()' before the system
    // can be constructed.

    // Generate the unparameterised system.
    System system = startSystem(map);

    // Initialise the parameter objects.
    QMultiHash<QString, CharmmParam> bond_params;
    QMultiHash<QString, CharmmParam> angle_params;
    QMultiHash<QString, CharmmParam> dihedral_params;
    QMultiHash<QString, CharmmParam> improper_params;
    QMultiHash<QString, CharmmParam> nonbonded_params;

    // Initialise a periodic box object.
    PeriodicBox box;
    bool has_box_params = false;

    // Parse and validate the parameter file.
    parseParameters(param_lines, bond_params, angle_params,
        dihedral_params, improper_params, nonbonded_params, box, has_box_params);

    // Early exit if parameters are missing.
    if ((nBonds()      > 0 and bond_params.count()     == 0) or
        (nAngles()     > 0 and angle_params.count()    == 0) or
        (nDihedrals()  > 0 and dihedral_params.count() == 0) or
        (nImpropers()  > 0 and improper_params.count() == 0))
    {
        throw SireError::incompatible_error(QObject::tr("The parameter file "
            "is missing information that is required for generating a system!"));
    }

    // Get the MolNums of each molecule in the System - this returns the
    // numbers in MolIdx order.
    const QVector<MolNum> molnums = system.getMoleculeNumbers().toVector();

    // Store the number of molecules.
    const int num_mols = molnums.count();

    // No molecules in the system.
    if (num_mols == 0)
        return system;

    // A vector of the updated molecules.
    QVector<Molecule> mols(num_mols);

    // Parameterise each molecule in the system.
    if (usesParallel())
    {
        tbb::parallel_for(tbb::blocked_range<int>(0, num_mols),
                           [&](const tbb::blocked_range<int> r)
        {
            for (int i=r.begin(); i<r.end(); ++i)
            {
                // Parameterise the molecule.
                mols[i] = parameteriseMolecule(i, system[molnums[i]].molecule(),
                    bond_params, angle_params, dihedral_params, improper_params,
                    nonbonded_params, map);
            }
        });
    }
    else
    {
        for (int i=0; i<num_mols; ++i)
        {
            // Parameterise the molecule.
            mols[i] = parameteriseMolecule(i, system[molnums[i]].molecule(),
                bond_params, angle_params, dihedral_params, improper_params,
                nonbonded_params, map);
        }
    }

    // Update the system.
    system.update(Molecules(mols));

    // Convert the vector of params to a single string.
    QString params_string;
	for (const auto &line : param_lines)
		params_string.append(line + "\n");

    // Add the CHARMM parameters as a property.
    system.setProperty(map["charmm-params"].source(), StringProperty(params_string));

    // Add the periodic box property.
    if (has_box_params)
        system.setProperty(map["space"].source(), box);

    // Return the parameterised system.
    return system;
}

/** Write the parsed data back to the file called 'filename'. This will
    overwrite the file if it exists already, so be careful! */
void CharmmPSF::writeToFile(const QString &filename) const
{
    QVector<QString> lines = toLines();

    if (lines.isEmpty())
        return;

    QFile f(filename);

    if (not f.open(QIODevice::WriteOnly | QIODevice::Text))
    {
        throw SireError::file_error(f, CODELOC);
    }

    // Write the PSF file.

    QTextStream ts(&f);

    for (const QString &line : lines)
    {
        ts << line << '\n';
    }

    f.close();

    // Write supplementary CHARMM parameters.

    if (not charmm_params.isEmpty())
    {
        QFileInfo fi(filename);

        // Create the name of the parameter file.
        QString param_filename = fi.absolutePath()
                               + '/'
                               + fi.completeBaseName()
                               + ".params";

        QFile f(param_filename);

        if (not f.open(QIODevice::WriteOnly | QIODevice::Text))
        {
            throw SireError::file_error(f, CODELOC);
        }

        QTextStream ts(&f);

        for (const QString &line : charmm_params)
        {
            ts << line << '\n';
        }

        f.close();
    }

    // Write periodic box information as a NAMD XSC file.

    if (has_box)
    {
        QFileInfo fi(filename);

        // Create the name of the XSC file.
        QString xsc_filename = fi.absolutePath()
                             + '/'
                             + fi.completeBaseName()
                             + ".xsc";

        QFile f(xsc_filename);

        if (not f.open(QIODevice::WriteOnly | QIODevice::Text))
        {
            throw SireError::file_error(f, CODELOC);
        }

        // Get the box dimensions.
        Vector box_size = box.dimensions();

        QTextStream ts(&f);

        QString l1 = "# NAMD extended system configuration output file\n";
        QString l2 = "#$LABELS step a_x a_y a_z b_x b_y b_z c_x c_y c_z\n";
        QString l3 = QString("0 %1 0 0 0 %2 0 0 0 %3 0 0 0\n")
            .arg(box_size[0]).arg(box_size[1]).arg(box_size[2]);

        // Write the box data to file.
        ts << l1 << l2 << l3;

        f.close();
    }
}

/** Internal function used to get the molecule structure for molecule 'imol'. */
MolStructureEditor CharmmPSF::getMolStructure(int imol, const PropertyName &cutting) const
{
    // Make sure the frame index is within range.
    if ((imol < 0) or (imol > molecules.count()))
    {
        throw SireError::program_bug(QObject::tr("The molecule index %1 is out of "
            "range, 0 - %2").arg(imol).arg(molecules.count()), CODELOC);
    }

    // Make sure that there are atoms in the molecule.
    if (molecules[imol].count() == 0)
    {
        throw SireError::program_bug(QObject::tr(
            "Strange - there are no atoms in molecule %1?")
                .arg(imol), CODELOC);
    }

    // First step is to build the structure of the molecule, i.e.
    // the layout of cutgroups, residues and atoms.
    MolStructureEditor mol;

    /* To do this we'll walk through all of the atoms in the molecule to
       work out what residues exist. Then we can add the residues to
       the molecules, then add the atoms, reparenting them to their residue.

       Note that PSF does not contain chain information.
     */

    // Mapping between residues and atoms.
    // There will be multiple atoms per residue.
    QMultiMap<int, int> res_to_atom;

    // Mapping between residue number and name.
    QMap<int, QString> num_to_name;

    // Loop through all the atoms in the molecule.
    for (int i=0; i<nAtoms(imol); ++i)
    {
        // Store the current atom ID.
        int atom_id = molecules[imol][i];

        // Store a reference to the current atom.
        const PSFAtom &atom = atoms[atom_id];

        // Store the atom number.
        const int res_num = atom.getResNum();

        // Map the atom to its residue.
        res_to_atom.insert(res_num, atom_id);

        // Make sure the residue number doesn't already exist
        // in the map. Residue numbers must be unique.
        if (not num_to_name.contains(res_num))
        {
            // Map the residue number to its name.
            num_to_name.insert(res_num, atom.getResName());
        }
        else
        {
            // This residue number has already appeared with a different name.
            if (atom.getResName() != num_to_name[res_num])
            {
                throw SireError::incompatible_error(QObject::tr("The parsed data "
                    "contains non-unique residue numbers! Residue number %1 "
                    "for atom number %2 is named %3, was previously assigned name %4")
                    .arg(res_num).arg(atom.getNumber()).arg(atom.getResName())
                    .arg(num_to_name[res_num]), CODELOC);
            }
        }
    }

    // Residue index.
    int ires = 0;

    // Loop over all unique residues in the molecule.
    for (auto res_num : res_to_atom.uniqueKeys())
    {
        // Get the residue name.
        QString res_name = num_to_name[res_num];

        // By default we will use one CutGroup per residue.
        // This may be changed later by the cutting system.
        auto cutgroup = mol.add(CGName(QString::number(ires)));

        // Get a sorted list of the atoms in the residue.
        QList<int> res_atoms = res_to_atom.values(res_num);
        qSort(res_atoms);

        // Add the residue to the molecule.
        auto res = mol.add(ResNum(res_num));
        res.rename(ResName(res_name.trimmed()));

        // Add each atom in the residue to the molecule.
        for (auto res_atom : res_atoms)
        {
            auto atom = cutgroup.add(AtomNum(atoms[res_atom].getNumber()));
            atom.rename(AtomName(atoms[res_atom].getName().trimmed()));

            // Reparent the atom to its residue.
            atom.reparent(ResNum(res_num));
        }

        ires++;
    }

    if (cutting.hasValue())
    {
        const CuttingFunction &cutfunc = cutting.value().asA<CuttingFunction>();

        if (not cutfunc.isA<ResidueCutting>())
        {
            mol = cutfunc(mol);
        }
    }

    return mol;
}

/** Internal function used to get the molecule structure for molecule 'imol'. */
MolEditor CharmmPSF::getMolecule(int imol, const PropertyMap &map) const
{
    // Make sure the molecule (model) index is within range.
    if ((imol < 0) or (imol > molecules.count()))
    {
        throw SireError::program_bug(QObject::tr("The frame index %1 is out of "
            "range, 0 - %2").arg(imol).arg(molecules.count()), CODELOC);
    }

    // Make sure that there are atoms in the frame.
    if (nAtoms(imol) == 0)
    {
        throw SireError::program_bug(QObject::tr(
            "Strange - there are no atoms in molecule %1?")
                .arg(imol), CODELOC);
    }

    // First, construct the layout of the molecule (sorting of atoms into residues and cutgroups).
    auto mol = this->getMolStructure(imol, map["cutting"]).commit().edit();

    // Set the molecule name.
    mol.setProperty(map["mol-name"], StringProperty(atoms[molecules[imol][0]].getSegment()));

    // Get the info object that can map between AtomNum to AtomIdx etc.
    const auto molinfo = mol.info();

    // Instantiate the atom property objects that we need.
    AtomStringProperty types(molinfo);
    AtomCharges        charges(molinfo);
    AtomMasses         masses(molinfo);

    // Now loop through the atoms in the molecule and set each property.
    for (int i=0; i<nAtoms(imol); ++i)
    {
        // Store a reference to the current atom.
        const auto &atom = atoms[molecules[imol][i]];

        // Determine the CGAtomIdx for this atom.
        auto cgatomidx = molinfo.cgAtomIdx(AtomNum(atom.getNumber()));

        // Set the properties.
        types.set(cgatomidx, atom.getType());
        masses.set(cgatomidx, atom.getMass() * SireUnits::g_per_mol);
        charges.set(cgatomidx, double(atom.getCharge()) * SireUnits::mod_electron);
    }

    // Now create the connectivity object for the molecular system. */

    // Connectivity object for bonded atoms.
    ConnectivityEditor connectivity(molinfo);

    // Loop over all bonds in the molecule.
    for (const auto &bond : mol_bonds[imol])
    {
        // Get the atom indices for the two bonds.
        int idx1 = num_to_idx[bonds[bond][0]];
        int idx2 = num_to_idx[bonds[bond][1]];

        // Add the bond to the connectivity object.
        connectivity.connect(AtomNum(atoms[idx1].getNumber()),
                             AtomNum(atoms[idx2].getNumber()));
    }

    return mol.setProperty(map["atomtype"], types)
              .setProperty(map["charge"], charges)
              .setProperty(map["mass"], masses)
              .setProperty(map["connectivity"], connectivity)
              .commit();
}

/** Internal function to break the parsed PSF record data into molecules. */
void CharmmPSF::findMolecules()
{
    // Clear any existing molecule data.
    molecules.clear();

    // Create a hash of the bonded atoms.
    QHash<int, int> bonded_atoms;

    for (int i=0; i<nBonds(); ++i)
    {
        bonded_atoms.insertMulti(bonds[i][0], bonds[i][1]);
        bonded_atoms.insertMulti(bonds[i][1], bonds[i][0]);
    }

    // Now recursively walk along each atom to find all the atoms that
    // are in the same molecule.

    // The number of molecules that are found.
    int nmols = 0;

    // A hash between atom and molecule indices.
    QHash<int, int> atom_to_mol;

    // Create the hash between atom number and index.
    // Clear any existing hash first.
    num_to_idx.clear();
    for(int i=0; i<nAtoms(); ++i)
        num_to_idx.insert(atoms[i].getNumber(), i);

    /*************************************************/
    /* Work out which atoms belong to this molecule. */
    /*************************************************/

    // Loop over all atoms by index.
    for (int i = 0; i<nAtoms(); ++i)
    {
        // Get the atom number.
        int num = atoms[i].getNumber();

        // The molecule doesn't already contain this atom.
        if (not atom_to_mol.contains(num))
        {
            // Initialise a set for atoms in this molecule.
            QSet<qint64> atoms_in_mol;

            nmols += 1;
            atom_to_mol[i] = nmols;
            atoms_in_mol.insert(num);

            // Recursive walk from this atom.
            findBondedAtoms(num, nmols, bonded_atoms, atom_to_mol, atoms_in_mol);

            // We've now found all of the atoms in this molecule!

            // Convert to a vector of atom numbers.
            QVector<qint64> mol_atoms = atoms_in_mol.toList().toVector();

            // Now convert the atom numbers to indices in the atoms vector
            // and set the molecule index for each atom.
            for (auto &atom : mol_atoms)
            {
                atom = num_to_idx[atom];
                atoms[atom].setMolIndex(nmols-1);
            }

            // Add the sorted atom indices.
            qSort(mol_atoms);
            molecules.append(mol_atoms);
        }
    }

    // Store the number of molecules.
    const auto num_mols = nMolecules();

    /*************************************************/
    /* Work out which bonds belong to this molecule. */
    /*************************************************/

    mol_bonds.resize(num_mols);

    // Loop over all of the bonds.
    for (int i=0; i<nBonds(); ++i)
    {
        // Get the molecule index.
        const int molidx = atoms[num_to_idx[bonds[i][0]]].getMolIndex();

        // Make sure the terminal atom is also in the molecule.
        if (atoms[num_to_idx[bonds[i][1]]].getMolIndex() == molidx)
        {
            mol_bonds[molidx].append(i);
        }
        else
        {
            throw SireError::program_bug(QObject::tr("The bonded atoms "
                "are not in the same atom: [ AtomNum(%1), MolIdx(%2) ] "
                "and [ AtomNum(%3), MolIdx(%4) ]")
                .arg(bonds[i][0]).arg(atoms[num_to_idx[bonds[i][0]]].getMolIndex())
                .arg(bonds[i][1]).arg(atoms[num_to_idx[bonds[i][1]]].getMolIndex()), CODELOC);
        }
    }

    /***************************************************/
    /*  Work out which angles belong to this molecule. */
    /***************************************************/

    mol_angles.resize(num_mols);

    // Loop over all of the angles.
    for (int i=0; i<nAngles(); ++i)
    {
        // Get the molecule index.
        const int molidx = atoms[num_to_idx[angles[i][0]]].getMolIndex();

        // Make sure the terminal atom is also in the molecule.
        if (atoms[num_to_idx[angles[i][2]]].getMolIndex() == molidx)
        {
            mol_angles[molidx].append(i);
        }
        else
        {
            throw SireError::program_bug(QObject::tr("The atoms involved in an "
                "angle interaction are not in the same molecule: [ AtomNum(%1), MolIdx(%2) ] "
                "and [ AtomNum(%3), MolIdx(%4) ]")
                .arg(angles[i][0]).arg(atoms[num_to_idx[angles[i][0]]].getMolIndex())
                .arg(angles[i][2]).arg(atoms[num_to_idx[angles[i][2]]].getMolIndex()), CODELOC);
        }
    }

    /*****************************************************/
    /* Work out which dihedrals belong to this molecule. */
    /*****************************************************/

    mol_dihedrals.resize(num_mols);

    // Loop over all of the dihedrals.
    for (int i=0; i<nDihedrals(); ++i)
    {
        // Get the molecule index.
        int molidx = atoms[num_to_idx[dihedrals[i][0]]].getMolIndex();

        // Make sure the terminal atom is also in the molecule.
        if (atoms[num_to_idx[dihedrals[i][3]]].getMolIndex() == molidx)
        {
            mol_dihedrals[molidx].append(i);
        }
        else
        {
            throw SireError::program_bug(QObject::tr("The atoms involved in a "
                "dihedral interaction are not in the same molecule: [ AtomNum(%1), MolIdx(%2) ] "
                "and [ AtomNum(%3), MolIdx(%4) ]")
                .arg(dihedrals[i][0]).arg(atoms[num_to_idx[dihedrals[i][0]]].getMolIndex())
                .arg(dihedrals[i][2]).arg(atoms[num_to_idx[dihedrals[i][3]]].getMolIndex()), CODELOC);
        }
    }

    /*****************************************************/
    /* Work out which impropers belong to this molecule. */
    /*****************************************************/

    mol_impropers.resize(num_mols);

    // Loop over all of the impropers.
    for (int i=0; i<nImpropers(); ++i)
    {
        // Get the molecule index.
        int molidx = atoms[num_to_idx[impropers[i][0]]].getMolIndex();

        // Make sure the terminal atom is also in the molecule.
        if (atoms[num_to_idx[impropers[i][3]]].getMolIndex() == molidx)
        {
            mol_impropers[molidx].append(i);
        }
        else
        {
            throw SireError::program_bug(QObject::tr("The atoms involved in an "
                "improper interaction are not in the same molecule: [ AtomNum(%1), MolIdx(%2) ] "
                "and [ AtomNum(%3), MolIdx(%4) ]")
                .arg(impropers[i][0]).arg(atoms[num_to_idx[impropers[i][0]]].getMolIndex())
                .arg(impropers[i][2]).arg(atoms[num_to_idx[impropers[i][3]]].getMolIndex()), CODELOC);
        }
    }

    /*******************************************************/
    /* Work out which cross-terms belong to this molecule. */
    /*******************************************************/

    mol_cross_terms.resize(num_mols);

    // Loop over all of the cross-terms.
    for (int i=0; i<nCrossTerms(); ++i)
    {
        // Get the molecule index.
        int molidx = atoms[num_to_idx[cross_terms[i][0]]].getMolIndex();

        // Make sure the terminal atom is also in the molecule.
        if (atoms[num_to_idx[cross_terms[i][3]]].getMolIndex() == molidx)
        {
            mol_cross_terms[molidx].append(i);
        }
        else
        {
            throw SireError::program_bug(QObject::tr("The atoms involved in a "
                "cross term are not in the same molecule: [ AtomNum(%1), MolIdx(%2) ] "
                "and [ AtomNum(%3), MolIdx(%4) ]")
                .arg(cross_terms[i][0]).arg(atoms[num_to_idx[cross_terms[i][0]]].getMolIndex())
                .arg(cross_terms[i][2]).arg(atoms[num_to_idx[cross_terms[i][3]]].getMolIndex()), CODELOC);
        }
    }
}

/** Helper function to recursively walk through bonded atoms in a molecule. */
void CharmmPSF::findBondedAtoms(int atom_num, int mol_idx, const QHash<int, int> &bonded_atoms,
    QHash<int, int> &atom_to_mol, QSet<qint64> &atoms_in_mol) const
{
    for (auto bonded_atom : bonded_atoms.values(atom_num))
    {
        // The molecule doesn't already contain this atom.
        if (not atoms_in_mol.contains(bonded_atom))
        {
            // Add the atom tho the molecule.
            atom_to_mol[bonded_atom] = mol_idx;
            atoms_in_mol.insert(bonded_atom);

            // Continue search from the next atom in the chain.
            findBondedAtoms(bonded_atom, mol_idx, bonded_atoms, atom_to_mol, atoms_in_mol);
        }
    }
}

/** Internal function used to grab the property, catching errors and signalling if
    the correct property has been found */
template<class T>
T CharmmPSF::getProperty(const PropertyName &prop, const MoleculeData &moldata, bool *found)
{
    if (moldata.hasProperty(prop))
    {
        const Property &p = moldata.property(prop);

        if (p.isA<T>())
        {
            *found = true;
            return p.asA<T>();
        }
    }

    *found = false;
    return T();
}

/** Construct PSF bond records from the set of two-atom functions. */
void CharmmPSF::getBondsFrom(const TwoAtomFunctions &funcs, const Molecule &sire_mol,
    QVector<QVector<qint64> > &local_bonds, QSet<QString> &bond_params, const PropertyMap &map)
{
    // Get the molecule info object.
    const auto molinfo = sire_mol.info();

    // Get the set of all bond functions.
    const auto potentials = funcs.potentials();

    // Store the number of bonds.
    const int num_bonds = potentials.count();

    // Resize the bonds vector.
    local_bonds.resize(num_bonds);

    // A multi-map between the start atom in a bond and the atom pair.
    // This allows us to reconstruct the correct bond order.
    QMultiMap<int, QVector<qint64> > bond_map;

    // Create the bond map.
    for (int i=0; i<num_bonds; ++i)
    {
        // Resize the bond vector.
        local_bonds[i].resize(2);

        // Get the potential.
        const auto potential = potentials.constData()[i];

        // Extract the cgAtomIdx for the two atoms.
        const auto idx0 = molinfo.cgAtomIdx(potential.atom0());
        const auto idx1 = molinfo.cgAtomIdx(potential.atom1());

        // Store the atom numbers.
        int atom0 = molinfo.number(idx0).value();
        int atom1 = molinfo.number(idx1).value();

        // If the atoms have an "atomtype" property then we can generate
        // CHARMM bond parameters.
        if ((sire_mol.atom(idx0).hasProperty(map["atomtype"])) and
            (sire_mol.atom(idx1).hasProperty(map["atomtype"])))
        {
            // Extract the atom types.
            QString type0 = sire_mol.atom(idx0).property<QString>(map["atomtype"]);
            QString type1 = sire_mol.atom(idx1).property<QString>(map["atomtype"]);

            // Create a string from the types.
            QString bond_atoms = QString("%1 %2").arg(type0, -4).arg(type1, -4);

            // Create the CHARMM bond parameters.
            const auto R = InternalPotential::symbols().bond().r();
            bond_params.insert(toHarmonicParameter(bond_atoms, potential.function(), R));
        }

        // Create the map value.
        QVector<qint64> value({atom0, atom1});

        // Insert the bond into the map.
        // Use negative atom number as the index so we can loop over in ascending order.
        bond_map.insert(-atom0, value);
    }

    // Now populate the bond vector.
    int i = 0;
    while (not bond_map.isEmpty())
    {
        local_bonds[i] = bond_map.take(bond_map.lastKey());
        i++;
    }
}

/** Construct PSF angle records from the set of three-atom functions. */
void CharmmPSF::getAnglesFrom(const ThreeAtomFunctions &funcs, const TwoAtomFunctions &ub_funcs,
    const Molecule &sire_mol, QVector<QVector<qint64> > &local_angles,
    QSet<QString> &angle_params, const PropertyMap &map)
{
    // Get the molecule info object.
    const auto molinfo = sire_mol.info();

    // Get the set of all angle functions.
    const auto potentials = funcs.potentials();

    // Get the set of all Urey-Bradley functions.
    const auto ub_potentials = ub_funcs.potentials();

    // Store the number of angles.
    const int num_angles = potentials.count();

    // Store the number of Urey-Bradley functions.
    const int num_ubs = ub_potentials.count();

    // Resize the angles vector.
    local_angles.resize(num_angles);

    // A multi-map between the start atom in an angle and the atom triplet.
    // This allows us to reconstruct the correct order for the PSF records.
    QMultiMap<int, QVector<qint64> > angle_map;

    // Populate the angle map.
    for (int i=0; i<num_angles; ++i)
    {
        // Resize the angles vector.
        local_angles[i].resize(3);

        // Get the potential.
        const auto potential = potentials.constData()[i];

        // Extract the cgAtomIdx for the three atoms.
        const auto idx0 = molinfo.cgAtomIdx(potential.atom0());
        const auto idx1 = molinfo.cgAtomIdx(potential.atom1());
        const auto idx2 = molinfo.cgAtomIdx(potential.atom2());

        // Store the atom numbers.
        int atom0 = molinfo.number(idx0).value();
        int atom1 = molinfo.number(idx1).value();
        int atom2 = molinfo.number(idx2).value();

        // If the atoms have an "atomtype" property then we can generate
        // CHARMM angle parameters.
        if ((sire_mol.atom(idx0).hasProperty(map["atomtype"])) and
            (sire_mol.atom(idx1).hasProperty(map["atomtype"])) and
            (sire_mol.atom(idx2).hasProperty(map["atomtype"])))
        {
            // Extract the atom types.
            QString type0 = sire_mol.atom(idx0).property<QString>(map["atomtype"]);
            QString type1 = sire_mol.atom(idx1).property<QString>(map["atomtype"]);
            QString type2 = sire_mol.atom(idx2).property<QString>(map["atomtype"]);

            // Create a string from the types.
            QString angle_atoms = QString("%1 %2 %3")
                .arg(type0, -4).arg(type1, -4).arg(type2, -4);

            // Create the CHARMM angle parameters.
            const auto Theta = InternalPotential::symbols().angle().theta();
            QString angle_param = toHarmonicParameter(angle_atoms, potential.function(), Theta, 3);

            // Now check whether there is an additional Urey-Bradley term for this angle.
            for (int j=0; j<num_ubs; ++j)
            {
                // Get the potential.
                const auto potential = ub_potentials.constData()[j];

                // Extract the cgAtomIdx for the three atoms.
                const auto ub_idx0 = molinfo.cgAtomIdx(potential.atom0());
                const auto ub_idx1 = molinfo.cgAtomIdx(potential.atom1());

                // If the atoms have an "atomtype" property then we can generate
                // CHARMM Urey-Bradley parameters.
                if ((sire_mol.atom(ub_idx0).hasProperty(map["atomtype"])) and
                    (sire_mol.atom(ub_idx1).hasProperty(map["atomtype"])))
                {
                    // Extract the atom types.
                    QString ub_type0 = sire_mol.atom(ub_idx0).property<QString>(map["atomtype"]);
                    QString ub_type1 = sire_mol.atom(ub_idx1).property<QString>(map["atomtype"]);

                    // If the index and type of the first and third atoms match, then add
                    // Urey-Bradley terms to the angle parameter string.
                    if ((ub_idx0  == idx0)  and (ub_idx1  == idx2) and
                        (ub_type0 == type0) and (ub_type1 == type2))
                    {
                        const auto R = InternalPotential::symbols().bond().r();
                        angle_param += toHarmonicParameter(QString(""), potential.function(), R);

                        break;
                    }
                }
            }

            // Insert the angle parameter string into the set.
            angle_params.insert(angle_param);
        }

        // Create the map value.
        QVector<qint64> value({atom0, atom1, atom2});

        // Insert the angle into the map.
        // Use negative atom number as the index so we can loop over in ascending order.
        angle_map.insert(-atom0, value);
    }

    // Now populate the angle vector.
    int i = 0;
    while (not angle_map.isEmpty())
    {
        local_angles[i] = angle_map.take(angle_map.lastKey());
        i++;
    }
}

/** Construct PSF four-atom records from the set of four-atom functions. */
void CharmmPSF::getFourAtomFrom(const FourAtomFunctions &funcs, const Molecule &sire_mol,
    QVector<QVector<qint64> > &four_atom, QSet<QString> &four_atom_params, const PropertyMap &map,
    bool is_improper)
{
    // Get the molecule info object.
    const auto molinfo = sire_mol.info();

    // Get the set of all four-atom functions.
    const auto potentials = funcs.potentials();

    // Store the number of functions.
    const int num_funcs = potentials.count();

    // Resize the fuctions vector.
    four_atom.resize(num_funcs);

    // A multi-map between the start atom in a four-atom function and the atom quartet.
    // This allows us to reconstruct the correct order for the PSF records.
    QMultiMap<int, QVector<qint64> > four_atom_map;

    // Populate the function map.
    for (int i=0; i<num_funcs; ++i)
    {
        // Resize the vector.
        four_atom[i].resize(4);

        // Get the potential.
        const auto potential = potentials.constData()[i];

        // Extract the cgAtomIdx for the four atoms.
        const auto idx0 = molinfo.cgAtomIdx(potential.atom0());
        const auto idx1 = molinfo.cgAtomIdx(potential.atom1());
        const auto idx2 = molinfo.cgAtomIdx(potential.atom2());
        const auto idx3 = molinfo.cgAtomIdx(potential.atom3());

        // Store the atom numbers.
        int atom0 = molinfo.number(idx0).value();
        int atom1 = molinfo.number(idx1).value();
        int atom2 = molinfo.number(idx2).value();
        int atom3 = molinfo.number(idx3).value();

        // If the atoms have an "atomtype" property then we can generate
        // CHARMM four-atom parameters.
        if ((sire_mol.atom(idx0).hasProperty(map["atomtype"])) and
            (sire_mol.atom(idx1).hasProperty(map["atomtype"])) and
            (sire_mol.atom(idx2).hasProperty(map["atomtype"])) and
            (sire_mol.atom(idx3).hasProperty(map["atomtype"])))
        {
            // Extract the atom types.
            QString type0 = sire_mol.atom(idx0).property<QString>(map["atomtype"]);
            QString type1 = sire_mol.atom(idx1).property<QString>(map["atomtype"]);
            QString type2 = sire_mol.atom(idx2).property<QString>(map["atomtype"]);
            QString type3 = sire_mol.atom(idx3).property<QString>(map["atomtype"]);

            // Create a string from the types.
            QString func_atoms = QString("%1 %2 %3 %4")
                .arg(type0, -4).arg(type1, -4).arg(type2, -4).arg(type3, -4);

            // Create the CHARMM parameters.
            if (is_improper)
            {
                const auto phi = InternalPotential::symbols().dihedral().phi();

                // Generate the improper parameter string.
                four_atom_params.insert(toHarmonicParameter(func_atoms, potential.function(), phi, 4));
            }
            else
            {
                // Generate the vector of dihedral parameter strings.
                QVector<QString> params = toFourAtomParameter(func_atoms, potential.function());

                for (const auto &param : params)
                    four_atom_params.insert(param);
            }
        }

        // Create the map value.
        QVector<qint64> value({atom0, atom1, atom2, atom3});

        // Insert into the map.
        // Use negative atom number as the index so we can loop over in ascending order.
        four_atom_map.insert(-atom0, value);
    }

    // Now populate the vector.
    int i = 0;
    while (not four_atom_map.isEmpty())
    {
        four_atom[i] = four_atom_map.take(four_atom_map.lastKey());
        i++;
    }
}

/** Convert a Sire two-atom function to a CHARMM bond parameter string. */
QString CharmmPSF::toHarmonicParameter(const QString &bond_atoms, const Expression &func,
    const Symbol &R, int num_atoms)
{
    // The function should be of the form "k(r - r0)^2".
    // We need to get the factors of R.
    const auto factors = func.expand(R);

    bool has_k = false;

    QStringList errors;

    double k = 0.0;
    double kr0_2 = 0.0;
    double kr0 = 0.0;

    for (const auto factor : factors)
    {
        if (factor.symbol() == R)
        {
            if (not factor.power().isConstant())
            {
                errors.append(QObject::tr("Power of R must be constant, not %1")
                                .arg(factor.power().toString()));
                continue;
            }

            if (not factor.factor().isConstant())
            {
                errors.append(QObject::tr("The value of K in K (R - R0)^2 must be constant. "
                                "Here it is %1").arg(factor.factor().toString()));
                continue;
            }

            double power = factor.power().evaluate(Values());

            if (power == 0.0)
            {
                // This is the constant.
                kr0_2 += factor.factor().evaluate(Values());
            }
            else if (power == 1.0)
            {
                // This is the -kR0 term.
                kr0 = factor.factor().evaluate(Values());
            }
            else if (power == 2.0)
            {
                // This is the R^2 term.
                if (has_k)
                {
                    // We cannot have two R2 factors?
                    errors.append(QObject::tr("Cannot have two R^2 factors!"));
                    continue;
                }

                k = factor.factor().evaluate(Values());
                has_k = true;
            }
            else
            {
                errors.append(QObject::tr("Power of R^2 must equal 2.0, 1.0 or 0.0, not %1")
                                .arg(power));
                continue;
            }
        }
        else
        {
            errors.append(QObject::tr("Cannot have a factor that does not include R. %1")
                        .arg(factor.symbol().toString()));
        }
    }

    if (kr0_2 < 0)
    {
        errors.append(QObject::tr("How can K R0^2 be negative? %1").arg(kr0_2));
    }

    double r0 = std::sqrt( kr0_2 / k );

    // kr0 should be equal to -2 k r0
    if (std::abs(k*r0 + 0.5 * kr0) > 0.001)
    {
        errors.append(QObject::tr("How can the power of R be %1. It should be 2 x %2 x %3 = %4.")
                        .arg(kr0).arg(k).arg(r0).arg(2*k*r0));
    }

    if (not errors.isEmpty())
    {
        throw SireError::incompatible_error(QObject::tr(
                "Cannot construct a CHARMM bond parameter K ( %1 - R0 )^2 from the "
                "expression %2, because\n%3")
                    .arg(R.toString()).arg(func.toString()).arg(errors.join("\n")), CODELOC );
    }

    // Improper function.
    if (num_atoms == 4)
    {
        // Convert out of plane angle to degrees.
        r0 = qRadiansToDegrees(r0);

        return QString("%1 %2 %3 %4")
            .arg(bond_atoms).arg(k, 10, 'f', 3).arg(0, 10).arg(r0, 10, 'f', 4);
    }

    else
    {
        // Angle function.
        // Convert equilibrium angle deviation to degrees.
        if (num_atoms == 3)
            r0 = qRadiansToDegrees(r0);

        return QString("%1 %2 %3").arg(bond_atoms).arg(k, 10, 'f', 3).arg(r0, 10, 'f', 4);
    }
}

/** Convert a Sire four-atom function to a set of CHARMM dihedral parameter strings. */
QVector<QString> CharmmPSF::toFourAtomParameter(const QString &dihedral_atoms, const Expression &func)
{
    const auto phi = InternalPotential::symbols().dihedral().phi();

    // This expression should be a sum of cos terms, plus constant terms.
    QList<Expression> cos_terms;
    double constant = 0.0;

    QStringList errors;

    if (func.base().isA<Sum>())
    {
        for (const auto child : func.base().asA<Sum>().children())
        {
            if (child.isConstant())
            {
                constant += func.factor() * child.evaluate(Values());
            }
            else if (child.base().isA<Cos>())
            {
                if (func.factor() == 1)
                {
                    cos_terms.append(child);
                }
                else
                {
                    cos_terms.append(func.factor() * child);
                }
            }
            else
            {
                errors.append(QObject::tr("Cannot interpret the dihedral expression "
                   "from '%1' as it should be a series of cosine terms involving %2.")
                        .arg(func.toString()).arg(phi.toString()));
            }
        }
    }
    else
    {
        if (func.isConstant())
        {
            constant += func.evaluate(Values());
        }
        else if (func.base().isA<Cos>())
        {
            cos_terms.append(func);
        }
        else
        {
            errors.append( QObject::tr("Cannot interpret the dihedral expression "
               "from '%1' as it should be a series of cosine terms involving %2.")
                    .arg(func.toString()).arg(phi.toString()));
        }
    }

    // Next extract all of the data from the cos terms.
    QList<std::tuple<double,double,double> > terms;

    double check_constant = 0;

    for (const auto cos_term : cos_terms)
    {
        // Term should be of the form 'k cos( periodicity * phi - phase )'.
        double k = cos_term.factor();
        check_constant += k;

        double periodicity = 0.0;
        double phase = 0.0;

        const auto factors = cos_term.base().asA<Cos>().argument().expand(phi);

        bool ok = true;

        for (const auto factor : factors)
        {
            if (not factor.power().isConstant())
            {
                errors.append(QObject::tr("Power of phi must be constant, not %1")
                                .arg(factor.power().toString()));
                ok = false;
                continue;
            }

            if (not factor.factor().isConstant())
            {
                errors.append(QObject::tr("The value of periodicity must be "
                                "constant. Here it is %1").arg(factor.factor().toString()));
                ok = false;
                continue;
            }

            double power = factor.power().evaluate(Values());

            if (power == 0.0)
            {
                // This is the constant phase.
                phase += factor.factor().evaluate(Values());
            }
            else if (power == 1.0)
            {
                // This is the periodicity * phi term
                periodicity = factor.factor().evaluate(Values());
            }
            else
            {
                errors.append(QObject::tr("Power of phi must equal 1.0 or 0.0, not %1")
                                    .arg(power));
                ok = false;
                continue;
            }
        }

        if (ok)
        {
            terms.append(std::make_tuple(k, periodicity, phase));
        }
    }

    // Now sum up the individual values of 'k'. These should equal the constant, which
    // is the aggregate of the "k [ 1 + cos ]" terms.
    if ( std::abs(constant - check_constant) > 0.001 )
    {
        errors.append(QObject::tr("The set of constants should sum up to the same total, "
            "i.e. should equal %1. Instead they equal %2. This is weird.")
                .arg(constant).arg(check_constant));
    }

    if (not errors.isEmpty())
    {
        throw SireError::incompatible_error(QObject::tr(
            "Cannot extract a CHARMM dihedral parameter from '%1' as "
            "the expression must be a series of terms of type "
            "'k{ 1 + cos[ per %2 - phase ] }'. Errors include\n%3")
                .arg(func.toString()).arg(phi.toString()).arg(errors.join("\n")),
                    CODELOC);
    }

    // The vector of dihedral parameter strings.
    QVector<QString> parameter_strings;

    if (not terms.isEmpty())
    {
        // Append each parameter string to the vector, converting the phase
        // shift to degrees.
        for (const auto term : terms)
        {
            parameter_strings.append(QString("%1 %2 %3 %4")
                .arg(dihedral_atoms)
                .arg(std::get<0>(term), 10, 'f', 4)
                .arg(std::get<1>(term), 3)
                .arg(-qRadiansToDegrees(std::get<2>(term)), 10, 'f', 2));
        }
    }

    return parameter_strings;
}
