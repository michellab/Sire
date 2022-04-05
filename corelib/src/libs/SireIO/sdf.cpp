/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2022  Christopher Woods
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


#include "SireIO/sdf.h"

#include "SireSystem/system.h"

#include "SireBase/parallel.h"
#include "SireBase/stringproperty.h"

#include "SireError/errors.h"
#include "SireIO/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireMol/atomcharges.h"
#include "SireMol/atomcoords.h"
#include "SireMol/atomelements.h"
#include "SireMol/atommasses.h"
#include "SireMol/errors.h"
#include "SireMol/molecule.h"
#include "SireMol/moleditor.h"
#include "SireMol/connectivity.h"
#include "SireMol/bondid.h"

#include "SireBase/propertylist.h"

#include "SireUnits/units.h"

#include "sire_version.h"

#include <QFile>
#include <QtMath>

using namespace SireBase;
using namespace SireIO;
using namespace SireMol;
using namespace SireStream;
using namespace SireSystem;
using namespace SireUnits;

namespace SireIO { namespace detail {

class SDFAtom
{
public:
    SDFAtom() : x(0), y(0), z(0), mass_difference(0), chg_difference(0)
    {}

    ~SDFAtom()
    {}

    QString toString() const
    {
        return QString("%1  %2 %3 %4  %5 %6  %7")
                .arg(name).arg(x).arg(y).arg(z)
                .arg(mass_difference).arg(chg_difference)
                .arg(fields.join(":"));
    }

    void completeFields()
    {
        while (fields.count() < 10)
        {
            fields.append("0");
        }
    }

    QString name;
    double x;
    double y;
    double z;
    qint32 mass_difference;
    qint32 chg_difference;
    QStringList fields;
};

class SDFBond
{
public:
    SDFBond() : atom0(0), atom1(0), typ(0), stereoscopy(0)
    {}

    ~SDFBond()
    {}

    QString toString() const
    {
        return QString("%1-%2  %3  %4  %5")
                    .arg(atom0).arg(atom1)
                    .arg(typ).arg(stereoscopy)
                    .arg(fields.join(":"));
    }

    void completeFields()
    {
        while (fields.count() < 3)
        {
            fields.append("0");
        }
    }

    qint32 atom0;
    qint32 atom1;
    qint32 typ;
    qint32 stereoscopy;
    QStringList fields;
};

class SDFMolecule
{
public:
    SDFMolecule()
    {}

    ~SDFMolecule()
    {}

    bool isValid() const
    {
        return atoms.count() > 0;
    }

    QString toString() const
    {
        QStringList lines;

        lines.append(QString("name = %1").arg(name));
        lines.append(QString("software = %1").arg(software));
        lines.append(QString("comment = %1").arg(comment));

        lines.append(QString("counts: %1").arg(counts.join(" ")));

        lines.append(QString("nAtoms == %1").arg(atoms.count()));

        for (auto atom : atoms)
        {
            lines.append(atom.toString());
        }

        lines.append(QString("nBonds == %1").arg(bonds.count()));

        for (auto bond : bonds)
        {
            lines.append(bond.toString());
        }

        for (auto key : properties.keys())
        {
            lines.append(QString("property %1").arg(key));

            for (auto s : properties[key])
            {
                lines.append(s);
            }
        }

        for (auto key : data.keys())
        {
            lines.append(QString("data %1").arg(key));

            for (auto s : data[key])
            {
                lines.append(s);
            }
        }

        return lines.join("\n");
    }

    int getCharge(int i) const
    {
        return 0;
    }

    Element getElement(int i) const
    {
        return Element(atoms[i].name);
    }

    double getMass(int i) const
    {
        return 0.0;
    }

    void completeCounts()
    {
        while (counts.count() < 10)
        {
            counts.append("0");
        }
    }

    QString name;
    QString software;
    QString comment;
    QStringList counts;

    QVector<SDFAtom> atoms;
    QVector<SDFBond> bonds;

    QHash<QString, QStringList> properties;
    QHash<QString, QStringList> data;
};

}} // end of namespace SireIO::detail

QDataStream &operator<<(QDataStream &ds, const SireIO::detail::SDFAtom &atom)
{
    ds << atom.name << atom.x << atom.y << atom.z
       << atom.chg_difference << atom.mass_difference
       << atom.fields;

    return ds;
}

QDataStream &operator>>(QDataStream &ds, SireIO::detail::SDFAtom &atom)
{
    ds >> atom.name >> atom.x >> atom.y >> atom.z
       >> atom.chg_difference >> atom.mass_difference
       >> atom.fields;

    return ds;
}

QDataStream &operator<<(QDataStream &ds, const SireIO::detail::SDFBond &bond)
{
    ds << bond.atom0 << bond.atom1 << bond.typ
       << bond.stereoscopy << bond.fields;

    return ds;
}

QDataStream &operator>>(QDataStream &ds, SireIO::detail::SDFBond &bond)
{
    ds >> bond.atom0 >> bond.atom1 >> bond.typ
       >> bond.stereoscopy >> bond.fields;

    return ds;
}

QDataStream &operator<<(QDataStream &ds, const SireIO::detail::SDFMolecule &mol)
{
    ds << mol.atoms << mol.bonds << mol.properties << mol.data;
    return ds;
}

QDataStream &operator>>(QDataStream &ds, SireIO::detail::SDFMolecule &mol)
{
    ds >> mol.atoms >> mol.bonds >> mol.properties >> mol.data;
    return ds;
}

using namespace SireIO::detail;

const RegisterParser<SDF> register_sdf;
static const RegisterMetaType<SDF> r_sdf;

QDataStream &operator<<(QDataStream &ds, const SDF &sdf)
{
    writeHeader(ds, r_sdf, 1);

    SharedDataStream sds(ds);

    sds << sdf.molecules << sdf.parse_warnings
        << static_cast<const MoleculeParser&>(sdf);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, SDF &sdf)
{
    VersionID v = readHeader(ds, r_sdf);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> sdf.molecules >> sdf.parse_warnings
            >> static_cast<MoleculeParser&>(sdf);
    }
    else
        throw version_error(v, "1", r_sdf, CODELOC);

    return ds;
}

/** Constructor */
SDF::SDF() : ConcreteProperty<SDF,MoleculeParser>()
{}

/** Construct to read in the data from the file called 'filename'. The
    passed property map can be used to pass extra parameters to control
    the parsing */
SDF::SDF(const QString &filename, const PropertyMap &map) :
    ConcreteProperty<SDF,MoleculeParser>(filename,map)
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
SDF::SDF(const QStringList &lines, const PropertyMap &map) :
    ConcreteProperty<SDF,MoleculeParser>(lines,map)
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

QStringList toLines(const SDFMolecule &molecule)
{
    QStringList lines;

    lines.append(molecule.name);
    lines.append(QString("  -Sire-%1").arg(SIRE_VERSION));
    lines.append(molecule.comment);

    if (molecule.counts.size() != 10)
    {
        throw SireError::program_bug(
            QObject::tr("Problem with the counts line! %1")
                    .arg(molecule.counts.join(",")), CODELOC);
    }

    QString count_line = QString("%1%2")
                            .arg(molecule.atoms.count(), 3)
                            .arg(molecule.bonds.count(), 3);

    for (int i=0; i<9; ++i)
    {
        QString count = molecule.counts[i].trimmed();
        count.truncate(3);
        count_line += QString("%1").arg(count, 3);
    }

    QString count = molecule.counts[9].trimmed();

    count.truncate(6);
    count_line += QString("%1").arg(count, 6);

    lines.append(count_line);

    for (const auto &atom : molecule.atoms)
    {
        QString name = atom.name.trimmed();
        name.truncate(3);

        QString atom_line = QString("%1%2%3 %4%5%6")
                                .arg(atom.x, 10, 'f', 4)
                                .arg(atom.y, 10, 'f', 4)
                                .arg(atom.z, 10, 'f', 4)
                                .arg(name, -3)
                                .arg(atom.mass_difference, 2)
                                .arg(atom.chg_difference, 3);

        if (atom.fields.count() != 10)
        {
            throw SireError::program_bug(
                QObject::tr("Problem with the atom line! %1")
                        .arg(atom.fields.join(",")), CODELOC);
        }

        for (int i=0; i<10; ++i)
        {
            QString f = atom.fields[i].trimmed();
            f.truncate(3);
            atom_line += QString("%1").arg(f, 3);
        }

        lines.append(atom_line);
    }

    for (const auto &bond : molecule.bonds)
    {
        QString bond_line = QString("%1%2%3%4")
                                .arg(bond.atom0, 3)
                                .arg(bond.atom1, 3)
                                .arg(bond.typ, 3)
                                .arg(bond.stereoscopy, 3);

        if (bond.fields.count() != 3)
        {
            throw SireError::program_bug(
                QObject::tr("Problem with the bond line! %1")
                        .arg(bond.fields.join(",")), CODELOC);
        }

        for (int i=0; i<3; ++i)
        {
            QString f = bond.fields[i].trimmed();
            f.truncate(3);

            bond_line += QString("%1").arg(f, 3);
        }

        lines.append(bond_line);
    }

    for (const auto &key : molecule.properties.keys())
    {
        for (const auto &value : molecule.properties[key])
        {
            QString k = key.trimmed();
            k.truncate(3);

            lines.append(QString("M  %1%2").arg(k).arg(value));
        }
    }

    lines.append("M  END");

    bool has_last_line = false;

    for (const auto &key : molecule.data.keys())
    {
        lines.append(QString("> <%1>").arg(key.trimmed()));

        for (const auto &value : molecule.data[key])
        {
            lines.append(value);
        }

        lines.append("");
        has_last_line = true;
    }

    if (not has_last_line)
        lines.append("");

    lines.append("$$$$");

    return lines;
}

QStringList toLines(const QList<SDFMolecule> &molecules,
                    bool uses_parallel=false)
{
    QStringList lines;

    if (uses_parallel and molecules.count() > 1)
    {
        QVector<QStringList> molecule_lines(molecules.count());
        QStringList *mollines_ptr = molecule_lines.data();

        tbb::parallel_for(tbb::blocked_range<int>(0, molecules.count()),
                          [&](const tbb::blocked_range<int> r)
        {
            for (int i=r.begin(); i<r.end(); ++i)
            {
                mollines_ptr[i] = ::toLines(molecules[i]);
            }
        });

        for (const auto &l : molecule_lines)
        {
            lines += l;
        }
    }
    else
    {
        for (const auto &molecule : molecules)
        {
            lines += ::toLines(molecule);
        }
    }

    return lines;
}

SDFMolecule parseMolecule(const Molecule &molecule,
                          QStringList &errors,
                          const PropertyMap &map)
{
    // Store the number of atoms in the molecule.
    int num_atoms = molecule.nAtoms();

    // Early exit.
    if (num_atoms == 0) return SDFMolecule();

    // TODO: Do we want a hard limit on the number of atoms?
    if (num_atoms > 999)
    {
        errors.append(QObject::tr("The number of atoms (%1) exceeds "
            " the SDF file format limit (999)!").arg(num_atoms));

        return SDFMolecule();
    }

    Connectivity connectivity;

    if (molecule.hasProperty(map["connectivity"]))
    {
        connectivity = molecule.property(
                        map["connectivity"]).asA<Connectivity>();
    }

    int num_bonds = connectivity.nConnections();

    if (num_bonds > 999)
    {
        errors.append(QObject::tr("The number of bonds (%1) exceeds "
            " the SDF file format limit (999)!").arg(num_bonds));
    }

    if (not molecule.hasProperty(map["coordinates"]))
    {
        errors.append(QObject::tr("The molecule is missing a coordinates "
              "property. This is needed for the SDF file!"));
        return SDFMolecule();
    }

    SDFMolecule sdfmol;

    if (molecule.hasProperty(map["name"]))
    {
        sdfmol.name = molecule.property(map["name"]).toString();
    }

    if (molecule.hasProperty(map["software"]))
    {
        sdfmol.software = molecule.property(map["software"]).toString();
    }

    if (molecule.hasProperty(map["comment"]))
    {
        sdfmol.comment = molecule.property(map["comment"]).toString();
    }

    for (int i=0; i<num_atoms; ++i)
    {
        const auto atom = molecule.atom(AtomIdx(i));

        const auto coords = atom.property<Vector>(map["coordinates"]);

        SDFAtom sdf_atom;

        sdf_atom.x = coords.x();
        sdf_atom.y = coords.y();
        sdf_atom.z = coords.z();

        if (atom.hasProperty(map["element"]))
        {
            sdf_atom.name = atom.property<Element>(map["element"]).symbol();
        }
        else
        {
            sdf_atom.name = atom.name();
        }

        sdf_atom.completeFields();

        sdfmol.atoms.append(sdf_atom);
    }

    for (auto bond : connectivity.getBonds())
    {
        qDebug() << bond.toString();
    }

    sdfmol.completeCounts();

    return sdfmol;
}

/** Construct this parser by extracting all necessary information from the
    passed SireSystem::System, looking for the properties that are specified
    in the passed property map */
SDF::SDF(const SireSystem::System &system, const PropertyMap &map) :
    ConcreteProperty<SDF,MoleculeParser>(map)
{
    // Get the MolNums of each molecule in the System - this returns the
    // numbers in MolIdx order.
    const QVector<MolNum> molnums = system.getMoleculeNumbers().toVector();

    // Store the number of molecules.
    const int nmols = molnums.count();

    // No molecules in the system.
    if (nmols == 0)
    {
        this->operator=(SDF());
        return;
    }

    QVector<SDFMolecule> mols(nmols);

    SDFMolecule *mols_ptr = mols.data();

    if (usesParallel())
    {
        QMutex mutex;

        tbb::parallel_for(tbb::blocked_range<int>(0, nmols),
                          [&](const tbb::blocked_range<int> r)
        {
            // Create local data objects.
            QStringList local_errors;

            for (int i=r.begin(); i<r.end(); ++i)
            {
                // Now parse the rest of the molecular data, i.e. atoms, residues, etc.
                mols_ptr[i] = parseMolecule(system[molnums[i]].molecule(),
                                            local_errors, map);
            }

            if (not local_errors.isEmpty())
            {
                // Acquire a lock.
                QMutexLocker lkr(&mutex);

                // Update the warning messages.
                parse_warnings += local_errors;
            }
        });
    }
    else
    {
        for (int i=0; i<nmols; ++i)
        {
            // Now parse the rest of the molecular data, i.e. atoms, residues, etc.
            mols_ptr[i] = parseMolecule(system[molnums[i]].molecule(),
                                        parse_warnings, map);
        }
    }

    QList<SDFMolecule> checked_mols;

    for (int i=0; i<nmols; ++i)
    {
        if (mols_ptr[i].isValid())
        {
            checked_mols.append(mols_ptr[i]);
        }
    }

    mols.clear();

    if (checked_mols.count() > 0)
    {
        QStringList lines = ::toLines(checked_mols, usesParallel());
        checked_mols.clear();

        // Reparse the lines as a self-consistency check.
        SDF parsed(lines, map);

        QStringList copy_warnings = parse_warnings;
        this->operator=(parsed);
        parse_warnings = copy_warnings;
    }
}

/** Copy constructor */
SDF::SDF(const SDF &other) :
    ConcreteProperty<SDF,MoleculeParser>(other),
    molecules(other.molecules),
    parse_warnings(other.parse_warnings)
{}

/** Destructor */
SDF::~SDF()
{}

/** Copy assignment operator */
SDF& SDF::operator=(const SDF &other)
{
    if (this != &other)
    {
        this->molecules = other.molecules;
        this->parse_warnings = other.parse_warnings;

        MoleculeParser::operator=(other);
    }

    return *this;
}

/** Comparison operator */
bool SDF::operator==(const SDF &other) const
{
    return MoleculeParser::operator==(other);
}

/** Comparison operator */
bool SDF::operator!=(const SDF &other) const
{
    return not operator==(other);
}

/** Return the C++ name for this class */
const char* SDF::typeName()
{
    return QMetaType::typeName( qMetaTypeId<SDF>() );
}

/** Return the C++ name for this class */
const char* SDF::what() const
{
    return SDF::typeName();
}

/** Return the parser that has been constructed by reading in the passed
    file using the passed properties */
MoleculeParserPtr SDF::construct(const QString &filename,
                                 const PropertyMap &map) const
{
    return SDF(filename,map);
}

/** Return the parser that has been constructed by reading in the passed
    text lines using the passed properties */
MoleculeParserPtr SDF::construct(const QStringList &lines,
                                 const PropertyMap &map) const
{
    return SDF(lines,map);
}

/** Return the parser that has been constructed by extract all necessary
    data from the passed SireSystem::System using the specified properties */
MoleculeParserPtr SDF::construct(const SireSystem::System &system,
                                 const PropertyMap &map) const
{
    return SDF(system,map);
}

/** Return a string representation of this parser */
QString SDF::toString() const
{
    if (lines().isEmpty())
        return QObject::tr("SDF::null");
    else
    {
        return QObject::tr("SDF( nMolecules() == %1 )")
                        .arg(this->nMolecules());
    }
}

/** Return the format name that is used to identify this file format within Sire */
QString SDF::formatName() const
{
    return "SDF";
}

/** Return any warnings raised when parsing this file */
QStringList SDF::parseWarnings() const
{
    return this->parse_warnings;
}

/** Return the number of molecules loaded in this file */
int SDF::nMolecules() const
{
    return this->molecules.count();
}

/** Return the total number of atoms. */
int SDF::nAtoms() const
{
    int num_atoms = 0;

    for (const auto &mol : this->molecules)
        num_atoms += mol.atoms.count();

    return num_atoms;
}

/** Return the number of atoms in molecule 'i'. */
int SDF::nAtoms(int i) const
{
    i = Index(i).map(this->molecules.count());
    return this->molecules[i].atoms.count();
}

/** Return a description of the file format */
QString SDF::formatDescription() const
{
    return QObject::tr("Structure Data File (SDF) format files.");
}

/** Return the suffixes that these files are normally associated with */
QStringList SDF::formatSuffix() const
{
    static const QStringList suffixes = { "SDF" };
    return suffixes;
}

/** Return whether or not this is a lead parser. The lead parser is responsible
    for starting the process of turning the parsed file into the System. There
    must be one and one-only lead parser in a set of parsers creating a System */
bool SDF::isLead() const
{
    return true;
}

/** The SDF cannot follow another lead parsers. */
bool SDF::canFollow() const
{
    return false;
}

/** Function that is called to assert that this object is sane. This
    should raise an exception if the parser is in an invalid state */
void SDF::assertSane() const
{
    if (this->nMolecules() == 0)
    {
        QStringList errors = this->parse_warnings;

        if (not errors.isEmpty())
        {
            throw SireIO::parse_error(QObject::tr(
                "There were errors reading the SDF format "
                "file:\n%1").arg(errors.join("\n\n")), CODELOC);
        }
    }
}

/** Internal function that is used to parse a single molecule's set
    of lines from the SDF file */
void SDF::parseMoleculeLines(const PropertyMap &map,
                             const QStringList &l)
{
    /* File format is decribed here:
        https://www.herongyang.com/Molecule/SDF-Format-Specification.html
        http://www.nonlinear.com/progenesis/sdf-studio/v0.9/faq/sdf-file-format-guidance.aspx
     */

    // The first three lines are the header block
    if (l.count() < 4)
    {
        // this is not a valid SDF file
        this->parse_warnings.append(
            QObject::tr("There are no enough lines for this "
                        "to be a valid SDF-formatted file.")
        );
        return;
    }

    // first line is the molecule name
    const auto molname = l[0].simplified();
    // then the software used to generate the SDF file
    const auto software = l[1].simplified();
    // then a user-supplied comment (if any)
    const auto comment = l[2].simplified();

    // now the counts line - 12 fixed-width fields. The first eleven are 3
    // characters long. The last is 6 characters long
    const auto counts_line = l[3];

    if (counts_line.size() < 39)
    {
        this->parse_warnings.append(
            QObject::tr("The counts line in this SDF file does not "
                        "have enough characters! '%1'. It should be "
                        "at least 39 characters wide.").arg(counts_line)
        );
        return;
    }

    QStringList counts;

    for (int i=0; i<11; ++i)
    {
        counts.append(counts_line.mid(i*3, 3));
    }

    // the last counts line item can also be a string!
    counts.append(counts_line.mid(33, 6));

    bool ok;

    // Atom counter.
    const int natoms = counts.takeFirst().toInt(&ok);

    if (not ok)
    {
        this->parse_warnings.append(
            QObject::tr("Cannot interpret the number of atoms from the "
                        "counts line: %1").arg(counts_line));
    }

    // Bonds counter
    const int nbonds = counts.takeFirst().toInt(&ok);

    if (not ok)
    {
        this->parse_warnings.append(
            QObject::tr("Cannot interpret the number of bonds from the "
                        "counts line: %1").arg(counts_line));
    }

    if (natoms == 0)
    {
        // nothing to read?
        this->parse_warnings.append(
            QObject::tr("The number of atoms to read is set to zero?")
        );
        return;
    }

    // next, read in the atoms
    if (l.count() < 4 + natoms + nbonds)
    {
        this->parse_warnings.append(
            QObject::tr("There aren't enough lines in this file to "
                        "contain all of the atoms and bonds. File is "
                        "corrupted?")
        );
        return;
    }

    QVector<SDFAtom> atoms(natoms);

    for (int i=0; i<natoms; ++i)
    {
        QString line = l[i+4];

        //000000000011111111112222222222333333333344444444445555555555666666666
        //012345678901234567890123456789012345678901234567890123456789012345678
        //    0.5369    0.9749    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
        while (line.size() < 69)
        {
            line.append(" ");
        }

        SDFAtom &atom = atoms[i];

        auto assert_ok = [&](bool is_ok, int line_num,
                             const QString &line, const QString &field)
                        {
                            if (not is_ok)
                            {
                                this->parse_warnings.append(
                                    QObject::tr("Atom line %1 has a problem "
                                      "with field '%2'. '%3'")
                                        .arg(line_num)
                                        .arg(field)
                                        .arg(line)
                                );

                                return false;
                            }

                            return true;
                        };

        atom.x = line.midRef(0,10).toDouble(&ok);

        if (not assert_ok(ok, i+1, line, "x"))
            return;

        atom.y = line.midRef(10, 10).toDouble(&ok);

        if (not assert_ok(ok, i+1, line, "y"))
            return;

        atom.z = line.midRef(20, 10).toDouble(&ok);

        if (not assert_ok(ok, i+1, line, "z"))
            return;

        atom.name = line.mid(31,3);

        atom.mass_difference = line.midRef(34, 2).toInt(&ok);

        if (not assert_ok(ok, i+1, line, "mass difference"))
            return;

        if (atom.mass_difference < -3 or atom.mass_difference > 4)
        {
            this->parse_warnings.append(QObject::tr(
                "Only mass differences between -3 and 4 are supported. "
                "Cannot have a difference of %1 on line %2. '%3'")
                    .arg(atom.mass_difference)
                    .arg(i+1).arg(line));
            return;
        }

        atom.chg_difference = line.midRef(36, 3).toInt(&ok);

        if (not assert_ok(ok, i+1, line, "charge difference"))
            return;

        if (atom.chg_difference < 0 or atom.chg_difference > 7)
        {
            this->parse_warnings.append(QObject::tr(
                "Only charge differences between 0 and 7 are supported. "
                "Cannot have a difference of %1 on line %2. '%3'")
                    .arg(atom.chg_difference)
                    .arg(i+1).arg(line));
            return;
        }

        // ten more fields of 3 characters each. We won't convert these
        // to numbers - just leave as strings
        for (int j=0; j<10; ++j)
        {
            atom.fields.append(line.mid(38+(3*j), 3));
        }
    }

    QVector<SDFBond> bonds(nbonds);

    for (int i=0; i<nbonds; ++i)
    {
        QString line = l[4 + natoms + i];

        //000000000011111111112
        //012345678901234567890
        //  1  2  1  0  0  0  0
        while (line.size() < 20)
        {
            line.append(" ");
        }

        SDFBond &bond = bonds[i];

        auto assert_ok = [&](bool is_ok, int line_num,
                             const QString &line, const QString &field)
                        {
                            if (not is_ok)
                            {
                                this->parse_warnings.append(
                                    QObject::tr("Bond line %1 has a problem "
                                      "with field '%2'. '%3'")
                                        .arg(line_num)
                                        .arg(field)
                                        .arg(line)
                                );

                                return false;
                            }

                            return true;
                        };

        bond.atom0 = line.midRef(0, 3).toInt(&ok);

        if (not assert_ok(ok, i+1, line, "atom0"))
            return;

        bond.atom1 = line.midRef(3, 3).toInt(&ok);

        if (not assert_ok(ok, i+1, line, "atom1"))
            return;

        bond.typ = line.midRef(6, 3).toInt(&ok);

        if (not assert_ok(ok, i+1, line, "bond type"))
            return;

        bond.stereoscopy = line.midRef(9, 3).toInt(&ok);

        if (not assert_ok(ok, i+1, line, "stereoscopy"))
            return;

        // now add on the three fields, which we will leave as strings
        for (int j=0; j<3; ++j)
        {
            bond.fields.append(line.mid(12+(3*j), 3));
        }
    }

    // now read in the properties...
    int linenum = 4 + natoms + nbonds;

    QHash<QString, QStringList> properties;

    while (linenum < l.count())
    {
        const auto line = l[linenum];

        if (line.startsWith("M  END"))
        {
            // end of the properties
            break;
        }
        else if (line.startsWith("M  "))
        {
            // this is a new property line
            if (line.size() >= 6)
            {
                QString key = line.mid(3, 3);
                properties[key].append(line.mid(6));
            }
        }

        linenum += 1;
    }

    // now read in the data...
    QString key = "";
    QHash<QString, QStringList> data;

    while (linenum < l.count())
    {
        QString line = l[linenum];

        if (line.startsWith("$$$$"))
        {
            // end of the data (and molecule)
            break;
        }
        else if (line.startsWith("> "))
        {
            line = line.mid(2);

            int start_idx = line.indexOf("<");
            int end_idx = line.indexOf(">");

            if (start_idx >= 0 and end_idx >= 0)
            {
                key = line.mid(start_idx+1, end_idx-start_idx-1);
            }
            else
            {
                key = "";
            }
        }
        else if (key.size() > 0)
        {
            if (line.simplified().size() > 0)
            {
                data[key].append(line);
            }
        }

        linenum += 1;
    }

    SDFMolecule molecule;
    molecule.name = molname;
    molecule.software = software;
    molecule.comment = comment;
    molecule.counts = counts;
    molecule.atoms = atoms;
    molecule.bonds = bonds;
    molecule.properties = properties;
    molecule.data = data;

    this->molecules.append(molecule);
}

/** Internal function that is used to actually parse the data contained
    in the lines of the file */
void SDF::parseLines(const PropertyMap &map)
{
    const auto &l = this->lines();

    // how many molecules are there? Each molecule is separated
    // by a '$$$$', so we count these lines
    int nmolecules = 0;

    for (const auto &line : l)
    {
        if (line == "$$$$")
        {
            nmolecules += 1;
        }
    }

    if (nmolecules == 1)
    {
        this->parseMoleculeLines(map, l.toList());
    }
    else
    {
        // break this into sets of lines and parse independently
        QStringList lines;

        for (const auto &line : l)
        {
            if (line == "$$$$")
            {
                lines.append(line);
                this->parseMoleculeLines(map, lines);
                lines = QStringList();
            }
            else
            {
                lines.append(line);
            }
        }
    }

    this->setScore(100 * this->molecules.count());
}

/** Use the data contained in this parser to create a new System of molecules,
    assigning properties based on the mapping in 'map'. */
System SDF::startSystem(const PropertyMap &map) const
{
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

/** Internal function used to get the molecule structure for molecule 'imol'. */
MolStructureEditor SDF::getMolStructure(const SDFMolecule &sdfmol,
                                        const PropertyName &cutting,
                                        const QString &resname) const
{
    // Make sure that there are atoms in the molecule.
    if (sdfmol.atoms.count() == 0)
    {
        throw SireError::program_bug(QObject::tr(
            "Strange - there are no atoms in molecule %1?")
                .arg(sdfmol.toString()), CODELOC);
    }

    // First step is to build the structure of the molecule, i.e.
    // the layout of cutgroups, residues and atoms.
    MolStructureEditor mol;

    // SDF files don't define resiues, so place all atoms in a single residue
    auto cutgroup = mol.add(CGName("1"));

    auto res = mol.add(ResNum(1));
    res.rename(ResName(resname));

    // Add each atom in the residue to the molecule.
    for (int i=0; i<sdfmol.atoms.count(); ++i)
    {
        const auto &sdfatom = sdfmol.atoms.at(i);

        auto atom = cutgroup.add(AtomNum(i+1));
        atom.rename(AtomName(sdfatom.name.trimmed()));
        // Reparent the atom to its residue.
        atom.reparent(ResIdx(0));
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
MolEditor SDF::getMolecule(int imol, const PropertyMap &map) const
{
    // Make sure the molecule (model) index is within range.
    imol = Index(imol).map(this->molecules.count());

    const auto sdfmol = this->molecules.at(imol);

    // Make sure that there are atoms in the molecule.
    if (sdfmol.atoms.count() == 0)
    {
        throw SireError::program_bug(QObject::tr(
            "Strange - there are no atoms in molecule %1?")
                .arg(imol), CODELOC);
    }

    QString resname = "MOL";

    if (map.specified("resname"))
    {
        resname = map["resname"].value().toString();
    }

    // First, construct the layout of the molecule (sorting of atoms into residues and cutgroups).
    auto mol = this->getMolStructure(sdfmol, map["cutting"],
                                     resname).commit().edit();

    // Get the info object that can map between AtomNum to AtomIdx etc.
    const auto molinfo = mol.info();

    // Atom property objects.
    AtomCoords         coords(molinfo);
    AtomCharges        charges(molinfo);
    AtomElements       elements(molinfo);
    AtomMasses         masses(molinfo);

    // Now loop through the atoms in the molecule and set each property.
    for (int i=0; i<sdfmol.atoms.count(); ++i)
    {
        // Store a reference to the current atom.
        const auto &atom = sdfmol.atoms.at(i);

        // Determine the CGAtomIdx for this atom.
        auto cgatomidx = molinfo.cgAtomIdx(AtomNum(i+1));

        // Set the properties.
        coords.set(cgatomidx, Vector(atom.x, atom.y, atom.z));
        charges.set(cgatomidx, int(sdfmol.getCharge(i)) * SireUnits::mod_electron);
        elements.set(cgatomidx, sdfmol.getElement(i));
        masses.set(cgatomidx, sdfmol.getMass(i) * SireUnits::g_per_mol);
    }

    // need to do the bonds...

    return mol.setProperty(map["coordinates"], coords)
              .setProperty(map["formal_charge"], charges)
              .setProperty(map["element"], elements)
              .setProperty(map["mass"], masses)
              .setProperty(map["software"], SireBase::wrap(sdfmol.software))
              .setProperty(map["name"], SireBase::wrap(sdfmol.name))
              .setProperty(map["comment"], SireBase::wrap(sdfmol.comment))
              .commit();
}
