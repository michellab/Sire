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

#include "SireError/errors.h"
#include "SireIO/errors.h"

#include "SireMol/atomcharges.h"
#include "SireMol/atomelements.h"
#include "SireMol/atommasses.h"
#include "SireMol/connectivity.h"
#include "SireMol/molecule.h"
#include "SireMol/moleditor.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireUnits/units.h"

#include <QDateTime>
#include <QtMath>

using namespace SireBase;
using namespace SireIO;
using namespace SireMol;
using namespace SireStream;
using namespace SireSystem;

const RegisterParser<CharmmPSF> register_psf;
static const RegisterMetaType<CharmmPSF> r_psf;
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

QDataStream SIREIO_EXPORT &operator<<(QDataStream &ds, const CharmmPSF &psf)
{
    writeHeader(ds, r_psf, 1);

    ds << psf.atoms << psf.bonds << psf.angles << psf.dihedrals << psf.impropers
       << psf.cross_terms << psf.molecules << psf.num_to_idx
       << static_cast<const MoleculeParser&>(psf);

    return ds;
}

QDataStream SIREIO_EXPORT &operator>>(QDataStream &ds, CharmmPSF &psf)
{
    VersionID v = readHeader(ds, r_psf);

    if (v == 1)
    {
        ds >> psf.atoms >> psf.bonds >> psf.angles >> psf.dihedrals >> psf.impropers
           >> psf.cross_terms >> psf.molecules >> psf.num_to_idx
           >> static_cast<MoleculeParser&>(psf);
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
PSFAtom::PSFAtom(const SireMol::Atom &atom, bool is_ter, QStringList &errors) :
    index(0),
    mol_idx(0),
    number(0),
    res_num(0),
    type("X"),
    charge(0),
    mass(0)
{
}

/** Generate a PSF record from the atom data. */
QString PSFAtom::toPSFRecord() const
{
    QString line;

    line.append(QString("%1").arg(number, 8));
    line.append(QString(" %1").arg(segment, -4));
    line.append(QString(" %1").arg(res_num, -4));
    line.append(QString(" %1").arg(res_name, -4));
    line.append(QString(" %1").arg(toPDBName(), 4));
    line.append(QString(" %1").arg(type, -4));
    line.append(QString(" %1").arg(charge, 10, 'f', 6));
    line.append(QString(" %1").arg(mass, 13, 'f', 4));
    line.append(QString("%1").arg("0", 12));

    return line;
}

/** Convert the name to PDB format. */
QString PSFAtom::toPDBName() const
{
    QString pdb_name = name;

    /* PDB atom names are a maximum of 4 characters wide and obey the
       following rules:

       All names start in the second position (i.e. start with a space)
       unless the first character of the name is a digit.

       3 character names:
         - Name starts with a digit
             --> append a space, e.g. "1HE "
         - Else...
             --> prepend a space, e.g. " NE1"

       2 character names:
         - Start in second position, e.g. " CA "

       1 character names:
         - Put in second position, e.g. " H  "
     */

    // Truncate the name to 4 characters.
    if (pdb_name.count() > 4) pdb_name = pdb_name.left(4);

    // Apply formatting rules.
    else if (pdb_name.count() < 4)
    {
        if (pdb_name.count() == 3)
        {
            if (pdb_name[0].isDigit()) pdb_name.append(" ");
            else                       pdb_name.prepend(" ");
        }
        else if (pdb_name.count() == 2)
        {
            pdb_name.prepend(" ");
            pdb_name.append(" ");
        }
        else if (pdb_name.count() == 1)
        {
            pdb_name.prepend(" ");
            pdb_name.append("  ");
        }
    }

    return pdb_name;
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

/** Constructor */
CharmmPSF::CharmmPSF() : ConcreteProperty<CharmmPSF,MoleculeParser>()
{}

/** Construct to read in the data from the file called 'filename'. The
passed property map can be used to pass extra parameters to control
the parsing */
CharmmPSF::CharmmPSF(const QString &filename, const PropertyMap &map)
     : ConcreteProperty<CharmmPSF,MoleculeParser>(filename,map)
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
CharmmPSF::CharmmPSF(const QStringList &lines, const PropertyMap &map)
     : ConcreteProperty<CharmmPSF,MoleculeParser>(lines,map)
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
CharmmPSF::CharmmPSF(const SireSystem::System &system, const PropertyMap &map)
     : ConcreteProperty<CharmmPSF,MoleculeParser>(map)
{
    //this->operator=(parsed);
}

/** Copy constructor */
CharmmPSF::CharmmPSF(const CharmmPSF &other) :
    ConcreteProperty<CharmmPSF,MoleculeParser>(other),
    atoms(other.atoms),
    bonds(other.bonds),
    angles(other.angles),
    dihedrals(other.dihedrals),
    impropers(other.impropers),
    cross_terms(other.cross_terms),
    num_to_idx(other.num_to_idx),
    molecules(other.molecules),
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
        angles = other.angles;
        dihedrals = other.dihedrals;
        impropers = other.impropers;
        cross_terms = other.cross_terms;
        num_to_idx = other.num_to_idx;
        molecules = other.molecules;
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
    // TODO: Add BioSimSpace version info.
    lines.append("PSF");
    lines.append("");
    lines.append("       2 !NTITLE");
    lines.append(QString(" REMARKS DATE:%1    created by BioSimSpace (v)")
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
            tbb::parallel_for( tbb::blocked_range<int>(0, num_records),
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
        record_lines[0] = QString("%1 !NBONDS: bonds").arg(num_records, 8);

        if (usesParallel())
        {
            tbb::parallel_for( tbb::blocked_range<int>(0, num_lines),
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
            tbb::parallel_for( tbb::blocked_range<int>(0, num_lines),
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
            tbb::parallel_for( tbb::blocked_range<int>(0, num_lines),
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
            tbb::parallel_for( tbb::blocked_range<int>(0, num_lines),
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
            tbb::parallel_for( tbb::blocked_range<int>(0, num_lines),
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
    return "CHARMMPSF";
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

/** Return the number of angle records. */
int CharmmPSF::nAngles() const
{
    return angles.count();
}

/** Return the number of dihedral records. */
int CharmmPSF::nDihedrals() const
{
    return dihedrals.count();
}

/** Return the number of improper records. */
int CharmmPSF::nImpropers() const
{
    return impropers.count();
}

/** Return the number of cross-term records. */
int CharmmPSF::nCrossTerms() const
{
    return cross_terms.count();
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
        // Tokenize the line, splitting using a single whitespace character.
        QStringList data = lines()[iline].simplified().split(QRegExp("\\s"));

        // Atom records.
        if (data.last() == "!NATOM")
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

                tbb::parallel_for( tbb::blocked_range<int>(0, num_atoms),
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
        else if (data.contains("!NBOND:"))
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

                tbb::parallel_for( tbb::blocked_range<int>(0, num_lines),
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
        else if (data.contains("!NTHETA:"))
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

                tbb::parallel_for( tbb::blocked_range<int>(0, num_lines),
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
        else if (data.contains("!NPHI:"))
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

                tbb::parallel_for( tbb::blocked_range<int>(0, num_lines),
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
        else if (data.contains("!NIMPHI:"))
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

                tbb::parallel_for( tbb::blocked_range<int>(0, num_lines),
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
        else if (data.contains("!NCRTERM:"))
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

                tbb::parallel_for( tbb::blocked_range<int>(0, num_lines),
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

    return system;
}

/** Use the data contained in this parser to add information from the file to
    the molecules that exist already in the passed System. For example, this
    may be used to add coordinate data from this file to the molecules in
    the passed System that are missing coordinate data. */
void CharmmPSF::addToSystem(System &system, const PropertyMap &map) const
{
    //you should loop through each molecule in the system and work out
    //which ones are described in the file, and then add data from the file
    //to thise molecules.
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

    // Get the info object that can map between AtomNum to AtomIdx etc.
    const auto molinfo = mol.info();

    // Instantiate the atom property objects that we need.
    AtomCharges charges(molinfo);
    AtomMasses  masses(molinfo);

    // Now loop through the atoms in the molecule and set each property.
    for (int i=0; i<nAtoms(imol); ++i)
    {
        // Store a reference to the current atom.
        const auto &atom = atoms[molecules[imol][i]];

        // Determine the CGAtomIdx for this atom.
        auto cgatomidx = molinfo.cgAtomIdx(AtomNum(atom.getNumber()));

        // Set the properties.
        masses.set(cgatomidx, atom.getMass() * SireUnits::g_per_mol);
        charges.set(cgatomidx, double(atom.getCharge()) * SireUnits::mod_electron);
    }

    // Now work out which bonds are part of this molecule.

    // The indices of the bonds in the molecule.
    QSet<int> mol_bonds;

    // Loop over all of the bonds.
    for (int i=0; i<nBonds(); ++i)
    {
        // The bond is part of this molecule.
        if (atoms[num_to_idx[bonds[i][0]]].getMolIndex() == imol)
        {
            // Make sure the terminal atom is also in the molecule.
            if (atoms[num_to_idx[bonds[i][1]]].getMolIndex() == imol)
            {
                mol_bonds.insert(i);
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
    }

    // Connectivity object for bonded atoms.
    ConnectivityEditor connectivity(molinfo);

    // Loop over all bonds in the molecule.
    for (const auto &bond : mol_bonds)
    {
        // Get the atom indices for the two bonds.
        int idx1 = num_to_idx[bonds[bond][0]];
        int idx2 = num_to_idx[bonds[bond][1]];

        // Add the bond to the connectivity object.
        connectivity.connect(AtomNum(atoms[idx1].getNumber()),
                             AtomNum(atoms[idx2].getNumber()));
    }

    return mol.setProperty(map["charge"], charges)
              .setProperty(map["mass"], masses)
              .setProperty(map["connectivity"], connectivity)
              .commit();
}

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
