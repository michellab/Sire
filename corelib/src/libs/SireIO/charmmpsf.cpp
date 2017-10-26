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
#include "SireMol/atomcoords.h"
#include "SireMol/atomelements.h"
#include "SireMol/molecule.h"
#include "SireMol/moleditor.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireUnits/units.h"

#include <QtMath>

using namespace SireIO;
using namespace SireMol;
using namespace SireBase;
using namespace SireSystem;
using namespace SireStream;

const RegisterParser<CharmmPSF> register_psf;
static const RegisterMetaType<CharmmPSF> r_psf;
static const RegisterMetaType<PSFAtom> r_psfatom(NO_ROOT);

QDataStream SIREIO_EXPORT &operator<<(QDataStream &ds, const PSFAtom &psfatom)
{
    writeHeader(ds, r_psfatom, 1);

    SharedDataStream sds(ds);

    sds << psfatom.number << psfatom.segment << psfatom.res_num << psfatom.res_name
        << psfatom.name << psfatom.type << psfatom.charge << psfatom.mass;

    return ds;
}

QDataStream SIREIO_EXPORT &operator>>(QDataStream &ds, PSFAtom &psfatom)
{
    VersionID v = readHeader(ds, r_psfatom);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> psfatom.number >> psfatom.segment >> psfatom.res_num >> psfatom.res_name
            >> psfatom.name >> psfatom.type >> psfatom.charge >> psfatom.mass;
    }
    else
        throw version_error(v, "1", r_psfatom, CODELOC);

    return ds;
}

QDataStream SIREIO_EXPORT &operator<<(QDataStream &ds, const CharmmPSF &psf)
{
    writeHeader(ds, r_psf, 1);

    ds << psf.atoms << psf.bonds << psf.angles << psf.dihedrals << psf.impropers
       << psf.cross_terms << psf.coords << static_cast<const MoleculeParser&>(psf);

    return ds;
}

QDataStream SIREIO_EXPORT &operator>>(QDataStream &ds, CharmmPSF &psf)
{
    VersionID v = readHeader(ds, r_psf);

    if (v == 1)
    {
        ds >> psf.atoms >> psf.bonds << psf.angles >> psf.dihedrals << psf.impropers
           >> psf.cross_terms >> psf.coords >> static_cast<MoleculeParser&>(psf);
    }
    else
        throw version_error(v, "1", r_psf, CODELOC);

    return ds;
}

/** Default constructor. */
PSFAtom::PSFAtom() :
    number(0),
    res_num(0),
    type("X"),
    charge(0),
    mass(0)
{

}

/** Constructor. */
PSFAtom::PSFAtom(const QString &line, QStringList &errors) :
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

    line.append(QString("%1").arg(number, 7, 10));
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

    // Store the name of the input file.
    this->filename = filename;

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

    // Store the name of the input file.
    this->filename = filename;

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
        /*return QObject::tr("CharmmPSF( nMolecules() = %1, "
            "nResidues() = %2, nAtoms() = %3 )")
            .arg(nMolecules()).arg(nSubstructures()).arg(nAtoms());*/
    }
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
    auto parse_line = [](const QString& line, int iline,
        int num_lines, int num_records, int width, QStringList &errors)
    {
        // Work out the number of records on this line.

        // The default number.
        int num_line = width;

        // We're on the last record line.
        if (iline == (num_lines -1))
            num_line = width - (width*num_lines - num_records);

        // Tokenize the line, splitting using a single whitespace character.
        QStringList data = line.simplified().split(QRegExp("\\s"));

        // Check that the line has the right number of records.
        if (data.count() != num_line)
        {
            errors.append(QObject::tr("The number of records on the line is "
                "incorrect. Expected %1, found %2: %3")
                .arg(width).arg(data.count()).arg(line));

            return;
        }

        // Work out the index of the first record.
        int istart = 4 * iline;

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
                        atoms[i] = PSFAtom(lines()[iline+i], local_errors);
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
                    atoms[i] = PSFAtom(lines()[++iline], parse_warnings);
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
            int num_record_lines = qCeil(num_bonds/4.0);
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
            int num_record_lines = qCeil(num_angles/3.0);
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
            int num_record_lines = qCeil(num_dihedrals/2.0);
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
            int num_record_lines = qCeil(num_impropers/2.0);
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
            int num_record_lines = qCeil(num_cross/2.0);
        }
    }

    this->setScore(num_atoms);
}

/** Use the data contained in this parser to create a new System of molecules,
    assigning properties based on the mapping in 'map' */
System CharmmPSF::startSystem(const PropertyMap &map) const
{
    /*const int nmols = nMolecules();

    if (nmols == 0)
        return System();

    QVector<Molecule> mols(nmols);
    Molecule *mols_array = mols.data();

    if (usesParallel())
    {
        tbb::parallel_for( tbb::blocked_range<int>(0, nmols),
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
    system.setProperty(map["filename"].source(), StringProperty(filename));
    system.setProperty(map["fileformat"].source(), StringProperty(this->formatName()));*/

    //return system;
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
    MolStructureEditor mol;

    return mol;
}

/** Internal function used to get the molecule structure for molecule 'imol'. */
MolEditor CharmmPSF::getMolecule(int imol, const PropertyMap &map) const
{
    /*return mol.setProperty(map["coordinates"], coords)
              .setProperty(map["charge"], charges)
              .setProperty(map["element"], elements)
              .commit();*/
}
