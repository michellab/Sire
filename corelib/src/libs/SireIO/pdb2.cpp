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

#include "SireIO/pdb2.h"

#include "SireMM/pdbparams.h"

#include "SireSystem/system.h"

#include "SireBase/parallel.h"
#include "SireBase/stringproperty.h"

#include "SireError/errors.h"
#include "SireIO/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireMol/atomcharges.h"
#include "SireMol/atomelements.h"
#include "SireMol/atomcoords.h"
#include "SireMol/molecule.h"
#include "SireMol/moleditor.h"

#include "SireUnits/units.h"

using namespace SireIO;
using namespace SireMM;
using namespace SireMol;
using namespace SireBase;
using namespace SireSystem;
using namespace SireStream;
using namespace SireUnits;

const RegisterParser<PDB2> register_pdb;
static const RegisterMetaType<PDB2> r_pdb2;
static const RegisterMetaType<PDBAtom> r_pdbatom(NO_ROOT);

QDataStream SIREIO_EXPORT &operator<<(QDataStream &ds, const PDBAtom &pdbatom)
{
    writeHeader(ds, r_pdbatom, 1);

    SharedDataStream sds(ds);

    sds << pdbatom.record << pdbatom.serial << pdbatom.name << pdbatom.alt_loc
        << pdbatom.res_name << pdbatom.chain_id << pdbatom.res_num << pdbatom.insert_code
        << pdbatom.coord << pdbatom.occupancy << pdbatom.temperature << pdbatom.element
        << pdbatom.charge << pdbatom.is_het << pdbatom.is_ter;

    return ds;
}

QDataStream SIREIO_EXPORT &operator>>(QDataStream &ds, PDBAtom &pdbatom)
{
    VersionID v = readHeader(ds, r_pdbatom);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> pdbatom.record >> pdbatom.serial >> pdbatom.name >> pdbatom.alt_loc
            >> pdbatom.res_name >> pdbatom.chain_id >> pdbatom.res_num >> pdbatom.insert_code
            >> pdbatom.coord >> pdbatom.occupancy >> pdbatom.temperature >> pdbatom.element
            >> pdbatom.charge >> pdbatom.is_het >> pdbatom.is_ter;
    }
    else
        throw version_error(v, "1", r_pdbatom, CODELOC);

    return ds;
}

QDataStream SIREIO_EXPORT &operator<<(QDataStream &ds, const PDB2 &pdb2)
{
    writeHeader(ds, r_pdb2, 1);

    SharedDataStream sds(ds);

    sds << pdb2.atoms << pdb2.residues << pdb2.chains << pdb2.num_ters
        << pdb2.parse_warnings << static_cast<const MoleculeParser&>(pdb2);

    return ds;
}

QDataStream SIREIO_EXPORT &operator>>(QDataStream &ds, PDB2 &pdb2)
{
    VersionID v = readHeader(ds, r_pdb2);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> pdb2.atoms >> pdb2.residues >> pdb2.chains >> pdb2.num_ters
            >> pdb2.parse_warnings >> static_cast<MoleculeParser&>(pdb2);
    }
    else
        throw version_error(v, "1", r_pdb2, CODELOC);

    return ds;
}

/** Default constructor. */
PDBAtom::PDBAtom() :
    occupancy(1.0),
    element("X"),
    charge(0),
    is_het(false),
    is_ter(false)
{
}

/** Constructor.
    @param line
        An ATOM record line from a PDB file.

    @param errors
        An array of error messages.
 */
PDBAtom::PDBAtom(const QString &line, QStringList &errors) :
    record(line),
    occupancy(1.0),
    element("X"),
    charge(0),
    is_het(false),
    is_ter(false)
{
    if (line.length() < 54)
    {
        errors.append(QObject::tr("Cannot parse ATOM record "
            " since it does not match the format! '%1'").arg(line));

        return;
    }

    // Flag that this is a HETATM record.
    if (line.leftRef(6) == "HETATM")
        is_het = true;

    // Extract the atom serial number.
    bool ok;
    int tmp_int = line.midRef(6,5).toInt(&ok);

    if (not ok)
    {
        errors.append(QObject::tr("Cannot extract the atom serial number "
            "from part (%1) from line '%2'").arg(line.mid(6,5)).arg(line));

        return;
    }
    serial = tmp_int;

    // Extract the atom name.
    // Somewhere we'll test this against a list of valid names.
    name = line.mid(12,4).simplified();

    // Extract the alternative atom location indicator.
    alt_loc = line[16];

    // Extract the residue name.
    res_name = line.mid(17,3).simplified();

    // Extract the chain ID.
    chain_id = line[21];

    // Extract the residue sequence number.
    tmp_int = line.midRef(22,4).toInt(&ok);

    if (not ok)
    {
        errors.append(QObject::tr("Cannot extract the residue sequence number "
            "from part (%1) from line '%2'").arg(line.mid(22,4)).arg(line));

        return;
    }
    res_num = tmp_int;

    // Extract the residue insertion code.
    insert_code = line[26];

    // Now try to extract the coordinate data.
    bool ok_x, ok_y, ok_z;
    double x = line.midRef(30,8).toDouble(&ok_x);
    double y = line.midRef(38,8).toDouble(&ok_y);
    double z = line.midRef(46,8).toDouble(&ok_z);

    if (not (ok_x and ok_y and ok_z))
    {
        errors.append(QObject::tr("There was a problem reading the coordinate "
            "values of x, y, and z from the data '%1' in line '%2'")
            .arg(line.mid(30,24)).arg(line));

        return;
    }

    // Store the atom coordinates.
    coord = Vector(x, y, z);

    // Now try to extract the "optional" data from the atom record.
    // We store this data if present, and throw warnings if it does not
    // match the format from the PDB standard (don't bail out).

    // Extract occupancy data.
    double tmp_dbl = line.midRef(54,6).toDouble(&ok);

    if (not ok)
    {
        errors.append(QObject::tr("There was a problem reading the occupancy "
            "value from the data '%1' in line '%2'")
            .arg(line.mid(54,6)).arg(line));
    }
    else
    {
        occupancy = tmp_dbl;

        // Check occupancy is valid.
        if ((occupancy < 1) or
            (occupancy > 1))
        {
            errors.append(QObject::tr("The occupancy (%1) was out of range! "
                "Setting to default value of 1.0.")
                .arg(occupancy));

            occupancy = 1;
            return;
        }
    }

    // Extract temperature data.
    tmp_dbl = line.midRef(60,6).toDouble(&ok);

    if (not ok)
    {
        errors.append(QObject::tr("There was a problem reading the temperature "
            "value from the data '%1' in line '%2'")
            .arg(line.mid(60,6)).arg(line));
    }
    else
    {
        temperature = tmp_dbl;

        // Check temperature is valid.
        if (temperature < 0)
        {
            errors.append(QObject::tr("The temperature factor (%1) is negative! "
                "Setting to the default value of 0.0.")
                .arg(temperature));

            temperature = 0;
            return;
        }
    }

    // Extract the element name.
    element = line.mid(76,2);

    // Extract the charge on the atom.
    QString chargeString = line.mid(78,2);

    // Format should be 2+ or 1- or something like that.
    int factor = 1;

    if (chargeString.contains("-"))
        factor = -1;

    chargeString.remove("-").remove("+");
    charge = factor * chargeString.toInt(&ok);

    // Something went wrong - ignore the charge.
    if (not ok)
        charge = 0;
}

/** Constructor.
    @param atom
        A reference to a SireMol::Atom object.

    @param is_ter
        Whether this is a terminal atom.

    @param errors
        An array of error messages.
 */
PDBAtom::PDBAtom(const SireMol::Atom &atom, bool is_ter, QStringList &errors) :
    serial(atom.number().value()),
    name(atom.name().value()),
    occupancy(1.0),
    temperature(0.0),
    element("X"),
    charge(0),
    is_het(false),
    is_ter(is_ter)
{
    // The atom must have atomic coordinates to be valid.
    if (not atom.hasProperty("coordinates"))
    {
        errors.append(QObject::tr("The atom does not have coordinates!"));

        return;
    }

    // Extract the atomic coordinates.
    coord = atom.property<SireMaths::Vector>("coordinates");

    // The atom is within a residue.
    if (atom.isWithinResidue())
    {
        res_name = atom.residue().name().value();
        res_num  = atom.residue().number().value();
    }

    // The atom is within a chain.
    if (atom.isWithinChain())
    {
        // TODO: Make sure we can handle situations where the
        // chain ID is a string.
        chain_id = atom.chain().name().value().at(0);
    }

    // Extract the occupancy.
    if (atom.hasProperty("occupancy"))
    {
        occupancy = atom.property<double>("occupancy");
    }

    // Extract the temperature factor.
    if (atom.hasProperty("beta-factor"))
    {
        temperature = atom.property<double>("beta-factor");
    }

    // Extract the element name.
    if (atom.hasProperty("element"))
    {
        element = atom.property<Element>("element").symbol();
    }
    else
    {
        // TODO: Infer the element...
    }

    // Extract the atomic charge.
    if (atom.hasProperty("formal-charge"))
    {
        // TODO: This doesn't seem to be working in all cases.
        charge = atom.property<SireUnits::Dimension::Charge>("formal-charge").value();
        charge /= SireUnits::mod_electron;
    }

    // Determine whether this is a HETATM.
    if (atom.hasProperty("is-het"))
    {
        if (atom.property<QString>("is-het") == "True")
            is_het = true;
    }
}

/** Convert the name to PDB format. */
QString PDBAtom::toPDBName() const
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

/** Set the terminal atom flag.
    @param is_ter
        Whether this is a terminal atom.
 */
void PDBAtom::setTerminal(bool _is_ter)
{
    is_ter = _is_ter;
}

/** Whether this is a HETATM. */
bool PDBAtom::isHet() const
{
    return is_het;
}

/** Whether this is a terminal atom. */
bool PDBAtom::isTer() const
{
    return is_ter;
}

/** Get the atom serial number. */
qint64 PDBAtom::getSerial() const
{
    return serial;
}

/** Set the atom serial number. */
void PDBAtom::setSerial(int serial)
{
    this->serial = serial;
}

/** Get the atom name. */
QString PDBAtom::getName() const
{
    return name;
}

/** Get the residue name. */
QString PDBAtom::getResName() const
{
    return res_name;
}

/** Get the chain id. */
QChar PDBAtom::getChainId() const
{
    return chain_id;
}

/** Get the residue sequence number. */
qint64 PDBAtom::getResNum() const
{
    return res_num;
}

/** Get the residue sequence number. */
SireMaths::Vector PDBAtom::getCoord() const
{
    return coord;
}

/** Get the occupancy. */
double PDBAtom::getOccupancy() const
{
    return occupancy;
}

/** Get the temperature factor. */
double PDBAtom::getTemperature() const
{
    return temperature;
}

/** Set the temperature factor. */
void PDBAtom::setTemperature(double temperature)
{
    this->temperature = temperature;
}

/** Get the element symbol. */
QString PDBAtom::getElement() const
{
    return element;
}

/** Get the charge on the atom. */
qint64 PDBAtom::getCharge() const
{
    return charge;
}

/** Return the C++ name for this class */
const char* PDBAtom::typeName()
{
    return QMetaType::typeName( qMetaTypeId<PDBAtom>() );
}

/** Return a PDB record line for this atom. */
QString PDBAtom::toPDBRecord() const
{
    QString line;

    if (is_het) line.append("HETATM");
    else        line.append("ATOM  ");

    // Append the atom serial number. Since we enforce a hard
    // limit of 99999 atoms, we can assume that the serial is
    // 5 characters or less.
    // TODO: Handle non-standard serial numbers.
    line.append(QString("%1").arg(serial, 5, 10));

    // Append the atom name, converting it to the correct PDB formatting.
    line.append(QString(" %1").arg(toPDBName(), 4));

    // Append the alternate location indicator.
    if (not alt_loc.isNull()) line.append(alt_loc);
    else                      line.append(" ");

    // Append the residue name, truncating if neccessary.
    line.append(res_name.left(3));

    // Append the chain ID.
    if (not chain_id.isNull()) line.append(QString(" %1").arg(chain_id));
    else                       line.append("  ");

    // Append the residue sequence number.
    line.append(QString("%1").arg(res_num, 4, 10));

    // Append the residue insertion code
    if (not insert_code.isNull()) line.append(insert_code);
    else                          line.append(" ");

    // Append the atomic coordinates.
    line.append(QString("   %1\%2\%3")
        .arg(coord[0], 8, 'f', 3)
        .arg(coord[1], 8, 'f', 3)
        .arg(coord[2], 8, 'f', 3));

    // Append the occupancy.
    line.append(QString("%1").arg(occupancy, 6, 'f', 2));

    // Append the beta factor.
    line.append(QString("%1").arg(temperature, 6, 'f', 2));

    // Append the element, truncating if neccessary.
    line.append(QString("          %1").arg(element.left(2)));

    // Apped the atomic charge, truncating as necessary.
    if (charge != 0)
    {
        line.append(QString("%1").arg(charge));
        if (charge < 0) line.append("-");
        else            line.append("+");
    }

    return line;
}

/** Return a string representation of this object */
QString PDBAtom::toString() const
{
    return QObject::tr("PDBAtom::null");
}

/** Constructor */
PDB2::PDB2() : ConcreteProperty<PDB2,MoleculeParser>()
{
}

/** Construct to read in the data from the file called 'filename'. The
    passed property map can be used to pass extra parameters to control
    the parsing */
PDB2::PDB2(const QString &filename, const PropertyMap &map) :
    ConcreteProperty<PDB2,MoleculeParser>(filename,map)
{
    //the file has been read into memory and is available via
    //the MoleculeParser::lines() function.

    //a parameter has also been read in MoleculeParser to say whether
    //we are allowed to use multiple cores to parse the file, e.g.
    //MoleculeParser::usesParallel() will be true

    // Set the file name.
    this->filename = filename;

    //parse the data in the parse function
    this->parseLines(map);

    //now make sure that everything is correct with this object
    this->assertSane();
}

/** Construct to read in the data from the passed text lines. The
    passed property map can be used to pass extra parameters to control
    the parsing */
PDB2::PDB2(const QStringList &lines, const PropertyMap &map) :
    ConcreteProperty<PDB2,MoleculeParser>(lines,map)
{
    //the file has been read into memory and is available via
    //the MoleculeParser::lines() function.

    //a parameter has also been read in MoleculeParser to say whether
    //we are allowed to use multiple cores to parse the file, e.g.
    //MoleculeParser::usesParallel() will be true

    // Set the file name.
    this->filename = filename;

    //parse the data in the parse function
    this->parseLines(map);

    //now make sure that everything is correct with this object
    this->assertSane();
}

/** Construct this parser by extracting all necessary information from the
    passed SireSystem::System, looking for the properties that are specified
    in the passed property map */
PDB2::PDB2(const SireSystem::System &system, const PropertyMap &map) :
    ConcreteProperty<PDB2,MoleculeParser>(map)
{
    // Get the MolNums of each molecule in the System - this returns the
    // numbers in MolIdx order.
    const QVector<MolNum> molnums = system.getMoleculeNumbers().toVector();

    // Store the number of molecules.
    const int nmols = molnums.count();

    // No molecules in the system.
    if (nmols == 0)
    {
        this->operator=(PDB2());
        return;
    }

    // The list of lines.
    QStringList lines;

    // Lines for different PDB data records (one for each molecule).
    QVector<QVector<QString> > atom_lines(nmols);

    // Set the name of the file from which the SireSystem was constructed.
    if (system.properties().hasProperty("filename"))
    {
        filename = system.property("filename").toString();
    }
    else filename.clear();

    if (usesParallel())
    {
        QMutex mutex;

        tbb::parallel_for( tbb::blocked_range<int>(0, nmols),
                           [&](const tbb::blocked_range<int> r)
        {
            // Create local data objects.
            QStringList local_errors;

            for (int i=r.begin(); i<r.end(); ++i)
            {
                // Now parse the rest of the molecular data, i.e. atoms, residues, etc.
                parseMolecule(system[molnums[i]].molecule(),
                    atom_lines[i], local_errors, map);
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
            parseMolecule(system[molnums[i]].molecule(),
                atom_lines[i], parse_warnings, map);
        }
    }

    // Now assemble the lines from the record data for each molecule.
    // We do this in serial since the order matters.
    for (int i=0; i<nmols; ++i)
    {
        // MODEL record.
        lines.append(QString("MODEL     %1").arg(i+1));

        // ATOM records.
        lines.append(atom_lines[i].toList());

        // ENDMDL record.
        lines.append("ENDMDL");
    }

    lines.append("END");

    // Reparse the lines as a self-consistency check.
    PDB2 parsed(lines, map);

    this->operator=(parsed);
}

/** Copy constructor */
PDB2::PDB2(const PDB2 &other) :
    ConcreteProperty<PDB2,MoleculeParser>(other),
    atoms(other.atoms),
    residues(other.residues),
    chains(other.chains),
    num_ters(other.num_ters),
    filename(other.filename),
    parse_warnings(other.parse_warnings)
{}

/** Destructor */
PDB2::~PDB2()
{}

/** Copy assignment operator */
PDB2& PDB2::operator=(const PDB2 &other)
{
    if (this != &other)
    {
        this->atoms = other.atoms;
        this->residues = other.residues;
        this->chains = other.chains;
        this->num_ters = other.num_ters;
        this->filename = other.filename;
        this->parse_warnings = other.parse_warnings;

        MoleculeParser::operator=(other);
    }

    return *this;
}

/** Comparison operator */
bool PDB2::operator==(const PDB2 &other) const
{
    return MoleculeParser::operator==(other);
}

/** Comparison operator */
bool PDB2::operator!=(const PDB2 &other) const
{
    return not operator==(other);
}

/** Return the C++ name for this class */
const char* PDB2::typeName()
{
    return QMetaType::typeName( qMetaTypeId<PDB2>() );
}

/** Return the C++ name for this class */
const char* PDB2::what() const
{
    return PDB2::typeName();
}

/** Return the parser that has been constructed by reading in the passed
    file using the passed properties */
MoleculeParserPtr PDB2::construct(const QString &filename,
                                  const PropertyMap &map) const
{
    return PDB2(filename,map);
}

/** Return the parser that has been constructed by reading in the passed
    text lines using the passed properties */
MoleculeParserPtr PDB2::construct(const QStringList &lines,
                                  const PropertyMap &map) const
{
    return PDB2(lines,map);
}

/** Return the parser that has been constructed by extract all necessary
    data from the passed SireSystem::System using the specified properties */
MoleculeParserPtr PDB2::construct(const SireSystem::System &system,
                                  const PropertyMap &map) const
{
    return PDB2(system,map);
}

/** Return a string representation of this parser */
QString PDB2::toString() const
{
    return QObject::tr("PDB2::null");
}

/** Return the format name that is used to identify this file format within Sire */
QString PDB2::formatName() const
{
    return "PDB";
}

/** Return a description of the file format */
QString PDB2::formatDescription() const
{
    return QObject::tr("Protein Data Bank (PDB) format files.");
}

/** Return the suffixes that these files are normally associated with */
QStringList PDB2::formatSuffix() const
{
    static const QStringList suffixes = { "PDB" };
    return suffixes;
}

/** Return the number of models (molecules). */
int PDB2::nModels() const
{
    return atoms.count();
}

/** Return the total number of atoms. */
int PDB2::nAtoms() const
{
    int num_atoms = 0;

    for (int i=0; i<atoms.count(); ++i)
        num_atoms += atoms[i].count();

    return num_atoms;
}

/** Return the number of atoms in model 'i'. */
int PDB2::nAtoms(int i) const
{
    return atoms[i].count();
}

/** Return the total number of TER records. */
int PDB2::nTers() const
{
    int n = 0;

    for (int i=0; i<atoms.count(); ++i)
        n += num_ters[i];

    return n;
}

/** Return the total number of TER records in model 'i'. */
int PDB2::nTers(int i) const
{
    return num_ters[i];
}

/** Return the atoms. */
QVector<QVector<PDBAtom> > PDB2::getAtoms() const
{
    return atoms;
}

/** Function that is called to assert that this object is sane. This
    should raise an exception if the parser is in an invalid state */
void PDB2::assertSane() const
{
    QStringList errors;

    if (not errors.isEmpty())
    {
        throw SireIO::parse_error(QObject::tr("There were errors reading the PDB format "
          "file:\n%1").arg(errors.join("\n\n")), CODELOC);
    }
}

/** Internal function that is used to actually parse the data contained
    in the lines of the file */
void PDB2::parseLines(const PropertyMap &map)
{
    /* File format is decribed here:
        http://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html

       Note that the PDB file format has been frozen since 21 November 2012.
     */

    // Atom counter.
    int nats = 0;

    // Model index.
    int imdl = 0;

    // Frame index;
    int iframe = 0;

    // Number of atoms in the previous model.
    // Each MODEL section should contain the same number of atoms.
    // For each MODEL record there should be a corresponding ENDMDL.
    int prev_nats = 0;

    // The indices for lines containing atom data.
    // For each frame we identify the lines containing atom data, then parse them.
    QVector<int> atom_lines;

    // The residue map for the frame.
    QMultiMap<QPair<qint64, QString>, qint64> frame_residues;

    // The chain identifier map for the frame.
    QMultiMap<QChar, qint64> frame_chains;

    // Internal function used to parse a single atom line in the file.
    auto parse_atoms = [&](const QString &line, int iatm, int iframe, int iline, int num_lines,
        PDBAtom &frame_atom, QMultiMap<QPair<qint64, QString>, qint64> &local_residues,
        QMultiMap<QChar, qint64> &local_chains, QStringList &errors)
    {
        // Populate atom data.
        frame_atom = PDBAtom(line, errors);

        // Check the next line for additional information.
        if (iline + 1 < num_lines)
        {
            // This is a terminal atom.
            if (lines()[iline + 1].leftRef(6) == "TER   ")
            {
                frame_atom.setTerminal(true);

                // Check whether this is a "real" TER record.
                // Many TER records are used to end group of atoms, e.g. a Tip4P molecule together,
                // rather than a chain.
                // Do we actually care about this? Probably not.
                if (lines()[iline + 1].length() < 27)
                {
                    // Perhaps some kind of warning.
                }
            }
        }

        // Validate the atom data against the first frame.
        if (iframe > 1)
        {
            // Make sure that the first frame has at least 'iatm' entries.
            if (iatm < atoms[0].count())
            {
                // Make sure the atom is consistent between frames.
                if (not validateAtom(frame_atom, atoms[0][iatm]))
                {
                    errors.append(QObject::tr("Invalid record for atom %1 "
                        "in model %2 on line %3. Record doesn't match does not "
                        " match previous model! '%4'")
                        .arg(iatm).arg(iframe).arg(iline).arg(line));
                }
            }
            else
            {
                // Only record this error once per frame.
                if ((iatm + 1) == atoms[0].count())
                {
                    errors.append(QObject::tr("The number of atoms in model "
                        "%1 is larger than the first model!")
                        .arg(iatm));
                }
            }
        }

        // Create residue <number, name> pair.
        QPair<qint64, QString> res(frame_atom.getResNum(), frame_atom.getResName());

        // Update the residue multi-map (residue <number, name> --> atom index).
        local_residues.insert(res, iatm);

        // Insert the chain identifier (ignore if it is blank).
        if (not frame_atom.getChainId().isSpace())
        {
            // Don't duplicate values, only keys.
            if (not frame_chains.contains(frame_atom.getChainId(), frame_atom.getResNum()))
            {
                local_chains.insert(frame_atom.getChainId(), frame_atom.getResNum());
            }
        }

        return frame_atom.isTer();
    };

    // Loop through all lines in the file.
    for (int iline=0; iline<lines().count(); ++iline)
    {
        // Store a reference to the line.
        const QString &line = lines()[iline];

        // Whether to parse atom data at the end of the current loop.
        bool isParse = false;

        // Extract the record type.
        // Could simplify this, i.e. remove whitespace.
        QString record = lines()[iline].left(6);

        // Start of a MODEL record.
        // These are used to define an atom configuratation, so can be used as
        // frames in a trajectory file. Each model entry must be consistent, i.e.
        // it must contain the same number and type of atoms.
        if (record.left(5) == "MODEL")
        {
            imdl++;

            // Extract the model entry number.
            // This should be 4 characters starting at column 11, but we'll assume
            // that it could be the entirety of the rest of the line.
            bool ok;
            int nmod = line.rightRef(line.count() - 6).toInt(&ok);

            // If this failed, try just extracting columns 11-14, as in the PDB standard.
            if (not ok)
                nmod = line.midRef(10,4).toInt(&ok);

            // Check that the model entry number is correct.
            // These must be in ascending order, starting at 1.
            if ((not ok) or (imdl != nmod))
            {
                parse_warnings.append(QObject::tr("Cannot parse the data "
                    "for MODEL %1 has incorrect MODEL entry number '%2'!")
                    .arg(imdl).arg(nmod));

                return;
            }
        }

        // End of a MODEL record.
        else if (record == "ENDMDL")
        {
            if (imdl > 1)
            {
                // Check that the atom number is consisent.
                if (nats != prev_nats)
                {
                    parse_warnings.append(QObject::tr("Cannot parse the data "
                    "for MODEL %1 is not the same size as MODEL %2!")
                    .arg(imdl).arg(imdl-1));
                }
            }

            // Record the number of atoms in the first model entry.
			if (imdl == 1) prev_nats = nats;

            // Flag that a model has been recorded and we can now parse
            // the atom records.
            isParse = true;
        }

        // An ATOM, or HETATM record.
        else if (record == "ATOM  " or
                 record == "HETATM")
        {
            // Store the line number of the atom record.
            atom_lines.append(iline);
            nats++;
        }

        // End of the file.
        if (iline + 1 == lines().count())
            isParse = true;

        // Parse the atom data.
        if (isParse)
        {
            // Don't proceed if there are no more ATOM records to parse.
            if (atom_lines.count() > 0)
            {
                // Initialise atom vector for the frame.
                QVector<PDBAtom> frame_atoms(nats);

                // The number of TER records in the current model.
                int num_ters_model = 0;

                if (usesParallel())
                {
                    QMutex mutex;

                    // Chain-residue multimap for parallel processing.
                    // There will likely be duplicate residue entries for each chain
                    // which will need to be corrected afterwards.
                    QMultiMap<QChar, qint64> temp_chains;

                    tbb::parallel_for( tbb::blocked_range<int>(0, nats),
                                    [&](const tbb::blocked_range<int> &r)
                    {
                        qint64 local_num_ters = 0;
                        QStringList local_errors;
                        QMultiMap<QPair<qint64, QString>, qint64> local_residues;
                        QMultiMap<QChar, qint64> local_chains;

                        for (int i=r.begin(); i<r.end(); ++i)
                        {
                            // Parse the atom record and determine whether it is a terminal record.
                            if (parse_atoms( lines().constData()[atom_lines[i]], i, iframe,
                                atom_lines[i], lines().count(), frame_atoms[i], local_residues,
                                local_chains, local_errors ))
                                    local_num_ters++;
                        }

                        QMutexLocker lkr(&mutex);

                        num_ters_model += local_num_ters;
                        frame_residues += local_residues;
                        temp_chains    += local_chains;
                        parse_warnings += local_errors;
                    });

                    // Remove duplicate residue records from the chains.
                    for (auto chain : temp_chains.uniqueKeys())
                    {
                        // Get a list of the residues that are part of this chain.
                        QList<qint64> chain_residues = temp_chains.values(chain);

                        // Loop over all of the residues.
                        for (auto residue : chain_residues)
                        {
                            // Only insert the residue once.
                            if (not frame_chains.contains(chain, residue))
                                frame_chains.insert(chain, residue);
                        }
                    }
                }
                else
                {
                    for (int i=0; i<nats; ++i)
                    {
                        // Parse the atom record and determine whether it is a terminal record.
                        if (parse_atoms( lines().constData()[atom_lines[i]], i, iframe,
                            atom_lines[i], lines().count(), frame_atoms[i], frame_residues,
                            frame_chains, parse_warnings ))
                                num_ters_model++;
                    }
                }

                // A list of the recorded atom numbers.
                // This is used to sort out issues where there are multiple PDB records
                // with the same atom serial number, which is a common problem.
                QVector<int> atom_numbers;

                // Whether there are duplicate atom numbers.
                bool has_duplicates = false;

                // Check whether there are duplicate atom numbers.
                // If there are, re-number them according to their indices.

                for (int i=0; i<nats; ++i)
                {
                    int number = frame_atoms[i].getSerial();

                    if (not atom_numbers.contains(number))
                    {
                        atom_numbers.append(number);
                    }
                    else
                    {
                        has_duplicates = true;
                        break;
                    }
                }

                // Re-number all of the atoms.
                if (has_duplicates)
                {
                    parse_warnings
                        .append(QObject::tr("Warning: There are duplicate atom "
                        "numbers in the PDB file. Atom numbers have been replaced "
                        "by their index."));

                    for (int i=0; i<nats; ++i)
                        frame_atoms[i].setSerial(i+1);
                }

                // Now check whether the temperature factors are sane.
                // If any exceed 100, then set all to the default of zero.
                bool is_valid = true;

                for (int i=0; i<nats; ++i)
                {
                    if (qAbs(frame_atoms[i].getTemperature()) > 100)
                    {
                        is_valid = false;
                        break;
                    }
                }

                // Zero all of the temperature factors.
                if (not is_valid)
                {
                    for (int i=0; i<nats; ++i)
                        frame_atoms[i].setTemperature(0);
                }

                // Append frame data and clear the vectors.

                atoms.append(frame_atoms);
                atom_lines.clear();

                residues.append(frame_residues);
                frame_residues.clear();

                chains.append(frame_chains);
                frame_chains.clear();

                num_ters.append(num_ters_model);

                iframe++;
                nats = 0;
            }
        }
    }

    this->setScore(nAtoms());
}

/** Helper function used to validate atom data from different model records */
bool PDB2::validateAtom(const PDBAtom &atom1, const PDBAtom &atom2) const
{
    // Different models are typically indexed by atom number, e.g.
    // 1 - nats, nats+1 - 2*nats, ..., or by using the alternate location
    // entry of the atom record, e.g. AALA in the first model, BALA in the second.

    // The following data must be the same for all models (I think...)
    if (atom1.getName()    != atom2.getName()    or
        atom1.getResName() != atom2.getResName() or
        atom1.getChainId() != atom2.getChainId() or
        atom1.getResNum()  != atom2.getResNum()) return false;
    else return true;
}

/** Use the data contained in this parser to create a new System of molecules,
    assigning properties based on the mapping in 'map'. */
System PDB2::startSystem(const PropertyMap &map) const
{
    const int nmols = nModels();

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
    system.setProperty(map["fileformat"].source(), StringProperty(this->formatName()));

    return system;
}

/** Use the data contained in this parser to add information from the file to
    the molecules that exist already in the passed System. For example, this
    may be used to add coordinate data from this file to the molecules in
    the passed System that are missing coordinate data. */
void PDB2::addToSystem(System &system, const PropertyMap &map) const
{
    //you should loop through each molecule in the system and work out
    //which ones are described in the file, and then add data from the file
    //to thise molecules.
}

/** Internal function used to get the molecule structure for molecule 'imol'. */
MolStructureEditor PDB2::getMolStructure(int imol, const PropertyName &cutting) const
{
    // Make sure the frame index is within range.
    if ((imol < 0) or (imol > atoms.count()))
    {
        throw SireError::program_bug(QObject::tr("The molecule index %1 is out of "
            "range, 0 - %2").arg(imol).arg(atoms.count()), CODELOC);
    }

    // Make sure that there are atoms in the molecule.
    if (atoms[imol].count() == 0)
    {
        throw SireError::program_bug(QObject::tr(
            "Strange - there are no atoms in molecule %1?")
                .arg(imol), CODELOC);
    }

    // First step is to build the structure of the molecule, i.e.
    // the layout of cutgroups, residues and atoms.
    MolStructureEditor mol;

    // To do this we'll walk through all of the residues in the frame,
    // adding each atom contained within the residue.

    // Residue index.
    int ires = 0;

    // Create a reverse mapping between residues and chains.
    QMap<qint64, QChar> res_to_chain;

    // Add any chains to the molecule.
    for (auto chain : chains[imol].uniqueKeys())
    {
        mol.add(ChainName(chain));

        // Get a list of the residues that are part of this chain.
        QList<qint64> chain_residues = chains[imol].values(chain);

        // Create the reverse mapping.
        for (auto residue : chain_residues)
            res_to_chain.insert(residue, chain);
    }

    // Loop over all unique residues in the frame.
    for (auto residue : residues[imol].uniqueKeys())
    {
        // Extract the residue number and name.
        auto res_num  = residue.first;
        auto res_name = residue.second;

        // By default we will use one CutGroup per residue.
        // This may be changed later by the cutting system.
        auto cutgroup = mol.add(CGName(QString::number(ires)));

        // Get a sorted list of the atoms that are part of the residue.
        QList<qint64> res_atoms = residues[imol].values(residue);
        qSort(res_atoms);

        // Add the residue to the molecule.
        auto res = mol.add(ResNum(res_num));
        res.rename(ResName(res_name.trimmed()));

        // Reparent the residue to its chain.
        if (res_to_chain.contains(res_num))
            res.reparent(ChainName(res_to_chain[res_num]));

        // Add each atom in the residue to the molecule.
        for (auto res_atom : res_atoms)
        {
            auto atom = cutgroup.add(AtomNum(atoms[imol][res_atom].getSerial()));
            atom.rename(AtomName(atoms[imol][res_atom].getName().trimmed()));

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
MolEditor PDB2::getMolecule(int imol, const PropertyMap &map) const
{
    // At the moment we'll assume that there is a single molecule. Once the molecule
    // is constructed we can use connectivity information to break it into sub-molecules.

    // Make sure the molecule (model) index is within range.
    if ((imol < 0) or (imol > atoms.count()))
    {
        throw SireError::program_bug(QObject::tr("The frame index %1 is out of "
            "range, 0 - %2").arg(imol).arg(atoms.count()), CODELOC);
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
    AtomCoords         coords(molinfo);
    AtomCharges        charges(molinfo);
    AtomElements       elements(molinfo);
    AtomFloatProperty  occupancies(molinfo);
    AtomFloatProperty  temperatures(molinfo);
    AtomStringProperty is_het_atom(molinfo);

    // Now loop through the atoms in the molecule and set each property.
    for (int i=0; i<nAtoms(imol); ++i)
    {
        // Store a reference to the current atom.
        const PDBAtom &atom = atoms[imol][i];

        // Determine the CGAtomIdx for this atom.
        auto cgatomidx = molinfo.cgAtomIdx(AtomNum(atom.getSerial()));

        // Set the properties.
        coords.set(cgatomidx, atom.getCoord());
        charges.set(cgatomidx, int(atom.getCharge()) * SireUnits::mod_electron);
        elements.set(cgatomidx, atom.getElement());
        occupancies.set(cgatomidx, atom.getOccupancy());
        temperatures.set(cgatomidx, atom.getTemperature());

        bool isHet = atom.isHet();

        if (isHet) is_het_atom.set(cgatomidx, "True");
        else       is_het_atom.set(cgatomidx, "False");
    }

    return mol.setProperty(map["coordinates"], coords)
              .setProperty(map["formal-charge"], charges)
              .setProperty(map["element"], elements)
              .setProperty(map["occupancy"], occupancies)
              .setProperty(map["beta-factor"], temperatures)
              .setProperty(map["is-het"], is_het_atom)
              .commit();
}

/** Internal function used to parse a Sire molecule view into a PDB ATOM records using
    the parameters in the property map. */
void PDB2::parseMolecule(const SireMol::Molecule &sire_mol, QVector<QString> &atom_lines,
    QStringList &errors, const SireBase::PropertyMap &map)
{
    // Convert the molecule view to an actual molecule object.
    auto mol = sire_mol.molecule();

    // Store the number of atoms in the molecule.
    int num_atoms = sire_mol.nAtoms();

    // Early exit.
    if (num_atoms == 0) return;

    // TODO: Do we want a hard limit on the number of atoms?
    if (num_atoms > 99999)
    {
        errors.append(QObject::tr("The number of atoms (%1) exceeds "
            " the PDB file format limit (99999)!").arg(num_atoms));

        return;
    }

    // Store the number of chains in the molecule.
    int num_chains = sire_mol.nChains();

    // Resize the data record containers (one TER record for each chain).
    atom_lines.resize(num_atoms + num_chains);

    // Whether each atom is a terminal atom, i.e. the end of a chain.
    QVector<bool> is_ter(num_atoms);
    is_ter.fill(false);

    // Loop over the chains.
    for (int i=0; i<num_chains; ++i)
    {
        // Extract the chain.
        auto chain = sire_mol.chain(ChainIdx(i));

        // Extract the atoms from the chain.
        auto atoms = chain.atoms();

        // Extract the number of the last atom in the chain
        int terminal_atom = atoms[atoms.count()-1]
                           .read()
                           .asA<SireMol::Atom>()
                           .index()
                           .value();

        // Set the terminal atom.
        is_ter[terminal_atom] = true;
    }

    if (usesParallel())
    {
        QMutex mutex;

        // Local data storage.
        QVector<PDBAtom> local_atoms(num_atoms);

        tbb::parallel_for( tbb::blocked_range<int>(0, num_atoms),
                        [&](const tbb::blocked_range<int> &r)
        {
            // Create local data objects.
            QStringList local_errors;

            // Convert each atom into a PDBAtom object
            // and generate a PDB data record.
            for (int i=r.begin(); i<r.end(); ++i)
            {
                local_atoms[i] = PDBAtom(sire_mol.atom(AtomIdx(i)), is_ter[i], local_errors);
                atom_lines[i] = local_atoms[i].toPDBRecord();
            }

            if (not local_errors.isEmpty())
            {
                // Acquire a lock.
                QMutexLocker lkr(&mutex);

                // Update the warning messages.
                errors += local_errors;
            }
        });

        // Create TER records if the system contains chains.
        if (num_chains > 0)
        {
            // Make a copy of the atom lines.
            auto lines = atom_lines;

            // Line index.
            int iline = 0;

            // Now loop through the atoms and insert TER records as needed.
            // This has to be done in serial.
            for (int i=0; i<num_atoms; ++i)
            {
                // Copy the atom record across.
                atom_lines[iline] = lines[i];
                iline++;

                // Add a TER record for this atom.
                if (is_ter[i])
                {
                    atom_lines[iline] = QString("TER   %1      %2 %3\%4\%5")
                                            .arg(lines[i].mid(6, 5))
                                            .arg(lines[i].mid(17, 3))
                                            .arg(lines[i].at(21))
                                            .arg(lines[i].mid(22, 4))
                                            .arg(lines[i].at(26));

                    iline++;
                }
            }
        }
    }
    else
    {
        // Line index.
        int iline = 0;

        // Loop over all of the atoms.
        for (int i=0; i<num_atoms; ++i)
        {
            // Initalise a PDBAtom.
            PDBAtom atom(sire_mol.atom(AtomIdx(i)), is_ter[i], errors);

            // Generate a PDB atom data record.
            atom_lines[iline] = atom.toPDBRecord();
            iline++;

            // Add a TER record for this atom.
            if (is_ter[i])
            {
                atom_lines[iline] = QString("TER   %1      %2 %3\%4\%5")
                                        .arg(atom_lines[iline-1].mid(6, 5))
                                        .arg(atom_lines[iline-1].mid(17, 3))
                                        .arg(atom_lines[iline-1].at(21))
                                        .arg(atom_lines[iline-1].mid(22, 4))
                                        .arg(atom_lines[iline-1].at(26));

                iline++;
            }
        }
    }
}
