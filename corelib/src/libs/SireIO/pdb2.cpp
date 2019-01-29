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
#include "SireMol/errors.h"
#include "SireMol/molecule.h"
#include "SireMol/moleditor.h"

#include "SireUnits/units.h"

#include <QFile>
#include <QtMath>

using namespace SireBase;
using namespace SireIO;
using namespace SireMol;
using namespace SireStream;
using namespace SireSystem;
using namespace SireUnits;

const RegisterParser<PDB2> register_pdb;
static const RegisterMetaType<PDB2> r_pdb2;
static const RegisterMetaType<PDBAtom> r_pdbatom(NO_ROOT);

QDataStream &operator<<(QDataStream &ds, const PDBAtom &pdbatom)
{
    writeHeader(ds, r_pdbatom, 1);

    SharedDataStream sds(ds);

    sds << pdbatom.record << pdbatom.serial << pdbatom.name << pdbatom.alt_loc
        << pdbatom.res_name << pdbatom.chain_id << pdbatom.res_num << pdbatom.insert_code
        << pdbatom.coord << pdbatom.occupancy << pdbatom.temperature << pdbatom.element
        << pdbatom.charge << pdbatom.is_het << pdbatom.is_ter;

    return ds;
}

QDataStream &operator>>(QDataStream &ds, PDBAtom &pdbatom)
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

QDataStream &operator<<(QDataStream &ds, const PDB2 &pdb2)
{
    writeHeader(ds, r_pdb2, 1);

    SharedDataStream sds(ds);

    sds << pdb2.atoms << pdb2.residues << pdb2.chains
        << pdb2.parse_warnings << static_cast<const MoleculeParser&>(pdb2);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, PDB2 &pdb2)
{
    VersionID v = readHeader(ds, r_pdb2);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> pdb2.atoms >> pdb2.residues >> pdb2.chains
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
    serial = line.midRef(6,5).toInt(&ok);

    if (not ok)
    {
        errors.append(QObject::tr("Cannot extract the atom serial number "
            "from part (%1) from line '%2'").arg(line.mid(6,5)).arg(line));

        return;
    }

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
    res_num = line.midRef(22,4).toInt(&ok);

    if (not ok)
    {
        errors.append(QObject::tr("Cannot extract the residue sequence number "
            "from part (%1) from line '%2'").arg(line.mid(22,4)).arg(line));

        return;
    }

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
        if ((occupancy < 0) or
            (occupancy > 1))
        {
            errors.append(QObject::tr("The occupancy (%1) was out of range! "
                "Setting to default value of 1.0.")
                .arg(occupancy));

            occupancy = 1;
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
        }
    }

    // Extract the element name.
    element = line.mid(76,2);

    // If the element is empty, try to guess from the atom name.
    if (element.simplified().isEmpty())
        element = Element(name).symbol()[0];

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
PDBAtom::PDBAtom(const SireMol::Atom &atom, bool is_ter, const PropertyMap &map, QStringList &errors) :
    serial(atom.number().value()),
    name(atom.name().value().toUpper()),
    occupancy(1.0),
    temperature(0.0),
    element("X"),
    charge(0),
    is_het(false),
    is_ter(is_ter)
{
    // The atom must have atomic coordinates to be valid.
    if (not atom.hasProperty(map["coordinates"]))
    {
        errors.append(QObject::tr("The atom does not have coordinates!"));

        return;
    }

    // Extract the atomic coordinates.
    coord = atom.property<SireMaths::Vector>(map["coordinates"]);

    // The atom is within a residue.
    if (atom.isWithinResidue())
    {
        res_name = atom.residue().name().value();
        res_num  = atom.residue().number().value();

        // Optional insertion code property.
        if (atom.residue().hasProperty("insert_code"))
            insert_code = atom.residue().property<QString>(map["insert_code"])[0];
    }

    // The atom is within a chain.
    if (atom.isWithinChain())
    {
        // TODO: Make sure we can handle situations where the
        // chain ID is a string.
        chain_id = atom.chain().name().value().at(0);
    }

    // Extract the occupancy.
    if (atom.hasProperty(map["occupancy"]))
    {
        occupancy = atom.property<double>(map["occupancy"]);
    }

    // Extract the temperature factor.
    if (atom.hasProperty(map["beta_factor"]))
    {
        temperature = atom.property<double>(map["beta_factor"]);
    }

    // Extract the element name.
    if (atom.hasProperty(map["element"]))
    {
        element = atom.property<Element>(map["element"]).symbol().toUpper();
    }
    // Otherwise, try to guess from the atom name.
    else
    {
        element = Element(name).symbol().toUpper();
    }

    // Extract the atomic charge.
    if (atom.hasProperty(map["formal_charge"]))
    {
        // TODO: This doesn't seem to be working in all cases.
        charge = atom.property<SireUnits::Dimension::Charge>(map["formal_charge"]).value();
        charge /= SireUnits::mod_electron;
    }

    // Determine whether this is a HETATM.
    if (atom.hasProperty(map["is_het"]))
    {
        if (atom.property<QString>(map["is_het"]) == "True")
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
void PDBAtom::setTerminal(bool is_ter)
{
    this->is_ter = is_ter;
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
qint64 PDBAtom::getNumber() const
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
QChar PDBAtom::getChainID() const
{
    return chain_id;
}

/** Set the chain id. */
void PDBAtom::setChainID(QChar id)
{
    chain_id = id;
}

/** Get the residue sequence number. */
qint64 PDBAtom::getResNum() const
{
    return res_num;
}

/** Set the residue sequence number. */
void PDBAtom::setResNum(int num)
{
    res_num = num;
}

/** Get the residue index. */
qint64 PDBAtom::getResIdx() const
{
    return res_idx;
}

/** Set the residue sequence index. */
void PDBAtom::setResIdx(int idx)
{
    res_idx = idx;
}

/** Get the residue insertion code. */
QChar PDBAtom::getInsertCode() const
{
    return insert_code;
}

/** Get the atomic coordinates. */
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

    // Write all records as ATOM, i.e. convert HETATM to ATOM.
    line.append("ATOM  ");

    // Append the atom serial number. Since we enforce a hard
    // limit of 99999 atoms, we can assume that the serial is
    // 5 characters or less.
    // TODO: Handle non-standard serial numbers.
    line.append(QString("%1").arg(serial, 5));

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

    // Make a copy of the res_num and insert_code member data.
    auto local_res_num = res_num;
    auto local_insert_code = insert_code;

    // Make sure the residue number doesn't exceed 9999.
    // If it does, the we set the number to mod(resnum, 9999) and use the insertion
    // code character to identify unique residues.
    if (local_res_num > 9999)
    {
        // Upper-case alphabet and digits 0 through 9.
        QString characters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";

        // Work out the insertion code index (the multiple of 9999).
        int index = qFloor(local_res_num / 9999.0) - 1;

        // Cap the character index.
        if (index > 35)
            index = 35;

        // Update the insertion code.
        local_insert_code = characters[index];

        // Wrap the residue number.
        local_res_num = local_res_num % 9999;
    }

    // Append the residue sequence number.
    line.append(QString("%1").arg(local_res_num, 4));

    // Append the residue insertion code
    if (not local_insert_code.isNull()) line.append(local_insert_code);
    else                                line.append(" ");

    // Append the atomic coordinates.
    line.append(QString("   %1\%2\%3")
        .arg(coord[0], 8, 'f', 3)
        .arg(coord[1], 8, 'f', 3)
        .arg(coord[2], 8, 'f', 3));

    // Append the occupancy.
    line.append(QString("%1").arg(occupancy, 6, 'f', 2));

    // Append the beta factor.
    line.append(QString("%1").arg(temperature, 6, 'f', 2));

    // Append the element (right-justified).
    auto el = element.left(2);
    if (el.count() == 1) el = QString(" %1").arg(el);
    line.append(QString("          %1").arg(el));

    // Apped the atomic charge, truncating as necessary.
    if (charge != 0)
    {
        if (charge < 0) line.append(QString("%1-").arg(qAbs(charge)));
        else            line.append(QString("%1+").arg(charge));
    }

    return line;
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

    // Check for velocities.
    // All molecules must have velocity properties to be able to generate a PDB
    // velocity file for NAMD.
    for (int i=0; i<nmols; ++i)
    {
        // Extract the molecule.
        auto mol = system[molnums[i]].molecule();

        // All molecules must have velocities.
        if (not mol.hasProperty(map["velocity"]))
        {
            velocities.clear();
            break;
        }

        // Loop over all atoms in the molecule.
        for (int j=0; j<mol.nAtoms(); j++)
        {
            // Extract the atom.
            auto atom = mol.atom(AtomIdx(j));

            // All atoms in the molecule must have velocities.
            if (not atom.hasProperty(map["velocity"]))
            {
                velocities.clear();
                i = nmols;
                break;
            }

            // Append the velocity to the vector.
            else
            {
                velocities.append(atom.property<SireMol::Velocity3D>(map["velocity"]));
            }
        }
    }

    // The list of lines.
    QStringList lines;

    // Lines for different PDB data records (one for each molecule).
    QVector<QVector<QString> > atom_lines(nmols);

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

    // Whether the object contains MODEL records (trajectory frames).
    bool is_model = isModel(system);

    // Now assemble the lines from the record data for each molecule.
    // We do this in serial since the order matters.
    for (int i=0; i<nmols; ++i)
    {
        // MODEL record.
        if (is_model)
            lines.append(QString("MODEL     %1").arg(i+1));

        // ATOM records.
        lines.append(atom_lines[i].toList());

        // Record the end of the molecule.
        if (is_model) lines.append("ENDMDL");
        else          lines.append("TER");
    }

    lines.append("END");

    // Reparse the lines as a self-consistency check.
    PDB2 parsed(lines, map);

    // Copy the velocities vector to the new object.
    parsed.velocities = velocities;

    this->operator=(parsed);
}

/** Copy constructor */
PDB2::PDB2(const PDB2 &other) :
    ConcreteProperty<PDB2,MoleculeParser>(other),
    atoms(other.atoms),
    chains(other.chains),
    residues(other.residues),
    velocities(other.velocities),
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
        this->chains = other.chains;
        this->residues = other.residues;
        this->velocities = other.velocities;
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
    if (lines().isEmpty())
        return QObject::tr("PDB2::null");
    else
    {
        return QObject::tr("PDB2( nMolecules() = %1, "
            "nChains() = %2, nResidues() = %3, nAtoms() = %4 )")
            .arg(nMolecules()).arg(nChains()).arg(nResidues()).arg(nAtoms());
    }
}

/** Convert the parsed data to a collection of PDB record lines. */
QVector<QString> PDB2::toLines(bool is_velocity) const
{
    // Store the number of molecules.
    const int nmols = nMolecules();

    // No molecules.
    if (nmols == 0)
        return QVector<QString>();

    // The list of lines.
    QVector<QString> lines;

    // Whether the object contains MODEL records (trajectory frames).
    bool is_model = isModel();

    // The index offset for the velocities vector.
    int offset = 0;

    // Now assemble the lines from the record data for each molecule.
    // We do this in serial since the order matters.
    for (int i=0; i<nmols; ++i)
    {
        // MODEL record.
        if (is_model)
            lines.append(QString("MODEL     %1").arg(i+1));

        // The number of atoms for this model.
        const int num_atoms = nAtoms(i);

        // The number of chains for this model.
        const int num_chains = nChains(i);

        // The atoms lines for this model.
        QVector<QString> atom_lines(num_atoms + num_chains);

        if (usesParallel())
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, num_atoms),
                            [&](const tbb::blocked_range<int> r)
            {
                for (int j=r.begin(); j<r.end(); ++j)
                {
                    atom_lines[j] = atoms[i][j].toPDBRecord();

                    // We are writing a PDB "velocity" file for NAMD.
                    if (is_velocity)
                    {
                        // Create the velocity string.
                        QString vel_string(QString("   %1\%2\%3")
                            .arg(velocities[offset+j][0].value(), 8, 'f', 3)
                            .arg(velocities[offset+j][1].value(), 8, 'f', 3)
                            .arg(velocities[offset+j][2].value(), 8, 'f', 3));

                        // Replace the coordinate record data with velocities.
                        atom_lines[j].replace(30, 24, vel_string);
                    }
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
                for (int j=0; j<num_atoms; ++j)
                {
                    // Copy the atom record across.
                    atom_lines[iline] = lines[j];
                    iline++;

                    // Add a TER record for this atom.
                    if (atoms[i][j].isTer())
                    {
                        atom_lines[iline] = QString("TER   %1      %2 %3\%4\%5")
                                                .arg(lines[j].mid(6, 5))
                                                .arg(lines[j].mid(17, 3))
                                                .arg(lines[j].at(21))
                                                .arg(lines[j].mid(22, 4))
                                                .arg(lines[j].at(26));

                        iline++;
                    }
                }
            }
        }
        else
        {
            // Line index.
            int iline = 0;

            for (int j=0; j<num_atoms; ++j)
            {
                atom_lines[iline] = atoms[i][j].toPDBRecord();

                // We are writing a PDB "velocity" file for NAMD.
                if (is_velocity)
                {
                    // Create the velocity string.
                    QString vel_string(QString("   %1\%2\%3")
                        .arg(velocities[offset+j][0].value(), 8, 'f', 3)
                        .arg(velocities[offset+j][1].value(), 8, 'f', 3)
                        .arg(velocities[offset+j][2].value(), 8, 'f', 3));

                    // Replace the coordinate record data with velocities.
                    atom_lines[j].replace(30, 24, vel_string);
                }

                iline++;

                // Add a TER record for this atom.
                if (atoms[i][j].isTer())
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

        // ATOM records.
        lines += atom_lines;

        // Record the end of this molecule.
        if (is_model) lines.append("ENDMDL");
        else          lines.append("TER");

        // Update the offset.
        offset += num_atoms;
    }

    lines.append("END");

    return lines;
}

/** Write a velocity file in PDB format. This can be used as a restart for NAMD simulations. */
bool PDB2::writeVelocityFile(const QString &filename) const
{
    if (velocities.isEmpty())
        return false;

    // Generate the vector of record lines.
    // Here atomic coordinates are replaced by velocities.
    QVector<QString> lines = toLines(true);

    if (lines.isEmpty())
        return false;

    QFile f(filename);

    if (not f.open( QIODevice::WriteOnly | QIODevice::Text ))
    {
        throw SireError::file_error(f, CODELOC);
    }

    QTextStream ts(&f);

    for (const QString &line : lines)
    {
        ts << line << '\n';
    }

    f.close();

    return true;
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

/** Return whether or not this is a lead parser. The lead parser is responsible
    for starting the process of turning the parsed file into the System. There
    must be one and one-only lead parser in a set of parsers creating a System */
bool PDB2::isLead() const
{
    return true;
}

/** Return the number of models (molecules). */
int PDB2::nMolecules() const
{
    return atoms.count();
}

/** Return the total number of residues. */
int PDB2::nResidues() const
{
    int num_residues = 0;

    for (int i=0; i<residues.count(); ++i)
        num_residues += residues[i].uniqueKeys().count();

    return num_residues;
}

/** Return the number of residues in molecule 'i'. */
int PDB2::nResidues(int i) const
{
    return residues[i].uniqueKeys().count();
}

/** Return the total number of chains. */
int PDB2::nChains() const
{
    int num_chains = 0;

    for (int i=0; i<chains.count(); ++i)
        num_chains += chains[i].uniqueKeys().count();

    return num_chains;
}

/** Return the number of chains in molecule 'i'. */
int PDB2::nChains(int i) const
{
    return chains[i].uniqueKeys().count();
}

/** Return the total number of atoms. */
int PDB2::nAtoms() const
{
    int num_atoms = 0;

    for (int i=0; i<atoms.count(); ++i)
        num_atoms += atoms[i].count();

    return num_atoms;
}

/** Return the number of atoms in molecule 'i'. */
int PDB2::nAtoms(int i) const
{
    return atoms[i].count();
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

    // Number of atoms in the previous model.
    // Each MODEL section should contain the same number of atoms.
    // For each MODEL record there should be a corresponding ENDMDL.
    int prev_nats = 0;

    // The indices for lines containing atom data.
    // For each molecule we identify the lines containing atom data, then parse them.
    QVector<int> atom_lines;

    // The chain mapping for the molecule.
    QMultiMap<QChar, qint64> mol_chains;

    // The residue mapping for the molecule.
    QMultiMap<qint64, qint64> mol_residues;

    // Internal function used to parse a single atom line in the file.
    auto parse_atoms = [&](const QString &line, int iatm,
        int iline, int num_lines, PDBAtom &atom, QStringList &errors)
    {
        // Populate atom data.
        atom = PDBAtom(line, errors);

        // Check the next line for additional information.
        if (iline + 1 < num_lines)
        {
            // Check whether this is a "real" TER record.
            // Standalone TER records are used to end molecules, e.g. a Tip4P molecule.
            if ((lines()[iline + 1].leftRef(6) == "TER   ") and
                (lines()[iline + 1].simplified() != "TER"))
            {
                atom.setTerminal(true);
            }
        }
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
                parse_warnings.append(QObject::tr("Incorrect MODEL entry "
                    "number. Found '%1', should be '%2'")
                    .arg(nmod).arg(imdl));

                nmod = imdl;
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

        // A standalone TER record.
        // This is used to flag the end of a molecule.
        else if (lines()[iline].simplified() == "TER")
            isParse = true;

        // End of the file.
        if (iline + 1 == lines().count())
            isParse = true;

        // Parse the atom data.
        if (isParse)
        {
            // Don't proceed if there are no more ATOM records to parse.
            if (atom_lines.count() > 0)
            {
                // Initialise atom vector for the molecule.
                QVector<PDBAtom> mol_atoms(nats);

                if (usesParallel())
                {
                    QMutex mutex;

                    tbb::parallel_for(tbb::blocked_range<int>(0, nats),
                                    [&](const tbb::blocked_range<int> &r)
                    {
                        QStringList local_errors;

                        for (int i=r.begin(); i<r.end(); ++i)
                        {
                            // Parse the atom record.
                            parse_atoms(lines().constData()[atom_lines[i]], i, atom_lines[i],
                                lines().count(), mol_atoms[i], local_errors);
                        }

                        if (not local_errors.isEmpty())
                        {
                            // Acquire a lock.
                            QMutexLocker lkr(&mutex);

                            // Update the global error messages.
                            parse_warnings += local_errors;
                        }
                    });
                }
                else
                {
                    for (int i=0; i<nats; ++i)
                    {
                        // Parse the atom record.
                        parse_atoms(lines().constData()[atom_lines[i]], i, atom_lines[i],
                            lines().count(), mol_atoms[i], parse_warnings);
                    }
                }

                /* We now attempt to the following common PDB errors:

                   1) Duplicate atom names.
                   1) Duplicate chain identifiers.
                   3) Duplicate residue numbers.
                   4) Non-sensical temperature factors.
                 */

                /**************** 1) FIX ATOM NUMBERS ****************/

                // Check whether there are duplicate atom numbers.
                // If there are, re-number them according to their indices.

                // A vector of the recorded atom numbers.
                QVector<int> atom_numbers;

                bool ok = true;
                for (int i=0; i<nats; ++i)
                {
                    int num = mol_atoms[i].getNumber();

                    if (not atom_numbers.contains(num))
                    {
                        atom_numbers.append(num);
                    }
                    else
                    {
                        ok = false;
                        break;
                    }
                }

                // There were duplicates: Re-number all of the atoms.
                if (not ok)
                {
                    parse_warnings.append(QObject::tr("Warning: There are duplicate atom "
                        "numbers in the PDB file. Atom numbers have been replaced "
                        "by their index."));

                    for (int i=0; i<nats; ++i)
                        mol_atoms[i].setSerial(i+1);
                }

                /**************** 2) FIX CHAIN IDS ****************/

                /* Check whether there are duplicate chain idenfitiers.
                   If there are, re-label them using alphanumberic characters,
                   first the alphabet in upper case, followed by the alphabet
                   in lower case, then digits 0-9.
                 */

                // A vector of the recorded chain identifiers.
                QVector<QChar> chain_ids;

                QChar curr_id = mol_atoms[0].getChainID();

                // Insert the first chain identifier.
                if (not curr_id.isSpace())
                    chain_ids.append(curr_id);

                ok = true;
                for (int i=1; i<nats; ++i)
                {
                    QChar id = mol_atoms[i].getChainID();

                    // We've reached the next chain record.
                    if ((not id.isSpace()) and (id != curr_id))
                    {
                        // Check whether we've already seen this chain.
                        if (not chain_ids.contains(id))
                        {
                            chain_ids.append(id);
                            curr_id = id;
                        }
                        else
                        {
                            // Often HETATM records can be at the end of the
                            // file, far away from the rest of the chain record.
                            if (not mol_atoms[i].isHet())
                            {
                                ok = false;
                                break;
                            }

                            /* TODO:
                               Because of the above, we would no longer be able
                               to write to file in order of chain, then residue.
                               This means that isolated HETATM records will appear
                               at the end of the file, which might be the TER
                               entry for their chain (which appears higher up in
                               the file).

                               The only way to fix this would be to insert HETATM
                               records into the appropriate place, i.e. at the end
                               of the chain where it first appears, although this
                               in turn can cause problems with atom numbering,
                               although these could be updated afterwards.
                             */
                        }
                    }
                }

                // There were duplicates: re-label all of the chains.
                if (not ok)
                {
                    parse_warnings.append(QObject::tr("Warning: There are duplicate chain "
                        "identifiers in the PDB file. Chains have been relabelled!"));

                    int num_chains = 0;

                    // The current chain identifier.
                    QChar curr_id = mol_atoms[0].getChainID();

                    // List of all possible chain identifiers.
                    QString identifiers("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789");

                    for (int i=0; i<nats; ++i)
                    {
                        if (not mol_atoms[i].getChainID().isSpace())
                        {
                            // Get the chain identifier for the atom.
                            QChar id = mol_atoms[i].getChainID();

                            if ((not id.isSpace()) and (id != curr_id))
                            {
                                // Increment the chain count and check that it
                                // doesn't exceed the PDF file limit.
                                num_chains++;

                                if (num_chains > 62)
                                {
                                    parse_warnings.append(QObject::tr("Warning: There number "
                                        "of chains exceeds the limit for the PDB format (62)."));

                                    return;
                                }

                                // Update the current chain identifier.
                                curr_id = mol_atoms[i].getChainID();
                            }

                            // Update the chain idenfifier.
                            mol_atoms[i].setChainID(identifiers[num_chains]);
                        }
                    }
                }

                /************* 3) FIX RESIDUE NUMBERS *************/

                /* Check whether there are duplicate residue numbers,
                   if so, we need to re-number all of the residues in
                   ascending order.

				   However, residue numbers are not unique across
                   chains, so we need to preserved the numbering if
                   the chain ID differs.
                 */

                /* A hash betwen residue numbers and a residue name, chain identifier pair.
                   Note that we create a string representation of the residue number by
                   appending the insertion code. This ensures a unique number if a residue
                   has been inserted into the PDB.
                 */
                QHash<QString, QPair<QString, QChar>> res_hash;

                // A multimap between duplicate residues and the atoms in those residues.
                QMultiMap<int, int> duplicates;

                // Initalise the maximum residue number.
                int max_res_num = -1000000;

                for (int i=0; i<nats; ++i)
                {
                    QString res_name = mol_atoms[i].getResName();
                    int     res_num  = mol_atoms[i].getResNum();
                    QChar   icode    = mol_atoms[i].getInsertCode();
                    QChar   chain_id = mol_atoms[i].getChainID();

                    // Create a string out of the residue number and insertion code.
                    QString num(QString("%1%2").arg(res_num).arg(icode));

                    // Check if this residue number exceeds the current maximum.
                    if (res_num > max_res_num)
                        max_res_num = res_num;

                    // We've already seen this residue number.
                    if (res_hash.contains(num))
                    {
                        // The previous name was different and the residue is
                        // part of the same chain.
                        if ((res_name != res_hash[num].first) and
                            (chain_id == res_hash[num].second))
                        {
                            // Insert the atom into the duplicate resiude multi-map.
                            duplicates.insert(res_num, i);
                        }

                        // Add the residue to the hash.
                        else
                        {
                            res_hash.insert(num,
                                QPair<QString, QChar>(res_name, chain_id));
                        }
                    }

                    // Add the residue to the hash.
                    else
                    {
                        res_hash.insert(num,
                            QPair<QString, QChar>(res_name, chain_id));
                    }
                }

                // Re-number all of the residues.
                if (not duplicates.isEmpty())
                {
                    // The incorrect residues are re-numbered starting at a value of one
                    // above the maximum residue number found in the PDB file. This means
                    // that we preserve the numbers of the correct residues. However, the
                    // residue numbers in the file will now be out of sequence.

                    parse_warnings.append(QObject::tr("Warning: There are duplicate residue "
                        "numbers in the PDB file. Residue numbers have been updated."));

                    // The residue counter. Start at one above the current maximum.
                    int num_res = max_res_num + 1;

                    // Loop over all duplicated residue numbers in the molecule.
                    for (const auto &res_num : duplicates.uniqueKeys())
                    {
                        // Loop over all of the atoms in the residue.
                        for (const auto &atom_num : duplicates.values(res_num))
                        {
                            // Update the residue number.
                            mol_atoms[atom_num].setResNum(num_res);
                        }

                        // Increment the residue number.
                        num_res++;
                    }
                }

                // Now we need to loop through all of the atoms and set a residue "index".
                // This will help with breaking the molecule up into its constituent parts.

                // A string identifying the current residue.
                // name + number + insert_code + chain.
                QString res_string(QString("%1%2%3%4")
                    .arg(mol_atoms[0].getResName())
                    .arg(mol_atoms[0].getResNum())
                    .arg(mol_atoms[0].getInsertCode())
                    .arg(mol_atoms[0].getChainID()));

                // The current residue index.
                int res_idx = 0;

                // Set the index for the residue of the first atom.
                mol_atoms[0].setResIdx(0);

                // A has between the residue string and its index.
                QHash<QString, int> res_indices;

                // Add the first residue to the hash.
                res_indices[res_string] = res_idx;

                // Loop through the rest of the residues.
                for (int i=1; i<nats; ++i)
                {
                    // A string identifying the current residue.
                    QString res_string(QString("%1%2%3%4")
                        .arg(mol_atoms[i].getResName())
                        .arg(mol_atoms[i].getResNum())
                        .arg(mol_atoms[i].getInsertCode())
                        .arg(mol_atoms[i].getChainID()));

                    // This residue has already been added.
                    if (res_indices.contains(res_string))
                    {
                        // Set the residue index to the hash value.
                        mol_atoms[i].setResIdx(res_indices[res_string]);
                    }

                    // This is a new residue.
                    else
                    {
                        // Increment the residue index.
                        res_idx++;

                        // Set the residue index.
                        mol_atoms[i].setResIdx(res_idx);

                        // Add the new residue to the hash.
                        res_indices[res_string] = res_idx;
                    }
                }

                // Now check whether the temperature factors are sane.
                // If any exceed 100, then set all to the default of zero.

                ok = true;
                for (int i=0; i<nats; ++i)
                {
                    if (qAbs(mol_atoms[i].getTemperature()) > 100)
                    {
                        ok = false;
                        break;
                    }
                }

                // Zero all of the temperature factors.
                if (not ok)
                {
                    parse_warnings.append(QObject::tr("Warning: Invalid temperature "
                        "factor found. All values have been zeroed."));

                    for (int i=0; i<nats; ++i)
                        mol_atoms[i].setTemperature(0);
                }

                /* Now used the parsed data to construct the residue and chain
                   information for the molecule.

                   Each chain contains residues, and the residues contain atoms.
                   However, residues within different chains may have the same
                   same name and number.
                 */

                if (usesParallel())
                {
                    QMutex mutex;

                    tbb::parallel_for(tbb::blocked_range<int>(0, nats),
                                    [&](const tbb::blocked_range<int> &r)
                    {
                        // Local data objects.
                        QMultiMap<QChar, qint64>  local_chains;
                        QMultiMap<qint64, qint64> local_residues;

                        for (int i=r.begin(); i<r.end(); ++i)
                        {
                            // Store a reference to the atom.
                            const PDBAtom &atom = mol_atoms[i];

                            // Map the chain identifier to the residue index.
                            if (not atom.getChainID().isSpace())
                                local_chains.insert(atom.getChainID(), atom.getResIdx());

                            // Map the residue index to the atom index.
                            local_residues.insert(atom.getResIdx(), i);
                        }

                        QMutexLocker lkr(&mutex);

                        mol_chains   += local_chains;
                        mol_residues += local_residues;
                    });
                }
                else
                {
                    for (int i=0; i<nats; ++i)
                    {
                        // Store a reference to the atom.
                        const PDBAtom &atom = mol_atoms[i];

                        // Map the chain identifier to the residue index.
                        if (not atom.getChainID().isSpace())
                            mol_chains.insert(atom.getChainID(), atom.getResIdx());

                        // Map the residue index to the atom index.
                        mol_residues.insert(atom.getResIdx(), i);
                    }
                }

                // Finally, append molecule data and clear the vectors.

                atoms.append(mol_atoms);
                atom_lines.clear();

                chains.append(mol_chains);
                mol_chains.clear();

                residues.append(mol_residues);
                mol_residues.clear();

                nats = 0;
            }
        }
    }

    this->setScore(nAtoms());
}

/** Helper function used to determine whether the object contains multiple models. */
bool PDB2::isModel() const
{
    // Different MODEL records are typically used to index frames in
    // a trajectory. In this case, the only difference between the records
    // in different models is the atom coordinates.

    // The number of atoms in the first molecule.
    int nats = nAtoms(0);

    // Extract the first atom from the first molecule.
    PDBAtom atom1 = atoms[0].at(0);

    for (int i=1; i<nMolecules(); ++i)
    {
        // Number of atoms don't match.
        if (nAtoms(i) != nats) return false;
        else
        {
            // Extract the first atom from the molecule.
            PDBAtom atom2 = atoms[i].at(0);

            // The following data must be the same for all models (I think...)
            if (atom1.getName()    != atom2.getName()    or
                atom1.getResName() != atom2.getResName() or
                atom1.getChainID() != atom2.getChainID() or
                atom1.getResNum()  != atom2.getResNum()) return false;
        }
    }

    // If we get this far, then we have a file with MODELS.
    return true;
}

/** Helper function used to determine whether the Sire system contains multiple models. */
bool PDB2::isModel(const SireSystem::System &system) const
{
    // Get the MolNums of each molecule in the System - this returns the
    // numbers in MolIdx order.
    const QVector<MolNum> molnums = system.getMoleculeNumbers().toVector();

    // The number of molecules in the system.
    int nmols = molnums.count();

    // The number of atoms in the first molecule.
    int nats = system[molnums[0]].molecule().nAtoms();

    // Extract the first atom from the first molecule.
    auto atom1 = system[molnums[0]]
                .molecule()
                .atoms()[0]
                .read()
                .asA<SireMol::Atom>();

    for (int i=1; i<nmols; ++i)
    {
        // Number of atoms don't match.
        if (system[molnums[i]].molecule().nAtoms() != nats) return false;
        else
        {
            // Extract the first atom from the molecule.
            auto atom2 = system[molnums[i]]
                        .atoms()[0]
                        .read()
                        .asA<SireMol::Atom>();

            // Check the atom name and number.
            if (atom1.name()    != atom2.name()     or
                atom1.number()  != atom2.number()   or
                atom1.residue() != atom2.residue()) return false;
        }
    }

    // If we get this far, then we have a system with MODELS.
    return true;
}

/** Use the data contained in this parser to create a new System of molecules,
    assigning properties based on the mapping in 'map'. */
System PDB2::startSystem(const PropertyMap &map) const
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

/** Use the data contained in this parser to add information from the file to
    the molecules that exist already in the passed System. For example, this
    may be used to add coordinate data from this file to the molecules in
    the passed System that are missing coordinate data. */
void PDB2::addToSystem(System &system, const PropertyMap &map) const
{
    /* Here we add information to an existing system that has been
       generated by a lead parser. In particular, this will be used
       to add coordinate information to a system generated from a
       CHARMM Protein Structure File (PSF). In this case, we'll need
       to make sure that the system has consistent atom information.
     */

    // Get the MolNums of each molecule in the System - this returns the
    // numbers in MolIdx order.
    const QVector<MolNum> molnums = system.getMoleculeNumbers().toVector();

    // Store the number of molecules.
    const int num_mols = molnums.count();

    // No molecules in the system.
    if (num_mols == 0)
        return;

    // A vector of the updated molecules.
    QVector<Molecule> mols(num_mols);

    // Loop over the molecules to work out the total number of atoms.
    int num_atoms = 0;
    for (int i=0; i<num_mols; ++i)
    {
        num_atoms += system[molnums[i]].molecule().nAtoms();
    }

    // The number of atoms must match.
    if (nAtoms() != num_atoms)
    {
        throw SireError::incompatible_error(QObject::tr("The number of atoms "
            "does not match the passed system. Found %1, expected %2")
            .arg(nAtoms()).arg(num_atoms), CODELOC);
    }


    // Whether each atom from the PDB record has been used as
    // a match for the system constructed by the lead parser.
    // This is used to avoid duplicate matches.
    QVector<bool> used_atoms(nAtoms(), false);

    // We now need to iterate over the molecules in the passed system,
    // find the matching atoms from the PDB data and add coordinate properties.
    if (usesParallel())
    {
        tbb::parallel_for(tbb::blocked_range<int>(0, num_mols),
                           [&](const tbb::blocked_range<int> r)
        {
            for (int i=r.begin(); i<r.end(); ++i)
            {
                // Add coordinate information to the molecule.
                mols[i] = updateMolecule(system[molnums[i]].molecule(), used_atoms, map);
            }
        });
    }
    else
    {
        for (int i=0; i<num_mols; ++i)
        {
            // Add coordinate information to the molecule.
            mols[i] = updateMolecule(system[molnums[i]].molecule(), used_atoms, map);
        }
    }

    // Update the system.
    system.update(Molecules(mols));

	// Update the System fileformat property to record that it includes
    // data from this file format.
    QString fileformat = this->formatName();

    PropertyName fileformat_property = map["fileformat"];

    try
    {
        QString last_format = system.property(fileformat_property)
                                    .asA<StringProperty>().value();

        fileformat = QString("%1,%2").arg(last_format,fileformat);
    }
    catch(...)
    {}

    if (fileformat_property.hasSource())
    {
        system.setProperty(fileformat_property.source(), StringProperty(fileformat));
    }
    else
    {
        system.setProperty("fileformat", StringProperty(fileformat));
    }
}

/** Internal function used to get the molecule structure for molecule 'imol'. */
MolStructureEditor PDB2::getMolStructure(int imol, const PropertyName &cutting) const
{
    // Make sure the molecule index is within range.
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

    // Now add any chains to the molecule.
    for (auto chain : chains[imol].uniqueKeys())
    {
        mol.add(ChainName(chain));
    }

    // Loop over all unique residues in the molecule.
    for (auto res_idx : residues[imol].uniqueKeys())
    {
        // By default we will use one CutGroup per residue.
        // This may be changed later by the cutting system.
        auto cutgroup = mol.add(CGName(QString::number(res_idx)));

        // Get a sorted list of the atoms that are part of the residue.
        QList<qint64> res_atoms = residues[imol].values(res_idx);
        qSort(res_atoms);

        // Use the first atom to get the name and number of this residue
        QString res_name = atoms[imol][res_atoms[0]].getResName();
        qint64  res_num  = atoms[imol][res_atoms[0]].getResNum();

        // Store the chain identifier and insertion code.
        QChar chain_id = atoms[imol][res_atoms[0]].getChainID();

        // Add the residue to the molecule.
        auto res = mol.add(ResNum(res_num));
        res.rename(ResName(res_name.trimmed()));

        // Reparent the residue to its chain.
        if (not chain_id.isSpace())
            res.reparent(ChainName(chain_id));

        // Add each atom in the residue to the molecule.
        for (auto res_atom : res_atoms)
        {
            auto atom = cutgroup.add(AtomNum(atoms[imol][res_atom].getNumber()));
            atom.rename(AtomName(atoms[imol][res_atom].getName().trimmed()));

            // Reparent the atom to its residue.
            atom.reparent(ResIdx(res_idx));
        }
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
    // Make sure the molecule (model) index is within range.
    if ((imol < 0) or (imol > atoms.count()))
    {
        throw SireError::program_bug(QObject::tr("The molecule index %1 is out of "
            "range, 0 - %2").arg(imol).arg(atoms.count()), CODELOC);
    }

    // Make sure that there are atoms in the molecule.
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

    // Atom property objects.
    AtomCoords         coords(molinfo);
    AtomCharges        charges(molinfo);
    AtomElements       elements(molinfo);
    AtomFloatProperty  occupancies(molinfo);
    AtomFloatProperty  temperatures(molinfo);
    AtomStringProperty is_het_atom(molinfo);

    // Residue property objects.
    ResStringProperty  insert_codes(molinfo);

    // The list of residues that have had properties updated.
    QVector<qint64> res_list;

    // Now loop through the atoms in the molecule and set each property.
    for (int i=0; i<nAtoms(imol); ++i)
    {
        // Store a reference to the current atom.
        const PDBAtom &atom = atoms[imol][i];

        // Determine the CGAtomIdx for this atom.
        auto cgatomidx = molinfo.cgAtomIdx(AtomNum(atom.getNumber()));

        // Set the properties.
        coords.set(cgatomidx, atom.getCoord());
        charges.set(cgatomidx, int(atom.getCharge()) * SireUnits::mod_electron);
        elements.set(cgatomidx, atom.getElement());
        occupancies.set(cgatomidx, atom.getOccupancy());
        temperatures.set(cgatomidx, atom.getTemperature());

        bool isHet = atom.isHet();

        if (isHet) is_het_atom.set(cgatomidx, "True");
        else       is_het_atom.set(cgatomidx, "False");

        // Store the residue index for the atom.
        qint64 res_idx = atom.getResIdx();

        // This residue hasn't already been processed.
        if (not res_list.contains(res_idx))
        {
            QChar icode = atom.getInsertCode();

            // Set the insert code property for this residue.
            if (not icode.isSpace())
            {
                insert_codes.set(ResIdx(res_idx), QString(icode));
            }

            res_list.append(res_idx);
        }
    }

    return mol.setProperty(map["coordinates"], coords)
              .setProperty(map["formal_charge"], charges)
              .setProperty(map["element"], elements)
              .setProperty(map["occupancy"], occupancies)
              .setProperty(map["beta_factor"], temperatures)
              .setProperty(map["is_het"], is_het_atom)
              .setProperty(map["insert_code"], insert_codes)
              .commit();
}

/** Internal function used to parse a Sire molecule view into a PDB ATOM records using
    the parameters in the property map. */
void PDB2::parseMolecule(const SireMol::Molecule &sire_mol, QVector<QString> &atom_lines,
    QStringList &errors, const SireBase::PropertyMap &map)
{
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
    QVector<bool> is_ter(num_atoms, false);

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

        tbb::parallel_for(tbb::blocked_range<int>(0, num_atoms),
                        [&](const tbb::blocked_range<int> &r)
        {
            // Create local data objects.
            QStringList local_errors;

            // Convert each atom into a PDBAtom object
            // and generate a PDB data record.
            for (int i=r.begin(); i<r.end(); ++i)
            {
                local_atoms[i] = PDBAtom(sire_mol.atom(AtomIdx(i)), is_ter[i], map, local_errors);
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
            PDBAtom atom(sire_mol.atom(AtomIdx(i)), is_ter[i], map, errors);

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

/** Internal function used to parse a add PDB coordinate data to an existing
    Sire Molecule. */
SireMol::Molecule PDB2::updateMolecule(const SireMol::Molecule &sire_mol,
    QVector<bool> &used_atoms, const SireBase::PropertyMap &map) const
{
    /* Now we need to match the atoms in the passed system to those
       in the PDB2 object.

       Annoyingly, the number and layout of molecules may not match,
       since PDB files that are used as a companion coordinate file
       for PSF often have useful molecule formatting information
       stripped, i.e. they are missing TER records.

       This means we need to matching atom numbers could be in any
       of the PDB molecules. As such, we loop over all PDB molecules
       to find matches. If we find a matching atom number we then
       check to see if the data matches that in the system. If not,
       we continue our search. If we find no valid matches then
       we throw an error and abort.
     */

    // Get the info object that can map between AtomNum to AtomIdx etc.
    const auto molinfo = sire_mol.info();

    // Initialise the property objects.
    AtomCoords coords(molinfo);
    AtomCharges        charges(molinfo);
    AtomElements       elements(molinfo);
    AtomFloatProperty  occupancies(molinfo);
    AtomFloatProperty  temperatures(molinfo);
    AtomStringProperty is_het_atom(molinfo);

    // Residue property objects.
    ResStringProperty  insert_codes(molinfo);

    // The list of residues that have had properties updated.
    QVector<qint64> res_list;

    // Internal function to update atom properties for AtomIdx 'i'
    auto update_atom = [&](int i, int mol_idx, int atom_idx)
    {
        // Determine the CGAtomIdx for this atom.
        auto cgatomidx = molinfo.cgAtomIdx(AtomIdx(i));

        // Store a reference to the atom.
        const PDBAtom &atom = atoms[mol_idx][atom_idx];

        // Set the properties.
        coords.set(cgatomidx, atom.getCoord());
        charges.set(cgatomidx, int(atom.getCharge()) * SireUnits::mod_electron);
        elements.set(cgatomidx, atom.getElement());
        occupancies.set(cgatomidx, atom.getOccupancy());
        temperatures.set(cgatomidx, atom.getTemperature());

        bool isHet = atom.isHet();

        if (isHet) is_het_atom.set(cgatomidx, "True");
        else       is_het_atom.set(cgatomidx, "False");

        // Store the residue index for the atom.
        qint64 res_idx = atom.getResIdx();

        // This residue hasn't already been processed.
        if (not res_list.contains(res_idx))
        {
            QChar icode = atom.getInsertCode();

            // Set the insert code property for this residue.
            if (not icode.isSpace())
            {
                insert_codes.set(ResIdx(res_idx), QString(icode));
            }

            res_list.append(res_idx);
        }
    };

    // Loop over all atoms in the molecule.
    // Here we are looping by AtomIdx.
    if (usesParallel())
    {
        QMutex mutex;

        tbb::parallel_for(tbb::blocked_range<int>(0, sire_mol.nAtoms()),
                           [&](const tbb::blocked_range<int> r)
        {
            for (int i=r.begin(); i<r.end(); ++i)
            {
                // The atom and molecule indices.
                int mol_idx, atom_idx;

                // Try to find a matching atom in the atoms vector.
                findAtom(sire_mol.atom(AtomIdx(i)), mol_idx, atom_idx, used_atoms);

                // Acquire a lock.
                QMutexLocker lkr(&mutex);

                // Update the atom properties.
                update_atom(i, mol_idx, atom_idx);
            }
        });
    }
    else
    {
        for (int i=0; i<sire_mol.nAtoms(); ++i)
        {
            // The atom and molecule indices.
            int mol_idx, atom_idx;

            // Try to find a matching atom in the atoms vector.
            findAtom(sire_mol.atom(AtomIdx(i)), mol_idx, atom_idx, used_atoms);

            // Update the atom properties.
            update_atom(i, mol_idx, atom_idx);
        }
    }

    // Set the additional properties.

    MolEditor edit_mol = sire_mol.edit();

    // Coordinates.
    if (not sire_mol.hasProperty(map["coordinates"]))
    {
        edit_mol.setProperty(map["coordinates"], coords);
    }
    else
    {
        if (not sire_mol.hasProperty(map["PDB.coordinates"]))
        {
            edit_mol.setProperty(map["PDB.coordinates"], coords);
        }
        else
        {
            edit_mol.setProperty(map["PDB.coordinates[2]"], coords);
        }
    }

    // Charges.
    if (not sire_mol.hasProperty(map["coordinates"]))
    {
        edit_mol.setProperty(map["formal_charge"], charges);
    }
    else
    {
        if (not sire_mol.hasProperty(map["PDB.formal_charge"]))
        {
            edit_mol.setProperty(map["PDB.formal_charge"], charges);
        }
        else
        {
            edit_mol.setProperty(map["PDB.formal_charge[2]"], charges);
        }
    }

    // Elements.
    if (not sire_mol.hasProperty(map["element"]))
    {
        edit_mol.setProperty(map["element"], elements);
    }
    else
    {
        if (not sire_mol.hasProperty(map["PDB.element"]))
        {
            edit_mol.setProperty(map["PDB.element"], elements);
        }
        else
        {
            edit_mol.setProperty(map["PDB.element[2]"], elements);
        }
    }

    // Occupancy.
    if (not sire_mol.hasProperty(map["occupancy"]))
    {
        edit_mol.setProperty(map["occupancy"], occupancies);
    }
    else
    {
        if (not sire_mol.hasProperty(map["PDB.occupancy"]))
        {
            edit_mol.setProperty(map["PDB.occupancy"], occupancies);
        }
        else
        {
            edit_mol.setProperty(map["PDB.occupancy[2]"], occupancies);
        }
    }

    // Temperature factors.
    if (not sire_mol.hasProperty(map["beta_factor"]))
    {
        edit_mol.setProperty(map["beta_factor"], temperatures);
    }
    else
    {
        if (not sire_mol.hasProperty(map["PDB.beta_factor"]))
        {
            edit_mol.setProperty(map["PDB.beta_factor"], temperatures);
        }
        else
        {
            edit_mol.setProperty(map["PDB.beta_factor[2]"], temperatures);
        }
    }

    // Hetero-atoms.
    if (not sire_mol.hasProperty(map["is_het"]))
    {
        edit_mol.setProperty(map["is_het"], is_het_atom);
    }
    else
    {
        if (not sire_mol.hasProperty(map["PDB.is_het"]))
        {
            edit_mol.setProperty(map["PDB.is_het"], is_het_atom);
        }
        else
        {
            edit_mol.setProperty(map["PDB.is_het[2]"], is_het_atom);
        }
    }

    // Residue insert codes.
    if (not sire_mol.hasProperty(map["insert_code"]))
    {
        edit_mol.setProperty(map["insert_code"], insert_codes);
    }
    else
    {
        if (not sire_mol.hasProperty(map["PDB.insert_code"]))
        {
            edit_mol.setProperty(map["PDB.insert_code"], insert_codes);
        }
        else
        {
            edit_mol.setProperty(map["PDB.insert_code[2]"], insert_codes);
        }
    }

    return edit_mol.commit();
}

/** Find the atom at index [ mol_idx ][ atom_idx ] in the atoms vector that
    matches the Sire Atom. Return false if no match is found. */
void PDB2::findAtom(const SireMol::Atom &sire_atom, int &mol_idx,
    int &atom_idx, QVector<bool> &used_atoms) const
{
    // Store the atom properties.
    const QString atom_name = sire_atom.name().value();
    const int atom_num = sire_atom.number().value();

    bool is_within_res = false;
    if (sire_atom.isWithinResidue())
        is_within_res = true;

    // Store residue properties.
    const QString res_name = is_within_res ?
        sire_atom.residue().name().value().left(3) : QString();

    const int res_num = is_within_res ?
        sire_atom.residue().number().value() : 0;

    // First, we'll locate the molecule that likely contains the
    // atom, i.e. the molecule where the last atom has a number
    // larger than the one that we're searching for.
    // This avoids doing an order N^2 search if we assume that
    // the atoms are likely in the same order in the lead file.

    // Molecule index.
    int imol = 0;
    while (atoms[imol][nAtoms(imol)-1].getNumber() < atom_num)
    {
        imol++;

        // Make sure the molecule index is in range.
        if (imol == nMolecules())
        {
            imol = 0;
            break;
        }
    }

    // Now search the molecule in question.

    // The index to start our search, assuming ascending atom numbers.
    int hint = atom_num - atoms[imol][0].getNumber();

    // Make sure the start index is in range.
    if (hint < 0 or hint >= nAtoms(imol))
        hint = 0;

    // A map of partial matches: score -> (mol_idx, atom_idx)
    QMap<int, QPair<int, int>> partial_matches;

    // Internal function to check whether we find a matching atom.
    auto is_match = [&](int i, int j)
    {
        bool same_atom_name = (atom_name == atoms[i][j].getName());
        bool same_atom_num  = (atom_num  == atoms[i][j].getNumber());

        // This is true by default, since some atoms aren't in residues.
        bool same_res_name = true;
        bool same_res_num  = true;

        if (is_within_res)
        {
            same_res_name = (res_name == atoms[i][j].getResName());
            same_res_num  = (res_num  == atoms[i][j].getResNum());
        }

        if (same_atom_name and same_res_name and j == hint)
        {
            return true;
        }
        if (same_atom_name and same_atom_num and
            same_res_name  and same_res_num)
        {
            return true;
        }
        else if (same_atom_name)
        {
            // Score a partial match. MUST match atom name, then ideally
            // residue name, atom number, residue number.
            partial_matches.insert((10 * same_res_name ) +
                                   ( 5 * same_atom_num ) +
                                   ( 2 * same_res_num  ) +
                                   ( 1 * (j == hint)   ),
                                   QPair<int, int>(i, j));
        }

        return false;
    };

    // Internal function to determine whether an atom has already been
    // used as a match for a different record.
    auto is_available = [&](int i, int j)
    {
        // Work out the index of this atom.

        int atom_idx = 0;

        for (int idx=0; idx<i; idx++)
            atom_idx += nAtoms(idx);

        atom_idx += j;

        if (used_atoms[atom_idx])
        {
            return false;
        }
        else
        {
            used_atoms[atom_idx] = true;
            return true;
        }
    };

    // Whether a matching atom was found.
    bool match_found = false;

    // Search all atoms in the suspected molecule.
    for (int i=hint; i<nAtoms(imol); i++)
    {
        if (is_match(imol, i) and
            is_available(imol, i))
        {
            mol_idx  = imol;
            atom_idx = i;

            match_found = true;
            break;
        }
    }

    // Didn't find a match, search all molecules.
    if (not match_found)
    {
        // Loop over all of the molecules.
        for (int i=0; i<nMolecules(); ++i)
        {
            // Loop over all of the atoms in the molecule.
            for (int j=0; j<nAtoms(i); ++j)
            {
                if (is_match(i, j) and
                    is_available(i, j))
                {
                    mol_idx  = i;
                    atom_idx = j;

                    match_found = true;
                    break;
                }
            }
        }
    }

    // Still no match, check for partial matches.
    if (not match_found)
    {
        // Try the partial matches in order of decreasing score.
        while (not partial_matches.isEmpty())
        {
            // Take the match with the highest score from the map.
            auto match = partial_matches.take(partial_matches.lastKey());

            mol_idx  = match.first;
            atom_idx = match.second;

            if (is_available(mol_idx, atom_idx))
            {
                match_found = true;
                break;
            }
        }

        if (not match_found)
        {
            if (is_within_res)
            {
                throw SireMol::missing_atom(QObject::tr("Could not find "
                    "a matching atom record: AtomName(\'%1\'), AtomNum(%2), ResName(\'%3\'), ResNum(%4)")
                    .arg(atom_name).arg(atom_num).arg(res_name).arg(res_num), CODELOC);
            }
            else
            {
                throw SireMol::missing_atom(QObject::tr("Could not find "
                    "a matching atom record: AtomName(\'%1\'), AtomNum(%2)")
                    .arg(atom_name).arg(atom_num), CODELOC);
            }
        }
    }
}
