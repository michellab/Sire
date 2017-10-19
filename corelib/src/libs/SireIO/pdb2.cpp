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
static const RegisterMetaType<PDBCrystal> r_pdbcrystal(NO_ROOT);
static const RegisterMetaType<PDBHelix> r_pdbhelix(NO_ROOT);
static const RegisterMetaType<PDBMaster> r_pdbmaster(NO_ROOT);
static const RegisterMetaType<PDBSheet> r_pdbsheet(NO_ROOT);
static const RegisterMetaType<PDBTitle> r_pdbtitle(NO_ROOT);
static const RegisterMetaType<PDBTransform> r_pdbtransform(NO_ROOT);

QDataStream SIREIO_EXPORT &operator<<(QDataStream &ds, const PDBAtom &pdbatom)
{
    writeHeader(ds, r_pdbatom, 1);

    SharedDataStream sds(ds);

    sds << pdbatom.record << pdbatom.serial << pdbatom.name << pdbatom.alt_loc
        << pdbatom.res_name << pdbatom.chain_id << pdbatom.res_num << pdbatom.insert_code
        << pdbatom.coord << pdbatom.occupancy << pdbatom.temperature << pdbatom.element
        << pdbatom.charge << pdbatom.is_het << pdbatom.is_ter << pdbatom.is_anis
        << pdbatom.anis_facts[0] << pdbatom.anis_facts[1] << pdbatom.anis_facts[2]
        << pdbatom.anis_facts[3] << pdbatom.anis_facts[4] << pdbatom.anis_facts[5];

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
            >> pdbatom.charge >> pdbatom.is_het >> pdbatom.is_ter >> pdbatom.is_anis
            >> pdbatom.anis_facts[0] >> pdbatom.anis_facts[1] >> pdbatom.anis_facts[2]
            >> pdbatom.anis_facts[3] >> pdbatom.anis_facts[4] >> pdbatom.anis_facts[5];
    }
    else
        throw version_error(v, "1", r_pdbatom, CODELOC);

    return ds;
}

QDataStream SIREIO_EXPORT &operator<<(QDataStream &ds, const PDBCrystal &pdbcrystal)
{
    writeHeader(ds, r_pdbcrystal, 1);

    SharedDataStream sds(ds);

    sds << pdbcrystal.record << pdbcrystal.base_lengths
        << pdbcrystal.angles << pdbcrystal.z;

    return ds;
}

QDataStream SIREIO_EXPORT &operator>>(QDataStream &ds, PDBCrystal &pdbcrystal)
{
    VersionID v = readHeader(ds, r_pdbcrystal);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> pdbcrystal.record >> pdbcrystal.base_lengths
            >> pdbcrystal.angles >> pdbcrystal.z;
    }
    else
        throw version_error(v, "1", r_pdbcrystal, CODELOC);

    return ds;
}

QDataStream SIREIO_EXPORT &operator<<(QDataStream &ds, const PDBHelix &pdbhelix)
{
    writeHeader(ds, r_pdbhelix, 1);

    SharedDataStream sds(ds);

    sds << pdbhelix.record << pdbhelix.serial << pdbhelix.id << pdbhelix.init_res_name
        << pdbhelix.init_chain_id << pdbhelix.init_res_num << pdbhelix.init_insert_code
        << pdbhelix.end_res_name << pdbhelix.end_chain_id << pdbhelix.end_insert_code
        << pdbhelix.helix_class << pdbhelix.comment << pdbhelix.length;

    return ds;
}

QDataStream SIREIO_EXPORT &operator>>(QDataStream &ds, PDBHelix &pdbhelix)
{
    VersionID v = readHeader(ds, r_pdbhelix);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> pdbhelix.record >> pdbhelix.serial >> pdbhelix.id >> pdbhelix.init_res_name
            >> pdbhelix.init_chain_id >> pdbhelix.init_res_num >> pdbhelix.init_insert_code
            >> pdbhelix.end_res_name >> pdbhelix.end_chain_id >> pdbhelix.end_insert_code
            >> pdbhelix.helix_class >> pdbhelix.comment << pdbhelix.length;
    }
    else
        throw version_error(v, "1", r_pdbhelix, CODELOC);

    return ds;
}

QDataStream SIREIO_EXPORT &operator<<(QDataStream &ds, const PDBMaster &pdbmaster)
{
    writeHeader(ds, r_pdbmaster, 1);

    SharedDataStream sds(ds);

    sds << pdbmaster.record << pdbmaster.num_remarks << pdbmaster.num_hets
        << pdbmaster.num_helices << pdbmaster.num_sheets << pdbmaster.num_sites
        << pdbmaster.num_transforms << pdbmaster.num_atoms << pdbmaster.num_ters
        << pdbmaster.num_connects << pdbmaster.num_sequences;

    return ds;
}

QDataStream SIREIO_EXPORT &operator>>(QDataStream &ds, PDBMaster &pdbmaster)
{
    VersionID v = readHeader(ds, r_pdbmaster);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> pdbmaster.record >> pdbmaster.num_remarks >> pdbmaster.num_hets
            >> pdbmaster.num_helices >> pdbmaster.num_sheets >> pdbmaster.num_sites
            >> pdbmaster.num_transforms >> pdbmaster.num_atoms >> pdbmaster.num_ters
            >> pdbmaster.num_connects >> pdbmaster.num_sequences;
    }
    else
        throw version_error(v, "1", r_pdbmaster, CODELOC);

    return ds;
}

QDataStream SIREIO_EXPORT &operator<<(QDataStream &ds, const PDBSheet &pdbsheet)
{
    writeHeader(ds, r_pdbsheet, 1);

    SharedDataStream sds(ds);

    sds << pdbsheet.record << pdbsheet.id << pdbsheet.num_strands
        << pdbsheet.init_res_name << pdbsheet.init_chain_id << pdbsheet.init_res_num
        << pdbsheet.init_insert_code << pdbsheet.end_res_name << pdbsheet.end_chain_id
        << pdbsheet.end_res_num << pdbsheet.end_insert_code << pdbsheet.sense
        << pdbsheet.curr_atm_name << pdbsheet.curr_res_name << pdbsheet.curr_chain_id
        << pdbsheet.curr_insert_code << pdbsheet.prev_atm_name << pdbsheet.prev_res_name
        << pdbsheet.prev_chain_id << pdbsheet.prev_insert_code;

    return ds;
}

QDataStream SIREIO_EXPORT &operator>>(QDataStream &ds, PDBSheet &pdbsheet)
{
    VersionID v = readHeader(ds, r_pdbsheet);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> pdbsheet.record >> pdbsheet.id >> pdbsheet.num_strands
            >> pdbsheet.init_res_name >> pdbsheet.init_chain_id >> pdbsheet.init_res_num
            >> pdbsheet.init_insert_code >> pdbsheet.end_res_name >> pdbsheet.end_chain_id
            >> pdbsheet.end_res_num >> pdbsheet.end_insert_code >> pdbsheet.sense
            >> pdbsheet.curr_atm_name >> pdbsheet.curr_res_name >> pdbsheet.curr_chain_id
            >> pdbsheet.curr_insert_code >> pdbsheet.prev_atm_name >> pdbsheet.prev_res_name
            >> pdbsheet.prev_chain_id >> pdbsheet.prev_insert_code;
    }
    else
        throw version_error(v, "1", r_pdbsheet, CODELOC);

    return ds;
}

QDataStream SIREIO_EXPORT &operator<<(QDataStream &ds, const PDBTitle &pdbtitle)
{
    writeHeader(ds, r_pdbtitle, 1);

    SharedDataStream sds(ds);

    sds << pdbtitle.records << pdbtitle.header << pdbtitle.obsoletes << pdbtitle.titles
        << pdbtitle.splits << pdbtitle.caveats << pdbtitle.compounds << pdbtitle.sources
        << pdbtitle.keywords << pdbtitle.experiments << pdbtitle.num_models
        << pdbtitle.model_types << pdbtitle.authors << pdbtitle.revisions
        << pdbtitle.supersedes << pdbtitle.journals << pdbtitle.remarks;

    return ds;
}

QDataStream SIREIO_EXPORT &operator>>(QDataStream &ds, PDBTitle &pdbtitle)
{
    VersionID v = readHeader(ds, r_pdbtitle);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> pdbtitle.records >> pdbtitle.header >> pdbtitle.obsoletes >> pdbtitle.titles
            >> pdbtitle.splits >> pdbtitle.caveats >> pdbtitle.compounds >> pdbtitle.sources
            >> pdbtitle.keywords >> pdbtitle.experiments >> pdbtitle.num_models
            >> pdbtitle.model_types >> pdbtitle.authors >> pdbtitle.revisions
            >> pdbtitle.supersedes >> pdbtitle.journals >> pdbtitle.remarks;
    }
    else
        throw version_error(v, "1", r_pdbtitle, CODELOC);

    return ds;
}

QDataStream SIREIO_EXPORT &operator<<(QDataStream &ds, const PDBTransform &pdbtransform)
{
    writeHeader(ds, r_pdbtransform, 1);

    SharedDataStream sds(ds);

    sds << pdbtransform.records << pdbtransform.serial << pdbtransform.isGiven
        << pdbtransform.transforms << pdbtransform.offsets;

    return ds;
}

QDataStream SIREIO_EXPORT &operator>>(QDataStream &ds, PDBTransform &pdbtransform)
{
    VersionID v = readHeader(ds, r_pdbtransform);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> pdbtransform.records >> pdbtransform.serial >> pdbtransform.isGiven
            >> pdbtransform.transforms >> pdbtransform.offsets;
    }
    else
        throw version_error(v, "1", r_pdbcrystal, CODELOC);

    return ds;
}

QDataStream SIREIO_EXPORT &operator<<(QDataStream &ds, const PDB2 &pdb2)
{
    writeHeader(ds, r_pdb2, 1);

    SharedDataStream sds(ds);

    sds << pdb2.title << pdb2.atoms << pdb2.residues << pdb2.chains << pdb2.segments
        << pdb2.connections << pdb2.helices << pdb2.sheets << pdb2.trans_orig
        << pdb2.trans_scale << pdb2.trans_matrix << pdb2.master << pdb2.num_ters
        << pdb2.invalid_records << pdb2.parse_warnings
        << static_cast<const MoleculeParser&>(pdb2);

    return ds;
}

QDataStream SIREIO_EXPORT &operator>>(QDataStream &ds, PDB2 &pdb2)
{
    VersionID v = readHeader(ds, r_pdb2);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> pdb2.title >> pdb2.atoms >> pdb2.residues >> pdb2.chains >> pdb2.segments
            >> pdb2.connections >> pdb2.helices >> pdb2.sheets >> pdb2.trans_orig
            >> pdb2.trans_scale >> pdb2.trans_matrix >> pdb2.master >> pdb2.num_ters
            >> pdb2.invalid_records >> pdb2.parse_warnings
            >> static_cast<MoleculeParser&>(pdb2);
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
    is_ter(false),
    is_anis(false)
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
    is_ter(false),
    is_anis(false)
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

    @param errors
        An array of error messages.
 */
PDBAtom::PDBAtom(const SireMol::Atom &atom, QStringList &errors) :
    serial(atom.number().value()),
    name(atom.name().value()),
    occupancy(1.0),
    temperature(0.0),
    element("X"),
    charge(0),
    is_het(false),
    is_ter(false),
    is_anis(false)
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

        // Now work out whether this is a terminal atom.

        // Extract the number of the last atom in the chain
        int terminal_atom = atom.chain()
                           .atoms()[atom.chain().atoms().count()-1]
                           .read()
                           .asA<SireMol::Atom>()
                           .number()
                           .value();

        // This is the last atom in the chain.
        if (terminal_atom == serial) is_ter = true;
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
    else if (atom.hasProperty("charge"))
    {
        // TODO: Some kind of conversion needed?
        charge = atom.property<SireUnits::Dimension::Charge>("charge").value();
    }
    else
    {
        // TODO: Infer the charge...
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

/** Set anisotropic temperature record data.
    @param line1
        An ATOM record line from a PDB file.

    @param line2
        The ANISOU record line for the atom.

    @param errors
        An array of error messages.
 */
void PDBAtom::setAnisTemp(const QString &line1, const QString &line2, QStringList &errors)
{
    is_anis = true;

    // Check that data from this record matches that of the atom.
    // Columns 6-26 and 72-80 should be identical.
    if (line1.midRef(6,21) == line2.midRef(6,21) and
        line1.midRef(72,8) == line2.midRef(72,8))
    {
        // Extract the six anisotropic temperature factors.
        for (int i=0; i<6; ++i)
        {
            bool ok;
            int tmp_int = line2.midRef(28 + i*7,7).toInt(&ok);

            if (not ok)
            {
                errors.append(QObject::tr("There was a problem reading the "
                    "anisotropic temperature factors from the data '%1' in "
                    "line '%2'").arg(line2.mid(28,42)).arg(line2));

                return;
            }

            // Store the temperature factor.
            anis_facts[i] = tmp_int;
        }
    }
    else
    {
        // Do we wan't to bail out here?
        // Probably just catch the error and continue.
        errors.append(QObject::tr("ANISOU record does not match the format! '%2' '%3'")
            .arg(line1).arg(line2));
    }
}

/** Get the atom serial number. */
qint64 PDBAtom::getSerial() const
{
    return serial;
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

/** Default constructor. */
PDBCrystal::PDBCrystal() : record("NULL")
{
}

/** Constructor.
    @param line
        An CRYST1 record line from a PDB file.

    @param errors
        An array of error messages.
 */
PDBCrystal::PDBCrystal(const QString &line,
                       QStringList &errors) : record(line)
{
    // Read unit cell base length records.

    bool ok_x, ok_y, ok_z;
    double x = line.midRef(6,9).toDouble(&ok_x);
    double y = line.midRef(15,9).toDouble(&ok_y);
    double z = line.midRef(24,9).toDouble(&ok_z);

    if (not (ok_x and ok_y and ok_z))
    {
        errors.append(QObject::tr("There was a problem reading the crystallographic "
            "values of a, c, and c from the data '%1' in line '%2'")
            .arg(line.mid(6,27)).arg(line));

        return;
    }
    base_lengths = Vector(x, y, z);

    // Read unit cell angle records.

    x = line.midRef(33,7).toDouble(&ok_x);
    y = line.midRef(40,7).toDouble(&ok_y);
    z = line.midRef(47,7).toDouble(&ok_z);

    if (not (ok_x and ok_y and ok_z))
    {
        errors.append(QObject::tr("There was a problem reading the crystallographic "
            "values of alpha, beta, and gamma from the data '%1' in line '%2'")
            .arg(line.mid(33,21)).arg(line));

        return;
    }
    angles = Vector(x, y, z);

    // Extract the crystallographic space group.
    space_group = line.mid(55,11).simplified();

    // Extract the Z value.
    int tmp_int = line.midRef(66,4).toInt(&ok_x);

    if (not ok_x)
    {
        errors.append(QObject::tr("There was a problem reading the crystallographic "
            "Z value from the data '%1' in line '%2'")
            .arg(line.mid(66,4)).arg(line));

        return;
    }
    z = tmp_int;
}

/** Return the C++ name for this class */
const char* PDBCrystal::typeName()
{
    return QMetaType::typeName( qMetaTypeId<PDBCrystal>() );
}

/** Return a PDB record line for this crystallographic object. */
QString PDBCrystal::toPDBRecord() const
{
    return record;
}

/** Return a string representation of this object */
QString PDBCrystal::toString() const
{
    return QObject::tr("PDBCrystal::null");
}

/** Whether the object contains a record. */
bool PDBCrystal::hasRecord() const
{
    return (record != "NULL");
}

/** Default constructor. */
PDBHelix::PDBHelix()
{
}

/** Constructor.
    @param line
        A HELIX record line from a PDB file.

    @param errors
        An array of error messages.
 */
PDBHelix::PDBHelix(const QString &line,
                   QStringList &errors) : record(line)
{
    // Extract the helix serial number.
    bool ok;
    int tmp_int = line.midRef(7,3).toInt(&ok);

    if (not ok)
    {
        errors.append(QObject::tr("Cannot extract the helix serial "
            "number from part (%1) from line '%2'")
            .arg(line.mid(7,3)).arg(line));

        return;
    }
    serial = tmp_int;

    // Extract the helix ID.
    id = line.mid(11,3).simplified();

    // Extract the name of the initial residue.
    init_res_name = line.mid(15,3).simplified();

    // Extract the ID of the initial chain.
    init_chain_id = line[19];

    // Extract the sequence number of the initial residue.
    tmp_int = line.midRef(21,4).toInt(&ok);

    if (not ok)
    {
        errors.append(QObject::tr("Cannot extract the initial residue sequence "
            "number from part (%1) from line '%2'")
            .arg(line.mid(21,4)).arg(line));

        return;
    }
    init_res_num = tmp_int;

    // Extract the insertion code of the end residue.
    end_insert_code = line[25];

    // Extract the name of the end residue.
    end_res_name = line.mid(27,3).simplified();

    // Extract the ID of the end chain.
    end_chain_id = line[31];

    // Extract the sequence number of the end residue.
    tmp_int = line.midRef(33,4).toInt(&ok);

    if (not ok)
    {
        errors.append(QObject::tr("Cannot extract the end residue sequence "
            "number from part (%1) from line '%2'")
            .arg(line.mid(33,4)).arg(line));

        return;
    }
    end_res_num = tmp_int;

    // Extract the insertion code of the end residue.
    end_insert_code = line[37];

    // Extract the helix class.
    tmp_int = line.midRef(38,2).toInt(&ok);

    if (not ok)
    {
        errors.append(QObject::tr("Cannot extract the helix class "
            "from part (%1) from line '%2'")
            .arg(line.mid(38,2)).arg(line));

        return;
    }
    helix_class = tmp_int;

    // Extract any comment about the helix.
    comment = line.mid(40,30).simplified();

    // Extract the helix length.
    tmp_int = line.midRef(71,6).toInt(&ok);

    if (not ok)
    {
        errors.append(QObject::tr("Cannot extract the helix length "
            "from part (%1) from line '%2'")
            .arg(line.mid(71,6)).arg(line));

        return;
    }
    length = tmp_int;
}

/** Return the C++ name for this class */
const char* PDBHelix::typeName()
{
    return QMetaType::typeName( qMetaTypeId<PDBHelix>() );
}

/** Return a PDB record line for this helix. */
QString PDBHelix::toPDBRecord() const
{
    return record;
}

/** Return a string representation of this object */
QString PDBHelix::toString() const
{
    return QObject::tr("PDBHelix::null");
}

/** Default constructor. */
PDBMaster::PDBMaster() :
    num_remarks(0),
    num_hets(0),
    num_helices(0),
    num_sheets(0),
    num_sites(0),
    num_transforms(0),
    num_atoms(0),
    num_ters(0),
    num_connects(0),
    num_sequences(0)
{
}

/** Constructor.
    @param line
        A MASTER record line from a PDB file.

    @param errors
        An array of error messages.
 */
PDBMaster::PDBMaster(const QString &line,
                     QStringList &errors) :
    record(line),
    num_remarks(0),
    num_hets(0),
    num_helices(0),
    num_sheets(0),
    num_sites(0),
    num_transforms(0),
    num_atoms(0),
    num_ters(0),
    num_connects(0),
    num_sequences(0)
{
    // Extract the number of remark records.
    bool ok;
    int tmp_int = line.midRef(10,5).toInt(&ok);

    if (not ok)
    {
        errors.append(QObject::tr("Cannot extract the number of REMARK records "
            "from part (%1) from line '%2'")
            .arg(line.mid(10,5)).arg(line));

        return;
    }
    num_remarks = tmp_int;

    // Extract the number of HET records.
    tmp_int = line.midRef(20,5).toInt(&ok);

    if (not ok)
    {
        errors.append(QObject::tr("Cannot extract the number of HET records "
            "from part (%1) from line '%2'")
            .arg(line.mid(20,5)).arg(line));

        return;
    }
    num_hets = tmp_int;

    // Extract the number of HELIX records.
    tmp_int = line.midRef(25,5).toInt(&ok);

    if (not ok)
    {
        errors.append(QObject::tr("Cannot extract the number of HELIX records "
            "from part (%1) from line '%2'")
            .arg(line.mid(25,5)).arg(line));

        return;
    }
    num_helices = tmp_int;

    // Extract the number of SHEET records.
    tmp_int = line.midRef(30,5).toInt(&ok);

    if (not ok)
    {
        errors.append(QObject::tr("Cannot extract the number of SHEET records "
            "from part (%1) from line '%2'")
            .arg(line.mid(30,5)).arg(line));

        return;
    }
    num_sheets = tmp_int;

    // Extract the number of SITE records.
    tmp_int = line.midRef(40,5).toInt(&ok);

    if (not ok)
    {
        errors.append(QObject::tr("Cannot extract the number of SITE records "
            "from part (%1) from line '%2'")
            .arg(line.mid(40,5)).arg(line));

        return;
    }
    num_sites = tmp_int;

    // Extract the number of transformation records.
    tmp_int = line.midRef(45,5).toInt(&ok);

    if (not ok)
    {
        errors.append(QObject::tr("Cannot extract the number of coordinate "
            "transformation records from part (%1) from line '%2'")
            .arg(line.mid(45,5)).arg(line));

        return;
    }
    num_transforms = tmp_int;

    // Extract the number of ATOM and HETATM records.
    tmp_int = line.midRef(50,5).toInt(&ok);

    if (not ok)
    {
        errors.append(QObject::tr("Cannot extract the number of ATOM and HETATM records "
            "from part (%1) from line '%2'")
            .arg(line.mid(50,5)).arg(line));

        return;
    }
    num_atoms = tmp_int;

    // Extract the number of TER records.
    tmp_int = line.midRef(55,5).toInt(&ok);

    if (not ok)
    {
        errors.append(QObject::tr("Cannot extract the number of TER records "
            "from part (%1) from line '%2'")
            .arg(line.mid(55,5)).arg(line));

        return;
    }
    num_ters = tmp_int;

    // Extract the number of CONECT records.
    tmp_int = line.midRef(60,5).toInt(&ok);

    if (not ok)
    {
        errors.append(QObject::tr("Cannot extract the number of CONECT records "
            "from part (%1) from line '%2'")
            .arg(line.mid(60,5)).arg(line));

        return;
    }
    num_connects = tmp_int;

    // Extract the number of SEQRES records.
    tmp_int = line.midRef(65,5).toInt(&ok);

    if (not ok)
    {
        errors.append(QObject::tr("Cannot extract the number of SEQRES records "
            "from part (%1) from line '%2'")
            .arg(line.mid(65,5)).arg(line));

        return;
    }
    num_sequences = tmp_int;
}

/** Return the C++ name for this class */
const char* PDBMaster::typeName()
{
    return QMetaType::typeName( qMetaTypeId<PDBMaster>() );
}

/** Return a PDB record line for this master object. */
QString PDBMaster::toPDBRecord() const
{
    return record;
}

/** Return a string representation of this object */
QString PDBMaster::toString() const
{
    return QObject::tr("PDBMaster::null");
}

/** Return the number of REMARK records. */
int PDBMaster::nRemarks() const
{
    return num_remarks;
}

/** Return the number of HET records. */
int PDBMaster::nHets() const
{
    return num_hets;
}

/** Return the number of HELIX records. */
int PDBMaster::nHelices() const
{
    return num_helices;
}

/** Return the number of SHEET records. */
int PDBMaster::nSheets() const
{
    return num_sheets;
}

/** Return the number of SITE records. */
int PDBMaster::nSites() const
{
    return num_sites;
}

/** Return the number of coordinate transformation records. */
int PDBMaster::nTransforms() const
{
    return num_transforms;
}

/** Return the number of ATOM and HETATM records. */
int PDBMaster::nAtoms() const
{
    return num_atoms;
}

/** Return the number of TER records. */
int PDBMaster::nTers() const
{
    return num_ters;
}

/** Return the number of CONECT records. */
int PDBMaster::nConnects() const
{
    return num_connects;
}

/** Return the number of SEQRES records. */
int PDBMaster::nSequences() const
{
    return num_sequences;
}

/** Default constructor. */
PDBSheet::PDBSheet()
{
}

/** Constructor.
    @param line
        A SHEET record line from a PDB file.

    @param errors
        An array of error messages.
 */
PDBSheet::PDBSheet(const QString &line,
                   QStringList &errors) : record(line)
{
    // Extract the strand number.
    bool ok;
    int tmp_int = line.midRef(7,3).toInt(&ok);

    if (not ok)
    {
        errors.append(QObject::tr("Cannot extract the strand "
            "number from part (%1) from line '%2'")
            .arg(line.mid(7,3)).arg(line));

        return;
    }
    strand = tmp_int;

    // Extract the sheet ID.
    id = line.mid(11,3).simplified();

    // Extract the number of strands in the sheet.
    tmp_int = line.midRef(14,2).toInt(&ok);

    if (not ok)
    {
        errors.append(QObject::tr("Cannot extract the number of strands "
            "from part (%1) from line '%2'")
            .arg(line.mid(14,2)).arg(line));

        return;
    }
    num_strands = tmp_int;

    // Extract the name of the initial residue.
    init_res_name = line.mid(17,3).simplified();

    // Extract the ID of the initial chain.
    init_chain_id = line[21];

    // Extract the sequence number of the initial residue.
    tmp_int = line.midRef(22,4).toInt(&ok);

    if (not ok)
    {
        errors.append(QObject::tr("Cannot extract the initial residue sequence "
            "number from part (%1) from line '%2'")
            .arg(line.mid(22,4)).arg(line));

        return;
    }
    init_res_num = tmp_int;

    // Extract the insertion code of the end residue.
    end_insert_code = line[26];

    // Extract the name of the end residue.
    end_res_name = line.mid(28,3).simplified();

    // Extract the ID of the end chain.
    end_chain_id = line[32];

    // Extract the sequence number of the end residue.
    tmp_int = line.midRef(33,4).toInt(&ok);

    if (not ok)
    {
        errors.append(QObject::tr("Cannot extract the end residue sequence "
            "number from part (%1) from line '%2'")
            .arg(line.mid(33,4)).arg(line));

        return;
    }
    end_res_num = tmp_int;

    // Extract the insertion code of the end residue.
    end_insert_code = line[37];

    // Extract the sense of the strand with respect to previous strand.
    tmp_int = line.midRef(38,2).toInt(&ok);

    if (not ok)
    {
        errors.append(QObject::tr("Cannot extract the strand sense "
            "from part (%1) from line '%2'")
            .arg(line.mid(38,2)).arg(line));

        return;
    }
    sense = tmp_int;

    // Parse stand data.
    if (strand > 1)
    {
        // DATA FOR CURRENT STRAND

        // Extract the atom name in the current strand.
        curr_atm_name = line.mid(41,4).simplified();

        // Extract the residue name in the current strand.
        curr_res_name = line.mid(45,3).simplified();

        // Extract the chain ID in the current strand.
        curr_chain_id = line[49];

        // Extract the sequence number in the current strand.
        tmp_int = line.midRef(50,4).toInt(&ok);

        if (not ok)
        {
            errors.append(QObject::tr("Cannot extract the residue sequence "
                "number from part (%1) from line '%2'")
                .arg(line.mid(50,4)).arg(line));

            return;
        }
        curr_res_num = tmp_int;

        // Extract the insertion code in the current strand.
        curr_insert_code = line[54];

        // DATA FOR PREVIOUS STRAND

        // Extract the atom name in the previous strand.
        prev_atm_name = line.mid(56,4).simplified();

        // Extract the residue name in the previous strand.
        prev_res_name = line.mid(60,3).simplified();

        // Extract the chain ID in the previous strand.
        prev_chain_id = line[64];

        // Extract the sequence number in the previous strand.
        tmp_int = line.midRef(65,4).toInt(&ok);

        if (not ok)
        {
            errors.append(QObject::tr("Cannot extract the residue sequence "
                "number from part (%1) from line '%2'")
                .arg(line.mid(65,4)).arg(line));

            return;
        }
        prev_res_num = tmp_int;

        // Extract the insertion code in the previous strand.
        prev_insert_code = line[69];
    }
}

/** Return the C++ name for this class */
const char* PDBSheet::typeName()
{
    return QMetaType::typeName( qMetaTypeId<PDBSheet>() );
}

/** Return a PDB record line for this sheet. */
QString PDBSheet::toPDBRecord() const
{
    return record;
}

/** Return a string representation of this object */
QString PDBSheet::toString() const
{
    return QObject::tr("PDBSheet::null");
}

/** Default constructor. */
PDBTitle::PDBTitle()
{
}

/** Append a PDB record.
    @param line
        The PDB record.

    @param record_type
        The title record type.

    @param errors
        An array of error messages.
 */
void PDBTitle::appendRecord(const QString &line,
                            RECORD_TYPE record_type,
                            QStringList &errors)
{
    // Store the orginal record string.
    records.append(line);

    switch(record_type)
    {
        case RECORD_TYPE::HEADER: header = line;
                                  break;

        case RECORD_TYPE::OBSLTE: obsoletes.append(line);
                                  break;

        case RECORD_TYPE::TITLE : titles.append(line);
                                  break;

        case RECORD_TYPE::SPLIT : splits.append(line);
                                  break;

        case RECORD_TYPE::CAVEAT: caveats.append(line);
                                  break;

        case RECORD_TYPE::COMPND: compounds.append(line);
                                  break;

        case RECORD_TYPE::SOURCE: sources.append(line);
                                  break;

        case RECORD_TYPE::KEYWDS: keywords.append(line);
                                  break;

        case RECORD_TYPE::EXPDTA: experiments.append(line);
                                  break;

        case RECORD_TYPE::MDLTYP: model_types.append(line);
                                  break;

        case RECORD_TYPE::AUTHOR: authors.append(line);
                                  break;

        case RECORD_TYPE::REVDAT: revisions.append(line);
                                  break;

        case RECORD_TYPE::SPRSDE: supersedes.append(line);
                                  break;

        case RECORD_TYPE::JRNL  : journals.append(line);
                                  break;

        case RECORD_TYPE::REMARK: remarks.append(line);
                                  break;

        case RECORD_TYPE::NUMMDL:
        {
            bool ok;
            int tmp_int = line.midRef(10,4).toInt(&ok);

            if (not ok)
            {
                errors.append(QObject::tr("Cannot extract the number of models "
                    "from part (%1) from line '%2'")
                    .arg(line.mid(10,4)).arg(line));

                return;
            }
            num_models = tmp_int;

            break;
        }

        default:
        {
            errors.append(QObject::tr("Title section record not recognised! "
                "'%2'").arg(line));

            return;
        }
    }
}

/** Return PDB records for the title section. */
QStringList PDBTitle::toPDBRecord() const
{
    return records;
}

/** Return a string representation of this object */
QString PDBTitle::toString() const
{
    return QObject::tr("PDBTitle::null");
}

/** Return the number of title section records. */
int PDBTitle::nRecords() const
{
    return records.count();
}

/** Return the number of REMARK records. */
int PDBTitle::nRemarks() const
{
    return remarks.count();
}

/** Return the number of models that should be in the PDB file. */
int PDBTitle::nModels() const
{
    return num_models;
}

/** Whether the object contains any records. */
bool PDBTitle::hasRecords() const
{
    return (records.count() > 0);
}

/** Return the C++ name for this class */
const char* PDBTitle::typeName()
{
    return QMetaType::typeName( qMetaTypeId<PDBTitle>() );
}

/** Default constructor. */
PDBTransform::PDBTransform() :
    serial(-1),
    isGiven(false),
    transforms(QVector<SireMaths::Vector>(3)),
    isDimension{false, false, false}
{
}

/** Append a record to the transformation object.
    @param line
        An transformation record line from a PDB file.

    @param dimension
        Whether this corresponds to an x, y, or z record.

    @param isNonCryst
        Whether this is a non-crystallographic transformation.

    @param errors
        An array of error messages.
 */
void PDBTransform::appendRecord(const QString &line,
                                int dimension,
                                bool isNonCryst,
                                QStringList &errors)
{
    // Store the orginal record string.
    records.append(line);

    // Flag that we've seen a record for this dimension.
    isDimension[dimension] = true;

    // Extract transformation matrix entry for dimension.

    bool ok_x, ok_y, ok_z;
    double x = line.midRef(10,10).toDouble(&ok_x);
    double y = line.midRef(20,10).toDouble(&ok_y);
    double z = line.midRef(30,10).toDouble(&ok_z);

    if (not (ok_x and ok_y and ok_z))
    {
        errors.append(QObject::tr("There was a problem reading the transformation "
            "values of x, y, and z from the data '%1' in line '%2'")
            .arg(line.mid(10,30)).arg(line));

        return;
    }
    transforms[dimension] = Vector(x, y, z);

    // Extract offset.

    bool ok;
    double offset = line.midRef(45,10).toDouble(&ok);

    if (not ok)
    {
        errors.append(QObject::tr("There was a problem reading the transformation "
            "offset value from the data '%1' in line '%2'")
            .arg(line.mid(45,10)).arg(line));

        return;
    }
    if      (dimension == 0) offsets.setX(offset);
    else if (dimension == 1) offsets.setY(offset);
    else if (dimension == 2) offsets.setZ(offset);

    // Read additional MTRIX record data.

    if (isNonCryst)
    {
        // Extract serial number.

        bool ok;
        int tmp_int = line.midRef(7,3).toInt(&ok);

        if (not ok)
        {
            errors.append(QObject::tr("There was a problem reading the transformation "
                "serial number from the data '%1' in line '%2'")
                .arg(line.mid(7,3)).arg(line));

            return;
        }
        serial = tmp_int;

        // Check whether the coordinates for this entry are approximated
        // by this transformation.

        tmp_int = line[60].toLatin1();

        if (tmp_int) isNonCryst = true;
        else         isNonCryst = false;
    }
}

/** Return PDB records for the transformation record. */
QStringList PDBTransform::toPDBRecord() const
{
    return records;
}

/** Return a string representation of this object */
QString PDBTransform::toString() const
{
    return QObject::tr("PDBTransform::null");
}

/** Whether there is a complete transformation record. */
bool PDBTransform::hasRecord() const
{
    return (isDimension[0] and
            isDimension[1] and
            isDimension[2]);
}

/** Return the C++ name for this class */
const char* PDBTransform::typeName()
{
    return QMetaType::typeName( qMetaTypeId<PDBTransform>() );
}

/** Constructor */
PDB2::PDB2() : ConcreteProperty<PDB2,MoleculeParser>(),
               num_ters(0),
               has_master(false)
{
}

/** Construct to read in the data from the file called 'filename'. The
    passed property map can be used to pass extra parameters to control
    the parsing */
PDB2::PDB2(const QString &filename, const PropertyMap &map) :
    ConcreteProperty<PDB2,MoleculeParser>(filename,map),
    num_ters(0),
    has_master(false)
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
    ConcreteProperty<PDB2,MoleculeParser>(lines,map),
    num_ters(0),
    has_master(false)
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
    ConcreteProperty<PDB2,MoleculeParser>(map),
    num_ters(0),
    has_master(false)
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
        lines.append("MODEL");
        lines.append(QString("%1").arg(i+1));

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
    ConcreteProperty<PDB2,MoleculeParser>(other)
{}

/** Destructor */
PDB2::~PDB2()
{}

/** Copy assignment operator */
PDB2& PDB2::operator=(const PDB2 &other)
{
    if (this != &other)
    {
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

/** Return the number of title section records. */
int PDB2::nTitles() const
{
    return title.nRecords();
}

/** Return the number of models (molecules). */
int PDB2::nModels() const
{
    return atoms.count();
}

/** Return the number of atoms. */
int PDB2::nAtoms() const
{
    return atoms[0].count();
}

/** Return the number of helices. */
int PDB2::nHelices() const
{
    return helices.count();
}

/** Return the number of sheets. */
int PDB2::nSheets() const
{
    return sheets.count();
}

/** Return the number of TER records. */
int PDB2::nTers() const
{
    return num_ters;
}

/** Whether the object contains crystallographic record data. */
bool PDB2::hasCrystal() const
{
    return crystal.hasRecord();
}

/** Whether the object contains a MASTER record. */
bool PDB2::hasMaster() const
{
    return has_master;
}

/** Whether the object contains an ORIGX transformation record. */
bool PDB2::hasTransOrig() const
{
    return trans_orig.hasRecord();
}

/** Whether the object contains an SCALE transformation record. */
bool PDB2::hasTransScale() const
{
    return trans_scale.hasRecord();
}

/** Whether the object contains an MTRIX transformation record. */
bool PDB2::hasTransMatrix() const
{
    return trans_matrix.hasRecord();
}

/** Return the title object. */
PDBTitle PDB2::getTitle() const
{
    return title;
}

/** Return the atoms. */
QVector<QVector<PDBAtom> > PDB2::getAtoms() const
{
    return atoms;
}

/** Return the helices. */
QVector<PDBHelix> PDB2::getHelices() const
{
    return helices;
}

/** Return the sheets. */
QVector<PDBSheet> PDB2::getSheets() const
{
    return sheets;
}

/** Return the master record object. */
PDBMaster PDB2::getMaster() const
{
    return master;
}

/** Return the map of invalid records. */
QMap<qint64, QString> PDB2::getInvalidRecords() const
{
    return invalid_records;
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
    QMultiMap<QChar, QString> frame_chains;

    // The residue sequence (segment) numbers for the frame.
    QSet<qint64> frame_segments;

    // The connectivity map for the frame.
    QMultiMap<qint64, qint64> frame_connections;

    // Internal function used to parse a single atom line in the file.
    auto parse_atoms = [&](const QString &line, int iatm, int iframe, int iline, int num_lines,
        PDBAtom &frame_atom, QMultiMap<QPair<qint64, QString>, qint64> &local_residues,
        QMultiMap<QChar, QString> &local_chains, QSet<qint64> &local_segments, QStringList &errors)
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

            // There are anisotropic temperature factors for this atom.
            else if (lines()[iline + 1].leftRef(6) == "ANISOU")
            {
                frame_atom.setAnisTemp(line, lines().constData()[iline + 1], errors);
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
            local_chains.insert(frame_atom.getChainId(), frame_atom.getResName());

        // Insert the residue sequence number.
        local_segments.insert(frame_atom.getResNum());

        return frame_atom.isTer();
    };

    // Loop through all lines in the file.
    for (int iline=0; iline<lines().count(); ++iline)
    {
        // Store a reference to the line.
        const QString &line = lines()[iline];

        // Whether to parse atom data at the end of the current loop.
        bool isParse = false;

        // Whether a model section has just been parsed.
        // This stops atom data being parsed twice if we've just recorded a model
        // then hit the end of the file.
        bool isModel = false;

        // Extract the record type.
        // Could simplify this, i.e. remove whitespace.
        QString record = lines()[iline].left(6);

        // Parse TITLE section records.
        if      (record == "HEADER") title.appendRecord(line, PDBTitle::RECORD_TYPE::HEADER, parse_warnings);
        else if (record == "OBSLTE") title.appendRecord(line, PDBTitle::RECORD_TYPE::OBSLTE, parse_warnings);
        else if (record == "TITLE ") title.appendRecord(line, PDBTitle::RECORD_TYPE::TITLE,  parse_warnings);
        else if (record == "SPLIT ") title.appendRecord(line, PDBTitle::RECORD_TYPE::SPLIT,  parse_warnings);
        else if (record == "CAVEAT") title.appendRecord(line, PDBTitle::RECORD_TYPE::CAVEAT, parse_warnings);
        else if (record == "COMPND") title.appendRecord(line, PDBTitle::RECORD_TYPE::COMPND, parse_warnings);
        else if (record == "SOURCE") title.appendRecord(line, PDBTitle::RECORD_TYPE::SOURCE, parse_warnings);
        else if (record == "KEYWDS") title.appendRecord(line, PDBTitle::RECORD_TYPE::KEYWDS, parse_warnings);
        else if (record == "EXPDTA") title.appendRecord(line, PDBTitle::RECORD_TYPE::EXPDTA, parse_warnings);
        else if (record == "MDLTYP") title.appendRecord(line, PDBTitle::RECORD_TYPE::MDLTYP, parse_warnings);
        else if (record == "AUTHOR") title.appendRecord(line, PDBTitle::RECORD_TYPE::AUTHOR, parse_warnings);
        else if (record == "REVDAT") title.appendRecord(line, PDBTitle::RECORD_TYPE::REVDAT, parse_warnings);
        else if (record == "SPRSDE") title.appendRecord(line, PDBTitle::RECORD_TYPE::SPRSDE, parse_warnings);
        else if (record == "JRNL  ") title.appendRecord(line, PDBTitle::RECORD_TYPE::JRNL,   parse_warnings);
        else if (record == "REMARK") title.appendRecord(line, PDBTitle::RECORD_TYPE::REMARK, parse_warnings);
        else if (record == "NUMMDL") title.appendRecord(line, PDBTitle::RECORD_TYPE::NUMMDL, parse_warnings);

        // Start of a MODEL record.
        // These are used to define an atom configuratation, so can be used as
        // frames in a trajectory file. Each model entry must be consistent, i.e.
        // it must contain the same number and type of atoms.
        else if (record == "MODEL ")
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

            // Flag that a model has been recorded.
            isModel = true;
        }

        // An ATOM, or HETATM record.
        else if (record == "ATOM  " or
                 record == "HETATM")
        {
            // Store the line number of the atom record.
            atom_lines.append(iline);
            nats++;
        }

        // A HELIX record.
        else if (record == "HELIX ")
        {
            // Create a helix object.
            PDBHelix helix(line, parse_warnings);

            // Append the helix.
            helices.append(helix);
        }

        // A chain terminus record.
        // These are processed by the atom parsing section so we just skip ahead.
        else if (record == "TER   ");

        // A SHEET record.
        else if (record == "SHEET ")
        {
            // Create a sheet object.
            PDBSheet sheet(line, parse_warnings);

            // Append the sheet.
            sheets.append(sheet);
        }

        // A CRYST1 record (crystallographic data).
        else if (record == "CRYST1")
        {
            // Create the crystallographic data object.
            crystal = PDBCrystal(line, parse_warnings);
        }

        // ORIGXn transformation records.
        else if (record == "ORIGX1") trans_orig.appendRecord(line, 0, false, parse_warnings);
        else if (record == "ORIGX2") trans_orig.appendRecord(line, 1, false, parse_warnings);
        else if (record == "ORIGX3") trans_orig.appendRecord(line, 2, false, parse_warnings);

        // SCALEn transformation records.
        else if (record == "SCALE1") trans_scale.appendRecord(line, 0, false, parse_warnings);
        else if (record == "SCALE2") trans_scale.appendRecord(line, 1, false, parse_warnings);
        else if (record == "SCALE3") trans_scale.appendRecord(line, 2, false, parse_warnings);

        // MTRIXn transformation records.
        else if (record == "MTRIX1") trans_matrix.appendRecord(line, 0, true, parse_warnings);
        else if (record == "MTRIX2") trans_matrix.appendRecord(line, 1, true, parse_warnings);
        else if (record == "MTRIX3") trans_matrix.appendRecord(line, 2, true, parse_warnings);

        // MASTER record.
        else if (record == "MASTER")
        {
            // There should only be a single MASTER record and it should only refer to
            // the first MODEL entry.
            if (has_master)
            {
                parse_warnings.append(QObject::tr("Ignoring invalid MASTER record "
                    "on line %1: '%2'")
                    .arg(iline).arg(line));
            }

            master = PDBMaster(line, parse_warnings);
            has_master = true;
        }

        // CONECT record.
        else if (record == "CONECT")
        {
            // There are a maximum of four bonds per line. If there are more
            // than four bonds per atom then a second CONECT entry can be used.

            // It could be hard to check for errors, since we don't know how many
            // bonds each atom should have. However, the connectivity information
            // is redundant, i.e. each bond entry appears twice, so we can validate
            // by constructing an adjacency matrix afterwards and ensuring that it
            // is symmetric.

            // First extract the atom serial number.
            bool ok;
            int atom = line.midRef(6,5).toInt(&ok);

            if (not ok)
            {
                parse_warnings.append(QObject::tr("Cannot extract the atom serial number "
                    "from part (%1) from line '%2'")
                    .arg(line.mid(6,5)).arg(line));

                return;
            }

            // Loop over all possible bond records.
            for (int i=0; i<4; ++i)
            {
                bool ok;
                int bond = line.midRef(11 + 5*i, 5).toInt(&ok);

                // Add the bond to the connectivity map.
                if (ok)
                {
                    frame_connections.insert(atom, bond);
                }
            }
        }

        // Invalid record
        else
        {
            parse_warnings.append(QObject::tr("Invalid PDB record found on "
                "line %1: '%2'").arg(iline).arg(line));

            invalid_records[iline] = line;
        }

        // End of the file.
        if (iline + 1 == lines().count())
            isParse = true;

        // Parse the atom data.
        if (isParse ^ isModel)
        {
            // Don't proceed if there are no more ATOM records to parse.
            if (atom_lines.count() > 0)
            {
                // Initialise atom vector for the frame.
                QVector<PDBAtom> frame_atoms(nats);

                if (usesParallel())
                {
                    QMutex mutex;

                    tbb::parallel_for( tbb::blocked_range<int>(0, nats),
                                    [&](const tbb::blocked_range<int> &r)
                    {
                        qint64 local_num_ters = 0;
                        QStringList local_errors;
                        QMultiMap<QPair<qint64, QString>, qint64> local_residues;
                        QMultiMap<QChar, QString> local_chains;
                        QSet<qint64> local_segments;

                        for (int i=r.begin(); i<r.end(); ++i)
                        {
                            // Parse the atom record and determine whether it is a terminal record.
                            if (parse_atoms( lines().constData()[atom_lines[i]], i, iframe,
                                atom_lines[i], lines().count(), frame_atoms[i], local_residues,
                                local_chains, local_segments, local_errors ))
                                    local_num_ters++;
                        }

                        QMutexLocker lkr(&mutex);

                        num_ters       += local_num_ters;
                        frame_residues += local_residues;
                        frame_chains   += local_chains;
                        frame_segments += local_segments;
                        parse_warnings += local_errors;
                    });
                }
                else
                {
                    for (int i=0; i<nats; ++i)
                    {
                        // Parse the atom record and determine whether it is a terminal record.
                        if (parse_atoms( lines().constData()[atom_lines[i]], i, iframe,
                            atom_lines[i], lines().count(), frame_atoms[i], frame_residues,
                            frame_chains, frame_segments, parse_warnings ))
                                num_ters++;
                    }
                }

                // Append frame data and clear the vectors.

                atoms.append(frame_atoms);
                atom_lines.clear();

                residues.append(frame_residues);
                frame_residues.clear();

                chains.append(frame_chains);
                frame_chains.clear();

                segments.append(frame_segments);
                frame_segments.clear();

                connections.append(frame_connections);
                frame_connections.clear();

                iframe++;
                nats = 0;
            }
        }
    }

    // Validate the parsed data against the MASTER record.
    if (has_master)
    {
        if (not (master.nRemarks() == title.nRemarks()))
        {
            parse_warnings.append(QObject::tr("Error in PDB file as the number of "
               "REMARK records (%1) read does not equal the number in "
               "the MASTER record (%2)!")
                .arg(title.nRemarks()).arg(master.nRemarks()));
        }

        if (not (master.nAtoms() == nAtoms()))
        {
            parse_warnings.append(QObject::tr("Error in PDB file as the number of "
               "ATOM and HETATM records (%1) read does not equal the number in "
               "the MASTER record (%2)!")
                .arg(nAtoms()).arg(master.nAtoms()));
        }

        if (not (master.nHelices() == nHelices()))
        {
            parse_warnings.append(QObject::tr("Error in PDB file as the number of "
               "HELIX records (%1) read does not equal the number in "
               "the MASTER record (%2)!")
                .arg(nHelices()).arg(master.nHelices()));
        }

        if (not (master.nSheets() == nSheets()))
        {
            parse_warnings.append(QObject::tr("Error in PDB file as the number of "
               "SHEET records (%1) read does not equal the number in "
               "the MASTER record (%2)!")
                .arg(nSheets()).arg(master.nSheets()));
        }

        if (not (master.nTers() == nTers()))
        {
            parse_warnings.append(QObject::tr("Error in PDB file as the number of "
               "TER records (%1) read does not equal the number in "
               "the MASTER record (%2)!")
                .arg(nTers()).arg(master.nTers()));
        }

        // Work out the number of coordinate transformation records.
        // Each object should contain three records, one for each dimension.
        int num_transforms = 3 * ( hasTransOrig() + hasTransScale() + hasTransMatrix() );

        if (not (master.nTransforms() == num_transforms ))
        {
            parse_warnings.append(QObject::tr("Error in PDB file as the number of coordinate "
               "transformation records (%1) read does not equal the number in "
               "the MASTER record (%2)!")
                .arg(num_transforms).arg(master.nTransforms()));
        }
    }

    this->setScore(nats);
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
    assigning properties based on the mapping in 'map', for the configuration
    'iframe' (a trajectory index). */
System PDB2::startSystem(int iframe, const PropertyMap &map) const
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
                mols_array[i] = this->getMolecule(i, iframe, map);
            }
        });
    }
    else
    {
        for (int i=0; i<nmols; ++i)
        {
            mols_array[i] = this->getMolecule(i, iframe, map);
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

/** Use the data contained in this parser to create a new System of molecules,
    assigning properties based on the mapping in 'map'. */
System PDB2::startSystem(const PropertyMap &map) const
{
    // Generate system for single model PDF files.
    return startSystem(0, map);
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

/** Internal function used to get the molecule structure for molecule 'imol'
    in the frame (model) 'iframe'. */
MolStructureEditor PDB2::getMolStructure(int imol,
    int iframe, const PropertyName &cutting) const
{
    // Make sure the frame index is within range.
    if ((iframe < 0) or (iframe > atoms.count()))
    {
        throw SireError::program_bug(QObject::tr("The frame index %1 is out of "
            "range, 0 - %2").arg(iframe).arg(atoms.count()), CODELOC);
    }

    // Make sure that there are atoms in the molecule.
    if (atoms[iframe].count() == 0)
    {
        throw SireError::program_bug(QObject::tr(
            "Strange - there are no atoms in molecule %1 in frame %2?")
                .arg(imol).arg(iframe), CODELOC);
    }

    // First step is to build the structure of the molecule, i.e.
    // the layout of cutgroups, residues and atoms.
    MolStructureEditor mol;

    // To do this we'll walk through all of the residues in the frame,
    // adding each atom contained within the residue.

    // Residue index.
    int ires = 0;

    // Create a reverse mapping between residues and chains.
    QMap<QString, QChar> res_to_chain;

    // Add any chains to the molecule.
    for (auto chain : chains[iframe].uniqueKeys())
    {
        mol.add(ChainName(chain));

        // Get a list of the residues that are part of this chain.
        QList<QString> chain_residues = chains[iframe].values(chain);

        // Create the reverse mapping.
        for (auto residue : chain_residues)
            res_to_chain.insert(residue, chain);
    }

    // Add any segments to the molecule.
    for (auto segment : segments[iframe])
        mol.add(SegName(QString::number(segment)));

    // Loop over all unique residues in the frame.
    for (auto residue : residues[iframe].uniqueKeys())
    {
        // Extract the residue number and name.
        auto res_num  = residue.first;
        auto res_name = residue.second;

        // By default we will use one CutGroup per residue.
        // This may be changed later by the cutting system.
        auto cutgroup = mol.add(CGName(QString::number(ires)));

        // Get a sorted list of the atoms that are part of the residue.
        QList<qint64> res_atoms = residues[iframe].values(residue);
        qSort(res_atoms);

        // Add the residue to the molecule.
        auto res = mol.add(ResNum(res_num));
        res.rename(ResName(res_name.trimmed()));

        // Reparent the residue to its chain.
        if (res_to_chain.contains(res_name))
            res.reparent(ChainName(res_to_chain[res_name]));

        // Add each atom in the residue to the molecule.
        for (auto res_atom : res_atoms)
        {
            auto atom = cutgroup.add(AtomNum(atoms[iframe][res_atom].getSerial()));
            atom.rename(AtomName(atoms[iframe][res_atom].getName().trimmed()));

            // Reparent the atom to its residue and segment.
            atom.reparent(ResNum(res_num));
            atom.reparent(SegName(QString::number(atoms[iframe][res_atom].getResNum())));
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

/** Internal function used to get the molecule structure for molecule 'imol'
    in the frame (model) 'iframe'. */
MolEditor PDB2::getMolecule(int imol, int iframe, const PropertyMap &map) const
{
    // At the moment we'll assume that there is a single molecule. Once the molecule
    // is constructed we can use connectivity information to break it into sub-molecules.

    // Make sure the frame index is within range.
    if ((iframe < 0) or (iframe > atoms.count()))
    {
        throw SireError::program_bug(QObject::tr("The frame index %1 is out of "
            "range, 0 - %2").arg(iframe).arg(atoms.count()), CODELOC);
    }

    // Make sure that there are atoms in the frame.
    if (atoms[iframe].count() == 0)
    {
        throw SireError::program_bug(QObject::tr(
            "Strange - there are no atoms in molecule %1 in frame %2?")
                .arg(imol).arg(iframe), CODELOC);
    }

    // First, construct the layout of the molecule (sorting of atoms into residues and cutgroups).
    auto mol = this->getMolStructure(imol, iframe, map["cutting"]).commit().edit();

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
    for (int i=0; i<nAtoms(); ++i)
    {
        // Store a reference to the current atom.
        const PDBAtom &atom = atoms[iframe][i];

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

    QVector<bool> is_ter(num_atoms);
    is_ter.fill(false);

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
                local_atoms[i] = PDBAtom(sire_mol.atom(AtomIdx(i)), local_errors);
                atom_lines[i] = local_atoms[i].toPDBRecord();

                // Flag whether this is a terminal atom.
                if (num_chains > 0) is_ter[i] = local_atoms[i].isTer();
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
            PDBAtom atom(sire_mol.atom(AtomIdx(i)), errors);

            // Generate a PDB atom data record.
            atom_lines[iline] = atom.toPDBRecord();
            iline++;

            // Add a TER record for this atom.
            if (atom.isTer())
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
