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

#include <iostream>

#include "SireIO/pdb2.h"

#include "SireMM/pdbparams.h"

#include "SireSystem/system.h"

#include "SireBase/parallel.h"

#include "SireError/errors.h"
#include "SireIO/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireIO;
using namespace SireMM;
using namespace SireMol;
using namespace SireBase;
using namespace SireSystem;
using namespace SireStream;

const RegisterParser<PDB2> register_pdb;
static const RegisterMetaType<PDB2> r_pdb2;
static const RegisterMetaType<PDBAtom> r_pdbatom(NO_ROOT);
static const RegisterMetaType<PDBHelix> r_pdbhelix(NO_ROOT);
static const RegisterMetaType<PDBSheet> r_pdbsheet(NO_ROOT);
static const RegisterMetaType<PDBTitle> r_pdbtitle(NO_ROOT);

QDataStream SIREIO_EXPORT &operator<<(QDataStream &ds, const PDBAtom &pdbatom)
{
    writeHeader(ds, r_pdbatom, 1);

    SharedDataStream sds(ds);

    sds << pdbatom.record << pdbatom.serial << pdbatom.name << pdbatom.alt_loc
        << pdbatom.res_name << pdbatom.chain_id << pdbatom.res_num << pdbatom.insert_code
        << pdbatom.coord << pdbatom.occupancy << pdbatom.temperature << pdbatom.element
        << pdbatom.charge << pdbatom.isHet << pdbatom.isTer << pdbatom.isAnis
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
            >> pdbatom.charge >> pdbatom.isHet >> pdbatom.isTer >> pdbatom.isAnis
            >> pdbatom.anis_facts[0] >> pdbatom.anis_facts[1] >> pdbatom.anis_facts[2]
            >> pdbatom.anis_facts[3] >> pdbatom.anis_facts[4] >> pdbatom.anis_facts[5];
    }
    else
        throw version_error(v, "1", r_pdbatom, CODELOC);

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

QDataStream SIREIO_EXPORT &operator<<(QDataStream &ds, const PDB2 &pdb2)
{
    writeHeader(ds, r_pdb2, 1);

    SharedDataStream sds(ds);

    sds << pdb2.title << pdb2.atoms << pdb2.helices << pdb2.sheets << pdb2.invalid_records
        << pdb2.parse_warnings;

    return ds;
}

QDataStream SIREIO_EXPORT &operator>>(QDataStream &ds, PDB2 &pdb2)
{
    VersionID v = readHeader(ds, r_pdb2);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> pdb2.title >> pdb2.atoms >> pdb2.helices >> pdb2.sheets >> pdb2.invalid_records
            >> pdb2.parse_warnings;
    }
    else
        throw version_error(v, "1", r_pdb2, CODELOC);

    return ds;
}

/** Default constructor. */
PDBAtom::PDBAtom() : occupancy(1.0),
                     element("X"),
                     charge("0"),
                     isHet(false),
                     isTer(false),
                     isAnis(false)
{
}

/** Constructor. */
PDBAtom::PDBAtom(const QString &line, QStringList &errors)
                    : record(line),
                      occupancy(1.0),
                      element("X"),
                      charge("0"),
                      isHet(false),
                      isTer(false),
                      isAnis(false)
{
    if (line.length() < 54)
    {
        errors.append(QObject::tr("Cannot parse ATOM record "
            " since it does not match the format! '%1'").arg(line));

        return;
    }

    // Flag that this is a HETATM record.
    if (line.leftRef(6) == "HETATM")
        isHet = true;

    // Extract the atom serial number.
    bool ok;
    int tmp_int = line.midRef(6,5).toInt(&ok);

    if (ok)
    {
        serial = tmp_int;
    }
    else
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
    tmp_int = line.midRef(22,4).toInt(&ok);

    if (ok)
    {
        res_num = tmp_int;
    }
    else
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
    charge = line.mid(78,2);
}

/** Set the terminal atom flag. */
void PDBAtom::setTerminal(bool _isTer)
{
    isTer = _isTer;
}

/** Set anisotropic temperature factor data. */
void PDBAtom::setAnisTemp(const QString &line1, const QString &line2, QStringList &errors)
{
    isAnis = true;

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
QString PDBAtom::getChainId() const
{
    return chain_id;
}

/** Get the residue sequence number. */
qint64 PDBAtom::getResNum() const
{
    return res_num;
}

/** Return the C++ name for this class */
const char* PDBAtom::typeName()
{
    return QMetaType::typeName( qMetaTypeId<PDBAtom>() );
}

/** Return a PDB record line for this atom. */
QString PDBAtom::toPDBLine() const
{
    return record;
}

/** Default constructor. */
PDBHelix::PDBHelix()
{
}

/** Constructor. */
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
QString PDBHelix::toPDBLine() const
{
    return record;
}

/** Default constructor. */
PDBSheet::PDBSheet()
{
}

/** Constructor. */
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
QString PDBSheet::toPDBLine() const
{
    return record;
}

/** Default constructor. */
PDBTitle::PDBTitle()
{
}

/** Append a record to the title object. */
void PDBTitle::appendRecord(const QString &record,
                            RECORD_TYPE record_type,
                            QStringList &errors)
{
    switch(record_type)
    {
        case RECORD_TYPE::HEADER: header = record;
        case RECORD_TYPE::OBSLTE: obsoletes.append(record);
        case RECORD_TYPE::TITLE : titles.append(record);
        case RECORD_TYPE::SPLIT : splits.append(record);
        case RECORD_TYPE::CAVEAT: caveats.append(record);
        case RECORD_TYPE::COMPND: compounds.append(record);
        case RECORD_TYPE::SOURCE: sources.append(record);
        case RECORD_TYPE::KEYWDS: keywords.append(record);
        case RECORD_TYPE::EXPDTA: experiments.append(record);
        case RECORD_TYPE::MDLTYP: model_types.append(record);
        case RECORD_TYPE::AUTHOR: authors.append(record);
        case RECORD_TYPE::REVDAT: revisions.append(record);
        case RECORD_TYPE::SPRSDE: supersedes.append(record);
        case RECORD_TYPE::JRNL  : journals.append(record);
        case RECORD_TYPE::REMARK: remarks.append(record);
        case RECORD_TYPE::NUMMDL:
        {
            bool ok;
            int tmp_int = record.midRef(10,4).toInt(&ok);

            if (not ok)
            {
                errors.append(QObject::tr("Cannot extract the number of models "
                    "from part (%1) from line '%2'")
                    .arg(record.mid(10,4)).arg(record));

                return;
            }
            num_models = tmp_int;
        }
        default:
        {
            errors.append(QObject::tr("Title section record not recognised! "
                "'%2'").arg(record));

            return;
        }
    }
}

/** Return the C++ name for this class */
const char* PDBTitle::typeName()
{
    return QMetaType::typeName( qMetaTypeId<PDBTitle>() );
}

/** Constructor */
PDB2::PDB2() : ConcreteProperty<PDB2,MoleculeParser>()
{
}

/** Construct to read in the data from the file called 'filename'. The
    passed property map can be used to pass extra parameters to control
    the parsing */
PDB2::PDB2(const QString &filename, const PropertyMap &map)
     : ConcreteProperty<PDB2,MoleculeParser>(filename,map)
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
PDB2::PDB2(const QStringList &lines, const PropertyMap &map)
     : ConcreteProperty<PDB2,MoleculeParser>(lines,map)
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
PDB2::PDB2(const SireSystem::System &system, const PropertyMap &map)
     : ConcreteProperty<PDB2,MoleculeParser>(map)
{
    //look through the system and extract out the information
    //that is needed to go into the file, based on the properties
    //found in 'map'. Do this to create the set of text lines that
    //will make up the file

    QStringList lines;  // = code used to generate the lines

    //now that you have the lines, reparse them back into a PDB2 object,
    //so that the information is consistent, and you have validated that
    //the lines you have written are able to be correctly read in by
    //this parser. This will also implicitly call 'assertSane()'
    PDB2 parsed( lines, map );

    this->operator=(parsed);
}

/** Copy constructor */
PDB2::PDB2(const PDB2 &other)
     : ConcreteProperty<PDB2,MoleculeParser>(other)
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

/** Function that is called to assert that this object is sane. This
    should raise an exception if the parser is in an invalid state */
void PDB2::assertSane() const
{
    //check state, raise SireError::program_bug if we are in an invalid state
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

    // Internal function used to parse a single atom line in the file.
    auto parse_atoms = [&](const QString &line, int iatm, int iframe,
        int iline, int num_lines, PDBAtom &frame_atom, QStringList &errors)
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

        // Validate the atom data.
        if (iframe > 1)
        {
            if (not validateAtom(frame_atom, atoms[iframe-1][iatm]))
            {
                errors.append( QObject::tr( "Invalid record for atom %1 "
                    "in frame %2 on line %3. Record doesn't match does not "
                    " match previous frame! '%4'")
                    .arg(iatm).arg(iframe).arg(iline).arg(line));
            }

        }
    };

    // Loop through all lines in the file.
    for (int iline=0; iline<lines().count(); ++iline)
    {
        // Whether to parse atom data at the end of the current loop.
        bool isParse = false;

        // Whether a model section has just been parsed.
        // This stops atom data being parsed twice if we've just recorded a model
        // then hit the end of the file.
        bool isModel = false;

        // Extract the record type.
        QString record = lines()[iline].left(6);

        // Parse TITLE section records.
        if      (record == "HEADER") title.appendRecord(lines()[iline], PDBTitle::RECORD_TYPE::HEADER, parse_warnings);
        else if (record == "OBSLTE") title.appendRecord(lines()[iline], PDBTitle::RECORD_TYPE::OBSLTE, parse_warnings);
        else if (record == "TITLE ") title.appendRecord(lines()[iline], PDBTitle::RECORD_TYPE::TITLE,  parse_warnings);
        else if (record == "SPLIT ") title.appendRecord(lines()[iline], PDBTitle::RECORD_TYPE::SPLIT,  parse_warnings);
        else if (record == "CAVEAT") title.appendRecord(lines()[iline], PDBTitle::RECORD_TYPE::CAVEAT, parse_warnings);
        else if (record == "COMPND") title.appendRecord(lines()[iline], PDBTitle::RECORD_TYPE::COMPND, parse_warnings);
        else if (record == "KEYWDS") title.appendRecord(lines()[iline], PDBTitle::RECORD_TYPE::KEYWDS, parse_warnings);
        else if (record == "EXPDTA") title.appendRecord(lines()[iline], PDBTitle::RECORD_TYPE::EXPDTA, parse_warnings);
        else if (record == "MDLTYP") title.appendRecord(lines()[iline], PDBTitle::RECORD_TYPE::MDLTYP, parse_warnings);
        else if (record == "AUTHOR") title.appendRecord(lines()[iline], PDBTitle::RECORD_TYPE::AUTHOR, parse_warnings);
        else if (record == "REVDAT") title.appendRecord(lines()[iline], PDBTitle::RECORD_TYPE::REVDAT, parse_warnings);
        else if (record == "SPRSDE") title.appendRecord(lines()[iline], PDBTitle::RECORD_TYPE::SPRSDE, parse_warnings);
        else if (record == "JRNL  ") title.appendRecord(lines()[iline], PDBTitle::RECORD_TYPE::JRNL,   parse_warnings);
        else if (record == "REMARK") title.appendRecord(lines()[iline], PDBTitle::RECORD_TYPE::REMARK, parse_warnings);
        else if (record == "NUMMDL") title.appendRecord(lines()[iline], PDBTitle::RECORD_TYPE::NUMMDL, parse_warnings);

        // Start of a MODEL record.
        // These are used to define an atom configuratation, so can be used as
        // frames in a trajectory file. Each model entry must be consistent, i.e.
        // it must contain the same number and type of atoms.
        if (record == "MODEL ")
        {
            imdl++;

            // Extract the model entry number.
            int nmod = lines()[iline].midRef(10,4).toInt();

            // Check that the model entry number is correct.
            // These must be in ascending order, starting at 1.
            if (imdl != nmod)
            {
                parse_warnings.append( QObject::tr( "Cannot parse the data "
                    "for MODEL %1 has incorrect MODEL entry number '%2'!")
                    .arg(imdl).arg(nmod) );

                return;
            }
        }

        // End of a MODEL record.
        else if (record == "ENDMDL ")
        {
            if (imdl > 1)
            {
                // Check that the atom number is consisent.
                if (nats != prev_nats)
                {
                    parse_warnings.append( QObject::tr( "Cannot parse the data "
                    "for MODEL %1 is not the same size as MODEL %2!")
                    .arg(imdl).arg(imdl-1) );

                    return;
                }
            }

            // Record the number of atoms in last model entry.
            prev_nats = nats;

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
            PDBHelix helix(lines().constData()[iline], parse_warnings);

            // Append the helix.
            helices.append(helix);
        }

        // A SHEET record.
        else if (record == "SHEET ")
        {
            // Create a sheet object.
            PDBSheet sheet(lines().constData()[iline], parse_warnings);

            // Append the sheet.
            sheets.append(sheet);
        }

        // Invalid record
        else
        {
            parse_warnings.append( QObject::tr( "Invalid PDB record found on "
                "line %1: '%2'").arg(iline).arg(lines()[iline]));

            invalid_records[iline] = lines()[iline];
        }

        // End of the file.
        if (iline + 1 == lines().count())
            isParse = true;

        // Parse the atom data.
        if (isParse or isModel)
        {
            if (nats > 99999)
            {
                parse_warnings.append( QObject::tr( "Number of atoms exceeds the PDB file limit "
                    "of 99999. Please split large entries into multiple files"));

                return;
            }

            // Initialise atom vector for the frame.
            QVector<PDBAtom> frame_atoms(nats);

            if (usesParallel())
            {
                QMutex mutex;

                tbb::parallel_for( tbb::blocked_range<int>(0,nats),
                                [&](const tbb::blocked_range<int> &r)
                {
                    QStringList local_errors;

                    for (int i=r.begin(); i<r.end(); ++i)
                    {
                        parse_atoms( lines().constData()[atom_lines[i]], i, iframe,
                            atom_lines[i], lines().count(), frame_atoms[i], local_errors );
                    }

                    if (not local_errors.isEmpty())
                    {
                        QMutexLocker lkr(&mutex);
                        parse_warnings += local_errors;
                    }
                });
            }
            else
            {
                for (int i=0; i<nats; ++i)
                {
                    parse_atoms( lines().constData()[atom_lines[i]], i, iframe,
                        atom_lines[i], lines().count(), frame_atoms[i], parse_warnings );
                }
            }

            // Append the atom data.
            atoms.append(frame_atoms);

            // Clear the atom line index vector.
            atom_lines.clear();

            iframe++;
            nats = 0;
        }
    }

    //for (int i=0; i<title.titles.size(); ++i)
        //std::cout << title.titles[i].toStdString() << '\n';

    /*for (int i=0; i<atoms[0].size(); ++i)
    {
        std::cout << atoms[0][i].serial << ' '
                  << atoms[0][i].name.toStdString() << ' '
                  << atoms[0][i].alt_loc.toLatin1() << ' '
                  << atoms[0][i].res_name.toStdString() << ' '
                  << atoms[0][i].chain_id.toLatin1() << ' '
                  << atoms[0][i].res_num << ' '
                  << atoms[0][i].insert_code.toLatin1() << ' '
                  << atoms[0][i].coord[0] << ' '
                  << atoms[0][i].coord[1] << ' '
                  << atoms[0][i].coord[2] << ' '
                  << atoms[0][i].occupancy << ' '
                  << atoms[0][i].temperature << ' '
                  << atoms[0][i].element.toStdString() << ' '
                  << atoms[0][i].charge.toStdString() << '\n';
        std::cin.get();
    }*/

    /*for (int i=0; i<helices.size(); ++i)
    {
        std::cout << helices[i].serial << ' '
                  << helices[i].id.toStdString() << ' '
                  << helices[i].init_res_name.toStdString() << ' '
                  << helices[i].init_chain_id.toLatin1() << ' '
                  << helices[i].init_res_num << ' '
                  << helices[i].init_insert_code.toLatin1() << ' '
                  << helices[i].end_res_name.toStdString() << ' '
                  << helices[i].end_chain_id.toLatin1() << ' '
                  << helices[i].end_res_num << ' '
                  << helices[i].end_insert_code.toLatin1() << ' '
                  << helices[i].helix_class << ' '
                  << helices[i].comment.toStdString() << ' '
                  << helices[i].length << '\n';
        std::cin.get();
    }*/

    /*for (int i=0; i<sheets.size(); ++i)
    {
        std::cout << sheets[i].strand << ' '
                  << sheets[i].id.toStdString() << ' '
                  << sheets[i].init_res_name.toStdString() << ' '
                  << sheets[i].init_chain_id.toLatin1() << ' '
                  << sheets[i].init_res_num << ' '
                  << sheets[i].init_insert_code.toLatin1() << ' '
                  << sheets[i].end_res_name.toStdString() << ' '
                  << sheets[i].end_chain_id.toLatin1() << ' '
                  << sheets[i].end_res_num << ' '
                  << sheets[i].end_insert_code.toLatin1() << ' '
                  << sheets[i].sense << ' ';
        if (i > 0)
        {
            std::cout << sheets[i].curr_atm_name.toStdString() << ' '
                      << sheets[i].curr_res_name.toStdString() << ' '
                      << sheets[i].curr_chain_id.toLatin1() << ' '
                      << sheets[i].curr_res_num << ' '
                      << sheets[i].curr_insert_code.toLatin1() << ' '
                      << sheets[i].prev_atm_name.toStdString() << ' '
                      << sheets[i].prev_res_name.toStdString() << ' '
                      << sheets[i].prev_chain_id.toLatin1() << ' '
                      << sheets[i].prev_res_num << ' '
                      << sheets[i].prev_insert_code.toLatin1() << ' ';
        }
        std::cout << '\n';
        std::cin.get();
    }*/

    num_atom = atoms[0].count();
    num_helix = helices.count();
    num_sheet = sheets.count();
    this->setScore(nats);
}

/** Helper function used to validate atom data from different model records */
bool PDB2::validateAtom(const PDBAtom &atom1, const PDBAtom &atom2)
{
    // Different models are typically indexed by atom number, e.g.
    // 1 - nats, nats+1 - 2*nats, ..., or by using the alternate location
    // entry of the atom record, e.g. AALA in the first model, BALA in the second.

    // The following data must be the same for all models (I think...)
    if (atom1.getName()    != atom2.getName()    or
        atom1.getResName() != atom2.getResName() or
        atom1.getChainId() != atom2.getChainId() or
        atom1.getResNum()  != atom2.getResNum())
    {
        return false;
    }
    else return true;
}

/** Use the data contained in this parser to create a new System of molecules,
    assigning properties based on the mapping in 'map' */
System PDB2::startSystem(const PropertyMap &map) const
{
    //you should now take the data that you have parsed into the intermediate
    //format and use it to start a new System and create all molecules,
    //which are added to this System

    System system;


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
