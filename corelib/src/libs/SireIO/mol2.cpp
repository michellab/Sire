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

#include "SireIO/mol2.h"

#include "SireMM/mol2params.h"

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

const RegisterParser<Mol2> register_mol2;
static const RegisterMetaType<Mol2> r_mol2;
static const RegisterMetaType<Mol2Atom> r_mol2atom(NO_ROOT);
static const RegisterMetaType<Mol2Bond> r_mol2bond(NO_ROOT);
static const RegisterMetaType<Mol2Molecule> r_mol2molecule(NO_ROOT);
static const RegisterMetaType<Mol2Substructure> r_mol2subst(NO_ROOT);

QDataStream SIREIO_EXPORT &operator<<(QDataStream &ds, const Mol2Atom &mol2atom)
{
    writeHeader(ds, r_mol2atom, 1);

    SharedDataStream sds(ds);

    sds << mol2atom.record << mol2atom.number << mol2atom.name << mol2atom.coord
        << mol2atom.type << mol2atom.subst_id << mol2atom.subst_name
        << mol2atom.charge << mol2atom.status_bit;

    return ds;
}

QDataStream SIREIO_EXPORT &operator>>(QDataStream &ds, Mol2Atom &mol2atom)
{
    VersionID v = readHeader(ds, r_mol2atom);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> mol2atom.record >> mol2atom.number >> mol2atom.name >> mol2atom.coord
            >> mol2atom.type >> mol2atom.subst_id >> mol2atom.subst_name
            >> mol2atom.charge >> mol2atom.status_bit;
    }
    else
        throw version_error(v, "1", r_mol2atom, CODELOC);

    return ds;
}

QDataStream SIREIO_EXPORT &operator<<(QDataStream &ds, const Mol2Bond &mol2bond)
{
    writeHeader(ds, r_mol2bond, 1);

    SharedDataStream sds(ds);

    sds << mol2bond.record << mol2bond.number << mol2bond.origin
        << mol2bond.target << mol2bond.type << mol2bond.subst_id
        << mol2bond.status_bit;

    return ds;
}

QDataStream SIREIO_EXPORT &operator>>(QDataStream &ds, Mol2Bond &mol2bond)
{
    VersionID v = readHeader(ds, r_mol2bond);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> mol2bond.record >> mol2bond.number >> mol2bond.origin
            >> mol2bond.target >> mol2bond.type >> mol2bond.subst_id
            >> mol2bond.status_bit;
    }
    else
        throw version_error(v, "1", r_mol2bond, CODELOC);

    return ds;
}

QDataStream SIREIO_EXPORT &operator<<(QDataStream &ds, const Mol2Molecule &mol2molecule)
{
    writeHeader(ds, r_mol2molecule, 1);

    SharedDataStream sds(ds);

    sds << mol2molecule.record << mol2molecule.name << mol2molecule.num_atoms
        << mol2molecule.num_bonds << mol2molecule.num_subst << mol2molecule.num_feats
        << mol2molecule.num_sets << mol2molecule.mol_type << mol2molecule.charge_type
        << mol2molecule.status_bit << mol2molecule.comment << mol2molecule.atoms
        << mol2molecule.bonds << mol2molecule.substructures;

    return ds;
}

QDataStream SIREIO_EXPORT &operator>>(QDataStream &ds, Mol2Molecule &mol2molecule)
{
    VersionID v = readHeader(ds, r_mol2molecule);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> mol2molecule.record >> mol2molecule.name >> mol2molecule.num_atoms
            >> mol2molecule.num_bonds >> mol2molecule.num_subst >> mol2molecule.num_feats
            >> mol2molecule.num_sets >> mol2molecule.mol_type >> mol2molecule.charge_type
            >> mol2molecule.status_bit >> mol2molecule.comment >> mol2molecule.atoms
            >> mol2molecule.bonds >> mol2molecule.substructures;
    }
    else
        throw version_error(v, "1", r_mol2molecule, CODELOC);

    return ds;
}

QDataStream SIREIO_EXPORT &operator<<(QDataStream &ds, const Mol2Substructure &mol2subst)
{
    writeHeader(ds, r_mol2subst, 1);

    SharedDataStream sds(ds);

    sds << mol2subst.record << mol2subst.number << mol2subst.name << mol2subst.root_atom
        << mol2subst.type << mol2subst.dict_type << mol2subst.chain << mol2subst.sub_type
        << mol2subst.num_inter_bonds << mol2subst.status_bit << mol2subst.comment;

    return ds;
}

QDataStream SIREIO_EXPORT &operator>>(QDataStream &ds, Mol2Substructure &mol2subst)
{
    VersionID v = readHeader(ds, r_mol2subst);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> mol2subst.record >> mol2subst.number >> mol2subst.name >> mol2subst.root_atom
            >> mol2subst.type >> mol2subst.dict_type >> mol2subst.chain >> mol2subst.sub_type
            >> mol2subst.num_inter_bonds >> mol2subst.status_bit >> mol2subst.comment;
    }
    else
        throw version_error(v, "1", r_mol2subst, CODELOC);

    return ds;
}

QDataStream SIREIO_EXPORT &operator<<(QDataStream &ds, const Mol2 &mol2)
{
    writeHeader(ds, r_mol2, 1);

    ds << mol2.molecules << static_cast<const MoleculeParser&>(mol2);

    return ds;
}

QDataStream SIREIO_EXPORT &operator>>(QDataStream &ds, Mol2 &mol2)
{
    VersionID v = readHeader(ds, r_mol2);

    if (v == 1)
    {
        ds >> mol2.molecules >> static_cast<MoleculeParser&>(mol2);
    }
    else
        throw version_error(v, "1", r_mol2, CODELOC);

    return ds;
}

/** Default constructor. */
Mol2Atom::Mol2Atom() :
    name("X"),
    charge(0)
{
}

/** Constructor. */
Mol2Atom::Mol2Atom(const QString &line, QStringList &errors) :
    record(line),
    name("X"),
    charge(0)
{
    // Tokenize the string, splitting using a single whitespace characters.
    QStringList data = line.simplified().split(QRegExp("\\s"));

    // The Mol2 record indicator. All record types start with this prefix.
    const QString record_indicator("@<TRIPOS>");

    // Check that we've not hit another TRIPOS record indicator.
    // This can happen with poorly formatted files, e.g. the number of atoms is incorrect.
    if (data[0].left(9) == record_indicator)
    {
        errors.append(QObject::tr("The @<TRIPOS>ATOM record is incorrectly formatted!"));

        return;
    }

    // Extract the atom ID.
    bool ok;
    number = data[0].toInt(&ok);

    if (not ok)
    {
        errors.append(QObject::tr("Cannot extract the atom ID number "
            "from part (%1) from line '%2'").arg(data[0]).arg(line));

        return;
    }

    // Extract the atom name.
    name = data[1].simplified();

    bool ok_x, ok_y, ok_z;
    double x = data[2].toDouble(&ok_x);
    double y = data[3].toDouble(&ok_y);
    double z = data[4].toDouble(&ok_z);

    if (not (ok_x and ok_y and ok_z))
    {
        errors.append(QObject::tr("There was a problem reading the coordinate "
            "values of x, y, and z from the data '%1' in line '%2'")
            .arg(data[2] + ' ' + data[3] + ' ' + data[4]).arg(line));

        return;
    }
    coord = Vector(x, y, z);

    // Extract the SYBYL atom type.
    type = data[5].simplified();

    // Extract the substructure ID.
    subst_id = data[6].toInt(&ok);

    if (not ok)
    {
        errors.append(QObject::tr("Cannot extract the substructure ID number "
            "from part (%1) from line '%2'").arg(data[6]).arg(line));

        return;
    }

    // Extract the substructure name.
    subst_name = data[7].simplified();

    // Extract the charge on the atom.
    charge = data[8].toDouble(&ok);

    if (not ok)
    {
        errors.append(QObject::tr("Cannot extract the atom charge "
            "from part (%1) from line '%2'").arg(data[8]).arg(line));

        return;
    }

    // There is a SYBYL status bit entry for this atom.
    if (data.count() == 10)
    {
        // List of valid SYBYL atom status bits.
        QStringList valid_bits;

        // Populate the list.
        valid_bits << "DSPMOD"
                   << "TPECOL"
                   << "CAP"
                   << "BACKBONE"
                   << "DICT"
                   << "ESSENTIAL"
                   << "WATER"
                   << "DIRECT";

        // Extract the SYBYL status bit.
        status_bit = data[9].simplified().toUpper();

        // Check that the status bit isn't empty.
        if (not status_bit.isEmpty())
        {
            // Check that the status bit is valid.
            if (not valid_bits.contains(status_bit))
            {
                errors.append(QObject::tr("Invalid SYBYL atom status bit "
                    "in part (%1) on line '%2'").arg(status_bit).arg(line));
            }
        }
    }
}

/** Generate a Mol2 record from the atom data. */
QString Mol2Atom::toMol2Record() const
{
    return record;
}

/** Generate a string representation of the object. */
QString Mol2Atom::toString() const
{
    return QObject::tr("Mol2Atom::null");
}

const char* Mol2Atom::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Mol2Atom>() );
}

SireMaths::Vector Mol2Atom::getCoord() const
{
    return coord;
}

double Mol2Atom::getCharge() const
{
    return charge;
}

/** Default constructor. */
Mol2Bond::Mol2Bond()
{
}

/** Constructor. */
Mol2Bond::Mol2Bond(const QString &line, QStringList &errors) :
    record(line)
{
    // Tokenize the string, splitting using a single whitespace characters.
    QStringList data = line.simplified().split(QRegExp("\\s"));

    // There must be at least four records.
    if (data.count() < 4)
    {
        errors.append(QObject::tr("The @<TRIPOS>BOND record "
            "is incorrectly formatted. Should have a minimum of 4 entries, has %1!")
            .arg(data.count()));
    }

    bool ok;

    // Extract the ID bond.
    number = data[0].toInt(&ok);

    if (not ok)
    {
        errors.append(QObject::tr("Cannot extract the bond ID "
            "from part (%1) from line '%2'")
            .arg(data[0]).arg(line));

        return;
    }

    // Extract the ID of the origin atom.
    origin = data[1].toInt(&ok);

    if (not ok)
    {
        errors.append(QObject::tr("Cannot extract the origin atom ID "
            "from part (%1) from line '%2'")
            .arg(data[1]).arg(line));

        return;
    }

    // Extract the ID of the target atom.
    target = data[2].toInt(&ok);

    if (not ok)
    {
        errors.append(QObject::tr("Cannot extract the target atom ID "
            "from part (%1) from line '%2'")
            .arg(data[2]).arg(line));

        return;
    }

    // List of valid SYBYL bond types.
    QStringList valid_bonds;

    // Populate the list.
    valid_bonds << "1"
                << "2"
                << "3"
                << "am"
                << "ar"
                << "du"
                << "un"
                << "nc";

    // Extract the bond type.
    type = data[3].toLower();

    if (not valid_bonds.contains(type))
    {
        errors.append(QObject::tr("Invalid SYBYL bond type "
            "in part (%1) on line '%2'").arg(type).arg(line));

        return;
    }

    // The record contains SYBYL status bit information.
    if (data.count() == 5)
    {
        // List of valid SYBYL molecule status bits.
        QStringList valid_bits;

        // Populate the list.
        valid_bits << "TYPECOL"
                   << "GROUP"
                   << "CAP"
                   << "BACKBONE"
                   << "DICT"
                   << "INTERRES";

        // Extract the SYBYL status bit.
        status_bit = data[4].simplified().toUpper();

        // Check that the status bit isn't empty.
        if (not status_bit.isEmpty())
        {
            // Check that the status bit is valid.
            if (not valid_bits.contains(status_bit))
            {
                errors.append(QObject::tr("Invalid SYBYL atom status bit "
                    "in part (%1) on line '%2'").arg(status_bit).arg(line));
            }
        }
    }
}

/** Generate a Mol2 record from the bond data. */
QString Mol2Bond::toMol2Record() const
{
    return record;
}

/** Generate a string representation of the object. */
QString Mol2Bond::toString() const
{
    return QObject::tr("Mol2Bond::null");
}

const char* Mol2Bond::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Mol2Bond>() );
}

/** Get the bond ID. */
qint64 Mol2Bond::getID() const
{
    return number;
}

/** Get the ID of the origin atom. */
qint64 Mol2Bond::getOrigin() const
{
    return origin;
}

/** Get the ID of the target atom. */
qint64 Mol2Bond::getTarget() const
{
    return target;
}

/** Get the bond type. */
QString Mol2Bond::getType() const
{
    return type;
}

/** Default constructor. */
Mol2Molecule::Mol2Molecule() :
    num_atoms(0),
    num_bonds(0),
    num_subst(0),
    num_feats(0),
    num_sets(0)
{
}

/** Constructor.
    @param lines
        A vector of strings containing the molecule record.

    @param errors
        An array of error messages.

    @param
        The number of records found.
 */
Mol2Molecule::Mol2Molecule(const QVector<QString> &lines,
    QStringList &errors, int &num_records) :
    record(lines),
    num_atoms(0),
    num_bonds(0),
    num_subst(0),
    num_feats(0),
    num_sets(0)
{
    // Extract the molecule name.
    name = lines[0].simplified();

    // Tokenize the string, splitting using a single whitespace characters.
    QStringList data = lines[1].simplified().split(QRegExp("\\s"));

    // List of "number" data strings.
    QStringList strings;
    strings << "atoms" << "bonds" << "substructures" << "features" << "sets";

    // Loop over all entries in the data list.
    for (int i=0; i<data.count(); ++i)
    {
        bool ok;

        // Extract the appropriate data record.
        if      (i == 0) num_atoms = data[0].toInt(&ok);
        else if (i == 1) num_bonds = data[1].toInt(&ok);
        else if (i == 2) num_subst = data[2].toInt(&ok);
        else if (i == 3) num_feats = data[3].toInt(&ok);
        else if (i == 4) num_sets  = data[4].toInt(&ok);

        // We currently store the number of SET and FEATURE records, even though
        // these are redundant for the purposes of constructing a molecule.
        // Any file written by the parser will have zero of these record types.
        // TODO: Maybe delete these data members.

        if (not ok)
        {
            errors.append(QObject::tr("Cannot extract the number of %1 "
                "from part (%2) from line '%3'")
                .arg(strings[i]).arg(data[0]).arg(lines[1]));

            return;
        }
    }

    // The Mol2 record indicator. All record types start with this prefix.
    const QString record_indicator("@<TRIPOS>");

    // List of valid molecule types.
    QStringList valid_mols;

    // Populate the list.
    valid_mols << "SMALL"
               << "BIOPOLYMER"
               << "PROTEIN"
               << "NUCLEIC_ACID"
               << "SACCHARIDE";

    // Extract the molecule type.
    mol_type = lines[2].simplified().toUpper();

    // Check that we've not hit another TRIPOS record indicator.
    // This can happen with poorly formatted files.
    if (mol_type.left(9) == record_indicator)
    {
        errors.append(QObject::tr("The @<TRIPOS>MOLECULE record "
            "is incorrectly formatted. Should have 6 lines, has %1!")
            .arg(2));

        num_records = 2;

        return;
    }

    // Check that the molecule type isn't empty.
    if (not mol_type.isEmpty())
    {
        // Check that the status bit is valid.
        if (not valid_mols.contains(mol_type))
        {
            errors.append(QObject::tr("Invalid molecule type "
                "in part (%1) on line '%2'").arg(mol_type).arg(lines[2]));
        }
    }

    // List of valid charge types.
    QStringList valid_chgs;

    // Populate the list.
    valid_chgs << "NO_CHARGES"
               << "DEL_RE"
               << "GASTEIGER"
               << "GAST_HUCK"
               << "HUCKEL"
               << "PULLMAN"
               << "GAUSS80_CHARGES"
               << "AMPAC_CHARGES"
               << "MULLIKEN_CHARGES"
               << "DICT_CHARGES"
               << "MMFF94_CHARGES"
               << "USER_CHARGES";

    // Extract the charge type.
    charge_type = lines[3].simplified().toUpper();

    // Check that we've not hit another TRIPOS record indicator.
    if (charge_type.left(9) == record_indicator)
    {
        errors.append(QObject::tr("The @<TRIPOS>MOLECULE record "
            "is incorrectly formatted. Should have 6 lines, has %1!")
            .arg(3));

        num_records = 3;

        return;
    }

    // Check that the charge type isn't empty.
    if (not charge_type.isEmpty())
    {
        // Check that the status bit is valid.
        if (not valid_chgs.contains(charge_type))
        {
            errors.append(QObject::tr("Invalid charge type "
                "in part (%1) on line '%2'").arg(charge_type).arg(lines[3]));
        }
    }

    // List of valid SYBYL molecule status bits.
    QStringList valid_bits;

    // Populate the list.
    valid_bits << "SYSTEM"
               << "INVALID_CHARGES"
               << "ANALYZED"
               << "SUBSTITUTED"
               << "ALTERED"
               << "REF_ANGLE"
               << "****";

    // Extract the SYBYL status bit.
    status_bit = lines[4].simplified().toUpper();
    bool has_status = true;

    // Check that we've not hit another TRIPOS record indicator.
    if (status_bit.left(9) == record_indicator)
    {
        errors.append(QObject::tr("The @<TRIPOS>MOLECULE record "
            "is incorrectly formatted. Should have 6 lines, has %1!")
            .arg(4));

        num_records = 4;

        return;
    }

    // Check that the status bit isn't empty.
    if (not status_bit.isEmpty())
    {
        // Check that the status bit is valid.
        if (not valid_bits.contains(status_bit))
        {
            errors.append(QObject::tr("Invalid SYBYL atom status bit "
                "in part (%1) on line '%2'").arg(status_bit).arg(lines[4]));

            has_status = false;
        }
    }

    // Extract the comment string.
    comment = lines[5];

    // Check that we've not hit another TRIPOS record indicator.
    if (comment.left(9) == record_indicator)
    {
        errors.append(QObject::tr("The @<TRIPOS>MOLECULE record "
            "is incorrectly formatted. Should have 6 lines, has %1!")
            .arg(5));

        num_records = 5;

        return;
    }

    // Check that the comment isn't empty.
    if (not comment.isEmpty())
    {
        // If we don't have a status bit, assume that the status bit
        // was the comment. This seems to be a common way of formatting
        // molecule records.

        if (not has_status)
        {
            comment = status_bit;
            status_bit = "";
        }
    }

    num_records = 6;
}

/** Generate a Mol2 record from the molecule data. */
QVector<QString> Mol2Molecule::toMol2Record() const
{
    return record;
}

/** Generate a string representation of the object. */
QString Mol2Molecule::toString() const
{
    return QObject::tr("Mol2Molecule::null");
}

const char* Mol2Molecule::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Mol2Molecule>() );
}

/** Get the number of atoms in the molecule. */
int Mol2Molecule::nAtoms() const
{
    return num_atoms;
}

/** Get the number bonds in the molecule. */
int Mol2Molecule::nBonds() const
{
    return num_bonds;
}

/** Get the number of substructures in the molecule. */
int Mol2Molecule::nSubstructures() const
{
    return num_subst;
}

/** Append an atom to the molecule. */
void Mol2Molecule::appendAtom(const Mol2Atom& atom)
{
    atoms.append(atom);
}

/** Append a bond to the molecule. */
void Mol2Molecule::appendBond(const Mol2Bond& bond)
{
    bonds.append(bond);
}

/** Append a substructure to the molecule. */
void Mol2Molecule::appendSubstructure(const Mol2Substructure& substructure)
{
    substructures.append(substructure);
}

/** Default constructor. */
Mol2Substructure::Mol2Substructure() :
    num_inter_bonds(0)
{
}

/** Constructor. */
Mol2Substructure::Mol2Substructure(const QString &line, QStringList &errors) :
    record(line),
    num_inter_bonds(0)
{
    // Tokenize the string, splitting using a single whitespace characters.
    QStringList data = line.simplified().split(QRegExp("\\s"));

    // There must be at least three records.
    if (data.count() < 3)
    {
        errors.append(QObject::tr("The @<TRIPOS>SUBSTRUCTURE record "
            "is incorrectly formatted. Should have a minimum of 3 entries, has %1!")
            .arg(data.count()));
    }

    // Extract the substructure ID number.
    bool ok;
    number = data[0].toInt(&ok);

    if (not ok)
    {
        errors.append(QObject::tr("Cannot extract the substructure ID number "
            "from part (%1) from line '%2'").arg(data[0]).arg(line));

        return;
    }

    // Extract the name of the substructure.
    name = data[1];

    // Extract the ID number of the root atom.
    root_atom = data[2].toInt(&ok);

    if (not ok)
    {
        errors.append(QObject::tr("Cannot extract the root atom ID number "
            "from part (%1) from line '%2'").arg(data[2]).arg(line));

        return;
    }

    // Parse optional records.
    if (data.count() > 3)
    {
        // List of valid substructure types.
        QStringList valid_types;

        // Populate the list.
        valid_types << "TEMP"
                    << "PERM"
                    << "RESIDUE"
                    << "GROUP"
                    << "DOMAIN"
                    << "****";

        // Extract the substructure type.
        type = data[3];

        // Check that the substructure type is valid.
        if (not valid_types.contains(type))
        {
            errors.append(QObject::tr("Invalid substructure type "
                "in part (%1) on line '%2'").arg(type).arg(line));

            return;
        }

        if (data.count() > 4)
        {
            // Extract the dictionary type.
            dict_type = data[4].toInt(&ok);

            if (not ok)
            {
                errors.append(QObject::tr("Cannot extract the dictionary type "
                    "from part (%1) from line '%2'").arg(data[4]).arg(line));

                return;
            }

            if (data.count() > 5)
            {
                // Extract the name of the chain to which the substructure belongs.
                chain = data[5];

                // If this is just a dummy entry, set the chain name to a NULL QString.
                if (chain == "****") chain = QString();

                if (data.count() > 6)
                {
                    // Extract the sub-type of the chain.
                    sub_type = data[6];

                    // If this is just a dummy entry, set the chain name to a NULL QString.
                    if (sub_type == "****") sub_type = QString();

                    if (data.count() > 7)
                    {
                        // Extract the number of inter substructure bonds.
                        num_inter_bonds = data[7].toInt(&ok);

                        if (not ok)
                        {
                            errors.append(QObject::tr("Cannot extract the number of inter "
                                "substructure bonds from part (%1) from line '%2'")
                                .arg(data[7]).arg(line));

                            return;
                        }

                        if (data.count() > 8)
                        {
                            // List of valid SYBYL status bits.
                            QStringList valid_bits;

                            // Populate the list.
                            valid_bits << "LEAF"
                                       << "ROOT"
                                       << "TYPECOL"
                                       << "DICT"
                                       << "BACKWARD"
                                       << "BLOCK"
                                       << "****";

                            // Extract the status bit.
                            status_bit = data[8];

                            // Check that the status bit is valid.
                            if (not valid_bits.contains(status_bit))
                            {
                                errors.append(QObject::tr("Invalid SYBYL status bit "
                                    "in part (%1) on line '%2'").arg(status_bit).arg(line));

                                return;
                            }

                            if (data.count() > 9)
                            {
                                // All other data records form a comment string.

                                // Extract the first string.
                                comment = data[9];

                                // Now concatenate with all other strings.
                                for (int i=10; i<data.count(); ++i)
                                    comment += (" " + data[i]);
                            }
                        }
                    }
                }
            }
        }
    }
}

/** Generate a Mol2 record from the substructure data. */
QString Mol2Substructure::toMol2Record() const
{
    return record;
}

/** Generate a string representation of the object. */
QString Mol2Substructure::toString() const
{
    return QObject::tr("Mol2Substructure::null");
}

const char* Mol2Substructure::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Mol2Substructure>() );
}

/** Constructor */
Mol2::Mol2() : ConcreteProperty<Mol2,MoleculeParser>()
{}

/** Construct to read in the data from the file called 'filename'. The
    passed property map can be used to pass extra parameters to control
    the parsing */
Mol2::Mol2(const QString &filename, const PropertyMap &map)
     : ConcreteProperty<Mol2,MoleculeParser>(filename,map)
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
Mol2::Mol2(const QStringList &lines, const PropertyMap &map)
     : ConcreteProperty<Mol2,MoleculeParser>(lines,map)
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
Mol2::Mol2(const SireSystem::System &system, const PropertyMap &map)
     : ConcreteProperty<Mol2,MoleculeParser>(map)
{
    //look through the system and extract out the information
    //that is needed to go into the file, based on the properties
    //found in 'map'. Do this to create the set of text lines that
    //will make up the file

    QStringList lines;  // = code used to generate the lines

    //now that you have the lines, reparse them back into a Mol2 object,
    //so that the information is consistent, and you have validated that
    //the lines you have written are able to be correctly read in by
    //this parser. This will also implicitly call 'assertSane()'
    Mol2 parsed( lines, map );

    this->operator=(parsed);
}

/** Copy constructor */
Mol2::Mol2(const Mol2 &other)
     : ConcreteProperty<Mol2,MoleculeParser>(other)
{}

/** Destructor */
Mol2::~Mol2()
{}

/** Copy assignment operator */
Mol2& Mol2::operator=(const Mol2 &other)
{
    if (this != &other)
    {
        MoleculeParser::operator=(other);
    }

    return *this;
}

/** Comparison operator */
bool Mol2::operator==(const Mol2 &other) const
{
    return MoleculeParser::operator==(other);
}

/** Comparison operator */
bool Mol2::operator!=(const Mol2 &other) const
{
    return not operator==(other);
}

/** Return the C++ name for this class */
const char* Mol2::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Mol2>() );
}

/** Return the C++ name for this class */
const char* Mol2::what() const
{
    return Mol2::typeName();
}

/** Return the parser that has been constructed by reading in the passed
    file using the passed properties */
MoleculeParserPtr Mol2::construct(const QString &filename,
                                  const PropertyMap &map) const
{
    return Mol2(filename,map);
}

/** Return the parser that has been constructed by reading in the passed
    text lines using the passed properties */
MoleculeParserPtr Mol2::construct(const QStringList &lines,
                                  const PropertyMap &map) const
{
    return Mol2(lines,map);
}

/** Return the parser that has been constructed by extract all necessary
    data from the passed SireSystem::System using the specified properties */
MoleculeParserPtr Mol2::construct(const SireSystem::System &system,
                                  const PropertyMap &map) const
{
    return Mol2(system,map);
}

/** Return a string representation of this parser */
QString Mol2::toString() const
{
    return QObject::tr("Mol2::null");
}

/** Return the format name that is used to identify this file format within Sire */
QString Mol2::formatName() const
{
    return "MOL2";
}

/** Return a description of the file format */
QString Mol2::formatDescription() const
{
    return QObject::tr("Sybyl Mol2 format files.");
}

/** Return the suffixes that these files are normally associated with */
QStringList Mol2::formatSuffix() const
{
    static const QStringList suffixes = { "mol2" };
    return suffixes;
}

/** Return the number of molecules in the system. */
int Mol2::nMols() const
{
    return molecules.count();
}

/** Return the number of atoms in each molecule. */
QVector<int> Mol2::nMolAtoms() const
{
    QVector<int> num_atoms(molecules.count());

    for (int i=0; i<molecules.count(); ++i)
        num_atoms[i] = molecules[i].nAtoms();

    return num_atoms;
}

/** Return the total number of atoms in all molecules. */
int Mol2::nAtoms() const
{
    int num_atoms = 0;

    for (int i=0; i<molecules.count(); ++i)
        num_atoms += molecules[i].nAtoms();

    return num_atoms;
}

/** Function that is called to assert that this object is sane. This
    should raise an exception if the parser is in an invalid state */
void Mol2::assertSane() const
{
    //check state, raise SireError::program_bug if we are in an invalid state
}

/** Internal function that is used to actually parse the data contained
    in the lines of the file */
void Mol2::parseLines(const PropertyMap &map)
{
    /* File format is described here:
         http://chemyang.ccnu.edu.cn/ccb/server/AIMMS/mol2.pdf

       Mol2 is a "free" file format, so records can potentially appear in any order.
       Typically there will be a MOLECULE record, which is followed by ATOM and, e.g.
       BOND records. For simplicity, we'll assume that there is a sensible ordering of
       data, then deal with any inconsistencies afterwards.
     */

    // The Mol2 record indicator. All record types start with this prefix.
    const QString record_indicator("@<TRIPOS>");

    // Molecule index.
    int imol = 0;

    // Internal function used to parse a single ATOM record.
    auto parse_atom = [&](Mol2Atom &local_atom,
        const QString &line, QStringList &local_errors)
    {
        local_atom = Mol2Atom(line, local_errors);
    };

    // Internal function used to parse a single BOND record.
    auto parse_bond = [&](Mol2Bond &local_bond,
        const QString &line, QStringList &local_errors)
    {
        local_bond = Mol2Bond(line, local_errors);
    };

    // Internal function used to parse a single SUBSTRUCTURE record.
    auto parse_subst= [&](Mol2Substructure &local_subst,
        const QString &line, QStringList &local_errors)
    {
        local_subst = Mol2Substructure(line, local_errors);
    };

    // Loop through all lines in the file.
    for (int iline=0; iline<lines().count(); ++iline)
    {
        // Store a reference to the current line.
        const QString &line = lines()[iline];

        // See if this the start of a Mol2 record.
        if (line.left(9) == record_indicator)
        {
            // Extract the record type.
            // Use "simplfied" method to clip trailing whitespace and formatting
            // characters, such as ^M carriage returns.
            QString record_type = line.section('>', 1).simplified();

            // Parse a MOLECULE record.
            if (record_type == "MOLECULE")
            {
                ++iline;

                int num_records = 0;

                // Create a molecule.
                Mol2Molecule mol(lines().mid(iline, iline+5), parse_warnings, num_records);

                // Append a new molecule.
                molecules.append(mol);

                // Fast-forward the line index, accounting for incorrectly formatted entries.
                iline += (num_records - 1);

                // Update the molecule index.
                ++imol;
            }

            // Parse an ATOM record section.
            else if (record_type == "ATOM")
            {
                // For correctly formatted files, the number of atoms should
                // be equal to "num_atoms" from the previous MOLECULE record.

                if (usesParallel())
                {
                    ++iline;

                    QMutex mutex;

                    tbb::parallel_for( tbb::blocked_range<int>(0,molecules[imol-1].nAtoms()),
                                    [&](const tbb::blocked_range<int> &r)
                    {
                        QStringList local_errors;
                        Mol2Atom local_atom;

                        for (int i=r.begin(); i<r.end(); ++i)
                            parse_atom(local_atom, lines()[iline+i], local_errors);

                        QMutexLocker lkr(&mutex);

                        molecules[imol-1].appendAtom(local_atom);
                        parse_warnings  += local_errors;
                    });

                    iline += (molecules[imol-1].nAtoms() - 1);
                }
                else
                {
                    for (int i=0; i<molecules[imol-1].nAtoms(); ++i)
                    {
                        // Create a new atom object.
                        Mol2Atom atom(lines()[++iline], parse_warnings);

                        // Insert the atom into the current molecule.
                        molecules[imol-1].appendAtom(atom);
                    }
                }
            }

            // Parse a BOND record section.
            else if (record_type == "BOND")
            {
                // For correctly formatted files, the number of bonds should
                // be equal to "num_bonds" from the previous MOLECULE record.

                if (usesParallel())
                {
                    ++iline;

                    QMutex mutex;

                    tbb::parallel_for( tbb::blocked_range<int>(0,molecules[imol-1].nBonds()),
                                    [&](const tbb::blocked_range<int> &r)
                    {
                        QStringList local_errors;
                        Mol2Bond local_bond;

                        for (int i=r.begin(); i<r.end(); ++i)
                            parse_bond(local_bond, lines()[iline+i], local_errors);

                        QMutexLocker lkr(&mutex);

                        molecules[imol-1].appendBond(local_bond);
                        parse_warnings  += local_errors;
                    });

                    iline += (molecules[imol-1].nBonds() - 1);
                }
                else
                {
                    for (int i=0; i<molecules[imol-1].nBonds(); ++i)
                    {
                        // Create a new bond object.
                        Mol2Bond bond(lines()[++iline], parse_warnings);

                        // Insert the bond into the current molecule.
                        molecules[imol-1].appendBond(bond);
                    }
                }
            }

            // Parse a SUBSTRUCTURE record.
            else if (record_type == "SUBSTRUCTURE")
            {
                // For correctly formatted files, the number of substructures should
                // be equal to "num_subst" from the previous MOLECULE record.

                if (usesParallel())
                {
                    ++iline;

                    QMutex mutex;

                    tbb::parallel_for( tbb::blocked_range<int>(0,molecules[imol-1].nSubstructures()),
                                    [&](const tbb::blocked_range<int> &r)
                    {
                        QStringList local_errors;
                        Mol2Substructure local_subst;

                        for (int i=r.begin(); i<r.end(); ++i)
                            parse_subst(local_subst, lines()[iline+i], local_errors);

                        QMutexLocker lkr(&mutex);

                        molecules[imol-1].appendSubstructure(local_subst);
                        parse_warnings  += local_errors;
                    });

                    iline += (molecules[imol-1].nSubstructures() - 1);
                }
                else
                {
                    for (int i=0; i<molecules[imol-1].nSubstructures(); ++i)
                    {
                        // Create a new substructure object.
                        Mol2Substructure substructure(lines()[++iline], parse_warnings);

                        // Insert the substructure into the current molecule.
                        molecules[imol-1].appendSubstructure(substructure);
                    }
                }
            }
        }
    }

    this->setScore(nAtoms());
}

/** Use the data contained in this parser to create a new System of molecules,
    assigning properties based on the mapping in 'map' */
System Mol2::startSystem(const PropertyMap &map) const
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
void Mol2::addToSystem(System &system, const PropertyMap &map) const
{
    //you should loop through each molecule in the system and work out
    //which ones are described in the file, and then add data from the file
    //to thise molecules.
}
