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

#include "SireIO/mol2.h"

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

using namespace SireBase;
using namespace SireIO;
using namespace SireMol;
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
        << mol2atom.charge << mol2atom.status_bits;

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
            >> mol2atom.charge >> mol2atom.status_bits;
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
        << mol2bond.status_bits;

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
            >> mol2bond.status_bits;
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
        << mol2molecule.status_bits << mol2molecule.comment << mol2molecule.atoms
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
            >> mol2molecule.status_bits >> mol2molecule.comment >> mol2molecule.atoms
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
        << mol2subst.num_inter_bonds << mol2subst.status_bits << mol2subst.comment;

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
            >> mol2subst.num_inter_bonds >> mol2subst.status_bits >> mol2subst.comment;
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
    type("Du"),
    charge(0)
{
}

/** Constructor. */
Mol2Atom::Mol2Atom(const QString &line, QStringList &errors) :
    record(line),
    name("X"),
    type("Du"),
    charge(0)
{
    // Tokenize the string, splitting using a single whitespace character.
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

    // Parse optional data records.
    if (data.count() > 6)
    {
        // Extract the substructure ID.
        subst_id = data[6].toInt(&ok);

        if (not ok)
        {
            errors.append(QObject::tr("Cannot extract the substructure ID number "
                "from part (%1) from line '%2'").arg(data[6]).arg(line));

            return;
        }

        if (data.count() > 7)
        {
            // Extract the substructure name.
            subst_name = data[7].simplified();

            if (data.count() > 8)
            {
                // Extract the charge on the atom.
                charge = data[8].toDouble(&ok);

                if (not ok)
                {
                    errors.append(QObject::tr("Cannot extract the atom charge "
                        "from part (%1) from line '%2'").arg(data[8]).arg(line));

                    return;
                }

                if (data.count() > 9)
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
                               << "DIRECT"
                               << "****";

                    // Extract the SYBYL status bit.
                    status_bits = data[9].simplified().toUpper();

                    // Check that the status bit isn't empty.
                    if (not status_bits.isEmpty())
                    {
                        // Tokenize the bits.
                        QStringList split_bits = status_bits.split('|');

                        // Check that the status bit is valid.
                        for (auto const &bit : split_bits)
                        {
                            // Check that the status bit is valid.
                            if (not valid_bits.contains(bit))
                            {
                                errors.append(QObject::tr("Invalid SYBYL substructure status bit "
                                    "in part (%1) on line '%2'").arg(bit).arg(line));
                            }
                        }
                    }
                }
            }
        }
    }
}

/** Constructor.
    @param atom
        A reference to a Sire Atom object.

    @param map
        A reference to the user parameter map.

    @param errors
        An array of error messages.

    @param is_idx
        Whether to number residues by their index. (optional)
 */
Mol2Atom::Mol2Atom(const SireMol::Atom &atom, const PropertyMap &map,
    QStringList &errors, bool is_idx) :
    number(atom.number().value()),
    name(atom.name().value()),
    type("Du"),
    charge(0)
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
        // Set the substructure name.
        subst_name = atom.residue().name().value();

        // Set the substructure number.
        if (is_idx) subst_id = atom.residue().index().value() + 1;
        else        subst_id = atom.residue().number().value();
    }

    // Extract the SYBYL atom type.
    if (atom.hasProperty(map["sybyl-atom-type"]))
    {
        type = atom.property<QString>(map["sybyl-atom-type"]);
    }
    else
    {
        // TODO: Some way of inferring the type...
    }

    // Extract the atomic charge.
    if (atom.hasProperty(map["charge"]))
    {
        charge = atom.property<SireUnits::Dimension::Charge>(map["charge"]).value();
    }
    else if (atom.hasProperty(map["formal-charge"]))
    {
        // TODO: Need some conversion here, I assume.
        charge = atom.property<SireUnits::Dimension::Charge>(map["formal-charge"]).value();
    }

    // Extract the SYBYL status bits.
    if (atom.hasProperty(map["atom-status-bits"]))
    {
        status_bits = atom.property<QString>(map["atom-status-bits"]);
    }
}

/** Generate a Mol2 record from the atom data. */
QString Mol2Atom::toMol2Record() const
{
    QString line = QString("%1 %2 %3 %4 %5 %6 %7 %8 %9")
        .arg(number, 7)
        .arg(name, -10)
        .arg(coord[0], 9, 'f', 4)
        .arg(coord[1], 9, 'f', 4)
        .arg(coord[2], 9, 'f', 4)
        .arg(type, -5)
        .arg(subst_id, 4)
        .arg(subst_name, -9)
        .arg(charge, 7, 'f', 4);

    if (not status_bits.isEmpty())
        line.append(QString(" %1").arg(status_bits));

    return line;
}

const char* Mol2Atom::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Mol2Atom>() );
}

/** Get the atom number. */
int Mol2Atom::getNumber() const
{
    return number;
}

/** Get the atom name. */
QString Mol2Atom::getName() const
{
    return name;
}

/** Get the atom coordinates. */
SireMaths::Vector Mol2Atom::getCoord() const
{
    return coord;
}

/** Get the SYBYL atom type. */
QString Mol2Atom::getType() const
{
    return type;
}

/** Get the number of the substructure containing the atom. */
qint64 Mol2Atom::getSubstructureNumber() const
{
    return subst_id;
}

/** Get the name of the substructure containing the atom. */
QString Mol2Atom::getSubstructureName() const
{
    return subst_name;
}

/** Get the atom charge. */
double Mol2Atom::getCharge() const
{
    return charge;
}

/** Get the SYBYL status bit. */
QString Mol2Atom::getStatusBits() const
{
    return status_bits;
}

/** Default constructor. */
Mol2Bond::Mol2Bond()
{
}

/** Constructor. */
Mol2Bond::Mol2Bond(const QString &line, QStringList &errors) :
    record(line)
{
    // Tokenize the string, splitting using a single whitespace character.
    QStringList data = line.simplified().split(QRegExp("\\s"));

    // There must be at least four records.
    if (data.count() < 4)
    {
        errors.append(QObject::tr("The @<TRIPOS>BOND record "
            "is incorrectly formatted. Should have a minimum of 4 entries, has %1!")
            .arg(data.count()));

        return;
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
                   << "INTERRES"
                   << "****";

        // Extract the SYBYL status bit.
        status_bits = data[4].simplified().toUpper();

        // Check that the status bit isn't empty.
        if (not status_bits.isEmpty())
        {
            // Tokenize the bits.
            QStringList split_bits = status_bits.split('|');

            // Check that the status bit is valid.
            for (auto const &bit : split_bits)
            {
                // Check that the status bit is valid.
                if (not valid_bits.contains(bit))
                {
                    errors.append(QObject::tr("Invalid SYBYL bond status bit "
                        "in part (%1) on line '%2'").arg(bit).arg(line));
                }
            }
        }
    }
}

/** Generate a Mol2 record from the bond data. */
QString Mol2Bond::toMol2Record() const
{
    return record;
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
    num_sets(0),
    mol_type("SMALL"),
    charge_type("GASTEIGER")
{
}

/** Constructor.
    @param lines
        A vector of strings containing the molecule record.

    @param errors
        An array of error messages.

    @param
        The number of records found.

    @param
        The name of the Mol2 file containing the data record (optional).

    @param
        The molecule index (optional).
 */
Mol2Molecule::Mol2Molecule(const QVector<QString> &lines,
    QStringList &errors, int &num_records, int imol) :
    record(lines),
    num_atoms(0),
    num_bonds(0),
    num_subst(0),
    num_feats(0),
    num_sets(0),
    mol_type("SMALL"),
    charge_type("GASTEIGER")
{
    // Extract the molecule name.
    name = lines[0].simplified();

    // If the name is blank, then name the molecule "Molecule".
    if (name.isEmpty())
    {
        name = QString("Molecule: ");

        // Append a molecule index (for multiple molecule records).
        if (imol != -1)
        {
            // Ensire indexing starts from one.
            if (imol == 0) imol = 1;

            name.append(QString(" %1").arg(imol));
        }
    }

    // Tokenize the string, splitting using a single whitespace character.
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
            "is incorrectly formatted. Should have at least 4 lines, has %1!")
            .arg(2));

        num_records = 2;

        return;
    }

    // Check that the status bit is valid.
    if (not valid_mols.contains(mol_type))
    {
        errors.append(QObject::tr("Invalid molecule type "
            "in part (%1) on line '%2'").arg(mol_type).arg(lines[2]));
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
            "is incorrectly formatted. Should have at least 4 lines, has %1!")
            .arg(3));

        num_records = 3;

        return;
    }

    // Check that the status bit is valid.
    if (not valid_chgs.contains(charge_type))
    {
        errors.append(QObject::tr("Invalid charge type "
            "in part (%1) on line '%2'").arg(charge_type).arg(lines[3]));
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
    status_bits = lines[4].simplified().toUpper();
    bool has_status = true;

    // Check that we've not hit another TRIPOS record indicator.
    if (status_bits.left(9) == record_indicator)
    {
        status_bits.clear();
        num_records = 4;

        return;
    }

    // Check that the status bit isn't empty.
    if (not status_bits.isEmpty())
    {
        // Tokenize the bits.
        QStringList split_bits = status_bits.split('|');

        // Check that the status bit is valid.
        for (auto const &bit : split_bits)
        {
            // Check that the status bit is valid.
            if (not valid_bits.contains(bit))
            {
                errors.append(QObject::tr("Invalid SYBYL molecule status bit "
                    "in part (%1) on line '%2'").arg(bit).arg(lines[4]));

                has_status = false;
            }
        }
    }

    // Extract the comment string.
    comment = lines[5];

    // Check that we've not hit another TRIPOS record indicator.
    if (comment.left(9) == record_indicator)
    {
        comment.clear();
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
            comment = status_bits;
            status_bits = "";
        }
    }

    num_records = 6;
}

/** Constructor.
    @param mol
        A reference to a Sire Molecule object.

    @param map
        A reference to the user parameter map.

    @param errors
        An array of error messages.

    @param imol
        The molecule index.
 */
Mol2Molecule::Mol2Molecule(const SireMol::Molecule &mol, const PropertyMap &map,
    QStringList &errors, int imol) :
    name(mol.name().value()),
    num_atoms(mol.nAtoms()),
    num_bonds(0),
    num_subst(mol.nResidues()),
    num_feats(0),
    num_sets(0),
    mol_type("SMALL"),
    charge_type("GASTEIGER")
{
    // If the name is blank, then name the molecule "Molecule".
    if (name.isEmpty())
    {
        name = QString("Molecule");

        // Append a molecule index (for multiple molecule records).
        if (imol != -1)
        {
            // Ensire indexing starts from one.
            if (imol == 0) imol = 1;

            name.append(QString(": %1").arg(imol));
        }
    }

    // Extract the molecule type.
    if (mol.hasProperty(map["mol-type"]))
    {
        // List of valid molecule types.
        QStringList valid_mols;

        // Populate the list.
        valid_mols << "SMALL"
                   << "BIOPOLYMER"
                   << "PROTEIN"
                   << "NUCLEIC_ACID"
                   << "SACCHARIDE";

        // Extract the molecule type.
        mol_type = mol.property(map["mol-type"]).toString().simplified().toUpper();

        // Check that the status bit is valid.
        if (not valid_mols.contains(mol_type))
        {
            errors.append(QObject::tr("Invalid molecule type: %1")
                .arg(mol_type));

            return;
        }
    }
    else
    {
        // TODO: Some way of inferring the type...
    }

    // Extract the charge type.
    if (mol.hasProperty(map["charge-type"]))
    {
        // List of valid molecule types.
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

        // Extract the molecule type.
        charge_type = mol.property(map["charge-type"]).toString().simplified().toUpper();

        // Check that the status bit is valid.
        if (not valid_chgs.contains(charge_type))
        {
            errors.append(QObject::tr("Invalid charge type: %1")
                .arg(charge_type));

            return;
        }
    }
    else
    {
        // TODO: Some way of inferring the charge...
    }

    // Extract the status bits.
    if (mol.hasProperty(map["mol-status-bits"]))
    {
        status_bits = mol.property(map["mol-status-bits"]).toString().simplified().toUpper();

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
                   << "DIRECT"
                   << "****";

        // Check that the status bit isn't empty.
        if (not status_bits.isEmpty())
        {
            // Tokenize the bits.
            QStringList split_bits = status_bits.split('|');

            // Check that the status bit is valid.
            for (auto const &bit : split_bits)
            {
                // Check that the status bit is valid.
                if (not valid_bits.contains(bit))
                {
                    errors.append(QObject::tr("Invalid SYBYL substructure "
                        "status bit: %1").arg(bit));
                }
            }
        }
    }

    // Extract the comment.
    if (mol.hasProperty(map["mol-comment"]))
    {
        comment = mol.property("mol-comment").toString();
    }
}

/** Generate a Mol2 record from the molecule data. */
QVector<QString> Mol2Molecule::toMol2Record() const
{
    QVector<QString> lines;

    lines.append(name);

    QString line = QString("%1 %2 %3 %4 %5")
                   .arg(nAtoms(), 5)
                   .arg(nBonds(), 5)
                   .arg(nSubstructures(), 5)
                   .arg(0, 5)
                   .arg(0, 5);

    lines.append(line);
    lines.append(mol_type);
    lines.append(charge_type);
    lines.append(status_bits);
    lines.append(comment);

    return lines;
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
void Mol2Molecule::appendAtom(const Mol2Atom &atom)
{
    atoms.append(atom);
}

/** Append a vector of atoms to the molecule. */
void Mol2Molecule::appendAtoms(const QVector<Mol2Atom> &atoms)
{
    for (int i=0; i<atoms.count(); ++i)
        this->atoms.append(atoms[i]);
}

/** Append a bond to the molecule. */
void Mol2Molecule::appendBond(const Mol2Bond &bond)
{
    bonds.append(bond);
}

/** Append a vector of bonds to the molecule. */
void Mol2Molecule::appendBonds(const QVector<Mol2Bond> &bonds)
{
    for (int i=0; i<bonds.count(); ++i)
        this->bonds.append(bonds[i]);
}

/** Append a substructure to the molecule. */
void Mol2Molecule::appendSubstructure(const Mol2Substructure &substructure)
{
    substructures.append(substructure);
}

/** Append a vector of substructures to the molecule. */
void Mol2Molecule::appendSubstructures(const QVector<Mol2Substructure> &substructures)
{
    for (int i=0; i<substructures.count(); ++i)
        this->substructures.append(substructures[i]);
}

/** Get the name of the molecule. */
QString Mol2Molecule::getName() const
{
    return name;
}

/** Get the molecule type. */
QString Mol2Molecule::getMolType() const
{
    return mol_type;
}

/** Get the charge type. */
QString Mol2Molecule::getChargeType() const
{
    return charge_type;
}

/** Get the status bits. */
QString Mol2Molecule::getStatusBits() const
{
    return status_bits;
}

/** Get the comments. */
QString Mol2Molecule::getComment() const
{
    return comment;
}

/** Get the atoms from the molecule. */
QVector<Mol2Atom> Mol2Molecule::getAtoms() const
{
    return atoms;
}

/** Get a specific atom from the molecule. */
Mol2Atom Mol2Molecule::getAtom(int i) const
{
    return atoms[i];
}

/** Get the bonds from the molecule. */
QVector<Mol2Bond> Mol2Molecule::getBonds() const
{
    return bonds;
}

/** Get a specific bond from the molecule. */
Mol2Bond Mol2Molecule::getBond(int i) const
{
    return bonds[i];
}

/** Get the substructures from the molecule. */
QVector<Mol2Substructure> Mol2Molecule::getSubstructures() const
{
    return substructures;
}

/** Get a specific substructure from the molecule. */
Mol2Substructure Mol2Molecule::getSubstructure(int i) const
{
    return substructures[i];
}

/** Default constructor. */
Mol2Substructure::Mol2Substructure() :
    root_atom(0),
    dict_type(0),
    num_inter_bonds(0)
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
Mol2Substructure::Mol2Substructure(const QString &line, QStringList &errors) :
    record(line),
    root_atom(0),
    dict_type(0),
    num_inter_bonds(0)
{
    // Tokenize the string, splitting using a single whitespace character.
    QStringList data = line.simplified().split(QRegExp("\\s"));

    // There must be at least three records.
    if (data.count() < 3)
    {
        errors.append(QObject::tr("The @<TRIPOS>SUBSTRUCTURE record "
            "is incorrectly formatted. Should have a minimum of 3 entries, has %1!")
            .arg(data.count()));

        return;
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
                            status_bits = data[8];

                            // Tokenize the bits.
                            QStringList split_bits = status_bits.split('|');

                            // Check that the status bit is valid.
                            for (auto const &bit : split_bits)
                            {
                                if (not valid_bits.contains(bit))
                                {
                                    errors.append(QObject::tr("Invalid SYBYL status bit "
                                        "in part (%1) on line '%2'").arg(bit).arg(line));

                                    return;
                                }
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

/** Constructor.
    @param res
        A reference to a Sire Residue object.

    @param map
        A reference to the user parameter map.

    @param errors
        An array of error messages.

    @param is_idx
        Whether to number residues by their index. (optional)
 */
Mol2Substructure::Mol2Substructure(const SireMol::Residue &res,
    const PropertyMap &map, QStringList &errors,
    bool is_idx) :
    name(res.name().value()),
    type("RESIDUE"),
    dict_type(0),
    num_inter_bonds(0)
{
    // Set the residue number.
    if (is_idx) number = res.index().value() + 1;
    else        number = res.number().value();

    // Set the chain.
    if (res.isWithinChain())
        chain = res.chain().name().value();

    // Extract the atoms from the residue.
    auto atoms = res.atoms();

    // Set the root atom of the residue.
    root_atom = atoms[0].read().asA<SireMol::Atom>().number().value();

    // Extract some optional data.

    // Substructure type.
    if (res.hasProperty(map["res-type"]))
        type = res.property<QString>(map["res-type"]);

    // Dictionary type.
    if (res.hasProperty(map["res-dict-type"]))
        dict_type = res.property<qint64>(map["res-dict-type"]);

    // Chain sub-type.
    if (res.hasProperty(map["chain-sub-type"]))
        sub_type = res.property<QString>(map["chain-sub-type"]);

    // Internal substructure bonds.
    if (res.hasProperty(map["res-inter-bonds"]))
        num_inter_bonds = res.property<qint64>(map["res-inter-bonds"]);

    // Status bits.
    if (res.hasProperty(map["res-status-bits"]))
        status_bits = res.property<QString>(map["res-status-bits"]);

    // Comments.
    if (res.hasProperty(map["res-comment"]))
        comment = res.property<QString>(map["res-comment"]);
}

/** Generate a Mol2 record from the substructure data. */
QString Mol2Substructure::toMol2Record() const
{
    // The first part of the record is mandatory.
    QString line = QString("%1 %2 %3")
        .arg(number, 6)
        .arg(name, -6)
        .arg(root_atom, 6);

    // The optional part of the line.
    QString optional;

    // Work backwards through the records and add each one.
    // If an entry down the line has been added, then we need to enter "****"
    // if the current string record is empty, or zero for an unspecified integer record.

    if (not comment.isEmpty())
    {
        optional.prepend(QString(" %1").arg(comment));
    }

    if (not status_bits.isEmpty())
    {
        optional.prepend(QString(" %1").arg(status_bits));
    }
    else
    {
        if (not optional.isEmpty())
        {
            optional.prepend(QString(" ****"));
        }
    }

    if (num_inter_bonds != 0)
    {
        optional.prepend(QString(" %1").arg(num_inter_bonds, 6));
    }
    else
    {
        if (not optional.isEmpty())
        {
            optional.prepend(QString("      0"));
        }
    }

    if (not sub_type.isEmpty())
    {
        optional.prepend(QString(" %1").arg(sub_type, -6));
    }
    else
    {
        if (not optional.isEmpty())
        {
            optional.prepend(QString(" ****  "));
        }
    }

    if (not chain.isEmpty())
    {
        optional.prepend(QString(" %1").arg(chain, -5));
    }
    else
    {
        if (not optional.isEmpty())
        {
            optional.prepend(QString(" **** "));
        }
    }

    if (dict_type != 0)
    {
        optional.prepend(QString(" %1").arg(dict_type, 3));
    }
    else
    {
        if (not optional.isEmpty())
        {
            optional.prepend(QString("   0"));
        }
    }

    if (not type.isEmpty())
    {
        optional.prepend(QString(" %1").arg(type, -6));
    }
    else
    {
        if (not optional.isEmpty())
        {
            optional.prepend(QString(" ****  "));
        }
    }

    line.append(QString(" %1").arg(optional));

    return line;
}

const char* Mol2Substructure::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Mol2Substructure>() );
}

/** Get the substructure number. */
qint64 Mol2Substructure::getNumber() const
{
    return number;
}

/** Get the substructure name. */
QString Mol2Substructure::getName() const
{
    return name;
}

/** Get the substructure type. */
QString Mol2Substructure::getType() const
{
    return type;
}

/** Get the dictionary type. */
qint64 Mol2Substructure::getDictType() const
{
    return dict_type;
}

/** Get the substructure chain. */
QString Mol2Substructure::getChain() const
{
    return chain;
}

/** Get the chain sub-type. */
QString Mol2Substructure::getChainSubType() const
{
    return sub_type;
}

/** Get the number of inter-substructure bonds. */
qint64 Mol2Substructure::getInterBonds() const
{
    return num_inter_bonds;
}

/** Get the status bits. */
QString Mol2Substructure::getStatusBits() const
{
    return status_bits;
}

/** Get the comment. */
QString Mol2Substructure::getComment() const
{
    return comment;
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
    // Get the MolNums of each molecule in the System - this returns the
    // numbers in MolIdx order.
    const QVector<MolNum> molnums = system.getMoleculeNumbers().toVector();

    // Store the number of molecules.
    const int nmols = molnums.count();

    // No molecules in the system.
    if (nmols == 0)
    {
        this->operator=(Mol2());
        return;
    }

    // The list of lines.
    QStringList lines;

    // Lines for different Mol2 data records (one for each molecule).
    QVector<QVector<QString> > molecule_lines(nmols);
    QVector<QVector<QString> > atom_lines(nmols);
    QVector<QVector<QString> > substructure_lines(nmols);

    // Resize the molecules vector.
    molecules.resize(nmols);

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
                // Parse the SireMolecule data into a Mol2Molecule.
                molecules[i] = Mol2Molecule(system[molnums[i]].molecule(),
                    map, local_errors, i);

                // Now parse the rest of the molecular data, i.e. atoms, residues, etc.
                parseMolecule(molecules[i], system[molnums[i]].molecule(), i,
                    atom_lines[i], substructure_lines[i], local_errors, map);

                // Generate the Mol2 data record lines.
                molecule_lines[i] = molecules[i].toMol2Record();
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
            // Parse the SireMolecule data into a Mol2Molecule.
            molecules[i] = Mol2Molecule(system[molnums[i]].molecule(), map, parse_warnings);

            // Now parse the rest of the molecular data, i.e. atoms, residues, etc.
            parseMolecule(molecules[i], system[molnums[i]].molecule(), i, atom_lines[i],
                substructure_lines[i], parse_warnings, map);

            // Generate the Mol2 data record lines.
            molecule_lines[i] = molecules[i].toMol2Record();
        }
    }

    // Now assemble the lines from the record data for each molecule.
    // We do this in serial since the order matters.
    for (int i=0; i<nmols; ++i)
    {
        // Molecule records.
        lines.append("@<TRIPOS>MOLECULE");
        lines.append(molecule_lines[i].toList());

        // Atom records.
        lines.append("@<TRIPOS>ATOM");
        lines.append(atom_lines[i].toList());

        // Substructure records.
        lines.append("@<TRIPOS>SUBSTRUCTURE");
        lines.append(substructure_lines[i].toList());
    }

    // Reparse the lines as a self-consistency check.
    Mol2 parsed(lines, map);

    this->operator=(parsed);
}

/** Copy constructor */
Mol2::Mol2(const Mol2 &other) :
    ConcreteProperty<Mol2,MoleculeParser>(other),
    molecules(other.molecules),
    parse_warnings(other.parse_warnings)
{}

/** Destructor */
Mol2::~Mol2()
{}

/** Copy assignment operator */
Mol2& Mol2::operator=(const Mol2 &other)
{
    if (this != &other)
    {
        molecules = other.molecules;
        parse_warnings = other.parse_warnings;

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
    if (lines().isEmpty())
        return QObject::tr("Mol2::null");
    else
    {
        return QObject::tr("Mol2( nMolecules() = %1, "
            "nResidues() = %2, nAtoms() = %3 )")
            .arg(nMolecules()).arg(nSubstructures()).arg(nAtoms());
    }
}

/** Convert the the parsed data to a collection of Mol2 record lines. */
QVector<QString> Mol2::toLines() const
{
    // Store the number of molecules.
    const int num_mols = nMolecules();

    // No molecules in the system.
    if (num_mols == 0)
        return QVector<QString>();

    // The vector of Mol2 record lines.
    QVector<QString> lines;

    // Loop over all molecules.
    for (int i=0; i<num_mols; ++i)
    {
        const int num_atoms = nAtoms(i);
        const int num_subst = nSubstructures(i);

        // Generate the Mol2 moleucle data record lines.
        QVector<QString> molecule_lines = molecules[i].toMol2Record();

        // Data record lines for the molecule.
        QVector<QString> atom_lines(num_atoms);
        QVector<QString> substructure_lines(num_subst);

        if (usesParallel())
        {
            tbb::parallel_for( tbb::blocked_range<int>(0, num_atoms),
                            [&](const tbb::blocked_range<int> r)
            {
                for (int j=r.begin(); j<r.end(); ++j)
                {
                    // Generate the Mol2 atom record lines.
                    atom_lines[j] = molecules[i].getAtom(j).toMol2Record();
                }
            });

            tbb::parallel_for( tbb::blocked_range<int>(0, num_subst),
                            [&](const tbb::blocked_range<int> r)
            {
                for (int j=r.begin(); j<r.end(); ++j)
                {
                    // Generate the Mol2 substructure record lines.
                    substructure_lines[j] = molecules[i].getSubstructure(j).toMol2Record();
                }
            });

        }
        else
        {
            for (int j=0; j<num_atoms; ++j)
            {
                // Generate the Mol2 atom record lines.
                atom_lines[j] = molecules[i].getAtom(j).toMol2Record();
            }

            for (int j=0; j<num_subst; ++j)
            {
                // Generate the Mol2 substructure record lines.
                substructure_lines[j] = molecules[i].getSubstructure(j).toMol2Record();
            }
        }

        // Molecule records.
        lines.append("@<TRIPOS>MOLECULE");
        lines += molecule_lines;

        // Atom records.
        lines.append("@<TRIPOS>ATOM");
        lines += atom_lines;

        // Substructure records.
        lines.append("@<TRIPOS>SUBSTRUCTURE");
        lines += substructure_lines;
    }

    return lines;
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
int Mol2::nMolecules() const
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

/** Return the number of atoms in a specific molecule. */
int Mol2::nAtoms(int i) const
{
    return molecules[i].nAtoms();
}

/** Return the total number of atoms in all molecules. */
int Mol2::nAtoms() const
{
    int num_atoms = 0;

    for (int i=0; i<molecules.count(); ++i)
        num_atoms += molecules[i].nAtoms();

    return num_atoms;
}

/** Return the number of bonds in each molecule. */
QVector<int> Mol2::nMolBonds() const
{
    QVector<int> num_bonds(molecules.count());

    for (int i=0; i<molecules.count(); ++i)
        num_bonds[i] = molecules[i].nBonds();

    return num_bonds;
}

/** Return the number of bonds in a specific molecule. */
int Mol2::nBonds(int i) const
{
    return molecules[i].nBonds();
}

/** Return the total number of bonds in all molecules. */
int Mol2::nBonds() const
{
    int num_bonds = 0;

    for (int i=0; i<molecules.count(); ++i)
        num_bonds += molecules[i].nBonds();

    return num_bonds;
}

/** Return the number of substructures in each molecule. */
QVector<int> Mol2::nMolSubstructures() const
{
    QVector<int> num_subst(molecules.count());

    for (int i=0; i<molecules.count(); ++i)
        num_subst[i] = molecules[i].nSubstructures();

    return num_subst;
}

/** Return the number of substructures in a specific molecule. */
int Mol2::nSubstructures(int i) const
{
    return molecules[i].nSubstructures();
}

/** Return the total number of substructures in all molecules. */
int Mol2::nSubstructures() const
{
    int num_subst = 0;

    for (int i=0; i<molecules.count(); ++i)
        num_subst += molecules[i].nSubstructures();

    return num_subst;
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

    // Number of lines in the file.
    int num_lines = lines().count();

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

                // Check that the file contains enough lines for the record.
                if ((iline == num_lines) or ((iline + 4) >= num_lines))
                {
                    parse_warnings.append(QObject::tr("We've unexpectedly hit the end "
                        "of the file parsing a MOLECULE record!"));

                    return;
                }

                int num_records = 0;

                // Create a molecule.
                Mol2Molecule mol(lines().mid(iline, iline+5),
                    parse_warnings, num_records, ++imol);

                // Append a new molecule.
                molecules.append(mol);

                // Fast-forward the line index, accounting for incorrectly formatted entries.
                iline += (num_records - 1);
            }

            // Parse an ATOM record section.
            else if (record_type == "ATOM")
            {
                // For correctly formatted files, the number of atoms should
                // be equal to "num_atoms" from the previous MOLECULE record.
                const int num_atoms = molecules[imol-1].nAtoms();

                // Check that the file contains enough lines for the record.
                if (((iline + 1) == num_lines) or
                   ((iline + num_atoms) >= num_lines))
                {
                    parse_warnings.append(QObject::tr("We've unexpectedly hit the end "
                        "of the file parsing an ATOM record!"));

                    return;
                }

                if (usesParallel())
                {
                    ++iline;

                    QMutex mutex;

                    // Local data storage for atom objects.
                    QVector<Mol2Atom> local_atoms(num_atoms);

                    tbb::parallel_for( tbb::blocked_range<int>(0, num_atoms),
                                    [&](const tbb::blocked_range<int> &r)
                    {
                        // Create local data objects.
                        QStringList local_errors;

                        for (int i=r.begin(); i<r.end(); ++i)
                        {
                            // Parse the data from the atom record.
                            local_atoms[i] = Mol2Atom(lines()[iline+i], local_errors);
                        }

                        if (not local_errors.isEmpty())
                        {
                            // Acquire a lock.
                            QMutexLocker lkr(&mutex);

                            // Update the warning messages.
                            parse_warnings += local_errors;
                        }
                    });

                    // Append the atoms to the molecule.
                    molecules[imol-1].appendAtoms(local_atoms);

                    // Fast-forward the line index.
                    iline += (num_atoms - 1);
                }
                else
                {
                    for (int i=0; i<num_atoms; ++i)
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
                const int num_bonds = molecules[imol-1].nBonds();

                // Check that the file contains enough lines for the record.
                if (((iline + 1) == num_lines) or
                   ((iline + num_bonds) >= num_lines))
                {
                    parse_warnings.append(QObject::tr("We've unexpectedly hit the end "
                        "of the file parsing a BOND record!"));

                    return;
                }

                if (usesParallel())
                {
                    ++iline;

                    QMutex mutex;

                    // Local data storage for bond objects.
                    QVector<Mol2Bond> local_bonds(num_bonds);

                    tbb::parallel_for( tbb::blocked_range<int>(0, num_bonds),
                                    [&](const tbb::blocked_range<int> &r)
                    {
                        // Create local data objects.
                        QStringList local_errors;

                        for (int i=r.begin(); i<r.end(); ++i)
                        {
                            // Parse the data from the bond record.
                            local_bonds[i] = Mol2Bond(lines()[iline+i], local_errors);
                        }

                        if (not local_errors.isEmpty())
                        {
                            // Acquire a lock.
                            QMutexLocker lkr(&mutex);

                            // Update the warning messages.
                            parse_warnings += local_errors;
                        }
                    });

                    // Append the bonds to the molecule.
                    molecules[imol-1].appendBonds(local_bonds);

                    // Fast-forward the line index.
                    iline += (num_bonds - 1);
                }
                else
                {
                    for (int i=0; i<num_bonds; ++i)
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
                const int num_subst = molecules[imol-1].nSubstructures();

                // Check that the file contains enough lines for the record.
                if (((iline + 1) == num_lines) or
                   ((iline + num_subst) >= num_lines))
                {
                    parse_warnings.append(QObject::tr("We've unexpectedly hit the end "
                        "of the file parsing a SUBSTRUCTURE record!"));

                    return;
                }

                if (usesParallel())
                {
                    ++iline;

                    QMutex mutex;

                    // Local data storage for substructure objects.
                    QVector<Mol2Substructure> local_subst(num_subst);

                    tbb::parallel_for( tbb::blocked_range<int>(0, num_subst),
                                    [&](const tbb::blocked_range<int> &r)
                    {
                        // Create local data objects.
                        QStringList local_errors;

                        for (int i=r.begin(); i<r.end(); ++i)
                        {
                            // Parse the data from the substructure record.
                            local_subst[i] = Mol2Substructure(lines()[iline+i], local_errors);
                        }

                        if (not local_errors.isEmpty())
                        {
                            // Acquire a lock.
                            QMutexLocker lkr(&mutex);

                            // Update the warning messages.
                            parse_warnings += local_errors;
                        }
                    });

                    // Append the substructures to the molecule.
                    molecules[imol-1].appendSubstructures(local_subst);

                    // Fast-forward the line index.
                    iline += (num_subst - 1);
                }
                else
                {
                    for (int i=0; i<num_subst; ++i)
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
    const int nmols = nMolecules();

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
    system.setProperty(map["fileformat"].source(), StringProperty(this->formatName()));

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

/** Internal function used to get the molecule structure for molecule 'imol'. */
MolStructureEditor Mol2::getMolStructure(int imol, const PropertyName &cutting) const
{
    // Get the number of atoms in the molecule.
    const int nats = nAtoms(imol);

    // Make sure that there are atoms in the molecule.
    if (nats == 0)
    {
        throw SireError::program_bug(QObject::tr(
            "Strange - there are no atoms in molecule %1?")
                .arg(imol), CODELOC);
    }

    /* First step is to build the structure of the molecule, i.e.
       the layout of cutgroups, residues and atoms.

       To do this we'll walk through all of the atoms in the molecule
       and work out which substructure they belong to. This is more
       robust than using the substructure record information, which
       is often incomplete, or missing entirely.

       Note that we are assuming that all substructure types can
       be represented as a residue in SIRE.

       If substructure records are present, then we can use them
       to assign the residues into chains.

       TODO: Check which substructure record types can be interpreted
       as a residue.
     */

    MolStructureEditor mol;

    // Create a multimap between a residue (number) and the atoms that it contains.
    QMultiMap<int, int> residues;

    // Create a map between residue number and name.
    QMap<int, QString> res_names;

    // Loop through the atoms in the molecule and store the residues.
    for (int i=0; i<nats; ++i)
    {
        // Get the atom.
        auto atom = molecules[imol].getAtom(i);

        // Insert the atom into the residue map.
        residues.insert(atom.getSubstructureNumber(), i);

        // Insert the residue into the residue name map.
        res_names.insert(atom.getSubstructureNumber(), atom.getSubstructureName());
    }

    // Create a multimap between residues and chains.
    // This information is contained in substructure records and
    // may be missing, or incomplete.
    QMultiMap<QString, int> chains;

    // Create a reverse mapping between residues and chains.
    QMap<int, QString> res_to_chain;

    // Residue index.
    int ires = 0;

    for (int i=0; i<nSubstructures(imol); ++i)
    {
        // Get the substructure.
        auto substructure = molecules[imol].getSubstructure(i);

        // Get the substructure's number.
        auto res_num = substructure.getNumber();

        // Get the name of the chain to which the substructure belongs.
        auto chain = substructure.getChain();

        // Make sure that this is a valid chain.
        if (not chain.isEmpty())
        {
            // Make sure that this substructure exists.
            if (residues.contains(res_num))
            {
                // Insert the substructure into the chain map.
                chains.insert(chain, res_num);
            }
        }
    }

    // Add any chains to the molecule.
    for (auto chain : chains.uniqueKeys())
    {
        mol.add(ChainName(chain));

        // Get a list of the residues that are part of this chain.
        QList<int> chain_residues = chains.values(chain);

        // Create the reverse mapping.
        for (auto residue : chain_residues)
            res_to_chain.insert(residue, chain);
    }

    // Loop over all unique residues in the frame (by number).
    for (auto res_num : residues.uniqueKeys())
    {
        // Extract the residue name.
        auto res_name = res_names[res_num];

        // By default we will use one CutGroup per residue.
        // This may be changed later by the cutting system.
        auto cutgroup = mol.add(CGName(QString::number(ires)));

        // Get a list of the atoms that are part of the residue.
        QList<int> res_atoms = residues.values(res_num);
        qSort(res_atoms);

        // Add the residue to the molecule.
        // Here we use the substructure name from the atom record,
        // so as to avoid any naming inconsistencies.
        auto res = mol.add(ResNum(res_num));
        res.rename(ResName(res_name.trimmed()));

        // Reparent the residue to its chain.
        if (res_to_chain.contains(res_num))
            res.reparent(ChainName(res_to_chain[res_num]));

        // Add each atom in the residue to the molecule.
        for (auto res_atom : res_atoms)
        {
            auto atom = cutgroup.add(AtomNum(molecules[imol].getAtom(res_atom).getNumber()));
            atom.rename(AtomName(molecules[imol].getAtom(res_atom).getName().trimmed()));

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
MolEditor Mol2::getMolecule(int imol, const PropertyMap &map) const
{
    // Get the number of atoms in the molecule.
    const int nats = nAtoms(imol);

    // Make sure that there are atoms in the frame.
    if (nats == 0)
    {
        throw SireError::program_bug(QObject::tr(
            "Strange - there are no atoms in molecule %1?")
                .arg(imol), CODELOC);
    }

    // First, construct the layout of the molecule (sorting of atoms into residues and cutgroups).
    auto mol = this->getMolStructure(imol, map["cutting"]).commit().edit();

    // Get the info object that can map between AtomNum to AtomIdx etc.
    const auto molinfo = mol.info();

    // Rename the molecule.
    mol.rename(MolName(molecules[imol].getName()));

    // Set properties for the molecule.
    mol.setProperty(map["mol-type"], StringProperty(molecules[imol].getMolType()))
       .setProperty(map["charge-type"], StringProperty(molecules[imol].getChargeType()))
       .commit();

    // Add status bits, if present.
    if (not molecules[imol].getStatusBits().isEmpty())
    {
        mol.setProperty(map["mol-status-bits"], StringProperty(molecules[imol].getStatusBits()))
           .commit();
    }

    // Add comments, if present.
    if (not molecules[imol].getComment().isEmpty())
    {
        mol.setProperty(map["mol-comment"], StringProperty(molecules[imol].getComment()))
           .commit();
    }

    // Instantiate the atom property objects that we need.
    AtomCoords         coords(molinfo);
    AtomCharges        charges(molinfo);
    AtomElements       elements(molinfo);
    AtomFloatProperty  occupancies(molinfo);
    AtomStringProperty types(molinfo);
    AtomStringProperty status_bits(molinfo);

    // Now loop through the atoms in the molecule and set each property.
    for (int i=0; i<nats; ++i)
    {
        // Get the current atom.
        auto atom = molecules[imol].getAtom(i);

        // Determine the CGAtomIdx for this atom.
        auto cgatomidx = molinfo.cgAtomIdx(AtomNum(atom.getNumber()));

        // Set the atom properties.
        coords.set(cgatomidx, atom.getCoord());
        charges.set(cgatomidx, double(atom.getCharge()) * SireUnits::mod_electron);
        types.set(cgatomidx, atom.getType());
        status_bits.set(cgatomidx, atom.getStatusBits());

        // The element is usuall the first character of the atom,
        // unless the name starts with a digit, in which case it's the second.
        auto name = atom.getName();
        if (name[0].isDigit()) elements.set(cgatomidx, Element(QString(name[1])));
        else                   elements.set(cgatomidx, Element(QString(name[0])));
    }

    // Instantiate the residue property objects that we need.
    ResStringProperty subst_types(molinfo);
    ResStringProperty chain_sub_types(molinfo);
    ResStringProperty subst_status_bits(molinfo);
    ResStringProperty comments(molinfo);
    ResIntProperty    dict_types(molinfo);
    ResIntProperty    inter_bonds(molinfo);

    // Now loop through the substructures in the molecule and set each property.
    for (int i=0; i<nSubstructures(imol); ++i)
    {
        // Get the current substructure.
        auto subst = molecules[imol].getSubstructure(i);

        // Determine the ResIdx for this residue.
        auto residx = molinfo.resIdx(ResNum(subst.getNumber()));

        // Set the residue properties.
        subst_types.set(residx, subst.getType());
        dict_types.set(residx, subst.getDictType());
        chain_sub_types.set(residx, subst.getChainSubType());
        inter_bonds.set(residx, subst.getInterBonds());
        subst_status_bits.set(residx, subst.getStatusBits());
        comments.set(residx, subst.getComment());
    }

    return mol.setProperty(map["coordinates"], coords)
              .setProperty(map["charge"], charges)
              .setProperty(map["element"], elements)
              .setProperty(map["sybyl-atom-type"], types)
              .setProperty(map["atom-status-bits"], status_bits)
              .setProperty(map["res-type"], subst_types)
              .setProperty(map["res-dict-type"], dict_types)
              .setProperty(map["chain-sub-type"], chain_sub_types)
              .setProperty(map["res-inter-bonds"], inter_bonds)
              .setProperty(map["res-status-bits"], subst_status_bits)
              .setProperty(map["res-comment"], comments)
              .commit();
}

/** Internal function used to parse a Sire molecule view into a Mol2 molecule using
    the parameters in the property map. */
void Mol2::parseMolecule(Mol2Molecule &mol2_mol, const SireMol::Molecule &sire_mol, int imol,
    QVector<QString> &atom_lines, QVector<QString> &substructure_lines, QStringList &errors,
    const SireBase::PropertyMap &map)
{
    // Store the number of atoms in the molecule.
    int num_atoms = sire_mol.nAtoms();

    // Early exit.
    if (num_atoms == 0) return;

    // Store the number of residues in the molecule.
    int num_res = sire_mol.nResidues();

    // Resize the data record containers.
    atom_lines.resize(num_atoms);
    substructure_lines.resize(num_res);

    /* We first need to check whether there are duplicate residue numbers.
       If so, then we need to number residues by index since Mol2 assumes unique numbering.
       Note that this means that you can't always perform a round-trip file conversion
       and retain the same residue numbers, e.g. PDB -> Mol2 -> PDB. (PDB can have
       duplicate residue numbers in different chains, or with different insertion codes.
     */

    // Whether to number by residue index.
    bool is_idx = false;

    // A vector of the residue numbers seen to date.
    QVector<int> res_nums;

    for (int i=0; i<num_res; ++i)
    {
        int res_num = sire_mol.residue(ResIdx(i)).number().value();

        if (not res_nums.contains(res_num))
        {
            res_nums.append(res_num);
        }
        else
        {
            errors.append(QObject::tr("Warning: there are duplicate residue "
                "numbers in molecule %1, converting to unique indices.").arg(imol));

            is_idx = true;
            break;
        }
    }

    if (usesParallel())
    {
        QMutex mutex;

        // Local data storage.
        QVector<Mol2Atom> local_atoms(num_atoms);
        QVector<Mol2Substructure> local_subst(num_res);

        tbb::parallel_for( tbb::blocked_range<int>(0, num_atoms),
                        [&](const tbb::blocked_range<int> &r)
        {
            // Create local data objects.
            QStringList local_errors;

            // Convert each atom into a Mol2Atom object
            // and generate a Mol2 data record.
            for (int i=r.begin(); i<r.end(); ++i)
            {
                local_atoms[i] = Mol2Atom(sire_mol.atom(AtomIdx(i)), map, local_errors, is_idx);
                atom_lines[i] = local_atoms[i].toMol2Record();
            }

            if (not local_errors.isEmpty())
            {
                // Acquire a lock.
                QMutexLocker lkr(&mutex);

                // Update the warning messages.
                errors += local_errors;
            }
        });

        tbb::parallel_for( tbb::blocked_range<int>(0, num_res),
                        [&](const tbb::blocked_range<int> &r)
        {
            // Create local data objects.
            QStringList local_errors;

            // Convert each residue into a Mol2Substructure object
            // and generate a Mol2 data record.
            for (int i=r.begin(); i<r.end(); ++i)
            {
                local_subst[i] = Mol2Substructure(sire_mol.residue(ResIdx(i)), map, local_errors, is_idx);
                substructure_lines[i] = local_subst[i].toMol2Record();
            }

            if (not local_errors.isEmpty())
            {
                // Acquire a lock.
                QMutexLocker lkr(&mutex);

                // Update the warning messages.
                errors += local_errors;
            }
        });
    }
    else
    {
        // Loop over all of the atoms.
        for (int i=0; i<num_atoms; ++i)
        {
            // Initalise a Mol2Atom.
            Mol2Atom atom(sire_mol.atom(AtomIdx(i)), map, errors, is_idx);

            // Generate a Mol2 atom data record.
            atom_lines[i] = atom.toMol2Record();
        }

        // Loop over all of the residues.
        for (int i=0; i<num_res; ++i)
        {
            // Initalise a Mol2Substructure.
            Mol2Substructure subst(sire_mol.residue(ResIdx(i)), map, errors, is_idx);

            // Generate a Mol2 substructure data record.
            substructure_lines[i] = subst.toMol2Record();
        }
    }
}
