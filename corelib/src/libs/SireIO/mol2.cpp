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

#include "SireMM/mol2params.h"

#include "SireSystem/system.h"

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
static const RegisterMetaType<Mol2Feature> r_mol2feature(NO_ROOT);
static const RegisterMetaType<Mol2Molecule> r_mol2molecule(NO_ROOT);
static const RegisterMetaType<Mol2Set> r_mol2set(NO_ROOT);
static const RegisterMetaType<Mol2SubStructure> r_mol2subst(NO_ROOT);

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

    sds << mol2bond.record << mol2bond.number << mol2bond.origin_atom
        << mol2bond.target_atom << mol2bond.type << mol2bond.subst_id
        << mol2bond.status_bit;

    return ds;
}

QDataStream SIREIO_EXPORT &operator>>(QDataStream &ds, Mol2Bond &mol2bond)
{
    VersionID v = readHeader(ds, r_mol2bond);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> mol2bond.record >> mol2bond.number >> mol2bond.origin_atom
            >> mol2bond.target_atom >> mol2bond.type >> mol2bond.subst_id
            >> mol2bond.status_bit;
    }
    else
        throw version_error(v, "1", r_mol2bond, CODELOC);

    return ds;
}

QDataStream SIREIO_EXPORT &operator<<(QDataStream &ds, const Mol2Feature &mol2feature)
{
    writeHeader(ds, r_mol2feature, 1);

    SharedDataStream sds(ds);

    //sds << ;

    return ds;
}

QDataStream SIREIO_EXPORT &operator>>(QDataStream &ds, Mol2Feature &mol2feature)
{
    VersionID v = readHeader(ds, r_mol2feature);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        //sds >> ;
    }
    else
        throw version_error(v, "1", r_mol2feature, CODELOC);

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
        << mol2molecule.bonds << mol2molecule.features << mol2molecule.sets
        << mol2molecule.substructures;

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
            >> mol2molecule.bonds >> mol2molecule.features >> mol2molecule.sets
            >> mol2molecule.substructures;
    }
    else
        throw version_error(v, "1", r_mol2molecule, CODELOC);

    return ds;
}

QDataStream SIREIO_EXPORT &operator<<(QDataStream &ds, const Mol2Set &mol2set)
{
    writeHeader(ds, r_mol2set, 1);

    SharedDataStream sds(ds);

    sds << mol2set.record << mol2set.name << mol2set.type << mol2set.sub_type
        << mol2set.status_bit << mol2set.comment << mol2set.num_members
        << mol2set.member_id << mol2set.rule;

    return ds;
}

QDataStream SIREIO_EXPORT &operator>>(QDataStream &ds, Mol2Set &mol2set)
{
    VersionID v = readHeader(ds, r_mol2set);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> mol2set.record >> mol2set.name >> mol2set.type >> mol2set.sub_type
            >> mol2set.status_bit >> mol2set.comment >> mol2set.num_members
            >> mol2set.member_id >> mol2set.rule;
    }
    else
        throw version_error(v, "1", r_mol2set, CODELOC);

    return ds;
}

QDataStream SIREIO_EXPORT &operator<<(QDataStream &ds, const Mol2SubStructure &mol2subst)
{
    writeHeader(ds, r_mol2subst, 1);

    SharedDataStream sds(ds);

    sds << mol2subst.record << mol2subst.number << mol2subst.name << mol2subst.root_atom
        << mol2subst.type << mol2subst.dict_type << mol2subst.chain << mol2subst.sub_type
        << mol2subst.num_inter_bonds << mol2subst.status_bit << mol2subst.comment;

    return ds;
}

QDataStream SIREIO_EXPORT &operator>>(QDataStream &ds, Mol2SubStructure &mol2subst)
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

    ds << static_cast<const MoleculeParser&>(mol2);

    return ds;
}

QDataStream SIREIO_EXPORT &operator>>(QDataStream &ds, Mol2 &mol2)
{
    VersionID v = readHeader(ds, r_mol2);

    if (v == 1)
    {
        ds >> static_cast<MoleculeParser&>(mol2);
    }
    else
        throw version_error(v, "1", r_mol2, CODELOC);

    return ds;
}

/** Default constructor. */
Mol2Atom::Mol2Atom()
{
}

/** Constructor. */
Mol2Atom::Mol2Atom(const QString &line, QStringList &errors)
{
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

/** Default constructor. */
Mol2Bond::Mol2Bond()
{
}

/** Constructor. */
Mol2Bond::Mol2Bond(const QString &line, QStringList &errors)
{
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

/** Default constructor. */
Mol2Feature::Mol2Feature()
{
}

/** Constructor. */
Mol2Feature::Mol2Feature(const QStringList &lines, QStringList &errors)
{
}

/** Generate a Mol2 record from the feature data. */
QStringList Mol2Feature::toMol2Record() const
{
    //return record;
}

/** Generate a string representation of the object. */
QString Mol2Feature::toString() const
{
    return QObject::tr("Mol2Feature::null");
}

const char* Mol2Feature::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Mol2Feature>() );
}

/** Default constructor. */
Mol2Molecule::Mol2Molecule()
{
}

/** Constructor. */
Mol2Molecule::Mol2Molecule(const QStringList &lines, QStringList &errors)
{
}

/** Generate a Mol2 record from the molecule data. */
QStringList Mol2Molecule::toMol2Record() const
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

/** Default constructor. */
Mol2Set::Mol2Set()
{
}

/** Constructor. */
Mol2Set::Mol2Set(const QStringList &lines, QStringList &errors)
{
}

/** Generate a Mol2 record from the set data. */
QStringList Mol2Set::toMol2Record() const
{
    return record;
}

/** Generate a string representation of the object. */
QString Mol2Set::toString() const
{
    return QObject::tr("Mol2Set::null");
}

const char* Mol2Set::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Mol2Set>() );
}

/** Default constructor. */
Mol2SubStructure::Mol2SubStructure()
{
}

/** Constructor. */
Mol2SubStructure::Mol2SubStructure(const QString &line, QStringList &errors)
{
}

/** Generate a Mol2 record from the sub-structure data. */
QString Mol2SubStructure::toMol2Record() const
{
    return record;
}

/** Generate a string representation of the object. */
QString Mol2SubStructure::toString() const
{
    return QObject::tr("Mol2SubStructure::null");
}

const char* Mol2SubStructure::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Mol2SubStructure>() );
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

       Mol2 is a "free" file format, so records could potentially appear in any order.
       Typically there will be a MOLECULE record, which is followed by ATOM and, e.g.
       BOND records. For simplicity, we'll assume that there is a sensible ordering of
       data, then deal with any inconsistencies afterwards.
     */

    // The Mol2 record indicator. All record types start with this string.
    const QString record_indicator("@<TRIPOS>");

    // Loop through all lines in the file.
    for (int iline=0; iline<lines().count(); ++iline)
    {
        // Store a reference to the line.
        const QString &line = lines()[iline];
    }

    this->setScore(0);
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
