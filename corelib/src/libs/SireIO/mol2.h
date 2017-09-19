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

#ifndef SIREIO_MOL2_H
#define SIREIO_MOL2_H

#include "moleculeparser.h"

#include "SireMaths/vector.h"

SIRE_BEGIN_HEADER

namespace SireIO
{
class Mol2Atom;
class Mol2Bond;
class Mol2Feature;
class Mol2Molecule;
class Mol2Set;
class Mol2SubStructure;
class Mol2;
}

QDataStream& operator<<(QDataStream&, const SireIO::Mol2Atom&);
QDataStream& operator>>(QDataStream&, SireIO::Mol2Atom&);

QDataStream& operator<<(QDataStream&, const SireIO::Mol2Bond&);
QDataStream& operator>>(QDataStream&, SireIO::Mol2Bond&);

QDataStream& operator<<(QDataStream&, const SireIO::Mol2Feature&);
QDataStream& operator>>(QDataStream&, SireIO::Mol2Feature&);

QDataStream& operator<<(QDataStream&, const SireIO::Mol2Molecule&);
QDataStream& operator>>(QDataStream&, SireIO::Mol2Molecule&);

QDataStream& operator<<(QDataStream&, const SireIO::Mol2Set&);
QDataStream& operator>>(QDataStream&, SireIO::Mol2Set&);

QDataStream& operator<<(QDataStream&, const SireIO::Mol2SubStructure&);
QDataStream& operator>>(QDataStream&, SireIO::Mol2SubStructure&);

QDataStream& operator<<(QDataStream&, const SireIO::Mol2&);
QDataStream& operator>>(QDataStream&, SireIO::Mol2&);

namespace SireIO
{

/** This class provides functionality for reading/writing
    Mol2 ATOM records.

    @author Lester Hedges
*/
class SIREIO_EXPORT Mol2Atom
{

friend QDataStream& ::operator<<(QDataStream&, const Mol2Atom&);
friend QDataStream& ::operator>>(QDataStream&, Mol2Atom&);

public:
    /** Default constructor. */
    Mol2Atom();

    /** Constructor. */
    Mol2Atom(const QString &line, QStringList &errors);

    /** Generate a Mol2 record from the atom data. */
    QString toMol2Record() const;

    /** Generate a string representation of the object. */
    QString toString() const;

    static const char* typeName();

private:
    /** The original Mol2 record used to instantiate the atom. */
    QString record;

    /** The ID number of the atom at the time the file was created. */
    qint64 number;

    /** The name of the atom. */
    QString name;

    /** Coordinates. */
    SireMaths::Vector coord;

    /** The SYBL atom type. */
    QString type;

    /** The ID number of the substructure containing the atom. */
    qint64 subst_id;

    /** The name of the substructure containing the atom. */
    QString subst_name;

    /** The charge on the atom. */
    double charge;

    /** The internal SYBL status bits. */
    QString status_bit;
};

/** This class provides functionality for reading/writing
    Mol2 BOND records.

    @author Lester Hedges
*/
class SIREIO_EXPORT Mol2Bond
{

friend QDataStream& ::operator<<(QDataStream&, const Mol2Bond&);
friend QDataStream& ::operator>>(QDataStream&, Mol2Bond&);

public:
    /** Default constructor. */
    Mol2Bond();

    /** Constructor. */
    Mol2Bond(const QString &line, QStringList &errors);

    /** Generate a Mol2 record from the bond data. */
    QString toMol2Record() const;

    /** Generate a string representation of the object. */
    QString toString() const;

    static const char* typeName();

private:
    /** The original Mol2 record used to instantiate the bond. */
    QString record;

    /** The ID number of the bond at the time the file was created. */
    qint64 number;

    /** The ID number of the atom at one end of the bond. */
    qint64 origin_atom;

    /** The ID number of the atom at the other end of the bond. */
    qint64 target_atom;

    /** The SYBL bond type. */
    QString type;

    /** The ID number of the substructure containing the atom. */
    qint64 subst_id;

    /** The internal SYBL status bits. */
    QString status_bit;
};

/** This class provides functionality for reading/writing
    Mol2 FEATURE records.

    @author Lester Hedges
*/
class SIREIO_EXPORT Mol2Feature
{

friend QDataStream& ::operator<<(QDataStream&, const Mol2Feature&);
friend QDataStream& ::operator>>(QDataStream&, Mol2Feature&);

public:
    /** Default constructor. */
    Mol2Feature();

    /** Constructor. */
    Mol2Feature(const QStringList &lines, QStringList &errors);

    /** Generate a Mol2 record from the feature data. */
    QStringList toMol2Record() const;

    /** Generate a string representation of the object. */
    QString toString() const;

    static const char* typeName();

private:
};

/** This class provides functionality for reading/writing
    Mol2 MOLECULE records.

    @author Lester Hedges
*/
class SIREIO_EXPORT Mol2Molecule
{

friend QDataStream& ::operator<<(QDataStream&, const Mol2Molecule&);
friend QDataStream& ::operator>>(QDataStream&, Mol2Molecule&);

public:
    /** Default constructor. */
    Mol2Molecule();

    /** Constructor. */
    Mol2Molecule(const QVector<QString> &lines, QStringList &errors);

    /** Generate a Mol2 record from the molecule data. */
    QStringList toMol2Record() const;

    /** Generate a string representation of the object. */
    QString toString() const;

    static const char* typeName();

private:
    // Record data.

    /** The original Mol2 record used to instantiate the molecule. */
    QStringList record;

    /** The name of the molecule. */
    QString name;

    /** The number of atoms in the molecule. */
    qint64 num_atoms;

    /** The number of bonds in the molecule. */
    qint64 num_bonds;

    /** The number of sub-structures in the molecule. */
    qint64 num_subst;

    /** The number of features in the molecule. */
    qint64 num_feats;

    /** The number of sets in the molecule. */
    qint64 num_sets;

    /** The molecule type. */
    QString mol_type;

    /** The charge type. */
    QString charge_type;

    /** The internal SYBL status bits. */
    QString status_bit;

    /** Comments about the molecule. */
    QString comment;

    // The objects that make up the molecule.

    /** Atom data. */
    QVector<Mol2Atom> atoms;

    /** Bond data. */
    QVector<Mol2Bond> bonds;

    /** Feature data. */
    QVector<Mol2Feature> features;

    /** Set data. */
    QVector<Mol2Set> sets;

    /** Sub-structure data. */
    QVector<Mol2SubStructure> substructures;
};

/** This class provides functionality for reading/writing
    Mol2 SET records.

    @author Lester Hedges
*/
class SIREIO_EXPORT Mol2Set
{

friend QDataStream& ::operator<<(QDataStream&, const Mol2Set&);
friend QDataStream& ::operator>>(QDataStream&, Mol2Set&);

public:
    /** Default constructor. */
    Mol2Set();

    /** Constructor. */
    Mol2Set(const QStringList &lines, QStringList &errors);

    /** Generate a Mol2 record from the feature data. */
    QStringList toMol2Record() const;

    /** Generate a string representation of the object. */
    QString toString() const;

    static const char* typeName();

private:
    /** The original Mol2 record used to instantiate the set. */
    QStringList record;

    /** The name of the set. */
    QString name;

    /** The set type (STATIC, or DYNAMIC). */
    QString type;

    /** The sub type of the set. */
    QString sub_type;

    /** The internal SYBL status bits. */
    QString status_bit;

    /** Comments about the set. */
    QString comment;

    /** The number of members of the set (STATIC only). */
    qint64 num_members;

    /** The ID of a member of the set (STATIC only). */
    qint64 member_id;

    /** The rule defining the DYNAMIC set. */
    QString rule;
};

/** This class provides functionality for reading/writing
    Mol2 SUBSTRUCTURE records.

    @author Lester Hedges
*/
class SIREIO_EXPORT Mol2SubStructure
{

friend QDataStream& ::operator<<(QDataStream&, const Mol2SubStructure&);
friend QDataStream& ::operator>>(QDataStream&, Mol2SubStructure&);

public:
    /** Default constructor. */
    Mol2SubStructure();

    /** Constructor. */
    Mol2SubStructure(const QString &line, QStringList &errors);

    /** Generate a Mol2 record from the feature data. */
    QString toMol2Record() const;

    /** Generate a string representation of the object. */
    QString toString() const;

    static const char* typeName();

private:
    /** The original Mol2 record used to instantiate the sub-structure. */
    QString record;

    /** The ID number of the sub-structure. */
    qint64 number;

    /** The name of the sub-structure. */
    QString name;

    /** The ID number of the root atom. */
    qint64 root_atom;

    /** The type of the sub-structure. */
    QString type;

    /** The dictionary type associated with the sub-structure. */
    qint64 dict_type;

    /** The chain to which the sub-structure belongs. */
    QString chain;

    /** The sub type of the sub-structure. */
    QString sub_type;

    /** The number of inter sub-structure bonds. */
    qint64 num_inter_bonds;

    /** The internal SYBL status bits. */
    QString status_bit;

    /** Comments about the set. */
    QString comment;
};

/** This class holds a parser for reading and writing Tripos Mol2 files.

    @author Lester Hedges
*/
class SIREIO_EXPORT Mol2 : public SireBase::ConcreteProperty<Mol2,MoleculeParser>
{

friend QDataStream& ::operator<<(QDataStream&, const Mol2&);
friend QDataStream& ::operator>>(QDataStream&, Mol2&);

public:
    Mol2();
    Mol2(const QString &filename,
         const PropertyMap &map = PropertyMap());

    Mol2(const QStringList &lines,
         const PropertyMap &map = PropertyMap());
    Mol2(const SireSystem::System &system,
         const PropertyMap &map = PropertyMap());

    Mol2(const Mol2 &other);

    ~Mol2();

    Mol2& operator=(const Mol2 &other);

    bool operator==(const Mol2 &other) const;
    bool operator!=(const Mol2 &other) const;

    static const char* typeName();

    const char* what() const;

    MoleculeParserPtr construct(const QString &filename,
                                const PropertyMap &map) const;

    MoleculeParserPtr construct(const QStringList &lines,
                                const PropertyMap &map) const;

    MoleculeParserPtr construct(const SireSystem::System &system,
                                const PropertyMap &map) const;

    QString toString() const;

    QString formatName() const;
    QString formatDescription() const;
    QStringList formatSuffix() const;

protected:
    SireSystem::System startSystem(const PropertyMap &map) const;
    void addToSystem(SireSystem::System &system, const PropertyMap &map) const;

private:
    void assertSane() const;
    void parseLines(const PropertyMap &map);

    /** The molecular data object. */
    QVector<Mol2Molecule> molecules;

    /** Any warnings that were raised when reading the file. */
    QStringList parse_warnings;
};

}

Q_DECLARE_METATYPE( SireIO::Mol2Atom )
Q_DECLARE_METATYPE( SireIO::Mol2Bond )
Q_DECLARE_METATYPE( SireIO::Mol2Feature )
Q_DECLARE_METATYPE( SireIO::Mol2Molecule )
Q_DECLARE_METATYPE( SireIO::Mol2Set )
Q_DECLARE_METATYPE( SireIO::Mol2SubStructure )
Q_DECLARE_METATYPE( SireIO::Mol2 )

SIRE_EXPOSE_CLASS( SireIO::Mol2Atom )
SIRE_EXPOSE_CLASS( SireIO::Mol2Bond )
SIRE_EXPOSE_CLASS( SireIO::Mol2Feature )
SIRE_EXPOSE_CLASS( SireIO::Mol2Molecule )
SIRE_EXPOSE_CLASS( SireIO::Mol2Set )
SIRE_EXPOSE_CLASS( SireIO::Mol2SubStructure )
SIRE_EXPOSE_CLASS( SireIO::Mol2 )

SIRE_END_HEADER

#endif