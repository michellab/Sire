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
class Mol2Molecule;
class Mol2Substructure;
class Mol2;
}

namespace SireMol
{
class Atom;
class MolEditor;
class MoleculeInfoData;
class MoleculeView;
class Residue;
}

SIREIO_EXPORT QDataStream& operator<<(QDataStream&, const SireIO::Mol2Atom&);
SIREIO_EXPORT QDataStream& operator>>(QDataStream&, SireIO::Mol2Atom&);

SIREIO_EXPORT QDataStream& operator<<(QDataStream&, const SireIO::Mol2Bond&);
SIREIO_EXPORT QDataStream& operator>>(QDataStream&, SireIO::Mol2Bond&);

SIREIO_EXPORT QDataStream& operator<<(QDataStream&, const SireIO::Mol2Molecule&);
SIREIO_EXPORT QDataStream& operator>>(QDataStream&, SireIO::Mol2Molecule&);

SIREIO_EXPORT QDataStream& operator<<(QDataStream&, const SireIO::Mol2Substructure&);
SIREIO_EXPORT QDataStream& operator>>(QDataStream&, SireIO::Mol2Substructure&);

SIREIO_EXPORT QDataStream& operator<<(QDataStream&, const SireIO::Mol2&);
SIREIO_EXPORT QDataStream& operator>>(QDataStream&, SireIO::Mol2&);

namespace SireIO
{

/** This class provides functionality for reading/writing
    Mol2 ATOM records.

    @author Lester Hedges
*/
class SIREIO_EXPORT Mol2Atom
{

friend SIREIO_EXPORT QDataStream& ::operator<<(QDataStream&, const Mol2Atom&);
friend SIREIO_EXPORT QDataStream& ::operator>>(QDataStream&, Mol2Atom&);

public:
    /** Default constructor. */
    Mol2Atom();

    /** Constructor. */
    Mol2Atom(const QString &line, QStringList &errors);

    /** Constructor. */
    Mol2Atom(const SireMol::Atom &atom, const PropertyMap &map,
        QStringList &errors, bool is_idx = false);

    /** Generate a Mol2 record from the atom data. */
    QString toMol2Record() const;

    static const char* typeName();

    /** Get the atom number. */
    int getNumber() const;

    /** Get the atom name. */
    QString getName() const;

    /** Get the atom coordinates. */
    SireMaths::Vector getCoord() const;

    /** Get the SYBYL atom type. */
    QString getType() const;

    /** Get the number of the substructure containing the atom. */
    qint64 getSubstructureNumber() const;

    /** Get the name of the substructure containing the atom. */
    QString getSubstructureName() const;

    /** Get the atom charge. */
    double getCharge() const;

    /** Get the status bits. */
    QString getStatusBits() const;

private:
    /** The original Mol2 record used to instantiate the atom. */
    QString record;

    /** The ID number of the atom at the time the file was created. */
    qint64 number;

    /** The name of the atom. */
    QString name;

    /** Coordinates. */
    SireMaths::Vector coord;

    /** The SYBYL atom type. */
    QString type;

    /** The ID number of the substructure containing the atom. */
    qint64 subst_id;

    /** The name of the substructure containing the atom. */
    QString subst_name;

    /** The charge on the atom. */
    double charge;

    /** The internal SYBYL status bits. */
    QString status_bits;
};

/** This class provides functionality for reading/writing
    Mol2 BOND records.

    @author Lester Hedges
*/
class SIREIO_EXPORT Mol2Bond
{

friend SIREIO_EXPORT QDataStream& ::operator<<(QDataStream&, const Mol2Bond&);
friend SIREIO_EXPORT QDataStream& ::operator>>(QDataStream&, Mol2Bond&);

public:
    /** Default constructor. */
    Mol2Bond();

    /** Constructor. */
    Mol2Bond(const QString &line, QStringList &errors);

    /** Generate a Mol2 record from the bond data. */
    QString toMol2Record() const;

    static const char* typeName();

    /** Get the bond ID. */
    qint64 getID() const;

    /** Get the ID of the origin atom. */
    qint64 getOrigin() const;

    /** Get the ID of the target atom. */
    qint64 getTarget() const;

    /** Get the bond type. */
    QString getType() const;

private:
    /** The original Mol2 record used to instantiate the bond. */
    QString record;

    /** The ID number of the bond at the time the file was created. */
    qint64 number;

    /** The ID number of the atom at one end of the bond. */
    qint64 origin;

    /** The ID number of the atom at the other end of the bond. */
    qint64 target;

    /** The SYBYL bond type. */
    QString type;

    /** The ID number of the substructure containing the atom. */
    qint64 subst_id;

    /** The internal SYBYL status bits. */
    QString status_bits;
};

/** This class provides functionality for reading/writing
    Mol2 MOLECULE records.

    @author Lester Hedges
*/
class SIREIO_EXPORT Mol2Molecule
{

friend SIREIO_EXPORT QDataStream& ::operator<<(QDataStream&, const Mol2Molecule&);
friend SIREIO_EXPORT QDataStream& ::operator>>(QDataStream&, Mol2Molecule&);

public:
    /** Default constructor. */
    Mol2Molecule();

    /** Constructor (from a Mol2 data record). */
    Mol2Molecule(const QVector<QString> &lines, QStringList &errors,
        int &num_records, int imol = -1);

    /** Constructor (from a Sire Molecule). */
    Mol2Molecule(const SireMol::Molecule &mol, const PropertyMap &map,
        QStringList &errors, int imol = -1);

    /** Generate a Mol2 record from the molecule data. */
    QVector<QString> toMol2Record() const;

    static const char* typeName();

    /** Get the number of atoms in the molecule. */
    int nAtoms() const;

    /** Get the number of bonds in the molecule. */
    int nBonds() const;

    /** Get the number of substructures in the molecule. */
    int nSubstructures() const;

    /** Append an atom to the molecule. */
    void appendAtom(const Mol2Atom &atom);

    /** Append a vector of atoms to the molecule. */
    void appendAtoms(const QVector<Mol2Atom> &atoms);

    /** Append a bond to the molecule. */
    void appendBond(const Mol2Bond &bond);

    /** Append a vector of bonds to the molecule. */
    void appendBonds(const QVector<Mol2Bond> &bonds);

    /** Append a substructure to the molecule. */
    void appendSubstructure(const Mol2Substructure &substructure);

    /** Append a vector of substructures to the molecule. */
    void appendSubstructures(const QVector<Mol2Substructure> &substructures);

    /** Get the name of the molecule. */
    QString getName() const;

    /** Get the molecule type. */
    QString getMolType() const;

    /** Get the charge type. */
    QString getChargeType() const;

    /** Get the status bits. */
    QString getStatusBits() const;

    /** Get the comment. */
    QString getComment() const;

    /** Get the atoms. */
    QVector<Mol2Atom> getAtoms() const;

    /** Get a specific atom. */
    Mol2Atom getAtom(int i) const;

    /** Get the bonds. */
    QVector<Mol2Bond> getBonds() const;

    /** Get a specific bond. */
    Mol2Bond getBond(int i) const;

    /** Get the substructures. */
    QVector<Mol2Substructure> getSubstructures() const;

    /** Get a specific substructure. */
    Mol2Substructure getSubstructure(int i) const;

private:
    // Record data.

    /** The original Mol2 record used to instantiate the molecule. */
    QVector<QString> record;

    /** The name of the molecule. */
    QString name;

    /** The number of atoms in the molecule. */
    qint64 num_atoms;

    /** The number of bonds in the molecule. */
    qint64 num_bonds;

    /** The number of substructures in the molecule. */
    qint64 num_subst;

    // TODO: Maybe delete this data member.
    /** The number of features in the molecule. */
    qint64 num_feats;

    // TODO: Maybe delete this data member.
    /** The number of sets in the molecule. */
    qint64 num_sets;

    /** The molecule type. */
    QString mol_type;

    /** The charge type. */
    QString charge_type;

    /** The internal SYBYL status bits. */
    QString status_bits;

    /** Comments about the molecule. */
    QString comment;

    // The objects that make up the molecule.

    /** Atom data. */
    QVector<Mol2Atom> atoms;

    /** Bond data. */
    QVector<Mol2Bond> bonds;

    /** substructure data. */
    QVector<Mol2Substructure> substructures;
};

/** This class provides functionality for reading/writing
    Mol2 SUBSTRUCTURE records.

    @author Lester Hedges
*/
class SIREIO_EXPORT Mol2Substructure
{

friend SIREIO_EXPORT QDataStream& ::operator<<(QDataStream&, const Mol2Substructure&);
friend SIREIO_EXPORT QDataStream& ::operator>>(QDataStream&, Mol2Substructure&);

public:
    /** Default constructor. */
    Mol2Substructure();

    /** Constructor. */
    Mol2Substructure(const QString &line, QStringList &errors);

    /** Constructor. */
    Mol2Substructure(const SireMol::Residue &res, const PropertyMap &map,
        QStringList &errors, bool is_idx = false);

    /** Generate a Mol2 record from the feature data. */
    QString toMol2Record() const;

    static const char* typeName();

    /** Get the substructure number. */
    qint64 getNumber() const;

    /** Get the substructure name. */
    QString getName() const;

    /** Get the root atom index. */
    qint64 getRootAtom() const;

    /** Get the substructure type. */
    QString getType() const;

    /** Get the dictionary type. */
    qint64 getDictType() const;

    /** Get the substructure chain. */
    QString getChain() const;

    /** Get the subtype of the chain. */
    QString getChainSubType() const;

    /** Get the number of inter substructure bonds. */
    qint64 getInterBonds() const;

    /** Get the status bits. */
    QString getStatusBits() const;

    /** Get the comment. */
    QString getComment() const;

private:
    /** The original Mol2 record used to instantiate the substructure. */
    QString record;

    /** The ID number of the substructure. */
    qint64 number;

    /** The name of the substructure. */
    QString name;

    /** The ID number of the root atom. */
    qint64 root_atom;

    /** The type of the substructure. */
    QString type;

    /** The dictionary type associated with the substructure. */
    qint64 dict_type;

    /** The chain to which the substructure belongs. */
    QString chain;

    /** The sub type of the chain. */
    QString sub_type;

    /** The number of inter substructure bonds. */
    qint64 num_inter_bonds;

    /** The internal SYBYL status bits. */
    QString status_bits;

    /** Comments about the substructure. */
    QString comment;
};

/** This class holds a parser for reading and writing Tripos Mol2 files.

    @author Lester Hedges
*/
class SIREIO_EXPORT Mol2 : public SireBase::ConcreteProperty<Mol2,MoleculeParser>
{

friend SIREIO_EXPORT QDataStream& ::operator<<(QDataStream&, const Mol2&);
friend SIREIO_EXPORT QDataStream& ::operator>>(QDataStream&, Mol2&);

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
    QVector<QString> toLines() const;

    QString formatName() const;
    QString formatDescription() const;
    QStringList formatSuffix() const;

    bool isTopology() const;
    bool isFrame() const;

    int nFrames() const;
    SireMol::Frame getFrame(int i) const;

    /** Return the number of molecules in the system. */
    int nMolecules() const;

    /** Return the number of atoms in each molecule. */
    QVector<int> nMolAtoms() const;

    /** Return the number of atoms in a specific molecule. */
    int nAtoms(int i) const;

    /** Return the total number of atoms in all molecules. */
    int nAtoms() const;

    /** Return the number of bonds in each molecule. */
    QVector<int> nMolBonds() const;

    /** Return the number of bonds in a specific molecule. */
    int nBonds(int i) const;

    /** Return the total number of bonds in all molecules. */
    int nBonds() const;

    /** Return the number of substructures in each molecule. */
    QVector<int> nMolSubstructures() const;

    /** Return the number of substructures in a specific molecule. */
    int nSubstructures(int i) const;

    /** Return the total number of substructures in all molecules. */
    int nSubstructures() const;

protected:
    SireSystem::System startSystem(const PropertyMap &map) const;
    void addToSystem(SireSystem::System &system, const PropertyMap &map) const;

private:
    void assertSane() const;
    void parseLines(const PropertyMap &map);
    void parseMolecule(Mol2Molecule &mol2_mol, const SireMol::Molecule &sire_mol, int imol,
        QVector<QString> &atom_lines, QVector<QString> &substructure_lines, QStringList &errors,
        const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SireMol::MolStructureEditor getMolStructure(int imol, QHash<int, int> &res_map,
        const SireBase::PropertyName &cutting) const;

    SireMol::MolEditor getMolecule(int imol,
        const PropertyMap &map = PropertyMap()) const;

    /** The molecular data object. */
    QVector<Mol2Molecule> molecules;

    /** Any warnings that were raised when reading the file. */
    QStringList parse_warnings;
};

}

Q_DECLARE_METATYPE( SireIO::Mol2Atom )
Q_DECLARE_METATYPE( SireIO::Mol2Bond )
Q_DECLARE_METATYPE( SireIO::Mol2Molecule )
Q_DECLARE_METATYPE( SireIO::Mol2Substructure )
Q_DECLARE_METATYPE( SireIO::Mol2 )

SIRE_EXPOSE_CLASS( SireIO::Mol2 )

SIRE_END_HEADER

#endif
