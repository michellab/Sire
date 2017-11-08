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

#ifndef SIREIO_PDB2_H
#define SIREIO_PDB2_H

#include "moleculeparser.h"

#include "SireMaths/vector.h"

SIRE_BEGIN_HEADER

namespace SireIO
{
class PDBAtom;
class PDB2;
}

namespace SireMol
{
class Atom;
class MolEditor;
class MoleculeInfoData;
}

QDataStream& operator<<(QDataStream&, const SireIO::PDBAtom&);
QDataStream& operator>>(QDataStream&, SireIO::PDBAtom&);

QDataStream& operator<<(QDataStream&, const SireIO::PDB2&);
QDataStream& operator>>(QDataStream&, SireIO::PDB2&);

namespace SireIO
{

/** This class provides functionality for reading/writing
    Protein Data Bank (PDB) ATOM records.

    @author Lester Hedges
*/
class SIREIO_EXPORT PDBAtom
{

friend QDataStream& ::operator<<(QDataStream&, const PDBAtom&);
friend QDataStream& ::operator>>(QDataStream&, PDBAtom&);

public:
    /** Default constructor. */
    PDBAtom();

    /** Constructor. */
    PDBAtom(const QString &line, QStringList &errors);

    /** Constructor. */
    PDBAtom(const SireMol::Atom &atom, bool is_ter, QStringList &errors);

    /** Generate a PDB record from the atom data. */
    QString toPDBRecord() const;

    static const char* typeName();

    /** Convert the atom name to PDB format. */
    QString toPDBName() const;

    /** Set the terminal atom flag. */
    void setTerminal(bool is_ter);

    /** Get the atom serial number. */
    qint64 getNumber() const;

    /** Set the atom serial number. */
    void setSerial(int serial);

    /** Get the atom name. */
    QString getName() const;

    /** Get the residue name. */
    QString getResName() const;

    /** Get the chain id. */
    QChar getChainID() const;

    /** Set the chain id. */
    void setChainID(QChar id);

    /** Get the residue sequence number. */
    qint64 getResNum() const;

    /** Set the residue sequence number. */
    void setResNum(int num);

    /** Get the residue index. */
    qint64 getResIdx() const;

    /** Set the residue index. */
    void setResIdx(int idx);

    /** Get the residue insertion code. */
    QChar getInsertCode() const;

    /** Get the atom coordinates. */
    SireMaths::Vector getCoord() const;

    /** Get the occupancy. */
    double getOccupancy() const;

    /** Get the atom temperature factor. */
    double getTemperature() const;

    /** Set the atom temperature factor. */
    void setTemperature(double temperature);

    /** Get the element symbol. */
    QString getElement() const;

    /** Get the charge on the atom. */
    qint64 getCharge() const;

    /** Whether this is a HETATM. */
    bool isHet() const;

    /** Whether this is a terminal atom. */
    bool isTer() const;

private:
    /** The original PDB record used to instantiate the atom. */
    QString record;

    /** Serial number. */
    qint64 serial;

    /** Name. */
    QString name;

    /** Alternate location indicator. */
    QChar alt_loc;

    /** Residue name. */
    QString res_name;

    /** Chain ID. */
    QChar chain_id;

    /** Residue sequence number. */
    qint64 res_num;

    /** Residue index. */
    qint64 res_idx;

    /** Residue insertion code. */
    QChar insert_code;

    /** Coordinates. */
    SireMaths::Vector coord;

    /** Occupancy. */
    double occupancy;

    /** Temperature factor. */
    double temperature;

    /** Element symbol. */
    QString element;

    /** Charge on the atom. */
    qint64 charge;

    /** Whether the atom is recorded as a HETAM entry. */
    bool is_het;

    /** Whether this is the last atom in a chain. */
    bool is_ter;
};

/** This class holds a parser for reading and writing
    Protein Data Bank (PDB) files

    @author Lester Hedges
*/
class SIREIO_EXPORT PDB2 : public SireBase::ConcreteProperty<PDB2,MoleculeParser>
{

friend QDataStream& ::operator<<(QDataStream&, const PDB2&);
friend QDataStream& ::operator>>(QDataStream&, PDB2&);

public:
    PDB2();
    PDB2(const QString &filename,
         const PropertyMap &map = PropertyMap());

    PDB2(const QStringList &lines,
         const PropertyMap &map = PropertyMap());
    PDB2(const SireSystem::System &system,
         const PropertyMap &map = PropertyMap());

    PDB2(const PDB2 &other);

    ~PDB2();

    PDB2& operator=(const PDB2 &other);

    bool operator==(const PDB2 &other) const;
    bool operator!=(const PDB2 &other) const;

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

    bool isLead() const;
    bool canFollow() const;

    int nMolecules() const;
    int nAtoms() const;
    int nResidues() const;
    int nResidues(int i) const;
    int nChains() const;
    int nChains(int i) const;
    int nAtoms(int i) const;

protected:
    SireSystem::System startSystem(const PropertyMap &map) const;
    void addToSystem(SireSystem::System &system, const PropertyMap &map) const;

private:
    void assertSane() const;
    void parseLines(const PropertyMap &map);

    void parseMolecule(const SireMol::Molecule &sire_mol, QVector<QString> &atom_lines,
        QStringList &errors, const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SireMol::Molecule updateMolecule(const SireMol::Molecule &sire_mol, QVector<bool> &used_atoms,
        const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;

    void findAtom(const SireMol::Atom &sire_atom, int &mol_idx,
        int &atom_idx, QVector<bool> &used_atoms) const;

    SireMol::MolStructureEditor getMolStructure(int imol,
        const SireBase::PropertyName &cutting) const;

    SireMol::MolEditor getMolecule(int imol, const PropertyMap &map = PropertyMap()) const;

    bool isModel() const;
    bool isModel(const SireSystem::System &system) const;

    //* Atom record data for each molecule. */
    QVector<QVector<PDBAtom> > atoms;

    //* Mapping between chain identifiers and residue index for each molecule. */
    QVector<QMultiMap<QChar, qint64> > chains;

    //* Mapping between residue and atom indices for each molecule. */
    QVector<QMultiMap<qint64, qint64> > residues;

    /** The name of the file from which data was parsed. */
    QString filename;

    /** Any warnings that were raised when reading the file. */
    QStringList parse_warnings;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** The PDB2 parser is a lead parser - it is capable alone
    of creating the System */
inline bool PDB2::isLead() const
{
    return true;
}

/** The PDB2 parser can follow another a lead parser. */
inline bool PDB2::canFollow() const
{
    return true;
}

#endif // SIRE_SKIP_INLINE_FUNCTIONS

}

Q_DECLARE_METATYPE( SireIO::PDBAtom )
Q_DECLARE_METATYPE( SireIO::PDB2 )

SIRE_EXPOSE_CLASS( SireIO::PDB2 )

SIRE_END_HEADER

#endif
