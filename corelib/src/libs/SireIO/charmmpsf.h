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

#ifndef SIREIO_CHARMMPSF_H
#define SIREIO_CHARMMPSF_H

#include "moleculeparser.h"

#include "SireMaths/vector.h"

SIRE_BEGIN_HEADER

namespace SireIO
{
class PSFAtom;
class CharmmPSF;
}

namespace SireMol
{
class Atom;
class MolEditor;
class MoleculeInfoData;
class MoleculeView;
class Residue;
}

QDataStream& operator<<(QDataStream&, const SireIO::PSFAtom&);
QDataStream& operator>>(QDataStream&, SireIO::PSFAtom&);

QDataStream& operator<<(QDataStream&, const SireIO::CharmmPSF&);
QDataStream& operator>>(QDataStream&, SireIO::CharmmPSF&);

namespace SireIO
{

/** This class provides functionality for reading/writing
    CHARMM PSF atom records (!NATOM).

    @author Lester Hedges
*/
class SIREIO_EXPORT PSFAtom
{

friend QDataStream& ::operator<<(QDataStream&, const PSFAtom&);
friend QDataStream& ::operator>>(QDataStream&, PSFAtom&);

public:
    /** Default constructor. */
    PSFAtom();

    /** Constructor. */
    PSFAtom(const QString &line, int index, QStringList &errors);

    /** Constructor. */
    PSFAtom(const SireMol::Atom &atom, bool is_ter, QStringList &errors);

    /** Generate a PSD record from the atom data. */
    QString toPSFRecord() const;

    /** Convert the atom name to PDB format. */
    QString toPDBName() const;

    static const char* typeName();

    /** Get the atom index. */
    int getIndex() const;

    /** Get the atom number. */
    int getNumber() const;

    /** Get the segment name. */
    QString getSegment() const;

    /** Get the residue number. */
    qint64 getResNum() const;

    /** Get the residue name. */
    QString getResName() const;

    /** Get the atom name. */
    QString getName() const;

    /** Get the atom type. */
    QString getType() const;

    /** Get the atom charge. */
    double getCharge() const;

    /** Get the atom mass. */
    double getMass() const;

private:
    /** The index in the atoms vector. */
    qint64 index;

    /** Serial number. */
    qint64 number;

    /** Segment name. */
    QString segment;

    /** Residue number. */
    qint64 res_num;

    /** Residue name. */
    QString res_name;

    /** Atom name. */
    QString name;

    /** Atom number. */
    QString type;

    /** Charge. */
    double charge;

    /** Mass. */
    double mass;
};

/** This class holds a parser for reading and writing CHARMM PSF files.

    @author Lester Hedges
*/
class SIREIO_EXPORT CharmmPSF : public SireBase::ConcreteProperty<CharmmPSF,MoleculeParser>
{

friend QDataStream& ::operator<<(QDataStream&, const CharmmPSF&);
friend QDataStream& ::operator>>(QDataStream&, CharmmPSF&);

public:
    CharmmPSF();
    CharmmPSF(const QString &filename,
         const PropertyMap &map = PropertyMap());

    CharmmPSF(const QStringList &lines,
         const PropertyMap &map = PropertyMap());
    CharmmPSF(const SireSystem::System &system,
         const PropertyMap &map = PropertyMap());

    CharmmPSF(const CharmmPSF &other);

    ~CharmmPSF();

    CharmmPSF& operator=(const CharmmPSF &other);

    bool operator==(const CharmmPSF &other) const;
    bool operator!=(const CharmmPSF &other) const;

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

    int nMolecules() const;
    int nAtoms() const;
    int nAtoms(int i) const;
    int nBonds() const;
    int nBonds(int i) const;
    int nAngles() const;
    int nDihedrals() const;
    int nImpropers() const;
    int nCrossTerms() const;

protected:
    SireSystem::System startSystem(const PropertyMap &map) const;
    void addToSystem(SireSystem::System &system, const PropertyMap &map) const;

private:
    void assertSane() const;
    void parseLines(const PropertyMap &map);

    SireMol::MolStructureEditor getMolStructure(int imol,
        const SireBase::PropertyName &cutting) const;

    SireMol::MolEditor getMolecule(int imol,
        const PropertyMap &map = PropertyMap()) const;

    void findMolecules();

    void findBondedAtoms(int atom_idx, int mol_idx, const QHash<int, int> &bonded_atoms,
        QHash<int, int> &atom_to_mol, QSet<qint64> &atoms_in_mol,
        QVector<QPair<qint64, qint64> >& bonds_in_mol) const;

    /** The indices of the atoms in each molecule. */
    QVector<QVector<qint64> > molecules;

    /** Pairs of bonds (atom-to-atom) for each molecule. */
    QVector<QVector<QPair<qint64, qint64> > > molecule_bonds;

    /** The atom record data (!NATOM). */
    QVector<PSFAtom> atoms;

    /** The bond record data (!NBOND). */
    QVector<QVector<qint64> > bonds;

    /** The angle record data (!NTHETA). */
    QVector<QVector<qint64> > angles;

    /** The dihedral record data (!NPHI). */
    QVector<QVector<qint64> > dihedrals;

    /** The improper record data (!NIMPHI). */
    QVector<QVector<qint64> > impropers;

    /** The cross term record data (!NCRTERM). */
    QVector<QVector<qint64> > cross_terms;

    /** The atomic coordinates. */
    QVector<SireMaths::Vector> coords;

    /** The name of the parsed file (if from a file). */
    QString filename;

    /** Any warnings that were raised when reading the file. */
    QStringList parse_warnings;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** The CharmmPSF parser is a lead parser - it is capable alone
    of creating the System */
inline bool CharmmPSF::isLead() const
{
    return true;
}

#endif // SIRE_SKIP_INLINE_FUNCTIONS

}

Q_DECLARE_METATYPE( SireIO::PSFAtom )
Q_DECLARE_METATYPE( SireIO::CharmmPSF )

SIRE_EXPOSE_CLASS( SireIO::CharmmPSF )

SIRE_END_HEADER

#endif
