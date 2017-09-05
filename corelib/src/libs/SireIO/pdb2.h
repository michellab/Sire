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

#include "SireMol/atomidx.h"

SIRE_BEGIN_HEADER

namespace SireIO
{
class PDBAtom;
class PDBCrystal;
class PDBHelix;
class PDBMaster;
class PDBSheet;
class PDBTitle;
class PDBTransform;
class PDB2;
}

namespace SireMol
{
class MolEditor;
class MoleculeInfoData;
}

QDataStream& operator<<(QDataStream&, const SireIO::PDBAtom&);
QDataStream& operator>>(QDataStream&, SireIO::PDBAtom&);

QDataStream& operator<<(QDataStream&, const SireIO::PDBCrystal&);
QDataStream& operator>>(QDataStream&, SireIO::PDBCrystal&);

QDataStream& operator<<(QDataStream&, const SireIO::PDBHelix&);
QDataStream& operator>>(QDataStream&, SireIO::PDBHelix&);

QDataStream& operator<<(QDataStream&, const SireIO::PDBMaster&);
QDataStream& operator>>(QDataStream&, SireIO::PDBMaster&);

QDataStream& operator<<(QDataStream&, const SireIO::PDBSheet&);
QDataStream& operator>>(QDataStream&, SireIO::PDBSheet&);

QDataStream& operator<<(QDataStream&, const SireIO::PDBTitle&);
QDataStream& operator>>(QDataStream&, SireIO::PDBTitle&);

QDataStream& operator<<(QDataStream&, const SireIO::PDBTransform&);
QDataStream& operator>>(QDataStream&, SireIO::PDBTransform&);

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

    /** Generate a PDB record from the atom data. */
    QString toPDBLine() const;

    /** Generate a string representation of the object. */
    QString toString() const;

    static const char* typeName();

    /** Set the terminal atom flag. */
    void setTerminal(bool _isTer);

    /** Get the atom serial number. */
    qint64 getSerial() const;

    /** Get the atom name. */
    QString getName() const;

    /** Get the residue name. */
    QString getResName() const;

    /** Get the chain id. */
    QString getChainId() const;

    /** Get the residue sequence number. */
    qint64 getResNum() const;

    /** Get the atom coordinates. */
    SireMaths::Vector getCoord() const;

    /** Get the occupancy. */
    double getOccupancy() const;

    /** Get the atom temperature factor. */
    double getTemperature() const;

    /** Get the element symbol. */
    QString getElement() const;

    /** Get the charge on the atom. */
    qint64 getCharge() const;

    /** Whether this is a terminal atom. */
    bool isTer() const;

    /** Set anisotropic temperature record data. */
    void setAnisTemp(const QString &line1, const QString &line2, QStringList& errors);

    /** Get the index of the atom within the molecule. */
    SireMol::AtomIdx getMolIndex() const;

    /** Set the index of the atom within the molecule. */
    void setMolIndex(SireMol::AtomIdx index);

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

    /** Whether there is an anisotropic temperature (ANISOU) record for the atom. */
    bool is_anis;

    /** Anisotropic temperature factors. */
    qint64 anis_facts[6];

    /** The index of the atom within the molecule. */
    SireMol::AtomIdx mol_index;
};

/** This class provides functionality for reading/writing
    Protein Data Bank (PDB) CRYST1 records.

    @author Lester Hedges
*/
class SIREIO_EXPORT PDBCrystal
{

friend QDataStream& ::operator<<(QDataStream&, const PDBCrystal&);
friend QDataStream& ::operator>>(QDataStream&, PDBCrystal&);

public:
    /** Default constructor. */
    PDBCrystal();

    /** Constructor. */
    PDBCrystal(const QString &line, QStringList &errors);

    /** Generate a PDB record from the crystallographic data. */
    QString toPDBLine() const;

    /** Generate a string representation of the object. */
    QString toString() const;

    /** Whether the object contains a record. */
    bool hasRecord() const;

    static const char* typeName();

private:
    /** The original PDB record used to instantiate the crsytal object. */
    QString record;

    /** Unit cell a, b, c base lengths. */
    SireMaths::Vector base_lengths;

    /** Unit cell alpha, beta, gamma angles. */
    SireMaths::Vector angles;

    /** Space group. */
    QString space_group;

    /** Z value. */
    qint64 z;
};

/** This class provides functionality for reading/writing
    Protein Data Bank (PDB) HELIX records.

    @author Lester Hedges
*/
class SIREIO_EXPORT PDBHelix
{

friend QDataStream& ::operator<<(QDataStream&, const PDBHelix&);
friend QDataStream& ::operator>>(QDataStream&, PDBHelix&);

public:
    /** Default constructor. */
    PDBHelix();

    /** Constructor. */
    PDBHelix(const QString &line, QStringList &errors);

    /** Generate a PDB record from the helix data. */
    QString toPDBLine() const;

    /** Generate a string representation of the object. */
    QString toString() const;

    static const char* typeName();

private:
    /** The original PDB record used to instantiate the helix. */
    QString record;

    /** Serial number. */
    qint64 serial;

    /** Helix ID. */
    QString id;

    /** Name of the initial residue. */
    QString init_res_name;

    /** ID for the initial chain containing the helix. */
    QChar init_chain_id;

    /** Sequence number of the initial residue. */
    qint64 init_res_num;

    /** Insertion code of the initial residue. */
    QChar init_insert_code;

    /** Name of the end residue. */
    QString end_res_name;

    /** ID for the end chain containing the helix. */
    QChar end_chain_id;

    /** Sequence number of the end residue. */
    qint64 end_res_num;

    /** Insertion code of the end residue. */
    QChar end_insert_code;

    /** Helix class. */
    qint64 helix_class;

    /** Comment. */
    QString comment;

    /** Length of the helix. */
    qint64 length;
};

/** This class provides functionality for reading/writing
    Protein Data Bank (PDB) MASTER records.

    If present, this record provides data for validating the parser.

    @author Lester Hedges
*/
class SIREIO_EXPORT PDBMaster
{

friend QDataStream& ::operator<<(QDataStream&, const PDBMaster&);
friend QDataStream& ::operator>>(QDataStream&, PDBMaster&);

public:
    /** Default constructor. */
    PDBMaster();

    /** Constructor. */
    PDBMaster(const QString &line, QStringList &errors);

    /** Generate a PDB record from the master data. */
    QString toPDBLine() const;

    /** Generate a string representation of the object. */
    QString toString() const;

    static const char* typeName();

    int nRemarks() const;
    int nHets() const;
    int nHelices() const;
    int nSheets() const;
    int nSites() const;
    int nTransforms() const;
    int nAtoms() const;
    int nTers() const;
    int nConnects() const;
    int nSequences() const;

private:
    /** The original PDB record used to instantiate the master object. */
    QString record;

    /** Number of REMARK records. */
    qint64 num_remarks;

    /** Number of HET records. */
    qint64 num_hets;

    /** Number of HELIX records. */
    qint64 num_helices;

    /** Number of SHEET records. */
    qint64 num_sheets;

    /** Number of SITE records. */
    qint64 num_sites;

    /** Number of tranformation records (ORIGX, SCALE, MTRIX). */
    qint64 num_transforms;

    /** Number of ATOM or HETATM records. */
    qint64 num_atoms;

    /** Number of TER records. */
    qint64 num_ters;

    /** Number of CONECT records. */
    qint64 num_connects;

    /** Number of SEQRES records. */
    qint64 num_sequences;
};

/** This class provides functionality for reading/writing
    Protein Data Bank (PDB) SHEET records.

    @author Lester Hedges
*/
class SIREIO_EXPORT PDBSheet
{

friend QDataStream& ::operator<<(QDataStream&, const PDBSheet&);
friend QDataStream& ::operator>>(QDataStream&, PDBSheet&);

public:
    /** Default constructor. */
    PDBSheet();

    /** Constructor. */
    PDBSheet(const QString &line, QStringList &errors);

    /** Generate a PDB record from the sheet data. */
    QString toPDBLine() const;

    /** Generate a string representation of the object. */
    QString toString() const;

    static const char* typeName();

private:
    /** The original PDB record used to instantiate the sheet. */
    QString record;

    /** Strand number. */
    qint64 strand;

    /** Sheet id. */
    QString id;

    /** Number of strands. */
    qint64 num_strands;

    /** Name of the initial residue. */
    QString init_res_name;

    /** ID for the initial chain containing the helix. */
    QChar init_chain_id;

    /** Sequence number of the initial residue. */
    qint64 init_res_num;

    /** Insertion code of the initial residue. */
    QChar init_insert_code;

    /** Name of the end residue. */
    QString end_res_name;

    /** ID for the terminal chain residue. */
    QChar end_chain_id;

    /** Sequence number of the end residue. */
    qint64 end_res_num;

    /** Insertion code of the end residue. */
    QChar end_insert_code;

    /** Sense of strand with respect to previous strand. */
    qint64 sense;

    /** Atom name in current strand. */
    QString curr_atm_name;

    /** Residue name in current strand. */
    QString curr_res_name;

    /** Chain ID in current strand. */
    QChar curr_chain_id;

    /** Sequence number in current strand. */
    qint64 curr_res_num;

    /** Insertion code in current strand. */
    QChar curr_insert_code;

    /** Atom name in previous strand. */
    QString prev_atm_name;

    /** Residue name in previous strand. */
    QString prev_res_name;

    /** Chain ID in previous strand. */
    QChar prev_chain_id;

    /** Sequence number in previous strand. */
    qint64 prev_res_num;

    /** Insertion code in previous strand. */
    QChar prev_insert_code;
};

/** This class provides functionality for reading/writing
    Protein Data Bank (PDB) title section records, e.g.
    HEADER, TITLE, REMARK, ...

    @author Lester Hedges
*/
class SIREIO_EXPORT PDBTitle
{

friend QDataStream& ::operator<<(QDataStream&, const PDBTitle&);
friend QDataStream& ::operator>>(QDataStream&, PDBTitle&);

public:
    enum RECORD_TYPE
    {
        HEADER = 0,
        OBSLTE = 1,
        TITLE  = 2,
        SPLIT  = 3,
        CAVEAT = 4,
        COMPND = 5,
        SOURCE = 6,
        KEYWDS = 7,
        EXPDTA = 8,
        NUMMDL = 9,
        MDLTYP = 10,
        AUTHOR = 11,
        REVDAT = 12,
        SPRSDE = 13,
        JRNL   = 14,
        REMARK = 15
    };

    /** Default constructor. */
    PDBTitle();

    /** Append a PDB record. */
    void appendRecord(const QString &line, RECORD_TYPE record_type, QStringList &errors);

    /** Generate PDB records from the atom data. */
    QStringList toPDBLines() const;

    /** Generate a string representation of the object. */
    QString toString() const;

    /** Return the number of title section records. */
    int nRecords() const;

    /** Return the number of REMARK records. */
    int nRemarks() const;

    /** Return the number of MODEL that should be in the PDB file. */
    int nModels() const;

    /** Whether the object contains any records. */
    bool hasRecords() const;

    static const char* typeName();

private:
    /** The original PDB records used to generate the title object. */
    QStringList records;

    /** Header record. */
    QString header;

    /** Obsolete records. */
    QStringList obsoletes;

    /** Title records. */
    QStringList titles;

    /** Split records. */
    QStringList splits;

    /** Caveat records. */
    QStringList caveats;

    /** Compound records. */
    QStringList compounds;

    /** Source records. */
    QStringList sources;

    /** Keyword records. */
    QStringList keywords;

    /** Experiment records. */
    QStringList experiments;

    /** The number of models. */
    qint64 num_models;

    /** Model type records. */
    QStringList model_types;

    /** Author records. */
    QStringList authors;

    /** Revision data records. */
    QStringList revisions;

    /** Superseded ID code records. */
    QStringList supersedes;

    /** Journal records. */
    QStringList journals;

    /** Remark records. */
    QStringList remarks;
};

/** This class provides functionality for reading/writing
    Protein Data Bank (PDB) ORIGXn transformation records.

    These records present the transformation from the orthogonal
    coordinates contained in the entry to the submitted coordinates.

    @author Lester Hedges
*/
class SIREIO_EXPORT PDBTransform
{

friend QDataStream& ::operator<<(QDataStream&, const PDBTransform&);
friend QDataStream& ::operator>>(QDataStream&, PDBTransform&);

public:
    /** Default constructor. */
    PDBTransform();

    /** Append a record to the transformation object. */
    void appendRecord(const QString &line, int dimension,
        bool isNonCryst, QStringList &errors);

    /** Generate a PDB record from the ORIGXn transformation data. */
    QStringList toPDBLines() const;

    /** Generate a string representation of the object. */
    QString toString() const;

    /** Whether this object contains a complete transformation record. */
    bool hasRecord() const;

    static const char* typeName();

private:
    /** The original PDB record used to instantiate the ORIGXn object. */
    QStringList records;

    /** Serial number (if present). */
    qint64 serial;

    /** Whether the coordinates in the PDB file are approximated
        by this transformation. */
    bool isGiven;

    /** Transformation vectors. */
    QVector<SireMaths::Vector> transforms;

    /** Offset vector. */
    SireMaths::Vector offsets;

    /** Whether there is a record for each dimension. */
    bool isDimension[3];
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

    QString formatName() const;
    QString formatDescription() const;
    QStringList formatSuffix() const;

    int nTitles() const;
    int nAtoms() const;
    int nHelices() const;
    int nSheets() const;
    int nTers() const;

    bool hasCrystal() const;
    bool hasMaster() const;
    bool hasTransOrig() const;
    bool hasTransScale() const;
    bool hasTransMatrix() const;

    PDBTitle getTitle() const;
    QVector<QVector<PDBAtom> > getAtoms() const;
    QVector<PDBHelix> getHelices() const;
    QVector<PDBSheet> getSheets() const;
    PDBMaster getMaster() const;
    QMap<qint64, QString> getInvalidRecords() const;

protected:
    SireSystem::System startSystem(int iframe, const PropertyMap &map);
    SireSystem::System startSystem(const PropertyMap &map);
    void addToSystem(SireSystem::System &system, const PropertyMap &map) const;

private:
    void assertSane() const;
    void parseLines(const PropertyMap &map);

    SireMol::MolStructureEditor getMolStructure(int imol,
        int iframe, const SireBase::PropertyName &cutting);

    SireMol::MolEditor getMolecule(int imol, int iframe,
        const PropertyMap &map = PropertyMap());

    /** Validate that the two atoms are the same (other than position). */
    bool validateAtom(const PDBAtom &atom1, const PDBAtom &atom2) const;

    //* Title record data. */
    PDBTitle title;

    //* Atom record data (possibly multiple frames). */
    QVector<QVector<PDBAtom> > atoms;

    //* Residue mapping. */
    QVector<QMultiMap<QPair<qint64, QString>, qint64> > residues;

    //* Connectivity data. */
    QVector<QMultiMap<qint64, qint64> > connections;

    //* Helix record data. */
    QVector<PDBHelix> helices;

    //* Sheet record data. */
    QVector<PDBSheet> sheets;

    //* Crystallographic record data. */
    PDBCrystal crystal;

    //* Transformation record: orthogonal to submitted coordinates. */
    PDBTransform trans_orig;

    //* Transformation record: orthogonal to fractional crystallographic. */
    PDBTransform trans_scale;

    //* Transformation record: non-crystallographic. */
    PDBTransform trans_matrix;

    //* Master record data. */
    PDBMaster master;

    /** Number of TER records. */
    qint64 num_ters;

    /** Whether a MASTER record was found (can be used for validation). */
    bool has_master;

    /** Invalid records */
    QMap<qint64, QString> invalid_records;

    /** Any warnings that were raised when reading the file */
    QStringList parse_warnings;
};

}

Q_DECLARE_METATYPE( SireIO::PDBAtom )
Q_DECLARE_METATYPE( SireIO::PDBCrystal )
Q_DECLARE_METATYPE( SireIO::PDBHelix )
Q_DECLARE_METATYPE( SireIO::PDBMaster )
Q_DECLARE_METATYPE( SireIO::PDBSheet )
Q_DECLARE_METATYPE( SireIO::PDBTitle )
Q_DECLARE_METATYPE( SireIO::PDBTransform )
Q_DECLARE_METATYPE( SireIO::PDB2 )

SIRE_EXPOSE_CLASS( SireIO::PDBAtom )
SIRE_EXPOSE_CLASS( SireIO::PDBCrystal )
SIRE_EXPOSE_CLASS( SireIO::PDBHelix )
SIRE_EXPOSE_CLASS( SireIO::PDBMaster )
SIRE_EXPOSE_CLASS( SireIO::PDBSheet )
SIRE_EXPOSE_CLASS( SireIO::PDBTitle )
SIRE_EXPOSE_CLASS( SireIO::PDBTransform )
SIRE_EXPOSE_CLASS( SireIO::PDB2 )

SIRE_END_HEADER

#endif
