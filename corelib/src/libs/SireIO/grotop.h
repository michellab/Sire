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

#ifndef SIREIO_GROTOP_H
#define SIREIO_GROTOP_H

#include "moleculeparser.h"

#include "SireMol/atomname.h"
#include "SireMol/resname.h"
#include "SireMol/bondid.h"
#include "SireMol/angleid.h"
#include "SireMol/dihedralid.h"

#include "SireMM/gromacsparams.h"
#include "SireMM/mmdetail.h"

#include <QMultiHash>

SIRE_BEGIN_HEADER

namespace SireIO
{
class GroTop;
class GroMolType;
class GroAtom;
class GroSystem;
}

QDataStream& operator<<(QDataStream&, const SireIO::GroTop&);
QDataStream& operator>>(QDataStream&, SireIO::GroTop&);

QDataStream& operator<<(QDataStream&, const SireIO::GroMolType&);
QDataStream& operator>>(QDataStream&, SireIO::GroMolType&);

QDataStream& operator<<(QDataStream&, const SireIO::GroAtom&);
QDataStream& operator>>(QDataStream&, SireIO::GroAtom&);

QDataStream& operator<<(QDataStream&, const SireIO::GroSystem&);
QDataStream& operator>>(QDataStream&, SireIO::GroSystem&);

namespace SireIO
{

using SireMM::GromacsBond;
using SireMM::GromacsAngle;
using SireMM::GromacsDihedral;

/** This class is used by GroTop to hold the intermediate representation of
    a Gromacs atom in a moleculetype

    @author Christopher Woods
*/
class SIREIO_EXPORT GroAtom
{

friend QDataStream& ::operator<<(QDataStream&, const GroAtom&);
friend QDataStream& ::operator>>(QDataStream&, GroAtom&);

public:
    GroAtom();

    GroAtom(const GroAtom &other);

    ~GroAtom();

    GroAtom& operator=(const GroAtom &other);

    bool operator==(const GroAtom &other) const;
    bool operator!=(const GroAtom &other) const;

    static const char* typeName();
    const char* what() const;

    QString toString() const;

    bool isNull() const;

    SireMol::AtomName name() const;
    SireMol::AtomNum number() const;

    SireMol::ResName residueName() const;
    SireMol::ResNum residueNumber() const;

    qint64 chargeGroup() const;

    QString atomType() const;
    QString bondType() const;

    SireUnits::Dimension::Charge charge() const;
    SireUnits::Dimension::MolarMass mass() const;

    void setName(const QString &name);
    void setNumber(qint64 number);

    void setResidueName(const QString &name);
    void setResidueNumber(qint64 number);

    void setChargeGroup(qint64 grp);

    void setAtomType(const QString &atomtype);
    void setBondType(const QString &bondtype);
    
    void setCharge(SireUnits::Dimension::Charge charge);
    void setMass(SireUnits::Dimension::MolarMass mass);

private:
    /** Name of the atom */
    QString atmname;

    /** Name of the residue */
    QString resname;

    /** Atom type */
    QString atmtyp;
    
    /** Bond type - normally the same as the atom type */
    QString bndtyp;

    /** Atom number */
    qint64 atmnum;

    /** Residue number */
    qint64 resnum;

    /** Charge group */
    qint64 chggrp;

    /** Charge */
    SireUnits::Dimension::Charge chg;

    /** Mass */
    SireUnits::Dimension::MolarMass mss;
};

/** This class is used by GroTop to hold an intermediate representation of a
    Gromacs moleculetype. This provides metadata about the molecule that is
    needed to construct the whole molecule.

    @author Christopher Woods
*/
class SIREIO_EXPORT GroMolType
{

friend QDataStream& ::operator<<(QDataStream&, const GroMolType&);
friend QDataStream& ::operator>>(QDataStream&, GroMolType&);

public:
    GroMolType();
    GroMolType(const SireMol::Molecule &mol, const PropertyMap &map=PropertyMap());

    GroMolType(const GroMolType &other);

    ~GroMolType();

    GroMolType& operator=(const GroMolType &other);

    bool operator==(const GroMolType &other) const;
    bool operator!=(const GroMolType &other) const;

    static const char* typeName();
    const char* what() const;

    QString toString() const;

    bool isNull() const;

    QString name() const;
    void setName(const QString &name);

    qint64 nExcludedAtoms() const;
    void setNExcludedAtoms(qint64 nexcl);

    SireMM::MMDetail forcefield() const;

    void addAtom(const GroAtom &atom);

    void addBond(const SireMol::BondID &bond, const GromacsBond &parm);
    void addAngle(const SireMol::AngleID &angle, const GromacsAngle &parm);
    void addDihedral(const SireMol::DihedralID &dihedral, const GromacsDihedral &parm);

    void addBonds(const QMultiHash<SireMol::BondID,GromacsBond> &bonds);
    void addAngles(const QMultiHash<SireMol::AngleID,GromacsAngle> &angles);
    void addDihedrals(const QMultiHash<SireMol::DihedralID,GromacsDihedral> &dihedrals);

    void sanitise(QString elecstyle, QString vdwstyle,
                  QString combrule, double elec14, double vdw14);

    void addWarning(const QString &warning);

    int nAtoms() const;
    int nResidues() const;

    GroAtom atom(const SireMol::AtomIdx &atomidx) const;
    GroAtom atom(const SireMol::AtomNum &atomnum) const;
    GroAtom atom(const SireMol::AtomName &atomnam) const;

    QVector<GroAtom> atoms() const;

    QVector<GroAtom> atoms(const SireMol::AtomName &atomnam) const;

    QVector<GroAtom> atoms(const SireMol::ResIdx &residx) const;
    QVector<GroAtom> atoms(const SireMol::ResNum &resnum) const;
    QVector<GroAtom> atoms(const SireMol::ResName &resnam) const;

    QMultiHash<SireMol::BondID,GromacsBond> bonds() const;
    QMultiHash<SireMol::AngleID,GromacsAngle> angles() const;
    QMultiHash<SireMol::DihedralID,GromacsDihedral> dihedrals() const;

    QStringList warnings() const;

    bool needsSanitising() const;

private:
    void _pvt_sanitise();

    /** The name of this moleculetype */
    QString nme;

    /** A set of warnings that are generated about this type */
    QStringList warns;

    /** Array of all of the atoms in this molecule */
    QVector<GroAtom> atms;

    /** Array giving the index of the first atom in each residue */
    QVector<qint64> first_atoms;

    /** Hash of all of the bonds */
    QMultiHash<SireMol::BondID,GromacsBond> bnds;

    /** Hash of all of the angles */
    QMultiHash<SireMol::AngleID,GromacsAngle> angs;

    /** Hash of all of the dihedrals */
    QMultiHash<SireMol::DihedralID,GromacsDihedral> dihs;

    /** The details about the forcefield used for this molecule */
    SireMM::MMDetail ffield;

    /** The number of excluded atoms */
    qint64 nexcl;
};

/** This class describes a Gromacs System */
class SIREIO_EXPORT GroSystem
{

friend QDataStream& ::operator<<(QDataStream&, const GroSystem&);
friend QDataStream& ::operator>>(QDataStream&, GroSystem&);

public:
    GroSystem();
    GroSystem(const QString &name);

    GroSystem(const GroSystem &other);

    ~GroSystem();

    GroSystem& operator=(const GroSystem &other);

    bool operator==(const GroSystem &other) const;
    bool operator!=(const GroSystem &other) const;

    QString operator[](int i) const;

    QString at(int i) const;

    int size() const;
    int count() const;
    int nMolecules() const;

    QStringList uniqueTypes() const;

    static const char* typeName();
    const char* what() const;

    QString name() const;
    void setName(QString name);

    QString toString() const;

    bool isNull() const;
    bool isEmpty() const;

    void add(QString moltype, int ncopies=1);

private:
    /** Name of the system */
    QString nme;

    /** The list of molecule types */
    QStringList moltypes;

    /** The number of each type of molecule */
    QList<qint64> nmols;

    /** The total number of molecules */
    qint64 total_nmols;
};

/** This class holds a parser for reading and writing Gromacs "top" topology files.

    @author Christopher Woods
*/
class SIREIO_EXPORT GroTop : public SireBase::ConcreteProperty<GroTop,MoleculeParser>
{

friend QDataStream& ::operator<<(QDataStream&, const GroTop&);
friend QDataStream& ::operator>>(QDataStream&, GroTop&);

public:
    GroTop();
    GroTop(const QString &filename,
           const PropertyMap &map = PropertyMap());

    GroTop(const QStringList &lines,
           const PropertyMap &map = PropertyMap());
    GroTop(const SireSystem::System &system,
           const PropertyMap &map = PropertyMap());

    GroTop(const GroTop &other);

    ~GroTop();

    GroTop& operator=(const GroTop &other);

    bool operator==(const GroTop &other) const;
    bool operator!=(const GroTop &other) const;

    static const char* typeName();

    const char* what() const;

    bool isLead() const;
    bool canFollow() const;

    QStringList includePath(bool absolute_paths=false) const;
    QStringList includedFiles(bool absolute_paths=false) const;

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

    int nonBondedFunctionType() const;
    int combiningRules() const;
    double fudgeLJ() const;
    double fudgeQQ() const;
    bool generateNonBondedPairs() const;

    SireMM::GromacsAtomType atomType(const QString &atm) const;

    SireMM::GromacsBond bond(const QString &atm0, const QString &atm1, int func) const;
    QList<SireMM::GromacsBond> bonds(const QString &atm0, const QString &atm1, int func) const;

    SireMM::GromacsAngle angle(const QString &atm0, const QString &atm1,
                               const QString &atm2, int func) const;

    QList<SireMM::GromacsAngle> angles(const QString &atm0, const QString &atm1,
                                       const QString &atm2, int func) const;

    SireMM::GromacsDihedral dihedral(const QString &atm0, const QString &atm1,
                                     const QString &atm2, const QString &atm3,
                                     int func) const;

    QList<SireMM::GromacsDihedral> dihedrals(const QString &atm0, const QString &atm1,
                                             const QString &atm2, const QString &atm3,
                                             int func) const;

    QHash<QString,SireMM::GromacsAtomType> atomTypes() const;

    QMultiHash<QString,SireMM::GromacsBond> bondPotentials() const;
    QMultiHash<QString,SireMM::GromacsAngle> anglePotentials() const;
    QMultiHash<QString,SireMM::GromacsDihedral> dihedralPotentials() const;

    GroMolType moleculeType(const QString &name) const;
    QVector<GroMolType> moleculeTypes() const;

    GroSystem groSystem() const;

    QStringList postprocessedLines() const;

    QStringList warnings() const;

protected:
    SireSystem::System startSystem(const PropertyMap &map) const;

private:
    void assertSane() const;
    void parseLines(const QString &path, const PropertyMap &map);

    void getIncludePath(const PropertyMap &map);

    QString searchForDihType(const QString &atom0, const QString &atom1,
                             const QString &atom2, const QString &atom3,
                             int func) const;

    QString findIncludeFile(QString filename, QString current_directory);

    QVector<QString> loadInclude(QString filename, QString current_directory);

    QVector<QString> preprocess(const QVector<QString> &lines,
                                QHash<QString,QString> defines,
                                const QString &current_directory,
                                const QString &parent_file);

    void interpret();

    QStringList processDirectives(const QMap<int,QString> &taglocs,
                                  const QHash<QString,int> &ntags);

    const QVector<QString>& expandedLines() const;
    
    SireMol::Molecule createMolecule(const GroMolType &moltype, QStringList &errors,
                                     const PropertyMap &map) const;
    SireMol::Molecule createMolecule(QString moltype, QStringList &errors,
                                     const PropertyMap &map) const;
    
    typedef std::tuple<SireBase::Properties,QStringList> PropsAndErrors;
    
    PropsAndErrors getAtomProperties(const SireMol::MoleculeInfo &molinfo,
                                     const GroMolType &moltype) const;
    PropsAndErrors getBondProperties(const SireMol::MoleculeInfo &molinfo,
                                     const GroMolType &moltype) const;
    PropsAndErrors getAngleProperties(const SireMol::MoleculeInfo &molinfo,
                                      const GroMolType &moltype) const;
    PropsAndErrors getDihedralProperties(const SireMol::MoleculeInfo &molinfo,
                                         const GroMolType &moltype) const;
        
    /** This is the full search path of all directories that should
        be searched for Gromacs include files */
    QStringList include_path;

    /** This is the set of files that had to be included as part of parsing
        this file, arranged as dependencies of the files */
    QHash<QString,QStringList> included_files;

    /** The post-processed lines */
    QVector<QString> expanded_lines;

    /** The database of atom types */
    QHash<QString,SireMM::GromacsAtomType> atom_types;

    /** The database of bond potentials */
    QMultiHash<QString,SireMM::GromacsBond> bond_potentials;

    /** The database of angle potentials */
    QMultiHash<QString,SireMM::GromacsAngle> ang_potentials;

    /** The database of dihedral potentials */
    QMultiHash<QString,SireMM::GromacsDihedral> dih_potentials;

    /** The list of all moleculetypes loaded from this file */
    QVector<GroMolType> moltypes;

    /** The GroSystem that describes the system in full */
    GroSystem grosys;

    /** The non-bonded function type to use for all molecules */
    qint32 nb_func_type;

    /** The combining rule used for all molecules */
    qint32 combining_rule;

    /** The fudge LJ value for all molecules */
    double fudge_lj;

    /** The fudge QQ value for all molecules */
    double fudge_qq;

    /** All of the parse warnings */
    QStringList parse_warnings;

    /** Whether or not to generate pairs for all molecules */
    bool generate_pairs;
};

}

Q_DECLARE_METATYPE( SireIO::GroTop )
Q_DECLARE_METATYPE( SireIO::GroMolType )
Q_DECLARE_METATYPE( SireIO::GroAtom )
Q_DECLARE_METATYPE( SireIO::GroSystem )

SIRE_EXPOSE_CLASS( SireIO::GroTop )
SIRE_EXPOSE_CLASS( SireIO::GroMolType )
SIRE_EXPOSE_CLASS( SireIO::GroAtom )
SIRE_EXPOSE_CLASS( SireIO::GroSystem )

SIRE_END_HEADER

#endif
