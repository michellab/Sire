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

#include "SireMM/gromacsparams.h"

SIRE_BEGIN_HEADER

namespace SireIO
{
class GroTop;
class GroMolType;
}

QDataStream& operator<<(QDataStream&, const SireIO::GroTop&);
QDataStream& operator>>(QDataStream&, SireIO::GroTop&);

QDataStream& operator<<(QDataStream&, const SireIO::GroMolType&);
QDataStream& operator>>(QDataStream&, SireIO::GroMolType&);

namespace SireIO
{

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
    
    GroMolType(const GroMolType &other);
    
    ~GroMolType();
    
    GroMolType& operator=(const GroMolType &other);
    
    bool operator==(const GroMolType &other) const;
    bool operator!=(const GroMolType &other) const;
    
    static const char* typeName();
    const char* what() const;
    
    QString toString() const;
    
    QString name() const;
    
    void addAtom(const QString &atomtype, int resnum, const QString &resname,
                 const QString &atomname, int chggroup, double chg, double mass);
    
    void checkSanity();
 
    void addWarning(const QString &warning);
    
    QStringList warnings() const;
    
private:
    /** The name of this moleculetype */
    QString nme;

    /** A set of warnings that are generated about this type */
    QStringList warns;
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

    SireMM::GromacsBond bond(const QString &atm0, const QString &atm1) const;
    QList<SireMM::GromacsBond> bonds(const QString &atm0, const QString &atm1) const;
    
    SireMM::GromacsAngle angle(const QString &atm0, const QString &atm1,
                               const QString &atm2) const;

    QList<SireMM::GromacsAngle> angles(const QString &atm0, const QString &atm1,
                                       const QString &atm2) const;

    SireMM::GromacsDihedral dihedral(const QString &atm0, const QString &atm1,
                                     const QString &atm2, const QString &atm3) const;

    QList<SireMM::GromacsDihedral> dihedrals(const QString &atm0, const QString &atm1,
                                             const QString &atm2, const QString &atm3) const;

    QHash<QString,SireMM::GromacsAtomType> atomTypes() const;
    
    QMultiHash<QString,SireMM::GromacsBond> bondPotentials() const;
    QMultiHash<QString,SireMM::GromacsAngle> anglePotentials() const;
    QMultiHash<QString,SireMM::GromacsDihedral> dihedralPotentials() const;

    GroMolType moleculeType(const QString &name) const;
    QVector<GroMolType> moleculeTypes() const;

protected:
    SireSystem::System startSystem(const PropertyMap &map) const;
    void addToSystem(SireSystem::System &system, const PropertyMap &map) const;

private:
    void assertSane() const;
    void parseLines(const QString &path, const PropertyMap &map);

    void getIncludePath(const PropertyMap &map);

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
    
    /** The non-bonded function type to use for all molecules */
    qint32 nb_func_type;
    
    /** The combining rule used for all molecules */
    qint32 combining_rule;
    
    /** The fudge LJ value for all molecules */
    double fudge_lj;
    
    /** The fudge QQ value for all molecules */
    double fudge_qq;
    
    /** Whether or not to generate pairs for all molecules */
    bool generate_pairs;
};

}

Q_DECLARE_METATYPE( SireIO::GroTop )
Q_DECLARE_METATYPE( SireIO::GroMolType )

SIRE_EXPOSE_CLASS( SireIO::GroTop )
SIRE_EXPOSE_CLASS( SireIO::GroMolType )

SIRE_END_HEADER

#endif
