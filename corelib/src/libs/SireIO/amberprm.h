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

#ifndef SIREIO_AMBERPRM_H
#define SIREIO_AMBERPRM_H

#include "moleculeparser.h"

#include "SireMM/mmdetail.h"

#include "SireMaths/vector.h"
#include "SireMM/ljparameter.h"

#include <QSet>
#include <QHash>

SIRE_BEGIN_HEADER

namespace SireIO
{
class AmberPrm;
class AmberRst7;
}

QDataStream& operator<<(QDataStream&, const SireIO::AmberPrm&);
QDataStream& operator>>(QDataStream&, SireIO::AmberPrm&);

namespace SireMol
{
class MolEditor;
class MoleculeInfoData;
}

namespace SireMM
{
class AmberParams;
}

namespace SireIO
{

/** This class represents an Amber-format parameter file, currently
    supporting top files produced from Amber7 until Amber16

    The format of this file is described here;

    http://ambermd.org/formats.html

    (specifically the "PARM" parameter/topology file specification)

    @author Christopher Woods
*/
class SIREIO_EXPORT AmberPrm : public SireBase::ConcreteProperty<AmberPrm,MoleculeParser>
{

friend QDataStream& ::operator<<(QDataStream&, const AmberPrm&);
friend QDataStream& ::operator>>(QDataStream&, AmberPrm&);

public:
    enum FLAG_TYPE { UNKNOWN = 0,
                     INTEGER = 1,
                     FLOAT = 2,
                     STRING = 3 };

    AmberPrm();

    AmberPrm(const QString &filename,
             const PropertyMap &map = PropertyMap());
    AmberPrm(const QStringList &lines,
             const PropertyMap &map = PropertyMap());

    AmberPrm(const SireSystem::System &system,
             const PropertyMap &map = PropertyMap());

    AmberPrm(const AmberPrm &other);

    ~AmberPrm();

    AmberPrm& operator=(const AmberPrm &other);

    bool operator==(const AmberPrm &other) const;
    bool operator!=(const AmberPrm &other) const;

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
    QStringList formatSuffix() const;

    QString formatDescription() const;

    bool isLead() const;
    bool canFollow() const;

    static AmberPrm parse(const QString &filename,
                          const PropertyMap &map = PropertyMap());

    SireMol::Molecule getMolecule(int i, const PropertyMap &map = PropertyMap()) const;

    SireMol::Molecule getMolecule(int i, const AmberRst7 &rst,
                                  const PropertyMap &map = PropertyMap()) const;

    QVector<QString> linesForFlag(const QString &flag) const;

    QStringList flags() const;

    FLAG_TYPE flagType(const QString &flag) const;

    QVector<qint64> intData(const QString &flag) const;
    QVector<double> floatData(const QString &flag) const;
    QVector<QString> stringData(const QString &flag) const;

    QString title() const;

    int nAtoms() const;

    int nTypes() const;

    int nBonds() const;
    int nBondsWithHydrogen() const;
    int nBondsNoHydrogen() const;

    int nAngles() const;
    int nAnglesWithHydrogen() const;
    int nAnglesNoHydrogen() const;

    int nDihedrals() const;
    int nDihedralsWithHydrogen() const;
    int nDihedralsNoHydrogen() const;

    int nExcluded() const;
    int nResidues() const;

    int nMolecules() const;

    int nAtoms(int molidx) const;

    void assertSane() const;

    SireMM::MMDetail forcefield() const;

protected:
    SireSystem::System startSystem(const PropertyMap &map) const;

private:
    void parse(const PropertyMap &map);

    void rebuildAfterReload();
    void rebuildLJParameters();
    void rebuildBADIndicies();
    void rebuildExcludedAtoms();
    void rebuildMolNumToAtomNums();

    SireMM::AmberParams getAmberParams(int imol, const SireMol::MoleculeInfoData &molinfo) const;

    SireMol::MolStructureEditor getMolStructure(int molidx,
                                                const SireBase::PropertyName &cutting) const;

    SireMol::MolEditor getMoleculeEditor(int molidx,
                                         const PropertyMap &map) const;

    QVector<int> getAtomIndexToMolIndex() const;

    /** Function to process all flags, returning the parsing score */
    double processAllFlags();

    /** A map showing the line number of all flags. This holds
        the start index and number of lines for each flag */
    QHash< QString,QPair<qint64,qint64> > flag_to_line;

    /** The raw int data for the integer flags */
    QHash< QString, QVector<qint64> > int_data;

    /** The raw float data for the float flags */
    QHash< QString, QVector<double> > float_data;

    /** The raw string data for the string flags */
    QHash< QString, QVector<QString> > string_data;

    /** All of the LJ parameters, indexed by atom type */
    QVector<SireMM::LJParameter> lj_data;

    /** The indicies of the bonds for each molecule */
    QVector< QVector<int> > bonds_inc_h, bonds_exc_h;

    /** The indicies of the angles for each molecule */
    QVector< QVector<int> > angs_inc_h, angs_exc_h;

    /** The indicies of the dihedrals for each molecule */
    QVector< QVector<int> > dihs_inc_h, dihs_exc_h;

    /** The excluded atoms for each atom of each molecule */
    QVector< QVector< QVector<int> > > excl_atoms;

    /** The AtomNums of each atom in each molecule (indexed by MolNum) */
    QVector< QVector<int> > molnum_to_atomnums;

    /** A copy of the POINTER data to prevent over-lookup */
    QVector<qint64> pointers;
    
    /** The forcefield for the molecules in this file */
    SireMM::MMDetail ffield;
};

}

Q_DECLARE_METATYPE( SireIO::AmberPrm )

SIRE_EXPOSE_CLASS( SireIO::AmberPrm )

SIRE_END_HEADER

#endif
