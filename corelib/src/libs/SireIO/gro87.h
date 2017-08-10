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

#ifndef SIREIO_GRO87_H
#define SIREIO_GRO87_H

#include "moleculeparser.h"

#include "SireMaths/vector.h"

SIRE_BEGIN_HEADER

namespace SireIO
{
class Gro87;
}

QDataStream& operator<<(QDataStream&, const SireIO::Gro87&);
QDataStream& operator>>(QDataStream&, SireIO::Gro87&);

namespace SireIO
{

/** This class holds a parser for reading and writing Gromacs Gro87 files.

    @author Christopher Woods
*/
class SIREIO_EXPORT Gro87 : public SireBase::ConcreteProperty<Gro87,MoleculeParser>
{

friend QDataStream& ::operator<<(QDataStream&, const Gro87&);
friend QDataStream& ::operator>>(QDataStream&, Gro87&);

public:
    Gro87();
    Gro87(const QString &filename,
          const PropertyMap &map = PropertyMap());

    Gro87(const QStringList &lines,
          const PropertyMap &map = PropertyMap());
    Gro87(const SireSystem::System &system,
          const PropertyMap &map = PropertyMap());
    
    Gro87(const Gro87 &other);
    
    ~Gro87();
    
    Gro87& operator=(const Gro87 &other);
    
    bool operator==(const Gro87 &other) const;
    bool operator!=(const Gro87 &other) const;
    
    static const char* typeName();
    
    const char* what() const;

    int count() const;
    int size() const;
    
    Gro87 operator[](int i) const;

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

    QString title() const;

    double time() const;
    double time(int frame) const;

    int nAtoms() const;
    int nResidues() const;
    
    bool hasCoordinates() const;
    bool hasVelocities() const;

    QVector<SireMaths::Vector> coordinates() const;
    QVector<SireMaths::Vector> velocities() const;

    QVector<qint64> atomNumbers() const;
    QVector<QString> atomNames() const;

    QVector<QString> residueNames() const;
    QVector<qint64> residueNumbers() const;

    int nFrames() const;
    
    QVector<SireMaths::Vector> coordinates(int frame) const;
    QVector<SireMaths::Vector> velocities(int frame) const;

    SireMaths::Vector boxV1() const;
    SireMaths::Vector boxV2() const;
    SireMaths::Vector boxV3() const;

    SireMaths::Vector boxV1(int frame) const;
    SireMaths::Vector boxV2(int frame) const;
    SireMaths::Vector boxV3(int frame) const;

    QStringList warnings() const;

protected:
    void addToSystem(SireSystem::System &system, const PropertyMap &map) const;

private:
    void assertSane() const;
    void parseLines(const PropertyMap &map);

    /** The title of the file */
    QString ttle;
    
    /** The current time of the simulation in picoseconds */
    QVector<double> current_time;
    
    /** The coordinate data in internal units (angstroms) */
    QVector< QVector<SireMaths::Vector> > coords;
    
    /** The velocity data in internal units (angstrom/picosecond) */
    QVector< QVector<SireMaths::Vector> > vels;
    
    /** The box dimensions - each of the three vectors */
    QVector<SireMaths::Vector> box_v1, box_v2, box_v3;
    
    /** The residue numbers */
    QVector<qint64> resnums;
    
    /** The residue names */
    QVector<QString> resnams;
    
    /** The atom names */
    QVector<QString> atmnams;
    
    /** The atom numbers */
    QVector<qint64> atmnums;
    
    /** Any warnings that were raised when reading the file */
    QStringList parse_warnings;
};

}

Q_DECLARE_METATYPE( SireIO::Gro87 )

SIRE_EXPOSE_CLASS( SireIO::Gro87 )

SIRE_END_HEADER

#endif
