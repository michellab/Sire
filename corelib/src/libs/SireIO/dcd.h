/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2022  Christopher Woods
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

#ifndef SIREIO_DCD_H
#define SIREIO_DCD_H

#include "moleculeparser.h"

#include "SireMol/trajectory.h"

#include "SireMaths/vector.h"

SIRE_BEGIN_HEADER

namespace SireIO
{
class DCD;
class DCDTrajectory;

namespace detail{ class DCDFile; }

}

SIREIO_EXPORT QDataStream& operator<<(QDataStream&, const SireIO::DCD&);
SIREIO_EXPORT QDataStream& operator>>(QDataStream&, SireIO::DCD&);

SIREIO_EXPORT QDataStream& operator<<(QDataStream&, const SireIO::DCDTrajectory&);
SIREIO_EXPORT QDataStream& operator>>(QDataStream&, SireIO::DCDTrajectory&);

namespace SireIO
{

class FortranFile;

namespace detail
{
/** This class provides a low-level interface to reading and writing
 *  a DCD file. It is designed to be used only with the
 *  DCD and DCDTrajectory classes
 */
class DCDFile
{
public:
    DCDFile();
    DCDFile(const QString &filename);
    ~DCDFile();

    void readHeader(FortranFile &file);

    QVector<double> readUnitCell(FortranFile &file, int frame) const;
    QVector<SireMaths::Vector> readCoordinates(FortranFile &file, int frame) const;

    QPair< QVector<double>,
           QVector<SireMaths::Vector> > readFrame(FortranFile &file, int frame) const;

    QString getTitle() const;

    double getTimeStep() const;
    qint64 getFrameStart() const;
    qint64 getFrameDelta() const;

    qint64 nFrames() const;
    qint64 nAtoms() const;

    double getTimeAtFrame(int frame) const;

private:
    QStringList title;

    QVector<qint32> fixed_atoms;
    QVector<SireMaths::Vector> first_frame;

    double timestep;

    qint64 istart;
    qint64 nsavc;
    qint64 nfixed;

    qint64 natoms;
    qint64 nframes;

    qint64 first_frame_line;

    bool CHARMM_FORMAT;
    bool HAS_EXTRA_BLOCK;
    bool HAS_FOUR_DIMS;
};
}

/** This class holds a complete DCD trajectory, and enables access
 *  to frames via the SireMol::Trajectory interface
*/
class SIREIO_EXPORT DCDTrajectory : public SireMol::TrajectoryData
{

friend QDataStream& ::operator<<(QDataStream&, const DCDTrajectory&);
friend QDataStream& ::operator>>(QDataStream&, DCDTrajectory&);

public:
    DCDTrajectory();
    DCDTrajectory(const QString &filename);
    DCDTrajectory(const DCDTrajectory &other);

    ~DCDTrajectory();

    const char* what() const;

    static const char* typeName();

    DCDTrajectory& operator=(const DCDTrajectory &other);

    DCDTrajectory* clone() const;

    int nFrames() const;
    int nAtoms() const;

    QStringList filenames() const;

    SireMol::Frame getFrame(int i) const;

protected:
    bool _equals(const TrajectoryData &other) const;

private:
    detail::DCDFile dcd;

    QString filename;
};

/** This class represents a DCD file reader. Note that this will
    only read the first frame from the file (it is for loading molecules).
    The whole trajectory will be added and parsed using
    the DCDTrajectory class.

    The format is described here

    https://www.ks.uiuc.edu/Research/vmd/plugins/molfile/dcdplugin.html

    @author Christopher Woods
*/
class SIREIO_EXPORT DCD : public SireBase::ConcreteProperty<DCD,MoleculeParser>
{

friend QDataStream& ::operator<<(QDataStream&, const DCD&);
friend QDataStream& ::operator>>(QDataStream&, DCD&);

public:
    DCD();
    DCD(const QString &filename,
        const PropertyMap &map = PropertyMap());
    DCD(const QStringList &lines,
        const PropertyMap &map = PropertyMap());
    DCD(const SireSystem::System &system,
        const PropertyMap &map = PropertyMap());

    DCD(const DCD &other);

    ~DCD();

    DCD& operator=(const DCD &other);

    bool operator==(const DCD &other) const;
    bool operator!=(const DCD &other) const;

    static const char* typeName();

    const char* what() const;

    MoleculeParserPtr construct(const QString &filename,
                                const PropertyMap &map) const;

    MoleculeParserPtr construct(const QStringList &lines,
                                const PropertyMap &map) const;

    MoleculeParserPtr construct(const SireSystem::System &system,
                                const PropertyMap &map) const;

    QString toString() const;

    bool isTextFile() const;

    QString formatName() const;
    QString formatDescription() const;
    QStringList formatSuffix() const;

    static DCD parse(const QString &filename);

    int nAtoms() const;

    bool hasCoordinates() const;
    bool hasVelocities() const;
    bool hasForces() const;

    QString title() const;
    double time() const;

    QVector<SireMaths::Vector> coordinates() const;
    QVector<SireMaths::Vector> velocities() const;
    QVector<SireMaths::Vector> forces() const;

    SireMaths::Vector boxDimensions() const;
    SireMaths::Vector boxAngles() const;

    QStringList warnings() const;

    void writeToFile(const QString &filename) const;

protected:
    void addToSystem(SireSystem::System &system, const PropertyMap &map) const;

private:
    void parse(const QString &filename, const PropertyMap &map);

    /** The title of the file */
    QString ttle;

    /** The current time of the simulation in picoseconds */
    double current_time;

    /** The coordinate data */
    QVector<SireMaths::Vector> coords;

    /** The velocity data */
    QVector<SireMaths::Vector> vels;

    /** The forces data */
    QVector<SireMaths::Vector> frcs;

    /** The box dimensions */
    SireMaths::Vector box_dims;

    /** The box angles */
    SireMaths::Vector box_angs;

    /** Any warnings that were raised when reading the file */
    QStringList parse_warnings;
};

}

Q_DECLARE_METATYPE( SireIO::DCD )
Q_DECLARE_METATYPE( SireIO::DCDTrajectory )

SIRE_EXPOSE_CLASS( SireIO::DCD )

SIRE_END_HEADER

#endif
