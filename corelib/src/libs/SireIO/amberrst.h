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

#ifndef SIREIO_AMBERRST_H
#define SIREIO_AMBERRST_H

#include "moleculeparser.h"

#include "SireMaths/vector.h"

SIRE_BEGIN_HEADER

namespace SireIO
{
class AmberRst;
}

SIREIO_EXPORT QDataStream& operator<<(QDataStream&, const SireIO::AmberRst&);
SIREIO_EXPORT QDataStream& operator>>(QDataStream&, SireIO::AmberRst&);

namespace SireIO
{

class NetCDFFile;

/** This class represents an Amber-format restart/coordinate file (binary),
    currently supporting these files from Amber7 to Amber16.

    This is a netcdf file format, which is described here;

    http://ambermd.org/netcdf/nctraj.xhtml

    @author Christopher Woods
*/
class SIREIO_EXPORT AmberRst : public SireBase::ConcreteProperty<AmberRst,MoleculeParser>
{

friend SIREIO_EXPORT QDataStream& ::operator<<(QDataStream&, const AmberRst&);
friend SIREIO_EXPORT QDataStream& ::operator>>(QDataStream&, AmberRst&);

public:
    AmberRst();
    AmberRst(const QString &filename,
             const PropertyMap &map = PropertyMap());
    AmberRst(const QStringList &lines,
             const PropertyMap &map = PropertyMap());
    AmberRst(const SireSystem::System &system,
             const PropertyMap &map = PropertyMap());

    AmberRst(const AmberRst &other);

    ~AmberRst();

    AmberRst& operator=(const AmberRst &other);

    bool operator==(const AmberRst &other) const;
    bool operator!=(const AmberRst &other) const;

    static const char* typeName();

    const char* what() const;

    int count() const;
    int size() const;

    AmberRst operator[](int i) const;

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

    static AmberRst parse(const QString &filename);

    QString title() const;

    double time() const;
    double time(int frame) const;

    int nAtoms() const;

    bool hasCoordinates() const;
    bool hasVelocities() const;
    bool hasForces() const;

    QVector<SireMaths::Vector> coordinates() const;
    QVector<SireMaths::Vector> velocities() const;
    QVector<SireMaths::Vector> forces() const;

    bool isFrame() const;

    int nFrames() const;

    int nBytes() const;

    SireMol::Frame getFrame(int i) const;

    QVector<SireMaths::Vector> coordinates(int frame) const;
    QVector<SireMaths::Vector> velocities(int frame) const;
    QVector<SireMaths::Vector> forces(int frame) const;

    SireMaths::Vector boxDimensions() const;
    SireMaths::Vector boxAngles() const;

    SireMaths::Vector boxDimensions(int frame) const;
    SireMaths::Vector boxAngles(int frame) const;

    QStringList warnings() const;

    QString creatorApplication() const;

    double formatVersion() const;

    bool createdFromRestart() const;
    bool createdFromTrajectory() const;

    bool isTextFile() const;

    void writeToFile(const QString &filename) const;

protected:
    void addToSystem(SireSystem::System &system, const PropertyMap &map) const;

private:
    void assertSane() const;
    void parse(const NetCDFFile &netcdf, const PropertyMap &map);

    SireMol::Frame _getFrame(int i) const;

    /** The title of the file */
    QString ttle;

    /** The current time of the simulation in picoseconds */
    QVector<double> current_time;

    /** The coordinate data */
    QVector< QVector<SireMaths::Vector> > coords;

    /** The velocity data in amber units (angstrom / 1/20.455 picosecond) */
    QVector< QVector<SireMaths::Vector> > vels;

    /** The force data in amber units (amu*angstrom/picosecond^2) */
    QVector< QVector<SireMaths::Vector> > frcs;

    /** The box dimensions */
    QVector<SireMaths::Vector> box_dims;

    /** The box angles */
    QVector<SireMaths::Vector> box_angs;

    /** The version of the file format */
    double convention_version;

    /** The name and version of the application that wrote this file */
    QString creator_app;

    /** Any warnings that were raised when reading the file */
    QStringList parse_warnings;

    /** The number of frames in this file */
    qint64 nframes;

    /** Whether or not this was read as a restart file */
    bool created_from_restart;
};

}

Q_DECLARE_METATYPE( SireIO::AmberRst )

SIRE_EXPOSE_CLASS( SireIO::AmberRst )

SIRE_END_HEADER

#endif
