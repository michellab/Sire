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

#ifndef SIREIO_DCDTRAJECTORY_H
#define SIREIO_DCDTRAJECTORY_H

#include "SireMol/trajectory.h"

SIRE_BEGIN_HEADER

namespace SireIO
{
class DCDTrajectory;
}

friend QDataStream& operator<<(QDataStream&, const SireIO::DCDTrajectory&);
friend QDataStream& operator>>(QDataStream&, SireIO::DCDTrajectory&);

namespace SireIO
{

/** This class provides a specialisation of TrajectoryData that can
 *  read coordinate data from a DCD file
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

    TrajectoryData* clone() const;

    int nFrames() const;
    SireUnits::Dimension::Time timestep() const;

    bool containsSpace() const;
    bool containsCoordinates() const;
    bool containsVelocities() const;
    bool containsForces() const;

    bool supportsSpace() const;
    bool supportsCoordinates() const;
    bool supportsVelocities() const;
    bool supportsForces() const;

    void assertValidIndex(qint64 index) const;

    int nFiles() const;

    bool isEmpty() const;

    QStringList filenames() const;

    const SireVol::Space& getSpace(int i) const;
    AtomCoords getFrame(int i, const AtomCoords &coords, qint64 idx) const;
    AtomVelocities getFrame(int i, const AtomVelocities &velocities, qint64 idx) const;
    AtomForces getFrame(int i, const AtomForces &forces, qint64 idx) const;
};

}

Q_DECLARE_METATYPE(SireIO::DCDTrajectory)

SIRE_END_HEADER

#endif
