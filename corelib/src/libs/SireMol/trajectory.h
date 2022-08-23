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

#ifndef SIREMOL_TRAJECTORY_H
#define SIREMOL_TRAJECTORY_H

#include "SireMol/molviewproperty.h"
#include "SireMol/atomcoords.h"
#include "SireMol/atomvelocities.h"
#include "SireMol/atomforces.h"

#include "SireBase/sharedpolypointer.hpp"
#include "SireBase/refcountdata.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class Trajectory;
class TrajectoryData;
}

QDataStream& operator<<(QDataStream&, const SireMol::Trajectory&);
QDataStream& operator>>(QDataStream&, SireMol::Trajectory&);

QDataStream& operator<<(QDataStream&, const SireMol::TrajectoryData&);
QDataStream& operator>>(QDataStream&, SireMol::TrajectoryData&);

namespace SireVol
{
class Space;
}

namespace SireMol
{

/** This is the virtual base class of all TrajectoryData objects.
 *  These objects hold handles to trajectory data. They do
 *  all of the work of loading trajectory data from the file formats
 *  they handle.
 */
class SIREMOL_EXPORT TrajectoryData : public SireBase::RefCountData
{

friend QDataStream& ::operator<<(QDataStream&, const TrajectoryData&);
friend QDataStream& ::operator>>(QDataStream&, TrajectoryData&);

public:
    typedef TrajectoryData ROOT;

    TrajectoryData();
    TrajectoryData(const TrajectoryData &other);

    virtual ~TrajectoryData();

    virtual const char* what() const=0;

    static const char* typeName();

    virtual TrajectoryData* clone() const=0;

    virtual int nFrames() const=0;
    virtual SireUnits::Dimension::Time timestep() const=0;

    virtual bool containsSpace() const=0;
    virtual bool containsCoordinates() const=0;
    virtual bool containsVelocities() const=0;
    virtual bool containsForces() const=0;

    virtual bool supportsSpace() const=0;
    virtual bool supportsCoordinates() const=0;
    virtual bool supportsVelocities() const=0;
    virtual bool supportsForces() const=0;

    void assertContainsSpace() const;
    void assertContainsCoordinates() const;
    void assertContainsVelocities() const;
    void assertContainsForces() const;

    void assertSupportsSpace() const;
    void assertSupportsCoordinates() const;
    void assertSupportsVelocities() const;
    void assertSupportsForces() const;

    virtual void assertValidIndex(qint64 index) const=0;

    int nFiles() const;

    bool isEmpty() const;

    virtual QStringList filenames() const=0;

    virtual const SireVol::Space& getSpace(int i) const=0;
    virtual AtomCoords getFrame(int i, const AtomCoords &coords, qint64 idx) const=0;
    virtual AtomVelocities getFrame(int i, const AtomVelocities &velocities, qint64 idx) const=0;
    virtual AtomForces getFrame(int i, const AtomForces &forces, qint64 idx) const=0;

    virtual void setFrame(int i, const AtomCoords &coords, qint64 idx)=0;
    virtual void setFrame(int i, const AtomVelocities &velocities, qint64 idx)=0;
    virtual void setFrame(int i, const AtomForces &forces, qint64 idx)=0;

    virtual void appendFrame(const AtomCoords &coords, qint64 idx);
    virtual void appendFrame(const AtomVelocities &velocities, qint64 idx);
    virtual void appendFrame(const AtomForces &forces, qint64 idx);

    virtual void insertFrame(int i, const AtomCoords &coords, qint64 idx);
    virtual void insertFrame(int i, const AtomVelocities &velocities, qint64 idx);
    virtual void insertFrame(int i, const AtomForces &forces, qint64 idx);

    virtual void deleteFrame(int i, qint64 idx)=0;

    virtual void insertFrame(int i)=0;
    virtual void deleteFrame(int i)=0;

protected:
    TrajectoryData& operator=(const TrajectoryData &other);
};


/** This is a molecular property that holds the handle to the
 *  trajectory data for that molecule. In addition to the
 *  handle, this also holds the index of the first atom
 *  in the underlying trajectory data (trajectory data is a
 *  vector of coordinates in atom index order for each molecule)
 */
class SIREMOL_EXPORT Trajectory
    : public SireBase::ConcreteProperty<Trajectory, SireMol::MoleculeProperty>
{

friend QDataStream& ::operator<<(QDataStream&, const Trajectory&);
friend QDataStream& ::operator>>(QDataStream&, Trajectory&);

public:
    Trajectory();
    Trajectory(const TrajectoryData &data, qint64 index);

    Trajectory(const Trajectory &other);

    ~Trajectory();

    Trajectory& operator=(const Trajectory &other);

    bool operator==(const Trajectory &other) const;
    bool operator!=(const Trajectory &other) const;

    static const char* typeName();

    const char* what() const;

    Trajectory* clone() const;

    bool isEmpty() const;

    int nFrames() const;
    SireUnits::Dimension::Time timestep() const;

    const SireVol::Space& getSpace(int i) const;
    AtomCoords getFrame(int i, const AtomCoords &coords) const;
    AtomVelocities getFrame(int i, const AtomVelocities &velocities) const;
    AtomForces getFrame(int i, const AtomForces &forces) const;

    void setFrame(int i, const AtomCoords &coords);
    void setFrame(int i, const AtomVelocities &velocities);
    void setFrame(int i, const AtomForces &forces);

    void appendFrame(const AtomCoords &coords);
    void appendFrame(const AtomVelocities &velocities);
    void appendFrame(const AtomForces &forces);

    void insertFrame(int i, const AtomCoords &coords);
    void insertFrame(int i, const AtomVelocities &velocities);
    void insertFrame(int i, const AtomForces &forces);

    void deleteFrame(int i);

    bool isCompatibleWith(const MoleculeInfoData &molinfo) const;

    bool containsSpace() const;
    bool containsCoordinates() const;
    bool containsVelocities() const;
    bool containsForces() const;

    bool supportsSpace() const;
    bool supportsCoordinates() const;
    bool supportsVelocities() const;
    bool supportsForces() const;

    void assertContainsSpace() const;
    void assertContainsCoordinates() const;
    void assertContainsVelocities() const;
    void assertContainsForces() const;

    void assertSupportsSpace() const;
    void assertSupportsCoordinates() const;
    void assertSupportsVelocities() const;
    void assertSupportsForces() const;

private:
    /** Handle to the underlying trajectory data */
    SireBase::SharedPolyPointer<TrajectoryData> d;

    /** Index of the molecule's first atom in the trajectory */
    qint64 idx;
};

}

Q_DECLARE_METATYPE(SireMol::Trajectory)

SIRE_EXPOSE_CLASS(SireMol::Trajectory)

SIRE_END_HEADER

#endif