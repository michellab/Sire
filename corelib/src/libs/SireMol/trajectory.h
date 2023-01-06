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

#include "SireVol/space.h"

#include "SireBase/sharedpolypointer.hpp"
#include "SireBase/refcountdata.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class Trajectory;
class Frame;
class TrajectoryData;
class MolTrajectoryData;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::Trajectory&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::Trajectory&);

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::Frame&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::Frame&);

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::TrajectoryData&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::TrajectoryData&);

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::MolTrajectoryData&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::MolTrajectoryData&);

namespace SireMol
{

typedef SireBase::SharedPolyPointer<TrajectoryData> TrajectoryDataPtr;

using SireMaths::Vector;

/** This is a single trajectory frame. */
class SIREMOL_EXPORT Frame
    : public SireBase::ConcreteProperty<Frame,MoleculeProperty>
{

friend QDataStream& ::operator<<(QDataStream&, const Frame&);
friend QDataStream& ::operator>>(QDataStream&, Frame&);

public:
    Frame();

    Frame(const Molecule &mol, const PropertyMap &map = PropertyMap());
    Frame(const MoleculeData &mol, const PropertyMap &map = PropertyMap());

    Frame(const QVector<Vector> &coordinates,
          const Space &space,
          SireUnits::Dimension::Time time);

    Frame(const QVector<Vector> &coordinates,
          const QVector<Velocity3D> &velocities,
          const Space &space,
          SireUnits::Dimension::Time time);

    Frame(const QVector<Vector> &coordinates,
          const QVector<Velocity3D> &velocites,
          const QVector<Force3D> &forces,
          const Space &space,
          SireUnits::Dimension::Time time);

    Frame(const Frame &other);

    ~Frame();

    Frame& operator=(const Frame &other);

    bool operator==(const Frame &other) const;
    bool operator!=(const Frame &other) const;

    static const char* typeName();

    const char* what() const;

    Frame* clone() const;

    QString toString() const;

    bool isEmpty() const;

    bool hasCoordinates() const;
    bool hasVelocities() const;
    bool hasForces() const;

    QVector<Vector> coordinates() const;
    QVector<Velocity3D> velocities() const;
    QVector<Force3D> forces() const;

    const SireVol::Space& space() const;
    SireUnits::Dimension::Time time() const;

    int nAtoms() const;

    int numBytes() const;

    Frame subset(int start_atom, int natoms) const;

    bool isCompatibleWith(const MoleculeInfoData &molinfo) const;

protected:
    void assertSane() const;

private:
    QVector<Vector> coords;
    QVector<Velocity3D> vels;
    QVector<Force3D> frcs;
    SireVol::SpacePtr spc;
    SireUnits::Dimension::Time t;
};

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

    bool operator==(const TrajectoryData &other) const;
    bool operator!=(const TrajectoryData &other) const;

    virtual int nFrames() const=0;
    virtual int nAtoms() const=0;

    int nFiles() const;

    bool isEmpty() const;

    virtual QStringList filenames() const=0;

    virtual Frame getFrame(int i) const=0;

    QList<Frame> getFrames() const;
    QList<Frame> getFrames(int start_atom, int natoms) const;

    Frame operator[](int i) const;

    virtual bool isEditable() const;
    virtual void assertIsEditable() const;

    virtual TrajectoryDataPtr makeEditable() const;
    virtual TrajectoryDataPtr makeSubsetEditable(int start_atom, int natoms) const;

    virtual void setFrame(int i, const Frame &frame);
    virtual void appendFrame(const Frame &frame);
    virtual void insertFrame(int i, const Frame &frame);
    virtual void deleteFrame(int i);

protected:
    TrajectoryData& operator=(const TrajectoryData &other);

    virtual bool _equals(const TrajectoryData &other) const=0;
};

/** This class holds an editable trajectory that has been extracted
 *  for a single molecule. This is what is used when parts of a trajectory
 *  are edited for a single molecule (or when new trajectory data is
 *  generated for a new molecule)
 */
class SIREMOL_EXPORT MolTrajectoryData : public TrajectoryData
{

friend QDataStream& ::operator<<(QDataStream&, const MolTrajectoryData&);
friend QDataStream& ::operator>>(QDataStream&, MolTrajectoryData&);

public:
    MolTrajectoryData();
    MolTrajectoryData(const Frame &frame);
    MolTrajectoryData(const QList<Frame> &frames);

    MolTrajectoryData(const MolTrajectoryData &other);

    ~MolTrajectoryData();

    const char* what() const;

    static const char* typeName();

    MolTrajectoryData& operator=(const MolTrajectoryData &other);

    MolTrajectoryData* clone() const;

    int nFrames() const;
    int nAtoms() const;

    QStringList filenames() const;

    Frame getFrame(int i) const;

    bool isEditable() const;

    TrajectoryDataPtr makeEditable() const;
    TrajectoryDataPtr makeSubsetEditable(int start_atom, int natoms) const;

    void setFrame(int i, const Frame &frame);
    void appendFrame(const Frame &frame);
    void insertFrame(int i, const Frame &frame);
    void deleteFrame(int i);

protected:
    bool _equals(const TrajectoryData &other) const;

private:
    QList<Frame> frames;
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

    Trajectory(const TrajectoryData &trajectory);
    Trajectory(const QList<TrajectoryDataPtr> &trajectories);

    Trajectory(const TrajectoryData &trajectory,
               int start_atom, int natoms);

    Trajectory(const QList<TrajectoryDataPtr> &trajectories,
               int start_atom, int natoms);

    Trajectory(const Trajectory &other);

    ~Trajectory();

    Trajectory& operator=(const Trajectory &other);

    bool operator==(const Trajectory &other) const;
    bool operator!=(const Trajectory &other) const;

    static const char* typeName();

    QString toString() const;

    const char* what() const;

    Trajectory* clone() const;

    bool isEmpty() const;

    int nFrames() const;

    int count() const;
    int size() const;

    int nAtoms() const;

    Frame getFrame(int i) const;

    Frame operator[](int i) const;
    QList<Frame> operator[](const QList<qint64> &idxs) const;
    QList<Frame> operator[](const SireBase::Slice &slice) const;

    void setFrame(int i, const Frame &frame);
    void appendFrame(const Frame &frame);
    void insertFrame(int i, const Frame &frame);
    void deleteFrame(int i);

    bool isCompatibleWith(const MoleculeInfoData &molinfo) const;

private:
    int _getIndexForFrame(int &frame) const;

    TrajectoryData& _makeEditable(int &frame);
    Frame _subset(const Frame &frame) const;

    /** Handles to the underlying trajectory data */
    QList<TrajectoryDataPtr> d;

    /** Index of the molecule's first atom in the trajectory */
    qint64 start_atom;

    /** The number of atoms in this molecule */
    qint64 natoms;
};

}

Q_DECLARE_METATYPE(SireMol::Trajectory)
Q_DECLARE_METATYPE(SireMol::Frame)
Q_DECLARE_METATYPE(SireMol::MolTrajectoryData)

SIRE_EXPOSE_CLASS(SireMol::Trajectory)
SIRE_EXPOSE_CLASS(SireMol::Frame)

SIRE_END_HEADER

#endif
