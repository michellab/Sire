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

#include "trajectory.h"

#include "SireID/index.h"

#include "SireVol/space.h"

#include "SireMol/core.h"

#include "SireUnits/dimensions.h"
#include "SireUnits/units.h"

#include "SireBase/generalunitproperty.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMol;
using namespace SireVol;
using namespace SireBase;
using namespace SireUnits;
using namespace SireUnits::Dimension;
using namespace SireStream;

////////
//////// Implementation of TrajectoryData
////////

static const RegisterMetaType<TrajectoryData> r_trajdata(MAGIC_ONLY, "SireMol::TrajectoryData");

SIREMOL_EXPORT QDataStream& operator<<(QDataStream &ds, const TrajectoryData &trajdata)
{
    writeHeader(ds, r_trajdata, 1);
    return ds;
}

SIREMOL_EXPORT QDataStream& operator>>(QDataStream &ds, TrajectoryData &trajdata)
{
    VersionID v = readHeader(ds, r_trajdata);

    if (v != 1)
        throw version_error(v, "1", r_trajdata, CODELOC);

    return ds;
}

TrajectoryData::TrajectoryData() : RefCountData()
{}

TrajectoryData::TrajectoryData(const TrajectoryData&) : RefCountData()
{}

TrajectoryData::~TrajectoryData()
{}

TrajectoryData& TrajectoryData::operator=(const TrajectoryData&)
{
    return *this;
}

const char* TrajectoryData::typeName()
{
    return "SireMol::TrajectoryData";
}

int TrajectoryData::nFiles() const
{
    return this->filenames().count();
}

bool TrajectoryData::isEmpty() const
{
    return this->nFrames() == 0;
}

QList<Frame> TrajectoryData::getFrames() const
{
    QList<Frame> frames;

    const int n = this->nFrames();

    frames.reserve(n);

    for (int i=0; i<n; ++i)
    {
        frames.append(this->getFrame(i));
    }

    return frames;
}

QList<Frame> TrajectoryData::getFrames(int start_atom, int natoms) const
{
    QList<Frame> frames;

    const int n = this->nFrames();

    frames.reserve(n);

    for (int i=0; i<n; ++i)
    {
        frames.append(this->getFrame(i).subset(start_atom, natoms));
    }

    return frames;
}

Frame TrajectoryData::operator[](int i) const
{
    return this->getFrame(i);
}

bool TrajectoryData::isEditable() const
{
    return false;
}

void TrajectoryData::assertIsEditable() const
{
    if (not this->isEditable())
        throw SireError::incompatible_error(QObject::tr(
            "You cannot edit a trajectory of type %1.")
                .arg(this->what()), CODELOC);
}

TrajectoryDataPtr TrajectoryData::makeEditable() const
{
    return TrajectoryDataPtr(new MolTrajectoryData(this->getFrames()));
}

TrajectoryDataPtr TrajectoryData::makeSubsetEditable(int start_atom, int natoms) const
{
    return TrajectoryDataPtr(new MolTrajectoryData(
                                    this->getFrames(start_atom, natoms)));
}

void TrajectoryData::setFrame(int, const Frame&)
{
    this->assertIsEditable();
}

void TrajectoryData::appendFrame(const Frame&)
{
    this->assertIsEditable();
}

void TrajectoryData::insertFrame(int, const Frame&)
{
    this->assertIsEditable();
}

void TrajectoryData::deleteFrame(int)
{
    this->assertIsEditable();
}

/////////////
///////////// Implementation of Trajectory
/////////////

static const RegisterMetaType<Trajectory> r_traj;

SIREMOL_EXPORT QDataStream& operator<<(QDataStream &ds, const Trajectory &traj)
{
    writeHeader(ds, r_traj, 1);

    SharedDataStream sds(ds);
    sds << traj.d << traj.start_atom << traj.natoms;

    return ds;
}

SIREMOL_EXPORT QDataStream& operator>>(QDataStream &ds, Trajectory &traj)
{
    VersionID v = readHeader(ds, r_traj);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> traj.d >> traj.start_atom >> traj.natoms;
    }
    else
        throw version_error(v, "1", r_traj, CODELOC);

    return ds;
}

Trajectory::Trajectory()
           : ConcreteProperty<Trajectory,MoleculeProperty>(),
             start_atom(0), natoms(0)
{}

Trajectory::Trajectory(const TrajectoryData &data)
           : ConcreteProperty<Trajectory,MoleculeProperty>(),
             start_atom(0), natoms(0)
{
    start_atom = 0;
    natoms = data.nAtoms();

    if (data.nFrames() > 0)
        d.append(TrajectoryDataPtr(data));
}

Trajectory::Trajectory(const QList<TrajectoryDataPtr> &data)
           : ConcreteProperty<Trajectory,MoleculeProperty>(),
             start_atom(0), natoms(0)
{
    start_atom = 0;
    natoms = 0;

    for (const auto &ptr : data)
    {
        if (ptr.constData() != 0)
        {
            if (ptr->nFrames() > 0)
            {
                if (natoms == 0)
                {
                    natoms = ptr->nAtoms();
                }
                else if (natoms != ptr->nAtoms())
                {
                    throw SireError::incompatible_error(
                        QObject::tr("Different number of atoms in different "
                                    "trajectories: %1 versus %2.")
                            .arg(natoms).arg(ptr->nAtoms()), CODELOC);
                }

                d.append(ptr);
            }
        }
    }
}

Trajectory::Trajectory(const TrajectoryData &data, int s, int n)
           : ConcreteProperty<Trajectory,MoleculeProperty>(),
             start_atom(s), natoms(n)
{
    if (natoms <= 0 or data.nFrames() <= 0)
        return;

    if (start_atom < 0 or (start_atom+natoms) >= data.nAtoms())
        throw SireError::incompatible_error(QObject::tr(
            "Cannot use start_atom %1 and natoms %2 for a trajectory with %3 atoms.")
                .arg(start_atom).arg(natoms).arg(data.nAtoms()),
                    CODELOC);

    d.append(TrajectoryDataPtr(data));
}

Trajectory::Trajectory(const QList<TrajectoryDataPtr> &data, int s, int n)
           : ConcreteProperty<Trajectory,MoleculeProperty>(),
             start_atom(s), natoms(n)
{
    if (natoms <= 0 or data.count() <= 0)
        return;

    n = 0;

    for (const auto &ptr : data)
    {
        if (ptr.constData() != 0)
        {
            if (ptr->nFrames() > 0)
            {
                if (n == 0)
                {
                    n = ptr->nAtoms();
                }
                else if (n != ptr->nAtoms())
                {
                    throw SireError::incompatible_error(
                        QObject::tr("Different number of atoms in different "
                                    "trajectories: %1 versus %2.")
                            .arg(natoms).arg(ptr->nAtoms()), CODELOC);
                }

                if (start_atom + natoms >= ptr->nAtoms())
                {
                    throw SireError::incompatible_error(
                        QObject::tr("Cannot use start_atom %1 and natoms %2 "
                                    "for a trajectory with %3 atoms.")
                            .arg(start_atom).arg(natoms)
                            .arg(ptr->nAtoms()), CODELOC);
                }

                d.append(ptr);
            }
        }
    }
}

Trajectory::Trajectory(const Trajectory &other)
           : ConcreteProperty<Trajectory,MoleculeProperty>(other),
             d(other.d), start_atom(other.start_atom), natoms(other.natoms)
{}

Trajectory::~Trajectory()
{}

Trajectory& Trajectory::operator=(const Trajectory &other)
{
    if (this != &other)
    {
        MoleculeProperty::operator=(other);
        d = other.d;
        start_atom = other.start_atom;
        natoms = other.natoms;
    }

    return *this;
}

bool Trajectory::operator==(const Trajectory &other) const
{
    return start_atom == other.start_atom and
           natoms == other.natoms and
           d == other.d;
}

bool Trajectory::operator!=(const Trajectory &other) const
{
    return not Trajectory::operator==(other);
}

const char* Trajectory::typeName()
{
    return QMetaType::typeName(qMetaTypeId<Trajectory>());
}

const char* Trajectory::what() const
{
    return Trajectory::typeName();
}

Trajectory* Trajectory::clone() const
{
    return new Trajectory(*this);
}

bool Trajectory::isEmpty() const
{
    return this->nFrames() == 0;
}

int Trajectory::nFrames() const
{
    int nframes = 0;

    for (const auto &ptr : d)
    {
        nframes += ptr->nFrames();
    }

    return nframes;
}

int Trajectory::nAtoms() const
{
    return natoms;
}

// THIS WILL CHANGE 'frame' SO IT HAS THE RIGHT VALUE
// FOR THE SPECIFIED TRAJECTORY
int Trajectory::_getIndexForFrame(int &frame) const
{
    frame = SireID::Index(frame).map(this->nFrames());

    for (int i=0; i<d.count(); ++i)
    {
        int n = d[i]->nFrames();

        if (frame < n)
            return i;
        else
            frame -= n;
    }

    // we should never get here
    throw SireError::program_bug(QObject::tr(
            "How did we get here? %1 %2").arg(frame).arg(this->nFrames()),
                CODELOC);

    return 0;
}

Frame Trajectory::getFrame(int i) const
{
    i = Index(i).map(this->nFrames());

    int idx = _getIndexForFrame(i);

    auto frame = d[idx]->getFrame(i);

    if (frame.nAtoms() == natoms)
        return frame;
    else
        return frame.subset(start_atom, natoms);
}

TrajectoryData& Trajectory::_makeEditable(int &frame)
{
    frame = Index(frame).map(this->nFrames());

    int idx = _getIndexForFrame(frame);

    if (d[idx]->nAtoms() != natoms or (not d[idx]->isEditable()))
    {
        if (d[idx]->nAtoms() == natoms)
        {
            d[idx] = d[idx]->makeEditable();
        }
        else
        {
            d[idx] = d[idx]->makeSubsetEditable(start_atom, natoms);

            if (d[idx]->nAtoms() != natoms)
                throw SireError::program_bug(QObject::tr(
                    "The result from making the subset editable should "
                    "be a data with natoms = %1 (not %2)")
                        .arg(natoms).arg(d[idx]->nAtoms()), CODELOC);
        }
    }

    return *(d[idx]);
}

Frame Trajectory::_subset(const Frame &frame) const
{
    if (frame.nAtoms() == natoms)
        return frame;
    else
        return frame.subset(start_atom, natoms);
}

void Trajectory::setFrame(int i, const Frame &frame)
{
    if (frame.isEmpty())
    {
        this->deleteFrame(i);
    }
    else
    {
        this->_makeEditable(i).setFrame(i, this->_subset(frame));
    }
}

void Trajectory::appendFrame(const Frame &frame)
{
    if (frame.isEmpty())
        return;

    auto subset = this->_subset(frame);

    if (d.isEmpty())
    {
        d.append(TrajectoryDataPtr(new MolTrajectoryData(subset)));
    }
    else if (d.last()->isEditable() and d.last()->nAtoms() == natoms)
    {
        d.last()->appendFrame(frame);
    }
    else
    {
        d.append(TrajectoryDataPtr(new MolTrajectoryData(subset)));
    }
}

void Trajectory::insertFrame(int i, const Frame &frame)
{
    if (frame.isEmpty())
        return;

    auto subset = this->_subset(frame);

    this->_makeEditable(i).insertFrame(i, subset);
}

void Trajectory::deleteFrame(int i)
{
    this->_makeEditable(i).deleteFrame(i);
}

bool Trajectory::isCompatibleWith(const MoleculeInfoData &molinfo) const
{
    return this->nAtoms() == molinfo.nAtoms();
}

///////
/////// Implementation of Frame
///////

static const RegisterMetaType<Frame> r_frame;

SIREMOL_EXPORT QDataStream& operator<<(QDataStream &ds, const Frame &frame)
{
    writeHeader(ds, r_frame, 1);

    SharedDataStream sds(ds);

    sds << frame.coords << frame.vels << frame.frcs
        << frame.spc << frame.t.to(picosecond);

    return ds;
}

SIREMOL_EXPORT QDataStream& operator>>(QDataStream &ds, Frame &frame)
{
    VersionID v = readHeader(ds, r_frame);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> frame.coords >> frame.vels >> frame.frcs >> frame.spc;

        double time;
        sds >> time;

        frame.t = time * picosecond;
    }
    else
        throw version_error(v, "1", r_frame, CODELOC);

    return ds;
}

Frame::Frame()
      : ConcreteProperty<Frame,MoleculeProperty>(),
        spc(Space::null()), t(0)
{}

Frame::Frame(const Molecule &mol, const PropertyMap &map)
      : ConcreteProperty<Frame,MoleculeProperty>(),
        spc(Space::null()), t(0)
{
    this->operator=(Frame(mol.data(), map));
}

Frame::Frame(const MoleculeData &mol, const PropertyMap &map)
      : ConcreteProperty<Frame,MoleculeProperty>(),
        spc(Space::null()), t(0)
{
    if (mol.hasProperty(map["coordinates"]))
    {
        try
        {
            coords = mol.property(map["coordinates"]).asA<AtomCoords>().toVector();
        }
        catch(...)
        {}
    }

    if (mol.hasProperty(map["velocities"]))
    {
        try
        {
            vels = mol.property(map["velocities"]).asA<AtomVelocities>().toVector();
        }
        catch(...)
        {}
    }

    if (mol.hasProperty(map["forces"]))
    {
        try
        {
            frcs = mol.property(map["forces"]).asA<AtomForces>().toVector();
        }
        catch(...)
        {}
    }

    if (mol.hasProperty(map["time"]))
    {
        try
        {
            t = mol.property(map["time"]).asA<GeneralUnitProperty>();
        }
        catch(...)
        {}
    }

    if (mol.hasProperty(map["space"]))
    {
        try
        {
            spc = mol.property(map["space"]).asA<Space>();
        }
        catch(...)
        {}
    }
}

Frame::Frame(const QVector<Vector> &coordinates,
             const QVector<Velocity3D> &velocities,
             const QVector<Force3D> &forces,
             const Space &space,
             SireUnits::Dimension::Time time)
      : ConcreteProperty<Frame,MoleculeProperty>(),
        coords(coordinates), vels(velocities), frcs(forces),
        spc(space), t(time)
{
    int natoms = coords.count();
    bool error = false;

    if (natoms == 0)
    {
        natoms = vels.count();
    }
    else if (vels.count() > 0 and natoms != vels.count())
    {
        error = true;
    }

    if (natoms == 0)
    {
        natoms = frcs.count();
    }
    else if (frcs.count() > 0 and natoms != frcs.count())
    {
        error = true;
    }

    if (error)
    {
        throw SireError::incompatible_error(QObject::tr(
            "Disagreement over the number of atoms! %1 vs %2 vs %3")
                .arg(coords.count()).arg(vels.count()).arg(frcs.count()),
                    CODELOC);
    }
}

Frame::Frame(const Frame &other)
      : ConcreteProperty<Frame,MoleculeProperty>(),
        coords(other.coords),
        vels(other.vels),
        frcs(other.frcs),
        spc(other.spc),
        t(other.t)
{}

Frame::~Frame()
{}

Frame& Frame::operator=(const Frame &other)
{
    if (this != &other)
    {
        coords = other.coords;
        vels = other.vels;
        frcs = other.frcs;
        spc = other.spc;
        t = other.t;
    }

    return *this;
}

bool Frame::operator==(const Frame &other) const
{
    return (this == &other) or
           (coords == other.coords and
            vels == other.vels and
            frcs == other.frcs and
            spc == other.spc and
            t == other.t);
}

bool Frame::operator!=(const Frame &other) const
{
    return not operator==(other);
}

const char* Frame::typeName()
{
    return QMetaType::typeName(qMetaTypeId<Frame>());
}

const char* Frame::what() const
{
    return Frame::typeName();
}

Frame* Frame::clone() const
{
    return new Frame(*this);
}

QString Frame::toString() const
{
    QStringList c;

    if (this->hasCoordinates())
        c.append("C");

    if (this->hasVelocities())
        c.append("V");

    if (this->hasForces())
        c.append("F");

    if (c.isEmpty())
        return QObject::tr("Frame::null");

    return QObject::tr("Frame( nAtoms=%1, time=%2, {%3} )")
                .arg(this->nAtoms()).arg(this->time().toString()).arg(c.join(","));
}

bool Frame::isEmpty() const
{
    return coords.isEmpty() and vels.isEmpty() and frcs.isEmpty();
}

bool Frame::hasCoordinates() const
{
    return not coords.isEmpty();
}

bool Frame::hasVelocities() const
{
    return not vels.isEmpty();
}

bool Frame::hasForces() const
{
    return not frcs.isEmpty();
}

QVector<Vector> Frame::coordinates() const
{
    return coords;
}

QVector<Velocity3D> Frame::velocities() const
{
    return vels;
}

QVector<Force3D> Frame::forces() const
{
    return frcs;
}

SireUnits::Dimension::Time Frame::time() const
{
    return t;
}

const Space& Frame::space() const
{
    return *spc;
}

int Frame::nAtoms() const
{
    int nats = coords.count();

    if (nats == 0)
    {
        nats = vels.count();
    }

    if (nats == 0)
    {
        nats = frcs.count();
    }

    return nats;
}

template<class T>
QVector<T> _mid(const QVector<T> &v, int start, int count)
{
    if (v.count() == 0)
        return v;

    if (start < 0 or count < 0 or start+count >= v.count())
        throw SireError::incompatible_error(QObject::tr(
            "Cannot subset a vector of length %1 using start %2, count %3.")
                .arg(v.count()).arg(start).arg(count), CODELOC);

    return v.mid(start, count);
}

Frame Frame::subset(int start_atom, int natoms) const
{
    Frame ret(*this);

    ret.coords = _mid(coords, start_atom, natoms);
    ret.vels = _mid(vels, start_atom, natoms);
    ret.frcs = _mid(frcs, start_atom, natoms);

    return ret;
}

bool Frame::isCompatibleWith(const MoleculeInfoData &molinfo) const
{
    return this->nAtoms() == molinfo.nAtoms();
}
