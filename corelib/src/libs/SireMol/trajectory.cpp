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

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMol;
using namespace SireVol;
using namespace SireBase;
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

void TrajectoryData::assertContainsSpace() const
{
    if (not this->containsSpace())
        throw SireError::incompatible_error(QObject::tr(
            "This trajectory does not contain any space information!"
        ), CODELOC);
}

void TrajectoryData::assertContainsCoordinates() const
{
    if (not this->containsCoordinates())
        throw SireError::incompatible_error(QObject::tr(
            "This trajectory does not contain any coordinates!"
        ), CODELOC);
}

void TrajectoryData::assertContainsVelocities() const
{
    if (not this->containsVelocities())
        throw SireError::incompatible_error(QObject::tr(
            "This trajectory does not contain any velocities!"
        ), CODELOC);
}

void TrajectoryData::assertContainsForces() const
{
    if (not this->containsForces())
        throw SireError::incompatible_error(QObject::tr(
            "This trajectory does not contain any forces!"
        ), CODELOC);
}

void TrajectoryData::assertSupportsSpace() const
{
    if (not this->supportsCoordinates())
        throw SireError::incompatible_error(QObject::tr(
            "This trajectory does not support space information!"
        ), CODELOC);
}

void TrajectoryData::assertSupportsCoordinates() const
{
    if (not this->supportsCoordinates())
        throw SireError::incompatible_error(QObject::tr(
            "This trajectory does not support coordinates!"
        ), CODELOC);
}

void TrajectoryData::assertSupportsVelocities() const
{
    if (not this->supportsVelocities())
        throw SireError::incompatible_error(QObject::tr(
            "This trajectory does not support velocities!"
        ), CODELOC);
}

void TrajectoryData::assertSupportsForces() const
{
    if (not this->supportsForces())
        throw SireError::incompatible_error(QObject::tr(
            "This trajectory does not support forces!"
        ), CODELOC);
}

int TrajectoryData::nFiles() const
{
    return this->filenames().count();
}

bool TrajectoryData::isEmpty() const
{
    return this->nFrames() == 0;
}

void TrajectoryData::appendFrame(const AtomCoords &coords, qint64 idx)
{
    this->insertFrame(this->nFrames(), coords, idx);
}

void TrajectoryData::appendFrame(const AtomVelocities &velocities, qint64 idx)
{
    this->insertFrame(this->nFrames(), velocities, idx);
}

void TrajectoryData::appendFrame(const AtomForces &forces, qint64 idx)
{
    this->insertFrame(this->nFrames(), forces, idx);
}

void TrajectoryData::insertFrame(int i, const AtomCoords &coords, qint64 idx)
{
    this->assertSupportsCoordinates();
    this->insertFrame(i);
    this->setFrame(i, coords, idx);
}

void TrajectoryData::insertFrame(int i, const AtomVelocities &velocities, qint64 idx)
{
    this->assertSupportsVelocities();
    this->insertFrame(i);
    this->setFrame(i, velocities, idx);
}

void TrajectoryData::insertFrame(int i, const AtomForces &forces, qint64 idx)
{
    this->assertSupportsForces();
    this->insertFrame(i);
    this->setFrame(i, forces, idx);
}

/////////////
///////////// Implementation of Trajectory
/////////////

static const RegisterMetaType<Trajectory> r_traj;

SIREMOL_EXPORT QDataStream& operator<<(QDataStream &ds, const Trajectory &traj)
{
    writeHeader(ds, r_traj, 1);

    SharedDataStream sds(ds);
    sds << traj.d << traj.idx;

    return ds;
}

SIREMOL_EXPORT QDataStream& operator>>(QDataStream &ds, Trajectory &traj)
{
    VersionID v = readHeader(ds, r_traj);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> traj.d >> traj.idx;
    }
    else
        throw version_error(v, "1", r_traj, CODELOC);

    return ds;
}

Trajectory::Trajectory()
           : ConcreteProperty<Trajectory,MoleculeProperty>(),
             idx(0)
{}

Trajectory::Trajectory(const TrajectoryData &data, qint64 index)
           : ConcreteProperty<Trajectory,MoleculeProperty>(),
             d(data), idx(index)
{
    data.assertValidIndex(index);
}

Trajectory::Trajectory(const Trajectory &other)
           : ConcreteProperty<Trajectory,MoleculeProperty>(other),
             d(other.d), idx(other.idx)
{}

Trajectory::~Trajectory()
{}

Trajectory& Trajectory::operator=(const Trajectory &other)
{
    if (this != &other)
    {
        MoleculeProperty::operator=(other);
        d = other.d;
        idx = other.idx;
    }

    return *this;
}

bool Trajectory::operator==(const Trajectory &other) const
{
    if (d.data() == 0)
        return other.d.data() == 0 and idx == other.idx;
    else
        return idx == other.idx and *d == *(other.d);
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
    if (d.data() == 0)
        return 1;
    else
        return d->nFrames();
}

SireUnits::Dimension::Time Trajectory::timestep() const
{
    if (d.data() == 0)
        return SireUnits::Dimension::Time(0);
    else
        return d->timestep();
}

const Space& Trajectory::getSpace(int i) const
{
    this->assertContainsSpace();
    return d->getSpace(i);
}

AtomCoords Trajectory::getFrame(int i, const AtomCoords &coords) const
{
    i = Index(i).map(this->nFrames());

    if (d.data() == 0)
        return coords;
    else
    {
        d->assertValidIndex(idx + coords.nAtoms());
        return d->getFrame(i, coords, idx);
    }
}

AtomVelocities Trajectory::getFrame(int i, const AtomVelocities &velocities) const
{
    i = Index(i).map(this->nFrames());

    if (d.data() == 0)
        return velocities;
    else
    {
        d->assertValidIndex(idx + velocities.nAtoms());
        return d->getFrame(i, velocities, idx);
    }
}

AtomForces Trajectory::getFrame(int i, const AtomForces &forces) const
{
    i = Index(i).map(this->nFrames());

    if (d.data() == 0)
        return forces;
    else
    {
        d->assertValidIndex(idx + forces.nAtoms());
        return d->getFrame(i, forces, idx);
    }
}

void Trajectory::setFrame(int i, const AtomCoords &coords)
{
    if (d.constData() == 0)
        throw SireError::incompatible_error(
            QObject::tr("You cannot edit a null trajectory."), CODELOC);

    this->assertSupportsCoordinates();
    d.constData()->assertValidIndex(idx + coords.nAtoms());

    d->setFrame(i, coords, idx);
}

void Trajectory::setFrame(int i, const AtomVelocities &velocities)
{
    if (d.constData() == 0)
        throw SireError::incompatible_error(
            QObject::tr("You cannot edit a null trajectory."), CODELOC);

    this->assertSupportsVelocities();
    d.constData()->assertValidIndex(idx + velocities.nAtoms());

    d->setFrame(i, velocities, idx);
}

void Trajectory::setFrame(int i, const AtomForces &forces)
{
    if (d.constData() == 0)
        throw SireError::incompatible_error(
            QObject::tr("You cannot edit a null trajectory."), CODELOC);

    this->assertSupportsForces();
    d.constData()->assertValidIndex(idx + forces.nAtoms());

    d->setFrame(i, forces, idx);
}

void Trajectory::appendFrame(const AtomCoords &coords)
{
    if (d.constData() == 0)
        throw SireError::incompatible_error(
            QObject::tr("You cannot edit a null trajectory."), CODELOC);

    this->assertSupportsCoordinates();
    d.constData()->assertValidIndex(idx + coords.nAtoms());

    d->appendFrame(coords, idx);
}

void Trajectory::appendFrame(const AtomVelocities &velocities)
{
    if (d.constData() == 0)
        throw SireError::incompatible_error(
            QObject::tr("You cannot edit a null trajectory."), CODELOC);

    this->assertSupportsVelocities();
    d.constData()->assertValidIndex(idx + velocities.nAtoms());

    d->appendFrame(velocities, idx);
}

void Trajectory::appendFrame(const AtomForces &forces)
{
    if (d.constData() == 0)
        throw SireError::incompatible_error(
            QObject::tr("You cannot edit a null trajectory."), CODELOC);

    this->assertSupportsForces();
    d.constData()->assertValidIndex(idx + forces.nAtoms());

    d->appendFrame(forces, idx);
}

void Trajectory::insertFrame(int i, const AtomCoords &coords)
{
    if (d.constData() == 0)
        throw SireError::incompatible_error(
            QObject::tr("You cannot edit a null trajectory."), CODELOC);

    this->assertSupportsCoordinates();
    d.constData()->assertValidIndex(idx + coords.nAtoms());

    d->insertFrame(i, coords, idx);
}

void Trajectory::insertFrame(int i, const AtomVelocities &velocities)
{
    if (d.constData() == 0)
        throw SireError::incompatible_error(
            QObject::tr("You cannot edit a null trajectory."), CODELOC);

    this->assertSupportsVelocities();
    d.constData()->assertValidIndex(idx + velocities.nAtoms());

    d->insertFrame(i, velocities, idx);
}

void Trajectory::insertFrame(int i, const AtomForces &forces)
{
    if (d.constData() == 0)
        throw SireError::incompatible_error(
            QObject::tr("You cannot edit a null trajectory."), CODELOC);

    this->assertSupportsForces();
    d.constData()->assertValidIndex(idx + forces.nAtoms());

    d->insertFrame(i, forces, idx);
}

void Trajectory::deleteFrame(int i)
{
    if (d.constData() == 0)
        throw SireError::incompatible_error(
            QObject::tr("You cannot edit a null trajectory."), CODELOC);

    i = Index(i).map(this->nFrames());
    d->deleteFrame(i);
}

bool Trajectory::isCompatibleWith(const MoleculeInfoData &molinfo) const
{
    if (d.data() == 0)
        return true;

    try
    {
        d->assertValidIndex(idx + molinfo.nAtoms());
        return true;
    }
    catch(...)
    {
        return false;
    }
}

bool Trajectory::containsSpace() const
{
    if (d.data() == 0)
        return false;
    else
        return d->containsCoordinates();
}

bool Trajectory::containsCoordinates() const
{
    if (d.data() == 0)
        return true;
    else
        return d->containsCoordinates();
}

bool Trajectory::containsVelocities() const
{
    if (d.data() == 0)
        return true;
    else
        return d->containsVelocities();
}

bool Trajectory::containsForces() const
{
    if (d.data() == 0)
        return true;
    else
        return d->containsForces();
}

bool Trajectory::supportsSpace() const
{
    if (d.data() == 0)
        return false;
    else
        return d->supportsCoordinates();
}

bool Trajectory::supportsCoordinates() const
{
    if (d.data() == 0)
        return true;
    else
        return d->supportsCoordinates();
}

bool Trajectory::supportsVelocities() const
{
    if (d.data() == 0)
        return true;
    else
        return d->supportsVelocities();
}

bool Trajectory::supportsForces() const
{
    if (d.data() == 0)
        return true;
    else
        return d->supportsForces();
}

void Trajectory::assertContainsSpace() const
{
    if (d.data() != 0)
        d->assertContainsSpace();
    else
        throw SireError::incompatible_error( QObject::tr(
            "The null trajectory does not contain space information."),
                CODELOC);
}

void Trajectory::assertContainsCoordinates() const
{
    if (d.data() != 0)
        d->assertContainsCoordinates();
}

void Trajectory::assertContainsVelocities() const
{
    if (d.data() != 0)
        d->assertContainsVelocities();
}

void Trajectory::assertContainsForces() const
{
    if (d.data() != 0)
        d->assertContainsForces();
}

void Trajectory::assertSupportsSpace() const
{
    if (d.data() != 0)
        d->assertContainsSpace();
    else
        throw SireError::incompatible_error( QObject::tr(
            "The null trajectory does not support space information."),
                CODELOC);
}

void Trajectory::assertSupportsCoordinates() const
{
    if (d.data() != 0)
        d->assertSupportsCoordinates();
}

void Trajectory::assertSupportsVelocities() const
{
    if (d.data() != 0)
        d->assertSupportsVelocities();
}

void Trajectory::assertSupportsForces() const
{
    if (d.data() != 0)
        d->assertSupportsForces();
}
