/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2007  Christopher Woods
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

#include "moleculeview.h"
#include "atomselection.h"
#include "atom.h"
#include "cutgroup.h"
#include "residue.h"
#include "chain.h"
#include "segment.h"
#include "selector.hpp"
#include "molecule.h"
#include "select.h"
#include "trajectory.h"

#include "SireVol/space.h"

#include "SireBase/slice.h"

#include "SireBase/generalunitproperty.h"

#include "SireVol/space.h"

#include "SireBase/errors.h"
#include "SireError/errors.h"
#include "SireMol/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>

using namespace SireStream;
using namespace SireBase;
using namespace SireVol;
using namespace SireMol;

RegisterMetaType<MoleculeView> r_molview( MAGIC_ONLY,
                                          "SireMol::MoleculeView" );

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                       const MoleculeView &molview)
{
    writeHeader(ds, r_molview, 2);

    SharedDataStream sds(ds);

    sds << molview.d
        << static_cast<const Property&>(molview);

    return ds;
}

/** Deserialise from a binary datastream */
QDataStream &operator>>(QDataStream &ds,
                                       MoleculeView &molview)
{
    VersionID v = readHeader(ds, r_molview);

    if (v == 2)
    {
        SharedDataStream sds(ds);

        sds >> molview.d >> static_cast<Property&>(molview);
    }
    else if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> molview.d;
    }
    else
        throw version_error(v, "1", r_molview, CODELOC);

    return ds;
}

/** Null constructor */
MoleculeView::MoleculeView() : Property(), d( MoleculeData::null() )
{}

/** Construct a view of the molecule whose data is in 'moldata' */
MoleculeView::MoleculeView(const MoleculeData &moldata)
             : Property(), d(moldata)
{}

/** Copy constructor */
MoleculeView::MoleculeView(const MoleculeView &other)
             : Property(other), d(other.d)
{}

/** Destructor */
MoleculeView::~MoleculeView()
{}

/** Copy assignment operator */
MoleculeView& MoleculeView::operator=(const MoleculeView &other)
{
    Property::operator=(other);
    d = other.d;
    return *this;
}

/** Comparison operator */
bool MoleculeView::operator==(const MoleculeView &other) const
{
    return d == other.d or *d == *(other.d);
}

/** Comparison operator */
bool MoleculeView::operator!=(const MoleculeView &other) const
{
    return d != other.d and *d != *(other.d);
}

/** Return whether or not this molecule view is null */
bool MoleculeView::isNull() const
{
    return *d == *MoleculeData::null();
}

int MoleculeView::nFrames() const
{
    return this->nFrames(PropertyMap());
}

int MoleculeView::nFrames(const SireBase::PropertyMap &map) const
{
    const auto traj_prop = map["trajectory"];

    if (d->hasProperty(traj_prop))
    {
        return d->property(traj_prop).asA<Trajectory>().nFrames();
    }
    else
    {
        return 1;
    }
}

void MoleculeView::_fromFrame(const Frame &frame,
                              const SireBase::PropertyMap &map)
{
    const auto coords_prop = map["coordinates"];
    const auto vels_prop = map["velocities"];
    const auto frcs_prop = map["forces"];
    const auto space_prop = map["space"];
    const auto time_prop = map["time"];

    if (frame.hasCoordinates() and coords_prop.hasSource())
    {
        auto coords = AtomCoords(d->info());
        coords.copyFrom(frame.coordinates());
        d->setProperty(coords_prop.source(), coords);
    }

    if (frame.hasVelocities() and vels_prop.hasSource())
    {
        auto vels = AtomVelocities(d->info());
        vels.copyFrom(frame.velocities());
        d->setProperty(vels_prop.source(), vels);
    }

    if (frame.hasForces() and frcs_prop.hasSource())
    {
        auto frcs = AtomForces(d->info());
        frcs.copyFrom(frame.forces());
        d->setProperty(frcs_prop.source(), frcs);
    }

    if (space_prop.hasSource())
    {
        d->setProperty(space_prop.source(), frame.space());
    }

    if (time_prop.hasSource())
    {
        d->setProperty(time_prop.source(), GeneralUnitProperty(frame.time()));
    }
}

Frame MoleculeView::_toFrame(const SireBase::PropertyMap &map) const
{
    const auto coords_prop = map["coordinates"];
    const auto vels_prop = map["velocities"];
    const auto frcs_prop = map["forces"];
    const auto space_prop = map["space"];
    const auto time_prop = map["time"];

    QVector<Vector> coords;
    QVector<Velocity3D> vels;
    QVector<Force3D> frcs;
    SpacePtr space;
    SireUnits::Dimension::Time time(0);

    if (d->hasProperty(coords_prop))
    {
        coords = d->property(coords_prop).asA<AtomCoords>().toVector();
    }

    if (d->hasProperty(vels_prop))
    {
        vels = d->property(vels_prop).asA<AtomVelocities>().toVector();
    }

    if (d->hasProperty(frcs_prop))
    {
        frcs = d->property(frcs_prop).asA<AtomForces>().toVector();
    }

    if (d->hasProperty(space_prop))
    {
        space = d->property(space_prop).asA<Space>();
    }

    if (d->hasProperty(time_prop))
    {
        time = d->property(time_prop).asA<GeneralUnitProperty>();
    }

    return Frame(coords, vels, frcs, space, time);
}

void MoleculeView::loadFrame(int frame)
{
    this->loadFrame(frame, PropertyMap());
}

void MoleculeView::saveFrame(int frame)
{
    this->saveFrame(frame, PropertyMap());
}

void MoleculeView::saveFrame()
{
    this->saveFrame(PropertyMap());
}

void MoleculeView::deleteFrame(int frame)
{
    this->deleteFrame(frame, PropertyMap());
}

void MoleculeView::loadFrame(int frame, const SireBase::PropertyMap &map)
{
    const auto traj_prop = map["trajectory"];

    if (frame == 0 and (not d->hasProperty(traj_prop)))
        return;

    auto traj = d->property(traj_prop).asA<Trajectory>();

    this->_fromFrame(traj[frame], map);
}

void MoleculeView::saveFrame(int frame, const SireBase::PropertyMap &map)
{
    const auto traj_prop = map["trajectory"];

    if (not (traj_prop.hasSource()))
        return;

    Trajectory traj;

    if (d->hasProperty(traj_prop))
    {
        traj = d->property(traj_prop.source()).asA<Trajectory>();
    }

    if (frame == traj.nFrames())
        this->saveFrame();
    else
    {
        traj.setFrame(frame, this->_toFrame(map));
        d->setProperty(traj_prop.source(), traj);
    }
}

void MoleculeView::saveFrame(const SireBase::PropertyMap &map)
{
    const auto traj_prop = map["trajectory"];

    if (not (traj_prop.hasSource()))
        return;

    Trajectory traj;

    if (d->hasProperty(traj_prop))
    {
        traj = d->property(traj_prop.source()).asA<Trajectory>();
    }

    traj.appendFrame(this->_toFrame(map));
    d->setProperty(traj_prop.source(), traj);
}

void MoleculeView::deleteFrame(int frame, const SireBase::PropertyMap &map)
{
    const auto traj_prop = map["trajectory"];

    if (not (traj_prop.hasSource() and d->hasProperty(traj_prop)))
        return;

    auto traj = d->property(traj_prop.source()).asA<Trajectory>();

    traj.deleteFrame(frame);

    d->setProperty(traj_prop.source(), traj);
}

/** Return whether or not this view is of the same molecule as 'other'
    (albeit perhaps a different version of the molecule) */
bool MoleculeView::isSameMolecule(const MoleculeData &other) const
{
    return d->number() == other.number();
}

/** Return whether or not this view is of the same molecule as 'other'
    (albeit perhaps a different version of the molecule) */
bool MoleculeView::isSameMolecule(const MoleculeView &other) const
{
    return this->isSameMolecule(other.data());
}

/** Assert that this view is looking at the molecule whose data is
    in 'other' (albeit perhaps a different version of that molecule)

    \throw SireError::incompatible_error
*/
void MoleculeView::assertSameMolecule(const MoleculeData &other) const
{
    if (d->number() != other.number())
        //these are different molecules!
        throw SireError::incompatible_error( QObject::tr(
            "The molecules \"%1\", number %2, and \"%3\", number %4, "
            "are different, and therefore incompatible.")
                .arg(d->name()).arg(d->number())
                .arg(other.name()).arg(other.number()),
                    CODELOC );
}

/** Assert that this is a view of the same molecule as 'other'
    (albeit at a different version)

    \throw SireError::incompatible_error
*/
void MoleculeView::assertSameMolecule(const MoleculeView &other) const
{
    this->assertSameMolecule(other.data());
}

/** Update this view with a new version of the molecule. You
    can only update the molecule if it has the same layout UID
    (so same atoms, residues, cutgroups etc.)

    \throw SireError::incompatible_error
*/
void MoleculeView::update(const MoleculeData &moldata)
{
    this->assertSameMolecule(moldata);
    d->info().assertEqualTo(moldata.info());

    d = moldata;
}

/** Synonym for MoleculeView::propertyKeys */
QStringList MoleculeView::keys() const
{
    return this->propertyKeys();
}

/** Return the type of the property at key 'key'

    \throw SireBase::missing_property
*/
const char* MoleculeView::propertyType(const PropertyName &key) const
{
    return d->property(key).what();
}

/** Return the type of the metadata at metakey 'metakey'

    \throw SireBase::missing_property
*/
const char* MoleculeView::metadataType(const PropertyName &metakey) const
{
    return d->metadata(metakey).what();
}

/** Return the type of the metadata at metakey 'metakey'
    for the property at key 'key'

    \throw SireBase::missing_property
*/
const char* MoleculeView::metadataType(const PropertyName &key,
                                       const PropertyName &metakey) const
{
    return d->metadata(key, metakey).what();
}

/** Assert that this view contains the atom at index 'atomidx'

    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
void MoleculeView::assertContains(AtomIdx atomidx) const
{
    if (not this->selection().selected(atomidx))
        throw SireMol::missing_atom( QObject::tr(
            "This view of the molecule \"%1\" (%2) does not "
            "contain the atom at index %3.")
                .arg(d->name()).arg(d->number())
                .arg(atomidx), CODELOC );
}

/** Assert that this contains a property at key 'key'

    \throw SireBase::missing_property
*/
void MoleculeView::assertHasProperty(const PropertyName &key) const
{
    if (not this->hasProperty(key))
        throw SireBase::missing_property( QObject::tr(
            "This view of the molecule \"%1\" (view type %2) "
            "does not have a valid property at key \"%3\".")
                .arg(d->name())
                .arg(this->what())
                .arg(key.toString()), CODELOC );
}

/** Assert that this contains some metadata at metakey 'metakey'

    \throw SireBase::missing_property
*/
void MoleculeView::assertHasMetadata(const PropertyName &metakey) const
{
    if (not this->hasMetadata(metakey))
        throw SireBase::missing_property( QObject::tr(
            "This view of the molecule \"%1\" (view type %2) "
            "does not have some valid metadata at metakey \"%3\".")
                .arg(d->name())
                .arg(this->what())
                .arg(metakey.toString()), CODELOC );
}

/** Assert that this contains some metadata at metakey 'metakey'
    for the property at key 'key'

    \throw SireBase::missing_property
*/
void MoleculeView::assertHasMetadata(const PropertyName &key,
                                     const PropertyName &metakey) const
{
    if (not this->hasMetadata(key,metakey))
        throw SireBase::missing_property( QObject::tr(
            "This view of the molecule \"%1\" (view type %2) "
            "does not have some valid metadata at metakey \"%3\" "
            "for the property at key \"%4\".")
                .arg(d->name())
                .arg(this->what())
                .arg(metakey.toString())
                .arg(key.toString()), CODELOC );
}

Atom MoleculeView::atom(int i, const PropertyMap &map) const
{
    auto s = this->selection();

    if (s.selectedAllAtoms())
    {
        return this->atom(AtomIdx(i), map);
    }
    else
    {
        auto atomidxs = s.selectedAtoms();
        return this->atom(atomidxs.at(Index(i).map(atomidxs.count())), map);
    }
}

Atom MoleculeView::atom(const QString &name, const PropertyMap &map) const
{
    try
    {
        return this->atom(AtomID::fromString(name), map);
    }
    catch(const SireMol::duplicate_atom &e)
    {
        throw e;
    }
    catch(const SireError::exception &e)
    {
        try
        {
            auto a = this->search(name).views().at(0).atom();
            return this->atom(a.index(), map);
        }
        catch(...)
        {
            if (name.length() < 5)
                //likely a name error
                e.throwSelf();
            else
                //likely a syntax error
                throw;
        }
    }

    return Atom();
}

/** Return the atom in this view that matches the ID 'atomid'

    \throw SireMol::missing_atom
    \throw SireError::invalid_index
    \throw SireMol::duplicate_atom
*/
Atom MoleculeView::atom(const AtomID &atomid, const PropertyMap &map) const
{
    return atomid.selectFrom(*this, map);
}

template<class T>
QList<qint64> _toIndicies(const QList<T> &ids)
{
    QList<qint64> idxs;

    for (const auto &id : ids)
    {
        idxs.append(id.value());
    }

    return idxs;
}

Selector<Atom> MoleculeView::atoms(const QString &name,
                                   const PropertyMap &map) const
{
    try
    {
        return this->atoms(AtomID::fromString(name), map);
    }
    catch(const SireError::exception &e)
    {
        try
        {
            const auto a = this->search(name).views();

            if (a.count() == 0)
                throw SireMol::missing_atom(QObject::tr(
                    "No atom matches '%1'").arg(name), CODELOC);

            return this->molecule().atoms(
                        _toIndicies(a[0].atoms().IDs()), map);
        }
        catch(...)
        {
            if (name.length() < 5)
                //likely a name error
                e.throwSelf();
            else
                //likely a syntax error
                throw;
        }
    }

    return Selector<Atom>();
}

Selector<Atom> MoleculeView::atoms(const QStringList &names,
                                   const PropertyMap &map) const
{
    if (names.count() == 0)
        throw SireMol::missing_atom(QObject::tr(
            "You must specify some names.."), CODELOC);

    auto s = this->atoms(names[0], map);

    for (int i=1; i<names.count(); ++i)
    {
        s = s + this->atoms(names[i], map);
    }

    return s;
}

Selector<Atom> MoleculeView::atoms(const QList<qint64> &values,
                                   const PropertyMap &map) const
{
    if (values.count() == 0)
        throw SireError::invalid_index(QObject::tr(
            "You must specify some indexes.."), CODELOC);

    if (this->isA< Selector<Atom> >())
    {
        auto IDs = this->asA< Selector<Atom> >().indexes();

        QList<AtomIdx> ret;
        ret.reserve(values.count());

        for (const auto &value : values)
        {
            ret.append(IDs.at(SireID::Index(value).map(IDs.count())));
        }

        return Selector<Atom>(this->data(), ret);
    }

    const auto s = this->selection();
    QList<AtomIdx> idxs;

    if (s.selectedAllAtoms())
    {
        for (auto value : values)
        {
            idxs.append(AtomIdx(value));
        }
    }
    else
    {
        auto atomidxs = s.selectedAtoms();

        for (auto value : values)
        {
            idxs.append(atomidxs.at(Index(value).map(atomidxs.count())));
        }
    }

    return Selector<Atom>(this->data(), idxs);
}

Selector<Atom> MoleculeView::atoms(const Slice &slice,
                                   const PropertyMap &map) const
{
    if (this->isA< Selector<Atom> >())
    {
        auto IDs = this->asA< Selector<Atom> >().indexes();

        QList<AtomIdx> ret;

        for (auto it = slice.begin(IDs.count()); not it.atEnd(); it.next())
        {
            ret.append(IDs.at(it.value()));
        }

        return Selector<Atom>(this->data(), ret);
    }

    const auto s = this->selection();
    QList<AtomIdx> idxs;

    if (s.selectedAllAtoms())
    {
        for (auto it = slice.begin(s.nAtoms()); not it.atEnd(); it.next())
        {
            idxs.append(AtomIdx(it.value()));
        }
    }
    else
    {
        auto atomidxs = s.selectedAtoms();

        for (auto it = slice.begin(atomidxs.count()); not it.atEnd(); it.next())
        {
            idxs.append(atomidxs.at(Index(it.value()).map(atomidxs.count())));
        }
    }

    return Selector<Atom>(this->data(), idxs);
}

/** Return the atoms from this view that match the ID 'atomid'

    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
Selector<Atom> MoleculeView::atoms(const AtomID &atomid,
                                   const PropertyMap &map) const
{
    return atomid.selectAllFrom(*this, map);
}

/** Return this view as a Atom - this will only work if
    this view contains only a single atom

    \throw SireMol::duplicate_atom
*/
Atom MoleculeView::atom() const
{
    QVector<AtomIdx> selected_atoms = this->selection().selectedAtoms();

    if (selected_atoms.isEmpty())
        throw SireMol::missing_atom( QObject::tr(
                "This view does not contain any atoms."), CODELOC );

    else if (selected_atoms.count() != 1)
        throw SireMol::duplicate_atom( QObject::tr(
                "Cannot convert this view (%1) into an Atom as "
                "we can only do this is just one atom is selected. "
                "The number of matching atoms is %2.")
                    .arg(this->toString()).arg(selected_atoms.count()),
                        CODELOC );

    return Atom(this->data(), selected_atoms.at(0));
}

/** Return all of the atoms in this view

    \throw SireMol::missing_atom
*/
Selector<Atom> MoleculeView::atoms() const
{
    if (this->isA< Selector<Atom> >())
        return this->asA< Selector<Atom> >();

    AtomSelection selected_atoms = this->selection();

    if (selected_atoms.selectedNone())
        throw SireMol::missing_atom( QObject::tr(
                "No atoms are available in this view (%1).")
                    .arg(this->toString()), CODELOC );

    return Selector<Atom>(this->data(), this->selection());
}

CutGroup MoleculeView::cutGroup(int i, const PropertyMap &map) const
{
    auto s = this->selection();

    if (s.selectedAllCutGroups())
    {
        return this->cutGroup(AtomIdx(i), map);
    }
    else
    {
        auto cgidxs = s.selectedCutGroups();
        return this->cutGroup(cgidxs.at(Index(i).map(cgidxs.count())), map);
    }
}

CutGroup MoleculeView::cutGroup(const QString &name, const PropertyMap &map) const
{
    return this->cutGroup(CGID::fromString(name), map);
}

/** Return the CutGroup whose atoms are in this view that matches
    the ID in 'cgid'

    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
    \throw SireMol::duplicate_cutgroup
*/
CutGroup MoleculeView::cutGroup(const CGID &cgid, const PropertyMap &map) const
{
    return cgid.selectFrom(*this, map);
}

Selector<CutGroup> MoleculeView::cutGroups(const QString &name,
                                           const PropertyMap &map) const
{
    return this->cutGroups(CGName(name), map);
}

Selector<CutGroup> MoleculeView::cutGroups(const Slice &slice,
                                           const PropertyMap &map) const
{
    if (this->isA< Selector<CutGroup> >())
    {
        auto IDs = this->asA< Selector<CutGroup> >().indexes();

        QList<CGIdx> ret;

        for (auto it = slice.begin(IDs.count()); not it.atEnd(); it.next())
        {
            ret.append(IDs.at(it.value()));
        }

        return Selector<CutGroup>(this->data(), ret);
    }

    const auto s = this->selection();
    QList<CGIdx> idxs;

    if (s.selectedAllCutGroups())
    {

        for (auto it = slice.begin(s.nCutGroups()); not it.atEnd(); it.next())
        {
            idxs.append(CGIdx(it.value()));
        }
    }
    else
    {
        auto cgidxs = s.selectedCutGroups();

        for (auto it = slice.begin(cgidxs.count()); not it.atEnd(); it.next())
        {
            idxs.append(cgidxs.at(it.value()));
        }
    }

    return Selector<CutGroup>(*this, idxs);
}

/** Return the CutGroups whose atoms are in this view that match
    the ID in 'cgid'

    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
*/
Selector<CutGroup> MoleculeView::cutGroups(const CGID &cgid,
                                           const PropertyMap &map) const
{
    return cgid.selectAllFrom(*this, map);
}

/** Return the CutGroup that contains the atom(s) in this view

    \throw SireMol::missing_cutgroup
    \throw SireMol::duplicate_cutgroup
*/
CutGroup MoleculeView::cutGroup() const
{
    QList<CGIdx> selected_cgs = this->selection().selectedCutGroups();

    if (selected_cgs.isEmpty())
        throw SireMol::missing_cutgroup( QObject::tr(
                "This view does not contain any CutGroups."), CODELOC );

    else if (selected_cgs.count() != 1)
        throw SireMol::duplicate_cutgroup( QObject::tr(
                "Cannot convert this view (%1) into a CutGroup as "
                "we can only do this is just one CutGroup is selected. "
                "These CutGroups are selected ; %2.")
                    .arg(this->toString(), Sire::toString(selected_cgs)),
                        CODELOC );

    return CutGroup(this->data(), selected_cgs.at(0));
}

/** Return all of the CutGroups that are involved in this view

    \throw SireMol::missing_cutgroup
*/
Selector<CutGroup> MoleculeView::cutGroups() const
{
    if (this->isA< Selector<CutGroup> >())
        return this->asA< Selector<CutGroup> >();

    QList<CGIdx> selected_cgs = this->selection().selectedCutGroups();

    if (selected_cgs.isEmpty())
        throw SireMol::missing_cutgroup( QObject::tr(
                "This view does not contain any CutGroups."), CODELOC );

    return Selector<CutGroup>(this->data(), selected_cgs);
}

Residue MoleculeView::residue(int i, const PropertyMap &map) const
{
    auto s = this->selection();

    if (s.selectedAllResidues())
    {
        return this->residue(ResIdx(i), map);
    }
    else
    {
        auto residxs = s.selectedResidues();
        return this->residue(residxs.at(Index(i).map(residxs.count())), map);
    }
}

Residue MoleculeView::residue(const QString &name, const PropertyMap &map) const
{
    try
    {
        return this->residue(ResID::fromString(name), map);
    }
    catch(const SireMol::duplicate_residue &e)
    {
        throw e;
    }
    catch(const SireError::exception &e)
    {
        try
        {
            auto a = this->search(name).views().at(0).residue();
            return this->residue(a.index(), map);
        }
        catch(...)
        {
            if (name.length() < 5)
                //likely a name error
                e.throwSelf();
            else
                //likely a syntax error
                throw;
        }
    }

    return Residue();
}

/** Return the residue from this view that matches the ID 'resid'

    \throw SireMol::missing_residue
    \throw SireError::invalid_index
    \throw SireMol::duplicate_residue
*/
Residue MoleculeView::residue(const ResID &resid, const PropertyMap &map) const
{
    return resid.selectFrom(*this, map);
}

Selector<Residue> MoleculeView::residues(const QString &name,
                                         const PropertyMap &map) const
{
    try
    {
        return this->residues(ResID::fromString(name), map);
    }
    catch(const SireError::exception &e)
    {
        try
        {
            const auto a = this->search(name).views();

            if (a.count() == 0)
                throw SireMol::missing_residue(QObject::tr(
                    "No residue matches '%1'").arg(name), CODELOC);

            return this->molecule().residues(
                        _toIndicies(a[0].residues().IDs()), map);
        }
        catch(...)
        {
            if (name.length() < 5)
                //likely a name error
                e.throwSelf();
            else
                //likely a syntax error
                throw;
        }
    }

    return Selector<Residue>();
}

Selector<Residue> MoleculeView::residues(const QStringList &names,
                                         const PropertyMap &map) const
{
    if (names.count() == 0)
        throw SireMol::missing_residue(QObject::tr(
            "You must specify some names.."), CODELOC);

    auto s = this->residues(names[0], map);

    for (int i=1; i<names.count(); ++i)
    {
        s = s + this->residues(names[i], map);
    }

    return s;
}

Selector<Residue> MoleculeView::residues(const QList<qint64> &values,
                                         const PropertyMap &map) const
{
    if (values.count() == 0)
        throw SireError::invalid_index(QObject::tr(
            "You must specify some indexes.."), CODELOC);

    if (this->isA< Selector<Residue> >())
    {
        auto IDs = this->asA< Selector<Residue> >().indexes();

        QList<ResIdx> ret;
        ret.reserve(values.count());

        for (const auto &value : values)
        {
            ret.append(IDs.at(SireID::Index(value).map(IDs.count())));
        }

        return Selector<Residue>(this->data(), ret);
    }

    const auto s = this->selection();
    QList<ResIdx> idxs;

    if (s.selectedAllResidues())
    {
        for (auto value : values)
        {
            idxs.append(ResIdx(value));
        }
    }
    else
    {
        auto residxs = s.selectedResidues();

        for (auto value : values)
        {
            idxs.append(residxs.at(Index(value).map(residxs.count())));
        }
    }

    return Selector<Residue>(this->data(), idxs);
}

Selector<Residue> MoleculeView::residues(const Slice &slice,
                                         const PropertyMap &map) const
{
    if (this->isA< Selector<Residue> >())
    {
        auto IDs = this->asA< Selector<Residue> >().indexes();

        QList<ResIdx> ret;

        for (auto it = slice.begin(IDs.count()); not it.atEnd(); it.next())
        {
            ret.append(IDs.at(it.value()));
        }

        return Selector<Residue>(this->data(), ret);
    }

    const auto s = this->selection();
    QList<ResIdx> idxs;

    if (s.selectedAllResidues())
    {
        for (auto it = slice.begin(s.nResidues()); not it.atEnd(); it.next())
        {
            idxs.append(ResIdx(it.value()));
        }
    }
    else
    {
        auto residxs = s.selectedResidues();

        for (auto it = slice.begin(residxs.count()); not it.atEnd(); it.next())
        {
            idxs.append(residxs.at(Index(it.value()).map(residxs.count())));
        }
    }

    return Selector<Residue>(*this, idxs);
}

/** Return the residues from this view that match the ID 'resid'

    \throw SireMol::missing_residue
    \throw SireError::invalid_index
*/
Selector<Residue> MoleculeView::residues(const ResID &resid,
                                         const PropertyMap &map) const
{
    return resid.selectAllFrom(*this, map);
}

/** Return the residue that is part of this view

    \throw SireMol::missing_residue
    \throw SireMol::duplicate_residue
*/
Residue MoleculeView::residue() const
{
    QList<ResIdx> selected_res = this->selection().selectedResidues();

    if (selected_res.isEmpty())
        throw SireMol::missing_residue( QObject::tr(
                "This view does not contain any residues."), CODELOC );

    else if (selected_res.count() != 1)
        throw SireMol::duplicate_residue( QObject::tr(
                "Cannot convert this view (%1) into a Residue as "
                "we can only do this is just one Residue is selected. "
                "These Residues are selected ; %2.")
                    .arg(this->toString(), Sire::toString(selected_res)),
                        CODELOC );

    return Residue(this->data(), selected_res.at(0));
}

/** Return all of the residues that are involved with this view

    \throw SireMol::missing_residue
*/
Selector<Residue> MoleculeView::residues() const
{
    if (this->isA< Selector<Residue> >())
        return this->asA< Selector<Residue> >();

    QList<ResIdx> selected_res = this->selection().selectedResidues();

    if (selected_res.isEmpty())
    {
        throw SireMol::missing_residue( QObject::tr(
                "This view does not contain any residues."), CODELOC );
    }

    return Selector<Residue>(this->data(), selected_res);
}

Chain MoleculeView::chain(int i, const PropertyMap &map) const
{
    auto s = this->selection();

    if (s.selectedAllChains())
    {
        return this->chain(ChainIdx(i), map);
    }
    else
    {
        auto cidxs = s.selectedChains();
        return this->chain(cidxs.at(Index(i).map(cidxs.count())), map);
    }
}

Chain MoleculeView::chain(const QString &name, const PropertyMap &map) const
{
    try
    {
        return this->chain(ChainID::fromString(name), map);
    }
    catch(const SireMol::duplicate_chain &e)
    {
        throw e;
    }
    catch(const SireError::exception &e)
    {
        try
        {
            auto a = this->search(name).views().at(0).chain();
            return this->chain(a.index(), map);
        }
        catch(...)
        {
            if (name.length() < 5)
                //likely a name error
                e.throwSelf();
            else
                //likely a syntax error
                throw;
        }
    }

    return Chain();
}

/** Return the chain that is involved with this view that matches
    the ID 'chainid'

    \throw SireMol::missing_chain
    \throw SireError::invalid_index
    \throw SireMol::duplicate_chain
*/
Chain MoleculeView::chain(const ChainID &chainid, const PropertyMap &map) const
{
    return chainid.selectFrom(*this, map);
}

Selector<Chain> MoleculeView::chains(const QString &name,
                                     const PropertyMap &map) const
{
    try
    {
        return this->chains(ChainID::fromString(name), map);
    }
    catch(const SireError::exception &e)
    {
        try
        {
            const auto a = this->search(name).views();

            if (a.count() == 0)
                throw SireMol::missing_chain(QObject::tr(
                    "No chain matches '%1'").arg(name), CODELOC);

            return this->molecule().chains(
                        _toIndicies(a[0].chains().IDs()), map);
        }
        catch(...)
        {
            if (name.length() < 5)
                //likely a name error
                e.throwSelf();
            else
                //likely a syntax error
                throw;
        }
    }

    return Selector<Chain>();
}

Selector<Chain> MoleculeView::chains(const QStringList &names,
                                     const PropertyMap &map) const
{
    if (names.count() == 0)
        throw SireMol::missing_atom(QObject::tr(
            "You must specify some names.."), CODELOC);

    auto s = this->chains(names[0], map);

    for (int i=1; i<names.count(); ++i)
    {
        s = s + this->chains(names[i], map);
    }

    return s;
}

Selector<Chain> MoleculeView::chains(const QList<qint64> &values,
                                     const PropertyMap &map) const
{
    if (values.count() == 0)
        throw SireError::invalid_index(QObject::tr(
            "You must specify some indexes.."), CODELOC);

    if (this->isA< Selector<Chain> >())
    {
        auto IDs = this->asA< Selector<Chain> >().indexes();

        QList<ChainIdx> ret;
        ret.reserve(values.count());

        for (const auto &value : values)
        {
            ret.append(IDs.at(SireID::Index(value).map(IDs.count())));
        }

        return Selector<Chain>(this->data(), ret);
    }

    const auto s = this->selection();
    QList<ChainIdx> idxs;

    if (s.selectedAllChains())
    {
        for (auto value : values)
        {
            idxs.append(ChainIdx(value));
        }
    }
    else
    {
        auto chainidxs = s.selectedChains();

        for (auto value : values)
        {
            idxs.append(chainidxs.at(Index(value).map(chainidxs.count())));
        }
    }

    return Selector<Chain>(this->data(), idxs);
}

Selector<Chain> MoleculeView::chains(const Slice &slice,
                                     const PropertyMap &map) const
{
    if (this->isA< Selector<Chain> >())
    {
        auto IDs = this->asA< Selector<Chain> >().indexes();

        QList<ChainIdx> ret;

        for (auto it = slice.begin(IDs.count()); not it.atEnd(); it.next())
        {
            ret.append(IDs.at(it.value()));
        }

        return Selector<Chain>(this->data(), ret);
    }

    const auto s = this->selection();
    QList<ChainIdx> idxs;

    if (s.selectedAllChains())
    {

        for (auto it = slice.begin(s.nChains()); not it.atEnd(); it.next())
        {
            idxs.append(ChainIdx(it.value()));
        }
    }
    else
    {
        auto cidxs = s.selectedChains();

        for (auto it = slice.begin(cidxs.count()); not it.atEnd(); it.next())
        {
            idxs.append(cidxs.at(it.value()));
        }
    }

    return Selector<Chain>(*this, idxs);
}

/** Return the chains that are involved with this view that match
    the ID 'chainid'

    \throw SireMol::missing_chain
    \throw SireError::invalid_index
    \throw SireMol::duplicate_chain
*/
Selector<Chain> MoleculeView::chains(const ChainID &chainid,
                                     const PropertyMap &map) const
{
    return chainid.selectAllFrom(*this, map);
}

/** Return the chain that is involved with this view

    \throw SireMol::missing_chain
    \throw SireMol::duplicate_chain
*/
Chain MoleculeView::chain() const
{
    QList<ChainIdx> selected_chn = this->selection().selectedChains();

    if (selected_chn.isEmpty())
        throw SireMol::missing_chain( QObject::tr(
                "This view does not contain any chains."), CODELOC );

    else if (selected_chn.count() != 1)
        throw SireMol::duplicate_chain( QObject::tr(
                "Cannot convert this view (%1) into a Chain as "
                "we can only do this is just one Chain is selected. "
                "These Chains are selected ; %2.")
                    .arg(this->toString(), Sire::toString(selected_chn)),
                        CODELOC );

    return Chain(this->data(), selected_chn.at(0));
}

/** Return the chains that are involved with this view

    \throw SireMol::missing_chain
*/
Selector<Chain> MoleculeView::chains() const
{
    if (this->isA< Selector<Chain> >())
        return this->asA< Selector<Chain> >();

    QList<ChainIdx> selected_chn = this->selection().selectedChains();

    if (selected_chn.isEmpty())
        throw SireMol::missing_chain( QObject::tr(
                "This view does not contain any chains."), CODELOC );

    return Selector<Chain>(this->data(), selected_chn);
}

Segment MoleculeView::segment(int i, const PropertyMap &map) const
{
    auto s = this->selection();

    if (s.selectedAllSegments())
    {
        return this->segment(SegIdx(i), map);
    }
    else
    {
        auto segidxs = s.selectedSegments();
        return this->segment(segidxs.at(Index(i).map(segidxs.count())), map);
    }
}

Segment MoleculeView::segment(const QString &name, const PropertyMap &map) const
{
    try
    {
        return this->segment(SegID::fromString(name), map);
    }
    catch(const SireMol::duplicate_segment &e)
    {
        throw e;
    }
    catch(const SireError::exception &e)
    {
        try
        {
            auto a = this->search(name).views().at(0).segment();
            return this->segment(a.index(), map);
        }
        catch(...)
        {
            if (name.length() < 5)
                //likely a name error
                e.throwSelf();
            else
                //likely a syntax error
                throw;
        }
    }

    return Segment();
}

/** Return the segment that is involved with this view that matches
    the ID 'segid'

    \throw SireMol::missing_segment
    \throw SireError::invalid_index
    \throw SireMol::duplicate_segment
*/
Segment MoleculeView::segment(const SegID &segid, const PropertyMap &map) const
{
    return segid.selectFrom(*this, map);
}

Selector<Segment> MoleculeView::segments(const QString &name,
                                         const PropertyMap &map) const
{
    try
    {
        return this->segments(SegID::fromString(name), map);
    }
    catch(const SireError::exception &e)
    {
        try
        {
            const auto a = this->search(name).views();

            if (a.count() == 0)
                throw SireMol::missing_segment(QObject::tr(
                    "No segment matches '%1'").arg(name), CODELOC);

            return this->molecule().segments(
                        _toIndicies(a[0].segments().IDs()), map);
        }
        catch(...)
        {
            if (name.length() < 5)
                //likely a name error
                e.throwSelf();
            else
                //likely a syntax error
                throw;
        }
    }

    return Selector<Segment>();
}

Selector<Segment> MoleculeView::segments(const QStringList &names,
                                         const PropertyMap &map) const
{
    if (names.count() == 0)
        throw SireMol::missing_atom(QObject::tr(
            "You must specify some names.."), CODELOC);

    auto s = this->segments(names[0], map);

    for (int i=1; i<names.count(); ++i)
    {
        s = s + this->segments(names[i], map);
    }

    return s;
}

Selector<Segment> MoleculeView::segments(const QList<qint64> &values,
                                         const PropertyMap &map) const
{
    if (values.count() == 0)
        throw SireError::invalid_index(QObject::tr(
            "You must specify some indexes.."), CODELOC);

    if (this->isA< Selector<Segment> >())
    {
        auto IDs = this->asA< Selector<Segment> >().indexes();

        QList<SegIdx> ret;
        ret.reserve(values.count());

        for (const auto &value : values)
        {
            ret.append(IDs.at(SireID::Index(value).map(IDs.count())));
        }

        return Selector<Segment>(this->data(), ret);
    }

    const auto s = this->selection();
    QList<SegIdx> idxs;

    if (s.selectedAllSegments())
    {
        for (auto value : values)
        {
            idxs.append(SegIdx(value));
        }
    }
    else
    {
        auto segidxs = s.selectedSegments();

        for (auto value : values)
        {
            idxs.append(segidxs.at(Index(value).map(segidxs.count())));
        }
    }

    return Selector<Segment>(this->data(), idxs);
}

Selector<Segment> MoleculeView::segments(const Slice &slice,
                                         const PropertyMap &map) const
{
    if (this->isA< Selector<Segment> >())
    {
        auto IDs = this->asA< Selector<Segment> >().indexes();

        QList<SegIdx> ret;

        for (auto it = slice.begin(IDs.count()); not it.atEnd(); it.next())
        {
            ret.append(IDs.at(it.value()));
        }

        return Selector<Segment>(this->data(), ret);
    }

    const auto s = this->selection();
    QList<SegIdx> idxs;

    if (s.selectedAllSegments())
    {

        for (auto it = slice.begin(s.nSegments()); not it.atEnd(); it.next())
        {
            idxs.append(SegIdx(it.value()));
        }
    }
    else
    {
        auto segidxs = s.selectedSegments();

        for (auto it = slice.begin(segidxs.count()); not it.atEnd(); it.next())
        {
            idxs.append(segidxs.at(it.value()));
        }
    }

    return Selector<Segment>(*this, idxs);
}

/** Return the segments that are involved with this view that match
    the ID 'segid'

    \throw SireMol::missing_segment
    \throw SireError::invalid_index
*/
Selector<Segment> MoleculeView::segments(const SegID &segid, const PropertyMap &map) const
{
    return segid.selectAllFrom(*this, map);
}

/** Return the segment that is involved with this view

    \throw SireMol::missing_segment
    \throw SireMol::duplicate_segment
*/
Segment MoleculeView::segment() const
{
    QList<SegIdx> selected_seg = this->selection().selectedSegments();

    if (selected_seg.isEmpty())
        throw SireMol::missing_segment( QObject::tr(
                "This view does not contain any segments."), CODELOC );

    else if (selected_seg.count() != 1)
        throw SireMol::duplicate_segment( QObject::tr(
                "Cannot convert this view (%1) into a Segment as "
                "we can only do this is just one Segment is selected. "
                "These Segments are selected ; %2.")
                    .arg(this->toString(), Sire::toString(selected_seg)),
                        CODELOC );

    return Segment(this->data(), selected_seg.at(0));
}

/** Return the segments that are involved with this view

    \throw SireMol::missing_segment
*/
Selector<Segment> MoleculeView::segments() const
{
    if (this->isA< Selector<Segment> >())
        return this->asA< Selector<Segment> >();

    QList<SegIdx> selected_seg = this->selection().selectedSegments();

    if (selected_seg.isEmpty())
        throw SireMol::missing_segment( QObject::tr(
                "This view does not contain any segments."), CODELOC );

    return Selector<Segment>(this->data(), selected_seg);
}

/** Return the molecule involved with this view */
Molecule MoleculeView::molecule() const
{
    return Molecule(this->data());
}

/** Return the CutGroup whose atoms are in this view that matches
    the ID in 'cgid'

    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
    \throw SireMol::duplicate_cutgroup
*/
CutGroup MoleculeView::select(const CGID &cgid, const PropertyMap &map) const
{
    return this->cutGroup(cgid, map);
}

/** Return the residue from this view that matches the ID 'resid'

    \throw SireMol::missing_residue
    \throw SireError::invalid_index
    \throw SireMol::duplicate_residue
*/
Residue MoleculeView::select(const ResID &resid, const PropertyMap &map) const
{
    return this->residue(resid, map);
}

/** Return the chain that is involved with this view that matches
    the ID 'chainid'

    \throw SireMol::missing_chain
    \throw SireError::invalid_index
    \throw SireMol::duplicate_chain
*/
Chain MoleculeView::select(const ChainID &chainid, const PropertyMap &map) const
{
    return this->chain(chainid, map);
}

/** Return the segment that is involved with this view that matches
    the ID 'segid'

    \throw SireMol::missing_segment
    \throw SireError::invalid_index
    \throw SireMol::duplicate_segment
*/
Segment MoleculeView::select(const SegID &segid, const PropertyMap &map) const
{
    return this->segment(segid, map);
}

/** Return the atom in this view that matches the ID 'atomid'

    \throw SireMol::missing_atom
    \throw SireError::invalid_index
    \throw SireMol::duplicate_atom
*/
Atom MoleculeView::select(const AtomID &atomid, const PropertyMap &map) const
{
    return this->atom(atomid, map);
}

/** Return the result of searching this molecule using the passed
    search string */
SelectResult MoleculeView::search(const QString &search_string) const
{
    return Select(search_string)(*this);
}

/** Return the atoms from this view that match the ID 'atomid'

    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
Selector<Atom> MoleculeView::selectAll(const AtomID &atomid,
                                       const PropertyMap &map) const
{
    return this->atoms(atomid, map);
}

/** Return all of the atoms in this view

    \throw SireMol::missing_atom
*/
Selector<Atom> MoleculeView::selectAll() const
{
    return this->atoms();
}

/** Return all of the atoms in this view

    \throw SireMol::missing_atom
*/
Selector<Atom> MoleculeView::selectAllAtoms() const
{
    return this->atoms();
}

/** Return the CutGroups whose atoms are in this view that match
    the ID in 'cgid'

    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
*/
Selector<CutGroup> MoleculeView::selectAll(const CGID &cgid,
                                           const PropertyMap &map) const
{
    return this->cutGroups(cgid, map);
}

/** Return all of the CutGroups that are involved in this view

    \throw SireMol::missing_cutgroup
*/
Selector<CutGroup> MoleculeView::selectAllCutGroups() const
{
    return this->cutGroups();
}

/** Return the residues from this view that match the ID 'resid'

    \throw SireMol::missing_residue
    \throw SireError::invalid_index
*/
Selector<Residue> MoleculeView::selectAll(const ResID &resid,
                                          const PropertyMap &map) const
{
    return this->residues(resid, map);
}

/** Return all of the residues that are involved with this view

    \throw SireMol::missing_residue
*/
Selector<Residue> MoleculeView::selectAllResidues() const
{
    return this->residues();
}

/** Return the chains that are involved with this view that match
    the ID 'chainid'

    \throw SireMol::missing_chain
    \throw SireError::invalid_index
    \throw SireMol::duplicate_chain
*/
Selector<Chain> MoleculeView::selectAll(const ChainID &chainid,
                                        const PropertyMap &map) const
{
    return this->chains(chainid, map);
}

/** Return the chains that are involved with this view

    \throw SireMol::missing_chain
*/
Selector<Chain> MoleculeView::selectAllChains() const
{
    return this->chains();
}

/** Return the segments that are involved with this view that match
    the ID 'segid'

    \throw SireMol::missing_segment
    \throw SireError::invalid_index
*/
Selector<Segment> MoleculeView::selectAll(const SegID &segid,
                                          const PropertyMap &map) const
{
    return this->segments(segid, map);
}

/** Return the segments that are involved with this view

    \throw SireMol::missing_segment
*/
Selector<Segment> MoleculeView::selectAllSegments() const
{
    return this->segments();
}

/** Return a completely null molecule */
const MoleculeView& MoleculeView::null()
{
    return *(create_shared_null<Molecule>());
}

/** Return the number of atoms in this view */
int MoleculeView::nAtoms() const
{
    return this->selection().nSelectedAtoms();
}

/** Return the number of residues in this view */
int MoleculeView::nResidues() const
{
    return this->selection().nSelectedResidues();
}

/** Return the number of CutGroups in this view */
int MoleculeView::nCutGroups() const
{
    return this->selection().nSelectedCutGroups();
}

/** Return the number of chains in this view */
int MoleculeView::nChains() const
{
    return this->selection().nSelectedChains();
}

/** Return the number of segments in this view */
int MoleculeView::nSegments() const
{
    return this->selection().nSelectedSegments();
}

/** Return the number of sub-views in this view. Most MoleculeViews are
    just a single view, but some (e.g. ViewsOfMol or Selector<T>) have
    multiple views */
int MoleculeView::nViews() const
{
    return 1;
}

/** Return the number of sub-views in this view. Most MoleculeViews are
    just a single view, but some (e.g. ViewsOfMol or Selector<T>) have
    multiple views */
int MoleculeView::size() const
{
    return this->nViews();
}

/** Return the number of sub-views in this view. Most MoleculeViews are
    just a single view, but some (e.g. ViewsOfMol or Selector<T>) have
    multiple views */
int MoleculeView::count() const
{
    return this->nViews();
}

/** Expand this into a list of unit classes. This will return the view itself if
 *  this is a unit class, e.g. Atom, Residue, Molecule etc.
 *  This will return a list of unit classes if this is a Selector<T> or
 *  equivalent type class.
*/
QList<MolViewPtr> MoleculeView::toList() const
{
    return QList<MolViewPtr>({MolViewPtr(this->clone())});
}

/** Return the ith view in this MoleculeView. */
MolViewPtr MoleculeView::operator[](int i) const
{
    return this->atom(i);
}

/** Return the child with the specified name from this MoleculeView */
MolViewPtr MoleculeView::operator[](const QString &name) const
{
    return this->atom(name);
}

/** Return the atom(s) that match 'atomid' in this view of the molecule */
MolViewPtr MoleculeView::operator[](const AtomID &atomid) const
{
    auto atoms = this->selectAll(atomid);

    if (atoms.count() == 1)
    {
        return atoms[0];
    }
    else
    {
        return atoms;
    }
}

/** Return the residue(s) that match 'resid' in this view of the molecule */
MolViewPtr MoleculeView::operator[](const ResID &resid) const
{
    auto residues = this->selectAll(resid);

    if (residues.count() == 1)
    {
        return residues[0];
    }
    else
    {
        return residues;
    }
}

/** Return the CutGroups(s) that match 'resid' in this view of the molecule */
MolViewPtr MoleculeView::operator[](const CGID &cgid) const
{
    auto cutgroups = this->selectAll(cgid);

    if (cutgroups.count() == 1)
    {
        return cutgroups[0];
    }
    else
    {
        return cutgroups;
    }
}

/** Return the residue(s) that match 'resid' in this view of the molecule */
MolViewPtr MoleculeView::operator[](const ChainID &chainid) const
{
    auto chains = this->selectAll(chainid);

    if (chains.count() == 1)
    {
        return chains[0];
    }
    else
    {
        return chains;
    }
}

/** Return the residue(s) that match 'resid' in this view of the molecule */
MolViewPtr MoleculeView::operator[](const SegID &segid) const
{
    auto segments = this->selectAll(segid);

    if (segments.count() == 1)
    {
        return segments[0];
    }
    else
    {
        return segments;
    }
}

/** This is an overload of operator[](int), allowing a SireID::Index to be used
    as the int */
MolViewPtr MoleculeView::operator[](const SireID::Index &idx) const
{
    return this->operator[](idx.value());
}

MolViewPtr MoleculeView::operator[](const Slice &slice) const
{
    return this->atoms(slice);
}

MolViewPtr MoleculeView::operator[](const QList<qint64> &idxs) const
{
    return this->atoms(idxs);
}

MolViewPtr MoleculeView::at(int i) const
{
    return this->operator[](i);
}

MolViewPtr MoleculeView::at(const AtomID &atomid) const
{
    return this->operator[](atomid);
}

MolViewPtr MoleculeView::at(const ResID &resid) const
{
    return this->operator[](resid);
}

MolViewPtr MoleculeView::at(const CGID &cgid) const
{
    return this->operator[](cgid);
}

MolViewPtr MoleculeView::at(const ChainID &chainid) const
{
    return this->operator[](chainid);
}

MolViewPtr MoleculeView::at(const SegID &segid) const
{
    return this->operator[](segid);
}

MolViewPtr MoleculeView::at(const SireID::Index &idx) const
{
    return this->operator[](idx);
}

namespace SireBase
{
    template class SireBase::PropPtr<SireMol::MoleculeView>;
}
