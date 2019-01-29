/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2014  Christopher Woods
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

#include "atom.h"
#include "atommatcher.h"
#include "atommatchers.h"
#include "atomidx.h"
#include "atomname.h"
#include "atomselection.h"
#include "atomidentifier.h"
#include "evaluator.h"
#include "moleculeinfodata.h"
#include "moleculeview.h"
#include "mover.h"
#include "selector.hpp"

#include "tostring.h"

#include "SireMaths/vector.h"

#include "SireUnits/units.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"

using namespace SireMol;
using namespace SireUnits;
using namespace SireBase;
using namespace SireStream;

/////////
///////// Implmentation of AtomIdxMatcher
/////////

static const RegisterMetaType<AtomIdxMatcher> r_idxmatcher;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const AtomIdxMatcher &idxmatcher)
{
    writeHeader(ds, r_idxmatcher, 1);
    ds << static_cast<const AtomMatcher&>(idxmatcher);

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, AtomIdxMatcher &idxmatcher)
{
    VersionID v = readHeader(ds, r_idxmatcher);

    if (v == 1)
    {
        ds >> static_cast<AtomMatcher&>(idxmatcher);
    }
    else
        throw version_error(v, "1", r_idxmatcher, CODELOC);

    return ds;
}

/** Constructor */
AtomIdxMatcher::AtomIdxMatcher() : ConcreteProperty<AtomIdxMatcher,AtomMatcher>()
{}

/** Copy constructor */
AtomIdxMatcher::AtomIdxMatcher(const AtomIdxMatcher &other)
               : ConcreteProperty<AtomIdxMatcher,AtomMatcher>(other)
{}

/** Destructor */
AtomIdxMatcher::~AtomIdxMatcher()
{}

/** Copy assignment operator */
AtomIdxMatcher& AtomIdxMatcher::operator=(const AtomIdxMatcher &other)
{
    return *this;
}

/** Comparison operator */
bool AtomIdxMatcher::operator==(const AtomIdxMatcher &other) const
{
    return true;
}

/** Comparison operator */
bool AtomIdxMatcher::operator!=(const AtomIdxMatcher &other) const
{
    return false;
}

QString AtomIdxMatcher::toString() const
{
    return QObject::tr("AtomIdxMatcher()");
}

/** Match the atoms in 'mol1' to the atoms in 'mol0' - this
    returns the AtomIdxs of the atoms in 'mol1' that are in
    'mol0', indexed by the AtomIdx of the atom in 'mol0'.

     This skips atoms in 'mol1' that are not in 'mol0'
*/
QHash<AtomIdx,AtomIdx> AtomIdxMatcher::pvt_match(const MoleculeView &mol0,
                                                 const PropertyMap &map0,
                                                 const MoleculeView &mol1,
                                                 const PropertyMap &map1) const
{
    const AtomSelection sel0 = mol0.selection();
    const AtomSelection sel1 = mol1.selection();

    QHash<AtomIdx,AtomIdx> map;

    if (sel0.selectedAll() and sel1.selectedAll())
    {
        const int nats = qMin(sel0.nSelectedAtoms(), sel1.nSelectedAtoms());

        map.reserve(nats);

        for (int i=0; i<nats; ++i)
        {
            map.insert( AtomIdx(i), AtomIdx(i) );
        }
    }
    else
    {
        const int nats1 = mol1.data().info().nAtoms();

        foreach (const AtomIdx atom, sel0.selectedAtoms())
        {
            if (atom.value() < nats1)
            {
                if (sel1.selected(atom))
                    map.insert(atom, atom);
            }
        }
    }

    return map;
}

/** Match the atoms in 'mol1' to the atoms in 'mol0' - this
    returns the AtomIdxs of the atoms in 'mol1' that are in
    'mol0', indexed by the AtomIdx of the atom in 'mol0'.

     This skips atoms in 'mol1' that are not in 'mol0'
*/
QHash<AtomIdx,AtomIdx> AtomIdxMatcher::pvt_match(const MoleculeInfoData &mol0,
                                                 const MoleculeInfoData &mol1) const
{
    QHash<AtomIdx,AtomIdx> map;

    const int nats = qMin(mol0.nAtoms(), mol1.nAtoms());

    map.reserve(nats);

    for (int i=0; i<nats; ++i)
    {
        map.insert( AtomIdx(i), AtomIdx(i) );
    }

    return map;
}

/** The AtomIdx matcher does not change the order of the atoms */
bool AtomIdxMatcher::pvt_changesOrder(const MoleculeInfoData &molinfo0,
                                      const MoleculeInfoData &molinfo1) const
{
    return false;
}

const char* AtomIdxMatcher::typeName()
{
    return QMetaType::typeName( qMetaTypeId<AtomIdxMatcher>() );
}

/////////
///////// Implementation of AtomNameMatcher
/////////

static const RegisterMetaType<AtomNameMatcher> r_namematcher;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                       const AtomNameMatcher &namematcher)
{
    writeHeader(ds, r_namematcher, 1);
    ds << static_cast<const AtomMatcher&>(namematcher);

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, AtomNameMatcher &namematcher)
{
    VersionID v = readHeader(ds, r_namematcher);

    if (v == 1)
    {
        ds >> static_cast<AtomMatcher&>(namematcher);
    }
    else
        throw version_error(v, "1", r_namematcher, CODELOC);

    return ds;
}

/** Constructor */
AtomNameMatcher::AtomNameMatcher() : ConcreteProperty<AtomNameMatcher,AtomMatcher>()
{}

/** Copy constructor */
AtomNameMatcher::AtomNameMatcher(const AtomNameMatcher &other)
               : ConcreteProperty<AtomNameMatcher,AtomMatcher>(other)
{}

/** Destructor */
AtomNameMatcher::~AtomNameMatcher()
{}

/** Copy assignment operator */
AtomNameMatcher& AtomNameMatcher::operator=(const AtomNameMatcher &other)
{
    return *this;
}

/** Comparison operator */
bool AtomNameMatcher::operator==(const AtomNameMatcher &other) const
{
    return true;
}

/** Comparison operator */
bool AtomNameMatcher::operator!=(const AtomNameMatcher &other) const
{
    return false;
}

QString AtomNameMatcher::toString() const
{
    return QObject::tr("AtomNameMatcher()");
}

/** Match the atoms in 'mol1' to the atoms in 'mol0' - this
    returns the AtomIdxs of the atoms in 'mol1' that are in
    'mol0', indexed by the AtomIdx of the atom in 'mol0'.

     This skips atoms in 'mol1' that are not in 'mol0'
*/
QHash<AtomIdx,AtomIdx> AtomNameMatcher::pvt_match(const MoleculeView &mol0,
                                                  const PropertyMap &map0,
                                                  const MoleculeView &mol1,
                                                  const PropertyMap &map1) const
{
    QHash<AtomIdx,AtomIdx> map;

    const AtomSelection sel0 = mol0.selection();
    const AtomSelection sel1 = mol1.selection();

    if (sel0.selectedAll() and sel1.selectedAll())
    {
        for (int i=0; i<mol0.data().info().nAtoms(); ++i)
        {
            const AtomIdx idx0(i);

            const AtomName name = mol0.data().info().name(idx0);

            try
            {
                AtomIdx idx1 = mol1.data().info().atomIdx(name);
                map.insert( idx0, idx1 );
            }
            catch(...)
            {}
        }
    }
    else
    {
        foreach (const AtomIdx idx0, sel0.selectedAtoms())
        {
            const AtomName name0 = mol0.data().info().name(idx0);

            // A list of matches.
            QList<AtomIdx> matches;

            foreach (const AtomIdx idx1, sel1.selectedAtoms())
            {
                const AtomName name1 = mol1.data().info().name(idx1);

                // Add the match.
                if (name0 == name1)
                    matches.append(idx1);
            }

            // Only insert unique matche into the map.
            if (matches.count() == 1)
                map.insert(idx0, matches[0]);
        }
    }

    return map;
}

/** Match the atoms in 'mol1' to the atoms in 'mol0' - this
    returns the AtomIdxs of the atoms in 'mol1' that are in
    'mol0', indexed by the AtomIdx of the atom in 'mol0'.

     This skips atoms in 'mol1' that are not in 'mol0'
*/
QHash<AtomIdx,AtomIdx> AtomNameMatcher::pvt_match(const MoleculeInfoData &mol0,
                                                  const MoleculeInfoData &mol1) const
{
    QHash<AtomIdx,AtomIdx> map;

    for (int i=0; i<mol0.nAtoms(); ++i)
    {
        const AtomIdx idx0(i);

        const AtomName name = mol0.name(idx0);

        try
        {
            AtomIdx idx1 = mol1.atomIdx(name);
            map.insert( idx0, idx1 );
        }
        catch(...)
        {}
    }

    return map;
}

const char* AtomNameMatcher::typeName()
{
    return QMetaType::typeName( qMetaTypeId<AtomNameMatcher>() );
}

/////////
///////// Implementation of AtomMCSMatcher
/////////

static const RegisterMetaType<AtomMCSMatcher> r_mcsmatcher;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                       const AtomMCSMatcher &mcsmatcher)
{
    writeHeader(ds, r_mcsmatcher, 2);

    SharedDataStream sds(ds);

    sds << mcsmatcher.prematcher << double(mcsmatcher.t.to(second))
        << mcsmatcher.match_light << mcsmatcher.verbose
        << static_cast<const AtomMatcher&>(mcsmatcher);

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, AtomMCSMatcher &mcsmatcher)
{
    VersionID v = readHeader(ds, r_mcsmatcher);

    if (v == 2)
    {
        double timeout;

        SharedDataStream sds(ds);

        sds >> mcsmatcher.prematcher >> timeout
            >> mcsmatcher.match_light >> mcsmatcher.verbose
            >> static_cast<AtomMatcher&>(mcsmatcher);

        mcsmatcher.t = timeout*second;
    }
    else if (v == 1)
    {
        double timeout;

        SharedDataStream sds(ds);

        sds >> mcsmatcher.prematcher >> timeout >> static_cast<AtomMatcher&>(mcsmatcher);

        mcsmatcher.t = timeout*second;
        mcsmatcher.match_light = false;
        mcsmatcher.verbose = false;
    }
    else
        throw version_error(v, "1,2", r_mcsmatcher, CODELOC);

    return ds;
}

/** Constructor */
AtomMCSMatcher::AtomMCSMatcher()
               : ConcreteProperty<AtomMCSMatcher,AtomMatcher>(), t(1*second),
                 match_light(false),
                 verbose(false)
{}

/** Constructor */
AtomMCSMatcher::AtomMCSMatcher(bool verbose)
               : ConcreteProperty<AtomMCSMatcher,AtomMatcher>(), t(1*second),
                 match_light(false),
                 verbose(verbose)
{}

/** Construct specifying the timeout for the MCS match */
AtomMCSMatcher::AtomMCSMatcher(const SireUnits::Dimension::Time &timeout,
                               bool verbose)
               : ConcreteProperty<AtomMCSMatcher,AtomMatcher>(), t(timeout),
                 match_light(false),
                 verbose(verbose)
{}

/** Construct specifying the prematcher for the MCS match */
AtomMCSMatcher::AtomMCSMatcher(const AtomMatcher &matcher,
                               bool verbose)
               : ConcreteProperty<AtomMCSMatcher,AtomMatcher>(),
                 prematcher(matcher), t(1*second),
                 match_light(false), verbose(verbose)
{}

/** Construct specifying the timeout and prematcher for the MCS match */
AtomMCSMatcher::AtomMCSMatcher(const AtomMatcher &matcher,
                               const SireUnits::Dimension::Time &timeout,
                               bool verbose)
               : ConcreteProperty<AtomMCSMatcher,AtomMatcher>(),
                 prematcher(matcher), t(timeout),
                 match_light(false), verbose(verbose)
{}

/** Constructor, specifying whether or not to match light atoms */
AtomMCSMatcher::AtomMCSMatcher(bool match_light_atoms,
                               bool verbose)
               : ConcreteProperty<AtomMCSMatcher,AtomMatcher>(), t(1*second),
                 match_light(match_light_atoms),
                 verbose(verbose)
{}

/** Construct specifying the timeout for the MCS match, and specifying whether or not
    to match light atoms */
AtomMCSMatcher::AtomMCSMatcher(const SireUnits::Dimension::Time &timeout,
                               bool match_light_atoms,
                               bool verbose)
               : ConcreteProperty<AtomMCSMatcher,AtomMatcher>(), t(timeout),
                 match_light(match_light_atoms),
                 verbose(verbose)
{}

/** Construct specifying the prematcher for the MCS match,
    and specifying whether or not to match light atoms
*/
AtomMCSMatcher::AtomMCSMatcher(const AtomMatcher &matcher,
                               bool match_light_atoms,
                               bool verbose)
               : ConcreteProperty<AtomMCSMatcher,AtomMatcher>(),
                 prematcher(matcher), t(1*second),
                 match_light(match_light_atoms), verbose(verbose)
{}

/** Construct specifying the timeout and prematcher for the MCS match,
    and specifying whether or not to match light atoms
*/
AtomMCSMatcher::AtomMCSMatcher(const AtomMatcher &matcher,
                               const SireUnits::Dimension::Time &timeout,
                               bool match_light_atoms,
                               bool verbose)
               : ConcreteProperty<AtomMCSMatcher,AtomMatcher>(),
                 prematcher(matcher), t(timeout),
                 match_light(match_light_atoms), verbose(verbose)
{}

/** Copy constructor */
AtomMCSMatcher::AtomMCSMatcher(const AtomMCSMatcher &other)
               : ConcreteProperty<AtomMCSMatcher,AtomMatcher>(other),
                 prematcher(other.prematcher), t(other.t),
                 match_light(other.match_light), verbose(other.verbose)
{}

/** Destructor */
AtomMCSMatcher::~AtomMCSMatcher()
{}

/** Copy assignment operator */
AtomMCSMatcher& AtomMCSMatcher::operator=(const AtomMCSMatcher &other)
{
    if (this != &other)
    {
        t = other.t;
        prematcher = other.prematcher;
        match_light = other.match_light;
        verbose = other.verbose;
    }

    return *this;
}

/** Comparison operator */
bool AtomMCSMatcher::operator==(const AtomMCSMatcher &other) const
{
    return prematcher == other.prematcher and t == other.t and match_light == other.match_light and verbose == other.verbose;
}

/** Comparison operator */
bool AtomMCSMatcher::operator!=(const AtomMCSMatcher &other) const
{
    return not operator==(other);
}

QString AtomMCSMatcher::toString() const
{
    if (prematcher.isNull() or prematcher.read().isNull())
    {
        return QObject::tr("AtomMCSMatcher( timeout() = %1 s, matchingLightAtoms() = %2, "
                "isVerbose() = %3 )")
                .arg(t.to(second)).arg(match_light).arg(verbose);
    }
    else
    {
        return QObject::tr("AtomMCSMatcher( preMatcher() = %1, timeout() = %2 s, "
                           "matchingLightAtoms() = %3, isVerbose() = %4 )")
                    .arg(prematcher.read().toString())
                    .arg(t.to(second))
                    .arg(match_light)
                    .arg(verbose);
    }
}

const char* AtomMCSMatcher::typeName()
{
    return QMetaType::typeName( qMetaTypeId<AtomMCSMatcher>() );
}

/** Return the prematcher (if any) that is used to pre-match atoms
    before the MCS match */
const AtomMatcher& AtomMCSMatcher::preMatcher() const
{
    return prematcher.read();
}

/** Return the timeout before the MCS match is abandoned */
SireUnits::Dimension::Time AtomMCSMatcher::timeout() const
{
    return t;
}

/** Return whether or not this will include light atoms (e.g. hydrogen)
    when searching for the maximum common substructure */
bool AtomMCSMatcher::matchingLightAtoms() const
{
    return match_light;
}

/** Return whether or not this will report progress to stdout. */
bool AtomMCSMatcher::isVerbose() const
{
    return verbose;
}

/** Match the atoms in 'mol1' to the atoms in 'mol0' - this
    returns the AtomIdxs of the atoms in 'mol1' that are in
    'mol0', indexed by the AtomIdx of the atom in 'mol0'.

     This skips atoms in 'mol1' that are not in 'mol0'
*/
QHash<AtomIdx,AtomIdx> AtomMCSMatcher::pvt_match(const MoleculeView &mol0,
                                                 const PropertyMap &map0,
                                                 const MoleculeView &mol1,
                                                 const PropertyMap &map1) const
{
    if (prematcher.isNull() or prematcher.read().isNull())
        return Evaluator(mol0).findMCS(mol1, t, match_light, map0, map1, this-verbose);
    else
        return Evaluator(mol0).findMCS(mol1, prematcher.read(), t, match_light, map0, map1, this-verbose);
}

/////////
///////// Implementation of AtomIDMatcher
/////////

static const RegisterMetaType<AtomIDMatcher> r_idmatcher;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                       const AtomIDMatcher &idmatcher)
{
    writeHeader(ds, r_idmatcher, 1);

    SharedDataStream sds(ds);
    sds << idmatcher.m << static_cast<const AtomMatcher&>(idmatcher);

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, AtomIDMatcher &idmatcher)
{
    VersionID v = readHeader(ds, r_idmatcher);

    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> idmatcher.m >> static_cast<AtomMatcher&>(idmatcher);
    }
    else
        throw version_error(v, "1", r_idmatcher, CODELOC);

    return ds;
}

/** Constructor */
AtomIDMatcher::AtomIDMatcher() : ConcreteProperty<AtomIDMatcher,AtomMatcher>()
{}

/** Construct to match atom names */
AtomIDMatcher::AtomIDMatcher(const QList< QPair<QString,QString> > &names)
              : ConcreteProperty<AtomIDMatcher,AtomMatcher>()
{
    for (int i=0; i<names.count(); ++i)
    {
        const QPair<QString,QString> &name = names[i];

        if (not (name.first.isEmpty() or name.second.isEmpty()))
        {
            m.append( QPair<AtomIdentifier,AtomIdentifier>( AtomName(name.first),
                                                            AtomName(name.second) ) );
        }
    }
}

/** Construct to match atom names */
AtomIDMatcher::AtomIDMatcher(const QHash<QString,QString> &names)
              : ConcreteProperty<AtomIDMatcher,AtomMatcher>()
{
    QSet<QString> matched_names;

    for (QHash<QString,QString>::const_iterator it = names.constBegin();
         it != names.constEnd();
         ++it)
    {
        if (not it.key().isEmpty() or it.value().isEmpty())
        {
            if (matched_names.contains(it.value()))
                throw SireError::invalid_arg( QObject::tr(
                        "You are trying to match multiple atoms (%1) to the same name (%2). "
                        "Please ensure that you have a unique name <=> name mapping.")
                            .arg(Sire::toString(names.keys(it.value())))
                            .arg(it.value()), CODELOC );

            matched_names.insert(it.value());
            m.append( QPair<AtomIdentifier,AtomIdentifier>( AtomName(it.key()),
                                                            AtomName(it.value()) ) );
        }
    }
}

/** Construct to match atom indexes */
AtomIDMatcher::AtomIDMatcher(const QList< QPair<int,int> > &idxs)
              : ConcreteProperty<AtomIDMatcher,AtomMatcher>()
{
    for (int i=0; i<idxs.count(); ++i)
    {
        const QPair<int,int> &idx = idxs[i];

        m.append( QPair<AtomIdentifier,AtomIdentifier>( AtomIdx(idx.first),
                                                        AtomIdx(idx.second) ) );
    }
}

/** Construct to match atom indexes */
AtomIDMatcher::AtomIDMatcher(const QHash<int,int> &idxs)
              : ConcreteProperty<AtomIDMatcher,AtomMatcher>()
{
    QSet<int> matched_idxs;

    for (QHash<int,int>::const_iterator it = idxs.constBegin();
         it != idxs.constEnd();
         ++it)
    {
        if (matched_idxs.contains(it.value()))
            throw SireError::invalid_arg( QObject::tr(
                    "You are trying to match multiple atoms (%1) to the same index (%2). "
                    "Please ensure that you have a unique index <=> index mapping.")
                        .arg(Sire::toString(idxs.keys(it.value())))
                        .arg(it.value()), CODELOC );

        matched_idxs.insert(it.value());

        m.append( QPair<AtomIdentifier,AtomIdentifier>( AtomIdx(it.key()),
                                                        AtomIdx(it.value()) ) );
    }
}

/** Construct to match specified AtomIdentifiers */
AtomIDMatcher::AtomIDMatcher(const QList< QPair<AtomIdentifier,AtomIdentifier> > &ids)
              : ConcreteProperty<AtomIDMatcher,AtomMatcher>()
{
    for (int i=0; i<ids.count(); ++i)
    {
        const QPair<AtomIdentifier,AtomIdentifier> id = ids[i];

        if (not (id.first.isNull() or id.second.isNull()))
            m.append(id);
    }
}

/** Construct to match specified AtomIdentifiers */
AtomIDMatcher::AtomIDMatcher(const QHash<AtomIdentifier,AtomIdentifier> &ids)
              : ConcreteProperty<AtomIDMatcher,AtomMatcher>()
{
    QSet<AtomIdentifier> matched_ids;

    for (QHash<AtomIdentifier,AtomIdentifier>::const_iterator it = ids.constBegin();
         it != ids.constEnd();
         ++it)
    {
        if (not it.key().isNull() or it.value().isNull())
        {
            if (matched_ids.contains(it.value()))
                throw SireError::invalid_arg( QObject::tr(
                        "You are trying to match multiple atoms (%1) to the same ID (%2). "
                        "Please ensure that you have a unique ID <=> ID mapping.")
                            .arg(Sire::toString(ids.keys(it.value())))
                            .arg(it.value().toString()), CODELOC );

            matched_ids.insert(it.value());

            m.append( QPair<AtomIdentifier,AtomIdentifier>(it.key(), it.value()) );
        }
    }
}

/** Construct to match atom names */
AtomIDMatcher::AtomIDMatcher(const QList< boost::tuple<QString,QString> > &names)
              : ConcreteProperty<AtomIDMatcher,AtomMatcher>()
{
    for (int i=0; i<names.count(); ++i)
    {
        const boost::tuple<QString,QString> &name = names[i];

        if (not (name.get<0>().isEmpty() or name.get<1>().isEmpty()))
        {
            m.append( QPair<AtomIdentifier,AtomIdentifier>( AtomName(name.get<0>()),
                                                            AtomName(name.get<1>()) ) );
        }
    }
}

/** Shorthand to construct to match atom names from the passed single string,
    with format 'atom0:atom1,atom2:atom3' etc. (i.e. comma separated pairs,
    each pair is colon separated to match atom to atom, e.g. this string
    matches atom0 to atom1, and atom2 to atom3) */
AtomIDMatcher::AtomIDMatcher(const QString &atom_names)
              : ConcreteProperty<AtomIDMatcher,AtomMatcher>()
{
    auto words = atom_names.split(",");

    if (not words.isEmpty())
    {
        QHash<QString,QString> match_names;
        match_names.reserve(words.count());

        for (auto word : words)
        {
            auto atoms = word.split(":");

            if (atoms.count() == 2)
            {
                match_names.insert( atoms[0].simplified(), atoms[1].simplified() );
            }
        }

        this->operator=( AtomIDMatcher(match_names) );
    }
}

/** Construct to match atom indexes */
AtomIDMatcher::AtomIDMatcher(const QList< boost::tuple<int,int> > &idxs)
              : ConcreteProperty<AtomIDMatcher,AtomMatcher>()
{
    for (int i=0; i<idxs.count(); ++i)
    {
        const boost::tuple<int,int> &idx = idxs[i];

        m.append( QPair<AtomIdentifier,AtomIdentifier>( AtomIdx(idx.get<0>()),
                                                        AtomIdx(idx.get<1>()) ) );
    }
}

/** Construct to match specified AtomIdentifiers */
AtomIDMatcher::AtomIDMatcher(const QList< boost::tuple<AtomIdentifier,AtomIdentifier> > &ids)
              : ConcreteProperty<AtomIDMatcher,AtomMatcher>()
{
    for (int i=0; i<ids.count(); ++i)
    {
        const boost::tuple<AtomIdentifier,AtomIdentifier> id = ids[i];

        if (not (id.get<0>().isNull() or id.get<1>().isNull()))
            m.append( QPair<AtomIdentifier,AtomIdentifier>(id.get<0>(),id.get<1>()) );
    }
}

/** Copy constructor */
AtomIDMatcher::AtomIDMatcher(const AtomIDMatcher &other)
              : ConcreteProperty<AtomIDMatcher,AtomMatcher>(other),
                m(other.m)
{}

/** Destructor */
AtomIDMatcher::~AtomIDMatcher()
{}

/** Copy assignment operator */
AtomIDMatcher& AtomIDMatcher::operator=(const AtomIDMatcher &other)
{
    m = other.m;
    return *this;
}

/** Comparison operator */
bool AtomIDMatcher::operator==(const AtomIDMatcher &other) const
{
    return m == other.m;
}

/** Comparison operator */
bool AtomIDMatcher::operator!=(const AtomIDMatcher &other) const
{
    return not operator==(*this);
}

/** Return whether or not this matcher is null (cannot be used for matching) */
bool AtomIDMatcher::isNull() const
{
    return m.isEmpty();
}

QString AtomIDMatcher::toString() const
{
    if (isNull())
        return QObject::tr("AtomIDMatcher::null");

    QStringList matches;

    for (QList< QPair<AtomIdentifier,AtomIdentifier> >::const_iterator it = m.constBegin();
         it != m.constEnd();
         ++it)
    {
        matches.append( QObject::tr(" %1 <=> %2").arg(it->first.toString(),
                                                      it->second.toString()) );
    }

    return QObject::tr("AtomIDMatcher( %1 )").arg( matches.join("\n") );
}

/** Match the atoms in 'mol1' to the atoms in 'mol0' - this
    returns the AtomIdxs of the atoms in 'mol1' that are in
    'mol0', indexed by the AtomIdx of the atom in 'mol0'.

     This skips atoms in 'mol1' that are not in 'mol0'
*/
QHash<AtomIdx,AtomIdx> AtomIDMatcher::pvt_match(const MoleculeView &mol0,
                                                const PropertyMap &map0,
                                                const MoleculeView &mol1,
                                                const PropertyMap &map1) const
{
    QHash<AtomIdx,AtomIdx> map;

    const AtomSelection sel0 = mol0.selection();
    const AtomSelection sel1 = mol1.selection();

    QSet<AtomIdx> found_1;

    const int nats0 = sel0.nSelectedAtoms();
    const int nats1 = sel1.nSelectedAtoms();

    for (QList< QPair<AtomIdentifier,AtomIdentifier> >::const_iterator it = m.constBegin();
         it != m.constEnd();
         ++it)
    {
        //find the atom in mol0
        try
        {
            AtomIdx idx0 = mol0.data().info().atomIdx(it->first);

            if (sel0.selected(idx0) and not map.contains(idx0))
            {
                try
                {
                    AtomIdx idx1 = mol1.data().info().atomIdx(it->second);
                    Atom atom = mol1.atom(idx1);

                    if (sel1.selected(idx1) and not found_1.contains(idx1))
                    {
                        map.insert(idx0, idx1);
                        found_1.insert(idx1);
                    }
                }
                catch(...)
                {}
            }
        }
        catch(...)
        {}

        if (map.count() >= nats0 or found_1.count() >= nats1)
            break;
    }

    return map;
}

/** Match the atoms in 'mol1' to the atoms in 'mol0' - this
    returns the AtomIdxs of the atoms in 'mol1' that are in
    'mol0', indexed by the AtomIdx of the atom in 'mol0'.

     This skips atoms in 'mol1' that are not in 'mol0'
*/
QHash<AtomIdx,AtomIdx> AtomIDMatcher::pvt_match(const MoleculeInfoData &mol0,
                                                const MoleculeInfoData &mol1) const
{
    QHash<AtomIdx,AtomIdx> map;


    QSet<AtomIdx> found_1;

    const int nats0 = mol0.nAtoms();
    const int nats1 = mol1.nAtoms();

    for (QList< QPair<AtomIdentifier,AtomIdentifier> >::const_iterator it = m.constBegin();
         it != m.constEnd();
         ++it)
    {
        //find the atom in mol0
        try
        {
            AtomIdx idx0 = mol0.atomIdx(it->first);

            if (not map.contains(idx0))
            {
                try
                {
                    AtomIdx idx1 = mol1.atomIdx(it->second);

                    if (not found_1.contains(idx1))
                    {
                        map.insert(idx0, idx1);
                        found_1.insert(idx1);
                    }
                }
                catch(...)
                {}
            }
        }
        catch(...)
        {}

        if (map.count() >= nats0 or found_1.count() >= nats1)
            break;
    }

    return map;
}

const char* AtomIDMatcher::typeName()
{
    return QMetaType::typeName( qMetaTypeId<AtomIDMatcher>() );
}

/////////
///////// Implementation of AtomMultiMatcher
/////////

static const RegisterMetaType<AtomMultiMatcher> r_multimatcher;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                       const AtomMultiMatcher &multimatcher)
{
    writeHeader(ds, r_multimatcher, 1);

    SharedDataStream sds(ds);
    sds << multimatcher.m << static_cast<const AtomMatcher&>(multimatcher);

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, AtomMultiMatcher &multimatcher)
{
    VersionID v = readHeader(ds, r_multimatcher);

    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> multimatcher.m >> static_cast<AtomMatcher&>(multimatcher);
    }
    else
        throw version_error(v, "1", r_multimatcher, CODELOC);

    return ds;
}

/** Constructor */
AtomMultiMatcher::AtomMultiMatcher() : ConcreteProperty<AtomMultiMatcher,AtomMatcher>()
{}

/** Construct from a single match */
AtomMultiMatcher::AtomMultiMatcher(const AtomMatcher &matcher)
                 : ConcreteProperty<AtomMultiMatcher,AtomMatcher>()
{
    if (matcher.isA<AtomMultiMatcher>())
    {
        this->operator=(matcher.asA<AtomMultiMatcher>());
    }
    else
    {
        m.append( matcher );
    }
}

/** Construct from a pair of matches */
AtomMultiMatcher::AtomMultiMatcher(const AtomMatcher &m0, const AtomMatcher &m1)
                 : ConcreteProperty<AtomMultiMatcher,AtomMatcher>()
{
    this->operator=( AtomMultiMatcher(m0) );

    if (m1.isA<AtomMultiMatcher>())
    {
        m += m1.asA<AtomMultiMatcher>().m;
    }
    else
    {
        m += m1;
    }
}

/** Copy constructor */
AtomMultiMatcher::AtomMultiMatcher(const AtomMultiMatcher &other)
                 : ConcreteProperty<AtomMultiMatcher,AtomMatcher>(other),
                   m(other.m)
{}

/** Destructor */
AtomMultiMatcher::~AtomMultiMatcher()
{}

/** Copy assignment operator */
AtomMultiMatcher& AtomMultiMatcher::operator=(const AtomMultiMatcher &other)
{
    m = other.m;
    return *this;
}

/** Comparison operator */
bool AtomMultiMatcher::operator==(const AtomMultiMatcher &other) const
{
    return m == other.m;
}

/** Comparison operator */
bool AtomMultiMatcher::operator!=(const AtomMultiMatcher &other) const
{
    return not operator==(*this);
}

/** Return whether or not this matcher is null (cannot be used for matching) */
bool AtomMultiMatcher::isNull() const
{
    return m.isEmpty();
}

QString AtomMultiMatcher::toString() const
{
    if (isNull())
        return QObject::tr("AtomMultiMatcher::null");

    QStringList matches;

    for (QList<AtomMatcherPtr>::const_iterator it = m.constBegin();
         it != m.constEnd();
         ++it)
    {
        matches.append( QObject::tr(" %1").arg(it->read().toString()) );
    }

    return QObject::tr("AtomMultiMatcher{ %1 }").arg( matches.join(" +\n") );
}

/** Match the atoms in 'mol1' to the atoms in 'mol0' - this
    returns the AtomIdxs of the atoms in 'mol1' that are in
    'mol0', indexed by the AtomIdx of the atom in 'mol0'.

     This skips atoms in 'mol1' that are not in 'mol0'
*/
QHash<AtomIdx,AtomIdx> AtomMultiMatcher::pvt_match(const MoleculeView &mol0,
                                                   const PropertyMap &map0,
                                                   const MoleculeView &mol1,
                                                   const PropertyMap &map1) const
{
    if (m.isEmpty())
        return QHash<AtomIdx,AtomIdx>();
    else if (m.count() == 1)
        return m.at(0).read().match(mol0,mol1);

    QHash<AtomIdx,AtomIdx> map;

    QSet<AtomIdx> found_1;
    const int nats0 = mol0.selection().nSelected();
    const int nats1 = mol1.selection().nSelected();

    for (QList<AtomMatcherPtr>::const_iterator it = m.constBegin();
         it != m.constEnd();
         ++it)
    {
        QHash<AtomIdx,AtomIdx> lmap = it->read().match(mol0, map0, mol1, map1);

        for (QHash<AtomIdx,AtomIdx>::const_iterator it2 = lmap.constBegin();
             it2 != lmap.constEnd();
             ++it2)
        {
            if (not (map.contains(it2.key()) or found_1.contains(it2.value())))
            {
                map.insert( it2.key(), it2.value() );
                found_1.insert( it2.value() );
            }
        }

        if (map.count() == nats0)
            break;

        if (found_1.count() == nats1)
            break;

        if (map.count() > nats0)
            throw SireError::program_bug( QObject::tr(
                    "Should not have excess matched atoms???"), CODELOC );

        if (found_1.count() > nats1)
            throw SireError::program_bug( QObject::tr(
                    "Should not have excess matched atoms???"), CODELOC );
    }

    return map;
}

/** Match the atoms in 'mol1' to the atoms in 'mol0' - this
    returns the AtomIdxs of the atoms in 'mol1' that are in
    'mol0', indexed by the AtomIdx of the atom in 'mol0'.

     This skips atoms in 'mol1' that are not in 'mol0'
*/
QHash<AtomIdx,AtomIdx> AtomMultiMatcher::pvt_match(const MoleculeInfoData &mol0,
                                                   const MoleculeInfoData &mol1) const
{
    if (m.isEmpty())
        return QHash<AtomIdx,AtomIdx>();
    else if (m.count() == 1)
        return m.at(0).read().match(mol0,mol1);

    QHash<AtomIdx,AtomIdx> map;

    QSet<AtomIdx> found_1;
    const int nats0 = mol0.nAtoms();
    const int nats1 = mol1.nAtoms();

    for (QList<AtomMatcherPtr>::const_iterator it = m.constBegin();
         it != m.constEnd();
         ++it)
    {
        QHash<AtomIdx,AtomIdx> lmap = it->read().match(mol0, mol1);

        for (QHash<AtomIdx,AtomIdx>::const_iterator it2 = lmap.constBegin();
             it2 != lmap.constEnd();
             ++it2)
        {
            if (not (map.contains(it2.key()) or found_1.contains(it2.value())))
            {
                map.insert( it2.key(), it2.value() );
                found_1.insert( it2.value() );
            }
        }

        if (map.count() == nats0)
            break;

        if (found_1.count() == nats1)
            break;

        if (map.count() > nats0)
            throw SireError::program_bug( QObject::tr(
                    "Should not have excess matched atoms???"), CODELOC );

        if (found_1.count() > nats1)
            throw SireError::program_bug( QObject::tr(
                    "Should not have excess matched atoms???"), CODELOC );
    }

    return map;
}

const char* AtomMultiMatcher::typeName()
{
    return QMetaType::typeName( qMetaTypeId<AtomMultiMatcher>() );
}

/////////
///////// Implementation of ResIdxAtomNameMatcher
/////////

static const RegisterMetaType<ResIdxAtomNameMatcher> r_residxatomnamematcher;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                       const ResIdxAtomNameMatcher &residxatomnamematcher)
{
    writeHeader(ds, r_residxatomnamematcher, 1);
    ds << static_cast<const ResIdxAtomNameMatcher&>(residxatomnamematcher);

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, ResIdxAtomNameMatcher &residxatomnamematcher)
{
    VersionID v = readHeader(ds, r_residxatomnamematcher);

    if (v == 1)
    {
        ds >> static_cast<ResIdxAtomNameMatcher&>(residxatomnamematcher);
    }
    else
        throw version_error(v, "1", r_residxatomnamematcher, CODELOC);

    return ds;
}

/** Constructor */
ResIdxAtomNameMatcher::ResIdxAtomNameMatcher() : ConcreteProperty<ResIdxAtomNameMatcher,AtomMatcher>()
{}

/** Copy constructor */
ResIdxAtomNameMatcher::ResIdxAtomNameMatcher(const ResIdxAtomNameMatcher &other)
               : ConcreteProperty<ResIdxAtomNameMatcher,AtomMatcher>(other)
{}

/** Destructor */
ResIdxAtomNameMatcher::~ResIdxAtomNameMatcher()
{}

/** Copy assignment operator */
ResIdxAtomNameMatcher& ResIdxAtomNameMatcher::operator=(const ResIdxAtomNameMatcher &other)
{
    return *this;
}

/** Comparison operator */
bool ResIdxAtomNameMatcher::operator==(const ResIdxAtomNameMatcher &other) const
{
    return true;
}

/** Comparison operator */
bool ResIdxAtomNameMatcher::operator!=(const ResIdxAtomNameMatcher &other) const
{
    return false;
}

QString ResIdxAtomNameMatcher::toString() const
{
    return QObject::tr("ResIdxAtomNameMatcher()");
}

/** Match the atoms in 'mol1' to the atoms in 'mol0' by name, searching each
    residue separately (by index) and combining the results. This returns the
    AtomIdxs of the atoms in 'mol1' that are in 'mol0', indexed by the AtomIdx
    of the atom in 'mol0'.

    This skips atoms in 'mol1' that are not in 'mol0'
*/
QHash<AtomIdx,AtomIdx> ResIdxAtomNameMatcher::pvt_match(const MoleculeView &mol0,
                                                        const PropertyMap &map0,
                                                        const MoleculeView &mol1,
                                                        const PropertyMap &map1) const
{
    const AtomSelection sel0 = mol0.selection();
    const AtomSelection sel1 = mol1.selection();

    if (sel0.selectedAll() and sel1.selectedAll())
    {
        return pvt_match(mol0.data().info(), mol1.data().info());
    }
    else
    {
        throw SireError::unsupported(QObject::tr("ResIdxAtomNameMatcher only works with "
            "full molecule selections."));
    }
}

/** Match the atoms in 'mol1' to the atoms in 'mol0' by name, searching each
    residue separately (by index) and combining the results. This returns the
    AtomIdxs of the atoms in 'mol1' that are in 'mol0', indexed by the AtomIdx
    of the atom in 'mol0'.

    This skips atoms in 'mol1' that are not in 'mol0'. Note that we only match
    the first unique atom within each residue.
*/
QHash<AtomIdx,AtomIdx> ResIdxAtomNameMatcher::pvt_match(const MoleculeInfoData &mol0,
                                                        const MoleculeInfoData &mol1) const
{
    QHash<AtomIdx,AtomIdx> map;

    // Get the list of residue indices from the reference molecule.
    auto resIdxs = mol0.getResidues();

    // Loop over all of the residues.
    for (const auto &resIdx : resIdxs)
    {
        // Get a list of atoms for the residue for both molecules.
        auto atoms0 = mol0.getAtomsIn(resIdx);
        auto atoms1 = mol1.getAtomsIn(resIdx);

        // The set of matched atoms from 'mol1' that belong to this residue.
        QSet<AtomIdx> matched;

        // For each atom from mol0, find those atoms in mol1 with the same name.
        for (const auto &idx0 : atoms0)
        {
            // Get the name of the atom.
            const auto name0 = mol0.name(idx0);

            // Loop over the atoms in mol1 until we find the first unique match.
            for (const auto &idx1 : atoms1)
            {
                // Get the name of the atom.
                const auto name1 = mol1.name(idx1);

                // This is a match, and it hasn't already been matched
                // to another atom in the residue.
                if ((name0 == name1) and (not matched.contains(idx1)))
                {
                    map.insert(idx0, idx1);
                    matched.insert(idx1);
                    break;
                }
            }
        }
    }

    return map;
}

const char* ResIdxAtomNameMatcher::typeName()
{
    return QMetaType::typeName( qMetaTypeId<ResIdxAtomNameMatcher>() );
}

/////////
///////// Implementation of ResIdxAtomMCSMatcher
/////////

static const RegisterMetaType<ResIdxAtomMCSMatcher> r_residxmcsmatcher;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                       const ResIdxAtomMCSMatcher &residxmcsmatcher)
{
    writeHeader(ds, r_mcsmatcher, 2);

    SharedDataStream sds(ds);

    sds << residxmcsmatcher.prematcher << double(residxmcsmatcher.t.to(second))
        << residxmcsmatcher.match_light << residxmcsmatcher.verbose
        << static_cast<const AtomMatcher&>(residxmcsmatcher);

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, ResIdxAtomMCSMatcher &residxmcsmatcher)
{
    VersionID v = readHeader(ds, r_residxmcsmatcher);

    if (v == 2)
    {
        double timeout;

        SharedDataStream sds(ds);

        sds >> residxmcsmatcher.prematcher >> timeout
            >> residxmcsmatcher.match_light >> residxmcsmatcher.verbose
            >> static_cast<AtomMatcher&>(residxmcsmatcher);

        residxmcsmatcher.t = timeout*second;
    }
    else if (v == 1)
    {
        double timeout;

        SharedDataStream sds(ds);

        sds >> residxmcsmatcher.prematcher >> timeout >> static_cast<AtomMatcher&>(residxmcsmatcher);

        residxmcsmatcher.t = timeout*second;
        residxmcsmatcher.match_light = false;
        residxmcsmatcher.verbose = false;
    }
    else
        throw version_error(v, "1,2", r_residxmcsmatcher, CODELOC);

    return ds;
}

/** Constructor */
ResIdxAtomMCSMatcher::ResIdxAtomMCSMatcher()
               : ConcreteProperty<ResIdxAtomMCSMatcher,AtomMatcher>(), t(1*second),
                 match_light(false),
                 verbose(false)
{}

/** Constructor */
ResIdxAtomMCSMatcher::ResIdxAtomMCSMatcher(bool verbose)
               : ConcreteProperty<ResIdxAtomMCSMatcher,AtomMatcher>(), t(1*second),
                 match_light(false),
                 verbose(verbose)
{}

/** Construct specifying the timeout for the MCS match */
ResIdxAtomMCSMatcher::ResIdxAtomMCSMatcher(const SireUnits::Dimension::Time &timeout,
                                           bool verbose)
               : ConcreteProperty<ResIdxAtomMCSMatcher,AtomMatcher>(), t(timeout),
                 match_light(false),
                 verbose(verbose)
{}

/** Construct specifying the prematcher for the MCS match */
ResIdxAtomMCSMatcher::ResIdxAtomMCSMatcher(const AtomMatcher &matcher,
                                           bool verbose)
               : ConcreteProperty<ResIdxAtomMCSMatcher,AtomMatcher>(),
                 prematcher(matcher), t(1*second),
                 match_light(false), verbose(verbose)
{}

/** Construct specifying the timeout and prematcher for the MCS match */
ResIdxAtomMCSMatcher::ResIdxAtomMCSMatcher(const AtomMatcher &matcher,
                                           const SireUnits::Dimension::Time &timeout,
                                           bool verbose)
               : ConcreteProperty<ResIdxAtomMCSMatcher,AtomMatcher>(),
                 prematcher(matcher), t(timeout),
                 match_light(false), verbose(verbose)
{}

/** Constructor, specifying whether or not to match light atoms */
ResIdxAtomMCSMatcher::ResIdxAtomMCSMatcher(bool match_light_atoms,
                                           bool verbose)
               : ConcreteProperty<ResIdxAtomMCSMatcher,AtomMatcher>(), t(1*second),
                 match_light(match_light_atoms),
                 verbose(verbose)
{}

/** Construct specifying the timeout for the MCS match, and specifying whether or not
    to match light atoms */
ResIdxAtomMCSMatcher::ResIdxAtomMCSMatcher(const SireUnits::Dimension::Time &timeout,
                                           bool match_light_atoms,
                                           bool verbose)
               : ConcreteProperty<ResIdxAtomMCSMatcher,AtomMatcher>(), t(timeout),
                 match_light(match_light_atoms),
                 verbose(verbose)
{}

/** Construct specifying the prematcher for the MCS match,
    and specifying whether or not to match light atoms
*/
ResIdxAtomMCSMatcher::ResIdxAtomMCSMatcher(const AtomMatcher &matcher,
                                           bool match_light_atoms,
                                           bool verbose)
               : ConcreteProperty<ResIdxAtomMCSMatcher,AtomMatcher>(),
                 prematcher(matcher), t(1*second),
                 match_light(match_light_atoms), verbose(verbose)
{}

/** Construct specifying the timeout and prematcher for the MCS match,
    and specifying whether or not to match light atoms
*/
ResIdxAtomMCSMatcher::ResIdxAtomMCSMatcher(const AtomMatcher &matcher,
                                           const SireUnits::Dimension::Time &timeout,
                                           bool match_light_atoms,
                                           bool verbose)
               : ConcreteProperty<ResIdxAtomMCSMatcher,AtomMatcher>(),
                 prematcher(matcher), t(timeout),
                 match_light(match_light_atoms), verbose(verbose)
{}

/** Copy constructor */
ResIdxAtomMCSMatcher::ResIdxAtomMCSMatcher(const ResIdxAtomMCSMatcher &other)
               : ConcreteProperty<ResIdxAtomMCSMatcher,AtomMatcher>(other),
                 prematcher(other.prematcher), t(other.t),
                 match_light(other.match_light), verbose(other.verbose)
{
}

/** Destructor */
ResIdxAtomMCSMatcher::~ResIdxAtomMCSMatcher()
{}

/** Copy assignment operator */
ResIdxAtomMCSMatcher& ResIdxAtomMCSMatcher::operator=(const ResIdxAtomMCSMatcher &other)
{
    if (this != &other)
    {
        t = other.t;
        prematcher = other.prematcher;
        match_light = other.match_light;
        verbose = other.verbose;
    }

    return *this;
}

/** Comparison operator */
bool ResIdxAtomMCSMatcher::operator==(const ResIdxAtomMCSMatcher &other) const
{
    return prematcher == other.prematcher and t == other.t and match_light == other.match_light and verbose == other.verbose;
}

/** Comparison operator */
bool ResIdxAtomMCSMatcher::operator!=(const ResIdxAtomMCSMatcher &other) const
{
    return not operator==(other);
}

QString ResIdxAtomMCSMatcher::toString() const
{
    if (prematcher.isNull() or prematcher.read().isNull())
    {
        return QObject::tr("ResIdxAtomMCSMatcher( timeout() = %1 s, matchingLightAtoms() = %2, "
                "isVerbose() = %3 )")
                .arg(t.to(second)).arg(match_light).arg(verbose);
    }
    else
    {
        return QObject::tr("ResIdxAtomMCSMatcher( preMatcher() = %1, timeout() = %2 s, "
                           "matchingLightAtoms() = %3, isVerbose() = %4 )")
                    .arg(prematcher.read().toString())
                    .arg(t.to(second))
                    .arg(match_light)
                    .arg(verbose);
    }
}

/** Return the prematcher (if any) that is used to pre-match atoms
    before the MCS match */
const AtomMatcher& ResIdxAtomMCSMatcher::preMatcher() const
{
    return prematcher.read();
}

/** Return the timeout before the MCS match is abandoned */
SireUnits::Dimension::Time ResIdxAtomMCSMatcher::timeout() const
{
    return t;
}

/** Return whether or not this will include light atoms (e.g. hydrogen)
    when searching for the maximum common substructure */
bool ResIdxAtomMCSMatcher::matchingLightAtoms() const
{
    return match_light;
}

/** Return whether or not this report progress to stdout. */
bool ResIdxAtomMCSMatcher::isVerbose() const
{
    return verbose;
}

/** Match the atoms in 'mol1' to the atoms in 'mol0' by MCS, searching each
    residue separately (by index) and combining the results. This returns the
    AtomIdxs of the atoms in 'mol1' that are in 'mol0', indexed by the AtomIdx
    of the atom in 'mol0'.

    This skips atoms in 'mol1' that are not in 'mol0'
*/
QHash<AtomIdx,AtomIdx> ResIdxAtomMCSMatcher::pvt_match(const MoleculeView &mol0,
                                                       const PropertyMap &map0,
                                                       const MoleculeView &mol1,
                                                       const PropertyMap &map1) const
{
    const AtomSelection sel0 = mol0.selection();
    const AtomSelection sel1 = mol1.selection();

    QHash<AtomIdx,AtomIdx> map;

    if (sel0.selectedAll() and sel1.selectedAll())
    {
        for (const auto &idx : mol0.data().info().getResidues())
        {
            // Perform an MCS match using the same residue in each molecule.
            QHash<AtomIdx,AtomIdx> local_map;
            if (prematcher.isNull() or prematcher.read().isNull())
            {
                local_map = Evaluator(mol0[idx]).findMCS(mol1[idx], t,
                    match_light, map0, map1, this->verbose);
            }
            else
            {
                local_map = Evaluator(mol0[idx]).findMCS(mol1[idx], prematcher.read(),
                    t, match_light, map0, map1, this->verbose);
            }

            // Update the atom index map.
            map.unite(local_map);
        }
    }
    else
    {
        throw SireError::unsupported(QObject::tr("ResIdxAtomMCSMatcher only works with "
            "full molecule selections."));
    }

    return map;
}

const char* ResIdxAtomMCSMatcher::typeName()
{
    return QMetaType::typeName( qMetaTypeId<ResIdxAtomNameMatcher>() );
}

/////////
///////// Implementation of ResIdxAtomCoordMatcher
/////////

static const RegisterMetaType<ResIdxAtomCoordMatcher> r_residxatomcoordmatcher;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                       const ResIdxAtomCoordMatcher &residxatomcoordmatcher)
{
    writeHeader(ds, r_residxatomcoordmatcher, 1);
    ds << static_cast<const ResIdxAtomCoordMatcher&>(residxatomcoordmatcher);

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, ResIdxAtomCoordMatcher &residxatomcoordmatcher)
{
    VersionID v = readHeader(ds, r_residxatomcoordmatcher);

    if (v == 1)
    {
        ds >> static_cast<ResIdxAtomCoordMatcher&>(residxatomcoordmatcher);
    }
    else
        throw version_error(v, "1", r_residxatomcoordmatcher, CODELOC);

    return ds;
}

/** Constructor */
ResIdxAtomCoordMatcher::ResIdxAtomCoordMatcher() : ConcreteProperty<ResIdxAtomCoordMatcher,AtomMatcher>()
{}

/** Copy constructor */
ResIdxAtomCoordMatcher::ResIdxAtomCoordMatcher(const ResIdxAtomCoordMatcher &other)
               : ConcreteProperty<ResIdxAtomCoordMatcher,AtomMatcher>(other)
{}

/** Destructor */
ResIdxAtomCoordMatcher::~ResIdxAtomCoordMatcher()
{}

/** Copy assignment operator */
ResIdxAtomCoordMatcher& ResIdxAtomCoordMatcher::operator=(const ResIdxAtomCoordMatcher &other)
{
    return *this;
}

/** Comparison operator */
bool ResIdxAtomCoordMatcher::operator==(const ResIdxAtomCoordMatcher &other) const
{
    return true;
}

/** Comparison operator */
bool ResIdxAtomCoordMatcher::operator!=(const ResIdxAtomCoordMatcher &other) const
{
    return false;
}

QString ResIdxAtomCoordMatcher::toString() const
{
    return QObject::tr("ResIdxAtomCoordMatcher()");
}

/** Match the atoms in 'mol1' to the atoms in 'mol0' by coordinates, searching each
    residue separately (by index) and combining the results. This returns the
    AtomIdxs of the atoms in 'mol1' that are in 'mol0', indexed by the AtomIdx
    of the atom in 'mol0'.

    This skips atoms in 'mol1' that are not in 'mol0'
*/
QHash<AtomIdx,AtomIdx> ResIdxAtomCoordMatcher::pvt_match(const MoleculeView &mol0,
                                                         const PropertyMap &map0,
                                                         const MoleculeView &mol1,
                                                         const PropertyMap &map1) const
{
    const AtomSelection sel0 = mol0.selection();
    const AtomSelection sel1 = mol1.selection();

    if (sel0.selectedAll() and sel1.selectedAll())
    {
        // A hash of mappings between atom indices in each molecule.
        QHash<AtomIdx,AtomIdx> matches;

        // Get the list of residue indices from the reference molecule.
        auto resIdxs = mol0.data().info().getResidues();

        // Vectors to store the centre of mass (CoM) of each molecule.
        Vector com0;
        Vector com1;

        // Work out the CoM of both molecules.

        // mol0
        for (int i=0; i<mol0.data().info().nAtoms(); ++i)
            com0 += mol0.atom(AtomIdx(i)).property<Vector>(map0["coordinates"]);
        com0 /= mol0.data().info().nAtoms();

        // mol1
        for (int i=0; i<mol1.data().info().nAtoms(); ++i)
            com1 += mol1.atom(AtomIdx(i)).property<Vector>(map1["coordinates"]);
        com1 /= mol1.data().info().nAtoms();

        // Loop over all of the residues.
        for (const auto &resIdx : resIdxs)
        {
            // Get a list of atoms for the residue for both molecules.
            auto atoms0 = mol0.data().info().getAtomsIn(resIdx);
            auto atoms1 = mol1.data().info().getAtomsIn(resIdx);

            // A set of matched atom indices.
            QSet<int> matched;

            // For each atom in atoms0, find the atom in atoms1 that is closest to it.
            // To account for possible coordinate frame translations, we shift the
            // coordinates of each atom by the CoM of its respective molecule.
            for (int i=0; i<atoms0.count(); ++i)
            {
                // Initialise the minimium difference to a large number.
                double min_diff = 1e6;

                // Initialise the match to an out-of-range number. This way we can tell
                // when no matches have been found.
                int match = -1;

                // Get the coordinates of atom0.
                auto coord0 = mol0.atom(atoms0[i])
                                  .property<Vector>(map0["coordinates"]);

                // Shift by the CoM.
                coord0 -= com0;

                // Loop over all of the atoms to match against.
                for (int j=0; j<atoms0.count(); ++j)
                {
                    // Get the coordinates of atom1.
                    auto coord1 = mol1.atom(atoms1[j])
                                      .property<Vector>(map0["coordinates"]);

                    // Shift by the CoM.
                    coord1 -= com1;

                    // Compute the separation between the atoms, accounting for
                    // the CoM. This avoids issues with coordinate frame translations.
                    double diff = qAbs((coord0 - coord1).magnitude());

                    // Is this the best match to date? If so, update the match and the
                    // minimum difference.
                    if (diff < min_diff)
                    {
                        // Has this atom already been matched?
                        if (not matched.contains(j))
                        {
                            min_diff = diff;
                            match = j;
                        }
                    }
                }

                // A match was found, store the best match and append to the list of
                // matched atoms.
                if (match != -1)
                {
                    matches.insert(atoms0[i], atoms1[match]);
                    matched.insert(match);
                }
            }
        }

        return matches;
    }
    else
    {
        throw SireError::unsupported(QObject::tr("ResIdxAtomCoordMatcher only works with "
            "full molecule selections."));
    }
}

const char* ResIdxAtomCoordMatcher::typeName()
{
    return QMetaType::typeName( qMetaTypeId<ResIdxAtomCoordMatcher>() );
}
