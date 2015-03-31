/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2008  Christopher Woods
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

#include "atommatcher.h"
#include "atommatchers.h"
#include "atomidx.h"
#include "atomname.h"
#include "atomselection.h"
#include "atomidentifier.h"
#include "evaluator.h"
#include "moleculeinfodata.h"
#include "moleculeview.h"

#include "tostring.h"

#include "SireUnits/units.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"

using namespace SireMol;
using namespace SireUnits;
using namespace SireBase;
using namespace SireStream;

/////////
///////// Implmentation of AtomMatcher
/////////

AtomMultiMatcher *null_matcher = 0;

const AtomMultiMatcher& AtomMatcher::null()
{
    if (not null_matcher)
        null_matcher = new AtomMultiMatcher();
    
    return *null_matcher;
}

static const RegisterMetaType<AtomMatcher> r_atommatcher( MAGIC_ONLY,
                                                          "SireMol::AtomMatcher" );

/** Serialise to a binary datastream */
QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds, const AtomMatcher &matcher)
{
    writeHeader(ds, r_atommatcher, 1);
    
    ds << static_cast<const Property&>(matcher);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds, AtomMatcher &matcher)
{
    VersionID v = readHeader(ds, r_atommatcher);
    
    if (v == 1)
    {
        ds >> static_cast<Property&>(matcher);
    }
    else
        throw version_error(v, "1", r_atommatcher, CODELOC);
        
    return ds;
}

/** Constructor */
AtomMatcher::AtomMatcher() : Property()
{}

/** Copy constructor */
AtomMatcher::AtomMatcher(const AtomMatcher &other) : Property(other)
{}

/** Destructor */
AtomMatcher::~AtomMatcher()
{}

/** Return whether or not this matcher is null (cannot be used for matching) */
bool AtomMatcher::isNull() const
{
    return false;
}

/** Return the matcher that matches using this matcher, and then 'other' (in that order) */
AtomMultiMatcher AtomMatcher::operator+(const AtomMatcher &other) const
{
    return AtomMultiMatcher(*this, other);
}
    
/** Return the matcher that matches using this matcher, and then 'other' (in that order) */
AtomMultiMatcher AtomMatcher::add(const AtomMatcher &other) const
{
    return AtomMultiMatcher(*this, other);
}

/** Return whether or not this match changes the order of number of atoms */
bool AtomMatcher::pvt_changesOrder(const MoleculeInfoData &mol0,
                                   const MoleculeInfoData &mol1) const
{
    if (mol0.nAtoms() != mol1.nAtoms())
        return true;
    
    QHash<AtomIdx,AtomIdx> map = this->match(mol0, mol1);
    
    if (map.count() != mol0.nAtoms())
        return true;
    
    for (QHash<AtomIdx,AtomIdx>::const_iterator it = map.constBegin();
         it != map.constEnd();
         ++it)
    {
        if (it.key() != it.value())
            return true;
    }
    
    return false;
}

/** Return whether or not this match changes the order or number of viewed atoms */
bool AtomMatcher::pvt_changesOrder(const MoleculeView &molview0,
                                   const PropertyMap &map0,
                                   const MoleculeView &molview1,
                                   const PropertyMap &map1) const
{
    const int nats0 = molview0.selection().nSelectedAtoms();
    const int nats1 = molview1.selection().nSelectedAtoms();
    
    if (nats0 != nats1)
        return true;

    QHash<AtomIdx,AtomIdx> map = this->match(molview0,map0,molview1,map1);
    
    if (map.count() != nats0)
        return true;
    
    for (QHash<AtomIdx,AtomIdx>::const_iterator it = map.constBegin();
         it != map.constEnd();
         ++it)
    {
        if (it.key() != it.value())
            return true;
    }
    
    return false;
}

bool AtomMatcher::pvt_changesOrder(const MoleculeView &molview0,
                                   const MoleculeView &molview1) const
{
    return this->changesOrder(molview0, PropertyMap(), molview1, PropertyMap());
}

bool AtomMatcher::pvt_changesOrder(const MoleculeView &molview0,
                                   const MoleculeView &molview1,
                                   const PropertyMap &map) const
{
    return this->changesOrder(molview0, map, molview1, map);
}

QHash<AtomIdx,AtomIdx> AtomMatcher::pvt_match(const MoleculeView &molview0,
                                              const MoleculeView &molview1) const
{
    return this->match(molview0, PropertyMap(), molview1, PropertyMap());
}

QHash<AtomIdx,AtomIdx> AtomMatcher::pvt_match(const MoleculeView &molview0,
                                              const MoleculeView &molview1,
                                              const PropertyMap &map) const
{
    return this->match(molview0, map, molview1, map);
}

/** Match atoms based only on the data in the MoleculeInfoData of the molecules. */
QHash<AtomIdx,AtomIdx> AtomMatcher::pvt_match(const MoleculeInfoData &mol0,
                                              const MoleculeInfoData &mol1) const
{
    throw SireError::unsupported( QObject::tr(
                "The AtomMatcher \"%1\" does not support matching using "
                "MoleculeInfoData objects only.")
                    .arg(this->toString()), CODELOC );

    return QHash<AtomIdx,AtomIdx>();
}

/** Return whether or not this match changes the order of number of atoms */
bool AtomMatcher::changesOrder(const MoleculeInfoData &mol0,
                               const MoleculeInfoData &mol1) const
{
    return this->pvt_changesOrder(mol0, mol1);
}

/** Return whether or not this match changes the order or number of viewed atoms */
bool AtomMatcher::changesOrder(const MoleculeView &molview0,
                               const PropertyMap &map0,
                               const MoleculeView &molview1,
                               const PropertyMap &map1) const
{
    return this->pvt_changesOrder(molview0,map0, molview1,map1);
}

bool AtomMatcher::changesOrder(const MoleculeView &molview0,
                               const MoleculeView &molview1) const
{
    return this->changesOrder(molview0, PropertyMap(), molview1, PropertyMap());
}

bool AtomMatcher::changesOrder(const MoleculeView &molview0,
                               const MoleculeView &molview1,
                               const PropertyMap &map) const
{
    return this->changesOrder(molview0, map, molview1, map);
}

QHash<AtomIdx,AtomIdx> AtomMatcher::match(const MoleculeView &molview0,
                                          const MoleculeView &molview1) const
{
    return this->pvt_match(molview0,molview1);
}

QHash<AtomIdx,AtomIdx> AtomMatcher::match(const MoleculeView &molview0,
                                          const MoleculeView &molview1,
                                          const PropertyMap &map) const
{
    return this->pvt_match(molview0,molview1,map);
}

QHash<AtomIdx,AtomIdx> AtomMatcher::match(const MoleculeView &molview0,
                                          const PropertyMap &map0,
                                          const MoleculeView &molview1,
                                          const PropertyMap &map1) const
{
    return this->pvt_match(molview0,map0,molview1,map1);
}

/** Match atoms based only on the data in the MoleculeInfoData of the molecules. */
QHash<AtomIdx,AtomIdx> AtomMatcher::match(const MoleculeInfoData &mol0,
                                          const MoleculeInfoData &mol1) const
{
    return this->pvt_match(mol0,mol1);
}

/////////
///////// Implmentation of AtomResultMatcher
/////////

static const RegisterMetaType<AtomResultMatcher> r_resmatcher;

/** Serialise to a binary datastream */
QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds, const AtomResultMatcher &resmatcher)
{
    writeHeader(ds, r_resmatcher, 1);
    SharedDataStream sds(ds);
    sds << resmatcher.m << static_cast<const AtomMatcher&>(resmatcher);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds, AtomResultMatcher &resmatcher)
{
    VersionID v = readHeader(ds, r_resmatcher);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> resmatcher.m >> static_cast<AtomMatcher&>(resmatcher);
    }
    else
        throw version_error(v, "1", r_resmatcher, CODELOC);

    return ds;
}

/** Constructor */
AtomResultMatcher::AtomResultMatcher() : ConcreteProperty<AtomResultMatcher,AtomMatcher>()
{}

/** Constructor */
AtomResultMatcher::AtomResultMatcher(const QHash<AtomIdx,AtomIdx> &results, bool invert)
                  : ConcreteProperty<AtomResultMatcher,AtomMatcher>(), m(results)
{
    if (invert and not results.isEmpty())
    {
        //invert the map (this allows reverse lookups)
        m.clear();
        m.reserve(results.count());
        
        for (QHash<AtomIdx,AtomIdx>::const_iterator it = results.constBegin();
             it != results.constEnd();
             ++it)
        {
            m.insert( it.value(), it.key() );
        }
    }
}

/** Copy constructor */
AtomResultMatcher::AtomResultMatcher(const AtomResultMatcher &other)
                  : ConcreteProperty<AtomResultMatcher,AtomMatcher>(other),
                    m(other.m)
{}

/** Destructor */
AtomResultMatcher::~AtomResultMatcher()
{}

/** Copy assignment operator */
AtomResultMatcher& AtomResultMatcher::operator=(const AtomResultMatcher &other)
{
    m = other.m;
    return *this;
}

/** Comparison operator */
bool AtomResultMatcher::operator==(const AtomResultMatcher &other) const
{
    return m == other.m;
}

/** Comparison operator */
bool AtomResultMatcher::operator!=(const AtomResultMatcher &other) const
{
    return not operator==(other);
}

bool AtomResultMatcher::isNull() const
{
    return m.isEmpty();
}

QString AtomResultMatcher::toString() const
{
    if (isNull())
        return QObject::tr("AtomResultMatcher::null");
    else
        return QObject::tr("AtomResultMatcher( %1 )").arg(Sire::toString(m));
}

/** Match the atoms in 'mol1' to the atoms in 'mol0' - this
    returns the AtomIdxs of the atoms in 'mol1' that are in
    'mol0', indexed by the AtomIdx of the atom in 'mol0'.
    
     This skips atoms in 'mol1' that are not in 'mol0'
*/
QHash<AtomIdx,AtomIdx> AtomResultMatcher::pvt_match(const MoleculeView &mol0,
                                                    const PropertyMap &map0,
                                                    const MoleculeView &mol1,
                                                    const PropertyMap &map1) const
{
    const AtomSelection sel0 = mol0.selection();
    const AtomSelection sel1 = mol1.selection();
    
    QHash<AtomIdx,AtomIdx> map;
    
    const int nats0 = mol0.data().info().nAtoms();
    const int nats1 = mol1.data().info().nAtoms();
    
    for (QHash<AtomIdx,AtomIdx>::const_iterator it = m.constBegin();
         it != m.constEnd();
         ++it)
    {
        if (it.key().value() >= 0 and it.key().value() < nats0)
        {
            if (it.value().value() >= 0 and it.value().value() < nats1)
            {
                if (sel0.selected(it.key()) and sel1.selected(it.value()))
                    map.insert(it.key(), it.value());
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
QHash<AtomIdx,AtomIdx> AtomResultMatcher::pvt_match(const MoleculeInfoData &mol0,
                                                    const MoleculeInfoData &mol1) const
{
    QHash<AtomIdx,AtomIdx> map;
    
    const int nats0 = mol0.nAtoms();
    const int nats1 = mol1.nAtoms();
    
    for (QHash<AtomIdx,AtomIdx>::const_iterator it = m.constBegin();
         it != m.constEnd();
         ++it)
    {
        if (it.key().value() >= 0 and it.key().value() < nats0)
        {
            if (it.value().value() >= 0 and it.value().value() < nats1)
            {
                map.insert(it.key(), it.value());
            }
        }
    }
    
    return map;
}

const char* AtomResultMatcher::typeName()
{
    return QMetaType::typeName( qMetaTypeId<AtomResultMatcher>() );
}

/////////
///////// Implmentation of AtomMatchInverter
/////////

static const RegisterMetaType<AtomMatchInverter> r_inverter;

/** Serialise to a binary datastream */
QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds, const AtomMatchInverter &inverter)
{
    writeHeader(ds, r_inverter, 1);
    SharedDataStream sds(ds);
    sds << inverter.m << static_cast<const AtomMatcher&>(inverter);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds, AtomMatchInverter &inverter)
{
    VersionID v = readHeader(ds, r_inverter);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> inverter.m >> static_cast<AtomMatcher&>(inverter);
    }
    else
        throw version_error(v, "1", r_inverter, CODELOC);

    return ds;
}

/** Constructor */
AtomMatchInverter::AtomMatchInverter() : ConcreteProperty<AtomMatchInverter,AtomMatcher>()
{}

/** Constructor */
AtomMatchInverter::AtomMatchInverter(const AtomMatcher &matcher)
                  : ConcreteProperty<AtomMatchInverter,AtomMatcher>()
{
    if (not matcher.isNull())
        m = matcher;
}

/** Copy constructor */
AtomMatchInverter::AtomMatchInverter(const AtomMatchInverter &other)
                  : ConcreteProperty<AtomMatchInverter,AtomMatcher>(other),
                    m(other.m)
{}

/** Destructor */
AtomMatchInverter::~AtomMatchInverter()
{}

/** Copy assignment operator */
AtomMatchInverter& AtomMatchInverter::operator=(const AtomMatchInverter &other)
{
    m = other.m;
    return *this;
}

/** Comparison operator */
bool AtomMatchInverter::operator==(const AtomMatchInverter &other) const
{
    return m == other.m;
}

/** Comparison operator */
bool AtomMatchInverter::operator!=(const AtomMatchInverter &other) const
{
    return not operator==(other);
}

bool AtomMatchInverter::isNull() const
{
    return m.constData() == 0 or m.read().isNull();
}

QString AtomMatchInverter::toString() const
{
    if (isNull())
        return QObject::tr("AtomMatchInverter::null");
    else
        return QObject::tr("AtomMatchInverter{ %1 }").arg(m.read().toString());
}

/** Match the atoms in 'mol1' to the atoms in 'mol0' - this
    returns the AtomIdxs of the atoms in 'mol1' that are in
    'mol0', indexed by the AtomIdx of the atom in 'mol0'.
    
     This skips atoms in 'mol1' that are not in 'mol0'
*/
QHash<AtomIdx,AtomIdx> AtomMatchInverter::pvt_match(const MoleculeView &mol0,
                                                    const PropertyMap &map0,
                                                    const MoleculeView &mol1,
                                                    const PropertyMap &map1) const
{
    if (isNull())
        return QHash<AtomIdx,AtomIdx>();
    
    //apply the match backwards, and then invert the result
    QHash<AtomIdx,AtomIdx> map = m.read().match(mol1,map1,mol0,map0);
    
    //invert the match
    if (not map.isEmpty())
    {
        QHash<AtomIdx,AtomIdx> invmap;
        invmap.reserve(map.count());
        
        for (QHash<AtomIdx,AtomIdx>::const_iterator it = map.constBegin();
             it != map.constEnd();
             ++it)
        {
            invmap.insert( it.value(), it.key() );
        }
        
        return invmap;
    }
    else
        return map;
}
/** Match the atoms in 'mol1' to the atoms in 'mol0' - this
    returns the AtomIdxs of the atoms in 'mol1' that are in
    'mol0', indexed by the AtomIdx of the atom in 'mol0'.
    
     This skips atoms in 'mol1' that are not in 'mol0'
*/
QHash<AtomIdx,AtomIdx> AtomMatchInverter::pvt_match(const MoleculeInfoData &mol0,
                                                    const MoleculeInfoData &mol1) const
{
    if (isNull())
        return QHash<AtomIdx,AtomIdx>();
    
    //apply the match backwards, and then invert the result
    QHash<AtomIdx,AtomIdx> map = m.read().match(mol1,mol0);
    
    //invert the match
    if (not map.isEmpty())
    {
        QHash<AtomIdx,AtomIdx> invmap;
        invmap.reserve(map.count());
        
        for (QHash<AtomIdx,AtomIdx>::const_iterator it = map.constBegin();
             it != map.constEnd();
             ++it)
        {
            invmap.insert( it.value(), it.key() );
        }
        
        return invmap;
    }
    else
        return map;
}

const char* AtomMatchInverter::typeName()
{
    return QMetaType::typeName( qMetaTypeId<AtomMatchInverter>() );
}
