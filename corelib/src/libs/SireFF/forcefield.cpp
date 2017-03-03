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

#include <QMutex>

#include "forcefield.h"

#include "SireMol/molecule.h"
#include "SireMol/mover.hpp"
#include "SireMol/viewsofmol.h"
#include "SireMol/partialmolecule.h"
#include "SireMol/moleculegroup.h"
#include "SireMol/mgnum.h"

#include "SireError/errors.h"
#include "SireBase/errors.h"
#include "SireMol/errors.h"

#include "SireStream/datastream.h"

#include <QDebug>

using namespace SireFF;
using namespace SireBase;
using namespace SireStream;

///////////
/////////// Implementation of NullFF
///////////

static const RegisterMetaType<NullFF> r_nullff;

/** Serialise to a binary datastream */
QDataStream SIREFF_EXPORT &operator<<(QDataStream &ds, const NullFF &nullff)
{
    writeHeader(ds, r_nullff, 0);
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREFF_EXPORT &operator>>(QDataStream &ds, NullFF &nullff)
{
    VersionID v = readHeader(ds, r_nullff);
    
    if (v != 0)
        throw version_error(v, "0", r_nullff, CODELOC);
        
    return ds;
}

/** Private constructor */
NullFF::NullFF(bool) : ConcreteProperty<NullFF,FF>()
{}

const NullFF& FF::null()
{
    return *(create_shared_null<NullFF>());
}

/** Constructor */
NullFF::NullFF() : ConcreteProperty<NullFF,FF>()
{}

/** Copy constructor */
NullFF::NullFF(const NullFF &other) : ConcreteProperty<NullFF,FF>(other)
{}

/** Destructor */
NullFF::~NullFF()
{}

/** Copy assignment */
NullFF& NullFF::operator=(const NullFF &other)
{
    FF::operator=(other);
    return *this;
}

/** All null forcefields are equal */
bool NullFF::operator==(const NullFF &other) const
{
    return true;
}

/** All null forcefields are equal */
bool NullFF::operator!=(const NullFF &other) const
{
    return false;
}

QString NullFF::toString() const
{
    return "NullFF()";
}

bool NullFF::setProperty(const QString &name, const Property &value)
{
    throw SireBase::missing_property( QObject::tr(
        "The null forcefield has no properties!"), CODELOC );

    return false;
}

const Property& NullFF::property(const QString &name) const
{
    throw SireBase::missing_property( QObject::tr(
        "The null forcefield has no properties!"), CODELOC );

    return props.property(name);
}

bool NullFF::containsProperty(const QString &name) const
{
    return false;
}

const Properties& NullFF::properties() const
{
    return props;
}

void NullFF::mustNowRecalculateFromScratch()
{}

const FFComponent& NullFF::components() const
{
    throw SireError::invalid_state( QObject::tr(
        "The null forcefield has no components!"), CODELOC );
}

void NullFF::recalculateEnergy()
{}

static void throwNoAdd()
{
    throw SireError::invalid_state( QObject::tr(
        "You cannot add anything to the null forcefield."), CODELOC );
}

void NullFF::group_add(quint32 i, const MoleculeView &molview,
                       const PropertyMap &map)
{
    throwNoAdd();
}

void NullFF::group_add(quint32 i, const ViewsOfMol &molviews, 
                       const PropertyMap &map)
{
    throwNoAdd();
}

void NullFF::group_add(quint32 i, const Molecules &molecules, 
                       const PropertyMap &map)
{
    throwNoAdd();
}

void NullFF::group_add(quint32 i, const MoleculeGroup &molgroup, 
                       const PropertyMap &map)
{
    throwNoAdd();
}

bool NullFF::group_addIfUnique(quint32 i, const MoleculeView &molview, 
                               const PropertyMap &map)
{
    throwNoAdd();
    return false;
}

ViewsOfMol NullFF::group_addIfUnique(quint32 i, const ViewsOfMol &molviews, 
                                     const PropertyMap &map)
{
    throwNoAdd();
    return ViewsOfMol();
}

QList<ViewsOfMol> NullFF::group_addIfUnique(quint32 i, const Molecules &molecules, 
                                            const PropertyMap &map)
{
    throwNoAdd();
    throw QList<ViewsOfMol>();
}

QList<ViewsOfMol> NullFF::group_addIfUnique(quint32 i, const MoleculeGroup &molgroup, 
                                            const PropertyMap &map)
{
    throwNoAdd();
    return QList<ViewsOfMol>();
}

static void throwNoRemove()
{
    throw SireError::invalid_state( QObject::tr(
        "You cannot remove anything from the null forcefield."), CODELOC );
}

bool NullFF::group_remove(quint32 i, const MoleculeView &molview)
{
    throwNoRemove();
    return false;
}

ViewsOfMol NullFF::group_remove(quint32 i, const ViewsOfMol &molviews)
{
    throwNoRemove();
    return ViewsOfMol();
}

QList<ViewsOfMol> NullFF::group_remove(quint32 i, const Molecules &molecules)
{
    throwNoRemove();
    return QList<ViewsOfMol>();
}

QList<ViewsOfMol> NullFF::group_remove(quint32 i, const MoleculeGroup &molgroup)
{
    throwNoRemove();
    return QList<ViewsOfMol>();
}

bool NullFF::group_removeAll(quint32 i, const MoleculeView &molview)
{
    throwNoRemove();
    return false;
}

ViewsOfMol NullFF::group_removeAll(quint32 i, const ViewsOfMol &molviews)
{
    throwNoRemove();
    return ViewsOfMol();
}

QList<ViewsOfMol> NullFF::group_removeAll(quint32 i, const Molecules &molecules)
{
    throwNoRemove();
    return QList<ViewsOfMol>();
}

QList<ViewsOfMol> NullFF::group_removeAll(quint32 i, const MoleculeGroup &molgroup)
{
    throwNoRemove();
    return QList<ViewsOfMol>();
}

ViewsOfMol NullFF::group_remove(quint32 i, MolNum molnum)
{
    throwNoRemove();
    return ViewsOfMol();
}

QList<ViewsOfMol> NullFF::group_remove(quint32 i, const QSet<MolNum> &molnums)
{
    throwNoRemove();
    return QList<ViewsOfMol>();
}

void NullFF::group_removeAll(quint32 i)
{
    throwNoRemove();
}

bool NullFF::group_update(quint32 i, const MoleculeData &moldata, bool auto_commit)
{
    return false;
}

QList<Molecule> NullFF::group_update(quint32 i, const Molecules &molecules, bool auto_commit)
{
    return QList<Molecule>();
}

QList<Molecule> NullFF::group_update(quint32 i, const MoleculeGroup &molgroup, bool auto_commmit)
{
    return QList<Molecule>();
}

static void throwNoSetContents()
{
    throw SireError::invalid_state( QObject::tr(
        "You cannot set the contents of a null forcefield!"), CODELOC );
}

bool NullFF::group_setContents(quint32 i, const MoleculeView &molview, 
                               const PropertyMap &map)
{
    throwNoSetContents();
    return false;
}

bool NullFF::group_setContents(quint32 i, const ViewsOfMol &molviews, 
                               const PropertyMap &map)
{
    throwNoSetContents();
    return false;
}

bool NullFF::group_setContents(quint32 i, const Molecules &molecules, 
                               const PropertyMap &map)
{
    throwNoSetContents();
    return false;
}

bool NullFF::group_setContents(quint32 i, const MoleculeGroup &molgroup, 
                               const PropertyMap &map)
{
    throwNoSetContents();
    return false;
}

void NullFF::_pvt_updateName()
{
    throw SireError::invalid_state( QObject::tr(
        "You cannot change the name of a null forcefield."), CODELOC );
}

static void throwNoGroups()
{
    throw SireMol::missing_group( QObject::tr(
        "There are no molecule groups in the null forcefield!"), CODELOC );
}

const MoleculeGroup& NullFF::at(MGNum mgnum) const
{
    throwNoGroups();
    return this->getGroup(mgnum);
}

const MoleculeGroup& NullFF::getGroup(MGNum mgnum) const
{
    return NullFF::at(mgnum);
}

void NullFF::getGroups(const QList<MGNum> &mgnums,
                       QVarLengthArray<const MoleculeGroup*,10> &groups) const
{
    throwNoGroups();
}

QHash<MGNum,const MoleculeGroup*> NullFF::getGroups() const
{
    return QHash<MGNum,const MoleculeGroup*>();
}

void NullFF::group_setName(quint32 i, const QString &new_name)
{
    throwNoGroups();
}

void NullFF::reindex()
{
    MolGroupsBase::clearIndex();
}

bool NullFF::needsAccepting() const
{
    return false;
}

void NullFF::accept()
{}

const char* NullFF::typeName()
{
    return QMetaType::typeName( qMetaTypeId<NullFF>() );
}
