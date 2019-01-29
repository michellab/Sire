/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2009  Christopher Woods
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

#include "mgidsandmaps.h"

#include "SireID/idorset.hpp"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMol;
using namespace SireBase;
using namespace SireID;
using namespace SireStream;

using boost::tuple;

static const RegisterMetaType<MGIDsAndMaps> r_mgidsandmaps;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, 
                                       const MGIDsAndMaps &mgidsandmaps)
{
    writeHeader(ds, r_mgidsandmaps, 1);
    
    SharedDataStream sds(ds);
    
    sds << mgidsandmaps.mgids << mgidsandmaps.maps
        << static_cast<const Property&>(mgidsandmaps);

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds,
                                       MGIDsAndMaps &mgidsandmaps)
{
    VersionID v = readHeader(ds, r_mgidsandmaps);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> mgidsandmaps.mgids >> mgidsandmaps.maps
            >> static_cast<Property&>(mgidsandmaps);
    }
    else
        throw version_error(v, "1", r_mgidsandmaps, CODELOC);
        
    return ds;
}

/** Null constructor */
MGIDsAndMaps::MGIDsAndMaps()
             : ConcreteProperty<MGIDsAndMaps,Property>()
{}

/** Construct from the passed groups / maps */
MGIDsAndMaps::MGIDsAndMaps(const MoleculeGroup &mgroup)
             : ConcreteProperty<MGIDsAndMaps,Property>()
{
    mgids.append( mgroup.number() );
    maps.append( PropertyMap() );
    
    mgids.squeeze();
    maps.squeeze();
}

/** Construct from the passed groups / maps */
MGIDsAndMaps::MGIDsAndMaps(const MoleculeGroup &mgroup, const PropertyMap &map)
             : ConcreteProperty<MGIDsAndMaps,Property>()
{
    mgids.append( mgroup.number() );
    maps.append(map);
    
    mgids.squeeze();
    maps.squeeze();
}

/** Construct from the passed groups / maps */
MGIDsAndMaps::MGIDsAndMaps(const tuple<MolGroupPtr,PropertyMap> group_and_map)
             : ConcreteProperty<MGIDsAndMaps,Property>()
{
    mgids.append( group_and_map.get<0>().read().number() );
    maps.append( group_and_map.get<1>() );
    
    mgids.squeeze();
    maps.squeeze();
}

static IDOrSet<MGID> getIDSet(const QList<MGNum> &mgnums)
{
    if (not mgnums.isEmpty())
        return IDOrSet<MGID>(mgnums);
    else
        return IDOrSet<MGID>();
}

static IDOrSet<MGID> getIDSet(const QList<MGIdentifier> &mgids)
{
    if (not mgids.isEmpty())
        return IDOrSet<MGID>(mgids);
    else
        return IDOrSet<MGID>();
}

static IDOrSet<MGID> getIDSet(const QList<MolGroupPtr> &mgroups)
{
    QList<MGNum> mgnums;
    
    for (QList<MolGroupPtr>::const_iterator it = mgroups.constBegin();
         it != mgroups.constEnd();
         ++it)
    {
        mgnums.append( it->read().number() );
    }
    
    return ::getIDSet(mgnums);
}

/** Construct from the passed groups / maps */
MGIDsAndMaps::MGIDsAndMaps(const QList<MolGroupPtr> &mgroups)
             : ConcreteProperty<MGIDsAndMaps,Property>()
{
    if (not mgroups.isEmpty())
    {
        mgids.append( ::getIDSet(mgroups) );
        maps.append( PropertyMap() );
        
        mgids.squeeze();
        maps.squeeze();
    }
}

/** Construct from the passed groups / maps */
MGIDsAndMaps::MGIDsAndMaps(const QList< tuple<MolGroupPtr,PropertyMap> > &groups_and_maps)
             : ConcreteProperty<MGIDsAndMaps,Property>()
{
    if (not groups_and_maps.isEmpty())
    {
        for (QList< tuple<MolGroupPtr,PropertyMap> >::const_iterator 
                                            it = groups_and_maps.constBegin();
             it != groups_and_maps.constEnd();
             ++it)
        {
            mgids.append( it->get<0>().read().number() );
            maps.append( it->get<1>() );
        }
        
        mgids.squeeze();
        maps.squeeze();
    }
}

/** Construct from the passed groups / maps */
MGIDsAndMaps::MGIDsAndMaps(const QList<MolGroupPtr> &mgroups, const PropertyMap &map)
             : ConcreteProperty<MGIDsAndMaps,Property>()
{
    if (not mgroups.isEmpty())
    {
        mgids.append( ::getIDSet(mgroups) );
        maps.append( map );
        
        mgids.squeeze();
        maps.squeeze();
    }
}

/** Construct from the passed groups / maps */
MGIDsAndMaps::MGIDsAndMaps(const tuple<QList<MolGroupPtr>,PropertyMap> &groups_and_maps)
             : ConcreteProperty<MGIDsAndMaps,Property>()
{
    if (not groups_and_maps.get<0>().isEmpty())
    {
        mgids.append( ::getIDSet(groups_and_maps.get<0>()) );
        maps.append( groups_and_maps.get<1>() );
        
        mgids.squeeze();
        maps.squeeze();
    }
}

/** Construct from the passed groups / maps */
MGIDsAndMaps::MGIDsAndMaps(const MGID &mgid)
             : ConcreteProperty<MGIDsAndMaps,Property>()
{
    mgids.append(mgid);
    maps.append(PropertyMap());
    
    mgids.squeeze();
    maps.squeeze();
}

/** Construct from the passed groups / maps */
MGIDsAndMaps::MGIDsAndMaps(const MGID &mgid, const PropertyMap &map)
             : ConcreteProperty<MGIDsAndMaps,Property>()
{
    mgids.append(mgid);
    maps.append(map);
    
    mgids.squeeze();
    maps.squeeze();
}

/** Construct from the passed groups / maps */
MGIDsAndMaps::MGIDsAndMaps(const tuple<MGIdentifier,PropertyMap> mgid_and_map)
             : ConcreteProperty<MGIDsAndMaps,Property>()
{
    mgids.append( mgid_and_map.get<0>() );
    maps.append( mgid_and_map.get<1>() );
    
    mgids.squeeze();
    maps.squeeze();
}

/** Construct from the passed groups / maps */
MGIDsAndMaps::MGIDsAndMaps(const QList<MGIdentifier> &mgs)
             : ConcreteProperty<MGIDsAndMaps,Property>()
{
    if (not mgs.isEmpty())
    {
        mgids.append( ::getIDSet(mgs) );
        maps.append( PropertyMap() );
        
        mgids.squeeze();
        maps.squeeze();
    }
}

/** Construct from the passed groups / maps */
MGIDsAndMaps::MGIDsAndMaps(const QList< tuple<MGIdentifier,PropertyMap> > &mgids_and_maps)
             : ConcreteProperty<MGIDsAndMaps,Property>()
{
    if (not mgids_and_maps.isEmpty())
    {
        for (QList< tuple<MGIdentifier,PropertyMap> >::const_iterator
                                                    it = mgids_and_maps.constBegin();
             it != mgids_and_maps.constEnd();
             ++it)
        {
            mgids.append( it->get<0>() );
            maps.append( it->get<1>() );
        }
        
        mgids.squeeze();
        maps.squeeze();
    }
}

/** Construct from the passed groups / maps */
MGIDsAndMaps::MGIDsAndMaps(const QList<MGIdentifier> &mgs, const PropertyMap &map)
             : ConcreteProperty<MGIDsAndMaps,Property>()
{
    if (not mgs.isEmpty())
    {
        mgids.append( ::getIDSet(mgs) );
        maps.append( map );
        
        mgids.squeeze();
        maps.squeeze();
    }
}

/** Construct from the passed groups / maps */
MGIDsAndMaps::MGIDsAndMaps(const tuple<QList<MGIdentifier>,PropertyMap> &mgids_and_maps)
             : ConcreteProperty<MGIDsAndMaps,Property>()
{
    if (not mgids_and_maps.get<0>().isEmpty())
    {
        mgids.append( ::getIDSet(mgids_and_maps.get<0>()) );
        maps.append( mgids_and_maps.get<1>() );
        
        mgids.squeeze();
        maps.squeeze();
    }
}

/** Construct from the passed groups / maps */
MGIDsAndMaps::MGIDsAndMaps(const QList<MGIDsAndMaps> &mgids_and_maps)
             : ConcreteProperty<MGIDsAndMaps,Property>()
{
    if (not mgids_and_maps.isEmpty())
    {
        for (QList<MGIDsAndMaps>::const_iterator it = mgids_and_maps.constBegin();
             it != mgids_and_maps.constEnd();
             ++it)
        {
            mgids += it->mgids;
            maps += it->maps;
        }
        
        mgids.squeeze();
        maps.squeeze();
    }
}

/** Copy constructor */
MGIDsAndMaps::MGIDsAndMaps(const MGIDsAndMaps &other)
             : ConcreteProperty<MGIDsAndMaps,Property>(other),
               mgids(other.mgids), maps(other.maps)
{}

/** Destructor */ 
MGIDsAndMaps::~MGIDsAndMaps()
{}

/** Copy assignment operator */
MGIDsAndMaps& MGIDsAndMaps::operator=(const MGIDsAndMaps &other)
{
    if (this != &other)
    {
        mgids = other.mgids;
        maps = other.maps;
        Property::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool MGIDsAndMaps::operator==(const MGIDsAndMaps &other) const
{
    return this == &other or
           (maps == other.maps and mgids == other.mgids);
}

/** Comparison operator */
bool MGIDsAndMaps::operator!=(const MGIDsAndMaps &other) const
{
    return not MGIDsAndMaps::operator==(other);
}

const char* MGIDsAndMaps::typeName()
{
    return QMetaType::typeName( qMetaTypeId<MGIDsAndMaps>() );
}

/** Return whether or not this is an empty group of groups */
bool MGIDsAndMaps::isEmpty() const
{
    return mgids.isEmpty();
}

/** Return the number of identifier groups in this set */
int MGIDsAndMaps::count() const
{
    return mgids.count();
}

/** Return a string representation of this set of IDs */
QString MGIDsAndMaps::toString() const
{
    QStringList ids;

    for (int i=0; i<mgids.count(); ++i)
    {
        if (maps[i].isDefault())
        {
            ids.append( mgids[i].toString() );
        }
        else 
        {
            ids.append( QString("(%1 : %2)")
                            .arg(mgids[i].toString(), maps[i].toString()) );
        }

    }
    
    return QString("MGIDsAndMaps[ %1 ]").arg(ids.join(", "));
}

/** Return the array of all of the IDs for this group of groups */
const QVector<MGIdentifier>& MGIDsAndMaps::mgIDs() const
{
    return mgids;
}

/** Return the array of the property maps to use for each group
    of IDs in this set (these are in the same order as the groups
    returned in MGIDsAndMaps::mgIDs()) */
const QVector<PropertyMap>& MGIDsAndMaps::propertyMaps() const
{
    return maps;
}
