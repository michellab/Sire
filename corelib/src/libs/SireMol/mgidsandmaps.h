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

#ifndef SIREMOL_MGIDS_AND_MAP_H
#define SIREMOL_MGIDS_AND_MAP_H

#include "mgidentifier.h"
#include "moleculegroup.h"

#include "SireBase/property.h"
#include "SireBase/propertymap.h"

#include <QVector>
#include <QList>

#include <boost/tuple/tuple.hpp>

SIRE_BEGIN_HEADER

namespace SireMol
{
class MGIDsAndMaps;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::MGIDsAndMaps&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::MGIDsAndMaps&);

namespace SireMol
{

/** This class holds a set of molecule group IDs, together
    with the property maps that should be used with the
    associated molecule groups. This provides a store for
    information used to identify groups of molecule groups,
    together with the properties needed to manipulate those
    groups
    
    @author Christopher Woods
*/
class SIREMOL_EXPORT MGIDsAndMaps 
        : public SireBase::ConcreteProperty<MGIDsAndMaps,SireBase::Property>
{

friend QDataStream& ::operator<<(QDataStream&, const MGIDsAndMaps&);
friend QDataStream& ::operator>>(QDataStream&, MGIDsAndMaps&);

public:
    MGIDsAndMaps();

    MGIDsAndMaps(const MoleculeGroup &mgroup);
    MGIDsAndMaps(const MoleculeGroup &mgroup, const PropertyMap &map);

    MGIDsAndMaps(const boost::tuple<MolGroupPtr,PropertyMap> group_and_map);

    MGIDsAndMaps(const QList<MolGroupPtr> &mgroups);
    MGIDsAndMaps(const QList< boost::tuple<MolGroupPtr,PropertyMap> > &groups_and_maps);

    MGIDsAndMaps(const QList<MolGroupPtr> &mgroups, const PropertyMap &map);
    MGIDsAndMaps(const boost::tuple<QList<MolGroupPtr>,PropertyMap> &groups_and_maps);
    
    MGIDsAndMaps(const MGID &mgid);
    MGIDsAndMaps(const MGID &mgid, const PropertyMap &map);
    
    MGIDsAndMaps(const boost::tuple<MGIdentifier,PropertyMap> mgid_and_map);
    
    MGIDsAndMaps(const QList<MGIdentifier> &mgids);
    MGIDsAndMaps(const QList< boost::tuple<MGIdentifier,PropertyMap> > &mgids_and_maps);

    MGIDsAndMaps(const QList<MGIdentifier> &mgids, const PropertyMap &map);
    MGIDsAndMaps(const boost::tuple<QList<MGIdentifier>,PropertyMap> &mgids_and_maps);
    
    MGIDsAndMaps(const QList<MGIDsAndMaps> &mgids_and_maps);

    MGIDsAndMaps(const MGIDsAndMaps &other);

    ~MGIDsAndMaps();
    
    MGIDsAndMaps& operator=(const MGIDsAndMaps &other);
    
    bool operator==(const MGIDsAndMaps &other) const;
    bool operator!=(const MGIDsAndMaps &other) const;
    
    static const char* typeName();
    
    bool isEmpty() const;
    
    int count() const;

    QString toString() const;
    
    const QVector<MGIdentifier>& mgIDs() const;
    const QVector<PropertyMap>& propertyMaps() const;
    
private:
    QVector<MGIdentifier> mgids;
    QVector<PropertyMap> maps;
};    

}

Q_DECLARE_METATYPE( SireMol::MGIDsAndMaps )

SIRE_EXPOSE_CLASS( SireMol::MGIDsAndMaps )

SIRE_END_HEADER

#endif
