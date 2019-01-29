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

#ifndef SIREVOL_VOLUMECHANGER_H
#define SIREVOL_VOLUMECHANGER_H

#include "SireBase/property.h"
#include "SireBase/propertymap.h"

#include "SireMol/mgidentifier.h"

#include "SireMaths/rangenerator.h"

#include "SireFF/point.h"

#include "SireUnits/dimensions.h"

SIRE_BEGIN_HEADER

namespace SireMove
{
class VolumeChanger;
class NullVolumeChanger;

class ScaleVolumeFromCenter;
}

SIREMOVE_EXPORT QDataStream& operator<<(QDataStream&, const SireMove::VolumeChanger&);
SIREMOVE_EXPORT QDataStream& operator>>(QDataStream&, SireMove::VolumeChanger&);

SIREMOVE_EXPORT QDataStream& operator<<(QDataStream&, const SireMove::NullVolumeChanger&);
SIREMOVE_EXPORT QDataStream& operator>>(QDataStream&, SireMove::NullVolumeChanger&);

SIREMOVE_EXPORT QDataStream& operator<<(QDataStream&, const SireMove::ScaleVolumeFromCenter&);
SIREMOVE_EXPORT QDataStream& operator>>(QDataStream&, SireMove::ScaleVolumeFromCenter&);

namespace SireSystem
{
class System;
}

namespace SireMol
{
class MoleculeGroup;
}

namespace SireMove
{

using SireMaths::RanGenerator;

using SireMol::MoleculeGroup;
using SireMol::MGID;

using SireFF::Point;
using SireFF::PointRef;

using SireBase::PropertyMap;

using SireSystem::System;

/** This is the base class of all volume changing function classes.
    These classes are used to change the volume of a system
    i.e. during a VolumeMove
    
    @author Christopher Woods
*/
class SIREMOVE_EXPORT VolumeChanger : public SireBase::Property
{

friend QDataStream& ::operator<<(QDataStream&, const VolumeChanger&);
friend QDataStream& ::operator>>(QDataStream&, VolumeChanger&);

public:
    VolumeChanger();
    VolumeChanger(const MGID &mgid);
    VolumeChanger(const MoleculeGroup &molgroup);
    
    VolumeChanger(const VolumeChanger &other);
    
    virtual ~VolumeChanger();
    
    VolumeChanger* clone() const=0;
    
    static const char* typeName()
    {
        return "SireMove::VolumeChanger";
    }

    void setGenerator(const RanGenerator &generator);
    const RanGenerator& generator() const;
    
    const MGID& groupID() const;
    
    void setGroup(const MGID &mgid);
    void setGroup(const MoleculeGroup &molgroup);
    
    virtual int setVolume(System &system,
                          const SireUnits::Dimension::Volume &volume,
                          const PropertyMap &map = PropertyMap()) const=0;
    
    virtual int changeVolume(System &system,
                             const SireUnits::Dimension::Volume &delta,
                             const PropertyMap &map = PropertyMap()) const;
    
    virtual int randomChangeVolume(System &system, 
                                   const SireUnits::Dimension::Volume &maxvolchange,
                                   double &new_bias, double &old_bias,
                                   const PropertyMap &map = PropertyMap()) const;
    
    static const NullVolumeChanger& null();
    
protected:
    VolumeChanger& operator=(const VolumeChanger &other);
    
    bool operator==(const VolumeChanger &other) const;
    bool operator!=(const VolumeChanger &other) const;
    
private:
    /** The random number generator used by this volume changer */
    RanGenerator rangen;

    /** The ID of the molecule group(s) that will be 
        moved by this volume changer */
    SireMol::MGIdentifier mgid;
};

/** This is a null volume changer that does nothing */
class SIREMOVE_EXPORT NullVolumeChanger
          : public SireBase::ConcreteProperty<NullVolumeChanger,VolumeChanger>
{

friend QDataStream& ::operator<<(QDataStream&, const NullVolumeChanger&);
friend QDataStream& ::operator>>(QDataStream&, NullVolumeChanger&);

public:
    NullVolumeChanger();
    
    NullVolumeChanger(const NullVolumeChanger &other);
    
    ~NullVolumeChanger();
    
    NullVolumeChanger& operator=(const NullVolumeChanger &other);
    
    bool operator==(const NullVolumeChanger &other) const;
    bool operator!=(const NullVolumeChanger &other) const;
    
    static const char* typeName();
    
    int setVolume(System &system,
                  const SireUnits::Dimension::Volume &volume,
                  const PropertyMap &map = PropertyMap()) const;
};

/** This is a volume changer that works by scaling the molecules 
    from a user-supplied center point, scaling those molecules closest
    to the center point the least, and those furthest from the 
    point the most
    
    @author Christopher Woods
*/
class SIREMOVE_EXPORT ScaleVolumeFromCenter
            : public SireBase::ConcreteProperty<ScaleVolumeFromCenter,VolumeChanger>
{

friend QDataStream& ::operator<<(QDataStream&, const ScaleVolumeFromCenter&);
friend QDataStream& ::operator>>(QDataStream&, ScaleVolumeFromCenter&);

public:
    ScaleVolumeFromCenter();
    
    ScaleVolumeFromCenter(const MGID &mgid);
    ScaleVolumeFromCenter(const MoleculeGroup &molgroup);
    
    ScaleVolumeFromCenter(const MGID &mgid, const PointRef &point);
    ScaleVolumeFromCenter(const MoleculeGroup &molgroup, const PointRef &point);
    
    ScaleVolumeFromCenter(const ScaleVolumeFromCenter &other);
    
    ~ScaleVolumeFromCenter();
    
    ScaleVolumeFromCenter& operator=(const ScaleVolumeFromCenter &other);
    
    bool operator==(const ScaleVolumeFromCenter &other) const;
    bool operator!=(const ScaleVolumeFromCenter &other) const;
    
    static const char* typeName();

    const Point& center() const;
    
    void setCenter(const PointRef &center);

    int setVolume(System &system,
                  const SireUnits::Dimension::Volume &volume,
                  const PropertyMap &map = PropertyMap()) const;
    
private:
    /** The point from which the volume will be scaled */
    SireFF::PointPtr scale_point;
};

typedef SireBase::PropPtr<VolumeChanger> VolumeChangerPtr;

}

Q_DECLARE_METATYPE( SireMove::NullVolumeChanger )
Q_DECLARE_METATYPE( SireMove::ScaleVolumeFromCenter )

SIRE_EXPOSE_CLASS( SireMove::VolumeChanger )
SIRE_EXPOSE_CLASS( SireMove::NullVolumeChanger )
SIRE_EXPOSE_CLASS( SireMove::ScaleVolumeFromCenter )

SIRE_END_HEADER

#endif
