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

#include "checkpoint.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireSystem;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<CheckPoint> r_ckpt;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const CheckPoint &ckpt)
{
    writeHeader(ds, r_ckpt, 1);
    
    SharedDataStream sds(ds);
    sds << ckpt.old_system;
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, CheckPoint &ckpt)
{
    VersionID v = readHeader(ds, r_ckpt);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> ckpt.old_system;
    }
    else
        throw version_error(v, "1", r_ckpt, CODELOC);
        
    return ds;
}

/** Create a null (empty) checkpoint */
CheckPoint::CheckPoint() : ConcreteProperty<CheckPoint,Property>(),
                           old_system( System::null() )
{}

/** Construct a checkpoint to hold the current state of the system 'system' */
CheckPoint::CheckPoint(const System &system)
           : ConcreteProperty<CheckPoint,Property>(),
             old_system(system)
{}

/** Copy constructor */
CheckPoint::CheckPoint(const CheckPoint &other)
           : ConcreteProperty<CheckPoint,Property>(other),
             old_system(other.old_system)
{}

/** Destructor */
CheckPoint::~CheckPoint()
{}

/** Set this equal to the checkpoint of the system 'system' */
CheckPoint& CheckPoint::operator=(const System &system)
{
    old_system = system;
    return *this;
}

/** Copy assignment operator */
CheckPoint& CheckPoint::operator=(const CheckPoint &other)
{
    old_system = other.old_system;
    Property::operator=(other);
    
    return *this;
}

/** Comparison operator */
bool CheckPoint::operator==(const CheckPoint &other) const
{
    return old_system == other.old_system;
}

/** Comparison operator */
bool CheckPoint::operator!=(const CheckPoint &other) const
{
    return old_system != other.old_system;
}

/** Allow automatic conversion of a CheckPoint to a System */
CheckPoint::operator System() const
{
    return old_system;
}

const char* CheckPoint::typeName()
{
    return QMetaType::typeName( qMetaTypeId<CheckPoint>() );
}
