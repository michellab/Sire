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

#include "suprasubmove.h"

#include "suprasubsystem.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMove;
using namespace SireBase;
using namespace SireStream;

//////////
////////// Implementation of SupraSubMove
//////////

static const RegisterMetaType<SupraSubMove> r_suprasubmove( MAGIC_ONLY,
                                                "SireMove::SupraSubMove" );
                                                
/** Serialise to a binary datastream */
QDataStream SIREMOVE_EXPORT &operator<<(QDataStream &ds, 
                                        const SupraSubMove &suprasubmove)
{
    writeHeader(ds, r_suprasubmove, 1);
    
    ds << static_cast<const Property&>(suprasubmove);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOVE_EXPORT &operator>>(QDataStream &ds, SupraSubMove &suprasubmove)
{
    VersionID v = readHeader(ds, r_suprasubmove);
    
    if (v == 1)
    {
        ds >> static_cast<Property&>(suprasubmove);
    }
    else
        throw version_error(v, "1", r_suprasubmove, CODELOC);
        
    return ds;
}

/** Constructor */
SupraSubMove::SupraSubMove() : Property()
{}

/** Copy constructor */
SupraSubMove::SupraSubMove(const SupraSubMove &other) : Property(other)
{}

/** Destructor */
SupraSubMove::~SupraSubMove()
{}

/** Copy assignment operator */
SupraSubMove& SupraSubMove::operator=(const SupraSubMove &other)
{
    return *this;
}

/** Comparison operator */
bool SupraSubMove::operator==(const SupraSubMove &other) const
{
    return true;
}

/** Comparison operator */
bool SupraSubMove::operator!=(const SupraSubMove &other) const
{
    return false;
}

/** Clear the move statistics */
void SupraSubMove::clearStatistics()
{}

Q_GLOBAL_STATIC( NullSupraSubMove, nullSupraSubMove )

/** Return the global null SupraSubMove */
const NullSupraSubMove& SupraSubMove::null()
{
    return *(nullSupraSubMove());
}

//////////
////////// Implementation of NullSupraSubMove
//////////

static const RegisterMetaType<NullSupraSubMove> r_nullsuprasubmove;

/** Serialise to a binary datastream */
QDataStream SIREMOVE_EXPORT &operator<<(QDataStream &ds,
                                        const NullSupraSubMove &nullsuprasubmove)
{
    writeHeader(ds, r_nullsuprasubmove, 1);
    
    ds << static_cast<const SupraSubMove&>(nullsuprasubmove);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOVE_EXPORT &operator>>(QDataStream &ds,
                                        NullSupraSubMove &nullsuprasubmove)
{
    VersionID v = readHeader(ds, r_nullsuprasubmove);
    
    if (v == 1)
    {
        ds >> static_cast<SupraSubMove&>(nullsuprasubmove);
    }
    else
        throw version_error(v, "1", r_nullsuprasubmove, CODELOC);
        
    return ds;
}

/** Constructor */
NullSupraSubMove::NullSupraSubMove()
                 : ConcreteProperty<NullSupraSubMove,SupraSubMove>()
{}

/** Copy constructor */
NullSupraSubMove::NullSupraSubMove(const NullSupraSubMove &other)
                 : ConcreteProperty<NullSupraSubMove,SupraSubMove>(other)
{}

/** Destructor */
NullSupraSubMove::~NullSupraSubMove()
{}

/** Copy assignment operator */
NullSupraSubMove& NullSupraSubMove::operator=(const NullSupraSubMove &other)
{
    SupraSubMove::operator=(other);
    return *this;
}

/** Comparison operator */
bool NullSupraSubMove::operator==(const NullSupraSubMove &other) const
{
    return true;
}

/** Comparison operator */
bool NullSupraSubMove::operator!=(const NullSupraSubMove &other) const
{
    return false;
}

/** Return a string representation of this move */
QString NullSupraSubMove::toString() const
{
    return QObject::tr("NullSupraSubMove()");
}

/** Null move, so doesn't do anything */
void NullSupraSubMove::move(SupraSubSystem&, int, int, bool)
{}

const char* NullSupraSubMove::typeName()
{
    return QMetaType::typeName( qMetaTypeId<NullSupraSubMove>() );
}
