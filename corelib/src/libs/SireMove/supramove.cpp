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

#include "supramove.h"

#include "suprasystem.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMove;
using namespace SireBase;
using namespace SireStream;

/////////////
///////////// Implementation of SupraMove
/////////////

static const RegisterMetaType<SupraMove> r_supramove(MAGIC_ONLY,
                                                "SireMove::SupraMove");
                                                
/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const SupraMove &supramove)
{
    writeHeader(ds, r_supramove, 1);
    
    ds << supramove.nmoves
       << static_cast<const Property&>(supramove);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, SupraMove &supramove)
{
    VersionID v = readHeader(ds, r_supramove);
    
    if (v == 1)
    {
        ds >> supramove.nmoves
           >> static_cast<Property&>(supramove);
    }
    else
        throw version_error(v, "1", r_supramove, CODELOC);
        
    return ds;
}

/** Constructor */
SupraMove::SupraMove() : Property(), nmoves(0)
{}

/** Copy constructor */
SupraMove::SupraMove(const SupraMove &other) 
          : Property(other), nmoves(other.nmoves)
{}

/** Destructor */
SupraMove::~SupraMove()
{}

/** Copy assignment operator */
SupraMove& SupraMove::operator=(const SupraMove &other)
{
    nmoves = other.nmoves;
    return *this;
}

/** Comparison operator */
bool SupraMove::operator==(const SupraMove &other) const
{
    return true;
}

/** Comparison operator */
bool SupraMove::operator!=(const SupraMove &other) const
{
    return false;
}

/** Clear all move statistics */
void SupraMove::clearStatistics()
{
    nmoves = 0;
}

/** Return the total number of supra-moves performed using this object */
int SupraMove::nMoves() const
{
    return nmoves;
}

/** Internal function used to increment the total number of moves performed */
void SupraMove::incrementNMoves(int nmore)
{
    if (nmore > 0)
        nmoves += nmore;
}

Q_GLOBAL_STATIC( NullSupraMove, nullSupraMove )

/** Return the global null move */
const NullSupraMove& SupraMove::null()
{
    return *(nullSupraMove());
}

/////////////
///////////// Implementation of NullSupraMove
/////////////

static const RegisterMetaType<NullSupraMove> r_nullsupramove;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                        const NullSupraMove &nullsupramove)
{
    writeHeader(ds, r_nullsupramove, 1);
    
    ds << static_cast<const SupraMove&>(nullsupramove);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds,
                                        NullSupraMove &nullsupramove)
{
    VersionID v = readHeader(ds, r_nullsupramove);
    
    if (v == 1)
    {
        ds >> static_cast<SupraMove&>(nullsupramove);
    }
    else
        throw version_error(v, "1", r_nullsupramove, CODELOC);
        
    return ds;
}

/** Constructor */
NullSupraMove::NullSupraMove() : ConcreteProperty<NullSupraMove,SupraMove>()
{}

/** Copy constructor */
NullSupraMove::NullSupraMove(const NullSupraMove &other)
              : ConcreteProperty<NullSupraMove,SupraMove>(other)
{}

/** Destructor */
NullSupraMove::~NullSupraMove()
{}

/** Copy assignment operator */
NullSupraMove& NullSupraMove::operator=(const NullSupraMove &other)
{
    SupraMove::operator=(other);
    return *this;
}

/** Comparison operator */
bool NullSupraMove::operator==(const NullSupraMove &other) const
{
    return true;
}

/** Comparison operator */
bool NullSupraMove::operator!=(const NullSupraMove &other) const
{
    return false;
}

/** Return a string representation of this move */
QString NullSupraMove::toString() const
{
    return QObject::tr("NullSupraMove()");
}

/** Perform the null move - this does nothing! */
void NullSupraMove::move(SupraSystem&, int, bool)
{}

const char* NullSupraMove::typeName()
{
    return QMetaType::typeName( qMetaTypeId<NullSupraMove>() );
}
