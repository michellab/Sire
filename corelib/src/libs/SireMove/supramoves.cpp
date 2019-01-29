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

#include "supramoves.h"

#include "suprasystem.h"

#include "SireID/index.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMove;
using namespace SireID;
using namespace SireBase;
using namespace SireStream;

/////////
///////// Implementation of SupraMoves
/////////

static const RegisterMetaType<SupraMoves> r_supramoves( MAGIC_ONLY,
                                                    "SireMove::SupraMoves" );
                                                    
/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const SupraMoves &supramoves)
{
    writeHeader(ds, r_supramoves, 1);
    
    ds << static_cast<const Property&>(supramoves);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, SupraMoves &supramoves)
{
    VersionID v = readHeader(ds, r_supramoves);
    
    if (v == 1)
    {
        ds >> static_cast<Property&>(supramoves);
    }
    else
        throw version_error(v, "1", r_supramoves, CODELOC);
        
    return ds;
}

/** Constructor */
SupraMoves::SupraMoves() : Property()
{}

/** Copy constructor */
SupraMoves::SupraMoves(const SupraMoves &other) : Property(other)
{}

/** Destructor */
SupraMoves::~SupraMoves()
{}

/** Copy assignment operator */
SupraMoves& SupraMoves::operator=(const SupraMoves &other)
{
    return *this;
}

/** Comparison operator */
bool SupraMoves::operator==(const SupraMoves &other) const
{
    return true;
}

/** Comparison operator */
bool SupraMoves::operator!=(const SupraMoves &other) const
{
    return false;
}

/** Return the ith move 

    \throw SireError::invalid_index
*/
const SupraMove& SupraMoves::operator[](int i) const
{
    QList<SupraMovePtr> mvs = this->moves();
    
    return mvs.at( Index(i).map(mvs.count()) ).read();
}

/** Return the number of different types of move in this set */
int SupraMoves::count() const
{
    return this->moves().count();
}

/** Return the number of different types of move in this set */
int SupraMoves::size() const
{
    return this->count();
}

/** Return the number of different types of move in this set */
int SupraMoves::nSubMoveTypes() const
{
    return this->count();
}

Q_GLOBAL_STATIC( SameSupraMoves, sameSupraMoves )

/** Return the global null SupraMoves object */
const SameSupraMoves& SupraMoves::null()
{
    return *(sameSupraMoves());
}

/////////
///////// Implementation of SameSupraMoves
/////////

static RegisterMetaType<SameSupraMoves> r_samesupramoves;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                        const SameSupraMoves &samesupramoves)
{
    writeHeader(ds, r_samesupramoves, 1);
    
    SharedDataStream sds(ds);
    
    sds << samesupramoves.mv
        << static_cast<const SupraMoves&>(samesupramoves);
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, SameSupraMoves &samesupramoves)
{
    VersionID v = readHeader(ds, r_samesupramoves);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> samesupramoves.mv
            >> static_cast<SupraMoves&>(samesupramoves);
    }
    else
        throw version_error(v, "1", r_samesupramoves, CODELOC);
        
    return ds;
}

/** Constructor */
SameSupraMoves::SameSupraMoves() : ConcreteProperty<SameSupraMoves,SupraMoves>()
{}

/** Construct to run the move 'move' repeatedly */
SameSupraMoves::SameSupraMoves(const SupraMove &move)
               : ConcreteProperty<SameSupraMoves,SupraMoves>(), mv(move)
{}

/** Copy constructor */
SameSupraMoves::SameSupraMoves(const SameSupraMoves &other)
               : ConcreteProperty<SameSupraMoves,SupraMoves>(other),
                 mv(other.mv)
{}

/** Destructor */
SameSupraMoves::~SameSupraMoves()
{}

/** Copy assignment operator */
SameSupraMoves& SameSupraMoves::operator=(const SameSupraMoves &other)
{
    SupraMoves::operator=(other);
    
    mv = other.mv;
    
    return *this;
}

/** Comparison operator */
bool SameSupraMoves::operator==(const SameSupraMoves &other) const
{
    return mv == other.mv and SupraMoves::operator==(other);
}

/** Comparison operator */
bool SameSupraMoves::operator!=(const SameSupraMoves &other) const
{
    return mv != other.mv or SupraMoves::operator!=(other);
}

/** Return a string representation of this moves set */
QString SameSupraMoves::toString() const
{
    return QObject::tr( "SameSupraMoves{ %1 }" ).arg(mv->toString());
}

/** Perform the moves 'nmoves' times */
void SameSupraMoves::move(SupraSystem &system, int nmoves, bool record_stats)
{
    if (nmoves == 0)
        return;
        
    mv.edit().move(system, nmoves, record_stats);
}

/** Clear all move statistics */
void SameSupraMoves::clearStatistics()
{
    mv.edit().clearStatistics();
}

/** Return a list of all of the moves */
QList<SupraMovePtr> SameSupraMoves::moves() const
{
    QList<SupraMovePtr> mvs;
    mvs.append(mv);
    
    return mvs;
}

/** Return the total number of moves that have been performed */
int SameSupraMoves::nMoves() const
{
    return mv.read().nMoves();
}

const char* SameSupraMoves::typeName()
{
    return QMetaType::typeName( qMetaTypeId<SameSupraMoves>() );
}
