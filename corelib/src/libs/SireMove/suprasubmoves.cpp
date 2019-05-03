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

#include "suprasubmoves.h"

#include "suprasubsystem.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireError/errors.h"

using namespace SireMove;
using namespace SireBase;
using namespace SireStream;

//////////
////////// Implementation of SupraSubMoves
//////////

static const RegisterMetaType<SupraSubMoves> r_suprasubmoves( MAGIC_ONLY,
                                                    "SireMove::SupraSubMoves" );
                                                    
/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                        const SupraSubMoves &suprasubmoves)
{
    writeHeader(ds, r_suprasubmoves, 1);
    
    ds << static_cast<const Property&>(suprasubmoves);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, SupraSubMoves &suprasubmoves)
{
    VersionID v = readHeader(ds, r_suprasubmoves);
    
    if (v == 1)
    {
        ds >> static_cast<Property&>(suprasubmoves);
    }
    else
        throw version_error(v, "1", r_suprasubmoves, CODELOC);
        
    return ds;
}

/** Constructor */
SupraSubMoves::SupraSubMoves() : Property()
{}

/** Copy constructor */
SupraSubMoves::SupraSubMoves(const SupraSubMoves &other) : Property(other)
{}

/** Destructor */
SupraSubMoves::~SupraSubMoves()
{}

/** Copy assignment operator */
SupraSubMoves& SupraSubMoves::operator=(const SupraSubMoves &other)
{
    return *this;
}

/** Comparison operator */
bool SupraSubMoves::operator==(const SupraSubMoves &other) const
{
    return true;
}

/** Comparison operator */
bool SupraSubMoves::operator!=(const SupraSubMoves &other) const
{
    return false;
}

/** Return the number of types of SupraSubMove objects in this set */
int SupraSubMoves::count() const
{
    return this->subMoves().count();
}

/** Return the number of types of SupraSubMove objects in this set */
int SupraSubMoves::size() const
{
    return this->count();
}

/** Return the number of types of SupraSubMove objects in this set */
int SupraSubMoves::nSubMoveTypes() const
{
    return this->count();
}

Q_GLOBAL_STATIC( SameSupraSubMoves, sameSupraSubMoves )

/** Return the global null SupraSubMoves object */
const SameSupraSubMoves& SupraSubMoves::null()
{
    return *(sameSupraSubMoves());
}

//////////
////////// Implementation of SameSupraSubMoves
//////////

static const RegisterMetaType<SameSupraSubMoves> r_samesuprasubmoves;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                        const SameSupraSubMoves &samesuprasubmoves)
{
    writeHeader(ds, r_samesuprasubmoves, 1);
    
    SharedDataStream sds(ds);
    
    sds << samesuprasubmoves.submove
        << static_cast<const SupraSubMoves&>(samesuprasubmoves);
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds,
                                        SameSupraSubMoves &samesuprasubmoves)
{
    VersionID v = readHeader(ds, r_samesuprasubmoves);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> samesuprasubmoves.submove
            >> static_cast<SupraSubMoves&>(samesuprasubmoves);
    }
    else
        throw version_error(v, "1", r_samesuprasubmoves, CODELOC);
        
    return ds;
}

/** Constructor */
SameSupraSubMoves::SameSupraSubMoves()
                  : ConcreteProperty<SameSupraSubMoves,SupraSubMoves>()
{}

/** Construct to perform the move 'move' repeatedly on a SupraSubSystem */
SameSupraSubMoves::SameSupraSubMoves(const SupraSubMove &move)
                  : ConcreteProperty<SameSupraSubMoves,SupraSubMoves>(),
                    submove(move)
{}

/** Copy constructor */
SameSupraSubMoves::SameSupraSubMoves(const SameSupraSubMoves &other)
                  : ConcreteProperty<SameSupraSubMoves,SupraSubMoves>(other),
                    submove(other.submove)
{}

/** Destructor */
SameSupraSubMoves::~SameSupraSubMoves()
{}

/** Copy assignment operator */
SameSupraSubMoves& SameSupraSubMoves::operator=(const SameSupraSubMoves &other)
{
    SupraSubMoves::operator=(other);
    submove = other.submove;
    
    return *this;
}

/** Comparison operator */
bool SameSupraSubMoves::operator==(const SameSupraSubMoves &other) const
{
    return submove == other.submove and SupraSubMoves::operator==(other);
}

/** Comparison operator */
bool SameSupraSubMoves::operator!=(const SameSupraSubMoves &other) const
{
    return submove != other.submove or SupraSubMoves::operator!=(other);
}

/** Return the ith type of move

    \throw SireError::invalid_index
*/
const SupraSubMove& SameSupraSubMoves::operator[](int i) const
{
    if (i != 0)
        throw SireError::invalid_index( QObject::tr(
            "There is only a single SupraSubMove in the SameSupraSubMoves, so "
            "only index 0 is valid (index %1 is thus invalid)")
                .arg(i), CODELOC );

    return submove.read();
}

/** Return a string representation of this set of moves */
QString SameSupraSubMoves::toString() const
{
    return QObject::tr("SameSupraSubMoves{ %1 }").arg( submove.read().toString() );
}

/** Perform the moves */
void SameSupraSubMoves::move(SupraSubSystem &system, int nsubmoves, 
                             int nsubmoves_per_block, bool record_substats)
{
    if (nsubmoves <= 0)
        return;
        
    submove.edit().move(system, nsubmoves, 
                        nsubmoves_per_block, record_substats);
}

/** Clear all of the move statistics */
void SameSupraSubMoves::clearStatistics()
{
    submove.edit().clearStatistics();
}

/** Return a list of all of the types of submove in this set */
QList<SupraSubMovePtr> SameSupraSubMoves::subMoves() const
{
    QList<SupraSubMovePtr> moves;
    moves.append(submove);
    
    return moves;
}

const char* SameSupraSubMoves::typeName()
{
    return QMetaType::typeName( qMetaTypeId<SameSupraSubMoves>() );
}
