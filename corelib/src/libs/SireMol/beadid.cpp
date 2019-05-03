/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2010  Christopher Woods
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

#include "beadid.h"
#include "beadidx.h"
#include "beadnum.h"

#include "SireStream/datastream.h"

using namespace SireMol;
using namespace SireID;
using namespace SireStream;

BeadID::BeadID() : ID()
{}

BeadID::BeadID(const BeadID &other) : ID(other)
{}

BeadID::~BeadID()
{}

/////////
///////// Implementation of BeadIdx
/////////

static const RegisterMetaType<BeadIdx> r_beadidx;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const BeadIdx &beadidx)
{
    writeHeader(ds, r_beadidx, 1);
    
    ds << static_cast<const SireID::Index_T_<BeadIdx>&>(beadidx);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, BeadIdx &beadidx)
{
    VersionID v = readHeader(ds, r_beadidx);
    
    if (v == 1)
    {
        ds >> static_cast<SireID::Index_T_<BeadIdx>&>(beadidx);
    }
    else
        throw version_error( v, "1", r_beadidx, CODELOC );
        
    return ds;
}

BeadIdx::BeadIdx() : SireID::Index_T_<BeadIdx>(), BeadID()
{}

BeadIdx::BeadIdx(qint32 idx) 
          : SireID::Index_T_<BeadIdx>(idx), BeadID()
{}

BeadIdx::BeadIdx(const BeadIdx &other) 
          : SireID::Index_T_<BeadIdx>(other), BeadID(other)
{}

BeadIdx::~BeadIdx()
{}

BeadIdx BeadIdx::null()
{
    return BeadIdx();
}

bool BeadIdx::isNull() const
{
    return SireID::Index_T_<BeadIdx>::isNull();
}

uint BeadIdx::hash() const
{
    return SireID::Index_T_<BeadIdx>::hash();
}

BeadIdx* BeadIdx::clone() const
{
    return new BeadIdx(*this);
}

QString BeadIdx::toString() const
{
    return QString("BeadIdx(%1)").arg(_idx);
}

BeadIdx& BeadIdx::operator=(const BeadIdx &other)
{
    SireID::IndexBase::operator=(other);
    BeadID::operator=(other);
    return *this;
}

bool BeadIdx::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare<BeadIdx>(*this, other);
}

const char* BeadIdx::typeName()
{
    return QMetaType::typeName( qMetaTypeId<BeadIdx>() );
}

///////
/////// Implementation of BeadNum
///////

static const RegisterMetaType<BeadNum> r_beadnum;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const BeadNum &beadnum)
{
    writeHeader(ds, r_beadnum, 1);
    
    ds << static_cast<const SireID::Number&>(beadnum);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, BeadNum &beadnum)
{
    VersionID v = readHeader(ds, r_beadnum);
    
    if (v == 1)
    {
        ds >> static_cast<SireID::Number&>(beadnum);
    }
    else
        throw version_error( v, "1", r_beadnum, CODELOC );
        
    return ds;
}

BeadNum::BeadNum() : SireID::Number(), BeadID()
{}

BeadNum::BeadNum(quint32 num) : SireID::Number(num), BeadID()
{}

BeadNum::BeadNum(const BeadNum &other) : SireID::Number(other), BeadID(other)
{}

BeadNum::~BeadNum()
{}

bool BeadNum::isNull() const
{
    return SireID::Number::isNull();
}

uint BeadNum::hash() const
{
    return ::qHash( static_cast<const SireID::Number&>(*this) );
}

QString BeadNum::toString() const
{
    return QString("BeadNum(%1)").arg(_num);
}

BeadNum& BeadNum::operator=(const BeadNum &other)
{
    SireID::Number::operator=(other);
    BeadID::operator=(other);
    return *this;
}

bool BeadNum::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare<BeadNum>(*this, other);
}

bool BeadNum::operator==(const BeadNum &other) const
{
    return _num == other._num;
}

bool BeadNum::operator!=(const BeadNum &other) const
{
    return _num != other._num;
}

bool BeadNum::operator<(const BeadNum &other) const
{
    return _num < other._num;
}

bool BeadNum::operator<=(const BeadNum &other) const
{
    return _num <= other._num;
}

bool BeadNum::operator>(const BeadNum &other) const
{
    return _num > other._num;
}

bool BeadNum::operator>=(const BeadNum &other) const
{
    return _num >= other._num;
}

BeadNum* BeadNum::clone() const
{
    return new BeadNum(*this);
}

const char* BeadNum::typeName()
{
    return QMetaType::typeName( qMetaTypeId<BeadNum>() );
}

