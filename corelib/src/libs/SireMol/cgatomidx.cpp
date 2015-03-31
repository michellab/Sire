/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2008  Christopher Woods
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

#include "cgatomidx.h"
#include "molinfo.h"

using namespace SireMol;

static const RegisterMetaType<CGAtomIdx> r_cgatomidx;

QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds, 
                                       const CGAtomIdx &cgatomidx)
{
    ds << cgatomidx._cgidx << cgatomidx._atmidx;
    return ds;
}

QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds,
                                       CGAtomIdx &cgatomidx)
{
    ds >> cgatomidx._cgidx >> cgatomidx._atmidx;
    return ds;
}

CGAtomIdx::CGAtomIdx() : AtomID()
{}

CGAtomIdx::CGAtomIdx(CGIdx cgid, SireID::Index atomid)
      : AtomID(), _cgidx(cgid), _atmidx(atomid)
{}

CGAtomIdx::CGAtomIdx(const CGAtomIdx &other)
      : AtomID(other), _cgidx(other._cgidx), _atmidx(other._atmidx)
{}

CGAtomIdx::~CGAtomIdx()
{}

CGAtomIdx CGAtomIdx::null()
{
    return CGAtomIdx( CGIdx::null(), SireID::Index::null() );
} 

bool CGAtomIdx::isNull() const
{
    return _cgidx.isNull() and _atmidx.isNull();
}

uint CGAtomIdx::hash() const
{
    return (::qHash(_cgidx) << 16) | (::qHash(_atmidx) & 0x0000FFFF);
}

QString CGAtomIdx::toString() const
{
    return QString("{%1,%2}").arg(_cgidx.toString(), _atmidx.toString());
}

CGAtomIdx& CGAtomIdx::operator=(const CGAtomIdx &other)
{
    _cgidx = other._cgidx;
    _atmidx = other._atmidx;
    
    AtomID::operator=(other);
    
    return *this;
}

bool CGAtomIdx::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare<CGAtomIdx>(*this, other);
}

bool CGAtomIdx::operator==(const CGAtomIdx &other) const
{
    return _cgidx == other._cgidx and _atmidx == other._atmidx;
}

bool CGAtomIdx::operator!=(const CGAtomIdx &other) const
{
    return _cgidx != other._cgidx or _atmidx != other._atmidx;
}

CGIdx CGAtomIdx::cutGroup() const
{
    return _cgidx;
}

SireID::Index CGAtomIdx::atom() const
{
    return _atmidx;
}

/** Combine a CGIdx with an Index */
CGAtomIdx CGIdx::operator+(const SireID::Index &other) const
{
    return CGAtomIdx(*this, other);
}

/** Combine a CGIdx with an Index */
namespace SireMol
{

CGAtomIdx SIREMOL_EXPORT operator+(const SireID::Index &idx, const CGIdx &cgidx)
{
    return CGAtomIdx(cgidx, idx);
}

}

QList<AtomIdx> CGAtomIdx::map(const MolInfo &molinfo) const
{
    QList<AtomIdx> atomidxs;
    atomidxs.append( molinfo.getAtom(_cgidx, _atmidx) );
    return atomidxs;
}

const char* CGAtomIdx::typeName()
{
    return QMetaType::typeName( qMetaTypeId<CGAtomIdx>() );
}

CGAtomIdx* CGAtomIdx::clone() const
{
    return new CGAtomIdx(*this);
}

