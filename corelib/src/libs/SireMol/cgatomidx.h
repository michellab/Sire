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

#ifndef SIREMOL_CGATOMIDX_H
#define SIREMOL_CGATOMIDX_H

#include "cgidx.h"
#include "atomidx.h"

#include "SireID/index.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class CGAtomIdx;
}

QDataStream& operator<<(QDataStream&, const SireMol::CGAtomIdx&);
QDataStream& operator>>(QDataStream&, SireMol::CGAtomIdx&);

namespace SireMol
{

/** This is the basic type used to ID atoms within a molecule. This
    provides the fastest way of indexing atoms and is the base
    type that all other AtomID classes map to.
    
    @author Christopher Woods
*/
class SIREMOL_EXPORT CGAtomIdx : public AtomID
{

friend QDataStream& ::operator<<(QDataStream&, const CGAtomIdx&);
friend QDataStream& ::operator>>(QDataStream&, CGAtomIdx&);

public:
    CGAtomIdx();
    
    CGAtomIdx(CGIdx cgid, SireID::Index atomid);
    
    CGAtomIdx(const CGAtomIdx &other);
    
    ~CGAtomIdx();
    
    static const char* typeName();
    
    const char* what() const
    {
        return CGAtomIdx::typeName();
    }
    
    CGAtomIdx* clone() const;
    
    static CGAtomIdx null();
    
    bool isNull() const;
    
    uint hash() const;
    
    QString toString() const;
    
    CGAtomIdx& operator=(const CGAtomIdx &other);
    
    bool operator==(const SireID::ID &other) const;
    
    bool operator==(const CGAtomIdx &other) const;
    
    bool operator!=(const CGAtomIdx &other) const;
    
    QList<AtomIdx> map(const MolInfo &molinfo) const;
    
    CGIdx cutGroup() const;
    
    SireID::Index atom() const;
    
private:
    /** The index of the CutGroup that contains the atom */
    CGIdx _cgidx;
    
    /** The index of the atom within the CutGroup */
    SireID::Index _atmidx;
};

}

Q_DECLARE_METATYPE(SireMol::CGAtomIdx);

SIRE_EXPOSE_CLASS( SireMol::CGAtomIdx )

SIRE_END_HEADER

#endif
