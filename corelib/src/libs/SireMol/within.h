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

#ifndef SIREMOL_WITHIN_H
#define SIREMOL_WITHIN_H

#include "atomid.h"
#include "atomidentifier.h"

#include "SireMaths/vector.h"

#include "SireUnits/dimensions.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class Within;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::Within&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::Within&);

namespace SireMol
{

using SireMaths::Vector;

/** This is an atom identifier that identifies atoms
    based on how far they are from other atoms, or
    from points in space
    
    @author Christopher Woods
*/
class SIREMOL_EXPORT Within : public AtomID
{

friend QDataStream& ::operator<<(QDataStream&, const Within&);
friend QDataStream& ::operator>>(QDataStream&, Within&);

public:
    Within();
    Within(SireUnits::Dimension::Length distance, 
           const Vector &point);
    
    Within(SireUnits::Dimension::Length distance,
           const AtomID &atomid);
    
    Within(const Within &other);
    
    ~Within();
    
    static const char* typeName();
    
    const char* what() const
    {
        return Within::typeName();
    }
    
    Within* clone() const;
    
    bool isNull() const;
    
    uint hash() const;
                
    QString toString() const;
    
    Within& operator=(const Within &other);
    
    bool operator==(const SireID::ID &other) const;
    using SireID::ID::operator!=;
   
    bool operator==(const Within &other) const;
    bool operator!=(const Within &other) const;
    
    QList<AtomIdx> map(const MolInfo &molinfo) const;

    QList<AtomIdx> map(const MoleculeView &molview, const PropertyMap &map) const;

private:
    /** The atom against which distance will be calculated */
    AtomIdentifier atomid;
    
    /** The point in space against which distance will be calculated */
    Vector point;
    
    /** The distance itself */
    SireUnits::Dimension::Length dist;
};

}

Q_DECLARE_METATYPE( SireMol::Within )

SIRE_EXPOSE_CLASS( SireMol::Within )

SIRE_END_HEADER

#endif
