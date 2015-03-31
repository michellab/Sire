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

#ifndef SIREVOL_COMBINESPACES_H
#define SIREVOL_COMBINESPACES_H

#include "SireBase/combineproperties.h"

#include "combinedspace.h"

SIRE_BEGIN_HEADER

namespace SireVol
{
class CombineSpaces;
}

QDataStream& operator<<(QDataStream&, const SireVol::CombineSpaces&);
QDataStream& operator>>(QDataStream&, SireVol::CombineSpaces&);

namespace SireVol
{

using SireBase::PropertyName;
using SireBase::Properties;

/** This is a property which creates a SireVol::CombinedSpace object
    of the specified properties (which must all be space objects
    themselves)
    
    @author Christopher Woods
*/
class SIREVOL_EXPORT CombineSpaces 
        : public SireBase::ConcreteProperty<CombineSpaces,SireBase::CombineProperties>
{

friend QDataStream& ::operator<<(QDataStream&, const CombineSpaces&);
friend QDataStream& ::operator>>(QDataStream&, CombineSpaces&);

public:
    CombineSpaces();

    CombineSpaces(const PropertyName &source);
    CombineSpaces(const PropertyName &source0, const PropertyName &source1);
    
    CombineSpaces(const QList<PropertyName> &sources);
    CombineSpaces(const QVector<PropertyName> &sources);
    
    CombineSpaces(const QList<QString> &sources);
    CombineSpaces(const QVector<QString> &sources);
    
    CombineSpaces(const CombineSpaces &other);
    
    ~CombineSpaces();
    
    CombineSpaces& operator=(const CombineSpaces &other);
    
    bool operator==(const CombineSpaces &other) const;
    bool operator!=(const CombineSpaces &other) const;
    
    static const char* typeName();
    
    void updateFrom(const Properties &properties);
};

}

Q_DECLARE_METATYPE( SireVol::CombineSpaces )

SIRE_EXPOSE_CLASS( SireVol::CombineSpaces )

SIRE_END_HEADER

#endif
