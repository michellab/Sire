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

#ifndef SIREFF_FFNAME_H
#define SIREFF_FFNAME_H

#include "SireID/name.h"

#include "ffid.h"

SIRE_BEGIN_HEADER

namespace SireFF
{
class FFName;
}

SIREFF_EXPORT QDataStream& operator<<(QDataStream&, const SireFF::FFName&);
SIREFF_EXPORT QDataStream& operator>>(QDataStream&, SireFF::FFName&);

namespace SireFF
{

/** This class holds the name of a forcefield
    
    @author Christopher Woods
*/
class SIREFF_EXPORT FFName : public SireID::Name, public FFID
{

friend SIREFF_EXPORT QDataStream& ::operator<<(QDataStream&, const FFName&);
friend SIREFF_EXPORT QDataStream& ::operator>>(QDataStream&, FFName&);

public:
    FFName();
    explicit FFName(const QString &name);
    
    FFName(const QString &name, SireID::CaseSensitivity case_sensitivity);
    
    FFName(const FFName &other);
    
    ~FFName();
    
    static const char* typeName();
    
    const char* what() const
    {
        return FFName::typeName();
    }
    
    FFName* clone() const;
    
    bool isNull() const;
    
    uint hash() const;
    
    QString toString() const;
    
    FFName& operator=(const FFName &other);
    
    bool operator==(const SireID::ID &other) const;
    
    bool operator==(const FFName &other) const;
    
    bool operator!=(const FFName &other) const;
    
    QList<FFIdx> map(const ForceFields &ffields) const;
};

}

Q_DECLARE_METATYPE(SireFF::FFName);

SIRE_EXPOSE_CLASS( SireFF::FFName )

SIRE_END_HEADER

#endif
