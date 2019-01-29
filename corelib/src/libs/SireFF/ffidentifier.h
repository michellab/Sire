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

#ifndef SIREFF_FFIDENTIFIER_H
#define SIREFF_FFIDENTIFIER_H

#include "ffid.h"

#include <boost/shared_ptr.hpp>

namespace SireFF
{
class FFIdentifier;
}

SIREFF_EXPORT QDataStream& operator<<(QDataStream&, const SireFF::FFIdentifier&);
SIREFF_EXPORT QDataStream& operator>>(QDataStream&, SireFF::FFIdentifier&);

namespace SireFF
{

/** This is a generic holder for any FFID class! 

    @author Christopher Woods
*/
class SIREFF_EXPORT FFIdentifier : public FFID
{

friend QDataStream& ::operator<<(QDataStream&, const FFIdentifier&);
friend QDataStream& ::operator>>(QDataStream&, FFIdentifier&);

public:
    FFIdentifier();
    FFIdentifier(const FFID &ffid);
    FFIdentifier(const FFIdentifier &other);
    
    ~FFIdentifier();
    
    static const char* typeName();
    
    const char* what() const
    {
        return FFIdentifier::typeName();
    }
    
    FFIdentifier* clone() const;
        
    bool isNull() const;
    
    uint hash() const;
                
    QString toString() const;
    
    const FFID& base() const;
    
    FFIdentifier& operator=(const FFIdentifier &other);
    FFIdentifier& operator=(const FFID &other);
    
    bool operator==(const SireID::ID &other) const;
    using SireID::ID::operator!=;
   
    bool operator==(const FFIdentifier &other) const;
    bool operator!=(const FFIdentifier &other) const;
    
    bool operator==(const FFID &other) const;
    bool operator!=(const FFID &other) const;

    QList<FFIdx> map(const ForceFields &ffields) const;

private:
    /** Pointer to the FFID */
    boost::shared_ptr<FFID> d;
};

inline uint qHash(const FFIdentifier &ffid)
{
    return ffid.hash();
}

}

#include "ffidx.h"

Q_DECLARE_METATYPE( SireID::Specify<SireFF::FFID> )
Q_DECLARE_METATYPE( SireID::IDAndSet<SireFF::FFID> )
Q_DECLARE_METATYPE( SireID::IDOrSet<SireFF::FFID> )

Q_DECLARE_METATYPE(SireFF::FFIdentifier);

#endif
