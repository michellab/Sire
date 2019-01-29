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

#ifndef SIREMOL_RESIDENTIFIER_H
#define SIREMOL_RESIDENTIFIER_H

#include "resid.h"

#include <boost/shared_ptr.hpp>

namespace SireMol
{
class ResIdentifier;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::ResIdentifier&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::ResIdentifier&);

namespace SireMol
{

/** This is the polymorphic holder of all residue IDs */
class SIREMOL_EXPORT ResIdentifier : public ResID
{

friend QDataStream& ::operator<<(QDataStream&, const ResIdentifier&);
friend QDataStream& ::operator>>(QDataStream&, ResIdentifier&);

public:
    ResIdentifier();
    ResIdentifier(const ResID &resid);
    ResIdentifier(const ResIdentifier &other);
    
    ~ResIdentifier();
    
    static const char* typeName();
    
    const char* what() const
    {
        return ResIdentifier::typeName();
    }
    
    ResIdentifier* clone() const;
    
    bool isNull() const;
    
    uint hash() const;
                
    QString toString() const;
    
    const ResID& base() const;
    
    ResIdentifier& operator=(const ResIdentifier &other);
    ResIdentifier& operator=(const ResID &other);
    
    bool operator==(const SireID::ID &other) const;
    using SireID::ID::operator!=;
   
    bool operator==(const ResIdentifier &other) const;
    bool operator!=(const ResIdentifier &other) const;
    
    bool operator==(const ResID &other) const;
    bool operator!=(const ResID &other) const;
    
    QList<ResIdx> map(const MolInfo &molinfo) const;

private:
    /** Pointer to the ResID */
    boost::shared_ptr<ResID> d;
};

inline uint qHash(const ResIdentifier &resid)
{
    return resid.hash();
}

}

#include "residx.h"
#include "atomidx.h"

Q_DECLARE_METATYPE( SireID::Specify<SireMol::ResID> );
Q_DECLARE_METATYPE( SireMol::AtomsIn<SireMol::ResID> );
Q_DECLARE_METATYPE( SireID::IDAndSet<SireMol::ResID> );
Q_DECLARE_METATYPE( SireID::IDOrSet<SireMol::ResID> );
Q_DECLARE_METATYPE( SireID::MatchAll<SireMol::ResID> );
Q_DECLARE_METATYPE( SireID::InvertMatch<SireMol::ResID> );

Q_DECLARE_METATYPE(SireMol::ResIdentifier);

#endif

