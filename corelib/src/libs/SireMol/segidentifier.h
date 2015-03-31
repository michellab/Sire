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

#ifndef SIREMOL_SEGIDENTIFIER_H
#define SIREMOL_SEGIDENTIFIER_H

#include "segid.h"

#include <boost/shared_ptr.hpp>

namespace SireMol
{
class SegIdentifier;
}

QDataStream& operator<<(QDataStream&, const SireMol::SegIdentifier&);
QDataStream& operator>>(QDataStream&, SireMol::SegIdentifier&);

namespace SireMol
{

class SegIdx;

/** This is the polymorphic holder of all Segment IDs */
class SIREMOL_EXPORT SegIdentifier : public SegID
{

friend QDataStream& ::operator<<(QDataStream&, const SegIdentifier&);
friend QDataStream& ::operator>>(QDataStream&, SegIdentifier&);

public:
    SegIdentifier();
    SegIdentifier(const SegID &segid);
    SegIdentifier(const SegIdentifier &other);
    
    ~SegIdentifier();
    
    static const char* typeName();
    
    const char* what() const
    {
        return SegIdentifier::typeName();
    }
    
    SegIdentifier* clone() const;
    
    bool isNull() const;
    
    uint hash() const;
                
    QString toString() const;
    
    const SegID& base() const;
    
    SegIdentifier& operator=(const SegIdentifier &other);
    SegIdentifier& operator=(const SegID &other);
    
    bool operator==(const SireID::ID &other) const;
    using SireID::ID::operator!=;
   
    bool operator==(const SegIdentifier &other) const;
    bool operator!=(const SegIdentifier &other) const;
    
    bool operator==(const SegID &other) const;
    bool operator!=(const SegID &other) const;
    
    QList<SegIdx> map(const MolInfo &molinfo) const;

private:
    /** Pointer to the SegID */
    boost::shared_ptr<SegID> d;
};

inline uint qHash(const SegIdentifier &segid)
{
    return segid.hash();
}

}

#include "segidx.h"
#include "atomidx.h"

Q_DECLARE_METATYPE( SireID::Specify<SireMol::SegID> );
Q_DECLARE_METATYPE( SireMol::AtomsIn<SireMol::SegID> );
Q_DECLARE_METATYPE( SireID::IDAndSet<SireMol::SegID> );
Q_DECLARE_METATYPE( SireID::IDOrSet<SireMol::SegID> );


Q_DECLARE_METATYPE(SireMol::SegIdentifier);

#endif
