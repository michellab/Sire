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

#ifndef SIREMOL_MOLIDENTIFIER_H
#define SIREMOL_MOLIDENTIFIER_H

#include "molid.h"

#include "atomsin.hpp"

#include "SireID/specify.hpp"
#include "SireID/idandset.hpp"
#include "SireID/idorset.hpp"

#include <boost/shared_ptr.hpp>

namespace SireMol
{
class MolIdentifier;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::MolIdentifier&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::MolIdentifier&);

namespace SireMol
{

/** This is a generic holder for any MolID class! 

    @author Christopher Woods
*/
class SIREMOL_EXPORT MolIdentifier : public MolID
{

friend SIREMOL_EXPORT QDataStream& ::operator<<(QDataStream&, const MolIdentifier&);
friend SIREMOL_EXPORT QDataStream& ::operator>>(QDataStream&, MolIdentifier&);

public:
    MolIdentifier();
    MolIdentifier(const MolID &atomid);
    MolIdentifier(const MolIdentifier &other);
    
    ~MolIdentifier();
    
    static const char* typeName();
    
    const char* what() const
    {
        return MolIdentifier::typeName();
    }
    
    MolIdentifier* clone() const;
    
    bool isNull() const;
    
    uint hash() const;
                
    QString toString() const;
    
    const MolID& base() const;
    
    MolIdentifier& operator=(const MolIdentifier &other);
    MolIdentifier& operator=(const MolID &other);
    
    bool operator==(const SireID::ID &other) const;
    using SireID::ID::operator!=;
   
    bool operator==(const MolIdentifier &other) const;
    bool operator!=(const MolIdentifier &other) const;
    
    bool operator==(const MolID &other) const;
    bool operator!=(const MolID &other) const;
    
    QList<MolNum> map(const Molecules &molecules) const;
    QList<MolNum> map(const MoleculeGroup &molgroup) const;
    QList<MolNum> map(const MolGroupsBase &molgroups) const;

private:
    /** Pointer to the MolID */
    boost::shared_ptr<MolID> d;
};

SIRE_ALWAYS_INLINE uint qHash(const MolIdentifier &molid)
{
    return molid.hash();
}

}

Q_DECLARE_METATYPE(SireMol::MolIdentifier);

#endif
