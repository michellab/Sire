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

#ifndef SIREMOL_MOLATOMID_H
#define SIREMOL_MOLATOMID_H

#include "atomidentifier.h"
#include "molidentifier.h"

#include <boost/tuple/tuple.hpp>

SIRE_BEGIN_HEADER

namespace SireMol
{
class MolAtomID;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::MolAtomID&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::MolAtomID&);

namespace SireMol
{

/** This class represents an ID that is used to identify
    a specific atom (or atoms) in a specific molecule
    (of molecules)
    
    @author Christopher Woods
*/
class SIREMOL_EXPORT MolAtomID : public AtomID
{

friend QDataStream& ::operator<<(QDataStream&, const MolAtomID&);
friend QDataStream& ::operator>>(QDataStream&, MolAtomID&);

public:
    MolAtomID();
    MolAtomID(const MolID &molid, const AtomID &atomid);
    MolAtomID(const AtomID &atomid, const MolID &molid);
    
    MolAtomID(const boost::tuple<MolIdentifier,AtomIdentifier> &molatomid);
    MolAtomID(const boost::tuple<AtomIdentifier,MolIdentifier> &molatomid);
    
    MolAtomID(const MolAtomID &other);
    
    ~MolAtomID();
    
    static const char* typeName();
    
    const char* what() const
    {
        return MolAtomID::typeName();
    }
    
    MolAtomID* clone() const;
    
    bool isNull() const;
    
    uint hash() const;
                
    QString toString() const;
    
    const AtomID& atomID() const;
    const MolID& molID() const;
    
    MolAtomID& operator=(const MolAtomID &other);
    
    bool operator==(const SireID::ID &other) const;
    using SireID::ID::operator!=;
   
    bool operator==(const MolAtomID &other) const;
    bool operator!=(const MolAtomID &other) const;
    
    QList<AtomIdx> map(const MolInfo &molinfo) const;

    QHash< MolNum,Selector<Atom> > selectAllFrom(const Molecules &molecules,
                                                 const PropertyMap &map = PropertyMap()) const;
    QHash< MolNum,Selector<Atom> > selectAllFrom(const MoleculeGroup &molgroup,
                                                 const PropertyMap &map = PropertyMap()) const;
    QHash< MolNum,Selector<Atom> > selectAllFrom(const MolGroupsBase &molgroups,
                                                 const PropertyMap &map = PropertyMap()) const;

private:
    void collapse();

    /** The molecule identifier */
    MolIdentifier molid;
    
    /** The atom identifier */
    AtomIdentifier atomid;
};

}

Q_DECLARE_METATYPE( SireMol::MolAtomID )

SIRE_EXPOSE_CLASS( SireMol::MolAtomID )

SIRE_END_HEADER

#endif
