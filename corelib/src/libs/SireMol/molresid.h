/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2014  Christopher Woods
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

#ifndef SIREMOL_MOLRESID_H
#define SIREMOL_MOLRESID_H

#include "residentifier.h"
#include "molidentifier.h"

#include "molnum.h"
#include "resnum.h"

#include <boost/tuple/tuple.hpp>

SIRE_BEGIN_HEADER

namespace SireMol
{
class MolResID;
class MolResNum;
}

QDataStream& operator<<(QDataStream&, const SireMol::MolResID&);
QDataStream& operator>>(QDataStream&, SireMol::MolResID&);

QDataStream& operator<<(QDataStream&, const SireMol::MolResNum&);
QDataStream& operator>>(QDataStream&, SireMol::MolResNum&);

namespace SireMol
{

/** This class represents an ID that is used to identify
    a specific residue (or residues) in a specific molecule
    (of molecules)
    
    @author Christopher Woods
*/
class SIREMOL_EXPORT MolResID : public ResID
{

friend QDataStream& ::operator<<(QDataStream&, const MolResID&);
friend QDataStream& ::operator>>(QDataStream&, MolResID&);

public:
    MolResID();
    MolResID(const MolID &molid, const ResID &resid);
    MolResID(const ResID &resid, const MolID &molid);
    
    MolResID(const boost::tuple<MolIdentifier,ResIdentifier> &molresid);
    MolResID(const boost::tuple<ResIdentifier,MolIdentifier> &molresid);
    
    MolResID(const MolResNum &molresnum);
    
    MolResID(const MolResID &other);
    
    ~MolResID();
    
    static const char* typeName();
    
    const char* what() const
    {
        return MolResID::typeName();
    }
    
    MolResID* clone() const;
    
    bool isNull() const;
    
    uint hash() const;
                
    QString toString() const;
    
    const ResID& resID() const;
    const MolID& molID() const;
    
    MolResID& operator=(const MolResID &other);
    
    bool operator==(const SireID::ID &other) const;
    using SireID::ID::operator!=;
   
    bool operator==(const MolResID &other) const;
    bool operator!=(const MolResID &other) const;
    
    QList<ResIdx> map(const MolInfo &molinfo) const;

    QHash< MolNum,Selector<Residue> > selectAllFrom(const Molecules &molecules,
                                                    const PropertyMap &map = PropertyMap()) const;
    QHash< MolNum,Selector<Residue> > selectAllFrom(const MoleculeGroup &molgroup,
                                                    const PropertyMap &map = PropertyMap()) const;
    QHash< MolNum,Selector<Residue> > selectAllFrom(const MolGroupsBase &molgroups,
                                                    const PropertyMap &map = PropertyMap()) const;

private:
    void collapse();

    /** The molecule identifier */
    MolIdentifier molid;
    
    /** The residue identifier */
    ResIdentifier resid;
};

/** This class represents an ID that is used to identify
    a specific residue using both the residue and 
    molecule number
    
    @author Christopher Woods
*/
class SIREMOL_EXPORT MolResNum : public ResID
{

friend QDataStream& ::operator<<(QDataStream&, const MolResNum&);
friend QDataStream& ::operator>>(QDataStream&, MolResNum&);

public:
    MolResNum();
    MolResNum(const MolNum &molnum, const ResNum &resnum);
    MolResNum(const ResNum &resnum, const MolNum &molnum);
    
    MolResNum(const MolResNum &other);
    
    ~MolResNum();
    
    static const char* typeName();
    
    const char* what() const
    {
        return MolResNum::typeName();
    }
    
    MolResNum* clone() const;
    
    bool isNull() const;
    
    uint hash() const;
                
    QString toString() const;
    
    const ResNum& resNum() const;
    const MolNum& molNum() const;
    
    MolResNum& operator=(const MolResNum &other);
    
    bool operator==(const SireID::ID &other) const;
    using SireID::ID::operator!=;
   
    bool operator==(const MolResNum &other) const;
    bool operator!=(const MolResNum &other) const;
    
    QList<ResIdx> map(const MolInfo &molinfo) const;

    QHash< MolNum,Selector<Residue> > selectAllFrom(const Molecules &molecules,
                                                    const PropertyMap &map = PropertyMap()) const;
    QHash< MolNum,Selector<Residue> > selectAllFrom(const MoleculeGroup &molgroup,
                                                    const PropertyMap &map = PropertyMap()) const;
    QHash< MolNum,Selector<Residue> > selectAllFrom(const MolGroupsBase &molgroups,
                                                    const PropertyMap &map = PropertyMap()) const;

private:
    /** The molecule number */
    MolNum molnum;
    
    /** The residue numner */
    ResNum resnum;
};

}

Q_DECLARE_METATYPE( SireMol::MolResID )
Q_DECLARE_METATYPE( SireMol::MolResNum )

SIRE_EXPOSE_CLASS( SireMol::MolResID )
SIRE_EXPOSE_CLASS( SireMol::MolResNum )

SIRE_END_HEADER

#endif
