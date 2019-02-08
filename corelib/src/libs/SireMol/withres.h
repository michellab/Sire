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

#ifndef SIREMOL_WITHRES_H
#define SIREMOL_WITHRES_H

#include "resid.h"
#include "chainid.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class ChainsWithRes;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::ChainsWithRes&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::ChainsWithRes&);

namespace SireMol
{

/** This ID class identifies chains that contain residues that
    match the passed ResID
    
    @author Christopher Woods
*/
class SIREMOL_EXPORT ChainsWithRes : public ChainID
{

friend SIREMOL_EXPORT QDataStream& ::operator<<(QDataStream&, const ChainsWithRes&);
friend SIREMOL_EXPORT QDataStream& ::operator>>(QDataStream&, ChainsWithRes&);

public:
    ChainsWithRes();
    ChainsWithRes(const ResID &resid);
    ChainsWithRes(const ChainsWithRes &other);
    
    ~ChainsWithRes();
    
    static const char* typeName();
    
    const char* what() const
    {
        return ChainsWithRes::typeName();
    }
    
    ChainsWithRes* clone() const;
    
    bool isNull() const;
    
    uint hash() const;
                
    QString toString() const;
    
    const ResID& resID() const;
    
    ChainsWithRes& operator=(const ChainsWithRes &other);
    
    bool operator==(const SireID::ID &other) const;
    using SireID::ID::operator!=;
   
    bool operator==(const ChainsWithRes &other) const;
    bool operator!=(const ChainsWithRes &other) const;
    
    QList<ChainIdx> map(const MolInfo &molinfo) const;

private:
    ResIdentifier resid;
};

}

Q_DECLARE_METATYPE( SireMol::ChainsWithRes )

SIRE_EXPOSE_CLASS( SireMol::ChainsWithRes )

SIRE_END_HEADER

#endif
