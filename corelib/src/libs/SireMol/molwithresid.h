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

#ifndef SIREMOL_MOLWITHRESID_H
#define SIREMOL_MOLWITHRESID_H

#include "molid.h"
#include "residentifier.h"
#include "resname.h"
#include "resnum.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class MolWithResID;
}

QDataStream& operator<<(QDataStream&, const SireMol::MolWithResID&);
QDataStream& operator>>(QDataStream&, SireMol::MolWithResID&);

namespace SireMol
{

/** This class is used to identify a molecule according to 
    whether or not it contains a residue with matching 
    residue ID
    
    @author Christopher Woods
*/
class SIREMOL_EXPORT MolWithResID : public MolID
{

friend QDataStream& ::operator<<(QDataStream&, const MolWithResID&);
friend QDataStream& ::operator>>(QDataStream&, MolWithResID&);

public:
    MolWithResID();
    explicit MolWithResID(const QString &resname);
    explicit MolWithResID(int resnum);
    explicit MolWithResID(const ResID &resid);
    
    MolWithResID(const QString &resname, SireID::CaseSensitivity case_sensitivity);
    
    MolWithResID(const MolWithResID &other);
    
    ~MolWithResID();
    
    static const char* typeName();
    
    const char* what() const
    {
        return MolWithResID::typeName();
    }
    
    MolWithResID* clone() const;
    
    bool isNull() const;
    
    uint hash() const;
    
    QString toString() const;
    
    MolWithResID& operator=(const MolWithResID &other);
    
    bool operator==(const SireID::ID &other) const;
    
    bool operator==(const MolWithResID &other) const;
    
    bool operator!=(const MolWithResID &other) const;
    
    const ResID& resID() const;
    
    QList<MolNum> map(const Molecules &molecules) const;
    QList<MolNum> map(const MoleculeGroup &molgroup) const;
    QList<MolNum> map(const MolGroupsBase &molgroups) const;

private:
    /** The ResID of the residue that any matching molecule
        must contain */
    ResIdentifier resid;
};

}

Q_DECLARE_METATYPE(SireMol::MolWithResID);

SIRE_EXPOSE_CLASS( SireMol::MolWithResID )

SIRE_END_HEADER

#endif

