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

#include "molwithresid.h"
#include "molecules.h"
#include "moleculegroup.h"
#include "moleculegroups.h"
#include "moleculeview.h"
#include "viewsofmol.h"

#include "SireMol/errors.h"

#include "SireStream/shareddatastream.h"

using namespace SireMol;
using namespace SireStream;

static const RegisterMetaType<MolWithResID> r_mol;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const MolWithResID &mol)
{
    writeHeader(ds, r_mol, 1);
    
    SharedDataStream sds(ds);
    sds << mol.resid;
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, MolWithResID &mol)
{
    VersionID v = readHeader(ds, r_mol);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> mol.resid;
    }
    else
        throw version_error( v, "1", r_mol, CODELOC );
        
    return ds;
}

/** Null constructor */
MolWithResID::MolWithResID() : MolID()
{}

/** Construct to match any molecule that contains the residue with name 'resname' */
MolWithResID::MolWithResID(const QString &resname)
             : MolID(), resid( ResName(resname) )
{}

/** Construct to match any molecule that contains the residue with name 'resname'. 
    (allowing you to specify whether the name match should be case sensitive or not) */
MolWithResID::MolWithResID(const QString &resname, SireID::CaseSensitivity case_sensitivity)
             : MolID(), resid( ResName(resname,case_sensitivity) )
{}

/** Construct to match any molecule that contains the residue with number 'resnum' */
MolWithResID::MolWithResID(int resnum)
             : MolID(), resid( ResNum(resnum) )
{}

/** Construct to match any molecule that contains the residue that matches ID 'resid' */
MolWithResID::MolWithResID(const ResID &id)
             : MolID(), resid(id)
{}

/** Copy constructor */
MolWithResID::MolWithResID(const MolWithResID &other) : MolID(other), resid(other.resid)
{}

/** Destructor */
MolWithResID::~MolWithResID()
{}

bool MolWithResID::isNull() const
{
    return resid.isNull();
}

uint MolWithResID::hash() const
{
    return resid.hash();
}

QString MolWithResID::toString() const
{
    if (isNull())
        return QObject::tr("MolWithResID::null");
    else
        return QObject::tr("MolWithResID{ residue = %1 }").arg(resid.toString());
}

MolWithResID& MolWithResID::operator=(const MolWithResID &other)
{
    resid = other.resid;
    MolID::operator=(other);
    return *this;
}

MolWithResID* MolWithResID::clone() const
{
    return new MolWithResID(*this);
}

bool MolWithResID::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare<MolWithResID>(*this, other);
}

bool MolWithResID::operator==(const MolWithResID &other) const
{
    return resid == other.resid;
}

bool MolWithResID::operator!=(const MolWithResID &other) const
{
    return not operator==(other);
}

/** Return the residue ID that is used to find the residue (and thus the molecule) */
const ResID& MolWithResID::resID() const
{
    return resid;
}

QList<MolNum> MolWithResID::map(const Molecules &molecules) const
{
    QList<MolNum> molnums;

    if (resid.isNull())
    {
        molnums = molecules.molNums().toList();
    }
    else
    {
        SireError::FastExceptionFlag f = SireError::exception::enableFastExceptions();
    
        for (Molecules::const_iterator it = molecules.constBegin();
             it != molecules.constEnd();
             ++it)
        {
            try
            {
                QList<ResIdx> residxs = resid.map(it.value().data().info());
                
                if (not it.value().selectedAll())
                {
                    AtomSelection selection = it.value().selection();
                    
                    bool is_selected = false;
                    
                    foreach (ResIdx residx, residxs)
                    {
                        if (selection.selected(residx))
                        {
                            is_selected = true;
                            break;
                        }
                    }
                    
                    if (is_selected)
                        molnums.append( it.key() );
                }
                else
                    molnums.append(it.key());
            }
            catch(...)
            {}
        }
    }
    
    if (molnums.isEmpty())
        throw SireMol::missing_molecule( QObject::tr(
            "There is no molecule with residue matching \"%1\" in the set of molecules.")
                .arg(resid.toString()), CODELOC );
                
    return molnums;
}

QList<MolNum> MolWithResID::map(const MoleculeGroup &molgroup) const
{
    return this->map(molgroup.molecules());
}

QList<MolNum> MolWithResID::map(const MolGroupsBase &molgroups) const
{
    return this->map(molgroups.molecules());
}

const char* MolWithResID::typeName()
{
    return QMetaType::typeName( qMetaTypeId<MolWithResID>() );
}

