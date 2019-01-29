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

#include "specifymol.h"

#include "molecules.h"
#include "moleculegroup.h"
#include "moleculegroups.h"

#include "SireMol/errors.h"
#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMol;
using namespace SireStream;

static const RegisterMetaType<SpecifyMol> r_specifymol;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const SpecifyMol &specifymol)
{
    writeHeader(ds, r_specifymol, 1);
    
    SharedDataStream sds(ds);
    sds << specifymol.molid << specifymol.strt << specifymol.end;
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, SpecifyMol &specifymol)
{
    VersionID v = readHeader(ds, r_specifymol);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> specifymol.molid >> specifymol.strt >> specifymol.end;
    }
    else
        throw version_error( v, "1", r_specifymol, CODELOC );
        
    return ds;
}

/** Constructor */
SpecifyMol::SpecifyMol() : MolID(), strt(0), end(-1)
{}

/** Construct to match all of the molecules that match the 
    ID 'molid' */
SpecifyMol::SpecifyMol(const MolID &mol_id)
           : MolID(), molid(mol_id), strt(0), end(-1)
{}

/** Construct to match the ith molecule that matches the ID 'molid' */
SpecifyMol::SpecifyMol(const MolID &mol_id, int i)
           : MolID(), molid(mol_id), strt(i), end(i)
{}

/** Construct to math the range of molecules from i to j that
    match the ID 'molid' */
SpecifyMol::SpecifyMol(const MolID &mol_id, int i, int j)
           : MolID(), molid(mol_id), strt(i), end(j)
{}

/** Copy constructor */
SpecifyMol::SpecifyMol(const SpecifyMol &other)
           : MolID(other), molid(other.molid),
             strt(other.strt), end(other.end)
{}

/** Destructor */
SpecifyMol::~SpecifyMol()
{}

/** Copy assignment operator */
SpecifyMol& SpecifyMol::operator=(const SpecifyMol &other)
{
    molid = other.molid;
    strt = other.strt;
    end = other.end;
    MolID::operator=(other);
    
    return *this;
}

/** Comparison operator */
bool SpecifyMol::operator==(const SpecifyMol &other) const
{
    return molid == other.molid and strt == other.strt
                  and end == other.end;
}

/** Comparison operator */
bool SpecifyMol::operator!=(const SpecifyMol &other) const
{
    return molid != other.molid or strt != other.strt
                       or end != other.end;
}

/** Comparison operator */
bool SpecifyMol::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare<SpecifyMol>(*this, other);
}

/** Comparison operator */
bool SpecifyMol::operator!=(const SireID::ID &other) const
{
    return not this->operator==(other);
}

/** Hash this ID */
uint SpecifyMol::hash() const
{
    return molid.hash() + strt + end;
}

/** Return a string representation of this ID */
QString SpecifyMol::toString() const
{
    if (strt == end)
        return QString("(%1)[%2]").arg(molid.toString()).arg(strt);
    else
        return QString("(%1)[%2:%3]").arg(molid.toString())
                                     .arg(strt).arg(end);
}

/** Return whether or not this ID is null */
bool SpecifyMol::isNull() const
{
    return molid.isNull() and strt == 0 and end == -1;
}

static QList<MolNum> get(const QList<MolNum> &molnums,
                         const Index &strt, const Index &end)
{
    //now get the specified matches
    int nmatches = molnums.count();

    int sane_strt = strt.map(nmatches);
    int sane_end = end.map(nmatches);
    
    if (sane_strt > sane_end)
        qSwap(sane_strt,sane_end);
    
    if (sane_end - sane_strt == nmatches)
        return molnums;
    else
    {
        QList<MolNum> specified_molnums;
    
        for (int i=strt; i<=end; ++i)
        {
            specified_molnums.append(molnums[i]);
        }
        
        return specified_molnums;
    }
}

/** Map this ID to the list of molecule numbers that match */
QList<MolNum> SpecifyMol::map(const Molecules &molecules) const
{
    QList<MolNum> molnums = molid.map(molecules);
    return ::get(molnums, strt, end);
}

/** Map this ID to the list of molecule numbers that match */
QList<MolNum> SpecifyMol::map(const MoleculeGroup &molgroup) const
{
    QList<MolNum> molnums = molid.map(molgroup);
    return ::get(molnums, strt, end);
}

/** Map this ID to the list of molecule numbers that match */
QList<MolNum> SpecifyMol::map(const MolGroupsBase &molgroups) const
{
    QList<MolNum> molnums = molid.map(molgroups);
    return ::get(molnums, strt, end);
}

const char* SpecifyMol::typeName()
{
    return QMetaType::typeName( qMetaTypeId<SpecifyMol>() );
}

SpecifyMol* SpecifyMol::clone() const
{
    return new SpecifyMol(*this);
}
