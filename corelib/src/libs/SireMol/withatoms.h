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

#ifndef SIREMOL_WITHATOMS_H
#define SIREMOL_WITHATOMS_H

#include "atomid.h"
#include "resid.h"
#include "cgid.h"
#include "chainid.h"
#include "segid.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class ResWithAtoms;
class CGsWithAtoms;
class ChainsWithAtoms;
class SegsWithAtoms;
}

QDataStream& operator<<(QDataStream&, const SireMol::ResWithAtoms&);
QDataStream& operator>>(QDataStream&, SireMol::ResWithAtoms&);

QDataStream& operator<<(QDataStream&, const SireMol::CGsWithAtoms&);
QDataStream& operator>>(QDataStream&, SireMol::CGsWithAtoms&);

QDataStream& operator<<(QDataStream&, const SireMol::ChainsWithAtoms&);
QDataStream& operator>>(QDataStream&, SireMol::ChainsWithAtoms&);

QDataStream& operator<<(QDataStream&, const SireMol::SegsWithAtoms&);
QDataStream& operator>>(QDataStream&, SireMol::SegsWithAtoms&);

namespace SireMol
{

/** This ID class identifies residues that contain atoms that
    match the passed AtomID
    
    @author Christopher Woods
*/
class SIREMOL_EXPORT ResWithAtoms : public ResID
{

friend QDataStream& ::operator<<(QDataStream&, const ResWithAtoms&);
friend QDataStream& ::operator>>(QDataStream&, ResWithAtoms&);

public:
    ResWithAtoms();
    ResWithAtoms(const AtomID &atomid);
    ResWithAtoms(const ResWithAtoms &other);
    
    ~ResWithAtoms();
    
    static const char* typeName();
    
    const char* what() const
    {
        return ResWithAtoms::typeName();
    }
    
    ResWithAtoms* clone() const;
    
    bool isNull() const;
    
    uint hash() const;
                
    QString toString() const;
    
    const AtomID& atomID() const;
    
    ResWithAtoms& operator=(const ResWithAtoms &other);
    
    bool operator==(const SireID::ID &other) const;
    using SireID::ID::operator!=;
   
    bool operator==(const ResWithAtoms &other) const;
    bool operator!=(const ResWithAtoms &other) const;
    
    QList<ResIdx> map(const MolInfo &molinfo) const;

private:
    AtomIdentifier atomid;
};

/** This ID class identifies CutGroups that contain atoms that
    match the passed AtomID
    
    @author Christopher Woods
*/
class SIREMOL_EXPORT CGsWithAtoms : public CGID
{

friend QDataStream& ::operator<<(QDataStream&, const CGsWithAtoms&);
friend QDataStream& ::operator>>(QDataStream&, CGsWithAtoms&);

public:
    CGsWithAtoms();
    CGsWithAtoms(const AtomID &atomid);
    CGsWithAtoms(const CGsWithAtoms &other);
    
    ~CGsWithAtoms();
    
    static const char* typeName();
    
    const char* what() const
    {
        return CGsWithAtoms::typeName();
    }
    
    CGsWithAtoms* clone() const;
    
    bool isNull() const;
    
    uint hash() const;
                
    QString toString() const;
    
    const AtomID& atomID() const;
    
    CGsWithAtoms& operator=(const CGsWithAtoms &other);
    
    bool operator==(const SireID::ID &other) const;
    using SireID::ID::operator!=;
   
    bool operator==(const CGsWithAtoms &other) const;
    bool operator!=(const CGsWithAtoms &other) const;
    
    QList<CGIdx> map(const MolInfo &molinfo) const;

private:
    AtomIdentifier atomid;
};

/** This ID class identifies chains that contain atoms that
    match the passed AtomID
    
    @author Christopher Woods
*/
class SIREMOL_EXPORT ChainsWithAtoms : public ChainID
{

friend QDataStream& ::operator<<(QDataStream&, const ChainsWithAtoms&);
friend QDataStream& ::operator>>(QDataStream&, ChainsWithAtoms&);

public:
    ChainsWithAtoms();
    ChainsWithAtoms(const AtomID &atomid);
    ChainsWithAtoms(const ChainsWithAtoms &other);
    
    ~ChainsWithAtoms();
    
    static const char* typeName();
    
    const char* what() const
    {
        return ChainsWithAtoms::typeName();
    }
    
    ChainsWithAtoms* clone() const;
    
    bool isNull() const;
    
    uint hash() const;
                
    QString toString() const;
    
    const AtomID& atomID() const;
    
    ChainsWithAtoms& operator=(const ChainsWithAtoms &other);
    
    bool operator==(const SireID::ID &other) const;
    using SireID::ID::operator!=;
   
    bool operator==(const ChainsWithAtoms &other) const;
    bool operator!=(const ChainsWithAtoms &other) const;
    
    QList<ChainIdx> map(const MolInfo &molinfo) const;

private:
    AtomIdentifier atomid;
};

/** This ID class identifies segments that contain atoms that
    match the passed AtomID
    
    @author Christopher Woods
*/
class SIREMOL_EXPORT SegsWithAtoms : public SegID
{

friend QDataStream& ::operator<<(QDataStream&, const SegsWithAtoms&);
friend QDataStream& ::operator>>(QDataStream&, SegsWithAtoms&);

public:
    SegsWithAtoms();
    SegsWithAtoms(const AtomID &atomid);
    SegsWithAtoms(const SegsWithAtoms &other);
    
    ~SegsWithAtoms();
    
    static const char* typeName();
    
    const char* what() const
    {
        return SegsWithAtoms::typeName();
    }
    
    SegsWithAtoms* clone() const;
    
    bool isNull() const;
    
    uint hash() const;
                
    QString toString() const;
    
    const AtomID& atomID() const;
    
    SegsWithAtoms& operator=(const SegsWithAtoms &other);
    
    bool operator==(const SireID::ID &other) const;
    using SireID::ID::operator!=;
   
    bool operator==(const SegsWithAtoms &other) const;
    bool operator!=(const SegsWithAtoms &other) const;
    
    QList<SegIdx> map(const MolInfo &molinfo) const;

private:
    AtomIdentifier atomid;
};

}

Q_DECLARE_METATYPE( SireMol::ResWithAtoms )
Q_DECLARE_METATYPE( SireMol::CGsWithAtoms )
Q_DECLARE_METATYPE( SireMol::ChainsWithAtoms )
Q_DECLARE_METATYPE( SireMol::SegsWithAtoms )

SIRE_EXPOSE_CLASS( SireMol::ResWithAtoms )
SIRE_EXPOSE_CLASS( SireMol::CGsWithAtoms )
SIRE_EXPOSE_CLASS( SireMol::ChainsWithAtoms )
SIRE_EXPOSE_CLASS( SireMol::SegsWithAtoms )

SIRE_END_HEADER

#endif
