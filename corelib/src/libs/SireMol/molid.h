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

#ifndef SIREMOL_MOLID_H
#define SIREMOL_MOLID_H

#include "SireID/id.h"

#include "SireID/idandset.hpp"
#include "SireID/idorset.hpp"

SIRE_BEGIN_HEADER

namespace SireMol
{

using SireID::IDAndSet;
using SireID::IDOrSet;

class MolIdx;
class MolIdentifier;
class MolNum;

class MolAtomID;

class SpecifyMol;

class AtomID;

class Molecules;
class MoleculeGroup;
class MolGroupsBase;

/** This is the base class of all identifiers that are used 
    to identify a Molecule

    @author Christopher Woods
*/
class SIREMOL_EXPORT MolID : public SireID::ID
{
public:
    typedef MolNum Index;
    typedef MolIdentifier Identifier;
    typedef Molecules SearchObject;

    MolID();
    
    MolID(const MolID &other);
    
    virtual ~MolID();

    static const char* typeName()
    {
        return "SireMol::MolID";
    }

    virtual MolID* clone() const=0;
    
    SpecifyMol operator[](int i) const;
    SpecifyMol operator()(int i) const;
    SpecifyMol operator()(int i, int j) const;
    
    IDAndSet<MolID> operator+(const MolID &other) const;
    MolAtomID operator+(const AtomID &other) const;

    IDOrSet<MolID> operator*(const MolID &other) const;

    IDAndSet<MolID> operator&&(const MolID &other) const;
    MolAtomID operator&&(const AtomID &other) const;

    IDAndSet<MolID> operator&(const MolID &other) const;
    MolAtomID operator&(const AtomID &other) const;
    
    IDOrSet<MolID> operator||(const MolID &other) const;
    IDOrSet<MolID> operator|(const MolID &other) const;
    
    IDOrSet<AtomID> operator*(const AtomID &other) const;
    IDOrSet<AtomID> operator|(const AtomID &other) const;
    IDOrSet<AtomID> operator||(const AtomID &other) const;
    
    /* TODO!!!
    SireID::InvertMatch<MolID> operator!() const;
    SireID::InvertMatch<MolID> invert() const;
    SireID::InvertMatch<MolID> inverse() const;
    
    static SireID::MatchAll<MolID> any();*/
    
    virtual QList<MolNum> map(const Molecules &molecules) const=0;
    virtual QList<MolNum> map(const MoleculeGroup &molgroup) const=0;
    virtual QList<MolNum> map(const MolGroupsBase &molgroupsbase) const=0;

protected:
    void processMatches(QList<MolNum> &matches, const Molecules &mols) const;
};

}

#include "molnum.h"
#include "molidentifier.h"

namespace SireID
{

using SireMol::MolID;
using SireMol::MolIdentifier;
using SireMol::MolNum;
using SireMol::Molecules;
using SireMol::MoleculeGroup;
using SireMol::MolGroupsBase;

template<>
class SIREMOL_EXPORT IDAndSet<MolID> : public MolID
{

friend QDataStream& ::operator<<<>(QDataStream&, const IDAndSet<MolID>&);
friend QDataStream& ::operator>><>(QDataStream&, IDAndSet<MolID>&);

public:
    IDAndSet();
    IDAndSet(const MolID &id);
    IDAndSet(const MolID &id0, const MolID &id1);
    
    IDAndSet(const QList<MolIdentifier> &ids);
    
    IDAndSet(const IDAndSet<MolID> &other);
    
    ~IDAndSet();
    
    static const char* typeName();
    
    const char* what() const
    {
        return IDAndSet<MolID>::typeName();
    }
    
    IDAndSet<MolID>* clone() const;
    
    bool isNull() const;
    
    uint hash() const;
                
    QString toString() const;

    const QSet<MolIdentifier>& IDs() const;
    
    IDAndSet<MolID>& operator=(const IDAndSet<MolID> &other);
    IDAndSet<MolID>& operator=(const MolID &other);
    
    bool operator==(const SireID::ID &other) const;
    bool operator!=(const SireID::ID &other) const;
   
    bool operator==(const IDAndSet<MolID> &other) const;
    bool operator!=(const IDAndSet<MolID> &other) const;
    
    bool operator==(const MolID &other) const;
    bool operator!=(const MolID &other) const;
    
    QList<MolNum> map(const Molecules &mols) const;
    QList<MolNum> map(const MoleculeGroup &molgroup) const;
    QList<MolNum> map(const MolGroupsBase &molgroups) const;

private:
    void add(const MolID &id);

    template<class T>
    QList<MolNum> _pvt_map(const T &group) const;
    
    QSet<MolIdentifier> ids;
};

template<>
class SIREMOL_EXPORT IDOrSet<MolID> : public MolID
{

friend QDataStream& ::operator<<<>(QDataStream&, const IDOrSet<MolID>&);
friend QDataStream& ::operator>><>(QDataStream&, IDOrSet<MolID>&);

public:
    IDOrSet();
    IDOrSet(const MolID &id);
    IDOrSet(const MolID &id0, const MolID &id1);
    
    IDOrSet(const QList<MolIdentifier> &ids);
    
    IDOrSet(const IDOrSet<MolID> &other);
    
    ~IDOrSet();
    
    static const char* typeName();
    
    const char* what() const
    {
        return IDOrSet<MolID>::typeName();
    }
    
    IDOrSet<MolID>* clone() const;
    
    bool isNull() const;
    
    uint hash() const;
                
    QString toString() const;

    const QSet<MolIdentifier>& IDs() const;
    
    IDOrSet<MolID>& operator=(const IDOrSet<MolID> &other);
    IDOrSet<MolID>& operator=(const MolID &other);
    
    bool operator==(const SireID::ID &other) const;
    bool operator!=(const SireID::ID &other) const;
   
    bool operator==(const IDOrSet<MolID> &other) const;
    bool operator!=(const IDOrSet<MolID> &other) const;
    
    bool operator==(const MolID &other) const;
    bool operator!=(const MolID &other) const;
    
    QList<MolNum> map(const Molecules &mols) const;
    QList<MolNum> map(const MoleculeGroup &molgroup) const;
    QList<MolNum> map(const MolGroupsBase &molgroups) const;

private:
    void add(const MolID &id);
    QList<MolNum> process(QList<MolNum> molnums) const;

    QSet<MolIdentifier> ids;
};

} // end of namespace SireID

Q_DECLARE_METATYPE( SireID::IDAndSet<SireMol::MolID> )
Q_DECLARE_METATYPE( SireID::IDOrSet<SireMol::MolID> )

#include "molidentifier.h"

SIRE_EXPOSE_CLASS( SireMol::MolID )
SIRE_EXPOSE_ALIAS( SireID::IDAndSet<SireMol::MolID>, SireMol::IDAndSet_MolID_ )
SIRE_EXPOSE_ALIAS( SireID::IDOrSet<SireMol::MolID>, SireMol::IDOrSet_MolID_ )

SIRE_END_HEADER

#endif
