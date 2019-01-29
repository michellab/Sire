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

#ifndef SIREMOL_BEADING_H
#define SIREMOL_BEADING_H

#include "molviewproperty.h"
#include "atombeads.h"

#include <boost/shared_ptr.hpp>

SIRE_BEGIN_HEADER

namespace SireMol
{
class Beading;
class MoleculeBeading;
class ResidueBeading;
class UserBeading;

class NullBeading;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::Beading&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::Beading&);

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::MoleculeBeading&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::MoleculeBeading&);

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::ResidueBeading&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::ResidueBeading&);

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::UserBeading&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::UserBeading&);

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::NullBeading&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::NullBeading&);

namespace SireID
{
class Index;
}

namespace SireBase
{
class PropertyName;
}

namespace SireMol
{

class AtomSelection;
class AtomIdx;
class Bead;
class Beads;
class BeadIdx;
class BeadProp;
class MoleculeInfoData;
class MoleculeData;

/** This is the virtual base class of the all beading properties.
    These are used to divide a molecule into beads, and are the 
    key classes used by the SireMol::Bead and SireMol::Beads
    molecule views
    
    @author Christopher Woods
*/
class SIREMOL_EXPORT Beading : public MolViewProperty
{

friend QDataStream& ::operator<<(QDataStream&, const Beading&);
friend QDataStream& ::operator>>(QDataStream&, Beading&);

public:
    Beading();
    Beading(const Beading &other);
    
    virtual ~Beading();
    
    virtual Beading* clone() const=0;
    
    static const char* typeName();
    
    static NullBeading null();
    
    bool isCompatibleWith(const MoleculeInfoData &molinfo) const;
    
protected:
    Beading& operator=(const Beading &other);
    bool operator==(const Beading &other) const;
    bool operator!=(const Beading &other) const;
    
    ////////////////////////////////////////////////////
    // Internal functions used only by Bead and Beads //
    // and BeadProp                                   //
    ////////////////////////////////////////////////////
    friend class Bead;
    friend class Beads;
    friend class BeadProp;

    virtual int nBeads(const MoleculeInfoData &moldata) const=0;

    virtual BeadNum beadNum(const MoleculeInfoData &moldata,
                            const BeadIdx &bead) const;

    virtual AtomIdx atomIdx(const MoleculeInfoData &moldata,
                            const BeadIdx &bead, int i) const=0;

    virtual SireBase::PropertyPtr atomProperty(const MoleculeData &moldata,
                                               const SireBase::PropertyName &key) const=0;
                             
    virtual AtomSelection selection(const MoleculeInfoData &moldata) const=0;
                             
    virtual AtomSelection selection(const MoleculeInfoData &moldata,
                                    const BeadIdx &bead) const=0;

    virtual QList<AtomIdx> atomIdxs(const MoleculeInfoData &moldata) const=0;
    virtual QList<AtomIdx> atomIdxs(const MoleculeInfoData &moldata,
                                    const BeadIdx &bead) const=0;

    void assertValidIndex(BeadIdx beadidx, const MoleculeInfoData &molinfo) const;
};

/** MoleculeBeading is used to create a single bead that contains
    all of the atoms of the molecule. This is typically used for 
    small molecules, such as small solvents, where multiple beads
    per molecule would not be sensible
    
    @author Christopher Woods
*/
class SIREMOL_EXPORT MoleculeBeading
        : public SireBase::ConcreteProperty<MoleculeBeading,Beading>
{

friend QDataStream& ::operator<<(QDataStream&, const MoleculeBeading&);
friend QDataStream& ::operator>>(QDataStream&, MoleculeBeading&);

public:
    MoleculeBeading();
    MoleculeBeading(const MoleculeBeading &other);
    
    ~MoleculeBeading();
    
    MoleculeBeading& operator=(const MoleculeBeading &other);
    
    bool operator==(const MoleculeBeading &other) const;
    bool operator!=(const MoleculeBeading &other) const;
    
    static const char* typeName();
    
protected:
    int nBeads(const MoleculeInfoData &moldata) const;

    AtomIdx atomIdx(const MoleculeInfoData &moldata,
                    const BeadIdx &bead, int i) const;

    SireBase::PropertyPtr atomProperty(const MoleculeData &moldata,
                                       const SireBase::PropertyName &key) const;
                             
    AtomSelection selection(const MoleculeInfoData &moldata) const;
                             
    AtomSelection selection(const MoleculeInfoData &moldata,
                            const BeadIdx &bead) const;

    QList<AtomIdx> atomIdxs(const MoleculeInfoData &moldata) const;
    QList<AtomIdx> atomIdxs(const MoleculeInfoData &moldata,
                            const BeadIdx &bead) const;
};

/** This is a beading function that breaks a molecule into beads
    with one bead per residue
    
    @author Christopher Woods
*/
class SIREMOL_EXPORT ResidueBeading
        : public SireBase::ConcreteProperty<ResidueBeading,Beading>
{

friend QDataStream& ::operator<<(QDataStream&, const ResidueBeading&);
friend QDataStream& ::operator>>(QDataStream&, ResidueBeading&);

public:
    ResidueBeading();
    ResidueBeading(const ResidueBeading &other);
    
    ~ResidueBeading();
    
    ResidueBeading& operator=(const ResidueBeading &other);
    
    bool operator==(const ResidueBeading &other) const;
    bool operator!=(const ResidueBeading &other) const;
    
    static const char* typeName();
    
protected:
    int nBeads(const MoleculeInfoData &moldata) const;

    BeadNum beadNum(const MoleculeInfoData &moldata,
                    const BeadIdx &bead) const;

    AtomIdx atomIdx(const MoleculeInfoData &moldata,
                    const BeadIdx &bead, int i) const;

    SireBase::PropertyPtr atomProperty(const MoleculeData &moldata,
                                       const SireBase::PropertyName &key) const;
                             
    AtomSelection selection(const MoleculeInfoData &moldata) const;
                             
    AtomSelection selection(const MoleculeInfoData &moldata,
                            const BeadIdx &bead) const;

    QList<AtomIdx> atomIdxs(const MoleculeInfoData &moldata) const;
    QList<AtomIdx> atomIdxs(const MoleculeInfoData &moldata,
                            const BeadIdx &bead) const;
};

namespace detail
{
class UserBeadingInfo;
class UserBeadingInfoRegistry;
}

/** This is a beading function that divides a molecule into 
    user-defined beads
    
    @author Christopher Woods
*/
class SIREMOL_EXPORT UserBeading
        : public SireBase::ConcreteProperty<UserBeading,Beading>
{

friend QDataStream& ::operator<<(QDataStream&, const UserBeading&);
friend QDataStream& ::operator>>(QDataStream&, UserBeading&);

public:
    UserBeading();
    UserBeading(const AtomBeads &beads);
    UserBeading(const UserBeading &other);
    
    ~UserBeading();
    
    UserBeading& operator=(const UserBeading &other);
    
    bool operator==(const UserBeading &other) const;
    bool operator!=(const UserBeading &other) const;
    
    static const char* typeName();
    
    const AtomBeads& atomBeads() const;

    bool isCompatibleWith(const MoleculeInfoData &molinfo) const;
    
protected:
    int nBeads(const MoleculeInfoData &moldata) const;

    BeadNum beadNum(const MoleculeInfoData &moldata,
                    const BeadIdx &bead) const;

    AtomIdx atomIdx(const MoleculeInfoData &moldata,
                    const BeadIdx &bead, int i) const;

    SireBase::PropertyPtr atomProperty(const MoleculeData &moldata,
                                       const SireBase::PropertyName &key) const;
                             
    AtomSelection selection(const MoleculeInfoData &moldata) const;
                             
    AtomSelection selection(const MoleculeInfoData &moldata,
                            const BeadIdx &bead) const;

    QList<AtomIdx> atomIdxs(const MoleculeInfoData &moldata) const;
    QList<AtomIdx> atomIdxs(const MoleculeInfoData &moldata,
                            const BeadIdx &bead) const;

private:
    const detail::UserBeadingInfo& getUserBeadingInfo(
                                        const MoleculeInfoData &moldata) const;
    
    /** Shared pointer to the UserBeadingInfo registry */
    boost::shared_ptr<detail::UserBeadingInfoRegistry> registry;
};

/** Null beading function */
class SIREMOL_EXPORT NullBeading
        : public SireBase::ConcreteProperty<NullBeading,Beading>
{

friend QDataStream& ::operator<<(QDataStream&, const NullBeading&);
friend QDataStream& ::operator>>(QDataStream&, NullBeading&);

public:
    NullBeading();
    NullBeading(const NullBeading &other);
    
    ~NullBeading();
    
    NullBeading& operator=(const NullBeading &other);
    
    bool operator==(const NullBeading &other) const;
    bool operator!=(const NullBeading &other) const;
    
    static const char* typeName();
    
protected:
    int nBeads(const MoleculeInfoData &moldata) const;

    AtomIdx atomIdx(const MoleculeInfoData &moldata,
                    const BeadIdx &bead, int i) const;

    SireBase::PropertyPtr atomProperty(const MoleculeData &moldata,
                                       const SireBase::PropertyName &key) const;
                             
    AtomSelection selection(const MoleculeInfoData &moldata) const;
                             
    AtomSelection selection(const MoleculeInfoData &moldata,
                            const BeadIdx &bead) const;

    QList<AtomIdx> atomIdxs(const MoleculeInfoData &moldata) const;
    QList<AtomIdx> atomIdxs(const MoleculeInfoData &moldata,
                            const BeadIdx &bead) const;
};


typedef SireBase::PropPtr<Beading> BeadingPtr;

} // end of namespace SireMol

Q_DECLARE_METATYPE( SireMol::MoleculeBeading )
Q_DECLARE_METATYPE( SireMol::ResidueBeading )
Q_DECLARE_METATYPE( SireMol::UserBeading )
Q_DECLARE_METATYPE( SireMol::NullBeading )

SIRE_EXPOSE_CLASS( SireMol::Beading )
SIRE_EXPOSE_CLASS( SireMol::MoleculeBeading )
SIRE_EXPOSE_CLASS( SireMol::ResidueBeading )
SIRE_EXPOSE_CLASS( SireMol::UserBeading )
SIRE_EXPOSE_CLASS( SireMol::NullBeading )

SIRE_EXPOSE_PROPERTY( SireMol::BeadingPtr, SireMol::Beading )

SIRE_END_HEADER

#endif
