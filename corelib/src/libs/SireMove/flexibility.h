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

#ifndef SIREMOVE_FLEXIBILITY_H
#define SIREMOVE_FLEXIBILITY_H

#include "SireBase/propertymap.h"
#include "SireBase/shareddatapointer.hpp"

#include "SireMol/atomidx.h"
#include "SireMol/molviewproperty.h"
#include "SireMol/mover.hpp"

#include "SireUnits/units.h"

SIRE_BEGIN_HEADER

namespace SireMove
{
class Flexibility;
class DofID;
}

SIREMOVE_EXPORT QDataStream& operator<<(QDataStream &ds, const SireMove::Flexibility&);
SIREMOVE_EXPORT QDataStream& operator>>(QDataStream &ds, SireMove::Flexibility&);

SIREMOVE_EXPORT QDataStream& operator<<(QDataStream&, const SireMove::DofID&);
SIREMOVE_EXPORT QDataStream& operator>>(QDataStream&, SireMove::DofID&);

namespace SireMol
{
class AtomID;
class AtomIdx;
class BondID;
class AngleID;
class DihedralID;
class Molecule;
class MoleculeInfoData;
class MoleculeData;
}


namespace SireMove
{

using SireBase::PropertyMap;
using SireUnits::Dimension::Length;
using SireUnits::Dimension::Angle;
using SireMol::AtomID;
using SireMol::AtomIdx;
using SireMol::BondID;
using SireMol::AngleID;
using SireMol::DihedralID;
using SireMol::Molecule;
using SireMol::MoleculeData;
using SireMol::MoleculeInfoData;


/** This class implements an unique label of degrees of freedom based 
    on atomic indices. It is based on the IDQuad class from SireMM/fouratomfunctions.h
    
    @author Julien Michel
*/
class SIREMOVE_EXPORT DofID
{

friend QDataStream& ::operator<<(QDataStream&, const DofID&);
friend QDataStream& ::operator>>(QDataStream&, DofID&);

public:
    DofID();

    DofID(const AtomIdx &atom0, const AtomIdx &atom1);

    DofID(const AtomIdx &atom0, const AtomIdx &atom1, 
          const AtomIdx &atom2);

    DofID(const AtomIdx &atom0, const AtomIdx &atom1, 
          const AtomIdx &atom2, const AtomIdx &atom3);

    DofID(const DofID &other);
    
    ~DofID();
    
    DofID& operator=(const DofID &other);
    
    bool operator==(const DofID &other) const;
    bool operator!=(const DofID &other) const;
    
    static const char* typeName();

    AtomIdx atom0() const;
    AtomIdx atom1() const;
    AtomIdx atom2() const;
    AtomIdx atom3() const;

    bool isNull() const;
    bool isBond() const;
    bool isAngle() const;
    bool isDihedral() const;

private:
    void sort();

    qint32 idx0;
    qint32 idx1;
    qint32 idx2;
    qint32 idx3; 

};

/** This class holds a the list of bonds, angles and dihedrals of a molecule that 
    can be moved by a MoverMove object, as well as the maximum rotations and 
    translations that are applied to a molecule by a rigid body move object
    
    @author Julien Michel
*/
class SIREMOVE_EXPORT Flexibility
    : public SireBase::ConcreteProperty<Flexibility,SireMol::MoleculeProperty>
{
    friend QDataStream& ::operator<<(QDataStream&, const Flexibility&);
    friend QDataStream& ::operator>>(QDataStream&, Flexibility&);
      
public:
    Flexibility();
    Flexibility(const MoleculeData &molecule);
    Flexibility(const Flexibility &other);

    ~Flexibility();
    
    static const char* typeName();
      
    Flexibility& operator=(const Flexibility &other);
    
    bool operator==(const Flexibility &other) const;
    bool operator!=(const Flexibility &other) const;
      
    const SireMol::MoleculeInfoData& info() const;

    QString toString() const;

    bool isCompatibleWith(const SireMol::MoleculeInfoData &molinfo) const;

    void setRotation(const Angle &rotation);
    void setTranslation(const Length &translation);
    void setMaximumBondVar(int maxvar);
    void setMaximumAngleVar(int maxvar);
    void setMaximumDihedralVar(int maxvar);

    Length translation() const;
    Angle rotation() const;
    int maximumBondVar() const;
    int maximumAngleVar() const;
    int maximumDihedralVar() const;   

    void add(const BondID &bond, const Length &delta);
    void add(const AngleID &angle, const Angle &delta);
    void add(const DihedralID &dihedral, const Angle &delta);
      
    void remove(const BondID &bond);
    void remove(const AngleID &angle);
    void remove(const DihedralID &dihedral);
      
    bool contains(const BondID &bond) const;
    bool contains(const AngleID &angle) const;
    bool contains(const DihedralID &dihedral) const;

    void setDelta(const BondID &bond, const Length &delta);
    void setDelta(const AngleID &angle, const Angle &delta);
    void setDelta(const DihedralID &dihedral, const Angle &delta);

    Length delta(const BondID &bond) const;
    Angle delta(const AngleID &angle) const;
    Angle delta(const DihedralID &dihedral) const;

    QList<BondID> flexibleBonds() const;
    QList<AngleID> flexibleAngles() const;
    QList<DihedralID> flexibleDihedrals() const;

private:
    /** The molecule that this flexibility operates on */
    SireBase::SharedDataPointer<SireMol::MoleculeInfoData> molinfo;

    /** The maximum translation for that molecule */
    Length maxtranslation;

    /** The maximum rotation for that molecule */
    Angle maxrotation;

    /** The maximum number of bond dofs to sample in one move */
    qint32 maxbondvar;

    /** The maximum number of angle dofs to sample in one move */
    qint32 maxanglevar;

    /** The maximum number of dihedral dofs to sample in one move */
    qint32 maxdihedralvar;
    
    /** The list of delta values for bonds*/
    QHash<DofID,Length> bond_deltas;

    /** The list of delta values for angle/dihedrals*/
    QHash<DofID,Angle> angle_deltas;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Return the index of the first atom */
inline AtomIdx DofID::atom0() const
{
    return AtomIdx(idx0);
}

/** Return the index of the second atom */
inline AtomIdx DofID::atom1() const
{
    return AtomIdx(idx1);
}

/** Return the index of the third atom */
inline AtomIdx DofID::atom2() const
{
    return AtomIdx(idx2);
}

/** Return the index of the fourth atom */
inline AtomIdx DofID::atom3() const
{
    return AtomIdx(idx3);
}

inline uint qHash(const DofID &dofid)
{
    return (dofid.atom0().value() << 24) |
           ( ( dofid.atom1().value() << 16) & 0x00FF0000) |
           ( ( dofid.atom2().value() << 8)  & 0x0000FF00) |
           ( dofid.atom3().value() & 0x000000FF); 
}

#endif // SIRE_SKIP_INLINE_FUNCTIONS

} // end of namespace SireMove

Q_DECLARE_METATYPE( SireMove::DofID )
Q_DECLARE_METATYPE( SireMove::Flexibility )

SIRE_EXPOSE_CLASS( SireMove::DofID )
SIRE_EXPOSE_CLASS( SireMove::Flexibility )

SIRE_END_HEADER

#endif
