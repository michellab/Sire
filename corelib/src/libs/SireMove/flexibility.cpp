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

#include "flexibility.h"

#include <QList>

#include "SireMol/molecule.h"
#include "SireMol/partialmolecule.h"
#include "SireMol/bondid.h"
#include "SireMol/angleid.h"
#include "SireMol/dihedralid.h"
#include "SireMol/atomidx.h"
#include "SireMol/mover.hpp"

#include "SireUnits/convert.h"
#include "SireUnits/units.h"

#include "SireMove/errors.h"
#include "SireMol/errors.h"
#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMove;
using namespace SireMol;
using namespace SireUnits;
using namespace SireUnits::Dimension;
using namespace SireStream;


//////
////// Implementation of DofID
//////

static const RegisterMetaType<DofID> r_dofid(NO_ROOT);

QDataStream &operator<<(QDataStream &ds, const DofID &dofid)
{
    writeHeader(ds, r_dofid, 2);

    SharedDataStream sds(ds);
  
    sds << dofid.idx0 << dofid.idx1
        << dofid.idx2 << dofid.idx3;
       
    return ds;
}

QDataStream &operator>>(QDataStream &ds, DofID &dofid)
{
    VersionID v = readHeader(ds, r_dofid);

    if (v == 2)
    {
        SharedDataStream sds(ds);
   
        sds >> dofid.idx0 >> dofid.idx1
            >> dofid.idx2 >> dofid.idx3;
    }
    else
        throw version_error(v, "2", r_dofid, CODELOC);

   return ds;
}

bool DofID::isBond() const
{
    return idx0 >=0 and idx1 >= 0 and idx2 < 0 and idx3 < 0;
}

bool DofID::isAngle() const
{
    return idx0 >=0 and idx1 >= 0 and idx2 >= 0 and idx3 < 0;
}

bool DofID::isDihedral() const
{
    return idx0 >=0 and idx1 >= 0 and idx2 >= 0 and idx3 >= 0;
}

bool DofID::isNull() const
{
    return idx0 < 0 and idx1 < 0 and idx2 < 0 and idx3 < 0;
}

void DofID::sort()
{
    if (isBond())
    {
        if (idx1 < idx0)
            qSwap(idx1, idx0);
    }
    else if (isAngle())
    {
        if (idx2 < idx0)
            qSwap(idx2, idx0);
    }
    else if (isDihedral())
    {
        if (idx3 < idx0)
        {
            qSwap(idx3, idx0);
            qSwap(idx2, idx1);
        }
    }
    else
    {
        idx0 = -1;
        idx1 = -1;
        idx2 = -1;
        idx3 = -1;
    }
}

/** Null constructor */
DofID::DofID() 
      : idx0(-1), idx1(-1), idx2(-1), idx3(-1)
{}

/** Constructor for a set of 2 AtomIdxs*/
DofID::DofID(const AtomIdx &atom0, const AtomIdx &atom1)
      : idx0(atom0.value()), idx1(atom1.value()),
        idx2(-1), idx3(-1)
{
    this->sort();
}

/** Constructor for a set of 3 AtomIdxs*/
DofID::DofID(const AtomIdx &atom0, const AtomIdx &atom1, const AtomIdx &atom2)
      : idx0(atom0.value()), idx1(atom1.value()),
        idx2(atom2.value()), idx3(-1)
{
    this->sort();
}

/** Constructor for a set of 4 AtomIdxs*/
DofID::DofID(const AtomIdx &atom0, const AtomIdx &atom1, 
             const AtomIdx &atom2, const AtomIdx &atom3)
      : idx0(atom0.value()), idx1(atom1.value()),
        idx2(atom2.value()), idx3(atom3.value())
{
    this->sort();
}

/** Copy constructor */
DofID::DofID(const DofID &other)
      : idx0(other.idx0), idx1(other.idx1), 
        idx2(other.idx2), idx3(other.idx3)
{}

/** Destructor */
DofID::~DofID()
{}

/** Copy assignment operator */
DofID& DofID::operator=(const DofID &other)
{
    idx0 = other.idx0;
    idx1 = other.idx1;
    idx2 = other.idx2;
    idx3 = other.idx3;
    
    return *this;
}

/** Comparison operator */
bool DofID::operator==(const DofID &other) const
{
    return idx0 == other.idx0 and idx1 == other.idx1 and
           idx2 == other.idx2 and idx3 == other.idx3;
}

/** Comparison operator */
bool DofID::operator!=(const DofID &other) const
{
    return not DofID::operator==(other);
}

const char* DofID::typeName()
{
    return QMetaType::typeName(qMetaTypeId<DofID>());
}

//
// Implementation of Flexibility
//
static const RegisterMetaType<Flexibility> r_flex;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const Flexibility &flex)
{
    writeHeader(ds, r_flex, 2);
    
    SharedDataStream sds(ds);
    
    sds << flex.molinfo << flex.maxtranslation 
        << flex.maxrotation << flex.maxbondvar 
	<< flex.maxanglevar << flex.maxdihedralvar
        << flex.bond_deltas << flex.angle_deltas 
        << static_cast<const MoleculeProperty&>(flex);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, Flexibility &flex)
{
    VersionID v = readHeader(ds, r_flex);
    
    if (v == 2)
    {
        SharedDataStream sds(ds);
        
        sds >> flex.molinfo >> flex.maxtranslation
            >> flex.maxrotation >> flex.maxbondvar 
	    >> flex.maxanglevar >> flex.maxdihedralvar
            >> flex.bond_deltas >> flex.angle_deltas
            >> static_cast<MoleculeProperty&>(flex);
    }
    else if (v == 1)
    {
        SharedDataStream sds(ds);

        qint32 maxvar;

        sds >> flex.molinfo >> flex.maxtranslation
            >> flex.maxrotation >> maxvar 
            >> flex.bond_deltas >> flex.angle_deltas
            >> static_cast<MoleculeProperty&>(flex);

        flex.maxbondvar = maxvar;
        flex.maxanglevar = maxvar;
        flex.maxdihedralvar = maxvar;
    }
    else
        throw version_error(v, "1,2", r_flex, CODELOC);
        
    return ds;
}

/** Null Constructor */
Flexibility::Flexibility() : ConcreteProperty<Flexibility,MoleculeProperty>()
{}

/** Constructor for the passed molecule*/
Flexibility::Flexibility(const MoleculeData &molecule)
            : ConcreteProperty<Flexibility,MoleculeProperty>(),
              molinfo(molecule.info())
{}

/** Copy constructor */
Flexibility::Flexibility(const Flexibility &other)
            : ConcreteProperty<Flexibility,MoleculeProperty>(),
              molinfo(other.molinfo),maxtranslation(other.maxtranslation),
              maxrotation(other.maxrotation),maxbondvar(other.maxbondvar),
	      maxanglevar(other.maxanglevar), maxdihedralvar(other.maxdihedralvar),
              bond_deltas(other.bond_deltas),angle_deltas(other.angle_deltas)
{}

/** Destructor */
Flexibility::~Flexibility()
{}

/** Copy assignment operator */
Flexibility& Flexibility::operator=(const Flexibility &other)
{
    if (this != &other)
    {
        MoleculeProperty::operator=(other);
        molinfo = other.molinfo;
        maxtranslation = other.maxtranslation;
        maxrotation = other.maxrotation;
        maxbondvar = other.maxbondvar;
	maxanglevar = other.maxanglevar;
	maxdihedralvar = other.maxdihedralvar;
        bond_deltas = other.bond_deltas;
        angle_deltas = other.angle_deltas;
    }
  
    return *this;
}

/** Comparison operator */
bool Flexibility::operator==(const Flexibility &other) const
{
    return (molinfo == other.molinfo and maxtranslation == other.maxtranslation and
	    maxrotation == other.maxrotation and maxbondvar == other.maxbondvar and
	    maxanglevar == other.maxanglevar and maxdihedralvar == other.maxdihedralvar and
	    bond_deltas == other.bond_deltas and angle_deltas == other.angle_deltas);
}

/** Comparison operator */
bool Flexibility::operator!=(const Flexibility &other) const
{
    return not Flexibility::operator==(other);
}

/** Return the layout of the molecule whose flexibility is contained
    in this object */
const MoleculeInfoData& Flexibility::info() const
{
    if (molinfo.constData() == 0)
        return MoleculeInfoData::null();
    else
        return *molinfo;
}

static QString getAtom(const AtomIdx &idx, const MoleculeInfoData &molinfo)
{
    return QString("%1:%2").arg( molinfo.name(idx) ).arg( molinfo.number(idx) );
}

/** Return a string representation of this flexibility */
QString Flexibility::toString() const
{
    QStringList lines;
    
    lines.append(QObject::tr("Flexibility ") );
    
    lines.append(QObject::tr("Rotation: %1 , Translation: %2 ")
                            .arg(maxrotation.toString())
                            .arg(maxtranslation.toString()));

    lines.append(QObject::tr("Maximum Bond Variables: %1")
                            .arg(maxbondvar));
    lines.append(QObject::tr("Maximum Angle Variables: %1")
                            .arg(maxanglevar));
    lines.append(QObject::tr("Maximum Dihedral Variables: %1")
                            .arg(maxdihedralvar));

    for (QHash<DofID,Length>::const_iterator it = bond_deltas.constBegin();
         it != bond_deltas.constEnd();
         ++it)
    {
        const DofID &bond = it.key();
    
        lines.append(QObject::tr("%2-%3 = %1")
                    .arg(it.value().toString())
                    .arg(getAtom(bond.atom0(),*molinfo), getAtom(bond.atom1(),*molinfo)));
    }
    
    for (QHash<DofID,Angle>::const_iterator it = angle_deltas.constBegin();
         it != angle_deltas.constEnd();
         ++it)
    {
        const DofID &angle = it.key();
    
        if (angle.atom3().isNull())
        {
            lines.append(QObject::tr("%1-%2-%3 = %4")
                    .arg(getAtom(angle.atom0(),*molinfo), getAtom(angle.atom1(),*molinfo))
                    .arg(getAtom(angle.atom2(),*molinfo))
                    .arg(it.value().toString()));
        }
    }
    
    for (QHash<DofID,Angle>::const_iterator it = angle_deltas.constBegin();
         it != angle_deltas.constEnd();
         ++it)
    {
        const DofID &angle = it.key();

        if (not angle.atom3().isNull())
        {
            lines.append(QObject::tr("%1-%2-%3-%4 = %5")
                    .arg(getAtom(angle.atom0(),*molinfo), getAtom(angle.atom1(),*molinfo))
                    .arg(getAtom(angle.atom2(),*molinfo), getAtom(angle.atom3(),*molinfo))
                    .arg(it.value().toString()));
        }
    }

    return lines.join("\n");
}

/** Return whether or not this flexibility is compatible with the molecule 
    whose info is in 'molinfo' */
bool Flexibility::isCompatibleWith(const SireMol::MoleculeInfoData &molinfo) const
{
    return info().UID() == molinfo.UID();
}

/** Set the maximum rotation of this flexibility */
void Flexibility::setRotation(const Angle &rotation)
{
    this->maxrotation = rotation;
}

/** Set the maximum translation of this flexibility */
void Flexibility::setTranslation(const Length &translation)
{
    this->maxtranslation = translation;
}

/** Set the maximum number of degrees of freedom that will be sampled in one move */
void Flexibility::setMaximumBondVar(int maxvar)
{
    this->maxbondvar = maxvar;
}

void Flexibility::setMaximumAngleVar(int maxvar)
{
    this->maxanglevar = maxvar;
}

void Flexibility::setMaximumDihedralVar(int maxvar)
{
    this->maxdihedralvar = maxvar;
}


/** Return the maximum rotation of this flexibility*/
Angle Flexibility::rotation() const
{
    return this->maxrotation;
}

/** Return the maximum translation of this flexibility */
Length Flexibility::translation() const
{
    return this->maxtranslation;
}

/** Return the maximum number of dofs that will be sampled in one move */
int Flexibility::maximumBondVar() const
{
    return this->maxbondvar;
}

int Flexibility::maximumAngleVar() const
{
    return this->maxanglevar;
}

int Flexibility::maximumDihedralVar() const
{
    return this->maxdihedralvar;
}


static DofID getBond(const BondID &bond, 
                     const SharedDataPointer<MoleculeInfoData> &molinfo)
{
    return DofID( molinfo->atomIdx(bond.atom0()), molinfo->atomIdx(bond.atom1()) );
}

static DofID getAngle(const AngleID &angle, 
                      const SharedDataPointer<MoleculeInfoData> &molinfo)
{
    return DofID( molinfo->atomIdx(angle.atom0()), molinfo->atomIdx(angle.atom1()),
                  molinfo->atomIdx(angle.atom2()) );
}

static DofID getDihedral(const DihedralID &dihedral, 
                         const SharedDataPointer<MoleculeInfoData> &molinfo)
{
    return DofID( molinfo->atomIdx(dihedral.atom0()),
                  molinfo->atomIdx(dihedral.atom1()),
                  molinfo->atomIdx(dihedral.atom2()), 
                  molinfo->atomIdx(dihedral.atom3()) );
}

/** Add bond with delta to this flexibility*/
void Flexibility::add(const BondID &bond, const Length &delta)
{
    bond_deltas.insert(::getBond(bond,molinfo), delta);
}

/** Add angle with delta to this flexibility*/
void Flexibility::add(const AngleID &angle, const Angle &delta)
{
    angle_deltas.insert(::getAngle(angle,molinfo), delta);
}

/** Add dihedral with delta to this flexibility*/
void Flexibility::add(const DihedralID &dihedral, const Angle &delta)
{
    angle_deltas.insert(::getDihedral(dihedral,molinfo), delta);
}

/** Remove bond from this flexibility*/
void Flexibility::remove(const BondID &bond)
{
    bond_deltas.remove(::getBond(bond,molinfo));
}

/** Remove angle from this flexibility*/
void Flexibility::remove(const AngleID &angle)
{
    angle_deltas.remove(::getAngle(angle,molinfo));
}

/** Remove dihedral from this flexibility*/
void Flexibility::remove(const DihedralID &dihedral)
{
    angle_deltas.remove(::getDihedral(dihedral,molinfo));
}

/** Check if bond is present in this flexibility */
bool Flexibility::contains(const BondID &bond) const
{
    return bond_deltas.contains(::getBond(bond,molinfo));
}

/** Check if angle is present in this flexibility */
bool Flexibility::contains(const AngleID &angle) const
{
    return angle_deltas.contains(::getAngle(angle,molinfo));
}

/** Check if angle is present in this flexibility */
bool Flexibility::contains(const DihedralID &dihedral) const
{
    return angle_deltas.contains(::getDihedral(dihedral,molinfo));
}

/** set the delta value of bond to delta*/
void Flexibility::setDelta(const BondID &bond, const Length &delta)
{
    this->add(bond, delta);
}

/** set the delta value of bond to delta*/
void Flexibility::setDelta(const AngleID &angle, const Angle &delta)
{
    this->add(angle, delta);
}

/** set the delta value of bond to delta*/
void Flexibility::setDelta(const DihedralID &dihedral, const Angle &delta)
{
    this->add(dihedral, delta);
}

/** Return the delta value of bond in this flexibility */
Length Flexibility::delta(const BondID &bond) const
{
    return bond_deltas.value( ::getBond(bond,molinfo), Length(0) );
}

/** Return the delta value of angle in this flexibility */
Angle Flexibility::delta(const AngleID &angle) const
{
    return angle_deltas.value( ::getAngle(angle,molinfo), Angle(0) );
}

/** Return the delta value of angle in this flexibility */
Angle Flexibility::delta(const DihedralID &dihedral) const
{
    return angle_deltas.value( ::getDihedral(dihedral,molinfo), Angle(0) );
}

/** Return the list of all flexible bonds */
QList<BondID> Flexibility::flexibleBonds() const
{
    QList<BondID> bonds;
    
    for (QHash<DofID,Length>::const_iterator it = bond_deltas.constBegin();
         it != bond_deltas.constEnd();
         ++it)
    {
        bonds.append( BondID(it.key().atom0(), it.key().atom1()) );
    }
    
    return bonds;
}

/** Return the list of all flexible angles */
QList<AngleID> Flexibility::flexibleAngles() const
{
    QList<AngleID> angles;
    
    for (QHash<DofID,Angle>::const_iterator it = angle_deltas.constBegin();
         it != angle_deltas.constEnd();
         ++it)
    {
        if (it.key().isAngle())
            angles.append( AngleID(it.key().atom0(), it.key().atom1(),
                                   it.key().atom2()) );
    }
    
    return angles;
}

/** Return the list of all flexible dihedrals */
QList<DihedralID> Flexibility::flexibleDihedrals() const
{
    QList<DihedralID> dihedrals;
    
    for (QHash<DofID,Angle>::const_iterator it = angle_deltas.constBegin();
         it != angle_deltas.constEnd();
         ++it)
    {
        if (it.key().isDihedral())
            dihedrals.append( DihedralID(it.key().atom0(), it.key().atom1(),
                                         it.key().atom2(), it.key().atom3()) );
    }
    
    return dihedrals;
}

const char* Flexibility::typeName()
{
    return QMetaType::typeName(qMetaTypeId<Flexibility>());
}
