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

#include "AmberParams.h"

#include "SireMol/molecule.h"
#include "SireMol/partialmolecule.h"
#include "SireMol/bondid.h"
#include "SireMol/angleid.h"
#include "SireMol/dihedralid.h"
#include "SireMol/improperid.h"
#include "SireMol/atomidx.h"

#include "SireCAS/expression.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMol;
using namespace SireCAS;
using namespace SireMM;
using namespace SireStream;

///////////
/////////// Implementation of AmberBond
///////////

QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, const AmberBond &bond)
{
    ds << bond._k << bond._r0;
    return ds;
}

QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, AmberBond &bond)
{
    ds >> bond._k >> bond._r0;
    return ds;
}

AmberBond::AmberBond(const Expression &f) : _k(0), _r0(0)
{}

AmberBond::~AmberBond()
{}

///////////
/////////// Implementation of AmberAngle
///////////

QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, const AmberAngle &angle)
{
    ds << angle._k << angle._theta0;
    return ds;
}

QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, AmberAngle &angle)
{
    ds >> angle._k >> angle._theta0;
    return ds;
}

AmberAngle::AmberAngle(const Expression &f) : _k(0), _theta0(0)
{}

AmberAngle::~AmberAngle()
{}

///////////
/////////// Implementation of AmberDihPart
///////////

QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, const AmberDihPart &dih)
{
    ds << dih._k << dih._periodicity << dih._phase;
    return ds;
}

QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, AmberDihPart &dih)
{
    ds >> dih._k >> dih._periodicity >> dih._phase;
    return ds;
}

AmberDihPart::~AmberDihPart()
{}

///////////
/////////// Implementation of AmberDihedral
///////////

QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, const AmberDihedral &dih)
{
    ds << dih._parts;
    return ds;
}

QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, AmberDihedral &dih)
{
    ds >> dih._parts;
    return ds;
}

AmberDihedral::AmberDihedral(const Expression &f)
{}

AmberDihedral::~AmberDihedral()
{}

///////////
/////////// Implementation of AmberNB14
///////////

QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, const AmberNB14 &nb)
{
    ds << nb._cscl << nb._ljscl;
    return ds;
}

QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, AmberNB14 &nb)
{
    ds >> nb._cscl >> nb._ljscl;
    return ds;
}

AmberNB14::~AmberNB14()
{}

///////////
/////////// Implementation of AmberParams
///////////

static const RegisterMetaType<AmberParams> r_amberparam;

/** Serialise to a binary datastream */
QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, const AmberParams &amberparam)
{
    writeHeader(ds, r_amberparam, 1);

    SharedDataStream sds(ds);

    sds << amberparam.molinfo << amberparam.bonds 
        << amberparam.angles << amberparam.dihedrals
        << amberparam.impropers << amberparam.nb14pairs
        << static_cast<const MoleculeProperty&>(amberparam);

    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, AmberParams &amberparam)
{
    VersionID v = readHeader(ds, r_amberparam);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> amberparam.molinfo >> amberparam.bonds
            >> amberparam.angles >> amberparam.dihedrals 
            >> amberparam.impropers >> amberparam.nb14pairs
            >> static_cast<MoleculeProperty&>(amberparam);
    }
    else
        throw version_error(v, "1", r_amberparam, CODELOC);
        
    return ds;
}

/** Null Constructor */
AmberParams::AmberParams() : ConcreteProperty<AmberParams,MoleculeProperty>()
{}

/** Constructor for the passed molecule*/
AmberParams::AmberParams(const MoleculeView &mol)
            : ConcreteProperty<AmberParams,MoleculeProperty>(),
              molinfo(mol.data().info())
{}

/** Constructor for the passed molecule*/
AmberParams::AmberParams(const MoleculeInfoData &info)
            : ConcreteProperty<AmberParams,MoleculeProperty>(),
              molinfo(info)
{}

/** Copy constructor */
AmberParams::AmberParams(const AmberParams &other)
            : ConcreteProperty<AmberParams,MoleculeProperty>(),
              bonds(other.bonds),angles(other.angles),dihedrals(other.dihedrals),
              impropers(other.impropers),nb14pairs(other.nb14pairs)
{}

/** Copy assignment operator */
AmberParams& AmberParams::operator=(const AmberParams &other)
{
    if (this != &other)
    {
        MoleculeProperty::operator=(other);
        molinfo = other.molinfo;
        bonds = other.bonds;
        angles = other.angles;
        dihedrals = other.dihedrals;
        impropers = other.impropers;
        nb14pairs = other.nb14pairs;
    }
  
    return *this;
}

/** Destructor */
AmberParams::~AmberParams()
{}

/** Comparison operator */
bool AmberParams::operator==(const AmberParams &other) const
{
  return (molinfo == other.molinfo and bonds == other.bonds and angles == other.angles 
          and dihedrals == other.dihedrals and impropers == other.impropers
          and nb14pairs == other.nb14pairs);
}

/** Comparison operator */
bool AmberParams::operator!=(const AmberParams &other) const
{
    return not AmberParams::operator==(other);
}

/** Return the layout of the molecule whose flexibility is contained
    in this object */
const MoleculeInfoData& AmberParams::info() const
{
    if (molinfo.constData() == 0)
        return MoleculeInfoData::null();
    else
        return *molinfo;
}

/** Return whether or not this flexibility is compatible with the molecule 
    whose info is in 'molinfo' */
bool AmberParams::isCompatibleWith(const SireMol::MoleculeInfoData &molinfo) const
{
    return info().UID() == molinfo.UID();
}

const char* AmberParams::typeName()
{
    return QMetaType::typeName(qMetaTypeId<AmberParams>());
}

void AmberParams::add(const BondID &bond, const double &k, const double &r0) 
{
    bonds.insert( bond, AmberBond(k,r0) );
}

void AmberParams::remove(const BondID &bond)
{
    bonds.remove( bond );
}

AmberBond AmberParams::getParams(const BondID &bond)
{
    return bonds.value(bond);
}

QList<BondID> AmberParams::getAllBonds()
{
  return bonds.keys();
}

void AmberParams::add(const AngleID &angle, const double &k, const double &theta0) 
{
    angles.insert( angle, AmberAngle(k,theta0) );
}

void AmberParams::remove(const AngleID &angle)
{
    angles.remove( angle );
}

AmberAngle AmberParams::getParams(const AngleID &angle)
{
    return angles.value(angle);
}

QList<AngleID> AmberParams::getAllAngles()
{
    return angles.keys();
}

void AmberParams::add(const DihedralID &dihedral, const double &v, const double &periodicity,
                      const double &phase)
{
    // If dihedral already exists, we will append parameters
    if (dihedrals.contains(dihedral))
    {
        dihedrals[dihedral] += AmberDihPart(v, periodicity, phase);
    }
    else if (dihedrals.contains(dihedral.mirror()))
    {
        dihedrals[dihedral.mirror()] += AmberDihPart(v, periodicity, phase);
    }
    else
    {
        dihedrals.insert(dihedral, AmberDihPart(v, periodicity, phase));
    }
}

void AmberParams::remove(const DihedralID &dihedral)
{
    dihedrals.remove( dihedral );
}

AmberDihedral AmberParams::getParams(const DihedralID &dihedral)
{
    return dihedrals.value(dihedral);
}

QList<DihedralID> AmberParams::getAllDihedrals()
{
    return dihedrals.keys();
}

void AmberParams::add(const ImproperID &improper, const double &v,
                      const double &periodicity, const double &phase)
{
    if (impropers.contains(improper))
    {
        impropers[improper] += AmberDihPart(v, periodicity, phase);
    }
    else
    {
        impropers.insert(improper, AmberDihPart(v, periodicity, phase));
    }
}

void AmberParams::remove(const ImproperID &improper)
{
    impropers.remove( improper );
}

AmberDihedral AmberParams::getParams(const ImproperID &improper)
{
    return impropers.value( improper );
}

QList<ImproperID> AmberParams::getAllImpropers()
{
    return impropers.keys();
}

void AmberParams::add14Pair(const BondID &pair, const double &cscl, const double &ljscl) 
{
    nb14pairs.insert( pair, AmberNB14(cscl,ljscl) );
}

void AmberParams::remove14Pair(const BondID &pair)
{
    nb14pairs.remove( pair );
}

AmberNB14 AmberParams::get14PairParams(const BondID &pair)
{
    return nb14pairs.value(pair);
}

QList<BondID> AmberParams::getAll14Pairs()
{
    return nb14pairs.keys();
}
