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

#include "amberparameters.h"

#include "SireMol/molecule.h"
#include "SireMol/partialmolecule.h"
#include "SireMol/bondid.h"
#include "SireMol/angleid.h"
#include "SireMol/dihedralid.h"
#include "SireMol/improperid.h"
#include "SireMol/atomidx.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMol;
using namespace SireStream;

///////////
/////////// Implementation of Amber Parameters
///////////

static const RegisterMetaType<AmberParameters> r_amberparam;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const AmberParameters &amberparam)
{
    //empty class so nothing to stream

    writeHeader(ds, r_amberparam, 2);

    SharedDataStream sds(ds);

    sds << amberparam.molinfo << amberparam.bonds 
	<< amberparam.angles << amberparam.dihedrals 
	<< amberparam.impropers << amberparam.nb14pairs << static_cast<const MoleculeProperty&>(amberparam); 

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, AmberParameters &amberparam)
{
    //empty class so nothing to stream

    VersionID v = readHeader(ds, r_amberparam);
    
    if (v != 1 and v != 2)
        throw version_error( v, "1", r_amberparam, CODELOC );
    else if (v == 2)
      {
	SharedDataStream sds(ds);

	sds >> amberparam.molinfo >> amberparam.bonds
	    >> amberparam.angles >> amberparam.dihedrals 
	    >> amberparam.impropers >> amberparam.nb14pairs >> static_cast<MoleculeProperty&>(amberparam); 
      }
    else
      {
	SharedDataStream sds(ds);

	sds >> amberparam.molinfo >> amberparam.bonds
	    >> amberparam.angles >> amberparam.dihedrals 
	    >> amberparam.impropers >> static_cast<MoleculeProperty&>(amberparam); 
      }
            
    return ds;
}

/** Null Constructor */
AmberParameters::AmberParameters() : ConcreteProperty<AmberParameters,MoleculeProperty>()
{}

/** Constructor for the passed molecule*/
AmberParameters::AmberParameters(const MoleculeData &molecule)
            : ConcreteProperty<AmberParameters,MoleculeProperty>(),
              molinfo(molecule.info())
{}


/** Copy constructor */
AmberParameters::AmberParameters(const AmberParameters &other)
                : ConcreteProperty<AmberParameters,MoleculeProperty>(),
		  molinfo(other.molinfo), bonds(other.bonds),angles(other.angles),
          dihedrals(other.dihedrals),impropers(other.impropers),nb14pairs(other.nb14pairs)
{}

/** Copy assignment operator */
AmberParameters& AmberParameters::operator=(const AmberParameters &other)
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
AmberParameters::~AmberParameters()
{}

/** Comparison operator */
bool AmberParameters::operator==(const AmberParameters &other) const
{
  return (molinfo == other.molinfo and bonds == other.bonds and angles == other.angles 
	  and dihedrals == other.dihedrals and impropers == other.impropers and nb14pairs == other.nb14pairs);
}

/** Comparison operator */
bool AmberParameters::operator!=(const AmberParameters &other) const
{
    return not AmberParameters::operator==(other);
}


/** Return the layout of the molecule whose flexibility is contained
    in this object */
const MoleculeInfoData& AmberParameters::info() const
{
    if (molinfo.constData() == 0)
        return MoleculeInfoData::null();
    else
        return *molinfo;
}

/** Return whether or not this flexibility is compatible with the molecule 
    whose info is in 'molinfo' */
bool AmberParameters::isCompatibleWith(const SireMol::MoleculeInfoData &molinfo) const
{
    return info().UID() == molinfo.UID();
}

const char* AmberParameters::typeName()
{
    return QMetaType::typeName(qMetaTypeId<AmberParameters>());
}

void AmberParameters::add(const BondID &bond, const double &k, const double &r0) 
{
  QList<double> params;
  params.append(k);
  params.append(r0);

  bonds.insert( bond, params);
}

void AmberParameters::remove(const BondID &bond)
{
  bonds.remove( bond );
}

QList<double> AmberParameters::getParams(const BondID &bond)
{
  return bonds.value(bond);
}

QList<BondID> AmberParameters::getAllBonds()
{
  return bonds.keys();
}

void AmberParameters::add(const AngleID &angle, const double &k, const double &theta0) 
{
  QList<double> params;
  params.append(k);
  params.append(theta0);

  angles.insert( angle, params);
}

void AmberParameters::remove(const AngleID &angle)
{
  angles.remove( angle );
}

QList<double> AmberParameters::getParams(const AngleID &angle)
{
  return angles.value(angle);
}

QList<AngleID> AmberParameters::getAllAngles()
{
  return angles.keys();
}

void AmberParameters::add(const DihedralID &dihedral, const double &v, const double &periodicity,
                          const double &phase)
{
    QList<double> params;

    // If dihedral already exists, we will append parameters
    QHash< DihedralID,QList<double> >::const_iterator it = dihedrals.constFind(dihedral);
    
    bool is_mirrored = false;
    
    if (it != dihedrals.constEnd())
    {
        params = *it;
    }
    else
    {
        it = dihedrals.constFind(dihedral.mirror());
        
        if (it != dihedrals.constEnd())
        {
            params = *it;
            is_mirrored = true;
        }
    }
    
    params.append(v);
    params.append(periodicity);
    params.append(phase);

    if (is_mirrored)
        dihedrals.insert( dihedral.mirror(), params );
    else
        dihedrals.insert( dihedral, params );
}

void AmberParameters::remove(const DihedralID &dihedral)
{
  dihedrals.remove( dihedral );
}

QList<double> AmberParameters::getParams(const DihedralID &dihedral)
{
  return dihedrals.value(dihedral);
}

QList<DihedralID> AmberParameters::getAllDihedrals()
{
  return dihedrals.keys();
}

void AmberParameters::add(const ImproperID &improper, const double &v,
                          const double &periodicity, const double &phase)
{
    QList<double> params;

    QHash< ImproperID, QList<double> >::const_iterator it = impropers.constFind(improper);
    
    if (it != impropers.constEnd())
    {
        params = *it;
    }

    params.append(v);
    params.append(periodicity);
    params.append(phase);

    impropers.insert( improper, params );
}

void AmberParameters::remove(const ImproperID &improper)
{
  impropers.remove( improper );
}

QList<double> AmberParameters::getParams(const ImproperID &improper)
{
  return impropers.value( improper );
}

QList<ImproperID> AmberParameters::getAllImpropers()
{
  return impropers.keys();
}

void AmberParameters::add14Pair(const BondID &pair, const double &cscl, const double &ljscl) 
{
  QList<double> params;
  params.append(cscl);
  params.append(ljscl);

  nb14pairs.insert( pair, params);
}

void AmberParameters::remove14Pair(const BondID &pair)
{
  nb14pairs.remove( pair );
}

QList<double> AmberParameters::get14PairParams(const BondID &pair)
{
  return nb14pairs.value(pair);
}

QList<BondID> AmberParameters::getAll14Pairs()
{
  return nb14pairs.keys();
}
