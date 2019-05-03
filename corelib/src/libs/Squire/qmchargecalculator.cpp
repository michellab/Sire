/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2009  Christopher Woods
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

#include "qmchargecalculator.h"

#include "SireMol/molecule.h"
#include "SireMol/partialmolecule.h"
#include "SireMol/molecules.h"
#include "SireMol/molnum.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace Squire;
using namespace SireMol;
using namespace SireBase;
using namespace SireStream;

//////////////
////////////// Implementation of ChargeCalculator
//////////////

static const RegisterMetaType<QMChargeCalculator> r_chgcalc( MAGIC_ONLY,
                                                    QMChargeCalculator::typeName() );

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                          const QMChargeCalculator &chgcalc)
{
    writeHeader(ds, r_chgcalc, 2);
    
    ds << chgcalc.sclfac
       << static_cast<const Property&>(chgcalc);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds,
                                          QMChargeCalculator &chgcalc)
{
    VersionID v = readHeader(ds, r_chgcalc);
    
    if (v == 2)
    {
        ds >> chgcalc.sclfac >> static_cast<Property&>(chgcalc);
    }
    else if (v == 1)
    {
        ds >> static_cast<Property&>(chgcalc);
        
        chgcalc.sclfac = 1;
    }
    else
        throw version_error(v, "1,2", r_chgcalc, CODELOC);
        
    return ds;
}

/** Constructor */
QMChargeCalculator::QMChargeCalculator() : Property(), sclfac(1)
{}

/** Copy constructor */
QMChargeCalculator::QMChargeCalculator(const QMChargeCalculator &other)
                   : Property(other), sclfac(other.sclfac)
{}

/** Destructor */
QMChargeCalculator::~QMChargeCalculator()
{}

/** Copy assignment operator */
QMChargeCalculator& QMChargeCalculator::operator=(const QMChargeCalculator &other)
{
    sclfac = other.sclfac;
    Property::operator=(other);
    return *this;
}

/** Comparison operator */
bool QMChargeCalculator::operator==(const QMChargeCalculator &other) const
{
    return sclfac == other.sclfac and Property::operator==(other);
}

/** Comparison operator */
bool QMChargeCalculator::operator!=(const QMChargeCalculator &other) const
{
    return sclfac != other.sclfac or Property::operator!=(other);
}

const char* QMChargeCalculator::typeName()
{
    return "Squire::QMChargeCalculator";
}

/** Return the partial charges on the molecule 'molecule' */
AtomCharges QMChargeCalculator::calculate(const PartialMolecule &molecule,
                                        const PropertyMap &map) const
{
    return this->operator()(molecule, map);
}

/** Return the scale factor for the charges */
double QMChargeCalculator::scaleFactor() const
{
    return sclfac;
}

/** Set the scale factor for the charges */
void QMChargeCalculator::setScaleFactor(double scl)
{
    sclfac = scl;
}

const NullQMChargeCalculator& QMChargeCalculator::null()
{
    return *(create_shared_null<NullQMChargeCalculator>());
}

//////////////
////////////// Implementation of NullQMChargeCalculator
//////////////

static const RegisterMetaType<NullQMChargeCalculator> r_nullchg;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                          const NullQMChargeCalculator &nullchg)
{
    writeHeader(ds, r_nullchg, 1);
    
    ds << static_cast<const QMChargeCalculator&>(nullchg);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds,
                                          NullQMChargeCalculator &nullchg)
{
    VersionID v = readHeader(ds, r_nullchg);
    
    if (v == 1)
    {
        ds >> static_cast<QMChargeCalculator&>(nullchg);
    }
    else
        throw version_error(v, "1", r_nullchg, CODELOC);
        
    return ds;
}

/** Constructor */
NullQMChargeCalculator::NullQMChargeCalculator()
                     : ConcreteProperty<NullQMChargeCalculator,QMChargeCalculator>()
{}

/** Copy constructor */
NullQMChargeCalculator::NullQMChargeCalculator(const NullQMChargeCalculator &other)
                     : ConcreteProperty<NullQMChargeCalculator,QMChargeCalculator>(other)
{}

/** Destructor */
NullQMChargeCalculator::~NullQMChargeCalculator()
{}

const char* NullQMChargeCalculator::typeName()
{
    return QMetaType::typeName( qMetaTypeId<NullQMChargeCalculator>() );
}

/** Copy assignment operator */
NullQMChargeCalculator& NullQMChargeCalculator::operator=(const NullQMChargeCalculator &other)
{
    QMChargeCalculator::operator=(other);
    return *this;
}

/** Comparison operator */
bool NullQMChargeCalculator::operator==(const NullQMChargeCalculator &other) const
{
    return QMChargeCalculator::operator==(other);
}

/** Comparison operator */
bool NullQMChargeCalculator::operator!=(const NullQMChargeCalculator &other) const
{
    return QMChargeCalculator::operator!=(other);
}
    
/** Return the charges calculate for the passed molecule - this returns
    zero charges! */
AtomCharges NullQMChargeCalculator::operator()(const PartialMolecule &molecule,
                                             const PropertyMap&) const
{
    AtomCharges charges(molecule.data().info());
    
    return charges;
}

/** Return whether or not the change from 'oldmol' to 'newmol'
    may change the charges. Use this function to see if 
    a charge calculation is required */
bool NullQMChargeCalculator::mayChangeCharges(const PartialMolecule&,
                                              const PartialMolecule&,
                                              const PropertyMap&) const
{
    return false;
}
