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

#include "atomfunctions.h"

#include "SireMol/moleculedata.h"
#include "SireMol/moleculeinfodata.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMM;
using namespace SireStream;

/////////
///////// Implementation of AtomFunction
/////////

QDataStream &operator<<(QDataStream &ds,
                                      const AtomFunction &atomfunc)
{
    SharedDataStream sds(ds);
    sds << atomfunc.func;
    
    return ds;
}

QDataStream &operator>>(QDataStream &ds,
                                      AtomFunction &atomfunc)
{
    SharedDataStream sds(ds);
    sds >> atomfunc.func;
    
    return ds;
}

/** Constructor */
AtomFunction::AtomFunction()
{}

/** Construct to hold the passed function */
AtomFunction::AtomFunction(const SireCAS::Expression &function)
             : func(function)
{}

/** Copy constructor */
AtomFunction::AtomFunction(const AtomFunction &other)
             : func(other.func)
{}

/** Destructor */
AtomFunction::~AtomFunction()
{}

/** Copy assignment operator */
AtomFunction& AtomFunction::operator=(const AtomFunction &other)
{
    func = other.func;
    return *this;
}

/** Comparison operator */
bool AtomFunction::operator==(const AtomFunction &other) const
{
    return func == other.func;
}

/** Comparison operator */
bool AtomFunction::operator!=(const AtomFunction &other) const
{
    return func != other.func;
}

/////////
///////// Implementation of AtomFunctions
/////////

static const RegisterMetaType<AtomFunctions> r_atomfuncs(MAGIC_ONLY,
                                                         "SireMM::AtomFunctions");

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                      const AtomFunctions &atomfuncs)
{
    writeHeader(ds, r_atomfuncs, 1);
    
    SharedDataStream sds(ds);
    sds << atomfuncs.molinfo << atomfuncs.symbls;
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds,
                                      AtomFunctions &atomfuncs)
{
    VersionID v = readHeader(ds, r_atomfuncs);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> atomfuncs.molinfo >> atomfuncs.symbls;
    }
    else
        throw version_error( v, "1", r_atomfuncs, CODELOC );
        
    return ds;
}

/** Constructor */
AtomFunctions::AtomFunctions() : MoleculeProperty()
{}

/** Construct to hold the atom functions for the molecule whose
    data is in 'moldata' */
AtomFunctions::AtomFunctions(const MoleculeData &moldata)
              : MoleculeProperty(),
                molinfo( moldata.info() )
{}

/** Construct to hold the atom functions for the molecule whose
    layout information is in 'molinfo' */
AtomFunctions::AtomFunctions(const MoleculeInfoData &info)
              : MoleculeProperty(),
                molinfo(info)
{}

/** Copy constructor */
AtomFunctions::AtomFunctions(const AtomFunctions &other)
              : MoleculeProperty(other),
                molinfo(other.molinfo), symbls(other.symbls)
{}

/** Destructor */
AtomFunctions::~AtomFunctions()
{}

/** Copy assignment operator */
AtomFunctions& AtomFunctions::operator=(const AtomFunctions &other)
{
    MoleculeProperty::operator=(other);
    molinfo = other.molinfo;
    symbls = other.symbls;
    
    return *this;
}

/** Comparison operator */
bool AtomFunctions::operator==(const AtomFunctions &other) const
{
    return info() == other.info() and MoleculeProperty::operator==(other);
}

/** Comparison operator */
bool AtomFunctions::operator!=(const AtomFunctions &other) const
{
    return info() != other.info() or MoleculeProperty::operator!=(other);
}

/** Return whether or not this property is compatible with the molecule
    whose layout information is in 'molinfo' */
bool AtomFunctions::isCompatibleWith(const MoleculeInfoData &molinfo) const
{
    return info() == molinfo;
}

/** Add the passed symbols to the index */
void AtomFunctions::addSymbols(const QSet<Symbol> &symbols)
{
    symbls += symbols;
}

/** Remove the passed symbols from the index */
void AtomFunctions::removeSymbols(const QSet<Symbol> &symbols)
{
    symbls -= symbols;
}

/** Clear all of the symbols */
void AtomFunctions::removeSymbols()
{
    symbls.clear();
}
