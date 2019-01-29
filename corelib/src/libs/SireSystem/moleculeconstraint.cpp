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

#include "moleculeconstraint.h"
#include "system.h"
#include "delta.h"

#include "SireMol/molecules.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <boost/shared_ptr.hpp>

#include <QDebug>

using namespace SireSystem;
using namespace SireMol;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<MoleculeConstraint> r_molconstraint( MAGIC_ONLY, 
                                                       MoleculeConstraint::typeName() );

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                          const MoleculeConstraint &molconstraint)
{
    writeHeader(ds, r_molconstraint, 2);
    
    ds << static_cast<const Constraint&>(molconstraint);
       
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, 
                                          MoleculeConstraint &molconstraint)
{
    VersionID v = readHeader(ds, r_molconstraint);
    
    if (v == 2)
    {
        ds >> static_cast<Constraint&>(molconstraint);
    }
    else if (v == 1)
    {
        QUuid sysuid;
    
        ds >> sysuid
           >> static_cast<Constraint&>(molconstraint);
    }
    else
        throw version_error( v, "1,2", r_molconstraint, CODELOC );

    return ds;
}

/** Constructor */
MoleculeConstraint::MoleculeConstraint() : Constraint()
{}

/** Copy constructor */
MoleculeConstraint::MoleculeConstraint(const MoleculeConstraint &other)
                   : Constraint(other)
{}

/** Destructor */
MoleculeConstraint::~MoleculeConstraint()
{}

const char* MoleculeConstraint::typeName()
{
    return "SireSystem::MoleculeConstraint";
}

/** Copy assignment operator */
MoleculeConstraint& MoleculeConstraint::operator=(const MoleculeConstraint &other)
{
    if (this != &other)
        Constraint::operator=(other);
    
    return *this;
}

/** Comparison operator */
bool MoleculeConstraint::operator==(const MoleculeConstraint &other) const
{
    return Constraint::operator==(other);
}

/** Comparison operator */
bool MoleculeConstraint::operator!=(const MoleculeConstraint &other) const
{
    return not this->operator==(other);
}
