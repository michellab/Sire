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

#include "integrationconstant.h"

#include "SireCAS/errors.h"

#include "SireStream/datastream.h"

using namespace SireStream;
using namespace SireCAS;

static const RegisterMetaType<IntegrationConstant> r_intconst;

/** Serialise to a binary datastream */
QDataStream SIRECAS_EXPORT &operator<<(QDataStream &ds, const IntegrationConstant &ic)
{
    writeHeader(ds, r_intconst, 1) << static_cast<const Symbol&>(ic);

    return ds;
}

/** Deserialise from a binary datastream */
QDataStream SIRECAS_EXPORT &operator>>(QDataStream &ds, IntegrationConstant &ic)
{
    VersionID v = readHeader(ds, r_intconst);

    if (v == 1)
    {
        ds >> static_cast<Symbol&>(ic);
    }
    else
        throw version_error(v, "1", r_intconst, CODELOC);

    return ds;
}

/** Constructor */
IntegrationConstant::IntegrationConstant()
                    : Symbol("C")
{}

/** Copy constructor */
IntegrationConstant::IntegrationConstant( const IntegrationConstant &other )
                    : Symbol(other)
{}

/** Destructor */
IntegrationConstant::~IntegrationConstant()
{}

/** Comparison operator */
bool IntegrationConstant::operator==(const ExBase &other) const
{
    const IntegrationConstant *other_ic = dynamic_cast<const IntegrationConstant*>(&other);

    return other_ic != 0 and typeid(other).name() == typeid(*this).name();
}

/** Return a hash for this object */
uint IntegrationConstant::hash() const
{
    return ( r_intconst.magicID() << 16 ) | ( r_intconst.magicID() << 16 );
}

/** Cannot integrate an expression containing an integration constant. This
    is to prevent integration constants from multiple integrations from
    appearing in the expression. */
Expression IntegrationConstant::integrate(const Symbol&) const
{
    throw SireCAS::unavailable_integral(QObject::tr(
        "Cannot integrate an IntegrationConstant. You must remove all "
        "integration constants from an expression before you can integrate "
        "it again (e.g. via ex.substitute(IntegrationConstant() == 0))."), CODELOC);

    return Expression();
}

const char* IntegrationConstant::typeName()
{
    return QMetaType::typeName( qMetaTypeId<IntegrationConstant>() );
}

IntegrationConstant* IntegrationConstant::clone() const
{
    return new IntegrationConstant(*this);
}

