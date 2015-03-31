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

#include "i.h"
#include "expression.h"
#include "complexvalues.h"

#include "SireMaths/errors.h"

#include "SireStream/datastream.h"

using namespace SireStream;
using namespace SireCAS;

static const RegisterMetaType<SireCAS::I> r_i;

/** Serialise to a binary datastream */
QDataStream SIRECAS_EXPORT &operator<<(QDataStream &ds, const I &i)
{
    writeHeader(ds, r_i, 1) << static_cast<const Constant&>(i);

    return ds;
}

/** Deserialise from a binary datastream */
QDataStream SIRECAS_EXPORT &operator>>(QDataStream &ds, I &i)
{
    VersionID v = readHeader(ds, r_i);

    if (v == 1)
    {
        ds >> static_cast<Constant&>(i);
    }
    else
        throw version_error(v, "1", r_i, CODELOC);

    return ds;
}

/** Constructor */
I::I() : Constant()
{}

/** Copy constructor */
I::I(const I &other) : Constant(other)
{}

/** Destructor */
I::~I()
{}

/** Comparison operator */
bool I::operator==(const ExBase &other) const
{
    const I *other_i = dynamic_cast<const I*>(&other);

    return other_i != 0 and typeid(other).name() == typeid(*this).name();
}

/** Return a hash of this expression */
uint I::hash() const
{
    return ( r_i.magicID() << 16 ) | ( r_i.magicID() & 0x0000FFFF );
}

/** Return a string representation */
QString I::toString() const
{
    return "i";
}

/** Cannot evaluate 'i' as a real number, so throw a domain error */
double I::evaluate(const Values&) const
{
    throw SireMaths::domain_error(QObject::tr(
        "Cannot evaluate 'i' as a real number"), CODELOC);
}

/** Evaluate this as a complex number - return 'i' */
Complex I::evaluate(const ComplexValues&) const
{
    return Complex(0,1);
}

/** Return the complex conjugate of 'i' (-i) */
Expression I::conjugate() const
{
    return -1 * I();
}

/** I is definitely complex :-) */
bool I::isComplex() const
{
    return true;
}

const char* I::typeName()
{
    return QMetaType::typeName( qMetaTypeId<I>() );
}

I* I::clone() const
{
    return new I(*this);
}
