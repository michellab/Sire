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

#include "abs.h"
#include "values.h"
#include "complexvalues.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <cmath>

using namespace SireCAS;
using namespace SireStream;

static const RegisterMetaType<Abs> r_abs;

/** Serialise to a binary datastream */
QDataStream SIRE_EXPORT &operator<<(QDataStream &ds, const Abs &abs)
{
    writeHeader(ds, r_abs, 1) << static_cast<const SingleFunc&>(abs);

    return ds;
}

/** Deserialise from a binary datastream */
QDataStream SIRE_EXPORT &operator>>(QDataStream &ds, Abs &abs)
{
    VersionID v = readHeader(ds, r_abs);

    if (v == 1)
    {
        ds >> static_cast<SingleFunc&>(abs);
    }
    else
        throw version_error(v, "1", r_abs, CODELOC);

    return ds;
}

/** Construct an empty Abs(0) */
Abs::Abs() : SingleFunc()
{}

/** Construct abs(expression) */
Abs::Abs(const Expression &expression) : SingleFunc(expression)
{}

/** Copy constructor */
Abs::Abs(const Abs &other) : SingleFunc(other)
{}

/** Destructor */
Abs::~Abs()
{}

/** Comparison operator */
bool Abs::operator==(const ExBase &other) const
{
    const Abs *other_abs = dynamic_cast<const Abs*>(&other);
    return other_abs != 0 and typeid(other).name() == typeid(*this).name()
             and this->argument() == other_abs->argument();
}

/** Return the magic for this function */
uint Abs::magic() const
{
    return r_abs.magicID();
}

/** Evaluate this function */
double Abs::evaluate(const Values &values) const
{
    return std::abs( x().evaluate(values) );
}

/** Complex evaluation */
Complex Abs::evaluate(const ComplexValues &values) const
{
    Complex arg = x().evaluate(values);
    
    return Complex( std::abs(arg.real()), std::abs(arg.imag()) );
}

/** d |x| / dx = x / |x|  (signum function) */
Expression Abs::diff() const
{
    return x() / *this;
}

/** Integral of |x| = 0.5( x |x| ) + C */
Expression Abs::integ() const
{
    return 0.5 * (x() * *this);
}

const char* Abs::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Abs>() );
}

Abs* Abs::clone() const
{
    return new Abs(*this);
}
