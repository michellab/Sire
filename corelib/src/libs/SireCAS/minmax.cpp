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

#include "minmax.h"
#include "identities.h"
#include "expression.h"
#include "complexvalues.h"

#include "SireMaths/errors.h"

#include "SireStream/datastream.h"

using namespace SireStream;
using namespace SireCAS;

////////
//////// Implementation of Min
////////

static const RegisterMetaType<Min> r_min;

/** Serialise to a binary datastream */
QDataStream SIRECAS_EXPORT &operator<<(QDataStream &ds, const Min &min)
{
    writeHeader(ds, r_min, 1) << static_cast<const DoubleFunc&>(min);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIRECAS_EXPORT &operator>>(QDataStream &ds, Min &min)
{
    VersionID v = readHeader(ds, r_min);
    
    if (v == 1)
    {
        ds >> static_cast<DoubleFunc&>(min);
    }
    else
        throw version_error(v, "1", r_min, CODELOC);
        
    return ds;
}

/** Constructor */
Min::Min() : DoubleFunc()
{}

/** Construct min(x(), y()) */
Min::Min(const Expression &x, const Expression &y) : DoubleFunc(x,y)
{}

/** Copy constructor */
Min::Min(const Min &other) : DoubleFunc(other)
{}

/** Destructor */
Min::~Min()
{}

/** Comparison operator */
bool Min::operator==(const ExBase &other) const
{
    const Min *other_min = dynamic_cast<const Min*>(&other);

    return other_min != 0 and typeid(other).name() == typeid(*this).name()
                 and this->x() == other_min->x() 
                 and this->y() == other_min->y();
}

/** Evaluate this function */
double Min::evaluate(const Values &values) const
{
    return qMin( x().evaluate(values), y().evaluate(values) );
}

/** Complex evaluation */
Complex Min::evaluate(const ComplexValues &values) const
{
    Complex xval = x().evaluate(values);
    Complex yval = y().evaluate(values);

    if (xval == yval)
        return xval;
        
    else if (xval.isReal() and yval.isReal())
        return Complex( qMin(xval.real(),yval.real()), 0 );
        
    else    
        throw SireMaths::domain_error( QObject::tr(
            "It is not possible to order the two complex numbers %1 and %2. "
            "It is thus not possible to find the minimum of %3.")
                .arg(xval.toString(), yval.toString(), this->toString()), CODELOC );

    return Complex(0,0);
}

/** Return this as a function of x and y */
Expression Min::functionOf(const Expression &ex, const Expression &ey) const
{
     if (ex == x() and ey == y())
     {
        return Expression(*this);
     }
     else
        return Expression(Min(ex, ey));
}

/** Return the magic */
uint Min::magic() const
{
    return r_min.magicID();
}

const char* Min::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Min>() );
}

////////
//////// Implementation of Max
////////

static const RegisterMetaType<Max> r_max;

/** Serialise to a binary datastream */
QDataStream SIRECAS_EXPORT &operator<<(QDataStream &ds, const Max &max)
{
    writeHeader(ds, r_max, 1) << static_cast<const DoubleFunc&>(max);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIRECAS_EXPORT &operator>>(QDataStream &ds, Max &max)
{
    VersionID v = readHeader(ds, r_max);
    
    if (v == 1)
    {
        ds >> static_cast<DoubleFunc&>(max);
    }
    else
        throw version_error(v, "1", r_max, CODELOC);
        
    return ds;
}

/** Constructor */
Max::Max() : DoubleFunc()
{}

/** Construct min(x(), y()) */
Max::Max(const Expression &x, const Expression &y) : DoubleFunc(x,y)
{}

/** Copy constructor */
Max::Max(const Max &other) : DoubleFunc(other)
{}

/** Destructor */
Max::~Max()
{}

/** Comparison operator */
bool Max::operator==(const ExBase &other) const
{
    const Max *other_max = dynamic_cast<const Max*>(&other);

    return other_max != 0 and typeid(other).name() == typeid(*this).name()
                 and this->x() == other_max->x() 
                 and this->y() == other_max->y();
}

/** Evaluate this function */
double Max::evaluate(const Values &values) const
{
    return qMax( x().evaluate(values), y().evaluate(values) );
}

/** Complex evaluation */
Complex Max::evaluate(const ComplexValues &values) const
{
    Complex xval = x().evaluate(values);
    Complex yval = y().evaluate(values);

    if (xval == yval)
        return xval;
        
    else if (xval.isReal() and yval.isReal())
        return Complex( qMax(xval.real(),yval.real()), 0 );
        
    else    
        throw SireMaths::domain_error( QObject::tr(
            "It is not possible to order the two complex numbers %1 and %2. "
            "It is thus not possible to find the minimum of %3.")
                .arg(xval.toString(), yval.toString(), this->toString()), CODELOC );

    return Complex(0,0);
}

/** Return this as a function of x and y */
Expression Max::functionOf(const Expression &ex, const Expression &ey) const
{
     if (ex == x() and ey == y())
     {
        return Expression(*this);
     }
     else
        return Expression(Max(ex, ey));
}

/** Return the magic */
uint Max::magic() const
{
    return r_max.magicID();
}

const char* Max::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Max>() );
}

Max* Max::clone() const
{
    return new Max(*this);
}


Min* Min::clone() const
{
    return new Min(*this);
}

