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

#include <QDataStream>

#include <boost/rational.hpp>

#include "rational.h"
#include "maths.h"

#include "SireStream/datastream.h"

using namespace SireStream;
using namespace SireMaths;

static const RegisterMetaType<Rational> r_rational(NO_ROOT);

/** Serialise a rational number to a binary datastream */
QDataStream SIREMATHS_EXPORT &operator<<(QDataStream &ds, const Rational &val)
{
    writeHeader(ds, r_rational, 1) << val.numerator() << val.denominator();
    return ds;
}

/** Deserialise a rational number from a binary datastream */
QDataStream SIREMATHS_EXPORT &operator>>(QDataStream &ds, Rational &val)
{
    VersionID v = readHeader(ds, r_rational);

    if (v == 1)
    {
        qint32 num,denom;
        ds >> num >> denom;

        val = Rational(num,denom);
    }
    else
        throw version_error(v, "1", r_rational, CODELOC);

    return ds;
}

Rational::~Rational()
{}

bool Rational::operator==(const Rational &other) const
{
    return boost::rational<qint32>(num,den) == boost::rational<qint32>(other.num,other.den);
}

bool Rational::operator!=(const Rational &other) const
{
    return boost::rational<qint32>(num,den) != boost::rational<qint32>(other.num,other.den);
}

bool Rational::operator>=(const Rational &other) const
{
    return boost::rational<qint32>(num,den) >= boost::rational<qint32>(other.num,other.den);
}

bool Rational::operator<=(const Rational &other) const
{
    return boost::rational<qint32>(num,den) <= boost::rational<qint32>(other.num,other.den);
}

bool Rational::operator>(const Rational &other) const
{
    return boost::rational<qint32>(num,den) > boost::rational<qint32>(other.num,other.den);
}

bool Rational::operator<(const Rational &other) const
{
    return boost::rational<qint32>(num,den) < boost::rational<qint32>(other.num,other.den);
}

bool Rational::operator==(qint32 number) const
{
    return boost::rational<qint32>(num,den) == boost::rational<qint32>(number);
}

bool Rational::operator!=(qint32 number) const
{
    return boost::rational<qint32>(num,den) != boost::rational<qint32>(number);
}

bool Rational::operator>=(qint32 number) const
{
    return boost::rational<qint32>(num,den) >= boost::rational<qint32>(number);
}

bool Rational::operator<=(qint32 number) const
{
    return boost::rational<qint32>(num,den) <= boost::rational<qint32>(number);
}

bool Rational::operator>(qint32 number) const
{
    return boost::rational<qint32>(num,den) > boost::rational<qint32>(number);
}

bool Rational::operator<(qint32 number) const
{
    return boost::rational<qint32>(num,den) < boost::rational<qint32>(number);
}

Rational& Rational::operator+=(const Rational &other)
{
    auto res = boost::rational<qint32>(num,den) + boost::rational<qint32>(other.num,other.den);
    num = res.numerator();
    den = res.denominator();
    return *this;
}

Rational& Rational::operator-=(const Rational &other)
{
    auto res = boost::rational<qint32>(num,den) - boost::rational<qint32>(other.num,other.den);
    num = res.numerator();
    den = res.denominator();
    return *this;
}

Rational& Rational::operator*=(const Rational &other)
{
    auto res = boost::rational<qint32>(num,den) * boost::rational<qint32>(other.num,other.den);
    num = res.numerator();
    den = res.denominator();
    return *this;
}

Rational& Rational::operator/=(const Rational &other)
{
    auto res = boost::rational<qint32>(num,den) / boost::rational<qint32>(other.num,other.den);
    num = res.numerator();
    den = res.denominator();
    return *this;
}

Rational& Rational::operator*=(qint32 number)
{
    auto res = boost::rational<qint32>(num,den) * boost::rational<qint32>(number);
    num = res.numerator();
    den = res.denominator();
    return *this;
}

Rational& Rational::operator/=(qint32 number)
{
    auto res = boost::rational<qint32>(num,den) / boost::rational<qint32>(number);
    num = res.numerator();
    den = res.denominator();
    return *this;
}

QString Rational::toString() const
{
    return SireMaths::toString(*this);
}

namespace SireMaths
{

/** Expose the boost gcd and lcm functions */
qint32 SIREMATHS_EXPORT gcd(qint32 n, qint32 m)
{
    return boost::gcd<qint32>(n,m);
}

/** Expose the boost gcd and lcm functions */
qint32 SIREMATHS_EXPORT lcm(qint32 n, qint32 m)
{
    return boost::lcm<qint32>(n,m);
}

QString SIREMATHS_EXPORT toString(const SireMaths::Rational &val)
{
    if (val.denominator() == 1)
        return QString::number(val.numerator());
    else
        return QString("%1/%2").arg(val.numerator()).arg(val.denominator());
}

/** Private function used by 'pow' to calculate 'x' raised to the positive
    fractional power 'power' */
double SIREMATHS_EXPORT pow_pvt(double x, const Rational &power)
{
    if ( x == 0 )
        return 0;
        
    else if ( x > 0.0 or SireMaths::isEven(power.numerator()) or 
              SireMaths::isOdd(power.denominator()) )
    {
        switch(power.denominator())
        {
            case 2:
                return std::sqrt( SireMaths::pow(x, power.numerator()) );
            default:
                return std::exp( std::log( SireMaths::pow(x, power.numerator()) ) 
                                              / power.denominator() );
        }
    }
    else 
        throw SireMaths::domain_error(
            QObject::tr("Cannot raise the negative number '%1' to a fractional "
                        "power (%2)").arg(x).arg(toString(power)), CODELOC);
}

/** Return x raised to the fractional power 'power' */
double SIREMATHS_EXPORT pow(double x, const Rational &power)
{
    if (power.denominator() == 1)
        return SireMaths::pow(x, power.numerator());
    else if (power > 0)
        return pow_pvt(x,power);
    else
        return double(1.0) / pow_pvt(x, -power);
}

/** Return whether this is a rational number (with maximum denominator = 'maxdenom') */
bool SIREMATHS_EXPORT isRational(double val, int maxdenom)
{
    for (int i=1; i<=maxdenom; ++i)
    {
        int ival = int( val * double(i) );
        
        double error = std::abs( val - (double(ival)/double(i)) );
        
        if (error < std::numeric_limits<double>::epsilon())
            return true;
    }
    
    return false;
}

/** Return 'val' converted into the best approximated rational number
    with maximum denominator 'maxdenom'. A perfect conversion will only 
    result if 'isRational(val,maxdenom)' returned true. */
Rational SIREMATHS_EXPORT toRational(double val, int maxdenom)
{
    Rational best_rational;
    double lowest_error = 0.0;

    //note, would be better if only tested primes...
    for (int i=1; i<=maxdenom; ++i)
    {
        int ival = int( val * double(i) );
        
        double error = std::abs( val - (double(ival)/double(i)) );
        
        if (error < std::numeric_limits<double>::epsilon())
            return Rational(ival, i);
        else if (i == 1 or error < lowest_error)
        {
            lowest_error = error;
            best_rational = Rational(ival, i);
        }
    }
    
    return best_rational;
}

/** Default value of maxdenom for toRational and isRational */
const int default_maxdenom = 500;

bool SIREMATHS_EXPORT isRational(double val)
{
    return isRational(val,default_maxdenom);
}

Rational SIREMATHS_EXPORT toRational(double val)
{
    return toRational(val,default_maxdenom);
}

/** Return 'val' converted to a double */
double SIREMATHS_EXPORT toDouble(const Rational &val)
{
    return double(val.numerator()) / double(val.denominator());
}

/** Return a hash of the rational number */
uint SIREMATHS_EXPORT qHash(const Rational &val)
{
    return (val.numerator()<<16) | (val.denominator() & 0x0000FFFF);
}

} // end of namespace SireMaths

