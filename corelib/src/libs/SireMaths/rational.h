/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2006  Christopher Woods
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

#ifndef SIREMATHS_RATIONAL_H
#define SIREMATHS_RATIONAL_H

#include <QString>
#include <QObject>

#include "SireMaths/maths.h"
#include "SireMaths/errors.h"

#include "sireglobal.h"

SIRE_BEGIN_HEADER

namespace SireMaths
{

/** Thin wrapper around boost::rational - cannot use directly as it kills
    compilers! */
class SIREMATHS_EXPORT Rational
{
public:
    Rational() : num(0), den(1)
    {}
    
    
    Rational(qint32 number) : num(number) , den(1)
    {}
    
    Rational(qint32 numerator, qint32 denominator)
        : num(numerator), den(denominator)
    {}
    
    Rational(const Rational &other) : num(other.num), den(other.den)
    {}
    
    ~Rational();
    
    Rational& operator=(const Rational &other)
    {
        num = other.num;
        den = other.den;
        return *this;
    }
    
    Rational& operator=(qint32 number)
    {
        num = number;
        den = 1;
        return *this;
    }
    
    bool operator==(const Rational &other) const;
    bool operator!=(const Rational &other) const;
    bool operator>=(const Rational &other) const;
    bool operator<=(const Rational &other) const;
    bool operator>(const Rational &other) const;
    bool operator<(const Rational &other) const;
    
    bool operator==(qint32 number) const;
    bool operator!=(qint32 number) const;
    bool operator>=(qint32 number) const;
    bool operator<=(qint32 number) const;
    bool operator>(qint32 number) const;
    bool operator<(qint32 number) const;
    
    Rational operator-() const
    {
        return Rational(-num,den);
    }
    
    Rational& operator+=(const Rational &other);
    Rational& operator-=(const Rational &other);
    Rational& operator*=(const Rational &other);
    Rational& operator/=(const Rational &other);
    
    Rational& operator+=(qint32 number)
    {
        num += number*den;
        return *this;
    }
    
    Rational& operator-=(qint32 number)
    {
        num -= number*den;
        return *this;
    }
    
    Rational& operator*=(qint32 number);
    Rational& operator/=(qint32 number);
    
    QString toString() const;
    
    const Rational& operator++()
    {
        num += den;
        return *this;
    }
    
    const Rational& operator--()
    {
        num -= den;
        return *this;
    }
    
    bool operator!() const
    {
        return !num;
    }

    qint32 numerator() const
    {
        return num;
    }
    
    qint32 denominator() const
    {
        return den;
    }

private:
    qint32 num, den;
};

/** Expose the boost gcd and lcm functions */
qint32 gcd(qint32 n, qint32 m);

/** Expose the boost gcd and lcm functions */
qint32 lcm(qint32 n, qint32 m);

}

class QDataStream;
QDataStream& operator<<(QDataStream&, const SireMaths::Rational&);
QDataStream& operator>>(QDataStream&, SireMaths::Rational&);

namespace SireMaths
{

/** Return a QString representation of the rational number */
QString toString(const SireMaths::Rational &val);

/** Private function used by 'pow' to calculate 'x' raised to the positive
    fractional power 'power' */
double pow_pvt(double x, const Rational &power);

/** Return x raised to the fractional power 'power' */
double pow(double x, const Rational &power);

/** Return whether this is a rational number (with maximum denominator = 'maxdenom') */
bool isRational(double val, int maxdenom);

/** Return 'val' converted into the best approximated rational number
    with maximum denominator 'maxdenom'. A perfect conversion will only 
    result if 'isRational(val,maxdenom)' returned true. */
Rational toRational(double val, int maxdenom);

bool isRational(double val);

Rational toRational(double val);

/** Return 'val' converted to a double */
double toDouble(const Rational &val);

/** Return a hash of the rational number */
uint qHash(const Rational &val);

}

Q_DECLARE_METATYPE(SireMaths::Rational)

SIRE_EXPOSE_CLASS(SireMaths::Rational)

SIRE_END_HEADER

#endif
