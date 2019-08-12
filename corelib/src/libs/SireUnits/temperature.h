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

#ifndef SIREUNITS_TEMPERATURE
#define SIREUNITS_TEMPERATURE

#include "dimensions.h"

SIRE_BEGIN_HEADER

namespace SireUnits
{

class Celsius;
class Fahrenheit;

namespace Dimension
{

//skip this completely when parsing with gccxml as it is broken!
#ifdef SKIP_BROKEN_GCCXML_PARTS

class Temperature
{
public:
    Temperature();
    Temperature(double);
    ~Temperature();

    operator double() const;
};

#endif // end of 'ifdef SKIP_BROKEN_GCCXML_PARTS'

class TempBase
{
friend class SireUnits::Celsius;
friend class SireUnits::Fahrenheit;

public:
    TempBase(double value = 0) : val(value)
    {}

    TempBase(const TempBase &other) : val(other.val)
    {}

    TempBase(const Temperature &temp) : val(temp)
    {}

    virtual ~TempBase()
    {}

    TempBase& operator=(const TempBase &other)
    {
        val = other.val;
        return *this;
    }

    TempBase& operator=(const Temperature &temp)
    {
        val = double(temp);
        return *this;
    }

    bool operator==(const TempBase &other) const
    {
        return val == other.val;
    }

    bool operator!=(const TempBase &other) const
    {
        return val != other.val;
    }

    bool operator==(const Temperature &temp) const
    {
        return val == double(temp);
    }

    bool operator!=(const Temperature &temp) const
    {
        return val != double(temp);
    }

    double value() const
    {
        return val;
    }

    QString toString() const
    {
        return QString("%1 %2").arg(this->convertFromInternal()).arg(this->unitString());
    }

    /** Convert this into a temperature object */
    operator Temperature() const
    {
        return Temperature(val);
    }

    operator double() const
    {
        return val;
    }

    double in(const TempBase &other) const
    {
        return other.convertFromInternal(val) / other.convertFromInternal();
    }

    double in(const Temperature &temp) const
    {
        return val * temp;
    }

    double to(const TempBase &other) const
    {
        return this->in(other);
    }

    virtual double convertToInternal(double value) const=0;
    virtual double convertFromInternal(double value) const=0;

    double convertFromInternal() const
    {
        return this->convertFromInternal(val);
    }

protected:
    virtual QString unitString() const
    {
        return "K";
    }

    /** This holds the temperature in internal units (K) */
    double val;
};

/** Construct a Unit from a TempBase */
SIRE_ALWAYS_INLINE Unit::Unit(const TempBase &temperature)
            : sclfac(temperature)
{}

} //end of namespace Dimension

class Celsius : public Dimension::TempBase
{

public:
    Celsius() : Dimension::TempBase(1)
    {}

    explicit Celsius(double value) : Dimension::TempBase()
    {
        val = convertToInternal(value);
    }

    Celsius(const Dimension::Temperature &temp) : Dimension::TempBase(temp)
    {}

    Celsius(const Dimension::TempBase &other) : Dimension::TempBase(other)
    {}

    Celsius(const Celsius &other)
          : Dimension::TempBase(other)
    {}

    ~Celsius()
    {}

    double convertToInternal(double value) const
    {
        return value + 273.15;
    }

    double convertFromInternal(double value) const
    {
        return value - 273.15;
    }

    double convertFromInternal() const
    {
        return Dimension::TempBase::convertFromInternal();
    }

    Celsius& operator=(const Celsius &other)
    {
        Dimension::TempBase::operator=(other);
        return *this;
    }

    Celsius& operator=(const Dimension::Temperature &temp)
    {
        Dimension::TempBase::operator=(temp);
        return *this;
    }

    Celsius operator-() const
    {
        return Celsius(-convertFromInternal());
    }

    Celsius operator+(const Celsius &other) const
    {
        return Celsius(convertFromInternal() + other.convertFromInternal());
    }

    Celsius operator-(const Celsius &other) const
    {
        return Celsius(convertFromInternal() - other.convertFromInternal());
    }

    Celsius& operator+=(const Celsius &other)
    {
        convertToInternal( convertFromInternal() + other.convertFromInternal() );
        return *this;
    }

    Celsius& operator-=(const Celsius &other)
    {
        convertToInternal( convertFromInternal() - other.convertFromInternal() );
        return *this;
    }

    Celsius operator+(const Dimension::Temperature &other) const
    {
        return *this + Celsius(other);
    }

    Celsius operator-(const Dimension::Temperature &other) const
    {
        return *this - Celsius(other);
    }

    Celsius& operator+=(const Dimension::Temperature &other)
    {
        return this->operator+=(Celsius(other));
    }

    Celsius& operator-=(const Dimension::Temperature &other)
    {
        return this->operator-=(Celsius(other));
    }

    Celsius operator*(double value) const
    {
        return Celsius(value * convertFromInternal());
    }

    Celsius operator/(double value) const
    {
        return Celsius(value / convertFromInternal());
    }

    Celsius operator*(int value) const
    {
        return Celsius(value * convertFromInternal());
    }

    Celsius operator/(int value) const
    {
        return Celsius(value / convertFromInternal());
    }

protected:
    QString unitString() const
    {
        return "C";
    }
};

SIRE_ALWAYS_INLINE Celsius operator*(double value, const Celsius &temp)
{
    return temp * value;
}

#ifndef SKIP_BROKEN_GCCXML_PARTS
SIRE_ALWAYS_INLINE Dimension::PhysUnit<0,0,0,0,-1,0,0> operator/(double value, const Celsius &temp)
{
    return Dimension::PhysUnit<0,0,0,0,-1,0,0>(value / temp.convertFromInternal());
}

SIRE_ALWAYS_INLINE Dimension::PhysUnit<0,0,0,0,-1,0,0> operator/(int value, const Celsius &temp)
{
    return Dimension::PhysUnit<0,0,0,0,-1,0,0>(value / temp.convertFromInternal());
}
#endif

SIRE_ALWAYS_INLINE Celsius operator*(int value, const Celsius &temp)
{
    return temp * value;
}

class Fahrenheit : public Dimension::TempBase
{

public:
    Fahrenheit() : Dimension::TempBase(1)
    {}

    explicit Fahrenheit(double value) : Dimension::TempBase()
    {
        val = convertToInternal(value);
    }

    Fahrenheit(const Dimension::Temperature &temp) : Dimension::TempBase(temp)
    {}

    Fahrenheit(const Dimension::TempBase &other) : Dimension::TempBase(other)
    {}

    Fahrenheit(const Fahrenheit &other)
          : Dimension::TempBase(other)
    {}

    ~Fahrenheit()
    {}

    double convertToInternal(double value) const
    {
        return (value + 459.67) / 1.8;
    }

    double convertFromInternal(double value) const
    {
        return (value * 1.8) - 459.67;
    }

    double convertFromInternal() const
    {
        return Dimension::TempBase::convertFromInternal();
    }

    Fahrenheit& operator=(const Fahrenheit &other)
    {
        Dimension::TempBase::operator=(other);
        return *this;
    }

    Fahrenheit& operator=(const Dimension::Temperature &temp)
    {
        Dimension::TempBase::operator=(temp);
        return *this;
    }

    Fahrenheit operator-() const
    {
        return Fahrenheit(-convertFromInternal());
    }

    Fahrenheit operator+(const Fahrenheit &other) const
    {
        return Fahrenheit(convertFromInternal() + other.convertFromInternal());
    }

    Fahrenheit operator-(const Fahrenheit &other) const
    {
        return Fahrenheit(convertFromInternal() - other.convertFromInternal());
    }

    Fahrenheit& operator+=(const Fahrenheit &other)
    {
        convertToInternal( convertFromInternal() + other.convertFromInternal() );
        return *this;
    }

    Fahrenheit& operator-=(const Fahrenheit &other)
    {
        convertToInternal( convertFromInternal() - other.convertFromInternal() );
        return *this;
    }

    Fahrenheit operator+(const Dimension::Temperature &other) const
    {
        return *this + Fahrenheit(other);
    }

    Fahrenheit operator-(const Dimension::Temperature &other) const
    {
        return *this - Fahrenheit(other);
    }

    Fahrenheit& operator+=(const Dimension::Temperature &other)
    {
        return this->operator+=(Fahrenheit(other));
    }

    Fahrenheit& operator-=(const Dimension::Temperature &other)
    {
        return this->operator-=(Fahrenheit(other));
    }

    Fahrenheit operator*(double value) const
    {
        return Fahrenheit(value * convertFromInternal());
    }

    Fahrenheit operator/(double value) const
    {
        return Fahrenheit(value / convertFromInternal());
    }

    Fahrenheit operator*(int value) const
    {
        return Fahrenheit(value * convertFromInternal());
    }

    Fahrenheit operator/(int value) const
    {
        return Fahrenheit(value / convertFromInternal());
    }

protected:
    QString unitString() const
    {
        return "F";
    }
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

SIRE_ALWAYS_INLINE Fahrenheit operator*(double value, const Fahrenheit &temp)
{
    return temp * value;
}

SIRE_ALWAYS_INLINE Dimension::PhysUnit<0,0,0,0,-1,0,0> operator/(double value, const Fahrenheit &temp)
{
    return Dimension::PhysUnit<0,0,0,0,-1,0,0>(value / temp.convertFromInternal());
}

SIRE_ALWAYS_INLINE Dimension::PhysUnit<0,0,0,0,-1,0,0> operator/(int value, const Fahrenheit &temp)
{
    return Dimension::PhysUnit<0,0,0,0,-1,0,0>(value / temp.convertFromInternal());
}

SIRE_ALWAYS_INLINE Fahrenheit operator*(int value, const Fahrenheit &temp)
{
    return temp * value;
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

const Celsius celsius(1);
const Fahrenheit fahrenheit(1);

}

SIRE_EXPOSE_CLASS( SireUnits::Dimension::TempBase )
SIRE_EXPOSE_CLASS( SireUnits::Celsius )
SIRE_EXPOSE_CLASS( SireUnits::Fahrenheit )

SIRE_END_HEADER

#endif
