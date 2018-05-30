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

#ifndef SIRECAS_COMPLEX_H
#define SIRECAS_COMPLEX_H

/**
The documentation for a lot of these functions is copied directly from the 
documentation of the gsl functions of the equivalent name. 

The copyright information for GSL is;

########################################################################
GSL

More information about GSL can be found at the project homepage, 
http://www.gnu.org/software/gsl/.

Printed copies of this manual can be purchased from Network Theory Ltd at
http://www.network-theory.co.uk/gsl/manual/. The money raised from sales of the 
manual helps support the development of GSL.

Copyright  1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006 The GSL Team.

Permission is granted to copy, distribute and/or modify this document 
under the terms of the GNU Free Documentation License, Version 1.2 or any 
later version published by the Free Software Foundation; with the Invariant 
Sections being "GNU General Public License" and "Free Software Needs Free Documentation", 
the Front-Cover text being "A GNU Manual", and with the Back-Cover Text 
being (a) (see below). A copy of the license is included in the section 
entitled ?GNU Free Documentation License?.

(a) The Back-Cover Text is: "You have freedom to copy and modify this GNU Manual, 
like GNU software."

GSL Complex

The functions described in this chapter provide support for complex numbers. 
The algorithms take care to avoid unnecessary intermediate underflows and overflows, 
allowing the functions to be evaluated over as much of the complex plane as possible.

For multiple-valued functions the branch cuts have been chosen to follow the 
conventions of Abramowitz and Stegun in the Handbook of Mathematical Functions. 
The functions return principal values which are the same as those in GNU Calc, 
which in turn are the same as those in Common Lisp, The Language (Second Edition)1 
and the HP-28/48 series of calculators.

The complex types are defined in the header file gsl_complex.h, while the 
corresponding complex functions and arithmetic operations are defined 
in gsl_complex_math.h.

#######################################################################
*/

#include <gsl/gsl_complex.h>

#include <QString>

#include <complex>

#include "SireMaths/maths.h"

SIRE_BEGIN_HEADER

namespace SireMaths
{
class Complex;
}

class QDataStream;
QDataStream& operator<<(QDataStream&, const SireMaths::Complex&);
QDataStream& operator>>(QDataStream&, SireMaths::Complex&);

namespace SireMaths
{

class Rational;

/**
This class represents a complex number to the same precision as 'double'. 
This is merely a thin wrapper around the gsl_complex struct, and the 
gsl_complex functions.

(indeed, this is publically derived from gsl_complex, so you can use 
this class whereever you would normally use a gsl_complex)

@author Christopher Woods
*/
class SIREMATHS_EXPORT Complex : public gsl_complex
{
public:

    Complex(double r=0.0, double i=0.0);
    Complex(const gsl_complex &complex);

    template<typename T>
    Complex(const std::complex<T> &stdcomplex);

    Complex(const Complex &other);

    ~Complex();

    static const char* typeName();
    
    const char* what() const
    {
        return Complex::typeName();
    }

    double real() const;
    double imag() const;
    
    template<typename T>
    operator std::complex<T>() const;

    bool isReal() const;
    bool isPurelyComplex() const;
    
    bool isZero() const;

    QString toString() const;

    static Complex rect(double x, double y);
    static Complex polar(double r, double theta);

    void setRectangular(double x, double y);
    void setPolar(double r, double theta);

    void setReal(double x);
    void setImag(double y);

    bool operator==(const Complex &other) const;
    bool operator!=(const Complex &other) const;

    template<typename T>
    bool operator==(const std::complex<T> &stdcomplex) const;

    template<typename T>
    bool operator!=(const std::complex<T> &stdcomplex) const;

    Complex& operator=(const Complex &other);

    template<typename T>
    Complex& operator=(const std::complex<T> &stdcomplex);

    Complex& operator+=(const Complex &other);
    Complex& operator-=(const Complex &other);
    Complex& operator*=(const Complex &other);
    Complex& operator/=(const Complex &other);

    Complex operator-() const;


    bool operator==(double r) const;
    bool operator!=(double r) const;

    Complex& operator=(double r);

    Complex& operator+=(double r);
    Complex& operator-=(double r);
    Complex& operator*=(double r);
    Complex& operator/=(double r);

    double arg() const;
    double abs() const;
    double abs2() const;

    double logAbs() const;

    Complex conjugate() const;

    Complex inverse() const;

    Complex negative() const;
};

Complex operator+(const Complex &z0, const Complex &z1);
Complex operator-(const Complex &z0, const Complex &z1);
Complex operator*(const Complex &z0, const Complex &z1);
Complex operator/(const Complex &z0, const Complex &z1);
Complex operator+(const Complex &z, double x);
Complex operator-(const Complex &z, double x);
Complex operator*(const Complex &z, double x);
Complex operator/(const Complex &z, double x);
Complex operator+(double x, const Complex &z);
Complex operator-(double x, const Complex &z);
Complex operator*(double x, const Complex &z);
Complex operator/(double x, const Complex &z);
Complex sqrt(const Complex &z);
Complex sqrt_real(double x);
Complex pow(const Complex &z, const Complex &a);
Complex pow(const Complex &z, double x);
Complex pow(const Complex &z, int n);
Complex pow(const Complex &z, const Rational &r);
Complex pow(double x, const Complex &z);
Complex exp(const Complex &z);
Complex log(const Complex &z);
Complex log10(const Complex &z);
Complex log_b(const Complex &z, const Complex &b);
Complex sin(const Complex &z);
Complex cos(const Complex &z);
Complex tan(const Complex &z);
Complex sec(const Complex &z);
Complex csc(const Complex &z);
Complex cot(const Complex &z);
Complex arcsin(const Complex &z);
Complex arcsin_real(double z);
Complex arccos(const Complex &z);
Complex arccos_real(double z);
Complex arctan(const Complex &z);
Complex arcsec(const Complex &z);
Complex arcsec_real(double z);
Complex arccsc(const Complex &z);
Complex arccsc_real(double z);
Complex arccot(const Complex &z);
Complex sinh(const Complex &z);
Complex cosh(const Complex &z);
Complex tanh(const Complex &z);
Complex sech(const Complex &z);
Complex csch(const Complex &z);
Complex coth(const Complex &z);
Complex arcsinh(const Complex &z);
Complex arccosh(const Complex &z);
Complex arccosh_real(double z);
Complex arctanh(const Complex &z);
Complex arctanh_real(double z);
Complex arcsech(const Complex &z);
Complex arccsch(const Complex &z);
Complex arccoth(const Complex &z);

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Construct from a std::complex */
template<typename T>
Complex::Complex(const std::complex<T> &stdcomplex)
        : gsl_complex( gsl_complex_rect(stdcomplex.real(),stdcomplex.imag()) )
{}

/** Assignment from a std::complex */
template<typename T>
Complex& Complex::operator=(const std::complex<T> &stdcomplex)
{
    setReal( stdcomplex.real() );
    setImag( stdcomplex.imag() );
    return *this;
}

/** Implicit conversion to a std::complex */
template<typename T>
Complex::operator std::complex<T>() const
{
    return std::complex<T>(real(),imag());
}

/** Comparison to std::complex */
template<typename T>
bool Complex::operator==(const std::complex<T> &stdcomplex) const
{
    return SireMaths::areEqual(real(), stdcomplex.real()) and
           SireMaths::areEqual(imag(), stdcomplex.imag());
}

/** Comparison to std::complex */
template<typename T>
bool Complex::operator!=(const std::complex<T> &stdcomplex) const
{
    return not operator==(stdcomplex);
}

/** Comparison with std::complex */
template<typename T>
SIRE_INLINE_TEMPLATE
bool operator==(const std::complex<T> &stdcomplex, const Complex &complex)
{
    return complex == stdcomplex;
}

/** Comparison with std::complex */
template<typename T>
SIRE_INLINE_TEMPLATE
bool operator!=(const std::complex<T> &stdcomplex, const Complex &complex)
{
    return complex != stdcomplex;
}

/** This is the std sqrt function. This helps resolve namespace issues... */
inline double sqrt(double x)
{
    return std::sqrt(x);
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

Q_DECLARE_METATYPE(SireMaths::Complex)
Q_DECLARE_TYPEINFO(SireMaths::Complex, Q_MOVABLE_TYPE);

SIRE_EXPOSE_CLASS( SireMaths::Complex )

SIRE_END_HEADER

#endif
