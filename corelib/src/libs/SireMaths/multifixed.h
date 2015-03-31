/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2013  Christopher Woods
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

#ifndef SIREMATHS_MULTIFIXED_H
#define SIREMATHS_MULTIFIXED_H

#include "SireMaths/multifloat.h"

SIRE_BEGIN_HEADER

namespace SireMaths
{

/** This class provides a vector of fixed point numbers which
    is built to complement MultiFloat. The number of fixed point
    numbers in the vector is equal to the number of floats in MultiFloat.
    
    Note that this is a 64bit balanced fixed point representation,
    i.e. 32bits of precision are available both before and after
    the decimal point (so about 9 significant figures on both sides,
    note that the minus sign takes one of the bits before the decimal point).
 
    Note also that there are no functions that let you retrieve 
    individual fixed point numbers from this vector. Fixed point
    numbers will automatically be converted to/from doubles when
    returned to the user.
 
    @author Christopher Woods
*/
class SIREMATHS_EXPORT MultiFixed
{
public:
    MultiFixed();
    
    MultiFixed(double value);
    
    MultiFixed(const double *array, int size);
    MultiFixed(const QVector<double> &array);
    
    MultiFixed(const MultiFloat &value);
    MultiFixed(const MultiFixed &other);
    
    ~MultiFixed();
    
    static QVector<MultiFixed> fromArray(const QVector<double> &array);
    
    static QVector<double> toArray(const QVector<MultiFixed> &array);
    
    MultiFixed& operator=(const MultiFixed &other);
    
    bool operator==(const MultiFixed &other) const;
    bool operator!=(const MultiFixed &other) const;
    
    bool operator<(const MultiFixed &other) const;
    bool operator>(const MultiFixed &other) const;
    
    bool operator<=(const MultiFixed &other) const;
    bool operator>=(const MultiFixed &other) const;
    
    MultiFixed compareEqual(const MultiFixed &other) const;
    MultiFixed compareNotEqual(const MultiFixed &other) const;

    MultiFixed compareLess(const MultiFixed &other) const;
    MultiFixed compareGreater(const MultiFixed &other) const;
    
    MultiFixed compareLessEqual(const MultiFixed &other) const;
    MultiFixed compareGreaterEqual(const MultiFixed &other) const;
    
    const char* what() const;
    static const char* typeName();
    
    QString toString() const;
    QString toBinaryString() const;
    
    static int size();
    static int count();
    
    double operator[](int i) const;
    
    void set(int i, double value);
    
    double get(int i) const;
    
    MultiFixed operator-() const;

    MultiFixed operator+(const MultiFixed &other) const;
    MultiFixed operator-(const MultiFixed &other) const;
    MultiFixed operator*(const MultiFixed &other) const;
    MultiFixed operator/(const MultiFixed &other) const;
    
    MultiFixed& operator+=(const MultiFixed &other);
    MultiFixed& operator-=(const MultiFixed &other);
    MultiFixed& operator*=(const MultiFixed &other);
    MultiFixed& operator/=(const MultiFixed &other);
    
    MultiFixed operator!() const;
    MultiFixed operator&(const MultiFixed &other) const;
    MultiFixed operator|(const MultiFixed &other) const;
    MultiFixed operator^(const MultiFixed &other) const;

    MultiFixed& operator&=(const MultiFixed &other);
    MultiFixed& operator|=(const MultiFixed &other);
    MultiFixed& operator^=(const MultiFixed &other);

    MultiFixed logicalNot() const;
    
    MultiFixed logicalAnd(const MultiFixed &other) const;
    MultiFixed logicalAndNot(const MultiFixed &other) const;
    
    MultiFixed logicalOr(const MultiFixed &other) const;
    MultiFixed logicalXor(const MultiFixed &other) const;
    
    MultiFixed& multiplyAdd(const MultiFixed &val0, const MultiFixed &val1);
    
    MultiFixed max(const MultiFixed &other) const;
    MultiFixed min(const MultiFixed &other) const;
    
    MultiFixed reciprocal() const;

    MultiFixed sqrt() const;
    MultiFixed rsqrt() const;
    
    MultiFixed rotate() const;
    
    double sum() const;

private:
    union
    {
        qint64 a[MULTIFLOAT_SIZE];
    } v;
};

}

SIRE_EXPOSE_CLASS( SireMaths::MultiFixed )

SIRE_END_HEADER

#endif
