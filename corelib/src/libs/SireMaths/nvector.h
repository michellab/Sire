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

#ifndef SIREMATHS_NVECTOR_H
#define SIREMATHS_NVECTOR_H

#include <QVector>

#include "sireglobal.h"

SIRE_BEGIN_HEADER

namespace SireMaths
{
class NVector;
}

QDataStream& operator<<(QDataStream&, const SireMaths::NVector&);
QDataStream& operator>>(QDataStream&, SireMaths::NVector&);

namespace SireMaths
{

class NMatrix;
class Vector;

/** This is a dense, double column Vector

    @author Christopher Woods
*/
class SIREMATHS_EXPORT NVector
{

friend QDataStream& ::operator<<(QDataStream&, const NVector&);
friend QDataStream& ::operator>>(QDataStream&, NVector&);

public:
    NVector();
    
    NVector(int nrows);
    NVector(int nrows, double initial_value);
    
    NVector(const Vector &vector);
    NVector(const QVector<double> &vector);
    
    NVector(const NVector &vector);
    
    ~NVector();
    
    static const char* typeName();

    const char* what() const;
    
    NVector& operator=(const NVector &other);
    
    bool operator==(const NVector &other) const;
    bool operator!=(const NVector &other) const;
    
    double& operator[](int i);

    double& operator()(int i);
    double& operator()(int i, int j);
    
    const double& operator[](int i) const;

    const double& operator()(int i) const;
    const double& operator()(int i, int j) const;
    
    NVector& operator+=(const NVector &other);
    NVector& operator-=(const NVector &other);
    
    NVector& operator*=(double scale);
    NVector& operator/=(double scale);
    
    NVector operator-() const;
    
    NVector operator+(const NVector &other) const;
    NVector operator-(const NVector &other) const;

    NVector operator*(double scale) const;
    NVector operator/(double scale) const;
    
    double* data();
    
    const double* data() const;
    const double* constData() const;
    
    void set(int i, double value);
    void set(int i, int j, double value);
    
    void setAll(double value);
    
    int count() const;
    int size() const;
    
    double length() const;
    double length2() const;
    
    NVector normalise() const;
    
    double sum() const;
    
    int nRows() const;
    int nColumns() const;
    
    bool isZero() const;
    
    QString toString() const;
    
    double dot(const NVector &other) const;
    
    NVector cross(const NVector &other) const;

    NMatrix transpose() const;

    void assertValidIndex(int i) const;
    void assertValidIndex(int i, int j) const;
    
    void assertNRows(int nrows) const;
    void assertNColumns(int ncolumns) const;

private:
    /** The raw data for the vector */
    QVector<double> array;
};

inline NVector operator*(double scale, const NVector &vector)
{
    return vector * scale;
}

}

Q_DECLARE_METATYPE( SireMaths::NVector )

SIRE_EXPOSE_CLASS( SireMaths::NVector )

SIRE_END_HEADER

#endif
