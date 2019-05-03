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

#ifndef SIREMATHS_TRIGMATRIX_H
#define SIREMATHS_TRIGMATRIX_H

#include <QVector>

#include "sireglobal.h"

#include <utility>

SIRE_BEGIN_HEADER

namespace SireMaths
{
class TrigMatrix;
}

SIREMATHS_EXPORT QDataStream& operator<<(QDataStream&, const SireMaths::TrigMatrix&);
SIREMATHS_EXPORT QDataStream& operator>>(QDataStream&, SireMaths::TrigMatrix&);

namespace SireMaths
{

class NMatrix;
class NVector;
class Vector;

/** This is a dense, double, general N*N 2-dimensional 
    triangular matrix, with m[i,j] == m[j,i].
    
    Only the upper diagonal half of the matrix is stored,
    in a packed format, which is not suitable for
    BLAS or LAPACK functions (a conversion to a NMatrix
    has to be performed first)
    
    The data is implicitly shared (copy on write), so 
    copying a matrix is very fast.
                  0 1 2 3 4 
    Matrix =  0 / a b c d e \    m[i,j] == m[j,i]
              1 | b f g h i |
              2 | c g j k l |    dimension == n == 5
              3 | d h k m n |
              4 \ e i l n o /
              
    This is stored in memory as;
           0 1 2 3 4 5 6 7 8 9 0 1 2 3 4
    v == [ a b c d e f g h i j k l m n o ]
    
    a == 0  (15-15), f == 5 (15-10), j == 9 (15-6), m == 12 (15-3), o == 14 (15-1)
    
    so m[i,j] == v[ S(n) - S(n-i) + j - i ]   if i >= j or
                 m[j,i]                       if j < i
    
    where S(k) == Sum_1^k (k) == (k^2 + k)/2
    
     S(1) == 1, S(2) == 3, S(3) == 6, S(4) == 10, S(5) == 15
    
    m[0,1] == 'b'   v[15 - 15 + 1 - 0 ] == v[1] == 'b'
    m[3,2] == 'k'   m[3,2] == m[2,3] == v[15 - 6 + 3 - 2] == v[10] == 'k'
    m[4,4] == 'o'   v[15 - 1 + 4 - 4 ] == v[14] == 'o'

    Simplified expression is m[i,j] = v[1/2 (-i - i^2 + 2 j + 2 i n)]
    
    @author Christopher Woods
*/
class SIREMATHS_EXPORT TrigMatrix
{

friend SIREMATHS_EXPORT QDataStream& ::operator<<(QDataStream&, const TrigMatrix&);
friend SIREMATHS_EXPORT QDataStream& ::operator>>(QDataStream&, TrigMatrix&);

public:
    TrigMatrix();
    
    TrigMatrix(int dimension);
    TrigMatrix(int dimension, double initial_value);
    
    TrigMatrix(const NMatrix &matrix, bool take_upper=true);
    
    TrigMatrix(const TrigMatrix &other);
    
    ~TrigMatrix();
    
    static const char* typeName();
    
    const char* what() const;
    
    TrigMatrix& operator=(const TrigMatrix &other);
    
    bool operator==(const TrigMatrix &other) const;
    bool operator!=(const TrigMatrix &other) const;
    
    double& operator()(int i, int j);
    const double& operator()(int i, int j) const;
    
    TrigMatrix& operator+=(const TrigMatrix &other);
    TrigMatrix& operator-=(const TrigMatrix &other);
    TrigMatrix& operator*=(const TrigMatrix &other);
    TrigMatrix& operator/=(const TrigMatrix &other);
    
    TrigMatrix& operator*=(double scale);
    TrigMatrix& operator/=(double scale);

    TrigMatrix operator-() const;
    
    TrigMatrix operator+(const TrigMatrix &other) const;
    TrigMatrix operator-(const TrigMatrix &other) const;
    TrigMatrix operator*(const TrigMatrix &other) const;
    TrigMatrix operator/(const TrigMatrix &other) const;
    
    TrigMatrix operator*(double scale) const;
    TrigMatrix operator/(double scale) const;
    
    NVector operator*(const NVector &vector) const;
    NVector operator*(const Vector &vector) const;
    
    int nRows() const;
    int nColumns() const;
    
    int count() const;
    int size() const;
    
    NVector row(int i) const;
    NVector column(int j) const;
    
    void set(int i, int j, double value);
    
    void setRow(int i, double value);
    void setRow(int i, const NVector &row);
    
    void setColumn(int j, double value);
    void setColumn(int j, const NVector &column);
    
    void setAll(double value);
    
    double* data();
    const double* data() const;
    const double* constData() const;
    
    QVector<double> memory() const;
    
    int offset(int i, int j) const;
    int checkedOffset(int i, int j) const;
    
    void redimension(int dimension);
    
    QString toString() const;

    double determinant() const;
    double trace() const;
    
    NVector diagonal() const;
    
    std::pair<NVector,NMatrix> diagonalise() const;

    TrigMatrix inverse() const;
    
    TrigMatrix transpose() const;
    TrigMatrix fullTranspose() const;

    bool isTransposed() const;

    void assertValidIndex(int i, int j) const;
    
    void assertValidRow(int i) const;
    void assertValidColumn(int j) const;
    
    void assertNRows(int nrows) const;
    void assertNColumns(int ncolumns) const;

    void assertSquare() const;

private:
    /** The raw data for the matrix */
    QVector<double> array;
    
    /** The number of rows in the matrix (square matrix!) */
    qint32 nrows;
};

/** Return the offset into the 1D array for the value
    at index [i,j] - note that this performs *NO* checking,
    and invalid input will result in invalid output. If you
    want to check the indicies, use checkOffset(int i, int j) */
inline int TrigMatrix::offset(int i, int j) const
{
    if (i <= j)
        return (2*(j + i*nrows) - i - i*i) / 2;
    else
        return (2*(i + j*nrows) - j - j*j) / 2;
}

/** Return the number of unique elements in this matrix
    (the size of the underlying 1D array) */
inline int TrigMatrix::count() const
{
	return array.count();
}

/** Return the number of unique elements in this matrix
    (the size of the underlying 1D array) */
inline int TrigMatrix::size() const
{
	return TrigMatrix::count();
}

}

Q_DECLARE_METATYPE( SireMaths::TrigMatrix )

SIRE_EXPOSE_CLASS( SireMaths::TrigMatrix )

SIRE_END_HEADER

#endif
