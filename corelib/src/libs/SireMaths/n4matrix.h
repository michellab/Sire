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

#ifndef SIREMATHS_N4MATRIX_H
#define SIREMATHS_N4MATRIX_H

#include <QVector>

#include "sireglobal.h"

SIRE_BEGIN_HEADER

namespace SireMaths
{
class N4Matrix;
}

SIREMATHS_EXPORT QDataStream& operator<<(QDataStream&, const SireMaths::N4Matrix&);
SIREMATHS_EXPORT QDataStream& operator>>(QDataStream&, SireMaths::N4Matrix&);

namespace SireBase
{
template<class T>
class Array2D;
}

namespace SireMaths
{

class NMatrix;
class NVector;

/** This is a dense, double, general N*M*L*K 4-dimensional matrix.
    The data is stored as a column-major 2D matrix of column-major
    2D matricies (so each 2D sub-matrix is suitable for
    use with Fortran BLAS or LAPACK functions). This is 
    designed for high speed.
    
    The data is implicitly shared (copy on write), so 
    copying a matrix is very fast.
    
    @author Christopher Woods
*/
class SIREMATHS_EXPORT N4Matrix
{

friend SIREMATHS_EXPORT QDataStream& ::operator<<(QDataStream&, const N4Matrix&);
friend SIREMATHS_EXPORT QDataStream& ::operator>>(QDataStream&, N4Matrix&);

public:
    N4Matrix();
    
    N4Matrix(int nbigrows, int nbigcolumns, int nrows, int columns);
    N4Matrix(int nbigrows, int nbigcolumn,
             int nrows, int ncolumns, double initial_value);
    
    N4Matrix(const NMatrix &matrix);
    N4Matrix(const SireBase::Array2D<NMatrix> &matrix);
    N4Matrix(const QVector< QVector< QVector< QVector<double> > > > &matrix);

    N4Matrix(const N4Matrix &other);
    
    ~N4Matrix();
    
    static const char* typeName();
    
    const char* what() const;
    
    N4Matrix& operator=(const N4Matrix &other);
    
    bool operator==(const N4Matrix &other) const;
    bool operator!=(const N4Matrix &other) const;
    
    double& operator()(int i, int j, int k, int l);
    const double& operator()(int i, int j, int k, int l) const;
    
    NMatrix operator()(int i, int j) const;
    
    N4Matrix& operator+=(const N4Matrix &other);
    N4Matrix& operator-=(const N4Matrix &other);
    
    N4Matrix& operator*=(double scale);
    N4Matrix& operator/=(double scale);

    N4Matrix operator-() const;
    
    N4Matrix operator+(const N4Matrix &other) const;
    N4Matrix operator-(const N4Matrix &other) const;
    
    N4Matrix operator*(double scale) const;
    N4Matrix operator/(double scale) const;
    
    void add(int i, int j, const NMatrix &matrix);
    void subtract(int i, int j, const NMatrix &matrix);
    
    int nBigRows() const;
    int nBigColumns() const;
    
    int nRows() const;
    int nColumns() const;
    
    NMatrix view(int i, int j) const;
    
    void set(int i, int j, const NMatrix &matrix);
    
    void set(int i, int j, int k, int l, double value);
    
    void setAll(double value);
    
    double* data();
    const double* data() const;
    const double* constData() const;
    
    QVector<double> memory() const;
    
    int offset(int i, int j, int k, int l) const;
    int checkedOffset(int i, int j, int k, int l) const;
    
    void redimension(int nbigrows, int nbigcolumns,
                     int nrows, int ncolumns);
    
    QString toString() const;

    void assertValidIndex(int i, int j, int k, int l) const;
    
    void assertValidBigRow(int i) const;
    void assertValidBigColumn(int j) const;
    
    void assertValidRow(int k) const;
    void assertValidColumn(int l) const;
    
    void assertNBigRows(int nbigrows) const;
    void assertNBigColumns(int nbigcolumns) const;
    
    void assertNRows(int nrows) const;
    void assertNColumns(int ncolumns) const;

private:
    /** The raw data for the matrix */
    QVector<double> array;
    
    /** The number of big rows, big columns, rows and columns in the matrix */
    qint32 nbigrows, nbigcolumns, nrows, ncolumns;
};

/** Return the offset into the 1D array for the value
    at index [i,j,k,l] - note that this performs *NO* checking,
    and invalid input will result in invalid output. If you
    want to check the indicies, use checkOffset(int i, int j, int k, int l) */
SIRE_ALWAYS_INLINE int N4Matrix::offset(int i, int j, int k, int l) const
{
    return k + (l*nrows) + (i*nrows*ncolumns) + (j*nrows*ncolumns*nbigrows);
}

}

Q_DECLARE_METATYPE( SireMaths::N4Matrix )

SIRE_EXPOSE_CLASS( SireMaths::N4Matrix )

SIRE_EXPOSE_ALIAS( SireBase::Array2D<SireMaths::NMatrix>, SireBase::Array2D_NMatrix_ )

#ifdef SIRE_INSTANTIATE_TEMPLATES
#include "SireBase/array2d.hpp"
#include "nmatrix.h"
template class SireBase::Array2D<SireMaths::NMatrix>;
#endif

SIRE_END_HEADER

#endif

