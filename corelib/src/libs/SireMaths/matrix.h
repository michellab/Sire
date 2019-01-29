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

#ifndef SIREMATHS_MATRIX_H
#define SIREMATHS_MATRIX_H

#include <QString>

#include <utility>

#include "vector.h"

#include <gsl/gsl_matrix.h>

#include <boost/tuple/tuple.hpp>

SIRE_BEGIN_HEADER

namespace SireMaths
{
class Matrix;
}

class QDataStream;
SIREMATHS_EXPORT QDataStream& operator<<(QDataStream&, const SireMaths::Matrix&);
SIREMATHS_EXPORT QDataStream& operator>>(QDataStream&, SireMaths::Matrix&);

namespace SireMaths
{

class Vector;
class NMatrix;

SIREMATHS_EXPORT const Matrix operator+(const Matrix &m1, const Matrix &m2);
SIREMATHS_EXPORT const Matrix operator-(const Matrix &m1, const Matrix &m2);
SIREMATHS_EXPORT const Matrix operator*(const Matrix &m1, const Matrix &m2);
SIREMATHS_EXPORT const Vector operator*(const Matrix &m, const Vector &p);
SIREMATHS_EXPORT const Matrix operator*(const Matrix &m, double c);
SIREMATHS_EXPORT const Matrix operator*(double c, const Matrix &m);

/**
This class represents a 3x3 square matrix, used to represent 3D transformations.

@author Christopher Woods
*/
class SIREMATHS_EXPORT Matrix
{

friend class Quaternion;

friend QDataStream& ::operator<<(QDataStream&, const Matrix&);
friend QDataStream& ::operator>>(QDataStream&, Matrix&);

public:
    Matrix();

    Matrix(double diagonal_value);

    Matrix(double xx, double xy, double xz,
           double yx, double yy, double yz,
           double zx, double zy, double zz);

    Matrix(const Vector& r1, const Vector& r2, const Vector& r3);
    Matrix(const tuple<Vector,Vector,Vector> &rows);

    Matrix(const NMatrix &m);
    Matrix(const gsl_matrix *m);

    Matrix(const Matrix& m);

    ~Matrix();

    static const char* typeName();
    
    const char* what() const
    {
        return Matrix::typeName();
    }

    const double& operator()(int i, int j) const;
    double& operator()(int i, int j);
    
    double* data();
    
    const double* data() const;
    const double* constData() const;

    int offset(int i, int j) const;
    int checkedOffset(int i, int j) const;

    double at(int i, int j) const;

    Matrix transpose() const;
    Matrix inverse() const;
    Vector trace() const;

    QString toString() const;

    bool isSymmetric() const;
    void enforceSymmetric();
    Matrix getPrincipalAxes() const;

    boost::tuple<Vector,Matrix> diagonalise() const;

    boost::tuple<Matrix,Matrix,Matrix> svd() const;
    boost::tuple<Matrix,Matrix,Matrix> singleValueDecomposition() const;

    double xx() const;
    double xy() const;
    double xz() const;
    
    double yx() const;
    double yy() const;
    double yz() const;
    
    double zx() const;
    double zy() const;
    double zz() const;

    Vector column0() const;
    Vector column1() const;
    Vector column2() const;

    Vector row0() const;
    Vector row1() const;
    Vector row2() const;

    double determinant() const;

    void setToIdentity();

    bool isIdentity() const;

    static Matrix identity();
    static Matrix zero();

    Matrix& operator=(const Matrix &other);

    bool operator==(const Matrix& m) const;
    bool operator!=(const Matrix& m) const;

    Matrix& operator+=(const Matrix &m);
    Matrix& operator-=(const Matrix &m);
    Matrix& operator*=(const Matrix &m);
    Matrix& operator*=(double c);
    Matrix& operator/=(double c);

    friend const Matrix operator+(const Matrix &m1, const Matrix &m2);
    friend const Matrix operator-(const Matrix &m1, const Matrix &m2);
    friend const Matrix operator*(const Matrix &m1, const Matrix &m2);
    friend const Vector operator*(const Matrix &m, const Vector &p);
    friend const Matrix operator*(const Matrix &m, double c);
    friend const Matrix operator*(double c, const Matrix &m);

    static Matrix covariance(const QVector<Vector> &p,
                             const QVector<Vector> &q,
                             int n=-1);

protected:
    /** The components of the matrix */
    double array[9];
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Return the offset into the array of the value at index [i,j] */
inline int Matrix::offset(int i, int j) const
{
    return 3*i + j;
}

/** Return the xx element (matrix[0,0]) */
inline double Matrix::xx() const
{
    return array[0];
}

/** Return the xy element (matrix[0,1]) */
inline double Matrix::xy() const
{
    return array[1];
}

/** Return the xz element (matrix[0,2]) */
inline double Matrix::xz() const
{
    return array[2];
}

/** Return the yx element (matrix[1,0]) */
inline double Matrix::yx() const
{
    return array[3];
}

/** Return the yy element (matrix[1,1]) */
inline double Matrix::yy() const
{
    return array[4];
}

/** Return the yz element (matrix[1,2]) */
inline double Matrix::yz() const
{
    return array[5];
}

/** Return the zx element (matrix[2,0]) */
inline double Matrix::zx() const
{
    return array[6];
}

/** Return the zy element (matrix[2,1]) */
inline double Matrix::zy() const
{
    return array[7];
}

/** Return the zz element (matrix[2,2]) */
inline double Matrix::zz() const
{
    return array[8];
}

#endif // #ifndef SIRE_SKIP_INLINE_FUNCTIONS

}

Q_DECLARE_METATYPE(SireMaths::Matrix)
Q_DECLARE_TYPEINFO(SireMaths::Matrix, Q_MOVABLE_TYPE);

SIRE_EXPOSE_CLASS( SireMaths::Matrix )

SIRE_END_HEADER

#endif
