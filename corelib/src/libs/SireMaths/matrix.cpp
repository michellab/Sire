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

#include <QString>
#include <boost/scoped_array.hpp>
#include "matrix.h"
#include "vector.h"
#include "maths.h"
#include "nmatrix.h"
#include "nvector.h"

#include "SireMaths/errors.h"
#include "SireError/errors.h"

#include "SireStream/datastream.h"

#include "third_party/eig3/eig3.h" // CONDITIONAL_INCLUDE

#include <gsl/gsl_eigen.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include <QDebug>

using namespace SireMaths;
using namespace SireStream;

using boost::tuple;

static const RegisterMetaType<Matrix> r_matrix(NO_ROOT);

/** Serialise to a binary data stream */
QDataStream SIREMATHS_EXPORT &operator<<(QDataStream &ds, const Matrix &matrix)
{
    writeHeader(ds, r_matrix, 1);
    
    for (int i=0; i<9; ++i)
    {
        ds << matrix.array[i];
    }

    return ds;
}

/** Deserialise from a binary data stream */
QDataStream SIREMATHS_EXPORT &operator>>(QDataStream &ds, Matrix &matrix)
{
    VersionID v = readHeader(ds, r_matrix);

    if (v == 1)
    {
        for (int i=0; i<9; ++i)
        {
            ds >> matrix.array[i];
        }
    }
    else
        throw version_error(v, "1", r_matrix, CODELOC);

    return ds;
}

/** Construct a default Matrix (identity matrix) */
Matrix::Matrix()
{
    for (int i=0; i<8; ++i)
    {
        array[i] = 0;
    }
    
    array[0] = 1;
    array[4] = 1;
    array[8] = 1;
}

/** Construct a matrix whose diagonal elements equal 'diagonal_value'
    and whose off-diagonal elements equal zero */
Matrix::Matrix(double diagonal_value)
{
    for (int i=0; i<8; ++i)
    {
        array[i] = 0;
    }
    
    array[0] = diagonal_value;
    array[4] = diagonal_value;
    array[8] = diagonal_value;
}

/** Construct a Matrix. Elements listed as row 1, then
    row 2, then row 3. */
Matrix::Matrix(double xx, double xy, double xz,
               double yx, double yy, double yz,
               double zx, double zy, double zz)
{
    array[0] = xx;
    array[1] = xy;
    array[2] = xz;
    
    array[3] = yx;
    array[4] = yy;
    array[5] = yz;
    
    array[6] = zx;
    array[7] = zy;
    array[8] = zz;
}

/** Construct from an NMatrix */
Matrix::Matrix(const NMatrix &m)
{
    if (m.nRows() != 3 or m.nColumns() != 3)
        throw SireError::incompatible_error( QObject::tr(
                "You cannot construct a 3x3 matrix from an NMatrix of dimension "
                "%1x%2.")
                    .arg(m.nRows()).arg(m.nColumns()), CODELOC );

    const double *d = m.constData();

    if (m.isTransposed())
    {
        array[0] = d[0]; array[3] = d[1]; array[6] = d[2];
        array[1] = d[3]; array[4] = d[4]; array[7] = d[5];
        array[2] = d[6]; array[5] = d[7]; array[8] = d[8];
    }
    else
    {
        for (int i=0; i<9; ++i)
            array[i] = d[i];
    }
}

/** Construct from a GSL matrix. This must obviously be a 3x3 matrix! */
Matrix::Matrix(const gsl_matrix *m)
{
    if (m->size1 > 3 or m->size2 > 3)
        throw SireError::incompatible_error( QObject::tr(
                    "SireMaths::Matrix is a 3x3 matrix class and cannot be initialised "
                    "from a gsl_matrix of size %1x%2.")
                        .arg(m->size1).arg(m->size2), CODELOC );

    for (int i=0; i<9; ++i)
        array[i] = 0;
    
    for (int i=0; i<m->size1; ++i)
    {
        for (int j=0; j<m->size2; ++j)
        {
            this->operator()(i,j) = gsl_matrix_get(m,i,j);
        }
    }
}

/** Copy constructor */
Matrix::Matrix(const Matrix &m)
{
    memcpy(this->data(), m.constData(), 9*sizeof(double));
}

/** Construct a matrix from three vectors - each vector is a row */
Matrix::Matrix(const Vector &r1, const Vector &r2, const Vector &r3)
{
    memcpy(this->data(), r1.constData(), 3*sizeof(double));
    memcpy(this->data()+3, r2.constData(), 3*sizeof(double));
    memcpy(this->data()+6, r3.constData(), 3*sizeof(double));
}

/** Construct a matrix from a tuple of three vectors - each
    vector is a row */
Matrix::Matrix(const tuple<Vector,Vector,Vector> &rows)
{
    const Vector &r1 = rows.get<0>();
    const Vector &r2 = rows.get<1>();
    const Vector &r3 = rows.get<2>();
    
    memcpy(this->data(), r1.constData(), 3*sizeof(double));
    memcpy(this->data()+3, r2.constData(), 3*sizeof(double));
    memcpy(this->data()+6, r3.constData(), 3*sizeof(double));
}

/** Destructor */
Matrix::~Matrix()
{}

/** Return the offset into the array of the value at index [i,j]

    \throw SireError::invalid_index
*/
int Matrix::checkedOffset(int i, int j) const
{
    if (i < 0 or i > 2 or j < 0 or j > 2)
        throw SireError::invalid_index( QObject::tr(    
                "Invalid index for 3x3 matrix - [%1,%2]")
                    .arg(i).arg(j), CODELOC );
                    
    return offset(i,j);
}

/** Return the value at index [i,j]

    \throw SireError::invalid_index
*/
const double& Matrix::operator()(int i, int j) const
{
    return array[ checkedOffset(i,j) ];
}

/** Return the value at index [i,j]

    \throw SireError::invalid_index
*/
double& Matrix::operator()(int i, int j)
{
    return array[ checkedOffset(i,j) ];
}

/** Return a raw pointer to the data of this matrix */
double* Matrix::data()
{
    return array;
}

/** Return a raw pointer to the data of this matrix */
const double* Matrix::data() const
{
    return array;
}

/** Return a raw pointer to the data of this matrix */
const double* Matrix::constData() const
{
    return array;
}

/** Return the element at index i,j */
double Matrix::at(int i, int j) const
{
    return this->operator()(i,j);
}

/** Return the determinant of the matrix */
double Matrix::determinant() const
{
    //  | a b c |
    //  | d e f |
    //  | g h i |

    // det = (aei+bfg+cdh)-(ceg+bdi+afh)

    const double a = array[offset(0,0)];
    const double b = array[offset(0,1)];
    const double c = array[offset(0,2)];

    const double d = array[offset(1,0)];
    const double e = array[offset(1,1)];
    const double f = array[offset(1,2)];

    const double g = array[offset(2,0)];
    const double h = array[offset(2,1)];
    const double i = array[offset(2,2)];

    return (a*e*i + b*f*g + c*d*h) - (c*e*g + b*d*i + a*f*h);
}

/** Return a QString representation of the matrix */
QString Matrix::toString() const
{
    return QObject::tr("/ %1, %2, %3 \\\n| %4, %5, %6 |\n\\ %7, %8, %9 /")
                  .arg(xx()).arg(xy()).arg(xz())
                  .arg(yx()).arg(yy()).arg(yz())
                  .arg(zx()).arg(zy()).arg(zz());
}

/** Return the trace of the matrix */
Vector Matrix::trace() const
{
    return Vector(xx(),yy(),zz());
}

/** Return each column */
Vector Matrix::column0() const
{
    return Vector(xx(),yx(),zx());
}

/** Return each column */
Vector Matrix::column1() const
{
    return Vector(xy(),yy(),zy());
}

/** Return each column */
Vector Matrix::column2() const
{
    return Vector(xz(),yz(),zz());
}

/** Return each row */
Vector Matrix::row0() const
{
    return Vector(xx(),xy(),xz());
}

/** Return each row */
Vector Matrix::row1() const
{
    return Vector(yx(),yy(),yz());
}

/** Return each row */
Vector Matrix::row2() const
{
    return Vector(zx(),zy(),zz());
}

/** Return the transpose of the matrix */
Matrix Matrix::transpose() const
{
    return Matrix(xx(),yx(),zx(),xy(),yy(),zy(),xz(),yz(),zz());
}

/** Set the matrix to identity */
void Matrix::setToIdentity()
{
    for (int i=1; i<8; ++i)
    {
        array[i] = 0;
    }

    array[0] = 1;
    array[4] = 1;
    array[8] = 1;
}

/** Return whether or not this matrix is equal to the identity matrix */
bool Matrix::isIdentity() const
{
    return xx() == 1 and yy() == 1 and zz() == 1 and
           xy() == 0 and xz() == 0 and yx() == 0 and
           yz() == 0 and zx() == 0 and zy() == 0;

}

/** Copy assignment operator */
Matrix& Matrix::operator=(const Matrix &other)
{
    if (this != &other)
    {
        memcpy(array, other.array, 9*sizeof(double));
    }
    
    return *this;
}

bool Matrix::operator==(const Matrix& m) const
{
    for (int i=0; i<9; ++i)
    {
        if (array[i] != m.array[i])
            return false;
    }
    
    return true;
}

bool Matrix::operator!=(const Matrix& m) const
{
    return not this->operator==(m);
}

/** Return the identity matrix */
Matrix Matrix::identity()
{
    return Matrix();
}

/** Return the null matrix */
Matrix Matrix::zero()
{
    return Matrix(0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0);
}

Matrix& Matrix::operator*=(const Matrix &m)
{
    //  xx xy xz    0 1 2
    //  yx yy yz    3 4 5
    //  zx zy zz    6 7 8

    const double *a = array;
    const double *o = m.array;

    const double sxx = a[0]*o[0] + a[1]*o[3] + a[2]*o[6];
    const double sxy = a[0]*o[1] + a[1]*o[4] + a[2]*o[7];
    const double sxz = a[0]*o[2] + a[1]*o[5] + a[2]*o[8];

    const double syx = a[3]*o[0] + a[4]*o[3] + a[5]*o[6];
    const double syy = a[3]*o[1] + a[4]*o[4] + a[5]*o[7];
    const double syz = a[3]*o[2] + a[4]*o[5] + a[5]*o[8];

    const double szx = a[6]*o[0] + a[7]*o[3] + a[8]*o[6];
    const double szy = a[6]*o[1] + a[7]*o[4] + a[8]*o[7];
    const double szz = a[6]*o[2] + a[7]*o[5] + a[8]*o[8];

    array[0] = sxx;
    array[1] = sxy;
    array[2] = sxz;

    array[3] = syx;
    array[4] = syy;
    array[5] = syz;

    array[6] = szx;
    array[7] = szy;
    array[8] = szz;

    return *this;
}

const Matrix SIREMATHS_EXPORT SireMaths::operator*(const Matrix &m1, const Matrix &m2)
{
    Matrix ret(m1);
    ret *= m2;
    
    return ret;
}

Matrix& Matrix::operator+=(const Matrix &m)
{
    for (int i=0; i<9; ++i)
    {
        array[i] += m.array[i];
    }

    return *this;
}

Matrix& Matrix::operator-=(const Matrix &m)
{
    for (int i=0; i<9; ++i)
    {
        array[i] -= m.array[i];
    }

    return *this;
}

Matrix& Matrix::operator*=(double c)
{
    for (int i=0; i<9; ++i)
    {
        array[i] *= c;
    }

    return *this;
}

Matrix& Matrix::operator/=(double c)
{
    if (isZero(c))
        throw SireMaths::math_error(QObject::tr(
                            "Cannot divide a matrix by 0"),CODELOC);

    return this->operator*=( 1 / c );
}

const Matrix SIREMATHS_EXPORT SireMaths::operator+(const Matrix &m1, const Matrix &m2)
{
    Matrix ret(m1);
    ret += m2;
    
    return ret;
}

const Matrix SIREMATHS_EXPORT SireMaths::operator-(const Matrix &m1, const Matrix &m2)
{
    Matrix ret(m1);
    ret -= m2;
    return ret;
}

const Matrix SIREMATHS_EXPORT SireMaths::operator*(const Matrix &m, double c)
{
    Matrix ret(m);
    ret *= c;
    return ret;
}

const Matrix SIREMATHS_EXPORT SireMaths::operator*(double c, const Matrix &m)
{
    Matrix ret(m);
    ret *= c;
    return ret;
}

/** Return the inverse of this matrix. Throws a math_error if this
    matrix cannot be inverted. */
Matrix Matrix::inverse() const
{
    //calculate the determinant of the matrix
    double det = this->determinant();

    //if the determinant is zero then this matrix cannot be inverted
    if (isZero(det))
    {
        throw SireMaths::math_error(QObject::tr(
                    "Matrix '%1' cannot be inverted!").arg(toString()),CODELOC);
    }

    //take the inverse of the determinant
    det = double(1.0) / det;

    //form the elements of the inverse matrix
    Matrix inv;

    inv.array[0] = det * (yy()*zz() - zy()*yz());
    inv.array[1] = det * (xz()*zy() - zz()*xy());
    inv.array[2] = det * (xy()*yz() - yy()*xz());

    inv.array[3] = det * (zx()*yz() - yx()*zz());
    inv.array[4] = det * (xx()*zz() - zx()*xz());
    inv.array[5] = det * (xz()*yx() - yz()*xx());

    inv.array[6] = det * (yx()*zy() - zx()*yy());
    inv.array[7] = det * (xy()*zx() - xx()*zy());
    inv.array[8] = det * (xx()*yy() - yx()*xy());

    return inv;
}

/** Return whether or not this is a symmetric matrix */
bool Matrix::isSymmetric() const
{
    return ( zx() == xz() and yx() == xy() and yz() == zy() );
}

/** Ensure that this matrix is symmetric - this is done by copying the upper-right
    diagonal to the lower-left diagonal. Note that you should only really use this
    function on matricies that you know are symmetric, but may have lost some of
    their symmetry due to numerical rounding error */
void Matrix::enforceSymmetric()
{
    /**
         xx yx zx
         xy yy zy
         xz yz zz  **/

    array[2] = zx();
    array[1] = yx();
    array[5] = zy();
}

Matrix convertGSLMatrix(gsl_matrix *mat)
{
    return Matrix(mat);
}

/** Obtain the principal axes of this matrix. This can only be performed if this
    matrix is symmetric. You should only call this function for matricies that
    you know are symmetric, as this function will assume that the matrix is
    symmetric, and will thus only use the upper-right diagonal of values.
    The returned principal axes will be sorted from the highest eigenvalue to
    the lowest. */
Matrix Matrix::getPrincipalAxes() const
{
    //assume that this matrix is symmetrical - we will thus
    //only look at the values in the upper-right diagonal

    //now use the GNU Scientific Library to solve the eigenvalue
    //problem for this matrix
    double new_array[9];
    memcpy(new_array, array, 9*sizeof(double));

    gsl_matrix_view m = gsl_matrix_view_array(new_array,3,3);

    //allocate space for the resulting eigenvectors and eigenvalues
    gsl_vector *eig_val = gsl_vector_alloc(3);
    gsl_matrix *eig_vec = gsl_matrix_alloc(3,3);

    //now allocate some workspace for the calculation...
    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(3);

    //perform the calculation
    gsl_eigen_symmv(&m.matrix, eig_val, eig_vec, w);

    //free the space used by the calculation
    gsl_eigen_symmv_free(w);

    //now sort the eigenvectors from the smallest eigenvalue to
    //the largest
    gsl_eigen_symmv_sort (eig_val, eig_vec, GSL_EIGEN_SORT_ABS_ASC);

    //now copy the results back into a new Matrix
    Matrix ret = convertGSLMatrix(eig_vec);

    //free up the memory used by the GSL data...
    gsl_vector_free(eig_val);
    gsl_matrix_free(eig_vec);

    //finally, return the matrix of principal components
    return ret;
}

/** Return the single value decomposition of this matrix.

    This calculates the decomposition of this matrix
    into U S V^T, returning U, S and V in the tuple
*/
boost::tuple<Matrix,Matrix,Matrix> Matrix::singleValueDecomposition() const
{
    //use GSL
    gsl_matrix *A = 0;
    gsl_matrix *W = 0;
    gsl_vector *S = 0;
    
    try
    {
        //copy this matrix into A
        A = gsl_matrix_alloc(3, 3);
        
        for (int i=0; i<3; ++i)
        {
            for (int j=0; j<3; ++j)
            {
                gsl_matrix_set( A, i, j, array[offset(i,j)] );
            }
        }
        
        //create space to hold the matrices for single value decomposition
        S = gsl_vector_alloc(3);
        W = gsl_matrix_alloc(3, 3);
            
        // calculate single value decomposition of A into V S W^T
        int ok = gsl_linalg_SV_decomp_jacobi(A, W, S);

        if (ok != 0)
            throw SireMaths::domain_error( QObject::tr(
                    "Could not calculate the single value decomposition of %1.")
                        .arg(this->toString()), CODELOC );
        
        //copy out the results...
        Matrix a(A);
        Matrix w = Matrix(W).transpose();
        Matrix s( gsl_vector_get(S,0), 0, 0,
                  0, gsl_vector_get(S,1), 0,
                  0, 0, gsl_vector_get(S,2) );
        
        gsl_matrix_free(A);
        gsl_vector_free(S);
        gsl_matrix_free(W);
        
        return boost::tuple<Matrix,Matrix,Matrix>(a,s,w);
    }
    catch(...)
    {
        gsl_matrix_free(A);
        gsl_vector_free(S);
        gsl_matrix_free(W);
        throw;

        return boost::tuple<Matrix,Matrix,Matrix>();
    }
}

/** Return the single value decomposition of this matrix.

    This calculates the decomposition of this matrix
    into U S V^T, returning U, S and V in the tuple
*/
boost::tuple<Matrix,Matrix,Matrix> Matrix::svd() const
{
    return this->singleValueDecomposition();
}

/** Return the eigenvectors and eigenvalues of this matrix */
boost::tuple<Vector,Matrix> Matrix::diagonalise() const
{
    if (this->isSymmetric())
    {
        //we can use the quick eig3 code
        double A[3][3], V[3][3], d[3];
        
        A[0][0] = array[0];
        A[0][1] = array[1];
        A[0][2] = array[2];
        A[1][0] = array[3];
        A[1][1] = array[4];
        A[1][2] = array[5];
        A[2][0] = array[6];
        A[2][1] = array[7];
        A[2][2] = array[8];
        
        eigen_decomposition(A, V, d);
        
        return boost::tuple<Vector,Matrix>(
                    Vector(d[0], d[1], d[2]),
                    Matrix(V[0][0], V[1][0], V[2][0],
                           V[0][1], V[1][1], V[2][1],
                           V[0][2], V[1][2], V[2][2]) );
    }
    else
    {
        //we need to use BLAS - via NMatrix
        std::pair<NVector,NMatrix> eigs = NMatrix(*this).diagonalise();
        
        return boost::tuple<Vector,Matrix>( Vector(eigs.first), Matrix(eigs.second) );
    }
}

/** Return the covariance matrix of the passed arrays of points. This
    matches point p[i] against point q[i], and only calculates up to either
    the specified number of points, if n > 0, or to min(len(p),len(q)) */
Matrix Matrix::covariance(const QVector<Vector> &p, const QVector<Vector> &q, int n)
{
    if (n < 0 or n > qMin(p.count(), q.count()))
    {
        n = qMin(p.count(), q.count());
    }

    gsl_matrix *P = 0;
    gsl_matrix *Q = 0;
    gsl_matrix *C = 0;

    try
    {
        //convert the two vectors of points into GSL matrices
        P = gsl_matrix_alloc(n, 3);
        Q = gsl_matrix_alloc(n, 3);
        
        for (int i=0; i<n; ++i)
        {
            gsl_matrix_set(P, i, 0, p[i].x());
            gsl_matrix_set(P, i, 1, p[i].y());
            gsl_matrix_set(P, i, 2, p[i].z());
            
            gsl_matrix_set(Q, i, 0, q[i].x());
            gsl_matrix_set(Q, i, 1, q[i].y());
            gsl_matrix_set(Q, i, 2, q[i].z());
        }
        
        //create space to hold the covariance matrix
        C = gsl_matrix_alloc(3, 3);
        
        for (int i=0; i<3; ++i)
        {
            for (int j=0; j<3; ++j)
            {
                gsl_matrix_set(C, i, j, 0);
            }
        }
        
        //compute the covariance matrix P^T Q
        int ok = gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, P, Q, 0.0, C);

        if (ok != 0)
            throw SireMaths::domain_error( QObject::tr(
                    "Something went wrong with the dgemm in covariance!"), CODELOC );

        Matrix c(C);
        
        gsl_matrix_free(P);
        gsl_matrix_free(Q);
        gsl_matrix_free(C);
        
        return c;
    }
    catch(...)
    {
        gsl_matrix_free(P);
        gsl_matrix_free(Q);
        gsl_matrix_free(C);
        throw;
        return Matrix();
    }
}

const char* Matrix::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Matrix>() );
}
