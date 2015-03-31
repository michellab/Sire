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

#include "sire_blas.h"
#include "nvector.h"
#include "nmatrix.h"

#include "SireError/errors.h"

#ifndef SIRE_DISABLE_FORTRAN

typedef int BLAS_INT;

#include "sire_blas_f.h" // CONDITIONAL_INCLUDE

// Here are all of the prototypes of the BLAS functions used
// by sire_blas
extern "C"
{
    /** This is dgemm - see BLAS API for documentation */
    void SireDGEMM(const char *TRANSA, const char *TRANSB,
                   const BLAS_INT *M, const BLAS_INT *N,
                   const BLAS_INT *K,
                   const double *ALPHA, 
                   const double *A, const BLAS_INT *LDA,
                   const double *B, const BLAS_INT *LDB,
                   const double *BETA,
                   const double *C, const BLAS_INT *LDC);

    /** This is dgemv - see BLAS API for documentation */
    void SireDGEMV(const char *TRANSA, const BLAS_INT *M,
                   const BLAS_INT *N, const double *ALPHA,
                   const double *A, const BLAS_INT *LDA,
                   const double *X, const BLAS_INT *INCX,
                   const double *BETA,
                   double *Y, const BLAS_INT *INCY);

} // end of extern "C"

#endif // SIRE_DISABLE_FORTRAN

namespace SireMaths
{

/** Return the product of matrix 'A' with vector 'X'

    Returns A \times X
*/
NVector SIREMATHS_EXPORT dgemv(const NMatrix &A, const NVector &X)
{
    #ifdef SIRE_DISABLE_FORTRAN
    throw SireError::unsupported( QObject::tr(
            "dgemv not available as BLAS does not work with this version of Sire."),
                    CODELOC );

    return NVector();

    #else

    A.assertNColumns( X.nRows() );

    char TRANSA;
    BLAS_INT M, N, LDA, INCX, INCY;
    double ALPHA, BETA;
    
    if (A.isTransposed())
    {
        TRANSA = 'T';
        M = A.nColumns();
        N = A.nRows();
    }
    else
    {
        TRANSA = 'N';
        M = A.nRows();
        N = A.nColumns();
    }
    
    BOOST_ASSERT( X.nRows() == A.nColumns() );
    
    LDA = M;
    INCX = 1;
    INCY = 1;
    
    NVector Y( X.nRows() );
    
    ALPHA = 1;
    BETA = 0;
    
    ::SireDGEMV( &TRANSA, &M, &N, &ALPHA, A.constData(),
                 &LDA, X.constData(), &INCX, &BETA,
                 Y.data(), &INCY );

    return Y;
    
    #endif // SIRE_DISABLE_FORTRAN
}

NVector SIREMATHS_EXPORT dgemv(double alpha, const NMatrix &A, const NVector &X,
                               double beta, const NVector &Y)
{
    #ifdef SIRE_DISABLE_FORTRAN
    throw SireError::unsupported( QObject::tr(
            "dgemv not available as BLAS does not work with this version of Sire."),
                    CODELOC );

    return NVector();

    #else

    A.assertNColumns( X.nRows() );
    Y.assertNRows( X.nRows() );

    char TRANSA;
    BLAS_INT M, N, LDA, INCX, INCY;
    double ALPHA, BETA;
    
    if (A.isTransposed())
    {
        TRANSA = 'T';
        M = A.nColumns();
        N = A.nRows();
    }
    else
    {
        TRANSA = 'N';
        M = A.nRows();
        N = A.nColumns();
    }
    
    BOOST_ASSERT( X.nRows() == A.nColumns() );
    BOOST_ASSERT( Y.nRows() == X.nRows() );
    
    LDA = M;
    INCX = 1;
    INCY = 1;
    
    NVector RESULT( Y );
    
    ALPHA = alpha;
    BETA = beta;
    
    ::SireDGEMV( &TRANSA, &M, &N, &ALPHA, A.constData(),
                 &LDA, X.constData(), &INCX, &BETA,
                 RESULT.data(), &INCY );

    return RESULT;
    
    #endif // SIRE_DISABLE_FORTRAN
}


/** Return the matrix product of A \times B */
NMatrix SIREMATHS_EXPORT dgemm(const NMatrix &A, const NMatrix &B)
{
    #ifdef SIRE_DISABLE_FORTRAN
    throw SireError::unsupported( QObject::tr(
            "dgemm not available as BLAS does not work with this version of Sire."),
                    CODELOC );

    return NMatrix();

    #else

    B.assertNRows( A.nColumns() );

    char TRANSA, TRANSB;
    BLAS_INT M, N, K, LDA, LDB, LDC;
    
    if (A.isTransposed())
    {
        TRANSA = 'T';
        M = A.nColumns();
        K = A.nRows();
    }
    else
    {
        TRANSA = 'N';
        M = A.nRows();
        K = A.nColumns();
    }
    
    if (B.isTransposed())
    {
        TRANSB = 'T';
        N = B.nRows();
        BOOST_ASSERT( K == B.nColumns() );
    }
    else
    {
        TRANSB = 'N';
        N = B.nColumns();
        BOOST_ASSERT( K == B.nRows() );
    }

    LDA = K;
    LDB = N;
    LDC = M;
    
    NMatrix C(LDC, N);
    
    double ALPHA = 1;
    double BETA = 0;

    ::SireDGEMM( &TRANSA, &TRANSB, &M, &N, &K, 
                 &ALPHA, A.constData(), &LDA,
                 B.constData(), &LDB, &BETA, C.data(), &LDC );
               
    return C;
    
    #endif // SIRE_DISABLE_FORTRAN
}

/** Return the matrix product of 

    \alpha \times A \times B + \beta \times C
*/
NMatrix dgemm(double alpha, const NMatrix &A, const NMatrix &B, 
              double beta, const NMatrix &C)
{
    #ifdef SIRE_DISABLE_FORTRAN
    throw SireError::unsupported( QObject::tr(
            "dgemm not available as BLAS does not work with this version of Sire."),
                    CODELOC );

    return NMatrix();

    #else
 
    A.assertNColumns( B.nRows() );
    C.assertNRows( A.nRows() );
    C.assertNColumns( B.nColumns() );

    char TRANSA, TRANSB;
    BLAS_INT M, N, K, LDA, LDB, LDC;
    
    if (A.isTransposed())
    {
        TRANSA = 'T';
        M = A.nColumns();
        K = A.nRows();
    }
    else
    {
        TRANSA = 'N';
        M = A.nRows();
        K = A.nColumns();
    }
    
    if (B.isTransposed())
    {
        TRANSB = 'T';
        N = B.nRows();
        BOOST_ASSERT( K == B.nColumns() );
    }
    else
    {
        TRANSB = 'N';
        N = B.nColumns();
        BOOST_ASSERT( K == B.nRows() );
    }

    LDA = K;
    LDB = N;
    LDC = C.nRows();

    BOOST_ASSERT( LDC >= M );
    BOOST_ASSERT( C.nColumns() == N );

    NMatrix RESULT(C);

    double ALPHA = alpha;
    double BETA = beta;

    ::SireDGEMM( &TRANSA, &TRANSB, &M, &N, &K, 
                 &ALPHA, A.constData(), &LDA,
                 B.constData(), &LDB, &BETA, RESULT.data(), &LDC );
               
    return RESULT;
    
    #endif // SIRE_DISABLE_FORTRAN
}

} // end of namespace SireMaths
