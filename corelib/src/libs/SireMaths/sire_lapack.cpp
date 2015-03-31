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

#include "sire_lapack.h"
#include "nvector.h"
#include "nmatrix.h"

#include <cmath>

#include "SireError/errors.h"
#include "SireMaths/errors.h"

#ifndef SIRE_DISABLE_FORTRAN

typedef int LAPACK_INT;

#include "sire_lapack_f.h" // CONDITIONAL_INCLUDE

// Here are all of the prototypes of the LAPACK functions used
// by sire_lapack
extern "C"
{
    /** This is dsyev - see LAPACK API for documentation */
    void SireDSYEV(const char *JOBZ, const char *UPLO, 
                   const LAPACK_INT *N, double *A, 
                   const LAPACK_INT *LDA, 
                   double *W, double *WORK, 
                   const LAPACK_INT *LWORK, 
                   LAPACK_INT *INFO);

} // end of extern "C"

#endif // SIRE_DISABLE_FORTRAN

namespace SireMaths
{

std::pair<NVector,NMatrix> SIREMATHS_EXPORT dsyev(const NMatrix &A, bool upper)
{
    #ifdef SIRE_DISABLE_FORTRAN
    throw SireError::unsupported( QObject::tr(
            "dsyev not available as LAPACK does not work with this version of Sire."),
                    CODELOC );

    return std::pair<NVector,NMatrix>();

    #else

    if (A.isTransposed())
    {
        //we can only process a column-major ordered matrix...
        return dsyev( A.transpose().fullTranspose(), upper );
    }

    char JOBZ, UPLO;
    LAPACK_INT N, LDA, LWORK, INFO;
    
    JOBZ = 'V';
    
    if (upper)
        UPLO = 'U';
    else
        UPLO = 'L';
        
    N = A.nRows();
    
    BOOST_ASSERT( A.nColumns() == N );
    
    LDA = N;
    
    NVector EIGVAL(N);
    
    QVector<double> WORK( 5*N );
    LWORK = WORK.count();
    
    INFO = 0;
    
    NMatrix EIGVEC( A );
    
    ::SireDSYEV(&JOBZ, &UPLO, &N, EIGVEC.data(),
                &LDA, EIGVAL.data(), WORK.data(), &LWORK, &INFO);
              
    if (INFO != 0)
        throw SireMaths::domain_error( QObject::tr(
                "There was a problem running dsyev - INFO == %1. A ==\n%2.")
                    .arg(INFO).arg(A.toString()), CODELOC );
    
    return std::pair<NVector,NMatrix>(EIGVAL, EIGVEC);

    #endif // SIRE_DISABLE_FORTRAN
}

NVector SIREMATHS_EXPORT dsyev_eigenvalues(const NMatrix &A, bool upper)
{
    #ifdef SIRE_DISABLE_FORTRAN
    
    throw SireError::unsupported( QObject::tr(
            "dsyev_eigenvalues not available as LAPACK does not work "
            "with this version of Sire."),
                    CODELOC );

    return NVector();

    #else

    if (A.isTransposed())
    {
        //we can only process a column-major ordered matrix...
        return dsyev_eigenvalues( A.transpose().fullTranspose(), upper );
    }

    char JOBZ, UPLO;
    LAPACK_INT N, LDA, LWORK, INFO;
    
    JOBZ = 'N';
    
    if (upper)
        UPLO = 'U';
    else
        UPLO = 'L';
        
    N = A.nRows();
    
    BOOST_ASSERT( A.nColumns() == N );
    
    LDA = N;
    
    NVector EIGVAL(N);
    
    QVector<double> WORK( 5*N );
    LWORK = WORK.count();
    
    INFO = 0;
    
    NMatrix A_COPY( A );
    
    ::SireDSYEV(&JOBZ, &UPLO, &N, A_COPY.data(),
                &LDA, EIGVAL.data(), WORK.data(), &LWORK, &INFO);
              
    if (INFO != 0)
        throw SireMaths::domain_error( QObject::tr(
                "There was a problem running dsyev - INFO == %1. A ==\n%2.")
                    .arg(INFO).arg(A.toString()), CODELOC );
    
    return EIGVAL;
    
    #endif // SIRE_DISABLE_FORTRAN
}

} // end of namespace SireMaths
