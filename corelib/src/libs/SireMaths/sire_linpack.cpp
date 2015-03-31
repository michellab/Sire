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

#include "sire_linpack.h"
#include "nvector.h"
#include "nmatrix.h"

#include <cmath>

#include "SireError/errors.h"
#include "SireMaths/errors.h"

#ifndef SIRE_DISABLE_FORTRAN

typedef int LINPACK_INT;

#include "sire_linpack_f.h" // CONDITIONAL_INCLUDE

// Here are all of the prototypes of the LINPACK functions used
// by sire_linpack
extern "C"
{
    /** This is dgeco - see LINPACK API for documentation */
    void SireDGECO(double *A, const LINPACK_INT *LDA, 
                   const LINPACK_INT *N, LINPACK_INT *IPVT, 
                   double *RCOND, double *Z);

    /** This is dgedi - see LINPACK API for documentation */
    void SireDGEDI(double *A, const LINPACK_INT *LDA, 
                   const LINPACK_INT *N, const LINPACK_INT *IPVT, 
                   double *DET, double *WORK, const LINPACK_INT *JOB);

} // end of extern "C"

#endif // SIRE_DISABLE_FORTRAN

namespace SireMaths
{

std::pair< NMatrix,QVector<int> > SIREMATHS_EXPORT dgeco(const NMatrix &A)
{
    #ifdef SIRE_DISABLE_FORTRAN

    throw SireError::unsupported( QObject::tr(
            "dgeco not available as LINPACK does not work with this version of Sire."),
                    CODELOC );

    return std::pair< NMatrix,QVector<int> >();

    #else

    LINPACK_INT LDA, N;
    
    if (A.isTransposed())
        return dgeco( A.transpose().fullTranspose() );
        
    N = A.nRows();
    BOOST_ASSERT( A.nColumns() == N );
    
    LDA = N;
    
    QVector<LINPACK_INT> IPVT(N);
    
    double RCOND;
    QVector<double> Z(N);
    
    NMatrix A_OUT(A);
    
    ::SireDGECO(A_OUT.data(), &LDA, &N, IPVT.data(), &RCOND, Z.data());
    
    return std::pair< NMatrix,QVector<int> >(A_OUT, IPVT);
    
    #endif // SIRE_DISABLE_FORTRAN
}

NMatrix SIREMATHS_EXPORT dgedi_inverse(const NMatrix &A, const QVector<int> &IPVT)
{
    #ifdef SIRE_DISABLE_FORTRAN

    throw SireError::unsupported( QObject::tr(
            "dgedi_inverse not available as LINPACK does not work "
            "with this version of Sire."),
                    CODELOC );

    return NMatrix();

    #else

    LINPACK_INT LDA, N, JOB;

    N = A.nRows();
    BOOST_ASSERT( A.nColumns() == N );
    BOOST_ASSERT( IPVT.count() == N );
    
    LDA = N;
    
    QVector<double> WORK(N);
    
    JOB = 01;
    
    NMatrix A_OUT( A );
    
    double DET[2];
    
    ::SireDGEDI(A_OUT.data(), &LDA, &N, IPVT.constData(), &(DET[0]), 
                WORK.data(), &JOB);
              
    return A_OUT;

    #endif // SIRE_DISABLE_FORTRAN
}

double SIREMATHS_EXPORT dgedi_determinant(const NMatrix &A, const QVector<int> &IPVT)
{
    #ifdef SIRE_DISABLE_FORTRAN

    throw SireError::unsupported( QObject::tr(
            "dgedi_determinant not available as LINPACK does not work "
            "with this version of Sire."),
                    CODELOC );

    return 0;

    #else
    
    LINPACK_INT LDA, N, JOB;

    N = A.nRows();
    BOOST_ASSERT( A.nColumns() == N );
    BOOST_ASSERT( IPVT.count() == N );
    
    LDA = N;
    
    QVector<double> WORK(N);
    
    JOB = 10;
    
    NMatrix A_OUT( A );
    
    double DET[2];
    
    ::SireDGEDI(A_OUT.data(), &LDA, &N, IPVT.constData(), &(DET[0]), 
                WORK.data(), &JOB);
              
    // DET contains DET[0] * 10^DET[1]
    return DET[0] * std::pow(10,DET[1]);
    
    #endif // SIRE_DISABLE_FORTRAN
}

std::pair<double,NMatrix> SIREMATHS_EXPORT dgedi(const NMatrix &A, 
                                                 const QVector<int> &IPVT)
{
    #ifdef SIRE_DISABLE_FORTRAN
    
    throw SireError::unsupported( QObject::tr(
            "dgedi not available as LINPACK does not work with this version of Sire."),
                    CODELOC );

    return std::pair<double,NMatrix>();

    #else

    LINPACK_INT LDA, N, JOB;

    N = A.nRows();
    BOOST_ASSERT( A.nColumns() == N );
    BOOST_ASSERT( IPVT.count() == N );
    
    LDA = N;
    
    QVector<double> WORK(N);
    
    JOB = 11;
    
    NMatrix A_OUT( A );
    
    double DET[2];
    
    ::SireDGEDI(A_OUT.data(), &LDA, &N, IPVT.constData(), &(DET[0]), 
                WORK.data(), &JOB);
              
    // DET contains DET[0] * 10^DET[1]
    return std::pair<double,NMatrix>( DET[0] * std::pow(10,DET[1]), A_OUT );

    #endif // SIRE_DISABLE_FORTRAN
}

} // end of namespace SireMaths
