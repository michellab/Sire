/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2010  Christopher Woods
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

#ifndef SIREMM_LJFUNCTION_H
#define SIREMM_LJFUNCTION_H

#include "ljpair.h"

#include "SireMaths/distvector.h"
#include "SireMaths/maths.h"

#ifdef SIRE_USE_SSE
    #ifdef __SSE__
        #include <emmintrin.h>   // CONDITIONAL_INCLUDE
    #else
        #undef SIRE_USE_SSE
    #endif
#endif

SIRE_BEGIN_HEADER

/** This file contains the functions used to calculate LJ energies. */

namespace SireMM
{

using SireMaths::DistVector;
using SireMaths::Vector;

/** Calculate the LJ energy from the passed
    distance^2 and LJ parameter pair */
SIRE_ALWAYS_INLINE double calcLJEnergy(const double r2, const LJPair &ljpair)
{
    if (r2 == 0)
        return 0;
    else
    {
        const double inv_r2 = 1 / r2;

        double sig_over_r6 = SireMaths::pow_3(ljpair.sigma()*ljpair.sigma()*inv_r2);
        double sig_over_r12 = SireMaths::pow_2(sig_over_r6);

        return 4 * ljpair.epsilon() * (sig_over_r12 - sig_over_r6);
    }
}

/** Calculate the LJ force from the passed distance vector and LJ parameter pair */
SIRE_ALWAYS_INLINE Vector calcLJForce(const DistVector &r, const LJPair &ljpair)
{
    if (r.length() == 0)
        return Vector(0);
    else
    {
        const double inv_r = 1 / r.length();

        const double sig_over_r6 = SireMaths::pow_6(ljpair.sigma()*inv_r);
        const double sig_over_r7 = sig_over_r6 * inv_r;
        double sig_over_r13 = SireMaths::pow_2(sig_over_r6) * inv_r;

        return (4 * ljpair.epsilon() 
                    * (6 * sig_over_r7 - 12*sig_over_r13)) * r.direction();
    }
}

/** Calculate the LJ force from the passed distance vector and LJ parameter pair,
    together with the LJ scaling factor and differential of the scaling factor
      - this is for feathered energies */
SIRE_ALWAYS_INLINE Vector calcLJForce(const DistVector &r, const LJPair &ljpair,
                          const double scl, const double dscl_dr)
{
    if (r.length() == 0)
        return Vector(0);
    else
    {
        const double inv_r = 1 / r.length();

        const double sig_over_r6 = SireMaths::pow_6(ljpair.sigma()*inv_r);
        const double sig_over_r7 = sig_over_r6 * inv_r;
        const double sig_over_r12 = SireMaths::pow_2(sig_over_r6);
        const double sig_over_r13 = sig_over_r12 * inv_r;

        //calculate the energy
        const double ljnrg = 4 * ljpair.epsilon() *
                              (sig_over_r12 - sig_over_r6);

        return (scl * 4 * ljpair.epsilon() * 
                    (6.0*sig_over_r7 - 12.0*sig_over_r13))
                       * r.direction()
                    
                     + (ljnrg * dscl_dr);
    }
}

#ifdef SIRE_USE_SSE

    /** SSE version of calcLJEnergy */
    SIRE_ALWAYS_INLINE __m128d calcLJEnergy(const double r0_2, const double r1_2,
                                const LJPair &lj0pair, const LJPair &lj1pair)
    {
        __m128d sse_nrg;
    
        if (r0_2 == 0)
        {
            sse_nrg = _mm_set_pd( 0, SireMM::calcLJEnergy(r1_2, lj1pair) );
        }
        else if (r1_2 == 0)
        {
            sse_nrg = _mm_set_pd( SireMM::calcLJEnergy(r0_2, lj0pair), 0 );
        }
        else
        {
            const __m128d sse_one = { 1.0, 1.0 };
            const __m128d sse_four = { 4.0, 4.0 };
        
            const __m128d sse_r2 = _mm_set_pd( r0_2, r1_2 );
                               
            const __m128d sse_sig = _mm_set_pd( lj0pair.sigma(), lj1pair.sigma() );
            const __m128d sse_eps = _mm_set_pd( lj0pair.epsilon(), 
                                                lj1pair.epsilon() );
            
            const __m128d sse_inv_r2 = _mm_div_pd(sse_one, sse_r2);
            
            //calculate (sigma/r)^6 and (sigma/r)^12
            const __m128d sse_sig2 = _mm_mul_pd(sse_sig,sse_sig);
            const __m128d sse_sig_over_r2 = _mm_mul_pd(sse_sig2, sse_inv_r2);
                                         
            __m128d sse_sig_over_r6 = _mm_mul_pd(sse_sig_over_r2,
                                                 sse_sig_over_r2);
                                                    
            sse_sig_over_r6 = _mm_mul_pd(sse_sig_over_r6,
                                         sse_sig_over_r2);
                                         
            const __m128d sse_sig_over_r12 = _mm_mul_pd(sse_sig_over_r6,
                                                        sse_sig_over_r6);
                                  
            __m128d sig12_sig6 = _mm_sub_pd(sse_sig_over_r12, sse_sig_over_r6);
            __m128d four_eps = _mm_mul_pd(sse_four, sse_eps);
            
            sse_nrg = _mm_mul_pd( four_eps, sig12_sig6 );
        }
        
        return sse_nrg;
    }

#endif //SIRE_USE_SSE

}

SIRE_END_HEADER

#endif
