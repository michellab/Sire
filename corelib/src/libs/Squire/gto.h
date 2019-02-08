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

#ifndef SQUIRE_GTO_H
#define SQUIRE_GTO_H

#include "orbital.h"

#include "SireMaths/vector.h"

SIRE_BEGIN_HEADER

namespace Squire
{
class GTO;
class GTOPair;
class S_GTO;
class SS_GTO;
}

SQUIRE_EXPORT QDataStream& operator<<(QDataStream&, const Squire::GTO&);
SQUIRE_EXPORT QDataStream& operator>>(QDataStream&, Squire::GTO&);

SQUIRE_EXPORT QDataStream& operator<<(QDataStream&, const Squire::GTOPair&);
SQUIRE_EXPORT QDataStream& operator>>(QDataStream&, Squire::GTOPair&);

namespace Squire
{

using SireMaths::Vector;

typedef SireBase::PropPtr<GTO> GTOPtr;

/** This is the base class of all single Gaussian Type Orbital shells (GTOs)
    (S_GTO (l==0), P_GTO (l==1), D_GTO (l==2), FPlus_GTO (l>=3))
    
    These orbitals are not used directly - rather they are combined into
    shell-pair orbitals, e.g. SS_GTO, SP_GTO etc. Integral functions then
    use these shell-pair orbital objects. The shell pair classes are;
    
    SS_GTO
    PP_GTO, PS_GTO
    DD_GTO, DP_GTO, DS_GTO
    FPlusFPlus_GTO, FPlusD_GTO, FPlusP_GTO, FPlusS_GTO
    
    An orbital shell contains all of the orbitals for a particular shell
    
    @author Christopher Woods
*/
class SQUIRE_EXPORT GTO : public OrbitalShell
{

friend SQUIRE_EXPORT QDataStream& ::operator<<(QDataStream&, const GTO&);
friend SQUIRE_EXPORT QDataStream& ::operator>>(QDataStream&, GTO&);

public:
    GTO();
    GTO(const GTO &other);
    
    virtual ~GTO();
    
    static const char* typeName();
    
    virtual GTO* clone() const=0;

	GTOPtr multiply(double coefficient) const;

    double alpha() const;
    double beta() const;
    
    double scale() const;

	bool isNull() const;

	static const GTO& null();

protected:
    GTO(double alpha, double scale);

    GTO& operator=(const GTO &other);
    
    bool operator==(const GTO &other) const;
    bool operator!=(const GTO &other) const;

private:
    /** The orbital exponent (alpha) */
    double alfa;
    
    /** The multiplication factor (including normalisation
        constant) */
    double scl;
};

/** This is the base class of all of the combined shell-pairs
    (e.g. SS_GTO, PS_GTO etc.)
    
    @author Christopher Woods
*/
class SQUIRE_EXPORT GTOPair : public ShellPair
{

friend SQUIRE_EXPORT QDataStream& ::operator<<(QDataStream&, const GTOPair&);
friend SQUIRE_EXPORT QDataStream& ::operator>>(QDataStream&, GTOPair&);

public:
    GTOPair();
    GTOPair(const GTOPair &other);
    
    virtual ~GTOPair();

    static const char* typeName();
    
    virtual GTOPair* clone() const=0;
    
    const Vector& P() const;
    const Vector& Q() const;
    
    double R2() const;
    
    double zeta() const;
    double eta() const;
    
    double xi() const;
    
    double K() const;
    double K_AB() const;
    double K_CD() const;

    double ss() const;
    
    double Q_AB() const;
    double Q_CD() const;
    
    bool isNull() const;
    
    static double T(const GTOPair &P, const GTOPair &Q);
    
    static double rho(const GTOPair &P, const GTOPair &Q);
    
    static double preFac(const GTOPair &P, const GTOPair &Q);
    
    static Vector W(const GTOPair &P, const GTOPair &Q);
    
    static const GTOPair& null();
    
protected:
    GTOPair(const Vector &A, const GTO &a,
            const Vector &B, const GTO &b);
    
    GTOPair& operator=(const GTOPair &other);
    
    bool operator==(const GTOPair &other) const;
    bool operator!=(const GTOPair &other) const;

	void setQ(double q);

private:
    /** The center of this combined SS shell pair - for 
        the orbitals a and b, with centers A and B and 
        exponents alpha and beta, we get;
        
        P = (alpha*A + beta*B) / (alpha + beta) */
    Vector _P;
    
    /** The distance squared between the two centers */
    double _R2;
    
    /** The zeta value of the combined gaussian. This is;
    
        zeta = alpha+beta
    */
    double _zeta;
    
    /** The xi value of the combined gaussian,
    
        xi = alpha*beta / (alpha + beta) 
    */
    double _xi;
    
    /** The K value for this combined gaussian. This is;
    
        K = sqrt(2) * pi^(5/4) * scl_a * scl_b * exp( (-alpha*beta/(alpha+beta))|A-B|^2 )
                / (alpha_beta)
    */
    double _K;
    
    /** The (s||s) overlap integral for this pair of orbitals */
    double _ss;
    
    /** The prescreening value - sqrt[ (ss|ss) ] - for this GTO pair */
    double _Q_AB;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Return the center of this shell-pair

    P = Q = (alpha*A + beta*B) / (alpha + beta)
*/
inline const Vector& GTOPair::P() const
{
    return _P;
}

/** Return the center of this shell-pair 

    Q = P = (alpha*A + beta*B) / (alpha + beta)
*/
inline const Vector& GTOPair::Q() const
{
    return GTOPair::P();
}

/** Return the distance-squared between the two orbitals
    that make up this shell-pair */
inline double GTOPair::R2() const
{
    return _R2;
}

/** Return the zeta value of this shell-pair 

    zeta = eta = alpha+beta
*/
inline double GTOPair::zeta() const
{
    return _zeta;
}

/** Return the eta value of this shell-pair 

    eta = zeta = alpha+beta
*/
inline double GTOPair::eta() const
{
    return GTOPair::zeta();
}

/** Return the xi value of this shell-pair.

    xi = (alpha * beta) / (alpha + beta)
*/
inline double GTOPair::xi() const
{
    return _xi;
}

/** Return the K value for this shell pair

    K = sqrt(2) * pi^(5/4) * scl_a * scl_b * exp( (-alpha*beta/(alpha+beta))|A-B|^2 )
                / (alpha_beta)

    (see Obara and Saika paper)
*/
inline double GTOPair::K() const
{
    return _K;
}

/** Synonym for GTOPair::K() */
inline double GTOPair::K_AB() const
{
    return GTOPair::K();
}

/** Synonym for GTOPair::K() */
inline double GTOPair::K_CD() const
{
    return GTOPair::K();
}

/** Return the value used for pre-screening this GTO pair. If this 
    pair is 'ab', then this is equal to max[ sqrt[ (a_i,b_j|a_i,b_j) ] ] */
inline double GTOPair::Q_AB() const
{
	return _Q_AB;
}

/** Synonym for GTOPair::Q() */
inline double GTOPair::Q_CD() const
{
	return GTOPair::Q_AB();
}

/** Return the (s||s) overlap integral for this pair of orbitals */
inline double GTOPair::ss() const
{
    return _ss;
}

#endif // SIRE_SKIP_INLINE_FUNCTIONS

}

SIRE_EXPOSE_CLASS( Squire::GTO )
SIRE_EXPOSE_CLASS( Squire::GTOPair )

SIRE_EXPOSE_PROPERTY( Squire::GTOPtr, Squire::GTO )

SIRE_END_HEADER

#endif
