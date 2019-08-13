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

#ifndef SIREMM_LJPAIR_H
#define SIREMM_LJPAIR_H

#include "ljparameter.h"
#include "ljparameterdb.h"

SIRE_BEGIN_HEADER

namespace SireMM
{
class LJPair;
}

class QDataStream;
QDataStream& operator<<(QDataStream&, const SireMM::LJPair&);
QDataStream& operator>>(QDataStream&, SireMM::LJPair&);

namespace SireMM
{

class LJPair;
class LJParameter;

/**
An LJPair holds a combined pair of Lennard Jones parameters
(as sigma and epsilon parameters)

@author Christopher Woods
*/
class SIREMM_EXPORT LJPair
{

friend SIREMM_EXPORT QDataStream& ::operator<<(QDataStream&, const LJPair&);
friend SIREMM_EXPORT QDataStream& ::operator>>(QDataStream&, LJPair&);

public:
    LJPair();
    LJPair(double sigma, double epsilon);
    LJPair(const LJParameter &ljparam);
    LJPair(const LJParameter &lj0, const LJParameter &lj1,
           LJParameterDB::CombiningRules combining_rules);       

    LJPair(const LJPair &other);

    ~LJPair();

    static const char* typeName();
    
    const char* what() const
    {
        return LJPair::typeName();
    }

    bool isDummy() const;
    bool zeroLJ() const;

    double sigma() const;
    double epsilon() const;
    double sqrtEpsilon() const;
    double sqrtSigma() const;

    double A() const;
    double B() const;

    double rmin() const;

    QString toString() const;

    bool operator==(const LJPair &other) const;
    bool operator!=(const LJPair &other) const;

    static LJPair dummy();
    static LJPair fromSigmaAndEpsilon(double sigma, double epsilon);
    static LJPair fromAAndB(double a, double b);
    static LJPair fromRMinAndEpsilon(double rmin, double epsilon);

    static LJPair arithmetic(const LJParameter &lj0, const LJParameter &lj1);
    static LJPair geometric(const LJParameter &lj0, const LJParameter &lj1);

private:
    /** the sigma parameter, in Angstroms */
    double sig;

    /** the epsilon parameter in kcal mol-1 */
    double eps;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Return the sigma value */
SIRE_ALWAYS_INLINE double LJPair::sigma() const
{
    return sig;
}

/** Return the epsilon value */
SIRE_ALWAYS_INLINE double LJPair::epsilon() const
{
    return eps;
}

/** Return an LJPair that represents the geometric combination of
    two LJParameters (sigma = sqrt(sig0*sig1), epsilon = sqrt(eps0*eps1) */
SIRE_ALWAYS_INLINE LJPair LJPair::geometric(const LJParameter &lj0, const LJParameter &lj1)
{
    LJPair ljpair;

    ljpair.sig = lj0.sqrtSigma() * lj1.sqrtSigma();
    ljpair.eps = lj0.sqrtEpsilon() * lj1.sqrtEpsilon();

    return ljpair;
}

/** Return an LJPair that represents the arithmetic combination
    of two LJParameters (Lorentz-Berthelot combining rules),
    sigma = 0.5(sig0+sig1), epsilon = sqrt(eps0*eps1) */
SIRE_ALWAYS_INLINE LJPair LJPair::arithmetic(const LJParameter &lj0, const LJParameter &lj1)
{
    LJPair ljpair;

    ljpair.sig = 0.5 * (lj0.sigma() + lj1.sigma());
    ljpair.eps = lj0.sqrtEpsilon() * lj1.sqrtEpsilon();

    return ljpair;
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

Q_DECLARE_TYPEINFO(SireMM::LJPair, Q_MOVABLE_TYPE);
Q_DECLARE_METATYPE(SireMM::LJPair);

SIRE_END_HEADER

#endif
