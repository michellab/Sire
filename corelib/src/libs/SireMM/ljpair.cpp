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

#include "ljpair.h"
#include "ljparameter.h"

#include "SireMaths/maths.h"

#include "SireStream/datastream.h"

using namespace SireMM;
using namespace SireMaths;
using namespace SireStream;

static const RegisterMetaType<LJPair> r_ljpair(NO_ROOT);

/** Serialise to a binary data stream */
QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, const LJPair &ljpair)
{
    writeHeader(ds, r_ljpair, 1) << ljpair.sig << ljpair.eps;
    
    return ds;
}

/** Deserialise from a binary data stream */
QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, LJPair &ljpair)
{
    VersionID v = readHeader(ds, r_ljpair);
    
    if (v == 1)
    {
        ds >> ljpair.sig >> ljpair.eps;
    }
    else
        throw version_error(v, "1", r_ljpair, CODELOC);
    
    return ds;
}

/** Construct a dummy LJPair */
LJPair::LJPair() : sig(0), eps(0)
{}

/** Construct from an LJParameter */
LJPair::LJPair(const LJParameter &ljparam)
            : sig(ljparam.sigma()),
              eps(ljparam.epsilon())
{}

/** Copy constructor */
LJPair::LJPair(const LJPair &other)
            : sig(other.sig), eps(other.eps)
{}

/** Construct a LJPair that has the specified sigma and epsilon */
LJPair::LJPair(double sigma, double epsilon)
            : sig( std::abs(sigma) ), eps( std::abs(epsilon) )
{
    if ( SireMaths::isZero(sig) or SireMaths::isZero(eps) )
    {
        sig = 0;
        eps = 0;
    }
}

/** Construct a LJPair that is the combination of lj0 and lj1 
    using the passed specified combining rules */
LJPair::LJPair(const LJParameter &lj0, const LJParameter &lj1,
               LJParameterDB::CombiningRules rules) : sig(0), eps(0)
{
    switch(rules)
    {
        case LJParameterDB::ARITHMETIC:
            this->operator=( LJPair::arithmetic(lj0, lj1) );
            break;
        
        case LJParameterDB::GEOMETRIC:
            this->operator=( LJPair::geometric(lj0, lj1) );
            break;
    }
}

/** Destructor */
LJPair::~LJPair()
{}

/** Return a dummy CLJPair */
LJPair LJPair::dummy()
{
    return LJPair(0,0);
}

/** Return whether or not two LJPairs are equal */
bool LJPair::operator==(const LJPair &other) const
{
    return sig == other.sig and eps == other.eps;
}

/** Return whether or not two LJPairs are different */
bool LJPair::operator!=(const LJPair &other) const
{
    return not operator==(other);
}

/** Return whether or not this is a dummy LJ parameter */
bool LJPair::isDummy() const
{
    //we only need to compare sigma as this will be set to zero if 
    //epsilon is zero
    return sig == 0.0;
}

/** Return whether or not this parameter has non-zero LJ parameters */
bool LJPair::zeroLJ() const
{
    //we only need to compare sigma as this will be set to zero if 
    //epsilon is zero
    return sig == 0.0;
}

/** Return sqrt(sigma) */
double LJPair::sqrtSigma() const
{
    return std::sqrt(sig);
}

/** Return sqrt(epsilon) */
double LJPair::sqrtEpsilon() const
{
    return std::sqrt(eps);
}

/** Return the LJ 'A' parameter */
double LJPair::A() const
{
    return double(4) * epsilon() * SireMaths::pow_12( sigma() );
}

/** Return the LJ 'B' parameter */
double LJPair::B() const
{
    return double(4) * epsilon() * SireMaths::pow_6( sigma() );
}

/** Return the LJ 'rmin' parameter - this is the location of the minimum.

    rmin = 2^(1/6) * sigma
*/
double LJPair::rmin() const
{
    // 2.0 ^ (1/6) = 1.122462048309372981439932526193103967671
    return sigma() * double(1.122462048309372981439932526193103967671);
}

/** Return a string representation of the CLJ parameter */
QString LJPair::toString() const
{
    return QString("LJ( sigma = %1 A, epsilon = %2 kcal mol-1 )").arg(sigma()).arg(epsilon());
}

/** Return a LJ parameter that corresponds to the passed values of sigma and epsilon,

    E(r) = 4 epsilon [ (sigma/r)^12 - (sigma/r)^6 ]
*/
LJPair LJPair::fromSigmaAndEpsilon(double sigma, double epsilon)
{
    return LJPair(sigma,epsilon);
}

/** Return a LJ parameter that corresponds to the passed LJ parameters A and B,

    E(r) = A r-12  - B r-6
*/
LJPair LJPair::fromAAndB(double a, double b)
{
    // A = 4 epsilon sigma^12,  B = 4 epsilon sigma^6

    // epsilon = A / (4 sigma^12) = B / (4 sigma^6)
    //           A / B  = (4 sigma^12) / (4 sigma^6)
    //           sigma = ( A / B )^(1/6)

    //           epsilon = B / (4 (A/B) )
    //                   = B^2 / 4A

    return LJPair( std::pow( (a/b), (1.0/6.0) ),
                        (b*b) / (4.0*a) );
}

/** Return a LJ parameter that corresponds to the curve that has a minimum at
    rmin, and a well-depth of epsilon.

    E(r) = 4 epilson [ (sigma/r)^12 - (sigma/r)^6 ], where

    rmin = 2^(1/6) sigma
*/
LJPair LJPair::fromRMinAndEpsilon(double rmin, double epsilon)
{
    //sigma = rmin / 2^(1/6) - 2^(1/6) = 1.122462048309372981439932526193103967671

    return LJPair( rmin / double(1.122462048309372981439932526193103967671),
                        epsilon );
}

const char* LJPair::typeName()
{
    return QMetaType::typeName( qMetaTypeId<LJPair>() );
}
