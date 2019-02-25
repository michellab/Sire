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

#include "ljparameter.h"
#include "ljpair.h"

#include "SireMaths/maths.h"

#include "SireStream/datastream.h"

using namespace SireStream;
using namespace SireUnits::Dimension;
using namespace SireMM;

static const RegisterMetaType<LJParameter> r_ljparam(NO_ROOT);

/** Serialise to a binary data stream */
QDataStream &operator<<(QDataStream &ds, const LJParameter &ljparam)
{
    writeHeader(ds, r_ljparam, 1) << ljparam.sqrtsig << ljparam.sqrteps;

    return ds;
}

/** Deserialise from a binary data stream */
QDataStream &operator>>(QDataStream &ds, LJParameter &ljparam)
{
    VersionID v = readHeader(ds, r_ljparam);

    if (v == 1)
    {
        ds >> ljparam.sqrtsig >> ljparam.sqrteps;
    }
    else
        throw version_error(v, "1", r_ljparam, CODELOC);

    return ds;
}

/** Construct a dummy LJ parameter */
LJParameter::LJParameter() : sqrtsig(0.0), sqrteps(0.0)
{}

/** Construct from an LJPair */
LJParameter::LJParameter(const LJPair &ljpair)
            : sqrtsig(ljpair.sqrtSigma()),
              sqrteps(ljpair.sqrtEpsilon())
{}

/** Copy constructor */
LJParameter::LJParameter(const LJParameter &other)
            : sqrtsig(other.sqrtsig), sqrteps(other.sqrteps)
{}

/** Construct a LJParameter that has the specified sigma and epsilon */
LJParameter::LJParameter(Length s, MolarEnergy e)
            : sqrtsig( std::sqrt(std::abs(s)) ), sqrteps( std::sqrt(std::abs(e)) )
{
    if ( SireMaths::isZero(sqrtsig) or SireMaths::isZero(sqrteps) )
    {
        sqrtsig = 0.0;
        sqrteps = 0.0;
    }
}

/** Destructor */
LJParameter::~LJParameter()
{}

/** Return the LJParameter that is this parameter combined with 'other'
    using arithmetic combining rules */
LJParameter LJParameter::combineArithmetic(const LJParameter &other) const
{
    LJPair pair(*this, other, LJParameterDB::ARITHMETIC);
    
    LJParameter ret;
    ret.sqrtsig = std::sqrt(pair.sigma());
    ret.sqrteps = std::sqrt(pair.epsilon());
    
    return ret;
}

/** Return the LJParameter that is this parameter combined with 'other'
    using geometric combining rules */
LJParameter LJParameter::combineGeometric(const LJParameter &other) const
{
    LJPair pair(*this, other, LJParameterDB::GEOMETRIC);
    
    LJParameter ret;
    ret.sqrtsig = std::sqrt(pair.sigma());
    ret.sqrteps = std::sqrt(pair.epsilon());
    
    return ret;
}

/** Return the LJParameter that is this parameter combined with 'other' according
    to the passed combining rules */
LJParameter LJParameter::combine(const LJParameter &other,
                                 LJParameter::CombiningRules rules) const
{
    switch(rules)
    {
        case LJParameter::ARITHMETIC:
            return this->combineArithmetic(other);
        case LJParameter::GEOMETRIC:
            return this->combineGeometric(other);
        default:
            return this->combineGeometric(other);
    }
}

/** Return a dummy CLJParameter */
LJParameter LJParameter::dummy()
{
    return LJParameter( Length(0), MolarEnergy(0) );
}

/** Return the LJ 'A' parameter */
double LJParameter::A() const
{
    return double(4.0) * double(epsilon()) 
                       * SireMaths::pow_12( double(sigma()) );
}

/** Return the LJ 'B' parameter */
double LJParameter::B() const
{
    return double(4.0) * double(epsilon()) 
                       * SireMaths::pow_6( double(sigma()) );
}

/** Return the LJ 'rmin' parameter - this is the location of the minimum.

    rmin = 2^(1/6) * sigma
*/
Length LJParameter::rmin() const
{
    // 2.0 ^ (1/6) = 1.122462048309372981439932526193103967671
    return sigma() * double(1.122462048309372981439932526193103967671);
}

/** Return a string representation of the CLJ parameter */
QString LJParameter::toString() const
{
    return QString("LJ( sigma = %1 A, epsilon = %2 kcal mol-1 )").arg(sigma()).arg(epsilon());
}

/** Return a LJ parameter that corresponds to the passed values of sigma and epsilon,

    E(r) = 4 epsilon [ (sigma/r)^12 - (sigma/r)^6 ]
*/
LJParameter LJParameter::fromSigmaAndEpsilon(Length sigma, MolarEnergy epsilon)
{
    return LJParameter(sigma,epsilon);
}

/** Return a LJ parameter that corresponds to the passed LJ parameters A and B,

    E(r) = A r-12  - B r-6
*/
LJParameter LJParameter::fromAAndB(double a, double b)
{
    // A = 4 epsilon sigma^12,  B = 4 epsilon sigma^6

    // epsilon = A / (4 sigma^12) = B / (4 sigma^6)
    //           A / B  = (4 sigma^12) / (4 sigma^6)
    //           sigma = ( A / B )^(1/6)

    //           epsilon = B / (4 (A/B) )
    //                   = B^2 / 4A

    return LJParameter( Length(std::pow( (a/b), (1.0/6.0)) ),
                        MolarEnergy((b*b) / (4.0*a)) );
}

/** Return a LJ parameter that corresponds to the curve that has a minimum at
    rmin, and a well-depth of epsilon.

    E(r) = 4 epilson [ (sigma/r)^12 - (sigma/r)^6 ], where

    rmin = 2^(1/6) sigma
*/
LJParameter LJParameter::fromRMinAndEpsilon(Length rmin, MolarEnergy epsilon)
{
    //sigma = rmin / 2^(1/6) - 2^(1/6) = 1.122462048309372981439932526193103967671

    return LJParameter( rmin / double(1.122462048309372981439932526193103967671),
                        epsilon );
}

const char* LJParameter::typeName()
{
    return QMetaType::typeName( qMetaTypeId<LJParameter>() );
}
