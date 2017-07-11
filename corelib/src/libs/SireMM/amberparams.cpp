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

#include "amberparams.h"

#include "SireMol/molecule.h"
#include "SireMol/partialmolecule.h"
#include "SireMol/bondid.h"
#include "SireMol/angleid.h"
#include "SireMol/dihedralid.h"
#include "SireMol/improperid.h"
#include "SireMol/atomidx.h"
#include "SireMol/connectivity.h"

#include "SireMM/twoatomfunctions.h"
#include "SireMM/threeatomfunctions.h"
#include "SireMM/fouratomfunctions.h"
#include "SireMM/cljnbpairs.h"

#include "SireCAS/expression.h"
#include "SireCAS/symbol.h"
#include "SireCAS/trigfuncs.h"
#include "SireCAS/values.h"
#include "SireCAS/sum.h"
#include "SireCAS/trigfuncs.h"

#include "SireBase/stringproperty.h"
#include "SireBase/parallel.h"

#include "SireError/errors.h"
#include "SireBase/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMol;
using namespace SireCAS;
using namespace SireMM;
using namespace SireBase;
using namespace SireStream;

///////////
/////////// Implementation of AmberBond
///////////

static const RegisterMetaType<AmberBond> r_bond(NO_ROOT);

QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, const AmberBond &bond)
{
    writeHeader(ds, r_bond, 1);

    ds << bond._k << bond._r0;
    return ds;
}

QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, AmberBond &bond)
{
    VersionID v = readHeader(ds, r_bond);
    
    if (v == 1)
    {
        ds >> bond._k >> bond._r0;
    }
    else
        throw version_error(v, "1", r_bond, CODELOC);

    return ds;
}

/** Construct with the passed bond constant and equilibrium bond length */
AmberBond::AmberBond(double k, double r0) : _k(k), _r0(r0)
{}

/** Construct from the passed expression */
AmberBond::AmberBond(const Expression &f, const Symbol &R) : _k(0), _r0(0)
{
    //expression should be of the form "k(r - r0)^2". We need to get the
    //factors of R
    const auto factors = f.expand(R);
    
    bool has_k = false;
    
    QStringList errors;
    
    double k = 0.0;
    double kr0_2 = 0.0;
    double kr0 = 0.0;
    
    for (const auto factor : factors)
    {
        if (factor.symbol() == R)
        {
            if (not factor.power().isConstant())
            {
                errors.append( QObject::tr("Power of R must be constant, not %1")
                                .arg(factor.power().toString()) );
                continue;
            }
        
            if (not factor.factor().isConstant())
            {
                errors.append( QObject::tr("The value of K in K (R - R0)^2 must be constant. "
                                "Here it is %1").arg(factor.factor().toString()) );
                continue;
            }

            double power = factor.power().evaluate(Values());
            
            if (power == 0.0)
            {
                //this is the constant
                kr0_2 += factor.factor().evaluate(Values());
            }
            else if (power == 1.0)
            {
                //this is the -kR0 term
                kr0 = factor.factor().evaluate(Values());
            }
            else if (power == 2.0)
            {
                //this is the R^2 term
                if (has_k)
                {
                    //we cannot have two R2 factors?
                    errors.append( QObject::tr("Cannot have two R^2 factors!") );
                    continue;
                }

                k = factor.factor().evaluate(Values());
                has_k = true;
            }
            else
            {
                errors.append( QObject::tr("Power of R^2 must equal 2.0, 1.0 or 0.0, not %1")
                                    .arg(power) );
                continue;
            }
        }
        else
        {
            errors.append( QObject::tr("Cannot have a factor that does not include R. %1")
                        .arg(factor.symbol().toString()) );
        }
    }
    
    if (kr0_2 < 0)
    {
        errors.append( QObject::tr("How can K R0^2 be negative? %1").arg(kr0_2) );
    }

    _k = k;
    _r0 = std::sqrt( kr0_2 / k );

    //kr0 should be equal to -2 k r0
    if ( std::abs(_k*_r0 + 0.5 * kr0) > 0.001 )
    {
        errors.append( QObject::tr("How can the power of R be %1. It should be 2 x %2 x %3 = %4.")
                        .arg(kr0).arg(_k).arg(_r0).arg(2*_k*_r0) );
    }
    
    
    if (not errors.isEmpty())
    {
        throw SireError::incompatible_error( QObject::tr(
                "Cannot extract an AmberBond with function K ( %1 - R0 )^2 from the "
                "expression %2, because\n%3")
                    .arg(R.toString()).arg(f.toString()).arg(errors.join("\n")), CODELOC );
    }
}

AmberBond::AmberBond(const AmberBond &other)
          : _k(other._k), _r0(other._r0)
{}

AmberBond::~AmberBond()
{}

double AmberBond::operator[](int i) const
{
    i = SireID::Index(i).map(2);
    
    if (i == 0)
        return _k;
    else
        return _r0;
}

AmberBond& AmberBond::operator=(const AmberBond &other)
{
    _k = other._k;
    _r0 = other._r0;
    return *this;
}

/** Comparison operator */
bool AmberBond::operator==(const AmberBond &other) const
{
    return _k == other._k and _r0 == other._r0;
}

/** Comparison operator */
bool AmberBond::operator!=(const AmberBond &other) const
{
    return not operator==(other);
}

/** Comparison operator */
bool AmberBond::operator<=(const AmberBond &other) const
{
    return (*this == other) or (*this < other);
}

/** Comparison operator */
bool AmberBond::operator>(const AmberBond &other) const
{
    return not (*this <= other);
}

/** Comparison operator */
bool AmberBond::operator>=(const AmberBond &other) const
{
    return not (*this < other);
}

const char* AmberBond::typeName()
{
    return QMetaType::typeName( qMetaTypeId<AmberBond>() );
}

const char* AmberBond::what() const
{
    return AmberBond::typeName();
}

/** Return the energy evaluated from this bond for the passed bond length */
double AmberBond::energy(double r) const
{
    return _k * SireMaths::pow_2(r - _r0);
}

/** Return an expression to evaluate the energy of this bond, using the passed
    symbol to represent the bond length */
Expression AmberBond::toExpression(const Symbol &R) const
{
    return _k * SireMaths::pow_2(R - _r0);
}

QString AmberBond::toString() const
{
    return QObject::tr("AmberBond( k = %1, r0 = %2 )").arg(_k).arg(_r0);
}

///////////
/////////// Implementation of AmberAngle
///////////

static const RegisterMetaType<AmberAngle> r_angle(NO_ROOT);

QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, const AmberAngle &angle)
{
    writeHeader(ds, r_angle, 1);
    ds << angle._k << angle._theta0;
    return ds;
}

QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, AmberAngle &angle)
{
    VersionID v = readHeader(ds, r_angle);
    
    if (v == 1)
    {
        ds >> angle._k >> angle._theta0;
    }
    else
        throw version_error(v, "1", r_angle, CODELOC);
    
    return ds;
}

AmberAngle::AmberAngle(double k, double theta0) : _k(k), _theta0(theta0)
{}

AmberAngle::AmberAngle(const Expression &f, const Symbol &theta) : _k(0), _theta0(0)
{
    //expression should be of the form "k(theta - theta0)^2". We need to get the
    //factors of theta
    const auto factors = f.expand(theta);
    
    bool has_k = false;
    
    QStringList errors;
    
    double k = 0.0;
    double ktheta0_2 = 0.0;
    double ktheta0 = 0.0;
    
    for (const auto factor : factors)
    {
        if (factor.symbol() == theta)
        {
            if (not factor.power().isConstant())
            {
                errors.append( QObject::tr("Power of theta must be constant, not %1")
                                .arg(factor.power().toString()) );
                continue;
            }
        
            if (not factor.factor().isConstant())
            {
                errors.append( QObject::tr("The value of K in K (theta - theta0)^2 must be "
                                "constant. Here it is %1").arg(factor.factor().toString()) );
                continue;
            }

            double power = factor.power().evaluate(Values());
            
            if (power == 0.0)
            {
                //this is the constant
                ktheta0_2 += factor.factor().evaluate(Values());
            }
            else if (power == 1.0)
            {
                //this is the -ktheta0 term
                ktheta0 = factor.factor().evaluate(Values());
            }
            else if (power == 2.0)
            {
                //this is the theta^2 term
                if (has_k)
                {
                    //we cannot have two R2 factors?
                    errors.append( QObject::tr("Cannot have two theta^2 factors!") );
                    continue;
                }

                k = factor.factor().evaluate(Values());
                has_k = true;
            }
            else
            {
                errors.append( QObject::tr("Power of theta^2 must equal 2.0, 1.0 or 0.0, not %1")
                                    .arg(power) );
                continue;
            }
        }
        else
        {
            errors.append( QObject::tr("Cannot have a factor that does not include theta. %1")
                        .arg(factor.symbol().toString()) );
        }
    }
    
    if (ktheta0_2 < 0)
    {
        errors.append( QObject::tr("How can K theta0^2 be negative? %1").arg(ktheta0_2) );
    }

    _k = k;
    _theta0 = std::sqrt( ktheta0_2 / k );

    //ktheta0 should be equal to -k theta0
    if ( std::abs(_k*_theta0 + 0.5*ktheta0) > 0.001 )
    {
        errors.append( QObject::tr(
                    "How can the power of theta be %1. It should be 2 x %2 x %3 = %4.")
                        .arg(ktheta0).arg(_k).arg(_theta0).arg(2*_k*_theta0) );
    }
    
    if (not errors.isEmpty())
    {
        throw SireError::incompatible_error( QObject::tr(
                "Cannot extract an AmberAngle with function K ( %1 - theta0 )^2 from the "
                "expression %2, because\n%3")
                    .arg(theta.toString()).arg(f.toString()).arg(errors.join("\n")), CODELOC );
    }
}

AmberAngle::AmberAngle(const AmberAngle &other)
           : _k(other._k), _theta0(other._theta0)
{}

AmberAngle::~AmberAngle()
{}

double AmberAngle::operator[](int i) const
{
    i = SireID::Index(i).map(2);
    
    if (i == 0)
        return _k;
    else
        return _theta0;
}

AmberAngle& AmberAngle::operator=(const AmberAngle &other)
{
    _k = other._k;
    _theta0 = other._theta0;
    return *this;
}

bool AmberAngle::operator==(const AmberAngle &other) const
{
    return _k == other._k and _theta0 == other._theta0;
}

bool AmberAngle::operator!=(const AmberAngle &other) const
{
    return not operator==(other);
}

/** Comparison operator */
bool AmberAngle::operator<=(const AmberAngle &other) const
{
    return (*this == other) or (*this < other);
}

/** Comparison operator */
bool AmberAngle::operator>(const AmberAngle &other) const
{
    return not (*this <= other);
}

/** Comparison operator */
bool AmberAngle::operator>=(const AmberAngle &other) const
{
    return not (*this < other);
}

const char* AmberAngle::typeName()
{
    return QMetaType::typeName( qMetaTypeId<AmberAngle>() );
}

const char* AmberAngle::what() const
{
    return AmberAngle::typeName();
}

double AmberAngle::energy(double theta) const
{
    return _k * SireMaths::pow_2(theta - _theta0);
}

Expression AmberAngle::toExpression(const Symbol &theta) const
{
    return _k * SireMaths::pow_2(theta - _theta0);
}

QString AmberAngle::toString() const
{
    return QObject::tr("AmberAngle( k = %1, theta0 = %2 )").arg(_k).arg(_theta0);
}

///////////
/////////// Implementation of AmberDihPart
///////////

static const RegisterMetaType<AmberDihPart> r_dihpart(NO_ROOT);

QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, const AmberDihPart &dih)
{
    writeHeader(ds, r_dihpart, 1);
    ds << dih._k << dih._periodicity << dih._phase;
    return ds;
}

QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, AmberDihPart &dih)
{
    VersionID v = readHeader(ds, r_dihpart);
    
    if (v == 1)
    {
        ds >> dih._k >> dih._periodicity >> dih._phase;
    }
    else
        throw version_error(v, "1", r_dihpart, CODELOC);
    
    return ds;
}

AmberDihPart::AmberDihPart(double k, double periodicity, double phase)
    : _k(k), _periodicity(periodicity), _phase(phase)
{}

AmberDihPart::AmberDihPart(const AmberDihPart &other)
             : _k(other._k), _periodicity(other._periodicity), _phase(other._phase)
{}

AmberDihPart::AmberDihPart::~AmberDihPart()
{}

double AmberDihPart::operator[](int i) const
{
    i = SireID::Index(i).map(3);
    
    if (i == 0)
        return _k;
    else if (i == 1)
        return _periodicity;
    else
        return _phase;
}

AmberDihPart& AmberDihPart::operator=(const AmberDihPart &other)
{
    _k = other._k;
    _periodicity = other._periodicity;
    _phase = other._phase;
    return *this;
}

bool AmberDihPart::operator==(const AmberDihPart &other) const
{
    return _k == other._k and _periodicity == other._periodicity
                 and _phase == other._phase;
}

bool AmberDihPart::operator!=(const AmberDihPart &other) const
{
    return not operator==(other);
}

/** Comparison operator */
bool AmberDihPart::operator<=(const AmberDihPart &other) const
{
    return (*this == other) or (*this < other);
}

/** Comparison operator */
bool AmberDihPart::operator>(const AmberDihPart &other) const
{
    return not (*this <= other);
}

/** Comparison operator */
bool AmberDihPart::operator>=(const AmberDihPart &other) const
{
    return not (*this < other);
}

const char* AmberDihPart::typeName()
{
    return QMetaType::typeName( qMetaTypeId<AmberDihPart>() );
}

const char* AmberDihPart::what() const
{
    return AmberDihPart::typeName();
}

double AmberDihPart::energy(double phi) const
{
    return _k * ( 1 + cos( (_periodicity * phi ) - _phase ) );
}

QString AmberDihPart::toString() const
{
    return QObject::tr("AmberDihPart( k = %1, periodicity = %2, phase = %3 )")
            .arg(_k).arg(_periodicity).arg(_phase);
}

///////////
/////////// Implementation of AmberDihedral
///////////

static const RegisterMetaType<AmberDihedral> r_dihedral(NO_ROOT);

QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, const AmberDihedral &dih)
{
    writeHeader(ds, r_dihedral, 1);
    ds << dih._parts;
    return ds;
}

QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, AmberDihedral &dih)
{
    VersionID v = readHeader(ds, r_dihedral);
    
    if (v == 1)
    {
        ds >> dih._parts;
    }
    else
        throw version_error(v, "1", r_dihedral, CODELOC);
    
    return ds;
}

AmberDihedral::AmberDihedral()
{}

AmberDihedral::AmberDihedral(AmberDihPart part)
{
    _parts = QVector<AmberDihPart>(1, part);
}

AmberDihedral::AmberDihedral(const Expression &f, const Symbol &phi)
{
    //this expression should be a sum of cos terms, plus constant terms
    QList<Expression> cos_terms;
    double constant = 0.0;
    
    QStringList errors;
    
    if (f.base().isA<Sum>())
    {
        for (const auto child : f.base().asA<Sum>().children())
        {
            if (child.isConstant())
            {
                constant += f.factor() * child.evaluate(Values());
            }
            else if (child.base().isA<Cos>())
            {
                if (f.factor() == 1)
                {
                    cos_terms.append(child);
                }
                else
                {
                    cos_terms.append( f.factor() * child );
                }
            }
            else
            {
                errors.append( QObject::tr( "Cannot interpret the dihedral expression "
                   "from '%1' as it should be a series of cosine terms involving %2.")
                        .arg(f.toString()).arg(phi.toString()) );
            }
        }
    }
    else
    {
        if (f.isConstant())
        {
            constant += f.evaluate(Values());
        }
        else if (f.base().isA<Cos>())
        {
            cos_terms.append(f);
        }
        else
        {
            errors.append( QObject::tr( "Cannot interpret the dihedral expression "
               "from '%1' as it should be a series of cosine terms involving %2.")
                    .arg(f.toString()).arg(phi.toString()) );
        }
    }

    //next extract all of the data from the cos terms
    QList< std::tuple<double,double,double> > terms;
    
    double check_constant = 0;
    
    for (const auto cos_term : cos_terms)
    {
        //term should be of the form 'k cos( periodicity * phi - phase )'
        double k = cos_term.factor();
        check_constant += k;
        
        double periodicity = 0.0;
        double phase = 0.0;
        
        const auto factors = cos_term.base().asA<Cos>().argument().expand(phi);
        
        bool ok = true;
        
        for (const auto factor : factors)
        {
            if (not factor.power().isConstant())
            {
                errors.append( QObject::tr("Power of phi must be constant, not %1")
                                .arg(factor.power().toString()) );
                ok = false;
                continue;
            }
        
            if (not factor.factor().isConstant())
            {
                errors.append( QObject::tr("The value of periodicity must be "
                                "constant. Here it is %1").arg(factor.factor().toString()) );
                ok = false;
                continue;
            }

            double power = factor.power().evaluate(Values());
            
            if (power == 0.0)
            {
                //this is the constant phase
                phase += factor.factor().evaluate(Values());
            }
            else if (power == 1.0)
            {
                //this is the periodicity * phi term
                periodicity = factor.factor().evaluate(Values());
            }
            else
            {
                errors.append( QObject::tr("Power of phi must equal 1.0 or 0.0, not %1")
                                    .arg(power) );
                ok = false;
                continue;
            }
        }
        
        if (ok)
        {
            terms.append( std::make_tuple(k, periodicity, phase) );
        }
    }
    
    //now sum up the individual values of 'k'. These should equal the constant, which
    //is the aggregate of the "k [ 1 + cos ]" terms
    if ( std::abs(constant - check_constant) > 0.001 )
    {
        errors.append( QObject::tr( "The set of constants should sum up to the same total, "
            "i.e. should equal %1. Instead they equal %2. This is weird.")
                .arg(constant).arg(check_constant) );
    }
    
    if (not errors.isEmpty())
    {
        throw SireError::incompatible_error( QObject::tr(
            "Cannot extract an Amber-format dihedral expression from '%1' as "
            "the expression must be a series of terms of type "
            "'k{ 1 + cos[ per %2 - phase ] }'. Errors include\n%3")
                .arg(f.toString()).arg(phi.toString()).arg(errors.join("\n")),
                    CODELOC );
    }

    //otherwise, add in all of the terms
    if (not terms.isEmpty())
    {
        _parts.reserve(terms.count());
        
        for (const auto term : terms)
        {
            //remember that the expression uses the negative of the phase ;-)
            _parts.append( AmberDihPart(std::get<0>(term), std::get<1>(term), -std::get<2>(term)) );
        }
    }
}

AmberDihedral::AmberDihedral(const AmberDihedral &other)
              : _parts(other._parts)
{}

AmberDihedral::~AmberDihedral()
{}

AmberDihedral& AmberDihedral::operator+=(const AmberDihPart &part)
{
    _parts.append(part);
    return *this;
}

AmberDihedral AmberDihedral::operator+(const AmberDihPart &part) const
{
    AmberDihedral ret(*this);
    ret += part;
    return *this;
}

AmberDihedral& AmberDihedral::operator=(const AmberDihedral &other)
{
    _parts = other._parts;
    return *this;
}

bool AmberDihedral::operator==(const AmberDihedral &other) const
{
    return _parts == other._parts;
}

bool AmberDihedral::operator!=(const AmberDihedral &other) const
{
    return not operator==(other);
}

const char* AmberDihedral::typeName()
{
    return QMetaType::typeName( qMetaTypeId<AmberDihedral>() );
}

const char* AmberDihedral::what() const
{
    return AmberDihedral::typeName();
}

AmberDihPart AmberDihedral::operator[](int i) const
{
    if (_parts.isEmpty())
    {
        //this is a zero dihedral
        i = SireID::Index(i).map(_parts.count());
        return AmberDihPart();
    }
    else
    {
        i = SireID::Index(i).map(_parts.count());
        return _parts[i];
    }
}

double AmberDihedral::energy(double phi) const
{
    double total = 0;
    for (int i=0; i<_parts.count(); ++i)
    {
        total += _parts.constData()[i].energy(phi);
    }
    return total;
}

Expression AmberDihedral::toExpression(const Symbol &phi) const
{
    Expression ret;
    
    for (auto part : _parts)
    {
        ret += part.k() * ( 1 + Cos( (part.periodicity() * phi ) - part.phase() ) );
    }

    return ret;
}

QString AmberDihedral::toString() const
{
    if (_parts.isEmpty())
    {
        return QObject::tr("AmberDihedral( 0 )");
    }

    QStringList s;
    for (int i=0; i<_parts.count(); ++i)
    {
        s.append( QObject::tr("k[%1] = %2, periodicity[%1] = %3, phase[%1] = %4")
                    .arg(i).arg(_parts[i].k()).arg(_parts[i].periodicity())
                    .arg(_parts[i].phase()) );
    }
    
    return QObject::tr("AmberDihedral( %1 )").arg(s.join(", "));
}

///////////
/////////// Implementation of AmberNB14
///////////

static const RegisterMetaType<AmberNB14> r_nb14(NO_ROOT);

QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, const AmberNB14 &nb)
{
    writeHeader(ds, r_nb14, 1);
    ds << nb._cscl << nb._ljscl;
    return ds;
}

QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, AmberNB14 &nb)
{
    VersionID v = readHeader(ds, r_nb14);
    
    if (v == 1)
    {
        ds >> nb._cscl >> nb._ljscl;
    }
    else
        throw version_error(v, "1", r_nb14, CODELOC);
    
    return ds;
}

AmberNB14::AmberNB14(double cscl, double ljscl) : _cscl(cscl), _ljscl(ljscl)
{}

AmberNB14::AmberNB14(const AmberNB14 &other)
          : _cscl(other._cscl), _ljscl(other._ljscl)
{}

AmberNB14::~AmberNB14()
{}

double AmberNB14::operator[](int i) const
{
    i = SireID::Index(i).map(2);
    
    if (i == 0)
        return _cscl;
    else
        return _ljscl;
}

AmberNB14& AmberNB14::operator=(const AmberNB14 &other)
{
    _cscl = other._cscl;
    _ljscl = other._ljscl;
    return *this;
}

/** Comparison operator */
bool AmberNB14::operator==(const AmberNB14 &other) const
{
    return _cscl == other._cscl and _ljscl == other._ljscl;
}

/** Comparison operator */
bool AmberNB14::operator!=(const AmberNB14 &other) const
{
    return not operator==(other);
}

/** Comparison operator */
bool AmberNB14::operator<(const AmberNB14 &other) const
{
    if (_cscl < other._cscl)
    {
        return true;
    }
    else if (_cscl == other._cscl)
    {
        return _ljscl < other._ljscl;
    }
    else
    {
        return false;
    }
}

/** Comparison operator */
bool AmberNB14::operator<=(const AmberNB14 &other) const
{
    return operator==(other) or operator<(other);
}

/** Comparison operator */
bool AmberNB14::operator>(const AmberNB14 &other) const
{
    return not operator<=(other);
}

/** Comparison operator */
bool AmberNB14::operator>=(const AmberNB14 &other) const
{
    return not operator<(other);
}

const char* AmberNB14::typeName()
{
    return QMetaType::typeName( qMetaTypeId<AmberNB14>() );
}

const char* AmberNB14::what() const
{
    return AmberNB14::typeName();
}

QString AmberNB14::toString() const
{
    return QObject::tr("AmberNB14( cscl = %1, ljscl = %2 )").arg(_cscl).arg(_ljscl);
}

/** Return the value converted to a CLJScaleFactor */
CLJScaleFactor AmberNB14::toScaleFactor() const
{
    return CLJScaleFactor(_cscl, _ljscl);
}

///////////
/////////// Implementation of AmberNBDihPart
///////////

static const RegisterMetaType<AmberNBDihPart> r_nbdihpart(NO_ROOT);

QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, const AmberNBDihPart &param)
{
    writeHeader(ds, r_nbdihpart, 1);
    ds << param.dih << param.nbscl;
    return ds;
}

QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, AmberNBDihPart &param)
{
    VersionID v = readHeader(ds, r_nbdihpart);
    
    if (v == 1)
    {
        ds >> param.dih >> param.nbscl;
    }
    else
        throw version_error(v, "1", r_nbdihpart, CODELOC);
    
    return ds;
}

/** Constructor */
AmberNBDihPart::AmberNBDihPart()
{}

/** Construct from the passed parameters */
AmberNBDihPart::AmberNBDihPart(const AmberDihPart &dihedral, const AmberNB14 &nb14)
               : dih(dihedral), nbscl(nb14)
{}

/** Copy constructor */
AmberNBDihPart::AmberNBDihPart(const AmberNBDihPart &other)
               : dih(other.dih), nbscl(other.nbscl)
{}

/** Destructor */
AmberNBDihPart::~AmberNBDihPart()
{}

/** Copy assignment operator */
AmberNBDihPart& AmberNBDihPart::operator=(const AmberNBDihPart &other)
{
    dih = other.dih;
    nbscl = other.nbscl;
    return *this;
}

/** Comparison operator */
bool AmberNBDihPart::operator==(const AmberNBDihPart &other) const
{
    return dih == other.dih and nbscl == other.nbscl;
}

/** Comparison operator */
bool AmberNBDihPart::operator!=(const AmberNBDihPart &other) const
{
    return not operator==(other);
}

/** Comparison operator */
bool AmberNBDihPart::operator<(const AmberNBDihPart &other) const
{
    if (nbscl < other.nbscl)
    {
        return true;
    }
    else if (nbscl == other.nbscl)
    {
        return dih < other.dih;
    }
    else
    {
        return false;
    }
}

/** Comparison operator */
bool AmberNBDihPart::operator<=(const AmberNBDihPart &other) const
{
    return this->operator==(other) or this->operator<(other);
}

/** Comparison operator */
bool AmberNBDihPart::operator>(const AmberNBDihPart &other) const
{
    return not this->operator<=(other);
}

/** Comparison operator */
bool AmberNBDihPart::operator>=(const AmberNBDihPart &other) const
{
    return not this->operator<(other);
}

const char* AmberNBDihPart::typeName()
{
    return QMetaType::typeName( qMetaTypeId<AmberNBDihPart>() );
}

const char* AmberNBDihPart::what() const
{
    return AmberNBDihPart::typeName();
}

QString AmberNBDihPart::toString() const
{
    return QObject::tr("AmberNBDihPart( param == %1, scl == %2 )")
                .arg(dih.toString()).arg(nbscl.toString());
}

///////////
/////////// Implementation of AmberParams
///////////

static const RegisterMetaType<AmberParams> r_amberparam;

/** Serialise to a binary datastream */
QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, const AmberParams &amberparam)
{
    writeHeader(ds, r_amberparam, 2);

    SharedDataStream sds(ds);

    sds << amberparam.molinfo
        << amberparam.amber_charges << amberparam.amber_ljs
        << amberparam.amber_masses << amberparam.amber_elements
        << amberparam.amber_types << amberparam.born_radii
        << amberparam.amber_screens << amberparam.amber_treechains
        << amberparam.exc_atoms
        << amberparam.amber_bonds
        << amberparam.amber_angles << amberparam.amber_dihedrals
        << amberparam.amber_impropers << amberparam.amber_nb14s
        << amberparam.radius_set
        << amberparam.propmap
        << static_cast<const MoleculeProperty&>(amberparam);

    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, AmberParams &amberparam)
{
    VersionID v = readHeader(ds, r_amberparam);
    
    if (v == 2)
    {
        SharedDataStream sds(ds);

        sds >> amberparam.molinfo
            >> amberparam.amber_charges >> amberparam.amber_ljs
            >> amberparam.amber_masses >> amberparam.amber_elements
            >> amberparam.amber_types >> amberparam.born_radii
            >> amberparam.amber_screens >> amberparam.amber_treechains
            >> amberparam.exc_atoms
            >> amberparam.amber_bonds
            >> amberparam.amber_angles >> amberparam.amber_dihedrals
            >> amberparam.amber_impropers >> amberparam.amber_nb14s
            >> amberparam.radius_set
            >> amberparam.propmap
            >> static_cast<MoleculeProperty&>(amberparam);
    }
    else
        throw version_error(v, "2", r_amberparam, CODELOC);
        
    return ds;
}

/** Null Constructor */
AmberParams::AmberParams() : ConcreteProperty<AmberParams,MoleculeProperty>()
{}

/** Constructor for the passed molecule*/
AmberParams::AmberParams(const MoleculeView &mol, const PropertyMap &map)
            : ConcreteProperty<AmberParams,MoleculeProperty>()
{
    const auto moldata = mol.data();
    
    const auto param_name = map["parameters"];
    
    //if possible, start from the existing parameters and update from there
    if (moldata.hasProperty(param_name))
    {
        const Property &param_prop = moldata.property(param_name);
        
        if (param_prop.isA<AmberParams>())
        {
            this->operator=(param_prop.asA<AmberParams>());
            
            if (propmap == map and this->isCompatibleWith(moldata.info()))
            {
                this->_pvt_updateFrom(moldata);
                return;
            }
        }
    }
    
    //otherwise construct this parameter from scratch
    this->operator=(AmberParams());
    
    molinfo = MoleculeInfo(moldata.info());
    propmap = map;
    this->_pvt_createFrom(moldata);
}

/** Constructor for the passed molecule*/
AmberParams::AmberParams(const MoleculeInfo &info)
            : ConcreteProperty<AmberParams,MoleculeProperty>(),
              molinfo(info)
{}

/** Constructor for the passed molecule*/
AmberParams::AmberParams(const MoleculeInfoData &info)
            : ConcreteProperty<AmberParams,MoleculeProperty>(),
              molinfo(info)
{}

/** Copy constructor */
AmberParams::AmberParams(const AmberParams &other)
            : ConcreteProperty<AmberParams,MoleculeProperty>(),
              molinfo(other.molinfo),
              amber_charges(other.amber_charges),amber_ljs(other.amber_ljs),
              amber_masses(other.amber_masses),amber_elements(other.amber_elements),
              amber_types(other.amber_types), born_radii(other.born_radii),
              amber_screens(other.amber_screens),
              amber_treechains(other.amber_treechains),
              exc_atoms(other.exc_atoms),
              amber_bonds(other.amber_bonds),
              amber_angles(other.amber_angles),
              amber_dihedrals(other.amber_dihedrals),
              amber_impropers(other.amber_impropers),
              amber_nb14s(other.amber_nb14s),
              radius_set(other.radius_set),
              propmap(other.propmap)
{}

/** Copy assignment operator */
AmberParams& AmberParams::operator=(const AmberParams &other)
{
    if (this != &other)
    {
        MoleculeProperty::operator=(other);
        molinfo = other.molinfo;
        amber_charges = other.amber_charges;
        amber_ljs = other.amber_ljs;
        amber_masses = other.amber_masses;
        amber_elements = other.amber_elements;
        amber_types = other.amber_types;
        born_radii = other.born_radii;
        amber_screens = other.amber_screens;
        amber_treechains = other.amber_treechains;
        exc_atoms = other.exc_atoms;
        amber_bonds = other.amber_bonds;
        amber_angles = other.amber_angles;
        amber_dihedrals = other.amber_dihedrals;
        amber_impropers = other.amber_impropers;
        amber_nb14s = other.amber_nb14s;
        radius_set = other.radius_set;
        propmap = other.propmap;
    }
  
    return *this;
}

/** Destructor */
AmberParams::~AmberParams()
{}

/** Comparison operator */
bool AmberParams::operator==(const AmberParams &other) const
{
  return (molinfo == other.molinfo
          and amber_charges == other.amber_charges
          and amber_ljs == other.amber_ljs
          and amber_masses == other.amber_masses
          and amber_elements == other.amber_elements
          and amber_types == other.amber_types
          and born_radii == other.born_radii
          and amber_screens == other.amber_screens
          and amber_treechains == other.amber_treechains
          and exc_atoms == other.exc_atoms
          and amber_bonds == other.amber_bonds
          and amber_angles == other.amber_angles
          and amber_dihedrals == other.amber_dihedrals
          and amber_impropers == other.amber_impropers
          and amber_nb14s == other.amber_nb14s
          and radius_set == other.radius_set
          and propmap == other.propmap);
}

/** Comparison operator */
bool AmberParams::operator!=(const AmberParams &other) const
{
    return not AmberParams::operator==(other);
}

/** Return the layout of the molecule whose flexibility is contained
    in this object */
MoleculeInfo AmberParams::info() const
{
    return molinfo;
}

/** Set the property map that should be used to find and update properties
    of the molecule */
void AmberParams::setPropertyMap(const PropertyMap &map)
{
    propmap = map;
}

/** Return the property map that is used to find and update properties
    of the molecule */
const PropertyMap& AmberParams::propertyMap() const
{
    return propmap;
}

/** Validate this set of parameters. This checks that all of the requirements
    for an Amber set of parameters are met, e.g. that all Atom indicies are 
    contiguous and in-order, and that all atoms contiguously fill all residues
    etc. This returns any errors as strings. An empty set of strings indicates
    that there are no errors */
QStringList AmberParams::validate() const
{
    return QStringList();
}

QString AmberParams::toString() const
{
    if (molinfo.nAtoms() == 0)
        return QObject::tr("AmberParams::null");

    return QObject::tr("AmberParams( nAtoms()=%6 nBonds=%1, nAngles=%2, nDihedrals=%3 "
                       "nImpropers=%4 n14s=%5 )")
                            .arg(amber_bonds.count())
                            .arg(amber_angles.count())
                            .arg(amber_dihedrals.count())
                            .arg(amber_impropers.count())
                            .arg(amber_nb14s.count())
                            .arg(molinfo.nAtoms());
}

/** Convert the passed BondID into AtomIdx IDs, sorted in index order */
BondID AmberParams::convert(const BondID &bond) const
{
    AtomIdx atom0 = info().atomIdx(bond.atom0());
    AtomIdx atom1 = info().atomIdx(bond.atom1());
    
    if (atom0.value() <= atom1.value())
        return BondID(atom0,atom1);
    else
        return BondID(atom1,atom0);
}

/** Convert the passed AngleID into AtomIdx IDs, sorted in index order */
AngleID AmberParams::convert(const AngleID &angle) const
{
    AtomIdx atom0 = info().atomIdx(angle.atom0());
    AtomIdx atom1 = info().atomIdx(angle.atom1());
    AtomIdx atom2 = info().atomIdx(angle.atom2());
    
    if (atom0.value() <= atom2.value())
        return AngleID(atom0,atom1,atom2);
    else
        return AngleID(atom2,atom1,atom0);
}

/** Convert the passed DihedralID into AtomIdx IDs, sorted in index order */
DihedralID AmberParams::convert(const DihedralID &dihedral) const
{
    AtomIdx atom0 = info().atomIdx(dihedral.atom0());
    AtomIdx atom1 = info().atomIdx(dihedral.atom1());
    AtomIdx atom2 = info().atomIdx(dihedral.atom2());
    AtomIdx atom3 = info().atomIdx(dihedral.atom3());
    
    if (atom0.value() < atom3.value())
        return DihedralID(atom0,atom1,atom2,atom3);
    else if (atom0.value() > atom3.value())
        return DihedralID(atom3,atom2,atom1,atom0);
    else if (atom1.value() <= atom2.value())
        return DihedralID(atom0,atom1,atom2,atom3);
    else
        return DihedralID(atom3,atom2,atom1,atom0);
}

/** Convert the passed ImproperID into AtomIdx IDs, sorted in index order */
ImproperID AmberParams::convert(const ImproperID &improper) const
{
    AtomIdx atom0 = info().atomIdx(improper.atom0());
    AtomIdx atom1 = info().atomIdx(improper.atom1());
    AtomIdx atom2 = info().atomIdx(improper.atom2());
    AtomIdx atom3 = info().atomIdx(improper.atom3());
    
    if (atom0.value() < atom3.value())
        return ImproperID(atom0,atom1,atom2,atom3);
    else if (atom0.value() > atom3.value())
        return ImproperID(atom3,atom2,atom1,atom0);
    else if (atom1.value() <= atom2.value())
        return ImproperID(atom0,atom1,atom2,atom3);
    else
        return ImproperID(atom3,atom2,atom1,atom0);
}

/** Return whether or not this flexibility is compatible with the molecule 
    whose info is in 'molinfo' */
bool AmberParams::isCompatibleWith(const SireMol::MoleculeInfoData &molinfo) const
{
    return info().UID() == molinfo.UID();
}

const char* AmberParams::typeName()
{
    return QMetaType::typeName(qMetaTypeId<AmberParams>());
}

/** Return the charges on the atoms */
AtomCharges AmberParams::charges() const
{
    return amber_charges;
}

/** Return the atom masses */
AtomMasses AmberParams::masses() const
{
    return amber_masses;
}

/** Return the atom elements */
AtomElements AmberParams::elements() const
{
    return amber_elements;
}

/** Return the atom LJ parameters */
AtomLJs AmberParams::ljs() const
{
    return amber_ljs;
}

/** Return all of the amber atom types */
AtomStringProperty AmberParams::amberTypes() const
{
    return amber_types;
}

/** Return all of the Born radii of the atoms */
AtomRadii AmberParams::gbRadii() const
{
    return born_radii;
}

/** Return all of the Born screening parameters for the atoms */
AtomFloatProperty AmberParams::gbScreening() const
{
    return amber_screens;
}

/** Return all of the Amber treechain classification for all of the atoms */
AtomStringProperty AmberParams::treeChains() const
{
    return amber_treechains;
}

/** Set the atom parameters for the specified atom to the provided values */
void AmberParams::add(const AtomID &atom,
                      SireUnits::Dimension::Charge charge,
                      SireUnits::Dimension::MolarMass mass,
                      const SireMol::Element &element,
                      const SireMM::LJParameter &ljparam,
                      const QString &amber_type,
                      SireUnits::Dimension::Length born_radius,
                      double screening_parameter,
                      const QString &treechain)
{
    CGAtomIdx idx = molinfo.cgAtomIdx(atom);
    
    if (amber_charges.isEmpty())
    {
        //set up the objects to hold these parameters
        amber_charges = AtomCharges(molinfo);
        amber_ljs = AtomLJs(molinfo);
        amber_masses = AtomMasses(molinfo);
        amber_elements = AtomElements(molinfo);
        amber_types = AtomStringProperty(molinfo);
        born_radii = AtomRadii(molinfo);
        amber_screens = AtomFloatProperty(molinfo);
        amber_treechains = AtomStringProperty(molinfo);
    }
    
    amber_charges.set(idx, charge);
    amber_ljs.set(idx, ljparam);
    amber_masses.set(idx, mass);
    amber_elements.set(idx, element);
    amber_types.set(idx, amber_type);
    born_radii.set(idx, born_radius);
    amber_screens.set(idx, screening_parameter);
    amber_treechains.set(idx, treechain);
}

/** Return the connectivity of the molecule implied by the
    the bonds */
Connectivity AmberParams::connectivity() const
{
    auto connectivity = Connectivity(molinfo).edit();
    
    for (auto it = amber_bonds.constBegin();
         it != amber_bonds.constEnd();
         ++it)
    {
        connectivity.connect( it.key().atom0(), it.key().atom1() );
    }
    
    return connectivity.commit();
}

/** Set the radius set used by LEAP to assign the Born radii 
    of the atoms. This is just a string that is used to label
    the radius set in the PRM file */
void AmberParams::setRadiusSet(const QString &rset)
{
    radius_set = rset;
}

/** Return the radius set used by LEAP to assign the Born radii */
QString AmberParams::radiusSet() const
{
    return radius_set;
}

/** Set the excluded atoms of the molecule. This should be a 
    CLJNBPairs with the value equal to 0 for atom0-atom1 pairs
    that are excluded, and 1 for atom0-atom1 pairs that are
    to be included in the non-bonded calculation */
void AmberParams::setExcludedAtoms(const CLJNBPairs &excluded_atoms)
{
    molinfo.assertCompatibleWith( excluded_atoms.info() );
    exc_atoms = excluded_atoms;
}

/** Return the excluded atoms of the molecule. The returned
    object has a matrix of all atom pairs, where the value
    is 0 for atom0-atom1 pairs that are to be excluded,
    and 1 for atom0-atom1 pairs that are to be included
    in the nonbonded calculation */
CLJNBPairs AmberParams::excludedAtoms() const
{
    if (exc_atoms.isEmpty())
    {
        if (molinfo.nAtoms() <= 3)
        {
            //everything is bonded, so scale factor is 0
            return CLJNBPairs(molinfo, CLJScaleFactor(0,0));
        }
        else
        {
            //nothing is explicitly excluded
            return CLJNBPairs(molinfo, CLJScaleFactor(1,1));
        }
    }
    else
        return exc_atoms;
}

/** Return the CLJ nonbonded 1-4 scale factors for the molecule */
CLJNBPairs AmberParams::cljScaleFactors() const
{
    //start from the set of excluded atoms
    CLJNBPairs nbpairs = this->excludedAtoms();
    
    //now add in all of the 1-4 nonbonded scale factors
    for (auto it = amber_nb14s.constBegin();
         it != amber_nb14s.constEnd();
         ++it)
    {
        nbpairs.set( it.key().atom0(), it.key().atom1(), it.value().toScaleFactor() );
    }
    
    return nbpairs;
}

void AmberParams::add(const BondID &bond, double k, double r0, bool includes_h)
{
    BondID b = convert(bond);
    amber_bonds.insert( this->convert(bond), qMakePair(AmberBond(k,r0),includes_h) );
}

void AmberParams::remove(const BondID &bond)
{
    amber_bonds.remove( this->convert(bond) );
}

AmberBond AmberParams::getParameter(const BondID &bond) const
{
    return amber_bonds.value(this->convert(bond)).first;
}

/** Return all of the bond parameters converted to a set of TwoAtomFunctions */
TwoAtomFunctions AmberParams::bondFunctions(const Symbol &R) const
{
    TwoAtomFunctions funcs(molinfo);
    
    for (auto it = amber_bonds.constBegin();
         it != amber_bonds.constEnd();
         ++it)
    {
        funcs.set( it.key(), it.value().first.toExpression(R) );
    }
    
    return funcs;
}

/** Return all of the bond parameters converted to a set of TwoAtomFunctions */
TwoAtomFunctions AmberParams::bondFunctions() const
{
    return bondFunctions( Symbol("r") );
}

void AmberParams::add(const AngleID &angle, double k, double theta0, bool includes_h)
{
    amber_angles.insert( this->convert(angle), qMakePair(AmberAngle(k,theta0),includes_h) );
}

void AmberParams::remove(const AngleID &angle)
{
    amber_angles.remove( this->convert(angle) );
}

AmberAngle AmberParams::getParameter(const AngleID &angle) const
{
    return amber_angles.value( this->convert(angle) ).first;
}

/** Return all of the angle parameters converted to a set of ThreeAtomFunctions */
ThreeAtomFunctions AmberParams::angleFunctions(const Symbol &THETA) const
{
    ThreeAtomFunctions funcs(molinfo);
    
    for (auto it = amber_angles.constBegin();
         it != amber_angles.constEnd();
         ++it)
    {
        funcs.set( it.key(), it.value().first.toExpression(THETA) );
    }
    
    return funcs;
}

/** Return all of the angle parameters converted to a set of ThreeAtomFunctions */
ThreeAtomFunctions AmberParams::angleFunctions() const
{
    return angleFunctions( Symbol("theta") );
}

void AmberParams::add(const DihedralID &dihedral, double k,
                      double periodicity, double phase, bool includes_h)
{
    //convert the dihedral into AtomIdx indicies
    DihedralID d = this->convert(dihedral);

    // If dihedral already exists, we will append parameters
    if (amber_dihedrals.contains(d))
    {
        amber_dihedrals[d].first += AmberDihPart(k, periodicity, phase);
    }
    else
    {
        amber_dihedrals.insert(d, qMakePair(AmberDihedral(AmberDihPart(k, periodicity, phase)),
                                            includes_h));
    }
}

void AmberParams::remove(const DihedralID &dihedral)
{
    amber_dihedrals.remove( this->convert(dihedral) );
}

AmberDihedral AmberParams::getParameter(const DihedralID &dihedral) const
{
    return amber_dihedrals.value( this->convert(dihedral) ).first;
}

/** Return all of the dihedral parameters converted to a set of FourAtomFunctions */
FourAtomFunctions AmberParams::dihedralFunctions(const Symbol &PHI) const
{
    FourAtomFunctions funcs(molinfo);
    
    for (auto it = amber_dihedrals.constBegin();
         it != amber_dihedrals.constEnd();
         ++it)
    {
        funcs.set( it.key(), it.value().first.toExpression(PHI) );
    }
    
    return funcs;
}

/** Return all of the dihedral parameters converted to a set of FourAtomFunctions */
FourAtomFunctions AmberParams::dihedralFunctions() const
{
    return dihedralFunctions( Symbol("phi") );
}

void AmberParams::add(const ImproperID &improper, double k,
                      double periodicity, double phase, bool includes_h)
{
    ImproperID imp = this->convert(improper);

    if (amber_impropers.contains(imp))
    {
        amber_impropers[imp].first += AmberDihPart(k, periodicity, phase);
    }
    else
    {
        amber_impropers.insert(imp, qMakePair(AmberDihedral(AmberDihPart(k, periodicity, phase)),
                                              includes_h));
    }
}

void AmberParams::remove(const ImproperID &improper)
{
    amber_impropers.remove( this->convert(improper) );
}

AmberDihedral AmberParams::getParameter(const ImproperID &improper) const
{
    return amber_impropers.value( this->convert(improper) ).first;
}

/** Return all of the improper parameters converted to a set of FourAtomFunctions */
FourAtomFunctions AmberParams::improperFunctions(const Symbol &PHI) const
{
    FourAtomFunctions funcs(molinfo);
    
    for (auto it = amber_impropers.constBegin();
         it != amber_impropers.constEnd();
         ++it)
    {
        funcs.set( it.key(), it.value().first.toExpression(PHI) );
    }
    
    return funcs;
}

/** Return all of the improper parameters converted to a set of FourAtomFunctions */
FourAtomFunctions AmberParams::improperFunctions() const
{
    return improperFunctions( Symbol("phi") );
}

void AmberParams::addNB14(const BondID &pair, double cscl, double ljscl)
{
    amber_nb14s.insert( this->convert(pair), AmberNB14(cscl,ljscl) );
}

void AmberParams::removeNB14(const BondID &pair)
{
    amber_nb14s.remove( this->convert(pair) );
}

AmberNB14 AmberParams::getNB14(const BondID &pair) const
{
    return amber_nb14s.value( this->convert(pair) );
}

/** Add the parameters from 'other' to this set */
AmberParams& AmberParams::operator+=(const AmberParams &other)
{
    if (not this->isCompatibleWith(other.info()) or propmap != other.propmap)
    {
        throw SireError::incompatible_error( QObject::tr(
                "Cannot combine Amber parameters, as the two sets are incompatible!"),
                    CODELOC );
    }

    if (not other.amber_charges.isEmpty())
    {
        //we overwrite these charges with 'other'
        amber_charges = other.amber_charges;
    }
    
    if (not other.exc_atoms.isEmpty())
    {
        //we overwrite our excluded atoms with 'other'
        exc_atoms = other.exc_atoms;
    }

    if (not other.amber_ljs.isEmpty())
    {
        //we overwrite these LJs with 'other'
        amber_ljs = other.amber_ljs;
    }

    if (not other.amber_masses.isEmpty())
    {
        //we overwrite these masses with 'other'
        amber_masses = other.amber_masses;
    }

    if (not other.amber_elements.isEmpty())
    {
        //we overwrite these elements with 'other'
        amber_elements = other.amber_elements;
    }

    if (not other.amber_types.isEmpty())
    {
        //we overwrite these types with 'other'
        amber_types = other.amber_types;
    }

    if (not other.born_radii.isEmpty())
    {
        //we overwrite these radii with 'other'
        born_radii = other.born_radii;
    }

    if (not other.amber_screens.isEmpty())
    {
        //we overwrite these screening parameters with 'other'
        amber_screens = other.amber_screens;
    }
    
    if (not other.amber_treechains.isEmpty())
    {
        //we overwrite these treechain classification with 'other'
        amber_treechains = other.amber_treechains;
    }

    if (amber_bonds.isEmpty())
    {
        amber_bonds = other.amber_bonds;
    }
    else if (not other.amber_bonds.isEmpty())
    {
        for (auto it = other.amber_bonds.constBegin(); it != other.amber_bonds.constEnd(); ++it)
        {
            amber_bonds.insert(it.key(), it.value());
        }
    }

    if (amber_angles.isEmpty())
    {
        amber_angles = other.amber_angles;
    }
    else if (not other.amber_angles.isEmpty())
    {
        for (auto it = other.amber_angles.constBegin(); it != other.amber_angles.constEnd(); ++it)
        {
            amber_angles.insert(it.key(), it.value());
        }
    }

    if (amber_dihedrals.isEmpty())
    {
        amber_dihedrals = other.amber_dihedrals;
    }
    else if (not other.amber_dihedrals.isEmpty())
    {
        for (auto it = other.amber_dihedrals.constBegin(); it != other.amber_dihedrals.constEnd();
             ++it)
        {
            amber_dihedrals.insert(it.key(), it.value());
        }
    }
    
    if (amber_impropers.isEmpty())
    {
        amber_impropers = other.amber_impropers;
    }
    else if (not other.amber_impropers.isEmpty())
    {
        for (auto it = other.amber_impropers.constBegin(); it != other.amber_impropers.constEnd();
             ++it)
        {
            amber_impropers.insert(it.key(), it.value());
        }
    }

    if (amber_nb14s.isEmpty())
    {
        amber_nb14s = other.amber_nb14s;
    }
    else if (not other.amber_nb14s.isEmpty())
    {
        for (auto it = other.amber_nb14s.constBegin(); it != other.amber_nb14s.constEnd(); ++it)
        {
            amber_nb14s.insert(it.key(), it.value());
        }
    }

    if (not other.radius_set.isEmpty())
    {
        //overwrite the radius set with other
        radius_set = other.radius_set;
    }

    return *this;
}

/** Return a combination of the two passed AmberParams */
AmberParams AmberParams::operator+(const AmberParams &other) const
{
    AmberParams ret(*this);
    
    ret += other;
    
    return ret;
}

/** Update these parameters from the contents of the passed molecule. This
    will only work if these parameters are compatible with this molecule */
void AmberParams::updateFrom(const MoleculeView &molview)
{
    this->assertCompatibleWith(molview);
    this->_pvt_updateFrom(molview.data());
}

/** Internal function used to grab the property, catching errors and signalling if
    the correct property has been found */
template<class T>
T getProperty(const PropertyName &prop, const MoleculeData &moldata, bool *found)
{
    if (moldata.hasProperty(prop))
    {
        const Property &p = moldata.property(prop);
        
        if (p.isA<T>())
        {
            *found = true;
            return p.asA<T>();
        }
    }
    
    *found = false;
    return T();
}

/** Internal function used to guess the masses of atoms based on their element */
void guessMasses(AtomMasses &masses, const AtomElements &elements, bool *has_masses)
{
    for (int i=0; i<elements.nCutGroups(); ++i)
    {
        const CGIdx cg(i);
    
        for (int j=0; i<elements.nAtoms(cg); ++j)
        {
            const CGAtomIdx idx( cg, Index(j) );
        
            masses.set(idx, elements[idx].mass());
        }
    }
    
    *has_masses = true;
}

/** Internal function used to guess the element of atoms based on their name */
AtomElements guessElements(const MoleculeInfoData &molinfo, bool *has_elements)
{
    AtomElements elements(molinfo);
    
    for (int i=0; i<elements.nCutGroups(); ++i)
    {
        const CGIdx cg(i);
    
        for (int j=0; j<elements.nAtoms(cg); ++j)
        {
            const CGAtomIdx idx( cg, Index(j) );
            
            elements.set(idx, Element::biologicalElement(molinfo.name(idx).value()));
        }
    }
    
    *has_elements = true;
    return elements;
}

/** Construct the hash of bonds */
void AmberParams::getAmberBondsFrom(const TwoAtomFunctions &funcs)
{
    // get the set of all bond functions
    const auto potentials = funcs.potentials();
    
    // create temporary space to hold the converted bonds
    QVector< std::tuple<BondID,AmberBond,bool> > bonds( potentials.count() );
    auto bonds_data = bonds.data();
    
    // convert each of these into an AmberBond
    tbb::parallel_for( tbb::blocked_range<int>(0,potentials.count()),
                       [&](const tbb::blocked_range<int> &r)
    {
        for (int i=r.begin(); i<r.end(); ++i)
        {
            const auto potential = potentials.constData()[i];
            
            //convert the atom IDs into a canonical form
            BondID bond = this->convert( BondID(potential.atom0(),potential.atom1()) );
            
            //does this bond involve hydrogen? - this relies on "AtomElements" being full
            bool contains_hydrogen = false;
            
            if (not amber_elements.isEmpty())
            {
                contains_hydrogen = (amber_elements.at(potential.atom0()).nProtons() < 2) or
                                    (amber_elements.at(potential.atom1()).nProtons() < 2);
            }
            
            bonds_data[i] = std::make_tuple(bond, AmberBond(potential.function(), Symbol("r")),
                                            contains_hydrogen);
        }
    });
    
    //finally add all of these into the amber_bonds hash
    amber_bonds.clear();
    amber_bonds.reserve(bonds.count());
    
    for (int i=0; i<bonds.count(); ++i)
    {
        amber_bonds.insert( std::get<0>(bonds_data[i]),
                            qMakePair(std::get<1>(bonds_data[i]),
                                      std::get<2>(bonds_data[i]) ) );
    }
}

/** Construct the hash of angles */
void AmberParams::getAmberAnglesFrom(const ThreeAtomFunctions &funcs)
{
    // get the set of all angle functions
    const auto potentials = funcs.potentials();
    
    // create temporary space to hold the converted angles
    QVector< std::tuple<AngleID,AmberAngle,bool> > angles( potentials.count() );
    auto angles_data = angles.data();
    
    // convert each of these into an AmberAngle
    tbb::parallel_for( tbb::blocked_range<int>(0,potentials.count()),
                       [&](const tbb::blocked_range<int> &r)
    {
        for (int i=r.begin(); i<r.end(); ++i)
        {
            const auto potential = potentials.constData()[i];
            
            //convert the atom IDs into a canonical form
            AngleID angle = this->convert( AngleID(potential.atom0(),
                                                   potential.atom1(),
                                                   potential.atom2()) );
            
            //does this angle involve hydrogen? - this relies on "AtomElements" being full
            bool contains_hydrogen = false;
            
            if (not amber_elements.isEmpty())
            {
                contains_hydrogen = (amber_elements.at(potential.atom0()).nProtons() < 2) or
                                    (amber_elements.at(potential.atom1()).nProtons() < 2) or
                                    (amber_elements.at(potential.atom2()).nProtons() < 2);
            }
            
            angles_data[i] = std::make_tuple(angle,
                                             AmberAngle(potential.function(), Symbol("theta")),
                                             contains_hydrogen);
        }
    });
    
    //finally add all of these into the amber_angles hash
    amber_angles.clear();
    amber_angles.reserve(angles.count());
    
    for (int i=0; i<angles.count(); ++i)
    {
        amber_angles.insert( std::get<0>(angles_data[i]),
                             qMakePair(std::get<1>(angles_data[i]),
                                       std::get<2>(angles_data[i]) ) );
    }
}

/** Construct the hash of dihedrals */
void AmberParams::getAmberDihedralsFrom(const FourAtomFunctions &funcs)
{
    // get the set of all dihedral functions
    const auto potentials = funcs.potentials();
    
    // create temporary space to hold the converted dihedrals
    QVector< std::tuple<DihedralID,AmberDihedral,bool> > dihedrals( potentials.count() );
    auto dihedrals_data = dihedrals.data();
    
    // convert each of these into an AmberDihedral
    tbb::parallel_for( tbb::blocked_range<int>(0,potentials.count()),
                       [&](const tbb::blocked_range<int> &r)
    {
        for (int i=r.begin(); i<r.end(); ++i)
        {
            const auto potential = potentials.constData()[i];
            
            //convert the atom IDs into a canonical form
            DihedralID dihedral = this->convert( DihedralID(potential.atom0(),
                                                            potential.atom1(),
                                                            potential.atom2(),
                                                            potential.atom3()) );
            
            //does this bond involve hydrogen? - this relies on "AtomElements" being full
            bool contains_hydrogen = false;
            
            if (not amber_elements.isEmpty())
            {
                contains_hydrogen = (amber_elements.at(potential.atom0()).nProtons() < 2) or
                                    (amber_elements.at(potential.atom1()).nProtons() < 2) or
                                    (amber_elements.at(potential.atom2()).nProtons() < 2) or
                                    (amber_elements.at(potential.atom3()).nProtons() < 2);
            }
            
            dihedrals_data[i] = std::make_tuple(dihedral,
                                                AmberDihedral(potential.function(), Symbol("phi")),
                                                contains_hydrogen);
        }
    });
    
    //finally add all of these into the amber_dihedrals hash
    amber_dihedrals.clear();
    amber_dihedrals.reserve(dihedrals.count());
    
    for (int i=0; i<dihedrals.count(); ++i)
    {
        amber_dihedrals.insert( std::get<0>(dihedrals_data[i]),
                                   qMakePair(std::get<1>(dihedrals_data[i]),
                                      std::get<2>(dihedrals_data[i]) ) );
    }
}

/** Construct the hash of impropers */
void AmberParams::getAmberImpropersFrom(const FourAtomFunctions &funcs)
{
    // get the set of all improper functions
    const auto potentials = funcs.potentials();
    
    // create temporary space to hold the converted dihedrals
    QVector< std::tuple<ImproperID,AmberDihedral,bool> > impropers( potentials.count() );
    auto impropers_data = impropers.data();
    
    // convert each of these into an AmberDihedral
    tbb::parallel_for( tbb::blocked_range<int>(0,potentials.count()),
                       [&](const tbb::blocked_range<int> &r)
    {
        for (int i=r.begin(); i<r.end(); ++i)
        {
            const auto potential = potentials.constData()[i];
            
            //convert the atom IDs into a canonical form
            ImproperID improper = this->convert( ImproperID(potential.atom0(),
                                                            potential.atom1(),
                                                            potential.atom2(),
                                                            potential.atom3()) );
            
            //does this bond involve hydrogen? - this relies on "AtomElements" being full
            bool contains_hydrogen = false;
            
            if (not amber_elements.isEmpty())
            {
                contains_hydrogen = (amber_elements.at(potential.atom0()).nProtons() < 2) or
                                    (amber_elements.at(potential.atom1()).nProtons() < 2) or
                                    (amber_elements.at(potential.atom2()).nProtons() < 2) or
                                    (amber_elements.at(potential.atom3()).nProtons() < 2);
            }
            
            impropers_data[i] = std::make_tuple(improper,
                                                AmberDihedral(potential.function(), Symbol("phi")),
                                                contains_hydrogen);
        }
    });
    
    //finally add all of these into the amber_dihedrals hash
    amber_impropers.clear();
    amber_impropers.reserve(impropers.count());
    
    for (int i=0; i<impropers.count(); ++i)
    {
        amber_impropers.insert( std::get<0>(impropers_data[i]),
                                   qMakePair(std::get<1>(impropers_data[i]),
                                      std::get<2>(impropers_data[i]) ) );
    }
}

/** Construct the excluded atom set and 14 NB parameters */
void AmberParams::getAmberNBsFrom(const CLJNBPairs &nbpairs,
                                  const FourAtomFunctions &dihedrals)
{
    //first, copy in the CLJNBPairs from the molecule
    exc_atoms = nbpairs;
    
    //now go through all dihedrals and get the 1-4 scale factors and
    //remove them from exc_atoms
    const auto potentials = dihedrals.potentials();
    
    //create new space to hold the 14 scale factors
    QHash<BondID,AmberNB14> new_nb14s;
    new_nb14s.reserve(potentials.count());
    
    for (int i=0; i<potentials.count(); ++i)
    {
        const auto potential = potentials.constData()[i];
        
        //convert the atom IDs into a canonical form
        BondID nb14pair = this->convert( BondID(potential.atom0(),
                                                potential.atom3()) );

        const auto molinfo = info();

        if (not new_nb14s.contains(nb14pair))
        {
            //extract the nb14 term from exc_atoms
            auto nbscl = nbpairs.get(nb14pair.atom0(),nb14pair.atom1());
            
            if (nbscl.coulomb() != 1.0 or nbscl.lj() != 1.0)
            {
                if (nbscl.coulomb() != 0.0 or nbscl.lj() != 0.0)
                {
                    //add them to the list of 14 scale factors
                    new_nb14s.insert(nb14pair, AmberNB14(nbscl.coulomb(),nbscl.lj()));
                
                    //and remove them from the excluded atoms list
                    exc_atoms.set(nb14pair.atom0(),nb14pair.atom1(), CLJScaleFactor(0));
                }
                else
                {
                    const auto tscl = nbpairs.get(nb14pair.atom0(),nb14pair.atom1());
                }
            }
            else
            {
                const auto tscl = nbpairs.get(nb14pair.atom0(),nb14pair.atom1());
            }
        }
    }
    
    amber_nb14s = new_nb14s;
}

/** Create this set of parameters from the passed object */
void AmberParams::_pvt_createFrom(const MoleculeData &moldata)
{
    //pull out all of the molecular properties
    const PropertyMap &map = propmap;
    
    //first, all of the atom-based properties
    bool has_charges, has_ljs, has_masses, has_elements, has_ambertypes;
    
    amber_charges = getProperty<AtomCharges>( map["charge"], moldata, &has_charges );
    amber_ljs = getProperty<AtomLJs>( map["LJ"], moldata, &has_ljs );
    amber_masses = getProperty<AtomMasses>( map["mass"], moldata, &has_masses );
    amber_elements = getProperty<AtomElements>( map["element"], moldata, &has_elements );
    amber_types = getProperty<AtomStringProperty>( map["ambertype"], moldata, &has_ambertypes );
    
    if (not has_elements)
    {
        //try to guess the elements from the names and/or masses
        amber_elements = guessElements(moldata.info(), &has_elements);
    }
    
    if (not has_masses)
    {
        //try to guess the masses from the elements
        if (has_elements)
        {
            amber_masses = AtomMasses(moldata.info());
            guessMasses(amber_masses, amber_elements, &has_masses);
        }
    }
    
    if (not (has_charges and has_ljs and has_masses and has_elements and has_ambertypes))
    {
        //it is not possible to create the parameter object if we don't have
        //these atom-based parameters
        throw SireBase::missing_property( QObject::tr(
                "Cannot create an AmberParams object for molecule %1 as it is missing "
                "necessary atom based properties: has_charges = %2, has_ljs = %3, "
                "has_masses = %4, has_elements = %5, has_ambertypes = %6.")
                    .arg( Molecule(moldata).toString() )
                    .arg(has_charges).arg(has_ljs).arg(has_masses)
                    .arg(has_elements).arg(has_ambertypes),
                        CODELOC );
    }
    
    // now see about the optional born radii and screening parameters
    bool has_radii, has_screening, has_treechains;
    
    born_radii = getProperty<AtomRadii>( map["gb_radii"], moldata, &has_radii );
    amber_screens = getProperty<AtomFloatProperty>( map["gb_screening"], moldata, &has_screening );
    amber_treechains = getProperty<AtomStringProperty>( map["treechain"],
                                                        moldata, &has_treechains );
    
    if (has_radii)
    {
        //see if there is a label for the source of the GB parameters
        bool has_source;
        
        radius_set = getProperty<StringProperty>( map["gb_radius_set"], moldata,
                                                  &has_source );
        
        if (not has_source)
        {
            radius_set = "unknown";
        }
    }
    else
    {
        radius_set = "unavailable";
    }

    // now lets get the bonded parameters (if they exist...)
    bool has_bonds, has_angles, has_dihedrals, has_impropers, has_nbpairs;
    
    const auto bonds = getProperty<TwoAtomFunctions>( map["bond"], moldata, &has_bonds );
    const auto angles = getProperty<ThreeAtomFunctions>( map["angle"], moldata, &has_angles );
    const auto dihedrals = getProperty<FourAtomFunctions>( map["dihedral"],
                                                           moldata, &has_dihedrals );
    const auto impropers = getProperty<FourAtomFunctions>( map["improper"],
                                                           moldata, &has_impropers );
    const auto nbpairs = getProperty<CLJNBPairs>( map["intrascale"],
                                                  moldata, &has_nbpairs );

    QVector< std::function<void()> > nb_functions;

    if (has_bonds)
    {
        nb_functions.append( [&](){ getAmberBondsFrom(bonds);} );
    }
    
    if (has_angles)
    {
        nb_functions.append( [&](){ getAmberAnglesFrom(angles);} );
    }
    
    if (has_dihedrals)
    {
        nb_functions.append( [&](){ getAmberDihedralsFrom(dihedrals);} );
    }
    
    if (has_impropers)
    {
        nb_functions.append( [&](){ getAmberImpropersFrom(impropers);} );
    }
    
    if (has_nbpairs)
    {
        nb_functions.append( [&](){ getAmberNBsFrom(nbpairs,dihedrals);} );
    }
    
    SireBase::parallel_invoke(nb_functions);
    
    //ensure that the resulting object is valid
    QStringList errors = this->validate();
    
    if (not errors.isEmpty())
    {
        throw SireError::io_error( QObject::tr(
                "Problem creating the AmberParams object for molecule %1 : %2. "
                "Errors include:\n%3")
                    .arg(moldata.name().value()).arg(moldata.number().value())
                    .arg(errors.join("\n")), CODELOC );
    }
}

/** Update this set of parameters from the passed object */
void AmberParams::_pvt_updateFrom(const MoleculeData &moldata)
{
    //for the moment we will just create everything from scratch.
    //However, one day we will optimise this and take existing
    //data that doesn't need to be regenerated.
    PropertyMap oldmap = propmap;
    const auto info = molinfo;
    
    this->operator=(AmberParams());

    propmap = oldmap;
    molinfo = info;

    this->_pvt_createFrom(moldata);
}
