/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2017  Christopher Woods
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

#include "gromacsparams.h"

#include "SireID/index.h"

#include "SireError/errors.h"

#include "SireUnits/units.h"

#include "SireCAS/exp.h"
#include "SireCAS/conditional.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMM;
using namespace SireMol;
using namespace SireCAS;
using namespace SireID;
using namespace SireUnits;
using namespace SireUnits::Dimension;
using namespace SireStream;

/////////
///////// Implementation of GromacsAtomType
/////////

static const RegisterMetaType<GromacsAtomType> r_atomtype(NO_ROOT);

QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, const GromacsAtomType &typ)
{
    writeHeader(ds, r_atomtype, 1);
    
    SharedDataStream sds(ds);
    
    sds << typ._typ << typ._mass.to(g_per_mol)
        << typ._chg.to(mod_electron) << typ.particleTypeString()
        << typ._lj << typ._elem;
    
    return ds;
}

QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, GromacsAtomType &typ)
{
    VersionID v = readHeader(ds, r_atomtype);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        double chg, mass;
        QString ptyp;
        
        sds >> typ._typ >> mass >> chg >> ptyp >> typ._lj >> typ._elem;
        
        typ._mass = mass * g_per_mol;
        typ._chg = chg * mod_electron;
        typ._ptyp = GromacsAtomType::toParticleType(ptyp);
    }
    else
        throw version_error(v, "1", r_atomtype, CODELOC);
    
    return ds;
}

/** Assert that this object is sane */
void GromacsAtomType::assertSane() const
{}

/** Null constructor */
GromacsAtomType::GromacsAtomType()
                : _typ(), _mass(0), _chg(0),
                  _ptyp(GromacsAtomType::UNKNOWN_TYPE),
                  _elem(0)
{}

/** Construct passing in all parameters */
GromacsAtomType::GromacsAtomType(QString atom_type,
                                 SireUnits::Dimension::MolarMass mass,
                                 SireUnits::Dimension::Charge charge,
                                 PARTICLE_TYPE particle_type,
                                 const LJParameter &ljparam,
                                 const Element &element)
                : _typ(atom_type), _mass(mass), _chg(charge),
                  _lj(ljparam), _ptyp(particle_type), _elem(element)
{
    assertSane();
}

/** Construct, specifying only the mass */
GromacsAtomType::GromacsAtomType(QString atom_type, SireUnits::Dimension::MolarMass mass)
                : _typ(atom_type), _mass(mass), _chg(0),
                  _ptyp(GromacsAtomType::UNKNOWN_TYPE),
                  _elem(0)
{}

/** Copy constructor */
GromacsAtomType::GromacsAtomType(const GromacsAtomType &other)
                : _typ(other._typ), _mass(other._mass), _chg(other._chg),
                  _lj(other._lj), _ptyp(other._ptyp), _elem(other._elem)
{}

/** Destructor */
GromacsAtomType::~GromacsAtomType()
{}

/** Copy assignment operator */
GromacsAtomType& GromacsAtomType::operator=(const GromacsAtomType &other)
{
    if (this != &other)
    {
        _typ = other._typ;
        _mass = other._mass;
        _chg = other._chg;
        _ptyp = other._ptyp;
        _lj = other._lj;
        _elem = other._elem;
    }
    
    return *this;
}

/** Comparison operator */
bool GromacsAtomType::operator==(const GromacsAtomType &other) const
{
    return _typ == other._typ and _mass == other._mass and _chg == other._chg and
           _ptyp == other._ptyp and _lj == other._lj and _elem == other._elem;
}

/** Comparison operator */
bool GromacsAtomType::operator!=(const GromacsAtomType &other) const
{
    return not operator==(other);
}

const char* GromacsAtomType::typeName()
{
    return QMetaType::typeName( qMetaTypeId<GromacsAtomType>() );
}

const char* GromacsAtomType::what() const
{
    return GromacsAtomType::typeName();
}

/** Return a string representation of this object */
QString GromacsAtomType::toString() const
{
    if (hasMassOnly())
    {
        return QObject::tr("GromacsAtomType( atomType() = %1, mass() = %2 g mol-1 )")
                .arg(atomType()).arg(mass().to(g_per_mol));
    }
    else
    {
        return QObject::tr("GromacsAtomType( atomType() = %1, mass() = %2 g mol-1, "
                              "charge() = %3 |e|, "
                              "particleType() = %4, ljParameter() = %5, element() = %6 )")
                              .arg(atomType()).arg(mass().to(g_per_mol))
                              .arg(charge().to(mod_electron))
                              .arg(particleTypeString())
                              .arg(ljParameter().toString())
                              .arg(element().toString());
    }
}

/** Return the single letter that represents the particle type */
QString GromacsAtomType::particleTypeLetter() const
{
    switch(_ptyp)
    {
    case ATOM:
        return "A";
    case SHELL:
        return "S";
    case VIRTUAL:
        return "V";
    default:
        throw SireError::unknown_type( QObject::tr(
            "The particle type for Gromacs atom type '%1' is unknown!")
                .arg(this->toString()), CODELOC );
    }
}

/** Return a string version of the particle type */
QString GromacsAtomType::particleTypeString() const
{
    switch(_ptyp)
    {
    case ATOM:
        return "atom";
    case SHELL:
        return "shell";
    case VIRTUAL:
        return "virtual";
    default:
        return "unknown";
    }
}

/** Convert the passed string to a Gromacs particle type. Use 'ok' to see if this
    worked correctly */
GromacsAtomType::PARTICLE_TYPE GromacsAtomType::toParticleType(const QString &word, bool *ok)
{
    QString lword = word.toLower();

    if (lword == "a" or word == "atom")
    {
        if (ok) *ok = true;
        return ATOM;
    }
    else if (lword == "s" or word == "shell")
    {
        if (ok) *ok = true;
        return SHELL;
    }
    else if (lword == "v" or lword == "virtual" or lword == "d" or lword == "dummy")
    {
        if (ok) *ok = true;
        return VIRTUAL;
    }
    else
    {
        if (ok) *ok = false;
        return UNKNOWN_TYPE;
    }
}

/////////
///////// Implementation of GromacsBond
/////////

static const RegisterMetaType<GromacsBond> r_bond(NO_ROOT);

QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, const GromacsBond &bond)
{
    writeHeader(ds, r_bond, 1);
    
    ds << bond.func_type;
    
    for (int i=0; i<bond.count(); ++i)
    {
        ds << bond.k[i];
    }
    
    return ds;
}

QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, GromacsBond &bond)
{
    VersionID v = readHeader(ds, r_bond);
    
    if (v == 1)
    {
        ds >> bond.func_type;
        
        for (int i=0; i<bond.count(); ++i)
        {
            ds >> bond.k[i];
        }
    }
    else
        throw version_error(v, "1", r_bond, CODELOC);
    
    return ds;
}

//the maximum number of parameters for any Gromacs bond function
static const int MAX_BOND_PARAMS = 4;

/** Null constructor */
GromacsBond::GromacsBond()
            : func_type(0)
{
    for (int i=0; i<MAX_BOND_PARAMS; ++i)
    {
        k[i] = 0;
    }
}

/** Construct from the passed 'bond', using 'R' as the symbol for the R value */
GromacsBond::GromacsBond(const SireCAS::Expression &bond, const SireCAS::Symbol &R)
            : func_type(0)
{
    for (int i=0; i<MAX_BOND_PARAMS; ++i)
    {
        k[i] = 0;
    }
    
    // a LOT of introspection will be needed to extract the function type
    // and parameters from a generic expression...
    throw SireError::incomplete_code( QObject::tr("THIS CODE NEEDS WRITING!"), CODELOC );
}

static void assert_valid_bond_function(int func_type)
{
    if (func_type < 1 or func_type > 10)
        throw SireError::invalid_arg( QObject::tr(
            "There is no Gromacs bond function with ID '%1'. The only valid IDs are "
            "the numbers 1-10.").arg(func_type), CODELOC );
}

/** Construct a bond of the specified function type with specified parameters
    (the order should be the same as in the Gromacs Manual, table 5.5) */
GromacsBond::GromacsBond(int function_type,
                         double k0, double k1, double k2, double k3)
            : func_type(function_type)
{
    assert_valid_bond_function(func_type);
    
    k[0] = k0;
    k[1] = k1;
    k[2] = k2;
    k[3] = k3;
    
    for (int i=count(); i<MAX_BOND_PARAMS; ++i)
    {
        k[i] = 0;
    }
}

/** Construct a bond of the specified function type by interpreting the parameter
    data from the passed list of parameter values. These should be in the 
    same order as in the Gromacs Manual, table 5.5 */
GromacsBond::GromacsBond(int function_type, const QList<double> &params)
            : func_type(function_type)
{
    assert_valid_bond_function(func_type);

    if (count() != params.count())
    {
        throw SireError::invalid_arg( QObject::tr(
            "Incorrect number of parameters (%1) passed for a Gromacs bond of type %2. "
            "You need to supply %3 parameters.").arg(params.count()).arg(function_type)
                    .arg(count()), CODELOC );
    }

    for (int i=0; i<count(); ++i)
    {
        k[i] = params[i];
    }

    for (int i=count(); i<MAX_BOND_PARAMS; ++i)
    {
        k[i] = 0;
    }
}

/** Copy constructor */
GromacsBond::GromacsBond(const GromacsBond &other)
            : func_type(other.func_type)
{
    for (int i=0; i<MAX_BOND_PARAMS; ++i)
    {
        k[i] = other.k[i];
    }
}

/** Destructor */
GromacsBond::~GromacsBond()
{}

/** Copy assigment operator */
GromacsBond& GromacsBond::operator=(const GromacsBond &other)
{
    if (this != &other)
    {
        func_type = other.func_type;
        for (int i=0; i<MAX_BOND_PARAMS; ++i)
        {
            k[i] = other.k[i];
        }
    }
    
    return *this;
}

/** Comparison operator */
bool GromacsBond::operator==(const GromacsBond &other) const
{
    if (func_type == other.func_type)
    {
        for (int i=0; i<MAX_BOND_PARAMS; ++i)
        {
            if (k[i] != other.k[i])
                return false;
        }

        return true;
    }
    else
        return false;
}

/** Comparison operator */
bool GromacsBond::operator!=(const GromacsBond &other) const
{
    return not operator==(other);
}

/** Comparison operator */
bool GromacsBond::operator<(const GromacsBond &other) const
{
    if (func_type < other.func_type)
        return true;
    else if (func_type == other.func_type)
    {
        if (k[0] < other.k[0])
            return true;
        else if (k[0] == other.k[0])
        {
            if (k[1] < other.k[1])
                return true;
            else if (k[1] == other.k[1])
            {
                if (k[2] < other.k[2])
                    return true;
                else if (k[2] == other.k[2])
                {
                    if (k[3] < other.k[3])
                        return true;
                }
            }
        }
    }

    return false;
}

/** Comparison operator */
bool GromacsBond::operator<=(const GromacsBond &other) const
{
    return *this == other or *this < other;
}

/** Comparison operator */
bool GromacsBond::operator>(const GromacsBond &other) const
{
    return not (*this <= other);
}

/** Comparison operator */
bool GromacsBond::operator>=(const GromacsBond &other) const
{
    return not (*this < other);
}

const char* GromacsBond::typeName()
{
    return QMetaType::typeName( qMetaTypeId<GromacsBond>() );
}

const char* GromacsBond::what() const
{
    return GromacsBond::typeName();
}

/** Return the ith parameter for this bond */
double GromacsBond::operator[](int i) const
{
    i = Index(i).map(count());
    return k[i];
}

/** Return the ith parameter for this bond */
double GromacsBond::at(int i) const
{
    return operator[](i);
}

/** Return the number of parameters associated with this bond type */
int GromacsBond::count() const
{
    switch(func_type)
    {
        case 1:
        case 2:
        case 6:
        case 7:
        case 8:
        case 9:
            return 2;
        case 3:
        case 4:
            return 3;
        case 5:
            return 0;
        case 10:
            return 4;
        default:
            return 0;
    }
}

/** Return the number of parameters associated with this bond type */
int GromacsBond::size() const
{
    return count();
}

/** Return the Gromacs ID number for the function type for this bond. See table
    5.5 in the Gromacs manual for information */
int GromacsBond::functionType() const
{
    return func_type;
}

/** Return the string description of the function type for this bond */
QString GromacsBond::functionTypeString() const
{
    switch (func_type)
    {
        case 1:
            return "bond";
        case 2:
            return "G96 bond";
        case 3:
            return "Morse bond";
        case 4:
            return "cubic bond";
        case 5:
            return "connection";
        case 6:
            return "harmonic potential";
        case 7:
            return "FENE bond";
        case 8:
            return "tabulated bond";
        case 9:
            return "tabulated bond (no exclusions)";
        case 10:
            return "restraint potential";
        default:
            return "unknown";
    }
}

/** Return all of the parameters for this bond */
QList<double> GromacsBond::parameters() const
{
    QList<double> params;
    
    for (int i=0; i<count(); ++i)
    {
        params.append( k[i] );
    }
    
    return params;
}

/** Return a string representation of this bond */
QString GromacsBond::toString() const
{
    if (func_type == 0)
        return QObject::tr("GromacsBond::null");
    else
    {
        QStringList params;
        
        for (int i=0; i<count(); ++i)
        {
            params.append( QString::number(k[i]) );
        }
        
        if (params.isEmpty())
            return QObject::tr("GromacsBond( functionType() = %1 )")
                    .arg(functionTypeString());
        else
            return QObject::tr("Gromacsbond( functionType() = %1, parameters() = [ %2 ] )")
                    .arg(params.join(", "));
    }
}

/** Return this function converted to a SireCAS::Expression using the passed symbol
    to represent the bond length */
SireCAS::Expression GromacsBond::toExpression(const SireCAS::Symbol &R) const
{
    const double kj_per_mol_per_nm2 = ((kJ_per_mol) / (nanometer*nanometer)).value();
    const double kj_per_mol_per_nm3 = ((kJ_per_mol) / (nanometer*nanometer*nanometer)).value();
    const double kj_per_mol_per_nm4 = ((kJ_per_mol) / (nanometer*nanometer*nanometer*nanometer))
                                                                .value();
    const double nm = nanometer.value();
    const double per_nm = (1.0 / nanometer).value();

    if (func_type == 1 or func_type == 6)
    {
        //standard bond : 0.5 k (r - r0)^2
        const double k0 = k[0] * kj_per_mol_per_nm2;
        const double r0 = k[1] * nm;
        
        return 0.5 * k0 * SireMaths::pow_2(R - r0);
    }
    else if (func_type == 2)
    {
        //gromos 96 bond : 0.25 k (r^2 - r0^2)^2
        const double k0 = k[0] * kj_per_mol_per_nm4;
        const double r0 = k[1] * nm;
        
        return 0.25 * k0 * SireMaths::pow_2(R*R - r0*r0);
    }
    else if (func_type == 3)
    {
        //morse potential : D[ 1 - exp{ -beta(r - r0) }]^2
        const double D = (k[1] * kJ_per_mol).value();
        const double beta = k[2] * per_nm;
        const double r0 = k[0] * nm;
        
        return D * SireMaths::pow_2( 1.0 - Exp( -beta*(R - r0) ) );
    }
    else if (func_type == 4)
    {
        //cubic bond : k1(r - r0)^2 + k1k2(r - r0)^3
        const double k1 = k[1] * kj_per_mol_per_nm2;
        const double k2 = k[2] * kj_per_mol_per_nm3;
        const double r0 = k[0] * nm;
        
        return k1*SireMaths::pow_2(R - r0) + k1*k2*SireMaths::pow_3(R - r0);
    }
    else if (func_type == 5)
    {
        //connection - zero interaction
        return SireCAS::Expression(0);
    }
    else if (func_type == 7)
    {
        //FENE bond : -0.5 k b log( 1 - (r^2/b^2) )
        const double b = k[0] * nm;
        const double k0 = k[1] * kj_per_mol_per_nm2;
        
        return -0.5 * k0 * Ln( 1.0 - ( (R*R)/(b*b)) );
    }
    else if (func_type == 8 or func_type == 9)
    {
        throw SireError::unsupported( QObject::tr(
            "It is not possible to convert a tabulated gromacs bond into a SireCAS::Expression!"),
                CODELOC );
    }
    else if (func_type == 10)
    {
        //restraint potential
        // if r < r0 : 0.5 k(r - r0)^2
        // if r0 <= r < r1 : 0
        // if r1 <= r < r2 : 0.5 k(r - r1)^2
        // else : 0.5 k (r2 - r1)(2r - r2 - r1)
        
        const double r0 = k[0] * nm;
        const double r1 = k[1] * nm;
        const double r2 = k[2] * nm;
        const double k0 = k[3] * kj_per_mol_per_nm2;
        
        return Conditional( LessThan(R, r0), 0.5*k0*SireMaths::pow_2(R - r0),
                  Conditional( LessThan(R, r1), Expression(0),
                      Conditional( LessThan(R, r2), 0.5*k0*SireMaths::pow_2(R - r1),
                          0.5*k0*(r2-r1)*(2.0*R - r2 - r1) )
                             )
                           );
    }
    else
        return SireCAS::Expression(0);
}

inline uint my_qHash(double key)
{
    return ::qHash( *(reinterpret_cast<const ulong*>(&key)) );
}

/** Return a hash for this bond */
uint GromacsBond::hash() const
{
    return  my_qHash(k[0]) | my_qHash(k[1]) | my_qHash(k[2]) | my_qHash(k[3]);
}
