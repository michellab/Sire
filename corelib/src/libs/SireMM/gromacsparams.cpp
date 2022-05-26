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
#include "amberparams.h"

#include "SireID/index.h"

#include "SireError/errors.h"

#include "SireUnits/units.h"

#include "SireCAS/exp.h"
#include "SireCAS/conditional.h"
#include "SireCAS/sum.h"
#include "SireCAS/trigfuncs.h"

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

QDataStream &operator<<(QDataStream &ds, const GromacsAtomType &typ)
{
    writeHeader(ds, r_atomtype, 2);

    SharedDataStream sds(ds);

    sds << typ._typ << typ._btyp << typ._mass.to(g_per_mol)
        << typ._chg.to(mod_electron) << typ.particleTypeString()
        << typ._lj << typ._elem;

    return ds;
}

QDataStream &operator>>(QDataStream &ds, GromacsAtomType &typ)
{
    VersionID v = readHeader(ds, r_atomtype);

    if (v == 2)
    {
        SharedDataStream sds(ds);

        double chg, mass;
        QString ptyp;

        sds >> typ._typ >> typ._btyp >> mass >> chg >> ptyp >> typ._lj >> typ._elem;

        typ._btyp = typ._typ;
        typ._mass = mass * g_per_mol;
        typ._chg = chg * mod_electron;
        typ._ptyp = GromacsAtomType::toParticleType(ptyp);
    }
    else if (v == 1)
    {
        SharedDataStream sds(ds);

        double chg, mass;
        QString ptyp;

        sds >> typ._typ >> mass >> chg >> ptyp >> typ._lj >> typ._elem;

        typ._btyp = typ._typ;
        typ._mass = mass * g_per_mol;
        typ._chg = chg * mod_electron;
        typ._ptyp = GromacsAtomType::toParticleType(ptyp);
    }
    else
        throw version_error(v, "1,2", r_atomtype, CODELOC);

    return ds;
}

/** Assert that this object is sane */
void GromacsAtomType::assertSane() const
{}

/** Null constructor */
GromacsAtomType::GromacsAtomType()
                : _typ(), _btyp(), _mass(0), _chg(0),
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
                : _typ(atom_type), _btyp(atom_type), _mass(mass), _chg(charge),
                  _lj(ljparam), _ptyp(particle_type), _elem(element)
{
    assertSane();
}

/** Construct passing in all parameters */
GromacsAtomType::GromacsAtomType(QString atom_type, QString bond_type,
                                 SireUnits::Dimension::MolarMass mass,
                                 SireUnits::Dimension::Charge charge,
                                 PARTICLE_TYPE particle_type,
                                 const LJParameter &ljparam,
                                 const Element &element)
                : _typ(atom_type), _btyp(bond_type), _mass(mass), _chg(charge),
                  _lj(ljparam), _ptyp(particle_type), _elem(element)
{
    assertSane();
}

/** Construct, specifying only the mass */
GromacsAtomType::GromacsAtomType(QString atom_type, SireUnits::Dimension::MolarMass mass)
                : _typ(atom_type), _btyp(atom_type), _mass(mass), _chg(0),
                  _ptyp(GromacsAtomType::UNKNOWN_TYPE),
                  _elem(0)
{}

/** Copy constructor */
GromacsAtomType::GromacsAtomType(const GromacsAtomType &other)
                : _typ(other._typ), _btyp(other._btyp), _mass(other._mass), _chg(other._chg),
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
        _btyp = other._btyp;
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
    return _typ == other._typ and _btyp == other._btyp and _mass == other._mass and
           _chg == other._chg and
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
    else if (_btyp == _typ)
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
    else
    {
        return QObject::tr("GromacsAtomType( atomType() = %1, bondType() = %7, "
                              "mass() = %2 g mol-1, "
                              "charge() = %3 |e|, "
                              "particleType() = %4, ljParameter() = %5, element() = %6 )")
                              .arg(atomType()).arg(mass().to(g_per_mol))
                              .arg(charge().to(mod_electron))
                              .arg(particleTypeString())
                              .arg(ljParameter().toString())
                              .arg(element().toString())
                              .arg(bondType());
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

/** Set the atom type to the passed value. This is useful if there are duplicate
    atom types with different parameters, as is the case for molecules parameterised
    using OpenFF, i.e. an atom of type "h1" will, in general, be different to an
    atom of type "h1" in a different molecule.
 */
void GromacsAtomType::setAtomType(const QString& atom_type)
{
    this->_typ = atom_type;
}

/** Set the element to the passed value. This is useful if the "atomtype" section has
    invalid mass informtion, as is the case for many topology files generted by acpype.
    This allows us to update the element of the type using the mass from the "atoms"
    section.
 */
void GromacsAtomType::setElement(SireMol::Element elem)
{
    this->_elem = elem;
}

/////////
///////// Implementation of GromacsBond
/////////

static const RegisterMetaType<GromacsBond> r_bond(NO_ROOT);

QDataStream &operator<<(QDataStream &ds, const GromacsBond &bond)
{
    writeHeader(ds, r_bond, 1);

    ds << bond.func_type;

    for (int i=0; i<bond.count(); ++i)
    {
        ds << bond.k[i];
    }

    return ds;
}

QDataStream &operator>>(QDataStream &ds, GromacsBond &bond)
{
    VersionID v = readHeader(ds, r_bond);

    if (v == 1)
    {
        bond = GromacsBond();

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

    //first, let's see if this is an Amber-style harmonic bond
    {
        AmberBond amberbond;
        bool is_amber = false;

        try
        {
            amberbond = AmberBond(bond, R);
            is_amber = true;
        }
        catch(...)
        {}

        if (is_amber)
        {
            //yes, this is a valid amber bond ( 0.5 kb (r - r0)^2 )
            double kb = 2.0 * amberbond.k();
            double r0 = amberbond.r0();

            if (kb == 0)
            {
                //this is a null bond (connection)
                func_type = 5;
            }
            else
            {
                //this is a harmonic bond
                func_type = 1;

                const double kj_per_mol_per_nm2 = ((kJ_per_mol) / (nanometer*nanometer)).value();
                const double nm = nanometer.value();

                k[0] = r0 / nm;
                k[1] = kb / kj_per_mol_per_nm2;
            }

            return;
        }
    }

    // a LOT of introspection will be needed to extract the function type
    // and parameters from a generic expression...
    throw SireError::incomplete_code( QObject::tr("Sire cannot yet interpret bonds "
       "that are not in a standard harmonic format! (%1)").arg(bond.toString()), CODELOC );
}

static void assert_valid_bond_function(int func_type)
{
    if (func_type < 1 or func_type > 10)
        throw SireError::invalid_arg( QObject::tr(
            "There is no Gromacs bond function with ID '%1'. The only valid IDs are "
            "the numbers 1-10.").arg(func_type), CODELOC );
}

//the value used to indicate that the parameter needs to be resolved
static double unresolved_parameter_value = std::numeric_limits<double>::infinity();

/** Construct a bond that is of the specified type, but the parameters have yet
    to be resolved. This is because Gromacs can indicate the required type of
    function in the molecule specification, without providing the parameters */
GromacsBond::GromacsBond(int function_type)
            : func_type(function_type)
{
    assert_valid_bond_function(func_type);

    for (int i=0; i<MAX_BOND_PARAMS; ++i)
    {
        k[i] = unresolved_parameter_value;
    }
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

    if (params.count() < count())
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

/** Return whether or not this parameter needs resolving */
bool GromacsBond::needsResolving() const
{
    //don't need to resolve the 'connection' function type (5), as it has no parameters
    return func_type != 5 and (k[0] == unresolved_parameter_value or func_type == 0);
}

/** Return whether or not the parameters for this bond are resolved */
bool GromacsBond::isResolved() const
{
    return not needsResolving();
}

/** Assert that the parameters for this bond have been resolved */
void GromacsBond::assertResolved() const
{
    if (needsResolving())
        throw SireError::invalid_state( QObject::tr(
            "The parameters for this GromacsBond have not been resolved! %1")
                .arg(this->toString()), CODELOC );
}

/** Return the ith parameter for this bond */
double GromacsBond::operator[](int i) const
{
    assertResolved();
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

/** All Gromacs bonds are simple (just a function of the bond length) */
bool GromacsBond::isSimple() const
{
    return true;
}

/** Return whether or not this is a harmonic bond */
bool GromacsBond::isHarmonic() const
{
    switch (func_type)
    {
    case 0: //null, so zero, which is a zero harmonic
    case 1: //harmonic
    case 5: //connection, so zero, which is a zero harmonic
    case 6: //harmonic
        return true;
    default:
        qDebug() << "NOT HARMONIC!" << func_type << this->toString();
        return false;
    }
}

/** Return all of the parameters for this bond */
QList<double> GromacsBond::parameters() const
{
    QList<double> params;

    if (isResolved())
    {
        for (int i=0; i<count(); ++i)
        {
            params.append( k[i] );
        }
    }

    return params;
}

/** Return a string representation of this bond */
QString GromacsBond::toString() const
{
    if (func_type == 0)
        return QObject::tr("GromacsBond::null");
    else if (needsResolving())
    {
        return QObject::tr("GromacsBond( functionType() = %1, needsResolving )")
                    .arg(functionTypeString());
    }
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
                    .arg(functionTypeString()).arg(params.join(", "));
    }
}

/** Return whether or not this GromacsBond implies that the atoms are actually
    bonded together */
bool GromacsBond::atomsAreBonded() const
{
    return true;
}

/** Return the equilibrium length of this bond */
SireUnits::Dimension::Length GromacsBond::equilibriumLength() const
{
    assertResolved();

    const double nm = nanometer.value();

    if (func_type == 1 or func_type == 6)
    {
        //standard bond : 0.5 k (r - r0)^2
        const double r0 = k[0] * nm;
        return SireUnits::Dimension::Length(r0);
    }
    else if (func_type == 2)
    {
        //gromos 96 bond : 0.25 k (r^2 - r0^2)^2
        const double r0 = k[0] * nm;
        return SireUnits::Dimension::Length(r0);
    }
    else if (func_type == 3)
    {
        //morse potential : D[ 1 - exp{ -beta(r - r0) }]^2
        const double r0 = k[0] * nm;
        return SireUnits::Dimension::Length(r0);
    }
    else if (func_type == 4)
    {
        //cubic bond : k1(r - r0)^2 + k1k2(r - r0)^3
        const double r0 = k[0] * nm;
        return SireUnits::Dimension::Length(r0);
    }
    else if (func_type == 5)
    {
        //connection - zero interaction
        return SireUnits::Dimension::Length(0);
    }
    else if (func_type == 7)
    {
        //FENE bond : -0.5 k b log( 1 - (r^2/b^2) )
        const double b = k[0] * nm;
        return SireUnits::Dimension::Length(b);
    }
    else if (func_type == 8 or func_type == 9)
    {
        throw SireError::unsupported( QObject::tr(
            "It is not possible to get an equilibrium length from a tabulated gromacs bond"),
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
        return SireUnits::Dimension::Length(r0);
    }
    else
        return SireUnits::Dimension::Length(0);
}

/** Return this function converted to a SireCAS::Expression using the passed symbol
    to represent the bond length */
SireCAS::Expression GromacsBond::toExpression(const SireCAS::Symbol &R) const
{
    assertResolved();

    const double kj_per_mol_per_nm2 = ((kJ_per_mol) / (nanometer*nanometer)).value();
    const double kj_per_mol_per_nm3 = ((kJ_per_mol) / (nanometer*nanometer*nanometer)).value();
    const double kj_per_mol_per_nm4 = ((kJ_per_mol) / (nanometer*nanometer*nanometer*nanometer))
                                                                .value();
    const double nm = nanometer.value();
    const double per_nm = (1.0 / nanometer).value();

    if (func_type == 1 or func_type == 6)
    {
        //standard bond : 0.5 k (r - r0)^2
        const double k0 = k[1] * kj_per_mol_per_nm2;
        const double r0 = k[0] * nm;

        return 0.5 * k0 * SireMaths::pow_2(R - r0);
    }
    else if (func_type == 2)
    {
        //gromos 96 bond : 0.25 k (r^2 - r0^2)^2
        const double k0 = k[1] * kj_per_mol_per_nm4;
        const double r0 = k[0] * nm;

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

SIRE_ALWAYS_INLINE uint my_qHash(double key)
{
    return ::qHash( *(reinterpret_cast<const ulong*>(&key)) );
}

/** Return a hash for this bond */
uint GromacsBond::hash() const
{
    return  my_qHash(k[0]) | my_qHash(k[1]) | my_qHash(k[2]) | my_qHash(k[3]);
}


/////////
///////// Implementation of GromacsAngle
/////////

static const RegisterMetaType<GromacsAngle> r_ang(NO_ROOT);

QDataStream &operator<<(QDataStream &ds, const GromacsAngle &ang)
{
    writeHeader(ds, r_ang, 1);

    ds << ang.func_type;

    for (int i=0; i<ang.count(); ++i)
    {
        ds << ang.k[i];
    }

    return ds;
}

QDataStream &operator>>(QDataStream &ds, GromacsAngle &ang)
{
    VersionID v = readHeader(ds, r_ang);

    if (v == 1)
    {
        ang = GromacsAngle();

        ds >> ang.func_type;

        for (int i=0; i<ang.count(); ++i)
        {
            ds >> ang.k[i];
        }
    }
    else
        throw version_error(v, "1", r_ang, CODELOC);

    return ds;
}

//the maximum number of parameters for any Gromacs angle function
static const int MAX_ANGLE_PARAMS = 6;

/** Null constructor */
GromacsAngle::GromacsAngle()
             : func_type(0)
{
    for (int i=0; i<MAX_ANGLE_PARAMS; ++i)
    {
        k[i] = 0;
    }
}

/** Construct from the passed 'angle', using 'theta' as the symbol for the theta value */
GromacsAngle::GromacsAngle(const SireCAS::Expression &angle, const SireCAS::Symbol &theta)
             : func_type(0)
{
    for (int i=0; i<MAX_ANGLE_PARAMS; ++i)
    {
        k[i] = 0;
    }

    //first, let's see if this is an Amber-style harmonic angle
    {
        AmberAngle amberangle;
        bool is_amber = false;

        try
        {
            amberangle = AmberAngle(angle, theta);
            is_amber = true;
        }
        catch(...)
        {}

        if (is_amber)
        {
            //yes, this is a valid amber angle ( 0.5 kb (theta - theta0)^2 )
            double kb = 2.0 * amberangle.k();
            double t0 = amberangle.theta0();

            //this is a harmonic angle
            func_type = 1;

            const double kj_per_mol_per_rad2 = ((kJ_per_mol) / (radian*radian)).value();
            const double deg = degree.value();

            k[0] = t0 / deg;
            k[1] = kb / kj_per_mol_per_rad2;

            return;
        }
    }

    // a LOT of introspection will be needed to extract the function type
    // and parameters from a generic expression...
    throw SireError::incomplete_code( QObject::tr("Sire cannot yet interpret angles "
       "that are not in a standard harmonic format! (%1)").arg(angle.toString()), CODELOC );
}

static void assert_valid_angle_function(int func_type)
{
    if (func_type < 1 or func_type > 10 or func_type == 7 or func_type == 9)
        throw SireError::invalid_arg( QObject::tr(
            "There is no Gromacs angle function with ID '%1'. The only valid IDs are "
            "the numbers 1-6,8,10.").arg(func_type), CODELOC );
}

/** Construct an angle that is of the specified type, but the parameters have yet
    to be resolved. This is because Gromacs can indicate the required type of
    function in the molecule specification, without providing the parameters */
GromacsAngle::GromacsAngle(int function_type)
             : func_type(function_type)
{
    assert_valid_angle_function(func_type);

    for (int i=0; i<MAX_ANGLE_PARAMS; ++i)
    {
        k[i] = unresolved_parameter_value;
    }
}

/** Construct an angle of the specified function type with specified parameters
    (the order should be the same as in the Gromacs Manual, table 5.5) */
GromacsAngle::GromacsAngle(int function_type,
                          double k0, double k1, double k2, double k3, double k4, double k5)
             : func_type(function_type)
{
    assert_valid_angle_function(func_type);

    k[0] = k0;
    k[1] = k1;
    k[2] = k2;
    k[3] = k3;
    k[4] = k4;
    k[5] = k5;

    for (int i=count(); i<MAX_ANGLE_PARAMS; ++i)
    {
        k[i] = 0;
    }
}

/** Construct an angle of the specified function type by interpreting the parameter
    data from the passed list of parameter values. These should be in the
    same order as in the Gromacs Manual, table 5.5 */
GromacsAngle::GromacsAngle(int function_type, const QList<double> &params)
             : func_type(function_type)
{
    assert_valid_angle_function(func_type);

    if (params.count() < count())
    {
        throw SireError::invalid_arg( QObject::tr(
            "Incorrect number of parameters (%1) passed for a Gromacs angle of type %2. "
            "You need to supply %3 parameters.").arg(params.count()).arg(function_type)
                    .arg(count()), CODELOC );
    }

    for (int i=0; i<count(); ++i)
    {
        k[i] = params[i];
    }

    for (int i=count(); i<MAX_ANGLE_PARAMS; ++i)
    {
        k[i] = 0;
    }
}

/** Copy constructor */
GromacsAngle::GromacsAngle(const GromacsAngle &other)
             : func_type(other.func_type)
{
    for (int i=0; i<MAX_ANGLE_PARAMS; ++i)
    {
        k[i] = other.k[i];
    }
}

/** Destructor */
GromacsAngle::~GromacsAngle()
{}

/** Copy assigment operator */
GromacsAngle& GromacsAngle::operator=(const GromacsAngle &other)
{
    if (this != &other)
    {
        func_type = other.func_type;
        for (int i=0; i<MAX_ANGLE_PARAMS; ++i)
        {
            k[i] = other.k[i];
        }
    }

    return *this;
}

/** Comparison operator */
bool GromacsAngle::operator==(const GromacsAngle &other) const
{
    if (func_type == other.func_type)
    {
        for (int i=0; i<MAX_ANGLE_PARAMS; ++i)
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
bool GromacsAngle::operator!=(const GromacsAngle &other) const
{
    return not operator==(other);
}

/** Comparison operator */
bool GromacsAngle::operator<(const GromacsAngle &other) const
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
                    else if (k[3] == other.k[3])
                    {
                        if (k[4] < other.k[4])
                            return true;
                        else if (k[4] == other.k[4])
                        {
                            if (k[5] < other.k[5])
                                return true;
                        }
                    }
                }
            }
        }
    }

    return false;
}

/** Comparison operator */
bool GromacsAngle::operator<=(const GromacsAngle &other) const
{
    return *this == other or *this < other;
}

/** Comparison operator */
bool GromacsAngle::operator>(const GromacsAngle &other) const
{
    return not (*this <= other);
}

/** Comparison operator */
bool GromacsAngle::operator>=(const GromacsAngle &other) const
{
    return not (*this < other);
}

const char* GromacsAngle::typeName()
{
    return QMetaType::typeName( qMetaTypeId<GromacsAngle>() );
}

const char* GromacsAngle::what() const
{
    return GromacsAngle::typeName();
}

/** Return the ith parameter for this angle */
double GromacsAngle::operator[](int i) const
{
    i = Index(i).map(count());
    return k[i];
}

/** Return the ith parameter for this angle */
double GromacsAngle::at(int i) const
{
    return operator[](i);
}

/** Return the number of parameters associated with this angle type */
int GromacsAngle::count() const
{
    switch(func_type)
    {
        case 1:
        case 2:
        case 8:
        case 10:
            return 2;
        case 3:
            return 3;
        case 4:
        case 5:
            return 4;
        case 6:
            return 5;
        default:
            return 0;
    }
}

/** Return the number of parameters associated with this angle type */
int GromacsAngle::size() const
{
    return count();
}

/** Return the Gromacs ID number for the function type for this angle. See table
    5.5 in the Gromacs manual for information */
int GromacsAngle::functionType() const
{
    return func_type;
}

/** Return whether or not this is a simple angle function, based only on the
    size of the angle */
bool GromacsAngle::isSimple() const
{
    switch (func_type)
    {
        case 1:
        case 2:
        case 6:
        case 8:
        case 10:
            return true;
        case 3:
        case 4:
            return false;
        case 5:
            return k[3] == 0;  // zero UB bond k value, so no bond term
        default:
            return false;
    }
}

/** Return whether or not this angle is really a mix of multiple bond terms */
bool GromacsAngle::isBondBondCrossTerm() const
{
    return (func_type == 3 or func_type == 4);
}

/** Return whether or not this angle is really a mix of bond and angle terms */
bool GromacsAngle::isBondAngleCrossTerm() const
{
    return (func_type == 5) and (k[3] != 0);
}

/** Return only the bond term part of this angle */
GromacsBond GromacsAngle::toBondTerm() const
{
    if (isBondAngleCrossTerm())
    {
        //Urey-Bradley bond term uses k[2] as r0 and k[3] as kb
        return GromacsBond(1, k[2], k[3]);
    }
    else
        return GromacsBond();
}

/** Return only the angle term part of this angle */
GromacsAngle GromacsAngle::toAngleTerm() const
{
    if (isBondBondCrossTerm())
    {
        return GromacsAngle();
    }
    else if (func_type == 5)
    {
        //Urey-Bradley angle term uses k[0] as theta0 and k[1] as k_t
        return GromacsAngle(1, k[0], k[1]);
    }
    else
        //this is already a simple angle
        return *this;
}

/** Return the string description of the function type for this angle */
QString GromacsAngle::functionTypeString() const
{
    switch (func_type)
    {
        case 1:
            return "angle";
        case 2:
            return "G96 angle";
        case 3:
            return "cross bond-bond";
        case 4:
            return "cross bond-angle";
        case 5:
            return "Urey-Bradley";
        case 6:
            return "quartic angle";
        case 8:
            return "tabulated angle";
        case 10:
            return "restricted bending potential";
        default:
            return "unknown";
    }
}

/** Return whether or not this is a harmonic angle */
bool GromacsAngle::isHarmonic() const
{
    switch (func_type)
    {
    case 0: //null, so zero, which is a zero harmonic
    case 1: //harmonic
        return true;
    default:
        return false;
    }
}

/** Return all of the parameters for this angle */
QList<double> GromacsAngle::parameters() const
{
    QList<double> params;

    for (int i=0; i<count(); ++i)
    {
        params.append( k[i] );
    }

    return params;
}

/** Return a string representation of this angle */
QString GromacsAngle::toString() const
{
    if (func_type == 0)
        return QObject::tr("GromacsAngle::null");
    else if (this->needsResolving())
    {
        return QObject::tr("GromacsAngle( functionType() = %1, needsResolving )")
                .arg(functionTypeString());
    }
    else
    {
        QStringList params;

        for (int i=0; i<count(); ++i)
        {
            params.append( QString::number(k[i]) );
        }

        if (params.isEmpty())
            return QObject::tr("GromacsAngle( functionType() = %1 )")
                    .arg(functionTypeString());
        else
            return QObject::tr("GromacsAngle( functionType() = %1, parameters() = [ %2 ] )")
                    .arg(functionTypeString()).arg(params.join(", "));
    }
}

/** Return this function converted to a SireCAS::Expression using the passed symbol
    to represent the angle size */
SireCAS::Expression GromacsAngle::toExpression(const SireCAS::Symbol &theta) const
{
    const double kj_per_mol = kJ_per_mol.value();
    const double kj_per_mol_per_rad = ((kJ_per_mol) / (radian)).value();
    const double kj_per_mol_per_rad2 = ((kJ_per_mol) / (radian*radian)).value();
    const double kj_per_mol_per_rad3 = ((kJ_per_mol) / (radian*radian*radian)).value();
    const double kj_per_mol_per_rad4 = ((kJ_per_mol) / (radian*radian*radian*radian)).value();
    const double kj_per_mol_per_nm2 = ((kJ_per_mol) / (nanometer*nanometer)).value();

    const double deg = degree.value();

    if (func_type == 1)
    {
        //standard angle : 0.5 k (theta - theta0)^2
        const double k0 = k[1] * kj_per_mol_per_rad2;
        const double theta0 = k[0] * deg;

        return 0.5 * k0 * SireMaths::pow_2(theta - theta0);
    }
    else if (func_type == 2)
    {
        //gromos 96 angle : 0.5 k (cos(theta) - cos(theta0))^2
        const double k0 = k[1] * kj_per_mol;
        const double theta0 = k[0] * deg;

        return 0.5 * k0 * SireMaths::pow_2( Cos(theta) - std::cos(theta0) );
    }
    else if (func_type == 3 or func_type == 4)
    {
        throw SireError::incompatible_error( QObject::tr(
            "Cannot convert the bond-bond type Gromacs angle '%1' to an expression "
            "using only theta.").arg(this->toString()), CODELOC );
    }
    else if (func_type == 5)
    {
        //Urey-Bradley : 0.5 k_t (theta - theta0)^2 + 0.5 k_b ( r - r0 )^2
        const double k0 = k[1] * kj_per_mol_per_rad2;
        const double theta0 = k[0] * deg;
        const double kb = k[3] * kj_per_mol_per_nm2;
        //const double r0 = k[2] * nm;

        if (kb != 0)
            throw SireError::incompatible_error( QObject::tr(
                "Cannot convert a Urey-Bradley Gromacs angle into an expression of only "
                "the angle size if the UB bond force constant is non-zero - %1")
                    .arg(this->toString()), CODELOC );

        return 0.5 * k0 * SireMaths::pow_2(theta - theta0);
    }
    else if (func_type == 6)
    {
        //Quartic angle : Sum_(n=0,5) C_n( theta - theta0 )^n
        const double k0 = k[1] * kj_per_mol;
        const double k1 = k[2] * kj_per_mol_per_rad;
        const double k2 = k[3] * kj_per_mol_per_rad2;
        const double k3 = k[4] * kj_per_mol_per_rad3;
        const double k4 = k[5] * kj_per_mol_per_rad4;
        const double theta0 = k[0] * deg;

        return (k0) + (k1 * (theta - theta0))
                    + (k2 * SireMaths::pow_2(theta - theta0))
                    + (k3 * SireMaths::pow_3(theta - theta0))
                    + (k4 * SireMaths::pow_4(theta - theta0));
    }
    else if (func_type == 8)
    {
        throw SireError::unsupported( QObject::tr(
            "It is not possible to convert a tabulated gromacs angle into a SireCAS::Expression!"),
                CODELOC );
    }
    else if (func_type == 10)
    {
        //restricted bending potential : 0.5 k [ (cos(theta) - cos(theta0))^2 / sin(theta)^2 ]
        const double k0 = k[1] * kj_per_mol;
        const double theta0 = k[0] * deg;

        return 0.5 * k0 * ( SireMaths::pow_2(Cos(theta) - std::cos(theta0)) /
                            SireMaths::pow_2(Sin(theta)) );
    }
    else
        return SireCAS::Expression(0);
}

/** Return this function converted to a SireCAS::Expression using the passed symbol
    to represent the bond length (r02) and angle size (t012) */
SireCAS::Expression GromacsAngle::toBondAngleExpression(const SireCAS::Symbol &r,
                                                        const SireCAS::Symbol &theta) const
{
    const double kj_per_mol_per_rad2 = ((kJ_per_mol) / (radian*radian)).value();
    const double deg = degree.value();

    const double kj_per_mol_per_nm2 = ((kJ_per_mol) / (nanometer*nanometer)).value();
    const double nm = nanometer.value();

    if (func_type == 5)
    {
        //Urey-Bradley : 0.5 k_t (theta - theta0)^2 + 0.5 k_b ( r - r0 )^2
        const double k0 = k[1] * kj_per_mol_per_rad2;
        const double theta0 = k[0] * deg;
        const double kb = k[3] * kj_per_mol_per_nm2;
        const double r0 = k[2] * nm;

        return 0.5 * k0 * SireMaths::pow_2(theta - theta0) +
               0.5 * kb * SireMaths::pow_2(r - r0);
    }
    else
    {
        //angle only
        return this->toExpression(theta);
    }
}

/** Return this function converted to a SireCAS::Expression using the passed symbol
    to represent the bond lengths r01, r12 and r02 */
SireCAS::Expression GromacsAngle::toBondBondExpression(const SireCAS::Symbol &r01,
                                                       const SireCAS::Symbol &r12,
                                                       const SireCAS::Symbol &r02) const
{
    const double kj_per_mol_per_nm2 = ((kJ_per_mol) / (nanometer*nanometer)).value();
    const double nm = nanometer.value();

    if (func_type == 3)
    {
        // cross bond-bond : kb (r01 - R1) (r12 - R2)
        const double kb = k[2] * kj_per_mol_per_nm2;
        const double R1 = k[0] * nm;
        const double R2 = k[1] * nm;

        return kb * (r01 - R1) * (r12 - R2);
    }
    else if (func_type == 4)
    {
        // cross bond-angle : kb (r02 - R3) ( r01 - R1 + r12 - R2 )
        const double kb = k[3] * kj_per_mol_per_nm2;
        const double R1 = k[0] * nm;
        const double R2 = k[1] * nm;
        const double R3 = k[2] * nm;

        return kb * (r02 - R3) * (r01 - R1 + r12 - R2);
    }
    else
    {
        throw SireError::incompatible_error( QObject::tr(
                "The Gromacs angle of type %1 cannot be expressed as an expression "
                "of bond lengths.")
                    .arg(this->functionTypeString()), CODELOC );
    }
}

/** Return a hash for this bond */
uint GromacsAngle::hash() const
{
    return  my_qHash(k[0]) | my_qHash(k[1]) | my_qHash(k[2]) | my_qHash(k[3]) |
            my_qHash(k[4]) | my_qHash(k[5]);
}

/** Return whether or not this parameter needs resolving */
bool GromacsAngle::needsResolving() const
{
    return k[0] == unresolved_parameter_value or func_type == 0;
}

/** Return whether or not the parameters for this angle are resolved */
bool GromacsAngle::isResolved() const
{
    return not needsResolving();
}

/** Assert that the parameters for this angle have been resolved */
void GromacsAngle::assertResolved() const
{
    if (needsResolving())
        throw SireError::invalid_state( QObject::tr(
            "The parameters for this GromacsAngle have not been resolved! %1")
                .arg(this->toString()), CODELOC );
}

/////////
///////// Implementation of GromacsDihedral
/////////

static const RegisterMetaType<GromacsDihedral> r_dih(NO_ROOT);

QDataStream &operator<<(QDataStream &ds, const GromacsDihedral &dih)
{
    writeHeader(ds, r_dih, 1);

    ds << dih.func_type;

    for (int i=0; i<dih.count(); ++i)
    {
        ds << dih.k[i];
    }

    return ds;
}

QDataStream &operator>>(QDataStream &ds, GromacsDihedral &dih)
{
    VersionID v = readHeader(ds, r_dih);

    if (v == 1)
    {
        dih = GromacsDihedral();

        ds >> dih.func_type;

        for (int i=0; i<dih.count(); ++i)
        {
            ds >> dih.k[i];
        }
    }
    else
        throw version_error(v, "1", r_dih, CODELOC);

    return ds;
}

//the maximum number of parameters for any Gromacs dihedral function
static const int MAX_DIHEDRAL_PARAMS = 6;

/** Null constructor */
GromacsDihedral::GromacsDihedral()
                : func_type(0)
{
    for (int i=0; i<MAX_DIHEDRAL_PARAMS; ++i)
    {
        k[i] = 0;
    }
}

/** Construct from the passed 'dihedral', using 'phi' as the symbol for the phi value */
QList<GromacsDihedral> GromacsDihedral::construct(const Expression &dihedral, const Symbol &phi)
{
    /* We first check whether the dihedral is of the Ryckaert-Bellemans form
       by repeatedly dividing the series terms by cos(phi - pi) to work out
       the constant prefactors. If we only have a single constant term, then
       this isn't the correct dihedral type and we check to see if it is an
       Amber-style cosine dihedral.
     */
    bool is_ryckaert_bellemans = true;
    auto exp = Expression(Cos(phi - SireMaths::pi));

    // Loop over all terms in the series.
    if (dihedral.base().isA<Sum>())
    {
        // Create a vector to hold the six possible factors. (Zero by default.)
        QVector<double> factors = QVector<double>(6);

        // Loop over all children.
        for (auto child : dihedral.base().asA<Sum>().children())
        {
            // Divide the child by cos(phi - pi) until it is a constant.
            // Ryckaert-Bellemans only goes up to cos(phi - pi)^6.
            int n = 0;
            while (not child.isConstant() and n < 6)
            {
                child = child.divide(exp);
                n++;
            }

            // We didn't find a constant factor, this doesn't look like the
            // correct dihedral type.
            if (n == 6)
            {
                is_ryckaert_bellemans = false;
                break;
            }
            else
            {
                // Store the factor.
                const double kj_per_mol = kJ_per_mol.value();
                factors[n] = child.factor() / kj_per_mol;
            }
        }

        if (is_ryckaert_bellemans)
        {
            QList<GromacsDihedral> dihs;
            dihs.append(GromacsDihedral(3, factors[0], factors[1], factors[2],
                                           factors[3], factors[4], factors[5]));

            return dihs;
        }
    }

    // Otherwise, we check whether this is an Amber-style cosine-based dihedral.
    AmberDihedral amberdihedral;
    bool is_amber = false;

    try
    {
        // Cast as an AmberDihedral, but don't try to re-cast as a
        // GromacsDihedral to avoid an infinite loop. The last parameter
        // (set to false) allows us to check whether the dihedral expression
        // is in the Gromacs Ryckaert-Bellemans form.
        amberdihedral = AmberDihedral(dihedral, phi, false);
        is_amber = true;
    }
    catch(...)
    {}

    if (is_amber)
    {
        // Yes, this is a valid amber dihedral. We need to create one
        // GromacsDihedral for each AmberDihPart.
        QList<GromacsDihedral> dihs;

        bool multiterm = amberdihedral.terms().count() > 1;

        for (const auto &amberdih : amberdihedral.terms())
        {
            double kb = amberdih.k();
            double per = amberdih.periodicity();
            double phase = amberdih.phase();

            // This is a cosine dihedral in from k [ 1 + cos(per phi - phase) ]
            // (will one day have to work out how to say this is an improper rather
            // than a dihedral...)
            int func_type = 1;

            // Multiple periodic dihedral.
            if (multiterm)
                func_type = 9;

            const double kj_per_mol = kJ_per_mol.value();
            const double deg = degree.value();

            phase = phase / deg;
            kb = kb / kj_per_mol;

            dihs.append( GromacsDihedral(func_type, phase, kb, per) );
        }

        return dihs;
    }

    // If we get here then this isn't a recognised AmberDihedral form. A LOT of introspection
    // will be needed to extract the function type and parameters from a generic expression...
    throw SireError::incomplete_code( QObject::tr("Sire cannot yet interpret dihedrals "
       "that are not in a standard cosine format! (%1)").arg(dihedral.toString()), CODELOC );
}

/** Construct from the passed 'improper', using 'phi' as the symbol for the phi value */
QList<GromacsDihedral> GromacsDihedral::constructImproper(
                                const Expression &dihedral, const Symbol &phi)
{
    auto parts = GromacsDihedral::construct(dihedral, phi);

    for (auto &part : parts)
    {
        if (part.functionType() == 1 or part.functionType() == 9)
        {
            part.func_type = 4;
        }
    }

    return parts;
}

/** Construct from the passed 'dihedral', using 'phi' as the symbol for the phi value */
GromacsDihedral::GromacsDihedral(const SireCAS::Expression &dihedral, const SireCAS::Symbol &phi)
                : func_type(0)
{
    for (int i=0; i<MAX_DIHEDRAL_PARAMS; ++i)
    {
        k[i] = 0;
    }

    auto parts = GromacsDihedral::construct(dihedral,phi);

    if (parts.count() == 1)
    {
        this->operator=(parts[0]);
        return;
    }
    else if (parts.count() > 1)
        throw SireError::incompatible_error( QObject::tr(
            "The passed expression (%1) is made up of multiple GromacsDihedral terms (%2). "
            "Please use GromacsDihedral::construct(...) to convert the expression into "
            "a valid set of GromacsDihedrals.")
                .arg(dihedral.toString()).arg(Sire::toString(parts)), CODELOC );
}

static void assert_valid_dihedral_function(int func_type)
{
    if (func_type < 1 or func_type > 11 or func_type == 6 or func_type == 7)
        throw SireError::invalid_arg( QObject::tr(
            "There is no Gromacs dihedral function with ID '%1'. The only valid IDs are "
            "the numbers 1-5,8-11.").arg(func_type), CODELOC );
}

/** Construct a dihedral that is of the specified type, but the parameters have yet
    to be resolved. This is because Gromacs can indicate the required type of
    function in the molecule specification, without providing the parameters */
GromacsDihedral::GromacsDihedral(int function_type)
                : func_type(function_type)
{
    assert_valid_dihedral_function(func_type);

    for (int i=0; i<MAX_DIHEDRAL_PARAMS; ++i)
    {
        k[i] = unresolved_parameter_value;
    }
}

/** Construct an dihedral of the specified function type with specified parameters
    (the order should be the same as in the Gromacs Manual, table 5.5) */
GromacsDihedral::GromacsDihedral(int function_type,
                                 double k0, double k1, double k2, double k3, double k4, double k5)
                : func_type(function_type)
{
    assert_valid_dihedral_function(func_type);

    k[0] = k0;
    k[1] = k1;
    k[2] = k2;
    k[3] = k3;
    k[4] = k4;
    k[5] = k5;

    for (int i=count(); i<MAX_DIHEDRAL_PARAMS; ++i)
    {
        k[i] = 0;
    }
}

/** Construct a dihedral of the specified function type by interpreting the parameter
    data from the passed list of parameter values. These should be in the
    same order as in the Gromacs Manual, table 5.5 */
GromacsDihedral::GromacsDihedral(int function_type, const QList<double> &params)
                : func_type(function_type)
{
    assert_valid_dihedral_function(func_type);

    if (params.count() < count())
    {
        throw SireError::invalid_arg( QObject::tr(
            "Incorrect number of parameters (%1) passed for a Gromacs dihedral of type %2. "
            "You need to supply %3 parameters.").arg(params.count()).arg(function_type)
                    .arg(count()), CODELOC );
    }

    for (int i=0; i<count(); ++i)
    {
        k[i] = params[i];
    }

    for (int i=count(); i<MAX_DIHEDRAL_PARAMS; ++i)
    {
        k[i] = 0;
    }
}

/** Copy constructor */
GromacsDihedral::GromacsDihedral(const GromacsDihedral &other)
                : func_type(other.func_type)
{
    for (int i=0; i<MAX_DIHEDRAL_PARAMS; ++i)
    {
        k[i] = other.k[i];
    }
}

/** Destructor */
GromacsDihedral::~GromacsDihedral()
{}

/** Copy assigment operator */
GromacsDihedral& GromacsDihedral::operator=(const GromacsDihedral &other)
{
    if (this != &other)
    {
        func_type = other.func_type;
        for (int i=0; i<MAX_DIHEDRAL_PARAMS; ++i)
        {
            k[i] = other.k[i];
        }
    }

    return *this;
}

/** Comparison operator */
bool GromacsDihedral::operator==(const GromacsDihedral &other) const
{
    if (func_type == other.func_type)
    {
        for (int i=0; i<MAX_DIHEDRAL_PARAMS; ++i)
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
bool GromacsDihedral::operator!=(const GromacsDihedral &other) const
{
    return not operator==(other);
}

/** Comparison operator */
bool GromacsDihedral::operator<(const GromacsDihedral &other) const
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
                    else if (k[3] == other.k[3])
                    {
                        if (k[4] < other.k[4])
                            return true;
                        else if (k[4] == other.k[4])
                        {
                            if (k[5] < other.k[5])
                                return true;
                        }
                    }
                }
            }
        }
    }

    return false;
}

/** Comparison operator */
bool GromacsDihedral::operator<=(const GromacsDihedral &other) const
{
    return *this == other or *this < other;
}

/** Comparison operator */
bool GromacsDihedral::operator>(const GromacsDihedral &other) const
{
    return not (*this <= other);
}

/** Comparison operator */
bool GromacsDihedral::operator>=(const GromacsDihedral &other) const
{
    return not (*this < other);
}

const char* GromacsDihedral::typeName()
{
    return QMetaType::typeName( qMetaTypeId<GromacsDihedral>() );
}

const char* GromacsDihedral::what() const
{
    return GromacsDihedral::typeName();
}

/** Return the ith parameter for this dihedral */
double GromacsDihedral::operator[](int i) const
{
    i = Index(i).map(count());
    return k[i];
}

/** Return the ith parameter for this dihedral */
double GromacsDihedral::at(int i) const
{
    return operator[](i);
}

/** Return the number of parameters associated with this dihedral type */
int GromacsDihedral::count() const
{
    switch(func_type)
    {
        case 1:
        case 4:
        case 9:
            return 3;
        case 2:
        case 8:
        case 10:
            return 2;
        case 3:
            return 6;
        case 5:
            return 4;
        case 11:
            return 5;
        default:
            return 0;
    }
}

/** Return the number of parameters associated with this dihedral type */
int GromacsDihedral::size() const
{
    return count();
}

/** Return the Gromacs ID number for the function type for this dihedral. See table
    5.5 in the Gromacs manual for information */
int GromacsDihedral::functionType() const
{
    return func_type;
}

/** Return whether or not this is a simple dihedral function, based only on the
    size of the torsion */
bool GromacsDihedral::isSimple() const
{
    switch (func_type)
    {
        case 1:
        case 3:
        case 4:
        case 5:
        case 8:
        case 9:
        case 10:
            return true;
        case 2:
        case 11:
            return false;
        default:
            return false;
    }
}

/** Return whether or not this dihedral is really an improper angle term */
bool GromacsDihedral::isImproperAngleTerm() const
{
    return func_type == 2;
}

/** Return whether or not this dihedral is a angle/torsion cross term */
bool GromacsDihedral::isAngleTorsionCrossTerm() const
{
    return func_type == 11;
}

/** Return the string description of the function type for this dihedral */
QString GromacsDihedral::functionTypeString() const
{
    switch (func_type)
    {
        case 1:
            return "proper dihedral";
        case 2:
            return "improper dihedral";
        case 3:
            return "Ryckaert-Bellemans dihedral";
        case 4:
            return "periodic improper dihedral";
        case 5:
            return "Fourier dihedral";
        case 8:
            return "tabulated dihedral";
        case 9:
            return "proper dihedral (multiple)";
        case 10:
            return "restricted dihedral";
        case 11:
            return "combined bending torsion potential";
        default:
            return "unknown";
    }
}

/** Return whether or not this is a cosine-series dihedral */
bool GromacsDihedral::isCosine() const
{
    switch (func_type)
    {
    case 0: //null, so zero, which is a zero cosine
    case 1: //proper dihedral
    case 3: //Ryckaert-Bellemans dihedral can be converted to cosine.
    case 4: //improper dihedral
    case 5: //fourier dihedral
    case 9: //proper dihedral (multiple)
        return true;
    default:
        return false;
    }
}

/** Return all of the parameters for this dihedral */
QList<double> GromacsDihedral::parameters() const
{
    QList<double> params;

    for (int i=0; i<count(); ++i)
    {
        params.append( k[i] );
    }

    return params;
}

/** Return a string representation of this dihedral */
QString GromacsDihedral::toString() const
{
    if (func_type == 0)
        return QObject::tr("GromacsDihedral::null");
    else if (this->needsResolving())
    {
        return QObject::tr("GromacsDihedral( functionType() = %1, needsResolving )")
                .arg(functionTypeString());
    }
    else
    {
        QStringList params;

        for (int i=0; i<count(); ++i)
        {
            params.append( QString::number(k[i]) );
        }

        if (params.isEmpty())
            return QObject::tr("GromacsDihedral( functionType() = %1 )")
                    .arg(functionTypeString());
        else
            return QObject::tr("GromacsDihedral( functionType() = %1, parameters() = [ %2 ] )")
                    .arg(functionTypeString()).arg(params.join(", "));
    }
}

/** Return this function converted to a SireCAS::Expression using the passed
    symbol to represent the improper angle eta */
SireCAS::Expression GromacsDihedral::toImproperExpression(const SireCAS::Symbol &eta) const
{
    const double kj_per_mol_per_rad2 = ((kJ_per_mol) / (radian*radian)).value();
    const double deg = degree.value();

    if (func_type == 2)
    {
        //improper : k0 * (eta - eta0)^2
        const double k0 = k[1] * kj_per_mol_per_rad2;
        const double eta0 = k[0] * deg;

        return k0 * SireMaths::pow_2( eta - eta0 );
    }
    else
        throw SireError::incompatible_error( QObject::tr(
            "Can only create an improper function from a Gromacs dihedral "
            "of type '2', not from this type '%1'.").arg(this->toString()), CODELOC );
}

/** Return this function converted to a SireCAS::Expression using the passed
    symbol to represent the torsion (phi) and the angles either side of the
    torsion (theta0 and theta1) */
SireCAS::Expression GromacsDihedral::toAngleTorsionExpression(const Symbol &theta0,
                                                              const Symbol &theta1,
                                                              const Symbol &phi) const
{
    const double kj_per_mol = kJ_per_mol.value();

    if (func_type == 11)
    {
        //bending / torsion potential
        // k sin(theta0)^3 sin(theta1)^3 sum_(n=0,4) a_n cos(phi)^n

        const double k0 = 1.0; //this doesn't appear to be set in the gromacs file...
        const double a0 = k[0] * kj_per_mol;
        const double a1 = k[1] * kj_per_mol;
        const double a2 = k[2] * kj_per_mol;
        const double a3 = k[3] * kj_per_mol;
        const double a4 = k[4] * kj_per_mol;

        return k0 * SireMaths::pow_3(Sin(theta0)) * SireMaths::pow_3(Sin(theta1)) *
                ( a0 + (a1 * Cos(phi)) +
                       (a2 * SireMaths::pow_2(Cos(phi))) +
                       (a3 * SireMaths::pow_3(Cos(phi))) +
                       (a4 * SireMaths::pow_4(Cos(phi))) );
    }
    else
        throw SireError::incompatible_error( QObject::tr(
            "Can only create an improper function from a Gromacs dihedral "
            "of type '2', not from this type '%1'.").arg(this->toString()), CODELOC );
}

/** Return this function converted to a SireCAS::Expression using the passed symbol
    to represent the torsion size */
SireCAS::Expression GromacsDihedral::toExpression(const SireCAS::Symbol &phi) const
{
    const double kj_per_mol = kJ_per_mol.value();
    const double deg = degree.value();

    if (func_type == 1 or func_type == 4 or func_type == 9)
    {
        //standard periodic / non-periodic proper / improper dihedral
        // k ( 1 + cos( n phi - phi_s ))
        const double k0 = k[1] * kj_per_mol;
        const double phi_s = k[0] * deg;
        const double n = k[2];

        return k0 * ( 1.0 + Cos( (n*phi) - phi_s ) );
    }
    else if (func_type == 2)
    {
        throw SireError::incompatible_error( QObject::tr(
            "Cannot convert the 'improper angle' type Gromacs dihedral '%1' to an expression "
            "using only the torsion angle 'phi'")
                .arg(this->toString()), CODELOC );
    }
    else if (func_type == 3)
    {
        //Ryckaert-Bellemans function
        // Sum_(n=0,5) C_n ( cos(psi) )^n   where psi = phi - 180
        const double c0 = k[0] * kj_per_mol;
        const double c1 = k[1] * kj_per_mol;
        const double c2 = k[2] * kj_per_mol;
        const double c3 = k[3] * kj_per_mol;
        const double c4 = k[4] * kj_per_mol;
        const double c5 = k[5] * kj_per_mol;

        //Gromacs calculates dihedrals with a shift of 180 degrees??? (not included)
        const auto cos_psi = Cos(phi - SireMaths::pi);


        auto f = c0 + (c1*cos_psi) + (c2*cos_psi*cos_psi) +
                                   (c3*cos_psi*cos_psi*cos_psi) +
                                   (c4*cos_psi*cos_psi*cos_psi*cos_psi) +
                                   (c5*cos_psi*cos_psi*cos_psi*cos_psi*cos_psi);

        return f;
    }
    else if (func_type == 5)
    {
        //fourier dihedral
        // 0.5 [ C1(1+cos(phi)) + C2(1-cos(2phi)) + C3(1+cos(3phi)) + C4(1+cos(4phi)) ]
        const double c1 = k[0] * kj_per_mol;
        const double c2 = k[1] * kj_per_mol;
        const double c3 = k[2] * kj_per_mol;
        const double c4 = k[4] * kj_per_mol;

        return 0.5 * ( c1*(1 + Cos(phi)) +
                       c2*(1 - Cos(2*phi)) +
                       c3*(1 + Cos(3*phi)) +
                       c4*(1 + Cos(4*phi)) );
    }
    else if (func_type == 8)
    {
        throw SireError::unsupported( QObject::tr(
            "It is not possible to convert a tabulated gromacs dihedral "
            "into a SireCAS::Expression!"),
                CODELOC );
    }
    else if (func_type == 10)
    {
        //restricted dihedral potential : 0.5 k [ (cos(phi) - cos(phi_0))^2 / sin(phi)^2 ]
        double k0 = k[1] * kj_per_mol;
        double phi_0 = k[0] * deg;

        return 0.5 * k0 * ( (Cos(phi) - SireMaths::pow_2(std::cos(phi_0))) /
                            SireMaths::pow_2(Sin(phi)) );
    }
    else if (func_type == 11)
    {
        throw SireError::incompatible_error( QObject::tr(
                "Cannot convert the 'torsion/angle' function '%1' to an expression "
                "that depends only on the torsion angle.")
                    .arg(this->toString()), CODELOC );
    }
    else
        return SireCAS::Expression(0);
}

/** Return a hash for this bond */
uint GromacsDihedral::hash() const
{
    return  my_qHash(k[0]) | my_qHash(k[1]) | my_qHash(k[2]) | my_qHash(k[3]) |
            my_qHash(k[4]) | my_qHash(k[5]);
}

/** Return whether or not this parameter needs resolving */
bool GromacsDihedral::needsResolving() const
{
    return k[0] == unresolved_parameter_value or func_type == 0;
}

/** Return whether or not the parameters for this dihedral are resolved */
bool GromacsDihedral::isResolved() const
{
    return not needsResolving();
}

/** Assert that the parameters for this dihedral have been resolved */
void GromacsDihedral::assertResolved() const
{
    if (needsResolving())
        throw SireError::invalid_state( QObject::tr(
            "The parameters for this GromacsBond have not been resolved! %1")
                .arg(this->toString()), CODELOC );
}
