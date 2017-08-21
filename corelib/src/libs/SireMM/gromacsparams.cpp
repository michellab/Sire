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

#include "SireError/errors.h"

#include "SireUnits/units.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMM;
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
        << typ._v << typ._w << typ.combiningRuleString();
    
    return ds;
}

QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, GromacsAtomType &typ)
{
    VersionID v = readHeader(ds, r_atomtype);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        double chg, mass;
        QString ptyp, rule;
        
        sds >> typ._typ >> mass >> chg >> ptyp >> typ._v >> typ._w >> rule;
        
        typ._mass = mass * g_per_mol;
        typ._chg = chg * mod_electron;
        typ._ptyp = GromacsAtomType::toParticleType(ptyp);
        typ._rule = GromacsAtomType::toCombiningRule(rule);
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
                  _v(0), _w(0), _rule(GromacsAtomType::UNKNOWN_RULE)
{}

/** Construct passing in all parameters */
GromacsAtomType::GromacsAtomType(QString atom_type,
                                 SireUnits::Dimension::MolarMass mass,
                                 SireUnits::Dimension::Charge charge,
                                 PARTICLE_TYPE particle_type,
                                 double v, double w, COMBINING_RULE combining_rule)
                : _typ(atom_type), _mass(mass), _chg(charge),
                  _ptyp(particle_type), _v(v), _w(w), _rule(combining_rule)
{
    assertSane();
}

/** Copy constructor */
GromacsAtomType::GromacsAtomType(const GromacsAtomType &other)
                : _typ(other._typ), _mass(other._mass), _chg(other._chg),
                  _ptyp(other._ptyp), _v(other._v), _w(other._w), _rule(other._rule)
{}

/** Destructor */
GromacsAtomType::~GromacsAtomType()
{}

/** Create this object from the passed line from a Gromacs topology file */
GromacsAtomType GromacsAtomType::fromGromacsTopLine(const QString &line)
{
    return GromacsAtomType();
}

/** Copy assignment operator */
GromacsAtomType& GromacsAtomType::operator=(const GromacsAtomType &other)
{
    if (this != &other)
    {
        _typ = other._typ;
        _mass = other._mass;
        _chg = other._chg;
        _ptyp = other._ptyp;
        _rule = other._rule;
        _v = other._v;
        _w = other._w;
    }
    
    return *this;
}

/** Comparison operator */
bool GromacsAtomType::operator==(const GromacsAtomType &other) const
{
    return _typ == other._typ and _mass == other._mass and _chg == other._chg and
           _ptyp == other._ptyp and _rule == other._rule and
           _v == other._v and _w == other._w;
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
    return QObject::tr("GromacsAtomType( atomType() = %1, mass() = %2 g mol-1, "
                          "charge() = %3 |e|, "
                          "particleType() = %4 )")
                          .arg(atomType()).arg(mass().to(g_per_mol))
                          .arg(charge().to(mod_electron))
                          .arg(particleTypeString());
}

/** Return this parameter as a line that is in the right format to go into
    a Gromacs topology file */
QString GromacsAtomType::toGromacsTopLine() const
{
    //need to be converted to the right units
    double mass = _mass;
    double chg = _chg;

    QString letter = this->particleTypeLetter();

    //need to be converted depending on the combining rule...
    double v = _v;
    double w = _w;

    return QString("%1 %2 %3 %4 %5 %6").arg(_typ).arg(mass).arg(chg).arg(letter).arg(v).arg(w);
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

/** Return a string version of the combining rule */
QString GromacsAtomType::combiningRuleString() const
{
    switch(_rule)
    {
    case ARITHMETIC:
        return "arithmetic";
    case GEOMETRIC:
        return "geometric";
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

/** Convert the passed string to a Gromacs combining rule type. Use 'ok' to see if this
    worked correctly */
GromacsAtomType::COMBINING_RULE GromacsAtomType::toCombiningRule(const QString &word, bool *ok)
{
    if (word == "1")
    {
        if (ok) *ok = true;
        return GEOMETRIC;
    }
    else if (word == "2")
    {
        if (ok) *ok = true;
        return ARITHMETIC;
    }
    
    QString lword = word.toLower();
    
    if (lword == "geometric")
    {
        if (ok) *ok = true;
        return GEOMETRIC;
    }
    else if (lword == "arithmetic")
    {
        if (ok) *ok = true;
        return ARITHMETIC;
    }
    else
    {
        if (ok) *ok = false;
        return UNKNOWN_RULE;
    }
}
