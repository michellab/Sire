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
using namespace SireMol;
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
