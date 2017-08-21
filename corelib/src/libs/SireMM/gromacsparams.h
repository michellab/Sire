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

#ifndef SIREMM_GROMACSPARAMS_H
#define SIREMM_GROMACSPARAMS_H

#include "SireMol/molviewproperty.h"

#include "SireUnits/dimensions.h"

SIRE_BEGIN_HEADER

namespace SireMM
{
class GromacsAtomType;
}

QDataStream& operator<<(QDataStream&, const SireMM::GromacsAtomType&);
QDataStream& operator>>(QDataStream&, SireMM::GromacsAtomType&);

namespace SireMM
{

/** This class represents a Gromacs format atom type. This
    combines particle type information with charge, mass,
    and LJ parameter info
    
    @author Christopher Woods
*/
class SIREMM_EXPORT GromacsAtomType
{

friend QDataStream& ::operator<<(QDataStream&, const GromacsAtomType&);
friend QDataStream& ::operator>>(QDataStream&, GromacsAtomType&);

public:
    enum PARTICLE_TYPE
    {
        UNKNOWN_TYPE = 0,
        ATOM = 1,
        SHELL = 2,
        VIRTUAL = 3
    };

    enum COMBINING_RULE
    {
        UNKNOWN_RULE = 0,
        ARITHMETIC = 1,
        GEOMETRIC = 2
    };

    GromacsAtomType();
    
    GromacsAtomType(QString atom_type,
                    SireUnits::Dimension::MolarMass mass,
                    SireUnits::Dimension::Charge charge,
                    PARTICLE_TYPE particle_type,
                    double v, double w, COMBINING_RULE combining_rule);
    
    GromacsAtomType(const GromacsAtomType &other);
    
    ~GromacsAtomType();

    static GromacsAtomType fromGromacsTopLine(const QString &line);
    
    GromacsAtomType& operator=(const GromacsAtomType &other);
    
    bool operator==(const GromacsAtomType &other) const;
    bool operator!=(const GromacsAtomType &other) const;
    
    static const char* typeName();
    const char* what() const;
    
    QString toString() const;
    
    QString toGromacsTopLine() const;
    
    QString atomType() const;
    SireUnits::Dimension::MolarMass mass() const;
    SireUnits::Dimension::Charge charge() const;
    
    PARTICLE_TYPE particleType() const;
    
    QString particleTypeString() const;
    QString particleTypeLetter() const;

    QString combiningRuleString() const;
    
    static PARTICLE_TYPE toParticleType(const QString &word, bool *ok=0);
    static COMBINING_RULE toCombiningRule(const QString &word, bool *ok=0);
    
    double V() const;
    double W() const;
    
    COMBINING_RULE combiningRule() const;

private:
    void assertSane() const;

    QString _typ;
    SireUnits::Dimension::MolarMass _mass;
    SireUnits::Dimension::Charge _chg;
    PARTICLE_TYPE _ptyp;
    double _v, _w;
    COMBINING_RULE _rule;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Return the atom type name */
inline QString GromacsAtomType::atomType() const
{
    return _typ;
}

/** Return the atom mass */
inline SireUnits::Dimension::MolarMass GromacsAtomType::mass() const
{
    return _mass;
}

/** Return the atomic charge */
inline SireUnits::Dimension::Charge GromacsAtomType::charge() const
{
    return _chg;
}

/** Return the type of the particle */
inline GromacsAtomType::PARTICLE_TYPE GromacsAtomType::particleType() const
{
    return _ptyp;
}

/** Return the Lennard Jones V value (meaning depends on the combining rules) */
inline double GromacsAtomType::V() const
{
    return _v;
}

/** Return the Lennard Jones W value (meaning depends on the combining rules) */
inline double GromacsAtomType::W() const
{
    return _w;
}

/** Return the combining rule to use with this parameter */
inline GromacsAtomType::COMBINING_RULE GromacsAtomType::combiningRule() const
{
    return _rule;
}

#endif

}

Q_DECLARE_METATYPE( SireMM::GromacsAtomType )

SIRE_EXPOSE_CLASS( SireMM::GromacsAtomType )

SIRE_END_HEADER

#endif
