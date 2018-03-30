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
#include "SireMol/element.h"

#include "SireMM/ljparameter.h"

#include "SireUnits/dimensions.h"

#include "SireCAS/symbol.h"
#include "SireCAS/expression.h"

SIRE_BEGIN_HEADER

namespace SireMM
{
class GromacsAtomType;
class GromacsBond;
class GromacsAngle;
class GromacsDihedral;
}

QDataStream& operator<<(QDataStream&, const SireMM::GromacsAtomType&);
QDataStream& operator>>(QDataStream&, SireMM::GromacsAtomType&);

QDataStream& operator<<(QDataStream&, const SireMM::GromacsBond&);
QDataStream& operator>>(QDataStream&, SireMM::GromacsBond&);

QDataStream& operator<<(QDataStream&, const SireMM::GromacsAngle&);
QDataStream& operator>>(QDataStream&, SireMM::GromacsAngle&);

QDataStream& operator<<(QDataStream&, const SireMM::GromacsDihedral&);
QDataStream& operator>>(QDataStream&, SireMM::GromacsDihedral&);

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

    GromacsAtomType();
    
    GromacsAtomType(QString atom_type,
                    SireUnits::Dimension::MolarMass mass);
    
    GromacsAtomType(QString atom_type,
                    SireUnits::Dimension::MolarMass mass,
                    SireUnits::Dimension::Charge charge,
                    PARTICLE_TYPE particle_type,
                    const LJParameter &ljparam,
                    const SireMol::Element &element = SireMol::Element(0));

    GromacsAtomType(QString atom_type, QString bond_type,
                    SireUnits::Dimension::MolarMass mass,
                    SireUnits::Dimension::Charge charge,
                    PARTICLE_TYPE particle_type,
                    const LJParameter &ljparam,
                    const SireMol::Element &element = SireMol::Element(0));
    
    GromacsAtomType(const GromacsAtomType &other);
    
    ~GromacsAtomType();
    
    GromacsAtomType& operator=(const GromacsAtomType &other);
    
    bool operator==(const GromacsAtomType &other) const;
    bool operator!=(const GromacsAtomType &other) const;
    
    static const char* typeName();
    const char* what() const;
    
    bool isNull() const;
    
    bool hasMassOnly() const;
    
    QString toString() const;
    
    QString atomType() const;
    QString bondType() const;
    SireUnits::Dimension::MolarMass mass() const;
    SireUnits::Dimension::Charge charge() const;
    SireMM::LJParameter ljParameter() const;
    
    SireMol::Element element() const;
    
    PARTICLE_TYPE particleType() const;
    
    QString particleTypeString() const;
    QString particleTypeLetter() const;
    
    static PARTICLE_TYPE toParticleType(const QString &word, bool *ok=0);

private:
    void assertSane() const;

    QString _typ;
    QString _btyp;
    SireUnits::Dimension::MolarMass _mass;
    SireUnits::Dimension::Charge _chg;
    SireMM::LJParameter _lj;
    PARTICLE_TYPE _ptyp;
    SireMol::Element _elem;
};

/** This class holds all of the information about a Gromacs Bond

    @author Christopher Woods
*/
class SIREMM_EXPORT GromacsBond
{

friend QDataStream& ::operator<<(QDataStream&, const GromacsBond&);
friend QDataStream& ::operator>>(QDataStream&, GromacsBond&);

public:
    GromacsBond();
    GromacsBond(int function_type);
    GromacsBond(int function_type,
                double k0, double k1=0, double k2=0, double k3=0);
    GromacsBond(int function_type, const QList<double> &params);
    
    GromacsBond(const SireCAS::Expression &bond, const SireCAS::Symbol &R);
    
    GromacsBond(const GromacsBond &other);
    
    ~GromacsBond();
    
    GromacsBond& operator=(const GromacsBond &other);
    
    bool operator==(const GromacsBond &other) const;
    bool operator!=(const GromacsBond &other) const;
    
    bool operator<(const GromacsBond &other) const;
    bool operator<=(const GromacsBond &other) const;
    
    bool operator>(const GromacsBond &other) const;
    bool operator>=(const GromacsBond &other) const;
    
    static const char* typeName();
    const char* what() const;
    
    double operator[](int i) const;
    
    double at(int i) const;
    
    int count() const;
    int size() const;
    
    int functionType() const;
    QString functionTypeString() const;
    
    bool needsResolving() const;
    bool isResolved() const;
    
    bool isSimple() const;
    
    bool isHarmonic() const;
    
    void assertResolved() const;
    
    QList<double> parameters() const;
    
    QString toString() const;
    SireCAS::Expression toExpression(const SireCAS::Symbol &R) const;
    
    bool atomsAreBonded() const;
    
    uint hash() const;

private:
    /** Space to hold up to 4 parameters */
    double k[4];

    /** The function type */
    qint32 func_type;
};

/** This class holds all of the information about a Gromacs Angle

    @author Christopher Woods
*/
class SIREMM_EXPORT GromacsAngle
{

friend QDataStream& ::operator<<(QDataStream&, const GromacsAngle&);
friend QDataStream& ::operator>>(QDataStream&, GromacsAngle&);

public:
    GromacsAngle();
    GromacsAngle(int function_type);
    GromacsAngle(int function_type,
                 double k0, double k1=0, double k2=0, double k3=0, double k4=0, double k5=0);
    GromacsAngle(int function_type, const QList<double> &params);
    
    GromacsAngle(const SireCAS::Expression &angle, const SireCAS::Symbol &theta);
    
    GromacsAngle(const GromacsAngle &other);
    
    ~GromacsAngle();
    
    GromacsAngle& operator=(const GromacsAngle &other);
    
    bool operator==(const GromacsAngle &other) const;
    bool operator!=(const GromacsAngle &other) const;
    
    bool operator<(const GromacsAngle &other) const;
    bool operator<=(const GromacsAngle &other) const;
    
    bool operator>(const GromacsAngle &other) const;
    bool operator>=(const GromacsAngle &other) const;
    
    static const char* typeName();
    const char* what() const;
    
    double operator[](int i) const;
    
    double at(int i) const;
    
    int count() const;
    int size() const;
    
    int functionType() const;
    QString functionTypeString() const;
    
    QList<double> parameters() const;
    
    bool isSimple() const;
    
    bool isBondBondCrossTerm() const;
    bool isBondAngleCrossTerm() const;
    
    GromacsBond toBondTerm() const;
    GromacsAngle toAngleTerm() const;
    
    QString toString() const;
    
    bool isHarmonic() const;
    
    bool isResolved() const;
    bool needsResolving() const;
    
    void assertResolved() const;
    
    SireCAS::Expression toExpression(const SireCAS::Symbol &theta) const;
    
    SireCAS::Expression toBondAngleExpression(const SireCAS::Symbol &r,
                                              const SireCAS::Symbol &theta) const;
    
    SireCAS::Expression toBondBondExpression(const SireCAS::Symbol &r01,
                                             const SireCAS::Symbol &r12,
                                             const SireCAS::Symbol &r02) const;
    
    uint hash() const;

private:
    /** Space to hold up to 6 parameters */
    double k[6];

    /** The function type */
    qint32 func_type;
};


/** This class holds all of the information about a Gromacs Dihedral

    @author Christopher Woods
*/
class SIREMM_EXPORT GromacsDihedral
{

friend QDataStream& ::operator<<(QDataStream&, const GromacsDihedral&);
friend QDataStream& ::operator>>(QDataStream&, GromacsDihedral&);

public:
    GromacsDihedral();
    GromacsDihedral(int function_type);
    GromacsDihedral(int function_type,
                    double k0, double k1=0, double k2=0, double k3=0, double k4=0, double k5=0);
    GromacsDihedral(int function_type, const QList<double> &params);
    
    GromacsDihedral(const SireCAS::Expression &dihedral, const SireCAS::Symbol &phi);
    
    static QList<GromacsDihedral> construct(const SireCAS::Expression &dihedral,
                                            const SireCAS::Symbol &phi);
    static QList<GromacsDihedral> constructImproper(
                                            const SireCAS::Expression &dihedral,
                                            const SireCAS::Symbol &phi);
    
    GromacsDihedral(const GromacsDihedral &other);
    
    ~GromacsDihedral();
    
    GromacsDihedral& operator=(const GromacsDihedral &other);
    
    bool operator==(const GromacsDihedral &other) const;
    bool operator!=(const GromacsDihedral &other) const;
    
    bool operator<(const GromacsDihedral &other) const;
    bool operator<=(const GromacsDihedral &other) const;
    
    bool operator>(const GromacsDihedral &other) const;
    bool operator>=(const GromacsDihedral &other) const;
    
    static const char* typeName();
    const char* what() const;
    
    double operator[](int i) const;
    
    double at(int i) const;
    
    int count() const;
    int size() const;
    
    int functionType() const;
    QString functionTypeString() const;
    
    QList<double> parameters() const;
    
    bool isCosine() const;
    
    bool isSimple() const;
    
    bool isImproperAngleTerm() const;
    bool isAngleTorsionCrossTerm() const;

    bool isResolved() const;
    bool needsResolving() const;
    
    void assertResolved() const;

    QString toString() const;
    
    SireCAS::Expression toExpression(const SireCAS::Symbol &phi) const;
    
    SireCAS::Expression toImproperExpression(const SireCAS::Symbol &eta) const;
    
    SireCAS::Expression toAngleTorsionExpression(const SireCAS::Symbol &theta0,
                                                 const SireCAS::Symbol &theta1,
                                                 const SireCAS::Symbol &phi) const;
    
    uint hash() const;

private:
    /** Space to hold up to 6 parameters */
    double k[6];

    /** The function type */
    qint32 func_type;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Return the atom type name */
inline QString GromacsAtomType::atomType() const
{
    return _typ;
}

/** Return the bond type name - this is normally the same as the atom type,
    but is different for some forcefields, like OPLS */
inline QString GromacsAtomType::bondType() const
{
    return _btyp;
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

/** Return the Lennard Jones parameter */
inline SireMM::LJParameter GromacsAtomType::ljParameter() const
{
    return _lj;
}

/** Return whether or not this is a null atom type */
inline bool GromacsAtomType::isNull() const
{
    return _typ.isNull();
}

/** Return whether or not this atom type has only the mass specified */
inline bool GromacsAtomType::hasMassOnly() const
{
    return _ptyp == UNKNOWN_TYPE and _chg.value() == 0 and _lj.isDummy();
}

/** Return the element type of this type */
inline SireMol::Element GromacsAtomType::element() const
{
    return _elem;
}

#endif

}

Q_DECLARE_METATYPE( SireMM::GromacsAtomType )
Q_DECLARE_METATYPE( SireMM::GromacsBond )
Q_DECLARE_METATYPE( SireMM::GromacsAngle )
Q_DECLARE_METATYPE( SireMM::GromacsDihedral )

SIRE_EXPOSE_CLASS( SireMM::GromacsAtomType )
SIRE_EXPOSE_CLASS( SireMM::GromacsBond )
SIRE_EXPOSE_CLASS( SireMM::GromacsAngle )
SIRE_EXPOSE_CLASS( SireMM::GromacsDihedral )

SIRE_END_HEADER

#endif
