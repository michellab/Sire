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

#ifndef SIREMM_AMBERPARAMS_H
#define SIREMM_AMBERPARAMS_H

#include "SireBase/propertymap.h"
#include "SireBase/shareddatapointer.hpp"

#include "SireID/index.h"

#include "SireMol/molecule.h"
#include "SireMol/partialmolecule.h"
#include "SireMol/bondid.h"
#include "SireMol/angleid.h"
#include "SireMol/dihedralid.h"
#include "SireMol/atomidx.h"
#include "SireMol/mover.hpp"
#include "SireMol/molviewproperty.h"

SIRE_BEGIN_HEADER

namespace SireMM
{
class AmberParams;
class AmberBond;
class AmberAngle;
class AmberDihPart;
class AmberDihedral;
class AmberNB14;
}

QDataStream& operator<<(QDataStream&, const SireMM::AmberParams&);
QDataStream& operator>>(QDataStream&, SireMM::AmberParams&);

QDataStream& operator<<(QDataStream&, const SireMM::AmberBond&);
QDataStream& operator>>(QDataStream&, SireMM::AmberBond&);

QDataStream& operator<<(QDataStream&, const SireMM::AmberAngle&);
QDataStream& operator>>(QDataStream&, SireMM::AmberAngle&);

QDataStream& operator<<(QDataStream&, const SireMM::AmberDihPart&);
QDataStream& operator>>(QDataStream&, SireMM::AmberDihPart&);

QDataStream& operator<<(QDataStream&, const SireMM::AmberDihedral&);
QDataStream& operator>>(QDataStream&, SireMM::AmberDihedral&);

QDataStream& operator<<(QDataStream&, const SireMM::AmberNB14&);
QDataStream& operator>>(QDataStream&, SireMM::AmberNB14&);

namespace SireMol
{
    class BondID;
    class AngleID;
    class DihedralID;
    class ImproperID;
    class Molecule;
    class MoleculeData;
    class MoleculeInfoData;
}

namespace SireCAS
{
    class Expression;
}

namespace SireMM
{
using SireBase::PropertyMap;

using SireMol::BondID;
using SireMol::AngleID;
using SireMol::DihedralID;
using SireMol::ImproperID;
using SireMol::Molecule;
using SireMol::MoleculeData;
using SireMol::MoleculeInfoData;

/** This simple class holds Amber parameters for a bond */
class SIREMM_EXPORT AmberBond
{
private:
    double _k, _r0;

public:
    friend QDataStream& ::operator<<(QDataStream&, const AmberBond&);
    friend QDataStream& ::operator>>(QDataStream&, AmberBond&);

    AmberBond(double k=0, double r0=0) : _k(k), _r0(r0)
    {}
    
    AmberBond(const SireCAS::Expression &eqn);
    
    ~AmberBond();
    
    double k() const
    {
        return _k;
    }
    
    double r0() const
    {
        return _r0;
    }
    
    double operator[](int i) const
    {
        i = SireID::Index(i).map(2);
        
        if (i == 0)
            return _k;
        else
            return _r0;
    }
    
    bool operator==(const AmberBond &other) const
    {
        return _k == other._k and _r0 == other._r0;
    }
    
    double energy(double r) const
    {
        return _k * SireMaths::pow_2(r - _r0);
    }
    
    QString toString() const
    {
        return QObject::tr("AmberBond( k = %1, r0 = %2 )").arg(_k).arg(_r0);
    }
};

/** This simple class holds Amber parameters for an angle */
class SIREMM_EXPORT AmberAngle
{
private:
    double _k, _theta0;

public:
    friend QDataStream& ::operator<<(QDataStream&, const AmberAngle&);
    friend QDataStream& ::operator>>(QDataStream&, AmberAngle&);

    AmberAngle(double k=0, double theta0=0) : _k(k), _theta0(theta0)
    {}
    
    AmberAngle(const SireCAS::Expression &eqn);
    
    ~AmberAngle();
    
    double k() const
    {
        return _k;
    }
    
    double theta0() const
    {
        return _theta0;
    }
    
    double operator[](int i) const
    {
        i = SireID::Index(i).map(2);
        
        if (i == 0)
            return _k;
        else
            return _theta0;
    }
    
    bool operator==(const AmberAngle &other) const
    {
        return _k == other._k and _theta0 == other._theta0;
    }
    
    double energy(double theta) const
    {
        return _k * SireMaths::pow_2(theta - _theta0);
    }
    
    QString toString() const
    {
        return QObject::tr("AmberAngle( k = %1, theta0 = %2 )").arg(_k).arg(_theta0);
    }
};

/** This simple class holds Amber dihedral or improper parameter parts */
class SIREMM_EXPORT AmberDihPart
{
private:
    double _k, _periodicity, _phase;
    
public:
    friend QDataStream& ::operator<<(QDataStream&, const AmberDihPart&);
    friend QDataStream& ::operator>>(QDataStream&, AmberDihPart&);

    AmberDihPart(double k=0, double periodicity=0, double phase=0)
        : _k(k), _periodicity(periodicity), _phase(phase)
    {}
    
    ~AmberDihPart();
    
    double k() const
    {
        return _k;
    }
    
    double periodicity() const
    {
        return _periodicity;
    }
    
    double phase() const
    {
        return _phase;
    }
    
    double operator[](int i) const
    {
        i = SireID::Index(i).map(3);
        
        if (i == 0)
            return _k;
        else if (i == 1)
            return _periodicity;
        else
            return _phase;
    }
    
    bool operator==(const AmberDihPart &other) const
    {
        return _k == other._k and _periodicity == other._periodicity
                     and _phase == other._phase;
    }
    
    double energy(double phi) const
    {
        return _k * ( 1 + cos( _periodicity * ( phi - 0 ) - _phase ) );
    }
    
    QString toString() const
    {
        return QObject::tr("AmberDihPart( k = %1, periodicity = %2, phase = %3 )")
                .arg(_k).arg(_periodicity).arg(_phase);
    }
};

/** This simple class holds Amber dihedral or improper parameter */
class SIREMM_EXPORT AmberDihedral
{
private:
    QVector<AmberDihPart> _parts;
    
public:
    friend QDataStream& ::operator<<(QDataStream&, const AmberDihedral&);
    friend QDataStream& ::operator>>(QDataStream&, AmberDihedral&);

    AmberDihedral()
    {}
    
    AmberDihedral(AmberDihPart part)
    {
        _parts = QVector<AmberDihPart>(1, part);
    }
    
    AmberDihedral(const SireCAS::Expression &f);
    
    ~AmberDihedral();
    
    AmberDihedral& operator+=(const AmberDihPart &part)
    {
        _parts.append(part);
        return *this;
    }
    
    AmberDihedral operator+(const AmberDihPart &part) const
    {
        AmberDihedral ret(*this);
        ret += part;
        return *this;
    }
    
    bool operator==(const AmberDihedral &other) const
    {
        return _parts == other._parts;
    }
    
    const AmberDihPart& operator[](int i) const
    {
        i = SireID::Index(i).map(_parts.count());
        return _parts[i];
    }
    
    double energy(double phi) const
    {
        double total = 0;
        for (int i=0; i<_parts.count(); ++i)
        {
            total += _parts.constData()[i].energy(phi);
        }
        return total;
    }
    
    QString toString() const
    {
        QStringList s;
        for (int i=0; i<_parts.count(); ++i)
        {
            s.append( QObject::tr("k[%1] = %2, periodicity[%1] = %3, phase[%1] = %4")
                        .arg(i).arg(_parts[i].k()).arg(_parts[i].periodicity())
                        .arg(_parts[i].phase()) );
        }
        
        return QObject::tr("AmberDihedral( %1 )").arg(s.join(", "));
    }
};

/** This simple class holds Amber parameters for a 1-4 scale factor */
class SIREMM_EXPORT AmberNB14
{
private:
    double _cscl, _ljscl;

public:
    friend QDataStream& ::operator<<(QDataStream&, const AmberNB14&);
    friend QDataStream& ::operator>>(QDataStream&, AmberNB14&);
    
    AmberNB14(double cscl=0, double ljscl=0) : _cscl(cscl), _ljscl(ljscl)
    {}
    
    ~AmberNB14();
    
    double cscl() const
    {
        return _cscl;
    }
    
    double ljscl() const
    {
        return _ljscl;
    }
    
    double operator[](int i) const
    {
        i = SireID::Index(i).map(2);
        
        if (i == 0)
            return _cscl;
        else
            return _ljscl;
    }
    
    bool operator==(const AmberNB14 &other) const
    {
        return _cscl == other._cscl and _ljscl == other._ljscl;
    }
    
    QString toString() const
    {
        return QObject::tr("AmberNB14( cscl = %1, ljscl = %2 )").arg(_cscl).arg(_ljscl);
    }
};


/** This class stores AMBER bonded force field parameters for 
    a collection of bonds, angles and dihedrals

    @author Julien Michel / Christopher Woods
 */
class SIREMOL_EXPORT AmberParams
    : public SireBase::ConcreteProperty<AmberParams,SireMol::MoleculeProperty>
{

friend QDataStream& ::operator<<(QDataStream&, const SireMM::AmberParams&);
friend QDataStream& ::operator>>(QDataStream&, SireMM::AmberParams&);

public:
    AmberParams();
    AmberParams(const SireMol::MoleculeView &molecule);
    AmberParams(const SireMol::MoleculeInfoData &molinfo);
    
    AmberParams(const AmberParams &other);

    ~AmberParams();

    static const char* typeName();

    AmberParams& operator=(const AmberParams &other);
    
    bool operator==(const AmberParams &other) const;
    bool operator!=(const AmberParams &other) const;

    const SireMol::MoleculeInfoData& info() const;

    bool isCompatibleWith(const SireMol::MoleculeInfoData &molinfo) const;

    void add(const BondID &bond, const double &k, const double &ro);
    void remove(const BondID &bond);
    
    AmberBond getParams(const BondID &bond);
    QList<BondID> getAllBonds();

    void add(const AngleID &angle, const double &k, const double &theta0);
    void remove(const AngleID &angle);

    AmberAngle getParams(const AngleID &angle);
    QList<AngleID> getAllAngles();

    void add(const DihedralID &dihedral, const double &v, const double &periodicity, const double &phase);
    void remove(const DihedralID &dihedral);

    AmberDihedral getParams(const DihedralID &dihedral);
    QList<DihedralID> getAllDihedrals();

    void add(const ImproperID &improper, const double &v, const double &periodicity, const double &phase);
    void remove(const ImproperID &improper);

    AmberDihedral getParams(const ImproperID &improper);
    QList<ImproperID> getAllImpropers();

    void add14Pair(const BondID &pair, const double &cscl, const double &ljscl);
    void remove14Pair(const BondID &pair);

    AmberNB14 get14PairParams(const BondID &pair);
    QList<BondID> getAll14Pairs();
  
 private:
    /** The molecule that this flexibility operates on */
    SireBase::SharedDataPointer<SireMol::MoleculeInfoData> molinfo;

    /**A Hash of force constants and equilibrium bond lengths for bonds **/
    QHash<BondID,AmberBond> bonds;

    /**A Hash of force constants and equilibrium bond angles for angles **/
    QHash<AngleID,AmberAngle> angles;

    /**A Hash of torsional parameters for dihedrals **/
    QHash<DihedralID,AmberDihedral> dihedrals;

    /**A Hash of torsional parameters for impropers **/
    QHash<ImproperID,AmberDihedral> impropers;

    /**A Hash of coulombic and lennard jones scale factors for 1,4 pairs**/
    QHash<BondID,AmberNB14> nb14pairs;
};

}

Q_DECLARE_METATYPE( SireMM::AmberParams )

SIRE_EXPOSE_CLASS( SireMM::AmberParams )
SIRE_EXPOSE_CLASS( SireMM::AmberBond )
SIRE_EXPOSE_CLASS( SireMM::AmberAngle )
SIRE_EXPOSE_CLASS( SireMM::AmberDihPart )
SIRE_EXPOSE_CLASS( SireMM::AmberDihedral )
SIRE_EXPOSE_CLASS( SireMM::AmberNB14 )

SIRE_END_HEADER

#endif
