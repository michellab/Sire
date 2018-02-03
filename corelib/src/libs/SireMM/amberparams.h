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
#include "SireMol/moleculeinfo.h"
#include "SireMol/partialmolecule.h"
#include "SireMol/bondid.h"
#include "SireMol/angleid.h"
#include "SireMol/dihedralid.h"
#include "SireMol/improperid.h"
#include "SireMol/atomidx.h"
#include "SireMol/mover.hpp"
#include "SireMol/molviewproperty.h"
#include "SireMol/atomcharges.h"
#include "SireMol/atommasses.h"
#include "SireMol/atomelements.h"
#include "SireMol/atomradii.h"
#include "SireMol/moleculeinfo.h"

#include "SireMM/cljnbpairs.h"
#include "SireMM/atomljs.h"

#include "SireCAS/symbol.h"

#include <QHash>
#include <QPair>

SIRE_BEGIN_HEADER

namespace SireMM
{
class AmberParams;
class AmberBond;
class AmberAngle;
class AmberDihPart;
class AmberDihedral;
class AmberNB14;
class AmberNBDihPart;
class TwoAtomFunctions;
class ThreeAtomFunctions;
class FourAtomFunctions;
class CLJNBPairs;
class CLJScaleFactor;
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

QDataStream& operator<<(QDataStream&, const SireMM::AmberNBDihPart&);
QDataStream& operator>>(QDataStream&, SireMM::AmberNBDihPart&);

namespace SireMol
{
    class BondID;
    class AngleID;
    class DihedralID;
    class ImproperID;
    class Molecule;
    class MoleculeData;
    class Connectivity;
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
public:
    friend QDataStream& ::operator<<(QDataStream&, const AmberBond&);
    friend QDataStream& ::operator>>(QDataStream&, AmberBond&);

    AmberBond(double k=0, double r0=0);
    
    AmberBond(const SireCAS::Expression &eqn, const SireCAS::Symbol &R);
    
    AmberBond(const AmberBond &other);
    
    ~AmberBond();
    
    static const char* typeName();
    const char* what() const;
    
    double k() const;
    double r0() const;
    
    double operator[](int i) const;
    
    AmberBond& operator=(const AmberBond &other);
    
    bool operator==(const AmberBond &other) const;
    bool operator!=(const AmberBond &other) const;
    
    bool operator<(const AmberBond &other) const;
    bool operator<=(const AmberBond &other) const;
    
    bool operator>(const AmberBond &other) const;
    bool operator>=(const AmberBond &other) const;
    
    double energy(double r) const;
    
    QString toString() const;
    SireCAS::Expression toExpression(const SireCAS::Symbol &R) const;

    uint hash() const;

private:
    double _k, _r0;
};

/** This simple class holds Amber parameters for an angle */
class SIREMM_EXPORT AmberAngle
{
public:
    friend QDataStream& ::operator<<(QDataStream&, const AmberAngle&);
    friend QDataStream& ::operator>>(QDataStream&, AmberAngle&);

    AmberAngle(double k=0, double theta0=0);
    
    AmberAngle(const SireCAS::Expression &eqn, const SireCAS::Symbol &THETA);
    
    AmberAngle(const AmberAngle &other);
    
    ~AmberAngle();
    
    static const char* typeName();
    const char* what() const;
    
    double k() const;
    
    double theta0() const;
    
    double operator[](int i) const;
    
    AmberAngle& operator=(const AmberAngle &other);
    
    bool operator==(const AmberAngle &other) const;
    bool operator!=(const AmberAngle &other) const;
    
    bool operator<(const AmberAngle &other) const;
    bool operator<=(const AmberAngle &other) const;
    
    bool operator>(const AmberAngle &other) const;
    bool operator>=(const AmberAngle &other) const;
    
    double energy(double theta) const;
    
    SireCAS::Expression toExpression(const SireCAS::Symbol &THETA) const;
    
    QString toString() const;

    uint hash() const;

private:
    double _k, _theta0;
};

/** This simple class holds Amber dihedral or improper parameter parts */
class SIREMM_EXPORT AmberDihPart
{
public:
    friend QDataStream& ::operator<<(QDataStream&, const AmberDihPart&);
    friend QDataStream& ::operator>>(QDataStream&, AmberDihPart&);

    AmberDihPart(double k=0, double periodicity=0, double phase=0);
    
    AmberDihPart(const AmberDihPart &other);
    
    ~AmberDihPart();
    
    static const char* typeName();
    const char* what() const;
    
    double k() const;
    
    double periodicity() const;
    
    double phase() const;
    
    double operator[](int i) const;
    
    AmberDihPart& operator=(const AmberDihPart &other);
    
    bool operator==(const AmberDihPart &other) const;
    bool operator!=(const AmberDihPart &other) const;
    
    bool operator<(const AmberDihPart &other) const;
    bool operator<=(const AmberDihPart &other) const;
    
    bool operator>(const AmberDihPart &other) const;
    bool operator>=(const AmberDihPart &other) const;
    
    double energy(double phi) const;
    
    QString toString() const;

    uint hash() const;

private:
    double _k, _periodicity, _phase;
};

/** This simple class holds Amber dihedral or improper parameter */
class SIREMM_EXPORT AmberDihedral
{
public:
    friend QDataStream& ::operator<<(QDataStream&, const AmberDihedral&);
    friend QDataStream& ::operator>>(QDataStream&, AmberDihedral&);

    AmberDihedral();
    
    AmberDihedral(AmberDihPart part);
    
    AmberDihedral(const SireCAS::Expression &f, const SireCAS::Symbol &PHI);
    
    AmberDihedral(const AmberDihedral &other);
    
    ~AmberDihedral();
    
    static const char* typeName();
    const char* what() const;
    
    AmberDihedral& operator+=(const AmberDihPart &part);
    AmberDihedral operator+(const AmberDihPart &part) const;
    
    AmberDihedral& operator=(const AmberDihedral &other);
    
    bool operator==(const AmberDihedral &other) const;
    bool operator!=(const AmberDihedral &other) const;
    
    AmberDihPart operator[](int i) const;
    
    QVector<AmberDihPart> terms() const;
    
    double energy(double phi) const;
    
    SireCAS::Expression toExpression(const SireCAS::Symbol &PHI) const;
    
    QString toString() const;

    uint hash() const;

private:
    QVector<AmberDihPart> _parts;
};

/** This simple class holds Amber parameters for a 1-4 scale factor */
class SIREMM_EXPORT AmberNB14
{
public:
    friend QDataStream& ::operator<<(QDataStream&, const AmberNB14&);
    friend QDataStream& ::operator>>(QDataStream&, AmberNB14&);
    
    AmberNB14(double cscl=0, double ljscl=0);
    
    AmberNB14(const AmberNB14 &other);
    
    ~AmberNB14();
    
    static const char* typeName();
    const char* what() const;
    
    double cscl() const;
    
    double ljscl() const;
    
    double operator[](int i) const;
    
    AmberNB14& operator=(const AmberNB14 &other);
    
    bool operator==(const AmberNB14 &other) const;
    bool operator!=(const AmberNB14 &other) const;
    
    bool operator<(const AmberNB14 &other) const;
    bool operator<=(const AmberNB14 &other) const;
    
    bool operator>(const AmberNB14 &other) const;
    bool operator>=(const AmberNB14 &other) const;
    
    QString toString() const;
    
    SireMM::CLJScaleFactor toScaleFactor() const;

private:
    double _cscl, _ljscl;
};

/** This simple class holds a combination of the AmberDihPart with AmberNB14,
    as these are combined together in some amber file types */
class SIREMM_EXPORT AmberNBDihPart
{
public:
    friend QDataStream& ::operator<<(QDataStream&, const AmberNBDihPart&);
    friend QDataStream& ::operator>>(QDataStream&, AmberNBDihPart&);

    AmberNBDihPart();
    AmberNBDihPart(const AmberDihPart &dih, const AmberNB14 &nb14);
    
    AmberNBDihPart(const AmberNBDihPart &other);
    
    ~AmberNBDihPart();
    
    static const char* typeName();
    const char* what() const;

    const AmberDihPart& parameter() const;
    const AmberNB14& scl() const;
    
    AmberNBDihPart& operator=(const AmberNBDihPart &other);
    
    bool operator==(const AmberNBDihPart &other) const;
    bool operator!=(const AmberNBDihPart &other) const;
    
    bool operator<(const AmberNBDihPart &other) const;
    bool operator<=(const AmberNBDihPart &other) const;
    
    bool operator>(const AmberNBDihPart &other) const;
    bool operator>=(const AmberNBDihPart &other) const;
    
    QString toString() const;

    uint hash() const;

private:
    /** The dihedral parameters */
    AmberDihPart dih;
    
    /** The 14 scale factor */
    AmberNB14 nbscl;
};

/** This class stores AMBER bonded force field parameters for 
    a collection of bonds, angles, dihedrals, impropers 
    and 1-4 scaling factors.

    @author Julien Michel / Christopher Woods
 */
class SIREMM_EXPORT AmberParams
    : public SireBase::ConcreteProperty<AmberParams,SireMol::MoleculeProperty>
{

friend QDataStream& ::operator<<(QDataStream&, const SireMM::AmberParams&);
friend QDataStream& ::operator>>(QDataStream&, SireMM::AmberParams&);

public:
    AmberParams();
    AmberParams(const SireMol::MoleculeView &molecule,
                const SireBase::PropertyMap &map = SireBase::PropertyMap());
                
    AmberParams(const SireMol::MoleculeInfo &molinfo);
    AmberParams(const SireMol::MoleculeInfoData &molinfo);
    
    AmberParams(const AmberParams &other);

    ~AmberParams();

    static const char* typeName();

    AmberParams& operator=(const AmberParams &other);
    
    bool operator==(const AmberParams &other) const;
    bool operator!=(const AmberParams &other) const;

    AmberParams& operator+=(const AmberParams &other);
    
    AmberParams operator+(const AmberParams &other) const;

    SireMol::MoleculeInfo info() const;

    bool isCompatibleWith(const SireMol::MoleculeInfoData &molinfo) const;

    QString toString() const;

    void add(const SireMol::AtomID &atom,
             SireUnits::Dimension::Charge charge,
             SireUnits::Dimension::MolarMass mass,
             const SireMol::Element &element,
             const SireMM::LJParameter &ljparam,
             const QString &amber_type,
             SireUnits::Dimension::Length born_radius,
             double screening_parameter,
             const QString &treechain);

    SireMol::AtomCharges charges() const;
    SireMol::AtomMasses masses() const;
    SireMol::AtomElements elements() const;
    SireMM::AtomLJs ljs() const;
    SireMol::AtomStringProperty amberTypes() const;
    SireMol::AtomRadii gbRadii() const;
    SireMol::AtomFloatProperty gbScreening() const;
    SireMol::AtomStringProperty treeChains() const;

    void setRadiusSet(const QString &radius_set);
    QString radiusSet() const;

    void setExcludedAtoms(const CLJNBPairs &excluded_atoms);
    CLJNBPairs excludedAtoms() const;

    SireMol::Connectivity connectivity() const;

    void add(const BondID &bond, double k, double r0, bool includes_hydrogen);
    void remove(const BondID &bond);
    
    AmberBond getParameter(const BondID &bond) const;
    QHash< BondID,QPair<AmberBond,bool> > bonds() const;
    TwoAtomFunctions bondFunctions() const;
    TwoAtomFunctions bondFunctions(const SireCAS::Symbol &R) const;

    void add(const AngleID &angle, double k, double theta0, bool includes_hydrogen);
    void remove(const AngleID &angle);

    AmberAngle getParameter(const AngleID &angle) const;
    QHash< AngleID,QPair<AmberAngle,bool> > angles() const;
    ThreeAtomFunctions angleFunctions() const;
    ThreeAtomFunctions angleFunctions(const SireCAS::Symbol &THETA) const;

    void add(const DihedralID &dihedral, double k, double periodicity, double phase,
             bool includes_hydrogen);
    void remove(const DihedralID &dihedral);

    AmberDihedral getParameter(const DihedralID &dihedral) const;
    QHash< DihedralID,QPair<AmberDihedral,bool> > dihedrals() const;
    FourAtomFunctions dihedralFunctions() const;
    FourAtomFunctions dihedralFunctions(const SireCAS::Symbol &PHI) const;

    void add(const ImproperID &improper, double v, double periodicity, double phase,
             bool includes_hydrogen);
    void remove(const ImproperID &improper);

    AmberDihedral getParameter(const ImproperID &improper) const;
    QHash< ImproperID,QPair<AmberDihedral,bool> > impropers() const;
    FourAtomFunctions improperFunctions() const;
    FourAtomFunctions improperFunctions(const SireCAS::Symbol &PHI) const;

    void addNB14(const BondID &pair, double cscl, double ljscl);
    void removeNB14(const BondID &pair);

    AmberNB14 getNB14(const BondID &pair) const;
    QHash<BondID,AmberNB14> nb14s() const;
    CLJNBPairs cljScaleFactors() const;
  
    QStringList validate() const;
    QStringList validateAndFix();
  
    void setPropertyMap(const PropertyMap &map);
    const PropertyMap& propertyMap() const;
  
    void updateFrom(const MoleculeView &molview);

    BondID convert(const BondID &bond) const;
    AngleID convert(const AngleID &angle) const;
    DihedralID convert(const DihedralID &dihedral) const;
    ImproperID convert(const ImproperID &improper) const;
  
private: 
    void _pvt_createFrom(const MoleculeData &moldata);
    void _pvt_updateFrom(const MoleculeData &moldata);

    void getAmberBondsFrom(const TwoAtomFunctions &funcs);
    void getAmberAnglesFrom(const ThreeAtomFunctions &funcs);
    void getAmberDihedralsFrom(const FourAtomFunctions &funcs);
    void getAmberImpropersFrom(const FourAtomFunctions &funcs);
    
    void getAmberNBsFrom(const CLJNBPairs &nbpairs, const FourAtomFunctions &funcs);
 
    /** The molecule that this flexibility operates on */
    SireMol::MoleculeInfo molinfo;

    /** All of the atom charges */
    SireMol::AtomCharges amber_charges;
    
    /** All of the atom LJs */
    SireMM::AtomLJs amber_ljs;
    
    /** All of the atom masses */
    SireMol::AtomMasses amber_masses;
    
    /** All of the atom elements */
    SireMol::AtomElements amber_elements;
    
    /** All of the amber atom types */
    SireMol::AtomStringProperty amber_types;

    /** All of the amber atom Born radii */
    SireMol::AtomRadii born_radii;
    
    /** All of the amber GB screening parameters */
    SireMol::AtomFloatProperty amber_screens;

    /** All of the amber treechain classifications for all atoms */
    SireMol::AtomStringProperty amber_treechains;

    /** The excluded atoms in the molecule (atom pairs between
        which a nonbonded calculation is not evaluated) */
    SireMM::CLJNBPairs exc_atoms;

    /** A hash of force constants and equilibrium bond lengths for bonds **/
    QHash< BondID,QPair<AmberBond,bool> > amber_bonds;

    /** A hash of force constants and equilibrium bond angles for angles **/
    QHash< AngleID,QPair<AmberAngle,bool> > amber_angles;

    /** A hash of torsional parameters for dihedrals **/
    QHash< DihedralID,QPair<AmberDihedral,bool> > amber_dihedrals;

    /** A hash of torsional parameters for impropers **/
    QHash< ImproperID,QPair<AmberDihedral,bool> > amber_impropers;

    /** A hash of coulombic and lennard jones scale factors for 1,4 pairs**/
    QHash<BondID,AmberNB14> amber_nb14s;
    
    /** The radius set chosen by LEAP when assigning the Born radii */
    QString radius_set;
    
    /** The property map used (if any) to identify the properties that
        hold the amber parameters */
    SireBase::PropertyMap propmap;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Return the force constant in kcal mol-1 A-2 */
inline double AmberBond::k() const
{
    return _k;
}

/** Return the equilibrium bond length, in angstroms */
inline double AmberBond::r0() const
{
    return _r0;
}

inline uint my_qHash(double key)
{
    return ::qHash( *(reinterpret_cast<const ulong*>(&key)) );
}

/** Return a hash for the bond */
inline uint AmberBond::hash() const
{
    return my_qHash(_k) | my_qHash(_r0);
}

/** Return a hash of this bond */
inline uint qHash(const AmberBond &param)
{
    return param.hash();
}

/** Comparison operator */
inline bool AmberBond::operator<(const AmberBond &other) const
{
    if (_k < other._k)
    {
        return true;
    }
    else if (_k == other._k)
    {
        return _r0 < other._r0;
    }
    else
        return false;
}

/** Return the force constant in kcal mol-1 radian-2 */
inline double AmberAngle::k() const
{
    return _k;
}

/** Return the equilibrium angle, in radians */
inline double AmberAngle::theta0() const
{
    return _theta0;
}

/** Return a hash for the angle */
inline uint AmberAngle::hash() const
{
    return my_qHash(_k) | my_qHash(_theta0);
}

/** Return a hash of this angle */
inline uint qHash(const AmberAngle &param)
{
    return param.hash();
}

/** Comparison operator */
inline bool AmberAngle::operator<(const AmberAngle &other) const
{
    if (_k < other._k)
    {
        return true;
    }
    else if (_k == other._k)
    {
        return _theta0 < other._theta0;
    }
    else
        return false;
}

/** Return the force constant, in kcal mol-1 */
inline double AmberDihPart::k() const
{
    return _k;
}

/** Return the periodicity */
inline double AmberDihPart::periodicity() const
{
    return _periodicity;
}

/** Return the phase */
inline double AmberDihPart::phase() const
{
    return _phase;
}

/** Return a hash for the dihedral part */
inline uint AmberDihPart::hash() const
{
    return my_qHash(_k) | my_qHash(_periodicity) | my_qHash(_phase);
}

/** Return a hash of this dihedral */
inline uint qHash(const AmberDihPart &param)
{
    return param.hash();
}

/** Return a hash for the dihedral */
inline uint AmberDihedral::hash() const
{
    uint h = 0;
    
    for (int i=0; i<_parts.count(); ++i)
    {
        h |= _parts[i].hash();
    }

    return h;
}

/** Comparison operator */
inline bool AmberDihPart::operator<(const AmberDihPart &other) const
{
    if (_k < other._k)
    {
        return true;
    }
    else if (_k == other._k)
    {
        if (_periodicity < other._periodicity)
        {
            return true;
        }
        else if (_periodicity == other._periodicity)
        {
            return _phase < other._phase;
        }
        else
            return false;
    }
    else
        return false;
}

/** Return a hash of the AmberNBDihPart */
inline uint AmberNBDihPart::hash() const
{
    return dih.hash();
}

inline uint qHash(const AmberNBDihPart &param)
{
    return param.hash();
}

/** Return an array of all of the terms in the dihedral */
inline QVector<AmberDihPart> AmberDihedral::terms() const
{
    if (_parts.isEmpty())
    {
        return QVector<AmberDihPart>( 1, AmberDihPart() );
    }
    else
    {
        return _parts;
    }
}

/** Return a hash of this bond */
inline uint qHash(const AmberDihedral &param)
{
    return param.hash();
}

/** Return the 14 electrostatic scaling factor */
inline double AmberNB14::cscl() const
{
    return _cscl;
}

/** Return the 14 LJ scaling factor */
inline double AmberNB14::ljscl() const
{
    return _ljscl;
}

/** Return a hash of all of the bond parameters */
inline QHash< BondID,QPair<AmberBond,bool> > AmberParams::bonds() const
{
    return amber_bonds;
}

/** Return a hash of all of the angle parameters */
inline QHash< AngleID,QPair<AmberAngle,bool> > AmberParams::angles() const
{
    return amber_angles;
}

/** Return a hash of all of the dihedral parameters */
inline QHash< DihedralID,QPair<AmberDihedral,bool> > AmberParams::dihedrals() const
{
    return amber_dihedrals;
}

/** Return a hash of all of the improper parameters */
inline QHash< ImproperID,QPair<AmberDihedral,bool> > AmberParams::impropers() const
{
    return amber_impropers;
}

/** Return a hash of all of the 1-4 nonbonded parameters */
inline QHash<BondID,AmberNB14> AmberParams::nb14s() const
{
    return amber_nb14s;
}

/** Return the parameter part of the class */
inline const AmberDihPart& AmberNBDihPart::parameter() const
{
    return dih;
}

/** Return the NB14 scale factor */
inline const AmberNB14& AmberNBDihPart::scl() const
{
    return nbscl;
}

#endif // SIRE_SKIP_INLINE_FUNCTIONS

}

Q_DECLARE_METATYPE( SireMM::AmberParams )
Q_DECLARE_METATYPE( SireMM::AmberBond )
Q_DECLARE_METATYPE( SireMM::AmberAngle )
Q_DECLARE_METATYPE( SireMM::AmberDihPart )
Q_DECLARE_METATYPE( SireMM::AmberDihedral )
Q_DECLARE_METATYPE( SireMM::AmberNB14 )
Q_DECLARE_METATYPE( SireMM::AmberNBDihPart )

SIRE_EXPOSE_CLASS( SireMM::AmberParams )
SIRE_EXPOSE_CLASS( SireMM::AmberBond )
SIRE_EXPOSE_CLASS( SireMM::AmberAngle )
SIRE_EXPOSE_CLASS( SireMM::AmberDihPart )
SIRE_EXPOSE_CLASS( SireMM::AmberDihedral )
SIRE_EXPOSE_CLASS( SireMM::AmberNB14 )

SIRE_END_HEADER

#endif
