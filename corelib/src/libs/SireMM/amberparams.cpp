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

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMol;
using namespace SireCAS;
using namespace SireMM;
using namespace SireStream;

///////////
/////////// Implementation of AmberBond
///////////

QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, const AmberBond &bond)
{
    ds << bond._k << bond._r0;
    return ds;
}

QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, AmberBond &bond)
{
    ds >> bond._k >> bond._r0;
    return ds;
}

/** Construct with the passed bond constant and equilibrium bond length */
AmberBond::AmberBond(double k, double r0) : _k(k), _r0(r0)
{}

/** Construct from the passed expression */
AmberBond::AmberBond(const Expression &f, const Symbol &R) : _k(0), _r0(0)
{}

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

QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, const AmberAngle &angle)
{
    ds << angle._k << angle._theta0;
    return ds;
}

QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, AmberAngle &angle)
{
    ds >> angle._k >> angle._theta0;
    return ds;
}

AmberAngle::AmberAngle(double k, double theta0) : _k(k), _theta0(theta0)
{}

AmberAngle::AmberAngle(const Expression &f, const Symbol &theta) : _k(0), _theta0(0)
{}

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

QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, const AmberDihPart &dih)
{
    ds << dih._k << dih._periodicity << dih._phase;
    return ds;
}

QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, AmberDihPart &dih)
{
    ds >> dih._k >> dih._periodicity >> dih._phase;
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

QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, const AmberDihedral &dih)
{
    ds << dih._parts;
    return ds;
}

QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, AmberDihedral &dih)
{
    ds >> dih._parts;
    return ds;
}

AmberDihedral::AmberDihedral()
{}

AmberDihedral::AmberDihedral(AmberDihPart part)
{
    _parts = QVector<AmberDihPart>(1, part);
}

AmberDihedral::AmberDihedral(const Expression &f, const Symbol &phi)
{}

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

const AmberDihPart& AmberDihedral::operator[](int i) const
{
    i = SireID::Index(i).map(_parts.count());
    return _parts[i];
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

QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, const AmberNB14 &nb)
{
    ds << nb._cscl << nb._ljscl;
    return ds;
}

QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, AmberNB14 &nb)
{
    ds >> nb._cscl >> nb._ljscl;
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

bool AmberNB14::operator==(const AmberNB14 &other) const
{
    return _cscl == other._cscl and _ljscl == other._ljscl;
}

bool AmberNB14::operator!=(const AmberNB14 &other) const
{
    return not operator==(other);
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
        << amberparam.amber_types << amberparam.exc_atoms
        << amberparam.amber_bonds
        << amberparam.amber_angles << amberparam.amber_dihedrals
        << amberparam.amber_impropers << amberparam.amber_nb14s
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
            >> amberparam.amber_types >> amberparam.exc_atoms
            >> amberparam.amber_bonds
            >> amberparam.amber_angles >> amberparam.amber_dihedrals
            >> amberparam.amber_impropers >> amberparam.amber_nb14s
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
AmberParams::AmberParams(const MoleculeView &mol)
            : ConcreteProperty<AmberParams,MoleculeProperty>(),
              molinfo(mol.data().info())
{}

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
              amber_types(other.amber_types), exc_atoms(other.exc_atoms),
              amber_bonds(other.amber_bonds),
              amber_angles(other.amber_angles),
              amber_dihedrals(other.amber_dihedrals),
              amber_impropers(other.amber_impropers),
              amber_nb14s(other.amber_nb14s)
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
        exc_atoms = other.exc_atoms;
        amber_bonds = other.amber_bonds;
        amber_angles = other.amber_angles;
        amber_dihedrals = other.amber_dihedrals;
        amber_impropers = other.amber_impropers;
        amber_nb14s = other.amber_nb14s;
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
          and exc_atoms == other.exc_atoms
          and amber_bonds == other.amber_bonds
          and amber_angles == other.amber_angles
          and amber_dihedrals == other.amber_dihedrals
          and amber_impropers == other.amber_impropers
          and amber_nb14s == other.amber_nb14s);
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

/** Set the atom parameters for the specified atom to the provided values */
void AmberParams::add(const AtomID &atom,
                      SireUnits::Dimension::Charge charge,
                      SireUnits::Dimension::MolarMass mass,
                      const SireMol::Element &element,
                      const SireMM::LJParameter &ljparam,
                      const QString &amber_type)
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
    }
    
    amber_charges.set(idx, charge);
    amber_ljs.set(idx, ljparam);
    amber_masses.set(idx, mass);
    amber_elements.set(idx, element);
    amber_types.set(idx, amber_type);
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
    if (not this->isCompatibleWith(other.info()))
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

    return *this;
}

/** Return a combination of the two passed AmberParams */
AmberParams AmberParams::operator+(const AmberParams &other) const
{
    AmberParams ret(*this);
    
    ret += other;
    
    return ret;
}
