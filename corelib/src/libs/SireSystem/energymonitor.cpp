/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  * 
  *  Copyright (C) 2012  Christopher Woods
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

#include "SireSystem/energymonitor.h"
#include "SireSystem/system.h"

#include "SireMol/partialmolecule.h"
#include "SireMol/atomcoords.h"
#include "SireMol/atomcharges.h"
#include "SireMM/atomljs.h"
#include "SireMM/ljpair.h"
#include "SireMol/mgname.h"

#include "SireBase/array2d.hpp"

#include "SireUnits/units.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireSystem;
using namespace SireMol;
using namespace SireMaths;
using namespace SireMM;
using namespace SireUnits;
using namespace SireUnits::Dimension;
using namespace SireBase;
using namespace SireStream;

using std::pair;

static const RegisterMetaType<EnergyMonitor> r_nrgmonitor;

QDataStream &operator<<(QDataStream &ds, 
                                          const EnergyMonitor &nrgmonitor)
{
    writeHeader(ds, r_nrgmonitor, 3);
    
    SharedDataStream sds(ds);
    
    sds << nrgmonitor.grp0 << nrgmonitor.grp1 
        << nrgmonitor.asgn0 << nrgmonitor.asgn1
        << nrgmonitor.accum
        << nrgmonitor.coul_nrgs << nrgmonitor.lj_nrgs
        << nrgmonitor.alpha_component << nrgmonitor.alfa
        << nrgmonitor.shift_delta
        << nrgmonitor.coulomb_power
        << static_cast<const SystemMonitor&>(nrgmonitor);
        
    return ds;
}

QDataStream &operator>>(QDataStream &ds,
                                          EnergyMonitor &nrgmonitor)
{
    VersionID v = readHeader(ds, r_nrgmonitor);
    
    if (v == 3)
    {
        SharedDataStream sds(ds);
        
        sds >> nrgmonitor.grp0 >> nrgmonitor.grp1 
            >> nrgmonitor.asgn0 >> nrgmonitor.asgn1
            >> nrgmonitor.accum
            >> nrgmonitor.coul_nrgs >> nrgmonitor.lj_nrgs
            >> nrgmonitor.alpha_component >> nrgmonitor.alfa
            >> nrgmonitor.shift_delta >> nrgmonitor.coulomb_power
            >> static_cast<SystemMonitor&>(nrgmonitor);
    }
    else if (v == 2)
    {
        SharedDataStream sds(ds);

        nrgmonitor.alfa = 0;
        nrgmonitor.shift_delta = 0;
        nrgmonitor.coulomb_power = 0;
        nrgmonitor.alpha_component = Symbol();        

        sds >> nrgmonitor.grp0 >> nrgmonitor.grp1 
            >> nrgmonitor.asgn0 >> nrgmonitor.asgn1
            >> nrgmonitor.accum
            >> nrgmonitor.coul_nrgs >> nrgmonitor.lj_nrgs
            >> static_cast<SystemMonitor&>(nrgmonitor);
    }
    else if (v == 1)
    {
        SharedDataStream sds(ds);

        nrgmonitor.alfa = 0;
        nrgmonitor.shift_delta = 0;
        nrgmonitor.coulomb_power = 0;
        nrgmonitor.alpha_component = Symbol();        

        sds >> nrgmonitor.grp0 >> nrgmonitor.grp1 >> nrgmonitor.accum
            >> nrgmonitor.coul_nrgs >> nrgmonitor.lj_nrgs
            >> static_cast<SystemMonitor&>(nrgmonitor);
    }
    else
        throw version_error(v, "1,2", r_nrgmonitor, CODELOC);
        
    return ds;
}

/** Null constructor */
EnergyMonitor::EnergyMonitor() 
              : ConcreteProperty<EnergyMonitor,SystemMonitor>(),
                alfa(0), shift_delta(0), coulomb_power(0)
{}

/** Construct to monitor the energies between all pairs of molecule views in the
    two passed groups. This will accumulate the average and standard deviation
    of each of the energies */
EnergyMonitor::EnergyMonitor(const MoleculeGroup &group0, 
                             const MoleculeGroup &group1)
              : ConcreteProperty<EnergyMonitor,SystemMonitor>(),
                grp0(group0), grp1(group1), accum( AverageAndStddev() ),
                alfa(0), shift_delta(0), coulomb_power(0)
{}
            
/** Construct to monitor the energies between all pairs of molecule views in 
    the two passed groups, accumulating the energies using the passed
    accumulator */
EnergyMonitor::EnergyMonitor(const MoleculeGroup &group0, 
                             const MoleculeGroup &group1,
                             const SireMaths::Accumulator &accumulator)
              : ConcreteProperty<EnergyMonitor,SystemMonitor>(),
                grp0(group0), grp1(group1), accum(accumulator),
                alfa(0), shift_delta(0), coulomb_power(0)
{}

/** Construct to monitor the energies between all pairs of molecule views in the
    two passed groups. This will accumulate the average and standard deviation
    of each of the energies */
EnergyMonitor::EnergyMonitor(const MoleculeGroup &group0,
                             const IDAssigner &group1)
              : ConcreteProperty<EnergyMonitor,SystemMonitor>(),
                grp0(group0), asgn1(group1), accum( AverageAndStddev() ),
                alfa(0), shift_delta(0), coulomb_power(0)
{}
            
/** Construct to monitor the energies between all pairs of molecule views in 
    the two passed groups, accumulating the energies using the passed
    accumulator */
EnergyMonitor::EnergyMonitor(const MoleculeGroup &group0, 
                             const IDAssigner &group1,
                             const SireMaths::Accumulator &accumulator)
              : ConcreteProperty<EnergyMonitor,SystemMonitor>(),
                grp0(group0), asgn1(group1), accum(accumulator),
                alfa(0), shift_delta(0), coulomb_power(0)
{}

/** Construct to monitor the energies between all pairs of molecule views in the
    two passed groups. This will accumulate the average and standard deviation
    of each of the energies */
EnergyMonitor::EnergyMonitor(const IDAssigner &group0,
                             const MoleculeGroup &group1)
              : ConcreteProperty<EnergyMonitor,SystemMonitor>(),
                grp1(group1), asgn0(group0), accum( AverageAndStddev() ),
                alfa(0), shift_delta(0), coulomb_power(0)
{}
            
/** Construct to monitor the energies between all pairs of molecule views in 
    the two passed groups, accumulating the energies using the passed
    accumulator */
EnergyMonitor::EnergyMonitor(const IDAssigner &group0, 
                             const MoleculeGroup &group1,
                             const SireMaths::Accumulator &accumulator)
              : ConcreteProperty<EnergyMonitor,SystemMonitor>(),
                grp1(group1), asgn0(group0), accum(accumulator),
                alfa(0), shift_delta(0), coulomb_power(0)
{}

/** Construct to monitor the energies between all pairs of molecule views in the
    two passed groups. This will accumulate the average and standard deviation
    of each of the energies */
EnergyMonitor::EnergyMonitor(const IDAssigner &group0,
                             const IDAssigner &group1)
              : ConcreteProperty<EnergyMonitor,SystemMonitor>(),
                asgn0(group0), asgn1(group1), accum( AverageAndStddev() ),
                alfa(0), shift_delta(0), coulomb_power(0)
{}
            
/** Construct to monitor the energies between all pairs of molecule views in 
    the two passed groups, accumulating the energies using the passed
    accumulator */
EnergyMonitor::EnergyMonitor(const IDAssigner &group0, 
                             const IDAssigner &group1,
                             const SireMaths::Accumulator &accumulator)
              : ConcreteProperty<EnergyMonitor,SystemMonitor>(),
                asgn0(group0), asgn1(group1), accum(accumulator),
                alfa(0), shift_delta(0), coulomb_power(0)
{}

/** Copy constructor */
EnergyMonitor::EnergyMonitor(const EnergyMonitor &other)
              : ConcreteProperty<EnergyMonitor,SystemMonitor>(other),
                grp0(other.grp0), grp1(other.grp1), 
                asgn0(other.asgn0), asgn1(other.asgn1),
                accum(other.accum),
                coul_nrgs(other.coul_nrgs), lj_nrgs(other.lj_nrgs),
                alpha_component(other.alpha_component),
                alfa(other.alfa), shift_delta(other.shift_delta),
                coulomb_power(other.coulomb_power) 
{}

/** Destructor */
EnergyMonitor::~EnergyMonitor()
{}

/** Copy assignment operator */
EnergyMonitor& EnergyMonitor::operator=(const EnergyMonitor &other)
{
    if (this != &other)
    {
        grp0 = other.grp0;
        grp1 = other.grp1;
        asgn0 = other.asgn0;
        asgn1 = other.asgn1;
        accum = other.accum;
        coul_nrgs = other.coul_nrgs;
        lj_nrgs = other.lj_nrgs;
        alpha_component = other.alpha_component;
        alfa = other.alfa;
        shift_delta = other.shift_delta;
        coulomb_power = other.coulomb_power;
        SystemMonitor::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool EnergyMonitor::operator==(const EnergyMonitor &other) const
{
    return this == &other or
           (grp0 == other.grp0 and grp1 == other.grp1 and
            asgn0 == other.asgn0 and asgn1 == other.asgn1 and
            accum == other.accum and 
            coul_nrgs == other.coul_nrgs and
            lj_nrgs == other.lj_nrgs and
            alpha_component == other.alpha_component and
            alfa == other.alfa and
            shift_delta == other.shift_delta and
            coulomb_power == other.coulomb_power and
            SystemMonitor::operator==(other));
}

/** Comparison operator */
bool EnergyMonitor::operator!=(const EnergyMonitor &other) const
{
    return not EnergyMonitor::operator==(other);
}

/** Return the typename of the class */
const char* EnergyMonitor::typeName()
{
    return QMetaType::typeName( qMetaTypeId<EnergyMonitor>() );
}

/** Return the array of all accumulated coulomb energies */
Array2D<AccumulatorPtr> EnergyMonitor::coulombEnergies() const
{
    return coul_nrgs;
}

/** Return the array of all accumulated LJ energies */
Array2D<AccumulatorPtr> EnergyMonitor::ljEnergies() const
{
    return lj_nrgs;
}

/** Return the array of the first group of molecule views in the same order as they
    appear in the arrays of energies */
QVector<SireMol::PartialMolecule> EnergyMonitor::views0() const
{
    if (asgn0.isNull())
    {
        int n = grp0->nViews();
    
        QVector<PartialMolecule> views(n);
    
        for (int i=0; i<n; ++i)
        {
            views[i] = grp0->viewAt(i);
        }
    
        return views;
    }
    else
    {
        return asgn0.read().asA<IDAssigner>().identifiedMolecules();
    }
}

/** Return the array of the second group of molecule views in the same order as they
    appear in the arrays of energies */
QVector<SireMol::PartialMolecule> EnergyMonitor::views1() const
{
    if (asgn1.isNull())
    {
        int n = grp1->nViews();
    
        QVector<PartialMolecule> views(n);
    
        for (int i=0; i<n; ++i)
        {
            views[i] = grp1->viewAt(i);
        }
    
        return views;
    }
    else
    {
        return asgn1.read().asA<IDAssigner>().identifiedMolecules();
    }
}

/** Return the molecule group from which "views0" are drawn. Note that this
    will return the molecule group used by "assigner0" if an assigner is
    used to choose views */
const SireMol::MoleculeGroup& EnergyMonitor::group0() const
{
    if (asgn0.isNull())
        return grp0.read();
    else
        return asgn0.read().asA<IDAssigner>().moleculeGroup();
}

/** Return the molecule group from which "views1" are drawn. Note that this
    will return the molecule group used by "assigner1" if an assigner is
    used to choose views */
const SireMol::MoleculeGroup& EnergyMonitor::group1() const
{
    if (asgn1.isNull())
        return grp1.read();
    else
        return asgn1.read().asA<IDAssigner>().moleculeGroup();
}

/** Return the assigner used to select "views0". Note that this will
    raise an exception if an assigner is not used to choose these views 
    
    \throw SireError::unavailable_resource
*/
const IDAssigner& EnergyMonitor::assigner0() const
{
    if (asgn0.isNull())
        throw SireError::unavailable_resource( QObject::tr(
                "An IDAssigner is not used to select the views for \"views0\"."),
                    CODELOC );
                    
    return asgn0.read().asA<IDAssigner>();
}

/** Return the assigner used to select "views1". Note that this will
    raise an exception if an assigner is not used to choose these views 
    
    \throw SireError::unavailable_resource
*/
const IDAssigner& EnergyMonitor::assigner1() const
{
    if (asgn1.isNull())
        throw SireError::unavailable_resource( QObject::tr(
                "An IDAssigner is not used to select the views for \"views1\"."),
                    CODELOC );
                    
    return asgn1.read().asA<IDAssigner>();
}

/** Clear all statistics */
void EnergyMonitor::clearStatistics()
{
    coul_nrgs = Array2D<AccumulatorPtr>();
    lj_nrgs = Array2D<AccumulatorPtr>();
}

inline pair<double,double> getCLJEnergy(
                const QVector<Vector> &coords0, const QVector<Charge> &chgs0,
                const QVector<LJParameter> &ljs0,
                const QVector<Vector> &coords1, const QVector<Charge> &chgs1,
                const QVector<LJParameter> &ljs1)
{
    const int nats0 = coords0.count();
    const int nats1 = coords1.count();
    
    double cnrg = 0;
    double ljnrg = 0;
    
    bool arithmetic_combining_rules = true;
    
    for (int i=0; i<nats0; ++i)
    {
        const Vector &coord0 = coords0.at(i);
        const Charge &chg0 = chgs0.at(i);
        const LJParameter &lj0 = ljs0.at(i);
        
        for (int j=0; j<nats1; ++j)
        {
            const Vector &coord1 = coords1.at(j);
            const Charge &chg1 = chgs1.at(j);
            const LJParameter &lj1 = ljs1.at(j);
            
            LJPair ljpair;

            if (arithmetic_combining_rules)
                ljpair = LJPair::arithmetic(lj0, lj1);
            else
                ljpair = LJPair::geometric(lj0, lj1);
            
            double one_over_r2 = Vector::invDistance2(coord0, coord1);
            double one_over_r = std::sqrt(one_over_r2);
            double one_over_r6 = one_over_r2 * one_over_r2 * one_over_r2;
            double one_over_r12 = one_over_r6 * one_over_r6;
            
            ljnrg += ljpair.A() * one_over_r12 - ljpair.B() * one_over_r6;
            
            cnrg += chg0.value() * chg1.value() * one_over_r 
                        * one_over_four_pi_eps0;
        }
    }
    
    return pair<double,double>(cnrg,ljnrg);
}

inline pair<double,double> getSoftCLJEnergy(
                const QVector<Vector> &coords0, const QVector<Charge> &chgs0,
                const QVector<LJParameter> &ljs0,
                const QVector<Vector> &coords1, const QVector<Charge> &chgs1,
                const QVector<LJParameter> &ljs1,
                double alpha, double shift_delta, int coulomb_power)
{
    //this uses the following potentials
    //           Zacharias and McCammon, J. Chem. Phys., 1994, and also,
    //           Michel et al., JCTC, 2007
    //
    //  V_{LJ}(r) = 4 epsilon [ ( sigma^12 / (delta*sigma + r^2)^6 ) - 
    //                          ( sigma^6  / (delta*sigma + r^2)^3 ) ]
    //
    //  delta = shift_delta * alpha
    //
    //  V_{coul}(r) = (1-alpha)^n q_i q_j / 4 pi eps_0 (alpha+r^2)^(1/2)
    //
    // This contrasts to Rich T's LJ softcore function, which was;
    //
    //  V_{LJ}(r) = 4 epsilon [ (sigma^12 / (alpha^m sigma^6 + r^6)^2) - 
    //                          (sigma^6  / (alpha^m sigma^6 + r^6) ) ]

    double one_minus_alfa_to_n = 1;
    const double delta = shift_delta * alpha;

    if (coulomb_power != 0)
        one_minus_alfa_to_n = SireMaths::pow(1 - alpha, coulomb_power);

    const int nats0 = coords0.count();
    const int nats1 = coords1.count();
    
    double cnrg = 0;
    double ljnrg = 0;
    
    bool arithmetic_combining_rules = true;
    
    for (int i=0; i<nats0; ++i)
    {
        const Vector &coord0 = coords0.at(i);
        const Charge &chg0 = chgs0.at(i);
        const LJParameter &lj0 = ljs0.at(i);
        
        for (int j=0; j<nats1; ++j)
        {
            const Vector &coord1 = coords1.at(j);
            const Charge &chg1 = chgs1.at(j);
            const LJParameter &lj1 = ljs1.at(j);
            
            LJPair ljpair;

            if (arithmetic_combining_rules)
                ljpair = LJPair::arithmetic(lj0, lj1);
            else
                ljpair = LJPair::geometric(lj0, lj1);
            
            double r2 = Vector::distance2(coord0, coord1);

            const double shift = ljpair.sigma() * delta;
            double lj_denom = r2 + shift;
            lj_denom = lj_denom * lj_denom * lj_denom;

            const double sig2 = ljpair.sigma() * ljpair.sigma();
            const double sig6 = sig2 * sig2 * sig2;
                        
            const double sig6_over_denom = sig6 / lj_denom;
            const double sig12_over_denom2 = sig6_over_denom *
                                             sig6_over_denom;
            
            ljnrg += ljpair.epsilon() * (sig12_over_denom2 - sig6_over_denom);
            
            cnrg += chg0.value() * chg1.value() * one_over_four_pi_eps0 / 
                          std::sqrt(alpha + r2);
        }
    }

    cnrg *= one_minus_alfa_to_n;
    ljnrg *= 4;
    
    return pair<double,double>(cnrg,ljnrg);
}

/** Return whether or not this monitor uses a soft-core potential to
    calculate the CLJ energy between the molecules in views0() and the
    molecules in views1() */
bool EnergyMonitor::usesSoftCore() const
{
    return (shift_delta != 0 or coulomb_power != 0) and
           (alfa != 0);
}

/** Set the system component symbol used to get the value of alpha 
    if using a soft-core potential. Note that this will overwrite
    any explicitly-set value of alpha */
void EnergyMonitor::setAlphaComponent(const Symbol &component)
{
    alpha_component = component;
    alfa = 1;
}

/** Explicitly set the value of alpha used if a soft-core potential is used.
    This clears any set alpha component symbol. */ 
void EnergyMonitor::setAlpha(double alpha)
{
    alfa = alpha;
    alpha_component = Symbol();
}

/** Set the shift delta parameter used by the soft-core potential */
void EnergyMonitor::setShiftDelta(double delta)
{
    shift_delta = delta;
}

/** Set the coulomb power parameter used by the soft-core potential */
void EnergyMonitor::setCoulombPower(int power)
{
    if (power >= 0)
        coulomb_power = power;
}

/** Return the shift delta parameter if a soft-core potential is used.
    This returns 0 if a LJ shifting term is not used */
double EnergyMonitor::shiftDelta() const
{
    return shift_delta;
}

/** Return the value of alpha (either the explicitly set value, or 
    the last value used when calculating the energy if an alpha
    component is used) */
double EnergyMonitor::alpha() const
{
    return alfa;
}

/** Return the coulomb power, if extra coulomb softening is used.
    This returns 0 if coulomb softening is not used */
int EnergyMonitor::coulombPower() const
{
    return coulomb_power;
}

/** Accumulate energies from the passed system */
void EnergyMonitor::monitor(System &system)
{
    EnergyMonitor old_state(*this);

    try
    {
        if (asgn0.isNull())
        {
            if (system.contains(grp0->number()))
                grp0 = system[grp0->number()];
        }
        else
        {
            asgn0.edit().asA<IDAssigner>().update(system);
        }
            
        if (asgn1.isNull())
        {
            if (system.contains(grp1->number()))
                grp1 = system[grp1->number()];
        }
        else
        {
            asgn1.edit().asA<IDAssigner>().update(system);
        }

        //see if we need to update the alpha value
        if (not alpha_component.isNull())
        {
            alfa = system.componentValue(alpha_component);
        }

        QVector<PartialMolecule> v0 = this->views0();
        QVector<PartialMolecule> v1 = this->views1();

        // extract the charge, LJ and coordinates of all of the views
        const PropertyName charge_prop("charge");
        const PropertyName lj_prop("LJ");
        const PropertyName coords_prop("coordinates");
        
        const int n0 = v0.count();
        const int n1 = v1.count();

        //has the number of views changed?
        if (n0 != coul_nrgs.nRows() or n1 != coul_nrgs.nColumns())
        {
            Array2D<AccumulatorPtr> coul_nrgs2(n0,n1,accum);
            Array2D<AccumulatorPtr> lj_nrgs2(n0,n1,accum);
            
            for (int i=0; i < qMin(n0,coul_nrgs.nRows()); ++i)
            {
                for (int j=0; j < qMin(n1,coul_nrgs.nColumns()); ++j)
                {
                    coul_nrgs2(i,j) = coul_nrgs(i,j);
                    lj_nrgs2(i,j) = lj_nrgs(i,j);
                }
            }
            
            coul_nrgs = coul_nrgs2;
            lj_nrgs = lj_nrgs2;
        }

        QVector< QVector<Vector> > grp0_coords(n0);
        QVector< QVector<Charge> > grp0_chgs(n0);
        QVector< QVector<LJParameter> > grp0_ljs(n0);
        
        QVector< QVector<Vector> > grp1_coords(n1);
        QVector< QVector<Charge> > grp1_chgs(n1);
        QVector< QVector<LJParameter> > grp1_ljs(n1);

        for (int i=0; i<n0; ++i)
        {
            const PartialMolecule &mol0 = v0.constData()[i];
            
            if (mol0.selectedAll())
            {
                grp0_coords[i] = mol0.property(coords_prop)
                                     .asA<AtomCoords>()
                                     .toVector();
            
                grp0_chgs[i] = mol0.property(charge_prop)
                                   .asA<AtomCharges>()
                                   .toVector();
                                                  
                grp0_ljs[i] = mol0.property(lj_prop)
                                  .asA<AtomLJs>()
                                  .toVector();
            }
            else
            {
                const AtomSelection selected_atoms = mol0.selection();

                grp0_coords[i] = mol0.property(coords_prop)
                                     .asA<AtomCoords>()
                                     .toVector(selected_atoms);
            
                grp0_chgs[i] = mol0.property(charge_prop)
                                   .asA<AtomCharges>()
                                   .toVector(selected_atoms);
                                                  
                grp0_ljs[i] = mol0.property(lj_prop)
                                  .asA<AtomLJs>()
                                  .toVector(selected_atoms);
            }
        }

        for (int i=0; i<n1; ++i)
        {
            const PartialMolecule &mol1 = v1.constData()[i];
            
            if (mol1.selectedAll())
            {
                grp1_coords[i] = mol1.property(coords_prop)
                                     .asA<AtomCoords>()
                                     .toVector();
            
                grp1_chgs[i] = mol1.property(charge_prop)
                                   .asA<AtomCharges>()
                                   .toVector();
                                                  
                grp1_ljs[i] = mol1.property(lj_prop)
                                  .asA<AtomLJs>()
                                  .toVector();
            }
            else
            {
                const AtomSelection selected_atoms = mol1.selection();

                grp1_coords[i] = mol1.property(coords_prop)
                                     .asA<AtomCoords>()
                                     .toVector(selected_atoms);
            
                grp1_chgs[i] = mol1.property(charge_prop)
                                   .asA<AtomCharges>()
                                   .toVector(selected_atoms);
                                                  
                grp1_ljs[i] = mol1.property(lj_prop)
                                  .asA<AtomLJs>()
                                  .toVector(selected_atoms);
            }
        }
    
        if (this->usesSoftCore())
        {
            //calculate the energy of each view-pair
            for (int i=0; i<n0; ++i)
            {
                const QVector<Vector> &mol0_coords = grp0_coords.at(i);
                const QVector<Charge> &mol0_chgs = grp0_chgs.at(i);
                const QVector<LJParameter> &mol0_ljs = grp0_ljs.at(i);
            
                for (int j=0; j<n1; ++j)
                {
                    const QVector<Vector> mol1_coords = grp1_coords.at(j);
                    const QVector<Charge> mol1_chgs = grp1_chgs.at(j);
                    const QVector<LJParameter> mol1_ljs = grp1_ljs.at(j);
                
                    pair<double,double> cljnrg 
                        = getSoftCLJEnergy(mol0_coords, mol0_chgs, mol0_ljs,
                                           mol1_coords, mol1_chgs, mol1_ljs,
                                           alfa, shift_delta, coulomb_power);
            
                    coul_nrgs(i,j).edit().accumulate(cljnrg.first);
                    lj_nrgs(i,j).edit().accumulate(cljnrg.second);
                }
            }
        }
        else
        {
            //calculate the energy of each view-pair
            for (int i=0; i<n0; ++i)
            {
                const QVector<Vector> &mol0_coords = grp0_coords.at(i);
                const QVector<Charge> &mol0_chgs = grp0_chgs.at(i);
                const QVector<LJParameter> &mol0_ljs = grp0_ljs.at(i);
            
                for (int j=0; j<n1; ++j)
                {
                    const QVector<Vector> mol1_coords = grp1_coords.at(j);
                    const QVector<Charge> mol1_chgs = grp1_chgs.at(j);
                    const QVector<LJParameter> mol1_ljs = grp1_ljs.at(j);
                
                    pair<double,double> cljnrg 
                        = getCLJEnergy(mol0_coords, mol0_chgs, mol0_ljs,
                                       mol1_coords, mol1_chgs, mol1_ljs);
            
                    coul_nrgs(i,j).edit().accumulate(cljnrg.first);
                    lj_nrgs(i,j).edit().accumulate(cljnrg.second);
                }
            }
        }
    }
    catch(...)
    {
        EnergyMonitor::operator=(old_state);
        throw;
    }
}
