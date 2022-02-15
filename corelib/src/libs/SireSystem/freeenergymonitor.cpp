/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2013  Christopher Woods
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

#include "SireSystem/freeenergymonitor.h"
#include "SireSystem/system.h"

#include "SireMol/partialmolecule.h"
#include "SireMol/atomcoords.h"
#include "SireMol/atomcharges.h"
#include "SireMM/atomljs.h"
#include "SireMM/ljpair.h"
#include "SireMol/mgname.h"

#include "SireMM/cljatoms.h"
#include "SireMM/cljshiftfunction.h"

#include "SireUnits/units.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QTime>

using namespace SireSystem;
using namespace SireMol;
using namespace SireMaths;
using namespace SireMM;
using namespace SireUnits;
using namespace SireUnits::Dimension;
using namespace SireBase;
using namespace SireStream;

using std::pair;

/////////
///////// Implementation of AssignerGroup
/////////

static const RegisterMetaType<AssignerGroup> r_assignergroup(NO_ROOT);

QDataStream &operator<<(QDataStream &ds, const AssignerGroup &group)
{
    writeHeader(ds, r_assignergroup, 1);

    SharedDataStream sds(ds);

    sds << group.molgroup << group.assgnr;

    return ds;
}

QDataStream &operator>>(QDataStream &ds, AssignerGroup &group)
{
    VersionID v = readHeader(ds, r_assignergroup);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> group.molgroup >> group.assgnr;
    }
    else
        throw version_error(v, "1", r_assignergroup, CODELOC);

    return ds;
}

/** Constructor */
AssignerGroup::AssignerGroup()
{}

/** Construct as a holder for a MoleculeGroup */
AssignerGroup::AssignerGroup(const MoleculeGroup &group)
              : molgroup(group)
{}

/** Construct as a holder for an IDAssigner */
AssignerGroup::AssignerGroup(const IDAssigner &assigner)
              : assgnr(assigner)
{}

/** Copy constructor */
AssignerGroup::AssignerGroup(const AssignerGroup &other)
              : molgroup(other.molgroup), assgnr(other.assgnr)
{}

/** Destructor */
AssignerGroup::~AssignerGroup()
{}

/** Copy assignment operator */
AssignerGroup& AssignerGroup::operator=(const AssignerGroup &other)
{
    molgroup = other.molgroup;
    assgnr = other.assgnr;
    return *this;
}

/** Comparison operator */
bool AssignerGroup::operator==(const AssignerGroup &other) const
{
    return molgroup == other.molgroup and assgnr == other.assgnr;
}

/** Comparison operator */
bool AssignerGroup::operator!=(const AssignerGroup &other) const
{
    return not AssignerGroup::operator==(other);
}

const char* AssignerGroup::typeName()
{
    return QMetaType::typeName( qMetaTypeId<AssignerGroup>() );
}

const char* AssignerGroup::what() const
{
    return AssignerGroup::typeName();
}

/** Return whether or not this is empty */
bool AssignerGroup::isEmpty() const
{
    if (molgroup.read().nViews() == 0 and assgnr.isNull())
        return true;

    else if (assgnr.isNull())
        return molgroup->isEmpty();

    else
    {
        return assgnr->asA<IDAssigner>().nPoints() == 0;
    }
}

/** Return whether or not this is a holding a MoleculeGroup */
bool AssignerGroup::isMoleculeGroup() const
{
    if (not assgnr.isNull())
        return false;
    else
        return molgroup.read().nViews() > 0;
}

/** Return whether or not this is holding an IDAssigner */
bool AssignerGroup::isAssigner() const
{
    return not assgnr.isNull();
}

/** Return the molecule group */
const MoleculeGroup& AssignerGroup::group() const
{
    if (molgroup.isNull())
        throw SireError::invalid_state( QObject::tr(
                "Cannot return a MoleculeGroup from an AssignerGroup that does not "
                "hold a MoleculeGroup object."), CODELOC );

    return molgroup.read();
}

/** Return the IDAssigner */
const IDAssigner& AssignerGroup::assigner() const
{
    if (assgnr.isNull())
        throw SireError::invalid_state( QObject::tr(
                "Cannot return an IDAssigner from an AssignerGroup that does not "
                "hold an IDAssigner object."), CODELOC );

    return assgnr.read().asA<IDAssigner>();
}

/** Return the molecule views contained in this group */
QVector<PartialMolecule> AssignerGroup::views() const
{
    if (this->isEmpty())
    {
        return QVector<PartialMolecule>();
    }
    else if (this->isMoleculeGroup())
    {
        const MoleculeGroup &group = this->group();

        int n = group.nViews();

        QVector<PartialMolecule> views(n);

        for (int i=0; i<n; ++i)
        {
            views[i] = group.viewAt(i);
        }

        return views;
    }
    else
    {
        return this->assigner().identifiedMolecules();
    }
}

/** Update the contained group or assigner to match the version
    in the passed system */
void AssignerGroup::update(const System &system)
{
    if (not assgnr.isNull())
    {
        assgnr.edit().asA<IDAssigner>().update(system);
    }
    else if (not molgroup.isNull())
    {
        if (system.contains(molgroup.read().number()))
            molgroup = system[molgroup.read().number()];
    }
}

/** Return whether or not this group is compatible with 'other'.
    Compatible means is the same type, refers to the same MoleculeGroup etc. */
bool AssignerGroup::isCompatible(const AssignerGroup &other) const
{
    if (this == &other)
        return true;

    if (assgnr.isNull())
    {
        if (not other.assgnr.isNull())
            return false;

        if (molgroup.isNull())
            return other.molgroup.isNull();

        //check that the molecule group number is the same and the number
        //of views is the same
        return this->group().number() == other.group().number() and
               this->group().nViews() == other.group().nViews();
    }
    else
    {
        if (other.assgnr.isNull())
            return false;

        // we cannot do much of a comparison of assigners?
        // Just check that the number of points is the same
        return assigner().nPoints() == other.assigner().nPoints();
    }
}

/////////
///////// Implementation of FreeEnergyMonitor
/////////

static const RegisterMetaType<FreeEnergyMonitor> r_nrgmonitor;

QDataStream &operator<<(QDataStream &ds,
                                          const FreeEnergyMonitor &nrgmonitor)
{
    writeHeader(ds, r_nrgmonitor, 1);

    SharedDataStream sds(ds);

    sds << nrgmonitor.refgroup
        << nrgmonitor.group_a << nrgmonitor.group_b
        << nrgmonitor.total_nrgs
        << nrgmonitor.coul_nrgs << nrgmonitor.lj_nrgs
        << nrgmonitor.nrg_template
        << nrgmonitor.lambda_symbol
        << nrgmonitor.shift_delta
        << nrgmonitor.lamval
        << nrgmonitor.delta_lambda
        << nrgmonitor.coulomb_power
        << static_cast<const SystemMonitor&>(nrgmonitor);

    return ds;
}

QDataStream &operator>>(QDataStream &ds,
                                          FreeEnergyMonitor &nrgmonitor)
{
    VersionID v = readHeader(ds, r_nrgmonitor);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> nrgmonitor.refgroup
            >> nrgmonitor.group_a >> nrgmonitor.group_b
            >> nrgmonitor.total_nrgs
            >> nrgmonitor.coul_nrgs >> nrgmonitor.lj_nrgs
            >> nrgmonitor.nrg_template
            >> nrgmonitor.lambda_symbol
            >> nrgmonitor.shift_delta
            >> nrgmonitor.lamval
            >> nrgmonitor.delta_lambda
            >> nrgmonitor.coulomb_power
            >> static_cast<SystemMonitor&>(nrgmonitor);
    }
    else
        throw version_error(v, "1", r_nrgmonitor, CODELOC);

    return ds;
}

/** Null constructor. By default, we don't collect a histogram of each of
    the components energies as this is too memory hungry */
FreeEnergyMonitor::FreeEnergyMonitor()
                  : ConcreteProperty<FreeEnergyMonitor,SystemMonitor>(),
                    nrg_template( FreeEnergyAverage(25*celsius, MolarEnergy(0)) ),
                    lambda_symbol("lambda"), shift_delta(0),
                    lamval(0), delta_lambda(0.001), coulomb_power(0)
{}

/** Construct to monitor the free energy difference of the reference group
    interacting with group A as it is perturbed into group B. By default,
    we don't collect a histogram of each of the components energies as this
    is too memory hungry */
FreeEnergyMonitor::FreeEnergyMonitor(const AssignerGroup &ref,
                                     const AssignerGroup &ga,
                                     const AssignerGroup &gb)
              : ConcreteProperty<FreeEnergyMonitor,SystemMonitor>(),
                refgroup(ref), group_a(ga), group_b(gb),
                nrg_template( FreeEnergyAverage(25*celsius, MolarEnergy(0)) ),
                lambda_symbol("lambda"), shift_delta(0), lamval(0), delta_lambda(0.001),
                coulomb_power(0)
{}

/** Copy constructor */
FreeEnergyMonitor::FreeEnergyMonitor(const FreeEnergyMonitor &other)
                  : ConcreteProperty<FreeEnergyMonitor,SystemMonitor>(other),
                    refgroup(other.refgroup), group_a(other.group_a), group_b(other.group_b),
                    total_nrgs(other.total_nrgs),
                    coul_nrgs(other.coul_nrgs), lj_nrgs(other.lj_nrgs),
                    nrg_template(other.nrg_template),
                    lambda_symbol(other.lambda_symbol),
                    shift_delta(other.shift_delta),
                    lamval(other.lamval),
                    delta_lambda(other.delta_lambda),
                    coulomb_power(other.coulomb_power)
{}

/** Destructor */
FreeEnergyMonitor::~FreeEnergyMonitor()
{}

/** Copy assignment operator */
FreeEnergyMonitor& FreeEnergyMonitor::operator=(const FreeEnergyMonitor &other)
{
    if (this != &other)
    {
        refgroup = other.refgroup;
        group_a = other.group_a;
        group_b = other.group_b;
        total_nrgs = other.total_nrgs;
        coul_nrgs = other.coul_nrgs;
        lj_nrgs = other.lj_nrgs;
        nrg_template = other.nrg_template;
        lambda_symbol = other.lambda_symbol;
        shift_delta = other.shift_delta;
        lamval = other.lamval;
        delta_lambda = other.delta_lambda;
        coulomb_power = other.coulomb_power;
        SystemMonitor::operator=(other);
    }

    return *this;
}

/** Comparison operator */
bool FreeEnergyMonitor::operator==(const FreeEnergyMonitor &other) const
{
    return this == &other or
           (refgroup == other.refgroup and
            group_a == other.group_a and group_b == other.group_b and
            total_nrgs == other.total_nrgs and
            coul_nrgs == other.coul_nrgs and
            lj_nrgs == other.lj_nrgs and
            nrg_template == other.nrg_template and
            lambda_symbol == other.lambda_symbol and
            shift_delta == other.shift_delta and
            lamval == other.lamval and
            delta_lambda == other.delta_lambda and
            coulomb_power == other.coulomb_power and
            SystemMonitor::operator==(other));
}

/** Comparison operator */
bool FreeEnergyMonitor::operator!=(const FreeEnergyMonitor &other) const
{
    return not FreeEnergyMonitor::operator==(other);
}

/** Return the typename of the class */
const char* FreeEnergyMonitor::typeName()
{
    return QMetaType::typeName( qMetaTypeId<FreeEnergyMonitor>() );
}

/** Return the array of all accumulated total free energies */
QVector<FreeEnergyAverage> FreeEnergyMonitor::freeEnergies() const
{
    return total_nrgs;
}

/** Return the array of all accumulated coulomb free energies */
QVector<FreeEnergyAverage> FreeEnergyMonitor::coulombFreeEnergies() const
{
    return coul_nrgs;
}

/** Return the array of all accumulated LJ free energies */
QVector<FreeEnergyAverage> FreeEnergyMonitor::ljFreeEnergies() const
{
    return lj_nrgs;
}

/** Return the array of the reference group molecule views in the same order as they
    appear in the arrays of free energies */
QVector<PartialMolecule> FreeEnergyMonitor::referenceViews() const
{
    return refgroup.views();
}

/** Return the reference group */
const AssignerGroup& FreeEnergyMonitor::referenceGroup() const
{
    return refgroup;
}

/** Return group A (the group interacting with the reference group at lambda = 0) */
const AssignerGroup& FreeEnergyMonitor::groupA() const
{
    return group_a;
}

/** Return group B (the group interacting with the reference group at lambda = 1) */
const AssignerGroup& FreeEnergyMonitor::groupB() const
{
    return group_b;
}

/** Set the temperature at which the free energies are accumulated */
void FreeEnergyMonitor::setTemperature(const SireUnits::Dimension::Temperature &temperature)
{
    nrg_template = FreeEnergyAverage(temperature,
                                     MolarEnergy(nrg_template.histogram().binWidth()));
}

/** Set the bin width for the histogram of recorded free energies.
    By default, we don't collect a histogram of each of the components energies
    as this is too memory hungry. */
void FreeEnergyMonitor::setBinWidth(const SireUnits::Dimension::MolarEnergy &binwidth)
{
    nrg_template = FreeEnergyAverage(nrg_template.temperature(), binwidth);
}

/** Return the temperature at which the free energies are calculated */
SireUnits::Dimension::Temperature FreeEnergyMonitor::temperature() const
{
    return nrg_template.temperature();
}

/** Set the value of delta lambda to use to calculate the free energy gradients */
void FreeEnergyMonitor::setDeltaLambda(double delta_lam)
{
    if (delta_lam > 0 and delta_lam < 1)
        delta_lambda = delta_lam;
}

/** Return the value of delta lambda used to calculate the free energy gradients */
double FreeEnergyMonitor::deltaLambda() const
{
    return delta_lambda;
}

/** Return the width of the histogram bins used to record the free energies */
SireUnits::Dimension::MolarEnergy FreeEnergyMonitor::binWidth() const
{
    return MolarEnergy(nrg_template.histogram().binWidth());
}

/** Clear all statistics */
void FreeEnergyMonitor::clearStatistics()
{
    total_nrgs = QVector<FreeEnergyAverage>();
    coul_nrgs = QVector<FreeEnergyAverage>();
    lj_nrgs = QVector<FreeEnergyAverage>();
}

SIRE_ALWAYS_INLINE pair<double,double> getCLJEnergy(
                const QVector<Vector> &coords, const QVector<Charge> &chgs,
                const QVector<LJParameter> &ljs,
                const QVector<Vector> &coords_a, const QVector<Charge> &chgs_a,
                const QVector<LJParameter> &ljs_a,
                const QVector<Vector> &coords_b, const QVector<Charge> &chgs_b,
                const QVector<LJParameter> &ljs_b,
                double lamval, double delta_lambda)
{
    const int nats = coords.count();
    const int nats_a = coords_a.count();
    const int nats_b = coords_b.count();

    double cnrg = 0;
    double ljnrg = 0;
    double cnrg_f = 0;
    double ljnrg_f = 0;

    if (lamval < 0)
        lamval = 0;

    if (lamval+delta_lambda > 1)
        lamval = 1 - delta_lambda;

    bool arithmetic_combining_rules = true;

    // total_nrg(lam) = (1-lam) * group:group_A + lam * group:group_B
    // delta_nrg = total_nrg(lam+delta_lambda) - total_nrg(lam)

    for (int i=0; i<nats; ++i)
    {
        const Vector &coord0 = coords.at(i);
        const Charge &chg0 = chgs.at(i);
        const LJParameter &lj0 = ljs.at(i);

        double icnrg = 0;
        double iljnrg = 0;

        for (int j=0; j<nats_a; ++j)
        {
            const Vector &coord1 = coords_a.at(j);
            const Charge &chg1 = chgs_a.at(j);
            const LJParameter &lj1 = ljs_a.at(j);

            LJPair ljpair;

            if (arithmetic_combining_rules)
                ljpair = LJPair::arithmetic(lj0, lj1);
            else
                ljpair = LJPair::geometric(lj0, lj1);

            double one_over_r2 = Vector::invDistance2(coord0, coord1);
            double one_over_r = std::sqrt(one_over_r2);
            double one_over_r6 = one_over_r2 * one_over_r2 * one_over_r2;
            double one_over_r12 = one_over_r6 * one_over_r6;

            iljnrg += ljpair.A() * one_over_r12 - ljpair.B() * one_over_r6;

            icnrg += chg0.value() * chg1.value() * one_over_r
                                  * one_over_four_pi_eps0;
        }

        cnrg += (1-lamval) * icnrg;
        cnrg_f += (1-lamval-delta_lambda) * icnrg;
        ljnrg += (1-lamval) * iljnrg;
        ljnrg_f += (1-lamval-delta_lambda) * iljnrg;

        icnrg = 0;
        iljnrg = 0;

        for (int j=0; j<nats_b; ++j)
        {
            const Vector &coord1 = coords_b.at(j);
            const Charge &chg1 = chgs_b.at(j);
            const LJParameter &lj1 = ljs_b.at(j);

            LJPair ljpair;

            if (arithmetic_combining_rules)
                ljpair = LJPair::arithmetic(lj0, lj1);
            else
                ljpair = LJPair::geometric(lj0, lj1);

            double one_over_r2 = Vector::invDistance2(coord0, coord1);
            double one_over_r = std::sqrt(one_over_r2);
            double one_over_r6 = one_over_r2 * one_over_r2 * one_over_r2;
            double one_over_r12 = one_over_r6 * one_over_r6;

            iljnrg += ljpair.A() * one_over_r12 - ljpair.B() * one_over_r6;

            icnrg += chg0.value() * chg1.value() * one_over_r
                                  * one_over_four_pi_eps0;
        }

        cnrg += (lamval) * icnrg;
        cnrg_f += (lamval+delta_lambda) * icnrg;
        ljnrg += (lamval) * iljnrg;
        ljnrg_f += (lamval+delta_lambda) * iljnrg;
    }

    return pair<double,double>(cnrg_f-cnrg,ljnrg_f-ljnrg);
}

SIRE_ALWAYS_INLINE pair<double,double> getSoftCLJEnergy(
                const CLJAtoms &refatoms,
                const CLJAtoms &atoms_a, const CLJAtoms &atoms_b,
                double lamval, double delta_lambda, double shift_delta, int coulomb_power)
{
    if (lamval < 0)
        lamval = 0;

    if (lamval+delta_lambda > 1)
        lamval = 1 - delta_lambda;

    bool arithmetic_combining_rules = true;

    const double alpha_a = lamval;
    const double alpha_a_f = lamval + delta_lambda;
    double one_minus_alpha_a_to_n = 1;
    double one_minus_alpha_a_to_n_f = 1;
    const double delta_a = shift_delta * alpha_a;
    const double delta_a_f = shift_delta * alpha_a_f;

    const double alpha_b = 1 - lamval;
    const double alpha_b_f = 1 - lamval - delta_lambda;
    double one_minus_alpha_b_to_n = 1;
    double one_minus_alpha_b_to_n_f = 1;
    const double delta_b = shift_delta * alpha_b;
    const double delta_b_f = shift_delta * alpha_b_f;

    if (coulomb_power != 0)
    {
        one_minus_alpha_a_to_n = SireMaths::pow(1 - alpha_a, coulomb_power);
        one_minus_alpha_a_to_n_f = SireMaths::pow(1 - alpha_a_f, coulomb_power);
        one_minus_alpha_b_to_n = SireMaths::pow(1 - alpha_b, coulomb_power);
        one_minus_alpha_b_to_n_f = SireMaths::pow(1 - alpha_b_f, coulomb_power);
    }

    CLJSoftShiftFunction cljfunc_a, cljfunc_a_f;

    cljfunc_a.setShiftDelta(delta_a);
    cljfunc_a_f.setShiftDelta(delta_a_f);
    cljfunc_a.setCoulombPower(coulomb_power);
    cljfunc_a_f.setCoulombPower(coulomb_power);
    cljfunc_a.setArithmeticCombiningRules(arithmetic_combining_rules);
    cljfunc_a_f.setArithmeticCombiningRules(arithmetic_combining_rules);

    CLJSoftShiftFunction cljfunc_b, cljfunc_b_f;

    cljfunc_b.setShiftDelta(delta_b);
    cljfunc_b_f.setShiftDelta(delta_b_f);
    cljfunc_b.setCoulombPower(coulomb_power);
    cljfunc_b_f.setCoulombPower(coulomb_power);
    cljfunc_b.setArithmeticCombiningRules(arithmetic_combining_rules);
    cljfunc_b_f.setArithmeticCombiningRules(arithmetic_combining_rules);

    double cnrg_a(0), cnrg_a_f(0), cnrg_b(0), cnrg_b_f(0);
    double ljnrg_a(0), ljnrg_a_f(0), ljnrg_b(0), ljnrg_b_f(0);

    cljfunc_a.total(refatoms, atoms_a, cnrg_a, ljnrg_a);
    cljfunc_a_f.total(refatoms, atoms_a, cnrg_a_f, ljnrg_a_f);
    cljfunc_b.total(refatoms, atoms_b, cnrg_b, ljnrg_b);
    cljfunc_b_f.total(refatoms, atoms_b, cnrg_b_f, ljnrg_b_f);

    cnrg_a *= (1 - lamval) * one_minus_alpha_a_to_n;
    ljnrg_a *= (1 - lamval) * one_minus_alpha_a_to_n;

    cnrg_a_f *= (1 - lamval - delta_lambda) * one_minus_alpha_a_to_n_f;
    ljnrg_a_f *= (1 - lamval - delta_lambda) * one_minus_alpha_a_to_n_f;

    cnrg_b *= (lamval) * one_minus_alpha_b_to_n;
    ljnrg_b *= (lamval) * one_minus_alpha_b_to_n;

    cnrg_b_f *= (lamval+delta_lambda) * one_minus_alpha_b_to_n_f;
    ljnrg_b_f *= (lamval+delta_lambda) * one_minus_alpha_b_to_n_f;

    return pair<double,double>(cnrg_a_f+cnrg_b_f-cnrg_a-cnrg_b,
                               ljnrg_a_f+ljnrg_b_f-ljnrg_a-ljnrg_b);
}

SIRE_ALWAYS_INLINE pair<double,double> getSoftCLJEnergy(
                const QVector<Vector> &coords, const QVector<Charge> &chgs,
                const QVector<LJParameter> &ljs,
                const QVector<Vector> &coords_a, const QVector<Charge> &chgs_a,
                const QVector<LJParameter> &ljs_a,
                const QVector<Vector> &coords_b, const QVector<Charge> &chgs_b,
                const QVector<LJParameter> &ljs_b,
                double lamval, double delta_lambda, double shift_delta, int coulomb_power)
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

    const int nats = coords.count();
    const int nats_a = coords_a.count();
    const int nats_b = coords_b.count();

    double cnrg = 0;
    double ljnrg = 0;
    double cnrg_f = 0;
    double ljnrg_f = 0;

    if (lamval < 0)
        lamval = 0;

    if (lamval+delta_lambda > 1)
        lamval = 1 - delta_lambda;

    bool arithmetic_combining_rules = true;

    const double alpha_a = lamval;
    const double alpha_a_f = lamval + delta_lambda;
    double one_minus_alpha_a_to_n = 1;
    double one_minus_alpha_a_to_n_f = 1;
    const double delta_a = shift_delta * alpha_a;
    const double delta_a_f = shift_delta * alpha_a_f;

    const double alpha_b = 1 - lamval;
    const double alpha_b_f = 1 - lamval - delta_lambda;
    double one_minus_alpha_b_to_n = 1;
    double one_minus_alpha_b_to_n_f = 1;
    const double delta_b = shift_delta * alpha_b;
    const double delta_b_f = shift_delta * alpha_b_f;

    if (coulomb_power != 0)
    {
        one_minus_alpha_a_to_n = SireMaths::pow(1 - alpha_a, coulomb_power);
        one_minus_alpha_a_to_n_f = SireMaths::pow(1 - alpha_a_f, coulomb_power);
        one_minus_alpha_b_to_n = SireMaths::pow(1 - alpha_b, coulomb_power);
        one_minus_alpha_b_to_n_f = SireMaths::pow(1 - alpha_b_f, coulomb_power);
    }

    // total_nrg(lam) = (1-lam) * group:group_A + lam * group:group_B
    // delta_nrg = total_nrg(lam+delta_lambda) - total_nrg(lam)

    for (int i=0; i<nats; ++i)
    {
        const Vector &coord0 = coords.at(i);
        const Charge &chg0 = chgs.at(i);
        const LJParameter &lj0 = ljs.at(i);

        double icnrg = 0;
        double iljnrg = 0;
        double icnrg_f = 0;
        double iljnrg_f = 0;

        for (int j=0; j<nats_a; ++j)
        {
            const Vector &coord1 = coords_a.at(j);
            const Charge &chg1 = chgs_a.at(j);
            const LJParameter &lj1 = ljs_a.at(j);

            LJPair ljpair;

            if (arithmetic_combining_rules)
                ljpair = LJPair::arithmetic(lj0, lj1);
            else
                ljpair = LJPair::geometric(lj0, lj1);

            double r2 = Vector::distance2(coord0, coord1);

            const double shift = ljpair.sigma() * delta_a;
            const double shift_f = ljpair.sigma() * delta_a_f;
            double lj_denom = r2 + shift;
            lj_denom = lj_denom * lj_denom * lj_denom;

            double lj_denom_f = r2 + shift_f;
            lj_denom_f = lj_denom_f * lj_denom_f * lj_denom_f;

            const double sig2 = ljpair.sigma() * ljpair.sigma();
            const double sig6 = sig2 * sig2 * sig2;

            const double sig6_over_denom = sig6 / lj_denom;
            const double sig6_over_denom_f = sig6 / lj_denom_f;
            const double sig12_over_denom2 = sig6_over_denom *
                                             sig6_over_denom;
            const double sig12_over_denom2_f = sig6_over_denom_f *
                                               sig6_over_denom_f;

            iljnrg += 4 * ljpair.epsilon() * (sig12_over_denom2 - sig6_over_denom);
            iljnrg_f += 4 * ljpair.epsilon() * (sig12_over_denom2_f - sig6_over_denom_f);

            icnrg += one_minus_alpha_a_to_n *
                        chg0.value() * chg1.value() * one_over_four_pi_eps0 /
                           std::sqrt(alpha_a + r2);

            icnrg_f += one_minus_alpha_a_to_n_f *
                        chg0.value() * chg1.value() * one_over_four_pi_eps0 /
                           std::sqrt(alpha_a_f + r2);
        }

        cnrg += (1-lamval) * icnrg;
        cnrg_f += (1-lamval-delta_lambda) * icnrg_f;
        ljnrg += (1-lamval) * iljnrg;
        ljnrg_f += (1-lamval-delta_lambda) * iljnrg_f;

        icnrg = 0;
        iljnrg = 0;
        icnrg_f = 0;
        iljnrg_f = 0;

        for (int j=0; j<nats_b; ++j)
        {
            const Vector &coord1 = coords_b.at(j);
            const Charge &chg1 = chgs_b.at(j);
            const LJParameter &lj1 = ljs_b.at(j);

            LJPair ljpair;

            if (arithmetic_combining_rules)
                ljpair = LJPair::arithmetic(lj0, lj1);
            else
                ljpair = LJPair::geometric(lj0, lj1);

            double r2 = Vector::distance2(coord0, coord1);

            const double shift = ljpair.sigma() * delta_b;
            const double shift_f = ljpair.sigma() * delta_b_f;
            double lj_denom = r2 + shift;
            lj_denom = lj_denom * lj_denom * lj_denom;

            double lj_denom_f = r2 + shift_f;
            lj_denom_f = lj_denom_f * lj_denom_f * lj_denom_f;

            const double sig2 = ljpair.sigma() * ljpair.sigma();
            const double sig6 = sig2 * sig2 * sig2;

            const double sig6_over_denom = sig6 / lj_denom;
            const double sig6_over_denom_f = sig6 / lj_denom_f;
            const double sig12_over_denom2 = sig6_over_denom *
                                             sig6_over_denom;
            const double sig12_over_denom2_f = sig6_over_denom_f *
                                               sig6_over_denom_f;

            iljnrg += 4 * ljpair.epsilon() * (sig12_over_denom2 - sig6_over_denom);
            iljnrg_f += 4 * ljpair.epsilon() * (sig12_over_denom2_f - sig6_over_denom_f);

            icnrg += one_minus_alpha_b_to_n *
                        chg0.value() * chg1.value() * one_over_four_pi_eps0 /
                           std::sqrt(alpha_b + r2);

            icnrg_f += one_minus_alpha_b_to_n_f *
                        chg0.value() * chg1.value() * one_over_four_pi_eps0 /
                           std::sqrt(alpha_b_f + r2);
        }

        cnrg += (lamval) * icnrg;
        cnrg_f += (lamval+delta_lambda) * icnrg_f;
        ljnrg += (lamval) * iljnrg;
        ljnrg_f += (lamval+delta_lambda) * iljnrg_f;
    }

    return pair<double,double>(cnrg_f-cnrg,ljnrg_f-ljnrg);
}

/** Return whether or not this monitor uses a soft-core potential to
    calculate the CLJ energy between the molecules in views0() and the
    molecules in views1() */
bool FreeEnergyMonitor::usesSoftCore() const
{
    return (shift_delta != 0 or coulomb_power != 0);
}

/** Set the system component symbol used to get the value of lambda */
void FreeEnergyMonitor::setLambdaComponent(const Symbol &component)
{
    lambda_symbol = component;
}

/** Set the shift delta parameter used by the soft-core potential */
void FreeEnergyMonitor::setShiftDelta(double delta)
{
    shift_delta = delta;
}

/** Set the coulomb power parameter used by the soft-core potential */
void FreeEnergyMonitor::setCoulombPower(int power)
{
    if (power >= 0)
        coulomb_power = power;
}

/** Return the shift delta parameter if a soft-core potential is used.
    This returns 0 if a LJ shifting term is not used */
double FreeEnergyMonitor::shiftDelta() const
{
    return shift_delta;
}

/** Return the coulomb power, if extra coulomb softening is used.
    This returns 0 if coulomb softening is not used */
int FreeEnergyMonitor::coulombPower() const
{
    return coulomb_power;
}

/** Return the lambda value at which the free energy components were collected */
double FreeEnergyMonitor::lambdaValue() const
{
    return lamval;
}

/** Return the symbol used to find the value of lambda from the system */
Symbol FreeEnergyMonitor::lambdaComponent() const
{
    return lambda_symbol;
}

/** Return whether this is empty (has no group data) */
bool FreeEnergyMonitor::isEmpty() const
{
    return refgroup.isEmpty();
}

/** Conserve memory by copying the molecule data etc. from 'other' into this monitor */
void FreeEnergyMonitor::conserveMemory(const FreeEnergyMonitor &other)
{
    if (this->isCompatibleExceptLambda(other))
    {
        refgroup = other.refgroup;
        group_a = other.group_a;
        group_b = other.group_b;
    }
}

/** Return whether or not this monitor is compatible with 'other'
    (have the same groups, soft-core parameters, delta lambda, temperature etc.) */
bool FreeEnergyMonitor::isCompatible(const FreeEnergyMonitor &other) const
{
    return this->lambdaValue() == other.lambdaValue() and
           this->lambdaComponent() == other.lambdaComponent() and
           this->usesSoftCore() == other.usesSoftCore() and
           this->shiftDelta() == other.shiftDelta() and
           this->coulombPower() == other.coulombPower() and
           this->temperature() == other.temperature() and
           this->binWidth() == other.binWidth() and
           refgroup.isCompatible(other.refgroup) and
           group_a.isCompatible(other.group_a) and
           group_b.isCompatible(other.group_b);
}

/** Return whether or not this monitor is compatible with other, ignoring that
    the monitors have different lambda values. This will let you know if it is sensible
    to construct PMFs from a combination of these monitors */
bool FreeEnergyMonitor::isCompatibleExceptLambda(const FreeEnergyMonitor &other) const
{
    return this->lambdaComponent() == other.lambdaComponent() and
           this->usesSoftCore() == other.usesSoftCore() and
           this->shiftDelta() == other.shiftDelta() and
           this->coulombPower() == other.coulombPower() and
           this->temperature() == other.temperature() and
           this->binWidth() == other.binWidth() and
           refgroup.isCompatible(other.refgroup) and
           group_a.isCompatible(other.group_a) and
           group_b.isCompatible(other.group_b);
}

/** Self-addition operator - you can only add two monitors together if they
    have the same groups, soft-core parameters, delta lambda and temperature etc.

    \throw SireError::incompatible_error
*/
FreeEnergyMonitor& FreeEnergyMonitor::operator+=(const FreeEnergyMonitor &other)
{
    if (this == &other)
    {
        this->operator+=( FreeEnergyMonitor(other) );
        return *this;
    }

    if (not this->isCompatible(other))
    {
        throw SireError::incompatible_error( QObject::tr(
                "Cannot add together two FreeEnergyMonitors as they are in some "
                "way incompatible."), CODELOC );
    }

    if (total_nrgs.isEmpty())
    {
        total_nrgs = other.total_nrgs;
        coul_nrgs = other.coul_nrgs;
        lj_nrgs = other.lj_nrgs;
    }
    else if (not other.total_nrgs.isEmpty())
    {
        if (total_nrgs.count() != other.total_nrgs.count())
            throw SireError::program_bug( QObject::tr(
                    "It should not be possible for two FreeEnergyMonitors to be compatible "
                    "but have different numbers of free energies (%1 vs. %2)")
                        .arg(total_nrgs.count())
                        .arg(other.total_nrgs.count()), CODELOC );

        for (int i=0; i<total_nrgs.count(); ++i)
        {
            total_nrgs[i] += other.total_nrgs[i];
        }

        for (int i=0; i<coul_nrgs.count(); ++i)
        {
            coul_nrgs[i] += other.coul_nrgs[i];
        }

        for (int i=0; i<lj_nrgs.count(); ++i)
        {
            lj_nrgs[i] += other.lj_nrgs[i];
        }
    }

    return *this;
}

/** Addition operator - you can only add two monitors together if they
    have the same groups, soft-core parameters, delta lambda, temperature etc.

    \throw SireError:incompatible_error
*/
FreeEnergyMonitor FreeEnergyMonitor::operator+(const FreeEnergyMonitor &other) const
{
    FreeEnergyMonitor ret(*this);
    ret += other;
    return ret;
}

/** Merge a whole set of free energy monitors together. Note that you can
    only merge them if they have the same groups, soft-core parameters, delta lambda,
    temperature etc.

    \throw SireError::incompatible_error
*/
FreeEnergyMonitor FreeEnergyMonitor::merge(const QList<FreeEnergyMonitor> &monitors)
{
    if (monitors.isEmpty())
        return FreeEnergyMonitor();

    FreeEnergyMonitor ret = monitors.at(0);

    for (int i=1; i<monitors.count(); ++i)
    {
        ret += monitors.at(i);
    }

    return ret;
}

/** Return the number of samples used to form all of the free energy averages
    in this monitor */
int FreeEnergyMonitor::nSamples() const
{
    int nsamples = 0;

    for (int i=0; i<total_nrgs.count(); ++i)
    {
        nsamples += total_nrgs.at(i).nSamples();
    }

    return nsamples;
}

/** Accumulate energies from the passed system */
void FreeEnergyMonitor::monitor(System &system)
{
    FreeEnergyMonitor old_state(*this);

    try
    {
        refgroup.update(system);
        group_a.update(system);
        group_b.update(system);

        //see if we need to update the lambda value
        if (not lambda_symbol.isNull())
        {
            double new_lamval = system.componentValue(lambda_symbol);

            if (lamval != new_lamval and not total_nrgs.isEmpty())
            {
                qDebug() << "WARNING: Throwing away component statistics (1).";
                this->clearStatistics();
            }

            lamval = new_lamval;
        }

        QVector<PartialMolecule> refviews = refgroup.views();
        QVector<PartialMolecule> views_a = group_a.views();
        QVector<PartialMolecule> views_b = group_b.views();

        if (group_a.isAssigner() or group_b.isAssigner())
            throw SireError::unsupported( QObject::tr(
                    "This is not supported"), CODELOC );

        QVector<CLJAtoms> refcljatoms;
        CLJAtoms cljatoms_a(group_a.group());
        CLJAtoms cljatoms_b(group_b.group());

        for (const PartialMolecule &view : refviews)
        {
            refcljatoms.append( CLJAtoms(view) );
        }

        // extract the charge, LJ and coordinates of all of the views
        const PropertyName charge_prop("charge");
        const PropertyName lj_prop("LJ");
        const PropertyName coords_prop("coordinates");

        const int nref = refviews.count();

        if (nref == 0)
        {
            this->clearStatistics();
            return;
        }
        else if (total_nrgs.isEmpty())
        {
            //this is the first time we are calculating energies
            total_nrgs = QVector<FreeEnergyAverage>(nref, nrg_template);
            coul_nrgs = total_nrgs;
            lj_nrgs = total_nrgs;
        }
        else if (total_nrgs.count() != nref)
        {
            qDebug() << "WARNING: Throwing away old component free energies! (2)";
            total_nrgs = QVector<FreeEnergyAverage>(nref, nrg_template);
            coul_nrgs = total_nrgs;
            lj_nrgs = total_nrgs;
        }

        //now get all of the coordinates and charges / LJ parameters of the groups
        QVector< QVector<Vector> > ref_coords(nref);
        QVector< QVector<Charge> > ref_chgs(nref);
        QVector< QVector<LJParameter> > ref_ljs(nref);

        QVector<Vector> a_coords = cljatoms_a.coordinates();
        QVector<Charge> a_chgs = cljatoms_a.charges();
        QVector<LJParameter> a_ljs = cljatoms_a.ljParameters();

        QVector<Vector> b_coords = cljatoms_b.coordinates();
        QVector<Charge> b_chgs = cljatoms_b.charges();
        QVector<LJParameter> b_ljs = cljatoms_b.ljParameters();

        for (int i=0; i<nref; ++i)
        {
            ref_coords[i] = refcljatoms[i].coordinates();
            ref_chgs[i] = refcljatoms[i].charges();
            ref_ljs[i] = refcljatoms[i].ljParameters();
        }

        if (this->usesSoftCore())
        {
            //calculate the energy difference for each view against group_a and group_b
            for (int i=0; i<nref; ++i)
            {
                const QVector<Vector> &mol_coords = ref_coords.at(i);
                const QVector<Charge> &mol_chgs = ref_chgs.at(i);
                const QVector<LJParameter> &mol_ljs = ref_ljs.at(i);

                pair<double,double> delta_cljnrg = getSoftCLJEnergy(
                            mol_coords, mol_chgs, mol_ljs,
                            a_coords, a_chgs, a_ljs,
                            b_coords, b_chgs, b_ljs,
                            lamval, delta_lambda,
                            shift_delta, coulomb_power);

                total_nrgs[i].accumulate(delta_cljnrg.first + delta_cljnrg.second);
                coul_nrgs[i].accumulate(delta_cljnrg.first);
                lj_nrgs[i].accumulate(delta_cljnrg.second);
            }
        }
        else
        {
            //calculate the energy difference for each view against group_a and group_b
            for (int i=0; i<nref; ++i)
            {
                const QVector<Vector> &mol_coords = ref_coords.at(i);
                const QVector<Charge> &mol_chgs = ref_chgs.at(i);
                const QVector<LJParameter> &mol_ljs = ref_ljs.at(i);

                pair<double,double> delta_cljnrg = getCLJEnergy(
                            mol_coords, mol_chgs, mol_ljs,
                            a_coords, a_chgs, a_ljs,
                            b_coords, b_chgs, b_ljs,
                            lamval, delta_lambda);

                total_nrgs[i].accumulate(delta_cljnrg.first + delta_cljnrg.second);
                coul_nrgs[i].accumulate(delta_cljnrg.first);
                lj_nrgs[i].accumulate(delta_cljnrg.second);
            }
        }
    }
    catch(...)
    {
        FreeEnergyMonitor::operator=(old_state);
        throw;
    }
}
