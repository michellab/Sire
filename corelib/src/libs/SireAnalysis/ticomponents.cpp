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

#include "ticomponents.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireError/errors.h"

#include "tostring.h"

using namespace SireAnalysis;
using namespace SireSystem;
using namespace SireMol;
using namespace SireBase;
using namespace SireUnits::Dimension;
using namespace SireStream;

//////////
////////// Implementation of ComponentGradients
//////////

static const RegisterMetaType<ComponentGradients> r_grads;

QDataStream &operator<<(QDataStream &ds, const ComponentGradients &grads)
{
    writeHeader(ds, r_grads, 1);

    SharedDataStream sds(ds);

    sds << grads.grads;

    return ds;
}

QDataStream &operator>>(QDataStream &ds, ComponentGradients &grads)
{
    VersionID v = readHeader(ds, r_grads);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> grads.grads;
    }
    else
        throw version_error(v, "1", r_grads, CODELOC);

    return ds;
}

//assert that the contained gradients are sane
void ComponentGradients::checkSane() const
{
    if (grads.isEmpty())
        return;

    const FreeEnergyMonitor &first = grads.constBegin().value();

    for (QMap<double,FreeEnergyMonitor>::const_iterator it = grads.constBegin();
         it != grads.constEnd();
         ++it)
    {
        if (it.key() != it.value().lambdaValue())
            throw SireError::program_bug( QObject::tr(
                    "Disagreement of lambda value: %1 vs. %2.")
                        .arg(it.key()).arg(it.value().lambdaValue()), CODELOC );

        if (not first.isCompatibleExceptLambda(it.value()))
            throw SireError::incompatible_error( QObject::tr(
                    "All of the FreeEnergyMonitors for each lambda value must be compatible. "
                    "The monitors at lambda=%1 and lambda=%2 are not!")
                        .arg(first.lambdaValue()).arg(it.key()), CODELOC );
    }
}

/** Constructor */
ComponentGradients::ComponentGradients() : ConcreteProperty<ComponentGradients,Property>()
{}

/** This function reduces the memory used by this object by ensuring that
    the FreeEnergyMonitor at each lambda value uses the copy of the
    molecules used at the first lambda value */
void ComponentGradients::conserveMemory()
{
    if (grads.count() <= 1)
        return;

    FreeEnergyMonitor first_monitor = grads[ lambdaValues().at(0) ];

    for (QMap<double,FreeEnergyMonitor>::iterator it = grads.begin();
         it != grads.end();
         ++it)
    {
        it.value().conserveMemory(first_monitor);
    }
}

/** This function conserves memory by copying in all of the shared molecule
    data etc. from 'other' into this object */
void ComponentGradients::conserveMemory(const ComponentGradients &other)
{
    if (grads.isEmpty() or other.grads.isEmpty())
        return;

    FreeEnergyMonitor first_monitor = other.grads[ other.lambdaValues().at(0) ];

    for (QMap<double,FreeEnergyMonitor>::iterator it = grads.begin();
         it != grads.end();
         ++it)
    {
        it.value().conserveMemory(first_monitor);
    }
}

/** Construct from the passed map of component monitors */
ComponentGradients::ComponentGradients(const QMap<double,FreeEnergyMonitor> &gradients,
                                       bool conserve_memory)
                   : ConcreteProperty<ComponentGradients,Property>()
{
    //are any empty?
    bool any_empty = false;
    for (QMap<double,FreeEnergyMonitor>::const_iterator it = gradients.constBegin();
         it != gradients.constEnd();
         ++it)
    {
        if (it.value().isEmpty() or it.value().nSamples() == 0)
        {
            any_empty = true;
            break;
        }
    }

    if (any_empty)
    {
        for (QMap<double,FreeEnergyMonitor>::const_iterator it = gradients.constBegin();
             it != gradients.constEnd();
             ++it)
        {
            if (not (it.value().isEmpty() or it.value().nSamples() == 0))
                grads.insert(it.key(), it.value());
        }
    }
    else if (not gradients.isEmpty())
    {
        grads = gradients;
    }

    checkSane();

    if (conserve_memory)
        conserveMemory();
}

/** Construct from the passed list of component monitors */
ComponentGradients::ComponentGradients(const QList<FreeEnergyMonitor> &gradients,
                                       bool conserve_memory)
                   : ConcreteProperty<ComponentGradients,Property>()
{
    foreach (const FreeEnergyMonitor &gradient, gradients)
    {
        if (not (gradient.isEmpty() or gradient.nSamples() == 0))
            grads.insert(gradient.lambdaValue(), gradient);
    }

    checkSane();

    if (conserve_memory)
        conserveMemory();
}

/** Copy constructor */
ComponentGradients::ComponentGradients(const ComponentGradients &other)
                   : ConcreteProperty<ComponentGradients,Property>(),
                     grads(other.grads)
{}

/** Destructor */
ComponentGradients::~ComponentGradients()
{}

/** Copy assignment operator */
ComponentGradients& ComponentGradients::operator=(const ComponentGradients &other)
{
    grads = other.grads;
    return *this;
}

/** Comparison operator */
bool ComponentGradients::operator==(const ComponentGradients &other) const
{
    return grads == other.grads;
}

/** Comparison operator */
bool ComponentGradients::operator!=(const ComponentGradients &other) const
{
    return not ComponentGradients::operator==(other);
}

const char* ComponentGradients::typeName()
{
    return QMetaType::typeName(qMetaTypeId<ComponentGradients>());
}

const char* ComponentGradients::what() const
{
    return ComponentGradients::typeName();
}

QString ComponentGradients::toString() const
{
    return QObject::tr("ComponentGradients( nComponents() == %1, nSamples() == %2, "
                                           "nLambdaValues() == %3 )")
                    .arg(nComponents())
                    .arg(nSamples())
                    .arg(nLambdaValues());
}

/** Return whether or not this set is empty */
bool ComponentGradients::isEmpty() const
{
    return grads.isEmpty();
}

/** Return whether or not this set of gradients is compatible with the ones
    provided in 'other' */
bool ComponentGradients::isCompatible(const ComponentGradients &other) const
{
    if (not (grads.isEmpty() or other.grads.isEmpty()))
    {
        return grads.constBegin().value().isCompatible( other.grads.constBegin().value() );
    }
    else
        return true;
}

/** Add all of the component gradients from 'other' onto this set. Note that
    this and 'other' must be compatible or else an error will be raised */
ComponentGradients& ComponentGradients::operator+=(const ComponentGradients &other)
{
    if (not this->isCompatible(other))
        throw SireError::incompatible_error( QObject::tr(
                "You cannot add the ComponentGradients %1 and %2 together as they "
                "are not compatible.")
                    .arg(this->toString(), other.toString()), CODELOC );

    if (this->isEmpty())
    {
        this->operator=(other);
    }
    else if (not other.isEmpty())
    {
        for (QMap<double,FreeEnergyMonitor>::const_iterator it = other.grads.constBegin();
             it != other.grads.constEnd();
             ++it)
        {
            if (grads.contains(it.key()))
            {
                grads[it.key()] += it.value();
            }
            else
            {
                grads.insert(it.key(), it.value());
            }
        }
    }

    return *this;
}

/** Addition operator. Note that you can only add to ComponentGradients objects
    together if they are compatible (this->isCompatible(other) is true). If they
    are not, then an exception will be raised */
ComponentGradients ComponentGradients::operator+(const ComponentGradients &other) const
{
    ComponentGradients ret(*this);
    ret += other;
    return ret;
}

/** Merge together all of the passed gradients. Note that they must all be compatible
    with one another, otherwise an exception will be raised */
ComponentGradients ComponentGradients::merge(const QList<ComponentGradients> &gradients)
{
    if (gradients.isEmpty())
        return ComponentGradients();
    else if (gradients.count() == 1)
        return gradients.at(0);
    else
    {
        ComponentGradients merged = gradients.at(0);

        for (int i=1; i<gradients.count(); ++i)
        {
            merged += gradients.at(i);
        }

        return merged;
    }
}

/** Return the temperature at which all of the components were collected */
Temperature ComponentGradients::temperature() const
{
    if (this->isEmpty())
        return Temperature(0);
    else
        return grads.constBegin().value().temperature();
}

/** Return the lambda values over which all of the components were collected */
QList<double> ComponentGradients::lambdaValues() const
{
    if (grads.isEmpty())
        return QList<double>();
    else
    {
        QList<double> lamvals = grads.keys();
        std::sort(lamvals.begin(), lamvals.end());

        return lamvals;
    }
}

/** Return the ith view that corresponds to the ith free energy component.
    Note that this returns the view in the numerically first (lowest) lambda
    value. Use viewAt(int i, double lamval) if you want to specify the lambda
    value from which you want to extract the view. */
PartialMolecule ComponentGradients::viewAt(int i) const
{
    if (grads.isEmpty())
        throw SireError::invalid_index( QObject::tr(
                "Cannot access the view for component %1, as this is an empty "
                "set of component gradients.")
                    .arg(i), CODELOC );

    QVector<PartialMolecule> views = grads.constBegin().value().referenceViews();

    return views.at( Index(i).map(views.count()) );
}

/** Return the ith view from lambda value 'lamval' that corresponds to the
    ith free energy component. */
PartialMolecule ComponentGradients::viewAt(int i, double lamval) const
{
    if (grads.isEmpty())
        throw SireError::invalid_index( QObject::tr(
                "Cannot access the view for component %1, as this is an empty "
                "set of component gradients.")
                    .arg(i), CODELOC );
    else if (not grads.contains(lamval))
        throw SireError::invalid_index( QObject::tr(
                "Cannot access the view for component %1 at lambda value %2, as this "
                "set of component gradient doesn't contain a view for this lambda value. "
                "(available lambda values are [ %3 ])")
                    .arg(i).arg(lamval).arg(Sire::toString(lambdaValues())), CODELOC );
    else
    {
        QVector<PartialMolecule> views = grads[lamval].referenceViews();

        return views.at( Index(i).map(views.count()) );
    }
}

/** Return the set of free energy gradients for the ith free energy component */
Gradients ComponentGradients::gradientsAt(int i) const
{
    if (grads.isEmpty())
        throw SireError::invalid_index( QObject::tr(
                "Cannot access the gradients for component %1, as this is an empty "
                "set of component gradients.")
                    .arg(i), CODELOC );

    QMap<double,FreeEnergyAverage> comp_grads;

    for (QMap<double,FreeEnergyMonitor>::const_iterator it = grads.constBegin();
         it != grads.constEnd();
         ++it)
    {
        QVector<FreeEnergyAverage> freenrgs = it.value().freeEnergies();
        comp_grads.insert(it.key(), freenrgs.at( Index(i).map(freenrgs.count()) ));
    }

    return Gradients(comp_grads,grads.constBegin().value().deltaLambda());
}

/** Return the set of coulomb free energy gradients for the ith free energy component */
Gradients ComponentGradients::coulombGradientsAt(int i) const
{
    if (grads.isEmpty())
        throw SireError::invalid_index( QObject::tr(
                "Cannot access the gradients for component %1, as this is an empty "
                "set of component gradients.")
                    .arg(i), CODELOC );

    QMap<double,FreeEnergyAverage> comp_grads;

    for (QMap<double,FreeEnergyMonitor>::const_iterator it = grads.constBegin();
         it != grads.constEnd();
         ++it)
    {
        QVector<FreeEnergyAverage> freenrgs = it.value().coulombFreeEnergies();
        comp_grads.insert(it.key(), freenrgs.at( Index(i).map(freenrgs.count()) ));
    }

    return Gradients(comp_grads,grads.constBegin().value().deltaLambda());
}

/** Return the set of LJ free energy gradients for the ith free energy component */
Gradients ComponentGradients::ljGradientsAt(int i) const
{
    if (grads.isEmpty())
        throw SireError::invalid_index( QObject::tr(
                "Cannot access the gradients for component %1, as this is an empty "
                "set of component gradients.")
                    .arg(i), CODELOC );

    QMap<double,FreeEnergyAverage> comp_grads;

    for (QMap<double,FreeEnergyMonitor>::const_iterator it = grads.constBegin();
         it != grads.constEnd();
         ++it)
    {
        QVector<FreeEnergyAverage> freenrgs = it.value().ljFreeEnergies();
        comp_grads.insert(it.key(), freenrgs.at( Index(i).map(freenrgs.count()) ));
    }

    return Gradients(comp_grads,grads.constBegin().value().deltaLambda());
}

/** Return the number of free energy components (number of molecule views whose
    free energy of interaction was recorded) */
int ComponentGradients::nComponents() const
{
    if (grads.isEmpty())
        return 0;
    else
        return grads.constBegin().value().referenceViews().count();
}

/** Return the number of lambda values over which the free energy components have
    been recorded */
int ComponentGradients::nLambdaValues() const
{
    return grads.count();
}

/** Return the number of samples used to form all of the average components */
qint64 ComponentGradients::nSamples() const
{
    qint64 n = 0;

    for (QMap<double,FreeEnergyMonitor>::const_iterator it = grads.constBegin();
         it != grads.constEnd();
         ++it)
    {
        QVector<FreeEnergyAverage> freenrgs = it.value().freeEnergies();

        foreach (const FreeEnergyAverage &freenrg, freenrgs)
        {
            n += freenrg.nSamples();
        }
    }

    return n;
}

/** Return the value of delta lambda used to approximate the free energy gradients */
double ComponentGradients::deltaLambda() const
{
    if (grads.isEmpty())
        return 0;
    else
        return grads.constBegin().value().deltaLambda();
}

/** Return the actual values of the free energy gradients of the ith component */
QVector<DataPoint> ComponentGradients::values(int i) const
{
    return gradientsAt(i).values();
}

/** Return the actual values of the coulomb free energy gradients of the ith component */
QVector<DataPoint> ComponentGradients::coulombValues(int i) const
{
    return coulombGradientsAt(i).values();
}

/** Return the actual values of the LJ free energy gradients of the ith component */
QVector<DataPoint> ComponentGradients::ljValues(int i) const
{
    return ljGradientsAt(i).values();
}

/** Return the raw data for all of the free energy components */
QMap<double,FreeEnergyMonitor> ComponentGradients::data() const
{
    return grads;
}

/** Integrate the free energy gradients of the ith component
    and return the resulting PMF */
TIPMF ComponentGradients::integrate(int i) const
{
    return gradientsAt(i).integrate();
}

/** Integrate the free energy gradients of the ith component to order 'order' and
    return the resulting PMF */
TIPMF ComponentGradients::integrate(int i, int order) const
{
    return gradientsAt(i).integrate(order);
}

/** Integrate the free energy gradients of the ith component between the range
    'range_min' to 'range_max', and return the resulting PMF */
TIPMF ComponentGradients::integrate(int i, double range_min, double range_max) const
{
    return gradientsAt(i).integrate(range_min, range_max);
}

/** Integrate the free energy gradients of the ith component to order 'order'
    between the range 'range_min' to 'range_max' and return the resulting PMF */
TIPMF ComponentGradients::integrate(int i, double range_min, double range_max, int order) const
{
    return gradientsAt(i).integrate(range_min, range_max, order);
}

/** Integrate the coulomb free energy gradients of the ith component
    and return the resulting PMF */
TIPMF ComponentGradients::integrateCoulomb(int i) const
{
    return coulombGradientsAt(i).integrate();
}

/** Integrate the coulomb free energy gradients of the ith component to order 'order' and
    return the resulting PMF */
TIPMF ComponentGradients::integrateCoulomb(int i, int order) const
{
    return coulombGradientsAt(i).integrate(order);
}

/** Integrate the coulomb free energy gradients of the ith component between the range
    'range_min' to 'range_max', and return the resulting PMF */
TIPMF ComponentGradients::integrateCoulomb(int i, double range_min, double range_max) const
{
    return coulombGradientsAt(i).integrate(range_min, range_max);
}

/** Integrate the coulomb free energy gradients of the ith component to order 'order'
    between the range 'range_min' to 'range_max' and return the resulting PMF */
TIPMF ComponentGradients::integrateCoulomb(int i, double range_min,
                                           double range_max, int order) const
{
    return coulombGradientsAt(i).integrate(range_min, range_max, order);
}

/** Integrate the LJ free energy gradients of the ith component
    and return the resulting PMF */
TIPMF ComponentGradients::integrateLJ(int i) const
{
    return ljGradientsAt(i).integrate();
}

/** Integrate the LJ free energy gradients of the ith component to order 'order' and
    return the resulting PMF */
TIPMF ComponentGradients::integrateLJ(int i, int order) const
{
    return ljGradientsAt(i).integrate(order);
}

/** Integrate the LJ free energy gradients of the ith component between the range
    'range_min' to 'range_max', and return the resulting PMF */
TIPMF ComponentGradients::integrateLJ(int i, double range_min, double range_max) const
{
    return ljGradientsAt(i).integrate(range_min, range_max);
}

/** Integrate the LJ free energy gradients of the ith component to order 'order'
    between the range 'range_min' to 'range_max' and return the resulting PMF */
TIPMF ComponentGradients::integrateLJ(int i, double range_min, double range_max, int order) const
{
    return ljGradientsAt(i).integrate(range_min, range_max, order);
}

//////////
////////// Implementation of TIComponents
//////////

static const RegisterMetaType<TIComponents> r_ti;

QDataStream &operator<<(QDataStream &ds, const TIComponents &ti)
{
    writeHeader(ds, r_ti, 2);

    SharedDataStream sds(ds);
    sds << ti.grads << ti.should_conserve_memory;

    return ds;
}

QDataStream &operator>>(QDataStream &ds, TIComponents &ti)
{
    VersionID v = readHeader(ds, r_ti);

    if (v == 2)
    {
        SharedDataStream sds(ds);
        sds >> ti.grads >> ti.should_conserve_memory;
    }
    else if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> ti.grads;
        ti.should_conserve_memory = false;
    }
    else
        throw version_error(v, "1,2", r_ti, CODELOC);

    return ds;
}

/** Constructor */
TIComponents::TIComponents(bool conserve_memory)
             : ConcreteProperty<TIComponents,Property>(),
               should_conserve_memory(conserve_memory)
{}

/** Construct from a single iteration's worth of gradients */
TIComponents::TIComponents(const QMap<double,FreeEnergyMonitor> &gradients,
                           bool conserve_memory)
             : ConcreteProperty<TIComponents,Property>(),
               should_conserve_memory(conserve_memory)
{
    grads.append( ComponentGradients(gradients,conserve_memory) );
}

/** Construct from a single iteration's worth of gradients */
TIComponents::TIComponents(const ComponentGradients &gradients,
                           bool conserve_memory)
             : ConcreteProperty<TIComponents,Property>(),
               should_conserve_memory(conserve_memory)
{
    if (conserve_memory)
    {
        ComponentGradients small_grads(gradients);
        small_grads.conserveMemory();
        grads.append(small_grads);
    }
    else
        grads.append(gradients);
}

/** Copy constructor */
TIComponents::TIComponents(const TIComponents &other)
             : ConcreteProperty<TIComponents,Property>(other),
               grads(other.grads), should_conserve_memory(other.should_conserve_memory)
{}

/** Destructor */
TIComponents::~TIComponents()
{}

/** Copy assignment operator */
TIComponents& TIComponents::operator=(const TIComponents &other)
{
    grads = other.grads;
    should_conserve_memory = other.should_conserve_memory;
    return *this;
}

/** Comparison operator */
bool TIComponents::operator==(const TIComponents &other) const
{
    return grads == other.grads and should_conserve_memory == other.should_conserve_memory;
}

/** Comparison operator */
bool TIComponents::operator!=(const TIComponents &other) const
{
    return not TIComponents::operator==(other);
}

const char* TIComponents::typeName()
{
    return QMetaType::typeName( qMetaTypeId<TIComponents>() );
}

const char* TIComponents::what() const
{
    return TIComponents::typeName();
}

QString TIComponents::toString() const
{
    return QObject::tr("TIComponents( nComponents() == %1, nIterations() == %2 "
                                     "nSamples() == %3, nLambdaValues() == %4 )")
                .arg(nComponents())
                .arg(nIterations())
                .arg(nSamples())
                .arg(nLambdaValues());
}

/** Conserve memory by sharing as much data as possible between the different iterations */
void TIComponents::conserveMemory()
{
    if (grads.isEmpty())
        return;

    //find the first non-empty set of gradients
    ComponentGradients first;

    for (int i=0; i<grads.count(); ++i)
    {
        if (not grads.at(i).isEmpty())
        {
            first = grads.at(i);
            break;
        }
    }

    if (first.isEmpty())
        return;

    //conserve memory in this first set
    first.conserveMemory();

    for (int i=0; i<grads.count(); ++i)
    {
        if (not grads.at(i).isEmpty())
            grads[i].conserveMemory(first);
    }

    should_conserve_memory = true;
}

/** Set the gradients for the ith iteration. Note that these must be compatible
    with the gradients of the other iterations */
void TIComponents::set(int i, const ComponentGradients &gradients)
{
    if (not gradients.isEmpty())
    {
        //make sure it is compatible
        for (int i=0; i<grads.count(); ++i)
        {
            if (not grads.at(i).isEmpty())
            {
                if (not grads.at(i).isCompatible(gradients))
                    throw SireError::incompatible_error( QObject::tr(
                            "Cannot add the passed set of gradients as they are incompatible "
                            "with those that have already been added."), CODELOC );

                break;
            }
        }
    }

    while (i >= grads.count())
    {
        grads.append( ComponentGradients() );
    }

    if (should_conserve_memory)
    {
        ComponentGradients first;

        for (int j=0; j<grads.count(); ++j)
        {
            if (not grads.at(j).isEmpty())
            {
                first = grads.at(j);
                break;
            }
        }

        ComponentGradients small_grads(gradients);

        if (first.isEmpty())
        {
            small_grads.conserveMemory();
        }
        else
        {
            small_grads.conserveMemory(first);
        }

        grads[i] = small_grads;
    }
    else
    {
        grads[i] = gradients;
    }
}

/** Set the gradients for the ith iteration. Note that these must be compatible
    with the gradients of the other iterations */
void TIComponents::set(int i, const QMap<double,FreeEnergyMonitor> &gradients)
{
    this->set(i, ComponentGradients(gradients));
}

/** Add the passed set of component gradients. Note that these must be compatible
    with any that are already in this set */
void TIComponents::add(const QMap<double,FreeEnergyMonitor> &gradients)
{
    return this->add( ComponentGradients(gradients) );
}

/** Add the passed set of component gradients. Note that these must be compatible
    with any that are already in this set */
void TIComponents::add(const ComponentGradients &gradients)
{
    return this->set( grads.count(), gradients );
}

/** Return the number of components in this collection (number of views) */
int TIComponents::nComponents() const
{
    if (grads.isEmpty())
        return 0;
    else
    {
        for (int i=0; i<grads.count(); ++i)
        {
            if (not grads.at(i).isEmpty())
                return grads.at(i).nComponents();
        }

        return 0;
    }
}

/** Return the number of iterations */
int TIComponents::nIterations() const
{
    return grads.count();
}

/** Return the number of lambda values */
int TIComponents::nLambdaValues() const
{
    return lambdaValues().count();
}

/** Return the total number of samples */
qint64 TIComponents::nSamples() const
{
    qint64 n = 0;

    for (int i=0; i<grads.count(); ++i)
    {
        n += grads.at(i).nSamples();
    }

    return n;
}

/** Return the number of iterations */
int TIComponents::count() const
{
    return grads.count();
}

/** Return the number of iterations */
int TIComponents::size() const
{
    return grads.count();
}

/** Whether or not this object conserves memory */
bool TIComponents::conservesMemory() const
{
    return should_conserve_memory;
}

/** Return a list of all lambda values that contain data */
QList<double> TIComponents::lambdaValues() const
{
    QMap<double,int> vals;

    for (int i=0; i<grads.count(); ++i)
    {
        foreach (double lam, grads.at(i).lambdaValues())
        {
            vals.insert(lam, 0);
        }
    }

    QList<double> lams = vals.keys();
    std::sort(lams.begin(), lams.end());

    return lams;
}

/** Return the ith set of ComponentGradients */
ComponentGradients TIComponents::operator[](int i) const
{
    return grads.at( Index(i).map(grads.count()) );
}

/** Return the ith set of ComponentGradients */
ComponentGradients TIComponents::at(int i) const
{
    return this->operator[](i);
}

/** Return the raw list of component gradients */
QList<ComponentGradients> TIComponents::gradients() const
{
    return grads;
}

/** Remove the data for iteration 'i'. This sets the data equal to ComponentGradients() */
void TIComponents::removeAt(int i)
{
    i = Index(i).map(grads.count());
    grads[i] = ComponentGradients();
}

/** Remove every iteration from 'start' to 'end' (inclusively). This sets
    the data equal to ComponentGradients() */
void TIComponents::removeRange(int start, int end)
{
    start = Index(start).map(grads.count());
    end = Index(end).map(grads.count());

    if (start > end)
        qSwap(start, end);

    for (int i = start; i <= end; ++i)
    {
        grads[i] = ComponentGradients();
    }
}

/** Remove all values from the histogram */
void TIComponents::clear()
{
    this->operator=( TIComponents() );
}

/** Merge (average) together the gradients from iteration "start" to iteration
    "end" inclusive */
ComponentGradients TIComponents::merge(int start, int end) const
{
    start = Index(start).map(grads.count());
    end = Index(end).map(grads.count());

    QList<ComponentGradients> set;

    for (int i=start; i<=end; ++i)
    {
        if (not grads.at(i).isEmpty())
            set.append( grads.at(i) );
    }

    return ComponentGradients::merge(set);
}

/** Merge together the gradients from the iterations with the passed indicies */
ComponentGradients TIComponents::merge(QList<int> indicies) const
{
    QList<ComponentGradients> set;

    foreach (int idx, indicies)
    {
        int i = Index(idx).map(grads.count());

        if (not grads.at(i).isEmpty())
            set.append( grads.at(i) );
    }

    return ComponentGradients::merge(set);
}

/** Return a list of Gradients that represents the rolling average over 'niterations'
    iterations over this TI data set. If this data set contains 100 iterations, and
    we calculate the rolling average over 50 iterations, then the returned Gradients
    will be the average from 1-50, then 2-51, 3-52.....51-100 */
QList<ComponentGradients> TIComponents::rollingAverage(int niterations) const
{
    if (grads.isEmpty())
        return QList<ComponentGradients>();

    QList<ComponentGradients> merged;

    if (niterations >= grads.count())
    {
        ComponentGradients m = this->merge(0,-1);

        if (not m.isEmpty())
            merged.append(m);
    }
    else if (niterations <= 1)
    {
        for (int i=0; i<grads.count(); ++i)
        {
            if (not grads.at(i).isEmpty())
                merged.append(grads.at(i));
        }
    }
    else
    {
        QList<ComponentGradients> set;

        int i = 0;

        for (i=0; i<grads.count(); ++i)
        {
            if (not grads.at(i).isEmpty())
            {
                set.append(grads.at(i));

                if (set.count() == niterations)
                    break;
            }
        }

        if (not set.isEmpty())
            merged.append( ComponentGradients::merge(set) );

        for (i = i+1; i<grads.count(); ++i)
        {
            if (not grads.at(i).isEmpty())
            {
                set.removeFirst();
                set.append(grads.at(i));
                merged.append( ComponentGradients::merge(set) );
            }
        }
    }

    return merged;
}
