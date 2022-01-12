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

#include "ti.h"

#include "SireID/index.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"
#include "SireStream/registeralternativename.h"

#include "third_party/regress.h"  // CONDITIONAL_INCLUDE

#include <cmath>

#include "tostring.h"

using namespace SireAnalysis;
using namespace SireMaths;
using namespace SireID;
using namespace SireBase;
using namespace SireStream;
using namespace SireUnits::Dimension;

/////////
///////// Implementation of TIPMF
/////////

static const RegisterMetaType<TIPMF> r_tipmf;
static const RegisterAlternativeName<TIPMF> r_alttipmf("Soiree::TIPMF");

QDataStream &operator<<(QDataStream &ds, const TIPMF &pmf)
{
    writeHeader(ds, r_tipmf, 1);

    SharedDataStream sds(ds);

    sds << pmf.grads << pmf.range_min << pmf.range_max << pmf.npoly
        << static_cast<const PMF&>(pmf);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, TIPMF &pmf)
{
    VersionID v = readHeader(ds, r_tipmf);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> pmf.grads >> pmf.range_min >> pmf.range_max >> pmf.npoly
            >> static_cast<PMF&>(pmf);

        pmf.recalculate();
    }
    else
        throw version_error(v, "1", r_tipmf, CODELOC);

    return ds;
}

/** Construct a PMF that will use 10 polynomials to fit
    and integrate the gradients between 0 and 1 */
TIPMF::TIPMF() : ConcreteProperty<TIPMF,PMF>(),
                 range_min(0), range_max(1), quad_value(0), npoly(10)
{}

/** Construct a PMF that will use the passed number of polynomials
    to fit and integrate the gradients between 0 and 1 */
TIPMF::TIPMF(int order) : ConcreteProperty<TIPMF,PMF>(),
                          range_min(0), range_max(1), quad_value(0), npoly(order)
{
    if (order <= 0)
        npoly = 1;
    else if (order > 100)
        npoly = 100;
}

/** Construct a PMF that will use 10 polynomials to fit and integrate
    the gradients in the passed range */
TIPMF::TIPMF(double rmin, double rmax)
      : ConcreteProperty<TIPMF,PMF>(),
        range_min(rmin), range_max(rmax), quad_value(0), npoly(10)
{
    if (range_min > range_max)
        qSwap(range_min, range_max);
}

/** Construct a PMF that will use the passed number of polynomials to fit and
    integrate the gradients in the passed range */
TIPMF::TIPMF(double min, double max, int order)
      : ConcreteProperty<TIPMF,PMF>(),
        range_min(min), range_max(max), quad_value(0), npoly(order)
{
    if (order <= 0)
        npoly = 1;
    else if (order > 100)
        npoly = 100;

    if (range_min > range_max)
        qSwap(range_min, range_max);
}

/** Set the order (number of polynomials) to fit the gradients for
    PMF integration */
void TIPMF::setOrder(int order)
{
    if (order <= 0)
        npoly = 1;
    else if (order > 100)
        npoly = 100;
    else
        npoly = order;

    recalculate();
}

/** Set the range of integration */
void TIPMF::setRange(double min_x, double max_x)
{
    range_min = min_x;
    range_max = max_x;

    if (range_min > range_max)
        qSwap(range_min, range_max);

    recalculate();
}

/** Set the raw gradients to be integrated */
void TIPMF::setGradients(const QVector<DataPoint> &gradients)
{
    grads = gradients;

    //ensure that the gradients are in sorted x numerical order
    bool sorted = false;

    while (not sorted)
    {
        sorted = true;

        //bubble sort - compare neighbours and swap if in the wrong order
        for (int i=1; i<grads.count(); ++i)
        {
            if (grads[i-1].x() > grads[i].x())
            {
                qSwap(grads[i-1], grads[i]);
                sorted = false;
            }
        }
    }

    recalculate();
}

/** Return the order (number of polynomials) used to integrate
    the gradients to get the PMF */
int TIPMF::order() const
{
    return npoly;
}

/** Return the minimum value of the range of integration */
double TIPMF::rangeMin() const
{
    return range_min;
}

/** Return the maximum value of the range of integration */
double TIPMF::rangeMax() const
{
    return range_max;
}

/** Internal function used to fit the raw gradients to a set
    of polynomials and to then integrate those to obtain the PMF */
void TIPMF::recalculate()
{
    if (grads.isEmpty())
    {
        smoothed_grads = QVector<DataPoint>();
        quad_value = 0;
        setValues(QVector<DataPoint>());
        return;
    }

    QVector<DataPoint> vals;

    //we will integrate these gradients using curve fitting to the underlying
    //gradients - this uses the "regress" code in third_party/regress.h
    // (author Conrad Shyu)
    std::list<stREGRESS> regress_grads;
    std::list<stREGRESS> regress_grads_plus_endpoints;
    std::list<stREGRESS> max_regress_grads;
    std::list<stREGRESS> min_regress_grads;
    std::list<stREGRESS> max_max_regress_grads;
    std::list<stREGRESS> min_min_regress_grads;
    {
        stREGRESS val;

        foreach (const DataPoint &grad, grads)
        {
            val.x = grad.x();
            val.y = grad.y();

            regress_grads.push_back(val);
            regress_grads_plus_endpoints.push_back(val);

            val.y = grad.y() + grad.yMinError();
            max_regress_grads.push_back(val);
            val.y = grad.y() + grad.yMaxError();
            max_max_regress_grads.push_back(val);

            val.y = grad.y() - grad.yMinError();
            min_regress_grads.push_back(val);
            val.y = grad.y() - grad.yMaxError();
            min_min_regress_grads.push_back(val);
        }

        //if the range extends before the first available gradient, then
        //we assume that this gradient is constant in this range
        if (grads.first().x() > range_min)
        {
            const DataPoint &grad = grads.first();

            val.x = range_min;
            val.y = grad.y();

            regress_grads_plus_endpoints.push_front(val);

            val.y = grad.y() + grad.yMinError();
            max_regress_grads.push_front(val);
            val.y = grad.y() + grad.yMaxError();
            max_max_regress_grads.push_front(val);

            val.y = grad.y() - grad.yMinError();
            min_regress_grads.push_front(val);
            val.y = grad.y() - grad.yMaxError();
            min_min_regress_grads.push_front(val);
        }

        //similarly, if the range extends beyond the last datapoint, then
        //continue this gradient
        if (grads.last().x() < range_max)
        {
            const DataPoint &grad = grads.last();

            val.x = range_max;
            val.y = grad.y();

            regress_grads_plus_endpoints.push_back(val);

            val.y = grad.y() + grad.yMinError();
            max_regress_grads.push_back(val);
            val.y = grad.y() + grad.yMaxError();
            max_max_regress_grads.push_back(val);

            val.y = grad.y() - grad.yMinError();
            min_regress_grads.push_back(val);
            val.y = grad.y() - grad.yMaxError();
            min_min_regress_grads.push_back(val);
        }
    }

    Regress regress(regress_grads, npoly);
    Regress regress_plus_endpoints(regress_grads_plus_endpoints, npoly);
    Regress max_regress(max_regress_grads, npoly);
    Regress max_max_regress(max_max_regress_grads, npoly);
    Regress min_regress(min_regress_grads, npoly);
    Regress min_min_regress(min_min_regress_grads, npoly);

    //get the coefficients of the polynomial
    std::vector<double> coeffs = regress.GetPolynomial();
    std::vector<double> max_coeffs = max_regress.GetPolynomial();
    std::vector<double> max_max_coeffs = max_max_regress.GetPolynomial();
    std::vector<double> min_coeffs = min_regress.GetPolynomial();
    std::vector<double> min_min_coeffs = min_min_regress.GetPolynomial();

    //now calculate the values of this polynomial
    //across lambda to get the smoothed gradients
    if (range_min != range_max)
    {
        smoothed_grads.clear();

        double x = range_min;
        double y = 0;
        double max_y = 0;
        double max_max_y = 0;
        double min_y = 0;
        double min_min_y = 0;
        double step = (range_max - range_min) * 0.01;

        while (x < range_max)
        {
            y = 0;
            max_y = 0;
            max_max_y = 0;
            min_y = 0;
            min_min_y = 0;

            for (size_t i=0; i<coeffs.size(); ++i)
            {
                double power = i;
                y += ( std::pow( x, power ) * coeffs[ i ] );
                max_y += ( std::pow( x, power ) * max_coeffs[ i ] );
                max_max_y += ( std::pow( x, power ) * max_max_coeffs[ i ] );
                min_y += ( std::pow( x, power ) * min_coeffs[ i ] );
                min_min_y += ( std::pow( x, power ) * min_min_coeffs[ i ] );
            }

            double err = qMax( std::abs(max_y-y), std::abs(y-min_y) );
            double max_err = qMax( std::abs(max_max_y-y), std::abs(y-min_min_y) );

            smoothed_grads.append( DataPoint(x,y,0,err,0,max_err) );
            x += step;
        }

        x = range_max;
        y = 0;
        max_y = 0;
        max_max_y = 0;
        min_y = 0;
        min_min_y = 0;

        for (size_t i=0; i<coeffs.size(); ++i)
        {
            double power = i;
            y += ( std::pow( x, power ) * coeffs[ i ] );
            max_y += ( std::pow( x, power ) * max_coeffs[ i ] );
            max_max_y += ( std::pow( x, power ) * max_max_coeffs[ i ] );
            min_y += ( std::pow( x, power ) * min_coeffs[ i ] );
            min_min_y += ( std::pow( x, power ) * min_min_coeffs[ i ] );
        }

        double err = qMax( std::abs(max_y-y), std::abs(y-min_y) );
        double max_err = qMax( std::abs(max_max_y-y), std::abs(y-min_min_y) );

        smoothed_grads.append( DataPoint(x,y,0,err,0,max_err) );
    }

    //now numerically integrate the gradients
    //across lambda to get the PMF
    if (range_min != range_max)
    {
        vals.clear();

        double x = range_min;
        double y = 0;
        double max_y = 0;
        double max_max_y = 0;
        double min_y = 0;
        double min_min_y = 0;
        double step = (range_max - range_min) * 0.01;

        vals.append( DataPoint(x,y) );

        double xmin, xmax;

        while (x < range_max)
        {
            xmin = x;
            xmax = x+step;

            double area = 0;
            double max_area = 0;
            double max_max_area = 0;
            double min_area = 0;
            double min_min_area = 0;

            for (size_t i=0; i<coeffs.size(); ++i)
            {
                double power = double(i + 1);

                area += ( ( pow( xmax, power ) / power ) * coeffs[ i ] -
                          ( pow( xmin, power ) / power ) * coeffs[ i ] );

                max_area += ( ( pow( xmax, power ) / power ) * max_coeffs[ i ] -
                              ( pow( xmin, power ) / power ) * max_coeffs[ i ] );

                max_max_area += ( ( pow( xmax, power ) / power ) * max_max_coeffs[ i ] -
                                  ( pow( xmin, power ) / power ) * max_max_coeffs[ i ] );

                min_area += ( ( pow( xmax, power ) / power ) * min_coeffs[ i ] -
                              ( pow( xmin, power ) / power ) * min_coeffs[ i ] );

                min_min_area += ( ( pow( xmax, power ) / power ) * min_min_coeffs[ i ] -
                                  ( pow( xmin, power ) / power ) * min_min_coeffs[ i ] );
            }

            y += area;
            max_y += max_area;
            max_max_y += max_max_area;
            min_y += min_area;
            min_min_y += min_min_area;

            double err = qMax( std::abs(max_y-y), std::abs(y-min_y) );
            double max_err = qMax( std::abs(max_max_y-y), std::abs(y-min_min_y) );

            vals.append( DataPoint(x,y,0,err,0,max_err) );

            x = xmax;
        }

        xmax = range_max;

        if (x != xmax)
        {
            xmin = x;

            double area = 0;
            double max_area = 0;
            double max_max_area = 0;
            double min_area = 0;
            double min_min_area = 0;

            for (size_t i=0; i<coeffs.size(); ++i)
            {
                double power = double(i + 1);

                area += ( ( pow( xmax, power ) / power ) * coeffs[ i ] -
                          ( pow( xmin, power ) / power ) * coeffs[ i ] );

                max_area += ( ( pow( xmax, power ) / power ) * max_coeffs[ i ] -
                              ( pow( xmin, power ) / power ) * max_coeffs[ i ] );

                max_max_area += ( ( pow( xmax, power ) / power ) * max_max_coeffs[ i ] -
                                  ( pow( xmin, power ) / power ) * max_max_coeffs[ i ] );

                min_area += ( ( pow( xmax, power ) / power ) * min_coeffs[ i ] -
                              ( pow( xmin, power ) / power ) * min_coeffs[ i ] );

                min_min_area += ( ( pow( xmax, power ) / power ) * min_min_coeffs[ i ] -
                                  ( pow( xmin, power ) / power ) * min_min_coeffs[ i ] );
            }

            y += area;
            max_y += max_area;
            max_max_y += max_max_area;
            min_y += min_area;
            min_min_y += min_min_area;

            double err = qMax( std::abs(max_y-y), std::abs(y-min_y) );
            double max_err = qMax( std::abs(max_max_y-y), std::abs(y-min_min_y) );

            vals.append( DataPoint(x,y,0,err,0,max_err) );
        }

        quad_value = regress_plus_endpoints.DoQuadrature();
    }

    setValues(vals);
}

/** Copy constructor */
TIPMF::TIPMF(const TIPMF &other)
      : ConcreteProperty<TIPMF,PMF>(other),
        grads(other.grads), smoothed_grads(other.smoothed_grads),
        range_min(other.range_min), range_max(other.range_max),
        quad_value(other.quad_value), npoly(other.npoly)
{}

/** Destructor */
TIPMF::~TIPMF()
{}

/** Copy assignment operator */
TIPMF& TIPMF::operator=(const TIPMF &other)
{
    if (this != &other)
    {
        grads = other.grads;
        smoothed_grads = other.smoothed_grads;
        range_min = other.range_min;
        range_max = other.range_max;
        quad_value = other.quad_value;
        npoly = other.npoly;
    }

    return *this;
}

/** Comparison operator */
bool TIPMF::operator==(const TIPMF &other) const
{
    return grads == other.grads and range_min == other.range_min and
           range_max == other.range_max and npoly == other.npoly and
           PMF::operator==(other);
}

/** Comparison operator */
bool TIPMF::operator!=(const TIPMF &other) const
{
    return not operator==(other);
}

const char* TIPMF::what() const
{
    return TIPMF::typeName();
}

const char* TIPMF::typeName()
{
    return QMetaType::typeName( qMetaTypeId<TIPMF>() );
}

/** Return the free energy calculated using integration of the
    polynomial fitted to the gradients */
double TIPMF::integral() const
{
    if (this->isEmpty())
        return 0;
    else
        return values().back().y();
}

/** Return the free energy calculated using trapezium quadrature
    from the raw gradients */
double TIPMF::quadrature() const
{
    return quad_value;
}

QString TIPMF::toString() const
{
    if (this->isEmpty())
        return QString("TIPMF()");
    else
        return QString("TIPMF( { integral() == %1, quadrature() == %2, "
                       "deltaG() == %3, error() == %4 } )")
                    .arg(integral()).arg(quadrature())
                    .arg(deltaG()).arg(error());
}

/** Return the raw gradients used to calculate the PMF */
QVector<DataPoint> TIPMF::gradients() const
{
    return grads;
}

/** Return the smoothed (fitted) gradients used to calculate the PMF */
QVector<DataPoint> TIPMF::smoothedGradients() const
{
    return smoothed_grads;
}

/** Return a copy of the PMF where the gradients at the end points
    (the first and last gradients) have been removed. This can be used
    to estimate the effect of end-point error */
TIPMF TIPMF::dropEndPoints() const
{
    TIPMF ret(*this);

    QVector<DataPoint> reduced_grads = grads;

    if (not reduced_grads.isEmpty())
        reduced_grads.pop_front();

    if (not reduced_grads.isEmpty())
        reduced_grads.pop_back();

    ret.setGradients(reduced_grads);

    return ret;
}

/////////
///////// Implementation of Gradients
/////////

static const RegisterMetaType<Gradients> r_grads;
static const RegisterAlternativeName<Gradients> r_altgrads("Soiree::Gradients");

QDataStream &operator<<(QDataStream &ds, const Gradients &grads)
{
    writeHeader(ds, r_grads, 2);

    SharedDataStream sds(ds);

    sds << grads.analytic << grads.fwds << grads.bwds << grads.delta_lam;

    return ds;
}

QDataStream &operator>>(QDataStream &ds, Gradients &grads)
{
    VersionID v = readHeader(ds, r_grads);

    if (v == 2)
    {
        SharedDataStream sds(ds);

        sds >> grads.analytic >> grads.fwds >> grads.bwds >> grads.delta_lam;
    }
    else if (v == 1)
    {
        SharedDataStream sds(ds);

        grads.analytic.clear();

        sds >> grads.fwds >> grads.bwds >> grads.delta_lam;
    }
    else
        throw version_error(v, "1,2", r_grads, CODELOC);

    return ds;
}

/** Constructor */
Gradients::Gradients() : ConcreteProperty<Gradients,Property>(), delta_lam(0)
{}

void Gradients::checkSane() const
{
    if (fwds.isEmpty() and bwds.isEmpty())
        return;

    Temperature t;
    bool have_first = false;

    for (QMap<double,FreeEnergyAverage>::const_iterator it = fwds.constBegin();
         it != fwds.constEnd();
         ++it)
    {
        if (not have_first)
        {
            t = it.value().temperature();
            have_first = true;
        }
        else if (it.value().temperature() != t)
        {
            throw SireError::invalid_arg( QObject::tr(
                "You cannot construct a set of TI gradients using free energy "
                "averages collected at different temperatures. %1 vs. %2")
                    .arg(t.toString()).arg(it.value().temperature().toString()),
                        CODELOC );
        }
    }

    for (QMap<double,FreeEnergyAverage>::const_iterator it = bwds.constBegin();
         it != bwds.constEnd();
         ++it)
    {
        if (not have_first)
        {
            t = it.value().temperature();
            have_first = true;
        }
        else if (it.value().temperature() != t)
        {
            throw SireError::invalid_arg( QObject::tr(
                "You cannot construct a set of TI gradients using free energy "
                "averages collected at different temperatures. %1 vs. %2")
                    .arg(t.toString()).arg(it.value().temperature().toString()),
                        CODELOC );
        }
    }
}

/** Construct from the passed analytic TI gradients */
Gradients::Gradients(const QMap<double,AverageAndStddev> &gradients)
          : ConcreteProperty<Gradients,Property>(),
            analytic(gradients), delta_lam(0)
{}

/** Construct from the passed full finite difference TI gradients */
Gradients::Gradients(const QMap<double,FreeEnergyAverage> &gradients)
          : ConcreteProperty<Gradients,Property>(),
            fwds(gradients), bwds(gradients), delta_lam(0)
{
    checkSane();
}

/** Construct from the passed finite difference TI gradients, using the passed
    value of delta lambda */
Gradients::Gradients(const QMap<double,FreeEnergyAverage> &gradients,
                     double delta_lambda)
          : ConcreteProperty<Gradients,Property>(),
            fwds(gradients), bwds(gradients), delta_lam(delta_lambda)
{
    if (delta_lam <= 0)
        throw SireError::invalid_arg( QObject::tr(
                "How can you have finite difference gradients with a value of "
                "delta lambda that is less than or equal to zero? %1")
                    .arg(delta_lam), CODELOC );

    checkSane();
}

/** Construct from the passed finite difference TI forwards and backwards
    gradients, using the passed value of delta lambda. Note that the
    forwards gradients should be the zwanzig free energies from
    lambda->lambda+delta_lambda, while the backwards gradients should
    be the zwanzig free energies from lambda-delta_lambda->lambda */
Gradients::Gradients(const QMap<double,FreeEnergyAverage> &forwards,
                     const QMap<double,FreeEnergyAverage> &backwards,
                     double delta_lambda)
          : ConcreteProperty<Gradients,Property>(),
            fwds(forwards), bwds(backwards), delta_lam(delta_lambda)
{
    if (delta_lam <= 0)
        throw SireError::invalid_arg( QObject::tr(
                "How can you have finite difference gradients with a value of "
                "delta lambda that is less than or equal to zero? %1")
                    .arg(delta_lam), CODELOC );

    //if the lambda values don't match, then copy the gradients from the
    //partners set
    if (backwards.isEmpty())
        bwds = fwds;

    else if (forwards.isEmpty())
        fwds = bwds;

    else
    {
        for (QMap<double,FreeEnergyAverage>::const_iterator it = bwds.constBegin();
             it != bwds.constEnd();
             ++it)
        {
            if (not fwds.contains(it.key()))
                fwds.insert( it.key(), it.value() );
        }

        for (QMap<double,FreeEnergyAverage>::const_iterator it = fwds.constBegin();
             it != fwds.constEnd();
             ++it)
        {
            if (not bwds.contains(it.key()))
                bwds.insert( it.key(), it.value() );
        }
    }

    checkSane();
}

/** Copy constructor */
Gradients::Gradients(const Gradients &other)
          : ConcreteProperty<Gradients,Property>(other),
            analytic(other.analytic),
            fwds(other.fwds), bwds(other.bwds), delta_lam(other.delta_lam)
{}

/** Destructor */
Gradients::~Gradients()
{}

/** Return whether or not this is empty */
bool Gradients::isEmpty() const
{
    return (analytic.isEmpty() and fwds.isEmpty() and bwds.isEmpty());
}

/** Copy assignment operator */
Gradients& Gradients::operator=(const Gradients &other)
{
    if (this != &other)
    {
        analytic = other.analytic;
        fwds = other.fwds;
        bwds = other.bwds;
        delta_lam = other.delta_lam;
    }

    return *this;
}

/** Comparison operator */
bool Gradients::operator==(const Gradients &other) const
{
    return analytic == other.analytic and fwds == other.fwds and
           bwds == other.bwds and delta_lam == other.delta_lam;
}

/** Comparison operator */
bool Gradients::operator!=(const Gradients &other) const
{
    return not operator==(other);
}

/** Return the value of the gradient at the passed lambda value */
MolarEnergy Gradients::operator[](double lam) const
{
    return gradient(lam);
}

/** Return the temperature at which the gradients were collected */
Temperature Gradients::temperature() const
{
    if (this->isEmpty())
        return Temperature(0);
    else
    {
        if (not analytic.isEmpty())
            return Temperature(0);
        else if (not fwds.isEmpty())
            return fwds.constBegin()->temperature();
        else
            return bwds.constBegin()->temperature();
    }
}

QString Gradients::toString() const
{
    if (not analytic.isEmpty())
    return QObject::tr("Gradients( nLambdaValues() == %1, nSamples() == %2 )")
                .arg(nLambdaValues()).arg(nSamples());
    else
        return QObject::tr("Gradients( nLambdaValues() == %1, nSamples() == %2, "
                           "temperature() == %3 )")
                .arg(nLambdaValues()).arg(nSamples()).arg(temperature().toString());
}

const char* Gradients::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Gradients>() );
}

const char* Gradients::what() const
{
    return Gradients::typeName();
}

/** Self-addition operator */
Gradients& Gradients::operator+=(const Gradients &other)
{
    if (this->isEmpty())
    {
        this->operator=(other);
        return *this;
    }
    else if (other.isEmpty())
    {
        return *this;
    }
    else
    {
        if (delta_lam != other.delta_lam)
            throw SireError::incompatible_error( QObject::tr(
                    "Cannot combine together these free energy Gradients "
                    "as the values of delta lambda are different. %1 vs. %2.")
                        .arg(delta_lam).arg(other.delta_lam), CODELOC );

        if (temperature() != other.temperature())
            throw SireError::incompatible_error( QObject::tr(
                "Cannot combine together these two free energy Gradients "
                "as the temperatures are different. %1 vs. %2.")
                    .arg(temperature().toString())
                    .arg(other.temperature().toString()), CODELOC );

        if (analytic.isEmpty())
        {
            QMap<double,FreeEnergyAverage> new_fwds = fwds;
            QMap<double,FreeEnergyAverage> new_bwds = bwds;

            for (QMap<double,FreeEnergyAverage>::const_iterator it = other.fwds.constBegin();
                 it != other.fwds.constEnd();
                 ++it)
            {
                double lam = it.key();
                const FreeEnergyAverage &grad = it.value();

                if (grad.nSamples() != 0)
                {
                    if (new_fwds.contains(lam))
                    {
                        new_fwds[lam] += grad;
                    }
                    else
                    {
                        new_fwds.insert(lam, grad);
                    }
                }
            }

            for (QMap<double,FreeEnergyAverage>::const_iterator it = other.bwds.constBegin();
                 it != other.bwds.constEnd();
                 ++it)
            {
                double lam = it.key();
                const FreeEnergyAverage &grad = it.value();

                if (grad.nSamples() != 0)
                {
                    if (new_bwds.contains(lam))
                    {
                        new_bwds[lam] += grad;
                    }
                    else
                    {
                        new_bwds.insert(lam, grad);
                    }
                }
            }

            if (new_fwds == new_bwds)
            {
                fwds = new_fwds;
                bwds = new_fwds;
            }
            else
            {
                fwds = new_fwds;
                bwds = new_bwds;
            }
        }
        else
        {
            QMap<double,AverageAndStddev> new_analytic = analytic;

            for (QMap<double,AverageAndStddev>::const_iterator it = other.analytic.constBegin();
                 it != other.analytic.constEnd();
                 ++it)
            {
                double lam = it.key();
                const AverageAndStddev &grad = it.value();

                if (grad.nSamples() != 0)
                {
                    if (new_analytic.contains(lam))
                    {
                        new_analytic[lam] += grad;
                    }
                    else
                    {
                        new_analytic.insert(lam, grad);
                    }
                }
            }

            analytic = new_analytic;
        }

        return *this;
    }
}

/** Addition operator */
Gradients Gradients::operator+(const Gradients &other) const
{
    Gradients ret(*this);
    ret += other;
    return ret;
}

/** Merge together the passed list of Gradients into a single object.
    Note that all of the passed gradients must be compatible, e.g.
    have the same temperature and delta lambda values */
Gradients Gradients::merge(const QList<Gradients> &gradients)
{
    if (gradients.isEmpty())
        return Gradients();
    else if (gradients.count() == 1)
        return gradients.at(0);

    Gradients ret = gradients.at(0);

    for (int i=1; i<gradients.count(); ++i)
    {
        ret += gradients.at(i);
    }

    return ret;
}

/** Return the (sorted) list of the lambda values */
QList<double> Gradients::lambdaValues() const
{
    if (not analytic.isEmpty())
    {
        QList<double> lams = analytic.keys();
        std::sort(lams.begin(), lams.end());
        return lams;
    }
    else if (delta_lam == 0 or fwds == bwds)
    {
        QList<double> lams = fwds.keys();
        std::sort(lams.begin(), lams.end());
        return lams;
    }
    else
    {
        QMap<double,int> lams;

        foreach (double lam, fwds.keys())
        {
            lams.insert(lam, 1);
        }

        foreach (double lam, bwds.keys())
        {
            lams.insert(lam, 1);
        }

        QList<double> l = lams.keys();
        std::sort(l.begin(), l.end());

        return l;
    }
}

/** Return the (sorted) list of all lambda values */
QList<double> Gradients::keys() const
{
    return lambdaValues();
}

/** Return the total number of lambda values */
int Gradients::nLambdaValues() const
{
    return lambdaValues().count();
}

/** Return the total number of samples used to calculate these gradients */
qint64 Gradients::nSamples() const
{
    qint64 n = 0;

    if (not analytic.isEmpty())
    {
        for (QMap<double,AverageAndStddev>::const_iterator it = analytic.constBegin();
             it != analytic.constEnd();
             ++it)
        {
            n += it.value().nSamples();
        }
    }
    else if (delta_lam == 0 or fwds == bwds)
    {
        for (QMap<double,FreeEnergyAverage>::const_iterator it = fwds.constBegin();
             it != fwds.constEnd();
             ++it)
        {
            n += it.value().nSamples();
        }
    }
    else
    {
        for (QMap<double,FreeEnergyAverage>::const_iterator it = fwds.constBegin();
             it != fwds.constEnd();
             ++it)
        {
            n += it.value().nSamples();
        }

        for (QMap<double,FreeEnergyAverage>::const_iterator it = bwds.constBegin();
             it != bwds.constEnd();
             ++it)
        {
            n += it.value().nSamples();
        }
    }

    return n;
}

/** Return the forwards gradient for the passed lambda value */
MolarEnergy Gradients::forwards(double lam) const
{
    if (not analytic.isEmpty())
    {
        if (not analytic.contains(lam))
            throw SireError::invalid_key( QObject::tr(
                    "There is no gradient value at lambda == %1. Available lambda values "
                    "are %2.").arg(lam).arg( Sire::toString(lambdaValues()) ), CODELOC );

        return MolarEnergy( analytic[lam].average() );
    }
    else
    {
        if (not fwds.contains(lam))
            throw SireError::invalid_key( QObject::tr(
                    "There is no gradient value at lambda == %1. Available lambda values "
                    "are %2.").arg(lam).arg( Sire::toString(lambdaValues()) ), CODELOC );

        if (delta_lam == 0)
            //pure TI, so just need the normal average energy
            return MolarEnergy( fwds[lam].histogram().mean() );
        else
        {
            //finite difference TI, so take the average of the forwards and backwards
            //values and divide by delta lambda
            return MolarEnergy( fwds[lam].fepFreeEnergy() / delta_lam );
        }
    }
}

/** Return the backwards gradient for the passed lambda value */
MolarEnergy Gradients::backwards(double lam) const
{
    if (not analytic.isEmpty())
    {
        if (not analytic.contains(lam))
            throw SireError::invalid_key( QObject::tr(
                    "There is no gradient value at lambda == %1. Available lambda values "
                    "are %2.").arg(lam).arg( Sire::toString(lambdaValues()) ), CODELOC );

        return MolarEnergy( analytic[lam].average() );
    }
    else
    {
        if (not fwds.contains(lam))
            throw SireError::invalid_key( QObject::tr(
                    "There is no gradient value at lambda == %1. Available lambda values "
                    "are %2.").arg(lam).arg( Sire::toString(lambdaValues()) ), CODELOC );

        if (delta_lam == 0)
            //pure TI, so just need the normal average energy
            return MolarEnergy( fwds[lam].histogram().mean() );
        else
        {
            //finite difference TI, so take the average of the forwards and backwards
            //values and divide by delta lambda
            return MolarEnergy( bwds[lam].fepFreeEnergy() / delta_lam );
        }
    }
}

/** Return the gradient at the passed lambda value. This is the
    average of the forwards and backwards gradient if finite difference
    is used */
MolarEnergy Gradients::gradient(double lam) const
{
    if (not analytic.isEmpty())
    {
        if (not analytic.contains(lam))
            throw SireError::invalid_key( QObject::tr(
                    "There is no gradient value at lambda == %1. Available lambda values "
                    "are %2.").arg(lam).arg( Sire::toString(lambdaValues()) ), CODELOC );

        return MolarEnergy( analytic[lam].average() );
    }
    else
    {
        if (not fwds.contains(lam))
            throw SireError::invalid_key( QObject::tr(
                    "There is no gradient value at lambda == %1. Available lambda values "
                    "are %2.").arg(lam).arg( Sire::toString(lambdaValues()) ), CODELOC );

        if (delta_lam == 0)
            //pure TI, so just need the normal average energy
            return MolarEnergy( fwds[lam].histogram().mean() );
        else
        {
            //finite difference TI, so take the average of the forwards and backwards
            //values and divide by delta lambda
            return MolarEnergy( 0.5 * (fwds[lam].fepFreeEnergy() + bwds[lam].fepFreeEnergy())
                                    / delta_lam );
        }
    }
}

/** Return the value of delta lambda. This will be zero if these are
    pure TI gradients */
double Gradients::deltaLambda() const
{
    return delta_lam;
}

/** Return the raw data for the analytic gradients */
QMap<double,AverageAndStddev> Gradients::analyticData() const
{
    return analytic;
}

/** Return the raw data for the forwards free energy gradients */
QMap<double,FreeEnergyAverage> Gradients::forwardsData() const
{
    return fwds;
}

/** Return the raw data for the backwards free energy gradients */
QMap<double,FreeEnergyAverage> Gradients::backwardsData() const
{
    return bwds;
}

/** Return the values of the gradients as data points. This returns the
    average of the forwards and backwards gradients, with errors calculated
    based on both the difference between the forwards and backwards values,
    and the 90% confidence level of the average of gradients */
QVector<DataPoint> Gradients::values() const
{
    QList<double> lamvals = this->lambdaValues();

    if (lamvals.isEmpty())
        return QVector<DataPoint>();

    QVector<DataPoint> points( lamvals.count() );

    if (not analytic.isEmpty())
    {
        for (int i=0; i<lamvals.count(); ++i)
        {
            double lam = lamvals[i];

            const AverageAndStddev &avg = *(analytic.constFind(lam));
            points[i] = DataPoint(lam, avg.average(), 0, avg.standardError(90));
        }
    }
    else
    {
        for (int i=0; i<lamvals.count(); ++i)
        {
            double lam = lamvals[i];

            if (delta_lam == 0)
            {
                //pure TI data
                if (fwds.contains(lam))
                {
                    const FreeEnergyAverage &avg = *(fwds.constFind(lam));
                    points[i] = DataPoint(lam, avg.histogram().mean(), 0,
                                               avg.histogram().standardError(90));
                }
            }
            else
            {
                const FreeEnergyAverage *fwdsavg = 0;
                const FreeEnergyAverage *bwdsavg = 0;

                if (fwds.contains(lam))
                    fwdsavg = &(*(fwds.constFind(lam)));

                if (bwds.contains(lam))
                    bwdsavg = &(*(bwds.constFind(lam)));

                if (bwdsavg == 0)
                    bwdsavg = fwdsavg;
                else if (fwdsavg == 0)
                    fwdsavg = bwdsavg;

                if (fwdsavg == bwdsavg or (*fwdsavg == *bwdsavg))
                {
                    double fwdsval = fwdsavg->average() / delta_lam;
                    double fwdserr = fwdsavg->histogram().standardError(90) / delta_lam;

                    double val = fwdsval;
                    double maxerr = fwdserr;

                    if (fwdserr > maxerr)
                        qSwap(fwdserr, maxerr);

                    points[i] = DataPoint(lam, val, 0, fwdserr, 0, maxerr);
                }
                else
                {
                    double fwdsval = fwdsavg->average() / delta_lam;
                    double bwdsval = bwdsavg->average() / delta_lam;

                    double fwdserr = fwdsavg->histogram().standardError(90) / delta_lam;
                    double bwdserr = bwdsavg->histogram().standardError(90) / delta_lam;

                    double val = 0.5 * (fwdsval + bwdsval);
                    double minerr = std::abs( fwdsval - bwdsval );
                    double maxerr = minerr + fwdserr + bwdserr;

                    points[i] = DataPoint(lam, val, 0, minerr, 0, maxerr);
                }
            }
        }
    }

    return points;
}

/** Return the values of the forwards gradients as data points. This returns the
    average forwards gradient for each lambda value, together with the
    standard error at the 90% confidence level */
QVector<DataPoint> Gradients::forwardsValues() const
{
    if (not analytic.isEmpty())
        return values();

    QList<double> lamvals = this->lambdaValues();

    if (lamvals.isEmpty())
        return QVector<DataPoint>();

    QVector<DataPoint> points( lamvals.count() );

    for (int i=0; i<lamvals.count(); ++i)
    {
        double lam = lamvals[i];

        if (fwds.contains(lam))
        {
            if (delta_lam == 0)
            {
                //pure TI data
                const FreeEnergyAverage &avg = *(fwds.constFind(lam));
                points[i] = DataPoint(lam, avg.histogram().mean(), 0,
                                           avg.histogram().standardError(90));
            }
            else
            {
                const FreeEnergyAverage &fwdsavg = *(fwds.constFind(lam));

                double fwdsval = fwdsavg.fepFreeEnergy() / delta_lam;
                double fwdserr = fwdsavg.histogram().standardError(90) / delta_lam;

                double val = fwdsval;
                double maxerr = fwdserr;

                if (fwdserr > maxerr)
                    qSwap(fwdserr, maxerr);

                points[i] = DataPoint(lam, val, 0, fwdserr, 0, maxerr);
            }
        }
    }

    return points;
}

/** Return the values of the backwards gradients as data points. This returns the
    average backwards gradient for each lambda value, together with the
    standard error at the 90% confidence level */
QVector<DataPoint> Gradients::backwardsValues() const
{
    if (not analytic.isEmpty())
        return values();

    QList<double> lamvals = this->lambdaValues();

    if (lamvals.isEmpty())
        return QVector<DataPoint>();

    QVector<DataPoint> points( lamvals.count() );

    for (int i=0; i<lamvals.count(); ++i)
    {
        double lam = lamvals[i];

        if (bwds.contains(lam))
        {
            if (delta_lam == 0)
            {
                //pure TI data
                const FreeEnergyAverage &avg = *(bwds.constFind(lam));
                points[i] = DataPoint(lam, avg.histogram().mean(), 0,
                                           avg.histogram().standardError(90));
            }
            else
            {
                const FreeEnergyAverage &bwdsavg = *(bwds.constFind(lam));

                double bwdsval = bwdsavg.fepFreeEnergy() / delta_lam;
                double bwdserr = bwdsavg.histogram().standardError(90) / delta_lam;

                double val = bwdsval;
                double maxerr = bwdserr;

                if (bwdserr > maxerr)
                    qSwap(bwdserr, maxerr);

                points[i] = DataPoint(lam, val, 0, bwdserr, 0, maxerr);
            }
        }
    }

    return points;
}

/** Integrate these gradients between 'range_min' to 'range_max' using
    a polynomial of passed order and return them as a potential of mean force (PMF) */
TIPMF Gradients::integrate(double range_min, double range_max, int order) const
{
    TIPMF pmf(range_min, range_max, order);

    if (not this->isEmpty())
        pmf.setGradients(this->values());

    return pmf;
}

/** Integrate these gradients between 'range_min' to 'range_max' using
    a polynomial of order 10 and return them as a potential of mean force (PMF) */
TIPMF Gradients::integrate(double range_min, double range_max) const
{
    return this->integrate(range_min, range_max, qMax(2, lambdaValues().count() - 2));
}

/** Integrate these gradients between 0 and 1 using a polynomial
    of passed order and return them as a potential of mean force (PMF) */
TIPMF Gradients::integrate(int order) const
{
    return this->integrate(0, 1, order);
}

/** Integrate these gradients between 0 and 1 using a polynomial of
    order ngradients-2 and return them as a potential of mean force (PMF) */
TIPMF Gradients::integrate() const
{
    return this->integrate(0, 1, qMax(2, lambdaValues().count() - 2));
}

/////////
///////// Implementation of TI
/////////

static const RegisterMetaType<TI> r_ti;
static const RegisterAlternativeName<TI> r_altti("Soiree::TI");

QDataStream &operator<<(QDataStream &ds, const TI &ti)
{
    writeHeader(ds, r_ti, 1);

    SharedDataStream sds(ds);
    sds << ti.grads;

    return ds;
}

QDataStream &operator>>(QDataStream &ds, TI &ti)
{
    VersionID v = readHeader(ds, r_ti);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> ti.grads;
    }
    else
        throw version_error(v, "1", r_ti, CODELOC);

    return ds;
}

/** Null constructor */
TI::TI() : ConcreteProperty<TI,Property>()
{}

/** Construct from the passed set of gradients */
TI::TI(const Gradients &gradients) : ConcreteProperty<TI,Property>()
{
    if (not gradients.isEmpty())
        grads.append(gradients);
}

/** Construct from the passed list of gradients from each iteration */
TI::TI(const QList<Gradients> &gradients) : ConcreteProperty<TI,Property>()
{
    foreach (const Gradients &g, gradients)
    {
        if (not g.isEmpty())
            grads.append(g);
    }
}

/** Copy constructor */
TI::TI(const TI &other) : ConcreteProperty<TI,Property>(other), grads(other.grads)
{}

/** Destructor */
TI::~TI()
{}

/** Copy assignment operator */
TI& TI::operator=(const TI &other)
{
    grads = other.grads;
    return *this;
}

/** Comparison operator */
bool TI::operator==(const TI &other) const
{
    return grads == other.grads;
}

/** Comparison operator */
bool TI::operator!=(const TI &other) const
{
    return not TI::operator==(other);
}

const char* TI::typeName()
{
    return QMetaType::typeName( qMetaTypeId<TI>() );
}

const char* TI::what() const
{
    return TI::typeName();
}

/** Return the raw list of gradients */
QList<Gradients> TI::gradients() const
{
    return grads;
}

/** Add the passed free energy gradients to the TI calculation. The gradients
    are in a dictionary, indexed by lambda value, and are analytic TI gradients */
void TI::add(const QMap<double,AverageAndStddev> &gradients)
{
    if (not gradients.isEmpty())
        grads.append( Gradients(gradients) );
}

/** Add the passed free energy gradients to the TI calculation. The gradients
    are in a dictionary, indexed by lambda value, and are "pure" TI gradients,
    i.e. they have been calculated exactly with an infinitesimal delta lambda */
void TI::add(const QMap<double,FreeEnergyAverage> &gradients)
{
    if (not gradients.isEmpty())
        grads.append( Gradients(gradients) );
}

/** Add the passed free energy gradients to the TI calcualtion. The gradients
    are in a dictionary, indexed by lambda value, and are the raw free energies
    calculated via the zwanzig equation as part of a finite-difference TI calculation.
    The value of delta lambda used must also be passed (so that we can then divide
    each gradient by delta lambda to get an approximation of the gradient) */
void TI::add(const QMap<double,FreeEnergyAverage> &gradients,
             double delta_lambda)
{
    if (not gradients.isEmpty())
        grads.append( Gradients(gradients, delta_lambda) );
}

/** Add the passed free energy gradients to the TI calculation. The gradients
    are in dictionaries, indexed by lambda values, and are the raw forwards and
    backwards free energies calculated via the zwanzig equation as part of
    a finite-difference TI calculation (backwards gradients calculated as the
    free energy from lambda-delta_lambda -> lambda, while forwards gradients calculated as the
    difference between lambda -> lambda+delta_lambda). The value of delta lambda
    must be passed (so that we can then divide each gradient by delta lambda to get
    an approximation of the true gradient) */
void TI::add(const QMap<double,FreeEnergyAverage> &forwards,
             const QMap<double,FreeEnergyAverage> &backwards,
             double delta_lambda)
{
    if (not (forwards.isEmpty() and backwards.isEmpty()))
    {
        grads.append( Gradients(forwards, backwards, delta_lambda) );
    }
}

void TI::add(const Gradients &gradients)
{
    if (not gradients.isEmpty())
        grads.append(gradients);
}

/** Return the number of iterations (number of sets of gradients that have been added) */
int TI::nIterations() const
{
    return grads.count();
}

/** Return all values of lambda that have data. The values are returned
    in numerical order */
QList<double> TI::lambdaValues() const
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

/** Return the number of lambda values */
int TI::nLambdaValues() const
{
    return lambdaValues().count();
}

/** Return the total number of samples in this calculation */
qint64 TI::nSamples() const
{
    qint64 n = 0;

    for (int i=0; i<grads.count(); ++i)
    {
        n += grads.at(i).nSamples();
    }

    return n;
}

/** Return the number of iterations */
int TI::count() const
{
    return nIterations();
}

/** Return the number of iterations */
int TI::size() const
{
    return nIterations();
}

/** Return the free energy gradient data for the ith iteration. This returns
    a tuple of the forwards gradients, backwards gradients and the value
    of delta lambda. Note that for pure TI calculations, the forwards and
    backwards gradients will be equal and the value of delta lambda will be 0 */
Gradients TI::operator[](int i) const
{
    i = Index(i).map(grads.count());
    return grads.at(i);
}

/** Return the free energy gradient data for the ith iteration. This returns
    a tuple of the forwards gradients, backwards gradients and the value
    of delta lambda. Note that for pure TI calculations, the forwards and
    backwards gradients will be equal and the value of delta lambda will be 0 */
Gradients TI::at(int i) const
{
    return operator[](i);
}

/** Set the gradients for the ith iteration equal to 'gradients'. These
    are analytic TI gradients */
void TI::set(int i, const QMap<double,AverageAndStddev> &gradients)
{
    while (i >= grads.count())
    {
        grads.append( Gradients() );
    }

    i = Index(i).map(grads.count());

    if (gradients.isEmpty())
        grads[i] = Gradients();
    else
        grads[i] = Gradients(gradients);
}

/** Set the gradients for the ith iteration equal to 'gradients'. These
    must be pure TI gradients, with no associated delta lambda value */
void TI::set(int i, const QMap<double,FreeEnergyAverage> &gradients)
{
    while (i >= grads.count())
    {
        grads.append( Gradients() );
    }

    i = Index(i).map(grads.count());

    if (gradients.isEmpty())
        grads[i] = Gradients();
    else
        grads[i] = Gradients(gradients);
}

/** Set the gradients for the ith iteration equal to 'gradients'. These
    are finite difference TI gradients, which are the raw zwanzig
    free energies together with the passed value of delta lambda */
void TI::set(int i, const QMap<double,FreeEnergyAverage> &gradients,
             double delta_lambda)
{
    while (i >= grads.count())
    {
        grads.append( Gradients() );
    }

    i = Index(i).map(grads.count());

    if (gradients.isEmpty())
        grads[i] = Gradients();
    else
        grads[i] = Gradients(gradients, delta_lambda);
}

/** Set the gradients for the ith iteration to the passed forwards
    and backwards finite difference TI values (together with associated
    delta lambda). The forwards gradients should be the raw zwanzig
    values from lambda -> lambda+delta_lambda, while the backwards
    gradients should be the raw zwanzig values from lambda-delta_lambda -> lambda */
void TI::set(int i, const QMap<double,FreeEnergyAverage> &forwards,
                    const QMap<double,FreeEnergyAverage> &backwards,
                    double delta_lambda)
{
    while (i >= grads.count())
    {
        grads.append( Gradients() );
    }

    i = Index(i).map(grads.count());

    if (forwards.isEmpty() and backwards.isEmpty())
        grads[i] = Gradients();
    else
        grads[i] = Gradients(forwards, backwards, delta_lambda);
}

void TI::set(int i, const Gradients &gradients)
{
    while (i >= grads.count())
    {
        grads.append( Gradients() );
    }

    i = Index(i).map(grads.count());

    if (gradients.isEmpty())
        grads[i] = Gradients();
    else
        grads[i] = gradients;
}

/** Remove the data for iteration 'i'. This sets the data equal to Gradients() */
void TI::removeAt(int i)
{
    i = Index(i).map(grads.count());
    grads[i] = Gradients();
}

/** Remove every iteration from 'start' to 'end' (inclusively). This sets
    the data equal to Gradients() */
void TI::removeRange(int start, int end)
{
    start = Index(start).map(grads.count());
    end = Index(end).map(grads.count());

    if (start > end)
        qSwap(start, end);

    for (int i = start; i <= end; ++i)
    {
        grads[i] = Gradients();
    }
}

/** Remove all values from the histogram */
void TI::clear()
{
    this->operator=( TI() );
}

/** Merge (average) together the gradients from iteration "start" to iteration
    "end" inclusive */
Gradients TI::merge(int start, int end) const
{
    start = Index(start).map(grads.count());
    end = Index(end).map(grads.count());

    QList<Gradients> set;

    for (int i=start; i<=end; ++i)
    {
        if (not grads.at(i).isEmpty())
            set.append( grads.at(i) );
    }

    return Gradients::merge(set);
}

/** Merge together the gradients from the iterations with the passed indicies */
Gradients TI::merge(QList<int> indicies) const
{
    QList<Gradients> set;

    foreach (int idx, indicies)
    {
        int i = Index(idx).map(grads.count());

        if (not grads.at(i).isEmpty())
            set.append( grads.at(i) );
    }

    return Gradients::merge(set);
}

/** Return a list of Gradients that represents the rolling average over 'niterations'
    iterations over this TI data set. If this data set contains 100 iterations, and
    we calculate the rolling average over 50 iterations, then the returned Gradients
    will be the average from 1-50, then 2-51, 3-52.....51-100 */
QList<Gradients> TI::rollingAverage(int niterations) const
{
    if (grads.isEmpty())
        return QList<Gradients>();

    QList<Gradients> merged;

    if (niterations >= grads.count())
    {
        Gradients m = this->merge(0,-1);

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
        QList<Gradients> set;

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
            merged.append( Gradients::merge(set) );

        for (i = i+1; i<grads.count(); ++i)
        {
            if (not grads.at(i).isEmpty())
            {
                set.removeFirst();
                set.append(grads.at(i));
                merged.append( Gradients::merge(set) );
            }
        }
    }

    return merged;
}

QString TI::toString() const
{
    return QObject::tr( "TI( nLambdaValues() == %1, nIterations() == %2, nSamples() == %3 )" )
                .arg(nLambdaValues()).arg(nIterations()).arg(nSamples());
}
