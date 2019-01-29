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

#include "fep.h"

#include "SireMaths/maths.h"

#include "SireID/index.h"

#include "SireError/errors.h"

#include "SireStream/shareddatastream.h"
#include "SireStream/registeralternativename.h"

#include "tostring.h"

using namespace SireAnalysis;
using namespace SireMaths;
using namespace SireBase;
using namespace SireID;
using namespace SireUnits::Dimension;
using namespace SireStream;


/////////
///////// Implementation of DataPoint
/////////

static const RegisterMetaType<DataPoint> r_dp(NO_ROOT);
static const RegisterAlternativeName<DataPoint> r_altdp("Soiree::DataPoint");

QDataStream &operator<<(QDataStream &ds, const DataPoint &dp)
{
    writeHeader(ds, r_dp, 1);
    
    ds << dp._x << dp._y
       << dp._xminerr << dp._yminerr
       << dp._xmaxerr << dp._ymaxerr;
    
    return ds;
}

QDataStream &operator>>(QDataStream &ds, DataPoint &dp)
{
    VersionID v = readHeader(ds, r_dp);
    
    if (v == 1)
    {
        ds >> dp._x >> dp._y
           >> dp._xminerr >> dp._yminerr
           >> dp._xmaxerr >> dp._ymaxerr;
    }
    else
        throw version_error(v, "1", r_dp, CODELOC);
    
    return ds;
}

/** Constructor. This constructs the point (0,0) with no error */
DataPoint::DataPoint()
          : _x(0), _y(0), _xminerr(0), _yminerr(0), _xmaxerr(0), _ymaxerr(0)
{}

/** Construct the point (x,y) with no error */
DataPoint::DataPoint(double x, double y)
          : _x(x), _y(y), _xminerr(0), _yminerr(0), _xmaxerr(0), _ymaxerr(0)
{}

/** Construct the point (x,y) with error (xerror,yerror) */
DataPoint::DataPoint(double x, double y, double xerror, double yerror)
          : _x(x), _y(y), _xminerr(std::abs(xerror)), _yminerr(std::abs(yerror)),
            _xmaxerr(std::abs(xerror)), _ymaxerr(std::abs(yerror))
{}

/** Construct the point (x,y) with a minimum error of (xminerror,yminerror)
    and a maximum error of (xmaxerror,ymaxerror). This is for situations where
    there may be multiple different error measures on a point and you want
    to store the range of errors (based on a range of error estimates) */
DataPoint::DataPoint(double x, double y, double xminerror, double yminerror,
                     double xmaxerror, double ymaxerror)
          : _x(x), _y(y),
            _xminerr(std::abs(xminerror)), _yminerr(std::abs(yminerror)),
            _xmaxerr(std::abs(xmaxerror)), _ymaxerr(std::abs(ymaxerror))
{
    if (_xminerr > _xmaxerr)
        qSwap(_xminerr, _xmaxerr);
    
    if (_yminerr > _ymaxerr)
        qSwap(_yminerr, _ymaxerr);
}

/** Copy constructor */
DataPoint::DataPoint(const DataPoint &other)
          : _x(other._x), _y(other._y),
            _xminerr(other._xminerr), _yminerr(other._yminerr),
            _xmaxerr(other._xmaxerr), _ymaxerr(other._ymaxerr)
{}

/** Destructor */
DataPoint::~DataPoint()
{}

/** Copy assignment operator */
DataPoint& DataPoint::operator=(const DataPoint &other)
{
    if (this != &other)
    {
        _x = other._x;
        _y = other._y;
        _xminerr = other._xminerr;
        _yminerr = other._yminerr;
        _xmaxerr = other._xmaxerr;
        _ymaxerr = other._ymaxerr;
    }
    
    return *this;
}

/** Comparison operator */
bool DataPoint::operator==(const DataPoint &other) const
{
    return _x == other._x and _y == other._y and
           _xminerr == other._xminerr and _yminerr == other._yminerr and
           _xmaxerr == other._xmaxerr and _ymaxerr == other._ymaxerr;
}

/** Comparison operator */
bool DataPoint::operator!=(const DataPoint &other) const
{
    return not operator==(other);
}

const char* DataPoint::what() const
{
    return DataPoint::typeName();
}

const char* DataPoint::typeName()
{
    return QMetaType::typeName( qMetaTypeId<DataPoint>() );
}

QString DataPoint::toString() const
{
    if (hasError())
    {
        QString xstr;
        QString ystr;
    
        if (hasXError())
        {
            if (xMinError() != xMaxError())
            {
                xstr = QString("%1 +/- %2 [%3]").arg(x()).arg(xMinError()).arg(xMaxError());
            }
            else
            {
                xstr = QString("%1 +/- %2").arg(x()).arg(xError());
            }
        }
        else
            xstr = QString("%1").arg(x());

        if (hasYError())
        {
            if (yMinError() != yMaxError())
            {
                ystr = QString("%1 +/- %2 [%3]").arg(y()).arg(yMinError()).arg(yMaxError());
            }
            else
            {
                ystr = QString("%1 +/- %2").arg(y()).arg(yError());
            }
        }
        else
            ystr = QString("%1").arg(y());
        
        return QString("DataPoint( %1, %2 )").arg(xstr, ystr);
    }
    else
        return QString("DataPoint( %1, %2 )").arg(x()).arg(y());
}

/** Return the x value of the point */
double DataPoint::x() const
{
    return _x;
}

/** Return the y value of the point */
double DataPoint::y() const
{
    return _y;
}

/** Return the error on the x value. This is the average
    of the minimum and maximum error */
double DataPoint::xError() const
{
    return 0.5 * (_xminerr + _xmaxerr);
}

/** Return the error on the y value. This is the average
    of the minimum and maximum error */
double DataPoint::yError() const
{
    return 0.5 * (_yminerr + _ymaxerr);
}

/** Return the minimum size of the error on the x value */
double DataPoint::xMinError() const
{
    return _xminerr;
}

/** Return the minimum size of the error on the y value */
double DataPoint::yMinError() const
{
    return _yminerr;
}

/** Return the maximum size of the error on the x value */
double DataPoint::xMaxError() const
{
    return _xmaxerr;
}

/** Return the maximum size of the error on the y value */
double DataPoint::yMaxError() const
{
    return _ymaxerr;
}

/** Return whether or not this data point has any error */
bool DataPoint::hasError() const
{
    return _xminerr > 0 or _yminerr > 0;
}

/** Return whether or not this data point has an error range */
bool DataPoint::hasErrorRange() const
{
    return _xminerr != _xmaxerr or _yminerr != _ymaxerr;
}

/** Return whether or not there is any error in the x value */
bool DataPoint::hasXError() const
{
    return _xminerr > 0;
}

/** Return whether or not there is an error range on the x value */
bool DataPoint::hasXErrorRange() const
{
    return _xminerr != _xmaxerr;
}

/** Return whether or not there is any error in the y value */
bool DataPoint::hasYError() const
{
    return _yminerr > 0;
}

/** Return whether or not there is an error range on the x value */
bool DataPoint::hasYErrorRange() const
{
    return _yminerr != _ymaxerr;
}

/** Return whether or not this data point is equal to the other, within
    the error range of the two points */
bool DataPoint::equalWithinError(const DataPoint &other) const
{
    return x() + xError() >= other.x() - other.xError() and
           x() - xError() <= other.x() + other.xError() and
    
           y() + yError() >= other.y() - other.yError() and
           y() - yError() <= other.y() + other.yError();
}

/** Return whether or not this data point in equal to the other, within
    the minimum error range of the two points */
bool DataPoint::equalWithinMinError(const DataPoint &other) const
{
    return x() + xMinError() >= other.x() - other.xMinError() and
           x() - xMinError() <= other.x() + other.xMinError() and
    
           y() + yMinError() >= other.y() - other.yMinError() and
           y() - yMinError() <= other.y() + other.yMinError();
}

/** Return whether or not this data point in equal to the other, within
    the maximum error range of the two points */
bool DataPoint::equalWithinMaxError(const DataPoint &other) const
{
    return x() + xMaxError() >= other.x() - other.xMaxError() and
           x() - xMaxError() <= other.x() + other.xMaxError() and
    
           y() + yMaxError() >= other.y() - other.yMaxError() and
           y() - yMaxError() <= other.y() + other.yMaxError();
}

/////////
///////// Implementation of PMF
/////////

static const RegisterMetaType<PMF> r_pmf;
static const RegisterAlternativeName<PMF> r_altpmf("Soiree::PMF");

QDataStream &operator<<(QDataStream &ds, const PMF &pmf)
{
    writeHeader(ds, r_pmf, 2);
    
    SharedDataStream sds(ds);
    
    sds << pmf.vals;
    
    return ds;
}

QDataStream &operator>>(QDataStream &ds, PMF &pmf)
{
    VersionID v = readHeader(ds, r_pmf);
    
    if (v == 2)
    {
        SharedDataStream sds(ds);
        
        sds >> pmf.vals;
    }
    else
        throw version_error(v, "2", r_pmf, CODELOC);
    
    return ds;
}

/** Construct an empty PMF */
PMF::PMF() : ConcreteProperty<PMF,Property>()
{}

void PMF::setValues(const QVector<DataPoint> &values)
{
    vals = values;
    
    //ensure that the points are in sorted x numerical order
    bool sorted = false;
    
    while (not sorted)
    {
        sorted = true;
    
        //bubble sort - compare neighbours and swap if in the wrong order
        for (int i=1; i<vals.count(); ++i)
        {
            if (vals[i-1].x() > vals[i].x())
            {
                qSwap(vals[i-1], vals[i]);
                sorted = false;
            }
        }
    }
}

/** Construct a from the passed data points */
PMF::PMF(const QVector<DataPoint> &values) : ConcreteProperty<PMF,Property>()
{
    setValues(values);
}

/** Copy constructor */
PMF::PMF(const PMF &other)
    : ConcreteProperty<PMF,Property>(other), vals(other.vals)
{}

/** Destructor */
PMF::~PMF()
{}

/** Copy assignment operator */
PMF& PMF::operator=(const PMF &other)
{
    if (this != &other)
    {
        vals = other.vals;
    }
    
    return *this;
}

/** Comparison operator */
bool PMF::operator==(const PMF &other) const
{
    return vals == other.vals;
}

/** Comparison operator */
bool PMF::operator!=(const PMF &other) const
{
    return not operator==(other);
}

const char* PMF::what() const
{
    return PMF::typeName();
}

const char* PMF::typeName()
{
    return QMetaType::typeName( qMetaTypeId<PMF>() );
}

/** Return the minimum x-value of the PMF */
double PMF::rangeMin() const
{
    if (vals.isEmpty())
        return 0;
    else
        return vals.first().x();
}

/** Return the maximum x-value of the PMF */
double PMF::rangeMax() const
{
    if (vals.isEmpty())
        return 0;
    else
        return vals.last().x();
}

/** Return the total free energy change along the PMF (difference in
    free energy of the end-points) */
double PMF::deltaG() const
{
    if (vals.isEmpty())
        return 0;
    else
        return vals.last().y() - vals.first().y();
}

/** Return the error on the total free energy calculation */
double PMF::error() const
{
    if (vals.isEmpty())
        return 0;

    else if (vals.count() == 1)
        return vals.at(0).yError();
    else
        return vals.last().yError() + vals.first().yError();
}

/** Return whether or not this PMF is empty (has not values) */
bool PMF::isEmpty() const
{
    return vals.isEmpty();
}

QString PMF::toString() const
{
    if (vals.isEmpty())
        return QString("PMF()");
    else
        return QString("PMF( { deltaG() == %1, error() == %2 } )")
                    .arg(deltaG()).arg(error());
}

/** Return the raw data for the PMF */
QVector<DataPoint> PMF::values() const
{
    return vals;
}

/////////
///////// Implementation of FEPDeltas
/////////

static const RegisterMetaType<FEPDeltas> r_deltas;
static const RegisterAlternativeName<FEPDeltas> r_altdeltas("Soiree::FEPDeltas");

QDataStream &operator<<(QDataStream &ds, const FEPDeltas &deltas)
{
    writeHeader(ds, r_deltas, 1);
    
    SharedDataStream sds(ds);
    
    sds << deltas.lamvals << deltas.fwds_deltas << deltas.bwds_deltas;
    
    return ds;
}

QDataStream &operator>>(QDataStream &ds, FEPDeltas &deltas)
{
    VersionID v = readHeader(ds, r_deltas);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> deltas.lamvals >> deltas.fwds_deltas >> deltas.bwds_deltas;
    }
    else
        throw version_error(v, "1", r_deltas, CODELOC);
    
    return ds;
}

void FEPDeltas::checkSane() const
{
    Temperature t;
    bool have_first = false;

    //make sure that there are no repeated lambda windows
    for (int i=0; i<lamvals.count()-1; ++i)
    {
        if (lamvals[i] == lamvals[i+1])
            throw SireError::invalid_arg( QObject::tr(
                    "You cannot have duplicate values of FEP windows. %1")
                        .arg(Sire::toString(lamvals)), CODELOC );
    }

    //check that all of the deltas match up with the supplied lambda values
    for (QMap<double,FreeEnergyAverage>::const_iterator it = fwds_deltas.constBegin();
         it != fwds_deltas.constEnd();
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
                "You cannot construct a set of FEP deltas using free energy "
                "averages collected at different temperatures. %1 vs. %2")
                    .arg(t.toString()).arg(it.value().temperature().toString()),
                        CODELOC );
        }
    
        int idx = lamvals.indexOf(it.key());
        
        if (idx == -1)
            throw SireError::invalid_arg( QObject::tr(
                    "All of the FEP deltas must correspond to one of the FEP windows. "
                    "The forwards delta with value %1 at window %2 is not in the list of "
                    "windows %3")
                        .arg(it.key()).arg(it.value().toString())
                        .arg(Sire::toString(lamvals)), CODELOC );
        
        if (idx == lamvals.count())
            //there should be no forwards delta for the last window
            throw SireError::invalid_arg( QObject::tr(
                    "There should be no forwards delta (%1) for the last FEP window (%2).")
                        .arg(it.value().toString()).arg(it.key()), CODELOC );
    }

    for (QMap<double,FreeEnergyAverage>::const_iterator it = bwds_deltas.constBegin();
         it != bwds_deltas.constEnd();
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
                "You cannot construct a set of FEP deltas using free energy "
                "averages collected at different temperatures. %1 vs. %2")
                    .arg(t.toString()).arg(it.value().temperature().toString()),
                        CODELOC );
        }

        int idx = lamvals.indexOf(it.key());
        
        if (idx == -1)
            throw SireError::invalid_arg( QObject::tr(
                    "All of the FEP deltas must correspond to one of the FEP windows. "
                    "The backwards delta with value %1 at window %2 is not in the list of "
                    "windows %3")
                        .arg(it.key()).arg(it.value().toString())
                        .arg(Sire::toString(lamvals)), CODELOC );
        
        if (idx == lamvals.count())
            //there should be no backwards delta for the first window
            throw SireError::invalid_arg( QObject::tr(
                    "There should be no backwards delta (%1) for the first FEP window (%2).")
                        .arg(it.value().toString()).arg(it.key()), CODELOC );
    }
}

/** Construct an empty set of deltas */
FEPDeltas::FEPDeltas() : ConcreteProperty<FEPDeltas,Property>()
{}

/** Construct the deltas as the deltas between each window and the window above */
FEPDeltas::FEPDeltas(const QList<double> &windows, const QMap<double,FreeEnergyAverage> &deltas)
          : ConcreteProperty<FEPDeltas,Property>(),
            lamvals(windows), fwds_deltas(deltas)
{
    qSort(lamvals);
    checkSane();
}

/** Construct the deltas as the deltas between each window and the windows above
    (forwards_deltas) and windows below (backwards deltas) */
FEPDeltas::FEPDeltas(const QList<double> &windows,
                     const QMap<double,FreeEnergyAverage> &forwards_deltas,
                     const QMap<double,FreeEnergyAverage> &backwards_deltas)
          : ConcreteProperty<FEPDeltas,Property>(),
            lamvals(windows), fwds_deltas(forwards_deltas), bwds_deltas(backwards_deltas)
{
    qSort(lamvals);
    checkSane();
}

/** Copy constructor */
FEPDeltas::FEPDeltas(const FEPDeltas &other)
          : ConcreteProperty<FEPDeltas,Property>(other),
            lamvals(other.lamvals), fwds_deltas(other.fwds_deltas), bwds_deltas(other.bwds_deltas)
{}

/** Destructor */
FEPDeltas::~FEPDeltas()
{}

/** Copy assignment operator */
FEPDeltas& FEPDeltas::operator=(const FEPDeltas &other)
{
    if (this != &other)
    {
        lamvals = other.lamvals;
        fwds_deltas = other.fwds_deltas;
        bwds_deltas = other.bwds_deltas;
    }
    
    return *this;
}

/** Comparison operator */
bool FEPDeltas::operator==(const FEPDeltas &other) const
{
    return lamvals == other.lamvals and
           fwds_deltas == other.fwds_deltas and
           bwds_deltas == other.bwds_deltas;
}

/** Comparison operator */
bool FEPDeltas::operator!=(const FEPDeltas &other) const
{
    return not operator==(other);
}

const char* FEPDeltas::what() const
{
    return FEPDeltas::typeName();
}

const char* FEPDeltas::typeName()
{
    return QMetaType::typeName( qMetaTypeId<FEPDeltas>() );
}

/** Return the temperature at which the FEP deltas were all collected */
Temperature FEPDeltas::temperature() const
{
    if (this->isEmpty())
        return Temperature(0);
    else
    {
        if (not fwds_deltas.isEmpty())
            return fwds_deltas.constBegin()->temperature();
        else
            return bwds_deltas.constBegin()->temperature();
    }
}

QString FEPDeltas::toString() const
{
    return QObject::tr("FEPDeltas( nWindows() == %1, nSamples() == %2, temperature() == %3 )")
                .arg(nWindows()).arg(nSamples()).arg(temperature().toString());
}

/** Return whether or not this is empty */
bool FEPDeltas::isEmpty() const
{
    return lamvals.isEmpty();
}

/** Self-addition operator */
FEPDeltas& FEPDeltas::operator+=(const FEPDeltas &other)
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
        if (lamvals != other.lamvals)
            throw SireError::incompatible_error( QObject::tr(
                "Cannot add together these two FEPDeltas as the lambda windows are different. "
                "%1 vs. %2.")
                    .arg(Sire::toString(lamvals)).arg(Sire::toString(other.lamvals)),
                        CODELOC);
        
        if (temperature() != other.temperature())
            throw SireError::incompatible_error( QObject::tr(
                "Cannot add together these two FEPDeltas as the temperature at which their "
                "free energies were collected are different. %1 vs. %2")
                    .arg(temperature().toString()).arg(other.temperature().toString()),
                        CODELOC );
        
        QMap<double,FreeEnergyAverage> new_fwds_deltas = fwds_deltas;
        QMap<double,FreeEnergyAverage> new_bwds_deltas = bwds_deltas;
        
        for (QMap<double,FreeEnergyAverage>::const_iterator it = other.fwds_deltas.constBegin();
             it != other.fwds_deltas.constEnd();
             ++it)
        {
            if (new_fwds_deltas.contains(it.key()))
                new_fwds_deltas[it.key()] += it.value();
            else
                new_fwds_deltas.insert(it.key(), it.value());
        }

        for (QMap<double,FreeEnergyAverage>::const_iterator it = other.bwds_deltas.constBegin();
             it != other.bwds_deltas.constEnd();
             ++it)
        {
            if (new_bwds_deltas.contains(it.key()))
                new_bwds_deltas[it.key()] += it.value();
            else
                new_bwds_deltas.insert(it.key(), it.value());
        }
        
        fwds_deltas = new_fwds_deltas;
        bwds_deltas = new_bwds_deltas;
        
        return *this;
    }
}

/** Addition operator */
FEPDeltas FEPDeltas::operator+(const FEPDeltas &other) const
{
    FEPDeltas ret(*this);
    ret += other;
    return *this;
}

/** Merge together all of the passed FEPDeltas into a single object */
FEPDeltas FEPDeltas::merge(const QList<FEPDeltas> &deltas)
{
    if (deltas.isEmpty())
        return FEPDeltas();
    
    else if (deltas.count() == 1)
        return deltas.at(0);
    
    else
    {
        FEPDeltas ret = deltas.at(0);

        for (int i=1; i<deltas.count(); ++i)
        {
            ret += deltas.at(i);
        }
        
        return ret;
    }
}

/** Return the lambda values for all of the windows */
QList<double> FEPDeltas::lambdaValues() const
{
    return lamvals;
}

/** Return the values of all of the windows */
QList<double> FEPDeltas::windows() const
{
    return lamvals;
}

/** Return the number of lambda values (windows) */
int FEPDeltas::nLambdaValues() const
{
    return lambdaValues().count();
}

/** Return the number of windows */
int FEPDeltas::nWindows() const
{
    return windows().count();
}

/** Return the total number of samples in the deltas */
qint64 FEPDeltas::nSamples() const
{
    quint64 n = 0;
    
    for (QMap<double,FreeEnergyAverage>::const_iterator it = fwds_deltas.constBegin();
         it != fwds_deltas.constEnd();
         ++it)
    {
        n += it.value().nSamples();
    }
    
    for (QMap<double,FreeEnergyAverage>::const_iterator it = bwds_deltas.constBegin();
         it != bwds_deltas.constEnd();
         ++it)
    {
        n += it.value().nSamples();
    }
    
    return n;
}

/** Return the values between windows. This returns the average of the 
    forwards and backwards values */
QVector<DataPoint> FEPDeltas::values() const
{
    QVector<DataPoint> points;
    
    for (int i=0; i<lamvals.count()-1; ++i)
    {
        const FreeEnergyAverage *fwds = 0;
        const FreeEnergyAverage *bwds = 0;
        
        double lam = lamvals.at(i);
        
        if (fwds_deltas.contains(lam))
            fwds = &(*(fwds_deltas.constFind(lam)));
        
        if (bwds_deltas.contains(lamvals.at(i+1)))
            bwds = &(*(bwds_deltas.constFind(lamvals.at(i+1))));
        
        if (fwds == 0)
        {
            double val = bwds->fepFreeEnergy();
            double minerr = bwds->histogram().standardError(90);
            double maxerr = minerr;
            
            points.append( DataPoint(lam, val, 0, minerr, 0, maxerr) );
        }
        else if (bwds == 0)
        {
            double val = fwds->fepFreeEnergy();
            double minerr = fwds->histogram().standardError(90);
            double maxerr = minerr;
            
            points.append( DataPoint(lam, val, 0, minerr, 0, maxerr) );
        }
        else
        {
            double val = 0.5 * (fwds->fepFreeEnergy() + bwds->fepFreeEnergy());
            double minerr = std::abs(fwds->fepFreeEnergy() - bwds->fepFreeEnergy());
            double maxerr = minerr + fwds->histogram().standardError(90)
                                   + bwds->histogram().standardError(90);
            
            points.append( DataPoint(lam, val, 0, minerr, 0, maxerr) );
        }
    }
    
    return points;
}

/** Return the forwards deltas. This returns the lambda value of the from window,
    together with the value of the free energy delta */
QVector<DataPoint> FEPDeltas::forwardsValues() const
{
    QVector<DataPoint> points;
    
    foreach (double lamval, lamvals)
    {
        if (fwds_deltas.contains(lamval))
        {
            const FreeEnergyAverage &fwds = *(fwds_deltas.constFind(lamval));
        
            double val = fwds.fepFreeEnergy();
            double maxerr = fwds.histogram().standardError(90);
            double minerr = maxerr;
            points.append( DataPoint(lamval, val, 0, minerr, 0, maxerr) );
        }
    }
    
    return points;
}

/** Return the backwards deltas. This returns the lambda value of the from window,
    together with the free energy delta */
QVector<DataPoint> FEPDeltas::backwardsValues() const
{
    QVector<DataPoint> points;
    
    for (int i=1; i<lamvals.count(); ++i)
    {
        if (bwds_deltas.contains(lamvals.at(i)))
        {
            const FreeEnergyAverage &bwds = *(bwds_deltas.constFind(lamvals.at(i)));
        
            double val = bwds.fepFreeEnergy();
            double maxerr = bwds.histogram().standardError(90);
            double minerr = maxerr;
            points.append( DataPoint(lamvals.at(i-1), val, 0, minerr, 0, maxerr) );
        }
    }
    
    return points;
}

/** Return the raw data for the fowards deltas */
QMap<double,FreeEnergyAverage> FEPDeltas::forwardsData() const
{
    return fwds_deltas;
}

/** Return the raw data for the backwards deltas */
QMap<double,FreeEnergyAverage> FEPDeltas::backwardsData() const
{
    return bwds_deltas;
}

/** Return the raw data for the fowards deltas */
QMap<double,FreeEnergyAverage> FEPDeltas::forwardsDeltas() const
{
    return fwds_deltas;
}

/** Return the raw data for the backwards deltas */
QMap<double,FreeEnergyAverage> FEPDeltas::backwardsDeltas() const
{
    return bwds_deltas;
}

/** Integrate (sum) the deltas across the windows to return the PMF */
PMF FEPDeltas::sum() const
{
    if (lamvals.isEmpty())
        return PMF();

    QVector<DataPoint> points;
    
    double total = 0;
    double total_minerr = 0;
    double total_maxerr = 0;
    
    points.append( DataPoint(lamvals.first(),0) );
    
    for (int i=0; i<lamvals.count()-1; ++i)
    {
        const FreeEnergyAverage *fwds = 0;
        const FreeEnergyAverage *bwds = 0;
        
        double lam = lamvals.at(i);
        
        double val = 0;
        double minerr = 0;
        double maxerr = 0;
        
        if (fwds_deltas.contains(lam))
            fwds = &(*(fwds_deltas.constFind(lam)));
        
        if (bwds_deltas.contains(lamvals.at(i+1)))
            bwds = &(*(bwds_deltas.constFind(lamvals.at(i+1))));
        
        if (fwds == 0)
        {
            val = bwds->fepFreeEnergy();
            minerr = bwds->histogram().standardError(90);
            maxerr = minerr;
        }
        else if (bwds == 0)
        {
            val = fwds->fepFreeEnergy();
            minerr = fwds->histogram().standardError(90);
            maxerr = minerr;
        }
        else
        {
            val = 0.5 * (fwds->fepFreeEnergy() + bwds->fepFreeEnergy());
            minerr = std::abs(fwds->fepFreeEnergy() - bwds->fepFreeEnergy());
            maxerr = minerr + fwds->histogram().standardError(90) +
                              bwds->histogram().standardError(90);
        }
        
        total += val;
        total_minerr += minerr;
        total_maxerr += maxerr;
        
        points.append( DataPoint(lamvals.at(i+1), total,
                                 0, total_minerr,
                                 0, total_maxerr ) );
    }
    
    return PMF(points);
}

/** Integrate (sum) the forwards deltas across the windows to return the PMF */
PMF FEPDeltas::sumForwards() const
{
    if (lamvals.isEmpty())
        return PMF();

    QVector<DataPoint> points;
    
    double total = 0;
    double total_err = 0;
    
    points.append( DataPoint(lamvals.first(),0) );
    
    for (int i=0; i<lamvals.count()-1; ++i)
    {
        const FreeEnergyAverage *fwds = 0;
        const FreeEnergyAverage *bwds = 0;
        
        double lam = lamvals.at(i);
        
        double val = 0;
        double err = 0;
        
        if (fwds_deltas.contains(lam))
            fwds = &(*(fwds_deltas.constFind(lam)));
        
        if (bwds_deltas.contains(lamvals.at(i+1)))
            bwds = &(*(bwds_deltas.constFind(lamvals.at(i+1))));
        
        if (fwds == 0)
        {
            err = bwds->histogram().standardError(90);
            val = bwds->average();
        }
        else
        {
            err = fwds->histogram().standardError(90);
            val = fwds->average();
        }
        
        total += val;
        total_err += err;
        
        points.append( DataPoint(lamvals.at(i+1), total, 0, total_err) );
    }
    
    return PMF(points);
}

/** Integrate (sum) the backwards deltas across the windows to return the PMF */
PMF FEPDeltas::sumBackwards() const
{
    if (lamvals.isEmpty())
        return PMF();

    QVector<DataPoint> points;
    
    double total = 0;
    double total_err = 0;
    
    points.append( DataPoint(lamvals.first(),0) );
    
    for (int i=0; i<lamvals.count()-1; ++i)
    {
        const FreeEnergyAverage *fwds = 0;
        const FreeEnergyAverage *bwds = 0;
        
        double lam = lamvals.at(i);
        
        double val = 0;
        double err = 0;
        
        if (fwds_deltas.contains(lam))
            fwds = &(*(fwds_deltas.constFind(lam)));
        
        if (bwds_deltas.contains(lamvals.at(i+1)))
            bwds = &(*(bwds_deltas.constFind(lamvals.at(i+1))));
        
        if (bwds == 0)
        {
            err = fwds->histogram().standardError(90);
            val = fwds->average();
        }
        else
        {
            err = bwds->histogram().standardError(90);
            val = bwds->average();
        }
        
        total += val;
        total_err += err;
        
        points.append( DataPoint(lamvals.at(i+1), total, 0, total_err) );
    }
    
    return PMF(points);
}

/** Integrate (sum) the forwards taylor expansions across the windows to return the PMF */
PMF FEPDeltas::sumForwardsTaylor() const
{
    if (lamvals.isEmpty())
        return PMF();

    QVector<DataPoint> points;
    
    double total = 0;
    double total_err = 0;
    
    points.append( DataPoint(lamvals.first(),0) );
    
    for (int i=0; i<lamvals.count()-1; ++i)
    {
        const FreeEnergyAverage *fwds = 0;
        const FreeEnergyAverage *bwds = 0;
        
        double lam = lamvals.at(i);
        
        double val = 0;
        double err = 0;
        
        if (fwds_deltas.contains(lam))
            fwds = &(*(fwds_deltas.constFind(lam)));
        
        if (bwds_deltas.contains(lamvals.at(i+1)))
            bwds = &(*(bwds_deltas.constFind(lamvals.at(i+1))));
        
        if (fwds == 0)
        {
            err = bwds->histogram().standardError(90);
            val = bwds->taylorExpansion();
        }
        else
        {
            err = fwds->histogram().standardError(90);
            val = fwds->taylorExpansion();
        }
        
        total += val;
        total_err += err;
        
        points.append( DataPoint(lamvals.at(i+1), total, 0, total_err) );
    }
    
    return PMF(points);
}

/** Integrate (sum) the backwards Taylor expansions across the windows to return the PMF */
PMF FEPDeltas::sumBackwardsTaylor() const
{
    if (lamvals.isEmpty())
        return PMF();

    QVector<DataPoint> points;
    
    double total = 0;
    double total_err = 0;
    
    points.append( DataPoint(lamvals.first(),0) );
    
    for (int i=0; i<lamvals.count()-1; ++i)
    {
        const FreeEnergyAverage *fwds = 0;
        const FreeEnergyAverage *bwds = 0;
        
        double lam = lamvals.at(i);
        
        double val = 0;
        double err = 0;
        
        if (fwds_deltas.contains(lam))
            fwds = &(*(fwds_deltas.constFind(lam)));
        
        if (bwds_deltas.contains(lamvals.at(i+1)))
            bwds = &(*(bwds_deltas.constFind(lamvals.at(i+1))));
        
        if (bwds == 0)
        {
            err = fwds->histogram().standardError(90);
            val = fwds->taylorExpansion();
        }
        else
        {
            err = bwds->histogram().standardError(90);
            val = bwds->taylorExpansion();
        }
        
        total += val;
        total_err += err;
        
        points.append( DataPoint(lamvals.at(i+1), total, 0, total_err) );
    }
    
    return PMF(points);
}

/** Integrate (sum) the deltas across the windows to return the PMF */
PMF FEPDeltas::integrate() const
{
    return sum();
}

//////////////
////////////// Implementation of FEP
//////////////

static const RegisterMetaType<FEP> r_fep;
static const RegisterAlternativeName<FEP> r_altfep("Soiree::FEP");

QDataStream &operator<<(QDataStream &ds, const FEP &fep)
{
    writeHeader(ds, r_fep, 1);
    
    SharedDataStream sds(ds);
    
    sds << fep.dltas;
    
    return ds;
}

QDataStream &operator>>(QDataStream &ds, FEP &fep)
{
    VersionID v = readHeader(ds, r_fep);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> fep.dltas;
    }
    else
        throw version_error(v, "1", r_fep, CODELOC);
    
    return ds;
}

/** Constructor */
FEP::FEP() : ConcreteProperty<FEP,Property>()
{}

/** Construct to use the passed set of windows, with the free energy deltas from
    each window to the window above */
FEP::FEP(const QList<double> &windows, const QMap<double,FreeEnergyAverage> &deltas)
    : ConcreteProperty<FEP,Property>()
{
    this->add( FEPDeltas(windows,deltas) );
}

/** Construct to use the passed windows, with the free energy deltas from 
    each window to the window above in 'forwards_deltas' and from the window
    below to each window in 'backwards_deltas' */
FEP::FEP(const QList<double> &windows,
         const QMap<double,FreeEnergyAverage> &forwards_deltas,
         const QMap<double,FreeEnergyAverage> &backwards_deltas)
    : ConcreteProperty<FEP,Property>()
{
    this->add( FEPDeltas(windows,forwards_deltas,backwards_deltas) );
}

/** Construct to use the passed FEP deltas */
FEP::FEP(const FEPDeltas &deltas)
    : ConcreteProperty<FEP,Property>()
{
    this->add(deltas);
}

/** Copy constructor */
FEP::FEP(const FEP &other) : ConcreteProperty<FEP,Property>(other), dltas(other.dltas)
{}

/** Destructor */
FEP::~FEP()
{}

/** Copy assignment operator */
FEP& FEP::operator=(const FEP &other)
{
    dltas = other.dltas;
    return *this;
}

/** Comparison operator */
bool FEP::operator==(const FEP &other) const
{
    return dltas == other.dltas;
}

/** Comparison operator */
bool FEP::operator!=(const FEP &other) const
{
    return not operator==(other);
}

const char* FEP::what() const
{
    return FEP::typeName();
}

const char* FEP::typeName()
{
    return QMetaType::typeName( qMetaTypeId<FEP>() );
}

QString FEP::toString() const
{
    return QObject::tr("FEP( nWindows() == %1, nIterations() == %2, nSamples() == %3 )")
                .arg(nWindows()).arg(nIterations()).arg(nSamples());
}

/** Add the data for the next iteration, which contains the deltas for the passed windows,
    with the free energy being for each window to the next window */
void FEP::add(const QList<double> &windows,
              const QMap<double,FreeEnergyAverage> &deltas)
{
    this->add( FEPDeltas(windows,deltas) );
}

/** Add the data for the next iteration, which contains the deltas for the passed windows,
    with forwards_deltas containing the free energy from each window to the next window,
    and backwards_deltas containing the free energy from the previous window to each window */
void FEP::add(const QList<double> &windows,
              const QMap<double,FreeEnergyAverage> &forwards_deltas,
              const QMap<double,FreeEnergyAverage> &backwards_deltas)
{
    this->add( FEPDeltas(windows,forwards_deltas,backwards_deltas) );
}

/** Add the data for the next iteration */
void FEP::add(const FEPDeltas &deltas)
{
    if (not deltas.isEmpty())
        dltas.append(deltas);
}

/** Return the number of iterations */
int FEP::nIterations() const
{
    return dltas.count();
}

/** Return the number of windows */
int FEP::nWindows() const
{
    return windows().count();
}

/** Return the number of lambda values (windows) */
int FEP::nLambdaValues() const
{
    return nWindows();
}

/** Return the total number of samples in the simulation */
qint64 FEP::nSamples() const
{
    quint64 n = 0;
    
    foreach (const FEPDeltas &delta, dltas)
    {
        n += delta.nSamples();
    }
    
    return n;
}

/** Return the number of iterations */
int FEP::count() const
{
    return dltas.count();
}

/** Return the number of iterations */
int FEP::size() const
{
    return dltas.count();
}

/** Return the values of all windows */
QList<double> FEP::lambdaValues() const
{
    return windows();
}

/** Return the value of all windows */
QList<double> FEP::windows() const
{
    QMap<double,int> vals;
    
    foreach (const FEPDeltas &delta, dltas)
    {
        foreach (double window, delta.windows())
        {
            vals.insert(window,1);
        }
    }
    
    QList<double> wdows = vals.keys();
    qSort(wdows);
    return wdows;
}

/** Return the deltas for the ith iteration */
FEPDeltas FEP::operator[](int i) const
{
    return dltas.at( Index(i).map(dltas.count()) );
}

/** Return the deltas for the ith iteration */
FEPDeltas FEP::at(int i) const
{
    return operator[](i);
}

/** Return the deltas for all iterations */
QList<FEPDeltas> FEP::deltas() const
{
    return dltas;
}

/** Set the deltas for the ith iteration */
void FEP::set(int i, const QList<double> &windows,
              const QMap<double,FreeEnergyAverage> &deltas)
{
    set(i, FEPDeltas(windows,deltas));
}

/** Set the deltas for the ith iteration */
void FEP::set(int i, const QList<double> &windows,
              const QMap<double,FreeEnergyAverage> &forwards_deltas,
              const QMap<double,FreeEnergyAverage> &backwards_deltas)
{
    set(i, FEPDeltas(windows,forwards_deltas,backwards_deltas));
}

/** Set the deltas for the ith iteration */
void FEP::set(int i, const FEPDeltas &deltas)
{
    while (i >= dltas.count())
    {
        dltas.append( FEPDeltas() );
    }

    i = Index(i).map(dltas.count());

    if (deltas.isEmpty())
        dltas[i] = FEPDeltas();
    else
        dltas[i] = deltas;
}

/** Merge the deltas for iterations start->end */
FEPDeltas FEP::merge(int start, int end) const
{
    start = Index(start).map(dltas.count());
    end = Index(end).map(dltas.count());
    
    QList<FEPDeltas> set;
    
    for (int i=start; i<=end; ++i)
    {
        if (not dltas.at(i).isEmpty())
            set.append( dltas.at(i) );
    }
    
    return FEPDeltas::merge(set);
}

/** Merge the deltas at the passed indicies */
FEPDeltas FEP::merge(QList<int> indicies) const
{
    QList<FEPDeltas> set;
    
    foreach (int idx, indicies)
    {
        int i = Index(idx).map(dltas.count());

        if (not dltas.at(i).isEmpty())
           set.append( dltas.at(i) );
    }
 
    return FEPDeltas::merge(set);
}

/** Return a list of Gradients that represents the rolling average over 'niterations'
    iterations over this TI data set. If this data set contains 100 iterations, and 
    we calculate the rolling average over 50 iterations, then the returned Gradients
    will be the average from 1-50, then 2-51, 3-52.....51-100 */
QList<FEPDeltas> FEP::rollingAverage(int niterations) const
{
    QList<FEPDeltas> merged;

    if (dltas.isEmpty())
        return merged;
    
    if (niterations >= dltas.count())
    {
        FEPDeltas d = this->merge(0,-1);
        if (not d.isEmpty())
            merged.append(d);
    }
    else if (niterations <= 1)
    {
        for (int i=0; i<dltas.count(); ++i)
        {
            if (not dltas.at(i).isEmpty())
                merged.append(dltas.at(i));
        }
    }
    else
    {
        QList<FEPDeltas> set;
        
        int i=0;
        
        for (i=0; i<dltas.count(); ++i)
        {
            if (not dltas.at(i).isEmpty())
            {
                set.append(dltas.at(i));
                if (set.count() == niterations)
                    break;
            }
        }
        
        merged.append( FEPDeltas::merge(set) );
        
        for (i=i+1; i<dltas.count(); ++i)
        {
            if (not dltas.at(i).isEmpty())
            {
                set.removeFirst();
                set.append(dltas.at(i));
                merged.append( FEPDeltas::merge(set) );
            }
        }
    }
    
    return merged;
}

/** Remove the data for iteration 'i' */
void FEP::removeAt(int i)
{
    i = Index(i).map(dltas.count());
    dltas[i] = FEPDeltas();
}

/** Remove every iteration from 'start' to 'end' (inclusively) */
void FEP::removeRange(int start, int end)
{
    start = Index(start).map(dltas.count());
    end = Index(end).map(dltas.count());
    
    if (start > end)
        qSwap(start, end);
    
    for (int i = start; i <= end; ++i)
    {
        dltas[i] = FEPDeltas();
    }
}

/** Remove all values from the histogram */
void FEP::clear()
{
    this->operator=( FEP() );
}
