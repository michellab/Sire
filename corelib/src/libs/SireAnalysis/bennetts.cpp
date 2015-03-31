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

#include "bennetts.h"

#include "SireMaths/maths.h"

#include "SireID/index.h"

#include "SireError/errors.h"

#include "SireUnits/units.h"
#include "SireUnits/temperature.h"

#include "SireStream/shareddatastream.h"
#include "SireStream/registeralternativename.h"

#include "tostring.h"

using namespace SireAnalysis;
using namespace SireMaths;
using namespace SireBase;
using namespace SireID;
using namespace SireUnits;
using namespace SireUnits::Dimension;
using namespace SireStream;

///////////
/////////// Implementation of BennettsRatios
///////////

static const RegisterMetaType<BennettsRatios> r_ratios;
static const RegisterAlternativeName<BennettsRatios> r_altratios("Soiree::BennettsRatios");

QDataStream SIREANALYSIS_EXPORT &operator<<(QDataStream &ds, const BennettsRatios &ratios)
{
    writeHeader(ds, r_ratios, 1);
    
    SharedDataStream sds(ds);
    
    sds << ratios.lamvals << ratios.fwds_ratios << ratios.bwds_ratios;
    
    return ds;
}

QDataStream SIREANALYSIS_EXPORT &operator>>(QDataStream &ds, BennettsRatios &ratios)
{
    VersionID v = readHeader(ds, r_ratios);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> ratios.lamvals >> ratios.fwds_ratios >> ratios.bwds_ratios;
    }
    else
        throw version_error(v, "1", r_ratios, CODELOC);
    
    return ds;
}

void BennettsRatios::checkSane() const
{
    Temperature t;
    bool have_first = false;

    //make sure that there are no repeated lambda windows
    for (int i=0; i<lamvals.count()-1; ++i)
    {
        if (lamvals[i] == lamvals[i+1])
            throw SireError::invalid_arg( QObject::tr(
                    "You cannot have duplicate values of Bennetts windows. %1")
                        .arg(Sire::toString(lamvals)), CODELOC );
    }

    //check that all of the deltas match up with the supplied lambda values
    for (QMap<double,BennettsFreeEnergyAverage>::const_iterator it = fwds_ratios.constBegin();
         it != fwds_ratios.constEnd();
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
                "You cannot construct a set of Bennetts acceptance ratios using free energy "
                "ratios collected at different temperatures. %1 vs. %2")
                    .arg(t.toString()).arg(it.value().temperature().toString()),
                        CODELOC );
        }
    
        int idx = lamvals.indexOf(it.key());
        
        if (idx == -1)
            throw SireError::invalid_arg( QObject::tr(
                    "All of the Bennetts ratios must correspond to one of the Bennetts windows. "
                    "The forwards ratio with value %1 at window %2 is not in the list of "
                    "windows %3")
                        .arg(it.key()).arg(it.value().toString())
                        .arg(Sire::toString(lamvals)), CODELOC );
        
        if (idx == lamvals.count())
            //there should be no forwards delta for the last window
            throw SireError::invalid_arg( QObject::tr(
                    "There should be no forwards ratio (%1) for the last Bennetts window (%2).")
                        .arg(it.value().toString()).arg(it.key()), CODELOC );
    }

    for (QMap<double,BennettsFreeEnergyAverage>::const_iterator it = bwds_ratios.constBegin();
         it != bwds_ratios.constEnd();
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
                "You cannot construct a set of Bennetts acceptance ratios using free energy "
                "ratios collected at different temperatures. %1 vs. %2")
                    .arg(t.toString()).arg(it.value().temperature().toString()),
                        CODELOC );
        }

        int idx = lamvals.indexOf(it.key());
        
        if (idx == -1)
            throw SireError::invalid_arg( QObject::tr(
                    "All of the Bennetts ratios must correspond to one of the Bennetts windows. "
                    "The backwards ratio with value %1 at window %2 is not in the list of "
                    "windows %3")
                        .arg(it.key()).arg(it.value().toString())
                        .arg(Sire::toString(lamvals)), CODELOC );
        
        if (idx == lamvals.count())
            //there should be no backwards delta for the first window
            throw SireError::invalid_arg( QObject::tr(
                    "There should be no backwards ratio (%1) for the first Bennetts window (%2).")
                        .arg(it.value().toString()).arg(it.key()), CODELOC );
    }
}

/** Construct an empty set of deltas */
BennettsRatios::BennettsRatios() : ConcreteProperty<BennettsRatios,Property>()
{}

/** Construct the ratios as the ratios between each window and the windows above
    (forwards_ratios) and windows below (backwards_ratios) */
BennettsRatios::BennettsRatios(const QList<double> &windows,
                               const QMap<double,BennettsFreeEnergyAverage> &forwards_ratios,
                               const QMap<double,BennettsFreeEnergyAverage> &backwards_ratios)
          : ConcreteProperty<BennettsRatios,Property>(),
            lamvals(windows), fwds_ratios(forwards_ratios), bwds_ratios(backwards_ratios)
{
    qSort(lamvals);
    checkSane();
}

/** Copy constructor */
BennettsRatios::BennettsRatios(const BennettsRatios &other)
          : ConcreteProperty<BennettsRatios,Property>(other),
            lamvals(other.lamvals), fwds_ratios(other.fwds_ratios), bwds_ratios(other.bwds_ratios)
{}

/** Destructor */
BennettsRatios::~BennettsRatios()
{}

/** Copy assignment operator */
BennettsRatios& BennettsRatios::operator=(const BennettsRatios &other)
{
    if (this != &other)
    {
        lamvals = other.lamvals;
        fwds_ratios = other.fwds_ratios;
        bwds_ratios = other.bwds_ratios;
    }
    
    return *this;
}

/** Comparison operator */
bool BennettsRatios::operator==(const BennettsRatios &other) const
{
    return lamvals == other.lamvals and
           fwds_ratios == other.fwds_ratios and
           bwds_ratios == other.bwds_ratios;
}

/** Comparison operator */
bool BennettsRatios::operator!=(const BennettsRatios &other) const
{
    return not operator==(other);
}

const char* BennettsRatios::what() const
{
    return BennettsRatios::typeName();
}

const char* BennettsRatios::typeName()
{
    return QMetaType::typeName( qMetaTypeId<BennettsRatios>() );
}

/** Return the temperature at which the Bennetts deltas were all collected */
Temperature BennettsRatios::temperature() const
{
    if (this->isEmpty())
        return Temperature(0);
    else
    {
        if (not fwds_ratios.isEmpty())
            return fwds_ratios.constBegin()->temperature();
        else
            return bwds_ratios.constBegin()->temperature();
    }
}

QString BennettsRatios::toString() const
{
    return QObject::tr("BennettsRatios( nWindows() == %1, nSamples() == %2, temperature() == %3 )")
                .arg(nWindows()).arg(nSamples()).arg(temperature().toString());
}

/** Return whether or not this is empty */
bool BennettsRatios::isEmpty() const
{
    return lamvals.isEmpty();
}

/** Self-addition operator */
BennettsRatios& BennettsRatios::operator+=(const BennettsRatios &other)
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
                "Cannot add together these two BennettsRatios as the lambda windows are different. "
                "%1 vs. %2.")
                    .arg(Sire::toString(lamvals)).arg(Sire::toString(other.lamvals)),
                        CODELOC);
        
        if (temperature() != other.temperature())
            throw SireError::incompatible_error( QObject::tr(
                "Cannot add together these two BennettsRatios as the temperature at which they "
                "were collected are different. %1 vs. %2")
                    .arg(temperature().toString()).arg(other.temperature().toString()),
                        CODELOC );
        
        QMap<double,BennettsFreeEnergyAverage> new_fwds_ratios = fwds_ratios;
        QMap<double,BennettsFreeEnergyAverage> new_bwds_ratios = bwds_ratios;
        
        for (QMap<double,BennettsFreeEnergyAverage>::const_iterator
                                    it = other.fwds_ratios.constBegin();
             it != other.fwds_ratios.constEnd();
             ++it)
        {
            if (new_fwds_ratios.contains(it.key()))
                new_fwds_ratios[it.key()] += it.value();
            else
                new_fwds_ratios.insert(it.key(), it.value());
        }

        for (QMap<double,BennettsFreeEnergyAverage>::const_iterator
                                    it = other.bwds_ratios.constBegin();
             it != other.bwds_ratios.constEnd();
             ++it)
        {
            if (new_bwds_ratios.contains(it.key()))
                new_bwds_ratios[it.key()] += it.value();
            else
                new_bwds_ratios.insert(it.key(), it.value());
        }
        
        fwds_ratios = new_fwds_ratios;
        bwds_ratios = new_bwds_ratios;
        
        return *this;
    }
}

/** Addition operator */
BennettsRatios BennettsRatios::operator+(const BennettsRatios &other) const
{
    BennettsRatios ret(*this);
    ret += other;
    return *this;
}

/** Merge together all of the passed BennettsRatios into a single object */
BennettsRatios BennettsRatios::merge(const QList<BennettsRatios> &deltas)
{
    if (deltas.isEmpty())
        return BennettsRatios();
    
    else if (deltas.count() == 1)
        return deltas.at(0);
    
    else
    {
        BennettsRatios ret = deltas.at(0);

        for (int i=1; i<deltas.count(); ++i)
        {
            ret += deltas.at(i);
        }
        
        return ret;
    }
}

/** Return the lambda values for all of the windows */
QList<double> BennettsRatios::lambdaValues() const
{
    return lamvals;
}

/** Return the values of all of the windows */
QList<double> BennettsRatios::windows() const
{
    return lamvals;
}

/** Return the number of lambda values (windows) */
int BennettsRatios::nLambdaValues() const
{
    return lambdaValues().count();
}

/** Return the number of windows */
int BennettsRatios::nWindows() const
{
    return windows().count();
}

/** Return the total number of samples in the deltas */
qint64 BennettsRatios::nSamples() const
{
    quint64 n = 0;
    
    for (QMap<double,BennettsFreeEnergyAverage>::const_iterator it = fwds_ratios.constBegin();
         it != fwds_ratios.constEnd();
         ++it)
    {
        n += it.value().nSamples();
    }
    
    for (QMap<double,BennettsFreeEnergyAverage>::const_iterator it = bwds_ratios.constBegin();
         it != bwds_ratios.constEnd();
         ++it)
    {
        n += it.value().nSamples();
    }
    
    return n;
}

const DataPoint& getPoint(const QVector<DataPoint> &points, double lam, bool *found)
{
    foreach (const DataPoint &point, points)
    {
        if (point.x() == lam)
        {
            *found = true;
            return point;
        }
    }
    
    static DataPoint empty;
    *found = false;
    return empty;
}

/** Return the values between windows. This returns the value of lambda of the from
    window, and the difference in free energy between this and the next window */
QVector<DataPoint> BennettsRatios::values() const
{
    QVector<DataPoint> consts = constants();
    QVector<DataPoint> nums = numerators();
    QVector<DataPoint> denoms = denominators();

    QVector<DataPoint> vals;

    for (int i=0; i<lamvals.count(); ++i)
    {
        double lamval = lamvals.at(i);
        
        bool found_const;
        bool found_num;
        bool found_denom;
        
        const DataPoint &c = getPoint(consts, lamval, &found_const);
        const DataPoint &num = getPoint(nums, lamval, &found_num);
        const DataPoint &denom = getPoint(denoms, lamval, &found_denom);
        
        if (found_const and found_num and found_denom)
        {
            //we have found a matched pair - calculate the ratio, and minimum
            //and maximum values
            double val = num.y() / denom.y();
            double minerr = (num.y()+num.yMinError())/(denom.y()-denom.yMinError());
            double maxerr = (num.y()+num.yMaxError())/(denom.y()-denom.yMaxError());
            
            // ratio = e^(-beta (dG - C)) so dG = -(1/beta) ln(ratio) + C
            
            val = -(k_boltz * temperature().to(kelvin) * std::log(val)) + c.y();
            minerr = -(k_boltz * temperature().to(kelvin) * std::log(minerr)) + c.y();
            maxerr = -(k_boltz * temperature().to(kelvin) * std::log(maxerr)) + c.y();
            
            vals.append( DataPoint(num.x(),val, 0, std::abs(val-minerr),
                                                0, std::abs(val-maxerr)) );
        }
        else if (found_num and found_denom)
        {
            //there are forwards and backwards FEP values that we can use to get dG
            //for this window
            if (not fwds_ratios.contains(lamval))
                throw SireError::program_bug( QObject::tr(
                        "No forwards value for lambda %1? %2")
                            .arg(lamval).arg(Sire::toString(fwds_ratios.keys())),
                                CODELOC );

            if (not bwds_ratios.contains( lamvals.at(i+1) ))
                throw SireError::program_bug( QObject::tr(
                        "No backwards value for lambda %1? %2")
                            .arg(lamvals.at(i+1)).arg(Sire::toString(bwds_ratios.keys())),
                                CODELOC );

            const BennettsFreeEnergyAverage &fwds = *(fwds_ratios.constFind(lamval));
            const BennettsFreeEnergyAverage &bwds = *(bwds_ratios.constFind(lamvals.at(i+1)));

            double fwdsval = fwds.fepFreeEnergy();
            double bwdsval = bwds.fepFreeEnergy();
            double val = 0.5 * (fwdsval + bwdsval);

            double minerr = 0.5 * std::abs(fwdsval - bwdsval);
            double maxerr = minerr + fwds.histogram().standardError(90) +
                                     bwds.histogram().standardError(90);
            
            vals.append( DataPoint(lamval,val, 0,minerr, 0,maxerr) );
        }
        else if (found_num)
        {
            //there is a forwards FEP value available from this lambda value to the next
            if (not fwds_ratios.contains(lamval))
                throw SireError::program_bug( QObject::tr(
                        "No forwards value for lambda %1? %2")
                            .arg(lamval).arg(Sire::toString(fwds_ratios.keys())),
                                CODELOC );

            const BennettsFreeEnergyAverage &fwds = *(fwds_ratios.constFind(lamval));
        
            double val = fwds.average();
            
            double err = fwds.histogram().standardError(90);
            
            vals.append( DataPoint(lamval, val, 0, err) );
        }
        else if (found_denom)
        {
            //there is a backwards FEP value available from the lambda value
            //above back down to this lambda value
            if (not bwds_ratios.contains( lamvals.at(i+1) ))
                throw SireError::program_bug( QObject::tr(
                        "No backwards value for lambda %1? %2")
                            .arg(lamvals.at(i+1)).arg(Sire::toString(bwds_ratios.keys())),
                                CODELOC );
            
            const FreeEnergyAverage &bwds = *(bwds_ratios.constFind(lamvals.at(i+1)));
        
            double val = bwds.average();
            double err = bwds.histogram().standardError(90);
            
            vals.append( DataPoint(lamval, val, 0, err) );
        }
        else
        {
            //no value for this lambda window is available
            vals.append( DataPoint(lamval, 0) );
        }
    }
    
    return vals;
}

/** Return the constants for each set of Bennetts acceptance ratios. This
    returns the lambda value of the from window, together with the constant
    used for the numerator from this window to the next window (which must
    be the same as the constant used for the denominator for the next window
    back to this window) */
QVector<DataPoint> BennettsRatios::constants() const
{
    QVector<DataPoint> points;
    
    for (int i=1; i<lamvals.count(); ++i)
    {
        if (fwds_ratios.contains(lamvals[i-1]) and bwds_ratios.contains(lamvals[i]))
        {
            const BennettsFreeEnergyAverage &fwds = *(fwds_ratios.constFind(lamvals[i-1]));
            const BennettsFreeEnergyAverage &bwds = *(bwds_ratios.constFind(lamvals[i]));
            
            if (fwds.constant() == bwds.constant())
                points.append( DataPoint(lamvals[i-1],fwds.constant().value()) );
            else
                qDebug() << "WARNING: Bennetts constants for numerator and denominator "
                         << "don't match!" << fwds.toString() << "vs." << bwds.toString();
            
        }
    }
    
    return points;
}

/** Return the numerators for the Bennetts acceptance ratio. This returns the 
    lambda value of the from window, together with the Bennetts ratio for the
    energy difference from this window to the next window */
QVector<DataPoint> BennettsRatios::numerators() const
{
    QVector<DataPoint> points;
    
    foreach (double lamval, lamvals)
    {
        if (fwds_ratios.contains(lamval))
        {
            const BennettsFreeEnergyAverage &fwds = *(fwds_ratios.constFind(lamval));
        
            double val = fwds.bennettsRatio();
            double minerr = fwds.bennettsStandardError(50);
            double maxerr = fwds.bennettsStandardError(95);
            
            if (maxerr < minerr)
                qSwap(maxerr, minerr);
            
            points.append( DataPoint(lamval, val, 0, minerr, 0, maxerr) );
        }
    }
    
    return points;
}

/** Return the denominators for the Bennetts acceptance ratio. This returns the 
    lambda value of the previous window, together with the Bennetts ratio for the
    energy difference from this window to the previous window */
QVector<DataPoint> BennettsRatios::denominators() const
{
    QVector<DataPoint> points;
    
    for (int i=1; i<lamvals.count(); ++i)
    {
        if (bwds_ratios.contains(lamvals[i]))
        {
            const BennettsFreeEnergyAverage &bwds = *(bwds_ratios.constFind(lamvals[i]));
        
            double val = bwds.bennettsRatio();
            double minerr = bwds.bennettsStandardError(50);
            double maxerr = bwds.bennettsStandardError(95);
            
            if (maxerr < minerr)
                qSwap(maxerr, minerr);
            
            points.append( DataPoint(lamvals[i-1], val, 0, minerr, 0, maxerr) );
        }
    }
    
    return points;
}

/** Return the raw data for the fowards ratios */
QMap<double,BennettsFreeEnergyAverage> BennettsRatios::forwardsData() const
{
    return fwds_ratios;
}

/** Return the raw data for the backwards ratios */
QMap<double,BennettsFreeEnergyAverage> BennettsRatios::backwardsData() const
{
    return bwds_ratios;
}

/** Return the raw data for the fowards ratios */
QMap<double,BennettsFreeEnergyAverage> BennettsRatios::forwardsRatios() const
{
    return fwds_ratios;
}

/** Return the raw data for the backwards ratios */
QMap<double,BennettsFreeEnergyAverage> BennettsRatios::backwardsRatios() const
{
    return bwds_ratios;
}

/** Integrate (sum) the deltas across the windows to return the PMF */
PMF BennettsRatios::sum() const
{
    if (lamvals.isEmpty())
        return PMF();

    QVector<DataPoint> vals = this->values();
    
    double total = 0;
    double total_minerr = 0;
    double total_maxerr = 0;
    
    QVector<DataPoint> points;
    points.append( DataPoint(lamvals.first(),0) );
    
    for (int i=0; i<lamvals.count()-1; ++i)
    {
        foreach (const DataPoint &val, vals)
        {
            if (val.x() == lamvals[i])
            {
                total += val.y();
                total_minerr += val.yMinError();
                total_maxerr += val.yMaxError();
                
                break;
            }
        }

        points.append( DataPoint(lamvals[i+1], total, 0, total_minerr, 0, total_maxerr) );
    }
    
    return PMF(points);
}

/** Integrate (sum) the deltas across the windows to return the PMF */
PMF BennettsRatios::integrate() const
{
    return sum();
}

///////////
/////////// Implementation of Bennetts
///////////

static const RegisterMetaType<Bennetts> r_bennetts;
static const RegisterAlternativeName<Bennetts> r_altbennets("Soiree::Bennetts");

QDataStream SIREANALYSIS_EXPORT &operator<<(QDataStream &ds, const Bennetts &bennetts)
{
    writeHeader(ds, r_bennetts, 1);
    
    SharedDataStream sds(ds);
    
    sds << bennetts.rtios;
    
    return ds;
}

QDataStream SIREANALYSIS_EXPORT &operator>>(QDataStream &ds, Bennetts &bennetts)
{
    VersionID v = readHeader(ds, r_bennetts);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> bennetts.rtios;
    }
    else
        throw version_error(v, "1", r_bennetts, CODELOC);
    
    return ds;
}

/** Constructor */
Bennetts::Bennetts() : ConcreteProperty<Bennetts,Property>()
{}

/** Construct to use the passed windows, with the free energy ratios from
    each window to the window above in 'forwards_ratios' and from the window
    below to each window in 'backwards_ratios' */
Bennetts::Bennetts(const QList<double> &windows,
                   const QMap<double,BennettsFreeEnergyAverage> &forwards_ratios,
                   const QMap<double,BennettsFreeEnergyAverage> &backwards_ratios)
    : ConcreteProperty<Bennetts,Property>()
{
    this->add( BennettsRatios(windows,forwards_ratios,backwards_ratios) );
}

/** Construct to use the passed Bennetts ratios */
Bennetts::Bennetts(const BennettsRatios &ratios)
         : ConcreteProperty<Bennetts,Property>()
{
    this->add(ratios);
}

/** Copy constructor */
Bennetts::Bennetts(const Bennetts &other)
         : ConcreteProperty<Bennetts,Property>(other), rtios(other.rtios)
{}

/** Destructor */
Bennetts::~Bennetts()
{}

/** Copy assignment operator */
Bennetts& Bennetts::operator=(const Bennetts &other)
{
    rtios = other.rtios;
    return *this;
}

/** Comparison operator */
bool Bennetts::operator==(const Bennetts &other) const
{
    return rtios == other.rtios;
}

/** Comparison operator */
bool Bennetts::operator!=(const Bennetts &other) const
{
    return not operator==(other);
}

const char* Bennetts::what() const
{
    return Bennetts::typeName();
}

const char* Bennetts::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Bennetts>() );
}

QString Bennetts::toString() const
{
    return QObject::tr("Bennetts( nWindows() == %1, nIterations() == %2, nSamples() == %3 )")
                .arg(nWindows()).arg(nIterations()).arg(nSamples());
}

/** Add the data for the next iteration, which contains the ratios for the passed windows,
    with forwards_ratios containing the free energy from each window to the next window,
    and backwards_ratios containing the free energy from the previous window to each window */
void Bennetts::add(const QList<double> &windows,
                   const QMap<double,BennettsFreeEnergyAverage> &forwards_ratios,
                   const QMap<double,BennettsFreeEnergyAverage> &backwards_ratios)
{
    this->add( BennettsRatios(windows,forwards_ratios,backwards_ratios) );
}

/** Add the data for the next iteration */
void Bennetts::add(const BennettsRatios &ratios)
{
    if (not ratios.isEmpty())
        rtios.append(ratios);
}

/** Return the number of iterations */
int Bennetts::nIterations() const
{
    return rtios.count();
}

/** Return the number of windows */
int Bennetts::nWindows() const
{
    return windows().count();
}

/** Return the number of lambda values (windows) */
int Bennetts::nLambdaValues() const
{
    return nWindows();
}

/** Return the total number of samples in the simulation */
qint64 Bennetts::nSamples() const
{
    quint64 n = 0;
    
    foreach (const BennettsRatios &ratio, rtios)
    {
        n += ratio.nSamples();
    }
    
    return n;
}

/** Return the number of iterations */
int Bennetts::count() const
{
    return rtios.count();
}

/** Return the number of iterations */
int Bennetts::size() const
{
    return rtios.count();
}

/** Return the values of all windows */
QList<double> Bennetts::lambdaValues() const
{
    return windows();
}

/** Return the value of all windows */
QList<double> Bennetts::windows() const
{
    QMap<double,int> vals;
    
    foreach (const BennettsRatios &ratio, rtios)
    {
        foreach (double window, ratio.windows())
        {
            vals.insert(window,1);
        }
    }
    
    QList<double> wdows = vals.keys();
    qSort(wdows);
    return wdows;
}

/** Return the deltas for the ith iteration */
BennettsRatios Bennetts::operator[](int i) const
{
    return rtios.at( Index(i).map(rtios.count()) );
}

/** Return the deltas for the ith iteration */
BennettsRatios Bennetts::at(int i) const
{
    return operator[](i);
}

/** Return the deltas for all iterations */
QList<BennettsRatios> Bennetts::ratios() const
{
    return rtios;
}

/** Set the deltas for the ith iteration */
void Bennetts::set(int i, const QList<double> &windows,
                   const QMap<double,BennettsFreeEnergyAverage> &forwards_ratios,
                   const QMap<double,BennettsFreeEnergyAverage> &backwards_ratios)
{
    set(i, BennettsRatios(windows,forwards_ratios,backwards_ratios));
}

/** Set the deltas for the ith iteration */
void Bennetts::set(int i, const BennettsRatios &ratios)
{
    while (i >= rtios.count())
    {
        rtios.append(BennettsRatios());
    }

    i = Index(i).map(rtios.count());

    if (ratios.isEmpty())
        rtios[i] = BennettsRatios();
    else
        rtios[i] = ratios;
}

/** Merge the deltas for iterations start->end */
BennettsRatios Bennetts::merge(int start, int end) const
{
    start = Index(start).map(rtios.count());
    end = Index(end).map(rtios.count());
    
    QList<BennettsRatios> set;
    
    for (int i=start; i<=end; ++i)
    {
        if (not rtios.at(i).isEmpty())
            set.append( rtios.at(i) );
    }
    
    return BennettsRatios::merge(set);
}

/** Merge the deltas at the passed indicies */
BennettsRatios Bennetts::merge(QList<int> indicies) const
{
    QList<BennettsRatios> set;
    
    foreach (int idx, indicies)
    {
        int i = Index(idx).map(rtios.count());
        
        if (not rtios.at(i).isEmpty())
            set.append( rtios.at(i) );
    }
 
    return BennettsRatios::merge(set);
}

/** Return a list of Gradients that represents the rolling average over 'niterations'
    iterations over this TI data set. If this data set contains 100 iterations, and 
    we calculate the rolling average over 50 iterations, then the returned Gradients
    will be the average from 1-50, then 2-51, 3-52.....51-100 */
QList<BennettsRatios> Bennetts::rollingAverage(int niterations) const
{
    QList<BennettsRatios> merged;

    if (rtios.isEmpty())
        return merged;
    
    else if (niterations >= rtios.count())
    {
        BennettsRatios r = this->merge(0,-1);
        
        if (not r.isEmpty())
            merged.append(r);
    }
    else if (niterations <= 1)
    {
        for (int i=0; i<rtios.count(); ++i)
        {
            if (not rtios.at(i).isEmpty())
                merged.append(rtios.at(i));
        }
    }
    else
    {
        QList<BennettsRatios> set;
        
        int i=0;
        
        for (i=0; i<rtios.count(); ++i)
        {
            if (not rtios.at(i).isEmpty())
            {
                set.append(rtios.at(i));
                if (set.count() == niterations)
                    break;
            }
        }
        
        merged.append( BennettsRatios::merge(set) );
        
        for (i=i+1; i<rtios.count(); ++i)
        {
            if (not rtios.at(i).isEmpty())
            {
                set.removeFirst();
                set.append(rtios.at(i));
                merged.append( BennettsRatios::merge(set) );
            }
        }
    }
    
    return merged;
}

/** Remove the data for iteration 'i' */
void Bennetts::removeAt(int i)
{
    i = Index(i).map(rtios.count());
    rtios[i] = BennettsRatios();
}

/** Remove every iteration from 'start' to 'end' (inclusively) */
void Bennetts::removeRange(int start, int end)
{
    start = Index(start).map(rtios.count());
    end = Index(end).map(rtios.count());
    
    if (start > end)
        qSwap(start, end);
    
    for (int i = start; i <= end; ++i)
    {
        rtios[i] = BennettsRatios();
    }
}

/** Remove all values from the histogram */
void Bennetts::clear()
{
    this->operator=( Bennetts() );
}
