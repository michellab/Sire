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

#ifndef SIREANALYSIS_FEP_H
#define SIREANALYSIS_FEP_H

#include "SireMaths/freeenergyaverage.h"

#include "SireUnits/dimensions.h"

SIRE_BEGIN_HEADER

namespace SireAnalysis
{
class FEP;
class FEPDeltas;
class DataPoint;
class PMF;
}

QDataStream& operator<<(QDataStream&, const SireAnalysis::FEP&);
QDataStream& operator>>(QDataStream&, SireAnalysis::FEP&);

QDataStream& operator<<(QDataStream&, const SireAnalysis::FEPDeltas&);
QDataStream& operator>>(QDataStream&, SireAnalysis::FEPDeltas&);

QDataStream& operator<<(QDataStream&, const SireAnalysis::DataPoint&);
QDataStream& operator>>(QDataStream&, SireAnalysis::DataPoint&);

QDataStream& operator<<(QDataStream&, const SireAnalysis::PMF&);
QDataStream& operator>>(QDataStream&, SireAnalysis::PMF&);

namespace SireAnalysis
{

using SireMaths::FreeEnergyAverage;

/** This class represents a single datapoint on an x,y graph. The point
    has associated errors (small and large) on both the x and y axes
    
    @author Christopher Woods
*/
class SIREANALYSIS_EXPORT DataPoint
{

friend QDataStream& ::operator<<(QDataStream&, const DataPoint&);
friend QDataStream& ::operator>>(QDataStream&, DataPoint&);

public:
    DataPoint();
    DataPoint(double x, double y);
    DataPoint(double x, double y,
              double xerror, double yerror);
    DataPoint(double x, double y,
              double xminerror, double yminerror,
              double xmaxerror, double ymaxerror);
    
    DataPoint(const DataPoint &other);
    
    ~DataPoint();
    
    DataPoint& operator=(const DataPoint &other);
    
    bool operator==(const DataPoint &other) const;
    bool operator!=(const DataPoint &other) const;
    
    const char* what() const;
    static const char* typeName();
   
    QString toString() const;
    
    double x() const;
    double y() const;
    
    double xError() const;
    double yError() const;
    
    double xMinError() const;
    double yMinError() const;
    
    double xMaxError() const;
    double yMaxError() const;
    
    bool hasError() const;
    bool hasErrorRange() const;
    
    bool hasXError() const;
    bool hasYError() const;
    
    bool hasXErrorRange() const;
    bool hasYErrorRange() const;
    
    bool equalWithinError(const DataPoint &other) const;
    bool equalWithinMinError(const DataPoint &other) const;
    bool equalWithinMaxError(const DataPoint &other) const;
    
private:
    /** The individual values of the data point */
    double _x, _y, _xminerr, _yminerr, _xmaxerr, _ymaxerr;
};

/** This class contains the complete potential of mean force
    
    @author Christopher Woods
*/
class SIREANALYSIS_EXPORT PMF : public SireBase::ConcreteProperty<PMF,SireBase::Property>
{

friend QDataStream& ::operator<<(QDataStream&, const PMF&);
friend QDataStream& ::operator>>(QDataStream&, PMF&);

public:
    PMF();
    PMF(const QVector<DataPoint> &values);
    
    PMF(const PMF &other);
    
    ~PMF();
    
    PMF& operator=(const PMF &other);
    
    bool operator==(const PMF &other) const;
    bool operator!=(const PMF &other) const;
    
    const char* what() const;
    static const char* typeName();

    bool isEmpty() const;

    virtual QString toString() const;

    virtual QVector<DataPoint> values() const;
    
    virtual double rangeMin() const;
    virtual double rangeMax() const;

    virtual double deltaG() const;
    virtual double error() const;

protected:
    void setValues(const QVector<DataPoint> &values);

private:
    /** The smoothed x,y points of the PMF, with errors */
    QVector<DataPoint> vals;
};

/** This class is used to hold the set of FEP deltas from a single 
    iteration 
    
    @author Christopher Woods
*/
class SIREANALYSIS_EXPORT FEPDeltas
      : public SireBase::ConcreteProperty<FEPDeltas,SireBase::Property>
{

friend QDataStream& ::operator<<(QDataStream&, const FEPDeltas&);
friend QDataStream& ::operator>>(QDataStream&, FEPDeltas&);

public:
    FEPDeltas();
    FEPDeltas(const QList<double> &windows, const QMap<double,FreeEnergyAverage> &deltas);
    FEPDeltas(const QList<double> &windows,
              const QMap<double,FreeEnergyAverage> &forwards_deltas,
              const QMap<double,FreeEnergyAverage> &backwards_deltas);
    
    FEPDeltas(const FEPDeltas &other);
    
    ~FEPDeltas();
    
    FEPDeltas& operator=(const FEPDeltas &other);
    
    bool operator==(const FEPDeltas &other) const;
    bool operator!=(const FEPDeltas &other) const;
    
    FEPDeltas& operator+=(const FEPDeltas &other);
    
    FEPDeltas operator+(const FEPDeltas &other) const;
    
    const char* what() const;
    static const char* typeName();
    
    QString toString() const;

    bool isEmpty() const;
    
    static FEPDeltas merge(const QList<FEPDeltas> &deltas);

    SireUnits::Dimension::Temperature temperature() const;

    QList<double> lambdaValues() const;
    QList<double> windows() const;

    int nLambdaValues() const;
    int nWindows() const;
    qint64 nSamples() const;

    QVector<DataPoint> values() const;

    QVector<DataPoint> forwardsValues() const;
    QVector<DataPoint> backwardsValues() const;

    QMap<double,FreeEnergyAverage> forwardsData() const;
    QMap<double,FreeEnergyAverage> backwardsData() const;

    QMap<double,FreeEnergyAverage> forwardsDeltas() const;
    QMap<double,FreeEnergyAverage> backwardsDeltas() const;

    PMF sum() const;
    PMF integrate() const;

    PMF sumForwards() const;
    PMF sumBackwards() const;
    
    PMF sumForwardsTaylor() const;
    PMF sumBackwardsTaylor() const;

private:
    void checkSane() const;

    /** The lambda values of all of the windows. The forwards deltas
        give the free energy between each window and the next one in the 
        list, while the backwards deltas give the free energy between
        each window and the previous one in the list */
    QList<double> lamvals;

    /** The forwards deltas */
    QMap<double,FreeEnergyAverage> fwds_deltas;

    /** The backwards deltas */
    QMap<double,FreeEnergyAverage> bwds_deltas;
};

/** This class is used to analyse the free energies that are
    calculated during a free energy perturbation (FEP) simulation
    
    @author Christopher Woods
*/
class SIREANALYSIS_EXPORT FEP : public SireBase::ConcreteProperty<FEP,SireBase::Property>
{

friend QDataStream& ::operator<<(QDataStream&, const FEP&);
friend QDataStream& ::operator>>(QDataStream&, FEP&);

public:
    FEP();
    
    FEP(const QList<double> &windows, const QMap<double,FreeEnergyAverage> &deltas);
    FEP(const QList<double> &windows,
        const QMap<double,FreeEnergyAverage> &forwards_deltas,
        const QMap<double,FreeEnergyAverage> &backwards_deltas);

    FEP(const FEPDeltas &deltas);
    
    FEP(const FEP &other);
    
    ~FEP();
    
    FEP& operator=(const FEP &other);
    
    bool operator==(const FEP &other) const;
    bool operator!=(const FEP &other) const;
    
    const char* what() const;
    static const char* typeName();
    
    QString toString() const;
    
    void add(const QList<double> &windows,
             const QMap<double,FreeEnergyAverage> &deltas);

    void add(const QList<double> &windows,
             const QMap<double,FreeEnergyAverage> &forwards_deltas,
             const QMap<double,FreeEnergyAverage> &backwards_deltas);
    
    void add(const FEPDeltas &deltas);
    
    int nIterations() const;
    int nLambdaValues() const;
    int nWindows() const;
    qint64 nSamples() const;
    
    int count() const;
    int size() const;
    
    QList<double> lambdaValues() const;
    QList<double> windows() const;
    
    FEPDeltas operator[](int i) const;
    FEPDeltas at(int i) const;
    
    QList<FEPDeltas> deltas() const;
    
    void set(int i, const QList<double> &windows,
             const QMap<double,FreeEnergyAverage> &deltas);
    void set(int i, const QList<double> &windows,
             const QMap<double,FreeEnergyAverage> &forwards_deltas,
             const QMap<double,FreeEnergyAverage> &backwards_deltas);
    
    void set(int i, const FEPDeltas &deltas);
    
    FEPDeltas merge(int start, int end) const;
    FEPDeltas merge(QList<int> indicies) const;
    
    QList<FEPDeltas> rollingAverage(int niterations) const;
    
    void removeAt(int i);
    void removeRange(int start, int end);
    
    void clear();
    
private:
    /** The set of differences between neighbouring lambda
        windows for each iteration */
    QList<FEPDeltas> dltas;
};

}

Q_DECLARE_METATYPE( SireAnalysis::FEP )
Q_DECLARE_METATYPE( SireAnalysis::FEPDeltas)
Q_DECLARE_METATYPE( SireAnalysis::DataPoint )
Q_DECLARE_METATYPE( SireAnalysis::PMF )

SIRE_EXPOSE_CLASS( SireAnalysis::FEP )
SIRE_EXPOSE_CLASS( SireAnalysis::FEPDeltas )
SIRE_EXPOSE_CLASS( SireAnalysis::DataPoint )
SIRE_EXPOSE_CLASS( SireAnalysis::PMF )

SIRE_END_HEADER

#endif
