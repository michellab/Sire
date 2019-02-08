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

#ifndef SIREANALYSIS_TI_H
#define SIREANALYSIS_TI_H

#include "SireAnalysis/fep.h"

SIRE_BEGIN_HEADER

namespace SireAnalysis
{
class TI;
class Gradients;
class TIPMF;
}

SIREANALYSIS_EXPORT QDataStream& operator<<(QDataStream&, const SireAnalysis::TI&);
SIREANALYSIS_EXPORT QDataStream& operator>>(QDataStream&, SireAnalysis::TI&);

SIREANALYSIS_EXPORT QDataStream& operator<<(QDataStream&, const SireAnalysis::Gradients&);
SIREANALYSIS_EXPORT QDataStream& operator>>(QDataStream&, SireAnalysis::Gradients&);

SIREANALYSIS_EXPORT QDataStream& operator<<(QDataStream&, const SireAnalysis::TIPMF&);
SIREANALYSIS_EXPORT QDataStream& operator>>(QDataStream&, SireAnalysis::TIPMF&);

namespace SireAnalysis
{

using SireMaths::AverageAndStddev;
using SireMaths::FreeEnergyAverage;
using SireUnits::Dimension::MolarEnergy;

/** This class contains the complete potential of mean force
    that has been created by integrating the TI gradients
  
    @author Christopher Woods
*/
class SIREANALYSIS_EXPORT TIPMF : public SireBase::ConcreteProperty<TIPMF,PMF>
{

friend SIREANALYSIS_EXPORT QDataStream& ::operator<<(QDataStream&, const TIPMF&);
friend SIREANALYSIS_EXPORT QDataStream& ::operator>>(QDataStream&, TIPMF&);

public:
    TIPMF();
    TIPMF(int order);
    
    TIPMF(double min_x, double max_x);
    TIPMF(double min_x, double max_x, int order);
    
    TIPMF(const TIPMF &other);
    
    ~TIPMF();
    
    TIPMF& operator=(const TIPMF &other);
    
    bool operator==(const TIPMF &other) const;
    bool operator!=(const TIPMF &other) const;
    
    const char* what() const;
    static const char* typeName();

    QString toString() const;

    QVector<DataPoint> gradients() const;
    QVector<DataPoint> smoothedGradients() const;

    void setOrder(qint32 order);
    void setRange(double min_x, double max_x);

    void setGradients(const QVector<DataPoint> &gradients);
    
    int order() const;
    
    double rangeMin() const;
    double rangeMax() const;
    
    double integral() const;
    double quadrature() const;

    TIPMF dropEndPoints() const;

private:
    void recalculate();

    /** The raw x,y points of the PMF gradients, complete with errors */
    QVector<DataPoint> grads;
    
    /** The smoothed x,y gradients, with errors */
    QVector<DataPoint> smoothed_grads;
    
    /** The minimum value of the range of the PMF */
    double range_min;
    
    /** The maximum value of the range of the PMF */
    double range_max;
    
    /** The free energy from quadrature */
    double quad_value;
    
    /** The number of polynomials to use to fit the gradients */
    qint32 npoly;
};

/** This class contains the free energy gradients from a TI simulation

    @author Christopher Woods
*/
class SIREANALYSIS_EXPORT Gradients
        : public SireBase::ConcreteProperty<Gradients,SireBase::Property>
{

friend SIREANALYSIS_EXPORT QDataStream& ::operator<<(QDataStream&, const Gradients&);
friend SIREANALYSIS_EXPORT QDataStream& ::operator>>(QDataStream&, Gradients&);

public:
    Gradients();
    Gradients(const QMap<double,AverageAndStddev> &gradients);
    Gradients(const QMap<double,FreeEnergyAverage> &gradients);
    Gradients(const QMap<double,FreeEnergyAverage> &gradients,
              double delta_lambda);
    Gradients(const QMap<double,FreeEnergyAverage> &forwards,
              const QMap<double,FreeEnergyAverage> &backwards,
              double delta_lambda);
    
    Gradients(const Gradients &other);
    
    ~Gradients();
    
    Gradients& operator=(const Gradients &other);
    
    bool operator==(const Gradients &other) const;
    bool operator!=(const Gradients &other) const;

    MolarEnergy operator[](double lam) const;

    static const char* typeName();
    
    const char* what() const;
    
    bool isEmpty() const;
    
    QString toString() const;
    
    Gradients& operator+=(const Gradients &other);
    Gradients operator+(const Gradients &other) const;
    
    static Gradients merge(const QList<Gradients> &gradients);

    SireUnits::Dimension::Temperature temperature() const;

    QList<double> lambdaValues() const;
    QList<double> keys() const;

    int nLambdaValues() const;
    qint64 nSamples() const;

    double deltaLambda() const;
    
    QVector<DataPoint> values() const;
    
    QVector<DataPoint> forwardsValues() const;
    QVector<DataPoint> backwardsValues() const;

    QMap<double,AverageAndStddev> analyticData() const;

    QMap<double,FreeEnergyAverage> forwardsData() const;
    QMap<double,FreeEnergyAverage> backwardsData() const;

    MolarEnergy forwards(double lam) const;
    MolarEnergy backwards(double lam) const;
    MolarEnergy gradient(double lam) const;

    TIPMF integrate() const;
    TIPMF integrate(int order) const;
    
    TIPMF integrate(double range_min, double range_max) const;
    TIPMF integrate(double range_min, double range_max, int order) const;

private:
    void checkSane() const;

    /** The analytic gradient data (from analytic TI) */
    QMap<double,AverageAndStddev> analytic;

    /** The forwards values (from finite difference TI) */
    QMap<double,FreeEnergyAverage> fwds;

    /** The backwards values */
    QMap<double,FreeEnergyAverage> bwds;

    /** The value of delta lambda */
    double delta_lam;
};

/** This class is used to analyse the free energies that are
    calculated during a thermodynamic integration simulation
    
    @author Christopher Woods
*/
class SIREANALYSIS_EXPORT TI : public SireBase::ConcreteProperty<TI,SireBase::Property>
{

friend SIREANALYSIS_EXPORT QDataStream& ::operator<<(QDataStream&, const TI&);
friend SIREANALYSIS_EXPORT QDataStream& ::operator>>(QDataStream&, TI&);

public:
    TI();
    TI(const Gradients &gradients);
    TI(const QList<Gradients> &gradients);
    
    TI(const TI &other);
    
    ~TI();
    
    TI& operator=(const TI &other);
    
    bool operator==(const TI &other) const;
    bool operator!=(const TI &other) const;
    
    const char* what() const;
    static const char* typeName();
    
    QString toString() const;
    
    void add(const QMap<double,AverageAndStddev> &gradients);
    void add(const QMap<double,FreeEnergyAverage> &gradients);
    void add(const QMap<double,FreeEnergyAverage> &gradients,
             double delta_lambda);
    
    void add(const QMap<double,FreeEnergyAverage> &forwards,
             const QMap<double,FreeEnergyAverage> &backwards,
             double delta_lambda);
    
    void add(const Gradients &gradients);
    
    int nIterations() const;
    int nLambdaValues() const;
    qint64 nSamples() const;
    
    int count() const;
    int size() const;
    
    QList<double> lambdaValues() const;
    
    Gradients operator[](int i) const;
    Gradients at(int i) const;
    
    QList<Gradients> gradients() const;
    
    void set(int i, const QMap<double,AverageAndStddev> &gradients);
    void set(int i, const QMap<double,FreeEnergyAverage> &gradients);
    void set(int i, const QMap<double,FreeEnergyAverage> &gradients,
             double delta_lambda);
    
    void set(int i, const QMap<double,FreeEnergyAverage> &forwards,
                    const QMap<double,FreeEnergyAverage> &backwards,
                    double delta_lambda);

    void set(int i, const Gradients &gradients);
    
    Gradients merge(int start, int end) const;
    Gradients merge(QList<int> indicies) const;
    
    QList<Gradients> rollingAverage(int niterations) const;
    
    void removeAt(int i);
    void removeRange(int start, int end);
    
    void clear();
    
private:
    /** The set of gradients for each iteration */
    QList<Gradients> grads;
};

}

Q_DECLARE_METATYPE( SireAnalysis::Gradients )
Q_DECLARE_METATYPE( SireAnalysis::TI )
Q_DECLARE_METATYPE( SireAnalysis::TIPMF )

SIRE_EXPOSE_CLASS( SireAnalysis::Gradients )
SIRE_EXPOSE_CLASS( SireAnalysis::TI )
SIRE_EXPOSE_CLASS( SireAnalysis::TIPMF )

SIRE_END_HEADER

#endif

