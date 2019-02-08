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

#ifndef SIREANALYSIS_BENNETS_H
#define SIREANALYSIS_BENNETS_H

#include "fep.h"

SIRE_BEGIN_HEADER

namespace SireAnalysis
{
class Bennetts;
class BennettsRatios;
}

SIREANALYSIS_EXPORT QDataStream& operator<<(QDataStream&, const SireAnalysis::Bennetts&);
SIREANALYSIS_EXPORT QDataStream& operator>>(QDataStream&, SireAnalysis::Bennetts&);

SIREANALYSIS_EXPORT QDataStream& operator<<(QDataStream&, const SireAnalysis::BennettsRatios&);
SIREANALYSIS_EXPORT QDataStream& operator>>(QDataStream&, SireAnalysis::BennettsRatios&);

namespace SireAnalysis
{

using SireMaths::BennettsFreeEnergyAverage;

/** This class is used to hold a set of Bennets acceptance ratios
    for a single iteration
    
    @author Christopher Woods
*/
class SIREANALYSIS_EXPORT BennettsRatios
        : public SireBase::ConcreteProperty<BennettsRatios,SireBase::Property>
{

friend SIREANALYSIS_EXPORT QDataStream& ::operator<<(QDataStream&, const BennettsRatios&);
friend SIREANALYSIS_EXPORT QDataStream& ::operator>>(QDataStream&, BennettsRatios&);

public:
    BennettsRatios();
    BennettsRatios(const QList<double> &windows,
                   const QMap<double,BennettsFreeEnergyAverage> &forwards_ratios,
                   const QMap<double,BennettsFreeEnergyAverage> &backwards_ratios);
    
    BennettsRatios(const BennettsRatios &other);
    
    ~BennettsRatios();
    
    BennettsRatios& operator=(const BennettsRatios &other);
    
    bool operator==(const BennettsRatios &other) const;
    bool operator!=(const BennettsRatios &other) const;
    
    BennettsRatios& operator+=(const BennettsRatios &other);
    
    BennettsRatios operator+(const BennettsRatios &other) const;
    
    const char* what() const;
    static const char* typeName();
    
    QString toString() const;

    bool isEmpty() const;
    
    static BennettsRatios merge(const QList<BennettsRatios> &deltas);

    SireUnits::Dimension::Temperature temperature() const;

    QList<double> lambdaValues() const;
    QList<double> windows() const;

    int nLambdaValues() const;
    int nWindows() const;
    qint64 nSamples() const;

    QVector<DataPoint> values() const;

    QVector<DataPoint> numerators() const;
    QVector<DataPoint> denominators() const;

    QVector<DataPoint> constants() const;

    QMap<double,BennettsFreeEnergyAverage> forwardsData() const;
    QMap<double,BennettsFreeEnergyAverage> backwardsData() const;

    QMap<double,BennettsFreeEnergyAverage> forwardsRatios() const;
    QMap<double,BennettsFreeEnergyAverage> backwardsRatios() const;

    PMF sum() const;
    PMF integrate() const;

private:
    void checkSane() const;

    /** The lambda values of all of the windows. The forwards ratios
        give the free energy between each window and the next one in the 
        list, while the backwards ratios give the free energy between
        each window and the previous one in the list */
    QList<double> lamvals;

    /** The forwards deltas */
    QMap<double,BennettsFreeEnergyAverage> fwds_ratios;

    /** The backwards deltas */
    QMap<double,BennettsFreeEnergyAverage> bwds_ratios;
};

/** This class is used to analyse the free energies that are
    calculated during a Bennetts Acceptance Ratio simulation
    
    @author Christopher Woods
*/
class SIREANALYSIS_EXPORT Bennetts : public SireBase::ConcreteProperty<Bennetts,SireBase::Property>
{

friend SIREANALYSIS_EXPORT QDataStream& ::operator<<(QDataStream&, const Bennetts&);
friend SIREANALYSIS_EXPORT QDataStream& ::operator>>(QDataStream&, Bennetts&);

public:
    Bennetts();
    
    Bennetts(const QList<double> &windows,
             const QMap<double,BennettsFreeEnergyAverage> &forwards_ratios,
             const QMap<double,BennettsFreeEnergyAverage> &backwards_ratios);

    Bennetts(const BennettsRatios &ratios);
    
    Bennetts(const Bennetts &other);
    
    ~Bennetts();
    
    Bennetts& operator=(const Bennetts &other);
    
    bool operator==(const Bennetts &other) const;
    bool operator!=(const Bennetts &other) const;
    
    const char* what() const;
    static const char* typeName();
    
    QString toString() const;
    
    void add(const QList<double> &windows,
             const QMap<double,BennettsFreeEnergyAverage> &forwards_ratios,
             const QMap<double,BennettsFreeEnergyAverage> &backwards_ratios);
    
    void add(const BennettsRatios &ratios);
    
    int nIterations() const;
    int nLambdaValues() const;
    int nWindows() const;
    qint64 nSamples() const;
    
    int count() const;
    int size() const;
    
    QList<double> lambdaValues() const;
    QList<double> windows() const;
    
    BennettsRatios operator[](int i) const;
    BennettsRatios at(int i) const;
    
    QList<BennettsRatios> ratios() const;
    
    void set(int i, const QList<double> &windows,
             const QMap<double,BennettsFreeEnergyAverage> &forwards_ratios,
             const QMap<double,BennettsFreeEnergyAverage> &backwards_ratios);
    
    void set(int i, const BennettsRatios &ratios);
    
    BennettsRatios merge(int start, int end) const;
    BennettsRatios merge(QList<int> indicies) const;
    
    QList<BennettsRatios> rollingAverage(int niterations) const;
    
    void removeAt(int i);
    void removeRange(int start, int end);
    
    void clear();
    
private:
    /** The set of Bennetts ratios for neighbouring lambda
        windows for each iteration */
    QList<BennettsRatios> rtios;
};

}

Q_DECLARE_METATYPE( SireAnalysis::BennettsRatios )
Q_DECLARE_METATYPE( SireAnalysis::Bennetts )

SIRE_EXPOSE_CLASS( SireAnalysis::BennettsRatios )
SIRE_EXPOSE_CLASS( SireAnalysis::Bennetts )

SIRE_END_HEADER

#endif
