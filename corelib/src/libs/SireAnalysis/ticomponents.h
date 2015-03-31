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

#ifndef SIREANALYSIS_TICOMPONENTS_H
#define SIREANALYSIS_TICOMPONENTS_H

#include "SireSystem/freeenergymonitor.h"
#include "SireAnalysis/ti.h"

SIRE_BEGIN_HEADER

namespace SireAnalysis
{
class TIComponents;
class ComponentGradients;
}

QDataStream& operator<<(QDataStream&, const SireAnalysis::TIComponents&);
QDataStream& operator>>(QDataStream&, SireAnalysis::TIComponents&);

QDataStream& operator<<(QDataStream&, const SireAnalysis::ComponentGradients&);
QDataStream& operator>>(QDataStream&, SireAnalysis::ComponentGradients&);

namespace SireAnalysis
{

using SireSystem::FreeEnergyMonitor;

/** This class is used to hold the individual free energy gradients
    for each of the components collected by the SireSystem::FreeEnergyMonitor
    class
    
    @author Christopher Woods
*/
class SIREANALYSIS_EXPORT ComponentGradients
            : public SireBase::ConcreteProperty<ComponentGradients,SireBase::Property>
{

friend QDataStream& ::operator<<(QDataStream&, const ComponentGradients&);
friend QDataStream& ::operator>>(QDataStream&, ComponentGradients&);

public:
    ComponentGradients();
    ComponentGradients(const QMap<double,FreeEnergyMonitor> &gradients,
                       bool conserve_memory=true);
    ComponentGradients(const QList<FreeEnergyMonitor> &gradients,
                       bool conserve_memory=true);
    
    ComponentGradients(const ComponentGradients &other);
    
    ~ComponentGradients();
    
    ComponentGradients& operator=(const ComponentGradients &other);
    
    bool operator==(const ComponentGradients &other) const;
    bool operator!=(const ComponentGradients &other) const;
    
    static const char* typeName();
    
    const char* what() const;
    
    QString toString() const;

    bool isEmpty() const;
    
    bool isCompatible(const ComponentGradients &other) const;
    
    ComponentGradients& operator+=(const ComponentGradients &other);
    ComponentGradients operator+(const ComponentGradients &other) const;
    
    static ComponentGradients merge(const QList<ComponentGradients> &gradients);

    SireUnits::Dimension::Temperature temperature() const;

    QList<double> lambdaValues() const;

    SireMol::PartialMolecule viewAt(int i) const;
    SireMol::PartialMolecule viewAt(int i, double lamval) const;

    Gradients gradientsAt(int i) const;
    Gradients coulombGradientsAt(int i) const;
    Gradients ljGradientsAt(int i) const;

    int nComponents() const;
    int nLambdaValues() const;
    qint64 nSamples() const;

    double deltaLambda() const;
    
    QVector<DataPoint> values(int i) const;
    QVector<DataPoint> coulombValues(int i) const;
    QVector<DataPoint> ljValues(int i) const;

    QMap<double,FreeEnergyMonitor> data() const;

    TIPMF integrate(int i) const;
    TIPMF integrate(int i, int order) const;
    
    TIPMF integrate(int i, double range_min, double range_max) const;
    TIPMF integrate(int i, double range_min, double range_max, int order) const;

    TIPMF integrateCoulomb(int i) const;
    TIPMF integrateCoulomb(int i, int order) const;
    
    TIPMF integrateLJ(int i) const;
    TIPMF integrateLJ(int i, int order) const;
    
    TIPMF integrateCoulomb(int i, double range_min, double range_max) const;
    TIPMF integrateCoulomb(int i, double range_min, double range_max, int order) const;
    
    TIPMF integrateLJ(int i, double range_min, double range_max) const;
    TIPMF integrateLJ(int i, double range_min, double range_max, int order) const;

    void conserveMemory();
    void conserveMemory(const ComponentGradients &other);

private:
    void checkSane() const;

    /** The set of free energy monitors for each lambda value */
    QMap<double,FreeEnergyMonitor> grads;
};

/** This class is used to analyse the free energy components that are
    collected by the SireSystem::FreeEnergyMonitor class
    
    @author Christopher Woods
*/
class SIREANALYSIS_EXPORT TIComponents
            : public SireBase::ConcreteProperty<TIComponents,SireBase::Property>
{

friend QDataStream& ::operator<<(QDataStream&, const TIComponents&);
friend QDataStream& ::operator>>(QDataStream&, TIComponents&);

public:
    TIComponents(bool conserve_memory = true);
    TIComponents(const QMap<double,FreeEnergyMonitor> &gradients,
                 bool conserve_memory = true);
    TIComponents(const ComponentGradients &gradients,
                 bool conserve_memory = true);
    
    TIComponents(const TIComponents &other);
    
    ~TIComponents();
    
    TIComponents& operator=(const TIComponents &other);
    
    bool operator==(const TIComponents &other) const;
    bool operator!=(const TIComponents &other) const;
    
    static const char* typeName();
    
    const char* what() const;
    
    QString toString() const;
    
    void add(const QMap<double,FreeEnergyMonitor> &gradients);
    void add(const ComponentGradients &gradients);
    
    void set(int i, const QMap<double,FreeEnergyMonitor> &gradients);
    void set(int i, const ComponentGradients &gradients);
    
    int nComponents() const;
    int nIterations() const;
    int nLambdaValues() const;
    qint64 nSamples() const;
    
    int count() const;
    int size() const;
    
    bool conservesMemory() const;
    
    QList<double> lambdaValues() const;
    
    ComponentGradients operator[](int i) const;
    ComponentGradients at(int i) const;
    
    QList<ComponentGradients> gradients() const;
    
    ComponentGradients merge(int start, int end) const;
    ComponentGradients merge(QList<int> indicies) const;
    
    QList<ComponentGradients> rollingAverage(int niterations) const;
    
    void removeAt(int i);
    void removeRange(int start, int end);
    
    void clear();

    void conserveMemory();

private:
    /** All of the free energy monitors from each iteration */
    QList<ComponentGradients> grads;
    
    /** Whether or not we should conserve memory */
    bool should_conserve_memory;
};

}

Q_DECLARE_METATYPE( SireAnalysis::TIComponents )
Q_DECLARE_METATYPE( SireAnalysis::ComponentGradients )

SIRE_EXPOSE_CLASS( SireAnalysis::TIComponents )
SIRE_EXPOSE_CLASS( SireAnalysis::ComponentGradients )

SIRE_END_HEADER

#endif
