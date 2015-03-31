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

#ifndef SIRESYSTEM_FREEENERGYMONITOR_H
#define SIRESYSTEM_FREEENERGYMONITOR_H

#include "SireSystem/systemmonitor.h"
#include "SireSystem/idassigner.h"

#include "SireMaths/freeenergyaverage.h"

#include "SireCAS/symbol.h"

#include "SireMol/moleculegroup.h"
#include "SireMol/partialmolecule.h"

#include "idassigner.h"

SIRE_BEGIN_HEADER

namespace SireSystem
{
class AssignerGroup;
class FreeEnergyMonitor;
}

QDataStream& operator<<(QDataStream&, const SireSystem::AssignerGroup&);
QDataStream& operator>>(QDataStream&, SireSystem::AssignerGroup&);

QDataStream& operator<<(QDataStream&, const SireSystem::FreeEnergyMonitor&);
QDataStream& operator>>(QDataStream&, SireSystem::FreeEnergyMonitor&);

namespace SireSystem
{

using SireMol::MoleculeGroup;
using SireMaths::Accumulator;

/** This is a simple class that holds either a MoleculeGroup or an IDAssigner

    @author Christopher Woods
*/
class SIRESYSTEM_EXPORT AssignerGroup
{

friend QDataStream& ::operator<<(QDataStream&, const AssignerGroup&);
friend QDataStream& ::operator>>(QDataStream&, AssignerGroup&);

public:
    AssignerGroup();
    AssignerGroup(const MoleculeGroup &molgroup);
    AssignerGroup(const IDAssigner &assigner);
    
    AssignerGroup(const AssignerGroup &other);
    
    ~AssignerGroup();
    
    AssignerGroup& operator=(const AssignerGroup &other);
    
    bool operator==(const AssignerGroup &other) const;
    bool operator!=(const AssignerGroup &other) const;
    
    static const char* typeName();
    
    const char* what() const;
    
    bool isEmpty() const;
    bool isMoleculeGroup() const;
    bool isAssigner() const;
    
    bool isCompatible(const AssignerGroup &other) const;
    
    const MoleculeGroup &group() const;
    const IDAssigner &assigner() const;

    QVector<PartialMolecule> views() const;
    
    void update(const System &system);

private:
    /** MoleculeGroup pointer, if this contains a molecule group */
    SireMol::MolGroupPtr molgroup;
    
    /** Pointer to the IDAssigner, if this contains one */
    SireBase::PropertyPtr assgnr;
};

/** This monitor is used to monitor the free energy difference of two
    molecule groups against every molecule view in the third. 
    This uses dual topology to calculate
    the free energy difference between the interaction of each view in the
    reference molecule group, and the group A and group B groups. This is intended
    to complement free energy calculations by letting you decompose
    the free energy difference into per-residue components. The coulomb
    and LJ components are also separately calculated and accumulated.
    
    @author Christopher Woods
*/
class SIRESYSTEM_EXPORT FreeEnergyMonitor
        : public SireBase::ConcreteProperty<FreeEnergyMonitor,SystemMonitor>
{

friend QDataStream& ::operator<<(QDataStream&, const FreeEnergyMonitor&);
friend QDataStream& ::operator>>(QDataStream&, FreeEnergyMonitor&);

public:
    FreeEnergyMonitor();
    FreeEnergyMonitor(const AssignerGroup &reference,
                      const AssignerGroup &groupA, const AssignerGroup &groupB);

    FreeEnergyMonitor(const FreeEnergyMonitor &other);
    
    ~FreeEnergyMonitor();
    
    FreeEnergyMonitor& operator=(const FreeEnergyMonitor &other);
    
    bool operator==(const FreeEnergyMonitor &other) const;
    bool operator!=(const FreeEnergyMonitor &other) const;
    
    FreeEnergyMonitor& operator+=(const FreeEnergyMonitor &other);
    
    FreeEnergyMonitor operator+(const FreeEnergyMonitor &other) const;

    static FreeEnergyMonitor merge(const QList<FreeEnergyMonitor> &monitors);
    
    static const char* typeName();
    
    bool isEmpty() const;
    bool isCompatible(const FreeEnergyMonitor &other) const;
    bool isCompatibleExceptLambda(const FreeEnergyMonitor &other) const;
    
    void clearStatistics();
    
    int nSamples() const;
    
    void monitor(System &system);

    void setLambdaComponent(const SireCAS::Symbol &component);
    
    void setShiftDelta(double delta);
    void setCoulombPower(int power);

    void setDeltaLambda(double delta_lam);

    void setTemperature(const SireUnits::Dimension::Temperature &temperature);
    void setBinWidth(const SireUnits::Dimension::MolarEnergy &binwidth);

    SireUnits::Dimension::Temperature temperature() const;
    SireUnits::Dimension::MolarEnergy binWidth() const;

    double shiftDelta() const;
    int coulombPower() const;

    double deltaLambda() const;

    bool usesSoftCore() const;
    
    SireCAS::Symbol lambdaComponent() const;
    double lambdaValue() const;

    QVector<SireMol::PartialMolecule> referenceViews() const;

    const AssignerGroup& referenceGroup() const;
    const AssignerGroup& groupA() const;
    const AssignerGroup& groupB() const;

    QVector<SireMaths::FreeEnergyAverage> freeEnergies() const;
    QVector<SireMaths::FreeEnergyAverage> coulombFreeEnergies() const;
    QVector<SireMaths::FreeEnergyAverage> ljFreeEnergies() const;

    void conserveMemory(const FreeEnergyMonitor &other);

private:
    /** The reference group that contains the groups against which
        the free energy will be calculated */
    AssignerGroup refgroup;
    
    /** The two groups over which the free energy will be calculated.
        The free energy difference of group A interacting with the 
        reference group, and group B interacting with the reference
        group will be evaluated */
    AssignerGroup group_a, group_b;
    
    /** The accumulated total free energy differences */
    QVector<SireMaths::FreeEnergyAverage> total_nrgs;
    
    /** The accumulated coulomb free energies */
    QVector<SireMaths::FreeEnergyAverage> coul_nrgs;
    
    /** The accumulated LJ free energies */
    QVector<SireMaths::FreeEnergyAverage> lj_nrgs;

    /** The template for all of the FreeEnergyAverages */
    SireMaths::FreeEnergyAverage nrg_template;

    /** The Symbol used to get the current value of lambda */
    SireCAS::Symbol lambda_symbol;

    /** The shift-delta parameter */
    double shift_delta;

    /** The value of lambda */
    double lamval;

    /** The value of delta lambda */
    double delta_lambda;

    /** The coulomb power parameter */
    quint32 coulomb_power;
};

}

Q_DECLARE_METATYPE( SireSystem::FreeEnergyMonitor )
Q_DECLARE_METATYPE( SireSystem::AssignerGroup )

SIRE_EXPOSE_CLASS( SireSystem::FreeEnergyMonitor )
SIRE_EXPOSE_CLASS( SireSystem::AssignerGroup )

SIRE_END_HEADER

#endif
