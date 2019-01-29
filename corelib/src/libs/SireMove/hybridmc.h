/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2010  Christopher Woods
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

#ifndef SIREMOVE_HYBRIDMC_H
#define SIREMOVE_HYBRIDMC_H

#include "moleculardynamics.h"
#include "montecarlo.h"
#include "velocitygenerator.h"

SIRE_BEGIN_HEADER

namespace SireMove
{
class HybridMC;
class HMCVelGen;
class HMCGenerator;
//class MEHMCGenerator;
}

SIREMOVE_EXPORT QDataStream& operator<<(QDataStream&, const SireMove::HybridMC&);
SIREMOVE_EXPORT QDataStream& operator>>(QDataStream&, SireMove::HybridMC&);

SIREMOVE_EXPORT QDataStream& operator<<(QDataStream&, const SireMove::HMCVelGen&);
SIREMOVE_EXPORT QDataStream& operator>>(QDataStream&, SireMove::HMCVelGen&);

SIREMOVE_EXPORT QDataStream& operator<<(QDataStream&, const SireMove::HMCGenerator&);
SIREMOVE_EXPORT QDataStream& operator>>(QDataStream&, SireMove::HMCGenerator&);

namespace SireMove
{

/** This class provides the base class of an extension 
    of the velocity generator that generates velocities 
    in such a way that the bias can be
    calculated and accounted for within the hybrid MC move
    
    @author Christopher Woods
*/
class SIREMOVE_EXPORT HMCVelGen : public VelocityGenerator
{

friend QDataStream& ::operator<<(QDataStream&, const HMCVelGen&);
friend QDataStream& ::operator>>(QDataStream&, HMCVelGen&);

public:
    HMCVelGen();

    HMCVelGen(const HMCVelGen &other);
    
    ~HMCVelGen();

    void setGenerator(const RanGenerator &generator);

    SireUnits::Dimension::Temperature temperature() const;

    virtual void setTemperature(SireUnits::Dimension::Temperature temperature);

    /** Generate the velocites in the passed MD object, 
        returning the biasing factor for the HMC algorithm */
    virtual double generate(const System &system, MolecularDynamics &md)=0;
    
    /** Return the bias for the velocities in the passed MD object */
    virtual double getBias(const MolecularDynamics &md)=0;

protected:
    HMCVelGen& operator=(const HMCVelGen &other);
    
    bool operator==(const HMCVelGen &other) const;
    bool operator!=(const HMCVelGen &other) const;

    const RanGenerator& rand() const;

private:
    /** The random number generator used to generate the random
        numbers needed by the velocities */
    RanGenerator rangen;

    /** The temperature used to control the generation of 
        the velocities - this is the temperature of the MC move */
    SireUnits::Dimension::Temperature temp;
};

/** This is the velocity generator used for the standard
    Hybrid Monte Carlo move. This generates velocities such
    that the acceptance test is just based on the change
    in total energy, which, for a good integrator, is 
    near zero, so nearly all moves should be accepted */
class SIREMOVE_EXPORT HMCGenerator 
        : public SireBase::ConcreteProperty<HMCGenerator,HMCVelGen>
{

friend QDataStream& ::operator<<(QDataStream&, const HMCGenerator&);
friend QDataStream& ::operator>>(QDataStream&, HMCGenerator&);

public:
    HMCGenerator();
    HMCGenerator(const HMCGenerator &other);
    
    ~HMCGenerator();
    
    static const char* typeName();
    
    HMCGenerator& operator=(const HMCGenerator &other);
    
    bool operator==(const HMCGenerator &other) const;
    bool operator!=(const HMCGenerator &other) const;
    
    double generate(const System &system, MolecularDynamics &md);
    double getBias(const MolecularDynamics &md);
    
    AtomVelocities generate(const MoleculeView &molview, 
                            const PropertyMap &map = PropertyMap()) const;
};

/** This class implements a hybrid Monte Carlo move. This is a 
    Monte Carlo move that uses a symplectic, time-reversible
    NVE integrator to run some molecular dynamics. The block
    of MD is accepted according to a Metropolis MC test
    on the change in total energy (kinetic+potential).
    
    @author Christopher Woods
*/
class SIREMOVE_EXPORT HybridMC : public SireBase::ConcreteProperty<HybridMC,MonteCarlo>
{

friend QDataStream& ::operator<<(QDataStream&, const HybridMC&);
friend QDataStream& ::operator>>(QDataStream&, HybridMC&);

public:
    HybridMC(const PropertyMap &map = PropertyMap());

    HybridMC(const MoleculeGroup &molgroup,
             const PropertyMap &map = PropertyMap());

    HybridMC(const MoleculeGroup &molgroup,
             const Integrator &integrator,
             const PropertyMap &map = PropertyMap());

    HybridMC(const MoleculeGroup &molgroup,
             SireUnits::Dimension::Time timestep,
             const PropertyMap &map = PropertyMap());

    HybridMC(const MoleculeGroup &molgroup,
             const Integrator &integrator,
             SireUnits::Dimension::Time timestep,
             const PropertyMap &map = PropertyMap());

    HybridMC(const MoleculeGroup &molgroup, int nsteps,
             const PropertyMap &map = PropertyMap());

    HybridMC(const MoleculeGroup &molgroup,
             const Integrator &integrator, int nsteps,
             const PropertyMap &map = PropertyMap());

    HybridMC(const MoleculeGroup &molgroup,
             SireUnits::Dimension::Time timestep, int nsteps,
             const PropertyMap &map = PropertyMap());

    HybridMC(const MoleculeGroup &molgroup,
             const Integrator &integrator,
             SireUnits::Dimension::Time timestep, int nsteps,
             const PropertyMap &map = PropertyMap());

    HybridMC(const HybridMC &other);

    ~HybridMC();
    
    static const char* typeName();
    
    HybridMC& operator=(const HybridMC &other);
    
    bool operator==(const HybridMC &other) const;
    bool operator!=(const HybridMC &other) const;
    
    QString toString() const;

    const MoleculeGroup& moleculeGroup() const;

    void setCoordinatesProperty(const PropertyName &coords_property);
    void setSpaceProperty(const PropertyName &space_property);

    void setGenerator(const RanGenerator &rangenerator);
    
    void setVelocityGenerator(const HMCVelGen &generator);
    
    const HMCVelGen& velocityGenerator() const;
    
    int nDynamicsSteps();

    void setNDynamicsSteps(int nsteps);

    SireUnits::Dimension::Time timeStep() const;
    
    void setTimeStep(SireUnits::Dimension::Time timestep);

    void move(System &system, int nmoves, bool record_stats=true);

private:
    void assertValidIntegrator(const Integrator &integrator) const;

    /** The MD object used to run the molecular dynamics moves */
    MolecularDynamics md;
    
    /** The generator used to generate the random velocities for each move */
    VelGenPtr velgen;
    
    /** The number of MD moves to perform per HybridMC move */
    quint32 nsteps;
};

}

Q_DECLARE_METATYPE( SireMove::HybridMC )
Q_DECLARE_METATYPE( SireMove::HMCGenerator )

SIRE_EXPOSE_CLASS( SireMove::HybridMC )
SIRE_EXPOSE_CLASS( SireMove::HMCVelGen )
SIRE_EXPOSE_CLASS( SireMove::HMCGenerator )

SIRE_END_HEADER

#endif
