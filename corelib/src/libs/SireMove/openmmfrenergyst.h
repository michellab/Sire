/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2009  Christopher Woods
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

#ifndef SIREMOVE_OPENMMFRENERGYST_H
#define SIREMOVE_OPENMMFRENERGYST_H

#include "integrator.h"

#ifdef SIRE_USE_OPENMM
  #include <OpenMM.h>   // CONDITIONAL_INCLUDE
#endif

#include <cstdio>
#include "SireUnits/temperature.h"
#include "SireSystem/system.h"
SIRE_BEGIN_HEADER

#ifdef SIRE_USE_OPENMM

namespace SireMove
{
class OpenMMFrEnergyST;
}

QDataStream& operator<<(QDataStream&, const SireMove::OpenMMFrEnergyST&);
QDataStream& operator>>(QDataStream&, SireMove::OpenMMFrEnergyST&);

namespace SireMove
{

/** This class implements single topology a free energy method using OpenMM. 
 
    @author Julien Michel,Gaetano Calabro and Antonia Mey
*/
class SIREMOVE_EXPORT OpenMMFrEnergyST
        : public SireBase::ConcreteProperty<OpenMMFrEnergyST,Integrator>
{

friend QDataStream& ::operator<<(QDataStream&, const OpenMMFrEnergyST&);
friend QDataStream& ::operator>>(QDataStream&, OpenMMFrEnergyST&);

public:
    OpenMMFrEnergyST(bool frequent_save_velocities = false);

    OpenMMFrEnergyST(const MoleculeGroup &molecule_group,
		     const MoleculeGroup &solutes,
		     const MoleculeGroup &solute_hard, 
		     const MoleculeGroup &solute_todummy, 
		     const MoleculeGroup & solute_fromdummy,
		     bool frequent_save_velocities = false);

    OpenMMFrEnergyST(const OpenMMFrEnergyST &other);

    ~OpenMMFrEnergyST();

    OpenMMFrEnergyST& operator=(const OpenMMFrEnergyST &other);

    bool operator==(const OpenMMFrEnergyST &other) const;
    bool operator!=(const OpenMMFrEnergyST &other) const;

    static const char* typeName();

    QString toString() const;

    Ensemble ensemble() const;

    bool isTimeReversible() const;

    void initialise();
    
    
    SireUnits::Dimension::MolarEnergy getPotentialEnergy(const System &system);
    
    System minimizeEnergy(System &system, double tolerance, int max_iteration); 
    
    System annealLambda(System &system, SireUnits::Dimension::Time timestep, 
                        int annealingSteps);

    void integrate(IntegratorWorkspace &workspace,
                   const Symbol &nrg_component,
                   SireUnits::Dimension::Time timestep,
                   int nmoves, bool record_stats) ;

    IntegratorWorkspacePtr createWorkspace(const PropertyMap &map = PropertyMap()) const;
    IntegratorWorkspacePtr createWorkspace(const MoleculeGroup &molgroup,const PropertyMap &map = PropertyMap()) const;

    QString getCutoffType(void);
    void setCutoffType(QString);

    SireUnits::Dimension::Length getCutoff_distance(void);
    void setCutoff_distance(SireUnits::Dimension::Length);

    double getField_dielectric(void);
    void setField_dielectric(double);

    bool getAndersen(void);
    void setAndersen(bool);

    double getAndersen_frequency(void);
    void setAndersen_frequency(double);

    bool getMCBarostat(void);
    void setMCBarostat(bool);

    void setMCBarostat_frequency(int);
    int getMCBarostat_frequency(void);

    QString getConstraintType(void);
    void setConstraintType(QString);

    SireUnits::Dimension::Pressure getPressure(void);
    void setPressure(SireUnits::Dimension::Pressure);

    SireUnits::Dimension::Temperature getTemperature(void);
    void setTemperature(SireUnits::Dimension::Temperature);

    QString getPlatform(void);
    void setPlatform(QString);

    bool getRestraint(void);
    void setRestraint(bool);

    int getCMMremoval_frequency(void);
    void setCMMremoval_frequency(int);

    int getBufferFrequency();
    void setBufferFrequency(int);

    int getEnergyFrequency();
    void setEnergyFrequency(int);

    void setDeviceIndex(QString);
    QString getDeviceIndex(void);

    void setPrecision(QString);
    QString getPrecision(void);

    double getAlchemical_value(void);
    void setAlchemical_value(double);

    float getCoulomb_power(void);
    void setCoulomb_power(float);

    double getShift_delta(void);
    void setShift_delta(double);

    double getDeltaAlchemical(void);
    void setDeltatAlchemical(double);

    QVector<double> getGradients(void);
    QVector<double> getEnergies(void);

    QString getIntegrator(void);
    void setIntegrator(QString);

    SireUnits::Dimension::Time getFriction(void);

    void setFriction(SireUnits::Dimension::Time);

    double getIntegration_tollerance(void);
    void setIntegration_tollerance(double tollerance);

    SireUnits::Dimension::Time getTimetoSkip(void);
    void setTimetoSkip(SireUnits::Dimension::Time);


    int getEquilib_iterations(void);
    void setEquilib_iterations(int);

    SireUnits::Dimension::Time getEquilib_time_step(void);
    void setEquilib_time_step(SireUnits::Dimension::Time);

    void setReinitializeContext(bool);

    int getRandomSeed(void);
    void setRandomSeed(int);

private:
    void createContext(IntegratorWorkspace &workspace,
                       SireUnits::Dimension::Time timestep);
    void destroyContext();
    
    /** Whether or not to save the velocities after every step, or to save them at the end of all of the steps */
    bool frequent_save_velocities;
    /** The Molecule Group on which the integrator operates */
    MolGroupPtr molgroup;
    /** The solute group on which the integrator operates */
    MolGroupPtr solute;
    /** The Solute hard Group on which the integrator operates */
    MolGroupPtr solutehard;
    /** The To Dummy Solute Group on which the integrator operates */
    MolGroupPtr solutetodummy;
    /** The From Dummy Solute Group on which the integrator operates */
    MolGroupPtr solutefromdummy;

    /**Try instead to...keep a copy of OpenMM::System */
    OpenMM::System* openmm_system;

    OpenMM::Context* openmm_context;

    /** Whether the openmm system and the context have been initialised*/
    bool isSystemInitialised;
    bool isContextInitialised;

    QString CutoffType;
    SireUnits::Dimension::Length cutoff_distance;
    double field_dielectric;

    bool Andersen_flag;
    double Andersen_frequency;

    bool MCBarostat_flag;
    int MCBarostat_frequency;

    QString ConstraintType;

    SireUnits::Dimension::Pressure Pressure;
    SireUnits::Dimension::Temperature Temperature;

    QString platform_type;

    bool Restraint_flag;

    int CMMremoval_frequency;
   
    int buffer_frequency;

    int energy_frequency;

    QString device_index;
    
    QString precision;

    double Alchemical_value;

    float coulomb_power;

    double shift_delta;

    double delta_alchemical;

    QVector<double> gradients;

    QVector<double> energies; 

    QVector<bool> perturbed_energies;

    QString Integrator_type;

    SireUnits::Dimension::Time friction;

    double integration_tol;

    SireUnits::Dimension::Time timeskip;

    int equilib_iterations;
    SireUnits::Dimension::Time equilib_time_step;

    bool reinetialize_context;

    double GF_acc;
    double GB_acc;
    
    bool Debug;

    int random_seed;

};


}

Q_DECLARE_METATYPE( SireMove::OpenMMFrEnergyST )

SIRE_EXPOSE_CLASS( SireMove::OpenMMFrEnergyST )

SIRE_END_HEADER

#else // SIRE_USE_OPENMM

namespace SireMove
{

    class OpenMMFrEnergyST{
        public:
            OpenMMFrEnergyST(){}
            ~OpenMMFrEnergyST(){}

            static const char* typeName(){ return "SireMM::OpenMMFrEnergyST"; }

    };

}

Q_DECLARE_METATYPE( SireMove::OpenMMFrEnergyST )

#endif // SIRE_USE_OPENMM

#endif
