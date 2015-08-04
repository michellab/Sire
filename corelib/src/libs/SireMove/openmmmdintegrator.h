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

#ifndef SIREMOVE_OPENMMMDINTEGRATOR_H
#define SIREMOVE_OPENMMMDINTEGRATOR_H

#include "integrator.h"

#ifdef SIRE_USE_OPENMM
#include <OpenMM.h>   // CONDITIONAL_INCLUDE
#endif

#include <cstdio>
#include "SireUnits/temperature.h"
#include "SireSystem/system.h"
SIRE_BEGIN_HEADER

#ifdef SIRE_USE_OPENMM

        namespace SireMove {
    class OpenMMMDIntegrator;
}

QDataStream& operator<<(QDataStream&, const SireMove::OpenMMMDIntegrator&);
QDataStream& operator>>(QDataStream&, SireMove::OpenMMMDIntegrator&);

namespace SireMove {

    /** This class implements a "pure" MD integrator using OpenMM. 
     * No free energy methods are supported.

        @author Julien Michel and Gaetano Calabro
     */
    class SIREMOVE_EXPORT OpenMMMDIntegrator
    : public SireBase::ConcreteProperty<OpenMMMDIntegrator, Integrator> {
        friend QDataStream& ::operator<<(QDataStream&, const OpenMMMDIntegrator&);
        friend QDataStream& ::operator>>(QDataStream&, OpenMMMDIntegrator&);

    public:
        OpenMMMDIntegrator(bool frequent_save_velocities = false);

        OpenMMMDIntegrator(const MoleculeGroup &molecule_group, bool frequent_save_velocities = false);

        OpenMMMDIntegrator(const OpenMMMDIntegrator &other);

        ~OpenMMMDIntegrator();

        OpenMMMDIntegrator& operator=(const OpenMMMDIntegrator &other);

        bool operator==(const OpenMMMDIntegrator &other) const;
        bool operator!=(const OpenMMMDIntegrator &other) const;

        static const char* typeName();

        QString toString() const;

        Ensemble ensemble() const;

        bool isTimeReversible() const;

        void initialise();

        SireUnits::Dimension::MolarEnergy getPotentialEnergy(const System &system);

        System minimiseEnergy(System &system, double tolerance, int max_iteration);

        System equilibrateSystem(System &system, SireUnits::Dimension::Time equib_time_step,
                int equib_steps);

        void integrate(IntegratorWorkspace &workspace,
                const Symbol &nrg_component,
                SireUnits::Dimension::Time timestep,
                int nmoves, bool record_stats);

        IntegratorWorkspacePtr createWorkspace(const PropertyMap &map = PropertyMap()) const;
        IntegratorWorkspacePtr createWorkspace(const MoleculeGroup &molgroup, const PropertyMap &map = PropertyMap()) const;

        QString getIntegrator(void);
        void setIntegrator(QString);

        SireUnits::Dimension::Time getFriction(void);
        void setFriction(SireUnits::Dimension::Time);

        QString getCutoffType(void);
        void setCutoffType(QString);

        SireUnits::Dimension::Length getCutoff_distance(void);
        void setCutoff_distance(SireUnits::Dimension::Length);

        double getField_dielectric(void);
        void setField_dielectric(double);

        double getToleranceEwaldPME(void);
        void setToleranceEwaldPME(double);

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

        void setDeviceIndex(QString);
        QString getDeviceIndex(void);

        void setLJDispersion(bool);
        bool getLJDispersion(void);

        void setPrecision(QString);
        QString getPrecision(void);

        void setReinitialiseContext(bool);

        double getIntegration_tolerance(void);
        void setIntegration_tolerance(double tollerance);

        SireUnits::Dimension::Time getTimetoSkip(void);
        void setTimetoSkip(SireUnits::Dimension::Time);

//        void setMinimization(bool);
//
//        double getMinimiseTol(void);
//        void setMinimiseTol(double);
//
//        int getMinimiseIterations(void);
//        void setMinimiseIterations(int);
//
//        int getEquilib_iterations(void);
//        void setEquilib_iterations(int);

//        SireUnits::Dimension::Time getEquilib_time_step(void);
//        void setEquilib_time_step(SireUnits::Dimension::Time);


    private:
        void createContext(IntegratorWorkspace &workspace,
                SireUnits::Dimension::Time timestep);
        void destroyContext();

        /** Whether or not to save the velocities after every step, or to save them at the end of all of the steps */
        bool frequent_save_velocities;
        /** The Molecule Group on which the integrator operates */
        MolGroupPtr molgroup;
        /** Pointer to OpenMM context that describes the desired simulation*/
        //OpenMM::Context* context;

        /**Try instead to...keep a copy of OpenMM::System */
        OpenMM::System* openmm_system;

        OpenMM::Context* openmm_context;

        /** Whether the openmm system has been initialised*/
        bool isSystemInitialised;
        bool isContextInitialised;

        QString Integrator_type;
        SireUnits::Dimension::Time friction;

        QString CutoffType;
        SireUnits::Dimension::Length cutoff_distance;
        double field_dielectric;

        double tolerance_ewald_pme;

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

        QString device_index;

        bool LJ_dispersion;

        QString precision;

        bool reinetialise_context;

        double integration_tol;

        SireUnits::Dimension::Time timeskip;

        //bool minimise;

        //double minimise_tol;
        //int minimise_iterations;


        //int equilib_iterations;
        //SireUnits::Dimension::Time equilib_time_step;

        bool is_periodic;

    };


}

Q_DECLARE_METATYPE(SireMove::OpenMMMDIntegrator)

SIRE_EXPOSE_CLASS(SireMove::OpenMMMDIntegrator)

SIRE_END_HEADER

#else // SIRE_USE_OPENMM

        namespace SireMove {

    class OpenMMMDIntegrator {
    public:

        OpenMMMDIntegrator() {
        }

        ~OpenMMMDIntegrator() {
        }

        static const char* typeName() {
            return "SireMM::OpenMMMDIntegrator";
        }

    };

}

Q_DECLARE_METATYPE(SireMove::OpenMMMDIntegrator)

#endif // SIRE_USE_OPENMM

#endif
