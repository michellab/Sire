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

SIREMOVE_EXPORT QDataStream& operator<<(QDataStream&, const SireMove::OpenMMMDIntegrator&);
SIREMOVE_EXPORT QDataStream& operator>>(QDataStream&, SireMove::OpenMMMDIntegrator&);

namespace SireMove {

    /** This class implements a "pure" MD integrator using OpenMM. 
     * No free energy methods are supported.

        @author Julien Michel and Gaetano Calabro
     */
    class SIREMOVE_EXPORT OpenMMMDIntegrator
    : public SireBase::ConcreteProperty<OpenMMMDIntegrator, Integrator> {
        friend SIREMOVE_EXPORT QDataStream& ::operator<<(QDataStream&, const OpenMMMDIntegrator&);
        friend SIREMOVE_EXPORT QDataStream& ::operator>>(QDataStream&, OpenMMMDIntegrator&);

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
        SireUnits::Dimension::MolarEnergy getKineticEnergy();

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

        SireUnits::Dimension::Length getCutoffDistance(void);
        void setCutoffDistance(SireUnits::Dimension::Length);

        double getFieldDielectric(void);
        void setFieldDielectric(double);

        double getToleranceEwaldPME(void);
        void setToleranceEwaldPME(double);

        bool getAndersen(void);
        void setAndersen(bool);

        double getAndersenFrequency(void);
        void setAndersenFrequency(double);

        bool getMCBarostat(void);
        void setMCBarostat(bool);

        void setMCBarostatFrequency(int);
        int getMCBarostatFrequency(void);

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

        int getCMMremovalFrequency(void);
        void setCMMremovalFrequency(int);

        int getBufferFrequency();
        void setBufferFrequency(int);

        void setDeviceIndex(QString);
        QString getDeviceIndex(void);

        void setLJDispersion(bool);
        bool getLJDispersion(void);

        void setPrecision(QString);
        QString getPrecision(void);

        void setReinitialiseContext(bool);

        double getIntegrationTolerance(void);
        void setIntegrationTolerance(double tollerance);

        SireUnits::Dimension::Time getTimetoSkip(void);
        void setTimetoSkip(SireUnits::Dimension::Time);

    private:
        void createContext(IntegratorWorkspace &workspace,
                SireUnits::Dimension::Time timestep);
        void destroyContext();
        void updateBoxDimensions(OpenMM::State &state_openmm, 
        QVector< Vector> &buffered_dimensions, bool Debug, 
        AtomicVelocityWorkspace &ws);

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

        bool is_periodic;
        
        SireUnits::Dimension::MolarEnergy openmmKineticEnergy;

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
