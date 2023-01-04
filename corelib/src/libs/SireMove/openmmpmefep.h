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

#ifndef SIREMOVE_OPENMMPMEFEP_H
#define SIREMOVE_OPENMMPMEFEP_H

#include <cstdio>
#include <list>
#include <utility>

#include <boost/tuple/tuple.hpp>

#ifdef SIRE_USE_OPENMM
    #include <OpenMM.h>
#endif

#include "integrator.h"
#include "SireUnits/temperature.h"
#include "SireSystem/system.h"

SIRE_BEGIN_HEADER

#ifdef SIRE_USE_OPENMM

namespace SireMove {
    class OpenMMPMEFEP;
}

SIREMOVE_EXPORT QDataStream& operator<<(QDataStream&, const SireMove::OpenMMPMEFEP&);
SIREMOVE_EXPORT QDataStream& operator>>(QDataStream&, SireMove::OpenMMPMEFEP&);

namespace SireMove {
    using tmpl_str = const QString;

    /** This class implements the single topology free energy method using
	OpenMM.

        @author Julien Michel, Gaetano Calabro, Antonia Mey, Hannes H Loeffler
     */
    class SIREMOVE_EXPORT OpenMMPMEFEP
    : public SireBase::ConcreteProperty<OpenMMPMEFEP, Integrator> {
        friend QDataStream& ::operator<<(QDataStream&, const OpenMMPMEFEP&);
        friend QDataStream& ::operator>>(QDataStream&, OpenMMPMEFEP&);

    public:
        OpenMMPMEFEP(bool frequent_save_velocities = false);

        OpenMMPMEFEP(const MoleculeGroup &molecule_group,
		     const MoleculeGroup &solutes,
		     const MoleculeGroup &solute_hard,
		     const MoleculeGroup &solute_todummy,
		     const MoleculeGroup &solute_fromdummy,
		     bool frequent_save_velocities = false
		     );

        OpenMMPMEFEP(const OpenMMPMEFEP &other);

        ~OpenMMPMEFEP();

        OpenMMPMEFEP& operator=(const OpenMMPMEFEP &other);

        bool operator==(const OpenMMPMEFEP &other) const;
        bool operator!=(const OpenMMPMEFEP &other) const;

        static const char* typeName();

        QString toString() const;

        Ensemble ensemble() const;

        bool isTimeReversible() const;

        void initialise(bool fullPME = false);

        SireUnits::Dimension::MolarEnergy getPotentialEnergy(const System &system);

        System minimiseEnergy(System &system, double tolerance, int max_iteration);

        System annealSystemToLambda(System &system, SireUnits::Dimension::Time anneal_step_size,
                int annealing_steps);

        void integrate(IntegratorWorkspace &workspace,
                const Symbol &nrg_component,
                SireUnits::Dimension::Time timestep,
                int nmoves, bool record_stats);

        IntegratorWorkspacePtr createWorkspace(const PropertyMap &map = PropertyMap()) const;
        IntegratorWorkspacePtr createWorkspace(const MoleculeGroup &molgroup, const PropertyMap &map = PropertyMap()) const;

	QString getCombiningRules(void);
	void setCombiningRules(QString);

        QString getCutoffType(void);

        SireUnits::Dimension::Length getCutoffDistance(void);
        void setCutoffDistance(SireUnits::Dimension::Length);

        double getFieldDielectric(void);
        void setFieldDielectric(double);

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

        int getEnergyFrequency();
        void setEnergyFrequency(int);

        void setDeviceIndex(QString);
        QString getDeviceIndex(void);

        void setPrecision(QString);
        QString getPrecision(void);

        double getAlchemicalValue(void);
        void setAlchemicalValue(double);

        void setAlchemicalArray(QVector<double>);

        float getCoulombPower(void);
        void setCoulombPower(float);

        double getShiftDelta(void);
        void setShiftDelta(double);

        double getDeltaAlchemical(void);
        void setDeltatAlchemical(double);

        QVector<double> getGradients(void);
        QVector<double> getEnergies(void);

        QVector<double> getForwardMetropolis(void);
        QVector<double> getBackwardMetropolis(void);

        QVector<QVector <double> > getReducedPerturbedEnergies(void);

        QString getIntegrator(void);
        void setIntegrator(QString);

        SireUnits::Dimension::Time getFriction(void);

        void setFriction(SireUnits::Dimension::Time);

        double getIntegrationTolerance(void);
        void setIntegrationTolerance(double tollerance);

        SireUnits::Dimension::Time getTimetoSkip(void);
        void setTimetoSkip(SireUnits::Dimension::Time);

        void setReinitialiseContext(bool);

        int getRandomSeed(void);
        void setRandomSeed(int);

	void setDebug(bool);

    private:
        void createContext(IntegratorWorkspace &workspace,
                SireUnits::Dimension::Time timestep);
        void destroyContext();
        void updateBoxDimensions(OpenMM::State &state_openmm,
				 QVector<QVector<Vector>> &buffered_dimensions,
				 AtomicVelocityWorkspace &ws);

        double getPotentialEnergyAtLambda(double lambda);
        void updateOpenMMContextLambda(double lambda);
        boost::tuples::tuple<double, double, double> calculateGradient(double increment_plus,
        double increment_minus, double potential_energy_lambda, double beta);
        QVector<double> computeReducedPerturbedEnergies(double);
        void emptyContainers(void);

        void addAndersenThermostat(OpenMM::System &system);
        void addMCBarostat(OpenMM::System &system);

        /** Whether or not to save the velocities after every step, or to save them
	    at the end of all of the steps */
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

	QString combiningRules;
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

        double current_lambda;

        float coulomb_power;

        double shift_delta;

        double delta_alchemical;

        QVector<double> alchemical_array;

        QVector<double> finite_diff_gradients;

        QVector<double> pot_energies;

        QVector<double> forward_Metropolis;

        QVector<double> backward_Metropolis;

        QVector<QVector <double> > reduced_perturbed_energies;

        QVector<bool> perturbed_energies;

        QString Integrator_type;

        SireUnits::Dimension::Time friction;

        double integration_tol;

        SireUnits::Dimension::Time timeskip;

        bool reinitialise_context;

        bool Debug;

        int random_seed;

	static tmpl_str GENERAL;
	static tmpl_str GENERAL_SIGMA[2];

	static tmpl_str TODUMMY;
	static tmpl_str TODUMMY_SIGMA[2];

	static tmpl_str FROMDUMMY;
	static tmpl_str FROMDUMMY_SIGMA[2];

	static tmpl_str FROMTODUMMY;
	static tmpl_str FROMTODUMMY_SIGMA[2];

	static tmpl_str INTRA_14_CLJ;
	static tmpl_str INTRA_14_CLJ_SIGMA[2];

	static tmpl_str CORR_RECIP;

	/* "Light" atoms are defined to have a mass of HMASS or smaller.  This
           ensures that hydrogens in the HMR scheme will be constraint.  The
           specific value of 5.0 assumes that the HMR factor does not exceed
           4.0 and the heavy atom in CH3 or NH3 does have a minimum mass (see
           OpenMMMD.py) larger than HMASS.  In this way hydrogens and heavier
           atoms (assuming no elements between H and C) should be cleanly
           separated by mass. */
        const double HMASS = 5.0; // g/mol
        const double SMALL = 0.0001;
    };
}

Q_DECLARE_METATYPE(SireMove::OpenMMPMEFEP)

SIRE_EXPOSE_CLASS(SireMove::OpenMMPMEFEP)

SIRE_END_HEADER

#else // SIRE_USE_OPENMM

namespace SireMove {

    class OpenMMPMEFEP {
    public:

        OpenMMPMEFEP() {}

        ~OpenMMPMEFEP() {}

        static const char* typeName()
	{
            return "SireMM::OpenMMPMEFEP";
        }
    };
}

Q_DECLARE_METATYPE(SireMove::OpenMMPMEFEP)

#endif // SIRE_USE_OPENMM

#endif
