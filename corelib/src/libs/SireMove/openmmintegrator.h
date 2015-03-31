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

#ifndef SIREMOVE_OPENMMINTEGRATOR_H
#define SIREMOVE_OPENMMINTEGRATOR_H

#include "integrator.h"

#ifdef SIRE_USE_OPENMM
  #include <OpenMM.h>   // CONDITIONAL_INCLUDE
#endif

#include <cstdio>
#include "SireUnits/temperature.h"

SIRE_BEGIN_HEADER

#ifdef SIRE_USE_OPENMM

namespace SireMove
{
class OpenMMIntegrator;
}

QDataStream& operator<<(QDataStream&, const SireMove::OpenMMIntegrator&);
QDataStream& operator>>(QDataStream&, SireMove::OpenMMIntegrator&);

namespace SireMove
{

/** This class implements an MD integrator using OpenMM 
 
    @author Julien Michel and Gaetano Calabro
*/
class SIREMOVE_EXPORT OpenMMIntegrator
          : public SireBase::ConcreteProperty<OpenMMIntegrator,Integrator>
{

friend QDataStream& ::operator<<(QDataStream&, const OpenMMIntegrator&);
friend QDataStream& ::operator>>(QDataStream&, OpenMMIntegrator&);

public:
    OpenMMIntegrator(bool frequent_save_velocities = false);
    
    OpenMMIntegrator(const OpenMMIntegrator &other);
    
    ~OpenMMIntegrator();
    
    OpenMMIntegrator& operator=(const OpenMMIntegrator &other);
    
    bool operator==(const OpenMMIntegrator &other) const;
    bool operator!=(const OpenMMIntegrator &other) const;
    
    static const char* typeName();
    
    QString toString() const;
    
    Ensemble ensemble() const;
    
    bool isTimeReversible() const;
    
    void integrate(IntegratorWorkspace &workspace,
                   const Symbol &nrg_component, 
                   SireUnits::Dimension::Time timestep,
                   int nmoves, bool record_stats);

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
	
	QVector<double> getAlchemical_values(void);
	void setAlchemical_values(QVector<double>);
	
	bool getRestraint(void);
	void setRestraint(bool);
	
	int getCMMremoval_frequency(void);
	void setCMMremoval_frequency(int);


private:
	/** Whether or not to save the velocities after every step, or to save them at the end of all of the steps */

	bool frequent_save_velocities;

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

	QVector<double> Alchemical_values;
	
	bool Restraint_flag;
	
	int CMMremoval_frequency;

};


}

Q_DECLARE_METATYPE( SireMove::OpenMMIntegrator )

SIRE_EXPOSE_CLASS( SireMove::OpenMMIntegrator )

SIRE_END_HEADER

#else // SIRE_USE_OPENMM

namespace SireMove
{

    class OpenMMIntegrator
    {
    public:
        OpenMMIntegrator(){}
        ~OpenMMIntegrator(){}

        static const char* typeName(){ return "SireMM::OpenMMIntegrator"; }

    };

}

Q_DECLARE_METATYPE( SireMove::OpenMMIntegrator )

#endif // SIRE_USE_OPENMM

#endif
