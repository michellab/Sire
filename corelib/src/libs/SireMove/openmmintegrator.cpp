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

#include "openmmintegrator.h"
#include "ensemble.h"

#include "SireMol/moleculegroup.h"
#include "SireMol/partialmolecule.h"
#include "SireMol/molecule.h"
#include "SireMol/atommasses.h"
#include "SireMol/atomcoords.h"
#include "SireMol/moleditor.h"

#include "SireMol/amberparameters.h"

#include "SireSystem/system.h"

#include "SireFF/forcetable.h"

#include "SireMaths/rangenerator.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireUnits/units.h"
#include "SireUnits/temperature.h"
#include "SireUnits/convert.h"

// ADDED BY JM 
#include "SireMol/connectivity.h"
#include "SireMol/bondid.h"
#include "SireMol/atomcharges.h"
#include "SireMM/internalff.h"
#include "SireIO/amber.h"
#include "SireMM/atomljs.h"

#include "SireVol/periodicbox.h"

#include "SireMove/flexibility.h"

//ADDED BY GAC
#include "SireMaths/vector.h"
#include "SireMol/mgname.h"
#include <iostream>
#include <QElapsedTimer>
#include <iomanip>
#include "fastio.h"
#include <QDebug>
#include <queue>

/* defines used by write_dcdstep */
#define NFILE_POS 8L
#define NSTEP_POS 20L

/* Define error codes that may be returned by the DCD routines */
#define DCD_SUCCESS      0  /* No problems                     */
#define DCD_EOF         -1  /* Normal EOF                      */
#define DCD_DNE         -2  /* DCD file does not exist         */
#define DCD_OPENFAILED  -3  /* Open of DCD file failed         */
#define DCD_BADREAD     -4  /* read call on DCD file failed    */
#define DCD_BADEOF      -5  /* premature EOF found in DCD file */
#define DCD_BADFORMAT   -6  /* format of DCD file is wrong     */
#define DCD_FILEEXISTS  -7  /* output file already exists      */
#define DCD_BADMALLOC   -8  /* malloc failed                   */
#define DCD_BADWRITE    -9  /* write call on DCD file failed   */

using namespace SireMove;
using namespace SireSystem;
using namespace SireMol;
using namespace SireFF;
using namespace SireCAS;
using namespace SireVol;
using namespace SireBase;
using namespace SireStream;
using namespace SireUnits;
using namespace SireUnits::Dimension;

//ADDED BY JM
using namespace SireMM;
using namespace SireIO;

//ADDED BY GAC
using namespace std;

typedef vector<pair<int,set<int> > > Vector_of_IntSet;
typedef vector<pair<int,int> > Vector_of_IntInt;

/* WRITE Macro to make porting easier */
#define WRITE(fd, buf, size) fio_fwrite(((void *) buf), (size), 1, (fd))


static int write_dcdheader(fio_fd fd, const char *remarks, int N, int ISTART, int NSAVC, double DELTA, int with_unitcell, int charmm);

static int write_dcdstep(fio_fd fd, int curframe, int curstep, int N, float *X, float *Y, float *Z, const double *unitcell, int charmm);

QString file_name(int i);

double gasdev(void);

void create_solute_solvent_lists(Vector_of_IntSet & solute_solvent_inter_lists, Vector_of_IntInt & solute_solute, Vector_of_IntInt & solvent_solvent, Vector_of_IntInt & solute_solvent); 

void create_intra_14_15_pairs(vector<pair<int,int> > & bond_pair, int N,  vector<pair<int,int> > & list_14, vector<pair<int,int> > & list_15); 

void BFS(vector<set<int> > & bonded, int root, int nbonds, set<int> & list);

void BFS_GE(vector<set<int> > & bonded, int root, int nbonds, set<int> & list);


enum {
	NOCUTOFF = 0,
	CUTOFFNONPERIODIC = 1,
	CUTOFFPERIODIC = 2
};

enum {

	NONE = 0,
	HBONDS = 1,
	ALLBONDS = 2,
	HANGLES = 3
	
};

class DCD {

	public:
		DCD(const char * file_name,std::vector<OpenMM::Vec3> & coordinates, std::vector<OpenMM::Vec3> & velocities,
			bool wrap_coord,int total_steps,int atoms,int frequency,int cutoff_type,int free_energy_nsamples,OpenMM::Context & context,
			OpenMM::VerletIntegrator & integrator,AtomicVelocityWorkspace & work_space ,std::vector<double> & center):
			positions_openmm(coordinates),velocities_openmm(velocities), context_openmm(context),integrator_openmm(integrator),ws(work_space),box_center(center)
			{
				filename=file_name;
				wrap=wrap_coord;
				md_steps=total_steps;
				nats=atoms;
				frequency_dcd=frequency;
				flag_cutoff=cutoff_type;
				free_energy_samples=free_energy_nsamples;
				fio_open(file_name,FIO_WRITE, &fd);

				X = new float[nats];
				Y = new float[nats];
				Z = new float[nats];
				
				
				for(int i=0; i<nats;i++){
					X[i]=0.0;
					Y[i]=0.0;
					Z[i]=0.0;

				}
				
				box_dims = new double[6];
				
				for(int i=0; i<5;i++){
					box_dims[i]=0.0;
				}


			}

	~DCD();

	void integrateMD(void);
	void integrateFreeEnergy(int current_frame);

	private:
		const char * filename;
		bool wrap;
		int md_steps;
		int nats;
		int frequency_dcd;
		int flag_cutoff;
		int free_energy_samples; 
		std::vector<OpenMM::Vec3> & positions_openmm;
		std::vector<OpenMM::Vec3> & velocities_openmm;
		OpenMM::Context & context_openmm;
		OpenMM::VerletIntegrator & integrator_openmm;
		fio_fd fd;
		float *X;
		float *Y;
		float *Z;
		double *box_dims;
		std::vector<double> box_center;
		AtomicVelocityWorkspace &ws;


};


static const RegisterMetaType<OpenMMIntegrator> r_openmmint;


/** Serialise to a binary datastream */
QDataStream SIREMOVE_EXPORT &operator<<(QDataStream &ds, const OpenMMIntegrator &velver)
{
    writeHeader(ds, r_openmmint, 1);
    
    SharedDataStream sds(ds);
    
    sds << velver.frequent_save_velocities << velver.CutoffType << velver.cutoff_distance << velver.field_dielectric
    	<< velver.Andersen_flag <<  velver.Andersen_frequency 
    	<< velver.MCBarostat_flag << velver.MCBarostat_frequency << velver.ConstraintType << velver.Pressure << velver.Temperature
    	<<velver.platform_type << velver.Alchemical_values << velver.Restraint_flag << velver.CMMremoval_frequency
    	<< static_cast<const Integrator&>(velver);
    
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOVE_EXPORT &operator>>(QDataStream &ds, OpenMMIntegrator &velver)
{
    VersionID v = readHeader(ds, r_openmmint);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> velver.frequent_save_velocities >> velver.CutoffType >> velver.cutoff_distance >> velver.field_dielectric
    		>> velver.Andersen_flag >>  velver.Andersen_frequency 
    		>> velver.MCBarostat_flag >> velver.MCBarostat_frequency >> velver.ConstraintType >> velver.Pressure >> velver.Temperature
    		>> velver.platform_type >> velver.Alchemical_values >> velver.Restraint_flag >> velver.CMMremoval_frequency
        	>> static_cast<Integrator&>(velver);
    }
    else
        throw version_error(v, "1", r_openmmint, CODELOC);
        
    return ds;
}

/** Constructor */
OpenMMIntegrator::OpenMMIntegrator(bool frequent_save) 
               : ConcreteProperty<OpenMMIntegrator,Integrator>(),
                 frequent_save_velocities(frequent_save), 
                 CutoffType("nocutoff"), cutoff_distance(1.0 * nanometer),field_dielectric(78.3),
                 Andersen_flag(false),Andersen_frequency(90.0), MCBarostat_flag(false),
                 MCBarostat_frequency(25),ConstraintType("none"),
                 Pressure(1.0 * bar),Temperature(300.0 * kelvin),platform_type("Reference"),Alchemical_values(),Restraint_flag(false),CMMremoval_frequency(0)
           
{}

/** Copy constructor */
OpenMMIntegrator::OpenMMIntegrator(const OpenMMIntegrator &other)
               : ConcreteProperty<OpenMMIntegrator,Integrator>(other),
                 frequent_save_velocities(other.frequent_save_velocities),
                 CutoffType(other.CutoffType),cutoff_distance(other.cutoff_distance),
                 field_dielectric(other.field_dielectric), Andersen_flag(other.Andersen_flag),
                 Andersen_frequency(other.Andersen_frequency), MCBarostat_flag(other.MCBarostat_flag),
                 MCBarostat_frequency(other.MCBarostat_frequency),ConstraintType(other.ConstraintType), 
                 Pressure(other.Pressure), Temperature(other.Temperature),platform_type(other.platform_type),
                 Alchemical_values(other.Alchemical_values),Restraint_flag(other.Restraint_flag),CMMremoval_frequency(other.CMMremoval_frequency)

{}

/** Destructor */
OpenMMIntegrator::~OpenMMIntegrator()
{}

/** Copy assignment operator */
OpenMMIntegrator& OpenMMIntegrator::operator=(const OpenMMIntegrator &other)
{
    Integrator::operator=(other);
    frequent_save_velocities = other.frequent_save_velocities;
    
    return *this;
}

/** Comparison operator */
bool OpenMMIntegrator::operator==(const OpenMMIntegrator &other) const
{
    return frequent_save_velocities == other.frequent_save_velocities and
           Integrator::operator==(other);
}

/** Comparison operator */
bool OpenMMIntegrator::operator!=(const OpenMMIntegrator &other) const
{
    return not OpenMMIntegrator::operator==(other);
}

/** Return a string representation of this integrator */
QString OpenMMIntegrator::toString() const
{
    return QObject::tr("OpenMMIntegrator()");
}
                                                       
/** Integrate the coordinates of the atoms in the molecules in 'molgroup'
    using the forces in 'forcetable', using the optionally supplied 
    property map to find the necessary molecular properties 
    
    \throw SireMol::missing_molecule
    \throw SireBase::missing_property
    \throw SireError:invalid_cast
    \throw SireError::incompatible_error
*/
void OpenMMIntegrator::integrate(IntegratorWorkspace &workspace, const Symbol &nrg_component, SireUnits::Dimension::Time timestep, int nmoves, bool record_stats) {
				       
	cout << "In OpenMMIntegrator::integrate()\n\n" ;
		
	
	// Initialise OpenMM

	// Convert the Sire data into data understood by OpenMM

	
	AtomicVelocityWorkspace &ws = workspace.asA<AtomicVelocityWorkspace>();
	
	const double dt = convertTo( timestep.value(), picosecond);
	
	cout << "time step = " << dt << " ps" << "\n";
	
	cout << "Number of moves = " << nmoves << "\n";

	const int nmols = ws.nMolecules();
    
	int nats = 0;
	
	int flag_cutoff;
	
	int flag_constraint;
	
	bool free_energy_calculation;
	
	Vector_of_IntSet solute_solvent_inter_lists;
	
	DCD * dcd_p=NULL;

		
	const System & ptr_sys = ws.system();
		
	cout << "Save frequency velocities = " << frequent_save_velocities << "\n";
	
	if (CutoffType == "nocutoff")
		flag_cutoff = NOCUTOFF;
	
	else if (CutoffType == "cutoffnonperiodic")
		flag_cutoff = CUTOFFNONPERIODIC;
	
	else if (CutoffType == "cutoffperiodic")
		flag_cutoff = CUTOFFPERIODIC;
	else
		throw SireError::program_bug(QObject::tr(
            "The CutOff method has not been specified. Possible choises: nocutoff, cutoffnonperiodic, cutoffperiodic"), CODELOC);
				

	if (ConstraintType == "none")
		flag_constraint = NONE;
	
	else if (ConstraintType == "hbonds")
		flag_constraint = HBONDS;
	
	else if(ConstraintType == "allbonds")
		flag_constraint = ALLBONDS;
		
	else if (ConstraintType == "hangles")
		flag_constraint = HANGLES;

	else
		throw SireError::program_bug(QObject::tr(
            "The Constraints method has not been specified. Possible choises: none, hbonds, allbonds, hangles"), CODELOC);
		
	
	cout << "\nConstraint Type = " << ConstraintType.toStdString() << "\n";
	
	
	// Free energy calculation flag
	
	if(Alchemical_values.size() == 0){
		free_energy_calculation = false;
		cout << "\nFree Energy calculation = OFF\n\n";
	}
	else{
		free_energy_calculation = true;
		cout << "\nFree Energy calculation = ON\n\n";
	}
	
	
	for (int i=0; i<nmols; ++i){
		nats = nats + ws.nAtoms(i);
	}

	cout << "There are " << nats << " atoms " << "There are " << nmols << " molecules" <<"\n" ;
	
	// INITIALIZE OpenMM
	
	const double Coulob14Scale = 1.0/1.2;
	
	const double LennardJones14Scale = 1.0/2.0;

	//Load Plugins from the OpenMM standard Plugin Directory 
	
	OpenMM::Platform::loadPluginsFromDirectory(OpenMM::Platform::getDefaultPluginsDirectory());
	
	//OpenMM system

	OpenMM::System system_openmm;
	
	//flag to extract information from the openmm system
	
	int infoMask = 0;

	infoMask = OpenMM::State::Positions;

	infoMask = infoMask + OpenMM::State::Velocities; 

	infoMask = infoMask +  OpenMM::State::Energy;
	

	
	//OpenMM non Bonded Forces
	

	OpenMM::NonbondedForce * nonbond_openmm = new OpenMM::NonbondedForce();
	
	
	nonbond_openmm->setUseDispersionCorrection(false);
	
	/*****************************************************************IMPORTANT********************************************************************/

	if(free_energy_calculation == false)
		system_openmm.addForce(nonbond_openmm);


	OpenMM::CustomNonbondedForce * custom_softcore_solute_solvent = NULL;
	
	
	OpenMM::CustomNonbondedForce * custom_solute_solute_solvent_solvent = NULL;
	
	
	OpenMM::CustomBondForce * custom_intra_14_15 = NULL;


	//Long Range interaction non bonded force method setting

	if(flag_cutoff == NOCUTOFF){
				
		nonbond_openmm->setNonbondedMethod(OpenMM::NonbondedForce::NoCutoff);
				
		if(free_energy_calculation == true){

			
			custom_softcore_solute_solvent = new OpenMM::CustomNonbondedForce("(10.0*Hls+100.0*Hcs)*ZeroOne;"
																			  "ZeroOne=issolute1*(1-issolute2)+issolute2*(1-issolute1);"
																			  "Hcs=((0.01*lam_cl)^(n+1))*138.935456*q_prod/sqrt(diff_cl+r^2);"
																			  "diff_cl=(1.0-lam_cl)*0.01;"
																			  "lam_cl=min(1,max(0,lambda-1));"
																			  "Hls=0.1*lam_lj*4.0*eps_avg*(LJ*LJ-LJ);"
																			  "LJ=((sigma_avg*sigma_avg)/soft)^3;"
																			  "soft=(diff_lj*delta*sigma_avg+r*r);"
																			  "diff_lj=(1.0-lam_lj)*0.1;"
																			  "lam_lj=max(0,min(1,lambda));"
																			  "q_prod=q1*q2;"
																			  "eps_avg=sqrt(eps1*eps2);"
																			  "sigma_avg=0.5*(sigma1+sigma2)");
																

			custom_softcore_solute_solvent->addGlobalParameter("lambda",Alchemical_values[0]);
			
			//int coulomb_Power = ptr_sys.property("coulombPower").toString().toInt();
			//double shift_Delta = ptr_sys.property("shiftDelta").toString().toDouble();
			
			int coulomb_Power = 0;
			double shift_Delta = 2.0;
			
			
			custom_softcore_solute_solvent->addGlobalParameter("delta",shift_Delta);
			custom_softcore_solute_solvent->addGlobalParameter("n",coulomb_Power);
			
			custom_softcore_solute_solvent->setNonbondedMethod(OpenMM::CustomNonbondedForce::NoCutoff);
			
			custom_solute_solute_solvent_solvent = new OpenMM::CustomNonbondedForce("(Hl+Hc)*ZeroOne;"
																					"ZeroOne=(1-issolute1)*(1-issolute2)+issolute1*issolute2;"
																					"Hl=4*eps_avg*((sigma_avg/r)^12-(sigma_avg/r)^6);"
																					"Hc=138.935456*q_prod/r;"
																					"q_prod=q1*q2;"
																					"eps_avg=sqrt(eps1*eps2);"
																					"sigma_avg=0.5*(sigma1+sigma2)");
			
			
			custom_solute_solute_solvent_solvent->setNonbondedMethod(OpenMM::CustomNonbondedForce::NoCutoff);
			
			
			custom_intra_14_15 = new OpenMM::CustomBondForce("Hl+Hc;"
													   		 "Hl=4*eps_avg*((sigma_avg/r)^12-(sigma_avg/r)^6);"
													   		 "Hc=138.935456*q_prod/r");
		
			
			cout << "Lambda = " << Alchemical_values[0] << " Coulomb Power = " << coulomb_Power << " Delta Shift = " << shift_Delta <<"\n";
			

		}

		cout << "\nCut off type = " << CutoffType.toStdString() << "\n";
	}


	if(flag_cutoff == CUTOFFNONPERIODIC || flag_cutoff == CUTOFFPERIODIC){

		const double converted_cutoff_distance = convertTo(cutoff_distance.value(), nanometer);
		
		if(flag_cutoff == CUTOFFNONPERIODIC)
			nonbond_openmm->setNonbondedMethod(OpenMM::NonbondedForce::CutoffNonPeriodic);

		else
			nonbond_openmm->setNonbondedMethod(OpenMM::NonbondedForce::CutoffPeriodic);


		nonbond_openmm->setCutoffDistance(converted_cutoff_distance);
		
		//Set Dielectric constant media
		nonbond_openmm->setReactionFieldDielectric(field_dielectric);
		
		//system_openmm.setDefaultPeriodicBoxVectors(OpenMM::Vec3(0,0,0),OpenMM::Vec3(0,0,0),OpenMM::Vec3(0,0,0));
		

		if(free_energy_calculation == true){

			/*custom_softcore_solute_solvent = new OpenMM::CustomNonbondedForce("(10.0*Hls+100.0*Hcs)*ZeroOne;"
																			  "ZeroOne=issolute1*(1-issolute2)+issolute2*(1-issolute1);"
																			  "Hcs=((0.01*lam_cl)^(n+1))*138.935456*q_prod/sqrt(diff_cl+r^2);"
																			  "diff_cl=(1.0-lam_cl)*0.01;"
																			  "lam_cl=min(1,max(0,lambda-1));"
																			  "Hls=0.1*lam_lj*4.0*eps_avg*(LJ*LJ-LJ);"
																			  "LJ=((sigma_avg*sigma_avg)/soft)^3;"
																			  "soft=(diff_lj*delta*sigma_avg+r*r);"
																			  "diff_lj=(1.0-lam_lj)*0.1;"
																			  "lam_lj=max(0,min(1,lambda));"
																			  "q_prod=q1*q2;"
																			  "eps_avg=sqrt(eps1*eps2);"
																			  "sigma_avg=0.5*(sigma1+sigma2)");*/
																			  
																			  
			custom_softcore_solute_solvent = new OpenMM::CustomNonbondedForce("(10.0*Hls+100.0*Hcs)*ZeroOne;"
																			  "ZeroOne=issolute1*(1-issolute2)+issolute2*(1-issolute1);"
																			  "Hcs=((0.01*lam_cl)^(n+1))*138.935456*q_prod*(1/sqrt(diff_cl+r*r) + krf*(diff_cl+r*r)-crf);"
																			  "diff_cl=(1.0-lam_cl)*0.01;"
																			  "lam_cl=min(1,max(0,lambda-1));"
																			  "Hls=0.1*lam_lj*4.0*eps_avg*(LJ*LJ-LJ);"
																			  "LJ=((sigma_avg*sigma_avg)/soft)^3;"
																			  "soft=(diff_lj*delta*sigma_avg+r*r);"
																			  "diff_lj=(1.0-lam_lj)*0.1;"
																			  "lam_lj=max(0,min(1,lambda));"
																			  "q_prod=q1*q2;"
																			  "eps_avg=sqrt(eps1*eps2);"
																			  "sigma_avg=0.5*(sigma1+sigma2)");


			custom_softcore_solute_solvent->addGlobalParameter("lambda",Alchemical_values[0]);

			int coulomb_Power = ptr_sys.property("coulombPower").toString().toInt();
			double shift_Delta = ptr_sys.property("shiftDelta").toString().toDouble();
			
			custom_softcore_solute_solvent->addGlobalParameter("delta",shift_Delta);
			custom_softcore_solute_solvent->addGlobalParameter("n",coulomb_Power);
			
			double eps2 = (field_dielectric - 1.0)/(2*field_dielectric+1.0);
			
			double kvalue = eps2/(converted_cutoff_distance * converted_cutoff_distance * converted_cutoff_distance);
			
			double cvalue = (1.0/converted_cutoff_distance)*(3.0*field_dielectric)/(2.0*field_dielectric+1.0);
			
			custom_softcore_solute_solvent->addGlobalParameter("krf",kvalue);

			custom_softcore_solute_solvent->addGlobalParameter("crf",cvalue);


			/*custom_solute_solute_solvent_solvent = new OpenMM::CustomNonbondedForce("(Hl+Hc)*ZeroOne;"
																					"ZeroOne=(1-issolute1)*(1-issolute2)+issolute1*issolute2;"
																					"Hl=4*eps_avg*((sigma_avg/r)^12-(sigma_avg/r)^6);"
																					"Hc=138.935456*q_prod/r;"
																					"q_prod=q1*q2;"
																					"eps_avg=sqrt(eps1*eps2);"
																					"sigma_avg=0.5*(sigma1+sigma2)");*/
																					
																					
			custom_solute_solute_solvent_solvent = new OpenMM::CustomNonbondedForce("(Hl+Hc)*ZeroOne;"
																					"ZeroOne=(1-issolute1)*(1-issolute2)+issolute1*issolute2;"
																					"Hl=4*eps_avg*((sigma_avg/r)^12-(sigma_avg/r)^6);"
																					"Hc=138.935456*q_prod*(1.0/r+(krf*r*r)-crf);"
																					"q_prod=q1*q2;"
																					"eps_avg=sqrt(eps1*eps2);"
																					"sigma_avg=0.5*(sigma1+sigma2)");


			
			custom_solute_solute_solvent_solvent->addGlobalParameter("krf",kvalue);

			custom_solute_solute_solvent_solvent->addGlobalParameter("crf",cvalue);


			
			if(flag_cutoff == CUTOFFNONPERIODIC){
				custom_softcore_solute_solvent->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffNonPeriodic);
				custom_solute_solute_solvent_solvent->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffNonPeriodic);
				
			}
			else{
				custom_softcore_solute_solvent->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffPeriodic);
				custom_solute_solute_solvent_solvent->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffPeriodic);
			}

			custom_softcore_solute_solvent->setCutoffDistance(converted_cutoff_distance);
			
			custom_solute_solute_solvent_solvent->setCutoffDistance(converted_cutoff_distance);


			custom_intra_14_15 = new OpenMM::CustomBondForce("withinCutoff*(Hl+Hc);"
															 "withinCutoff=step(cutoff-r);"
															 "Hl=4*eps_avg*((sigma_avg/r)^12-(sigma_avg/r)^6);"
														 	 "Hc=138.935456*q_prod/r");
																  
			custom_intra_14_15->addGlobalParameter("cutoff",converted_cutoff_distance);


			cout << "Lambda = " << Alchemical_values[0] << " Coulomb Power = " << coulomb_Power << " Delta Shift = " << shift_Delta <<"\n";
		
		}


		cout << "\nCut off type = " << CutoffType.toStdString() << "\n";
		cout << "CutOff distance = " << converted_cutoff_distance  << " Nm" << "\n";
		cout << "Dielectric constant= " << field_dielectric << "\n\n";
	
	}
	
	/*****************************************************************IMPORTANT********************************************************************/
	
	if(free_energy_calculation == true){
		system_openmm.addForce(custom_softcore_solute_solvent);
		system_openmm.addForce(custom_solute_solute_solvent_solvent);
		system_openmm.addForce(custom_intra_14_15);
	}


	//Andersen thermostat
	
	if (Andersen_flag == true){
	
		const double converted_Temperature = convertTo(Temperature.value(), kelvin);
		
		OpenMM::AndersenThermostat * thermostat = new OpenMM::AndersenThermostat(converted_Temperature, Andersen_frequency);
		system_openmm.addForce(thermostat);
		
		cout << "\nAndersen Thermostat set\n";
		cout << "Temperature = " << converted_Temperature << " K\n";
		cout << "Frequency collisions = " << Andersen_frequency << " 1/ps\n";
	
	}
	
	//Monte Carlo Barostat
	
	if (MCBarostat_flag == true) {
	
		const double converted_Temperature = convertTo(Temperature.value(), kelvin);
		const double converted_Pressure = convertTo(Pressure.value(), bar);
		
		
		OpenMM::MonteCarloBarostat * barostat = new OpenMM::MonteCarloBarostat(converted_Pressure, converted_Temperature, MCBarostat_frequency);
		system_openmm.addForce(barostat);
		
		cout << "\nMonte Carlo Barostat set\n";
		cout << "Temperature = " << converted_Temperature << " K\n";
		cout << "Pressure = " << converted_Pressure << " bar\n";
		cout << "Frequency every " << MCBarostat_frequency << " steps\n";
		
	}
	
	
	//OpenMM Bonded Forces

	OpenMM::HarmonicBondForce * bondStretch_openmm = new OpenMM::HarmonicBondForce();

	OpenMM::HarmonicAngleForce * bondBend_openmm = new OpenMM::HarmonicAngleForce();

	OpenMM::PeriodicTorsionForce * bondTorsion_openmm = new OpenMM::PeriodicTorsionForce();
	

	/*****************************************************************IMPORTANT********************************************************************/
	
	system_openmm.addForce(bondStretch_openmm);

	system_openmm.addForce(bondBend_openmm);

	system_openmm.addForce(bondTorsion_openmm);
	
	// Check whether positional restraints have been defined for a set of atoms in that molecule.
	// You can get the information out by getting the property and casting to VariantProperty
	//From VariantProperty you have the QVariant, so you can call .toDouble() and .toInt() there
	//so VariantProperty num = mol.property(QString("AtomNum(%1)").arg(i)).asA<VariantProperty>();
	//AtomNum atomnum( num.toInt() );
	//double x = mol.property(QString("x(%1)").arg(i)).asA<VariantProperty>().toDouble();
	//double y = ...; double z = ...;
	//QVector< QPair<AtomNum,Vector> > vals;
	//vals.append( QPair<AtomNum,Vector>(AtomNum(num.toInt()), Vector(x,y,z) ) );


	OpenMM::CustomExternalForce * positionalRestraints_openmm = NULL;

	if (Restraint_flag == true){
	
		positionalRestraints_openmm = new OpenMM::CustomExternalForce("k*( (x-xref)^2 + (y-yref)^2 + (z-zref)^2 )");
		
		positionalRestraints_openmm->addPerParticleParameter("xref");
		positionalRestraints_openmm->addPerParticleParameter("yref");
		positionalRestraints_openmm->addPerParticleParameter("zref");
		positionalRestraints_openmm->addPerParticleParameter("k");

		system_openmm.addForce(positionalRestraints_openmm);
		
		cout << "\nRestraint = ON\n\n";
	
	}


	std::vector<std::pair<int,int> > bondPairs;


	//OpenMM vector coordinate

	std::vector<OpenMM::Vec3> positions_openmm(nats);
	
	
	//OpenMM vector momenta
	
	std::vector<OpenMM::Vec3> velocities_openmm(nats);
	
	int system_index = 0;

	//cout << "INITAL COORDINATES AND VELOCITIES\n\n";
	
	QElapsedTimer timer_IN;
	
	timer_IN.start();
	
	MoleculeGroup molgroup = ws.moleculeGroup();
	
	
	// To avoid possible mismatch between the index in which atoms are added to the openmm system arrays and 
	// their atomic numbers in sire, one array is populated while filling up the openmm global arrays
	//  AtomNumtoopenmmIndex
	QHash<int, int> AtomNumToOpenMMIndex;

	 // Conversion factor because sire units of time are in AKMA, whereas OpenMM uses picoseconds

	double AKMAPerPs = 0.04888821;
	double PsPerAKMA = 1 / AKMAPerPs;

	for (int i=0; i < nmols; ++i){

		const int nats_mol = ws.nAtoms(i);

		Vector *c = ws.coordsArray(i);
		
		Vector *p = ws.momentaArray(i);
	
		const double *m = ws.massArray(i);
		
		//cout << "Molecule = " << i <<"\n";
		
		
		MolNum molnum = molgroup.molNumAt(i);
		const ViewsOfMol &molview = molgroup[molnum].data();
		const Molecule &mol = molview.molecule();
		Selector<Atom> molatoms = mol.atoms();
		
		
	/*****************************************************************IMPORTANT********************************************************************/

		bool boltz = true;//Generate velocities using Boltzmann Distribution (variance = 1 average = 0)
		
		
		for (int j=0; j < nats_mol; ++j){

			positions_openmm[system_index] = OpenMM::Vec3(c[j].x() * (OpenMM::NmPerAngstrom) ,
							c[j].y() * (OpenMM::NmPerAngstrom),c[j].z() * (OpenMM::NmPerAngstrom));

			if(boltz == false){
		
				velocities_openmm[system_index] = OpenMM::Vec3(p[j].x()/m[j] * (OpenMM::NmPerAngstrom) * PsPerAKMA,
														   p[j].y()/m[j] * (OpenMM::NmPerAngstrom) * PsPerAKMA,
														   p[j].z()/m[j] * (OpenMM::NmPerAngstrom) * PsPerAKMA);

			}
			else{

				double vx = gasdev();
				double vy = gasdev();
				double vz = gasdev();
				velocities_openmm[system_index] = OpenMM::Vec3(vx,vy,vz);

			}
			
			/*cout << "\natom = " << system_index  << " VELOCITY X = " << velocities_openmm[system_index][0] 
												 << " VELOCITY Y = " << velocities_openmm[system_index][1] 
												 << " VELOCITY Z = " << velocities_openmm[system_index][2];*/


			system_openmm.addParticle(m[j]) ;

			Atom at = molatoms.at(j);
			AtomNum atnum = at.number();

			//qDebug() << " openMM_index " << system_index << " Sire Atom Number " << atnum.value();

			AtomNumToOpenMMIndex[atnum.value()] = system_index;

			


			system_index = system_index + 1;
			
			
			/*cout << "\natom = " << system_index - 1  << " COORD X = " << c[j].x() 
							    					 << " COORD Y = " << c[j].y() 
							  						 << " COORD Z = " << c[j].z()<<"\n" ;

			cout << "atom = " << system_index - 1  << " MOMENTA X = " << p[j].x() 
												   << " MOMENTA Y = " << p[j].y() 
												   << " MOMENTA Z = " << p[j].z()<<"\n" ;

			cout << "atom = " << j << " MASS = " << m[j]<<"\n" ;*/


		}
		

	}


	int num_atoms_till_i = 0;
	
	
	MoleculeGroup solute_group;

	
	
	if(free_energy_calculation == true){

		custom_softcore_solute_solvent->addPerParticleParameter("q");
		custom_softcore_solute_solvent->addPerParticleParameter("sigma");
		custom_softcore_solute_solvent->addPerParticleParameter("eps");
		custom_softcore_solute_solvent->addPerParticleParameter("issolute");
		
		custom_solute_solute_solvent_solvent->addPerParticleParameter("q");
		custom_solute_solute_solvent_solvent->addPerParticleParameter("sigma");
		custom_solute_solute_solvent_solvent->addPerParticleParameter("eps");
		custom_solute_solute_solvent_solvent->addPerParticleParameter("issolute");
		
		
		custom_intra_14_15->addPerBondParameter("q_prod");
		custom_intra_14_15->addPerBondParameter("sigma_avg");
		custom_intra_14_15->addPerBondParameter("eps_avg");
		
		const QString solute = "solute";

		MGName Solute(solute);
		
		solute_group = ptr_sys[Solute];
	
	}

	int nions = 0;

	for (int i=0; i < nmols ; i++){
	
		const Vector *c = ws.coordsArray(i);
	
		Molecule molecule = molgroup.moleculeAt(i).molecule();
		
		int num_atoms_molecule = molecule.nAtoms();
		
		
		
		/*cout << "\nMOLECULE NUM = "<< i <<" N_ATOMS MOLECULE = " 
									 << num_atoms_molecule 
									 << " NUM ATOMS UNTILL CURRENT MOLUCULE = "
									 << num_atoms_till_i << "\n\n" ; */
		
	      
		//NON BONDED TERMS  
		
		//Lennard Jones
	  
		AtomLJs atomvdws = molecule.property("LJ").asA<AtomLJs>();
	
		AtomCharges atomcharges = molecule.property("charge").asA<AtomCharges>();
	
		QVector<SireMM::LJParameter> ljparameters = atomvdws.toVector();
	
		QVector<SireUnits::Dimension::Charge> charges = atomcharges.toVector();
		
	
		set<int> solvent;
		
		set<int> solute;
		
	
		for(int j=0; j< ljparameters.size();j++){
	
			double sigma = ljparameters[j].sigma();
		
			double epsilon = ljparameters[j].epsilon();
		
			double charge = charges[j].value();
			
			nonbond_openmm->addParticle(charge, sigma * OpenMM::NmPerAngstrom, epsilon * OpenMM::KJPerKcal);
			

			Atom atom = molecule.molecule().atoms()[j];
			
			AtomNum atnum = atom.number();
			
			vector<double> params(4);

			
			if(solute_group.contains(atom.molecule())){//solute atom
			
				if(free_energy_calculation == true){
	
					params[0]=charge;
					params[1]=sigma * OpenMM::NmPerAngstrom;
					params[2]=epsilon * OpenMM::KJPerKcal;
					params[3]=1.0;

					custom_softcore_solute_solvent->addParticle(params);
					custom_solute_solute_solvent_solvent->addParticle(params);
				}
				
				solute.insert(AtomNumToOpenMMIndex[atnum.value()]);

			}
			else{//solvent atom
				
				if(free_energy_calculation == true){

					params[0]=charge;
					params[1]=sigma * OpenMM::NmPerAngstrom;
					params[2]=epsilon * OpenMM::KJPerKcal;
					params[3]=0.0;

					custom_softcore_solute_solvent->addParticle(params);
					custom_solute_solute_solvent_solvent->addParticle(params);
				}
				
				
				solvent.insert(AtomNumToOpenMMIndex[atnum.value()]);
			
			}


			/*cout << "Sire Atom number = " << atnum.value() << " OpenMM index = " << AtomNumToOpenMMIndex[atnum.value()] << "\n";
			
			
			cout << "Is Solute = "<< (double) solute_group.contains(atom.molecule()) <<" J = "<< j << "\n";
			

			//qDebug()<< atomvdws.toString();

			cout << "sigma :" << sigma <<" A\n";
		
			cout << "epsilon :" << epsilon << " kcal/mol\n";

			cout << "charges : " << charge << " |e|\n";

			cout << "\n";*/
		}
		
		if(solute.size() != 0)
			solute_solvent_inter_lists.push_back(make_pair(1,solute));
		else
			solute_solvent_inter_lists.push_back(make_pair(0,solvent));


		//BONDED TERMS


		// The connectivity...
		bool hasConnectivity = molecule.hasProperty("connectivity");

		//Connectivity connectivity = molecule.property("connectivity").asA<Connectivity>();
		
		/* If the molecule does not have a connectivity, then we cannot get bonds/angles/dihedrals (effectively assuming it is a monoatomic molecule)*/

		if ( !hasConnectivity ){
			num_atoms_till_i = num_atoms_till_i + num_atoms_molecule ;
			/*cout << "\nAtoms = " <<  num_atoms_molecule << " Num atoms till i =" << num_atoms_till_i <<"\n";*/
			//cout << "\nION DETECTED\n";
			nions=nions+1;
			continue;
		}
		
		
		
		// The bonded parameters are stored in "amberparameters"
	
		AmberParameters amber_params = molecule.property("amberparameters").asA<AmberParameters>();
		
		//Bonds
		
		QList<BondID> bonds_ff = amber_params.getAllBonds();
		
		QVector<BondID> bonds = bonds_ff.toVector();


		for (int j=0; j < bonds_ff.length() ; j++){

			BondID bond_ff = bonds_ff[j];

			QList<double> bond_params = amber_params.getParams(bond_ff);

			double k = bond_params[0];
			double r0 = bond_params[1];
			
			int idx0 = bonds[j].atom0().asA<AtomIdx>().value();
			int idx1 = bonds[j].atom1().asA<AtomIdx>().value();
			
			//Select the atom type
			QString atom0 =  molecule.atom(AtomIdx(idx0)).toString();
			QString atom1 =  molecule.atom(AtomIdx(idx1)).toString();


			idx0 = idx0 + num_atoms_till_i;
			idx1 = idx1 + num_atoms_till_i;
			
			
			
			if(flag_constraint == NONE ){
				
				bondStretch_openmm->addBond(idx0,idx1,r0 * OpenMM::NmPerAngstrom, 
				k * 2.0 * OpenMM::KJPerKcal * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm); 
				
				//cout << "\nBOND ADDED TO "<< atom0.toStdString() << " AND " << atom1.toStdString() << "\n";
				
			}
			
			else if ( flag_constraint == ALLBONDS || flag_constraint == HANGLES ) {
			
				system_openmm.addConstraint(idx0,idx1, r0 *  OpenMM::NmPerAngstrom);
				
				//cout << "\nALLBONDS or HANGLES ADDED BOND CONSTRAINT TO " << atom0.toStdString() << " AND " << atom1.toStdString() << "\n";
				
			}
			
			else if ( flag_constraint == HBONDS ) {
			
				if ( (atom0[6] == 'H') || (atom1[6] == 'H') ){
					
					system_openmm.addConstraint(idx0,idx1, r0 *  OpenMM::NmPerAngstrom);
					
					//cout << "\nHBONDS ADDED BOND CONSTRAINT TO " << atom0.toStdString() << " AND " << atom1.toStdString() << "\n";
					
				}
					
				else {
					
					bondStretch_openmm->addBond(idx0,idx1,r0 * OpenMM::NmPerAngstrom, 
												k * 2.0 * OpenMM::KJPerKcal * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm);
					
					//cout << "\nHBONDS ADDED BOND TO " << atom0.toStdString() << " AND " << atom1.toStdString() << "\n";
				
				}	
			
			}
			
			
			
			
			//Bond exclusion List
			
			bondPairs.push_back(std::make_pair(idx0,idx1));

			/*cout << "Bond between atom global index " << idx0 << " and " << idx1 <<"\n"; 
			
			cout << "k_bond = " << k << " r0 = " << r0<<"\n";*/

			//cout << "Bond local" << bond_ff.toString() << " has k " << k << " r0 " << r0; 

			//cout << "\n";
		}

	
		//Angles
	
		QList<AngleID> angles_ff = amber_params.getAllAngles();
		
		QVector<AngleID> angles = angles_ff.toVector();
		
	
		for (int j=0; j < angles_ff.length() ; j++){
		
			AngleID angle_ff = angles_ff[j];

			QList<double> angle_params = amber_params.getParams(angle_ff);

			double k = angle_params[0];

			double theta0 = angle_params[1];// It is already in radiant

			int idx0 = angles[j].atom0().asA<AtomIdx>().value();

			int idx1 = angles[j].atom1().asA<AtomIdx>().value();

			int idx2= angles[j].atom2().asA<AtomIdx>().value();
			
			
			QString atom0 =  molecule.atom(AtomIdx(idx0)).toString();
			QString atom1 =  molecule.atom(AtomIdx(idx1)).toString();
			QString atom2 =  molecule.atom(AtomIdx(idx2)).toString();
			
			Vector diff = c[idx2] - c[idx0];
			
			/*cout << "\nATOM0 = " << atom0.toStdString() << " coords = [A] (" << c[idx0].x() << " ," << c[idx0].y() << " ,"<<  c[idx0].z()<< ")\n";
			cout << "\nATOM1 = " << atom1.toStdString() << " coords = [A] (" << c[idx1].x() << " ," << c[idx1].y() << " ,"<<  c[idx1].z()<< ")\n";
			cout << "\nATOM2 = " << atom2.toStdString() << " coords = [A] (" << c[idx2].x() << " ," << c[idx2].y() << " ,"<<  c[idx2].z()<< ")\n\n";*/
			 
		
			idx0 = idx0 + num_atoms_till_i;
			idx1 = idx1 + num_atoms_till_i;
			idx2 = idx2 + num_atoms_till_i;
			
			if ( flag_constraint == HANGLES ){
			
				if( ((atom0[6] == 'H') && (atom2[6] == 'H')) ){
								
					system_openmm.addConstraint(idx0,idx2, diff.length() *  OpenMM::NmPerAngstrom);
					
					/*cout << "Diff coords = [A] (" << diff.x() << " ," << diff.y() << " ,"<<  diff.z()<< ")\n";
					
					cout << "\nHANGLES H-X-H ADDED ANGLE CONSTRAINT TO " << atom0.toStdString() << " AND " 
														   << atom1.toStdString() << " AND " 
														   << atom2.toStdString() << " Distance = " << diff.length() << " A \n";*/
					
				}
				
				else if ( ((atom0[6] == 'H') && (atom1[6] == 'O')) || ((atom1[6] == 'O') && (atom2[6] == 'H')) ) {
				
					system_openmm.addConstraint(idx0,idx2, diff.length() *  OpenMM::NmPerAngstrom);
					
					/*cout << "\nHANGLES H-O-X or X-O-H ADDED ANGLE CONSTRAINT TO " << atom0.toStdString() << " AND " 
														   << atom1.toStdString() << " AND " 
														   << atom2.toStdString() << " Distance = " << diff.length() << " A \n";*/
				
				}
				
				else{
				
					bondBend_openmm->addAngle(idx0,idx1,idx2, theta0 , k * 2.0 * OpenMM::KJPerKcal);
					
					/*cout << "\nHANLGES ADDED ANGLE BOND TO " << atom0.toStdString() << " AND " 
													 << atom1.toStdString() << " AND "
													 << atom2.toStdString() << "\n";*/
				}
				
			}
			
			else {
				
				bondBend_openmm->addAngle(idx0,idx1,idx2, theta0 , k * 2.0 * OpenMM::KJPerKcal);
				
				/*cout << "\nADDED ANGLE BOND TO " << atom0.toStdString() << " AND " 
													 << atom1.toStdString() << " AND "
													 << atom2.toStdString() << "\n";*/
				
			}

			/*cout << "Angle between atom global index " << idx0 << " and " <<idx1 <<" and "<< idx2<<"\n";
			
			cout << "k_angle = " << k << " theta0 = " << theta0<<"\n"; */

			//cout << "Angle local " << angle_ff.toString() << " has k " << k << " theta0 " << theta0; 

			//cout << "\n";

		}

	
		//Dihedrals

		QList<DihedralID> dihedrals_ff = amber_params.getAllDihedrals();
		
		QVector<DihedralID> dihedrals = dihedrals_ff.toVector();
		

		for (int j=0; j < dihedrals_ff.length() ; j++ ){

			DihedralID dihedral_ff = dihedrals_ff[j];

			QList<double> dihedral_params = amber_params.getParams(dihedral_ff);
		
			int idx0 = dihedrals[j].atom0().asA<AtomIdx>().value() + num_atoms_till_i;
		
			int idx1 = dihedrals[j].atom1().asA<AtomIdx>().value() + num_atoms_till_i;
		
			int idx2 = dihedrals[j].atom2().asA<AtomIdx>().value() + num_atoms_till_i;
		
			int idx3 = dihedrals[j].atom3().asA<AtomIdx>().value() + num_atoms_till_i;

			// Variable number of parameters

			for (int k=0 ; k < dihedral_params.length() ; k = k + 3 ){
			
				double v = dihedral_params[ k ];
			
				int periodicity = dihedral_params[ k + 1 ];
			
				double phase = dihedral_params[ k + 2 ];
			
				bondTorsion_openmm->addTorsion(idx0, idx1, idx2, idx3, periodicity, phase , v * OpenMM::KJPerKcal);
			
				/*cout << "Dihedral between atom global index " << idx0 << " and " << idx1 << " and " << idx2 << " and " << idx3<<"\n";
				
				cout << "Amplitude_dih = " << v << " periodicity " << periodicity << " phase " << phase<<"\n";*/

				//cout << "Dihedral local" << dihedral_ff.toString() << " v " << v << " periodicity " << periodicity << " phase " << phase;
				
				//cout << "\n";
				
			}
		}


		//Improper Dihedrals

		QList<ImproperID> impropers_ff = amber_params.getAllImpropers();
		
		QVector<ImproperID> impropers = impropers_ff.toVector();
	

		for (int j=0; j < impropers_ff.length() ; j++ ){
		
			ImproperID improper_ff = impropers_ff[j];

			QList<double> improper_params = amber_params.getParams(improper_ff);

			// Variable number of parameters
			
			int idx0 = impropers[j].atom0().asA<AtomIdx>().value() + num_atoms_till_i;
		
			int idx1 = impropers[j].atom1().asA<AtomIdx>().value() + num_atoms_till_i;
		
			int idx2 = impropers[j].atom2().asA<AtomIdx>().value() + num_atoms_till_i;
		
			int idx3 = impropers[j].atom3().asA<AtomIdx>().value() + num_atoms_till_i;


			for (int k=0 ; k < improper_params.length() ; k = k + 3 ){
				
				double v = improper_params[ k ];
				
				double periodicity = improper_params[ k + 1 ];
				
				double phase = improper_params[ k + 2 ];
				
				bondTorsion_openmm->addTorsion(idx0, idx1, idx2, idx3, periodicity, phase , v * OpenMM::KJPerKcal);
				
				/*cout << "Improper Dihedral between atom global index " << idx0 << " and " << idx1 << " and " << idx2 << " and " << idx3<<"\n";
			
				cout << "Amplitude_impr = " << v << " periodicity " << periodicity << " phase " << phase <<"\n";*/
				
				//cout << "\n";
			}
		}
		
		if(Restraint_flag == true){
		
			bool hasRestrainedAtoms = molecule.hasProperty("restrainedatoms");
		
			if(hasRestrainedAtoms){

				Properties restrainedAtoms = molecule.property("restrainedatoms").asA<Properties>();

				int nrestrainedatoms = restrainedAtoms.property(QString("nrestrainedatoms")).asA<VariantProperty>().toInt();

				cout << " nrestrainedatoms " << nrestrainedatoms ;

				for (int i=0; i < nrestrainedatoms ; i++){

					int atomnum = restrainedAtoms.property(QString("AtomNum(%1)").arg(i)).asA<VariantProperty>().toInt();
				
					double xref = restrainedAtoms.property(QString("x(%1)").arg(i)).asA<VariantProperty>().toDouble();
					double yref = restrainedAtoms.property(QString("y(%1)").arg(i)).asA<VariantProperty>().toDouble();
					double zref = restrainedAtoms.property(QString("z(%1)").arg(i)).asA<VariantProperty>().toDouble();
					double k = restrainedAtoms.property(QString("k(%1)").arg(i)).asA<VariantProperty>().toDouble();

					int openmmindex = AtomNumToOpenMMIndex[atomnum];

					cout << " atomnum " << atomnum << " openmmindex " << openmmindex << " x " << xref << " y " << yref << " z " << zref << " k " << k;

					int posrestrdim = 4;

					std::vector<double> params(posrestrdim);

					params[0] = xref * OpenMM::NmPerAngstrom;
					params[1] = yref * OpenMM::NmPerAngstrom;
					params[2] = zref * OpenMM::NmPerAngstrom;
					params[3] = k  * ( OpenMM::KJPerKcal * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm );

					positionalRestraints_openmm->addParticle(openmmindex, params);

				}

			}
		}

		num_atoms_till_i = num_atoms_till_i + num_atoms_molecule ;


	}// end of loop over mols
	
	if(nions!=0)
		cout << "\n\nNumber of ions = " << nions << "\n\n";
	
	//Exclude the 1-2, 1-3 bonded atoms from nonbonded forces, and scale down 1-4 bonded atoms
	
	nonbond_openmm->createExceptionsFromBonds(bondPairs, Coulob14Scale, LennardJones14Scale);

	/*for(unsigned int i=0; i<bondPairs.size();i++){

		cout << "Bond excluded = ( " << bondPairs[i].first << " , " << bondPairs[i].second << " )";
		cout << "\n";
	}*/
	
	//Set Center of Mass motion removal 
	
	OpenMM::CMMotionRemover * cmmotionremover = NULL;
	
	if(CMMremoval_frequency > 0){
	
		cmmotionremover = new OpenMM::CMMotionRemover(CMMremoval_frequency);

		system_openmm.addForce(cmmotionremover);

		cout << "\n\nWill remove Center of Mass motion every " << CMMremoval_frequency << " steps\n\n";
	}


	
	if(free_energy_calculation == true){
	
		/*Vector_of_IntInt solute_solute;
		Vector_of_IntInt solvent_solvent;
		Vector_of_IntInt solute_solvent;
		
		create_solute_solvent_lists(solute_solvent_inter_lists, solute_solute, solvent_solvent, solute_solvent);
		
		
		cout << "Solute Solute list  SIZE = " << solute_solute.size() << "\n";
		cout << "Solvent Solvent list SIZE = " << solvent_solvent.size() << "\n";
		cout << "Solute Solvent list SIZE = " << solute_solvent.size() << "\n";


		cout << "\n\nSolute - Solute\n\n";
	
		for(unsigned int i=0;i<solute_solute.size();i++){
	
			cout << "( " << solute_solute[i].first << " , " << solute_solute[i].second << " )\n";

		}
		
		cout << "\n\nSolvent - Solvent\n\n";
	
		for(unsigned int i=0;i<solvent_solvent.size();i++){
	
			cout << "( " << solvent_solvent[i].first << " , " << solvent_solvent[i].second << " )\n";

		}
		
		cout << "\n\nSolute - Solvent\n\n";
		
		for(unsigned int i=0;i<solute_solvent.size();i++){
	
			cout << "( " << solute_solvent[i].first << " , " << solute_solvent[i].second << " )\n";

		}
		
		
		vector<pair<int,int> > list_14;
		
		vector<pair<int,int> > list_15;
		
		
		create_intra_14_15_pairs(bondPairs, nats,list_14, list_15);
		
	
		for(unsigned int i=0; i<list_14.size();i++){

			cout << "14 List = ( " << list_14[i].first << " , " <<list_14[i].second << " )";
			cout << "\n";
		}
		

		for(unsigned int i=0; i<list_15.size();i++){

			cout << "15 List = ( " << list_15[i].first << " , " <<list_15[i].second << " )";
			cout << "\n";
		}
	
		
		for(unsigned int i=0; i<list_14.size();i++){
		
			int p1= list_14[i].first;
			int p2= list_14[i].second;
			
			double  charge_p1=0;
			double  sigma_p1=0;
			double  epsilon_p1=0;
			
			double  charge_p2=0;
			double  sigma_p2=0;
			double  epsilon_p2=0;
		
			nonbond_openmm->getParticleParameters(p1,charge_p1,sigma_p1,epsilon_p1);
			
			nonbond_openmm->getParticleParameters(p2,charge_p2,sigma_p2,epsilon_p2);
		
			//cout << "p1 = " << p1 << " charge  p1 = " << charge_p1 << " sigma p1 = " << sigma_p1 << " epsilon p1 = " << epsilon_p1 <<"\n";
			//cout << "p2 = " << p2 << " charge  p2 = " << charge_p2 << " sigma p2 = " << sigma_p2 << " epsilon p2 = " << epsilon_p1 <<"\n\n";
			
			
			double tmp[]={Coulob14Scale*charge_p1*charge_p2,(sigma_p1+sigma_p2)*0.5, LennardJones14Scale*sqrt(epsilon_p1*epsilon_p2)};

			const std::vector<double> params(tmp,tmp+3);

			custom_intra_14_15->addBond(p1,p2,params);
			
			//nonbond_openmm->addException(p1,p2,0.0,1.0,0.0,true);
		
		}
		
		if(list_14.size() !=0 )
			cout << "\n\n14 list has been added\n\n";

		for(unsigned int i=0; i<list_15.size();i++){

			int p1= list_15[i].first;
			int p2= list_15[i].second;

			double  charge_p1=0;
			double  sigma_p1=0;
			double  epsilon_p1=0;

			double  charge_p2=0;
			double  sigma_p2=0;
			double  epsilon_p2=0;

			nonbond_openmm->getParticleParameters(p1,charge_p1,sigma_p1,epsilon_p1);
			
			nonbond_openmm->getParticleParameters(p2,charge_p2,sigma_p2,epsilon_p2);

			//cout << "p1 = " << p1 << " charge  p1 = " << charge_p1 << " sigma p1 = " << sigma_p1 << " epsilon p1 = " << epsilon_p1 <<"\n";
			//cout << "p2 = " << p2 << " charge  p2 = " << charge_p2 << " sigma p2 = " << sigma_p2 << " epsilon p2 = " << epsilon_p1 <<"\n\n";

			double tmp[]={charge_p1*charge_p2,(sigma_p1+sigma_p2)*0.5, sqrt(epsilon_p1*epsilon_p2)};

			const std::vector<double> params(tmp,tmp+3);

			custom_intra_14_15->addBond(p1,p2,params);
		
			nonbond_openmm->addException(p1,p2,0.0,1.0,0.0,true);

		}
		
		if(list_15.size() !=0 )
			cout << "\n\n15 list has been added\n\n";*/


		int num_exceptions = nonbond_openmm->getNumExceptions();
		
		//cout << "NUM EXCEPTIONS = " << num_exceptions << "\n";
		
		
		for(int i=0;i<num_exceptions;i++){

			int p1,p2;
			
			double charge_prod,sigma_avg,epsilon_avg;


			nonbond_openmm->getExceptionParameters(i,p1,p2,charge_prod,sigma_avg,epsilon_avg);

			/*cout << "Exception = " << i << " p1 = " << p1 << " p2 = " << p2 << " charge prod = " << charge_prod << 
				    " sigma avg = " << sigma_avg << " epsilon_avg = " << epsilon_avg << "\n";*/

			if(charge_prod!=0 && sigma_avg!=1 && epsilon_avg!=0){//1-4 interactions
			
				double tmp[]={charge_prod,sigma_avg, epsilon_avg};

				const std::vector<double> params(tmp,tmp+3);

				custom_intra_14_15->addBond(p1,p2,params);

			}

			custom_softcore_solute_solvent->addExclusion(p1,p2);
			custom_solute_solute_solvent_solvent->addExclusion(p1,p2);

		}

	}


	//OpenMM Integrator
	
	OpenMM::VerletIntegrator integrator_openmm(dt);//dt in pico seconds
	

	//OpenMM Context
	
	
	
	OpenMM::Platform& platform_openmm = OpenMM::Platform::getPlatformByName(platform_type.toStdString()); 

	OpenMM::Context context_openmm(system_openmm,integrator_openmm,platform_openmm);  
	
	
	if(flag_cutoff == CUTOFFPERIODIC){
		
		const System & ptr_sys = ws.system();
	
		const PropertyName &space_property = PropertyName("space");
	
		const PeriodicBox &space = ptr_sys.property(space_property).asA<PeriodicBox>();

		const double Box_x_Edge_Length = space.dimensions()[0] * OpenMM::NmPerAngstrom; //units in nm

		const double Box_y_Edge_Length = space.dimensions()[1] * OpenMM::NmPerAngstrom; //units in nm

		const double Box_z_Edge_Length = space.dimensions()[2] * OpenMM::NmPerAngstrom; //units in nm

		cout << "\nBOX SIZE [A] = (" << space.dimensions()[0] << " , " << space.dimensions()[1] << " ,  " << space.dimensions()[2] << ")\n\n";

		//Set Periodic Box Condition

		context_openmm.setPeriodicBoxVectors(OpenMM::Vec3(Box_x_Edge_Length,0,0),
											 OpenMM::Vec3(0,Box_y_Edge_Length,0),
											 OpenMM::Vec3(0,0,Box_z_Edge_Length));

	}

	//Add the coordinates  and velocities of the atoms to the OpenMM context

	context_openmm.setPositions(positions_openmm);  

	context_openmm.setVelocities(velocities_openmm);


	cout << "\nCOPY IN Simulation time = " << timer_IN.elapsed() / 1000.0 << " s"<<"\n\n";

	cout << "\n\nREMARK  Using OpenMM platform = " <<context_openmm.getPlatform().getName().c_str()<<"\n";

	
	OpenMM::State state_openmm;
	
	
	//Time benchmark
	
	QElapsedTimer timer_MD;
	
	QElapsedTimer timer_OUT;


	timer_MD.start();

	//MD SIMULATION FOR NMOVES

	//OpenMM::CMMotionRemover(1);


	if(free_energy_calculation == true){

		double delta = 0.001;

		const double beta = 1.0 / (0.0083144621 * convertTo(Temperature.value(), kelvin));

		cout << "\nBETA = " << beta <<" mol/kJ\n";

		int frequency_energy = 100; 

		int n_samples = nmoves/frequency_energy;
		
		int md_steps_per_sample = 100;

		cout << "Number Energy Samples = "<< n_samples << "\n\n";

		

		bool dcd = false;
		
		bool wrap = false;
		
		vector<double> center (3,0.0);  

		int frequency_dcd = 10;

		if(dcd == true)
			dcd_p = new DCD("dynamic.dcd", positions_openmm, velocities_openmm, wrap, md_steps_per_sample, nats, 
							frequency_dcd, flag_cutoff, n_samples, context_openmm, integrator_openmm,ws,center);

		for (int i=0; i<Alchemical_values.size();i++){

			double GF_acc = 0.0;

			double GB_acc = 0.0;

			double lam_val = 0.0;

			context_openmm.setParameter("lambda",Alchemical_values[i]);

			state_openmm=context_openmm.getState(infoMask);

			lam_val = context_openmm.getParameter("lambda");

			
			cout << "Start - Lambda = " << lam_val << " Potential energy lambda  = " << state_openmm.getPotentialEnergy() * OpenMM::KcalPerKJ << " [A + A^2] kcal/mol " << "\n";
	
			double j=0.0;



			while(j<n_samples){

				/*state_openmm=context_openmm.getState(infoMask);
				
				double potential_energy = state_openmm.getPotentialEnergy();

				cout << " Energy P  >>>>>>>>>>>>>>>>>>>> " << (potential_energy )  * OpenMM::KcalPerKJ << " kcal/mol" << "\n";*/
				


				if(dcd == true)
					dcd_p->integrateFreeEnergy(j);
				else
					integrator_openmm.step(md_steps_per_sample);



				state_openmm=context_openmm.getState(infoMask);

				cout<< "\nTotal Time = " << state_openmm.getTime() << " ps"<<"\n\n";


				lam_val = context_openmm.getParameter("lambda");

				double potential_energy_lambda = state_openmm.getPotentialEnergy();
				
				double potential_energy_lambda_plus_delta;
				
				double potential_energy_lambda_minus_delta;
				
				double plus;
				
				double minus;

				cout << "Lambda = " << lam_val << " Potential energy lambda MD = " << potential_energy_lambda  * OpenMM::KcalPerKJ << " kcal/mol" << "\n";

				if((Alchemical_values[i]+delta)>2.0){
					
					context_openmm.setParameter("lambda",Alchemical_values[i]-delta);
					
					state_openmm=context_openmm.getState(infoMask);

					potential_energy_lambda_minus_delta = state_openmm.getPotentialEnergy();
					
					minus =  exp(-beta * potential_energy_lambda_minus_delta) * exp(beta * potential_energy_lambda);
					
					plus = exp(beta * potential_energy_lambda_minus_delta) * exp(-beta*potential_energy_lambda);
					
					lam_val = context_openmm.getParameter("lambda");
					
					cout << "Lambda + delta > 2.0\n";
					
					cout << "Lambda - delta = " << lam_val << " Potential energy minus  = " << potential_energy_lambda_minus_delta * OpenMM::KcalPerKJ  << " kcal/mol" << "\n"; 
				}
				
				else if((Alchemical_values[i]-delta)<0.0){
					
					context_openmm.setParameter("lambda",Alchemical_values[i]+delta);
					
					state_openmm=context_openmm.getState(infoMask);
					
					potential_energy_lambda_plus_delta = state_openmm.getPotentialEnergy();
					
					plus = exp(-beta * potential_energy_lambda_plus_delta) * exp(beta * potential_energy_lambda);
					
					minus = exp(beta * potential_energy_lambda_plus_delta) * exp(-beta * potential_energy_lambda);
					
					lam_val = context_openmm.getParameter("lambda");
					
					cout << "Lambda + delta = " << lam_val << " Potential energy plus  = " << potential_energy_lambda_plus_delta * OpenMM::KcalPerKJ << " kcal/mol" << "\n";
				
					cout << "Lambda - delta < 0.0\n";
				}
				
				else{
				
					context_openmm.setParameter("lambda",Alchemical_values[i]+delta);

					state_openmm=context_openmm.getState(infoMask);

					potential_energy_lambda_plus_delta = state_openmm.getPotentialEnergy();

					lam_val = context_openmm.getParameter("lambda");

					cout << "Lambda + delta = " << lam_val << " Potential energy plus  = " << potential_energy_lambda_plus_delta * OpenMM::KcalPerKJ << " kcal/mol" << "\n";

					context_openmm.setParameter("lambda",Alchemical_values[i]-delta);

					state_openmm=context_openmm.getState(infoMask);

					potential_energy_lambda_minus_delta = state_openmm.getPotentialEnergy();

					plus = exp(-beta * potential_energy_lambda_plus_delta) * exp(beta * potential_energy_lambda);

					minus =  exp(-beta * potential_energy_lambda_minus_delta) * exp(beta * potential_energy_lambda);

					lam_val = context_openmm.getParameter("lambda");

					cout << "Lambda - delta = " << lam_val << " Potential energy minus  = " << potential_energy_lambda_minus_delta * OpenMM::KcalPerKJ  << " kcal/mol" << "\n"; 
				}


				GF_acc = GF_acc + plus;

				GB_acc = GB_acc + minus;


				cout << "\nplus = " << plus << " # minus = " << minus << "\n\n";

				//cout << "\nDifference +Delta= " << potential_energy_lambda_plus_delta - potential_energy_lambda  << " Difference -Delta=  " << potential_energy_lambda_minus_delta - potential_energy_lambda  << "\n\n";

				if(isnormal(GF_acc==0) || isnormal(GB_acc==0)){ 
					cout << "\n\n ********************** ERROR NAN!! ****************************\n\n";
					exit(-1); 
				}

				double avg_GF = GF_acc /(j+1.0);

				double avg_GB = GB_acc /(j+1.0);


				double Energy_GF = -(1.0/beta)*log(avg_GF);

				double Energy_GB = -(1.0/beta)*log(avg_GB);

				double Energy_Gradient_lamda = (Energy_GF - Energy_GB) / (2 * delta);

				cout << "\n\n*Energy Gradient = " << Energy_Gradient_lamda * OpenMM::KcalPerKJ << " kcal/(mol lambda)" << "\n\n";


				context_openmm.setParameter("lambda",Alchemical_values[i]);

				state_openmm=context_openmm.getState(infoMask);
				double dummy = state_openmm.getPotentialEnergy();
				cout << "\nDifference dummy = " << dummy - potential_energy_lambda<< "\n\n";
				
				j=j+1.0;
			}

		}

		state_openmm=context_openmm.getState(infoMask);

		positions_openmm = state_openmm.getPositions();
		
		velocities_openmm = state_openmm.getVelocities();

		cout<< "Total Time = " << state_openmm.getTime() << " ps"<<"\n\n";
		
		double kinetic_energy = 0.0; 
		double potential_energy = 0.0; 

		kinetic_energy = state_openmm.getKineticEnergy(); 

		potential_energy = state_openmm.getPotentialEnergy(); 

		
		cout<< "After MD" <<"\n";

		cout <<"Total Energy = " << (kinetic_energy + potential_energy) * OpenMM::KcalPerKJ << " Kcal/mol "
		 	 << " Kinetic Energy = " << kinetic_energy  * OpenMM::KcalPerKJ << " Kcal/mol " 
		 	 << " Potential Energy = " << potential_energy * OpenMM::KcalPerKJ << " Kcal/mol";

		cout<<"\n";

		cout << "\nMD Simulation time = " << timer_MD.elapsed() / 1000.0 << " s"<<"\n\n";


	}else{/******************** MD ***********************/

		timer_OUT.start();


		state_openmm=context_openmm.getState(infoMask);

		double kinetic_energy = 0.0; 
		double potential_energy = 0.0; 

		kinetic_energy = state_openmm.getKineticEnergy(); 

		potential_energy = state_openmm.getPotentialEnergy(); 

		cout<< "Before MD " <<"\n";

		cout <<"*Total Energy = " << (kinetic_energy + potential_energy) * OpenMM::KcalPerKJ << " Kcal/mol "
		 	 << " Kinetic Energy = " << kinetic_energy  * OpenMM::KcalPerKJ << " Kcal/mol " 
			 << " Potential Energy = " << potential_energy * OpenMM::KcalPerKJ << " Kcal/mol";

		cout<<"\n";

		int frequency_dcd = 10;//save every frequency_dcd

		bool dcd = false;
		
		bool wrap = false;
		
		vector<double> center;
		
		center.push_back(0.0);
		center.push_back(0.0);
		center.push_back(0.0); 
		
		if(dcd == true){
			dcd_p = new DCD("dynamic.dcd",positions_openmm,velocities_openmm,wrap,nmoves,nats,
							frequency_dcd,flag_cutoff,1,context_openmm,integrator_openmm,ws,center);
			
			dcd_p->integrateMD();

		}
		else
			integrator_openmm.step(nmoves);


		cout << "\nMD Simulation time = " << timer_MD.elapsed() / 1000.0 << " s"<<"\n\n";

		state_openmm=context_openmm.getState(infoMask);

		positions_openmm = state_openmm.getPositions();
		
		velocities_openmm = state_openmm.getVelocities();

		cout<< "Total Time = " << state_openmm.getTime() << " ps"<<"\n\n";

		kinetic_energy = state_openmm.getKineticEnergy(); 
	
		potential_energy = state_openmm.getPotentialEnergy(); 

		cout<< "After MD" <<"\n";
	
		cout <<"Total Energy = " << (kinetic_energy + potential_energy) * OpenMM::KcalPerKJ << " Kcal/mol "

		 	 << " Kinetic Energy = " << kinetic_energy  * OpenMM::KcalPerKJ << " Kcal/mol " 
		 	 << " Potential Energy = " << potential_energy * OpenMM::KcalPerKJ << " Kcal/mol";
	
		cout<<"\n";

		cout << "\nMD Simulation time = " << timer_MD.elapsed() / 1000.0 << " s"<<"\n\n";

	}

	timer_OUT.start();

	//Copy back to Sire positions and velocities 
	
	int k=0;


	for(int i=0; i<nmols;i++){

		Vector *sire_coords = ws.coordsArray(i);


		Vector *sire_momenta = ws.momentaArray(i);

		const double *m = ws.massArray(i);

		for(int j=0; j<ws.nAtoms(i);j++){
		
			sire_coords[j] = Vector(positions_openmm[j+k][0] * (OpenMM::AngstromsPerNm),
									positions_openmm[j+k][1] * (OpenMM::AngstromsPerNm),
									positions_openmm[j+k][2] * (OpenMM::AngstromsPerNm));

			sire_momenta[j] = Vector(velocities_openmm[j+k][0] * m[j] * (OpenMM::AngstromsPerNm) * AKMAPerPs,
									 velocities_openmm[j+k][1] * m[j] * (OpenMM::AngstromsPerNm) * AKMAPerPs,
									 velocities_openmm[j+k][2] * m[j] * (OpenMM::AngstromsPerNm) * AKMAPerPs);
		}

		k=k+ws.nAtoms(i);

	}

	ws.commitCoordinates();
	ws.commitVelocities();

	if (MCBarostat_flag == true) {

		OpenMM::Vec3 a;
		OpenMM::Vec3 b;
		OpenMM::Vec3 c;

		state_openmm.getPeriodicBoxVectors(a,b,c);

		cout << "\n\nNEW BOX DIMENSIONS [A] = (" << a[0] * OpenMM::AngstromsPerNm << ", " 
										   		 << b[1] * OpenMM::AngstromsPerNm << ", " 
										   		 << c[2] * OpenMM::AngstromsPerNm << ")\n\n";
										   		 

		Vector new_dims = Vector(a[0] * OpenMM::AngstromsPerNm, b[1] * OpenMM::AngstromsPerNm, c[2] * OpenMM::AngstromsPerNm);

		System & ptr_sys = ws.nonConstsystem();

		PeriodicBox  sp = ptr_sys.property("space").asA<PeriodicBox>();

		sp.setDimensions(new_dims);

		//cout << "\nBOX SIZE GAC = (" << sp.dimensions()[0] << " , " << sp.dimensions()[1] << " ,  " << sp.dimensions()[2] << ")\n\n";

		const QString string = "space" ;

		ptr_sys.setProperty(string,sp);

	}

	//cout << " record_stats " << record_stats ;

	if (record_stats)
		ws.collectStatistics();


	cout << "\n\nCOPY OUT Simulation time = " << timer_OUT.elapsed() / 1000.0 << " s"<<"\n\n";    

	/*cout << "FINAL COORDINATES AND VELOCITIES AFTER COMMIT" <<"\n\n";

	const int nmols_gac = ws.nMolecules();

	int system_index_gac = 0;

	for (int i=0; i < nmols_gac; ++i){

		const int nats_gac = ws.nAtoms(i);

		Vector *c_gac = ws.coordsArray(i);

		Vector *p_gac = ws.momentaArray(i);

		for (int j=0; j < nats_gac; ++j){

			cout <<"\n";
			cout << "atom = " << system_index_gac  << " COORD X = " << c_gac[j].x() << " COORD Y = " << c_gac[j].y() << " COORD Z = " << c_gac[j].z()<<"\n" ;
			cout << "atom = " << system_index_gac  << " MOMENTA X = " << p_gac[j].x() << " MOMENTA Y = " << p_gac[j].y() << " MOMENTA Z = " << p_gac[j].z()<<"\n" ;

			system_index_gac = system_index_gac + 1;

		}
		

	}*/

	
	
	cout << "\n---------------------------------------------\n";
	

	return;

}


/** Get the cufott type: nocutoff, cutoffnonperiodic, cutoffperiodic */
QString OpenMMIntegrator::getCutoffType(void){

	return CutoffType;

}

/** Set the cufott type: nocutoff, cutoffnonperiodic, cutoffperiodic */
void OpenMMIntegrator::setCutoffType(QString cutoff_type){

	CutoffType = cutoff_type;

}


/** Get the cutoff distance in A */
SireUnits::Dimension::Length OpenMMIntegrator::getCutoff_distance(void){

	return cutoff_distance;

}

/** Set the cutoff distance in A */
void OpenMMIntegrator::setCutoff_distance(SireUnits::Dimension::Length distance){

	cutoff_distance = distance;	

}

/** Get the dielectric constant */
double OpenMMIntegrator::getField_dielectric(void){

	return field_dielectric;
}

/** Set the dielectric constant */
void OpenMMIntegrator::setField_dielectric(double dielectric){
	
	field_dielectric=dielectric;

}

/** Set Andersen thermostat */

void OpenMMIntegrator::setAndersen(bool andersen){
	Andersen_flag = andersen;	
}

/** Get Andersen thermostat status on/off */
bool OpenMMIntegrator::getAndersen(void){
	
	return 	Andersen_flag;
	
}



/** Get the Andersen Thermostat frequency collision */
double OpenMMIntegrator::getAndersen_frequency(void){
	
	return Andersen_frequency;
	
}

/** Set the Andersen Thermostat frequency collision */
void OpenMMIntegrator::setAndersen_frequency(double freq){
	
	Andersen_frequency=freq;
	
}


/** Get the bath Temperature */
SireUnits::Dimension::Temperature OpenMMIntegrator::getTemperature(void){
		
	return Temperature;
}

/** Set the Temperature */
void OpenMMIntegrator::setTemperature(SireUnits::Dimension::Temperature temperature){
		
	Temperature = temperature;
}

/** Set Monte Carlo Barostat on/off */
void OpenMMIntegrator::setMCBarostat(bool MCBarostat){
	MCBarostat_flag = MCBarostat;
}

/** Get Andersen thermostat status on/off */
bool OpenMMIntegrator::getMCBarostat(void){
	
	return 	MCBarostat_flag;
	
}

/** Get the Monte Carlo Barostat frequency in time speps */
int OpenMMIntegrator::getMCBarostat_frequency(void){
	
	return  MCBarostat_frequency;
	
}

/** Set the Monte Carlo Barostat frequency in time speps */
void OpenMMIntegrator::setMCBarostat_frequency(int freq){
	
	MCBarostat_frequency=freq;
	
}

/** Get the Presure */
SireUnits::Dimension::Pressure OpenMMIntegrator::getPressure(void){
		
	return Pressure;
}

/** Set the Pressure */
void OpenMMIntegrator::setPressure(SireUnits::Dimension::Pressure pressure){

	Pressure = pressure;
}


/** Get the Constraint type: none, hbonds, allbonds, hangles */
QString OpenMMIntegrator::getConstraintType(void){

	return ConstraintType;

}

/** Set the Constraint type: none, hbonds, allbonds, hangles */
void OpenMMIntegrator::setConstraintType(QString constrain){

	ConstraintType = constrain;

}

/** Get the OpenMM Platform: CUDA, OpenCL, CPU */
QString OpenMMIntegrator::getPlatform(void){

	return platform_type;

}

/** Set the OpenMM Platform: CUDA, OpenCL, CPU */
void OpenMMIntegrator::setPlatform(QString platform){

	platform_type = platform;

}

/** Get the alchemical values used to calculate the free energy change via TI method*/
QVector<double> OpenMMIntegrator::getAlchemical_values(void){

	return Alchemical_values;

}

/** Set the alchemical values used to calculate the free energy change via TI method*/
void OpenMMIntegrator::setAlchemical_values(QVector<double> lambda_values){

	Alchemical_values = lambda_values;

}


/** Get the Restaint mode*/
bool OpenMMIntegrator::getRestraint(void){

	return Restraint_flag;

}

/** Set the Retraint mode */
void OpenMMIntegrator::setRestraint(bool Restraint){

	Restraint_flag = Restraint;
}

/** Get the Center of Mass motion removal frequency */
int OpenMMIntegrator::getCMMremoval_frequency(void){

	return CMMremoval_frequency;
}

/** Set the Center of Mass motion removal frequency */
void OpenMMIntegrator::setCMMremoval_frequency(int frequency){

	CMMremoval_frequency = frequency;
}


/** Create an empty workspace */
IntegratorWorkspacePtr OpenMMIntegrator::createWorkspace(const PropertyMap &map) const
{
	return IntegratorWorkspacePtr( new AtomicVelocityWorkspace(map) );
}

/** Return the ensemble of this integrator */
Ensemble OpenMMIntegrator::ensemble() const
{
	return Ensemble::NVE();
}

/** Return whether or not this integrator is time-reversible */
bool OpenMMIntegrator::isTimeReversible() const
{
	return true;
}

/** Create a workspace for this integrator for the molecule group 'molgroup' */
IntegratorWorkspacePtr OpenMMIntegrator::createWorkspace(const MoleculeGroup &molgroup, const PropertyMap &map) const
{
	return IntegratorWorkspacePtr( new AtomicVelocityWorkspace(molgroup,map) );
}

const char* OpenMMIntegrator::typeName()
{
	return QMetaType::typeName( qMetaTypeId<OpenMMIntegrator>() );
}




/*
 * Write a header for a new dcd file
 * Input: fd - file struct opened for binary writing
 *        remarks - string to be put in the remarks section of the header.  
 *                  The string will be truncated to 70 characters.
 *        natoms, istart, nsavc, delta - see comments in read_dcdheader
 * Output: 0 on success, negative error code on failure.
 * Side effects: Header information is written to the dcd file.
 */
static int write_dcdheader(fio_fd fd, const char *remarks, int N, 
                    int ISTART, int NSAVC, double DELTA, int with_unitcell,
                    int charmm) {
  int out_integer;
  float out_float;
  char title_string[200];
  time_t cur_time;
  struct tm *tmbuf;
  char time_str[81];

  out_integer = 84;
  WRITE(fd, (char *) & out_integer, sizeof(int));
  strcpy(title_string, "CORD");
  WRITE(fd, title_string, 4);
  fio_write_int32(fd, 0);      /* Number of frames in file, none written yet   */
  fio_write_int32(fd, ISTART); /* Starting timestep                            */
  fio_write_int32(fd, NSAVC);  /* Timesteps between frames written to the file */
  fio_write_int32(fd, 0);      /* Number of timesteps in simulation            */
  fio_write_int32(fd, 0);      /* NAMD writes NSTEP or ISTART - NSAVC here?    */
  fio_write_int32(fd, 0);
  fio_write_int32(fd, 0);
  fio_write_int32(fd, 0);
  fio_write_int32(fd, 0);
  if (charmm) {
    out_float = DELTA;
    WRITE(fd, (char *) &out_float, sizeof(float));
    if (with_unitcell) {
      fio_write_int32(fd, 1);
    } else {
      fio_write_int32(fd, 0);
    }
  } else {
    WRITE(fd, (char *) &DELTA, sizeof(double));
  }
  fio_write_int32(fd, 0);
  fio_write_int32(fd, 0);
  fio_write_int32(fd, 0);
  fio_write_int32(fd, 0);
  fio_write_int32(fd, 0);
  fio_write_int32(fd, 0);
  fio_write_int32(fd, 0);
  fio_write_int32(fd, 0);
  if (charmm) {
    fio_write_int32(fd, 24); /* Pretend to be Charmm version 24 */
  } else {
    fio_write_int32(fd, 0);
  }
  fio_write_int32(fd, 84);
  fio_write_int32(fd, 164);
  fio_write_int32(fd, 2);

  strncpy(title_string, remarks, 80);
  title_string[79] = '\0';
  WRITE(fd, title_string, 80);

  cur_time=time(NULL);
  tmbuf=localtime(&cur_time);
  strftime(time_str, 80, "REMARKS Created %d %B, %Y at %R", tmbuf);
  WRITE(fd, time_str, 80);

  fio_write_int32(fd, 164);
  fio_write_int32(fd, 4);
  fio_write_int32(fd, N);
  fio_write_int32(fd, 4);

  return DCD_SUCCESS;
}


/* 
 * Write a timestep to a dcd file
 * Input: fd - a file struct for which a dcd header has already been written
 *       curframe: Count of frames written to this file, starting with 1.
 *       curstep: Count of timesteps elapsed = istart + curframe * nsavc.
 *        natoms - number of elements in x, y, z arrays
 *        x, y, z: pointers to atom coordinates
 * Output: 0 on success, negative error code on failure.
 * Side effects: coordinates are written to the dcd file.
 */
static int write_dcdstep(fio_fd fd, int curframe, int curstep, int N, 
                  float *X, float *Y, float *Z, 
                  const double *unitcell, int charmm) {
  int out_integer;
  //int rc;

  if (charmm) {
    /* write out optional unit cell */
    if (unitcell != NULL) {
      out_integer = 48; /* 48 bytes (6 doubles) */
      fio_write_int32(fd, out_integer);
      WRITE(fd, unitcell, out_integer);
      fio_write_int32(fd, out_integer);
    }
  }

  /* write out coordinates */
  out_integer = N*4; /* N*4 bytes per X/Y/Z array (N floats per array) */
  fio_write_int32(fd, out_integer);
  if (fio_fwrite((void *) X, out_integer, 1, fd) != 1) return DCD_BADWRITE;
  fio_write_int32(fd, out_integer);
  fio_write_int32(fd, out_integer);
  if (fio_fwrite((void *) Y, out_integer, 1, fd) != 1) return DCD_BADWRITE;
  fio_write_int32(fd, out_integer);
  fio_write_int32(fd, out_integer);
  if (fio_fwrite((void *) Z, out_integer, 1, fd) != 1) return DCD_BADWRITE;
  fio_write_int32(fd, out_integer);

  /* update the DCD header information */
  fio_fseek(fd, NFILE_POS, FIO_SEEK_SET);
  fio_write_int32(fd, curframe);
  fio_fseek(fd, NSTEP_POS, FIO_SEEK_SET);
  fio_write_int32(fd, curstep);
  fio_fseek(fd, 0, FIO_SEEK_END);

  return DCD_SUCCESS;
}



QString file_name(int i){

	QString num = QString::number(i);

	QString str;

	if(i<10)
		str = "0000";
	if((i>=10) && (i<100))
		str = "000";
	if((i>=100) && (i<1000))
		str = "00";
	if((i>=1000) & (i<10000))
		str = "0";

	str = str + num;

	str= str + ".dcd\0";

	return str;

}

/*Box Muller algorithm */
double gasdev(void) {

	static bool available = false;

	static double gset;

	double fac, rsq, v1, v2;

	srand((unsigned) time ( NULL ));


	if (!available) {
		do {
			v1 = 2.0 * rand() / double(RAND_MAX) - 1.0;
			v2 = 2.0 * rand() / double(RAND_MAX) - 1.0;
			rsq = v1 * v1 + v2 * v2;
		} while (rsq >= 1.0 || rsq == 0.0);

		fac = sqrt(-2.0 * log(rsq) / rsq);
		gset = v1 * fac;
		available = true;
		return v2*fac;
	} 
	else {
		available = false;
		return gset;
	}
}




typedef  std::pair<int, int> PairInt;

bool compare(const PairInt &l,const PairInt &r)
{
  int lfirst = std::min(l.first,l.second);
  int rfirst = std::min(r.first,r.second);
  if (lfirst<rfirst) return true;
  if (rfirst<lfirst) return false;
  return std::max(l.first,l.second)<std::max(r.first,r.second);
}


void create_intra_14_15_pairs(vector<pair<int,int> > & bondPairs, int Natoms, vector<pair<int,int> > & list_pairs_14, vector<pair<int,int> > & list_pairs_15){

	typedef std::set<PairInt,bool (*)(const PairInt &,const PairInt &)> IntSet;
  	
  	IntSet list14_pairs(compare);

	IntSet list15_pairs(compare);
	
	
	vector<set<int> > bonded(Natoms);

	for (int i = 0; i < (int) bondPairs.size(); i++) {
		bonded[bondPairs[i].first].insert(bondPairs[i].second);
		bonded[bondPairs[i].second].insert(bondPairs[i].first);
	}

	/*for( int i=0;i<(int) bonded.size(); i++){
	
			cout << "\nParticle = " << i << "\n";
			
			for (set<int>::const_iterator iter = bonded[i].begin(); iter != bonded[i].end(); ++iter){
							
					cout << " " << (*iter);								
			}				
	}*/
	
	
	for(int i=0;i<Natoms;i++){
	
		set<int> list12;
	
		BFS(bonded, i, 1, list12);
	
		/*for (set<int>::const_iterator iter = list12.begin(); iter != list12.end(); ++iter){

			cout << "\n12 = " << (*iter) << "\n";
	
		}*/
	
		set<int> list13;
	
		BFS(bonded, i, 2, list13);
	
	
		for (set<int>::const_iterator iter = list13.begin(); iter != list13.end(); ++iter){
	
			set<int>::iterator it;

			it = list12.find(*iter);
		
			if(it != list12.end())
				list13.erase((*it));
		
		}
	
	
		/*for (set<int>::const_iterator iter = list13.begin(); iter != list13.end(); ++iter){

			cout << "\n13 = " << (*iter) << "\n";
	
		}*/
	
		set<int> list14;
	
		BFS(bonded, i, 3, list14);
	
	
		for (set<int>::const_iterator iter = list14.begin(); iter != list14.end(); ++iter){
	
			set<int>::iterator it;

			it = list12.find(*iter);
		
			if(it != list12.end())
				list14.erase((*it));
		}
	

		for (set<int>::const_iterator iter = list14.begin(); iter != list14.end(); ++iter){

			set<int>::iterator it;

			it = list13.find(*iter);
		
			if(it != list13.end())
				list14.erase((*it));
	
		}
	
	
		for (set<int>::const_iterator iter = list14.begin(); iter != list14.end(); ++iter){

			list14_pairs.insert(PairInt(i,*iter));
			
			
			//cout << "\n14 = " << (*iter) << "\n";
	
		}
	
	
		
	
	
		set<int> list15;
	
		BFS_GE(bonded, i, 4, list15);
	

	
    	for (set<int>::const_iterator iter = list15.begin(); iter != list15.end(); ++iter){
	
			set<int>::iterator it;

			it = list12.find(*iter);
		
			if(it != list12.end())
				list15.erase((*it));
		}
	

		for (set<int>::const_iterator iter = list15.begin(); iter != list15.end(); ++iter){

			set<int>::iterator it;

			it = list13.find(*iter);
		
			if(it != list13.end())
				list15.erase((*it));
		}

		for (set<int>::const_iterator iter = list15.begin(); iter != list15.end(); ++iter){

			set<int>::iterator it;

			it = list14.find(*iter);
		
			if(it != list14.end())
				list15.erase((*it));
	
		}
		
	
	
	
		for (set<int>::const_iterator iter = list15.begin(); iter != list15.end(); ++iter){

			
			list15_pairs.insert(PairInt(i,*iter));
			
			//cout << "\n15 = " << (*iter) << "\n";
	
		}
	
	
	}	
	
	
	for (IntSet::const_iterator i=list14_pairs.begin(); i!=list14_pairs.end(); ++i) { 
    	
    	list_pairs_14.push_back(make_pair((*i).first,(*i).second));
    	
    	//cout << "( " << (*i).first << "," << (*i).second << " )" << "\n";
  	
  	
  	} 
  	
  	
  	for (IntSet::const_iterator i=list15_pairs.begin(); i!=list15_pairs.end(); ++i) { 
  		
  		list_pairs_15.push_back(make_pair((*i).first,(*i).second));
    	
    	//cout << "( " << (*i).first << "," << (*i).second << " )" << "\n";
  	}
  	

	

}



void BFS(vector<set<int> > & bonded, int root, int nbonds, set<int> & list){

	vector<int> visited;
	
	vector<int> distance;
	
	queue<int> Q;

	for(unsigned int i=0; i< bonded.size(); i++){
		visited.push_back(0);
		distance.push_back(0);
	}
	
	visited[root] = 1;
	distance[root] = 0;
	
	Q.push(root);
	
	while(!Q.empty()){
	
		int u = Q.front();
		Q.pop();
		
		for(set<int>::const_iterator v = bonded[u].begin(); v != bonded[u].end(); ++v){
			
			if(visited[*v] == 0){
			
				visited[*v] = 1;
				
				distance[*v] = distance[u]+1;
				
				if(distance[*v] == nbonds)
					list.insert(*v);
				
				Q.push(*v);

			}
		
		}

	}
	
	/*cout << "\n\nNODE = " << root << "\n\n";
	
	for(unsigned int i=0; i< distance.size(); i++){
		cout << "distance (" << root << " , " << i << ") = "<< distance[i] << "\n";
	}*/
	
}


void BFS_GE(vector<set<int> > & bonded, int root, int nbonds, set<int> & list){

	vector<int> visited;
	
	vector<int> distance;
	
	queue<int> Q;

	for(unsigned int i=0; i< bonded.size(); i++){
		visited.push_back(0);
		distance.push_back(0);
	}
	
	visited[root] = 1;
	distance[root] = 0;
	
	Q.push(root);
	
	while(!Q.empty()){
	
		int u = Q.front();
		Q.pop();
		
		for(set<int>::const_iterator v = bonded[u].begin(); v != bonded[u].end(); ++v){
			
			if(visited[*v] == 0){
			
				visited[*v] = 1;
				
				distance[*v] = distance[u]+1;
				
				if(distance[*v] >= nbonds)
					list.insert(*v);
				
				Q.push(*v);
						
			}
		
		} 	

	}
	
	/*cout << "\n\nNODE = " << root << "\n\n";
	
	for(unsigned int i=0; i< distance.size(); i++){
		cout << "distance (" << root << " , " << i << ") = "<< distance[i] << "\n";
	}*/
	
}

void create_solute_solvent_lists(Vector_of_IntSet & solute_solvent_inter_lists, Vector_of_IntInt & solute_solute, Vector_of_IntInt & solvent_solvent, Vector_of_IntInt & solute_solvent){

	unsigned int Nmolecules = solute_solvent_inter_lists.size();
	
	for(unsigned int i=0; i<Nmolecules; i++){
	
		int type_i = solute_solvent_inter_lists[i].first;
		
		for(unsigned int j=i+1;j<Nmolecules;j++){

				int type_j = solute_solvent_inter_lists[j].first;
				
				if(type_i == type_j){

					if(type_i == 0){//solvent-solvent
					
						for (set<int>::const_iterator k = solute_solvent_inter_lists[i].second.begin(); k != solute_solvent_inter_lists[i].second.end(); ++k){
							
							for (set<int>::const_iterator l = solute_solvent_inter_lists[j].second.begin(); l != solute_solvent_inter_lists[j].second.end(); ++l){
							
									solvent_solvent.push_back(make_pair(*k,*l));

							}

						}

					}
					
					if(type_i == 1){//solute-solute
						
						for (set<int>::const_iterator k = solute_solvent_inter_lists[i].second.begin(); k != solute_solvent_inter_lists[i].second.end(); ++k){
							
							for (set<int>::const_iterator l = solute_solvent_inter_lists[j].second.begin(); l != solute_solvent_inter_lists[j].second.end(); ++l){
							
									solute_solute.push_back(make_pair(*k,*l));

							}

						}
					}

				}
				else{//solute-solvent
					
					for (set<int>::const_iterator k = solute_solvent_inter_lists[i].second.begin(); k != solute_solvent_inter_lists[i].second.end(); ++k){
							
							for (set<int>::const_iterator l = solute_solvent_inter_lists[j].second.begin(); l != solute_solvent_inter_lists[j].second.end(); ++l){
									
									solute_solvent.push_back(make_pair(*k,*l));

							}
						}
				}
				
		}
	
	}

}

void DCD::integrateMD(void){

	int cycles = 0;

	int steps;

	const int nmols = ws.nMolecules();
	
	MoleculeGroup molgroup = ws.moleculeGroup();

	double COG[3] = {0.0,0.0,0.0};//center of geometry
	
	double COT[3] = {0.0,0.0,0.0};//center of translation
	
	double dt = integrator_openmm.getStepSize();

	int infoMask = 0;

	infoMask = OpenMM::State::Positions;

	infoMask = infoMask + OpenMM::State::Velocities; 

	OpenMM::State state_openmm = context_openmm.getState(infoMask);


	OpenMM::Vec3 a;
	OpenMM::Vec3 b;
	OpenMM::Vec3 c;


	double delta = dt/0.0488821;

	cycles = md_steps/frequency_dcd;

	steps=frequency_dcd;

	int box = 0;

	if(flag_cutoff == CUTOFFPERIODIC){
		box=1;
	}

	write_dcdheader(fd, "Created by OpenMM", nats,0,frequency_dcd, delta, box,1);

	if(wrap == true){
		for(unsigned int i=0;i<positions_openmm.size();i++){

			COG[0] = COG[0] + positions_openmm[i][0]*(OpenMM::AngstromsPerNm); //X Cennter of Geometry
			COG[1] = COG[1] + positions_openmm[i][1]*(OpenMM::AngstromsPerNm); //Y Cennter of Geometry
			COG[2] = COG[2] + positions_openmm[i][2]*(OpenMM::AngstromsPerNm); //Z Cennter of Geometry
		}

		COG[0] = COG[0]/nats;
		COG[1] = COG[1]/nats;
		COG[2] = COG[2]/nats;
		
		COT[0]=box_center[0]-COG[0];
		COT[1]=box_center[1]-COG[1];
		COT[2]=box_center[2]-COG[2];


		

	}


	
	/*cout << "\n\nBox center[0] = " << box_center[0] << " Box center[1] = " << box_center[1] << " Box center[2] = " << box_center[2] << "\n\n";

	cout << "COG[0] = " << COG[0] << " COG[1] = " << COG[1] << " COG[2] = " << COG[2] << "\n\n";

	cout << "COT[0] = " << COT[0] << " COT[1] = " << COT[1] << " COT[2] = " << COT[2] << "\n\n";*/

	for(int i=0;i<cycles;i++){

		integrator_openmm.step(steps);/******************************IMPORTANT*****************************************/

		state_openmm=context_openmm.getState(infoMask);	

		positions_openmm = state_openmm.getPositions();

		velocities_openmm = state_openmm.getVelocities();


		if(flag_cutoff == CUTOFFPERIODIC){

			state_openmm.getPeriodicBoxVectors(a,b,c);

			box_dims[0]=a[0] * OpenMM::AngstromsPerNm;
			box_dims[1]=0.0;
			box_dims[2]=b[1] * OpenMM::AngstromsPerNm;
			box_dims[3]=0.0;
			box_dims[4]=0.0;
			box_dims[5]=c[2] * OpenMM::AngstromsPerNm;

		}

		int num_atoms_till_l=0;

		for (int l=0; l < nmols; ++l){

			const int nats_mol = ws.nAtoms(l);
			
			double molecule_COG[3] = {0.0,0.0,0.0};
			double T[3] = {0.0,0.0,0.0};//Translation vector

			if((wrap == true) && (flag_cutoff == CUTOFFPERIODIC)){

				for (int m=0; m < nats_mol; m++){
					
					int openmmindex = num_atoms_till_l+m;

					//cout << "Molecule index = " << l << " Atom index = " << m <<" Openmm Index = " << openmmindex <<"\n\n";

					molecule_COG[0] = molecule_COG[0] + positions_openmm[openmmindex][0]*OpenMM::AngstromsPerNm;;
					molecule_COG[1] = molecule_COG[1] + positions_openmm[openmmindex][1]*OpenMM::AngstromsPerNm;;
					molecule_COG[2] = molecule_COG[2] + positions_openmm[openmmindex][2]*OpenMM::AngstromsPerNm;;

				}

				molecule_COG[0] = molecule_COG[0]/nats_mol;
				molecule_COG[1] = molecule_COG[1]/nats_mol;
				molecule_COG[2] = molecule_COG[2]/nats_mol;
				
				//cout << "molecule_COG[0] = " << molecule_COG[0] << " molecule_COG[1] = " << molecule_COG[1] << " molecule_COG[2] = " << molecule_COG[2] << "\n\n";
				
				molecule_COG[0] = molecule_COG[0] + COT[0];
				molecule_COG[1] = molecule_COG[1] + COT[1];
				molecule_COG[2] = molecule_COG[2] + COT[2];
				
				//cout << "molecule_COG_T[0] = " << molecule_COG[0] << " molecule_COG_T[1] = " << molecule_COG[1] << " molecule_COG_T[2] = " << molecule_COG[2] << "\n\n";
				
				T[0] = molecule_COG[0] - box_dims[0]*round(molecule_COG[0]/box_dims[0]);
				T[1] = molecule_COG[1] - box_dims[2]*round(molecule_COG[1]/box_dims[2]);
				T[2] = molecule_COG[2] - box_dims[5]*round(molecule_COG[2]/box_dims[5]);
				
				T[0] = T[0] - molecule_COG[0];
				T[1] = T[1] - molecule_COG[1];
				T[2] = T[2] - molecule_COG[2];
				
				//cout << "T[0] = " << T[0] << " T[1] = " << T[1] << " T[2] = " << T[2] << "\n\n";

			}//end if wrap


			for (int m=0; m < nats_mol; m++){

				int openmmindex = num_atoms_till_l+m;

				X[openmmindex] = positions_openmm[openmmindex][0]*OpenMM::AngstromsPerNm;
				Y[openmmindex] = positions_openmm[openmmindex][1]*OpenMM::AngstromsPerNm;
				Z[openmmindex] = positions_openmm[openmmindex][2]*OpenMM::AngstromsPerNm;

				if((wrap == true) && (flag_cutoff == CUTOFFPERIODIC)){

					X[openmmindex] = X[openmmindex] + COT[0];
					Y[openmmindex] = Y[openmmindex] + COT[1];
					Z[openmmindex] = Z[openmmindex] + COT[2];
					
					/*cout << "X[" << openmmindex << "] =" << X[openmmindex] <<
						    " Y[" << openmmindex << "] = "<< Y[openmmindex] <<
						    " Z[" << openmmindex << "] = " << Z[openmmindex] << "\n";*/

					X[openmmindex] = X[openmmindex] + T[0];
					Y[openmmindex] = Y[openmmindex] + T[1];
					Z[openmmindex] = Z[openmmindex] + T[2];
					
					/*cout << "X[" << openmmindex << "] =" << X[openmmindex] <<
						    " Y[" << openmmindex << "] = "<< Y[openmmindex] <<
						    " Z[" << openmmindex << "] = " << Z[openmmindex] << "\n";*/

				}

			}

			num_atoms_till_l=num_atoms_till_l+nats_mol;
		
		}//end for molecules
		
		if(flag_cutoff == CUTOFFPERIODIC)
			write_dcdstep(fd, i+1, (i+1)*frequency_dcd, nats, X, Y, Z,box_dims, 1);

		else
			write_dcdstep(fd, i+1, (i+1)*frequency_dcd, nats, X, Y, Z,NULL, 1);

	}//end for cycles

	fio_fclose(fd);

	cout << "\nTrajectory file = " << filename << " has been written...." << "\n\n";

}


void DCD::integrateFreeEnergy(int current_frame){

	int cycles = 0;

	int steps;

	static double COG[3] = {0.0,0.0,0.0};
	
	static double COT[3] = {0.0,0.0,0.0};//center of translation

	double dt = integrator_openmm.getStepSize();
	
	const int nmols = ws.nMolecules();
	
	MoleculeGroup molgroup = ws.moleculeGroup();

	int infoMask = 0;

	infoMask = OpenMM::State::Positions;

	infoMask = infoMask + OpenMM::State::Velocities; 

	OpenMM::State state_openmm = context_openmm.getState(infoMask);


	OpenMM::Vec3 a;
	OpenMM::Vec3 b;
	OpenMM::Vec3 c;


	double delta = dt/0.0488821;

	cycles = md_steps/frequency_dcd;

	steps=frequency_dcd;

	int box = 0;

	if(flag_cutoff == CUTOFFPERIODIC){
		box=1;
	}

	if(current_frame == 0)
		write_dcdheader(fd, "Created by OpenMM", nats,0,frequency_dcd, delta, box,1);

	if(wrap == true && current_frame == 0){
		for(unsigned int i=0;i<positions_openmm.size();i++){

			COG[0] = COG[0] + positions_openmm[i][0]*(OpenMM::AngstromsPerNm); //X Cennter of Geometry
			COG[1] = COG[1] + positions_openmm[i][1]*(OpenMM::AngstromsPerNm); //Y Cennter of Geometry
			COG[2] = COG[2] + positions_openmm[i][2]*(OpenMM::AngstromsPerNm); //Z Cennter of Geometry
		}

		COG[0] = COG[0]/nats;
		COG[1] = COG[1]/nats;
		COG[2] = COG[2]/nats;
		
		COT[0]=box_center[0]-COG[0];
		COT[1]=box_center[1]-COG[1];
		COT[2]=box_center[2]-COG[2];
		
		/*cout << "\n\nBox center[0] = " << box_center[0] << " Box center[1] = " << box_center[1] << " Box center[2] = " << box_center[2] << "\n\n";

		cout << "COG[0] = " << COG[0] << " COG[1] = " << COG[1] << " COG[2] = " << COG[2] << "\n\n";

		cout << "COT[0] = " << COT[0] << " COT[1] = " << COT[1] << " COT[2] = " << COT[2] << "\n\n";*/

	}


	for(int i=0;i<cycles;i++){

		integrator_openmm.step(steps);/******************************IMPORTANT*****************************************/

		state_openmm=context_openmm.getState(infoMask);

		positions_openmm = state_openmm.getPositions();

		velocities_openmm = state_openmm.getVelocities();


		if(flag_cutoff == CUTOFFPERIODIC){

			state_openmm.getPeriodicBoxVectors(a,b,c);

			box_dims[0]=a[0] * OpenMM::AngstromsPerNm;
			box_dims[1]=0.0;
			box_dims[2]=b[1] * OpenMM::AngstromsPerNm;
			box_dims[3]=0.0;
			box_dims[4]=0.0;
			box_dims[5]=c[2] * OpenMM::AngstromsPerNm;

		}

		int num_atoms_till_l=0;

		for (int l=0; l < nmols; ++l){

			const int nats_mol = ws.nAtoms(l);
			
			double molecule_COG[3] = {0.0,0.0,0.0};
			double T[3] = {0.0,0.0,0.0};//Translation vector

			if((wrap == true) && (flag_cutoff == CUTOFFPERIODIC)){

				for (int m=0; m < nats_mol; m++){
					
					int openmmindex = num_atoms_till_l+m;

					//cout << "Molecule index = " << l << " Atom index = " << m <<" Openmm Index = " << openmmindex <<"\n\n";

					molecule_COG[0] = molecule_COG[0] + positions_openmm[openmmindex][0]*OpenMM::AngstromsPerNm;;
					molecule_COG[1] = molecule_COG[1] + positions_openmm[openmmindex][1]*OpenMM::AngstromsPerNm;;
					molecule_COG[2] = molecule_COG[2] + positions_openmm[openmmindex][2]*OpenMM::AngstromsPerNm;;

				}

				molecule_COG[0] = molecule_COG[0]/nats_mol;
				molecule_COG[1] = molecule_COG[1]/nats_mol;
				molecule_COG[2] = molecule_COG[2]/nats_mol;
				
				//cout << "molecule_COG[0] = " << molecule_COG[0] << " molecule_COG[1] = " << molecule_COG[1] << " molecule_COG[2] = " << molecule_COG[2] << "\n\n";
				
				molecule_COG[0] = molecule_COG[0] + COT[0];
				molecule_COG[1] = molecule_COG[1] + COT[1];
				molecule_COG[2] = molecule_COG[2] + COT[2];
				
				//cout << "molecule_COG_T[0] = " << molecule_COG[0] << " molecule_COG_T[1] = " << molecule_COG[1] << " molecule_COG_T[2] = " << molecule_COG[2] << "\n\n";
				
				T[0] = molecule_COG[0] - box_dims[0]*round(molecule_COG[0]/box_dims[0]);
				T[1] = molecule_COG[1] - box_dims[2]*round(molecule_COG[1]/box_dims[2]);
				T[2] = molecule_COG[2] - box_dims[5]*round(molecule_COG[2]/box_dims[5]);
				
				T[0] = T[0] - molecule_COG[0];
				T[1] = T[1] - molecule_COG[1];
				T[2] = T[2] - molecule_COG[2];
				
				//cout << "T[0] = " << T[0] << " T[1] = " << T[1] << " T[2] = " << T[2] << "\n\n";

			}//end if wrap


			for (int m=0; m < nats_mol; m++){

				int openmmindex = num_atoms_till_l+m;

				X[openmmindex] = positions_openmm[openmmindex][0]*OpenMM::AngstromsPerNm;
				Y[openmmindex] = positions_openmm[openmmindex][1]*OpenMM::AngstromsPerNm;
				Z[openmmindex] = positions_openmm[openmmindex][2]*OpenMM::AngstromsPerNm;

				if((wrap == true) && (flag_cutoff == CUTOFFPERIODIC)){

					X[openmmindex] = X[openmmindex] + COT[0];
					Y[openmmindex] = Y[openmmindex] + COT[1];
					Z[openmmindex] = Z[openmmindex] + COT[2];
					
					/*cout << "X[" << openmmindex << "] =" << X[openmmindex] <<
						    " Y[" << openmmindex << "] = "<< Y[openmmindex] <<
						    " Z[" << openmmindex << "] = " << Z[openmmindex] << "\n";*/

					X[openmmindex] = X[openmmindex] + T[0];
					Y[openmmindex] = Y[openmmindex] + T[1];
					Z[openmmindex] = Z[openmmindex] + T[2];
					
					/*cout << "X[" << openmmindex << "] =" << X[openmmindex] <<
						    " Y[" << openmmindex << "] = "<< Y[openmmindex] <<
						    " Z[" << openmmindex << "] = " << Z[openmmindex] << "\n";*/

				}

			}

			num_atoms_till_l=num_atoms_till_l+nats_mol;
		
		}//end for molecules
		
		if(flag_cutoff == CUTOFFPERIODIC)
			write_dcdstep(fd, current_frame*cycles + (i+1), current_frame*cycles + (i+1)*frequency_dcd, nats, X, Y, Z,box_dims, 1);
		else
			write_dcdstep(fd, current_frame*cycles + (i+1), current_frame*cycles + (i+1)*frequency_dcd, nats, X, Y, Z,NULL, 1);

	}//end for cycles

	
	if(current_frame == free_energy_samples-1){
		fio_fclose(fd);
		cout << "\nTrajectory file = " << filename << " has been written...." << "\n\n";
	}

}
