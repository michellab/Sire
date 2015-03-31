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

#include "openmmacceleratedmdintegrator.h"
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

#include "SireMaths/constants.h"

//ADDED BY GAC
#include "SireMaths/vector.h"
#include "SireMol/mgname.h"
#include <iostream>
/* include <QElapsedTimer> */
#include <QDebug>
#include <QTime>

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
//using namespace SireMaths;

//ADDED BY GAC
using namespace std;


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

static const RegisterMetaType<OpenMMAMDIntegrator> r_openmmint;


/** Serialise to a binary datastream */
QDataStream SIREMOVE_EXPORT &operator<<(QDataStream &ds, const OpenMMAMDIntegrator &velver)
{
	writeHeader(ds, r_openmmint, 1);

	SharedDataStream sds(ds);

	sds << velver.frequent_save_velocities << velver.molgroup 
	<< velver.Integrator_type << velver.friction
	<< velver.CutoffType << velver.cutoff_distance << velver.field_dielectric
	<< velver.Andersen_flag <<  velver.Andersen_frequency 
	<< velver.MCBarostat_flag << velver.MCBarostat_frequency << velver.ConstraintType << velver.Pressure << velver.Temperature
	<<velver.platform_type << velver.Restraint_flag << velver.CMMremoval_frequency << velver.buffer_frequency
	<< velver.device_index
	<< velver.boosted_torsions << velver.etorsion << velver.alphatorsion 
	<< static_cast<const Integrator&>(velver);

	// Free OpenMM pointers??

	return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOVE_EXPORT &operator>>(QDataStream &ds, OpenMMAMDIntegrator &velver)
{
	VersionID v = readHeader(ds, r_openmmint);

	if (v == 1){

		SharedDataStream sds(ds);

		sds >> velver.frequent_save_velocities >> velver.molgroup 
		>> velver.Integrator_type >> velver.friction
		>> velver.CutoffType >> velver.cutoff_distance >> velver.field_dielectric
		>> velver.Andersen_flag >>  velver.Andersen_frequency 
		>> velver.MCBarostat_flag >> velver.MCBarostat_frequency >> velver.ConstraintType >> velver.Pressure >> velver.Temperature
		>> velver.platform_type >> velver.Restraint_flag >> velver.CMMremoval_frequency >> velver.buffer_frequency
		>> velver.device_index 
		>> velver.boosted_torsions >> velver.etorsion >> velver.alphatorsion
		>> static_cast<Integrator&>(velver);

	// Maybe....need to reinitialise from molgroup because openmm system was not serialised...
	velver.isInitialised = false;
	
	//qDebug() << " Re-initialisation of openmmmdintegrator from datastream";

	velver.initialise();
	}
	else
		throw version_error(v, "1", r_openmmint, CODELOC);

	return ds;
}

/** Constructor*/
OpenMMAMDIntegrator::OpenMMAMDIntegrator(bool frequent_save) : ConcreteProperty<OpenMMAMDIntegrator,Integrator>(),
					frequent_save_velocities(frequent_save), molgroup(MoleculeGroup()),
					openmm_system(0),isInitialised(false), Integrator_type("leapfrogverlet"),friction(1.0 / picosecond ),
					CutoffType("nocutoff"), cutoff_distance(1.0 * nanometer),field_dielectric(78.3),
					Andersen_flag(false),Andersen_frequency(90.0), MCBarostat_flag(false),
					MCBarostat_frequency(25),ConstraintType("none"),
					Pressure(1.0 * bar),Temperature(300.0 * kelvin),platform_type("Reference"),Restraint_flag(false),
					CMMremoval_frequency(0), buffer_frequency(0),device_index("0"),
					boosted_torsions(QList<DihedralID>()), etorsion( 0*kcal_per_mol ), alphatorsion(99999.0)
{}

/** Constructor using the passed molecule group */
OpenMMAMDIntegrator::OpenMMAMDIntegrator(const MoleculeGroup &molecule_group, bool frequent_save) 
					: ConcreteProperty<OpenMMAMDIntegrator,Integrator>(),
					frequent_save_velocities(frequent_save), 
					molgroup(molecule_group),
					openmm_system(0),isInitialised(false),
					Integrator_type("leapfrogverlet"),friction(1.0 / picosecond ),
					CutoffType("nocutoff"), cutoff_distance(1.0 * nanometer),field_dielectric(78.3),
					Andersen_flag(false),Andersen_frequency(90.0), MCBarostat_flag(false),
					MCBarostat_frequency(25),ConstraintType("none"),
					Pressure(1.0 * bar),Temperature(300.0 * kelvin),platform_type("Reference"),Restraint_flag(false),
					CMMremoval_frequency(0), buffer_frequency(0),device_index("0"),
					boosted_torsions(QList<DihedralID>()), etorsion( 0*kcal_per_mol ), alphatorsion(99999.0)
{}

/** Copy constructor */
OpenMMAMDIntegrator::OpenMMAMDIntegrator(const OpenMMAMDIntegrator &other)
					: ConcreteProperty<OpenMMAMDIntegrator,Integrator>(other),
					frequent_save_velocities(other.frequent_save_velocities),
					molgroup(other.molgroup),
					openmm_system(other.openmm_system),isInitialised(other.isInitialised),
					Integrator_type(other.Integrator_type),friction(other.friction),
					CutoffType(other.CutoffType),cutoff_distance(other.cutoff_distance),
					field_dielectric(other.field_dielectric), Andersen_flag(other.Andersen_flag),
					Andersen_frequency(other.Andersen_frequency), MCBarostat_flag(other.MCBarostat_flag),
					MCBarostat_frequency(other.MCBarostat_frequency),ConstraintType(other.ConstraintType), 
					Pressure(other.Pressure), Temperature(other.Temperature),platform_type(other.platform_type),
					Restraint_flag(other.Restraint_flag),CMMremoval_frequency(other.CMMremoval_frequency),
					buffer_frequency(other.buffer_frequency),device_index(other.device_index),
					boosted_torsions(other.boosted_torsions), etorsion(other.etorsion), alphatorsion(other.alphatorsion)
{}

/** Destructor */
OpenMMAMDIntegrator::~OpenMMAMDIntegrator()
{
  //delete openmm_system;
}

/** Copy assignment operator */
OpenMMAMDIntegrator& OpenMMAMDIntegrator::operator=(const OpenMMAMDIntegrator &other){
	Integrator::operator=(other);
	frequent_save_velocities = other.frequent_save_velocities;
	molgroup = other.molgroup; 
	openmm_system = other.openmm_system;
	isInitialised = other.isInitialised;
	Integrator_type = other.Integrator_type;
	friction = other.friction;
	CutoffType = other.CutoffType;
	cutoff_distance = other.cutoff_distance;
	field_dielectric = other.field_dielectric;
	Andersen_flag = other.Andersen_flag;
	Andersen_frequency = other.Andersen_frequency;
	MCBarostat_flag = other.MCBarostat_flag;
	MCBarostat_frequency = other.MCBarostat_frequency;
	ConstraintType = other.ConstraintType;
	Pressure = other.Pressure;
	Temperature = other.Temperature;
	platform_type = other.platform_type;
	Restraint_flag = other.Restraint_flag;
	CMMremoval_frequency = other.CMMremoval_frequency;
	buffer_frequency = other.buffer_frequency;
	device_index = other.device_index;
	boosted_torsions = other.boosted_torsions;
	etorsion = other.etorsion;
	alphatorsion = other.alphatorsion;

	return *this;
}

/** Comparison operator */
bool OpenMMAMDIntegrator::operator==(const OpenMMAMDIntegrator &other) const{
	return frequent_save_velocities == other.frequent_save_velocities and Integrator::operator==(other);
}

/** Comparison operator */
bool OpenMMAMDIntegrator::operator!=(const OpenMMAMDIntegrator &other) const
{
	return not OpenMMAMDIntegrator::operator==(other);
}

/** Return a string representation of this integrator */
QString OpenMMAMDIntegrator::toString() const
{
	return QObject::tr("OpenMMAMDIntegrator()");
}


void OpenMMAMDIntegrator::initialise(){

	bool Debug = true;

	if (Debug)
		qDebug() << " initialising OpenMMAMDIntegrator";

	// Create a workspace using the stored molgroup

	const MoleculeGroup moleculegroup = this->molgroup.read();

	if ( moleculegroup.isEmpty() ){
		throw SireError::program_bug(QObject::tr(
			"Cannot initialise OpenMMAMDIntegrator because molgroup has not been defined"), CODELOC);
	}

	AtomicVelocityWorkspace ws = this->createWorkspace(moleculegroup).read().asA<AtomicVelocityWorkspace>();

	const int nmols = ws.nMolecules();

	int nats = 0;

	for (int i=0; i<nmols; ++i){
		nats = nats + ws.nAtoms(i);
	}

	if (Debug)
		qDebug() << "There are " << nats << " atoms " << "There are " << nmols << " molecules" <<"\n" ;

	int flag_cutoff;
	int flag_constraint;

	if (CutoffType == "nocutoff")
		flag_cutoff = NOCUTOFF;
	else if (CutoffType == "cutoffnonperiodic")
		flag_cutoff = CUTOFFNONPERIODIC;
	else if (CutoffType == "cutoffperiodic")
		flag_cutoff = CUTOFFPERIODIC;
	else
		throw SireError::program_bug(QObject::tr("The CutOff method has not been specified. Possible choises: nocutoff, cutoffnonperiodic, cutoffperiodic"), CODELOC);

	if (Debug)
		qDebug() << "\nCutoffType = " << CutoffType << "\n";

	if (ConstraintType == "none")
		flag_constraint = NONE;
	else if (ConstraintType == "hbonds")
		flag_constraint = HBONDS;
	else if(ConstraintType == "allbonds")
		flag_constraint = ALLBONDS;
	else if (ConstraintType == "hangles")
		flag_constraint = HANGLES;
	else
		throw SireError::program_bug(QObject::tr("The Constraints method has not been specified. Possible choises: none, hbonds, allbonds, hangles"), CODELOC);

	if (Debug)
		qDebug() << "\nConstraint Type = " << ConstraintType<< "\n";

	//Load Plugins from the OpenMM standard Plugin Directory 
	OpenMM::Platform::loadPluginsFromDirectory(OpenMM::Platform::getDefaultPluginsDirectory());

	OpenMM::System *system_openmm = new OpenMM::System();

	OpenMM::NonbondedForce * nonbond_openmm = new OpenMM::NonbondedForce();

	nonbond_openmm->setUseDispersionCorrection(false);

	system_openmm->addForce(nonbond_openmm);

	if ( flag_cutoff == NOCUTOFF ){
		nonbond_openmm->setNonbondedMethod(OpenMM::NonbondedForce::NoCutoff);

		if (Debug)
			qDebug() << "\nCut off type = " << CutoffType << "\n";
	}
	else {
		const double converted_cutoff_distance = convertTo(cutoff_distance.value(), nanometer);

		if( flag_cutoff == CUTOFFNONPERIODIC )
			nonbond_openmm->setNonbondedMethod(OpenMM::NonbondedForce::CutoffNonPeriodic);
		else
			nonbond_openmm->setNonbondedMethod(OpenMM::NonbondedForce::CutoffPeriodic);

			nonbond_openmm->setCutoffDistance(converted_cutoff_distance);

		//Set Dielectric constant media
		nonbond_openmm->setReactionFieldDielectric(field_dielectric);

		if (Debug) {
			qDebug() << "\nCut off type = " << CutoffType << "\n";
			qDebug() << "CutOff distance = " << converted_cutoff_distance  << " Nm" << "\n";
			qDebug() << "Dielectric constant= " << field_dielectric << "\n\n";
		}
	}


	// Andersen thermostat. Complain if NOT using Verlet
	if (Andersen_flag == true){
		if (Integrator_type != "leapfrogverlet" or Integrator_type != "variableleapfrogverlet") {
			throw SireError::program_bug(QObject::tr("The Andersen thermostat can only be used with the leapfrogverlet or variableleapfrogverlet integrators"), CODELOC);
		}

	const double converted_Temperature = convertTo(Temperature.value(), kelvin);

	OpenMM::AndersenThermostat * thermostat = new OpenMM::AndersenThermostat(converted_Temperature, Andersen_frequency);

	system_openmm->addForce(thermostat);

		if (Debug) {
			qDebug() << "\nAndersen Thermostat set\n";
			qDebug() << "Temperature = " << converted_Temperature << " K\n";
			qDebug() << "Frequency collisions = " << Andersen_frequency << " 1/ps\n";
		}
	}

	// Monte Carlo Barostat
	if (MCBarostat_flag == true) {
		const double converted_Temperature = convertTo(Temperature.value(), kelvin);
		const double converted_Pressure = convertTo(Pressure.value(), bar);

		OpenMM::MonteCarloBarostat * barostat = new OpenMM::MonteCarloBarostat(converted_Pressure, converted_Temperature, MCBarostat_frequency);
		system_openmm->addForce(barostat);

		if (Debug) {
			qDebug() << "\nMonte Carlo Barostat set\n";
			qDebug() << "Temperature = " << converted_Temperature << " K\n";
			qDebug() << "Pressure = " << converted_Pressure << " bar\n";
			qDebug() << "Frequency every " << MCBarostat_frequency << " steps\n";
		}
	}


	//OpenMM Bonded Forces

	OpenMM::HarmonicBondForce * bondStretch_openmm = new OpenMM::HarmonicBondForce();

	OpenMM::HarmonicAngleForce * bondBend_openmm = new OpenMM::HarmonicAngleForce();

	OpenMM::PeriodicTorsionForce * bondTorsion_openmm = new OpenMM::PeriodicTorsionForce();

	system_openmm->addForce(bondStretch_openmm);

	system_openmm->addForce(bondBend_openmm);

	system_openmm->addForce(bondTorsion_openmm);

	OpenMM::CustomTorsionForce * TorsionBoost = NULL;

	TorsionBoost = new OpenMM::CustomTorsionForce("on_off * Hb *step(Eth - Ht);"
												  "Hb=((Eth - Ht)^2)/(Ath + Eth - Ht);"
												  "Ht=k*(1 + cos(n*theta - theta0))");

	double eth = etorsion.value();

	double ath = alphatorsion.value();

	TorsionBoost->addGlobalParameter("Eth", eth * OpenMM::KJPerKcal);

	TorsionBoost->addGlobalParameter("Ath", ath * OpenMM::KJPerKcal);

	TorsionBoost->addGlobalParameter("on_off", 1.0);

	system_openmm->addForce(TorsionBoost);


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

		system_openmm->addForce(positionalRestraints_openmm);

		if (Debug)
			qDebug() << "\nRestraint = ON\n\n";

	}

	//OpenMM vector coordinate
	std::vector<OpenMM::Vec3> positions_openmm(nats);

	//OpenMM vector momenta
	std::vector<OpenMM::Vec3> velocities_openmm(nats);

	std::vector<std::pair<int,int> > bondPairs;

	int system_index = 0;

	// To avoid possible mismatch between the index in which atoms are added to the openmm system arrays and 
	// their atomic numbers in sire, one array is populated while filling up the openmm global arrays
	//  AtomNumtoopenmmIndex
	QHash<int, int> AtomNumToOpenMMIndex;

	for (int i=0; i < nmols; ++i){
		const int nats_mol = ws.nAtoms(i);
		const double *m = ws.massArray(i);

		MolNum molnum = moleculegroup.molNumAt(i);
		const ViewsOfMol &molview = moleculegroup[molnum].data();
		const Molecule &mol = molview.molecule();
		Selector<Atom> molatoms = mol.atoms();

		for (int j=0; j < nats_mol; ++j){

			system_openmm->addParticle(m[j]) ;

			Atom at = molatoms.at(j);
			AtomNum atnum = at.number();

			if (Debug)
				qDebug() << " openMM_index " << system_index << " Sire Atom Number " << atnum.toString();

			AtomNumToOpenMMIndex[atnum.value()] = system_index;

			// JM Nov 12
			// The code below implements a ThreeParticleAverageSite for virtual sites for EPW atoms present in a WAT residue
			// This is very AMBER specific. 

			AtomName atname = at.name();

			if (Debug)
				qDebug() << " atname " << atname.value() << " mol " << i;

			if ( atname == AtomName("EPW") ) {
				ResName resname = at.residue().name();
				if ( resname == ResName("WAT") ){
					Atom oatom = molatoms.select( AtomName("O") );
					Atom h1atom = molatoms.select( AtomName("H1") );
					Atom h2atom = molatoms.select( AtomName("H2") );

					AmberParameters amber_params = mol.property("amberparameters").asA<AmberParameters>();
					QList<BondID> bonds_ff = amber_params.getAllBonds();

					double distoh = -1.0;
					double disthh = -1.0;
					double distoe = -1.0;

					for (int k=0; k < bonds_ff.length() ; k++) {
						BondID bond_ff = bonds_ff[k];
						QList<double> bond_params = amber_params.getParams(bond_ff);

						double r0 = bond_params[1];

					AtomName at0name = mol.select( bond_ff.atom0() ).name();
					AtomName at1name = mol.select( bond_ff.atom1() ).name();

					// qDebug() << " at0name " << at0name.toString() << " at1name " << at1name.toString();

						if ( ( at0name == AtomName("O") and at1name == AtomName("H1") ) or ( at0name == AtomName("H1") and at1name == AtomName("O") ) ) {
							distoh = r0;
						}
						else if ( ( at0name == AtomName("H1") and at1name == AtomName("H2") ) or ( at0name == AtomName("H2") and at1name == AtomName("H1") ) ) {
							disthh = r0;
						}
						else if ( ( at0name == AtomName("EPW") and at1name== AtomName("O") ) or ( at0name == AtomName("O") and at1name == AtomName("EPW") ) ) {
							distoe = r0;
						}
					}

					if ( distoh < 0 or disthh < 0 or distoe < 0){
						throw SireError::program_bug(QObject::tr("Could not find expected atoms in TIP4P water molecule."), CODELOC);
					}

					//qDebug() << " distoe " << distoe << " distoh " << distoh << " disthh " << disthh;

					double weightH = distoe / sqrt( (distoh*distoh) - ( 0.25 * disthh*disthh) );

					int o_index = AtomNumToOpenMMIndex[oatom.number().value()];
					int h1_index = AtomNumToOpenMMIndex[h1atom.number().value()];
					int h2_index = AtomNumToOpenMMIndex[h2atom.number().value()];

					if (Debug)
						qDebug() << "virtual site " << system_index << " o " << o_index << " h1 " << h1_index << " h2 " << h2_index << " 1 - weightH " << 1 - weightH << " weightH/2 " << weightH/2 ;

					OpenMM::ThreeParticleAverageSite * vsite =  new OpenMM::ThreeParticleAverageSite(o_index, h1_index, h2_index, 1-weightH, weightH/2, weightH/2);

					system_openmm->setVirtualSite( system_index, vsite );

				}
			}

			system_index = system_index + 1;

		}// end of loop on atoms in molecule
	}//end of loop on molecules in workspace

	int num_atoms_till_i = 0;

	for (int i=0; i < nmols ; i++){
		const Vector *c = ws.coordsArray(i);

		Molecule molecule = moleculegroup.moleculeAt(i).molecule();
		int num_atoms_molecule = molecule.nAtoms();

		// The atomic parameters
		AtomLJs atomvdws = molecule.property("LJ").asA<AtomLJs>();
		AtomCharges atomcharges = molecule.property("charge").asA<AtomCharges>();
		QVector<SireMM::LJParameter> ljparameters = atomvdws.toVector();
		QVector<SireUnits::Dimension::Charge> charges = atomcharges.toVector();

		for(int j=0; j< ljparameters.size();j++){
			double sigma = ljparameters[j].sigma();
			double epsilon = ljparameters[j].epsilon();
			double charge = charges[j].value();

			nonbond_openmm->addParticle(charge, sigma * OpenMM::NmPerAngstrom, epsilon * OpenMM::KJPerKcal);
		}


		if(Restraint_flag == true){
			bool hasRestrainedAtoms = molecule.hasProperty("restrainedatoms");

			if(hasRestrainedAtoms){
				Properties restrainedAtoms = molecule.property("restrainedatoms").asA<Properties>();

				int nrestrainedatoms = restrainedAtoms.property(QString("nrestrainedatoms")).asA<VariantProperty>().toInt();
				if (Debug)
					qDebug() << " nrestrainedatoms " << nrestrainedatoms ;

				for (int i=0; i < nrestrainedatoms ; i++){
					int atomnum = restrainedAtoms.property(QString("AtomNum(%1)").arg(i)).asA<VariantProperty>().toInt();
					double xref = restrainedAtoms.property(QString("x(%1)").arg(i)).asA<VariantProperty>().toDouble();
					double yref = restrainedAtoms.property(QString("y(%1)").arg(i)).asA<VariantProperty>().toDouble();
					double zref = restrainedAtoms.property(QString("z(%1)").arg(i)).asA<VariantProperty>().toDouble();
					double k = restrainedAtoms.property(QString("k(%1)").arg(i)).asA<VariantProperty>().toDouble();

					int openmmindex = AtomNumToOpenMMIndex[atomnum];

					if (Debug)
						qDebug() << " atomnum " << atomnum << " openmmindex " << openmmindex << " x " << xref << " y " << yref << " z " << zref << " k " << k;

					int posrestrdim = 4;
					std::vector<double> params(posrestrdim);

					params[0] = xref * OpenMM::NmPerAngstrom;
					 params[1] = yref * OpenMM::NmPerAngstrom;
					params[2] = zref * OpenMM::NmPerAngstrom;
					params[3] = k  * ( OpenMM::KJPerKcal * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm );

					positionalRestraints_openmm->addParticle(openmmindex, params);
				}
			}
		}//end of restraint flag


		// The bonded parameters
		bool hasConnectivity = molecule.hasProperty("connectivity");

		if ( !hasConnectivity ){
			num_atoms_till_i = num_atoms_till_i + num_atoms_molecule ;
			if (Debug) {
				qDebug() << "\nAtoms = " <<  num_atoms_molecule << " Num atoms till i =" << num_atoms_till_i <<"\n";
				qDebug() << "\n*********************MONOATOMIC MOLECULE DETECTED**************************\n";
			}
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
				bondStretch_openmm->addBond(idx0,idx1,r0 * OpenMM::NmPerAngstrom, k * 2.0 * OpenMM::KJPerKcal * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm); 

			//cout << "\nBOND ADDED TO "<< atom0.toStdString() << " AND " << atom1.toStdString() << "\n";
			}
			else if ( flag_constraint == ALLBONDS || flag_constraint == HANGLES ) {
				system_openmm->addConstraint(idx0,idx1, r0 *  OpenMM::NmPerAngstrom);
				//cout << "\nALLBONDS or HANGLES ADDED BOND CONSTRAINT TO " << atom0.toStdString() << " AND " << atom1.toStdString() << "\n";
			}
			else if ( flag_constraint == HBONDS ){
				if ( (atom0[6] == 'H') || (atom1[6] == 'H') ){
					system_openmm->addConstraint(idx0,idx1, r0 *  OpenMM::NmPerAngstrom);
					//cout << "\nHBONDS ADDED BOND CONSTRAINT TO " << atom0.toStdString() << " AND " << atom1.toStdString() << "\n";
				}
				else {
					bondStretch_openmm->addBond(idx0,idx1,r0 * OpenMM::NmPerAngstrom,k * 2.0 * OpenMM::KJPerKcal * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm);
					//cout << "\nHBONDS ADDED BOND TO " << atom0.toStdString() << " AND " << atom1.toStdString() << "\n";
				}
			}
			//Bond exclusion List
			bondPairs.push_back(std::make_pair(idx0,idx1));
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

			idx0 = idx0 + num_atoms_till_i;
			idx1 = idx1 + num_atoms_till_i;
			idx2 = idx2 + num_atoms_till_i;

			if ( flag_constraint == HANGLES ){
				if( ((atom0[6] == 'H') && (atom2[6] == 'H')) ){
					system_openmm->addConstraint(idx0,idx2, diff.length() *  OpenMM::NmPerAngstrom);
				}
				else if ( ((atom0[6] == 'H') && (atom1[6] == 'O')) || ((atom1[6] == 'O') && (atom2[6] == 'H')) ) {
					system_openmm->addConstraint(idx0,idx2, diff.length() *  OpenMM::NmPerAngstrom);
				}
				else{
					bondBend_openmm->addAngle(idx0,idx1,idx2, theta0 , k * 2.0 * OpenMM::KJPerKcal);
				}
			}
			else {
				bondBend_openmm->addAngle(idx0,idx1,idx2, theta0 , k * 2.0 * OpenMM::KJPerKcal);
			}
		}//end of angles

		//Dihedrals
		QList<DihedralID> dihedrals_ff = amber_params.getAllDihedrals();
		QVector<DihedralID> dihedrals = dihedrals_ff.toVector();


		TorsionBoost->addPerTorsionParameter("theta0");
		TorsionBoost->addPerTorsionParameter("n");
		TorsionBoost->addPerTorsionParameter("k");


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

				if ( boosted_torsions.contains(dihedral_ff) ){

					vector<double> parameters(3);
					parameters[0] = phase;
					parameters[1] = periodicity;
					parameters[2] = v * OpenMM::KJPerKcal;

					TorsionBoost->addTorsion(idx0, idx1, idx2, idx3, parameters);

					qDebug() << " The selected Dihedral is in the Boost List ";

				}

				if (Debug){
					qDebug() << "Dihedral between atom global index " << idx0 << " and " << idx1 << " and " << idx2 << " and " << idx3;
					qDebug() << "Amplitude_dih = " << v << " periodicity " << periodicity << " phase " << phase;
					qDebug() << "Dihedral local" << dihedral_ff.toString() << " v " << v << " periodicity " << periodicity << " phase " << phase;
				}
			}

		} // end of dihedrals


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
		}//end of impropers

		num_atoms_till_i = num_atoms_till_i + num_atoms_molecule ;
	}// end of loop over molecules

	//Exclude the 1-2, 1-3 bonded atoms from nonbonded forces, and scale down 1-4 bonded atoms
	const double Coulomb14Scale = 1.0/1.2;
	const double LennardJones14Scale = 1.0/2.0;
	nonbond_openmm->createExceptionsFromBonds(bondPairs, Coulomb14Scale, LennardJones14Scale);

	if(CMMremoval_frequency > 0){
		OpenMM::CMMotionRemover * cmmotionremover = new OpenMM::CMMotionRemover(CMMremoval_frequency);

		system_openmm->addForce(cmmotionremover);

		if (Debug)
			qDebug() << "\n\nWill remove Center of Mass motion every " << CMMremoval_frequency << " steps\n\n";
	}

	this->openmm_system = system_openmm;

	this->isInitialised = true;
}

void OpenMMAMDIntegrator::integrate(IntegratorWorkspace &workspace, const Symbol &nrg_component, SireUnits::Dimension::Time timestep, int nmoves, bool record_stats)  {

	bool Debug = false; 

	QTime timer;

	timer.start();

	if (Debug)
		qDebug() << "In OpenMMAMDIntegrator::integrate()\n\n" ;

	// Check that the openmm system has been initialised
	// !! Should check that the workspace is compatible with molgroup
	if ( not this->isInitialised){
			qDebug() << "Not initialised ! ";
			throw SireError::program_bug(QObject::tr("OpenMMMDintegrator should have been initialised before calling integrate."), CODELOC);
	}

	OpenMM::System *system_openmm = openmm_system;

	int nats = system_openmm->getNumParticles();

	if (Debug)
		qDebug() << " openmm nats " << nats;

	// Create here integrator + Platform, then context
	const double dt = convertTo( timestep.value(), picosecond);
	const double converted_Temperature = convertTo(Temperature.value(), kelvin);
	const double converted_friction = convertTo( friction.value(), picosecond);
	// Currently not exposed to Sire API
	const double integration_tol = 0.001;

	OpenMM::Integrator * integrator_openmm = NULL;

	if (Integrator_type == "leapfrogverlet")
		integrator_openmm = new OpenMM::VerletIntegrator(dt);//dt in picosecond
		//OpenMM::VerletIntegrator integrator_openmm(dt);
	else if (Integrator_type == "variableleapfrogverlet")
		integrator_openmm = new OpenMM::VariableVerletIntegrator(integration_tol);//integration tolerance error unitless
		//OpenMM::VariableVerletIntegrator integrator_openmm(integration_tol);// integrator tolerance error unitless
	else if (Integrator_type == "langevin")
		integrator_openmm = new OpenMM::LangevinIntegrator(converted_Temperature, converted_friction, dt);
		//OpenMM::LangevinIntegrator integrator_openmm(converted_Temperature, converted_friction, dt);
	else if (Integrator_type == "variablelangevin")
		integrator_openmm = new OpenMM::VariableLangevinIntegrator(converted_Temperature, converted_friction, integration_tol);
		//OpenMM::VariableLangevinIntegrator integrator_openmm(converted_Temperature, converted_friction, integration_tol);
	else if (Integrator_type == "brownian")
		integrator_openmm = new OpenMM::BrownianIntegrator(converted_Temperature, converted_friction, dt);
		//OpenMM::BrownianIntegrator(converted_Temperature, converted_friction, dt);
	else 
		throw SireError::program_bug(QObject::tr("The user defined Integrator type is not supported. Available types are leapfrogverlet, variableleapfrogverlet, langevin, variablelangevin, brownian"), CODELOC);

	if (Debug)
		qDebug() << "Using Integrator; " << Integrator_type;

	OpenMM::Platform& platform_openmm = OpenMM::Platform::getPlatformByName(platform_type.toStdString()); 	

	if (platform_type == "OpenCL"){
		const std::string prop = std::string("OpenCLDeviceIndex");
		platform_openmm.setPropertyDefaultValue(prop, device_index.toStdString() );
		//qDebug() << " setting up OpenCL default Index to " << device_index;
	}
	else if (platform_type == "Cuda"){
		const std::string prop = std::string("CudaDeviceIndex");
		platform_openmm.setPropertyDefaultValue(prop, device_index.toStdString() );   
	}
	// JM Dec 12. Do we need to set for Reference?

	// Creating a context is the bottleneck in the setup step 
	// Another implementation could have the context created once during initialisation. 
	// But then have to figure out how to properly allocate/free context on the heap and make it 
	// compatible with sire objects
	OpenMM::Context context_openmm( *system_openmm, *integrator_openmm, platform_openmm);  

	if (Debug)
		qDebug() << "\n Using OpenMM platform = " <<context_openmm.getPlatform().getName().c_str()<<"\n";

	// Now update coordinates / velocities / dimensions with sire data

	//OpenMM vector coordinate
	std::vector<OpenMM::Vec3> positions_openmm(nats);
	//OpenMM vector momenta
	std::vector<OpenMM::Vec3> velocities_openmm(nats);

	AtomicVelocityWorkspace &ws = workspace.asA<AtomicVelocityWorkspace>();

	// Conversion factor because sire units of time are in AKMA, whereas OpenMM uses picoseconds
	double AKMAPerPs = 0.04888821;
	double PsPerAKMA = 1 / AKMAPerPs;

	const int nmols = ws.nMolecules();

	int system_index = 0;

	for (int i=0; i < nmols; ++i){
		const int nats_mol = ws.nAtoms(i);
		Vector *c = ws.coordsArray(i);
		Vector *p = ws.momentaArray(i);
		const double *m = ws.massArray(i);

		for (int j=0; j < nats_mol; ++j){
			positions_openmm[system_index] = OpenMM::Vec3(c[j].x() * (OpenMM::NmPerAngstrom),c[j].y() * (OpenMM::NmPerAngstrom),c[j].z() * (OpenMM::NmPerAngstrom));
			if (m[j] > SireMaths::small) {
				velocities_openmm[system_index] = OpenMM::Vec3(p[j].x()/m[j] * (OpenMM::NmPerAngstrom) * PsPerAKMA, p[j].y()/m[j] * (OpenMM::NmPerAngstrom) * PsPerAKMA,p[j].z()/m[j] * (OpenMM::NmPerAngstrom) * PsPerAKMA);
			}
			else{
				velocities_openmm[system_index] =  OpenMM::Vec3(0.0, 0.0, 0.0);
			}
			system_index++;
		}
	}

	if ( system_index != nats ){
		if (Debug)
			qDebug() << " system_index " << system_index << " nats " << nats;
			throw SireError::program_bug(QObject::tr("The number of atoms in the openmm system does not match the number of atoms in the sire workspace"), CODELOC);
	}

	context_openmm.setPositions(positions_openmm);  
	context_openmm.setVelocities(velocities_openmm);

	bool isperiodic = false;

	if ( CutoffType == "cutoffperiodic" ){
		const System & ptr_sys = ws.system();
		const PropertyName &space_property = PropertyName("space");
		const PeriodicBox &space = ptr_sys.property(space_property).asA<PeriodicBox>();

		const double Box_x_Edge_Length = space.dimensions()[0] * OpenMM::NmPerAngstrom; //units in nm
		const double Box_y_Edge_Length = space.dimensions()[1] * OpenMM::NmPerAngstrom; //units in nm
		const double Box_z_Edge_Length = space.dimensions()[2] * OpenMM::NmPerAngstrom; //units in nm

		if (Debug)
			qDebug() << "\nBOX SIZE [A] = (" << space.dimensions()[0] << " , " << space.dimensions()[1] << " ,  " << space.dimensions()[2] << ")\n\n";

		//Set Periodic Box Condition

		context_openmm.setPeriodicBoxVectors( OpenMM::Vec3(Box_x_Edge_Length,0,0), OpenMM::Vec3(0,Box_y_Edge_Length,0),OpenMM::Vec3(0,0,Box_z_Edge_Length) );
		isperiodic = true;
	}

	int infoMask = 0;
	infoMask = OpenMM::State::Positions;
	infoMask = infoMask + OpenMM::State::Velocities; 
	infoMask = infoMask +  OpenMM::State::Energy;

	if (Debug)
		qDebug() << " Setup dynamics, time elapsed ms " << timer.elapsed() << " ms ";

	/** Now perform some steps of dynamics */
	timer.restart();

	if (Debug)
		qDebug() << " Doing " << nmoves << " steps of dynamics ";

	// Coordinates are buffered every coord_freq steps
	int coord_freq = buffer_frequency;

	int nframes;
	int MAXFRAMES = 1000;

	// Limit excessive internal buffering
	if ( coord_freq > 0 ){
		nframes = ( nmoves / coord_freq ) ;

		if  ( nframes > MAXFRAMES ) {
			throw SireError::program_bug(QObject::tr("You are requesting to buffer %1 frames, which is above the hardcoded limit of %2.").arg(nframes, MAXFRAMES), CODELOC);
		}
	}
	else{
		nframes = 0;
	}

	QVector< std::vector<OpenMM::Vec3> > buffered_positions;
	QVector< Vector> buffered_dimensions;

	OpenMM::State state_openmm;

	OpenMM::Vec3 a;
	OpenMM::Vec3 b;
	OpenMM::Vec3 c;


	OpenMM::Integrator& context_integrator = context_openmm.getIntegrator();

	if ( coord_freq > 0 ){/** Break nmoves in several steps to buffer coordinates*/
		for (int i=0; i < nmoves ; i = i + coord_freq){
			context_integrator.step(coord_freq);
			if (Debug)
				qDebug() << " i now " << i;

			state_openmm = context_openmm.getState(infoMask);
			positions_openmm = state_openmm.getPositions();
			buffered_positions.append( positions_openmm );

			state_openmm.getPeriodicBoxVectors(a,b,c);
			Vector dims = Vector( a[0] * OpenMM::AngstromsPerNm, b[1] * OpenMM::AngstromsPerNm, c[2] * OpenMM::AngstromsPerNm);
			buffered_dimensions.append( dims );
			}
	}
	else{/** No buffering*/
		context_integrator.step(nmoves);
	}

	if (Debug)
		qDebug() << " Done dynamics, time elapsed ms " << timer.elapsed() << " ms ";

	/** Now update the sire coordinates/velocities and box dimensions */
	timer.restart();

	if (Debug) {
		double kinetic_energy = state_openmm.getKineticEnergy(); 
		double potential_energy = state_openmm.getPotentialEnergy(); 
		qDebug() << " After MD kinetic energy " << kinetic_energy  * OpenMM::KcalPerKJ  << " kcal/mol potential " << potential_energy  * OpenMM::KcalPerKJ << " kcal/mol ";
	}

	state_openmm = context_openmm.getState(infoMask);	
	positions_openmm = state_openmm.getPositions();
	velocities_openmm = state_openmm.getVelocities();

	// Vector of Vector of molecules that are vector of atomic coordinates...
	QVector< QVector< QVector< Vector > > > buffered_workspace(nframes);
	for (int i=0; i < buffered_workspace.size() ; i++){
		buffered_workspace[i].resize(nmols);
		for (int j=0 ; j < nmols; j++){
			int nats = ws.nAtoms(j);
			buffered_workspace[i][j].resize(nats);
		}
	}

	int k=0;

	for(int i=0; i<nmols;i++){
		Vector *sire_coords = ws.coordsArray(i);
		Vector *sire_momenta = ws.momentaArray(i);
		const double *m = ws.massArray(i);

		for(int j=0; j < ws.nAtoms(i) ; j++){

			sire_coords[j] = Vector(positions_openmm[j+k][0] * (OpenMM::AngstromsPerNm), positions_openmm[j+k][1] * (OpenMM::AngstromsPerNm), positions_openmm[j+k][2] * (OpenMM::AngstromsPerNm));

			for (int l=0; l < nframes ; l++){
				//qDebug() << " i " << i << " j " << j << " k " << k << " l " << l;

				Vector buffered_atcoord = Vector(  buffered_positions[l][j+k][0] * (OpenMM::AngstromsPerNm),buffered_positions[l][j+k][1] * (OpenMM::AngstromsPerNm),buffered_positions[l][j+k][2] * (OpenMM::AngstromsPerNm) );
				buffered_workspace[l][i][j] = buffered_atcoord;

			}

			sire_momenta[j] = Vector(velocities_openmm[j+k][0] * m[j] * (OpenMM::AngstromsPerNm) * AKMAPerPs,velocities_openmm[j+k][1] * m[j] * (OpenMM::AngstromsPerNm) * AKMAPerPs,velocities_openmm[j+k][2] * m[j] * (OpenMM::AngstromsPerNm) * AKMAPerPs);
		}

		k= k + ws.nAtoms(i);
	}

	if ( nframes <= 0 )
		ws.commitCoordinatesAndVelocities();
	else
		ws.commitBufferedCoordinatesAndVelocities( buffered_workspace );

	/** Now the box dimensions (if the simulation used a periodic space) */
	if (isperiodic){
		state_openmm.getPeriodicBoxVectors(a,b,c);

		Vector new_dims = Vector(a[0] * OpenMM::AngstromsPerNm, b[1] * OpenMM::AngstromsPerNm, c[2] * OpenMM::AngstromsPerNm);

		if (Debug)
			qDebug() << " a " << a[0] << " b " << b[1] << " c " << c[2];

		System & ptr_sys = ws.nonConstsystem();
		PeriodicBox  sp = ptr_sys.property("space").asA<PeriodicBox>();

		sp.setDimensions(new_dims);
		const QString string = "space" ;
		ptr_sys.setProperty(string, sp);

		/** Buffer dimensions if necessary */
		for (int k=0; k < buffered_dimensions.size() ; k++){
			const QString buffered_space = "buffered_space_" + QString::number(k) ;
			PeriodicBox buff_space = PeriodicBox( buffered_dimensions[k] );
			ptr_sys.setProperty( buffered_space, buff_space);
		}
	}

	/** Clear all buffers */
	buffered_workspace.clear();
	buffered_dimensions.clear();

	if (Debug)
		qDebug() << " Updating system coordinates, time elapsed ms " << timer.elapsed() << " ms ";

	return;
}


QString OpenMMAMDIntegrator::getIntegrator(void) {
  return Integrator_type;
}

void OpenMMAMDIntegrator::setIntegrator(QString intgrator) {
  Integrator_type = intgrator;
}

SireUnits::Dimension::Time OpenMMAMDIntegrator::getFriction(void) {
  return friction;
}

void OpenMMAMDIntegrator::setFriction(SireUnits::Dimension::Time thefriction) {
  friction = thefriction;
} 

/** Get the cufott type: nocutoff, cutoffnonperiodic, cutoffperiodic */
QString OpenMMAMDIntegrator::getCutoffType(void){

	return CutoffType;

}

/** Set the cufott type: nocutoff, cutoffnonperiodic, cutoffperiodic */
void OpenMMAMDIntegrator::setCutoffType(QString cutoff_type){

	CutoffType = cutoff_type;

}


/** Get the cutoff distance in A */
SireUnits::Dimension::Length OpenMMAMDIntegrator::getCutoff_distance(void){

	return cutoff_distance;

}

/** Set the cutoff distance in A */
void OpenMMAMDIntegrator::setCutoff_distance(SireUnits::Dimension::Length distance){

	cutoff_distance = distance;	

}

/** Get the dielectric constant */
double OpenMMAMDIntegrator::getField_dielectric(void){

	return field_dielectric;
}

/** Set the dielectric constant */
void OpenMMAMDIntegrator::setField_dielectric(double dielectric){
	
	field_dielectric=dielectric;

}

/** Set Andersen thermostat */

void OpenMMAMDIntegrator::setAndersen(bool andersen){
	Andersen_flag = andersen;	
}

/** Get Andersen thermostat status on/off */
bool OpenMMAMDIntegrator::getAndersen(void){
	
	return 	Andersen_flag;
	
}


/** Get the Andersen Thermostat frequency collision */
double OpenMMAMDIntegrator::getAndersen_frequency(void){
	
	return Andersen_frequency;
	
}

/** Set the Andersen Thermostat frequency collision */
void OpenMMAMDIntegrator::setAndersen_frequency(double freq){
	
	Andersen_frequency=freq;
	
}


/** Get the bath Temperature */
SireUnits::Dimension::Temperature OpenMMAMDIntegrator::getTemperature(void){
		
	return Temperature;
}

/** Set the Temperature */
void OpenMMAMDIntegrator::setTemperature(SireUnits::Dimension::Temperature temperature){
		
	Temperature = temperature;
}

/** Set Monte Carlo Barostat on/off */
void OpenMMAMDIntegrator::setMCBarostat(bool MCBarostat){
	MCBarostat_flag = MCBarostat;
}

/** Get Andersen thermostat status on/off */
bool OpenMMAMDIntegrator::getMCBarostat(void){
	
	return 	MCBarostat_flag;
	
}

/** Get the Monte Carlo Barostat frequency in time speps */
int OpenMMAMDIntegrator::getMCBarostat_frequency(void){
	
	return  MCBarostat_frequency;
	
}

/** Set the Monte Carlo Barostat frequency in time speps */
void OpenMMAMDIntegrator::setMCBarostat_frequency(int freq){
	
	MCBarostat_frequency=freq;
	
}

/** Get the Presure */
SireUnits::Dimension::Pressure OpenMMAMDIntegrator::getPressure(void){
		
	return Pressure;
}

/** Set the Pressure */
void OpenMMAMDIntegrator::setPressure(SireUnits::Dimension::Pressure pressure){

	Pressure = pressure;
}


/** Get the Constraint type: none, hbonds, allbonds, hangles */
QString OpenMMAMDIntegrator::getConstraintType(void){

	return ConstraintType;

}

/** Set the Constraint type: none, hbonds, allbonds, hangles */
void OpenMMAMDIntegrator::setConstraintType(QString constrain){

	ConstraintType = constrain;

}

/** Get the OpenMMMD Platform: CUDA, OpenCL, CPU */
QString OpenMMAMDIntegrator::getPlatform(void){

	return platform_type;

}

/** Set the OpenMM Platform: CUDA, OpenCL, CPU */
void OpenMMAMDIntegrator::setPlatform(QString platform){

	platform_type = platform;

}

/** Get the OpenMMMD Platform: CUDA, OpenCL, CPU */
QString OpenMMAMDIntegrator::getDeviceIndex(void){

	return device_index;

}

/** Set the OpenMM Platform: CUDA, OpenCL, CPU */
void OpenMMAMDIntegrator::setDeviceIndex(QString deviceidx){

	device_index = deviceidx;

}

/** Get the Restaint mode*/
bool OpenMMAMDIntegrator::getRestraint(void){

	return Restraint_flag;

}

/** Set the Retraint mode */
void OpenMMAMDIntegrator::setRestraint(bool Restraint){

	Restraint_flag = Restraint;
}

/** Get the Center of Mass motion removal frequency */
int OpenMMAMDIntegrator::getCMMremoval_frequency(void){

	return CMMremoval_frequency;
}

/** Set the Center of Mass motion removal frequency */
void OpenMMAMDIntegrator::setCMMremoval_frequency(int frequency){

	CMMremoval_frequency = frequency;
}

/** Get the frequency of buffering coordinates */
int OpenMMAMDIntegrator::getBufferFrequency(){

	return buffer_frequency;
}

/** Set the Center of Mass motion removal frequency */
void OpenMMAMDIntegrator::setBufferFrequency(int frequency){

	buffer_frequency = frequency;
}

/** Set the torsions to boost **/
void OpenMMAMDIntegrator::setBoostedTorsions( QList<DihedralID> torsions )
{
  boosted_torsions = torsions;
}

/** Get the torsions that are being boosted **/
QList<DihedralID> OpenMMAMDIntegrator::getBoostedTorsions(void)
{
  return boosted_torsions;
}

/** Set the level of boosting (Energy)**/
void OpenMMAMDIntegrator::setEtorsion( SireUnits::Dimension::MolarEnergy energy)
{
  etorsion = energy;
}

/** Get the level of boosting (energy) **/
SireUnits::Dimension::MolarEnergy OpenMMAMDIntegrator::getEtorsion(void)
{
  return etorsion;
}

/** Set the boost acceleration (alpha) **/
void OpenMMAMDIntegrator::setAlphatorsion(SireUnits::Dimension::MolarEnergy alpha)
{
  alphatorsion = alpha;
}

/** Get the boost acceleration (alpha) **/
SireUnits::Dimension::MolarEnergy OpenMMAMDIntegrator::getAlphatorsion(void)
{
  return alphatorsion;
}


/** Create an empty workspace */
IntegratorWorkspacePtr OpenMMAMDIntegrator::createWorkspace(const PropertyMap &map) const
{
	return IntegratorWorkspacePtr( new AtomicVelocityWorkspace(map) );
}

/** Return the ensemble of this integrator */
Ensemble OpenMMAMDIntegrator::ensemble() const
{
	return Ensemble::NVE();
}

/** Return whether or not this integrator is time-reversible */
bool OpenMMAMDIntegrator::isTimeReversible() const
{
	return true;
}

/** Create a workspace for this integrator for the molecule group 'molgroup' */
IntegratorWorkspacePtr OpenMMAMDIntegrator::createWorkspace(const MoleculeGroup &molgroup, const PropertyMap &map) const
{
	return IntegratorWorkspacePtr( new AtomicVelocityWorkspace(molgroup,map) );
}

const char* OpenMMAMDIntegrator::typeName()
{
	return QMetaType::typeName( qMetaTypeId<OpenMMAMDIntegrator>() );
}













