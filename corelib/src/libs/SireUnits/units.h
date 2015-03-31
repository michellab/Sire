/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2006  Christopher Woods
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

#ifndef SIREUNITS_UNITS_H
#define SIREUNITS_UNITS_H

//skip this completely when parsing with gccxml as it is broken!
#ifndef SKIP_BROKEN_GCCXML_PARTS

#include <limits>
#include <cmath>

#include "dimensions.h"

SIRE_BEGIN_HEADER

namespace SireUnits
{

/** This file defines physical constants, in internal units of this program

    We use the AKMA units (same as charmm)
        Angstroms, Kilocalories per Mole, Atomic mass units

    energy = kcal mol-1 (thermal)     (really MolarEnergy)
    length = angstrom
    mass = g mol-1                    (really MolarMass)
    time = AKMA time == 48.88821 fs == 0.04888821 ps
    charge = unit electrons

    Where necessary, physical constants were downloaded from
    the NIST website (web pages referenced where appropriate).
    
    Physical constants were last checked on 28/10/2008
*/

/** Avogadro's number */
//http://physics.nist.gov/cgi-bin/cuu/Value?na|search_for=physchem_in!
const Dimension::Quantity mole( 6.02214179e23 );

const Dimension::Quantity dozen( 12 );

/////////////////////////////////////////////////
// Units of angle. Internal unit = radians     //
/////////////////////////////////////////////////

const Dimension::Angle radians( 1 );
const Dimension::Angle radian( 1 );

const Dimension::Angle degrees = radians * pi / 180.0;
const Dimension::Angle degree = degrees;

const Dimension::Angle angle_minute = degree / 60;
const Dimension::Angle angle_minutes = angle_minute;

const Dimension::Angle angle_second = angle_minute / 60;
const Dimension::Angle angle_seconds = angle_second;

const Dimension::Angle octant = 45 * degrees;
const Dimension::Angle octants = octant;

const Dimension::Angle sextant = 60 * degrees;
const Dimension::Angle sextants = sextant;

const Dimension::Angle quadrant = 90 * degrees;
const Dimension::Angle quadrants = quadrant;

const Dimension::Angle gradian = quadrant / 100;
const Dimension::Angle gradians = gradian;
const Dimension::Angle grad = gradian;
const Dimension::Angle gon = gradian;

const Dimension::Angle revolution = 360 * degrees;
const Dimension::Angle revolutions = revolution;
const Dimension::Angle revs = revolution;

const Dimension::Angle circumference = revolution;

/////////////////////////////////////////////////
// Units of length. Internal unit = Angstroms  //
/////////////////////////////////////////////////

const Dimension::Length angstrom(1);
const Dimension::Length angstroms = angstrom;

const Dimension::Length picometer( 0.01 * angstrom );
const Dimension::Length nanometer( 1000 * picometer );
const Dimension::Length micrometer( 1000 * nanometer );
const Dimension::Length millimeter( 1000 * micrometer );
const Dimension::Length centimeter( 10 * millimeter );
const Dimension::Length meter( 100 * centimeter );
const Dimension::Length kilometer( 1000 * meter );

const Dimension::Length picometers = picometer;
const Dimension::Length nanometers = nanometer;
const Dimension::Length micrometers = micrometer;
const Dimension::Length millimeters = millimeter;
const Dimension::Length centimeters = centimeter;
const Dimension::Length meters = meter;
const Dimension::Length kilometers = kilometer;

//http://physics.nist.gov/cgi-bin/cuu/Value?bohrrada0|search_for=bohr
const Dimension::Length bohr_radii( 0.52917720859 * angstrom );

const Dimension::Length inch( 2.54 * centimeter );
const Dimension::Length foot( 12 * inch );
const Dimension::Length yard( 3 * foot );
const Dimension::Length mile( 1760 * yard );

const Dimension::Length inches = inch;
const Dimension::Length feet = foot;
const Dimension::Length yards = yard;
const Dimension::Length miles = mile;

/////////////////////////////////////////////////
// Units of area. Internal unit = Angstroms^2  //
/////////////////////////////////////////////////

const Dimension::Area angstrom2 = angstrom * angstrom;
const Dimension::Area picometer2 = picometer * picometer;
const Dimension::Area nanometer2 = nanometer * nanometer;

const Dimension::Area meter2 = meter * meter;

/////////////////////////////////////////////////
// Units of volume. Internal unit = Angstroms^3  //
/////////////////////////////////////////////////

const Dimension::Volume angstrom3 = angstrom * angstrom * angstrom;
const Dimension::Volume picometer3 = picometer * picometer * picometer;
const Dimension::Volume nanometer3 = nanometer * nanometer * nanometer;

const Dimension::Volume meter3 = meter * meter * meter;

///////////////////////////////////////////////////////
// Units of mass. Internal unit = g mol-1            //
///////////////////////////////////////////////////////

const Dimension::MolarMass g_per_mol(1);

const Dimension::Mass gram( mole * g_per_mol );
const Dimension::Mass kilogram( 1000 * gram );
const Dimension::Mass tonne( 1000 * kilogram );

const Dimension::Mass milligram( 0.001 * gram );
const Dimension::Mass microgram( 0.001 * milligram );
const Dimension::Mass nanogram( 0.001 * microgram );
const Dimension::Mass picogram( 0.001 * nanogram );
const Dimension::Mass femtogram( 0.001 * picogram );

const Dimension::MolarMass kg_per_mol( 1000 * g_per_mol );
const Dimension::MolarMass tonne_per_mol( 1000 * kg_per_mol );
const Dimension::MolarMass mg_per_mol( 0.001 * g_per_mol );
const Dimension::MolarMass ug_per_mol( 0.001 * mg_per_mol );
const Dimension::MolarMass ng_per_mol( 0.001 * ug_per_mol );
const Dimension::MolarMass pg_per_mol( 0.001 * ng_per_mol );
const Dimension::MolarMass fg_per_mol( 0.001 * pg_per_mol );

///////////////////////////////////////////////////////
// Units of Charge. Internal unit = |e|              //
///////////////////////////////////////////////////////

const Dimension::Charge mod_electron(1);
const Dimension::MolarCharge faraday(1);

//http://physics.nist.gov/cgi-bin/cuu/Value?e|search_for=elecmag_in!
const Dimension::Charge coulomb = mod_electron / 1.602176487e-19;
const Dimension::MolarCharge coulomb_per_mol = coulomb / mole;

const Dimension::Charge e_charge = -mod_electron;

/////////////////////////////////////////////////
// Units of Energy. Internal unit = kcal mol-1 //
/////////////////////////////////////////////////

const Dimension::MolarEnergy kcal_per_mol(1);
const Dimension::Energy kcal = mole * kcal_per_mol;

const Dimension::MolarEnergy cal_per_mol = 0.001 * kcal_per_mol;
const Dimension::Energy cal = 0.001 * kcal;

const Dimension::MolarEnergy kJ_per_mol = kcal_per_mol / 4.184;
const Dimension::Energy kilojoule = mole * kJ_per_mol;

const Dimension::MolarEnergy MJ_per_mol = 1000 * kJ_per_mol;
const Dimension::Energy megajoule = 1000 * kilojoule;

const Dimension::MolarEnergy J_per_mol = 0.001 * kJ_per_mol;
const Dimension::Energy joule = 0.001 * kilojoule;

/** Conversion factor from international kcal mol-1 to internal units  */
const Dimension::MolarEnergy int_kcal_per_mol( 4.1868 * kJ_per_mol );
const Dimension::MolarEnergy int_cal_per_mol( 0.001 * int_kcal_per_mol );

const Dimension::Energy int_kcal( mole * int_kcal_per_mol );
const Dimension::Energy int_cal( 0.001 * int_kcal );

//http://physics.nist.gov/cgi-bin/cuu/Value?hr|search_for=hartree
const Dimension::Energy hartree(4.35974394e-18 * joule);

////////////////////////////////////////////////////////////
// Units of time. Internal unit = akma_time == 48.8882 fs //
////////////////////////////////////////////////////////////

const Dimension::Time akma_time(1);

const Dimension::Time second( std::sqrt( (kg_per_mol * meter * meter) / J_per_mol ) );
const Dimension::Time millisecond = 0.001 * second;
const Dimension::Time microsecond = 0.001 * millisecond;
const Dimension::Time nanosecond = 0.001 * microsecond;
const Dimension::Time picosecond = 0.001 * nanosecond;
const Dimension::Time femtosecond = 0.001 * picosecond;

const Dimension::Time minute( 60 * second );
const Dimension::Time hour( 60 * minute );
const Dimension::Time day( 24 * hour );
const Dimension::Time week( 7 * day );
const Dimension::Time fortnight( 2 * week );

/////////////////////////////////////////////////////////////
// Units of velocity. Internal unit = Angstrom / AKMA time //
/////////////////////////////////////////////////////////////

const Dimension::Velocity akma_velocity(1);

const Dimension::Velocity angstroms_per_fs( angstrom / femtosecond );
const Dimension::Velocity meters_per_second( meter / second );
const Dimension::Velocity kilometers_per_hour( kilometer / hour );
const Dimension::Velocity miles_per_hour( mile / hour );
const Dimension::Velocity mph = miles_per_hour;
const Dimension::Velocity kph = kilometers_per_hour;

//////////////////////////////////////////////////////////
// Units of force. Internal units = kcal mol-1 A-1      //
//////////////////////////////////////////////////////////

/** Convert a force in Newtons to internal units */
const Dimension::Force newton( joule / meter );

/** Weights */
const Dimension::Force ounce = 0.27801385095 * newton;
const Dimension::Force pound = 16 * ounce;
const Dimension::Force stone = 14 * pound;
const Dimension::Force hundredweight = 8 * stone;

//////////////////////////////////////////////////////////
// Units of pressure. Internal units = kcal mol-1 A-2   //
//////////////////////////////////////////////////////////

const Dimension::Pressure pascal = newton / (meter*meter);

const Dimension::Pressure bar = 100000 * pascal;
const Dimension::Pressure atm = 101325 * pascal;

const Dimension::Pressure psi = pound / (inch*inch);
const Dimension::Pressure mmHg = 133.322 * pascal;

//////////////////////////////////////////////////////////
// Units of temperature. Internal units = Kelvin        //
//////////////////////////////////////////////////////////

const Dimension::Temperature kelvin(1);

// other temperature units defined in temperature.h

//////////////////////////////////////////////////////////
// Now some miscellaneous units                         //
//////////////////////////////////////////////////////////

/** Convert the units of current (amps) */
const Dimension::Current amp = coulomb / second;

/** Volts */
const Dimension::Potential volt = joule / coulomb;

/** Convert the units of capacitance (farads) */
const Dimension::Capacitance farad = coulomb / volt;

/** Convert power in Watts */
const Dimension::Power watt = joule / second;
const Dimension::MolarPower watt_per_mol = J_per_mol / second;

///////////////////////////////////////////////////////////
// Now physical constants converted into internal units. //
//  The values of these are taken from the 1998 CODATA   //
//  values - see fundemental_constants.pdf in techdocs   //
///////////////////////////////////////////////////////////

/** Speed of light in a vacuum */
//http://physics.nist.gov/cgi-bin/cuu/Value?c|search_for=c
const Dimension::Velocity c = 299792458 * (meter / second);

/** Epsilon_0 (electrostatic constant) 8.854187817e-12 F m-1 */
//http://physics.nist.gov/cgi-bin/cuu/Value?ep0|search_for=permittivity
const double epsilon0 = 8.854187817e-12 * (farad / meter);

/** 4 * pi * epsilon_0 */
const double four_pi_eps0 = 4.0 * SireMaths::pi * epsilon0;

/** 1 / (4 * pi * epsilon0) */
const double one_over_four_pi_eps0 = 1.0 / four_pi_eps0;

/** Gas constant (8.314472 J mol-1 K-1) */
//http://physics.nist.gov/cgi-bin/cuu/Value?r|search_for=gas
const double gasr = 8.314472 * (J_per_mol / kelvin);

/** Boltzmann constant J K-1 (is equal to gasr in internal units of kcal mol-1 K-1) */
const double k_boltz = gasr;

/** Magnetic constant, mu0, 4pi * 10-7 N A-2 */
const double mu0 = 4.0e-7 * pi * (newton / (amp*amp));

/** Newton's gravitational constant */
//http://physics.nist.gov/cgi-bin/cuu/Value?bg|search_for=gravitational
const double G_newton = 6.67428e-11 * ((meter*meter*meter) / (kilogram * second * second));

/** Acceleration due to gravity on Earth */
const Dimension::Acceleration g_accel = 9.8 * meter / (second*second);

/** Planck's constant */
//http://physics.nist.gov/cgi-bin/cuu/Value?h|search_for=planck
const double h_planck = 6.62606896e-34 * (joule * second);

/** Plank / 2pi */
const double h_slash = h_planck / (2.0*pi);

/** Mass of an electron */
//http://physics.nist.gov/cgi-bin/cuu/Value?me|search_for=mass
const Dimension::Mass electron_mass = 9.10938215e-31 * kilogram;

/** Mass of a proton */
//http://physics.nist.gov/cgi-bin/cuu/Value?mp|search_for=mass
const Dimension::Mass proton_mass = 1.672621637e-27 * kilogram;

/** Mass of a neutron */
//http://physics.nist.gov/cgi-bin/cuu/Value?mn|search_for=mass
const Dimension::Mass neutron_mass = 1.674927211e-27 * kilogram;

/** Atomic mass constant */
//http://physics.nist.gov/cgi-bin/cuu/Value?u|search_for=mass
const Dimension::Mass atomic_mass_constant = 1.660538782e-27 * kilogram;

/** Molar volume of an ideal gas  (273.15 K, 101.325 kPa) */
//http://physics.nist.gov/cgi-bin/cuu/Value?mvolstd|search_for=molar+volume
const Dimension::MolarVolume molar_volume = 22.413996e-3 * (meter*meter*meter) / mole;

}

SIRE_END_HEADER

#endif // end of 'ifndef SKIP_BROKEN_GCCXML_PARTS'

#endif
