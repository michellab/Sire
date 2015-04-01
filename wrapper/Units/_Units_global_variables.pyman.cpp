
#include "_Units_global_variables.pyman.hpp"
#include <boost/python.hpp>
#include "SireUnits/units.h"
#include "SireUnits/temperature.h"

using namespace boost::python;
using namespace SireUnits;
using namespace SireUnits::Dimension;

void register_man_global_variables()
{
    scope().attr("mole") = mole;

    scope().attr("dozen") = dozen;

    scope().attr("radians") = radians;

    scope().attr("radian") = radian;

    scope().attr("degrees") = degrees;

    scope().attr("degree") = degree;

    scope().attr("angle_minute") = angle_minute;

    scope().attr("angle_minutes") = angle_minutes;

    scope().attr("angle_second") = angle_second;

    scope().attr("angle_seconds") = angle_seconds;

    scope().attr("octant") = octant;

    scope().attr("octants") = octants;

    scope().attr("sextant") = sextant;

    scope().attr("sextants") = sextants;

    scope().attr("quadrant") = quadrant;

    scope().attr("quadrants") = quadrants;

    scope().attr("gradian") = gradian;

    scope().attr("gradians") = gradians;

    scope().attr("grad") = grad;

    scope().attr("gon") = gon;

    scope().attr("revolution") = revolution;

    scope().attr("revolutions") = revolutions;

    scope().attr("revs") = revs;

    scope().attr("circumference") = circumference;

    scope().attr("angstrom") = angstrom;

    scope().attr("angstroms") = angstroms;

    scope().attr("picometer") = picometer;

    scope().attr("nanometer") = nanometer;

    scope().attr("micrometer") = micrometer;

    scope().attr("millimeter") = millimeter;

    scope().attr("centimeter") = centimeter;

    scope().attr("meter") = meter;

    scope().attr("kilometer") = kilometer;

    scope().attr("picometers") = picometers;

    scope().attr("nanometers") = nanometers;

    scope().attr("micrometers") = micrometers;

    scope().attr("millimeters") = millimeters;

    scope().attr("centimeters") = centimeters;

    scope().attr("meters") = meters;

    scope().attr("kilometers") = kilometers;

    scope().attr("bohr_radii") = bohr_radii;

    scope().attr("inch") = inch;

    scope().attr("foot") = foot;

    scope().attr("yard") = yard;

    scope().attr("mile") = mile;

    scope().attr("inches") = inches;

    scope().attr("feet") = feet;

    scope().attr("yards") = yards;

    scope().attr("miles") = miles;

    scope().attr("angstrom2") = angstrom2;

    scope().attr("picometer2") = picometer2;

    scope().attr("nanometer2") = nanometer2;

    scope().attr("meter2") = meter2;

    scope().attr("angstrom3") = angstrom3;

    scope().attr("picometer3") = picometer3;

    scope().attr("nanometer3") = nanometer3;

    scope().attr("meter3") = meter3;

    scope().attr("g_per_mol") = g_per_mol;

    scope().attr("gram") = gram;

    scope().attr("kilogram") = kilogram;

    scope().attr("tonne") = tonne;

    scope().attr("milligram") = milligram;

    scope().attr("microgram") = microgram;

    scope().attr("nanogram") = nanogram;

    scope().attr("picogram") = picogram;

    scope().attr("femtogram") = femtogram;

    scope().attr("kg_per_mol") = kg_per_mol;

    scope().attr("tonne_per_mol") = tonne_per_mol;

    scope().attr("mg_per_mol") = mg_per_mol;

    scope().attr("ug_per_mol") = ug_per_mol;

    scope().attr("ng_per_mol") = ng_per_mol;

    scope().attr("pg_per_mol") = pg_per_mol;

    scope().attr("fg_per_mol") = fg_per_mol;

    scope().attr("mod_electron") = mod_electron;

    scope().attr("faraday") = faraday;

    scope().attr("coulomb") = coulomb;

    scope().attr("coulomb_per_mol") = coulomb_per_mol;

    scope().attr("e_charge") = e_charge;

    scope().attr("kcal_per_mol") = kcal_per_mol;

    scope().attr("kcal") = kcal;

    scope().attr("cal_per_mol") = cal_per_mol;

    scope().attr("cal") = cal;

    scope().attr("kJ_per_mol") = kJ_per_mol;

    scope().attr("kilojoule") = kilojoule;

    scope().attr("MJ_per_mol") = MJ_per_mol;

    scope().attr("megajoule") = megajoule;

    scope().attr("J_per_mol") = J_per_mol;

    scope().attr("joule") = joule;

    scope().attr("int_kcal_per_mol") = int_kcal_per_mol;

    scope().attr("int_cal_per_mol") = int_cal_per_mol;

    scope().attr("int_kcal") = int_kcal;

    scope().attr("int_cal") = int_cal;

    scope().attr("hartree") = hartree;

    scope().attr("akma_time") = akma_time;

    scope().attr("second") = second;

    scope().attr("millisecond") = millisecond;

    scope().attr("microsecond") = microsecond;

    scope().attr("nanosecond") = nanosecond;

    scope().attr("picosecond") = picosecond;

    scope().attr("femtosecond") = femtosecond;

    scope().attr("minute") = minute;

    scope().attr("hour") = hour;

    scope().attr("day") = day;

    scope().attr("week") = week;

    scope().attr("fortnight") = fortnight;

    scope().attr("akma_velocity") = akma_velocity;

    scope().attr("angstroms_per_fs") = angstroms_per_fs;

    scope().attr("meters_per_second") = meters_per_second;

    scope().attr("kilometers_per_hour") = kilometers_per_hour;

    scope().attr("miles_per_hour") = miles_per_hour;

    scope().attr("mph") = mph;

    scope().attr("kph") = kph;

    scope().attr("newton") = newton;

    scope().attr("ounce") = ounce;

    scope().attr("pound") = pound;

    scope().attr("stone") = stone;

    scope().attr("hundredweight") = hundredweight;

    scope().attr("pascal") = pascal;

    scope().attr("bar") = bar;

    scope().attr("atm") = atm;

    scope().attr("psi") = psi;

    scope().attr("mmHg") = mmHg;

    scope().attr("kelvin") = kelvin;

    scope().attr("amp") = amp;

    scope().attr("volt") = volt;

    scope().attr("farad") = farad;

    scope().attr("watt") = watt;

    scope().attr("watt_per_mol") = watt_per_mol;

    scope().attr("c") = c;

    scope().attr("epsilon0") = epsilon0;

    scope().attr("four_pi_eps0") = four_pi_eps0;

    scope().attr("one_over_four_pi_eps0") = one_over_four_pi_eps0;

    scope().attr("gasr") = gasr;

    scope().attr("k_boltz") = k_boltz;

    scope().attr("mu0") = mu0;

    scope().attr("G_newton") = G_newton;

    scope().attr("g_accel") = g_accel;

    scope().attr("h_planck") = h_planck;

    scope().attr("h_slash") = h_slash;

    scope().attr("electron_mass") = electron_mass;

    scope().attr("proton_mass") = proton_mass;

    scope().attr("neutron_mass") = neutron_mass;

    scope().attr("atomic_mass_constant") = atomic_mass_constant;

    scope().attr("molar_volume") = molar_volume;

    scope().attr("celsius") = celsius;

    scope().attr("fahrenheit") = fahrenheit;

}

