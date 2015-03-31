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

#ifndef SIREMOL_ELEMENT_DATA_H
#define SIREMOL_ELEMENT_DATA_H

void ElementDB::populate()
{
this->import( new ElementData(0, QObject::tr("dummy"), "Xx", 0, 0, 0.000, 0.000, 0.000, 0, 0.000, 0.00, 0.07, 0.50, 0.70) );
this->import( new ElementData(1, QObject::tr("Hydrogen"), "H", 1, 1, 0.230, 0.330, 1.200, 1, 1.00794, 2.10, 1.00, 1.00, 1.00) );
this->import( new ElementData(2, QObject::tr("Helium"), "He", 18, 1, 0.930, 0.700, 1.400, 0, 4.002602, 0.00, 0.85, 1.00, 1.00) );
this->import( new ElementData(3, QObject::tr("Lithium"), "Li", 1, 2, 0.680, 1.230, 1.820, 1, 6.941, 0.98, 0.80, 0.50, 1.00) );
this->import( new ElementData(4, QObject::tr("Beryllium"), "Be", 2, 2, 0.350, 0.900, 1.700, 2, 9.012182, 1.57, 0.76, 1.00, 0.00) );
this->import( new ElementData(5, QObject::tr("Boron"), "B", 13, 2, 0.830, 0.820, 2.080, 3, 10.811, 2.04, 1.00, 0.71, 0.71) );
this->import( new ElementData(6, QObject::tr("Carbon"), "C", 14, 2, 0.680, 0.770, 1.950, 4, 12.0107, 2.55, 0.50, 0.50, 0.50) );
this->import( new ElementData(7, QObject::tr("Nitrogen"), "N", 15, 2, 0.680, 0.700, 1.850, 4, 14.0067, 3.04, 0.05, 0.05, 1.00) );
this->import( new ElementData(8, QObject::tr("Oxygen"), "O", 16, 2, 0.680, 0.660, 1.700, 2, 15.9994, 3.44, 1.00, 0.05, 0.05) );
this->import( new ElementData(9, QObject::tr("Fluorine"), "F", 17, 2, 0.640, 0.611, 1.730, 1, 18.9984032, 3.98, 0.70, 1.00, 1.00) );
this->import( new ElementData(10, QObject::tr("Neon"), "Ne", 18, 2, 1.120, 0.700, 1.540, 0, 20.1797, 0.00, 0.70, 0.89, 0.96) );
this->import( new ElementData(11, QObject::tr("Sodium"), "Na", 1, 3, 0.970, 1.540, 2.270, 1, 22.989770, 0.93, 0.67, 0.36, 0.95) );
this->import( new ElementData(12, QObject::tr("Magnesium"), "Mg", 2, 3, 1.100, 1.360, 1.730, 2, 24.3050, 1.31, 0.54, 1.00, 0.00) );
this->import( new ElementData(13, QObject::tr("Aluminium"), "Al", 13, 3, 1.350, 1.180, 2.050, 6, 26.981538, 1.61, 0.75, 0.65, 0.65) );
this->import( new ElementData(14, QObject::tr("Silicon"), "Si", 14, 3, 1.200, 0.937, 2.100, 6, 28.0855, 1.90, 0.50, 0.60, 0.60) );
this->import( new ElementData(15, QObject::tr("Phosphorus"), "P", 15, 3, 1.100, 0.890, 2.080, 5, 30.973761, 2.19, 1.00, 0.50, 0.00) );
this->import( new ElementData(16, QObject::tr("Sulfur"), "S", 16, 3, 1.020, 1.040, 2.000, 6, 32.065, 2.58, 1.00, 1.00, 0.19) );
this->import( new ElementData(17, QObject::tr("Chlorine"), "Cl", 17, 3, 0.990, 0.997, 1.970, 1, 35.453, 3.16, 0.12, 0.94, 0.12) );
this->import( new ElementData(18, QObject::tr("Argon"), "Ar", 18, 3, 1.570, 1.740, 1.880, 0, 39.948, 0.00, 0.50, 0.82, 0.89) );
this->import( new ElementData(19, QObject::tr("Potassium"), "K", 1, 4, 1.330, 2.030, 2.750, 1, 39.0983, 0.82, 0.56, 0.25, 0.83) );
this->import( new ElementData(20, QObject::tr("Calcium"), "Ca", 2, 4, 0.990, 1.740, 1.973, 2, 40.078, 1.00, 0.24, 1.00, 0.00) );
this->import( new ElementData(21, QObject::tr("Scandium"), "Sc", 3, 4, 1.440, 1.440, 1.700, 6, 44.955910, 1.36, 0.90, 0.90, 0.90) );
this->import( new ElementData(22, QObject::tr("Titanium"), "Ti", 4, 4, 1.470, 1.320, 1.700, 6, 47.867, 1.54, 0.75, 0.76, 0.78) );
this->import( new ElementData(23, QObject::tr("Vanadium"), "V", 5, 4, 1.330, 1.220, 1.700, 6, 50.9415, 1.63, 0.65, 0.65, 0.67) );
this->import( new ElementData(24, QObject::tr("Chromium"), "Cr", 6, 4, 1.350, 1.180, 1.700, 6, 51.9961, 1.66, 0.54, 0.60, 0.78) );
this->import( new ElementData(25, QObject::tr("Manganese"), "Mn", 7, 4, 1.350, 1.170, 1.700, 8, 54.938049, 1.55, 0.61, 0.48, 0.78) );
this->import( new ElementData(26, QObject::tr("Iron"), "Fe", 8, 4, 1.340, 1.170, 1.700, 6, 55.845, 1.83, 0.50, 0.48, 0.78) );
this->import( new ElementData(27, QObject::tr("Cobalt"), "Co", 9, 4, 1.330, 1.160, 1.700, 6, 58.933200, 1.88, 0.44, 0.48, 0.78) );
this->import( new ElementData(28, QObject::tr("Nickel"), "Ni", 10, 4, 1.500, 1.150, 1.630, 6, 58.6934, 1.91, 0.36, 0.48, 0.76) );
this->import( new ElementData(29, QObject::tr("Copper"), "Cu", 11, 4, 1.520, 1.170, 1.400, 6, 63.546, 1.90, 1.00, 0.48, 0.38) );
this->import( new ElementData(30, QObject::tr("Zinc"), "Zn", 12, 4, 1.450, 1.250, 1.390, 6, 65.409, 1.65, 0.49, 0.50, 0.69) );
this->import( new ElementData(31, QObject::tr("Gallium"), "Ga", 13, 4, 1.220, 1.260, 1.870, 3, 69.723, 1.81, 0.76, 0.56, 0.56) );
this->import( new ElementData(32, QObject::tr("Germanium"), "Ge", 14, 4, 1.170, 1.188, 1.700, 4, 72.64, 2.01, 0.40, 0.56, 0.56) );
this->import( new ElementData(33, QObject::tr("Arsenic"), "As", 15, 4, 1.210, 1.200, 1.850, 3, 74.92160, 2.18, 0.74, 0.50, 0.89) );
this->import( new ElementData(34, QObject::tr("Selenium"), "Se", 16, 4, 1.220, 1.170, 1.900, 2, 78.96, 2.55, 1.00, 0.63, 0.00) );
this->import( new ElementData(35, QObject::tr("Bromine"), "Br", 17, 4, 1.210, 1.167, 2.100, 1, 79.904, 2.96, 0.65, 0.16, 0.16) );
this->import( new ElementData(36, QObject::tr("Krypton"), "Kr", 18, 4, 1.910, 1.910, 2.020, 0, 83.798, 0.00, 0.36, 0.72, 0.82) );
this->import( new ElementData(37, QObject::tr("Rubidium"), "Rb", 1, 5, 1.470, 2.160, 1.700, 1, 85.4678, 0.82, 0.44, 0.18, 0.69) );
this->import( new ElementData(38, QObject::tr("Strontium"), "Sr", 2, 5, 1.120, 1.910, 1.700, 2, 87.62, 0.95, 0.00, 1.00, 0.00) );
this->import( new ElementData(39, QObject::tr("Yttrium"), "Y", 3, 5, 1.780, 1.620, 1.700, 6, 88.90585, 1.22, 0.58, 1.00, 1.00) );
this->import( new ElementData(40, QObject::tr("Zirconium"), "Zr", 4, 5, 1.560, 1.450, 1.700, 6, 91.224, 1.33, 0.58, 0.88, 0.88) );
this->import( new ElementData(41, QObject::tr("Niobium"), "Nb", 5, 5, 1.480, 1.340, 1.700, 6, 92.90638, 1.6, 0.45, 0.76, 0.79) );
this->import( new ElementData(42, QObject::tr("Molybdenum"), "Mo", 6, 5, 1.470, 1.300, 1.700, 6, 95.94, 2.16, 0.33, 0.71, 0.71) );
this->import( new ElementData(43, QObject::tr("Technetium"), "Tc", 7, 5, 1.350, 1.270, 1.700, 6, 98, 1.9, 0.23, 0.62, 0.62) );
this->import( new ElementData(44, QObject::tr("Ruthenium"), "Ru", 8, 5, 1.400, 1.250, 1.700, 6, 101.07, 2.2, 0.14, 0.56, 0.56) );
this->import( new ElementData(45, QObject::tr("Rhodium"), "Rh", 9, 5, 1.450, 1.250, 1.700, 6, 102.90550, 2.28, 0.04, 0.49, 0.55) );
this->import( new ElementData(46, QObject::tr("Palladium"), "Pd", 10, 5, 1.500, 1.280, 1.630, 6, 106.42, 2.20, 0.00, 0.41, 0.52) );
this->import( new ElementData(47, QObject::tr("Silver"), "Ag", 11, 5, 1.590, 1.340, 1.720, 6, 107.8682, 1.93, 0.88, 0.88, 1.00) );
this->import( new ElementData(48, QObject::tr("Cadmium"), "Cd", 12, 5, 1.690, 1.480, 1.580, 6, 112.411, 1.69, 1.00, 0.85, 0.56) );
this->import( new ElementData(49, QObject::tr("Indium"), "In", 13, 5, 1.630, 1.440, 1.930, 3, 114.818, 1.78, 0.65, 0.46, 0.45) );
this->import( new ElementData(50, QObject::tr("Tin"), "Sn", 14, 5, 1.460, 1.385, 2.170, 4, 118.710, 1.96, 0.40, 0.50, 0.50) );
this->import( new ElementData(51, QObject::tr("Antimony"), "Sb", 15, 5, 1.460, 1.400, 2.200, 3, 121.760, 2.05, 0.62, 0.39, 0.71) );
this->import( new ElementData(52, QObject::tr("Tellurium"), "Te", 16, 5, 1.470, 1.378, 2.060, 2, 127.60, 2.1, 0.83, 0.48, 0.00) );
this->import( new ElementData(53, QObject::tr("Iodine"), "I", 17, 5, 1.400, 1.387, 2.150, 1, 126.90447, 2.66, 0.58, 0.00, 0.58) );
this->import( new ElementData(54, QObject::tr("Xenon"), "Xe", 18, 5, 1.980, 1.980, 2.160, 0, 131.293, 2.6, 0.26, 0.62, 0.69) );
this->import( new ElementData(55, QObject::tr("Caesium"), "Cs", 1, 6, 1.670, 2.350, 1.700, 1, 132.90545, 0.79, 0.34, 0.09, 0.56) );
this->import( new ElementData(56, QObject::tr("Barium"), "Ba", 2, 6, 1.340, 1.980, 1.700, 2, 137.327, 0.89, 0.00, 0.79, 0.00) );
this->import( new ElementData(57, QObject::tr("Lanthanum"), "La", 3, 6, 1.870, 1.690, 1.700, 12, 138.9055, 1.10, 0.44, 0.83, 1.00) );
this->import( new ElementData(58, QObject::tr("Cerium"), "Ce", 0, 6, 1.830, 1.830, 1.700, 6, 140.116, 1.12, 1.00, 1.00, 0.78) );
this->import( new ElementData(59, QObject::tr("Praseodymium"), "Pr", 0, 6, 1.820, 1.820, 1.700, 6, 140.90765, 1.13, 0.85, 1.00, 0.78) );
this->import( new ElementData(60, QObject::tr("Neodymium"), "Nd", 0, 6, 1.810, 1.810, 1.700, 6, 144.24, 1.14, 0.78, 1.00, 0.78) );
this->import( new ElementData(61, QObject::tr("Promethium"), "Pm", 0, 6, 1.800, 1.800, 1.700, 6, 145, 1.13, 0.64, 1.00, 0.78) );
this->import( new ElementData(62, QObject::tr("Samarium"), "Sm", 0, 6, 1.800, 1.800, 1.700, 6, 150.36, 1.17, 0.56, 1.00, 0.78) );
this->import( new ElementData(63, QObject::tr("Europium"), "Eu", 0, 6, 1.990, 1.990, 1.700, 6, 151.964, 1.2, 0.38, 1.00, 0.78) );
this->import( new ElementData(64, QObject::tr("Gadolinium"), "Gd", 0, 6, 1.790, 1.790, 1.700, 6, 157.25, 1.20, 0.27, 1.00, 0.78) );
this->import( new ElementData(65, QObject::tr("Terbium"), "Tb", 0, 6, 1.760, 1.760, 1.700, 6, 158.92534, 1.1, 0.19, 1.00, 0.78) );
this->import( new ElementData(66, QObject::tr("Dysprosium"), "Dy", 0, 6, 1.750, 1.750, 1.700, 6, 162.500, 1.22, 0.12, 1.00, 0.78) );
this->import( new ElementData(67, QObject::tr("Holmium"), "Ho", 0, 6, 1.740, 1.740, 1.700, 6, 164.93032, 1.23, 0.00, 1.00, 0.61) );
this->import( new ElementData(68, QObject::tr("Erbium"), "Er", 0, 6, 1.730, 1.730, 1.700, 6, 167.259, 1.24, 0.00, 0.90, 0.46) );
this->import( new ElementData(69, QObject::tr("Thulium"), "Tm", 0, 6, 1.720, 1.720, 1.700, 6, 168.93421, 1.25, 0.00, 0.83, 0.32) );
this->import( new ElementData(70, QObject::tr("Ytterbium"), "Yb", 0, 6, 1.940, 1.940, 1.700, 6, 173.04, 1.1, 0.00, 0.75, 0.22) );
this->import( new ElementData(71, QObject::tr("Lutetium"), "Lu", 0, 6, 1.720, 1.720, 1.700, 6, 174.967, 1.27, 0.00, 0.67, 0.14) );
this->import( new ElementData(72, QObject::tr("Hafnium"), "Hf", 4, 6, 1.570, 1.440, 1.700, 6, 178.49, 1.3, 0.30, 0.76, 1.00) );
this->import( new ElementData(73, QObject::tr("Tantalum"), "Ta", 5, 6, 1.430, 1.340, 1.700, 6, 180.9479, 1.5, 0.30, 0.65, 1.00) );
this->import( new ElementData(74, QObject::tr("Tungsten"), "W", 6, 6, 1.370, 1.300, 1.700, 6, 183.84, 2.36, 0.13, 0.58, 0.84) );
this->import( new ElementData(75, QObject::tr("Rhenium"), "Re", 7, 6, 1.350, 1.280, 1.700, 6, 186.207, 1.9, 0.15, 0.49, 0.67) );
this->import( new ElementData(76, QObject::tr("Osmium"), "Os", 8, 6, 1.370, 1.260, 1.700, 6, 190.23, 2.2, 0.15, 0.40, 0.59) );
this->import( new ElementData(77, QObject::tr("Iridium"), "Ir", 9, 6, 1.320, 1.270, 1.700, 6, 192.217, 2.20, 0.09, 0.33, 0.53) );
this->import( new ElementData(78, QObject::tr("Platinum"), "Pt", 10, 6, 1.500, 1.300, 1.720, 6, 195.078, 2.28, 0.96, 0.93, 0.82) );
this->import( new ElementData(79, QObject::tr("Gold"), "Au", 11, 6, 1.500, 1.340, 1.660, 6, 196.96655, 2.54, 0.80, 0.82, 0.12) );
this->import( new ElementData(80, QObject::tr("Mercury"), "Hg", 12, 6, 1.700, 1.490, 1.550, 6, 200.59, 2.00, 0.71, 0.71, 0.76) );
this->import( new ElementData(81, QObject::tr("Thallium"), "Tl", 13, 6, 1.550, 1.480, 1.960, 3, 204.3833, 2.04, 0.65, 0.33, 0.30) );
this->import( new ElementData(82, QObject::tr("Lead"), "Pb", 14, 6, 1.540, 1.480, 2.020, 4, 207.2, 2.33, 0.34, 0.35, 0.38) );
this->import( new ElementData(83, QObject::tr("Bismuth"), "Bi", 15, 6, 1.540, 1.450, 1.700, 3, 208.98038, 2.02, 0.62, 0.31, 0.71) );
this->import( new ElementData(84, QObject::tr("Polonium"), "Po", 16, 6, 1.680, 1.460, 1.700, 2, 209, 2.0, 0.67, 0.36, 0.00) );
this->import( new ElementData(85, QObject::tr("Astatine"), "At", 17, 6, 1.700, 1.450, 1.700, 1, 210, 2.2, 0.46, 0.31, 0.27) );
this->import( new ElementData(86, QObject::tr("Radon"), "Rn", 18, 6, 2.400, 2.400, 1.700, 0, 222, 0.0, 0.26, 0.51, 0.59) );
this->import( new ElementData(87, QObject::tr("Francium"), "Fr", 1, 7, 2.000, 2.000, 1.700, 1, 223, 0.7, 0.26, 0.00, 0.40) );
this->import( new ElementData(88, QObject::tr("Radium"), "Ra", 2, 7, 1.900, 1.900, 1.700, 2, 226, 0.89, 0.00, 0.49, 0.00) );
this->import( new ElementData(89, QObject::tr("Actinium"), "Ac", 3, 7, 1.880, 1.880, 1.700, 6, 227, 1.1, 0.44, 0.67, 0.98) );
this->import( new ElementData(90, QObject::tr("Thorium"), "Th", 0, 7, 1.790, 1.790, 1.700, 6, 232.0381, 1.3, 0.00, 0.73, 1.00) );
this->import( new ElementData(91, QObject::tr("Protactinium"), "Pa", 0, 7, 1.610, 1.610, 1.700, 6, 231.03588, 1.5, 0.00, 0.63, 1.00) );
this->import( new ElementData(92, QObject::tr("Uranium"), "U", 0, 7, 1.580, 1.580, 1.860, 6, 238.02891, 1.38, 0.00, 0.56, 1.00) );
this->import( new ElementData(93, QObject::tr("Neptunium"), "Np", 0, 7, 1.550, 1.550, 1.700, 6, 237, 1.36, 0.00, 0.50, 1.00) );
this->import( new ElementData(94, QObject::tr("Plutonium"), "Pu", 0, 7, 1.530, 1.530, 1.700, 6, 244, 1.28, 0.00, 0.42, 1.00) );
this->import( new ElementData(95, QObject::tr("Americium"), "Am", 0, 7, 1.510, 1.070, 1.700, 6, 243, 1.3, 0.33, 0.36, 0.95) );
this->import( new ElementData(96, QObject::tr("Curium"), "Cm", 0, 7, 1.500, 0.000, 1.700, 6, 247, 1.3, 0.47, 0.36, 0.89) );
this->import( new ElementData(97, QObject::tr("Berkelium"), "Bk", 0, 7, 1.500, 0.000, 1.700, 6, 247, 1.3, 0.54, 0.31, 0.89) );
this->import( new ElementData(98, QObject::tr("Californium"), "Cf", 0, 7, 1.500, 0.000, 1.700, 6, 251, 1.3, 0.63, 0.21, 0.83) );
this->import( new ElementData(99, QObject::tr("Einsteinium"), "Es", 0, 7, 1.500, 0.000, 1.700, 6, 252, 1.3, 0.70, 0.12, 0.83) );
this->import( new ElementData(100, QObject::tr("Fermium"), "Fm", 0, 7, 1.500, 0.000, 1.700, 6, 257, 1.3, 0.70, 0.12, 0.73) );
this->import( new ElementData(101, QObject::tr("Mendelevium"), "Md", 0, 7, 1.500, 0.000, 1.700, 6, 258, 1.3, 0.70, 0.05, 0.65) );
this->import( new ElementData(102, QObject::tr("Nobelium"), "No", 0, 7, 1.500, 0.000, 1.700, 6, 259, 1.3, 0.74, 0.05, 0.53) );
this->import( new ElementData(103, QObject::tr("Lawrencium"), "Lr", 0, 7, 1.500, 0.000, 1.700, 6, 262, 1.3, 0.78, 0.00, 0.40) );
this->import( new ElementData(104, QObject::tr("Rutherfordium"), "Rf", 4, 7, 1.600, 0.000, 1.700, 6, 261, 0.0, 0.80, 0.00, 0.35) );
this->import( new ElementData(105, QObject::tr("Dubnium"), "Db", 5, 7, 1.600, 0.000, 1.700, 6, 262, 0.0, 0.82, 0.00, 0.31) );
this->import( new ElementData(106, QObject::tr("Seaborgium"), "Sg", 6, 7, 1.600, 0.000, 1.700, 6, 263, 0.0, 0.85, 0.00, 0.27) );
this->import( new ElementData(107, QObject::tr("Bohrium"), "Bh", 7, 7, 1.600, 0.000, 1.700, 6, 264, 0.0, 0.88, 0.00, 0.22) );
this->import( new ElementData(108, QObject::tr("Hassium"), "Hs", 8, 7, 1.600, 0.000, 1.700, 6, 265, 0.0, 0.90, 0.00, 0.18) );
this->import( new ElementData(109, QObject::tr("Meitnerium"), "Mt", 9, 7, 1.600, 0.000, 1.700, 6, 268, 0.0, 0.92, 0.00, 0.15) );
this->import( new ElementData(110, QObject::tr("Darmstadtium"), "Ds", 10, 7, 1.600, 0.000, 1.700, 6, 271, 0.0, 0.95, 0.00, 0.10) );
this->import( new ElementData(111, QObject::tr("Roentgenium"), "Rg", 11, 7, 1.600, 0.000, 1.700, 6, 272, 0.0, 0.98, 0.00, 0.02) );
}

#endif

