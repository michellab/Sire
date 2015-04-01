/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2007  Christopher Woods
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

#include <Python.h>
#include <boost/python.hpp>

#include "SireUnits/units.h"
#include "SireUnits/dimensions.h"
#include "SireUnits/temperature.h"

#include "generalunit.h"

#include "SireMM/distancerestraint.h"

using namespace SireUnits;
using namespace SireUnits::Dimension;

using namespace boost::python;

void register_SireUnits_dimensions()
{
    register_dimension< Dimensionless >();
    register_dimension< Mass >();
    register_dimension< MolarMass >();
    register_dimension< Length >();
    register_dimension< Time >();
    register_dimension< Charge >();
    register_dimension< MolarCharge >();
    register_dimension< Temperature >();
    register_dimension< Quantity >();
    register_dimension< Angle >();
    register_dimension< Area >();
    register_dimension< Volume >();
    register_dimension< MolarVolume >();
    register_dimension< Velocity >();
    register_dimension< Acceleration >();
    register_dimension< AngularVelocity >();
    register_dimension< AngularAcceleration >();
    register_dimension< Energy >();
    register_dimension< MolarEnergy >();
    register_dimension< Power >();
    register_dimension< MolarPower >();
    register_dimension< Density >();
    register_dimension< MolarDensity >();
    register_dimension< Force >();
    register_dimension< Pressure >();
    register_dimension< Current >();
    register_dimension< Potential >();
    register_dimension< Capacitance >();

    register_dimension< SireMM::HarmonicDistanceForceConstant >();

}
