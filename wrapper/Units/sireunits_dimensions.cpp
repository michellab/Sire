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

#include "SireUnits/generalunit.h"

#include "SireMM/distancerestraint.h"

using namespace SireUnits;
using namespace SireUnits::Dimension;

using namespace boost::python;

template<class D>
struct from_general_unit
{
    /** Constructor - register the conversion functions
        for this type */
    from_general_unit()
    {
        boost::python::converter::registry::push_back(
            &convertible,
            &construct,
            type_id< D >());
    }

    /** Test whether or not it is possible to convert the PyObject
        holding a GeneralUnit to this specific PhysUnit*/
    static void* convertible(PyObject* obj_ptr)
    {
        object obj( handle<>(borrowed(obj_ptr)) );

        extract<const Dimension::GeneralUnit&> x(obj);

        //is this a GeneralUnit?
        if ( x.check() )
        {
            //it is ;-)  Get a reference to it and make sure
            //that it is of the right dimension
            const Dimension::GeneralUnit &gen_unit = x();

            if (gen_unit.isZero())
            {
                // zero is compatible with everything
                return obj_ptr;
            }
            else if ( gen_unit.MASS() == D::MASS() and
                 gen_unit.LENGTH() == D::LENGTH() and
                 gen_unit.TIME() == D::TIME() and
                 gen_unit.CHARGE() == D::CHARGE() and
                 gen_unit.TEMPERATURE() == D::TEMPERATURE() and
                 gen_unit.QUANTITY() == D::QUANTITY() and
                 gen_unit.ANGLE() == D::ANGLE() )
            {
                //this has the right dimension :-)
                return obj_ptr;
            }
        }

        //could not recognise the type or the dimension was wrong
        return 0;
    }

    /** Construct a PhysUnit from the passed GeneralUnit */
    static void construct(
        PyObject* obj_ptr,
        boost::python::converter::rvalue_from_python_stage1_data* data)
    {
        object obj( handle<>(borrowed(obj_ptr)) );

        extract<const Dimension::GeneralUnit&> x(obj);

        //is this a GeneralUnit?
        if ( x.check() )
        {
            //it is ;-)  Get a reference to it and make sure
            //that it is of the right dimension
            const Dimension::GeneralUnit &gen_unit = x();

            if ( gen_unit.isZero() )
            {
                // zero is compatible with everything
                void* storage =
                    ( (converter::rvalue_from_python_storage<D>*)data )->storage.bytes;

                //create the T container
                new (storage) D(0.0);

                data->convertible = storage;
            }
            else if ( gen_unit.MASS() == D::MASS() and
                 gen_unit.LENGTH() == D::LENGTH() and
                 gen_unit.TIME() == D::TIME() and
                 gen_unit.CHARGE() == D::CHARGE() and
                 gen_unit.TEMPERATURE() == D::TEMPERATURE() and
                 gen_unit.QUANTITY() == D::QUANTITY() and
                 gen_unit.ANGLE() == D::ANGLE() )
            {
                //locate the storage space for the result
                void* storage =
                    ( (converter::rvalue_from_python_storage<D>*)data )->storage.bytes;

                //create the T container
                new (storage) D( gen_unit.scaleFactor() );

                data->convertible = storage;
            }
        }
    }
};

template<class D>
struct to_general_unit
{
    static PyObject* convert(const D &unit)
    {
        return incref( object(Dimension::GeneralUnit(unit)).ptr() );
    }
};

template<class D>
void register_dimension()
{
    to_python_converter< D, to_general_unit<D> >();

    converter::registry::push_back(
          &from_general_unit<D>::convertible,
          &from_general_unit<D>::construct,
          type_id<D>() );
}

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
