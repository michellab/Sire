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

#include <Python.h>

#include <boost/python.hpp>

#include "generalunit.h"

#include "SireError/errors.h"

#include "Helpers/release_gil_policy.hpp"

#include <QDebug>

using namespace boost::python;

/** This function converts a GeneralUnitProperty to a python GeneralUnit

*/
struct GeneralUnitProperty_to_python
{
    static PyObject* convert(const SireUnits::Dimension::GeneralUnitProperty &prop)
    {
        qDebug() << "CONVERT!";
        auto raii = boost::python::release_gil_policy::acquire_gil();
        return bp::incref(bp::object(prop.value()).ptr());
    }
};

void autoconvert_GeneralUnitProperty()
{
    qDebug() << "REGISTER";
    //code to get a Python object from a GeneralUnitProperty
    bp::to_python_converter< SireUnits::Dimension::GeneralUnitProperty, GeneralUnitProperty_to_python >();
}
