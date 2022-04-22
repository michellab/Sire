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

#include "SireBase/slice.h"

#include "SireError/errors.h"

#include "Helpers/release_gil_policy.hpp"

#include <QDebug>

using namespace boost::python;

/** This function convert a python string or unicode to a QChar */
void Slice_from_python_slice(PyObject* obj_ptr,
                             converter::rvalue_from_python_stage1_data* data)
{
    auto raii = boost::python::release_gil_policy::acquire_gil();
    if ( PySlice_Check(obj_ptr) )
    {
        Py_ssize_t start, stop, step;
        int ok = PySlice_Unpack(obj_ptr, &start, &stop, &step);

        if (ok != 0)
            throw SireError::assertation_failed(
                QObject::tr("Cannot unpack the slice! %1").arg(ok), CODELOC );

        void* storage = ((converter::rvalue_from_python_storage<SireBase::Slice>*) data)->storage.bytes;

        if (stop == 9223372036854775807L)
        {
            // magic number used when the end value has not been set
            new (storage) SireBase::Slice( SireBase::Slice::fromStart(start, step));
        }
        else
        {
            if (stop == -9223372036854775808L)
            {
                stop = 0;
            }

            new (storage) SireBase::Slice( SireBase::Slice::fromStartStop(start, stop, step));
        }

        data->convertible = storage;
    }
}

/** The actual struct used to control the conversion */
struct Slice_from_python
{
    Slice_from_python()
    {
        converter::registry::push_back(  &convertible,
                                         &construct,
                                         type_id<SireBase::Slice>() );
    }

    /** Can the python object pointed to by 'obj_ptr' be converted
        to a Slice?
    */
    static void* convertible(PyObject* obj_ptr)
    {
        auto raii = boost::python::release_gil_policy::acquire_gil();
        if ( PySlice_Check(obj_ptr) )
        {
             return obj_ptr;
        }
        else
            return 0;
    }

    /** Perform the actual conversion */
    static void construct(  PyObject* obj_ptr,
                            converter::rvalue_from_python_stage1_data* data)
    {
        Slice_from_python_slice(obj_ptr, data);
    }
};

void autoconvert_Slice()
{
    //code to get a Slice from a python object
    Slice_from_python();
}
