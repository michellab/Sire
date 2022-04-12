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

using namespace boost::python;

/** This function convert a python string or unicode to a QChar */
void Slice_from_python_slice(PyObject* obj_ptr,
                             converter::rvalue_from_python_stage1_data* data)
{
    Py_ssize_t start, stop, step;
    int ok = PySlice_Unpack(obj_ptr, &start, &stop, &step);

    if (not ok)
        boost::python::throw_error_already_set();

    void* storage = ((converter::rvalue_from_python_storage<SireBase::Slice>*) data)->storage.bytes;

    new (storage) SireBase::Slice( SireBase::Slice::fromStartStop(start, stop, step));
    data->convertible = storage;
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
