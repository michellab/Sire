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

#ifndef PYWRAP_SIREPY_CONVERTSET_HPP
#define PYWRAP_SIREPY_CONVERTSET_HPP

#include "sireglobal.h"

#include <boost/python.hpp>
#include <boost/tuple/tuple.hpp>

namespace bp = boost::python;

SIRE_BEGIN_HEADER

/** This struct provides the from-Python conversion from a list or
    tuple to a list-like container of type 'C' (e.g. QSet) */
template<class C, class T>
struct from_py_set
{
    /** Constructor - register the conversion functions
        for this type */
    from_py_set()
    {
        bp::converter::registry::push_back(
            &convertible,
            &construct,
            bp::type_id< C >());
    }

    /** Test whether or not it is possible to convert the PyObject
        to a QVector where all of the elements are of type 'T' */
    static void* convertible(PyObject* obj_ptr)
    {
        //is this a tuple type?
        if ( PyTuple_Check(obj_ptr) )
        {
            //check the tuple elements... - convert to a boost::tuple object
            bp::tuple t( bp::handle<>(bp::borrowed(obj_ptr)) );

            //how many elements are there?
            int n = PyTuple_Size(obj_ptr);

            //can they all be converted to type 'T'?
            for (int i=0; i<n; ++i)
            {
                if (not bp::extract<T>(t[i]).check())
                    return 0;
            }

            //the tuple is ok!
            return obj_ptr;
        }
        //is this a list type?
        else if ( PyList_Check(obj_ptr) )
        {
            //check that all of the list elements can be converted to the right type
            bp::list l( bp::handle<>(bp::borrowed(obj_ptr)) );

            //how many elements are there?
            int n = PyList_Size(obj_ptr);

            //can all of the elements be converted to type 'T'?
            for (int i=0; i<n; ++i)
            {
                if (not bp::extract<T>(l[i]).check())
                    return 0;
            }

            //the list is ok!
            return obj_ptr;
        }
        else
            //could not recognise the type...
            return 0;
    }

    /** Construct a container of type T from the PyObject pointed to
        by 'obj_ptr' */
    static void construct(
        PyObject* obj_ptr,
        bp::converter::rvalue_from_python_stage1_data* data)
    {
        // need to re-acquire the GIL when creating new list objects
        PyGILState_STATE gstate;
        gstate = PyGILState_Ensure();

        if (PyTuple_Check(obj_ptr))
        {
            //convert the PyObject to a boost::python::object
            bp::tuple t( bp::handle<>(bp::borrowed(obj_ptr)) );

            //locate the storage space for the result
            void* storage =
                ( (bp::converter::rvalue_from_python_storage<C>*)data )->storage.bytes;

            //create the T container
            new (storage) C();

            C *container = static_cast<C*>(storage);

            //add all of the elements from the tuple - get the number of elements in the tuple
            int n = PyTuple_Size(obj_ptr);

            for (int i=0; i<n; ++i)
            {
                container->insert( bp::extract<T>(t[i])() );
            }

            data->convertible = storage;
        }
        else if (PyList_Check(obj_ptr))
        {
            //convert the PyObject to a boost::python::object
            bp::list l( bp::handle<>(bp::borrowed(obj_ptr)) );

            //locate the storage space for the result
            void* storage =
                ( (bp::converter::rvalue_from_python_storage<C>*)data )->storage.bytes;

            //create the T container
            new (storage) C();

            C *container = static_cast<C*>(storage);

            //add all of the elements from the tuple - get the number of elements in the tuple
            int n = PyList_Size(obj_ptr);

            for (int i=0; i<n; ++i)
            {
                container->insert( bp::extract<T>(l[i])() );
            }

            data->convertible = storage;
        }

        PyGILState_Release(gstate);
    }
};

template<class C>
struct to_py_set
{
    static PyObject* convert(const C &cpp_set)
    {
        bp::list python_set;

        //add all items to the python dictionary
        for (typename C::const_iterator it = cpp_set.begin();
             it != cpp_set.end();
             ++it)
        {
            python_set.append(*it);
        }

        return bp::incref( python_set.ptr() );
    }
};

template<class C>
void register_set()
{
    bp::to_python_converter< C, to_py_set<C> >();

    typedef typename C::value_type value_type;

    bp::converter::registry::push_back( &from_py_set<C,value_type>::convertible,
                                        &from_py_set<C,value_type>::construct,
                                        bp::type_id<C>() );
}

template<class C, class value_type>
void register_set()
{
    bp::to_python_converter< C, to_py_set<C> >();

    bp::converter::registry::push_back( &from_py_set<C,value_type>::convertible,
                                        &from_py_set<C,value_type>::construct,
                                        bp::type_id<C>() );
}

SIRE_END_HEADER

#endif
