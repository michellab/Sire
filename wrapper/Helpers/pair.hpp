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

#ifndef PYWRAP_SIREPY_PAIR_HPP
#define PYWRAP_SIREPY_PAIR_HPP

#include "sireglobal.h"

#include <boost/python.hpp>
#include <utility>

#include "Helpers/release_gil_policy.hpp"

namespace bp = boost::python;

SIRE_BEGIN_HEADER

/** This struct provides the from-Python conversion from a list or
    tuple to a pair of type std::pair<S,T> */
template<class S, class T>
struct from_py_pair
{
    /** Constructor - register the conversion functions
        for this type */
    from_py_pair()
    {
        boost::python::converter::registry::push_back(
            &convertible,
            &construct,
            bp::type_id< std::pair<S,T> >());
    }

    /** Test whether or not it is possible to convert the PyObject
        to a QVector where all of the elements are of type 'T' */
    static void* convertible(PyObject* obj_ptr)
    {
        auto raii = boost::python::release_gil_policy::acquire_gil();

        //is this a tuple type?
        if ( PyTuple_Check(obj_ptr) )
        {
            //check the tuple elements... - convert to a boost::tuple object
            bp::tuple t( bp::handle<>(bp::borrowed(obj_ptr)) );

            //how many elements are there?
            int n = PyTuple_Size(obj_ptr);

            if (n != 2)
                return 0;

            if (not bp::extract<S>(t[0]).check())
                return 0;

            if (not bp::extract<T>(t[1]).check())
                return 0;

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

            if (n != 2)
                return 0;

            if (not bp::extract<S>(l[0]).check())
                return 0;

            if (not bp::extract<T>(l[0]).check())
                return 0;

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
        auto raii = boost::python::release_gil_policy::acquire_gil();

        if (PyTuple_Check(obj_ptr))
        {
            //convert the PyObject to a boost::python::object
            bp::tuple t( bp::handle<>(bp::borrowed(obj_ptr)) );

            //locate the storage space for the result
            void* storage =
                ( (bp::converter::rvalue_from_python_storage< std::pair<S,T> >*)data )->storage.bytes;

            //create the T container
            new (storage) std::pair<S,T>();

            std::pair<S,T> *container = static_cast<std::pair<S,T>*>(storage);

            //add all of the elements from the tuple - get the number of elements in the tuple
            container->first = bp::extract<S>(t[0])();
            container->second = bp::extract<T>(t[1])();

            data->convertible = storage;
        }
        else if (PyList_Check(obj_ptr))
        {
            //convert the PyObject to a boost::python::object
            bp::list l( bp::handle<>(bp::borrowed(obj_ptr)) );

            //locate the storage space for the result
            void* storage =
                ( (bp::converter::rvalue_from_python_storage< std::pair<S,T> >*)data )->storage.bytes;

            //create the T container
            new (storage) std::pair<S,T>();

            std::pair<S,T> *container = static_cast<std::pair<S,T>*>(storage);

            //add all of the elements from the tuple - get the number of elements in the tuple
            container->first = bp::extract<S>(l[0])();
            container->second = bp::extract<T>(l[1])();

            data->convertible = storage;
        }
    }
};

template<class S, class T>
struct to_py_pair
{
    static PyObject* convert(const std::pair<S,T> &cpp_pair)
    {
        auto raii = boost::python::release_gil_policy::acquire_gil();
        {
            bp::tuple python_tuple = bp::make_tuple( cpp_pair.first, cpp_pair.second );

            return bp::incref( python_tuple.ptr() );
        }
    }
};

template<class S, class T>
void register_pair()
{
    bp::to_python_converter< std::pair<S,T>, to_py_pair<S,T> >();

    bp::converter::registry::push_back( &from_py_pair<S,T>::convertible,
                                        &from_py_pair<S,T>::construct,
                                        bp::type_id< std::pair<S,T> >() );
}

SIRE_END_HEADER

#endif
