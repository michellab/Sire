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

#ifndef PYWRAP_SIREPY_CONVERTDICT_HPP
#define PYWRAP_SIREPY_CONVERTDICT_HPP

#include "sireglobal.h"

#include <boost/python.hpp>
#include <boost/tuple/tuple.hpp>

#include "Helpers/release_gil_policy.hpp"

namespace bp = boost::python;

SIRE_BEGIN_HEADER

/** This struct provides the from-Python conversion from a dict
    to a dict or hash-like container of type 'C' (e.g. QHash, QMap) */
template<class C, class key_type, class mapped_type>
struct from_py_dict
{
    /** Constructor - register the conversion functions
        for this type */
    from_py_dict()
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
        auto raii = boost::python::release_gil_policy::acquire_gil();

        //is this a dict type?
        if ( PyDict_Check(obj_ptr) )
        {
            //check the tuple elements... - convert to a boost::tuple object
            bp::dict d( bp::handle<>(bp::borrowed(obj_ptr)) );

            //check the items in the dict (items is a list of 2-tuples key-value
            bp::list items = d.items();

            int nitems = bp::extract<int>(items.attr("__len__")())();

            for (int i=0; i<nitems; ++i)
            {
                bp::tuple item = bp::extract<bp::tuple>(items[i])();

                if ( not (bp::extract<key_type>(item[0]).check() and
                        bp::extract<mapped_type>(item[1]).check()) )
                {
                    //either the key of value is wrong
                    return 0;
                }
            }

            //the tuple is ok!
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
        {
            //convert the PyObject to a boost::python::dict
            bp::dict d( bp::handle<>(bp::borrowed(obj_ptr)) );

            //locate the storage space for the result
            void* storage =
                ( (bp::converter::rvalue_from_python_storage<C>*)data )->storage.bytes;

            //create the container
            new (storage) C();

            C *container = static_cast<C*>(storage);

            //add all of the elements from the dict - do this by converting
            //to a list and then extracting each item
            bp::list items = d.items();

            int nitems = bp::extract<int>(items.attr("__len__")())();

            for (int i=0; i<nitems; ++i)
            {
                bp::tuple item = bp::extract<bp::tuple>(items[i])();

                container->insert( bp::extract<key_type>(item[0])(),
                                bp::extract<mapped_type>(item[1])() );
            }

            data->convertible = storage;
        }
    }
};

template<class C>
struct to_py_dict
{
    static PyObject* convert(const C &cpp_dict)
    {
        auto raii = boost::python::release_gil_policy::acquire_gil();
        {
            bp::dict python_dict;

            //add all items to the python dictionary
            for (typename C::const_iterator it = cpp_dict.begin();
                it != cpp_dict.end();
                ++it)
            {
                python_dict[it.key()] = it.value();
            }

            return bp::incref( python_dict.ptr() );
        }
    }
};

template<class C>
void register_dict()
{
    typedef typename C::key_type key_type;
    typedef typename C::mapped_type mapped_type;

    bp::to_python_converter< C, to_py_dict<C> >();

    bp::converter::registry::push_back( &from_py_dict<C,key_type,mapped_type>::convertible,
                                        &from_py_dict<C,key_type,mapped_type>::construct,
                                        bp::type_id<C>() );
}

template<class C, class key_type, class mapped_type>
void register_dict()
{
    bp::to_python_converter< C, to_py_dict<C> >();

    bp::converter::registry::push_back( &from_py_dict<C,key_type,mapped_type>::convertible,
                                        &from_py_dict<C,key_type,mapped_type>::construct,
                                        bp::type_id<C>() );
}

SIRE_END_HEADER

#endif
