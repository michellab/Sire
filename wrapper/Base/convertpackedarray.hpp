/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2008  Christopher Woods
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

#ifndef PYWRAP_SIREPY_CONVERTPACKEDARRAY_HPP
#define PYWRAP_SIREPY_CONVERTPACKEDARRAY_HPP

#include "sireglobal.h"

#include <boost/python.hpp>
#include <boost/tuple/tuple.hpp>

namespace bp = boost::python;

SIRE_BEGIN_HEADER

/** This struct provides the from-Python conversion from a list or
    tuple to a PackedArray of type 'C' */
template<class C>
struct from_py_PackedArray
{
    typedef typename C::value_type T;

    /** Constructor - register the conversion functions
        for this type */
    from_py_PackedArray()
    {
        boost::python::converter::registry::push_back(
            &convertible,
            &construct,
            bp::type_id< C >());
    }

    static bool is_list(PyObject* obj_ptr)
    {
        if (PyTuple_Check(obj_ptr))
        {
            bp::tuple t( bp::handle<>(bp::borrowed(obj_ptr)) );

            int n = PyTuple_Size(obj_ptr);

            for (int i=0; i<n; ++i)
            {
                if (not bp::extract<T>(t[i]).check())
                    return false; 
            }

            return true;
        }
        else if (PyList_Check(obj_ptr))
        {
            bp::list l( bp::handle<>(bp::borrowed(obj_ptr)) );

            int n = PyList_Size(obj_ptr);
     
            for (int i=0; i<n; ++i)
            {
                if (not bp::extract<T>(l[i]).check())
                    return false;
            }

            return true;
        }
        else
            return false;
    }

    static bool is_list_of_list(PyObject* obj_ptr)
    {
        if (PyTuple_Check(obj_ptr))
        {
            bp::tuple t( bp::handle<>(bp::borrowed(obj_ptr)) );

            int n = PyTuple_Size(obj_ptr);

            for (int i=0; i<n; ++i)
            {
                if (not from_py_PackedArray<C>::is_list( bp::object(t[i]).ptr()) )
                    return false;
            }

            return true;
        }
        else if (PyList_Check(obj_ptr))
        {
            bp::list l( bp::handle<>(bp::borrowed(obj_ptr)) );

            int n = PyList_Size(obj_ptr);

            for (int i=0; i<n; ++i)
            {
                if (not from_py_PackedArray<C>::is_list( bp::object(l[i]).ptr()) )
                    return false;
            }

            return true;
        }
        else
            return false;
    }

    /** Test whether or not it is possible to convert the PyObject
        to a PackedArray where all of the elements are of type 'T' */
    static void* convertible(PyObject* obj_ptr)
    {
        if (from_py_PackedArray<C>::is_list(obj_ptr))
            return obj_ptr;

        else if (from_py_PackedArray<C>::is_list_of_list(obj_ptr))
            return obj_ptr;

        else
            return 0;
    }

    static QVector<T> convertToVector(PyObject* obj_ptr)
    {
        //this is a single list
        if (PyTuple_Check(obj_ptr))
        {
            //this is a tuple
            bp::tuple t( bp::handle<>(bp::borrowed(obj_ptr)) );

            int n = PyTuple_Size(obj_ptr);

            QVector<T> array(n);
 
            for (int i=0; i<n; ++i)
            {
                array[i] = bp::extract<T>(t[i])();
            }

            return array;
        }
        else if (PyList_Check(obj_ptr))
        {
            //this is a list
            bp::list l( bp::handle<>(bp::borrowed(obj_ptr)) );

            int n = PyList_Size(obj_ptr);

            QVector<T> array(n);

            for (int i=0; i<n; ++i)
            {
                array[i] = bp::extract<T>(l[i])();
            }

            return array;
        }

        return QVector<T>();
    }

    /** Construct a container of type T from the PyObject pointed to
        by 'obj_ptr' */
    static void construct(
        PyObject* obj_ptr,
        bp::converter::rvalue_from_python_stage1_data* data)
    {
        QVector< QVector<T> > packed_array;

        if (from_py_PackedArray<C>::is_list(obj_ptr))
        {
            packed_array.append( from_py_PackedArray<C>::convertToVector(obj_ptr) );
        }
        else if (from_py_PackedArray<C>::is_list_of_list(obj_ptr))
        {
            if (PyTuple_Check(obj_ptr))
            {
                bp::tuple t( bp::handle<>(bp::borrowed(obj_ptr)) );

                int n = PyTuple_Size(obj_ptr);

                for (int i=0; i<n; ++i)
                {
                    packed_array.append( from_py_PackedArray<C>::convertToVector(
                                                                  bp::object(t[i]).ptr()) );
                }
            }
            else if (PyList_Check(obj_ptr))
            {
                bp::list l( bp::handle<>(bp::borrowed(obj_ptr)) );
       
                int n = PyList_Size(obj_ptr);
    
                for (int i=0; i<n; ++i)
                {
                    packed_array.append( from_py_PackedArray<C>::convertToVector(
                                                                  bp::object(l[i]).ptr()) );
                }
            }
        }

        //locate the storage space for the result
        void* storage =
            ( (bp::converter::rvalue_from_python_storage<C>*)data )->storage.bytes;

        //create the T container
        new (storage) C( packed_array );

        data->convertible = storage;
    }
};

template<class C>
struct to_py_PackedArray
{
    typedef typename C::value_type T;

    static PyObject* convert(const C &packed_array)
    {
        bp::list python_list;

        if (packed_array.nArrays() == 1)
        {
            //just a single array
            const T *data = packed_array.constData(0);
            int nvalues = packed_array.nValues(0);

            for (int i=0; i<nvalues; ++i)
            {
                python_list.append( data[i] );
            }
        }
        else
        {
            for (int i=0; i<packed_array.nArrays(); ++i)
            {
                bp::list array;

                const T *data = packed_array.constData(i);
                int nvalues = packed_array.nValues(i);

                for (int j=0; j<nvalues; ++j)
                {
                    array.append( data[j] );
                }

                python_list.append( array );
            }
        }

        return bp::incref( python_list.ptr() );
    }
};

template<class C>
void register_PackedArray()
{
    bp::to_python_converter< C, to_py_PackedArray<C> >();

    bp::converter::registry::push_back( &from_py_PackedArray<C>::convertible,
                                        &from_py_PackedArray<C>::construct,
                                        bp::type_id<C>() );
}

SIRE_END_HEADER

#endif
