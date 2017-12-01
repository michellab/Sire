/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2007   Christopher Woods
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

#include <QVector>
#include <QSet>

#include <boost/tuple/tuple.hpp>

#include "Helpers/convertlist.hpp"
#include "Helpers/convertdict.hpp"
#include "Helpers/convertset.hpp"
#include "Helpers/tuples.hpp"
#include "Base/convertpackedarray.hpp"

#include "SireBase/propertymap.h"
#include "SireBase/unittest.h"
#include "SireBase/propertylist.h"
#include "SireBase/numberproperty.h"
#include "SireBase/stringproperty.h"

#include "SireCAS/expressionproperty.h"

using namespace SireBase;

using boost::python::register_tuple;

/** This struct provides the from-Python conversion from a list or
    tuple of Properties / strings / numbers to a PropertyList */
struct from_py_PropertyList
{
    /** Constructor - register the conversion functions
        for this type */
    from_py_PropertyList()
    {
        boost::python::converter::registry::push_back(
            &convertible,
            &construct,
            bp::type_id< PropertyList >());
    }

    template<class T>
    static bool can_convert_to_property(const T &val)
    {
        return bp::extract<double>(val).check() or
               bp::extract<QString>(val).check() or
               bp::extract<const Property*>(val).check() or
               bp::extract<const SireCAS::ExBase*>(val).check() or
               bp::extract<const SireCAS::Expression*>(val).check();
    }

    template<class T>
    static PropertyPtr convert_to_property(const T &val)
    {
        //number
        {
            bp::extract<double> e(val);
            if (e.check())
                return NumberProperty(e());
        }
        //string
        {
            bp::extract<QString> e(val);
            if (e.check())
                return StringProperty(e());
        }
        //Property class
        {
            bp::extract<const Property*> e(val) BOOST_EXTRACT_WORKAROUND;
            if (e.check())
                return PropertyPtr( *(e()) );
        }
        //SireCAS::ExBase class
        {
            bp::extract<const SireCAS::ExBase*> e(val) BOOST_EXTRACT_WORKAROUND;
            if (e.check())
                return SireCAS::ExpressionProperty( *(e()) );
        }
        //SireCAS::Expression class
        {
            bp::extract<const SireCAS::Expression*> e(val) BOOST_EXTRACT_WORKAROUND;
            if (e.check())
                return SireCAS::ExpressionProperty( *(e()) );
        }

        return PropertyPtr();
    }

    static bool is_list(PyObject* obj_ptr)
    {
        if (PyTuple_Check(obj_ptr))
        {
            bp::tuple t( bp::handle<>(bp::borrowed(obj_ptr)) );

            int n = PyTuple_Size(obj_ptr);

            for (int i=0; i<n; ++i)
            {
                if (not can_convert_to_property(t[i]))
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
                if (not can_convert_to_property(l[i]))
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
                if ( from_py_PropertyList::is_list( bp::object(t[i]).ptr()) )
                    return true;
            }

            return false;
        }
        else if (PyList_Check(obj_ptr))
        {
            bp::list l( bp::handle<>(bp::borrowed(obj_ptr)) );

            int n = PyList_Size(obj_ptr);

            for (int i=0; i<n; ++i)
            {
                if ( from_py_PropertyList::is_list( bp::object(l[i]).ptr()) )
                    return true;
            }

            return false;
        }
        else
            return false;
    }

    /** Test whether or not it is possible to convert the PyObject
        to a PropertyList */
    static void* convertible(PyObject* obj_ptr)
    {
        if (from_py_PropertyList::is_list(obj_ptr))
            return obj_ptr;

        else if (from_py_PropertyList::is_list_of_list(obj_ptr))
            return obj_ptr;

        else
            return 0;
    }

    static PropertyList convertToList(PyObject* obj_ptr)
    {
        //this is a single list
        if (PyTuple_Check(obj_ptr))
        {
            //this is a tuple
            bp::tuple t( bp::handle<>(bp::borrowed(obj_ptr)) );

            int n = PyTuple_Size(obj_ptr);

            PropertyList list;
 
            for (int i=0; i<n; ++i)
            {
                list.append( convert_to_property(t[i]) );
            }

            return list;
        }
        else if (PyList_Check(obj_ptr))
        {
            //this is a list
            bp::list l( bp::handle<>(bp::borrowed(obj_ptr)) );

            int n = PyList_Size(obj_ptr);

            PropertyList list;

            for (int i=0; i<n; ++i)
            {
                list.append( convert_to_property(l[i]) );
            }

            return list;
        }

        return PropertyList();
    }

    /** Construct a PropertyList from the PyObject pointed to
        by 'obj_ptr' */
    static void construct(
        PyObject* obj_ptr,
        bp::converter::rvalue_from_python_stage1_data* data)
    {
        PropertyList list;

        if (from_py_PropertyList::is_list_of_list(obj_ptr))
        {
            if (PyTuple_Check(obj_ptr))
            {
                bp::tuple t( bp::handle<>(bp::borrowed(obj_ptr)) );

                int n = PyTuple_Size(obj_ptr);

                for (int i=0; i<n; ++i)
                {
                    PropertyList sublist = from_py_PropertyList::convertToList(
                                                          bp::object(t[i]).ptr());

                    if (sublist.count() == 0)
                        list.append( convert_to_property(t[i]) );
                    else if (sublist.count() == 1)
                        list.append(sublist.at(0));
                    else
                        list.append(sublist);                   
                }
            }
            else if (PyList_Check(obj_ptr))
            {
                bp::list l( bp::handle<>(bp::borrowed(obj_ptr)) );
       
                int n = PyList_Size(obj_ptr);

                for (int i=0; i<n; ++i)
                {
                    PropertyList sublist = from_py_PropertyList::convertToList(
                                                        bp::object(l[i]).ptr());

                    if (sublist.count() == 0)
                        list.append( convert_to_property(l[i]) );
                    else if (sublist.count() == 1)
                        list.append(sublist.at(0));
                    else
                        list.append(sublist);
                }
            }
        }
        else if (from_py_PropertyList::is_list(obj_ptr))
        {
            list = from_py_PropertyList::convertToList(obj_ptr);
        }

        //locate the storage space for the result
        void* storage =
            ( (bp::converter::rvalue_from_python_storage<PropertyList>*)data )->storage.bytes;

        //create the PropertyList
        new (storage) PropertyList( list );

        data->convertible = storage;
    }
};

void register_PropertyList()
{
    //bp::to_python_converter< C, to_py_PackedArray<C> >();

    bp::converter::registry::push_back( &from_py_PropertyList::convertible,
                                        &from_py_PropertyList::construct,
                                        bp::type_id<PropertyList>() );
}

void register_SireBase_containers()
{
    register_list< QList< boost::shared_ptr<UnitTest> > >();
    register_list< QList< PropertyPtr > >();

    register_PropertyList();

    #if QT_VERSION >= 0x402000
    register_dict< QHash<QString,PropertyName> >();

    #else
    register_dict< QHash<QString,PropertyName>, QString, PropertyName>();

    #endif    
}
