#ifndef CONVERTPROPERTY_HPP
#define CONVERTPROPERTY_HPP

#include <boost/python.hpp>
#include <boost/python/detail/indirect_traits.hpp>
#include <boost/mpl/if.hpp>
#include <boost/python/to_python_indirect.hpp>

#include "SireBase/property.h"

template<class P, class Base>
struct property_to_base_object
{
    static PyObject* convert(const P &property_container)
    {
        try
        {
            //try the first route, which will downcast to the derived type
            return boost::python::detail::make_owning_holder::execute<Base>
                   ( property_container.read().clone() );
        }
        catch(...)
        {
            PyErr_Clear();

            //failed as the derived type is not wrapped. Instead, 
            //return a python wrapper to the base type
            boost::python::object obj( property_container.read().clone() );
            return boost::python::incref( obj.ptr() );
        }
    }
};

template<class P, class Base>
void register_property_container()
{
    boost::python::to_python_converter< P, property_to_base_object<P,Base> >();
    boost::python::implicitly_convertible< Base, P >();
}

#endif

