#ifndef CONVERTPROPERTY_HPP
#define CONVERTPROPERTY_HPP

#include "SireBase/sharedpolypointer.hpp"
#include "SireBase/shareddatapointer.hpp"
#include "SireBase/property.h"

namespace SireBase
{
    template<class T> 
    const T * get_pointer(const SharedPolyPointer<T> & p)
    {
        return p.constData();
    }

    template<class T> T * get_pointer(SharedPolyPointer<T> & p )
    {
        return p.data();
    }

    template<class T>
    const T * get_pointer(const SharedDataPointer<T> & p)
    {
        return p.constData();
    }

    template<class T> T * get_pointer(SharedDataPointer<T> & p )     
    {
        return p.data();
    }
}

template<class P, class Base>
struct property_to_base_object
{
    static PyObject* convert(const P &property_container)
    {
        return boost::python::incref( boost::python::object( 
                         SireBase::SharedPolyPointer<Base>(property_container.read()) ).ptr() );
    }
};

template<class P, class Base>
void register_property_container()
{
    boost::python::to_python_converter< P, property_to_base_object<P,Base> >();

    boost::python::implicitly_convertible< Base, P >();

    boost::python::register_ptr_to_python< SireBase::SharedPolyPointer<Base> >();
}

#endif

