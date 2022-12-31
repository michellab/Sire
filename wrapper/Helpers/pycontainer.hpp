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

#ifndef SIRE_PYWRAP_PYCONTAINER_HPP
#define SIRE_PYWRAP_PYCONTAINER_HPP

#include "Python.h"
#include <boost/python.hpp>
#include <boost/python/slice.hpp>

#include <boost/type_traits/is_pod.hpp>
#include <boost/mpl/if.hpp>

#include <exception>

#include <QObject>
#include <QDebug>

#include "sireglobal.h"

using namespace boost::python;

/** Struct that returns the return policy for the reference of type 'V'.
    This returns 'copy_non_const_reference' for POD types, and
    return_internal_reference< 1, with_custodian_and_ward<1,2> > for non-POD
    types */
template<class V>
struct get_return_reference_policy
{
    typedef typename boost::mpl::if_< boost::is_pod<V>,
                             return_value_policy<copy_non_const_reference>,
                             return_internal_reference< 1,
                                with_custodian_and_ward<1,2> > >::type type;
};

template<class C, class VALUE>
VALUE& __getitem__list(C &container, int index)
{
    int container_size = container.size();

    //implement negative indicies
    if (index < 0) index += container_size;

    if ( index < 0 or index >= container_size )
    {
        PyErr_SetString( PyExc_IndexError,
                qPrintable(QObject::tr(
                    "Index out of range (%1 from %2)")
                        .arg(index).arg(container_size)) );

        qDebug() << CODELOC;
        PyErr_Print();

        throw_error_already_set();
    }

    return container[index];
}

template<class C, class VALUE>
VALUE __getitem__list_const(const C &container, int index)
{
    int container_size = container.size();

    //implement negative indicies
    if (index < 0) index += container_size;

    if ( index < 0 or index >= container_size )
    {
        PyErr_SetString( PyExc_IndexError,
                qPrintable(QObject::tr(
                    "Index out of range (%1 from %2)")
                        .arg(index).arg(container_size)) );

        qDebug() << CODELOC;
        PyErr_Print();

        throw_error_already_set();
    }

    return container[index];
}

template<class C>
C __getitem__slice(const C &container, const slice &index)
{
    slice::range<typename C::const_iterator> bounds;

    try
    {
        bounds = index.get_indicies<>(container.begin(), container.end());
    }
    catch (std::invalid_argument)
    {
        //invalid slice - return an empty container
        return C();
    }

    //extract all of the elements
    C returned_container;

    while (bounds.start != bounds.stop)
    {
        returned_container.append( *bounds.start );
        std::advance( bounds.start, bounds.step);
    }
    returned_container.append( *bounds.start );

    return returned_container;
}

template<class C, class KEY, class VALUE>
VALUE& __getitem__dict(C &container, const KEY &key)
{
    return container[key];
}

template<class C, class KEY, class VALUE>
VALUE __getitem__dict_const(const C &container, const KEY &key)
{
    if (not container.contains(key))
    {
        PyErr_SetString( PyExc_KeyError,
                qPrintable(QObject::tr(
                    "Key not in container.")));

        qDebug() << CODELOC;
        PyErr_Print();

        boost::python::throw_error_already_set();
    }

    return container[key];
}

template<class C, class VALUE>
void __setitem__list(C &container, int index, const VALUE &value)
{
    int container_size = container.size();

    //implement negative indicies
    if (index < 0) index += container_size;

    if ( index < 0 or index >= container_size )
    {
        PyErr_SetString( PyExc_IndexError,
                qPrintable(QObject::tr(
                    "Index out of range (%1 from %2)")
                        .arg(index).arg(container_size)) );

        qDebug() << CODELOC;
        PyErr_Print();

        boost::python::throw_error_already_set();
    }

    container[index] = value;
}

template<class C, class KEY, class VALUE>
void __setitem__dict(C &container, const KEY &key, const VALUE &value)
{
    container[key] = value;
}

template<class C>
void __delitem__list(C &container, int index)
{
    int container_size = container.size();

    //implement negative indicies
    if (index < 0) index += container_size;

    if ( index < 0 or index >= container_size )
    {
        PyErr_SetString( PyExc_IndexError,
                qPrintable(QObject::tr(
                    "Index out of range (%1 from %2)")
                        .arg(index).arg(container_size)) );

        qDebug() << CODELOC;
        PyErr_Print();

        boost::python::throw_error_already_set();
    }

    container.removeAt(index);
}

template<class C>
void __delitem__list_remove(C &container, int index)
{
    int container_size = container.size();

    //implement negative indicies
    if (index < 0) index += container_size;

    if ( index < 0 or index >= container_size )
    {
        PyErr_SetString( PyExc_IndexError,
                qPrintable(QObject::tr(
                    "Index out of range (%1 from %2)")
                        .arg(index).arg(container_size)) );

        qDebug() << CODELOC;
        PyErr_Print();

        boost::python::throw_error_already_set();
    }

    container.remove(index);
}

template<class C>
int __len__(const C &container)
{
    return container.size();
}

template<class C, class KEY>
void __delitem__dict(C &container, const KEY &key)
{
    container.remove(key);
}

template<class C, class KEY>
bool __contains__(const C &container, const KEY &key)
{
    return container.contains(key);
}

template<class C, class VALUE>
void __lshift__(C &container, const VALUE &value)
{
    container << value;
}

template<class C, class VALUE>
void __rshift__(C &container, VALUE &value)
{
    container >> value;
}

#endif
