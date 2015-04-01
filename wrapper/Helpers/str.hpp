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

#ifndef PYWRAP_SIREPY_STR_HPP
#define PYWRAP_SIREPY_STR_HPP

#include <Python.h>
#include <boost/python.hpp>
#include <boost/mpl/if.hpp>
#include <boost/type_traits/is_pod.hpp>

#include <QString>

#include "sireglobal.h"

using namespace boost::python;

SIRE_BEGIN_HEADER

/** The function used to get a string representation of a POD type */
template<class T>
struct __str__pod
{
    static object __str__(const T &object)
    {
        return str(object);
    }
};

/** Template function used to get a string representation of 'obj' */
template<class T>
QString get_string_representation(const T &obj)
{
    return obj.toString();
}

/** The function used to get a string representation of a non-POD type - 
    the non-POD type must implement a QString conversion function */
template<class T>
struct __str__nonpod
{
    static object __str__(const T &object)
    {
        return str( get_string_representation<T>(object) );
    }
};

/** The function used to get a string representation of an object for 
    a non-pod type - it returns the string representation in quotes
    unless the object is a container type */
template<class T>
struct __str__nonpod_list
{
    static object __str__(const T &obj)
    {
        //convert 'obj' to an object
        object pyobj( obj );
    
        //is this a container object?
        //how do I do this?
      
        //add quotes to the string
        return "'" + str(pyobj) + "'";
    }
};

/** The function used to return a string representation 
    of an object. */
template<class T>
object __str__(const T &object)
{
    return boost::mpl::if_< boost::is_pod<T>,
                            __str__pod<T>,
                            __str__nonpod<T> >::type::__str__( object );
}

/** Return the string form of the object 'obj' */
template<class T>
object get_list_string( const T &obj )
{
    return boost::mpl::if_< boost::is_pod<T>,
                            __str__pod<T>,
                            __str__nonpod_list<T> >::type::__str__( obj );
}

/** The function used to return a string representation of
    a list of objects - the contents of the list must all have been 
    exposed to python for this to work! */
template<class C>
object __str__list(const C &container)
{
    if (container.isEmpty())
        return str("[]");
    else
    {
        //create a list to hold the string representations of the 
        //contents of this list
        list strlist;
    
        //append a string version of the object to a python list...
        for (typename C::const_iterator it = container.begin();
             it != container.end();
             ++it)
        {
            strlist.append( get_list_string( *it ) );
        }
    
        //now join the items together and bracket the resulting string
        return str("[") + str(", ").join(strlist) + str("]");
    }
}

/** The function used to return a string representation of 
    a hash or map of objects - the contents and keys of the 
    hash or map must have been exposed to python for this to work! */
template<class C>
object __str__dict(const C &container)
{
    if (container.isEmpty())
        return str("{}");
    else
    {
        //create a list to hold the string representation of the 
        //key-value pairs
        list strlist;
        
        //append a string version of the key-value pair to a python list...
        for (typename C::const_iterator it = container.begin();
             it != container.end();
             ++it)
        {
            strlist.append( get_list_string(it.key()) + ": " + 
                            get_list_string(it.value()) );
        }
    
        //now join the items together and bracket the resulting string
        return str("{") + str(", ").join(strlist) + str("}");
    }
}

SIRE_END_HEADER

#endif
