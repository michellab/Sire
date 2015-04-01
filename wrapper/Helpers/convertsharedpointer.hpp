/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2007  Christopher Woods
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

#ifndef PYWRAP_SIREPY_CONVERTSHAREDPOINTER_HPP
#define PYWRAP_SIREPY_CONVERTSHAREDPOINTER_HPP

namespace bp = boost::python;

#include "SireBase/sharedpolypointer.hpp"

SIRE_BEGIN_HEADER

using SireBase::SharedPolyPointer;

template<class P, class Base>
struct to_base_object
{
    static PyObject* convert(const P &object_container)
    {
        return bp::incref( object( SireBase::SharedPolyPointer<Base>(object_container.base()) ).ptr() );
    }
};

template<class P, class Base>
void register_container()
{
    bp::to_python_converter< P, to_base_object<P,Base> >();

    bp::implicitly_convertible< Base, P >();

    bp::register_ptr_to_python< SharedPolyPointer<Base> >();
}

SIRE_END_HEADER

#endif
