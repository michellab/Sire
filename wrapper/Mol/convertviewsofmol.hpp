/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2017  Christopher Woods
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

#ifndef PYWRAP_SIREMOL_CONVERTVIEWSOFMOL_HPP
#define PYWRAP_SIREMOL_CONVERTVIEWSOFMOL_HPP

#include "sireglobal.h"

#include <boost/python.hpp>
#include <boost/tuple/tuple.hpp>

#include "SireMol/viewsofmol.h"

#include "Helpers/release_gil_policy.hpp"

#include <QDebug>

namespace bp = boost::python;

SIRE_BEGIN_HEADER

/** This struct provides the from-Python conversion from a list of
    the views of a single molecule to a ViewsOfMol object */
struct viewsofmol_from_py_list
{
    /** Constructor - register the conversion functions
        for this type */
    viewsofmol_from_py_list()
    {
        boost::python::converter::registry::push_back(
            &convertible,
            &construct,
            bp::type_id< SireMol::ViewsOfMol >());
    }

    /** Test whether or not it is possible to convert the PyObject
        to a ViewsOfMol */
    static void* convertible(PyObject* obj_ptr)
    {
        auto raii = bp::release_gil_policy::acquire_gil();

        //is this a tuple type?
        if ( PyTuple_Check(obj_ptr) )
        {
            //check the tuple elements... - convert to a boost::tuple object
            bp::tuple t( bp::handle<>(bp::borrowed(obj_ptr)) );

            //how many elements are there?
            int n = PyTuple_Size(obj_ptr);

            //can they all be converted to type MoleculeView
            //and they all have the same molecule number
            SireMol::MolNum molnum;

            for (int i=0; i<n; ++i)
            {
                bp::extract<SireMol::MoleculeView*> get_molview(t[i]);

                if (get_molview.check())
                {
                    if (molnum.isNull())
                    {
                        molnum = get_molview()->data().number();
                    }
                    else if (molnum != get_molview()->data().number())
                    {
                        return 0;
                    }
                }
                else
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

            //can all of the elements be converted to type MoleculeView
            //and all have the same molecule number
            SireMol::MolNum molnum;

            for (int i=0; i<n; ++i)
            {
                bp::extract<SireMol::MoleculeView*> get_molview(l[i]);

                if (get_molview.check())
                {
                    if (molnum.isNull())
                    {
                        molnum = get_molview()->data().number();
                    }
                    else if (molnum != get_molview()->data().number())
                    {
                        return 0;
                    }
                }
                else
                    return 0;
            }

            //the list is ok!
            return obj_ptr;
        }
        else
            //could not recognise the type...
            return 0;
    }

    /** Construct a container of type ViewsOfMol from the PyObject pointed to
        by 'obj_ptr' */
    static void construct(
        PyObject* obj_ptr,
        bp::converter::rvalue_from_python_stage1_data* data)
    {
        auto raii = bp::release_gil_policy::acquire_gil();

        if (PyTuple_Check(obj_ptr))
        {
            //convert the PyObject to a boost::python::object
            bp::tuple t( bp::handle<>(bp::borrowed(obj_ptr)) );

            //locate the storage space for the result
            void* storage =
                ( (bp::converter::rvalue_from_python_storage<SireMol::ViewsOfMol>*)data
                             )->storage.bytes;

            //create the ViewsOfMol container
            new (storage) SireMol::ViewsOfMol();

            SireMol::ViewsOfMol *container = static_cast<SireMol::ViewsOfMol*>(storage);

            //add all of the elements from the tuple - get the number of elements in the tuple
            int n = PyTuple_Size(obj_ptr);

            if (n > 0)
            {
                *container = SireMol::ViewsOfMol(
                         *(bp::extract<SireMol::MoleculeView*>(t[0])()) );
            }

            for (int i=1; i<n; ++i)
            {
                container->add( bp::extract<SireMol::MoleculeView*>(t[i])()
                                     ->selection() );
            }

            data->convertible = storage;
        }
        else if (PyList_Check(obj_ptr))
        {
            //convert the PyObject to a boost::python::object
            bp::list l( bp::handle<>(bp::borrowed(obj_ptr)) );

            //locate the storage space for the result
            void* storage =
                ( (bp::converter::rvalue_from_python_storage<SireMol::ViewsOfMol>*)data
                       )->storage.bytes;

            //create the T container
            new (storage) SireMol::ViewsOfMol();

            SireMol::ViewsOfMol *container = static_cast<SireMol::ViewsOfMol*>(storage);

            //add all of the elements from the tuple - get the number of elements in the tuple
            int n = PyList_Size(obj_ptr);

            if (n > 0)
            {
                *container = SireMol::ViewsOfMol(
                         *(bp::extract<SireMol::MoleculeView*>(l[0])()) );
            }

            for (int i=1; i<n; ++i)
            {
                container->add( bp::extract<SireMol::MoleculeView*>(l[i])()
                                     ->selection() );
            }

            data->convertible = storage;
        }
    }
};

/** Function that returns the passed molecule view converted into
    a boost python object of the appropriate type, e.g. if the
    view is a single atom, it is an Atom, if it is a single residue,
    it is a Residue etc. etc. */
static bp::object get_molview(const SireMol::PartialMolecule &mol)
{
    if (mol.selectedAll())
    {
        return bp::object( mol.molecule() );
    }
    else if (mol.nAtoms() == 1)
    {
        return bp::object( mol.atom() );
    }

    if (mol.nResidues() == 1)
    {
        const auto res = mol.residue();
        if (mol.selection().selectedAll(res.number()))
        {
            return bp::object(res);
        }
    }

    if (mol.nCutGroups() == 1)
    {
        const auto cg = mol.cutGroup();
        if (mol.selection().selectedAll(cg.name()))
        {
            return bp::object(cg);
        }
    }

    if (mol.nChains() == 1)
    {
        const auto chain = mol.chain();
        if (mol.selection().selectedAll(chain.name()))
        {
            return bp::object(chain);
        }
    }

    if (mol.nSegments() == 1)
    {
        const auto segment = mol.segment();
        if (mol.selection().selectedAll(segment.name()))
        {
            return bp::object(segment);
        }
    }

    // this is a mixed view
    return bp::object(mol);
}

struct viewsofmol_to_py_list
{
    //this converts a ViewsOfMol into the most appropriate type,
    //depending on the views. If there is a single view, then
    //this will be an Atom, Residue, Molecule, PartialMolecule.
    //If this is multiple views, then this will be a list of
    //Atom, Residue, Molecule etc. views
    static PyObject* convert(const SireMol::ViewsOfMol &views)
    {
        auto raii = bp::release_gil_policy::acquire_gil();

        if (views.nViews() == 0)
        {
            //return None
            return 0;
        }
        else if (views.nViews() == 1)
        {
            //return the view itself
            bp::object obj = get_molview(views.valueAt(0));
            return bp::incref( obj.ptr() );
        }
        else
        {
            // are the views all of one type?
            QString typ;

            try
            {
                typ = views.getCommonType();
            }
            catch(...)
            {}

            if (typ == "SireMol::Atom")
            {
                return bp::incref(bp::object(views.atoms()).ptr());
            }
            else if (typ == "SireMol::Residue")
            {
                return bp::incref(bp::object(views.residues()).ptr());
            }
            else if (typ == "SireMol::Chain")
            {
                return bp::incref(bp::object(views.chains()).ptr());
            }
            else if (typ == "SireMol::CutGroup")
            {
                return bp::incref(bp::object(views.cutGroups()).ptr());
            }
            else if (typ == "SireMol::Segment")
            {
                return bp::incref(bp::object(views.segments()).ptr());
            }

            // this is a mixture of types, so return this as a
            // python list of views
            bp::list python_list;

            //add all items to the python list
            for (int i=0; i<views.nViews(); ++i)
            {
                python_list.append( get_molview(views.valueAt(i)) );
            }

            return bp::incref( python_list.ptr() );
        }
    }
};

void register_viewsofmol_list()
{
    bp::to_python_converter< SireMol::ViewsOfMol, viewsofmol_to_py_list >();

    bp::converter::registry::push_back( &viewsofmol_from_py_list::convertible,
                                        &viewsofmol_from_py_list::construct,
                                        bp::type_id<SireMol::ViewsOfMol>() );
}

SIRE_END_HEADER

#endif
