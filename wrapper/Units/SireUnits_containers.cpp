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

#include "SireUnits/dimensions.h"

#include "SireBase/packedarray2d.hpp"

using namespace SireUnits;

using boost::python::register_tuple;

void register_SireUnits_containers()
{
    register_list< QVector<Dimension::Charge> >();
    register_list< QVector<Dimension::Mass> >();
    register_list< QVector<Dimension::MolarMass> >();
    register_list< QVector<Dimension::Length> > ();
    register_list< QVector<Dimension::MolarEnergy> > ();

    register_list< QList<Dimension::Charge> >();
    register_list< QList<Dimension::Mass> >();
    register_list< QList<Dimension::MolarMass> >();
    register_list< QList<Dimension::Length> > ();
    register_list< QList<Dimension::MolarEnergy> > ();

    register_PackedArray< SireBase::PackedArray2D<Dimension::Charge> >();
    register_PackedArray< SireBase::PackedArray2D<Dimension::Mass> >();
    register_PackedArray< SireBase::PackedArray2D<Dimension::MolarMass> >();
}
