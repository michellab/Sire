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

#include "SireFF/ffname.h"
#include "SireFF/forcefield.h"
#include "SireFF/forcefields.h"
#include "SireFF/point.h"

#include "SireBase/properties.h"

using namespace SireFF;
using namespace SireBase;

using boost::python::register_tuple;

void register_SireFF_containers()
{
    register_list< QVector<ForceField> >();
    register_list< QVector<ForceFields> >();
    register_list< QVector<PointPtr> >();
    register_list< QList<PointPtr> >();

    #if QT_VERSION >= QT_VERSION_CHECK(4, 2, 0)
    register_dict< QHash<FFName,Properties> >();

    #else
    register_dict< QHash<FFName,Properties>, FFName, Properties > >();

    #endif
}
