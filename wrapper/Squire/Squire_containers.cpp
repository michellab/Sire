/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) <year>  <name of author>
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

#include "SireBase/property.h"

#include "Squire/qmprogram.h"
#include "Squire/pointcharge.h"
#include "Squire/gto.h"
#include "Squire/sgto.h"
#include "Squire/pgto.h"

using boost::python::register_tuple;

using namespace Squire;

void register_Squire_containers()
{
    register_tuple< boost::tuple<double,GTOPtr> >();
    register_tuple< boost::tuple<GTOPtr,double> >();

    register_list< QVector< boost::tuple<double,GTOPtr> > >();
    register_list< QVector< boost::tuple<GTOPtr,double> > >();

    register_list< QVector<PointCharge> >();

    register_list< QVector<S_GTO> >();
    register_list< QVector<P_GTO> >();

    register_list< QVector<GTOPtr> >();
}
