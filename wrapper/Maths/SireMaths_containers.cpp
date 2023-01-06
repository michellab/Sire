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
#include "Helpers/pair.hpp"

#include "SireMaths/vector.h"
#include "SireMaths/matrix.h"
#include "SireMaths/nmatrix.h"
#include "SireMaths/nvector.h"
#include "SireMaths/accumulator.h"
#include "SireMaths/multifloat.h"
#include "SireMaths/multifixed.h"
#include "SireMaths/multidouble.h"
#include "SireMaths/multiuint.h"
#include "SireMaths/histogram.h"
#include "SireMaths/sphere.h"
#include "SireMaths/freeenergyaverage.h"

using namespace SireMaths;

using boost::python::register_tuple;

void register_SireMaths_containers()
{
    register_list< QList<Vector> >();
    register_list< QVector<Vector> >();
    register_list< QVector<Sphere> >();

    register_list< QVector<MultiFloat> >();
    register_list< QVector<MultiFixed> >();
    register_list< QVector<MultiDouble> >();
    register_list< QVector<MultiUInt> >();

    register_list< QVector<AverageAndStddev> >();
    register_list< QVector<FreeEnergyAverage> >();

    register_list< QVector<HistogramValue> >();

    register_pair< NVector,NMatrix >();

    register_tuple< boost::tuple<Vector,Vector,Vector> >();
    register_tuple< boost::tuple<Matrix,Matrix,Matrix> >();
}

