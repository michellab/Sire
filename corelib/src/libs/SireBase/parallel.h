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

#ifndef SIREBASE_PARALLEL_H
#define SIREBASE_PARALLEL_H

#include "sireglobal.h"

#include "SireBase/propertymap.h"
#include "SireBase/booleanproperty.h"

SIRE_BEGIN_HEADER

#ifndef GCCXML_PARSE

#include <QVector>
#include <QMutex>

// We have to undef the 'emit' from Qt as this is a function
// name used in TBB! This should be safe for Qt as that
// code uses Q_EMIT, and we don't use signals and slots in Sire
#ifdef emit
#undef emit
#endif

#include <tbb/parallel_for.h>
#include <tbb/parallel_for_each.h>
#include <tbb/parallel_reduce.h>
#include <tbb/parallel_invoke.h>

#include <memory>

namespace SireBase
{
    inline bool should_run_in_parallel(int count, const PropertyMap &map=PropertyMap())
    {
        if (count < 8)
            return false;
        else if (map["parallel"].hasValue())
            return map["parallel"].value().asA<BooleanProperty>().value();
        else
            return true;
    }


    /** This function runs the passed array T of functions in parallel, if
        the optional 'run_parallel' is true. Otherwise, it runs the functions
        serially, one after another */
    template<class T>
    void parallel_invoke( const T &functions, bool run_parallel=true )
    {
        if (run_parallel)
        {
            tbb::parallel_for( tbb::blocked_range<int>(0, functions.count()),
                               [&](const tbb::blocked_range<int> &r)
            {
                for (int i=r.begin(); i<r.end(); ++i)
                {
                    functions[i]();
                }
            });
        }
        else
        {
            for (int i=0; i<functions.count(); ++i)
            {
                functions[i]();
            }
        }
    }

} // end of namespace SireBase

#endif

SIRE_END_HEADER

#endif

