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

#include <QVector>
#include <QMutex>

#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/parallel_invoke.h>
#include <tbb/tbb_exception.h>

#include <memory>

namespace SireBase
{
    namespace parallel
    {
        namespace detail
        {
            template<class ARG>
            size_t get_min_container_size(const std::vector<ARG> &arg)
            {
                return arg.size();
            }

            template<class ARG1, class ARG2>
            size_t get_min_container_size(const std::vector<ARG1> &arg1,
                                          const std::vector<ARG2> &arg2)
            {
                return std::min(arg1.size(), arg2.size());
            }

            template<class ARG1, class ARG2, class... ARGS>
            size_t get_min_container_size(const std::vector<ARG1> &arg1,
                                          const std::vector<ARG2> &arg2,
                                          const std::vector<ARGS>&... args)
            {
                size_t minsize = get_min_container_size(args...);
                return std::min( minsize, get_min_container_size(arg1,arg2) );
            }
        }

        /** This will map the function 'func' against the array(s) of argument(s) in args,
            returning a vector of results */
        template<class FUNC, class... ARGS>
        auto map(FUNC func, const QVector<ARGS>&... args)
        {
            typedef typename std::result_of<FUNC(ARGS...)>::type RETURN_TYPE;

            int nvals=detail::get_min_container_size(args...);

            QVector<RETURN_TYPE> result(nvals);

            RETURN_TYPE *result_data = result.data();

            tbb::parallel_for( tbb::blocked_range<int>(0,nvals),
                               [&](tbb::blocked_range<int> r)
            {
                for (int i=r.begin(); i<r.end(); ++i)
                {
                     result_data[i] = mapfunc(args[i]...);
                }
            });

            return result;
        }

        /** This will reduce the passed array of values using the function 'func',
            starting from the initial value 'initial' */
        template<class FUNC, class T>
        T reduce(FUNC func, const QVector<T> &values, const T &initial)
        {
            if (values.empty())
            {
                return initial;
            }
            else
            {
                const T *data = values.constData();
            
                const int nvals = values.count();
            
                T result =
                    tbb::parallel_reduce( tbb::blocked_range<int>(0,nvals),
                                          initial,
                        [&](tbb::blocked_range<int> r, T task_result)
                    {
                        for (int i=r.begin(); i<r.end(); ++i)
                        {
                            task_result = func(task_result,data[i]);
                        }

                        return task_result;

                    }, func );

                return result;
            }
        }

        /** This will reduce the passed array of values using the function 'func' */
        template<class FUNC, class T>
        T reduce(FUNC func, const QVector<T> &values)
        {
            return reduce(func, values, T(0));
        }

        template<class MAPFUNC, class REDFUNC, class... ARGS>
        auto mapReduce(MAPFUNC mapfunc, REDFUNC redfunc, const QVector<ARGS>&... args)
        {
            typedef typename std::result_of<MAPFUNC(ARGS...)>::type RETURN_TYPE;

            int nvals=detail::get_min_container_size(args...);

            RETURN_TYPE result =
                 tbb::parallel_reduce( tbb::blocked_range<int>(0,nvals),
                                       RETURN_TYPE(0),
                       [&](tbb::blocked_range<int> r, RETURN_TYPE task_result)
                 {
                     for (int i=r.begin(); i<r.end(); ++i)
                     {
                         task_result = redfunc(task_result, mapfunc(args[i]...) );
                     }

                     return task_result;

                 }, redfunc );

            return result;                               
        }
    }
}

#endif

