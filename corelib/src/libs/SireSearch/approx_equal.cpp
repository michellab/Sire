/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2022  Christopher Woods
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

#include "approx_equal.h"

#include <QReadWriteLock>

#include <QDebug>

struct Epsilon
{
    double epsilon = 1e-6;
    QReadWriteLock lock;
};

Q_GLOBAL_STATIC(Epsilon, epsilon);

namespace SireSearch
{
    double SIRESEARCH_EXPORT get_approx_epsilon()
    {
        auto e = epsilon();
        QReadLocker lkr(&(e->lock));
        return e->epsilon;
    }

    void SIRESEARCH_EXPORT set_approx_epsilon(double eps)
    {
        auto e = epsilon();
        QWriteLocker lkr(&(e->lock));
        e->epsilon = abs(eps);
    }

    bool SIRESEARCH_EXPORT approx_equal(double val0, double val1)
    {
        double delta = abs(val1 - val0);
        double eps = get_approx_epsilon();

        //check absolute difference is less than epsilon**2
        if (delta <= (eps*eps))
            return true;

        //check relative difference is less than eps
        return delta <= 0.5 * eps * (abs(val0) + abs(val1));
    }
}
