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

#include "refcountdata.h"

#include "SireError/getbacktrace.h"

#include <QDebug>

using namespace SireBase;

Q_GLOBAL_STATIC( tbb::spin_mutex, get_spin_mutex );

namespace SireBase
{
    namespace detail
    {
        tbb::spin_mutex SIREBASE_EXPORT *get_shared_null_mutex()
        {
            return get_spin_mutex();
        }
    }
}

/** Constructor */
RefCountData::RefCountData()
{}

/** Copying is forbidden */
RefCountData::RefCountData(const RefCountData&o)
{}

/** Destructor */
RefCountData::~RefCountData()
{}

/** Copy assignment is forbidden */
RefCountData& RefCountData::operator=(RefCountData&)
{
    return *this;
}

/** Only the same if exactly the same object */
bool RefCountData::operator==(const RefCountData &other) const
{
    return this == &other;
}

/** Only the same if exactly the same object */
bool RefCountData::operator!=(const RefCountData &other) const
{
    return not operator==(other);
}

/** Constructor */
RefCountData::Counter::Counter()
{
    refcount = 0;
}

/** Destructor */
RefCountData::Counter::~Counter()
{}

/** Error function called when a counter is dereferenced below 0 */
void RefCountData::Counter::doubleDereferenced() const
{
    qDebug() << "WARNING - MEMORY CORRUPTION: PROGRAM BUG: REFCOUNT" << qintptr(this)
             << "HAS BEEN DEREFERENCED BELOW ZERO - DOUBLE FREE'D. BACKTRACE\n"
             << SireError::getBackTrace().join("\n");
}
