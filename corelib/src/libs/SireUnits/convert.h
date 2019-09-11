/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2006  Christopher Woods
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

#ifndef SIREUNITS_CONVERT_H
#define SIREUNITS_CONVERT_H

#include "SireUnits/units.h"
#include "SireUnits/temperature.h"

#include <QDebug>

SIRE_BEGIN_HEADER

namespace SireUnits
{

/////////////////////////////////////////////////
// Conversion functions for derived units      //
/////////////////////////////////////////////////

SIRE_ALWAYS_INLINE double convertFrom(double val, const Dimension::TempBase &from_units)
{
    return from_units.convertToInternal(val);
}

SIRE_ALWAYS_INLINE double convertFrom(double val, const Dimension::Unit &from_units)
{
    return from_units.convertToInternal(val);
}

SIRE_ALWAYS_INLINE double convertTo(double val, const Dimension::TempBase &to_units)
{
    return to_units.convertFromInternal(val);
}

SIRE_ALWAYS_INLINE double convertTo(double val, const Dimension::Unit &to_units)
{
    return to_units.convertFromInternal(val);
}

SIRE_ALWAYS_INLINE double convert(double val, const Dimension::TempBase &from_units,
                                  const Dimension::TempBase &to_units)
{
    return convertTo( convertFrom(val,from_units), to_units );
}

SIRE_ALWAYS_INLINE double convert(double val, const Dimension::Unit &from_units,
                                  const Dimension::TempBase &to_units)
{
    return convertTo( convertFrom(val,from_units), to_units );
}

SIRE_ALWAYS_INLINE double convert(double val, const Dimension::Unit &from_units,
                                  const Dimension::Unit &to_units)
{
    return convertTo( convertFrom(val,from_units), to_units );
}

SIRE_ALWAYS_INLINE double convert(double val, const Dimension::TempBase &to_units)
{
    return convertTo(val, to_units);
}

SIRE_ALWAYS_INLINE double convert(double val, const Dimension::Unit &to_units)
{
    return convertTo(val, to_units);
}

}

SIRE_EXPOSE_FUNCTION( SireUnits::convert )
SIRE_EXPOSE_FUNCTION( SireUnits::convertTo )
SIRE_EXPOSE_FUNCTION( SireUnits::convertFrom )

SIRE_END_HEADER

#endif
