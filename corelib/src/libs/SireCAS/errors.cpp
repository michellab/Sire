/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2007  Christopher Woods
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

#include "SireCAS/errors.h"

using namespace SireCAS;

const char* unregistered_expression::typeName()
{
    return QMetaType::typeName( qMetaTypeId<unregistered_expression>() );
}

const char* unavailable_differential::typeName()
{
    return QMetaType::typeName( qMetaTypeId<unavailable_differential>() );
}

const char* unavailable_integral::typeName()
{
    return QMetaType::typeName( qMetaTypeId<unavailable_integral>() );
}

const char* rearrangement_error::typeName()
{
    return QMetaType::typeName( qMetaTypeId<rearrangement_error>() );
}

const char* invalid_symbol::typeName()
{
    return QMetaType::typeName( qMetaTypeId<invalid_symbol>() );
}

const char* missing_symbol::typeName()
{
    return QMetaType::typeName( qMetaTypeId<missing_symbol>() );
}

static const
RegisterMetaType<unregistered_expression> r_unreg;

static const
RegisterMetaType<unavailable_differential> r_diff;

static const
RegisterMetaType<unavailable_integral> r_int;

static const
RegisterMetaType<rearrangement_error> r_rearrange;

static const
RegisterMetaType<invalid_symbol> r_invalid_symbol;

static const
RegisterMetaType<missing_symbol> r_missing_symbol;
