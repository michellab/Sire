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

#include "SireSystem/errors.h"

using namespace SireSystem;

const char* missing_monitor::typeName()
{
    return QMetaType::typeName( qMetaTypeId<missing_monitor>() );
}

const char* duplicate_monitor::typeName()
{
    return QMetaType::typeName( qMetaTypeId<duplicate_monitor>() );
}

const char* missing_system::typeName()
{
    return QMetaType::typeName( qMetaTypeId<missing_system>() );
}

const char* duplicate_system::typeName()
{
    return QMetaType::typeName( qMetaTypeId<duplicate_system>() );
}

const char* constraint_error::typeName()
{
    return QMetaType::typeName( qMetaTypeId<constraint_error>() );
}

static const RegisterMetaType<missing_monitor> r_missmonitor;
static const RegisterMetaType<duplicate_monitor> r_dupmonitor;
static const RegisterMetaType<missing_system> r_misssystem;
static const RegisterMetaType<duplicate_system> r_dupsystem;
static const RegisterMetaType<constraint_error> r_constraint;
