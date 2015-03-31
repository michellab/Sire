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

#include "SireFF/errors.h"

#include "SireStream/datastream.h"

using namespace SireStream;
using namespace SireFF;

const char* missing_component::typeName()
{
    return QMetaType::typeName( qMetaTypeId<missing_component>() );
}

const char* missing_function::typeName()
{
    return QMetaType::typeName( qMetaTypeId<missing_function>() );
}

const char* missing_forcefield::typeName()
{
    return QMetaType::typeName( qMetaTypeId<missing_forcefield>() );
}

const char* missing_derivative::typeName()
{
    return QMetaType::typeName( qMetaTypeId<missing_derivative>() );
}

const char* missing_parameter::typeName()
{
    return QMetaType::typeName( qMetaTypeId<missing_parameter>() );
}

const char* duplicate_component::typeName()
{
    return QMetaType::typeName( qMetaTypeId<duplicate_component>() );
}

const char* duplicate_function::typeName()
{
    return QMetaType::typeName( qMetaTypeId<duplicate_function>() );
}

const char* duplicate_forcefield::typeName()
{
    return QMetaType::typeName( qMetaTypeId<duplicate_forcefield>() );
}

const char* invalid_group::typeName()
{
    return QMetaType::typeName( qMetaTypeId<invalid_group>() );
}

static const RegisterMetaType<missing_component> r_mcomp;
static const RegisterMetaType<missing_function> r_mfunc;
static const RegisterMetaType<missing_forcefield> r_mff;
static const RegisterMetaType<missing_derivative> r_mderiv;
static const RegisterMetaType<missing_parameter> r_mparam;
static const RegisterMetaType<duplicate_component> r_dupcomp;
static const RegisterMetaType<duplicate_function> r_dupfunc;
static const RegisterMetaType<duplicate_forcefield> r_dupff;
static const RegisterMetaType<invalid_group> r_group;
