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

#include "SireMol/errors.h"

using namespace SireMol;

const char* missing_atom::typeName()
{
    return QMetaType::typeName( qMetaTypeId<missing_atom>() );
}

const char* duplicate_atom::typeName()
{
    return QMetaType::typeName( qMetaTypeId<duplicate_atom>() );
}

const char* missing_residue::typeName()
{
    return QMetaType::typeName( qMetaTypeId<missing_residue>() );
}

const char* duplicate_residue::typeName()
{
    return QMetaType::typeName( qMetaTypeId<duplicate_residue>() );
}

const char* missing_cutgroup::typeName()
{
    return QMetaType::typeName( qMetaTypeId<missing_cutgroup>() );
}

const char* duplicate_cutgroup::typeName()
{
    return QMetaType::typeName( qMetaTypeId<duplicate_cutgroup>() );
}

const char* missing_chain::typeName()
{
    return QMetaType::typeName( qMetaTypeId<missing_chain>() );
}

const char* duplicate_chain::typeName()
{
    return QMetaType::typeName( qMetaTypeId<duplicate_chain>() );
}

const char* missing_segment::typeName()
{
    return QMetaType::typeName( qMetaTypeId<missing_segment>() );
}

const char* duplicate_segment::typeName()
{
    return QMetaType::typeName( qMetaTypeId<duplicate_segment>() );
}

const char* missing_group::typeName()
{
    return QMetaType::typeName( qMetaTypeId<missing_group>() );
}

const char* duplicate_group::typeName()
{
    return QMetaType::typeName( qMetaTypeId<duplicate_group>() );
}

const char* missing_molecule::typeName()
{
    return QMetaType::typeName( qMetaTypeId<missing_molecule>() );
}

const char* duplicate_molecule::typeName()
{
    return QMetaType::typeName( qMetaTypeId<duplicate_molecule>() );
}

const char* template_error::typeName()
{
    return QMetaType::typeName( qMetaTypeId<template_error>() );
}

const char* anchor_error::typeName()
{
    return QMetaType::typeName( qMetaTypeId<anchor_error>() );
}

const char* ring_error::typeName()
{
    return QMetaType::typeName( qMetaTypeId<ring_error>() );
}

const char* incompatible_molecule::typeName()
{
    return QMetaType::typeName( qMetaTypeId<incompatible_molecule>() );
}

const char* missing_bond::typeName()
{
    return QMetaType::typeName( qMetaTypeId<missing_bond>() );
}

const char* missing_angle::typeName()
{
    return QMetaType::typeName( qMetaTypeId<missing_angle>() );
}

const char* missing_dihedral::typeName()
{
    return QMetaType::typeName( qMetaTypeId<missing_dihedral>() );
}

static const RegisterMetaType<missing_atom> r_matom;
static const RegisterMetaType<missing_group> r_mgrp;
static const RegisterMetaType<duplicate_atom> r_datom;
static const RegisterMetaType<missing_residue> r_mres;
static const RegisterMetaType<duplicate_residue> r_dres;
static const RegisterMetaType<missing_cutgroup> r_mcgrp;
static const RegisterMetaType<duplicate_cutgroup> r_dcgrp;
static const RegisterMetaType<missing_chain> r_mchn;
static const RegisterMetaType<duplicate_chain> r_dchn;
static const RegisterMetaType<missing_segment> r_mseg;
static const RegisterMetaType<duplicate_segment> r_dseg;
static const RegisterMetaType<missing_molecule> r_mmol;
static const RegisterMetaType<duplicate_molecule> r_dmol;
static const RegisterMetaType<template_error> r_tmplerr;
static const RegisterMetaType<anchor_error> r_ancerr;
static const RegisterMetaType<ring_error> r_ringerr;
static const RegisterMetaType<incompatible_molecule> r_incompat;
static const RegisterMetaType<missing_bond> r_mbond;
static const RegisterMetaType<missing_angle> r_mangle;
static const RegisterMetaType<missing_angle> r_mdihedral;
