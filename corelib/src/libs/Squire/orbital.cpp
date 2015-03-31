/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2009  Christopher Woods
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

#include "orbital.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace Squire;
using namespace SireBase;
using namespace SireStream;

////////////
//////////// Implementation of Orbital
////////////

static const RegisterMetaType<Orbital> r_orbital( MAGIC_ONLY, Orbital::typeName() );

/** Serialise to a binary datastream */
QDataStream SQUIRE_EXPORT &operator<<(QDataStream &ds, const Orbital &orbital)
{
    writeHeader(ds, r_orbital, 1);
    
    ds << static_cast<const Property&>(orbital);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SQUIRE_EXPORT &operator>>(QDataStream &ds, Orbital &orbital)
{
    VersionID v = readHeader(ds, r_orbital);
    
    if (v == 1)
    {
        ds >> static_cast<Property&>(orbital);
    }
    else
        throw version_error(v, "1", r_orbital, CODELOC);

    return ds;
}

/** Constructor */
Orbital::Orbital() : Property()
{}

/** Copy constructor */
Orbital::Orbital(const Orbital &other) : Property(other)
{}

/** Destructor */
Orbital::~Orbital()
{}

/** Copy assignment operator */
Orbital& Orbital::operator=(const Orbital &other)
{
    Property::operator=(other);
    return *this;
}

/** Comparison operator */
bool Orbital::operator==(const Orbital &other) const
{
    return Property::operator==(other);
}

/** Comparison operator */
bool Orbital::operator!=(const Orbital &other) const
{
    return Property::operator!=(other);
}

const char* Orbital::typeName()
{
    return "Squire::Orbital";
}

////////////
//////////// Implementation of OrbitalShell
////////////

static const RegisterMetaType<OrbitalShell> r_orbshell( MAGIC_ONLY, 
                                                        OrbitalShell::typeName() );
                                                        
/** Serialise to a binary datastream */
QDataStream SQUIRE_EXPORT &operator<<(QDataStream &ds, const OrbitalShell &orbshell)
{
    writeHeader(ds, r_orbshell, 1);
    
    ds << static_cast<const Orbital&>(orbshell);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SQUIRE_EXPORT &operator>>(QDataStream &ds, OrbitalShell &orbshell)
{
    VersionID v = readHeader(ds, r_orbshell);
    
    if (v == 1)
    {
        ds >> static_cast<Orbital&>(orbshell);
    }
    else
        throw version_error(v, "1", r_orbshell, CODELOC);

    return ds;
}

/** Constructor */
OrbitalShell::OrbitalShell() : Orbital()
{}

/** Copy constructor */
OrbitalShell::OrbitalShell(const OrbitalShell &other) : Orbital(other)
{}

/** Destructor */
OrbitalShell::~OrbitalShell()
{}

/** Copy assignment operator */
OrbitalShell& OrbitalShell::operator=(const OrbitalShell &other)
{
    Orbital::operator=(other);
    return *this;
}

/** Comparison operator */
bool OrbitalShell::operator==(const OrbitalShell &other) const
{
    return Orbital::operator==(other);
}

/** Comparison operator */
bool OrbitalShell::operator!=(const OrbitalShell &other) const
{
    return Orbital::operator!=(other);
}

const char* OrbitalShell::typeName()
{
    return "Squire::OrbitalShell";
}

////////////
//////////// Implementation of ShellPair
////////////

static const RegisterMetaType<ShellPair> r_shellpair( MAGIC_ONLY, ShellPair::typeName() );

/** Serialise to a binary datastream */
QDataStream SQUIRE_EXPORT &operator<<(QDataStream &ds, const ShellPair &pair)
{
    writeHeader(ds, r_shellpair, 1);
    
    ds << static_cast<const Property&>(pair);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SQUIRE_EXPORT &operator>>(QDataStream &ds, ShellPair &pair)
{
    VersionID v = readHeader(ds, r_shellpair);
    
    if (v == 1)
    {
        ds >> static_cast<Property&>(pair);
    }
    else
        throw version_error(v, "1", r_shellpair, CODELOC);
        
    return ds;
}

/** Constructor */
ShellPair::ShellPair() : Property()
{}

/** Copy constructor */
ShellPair::ShellPair(const ShellPair &other) : Property(other)
{}

/** Destructor */
ShellPair::~ShellPair()
{}

/** Copy assignment operator */
ShellPair& ShellPair::operator=(const ShellPair &other)
{
    Property::operator=(other);
    return *this;
}

/** Comparison operator */
bool ShellPair::operator==(const ShellPair &other) const
{
    return Property::operator==(other);
}

/** Comparison operator */
bool ShellPair::operator!=(const ShellPair &other) const
{
    return Property::operator!=(other);
}

const char* ShellPair::typeName()
{
    return "Squire::ShellPair";
}
