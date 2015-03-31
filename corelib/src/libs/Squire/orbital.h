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

#ifndef SQUIRE_ORBITAL_H
#define SQUIRE_ORBITAL_H

#include "SireBase/property.h"

SIRE_BEGIN_HEADER

namespace Squire
{
class Orbital;
class OrbitalShell;
class ShellPair;
}

QDataStream& operator<<(QDataStream&, const Squire::Orbital&);
QDataStream& operator>>(QDataStream&, Squire::Orbital&);

QDataStream& operator<<(QDataStream&, const Squire::OrbitalShell&);
QDataStream& operator>>(QDataStream&, Squire::OrbitalShell&);

QDataStream& operator<<(QDataStream&, const Squire::ShellPair&);
QDataStream& operator>>(QDataStream&, Squire::ShellPair&);

namespace Squire
{

/** This is the base class of all orbitals.

    Orbital classes (like OrbitalShell) are designed to hold
    information about an orbital, and are not meant to be
    used directly in an integral program (indeed, an orbital
    normally doesn't even have coordinate information, as
    this is provided by the molecule that holds the orbital).
    
    Orbital classes are thus virtual and not designed for speed.

    @author Christopher Woods
*/
class SQUIRE_EXPORT Orbital : public SireBase::Property
{

friend QDataStream& ::operator<<(QDataStream&, const Orbital&);
friend QDataStream& ::operator>>(QDataStream&, Orbital&);

public:
    Orbital();
    
    Orbital(const Orbital &other);
    
    virtual ~Orbital();
    
    static const char* typeName();
    
    virtual Orbital* clone() const=0;
    
    virtual QString toString() const=0;

protected:
    Orbital& operator=(const Orbital &other);
    
    bool operator==(const Orbital &other) const;
    bool operator!=(const Orbital &other) const;
};

/** This is the base class of a shell of orbitals. An orbital
    shell contains the typical atomic orbitals (e.g. a p-shell
    contains 3 orbitals, p_x, p_y and p_z)

    OrbitalShell classes (like Orbital) are designed to hold
    information about an orbital, and are not meant to be
    used directly in an integral program (indeed, an orbital
    normally doesn't even have coordinate information, as
    this is provided by the molecule that holds the orbital).
    
    Orbital classes are thus virtual and not designed for speed.
    
    @author Christopher Woods
*/
class SQUIRE_EXPORT OrbitalShell : public Orbital
{

friend QDataStream& ::operator<<(QDataStream&, const OrbitalShell&);
friend QDataStream& ::operator>>(QDataStream&, OrbitalShell&);

public:
    OrbitalShell();
    OrbitalShell(const OrbitalShell &other);
    
    virtual ~OrbitalShell();
    
    static const char* typeName();
    
    virtual OrbitalShell* clone() const=0;
    
    virtual int angularMomentum() const=0;
    virtual int nOrbitals() const=0;

protected:
    OrbitalShell& operator=(const OrbitalShell &other);
    
    bool operator==(const OrbitalShell &other) const;
    bool operator!=(const OrbitalShell &other) const;
};

/** The base class of all combined pair orbitals */
class SQUIRE_EXPORT ShellPair : public SireBase::Property
{

friend QDataStream& ::operator<<(QDataStream&, const ShellPair&);
friend QDataStream& ::operator>>(QDataStream&, ShellPair&);

public:
    ShellPair();
    ShellPair(const ShellPair &other);
    
    virtual ~ShellPair();
    
    static const char* typeName();

    virtual ShellPair* clone() const=0;
    
    virtual int angularMomentum0() const=0;
    virtual int angularMomentum1() const=0;
    
    virtual int nOrbitals0() const=0;
    virtual int nOrbitals1() const=0;

protected:
    ShellPair& operator=(const ShellPair &other);
    
    bool operator==(const ShellPair &other) const;
    bool operator!=(const ShellPair &other) const;
};

}

SIRE_EXPOSE_CLASS( Squire::Orbital )
SIRE_EXPOSE_CLASS( Squire::OrbitalShell )

SIRE_EXPOSE_CLASS( Squire::ShellPair )

SIRE_END_HEADER

#endif
