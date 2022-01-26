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

#ifndef SIREMOL_ERRORS_H
#define SIREMOL_ERRORS_H

#include "SireError/exception.h"

SIRE_BEGIN_HEADER

/**
This file contains the exceptions that can be thrown by the SireMol library.

@author Christopher Woods
*/

namespace SireMol
{

/** This is the base class of all SireMol errors */
class SIREMOL_EXPORT siremol_error : public SireError::exception
{
public:
    siremol_error() : exception()
    {}

    siremol_error(QString err, QString place = QString()) : exception(err,place)
    {}

    siremol_error(const siremol_error &other) : exception(other)
    {}

    ~siremol_error() throw()
    {}

    static const char* typeName()
    {
        return "SireMol::siremol_error";
    }
};


/** This exception is thrown when a request is made of a non-existant atom

    @author Christopher Woods
*/
class SIREMOL_EXPORT missing_atom : public siremol_error
{
public:
    missing_atom() : siremol_error()
    {}

    missing_atom(QString err, QString place = QString())
              : siremol_error(err,place)
    {}

    missing_atom(const missing_atom &other) : siremol_error(other)
    {}

    ~missing_atom() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return missing_atom::typeName();
    }

    void throwSelf() const
    {
        throw missing_atom(*this);
    }
};

/** This exception is thrown whenever multiple atoms are found,
    but we want only one

    @author Christopher Woods
*/
class SIREMOL_EXPORT duplicate_atom : public siremol_error
{
public:
    duplicate_atom() : siremol_error()
    {}

    duplicate_atom(QString err, QString place = QString())
                : siremol_error(err,place)
    {}

    duplicate_atom(const duplicate_atom &other) : siremol_error(other)
    {}

    ~duplicate_atom() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return duplicate_atom::typeName();
    }

    void throwSelf() const
    {
        throw duplicate_atom(*this);
    }
};

/** This exception is thrown when a request is made of a non-existant residue

    @author Christopher Woods
*/
class SIREMOL_EXPORT missing_residue : public siremol_error
{
public:
    missing_residue() : siremol_error()
    {}

    missing_residue(QString err, QString place = QString())
                    : siremol_error(err,place)
    {}

    missing_residue(const missing_residue &other) : siremol_error(other)
    {}

    ~missing_residue() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return missing_residue::typeName();
    }

    void throwSelf() const
    {
        throw missing_residue(*this);
    }
};

/** This exception is thrown when we get more than one
    residue, but we want only one

    @author Christopher Woods
*/
class SIREMOL_EXPORT duplicate_residue : public siremol_error
{
public:
    duplicate_residue() : siremol_error()
    {}

    duplicate_residue(QString err, QString place = QString())
                  : siremol_error(err,place)
    {}

    duplicate_residue(const duplicate_residue &other) : siremol_error(other)
    {}

    ~duplicate_residue() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return duplicate_residue::typeName();
    }

    void throwSelf() const
    {
        throw duplicate_residue(*this);
    }
};

/** This exception is thrown when an action is requested of a non-existant CutGroup

    @author Christopher Woods
*/
class SIREMOL_EXPORT missing_cutgroup : public siremol_error
{
public:
    missing_cutgroup() : siremol_error()
    {}

    missing_cutgroup(QString err, QString place = QString())
                    : siremol_error(err,place)
    {}

    missing_cutgroup(const missing_cutgroup &other) : siremol_error(other)
    {}

    ~missing_cutgroup() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return missing_cutgroup::typeName();
    }

    void throwSelf() const
    {
        throw missing_cutgroup(*this);
    }
};

/** This exception is thrown when we get more than one
    cutgroup, but we want only one

    @author Christopher Woods
*/
class SIREMOL_EXPORT duplicate_cutgroup : public siremol_error
{
public:
    duplicate_cutgroup() : siremol_error()
    {}

    duplicate_cutgroup(QString err, QString place = QString())
                  : siremol_error(err,place)
    {}

    duplicate_cutgroup(const duplicate_cutgroup &other) : siremol_error(other)
    {}

    ~duplicate_cutgroup() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return duplicate_cutgroup::typeName();
    }

    void throwSelf() const
    {
        throw duplicate_cutgroup(*this);
    }
};

/** This exception is thrown when an action is requested of a non-existant chain

    @author Christopher Woods
*/
class SIREMOL_EXPORT missing_chain : public siremol_error
{
public:
    missing_chain() : siremol_error()
    {}

    missing_chain(QString err, QString place = QString())
                    : siremol_error(err,place)
    {}

    missing_chain(const missing_chain &other) : siremol_error(other)
    {}

    ~missing_chain() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return missing_chain::typeName();
    }

    void throwSelf() const
    {
        throw missing_chain(*this);
    }
};

/** This exception is thrown when we get more than one
    chain, but we want only one

    @author Christopher Woods
*/
class SIREMOL_EXPORT duplicate_chain : public siremol_error
{
public:
    duplicate_chain() : siremol_error()
    {}

    duplicate_chain(QString err, QString place = QString())
                  : siremol_error(err,place)
    {}

    duplicate_chain(const duplicate_chain &other) : siremol_error(other)
    {}

    ~duplicate_chain() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return duplicate_chain::typeName();
    }

    void throwSelf() const
    {
        throw duplicate_chain(*this);
    }
};

/** This exception is thrown when an action is requested of a non-existant segment

    @author Christopher Woods
*/
class SIREMOL_EXPORT missing_segment : public siremol_error
{
public:
    missing_segment() : siremol_error()
    {}

    missing_segment(QString err, QString place = QString())
                    : siremol_error(err,place)
    {}

    missing_segment(const missing_segment &other) : siremol_error(other)
    {}

    ~missing_segment() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return missing_segment::typeName();
    }

    void throwSelf() const
    {
        throw missing_segment(*this);
    }
};

/** This exception is thrown when we get more than one
    segment, but we want only one

    @author Christopher Woods
*/
class SIREMOL_EXPORT duplicate_segment : public siremol_error
{
public:
    duplicate_segment() : siremol_error()
    {}

    duplicate_segment(QString err, QString place = QString())
                  : siremol_error(err,place)
    {}

    duplicate_segment(const duplicate_segment &other) : siremol_error(other)
    {}

    ~duplicate_segment() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return duplicate_segment::typeName();
    }

    void throwSelf() const
    {
        throw duplicate_segment(*this);
    }
};

/** This exception is thrown when a request is made of a non-existant group

    @author Christopher Woods
*/
class SIREMOL_EXPORT missing_group : public siremol_error
{
public:
    missing_group() : siremol_error()
    {}

    missing_group(QString err, QString place = QString())
              : siremol_error(err,place)
    {}

    missing_group(const missing_group &other) : siremol_error(other)
    {}

    ~missing_group() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return missing_group::typeName();
    }

    void throwSelf() const
    {
        throw missing_group(*this);
    }
};

/** This exception is thrown when a new molecule group is added
    to a set that already contains a group with the same ID, or
    when multiple groups match the same ID, but only one group
    is required.

    @author Christopher Woods
*/
class SIREMOL_EXPORT duplicate_group : public siremol_error
{
public:
    duplicate_group() : siremol_error()
    {}

    duplicate_group(QString err, QString place = QString())
                  : siremol_error(err,place)
    {}

    duplicate_group(const duplicate_group &other) : siremol_error(other)
    {}

    ~duplicate_group() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return duplicate_group::typeName();
    }

    void throwSelf() const
    {
        throw duplicate_group(*this);
    }
};

/** This exception is thrown when an action is requested of a non-existant molecule

    @author Christopher Woods
*/
class SIREMOL_EXPORT missing_molecule : public siremol_error
{
public:
    missing_molecule() : siremol_error()
    {}

    missing_molecule(QString err, QString place = QString())
                    : siremol_error(err,place)
    {}

    missing_molecule(const missing_molecule &other) : siremol_error(other)
    {}

    ~missing_molecule() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return missing_molecule::typeName();
    }

    void throwSelf() const
    {
        throw missing_molecule(*this);
    }
};

/** This exception is thrown when a new molecule is added to a set that already
    contains a molecule with the same ID

    @author Christopher Woods
*/
class SIREMOL_EXPORT duplicate_molecule : public siremol_error
{
public:
    duplicate_molecule() : siremol_error()
    {}

    duplicate_molecule(QString err, QString place = QString())
                  : siremol_error(err,place)
    {}

    duplicate_molecule(const duplicate_molecule &other) : siremol_error(other)
    {}

    ~duplicate_molecule() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return duplicate_molecule::typeName();
    }

    void throwSelf() const
    {
        throw duplicate_molecule(*this);
    }
};

/** This exception is thrown when a problem is detected when applying
    a molecule or residue template.

    @author Christopher Woods
*/
class SIREMOL_EXPORT template_error : public siremol_error
{
public:
    template_error() : siremol_error()
    {}

    template_error(QString err, QString place = QString())
                  : siremol_error(err,place)
    {}

    template_error(const template_error &other) : siremol_error(other)
    {}

    ~template_error() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return template_error::typeName();
    }

    void throwSelf() const
    {
        throw template_error(*this);
    }
};

/** This exception is thrown when a problem is detected with anchors on moving
    atoms

    @author Christopher Woods
*/
class SIREMOL_EXPORT anchor_error : public siremol_error
{
public:
    anchor_error() : siremol_error()
    {}

    anchor_error(QString err, QString place = QString())
                  : siremol_error(err,place)
    {}

    anchor_error(const anchor_error &other) : siremol_error(other)
    {}

    ~anchor_error() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return anchor_error::typeName();
    }

    void throwSelf() const
    {
        throw anchor_error(*this);
    }
};

/** This exception is thrown when a ring is detected when trying to split
    a molecule into two

    @author Christopher Woods
*/
class SIREMOL_EXPORT ring_error : public siremol_error
{
public:
    ring_error() : siremol_error()
    {}

    ring_error(QString err, QString place = QString())
                  : siremol_error(err,place)
    {}

    ring_error(const ring_error &other) : siremol_error(other)
    {}

    ~ring_error() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return ring_error::typeName();
    }

    void throwSelf() const
    {
        throw ring_error(*this);
    }
};

/** This exception is thrown when an incompatible molecule is used

    @author Christopher Woods
*/
class SIREMOL_EXPORT incompatible_molecule : public siremol_error
{
public:
    incompatible_molecule() : siremol_error()
    {}

    incompatible_molecule(QString err, QString place = QString())
                          : siremol_error(err,place)
    {}

    incompatible_molecule(const incompatible_molecule &other)
                          : siremol_error(other)
    {}

    ~incompatible_molecule() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return incompatible_molecule::typeName();
    }

    void throwSelf() const
    {
        throw incompatible_molecule(*this);
    }
};

/** This exception is thrown when a request is made of a non-existant bond

    @author Julien Michel
*/
class SIREMOL_EXPORT missing_bond : public siremol_error
{
public:
    missing_bond() : siremol_error()
    {}

    missing_bond(QString err, QString place = QString())
              : siremol_error(err,place)
    {}

    missing_bond(const missing_bond &other) : siremol_error(other)
    {}

    ~missing_bond() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return missing_bond::typeName();
    }

    void throwSelf() const
    {
        throw missing_bond(*this);
    }
};

/** This exception is thrown when a request is made of a non-existant angle

    @author Julien Michel
*/
class SIREMOL_EXPORT missing_angle : public siremol_error
{
public:
    missing_angle() : siremol_error()
    {}

    missing_angle(QString err, QString place = QString())
              : siremol_error(err,place)
    {}

    missing_angle(const missing_angle &other) : siremol_error(other)
    {}

    ~missing_angle() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return missing_angle::typeName();
    }

    void throwSelf() const
    {
        throw missing_angle(*this);
    }
};

/** This exception is thrown when a request is made of a non-existant dihedral

    @author Julien Michel
*/
class SIREMOL_EXPORT missing_dihedral : public siremol_error
{
public:
    missing_dihedral() : siremol_error()
    {}

    missing_dihedral(QString err, QString place = QString())
              : siremol_error(err,place)
    {}

    missing_dihedral(const missing_dihedral &other) : siremol_error(other)
    {}

    ~missing_dihedral() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return missing_dihedral::typeName();
    }

    void throwSelf() const
    {
        throw missing_dihedral(*this);
    }
};

}

Q_DECLARE_METATYPE(SireMol::missing_atom)
Q_DECLARE_METATYPE(SireMol::duplicate_atom)
Q_DECLARE_METATYPE(SireMol::missing_residue)
Q_DECLARE_METATYPE(SireMol::duplicate_residue)
Q_DECLARE_METATYPE(SireMol::missing_cutgroup)
Q_DECLARE_METATYPE(SireMol::duplicate_cutgroup)
Q_DECLARE_METATYPE(SireMol::missing_chain)
Q_DECLARE_METATYPE(SireMol::duplicate_chain)
Q_DECLARE_METATYPE(SireMol::missing_segment)
Q_DECLARE_METATYPE(SireMol::duplicate_segment)
Q_DECLARE_METATYPE(SireMol::missing_group)
Q_DECLARE_METATYPE(SireMol::duplicate_group)
Q_DECLARE_METATYPE(SireMol::missing_molecule)
Q_DECLARE_METATYPE(SireMol::duplicate_molecule)
Q_DECLARE_METATYPE(SireMol::template_error)
Q_DECLARE_METATYPE(SireMol::anchor_error)
Q_DECLARE_METATYPE(SireMol::ring_error)
Q_DECLARE_METATYPE(SireMol::incompatible_molecule)
Q_DECLARE_METATYPE(SireMol::missing_bond)
Q_DECLARE_METATYPE(SireMol::missing_angle)
Q_DECLARE_METATYPE(SireMol::missing_dihedral)

SIRE_END_HEADER

#endif
