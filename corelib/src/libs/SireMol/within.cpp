/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2010  Christopher Woods
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

#include "within.h"
#include "atomcoords.h"
#include "moleculeview.h"
#include "molecule.h"
#include "atom.h"
#include "mover.hpp"
#include "selector.hpp"

#include "SireUnits/units.h"

#include "SireMol/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMol;
using namespace SireMaths;
using namespace SireID;
using namespace SireUnits;
using namespace SireUnits::Dimension;
using namespace SireStream;

static const RegisterMetaType<Within> r_within;

/** Serialise to a binary datastream */
QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds, const Within &within)
{
    writeHeader(ds, r_within, 1);
    
    SharedDataStream sds(ds);
    
    sds << within.atomid << within.point << within.dist.to(angstrom);

    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds, Within &within)
{
    VersionID v = readHeader(ds, r_within);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        double dist;
        
        sds >> within.atomid >> within.point >> dist;
            
        within.dist = dist*angstrom;
    }
    else
        throw version_error( v, "1", r_within, CODELOC );
        
    return ds;
}

/** Null constructor */
Within::Within() : AtomID(), point(0), dist(0)
{}

/** Construct to match atoms within 'dist' of 'point' */
Within::Within(Length distance, const Vector &p)
       : AtomID(), point(p), dist(distance)
{}

/** Construct to match atoms within 'dist' of the atoms
    from the molecule that match 'atomid' */
Within::Within(Length distance, const AtomID &id)
       : AtomID(), atomid(id), point(0), dist(distance)
{}

/** Copy constructor */
Within::Within(const Within &other)
       : AtomID(other), atomid(other.atomid), point(other.point),
         dist(other.dist)
{}

/** Destructor */
Within::~Within()
{}

/** Is this selection null? */
bool Within::isNull() const
{
    return atomid.isNull() and dist.value() == 0 and point == Vector(0);
}

/** Return a hash of this identifier */
uint Within::hash() const
{
    return atomid.hash();
}
            
/** Return a string representatio of this ID */
QString Within::toString() const
{
    if (this->isNull())
        return QObject::tr("Within::null");

    else if (atomid.isNull())
        return QObject::tr("Within( %1 A from %2 )")
                    .arg(dist.to(angstrom)).arg(point.toString());

    else
        return QObject::tr("Within( %1 A from %2 )")
                    .arg(dist.to(angstrom)).arg(atomid.toString());
}

/** Copy assignment operator */
Within& Within::operator=(const Within &other)
{
    if (this != &other)
    {
        atomid = other.atomid;
        point = other.point;
        dist = other.dist;
    }
    
    return *this;
}

/** Comparison operator */
bool Within::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare<Within>(*this, other);
}

/** Comparison operator */
bool Within::operator==(const Within &other) const
{
    return atomid == other.atomid and point == other.point
               and dist == other.dist;
}

/** Comparison operator */
bool Within::operator!=(const Within &other) const
{
    return not Within::operator==(other);
}

/** Map this ID to the list of indicies of atoms that match this ID

    Note that this function is not valid for this ID class, as
    we need to have access to the molecular coordinates

    \throw SireError::incompatible_error
*/
QList<AtomIdx> Within::map(const MolInfo&) const
{
    throw SireError::incompatible_error( QObject::tr(
            "The ID %1 cannot be used in this context. It has to be "
            "used to select atoms in a molecule.")
                .arg(this->toString()), CODELOC );

    return QList<AtomIdx>();
}

/** Map this ID to the list of atomidxs of specified atoms 
    in the passed molecule
    
    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
QList<AtomIdx> Within::map(const MoleculeView &molview, const PropertyMap &map) const
{
    QList<AtomIdx> atomidxs;

    const double dist2 = dist.value() * dist.value();

    const MoleculeData &moldata = molview.data();

    const AtomCoords &coords = moldata.property(map["coordinates"])
                                      .asA<AtomCoords>();

    if (atomid.isNull())
    {
        if (molview.selectedAll())
        {
            for (AtomIdx i(0); i<coords.nAtoms(); ++i)
            {
                if (Vector::distance2(point, coords[moldata.info().cgAtomIdx(i)]) < dist2)
                    atomidxs.append(i);
            }
        }
        else
        {
            foreach (AtomIdx i, molview.selection().selectedAtoms())
            {
                if (Vector::distance2(point, coords[moldata.info().cgAtomIdx(i)]) < dist2)
                    atomidxs.append(i);
            }
        }
    }
    else
    {
        Selector<Atom> atoms = molview.molecule().selectAll(atomid,map);
        
        QVector<Vector> points;
        
        
        for (int i=0; i<atoms.count(); ++i)
        {
            points.append( coords[atoms[i].cgAtomIdx()] );
        }
        
        if (molview.selectedAll())
        {
            for (AtomIdx i(0); i<coords.nAtoms(); ++i)
            {
                foreach (const Vector &p, points)
                {
                    if (Vector::distance2(p, coords[moldata.info().cgAtomIdx(i)]) < dist2)
                    {
                        atomidxs.append(i);
                        break;
                    }
                }
            }
        }
        else
        {
            foreach (AtomIdx i, molview.selection().selectedAtoms())
            {
                foreach (const Vector &p, points)
                {
                    if (Vector::distance2(p, coords[moldata.info().cgAtomIdx(i)]) < dist2)
                    {
                        atomidxs.append(i);
                        break;
                    }
                }
            }
        }
    }
    
    if (atomidxs.isEmpty())
        throw SireMol::missing_atom( QObject::tr(
                "There is no atom that matches %1.").arg(this->toString()),
                    CODELOC );
                    
    return atomidxs;
}

const char* Within::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Within>() );
}

Within* Within::clone() const
{
    return new Within(*this);
}
