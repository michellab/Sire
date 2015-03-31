/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2008  Christopher Woods
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

#include <QSet>

#include "atomiccoords3d.h"
#include "atomicparameters.hpp"
#include "atomicparameters3d.hpp"

#include "SireMol/partialmolecule.h"
#include "SireMol/atomcoords.h"

#include "SireMol/mover.hpp"

#include "SireStream/shareddatastream.h"

#include "SireError/errors.h"

#include <QDebug>

using namespace SireFF;
using namespace SireFF::detail;

using namespace SireVol;
using namespace SireMol;

using namespace SireStream;

/** Serialise to a binary datastream */
QDataStream SIREFF_EXPORT &operator<<(QDataStream &ds,
                                      const AtomicCoords3D &coords)
{
    SharedDataStream sds(ds);

    sds << coords.coords;
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREFF_EXPORT &operator>>(QDataStream &ds,
                                      AtomicCoords3D &coords)
{
    SharedDataStream sds(ds);

    sds >> coords.coords;
    
    return ds;
}

bool SIREFF_EXPORT SireFF::detail::selectedAll(const QSet<quint32> &changed_groups, 
                                               quint32 n)
{
    if (quint32(changed_groups.count()) >= n)
    {
        quint32 got = 0;
    
        for (QSet<quint32>::const_iterator it = changed_groups.begin();
             it != changed_groups.end();
             ++it)
        {
            if (*it < n)
            {
                ++got;
                
                if (got == n)
                    return true;
            }
        }
    }
    
    return false;
}

/** Null constructor */
AtomicCoords3D::AtomicCoords3D()
{}

/** Construct to hold the passed coordinates */
AtomicCoords3D::AtomicCoords3D(const CoordGroupArray &coordinates)
               : coords(coordinates)
{}

/** Construct from the atomic coordinates of the selected view
    of the molecule 'molecule', using the property 'coords_property'
    to find the coordinates property
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
AtomicCoords3D::AtomicCoords3D(const PartialMolecule &molecule,
                               const PropertyName &coords_property)
{
    //get the coordinates of the atoms
    coords = molecule.property(coords_property).asA<AtomCoords>().array();
    
    //these are for all of the atoms - we now need to mask this
    //so that only the coordinates for atoms in CutGroups that
    //contain at least one selected atom are stored
    if (not molecule.selection().selectedAllCutGroups())
    {
        //pick out the selected CoordGroups
        QList<CGIdx> selected_cgroups = molecule.selection().selectedCutGroups();
        
        QVector<CoordGroup> cgroups;
        cgroups.reserve(selected_cgroups.count());
        
        foreach (CGIdx selected_cgroup, selected_cgroups)
        {
            cgroups.append( coords.constData()[selected_cgroup] );
        }
        
        //rebuild the CoordGroupArray
        coords = CoordGroupArray(cgroups);
    }
}
               
/** Copy constructor */
AtomicCoords3D::AtomicCoords3D(const AtomicCoords3D &other)
               : coords(other.coords)
{}

/** Destructor */
AtomicCoords3D::~AtomicCoords3D()
{}

/** Copy assignment operator */
AtomicCoords3D& AtomicCoords3D::operator=(const AtomicCoords3D &other)
{
    coords = other.coords;
    return *this;
}

/** Comparison operator */
bool AtomicCoords3D::operator==(const AtomicCoords3D &other) const
{
    return coords == other.coords;
}

/** Comparison operator */
bool AtomicCoords3D::operator!=(const AtomicCoords3D &other) const
{
    return coords != other.coords;
}

/** Set the coordinates */
void AtomicCoords3D::setAtomicCoordinates(const AtomicCoords3D &other)
{
    coords = other.coords;
}

/** Have all of the parameters in this group changed? */
bool AtomicCoords3D::changedAllGroups(const AtomicCoords3D &other) const
{
    int ngroups = qMin(coords.count(), other.coords.count());
    
    const CoordGroup *this_array = coords.constData();
    const CoordGroup *other_array = other.coords.constData();
    
    for (int i=0; i<ngroups; ++i)
    {
        if (this_array[i] == other_array[i])
        {
            return false;
        }
    }
            
    return true;
}

/** Add the indicies of CoordGroups in this array that are different
    to those in 'params' */
void AtomicCoords3D::addChangedGroups(const AtomicCoords3D &params,
                                      QSet<quint32> &changed_groups) const
{
    quint32 ngroups = coords.count();

    if ( selectedAll(changed_groups, ngroups) )
        return;
        
    quint32 nsharedgroups = qMin(ngroups, quint32(params.coords.count()));
    
    const CoordGroup *this_cgroups = coords.constData();
    const CoordGroup *other_cgroups = params.coords.constData();
    
    for (quint32 i=0; i<nsharedgroups; ++i)
    {
        if (not changed_groups.contains(i)
            and this_cgroups[i] != other_cgroups[i])
        {
            changed_groups.insert(i);
            
            if (selectedAll(changed_groups, ngroups))
                return;
        }
    }
}

/** Return the indicies of CoordGroups in this set that are different to
    those in 'params'. */
QSet<quint32> AtomicCoords3D::getChangedGroups(const AtomicCoords3D &params) const
{
    QSet<quint32> diffgroups;
    diffgroups.reserve(coords.count());
    
    this->addChangedGroups(params, diffgroups);
    return diffgroups;
}
    
/** Mask this set of coordinates so that only the CoordGroups at indicies
    'idxs' are present */
AtomicCoords3D AtomicCoords3D::applyMask(const QSet<quint32> &idxs) const
{
    if (selectedAll(idxs, coords.count()))
        return *this;
    
    else if (idxs.isEmpty())
        return AtomicCoords3D();
        
    else if (idxs.count() == 1)
        return AtomicCoords3D( coords[*(idxs.constBegin())] );

    //mask by the indicies
    quint32 ngroups = coords.count();
    
    QVector<CoordGroup> cgroups;
    cgroups.reserve(idxs.count());
    
    const CoordGroup *cgroups_array = coords.constData();
    
    for (quint32 i=0; i<ngroups; ++i)
    {
        if (idxs.contains(i))
            cgroups.append( cgroups_array[i] );
    }
    
    return AtomicCoords3D( CoordGroupArray(cgroups) );
}
