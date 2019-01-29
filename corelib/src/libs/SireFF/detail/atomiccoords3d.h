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

#ifndef SIREFF_DETAIL_ATOMICCOORDS3D_H
#define SIREFF_DETAIL_ATOMICCOORDS3D_H

#include "SireVol/coordgroup.h"

SIRE_BEGIN_HEADER

namespace SireFF { namespace detail { class AtomicCoords3D; } }

SIREFF_EXPORT QDataStream& operator<<(QDataStream&, const SireFF::detail::AtomicCoords3D&);
SIREFF_EXPORT QDataStream& operator>>(QDataStream&, SireFF::detail::AtomicCoords3D&);

namespace SireVol
{
class CoordGroup;
class CoordGroupArray;
class CoordGroupArrayArray;
}

namespace SireMol
{
class PartialMolecule;
}

namespace SireBase
{
class PropertyName;
}

namespace SireFF
{

using SireVol::CoordGroup;
using SireVol::CoordGroupArray;
using SireVol::CoordGroupArrayArray;

using SireMol::PartialMolecule;
using SireBase::PropertyName;

namespace detail
{
    
SIREFF_EXPORT bool selectedAll(const QSet<quint32> &idxs, quint32 n);

/** This class holds the 3D coordinates of all of the atoms,
    optimised for access by the forcefield
    
    @author Christopher Woods
*/
class SIREFF_EXPORT AtomicCoords3D
{

friend QDataStream& ::operator<<(QDataStream&, const AtomicCoords3D&);
friend QDataStream& ::operator>>(QDataStream&, AtomicCoords3D&);

public:
    AtomicCoords3D();
    
    AtomicCoords3D(const PartialMolecule &molecule,
                   const PropertyName &coords_property);
    
    AtomicCoords3D(const CoordGroupArray &coords);
                   
    AtomicCoords3D(const AtomicCoords3D &other);
    
    ~AtomicCoords3D();
    
    AtomicCoords3D& operator=(const AtomicCoords3D &other);
    
    bool operator==(const AtomicCoords3D &other) const;
    bool operator!=(const AtomicCoords3D &other) const;
    
    const CoordGroupArray& atomicCoordinates() const;
    
    void setAtomicCoordinates(const AtomicCoords3D &other);
    
    bool changedAllGroups(const AtomicCoords3D &params) const;
    
    QSet<quint32> getChangedGroups(const AtomicCoords3D &params) const;
    
    void addChangedGroups(const AtomicCoords3D &params,
                          QSet<quint32> &changed_groups) const;
    
    AtomicCoords3D applyMask(const QSet<quint32> &idxs) const;
    
private:
    /** The 3D coordinates of all of the atoms in this forcefield,
        arranged by CutGroup */
    CoordGroupArray coords;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Return the 3D coordinates of the atoms in the forcefield */
inline const CoordGroupArray& AtomicCoords3D::atomicCoordinates() const
{
    return coords;
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

} // end of namespace detail

} // end of namespace SireFF

SIRE_END_HEADER

#endif
