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

#ifndef SIREFF_PATCH_H
#define SIREFF_PATCH_H

#include "ffparameters.h"

#include "SireVol/coordgroup.h"

#include <QVarLengthArray>

SIRE_BEGIN_HEADER

namespace SireFF
{
class Patch;
}

QDataStream& operator<<(QDataStream&, const SireFF::Patch&);
QDataStream& operator>>(QDataStream&, SireFF::Patch&);

namespace SireFF
{

using SireVol::CoordGroupArray;
using SireVol::CoordGroup;
using SireVol::AABox;

class Patches;
class FFBead;
class FFBeadChange;

/** A patch is used by the forcefield classes to organise the
    beads into groups (normally spatially grouped, e.g. into
    patches of neighbouring groups). This allows the forcefield
    classes to optimise intermolecular energy evaluations by
    performing a domain decomposition.
    
    This class is publically read-only, though can be modified
    by a Patches class (as the Patches class is designed to
    modify the Patch objects that it contains). This is because
    the Patches class contains the Patcher object that contains
    the geometric information about this patch, and also because
    the Patches object manages the bead ID numbers used to rapidly
    identify each bead in the set of patches.
    
    @author Christopher Woods
*/
class SIREFF_EXPORT Patch 
        : public SireBase::ConcreteProperty<Patch,SireBase::Property>
{

friend SIREFF_EXPORT QDataStream& ::operator<<(QDataStream&, const Patch&);
friend SIREFF_EXPORT QDataStream& ::operator>>(QDataStream&, Patch&);

friend class Patches;

public:
    Patch();
    Patch(const Patch &other);
    
    ~Patch();
    
    static const char* typeName();
    
    Patch& operator=(const Patch &other);
    
    bool operator==(const Patch &other) const;
    bool operator!=(const Patch &other) const;
    
    QString toString() const;
    
    bool isEmpty() const;
    
    int nBeads() const;
    
    const AABox& aaBox() const;

    const QVector<quint32> beadIDs() const;
    const CoordGroupArray& coordinates() const;
    const FFParametersArray& parameters() const;

    int getLocation(quint32 beadid) const;

protected:
    /////////////////////////////////////////////////////////////////////////
    //Functions to modify the Patch, called only by the parent Patches class

    FFBead add(quint32 beadid, const CoordGroup &coords, const FFParameters &params);
    QHash<quint32,FFBead> add(const QVarLengthArray<quint32> &beadids, 
                              const CoordGroupArray &coords,
                              const FFParametersArray &params);
             
    FFBeadChange update(quint32 beadid, const CoordGroup &coords);
    FFBeadChange update(quint32 beadid, const FFParameters &params);
    FFBeadChange update(quint32 beadid, const CoordGroup &coords, 
                                        const FFParameters &params);
                
    QHash<quint32,FFBeadChange> update(const QVarLengthArray<quint32> &beadids, 
                                       const CoordGroupArray &coords);
    QHash<quint32,FFBeadChange> update(const QVarLengthArray<quint32> &beadids, 
                                       const FFParametersArray &params);
    QHash<quint32,FFBeadChange> update(const QVarLengthArray<quint32> &beadids, 
                                       const CoordGroupArray &coords,
                                       const FFParametersArray &params);
                
    FFBead remove(quint32 beadid);
    QHash<quint32,FFBead> remove(const QVarLengthArray<quint32> &beadids);

    void removeAll();

    //
    //////////////////////////////////////////////////////////////////////////
private:
    int getBeadIdx(quint32 beadid) const;

    /** The coordinates of all of the atoms in the patch,
        arranged into beads (CoordGroups). This object also
        contains the AABox that completely encompasses the 
        beads in the patch */
    CoordGroupArray coords;
    
    /** The parameters for all of the beads, in the same order
        as the coordinates */
    FFParametersArrayPtr params;
    
    /** Internal ID number of each bead */
    QVector<quint32> idx_to_beadid;
    
    /** Index of each bead in the arrays */
    QHash<quint32,int> beadid_to_idx;
    
    /** The AABox that completely encloses all of the 
        atoms in this patch */
    AABox aabox;
};

}

Q_DECLARE_METATYPE( SireFF::Patch )

SIRE_END_HEADER

#endif
