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

#ifndef SIREFF_PATCHES_H
#define SIREFF_PATCHES_H

#include "patch.h"

#include "SireVol/patching.h"

SIRE_BEGIN_HEADER

namespace SireFF
{
class Patches;
class FFBead;
class FFBeadChange;
}

QDataStream& operator<<(QDataStream&, const SireFF::Patches&);
QDataStream& operator>>(QDataStream&, SireFF::Patches&);

QDataStream& operator<<(QDataStream&, const SireFF::FFBead&);
QDataStream& operator>>(QDataStream&, SireFF::FFBead&);

QDataStream& operator<<(QDataStream&, const SireFF::FFBeadChange&);
QDataStream& operator>>(QDataStream&, SireFF::FFBeadChange&);

namespace SireFF
{

using SireVol::Patching;
using SireVol::Space;
using SireVol::CoordGroup;
using SireVol::AABox;

/** This simple class holds the coordinates and forcefield parameters 
    for a bead
    
    @author Christopher Woods
*/
class SIREFF_EXPORT FFBead
{

friend SIREFF_EXPORT QDataStream& ::operator<<(QDataStream&, const FFBead&);
friend SIREFF_EXPORT QDataStream& ::operator>>(QDataStream&, FFBead&);

public:
    FFBead();
    FFBead(const CoordGroup &coordinates, const FFParameters &parameters);
    
    FFBead(const FFBead &other);
    
    ~FFBead();
    
    FFBead& operator=(const FFBead &other);
    
    bool operator==(const FFBead &other) const;
    bool operator!=(const FFBead &other) const;
    
    static const char* typeName();
    
    const CoordGroup& coordinates() const;
    const FFParameters& parameters() const;
    
    bool isEmpty() const;
    
private:
    /** The bead coordinates */
    CoordGroup coords;
    
    /** The bead parameters */
    FFParametersPtr params;
};

/** This simple class holds the data the describes a change
    in a bead
    
    @author Christopher Woods
*/
class SIREFF_EXPORT FFBeadChange
{

friend SIREFF_EXPORT QDataStream& ::operator<<(QDataStream&, const FFBeadChange&);
friend SIREFF_EXPORT QDataStream& ::operator>>(QDataStream&, FFBeadChange&);

public:
    FFBeadChange();
    FFBeadChange(const FFBead &old_bead, const FFBead &new_bead);
    
    FFBeadChange(const FFBeadChange &other);
    
    ~FFBeadChange();
    
    FFBeadChange& operator=(const FFBeadChange &other);
    
    bool operator==(const FFBeadChange &other) const;
    bool operator!=(const FFBeadChange &other) const;
    
    static const char* typeName();
    
    bool isEmpty() const;
    
    const FFBead& oldBead() const;
    const FFBead& newBead() const;
    
    FFBeadChange update(const FFBead &newer_bead) const;
    
private:
    /** The old state of the bead */
    FFBead old_bead;
    
    /** The new state of the bead */
    FFBead new_bead;
};

/** This class holds a collection of Patch objects. All of the beads
    in a forcefield are arranged into Patches, with each patch containing
    a neighbouring group of beads.
    
    @author Christopher Woods
*/
class SIREFF_EXPORT Patches
        : public SireBase::ConcreteProperty<Patches,SireBase::Property>
{

friend SIREFF_EXPORT QDataStream& ::operator<<(QDataStream&, const Patches&);
friend SIREFF_EXPORT QDataStream& ::operator>>(QDataStream&, Patches&);

public:
    Patches();
    Patches(const Space &space, const Patching &patching);
    
    Patches(const Patching &patching);

    Patches(const Patches &other);

    ~Patches();
    
    Patches& operator=(const Patches &other);
    
    bool operator==(const Patches &other) const;
    bool operator!=(const Patches &other) const;
    
    QString toString() const;
    
    const Patching& patching() const;
    const Space& space() const;
    
    void repatch(const Patching &patching);
    void repatch(const Space &space);
   
    int nPatches() const;
    int size() const;
    int count() const;
    
    bool isEmpty() const;
    
    int nBeads() const;
    
    const Patch& operator[](int i) const;
    const Patch& at(int i) const;
    
    const Patch* data() const;
    const Patch* constData() const;

    QPair<int,int> getLocation(quint32 beadid) const;

    QPair<quint32,FFBead> add(const CoordGroup &coords, const FFParameters &params);
    QPair<quint32,FFBead> add(const FFBead &bead);
    
    QHash<quint32,FFBead> add(const CoordGroupArray &coords, 
                              const FFParametersArray &params);
    
    QHash<quint32,FFBead> add(const QVector<CoordGroup> &coords,
                              const QVector<FFParametersPtr> &params);

    QHash<quint32,FFBead> add(const QVector<FFBead> &beads);
    
    FFBeadChange update(quint32 beadid, const CoordGroup &coords);
    FFBeadChange update(quint32 beadid, const FFParameters &params);
    FFBeadChange update(quint32 beadid, const CoordGroup &coords, 
                        const FFParameters &params);
    FFBeadChange update(quint32 beadid, const FFBead &bead);

    QHash<quint32,FFBeadChange> update(const QVector<quint32> &beadids, 
                                       const CoordGroupArray &coords);
                                 
    QHash<quint32,FFBeadChange> update(const QVector<quint32> &beadids, 
                                       const FFParametersArray &params);
                                 
    QHash<quint32,FFBeadChange> update(const QVector<quint32> &beadids, 
                                       const CoordGroupArray &coords,
                                       const FFParametersArray &params);
    
    FFBeadChange remove(quint32 beadid);
    QHash<quint32,FFBeadChange> remove(const QVector<quint32> &beadids);

    void removeAll();
    
private:
    int getIdx(quint32 beadid) const;

    /** The current patching scheme (contains the space as well) */
    SireVol::PatchingPtr ptchng;

    /** All of the patches */
    QVector<Patch> ptchs;
    
    /** Index of the patch that contains each bead */
    QHash<quint32,quint32> beadid_to_patch;
    
    /** The last assigned bead ID */
    quint32 last_beadid;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Return the coordinates of this bead */
inline const CoordGroup& FFBead::coordinates() const
{
    return coords;
}

/** Return the parameters of the bead */
inline const FFParameters& FFBead::parameters() const
{
    return params.read();
}

/** Return whether or not this bead is empty */
inline bool FFBead::isEmpty() const
{
    return coords.isEmpty();
}

/** Return whether the change is empty (represents no change) */    
inline bool FFBeadChange::isEmpty() const
{
    return old_bead.isEmpty() and new_bead.isEmpty();
}

/** Return the state of the bead before the change */
inline const FFBead& FFBeadChange::oldBead() const
{
    return old_bead;
}

/** Return the state of the bead after the change */
inline const FFBead& FFBeadChange::newBead() const
{
    return new_bead;
}

#endif // SIRE_SKIP_INLINE_FUNCTIONS

}

Q_DECLARE_METATYPE( SireFF::FFBead )
Q_DECLARE_METATYPE( SireFF::FFBeadChange )
Q_DECLARE_METATYPE( SireFF::Patches )

SIRE_END_HEADER

#endif
