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

#ifndef SIREMOVE_SUPRASYSTEM_H
#define SIREMOVE_SUPRASYSTEM_H

#include "SireBase/property.h"

#include "suprasubsystem.h"

SIRE_BEGIN_HEADER

namespace SireMove
{
class SupraSystem;
}

SIREMOVE_EXPORT QDataStream& operator<<(QDataStream&, const SireMove::SupraSystem&);
SIREMOVE_EXPORT QDataStream& operator>>(QDataStream&, SireMove::SupraSystem&);

namespace SireSystem
{
class System;
}

namespace SireMove
{

using SireSystem::System;
using SireSystem::SystemMonitor;
using SireSystem::SystemMonitors;

class SimStore;
class Moves;

class SupraSubSystem;

/** This is the base class of all supra-systems. A supra-system is a 
    collection of systems, used to perform a simulation that
    involves making moves on lots of connected sub-systems. A good
    example is a replica exchange simulation, where each replica
    is a sub-system, and the ensemble of systems is the supra-system.
    
    @author Christopher Woods
*/
class SIREMOVE_EXPORT SupraSystem
       : public SireBase::ConcreteProperty<SupraSystem,SireBase::Property>
{

friend SIREMOVE_EXPORT QDataStream& ::operator<<(QDataStream&, const SupraSystem&);
friend SIREMOVE_EXPORT QDataStream& ::operator>>(QDataStream&, SupraSystem&);

public:
    SupraSystem();
    SupraSystem(int n);
    SupraSystem(const System &system, int n=1);
    SupraSystem(const QVector<System> &systems);
    
    SupraSystem(const SupraSubSystem &subsystem, int n=1);
    SupraSystem(const SupraSystem &other);
    
    virtual ~SupraSystem();
    
    SupraSystem& operator=(const SupraSystem &other);
    
    bool operator==(const SupraSystem &other) const;
    bool operator!=(const SupraSystem &other) const;
    
    static const char* typeName();
    
    virtual SupraSystem* clone() const;
    
    const SupraSubSystem& operator[](int i) const;
    
    const SupraSubSystem& at(int i) const;
    
    bool isEmpty() const;
    
    int nSubSystems() const;
    int count() const;
    int size();
    
    void clearStatistics();
    void clearSubStatistics();
    
    void clearAllStatistics();
    
    void mustNowRecalculateFromScratch();

    bool isPacked() const;
    bool isPackedToMemory() const;
    bool isPackedToDisk() const;
    
    bool anyPacked() const;
    bool anyPackedToMemory() const;
    bool anyPackedToDisk() const;
    
    void pack();
    void unpack();
    
    void packToDisk();
    void packToDisk(const QString &tempdir);
    
    void packToMemory();
    
    void pack(int i);
    void unpack(int i);
    
    void packToDisk(int i);
    void packToDisk(int i, const QString &tempdir);
    
    void packToMemory(int i);
    
    void setSubSystems(const SupraSystem &system);
    
    virtual void setSubSystem(int i, const SupraSubSystem &subsystem);
    void setSubSystem(const SupraSubSystem &subsystem);
    
    virtual void setSubSystem(int i, const System &system);
    void setSubSystem(const System &system);
    
    virtual void setSubMoves(int i, const Moves &moves);
    void setSubMoves(const Moves &moves);
    
    void setSubSystemAndMoves(int i, const System &system, const Moves &moves);
    void setSubSystemAndMoves(const System &system, const Moves &moves);
    
    virtual void setSubSystemAndMoves(int i, const SimStore &simstore);
    void setSubSystemAndMoves(const SimStore &simstore);
    
    void setSubMonitors(const SystemMonitors &monitors, int frequency=1);
    void setSubMonitors(int i, const SystemMonitors &monitors, int frequency=1);
    
    void add(const QString &name, const SystemMonitor &monitor, int frequency=1);
                     
    void add(int i, const QString &name, const SystemMonitor &monitor, int frequency=1);
                     
    void add(const SystemMonitors &monitors, int frequency=1);
    void add(int i, const SystemMonitors &monitors, int frequency=1);
    
    virtual void setNSubMoves(int i, int nmoves);
    void setNSubMoves(int nmoves);
    
    virtual void setRecordStatistics(int i, bool record_stats);
    void setRecordStatistics(bool record_stats);
    
    virtual void setRecordSubStatistics(int i, bool record_stats);
    void setRecordSubStatistics(bool record_stats);

    void setRecordAllStatistics(int i, bool record_stats);
    void setRecordAllStatistics(bool record_stats);

    virtual void collectSupraStats();

    static const SupraSystem& null();

protected:
    virtual void _pre_pack();
    virtual void _post_pack();
    
    virtual void _pre_unpack();
    virtual void _post_unpack();
    
    const SupraSubSystem& _pvt_subSystem(int i) const;
    SupraSubSystem& _pvt_subSystem(int i);

    void updateSubSystems(const SupraSubSystem *old_system,
                          const SupraSubSystem *new_system,
                          QSet<int> &done_systems);

private:
    /** The array of all sub systems */
    QVector<SupraSubSystemPtr> subsystems;
};

typedef SireBase::PropPtr<SupraSystem> SupraSystemPtr;

}

Q_DECLARE_METATYPE( SireMove::SupraSystem )

SIRE_EXPOSE_CLASS( SireMove::SupraSystem )

SIRE_EXPOSE_PROPERTY( SireMove::SupraSystemPtr, SireMove::SupraSystem )

SIRE_END_HEADER

#endif
