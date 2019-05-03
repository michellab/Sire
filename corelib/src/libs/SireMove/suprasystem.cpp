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

#include "suprasystem.h"

#include "simstore.h"
#include "moves.h"

#include "SireSystem/system.h"
#include "SireSystem/systemmonitor.h"

#include "SireID/index.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireError/errors.h"

using namespace SireMove;
using namespace SireSystem;
using namespace SireID;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<SupraSystem> r_suprasystem;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const SupraSystem &suprasystem)
{
    writeHeader(ds, r_suprasystem, 1);
    
    SharedDataStream sds(ds);
    
    sds << suprasystem.subsystems
        << static_cast<const Property&>(suprasystem);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, SupraSystem &suprasystem)
{
    VersionID v = readHeader(ds, r_suprasystem);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> suprasystem.subsystems
            >> static_cast<Property&>(suprasystem);
    }
    else
        throw version_error(v, "1", r_suprasystem, CODELOC);
        
    return ds;
}

/** Null Constructor */
SupraSystem::SupraSystem() : ConcreteProperty<SupraSystem,Property>()
{}

/** Construct a supra-system consisting of 'n' sub systems */
SupraSystem::SupraSystem(int n)
            : ConcreteProperty<SupraSystem,Property>(),
              subsystems(n)
{
    subsystems.squeeze();
}

/** Construct a supra-system that consists of 'n' sub systems,
    which are initially created to be equal to 'system' */
SupraSystem::SupraSystem(const System &system, int n)
            : ConcreteProperty<SupraSystem,Property>(),
              subsystems(n)
{
    subsystems.squeeze();
    
    if (n > 0)
    {
        SupraSubSystem &subsys = subsystems[0].edit();
        
        subsys.setSubSystem(system);
        
        for (int i=1; i<n; ++i)
        {
            subsystems[i] = subsys;
        }
    }
}

/** Construct a supra-system that consists of the sub systems
    that are copies of those in 'systems' */
SupraSystem::SupraSystem(const QVector<System> &systems)
            : ConcreteProperty<SupraSystem,Property>(),
              subsystems(systems.count())
{
    subsystems.squeeze();
    
    int nsystems = systems.count();
    
    for (int i=0; i<nsystems; ++i)
    {
        subsystems[i].edit().setSubSystem( systems.at(i) );
    }
}

/** Construct a SupraSystem that consists of 'n' copies of
    the passed subsystem */
SupraSystem::SupraSystem(const SupraSubSystem &subsystem, int n)
            : ConcreteProperty<SupraSystem,Property>(),
              subsystems(n)
{
    subsystems.squeeze();
    
    if (n > 0)
    {
        for (int i=0; i<n; ++i)
        {
            subsystems[i] = subsystem;
        }
    }
}


/** Copy constructor */
SupraSystem::SupraSystem(const SupraSystem &other)
            : ConcreteProperty<SupraSystem,Property>(other), 
              subsystems(other.subsystems)
{}

/** Destructor */
SupraSystem::~SupraSystem()
{}

/** Copy assignment operator */
SupraSystem& SupraSystem::operator=(const SupraSystem &other)
{
    Property::operator=(other);
    subsystems = other.subsystems;
    
    return *this;
}

/** Comparison operator */
bool SupraSystem::operator==(const SupraSystem &other) const
{
    return (this == &other) or
           (subsystems.constData() == other.subsystems.constData()) or
           (subsystems == other.subsystems);
}

/** Comparison operator */
bool SupraSystem::operator!=(const SupraSystem &other) const
{
    return not this->operator==(other);
}

Q_GLOBAL_STATIC( SupraSystem, supraSystem )

/** Return the global null SupraSystem */
const SupraSystem& SupraSystem::null()
{
    return *(supraSystem());
}

/** This function is called after each SupraMove to collect any statistics
    about the SupraSystem */
void SupraSystem::collectSupraStats()
{}

/** This function is called just before all of the sub-systems are packed */
void SupraSystem::_pre_pack()
{}

/** This function is called just after all of the sub-systems are packed */
void SupraSystem::_post_pack()
{}

/** This function is called just before all of the sub-systems are unpacked */
void SupraSystem::_pre_unpack()
{}

/** This function is called just after all of the sub-systems are unpacked */
void SupraSystem::_post_unpack()
{}

/** Internal function used to return the 'ith' sub-system. No bounds
    checking is performed in 'i' */
const SupraSubSystem& SupraSystem::_pvt_subSystem(int i) const
{
    if (i < 0 or i >= subsystems.count())
        throw SireError::program_bug( QObject::tr(
            "Should not access system %1 as count == %2")
                .arg(i).arg(subsystems.count()), CODELOC );

    return subsystems.constData()[i].read();
}

/** Internal function used to return the 'ith' sub-system. No bounds
    checking is performed in 'i' */
SupraSubSystem& SupraSystem::_pvt_subSystem(int i)
{
    if (i < 0 or i >= subsystems.count())
        throw SireError::program_bug( QObject::tr(
            "Should not access system %1 as count == %2")
                .arg(i).arg(subsystems.count()), CODELOC );

    return subsystems.data()[i].edit();
}

/** Return the ith sub-system in this supra-system

    \throw SireError::invalid_index
*/
const SupraSubSystem& SupraSystem::operator[](int i) const
{
    return this->_pvt_subSystem( Index(i).map( subsystems.count() ) );
}

/** Return the ith sub-system in this supra-system

    \throw SireError::invalid_index
*/
const SupraSubSystem& SupraSystem::at(int i) const
{
    return this->operator[](i);
}

/** Return whether or not this SupraSystem contains no subsystems */
bool SupraSystem::isEmpty() const
{
    return subsystems.isEmpty();
}

/** Return the number of sub-systems in this supra-system */
int SupraSystem::nSubSystems() const
{
    return subsystems.count();
}

/** Return the number of sub-systems in this supra-system */
int SupraSystem::count() const
{
    return this->nSubSystems();
}

/** Return the number of sub-systems in this supra-system */
int SupraSystem::size()
{
    return this->nSubSystems();
}

/** Clear the statistics that are collected between blocks of 
    sub-system moves */
void SupraSystem::clearStatistics()
{
    int n = subsystems.count();
    
    for (int i=0; i<n; ++i)
    {
        subsystems[i].edit().clearStatistics();
    }
}

/** Clear the statistics that are collected within blocks of
    sub-system moves */
void SupraSystem::clearSubStatistics()
{
    int n = subsystems.count();
    
    for (int i=0; i<n; ++i)
    {
        subsystems[i].edit().clearSubStatistics();
    }
}

/** Clear all statistics collected during the moves */
void SupraSystem::clearAllStatistics()
{
    int n = subsystems.count();
    
    for (int i=0; i<n; ++i)
    {
        subsystems[i].edit().clearAllStatistics();
    }
}

/** Tell all sub-systems that the next energy calculate must
    be performed from scratch - this is useful for debugging */
void SupraSystem::mustNowRecalculateFromScratch()
{
    int n = subsystems.count();
    
    for (int i=0; i<n; ++i)
    {
        subsystems[i].edit().mustNowRecalculateFromScratch();
    }
}

/** Return whether or not *all* sub-systems are packed */
bool SupraSystem::isPacked() const
{
    int n = subsystems.count();
    
    for (int i=0; i<n; ++i)
    {
        if (not subsystems.at(i)->isPacked())
            return false;
    }
    
    return true;
}

/** Return whether or not *all* sub-systems are packed to memory */
bool SupraSystem::isPackedToMemory() const
{
    int n = subsystems.count();
    
    for (int i=0; i<n; ++i)
    {
        if (not subsystems.at(i)->isPackedToMemory())
            return false;
    }
    
    return true;
}

/** Return whether or not *all* sub-systems are packed to disk */
bool SupraSystem::isPackedToDisk() const
{
    int n = subsystems.count();
    
    for (int i=0; i<n; ++i)
    {
        if (not subsystems.at(i)->isPackedToDisk())
            return false;
    }
    
    return true;
}

/** Return whether or not any of the sub-systems are packed */
bool SupraSystem::anyPacked() const
{
    int n = subsystems.count();
    
    for (int i=0; i<n; ++i)
    {
        if (subsystems.at(i)->isPacked())
            return true;
    }
    
    return false;
}

/** Return whether or not any of the sub-systems are packed to memory */
bool SupraSystem::anyPackedToMemory() const
{
    int n = subsystems.count();
    
    for (int i=0; i<n; ++i)
    {
        if (subsystems.at(i)->isPackedToMemory())
            return true;
    }
    
    return false;
}

/** Return whether or not any of the sub-systems are packed to disk */
bool SupraSystem::anyPackedToDisk() const
{
    int n = subsystems.count();
    
    for (int i=0; i<n; ++i)
    {
        if (subsystems.at(i)->isPackedToDisk())
            return true;
    }
    
    return false;
}

/** Pack all sub-systems */
void SupraSystem::pack()
{
    int n = subsystems.count();

    if (n == 0)
        return;

    bool being_packed = not this->isPacked();

    if (being_packed)
        this->_pre_pack();
    
    for (int i=0; i<n; ++i)
    {
        if (not subsystems.at(i)->isPacked())
        {
            SupraSubSystemPtr unpacked_system = subsystems.at(i);
            
            subsystems[i].edit().pack();
            
            //see if there are any copies of this system (so that
            //we can copy the already-packed version)
            for (int j=i+1; j<n; ++j)
            {
                if (unpacked_system.constData() == subsystems.at(j).constData())
                {
                    subsystems[j] = subsystems.at(i);
                }
            }
        }
    }
    
    if (being_packed)
        this->_post_pack();
}

/** Unpack all sub-systems */
void SupraSystem::unpack()
{
    int n = subsystems.count();
    
    if (n == 0)
        return;
    
    bool being_unpacked = this->anyPacked();
    
    if (being_unpacked)
        this->_pre_unpack();
    
    for (int i=0; i<n; ++i)
    {
        if (subsystems.at(i)->isPacked())
        {
            SupraSubSystemPtr packed_system = subsystems.at(i);
            subsystems[i].edit().unpack();
            
            for (int j=i+1; j<n; ++j)
            {
                if (packed_system.constData() == subsystems.at(j).constData())
                {
                    subsystems[j] = subsystems.at(i);
                }
            }
        }
    }
    
    if (being_unpacked)
        this->_post_unpack();
}

/** Pack all sub-systems to disk */
void SupraSystem::packToDisk()
{
    int n = subsystems.count();
    
    if (n == 0)
        return;
        
    bool being_packed = not this->isPacked();
    
    if (being_packed)
        this->_pre_pack();
    
    for (int i=0; i<n; ++i)
    {
        if (not subsystems.at(i)->isPackedToDisk())
        {
            SupraSubSystemPtr old_system = subsystems.at(i);

            subsystems[i].edit().packToDisk();
            
            for (int j=i+1; j<n; ++j)
            {
                if (old_system.constData() == subsystems.at(j).constData())
                {
                    subsystems[j] = subsystems.at(i);
                }
            }
        }
    }
    
    if (being_packed)
        this->_post_pack();
}

/** Pack all sub-systems to disk, in the directory 'tempdir'. Note
    that this will not move sub-systems that are already packed into
    a different directory */
void SupraSystem::packToDisk(const QString &tempdir)
{
    int n = subsystems.count();
    
    if (n == 0)
        return;
        
    bool being_packed = not this->isPacked();
    
    if (being_packed)
        this->_pre_pack();
    
    for (int i=0; i<n; ++i)
    {
        if (not subsystems.at(i)->isPackedToDisk())
        {
            SupraSubSystemPtr old_system = subsystems.at(i);

            subsystems[i].edit().packToDisk(tempdir);
            
            for (int j=i+1; j<n; ++j)
            {
                if (old_system.constData() == subsystems.at(j).constData())
                {
                    subsystems[j] = subsystems.at(i);
                }
            }
        }
    }
    
    if (being_packed)
        this->_post_pack();
}

/** Pack all sub-systems to memory */
void SupraSystem::packToMemory()
{
    int n = subsystems.count();
    
    if (n == 0)
        return;
    
    bool being_packed = not this->isPacked();
    
    if (being_packed)
        this->_pre_pack();
    
    for (int i=0; i<n; ++i)
    {
        if (not subsystems.at(i)->isPackedToMemory())
        {
            SupraSubSystemPtr old_system = subsystems.at(i);
            subsystems[i].edit().packToMemory();
            
            for (int j=i+1; j<n; ++j)
            {
                if (old_system.constData() == subsystems.at(j).constData())
                {
                    subsystems[j] = subsystems.at(i);
                }
            }
        }
    }
    
    if (being_packed)
        this->_post_pack();
}

/** Pack the ith sub-system 

    \throw SireError::invalid_index
*/
void SupraSystem::pack(int i)
{
    i = Index(i).map( subsystems.count() );

    if (subsystems.at(i)->isPacked())
        return;
    
    this->_pre_pack();
    subsystems[i].edit().pack();
    this->_post_pack();
}

/** Unpack the ith sub-system 

    \throw SireError::invalid_index
*/
void SupraSystem::unpack(int i)
{
    i = Index(i).map( subsystems.count() );

    if (not subsystems.at(i)->isPacked())
        return;
    
    this->_pre_unpack();
    subsystems[i].edit().unpack();
    this->_post_unpack();
}

/** Pack the ith sub-system to disk

    \throw SireError::invalid_index
*/
void SupraSystem::packToDisk(int i)
{
    i = Index(i).map( subsystems.count() );

    if (subsystems.at(i)->isPackedToDisk())
        return;
    
    bool being_packed = not subsystems.at(i)->isPacked();
    
    if (being_packed)
        this->_pre_pack();
    
    subsystems[i].edit().packToDisk();
    
    if (being_packed)
        this->_post_pack();
}

/** Pack the ith sub-system to disk in the directory 'tempdir'

    \throw SireError::invalid_index
*/
void SupraSystem::packToDisk(int i, const QString &tempdir)
{
    i = Index(i).map( subsystems.count() );

    if (subsystems.at(i)->isPackedToDisk())
        return;
    
    bool being_packed = not subsystems.at(i)->isPacked();
    
    if (being_packed)
        this->_pre_pack();
        
    subsystems[i].edit().packToDisk(tempdir);

    if (being_packed)
        this->_post_pack();
}

/** Pack the ith sub-system to memory

    \throw SireError::invalid_index
*/
void SupraSystem::packToMemory(int i)
{
    i = Index(i).map( subsystems.count() );

    if (subsystems.at(i)->isPackedToMemory())
        return;
    
    bool being_packed = not subsystems.at(i)->isPacked();
    
    if (being_packed)
        this->_pre_pack();
        
    subsystems[i].edit().packToMemory();

    if (being_packed)
        this->_post_pack();
}

/** Internal function used to update all copies of 'old_subsystem' 
    so that they are equal to 'new_subsystem'. Use this function
    when you are performing an operation on all sub-systems and
    you want to ensure that any shared sub-systems are 
    still shared */
void SupraSystem::updateSubSystems(const SupraSubSystem *old_subsystem,
                                   const SupraSubSystem *new_subsystem,
                                   QSet<int> &done_systems)
{
    if (old_subsystem == 0 or new_subsystem == 0 or
        old_subsystem == new_subsystem)
    {
        return;
    }
        
    int n = subsystems.count();
    
    for (int i=0; i<n; ++i)
    {
        if (not done_systems.contains(i))
        {
            if ( &(subsystems.at(i).read()) == old_subsystem )
            {   
                subsystems[i] = *(new_subsystem);
                done_systems.insert(i);
            }
        }
    }
}

/** Set the ith sub-system equal to 'subsystem'

    \throw SireError::invalid_index
*/
void SupraSystem::setSubSystem(int i, const SupraSubSystem &subsystem)
{
    subsystems[ Index(i).map(subsystems.count()) ] = subsystem;
}

/** Set all sub-systems equal to 'subsystem' */
void SupraSystem::setSubSystem(const SupraSubSystem &subsystem)
{
    SupraSystemPtr old_state = this->clone();
    
    try
    {
        int n = subsystems.count();
        
        if (n == 0)
            return;
        
        QSet<int> done_systems;
        done_systems.reserve(n);
    
        //do the first subsystem
        this->setSubSystem(0, subsystem);
        
        //now set the rest from the copy of the first
        const SupraSubSystem &subsystem0 = *(subsystems.at(0));
    
        done_systems.insert(0);
    
        for (int i=1; i<n; ++i)
        {
            if ( not done_systems.contains(i) )
            {
                const SupraSubSystem *old_subsystem = subsystems.at(i).constData();
                this->setSubSystem(i, subsystem0);
                
                this->updateSubSystems(old_subsystem, subsystems.at(i).constData(), 
                                       done_systems);
            }
        }
    }
    catch(...)
    {
        this->copy(*old_state);
        throw;
    }
}

/** Set all sub-systems equal to a copy of those in 'system'. Note
    that both this SupraSystem and 'other' must have the same
    number of sub-systems 
    
    \throw SireError::incompatible_error
*/
void SupraSystem::setSubSystems(const SupraSystem &system)
{
    int n = subsystems.count();
    
    if (system.nSubSystems() != n)
    {
        throw SireError::incompatible_error( QObject::tr(
            "You cannot set the sub-systems for this supra-system from the "
            "passed supra-system, as they have different numbers of "
            "sub-systems (%1 vs. %2)")
                .arg(n).arg(system.nSubSystems()), CODELOC );
    }

    SupraSystemPtr old_state = this->clone();
    
    try
    {
        for (int i=0; i<n; ++i)
        {
            this->setSubSystem(i, system[i]);
        }
    }
    catch(...)
    {
        this->copy( *old_state );
        throw;
    }
}

/** Set the SireSystem::System used by the ith sub-system equal to 'system'

    \throw SireError::invalid_index
*/
void SupraSystem::setSubSystem(int i, const System &system)
{
    subsystems[ Index(i).map(subsystems.count()) ].edit().setSubSystem(system);
}

/** Set the SireSystem::System used by all sub-systems equal to 'system' */
void SupraSystem::setSubSystem(const System &system)
{
    SupraSystemPtr old_state = this->clone();
    
    try
    {
        int n = subsystems.count();
        QSet<int> done_systems;
        done_systems.reserve(n);
    
        for (int i=0; i<n; ++i)
        {
            if (not done_systems.contains(i))
            {
                const SupraSubSystem *old_subsystem = subsystems.at(i).constData();
            
                this->setSubSystem(i, system);
                
                this->updateSubSystems(old_subsystem, subsystems.at(i).constData(),
                                       done_systems);
            }
        }
    }
    catch(...)
    {
        this->copy( *old_state );
        throw;
    }
}

/** Set the moves that will be applied to the ith sub-system during
    each sub-move
    
    \throw SireError::invalid_index
*/
void SupraSystem::setSubMoves(int i, const Moves &moves)
{
    subsystems[ Index(i).map(subsystems.count()) ].edit().setSubMoves(moves);
}

/** Set the moves that will be applied to all sub-systems during
    each sub-move
    
    \throw SireError::invalid_index
*/
void SupraSystem::setSubMoves(const Moves &moves)
{
    SupraSystemPtr old_state = this->clone();
    
    try
    {
        int n = subsystems.count();
        QSet<int> done_systems;
        done_systems.reserve(n);
    
        for (int i=0; i<n; ++i)
        {
            if ( not done_systems.contains(i) )
            {
                const SupraSubSystem *old_subsystem = subsystems.at(i).constData();
            
                this->setSubMoves(i, moves);
                
                this->updateSubSystems(old_subsystem, subsystems.at(i).constData(),
                                       done_systems);
            }
        }
    }
    catch(...)
    {
        this->copy( *old_state );
        throw;
    }
}

/** Set the system and moves used for the ith sub-system to those
    contained in 'simstore'
    
    \throw SireError::invalid_index
*/
void SupraSystem::setSubSystemAndMoves(int i, const SimStore &simstore)
{
    subsystems[ Index(i).map(subsystems.count()) ].edit().setSubSystemAndMoves(simstore);
}

/** Set the system and moves used for all sub-systems to those
    contained in 'simstore' */
void SupraSystem::setSubSystemAndMoves(const SimStore &simstore)
{
    SupraSystemPtr old_state = this->clone();
    
    try
    {
        int n = subsystems.count();
        QSet<int> done_systems;
        done_systems.reserve(n);
    
        for (int i=0; i<n; ++i)
        {
            if (not done_systems.contains(i))
            {
                const SupraSubSystem *old_subsystem = subsystems.at(i).constData();
            
                this->setSubSystemAndMoves(i, simstore);
                
                this->updateSubSystems(old_subsystem, subsystems.at(i).constData(),
                                       done_systems);
            }
        }
    }
    catch(...)
    {
        this->copy(*old_state);
        throw;
    }
}

/** Set the system and moves that will be applied to the ith sub-system
    to 'system' and 'moves' 
    
    \throw SireError::invalid_index
*/
void SupraSystem::setSubSystemAndMoves(int i, const System &system, const Moves &moves)
{
    this->setSubSystemAndMoves( i, SimStore(system,moves) );
}

/** Set the system and moves used by all sub-systems to 
    'system' and 'moves' */
void SupraSystem::setSubSystemAndMoves(const System &system, const Moves &moves)
{
    this->setSubSystemAndMoves( SimStore(system,moves) );
}

/** Set the monitors used for all sub-systems to 'monitors', and set
    them all to update with a frequency of 'frequency' */
void SupraSystem::setSubMonitors(const SystemMonitors &monitors, int frequency)
{
    QVector<SupraSubSystemPtr> old_state = subsystems;
    
    try
    {
        SystemMonitors new_monitors = monitors;
        new_monitors.setAllFrequency(frequency);
    
        int n = subsystems.count();
    
        QSet<int> done_systems;
        done_systems.reserve(n);
    
        for (int i=0; i<n; ++i)
        {
            if (not done_systems.contains(i))
            {
                const SupraSubSystem *old_subsystem = subsystems.at(i).constData();
                subsystems[i].edit().setMonitors(monitors);
                this->updateSubSystems(old_subsystem, subsystems.at(i).constData(),
                                       done_systems);
            }
        }
    }
    catch(...)
    {
        subsystems = old_state;
        throw;
    }
}

/** Set the monitors used for the ith sub-system to 'monitors', and set
    them to update with a frequency of 'frequency' 
    
    \throw SireError::invalid_index
*/
void SupraSystem::setSubMonitors(int i, const SystemMonitors &monitors, int frequency)
{
    SystemMonitors new_monitors = monitors;
    new_monitors.setAllFrequency(frequency);
    
    subsystems[ Index(i).map(subsystems.count()) ].edit().setMonitors(monitors);
}

/** Add the monitor 'monitor' to all of the sub-systems, with the name 'name', 
    set to update with a frequency 'frequency' */
void SupraSystem::add(const QString &name, const SystemMonitor &monitor, int frequency)
{
    QVector<SupraSubSystemPtr> old_state = subsystems;
    
    try
    {
        int n = subsystems.count();
        
        QSet<int> done_systems;
        done_systems.reserve(n);
        
        for (int i=0; i<n; ++i)
        {
            if (not done_systems.contains(i))
            {
                const SupraSubSystem *old_subsystem = subsystems.at(i).constData();
            
                SupraSubSystem &system = subsystems[i].edit();
            
                SystemMonitors monitors = system.monitors();
            
                monitors.add(name, monitor, frequency);
            
                system.setMonitors(monitors);
                
                this->updateSubSystems(old_subsystem, subsystems.at(i).constData(),
                                       done_systems);
            }
        }
    }
    catch(...)
    {
        subsystems = old_state;
        throw;
    }
}

/** Add the monitor 'monitor' to the ith sub-system, with the name 'name',
    set to update with a frequency of 'frequency' 
    
    \throw SireError::invalid_index
*/
void SupraSystem::add(int i, const QString &name, const SystemMonitor &monitor,
                      int frequency)
{
    SupraSubSystem &system = subsystems[ Index(i).map(subsystems.count()) ].edit();
    
    SystemMonitors monitors = system.monitors();
    
    monitors.add(name, monitor, frequency);
    
    system.setMonitors(monitors);
}

/** Add the monitors 'monitors' to all of the sub-systems, set to update
    with a frequency of 'frequency' */             
void SupraSystem::add(const SystemMonitors &monitors, int frequency)
{
    QVector<SupraSubSystemPtr> old_state = subsystems;
    
    try
    {
        int n = subsystems.count();
        
        QSet<int> done_systems;
        done_systems.reserve(n);
        
        for (int i=0; i<n; ++i)
        {
            if (not done_systems.contains(i))
            {
                const SupraSubSystem *old_subsystem = subsystems.at(i).constData();
            
                SupraSubSystem &system = subsystems[i].edit();
            
                SystemMonitors sysmonitors = system.monitors();

                sysmonitors.add(monitors, frequency);
            
                system.setMonitors(sysmonitors);
                
                this->updateSubSystems(old_subsystem, subsystems.at(i).constData(),
                                       done_systems);
            }
        }
    }
    catch(...)
    {
        subsystems = old_state;
        throw;
    }
}

/** Add the monitors 'monitors' to the ith system, set to update with
    a frequency of 'frequency'
    
    \throw SireError::invalid_index
*/
void SupraSystem::add(int i, const SystemMonitors &monitors, int frequency)
{
    SupraSubSystem &system = subsystems[ Index(i).map(subsystems.count()) ].edit();
    
    SystemMonitors sysmonitors = system.monitors();
    
    sysmonitors.add(monitors, frequency);
    
    system.setMonitors(sysmonitors);
}

/** Set the number of moves to perform per block of sub-moves for 
    the ith sub-system
    
    \throw SireError::invalid_index
*/
void SupraSystem::setNSubMoves(int i, int nmoves)
{
    subsystems[ Index(i).map(subsystems.count()) ].edit().setNSubMoves(nmoves);
}

/** Set the number of moves to perform per block of sub-moves for
    every sub-system in this supra-system */
void SupraSystem::setNSubMoves(int nmoves)
{
    int n = subsystems.count();
    
    QSet<int> done_systems;
    done_systems.reserve(n);
    
    for (int i=0; i<n; ++i)
    {
        if (not done_systems.contains(i))
        {
            const SupraSubSystem *old_subsystem = subsystems.at(i).constData();
            this->setNSubMoves(i, nmoves);
            this->updateSubSystems(old_subsystem, subsystems.at(i).constData(),
                                   done_systems);
        }
    }
}

/** Set whether or not to record statistics between blocks of sub-moves
    at the ith sub-system
    
    \throw SireError::invalid_index
*/
void SupraSystem::setRecordStatistics(int i, bool record_stats)
{
    subsystems[ Index(i).map(subsystems.count()) ]
            .edit().setRecordStatistics(record_stats);
}

/** Set whether or not to record statistics between blocks of sub-moves
    for all sub-systems */
void SupraSystem::setRecordStatistics(bool record_stats)
{
    int n = subsystems.count();
    
    QSet<int> done_systems;
    done_systems.reserve(n);
    
    for (int i=0; i<n; ++i)
    {
        if (not done_systems.contains(i))
        {
            const SupraSubSystem *old_subsystem = subsystems.at(i).constData();
            this->setRecordStatistics(i, record_stats);
            this->updateSubSystems(old_subsystem, subsystems.at(i).constData(),
                                   done_systems);
        }
    }
}

/** Set whether or not to record statistics within blocks of sub-moves
    for the ith sub-system
    
    \throw SireError::invalid_index
*/
void SupraSystem::setRecordSubStatistics(int i, bool record_stats)
{
    subsystems[ Index(i).map(subsystems.count()) ]
            .edit().setRecordSubStatistics(record_stats);
}

/** Set whether or not to record statistics within blocks of sub-moves
    for all sub-systems */
void SupraSystem::setRecordSubStatistics(bool record_stats)
{
    int n = subsystems.count();
    
    QSet<int> done_systems;
    done_systems.reserve(n);
    
    for (int i=0; i<n; ++i)
    {
        if (not done_systems.contains(i))
        {
            const SupraSubSystem *old_subsystem = subsystems.at(i).constData();
            
            this->setRecordSubStatistics(i, record_stats);
            
            this->updateSubSystems(old_subsystem, subsystems.at(i).constData(),
                                   done_systems);
        }
    }
}

/** Set whether or not to record all statistics (both within and
    between sub-moves) for the ith sub-system
    
    \throw SireError::invalid_index
*/
void SupraSystem::setRecordAllStatistics(int i, bool record_stats)
{
    this->setRecordStatistics(i, record_stats);
    this->setRecordSubStatistics(i, record_stats);
}

/** Set whether or not to record all statistics (both within and
    between sub-moves) for all sub-systems */
void SupraSystem::setRecordAllStatistics(bool record_stats)
{
    int n = subsystems.count();
    
    QSet<int> done_systems;
    done_systems.reserve(n);
    
    for (int i=0; i<n; ++i)
    {
        if (not done_systems.contains(i))
        {
            const SupraSubSystem *old_subsystem = subsystems.at(i).constData();
            this->setRecordStatistics(i, record_stats);
            this->setRecordSubStatistics(i, record_stats);
            this->updateSubSystems( old_subsystem, subsystems.at(i).constData(),
                                    done_systems );
        }
    }
}

const char* SupraSystem::typeName()
{
    return QMetaType::typeName( qMetaTypeId<SupraSystem>() );
}

SupraSystem* SupraSystem::clone() const
{
    return new SupraSystem(*this);
}
