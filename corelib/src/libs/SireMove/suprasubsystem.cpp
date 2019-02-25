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

#include "suprasubsystem.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMove;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<SupraSubSystem> r_suprasubsystem;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                        const SupraSubSystem &suprasubsystem)
{
    writeHeader(ds, r_suprasubsystem, 1);
    
    SharedDataStream sds(ds);
    
    sds << suprasubsystem.simstore 
        << suprasubsystem.sys_monitors
        << suprasubsystem.nsubmoves
        << suprasubsystem.record_stats
        << suprasubsystem.record_sub_stats
        << suprasubsystem.recalc_next_from_scratch
        << suprasubsystem.clear_subsys_stats
        << static_cast<const Property&>(suprasubsystem);
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds,
                                        SupraSubSystem &suprasubsystem)
{
    VersionID v = readHeader(ds, r_suprasubsystem);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> suprasubsystem.simstore 
            >> suprasubsystem.sys_monitors
            >> suprasubsystem.nsubmoves
            >> suprasubsystem.record_stats
            >> suprasubsystem.record_sub_stats
            >> suprasubsystem.recalc_next_from_scratch
            >> suprasubsystem.clear_subsys_stats
            >> static_cast<Property&>(suprasubsystem);
    }
    else
        throw version_error(v, "1", r_suprasubsystem, CODELOC);
        
    return ds;
}

/** Constructor */
SupraSubSystem::SupraSubSystem() 
               : ConcreteProperty<SupraSubSystem,Property>(),
                 record_stats(true),
                 record_sub_stats(true),
                 recalc_next_from_scratch(false),
                 clear_subsys_stats(false)
{}

/** Copy constructor */
SupraSubSystem::SupraSubSystem(const SupraSubSystem &other)
               : ConcreteProperty<SupraSubSystem,Property>(other),
                 simstore(other.simstore), sys_monitors(other.sys_monitors),
                 nsubmoves(other.nsubmoves), 
                 record_stats(other.record_stats),
                 record_sub_stats(other.record_sub_stats),
                 recalc_next_from_scratch(other.recalc_next_from_scratch),
                 clear_subsys_stats(other.clear_subsys_stats)
{}

/** Destructor */
SupraSubSystem::~SupraSubSystem()
{}

/** Copy assignment operator */
SupraSubSystem& SupraSubSystem::operator=(const SupraSubSystem &other)
{
    Property::operator=(other);
    
    simstore = other.simstore;
    sys_monitors = other.sys_monitors;
    nsubmoves = other.nsubmoves;
    record_stats = other.record_stats;
    record_sub_stats = other.record_sub_stats;
    recalc_next_from_scratch = other.recalc_next_from_scratch;
    clear_subsys_stats = other.clear_subsys_stats;
    
    return *this;
}

/** Comparison operator */
bool SupraSubSystem::operator==(const SupraSubSystem &other) const
{
    return (this == &other) or
           (nsubmoves == other.nsubmoves and 
            record_stats == other.record_stats and
            record_sub_stats == other.record_sub_stats and
            recalc_next_from_scratch == other.recalc_next_from_scratch and
            clear_subsys_stats == other.clear_subsys_stats and
            sys_monitors == other.sys_monitors and
            simstore == other.simstore);
}

/** Comparison operator */
bool SupraSubSystem::operator!=(const SupraSubSystem &other) const
{
    return not this->operator==(other);
}

Q_GLOBAL_STATIC( SupraSubSystem, supraSubSystem )

/** Return the global null SupraSubSystem */
const SupraSubSystem& SupraSubSystem::null()
{
    return *(supraSubSystem());
}

/** Return the monitors that are applied to this sub-system. These
    are the monitors that are applied between blocks of sub-moves.
    e.g.
    
    subsys.monitors()   // Monitors applied between blocks of sub-moves.
                        // These monitors stay with the SupraSubSystem, and
                        // are collected using "subsys.collectStats()" and
                        // are cleared using "subsys.clearStatistics()"
                        // These are collected if record_stats is true
                        
    subsys.system().monitors()   // Monitors applied with each block of sub-moves.
                                 // These monitors stay with the System within each
                                 // SupraSubSystem, are collected within a block
                                 // of sub-system moves if record_sub_stats is true,
                                 // and are cleared using "subsys.clearSubStatistics()"
                                 // These are collected if record_sub_stats is true
*/
const SystemMonitors& SupraSubSystem::monitors() const
{
    return sys_monitors;
}

/** Return the system that is part of this sub-system 

    \throw SireError::invalid_state
*/
const System& SupraSubSystem::subSystem() const
{
    return simstore.system();
}

/** Return the moves that will be applied to the sub-system 

    \throw SireError::invalid_state
*/
const Moves& SupraSubSystem::subMoves() const
{
    return simstore.moves();
}

/** Return both the system and moves that are part of this sub-system */
const SimStore& SupraSubSystem::subSystemAndMoves() const
{
    return simstore;
}

/** Return the number of moves to be applied to this system */
int SupraSubSystem::nSubMoves() const
{
    return nsubmoves;
}

/** Return whether or not we are recording statistics 
    between blocks of sub-moves */
bool SupraSubSystem::recordingStatistics() const
{
    return record_stats;
}

/** Return whether or not we are recording statistics within 
    the sub-system */
bool SupraSubSystem::recordingSubStatistics() const
{
    return record_sub_stats;
}

/** Return whether or not the system is packed */
bool SupraSubSystem::isPacked() const
{
    return simstore.isPacked();
}

/** Return whether or not this system is packed to disk */
bool SupraSubSystem::isPackedToDisk() const
{
    return simstore.isPackedToDisk();
}

/** Return whether or not this system is packed to memory */
bool SupraSubSystem::isPackedToMemory() const
{
    return simstore.isPackedToMemory();
}

/** This function is called just before the system is packed */
void SupraSubSystem::_pre_pack()
{}

/** This function is called just after the system is packed */
void SupraSubSystem::_post_pack()
{}

/** This function is called just before the system is unpacked */
void SupraSubSystem::_pre_unpack()
{}

/** This function is called just after whenever the 
    simstore is unpacked */
void SupraSubSystem::_post_unpack()
{
    if (clear_subsys_stats or recalc_next_from_scratch)
    {
        System sys = simstore.system();
        
        if (clear_subsys_stats)
            sys.clearStatistics();
            
        if (recalc_next_from_scratch)
            sys.mustNowRecalculateFromScratch();
            
        simstore.setSystem(sys);
        
        clear_subsys_stats = false;
        recalc_next_from_scratch = false;
    }
}

/** Pack the system. This does nothing if the system is already packed */
void SupraSubSystem::pack()
{
    if (simstore.isPacked())
        return;
        
    this->_pre_pack();
    simstore.pack();
    this->_post_pack();
}

/** Unpack the system. This does nothing if the system is already unpacked */
void SupraSubSystem::unpack()
{
    if (not simstore.isPacked())
        return;
        
    this->_pre_unpack();
    simstore.unpack();
    this->_post_unpack();
}

/** Pack the system to disk. This does nothing if the system
    is already packed to disk */
void SupraSubSystem::packToDisk()
{
    if (simstore.isPackedToDisk())
        return;
   
    else if (not simstore.isPacked())
    {
        this->_pre_pack();
        simstore.packToDisk();
        this->_post_pack();
    }
    else
        simstore.packToDisk();
}

/** Pack the system to disk, into the directory 'tempdir'.
    This does nothing if the system is already packed to
    disk (even if the system is packed into a different
    directory) */
void SupraSubSystem::packToDisk(const QString &tempdir)
{
    if (simstore.isPackedToDisk())
        return;
        
    else if (not simstore.isPacked())
    {
        this->_pre_pack();
        simstore.packToDisk(tempdir);
        this->_post_pack();
    }
    else
        simstore.packToDisk(tempdir);
}

/** Pack the system to memory. This does nothing if the system
    is already packed to memory */
void SupraSubSystem::packToMemory()
{
    if (simstore.isPackedToMemory())
        return;
        
    else if (not simstore.isPacked())
    {
        this->_pre_pack();
        simstore.packToMemory();
        this->_post_pack();
    }
    else
        simstore.packToMemory();
}

/** Call this function to collect sub-system level statistics. This is
    called between blocks of sub-moves. */
void SupraSubSystem::collectStats()
{
    if (sys_monitors.isEmpty() or not record_stats)
        return;
        
    System system = simstore.system();
    
    sys_monitors.monitor(system);
    
    //copy the system back - this is because some of 
    //the monitors may change the system
    simstore.setSystem(system);
}

/** Perform the sub-system moves on this sub-system. This performs
    the moves in this->subMoves() on this->subSystem(), collecting
    statistics both during the block of sub-moves, and then after
    the block of sub-moves if 'recording_statistics' is true,
    and if 'sysmon.recordingSubStatistics()' and 
    'sysmon.recordingStatistics()' are true respectively.
    
    \throw SireError::invalid_state
*/
void SupraSubSystem::subMove(bool recording_statistics)
{
    if (nsubmoves <= 0)
        return;

    MovesPtr moves = simstore.moves();
    System system = simstore.system();
    
    system = moves.edit().move(system, nsubmoves, 
                               recording_statistics and record_sub_stats);
    
    simstore.setSystemAndMoves(system, moves);
    
    if (recording_statistics)
        this->collectStats();
}

/** Clear the SupraSubSystem level statistics (this clears the monitors
    that are applied at the end of blocks of sub-moves) */
void SupraSubSystem::clearStatistics()
{
    sys_monitors.clearStatistics();
}

/** Clear the system level statistics (this clears the monitors
    that are applied to the sub-system as it is being acted on
    by the sub-moves) */
void SupraSubSystem::clearSubStatistics()
{
    if (simstore.isPacked())
    {
        clear_subsys_stats = true;
    }
    else
    {
        System system = simstore.system();
        system.clearStatistics();
        simstore.setSystem(system);
    }
}

/** Completely clear all statistics from this system (this
    calls both clearStatistics and clearSubStatistics) */
void SupraSubSystem::clearAllStatistics()
{
    this->clearStatistics();
    this->clearSubStatistics();
}

/** Tell the system that the next energy to be calculated
    should be recalculated from scratch. This is mainly
    useful for debugging */
void SupraSubSystem::mustNowRecalculateFromScratch()
{
    if (simstore.isPacked())
    {
        recalc_next_from_scratch = true;
    }
    else
    {
        System system = simstore.system();
        system.mustNowRecalculateFromScratch();
        simstore.setSystem(system);
    }
}

/** Set the sub-system that will be acted on by the set of sub-moves */
void SupraSubSystem::setSubSystem(const System &subsystem)
{
    simstore.setSystem(subsystem);
    
    recalc_next_from_scratch = false;
    clear_subsys_stats = false;
}

/** Set the moves that will be used to move the sub-system in each
    block of sub-moves */
void SupraSubSystem::setSubMoves(const Moves &submoves)
{
    simstore.setMoves(submoves);
}

/** Set the sub-system and sub-moves together */
void SupraSubSystem::setSubSystemAndMoves(const SimStore &new_simstore)
{
    simstore = new_simstore;
    
    recalc_next_from_scratch = false;
    clear_subsys_stats = false;
}

/** Set the monitors that will be applied between blocks of sub-moves */
void SupraSubSystem::setMonitors(const SystemMonitors &monitors)
{
    sys_monitors = monitors;
}

/** Set the number of moves to apply to the sub-system within each
    block of sub-moves */
void SupraSubSystem::setNSubMoves(int n)
{
    if (n <= 0)
        nsubmoves = 0;
    else
        nsubmoves = n;
}

/** Set whether or not to record statistics within each block of
    sub-moves - this affects the monitors in subsys.system().monitors() */
void SupraSubSystem::setRecordSubStatistics(bool record_stats)
{
    record_sub_stats = record_stats;
}

/** Set whether or not to record statistics between blocks of sub-moves.
    This affects the subsys.monitors() */
void SupraSubSystem::setRecordStatistics(bool recording)
{
    record_stats = recording;
}

const char* SupraSubSystem::typeName()
{
    return QMetaType::typeName( qMetaTypeId<SupraSubSystem>() );
}

SupraSubSystem* SupraSubSystem::clone() const
{
    return new SupraSubSystem(*this);
}
