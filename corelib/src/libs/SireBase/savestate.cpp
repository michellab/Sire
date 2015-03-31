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

#include <QThreadStorage>
#include <QMutex>

#include "savestate.h"

#include "SireError/errors.h"

using namespace SireBase;

typedef QThreadStorage< QSet<const void*>* > StateRegistry;

static StateRegistry *state_registry(0);

Q_GLOBAL_STATIC( QMutex, globalMutex )

static QSet<const void*>& stateRegistry()
{
    if (state_registry == 0)
    {
        QMutexLocker lkr( globalMutex() );
        
        if (state_registry == 0)
            state_registry = new StateRegistry();
    }
    
    if (not state_registry->hasLocalData())
        state_registry->setLocalData( new QSet<const void*>() );

    return *(state_registry->localData());
}

/** Null constructor */
SaveState::SaveState() : current_state(0), forced_state(false)
{}

/** Internal constructor */
SaveState::SaveState(const Property &object, bool force)
          : old_state( object.clone() ), current_state(&object),
            forced_state(force)
{}

/** Copy constructor */
SaveState::SaveState(const SaveState &other)
          : old_state(other.old_state), current_state(other.current_state),
            forced_state(other.forced_state)
{}

/** Destructor */
SaveState::~SaveState()
{
    //if the state will be destroyed, then remove the
    //record from the thread-global dictionary
    if (current_state and not forced_state)
    {
        if (old_state.unique())
            ::stateRegistry().remove(current_state);
    }
}
    
/** Copy assignment operator */
SaveState& SaveState::operator=(const SaveState &other)
{
    if (this != &other)
    {
        if (current_state and not forced_state)
        {
            if (old_state.unique())
                ::stateRegistry().remove(current_state);
        }
    
        old_state = other.old_state;
        current_state = other.current_state;
        forced_state = other.forced_state;
    }
    
    return *this;
}
    
/** Restore the state of the object 'object' - this does nothing
    if this is a null state */
void SaveState::restore(Property &object)
{
    if (current_state)
    {
        if (&object != current_state)
            throw SireError::program_bug( QObject::tr(
                "An attempt is being made to restore the state of the wrong object! "
                "%1 vs. %2")
                    .arg(toInt(&object)).arg(toInt(current_state)), CODELOC );

        object.copy( *old_state );
    }
}

/** Return whether or not this is a null state */
bool SaveState::isNull()
{
    return current_state == 0;
}

/** Save the state of the object 'object' - this returns a null state
    if the state of this object has already been saved */
SaveState SaveState::save(const Property &object)
{
    QSet<const void*> &registry = ::stateRegistry();
    
    if (not registry.contains(&object))
    {
        registry.insert( &object );
        return SaveState(object);
    }
    else
        return SaveState();
}

/** Force the saving of the state of the object 'object' - this
    saves the state even if an existing state of this object has
    been saved */
SaveState SaveState::forceSave(const Property &object)
{
    QSet<const void*> &registry = ::stateRegistry();
    
    if (not registry.contains(&object))
    {
        registry.insert( &object );
        return SaveState(object);
    }
    else
        return SaveState(object, true);
}
