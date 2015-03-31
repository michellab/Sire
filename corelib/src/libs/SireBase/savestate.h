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

#ifndef SIREBASE_SAVESTATE_H
#define SIREBASE_SAVESTATE_H

#include "property.h"

SIRE_BEGIN_HEADER

namespace SireBase
{

/** This class is used to save the state of an object derived
    from SireBase::Property. This is used
    when the state of an object is repeatedly saved in nested
    function calls (to protect against an exception being thrown)
    but only the state in the top-level function call will actually
    be restored if there is a problem (so the other saved states
    just waste time and resources)
  
    e.g. A function that changes 10 things may save state before any
         of the 10 are changed, then call functions that change each
         of the 10 things individually. Each of the individual functions
         may also save state, but this is unnecessary as the 10-change function
         is the one that will restore the state if something goes wrong.
         The 10 saved states from the individual-change functions are thus
         unnecessary, and waste copying time and memory. However, those
         10 individual-change functions still need to save state, as 
         they may be called individually by the user. By using this class
         to save state, we can use the same code to save state in both
         the 10-change and the individual-change functions, and the
         right thing is done in both cases automatically (rather than
         us having to write an individual-change function that 
         saves state (to be called from top-level) then a separate
         individual-change function that doesn't save state (to be 
         called by the 10-change function)

    Usage:
    
    \code
    class Foo : public ConcreteProperty<Foo,Property>
    {
    public:
        Foo();
        
        void bar1();
        void bar2();
        void allBar();
    };
    
    void Foo::bar1()
    {
        SaveState old_state = SaveState::save(*this);
        
        try
        {
            //do stuff
        }
        catch(...)
        {
            old_state.restore(*this);
            throw;
        }
    }
    
    void Foo::bar2()
    {
        SaveState old_state = SaveState::save(*this);
        
        try
        {
            //do stuff
        }
        catch(...)
        {
            old_state.restore(*this);
            throw;
        }
    }
        
    void Foo::allBar()
    {
        SaveState old_state = SaveState::save(*this);
        
        try
        {
            this->bar1();   // state now won't be saved in bar1()
            this->bar2();   // state now won't be saved in bar2()
        }
        catch(...)
        {
            old_state.restore(*this);
            throw;
        }
    }
        
    \endcode    
    
    @author Christopher Woods
*/
class SIREBASE_EXPORT SaveState
{
public:
    SaveState();
    SaveState(const SaveState &other);
    ~SaveState();
        
    SaveState& operator=(const SaveState &other);
        
    void restore(Property &object);

    bool isNull();

    static SaveState save(const Property &object);
    static SaveState forceSave(const Property &object);

private:
    SaveState(const Property &object, bool force=false);

    /** Pointer to the old state of the object */
    PropertyPtr old_state;

    /** Pointer to the current state of the object */
    const Property *current_state;

    /** Whether or not this is a forced state */
    bool forced_state;
};

}

SIRE_END_HEADER

#endif
