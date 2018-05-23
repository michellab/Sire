/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2018  Christopher Woods
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

#include "idengine.h"

#include "tostring.h"

#include <QRegExp>

using namespace SireMol;
using namespace parser_idengine;

////////
//////// Implementation of the IDNameEngine
////////

IDNameEngine::IDNameEngine(IDObject o, NameValues v)
             : SelectEngine(), obj(o), vals(v)
{}

IDNameEngine::~IDNameEngine()
{}

SelectResult selectAtom(const NameValues &vals, const SelectResult &mols)
{
    SelectResult::Container result;
    
    return result;
}

SelectResult selectCutGroup(const NameValues &vals, const SelectResult &mols)
{
    SelectResult::Container result;
    
    return result;
}

SelectResult selectResidue(const NameValues &vals, const SelectResult &mols)
{
    SelectResult::Container result;
    
    return result;
}

SelectResult selectChain(const NameValues &vals, const SelectResult &mols)
{
    SelectResult::Container result;
    
    return result;
}

SelectResult selectSegment(const NameValues &vals, const SelectResult &mols)
{
    SelectResult::Container result;
    
    return result;
}

SelectResult selectMolecule(const NameValues &vals, const SelectResult &mols)
{
    SelectResult::Container result;
    
    return result;
}

SelectResult IDNameEngine::select(const SelectResult &mols, const PropertyMap&) const
{
    SelectResult::Container result;
    
    qDebug() << CODELOC;
    
    switch(obj)
    {
    case ATOM:
        return selectAtom(vals, mols);
    case CUTGROUP:
        return selectCutGroup(vals, mols);
    case RESIDUE:
        return selectResidue(vals, mols);
    case CHAIN:
        return selectChain(vals, mols);
    case SEGMENT:
        return selectSegment(vals, mols);
    case MOLECULE:
        return selectMolecule(vals, mols);
    default:
        return SelectResult();
    }
}

////////
//////// Implementation of the IDNumberEngine
////////

IDNumberEngine::IDNumberEngine( IDObject o, RangeValues v )
               : obj(o), vals(v)
{}

IDNumberEngine::~IDNumberEngine()
{}

SelectResult IDNumberEngine::select(const SelectResult &mols, const PropertyMap &map) const
{
    SelectResult::Container result;
    
    return result;
}

////////
//////// Implementation of the IDIndexEngine
////////

IDIndexEngine::IDIndexEngine( IDObject o, RangeValues v )
              : obj(o), vals(v)
{}

IDIndexEngine::~IDIndexEngine()
{}

SelectResult IDIndexEngine::select(const SelectResult &mols, const PropertyMap &map) const
{
    SelectResult::Container result;
    
    return result;
}

////////
//////// Implementation of the IDAndEngine
////////

IDAndEngine::IDAndEngine( const SelectEnginePtr &p0, const SelectEnginePtr &p1 )
            : part0(p0), part1(p1)
{}

IDAndEngine::~IDAndEngine()
{}

SelectResult IDAndEngine::select(const SelectResult &mols, const PropertyMap &map) const
{
    SelectResult::Container result;
    
    return result;
}

////////
//////// Implementation of the IDOrEngine
////////

IDOrEngine::IDOrEngine( const SelectEnginePtr &part0, const SelectEnginePtr &part1 )
{
    parts.append(part0);
    parts.append(part1);
}

IDOrEngine::IDOrEngine( const QList<SelectEnginePtr> &p ) : parts(p)
{}

IDOrEngine::~IDOrEngine()
{}

SelectResult IDOrEngine::select(const SelectResult &mols, const PropertyMap &map) const
{
    SelectResult::Container result;
    
    return result;
}

////////
//////// Implementation of the IDNotEngine
////////

IDNotEngine::IDNotEngine( const SelectEnginePtr &p ) : part(p)
{}

IDNotEngine::~IDNotEngine()
{}

SelectResult IDNotEngine::select(const SelectResult &mols, const PropertyMap &map) const
{
    SelectResult::Container result;
    
    return result;
}

////////
//////// Implementation of the IDJoinEngine
////////

IDJoinEngine::IDJoinEngine( const SelectEnginePtr &p ) : part(p)
{}

IDJoinEngine::~IDJoinEngine()
{}

SelectResult IDJoinEngine::select(const SelectResult &mols, const PropertyMap &map) const
{
    SelectResult::Container result;
    
    return result;
}

////////
//////// Implementation of the IDSubScriptEngine
////////

IDSubScriptEngine::IDSubScriptEngine( const SelectEnginePtr &p, const RangeValue &v )
                  : part(p), val(v)
{}

IDSubScriptEngine::~IDSubScriptEngine()
{}

SelectResult IDSubScriptEngine::select(const SelectResult &mols, const PropertyMap &map) const
{
    SelectResult::Container result;
    
    return result;
}

////////
//////// Implementation of the IDWithEngine
////////

IDWithEngine::IDWithEngine( IDObject o, IDToken t, const SelectEnginePtr &p )
             : obj(o), token(t), part(p)
{}

IDWithEngine::~IDWithEngine()
{}

SelectResult IDWithEngine::select(const SelectResult &mols, const PropertyMap &map) const
{
    SelectResult::Container result;
    
    return result;
}
