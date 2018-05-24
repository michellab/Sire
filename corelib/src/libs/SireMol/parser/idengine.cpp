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

IDNameEngine::IDNameEngine() : SelectEngine()
{}

SelectEnginePtr IDNameEngine::construct(IDObject o, NameValues vals)
{
    IDNameEngine *ptr = new IDNameEngine();
    ptr->obj = o;

    try
    {
        for (const auto val : vals)
        {
            if (val.value.which() == 0)
            {
                RegExpValue v = boost::get<RegExpValue>(val.value);
                QString r = QString::fromStdString(v.value);
                
                if (v.is_case_sensitive)
                    ptr->regexps.append( QRegExp(r, Qt::CaseSensitive) );
                else
                    ptr->regexps.append( QRegExp(r, Qt::CaseInsensitive) );
            }
            else if (val.value.which() == 1)
                ptr->names.append( QString::fromStdString(boost::get<std::string>(val.value)) );
        }
    }
    catch(...)
    {
        delete ptr;
        throw;
    }
    
    return makePtr(ptr);
}

IDNameEngine::~IDNameEngine()
{}

SelectResult IDNameEngine::selectAtoms(const SelectResult &mols) const
{
    SelectResult::Container result;
    
    for (const auto mol : mols)
    {
        auto selected_atoms = mol.selection();

        for (const auto name : names)
        {
            
        }
    }
    
    return result;
}

SelectResult IDNameEngine::selectCutGroups(const SelectResult &mols) const
{
    SelectResult::Container result;
    
    return result;
}

SelectResult IDNameEngine::selectResidues(const SelectResult &mols) const
{
    SelectResult::Container result;
    
    return result;
}

SelectResult IDNameEngine::selectChains(const SelectResult &mols) const
{
    SelectResult::Container result;
    
    return result;
}

SelectResult IDNameEngine::selectSegments(const SelectResult &mols) const
{
    SelectResult::Container result;
    
    return result;
}

SelectResult IDNameEngine::selectMolecules(const SelectResult &mols) const
{
    SelectResult::Container result;
    
    return result;
}

SelectResult IDNameEngine::select(const SelectResult &mols, const PropertyMap&) const
{
    SelectResult::Container result;
    
    switch(obj)
    {
    case AST::ATOM:
        return selectAtoms(mols);
    case AST::CUTGROUP:
        return selectCutGroups(mols);
    case AST::RESIDUE:
        return selectResidues(mols);
    case AST::CHAIN:
        return selectChains(mols);
    case AST::SEGMENT:
        return selectSegments(mols);
    case AST::MOLECULE:
        return selectMolecules(mols);
    default:
        return SelectResult();
    }
}

SelectEngine::ObjType IDNameEngine::objectType() const
{
    switch(obj)
    {
    case AST::ATOM:
        return SelectEngine::ATOM;
    case AST::CUTGROUP:
        return SelectEngine::CUTGROUP;
    case AST::RESIDUE:
        return SelectEngine::RESIDUE;
    case AST::CHAIN:
        return SelectEngine::CHAIN;
    case AST::SEGMENT:
        return SelectEngine::SEGMENT;
    case AST::MOLECULE:
        return SelectEngine::MOLECULE;
    default:
        return SelectEngine::COMPLEX;
    }
}

////////
//////// Implementation of the IDNumberEngine
////////

IDNumberEngine::IDNumberEngine()
{}

SelectEnginePtr IDNumberEngine::construct( IDObject o, RangeValues v )
{
    IDNumberEngine *ptr = new IDNumberEngine();
    ptr->obj = o;
    ptr->vals = v;

    return makePtr(ptr);
}

IDNumberEngine::~IDNumberEngine()
{}

SelectResult IDNumberEngine::select(const SelectResult &mols, const PropertyMap &map) const
{
    SelectResult::Container result;
    
    return result;
}

SelectEngine::ObjType IDNumberEngine::objectType() const
{
    switch(obj)
    {
    case AST::ATOM:
        return SelectEngine::ATOM;
    case AST::CUTGROUP:
        return SelectEngine::CUTGROUP;
    case AST::RESIDUE:
        return SelectEngine::RESIDUE;
    case AST::CHAIN:
        return SelectEngine::CHAIN;
    case AST::SEGMENT:
        return SelectEngine::SEGMENT;
    case AST::MOLECULE:
        return SelectEngine::MOLECULE;
    default:
        return SelectEngine::COMPLEX;
    }
}

////////
//////// Implementation of the IDIndexEngine
////////

IDIndexEngine::IDIndexEngine()
{}

SelectEnginePtr IDIndexEngine::construct( IDObject o, RangeValues v )
{
    IDIndexEngine *ptr = new IDIndexEngine();
    ptr->obj = o;
    ptr->vals = v;
    
    return makePtr(ptr);
}

IDIndexEngine::~IDIndexEngine()
{}

SelectResult IDIndexEngine::select(const SelectResult &mols, const PropertyMap &map) const
{
    SelectResult::Container result;
    
    return result;
}

SelectEngine::ObjType IDIndexEngine::objectType() const
{
    switch(obj)
    {
    case AST::ATOM:
        return SelectEngine::ATOM;
    case AST::CUTGROUP:
        return SelectEngine::CUTGROUP;
    case AST::RESIDUE:
        return SelectEngine::RESIDUE;
    case AST::CHAIN:
        return SelectEngine::CHAIN;
    case AST::SEGMENT:
        return SelectEngine::SEGMENT;
    case AST::MOLECULE:
        return SelectEngine::MOLECULE;
    default:
        return SelectEngine::COMPLEX;
    }
}

////////
//////// Implementation of the IDAndEngine
////////

IDAndEngine::IDAndEngine()
{}

SelectEnginePtr IDAndEngine::construct(SelectEnginePtr p0, SelectEnginePtr p1)
{
    IDAndEngine *ptr = new IDAndEngine();
    
    auto p = makePtr(ptr);
    
    if (p0)
        p0->setParent(p);

    if (p1)
        p1->setParent(p);

    ptr->part0 = p0;
    ptr->part1 = p1;
    
    return p;
}

IDAndEngine::~IDAndEngine()
{}

SelectResult IDAndEngine::select(const SelectResult &mols, const PropertyMap &map) const
{
    SelectResult::Container result;
    
    return result;
}

SelectEnginePtr IDAndEngine::simplify()
{
    if (part0.get())
        part0 = part0->simplify();
    
    if (part1.get())
        part1 = part1->simplify();
    
    if (part0.get() == 0)
        return part1;
    else if (part1.get() == 0)
        return part0;
    else
        return selfptr.lock();
}

SelectEngine::ObjType IDAndEngine::objectType() const
{
    if (part0.get() and part1.get())
    {
        auto o0 = part0->objectType();
        auto o1 = part1->objectType();
        
        if (o0 == o1)
            return o0;
        else
        {
            //the object type is always the smallest, e.g. atom and residue == atom
            return qMin(o0,o1);
        }
    }
    else if (part0.get())
    {
        return part0->objectType();
    }
    else if (part1.get())
    {
        return part1->objectType();
    }
    else
        return SelectEngine::COMPLEX;
}

////////
//////// Implementation of the IDOrEngine
////////

IDOrEngine::IDOrEngine()
{}

SelectEnginePtr IDOrEngine::construct(SelectEnginePtr part0, SelectEnginePtr part1)
{
    IDOrEngine *ptr = new IDOrEngine();
    
    auto p = makePtr(ptr);
    
    if (part0)
        part0->setParent(p);

    if (part1)
        part1->setParent(p);
    
    ptr->parts.append(part0);
    ptr->parts.append(part1);

    return p;
}

SelectEnginePtr IDOrEngine::construct(QList<SelectEnginePtr> parts)
{
    IDOrEngine *ptr = new IDOrEngine();
    
    auto p = makePtr(ptr);
    
    for (auto &part : parts)
    {
        if (part)
            part->setParent(p);
    }
    
    ptr->parts = parts;
    
    return p;
}

IDOrEngine::~IDOrEngine()
{}

SelectResult IDOrEngine::select(const SelectResult &mols, const PropertyMap &map) const
{
    SelectResult::Container result;
    
    return result;
}

SelectEnginePtr IDOrEngine::simplify()
{
    for (auto &part : parts)
    {
        part = part->simplify();
    }
    
    return selfptr.lock();
}

SelectEngine::ObjType IDOrEngine::objectType() const
{
    bool set = false;
    auto o = SelectEngine::COMPLEX;

    for (auto part : parts)
    {
        if (part.get())
        {
            if (not set)
            {
                o = part->objectType();
                set = true;
            }
            else if (o != part->objectType())
            {
                return SelectEngine::COMPLEX;
            }
        }
    }
    
    return o;
}

////////
//////// Implementation of the IDNotEngine
////////

IDNotEngine::IDNotEngine()
{}

SelectEnginePtr IDNotEngine::construct(SelectEnginePtr part)
{
    IDNotEngine *ptr = new IDNotEngine();
    auto p = makePtr(ptr);
    
    if (part)
        part->setParent(p);
    
    ptr->part = part;
    
    return p;
}

IDNotEngine::~IDNotEngine()
{}

SelectResult IDNotEngine::select(const SelectResult &mols, const PropertyMap &map) const
{
    SelectResult::Container result;
    
    return result;
}

SelectEnginePtr IDNotEngine::simplify()
{
    if (part.get())
        part = part->simplify();
    
    return selfptr.lock();
}

SelectEngine::ObjType IDNotEngine::objectType() const
{
    if (part)
        return part->objectType();
    else
        return SelectEngine::COMPLEX;
}

////////
//////// Implementation of the IDJoinEngine
////////

IDJoinEngine::IDJoinEngine()
{}

SelectEnginePtr IDJoinEngine::construct(SelectEnginePtr part)
{
    IDJoinEngine *ptr = new IDJoinEngine();
    auto p = makePtr(ptr);
    
    if (part)
        part->setParent(p);
    
    ptr->part = part;
    
    return p;
}

IDJoinEngine::~IDJoinEngine()
{}

SelectResult IDJoinEngine::select(const SelectResult &mols, const PropertyMap &map) const
{
    SelectResult::Container result;
    
    return result;
}

SelectEnginePtr IDJoinEngine::simplify()
{
    if (part.get())
        part = part->simplify();
    
    return selfptr.lock();
}

SelectEngine::ObjType IDJoinEngine::objectType() const
{
    return SelectEngine::COMPLEX;
}

////////
//////// Implementation of the IDSubScriptEngine
////////

IDSubScriptEngine::IDSubScriptEngine()
{}

SelectEnginePtr IDSubScriptEngine::construct(SelectEnginePtr part, const RangeValue &val)
{
    IDSubScriptEngine *ptr = new IDSubScriptEngine();
    auto p = makePtr(ptr);
    
    if (part)
        part->setParent(p);
    
    ptr->part = part;
    ptr->val = val;
    
    return p;
}

IDSubScriptEngine::~IDSubScriptEngine()
{}

SelectResult IDSubScriptEngine::select(const SelectResult &mols, const PropertyMap &map) const
{
    SelectResult::Container result;
    
    return result;
}

SelectEnginePtr IDSubScriptEngine::simplify()
{
    if (part.get())
        part = part->simplify();
    
    return selfptr.lock();
}

SelectEngine::ObjType IDSubScriptEngine::objectType() const
{
    if (part)
        return part->objectType();
    else
        return SelectEngine::COMPLEX;
}

////////
//////// Implementation of the IDWithEngine
////////

IDWithEngine::IDWithEngine()
{}

SelectEnginePtr IDWithEngine::construct( IDObject obj, IDToken token, SelectEnginePtr part)
{
    IDWithEngine *ptr = new IDWithEngine();
    auto p = makePtr(ptr);
    
    if (part)
        part->setParent(p);
    
    ptr->part = part;
    ptr->obj = obj;
    ptr->token = token;
    
    return p;
}

IDWithEngine::~IDWithEngine()
{}

SelectResult IDWithEngine::select(const SelectResult &mols, const PropertyMap &map) const
{
    SelectResult::Container result;
    
    return result;
}

SelectEnginePtr IDWithEngine::simplify()
{
    if (part.get())
        part = part->simplify();
    
    return selfptr.lock();
}

SelectEngine::ObjType IDWithEngine::objectType() const
{
    switch(obj)
    {
    case AST::ATOM:
        return SelectEngine::ATOM;
    case AST::CUTGROUP:
        return SelectEngine::CUTGROUP;
    case AST::RESIDUE:
        return SelectEngine::RESIDUE;
    case AST::CHAIN:
        return SelectEngine::CHAIN;
    case AST::SEGMENT:
        return SelectEngine::SEGMENT;
    case AST::MOLECULE:
        return SelectEngine::MOLECULE;
    default:
        return SelectEngine::COMPLEX;
    }
}
