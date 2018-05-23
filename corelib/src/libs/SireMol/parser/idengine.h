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

#ifndef SIREMOL_PARSER_IDENGINE_H
#define SIREMOL_PARSER_IDENGINE_H

#include "SireMol/select.h"

#include "ast.h"

SIRE_BEGIN_HEADER

namespace parser_idengine
{

using SireMol::parser::SelectEngine;

using SireMol::SelectResult;
using SireBase::PropertyMap;

using namespace AST;

/** Internal class providing the SelectEngine for objects 
    based on their names
    
    @author Christopher Woods
*/
class IDNameEngine : public SelectEngine
{
public:
    IDNameEngine( IDObject obj, NameValues vals );
    ~IDNameEngine();
    
protected:
    SelectResult select(const SelectResult &mols, const PropertyMap &map) const;
    
private:
    IDObject obj;
    NameValues vals;
};

/** Internal class providing the SelectEngine for objects 
    based on their numbers
    
    @author Christopher Woods
*/
class IDNumberEngine : public SelectEngine
{
public:
    IDNumberEngine( IDObject obj, RangeValues vals );
    ~IDNumberEngine();
    
protected:
    SelectResult select(const SelectResult &mols, const PropertyMap &map) const;
    
private:
    IDObject obj;
    RangeValues vals;
};

/** Internal class providing the SelectEngine for objects 
    based on their indicies (index)
    
    @author Christopher Woods
*/
class IDIndexEngine : public SelectEngine
{
public:
    IDIndexEngine( IDObject obj, RangeValues vals );
    ~IDIndexEngine();
    
protected:
    SelectResult select(const SelectResult &mols, const PropertyMap &map) const;
    
private:
    IDObject obj;
    RangeValues vals;
};

/** Internal class providing the SelectEngine for objects 
    in an "and" expression
    
    @author Christopher Woods
*/
class IDAndEngine : public SelectEngine
{
public:
    IDAndEngine( const SelectEnginePtr &part0, const SelectEnginePtr &part1 );
    ~IDAndEngine();
    
protected:
    SelectResult select(const SelectResult &mols, const PropertyMap &map) const;
    
private:
    SelectEnginePtr part0, part1;
};

/** Internal class providing the SelectEngine for objects 
    in an "or" expression
    
    @author Christopher Woods
*/
class IDOrEngine : public SelectEngine
{
public:
    IDOrEngine( const SelectEnginePtr &part0, const SelectEnginePtr &part1 );
    IDOrEngine( const QList<SelectEnginePtr> &parts );
    ~IDOrEngine();
    
protected:
    SelectResult select(const SelectResult &mols, const PropertyMap &map) const;
    
private:
    QList<SelectEnginePtr> parts;
};

/** Internal class providing the SelectEngine for objects 
    in a "not" expression
    
    @author Christopher Woods
*/
class IDNotEngine : public SelectEngine
{
public:
    IDNotEngine( const SelectEnginePtr &part );
    ~IDNotEngine();
    
protected:
    SelectResult select(const SelectResult &mols, const PropertyMap &map) const;
    
private:
    SelectEnginePtr part;
};

/** Internal class providing the SelectEngine for objects 
    in a "join" expression
    
    @author Christopher Woods
*/
class IDJoinEngine : public SelectEngine
{
public:
    IDJoinEngine( const SelectEnginePtr &part );
    ~IDJoinEngine();
    
protected:
    SelectResult select(const SelectResult &mols, const PropertyMap &map) const;
    
private:
    SelectEnginePtr part;
};

/** Internal class providing the SelectEngine for objects 
    in a "subscript" expression
    
    @author Christopher Woods
*/
class IDSubScriptEngine : public SelectEngine
{
public:
    IDSubScriptEngine( const SelectEnginePtr &part, const RangeValue &val );
    ~IDSubScriptEngine();
    
protected:
    SelectResult select(const SelectResult &mols, const PropertyMap &map) const;
    
private:
    SelectEnginePtr part;
    RangeValue val;
};

/** Internal class providing the SelectEngine for objects 
    in a "with" expression
    
    @author Christopher Woods
*/
class IDWithEngine : public SelectEngine
{
public:
    IDWithEngine( IDObject obj, IDToken token, const SelectEnginePtr &part );
    ~IDWithEngine();
    
protected:
    SelectResult select(const SelectResult &mols, const PropertyMap &map) const;
    
private:
    IDObject obj;
    IDToken token;
    SelectEnginePtr part;
};

}

SIRE_END_HEADER

#endif
