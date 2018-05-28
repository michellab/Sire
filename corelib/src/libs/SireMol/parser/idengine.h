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

#include <QString>
#include <QRegularExpression>

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
    static SelectEnginePtr construct( IDObject obj, NameValues vals );
    ~IDNameEngine();
    
    ObjType objectType() const;
    
protected:
    IDNameEngine();
    SelectResult select(const SelectResult &mols, const PropertyMap &map) const;
    
private:
    SelectResult selectAtoms(const SelectResult &mols, bool use_parallel) const;
    SelectResult selectCutGroups(const SelectResult &mols, bool use_parallel) const;
    SelectResult selectResidues(const SelectResult &mols, bool use_parallel) const;
    SelectResult selectChains(const SelectResult &mols, bool use_parallel) const;
    SelectResult selectSegments(const SelectResult &mols, bool use_parallel) const;
    SelectResult selectMolecules(const SelectResult &mols, bool use_parallel) const;

    bool match(const QString &name) const;

    IDObject obj;
    QStringList names;
    QList<QRegularExpression> regexps;
};

/** Internal class providing the SelectEngine for objects 
    based on their numbers
    
    @author Christopher Woods
*/
class IDNumberEngine : public SelectEngine
{
public:
    static SelectEnginePtr construct( IDObject obj, RangeValues vals );
    ~IDNumberEngine();

    ObjType objectType() const;
    
protected:
    IDNumberEngine();
    SelectResult select(const SelectResult &mols, const PropertyMap &map) const;
    
private:
    SelectResult selectAtoms(const SelectResult &mols, bool use_parallel) const;
    SelectResult selectResidues(const SelectResult &mols, bool use_parallel) const;
    SelectResult selectMolecules(const SelectResult &mols, bool use_parallel) const;

    bool match(int val) const;

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
    static SelectEnginePtr construct( IDObject obj, RangeValues vals );
    ~IDIndexEngine();
    
    ObjType objectType() const;

protected:
    IDIndexEngine();
    SelectResult select(const SelectResult &mols, const PropertyMap &map) const;
    
private:
    SelectResult selectAtoms(const SelectResult &mols, bool use_parallel) const;
    SelectResult selectCutGroups(const SelectResult &mols, bool use_parallel) const;
    SelectResult selectResidues(const SelectResult &mols, bool use_parallel) const;
    SelectResult selectChains(const SelectResult &mols, bool use_parallel) const;
    SelectResult selectSegments(const SelectResult &mols, bool use_parallel) const;
    SelectResult selectMolecules(const SelectResult &mols, bool use_parallel) const;

    bool match(int val, int count) const;

    IDObject obj;
    RangeValues vals;
};

/** Internal class providing the SelectEngine for selecting by chemical element

    @author Christopher Woods
*/
class IDElementEngine : public SelectEngine
{
public:
    static SelectEnginePtr construct(const std::vector<SireMol::Element> &values);
    ~IDElementEngine();
    
    ObjType objectType() const;

protected:
    IDElementEngine();
    SelectResult select(const SelectResult &mols, const PropertyMap &map) const;
    
private:
    QSet<SireMol::Element> elements;
};

/** Internal class providing the SelectEngine for objects
    in an "and" expression
    
    @author Christopher Woods
*/
class IDAndEngine : public SelectEngine
{
public:
    static SelectEnginePtr construct(SelectEnginePtr part0, SelectEnginePtr part1);
    ~IDAndEngine();

    ObjType objectType() const;
    
    SelectEnginePtr simplify();
    
protected:
    IDAndEngine();
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
    static SelectEnginePtr construct(SelectEnginePtr part0, SelectEnginePtr part1);
    static SelectEnginePtr construct(QList<SelectEnginePtr> parts);
    ~IDOrEngine();

    ObjType objectType() const;

    SelectEnginePtr simplify();

protected:
    IDOrEngine();
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
    static SelectEnginePtr construct(SelectEnginePtr part);
    ~IDNotEngine();

    ObjType objectType() const;

    SelectEnginePtr simplify();
    
protected:
    IDNotEngine();
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
    static SelectEnginePtr construct(SelectEnginePtr part);
    ~IDJoinEngine();

    ObjType objectType() const;

    SelectEnginePtr simplify();
    
protected:
    IDJoinEngine();
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
    static SelectEnginePtr construct(SelectEnginePtr part, const RangeValue &val);
    ~IDSubScriptEngine();

    ObjType objectType() const;

    SelectEnginePtr simplify();
    
protected:
    IDSubScriptEngine();
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
    static SelectEnginePtr construct( IDObject obj, IDToken token, SelectEnginePtr part);
    ~IDWithEngine();

    ObjType objectType() const;

    SelectEnginePtr simplify();

protected:
    IDWithEngine();
    SelectResult select(const SelectResult &mols, const PropertyMap &map) const;
    
private:
    IDObject obj;
    IDToken token;
    SelectEnginePtr part;
};

/** Internal class used to select objects that are within a certain
    distance of other objects
    
    @author Christopher Woods
*/
class IDDistanceEngine : public SelectEngine
{
public:
    static SelectEnginePtr construct( IDObject obj, SireUnits::Dimension::Length distance,
                                      SelectEnginePtr part );
    
    static SelectEnginePtr construct( IDObject obj, IDCoordType typ,
                                      SireUnits::Dimension::Length distance,
                                      SelectEnginePtr part );
    
    ~IDDistanceEngine();
    
    ObjType objectType() const;
    
    SelectEnginePtr simplify();
    
protected:
    IDDistanceEngine();
    SelectResult select(const SelectResult &mols, const PropertyMap &map) const;
    
private:
    IDObject obj;
    IDCoordType typ;
    SelectEnginePtr part;
    double distance;
};

}

SIRE_END_HEADER

#endif
