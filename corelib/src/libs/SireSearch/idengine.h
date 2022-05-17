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

#ifndef SIRESEARCH_IDENGINE_H
#define SIRESEARCH_IDENGINE_H

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

    QString toString() const;

protected:
    IDNameEngine();
    SelectResult select(const SelectResult &mols, const PropertyMap &map) const;

private:
    template<class T>
    SelectResult searchName(const SelectResult &mols, bool use_parallel) const;

    SelectResult searchMolName(const SelectResult &mols, bool use_parallel) const;

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
    template<class T>
    SelectResult searchNum(const SelectResult &mols, bool use_parallel) const;

    SelectResult searchMolNum(const SelectResult &mols, bool use_parallel) const;

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
    template<class T>
    SelectResult searchIdx(const SelectResult &mols, bool use_parallel) const;

    SelectResult searchMolIdx(const SelectResult &mols, bool use_parallel) const;

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
 *  in a "bonds" expression
*/
class IDBondEngine : public SelectEngine
{
public:
    static SelectEnginePtr construct( IDBondToken from_token,
                                      SelectEnginePtr from_value,
                                      IDBondToken to_token,
                                      SelectEnginePtr to_value);
    ~IDBondEngine();

    ObjType objectType() const;

    SelectEnginePtr simplify();

protected:
    IDBondEngine();
    SelectResult select(const SelectResult &mols, const PropertyMap &map) const;

private:
    IDBondToken from_token;
    SelectEnginePtr from_value;
    IDBondToken to_token;
    SelectEnginePtr to_value;
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

/** Internal class used to find atoms by their mass */
class IDMassEngine : public SelectEngine
{
public:
    static SelectEnginePtr construct(IDObject obj, IDComparison compare, SireUnits::Dimension::MolarMass mass);
    ~IDMassEngine();

    ObjType objectType() const;

protected:
    IDMassEngine();
    SelectResult select(const SelectResult &mols, const PropertyMap &map) const;

    template<class T>
    SelectResult select_t(const SelectResult &mols, const PropertyMap &map) const;

    SelectResult select_mols(const SelectResult &mols, const PropertyMap &map) const;

    SelectResult select_bonds(const SelectResult &mols, const PropertyMap &map) const;

private:
    IDObject obj;
    IDComparison compare;
    SireUnits::Dimension::MolarMass value;
};

/** Internal class used to find atoms by their charge */
class IDChargeEngine : public SelectEngine
{
public:
    static SelectEnginePtr construct(IDObject obj, IDComparison compare, SireUnits::Dimension::Charge mass);
    ~IDChargeEngine();

    ObjType objectType() const;

protected:
    IDChargeEngine();
    SelectResult select(const SelectResult &mols, const PropertyMap &map) const;

    template<class T>
    SelectResult select_t(const SelectResult &mols, const PropertyMap &map) const;

    SelectResult select_mols(const SelectResult &mols, const PropertyMap &map) const;

    SelectResult select_bonds(const SelectResult &mols, const PropertyMap &map) const;

private:
    IDObject obj;
    IDComparison compare;
    SireUnits::Dimension::Charge value;
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

/** Internal class used to select objects that are within a certain
    distance of a 3D position vector.

    @author Lester Hedges
*/
class IDDistanceVectorEngine : public SelectEngine
{
public:
    static SelectEnginePtr construct( IDObject obj, SireUnits::Dimension::Length distance,
                                      VectorValue position );

    static SelectEnginePtr construct( IDObject obj, IDCoordType typ,
                                      SireUnits::Dimension::Length distance,
                                      VectorValue position );

    ~IDDistanceVectorEngine();

    ObjType objectType() const;

    SelectEnginePtr simplify();

protected:
    IDDistanceVectorEngine();
    SelectResult select(const SelectResult &mols, const PropertyMap &map) const;

private:
    IDObject obj;
    IDCoordType typ;
    VectorValue position;
    double distance;
};

/** Internal class used to select all objects

    @author Christopher Woods
*/
class IDAllEngine : public SelectEngine
{
public:
    static SelectEnginePtr construct();

    ~IDAllEngine();

    ObjType objectType() const;

protected:
    IDAllEngine();
    SelectResult select(const SelectResult &mols, const PropertyMap &map) const;
};

/** Internal class used to select all water molecules.

    @author Lester Hedges
*/
class IDWaterEngine : public SelectEngine
{
public:
    static SelectEnginePtr construct();

    ~IDWaterEngine();

    ObjType objectType() const;

protected:
    IDWaterEngine();
    SelectResult select(const SelectResult &mols, const PropertyMap &map) const;
};

/** Internal class used to select all protein molecules */
class IDProteinEngine : public SelectEngine
{
public:
    static SelectEnginePtr construct();

    ~IDProteinEngine();

    ObjType objectType() const;

protected:
    IDProteinEngine();
    SelectResult select(const SelectResult &mols, const PropertyMap &map) const;
};

/** Internal class used to select all perturbable molecules.

    @author Lester Hedges
*/
class IDPerturbableEngine : public SelectEngine
{
public:
    static SelectEnginePtr construct();

    ~IDPerturbableEngine();

    ObjType objectType() const;

protected:
    IDPerturbableEngine();
    SelectResult select(const SelectResult &mols, const PropertyMap &map) const;
};

}

SIRE_END_HEADER

#endif
