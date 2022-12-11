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

#ifndef SIRESEARCH_PARSER_AST_H
#define SIRESEARCH_PARSER_AST_H

#include "SireUnits/dimensions.h"
#include "SireUnits/units.h"
#include "SireMol/select.h"
#include "SireMol/element.h"

#include "SireError/errors.h"

#include "parser.h"

namespace SireBase
{
class Slice;
}

SIRE_BEGIN_HEADER

// need to use this to increase the number of variants that can be
// held in ExpressionVariant from the default of 20. This is
// not a long term solution!
#define BOOST_MPL_CFG_NO_PREPROCESSED_HEADERS
#define BOOST_MPL_LIMIT_LIST_SIZE 30

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/support_line_pos_iterator.hpp>

#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/spirit/include/phoenix_object.hpp>

#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/io.hpp>

#include <boost/shared_ptr.hpp>

// A lot of the below code is heavily inspired by
// https://medium.com/@alinakipoglu/parsing-with-spirit-qi-fcaeaf4357b3

/** Namespace holding the objects used in the abstract syntax tree */
namespace AST
{
    using SireMol::parser::SelectEnginePtr;

    /** The different objects that can be identified */
    enum IDObject { ID_UNKNOWN = 0,
                    ATOM = 1,
                    BOND = 2,
                    ANGLE = 3,
                    DIHEDRAL = 4,
                    IMPROPER = 5,
                    RESIDUE = 6,
                    CHAIN = 7,
                    SEGMENT = 8,
                    CUTGROUP = 9,
                    MOLECULE = 10,
                    VIEW = 99 };

    QString idobject_to_string(IDObject obj);

    /** The different types of number that can be identified */
    enum IDNumType { ID_TYP_UNKNOWN = 0, ID_NUMBER = 1, ID_INDEX = 2 };

    QString idnumtype_to_string(IDNumType typ);

    /** The different types of logical operation */
    enum IDOperation { ID_OP_UNKNOWN = 0, ID_AND = 1, ID_OR = 2 };

    QString idoperation_to_string(IDOperation op);

    /** The different types of value comparison */
    enum IDComparison { ID_CMP_UNKNOWN = 0, ID_CMP_LT = 1, ID_CMP_LE = 2,
                        ID_CMP_EQ = 3, ID_CMP_NE = 4, ID_CMP_GT = 5, ID_CMP_GE = 6,
                        ID_CMP_AE = 7 };

    QString idcomparison_to_string(IDComparison cmp);

    /** The different miscellaneous tokens */
    enum IDToken { ID_TOKEN_UNKNOWN = 0, ID_WHERE = 1, ID_WITH = 2, ID_IN = 3 };

    QString idtoken_to_string(IDToken token);

    /** The different bond tokens */
    enum IDBondToken { ID_BOND_UNKNOWN = 0, ID_BOND_FROM = 1, ID_BOND_TO = 2 };

    QString idbondtoken_to_string(IDBondToken token);

    /** The different types of coordinate */
    enum IDCoordType { ID_COORD_UNKNOWN = 0, ID_COORD_CENTER = 1, ID_COORD_CENTER_X = 2,
                       ID_COORD_CENTER_Y = 3, ID_COORD_CENTER_Z = 4, ID_COORD_MAX = 5,
                       ID_COORD_MAX_X = 6, ID_COORD_MAX_Y = 7, ID_COORD_MAX_Z = 8,
                       ID_COORD_MIN = 9, ID_COORD_MIN_X = 10, ID_COORD_MIN_Y = 11,
                       ID_COORD_MIN_Z = 12, ID_COORD_X = 13, ID_COORD_Y = 14, ID_COORD_Z = 15,
                       ID_COORD_CLOSEST = 16 };

    QString idcoordtype_to_string(IDCoordType typ);

    struct NameValue;
    struct RangeValue;
    struct CompareValue;
    struct RegExpValue;
    struct LengthValue;
    struct VectorValue;

    struct IDNull;
    struct IDName;
    struct IDNumber;
    struct IDBinary;
    struct IDWith;
    struct IDWhere;
    struct IDCount;
    struct IDNot;
    struct IDSubscript;
    struct IDWithin;
    struct IDWithinVector;
    struct IDUser;
    struct IDJoin;
    struct IDAll;
    struct IDWater;
    struct IDProtein;
    struct IDPerturbable;
    struct IDBond;
    struct IDProperty;
    struct IDMass;
    struct IDCharge;
    struct IDCmpMass;
    struct IDCmpCharge;
    struct IDObjMass;
    struct IDObjCmpMass;
    struct IDObjCharge;
    struct IDObjCmpCharge;

    struct IDWhereCompare;
    struct IDWhereWithin;

    struct IDElement;

    struct Expression;
    struct ExpressionPart;
    struct Node;

    /** Base holder for all of the different ID expressions */
    using ExpressionVariant = boost::variant<boost::recursive_wrapper<IDNull>,
                                             boost::recursive_wrapper<IDName>,
                                             boost::recursive_wrapper<IDNumber>,
                                             boost::recursive_wrapper<IDElement>,
                                             boost::recursive_wrapper<IDBinary>,
                                             boost::recursive_wrapper<IDWith>,
                                             boost::recursive_wrapper<IDWhere>,
                                             boost::recursive_wrapper<IDCount>,
                                             boost::recursive_wrapper<IDNot>,
                                             boost::recursive_wrapper<IDSubscript>,
                                             boost::recursive_wrapper<IDWithin>,
                                             boost::recursive_wrapper<IDWithinVector>,
                                             boost::recursive_wrapper<IDUser>,
                                             boost::recursive_wrapper<IDJoin>,
                                             boost::recursive_wrapper<IDAll>,
                                             boost::recursive_wrapper<IDBond>,
                                             boost::recursive_wrapper<IDProperty>,
                                             boost::recursive_wrapper<IDWater>,
                                             boost::recursive_wrapper<IDProtein>,
                                             boost::recursive_wrapper<IDPerturbable>,
                                             boost::recursive_wrapper<IDMass>,
                                             boost::recursive_wrapper<IDCharge>,
                                             boost::recursive_wrapper<IDCmpMass>,
                                             boost::recursive_wrapper<IDCmpCharge>,
                                             boost::recursive_wrapper<IDObjMass>,
                                             boost::recursive_wrapper<IDObjCharge>,
                                             boost::recursive_wrapper<IDObjCmpMass>,
                                             boost::recursive_wrapper<IDObjCmpCharge>,
                                             boost::recursive_wrapper<ExpressionPart> >;

    QString expression_to_string(const ExpressionVariant &expression);

    /** Base holder for strings or regular expressions */
    using NameVariant = boost::variant<boost::recursive_wrapper<RegExpValue>,
                                       std::string>;

    /** Base holder for ranges of value comparisons */
    using RangeVariant = boost::variant<boost::recursive_wrapper<RangeValue>,
                                        boost::recursive_wrapper<CompareValue> >;

    /** Base holder for different methods of holding a "where" expression */
    using IDWhereVariant = boost::variant<boost::recursive_wrapper<IDWhereCompare>,
                                          boost::recursive_wrapper<IDWhereWithin> >;

    using IDNames = std::vector<IDName>;
    using Expressions = std::vector<Expression>;
    using NameValues = std::vector<NameValue>;
    using RangeValues = std::vector<RangeVariant>;

    /** Visitor used to get string representations from a boost::variant */
    class qstring_visitor : public boost::static_visitor<QString>
    {
    public:
        /** In general, use the .toString() function from a class */
        template<class T>
        QString operator()(const T &value) const
        {
            return value.toString();
        }

        /** Convert a std::string into a QString */
        QString operator()(const std::string &string) const
        {
            return QString("'%1'").arg(QString::fromStdString(string));
        }
    };

    /** Visitorused to get SelectEnginePtr values from a boost::variant */
    class engine_visitor : public boost::static_visitor<SelectEnginePtr>
    {
    public:
        /** In general, use the .toEngine() function from a class */
        template<class T>
        SelectEnginePtr operator()(const T &value) const
        {
            return value.toEngine();
        }
    };

    /** Struct that holds a length or distance, with units */
    struct LengthValue
    {
        LengthValue() : value(0), unit(1.0)
        {}

        LengthValue& operator+=(double v)
        {
            value = v;
            return *this;
        }

        LengthValue& operator+=(SireUnits::Dimension::Length v)
        {
            unit = v;
            return *this;
        }

        double value;
        SireUnits::Dimension::Length unit;

        QString toString() const;
    };

    /** Struct that holds a vector or point in space (within units) */
    struct VectorValue
    {
        LengthValue x;
        LengthValue y;
        LengthValue z;

        VectorValue() : _c(0)
        {}

        VectorValue& operator+=(const LengthValue &val)
        {
            if (_c == 0)
            {
                x = val;
                y = val;
                z = val;
                _c += 1;
            }
            else if (_c == 1)
            {
                y = val;
                z = val;
                _c += 1;
            }
            else if (_c == 2)
            {
                z = val;
                _c += 1;
            }
            else
                throw SireMol::parse_error( QObject::tr(
                    "Cannot add more than there points to a vector"), CODELOC );

            return *this;
        }

        QString toString() const;

    private:
        int _c;

    };

    /** Struct that holds a regular expression */
    struct RegExpValue
    {
        RegExpValue() : is_case_sensitive(true)
        {}

        template<class T>
        RegExpValue& operator+=(const T &val)
        {
            value = val;
            return *this;
        }

        template<class T>
        RegExpValue& operator*=(const T &val)
        {
            is_case_sensitive = false;
            return *this;
        }

        QString toString() const;

        std::string value;
        bool is_case_sensitive;
    };

    /** Struct that holds a name (either a string or regular expression) */
    struct NameValue
    {
        NameVariant value;

        QString toString() const
        {
            return boost::apply_visitor( qstring_visitor(), value );
        }
    };

    /** Struct that holds a comparison against an integer */
    struct CompareValue
    {
        IDComparison compare;
        int value;

        QString toString() const;
    };

    /** Struct that holds a range of integers */
    struct RangeValue
    {
        RangeValue() : _c(0)
        {}

        std::shared_ptr<int> start;
        std::shared_ptr<int> stop;
        std::shared_ptr<int> step;

        RangeValue& operator*=(const int val)
        {
            _c += 1;
            return *this;
        }

        RangeValue& operator+=(const int val)
        {
            if (_c == 0)
            {
                start.reset(new int(val));
            }
            else if (_c == 1)
            {
                stop.reset(new int(val));
            }
            else if (_c == 2)
            {
                step.reset(new int(val));
            }
            else
            {
                qDebug() << "extra +=" << val;
            }

            return *this;
        }

        QString toString() const;

        SireBase::Slice toSlice() const;

    private:
        int _c;
    };

    /** Null struct for empty values */
    struct IDNull
    {
        IDObject name;
        NameValues values;

        QString toString() const;
        SelectEnginePtr toEngine() const;
    };

    /** Struct that holds a general selection expression */
    struct Expression
    {
        ExpressionVariant value = IDNull();

        QString toString() const;

        SelectEnginePtr toEngine() const;
    };

    /** Struct that holds an atomic part of a selection expression */
    struct ExpressionPart
    {
        ExpressionVariant value;

        QString toString() const;

        SelectEnginePtr toEngine() const;
    };

    /** Struct that holds an ID token that matches everything */
    struct IDAll
    {
        IDAll(IDObject object = MOLECULE) : name(object)
        {}

        IDObject name;

        QString toString() const;
        SelectEnginePtr toEngine() const;
    };

    /** Struct to hold an ID token that matches water molecules */
    struct IDWater
    {
        QString toString() const;

        SelectEnginePtr toEngine() const;
    };

    /** Struct to hold an ID token that matches protein molecules */
    struct IDProtein
    {
        QString toString() const;

        SelectEnginePtr toEngine() const;
    };

    /** Struct to hold an ID token that matches perturbable molecules */
    struct IDPerturbable
    {
        QString toString() const;

        SelectEnginePtr toEngine() const;
    };

    /** Struct that holds an ID token that represents a user-supplied selection */
    struct IDUser
    {
        std::string token;
        ExpressionVariant value;

        IDUser()
        {}

        IDUser(const std::string t, const ExpressionVariant &s)
            : token(t), value(s)
        {}

        QString toString() const;

        SelectEnginePtr toEngine() const;
    };

    /** The root node of the AST - this holds a set of Expressions */
    struct Node
    {
        Expressions values;

        QString toString() const;

        SelectEnginePtr toEngine() const;
    };
    /** Struct that holds a name and associated values */
    struct IDName
    {
        IDObject name;
        NameValues values;

        QString toString() const;

        SelectEnginePtr toEngine() const;
    };

    /** Struct that holds a number and associated values */
    struct IDNumber
    {
        IDObject name;
        IDNumType numtype;
        RangeValues values;

        IDNumber() : name(ID_UNKNOWN), numtype(ID_TYP_UNKNOWN)
        {}

        IDNumber& operator+=(const QPair<IDObject,IDNumType> &v)
        {
            name = v.first;
            numtype = v.second;
            return *this;
        }

        IDNumber& operator+=(const RangeValues &vals)
        {
            values = vals;
            return *this;
        }

        QString toString() const;

        SelectEnginePtr toEngine() const;
    };

    /** Structs that holds a list of elements */
    struct IDElement
    {
        std::vector<QString> values;

        QString toString() const;

        SelectEnginePtr toEngine() const;
    };

    /** Struct that holds a binary expression, e.g. something and other */
    struct IDBinary
    {
        Expression part0;
        IDOperation operation;
        Expression part1;

        QString toString() const;

        SelectEnginePtr toEngine() const;
    };

    /** Struct that holds an ID token that represents a Property expression,
     * e.g. property charge, property charge > 0, property perturbable == True
     */
    struct IDProperty
    {
        IDObject name;
        std::string property;
        IDComparison compare;
        std::string value;

        IDProperty& operator/=(int)
        {
            //reset the search
            property.clear();
            value.clear();
            name = AST::VIEW;
            compare = AST::ID_CMP_EQ;
            return *this;
        }

        IDProperty& operator+=(const std::string &p)
        {
            property = p;
            return *this;
        }

        IDProperty& operator+=(const IDObject &n)
        {
            this->operator/=(1);
            name = n;
            return *this;
        }

        IDProperty& operator+=(const IDComparison &c)
        {
            compare = c;
            return *this;
        }

        IDProperty& operator*=(const std::string &v)
        {
            value = v;
            return *this;
        }

        QString toString() const;

        SelectEnginePtr toEngine() const;
    };

    /** Struct that holds an ID token that represents a Bond expression,
     *  e.g. bonds within resnum 1, bonds in resnum 1,
     *  bonds from atomnum 1 to atomnum 2
    */
    struct IDBond
    {
        IDBondToken from_token;
        Expression from_value;
        IDBondToken to_token;
        Expression to_value;

        QString toString() const;

        SelectEnginePtr toEngine() const;
    };

    /** Struct that holds a "with" expression, e.g. molecules with resname ala */
    struct IDWith
    {
        IDWith()
        {}

        Expression value0;
        IDToken token;
        Expression value1;

        IDWith& operator+=(const Expression &v)
        {
            value0 = v;
            return *this;
        }

        IDWith& operator+=(const IDToken &t)
        {
            token = t;
            return *this;
        }

        IDWith& operator*=(const Expression &v)
        {
            value1 = v;
            return *this;
        }

        QString toString() const;

        SelectEnginePtr toEngine() const;
    };

    /** Struct that holds a "mass" expression, e.g. 1 g_per_mol */
    struct IDMass
    {
        double value;
        SireUnits::Dimension::MolarMass units;

        IDMass() : value(0.0), units(1.0)
        {}

        IDMass(const IDMass &other) : value(other.value), units(other.units)
        {}

        IDMass& operator+=(double v)
        {
            value = v;
            return *this;
        }

        IDMass& operator+=(SireUnits::Dimension::MolarMass u)
        {
            units = u;
            return *this;
        }

        QString toString() const;

        SelectEnginePtr toEngine() const;
    };

    /** Struct that holds a mass comparison, e.g. mass > 1 g_per_mol */
    struct IDCmpMass
    {
        IDComparison compare;
        IDMass value;

        QString toString() const;

        SelectEnginePtr toEngine() const;
    };

    /** Struct that holds an object mass, e.g. atom mass 1 g_per_mol */
    struct IDObjMass
    {
        IDObject name;
        IDMass value;

        QString toString() const;

        SelectEnginePtr toEngine() const;
    };

    /** Struct that holds an object mass comparison, e.g. atom mass > 1 g_per_mol */
    struct IDObjCmpMass
    {
        IDObject name;
        IDComparison compare;
        IDMass value;

        QString toString() const;

        SelectEnginePtr toEngine() const;
    };

    /** Struct that holds a "charge" expression, e.g. 1 e */
    struct IDCharge
    {
        double value;
        SireUnits::Dimension::Charge units;

        IDCharge() : value(0.0), units(1.0)
        {}

        IDCharge(const IDCharge &other) : value(other.value), units(other.units)
        {}

        IDCharge& operator+=(double v)
        {
            value = v;
            return *this;
        }

        IDCharge& operator+=(SireUnits::Dimension::Charge u)
        {
            units = u;
            return *this;
        }

        QString toString() const;

        SelectEnginePtr toEngine() const;
    };

    /** Struct that holds a charge comparison, e.g. charge > 1 e */
    struct IDCmpCharge
    {
        IDComparison compare;
        IDCharge value;

        QString toString() const;

        SelectEnginePtr toEngine() const;
    };

    /** Struct that holds an object charge, e.g. atom char 1 e */
    struct IDObjCharge
    {
        IDObject name;
        IDCharge value;

        QString toString() const;

        SelectEnginePtr toEngine() const;
    };

    /** Struct that holds an object charge comparison, e.g. atom charge > 1 e */
    struct IDObjCmpCharge
    {
        IDObject name;
        IDComparison compare;
        IDCharge value;

        QString toString() const;

        SelectEnginePtr toEngine() const;
    };

    /** Struct that holds "count(*) < X" expressions */
    struct IDCount
    {
        Expression object;
        IDComparison compare;
        int value;

        IDCount& operator+=(const Expression &o)
        {
            object = o;
            return *this;
        }

        IDCount& operator+=(const IDComparison &c)
        {
            compare = c;
            return *this;
        }

        IDCount& operator+=(int v)
        {
            value = v;
            return *this;
        }

        QString toString() const;

        SelectEnginePtr toEngine() const;
    };

    /** Struct that holds a "where within"
        expression, e.g. residues where center is within 5 A of resname /lig/i */
    struct IDWhereWithin
    {
        LengthValue distance;
        Expression value;

        QString toString() const;

        SelectEnginePtr toEngine(IDObject name, IDCoordType typ) const;
    };

    /** Struct that holds a "where" comparison expression, e.g.
        residues where center.x >= 5 */
    struct IDWhereCompare
    {
        IDComparison compare;
        VectorValue value;

        QString toString() const;

        SelectEnginePtr toEngine(IDObject name, IDCoordType typ) const;
    };

    /** Struct that holds a general "where" expression */
    struct IDWhere
    {
        IDObject name;
        IDCoordType typ;
        IDWhereVariant value;

        QString toString() const;

        SelectEnginePtr toEngine() const;
    };

    /** Struct to hold a negated (not) expression */
    struct IDNot
    {
        Expression value;

        QString toString() const;

        SelectEnginePtr toEngine() const;
    };

    /** Struct to hold a join expression */
    struct IDJoin
    {
        Expression value;

        QString toString() const;

        SelectEnginePtr toEngine() const;
    };

    /** Struct to hold a subscripted expression, e.g. {something}[0:10:2] */
    struct IDSubscript
    {
        Expression value;
        RangeValue range;

        QString toString() const;

        SelectEnginePtr toEngine() const;
    };

    /** Struct to hold expressions that select based on being within a distance */
    struct IDWithin
    {
        IDObject name;
        LengthValue distance;
        Expression value;

        QString toString() const;

        SelectEnginePtr toEngine() const;
    };

    /** Struct to hold expressions that select based on being within a distance to a point*/
    struct IDWithinVector
    {
        IDObject name;
        LengthValue distance;
        VectorValue value;

        QString toString() const;

        SelectEnginePtr toEngine() const;
    };
}

BOOST_FUSION_ADAPT_STRUCT( AST::IDWithin,
                           (AST::IDObject,name)
                           (AST::LengthValue,distance)
                           (AST::Expression,value)
                         )

BOOST_FUSION_ADAPT_STRUCT( AST::IDCharge,
                           (double,value)
                           (SireUnits::Dimension::Charge,units)
                         )

BOOST_FUSION_ADAPT_STRUCT( AST::IDMass,
                           (double,value)
                           (SireUnits::Dimension::MolarMass,units)
                         )

BOOST_FUSION_ADAPT_STRUCT( AST::IDCmpCharge,
                           (AST::IDComparison,compare)
                           (AST::IDCharge,value)
                         )

BOOST_FUSION_ADAPT_STRUCT( AST::IDCmpMass,
                           (AST::IDComparison,compare)
                           (AST::IDMass,value)
                         )

BOOST_FUSION_ADAPT_STRUCT( AST::IDObjCharge,
                           (AST::IDObject, name)
                           (AST::IDCharge,value)
                         )

BOOST_FUSION_ADAPT_STRUCT( AST::IDObjMass,
                           (AST::IDObject, name)
                           (AST::IDMass,value)
                         )

BOOST_FUSION_ADAPT_STRUCT( AST::IDObjCmpCharge,
                           (AST::IDObject, name)
                           (AST::IDComparison,compare)
                           (AST::IDCharge,value)
                         )

BOOST_FUSION_ADAPT_STRUCT( AST::IDObjCmpMass,
                           (AST::IDObject, name)
                           (AST::IDComparison,compare)
                           (AST::IDMass,value)
                         )

BOOST_FUSION_ADAPT_STRUCT( AST::IDWithinVector,
                           (AST::IDObject,name)
                           (AST::LengthValue,distance)
                           (AST::VectorValue,value)
                         )

BOOST_FUSION_ADAPT_STRUCT( AST::LengthValue,
                           (double,value)
                           (SireUnits::Dimension::Length,unit)
                         )

BOOST_FUSION_ADAPT_STRUCT( AST::IDCount,
                           (AST::Expression,object)
                           (AST::IDComparison,compare)
                           (int,value)
                         )

BOOST_FUSION_ADAPT_STRUCT( AST::IDWhereWithin,
                           (AST::LengthValue,distance)
                           (AST::Expression,value)
                         )

BOOST_FUSION_ADAPT_STRUCT( AST::IDWhereCompare,
                           (AST::IDComparison,compare)
                           (AST::VectorValue,value)
                         )

BOOST_FUSION_ADAPT_STRUCT( AST::IDWhere,
                           (AST::IDObject,name)
                           (AST::IDCoordType,typ)
                           (AST::IDWhereVariant,value)
                         )

BOOST_FUSION_ADAPT_STRUCT( AST::NameValue,
                           (AST::NameVariant,value)
                         )

BOOST_FUSION_ADAPT_STRUCT( AST::CompareValue,
                           (AST::IDComparison,compare)
                           (int,value)
                         )

BOOST_FUSION_ADAPT_STRUCT( AST::IDNot,
                           (AST::Expression,value)
                         )

BOOST_FUSION_ADAPT_STRUCT( AST::IDJoin,
                           (AST::Expression,value)
                         )

BOOST_FUSION_ADAPT_STRUCT( AST::IDSubscript,
                           (AST::Expression,value)
                           (AST::RangeValue,range)
                         )

BOOST_FUSION_ADAPT_STRUCT( AST::RegExpValue,
                           (std::string,value)
                           (bool,is_case_sensitive)
                         )

BOOST_FUSION_ADAPT_STRUCT( AST::Node,
                           (AST::Expressions,values)
                         )

BOOST_FUSION_ADAPT_STRUCT( AST::IDNull,
                           (AST::IDObject, name),
                           (AST::NameValues, values)
                         )

BOOST_FUSION_ADAPT_STRUCT( AST::IDName,
                           (AST::IDObject, name),
                           (AST::NameValues, values)
                         )

BOOST_FUSION_ADAPT_STRUCT( AST::IDElement,
                           (std::vector<QString>, values)
                         )

BOOST_FUSION_ADAPT_STRUCT( AST::IDNumber,
                           (AST::IDObject, name),
                           (AST::IDNumType, numtype),
                           (AST::RangeValues, values)
                         )

BOOST_FUSION_ADAPT_STRUCT( AST::IDBinary,
                           (AST::Expression, part0),
                           (AST::IDOperation, operation),
                           (AST::Expression, part1)
                         )

BOOST_FUSION_ADAPT_STRUCT( AST::IDWith,
                           (AST::Expression,value0)
                           (AST::IDToken,token)
                           (AST::Expression,value1)
                         )

BOOST_FUSION_ADAPT_STRUCT( AST::IDBond,
                           (AST::IDBondToken,from_token)
                           (AST::Expression,from_value)
                           (AST::IDBondToken,to_token)
                           (AST::Expression,to_value)
                         )

BOOST_FUSION_ADAPT_STRUCT( AST::Expression,
                           (AST::ExpressionVariant, value)
                         )

BOOST_FUSION_ADAPT_STRUCT( AST::ExpressionPart,
                           (AST::ExpressionVariant, value)
                         )

SIRE_END_HEADER

#endif
