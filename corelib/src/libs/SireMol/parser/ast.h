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

#ifndef SIREMOL_PARSER_AST_H
#define SIREMOL_PARSER_AST_H

#include "SireUnits/dimensions.h"
#include "SireUnits/units.h"
#include "SireMol/select.h"
#include "SireMol/element.h"

#include "SireError/errors.h"

#include "parser.h"

SIRE_BEGIN_HEADER

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
    enum IDObject { ID_UNKNOWN = 0, ATOM = 1, CUTGROUP = 2,
                    RESIDUE = 3, CHAIN = 4, SEGMENT = 5, MOLECULE = 6 };
    
    QString idobject_to_string(IDObject obj);

    /** The different types of number that can be identified */
    enum IDNumType { ID_TYP_UNKNOWN = 0, ID_NUMBER = 1, ID_INDEX = 2 };
    
    QString idnumtype_to_string(IDNumType typ);

    /** The different types of logical operation */
    enum IDOperation { ID_OP_UNKNOWN = 0, ID_AND = 1, ID_OR = 2 };
    
    QString idoperation_to_string(IDOperation op);

    /** The different types of value comparison */
    enum IDComparison { ID_CMP_UNKNOWN = 0, ID_CMP_LT = 1, ID_CMP_LE = 2,
                        ID_CMP_EQ = 3, ID_CMP_NE = 4, ID_CMP_GT = 5, ID_CMP_GE = 6 };
    
    QString idcomparison_to_string(IDComparison cmp);

    /** The different miscellaneous tokens */
    enum IDToken { ID_TOKEN_UNKNOWN = 0, ID_WHERE = 1, ID_WITH = 2, ID_IN = 3 };
    
    QString idtoken_to_string(IDToken token);
    
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

    struct IDName;          
    struct IDNumber;        
    struct IDBinary;        
    struct IDWith;          
    struct IDWhere;         
    struct IDNot;           
    struct IDSubscript;     
    struct IDWithin;        
    struct IDUser;
    struct IDJoin;

    struct IDWhereCompare;
    struct IDWhereWithin;

    struct IDElement;

    struct Expression;
    struct ExpressionPart;
    struct Node;
    
    /** Base holder for all of the different ID expressions */
    using ExpressionVariant = boost::variant<boost::recursive_wrapper<IDName>,
                                             boost::recursive_wrapper<IDNumber>,
                                             boost::recursive_wrapper<IDElement>,
                                             boost::recursive_wrapper<IDBinary>,
                                             boost::recursive_wrapper<IDWith>,
                                             boost::recursive_wrapper<IDWhere>,
                                             boost::recursive_wrapper<IDNot>,
                                             boost::recursive_wrapper<IDSubscript>,
                                             boost::recursive_wrapper<IDWithin>,
                                             boost::recursive_wrapper<IDUser>,
                                             boost::recursive_wrapper<IDJoin>,
                                             boost::recursive_wrapper<ExpressionPart> >;
    
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
        RangeValue() : start(0), end(0), step(1), _c(0)
        {}
    
        int start;
        int end;
        int step;
        
        RangeValue& operator+=(const int val)
        {
            if (_c == 0)
            {
                start = val;
                end = val;
                step = 1;
                _c += 1;
            }
            else if (_c == 1)
            {
                end = val;
                step = 1;
                _c += 1;
            }
            else if (_c == 2)
            {
                step = val;
                _c += 1;
            }
            else
                throw SireError::program_bug( QObject::tr(
                    "Should not add more than three values to a range..."), CODELOC );
        
            return *this;
        }
        
        QString toString() const;

    private:
        int _c;
    };

    /** Struct that holds a general selection expression */
    struct Expression
    {
        ExpressionVariant value;
        
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
        std::vector<SireMol::Element> values;
        
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
    
    /** Struct that holds a "with" expression, e.g. molecules with resname ala */
    struct IDWith
    {
        IDObject name;
        IDToken token;
        Expression value;
        
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
}

BOOST_FUSION_ADAPT_STRUCT( AST::IDWithin,
                           (AST::IDObject,name)
                           (AST::LengthValue,distance)
                           (AST::Expression,value)
                         )

BOOST_FUSION_ADAPT_STRUCT( AST::LengthValue,
                           (double,value)
                           (SireUnits::Dimension::Length,unit)
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

BOOST_FUSION_ADAPT_STRUCT( AST::IDName,
                           (AST::IDObject, name),
                           (AST::NameValues, values)
                         )

BOOST_FUSION_ADAPT_STRUCT( AST::IDElement,
                           (std::vector<SireMol::Element>, values)
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
                           (AST::IDObject,name)
                           (AST::IDToken,token)
                           (AST::Expression,value)
                         )

BOOST_FUSION_ADAPT_STRUCT( AST::Expression,
                           (AST::ExpressionVariant, value)
                         )

BOOST_FUSION_ADAPT_STRUCT( AST::ExpressionPart,
                           (AST::ExpressionVariant, value)
                         )

SIRE_END_HEADER

#endif
