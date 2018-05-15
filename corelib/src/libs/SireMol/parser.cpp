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

#include "parser.h"

#include <QDebug>

using namespace SireMol;

const char* parse_error::typeName()
{
    return QMetaType::typeName( qMetaTypeId<parse_error>() );
}

static const RegisterMetaType<parse_error> r_parse;

// include boost::spirit::qi for parsing
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/support_line_pos_iterator.hpp>

#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/spirit/include/phoenix_object.hpp>

#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/io.hpp>

// A lot of the below code is heavily inspired by
// https://medium.com/@alinakipoglu/parsing-with-spirit-qi-fcaeaf4357b3

///////////
/////////// First define the objects used for the abstract syntax tree
///////////

namespace AST
{
    enum IDObject { ID_UNKNOWN = 0, ATOM = 1, CUTGROUP = 2,
                    RESIDUE = 3, CHAIN = 4, SEGMENT = 5, MOLECULE = 6 };
    
    QString idobject_to_string(IDObject obj)
    {
        switch(obj)
        {
        case ID_UNKNOWN:
            return "unknown";
        case ATOM:
            return "atom";
        case CUTGROUP:
            return "cutgroup";
        case RESIDUE:
            return "residue";
        case CHAIN:
            return "chain";
        case SEGMENT:
            return "segment";
        case MOLECULE:
            return "molecule";
        default:
            return "unknown";
        }
    }
    
    enum IDOperation { ID_OP_UNKNOWN = 0, ID_AND = 1, ID_OR = 2 };
    
    QString idoperation_to_string(IDOperation op)
    {
        switch(op)
        {
        case ID_AND:
            return "and";
        case ID_OR:
            return "or";
        default:
            return "unknown";
        }
    }
    
    struct NameValue;       // holder for a generic name value
    struct RegExpValue;     // holder for a regular expression value
    struct IDName;          // a named generic value
    struct IDBinary;        // a binary ID expression
    struct Expression;  // holder for a generic expression
    struct ExpressionPart;  //holder for a generic part of an expression
    struct Node;       // a node in the tree
    
    using ExpressionVariant = boost::variant<boost::recursive_wrapper<IDName>,
                                             boost::recursive_wrapper<IDBinary>,
                                             boost::recursive_wrapper<ExpressionPart> >;
    
    using NameVariant = boost::variant<boost::recursive_wrapper<RegExpValue>,
                                       std::string>;
    
    using IDNames = std::vector<IDName>;
    using Expressions = std::vector<Expression>;
    using NameValues = std::vector<NameValue>;

    // a holder for a regular expression
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
    
        QString toString() const
        {
            QString qstr = QString::fromStdString(value);
        
            if (is_case_sensitive)
                return QObject::tr("/%1/").arg(qstr);
            else
                return QObject::tr("/%1/i").arg(qstr);
        }
    
        std::string value;
        bool is_case_sensitive;
    };

    class qstring_visitor : public boost::static_visitor<QString>
    {
    public:
        template<class T>
        QString operator()(const T &value) const
        {
            return value.toString();
        }
    
        QString operator()(const std::string &string) const
        {
            return QString("'%1'").arg(QString::fromStdString(string));
        }
    };

    // a single name value. Holds a string or a regular expression
    struct NameValue
    {
        NameVariant value;
        
        QString toString() const
        {
            return boost::apply_visitor( qstring_visitor(), value );
        }
    };

    // an Expression
    struct Expression
    {
        ExpressionVariant value;
        
        QString toString() const
        {
            return boost::apply_visitor( qstring_visitor(), value );
        }
    };

    //part of an expression
    struct ExpressionPart
    {
        ExpressionVariant value;

        QString toString() const
        {
            return boost::apply_visitor( qstring_visitor(), value );
        }
    };

    // a Node contains an array of attributes (which are name-value pairs)
    struct Node
    {
        Expressions values;
        
        QString toString() const
        {
            QStringList lines;
            
            for (const auto value : values)
            {
                lines.append( value.toString() );
            }
            
            return lines.join("; ");
        }
    };

    //a holder for a name of an item
    struct IDName
    {
        IDObject name;
        NameValues values;
        
        QString toString() const
        {
            QStringList lines;
            for (const auto value : values)
            {
                lines.append( value.toString() );
            }
        
            return QObject::tr("%1name %2")
                            .arg( idobject_to_string(name) )
                            .arg(lines.join(","));
        }
    };

    // a binary ID expression, e.g. something AND something
    struct IDBinary
    {
        Expression part0;
        IDOperation operation;
        Expression part1;
        
        QString toString() const
        {
            return QObject::tr("(%1 %2 %3)")
                        .arg(part0.toString())
                        .arg(idoperation_to_string(operation))
                        .arg(part1.toString());
        }
    };
}

BOOST_FUSION_ADAPT_STRUCT( AST::NameValue,
                           (AST::NameVariant,value)
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

BOOST_FUSION_ADAPT_STRUCT( AST::IDBinary,
                           (AST::Expression, part0),
                           (AST::IDOperation, operation),
                           (AST::Expression, part1)
                         )

BOOST_FUSION_ADAPT_STRUCT( AST::Expression,
                           (AST::ExpressionVariant, value)
                         )

BOOST_FUSION_ADAPT_STRUCT( AST::ExpressionPart,
                           (AST::ExpressionVariant, value)
                         )

namespace spirit  = boost::spirit;
namespace qi      = spirit::qi;
namespace phoenix = boost::phoenix;

/** This is the grammar that enables skipping of spaces, newlines and comments */
template<typename IteratorT>
class SkipperGrammar : public qi::grammar<IteratorT>
{
public:
    SkipperGrammar() : SkipperGrammar::base_type( rule )
    {
        lineCommentRule  = qi::lit( "//" ) >>
                           *(qi::char_ -qi::eol) >> 
                           qi::eol;
        blockCommentRule = qi::lit( "/*" ) >> 
                           *(qi::char_ -qi::lit( "*/" ) ) >> 
                           qi::lit( "*/" );
        spaceRule        = qi::space;
        rule             = spaceRule | lineCommentRule | blockCommentRule;
    }

    qi::rule<IteratorT> lineCommentRule;
    qi::rule<IteratorT> blockCommentRule;
    qi::rule<IteratorT> spaceRule;
    qi::rule<IteratorT> rule;
};

/** This is a quoted string grammar that will parse quoted strings and also
    auto-escape characters */
template<typename IteratorT, typename SkipperT>
class ValueGrammar : public qi::grammar<IteratorT, std::string(), SkipperT>
{
public:
    ValueGrammar() : ValueGrammar::base_type( rule, "String" )
    {
        escapedStringRule %= qi::lexeme[
             qi::lit( "'" ) >>
             *( escapeCharSymbols | ( qi::char_ - qi::char_( "'" ) ) ) >
             qi::lit( "'" ) ];
        
        rawStringRule %= qi::lexeme[
                    +( qi::alnum |
                       qi::char_( '.' ) |
                       qi::char_( '/' ) |
                       qi::char_( '_' ) |
                       qi::char_( '-' )
                      ) ];
        
        rule %= rawStringRule | escapedStringRule;
        
        escapeCharSymbols.add( "\\a", '\a' )
                             ( "\\b", '\b' )
                             ( "\\f", '\f' )
                             ( "\\n", '\n' )
                             ( "\\r", '\r' )
                             ( "\\t", '\t' )
                             ( "\\v", '\v' )
                             ( "\\\\", '\\' )
                             ( "\\\'", '\'' )
                             ( "\\\"", '\"' );

        escapedStringRule.name( "Escaped String" );
        rawStringRule.name( "Escaped String" );
        
        escapeCharSymbols.name( "Escaped Chars" );
    }
    
    qi::rule<IteratorT, std::string(), SkipperT>   escapedStringRule;
    qi::rule<IteratorT, std::string(), SkipperT>   rawStringRule;
    qi::rule<IteratorT, std::string(), SkipperT>   rule;
    qi::symbols<const char, const char>            escapeCharSymbols;
};

template<typename IteratorT, typename SkipperT>
class Grammar : public qi::grammar<IteratorT, AST::Node(), SkipperT>
{
public:
    Grammar() : Grammar::base_type( nodeRule, "Node" )
    {
        using qi::lit;
        using qi::lexeme;
        using qi::eps;
        using qi::_1;
        using qi::on_error;
        using qi::fail;
        using namespace qi::labels;
        using qi::as_string;

        using phoenix::construct;
        using phoenix::val;
    
        using boost::spirit::ascii::char_;
    
        nodeRule %= expressionsRule;

        //must be first so that we greedily parse as much as we can
        expressionsRule %= ( expressionRule % qi::lit( ';' ) );

        //must be first so that we greedily parse as much as we can
        expressionRule %= binaryRule2 | binaryRule | expressionPartRule;

        idNameRule  %= name_token >> nameValuesRule;
        
        binaryRule %= (expressionPartRule >> op_token >> expressionPartRule) |
                      ( qi::lit('(') >> binaryRule >> qi::lit(')') );
        binaryRule2 %= binaryRule >> op_token >> binaryRule |
                       binaryRule >> op_token >> expressionPartRule |
                       (qi::lit('(') >> binaryRule2 >> qi::lit(')') );

        expressionPartRule %= ( idNameRule ) |
                              ( qi::lit('(') >> expressionPartRule >> qi::lit(')') );

        name_token.add( "resnam", AST::RESIDUE )
                      ( "resname", AST::RESIDUE )
                      ( "atomnam", AST::ATOM )
                      ( "atomname", AST::ATOM )
                      ( "chainnam", AST::CHAIN )
                      ( "chainname", AST::CHAIN );
        
        op_token.add( "and", AST::ID_AND )
                    ( "AND", AST::ID_AND )
                    ( "or", AST::ID_OR )
                    ( "OR", AST::ID_OR );
        
        nameValuesRule     %= ( nameValueRule % qi::lit( ',' ) );
        
        nameValueRule      %= regExpRule | stringRule;
        
        regExpRule = eps [ _val = AST::RegExpValue() ] >>
                     (
                        lexeme[ "/" >> as_string[+(char_ - "/")][ _val += _1 ] >> "/" ]
                        >> -qi::lit("i")[ _val *= 1 ]
                     )
                     ;
        
        nodeRule.name( "Node" );
        expressionsRule.name( "Expressions" );
        expressionRule.name( "Expression" );
        stringRule.name( "String" );
        regExpRule.name( "RegExp" );
        
        on_error<fail>
        (
            nodeRule
          , std::cout
                << val("Error! Expecting ")
                << _4                               // what failed?
                << val(" here: \"")
                << construct<std::string>(_3, _2)   // iterators to error-pos, end
                << val("\"")
                << std::endl
        );
    }
    
    qi::rule<IteratorT, AST::Node(), SkipperT> nodeRule;
    qi::rule<IteratorT, AST::IDName(), SkipperT> idNameRule;
    qi::rule<IteratorT, AST::IDBinary(), SkipperT> binaryRule;
    qi::rule<IteratorT, AST::IDBinary(), SkipperT> binaryRule2;

    qi::rule<IteratorT, AST::Expressions(), SkipperT> expressionsRule;
    qi::rule<IteratorT, AST::Expression(), SkipperT> expressionRule;
    
    qi::rule<IteratorT, AST::ExpressionPart(), SkipperT> expressionPartRule;
    
    qi::rule<IteratorT, AST::NameValues(), SkipperT> nameValuesRule;
    qi::rule<IteratorT, AST::NameValue(), SkipperT> nameValueRule;
    
    qi::symbols<char,AST::IDObject> name_token;
    qi::symbols<char,AST::IDOperation> op_token;

    ValueGrammar<IteratorT, SkipperT> stringRule;
    qi::rule<IteratorT, AST::RegExpValue(), SkipperT> regExpRule;
};

template<typename IteratorT>
QString toString(IteratorT begin, IteratorT end)
{
    QStringList lines;
    for (; begin != end; ++begin)
    {
        lines.append( QString( *begin ) );
    }
    
    return lines.join("");
}

template<typename IteratorT>
AST::Node parse(const IteratorT & begin, const IteratorT & end)
{
    using LinePosIteratorT  = spirit::line_pos_iterator<IteratorT>;
  
    using SkipperGrammarT   = SkipperGrammar<LinePosIteratorT>;
    using ParserGrammarT    = Grammar<LinePosIteratorT, SkipperGrammarT>;

    SkipperGrammarT  skipper;
    ParserGrammarT   grammar;
    LinePosIteratorT posIterBegin( begin );
    LinePosIteratorT posIterEnd( end );
  
    AST::Node result;

    const bool parseResult = qi::phrase_parse( posIterBegin,
                                               posIterEnd,
                                               grammar,
                                               skipper,
                                               result );
    
    if( not (parseResult && posIterBegin == posIterEnd) )
    {
        QString line = toString( LinePosIteratorT(begin), LinePosIteratorT(end) );
        QString left = toString( posIterBegin, posIterEnd );

        throw SireMol::parse_error( QObject::tr("Failed to parse the selection '%1'. "
          "Successfully parsed the beginning, but failed to parse '%2'")
            .arg(line).arg(left), CODELOC );
    }
    
    return result;
}

int parse_main(const std::string &str)
{
    // Read file contents.
    auto result = parse( str.begin(), str.end() );

    qDebug() << result.toString();
    
    return 0;
}

namespace SireMol
{
    /** Internal function used to parse the passed string and convert
        it into a SelectEngine object */
    boost::shared_ptr<parser::SelectEngine> parse( const QString &str )
    {
        parse_main( str.toStdString() );
        return boost::shared_ptr<parser::SelectEngine>();
    }
}
