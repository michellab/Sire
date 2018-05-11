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

// include boost::spirit::qi for parsing
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/support_line_pos_iterator.hpp>

// A lot of the below code is heavily inspired by
// https://medium.com/@alinakipoglu/parsing-with-spirit-qi-fcaeaf4357b3

///////////
/////////// First define the objects used for the abstract syntax tree
///////////

namespace AST
{
    struct Value;      // holder for a generic value
    struct NameValue;  // combination of a name with a value
    struct Node;       // a node in the tree
    struct Array;      // an array of nodes

    // a Variant can hold a single node, an array of values or a string
    using Variant = boost::variant<boost::recursive_wrapper<Node>,
                                   boost::recursive_wrapper<Array>,
                                   std::string>;
    
    // the attributes are the set of name-value pairs
    using Attributes = std::vector<NameValue>;
    
    // an array of generic value objects
    using Values = std::vector<Value>;

    // a value holds a single variant (which can be a Node, Array of Nodes or a string)
    struct Value
    {
        Variant value;
    };

    // a Node contains an array of attributes (which are name-value pairs)
    struct Node
    {
        Attributes attributes;
    };

    // an Array contains an array of generic value objects
    struct Array
    {
        Values values;
    };

    // a NameValue provides both a name and a value
    struct NameValue
    {
        std::string name;
        Value value;
    };

    void to_string(const NameValue &val);
    void to_string(const Node &node);
    void to_string(const Value &val);

    class print_visitor : public boost::static_visitor<int>
    {
    public:
        int operator()(const Node &node) const
        {
            std::cout << "node--\n";
            to_string(node);
            return 0;
        }
        
        int operator()(const Array &array) const
        {
            std::cout << "Array--\n";
            
            int i = 0;
            for (auto value : array.values)
            {
                std::cout << "value" << i << std::endl;
                to_string(value);
                std::cout << "\n";
                i += 1;
            }
            
            return 0;
        }
        
        int operator()(const std::string &str) const
        {
            std::cout << "value: " << str << " ";
            return 0;
        }
    };

    void to_string(const Value &value)
    {
        boost::apply_visitor( print_visitor(), value.value );
    }

    void to_string(const Node &node)
    {
        std::cout << "attributes\n";
        int i = 0;
        for (auto attribute : node.attributes)
        {
            std::cout << "attribute " << i << std::endl;
            std::cout << attribute.name << " = ";
            to_string( attribute.value );
            std::cout << "\n";
            i += 1;
        }
    }
}

BOOST_FUSION_ADAPT_STRUCT( AST::Value,
                           (AST::Variant,value)
                         )

BOOST_FUSION_ADAPT_STRUCT( AST::Node,
                           (AST::Attributes,attributes)
                         )

BOOST_FUSION_ADAPT_STRUCT( AST::Array,
                           (AST::Values,values)
                         )

BOOST_FUSION_ADAPT_STRUCT( AST::NameValue,
                           (std::string, name),
                           (AST::Value, value)
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
             qi::lit( '"' ) >>
             *( escapeCharSymbols | ( qi::char_ - qi::char_( '"' ) ) ) >
             qi::lit( '"') ];
        
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
    Grammar() : Grammar::base_type( objectRule, "Object" )
    {
        objectRule %= qi::lit( '{' ) >> -attributesRule > qi::lit( '}' );
        
        arrayRule %= qi::lit( '(' ) >>
                       -valuesRule >>
                       -qi::lit( ',' ) >
                        qi::lit( ')' );
        
        attributesRule %= attributeRule >>
                          *( qi::lit( ';' ) >> attributeRule ) >> 
                          -qi::lit( ';' );
        
        attributeRule  %= nameRule > qi::lit( '=' ) > valueRule;
        
        nameRule       %= qi::lexeme[ +( qi::alnum | qi::char_( '_' ) ) ];
        
        valuesRule     %= ( objectRule % qi::lit( ',' ) ) |
                          ( arrayRule % qi::lit( ',' ) ) |
                          ( stringRule % qi::lit( ',' ) );
        
        valueRule      %= objectRule | arrayRule | stringRule;
        
        objectRule.name( "Object" );
        arrayRule.name( "Array" );
        attributesRule.name( "Attributes" );
        attributeRule.name( "Attribute" );
        nameRule.name( "Name" );
        valuesRule.name( "Values" );
        valueRule.name( "Value" );
        stringRule.name( "String" );
    }
    
    qi::rule<IteratorT, AST::Node(), SkipperT> objectRule;
    qi::rule<IteratorT, AST::Array(), SkipperT> arrayRule;
    qi::rule<IteratorT, AST::Attributes(), SkipperT> attributesRule;
    qi::rule<IteratorT, AST::NameValue(), SkipperT> attributeRule;
    qi::rule<IteratorT, std::string(), SkipperT> nameRule;
    qi::rule<IteratorT, AST::Values(), SkipperT> valuesRule;
    qi::rule<IteratorT, AST::Value(), SkipperT> valueRule;
    ValueGrammar<IteratorT, SkipperT> stringRule;
};

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
  
    try
    {
        AST::Node result;
  
        const bool parseResult = qi::phrase_parse( posIterBegin,
                                                   posIterEnd,
                                                   grammar,
                                                   skipper,
                                                   result );
        
        if( parseResult && posIterBegin == posIterEnd )
        {
            qDebug() << "PARSED USING AST!";
            return result;
        }
        else
        {
            throw std::runtime_error{ "Parser Failed" };
        }
    }
    catch ( ... )
    {
        throw std::runtime_error{ "Parser Failed" };
    }
    
    return AST::Node{};
}

int parse_main(const std::string &str)
{
    // Read file contents.
    AST::Node result;

    try
    {
        result = parse( str.begin(), str.end() );
    }
    catch(...)
    {
        qDebug() << "ERROR IN PARSING";
    }

    AST::to_string(result);
    
    return 0;
}

namespace SireMol
{
    /** Internal function used to parse the passed string and convert
        it into a SelectEngine object */
    boost::shared_ptr<parser::SelectEngine> parse( const QString &str )
    {
        const auto cstr = str.toStdString();

        parse_main(cstr);
        
        return boost::shared_ptr<parser::SelectEngine>();
    }
}
