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
#include "parser/ast.h"

#include "SireError/errors.h"

#include <QMutex>
#include <QHash>
#include <QDebug>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/support_line_pos_iterator.hpp>

#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/spirit/include/phoenix_object.hpp>

#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/io.hpp>

using namespace SireMol;

namespace spirit  = boost::spirit;
namespace qi      = spirit::qi;
namespace phoenix = boost::phoenix;

////////
//////// implementation of parse_error
////////

const char* parse_error::typeName()
{
    return QMetaType::typeName( qMetaTypeId<parse_error>() );
}

static const RegisterMetaType<parse_error> r_parse;

////////
//////// implementation of the user tokens registry
////////

typedef qi::symbols<char,AST::IDUser> UserTokens;
static UserTokens *_user_tokens = 0;

Q_GLOBAL_STATIC( QMutex, tokensMutex )

/** Get the set of user-supplied tokens */
UserTokens getUserTokens()
{
    QMutexLocker lkr(tokensMutex());
    
    if (_user_tokens == 0)
        _user_tokens = new UserTokens();
    
    return *_user_tokens;
}

/** Clear all of the user-supplied tokens */
static void reset_tokens()
{
    QMutexLocker lkr(tokensMutex());
    
    delete _user_tokens;
    _user_tokens = 0;
}

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

#include "parser/grammar.h" //file containing the actual grammar - separated to ease reading

/** Function that parses the passed string (represented via its iterators) into an AST::Node */
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

/** Function used internally to associate a user token with a user-specified selection */
static void set_token(const std::string &token, const std::string &str)
{
    //first parse this into an AST::Node
    auto node = parse( str.begin(), str.end() );
    
    if (node.values.size() != 1)
        throw SireMol::parse_error( QObject::tr(
            "Cannot set a token based on a multi-line selection!"), CODELOC );
    
    QMutexLocker lkr(tokensMutex());
    
    if (_user_tokens == 0)
        _user_tokens = new UserTokens();
    
    _user_tokens->add(token, AST::IDUser(token,node.values[0].value));
}

/** Function used internally to parse a string into an AST::Node */
static AST::Node parse_main(const std::string &str)
{
    // Read file contents.
    return parse( str.begin(), str.end() );
}

namespace SireMol
{
    namespace parser
    {
        /** Internal function used to parse the passed string and convert
            it into a SelectEngine object */
        SelectEnginePtr parse( const QString &str )
        {
            auto ast = ::parse_main( str.toStdString() );
            
            auto engine = ast.toEngine();
            
            if (engine.get())
                engine = engine->simplify();
            
            return engine;
        }
        
        /** Internal function used to associate a named token with the
            passed user-specified selection string */
        void set_token(const QString &token, const QString &selection)
        {
            ::set_token(token.toStdString(), selection.toStdString());
        }
        
        /** Internal function used to clear all user-specified tokens */
        void reset_tokens()
        {
            ::reset_tokens();
        }
    }
}
