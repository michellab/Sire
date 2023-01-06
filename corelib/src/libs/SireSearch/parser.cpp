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

#define BOOST_SPIRIT_USE_PHOENIX_V3
#define BOOST_SPIRIT_UNICODE

#include "parser.h"
#include "ast.h"

#include "SireError/errors.h"

#include "SireMol/core.h"

#include "SireSearch/ast.h"

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

using namespace SireSearch;
using namespace SireMol;

namespace spirit  = boost::spirit;
namespace qi      = spirit::qi;
namespace unicode = boost::spirit::unicode;
namespace phoenix = boost::phoenix;

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

#include "grammar.h" //file containing the actual grammar - separated to ease reading

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

        if (line == left)
            throw SireMol::parse_error( QObject::tr(
                "Failed to parse any of the selection '%1'."
            ).arg(line), CODELOC);
        else
            throw SireMol::parse_error( QObject::tr(
              "Failed to parse the selection '%1'. "
              "Successfully parsed the beginning, but failed to parse '%2'.")
                .arg(line).arg(left), CODELOC );
    }

    return result;
}

/** Function used internally to associate a user token with a user-specified selection */
static void set_token(const std::string &token, const std::string &str)
{
    QString str_token = QString::fromStdString(token);
    str_token = str_token.simplified();

    // token names cannot contain spaces
    for (const auto &token_char : str_token)
    {
        if (token_char.isSpace())
        {
            throw SireError::invalid_key(QObject::tr(
                "You cannot use '%1' as the name of a search token as it "
                "contains a space (or space-like) character. Please replace "
                "all spaces with, e.g. underscores or hyphens.")
                    .arg(str_token), CODELOC);
        }
    }

    // they must also have a size...
    if (str_token.isEmpty())
    {
        throw SireError::invalid_key(QObject::tr(
            "You can't use an empty string as the name of a search token."),
                CODELOC);
    }

    // they also can't start with a number (or maths) character
    if (not str_token[0].isLetter())
    {
        throw SireError::invalid_key(QObject::tr(
            "You cannot use '%1' as the name of a search token as it "
            "starts with a non-letter character. Please "
            "use a name that starts with a letter.")
                .arg(str_token), CODELOC);
    }

    auto cleaned_token = str_token.toStdString();

    //first check that we *can't* parse the name of the token. We should
    //not be able to set token names that are valid search terms
    bool token_is_parseable = false;
    QString part_parseable;

    try
    {
        parse(cleaned_token.begin(), cleaned_token.end());
        token_is_parseable = true;
    }
    catch(const SireMol::parse_error &e)
    {
        // we want an exception to be raised as this
        // shows that the token is not parseable

        // but we need it to be entirely unparseable, else
        // the token may confuse the parse engine. Tokens that
        // start with 'protein' or 'resnum' confuse it

        // This is a really annoying restriction...
        if (e.why().indexOf("Failed to parse any of the selection") == -1)
        {
            token_is_parseable = true;
            part_parseable = e.why();
        }
    }

    if (token_is_parseable)
    {
        if (part_parseable.isEmpty())
        {
            throw SireError::invalid_key(QObject::tr(
                "You cannot use '%1' as the name of a search token as this is, "
                "itself, a valid search phrase. Please choose a name that is not "
                "an existing search token.").arg(str_token), CODELOC);
        }
        else
        {
            throw SireError::invalid_key(QObject::tr(
                "You cannot use '%1' as the name of a search token as it starts "
                "with an existing search token. This would confuse the parser. "
                "Please choose a name that does not start with any existing "
                "search term. To help, here is the partial-parse message: %2")
                    .arg(str_token).arg(part_parseable), CODELOC);
        }
    }

    //first parse this into an AST::Node
    auto node = parse( str.begin(), str.end() );

    QMutexLocker lkr(tokensMutex());

    if (_user_tokens == 0)
        _user_tokens = new UserTokens();

    _user_tokens->add(cleaned_token, AST::IDUser(token,node.values.value));
}

static bool has_token(const std::string &token)
{
    QMutexLocker lkr(tokensMutex());

    if (_user_tokens == 0)
        return false;

    return _user_tokens->find(token) != 0;
}

static QString get_token(const std::string &token)
{
    QMutexLocker lkr(tokensMutex());

    if (_user_tokens == 0)
        throw SireError::invalid_key( QObject::tr(
            "There is no user token '%1'").arg(QString::fromStdString(token)), CODELOC);

    auto p = _user_tokens->find(token);

    if (p == 0)
        throw SireError::invalid_key( QObject::tr(
            "There is no user token '%1'").arg(QString::fromStdString(token)), CODELOC);

    return expression_to_string(p->value);
}

static void delete_token(const std::string &token)
{
    QMutexLocker lkr(tokensMutex());

    if (_user_tokens == 0)
        return;

    auto p = _user_tokens->find(token);

    if (p)
        _user_tokens->remove(token);
}

/** Function used internally to parse a string into an AST::Node */
static AST::Node parse_main(const std::string &str)
{
    // Read file contents.
    return parse( str.begin(), str.end() );
}

namespace SireSearch
{
    namespace parser
    {
        SearchParser::SearchParser() : SireMol::parser::Parser()
        {}

        SearchParser::~SearchParser()
        {}

        void SearchParser::install()
        {
            SireMol::parser::Parser::install_parser(new SearchParser());
        }

        void SearchParser::set_token(const QString &token, const QString &selection)
        {
            ::set_token(token.toStdString(), selection.toStdString());
        }

        bool SearchParser::has_token(const QString &token)
        {
            return ::has_token(token.toStdString());
        }

        QString SearchParser::get_token(const QString &token)
        {
            return ::get_token(token.toStdString());
        }

        void SearchParser::delete_token(const QString &token)
        {
            ::delete_token(token.toStdString());
        }

        void SearchParser::delete_all_tokens()
        {
            ::reset_tokens();
        }

        void SearchParser::reset_tokens()
        {
            ::reset_tokens();
        }

        SireMol::parser::SelectEnginePtr SearchParser::parse(const QString &str)
        {
            auto ast = ::parse_main( str.toStdString() );

            auto engine = ast.toEngine();

            if (engine.get())
                engine = engine->simplify();

            return engine;
        }
    }
}
