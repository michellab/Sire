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
using namespace SireMol::parser;
using namespace SireBase;

//the parsers used to select by atom/chain/etc ID (name, number, index)
#include "parser/idengine.h"

// word count lexer
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/lex_lexertl.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_statement.hpp>
#include <boost/spirit/include/phoenix_container.hpp>

using namespace boost::spirit;
using namespace boost::spirit::ascii;

enum ID_TOKENS { ID_ANY = 1,
                 ID_OPEN_BRACKET = 2,
                 ID_CLOSE_BRACKET = 3,
                 ID_SINGLE_QUOTE = 4,
                 ID_SPACE = 5 };

template <typename Lexer>
struct word_count_tokens : lex::lexer<Lexer>
{
    word_count_tokens()
    {
        // define patterns (lexer macros) to be used during token definition 
        // below
        this->self.add_pattern("WORD", "[^ \t\n]+");

        // define tokens and associate them with the lexer
        word = "{WORD}";    // reference the pattern 'WORD' as defined above

        // this lexer will recognize 3 token types: words, newlines, and 
        // everything else
        this->self.add
            ('(', ID_OPEN_BRACKET)
            (')', ID_CLOSE_BRACKET)
            ('\'', ID_SINGLE_QUOTE)
            (word)
            ('\n', ID_SPACE)
            (".", ID_ANY)
        ;
    }

    // the token 'word' exposes the matched string as its parser attribute
    lex::token_def<std::string> word;
};

class ExpressionStack
{
public:
    ExpressionStack() : depth(0), is_quote_open(false)
    {}
    
    ~ExpressionStack()
    {}
    
    void pushText(const QString &text)
    {
        if (expressions.isEmpty())
        {
            expressions.push_back( QString() );
        }

        expressions.last() += text;
    }
    
    void pushBracket()
    {
        expressions.push_back( QString() );
        depth += 1;
    }
    
    void popBracket()
    {
        if (depth < 1)
            qDebug() << "SOMETHING WRONG WITH BRACKET DEPTH";
        
        QString line = expressions.last();
        expressions.pop_back();
        qDebug() << "NOW TIME TO PARSE" << line;
        depth -= 1;
    }
    
    ExpressionStack& operator+=(ID_TOKENS token)
    {
        if (is_quote_open)
        {
            switch(token)
            {
            case ID_OPEN_BRACKET:
                pushText("(");
                break;
            case ID_CLOSE_BRACKET:
                pushText(")");
                break;
            case ID_SINGLE_QUOTE:
                pushText("'");
                is_quote_open = false;
                break;
            case ID_ANY:
            case ID_SPACE:
                pushText(" ");
                break;
            };
        }
        else
        {
            switch(token)
            {
            case ID_OPEN_BRACKET:
                this->pushBracket();
                break;
            case ID_CLOSE_BRACKET:
                this->popBracket();
                break;
            case ID_SINGLE_QUOTE:
                pushText("'");
                is_quote_open = true;
                break;
            case ID_ANY:
            case ID_SPACE:
                pushText(" ");
                break;
            };
        }
        return *this;
    }
    
    ExpressionStack& operator+=(const std::string &word)
    {
        pushText( QString::fromStdString(word) );
        return *this;
    }
    
    QVector<QString> expressions;
    int depth;
    bool is_quote_open;
};

template <typename Iterator>
struct word_count_grammar : qi::grammar<Iterator>
{
    template <typename TokenDef>
    word_count_grammar(TokenDef const& tok)
      : word_count_grammar::base_type(start)
    {
        using boost::phoenix::ref;
        using boost::phoenix::size;

        start =  *(   qi::token(ID_OPEN_BRACKET)  [ref(stack) += ID_OPEN_BRACKET]
                  |   qi::token(ID_CLOSE_BRACKET) [ref(stack) += ID_CLOSE_BRACKET]
                  |   qi::token(ID_SINGLE_QUOTE)  [ref(stack) += ID_SINGLE_QUOTE]
                  |   tok.word                    [ref(stack) += _1]
                  |   lit('\n')                   [ref(stack) += ID_SPACE]
                  |   qi::token(ID_ANY)           [ref(stack) += ID_ANY]
                  |   qi::token(ID_SPACE)         [ref(stack) += ID_SPACE]
                  )
              ;
    }

    ExpressionStack stack;
    qi::rule<Iterator> start;
};

namespace SireMol
{
    /** Internal function used to parse the passed string and convert
        it into a SelectEngine object */
    boost::shared_ptr<SelectEngine> parse( const QString &str )
    {
        //get a standard wstring from the QString
        const auto s = str.toStdWString();
     
        typedef lex::lexertl::token<
                char const*, boost::mpl::vector<std::string>
            > token_type;
        
        typedef lex::lexertl::lexer<token_type> lexer_type;
        
        typedef word_count_tokens<lexer_type>::iterator_type iterator_type;
        
        // now we use the types defined above to create the lexer and grammar
        // object instances needed to invoke the parsing process
        word_count_tokens<lexer_type> word_count;          // Our lexer
        word_count_grammar<iterator_type> gw (word_count);  // Our parser
        
        const auto cstr = str.toStdString();

        char const* first = cstr.c_str();
        char const* last = &first[cstr.size()];
        
        bool r = lex::tokenize_and_parse(first, last, word_count, gw);
        
        // print results
        if (r) {
            std::cout << "successs\n";
        }
        else {
            std::string rest(first, last);
            std::cout << "Lexical analysis failed\n" << "stopped at: \""
                      << rest << "\"\n";
        }
    
        parser_idengine::IDEngine eng;
        parser_idengine::idengine_parser g;
        
        r = phrase_parse(s.begin(), s.end(), g, boost::spirit::ascii::space, eng);

        if (r)
        {
            qDebug() << "PARSING SUCCESSFUL";
            qDebug() << eng.toString();
        }
        else
        {
            qDebug() << "PARSING FAILED!";
        }
        
        return boost::shared_ptr<SelectEngine>(new parser_idengine::IDEngine(eng));
    }
}
