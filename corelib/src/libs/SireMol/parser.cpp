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

/*
class NGin
{
public:
    NGin()
    {}
    
    ~NGin()
    {}
    
    template<class T>
    NGin& operator+=(const T &ngin)
    {
        engine.reset( new T(ngin) );
        return *this;
    }
    
    boost::shared_ptr<SelectEngine> engine;
};

struct select_parser : qi::grammar<std::wstring::const_iterator, NGin(), ascii::space_type>
{
    typedef std::wstring::const_iterator Iterator;

    select_parser()
    {
        parser_
    }
    
    qi::rule<Iterator, NGin(), ascii::space_type> start;
};
*/

namespace SireMol
{
    /** Internal function used to parse the passed string and convert
        it into a SelectEngine object */
    boost::shared_ptr<SelectEngine> parse( const QString &str )
    {
        //get a standard wstring from the QString
        const auto s = str.toStdWString();
        
        parser_idengine::IDEngine eng;
        parser_idengine::idengine_parser g;
        
        bool r = phrase_parse(s.begin(), s.end(), g, boost::spirit::ascii::space, eng);

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
