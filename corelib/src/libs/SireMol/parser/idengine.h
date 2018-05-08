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

#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_object.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/io.hpp>

#include "SireMol/select.h"

namespace parser_idengine
{

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;
namespace phoenix = boost::phoenix;

using SireMol::parser::SelectEngine;

/** Internal class used to parse strings or regexps */
class StringsOrRegexps
{
public:
    StringsOrRegexps();
    ~StringsOrRegexps();
    
    StringsOrRegexps& operator+=(const std::string &s);
    StringsOrRegexps& operator+=(const std::wstring &s);
    
    QStringList strings;
    QList<QRegExp> regexps;
};

/** Internal class used to parse a range of integers */
class IDRange
{
public:
    IDRange();
    ~IDRange();
    
    IDRange& operator+=(int val);
    
    int start;
    int end;
    int step;

private:
    int c;
};

/** Internal class used to parse numbers or ranges */
class NumbersOrRanges
{
public:
    NumbersOrRanges();
    ~NumbersOrRanges();
    
    NumbersOrRanges& operator+=(const IDRange &range);
    
    QList<int> numbers;
};

/** Internal class providing the SelectEngine for objects 
    based on their IDs (names, numbers, indicies etc.) 
    
    @author Christopher Woods
*/
class IDEngine : public SelectEngine
{
public:
    enum Obj { ATOM = 1,
               CUTGROUP = 2,
               RESIDUE = 3,
               CHAIN = 4,
               SEGMENT = 5,
               MOLECULE = 6 };
    
    enum Typ { NAME = 1,
               NUMBER = 2,
               INDEX = 3 };

    IDEngine();
    ~IDEngine();
    
    IDEngine& operator+=(const QPair<int,int> &ids);
    IDEngine& operator+=(const StringsOrRegexps &names);
    IDEngine& operator+=(const NumbersOrRanges &numbers);
    
    QString toString() const;
    
    Obj obj;
    Typ typ;
    
    StringsOrRegexps idnams;
    NumbersOrRanges idnums;
};

struct idengine_parser : qi::grammar<std::wstring::const_iterator, IDEngine(), ascii::space_type>
{
    typedef std::wstring::const_iterator Iterator;

    idengine_parser();
    
    qi::rule<Iterator, IDEngine(), ascii::space_type> start;
    qi::rule<Iterator, IDRange(), ascii::space_type> id_range;
    qi::rule<Iterator, StringsOrRegexps(), ascii::space_type> strings_or_regexps;
    qi::rule<Iterator, NumbersOrRanges(), ascii::space_type> numbers_or_ranges;
};

}

#endif
