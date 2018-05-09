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

#include "idengine.h"

#include "tostring.h"

using namespace parser_idengine;

////////
//////// Implementation of the IDEngine that actually performs the selection
////////

IDEngine::IDEngine() : SelectEngine()
{}

IDEngine::~IDEngine()
{}

IDEngine& IDEngine::operator+=(const QPair<int,int> &ids)
{
    obj = IDEngine::Obj(ids.first);
    typ = IDEngine::Typ(ids.second);
    
    return *this;
}

IDEngine& IDEngine::operator+=(const StringsOrRegexps &names)
{
    idnams = names;
    return *this;
}

IDEngine& IDEngine::operator+=(const NumbersOrRanges &numbers)
{
    idnums = numbers;
    return *this;
}

QString IDEngine::toString() const
{
    QString objstr, typstr;
    
    switch (obj)
    {
    case ATOM:
        objstr = "atom";
        break;
    case CUTGROUP:
        objstr = "cg";
        break;
    case RESIDUE:
        objstr = "res";
        break;
    case CHAIN:
        objstr = "chain";
        break;
    case SEGMENT:
        objstr = "seg";
        break;
    case MOLECULE:
        objstr = "mol";
        break;
    }
    
    switch (typ)
    {
    case NAME:
        typstr = QString("name = %1").arg(idnams.strings.join(","));
        break;
    case NUMBER:
        typstr = QString("num = %1").arg(Sire::toString(idnums.numbers));
        break;
    case INDEX:
        typstr = QString("idx = %1").arg(Sire::toString(idnums.numbers));
    }
    
    return QString("%1%2").arg(objstr,typstr);
}

/////////
///////// Implementation of StringsOrRegexps
/////////

StringsOrRegexps::StringsOrRegexps()
{}

StringsOrRegexps::~StringsOrRegexps()
{}

StringsOrRegexps& StringsOrRegexps::operator+=(const QString &s)
{
    strings.append(s);
    qSort(strings);
    return *this;
}

StringsOrRegexps& StringsOrRegexps::operator+=(const std::wstring &s)
{
    return this->operator+=( QString::fromStdWString(s) );
}

StringsOrRegexps& StringsOrRegexps::operator+=(const std::string &s)
{
    return this->operator+=( QString::fromStdString(s) );
}

StringsOrRegexps& StringsOrRegexps::operator+=(const std::vector<char> &s)
{
    return this->operator+=( QString::fromUtf8(s.data(), s.size()) );
}

/////////
///////// Implementation of IDRange
/////////

IDRange::IDRange() : start(0), end(0), step(1), c(0)
{}

IDRange::~IDRange()
{}

IDRange& IDRange::operator+=(int val)
{
    if (c == 0)
    {
        start = val;
        end = val;
        step = 1;
        c += 1;
    }
    else if (c == 1)
    {
        end = val;
        step = 1;
        
        c += 1;
    }
    else if (c == 2)
    {
        step = val;
        
        if (step == 0)
            step = 1;
        else if (step < 0)
            step = -val;
        
        c += 1;
    }
    else
        qDebug() << "PARSE RANGE ERROR!";
    
    return *this;
}

/////////
///////// Implementation of NumbersOrRanges
/////////

NumbersOrRanges::NumbersOrRanges()
{}

NumbersOrRanges::~NumbersOrRanges()
{}

NumbersOrRanges& NumbersOrRanges::operator+=(const IDRange &range)
{
    if (range.start == range.end)
        numbers.append(range.start);
    else
    {
        if (range.end > range.start)
        {
            for (int i=range.start; i<=range.end; i+=range.step)
            {
                numbers.append(i);
            }
        }
        else
        {
            for (int i=range.start; i>=range.end; i-=range.step)
            {
                numbers.append(i);
            }
        }
    }

    qSort(numbers);
    return *this;
}

/////////
///////// Implementation of the spirit parser that describes the grammar
/////////

struct name_objects_ : qi::symbols< char, QPair<int,int> >
{
    name_objects_()
    {
        add
            ("atomname"    ,    QPair<int,int>(IDEngine::ATOM, IDEngine::NAME))
            ("groupname"    ,    QPair<int,int>(IDEngine::CUTGROUP, IDEngine::NAME))
            ("resname"    ,    QPair<int,int>(IDEngine::RESIDUE, IDEngine::NAME))
            ("chainname"    ,    QPair<int,int>(IDEngine::CHAIN, IDEngine::NAME))
            ("segname"    ,    QPair<int,int>(IDEngine::SEGMENT, IDEngine::NAME))
            ("molname"    ,    QPair<int,int>(IDEngine::MOLECULE, IDEngine::NAME))
            ("atomnam"    ,    QPair<int,int>(IDEngine::ATOM, IDEngine::NAME))
            ("groupnam"    ,    QPair<int,int>(IDEngine::CUTGROUP, IDEngine::NAME))
            ("resnam"    ,    QPair<int,int>(IDEngine::RESIDUE, IDEngine::NAME))
            ("chainnam"    ,    QPair<int,int>(IDEngine::CHAIN, IDEngine::NAME))
            ("segnam"    ,    QPair<int,int>(IDEngine::SEGMENT, IDEngine::NAME))
            ("molnam"    ,    QPair<int,int>(IDEngine::MOLECULE, IDEngine::NAME))
        ;
    }
} name_objects;

struct num_objects_ : qi::symbols< char, QPair<int,int> >
{
    num_objects_()
    {
        add
            ("atomnum"    ,    QPair<int,int>(IDEngine::ATOM, IDEngine::NUMBER))
            ("resnum"    ,    QPair<int,int>(IDEngine::RESIDUE, IDEngine::NUMBER))
            ("molnum"    ,    QPair<int,int>(IDEngine::MOLECULE, IDEngine::NUMBER))
            ("atomidx"    ,    QPair<int,int>(IDEngine::ATOM, IDEngine::INDEX))
            ("groupidx"    ,    QPair<int,int>(IDEngine::CUTGROUP, IDEngine::INDEX))
            ("residx"    ,    QPair<int,int>(IDEngine::RESIDUE, IDEngine::INDEX))
            ("chainidx"    ,    QPair<int,int>(IDEngine::CHAIN, IDEngine::INDEX))
            ("segidx"    ,    QPair<int,int>(IDEngine::SEGMENT, IDEngine::INDEX))
            ("molidx"    ,    QPair<int,int>(IDEngine::MOLECULE, IDEngine::INDEX))
            ("atomnumber"    ,    QPair<int,int>(IDEngine::ATOM, IDEngine::NUMBER))
            ("resnumber"    ,    QPair<int,int>(IDEngine::RESIDUE, IDEngine::NUMBER))
            ("molnumber"    ,    QPair<int,int>(IDEngine::MOLECULE, IDEngine::NUMBER))
            ("atomindex"    ,    QPair<int,int>(IDEngine::ATOM, IDEngine::INDEX))
            ("groupindex"    ,    QPair<int,int>(IDEngine::CUTGROUP, IDEngine::INDEX))
            ("resindex"    ,    QPair<int,int>(IDEngine::RESIDUE, IDEngine::INDEX))
            ("chainindex"    ,    QPair<int,int>(IDEngine::CHAIN, IDEngine::INDEX))
            ("segindex"    ,    QPair<int,int>(IDEngine::SEGMENT, IDEngine::INDEX))
            ("molindex"    ,    QPair<int,int>(IDEngine::MOLECULE, IDEngine::INDEX))
        ;
    }
} num_objects;

/** The actual grammar for the IDEngine parser */
idengine_parser::idengine_parser() : idengine_parser::base_type(start)
{
    using qi::int_;
    using qi::_val;
    using qi::eps;
    using qi::lit;
    using qi::double_;
    using qi::_1;
    using qi::_2;
    using qi::lexeme;
    using ascii::char_;

    strings_or_regexps = eps [ _val = StringsOrRegexps() ] >>
        (
            lexeme[+char_][ _val += _1 ]
        )
        ;
    
    id_range = eps [ _val = IDRange() ] >>
        (
            int_[ _val += _1 ] >>
            *( ':' >> int_[ _val += _1 ] )
        )
        ;
    
    numbers_or_ranges = eps [ _val = NumbersOrRanges() ] >>
        (
            (id_range[ _val += _1 ]) >>
            *( ',' >> id_range[ _val += _1 ] )
        )
        ;

    start = eps [ _val = IDEngine() ] >>
        (
            name_objects[ _val += _1 ] >>
            strings_or_regexps[ _val += _1 ]
        ) |
        (
            num_objects[ _val += _1 ] >>
            numbers_or_ranges[ _val += _1 ]
        )
        ;
}
