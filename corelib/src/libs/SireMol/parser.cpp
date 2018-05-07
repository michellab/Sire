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

#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/bind.hpp>

#include "parser.h"

#include <QDebug>

using namespace SireMol;

namespace parser
{
    namespace qi = boost::spirit::qi;
    using boost::spirit::qi::int_;

    void print(const int &i)
    {
        qDebug() << "FOUND" << i;
    }
    
    
    /** Internal function that uses boost::spirit to do the parsing */
    QList<SelectPtr> parse_select(const QString &str)
    {
        //get a standard wstring from the QString
        const auto wstring = str.toStdWString();
        
        boost::spirit::qi::parse(wstring.begin(), wstring.end(),
                    '{' >> int_[&print] >> '}');
        
        return QList<SelectPtr>();
    }
}

namespace SireMol
{
    /** Internal function used to parse the passed string and convert
        it into a set of Select objects */
    QList<SelectPtr> parse( const QString &str )
    {
        qDebug() << "Parsing" << str;
        return parser::parse_select(str);
    }
}
