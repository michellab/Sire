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

#ifndef SIREMOL_PARSER_H
#define SIREMOL_PARSER_H

#include "select.h"
#include "SireMol/errors.h"

SIRE_BEGIN_HEADER

namespace SireMol
{

/** This exception is thrown when there was an error parsing a selection

    @author Christopher Woods
*/
class SIREMOL_EXPORT parse_error : public siremol_error
{
public:
    parse_error() : siremol_error()
    {}

    parse_error(QString err, QString place = QString::null)
              : siremol_error(err,place)
    {}

    parse_error(const parse_error &other) : siremol_error(other)
    {}

    ~parse_error() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return parse_error::typeName();
    }
    
    void throwSelf() const
    {
        throw parse_error(*this);
    }
};

namespace parser
{
    using SireMol::parser::SelectEnginePtr;

    void set_token(const QString &token, const QString &selection);
    void reset_tokens();

    SelectEnginePtr parse(const QString &str);
}

}

Q_DECLARE_METATYPE(SireMol::parse_error)

SIRE_END_HEADER

#endif
