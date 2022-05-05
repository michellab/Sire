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

#include "SireError/errors.h"

namespace SireMol
{
    namespace parser
    {
        Parser* Parser::global_parser(0);

        void Parser::install_parser(Parser *parser)
        {
            if (global_parser != 0)
                delete global_parser;

            global_parser = parser;
        }

        Parser& Parser::globalParser()
        {
            if (global_parser == 0)
                throw SireError::program_bug( QObject::tr(
                    "No parser has been loaded!"), CODELOC );

            return *global_parser;
        }

        void set_token(const QString &token, const QString &selection)
        {
            Parser::globalParser().set_token(token, selection);
        }

        void reset_tokens()
        {
            Parser::globalParser().reset_tokens();
        }

        /** Internal function used to parse the passed string and convert
            it into a SelectEngine object */
        SelectEnginePtr parse( const QString &str )
        {
            return Parser::globalParser().parse(str);
        }
    }
}
