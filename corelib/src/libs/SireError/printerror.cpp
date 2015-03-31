/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2008  Christopher Woods
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

#include <QTextStream>

#include "printerror.h"
#include "SireError/exception.h"

namespace SireError
{

void printAaargh(QTextStream &ts)
{
    ts << " ^   ^  \n"
       << " O   O  \n"
       << "   \"   \n"
       << " /\"\"\"\\   AAAAAARRRRGGGHHHH!!!!!!!\n"
       << " \\___/\n\n";
}

/** Print the error 'text' to standard output */
void SIREERROR_EXPORT printError(const QString &text)
{
    QTextStream ts(stdout);
    
    ts << "\n**********************************************************\n"
       << "Something went wrong with the simulation. Here's the error\n"
       << "**********************************************************\n\n"
       << text
       << "\n\nSorry - your simulation terminated with an error. Scroll back up\n"
       << "to the top of the error to work out what has gone wrong.\n\n";
    
    printAaargh(ts);
}

/** Print the error 'e' to standard output */
void SIREERROR_EXPORT printError(const SireError::exception &e)
{
    SireError::printError( e.toString() );
}

}
