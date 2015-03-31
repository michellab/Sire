/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2009  Christopher Woods
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

#include "sharestrings.h"

#include <QSet>
#include <QMutex>

Q_GLOBAL_STATIC( QMutex, stringRegistryMutex )

static QSet<QString> *string_registry(0);

namespace SireStream
{

/** This function adds the string 'string' to shared storage, and returns 
    a copy of the shared-stored value. Use this function to ensure that
    there is only one copy of duplicate strings */
QString SIRESTREAM_EXPORT shareString(const QString &string)
{
    QMutex *mutex = stringRegistryMutex();
    
    if (not mutex)
        return string;
        
    QMutexLocker lkr(mutex);
    
    if (not string_registry)
    {
        string_registry = new QSet<QString>();
    }
    
    QSet<QString>::const_iterator it = string_registry->constFind(string);
    
    if (it != string_registry->constEnd())
        return *it;
    else
    {
        string_registry->insert(string);
        return string;
    }
}

/** This function adds all of the strings in 'strings' to shared storage */
void SIRESTREAM_EXPORT shareStrings(QStringList &strings)
{
    for (QStringList::iterator it = strings.begin();
         it != strings.end();
         ++it)
    {
        *it = shareString(*it);
    }
}

}
