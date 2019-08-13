/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2006  Christopher Woods
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

#ifndef SIREBASE_INCREMENT_H
#define SIREBASE_INCREMENT_H

#include <QString>
#include <QRegExp>

SIRE_BEGIN_HEADER

namespace SireBase
{

QString increment(QString name);

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** A small function that increments a name, e.g. the name foobar_1 will be incremented
    to foobar_2, or foobar will be incremented to foobar_1 */
SIRE_ALWAYS_INLINE QString increment(QString name)
{
    QRegExp nummatch("_(\\d+)$");
    
    int pos = nummatch.indexIn(name);
    if (pos != -1)
    {
        QString numstring = nummatch.cap(1);
        int lgth = numstring.size();
        int num = numstring.toInt();
        name.remove(pos,name.length());
        num++;
        return QString("%1_%2").arg(name).arg(QString::number(num),lgth,'0');
    }
    else
      return QString("%1_1").arg(name);
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

SIRE_EXPOSE_FUNCTION( SireBase::increment );

SIRE_END_HEADER

#endif
