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

#ifndef SIRECAS_FUNCTIONS_H
#define SIRECAS_FUNCTIONS_H

#include <QSet>

#include "function.h"

#include "sireglobal.h"

SIRE_BEGIN_HEADER

namespace SireCAS
{

/** Trival derivation of QSet<Function> that adds a constructor that
    automatically adds the passed Function. There are additional functions
    that retrieve all of the function names. 
    
    @author Christopher Woods
*/
class SIRECAS_EXPORT Functions : public QSet<Function>
{
public:
    Functions() : QSet<Function>()
    {}
    
    Functions(const Function &func) : QSet<Function>()
    {
        this->insert(func);
    }
    
    Functions(const QSet<Function> &other) : QSet<Function>(other)
    {}
    
    Functions(const QList<Function> &other) : QSet<Function>()
    {
        int n = other.count();
        for (int i=0; i<n; ++i)
            insert( other.at(i) );
    }
    
    ~Functions()
    {}
    
    void insert(const Function &func)
    {
        QSet<Function>::insert(func);
    }
    
    void insert(const Functions &funcs)
    {
        for (Functions::const_iterator it = funcs.begin();
             it != funcs.end();
             ++it)
        {
            this->insert(*it);
        }
    }
    
    QSet<QString> functionNames() const;
};

}

SIRE_END_HEADER

#endif
