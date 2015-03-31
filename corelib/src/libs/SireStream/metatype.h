/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2013  Christopher Woods
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

#ifndef SIRESTREAM_METATYPE_H
#define SIRESTREAM_METATYPE_H

#include <QSet>
#include <QString>

#include "sireglobal.h"

SIRE_BEGIN_HEADER

namespace SireStream
{
    bool isLeafOf(const char *base, const QString &type_name);
    bool isBranchOf(const char *base, const QString &type_name);
    bool isDerivedFrom(const char *base, const QString &type_name);
    bool isRootless(const QString &type_name);
    bool isRegistered(const QString &type_name);
    
    QString registeredRoot(const QString &type_name);
    
    QSet<QString> leafClassesOf(const char *base);
    QSet<QString> branchClassesOf(const char *base);
    QSet<QString> derivedClassesOf(const char *base);

    QSet<QString> rootlessClasses();
    QSet<QString> registeredClasses();

    /** Return whether the class called "type_name" is derived from the root (base) class T */
    template<class T>
    bool isDerivedFrom(const QString &type_name)
    {
        return isDerivedFrom(T::typeName(), type_name);
    }
    
    /** Return whether the class "T" is derived from the root (base) class "S" */
    template<class S, class T>
    bool isDerivedFrom()
    {
        return isDerivedFrom(S::typeName(), QLatin1String(T::typeName()));
    }

    /** Return whether the class called "type_name" is a leaf of the root (base) class T */
    template<class T>
    bool isLeafOf(const QString &type_name)
    {
        return isLeafOf(T::typeName(), type_name);
    }
    
    /** Return whether the class "T" is a leaf of the root (base) class "S" */
    template<class S, class T>
    bool isLeafOf()
    {
        return isLeafOf(S::typeName(), QLatin1String(T::typeName()));
    }

    /** Return whether the class called "type_name" is a branch of the root (base) class T */
    template<class T>
    bool isBranchOf(const QString &type_name)
    {
        return isBranchOf(T::typeName(), type_name);
    }
    
    /** Return whether the class "T" is a branch of the root (base) class "S" */
    template<class S, class T>
    bool isBranchOf()
    {
        return isBranchOf(S::typeName(), QLatin1String(T::typeName()));
    }
    
    /** Return whether or not class T is rootless */
    template<class T>
    bool isRootless()
    {
        return isRootless(T::typeName());
    }
    
    /** Return the registered root of class T. This returns an empty
        string if the class has not been registered */
    template<class T>
    QString registeredRoot()
    {
        return registeredRoot( T::typeName() );
    }
    
    /** Return whether or not class T is registered */
    template<class T>
    bool isRegistered()
    {
        return isRegistered(T::typeName());
    }
    
    /** Return all of the leaf classes derived from root/base class T */
    template<class T>
    QSet<QString> leafClassesOf()
    {
        return leafClassesOf(T::typeName());
    }
    
    /** Return all of the branch classes derived from root/base class T */
    template<class T>
    QSet<QString> branchClassesOf()
    {
        return branchClassesOf(T::typeName());
    }
    
    /** Return all of the classes derived from root/base class T */
    template<class T>
    QSet<QString> derivedClassesOf()
    {
        return derivedClassesOf(T::typeName());
    }
}

SIRE_END_HEADER

#endif
