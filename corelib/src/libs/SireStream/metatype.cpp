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

#include "metatype.h"

namespace Sire
{
    namespace detail
    {
        const QHash< QString, QSet<QString> > branchClasses();
        const QHash< QString, QSet<QString> > leafClasses();
        const QSet<QString> rootlessClasses();
    }
}

namespace SireStream
{
    /** Return whether the class 'type_name' is registered as a leaf of
        the base 'base' */
    bool SIRESTREAM_EXPORT isLeafOf(const char *base, const QString &type_name)
    {
        QHash< QString, QSet<QString> >::const_iterator
                    it = Sire::detail::leafClasses().constFind( QLatin1String(base) );
        
        if (it != Sire::detail::leafClasses().constEnd())
        {
            return it.value().contains(type_name);
        }
        else
            return false;
    }

    /** Return whether the class 'type_name' is registered as a branch of
        the base 'base' */
    bool SIRESTREAM_EXPORT isBranchOf(const char *base, const QString &type_name)
    {
        QHash< QString, QSet<QString> >::const_iterator
                    it = Sire::detail::branchClasses().constFind( QLatin1String(base) );
        
        if (it != Sire::detail::branchClasses().constEnd())
        {
            return it.value().contains(type_name);
        }
        else
            return false;

    }

    /** Return whether the class 'type_name' is registered to be derived
        from the base 'base' */
    bool SIRESTREAM_EXPORT isDerivedFrom(const char *base, const QString &type_name)
    {
        return isLeafOf(base, type_name) or isBranchOf(base, type_name);
    }
    
    /** Return whether or not the class 'type_name' has been registered
        as rootless */
    bool SIRESTREAM_EXPORT isRootless(const QString &type_name)
    {
        return Sire::detail::rootlessClasses().contains(type_name);
    }
    
    /** Return whether or not the class with name 'type_name' is registered
        with the metatype system */
    bool SIRESTREAM_EXPORT isRegistered(const QString &type_name)
    {
        if ( Sire::detail::rootlessClasses().contains(type_name) )
            return true;
        
        else
        {
            for (QHash< QString,QSet<QString> >::const_iterator
                        it = Sire::detail::leafClasses().constBegin();
                 it != Sire::detail::leafClasses().constEnd();
                 ++it)
            {
                if (it.value().contains(type_name))
                    return true;
            }
            
            for (QHash< QString,QSet<QString> >::const_iterator
                        it = Sire::detail::branchClasses().constBegin();
                 it != Sire::detail::branchClasses().constEnd();
                 ++it)
            {
                if (it.value().contains(type_name))
                    return true;
            }
            
            return false;
        }
    }
    
    /** Return the name of the class registered as the root of class 'type_name'.
        This returns an empty string if there is no registered root,
        or this class has not been registered */
    QString SIRESTREAM_EXPORT registeredRoot(const QString &type_name)
    {
        if (Sire::detail::rootlessClasses().contains(type_name))
        {
            return QString::null;
        }
        else
        {
            for (QHash< QString,QSet<QString> >::const_iterator
                        it = Sire::detail::leafClasses().constBegin();
                 it != Sire::detail::leafClasses().constEnd();
                 ++it)
            {
                if (it.value().contains(type_name))
                    return it.key();
            }

            for (QHash< QString,QSet<QString> >::const_iterator
                        it = Sire::detail::branchClasses().constBegin();
                 it != Sire::detail::branchClasses().constEnd();
                 ++it)
            {
                if (it.value().contains(type_name))
                    return it.key();
            }
            
            return QString::null;
        }
    }
    
    /** Return all of the leaf classes of root class 'base' */
    QSet<QString> SIRESTREAM_EXPORT leafClassesOf(const char *base)
    {
        return Sire::detail::leafClasses().value( QLatin1String(base) );
    }
    
    /** Return all of the branch classes of root class 'base' */
    QSet<QString> SIRESTREAM_EXPORT branchClassesOf(const char *base)
    {
        return Sire::detail::branchClasses().value( QLatin1String(base) );
    }
    
    /** Return all of the classes registered as derived from root class 'base' */
    QSet<QString> SIRESTREAM_EXPORT derivedClassesOf(const char *base)
    {
        return leafClassesOf(base).unite(branchClassesOf(base));
    }
    
    /** Return the set of all rootless registered classes */
    QSet<QString> SIRESTREAM_EXPORT rootlessClasses()
    {
        return Sire::detail::rootlessClasses();
    }
    
    /** Return the set of all registered classes */
    QSet<QString> SIRESTREAM_EXPORT registeredClasses()
    {
        QSet<QString> classes = rootlessClasses();
        
        for (QHash< QString,QSet<QString> >::const_iterator
                    it = Sire::detail::leafClasses().constBegin();
             it != Sire::detail::leafClasses().constEnd();
             ++it)
        {
            classes.unite(it.value());
        }

        for (QHash< QString,QSet<QString> >::const_iterator
                    it = Sire::detail::branchClasses().constBegin();
             it != Sire::detail::branchClasses().constEnd();
             ++it)
        {
            classes.unite(it.value());
        }
        
        return classes;
    }
}
