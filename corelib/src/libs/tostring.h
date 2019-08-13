/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2007  Christopher Woods
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

#ifndef SIRE_TOSTRING_H
#define SIRE_TOSTRING_H

#include <QString>
#include <QStringList>

#include <QVector>
#include <QList>
#include <QSet>
#include <QHash>
#include <QMultiHash>
#include <QMap>

#ifndef GCCXML_PARSE
   // GCCXML chokes on these files
   #include <boost/type_traits/is_pod.hpp>
   #include <boost/mpl/if.hpp>
#endif

SIRE_BEGIN_HEADER

namespace SireBase
{
template<class T>
class PropPtr;
}

namespace Sire
{

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** @file tostring.h

    This file contains the functions necessary to implement a generic "toString"
    function that obtains a string representation of most objects (including
    containers)
    
    @author Christopher Woods
*/

/** Return a string representation of a POD type */
template<typename T>
struct qstr_numtype
{
    static QString qstr(const T &val)
    {
        return QString::number(val);
    }
};

/** Return a string representation of a class type (must have 
    implemented a .toString() member function, or else you
    must supply a specific instantiation) */
template<class T>
struct qstr_classtype
{
    static QString qstr(const T &obj)
    {
        return obj.toString();
    }
};

/** Used to return the string of a non-container type */
template<class T>
QString qstr(const T &obj)
{
    return boost::mpl::if_<boost::is_pod<T>,
                     qstr_numtype<T>, qstr_classtype<T> >
                         ::type::qstr(obj);
}

/** Used to return a string representation of a QString! */
template<>
SIRE_ALWAYS_INLINE QString qstr(const QString &string)
{
    return string;
}

/** Used to return a string representation of a QStringList */
template<>
SIRE_ALWAYS_INLINE QString qstr(const QStringList &strings)
{
    if (strings.isEmpty())
        return QString("[ ]");
    else if (strings.count() == 1)
        return strings.at(0);
    else
        return QString("[ %1 ]").arg( strings.join(", ") );
}

/** Used to return a string representation of a SireBase::PropPtr */
template<class T>
QString qstr(const SireBase::PropPtr<T> &ptr)
{
    return qstr( ptr.read() );
}

/** Used to return a string representation of a QList */
template<class T>
QString qstr(const QList<T> &objs)
{
    QStringList strngs;
    
    for (typename QList<T>::const_iterator it = objs.begin();
         it != objs.end();
         ++it)
    {
        strngs.append( qstr(*it) );
    }
    
    return QString("[ %1 ]").arg( strngs.join(","));
}

/** Used to return a string representation of a QVector */
template<class T>
QString qstr(const QVector<T> &objs)
{
    QStringList strngs;
    
    for (typename QVector<T>::const_iterator it = objs.begin();
         it != objs.end();
         ++it)
    {
        strngs.append( qstr(*it) );
    }
    
    return QString("[ %1 ]").arg( strngs.join(","));
}

/** Used to return a string representation of a QSet */
template<class T>
QString qstr(const QSet<T> &objs)
{
    QStringList strngs;
    
    for (typename QSet<T>::const_iterator it = objs.begin();
         it != objs.end();
         ++it)
    {
        strngs.append( qstr(*it) );
    }
    
    return QString("{ %1 }").arg( strngs.join(","));
}

/** Used to return a string representation of a QHash */
template<class S, class T>
QString qstr(const QHash<S,T> &objs)
{
    QStringList strngs;
    
    for (typename QHash<S,T>::const_iterator it = objs.begin();
         it != objs.end();
         ++it)
    {
        strngs.append( QString(" %1 : %2").arg( qstr(it.key()), qstr(it.value()) ) );
    }
    
    return QString("{ %1 }").arg( strngs.join(","));
}

/** Used to return a string representation of a QMultiHash */
template<class S, class T>
QString qstr(const QMultiHash<S,T> &objs)
{
    QStringList strngs;
    
    for (typename QMultiHash<S,T>::const_iterator it = objs.begin();
         it != objs.end();
         ++it)
    {
        strngs.append( QString(" %1 : %2").arg( qstr(it.key()), qstr(it.value()) ) );
    }
    
    return QString("{ %1 }").arg( strngs.join(","));
}

/** Used to return a string representation of a QMap */
template<class S, class T>
QString qstr(const QMap<S,T> &objs)
{
    QStringList strngs;
    
    for (typename QMap<S,T>::const_iterator it = objs.begin();
         it != objs.end();
         ++it)
    {
        strngs.append( QString(" %1 : %2").arg( qstr(it.key()), qstr(it.value()) ) );
    }
    
    return QString("{ %1 }").arg( strngs.join(","));
}

/** Return a QString representation of the object 'obj' */
template<class T>
QString toString(const T &obj)
{
    return qstr(obj);
}

#else // SIRE_SKIP_INLINE_FUNCTIONS

template<class T>
QString qstr(const T &obj);

template<class T>
QString toString(const T &obj);

#endif // SIRE_SKIP_INLINE_FUNCTIONS

}

SIRE_END_HEADER

#endif
