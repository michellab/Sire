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

#include "SireError/exception.h"
#include "SireStream/datastream.h"

#include "SireError/printerror.h"
#include "SireError/errors.h"

#include "getbacktrace.h"

#include <QDataStream>
#include <QThreadStorage>

#include <QDebug>

#ifdef _HAVE_BOOST_STACKTRACE_HPP_
#define BOOST_STACKTRACE_GNU_SOURCE_NOT_REQUIRED 1
#include <boost/stacktrace.hpp>
#endif

using namespace SireError;
using namespace SireStream;

Q_GLOBAL_STATIC( QThreadStorage<QString*>, pidStrings );
Q_GLOBAL_STATIC( QString, processString );

namespace Sire
{
namespace detail
{
    QHash< QString, QSet<QString> > branch_classes;
    QHash< QString, QSet<QString> > leaf_classes;
    QSet<QString> rootless_classes;

    const QHash< QString, QSet<QString> > branchClasses()
    {
        return branch_classes;
    }

    const QHash< QString, QSet<QString> > leafClasses()
    {
        return leaf_classes;
    }

    const QSet<QString> rootlessClasses()
    {
        return rootless_classes;
    }

    void registerLeaf(const QString &type_name, const char *root)
    {
        QLatin1String r(root);
        if (not leaf_classes.contains(r))
        {
            leaf_classes.insert(r, QSet<QString>());
        }
        
        leaf_classes[r].insert(type_name);
    }
    
    void registerBranch(const QString &type_name, const char *root)
    {
        QLatin1String r(root);
        if (not branch_classes.contains(r))
        {
            branch_classes.insert(r, QSet<QString>());
        }
        
        branch_classes[r].insert(type_name);
    }
    
    void registerRootless(const QString &type_name)
    {
        rootless_classes.insert(type_name);
    }

} // end of namespace detail
} // end of namespace Sire

namespace SireError
{

/** Set the string that SireError will use to identify this process */
void setProcessString(const QString &s)
{
    *(processString()) = s;
}

/** Set the string that SireError will used to identify this thread
    within the process */
void setThreadString(const QString &s)
{
    pidStrings()->setLocalData( new QString(s) );
}

/** Return the string used by SireError to identify the process */
QString getProcessString()
{
    QString *s = processString();
    
    if (s->isEmpty())
    {
        return *s = QObject::tr("master");
    }
    
    return *s;
}

/** Return the string used to identify a particular thread */
QString getThreadString()
{
    QThreadStorage<QString*> *store = pidStrings();
    
    if (store->hasLocalData())
    {
        return *(store->localData());
    }
    else
    {
        return QString::null;
    }
}

/** Return the string used by SireError to identify a particular
    thread within a process */
QString getPIDString()
{
    QThreadStorage<QString*> *store = pidStrings();
    
    if (store->hasLocalData())
    {
        return QString("%1:%2").arg( getProcessString(),
                                     *(store->localData()) );
    }
    else
    {
        return getProcessString();
    }
}

} // end of namespace SireError

bool FastExceptionFlag::enable_fast_exceptions = false;

FastExceptionFlag::FastExceptionFlag()
{}

FastExceptionFlag::FastExceptionFlag(const FastExceptionFlag &other) : d(other.d)
{}

FastExceptionFlag::~FastExceptionFlag()
{}

FastExceptionFlag& FastExceptionFlag::operator=(const FastExceptionFlag &other)
{
    d = other.d;
    return *this;
}

FastExceptionFlag FastExceptionFlag::construct()
{
    FastExceptionFlag f;
    f.d.reset( new FastExceptionFlagData() );
    return f;
}

void FastExceptionFlag::disable()
{
    enable_fast_exceptions = false;
    d.reset();
}

FastExceptionFlag::FastExceptionFlagData::FastExceptionFlagData()
{
    enable_fast_exceptions = true;
}

FastExceptionFlag::FastExceptionFlagData::~FastExceptionFlagData()
{
    enable_fast_exceptions = false;
}

/** Switch on fast exceptions. These are used, e.g. when you know that
    you are going to be handling exceptions yourself, are throwing thousands
    of them, and thus don't want to slow the code down generating
    extra backtrace or other information */
FastExceptionFlag exception::enableFastExceptions()
{
    return FastExceptionFlag::construct();
}

/** Construct a null exception */
exception::exception()
{
    if (FastExceptionFlag::enable_fast_exceptions)
        return;

    //pidstr = getPIDString();
}

/** Constructor.
    \param error The error message associated with this exception.
    \param place From where in the code this exception was thrown. This
                 is supplied automatically via the 'CODELOC' macro.
*/
exception::exception(QString error, QString place) : err(error), plce(place)
{
    if (FastExceptionFlag::enable_fast_exceptions)
        return;

    #if defined(SIRE_ENABLE_BACKTRACE) || defined(SIRE_ENABLE_BOOST_BACKTRACE)
    #ifdef SIRE_ENABLE_BACKTRACE
        bt = getBackTrace();
    #else
        std::stringstream ss;
        ss << boost::stacktrace::stacktrace();
        bt = QString::fromStdString(ss.str()).split("\n");
    #endif
        pidstr = getPIDString();
    #endif
}

/** Copy constructor */
exception::exception(const exception &other)
          : std::exception(other), err(other.err), plce(other.plce), 
                                   bt(other.bt), pidstr(other.pidstr)
{}

/** Destructor */
exception::~exception() throw()
{}

/** Return a clone of this exception */
exception* exception::clone() const
{
    //get the ID number of this type
    int id = QMetaType::type( this->what() );

    if ( id == 0 or not QMetaType::isRegistered(id) )
        throw SireError::unknown_type(QObject::tr(
            "The exception with type \"%1\" does not appear to have been "
            "registered with QMetaType. It cannot be cloned! (%2, %3)")
                .arg(this->what()).arg(id).arg(QMetaType::isRegistered(id)), CODELOC);

    return static_cast<exception*>( QMetaType::create(id,this) );
}

/** Pack this exception into a binary QByteArray - this packs the exception
    with its type information */
QByteArray exception::pack() const
{
    QByteArray data;
    
    //reserve 128K of space for the exception (should be way
    //more than enough!)
    data.reserve( 128 * 1024 );
    
    QDataStream ds(&data, QIODevice::WriteOnly);
    
    //get the ID number of this type
    int id = QMetaType::type( this->what() );

    if ( id == 0 or not QMetaType::isRegistered(id) )
        throw SireError::unknown_type(QObject::tr(
            "The exception with type \"%1\" does not appear to have been "
            "registered with QMetaType. It cannot be streamed! (%2, %3)")
                .arg(this->what()).arg(id).arg(QMetaType::isRegistered(id)), CODELOC);

    //save the object type name
    ds << QString(this->what());

    //use the QMetaType streaming function to save this table
    if (not QMetaType::save(ds, id, this))
        throw SireError::program_bug(QObject::tr(
            "There was an error saving the exception of type \"%1\". "
            "Has the programmer added a RegisterMetaType for this exception?")
                .arg(this->what()), CODELOC);

    return data;
}

/** Unpack an exception from the raw data in the passed bytearray */
boost::shared_ptr<SireError::exception> exception::unpack(const QByteArray &data)
{
    QDataStream ds(data);
    
    //read the type of exception first
    QString type_name;
    
    ds >> type_name;
    
    //get the type that represents this name
    int id = QMetaType::type( type_name.toLatin1().constData() );

    if ( id == 0 or not QMetaType::isRegistered(id) )
        throw SireError::unknown_type( QObject::tr(
              "Cannot deserialise an exception of type \"%1\". "
              "Ensure that the library or module containing "
              "this exception has been loaded and that it has been registered "
              "with QMetaType.").arg(type_name), CODELOC );
    
    //construct an exception of this type
    boost::shared_ptr<exception> ptr(  
                        static_cast<exception*>(QMetaType::create(id,0)) );
                        
    //load the object from the datastream
    if ( not QMetaType::load(ds, id, ptr.get()) )
        throw SireError::program_bug(QObject::tr(
            "There was an error loading the exception of type \"%1\"")
                 .arg(type_name), CODELOC);

    return ptr;
}

/** Unpack the exception contained in the raw binary data 'data', and then
    throw that exception */
void exception::unpackAndThrow(const QByteArray &data)
{
    boost::shared_ptr<exception> e = exception::unpack(data);
    e->throwSelf();
}

/** Return a pretty string representation of all of the information stored
    in the exception, suitable for printing to the screen or to a log file. */
QString exception::toString() const throw()
{
    return QObject::tr("Exception '%1' thrown by the thread '%2'.\n"
                       "%3\n"
                       "Thrown from %4\n"
                       "__Backtrace__\n"
                       "%5\n"
                       "__EndTrace__\n"
                       "Exception '%1' thrown by the thread '%2'.\n"
                       "%3\n"
                       "Thrown from %4")
             .arg(what()).arg(pid()).arg(why()).arg(where()).arg(trace().join("\n"));
}

/** Return the error message associated with this exception.
    \return The error message associated with this exception.
*/
QString exception::error() const throw()
{
    return err;
}

/** Return the location (function, line, file) from which this exception
    was originally thrown. The amount of information available will depend
    on the compiler used to compile the program, and come from the 'CODELOC'
    macro.
    \return Description of from where the exception was thrown.
*/
QString exception::from() const throw()
{
    return plce;
}

/** Return the function backtrace when the exception was constructed */
QStringList exception::trace() const throw()
{
    #if defined(SIRE_ENABLE_BACKTRACE) || defined(SIRE_ENABLE_BOOST_BACKTRACE)
        return bt;
    #else
        QStringList btrace;
        btrace.append( QObject::tr("No backtrace available. "
                "Recompile the SireError library with -DSIRE_ENABLE_BACKTRACE "
                "to enable backtraces.") );
    
        return btrace;
    #endif
}

/** Overloaded functions to give logical names... */
QString exception::where() const throw()
{
    return this->from();
}

/** Overloaded functions to give logical names... */
QString exception::why() const throw()
{
    return this->error();
}

/** Return the identifying string for the thread on the process that threw
    this exception. You should call the static SireError::exception::setPID
    function in each thread for the result of this call to be intelligable */
QString exception::pid() const throw()
{
    return pidstr;
}
