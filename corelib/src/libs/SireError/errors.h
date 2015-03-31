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

#ifndef SIREERROR_ERRORS_H
#define SIREERROR_ERRORS_H

#include "sireglobal.h"

#include <boost/assert.hpp>

//#ifndef BOOST_ENABLE_ASSERT_HANDLER
//#define BOOST_ENABLE_ASSERT_HANDLER 1
//#endif

#include "SireError/exception.h"
#include "SireError/version_error.h"

#include <stdexcept>
#include <QObject>

class QFile;
class QProcess;

SIRE_BEGIN_HEADER

namespace SireError
{

QString fileError(const QFile &f);

QString processError(const QString &executable, const QProcess &p);

/** This exception is thrown when the code detects a programming bug. I hope that
    eventually Sire will not throw too many of these!

    @author Christopher oods
*/
class SIREERROR_EXPORT program_bug : public SireError::exception
{
public:
    program_bug() : exception()
    {}

    program_bug(QString err, QString place = QString::null) : exception(err,place)
    {}

    program_bug(const program_bug &other) : exception(other)
    {}

    ~program_bug() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return program_bug::typeName();
    }
    
    void throwSelf() const
    {
        throw program_bug(*this);
    }
};

/** This exception is thrown when unsupported hardware is detected (e.g. unsupported
    openGL system)

    @author Christopher Woods
*/
class SIREERROR_EXPORT unsupported : public SireError::exception
{
public:
    unsupported() : exception()
    {}

    unsupported(QString err, QString place = QString::null) : exception(err,place)
    {}

    unsupported(const unsupported &other) : exception(other)
    {}

    ~unsupported() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return unsupported::typeName();
    }
    
    void throwSelf() const
    {
        throw unsupported(*this);
    }
};

/** This exception is thrown when an invalid key is used, e.g. a hash or set is
    used with a non-existant key.
    @author Christopher Woods
*/
class SIREERROR_EXPORT invalid_key : public SireError::exception
{
public:
    invalid_key() : exception()
    {}

    invalid_key(QString err, QString place = QString::null) : exception(err,place)
    {}

    invalid_key(const invalid_key &other) : exception(other)
    {}

    ~invalid_key() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return invalid_key::typeName();
    }
    
    void throwSelf() const
    {
        throw invalid_key(*this);
    }
};

/** This exception is thrown when an invalid index is used, e.g. an invalid index
    is used to access an array.
    @author Christopher Woods
*/
class SIREERROR_EXPORT invalid_index : public SireError::exception
{
public:
    invalid_index() : exception()
    {}

    invalid_index(QString err, QString place=QString::null) : exception(err,place)
    {}

    invalid_index(const invalid_index &other) : exception(other)
    {}

    ~invalid_index() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return invalid_index::typeName();
    }
    
    void throwSelf() const
    {
        throw invalid_index(*this);
    }
};

/** This exception is thrown when an invalid cast is attempted, e.g. a dynamic_cast
    that fails.
    @author Christopher Woods
*/
class SIREERROR_EXPORT invalid_cast : public SireError::exception
{
public:
    invalid_cast() : exception()
    {}

    invalid_cast(QString err, QString place = QString::null) : exception(err,place)
    {}

    invalid_cast(const invalid_cast &other) : exception(other)
    {}

    ~invalid_cast() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return invalid_cast::typeName();
    }
    
    void throwSelf() const
    {
        throw invalid_cast(*this);
    }
};

/** This exception is thrown whenever something incompatible is requested, e.g.
    trying to combine two incompatible objects.

    @author Christopher Woods
*/
class SIREERROR_EXPORT incompatible_error : public SireError::exception
{
public:
    incompatible_error() : exception()
    {}

    incompatible_error(QString err, QString place = QString::null) : exception(err,place)
    {}

    incompatible_error(const incompatible_error &other) : exception(other)
    {}

    ~incompatible_error() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return incompatible_error::typeName();
    }
    
    void throwSelf() const
    {
        throw incompatible_error(*this);
    }
};

/** This exception is thrown when request is made of an unknown type, e.g.
    when using QMetaType

    @author Christopher Woods
*/
class SIREERROR_EXPORT unknown_type : public SireError::exception
{
public:
    unknown_type() : exception()
    {}

    unknown_type(QString err, QString place = QString::null) : exception(err,place)
    {}

    unknown_type(const unknown_type &other) : exception(other)
    {}

    ~unknown_type() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return unknown_type::typeName();
    }
    
    void throwSelf() const
    {
        throw unknown_type(*this);
    }
};

/** This exception is thrown whenever an attempt is made to copy an uncopyable object.
    This will normally be a result of a programming bug!
    @author Christopher Woods
*/
class SIREERROR_EXPORT noncopyable_error : public SireError::exception
{
public:
    noncopyable_error() : exception()
    {}

    noncopyable_error(QString err, QString place=QString::null) : exception(err,place)
    {}

    noncopyable_error(const noncopyable_error &other) : exception(other)
    {}

    ~noncopyable_error() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return noncopyable_error::typeName();
    }
    
    void throwSelf() const
    {
        throw noncopyable_error(*this);
    }
};

/** This exception is thrown when an attempt is made to access a null pointer.
    @author Christopher Woods
*/
class SIREERROR_EXPORT nullptr_error : public SireError::exception
{
public:
    nullptr_error() : exception()
    {}

    nullptr_error(QString err, QString place=QString::null) : exception(err,place)
    {}

    nullptr_error(const nullptr_error &other) : exception(other)
    {}

    ~nullptr_error() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return nullptr_error::typeName();
    }
    
    void throwSelf() const
    {
        throw nullptr_error(*this);
    }
};

/** This exception is thrown whenever there is an error with locking a resource,
    e.g. a file lock, or mutex lock.
    @author Christopher Woods
*/
class SIREERROR_EXPORT lock_error : public SireError::exception
{
public:
    lock_error() : exception()
    {}

    lock_error(QString err, QString place=QString::null) : exception(err,place)
    {}

    lock_error(const lock_error &other) : exception(other)
    {}

    ~lock_error() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return lock_error::typeName();
    }
    
    void throwSelf() const
    {
        throw lock_error(*this);
    }
};

/** This exception is thrown whenever an assertation fails. This is normally
    a result of a coding bug where a function is being used in an unintended
    and unsupported way.
    @author Christopher Woods
*/
class SIREERROR_EXPORT assertation_failed : public SireError::exception
{
public:
    assertation_failed() : exception()
    {}

    assertation_failed(QString err, QString place=QString::null)
          : exception(err,place)
    {}

    assertation_failed(const invalid_key &other) : exception(other)
    {}

    ~assertation_failed() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return assertation_failed::typeName();
    }
    
    void throwSelf() const
    {
        throw assertation_failed(*this);
    }
};

/** This exception is thrown whenever there is an file handling error */
class SIREERROR_EXPORT file_error : public SireError::exception
{
public:
    file_error() : exception()
    {}

    file_error(const QString &err, const QString &place=QString::null)
              : exception(err, place)
    {}

    file_error(const QFile &file, QString place=QString::null)
              : exception(fileError(file), place)
    {}

    file_error(const file_error &other) : exception(other)
    {}

    ~file_error() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return file_error::typeName();
    }
    
    void throwSelf() const
    {
        throw file_error(*this);
    }
};

/** This exception is thrown whenever there is an error with a child process */
class SIREERROR_EXPORT process_error : public SireError::exception
{
public:
    process_error() : exception()
    {}

    process_error(const QString &err, const QString place=QString::null)
              : exception(err, place)
    {}

    process_error(const QString &executable, const QProcess &process,
                  QString place=QString::null)
              : exception(processError(executable,process), place)
    {}

    process_error(const process_error &other) : exception(other)
    {}

    ~process_error() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return process_error::typeName();
    }
    
    void throwSelf() const
    {
        throw process_error(*this);
    }
};

/** This exception is thrown whenever a request is made of a resource
    that is just not available */
class SIREERROR_EXPORT unavailable_resource : public SireError::exception
{
public:
    unavailable_resource() : exception()
    {}

    unavailable_resource(QString err, QString place=QString::null)
              : exception(err,place)
    {}

    unavailable_resource(const unavailable_resource &other) : exception(other)
    {}

    ~unavailable_resource() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return unavailable_resource::typeName();
    }
    
    void throwSelf() const
    {
        throw unavailable_resource(*this);
    }
};

/** This exception is thrown whenever there is an input/output error */
class SIREERROR_EXPORT io_error : public SireError::exception
{
public:
    io_error() : exception()
    {}

    io_error(QString err, QString place=QString::null)
              : exception(err,place)
    {}

    io_error(const io_error &other) : exception(other)
    {}

    ~io_error() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return io_error::typeName();
    }
    
    void throwSelf() const
    {
        throw io_error(*this);
    }
};

/** This exception is thrown whenever a function is called with invalid arguments.
    This is normally as a result of a programming bug, or an attempt to use a function
    in an unsupported way.
    @author Christopher Woods
*/
class SIREERROR_EXPORT invalid_arg : public SireError::exception
{
public:
    invalid_arg() : exception()
    {}

    invalid_arg(QString err, QString place=QString::null) : exception(err,place)
    {}

    invalid_arg(const invalid_arg &other) : exception(other)
    {}

    ~invalid_arg() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return invalid_arg::typeName();
    }
    
    void throwSelf() const
    {
        throw invalid_arg(*this);
    }
};

/** This exception is thrown whenever the code detects that it has entered an
    invalid state. This will be as a result of a programming bug, or undocumented
    use of the code.
    @author Christopher Woods
*/
class SIREERROR_EXPORT invalid_state : public SireError::exception
{
public:
    invalid_state() : exception()
    {}

    invalid_state(QString err, QString place=QString::null) : exception(err,place)
    {}

    invalid_state(const invalid_state &other) : exception(other)
    {}

    ~invalid_state() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return invalid_state::typeName();
    }
    
    void throwSelf() const
    {
        throw invalid_state(*this);
    }
};

/** This exception is thrown whenever an invalid operation is requested

    @author Christopher Woods
*/
class SIREERROR_EXPORT invalid_operation : public SireError::exception
{
public:
    invalid_operation() : exception()
    {}

    invalid_operation(QString err, QString place=QString::null) : exception(err,place)
    {}

    invalid_operation(const invalid_operation &other) : exception(other)
    {}

    ~invalid_operation() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return invalid_operation::typeName();
    }
    
    void throwSelf() const
    {
        throw invalid_operation(*this);
    }
};

/** This exception is thrown whenever the program tries to use code that has yet
    to be written. This is definitely a programming bug, and results from the use
    of stubs that I have forgotton to expand!
    @author Christopher Woods
*/
class SIREERROR_EXPORT incomplete_code : public SireError::exception
{
public:
    incomplete_code() : exception()
    {}

    incomplete_code(QString err, QString place=QString::null) : exception(err,place)
    {}

    incomplete_code(const incomplete_code &other) : exception(other)
    {}

    ~incomplete_code() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return incomplete_code::typeName();
    }
    
    void throwSelf() const
    {
        throw incomplete_code(*this);
    }
};

/** This exception is used to translate a std::exception into a SireError::exception */
class SIREERROR_EXPORT std_exception : public SireError::exception
{
public:
    std_exception();
    std_exception(const std::logic_error &error);
    std_exception(const std::domain_error &error);
    std_exception(const std::invalid_argument &error);
    std_exception(const std::length_error &error);
    std_exception(const std::out_of_range &error);
    std_exception(const std::runtime_error &error);
    std_exception(const std::range_error &error);
    std_exception(const std::overflow_error &error);
    std_exception(const std::underflow_error &error);
    std_exception(const std::bad_alloc &error);
    std_exception(const std::exception &error);

    ~std_exception() throw();

    static const char* typeName();

    const char* what() const throw()
    {
        return std_exception::typeName();
    }
    
    void throwSelf() const
    {
        throw std_exception(*this);
    }

protected:
    static QString getWhatString(QString typstring, const std::exception &error);
};

/** This exception is thrown whenever an unidentified exception needs translating into
    a SireError::exception, or when an unidentified SireError::exception is detected.
    
    @author Christopher Woods
*/
class SIREERROR_EXPORT unknown_exception : public SireError::exception
{
public:
    unknown_exception() : exception()
    {}

    unknown_exception(QString err, QString place=QString::null) : exception(err,place)
    {}

    unknown_exception(const unknown_exception &other) : exception(other)
    {}

    ~unknown_exception() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return unknown_exception::typeName();
    }
    
    void throwSelf() const
    {
        throw unknown_exception(*this);
    }
};

/** This exception is thrown whenever a problem is detected
    that involves dependencies, e.g. a dependency is missing,
    or a dependency may be about to be broken.
    
    @author Christopher Woods
*/
class SIREERROR_EXPORT dependency_error : public SireError::exception
{
public:
    dependency_error() : exception()
    {}

    dependency_error(QString err, QString place=QString::null) : exception(err,place)
    {}

    dependency_error(const dependency_error &other) : exception(other)
    {}

    ~dependency_error() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return dependency_error::typeName();
    }
    
    void throwSelf() const
    {
        throw dependency_error(*this);
    }
};

/** This exception is thrown whenever there is an error with an
    ID number

    @author Christopher Woods
*/
class SIREERROR_EXPORT id_error : public exception
{
public:
    id_error() : exception()
    {}

    id_error(QString err, QString place = QString::null)
                  : exception(err,place)
    {}

    id_error(const id_error &other) : exception(other)
    {}

    ~id_error() throw()
    {}

    static const char* typeName();

    const char* what() const throw()
    {
        return id_error::typeName();
    }
    
    void throwSelf() const
    {
        throw id_error(*this);
    }
};

}

Q_DECLARE_METATYPE(SireError::program_bug)
Q_DECLARE_METATYPE(SireError::unsupported)
Q_DECLARE_METATYPE(SireError::id_error)
Q_DECLARE_METATYPE(SireError::invalid_key)
Q_DECLARE_METATYPE(SireError::invalid_index)
Q_DECLARE_METATYPE(SireError::invalid_cast)
Q_DECLARE_METATYPE(SireError::noncopyable_error)
Q_DECLARE_METATYPE(SireError::nullptr_error)
Q_DECLARE_METATYPE(SireError::lock_error)
Q_DECLARE_METATYPE(SireError::assertation_failed)
Q_DECLARE_METATYPE(SireError::incompatible_error)
Q_DECLARE_METATYPE(SireError::file_error)
Q_DECLARE_METATYPE(SireError::io_error)
Q_DECLARE_METATYPE(SireError::process_error)
Q_DECLARE_METATYPE(SireError::invalid_arg)
Q_DECLARE_METATYPE(SireError::invalid_state)
Q_DECLARE_METATYPE(SireError::invalid_operation)
Q_DECLARE_METATYPE(SireError::unavailable_resource)
Q_DECLARE_METATYPE(SireError::incomplete_code)
Q_DECLARE_METATYPE(SireError::std_exception)
Q_DECLARE_METATYPE(SireError::unknown_exception)
Q_DECLARE_METATYPE(SireError::unknown_type)
Q_DECLARE_METATYPE(SireError::dependency_error)

SIRE_END_HEADER

#endif
