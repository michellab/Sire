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

#ifndef SIREERROR_EXCEPTION_H
#define SIREERROR_EXCEPTION_H

#include <exception>

#include <QString>
#include <QStringList>

#include <boost/shared_ptr.hpp>

#include "sireglobal.h"

SIRE_BEGIN_HEADER

namespace SireError
{
class exception;
}

class QDataStream;
SIRESTREAM_EXPORT QDataStream& operator<<(QDataStream&, const SireError::exception&);
SIRESTREAM_EXPORT QDataStream& operator>>(QDataStream&, SireError::exception&);

namespace SireError
{

/** This is a small class that is used to ensure
    that fast exceptions are switched off when no longer
    needed */
class SIREERROR_EXPORT FastExceptionFlag
{
public:
    FastExceptionFlag();
    FastExceptionFlag(const FastExceptionFlag &other);
    ~FastExceptionFlag();

    FastExceptionFlag& operator=(const FastExceptionFlag &other);

    void disable();

private:
    friend class exception;
    static FastExceptionFlag construct();

    class FastExceptionFlagData
    {
    public:
        FastExceptionFlagData();
        ~FastExceptionFlagData();
    };

    boost::shared_ptr<FastExceptionFlagData> d;

    /** Whether or not fast exceptions are enabled */
    static bool enable_fast_exceptions;
};

/** This is the base class of all Sire specific exceptions. The python wrapping
    allows for automatic conversion of any exception derived from this base.
    The exception system in Sire is very basic, namely there are a collection
    of named exceptions. Each exception only has basic information such as
    from where is was thrown, and an optional message to say why is was thrown.

    @author Christopher Woods
*/
class SIREERROR_EXPORT exception : public std::exception
{

friend SIRESTREAM_EXPORT QDataStream& ::operator<<(QDataStream&, const exception&);
friend SIRESTREAM_EXPORT QDataStream& ::operator>>(QDataStream&, exception&);

public:
    typedef SireError::exception ROOT;

    exception();
    exception(QString error, QString place = QString());

    exception(const exception &other);

    virtual ~exception() throw();

    static const char* typeName()
    {
        return "SireError::exception";
    }

    virtual const char* what() const throw()=0;

    exception* clone() const;

    QByteArray pack() const;
    static boost::shared_ptr<SireError::exception> unpack(const QByteArray &data);

    static void unpackAndThrow(const QByteArray &errordata);

    static FastExceptionFlag enableFastExceptions();

    QString error() const throw();
    QString from() const throw();
    QStringList trace() const throw();
    QString where() const throw();
    QString why() const throw();
    QString pid() const throw();

    QString toString() const throw();

    virtual void throwSelf() const=0;

protected:
    QString err;  ///< The error associated with the exception.
    QString plce; ///< Description of from where the exception was thrown.
    QStringList bt; ///< Backtrace at the point that the exception was constructed
    QString pidstr; /**< String identifying the process from which
                             process/thread this exception was thrown */
};

}

SIRE_END_HEADER

#endif
