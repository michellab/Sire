/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2016  Christopher Woods
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

#ifndef SIREBASE_UNITTEST_H
#define SIREBASE_UNITTEST_H

#include "SireError/errors.h"
#include "tostring.h"

#ifdef SIRE_HAS_CPP_11
    #include <functional>
#endif

#include <QList>
#include <QMutex>
#include <QString>

#include <boost/shared_ptr.hpp>

namespace SireBase
{

/** This class is used to register a unit test 

    @author Christopher Woods
*/
class SIREBASE_EXPORT UnitTest
{
public:
    UnitTest();
    UnitTest(const QString &name, void (*test_func)(bool));
    UnitTest(const UnitTest &other);
    
    ~UnitTest();
    
    bool run(int nrepeats=1, bool verbose=false);

    bool wasSuccessful();
    bool wasError();
    
    QString name() const;
    
    QString errorString();
    
    quint64 runTime();

    static QList< boost::shared_ptr<UnitTest> > tests();

    static int runAll(bool verbose=false);

private:
    /** Mutex used to protect access to this test */
    QMutex run_mutex;

    /** The name of the test */
    QString test_name;

    /** Pointer to the function to run */
    void (*test_function)(bool);
    
    /** The error string. This is empty if the test has
        not been run or was successful */
    QString error_string;
    
    /** How long the test took to run (on average across all runs) */
    quint64 run_time;
    
    /** The set of all tests that have been registered */
    static QList< boost::shared_ptr<UnitTest> > all_tests;
};

#define SIRE_UNITTEST(X) static SireBase::UnitTest _SireBase_UnitTest_register_##X( #X, &X );

/** Assert that the passed two objects are equal */
template<class T>
void assert_equal( const T &a, const T &b, const QString &code_location )
{
    if (a != b)
    {
        throw SireError::assertation_failed( QObject::tr("NOT EQUAL!\n\n%1\n\n!=\n\n%2")
			.arg(Sire::toString(a)).arg(Sire::toString(b)), code_location );
    }
}

/** Assert that the passed two objects are not equal */
template<class T>
void assert_not_equal( const T &a, const T &b, const QString &code_location )
{
    if (a == b)
    {
        throw SireError::assertation_failed( QObject::tr("EQUAL!\n\n%1\n\n!=\n\n%2")
			.arg(Sire::toString(a)).arg(Sire::toString(b)), code_location );
    }
}

/** Assert that hte passed two objects are equal to within the specified range */
template<class T>
void assert_nearly_equal( const T &a, const T &b, const T &range, const QString &code_location )
{
    if ((a - b < -range) or (a - b > range))
    {                                    
        throw SireError::assertation_failed( QObject::tr("NOT NEARLY EQUAL\n\n%1\n\n!=\n\n%2\n\n"
                                                         "TO WITHIN %3")
                        .arg(Sire::toString(a)).arg(Sire::toString(b))
                        .arg(Sire::toString(range)), code_location );
    }
}

/** Assert that the passed expression is true */
SIRE_ALWAYS_INLINE void assert_true( bool result, const QString &code_location )
{
    if (not result)
        throw SireError::assertation_failed( QObject::tr(
                "The expression is not true"), code_location );
}

/** Assert that the passed expression is false */
SIRE_ALWAYS_INLINE void assert_false( bool result, const QString &code_location )
{
    if (result)
        throw SireError::assertation_failed( QObject::tr(
                "The expression is not false"), code_location );
}

#ifdef SIRE_HAS_CPP_11
    /** Assert that calling the passed function will raise the supplied exception */
    SIRE_ALWAYS_INLINE void assert_throws( std::function<void ()> func, const SireError::exception &e,
                               const QString &code_location )
    {
        try
        {
            func();
        }
        catch(const SireError::exception &thrown)
        {
            if ( QLatin1String(thrown.what()) != QLatin1String(e.what()) )
            {
                throw SireError::assertation_failed( QObject::tr(
                        "The expected exception %1 was not thrown. Instead the exception "
                        "%2 was thrown.").arg(e.what()).arg(thrown.what()), code_location );
            }
            
            return;
        }
        catch(const std::exception &thrown)
        {
            throw SireError::assertation_failed( QObject::tr(
                    "The expected exception %1 was not thrown. Instead the exception "
                    "%2 was thrown.").arg(e.what()).arg(thrown.what()), code_location );
        }
        catch(...)
        {
                throw SireError::assertation_failed( QObject::tr(
                        "The expected exception %1 was not thrown. Instead an unknown exception "
                        "was thrown.").arg(e.what()), code_location );
        }
        
        throw SireError::assertation_failed( QObject::tr(
                    "The expected exception %1 was not thrown, as no exception was thrown!")
                        .arg(e.what()), code_location );
    }
#endif


}

SIRE_EXPOSE_CLASS( SireBase::UnitTest )

#endif
