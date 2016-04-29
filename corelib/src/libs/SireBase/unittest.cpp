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

#include "unittest.h"

#include "SireError/exception.h"

#include <QElapsedTimer>
#include <QTextStream>

using namespace SireBase;

static QMutex global_mutex;

QList< boost::shared_ptr<UnitTest> > UnitTest::all_tests;

/** Null constructor */
UnitTest::UnitTest() : test_function(0), run_time(0)
{}

/** Construct a test called 'name' that will call the passed function. 
    The test function must have the signature
    
    void test_function(bool verbose)
    {}
    
    where 'verbose' tells the function whether or not to be verbose
    or silent. The test should signal failure by throwing
    an exception
*/
UnitTest::UnitTest(const QString &name, void (*test_func)(bool))
         : test_name(name), test_function(test_func), run_time(0)
{
    //register this test with the global list
    if (test_func != 0)
    {
        QMutexLocker lkr(&global_mutex);
        all_tests.append( boost::shared_ptr<UnitTest>(new UnitTest(*this)) );
    }
}

/** Copy constructor */
UnitTest::UnitTest(const UnitTest &other)
         : test_name(other.test_name), test_function(other.test_function),
           error_string(other.error_string), run_time(other.run_time)
{}

/** Destructor */
UnitTest::~UnitTest()
{}

/** Return the name of the test */
QString UnitTest::name() const
{
    return test_name;
}

static QTextStream cerr( stderr, QIODevice::WriteOnly | QIODevice::Unbuffered );

/** Run this test a total of 'nrepeats' times, printing out information
    if 'verbose' is true */
bool UnitTest::run(int nrepeats, bool verbose)
{
    QMutexLocker lkr(&run_mutex);

    if (nrepeats <= 0)
        return true;
    
    try
    {
        QElapsedTimer t;
        t.start();
        for (int i=0; i<nrepeats; ++i)
        {
            (*test_function)(verbose);
            cerr << ".";
            cerr.flush();
        }
        quint64 ns = t.nsecsElapsed();
        
        run_time = ns / nrepeats;
        error_string = QString::null;
        return true;
    }
    catch(const SireError::exception &e)
    {
        error_string = e.toString();
    }
    catch(const std::exception &e)
    {
        error_string = QObject::tr("std::exception( %1 )").arg(e.what());
    }
    catch(...)
    {
        error_string = QObject::tr("An unknown error occured!");
    }
    
    return false;
}

/** Return whether or not the test was successful */
bool UnitTest::wasSuccessful()
{
    QMutexLocker lkr(&run_mutex);
    return error_string.isEmpty();
}

/** Return whether or not the test was an error */
bool UnitTest::wasError()
{
    return not wasSuccessful();
}

/** Return the error string */
QString UnitTest::errorString()
{
    QMutexLocker lkr(&run_mutex);
    return error_string;
}

/** Return the time taken to run the test, in nanoseconds */
quint64 UnitTest::runTime()
{
    QMutexLocker lkr(&run_mutex);
    return run_time;
}

/** Return all of the tests that have been registered */
QList< boost::shared_ptr<UnitTest> > UnitTest::tests()
{
    QMutexLocker lkr(&global_mutex);
    return all_tests;
}

/** Run all of the tests one after another, printing out the results to the
    screen and returning the number of tests that failed */
int UnitTest::runAll(bool verbose)
{
    QMutexLocker lkr(&global_mutex);
    
    cerr << QObject::tr("\n== Running Unit Tests ==\n\n");
    cerr.flush();
    
    int i = 1;
    const int n = all_tests.count();
    int nfailed = 0;
    
    foreach( boost::shared_ptr<UnitTest> test, all_tests )
    {
        cerr << QObject::tr("(%1/%2) %3 ").arg(i).arg(n).arg(test->name());
        cerr.flush();
        bool passed = test->run(1, verbose);

        if (passed and not verbose)
        {
            //if the test was really quick, then run it several times to try different inputs
            if (test->runTime() < 500000000)
            {
                passed = test->run( qMin( int(500000000 / test->runTime()), 10 ), verbose );
            }
        }
        
        if (passed)
        {
            cerr << QObject::tr("(passed - took %1 ms)\n").arg( 0.000001 * test->runTime() );
            cerr.flush();
        }
        else
        {
            cerr << QObject::tr("(failed)\n");
            cerr.flush();
            nfailed += 1;
        }
        
        i += 1;
    }
    
    if (nfailed > 0)
    {
        cerr << QObject::tr("\n\nA number of tests (%1) failed!!!\n").arg(nfailed);
        
        foreach( boost::shared_ptr<UnitTest> test, all_tests )
        {
            if (test->wasError())
            {
                cerr << QObject::tr("\n== FAILED TEST | %1 ==\n%2\n")
                            .arg(test->name()).arg(test->errorString());
            }
        }
        
        cerr << QObject::tr("\nSOME TESTS FAILED!!! (%1)\n").arg(nfailed);
    }
    else
        cerr << QObject::tr("\nALL TESTS PASSED :-)\n");
    
    return nfailed;
}
