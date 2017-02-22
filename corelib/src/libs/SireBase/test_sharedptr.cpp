/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2017  Christopher Woods
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

#include "SireBase/parallel.h"
#include "SireBase/unittest.h"

#include "SireBase/shareddatapointer.hpp"
#include "SireBase/sharedpolypointer.hpp"

#include <cmath>
#include <QSharedData>

using namespace SireBase;

class Foo : public QSharedData
{
public:
    Foo() : QSharedData()
    {
        for (int i=0; i<size(); ++i)
        {
            data[i] = std::rand();
        }
    }
    
    Foo(const Foo &other) : QSharedData()
    {
        for (int i=0; i<size(); ++i)
        {
            data[i] = other.data[i];
        }
    }
    
    virtual ~Foo()
    {}
    
    Foo& operator=(const Foo &other)
    {
        for (int i=0; i<size(); ++i)
        {
            data[i] = other.data[i];
        }
        
        return *this;
    }
    
    bool operator==(const Foo &other) const
    {
        for (int i=0; i<size(); ++i)
        {
            if (data[i] != other.data[i])
                return false;
        }
        
        return true;
    }
    
    bool operator!=(const Foo &other) const
    {
        return not operator==(other);
    }
    
    int pointer() const
    {
        return qintptr(this);
    }
    
    static int size()
    {
        return 24;
    }
    
    virtual Foo* clone() const
    {
        return new Foo(*this);
    }
    
    QString toString() const
    {
        return QObject::tr("Foo( ptr=%1 )").arg(qintptr(this));
    }
    
private:
    float data[24];
};

typedef SharedDataPointer<Foo> FooPtr;

void test_sharedptr(bool verbose)
{
    FooPtr f1 = Foo();
    
    FooPtr f2 = f1;
    
    assert_equal( f1.read(), f2.read(), CODELOC );
    assert_equal( f1.read().pointer(), f2.read().pointer(), CODELOC );
    
    f2.detach();
    
    assert_equal( f1.read(), f2.read(), CODELOC );
    assert_not_equal( f1.read().pointer(), f2.read().pointer(), CODELOC );

    const int nfoo = 2048;
    
    QVector<FooPtr> foos(nfoo);
    
    tbb::parallel_for( tbb::blocked_range<int>(0,nfoo),
                       [&](tbb::blocked_range<int> r)
    {
        for (int i=r.begin(); i<r.end(); ++i)
        {
            foos[i] = Foo();
        }
    });

    QVector<FooPtr> foos2(foos);
    foos2.detach();
    
    for (int i=0; i<nfoo; ++i)
    {
        assert_equal( foos[i].read(), foos2[i].read(), CODELOC );
        assert_equal( foos[i].read().pointer(), foos2[i].read().pointer(), CODELOC );
    }
    
    QVector<FooPtr> foos3(foos);
    
    for (int i=0; i<nfoo; ++i)
    {
        assert_equal( foos[i].read(), foos3[i].read(), CODELOC );
        assert_equal( foos3[i].read().pointer(), foos3[i].read().pointer(), CODELOC );
    }
    
    foos.detach();
    
    tbb::parallel_for( tbb::blocked_range<int>(0,nfoo),
                        [&](tbb::blocked_range<int> r)
    {
        for (int i=r.begin(); i<r.end(); ++i)
        {
            foos[i].write();
        }
    });
    
    for (int i=0; i<nfoo; ++i)
    {
        assert_equal( foos[i].read(), foos3[i].read(), CODELOC );
        assert_not_equal( foos[i].read().pointer(), foos3[i].read().pointer(), CODELOC );
    }
    
}

SIRE_UNITTEST( test_sharedptr )
