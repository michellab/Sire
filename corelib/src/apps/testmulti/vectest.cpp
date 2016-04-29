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

#include "SireMaths/multivector.h"
#include "SireMaths/multiquaternion.h"

#include "SireMaths/rangenerator.h"
#include "SireMaths/quaternion.h"
#include "SireError/errors.h"
#include "tostring.h"

#include <QVector>
#include <QDebug>
#include <iostream>

using std::cerr;
using std::endl;

using namespace SireMaths;
using namespace SireError;

static RanGenerator rangen;

template<class T>
void assert_equal( const T &a, const T &b )
{
    if (a != b)
    {
        qDebug() << "\n\nASSERT EQUAL FAILED!" << Sire::toString(a)
                 << Sire::toString(b);
        throw assertation_failed( QObject::tr("NOT EQUAL! %1 != %2")
			.arg(Sire::toString(a)).arg(Sire::toString(b)), CODELOC );
    }
}

template<class T>
void assert_nearly_equal( const T &a, const T &b, const T &range )
{
    if ((a - b < -range) or (a - b > range))
    {                                    
        qDebug() << "ASSERT NEARLY EQUAL FAILED!" << Sire::toString(a)
                 << Sire::toString(b);
        throw assertation_failed( QObject::tr("NOT NEARLY EQUAL %1 != %2 TO WITHIN %3")
                        .arg(Sire::toString(a)).arg(Sire::toString(b))
                        .arg(Sire::toString(range)), CODELOC );
    }
}

void test_vecarray()
{
    //generate 1000 vectors
    QVector<Vector> vecs(1000);

    for (int i=0; i<1000; ++i)
    {
        vecs[i] = rangen.vectorOnSphere(10);
    }

    //convert to a MultiVector
    QVector<MultiVector> mvecs = MultiVector::fromArray(vecs);

    assert_equal( mvecs.count(), (vecs.count() / MultiVector::count() ) );

    for (int i=0; i<1000; ++i)
    {
        assert_equal( vecs.at(i), mvecs.at(i / MultiVector::count())
                                       .at(i % MultiVector::count()) );
    }
}

void test_length()
{
    QVector<Vector> v(MultiVector::count());

    for (int i=0; i<MultiVector::count(); ++i)
    {
        v[i] = Vector( rangen.rand(-10,10), rangen.rand(-10,10),
                       rangen.rand(-10,10) );
    }

    MultiVector mv(v);
    MultiDouble l = mv.length();

    for (int i=0; i<MultiVector::count(); ++i)
    {
        assert_equal( l.at(i), v[i].length() );
    }

    mv = mv.normalise();
    l = mv.length();

    for (int i=0; i<MultiVector::count(); ++i)
    {
        assert_nearly_equal( l.at(i), double(1.0), 0.0001 );
    }
}

void test_cross()
{
    QVector<Vector> v0(MultiVector::count()), v1(MultiVector::count());
    QVector<Vector> cross(MultiVector::count());

    for (int i=0; i<MultiVector::count(); ++i)
    {
        v0[i] = rangen.vectorOnSphere(5);
        v1[i] = rangen.vectorOnSphere(5);
        cross[i] = Vector::cross(v0[i], v1[i]);
    }

    MultiVector mv0(v0);
    MultiVector mv1(v1);
    MultiVector mcross = MultiVector::cross(mv0, mv1);

    for (int i=0; i<MultiVector::count(); ++i)
    {
        assert_equal( mcross.at(i), cross.at(i) );
    }
}

void test_getset()
{
    QVector<Vector> v(MultiVector::count());

    for (int i=0; i<MultiVector::count(); ++i)
    {
        v[i] = Vector( rangen.rand(-10,10), rangen.rand(-10,10),
                       rangen.rand(-10,10) );
    }

    MultiVector mv;

    for (int i=0; i<MultiVector::count(); ++i)
    {
        mv.set(i, v[i]);
        assert_equal( mv.getitem(i), v[i] );
    }

    assert_equal( MultiVector(v), mv );
}

void test_rotate()
{
    QVector<Vector> axis(MultiVector::count()), vec(MultiVector::count());
    QVector<double> angle(MultiVector::count());
    QVector<Vector> result(MultiVector::count());

    for (int i=0; i<MultiVector::count(); ++i)
    {
        axis[i] = rangen.vectorOnSphere(1);
        vec[i] = rangen.vectorOnSphere(5);
        angle[i] = rangen.rand(-2.0*SireMaths::pi, 2.0*SireMaths::pi);
        result[i] = Quaternion(Angle(angle[i]),axis[i]).rotate(vec[i]);
    }

    MultiVector maxis(axis);
    MultiVector mvec(vec);
    MultiDouble mangle(angle);

    MultiVector mresult = MultiQuaternion(mangle, maxis).rotate(mvec);

    for (int i=0; i<MultiVector::count(); ++i)
    {
        for (int j=0; j<3; ++j)
        {
            assert_nearly_equal( mresult.at(i)[j], result[i][j], 1e-5 );
        }
    }
}

int main(int argc, const char **argv)
{
    try
    {
        cerr << "test_getset";
        for (int i=0; i<100; ++i){ test_getset(); cerr << "."; }
        cerr << "(passed)\n";

        cerr << "test_vecarray...";
        test_vecarray();
        cerr << "(passed)\n";

        cerr << "test_cross";
        for (int i=0; i<100; ++i){ test_cross(); cerr << "."; }
        cerr << "(passed)\n";

        cerr << "test_length";
        for (int i=0; i<100; ++i){ test_length(); cerr << "."; }
        cerr << "(passed)\n";

        cerr << "test_rotate";
        for (int i=0; i<100; ++i){ test_rotate(); cerr << "."; }
        cerr << "(passed)\n";
    }
    catch(const SireError::exception &e)
    {
        qDebug() << e.toString();
        return -1;
    }

    return 0;
}
