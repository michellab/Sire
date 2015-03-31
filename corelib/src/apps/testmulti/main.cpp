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

#include "SireMaths/multifloat.h"
#include "SireMaths/multidouble.h"
#include "SireMaths/multiint.h"

#include <QDebug>

//#include <type_traits>

using namespace SireMaths;

void intTests()
{
    const int size = 100;

    QVector<qint32> vals(size);

    for (int i=0; i<size; ++i)
    {
        vals[i] = i;
    }

    QVector<MultiInt> mi = MultiInt::fromArray(vals);

    for (int i=0; i<mi.count(); ++i)
    {
        qDebug() << mi[i].toString() << mi[i].toBinaryString();
    }

    MultiInt m = mi[0];

    qDebug() << "ROTATE";
    for (int i=0; i<MultiInt::size(); ++i)
    {
        qDebug() << i << m.toString();
        m = m.rotate();
    }

    qDebug() << "CAST";
    qDebug() << m.reinterpretCastToFloat().toString();
    qDebug() << m.reinterpretCastToFloat().toBinaryString();

    qDebug() << "ADD";
    qDebug() << (m+m).toString();

    qDebug() << "SUB";
    qDebug() << (m-m).toString();

    qDebug() << "COMPARE ==";
    qDebug() << (m.compareEqual(m)).toBinaryString();

    qDebug() << "COMPARE !=";
    qDebug() << (m.compareNotEqual(m)).toBinaryString();

    qDebug() << "COMPARE >";
    qDebug() << (m.compareGreater(m.rotate())).toBinaryString();

    qDebug() << "COMPARE <";
    qDebug() << (m.rotate().compareLess(m)).toBinaryString();

    qDebug() << "AND";
    qDebug() << m.compareEqual(m).logicalAnd(m).toBinaryString();

    qDebug() << "OR";
    qDebug() << m.compareEqual(m).logicalOr(m).toBinaryString();
}

int main(int argc, const char **argv)
{
    QVector<float> vals;

    for (int i=1; i<=MultiFloat::size(); ++i)
    {
        vals.append( i );
    }

    MultiFloat mf(vals);
    MultiDouble md(vals);

    MultiFloat mf2(md);
    MultiDouble md2(mf);

    qDebug() << mf.toString();
    qDebug() << md.toString();    
    qDebug() << mf2.toString();
    qDebug() << md2.toString();

    mf2 = md;
    md2 = mf;

    qDebug() << mf2.toString();
    qDebug() << md2.toString();

    qDebug() << "ROTATE";
    for (int i=1; i<=MultiFloat::size(); ++i)
    {
        mf = mf.rotate();
        md = md.rotate();
        qDebug() << i << mf.toString();
        qDebug() << i << md.toString();
    }

    qDebug() << "SQRT";
    qDebug() << mf.sqrt().toString();
    qDebug() << md.sqrt().toString();

    for (int i=MultiFloat::count()+1; i<= 101; ++i)
    {
        vals.append(i);
    }

    QVector<MultiFloat> mfarray = MultiFloat::fromArray(vals);
    QVector<MultiDouble> mdarray = MultiDouble::fromArray(vals);

    for (int i=0; i<mfarray.count(); ++i)
    {
        qDebug() << mfarray[i].toString();
        qDebug() << mdarray[i].toString();
    }

    qDebug() << "\nMultiFloat tests";
    MultiFloat a(1);
    MultiFloat b(2);

    qDebug() << a.compareNotEqual(a).toBinaryString();
    qDebug() << a.compareNotEqual(b).toBinaryString();

    qDebug() << "\nAbs tests";

    for (int i=0; i<MultiFloat::size(); ++i)
    {
        vals[i] = 0.1 * i;
    }

    a = MultiFloat(vals.constData(), MultiFloat::size());

    for (int i=0; i<MultiFloat::size(); ++i)
    {
        vals[i] = -0.1 * i;
    }

    b = MultiFloat(vals.constData(), MultiFloat::size());

    qDebug() << a.toString() << a.abs().toString();
    qDebug() << b.toString() << b.abs().toString();

    /*
    qDebug() << "\nMultiInt tests";
    intTests();
    */

    return 0;

}
