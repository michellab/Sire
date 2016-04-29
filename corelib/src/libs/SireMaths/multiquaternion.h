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

#ifndef SIREMATHS_MULTIQUATERNION_H
#define SIREMATHS_MULTIQUATERNION_H

#include "quaternion.h"
#include "multivector.h"

SIRE_BEGIN_HEADER

namespace SireMaths
{

class MultiVector;
class MultiQuaternion;

MultiQuaternion operator+(const MultiQuaternion &p1, const MultiQuaternion &p2);
MultiQuaternion operator-(const MultiQuaternion &p1, const MultiQuaternion &p2);
MultiQuaternion operator*(const MultiQuaternion &p1, const MultiQuaternion &p2);

/**
This is the MultiX version of Quaternion
 
@author Christopher Woods
*/
class SIREMATHS_EXPORT MultiQuaternion
{
public:
    MultiQuaternion();
    MultiQuaternion(const MultiQuaternion& p1);
    
    MultiQuaternion(const MultiDouble &angle, const MultiVector &axis);
    MultiQuaternion(const MultiFloat &angle, const MultiVector &axis);
    MultiQuaternion(const MultiDouble &x, const MultiDouble &y,
                    const MultiDouble &z, const MultiDouble &w);
    
    ~MultiQuaternion();
    
    static const char* typeName();
    
    const char* what() const
    {
        return MultiQuaternion::typeName();
    }

    const MultiDouble& x() const;
    const MultiDouble& y() const;
    const MultiDouble& z() const;
    const MultiDouble& w() const;

    MultiQuaternion inverse() const;
    MultiQuaternion conjugate() const;

    MultiDouble dot(const MultiQuaternion &q) const;

    QString toString() const;

    MultiVector rotate(const MultiVector &p) const;
    QVector<MultiVector> rotate(const QVector<MultiVector> &points) const;

    MultiQuaternion slerp(const MultiQuaternion &q, const MultiDouble &lam) const;

    MultiQuaternion pow(const MultiDouble &n) const;

    static MultiQuaternion identity();
    
    void renormalise();

    bool operator==(const MultiQuaternion &p1) const;
    bool operator!=(const MultiQuaternion &p1) const;

    static int size();
    static int count();

    void set(int i, const Quaternion &val);
    void quickSet(int i, const Quaternion &val);

    Quaternion operator[](int i) const;
    Quaternion at(int i) const;
    Quaternion getitem(int i) const;

    MultiQuaternion& operator=(const MultiQuaternion &p);

    MultiQuaternion& operator+=(const MultiQuaternion &p);
    MultiQuaternion& operator-=(const MultiQuaternion &p);
    MultiQuaternion& operator*=(const MultiQuaternion &p);
    MultiQuaternion& operator*=(const MultiVector &p);

    friend MultiQuaternion operator+(const MultiQuaternion &p1, const MultiQuaternion &p2);
    friend MultiQuaternion operator-(const MultiQuaternion &p1, const MultiQuaternion &p2);
    friend MultiQuaternion operator*(const MultiQuaternion &p1, const MultiQuaternion &p2);
    friend MultiQuaternion operator*(const MultiQuaternion &p1, const MultiVector &p2);
    friend MultiQuaternion operator*(const MultiVector &p1, const MultiQuaternion &p2);

    void swap(MultiQuaternion &q0, int idx0, MultiQuaternion &q1, int idx1);

private:
    /** The x,y,z,w of the quaternion */
    MultiDouble sc[4];
};

}

SIRE_END_HEADER

#endif

