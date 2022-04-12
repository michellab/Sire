/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2022  Christopher Woods
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

#ifndef SIREBASE_SLICE_H
#define SIREBASE_SLICE_H

#include "sireglobal.h"

SIRE_BEGIN_HEADER

namespace SireBase
{

class SliceIterator;

/** Simple class that simplifies accessing Sire containers
    via Python slices
*/
class SIREBASE_EXPORT Slice
{

friend class SliceIterator;

public:
    Slice();
    Slice(const Slice &other);

    ~Slice();

    static Slice fromStartStop(int start, int stop, int step=1);
    static Slice fromStart(int start, int step=1);

    Slice& operator=(const Slice &other);

    bool operator==(const Slice &other) const;
    bool operator!=(const Slice &other) const;

    const char* what() const;
    static const char* typeName();

    QString toString() const;

    SliceIterator begin(int n) const;

private:
    int start;
    int stop;
    int step;
};

/** A simple iterator for the slice */
class SIREBASE_EXPORT SliceIterator
{
public:
    SliceIterator();
    SliceIterator(const Slice &slice, int n);
    SliceIterator(const SliceIterator &other);

    ~SliceIterator();

    SliceIterator& next();

    bool atEnd() const;

    int value() const;

private:
    int i;
    int start;
    int stop;
    int step;
};

}

SIRE_END_HEADER

#endif
