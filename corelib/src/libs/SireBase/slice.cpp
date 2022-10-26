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

#include "slice.h"

#include "SireError/errors.h"

#include <limits>

#include <QDebug>

using namespace SireBase;

const int _unset = std::numeric_limits<int>::max();

int _map(int val, int n, bool auto_fix)
{
    if (val < 0)
        val = n + val;

    if (auto_fix)
    {
        if (val >= n)
            val = n-1;

        if (val < 0)
            val = 0;
    }
    else if (val < 0 or val >= n)
        throw SireError::invalid_index(QObject::tr(
            "Invalid index for a container with count == %1"
                    ).arg(n), CODELOC);

    return val;
}

SliceIterator::SliceIterator()
              : i(0), start(0), stop(0), step(0)
{}

SliceIterator::SliceIterator(const Slice &slice, int n, bool auto_fix)
              : i(slice.start), start(slice.start),
                stop(slice.stop), step(slice.step)
{
    if (n <= 0)
        throw SireError::invalid_index(
            QObject::tr("Cannot slice an empty container!"), CODELOC);

    start = slice.start;

    if (start == _unset)
        start = 0;

    start = _map(start, n, auto_fix);
    stop = slice.stop;
    step = slice.step;

    if (step == _unset)
    {
        step = 1;
    }

    if (stop == _unset)
    {
        if (step >= 0)
            stop = n-1;
        else
            stop = 0;
    }
    else
    {
        if (step >= 0)
        {
            if (slice.start == slice.stop)
            {
                // this is a single-value slice
                stop = start;
            }
            else if (stop != 0)
                stop = _map(stop-1, n, auto_fix);
        }
        else
        {
            if (slice.start == slice.stop)
            {
                // this is a single-value slice
                stop = start;
            }
            else
                stop = _map(stop+1, n, auto_fix);
        }
    }

    if (step == 0)
        step = 1;

    if (start <= stop and step <= 0)
    {
        step *= -1;
    }
    else if (start > stop and step >= 0)
    {
        step *= -1;
    }

    i = start;
}

SliceIterator::SliceIterator(const SliceIterator &other)
              : i(other.i), start(other.start),
                stop(other.stop), step(other.step)
{}

SliceIterator::~SliceIterator()
{}

SliceIterator& SliceIterator::next()
{
    i += step;
    return *this;
}

bool SliceIterator::atEnd() const
{
    if (step > 0)
        return i > stop;
    else
        return i < stop;
}

int SliceIterator::value() const
{
    return i;
}

Slice::Slice() : start(_unset), stop(_unset), step(_unset)
{}

Slice::Slice(const Slice &other)
      : start(other.start), stop(other.stop), step(other.step)
{}

Slice::~Slice()
{}

Slice Slice::fromStartStop(int _start, int _stop, int _step)
{
    Slice s;
    s.start = _start;
    s.stop = _stop;
    s.step = _step;

    return s;
}

Slice Slice::fromStart(int _start, int _step)
{
    Slice s;
    s.start = _start;
    s.stop = _unset;
    s.step = _step;

    return s;
}

Slice& Slice::operator=(const Slice &other)
{
    start = other.start;
    stop = other.stop;
    step = other.step;
    return *this;
}

bool Slice::operator==(const Slice &other) const
{
    return start == other.start and stop == other.stop and step == other.step;
}

bool Slice::operator!=(const Slice &other) const
{
    return not operator==(other);
}

const char* Slice::what() const
{
    return Slice::typeName();
}

const char* Slice::typeName()
{
    return "SireBase::Slice";
}

int Slice::unset()
{
    return _unset;
}

QString Slice::toString() const
{
    if (start == _unset)
    {
        if (stop == _unset)
        {
            if (step == _unset)
                return QObject::tr("slice[::]");
            else
                return QObject::tr("slice[::%1]").arg(step);
        }
        else if (step == _unset)
        {
            return QObject::tr("slice[:%1]").arg(stop);
        }
        else
        {
            return QObject::tr("slice[:%1:%2]").arg(stop).arg(step);
        }
    }
    else if (stop == _unset)
    {
        if (step == _unset)
            return QObject::tr("slice[%1:]").arg(start);
        else
            return QObject::tr("slice[%1::%2]").arg(start).arg(step);
    }
    else
    {
        if (step == _unset)
            return QObject::tr("slice[%1:%2]").arg(start).arg(stop);
        else
            return QObject::tr("slice[%1:%2:%3]").arg(start).arg(stop).arg(step);
    }
}

SliceIterator Slice::begin(int n, bool auto_fix) const
{
    return SliceIterator(*this, n, auto_fix);
}
