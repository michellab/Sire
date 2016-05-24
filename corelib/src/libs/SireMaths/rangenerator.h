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

#ifndef SIREMATHS_RANGENERATOR_H
#define SIREMATHS_RANGENERATOR_H

#include <boost/shared_ptr.hpp>
#include <boost/noncopyable.hpp>

#include <QVector>

#include "sireglobal.h"

SIRE_BEGIN_HEADER

namespace SireMaths
{
class RanGenerator;
}

QDataStream& operator<<(QDataStream&, const SireMaths::RanGenerator&);
QDataStream& operator>>(QDataStream&, SireMaths::RanGenerator&);

namespace SireMaths
{

class Vector;

namespace detail
{
class RanGeneratorPvt;
}

void seed_qrand();

/** This class provides a thread-safe, copyable and streamable
    random number generator. Copies are guaranteed to produce
    different random number sequences (thus the possibility
    of accidental repeat random numbers is removed).

    @author Christopher Woods
*/
class SIREMATHS_EXPORT RanGenerator
{

friend QDataStream& ::operator<<(QDataStream&, const RanGenerator&);
friend QDataStream& ::operator>>(QDataStream&, RanGenerator&);

public:
    RanGenerator();
    RanGenerator(quint32 seed);
    RanGenerator(const QVector<quint32> &seed);

    RanGenerator(const RanGenerator &other);

    ~RanGenerator();

    static const char* typeName();
    
    const char* what() const
    {
        return RanGenerator::typeName();
    }

    RanGenerator& operator=(const RanGenerator &other);

    bool operator==(const RanGenerator &other) const;
    bool operator!=(const RanGenerator &other) const;

    void detach();

    void seed();
    void seed(quint32 seed);
    void seed(const QVector<quint32> &seed);
    void seed(const RanGenerator &other);

    double rand() const;
    double rand(double maxval) const;
    double rand(double minval, double maxval) const;

    double locked_rand() const;
    double locked_rand(double maxval) const;
    double locked_rand(double minval, double maxval) const;

    QVector<double> nrand(int n) const;
    QVector<double> nrand(int n, double maxval) const;
    QVector<double> nrand(int n, double minval, double maxval) const;
    
    void nrand(QVector<double> &result) const;
    void nrand(QVector<double> &result, double maxval) const;
    void nrand(QVector<double> &result, double minval, double maxval) const;

    double rand53() const;
    double rand53(double maxval) const;
    double rand53(double minval, double maxval) const;

    double randNorm() const;
    double randNorm(double mean, double variance) const;

    double locked_randNorm() const;
    double locked_randNorm(double minval, double maxval) const;

    void nrandNorm(QVector<double> &result, double mean, double variance) const;
    QVector<double> nrandNorm(int n, double mean, double variance) const;

    Vector vectorOnSphere() const;
    Vector vectorOnSphere(double radius) const;

    Vector locked_vectorOnSphere() const;
    Vector locked_vectorOnSphere(double radius) const;

    void nvectorOnSphere(QVector<Vector> &result) const;
    void nvectorOnSphere(QVector<Vector> &result, double radius) const;
    
    QVector<Vector> nvectorOnSphere(int n) const;
    QVector<Vector> nvectorOnSphere(int n, double radius) const;

    bool randBool() const;

    quint32 randInt() const;
    quint32 randInt(quint32 maxval) const;
    qint32 randInt(qint32 minval, qint32 maxval) const;

    quint64 randInt64() const;
    quint64 randInt64(quint64 maxval) const;
    qint64 randInt64(qint64 minval, qint64 maxval) const;

    QVector<quint32> getState() const;
    void setState(const QVector<quint32> &state);
    
    void lock() const;
    void unlock() const;
    
    static const RanGenerator& global();
    
    static void seedGlobal();
    static void seedGlobal(quint32 seed);
    static void seedGlobal(const QVector<quint32> &seed);
    static void seedGlobal(const RanGenerator &other);

private:
    detail::RanGeneratorPvt& nonconst_d() const;

    /** Shared pointer to the actual generator */
    boost::shared_ptr<detail::RanGeneratorPvt> d;
};


/** This is a small locker class that holds a lock on the random
    number generator */
class RanGeneratorLocker : boost::noncopyable
{
public:
    RanGeneratorLocker() : boost::noncopyable(), rangen(0)
    {}
    
    RanGeneratorLocker(const RanGenerator *r) : boost::noncopyable(), rangen(r)
    {
        if (rangen)
            rangen->lock();
    }
    
    ~RanGeneratorLocker()
    {
        if (rangen)
        {
            rangen->unlock();
            rangen = 0;
        }
    }
    
    void unlock()
    {
        if (rangen)
        {
            rangen->unlock();
            rangen = 0;
        }
    }
    
private:
    const RanGenerator *rangen;
};

}

Q_DECLARE_METATYPE(SireMaths::RanGenerator);

SIRE_EXPOSE_CLASS( SireMaths::RanGenerator )

SIRE_END_HEADER

#endif
