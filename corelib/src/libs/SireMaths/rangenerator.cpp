/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2007  Christopher Woods
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

#include <QMutex>
#include <QVector>
#include <QUuid>

#include <QDebug>

#include <boost/noncopyable.hpp>
#include <boost/scoped_array.hpp>

#include <limits>

#include "rangenerator.h"
#include "vector.h"

#include "ThirdParty/MersenneTwister.h"       // CONDITIONAL_INCLUDE

#include "SireBase/refcountdata.h"

#include "SireError/errors.h"
#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMaths;
using namespace SireMaths::detail;

using namespace SireStream;

namespace SireMaths
{
namespace detail
{
class RanGeneratorPvt;
}
}

QDataStream& operator<<(QDataStream&, const SireMaths::detail::RanGeneratorPvt&);
QDataStream& operator>>(QDataStream&, SireMaths::detail::RanGeneratorPvt&);

namespace SireMaths
{

namespace detail
{

/** This class provides the private implementation of RanGenerator.
    This class is explicitly shared by instances of RanGenerator,
    and is therefore thread-safe.

    @author Christopher Woods
*/
class RanGeneratorPvt
{

friend QDataStream& ::operator<<(QDataStream&, const RanGeneratorPvt&);
friend QDataStream& ::operator>>(QDataStream&, RanGeneratorPvt&);

typedef MTRand::uint32 MTUInt32;

public:
    /** Construct a generator with a random seed */
    RanGeneratorPvt() : mutex(QMutex::NonRecursive)
    {}

    /** Construct a generator with a specified seed */
    RanGeneratorPvt(quint32 seed) : mutex(QMutex::NonRecursive),
                                    mersenne_generator(seed)
    {}

    /** Construct a generator with a specified seed */
    RanGeneratorPvt(const QVector<quint32> &s) : mutex(QMutex::NonRecursive)
    {
        this->seed(s);
    }

    /** Copy constructor */
    RanGeneratorPvt(const RanGeneratorPvt &other) : mutex(QMutex::NonRecursive)
    {
        QMutexLocker lkr( const_cast<QMutex*>( &(other.mutex) ) );
        
        mersenne_generator = other.mersenne_generator;
    }

    /** Destructor */
    ~RanGeneratorPvt()
    {}

    /** Mutex to serialise access to the generator */
    QMutex mutex;

    /** The actual generator (Mersenne Twister) */
    MTRand mersenne_generator;

    /** Randomly seed the generator */
    void seed()
    {
        QMutexLocker lkr(&mutex);

        mersenne_generator.seed();
    }

    /** Reseed the generator from an array of uints */
    void seed(const QVector<quint32> &s)
    {
        if (s.isEmpty())
        {
            this->seed();
            return;
        }

        QMutexLocker lkr(&mutex);

        //we need to convert this into an array of MTRand::uint32
        int sz = s.count();

        boost::scoped_array<MTUInt32> array( new MTUInt32[sz] );

        for (int i=0; i<sz; ++i)
            array[i] = s.constData()[i];

        mersenne_generator.seed( array.get(),  sz );
    }

    /** Reseed the generator */
    void seed(quint32 s)
    {
        QMutexLocker lkr(&mutex);

        mersenne_generator.seed(s);
    }

    /** Return an array containing the state of the random generator */
    QVector<quint32> getState()
    {
        mutex.lock();

        //create an array to hold the state
        int state_size = MTRand::N + 1;
        boost::scoped_array<MTUInt32> array( new MTUInt32[state_size] );

        mersenne_generator.save(array.get());

        mutex.unlock();

        //copy the array into a QVector<quint32>
        QVector<quint32> ret(state_size);

        quint32 *ret_array = ret.data();

        for (int i=0; i<state_size; ++i)
            ret_array[i] = array[i];

        return ret;
    }

    /** Load the state from an array - the array must have size
        MTRand::N + 1 (625)

        \throw SireError::incompatible_error
    */
    void loadState(const QVector<quint32> &state)
    {
        int state_size = MTRand::N + 1;
    
        //check that the array is of the right size...
        if (state.count() != state_size)
            throw SireError::incompatible_error( QObject::tr(
                "Can only restore the state from an array of size %1, "
                "while you have provided an array of size %2.")
                    .arg(state_size).arg(state.count()), CODELOC );

        //convert the array to the right type...
        boost::scoped_array<MTUInt32> array( new MTUInt32[state_size] );

        const quint32 *state_array = state.constData();

        for (int i=0; i<state_size; ++i)
            array[i] = state_array[i];

        QMutexLocker lkr(&mutex);

        mersenne_generator.load(array.get());
    }
};

} // end of namespace detail

} // end of namespace SireMaths

////////////
//////////// Implementation of RanGeneratorPvt
////////////

QDataStream& operator<<(QDataStream &ds, 
                        const SireMaths::detail::RanGeneratorPvt &rangen)
{
    ds << const_cast<RanGeneratorPvt*>(&rangen)->getState();
    return ds;
}

QDataStream& operator>>(QDataStream &ds,
                        SireMaths::detail::RanGeneratorPvt &rangen)
{
    QVector<quint32> state;
    ds >> state;
     
    rangen.loadState(state);
    
    return ds;
}

////////////
//////////// Implementation of RanGenerator
////////////

static const RegisterMetaType<RanGenerator> r_rangen(NO_ROOT);

/** Serialise to a binary data stream */
QDataStream &operator<<(QDataStream &ds, const RanGenerator &rangen)
{
    writeHeader(ds, r_rangen, 1);

    SharedDataStream sds(ds);
    sds << rangen.d;

    return ds;
}

/** Deserialise from a binary data stream */
QDataStream &operator>>(QDataStream &ds, RanGenerator &rangen)
{
    VersionID v = readHeader(ds, r_rangen);

    if (v == 1)
    {
        //I need to detach from shared storage
        rangen.d.reset(new RanGeneratorPvt());
    
        SharedDataStream sds(ds);
        sds >> rangen.d;
    }
    else
        throw version_error(v, "1", r_rangen, CODELOC);

    return ds;
}

static RanGeneratorPvt* createSharedNull()
{
    auto mutex = SireBase::detail::get_shared_null_mutex();
    
    tbb::spin_mutex::scoped_lock lock(*mutex);

    RanGeneratorPvt *gen = new RanGeneratorPvt();
    
    int seed = gen->mersenne_generator.randInt();
    
    //QUuid *very annoyingly* calls qsrand using the current time.
    // THIS IS REALLY ANNOYING WHEN USING QUuid IN AN MPI PROGRAM!!!!
    
    // However, it only calls qsrand on the first QUuid - I'll thus
    // call qsrand after the first...
    QUuid::createUuid();
             
    qsrand(seed);

    //lets dispose of another QUuid while we're at it
    QUuid::createUuid();
    
    //what is even more annoying, is that qsrand is a *per-thread*
    //seed, so this has only fixed this thread. Other threads
    //will use a seed of 1, so will all create the same sequence
    //of QUuids - the only way to fix this is for all new threads
    //to seed qsrand
    
    return gen;
}

static boost::shared_ptr<RanGeneratorPvt> shared_null( ::createSharedNull() );

/** Create a randomly seeded generator
    (actually a copy of the global, random generator) */
RanGenerator::RanGenerator() : d(shared_null)
{}

/** Create a generator seeded with 'seed' */
RanGenerator::RanGenerator(quint32 seed)
             : d( new RanGeneratorPvt(seed) )
{}

/** Create a generator seeded with 'seed' */
RanGenerator::RanGenerator(const QVector<quint32> &seed)
             : d( new RanGeneratorPvt(seed) )
{}

/** Copy constructor - this takes an explicitly shared
    copy of 'other' (this is to prevent repeat random numbers
    from being generated by implicit copies!) */
RanGenerator::RanGenerator(const RanGenerator &other)
             : d( other.d )
{}

/** Destructor */
RanGenerator::~RanGenerator()
{}

/** Detach from shared storage */
void RanGenerator::detach()
{
    if (not d.unique())
    {
        d.reset( new RanGeneratorPvt(*d) );
    }
}

/** Copy assignment */
RanGenerator& RanGenerator::operator=(const RanGenerator &other)
{
    d = other.d;
    return *this;
}

/** Comparison operator - two generators are equal if they use the
    same underlying generator */
bool RanGenerator::operator==(const RanGenerator &other) const
{
    return d == other.d;
}

/** Comparison operator - two generators are equal if they use the
    same underlying generator */
bool RanGenerator::operator!=(const RanGenerator &other) const
{
    return d != other.d;
}

RanGeneratorPvt& RanGenerator::nonconst_d() const
{
    return const_cast<RanGeneratorPvt&>(*d);
}

/** See the generator with a new, random seed - this will detach
    this explicitly shared copy of the generator */
void RanGenerator::seed()
{
    if (d.unique())
        d->seed();
    else
        d.reset( new RanGeneratorPvt() );
}

/** Seed the generator with 's'  - this will detach
    this explicitly shared copy of the generator */
void RanGenerator::seed(quint32 s)
{
    if (d.unique())
        d->seed(s);
    else
        d.reset( new RanGeneratorPvt(s) );
}

/** Seed the generator with 'seed' - this will detach
    this explicitly shared copy of the generator */
void RanGenerator::seed(const QVector<quint32> &s)
{
    if (d.unique())
        d->seed(s);
    else
        d.reset( new RanGeneratorPvt(s) );
}

/** Seed the generator with another generator - this
    really just copies the generator as they are
    all explicit copies of one another! */
void RanGenerator::seed(const RanGenerator &other)
{
    d = other.d;
}

/** Call this function to seed the qrand generator for this thread */
namespace SireMaths
{
    void seed_qrand()
    {
        RanGenerator rand;
        qsrand( rand.randInt() );
    }
}

/** Return a random real number on [0,1] */
double RanGenerator::rand() const
{
    QMutexLocker lkr( &(nonconst_d().mutex) );
    return nonconst_d().mersenne_generator.rand();
}

/** Return a random real number on [0,maxval] */
double RanGenerator::rand(double maxval) const
{
    return maxval * rand();
}

/** Return a random real number on [minval,maxval] */
double RanGenerator::rand(double minval, double maxval) const
{
    if (minval > maxval)
        return rand(maxval, minval);
    else
        return minval + rand() * (maxval-minval);
}

/** Take hold of the generator lock. Only you can now generate
    random numbers while this lock is held */
void RanGenerator::lock() const
{
    nonconst_d().mutex.lock();
}

/** Release the generator lock */
void RanGenerator::unlock() const
{
    nonconst_d().mutex.unlock();
}

/** Return a random real number on [0,1]. Should only be called while
    you hold the generator lock */
double RanGenerator::locked_rand() const
{
    return nonconst_d().mersenne_generator.rand();
}

/** Return a random real number on [0,maxval]. Should only be called while
    you hold the generator lock */
double RanGenerator::locked_rand(double maxval) const
{
    return maxval * locked_rand();
}

/** Return a random real number on [minval,maxval]. Should only be called while
    you hold the generator lock */
double RanGenerator::locked_rand(double minval, double maxval) const
{
    if (minval > maxval)
        return locked_rand(maxval, minval);
    else
        return minval + locked_rand() * (maxval-minval);
}

/** Fill the passed array of doubles with random numbers. This replaces each
    value in the array with a random number on [0,1] */
void RanGenerator::nrand(QVector<double> &result) const
{
    const int n = result.count();

    if (n > 0)
    {
        QMutexLocker lkr( &(nonconst_d().mutex) );
    
        double *d = result.data();
    
        for (int i=0; i<n; ++i)
        {
            d[i] = nonconst_d().mersenne_generator.rand();
        }
    }
}

/** Fill the passed array of doubles with random numbers. This replaces each
    value in the array with a random number on [0,maxval] */
void RanGenerator::nrand(QVector<double> &result, double maxval) const
{
    const int n = result.count();

    if (n > 0)
    {
        double *d = result.data();

        if (maxval == 0)
        {
            for (int i=0; i<n; ++i)
            {
                d[i] = 0;
            }
        }
        else
        {
            QMutexLocker lkr( &(nonconst_d().mutex) );

            for (int i=0; i<n; ++i)
            {
                d[i] = maxval * nonconst_d().mersenne_generator.rand();
            }
        }
    }
}

/** Fill the passed array of doubles with random numbers. This replaces each
    value in the array with a random number on [minval,maxval] */
void RanGenerator::nrand(QVector<double> &result, double minval, double maxval) const
{
    if (minval < maxval)
    {
        nrand(result, maxval, minval);
        return;
    }

    const int n = result.count();

    if (n > 0)
    {
        double *d = result.data();
        
        const double diff = maxval - minval;
        
        if (diff == 0)
        {
            for (int i=0; i<n; ++i)
            {
                d[i] = 0;
            }
        }
        else
        {
            QMutexLocker lkr( &(nonconst_d().mutex) );

            for (int i=0; i<n; ++i)
            {
                d[i] = minval + (diff * nonconst_d().mersenne_generator.rand());
            }
        }
    }
}

/** Return an array of 'n' random numbers on [0,1] */
QVector<double> RanGenerator::nrand(int n) const
{
    if (n > 0)
    {
        QVector<double> result(n);
        this->nrand(result);
        return result;
    }
    else
        return QVector<double>();
}

/** Return an array of 'n' random numbers on [0,maxval] */
QVector<double> RanGenerator::nrand(int n, double maxval) const
{
    if (n > 0)
    {
        QVector<double> result(n);
        this->nrand(result, maxval);
        return result;
    }
    else
        return QVector<double>();
}

/** Return an array of 'n' random numbers on [minval,maxval] */
QVector<double> RanGenerator::nrand(int n, double minval, double maxval) const
{
    if (n > 0)
    {
        QVector<double> result(n);
        this->nrand(result, minval, maxval);
        return result;
    }
    else
        return QVector<double>();
}

/** Return a high-precision random real number on [0,1) */
double RanGenerator::rand53() const
{
    QMutexLocker lkr( &(nonconst_d().mutex) );
    return nonconst_d().mersenne_generator.rand53();
}

/** Return a high-precision random real number on [0,1) */
double RanGenerator::rand53(double maxval) const
{
    return maxval * rand53();
}

/** Return a high-precision random real number on [minval,maxval) */
double RanGenerator::rand53(double minval, double maxval) const
{
    return minval + rand53()*(maxval-minval);
}

/** Return a random number from the normal distribution
    with supplied mean and variance. */
double RanGenerator::randNorm(double mean, double variance) const
{
    QMutexLocker lkr( &(nonconst_d().mutex) );
    return nonconst_d().mersenne_generator.randNorm(mean, variance);
}

/** Return a random number generated from the normal distribution 
    with mean 0 and standard deviation 1 */
double RanGenerator::randNorm() const
{
    QMutexLocker lkr( &(nonconst_d().mutex) );
    return nonconst_d().mersenne_generator.randNorm(1,0);
}

/** Return a random number from the normal distribution
    with supplied mean and variance. You must hold the generator
    lock when calling this function */
double RanGenerator::locked_randNorm(double mean, double variance) const
{
    return nonconst_d().mersenne_generator.randNorm(mean, variance);
}

/** Return a random number generated from the normal distribution 
    with mean 0 and standard deviation 1. You must hold the generator
    lock when calling this function */
double RanGenerator::locked_randNorm() const
{
    return nonconst_d().mersenne_generator.randNorm(1,0);
}

/** Fill the passed array with random numbers drawn from the normal
    distribution with supplied mean and variance */
void RanGenerator::nrandNorm(QVector<double> &result, double mean, double variance) const
{
    const int n = result.count();
    
    if (n > 0)
    {
        double *d = result.data();
        
        QMutexLocker lkr( &(nonconst_d().mutex) );
        
        for (int i=0; i<n; ++i)
        {
            d[i] = nonconst_d().mersenne_generator.randNorm(mean, variance);
        }
    }
}

/** Return an array of 'N' random numbers drawn from the normal distribution with
    supplied mean and variance */
QVector<double> RanGenerator::nrandNorm(int n, double mean, double variance) const
{
    if (n > 0)
    {
        QVector<double> result(n);
        this->nrandNorm(result, mean, variance);
        return result;
    }
    else
        return QVector<double>();
}

/** Return a random vector on the unit sphere. You must hold the generator
    lock when calling this function */
Vector RanGenerator::locked_vectorOnSphere() const
{
    while (true)
    {
        Vector v;
        
        v.setX( 1.0 - 2.0 * nonconst_d().mersenne_generator.rand() );
        v.setY( 1.0 - 2.0 * nonconst_d().mersenne_generator.rand() );
        v.setZ( 1.0 - 2.0 * nonconst_d().mersenne_generator.rand() );
        
        const double lgth2 = v.length2();
        
        if (lgth2 < 1)
        {
            v /= std::sqrt(lgth2);
            return v;
        }
    }
}

/** Return a random vector on the unit sphere */
Vector RanGenerator::vectorOnSphere() const
{
    QMutexLocker lkr( &(nonconst_d().mutex) );
    return locked_vectorOnSphere();
}

/** Fill the passed array with random vectors on a unit sphere */
void RanGenerator::nvectorOnSphere(QVector<Vector> &result) const
{
    const int n = result.count();
    
    if (n > 0)
    {
        QMutexLocker lkr( &(nonconst_d().mutex) );
        
        int i=0;
        
        while (i < n)
        {
            Vector &v = result[i];
            
            v.setX( 1.0 - 2.0 * nonconst_d().mersenne_generator.rand() );
            v.setY( 1.0 - 2.0 * nonconst_d().mersenne_generator.rand() );
            v.setZ( 1.0 - 2.0 * nonconst_d().mersenne_generator.rand() );
            
            const double lgth2 = v.length2();
            
            if (lgth2 < 1)
            {
                v /= std::sqrt(lgth2);
                i += 1;
            }
        }
    }
}

/** Return an array of 'n' random vectors on a unit sphere */
QVector<Vector> RanGenerator::nvectorOnSphere(int n) const
{
    if (n > 0)
    {
        QVector<Vector> result(n);
        this->nvectorOnSphere(result);
        return result;
    }
    else
        return QVector<Vector>();
}

/** Return a random vector on the sphere with radius 'radius' */
Vector RanGenerator::vectorOnSphere(double radius) const
{
    return radius * this->vectorOnSphere();
}

/** Return a random vector on the sphere with radius 'radius'.
    You must hold the generator lock when calling this function */
Vector RanGenerator::locked_vectorOnSphere(double radius) const
{
    return radius * this->locked_vectorOnSphere();
}

/** Fill the passed array with random vectors on a sphere with radius 'radius' */
void RanGenerator::nvectorOnSphere(QVector<Vector> &result, double radius) const
{
    const int n = result.count();
    
    if (n > 0)
    {
        QMutexLocker lkr( &(nonconst_d().mutex) );
        
        int i=0;
        
        while (i < n)
        {
            Vector &v = result[i];
            
            v.setX( 1.0 - 2.0 * nonconst_d().mersenne_generator.rand() );
            v.setY( 1.0 - 2.0 * nonconst_d().mersenne_generator.rand() );
            v.setZ( 1.0 - 2.0 * nonconst_d().mersenne_generator.rand() );
            
            const double lgth2 = v.length2();
            
            if (lgth2 < 1)
            {
                v *= (radius * std::sqrt(lgth2));
                i += 1;
            }
        }
    }
}

/** Return an array of 'n' random vectors on a sphere of radius 'radius' */
QVector<Vector> RanGenerator::nvectorOnSphere(int n, double radius) const
{
    if (n > 0)
    {
        QVector<Vector> result(n);
        this->nvectorOnSphere(result, radius);
        return result;
    }
    else
        return QVector<Vector>();
}

/** Return a random 32bit unsigned integer in [0,2^32 - 1] */
quint32 RanGenerator::randInt() const
{
    QMutexLocker lkr( &(nonconst_d().mutex) );
    return nonconst_d().mersenne_generator.randInt();
}

/** Return a random true or false value */
bool RanGenerator::randBool() const
{
    return this->randInt() & 0x0001;
}

/** Return a random 32bit unsigned integer in [0,maxval] */
quint32 RanGenerator::randInt(quint32 maxval) const
{
    QMutexLocker lkr( &(nonconst_d().mutex) );
    return nonconst_d().mersenne_generator.randInt(maxval);
}

/** Return a random 32bit integer in [minval,maxval] */
qint32 RanGenerator::randInt(qint32 minval, qint32 maxval) const
{
    if (maxval == minval)
        return maxval;
    else if (maxval < minval)
        qSwap(minval,maxval);

    return minval + randInt(maxval-minval);
}

static quint64 randInt64(MTRand &mersenne_generator)
{
    quint64 ran0 = mersenne_generator.randInt();
    quint64 ran1 = mersenne_generator.randInt();

    return (ran0 << 32) | ran1;
}

/** Return a random 64bit unsigned integer on [0,2^64 - 1] */
quint64 RanGenerator::randInt64() const
{
    QMutexLocker lkr( &(nonconst_d().mutex) );
    return ::randInt64(nonconst_d().mersenne_generator);
}

/** Return a random 64bit unsigned integer on [0,maxval] */
quint64 RanGenerator::randInt64(quint64 maxval) const
{
    QMutexLocker lkr( &(nonconst_d().mutex) );

    if (maxval <= std::numeric_limits<quint32>::max())
        //maxval can fit into a 32bit int - there is no
        //point using a 64bit generator!
        return nonconst_d().mersenne_generator.randInt( quint32(maxval) );

    //use same algorithm in MersenneTwister.h
    quint64 used = maxval;

    used |= used >> 1;
    used |= used >> 2;
    used |= used >> 4;
    used |= used >> 8;
    used |= used >> 16;
    used |= used >> 32;

    quint64 i;

    do
    {
        i = ::randInt64(nonconst_d().mersenne_generator) & used;
    } while ( i > maxval );

    return i;
}

/** Return a random 64bit integer on [minval,maxval] */
qint64 RanGenerator::randInt64(qint64 minval, qint64 maxval) const
{
    if (maxval == minval)
        return maxval;
    else if (maxval < minval)
        qSwap(minval,maxval);

    return minval + randInt(maxval-minval);
}

/** Return the current state of the random number generator.
    Use this if you truly wish to get reproducible sequences
    of random numbers */
QVector<quint32> RanGenerator::getState() const
{
    return nonconst_d().getState();
}

/** Load the state into this generator - the state must have
    been produced by the getState() function above.

    This will detach this copy from shared storage.

    \throw SireError::incompatible_error
*/
void RanGenerator::setState(const QVector<quint32> &state)
{
    if (d.unique())
    {
        d->loadState(state);
    }
    else
    {
        d.reset( new RanGeneratorPvt() );
        d->loadState(state);
    }
}

Q_GLOBAL_STATIC( RanGenerator, globalGenerator );

/** Return a reference to the global random number generator 
    (shared between all threads) */
const RanGenerator& RanGenerator::global()
{
    return *(globalGenerator());
}

/** Seed the global random number generator */
void RanGenerator::seedGlobal()
{
    globalGenerator()->seed();
}

/** Seed the global random number generator */
void RanGenerator::seedGlobal(quint32 seed)
{
    globalGenerator()->seed(seed);
}

/** Seed the global random number generator */
void RanGenerator::seedGlobal(const QVector<quint32> &seed)
{
    globalGenerator()->seed(seed);
}

/** Seed the global random number generator */
void RanGenerator::seedGlobal(const RanGenerator &rangen)
{
    globalGenerator()->d.reset( new RanGeneratorPvt(*(rangen.d)) );
}

const char* RanGenerator::typeName()
{
    return QMetaType::typeName( qMetaTypeId<RanGenerator>() );
}
