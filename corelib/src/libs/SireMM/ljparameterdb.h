/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2008  Christopher Woods
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

#ifndef SIREMM_LJPARAMETERDB_H
#define SIREMM_LJPARAMETERDB_H

#include "ljparameter.h"

#include "SireBase/array2d.hpp"

#include <QVector>
#include <QHash>
#include <QReadWriteLock>

#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>

SIRE_BEGIN_HEADER

namespace SireStream
{
class SharedDataStream;
}

namespace SireMM
{

class LJPair;
typedef SireBase::Array2D<LJPair> LJPairMatrix;

namespace detail
{
    class LJDBIOLockData : public boost::noncopyable
    {
    public:
        LJDBIOLockData();
        ~LJDBIOLockData();
    };
}

/** Class used to help the serialisation of the LJParameterDB */
class SIREMM_EXPORT LJDBIOLock
{
public:
    LJDBIOLock();
    LJDBIOLock(const LJDBIOLock &other);
    ~LJDBIOLock();
    
private:
    friend class LJParameterDB;
    LJDBIOLock( const boost::shared_ptr<detail::LJDBIOLockData> &p );

    boost::shared_ptr<detail::LJDBIOLockData> d;
};

/** This static singleton class holds a complete database of 
    all of the LJ parameters used during the simulation. This
    provides rapid access to the combined parameters during the
    pair loop, removing the need to combine the LJ parameters
    together at simulation time (as it can all be done during
    setup time). A singleton class is used to save memory - 
    it is not necessary for each LJ potential to keep a copy
    of its own database, as in the normal case there will be
    several common LJ parameters.
  
    This class is both reentrant and thread-safe
        
    @author Christopher Woods
*/
class SIREMM_EXPORT LJParameterDB
{
public:
    enum CombiningRules { ARITHMETIC, GEOMETRIC };

    static LJDBIOLock saveParameters(SireStream::SharedDataStream &sds);
    static LJDBIOLock loadParameters(SireStream::SharedDataStream &sds);

    static LJPairMatrix getLJPairs(CombiningRules type);
    static quint32 addLJParameter(const LJParameter &ljparam);
    static LJParameter getLJParameter(quint32 id);

    static CombiningRules interpret(const QString &rule);
    static const QString& toString(CombiningRules rule);

    static void lock();
    
    static quint32 _locked_addLJParameter(const LJParameter &ljparam);
    static LJParameter _locked_getLJParameter(quint32 id);
    
    static void unlock();

private:
    friend class ::SireMM::detail::LJDBIOLockData;
    static void finishedIO();

    class SIREMM_EXPORT LJParameterDBData
    {
    public:
        LJParameterDBData();
        ~LJParameterDBData();
    
        const LJPairMatrix& getLJPairs(CombiningRules type);
        quint32 addLJParameter(const LJParameter &ljparam);
        LJParameter getLJParameter(quint32 id);
    
        void lock();
        quint32 _locked_addLJParameter(const LJParameter &ljparam);
        LJParameter _locked_getLJParameter(quint32 id) const;
        void unlock();
                
    private:
        LJPairMatrix combineArithmetic() const;
        LJPairMatrix combineGeometric() const;
    
        friend class LJParameterDB;
    
        /** Read-write lock used to control access to shared resources */
        QReadWriteLock db_lock;
    
        /** Any requested LJPair arrays, indexed by CombiningRule type */
        QHash<int, LJPairMatrix> ljpair_arrays;
    
        /** All of the LJ parameters indexed by their ID */
        QVector<LJParameter> ljparams_by_idx;
    
        /** Index allowing reverse lookup of a LJParameter's ID number */
        QHash<LJParameter,quint32> ljparams_by_value;
        
        /** Whether or not we are in the process of saving or loading the data */
        int dbio_count;
    };
    
    static LJParameterDBData ljdb;
};

} // end of namespace SireMM

#include "ljpair.h"

namespace SireMM
{

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Return the matrix of LJ pairs for the given combining rules */
inline LJPairMatrix LJParameterDB::getLJPairs(CombiningRules type)
{
    return ljdb.getLJPairs(type);
}

/** Add a new LJParameter to the database, returning the ID of the 
    parameter. This has to lock and unlock the database, so it may
    be slow if you are adding large numbers of parameters */
inline quint32 LJParameterDB::addLJParameter(const LJParameter &ljparam)
{
    return ljdb.addLJParameter(ljparam);
}

/** Get the ID number of a LJ parameter

    \throw SireFF::missing_parameter
*/
inline LJParameter LJParameterDB::getLJParameter(quint32 id)
{
    return ljdb.getLJParameter(id);
}

/** Lock the database - use this if you are going to add lots of 
    parameters via the '_locked_addLJParameter' function */
inline void LJParameterDB::lock()
{
    ljdb.lock();
}

/** Add a new LJParameter to the database, returning the ID of the 
    parameter. You can only call this function if you have manually
    locked the database via the lock() function, and you must
    unlock the database via the unlock() function once you have
    finished adding all of the parameters */
inline quint32 LJParameterDB::_locked_addLJParameter(const LJParameter &ljparam)
{
    return ljdb._locked_addLJParameter(ljparam);
}

/** Get the LJ parameter corresponding to the passed ID number
    - you can only call this function
    if you have manually locked the database via the lock() function

    \throw SireFF::missing_parameter
*/
inline LJParameter LJParameterDB::_locked_getLJParameter(quint32 id)
{
    return ljdb._locked_getLJParameter(id);
}

/** Unlock the database - ensure that you do this after you have finished
    adding parameter via the '_locked_addLJParameter' function, or else
    weird things may happen! */
inline void LJParameterDB::unlock()
{
    ljdb.unlock();
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

SIRE_END_HEADER

#endif
