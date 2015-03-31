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

#include <QStringList>

#include "ljparameterdb.h"

#include "SireFF/errors.h"
#include "SireError/errors.h"

using namespace SireMM;
using namespace SireBase;

static QHash<QString,LJParameterDB::CombiningRules> rule_types;

LJParameterDB::LJParameterDBData LJParameterDB::ljdb;

/** Convert the passed string into the CombiningRules ID */
LJParameterDB::CombiningRules LJParameterDB::interpret(const QString &rule)
{
    if (rule_types.isEmpty())
    {
        LJParameterDB::lock();
        
        rule_types.insert("arithmetic", ARITHMETIC);
        rule_types.insert("geometric", GEOMETRIC);
        
        LJParameterDB::unlock();
    }
    
    QHash<QString,CombiningRules>::const_iterator it = rule_types.constFind(rule);
    
    if (it == rule_types.constEnd())
        throw SireError::invalid_arg( QObject::tr(
            "There is no combining rule available that matches the ID \"%1\". "
            "Available rules are [ %2 ].")
                .arg(rule, QStringList(rule_types.keys()).join(", ")), CODELOC );
                
    return *it;
}

/** Convert the passed rule into the string representing that rule */
const QString& LJParameterDB::toString(LJParameterDB::CombiningRules rule)
{
    switch (rule)
    {
        case ARITHMETIC:
            return rule_types.constFind("arithmetic").key();
        case GEOMETRIC:
            return rule_types.constFind("geometric").key();
    }
    
    return rule_types.constFind("arithmetic").key();
}

using namespace SireMM::detail;

LJDBIOLockData::LJDBIOLockData() : boost::noncopyable()
{}

LJDBIOLockData::~LJDBIOLockData()
{
    LJParameterDB::finishedIO();
}

LJDBIOLock::LJDBIOLock()
{}

LJDBIOLock::LJDBIOLock(const boost::shared_ptr<LJDBIOLockData> &ptr) : d(ptr)
{}

LJDBIOLock::LJDBIOLock(const LJDBIOLock &other) : d(other.d)
{}

LJDBIOLock::~LJDBIOLock()
{}

/** Save the LJParameter database to the datastream */
LJDBIOLock LJParameterDB::saveParameters(SireStream::SharedDataStream &sds)
{
    ljdb.dbio_count += 1;
    
    if (ljdb.dbio_count == 1)
    {
        QReadLocker lkr( &(ljdb.db_lock) );
        sds << ljdb.ljparams_by_idx;
    }
    
    return LJDBIOLock( boost::shared_ptr<detail::LJDBIOLockData>( new detail::LJDBIOLockData() ) );
}

/** Load the LJParameter database from the datastream */
LJDBIOLock LJParameterDB::loadParameters(SireStream::SharedDataStream &sds)
{
    ljdb.dbio_count += 1;
    
    if (ljdb.dbio_count == 1)
    {
        QVector<LJParameter> ljparams;
        sds >> ljparams;
        
        ljdb.lock();
        foreach (const LJParameter &ljparam, ljparams)
        {
            ljdb._locked_addLJParameter(ljparam);
        }
        ljdb.unlock();
    }
    
    return LJDBIOLock( boost::shared_ptr<detail::LJDBIOLockData>( new detail::LJDBIOLockData() ) );
}

void LJParameterDB::finishedIO()
{
    ljdb.dbio_count -= 1;
    
    if (ljdb.dbio_count < 0){ ljdb.dbio_count = 0; }
}

/** Constructor */
LJParameterDB::LJParameterDBData::LJParameterDBData()
{
    //add the null LJ parameter - this has ID == 0
    ljparams_by_idx.append( LJParameter::dummy() );
    ljparams_by_value.insert( LJParameter::dummy(), 0 );
    dbio_count = 0;
}

/** Destructor */
LJParameterDB::LJParameterDBData::~LJParameterDBData()
{}

/** Return a square matrix of all of the combined LJ parameters, combined
    using geometric combining rules */
LJPairMatrix LJParameterDB::LJParameterDBData::combineGeometric() const
{
    int nlj = ljparams_by_idx.count();

    LJPairMatrix ljpairs( nlj, nlj );
    
    LJPair *ljpairs_array = ljpairs.data();
    
    const LJParameter *ljparams_array = ljparams_by_idx.constData();
    
    for (int i=0; i<nlj; ++i)
    {
        const LJParameter &ljparam0 = ljparams_array[i];
        LJPair *ljpairs_row = ljpairs_array + ljpairs.map(i,0);
    
        for (int j=0; j<nlj; ++j)
        {
            ljpairs_row[j] = LJPair::geometric(ljparam0, ljparams_array[j]);
        }
    }
    
    return ljpairs;
}

/** Return a square matrix of all of the combined LJ parameters, combined
    using arithmetic combining rules */
LJPairMatrix LJParameterDB::LJParameterDBData::combineArithmetic() const
{
    int nlj = ljparams_by_idx.count();

    LJPairMatrix ljpairs( nlj, nlj );
    
    LJPair *ljpairs_array = ljpairs.data();
    
    const LJParameter *ljparams_array = ljparams_by_idx.constData();
    
    for (int i=0; i<nlj; ++i)
    {
        const LJParameter &ljparam0 = ljparams_array[i];
        LJPair *ljpairs_row = ljpairs_array + ljpairs.map(i,0);
    
        for (int j=0; j<nlj; ++j)
        {
            ljpairs_row[j] = LJPair::arithmetic(ljparam0, ljparams_array[j]);
        }
    }
    
    return ljpairs;
}

/** Return the matrix of LJ pairs for the given combining rules */
const LJPairMatrix& LJParameterDB::LJParameterDBData::getLJPairs(CombiningRules type)
{
    //use a scope so that variables below are private - this allows
    //me to use a QReadLocker to safely lock the database
    {
        QReadLocker lkr(&db_lock);
        
        QHash<int, LJPairMatrix>::const_iterator it = ljpair_arrays.constFind(type);
                                    
        if (it != ljpair_arrays.constEnd())
            return it.value();
    
        //read lock is released at exit of this scope
    }
    
    //the array wasn't found - it needs to be constructed
    QWriteLocker lkr(&db_lock);
    
    switch (type)
    {
        case ARITHMETIC:
            ljpair_arrays.insert( type, this->combineArithmetic() );
            break;
        
        case GEOMETRIC:
            ljpair_arrays.insert( type, this->combineGeometric() );
            break;
    }
    
    return *(ljpair_arrays.constFind(type));
}

/** Lock the database - use this if you are going to add lots of 
    parameters via the '_locked_addLJParameter' function */
void LJParameterDB::LJParameterDBData::lock()
{
    db_lock.lockForWrite();
}

/** Add a new LJParameter to the database, returning the ID of the 
    parameter. You can only call this function if you have manually
    locked the database via the lock() function, and you must
    unlock the database via the unlock() function once you have
    finished adding all of the parameters */
quint32 LJParameterDB::LJParameterDBData::_locked_addLJParameter(
                                                    const LJParameter &ljparam)
{
    //does the parameter exist already in the database?
    QHash<LJParameter,quint32>::const_iterator 
                                    it = ljparams_by_value.constFind(ljparam);
                                    
    if (it != ljparams_by_value.constEnd())
    {
        return it.value();
    }
    
    //the parameter needs to be added
    quint32 idx = ljparams_by_idx.count();
    
    ljparams_by_idx.append(ljparam);
    ljparams_by_value.insert(ljparam, idx);
    
    //clear the existing LJ pair matricies, as these are now invalid
    ljpair_arrays.clear();
    
    return idx;
}

/** Unlock the database - ensure that you do this after you have finished
    adding parameter via the '_locked_addLJParameter' function, or else
    weird things may happen! */
void LJParameterDB::LJParameterDBData::unlock()
{
    db_lock.unlock();
}

/** Add a new LJParameter to the database, returning the ID of the 
    parameter. This has to lock and unlock the database, so it may
    be slow if you are adding large numbers of parameters */
quint32 LJParameterDB::LJParameterDBData::addLJParameter(const LJParameter &ljparam)
{
    QWriteLocker lkr(&db_lock);
    return this->_locked_addLJParameter(ljparam);
}

/** Get the ID number of a LJ parameter

    \throw SireFF::missing_parameter
*/
LJParameter LJParameterDB::LJParameterDBData::_locked_getLJParameter(quint32 id) const
{
    //does the parameter exist already in the database?
    if (id > quint32(ljparams_by_idx.count()))
        throw SireFF::missing_parameter( QObject::tr(
            "Could not find the LJ parameter with ID %1 in the global "
            "LJ parameter database.").arg(id), CODELOC );
            
    return ljparams_by_idx.constData()[id];
}

/** Get the ID number of a LJ parameter

    \throw SireFF::missing_parameter
*/
LJParameter LJParameterDB::LJParameterDBData::getLJParameter(quint32 id)
{
    QReadLocker lkr(&db_lock);
    return this->_locked_getLJParameter(id);
}
