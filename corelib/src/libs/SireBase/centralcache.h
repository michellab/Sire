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

#ifndef SIREBASE_CENTRALCACHE_H
#define SIREBASE_CENTRALCACHE_H

#include "sireglobal.h"

#include <QVariant>
#include <QHash>
#include <QDateTime>
#include <QMutex>

#include <memory>

#include <boost/noncopyable.hpp>

SIRE_BEGIN_HEADER

namespace SireBase
{

class CentralCache;

/** This is an item of cached data. Note that you
 *  must hold the QMutex in this data whenever you
 *  try to access or use this data
 */
class SIREBASE_EXPORT CachedData : public boost::noncopyable
{
public:
    CachedData();
    ~CachedData();

    bool isEmpty() const;

    QMutex* mutex();

    const QVariant& data() const;

    int nBytes() const;

    void setData(const QVariant &data, int num_bytes);

private:
    QVariant d;
    QMutex m;
    int num_bytes;
};


/** This provides a centralised cache, which evicts items
 *  in reverse order of last accessed
 */
class SIREBASE_EXPORT CentralCache : public boost::noncopyable
{
public:
    CentralCache();
    ~CentralCache();

    static std::shared_ptr<CachedData> get(const QString &key);

    static void setMaxCacheSize(qint64 size_in_bytes);
    static qint64 getMaxCacheSize();

    static qint64 getCacheSize();

    static void empty();

private:
    void clean();

    static CentralCache global_cache;

    QMutex mutex;
    QHash<QString, std::shared_ptr<CachedData>> cache;
    QHash<QString,qint64> last_access_times;

    qint64 last_clean_time;

    qint64 current_cache_size;
    qint64 max_cache_size;
};

}

SIRE_END_HEADER

#endif
