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

#include "centralcache.h"

#include <QDebug>

using namespace SireBase;

CachedData::CachedData() : boost::noncopyable(), num_bytes(0)
{}

CachedData::~CachedData()
{}

bool CachedData::isEmpty() const
{
    return d.isNull();
}

QMutex* CachedData::mutex()
{
    return &m;
}

const QVariant& CachedData::data() const
{
    return d;
}

int CachedData::nBytes() const
{
    return num_bytes;
}

void CachedData::setData(const QVariant &data, int n)
{
    d = data;
    num_bytes = n;

    if (num_bytes < 0)
        num_bytes = 0;
}

CentralCache CentralCache::global_cache;

// Construct, setting the initial cache size to 512 MB
CentralCache::CentralCache()
             : last_clean_time(-1), max_cache_size(512*1024*1024)
{
    // reserve space to cut down on mallocs
    const int reserve_size = 4097;

    this->cache.reserve(reserve_size);
    this->last_access_times.reserve(reserve_size);
}

CentralCache::~CentralCache()
{}

void CentralCache::setMaxCacheSize(qint64 v)
{
    if (v < 1000)
        return;

    QMutexLocker lkr(&(global_cache.mutex));
    global_cache.max_cache_size = v;
}

qint64 CentralCache::getMaxCacheSize()
{
    QMutexLocker lkr(&(global_cache.mutex));
    return global_cache.max_cache_size;
}

std::shared_ptr<CachedData> CentralCache::get(const QString &key)
{
    auto now = QDateTime::currentSecsSinceEpoch();

    QMutexLocker lkr(&(global_cache.mutex));

    if (not global_cache.cache.contains(key))
    {
        global_cache.cache.insert(key, std::shared_ptr<CachedData>(new CachedData()));
    }

    global_cache.last_access_times[key] = now;

    auto cached_data = global_cache.cache[key];

    if (global_cache.last_clean_time == -1)
    {
        global_cache.last_clean_time = now;
    }
    else if (now - global_cache.last_clean_time > 30)
    {
        global_cache.clean();
        global_cache.last_clean_time = QDateTime::currentSecsSinceEpoch();
    }

    return cached_data;
}

void CentralCache::empty()
{
    QMutexLocker lkr(&(global_cache.mutex));
    global_cache.cache.clear();
    global_cache.last_access_times.clear();
}

void CentralCache::clean()
{
    // calculate the total size of the cache
    qint64 total_size = 0;

    for (const auto &key : this->cache.keys())
    {
        auto cached_data = this->cache.value(key);

        if (cached_data->mutex()->tryLock())
        {
            total_size += cached_data->nBytes();
            cached_data->mutex()->unlock();
        }
        else
        {
            // this item is still being used, so update the
            // last access time
            this->last_access_times[key] = QDateTime::currentSecsSinceEpoch();
        }
    }

    while (total_size > this->max_cache_size)
    {
        // we need to remove items in last-accessed order
        qint64 oldest_time = -1;
        QString oldest_key;
        qint64 oldest_size = 0;

        for (const auto &key : this->last_access_times.keys())
        {
            const auto &t = this->last_access_times[key];

            if (oldest_time == -1 or t < oldest_time)
            {
                auto cached_data = this->cache.value(key);

                if (cached_data->mutex()->tryLock())
                {
                    if (cached_data->nBytes() > 0)
                    {
                        oldest_time = t;
                        oldest_key = key;
                        oldest_size = cached_data->nBytes();
                    }

                    cached_data->mutex()->unlock();
                }
                else
                {
                    cached_data->mutex()->unlock();
                }
            }
        }

        if (oldest_time == -1 or oldest_key.isEmpty())
            break;

        this->last_access_times.remove(oldest_key);
        this->cache.remove(oldest_key);
        total_size -= oldest_size;
    }
}

qint64 CentralCache::getCacheSize()
{
    QMutexLocker lkr(&(global_cache.mutex));
    return global_cache.current_cache_size;
}
