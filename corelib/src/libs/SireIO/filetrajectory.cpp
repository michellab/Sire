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

#include "filetrajectory.h"

#include <QMutex>
#include <QHash>

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireIO;
using namespace SireMol;
using namespace SireStream;

static const RegisterMetaType<FileTrajectory> r_traj;

SIREIO_EXPORT QDataStream& operator<<(QDataStream &ds, const FileTrajectory &file)
{
    writeHeader(ds, r_traj, 1);

    SharedDataStream sds(ds);

    QString filename;

    if (not file.parser.isNull())
        filename = file.parser.read().filename();

    sds << filename;

    return ds;
}

SIREIO_EXPORT QDataStream& operator>>(QDataStream &ds, FileTrajectory &file)
{
    VersionID v = readHeader(ds, r_traj);

    if (v == 1)
    {
        QString filename;
        SharedDataStream sds(ds);

        sds >> filename;

        if (filename.isEmpty())
        {
            qDebug() << "Cannot load the trajectory as the filename is NULL!";
            return ds;
        }

        try
        {
            file.parser = MoleculeParser::parse(filename);
        }
        catch(const SireError::exception &e)
        {
            qDebug() << "Cannot load the trajectory as something went wrong loading the file!";
            qDebug() << e.what() << e.error();
        }

        if (file.parser.read().isBroken())
        {
            qDebug() << "Cannot load the trajectory as the file was broken.";
            qDebug() << file.parser.read().errorReport();
            file.parser = MoleculeParserPtr();
        }
    }
    else
        throw version_error(v, "1", r_traj, CODELOC);

    return ds;
}

class TrajectoryCache
{
public:
    TrajectoryCache() : max_frames(100)
    {}

    ~TrajectoryCache()
    {}

    bool contains(int idx) const
    {
        return cache.contains(idx);
    }

    Frame get(int idx) const
    {
        return cache[idx];
    }

    void save(int idx, const Frame &frame)
    {
        if (order.isEmpty())
        {
            const int nbytes = frame.numBytes();

            if (nbytes == 0)
                return;

            // calculate an appropriate number of frames to cache,
            // assuming we don't want to consume more than ~512 MB
            max_frames = qMax(1, int(512*1024*1024 / nbytes) + 1);
        }

        while (order.count() >= max_frames-1)
        {
            if (order.isEmpty())
                break;

            cache.remove(order.takeFirst());
        }

        if (not cache.contains(idx))
        {
            order.append(idx);
        }

        cache[idx] = frame;
    }

    // you must hold this mutex when calling any of the
    // above functions!
    QMutex mutex;

private:
    QHash<int, Frame> cache;
    QList<int> order;

    int max_frames;
};

class TrajectoriesCache
{
public:
    static std::shared_ptr<TrajectoryCache> get(const void *ptr)
    {
        QMutexLocker lkr(&TrajectoriesCache::mutex);

        if (not TrajectoriesCache::cache.contains(ptr))
        {
            TrajectoriesCache::cache.insert(ptr, std::shared_ptr<TrajectoryCache>(new TrajectoryCache()));
        }

        return TrajectoriesCache::cache[ptr];
    }

    static void deleted(const void *ptr)
    {
        QMutexLocker lkr(&TrajectoriesCache::mutex);
        TrajectoriesCache::cache.remove(ptr);
    }

private:
    static QMutex mutex;
    static QHash<const void*, std::shared_ptr<TrajectoryCache> > cache;
};

QMutex TrajectoriesCache::mutex;
QHash<const void*, std::shared_ptr<TrajectoryCache> > TrajectoriesCache::cache;

FileTrajectory::FileTrajectory() : TrajectoryData()
{}

FileTrajectory::FileTrajectory(const MoleculeParser &p)
               : TrajectoryData(), parser(p)
{}

FileTrajectory::FileTrajectory(const FileTrajectory &other)
               : TrajectoryData(other), parser(other.parser)
{}

FileTrajectory::~FileTrajectory()
{
    TrajectoriesCache::deleted(this);
}

const char* FileTrajectory::what() const
{
    return FileTrajectory::typeName();
}

const char* FileTrajectory::typeName()
{
    return QMetaType::typeName(qMetaTypeId<FileTrajectory>());
}

FileTrajectory* FileTrajectory::clone() const
{
    return new FileTrajectory(*this);
}

bool FileTrajectory::_equals(const TrajectoryData &other) const
{
    const FileTrajectory *ptr = dynamic_cast<const FileTrajectory*>(&other);

    if (ptr)
    {
        return ptr->filenames() == this->filenames();
    }
    else
        return false;
}

int FileTrajectory::nFrames() const
{
    if (parser.isNull())
        return 0;
    else
        return parser.read().nFrames();
}

int FileTrajectory::nAtoms() const
{
    if (parser.isNull())
        return 0;
    else
        return parser.read().nAtoms();
}

QStringList FileTrajectory::filenames() const
{
    if (parser.isNull())
        return QStringList();
    else
    {
        QString filename = parser.read().filename();

        QStringList f;

        if (not filename.isEmpty())
            f.append(filename);

        return f;
    }
}

Frame FileTrajectory::getFrame(int i) const
{
    auto cache = TrajectoriesCache::get(this);

    QMutexLocker lkr(&(cache->mutex));

    if (cache->contains(i))
    {
        return cache->get(i);
    }

    i = SireID::Index(i).map(this->nFrames());

    auto frame = parser.read().getFrame(i);

    cache->save(i, frame);

    return frame;
}

bool FileTrajectory::isEditable() const
{
    return false;
}
