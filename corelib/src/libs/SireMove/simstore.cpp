/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2009  Christopher Woods
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

#include <QDir>
#include <QFileInfo>

#include "simstore.h"

#include "SireError/errors.h"
#include "SireStream/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>

using namespace SireMove;
using namespace SireSystem;
using namespace SireStream;

static const RegisterMetaType<SimStore> r_simstore(NO_ROOT);

/** Serialise to a binary datastream */
QDataStream SIREMOVE_EXPORT &operator<<(QDataStream &ds, const SimStore &simstore)
{
    writeHeader(ds, r_simstore, 5);
    
    SharedDataStream sds(ds);

    if (simstore.isPackedToDisk())
    {
        SimStore memory_packed = simstore;
        memory_packed.packToMemory();

        sds << quint32(SimStore::TO_DISK) << memory_packed.compressed_data;
    }
    else if (simstore.isPackedToMemory())
    {
        sds << quint32(SimStore::TO_MEMORY) << simstore.compressed_data;
    }
    else
    {
        sds << quint32(SimStore::NOT_PACKED)
            << simstore.sim_system << simstore.sim_moves
            << quint32(simstore.last_packing_state);
    }
    
    sds << simstore.packed_dir;
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOVE_EXPORT &operator>>(QDataStream &ds, SimStore &simstore)
{
    VersionID v = readHeader(ds, r_simstore);
   
    SimStore new_store;
   
    if (v == 5)
    {
        quint32 packed_state( SimStore::NOT_PACKED );
        
        SharedDataStream sds(ds);
        sds >> packed_state;
        
        switch (packed_state)
        {
            case SimStore::TO_DISK:
            {
                sds >> new_store.compressed_data
                    >> new_store.packed_dir;
                
                new_store.last_packing_state = SimStore::TO_DISK;
                
                new_store.packToDisk();
            }
            break;
            
            case SimStore::TO_MEMORY:
            {
                sds >> new_store.compressed_data
                    >> new_store.packed_dir;
                    
                new_store.last_packing_state = SimStore::TO_MEMORY;
            }
            break;
            
            case SimStore::NOT_PACKED:
            {
                sds >> new_store.sim_system >> new_store.sim_moves;
                
                quint32 last_state( SimStore::TO_MEMORY );
                
                sds >> last_state;
                
                switch (last_state)
                {
                    case SimStore::TO_MEMORY:
                        new_store.last_packing_state = SimStore::TO_MEMORY;
                        break;
                        
                    case SimStore::TO_DISK:
                        new_store.last_packing_state = SimStore::TO_DISK;
                        break;
                        
                    default:
                        new_store.last_packing_state = SimStore::TO_MEMORY;
                }
                
                sds >> new_store.packed_dir;
            }
            break;
            
            default:
                throw SireStream::corrupted_data( QObject::tr(
                    "There was an error reading the state of the SimStore - "
                    "a value of %1 is not valid!").arg(packed_state), CODELOC );
        }
    }
    else if (v == 4)
    {
        SharedDataStream sds(ds);
        
        QString packed_dir;
        
        sds >> new_store.sim_system >> new_store.sim_moves
            >> new_store.compressed_data >> packed_dir;

        if (not packed_dir.isEmpty())
            new_store.packToDisk(packed_dir);
    }
    else if (v == 3)
    {
        SharedDataStream sds(ds);
        
        sds >> new_store.sim_system >> new_store.sim_moves >> new_store.compressed_data;
        
        new_store.packed_file.reset();
        new_store.packed_filename = QString::null;
    }
    else if (v == 2)
    {
        SharedDataStream sds(ds);
        
        sds >> new_store.compressed_data;
        
        new_store.packed_file.reset();
        new_store.packed_filename = QString::null;
    }
    else if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> new_store.sim_system >> new_store.sim_moves;
        
        new_store.compressed_data = QByteArray();
        new_store.packed_file.reset();
        new_store.packed_filename = QString::null;
    }
    else
        throw version_error( v, "1-5", r_simstore, CODELOC );

    simstore = new_store;

    return ds;
}

/** Null constructor */
SimStore::SimStore() : last_packing_state( SimStore::TO_MEMORY )
{}

/** Construct from the passed system and moves - optionally specify whether
    or not to compress the data now */
SimStore::SimStore(const System &system, const Moves &moves, bool compress)
         : sim_system(system), sim_moves(moves),
           last_packing_state( SimStore::TO_MEMORY )
{
    if (compress)
        this->pack();
}

/** Copy constructor */
SimStore::SimStore(const SimStore &other)
         : sim_system(other.sim_system), sim_moves(other.sim_moves),
           compressed_data(other.compressed_data),
           packed_dir(other.packed_dir),
           packed_file(other.packed_file),
           packed_filename(other.packed_filename),
           last_packing_state(other.last_packing_state)
{}

/** Destructor */
SimStore::~SimStore()
{}

/** Copy assignment operator */
SimStore& SimStore::operator=(const SimStore &other)
{
    if (this == &other)
        return *this;
    
    sim_system = other.sim_system;
    sim_moves = other.sim_moves;
    compressed_data = other.compressed_data;
    packed_dir = other.packed_dir;
    packed_file = other.packed_file;
    packed_filename = other.packed_filename;
    last_packing_state = other.last_packing_state;
       
    return *this;
}

/** Comparison operator */
bool SimStore::operator==(const SimStore &other) const
{
    if (this == &other)
        return true;
        
    return sim_system == other.sim_system and 
           sim_moves == other.sim_moves and
           compressed_data == other.compressed_data and
           packed_file.get() == other.packed_file.get();
}

/** Comparison operator */
bool SimStore::operator!=(const SimStore &other) const
{
    return not this->operator==(other);
}

/** Return whether or not the data is packed into a compressed
    binary array */
bool SimStore::isPackedToMemory() const
{
    return not compressed_data.isEmpty();
}

/** Return whether or not the data is packed to disk */
bool SimStore::isPackedToDisk() const
{
    return packed_file.get() != 0;
}

/** Return whether or not the data is packed into memory or 
    is packed onto disk */
bool SimStore::isPacked() const
{
    return this->isPackedToMemory() or this->isPackedToDisk();
}

/** Internal function used to move the data from the packed file
    to memory */
void SimStore::_pvt_moveFromDiskToMemory()
{
    BOOST_ASSERT( packed_file.get() != 0 );
    BOOST_ASSERT( not packed_filename.isEmpty() );

    //open the packed data file using a separate file handle
    QFile f( packed_filename );
    
    if (not f.open(QIODevice::ReadOnly | QIODevice::Unbuffered))
    {
        throw SireError::file_error(f, CODELOC);
    }

    //now read all of the data
    compressed_data = f.readAll();
    
    f.close();

    //lose the file (this will auto-delete the temporary file
    //if no other SimStores are using it)
    packed_file.reset();
    packed_filename = QString::null;
}

/** Pack the system and moves to memory - this will compress
    them to a compressed binary array if they are not packed,
    will do nothing if they are already packed to memory,
    or will move them from disk to memory if they are already
    packed to disk */
void SimStore::packToMemory()
{
    if (this->isPackedToMemory())
    {
        return;
    }
    else if (this->isPackedToDisk())
    {
        this->_pvt_moveFromDiskToMemory();
    }
    else
    {
        QByteArray data;
        
        //start by reserving 32 MB
        data.reserve( 32L*1024L*1024L );
        
        QDataStream ds( &data, QIODevice::WriteOnly );
        
        SharedDataStream sds(ds);
        
        sds << sim_system << sim_moves;
        
        sim_system = System();
        sim_moves = MovesPtr();

        compressed_data = qCompress(data);
        
        packed_file.reset();
        packed_filename = QString::null;
    }
    
    last_packing_state = TO_MEMORY;
}

static QString getUserName()
{
    #ifdef Q_OS_UNIX
        return std::getenv("USER");
    #else
        return "USER";
    #endif
}

/** Pack the system and moves to disk - this places the data
    into a temporary file in 'tempdir' */
void SimStore::packToDisk(const QString &tempdir)
{
    if (this->isPackedToDisk())
        return;

    if (tempdir != QDir::tempPath())
        //save the desired packing directory
        packed_dir = tempdir;

    QDir dir(tempdir);
    
    if (not dir.exists())
        dir = QDir::temp();
        
    boost::shared_ptr<QTemporaryFile> tmp;
    
    tmp.reset( new QTemporaryFile( QString("%1/%2_sire_simstore_XXXXXX.data")
                            .arg(dir.absolutePath(), getUserName()) ) );
                            
    if (not tmp->open())
    {
        if (tempdir == QDir::tempPath())
            throw SireError::file_error( QObject::tr(
                "Cannot pack a SimStore to disk as the temporary directory "
                "(%1) is not writable. Check that you have permission and "
                "there is sufficient space.")
                    .arg(tempdir), CODELOC );

        this->packToDisk();
        return;
    }
    
    //pack the data to memory (if it isn't already)
    this->packToMemory();

    //now write the data to disk
    qint64 nbytes = tmp->write(compressed_data);

    //save the filename as some Qt versions lose the
    //filename when the QTemporaryFile is closed
    packed_filename = tmp->fileName();

    tmp->close();
    
    if (nbytes == -1)
    {
        throw SireError::file_error( QObject::tr(
            "There was an error writing the packed SimStore to disk. "
            "Maybe you don't have write permission or maybe there is not "
            "enough disk space (we require at least %1 MB)")
                .arg( compressed_data.count() / (1024.0*1024.0) ), CODELOC );
    }
    else if (nbytes != compressed_data.count())
    {
        throw SireError::file_error( QObject::tr(
            "There was a problem writing the packed SimStore to disk. "
            "Only %1 bytes were written, out of the necessary %2 bytes. "
            "Maybe the disk became full while it was being written?")
                .arg(nbytes).arg(compressed_data.count()), CODELOC );
    }
    
    packed_file = tmp;
    compressed_data = QByteArray();
    
    last_packing_state = TO_DISK;
}

/** Pack the system and moves to disk - this places the data in 
    a temporary file in QDir::tempPath() */
void SimStore::packToDisk()
{
    if (packed_dir.isEmpty())
        this->packToDisk(QDir::tempPath());
    else
        this->packToDisk(packed_dir);
}

/** Pack the system - this packs the system to the same state
    as the last time this function was called (so if it was previously
    packed to disk, then this will pack it to disk). This allows 
    you to call SimStore::unpack(), knowing that SimStore::pack()
    will restore the packing state */
void SimStore::pack()
{
    if (this->isPacked())
        //the data is already packed
        return;
    
    switch (last_packing_state)
    {
        case NOT_PACKED:
        case TO_MEMORY:
            this->packToMemory();
            break;
            
        case TO_DISK:
            this->packToDisk();
            break;
    }
}

/** Unpack the system and move from the compressed binary array */
void SimStore::unpack()
{
    if (this->isPackedToDisk())
        this->_pvt_moveFromDiskToMemory();

    if (this->isPackedToMemory())
    {
        System new_system;
        MovesPtr new_moves;

        //use a local scope so that the uncompressed data is deleted
        //as soon as possible
        {
            QByteArray data = qUncompress( compressed_data );
            
            QDataStream ds(data);
            SharedDataStream sds(ds);
        
            sds >> new_system >> new_moves;
        }

        compressed_data = QByteArray();
        sim_system = new_system;
        sim_moves = new_moves;
        packed_file.reset();
        packed_filename = QString::null;
    }
}

/** Set the system to be stored */
void SimStore::setSystem(const System &system)
{
    if (this->isPacked())
        throw SireError::invalid_state( QObject::tr(
            "Cannot set the system in a SimStore that is packed. You must "
            "unpack it first."), CODELOC );
    
    sim_system = system;
}

/** Set the moves to be stored */
void SimStore::setMoves(const Moves &moves)
{
    if (this->isPacked())
        throw SireError::invalid_state( QObject::tr(
            "Cannot set the moves in a SimStore that is packed. You must "
            "unpack it first."), CODELOC );
        
    sim_moves = moves;
}

/** Set both the system and moves to be stored */
void SimStore::setSystemAndMoves(const System &system, const Moves &moves)
{
    if (this->isPacked())
        throw SireError::invalid_state( QObject::tr(
            "Cannot set the system and moves in a SimStore that is packed. You must "
            "unpack it first."), CODELOC );

    sim_system = system;
    sim_moves = moves;
}

/** Return a copy of the system being stored */
const System& SimStore::system() const
{
    if (this->isPacked())
        throw SireError::invalid_state( QObject::tr(
            "You cannot get the system from a SimStore that is packed. You must "
            "unpack it first."), CODELOC );

    return sim_system;
}

/** Return a copy of the moves being stored */
const Moves& SimStore::moves() const
{
    if (this->isPacked())
        throw SireError::invalid_state( QObject::tr(
            "You cannot get the moves from a SimStore that is packed. You must "
            "unpack it first."), CODELOC );

    return sim_moves;
}

const char* SimStore::typeName()
{
    return QMetaType::typeName( qMetaTypeId<SimStore>() );
}

SimStore* SimStore::clone() const
{
    return new SimStore(*this);
}
