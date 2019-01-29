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

#ifndef SIREIO_TRAJECTORYMONITOR_H
#define SIREIO_TRAJECTORYMONITOR_H

#include <QList>

#include <boost/shared_ptr.hpp>

#include "SireSystem/systemmonitor.h"

#include "SireVol/space.h"

#include "SireMol/mgidentifier.h"

#include "SireBase/propertymap.h"

#include "iobase.h"

SIRE_BEGIN_HEADER

namespace SireIO
{
class TrajectoryMonitor;
}

SIREIO_EXPORT QDataStream& operator<<(QDataStream&, const SireIO::TrajectoryMonitor&);
SIREIO_EXPORT QDataStream& operator>>(QDataStream&, SireIO::TrajectoryMonitor&);

class QTemporaryFile;

namespace SireIO
{

using SireSystem::System;

using SireMol::MoleculeGroup;
using SireMol::MGIdentifier;
using SireMol::MGID;

using SireBase::PropertyMap;

/** This is a monitor that can be used to save a trajectory
    of an arbitrary collection of molecules from the system
    
    @author Christopher Woods
*/
class SIREIO_EXPORT TrajectoryMonitor 
    : public SireBase::ConcreteProperty<TrajectoryMonitor,SireSystem::SystemMonitor>
{

friend QDataStream& ::operator<<(QDataStream&, const TrajectoryMonitor&);
friend QDataStream& ::operator>>(QDataStream&, TrajectoryMonitor&);

public:
    TrajectoryMonitor();
    
    TrajectoryMonitor(const MoleculeGroup &molgroup,
                      const PropertyMap &map = PropertyMap());
    
    TrajectoryMonitor(const MoleculeGroup &molgroup, const IOBase &writer,
                      const PropertyMap &map = PropertyMap());
    
    TrajectoryMonitor(const MGID &mgid, const PropertyMap &map = PropertyMap());
    
    TrajectoryMonitor(const MGID &mgid, const IOBase &writer,
                      const PropertyMap &map = PropertyMap());
    
    TrajectoryMonitor(const TrajectoryMonitor &other);
    
    ~TrajectoryMonitor();
    
    TrajectoryMonitor& operator=(const TrajectoryMonitor &other);
    
    bool operator==(const TrajectoryMonitor &other) const;
    bool operator!=(const TrajectoryMonitor &other) const;
    
    static const char* typeName();

    void clearStatistics();
    
    void monitor(System &system);
    
    void setTempDir(const QString &tempdir);
    
    void writeToDisk(const QString &file_template) const;
    
private:
    /** The IO object used to create the coordinates file(s) 
        from the system */
    IOPtr io_writer;
    
    /** The ID of the molecule group that is being written */
    MGIdentifier mgid;
    
    /** Temporary files containing each frame of the animation */
    QList< QPair< QString,boost::shared_ptr<QTemporaryFile> > > traj_frames;
    
    /** The system space for each frame of the trajectory */
    QList<SireVol::SpacePtr> space_frames;
    
    /** The property map used to find the properties that will
        be written to the trajectory frame */
    PropertyMap mol_properties;
    
    /** Name of the directory in which to save the temporary files */
    QString temp_dir;
};

}

Q_DECLARE_METATYPE( SireIO::TrajectoryMonitor )

SIRE_EXPOSE_CLASS( SireIO::TrajectoryMonitor )

SIRE_END_HEADER

#endif
