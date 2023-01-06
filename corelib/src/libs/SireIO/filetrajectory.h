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

#ifndef SIREIO_FILETRAJECTORY_H
#define SIREIO_FILETRAJECTORY_H

#include "moleculeparser.h"

#include "SireMol/trajectory.h"

SIRE_BEGIN_HEADER

namespace SireIO
{
class FileTrajectory;
}

SIREIO_EXPORT QDataStream& operator<<(QDataStream&, const SireIO::FileTrajectory&);
SIREIO_EXPORT QDataStream& operator>>(QDataStream&, SireIO::FileTrajectory&);

namespace SireIO
{

/** This is a trajectory that loads frames dynamically from the
 *  filesystem using MoleculeParser objects that implement
 *  the trajectory API
 */
class SIREIO_EXPORT FileTrajectory : public SireMol::TrajectoryData
{

friend QDataStream& ::operator<<(QDataStream&, const FileTrajectory&);
friend QDataStream& ::operator>>(QDataStream&, FileTrajectory&);

public:
    FileTrajectory();
    FileTrajectory(const MoleculeParser &parser);

    FileTrajectory(const FileTrajectory &other);

    ~FileTrajectory();

    const char* what() const;

    static const char* typeName();

    FileTrajectory* clone() const;

    int nFrames() const;
    int nAtoms() const;

    QStringList filenames() const;

    SireMol::Frame getFrame(int i) const;

    bool isEditable() const;

protected:
    bool _equals(const TrajectoryData &other) const;

private:
    MoleculeParserPtr parser;
};

}

Q_DECLARE_METATYPE(SireIO::FileTrajectory)

SIRE_END_HEADER

#endif
