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

#ifndef SIREIO_DCD_H
#define SIREIO_DCD_H

#include "moleculeparser.h"

#include "SireMol/trajectory.h"

#include "SireMaths/vector.h"

SIRE_BEGIN_HEADER

namespace SireIO
{
class DCD;

namespace detail{ class DCDFile; }

}

SIREIO_EXPORT QDataStream& operator<<(QDataStream&, const SireIO::DCD&);
SIREIO_EXPORT QDataStream& operator>>(QDataStream&, SireIO::DCD&);

namespace SireIO
{

class FortranFile;

namespace detail
{
/** This class provides a low-level interface to reading and writing
 *  a DCD file. It is designed to be used only with the
 *  DCD class
 */
class DCDFile
{
public:
    DCDFile();
    DCDFile(const QString &filename);
    ~DCDFile();

    void readHeader(FortranFile &file);

    SireVol::SpacePtr readSpace(FortranFile &file, int frame) const;
    QVector<SireMaths::Vector> readCoordinates(FortranFile &file, int frame) const;

    SireMol::Frame readFrame(FortranFile &file, int frame) const;

    void setTitle(QString title);
    QString getTitle() const;

    double getTimeStep() const;
    qint64 getFrameStart() const;
    qint64 getFrameDelta() const;

    double getCurrentTime() const;
    void setCurrentTime(double time);

    void setSpace(const SireVol::Space &space);
    const SireVol::Space& getSpace() const;

    qint64 nFrames() const;
    qint64 nAtoms() const;

    double getTimeAtFrame(int frame) const;

private:
    QStringList title;

    QVector<qint32> fixed_atoms;
    QVector<SireMaths::Vector> first_frame;

    double timestep;

    SireVol::SpacePtr spc;

    qint64 istart;
    qint64 nsavc;
    qint64 nfixed;

    qint64 natoms;
    qint64 nframes;

    qint64 first_frame_line;

    bool CHARMM_FORMAT;
    bool HAS_EXTRA_BLOCK;
    bool HAS_FOUR_DIMS;
};
}

/** This class represents a DCD file reader.

    The format is described here

    https://www.ks.uiuc.edu/Research/vmd/plugins/molfile/dcdplugin.html

    @author Christopher Woods
*/
class SIREIO_EXPORT DCD : public SireBase::ConcreteProperty<DCD,MoleculeParser>
{

friend QDataStream& ::operator<<(QDataStream&, const DCD&);
friend QDataStream& ::operator>>(QDataStream&, DCD&);

public:
    DCD();
    DCD(const QString &filename,
        const PropertyMap &map = PropertyMap());
    DCD(const QStringList &lines,
        const PropertyMap &map = PropertyMap());
    DCD(const SireSystem::System &system,
        const PropertyMap &map = PropertyMap());

    DCD(const DCD &other);

    ~DCD();

    DCD& operator=(const DCD &other);

    bool operator==(const DCD &other) const;
    bool operator!=(const DCD &other) const;

    static const char* typeName();

    const char* what() const;

    MoleculeParserPtr construct(const QString &filename,
                                const PropertyMap &map) const;

    MoleculeParserPtr construct(const QStringList &lines,
                                const PropertyMap &map) const;

    MoleculeParserPtr construct(const SireSystem::System &system,
                                const PropertyMap &map) const;

    QString toString() const;

    bool isTextFile() const;

    QString formatName() const;
    QString formatDescription() const;
    QStringList formatSuffix() const;

    static DCD parse(const QString &filename);

    int nAtoms() const;

    bool isFrame() const;

    int nFrames() const;
    SireMol::Frame getFrame(int i) const;

    QString title() const;

    QVector<SireMaths::Vector> coordinates() const;

    const SireVol::Space& space() const;

    SireUnits::Dimension::Time time() const;

    QStringList warnings() const;

    void writeToFile(const QString &filename) const;

protected:
    void addToSystem(SireSystem::System &system, const PropertyMap &map) const;

private:
    void parse(const QString &filename, const PropertyMap &map);

    /** The coordinate data */
    QVector<SireMaths::Vector> coords;

    /** Handle to read in the DCD data */
    detail::DCDFile dcd;

    /** Any warnings that were raised when reading the file */
    QStringList parse_warnings;
};

}

Q_DECLARE_METATYPE( SireIO::DCD )

SIRE_EXPOSE_CLASS( SireIO::DCD )

SIRE_END_HEADER

#endif
