/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2017  Christopher Woods
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

#ifndef SIREIO_AMBERRST7_H
#define SIREIO_AMBERRST7_H

#include "moleculeparser.h"

#include "SireMaths/vector.h"

SIRE_BEGIN_HEADER

namespace SireIO
{
class AmberRst7;
}

QDataStream& operator<<(QDataStream&, const SireIO::AmberRst7&);
QDataStream& operator>>(QDataStream&, SireIO::AmberRst7&);

namespace SireIO
{

/** This class represents an Amber-format restart/coordinate file (ascii),
    currently supporting these files from Amber7 to Amber16.
    
    The format of this file is described here;
    
    http://ambermd.org/formats.html
    
    (specifically the "AMBER coordinate/restart" file specification
 
    @author Christopher Woods
*/
class SIREIO_EXPORT AmberRst7 : public SireBase::ConcreteProperty<AmberRst7,MoleculeParser>
{

friend QDataStream& ::operator<<(QDataStream&, const AmberRst7&);
friend QDataStream& ::operator>>(QDataStream&, AmberRst7&);

public:
    AmberRst7();
    AmberRst7(const QString &filename,
              const PropertyMap &map = PropertyMap());
    AmberRst7(const QStringList &lines,
              const PropertyMap &map = PropertyMap());
    AmberRst7(const SireSystem::System &system,
              const PropertyMap &map = PropertyMap());
    
    AmberRst7(const AmberRst7 &other);
    
    ~AmberRst7();
    
    AmberRst7& operator=(const AmberRst7 &other);
    
    bool operator==(const AmberRst7 &other) const;
    bool operator!=(const AmberRst7 &other) const;
    
    static const char* typeName();
    
    const char* what() const;
    
    MoleculeParserPtr construct(const QString &filename,
                                const PropertyMap &map) const;

    MoleculeParserPtr construct(const QStringList &lines,
                                const PropertyMap &map) const;
    
    MoleculeParserPtr construct(const SireSystem::System &system,
                                const PropertyMap &map) const;

    QString toString() const;

    QString formatName() const;
    QString formatDescription() const;
    QStringList formatSuffix() const;

    static AmberRst7 parse(const QString &filename);

    QString title() const;
    double time() const;

    int nAtoms() const;
    
    bool hasVelocities() const;

    QVector<SireMaths::Vector> coordinates() const;
    QVector<SireMaths::Vector> velocities() const;

    SireMaths::Vector boxDimensions() const;
    SireMaths::Vector boxAngles() const;

protected:
    void addToSystem(SireSystem::System &system, const PropertyMap &map) const;

private:
    void parse(const PropertyMap &map);

    void readBoxInfo(int boxidx);

    /** The title of the file */
    QString ttle;
    
    /** The current time of the simulation in picoseconds */
    double current_time;
    
    /** The coordinate data */
    QVector<SireMaths::Vector> coords;
    
    /** The velocity data in amber units (angstroms / 1/20.455 ps) */
    QVector<SireMaths::Vector> vels;
    
    /** The box dimensions */
    SireMaths::Vector box_dims;
    
    /** The box angles */
    SireMaths::Vector box_angs;
};

}

Q_DECLARE_METATYPE( SireIO::AmberRst7 )

SIRE_EXPOSE_CLASS( SireIO::AmberRst7 )

SIRE_END_HEADER

#endif
