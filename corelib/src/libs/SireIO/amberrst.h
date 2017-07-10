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

#ifndef SIREIO_AMBERRST_H
#define SIREIO_AMBERRST_H

#include "moleculeparser.h"

#include "SireMaths/vector.h"

SIRE_BEGIN_HEADER

namespace SireIO
{
class AmberRst;
}

QDataStream& operator<<(QDataStream&, const SireIO::AmberRst&);
QDataStream& operator>>(QDataStream&, SireIO::AmberRst&);

namespace SireIO
{

/** This class represents an Amber-format restart/coordinate file (binary),
    currently supporting these files from Amber7 to Amber16.
    
    This is a netcdf file format, which is described here;
    
    http://ambermd.org/netcdf/nctraj.xhtml
 
    @author Christopher Woods
*/
class SIREIO_EXPORT AmberRst : public SireBase::ConcreteProperty<AmberRst,MoleculeParser>
{

friend QDataStream& ::operator<<(QDataStream&, const AmberRst&);
friend QDataStream& ::operator>>(QDataStream&, AmberRst&);

public:
    AmberRst();
    AmberRst(const QString &filename,
             const PropertyMap &map = PropertyMap());
    AmberRst(const QStringList &lines,
             const PropertyMap &map = PropertyMap());
    AmberRst(const SireSystem::System &system,
             const PropertyMap &map = PropertyMap());
    
    AmberRst(const AmberRst &other);
    
    ~AmberRst();
    
    AmberRst& operator=(const AmberRst &other);
    
    bool operator==(const AmberRst &other) const;
    bool operator!=(const AmberRst &other) const;
    
    static const char* typeName();
    
    const char* what() const;

    QString toString() const;

    QString formatName() const;

    static AmberRst parse(const QString &filename);

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

    /** The title of the file */
    QString ttle;
    
    /** The current time of the simulation */
    double current_time;
    
    /** The coordinate data */
    QVector<SireMaths::Vector> coords;
    
    /** The velocity data in amber units */
    QVector<SireMaths::Vector> vels;
    
    /** The box dimensions */
    SireMaths::Vector box_dims;
    
    /** The box angles */
    SireMaths::Vector box_angs;
};

}

Q_DECLARE_METATYPE( SireIO::AmberRst )

SIRE_EXPOSE_CLASS( SireIO::AmberRst )

SIRE_END_HEADER

#endif
