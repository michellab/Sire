/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2010  Christopher Woods
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

#ifndef SIREIO_CUBE_H
#define SIREIO_CUBE_H

#include "SireFF/potentialtable.h"
#include "SireMol/moleculegroups.h"
#include "SireMol/moleculegroup.h"

#include "SireUnits/dimensions.h"

SIRE_BEGIN_HEADER

namespace SireIO
{

using SireBase::PropertyMap;
using SireFF::PotentialTable;
using SireMol::MolGroupsBase;
using SireMol::MoleculeGroup;

/** This class is used to write a PotentialTable as a 
    Gaussian cube file
    
    @author Christopher Woods
*/
class SIREIO_EXPORT Cube
{
public:
    Cube();
    Cube(SireUnits::Dimension::MolarEnergy cutoff);
    Cube(const Cube &other);
    
    ~Cube();
    
    Cube& operator=(const Cube &other);
    
    bool operator==(const Cube &other) const;
    bool operator!=(const Cube &other) const;

    void write(const PotentialTable &table,
               const QString &filename, const PropertyMap &map = PropertyMap()) const;
    
    void write(const PotentialTable &table, const MoleculeGroup &molgroup,
               const QString &filename, const PropertyMap &map = PropertyMap()) const;
    
    void write(const PotentialTable &table, const MolGroupsBase &molecules,
               const QString &filename, const PropertyMap &map = PropertyMap()) const;

private:
    double cutoff;
};

}

SIRE_EXPOSE_CLASS( SireIO::Cube )

SIRE_END_HEADER

#endif
