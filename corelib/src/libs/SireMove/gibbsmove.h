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

#ifndef SIREMOVE_GIBBSMOVE_H
#define SIREMOVE_GIBBSMOVE_H

#include "montecarlo.h"

SIRE_BEGIN_HEADER

namespace SireMove
{
class GibbsMove;
}

QDataStream& operator<<(QDataStream&, const SireMove::GibbsMove&);
QDataStream& operator>>(QDataStream&, SireMove::GibbsMove&);

namespace SireMove
{

/** This is a Monte Carlo move that is used to perform particle
    deletion from one box, and particle insertion into another box,
    so sampling from the Gibbs ensemble;
    
    "Phase equilibria by simulation in the Gibbs ensemble"
    A. Z. Panagiotopoulos, N. Quirke, M. Stapleton, D. J. Tildesley
    
    Mol. Phys., 1988, VOL. 63, NO. 4, 527-545

    The move will delete a molecule from one set of molecule groups
    and will add it with a random position and orientation to another
    set of molecule groups. The two sets of molecule groups must
    correspond to the two simulation boxes use to generate the 
    Gibbs ensemble, and the energy component this move samples
    must equal the sum of the energies of the two boxes, and the
    space property this move uses must be a SireVol::CombinedSpace
    which combines the spaces used by the two boxes (with the
    first space in the combined space representing the first box,
    and the second space in the combined space representing
    the second box). The two spaces in the combined spaces
    must be periodic (as we need a fixed space in which to 
    add new molecules, and to maintain constant volume
    or pressure).

	@author Christopher Woods
*/
class SIREMOVE_EXPORT GibbsMove
          : public SireBase::ConcreteProperty<GibbsMove,MonteCarlo>
{

friend SIREMOVE_EXPORT QDataStream& ::operator<<(QDataStream&, const GibbsMove&);
friend SIREMOVE_EXPORT QDataStream& ::operator>>(QDataStream&, GibbsMove&);

public:
    GibbsMove(const PropertyMap &map = PropertyMap());
    
    GibbsMove(const MGIDsAndMaps &group0, const MGIDsAndMaps &group1,
              const PropertyMap &map = PropertyMap());
              
	GibbsMove(const GibbsMove &other);
    
    ~GibbsMove();
    
	GibbsMove& operator=(const GibbsMove &other);
    
	bool operator==(const GibbsMove &other) const;
    bool operator!=(const GibbsMove &other) const;
    
    static const char* typeName();
        
    bool operator==(const VolumeMove &other) const;
    bool operator!=(const VolumeMove &other) const;
    
    QString toString() const;

    const MGIDsAndMaps& group0() const;
    const MGIDsAndMaps& group1() const;
    
    void setGroup0(const MGIDsAndMaps &group0);
    void setGroup1(const MGIDsAndMaps &group1);
    
    void setGenerator(const RanGenerator &rangenerator);
    
    void move(System &system, int nmoves, bool record_stats=true);

protected:
    void _pvt_setTemperature(const SireUnits::Dimension::Temperature &temperature);
    void _pvt_setPressure(const SireUnits::Dimension::Pressure &pressure);

private:
    /** The molecule group IDs and associated property maps
        for the two groups */
    MGIDsAndMaps g0, g1;
    
    /** The set of coordinate properties used by this move
        for each of the two groups - this is used as random
        coordinates need to be generated for each property */
    QVector<PropertyName> coord_props0, coord_props1;
};

}

SIRE_END_HEADER

#endif
