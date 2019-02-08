/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2008  Christopher Woods
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

#ifndef SIREMOVE_ZMATMOVE_H
#define SIREMOVE_ZMATMOVE_H

#include "montecarlo.h"
#include "sampler.h"

SIRE_BEGIN_HEADER

namespace SireMove
{
class ZMatMove;
}

SIREMOVE_EXPORT QDataStream& operator<<(QDataStream&, const SireMove::ZMatMove&);
SIREMOVE_EXPORT QDataStream& operator>>(QDataStream&, SireMove::ZMatMove&);

namespace SireMol
{
class MoleculeGroup;
class PartialMolecule;
class AtomIdx;
}

namespace SireMove
{

class Sampler;
class ZMatrixCoords;

using SireMol::MoleculeGroup;
using SireMol::PartialMolecule;

/** This class implements a z-matrix based intramolecular Monte Carlo
    move that may be applied to a random molecule (or part of a molecule)
    within a MoleculeGroup
    
    @author Christopher Woods
*/
class SIREMOVE_EXPORT ZMatMove 
            : public SireBase::ConcreteProperty<ZMatMove,MonteCarlo>
{

friend SIREMOVE_EXPORT QDataStream& ::operator<<(QDataStream&, const ZMatMove&);
friend SIREMOVE_EXPORT QDataStream& ::operator>>(QDataStream&, ZMatMove&);

public:
    ZMatMove(const PropertyMap &map = PropertyMap());
    
    ZMatMove(const MoleculeGroup &molgroup,
             const PropertyMap &map = PropertyMap());
    ZMatMove(const Sampler &sampler,
             const PropertyMap &map = PropertyMap());
    
    ZMatMove(const ZMatMove &other);
    
    ~ZMatMove();
    
    ZMatMove& operator=(const ZMatMove &other);
    
    static const char* typeName();

    bool operator==(const ZMatMove &other) const;
    bool operator!=(const ZMatMove &other) const;

    QString toString() const;

    void setSampler(const Sampler &sampler);
    void setSampler(const MoleculeGroup &molgroup);

    const Sampler& sampler() const;
    const MoleculeGroup& moleculeGroup() const;

    const PropertyName& zmatrixProperty() const;

    void setZMatrixProperty(const PropertyName &property);

    void setSynchronisedMotion(bool on);
    void setSynchronisedBonds(bool on);
    void setSynchronisedAngles(bool on);
    void setSynchronisedDihedrals(bool on);
    
    bool synchronisedMotion() const;
    bool synchronisedBonds() const;
    bool synchronisedAngles() const;
    bool synchronisedDihedrals() const;
    
    void setGenerator(const RanGenerator &rangenerator);

    void move(System &system, int nmoves, bool record_stats=true);

protected:
    void _pvt_setTemperature(const SireUnits::Dimension::Temperature &temperature);

private:
    void move(SireMol::AtomIdx atom, ZMatrixCoords &zmatrix,
              QHash< SireMol::AtomIdx, boost::tuple<SireUnits::Dimension::Length,
                     SireUnits::Dimension::Angle,SireUnits::Dimension::Angle> > &deltas);

    /** The sampler used to select random molecules for the move */
    SamplerPtr smplr;
    
    /** The name of the property that contains the z-matrix */
    PropertyName zmatrix_property;
    
    /** Whether or not to move equivalent bonds in different
        molecules by the same amount (equivalent is judged by 
        AtomIdx) */
    bool sync_bonds;
    
    /** Whether or not to move equivalent angles in different
        molecules by the same amount */
    bool sync_angles;
    
    /** Whether or not to move equivalent dihedrals in different
        molecules by the same amount */
    bool sync_dihedrals;
};

}

Q_DECLARE_METATYPE( SireMove::ZMatMove )

SIRE_EXPOSE_CLASS( SireMove::ZMatMove )

SIRE_END_HEADER

#endif
