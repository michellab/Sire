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

#ifndef SIREVOL_PATCHING_H
#define SIREVOL_PATCHING_H

#include "space.h"

#include <QPair>

namespace SireVol
{
class Patching;
class NullPatching;

class BoxPatching;
}

SIREVOL_EXPORT QDataStream& operator<<(QDataStream&, const SireVol::Patching&);
SIREVOL_EXPORT QDataStream& operator>>(QDataStream&, SireVol::Patching&);

SIREVOL_EXPORT QDataStream& operator<<(QDataStream&, const SireVol::NullPatching&);
SIREVOL_EXPORT QDataStream& operator>>(QDataStream&, SireVol::NullPatching&);

SIREVOL_EXPORT QDataStream& operator<<(QDataStream&, const SireVol::BoxPatching&);
SIREVOL_EXPORT QDataStream& operator>>(QDataStream&, SireVol::BoxPatching&);

namespace SireVol
{

typedef SireBase::PropPtr<Patching> PatchingPtr;

/** This is the base class of all Patching classes. Patching
    represents a scheme for decomposing a space into a set
    of regions (domain decomposition) that contain neighbouring
    CoordGroups - this allows inter-CoordGroup calculations to
    be accelerated as cutoff tests can be applied to the patches
    to eliminate tests of the contained CoordGroups.
    
    @author Christopher Woods
*/
class SIREVOL_EXPORT Patching : public SireBase::Property
{

friend QDataStream& ::operator<<(QDataStream&, const Patching&);
friend QDataStream& ::operator>>(QDataStream&, Patching&);

public:
    Patching();
    Patching(const Patching &other);
    
    virtual ~Patching();
    
    static const char* typeName();
    
    virtual Patching* clone() const=0;
    
    const Space& space() const;
    
    virtual int nPatches() const=0;
    
    virtual int patchIndex(const Vector &point) const=0;

    virtual QPair<int,Vector> patchIndexAndCenter(const Vector &point) const=0;

    virtual PatchingPtr repatch(const Space &new_space) const=0;

    virtual PatchingPtr rebalance(const Space &space,
                                  const QVector<CoordGroupArray> &patchcoords) const;

    static NullPatching null();

protected:
    Patching(const Space &space);
    
    Patching& operator=(const Patching &other);
    
    bool operator==(const Patching &other) const;
    bool operator!=(const Patching &other) const;

private:
    /** The space from which the patching was derived */
    SpacePtr spce;
};

/** Null patching */
class SIREVOL_EXPORT NullPatching 
        : public SireBase::ConcreteProperty<NullPatching,Patching>
{

friend QDataStream& ::operator<<(QDataStream&, const NullPatching&);
friend QDataStream& ::operator>>(QDataStream&, NullPatching&);

public:
    NullPatching();
    NullPatching(const Space &space);
    NullPatching(const NullPatching &other);
    
    ~NullPatching();
    
    static const char* typeName();
    
    NullPatching& operator=(const NullPatching &other);
    
    bool operator==(const NullPatching &other) const;
    bool operator!=(const NullPatching &other) const;
    
    int nPatches() const;
    
    int patchIndex(const Vector &point) const;

    QPair<int,Vector> patchIndexAndCenter(const Vector &point) const;
    
    PatchingPtr repatch(const Space &new_space) const;
};

/** This is a simple patching scheme that divides space up into a 
    series of cuboidal boxes (3D grid!)
    
    @author Christopher Woods
*/
class SIREVOL_EXPORT BoxPatching
        : public SireBase::ConcreteProperty<BoxPatching,Patching>
{

friend QDataStream& ::operator<<(QDataStream&, const BoxPatching&);
friend QDataStream& ::operator>>(QDataStream&, BoxPatching&);

public:
    BoxPatching();
    BoxPatching(const Space &space);
    BoxPatching(const Space &space, const Vector &center);
    
    BoxPatching(const Space &space, SireUnits::Dimension::Length patch_size);
    BoxPatching(const Space &space, SireUnits::Dimension::Length patch_size,
                const Vector &center);
    
    BoxPatching(const BoxPatching &other);
    
    ~BoxPatching();
    
    static const char* typeName();
    
    BoxPatching& operator=(const BoxPatching &other);
    
    bool operator==(const BoxPatching &other) const;
    bool operator!=(const BoxPatching &other) const;

    QString toString() const;

    Vector center() const;
    SireUnits::Dimension::Length patchSize() const;
    
    Vector patchDimension() const;
    
    AABox patchBox(int i) const;
    AABox patchBox(const Vector &point) const;
    
    int nPatches() const;
    
    int patchIndex(const Vector &point) const;

    QPair<int,Vector> patchIndexAndCenter(const Vector &point) const;
    
    PatchingPtr repatch(const Space &new_space) const;
    
private:
    int getIndex(const Vector &point) const;

    /** The initial dimension used to create the patch */
    SireUnits::Dimension::Length patch_size;

    /** The origin of the 3D grid */
    Vector orgn;
    
    /** The inverse unit vector of the 3D grid (the grid is centered
        around (0,0,0) */
    Vector inv_gridvec;
        
    /** The number of patches in the X, Y and Z dimensions */
    int nx, ny, nz;
};

}

Q_DECLARE_METATYPE( SireVol::NullPatching )
Q_DECLARE_METATYPE( SireVol::BoxPatching )

SIRE_EXPOSE_CLASS( SireVol::Patching )
SIRE_EXPOSE_CLASS( SireVol::NullPatching )
SIRE_EXPOSE_CLASS( SireVol::BoxPatching )

SIRE_EXPOSE_PROPERTY( SireVol::PatchingPtr, SireVol::Patching )

#endif
