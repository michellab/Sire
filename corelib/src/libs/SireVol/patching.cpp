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

#include "patching.h"

#include "SireVol/periodicbox.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireUnits/units.h"

#include "SireError/errors.h"

using namespace SireVol;
using namespace SireBase;
using namespace SireUnits;
using namespace SireUnits::Dimension;
using namespace SireStream;

/////////////
///////////// Implementation of Patching
/////////////

static const RegisterMetaType<Patching> r_patching( MAGIC_ONLY,
                                                    Patching::typeName() );
                                                    
QDataStream SIREVOL_EXPORT &operator<<(QDataStream &ds, const Patching &patching)
{
    writeHeader(ds, r_patching, 1);
    
    SharedDataStream sds(ds);
    
    sds << patching.spce << static_cast<const Property&>(patching);
    
    return ds;
}

QDataStream SIREVOL_EXPORT &operator>>(QDataStream &ds, Patching &patching)
{
    VersionID v = readHeader(ds, r_patching);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> patching.spce >> static_cast<Property&>(patching);
    }
    else
        throw version_error(v, "1", r_patching, CODELOC);
        
    return ds;
}

/** Null constructor */
Patching::Patching() : Property()
{}

/** Internal constructor used to set the space */
Patching::Patching(const Space &space) : Property(), spce(space)
{}

/** Copy constructor */
Patching::Patching(const Patching &other)
         : Property(other), spce(other.spce)
{}

/** Destructor */
Patching::~Patching()
{}

/** Copy assignment operator */
Patching& Patching::operator=(const Patching &other)
{
    if (this != &other)
    {
        spce = other.spce;
        Property::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool Patching::operator==(const Patching &other) const
{
    return spce == other.spce;
}

/** Comparison operator */
bool Patching::operator!=(const Patching &other) const
{
    return not Patching::operator==(other);
}

const char* Patching::typeName()
{
    return "SireVol::Patching";
}
    
/** Return the space used to create this patching scheme */
const Space& Patching::space() const
{
    return spce.read();
}

/** Rebalance the patching so that the patches for the passed space contain
    roughly equal numbers of CoordGroups */
PatchingPtr Patching::rebalance(const Space &space,
                                const QVector<CoordGroupArray> &patchcoords) const
{
    return PatchingPtr(*this);
}

/** Retunr the null patching object */
NullPatching Patching::null()
{
    return NullPatching();
}

/////////////
///////////// Implementation of NullPatching
/////////////

static const RegisterMetaType<NullPatching> r_nullpatching;

QDataStream SIREVOL_EXPORT &operator<<(QDataStream &ds, const NullPatching &nullpatching)
{
    writeHeader(ds, r_nullpatching, 1);
    
    ds << static_cast<const Patching&>(nullpatching);
    
    return ds;
}

QDataStream SIREVOL_EXPORT &operator>>(QDataStream &ds, NullPatching &nullpatching)
{
    VersionID v = readHeader(ds, r_nullpatching);
    
    if (v == 1)
    {
        ds >> static_cast<Patching&>(nullpatching);
    }
    else
        throw version_error(v, "1", r_nullpatching, CODELOC);
        
    return ds;
}

/** Constructor */
NullPatching::NullPatching() : ConcreteProperty<NullPatching,Patching>()
{}

/** Construct with the passed space */
NullPatching::NullPatching(const Space &space)
             : ConcreteProperty<NullPatching,Patching>(space)
{}

/** Copy constructor */
NullPatching::NullPatching(const NullPatching &other)
             : ConcreteProperty<NullPatching,Patching>(other)
{}

/** Destructor */
NullPatching::~NullPatching()
{}

const char* NullPatching::typeName()
{
    return QMetaType::typeName( qMetaTypeId<NullPatching>() );
}

/** Copy assignment operator */
NullPatching& NullPatching::operator=(const NullPatching &other)
{
    Patching::operator=(other);
    return *this;
}

/** Comparison operator */
bool NullPatching::operator==(const NullPatching &other) const
{
    return Patching::operator==(other);
}

/** Comparison operator */
bool NullPatching::operator!=(const NullPatching &other) const
{
    return Patching::operator!=(other);
}

/** Return the number of patches */
int NullPatching::nPatches() const
{
    return 1;
}

/** Return the patch index of the passed point */
int NullPatching::patchIndex(const Vector &point) const
{
    return 0;
}

/** Return the patch index and the center of the patch */
QPair<int,Vector> NullPatching::patchIndexAndCenter(const Vector &point) const
{
    return QPair<int,Vector>(0, Vector(0));
}

/** Repatch this patching for the passed space */
PatchingPtr NullPatching::repatch(const Space &new_space) const
{
    return NullPatching(new_space);
}

/////////////
///////////// Implementation of BoxPatching
/////////////

static const RegisterMetaType<BoxPatching> r_boxpatching;

QDataStream SIREVOL_EXPORT &operator<<(QDataStream &ds, const BoxPatching &boxpatching)
{
    writeHeader(ds, r_boxpatching, 1);
    
    SharedDataStream sds(ds);
    
    sds << boxpatching.patch_size.to(angstrom)
        << boxpatching.orgn << boxpatching.inv_gridvec
        << boxpatching.nx << boxpatching.ny << boxpatching.nz
        << static_cast<const Patching&>(boxpatching);
        
    return ds;
}

QDataStream SIREVOL_EXPORT &operator>>(QDataStream &ds, BoxPatching &boxpatching)
{
    VersionID v = readHeader(ds, r_boxpatching);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        double patch_size;
        
        sds >> patch_size >> boxpatching.orgn >> boxpatching.inv_gridvec
            >> boxpatching.nx >> boxpatching.ny >> boxpatching.nz
            >> static_cast<Patching&>(boxpatching);
           
        boxpatching.patch_size = patch_size*angstrom;
              
        return ds;
    }
    else
        throw version_error(v, "1", r_boxpatching, CODELOC);
        
    return ds;
}

/** Constructor */
BoxPatching::BoxPatching()
            : ConcreteProperty<BoxPatching,Patching>(),
              patch_size(0), orgn(0), inv_gridvec(0), nx(0), ny(0), nz(0)
{}

/** Construct for the passed space, placing the center of the grid at "center",
    trying to construct a grid that divides space using a patch size of
    approximately "patch_size" 
    
    Note that this patching is only compatible with cartesian spaces.
    
    \throw SireError::incompatible_error
*/
BoxPatching::BoxPatching(const Space &space, Length size, const Vector &center)
            : ConcreteProperty<BoxPatching,Patching>(space),
              patch_size(size), orgn(0), inv_gridvec(0), nx(0), ny(0), nz(0)
{
    if (not space.isCartesian())
        throw SireError::incompatible_error( QObject::tr(
                "BoxPatching is only compatible with cartesian spaces "
                "(spaces for which space.isCartesian() is true). The passed "
                "space (%1) is not a cartesian space.")
                    .arg(space.toString()), CODELOC );

    if (space.isPeriodic())
    {
        //need a virtual function call here - as at the moment it
        //depends on all cartesian periodic spaces being PeriodicBox...
        Vector dimensions = space.asA<PeriodicBox>().dimensions();
        
        nx = int( dimensions.x() / size.value() ) + 1;
        ny = int( dimensions.y() / size.value() ) + 1;
        nz = int( dimensions.z() / size.value() ) + 1;
        
        BOOST_ASSERT( nx > 0 and ny > 0 and nz > 0 );
        
        if (nx > 64)
        {
            qDebug() << "Limiting patching in x dimension to a maximum of 64:" << nx;
            nx = 64;
        }

        if (ny > 64)
        {
            qDebug() << "Limiting patching in y dimension to a maximum of 64:" << ny;
            ny = 64;
        }

        if (nz > 64)
        {
            qDebug() << "Limiting patching in z dimension to a maximum of 64:" << nz;
            nz = 64;
        }
        
        inv_gridvec = Vector( nx / dimensions.x(),
                              ny / dimensions.y(),
                              nz / dimensions.z() );
                              
        orgn = center - 0.5*dimensions;
    }
    else
    {
        //we will limit ourselves ot a 16*16*16 box (that is 4096 patches!)
        nx = 16;
        ny = 16;
        nz = 16;
        
        inv_gridvec = Vector( 1 / patch_size.value() );
        
        orgn = center - Vector(8 * patch_size.value());
    }
}

/** Construct for the passed space, placing the center of the grid at "center".
    This tries to divide space using a patch size of 8 A */
BoxPatching::BoxPatching(const Space &space, const Vector &center)
            : ConcreteProperty<BoxPatching,Patching>(space),
              patch_size(0), orgn(0), inv_gridvec(0), nx(0), ny(0), nz(0)
{
    this->operator=( BoxPatching(space, 8*angstrom, center) );
}

/** Construct for the passed space, using the passed patch size. This will try 
    to build a cubic grid of patches where the grid dimension is approximately
    'patch_size', with the center of the grid at (0,0,0) */
BoxPatching::BoxPatching(const Space &space, Length patch_size)
            : ConcreteProperty<BoxPatching,Patching>(space),
              patch_size(0), orgn(0), inv_gridvec(0), nx(0), ny(0), nz(0)
{
    this->operator=( BoxPatching(space, patch_size, Vector(0)) );
}

/** Construct for the passed space - this tries to divide the space
    using a patch size of 8 A, placing the center of the grid at (0,0,0) */
BoxPatching::BoxPatching(const Space &space)
            : ConcreteProperty<BoxPatching,Patching>()
{
    this->operator=( BoxPatching(space, 8*angstrom, Vector(0)) );
}

/** Copy constructor */
BoxPatching::BoxPatching(const BoxPatching &other)
            : ConcreteProperty<BoxPatching,Patching>(other),
              patch_size(other.patch_size), orgn(other.orgn),
              inv_gridvec(other.inv_gridvec),
              nx(other.nx), ny(other.ny), nz(other.nz)
{}

/** Destructor */
BoxPatching::~BoxPatching()
{}

const char* BoxPatching::typeName()
{
    return QMetaType::typeName( qMetaTypeId<BoxPatching>() );
}

/** Copy assignment operator */
BoxPatching& BoxPatching::operator=(const BoxPatching &other)
{
    if (this != &other)
    {
        patch_size = other.patch_size;
        orgn = other.orgn;
        inv_gridvec = other.inv_gridvec;
        nx = other.nx;
        ny = other.ny;
        nz = other.nz;
        Patching::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool BoxPatching::operator==(const BoxPatching &other) const
{
    return patch_size == other.patch_size and
           orgn == other.orgn and
           inv_gridvec == other.inv_gridvec and
           nx == other.nx and ny == other.ny and nz == other.nz and
           Patching::operator==(other);
}

/** Comparison operator */
bool BoxPatching::operator!=(const BoxPatching &other) const
{
    return not BoxPatching::operator==(other);
}

/** Return the patch size used to control the rough size of the patches */
Length BoxPatching::patchSize() const
{
    return patch_size;
}

/** Return the center of the patching grid */
Vector BoxPatching::center() const
{
    if (nx == 0 or ny == 0 or nz == 0)
        return Vector(0);

    else
        return orgn + 0.5*Vector( nx / inv_gridvec.x(),
                                  ny / inv_gridvec.y(),
                                  nz / inv_gridvec.z() );
}

/** Return the dimensions of each path (the lengths of each side of the box) */
Vector BoxPatching::patchDimension() const
{
    if (nx == 0 or ny == 0 or nz == 0)
        return Vector(0);
    
    else
        return Vector( 1.0 / inv_gridvec.x(),
                       1.0 / inv_gridvec.y(),
                       1.0 / inv_gridvec.z() );
}

/** Return a string representation of the patching */
QString BoxPatching::toString() const
{
    if (nx == 0 or ny == 0 or nz == 0)
        return QObject::tr("BoxPatching::null");
    else
        return QObject::tr("BoxPatching( nPatches()=[%1,%2,%3]=%4, center()=%5, "
                                        "patchDimension()=%6 )")
                    .arg(nx).arg(ny).arg(nz).arg(nx*ny*nz+1)
                    .arg(center().toString()).arg(patchDimension().toString());
}

/** Return the number of patches */
int BoxPatching::nPatches() const
{
    //there is one patch for each box in the grid, plus one
    //patch for all points that are outside the grid
    return nx*ny*nz + 1;
}

/** Return the index of the passed point - this assumes that the
    point has already been mapped into the correct space */
int BoxPatching::getIndex(const Vector &point) const
{
    //get the (i,j,k) indicies of the patch that contains the point
    //(each will space from -(n?/2) to (n?/2)-1

    //translate the point so that it is relative to the origin of the grid
    Vector p = (point - orgn);
    
    p = Vector( p.x() * inv_gridvec.x(),
                p.y() * inv_gridvec.y(),
                p.z() * inv_gridvec.z() );

    int i = int(p.x());
    int j = int(p.y());
    int k = int(p.z());
    
    if (i < 0 or i >= nx or j < 0 or j >= ny or k < 0 or k >= nz)
        //this point is outside of the grid
        return nx*ny*nz;
    else
        return i + j*nx + k*nx*ny;
}

/** Return the AABox that completely encloses the ith patch box */
AABox BoxPatching::patchBox(int ith) const
{
    if (ith < 0 or ith >= nx*ny*nz)
        return AABox(Vector(0), Vector(1.0e50));
        
    else
    {
        //decompose this back into i,j,k
        int k = ith / (nx*ny);
        ith -= k*nx*ny;

        int j = ith / nx;
        int i = ith - j*nx;
        
        return AABox( Vector( orgn.x() + (i+0.5) / inv_gridvec.x(),
                              orgn.y() + (j+0.5) / inv_gridvec.y(),
                              orgn.z() + (k+0.5) / inv_gridvec.z() ),
                      Vector( 0.5 / inv_gridvec.x(),
                              0.5 / inv_gridvec.y(),
                              0.5 / inv_gridvec.z() ) );
    }
}

/** Return the index of the patch that contains the passed point */
int BoxPatching::patchIndex(const Vector &point) const
{
    if ( space().isPeriodic() )
        return BoxPatching::getIndex( space().getMinimumImage(point, Vector(0)) );
    else
        return BoxPatching::getIndex( point );
}

/** Return the AABox of that completely encloses the patch that contains the 
    point 'point' */
AABox BoxPatching::patchBox(const Vector &point) const
{
    return patchBox( patchIndex(point) );
}

/** Return the index of the patch that contains the passed point, together
    with the point mapped into that patch (if the space is periodic) */
QPair<int,Vector> BoxPatching::patchIndexAndCenter(const Vector &point) const
{
    if ( space().isPeriodic() )
    {
        Vector mapped_point = space().getMinimumImage(point, Vector(0));
        
        return QPair<int,Vector>( BoxPatching::getIndex(mapped_point), mapped_point );
    }
    else
        return QPair<int,Vector>( BoxPatching::getIndex(point), point );
}

/** Recreate the patching for the passed space */
PatchingPtr BoxPatching::repatch(const Space &new_space) const
{
    return BoxPatching(new_space, patch_size, center());
}
