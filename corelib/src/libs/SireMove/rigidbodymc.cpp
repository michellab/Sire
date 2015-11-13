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

#include "rigidbodymc.h"
#include "uniformsampler.h"
#include "ensemble.h"

#include "SireSystem/system.h"

#include "SireMol/partialmolecule.h"
#include "SireMol/molecule.h"
#include "SireMol/moleculedata.h"
#include "SireMol/atomcoords.h"
#include "SireMol/atomelements.h"
#include "SireMol/atommasses.h"
#include "SireMol/atomselection.h"

#include "SireVol/space.h"

#include "SireMaths/quaternion.h"
#include "SireMaths/vectorproperty.h"
#include "SireMaths/sphere.h"

#include "SireUnits/units.h"
#include "SireUnits/temperature.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>
#include <QTime>
#include <QElapsedTimer>

using namespace SireMove;
using namespace SireSystem;
using namespace SireMol;
using namespace SireUnits;
using namespace SireVol;
using namespace SireStream;
using namespace SireMaths;

static const RegisterMetaType<RigidBodyMC> r_rbmc;

/** Serialise to a binary datastream */
QDataStream SIREMOVE_EXPORT &operator<<(QDataStream &ds,
                                        const RigidBodyMC &rbmc)
{
    writeHeader(ds, r_rbmc, 7);

    SharedDataStream sds(ds);

    sds << rbmc.smplr << rbmc.center_function
        << rbmc.adel << rbmc.rdel
        << rbmc.reflect_radius
        << rbmc.reflect_points
        << rbmc.reflect_moves
        << rbmc.sync_trans << rbmc.sync_rot << rbmc.common_center;
    
    sds << quint32( rbmc.mol_reflectors.count() );
    
    for (QHash< MolNum,QPair<Vector,double> >::const_iterator it = rbmc.mol_reflectors.constBegin();
         it != rbmc.mol_reflectors.constEnd();
         ++it)
    {
        sds << it.key() << it.value().first << it.value().second;
    }
    
    sds << static_cast<const MonteCarlo&>(rbmc);

    return ds;
}

/** Deserialise from a binary datastream */
QDataStream SIREMOVE_EXPORT &operator>>(QDataStream &ds, RigidBodyMC &rbmc)
{
    VersionID v = readHeader(ds, r_rbmc);

    rbmc = RigidBodyMC();

    rbmc.center_function = GetCOGPoint();
    rbmc.reflect_points.clear();
    rbmc.sync_trans = false;
    rbmc.sync_rot = false;
    rbmc.common_center = false;
    rbmc.reflect_moves = false;

    if (v == 7)
    {
        SharedDataStream sds(ds);

        sds >> rbmc.smplr >> rbmc.center_function
            >> rbmc.adel >> rbmc.rdel
            >> rbmc.reflect_radius
            >> rbmc.reflect_points
            >> rbmc.reflect_moves
            >> rbmc.sync_trans >> rbmc.sync_rot
            >> rbmc.common_center;
        
        quint32 nreflect;
        
        sds >> nreflect;
        
        if (nreflect > 0)
            rbmc.mol_reflectors.reserve(nreflect);
        else
            rbmc.mol_reflectors.clear();
        
        for (int i=0; i<nreflect; ++i)
        {
            MolNum molnum;
            Vector center;
            double radius;
            
            sds >> molnum >> center >> radius;
            
            rbmc.mol_reflectors.insert(molnum, QPair<Vector,double>(center,radius));
        }
        
        sds >> static_cast<MonteCarlo&>(rbmc);
    }
    else if (v == 6)
    {
        SharedDataStream sds(ds);

        Vector reflect_center;

        sds >> rbmc.smplr >> rbmc.center_function
            >> rbmc.adel >> rbmc.rdel
            >> reflect_center >> rbmc.reflect_radius
            >> rbmc.reflect_moves
            >> rbmc.sync_trans >> rbmc.sync_rot
            >> rbmc.common_center;
        
        if (rbmc.reflect_moves)
            rbmc.reflect_points.append(reflect_center);
        
        quint32 nreflect;
        
        sds >> nreflect;
        
        if (nreflect > 0)
            rbmc.mol_reflectors.reserve(nreflect);
        else
            rbmc.mol_reflectors.clear();
        
        for (int i=0; i<nreflect; ++i)
        {
            MolNum molnum;
            Vector center;
            double radius;
            
            sds >> molnum >> center >> radius;
            
            rbmc.mol_reflectors.insert(molnum, QPair<Vector,double>(center,radius));
        }
        
        sds >> static_cast<MonteCarlo&>(rbmc);
    }
    else if (v == 5)
    {
        SharedDataStream sds(ds);

        Vector reflect_center;

        sds >> rbmc.smplr >> rbmc.center_function
            >> rbmc.adel >> rbmc.rdel
            >> reflect_center >> rbmc.reflect_radius
            >> rbmc.reflect_moves
            >> rbmc.sync_trans >> rbmc.sync_rot
            >> rbmc.common_center
            >> static_cast<MonteCarlo&>(rbmc);

        if (rbmc.reflect_moves)
            rbmc.reflect_points.append(reflect_center);
    }
    else if (v == 4)
    {
        SharedDataStream sds(ds);
        
        sds >> rbmc.smplr >> rbmc.center_function
            >> rbmc.adel >> rbmc.rdel
            >> rbmc.sync_trans >> rbmc.sync_rot
            >> rbmc.common_center
            >> static_cast<MonteCarlo&>(rbmc);
    }
    else if (v == 3)
    {
        SharedDataStream sds(ds);
        
        sds >> rbmc.smplr
            >> rbmc.adel >> rbmc.rdel
            >> rbmc.sync_trans >> rbmc.sync_rot
            >> rbmc.common_center
            >> static_cast<MonteCarlo&>(rbmc);
    }
    else if (v == 2)
    {
        SharedDataStream sds(ds);
        
        sds >> rbmc.smplr
            >> rbmc.adel >> rbmc.rdel
            >> rbmc.sync_trans >> rbmc.sync_rot
            >> static_cast<MonteCarlo&>(rbmc);
            
        rbmc.common_center = false;
    }
    else if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> rbmc.smplr
            >> rbmc.adel >> rbmc.rdel
            >> static_cast<MonteCarlo&>(rbmc);
    }
    else
        throw version_error(v, "1-5", r_rbmc, CODELOC);

    return ds;
}

/** Null constructor */
RigidBodyMC::RigidBodyMC(const PropertyMap &map) 
            : ConcreteProperty<RigidBodyMC,MonteCarlo>(map),
              center_function( GetCOGPoint() ),
              adel( 0.15 * angstrom ), rdel( 15 * degrees ),
              reflect_radius(0), reflect_moves(false),
              sync_trans(false), sync_rot(false), common_center(false)
{
    MonteCarlo::setEnsemble( Ensemble::NVT(25*celsius) );
}

/** Construct a move that moves molecules returned by the sampler 'sampler' */
RigidBodyMC::RigidBodyMC(const Sampler &sampler, const PropertyMap &map)
            : ConcreteProperty<RigidBodyMC,MonteCarlo>(map),
              smplr(sampler), center_function( GetCOGPoint() ),
              adel( 0.15 * angstrom ),
              rdel( 15 * degrees ), 
              reflect_radius(0), reflect_moves(false),
              sync_trans(false), sync_rot(false), common_center(false)
{
    MonteCarlo::setEnsemble( Ensemble::NVT(25*celsius) );
    smplr.edit().setGenerator( this->generator() );
}

/** Construct a move that moves molecule views from the molecule group 'molgroup',
    selecting views randomly, with each view having an equal chance of
    being chosen */
RigidBodyMC::RigidBodyMC(const MoleculeGroup &molgroup, 
                         const PropertyMap &map)
            : ConcreteProperty<RigidBodyMC,MonteCarlo>(map), 
              smplr( UniformSampler(molgroup) ),
              center_function( GetCOGPoint() ),
              adel( 0.15 * angstrom ), rdel( 15 * degrees ),
              reflect_radius(0), reflect_moves(false),
              sync_trans(false), sync_rot(false), common_center(false)
{
    MonteCarlo::setEnsemble( Ensemble::NVT(25*celsius) );
    smplr.edit().setGenerator( this->generator() );
}

/** Copy constructor */
RigidBodyMC::RigidBodyMC(const RigidBodyMC &other)
            : ConcreteProperty<RigidBodyMC,MonteCarlo>(other), 
              smplr(other.smplr), center_function(other.center_function),
              adel(other.adel), rdel(other.rdel),
              reflect_radius(other.reflect_radius),
              reflect_points(other.reflect_points),
              mol_reflectors(other.mol_reflectors),
              reflect_moves(other.reflect_moves),
              sync_trans(other.sync_trans), sync_rot(other.sync_rot),
              common_center(other.common_center)
{}

/** Destructor */
RigidBodyMC::~RigidBodyMC()
{}

void RigidBodyMC::_pvt_setTemperature(const Temperature &temperature)
{
    MonteCarlo::setEnsemble( Ensemble::NVT(temperature) );
}

/** Copy assignment operator */
RigidBodyMC& RigidBodyMC::operator=(const RigidBodyMC &other)
{
    if (this != &other)
    {
        smplr = other.smplr;
        center_function = other.center_function;
        adel = other.adel;
        rdel = other.rdel;
        reflect_radius = other.reflect_radius;
        reflect_points = other.reflect_points;
        mol_reflectors = other.mol_reflectors;
        reflect_moves = other.reflect_moves;
        sync_trans = other.sync_trans;
        sync_rot = other.sync_rot;
        common_center = other.common_center;
        MonteCarlo::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool RigidBodyMC::operator==(const RigidBodyMC &other) const
{
    return smplr == other.smplr and center_function == other.center_function and
           adel == other.adel and rdel == other.rdel and
           reflect_radius == other.reflect_radius and
           reflect_points == other.reflect_points and
           mol_reflectors == other.mol_reflectors and
           reflect_moves == other.reflect_moves and 
           sync_trans == other.sync_trans and sync_rot == other.sync_rot and
           common_center == other.common_center and
           MonteCarlo::operator==(other);
}

/** Comparison operator */
bool RigidBodyMC::operator!=(const RigidBodyMC &other) const
{
    return not this->operator==(other);
}

/** Return a string representation of this move */
QString RigidBodyMC::toString() const
{
    return QObject::tr("RigidBodyMC( maximumTranslation() = %1 A, "
                       "maximumRotation() = %2 degrees "
                       "nAccepted() = %3 nRejected() = %4 )")
                  .arg(this->maximumTranslation().to(angstrom))
                  .arg(this->maximumRotation().to(degrees))
                  .arg(this->nAccepted())
                  .arg(this->nRejected());
}

/** Set the maximum delta for any translation */
void RigidBodyMC::setMaximumTranslation(Dimension::Length max_translation)
{
    adel = max_translation;
}

/** Set the maximum delta for any rotation */
void RigidBodyMC::setMaximumRotation(Dimension::Angle max_rotation)
{
    rdel = max_rotation;
}

/** Set the function used to get the center of rotation for each molecule */
void RigidBodyMC::setCenterOfRotation(const GetPoint &func)
{
    center_function = func;
}

/** Return the maximum translation for each move */
Dimension::Length RigidBodyMC::maximumTranslation() const
{
    return Dimension::Length(adel);
}

/** Return the maximum rotation for each move */
Dimension::Angle RigidBodyMC::maximumRotation() const
{
    return rdel;
}

/** Return the function used to get the center of rotation of each molecule */
const GetPoint& RigidBodyMC::centerOfRotation() const
{
    return center_function.read();
}

/** Set the sampler (and contained molecule group) that provides
    the random molecules to be moved. This gives the sampler the
    same random number generator that is used by this move */
void RigidBodyMC::setSampler(const Sampler &sampler)
{
    smplr = sampler;
    smplr.edit().setGenerator( this->generator() );
}

/** Set the sampler to be one that selects views at random 
    from the molecule group 'molgroup' (each view has an
    equal chance of being chosen) */
void RigidBodyMC::setSampler(const MoleculeGroup &molgroup)
{
    smplr = UniformSampler(molgroup);
}

/** Return the sampler that is used to draw random molecules */
const Sampler& RigidBodyMC::sampler() const
{
    return smplr;
}

/** Return the molecule group from which random molecule views 
    are drawn */
const MoleculeGroup& RigidBodyMC::moleculeGroup() const
{
    return smplr->group();
}

/** Set the random number generator used by this move */
void RigidBodyMC::setGenerator(const RanGenerator &rangenerator)
{
    MonteCarlo::setGenerator(rangenerator);
    smplr.edit().setGenerator(this->generator());
}

/** Completely switch off use of the reflection sphere or volume */
void RigidBodyMC::disableReflectionVolume()
{
    reflect_moves = false;
    reflect_points.clear();
    reflect_radius = 0;
    mol_reflectors.clear();
}

/** Completely switch off use of the reflection sphere or volume */
void RigidBodyMC::disableReflectionSphere()
{
    disableReflectionVolume();
}

/** Set whether or not to synchronise translation of all of the views */
void RigidBodyMC::setSynchronisedTranslation(bool on)
{
    sync_trans = on;
    
    if (sync_trans)
        disableReflectionVolume();
}

/** Set whether or not to synchronise rotation of all of the views */
void RigidBodyMC::setSynchronisedRotation(bool on)
{
    sync_rot = on;
    
    if (sync_rot)
        disableReflectionVolume();
}

/** Set whether or not to use the same rotation center for all
    synchronised molecules */
void RigidBodyMC::setSharedRotationCenter(bool on)
{
    common_center = on;
}

/** Turn on rigid body move reflections. This makes sure that only
    molecules within the specified sphere can be moved, and any moves
    are reflected so that these molecules will always remain within 
    the sphere */
void RigidBodyMC::setReflectionSphere(Vector sphere_center,
                                      SireUnits::Dimension::Length sphere_radius)
{
    if (sphere_radius.value() < 0.01)
    {
        reflect_moves = false;
        reflect_points.clear();
        reflect_radius = 0;
    }
    else
    {
        reflect_moves = true;
        reflect_radius = sphere_radius.value();
        reflect_points.clear();
        reflect_points.append(sphere_center);
        
        sync_trans = false;
        sync_rot = false;
    }
}

/** Return whether or not these moves use a reflection sphere */
bool RigidBodyMC::usesReflectionMoves() const
{
    return reflect_moves;
}

/** Return the center of the reflection sphere. Returns a zero vector
    if a reflection sphere is not being used */
Vector RigidBodyMC::reflectionSphereCenter() const
{
    if (reflect_points.count() == 1)
        return reflect_points.at(0);
    else
        return Vector(0);
}

/** Return the radius of the reflection sphere. This returns zero
    if the reflection sphere is not being used */
SireUnits::Dimension::Length RigidBodyMC::reflectionSphereRadius() const
{
    return SireUnits::Dimension::Length( reflect_radius );
}

/** Turn on rigid body move reflections for molecule 'molnum'. This makes sure
    this molecule is moved only within the sphere centered at 'sphere_center'
    with radius 'sphere_radius'. Any moves are reflected so that 
    this molecule will always remain within the sphere */
void RigidBodyMC::setReflectionSphere(MolNum molnum, Vector sphere_center,
                                      SireUnits::Dimension::Length sphere_radius)
{
    if (molnum.isNull())
        return;

    if (sphere_radius.value() < 0)
        sphere_radius = SireUnits::Dimension::Length(0);

    mol_reflectors.insert( molnum, QPair<Vector,double>(sphere_center,sphere_radius.value()) );
    
    sync_trans = false;
    sync_rot = false;
}

/** Turn on rigid body move reflections for molecule 'mol'. This makes sure
    this molecule is moved only within the sphere centered at 'sphere_center'
    with radius 'sphere_radius'. Any moves are reflected so that 
    this molecule will always remain within the sphere */
void RigidBodyMC::setReflectionSphere(const MoleculeView &mol, Vector sphere_center,
                                      SireUnits::Dimension::Length sphere_radius)
{
    this->setReflectionSphere(mol.data().number(), sphere_center, sphere_radius);
}

/** Return whether or not these moves use a reflection sphere for molecule molnum */
bool RigidBodyMC::usesReflectionMoves(MolNum molnum) const
{
    return reflect_moves or mol_reflectors.contains(molnum);
}

/** Return whether or not these moves use a reflection sphere for molecule 'mol' */
bool RigidBodyMC::usesReflectionMoves(const MoleculeView &mol) const
{
    return this->usesReflectionMoves(mol.data().number());
}

/** Return the center of the reflection sphere for molecule 'molnum'. Returns a zero vector
    if a reflection sphere is not being used */
Vector RigidBodyMC::reflectionSphereCenter(MolNum molnum) const
{
    if (mol_reflectors.contains(molnum))
        return mol_reflectors.value(molnum).first;
    else
        return reflectionSphereCenter();
}

/** Return the center of the reflection sphere for molecule 'mol'. Returns a null vector
    if a reflection sphere is not being used */
Vector RigidBodyMC::reflectionSphereCenter(const MoleculeView &mol) const
{
    return this->reflectionSphereCenter(mol.data().number());
}

/** Return the radius of the reflection sphere for molecule 'molnum'. This returns zero
    if the reflection sphere is not being used */
SireUnits::Dimension::Length RigidBodyMC::reflectionSphereRadius(MolNum molnum) const
{
    if (mol_reflectors.contains(molnum))
        return SireUnits::Dimension::Length( mol_reflectors.value(molnum).second );
    else
        return reflectionSphereRadius();
}

/** Return the radius of the reflection sphere for molecule 'mol'. This returns zero
    if the reflection sphere is not being used */
SireUnits::Dimension::Length RigidBodyMC::reflectionSphereRadius(const MoleculeView &mol) const
{
    return this->reflectionSphereRadius(mol.data().number());
}

/** Turn on and specify a reflection volume. This is like a reflection sphere,
    except now the reflection volume is formed as the intersection of the spheres
    whose centers are in 'points', all with radii 'radius'. This replaces any
    reflection sphere set (although not molecule-specific reflection spheres) */
void RigidBodyMC::setReflectionVolume(const QVector<Vector> &points,
                                      SireUnits::Dimension::Length radius)
{
    if (radius.value() < 0.01 or points.isEmpty())
    {
        reflect_moves = false;
        reflect_points.clear();
        reflect_radius = 0;
    }
    else
    {
        reflect_moves = true;
        reflect_points = points;
        reflect_radius = radius.value();
        
        sync_trans = false;
        sync_rot = false;
    }
}

/** Turn on and specify the reflection volume by using the coordinates of 
    all of the atoms in the passed molecule view (excluding light atoms
    if "heavy_atoms_only" is true) */
void RigidBodyMC::setReflectionVolume(const MoleculeView &mol,
                                      SireUnits::Dimension::Length radius,
                                      bool heavy_atoms_only,
                                      const PropertyMap &map)
{
    //get all of the coordinates
    QVector<Vector> points;
    
    const AtomCoords &coords = mol.data().property( map["coordinates"] ).asA<AtomCoords>();
    
    if (heavy_atoms_only)
    {
        if (mol.data().hasProperty( map["element"] ) and
            mol.data().property( map["element"] ).isA<AtomElements>())
        {
            const AtomElements &elems = mol.data().property( map["element"] ).asA<AtomElements>();
            
            if (mol.selectedAll())
            {
                for (int i=0; i<coords.nCutGroups(); ++i)
                {
                    CGIdx ci(i);
                
                    const Vector *coords_array = coords.constData(ci);
                    const Element *elems_array = elems.constData(ci);
                
                    for (int j=0; j<coords.nAtoms(ci); ++j)
                    {
                        if (elems_array[j].nProtons() > 5)
                            points.append(coords_array[j]);
                    }
                }
            }
            else
            {
                AtomSelection selected_atoms = mol.selection();
                
                foreach (CGIdx ci, selected_atoms.selectedCutGroups())
                {
                    const Vector *coords_array = coords.constData(ci);
                    const Element *elems_array = elems.constData(ci);
                    
                    foreach (Index j, selected_atoms.selectedAtoms(ci))
                    {
                        if (elems_array[j].nProtons() > 5)
                            points.append(coords_array[j]);
                    }
                }
            }
        }
        else
        {
            const AtomMasses &masses = mol.data().property( map["mass"] ).asA<AtomMasses>();
            
            if (mol.selectedAll())
            {
                for (int i=0; i<coords.nCutGroups(); ++i)
                {
                    CGIdx ci(i);
                
                    const Vector *coords_array = coords.constData(ci);
                    const SireUnits::Dimension::MolarMass *masses_array = masses.constData(ci);
                
                    for (int j=0; j<coords.nAtoms(ci); ++j)
                    {
                        if (masses_array[j] > Element(5).mass())
                            points.append(coords_array[j]);
                    }
                }
            }
            else
            {
                AtomSelection selected_atoms = mol.selection();
                
                foreach (CGIdx ci, selected_atoms.selectedCutGroups())
                {
                    const Vector *coords_array = coords.constData(ci);
                    const SireUnits::Dimension::MolarMass *masses_array = masses.constData(ci);
                    
                    foreach (Index j, selected_atoms.selectedAtoms(ci))
                    {
                        if (masses_array[j] > Element(5).mass())
                            points.append(coords_array[j]);
                    }
                }
            }
        }
    }
    else
    {
        if (mol.selectedAll())
        {
            for (int i=0; i<coords.nCutGroups(); ++i)
            {
                CGIdx ci(i);
            
                const Vector *coords_array = coords.constData(ci);
            
                for (int j=0; j<coords.nAtoms(ci); ++j)
                {
                    points.append(coords_array[j]);
                }
            }
        }
        else
        {
            AtomSelection selected_atoms = mol.selection();
            
            foreach (CGIdx ci, selected_atoms.selectedCutGroups())
            {
                const Vector *coords_array = coords.constData(ci);
                
                foreach (Index j, selected_atoms.selectedAtoms(ci))
                {
                    points.append(coords_array[j]);
                }
            }
        }
    }

    this->setReflectionVolume(points, radius);
}

/** Turn on and specify the reflection volume by using the coordinates of 
    all of the heavy atoms in the passed molecule view */
void RigidBodyMC::setReflectionVolume(const MoleculeView &mol,
                                      SireUnits::Dimension::Length radius,
                                      const PropertyMap &map)
{
    this->setReflectionVolume(mol, radius, true, map);
}

/** Return whether or not rigid body moves are confined to a reflection volume.
    Note that this also returns true if we are using a reflection sphere
    (as a sphere is just a special case of a reflection volume) */
bool RigidBodyMC::usesReflectionVolume() const
{
    return not reflect_points.isEmpty();
}

/** Return the points used to define the reflection volume. This will
    be a single point if a reflection sphere is used. If a reflection
    sphere or volume is not used then this returns an empty list */
QVector<Vector> RigidBodyMC::reflectionVolumePoints() const
{
    return reflect_points;
}

/** Return the set of spheres that are used to define the reflection volume.
    This returns an empty list if the reflection volume / sphere is not used */
QVector<Sphere> RigidBodyMC::reflectionVolume() const
{
    QVector<Sphere> spheres;
    
    for (int i=0; i<reflect_points.count(); ++i)
    {
        spheres.append( Sphere(reflect_points.at(i), reflect_radius) );
    }
    
    return spheres;
}

/** Return the volume of space occupied by the reflection volume */
double RigidBodyMC::reflectedVolume() const
{
    if (reflect_points.isEmpty())
        return 0;
    else
        return Sphere::combinedVolume( this->reflectionVolume() );
}

/** Return the reflection volume radius (same as the reflection sphere radius) */
SireUnits::Dimension::Length RigidBodyMC::reflectionVolumeRadius() const
{
    if (usesReflectionVolume())
        return SireUnits::Dimension::Length(reflect_radius);
    else
        return SireUnits::Dimension::Length(0);
}

/** Return whether or not translation of all molecules is synchronised */
bool RigidBodyMC::synchronisedTranslation() const
{
    return sync_trans;
}

/** Return whether or not rotation of all molecules is synchronised */
bool RigidBodyMC::synchronisedRotation() const
{
    return sync_rot;
}

/** Return whether or not synchronised rotation uses the same
    center of rotation for all molecules */
bool RigidBodyMC::sharedRotationCenter() const
{
    return common_center;
}

Vector getIntersectionPointWithSphere(const Vector &origin, const Vector &direction,
                                      const Vector &sphere_center, const double sphere_radius,
                                      bool *ok)
{
    *ok = true;

    //first, find the intersection of the delta vector with the sphere.
    // The delta vector has origin at O, direction D. The sphere
    // is at origin C, with radius R
    const Vector D = direction.normalise();
    const Vector O = origin;
    const Vector C = sphere_center;
    const double R2 = sphere_radius*sphere_radius;

    //a point P is on the surface of the sphere if (P-C).(P-C) = R^2
    //this means that the intersection of the vector with the sphere
    //must satisfy ( (O + xD) - C ).( (O + xD) - C ) = R^2
    // This gives;
    // (D.D) x^2 + 2 ( O-C ).D x + (O-C).(O-C) - R^2 = 0
    // which is A x^2 + B x + C = 0, where
    //
    // A = D.D
    // B = 2 (O-C).D
    // C = (O-C).(O-C) - R^2
    //
    // which can be solved using the quadratic formula
    //
    // roots = [-B - sqrt(B^2 - 4AC)] / 2A
    //       = [-B + sqrt(B^2 - 4AC)] / 2A
    //
    // To avoid numerical instability, we use;
    //
    // roots = Q / A and C / Q where
    //
    // Q = [-B + sqrt(B^2 - 4AC)] / 2   if B < 0
    // Q = [-B - sqrt(B^2 - 4AC)] / 2   otherwise

    double QA = Vector::dot(D,D);
    double QB = 2.0 * Vector::dot( O-C, D );
    double QC = Vector::dot(O-C, O-C) - R2;

    double B2_minus_4AC = QB*QB - 4*QA*QC;

    if (B2_minus_4AC < 0)
    {
        //the move does not intersect with the sphere... weird...
        *ok = false;
        qDebug() << "WEIRD: VECTOR DOES NOT INTERSECT WITH SPHERE";
        return Vector(0);
    }

    double Q;

    if (QB < 0)
    {
        Q = (-QB - std::sqrt(B2_minus_4AC)) * 0.5;
    }
    else
    {
        Q = (-QB + std::sqrt(B2_minus_4AC)) * 0.5;
    }

    double x0 = Q / QA;
    double x1 = QC / Q;

    if (x0 > x1){ qSwap(x0,x1); }

    if (x1 < 0)
    {
        //the intersection is behind us...
        qDebug() << "Intersection behind us..." << x1;
        *ok = false;
        return Vector(0);
    }

    double x = x0;

    if (x0 < 0){ x = x1; }

    //the intersection point, X, is O + xD
    Vector X = O + x*D;
    
    return X;
}

PartialMolecule reflectMolecule(const PartialMolecule &oldmol, PartialMolecule newmol,
                                const Vector reflect_cent, const double reflect_rad,
                                Vector old_center, Vector new_center,
                                bool has_center_property,
                                const PropertyName &center_property,
                                const GetPoint &center_function,
                                const PropertyMap &map,
                                bool *ok)
{
    *ok = true;

    double dist = (new_center - reflect_cent).length();

    if (dist > reflect_rad)
    {
        Vector delta = new_center - old_center;
    
        //this would move the molecule out of the sphere. We need
        //to work out where the molecule would intersect the surface
        //of the sphere, and then reflect the molecule back inside

        //first, find the intersection of the delta vector with the sphere.
        // The delta vector has origin at O, direction D. The sphere
        // is at origin C, with radius R
        Vector D = delta.normalise();
        Vector O = old_center;
        Vector C = reflect_cent;
        double R2 = reflect_rad*reflect_rad;

        //a point P is on the surface of the sphere if (P-C).(P-C) = R^2
        //this means that the intersection of the vector with the sphere
        //must satisfy ( (O + xD) - C ).( (O + xD) - C ) = R^2
        // This gives;
        // (D.D) x^2 + 2 ( O-C ).D x + (O-C).(O-C) - R^2 = 0
        // which is A x^2 + B x + C = 0, where
        //
        // A = D.D
        // B = 2 (O-C).D
        // C = (O-C).(O-C) - R^2
        //
        // which can be solved using the quadratic formula
        //
        // roots = [-B - sqrt(B^2 - 4AC)] / 2A
        //       = [-B + sqrt(B^2 - 4AC)] / 2A
        //
        // To avoid numerical instability, we use;
        //
        // roots = Q / A and C / Q where
        //
        // Q = [-B + sqrt(B^2 - 4AC)] / 2   if B < 0
        // Q = [-B - sqrt(B^2 - 4AC)] / 2   otherwise

        double QA = Vector::dot(D,D);
        double QB = 2.0 * Vector::dot( O-C, D );
        double QC = Vector::dot(O-C, O-C) - R2;

        double B2_minus_4AC = QB*QB - 4*QA*QC;

        if (B2_minus_4AC < 0)
        {
            //the move does not intersect with the sphere... weird...
            *ok = false;
            qDebug() << "WEIRD: MOVE DOES NOT INTERSECT WITH SPHERE" << B2_minus_4AC;
            return oldmol;
        }

        double Q;

        if (QB < 0)
        {
            Q = (-QB - std::sqrt(B2_minus_4AC)) * 0.5;
        }
        else
        {
            Q = (-QB + std::sqrt(B2_minus_4AC)) * 0.5;
        }

        double x0 = Q / QA;
        double x1 = QC / Q;

        if (x0 > x1){ qSwap(x0,x1); }

        if (x1 < 0)
        {
            //the intersection is behind us...
            qDebug() << "Intersection behind us..." << x1;
            *ok = false;
            return oldmol;
        }

        double x = x0;

        if (x0 < 0){ x = x1; }

        //the intersection point, X, is O + xD
        Vector X = O + x*D;

        //qDebug() << "Reflect at " << X.toString() << (X-C).length() << std::sqrt(R2);

        //ok - now we have the intersection point, the next step is to 
        //get the normal (N) to the sphere at this point, as this defines the
        //reflection plane
        Vector N = (X - C).normalise();

        //We want to reflect the unit vector that intersects with the 
        //sphere about this normal
        Vector X1 = X - D;
        Vector X2 = X1 + 2 * ( (X - Vector::dot(D,N)*N) - X1 );

        //X2 is the reflected vector. Now work out how much we have
        //moved along X1, and then move that same amount along X2
        double dist_x1 = (X-O).length();
        double dist_x2 = delta.length() - dist_x1;

        if (dist_x2 < 0)
        {
            qDebug() << "WEIRD. NEGATIVE REFLECTION DISTANCE??? " 
                     << dist_x2;
            
            *ok = false;
            return oldmol;
        }

        //work out where the new center of the molecule should lie
        Vector new_new_center = X + dist_x2*((X2-X).normalise());

        //now translate the molecule to the new position
        newmol = newmol.move().translate(new_new_center-new_center,map).commit();

        if (has_center_property)
        {
            new_center = newmol.property(center_property).asA<VectorProperty>();
        }
        else
        {
            new_center = center_function(newmol,map);
        }

        double dist = (new_center - reflect_cent).length();
        
        int check_count = 0;
        
        while (dist > reflect_rad)
        {
            if ( dist - reflect_rad > 0.2 )
            {
                qDebug() << "MOVED MOLECULE OUTSIDE SPHERE" << dist << reflect_rad;
                qDebug()
                    << "FIXING THE PROBLEM (MOSTLY CAUSED BY NUMERICAL IMPRECISION)";
            }
            
            //this will be due to a little numerical imprecision
            newmol = newmol.move().translate( 
                    (1.01*(dist-reflect_rad))
                        * ((reflect_cent-new_center).normalise()) ).commit();
            
            if (has_center_property)
            {
                new_center = newmol.property(center_property).asA<VectorProperty>();
            }
            else
            {
                new_center = center_function(newmol,map);
            }
            
            dist = (new_center - reflect_cent).length();
            
            check_count += 1;
            
            if (check_count > 10)
            {
                qDebug() << "WARNING: SOMETHING WEIRD GOING ON."
                         << "SKIPPING THIS MOVE.";
                         
                *ok = false;
                return oldmol;
            }
        }
    }
    
    return newmol;
}

bool inVolume(const Vector &c, const QVector<Vector> &points, const double radius)
{
    for (int i=0; i<points.count(); ++i)
    {
        if (Vector::distance( c, points.constData()[i] ) <= radius)
            return true;
    }
    
    return false;
}

/** This internal function is used to actually move the molecule(s) */
void RigidBodyMC::performMove(System &system,
                              double &old_bias, double &new_bias,
                              const PropertyMap &map)
{
    const PropertyName &center_property = map["center"];

    //update the sampler with the latest version of the molecules
    smplr.edit().updateFrom(system);

    //get the random amounts by which to translate and
    //rotate the molecule(s)
    Vector delta = generator().vectorOnSphere(adel);

    Quaternion rotdelta( rdel * generator().rand(),
                         generator().vectorOnSphere() );

    old_bias = 1;
    new_bias = 1;

    if (reflect_moves and (sync_trans or sync_rot))
        throw SireError::incomplete_code( QObject::tr(
                "Sire does not yet support using the reflection sphere together with "
                "synchronised translation or rotation of molecules."), CODELOC );

    if ( not (sync_trans or sync_rot) )
    {
        //randomly select a molecule to move
        tuple<PartialMolecule,double> mol_and_bias = smplr.read().sample();

        const PartialMolecule &oldmol = mol_and_bias.get<0>();
        
        if (smplr.read().isBiased())
            old_bias = mol_and_bias.get<1>();
        
        if (oldmol.isEmpty())
        {
            qDebug() << "Sampler returned an empty molecule in RigidBodyMC" << this->toString()
                     << this->moleculeGroup().toString()
                     << this->moleculeGroup().nMolecules() << smplr.read().toString();
            return;
        }
        
        const bool has_center_property = (oldmol.selectedAll() and
                                          oldmol.hasProperty(center_property));

        //move the molecule
        PartialMolecule newmol;

        // only use a reflection volume if there are more than one reflection points
        // and the molecule does not have its own reflector
        if (reflect_moves and reflect_points.count() > 1
                    and not mol_reflectors.contains(oldmol.number()))
        {
            if (reflect_radius == 0)
            {
                qDebug() << "CANNOT MOVE MOLECULE AS RESTRICT VOLUME RADIUS IS ZERO"
                         << oldmol.number().toString();
                return;
            }
        
            const double reflect_rad = reflect_radius;
        
            //make sure that the molecule starts from within this volume
            Vector old_center;
            
            if (has_center_property)
            {
                old_center = oldmol.property(center_property).asA<VectorProperty>();
            }
            else
            {
                old_center = center_function.read()(oldmol,map);
            }

            if (not ::inVolume(old_center,reflect_points,reflect_rad))
            {
                //the molecule is already outside the volume, so cannot be moved
                qDebug() << "HOW IS THE MOLECULE OUTSIDE THE VOLUME?"
                         << oldmol.number().toString();
                return;
            }

            int nattempts = 0;
        
            while (true)
            {
                Vector new_center;

                if (has_center_property)
                {
                    newmol = oldmol.move()
                                   .rotate(rotdelta,
                                           oldmol.property(center_property).asA<VectorProperty>(),
                                           map)
                                   .translate(delta, map)
                                   .commit();

                    new_center = newmol.property(center_property).asA<VectorProperty>();
                }
                else
                {
                    newmol = oldmol.move()
                                   .rotate(rotdelta,
                                           center_function.read()(oldmol,map),
                                           map)
                                   .translate(delta, map)
                                   .commit();

                    new_center = center_function.read()(newmol,map);
                }
                
                if (::inVolume(new_center,reflect_points,reflect_rad))
                {
                    //we have successfully moved the molecule :-)
                    break;
                }
                else
                {
                    //try to use the reflection sphere move to reflect the move back
                    //into the sphere
                    
                    //work out which spheres contained the molecule before the move
                    QVarLengthArray<int> spheres_with_point;
                    
                    for (int i=0; i<reflect_points.count(); ++i)
                    {
                        if (Vector::distance(reflect_points.constData()[i], old_center)
                                                                                <= reflect_rad)
                        {
                            spheres_with_point.append(i);
                        }
                    }
                    
                    if (spheres_with_point.isEmpty())
                        throw SireError::program_bug( QObject::tr(
                                "How can the molecule not be in a sphere???"), CODELOC );
                    
                    //now find the point of intersection between the move and each sphere
                    QVarLengthArray<Vector> intersect_points;
                    
                    for (int i=0; i<spheres_with_point.count(); ++i)
                    {
                        const Vector reflect_cent = reflect_points.at( spheres_with_point[i] );
                        
                        bool ok;
                        
                        Vector intersect = ::getIntersectionPointWithSphere(
                                                old_center, new_center-old_center,
                                                reflect_cent, reflect_rad, &ok);
                        
                        if (ok)
                        {
                            intersect_points.append(intersect);
                        }
                        else
                        {
                            spheres_with_point.remove(i);
                        }
                    }
                    
                    //now, hopefully, only one of these points of intersection is at the
                    //boundary of all of these spheres (if not, we will try again as I
                    //don't want to work out the maths of bouncing off the intersection
                    //point of more than one sphere!)
                    int reflect_sphere_id = -1;
                    
                    bool single_sphere = false;
                    
                    for (int i=0; i<spheres_with_point.count(); ++i)
                    {
                        const Vector intersect = intersect_points[i];
                        bool in_spheres = false;

                        for (int j=0; j<spheres_with_point.count(); ++j)
                        {
                            if (i != j)
                            {
                                Vector cent = reflect_points.constData()[spheres_with_point[j]];
                                
                                if (Vector::distance(cent,intersect) <= reflect_rad)
                                {
                                    in_spheres = true;
                                    break;
                                }
                            }
                        }

                        if (not in_spheres)
                        {
                            //this intersection point is not inside any of the other spheres
                            //This must mean that it is the point at which we bounce off of
                            //the boundary formed by all of the spheres.
                            if (reflect_sphere_id != -1)
                            {
                                //oh dear - we are bouncing off the boundary of more
                                //then one sphere!
                                qDebug() << "Bouncing off the boundary of more than one sphere";
                                single_sphere = false;
                                break;
                            }
                            else
                            {
                                reflect_sphere_id = spheres_with_point[i];
                                single_sphere = true;
                            }
                        }
                    }
                    
                    bool ok = false;

                    if (single_sphere)
                    {
                        //now we know from which sphere we bounced, work out the new
                        //position
                        newmol = ::reflectMolecule(oldmol, newmol,
                                                   reflect_points.at(reflect_sphere_id),
                                                   reflect_rad,
                                                   old_center, new_center, has_center_property,
                                                   center_property, center_function, map, &ok);
                    }
                
                    if (ok)
                    {
                        //now check that this new point is inside at least one of the spheres
                        //(if not, try again as I don't want to work out the maths of
                        //double richochets!)
                        if (has_center_property)
                        {
                            new_center = newmol.property(center_property).asA<VectorProperty>();
                        }
                        else
                        {
                            new_center = center_function.read()(newmol,map);
                        }
                        
                        if (::inVolume(new_center,reflect_points,reflect_rad))
                        {
                            //yes - everything is ok
                            break;
                        }
                    }
                    
                    nattempts += 1;

                    //randomly generate new amounts by which to translate and
                    //rotate the molecule
                    delta = generator().vectorOnSphere(adel);

                    rotdelta = Quaternion( rdel * generator().rand(),
                                           generator().vectorOnSphere() );
                    
                    if (nattempts > 50)
                    {
                        qDebug() << "Cannot move molecule as the number of volume reflection "
                                    "attempts has exceeded 50." << oldmol.number().toString();
                        
                        return;
                    }
                }
            }
        }
        else
        {
            if (has_center_property)
            {
                newmol = oldmol.move()
                               .rotate(rotdelta,
                                       oldmol.property(center_property).asA<VectorProperty>(),
                                       map)
                               .translate(delta, map)
                               .commit();
            }
            else
            {
                newmol = oldmol.move()
                               .rotate(rotdelta,
                                       center_function.read()(oldmol,map),
                                       map)
                               .translate(delta, map)
                               .commit();
            }

            //if we are reflecting moves in a sphere, then check that this move
            //won't take us out of the sphere.
            if (reflect_moves or mol_reflectors.contains(oldmol.number()))
            {
                //moves are constrained into a sphere of radius reflect_radius
                //around reflect_center. If the center of geometry moves outside
                //the sphere, then the molecule will bounce off the edge of 
                //the sphere and back into the sphere volume

                Vector reflect_cent;
                double reflect_rad = 0;
                
                if (mol_reflectors.contains(oldmol.number()))
                {
                    reflect_cent = mol_reflectors.value(oldmol.number()).first;
                    reflect_rad = mol_reflectors.value(oldmol.number()).second;
                }
                else if (reflect_points.count() == 1)
                {
                    reflect_cent = reflect_points.at(0);
                    reflect_rad = reflect_radius;
                }
                else
                    throw SireError::program_bug( QObject::tr(
                            "Should not be possible to get here? %1 %2")
                                .arg(reflect_moves).arg( Sire::toString(reflect_points) ),
                                    CODELOC );

                Vector old_center;
                
                if (has_center_property)
                {
                    old_center = oldmol.property(center_property).asA<VectorProperty>();
                }
                else
                {
                    old_center = center_function.read()(oldmol,map);
                }

                if ( (old_center-reflect_cent).length() > reflect_rad )
                {
                    //the molecule is already outside the sphere, so cannot be moved
                    qDebug() << "HOW IS THE MOLECULE OUTSIDE THE SPHERE?";
                    qDebug() << (old_center-reflect_cent).length() << reflect_rad;
                    return;
                }

                Vector new_center;
                
                if (has_center_property)
                {
                    new_center = newmol.property(center_property).asA<VectorProperty>();
                }
                else
                {
                    new_center = center_function.read()(newmol,map);
                }

                bool ok;

                newmol = ::reflectMolecule(oldmol, newmol, reflect_cent, reflect_rad,
                                           old_center, new_center, has_center_property,
                                           center_property, center_function, map, &ok);
            
                if (not ok)
                {
                    qDebug() << "Something went wrong with the reflection move...";
                    return;
                }
            }
        }

        //update the system with the new coordinates
        if (MonteCarlo::usingOptimisedMoves())
            system.update(newmol, false);
        else
            system.update(newmol, true);

        //get the new bias on this molecule
        if (smplr.read().isBiased())
            new_bias = smplr.read().probabilityOf(newmol);
    }
    else if (sync_trans)
    {
        if (sync_rot)
        {
            //translate and rotate all molecules
            const Molecules &molecules = smplr.read().group().molecules();

            Molecules new_molecules = molecules;

            if (common_center)
            {
                AABox box;
            
                //rotate all molecules around the same center
                for (Molecules::const_iterator it = molecules.constBegin();
                     it != molecules.constEnd();
                     ++it)
                {
                    if (it->selectedAll() and it->molecule().hasProperty(center_property))
                    {
                        box += it->molecule().property(center_property).asA<VectorProperty>();
                    }
                    else
                    {
                        box += center_function.read()(*it,map);
                    }
                }
                
                for (Molecules::const_iterator it = molecules.constBegin();
                     it != molecules.constEnd();
                     ++it)
                {
                    PartialMolecule newmol = it->move()
                                                .rotate(rotdelta, box.center(), map)
                                                .translate(delta, map)
                                                .commit();
                    
                    new_molecules.update(newmol);
                }            
                
            }
            else
            {
                for (Molecules::const_iterator it = molecules.constBegin();
                     it != molecules.constEnd();
                     ++it)
                {
                    PartialMolecule newmol = it->move()
                                                .rotate(rotdelta,
                                                        center_function.read()(*it,map),
                                                        map)
                                                .translate(delta, map)
                                                .commit();
                    
                    new_molecules.update(newmol);
                }            
            }
            
            system.update(new_molecules);            
        }
        else
        {
            //translate all molecules
            const Molecules &molecules = smplr.read().group().molecules();

            Molecules new_molecules = molecules;

            for (Molecules::const_iterator it = molecules.constBegin();
                 it != molecules.constEnd();
                 ++it)
            {
                PartialMolecule newmol = it->move()
                                            .translate(delta, map)
                                            .commit();

                new_molecules.update(newmol);
            }

            system.update(new_molecules);

            //then rotate a single random molecule
            smplr.edit().updateFrom(system);

            tuple<PartialMolecule,double> mol_and_bias = smplr.read().sample();

            const PartialMolecule &oldmol = mol_and_bias.get<0>();
            
            if (smplr.read().isBiased())
                old_bias = mol_and_bias.get<1>();

            PartialMolecule newmol;
            
            if (oldmol.selectedAll() and oldmol.hasProperty(center_property))
            {
                newmol = oldmol.move()
                               .rotate(rotdelta,
                                       oldmol.property(center_property).asA<VectorProperty>(),
                                       map)
                               .commit();
            }
            else
            {
                newmol = oldmol.move()
                               .rotate(rotdelta,
                                       center_function.read()(oldmol,map),
                                       map)
                               .commit();
            }

            //update the system with the new coordinates
            system.update(newmol);

            //get the new bias on this molecule
            if (smplr.read().isBiased())
                new_bias = smplr.read().probabilityOf(newmol);
        }
    }
    else if (sync_rot)
    {
        //rotate all of the molecules
        const Molecules &molecules = smplr.read().group().molecules();

        Molecules new_molecules = molecules;

        bool perform_translation = true;

        if (common_center)
        {
            //we cannot perform a move with synchronised rotation and
            //individual rotation if a common center is used, as it 
            //would be very hard to work out the probability of the 
            //reverse move... So we will instead have a 50/50 chance
            //of performing rotation only or translation only
            perform_translation = generator().randBool();
            
            if (not perform_translation)
            {
                AABox box;
                
                for (Molecules::const_iterator it = molecules.constBegin();
                     it != molecules.constEnd();
                     ++it)
                {
                    if (it->selectedAll() and it->molecule().hasProperty(center_property))
                    {
                        box += it->molecule().property(center_property).asA<VectorProperty>();
                    }
                    else
                    {
                        box += center_function.read()(*it,map);
                    }
                }
            
                for (Molecules::const_iterator it = molecules.constBegin();
                     it != molecules.constEnd();
                     ++it)
                {
                    PartialMolecule newmol = it->move()
                                                .rotate(rotdelta, box.center(), map)
                                                .commit();

                    new_molecules.update(newmol);
                }
            }
        }
        else
        {
            for (Molecules::const_iterator it = molecules.constBegin();
                 it != molecules.constEnd();
                 ++it)
            {
                if (it->selectedAll() and it->molecule().hasProperty(center_property))
                {
                    PartialMolecule newmol = it->move()
                                                .rotate(rotdelta,
                                                        it->molecule().property(center_property)
                                                                .asA<VectorProperty>(),
                                                        map)
                                                .commit();
                    
                    new_molecules.update(newmol);
                }
                else
                {
                    PartialMolecule newmol = it->move()
                                                .rotate(rotdelta,
                                                        center_function.read()(*it,map),
                                                        map)
                                                .commit();
                    
                    new_molecules.update(newmol);
                }
            }
        }

        system.update(new_molecules);

        if (perform_translation)
        {
            //then translate a single random molecule
            smplr.edit().updateFrom(system);

            tuple<PartialMolecule,double> mol_and_bias = smplr.read().sample();

            const PartialMolecule &oldmol = mol_and_bias.get<0>();
            
            if (smplr.read().isBiased())
                old_bias = mol_and_bias.get<1>();

            PartialMolecule newmol = oldmol.move()
                                           .translate(delta, map)
                                           .commit();

            //update the system with the new coordinates
            system.update(newmol);

            //get the new bias on this molecule
            if (smplr.read().isBiased())
                new_bias = smplr.read().probabilityOf(newmol);
        }
    }
}

/** Extract from 'mols' all molecules that are within the reflection sphere / volume
    described by this move. This allows the molecules that will be affected by this
    move to be separated out from the rest. 
    
    If 'buffer' is greater than zero, then this will also extract molecules that are
    within 'buffer' of the reflection volume, and will additionally translate those
    molecules so that they are moved into the volume. The buffer can be used
    to pack more molecules into the volume than would be achieved just with a simple
    selection.
    
    Note that this returns all molecules if a reflection sphere or volume is not used
*/
Molecules RigidBodyMC::extract(const Molecules &mols, SireUnits::Dimension::Length buffer) const
{
    if (reflect_points.isEmpty())
        return mols;
    
    Molecules extracted;

    const PropertyMap &map = propertyMap();
    const PropertyName center_property = map["center"];

    if (buffer.value() <= 0)
    {
        for (Molecules::const_iterator it = mols.constBegin();
             it != mols.constEnd();
             ++it)
        {
            Vector center;
        
            if (it.value().hasProperty(center_property))
            {
                center = it.value().data().property(center_property).asA<VectorProperty>();
            }
            else
            {
                center = center_function.read()(it.value(),map);
            }
            
            //is this molecule within any of the spheres
            bool in_sphere = false;
            
            for (int i=0; i<reflect_points.count(); ++i)
            {
                if (Vector::distance(reflect_points.at(i),center) <= reflect_radius)
                {
                    in_sphere = true;
                    break;
                }
            }
            
            if (in_sphere)
            {
                extracted.add(it.value());
            }
        }
    }
    else
    {
        //need to add on the buffer and move molecules if necessary
        double reflect_plus_buffer = reflect_radius + buffer.value();
        
        for (Molecules::const_iterator it = mols.constBegin();
             it != mols.constEnd();
             ++it)
        {
            Vector center;
        
            if (it.value().hasProperty(center_property))
            {
                center = it.value().data().property(center_property).asA<VectorProperty>();
            }
            else
            {
                center = center_function.read()(it.value(),map);
            }
            
            //is this molecule within any of the spheres
            bool in_buffer = false;
            bool in_sphere = false;
            
            for (int i=0; i<reflect_points.count(); ++i)
            {
                double dist = Vector::distance(reflect_points.at(i),center);
            
                if (dist <= reflect_plus_buffer)
                {
                    in_buffer = true;
                    
                    if (dist <= reflect_radius)
                    {
                        in_sphere = true;
                        break;
                    }
                }
            }
            
            if (in_sphere)
            {
                extracted.add(it.value());
            }
            else if (in_buffer)
            {
                //we need to translate this molecule into the volume - find the closest
                //reflection sphere and translate the molecule along the vector from the
                //sphere center to molecule center
                double shortest_dist = -1;
                int closest_sphere = -1;
                
                for (int i=0; i<reflect_points.count(); ++i)
                {
                    double dist = Vector::distance(reflect_points.at(i),center);
                    
                    if (closest_sphere == -1 or dist < shortest_dist)
                    {
                        shortest_dist = dist;
                        closest_sphere = i;
                    }
                }
                
                //to move the molecule into the sphere, we need to translate by
                // (reflect_radius - shortest_dist - 0.05) * vector(sphere_center -> mol_center)
                // (we use 0.05 so that the molecule is placed 'just' inside the sphere)
                Vector delta = (reflect_radius - shortest_dist - 0.05) *
                                (center - reflect_points.at(closest_sphere)).normalise();
                
                ViewsOfMol mol = it.value().move().translate(delta).commit();
                
                extracted.add(mol);
            }
        }
    }
    
    return extracted;
}

/** Extract from 'mols' all molecules that are within the reflection sphere / volume
    described by this move. This allows the molecules that will be affected by this
    move to be separated out from the rest. 
 
    Note that this returns all molecules if a reflection sphere or volume is not used
*/
Molecules RigidBodyMC::extract(const Molecules &mols) const
{
    return this->extract(mols, SireUnits::Dimension::Length(0));
}

/** Attempt 'n' rigid body moves of the views of the system 'system' */
void RigidBodyMC::move(System &system, int nmoves, bool record_stats)
{
    if (nmoves <= 0)
        return;

    QElapsedTimer t, t2;
    
    qint64 old_ns = 0;
    qint64 copy_ns = 0;
    qint64 nrg_ns = 0;
    qint64 move_ns = 0;
    qint64 test_ns = 0;
    qint64 reject_ns = 0;
    qint64 accept_ns = 0;

    const PropertyMap &map = Move::propertyMap();
    
    if (nmoves > 1)
        t2.start();
    
    for (int i=0; i<nmoves; ++i)
    {
        //get the old total energy of the system
        if (nmoves > 1)
            if (nmoves > 1)
                t.start();
            

        double old_nrg = system.energy( this->energyComponent() );

        if (nmoves > 1)

        //save the old system
        if (nmoves > 1)
            t.start();
        
        System old_system = system;

        if (nmoves > 1)
            copy_ns += t.nsecsElapsed();

        double old_bias = 1;
        double new_bias = 1;

        if (nmoves > 1)
            t.start();

        this->performMove(system, old_bias, new_bias, map);

        if (nmoves > 1)

            if (nmoves > 1)
                move_ns += t.nsecsElapsed();

        //calculate the energy of the system
        if (nmoves > 1)
            if (nmoves > 1)
                t.start();


        double new_nrg = system.energy( this->energyComponent() );

        if (nmoves > 1)

        //accept or reject the move based on the change of energy
        //and the biasing factors
        if (nmoves > 1)
            t.start();

        const bool accept_move = this->test(new_nrg, old_nrg, new_bias, old_bias);
        
        if (nmoves > 1)

            if (nmoves > 1)
                test_ns += t.nsecsElapsed();

        if (accept_move)
        {
            //the move has been rejected. Destroy the old state and accept the move
            if (nmoves > 1)
                t.start();
            
            old_system = System();
            system.accept();
            
            if (nmoves > 1)
                accept_ns += t.nsecsElapsed();
        }
        else
        {
            //the move has been rejected - reset the state
            if (nmoves > 1)
                t.start();
            
            system = old_system;
            
            if (nmoves > 1)
                reject_ns += t.nsecsElapsed();
        }

        if (record_stats)
        {
            system.collectStats();
        }
    }
    
    qint64 ns = t2.nsecsElapsed();
    
    /*if (nmoves > 1)
    {
        qDebug() << "Timing for" << nmoves << "(" << (0.000001*ns) << ")";
        qDebug() << "OLD:" << (0.000001*old_ns) << "COPY:" << (0.000001*copy_ns)
                 << "MOVE:" << (0.000001*move_ns) << "ENERGY:" << (0.000001*nrg_ns)
                 << "TEST:" << (0.000001*test_ns) << "ACCEPT:" << (0.000001*accept_ns)
                 << "REJECT:" << (0.000001*reject_ns);
    }*/
}

const char* RigidBodyMC::typeName()
{
    return QMetaType::typeName( qMetaTypeId<RigidBodyMC>() );
}
