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

#include "sireglobal.h"

#include <algorithm>

#include "internalparameters.h"

#include "SireMol/molecule.h"
#include "SireMol/partialmolecule.h"
#include "SireMol/mover.hpp"
#include "SireMol/cgidx.h"

#include "SireVol/coordgroup.h"

#include "SireFF/errors.h"

#include "SireError/errors.h"

#include "tostring.h"

using namespace SireMM;
using namespace SireMM::detail;

using namespace SireBase;
using namespace SireCAS;
using namespace SireVol;
using namespace SireMol;
using namespace SireFF::detail;

using namespace SireStream;

//////////
////////// Implementation of SireMM::detail::CGIDQuad
//////////

/** Serialise to a binary datastream */
QDataStream& operator<<(QDataStream &ds, const CGIDQuad &quad)
{
    ds << quad.cgIdx0() << quad.cgIdx1()
       << quad.cgIdx2() << quad.cgIdx3();
       
    return ds;
}

/** Extract from a binary datastream */
QDataStream& operator>>(QDataStream &ds, CGIDQuad &quad)
{
    CGIdx cgidx0, cgidx1, cgidx2, cgidx3;
    
    ds >> cgidx0 >> cgidx1 >> cgidx2 >> cgidx3;
    
    quad = CGIDQuad(cgidx0, cgidx1, cgidx2, cgidx3);
    
    return ds;
}

/** Constructor */
CGIDQuad::CGIDQuad()
{}

/** Construct the ID for a single CutGroup */
CGIDQuad::CGIDQuad(CGIdx cg0)
{
    cgidxs[0] = cg0;
    cgidxs[1] = cg0;
    cgidxs[2] = cg0;
    cgidxs[3] = cg0;
}

/** Construct the ID for a pair of CutGroups */
CGIDQuad::CGIDQuad(CGIdx cg0, CGIdx cg1)
{
    if (cg0 <= cg1)
    {
        cgidxs[0] = cg0;
        cgidxs[1] = cg1;
    }
    else
    {
        cgidxs[0] = cg1;
        cgidxs[1] = cg0;
    }

    cgidxs[2] = cgidxs[1];
    cgidxs[3] = cgidxs[1];
}

/** Construct the ID for a triple of CutGroups */
CGIDQuad::CGIDQuad(CGIdx cg0, CGIdx cg1, CGIdx cg2)
{
    cgidxs[0] = cg0;
    cgidxs[1] = cg1;
    cgidxs[2] = cg2;
    
    qSort(cgidxs, cgidxs+3);
    
    cgidxs[3] = cgidxs[2];
    
    if (cgidxs[0] == cgidxs[1])
        cgidxs[1] = cgidxs[2];
}

/** Construct the ID for a quadruplet of CutGroups */
CGIDQuad::CGIDQuad(CGIdx cg0, CGIdx cg1, CGIdx cg2, CGIdx cg3)
{
    cgidxs[0] = cg0;
    cgidxs[1] = cg1;
    cgidxs[2] = cg2;
    cgidxs[3] = cg3;
    
    qSort(cgidxs, cgidxs+4);
    
    if (cgidxs[0] == cgidxs[1])
    {
        cgidxs[1] = cgidxs[2];
        cgidxs[2] = cgidxs[3];
    }

    if (cgidxs[0] == cgidxs[1])
    {
        cgidxs[1] = cgidxs[2];
        cgidxs[2] = cgidxs[3];
    }

    if (cgidxs[1] == cgidxs[2])
        cgidxs[2] = cgidxs[3];
}

/** Copy constructor */
CGIDQuad::CGIDQuad(const CGIDQuad &other)
{
    cgidxs[0] = other.cgidxs[0];
    cgidxs[1] = other.cgidxs[1];
    cgidxs[2] = other.cgidxs[2];
    cgidxs[3] = other.cgidxs[3];
}

/** Destructor */
CGIDQuad::~CGIDQuad()
{}

/** Copy assignment operator */
CGIDQuad& CGIDQuad::operator=(const CGIDQuad &other)
{
    if (this != &other)
    {
        cgidxs[0] = other.cgidxs[0];
        cgidxs[1] = other.cgidxs[1];
        cgidxs[2] = other.cgidxs[2];
        cgidxs[3] = other.cgidxs[3];
    }

    return *this;
}

/** Comparison operator */
bool CGIDQuad::operator==(const CGIDQuad &other) const
{
    return cgidxs[0] == other.cgidxs[0] and cgidxs[1] == other.cgidxs[1] and
           cgidxs[2] == other.cgidxs[2] and cgidxs[3] == other.cgidxs[3];
}

/** Comparison operator */
bool CGIDQuad::operator!=(const CGIDQuad &other) const
{
    return cgidxs[0] != other.cgidxs[0] or cgidxs[1] != other.cgidxs[1] or
           cgidxs[2] != other.cgidxs[2] or cgidxs[3] != other.cgidxs[3];
}

/** Return whether this IDs only a single CutGroup */
bool CGIDQuad::isSingleCutGroup() const
{
    return cgidxs[0] == cgidxs[1];
}

/** Return whether this IDs a pair of CutGroups */
bool CGIDQuad::isDoubleCutGroup() const
{
    return cgidxs[1] == cgidxs[2];
}

/** Return whether this IDs a triple of CutGroups */
bool CGIDQuad::isTripleCutGroup() const
{
    return cgidxs[2] == cgidxs[3];
}

/** Return whether this IDs four distinct CutGroups */
bool CGIDQuad::isQuadrupleCutGroup() const
{
    return cgidxs[2] != cgidxs[3];
}

/** Return the index of the first identified CutGroup */
CGIdx CGIDQuad::cgIdx0() const
{
    return cgidxs[0];
}

/** Return the index of the first identified CutGroup */
CGIdx CGIDQuad::cgIdx1() const
{
    return cgidxs[1];
}

/** Return the index of the first identified CutGroup */
CGIdx CGIDQuad::cgIdx2() const
{
    return cgidxs[2];
}

/** Return the index of the first identified CutGroup */
CGIdx CGIDQuad::cgIdx3() const
{
    return cgidxs[3];
}

/** Return a hash of this ID */
uint CGIDQuad::hash() const
{
    return ( quint32(cgidxs[0]) << 24 ) |
           ( (quint32(cgidxs[1]) << 16) & 0x00FF0000 ) |
           ( (quint32(cgidxs[2]) <<  8) & 0x0000FF00 ) |
           ( quint32(cgidxs[3]) & 0x000000FF );
}

//////////
////////// Implementation of GroupInternalNonPhysParameters
//////////

/** Serialise to a binary datastream */
QDataStream& operator<<(QDataStream &ds,
                        const GroupInternalNonPhysParameters &group)
{
    SharedDataStream sds(ds);
    
    sds << group.improper_params
        << group.improper_theta_forces << group.improper_phi_forces
        << group.ub_params << group.ub_forces;
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream& operator>>(QDataStream &ds,
                        GroupInternalNonPhysParameters &group)
{
    SharedDataStream sds(ds);
    
    sds >> group.improper_params
        >> group.improper_theta_forces >> group.improper_phi_forces
        >> group.ub_params >> group.ub_forces;
        
    return ds;
}

/** Constructor */
GroupInternalNonPhysParameters::GroupInternalNonPhysParameters()
                               : RefCountData()
{}

GroupInternalNonPhysParameters::GroupInternalNonPhysParameters(
                                    const GroupInternalNonPhysParameters &other)
        : RefCountData(),
          improper_params(other.improper_params),
          improper_theta_forces(other.improper_theta_forces),
          improper_phi_forces(other.improper_phi_forces),
          ub_params(other.ub_params), ub_forces(other.ub_forces)
{}

/** Destructor */
GroupInternalNonPhysParameters::~GroupInternalNonPhysParameters()
{}

/** Copy assignment operator */
GroupInternalNonPhysParameters& GroupInternalNonPhysParameters::operator=(
                                        const GroupInternalNonPhysParameters &other)
{
    if (this != &other)
    {
        improper_params = other.improper_params;
        improper_theta_forces = other.improper_theta_forces;
        improper_phi_forces = other.improper_phi_forces;
        
        ub_params = other.ub_params;
        ub_forces = other.ub_forces;
    }
    
    return *this;
}

/** Comparison operator */
bool GroupInternalNonPhysParameters::operator==(
                            const GroupInternalNonPhysParameters &other) const
{
    return improper_params == other.improper_params and
           ub_params == other.ub_params;
}

/** Comparison operator */
bool GroupInternalNonPhysParameters::operator!=(
                            const GroupInternalNonPhysParameters &other) const
{
    return not this->operator==(other);
}

//////////
////////// Implementation of GroupInternalCrossParameters
//////////

/** Serialise to a binary datastream */
QDataStream& operator<<(QDataStream &ds,
                        const GroupInternalCrossParameters &group)
{
    SharedDataStream sds(ds);
    
    sds << group.stretch_stretch_params << group.stretch_stretch_r01_forces
        << group.stretch_stretch_r21_forces

        << group.stretch_bend_params << group.stretch_bend_theta_forces
        << group.stretch_bend_r01_forces << group.stretch_bend_r21_forces
    
        << group.bend_bend_params << group.bend_bend_theta012_forces
        << group.bend_bend_theta213_forces << group.bend_bend_theta310_forces

        << group.stretch_bend_torsion_params
        << group.stretch_bend_torsion_phi_forces
        << group.stretch_bend_torsion_r01_forces
        << group.stretch_bend_torsion_r12_forces
        << group.stretch_bend_torsion_r32_forces
        << group.stretch_bend_torsion_r03_forces
        << group.stretch_bend_torsion_theta012_forces
        << group.stretch_bend_torsion_theta321_forces;

    return ds;
}

/** Extract from a binary datastream */
QDataStream& operator>>(QDataStream &ds,
                        GroupInternalCrossParameters &group)
{
    SharedDataStream sds(ds);
    
    sds >> group.stretch_stretch_params >> group.stretch_stretch_r01_forces
        >> group.stretch_stretch_r21_forces

        >> group.stretch_bend_params >> group.stretch_bend_theta_forces
        >> group.stretch_bend_r01_forces >> group.stretch_bend_r21_forces
    
        >> group.bend_bend_params >> group.bend_bend_theta012_forces
        >> group.bend_bend_theta213_forces >> group.bend_bend_theta310_forces

        >> group.stretch_bend_torsion_params
        >> group.stretch_bend_torsion_phi_forces
        >> group.stretch_bend_torsion_r01_forces
        >> group.stretch_bend_torsion_r12_forces
        >> group.stretch_bend_torsion_r32_forces
        >> group.stretch_bend_torsion_r03_forces
        >> group.stretch_bend_torsion_theta012_forces
        >> group.stretch_bend_torsion_theta321_forces;

    return ds;
}

/** Constructor */
GroupInternalCrossParameters::GroupInternalCrossParameters()
                             : RefCountData()
{}

/** Copy constructor */
GroupInternalCrossParameters::GroupInternalCrossParameters(
                                    const GroupInternalCrossParameters &other)
    : RefCountData(),
      stretch_stretch_params(other.stretch_stretch_params),
      stretch_stretch_r01_forces(other.stretch_stretch_r01_forces),
      stretch_stretch_r21_forces(other.stretch_stretch_r21_forces),
      stretch_bend_params(other.stretch_bend_params),
      stretch_bend_theta_forces(other.stretch_bend_theta_forces),
      stretch_bend_r01_forces(other.stretch_bend_r01_forces),
      stretch_bend_r21_forces(other.stretch_bend_r21_forces),
      bend_bend_params(other.bend_bend_params),
      bend_bend_theta012_forces(other.bend_bend_theta012_forces),
      bend_bend_theta213_forces(other.bend_bend_theta213_forces),
      bend_bend_theta310_forces(other.bend_bend_theta310_forces),
      stretch_bend_torsion_params(other.stretch_bend_torsion_params),
      stretch_bend_torsion_phi_forces(other.stretch_bend_torsion_phi_forces),
      stretch_bend_torsion_r01_forces(other.stretch_bend_torsion_r01_forces),
      stretch_bend_torsion_r12_forces(other.stretch_bend_torsion_r12_forces),
      stretch_bend_torsion_r32_forces(other.stretch_bend_torsion_r32_forces),
      stretch_bend_torsion_r03_forces(other.stretch_bend_torsion_r03_forces),
      stretch_bend_torsion_theta012_forces(other.stretch_bend_torsion_theta012_forces),
      stretch_bend_torsion_theta321_forces(other.stretch_bend_torsion_theta321_forces)
{}

/** Destructor */    
GroupInternalCrossParameters::~GroupInternalCrossParameters()
{}
  
/** Copy assignment operator */  
GroupInternalCrossParameters& GroupInternalCrossParameters::operator=(
                                        const GroupInternalCrossParameters &other)
{
    if (this != &other)
    {
      stretch_stretch_params = other.stretch_stretch_params;
      stretch_stretch_r01_forces = other.stretch_stretch_r01_forces;
      stretch_stretch_r21_forces = other.stretch_stretch_r21_forces;
      stretch_bend_params = other.stretch_bend_params;
      stretch_bend_theta_forces = other.stretch_bend_theta_forces;
      stretch_bend_r01_forces = other.stretch_bend_r01_forces;
      stretch_bend_r21_forces = other.stretch_bend_r21_forces;
      bend_bend_params = other.bend_bend_params;
      bend_bend_theta012_forces = other.bend_bend_theta012_forces;
      bend_bend_theta213_forces = other.bend_bend_theta213_forces;
      bend_bend_theta310_forces = other.bend_bend_theta310_forces;
      stretch_bend_torsion_params = other.stretch_bend_torsion_params;
      stretch_bend_torsion_phi_forces = other.stretch_bend_torsion_phi_forces;
      stretch_bend_torsion_r01_forces = other.stretch_bend_torsion_r01_forces;
      stretch_bend_torsion_r12_forces = other.stretch_bend_torsion_r12_forces;
      stretch_bend_torsion_r32_forces = other.stretch_bend_torsion_r32_forces;
      stretch_bend_torsion_r03_forces = other.stretch_bend_torsion_r03_forces;
      stretch_bend_torsion_theta012_forces = other.stretch_bend_torsion_theta012_forces;
      stretch_bend_torsion_theta321_forces = other.stretch_bend_torsion_theta321_forces;
    }
    
    return *this;
}

/** Comparison operator */
bool GroupInternalCrossParameters::operator==(
                                    const GroupInternalCrossParameters &other) const
{
    return stretch_stretch_params == other.stretch_stretch_params and
           stretch_bend_params == other.stretch_bend_params and
           bend_bend_params == other.bend_bend_params and
           stretch_bend_torsion_params == other.stretch_bend_torsion_params;
}

//////////
////////// Implementation of GroupInternalParametersData
//////////

/** Serialise to a binary datastream */
QDataStream& operator<<(QDataStream &ds,
                        const GroupInternalParametersData &group)
{
    SharedDataStream sds(ds);
    
    sds << group.idquad 
        << group.bond_params << group.bond_forces
        << group.angle_params << group.angle_forces
        << group.dihedral_params << group.dihedral_forces
        << group.nonphys_terms << group.cross_terms;
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream& operator>>(QDataStream &ds,
                        GroupInternalParametersData &group)
{
    SharedDataStream sds(ds);
    
    sds >> group.idquad
        >> group.bond_params >> group.bond_forces
        >> group.angle_params >> group.angle_forces
        >> group.dihedral_params >> group.dihedral_forces
        >> group.nonphys_terms >> group.cross_terms;
        
    return ds;
}

static SharedDataPointer<GroupInternalNonPhysParameters>
                    null_nonphys( new GroupInternalNonPhysParameters() );

static SharedDataPointer<GroupInternalCrossParameters>
                    null_cross( new GroupInternalCrossParameters() );

/** Constructor */
GroupInternalParametersData::GroupInternalParametersData()
                            : RefCountData(),
                              nonphys_terms( null_nonphys ),
                              cross_terms( null_cross )
{}

/** Construct for the specified combination of CutGroups */
GroupInternalParametersData::GroupInternalParametersData(const CGIDQuad &cgids)
                            : RefCountData(),
                              idquad(cgids),
                              nonphys_terms( null_nonphys ),
                              cross_terms( null_cross )
{}                              

/** Copy constructor */
GroupInternalParametersData::GroupInternalParametersData(
                                       const GroupInternalParametersData &other)
      : RefCountData(),
        idquad(other.idquad),
        bond_params(other.bond_params), bond_forces(other.bond_forces),
        angle_params(other.angle_params), angle_forces(other.angle_forces),
        dihedral_params(other.dihedral_params), dihedral_forces(other.dihedral_forces),
        nonphys_terms(other.nonphys_terms),
        cross_terms(other.cross_terms)
{}

/** Destructor */
GroupInternalParametersData::~GroupInternalParametersData()
{}

/** Copy assignment operator */
GroupInternalParametersData& GroupInternalParametersData::operator=(
                                        const GroupInternalParametersData &other)
{
    if (this != &other)
    {
        idquad = other.idquad;
        bond_params = other.bond_params;
        bond_forces = other.bond_forces;
        angle_params = other.angle_params;
        angle_forces = other.angle_forces;
        dihedral_params = other.dihedral_params;
        dihedral_forces = other.dihedral_forces;
        nonphys_terms = other.nonphys_terms;
        cross_terms = other.cross_terms;
    }
    
    return *this;
}

/** Comparison operator */
bool GroupInternalParametersData::operator==(
                                    const GroupInternalParametersData &other) const
{
    return idquad == other.idquad and 
           bond_params == other.bond_params and
           angle_params == other.angle_params and
           dihedral_params == other.dihedral_params and
           
           (nonphys_terms.constData() == other.nonphys_terms.constData() or
            *(nonphys_terms) == *(other.nonphys_terms)) and
            
           (cross_terms.constData() == other.cross_terms.constData() or
            *(cross_terms) == *(other.cross_terms));
}

/** Comparison operator */
bool GroupInternalParametersData::operator!=(
                            const GroupInternalParametersData &other) const
{
    return not this->operator==(other);
}

/** Return whether or not this group has any cross terms */
bool GroupInternalParametersData::hasCrossTerms() const
{
    return cross_terms.constData() != null_cross.constData();
}

/** Return whether or not this group has any non-physical parameters */
bool GroupInternalParametersData::hasNonPhysicalParameters() const
{
    return nonphys_terms.constData() != null_nonphys.constData();
}

//////////
////////// Implementation of GroupInternalParameters
//////////

/** Serialise to a binary datastream */
QDataStream SIREMM_EXPORT  &operator<<(QDataStream &ds,
                                       const GroupInternalParameters &group)
{
    SharedDataStream sds(ds);
    sds << group.d;
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds,
                                      GroupInternalParameters &group)
{
    SharedDataStream sds(ds);
    sds >> group.d;
    
    return ds;
}

static SharedDataPointer<GroupInternalParametersData>
                            null_data( new GroupInternalParametersData() );

/** Constructor */
GroupInternalParameters::GroupInternalParameters()
                        : d(null_data)
{}

/** Construct for the specified CutGroups */
GroupInternalParameters::GroupInternalParameters(const detail::CGIDQuad &cgids)
                        : d( new GroupInternalParametersData(cgids) )
{}

/** Copy constructor */
GroupInternalParameters::GroupInternalParameters(const GroupInternalParameters &other)
                        : d(other.d)
{}

/** Destructor */
GroupInternalParameters::~GroupInternalParameters()
{}

/** Copy assignment operator */
GroupInternalParameters& GroupInternalParameters::operator=(
                                                const GroupInternalParameters &other)
{
    d = other.d;
    
    return *this;
}

/** Comparison operator */
bool GroupInternalParameters::operator==(const GroupInternalParameters &other) const
{
    return d.constData() == other.d.constData() or
           *d == *(other.d);
}

/** Comparison operator */
bool GroupInternalParameters::operator!=(const GroupInternalParameters &other) const
{
    return d.constData() != other.d.constData() and
           *d != *(other.d);
}

/** Return whether or not this has any physical parameters
    (bond, angle or dihedral) */
bool GroupInternalParameters::hasPhysicalParameters() const
{
    return not (d->bond_params.isEmpty() and d->angle_params.isEmpty() and
                d->dihedral_params.isEmpty());
}

/** Return whether or not this has any non-physical parameters
    (Urey-Bradley or improper terms) */
bool GroupInternalParameters::hasNonPhysicalParameters() const
{
    return d->hasNonPhysicalParameters();
}

/** Return whether or not this has any cross terms
    (stretch-stretch, stretch-bend, bend-bend, stretch-bend-torsion) */
bool GroupInternalParameters::hasCrossTerms() const
{
    return d->hasCrossTerms();
}

/** Return whether this group is empty (contains no parameters) */
bool GroupInternalParameters::isEmpty() const
{
    return d.constData() == null_data.constData() or
           not (this->hasPhysicalParameters() or
                this->hasNonPhysicalParameters() or
                this->hasCrossTerms());
}

/** Return whether or not this group contains parameters from only
    a single CutGroup */
bool GroupInternalParameters::isSingleCutGroup() const
{
    return d->idquad.isSingleCutGroup();
}

/** Return whether or not this group contains parameters from only
    two CutGroups */
bool GroupInternalParameters::isDoubleCutGroup() const
{
    return d->idquad.isDoubleCutGroup();
}

/** Return whether or not this group contains parameters from only
    three CutGroups */
bool GroupInternalParameters::isTripleCutGroup() const
{
    return d->idquad.isTripleCutGroup();
}

/** Return whether or not this group contains parameters from
    four CutGroups */
bool GroupInternalParameters::isQuadrupleCutGroup() const
{
    return d->idquad.isQuadrupleCutGroup();
}

/** Return whether or not this group contains parameters from
    only a single CutGroup, with index 'cgidx0' */
bool GroupInternalParameters::isSingleCutGroup(CGIdx cgidx0) const
{
    return d->idquad == CGIDQuad(cgidx0);
}

/** Return whether or not this group contains parameters from
    only two CutGroups, with indicies 'cgidx0' and 'cgidx1' */
bool GroupInternalParameters::isDoubleCutGroup(CGIdx cgidx0, CGIdx cgidx1) const
{
    return d->idquad == CGIDQuad(cgidx0,cgidx1);
}

/** Return whether or not this group contains parameters from
    three CutGroups, with indicies 'cgidx0', 'cgidx1' and 'cgidx2' */
bool GroupInternalParameters::isTripleCutGroup(CGIdx cgidx0, CGIdx cgidx1, 
                                               CGIdx cgidx2) const
{
    return d->idquad == CGIDQuad(cgidx0,cgidx1,cgidx2);
}

/** Return whether or not this group contains parameters from
    four CutGroups, with indicies 'cgidx0', 'cgidx1', 'cgidx2' and 'cgidx3' */
bool GroupInternalParameters::isQuadrupleCutGroup(CGIdx cgidx0, CGIdx cgidx1,
                                                  CGIdx cgidx2, CGIdx cgidx3) const
{
    return d->idquad == CGIDQuad(cgidx0,cgidx1,cgidx2,cgidx3);
}

/** Return whether or not this group contains parameters 
    that involve atoms in the CutGroup at index 'cgidx' */
bool GroupInternalParameters::refersTo(CGIdx cgidx) const
{
    return d->idquad.cgIdx0() == cgidx or
           d->idquad.cgIdx1() == cgidx or
           d->idquad.cgIdx2() == cgidx or
           d->idquad.cgIdx3() == cgidx;
}

/** Return whether or not this group contains parameters
    that involve atoms in any of the CutGroups whose indicies are
    in 'cgidxs' */
bool GroupInternalParameters::refersTo(const QSet<CGIdx> &cgidxs) const
{
    foreach (CGIdx cgidx, cgidxs)
    {
        if (this->refersTo(cgidx))
            return true;
    }
    
    return false;
}

/** Return the index of the first group */
CGIdx GroupInternalParameters::cgIdx0() const
{
    return d->idquad.cgIdx0();
}

/** Return the index of the first group */
CGIdx GroupInternalParameters::cgIdx1() const
{
    return d->idquad.cgIdx1();
}

/** Return the index of the first group */
CGIdx GroupInternalParameters::cgIdx2() const
{
    return d->idquad.cgIdx2();
}

/** Return the index of the first group */
CGIdx GroupInternalParameters::cgIdx3() const
{
    return d->idquad.cgIdx3();
}

/** Return all of the bond potentials for this group */
const QVector<TwoAtomFunction>& GroupInternalParameters::bondPotential() const
{
    return d->bond_params;
}

/** Return the bond force ( -dE/dr ) */
const QVector<TwoAtomFunction>& GroupInternalParameters::bondForces() const
{
    return d->bond_forces;
}

/** Return the angle potentials for this group */
const QVector<ThreeAtomFunction>& GroupInternalParameters::anglePotential() const
{
    return d->angle_params;
}

/** Return the angle force ( -dE/dtheta ) */
const QVector<ThreeAtomFunction>& GroupInternalParameters::angleForces() const
{
    return d->angle_forces;
}

/** Return the dihedral potentials for this group */
const QVector<FourAtomFunction>& GroupInternalParameters::dihedralPotential() const
{
    return d->dihedral_params;
}

/** Return the dihedral force ( -dE/dphi ) */
const QVector<FourAtomFunction>& GroupInternalParameters::dihedralForces() const
{
    return d->dihedral_forces;
}

/** Return the improper potentials for this group */
const QVector<FourAtomFunction>& GroupInternalParameters::improperPotential() const
{
    return d->nonphys_terms->improper_params;
}

/** Return the improper force ( -dE/dtheta ) */
const QVector<FourAtomFunction>& GroupInternalParameters::improper_Theta_Forces() const
{
    return d->nonphys_terms->improper_theta_forces;
}

/** Return the improper force ( -dE/dphi ) */
const QVector<FourAtomFunction>& GroupInternalParameters::improper_Phi_Forces() const
{
    return d->nonphys_terms->improper_phi_forces;
}

/** Return the Urey-Bradley potentials for this group */
const QVector<TwoAtomFunction>& GroupInternalParameters::ureyBradleyPotential() const
{
    return d->nonphys_terms->ub_params;
}

/** Return the Urey-Bradley force ( -dE/dr ) */
const QVector<TwoAtomFunction>& GroupInternalParameters::ureyBradleyForces() const
{
    return d->nonphys_terms->ub_forces;
}

/** Return the stretch-stretch potentials for this group */
const QVector<ThreeAtomFunction>& 
GroupInternalParameters::stretchStretchPotential() const
{
    return d->cross_terms->stretch_stretch_params;
}

/** Return the stretch-stretch force ( -dE/dr_01 ) */
const QVector<ThreeAtomFunction>& 
GroupInternalParameters::stretchStretch_R01_Forces() const
{
    return d->cross_terms->stretch_stretch_r01_forces;
}

/** Return the stretch-stretch force ( -dE/dr_01 ) */
const QVector<ThreeAtomFunction>& 
GroupInternalParameters::stretchStretch_R21_Forces() const
{
    return d->cross_terms->stretch_stretch_r01_forces;
}

/** Return the stretch-bend potentials for this group */
const QVector<ThreeAtomFunction>& 
GroupInternalParameters::stretchBendPotential() const
{
    return d->cross_terms->stretch_bend_params;
}

/** Return the stretch-bend force ( -dE/dtheta ) */
const QVector<ThreeAtomFunction>& 
GroupInternalParameters::stretchBend_Theta_Forces() const
{
    return d->cross_terms->stretch_bend_theta_forces;
}

/** Return the stretch-bend force ( -dE/dtheta ) */
const QVector<ThreeAtomFunction>& 
GroupInternalParameters::stretchBend_R01_Forces() const
{
    return d->cross_terms->stretch_bend_r01_forces;
}

/** Return the stretch-bend force ( -dE/dtheta ) */
const QVector<ThreeAtomFunction>& 
GroupInternalParameters::stretchBend_R21_Forces() const
{
    return d->cross_terms->stretch_bend_r21_forces;
}

/** Return the bend-bend potentials for this group */
const QVector<FourAtomFunction>& 
GroupInternalParameters::bendBendPotential() const
{
    return d->cross_terms->bend_bend_params;
}

/** Return bend-bend force ( -dE/dtheta_012 ) */
const QVector<FourAtomFunction>& 
GroupInternalParameters::bendBend_Theta012_Forces() const
{
    return d->cross_terms->bend_bend_theta012_forces;
}

/** Return bend-bend force ( -dE/dtheta_012 ) */
const QVector<FourAtomFunction>& 
GroupInternalParameters::bendBend_Theta213_Forces() const
{
    return d->cross_terms->bend_bend_theta213_forces;
}

/** Return bend-bend force ( -dE/dtheta_012 ) */
const QVector<FourAtomFunction>& 
GroupInternalParameters::bendBend_Theta310_Forces() const
{
    return d->cross_terms->bend_bend_theta310_forces;
}

/** Return the stretch-bend-torsion potentials for this group */
const QVector<FourAtomFunction>& 
GroupInternalParameters::stretchBendTorsionPotential() const
{
    return d->cross_terms->stretch_bend_torsion_params;
}

/** Return the stretch-bend-torsion force ( -dE/dphi ) */
const QVector<FourAtomFunction>& 
GroupInternalParameters::stretchBendTorsion_Phi_Forces() const
{
    return d->cross_terms->stretch_bend_torsion_phi_forces;
}

/** Return the stretch-bend-torsion force ( -dE/dr_01 ) */
const QVector<FourAtomFunction>& 
GroupInternalParameters::stretchBendTorsion_R01_Forces() const
{
    return d->cross_terms->stretch_bend_torsion_r01_forces;
}

/** Return the stretch-bend-torsion force ( -dE/dr_12 ) */
const QVector<FourAtomFunction>& 
GroupInternalParameters::stretchBendTorsion_R12_Forces() const
{
    return d->cross_terms->stretch_bend_torsion_r12_forces;
}

/** Return the stretch-bend-torsion force ( -dE/dr_32 ) */
const QVector<FourAtomFunction>& 
GroupInternalParameters::stretchBendTorsion_R32_Forces() const
{
    return d->cross_terms->stretch_bend_torsion_r32_forces;
}

/** Return the stretch-bend-torsion force ( -dE/dr_03 ) */
const QVector<FourAtomFunction>& 
GroupInternalParameters::stretchBendTorsion_R03_Forces() const
{
    return d->cross_terms->stretch_bend_torsion_r03_forces;
}

/** Return the stretch-bend-torsion force ( -dE/dtheta_012 ) */
const QVector<FourAtomFunction>& 
GroupInternalParameters::stretchBendTorsion_Theta012_Forces() const
{
    return d->cross_terms->stretch_bend_torsion_theta012_forces;
}

/** Return the stretch-bend-torsion force ( -dE/dtheta_321 ) */
const QVector<FourAtomFunction>& 
GroupInternalParameters::stretchBendTorsion_Theta321_Forces() const
{
    return d->cross_terms->stretch_bend_torsion_theta321_forces;
}

/** Internal function used to set the bond parameters */
void GroupInternalParameters::setBondPotential(
                                    const QVector<TwoAtomFunction> &potential,
                                    const QVector<TwoAtomFunction> &forces)
{
    d->bond_params = potential;
    d->bond_forces = forces;
}
                       
/** Internal function used to set the angle parameters */
void GroupInternalParameters::setAnglePotential(
                                    const QVector<ThreeAtomFunction> &potential,
                                    const QVector<ThreeAtomFunction> &forces)
{
    d->angle_params = potential;
    d->angle_forces = forces;
}

/** Internal function used to set the dihedral parameters */
void GroupInternalParameters::setDihedralPotential(
                                    const QVector<FourAtomFunction> &potential,
                                    const QVector<FourAtomFunction> &forces)
{
    d->dihedral_params = potential;
    d->dihedral_forces = forces;
}

/** Internal function used to set the improper parameters */
void GroupInternalParameters::setImproperPotential(
                                const QVector<FourAtomFunction> &potential,
                                const QVector<FourAtomFunction> &theta_forces,
                                const QVector<FourAtomFunction> &phi_forces)
{
    d->nonphys_terms->improper_params = potential;
    d->nonphys_terms->improper_theta_forces = theta_forces;
    d->nonphys_terms->improper_phi_forces = phi_forces;
}
                           
/** Internal function used to set the Urey-Bradley parameters */
void GroupInternalParameters::setUreyBradleyPotential(
                                const QVector<TwoAtomFunction> &potential,
                                const QVector<TwoAtomFunction> &forces)
{
    d->nonphys_terms->ub_params = potential;
    d->nonphys_terms->ub_forces = forces;
}

/** Internal function used to set the stretch-stretch parameters */
void GroupInternalParameters::setStretchStretchPotential(
                                const QVector<ThreeAtomFunction> &potential,
                                const QVector<ThreeAtomFunction> &r01_forces,
                                const QVector<ThreeAtomFunction> &r21_forces)
{
    d->cross_terms->stretch_stretch_params = potential;
    d->cross_terms->stretch_stretch_r01_forces = r01_forces;
    d->cross_terms->stretch_stretch_r21_forces = r21_forces;
}
                                 
/** Internal function used to set the stretch-bend parameters */
void GroupInternalParameters::setStretchBendPotential(
                                const QVector<ThreeAtomFunction> &potential,
                                const QVector<ThreeAtomFunction> &theta_forces,
                                const QVector<ThreeAtomFunction> &r01_forces,
                                const QVector<ThreeAtomFunction> &r21_forces)
{
    d->cross_terms->stretch_bend_params = potential;
    d->cross_terms->stretch_bend_theta_forces = theta_forces;
    d->cross_terms->stretch_bend_r01_forces = r01_forces;
    d->cross_terms->stretch_bend_r21_forces = r21_forces;
}
                              
/** Internal function used to set the bend-bend parameters */
void GroupInternalParameters::setBendBendPotential(
                                const QVector<FourAtomFunction> &potential,
                                const QVector<FourAtomFunction> &theta012_forces,
                                const QVector<FourAtomFunction> &theta213_forces,
                                const QVector<FourAtomFunction> &theta310_forces)
{
    d->cross_terms->bend_bend_params = potential;
    d->cross_terms->bend_bend_theta012_forces = theta012_forces;
    d->cross_terms->bend_bend_theta213_forces = theta213_forces;
    d->cross_terms->bend_bend_theta310_forces = theta310_forces;
}

/** Internal function used to set the stretch-bend-torsion parameters */
void GroupInternalParameters::setStretchBendTorsionPotential(
                                 const QVector<FourAtomFunction> &potential,
                                 const QVector<FourAtomFunction> &phi_forces,
                                 const QVector<FourAtomFunction> &r01_forces,
                                 const QVector<FourAtomFunction> &r12_forces,
                                 const QVector<FourAtomFunction> &r32_forces,
                                 const QVector<FourAtomFunction> &theta012_forces,
                                 const QVector<FourAtomFunction> &theta321_forces)
{
    d->cross_terms->stretch_bend_torsion_params = potential;
    d->cross_terms->stretch_bend_torsion_phi_forces = phi_forces;
    d->cross_terms->stretch_bend_torsion_r01_forces = r01_forces;
    d->cross_terms->stretch_bend_torsion_r12_forces = r12_forces;
    d->cross_terms->stretch_bend_torsion_r32_forces = r32_forces;
    d->cross_terms->stretch_bend_torsion_theta012_forces = theta012_forces;
    d->cross_terms->stretch_bend_torsion_theta321_forces = theta321_forces;
}

//////////
////////// Implementation of InternalSymbolsBase
//////////

InternalSymbolsBase::InternalSymbolsBase()
{}

InternalSymbolsBase::~InternalSymbolsBase()
{}

//////////
////////// Implementation of BondSymbols
//////////

BondSymbols::BondSymbols() : r_("r")
{
    symbols.insert(r_);
}

BondSymbols::~BondSymbols()
{}
    
/** Return the symbol representing the vector along the bond (r) */
const Symbol& BondSymbols::r() const
{
    return r_;
}

//////////
////////// Implementation of AngleSymbols
//////////

AngleSymbols::AngleSymbols() : theta_("theta")
{
    symbols.insert(theta_);
}

AngleSymbols::~AngleSymbols()
{}
    
/** Return the symbols representing the angle (theta) */
const Symbol& AngleSymbols::theta() const
{
    return theta_;
}

//////////
////////// Implementation of DihedralSymbols
//////////

DihedralSymbols::DihedralSymbols() : phi_("phi")
{
    symbols.insert(phi_);
}

DihedralSymbols::~DihedralSymbols()
{}
    
/** Return the symbol representing the torsion (phi) */
const Symbol& DihedralSymbols::phi() const
{
    return phi_;
}

//////////
////////// Implementation of ImproperSymbols
//////////

ImproperSymbols::ImproperSymbols() : theta_("theta"), phi_("phi")
{
    symbols.insert(theta_);
    symbols.insert(phi_);
}

ImproperSymbols::~ImproperSymbols()
{}
  
/** Return the symbol representing the angle between the improper
    and the plane formed by atoms 1-3 */  
const Symbol& ImproperSymbols::theta() const
{
    return theta_;
}

/** Return the symbol representing the torsion 0-1-2-3 */
const Symbol& ImproperSymbols::phi() const
{
    return phi_;
}

//////////
////////// Implementation of StretchStretchSymbols
//////////

StretchStretchSymbols::StretchStretchSymbols()
                      : r01_("r_{01}"), r21_("r_{21}"), r12_("r_{12}")
{
    symbols.insert(r01_);
    symbols.insert(r21_);
    symbols.insert(r12_);
}

StretchStretchSymbols::~StretchStretchSymbols()
{}
    
/** Return the symbol representing the bond length r_{01} */
const Symbol& StretchStretchSymbols::r01() const
{
    return r01_;
}

/** Return the symbol representing the bond length r_{21} */
const Symbol& StretchStretchSymbols::r21() const
{
    return r21_;
}

/** Return the symbol representing the bond length r_{12} */
const Symbol& StretchStretchSymbols::r12() const
{
    return r12_;
}

//////////
////////// Implementation of StretchBendSymbols
//////////

StretchBendSymbols::StretchBendSymbols()
                   : theta_("theta"), r01_("r_{01}"), r21_("r_{21}"), r12_("r_{12}")
{
    symbols.insert(theta_);
    symbols.insert(r01_);
    symbols.insert(r21_);
    symbols.insert(r12_);
}

StretchBendSymbols::~StretchBendSymbols()
{}

/** Return the symbol representing the angle, theta */  
const Symbol& StretchBendSymbols::theta() const
{
    return theta_;
}

/** Return the symbol representing the bond length, r_{01} */  
const Symbol& StretchBendSymbols::r01() const
{
    return r01_;
}

/** Return the symbol representing the bond length, r_{21} */  
const Symbol& StretchBendSymbols::r21() const
{
    return r21_;
}

/** Return the symbol representing the bond length r_{12} */
const Symbol& StretchBendSymbols::r12() const
{
    return r12_;
}

//////////
////////// Implementation of BendBendSymbols
//////////

BendBendSymbols::BendBendSymbols()
                : theta012_("theta_{012}"), theta213_("theta_{213}"),
                  theta310_("theta_{310}")
{
    symbols.insert(theta012_);
    symbols.insert(theta213_);
    symbols.insert(theta310_);
}

BendBendSymbols::~BendBendSymbols()
{}

/** Return the symbol representing the angle between atoms 0-1-2, theta_{012} */
const Symbol& BendBendSymbols::theta012() const
{
    return theta012_;
}

/** Return the symbol representing the angle between atoms 2-1-3, theta_{213} */
const Symbol& BendBendSymbols::theta213() const
{
    return theta213_;
}

/** Return the symbol representing the angle between atoms 3-1-0, theta_{310} */
const Symbol& BendBendSymbols::theta310() const
{
    return theta310_;
}

//////////
////////// Implementation of StretchBendTorsionSymbols
//////////

StretchBendTorsionSymbols::StretchBendTorsionSymbols()
                          : phi_("phi"), theta012_("theta_{012}"),
                            theta321_("theta_{321}"), r01_("r_{01}"),
                            r12_("r_{12}"), r32_("r_{32}"), r03_("r_{03}")
{
    symbols.insert(phi_);
    symbols.insert(theta012_);
    symbols.insert(theta321_);
    symbols.insert(r01_);
    symbols.insert(r12_);
    symbols.insert(r32_);
    symbols.insert(r03_);
}

StretchBendTorsionSymbols::~StretchBendTorsionSymbols()
{}

/** Return the symbol representing the torsion, phi */
const Symbol& StretchBendTorsionSymbols::phi() const
{
    return phi_;
}

/** Return the symbol representing the angle between atoms 0-1-2, theta_{012} */
const Symbol& StretchBendTorsionSymbols::theta012() const
{
    return theta012_;
}

/** Return the symbol representing the angle between atoms 3-2-1, theta_{321} */
const Symbol& StretchBendTorsionSymbols::theta321() const
{
    return theta321_;
}

/** Return the symbol representing the bond between atoms 0-1, r_{01} */
const Symbol& StretchBendTorsionSymbols::r01() const
{
    return r01_;
}

/** Return the symbol representing the bond between atoms 1-2, r_{12} */
const Symbol& StretchBendTorsionSymbols::r12() const
{
    return r12_;
}

/** Return the symbol representing the bond between atoms 3-2, r_{32} */
const Symbol& StretchBendTorsionSymbols::r32() const
{
    return r32_;
}

/** Return the symbol representing the distance from atom 0 to 3, r_{03} */
const Symbol& StretchBendTorsionSymbols::r03() const
{
    return r03_;
}

//////////
////////// Implementation of InternalSymbols
//////////

InternalSymbols::InternalSymbols()
{
    symbols += bond_;
    symbols += angle_;
    symbols += dihedral_;
    symbols += improper_;
    symbols += ureybradley_;
    symbols += stretchstretch_;
    symbols += stretchbend_;
    symbols += bendbend_;
    symbols += stretchbendtorsion_; 
}

InternalSymbols::~InternalSymbols()
{}

/** Return all of the symbols used in the bond parameters */
const BondSymbols& InternalSymbols::bond() const
{
    return bond_;
}

/** Return all of the symbols used in the angle parameters */
const AngleSymbols& InternalSymbols::angle() const
{
    return angle_;
}

/** Return all of the symbols used in the dihedral parameters */
const DihedralSymbols& InternalSymbols::dihedral() const
{
    return dihedral_;
}

/** Return all of the symbols used in the improper parameters */
const ImproperSymbols& InternalSymbols::improper() const
{
    return improper_;
}

/** Return all of the symbols used in the Urey-Bradley parameters */
const BondSymbols& InternalSymbols::ureyBradley() const
{
    return ureybradley_;
}

/** Return all of the symbols used in the stretch-stretch parameters */
const StretchStretchSymbols& InternalSymbols::stretchStretch() const
{
    return stretchstretch_;
}

/** Return all of the symbols used in the stretch-bend parameters */
const StretchBendSymbols& InternalSymbols::stretchBend() const
{
    return stretchbend_;
}

/** Return all of the symbols used in the bend-bend parameters */
const BendBendSymbols& InternalSymbols::bendBend() const
{
    return bendbend_;
}

/** Return all of the symbols used in the stretch-bend-torsion parameters */
const StretchBendTorsionSymbols& InternalSymbols::stretchBendTorsion() const
{
    return stretchbendtorsion_;
}

//////////
////////// Implementation of InternalParameters
//////////

static const RegisterMetaType<InternalParameters> r_params(NO_ROOT);

/** Serialise to a binary datastream */
QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds,
                                      const InternalParameters &params)
{
    writeHeader(ds, r_params, 1);
    
    SharedDataStream sds(ds);
    
    sds << params.group_params << params.groups_by_cgidx;
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds,
                                      InternalParameters &params)
{
    VersionID v = readHeader(ds, r_params);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> params.group_params >> params.groups_by_cgidx;
        
        params.updateState();
    }
    else
        throw version_error( v, "1", r_params, CODELOC );
        
    return ds;
}

InternalSymbols InternalParameters::function_symbols;

/** Return the index of the group that contains the potential and forces
    for the internals that act only within the CutGroup with index 'cgidx0'.
    
    This returns -1 if there is no such group 
*/
qint32 InternalParameters::getIndex(CGIdx cgidx0) const
{
    QHash< CGIdx, QSet<qint32> >::const_iterator 
                                    it = groups_by_cgidx.constFind(cgidx0);
                                    
    if (it != groups_by_cgidx.constEnd())
    {
        const GroupInternalParameters *group_params_array = group_params.constData();
        
        foreach (qint32 i, it.value())
        {
            if (group_params_array[i].isSingleCutGroup(cgidx0))
                return i;
        }
    }
    
    return -1;
}

/** Return the index of the group that contains the potential and forces
    for the internals that act only between the CutGroups with indicies 
    'cgidx0' and 'cgidx1'.
    
    This returns -1 if there is no such group 
*/
qint32 InternalParameters::getIndex(CGIdx cgidx0, CGIdx cgidx1) const
{
    QHash< CGIdx, QSet<qint32> >::const_iterator 
                                    it0 = groups_by_cgidx.constFind(cgidx0);

    QHash< CGIdx, QSet<qint32> >::const_iterator 
                                    it1 = groups_by_cgidx.constFind(cgidx1);
                                    
    if (it0 != groups_by_cgidx.constEnd() and
        it1 != groups_by_cgidx.constEnd())
    {
        const GroupInternalParameters *group_params_array = group_params.constData();
        
        QSet<qint32> shared_groups = it0.value();
        shared_groups.intersect(it1.value());
        
        foreach (qint32 i, shared_groups)
        {
            if (group_params_array[i].isDoubleCutGroup(cgidx0,cgidx1))
                return i;
        }
    }
    
    return -1;
}

/** Return the index of the group that contains the potential and forces
    for the internals that act only between the three CutGroups with indicies 
    'cgidx0', 'cgidx1' and 'cgidx2'.
    
    This returns -1 if there is no such group 
*/
qint32 InternalParameters::getIndex(CGIdx cgidx0, CGIdx cgidx1,
                                    CGIdx cgidx2) const
{
    QHash< CGIdx, QSet<qint32> >::const_iterator 
                                    it0 = groups_by_cgidx.constFind(cgidx0);

    QHash< CGIdx, QSet<qint32> >::const_iterator 
                                    it1 = groups_by_cgidx.constFind(cgidx1);

    QHash< CGIdx, QSet<qint32> >::const_iterator 
                                    it2 = groups_by_cgidx.constFind(cgidx2);
                                    
    if (it0 != groups_by_cgidx.constEnd() and
        it1 != groups_by_cgidx.constEnd() and
        it2 != groups_by_cgidx.constEnd())
    {
        const GroupInternalParameters *group_params_array = group_params.constData();
        
        QSet<qint32> shared_groups = it0.value();
        shared_groups.intersect(it1.value());
        shared_groups.intersect(it2.value());
        
        foreach (qint32 i, shared_groups)
        {
            if (group_params_array[i].isTripleCutGroup(cgidx0,cgidx1,cgidx2))
                return i;
        }
    }
    
    return -1;
}

/** Return the index of the group that contains the potential and forces
    for the internals that act only between the four CutGroups with indicies 
    'cgidx0', 'cgidx1', 'cgidx2' and 'cgidx3'.
    
    This returns -1 if there is no such group 
*/
qint32 InternalParameters::getIndex(CGIdx cgidx0, CGIdx cgidx1,
                                    CGIdx cgidx2, CGIdx cgidx3) const
{
    QHash< CGIdx, QSet<qint32> >::const_iterator 
                                    it0 = groups_by_cgidx.constFind(cgidx0);

    QHash< CGIdx, QSet<qint32> >::const_iterator 
                                    it1 = groups_by_cgidx.constFind(cgidx1);

    QHash< CGIdx, QSet<qint32> >::const_iterator 
                                    it2 = groups_by_cgidx.constFind(cgidx2);

    QHash< CGIdx, QSet<qint32> >::const_iterator 
                                    it3 = groups_by_cgidx.constFind(cgidx3);
                                    
    if (it0 != groups_by_cgidx.constEnd() and
        it1 != groups_by_cgidx.constEnd() and
        it2 != groups_by_cgidx.constEnd() and
        it3 != groups_by_cgidx.constEnd())
    {
        const GroupInternalParameters *group_params_array = group_params.constData();
        
        QSet<qint32> shared_groups = it0.value();
        shared_groups.intersect(it1.value());
        shared_groups.intersect(it2.value());
        shared_groups.intersect(it3.value());
        
        foreach (qint32 i, shared_groups)
        {
            if (group_params_array[i].isQuadrupleCutGroup(cgidx0,cgidx1,
                                                          cgidx2,cgidx3))
                return i;
        }
    }
    
    return -1;
}

/** Return the index of the group with ID 'idquad'. This returns -1
    if there is no group with this ID */
qint32 InternalParameters::getIndex(const CGIDQuad &idquad) const
{
    if (idquad.isSingleCutGroup())
        return this->getIndex(idquad.cgIdx0());
    
    else if (idquad.isDoubleCutGroup())
        return this->getIndex(idquad.cgIdx0(), idquad.cgIdx1());
        
    else if (idquad.isTripleCutGroup())
        return this->getIndex(idquad.cgIdx0(), idquad.cgIdx1(),
                              idquad.cgIdx2());
                              
    else
        return this->getIndex(idquad.cgIdx0(), idquad.cgIdx1(),
                              idquad.cgIdx2(), idquad.cgIdx3());
}

static const GroupInternalParameters shared_empty_group_params;

/** Return the potential and force parameters for all of the internals
    that involve only atoms in the CutGroup at index 'cgidx0'
    
    This returns an empty set of parameters if there are no internals
    that involve this CutGroup
*/
const GroupInternalParameters& InternalParameters::getGroup(CGIdx cgidx0) const
{
    qint32 i = this->getIndex(cgidx0);
    
    if (i >= 0)
        return group_params.at(i);
    else
        return shared_empty_group_params;
}

/** Return the potential and force parameters for all of the internals
    that involve the two CutGroups at indicies 'cgidx0' and 'cgidx1'
    
    This returns an empty set of parameters if there are no internals
    that involve this CutGroup
*/
const GroupInternalParameters& InternalParameters::getGroup(CGIdx cgidx0, 
                                                            CGIdx cgidx1) const
{
    qint32 i = this->getIndex(cgidx0, cgidx1);
    
    if (i >= 0)
        return group_params.at(i);
    else
        return shared_empty_group_params;
}

/** Return the potential and force parameters for all of the internals
    that involve the three CutGroups at indicies 'cgidx0', 'cgidx1'
    and 'cgidx2'.
    
    This returns an empty set of parameters if there are no internals
    that involve this CutGroup
*/
const GroupInternalParameters& InternalParameters::getGroup(CGIdx cgidx0, 
                                                            CGIdx cgidx1,
                                                            CGIdx cgidx2) const
{
    qint32 i = this->getIndex(cgidx0, cgidx1, cgidx2);
    
    if (i >= 0)
        return group_params.at(i);
    else
        return shared_empty_group_params;
}

/** Return the potential and force parameters for all of the internals
    that involve the four CutGroups at indicies 'cgidx0', 'cgidx1',
    'cgidx2' and 'cgidx3'.
    
    This returns an empty set of parameters if there are no internals
    that involve this CutGroup
*/
const GroupInternalParameters& InternalParameters::getGroup(CGIdx cgidx0, 
                                                            CGIdx cgidx1,
                                                            CGIdx cgidx2, 
                                                            CGIdx cgidx3) const
{
    qint32 i = this->getIndex(cgidx0, cgidx1, cgidx2, cgidx3);
    
    if (i >= 0)
        return group_params.at(i);
    else
        return shared_empty_group_params;
}

/** Add the group that contains the potential and forces for the 
    internals that act only between the atoms between the four
    CutGroups at indicies 'cgidx0', 'cgidx1', 'cgidx2' and 'cgidx3',
    and then return the index of the added group.
*/
qint32 InternalParameters::addGroup(const CGIDQuad &idquad)
{
    qint32 i = this->getIndex(idquad);
    
    if (i < 0)
    {
        //the group doesn't exist - it must be added
        group_params.append( GroupInternalParameters(idquad) );
        
        i = group_params.count() - 1;
        
        if (idquad.isSingleCutGroup())
        {
            groups_by_cgidx[idquad.cgIdx0()].insert(i);
        }
        else if (idquad.isDoubleCutGroup())
        {
            groups_by_cgidx[idquad.cgIdx0()].insert(i);
            groups_by_cgidx[idquad.cgIdx1()].insert(i);
        }
        else if (idquad.isTripleCutGroup())
        {
            groups_by_cgidx[idquad.cgIdx0()].insert(i);
            groups_by_cgidx[idquad.cgIdx1()].insert(i);
            groups_by_cgidx[idquad.cgIdx2()].insert(i);
        }
        else
        {
            groups_by_cgidx[idquad.cgIdx0()].insert(i);
            groups_by_cgidx[idquad.cgIdx1()].insert(i);
            groups_by_cgidx[idquad.cgIdx2()].insert(i);
            groups_by_cgidx[idquad.cgIdx3()].insert(i);
        }
    }
    
    return i;
}

GroupInternalParameters&
InternalParameters::getGroup(const CGIDQuad &idquad, 
                             QHash<CGIDQuad,qint32> &cached_groups)
{
    qint32 i = cached_groups.value(idquad, -1);
    
    if (i < 0)
    {
        i = this->addGroup(idquad);

        cached_groups.insert(idquad, i);
    }
    
    return group_params.data()[i];
}

/** Assert that the symbols in 'test_symbols' are all present in 'have_symbols'.
    This is to make sure that the function of 'test_symbols' is not using
    anything that isn't provided (only symbols in 'have_symbols are provided)
    
    \throw SireError::incompatible_error
*/
void InternalParameters::assertContainsOnly(const QSet<Symbol> &have_symbols,
                                            const QSet<Symbol> &test_symbols) const
{
    QList<Symbol> extra_symbols;

    for (QSet<Symbol>::const_iterator it = test_symbols.constBegin();
         it != test_symbols.constEnd();
         ++it)
    {
        if (not have_symbols.contains(*it))
            extra_symbols.append(*it);
    }
    
    if (not extra_symbols.isEmpty())
        throw SireError::incompatible_error( QObject::tr(
            "This is an incompatible function, as it requires symbols "
            "that are not provided by the internal function, %1")
                .arg( Sire::toString(extra_symbols) ), CODELOC );
}

/** Return the symbols used by the internal functions */
const InternalSymbols& InternalParameters::symbols() const
{
    return function_symbols;
}

/** This adds all of the bond parameters and forces in 'bondparams'
    to the list of parameters */
void InternalParameters::addBonds(const TwoAtomFunctions &bondparams,
                                  QHash<CGIDQuad,qint32> &cached_groups)
{
    if (bondparams.isEmpty())
        return;

    //assert that these are functions only of the bond length
    this->assertContainsOnly(this->symbols().bond(), bondparams.symbols());
    
    //get the potential and forces
    QVector<TwoAtomFunction> potentials = bondparams.potentials();
    QVector<TwoAtomFunction> forces = bondparams.forces(this->symbols().bond().r());

    //sort the internals into groups
    QHash< CGIDQuad, QVector<TwoAtomFunction> > group_potentials;
    QHash< CGIDQuad, QVector<TwoAtomFunction> > group_forces;
    
    const TwoAtomFunction *potentials_array = potentials.constData();
    int n = potentials.count();
    
    for (int i=0; i<n; ++i)
    {
        const TwoAtomFunction &potential = potentials_array[i];
    
        CGIDQuad idquad(potential.atom0().cutGroup(),
                        potential.atom1().cutGroup());
                                         
        group_potentials[idquad].append(potential);
    }

    //do forces separately as not all potentials may have a force
    //(as potential could be constant with r)
    const TwoAtomFunction *forces_array = forces.constData();
    n = forces.count();

    for (int i=0; i<n; ++i)
    {
        const TwoAtomFunction &force = forces_array[i];
        
        CGIDQuad idquad(force.atom0().cutGroup(),
                        force.atom1().cutGroup());
                        
        group_forces[idquad].append(force);
    }
    
    //add the groups into this set
    for (QHash< CGIDQuad,QVector<TwoAtomFunction> >::iterator
                                        it = group_potentials.begin();
         it != group_potentials.end();
         ++it)
    {
        it.value().squeeze();
        
        QVector<TwoAtomFunction> &group_force = group_forces[it.key()];
        group_force.squeeze();
        
        this->getGroup(it.key(), cached_groups)
                        .setBondPotential( it.value(), group_force );
    }
}

/** This adds all of the angle parameters and forces in 'angleparams'
    to the list of parameters */
void InternalParameters::addAngles(const ThreeAtomFunctions &angleparams,
                                   QHash<CGIDQuad,qint32> &cached_groups)
{
    if (angleparams.isEmpty())
        return;

    //assert that these are functions only of the angle (theta)
    this->assertContainsOnly(this->symbols().angle(), angleparams.symbols());
    
    //get the potential and forces
    QVector<ThreeAtomFunction> potentials = angleparams.potentials();
    QVector<ThreeAtomFunction> forces = angleparams.forces(this->symbols()
                                                            .angle().theta());
    
    //sort the internals into groups
    QHash< CGIDQuad, QVector<ThreeAtomFunction> > group_potentials;
    QHash< CGIDQuad, QVector<ThreeAtomFunction> > group_forces;
    
    const ThreeAtomFunction *potentials_array = potentials.constData();
    int n = potentials.count();
    
    for (int i=0; i<n; ++i)
    {
        const ThreeAtomFunction &potential = potentials_array[i];
    
        CGIDQuad idquad(potential.atom0().cutGroup(),
                        potential.atom1().cutGroup(),
                        potential.atom2().cutGroup());
                                         
        group_potentials[idquad].append(potential);
    }

    const ThreeAtomFunction *forces_array = forces.constData();
    n = forces.count();
    
    for (int i=0; i<n; ++i)
    {
        const ThreeAtomFunction &force = forces_array[i];
        
        CGIDQuad idquad(force.atom0().cutGroup(),
                        force.atom1().cutGroup(),
                        force.atom2().cutGroup());
                        
        group_forces[idquad].append(force);
    }
    
    //add the groups into this set
    for (QHash< CGIDQuad,QVector<ThreeAtomFunction> >::iterator
                                        it = group_potentials.begin();
         it != group_potentials.end();
         ++it)
    {
        it.value().squeeze();
        
        QVector<ThreeAtomFunction> &group_force = group_forces[it.key()];
        group_force.squeeze();
        
        this->getGroup(it.key(), cached_groups)
                        .setAnglePotential( it.value(), group_force );
    }
}

/** This adds all of the dihedral parameters and forces in 'dihedralparams'
    to the list of parameters */
void InternalParameters::addDihedrals(const FourAtomFunctions &dihedralparams,
                                      QHash<CGIDQuad,qint32> &cached_groups)
{
    if (dihedralparams.isEmpty())
        return;

    //assert that these are functions only of the torsion (phi)
    this->assertContainsOnly(this->symbols().dihedral(), dihedralparams.symbols());
    
    //get the potential and forces
    QVector<FourAtomFunction> potentials = dihedralparams.potentials();
    QVector<FourAtomFunction> forces = dihedralparams.forces(this->symbols()
                                                            .dihedral().phi());
    
    //sort the internals into groups
    QHash< CGIDQuad, QVector<FourAtomFunction> > group_potentials;
    QHash< CGIDQuad, QVector<FourAtomFunction> > group_forces;
    
    const FourAtomFunction *potentials_array = potentials.constData();
    int n = potentials.count();
    
    for (int i=0; i<n; ++i)
    {
        const FourAtomFunction &potential = potentials_array[i];
    
        CGIDQuad idquad(potential.atom0().cutGroup(),
                        potential.atom1().cutGroup(),
                        potential.atom2().cutGroup(),
                        potential.atom3().cutGroup());
                                         
        group_potentials[idquad].append(potential);
    }

    const FourAtomFunction *forces_array = forces.constData();
    n = forces.count();
    
    for (int i=0; i<n; ++i)
    {
        const FourAtomFunction &force = forces_array[i];
        
        CGIDQuad idquad(force.atom0().cutGroup(),
                        force.atom1().cutGroup(),
                        force.atom2().cutGroup(),
                        force.atom3().cutGroup());
                        
        group_forces[idquad].append(force);                        
    }
    
    //add the groups into this set
    for (QHash< CGIDQuad,QVector<FourAtomFunction> >::iterator
                                        it = group_potentials.begin();
         it != group_potentials.end();
         ++it)
    {
        it.value().squeeze();
        
        QVector<FourAtomFunction> &group_force = group_forces[it.key()];
        group_force.squeeze();
        
        this->getGroup(it.key(), cached_groups)
                        .setDihedralPotential( it.value(), group_force );
    }
}

/** This adds all of the improper parameters and forces in 'improperparams'
    to the list of parameters */
void InternalParameters::addImpropers(const FourAtomFunctions &improperparams,
                                      QHash<CGIDQuad,qint32> &cached_groups)
{
    if (improperparams.isEmpty())
        return;

    //assert that these are functions only of the angle (theta or phi)
    this->assertContainsOnly(this->symbols().improper(), improperparams.symbols());
    
    //get the potential and forces
    QVector<FourAtomFunction> potentials = improperparams.potentials();
    
    QVector<FourAtomFunction> theta_forces = improperparams.forces(this->symbols()
                                                            .improper().theta());
    QVector<FourAtomFunction> phi_forces = improperparams.forces(this->symbols()
                                                            .improper().phi());
    
    //sort the internals into groups
    QHash< CGIDQuad, QVector<FourAtomFunction> > group_potentials;
    
    QHash< CGIDQuad, QVector<FourAtomFunction> > group_theta_forces;
    QHash< CGIDQuad, QVector<FourAtomFunction> > group_phi_forces;
    
    const FourAtomFunction *potentials_array = potentials.constData();
    int n = potentials.count();
    
    for (int i=0; i<n; ++i)
    {
        const FourAtomFunction &potential = potentials_array[i];
    
        CGIDQuad idquad(potential.atom0().cutGroup(),
                        potential.atom1().cutGroup(),
                        potential.atom2().cutGroup(),
                        potential.atom3().cutGroup());
                        
        group_potentials[idquad].append(potential);
    }

    const FourAtomFunction *theta_forces_array = theta_forces.constData();
    n = theta_forces.count();
    
    for (int i=0; i<n; ++i)
    {
        const FourAtomFunction &force = theta_forces_array[i];
        
        CGIDQuad idquad(force.atom0().cutGroup(),
                        force.atom1().cutGroup(),
                        force.atom2().cutGroup(),
                        force.atom3().cutGroup());
                        
        group_theta_forces[idquad].append(force);
    }
    
    const FourAtomFunction *phi_forces_array = phi_forces.constData();
    n = phi_forces.count();
    
    for (int i=0; i<n; ++i)
    {
        const FourAtomFunction &force = phi_forces_array[i];
        
        CGIDQuad idquad(force.atom0().cutGroup(),
                        force.atom1().cutGroup(),
                        force.atom2().cutGroup(),
                        force.atom3().cutGroup());
                        
        group_phi_forces[idquad].append(force);
    }
    
    //add the groups into this set
    for (QHash< CGIDQuad,QVector<FourAtomFunction> >::iterator
                                        it = group_potentials.begin();
         it != group_potentials.end();
         ++it)
    {
        it.value().squeeze();
        
        QVector<FourAtomFunction> &group_theta_force = group_theta_forces[it.key()];
        QVector<FourAtomFunction> &group_phi_force = group_phi_forces[it.key()];
        
        group_theta_force.squeeze();
        group_phi_force.squeeze();
        
        this->getGroup(it.key(), cached_groups)
                        .setImproperPotential( it.value(), 
                                               group_theta_force, group_phi_force );
    }
}

/** This adds all of the bond parameters and forces in 'bondparams'
    to the list of parameters */
void InternalParameters::addUBs(const TwoAtomFunctions &ubparams,
                                QHash<CGIDQuad,qint32> &cached_groups)
{
    if (ubparams.isEmpty())
        return;

    //assert that these are functions only of the bond length
    this->assertContainsOnly(this->symbols().ureyBradley(), ubparams.symbols());
    
    //get the potential and forces
    QVector<TwoAtomFunction> potentials = ubparams.potentials();
    QVector<TwoAtomFunction> forces = ubparams.forces(this->symbols()
                                                    .ureyBradley().r());
    
    //sort the internals into groups
    QHash< CGIDQuad, QVector<TwoAtomFunction> > group_potentials;
    QHash< CGIDQuad, QVector<TwoAtomFunction> > group_forces;
    
    const TwoAtomFunction *potentials_array = potentials.constData();
    int n = potentials.count();
    
    for (int i=0; i<n; ++i)
    {
        const TwoAtomFunction &potential = potentials_array[i];
    
        CGIDQuad idquad(potential.atom0().cutGroup(),
                        potential.atom1().cutGroup());
                                         
        group_potentials[idquad].append(potential);
    }

    //do forces separately as not all potentials may have a force
    //(as potential could be constant with r)
    const TwoAtomFunction *forces_array = forces.constData();
    n = forces.count();

    for (int i=0; i<n; ++i)
    {
        const TwoAtomFunction &force = forces_array[i];
        
        CGIDQuad idquad(force.atom0().cutGroup(),
                        force.atom1().cutGroup());
                        
        group_forces[idquad].append(force);
    }
    
    //add the groups into this set
    for (QHash< CGIDQuad,QVector<TwoAtomFunction> >::iterator
                                        it = group_potentials.begin();
         it != group_potentials.end();
         ++it)
    {
        it.value().squeeze();
        
        QVector<TwoAtomFunction> &group_force = group_forces[it.key()];
        group_force.squeeze();
        
        this->getGroup(it.key(), cached_groups)
                        .setUreyBradleyPotential( it.value(), group_force );
    }
}

/** This adds all of the stretch-stretch parameters and forces in 'ssparams'
    to the list of parameters */
void InternalParameters::addSSs(const ThreeAtomFunctions &ssparams,
                                QHash<CGIDQuad,qint32> &cached_groups)
{
    if (ssparams.isEmpty())
        return;

    //assert that these are functions only of the angle (theta)
    this->assertContainsOnly(this->symbols().stretchStretch(), 
                             ssparams.symbols());
    
    //get the potential and forces
    QVector<ThreeAtomFunction> potentials = ssparams.potentials();

    QVector<ThreeAtomFunction> r01_forces = ssparams.forces(this->symbols()
                                                      .stretchStretch().r01());
    QVector<ThreeAtomFunction> r21_forces = ssparams.forces(this->symbols()
                                                      .stretchStretch().r21());
    
    //sort the internals into groups
    QHash< CGIDQuad, QVector<ThreeAtomFunction> > group_potentials;

    QHash< CGIDQuad, QVector<ThreeAtomFunction> > group_r01_forces;
    QHash< CGIDQuad, QVector<ThreeAtomFunction> > group_r21_forces;
    
    const ThreeAtomFunction *potentials_array = potentials.constData();
    int n = potentials.count();
    
    for (int i=0; i<n; ++i)
    {
        const ThreeAtomFunction &potential = potentials_array[i];
    
        CGIDQuad idquad(potential.atom0().cutGroup(),
                        potential.atom1().cutGroup(),
                        potential.atom2().cutGroup());
                                         
        group_potentials[idquad].append(potential);
    }

    const ThreeAtomFunction *r01_forces_array = r01_forces.constData();
    n = r01_forces.count();
    
    for (int i=0; i<n; ++i)
    {
        const ThreeAtomFunction &force = r01_forces_array[i];
        
        CGIDQuad idquad(force.atom0().cutGroup(),
                        force.atom1().cutGroup(),
                        force.atom2().cutGroup());
                        
        group_r01_forces[idquad].append(force);
    }

    const ThreeAtomFunction *r21_forces_array = r21_forces.constData();
    n = r21_forces.count();
    
    for (int i=0; i<n; ++i)
    {
        const ThreeAtomFunction &force = r21_forces_array[i];
        
        CGIDQuad idquad(force.atom0().cutGroup(),
                        force.atom1().cutGroup(),
                        force.atom2().cutGroup());
                        
        group_r21_forces[idquad].append(force);
    }
    
    //add the groups into this set
    for (QHash< CGIDQuad,QVector<ThreeAtomFunction> >::iterator
                                        it = group_potentials.begin();
         it != group_potentials.end();
         ++it)
    {
        it.value().squeeze();
        
        QVector<ThreeAtomFunction> &group_r01_force = group_r01_forces[it.key()];
        QVector<ThreeAtomFunction> &group_r21_force = group_r21_forces[it.key()];
        
        group_r01_force.squeeze();
        group_r21_force.squeeze();
        
        this->getGroup(it.key(), cached_groups)
                        .setStretchStretchPotential( it.value(), 
                                            group_r01_force, group_r21_force );
    }
}

/** This adds all of the stretch-bend parameters and forces in 'sbparams'
    to the list of parameters */
void InternalParameters::addSBs(const ThreeAtomFunctions &sbparams,
                                QHash<CGIDQuad,qint32> &cached_groups)
{
    if (sbparams.isEmpty())
        return;

    //assert that these are functions only of the angle (theta)
    this->assertContainsOnly(this->symbols().stretchBend(), 
                             sbparams.symbols());
    
    //get the potential and forces
    QVector<ThreeAtomFunction> potentials = sbparams.potentials();

    QVector<ThreeAtomFunction> theta_forces = sbparams.forces(this->symbols()
                                                      .stretchBend().theta());
    QVector<ThreeAtomFunction> r01_forces = sbparams.forces(this->symbols()
                                                      .stretchBend().r01());
    QVector<ThreeAtomFunction> r21_forces = sbparams.forces(this->symbols()
                                                      .stretchBend().r21());
    
    //sort the internals into groups
    QHash< CGIDQuad, QVector<ThreeAtomFunction> > group_potentials;

    QHash< CGIDQuad, QVector<ThreeAtomFunction> > group_theta_forces;
    QHash< CGIDQuad, QVector<ThreeAtomFunction> > group_r01_forces;
    QHash< CGIDQuad, QVector<ThreeAtomFunction> > group_r21_forces;
    
    const ThreeAtomFunction *potentials_array = potentials.constData();
    int n = potentials.count();
    
    for (int i=0; i<n; ++i)
    {
        const ThreeAtomFunction &potential = potentials_array[i];
    
        CGIDQuad idquad(potential.atom0().cutGroup(),
                        potential.atom1().cutGroup(),
                        potential.atom2().cutGroup());
                                         
        group_potentials[idquad].append(potential);
    }

    const ThreeAtomFunction *theta_forces_array = theta_forces.constData();
    n = theta_forces.count();
    
    for (int i=0; i<n; ++i)
    {
        const ThreeAtomFunction &force = theta_forces_array[i];
        
        CGIDQuad idquad(force.atom0().cutGroup(),
                        force.atom1().cutGroup(),
                        force.atom2().cutGroup());
                        
        group_theta_forces[idquad].append(force);
    }

    const ThreeAtomFunction *r01_forces_array = r01_forces.constData();
    n = r01_forces.count();
    
    for (int i=0; i<n; ++i)
    {
        const ThreeAtomFunction &force = r01_forces_array[i];
        
        CGIDQuad idquad(force.atom0().cutGroup(),
                        force.atom1().cutGroup(),
                        force.atom2().cutGroup());
                        
        group_r01_forces[idquad].append(force);
    }

    const ThreeAtomFunction *r21_forces_array = r21_forces.constData();
    n = r21_forces.count();
    
    for (int i=0; i<n; ++i)
    {
        const ThreeAtomFunction &force = r21_forces_array[i];
        
        CGIDQuad idquad(force.atom0().cutGroup(),
                        force.atom1().cutGroup(),
                        force.atom2().cutGroup());
                        
        group_r21_forces[idquad].append(force);
    }
    
    //add the groups into this set
    for (QHash< CGIDQuad,QVector<ThreeAtomFunction> >::iterator
                                        it = group_potentials.begin();
         it != group_potentials.end();
         ++it)
    {
        it.value().squeeze();
        
        QVector<ThreeAtomFunction> &group_theta_force = group_theta_forces[it.key()];
        QVector<ThreeAtomFunction> &group_r01_force = group_r01_forces[it.key()];
        QVector<ThreeAtomFunction> &group_r21_force = group_r21_forces[it.key()];
        
        group_theta_force.squeeze();
        group_r01_force.squeeze();
        group_r21_force.squeeze();
        
        this->getGroup(it.key(), cached_groups)
                        .setStretchBendPotential( it.value(), group_theta_force,
                                            group_r01_force, group_r21_force );
    }
}

/** This adds all of the bend-bend parameters and forces in 'bbparams'
    to the list of parameters */
void InternalParameters::addBBs(const FourAtomFunctions &bbparams,
                                QHash<CGIDQuad,qint32> &cached_groups)
{
    if (bbparams.isEmpty())
        return;

    //assert that these are functions only of the angle (theta)
    this->assertContainsOnly(this->symbols().bendBend(), 
                             bbparams.symbols());
    
    //get the potential and forces
    QVector<FourAtomFunction> potentials = bbparams.potentials();

    QVector<FourAtomFunction> theta012_forces = bbparams.forces(this->symbols()
                                                      .bendBend().theta012());
    QVector<FourAtomFunction> theta213_forces = bbparams.forces(this->symbols()
                                                      .bendBend().theta213());
    QVector<FourAtomFunction> theta310_forces = bbparams.forces(this->symbols()
                                                      .bendBend().theta310());
    
    //sort the internals into groups
    QHash< CGIDQuad, QVector<FourAtomFunction> > group_potentials;

    QHash< CGIDQuad, QVector<FourAtomFunction> > group_theta012_forces;
    QHash< CGIDQuad, QVector<FourAtomFunction> > group_theta213_forces;
    QHash< CGIDQuad, QVector<FourAtomFunction> > group_theta310_forces;
    
    const FourAtomFunction *potentials_array = potentials.constData();
    int n = potentials.count();
    
    for (int i=0; i<n; ++i)
    {
        const FourAtomFunction &potential = potentials_array[i];
    
        CGIDQuad idquad(potential.atom0().cutGroup(),
                        potential.atom1().cutGroup(),
                        potential.atom2().cutGroup(),
                        potential.atom3().cutGroup());
                                         
        group_potentials[idquad].append(potential);
    }

    const FourAtomFunction *theta012_forces_array = theta012_forces.constData();
    n = theta012_forces.count();
    
    for (int i=0; i<n; ++i)
    {
        const FourAtomFunction &force = theta012_forces_array[i];
        
        CGIDQuad idquad(force.atom0().cutGroup(),
                        force.atom1().cutGroup(),
                        force.atom2().cutGroup(),
                        force.atom3().cutGroup());
                        
        group_theta012_forces[idquad].append(force);
    }

    const FourAtomFunction *theta213_forces_array = theta213_forces.constData();
    n = theta213_forces.count();
    
    for (int i=0; i<n; ++i)
    {
        const FourAtomFunction &force = theta213_forces_array[i];
        
        CGIDQuad idquad(force.atom0().cutGroup(),
                        force.atom1().cutGroup(),
                        force.atom2().cutGroup(),
                        force.atom3().cutGroup());
                        
        group_theta213_forces[idquad].append(force);
    }

    const FourAtomFunction *theta310_forces_array = theta310_forces.constData();
    n = theta310_forces.count();
    
    for (int i=0; i<n; ++i)
    {
        const FourAtomFunction &force = theta310_forces_array[i];
        
        CGIDQuad idquad(force.atom0().cutGroup(),
                        force.atom1().cutGroup(),
                        force.atom2().cutGroup(),
                        force.atom3().cutGroup());
                        
        group_theta310_forces[idquad].append(force);
    }
    
    //add the groups into this set
    for (QHash< CGIDQuad,QVector<FourAtomFunction> >::iterator
                                        it = group_potentials.begin();
         it != group_potentials.end();
         ++it)
    {
        it.value().squeeze();
        
        QVector<FourAtomFunction> &group_theta012_force 
                                        = group_theta012_forces[it.key()];
        QVector<FourAtomFunction> &group_theta213_force 
                                        = group_theta213_forces[it.key()];
        QVector<FourAtomFunction> &group_theta310_force 
                                        = group_theta310_forces[it.key()];
        
        group_theta012_force.squeeze();
        group_theta213_force.squeeze();
        group_theta310_force.squeeze();
        
        this->getGroup(it.key(), cached_groups)
                        .setBendBendPotential( it.value(), group_theta012_force,
                                        group_theta213_force, group_theta310_force );
    }
}

/** This adds all of the stretch-bend-torsion parameters and forces in 'sbtparams'
    to the list of parameters */
void InternalParameters::addSBTs(const FourAtomFunctions &sbtparams,
                                QHash<CGIDQuad,qint32> &cached_groups)
{
    if (sbtparams.isEmpty())
        return;

    //assert that these are functions only of the angle (theta)
    this->assertContainsOnly(this->symbols().stretchBendTorsion(), 
                             sbtparams.symbols());
    
    //get the potential and forces
    QVector<FourAtomFunction> potentials = sbtparams.potentials();

    QVector<FourAtomFunction> phi_forces = sbtparams.forces(this->symbols()
                                                  .stretchBendTorsion().phi());

    QVector<FourAtomFunction> theta012_forces = sbtparams.forces(this->symbols()
                                                  .stretchBendTorsion().theta012());
    QVector<FourAtomFunction> theta321_forces = sbtparams.forces(this->symbols()
                                                  .stretchBendTorsion().theta321());

    QVector<FourAtomFunction> r01_forces = sbtparams.forces(this->symbols()
                                                  .stretchBendTorsion().r01());
    QVector<FourAtomFunction> r12_forces = sbtparams.forces(this->symbols()
                                                  .stretchBendTorsion().r12());
    QVector<FourAtomFunction> r32_forces = sbtparams.forces(this->symbols()
                                                  .stretchBendTorsion().r32());

    //sort the internals into groups
    QHash< CGIDQuad, QVector<FourAtomFunction> > group_potentials;

    QHash< CGIDQuad, QVector<FourAtomFunction> > group_phi_forces;

    QHash< CGIDQuad, QVector<FourAtomFunction> > group_theta012_forces;
    QHash< CGIDQuad, QVector<FourAtomFunction> > group_theta321_forces;

    QHash< CGIDQuad, QVector<FourAtomFunction> > group_r01_forces;
    QHash< CGIDQuad, QVector<FourAtomFunction> > group_r12_forces;
    QHash< CGIDQuad, QVector<FourAtomFunction> > group_r32_forces;
    
    const FourAtomFunction *potentials_array = potentials.constData();
    int n = potentials.count();
    
    for (int i=0; i<n; ++i)
    {
        const FourAtomFunction &potential = potentials_array[i];
    
        CGIDQuad idquad(potential.atom0().cutGroup(),
                        potential.atom1().cutGroup(),
                        potential.atom2().cutGroup(),
                        potential.atom3().cutGroup());
                                         
        group_potentials[idquad].append(potential);
    }

    const FourAtomFunction *phi_forces_array = phi_forces.constData();
    n = phi_forces.count();
    
    for (int i=0; i<n; ++i)
    {
        const FourAtomFunction &force = phi_forces_array[i];
        
        CGIDQuad idquad(force.atom0().cutGroup(),
                        force.atom1().cutGroup(),
                        force.atom2().cutGroup(),
                        force.atom3().cutGroup());
                        
        group_phi_forces[idquad].append(force);
    }

    const FourAtomFunction *theta012_forces_array = theta012_forces.constData();
    n = theta012_forces.count();
    
    for (int i=0; i<n; ++i)
    {
        const FourAtomFunction &force = theta012_forces_array[i];
        
        CGIDQuad idquad(force.atom0().cutGroup(),
                        force.atom1().cutGroup(),
                        force.atom2().cutGroup(),
                        force.atom3().cutGroup());
                        
        group_theta012_forces[idquad].append(force);
    }

    const FourAtomFunction *theta321_forces_array = theta321_forces.constData();
    n = theta321_forces.count();
    
    for (int i=0; i<n; ++i)
    {
        const FourAtomFunction &force = theta321_forces_array[i];
        
        CGIDQuad idquad(force.atom0().cutGroup(),
                        force.atom1().cutGroup(),
                        force.atom2().cutGroup(),
                        force.atom3().cutGroup());
                        
        group_theta321_forces[idquad].append(force);
    }

    const FourAtomFunction *r01_forces_array = r01_forces.constData();
    n = r01_forces.count();
    
    for (int i=0; i<n; ++i)
    {
        const FourAtomFunction &force = r01_forces_array[i];
        
        CGIDQuad idquad(force.atom0().cutGroup(),
                        force.atom1().cutGroup(),
                        force.atom2().cutGroup(),
                        force.atom3().cutGroup());
                        
        group_r01_forces[idquad].append(force);
    }

    const FourAtomFunction *r12_forces_array = r12_forces.constData();
    n = r12_forces.count();
    
    for (int i=0; i<n; ++i)
    {
        const FourAtomFunction &force = r12_forces_array[i];
        
        CGIDQuad idquad(force.atom0().cutGroup(),
                        force.atom1().cutGroup(),
                        force.atom2().cutGroup(),
                        force.atom3().cutGroup());
                        
        group_r12_forces[idquad].append(force);
    }

    const FourAtomFunction *r32_forces_array = r32_forces.constData();
    n = r32_forces.count();
    
    for (int i=0; i<n; ++i)
    {
        const FourAtomFunction &force = r32_forces_array[i];
        
        CGIDQuad idquad(force.atom0().cutGroup(),
                        force.atom1().cutGroup(),
                        force.atom2().cutGroup(),
                        force.atom3().cutGroup());
                        
        group_r32_forces[idquad].append(force);
    }
    
    //add the groups into this set
    for (QHash< CGIDQuad,QVector<FourAtomFunction> >::iterator
                                        it = group_potentials.begin();
         it != group_potentials.end();
         ++it)
    {
        it.value().squeeze();
        
        QVector<FourAtomFunction> &group_phi_force = group_phi_forces[it.key()];
        
        QVector<FourAtomFunction> &group_theta012_force 
                                        = group_theta012_forces[it.key()];
        QVector<FourAtomFunction> &group_theta321_force 
                                        = group_theta321_forces[it.key()];

        QVector<FourAtomFunction> &group_r01_force = group_r01_forces[it.key()];
        QVector<FourAtomFunction> &group_r12_force = group_r12_forces[it.key()];
        QVector<FourAtomFunction> &group_r32_force = group_r32_forces[it.key()];
        
        group_phi_force.squeeze();

        group_theta012_force.squeeze();
        group_theta321_force.squeeze();

        group_r01_force.squeeze();
        group_r12_force.squeeze();
        group_r32_force.squeeze();
        
        this->getGroup(it.key(), cached_groups)
                        .setStretchBendTorsionPotential( it.value(), 
                                        group_phi_force,
                                        group_theta012_force, group_theta321_force,
                                        group_r01_force, group_r12_force,
                                        group_r32_force );
    }
}

/** Construct the parameters from the specified properties of 
    the molecule 'molecule'. If 'isstrict' is true, then this
    includes only internals that are wholely contained within
    the selected atoms of the molecule. Otherwise, this contains
    internals from which at least one atom is in the selected
    atoms of the molecule */
InternalParameters::InternalParameters(const PartialMolecule &molecule,
                                       const PropertyName &bond_params,
                                       const PropertyName &angle_params,
                                       const PropertyName &dihedral_params,
                                       const PropertyName &improper_params,
                                       const PropertyName &ub_params,
                                       const PropertyName &ss_params,
                                       const PropertyName &sb_params,
                                       const PropertyName &bb_params,
                                       const PropertyName &sbt_params,
                                       bool isstrict
                                       )
                   : state(EMPTY)
{
    QHash<CGIDQuad,qint32> cached_groups;
    
    const Property &bond_property = molecule.property(bond_params);
    
    if (not bond_property.isA<NullProperty>())
    {
        const TwoAtomFunctions &bondparams = bond_property.asA<TwoAtomFunctions>();

        if (molecule.selection().selectedAll())
            this->addBonds(bondparams, cached_groups);
        else
            this->addBonds(bondparams.includeOnly(molecule.selection(), isstrict),
                           cached_groups);
    }
    
    const Property &angle_property = molecule.property(angle_params);
    
    if (not angle_property.isA<NullProperty>())
    {
        const ThreeAtomFunctions &angleparams = angle_property.asA<ThreeAtomFunctions>();
                                                        
        if (molecule.selection().selectedAll())
            this->addAngles(angleparams, cached_groups);
        else
            this->addAngles(angleparams.includeOnly(molecule.selection(), isstrict),
                            cached_groups);
    }
    
    const Property &dihedral_property = molecule.property(dihedral_params);
    
    if (not dihedral_property.isA<NullProperty>())
    {
        const FourAtomFunctions &dihedralparams 
                                        = dihedral_property.asA<FourAtomFunctions>();
                                                        
        if (molecule.selection().selectedAll())
            this->addDihedrals(dihedralparams, cached_groups);
        else
            this->addDihedrals(dihedralparams.includeOnly(molecule.selection(), 
                                                          isstrict), cached_groups);
    }
    
    const Property &improper_property = molecule.property(improper_params);
    
    if (not improper_property.isA<NullProperty>())
    {
        const FourAtomFunctions &improperparams 
                                        = improper_property.asA<FourAtomFunctions>();
                                                        
        if (molecule.selection().selectedAll())
            this->addImpropers(improperparams, cached_groups);
        else
            this->addImpropers(improperparams.includeOnly(molecule.selection(), 
                                                          isstrict), cached_groups);
    }
    
    const Property &ub_property = molecule.property(ub_params);
    
    if (not ub_property.isA<NullProperty>())
    {
        const TwoAtomFunctions &ubparams = ub_property.asA<TwoAtomFunctions>();
                                                        
        if (molecule.selection().selectedAll())
            this->addUBs(ubparams, cached_groups);
        else
            this->addUBs(ubparams.includeOnly(molecule.selection(), isstrict),
                         cached_groups);
    }
    
    const Property &ss_property = molecule.property(ss_params);
    
    if (not ss_property.isA<NullProperty>())
    {
        const ThreeAtomFunctions &ssparams = ss_property.asA<ThreeAtomFunctions>();
                                                        
        if (molecule.selection().selectedAll())
            this->addSSs(ssparams, cached_groups);
        else
            this->addSSs(ssparams.includeOnly(molecule.selection(), isstrict),
                         cached_groups);
    }
    
    const Property &sb_property = molecule.property(sb_params);
    
    if (not sb_property.isA<NullProperty>())
    {
        const ThreeAtomFunctions &sbparams = sb_property.asA<ThreeAtomFunctions>();
                                                        
        if (molecule.selection().selectedAll())
            this->addSBs(sbparams, cached_groups);
        else
            this->addSBs(sbparams.includeOnly(molecule.selection(), isstrict),
                         cached_groups);
    }
        
    const Property &bb_property = molecule.property(bb_params);
        
    if (not bb_property.isA<NullProperty>())
    {
        const FourAtomFunctions &bbparams = bb_property.asA<FourAtomFunctions>();
                                                        
        if (molecule.selection().selectedAll())
            this->addBBs(bbparams, cached_groups);
        else
            this->addBBs(bbparams.includeOnly(molecule.selection(), isstrict),
                         cached_groups);
    }
        
    const Property &sbt_property = molecule.property(sbt_params);
        
    if (not sbt_property.isA<NullProperty>())
    {
        const FourAtomFunctions &sbtparams = sbt_property.asA<FourAtomFunctions>();
                                                        
        if (molecule.selection().selectedAll())
            this->addSBTs(sbtparams, cached_groups);
        else
            this->addSBTs(sbtparams.includeOnly(molecule.selection(), isstrict),
                          cached_groups);
    }

    //squeeze all of the data to minimise memory use
    group_params.squeeze();
    
    for (QHash< CGIdx,QSet<qint32> >::iterator it = groups_by_cgidx.begin();
         it != groups_by_cgidx.end();
         ++it)
    {
        it->squeeze();
    }
    
    groups_by_cgidx.squeeze();
    
    this->updateState();
}

/** Null constructor */
InternalParameters::InternalParameters() : state(EMPTY)
{}

/** Copy constructor */
InternalParameters::InternalParameters(const InternalParameters &other)
                   : state(other.state),
                     group_params(other.group_params),
                     groups_by_cgidx(other.groups_by_cgidx)
{}

/** Destructor */
InternalParameters::~InternalParameters()
{}

/** Copy assignment operator */
InternalParameters& InternalParameters::operator=(const InternalParameters &other)
{
    state = other.state;
    group_params = other.group_params;
    groups_by_cgidx = other.groups_by_cgidx;
    
    return *this;
}

/** Comparison operator */
bool InternalParameters::operator==(const InternalParameters &other) const
{
    return state == other.state and
           group_params == other.group_params;
}

/** Comparison operator */
bool InternalParameters::operator!=(const InternalParameters &other) const
{
    return state != other.state or
           group_params != other.group_params;
}

/** Return whether all of the parameters for all CutGroups have changed
    compared to 'other' */
bool InternalParameters::changedAllGroups(const InternalParameters &other) const
{
    if (groups_by_cgidx == other.groups_by_cgidx)
    {
        BOOST_ASSERT( group_params.count() == other.group_params.count() );
        
        int n = group_params.count();
        const GroupInternalParameters *this_array = group_params.constData();
        const GroupInternalParameters *other_array = other.group_params.constData();
        
        if (this_array == other_array)
            return false;
            
        for (int i=0; i<n; ++i)
        {
            if (this_array[i] == other_array[i])
                return false;
        }
        
        return true;
    }
    else
        return true;
}

/** Add the indicies of CutGroups that have changed on to 'changed_groups' */
void InternalParameters::addChangedGroups(const InternalParameters &other, 
                                          QSet<quint32> &changed_groups) const
{
    if (this->changedAllGroups(other))
    {
        for (QHash< CGIdx,QSet<qint32> >::const_iterator 
                                                    it = groups_by_cgidx.constBegin();
             it != groups_by_cgidx.constEnd();
             ++it)
        {
            changed_groups.insert(it.key());
        }
        
        return;
    }
    else
    {
        //the groups_by_cgidx objects for both set of parameters are now
        //guaranteed to be equal, so the group_params arrays are directly
        //comparable
        int n = group_params.count();
        
        BOOST_ASSERT(group_params.count() == other.group_params.count());
        
        const GroupInternalParameters *this_array = group_params.constData();
        const GroupInternalParameters *other_array = other.group_params.constData();
        
        if (this_array == other_array)
            return;
            
        for (int i=0; i<n; ++i)
        {
            if (this_array[i] != other_array[i])
            {
                const GroupInternalParameters &group = this_array[i];
                
                if (group.isSingleCutGroup())
                {
                    changed_groups.insert(group.cgIdx0());
                }
                else if (group.isDoubleCutGroup())
                {
                    changed_groups.insert(group.cgIdx0());
                    changed_groups.insert(group.cgIdx1());
                }
                else if (group.isTripleCutGroup())
                {
                    changed_groups.insert(group.cgIdx0());
                    changed_groups.insert(group.cgIdx1());
                    changed_groups.insert(group.cgIdx2());
                }
                else
                {
                    changed_groups.insert(group.cgIdx0());
                    changed_groups.insert(group.cgIdx1());
                    changed_groups.insert(group.cgIdx2());
                    changed_groups.insert(group.cgIdx3());
                }
            }
        }
        
        return;
    }
}
        
/** Return the CGIdxs of the CutGroups that have changed from this group
    compared to the paramters held in 'other' */              
QSet<quint32> InternalParameters::getChangedGroups(const InternalParameters &other) const
{
    QSet<quint32> changed_groups;
    
    this->addChangedGroups(other, changed_groups);
    
    return changed_groups;
}

/** Return whether or not this contains only parameters that affect the atoms
    in the CutGroups whose indicies are in 'cgidxs' */
bool InternalParameters::containsOnly(const QSet<quint32> &cgidxs) const
{
    if (cgidxs.count() >= groups_by_cgidx.count())
    {
        //do we only contain valid groups?
        for (QHash< CGIdx, QSet<qint32> >::const_iterator
                                            it = groups_by_cgidx.constBegin();
             it != groups_by_cgidx.constEnd();
             ++it)
        {
            if (not cgidxs.contains(it.key()))
                return false;
        }
        
        return true;
    }
    
    return false;
}

/** Update the state */
void InternalParameters::updateState()
{
    state = EMPTY;
    
    if (group_params.isEmpty())
        return;
        
    int ngroups = group_params.count();
    const GroupInternalParameters *groups_array = group_params.constData();
    
    for (int i=0; i<ngroups; ++i)
    {
        const GroupInternalParameters &group = groups_array[i];
        
        if (group.hasPhysicalParameters())
        {
            if (not group.bondPotential().isEmpty())
                state |= HAS_BOND;
                
            if (not group.anglePotential().isEmpty())
                state |= HAS_ANGLE;
                
            if (not group.dihedralPotential().isEmpty())
                state |= HAS_DIHEDRAL;
        }
        
        if (group.hasNonPhysicalParameters())
        {
            if (not group.improperPotential().isEmpty())
                state |= HAS_IMPROPER;
            
            if (not group.ureyBradleyPotential().isEmpty())
                state |= HAS_UB;
        }
        
        if (group.hasCrossTerms())
        {
            if (not group.stretchStretchPotential().isEmpty())
                state |= HAS_SS;
                
            if (not group.stretchBendPotential().isEmpty())
                state |= HAS_SB;
                
            if (not group.bendBendPotential().isEmpty())
                state |= HAS_BB;
                
            if (not group.stretchBendTorsionPotential().isEmpty())
                state |= HAS_SBT;
        }
    }
}

/** Reindex the parameters */
void InternalParameters::reindex()
{
    if (group_params.isEmpty())
    {
        groups_by_cgidx.clear();
        return;
    }
    
    qint32 n = group_params.count();
    const GroupInternalParameters *groups_array = group_params.constData();
    
    for (qint32 i=0; i<n; ++i)
    {
        const GroupInternalParameters &group = groups_array[i];
        
        if (group.isSingleCutGroup())
        {
            groups_by_cgidx[group.cgIdx0()].insert(i);
        }
        else if (group.isDoubleCutGroup())
        {
            groups_by_cgidx[group.cgIdx0()].insert(i);
            groups_by_cgidx[group.cgIdx1()].insert(i);
        }
        else if (group.isTripleCutGroup())
        {
            groups_by_cgidx[group.cgIdx0()].insert(i);
            groups_by_cgidx[group.cgIdx1()].insert(i);
            groups_by_cgidx[group.cgIdx2()].insert(i);
        }
        else
        {
            groups_by_cgidx[group.cgIdx0()].insert(i);
            groups_by_cgidx[group.cgIdx1()].insert(i);
            groups_by_cgidx[group.cgIdx2()].insert(i);
            groups_by_cgidx[group.cgIdx3()].insert(i);
        }
    }
    
    for (QHash< CGIdx,QSet<qint32> >::iterator it = groups_by_cgidx.begin();
         it != groups_by_cgidx.end();
         ++it)
    {
        it->squeeze();
    }
    
    groups_by_cgidx.squeeze();
    
    this->updateState();
}

/** Return all of the parameters that involve any of the CutGroups whose
    indicies are in 'cgidxs' */
QVector<GroupInternalParameters> InternalParameters::groupParameters(
                                                    const QSet<quint32> &cgidxs) const
{
    //build up the set of groups to include...
    QSet<qint32> idxs;
    idxs.reserve(group_params.count());
    
    foreach (quint32 cgidx, cgidxs)
    {
        idxs += groups_by_cgidx.value( CGIdx(cgidx) );
    }
    
    if (idxs.isEmpty())
        return QVector<GroupInternalParameters>();

    else if (idxs.count() == group_params.count())
        return group_params;
        
    QVector<GroupInternalParameters> params(idxs.count());
    
    const GroupInternalParameters *params_array = group_params.constData();
    GroupInternalParameters *new_array = params.data();
    
    int i = 0;
    
    foreach (qint32 idx, idxs)
    {
        new_array[i] = params_array[idx];
        ++i;
    }
    
    return params;
}

/** Mask this set so that only the parameters for the specified CutGroups are
    included in the returned group. */
InternalParameters InternalParameters::applyMask(const QSet<quint32> &cgidxs) const
{
    if (this->containsOnly(cgidxs))
        return *this;
        
    InternalParameters ret;
    
    ret.group_params = this->groupParameters(cgidxs);
    ret.reindex();
    
    return ret;
}

/** Return the array of all of the internal parameters */
const QVector<GroupInternalParameters>& InternalParameters::groupParameters() const
{
    return group_params;
}

/** Return the array of all of the parameters that involve the CutGroup with index
    'cgidx'. This returns an empty array if there are no parameters for this CutGroup */
QVector<GroupInternalParameters> InternalParameters::groupParameters(quint32 cgidx) const
{
    QSet<qint32> idxs = groups_by_cgidx.value( CGIdx(cgidx) );
    
    if (idxs.isEmpty())
        return QVector<GroupInternalParameters>();
        
    QVector<GroupInternalParameters> params(idxs.count());
    GroupInternalParameters *new_array = params.data();
    
    const GroupInternalParameters *groups_array = group_params.constData();
    
    int i = 0;
    
    foreach (qint32 idx, idxs)
    {
        new_array[i] = groups_array[idx];
        ++i;
    }
    
    return params;
}

const char* InternalParameters::typeName()
{
    return QMetaType::typeName( qMetaTypeId<InternalParameters>() );
}

//////////
////////// Implementation of InternalParameters3D
//////////

static const RegisterMetaType<InternalParameters3D> r_params3d(NO_ROOT);

/** Serialise to a binary datastream */
QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds,
                                      const InternalParameters3D &params3d)
{
    writeHeader(ds, r_params3d, 1);
    
    SharedDataStream sds(ds);
    
    sds << static_cast<const InternalParameters&>(params3d)
        << static_cast<const AtomicCoords3D&>(params3d);
        
    return ds;
}              

/** Extract from a binary datastream */
QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds,
                                      InternalParameters3D &params3d)
{
    VersionID v = readHeader(ds, r_params3d);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> static_cast<InternalParameters&>(params3d)
            >> static_cast<AtomicCoords3D&>(params3d);
    }
    else
        throw version_error(v, "1", r_params3d, CODELOC);
        
    return ds;
}                        

/** Null constructor */
InternalParameters3D::InternalParameters3D()
                     : InternalParameters(), AtomicCoords3D()
{}

/** Internal constructor */
InternalParameters3D::InternalParameters3D(const AtomicCoords3D &coords,
                                           const InternalParameters &params)
                     : InternalParameters(params),
                       AtomicCoords3D(coords)
{}

/** Construct, creating the parameters from the passed molecule 
    using the supplied property names
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
InternalParameters3D::InternalParameters3D(const PartialMolecule &molecule,
                                           const PropertyName &coords_property,
                                           const PropertyName &bond_params,
                                           const PropertyName &angle_params,
                                           const PropertyName &dihedral_params,
                                           const PropertyName &improper_params,
                                           const PropertyName &ub_params,
                                           const PropertyName &ss_params,
                                           const PropertyName &sb_params,
                                           const PropertyName &bb_params,
                                           const PropertyName &sbt_params,
                                           bool isstrict)
     : InternalParameters(molecule, bond_params, angle_params,
                          dihedral_params, improper_params,
                          ub_params, ss_params, sb_params,
                          bb_params, sbt_params, isstrict),
       AtomicCoords3D( molecule.molecule(), coords_property )
{}

/** Copy constructor */
InternalParameters3D::InternalParameters3D(const InternalParameters3D &other)
                     : InternalParameters(other), AtomicCoords3D(other)
{}

/** Destructor */
InternalParameters3D::~InternalParameters3D()
{}

/** Copy assignment operator */
InternalParameters3D& InternalParameters3D::operator=(const InternalParameters3D &other)
{
    InternalParameters::operator=(other);
    AtomicCoords3D::operator=(other);
    
    return *this;
}

/** Comparison operator */
bool InternalParameters3D::operator==(const InternalParameters3D &other) const
{
    return AtomicCoords3D::operator==(other) and 
           InternalParameters::operator==(other);
}

/** Comparison operator */
bool InternalParameters3D::operator!=(const InternalParameters3D &other) const
{
    return AtomicCoords3D::operator!=(other) or
           InternalParameters::operator!=(other);
}

/** Return the coordinates */
const CoordGroupArray& InternalParameters3D::atomicCoordinates() const
{
    return AtomicCoords3D::atomicCoordinates();
}

/** Set the coordinates used by these parameters */
void InternalParameters3D::setAtomicCoordinates(const AtomicCoords3D &coords)
{
    const int ngroups = coords.atomicCoordinates().count();

    //there must be the same number of CutGroups as in the parameters!
    if (ngroups != atomicCoordinates().count())
    {
        throw SireError::program_bug( QObject::tr(
            "Error setting incompatible coordinates!!! %1 vs. %2")
                .arg(ngroups)
                .arg(atomicCoordinates().count()),
                    CODELOC );
    }
    
    const CoordGroup *this_array = atomicCoordinates().constData();
    const CoordGroup *other_array = coords.atomicCoordinates().constData();
    
    for (int i=0; i<ngroups; ++i)
    {
        if (this_array[i].count() != other_array[i].count())
            throw SireError::program_bug( QObject::tr(
                "Error setting incompatible coordinates in group %1. %2 vs. %3.")
                    .arg(i).arg(this_array[i].count())
                    .arg(other_array[i].count()),
                        CODELOC );
    }
    
    //ok, the coordinates are compatible
    AtomicCoords3D::operator=(coords);
}

/** Return the number of CutGroups in the molecule whose parameters are
    contained in this object */
int InternalParameters3D::nCutGroups() const
{
    return AtomicCoords3D::atomicCoordinates().count();
}

/** Return whether or not all of the CutGroup have changed compared to 'other' */
bool InternalParameters3D::changedAllGroups(const InternalParameters3D &other) const
{
    if (AtomicCoords3D::changedAllGroups(other) or
        InternalParameters::changedAllGroups(other))
    {
        return true;
    }
    else
    {
        //maybe some coordinates have changed, and some internal parameters.
        // However it takes too long to test for this, so we won't!
        return false;
    }
}

/** Add the changed groups that are different in 'other' compared to this
    to 'changed_groups' */
void InternalParameters3D::addChangedGroups(const InternalParameters3D &other, 
                                            QSet<quint32> &changed_groups) const
{
    const CoordGroup *this_coords_array = this->atomicCoordinates().constData();

    const CoordGroup *other_coords_array = other.atomicCoordinates().constData();

    if (this_coords_array != other_coords_array)
    {
        quint32 ngroups = qMin(this->atomicCoordinates().count(),
                               other.atomicCoordinates().count());
                       
        for (quint32 i=0; i<ngroups; ++i)
        {
            if ( this_coords_array[i] != other_coords_array[i] )
                changed_groups.insert(i);
        }
    }
    
    //check to see whether all of the groups have changed
    if (changed_groups.count() >= this->nCutGroups())
    {
        int count = 0;
    
        foreach (quint32 cgidx, changed_groups)
        {
            if ( cgidx < quint32(this->nCutGroups()) )
            {
                ++count;
            }
        }
        
        if (count >= this->nCutGroups())
            return;
    }
    
    //now add on changed parameters
    InternalParameters::addChangedGroups(other, changed_groups);
}
                      
/** Return the indicies of the CutGroups that have changed in 'other' compared
    to this set of parameters */
QSet<quint32> InternalParameters3D::getChangedGroups(
                                        const InternalParameters3D &other) const
{
    QSet<quint32> changed_groups;
    
    this->addChangedGroups(other, changed_groups);
    
    return changed_groups;
}

/** Mask these parameters so that only the parameters for the CutGroups
    whose indicies are in 'cgidxs' are contained. */
InternalParameters3D InternalParameters3D::applyMask(const QSet<quint32> &cgidxs) const
{
    return InternalParameters3D( *this, InternalParameters::applyMask(cgidxs) );
}

const char* InternalParameters3D::typeName()
{
    return QMetaType::typeName( qMetaTypeId<InternalParameters3D>() );
}
