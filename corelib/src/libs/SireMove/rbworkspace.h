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

#ifndef SIREMOVE_RBWORKSPACE_H
#define SIREMOVE_RBWORKSPACE_H

#include "integratorworkspace.h"

#include "SireMaths/matrix.h"
#include "SireMaths/quaternion.h"

SIRE_BEGIN_HEADER

namespace SireMove
{
class RBWorkspace;
}

QDataStream& operator<<(QDataStream&, const SireMove::RBWorkspace&);
QDataStream& operator>>(QDataStream&, SireMove::RBWorkspace&);

namespace SireMove
{

using SireMaths::Matrix;
using SireMaths::Quaternion;

using SireMol::ViewsOfMol;

/** This class provides a workspace for integrators that perform
    rigid body integration of atomic velocities and coordinates
    
    @author Christopher Woods
*/
class SIREMOVE_EXPORT RBWorkspace
       : public SireBase::ConcreteProperty<RBWorkspace,IntegratorWorkspace>
{

friend QDataStream& ::operator<<(QDataStream&, const RBWorkspace&);
friend QDataStream& ::operator>>(QDataStream&, RBWorkspace&);

public:
    RBWorkspace(const PropertyMap &map = PropertyMap());
    RBWorkspace(const MoleculeGroup &molgroup,
                const PropertyMap &map = PropertyMap());
    
    RBWorkspace(const RBWorkspace &other);
    
    ~RBWorkspace();

    RBWorkspace& operator=(const RBWorkspace &other);
    
    bool operator==(const RBWorkspace &other) const;
    bool operator!=(const RBWorkspace &other) const;
    
    static const char* typeName();
    
    SireUnits::Dimension::MolarEnergy kineticEnergy() const;
    SireUnits::Dimension::MolarEnergy kineticEnergy(MolNum molnum) const;
    SireUnits::Dimension::MolarEnergy kineticEnergy(const MoleculeView &molview) const;

    PropertyName beadingProperty() const;

    int nBeads() const;

    int nAtoms(int ibead) const;    
    
    bool setSystem(const System &system);

    bool calculateForces(const Symbol &nrg_component);

    void regenerateVelocities(const VelocityGenerator &generator);
    
    void commitCoordinates();
    void commitVelocities();
    
    void commitCoordinatesAndVelocities();

    Vector* beadCoordsArray();
    Quaternion* beadOrientationArray();
    
    Vector* beadLinearMomentaArray();
    Vector* beadAngularMomentaArray();

    const double* beadMassesArray() const;
    const Vector* beadInertiasArray() const;
        
    const Vector* beadForcesArray() const;
    const Vector* beadTorquesArray() const;

protected:
    void changedProperty(const QString &property);

private:
    void rebuildFromScratch();

    /** The coordinates of the atoms of each bead in the center of
        mass / principle inertia tensor frame. These
        will be constant, as the atoms don't move relative
        to the center of mass / inertia frame - they 
        are used to convert from that frame to the
        world cartesian frame */
    QVector< QVector<Vector> > atom_int_coords;

    /** The index of the bead in which each atom exists */
    QVector< QPair< qint32,QVector<qint32> > > atoms_to_beads;

    /** The center of mass coordinates of each bead */
    QVector<Vector> bead_coordinates;
    
    /** The matrix to map from bead internal coordinates
        to World cartesian coordinates */
    QVector<Matrix> bead_to_world;
    
    /** The orientation quaternion for each bead */
    QVector<Quaternion> bead_orientations;
    
    /** All of the beads' linear momenta */
    QVector<Vector> bead_linear_momenta;
    
    /** All of the beads' angular momenta */
    QVector<Vector> bead_angular_momenta;
    
    /** All of the forces for all of the beads */
    QVector<Vector> bead_forces;
    
    /** All of the torques acting on all of the beads */
    QVector<Vector> bead_torques;
    
    /** The mass of each bead */
    QVector<double> bead_masses;
    
    /** The diagonals of the inertia tensor for each bead */
    QVector<Vector> bead_inertia;
    
    /** The generator used to get the initial 
        linear and angular velocities */
    VelGenPtr vel_generator;
};

}

Q_DECLARE_METATYPE( SireMove::RBWorkspace )

SIRE_END_HEADER

#endif
