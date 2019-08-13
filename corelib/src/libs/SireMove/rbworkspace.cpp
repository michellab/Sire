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

#include "stdio.h"

#include "rbworkspace.h"

#include "SireMol/atomelements.h"
#include "SireMol/atommasses.h"
#include "SireMol/atomcoords.h"
#include "SireMol/molecule.h"
#include "SireMol/moleditor.h"

#include "SireMaths/axisset.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMove;
using namespace SireFF;
using namespace SireMol;
using namespace SireSystem;
using namespace SireMaths;
using namespace SireBase;
using namespace SireUnits::Dimension;
using namespace SireStream;

static const RegisterMetaType<RBWorkspace> r_rbws;

QDataStream &operator<<(QDataStream &ds, const RBWorkspace &rbws)
{
    writeHeader(ds, r_rbws, 1);
    
    SharedDataStream sds(ds);
    
    sds << static_cast<const IntegratorWorkspace&>(rbws);
        
    return ds;
}

QDataStream &operator>>(QDataStream &ds, RBWorkspace &rbws)
{
    VersionID v = readHeader(ds, r_rbws);
    
    if (v == 1)
    {
        RBWorkspace rb;
    
        SharedDataStream sds(ds);
    
        sds >> static_cast<IntegratorWorkspace&>(rb);
        
        rb.rebuildFromScratch();
        
        rbws = rb;
    }
    else
        throw version_error(v, "1", r_rbws, CODELOC);
        
    return ds;
}

SIRE_ALWAYS_INLINE double getMass(const MolarMass &mass)
{
    return mass.value();
}

SIRE_ALWAYS_INLINE double getMass(const Element &element)
{
    return element.mass().value();
}

static Vector getPrincipalAxes(Matrix &inertia)
{
    double *inertia_array = inertia.data();

    //remove near-zero elements
    for (int i=0; i<9; ++i)
    {
        if (inertia_array[i] < 1e-6 and inertia_array[i] > -1e-6)
            inertia_array[i] = 0;
    }
    
    //symmetric matrix
    inertia_array[inertia.offset(1,0)] = inertia_array[inertia.offset(0,1)];
    inertia_array[inertia.offset(2,0)] = inertia_array[inertia.offset(0,2)];
    inertia_array[inertia.offset(2,1)] = inertia_array[inertia.offset(1,2)];

    boost::tuple<Vector,Matrix> eigs = inertia.diagonalise();

    Vector &principle_inertia = eigs.get<0>();
    Matrix &orientation = eigs.get<1>();
    
    //if one or more of the eigenvalues is zero then we may have a problem
    //because the wrong eigenvector direction may be chosen - in this case,
    //we will build this eigenvector using a cross product to ensure that 
    //the right-hand-rule definition of our axes is maintained
    //
    // Also, even if we have three eigenvalues, we still need to make sure
    // that a right-hand-rule set is chosen, rather than the left-hand set
    bool zero_x = std::abs(principle_inertia[0]) < 1e-6;
    bool zero_y = std::abs(principle_inertia[1]) < 1e-6;
    bool zero_z = std::abs(principle_inertia[2]) < 1e-6;
    
    int n_zeroes = int(zero_x) + int(zero_y) + int(zero_z);
    
    if (n_zeroes == 3)
    {
        //no axes!
        orientation = Matrix(1);
    }
    else if (n_zeroes == 2)
    {
        //just one well-defined axis - I don't know how to handle this...
        throw SireError::incompatible_error( QObject::tr(
                "Sire cannot yet handle rigid body systems with only "
                "a single non-zero eigenvalue - %1, moment of interia equals\n%2")
                    .arg(principle_inertia.toString(),
                         orientation.toString()), CODELOC );
    }
    else if (n_zeroes == 1)
    {
        Vector r0 = orientation.row0();
        Vector r1 = orientation.row1();
        Vector r2 = orientation.row2();
        
        if (zero_x)
            r0 = Vector::cross(r1,r2);
        else if (zero_y)
            r1 = Vector::cross(r2,r0);
        else if (zero_z)
            r2 = Vector::cross(r0,r1);
        
        orientation = Matrix(r0, r1, r2);
    }
    else
    {
        Vector r0 = orientation.row0();
        Vector r1 = orientation.row1();
                 
        orientation = Matrix( r0, r1, Vector::cross(r0,r1) );
    }

    inertia = orientation;
    return principle_inertia;
}

template<class T>
static QVector<Vector> buildBead(const ViewsOfMol &mol, 
                                 const AtomCoords &coords,
                                 const AtomProperty<T> &masses,
                                 const QVector<qint32> &beading,
                                 Vector *bead_coords, Matrix *beads_to_world,
                                 Quaternion *bead_orients, double *bead_masses,
                                 Vector *bead_inertias)
{
    QVector<Vector> int_atom_coords;

    const MoleculeData &moldata = mol.data();

    if (mol.selectedAll())
    {
        int nats = moldata.info().nAtoms();
    
        int_atom_coords = QVector<Vector>(nats);
        int_atom_coords.squeeze();
        
        Vector *int_coords_array = int_atom_coords.data();
        
        const Vector *coords_array = coords.array().constCoordsData();
        const T *masses_array = masses.array().constValueData();
        
        if (beading.isEmpty())
        {
            //all atoms are part of the same bead
            Vector &bead_com = bead_coords[0];
            double &bead_mass = bead_masses[0];
            Matrix &bead_to_world = beads_to_world[0];
            Vector &bead_inertia = bead_inertias[0];
            
            //calculate the COM of the bead
            for (int i=0; i<nats; ++i)
            {
                const double mass = ::getMass(masses_array[i]);
            
                bead_com += mass * coords_array[i];
                bead_mass += mass;
            }
            
            bead_com /= bead_mass;
            bead_orients[0] = Quaternion();
            
            bead_to_world = Matrix( double(0) );
            
            //now calculate moments of inertia
            for (int i=0; i<nats; ++i)
            {
                double *inertia_array = bead_to_world.data();

                Vector d = coords_array[i] - bead_com;

                const double mass = ::getMass(masses_array[i]);
        
                inertia_array[ bead_to_world.offset(0,0) ] += 
                                                mass * (d.y()*d.y() + d.z()*d.z());
                inertia_array[ bead_to_world.offset(1,1) ] += 
                                                mass * (d.x()*d.x() + d.z()*d.z());
                inertia_array[ bead_to_world.offset(2,2) ] += 
                                                mass * (d.x()*d.x() + d.y()*d.y());
        
                inertia_array[ bead_to_world.offset(0,1) ] -= mass * d.x() * d.y();
                inertia_array[ bead_to_world.offset(0,2) ] -= mass * d.x() * d.z();
                inertia_array[ bead_to_world.offset(1,2) ] -= mass * d.y() * d.z();
            }

            bead_inertia = ::getPrincipalAxes(bead_to_world);
            
            //now calculate the coordinates of all of the atoms in terms
            //of the center of mass / orientaton frame
            Matrix inv_matrix = bead_to_world.inverse();

            for (int i=0; i<nats; ++i)
            {
                int_coords_array[i] = inv_matrix * (coords_array[i] - bead_com);
            }
        }
        else
        {
            //the molecule is split into several beads
            const qint32 *beading_array = beading.constData();
            int nbeads = 0;
        
            //calculate the COM of each bead
            for (int i=0; i<nats; ++i)
            {
                qint32 bead_idx = beading_array[i];
            
                if (bead_idx == -1)
                    //this atom is not part of any bead
                    continue;
                
                if (bead_idx+1 > nbeads)
                    nbeads = bead_idx + 1;
                
                const double mass = ::getMass(masses_array[i]);
                
                bead_coords[bead_idx] += mass * coords_array[i];
                bead_masses[bead_idx] += mass;
            }
        
            for (int i=0; i<nbeads; ++i)
            {
                if (bead_masses[i] != 0)
                    bead_coords[i] /= bead_masses[i];
                else
                    bead_coords[i] = Vector(0);
                    
                bead_orients[i] = Quaternion();
            }
        
            //now calculate moments of inertia for each bead
            for (int i=0; i<nats; ++i)
            {
                qint32 bead_idx = beading_array[i];
                
                if (bead_idx == -1)
                    continue;
                    
                Matrix &bead_to_world = beads_to_world[bead_idx];
                double *inertia_array = bead_to_world.data();

                Vector d = coords_array[i] - bead_coords[bead_idx];
        
                const double mass = ::getMass(masses_array[i]);
        
                inertia_array[ bead_to_world.offset(0,0) ] += 
                                                mass * (d.y()*d.y() + d.z()*d.z());
                inertia_array[ bead_to_world.offset(1,1) ] += 
                                                mass * (d.x()*d.x() + d.z()*d.z());
                inertia_array[ bead_to_world.offset(2,2) ] += 
                                                mass * (d.x()*d.x() + d.y()*d.y());
        
                inertia_array[ bead_to_world.offset(0,1) ] -= mass * d.x() * d.y();
                inertia_array[ bead_to_world.offset(0,2) ] -= mass * d.x() * d.z();
                inertia_array[ bead_to_world.offset(1,2) ] -= mass * d.y() * d.z();
            }

            for (int i=0; i<nbeads; ++i)
            {
                bead_inertias[i] = ::getPrincipalAxes(beads_to_world[i]);
            }
            
            //now calculate the coordinates of all of the atoms in terms
            //of the center of mass / orientaton frame
            qint32 last_idx = -1;
            Matrix inv_matrix;

            for (int i=0; i<nats; ++i)
            {
                qint32 bead_idx = beading_array[i];
        
                if (bead_idx == -1)
                    continue;
    
                if (bead_idx != last_idx)
                {
                    last_idx = bead_idx;
                    inv_matrix = beads_to_world[bead_idx].inverse();
                }
    
                int_coords_array[i] = inv_matrix
                                        * (coords_array[i] - bead_coords[bead_idx]);
            }
        }
    }
    else
        throw SireError::unsupported( QObject::tr(
                "The code to support rigid body dynamics of only part "
                "of a molecule has yet to be written..."), CODELOC );

    return int_atom_coords;
}

static QVector<Vector> buildBead(const ViewsOfMol &mol, const QVector<qint32> &beading,
                                 const PropertyName &coords_property,
                                 const PropertyName &masses_property,
                                 const PropertyName &elements_property,
                                 Vector *bead_coords, Matrix *bead_to_world,
                                 Quaternion *bead_orients, double *bead_masses,
                                 Vector *bead_inertia)
{
    const MoleculeData &moldata = mol.data();
    
    if (moldata.hasProperty(masses_property))
    {
        return ::buildBead(mol, moldata.property(coords_property).asA<AtomCoords>(),
                                moldata.property(masses_property).asA<AtomMasses>(),
                                beading, bead_coords, bead_to_world,
                                bead_orients, bead_masses, bead_inertia );
    }
    else
    {
        return ::buildBead(mol, moldata.property(coords_property).asA<AtomCoords>(),
                                moldata.property(elements_property).asA<AtomElements>(),
                                beading, bead_coords, bead_to_world,
                                bead_orients, bead_masses, bead_inertia);
    }
}

/** Return the property used to bead up a molecule */
PropertyName RBWorkspace::beadingProperty() const
{
    return propertyMap()["beading"];
}

static void getBeading(const ViewsOfMol &mol, const PropertyName &beading_property,
                       QPair< qint32,QVector<qint32> > &beading, qint32 &nbeads)
{
    const MoleculeData &moldata = mol.data();
    
    if (moldata.hasProperty(beading_property))
    {
        qDebug() << CODELOC;
    
        if (mol.selectedAll())
        {
            throw SireError::incomplete_code( QObject::tr("TODO"), CODELOC );
        }
        else
        {
            throw SireError::unsupported( QObject::tr(
                    "The code to support rigid body dynamics of partial "
                    "molecules has yet to be written."), CODELOC );
        }
    }
    else
    {
        beading.first = nbeads;
        beading.second = QVector<qint32>();
        nbeads += 1;
    }
}

/** Rebuild all of the data array from the current state of the system */
void RBWorkspace::rebuildFromScratch()
{
    const System &sys = this->system();
    
    PropertyName coords_property = this->coordinatesProperty();
    PropertyName masses_property = this->massesProperty();
    PropertyName elements_property = this->elementsProperty();
    PropertyName velgen_property = this->velocityGeneratorProperty();
    PropertyName beading_property = this->beadingProperty();
    
    const MoleculeGroup &molgroup = this->moleculeGroup();
    
    int nmols = molgroup.nMolecules();
    
    if (sys.containsProperty(velgen_property))
        vel_generator = sys.property(velgen_property).asA<VelocityGenerator>();
    else
        vel_generator = NullVelocityGenerator();
    
    atom_int_coords = QVector< QVector<Vector> >(nmols);
    atoms_to_beads = QVector< QPair< qint32,QVector<qint32> > >(nmols,
                                 QPair< qint32,QVector<qint32> >(-1,QVector<qint32>()) );
    
    atom_int_coords.squeeze();
    atoms_to_beads.squeeze();
    
    QVector<Vector> *atom_int_coords_array = atom_int_coords.data();
    QPair< qint32,QVector<qint32> > *atoms_to_beads_array = atoms_to_beads.data();

    //bead up the molecules
    qint32 nbeads = 0;
    
    for (int i=0; i<nmols; ++i)
    {
        MolNum molnum = molgroup.molNumAt(i);
        ::getBeading(molgroup[molnum], beading_property, atoms_to_beads_array[i], nbeads);
    }
    
    if (nbeads == 0)
    {
        bead_coordinates.clear();
        bead_orientations.clear();
        bead_to_world.clear();
        bead_masses.clear();
        bead_inertia.clear();
        return;
    }

    //now build all of the beads
    bead_coordinates = QVector<Vector>(nbeads, Vector(0));
    bead_orientations = QVector<Quaternion>(nbeads);
    bead_to_world = QVector<Matrix>(nbeads, Matrix(double(0)));
    bead_masses = QVector<double>(nbeads, 0.0);
    bead_inertia = QVector<Vector>(nbeads, Vector(0));
    
    bead_coordinates.squeeze();
    bead_orientations.squeeze();
    bead_to_world.squeeze();
    bead_masses.squeeze();
    bead_inertia.squeeze();
    
    Vector *bead_coords_array = bead_coordinates.data();
    Matrix *bead_to_world_array = bead_to_world.data();
    Quaternion *bead_orients_array = bead_orientations.data();
    double *bead_masses_array = bead_masses.data();
    Vector *bead_inertia_array = bead_inertia.data();
    
    for (int i=0; i<nmols; ++i)
    {
        const QPair< qint32,QVector<qint32> > &beading = atoms_to_beads_array[i];

        if (beading.first == -1)
            //there are no beads for this molecule
            continue;

        atom_int_coords_array[i] = ::buildBead(molgroup[molgroup.molNumAt(i)], 
                                               beading.second, 
                                               coords_property, masses_property,
                                               elements_property,
                                               bead_coords_array + beading.first,
                                               bead_to_world_array + beading.first,
                                               bead_orients_array + beading.first,
                                               bead_masses_array + beading.first,
                                               bead_inertia_array + beading.first);
    }
    
    //create space for the forces, torques and momenta
    bead_forces = QVector<Vector>(nbeads, Vector(0));
    bead_torques = QVector<Vector>(nbeads, Vector(0));
    bead_linear_momenta = QVector<Vector>(nbeads, Vector(0));
    bead_angular_momenta = QVector<Vector>(nbeads, Vector(0));
    
    bead_forces.squeeze();
    bead_torques.squeeze();
    bead_linear_momenta.squeeze();
    bead_angular_momenta.squeeze();
}

static Vector cross(const Vector &v0, const Vector &v1)
{
    return Vector( v0.y()*v1.z() - v0.z()*v1.y(),
                   v0.z()*v1.x() - v0.x()*v1.z(),
                   v0.x()*v1.y() - v0.y()*v1.x() );
}

static void calculateForces(const ViewsOfMol &mol, const MolForceTable &forces,
                            const QVector<Vector> &atom_int_coords,
                            const QVector<qint32> &beading,
                            const Matrix *beads_to_world,
                            const Quaternion *bead_orients,
                            Vector *bead_forces, Vector *bead_torques)
{
    if (mol.selectedAll())
    {
        const Vector *atomforces = forces.constValueData();
        const int nats = atom_int_coords.count();
        
        BOOST_ASSERT( forces.nValues() == nats );
        
        if (beading.isEmpty())
        {
            //there is only a single bead for the molecule
            const Matrix &bead_to_world = beads_to_world[0];
            const Quaternion &bead_orient = bead_orients[0];
            
            Vector &bead_force = bead_forces[0];
            Vector &bead_torque = bead_torques[0];
            
            if (nats == 1)
            {
                //and this is made of a single atom - force and no torque
                bead_force = atomforces[0];
                bead_torque = Vector(0);
            }
            else
            {
                Matrix orient = bead_orient.toMatrix() * bead_to_world;
                const Vector *int_coords = atom_int_coords.constData();
            
                for (int i=0; i<nats; ++i)
                {
                    bead_force += atomforces[i];
                
                    //calculate the vector from the center of mass to 
                    //the atom, in the World cartesian frame
                    Vector r = orient * int_coords[i];

                    //the torque is r cross force (need unnormalised cross product)
                    bead_torque -= ::cross(r, atomforces[i]);
                }

                //map the torque back from the cartesian frame to the 
                //internal frame
                bead_torque = orient.inverse() * bead_torque;
            }
        }
        else
        {
            const qint32 *beading_array = beading.constData();
            const Vector *int_coords = atom_int_coords.constData();
        
            qint32 last_idx = -1;
            Matrix orient;
            
            int nbeads = 0;
            
            for (int i=0; i<nats; ++i)
            {
                qint32 bead_idx = beading_array[i];
                
                if (i == -1)
                    //this atom is not part of a bead
                    continue;
                    
                if (bead_idx != last_idx)
                {
                    last_idx = bead_idx;
                    orient = bead_orients[bead_idx].toMatrix() * beads_to_world[bead_idx];
                    
                    if (bead_idx+1 > nbeads)
                        nbeads = bead_idx + 1;
                }
            
                bead_forces[bead_idx] += atomforces[i];
                
                //calculate the vector from the center of mass to 
                //the atom, in the World cartesian frame
                Vector r = orient * int_coords[i];

                //the torque is r cross force (need unnormalised cross product)
                bead_torques[bead_idx] -= ::cross(r, atomforces[i]);
            }

            for (int i=0; i<nbeads; ++i)
            {
                //map the torque back from the cartesian frame to the 
                //internal frame
                Matrix orient = bead_orients[i].toMatrix() * beads_to_world[i];

                bead_torques[i] = orient.inverse() * bead_torques[i];
            }
        }
    }
    else
        throw SireError::unsupported( QObject::tr(
                "The code to support rigid body dynamics of partial molecules "
                "has yet to be written."), CODELOC );
}

/** Calculate the forces and torques */
bool RBWorkspace::calculateForces(const Symbol &nrg_component)
{
    const int nbeads = bead_coordinates.count();
    
    if (nbeads == 0)
        //no beads, so no need to calculate forces
        return false;

    else if (not IntegratorWorkspace::calculateForces(nrg_component))
        return false;

    //forces have changed on the atoms, so recalculate all of
    //the forces and torques on all of the beads

    const MoleculeGroup &molgroup = moleculeGroup();
    const int nmols = molgroup.nMolecules();
    
    Vector *bead_forces_array = bead_forces.data();
    Vector *bead_torques_array = bead_torques.data();

    const Matrix *bead_to_world_array = bead_to_world.constData();
    const Quaternion *bead_orients_array = bead_orientations.constData();

    const QPair< qint32,QVector<qint32> > *atoms_to_beads_array 
                                                    = atoms_to_beads.constData();
    const QVector<Vector> *int_coords_array = atom_int_coords.constData();
    
    const ForceTable &forcetable = forceTable();
    
    for (int i=0; i<nbeads; ++i)
    {
        bead_forces_array[i] = Vector(0);
        bead_torques_array[i] = Vector(0);
    }
    
    for (int i=0; i<nmols; ++i)
    {
        const QPair< qint32,QVector<qint32> > &beading = atoms_to_beads_array[i];

        if (beading.first == -1)
            //there are no beads for this molecule, so no forces or torques
            continue;

        MolNum molnum = molgroup.molNumAt(i);

        ::calculateForces(molgroup[molnum], forcetable.getTable(molnum), 
                          int_coords_array[i], beading.second, 
                          bead_to_world_array + beading.first,
                          bead_orients_array + beading.first,
                          bead_forces_array + beading.first,
                          bead_torques_array + beading.first);
    }
    
    return true;
}

/** Construct an empty workspace */
RBWorkspace::RBWorkspace(const PropertyMap &map)
            : ConcreteProperty<RBWorkspace,IntegratorWorkspace>(map)
{}

/** Construct a workspace for the passed molecule group */
RBWorkspace::RBWorkspace(const MoleculeGroup &molgroup, const PropertyMap &map)
            : ConcreteProperty<RBWorkspace,IntegratorWorkspace>(molgroup, map)
{
    this->rebuildFromScratch();
}

/** Copy constructor */
RBWorkspace::RBWorkspace(const RBWorkspace &other)
            : ConcreteProperty<RBWorkspace,IntegratorWorkspace>(other),
              atom_int_coords(other.atom_int_coords),
              atoms_to_beads(other.atoms_to_beads),
              bead_coordinates(other.bead_coordinates),
              bead_to_world(other.bead_to_world),
              bead_orientations(other.bead_orientations),
              bead_linear_momenta(other.bead_linear_momenta),
              bead_angular_momenta(other.bead_angular_momenta),
              bead_forces(other.bead_forces),
              bead_torques(other.bead_torques),
              bead_masses(other.bead_masses),
              bead_inertia(other.bead_inertia),
              vel_generator(other.vel_generator)
{}

/** Destructor */
RBWorkspace::~RBWorkspace()
{}

/** Copy assignment operator */
RBWorkspace& RBWorkspace::operator=(const RBWorkspace &other)
{
    if (this != &other)
    {
        IntegratorWorkspace::operator=(other);
     
        atom_int_coords = other.atom_int_coords;
        atoms_to_beads = other.atoms_to_beads;
        bead_coordinates = other.bead_coordinates;
        bead_to_world = other.bead_to_world;
        bead_orientations = other.bead_orientations;
        bead_linear_momenta = other.bead_linear_momenta;
        bead_angular_momenta = other.bead_angular_momenta;
        bead_forces = other.bead_forces;
        bead_torques = other.bead_torques;
        bead_masses = other.bead_masses;
        bead_inertia = other.bead_inertia;
        vel_generator = other.vel_generator;
    }
    
    return *this;
}

/** Comparison operator */
bool RBWorkspace::operator==(const RBWorkspace &other) const
{
    return this == &other or 
           (IntegratorWorkspace::operator==(other) and
            atom_int_coords == other.atom_int_coords and 
            atoms_to_beads == other.atoms_to_beads and
            bead_coordinates == other.bead_coordinates and
            bead_to_world == other.bead_to_world and
            bead_orientations == other.bead_orientations and
            bead_linear_momenta == other.bead_linear_momenta and
            bead_angular_momenta == other.bead_angular_momenta and
            bead_masses == other.bead_masses and
            bead_inertia == other.bead_inertia and
            vel_generator == other.vel_generator);
            
    //don't need forces or torques as these are implied by the
    //forcetable in IntegratorWorkspace
}

/** Comparison operator */
bool RBWorkspace::operator!=(const RBWorkspace &other) const
{
    return not RBWorkspace::operator==(other);
}

const char* RBWorkspace::typeName()
{
    return QMetaType::typeName( qMetaTypeId<RBWorkspace>() );
}

/** Return the kinetic energy of all of the molecules being integrated */
MolarEnergy RBWorkspace::kineticEnergy() const
{
    //sum together the linear kinetic energy of the beads...
    double nrg = 0;
    
    int nbeads = bead_linear_momenta.count();
    
    const Vector *p = bead_linear_momenta.constData();
    const double *m = bead_masses.constData();
    
    for (int i=0; i<nbeads; ++i)
    {
        if (m[i] != 0)
            nrg += p[i].length2() / (2 * m[i]);
    }
    
    //now the angular kinetic energy
    const Vector *q = bead_angular_momenta.constData();
    const Vector *I = bead_inertia.constData();
    
    for (int i=0; i<nbeads; ++i)
    {
        if (I[i].x() != 0)
            nrg += pow_2(q[i].x()) / (2 * I[i].x());
            
        if (I[i].y() != 0)
            nrg += pow_2(q[i].y()) / (2 * I[i].y());
            
        if (I[i].z() != 0)
            nrg += pow_2(q[i].z()) / (2 * I[i].z());
    }
    
    return MolarEnergy(nrg);
}

/** Return the kinetic energy of the molecule with number 'molnum'

    \throw SireMol::missing_molecule
*/
MolarEnergy RBWorkspace::kineticEnergy(MolNum molnum) const
{
    throw SireError::incomplete_code( QObject::tr("Need to write!"), CODELOC );
    return MolarEnergy(0);
}

/** Return the kinetic energy of the atoms in the view 'molview'

    \throw SireMol::missing_molecule
*/
MolarEnergy RBWorkspace::kineticEnergy(const MoleculeView &molview) const
{
    throw SireError::incomplete_code( QObject::tr("Need to write!"), CODELOC );
    return MolarEnergy(0);
}

/** Return the number of rigid body beads to be integrated */
int RBWorkspace::nBeads() const
{
    return bead_coordinates.count();
}

/** Return the number of atoms in the ith bead */
int RBWorkspace::nAtoms(int i) const
{
    return atom_int_coords.constData()[i].count();
}

/** Set the system to be integrated */
bool RBWorkspace::setSystem(const System &system)
{
    if (IntegratorWorkspace::setSystem(system))
    {
        this->rebuildFromScratch();
        return true;
    }
    else
        return false;
}

/** Return the array of coordinates of the center of masses of the beads */
Vector* RBWorkspace::beadCoordsArray()
{
    return bead_coordinates.data();
}

/** Return the array of orientations of the beads (this is the rotation
    to be applied to the matrix that maps from the cartesian frame
    to the internal principle inertia frame) */
Quaternion* RBWorkspace::beadOrientationArray()
{
    return bead_orientations.data();
}

/** Return the array of bead linear momenta */
Vector* RBWorkspace::beadLinearMomentaArray()
{
    return bead_linear_momenta.data();
}

/** Return the array of bead angular momenta */
Vector* RBWorkspace::beadAngularMomentaArray()
{
    return bead_angular_momenta.data();
}

/** Return the array of bead masses */
const double* RBWorkspace::beadMassesArray() const
{
    return bead_masses.constData();
}

/** Retunr the array of bead inertia */
const Vector* RBWorkspace::beadInertiasArray() const
{
    return bead_inertia.constData();
}
    
/** Return the array of forces acting on the center of mass
    of each bead */
const Vector* RBWorkspace::beadForcesArray() const
{
    return bead_forces.data();
}

/** Return the array of bead torques */
const Vector* RBWorkspace::beadTorquesArray() const
{
    return bead_torques.data();
}

/** Regenerate all of the linear and angular velocities using the passed generator */
void RBWorkspace::regenerateVelocities(const VelocityGenerator &generator)
{
    throw SireError::incomplete_code( QObject::tr("Need to write!"), CODELOC );
}

static AtomCoords updateCoordinates(const ViewsOfMol &mol,
                                    const PropertyName &coords_property,
                                    const QVector<qint32> &beading,
                                    const QVector<Vector> &int_coords,
                                    const Vector *bead_coords,
                                    const Matrix *beads_to_world,
                                    const Quaternion *bead_orients)
{
    AtomCoords coords = mol.data().property(coords_property).asA<AtomCoords>();

    if (mol.selectedAll())
    {
        //calculate the world-cartesian coordinates of the beads
        //from the current orientation and internal coordinates
        const int nats = int_coords.count();
        const Vector *int_coords_array = int_coords.constData();
        QVector<Vector> new_coords(nats);
        Vector *new_coords_array = new_coords.data();
        
        if (beading.isEmpty())
        {
            //just one bead for the whole molecule
            const Vector &com = bead_coords[0];
            Matrix orient = bead_orients[0].toMatrix() * beads_to_world[0];
        
            for (int i=0; i<nats; ++i)
            {
                new_coords_array[i] = com + (orient * int_coords_array[i]);
            }
        }
        else
        {
            //lots of beads in the molecule
            qint32 last_idx = -1;
            Vector com;
            Matrix orient;
            
            const Vector *old_coords_array = coords.array().constCoordsData();
            
            const qint32 *beading_array = beading.constData();
            
            for (int i=0; i<nats; ++i)
            {
                qint32 bead_idx = beading_array[i];
                
                if (bead_idx == -1)
                {
                    //this atom is not part of a bead
                    new_coords_array[i] = old_coords_array[i];
                    continue;
                }
                
                if (bead_idx != last_idx)
                {
                    last_idx = bead_idx;
                    com = bead_coords[bead_idx];
                    orient = bead_orients[bead_idx].toMatrix() * beads_to_world[bead_idx];                     
                }
                
                new_coords_array[i] = com + (orient * int_coords_array[i]);
            }
        }
            
        coords.copyFrom(new_coords);
    }
    else
        throw SireError::unsupported( QObject::tr(
                "The code to support rigid body dynamics of partial "
                "molecules has yet to be written."), CODELOC );
                
    return coords;
}

/** Commit the coordinates back to the system. This maps the bead coordinates
    and orientations back to atomic coordinates and position and 
    updates the system with these */
void RBWorkspace::commitCoordinates()
{
    int nmols = atom_int_coords.count();
    
    const MoleculeGroup &molgroup = moleculeGroup();
    
    BOOST_ASSERT( molgroup.nMolecules() == nmols );
    
    const QVector<Vector> *int_coords_array = atom_int_coords.constData();
    const QPair< qint32,QVector<qint32> > 
                *atoms_to_beads_array = atoms_to_beads.constData();
    
    const Vector *bead_coords_array = bead_coordinates.constData();
    const Quaternion *bead_orients_array = bead_orientations.constData();
    const Matrix *bead_to_world_array = bead_to_world.constData();
    
    PropertyName coords_property = coordinatesProperty();
    
    Molecules changed_mols;
    changed_mols.reserve(nmols);
    
    for (int i=0; i<nmols; ++i)
    {
        const QPair< qint32,QVector<qint32> > &beading = atoms_to_beads_array[i];

        if (beading.first == -1)
            //there are no beads for this molecule, so no new coordinates
            continue;

        MolNum molnum = molgroup.molNumAt(i);
        const ViewsOfMol &mol = molgroup[molnum];

        AtomCoords new_coords = ::updateCoordinates(mol, coords_property,  
                                                    beading.second,
                                                    int_coords_array[i],
                                                    bead_coords_array + beading.first,
                                                    bead_to_world_array + beading.first,
                                                    bead_orients_array + beading.first);

        changed_mols.add( mol.molecule().edit()
                             .setProperty(coords_property,new_coords).commit() );
    }
    
    IntegratorWorkspace::pvt_update(changed_mols);
}

/** Commit the linear and angular velocities back to the system. This saves
    the velocities as bead properties */
void RBWorkspace::commitVelocities()
{
    //throw SireError::incomplete_code( QObject::tr("Need to write!"), CODELOC );
}

/** Commit both the coordinates and velocities - this performs the 
    equivalent of commitCoordinates() and commitVelocities() in
    a single call */
void RBWorkspace::commitCoordinatesAndVelocities()
{
    RBWorkspace::commitCoordinates();
    //throw SireError::incomplete_code( QObject::tr("Need to write!"), CODELOC );
}

/** This internal function is called whenever a property is changed. 
    This is used to see if the data has to be regenerated */
void RBWorkspace::changedProperty(const QString &property)
{
    this->rebuildFromScratch();
}
