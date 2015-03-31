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

#ifndef SIREMOVE_INTEGRATORWORKSPACEJM_H
#define SIREMOVE_INTEGRATORWORKSPACEJM_H

#include <QUuid>

#include "SireBase/property.h"
#include "SireBase/properties.h"
#include "SireBase/majorminorversion.h"

#include "SireCAS/symbol.h"

#include "SireMol/moleculegroup.h"
#include "SireMol/atomvelocities.h"
#include "SireMol/atomforces.h"
#include "SireMol/atommasses.h"

#include "SireFF/forcetable.h"
#include "SireFF/energytable.h"

#include "SireSystem/system.h"

#include "velocitygenerator.h"

SIRE_BEGIN_HEADER

namespace SireMove
{
class IntegratorWorkspaceJM;
class NullIntegratorWorkspaceJM;
class AtomicVelocityWorkspaceJM;
}

QDataStream& operator<<(QDataStream&, const SireMove::IntegratorWorkspaceJM&);
QDataStream& operator>>(QDataStream&, SireMove::IntegratorWorkspaceJM&);

QDataStream& operator<<(QDataStream&, const SireMove::NullIntegratorWorkspaceJM&);
QDataStream& operator>>(QDataStream&, SireMove::NullIntegratorWorkspaceJM&);

QDataStream& operator<<(QDataStream&, const SireMove::AtomicVelocityWorkspaceJM&);
QDataStream& operator>>(QDataStream&, SireMove::AtomicVelocityWorkspaceJM&);

namespace SireMol
{
class MoleculeView;
}

namespace SireMove
{

class VelocityGenerator;

using SireMol::MoleculeView;
using SireMol::MoleculeGroup;
using SireMol::MolGroupPtr;
using SireMol::MGNum;
using SireMol::MolNum;
using SireMol::MolID;
using SireMol::AtomForces;
using SireMol::AtomMasses;
using SireMol::AtomVelocities;
using SireMol::Velocity3D;
using SireMol::Molecules;

using SireSystem::System;

using SireFF::EnergyTable;
using SireFF::ForceTable;

using SireMaths::Vector;

using SireCAS::Symbol;

using SireBase::PropertyMap;
using SireBase::PropertyName;

/** This is the base class of the workspaces which are used to 
    hold the intermediate values used when integrating the 
    dynamics of a system
    
    @author Christopher Woods
*/
class SIREMOVE_EXPORT IntegratorWorkspaceJM : public SireBase::Property
{

friend QDataStream& ::operator<<(QDataStream&, const IntegratorWorkspaceJM&);
friend QDataStream& ::operator>>(QDataStream&, IntegratorWorkspaceJM&);

public:
    IntegratorWorkspaceJM(const PropertyMap &map = PropertyMap());
    IntegratorWorkspaceJM(const MoleculeGroup &molgroup,
                        const PropertyMap &map = PropertyMap());
    
    IntegratorWorkspaceJM(const IntegratorWorkspaceJM &other);
    
    virtual ~IntegratorWorkspaceJM();
    
    static const char* typeName()
    {
        return "SireBase::IntegratorWorkspaceJM";
    }

    virtual IntegratorWorkspaceJM* clone() const=0;

    const MoleculeGroup& moleculeGroup() const;

    const ForceTable& forceTable() const;
   
    const EnergyTable& energyTable() const;

    virtual bool setSystem(const System &system);

    const System& system() const;

    const PropertyMap& propertyMap() const;
    
    virtual void setPropertyMap(const PropertyMap &map);
    
    virtual void setGenerator(const RanGenerator &generator);

    virtual void setCoordinatesProperty(const PropertyName &source);
    virtual void setSpaceProperty(const PropertyName &source);
    virtual void setVelocitiesProperty(const PropertyName &source);
    virtual void setMassesProperty(const PropertyName &source);
    virtual void setElementsProperty(const PropertyName &source);

    virtual void setVelocityGeneratorProperty(const PropertyName &source);

    PropertyName coordinatesProperty() const;
    PropertyName spaceProperty() const;
    PropertyName velocitiesProperty() const;
    PropertyName massesProperty() const;
    PropertyName elementsProperty() const;
    PropertyName velocityGeneratorProperty() const;

    virtual bool calculateForces(const Symbol &nrg_component);

    //    virtual bool calculateEnergies(const Symbol &nrg_component);

    bool forcesNeedCalculating(const Symbol &nrg_component) const;
    
    void mustNowRecalculateFromScratch();

    void collectStatistics();

    virtual void regenerateVelocities(const VelocityGenerator &generator)=0;

    virtual SireUnits::Dimension::MolarEnergy kineticEnergy() const=0;
    virtual SireUnits::Dimension::MolarEnergy
                            kineticEnergy(MolNum molnum) const=0;
    virtual SireUnits::Dimension::MolarEnergy 
                            kineticEnergy(const MoleculeView &molview) const=0;
    
    static const NullIntegratorWorkspaceJM& null();

protected:
    IntegratorWorkspaceJM& operator=(const IntegratorWorkspaceJM &other);
    
    bool operator==(const IntegratorWorkspaceJM &other) const;
    bool operator!=(const IntegratorWorkspaceJM &other) const;

    virtual void changedProperty(const QString &property);

    void pvt_update(const Molecules &changed_mols);

private:
    /** The system being integrated */
    System sys;

    /** The molecule group containing the molecules being integrated */
    MolGroupPtr molgroup;

    /** The current forces acting on the molecules */
    ForceTable molforces;

    /** The current energies of the molecules */
    EnergyTable molenergies;

    /** The energy component used when we last got the forces */
    SireCAS::Symbol last_nrg_component;

    /** The property map used to find the sources of required properties */
    PropertyMap map;

    /** Whether or not the forces need to be recalculated */
    bool need_new_forces;



};

/** This is the null integrator workspace */
class SIREMOVE_EXPORT NullIntegratorWorkspaceJM
        : public SireBase::ConcreteProperty<NullIntegratorWorkspaceJM,IntegratorWorkspaceJM>
{

friend QDataStream& ::operator<<(QDataStream&, const NullIntegratorWorkspaceJM&);
friend QDataStream& ::operator>>(QDataStream&, NullIntegratorWorkspaceJM&);

public:
    NullIntegratorWorkspaceJM();
    NullIntegratorWorkspaceJM(const NullIntegratorWorkspaceJM &other);
    
    ~NullIntegratorWorkspaceJM();
    
    NullIntegratorWorkspaceJM& operator=(const NullIntegratorWorkspaceJM &other);
    
    bool operator==(const NullIntegratorWorkspaceJM &other) const;
    bool operator!=(const NullIntegratorWorkspaceJM &other) const;

    static const char* typeName();

    void regenerateVelocities(const VelocityGenerator &generator);

    SireUnits::Dimension::MolarEnergy kineticEnergy() const;
    SireUnits::Dimension::MolarEnergy kineticEnergy(MolNum molnum) const;
    SireUnits::Dimension::MolarEnergy kineticEnergy(const MoleculeView &molview) const;
};

/** This class provides a workspace for integrators that make use
    of atomic forces and velocities 
    
    @author Christopher Woods
*/
class SIREMOVE_EXPORT AtomicVelocityWorkspaceJM
       : public SireBase::ConcreteProperty<AtomicVelocityWorkspaceJM,IntegratorWorkspaceJM>
{

friend QDataStream& ::operator<<(QDataStream&, const AtomicVelocityWorkspaceJM&);
friend QDataStream& ::operator>>(QDataStream&, AtomicVelocityWorkspaceJM&);

public:
    AtomicVelocityWorkspaceJM(const PropertyMap &map = PropertyMap());
    AtomicVelocityWorkspaceJM(const MoleculeGroup &molgroup,
                            const PropertyMap &map = PropertyMap());
    
    AtomicVelocityWorkspaceJM(const AtomicVelocityWorkspaceJM &other);
    
    ~AtomicVelocityWorkspaceJM();

    AtomicVelocityWorkspaceJM& operator=(const AtomicVelocityWorkspaceJM &other);
    
    bool operator==(const AtomicVelocityWorkspaceJM &other) const;
    bool operator!=(const AtomicVelocityWorkspaceJM &other) const;
    
    static const char* typeName();
    
    SireUnits::Dimension::MolarEnergy kineticEnergy() const;
    SireUnits::Dimension::MolarEnergy kineticEnergy(MolNum molnum) const;
    SireUnits::Dimension::MolarEnergy kineticEnergy(const MoleculeView &molview) const;

    int nMolecules() const;

    int nAtoms(int i) const;
    
    Vector* coordsArray(int i);
    Vector* momentaArray(int i);

    const Vector* coordsArray(int i) const;
    const Vector* forceArray(int i) const;
    const Vector* energyArray(int i) const;   

    const Vector* momentaArray(int i) const;

    const double* massArray(int i) const;
    
    const Vector* constCoordsArray(int i) const;
    const Vector* constForceArray(int i) const;
    const Vector* constEnergyArray(int i) const;
    const Vector* constMomentaArray(int i) const;
    
    const double* constMassArray(int i) const;
    
    bool calculateForces(const Symbol &nrg_component);
    

    bool setSystem(const System &system);

    void regenerateVelocities(const VelocityGenerator &generator);
    
    void commitCoordinates();
    void commitVelocities();
    
    void commitCoordinatesAndVelocities();

protected:
    void changedProperty(const QString &property);

private:
    void rebuildFromScratch();

    /** All of the atomic coordinates */
    QVector< QVector<Vector> > atom_coords;
    
    /** All of the atomic momenta */
    QVector< QVector<Vector> > atom_momenta;
    
    /** All of the forces for molecules that are not 
        fully selected */
    QVector< QVector<Vector> > atom_forces; 

    /** All of the energies for molecules that are not 
	fully selected */
    QVector< QVector<Vector> > atom_energies; 

    /** All of the atom masses */
    QVector< QVector<double> > atom_masses;
    
    /** The generator used to get the initial velocities */
    VelGenPtr vel_generator;
};

typedef SireBase::PropPtr<IntegratorWorkspaceJM> IntegratorWorkspaceJMPtr;

}

Q_DECLARE_METATYPE( SireMove::NullIntegratorWorkspaceJM )
Q_DECLARE_METATYPE( SireMove::AtomicVelocityWorkspaceJM )

SIRE_END_HEADER

#endif
