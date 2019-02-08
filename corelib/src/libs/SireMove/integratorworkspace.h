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

#ifndef SIREMOVE_INTEGRATORWORKSPACE_H
#define SIREMOVE_INTEGRATORWORKSPACE_H

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

#include "SireSystem/system.h"

#include "velocitygenerator.h"

SIRE_BEGIN_HEADER

namespace SireMove
{
class IntegratorWorkspace;
class NullIntegratorWorkspace;
class AtomicVelocityWorkspace;
}

QDataStream& operator<<(QDataStream&, const SireMove::IntegratorWorkspace&);
QDataStream& operator>>(QDataStream&, SireMove::IntegratorWorkspace&);

QDataStream& operator<<(QDataStream&, const SireMove::NullIntegratorWorkspace&);
QDataStream& operator>>(QDataStream&, SireMove::NullIntegratorWorkspace&);

QDataStream& operator<<(QDataStream&, const SireMove::AtomicVelocityWorkspace&);
QDataStream& operator>>(QDataStream&, SireMove::AtomicVelocityWorkspace&);

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
class SIREMOVE_EXPORT IntegratorWorkspace : public SireBase::Property
{

friend SIREMOVE_EXPORT QDataStream& ::operator<<(QDataStream&, const IntegratorWorkspace&);
friend SIREMOVE_EXPORT QDataStream& ::operator>>(QDataStream&, IntegratorWorkspace&);

public:
    IntegratorWorkspace(const PropertyMap &map = PropertyMap());
    IntegratorWorkspace(const MoleculeGroup &molgroup,
                        const PropertyMap &map = PropertyMap());
    
    IntegratorWorkspace(const IntegratorWorkspace &other);
    
    virtual ~IntegratorWorkspace();
    
    static const char* typeName()
    {
        return "SireBase::IntegratorWorkspace";
    }

    virtual IntegratorWorkspace* clone() const=0;

    const MoleculeGroup& moleculeGroup() const;

    const ForceTable& forceTable() const;

    virtual bool setSystem(const System &system);

    const System& system() const;

    System& nonConstsystem();
          
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

    bool forcesNeedCalculating(const Symbol &nrg_component) const;
    
    void mustNowRecalculateFromScratch();

    void collectStatistics();

    virtual void regenerateVelocities(const VelocityGenerator &generator)=0;

    virtual SireUnits::Dimension::MolarEnergy kineticEnergy() const=0;
    virtual SireUnits::Dimension::MolarEnergy
                            kineticEnergy(MolNum molnum) const=0;
    virtual SireUnits::Dimension::MolarEnergy 
                            kineticEnergy(const MoleculeView &molview) const=0;
    
    static const NullIntegratorWorkspace& null();

protected:
    IntegratorWorkspace& operator=(const IntegratorWorkspace &other);
    
    bool operator==(const IntegratorWorkspace &other) const;
    bool operator!=(const IntegratorWorkspace &other) const;

    virtual void changedProperty(const QString &property);

    void pvt_update(const Molecules &changed_mols);

private:
    /** The system being integrated */
    System sys;

    /** The molecule group containing the molecules being integrated */
    MolGroupPtr molgroup;

    /** The current forces acting on the molecules */
    ForceTable molforces;
    
    /** The energy component used when we last got the forces */
    SireCAS::Symbol last_nrg_component;

    /** The property map used to find the sources of required properties */
    PropertyMap map;

    /** Whether or not the forces need to be recalculated */
    bool need_new_forces;
};

/** This is the null integrator workspace */
class SIREMOVE_EXPORT NullIntegratorWorkspace
        : public SireBase::ConcreteProperty<NullIntegratorWorkspace,IntegratorWorkspace>
{

friend SIREMOVE_EXPORT QDataStream& ::operator<<(QDataStream&, const NullIntegratorWorkspace&);
friend SIREMOVE_EXPORT QDataStream& ::operator>>(QDataStream&, NullIntegratorWorkspace&);

public:
    NullIntegratorWorkspace();
    NullIntegratorWorkspace(const NullIntegratorWorkspace &other);
    
    ~NullIntegratorWorkspace();
    
    NullIntegratorWorkspace& operator=(const NullIntegratorWorkspace &other);
    
    bool operator==(const NullIntegratorWorkspace &other) const;
    bool operator!=(const NullIntegratorWorkspace &other) const;

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
class SIREMOVE_EXPORT AtomicVelocityWorkspace
       : public SireBase::ConcreteProperty<AtomicVelocityWorkspace,IntegratorWorkspace>
{

friend SIREMOVE_EXPORT QDataStream& ::operator<<(QDataStream&, const AtomicVelocityWorkspace&);
friend SIREMOVE_EXPORT QDataStream& ::operator>>(QDataStream&, AtomicVelocityWorkspace&);

public:
    AtomicVelocityWorkspace(const PropertyMap &map = PropertyMap());
    AtomicVelocityWorkspace(const MoleculeGroup &molgroup,
                            const PropertyMap &map = PropertyMap());
    
    AtomicVelocityWorkspace(const AtomicVelocityWorkspace &other);
    
    ~AtomicVelocityWorkspace();

    AtomicVelocityWorkspace& operator=(const AtomicVelocityWorkspace &other);
    
    bool operator==(const AtomicVelocityWorkspace &other) const;
    bool operator!=(const AtomicVelocityWorkspace &other) const;
    
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
    const Vector* momentaArray(int i) const;

    const double* massArray(int i) const;
    
    const Vector* constCoordsArray(int i) const;
    const Vector* constForceArray(int i) const;
    const Vector* constMomentaArray(int i) const;
    
    const double* constMassArray(int i) const;
    
    bool calculateForces(const Symbol &nrg_component);
    
    bool setSystem(const System &system);

    void regenerateVelocities(const VelocityGenerator &generator);
    
    void commitCoordinates();
    void commitVelocities();
    
    void commitCoordinatesAndVelocities();

    void commitBufferedCoordinatesAndVelocities(  QVector < QVector< QVector < Vector > > > &buffered_coords);

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
    
    /** All of the atom masses */
    QVector< QVector<double> > atom_masses;
    
    /** The generator used to get the initial velocities */
    VelGenPtr vel_generator;
};

typedef SireBase::PropPtr<IntegratorWorkspace> IntegratorWorkspacePtr;

}

Q_DECLARE_METATYPE( SireMove::NullIntegratorWorkspace )
Q_DECLARE_METATYPE( SireMove::AtomicVelocityWorkspace )

SIRE_END_HEADER

#endif
