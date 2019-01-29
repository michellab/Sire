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

#ifndef SIREMOVE_DLMRIGIDBODY_H
#define SIREMOVE_DLMRIGIDBODY_H

#include "integrator.h"

SIRE_BEGIN_HEADER

namespace SireMove
{
class DLMRigidBody;
}

SIREMOVE_EXPORT QDataStream& operator<<(QDataStream&, const SireMove::DLMRigidBody&);
SIREMOVE_EXPORT QDataStream& operator>>(QDataStream&, SireMove::DLMRigidBody&);

namespace SireMove
{

/** This class implements the rigid body dynamics integrator
    as described in;
    
    "Symplectic splitting methods for rigid body molecular dynamics"
    Andreas Dullweber, Benedict Leimkuhler, Robert McLachlan
    
    J. Chem. Phys., 107:15, 5840-5851, 1997
    
    This is a time-reversible, symplectic integrator for 
    general rigid bodies
    
    @author Christopher Woods
*/
class SIREMOVE_EXPORT DLMRigidBody 
        : public SireBase::ConcreteProperty<DLMRigidBody,Integrator>
{

friend QDataStream& ::operator<<(QDataStream&, const DLMRigidBody&);
friend QDataStream& ::operator>>(QDataStream&, DLMRigidBody&);

public:
    DLMRigidBody(bool frequent_save_velocities = false);
    
    DLMRigidBody(const DLMRigidBody &other);
    
    ~DLMRigidBody();
    
    DLMRigidBody& operator=(const DLMRigidBody &other);
    
    bool operator==(const DLMRigidBody &other) const;
    bool operator!=(const DLMRigidBody &other) const;
    
    static const char* typeName();
    
    QString toString() const;
    
    Ensemble ensemble() const;
    
    bool isTimeReversible() const;
    
    void integrate(IntegratorWorkspace &workspace,
                   const Symbol &nrg_component, 
                   SireUnits::Dimension::Time timestep,
                   int nmoves, bool record_stats);

    IntegratorWorkspacePtr createWorkspace(const PropertyMap &map = PropertyMap()) const;
    IntegratorWorkspacePtr createWorkspace(const MoleculeGroup &molgroup,
                                           const PropertyMap &map = PropertyMap()) const;

private:
    /** Whether or not to save velocities at every step
        (or just at the end of all steps) */
    bool frequent_save_velocities;
};

}

Q_DECLARE_METATYPE( SireMove::DLMRigidBody )

SIRE_EXPOSE_CLASS( SireMove::DLMRigidBody )

SIRE_END_HEADER

#endif
