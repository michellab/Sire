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

#ifndef SIREMOVE_INTEGRATOR_H
#define SIREMOVE_INTEGRATOR_H

#include "SireBase/property.h"
#include "SireBase/propertymap.h"

#include "SireUnits/dimensions.h"

#include "SireMol/moleculegroup.h"

#include "integratorworkspace.h"

SIRE_BEGIN_HEADER

namespace SireMove
{
class Integrator;
class NullIntegrator;
}

SIREMOVE_EXPORT QDataStream& operator<<(QDataStream&, const SireMove::Integrator&);
SIREMOVE_EXPORT QDataStream& operator>>(QDataStream&, SireMove::Integrator&);

SIREMOVE_EXPORT QDataStream& operator<<(QDataStream&, const SireMove::NullIntegrator&);
SIREMOVE_EXPORT QDataStream& operator>>(QDataStream&, SireMove::NullIntegrator&);

namespace SireFF
{
class ForceTable;
}

namespace SireMol
{
class MoleculeGroup;
class MoleculeView;
}

namespace SireMaths
{
class RanGenerator;
}

namespace SireSystem
{
class System;
}

namespace SireCAS
{
class Symbol;
}

namespace SireMove
{

class Ensemble;

using SireMol::MoleculeGroup;
using SireMol::MoleculeView;

using SireSystem::System;

using SireFF::ForceTable;

using SireCAS::Symbol;

using SireBase::PropertyMap;

using SireMaths::RanGenerator;

/** This is the virtual base class of all dynamics integrators. An integrator
    is the kernel used to advance the coordinates of the system from one
    timestep to the next
    
    @author Christopher Woods
*/
class SIREMOVE_EXPORT Integrator : public SireBase::Property
{

friend SIREMOVE_EXPORT QDataStream& ::operator<<(QDataStream&, const Integrator&);
friend SIREMOVE_EXPORT QDataStream& ::operator>>(QDataStream&, Integrator&);

public:
    Integrator();
    
    Integrator(const Integrator &other);
    
    virtual ~Integrator();
    
    static const char* typeName()
    {
        return "SireMove::Integrator";
    }

    virtual Integrator* clone() const=0;
    
    static const NullIntegrator& null();
    
    virtual QString toString() const=0;
    
    virtual Ensemble ensemble() const=0;
    
    virtual bool isTimeReversible() const=0;
    
    virtual void integrate(IntegratorWorkspace &workspace, 
                           const Symbol &nrg_component,
                           SireUnits::Dimension::Time timestep,
                           int nmoves, bool record_stats)=0;
    
    virtual IntegratorWorkspacePtr createWorkspace(
                                        const PropertyMap &map = PropertyMap()) const=0;
    
    virtual IntegratorWorkspacePtr 
                        createWorkspace(const MoleculeGroup &molgroup,
                                        const PropertyMap &map = PropertyMap()) const=0;
    
protected:
    Integrator& operator=(const Integrator &other);
    
    bool operator==(const Integrator &other) const;
    bool operator!=(const Integrator &other) const;
};

/** This class holds a null integrator, which doesn't advance anything anywhere

    @author Christopher Woods
*/
class SIREMOVE_EXPORT NullIntegrator 
            : public SireBase::ConcreteProperty<NullIntegrator,Integrator>
{

friend SIREMOVE_EXPORT QDataStream& ::operator<<(QDataStream&, const NullIntegrator&);
friend SIREMOVE_EXPORT QDataStream& ::operator>>(QDataStream&, NullIntegrator&);

public:
    NullIntegrator();
    
    NullIntegrator(const NullIntegrator &other);
    
    ~NullIntegrator();
    
    NullIntegrator& operator=(const NullIntegrator &other);
    
    static const char* typeName();
    
    bool operator==(const NullIntegrator &other) const;
    bool operator!=(const NullIntegrator &other) const;
    
    Ensemble ensemble() const;
    
    bool isTimeReversible() const;
    
    QString toString() const;
    
    void integrate(IntegratorWorkspace &workspace, const Symbol &nrg_component, 
                   SireUnits::Dimension::Time timestep, 
                   int nmoves, bool record_stats);
    
    IntegratorWorkspacePtr createWorkspace(const PropertyMap &map = PropertyMap()) const;
    
    IntegratorWorkspacePtr createWorkspace(const MoleculeGroup &molgroup,
                                           const PropertyMap &map = PropertyMap()) const;
};

typedef SireBase::PropPtr<Integrator> IntegratorPtr;

}

Q_DECLARE_METATYPE( SireMove::NullIntegrator )

SIRE_EXPOSE_CLASS( SireMove::Integrator )
SIRE_EXPOSE_CLASS( SireMove::NullIntegrator )

SIRE_EXPOSE_PROPERTY( SireMove::IntegratorPtr, SireMove::Integrator )

SIRE_END_HEADER

#endif
