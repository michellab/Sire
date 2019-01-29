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

#include "integrator.h"
#include "integratorworkspace.h"
#include "ensemble.h"

#include "SireFF/forcetable.h"
#include "SireFF/forcefields.h"

#include "SireSystem/system.h"

#include "SireCAS/symbol.h"

#include "SireMaths/rangenerator.h"

#include "SireMol/moleculegroup.h"
#include "SireMol/moleculeview.h"

#include "SireStream/datastream.h"

using namespace SireMove;
using namespace SireFF;
using namespace SireSystem;
using namespace SireCAS;
using namespace SireBase;
using namespace SireStream;
using namespace SireUnits::Dimension;

////////////
//////////// Implementation of Integrator
////////////

static const RegisterMetaType<Integrator> r_integrator(MAGIC_ONLY,
                                               "SireMove::Integrator");

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const Integrator &integrator)
{
    writeHeader(ds, r_integrator, 1);
    
    ds << static_cast<const Property&>(integrator);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, Integrator &integrator)
{
    VersionID v = readHeader(ds, r_integrator);
    
    if (v == 1)
    {
        ds >> static_cast<Property&>(integrator);
    }
    else
        throw version_error(v, "1", r_integrator, CODELOC);
        
    return ds;
}

/** Constructor */
Integrator::Integrator() : Property()
{}

/** Copy constructor */
Integrator::Integrator(const Integrator &other) : Property(other)
{}

/** Destructor */
Integrator::~Integrator()
{}

/** Copy assignment operator */
Integrator& Integrator::operator=(const Integrator &other)
{
    return *this;
}

/** Comparison operator */
bool Integrator::operator==(const Integrator &other) const
{
    return true;
}

/** Comparison operator */
bool Integrator::operator!=(const Integrator &other) const
{
    return false;
}

Q_GLOBAL_STATIC( NullIntegrator, getNullIntegrator );

/** Return a NullIntegrator */
const NullIntegrator& Integrator::null()
{
    return *(getNullIntegrator());
}

////////////
//////////// Implementation of NullIntegrator
////////////

static const RegisterMetaType<NullIntegrator> r_nullint;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const NullIntegrator &nullint)
{
    writeHeader(ds, r_nullint, 1);
    
    ds << static_cast<const Integrator&>(nullint);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, NullIntegrator &nullint)
{
    VersionID v = readHeader(ds, r_nullint);
    
    if (v == 1)
    {
        ds >> static_cast<Integrator&>(nullint);
    }
    else
        throw version_error(v, "1", r_nullint, CODELOC);
        
    return ds;
}

/** Constructor */
NullIntegrator::NullIntegrator() : ConcreteProperty<NullIntegrator,Integrator>()
{}

/** Copy constructor */
NullIntegrator::NullIntegrator(const NullIntegrator &other)
               : ConcreteProperty<NullIntegrator,Integrator>(other)
{}

/** Destructor */
NullIntegrator::~NullIntegrator()
{}

/** Copy assignment operator */
NullIntegrator& NullIntegrator::operator=(const NullIntegrator &other)
{
    Integrator::operator=(other);
    return *this;
}

/** Comparison operator - all NullIntegrators are the same */
bool NullIntegrator::operator==(const NullIntegrator &other) const
{
    return Integrator::operator==(other);
}

/** Comparison operator - all NullIntegrators are the same */
bool NullIntegrator::operator!=(const NullIntegrator &other) const
{
    return Integrator::operator!=(other);
}

/** Return a string representation of this integrator */
QString NullIntegrator::toString() const
{
    return QObject::tr("NullIntegrator");
}

/** Return the ensemble of this integrator */
Ensemble NullIntegrator::ensemble() const
{
    return Ensemble::NVE();
}

/** Return whether or not this integrator is time-reversible */
bool NullIntegrator::isTimeReversible() const
{
    return true;
}

/** The null integrator does nothing */
void NullIntegrator::integrate(IntegratorWorkspace&, 
                               const Symbol&, 
                               SireUnits::Dimension::Time,
                               int, bool)
{}

/** This returns a null workspace */
IntegratorWorkspacePtr NullIntegrator::createWorkspace(const PropertyMap&) const
{
    return IntegratorWorkspacePtr();
}

/** This returns a null workspace */
IntegratorWorkspacePtr NullIntegrator::createWorkspace(const MoleculeGroup&,
                                                       const PropertyMap&) const
{
    return IntegratorWorkspacePtr();
}

const char* NullIntegrator::typeName()
{
    return QMetaType::typeName( qMetaTypeId<NullIntegrator>() );
}
