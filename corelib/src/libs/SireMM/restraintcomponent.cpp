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

#include "restraintcomponent.h"

#include "SireStream/datastream.h"

using namespace SireMM;
using namespace SireFF;
using namespace SireCAS;
using namespace SireStream;

static const RegisterMetaType<RestraintComponent> r_rest;

/** Serialise to a binary datastream */
QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, const RestraintComponent &rest)
{
    writeHeader(ds, r_rest, 1);
    ds << static_cast<const FFComponent&>(rest);
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, RestraintComponent &rest)
{
    VersionID v = readHeader(ds, r_rest);
    
    if (v == 1)
    {
        ds >> static_cast<FFComponent&>(rest);
    }
    else
        throw version_error(v, "1", r_rest, CODELOC);
        
    return ds;
}

/** Constructor */
RestraintComponent::RestraintComponent(const FFName &ffname)
                   : FFComponent(ffname, QLatin1String("restraint"))
{}

/** Construct using the passed forcefield name and suffix */
RestraintComponent::RestraintComponent(const FFName &ffname, const QString &suffix)
                   : FFComponent(ffname, QString("restraint_{%1}").arg(suffix))
{}

/** Construct from a symbol

    \throw SireError::incompatible_error
*/
RestraintComponent::RestraintComponent(const SireCAS::Symbol &symbol)
                   : FFComponent(symbol, QLatin1String("restraint"))
{}

/** Copy constructor */  
RestraintComponent::RestraintComponent(const RestraintComponent &other)
                   : FFComponent(other)
{}
  
/** Destructor */  
RestraintComponent::~RestraintComponent()
{}

/** Set the restraint component of the energy in the forcefield 'ff'
    to equal to the passed RestraintEnergy */
void RestraintComponent::setEnergy(FF &ff, const RestraintEnergy &nrg) const
{
    FFComponent::setEnergy(ff, this->total(), nrg);
}

/** Change the restraint component of the energy in the forcefield 'ff'
    by 'delta' */
void RestraintComponent::changeEnergy(FF &ff, const RestraintEnergy &delta) const
{
    FFComponent::changeEnergy(ff, this->total(), delta);
}

const char* RestraintComponent::typeName()
{
    return QMetaType::typeName( qMetaTypeId<RestraintComponent>() );
}
