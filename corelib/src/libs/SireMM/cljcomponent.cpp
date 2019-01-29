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

#include "cljcomponent.h"

#include "SireStream/datastream.h"

using namespace SireMM;
using namespace SireFF;
using namespace SireCAS;
using namespace SireStream;

//////
////// Implementation of CoulombComponent
//////

static const RegisterMetaType<CoulombComponent> r_coul;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const CoulombComponent &coul)
{
    writeHeader(ds, r_coul, 1);
    ds << static_cast<const FFComponent&>(coul);
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, CoulombComponent &coul)
{
    VersionID v = readHeader(ds, r_coul);
    
    if (v == 1)
    {
        ds >> static_cast<FFComponent&>(coul);
    }
    else
        throw version_error(v, "1", r_coul, CODELOC);
        
    return ds;
}

/** Constructor */
CoulombComponent::CoulombComponent(const FFName &ffname)
                 : FFComponent(ffname, QLatin1String("coulomb"))
{}

/** Construct using the passed forcefield name and suffix */
CoulombComponent::CoulombComponent(const FFName &ffname, const QString &suffix)
                 : FFComponent(ffname, QString("coulomb_{%1}").arg(suffix))
{}

/** Construct from a symbol

    \throw SireError::incompatible_error
*/
CoulombComponent::CoulombComponent(const SireCAS::Symbol &symbol)
                 : FFComponent(symbol, QLatin1String("coulomb"))
{}

/** Copy constructor */  
CoulombComponent::CoulombComponent(const CoulombComponent &other)
                 : FFComponent(other)
{}
  
/** Destructor */  
CoulombComponent::~CoulombComponent()
{}

/** Set the coulomb component of the energy in the forcefield 'ff'
    to equal to the passed CoulombEnergy */
void CoulombComponent::setEnergy(FF &ff, const CoulombEnergy &cnrg) const
{
    FFComponent::setEnergy(ff, this->total(), cnrg);
}

/** Change the coulomb component of the energy in the forcefield 'ff'
    by 'delta' */
void CoulombComponent::changeEnergy(FF &ff, const CoulombEnergy &delta) const
{
    FFComponent::changeEnergy(ff, this->total(), delta);
}

const char* CoulombComponent::typeName()
{
    return QMetaType::typeName( qMetaTypeId<CoulombComponent>() );
}

//////
////// Implementation of LJComponent
//////

static const RegisterMetaType<LJComponent> r_lj;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const LJComponent &lj)
{
    writeHeader(ds, r_lj, 1);
    ds << static_cast<const FFComponent&>(lj);
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, LJComponent &lj)
{
    VersionID v = readHeader(ds, r_lj);
    
    if (v == 1)
    {
        ds >> static_cast<FFComponent&>(lj);
    }
    else
        throw version_error(v, "1", r_lj, CODELOC);
        
    return ds;
}

/** Constructor */
LJComponent::LJComponent(const FFName &ffname)
            : FFComponent(ffname, QLatin1String("LJ"))
{}

/** Construct using the name of the forcefield, and the passed suffix */
LJComponent::LJComponent(const FFName &ffname, const QString &suffix)
            : FFComponent(ffname, QString("LJ_{%1}").arg(suffix))
{}

/** Construct from a symbol

    \throw SireError::incompatible_error
*/
LJComponent::LJComponent(const SireCAS::Symbol &symbol)
            : FFComponent(symbol, QLatin1String("LJ"))
{}

/** Copy constructor */  
LJComponent::LJComponent(const LJComponent &other)
            : FFComponent(other)
{}
  
/** Destructor */  
LJComponent::~LJComponent()
{}

/** Set the LJ component of the energy in the forcefield 'ff'
    to equal to the passed LJEnergy */
void LJComponent::setEnergy(FF &ff, const LJEnergy &ljnrg) const
{
    FFComponent::setEnergy(ff, this->total(), ljnrg);
}

/** Change the LJ component of the energy in the forcefield 'ff'
    by 'delta' */
void LJComponent::changeEnergy(FF &ff, const LJEnergy &delta) const
{
    FFComponent::changeEnergy(ff, this->total(), delta);
}

const char* LJComponent::typeName()
{
    return QMetaType::typeName( qMetaTypeId<LJComponent>() );
}

//////
////// Implementation of CLJComponent
//////

static const RegisterMetaType<CLJComponent> r_clj;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const CLJComponent &clj)
{
    writeHeader(ds, r_clj, 1);
    ds << static_cast<const FFComponent&>(clj)
       << clj.coul_component << clj.lj_component;
       
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, CLJComponent &clj)
{
    VersionID v = readHeader(ds, r_clj);
    
    if (v == 1)
    {
        ds >> static_cast<FFComponent&>(clj)
           >> clj.coul_component >> clj.lj_component;
    }
    else
        throw version_error(v, "1", r_clj, CODELOC);
        
    return ds;
}

/** Constructor */
CLJComponent::CLJComponent(const FFName &ffname)
            : FFComponent(ffname, QLatin1String("CLJ")),
               coul_component(ffname), lj_component(ffname)
{}

/** Construct from the passed forcefield name and suffix */
CLJComponent::CLJComponent(const FFName &ffname, const QString &suffix)
             : FFComponent(ffname, QString("CLJ_{%1}").arg(suffix)),
               coul_component(ffname,suffix), lj_component(ffname,suffix)
{}

/** Construct from a symbol

    \throw SireError::incompatible_error
*/
CLJComponent::CLJComponent(const SireCAS::Symbol &symbol)
            : FFComponent(symbol, QLatin1String("CLJ"))
{
    coul_component = CoulombComponent( this->forceFieldName() );
    lj_component = LJComponent( this->forceFieldName() );
}

/** Copy constructor */  
CLJComponent::CLJComponent(const CLJComponent &other)
            : FFComponent(other), coul_component(other.coul_component), 
              lj_component(other.lj_component)
{}
  
/** Destructor */  
CLJComponent::~CLJComponent()
{}

/** Set the CLJ components of the forcefield 'ff' to the passed values */
void CLJComponent::setEnergy(FF &ff, const CLJEnergy &value) const
{
    FFComponent::setEnergy(ff, this->total(), value.total());
    FFComponent::setEnergy(ff, this->coulomb(), value.coulomb());
    FFComponent::setEnergy(ff, this->lj(), value.lj());
}

/** Change the CLJ components of the forcefield 'ff' by 'delta' */
void CLJComponent::changeEnergy(FF &ff, const CLJEnergy &delta) const
{
    FFComponent::changeEnergy(ff, this->total(), delta.total());
    FFComponent::changeEnergy(ff, this->coulomb(), delta.coulomb());
    FFComponent::changeEnergy(ff, this->lj(), delta.lj());
}

Symbols CLJComponent::symbols() const
{
    Symbols symbls;
    symbls.reserve(3);
    
    symbls.insert(coul_component);
    symbls.insert(lj_component);
    symbls.insert(*this);
    
    return symbls;
}

const char* CLJComponent::typeName()
{
    return QMetaType::typeName( qMetaTypeId<CLJComponent>() );
}
