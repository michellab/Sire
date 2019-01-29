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

#include "internalcomponent.h"

#include "SireFF/ff.h"

#include "SireStream/datastream.h"

using namespace SireMM;
using namespace SireFF;
using namespace SireCAS;
using namespace SireStream;

//////
////// Implementation of BondComponent
//////

static const RegisterMetaType<BondComponent> r_bond;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const BondComponent &bond)
{
    writeHeader(ds, r_bond, 1);
    ds << static_cast<const FFComponent&>(bond);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, BondComponent &bond)
{
    VersionID v = readHeader(ds, r_bond);
    
    if (v == 1)
    {
        ds >> static_cast<FFComponent&>(bond);
    }
    else
        throw version_error(v, "1", r_bond, CODELOC);
        
    return ds;
}

/** Constructor */
BondComponent::BondComponent(const FFName &ffname)
              : FFComponent(ffname, QLatin1String("bond"))
{}

/** Construct from a symbol

    \throw SireError::incompatible_error
*/
BondComponent::BondComponent(const SireCAS::Symbol &symbol)
              : FFComponent(symbol, QLatin1String("bond"))
{}

/** Copy constructor */  
BondComponent::BondComponent(const BondComponent &other)
              : FFComponent(other)
{}
  
/** Destructor */  
BondComponent::~BondComponent()
{}

/** Set the component of the energy in the forcefield 'ff'
    to be equal to the passed energy */
void BondComponent::setEnergy(FF &ff, const BondEnergy &nrg) const
{
    FFComponent::setEnergy(ff, this->total(), nrg);
}

/** Change the component of the energy in the forcefield 'ff'
    by 'delta' */
void BondComponent::changeEnergy(FF &ff, const BondEnergy &delta) const
{
    FFComponent::changeEnergy(ff, this->total(), delta);
}

//////
////// Implementation of AngleComponent
//////

static const RegisterMetaType<AngleComponent> r_angle;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const AngleComponent &angle)
{
    writeHeader(ds, r_angle, 1);
    ds << static_cast<const FFComponent&>(angle);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, AngleComponent &angle)
{
    VersionID v = readHeader(ds, r_angle);
    
    if (v == 1)
    {
        ds >> static_cast<FFComponent&>(angle);
    }
    else
        throw version_error(v, "1", r_angle, CODELOC);
        
    return ds;
}

/** Constructor */
AngleComponent::AngleComponent(const FFName &ffname)
               : FFComponent(ffname, QLatin1String("angle"))
{}

/** Construct from a symbol

    \throw SireError::incompatible_error
*/
AngleComponent::AngleComponent(const SireCAS::Symbol &symbol)
               : FFComponent(symbol, QLatin1String("angle"))
{}

/** Copy constructor */  
AngleComponent::AngleComponent(const AngleComponent &other)
              : FFComponent(other)
{}
  
/** Destructor */  
AngleComponent::~AngleComponent()
{}

/** Set the component of the energy in the forcefield 'ff'
    to be equal to the passed energy */
void AngleComponent::setEnergy(FF &ff, const AngleEnergy &nrg) const
{
    FFComponent::setEnergy(ff, this->total(), nrg);
}

/** Change the component of the energy in the forcefield 'ff'
    by 'delta' */
void AngleComponent::changeEnergy(FF &ff, const AngleEnergy &delta) const
{
    FFComponent::changeEnergy(ff, this->total(), delta);
}

//////
////// Implementation of DihedralComponent
//////

static const RegisterMetaType<DihedralComponent> r_dihedral;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const DihedralComponent &dihedral)
{
    writeHeader(ds, r_dihedral, 1);
    ds << static_cast<const FFComponent&>(dihedral);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, DihedralComponent &dihedral)
{
    VersionID v = readHeader(ds, r_dihedral);
    
    if (v == 1)
    {
        ds >> static_cast<FFComponent&>(dihedral);
    }
    else
        throw version_error(v, "1", r_dihedral, CODELOC);
        
    return ds;
}

/** Constructor */
DihedralComponent::DihedralComponent(const FFName &ffname)
                  : FFComponent(ffname, QLatin1String("dihedral"))
{}

/** Construct from a symbol

    \throw SireError::incompatible_error
*/
DihedralComponent::DihedralComponent(const SireCAS::Symbol &symbol)
                  : FFComponent(symbol, QLatin1String("dihedral"))
{}

/** Copy constructor */  
DihedralComponent::DihedralComponent(const DihedralComponent &other)
              : FFComponent(other)
{}
  
/** Destructor */  
DihedralComponent::~DihedralComponent()
{}

/** Set the component of the energy in the forcefield 'ff'
    to be equal to the passed energy */
void DihedralComponent::setEnergy(FF &ff, const DihedralEnergy &nrg) const
{
    FFComponent::setEnergy(ff, this->total(), nrg);
}

/** Change the component of the energy in the forcefield 'ff'
    by 'delta' */
void DihedralComponent::changeEnergy(FF &ff, const DihedralEnergy &delta) const
{
    FFComponent::changeEnergy(ff, this->total(), delta);
}

//////
////// Implementation of ImproperComponent
//////

static const RegisterMetaType<ImproperComponent> r_improper;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const ImproperComponent &improper)
{
    writeHeader(ds, r_improper, 1);
    ds << static_cast<const FFComponent&>(improper);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, ImproperComponent &improper)
{
    VersionID v = readHeader(ds, r_improper);
    
    if (v == 1)
    {
        ds >> static_cast<FFComponent&>(improper);
    }
    else
        throw version_error(v, "1", r_improper, CODELOC);
        
    return ds;
}

/** Constructor */
ImproperComponent::ImproperComponent(const FFName &ffname)
                  : FFComponent(ffname, QLatin1String("improper"))
{}

/** Construct from a symbol

    \throw SireError::incompatible_error
*/
ImproperComponent::ImproperComponent(const SireCAS::Symbol &symbol)
                  : FFComponent(symbol, QLatin1String("improper"))
{}

/** Copy constructor */  
ImproperComponent::ImproperComponent(const ImproperComponent &other)
                  : FFComponent(other)
{}
  
/** Destructor */  
ImproperComponent::~ImproperComponent()
{}

/** Set the component of the energy in the forcefield 'ff'
    to be equal to the passed energy */
void ImproperComponent::setEnergy(FF &ff, const ImproperEnergy &nrg) const
{
    FFComponent::setEnergy(ff, this->total(), nrg);
}

/** Change the component of the energy in the forcefield 'ff'
    by 'delta' */
void ImproperComponent::changeEnergy(FF &ff, const ImproperEnergy &delta) const
{
    FFComponent::changeEnergy(ff, this->total(), delta);
}

//////
////// Implementation of UreyBradleyComponent
//////

static const RegisterMetaType<UreyBradleyComponent> r_ub;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const UreyBradleyComponent &ub)
{
    writeHeader(ds, r_ub, 1);
    ds << static_cast<const FFComponent&>(ub);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, UreyBradleyComponent &ub)
{
    VersionID v = readHeader(ds, r_ub);
    
    if (v == 1)
    {
        ds >> static_cast<FFComponent&>(ub);
    }
    else
        throw version_error(v, "1", r_ub, CODELOC);
        
    return ds;
}

/** Constructor */
UreyBradleyComponent::UreyBradleyComponent(const FFName &ffname)
                     : FFComponent(ffname, QLatin1String("urey_bradley"))
{}

/** Construct from a symbol

    \throw SireError::incompatible_error
*/
UreyBradleyComponent::UreyBradleyComponent(const SireCAS::Symbol &symbol)
                     : FFComponent(symbol, QLatin1String("urey_bradley"))
{}

/** Copy constructor */  
UreyBradleyComponent::UreyBradleyComponent(const UreyBradleyComponent &other)
                     : FFComponent(other)
{}
  
/** Destructor */  
UreyBradleyComponent::~UreyBradleyComponent()
{}

/** Set the component of the energy in the forcefield 'ff'
    to be equal to the passed energy */
void UreyBradleyComponent::setEnergy(FF &ff, const UreyBradleyEnergy &nrg) const
{
    FFComponent::setEnergy(ff, this->total(), nrg);
}

/** Change the component of the energy in the forcefield 'ff'
    by 'delta' */
void UreyBradleyComponent::changeEnergy(FF &ff, const UreyBradleyEnergy &delta) const
{
    FFComponent::changeEnergy(ff, this->total(), delta);
}

//////
////// Implementation of StretchStretchComponent
//////

static const RegisterMetaType<StretchStretchComponent> r_ss;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const StretchStretchComponent &ss)
{
    writeHeader(ds, r_ss, 1);
    ds << static_cast<const FFComponent&>(ss);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, StretchStretchComponent &ss)
{
    VersionID v = readHeader(ds, r_ss);
    
    if (v == 1)
    {
        ds >> static_cast<FFComponent&>(ss);
    }
    else
        throw version_error(v, "1", r_ss, CODELOC);
        
    return ds;
}

/** Constructor */
StretchStretchComponent::StretchStretchComponent(const FFName &ffname)
                        : FFComponent(ffname, QLatin1String("stretch-stretch"))
{}

/** Construct from a symbol

    \throw SireError::incompatible_error
*/
StretchStretchComponent::StretchStretchComponent(const SireCAS::Symbol &symbol)
                        : FFComponent(symbol, QLatin1String("stretch-stretch"))
{}

/** Copy constructor */  
StretchStretchComponent::StretchStretchComponent(const StretchStretchComponent &other)
                        : FFComponent(other)
{}
  
/** Destructor */  
StretchStretchComponent::~StretchStretchComponent()
{}

/** Set the component of the energy in the forcefield 'ff'
    to be equal to the passed energy */
void StretchStretchComponent::setEnergy(FF &ff, const StretchStretchEnergy &nrg) const
{
    FFComponent::setEnergy(ff, this->total(), nrg);
}

/** Change the component of the energy in the forcefield 'ff'
    by 'delta' */
void StretchStretchComponent::changeEnergy(FF &ff, 
                                           const StretchStretchEnergy &delta) const
{
    FFComponent::changeEnergy(ff, this->total(), delta);
}

//////
////// Implementation of StretchBendComponent
//////

static const RegisterMetaType<StretchBendComponent> r_sb;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const StretchBendComponent &sb)
{
    writeHeader(ds, r_sb, 1);
    ds << static_cast<const FFComponent&>(sb);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, StretchBendComponent &sb)
{
    VersionID v = readHeader(ds, r_sb);
    
    if (v == 1)
    {
        ds >> static_cast<FFComponent&>(sb);
    }
    else
        throw version_error(v, "1", r_sb, CODELOC);
        
    return ds;
}

/** Constructor */
StretchBendComponent::StretchBendComponent(const FFName &ffname)
                     : FFComponent(ffname, QLatin1String("stretch-bend"))
{}

/** Construct from a symbol

    \throw SireError::incompatible_error
*/
StretchBendComponent::StretchBendComponent(const SireCAS::Symbol &symbol)
                     : FFComponent(symbol, QLatin1String("stretch-bend"))
{}

/** Copy constructor */  
StretchBendComponent::StretchBendComponent(const StretchBendComponent &other)
                     : FFComponent(other)
{}
  
/** Destructor */  
StretchBendComponent::~StretchBendComponent()
{}

/** Set the component of the energy in the forcefield 'ff'
    to be equal to the passed energy */
void StretchBendComponent::setEnergy(FF &ff, const StretchBendEnergy &nrg) const
{
    FFComponent::setEnergy(ff, this->total(), nrg);
}

/** Change the component of the energy in the forcefield 'ff'
    by 'delta' */
void StretchBendComponent::changeEnergy(FF &ff, const StretchBendEnergy &delta) const
{
    FFComponent::changeEnergy(ff, this->total(), delta);
}

//////
////// Implementation of BendBendComponent
//////

static const RegisterMetaType<BendBendComponent> r_bb;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const BendBendComponent &bb)
{
    writeHeader(ds, r_bb, 1);
    ds << static_cast<const FFComponent&>(bb);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, BendBendComponent &bb)
{
    VersionID v = readHeader(ds, r_bb);
    
    if (v == 1)
    {
        ds >> static_cast<FFComponent&>(bb);
    }
    else
        throw version_error(v, "1", r_bb, CODELOC);
        
    return ds;
}

/** Constructor */
BendBendComponent::BendBendComponent(const FFName &ffname)
                  : FFComponent(ffname, QLatin1String("bend-bend"))
{}

/** Construct from a symbol

    \throw SireError::incompatible_error
*/
BendBendComponent::BendBendComponent(const SireCAS::Symbol &symbol)
                  : FFComponent(symbol, QLatin1String("bend-bend"))
{}

/** Copy constructor */  
BendBendComponent::BendBendComponent(const BendBendComponent &other)
                  : FFComponent(other)
{}
  
/** Destructor */  
BendBendComponent::~BendBendComponent()
{}

/** Set the component of the energy in the forcefield 'ff'
    to be equal to the passed energy */
void BendBendComponent::setEnergy(FF &ff, const BendBendEnergy &nrg) const
{
    FFComponent::setEnergy(ff, this->total(), nrg);
}

/** Change the component of the energy in the forcefield 'ff'
    by 'delta' */
void BendBendComponent::changeEnergy(FF &ff, const BendBendEnergy &delta) const
{
    FFComponent::changeEnergy(ff, this->total(), delta);
}

//////
////// Implementation of StretchBendTorsionComponent
//////

static const RegisterMetaType<StretchBendTorsionComponent> r_sbt;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, 
                                      const StretchBendTorsionComponent &sbt)
{
    writeHeader(ds, r_sbt, 1);
    ds << static_cast<const FFComponent&>(sbt);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, StretchBendTorsionComponent &sbt)
{
    VersionID v = readHeader(ds, r_sbt);
    
    if (v == 1)
    {
        ds >> static_cast<FFComponent&>(sbt);
    }
    else
        throw version_error(v, "1", r_bond, CODELOC);
        
    return ds;
}

/** Constructor */
StretchBendTorsionComponent::StretchBendTorsionComponent(const FFName &ffname)
                  : FFComponent(ffname, QLatin1String("stretch-bend-torsion"))
{}

/** Construct from a symbol

    \throw SireError::incompatible_error
*/
StretchBendTorsionComponent::StretchBendTorsionComponent(const SireCAS::Symbol &symbol)
                  : FFComponent(symbol, QLatin1String("stretch-bend-torsion"))
{}

/** Copy constructor */  
StretchBendTorsionComponent::StretchBendTorsionComponent(
                                        const StretchBendTorsionComponent &other)
              : FFComponent(other)
{}
  
/** Destructor */  
StretchBendTorsionComponent::~StretchBendTorsionComponent()
{}

/** Set the component of the energy in the forcefield 'ff'
    to be equal to the passed energy */
void StretchBendTorsionComponent::setEnergy(FF &ff, 
                                        const StretchBendTorsionEnergy &nrg) const
{
    FFComponent::setEnergy(ff, this->total(), nrg);
}

/** Change the component of the energy in the forcefield 'ff'
    by 'delta' */
void StretchBendTorsionComponent::changeEnergy(FF &ff, 
                                        const StretchBendTorsionEnergy &delta) const
{
    FFComponent::changeEnergy(ff, this->total(), delta);
}

//////
////// Implementation of Intra14CoulombComponent
//////

static const RegisterMetaType<Intra14CoulombComponent> r_14coul;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const Intra14CoulombComponent &comp)
{
    writeHeader(ds, r_14coul, 1);
    ds << static_cast<const FFComponent&>(comp);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, Intra14CoulombComponent &comp)
{
    VersionID v = readHeader(ds, r_14coul);
    
    if (v == 1)
    {
        ds >> static_cast<FFComponent&>(comp);
    }
    else
        throw version_error(v, "1", r_14coul, CODELOC);
        
    return ds;
}

/** Constructor */
Intra14CoulombComponent::Intra14CoulombComponent(const FFName &ffname)
                        : FFComponent(ffname, QLatin1String("1-4[coulomb]"))
{}

/** Construct from a symbol

    \throw SireError::incompatible_error
*/
Intra14CoulombComponent::Intra14CoulombComponent(const SireCAS::Symbol &symbol)
                        : FFComponent(symbol, QLatin1String("1-4[coulomb]"))
{}

/** Copy constructor */  
Intra14CoulombComponent::Intra14CoulombComponent(const Intra14CoulombComponent &other)
                        : FFComponent(other)
{}
  
/** Destructor */  
Intra14CoulombComponent::~Intra14CoulombComponent()
{}

/** Set the component of the energy in the forcefield 'ff'
    to be equal to the passed energy */
void Intra14CoulombComponent::setEnergy(FF &ff, const Intra14CoulombEnergy &nrg) const
{
    FFComponent::setEnergy(ff, this->total(), nrg);
}

/** Change the component of the energy in the forcefield 'ff'
    by 'delta' */
void Intra14CoulombComponent::changeEnergy(FF &ff, const Intra14CoulombEnergy &delta) const
{
    FFComponent::changeEnergy(ff, this->total(), delta);
}

//////
////// Implementation of Intra14LJComponent
//////

static const RegisterMetaType<Intra14LJComponent> r_14lj;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const Intra14LJComponent &comp)
{
    writeHeader(ds, r_14lj, 1);
    ds << static_cast<const FFComponent&>(comp);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, Intra14LJComponent &comp)
{
    VersionID v = readHeader(ds, r_14lj);
    
    if (v == 1)
    {
        ds >> static_cast<FFComponent&>(comp);
    }
    else
        throw version_error(v, "1", r_14lj, CODELOC);
        
    return ds;
}

/** Constructor */
Intra14LJComponent::Intra14LJComponent(const FFName &ffname)
                   : FFComponent(ffname, QLatin1String("1-4[LJ]"))
{}

/** Construct from a symbol

    \throw SireError::incompatible_error
*/
Intra14LJComponent::Intra14LJComponent(const SireCAS::Symbol &symbol)
                   : FFComponent(symbol, QLatin1String("1-4[LJ]"))
{}

/** Copy constructor */  
Intra14LJComponent::Intra14LJComponent(const Intra14LJComponent &other)
                   : FFComponent(other)
{}
  
/** Destructor */  
Intra14LJComponent::~Intra14LJComponent()
{}

/** Set the component of the energy in the forcefield 'ff'
    to be equal to the passed energy */
void Intra14LJComponent::setEnergy(FF &ff, const Intra14LJEnergy &nrg) const
{
    FFComponent::setEnergy(ff, this->total(), nrg);
}

/** Change the component of the energy in the forcefield 'ff'
    by 'delta' */
void Intra14LJComponent::changeEnergy(FF &ff, const Intra14LJEnergy &delta) const
{
    FFComponent::changeEnergy(ff, this->total(), delta);
}

//////
////// Implementation of Intra14Component
//////

static const RegisterMetaType<Intra14Component> r_14;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const Intra14Component &comp)
{
    writeHeader(ds, r_14, 1);
    
    ds << static_cast<const FFComponent&>(comp)
       << comp.coul_component << comp.lj_component;
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, Intra14Component &comp)
{
    VersionID v = readHeader(ds, r_14);
    
    if (v == 1)
    {
        ds >> static_cast<FFComponent&>(comp)
           >> comp.coul_component >> comp.lj_component;
    }
    else
        throw version_error(v, "1", r_14, CODELOC);
        
    return ds;
}

/** Constructor */
Intra14Component::Intra14Component(const FFName &ffname)
                  : FFComponent(ffname, QLatin1String("1-4")),
                    coul_component(ffname), lj_component(ffname)
{}

/** Construct from a symbol

    \throw SireError::incompatible_error
*/
Intra14Component::Intra14Component(const SireCAS::Symbol &symbol)
                 : FFComponent(symbol, QLatin1String("1-4"))
{
    coul_component = Intra14CoulombComponent(this->forceFieldName());
    lj_component = Intra14LJComponent(this->forceFieldName());
}

/** Copy constructor */  
Intra14Component::Intra14Component(const Intra14Component &other)
                 : FFComponent(other),
                   coul_component(other.coul_component),
                   lj_component(other.lj_component)
{}

/** Destructor */  
Intra14Component::~Intra14Component()
{}

/** Return all of the components in this set */
Symbols Intra14Component::symbols() const
{
    Symbols symbls;
    
    symbls.reserve(3);
    
    symbls.insert(*this);
    symbls.insert(coul_component);
    symbls.insert(lj_component);
    
    return symbls;
}

/** Set the internal components of the forcefield 'ff' to the passed values */
void Intra14Component::setEnergy(FF &ff, const Intra14Energy &value) const
{
    FFComponent::setEnergy(ff, this->total(), value.total());
    FFComponent::setEnergy(ff, this->coulomb(), value.coulomb());
    FFComponent::setEnergy(ff, this->lj(), value.lj());
}

/** Change the internal components of the forcefield 'ff' by 'delta' */
void Intra14Component::changeEnergy(FF &ff, const Intra14Energy &delta) const
{
    FFComponent::changeEnergy(ff, this->total(), delta.total());
    FFComponent::changeEnergy(ff, this->coulomb(), delta.coulomb());
    FFComponent::changeEnergy(ff, this->lj(), delta.lj());
}

//////
////// Implementation of InternalComponent
//////

static const RegisterMetaType<InternalComponent> r_internal;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const InternalComponent &internal)
{
    writeHeader(ds, r_internal, 2);
    
    ds << static_cast<const FFComponent&>(internal)
       << internal.bond_component
       << internal.angle_component
       << internal.dihedral_component
       << internal.improper_component
       << internal.ub_component
       << internal.ss_component
       << internal.sb_component
       << internal.bb_component
       << internal.sbt_component
       << internal.nb_component;
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, InternalComponent &internal)
{
    VersionID v = readHeader(ds, r_internal);
    
    if (v == 2)
    {
        ds >> static_cast<FFComponent&>(internal)
           >> internal.bond_component
           >> internal.angle_component
           >> internal.dihedral_component
           >> internal.improper_component
           >> internal.ub_component
           >> internal.ss_component
           >> internal.sb_component
           >> internal.bb_component
           >> internal.sbt_component
           >> internal.nb_component;
    }
    else if (v == 1)
    {
        ds >> static_cast<FFComponent&>(internal)
           >> internal.bond_component
           >> internal.angle_component
           >> internal.dihedral_component
           >> internal.improper_component
           >> internal.ub_component
           >> internal.ss_component
           >> internal.sb_component
           >> internal.bb_component
           >> internal.sbt_component;
        
        internal.nb_component = Intra14Component( internal.forceFieldName() );
    }
    else
        throw version_error(v, "1,2", r_internal, CODELOC);
        
    return ds;
}

/** Constructor */
InternalComponent::InternalComponent(const FFName &ffname)
                  : FFComponent(ffname, QLatin1String("internal")),
                    bond_component(ffname), angle_component(ffname),
                    dihedral_component(ffname), 
                    improper_component(ffname), ub_component(ffname),
                    ss_component(ffname), sb_component(ffname),
                    bb_component(ffname), sbt_component(ffname),
                    nb_component(ffname)
{}

/** Construct from a symbol

    \throw SireError::incompatible_error
*/
InternalComponent::InternalComponent(const SireCAS::Symbol &symbol)
                  : FFComponent(symbol, QLatin1String("internal"))
{
    bond_component = BondComponent( this->forceFieldName() );
    angle_component = AngleComponent( this->forceFieldName() );
    dihedral_component = DihedralComponent( this->forceFieldName() );
    
    improper_component = ImproperComponent( this->forceFieldName() );
    ub_component = UreyBradleyComponent( this->forceFieldName() );
    
    ss_component = StretchStretchComponent( this->forceFieldName() );
    sb_component = StretchBendComponent( this->forceFieldName() );
    bb_component = BendBendComponent( this->forceFieldName() );
    sbt_component = StretchBendTorsionComponent( this->forceFieldName() );
    
    nb_component = Intra14Component( this->forceFieldName() );
}

/** Copy constructor */  
InternalComponent::InternalComponent(const InternalComponent &other)
                  : FFComponent(other), 
                    bond_component(other.bond_component),
                    angle_component(other.angle_component),
                    dihedral_component(other.dihedral_component),
                    improper_component(other.improper_component),
                    ub_component(other.ub_component),
                    ss_component(other.ss_component),
                    sb_component(other.sb_component),
                    bb_component(other.bb_component),
                    sbt_component(other.sbt_component),
                    nb_component(other.nb_component)
{}
  
/** Destructor */  
InternalComponent::~InternalComponent()
{}

/** Return all of the components in this set */
Symbols InternalComponent::symbols() const
{
    Symbols symbls;
    
    symbls.reserve(13);
    
    symbls.insert(*this);
    
    symbls.insert(bond_component);
    symbls.insert(angle_component);
    symbls.insert(dihedral_component);
    
    symbls.insert(improper_component);
    symbls.insert(ub_component);
    
    symbls.insert(ss_component);
    symbls.insert(sb_component);
    symbls.insert(bb_component);
    symbls.insert(sbt_component);
    
    symbls.insert(nb_component);
    symbls.insert(nb_component.coulomb());
    symbls.insert(nb_component.lj());
    
    return symbls;
}

/** Set the internal components of the forcefield 'ff' to the passed values */
void InternalComponent::setEnergy(FF &ff, const InternalEnergy &value) const
{
    FFComponent::setEnergy(ff, this->total(), value.total());
    FFComponent::setEnergy(ff, this->bond(), value.bond());
    FFComponent::setEnergy(ff, this->angle(), value.angle());
    FFComponent::setEnergy(ff, this->dihedral(), value.dihedral());
    
    FFComponent::setEnergy(ff, this->improper(), value.improper());
    FFComponent::setEnergy(ff, this->ureyBradley(), value.ureyBradley());
    
    FFComponent::setEnergy(ff, this->stretchStretch(), value.stretchStretch());
    FFComponent::setEnergy(ff, this->stretchBend(), value.stretchBend());
    FFComponent::setEnergy(ff, this->bendBend(), value.bendBend());
    FFComponent::setEnergy(ff, this->stretchBendTorsion(), value.stretchBendTorsion());
    
    FFComponent::setEnergy(ff, this->intra14(), value.intra14());
    FFComponent::setEnergy(ff, this->intra14Coulomb(), value.intra14Coulomb());
    FFComponent::setEnergy(ff, this->intra14LJ(), value.intra14LJ());
}

/** Change the internal components of the forcefield 'ff' by 'delta' */
void InternalComponent::changeEnergy(FF &ff, const InternalEnergy &delta) const
{
    FFComponent::changeEnergy(ff, this->total(), delta.total());
    FFComponent::changeEnergy(ff, this->bond(), delta.bond());
    FFComponent::changeEnergy(ff, this->angle(), delta.angle());
    FFComponent::changeEnergy(ff, this->dihedral(), delta.dihedral());
    
    FFComponent::changeEnergy(ff, this->improper(), delta.improper());
    FFComponent::changeEnergy(ff, this->ureyBradley(), delta.ureyBradley());
    
    FFComponent::changeEnergy(ff, this->stretchStretch(), delta.stretchStretch());
    FFComponent::changeEnergy(ff, this->stretchBend(), delta.stretchBend());
    FFComponent::changeEnergy(ff, this->bendBend(), delta.bendBend());
    FFComponent::changeEnergy(ff, this->stretchBendTorsion(), 
                              delta.stretchBendTorsion());

    FFComponent::changeEnergy(ff, this->intra14(), delta.intra14());
    FFComponent::changeEnergy(ff, this->intra14Coulomb(), delta.intra14Coulomb());
    FFComponent::changeEnergy(ff, this->intra14LJ(), delta.intra14LJ());
}

/////////
///////// Implementation of InternalEnergy
/////////

/** Constructor */
InternalEnergy::InternalEnergy(double bondnrg, double anglenrg,
                               double dihedralnrg, 
                               double impropernrg, double ubnrg,
                               double ssnrg, double sbnrg,
                               double bbnrg, double sbtnrg,
                               Intra14Energy inrg)
               : i14nrg(inrg),
                 ibndnrg(bondnrg), iangnrg(anglenrg),
                 idihnrg(dihedralnrg), 
                 iimpnrg(impropernrg), iubnrg(ubnrg),
                 issnrg(ssnrg), isbnrg(sbnrg),
                 ibbnrg(bbnrg), isbtnrg(sbtnrg)
{}
  
/** Copy constructor */  
InternalEnergy::InternalEnergy(const InternalEnergy &other)
               : i14nrg(other.i14nrg),
                 ibndnrg(other.ibndnrg), iangnrg(other.iangnrg),
                 idihnrg(other.idihnrg), 
                 iimpnrg(other.iimpnrg), iubnrg(other.iubnrg),
                 issnrg(other.issnrg), isbnrg(other.isbnrg),
                 ibbnrg(other.ibbnrg), isbtnrg(other.isbtnrg)
{}
  
/** Destructor */  
InternalEnergy::~InternalEnergy()
{}

/////////
///////// Implementation of Intra14Energy
/////////

/** Constructor */
Intra14Energy::Intra14Energy(double coul_nrg, double lj_nrg)
               : cnrg(coul_nrg), ljnrg(lj_nrg)
{}
  
/** Copy constructor */  
Intra14Energy::Intra14Energy(const Intra14Energy &other)
              : cnrg(other.cnrg), ljnrg(other.ljnrg)
{}
  
/** Destructor */  
Intra14Energy::~Intra14Energy()
{}

/// typename functions

const char* BondComponent::typeName()
{
    return QMetaType::typeName( qMetaTypeId<BondComponent>() );
}

const char* AngleComponent::typeName()
{
    return QMetaType::typeName( qMetaTypeId<AngleComponent>() );
}

const char* DihedralComponent::typeName()
{
    return QMetaType::typeName( qMetaTypeId<DihedralComponent>() );
}

const char* ImproperComponent::typeName()
{
    return QMetaType::typeName( qMetaTypeId<ImproperComponent>() );
}

const char* UreyBradleyComponent::typeName()
{
    return QMetaType::typeName( qMetaTypeId<UreyBradleyComponent>() );
}

const char* StretchStretchComponent::typeName()
{
    return QMetaType::typeName( qMetaTypeId<StretchStretchComponent>() );
}

const char* StretchBendComponent::typeName()
{
    return QMetaType::typeName( qMetaTypeId<StretchBendComponent>() );
}

const char* BendBendComponent::typeName()
{
    return QMetaType::typeName( qMetaTypeId<BendBendComponent>() );
}

const char* StretchBendTorsionComponent::typeName()
{
    return QMetaType::typeName( qMetaTypeId<StretchBendTorsionComponent>() );
}

const char* Intra14CoulombComponent::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Intra14CoulombComponent>() );
}

const char* Intra14LJComponent::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Intra14LJComponent>() );
}

const char* Intra14Component::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Intra14Component>() );
}

const char* InternalComponent::typeName()
{
    return QMetaType::typeName( qMetaTypeId<InternalComponent>() );
}
