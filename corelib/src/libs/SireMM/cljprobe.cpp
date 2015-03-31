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

#include "cljprobe.h"

#include "SireUnits/units.h"

#include "SireStream/datastream.h"

using namespace SireMM;
using namespace SireFF;
using namespace SireUnits;
using namespace SireUnits::Dimension;
using namespace SireBase;
using namespace SireStream;

/////////
///////// Implementation of CoulombProbe
/////////

static const RegisterMetaType<CoulombProbe> r_cprobe;

QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, 
                                      const CoulombProbe &cprobe)
{
    writeHeader(ds, r_cprobe, 1);
    
    ds << cprobe.chg.to( mod_electron )
       << static_cast<const Probe&>(cprobe);
       
    return ds;
}

QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, CoulombProbe &cprobe)
{
    VersionID v = readHeader(ds, r_cprobe);
    
    if (v == 1)
    {
        double chg;
        
        ds >> chg;
        
        cprobe = CoulombProbe( chg*mod_electron );
        
        ds >> static_cast<Probe&>(cprobe);
    }
    else
        throw version_error(v, "1", r_cprobe, CODELOC);
        
    return ds;
}

static double getCharge(const Charge &charge)
{
    const double sqrt_4pieps0 = std::sqrt(SireUnits::one_over_four_pi_eps0);
    
    return charge.value() * sqrt_4pieps0;
}

/** Construct a default probe (+1 unit charge) */
CoulombProbe::CoulombProbe() 
             : ConcreteProperty<CoulombProbe,Probe>(), chg(1 * mod_electron)
{
    reduced_chg = ::getCharge(chg);
}

/** Construct a probe with charge 'charge' */
CoulombProbe::CoulombProbe(Charge charge)
             : ConcreteProperty<CoulombProbe,Probe>(), chg(charge)
{
    reduced_chg = ::getCharge(chg);
}

/** Construct a probe with charge taken from 'cljprobe' */
CoulombProbe::CoulombProbe(const CLJProbe &cljprobe)
             : ConcreteProperty<CoulombProbe,Probe>(), chg(cljprobe.charge())
{
    reduced_chg = ::getCharge(chg);
}

/** Construct from the passed probe */
CoulombProbe::CoulombProbe(const Probe &probe)
             : ConcreteProperty<CoulombProbe,Probe>(), chg(0)
{
    if (probe.isA<CoulombProbe>())
        chg = probe.asA<CoulombProbe>().chg;

    else if (probe.isA<CLJProbe>())
        chg = probe.asA<CLJProbe>().charge();
        
    reduced_chg = ::getCharge(chg);
}

/** Copy constructor */
CoulombProbe::CoulombProbe(const CoulombProbe &other)
             : ConcreteProperty<CoulombProbe,Probe>(other),
               chg(other.chg), reduced_chg(other.reduced_chg)
{}

/** Destructor */
CoulombProbe::~CoulombProbe()
{}

/** Copy assignment operator */
CoulombProbe& CoulombProbe::operator=(const CoulombProbe &other)
{
    if (this != &other)
    {
        Probe::operator=(other);
        chg = other.chg;
        reduced_chg = other.reduced_chg;
    }
    
    return *this;
}

/** Comparison operator */
bool CoulombProbe::operator==(const CoulombProbe &other) const
{
    return reduced_chg == other.reduced_chg and Probe::operator==(other);
}

/** Comparison operator */
bool CoulombProbe::operator!=(const CoulombProbe &other) const
{
    return not CoulombProbe::operator==(other);
}

const char* CoulombProbe::typeName()
{
    return QMetaType::typeName( qMetaTypeId<CoulombProbe>() );
}

/////////
///////// Implementation of LJProbe
/////////

static const RegisterMetaType<LJProbe> r_ljprobe;

QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, const LJProbe &ljprobe)
{
    writeHeader(ds, r_ljprobe, 1);
    
    ds << ljprobe.ljparam << static_cast<const Probe&>(ljprobe);
    
    return ds;
}

QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, LJProbe &ljprobe)
{
    VersionID v = readHeader(ds, r_ljprobe);
    
    if (v == 1)
    {
        ds >> ljprobe.ljparam >> static_cast<Probe&>(ljprobe);
    }
    else
        throw version_error(v, "1", r_ljprobe, CODELOC);
        
    return ds;
}

static const LJParameter ua_methane( 3.73*angstrom, 0.294*kcal_per_mol );

/** Constructor - this makes a probe that is the equivalent
    of an OPLS united atom methane molecule */
LJProbe::LJProbe() : ConcreteProperty<LJProbe,Probe>(), ljparam(ua_methane)
{}

/** Construct a probe with parameters in 'ljparam' */
LJProbe::LJProbe(const LJParameter &lj)
        : ConcreteProperty<LJProbe,Probe>(), ljparam(lj)
{}

/** Construct to take the LJ probe from the passed CLJProbe */
LJProbe::LJProbe(const CLJProbe &cljprobe)
        : ConcreteProperty<LJProbe,Probe>(), ljparam(cljprobe.lj())
{}

/** Construct from the passed probe */
LJProbe::LJProbe(const Probe &probe)
        : ConcreteProperty<LJProbe,Probe>()
{
    if (probe.isA<LJProbe>())
        ljparam = probe.asA<LJProbe>().ljparam;

    else if (probe.isA<CLJProbe>())
        ljparam = probe.asA<CLJProbe>().lj();
}

/** Copy constructor */
LJProbe::LJProbe(const LJProbe &other)
        : ConcreteProperty<LJProbe,Probe>(other), ljparam(other.ljparam)
{}

/** Destructor */
LJProbe::~LJProbe()
{}

/** Copy assignment operator */
LJProbe& LJProbe::operator=(const LJProbe &other)
{
    if (this != &other)
    {
        ljparam = other.ljparam;
        Probe::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool LJProbe::operator==(const LJProbe &other) const
{
    return ljparam == other.ljparam and Probe::operator==(other);
}

/** Comparison operator */
bool LJProbe::operator!=(const LJProbe &other) const
{
    return not LJProbe::operator==(other);
}

const char* LJProbe::typeName()
{
    return QMetaType::typeName( qMetaTypeId<LJProbe>() );
}

/////////
///////// Implementation of CLJProbe
/////////

static const RegisterMetaType<CLJProbe> r_cljprobe;

QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, const CLJProbe &cljprobe)
{
    writeHeader(ds, r_cljprobe, 1);
    
    ds << cljprobe.chg.to(mod_electron) << cljprobe.ljparam
       << static_cast<const Probe&>(cljprobe);
       
    return ds;
}

QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, CLJProbe &cljprobe)
{
    VersionID v = readHeader(ds, r_cljprobe);
    
    if (v == 1)
    {
        double chg;
        ds >> chg;
        
        cljprobe = CLJProbe(chg*mod_electron);
        
        ds >> cljprobe.ljparam >> static_cast<Probe&>(cljprobe);
    }
    else
        throw version_error(v, "1", r_cljprobe, CODELOC);
        
    return ds;
}

/** Construct a default probe the represents a unit positive charge,
    and an OPLS united atom methane */
CLJProbe::CLJProbe() : ConcreteProperty<CLJProbe,Probe>(),
                       ljparam(ua_methane), chg(1 * mod_electron)
{
    reduced_chg = ::getCharge(chg);
}

/** Construct a probe that is just the passed charge (with dummy LJ parameters) */
CLJProbe::CLJProbe(SireUnits::Dimension::Charge charge)
         : ConcreteProperty<CLJProbe,Probe>(),
           ljparam( LJParameter::dummy() ), chg(charge)
{
    reduced_chg = ::getCharge(chg);
}

/** Construct a probe that is just the passed LJ parameter (zero charge) */
CLJProbe::CLJProbe(const LJParameter &lj)
         : ConcreteProperty<CLJProbe,Probe>(),
           ljparam(lj), chg(0), reduced_chg(0)
{}

/** Construct a probe with the passed charge and LJ parameter */
CLJProbe::CLJProbe(Charge charge, const LJParameter &lj)
         : ConcreteProperty<CLJProbe,Probe>(),
           ljparam(lj), chg(charge)
{
    reduced_chg = ::getCharge(chg);
}

/** Construct to hold just the charge from 'probe' (dummy LJ) */
CLJProbe::CLJProbe(const CoulombProbe &probe)
         : ConcreteProperty<CLJProbe,Probe>(),
           ljparam( LJParameter::dummy() ), chg(probe.charge())
{
    reduced_chg = ::getCharge(chg);
}

/** Construct to hold just the LJ parameters from 'probe' (zero charge) */
CLJProbe::CLJProbe(const LJProbe &probe)
         : ConcreteProperty<CLJProbe,Probe>(),
           ljparam(probe.lj()), chg(0), reduced_chg(0)
{}

/** Copy constructor */
CLJProbe::CLJProbe(const CLJProbe &other)
         : ConcreteProperty<CLJProbe,Probe>(other),
           ljparam(other.ljparam), chg(other.chg),
           reduced_chg(other.reduced_chg)
{}

/** Destructor */
CLJProbe::~CLJProbe()
{}

/** Construct from the passed probe */
CLJProbe::CLJProbe(const Probe &other)
         : ConcreteProperty<CLJProbe,Probe>(), chg(0), reduced_chg(0)
{
    if (other.isA<CLJProbe>())
        this->operator=(other.asA<CLJProbe>());
        
    else if (other.isA<CoulombProbe>())
        this->operator=(other.asA<CoulombProbe>());
        
    else if (other.isA<LJProbe>())
        this->operator=(other.asA<LJProbe>());
}

/** Copy assignment operator */
CLJProbe& CLJProbe::operator=(const CLJProbe &other)
{
    if (this != &other)
    {
        Probe::operator=(other);
        ljparam = other.ljparam;
        chg = other.chg;
        reduced_chg = other.reduced_chg;
    }
    
    return *this;
}

/** Comparison operator */
bool CLJProbe::operator==(const CLJProbe &other) const
{
    return reduced_chg == other.reduced_chg and 
           ljparam == other.ljparam and Probe::operator==(other);
}

/** Comparison operator */
bool CLJProbe::operator!=(const CLJProbe &other) const
{
    return not CLJProbe::operator==(other);
}

const char* CLJProbe::typeName()
{
    return QMetaType::typeName( qMetaTypeId<CLJProbe>() );
}
