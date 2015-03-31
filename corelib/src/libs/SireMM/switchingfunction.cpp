/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2007  Christopher Woods
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

#include <numeric>
#include <cmath>

#include <QMutex>

#include "switchingfunction.h"

#include "SireMaths/maths.h"

#include "SireUnits/units.h"

#include "SireFF/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMM;
using namespace SireStream;
using namespace SireBase;
using namespace SireUnits;

/////////////
///////////// Implementation of SwitchingFunction
/////////////

static const RegisterMetaType<SwitchingFunction> r_switchbase(MAGIC_ONLY,
                                                           "SireMM::SwitchingFunction");

/** Serialise to a binary data stream */
QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, const SwitchingFunction &switchfunc)
{
    writeHeader(ds, r_switchbase, 1)
        << static_cast<const Property&>(switchfunc)
        << switchfunc.cut_elec << switchfunc.feather_elec
        << switchfunc.cut_vdw << switchfunc.feather_vdw;

    return ds;
}

/** Deserialise from a binary data stream */
QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, SwitchingFunction &switchfunc)
{
    VersionID v = readHeader(ds, r_switchbase);

    if (v == 1)
    {
        ds >> static_cast<Property&>(switchfunc)
           >> switchfunc.cut_elec >> switchfunc.feather_elec
           >> switchfunc.cut_vdw >> switchfunc.feather_vdw;
           
        switchfunc.cutdist = qMax(switchfunc.cut_elec, switchfunc.cut_vdw);
        switchfunc.featherdist = qMax(switchfunc.feather_elec,
                                      switchfunc.feather_vdw);
    }
    else
        throw version_error(v, "1", r_switchbase, CODELOC);

    return ds;
}

/** Null constructor - this places the cutoff distance at
    the maximum value of a double (e.g. on cutoff) */
SwitchingFunction::SwitchingFunction()
                  : Property(),
                    cutdist( std::numeric_limits<double>::max() ),
                    featherdist( std::numeric_limits<double>::max() ),
                    cut_elec( std::numeric_limits<double>::max() ),
                    feather_elec( std::numeric_limits<double>::max() ),
                    cut_vdw( std::numeric_limits<double>::max() ),
                    feather_vdw( std::numeric_limits<double>::max() )
{}

/** Construct, placing the ultimate cutoff distance at 'cutdistance' */
SwitchingFunction::SwitchingFunction(Length cutdistance)
                  : Property(),
                    cutdist( std::abs(cutdistance) ),
                    featherdist(cutdist),
                    cut_elec(cutdist),
                    feather_elec(cutdist),
                    cut_vdw(cutdist),
                    feather_vdw(cutdist)
{}

/** Construct placing the cutoff distance at 'cutdistance' and 
    start feathering at 'featherdistance' */
SwitchingFunction::SwitchingFunction(Length cutdistance, Length featherdistance)
                  : Property(),
                    cutdist( std::abs(cutdistance) ),
                    featherdist( qMin( cutdist, std::abs(featherdistance) ) ),
                    cut_elec(cutdist),
                    feather_elec(featherdist),
                    cut_vdw(cutdist),
                    feather_vdw(featherdist)
{}

/** Construct placing the electrostatic cutoff at 'eleccut' with 
    the feather at 'elecfeather', while the vdw cutoff is at
    'vdwcut' with the feather at 'vdwfeather' */
SwitchingFunction::SwitchingFunction(Length eleccut, Length elecfeather,
                                     Length vdwcut, Length vdwfeather)
                  : Property()
{
    cut_elec = std::abs(eleccut);
    cut_vdw = std::abs(vdwcut);
    
    cutdist = qMax(cut_elec, cut_vdw);
    
    feather_elec = qMin( cut_elec, std::abs(elecfeather) );
    feather_vdw = qMin( cut_vdw, std::abs(vdwfeather) );
    
    featherdist = qMax(feather_elec, feather_vdw);
}

/** Copy constructor */
SwitchingFunction::SwitchingFunction(const SwitchingFunction &other)
                  : Property(other), 
                    cutdist(other.cutdist),
                    featherdist(other.featherdist),
                    cut_elec(other.cut_elec),
                    feather_elec(other.feather_elec),
                    cut_vdw(other.cut_vdw),
                    feather_vdw(other.feather_vdw)
{}

/** Destructor */
SwitchingFunction::~SwitchingFunction()
{}

/** Copy assignment operator */
SwitchingFunction& SwitchingFunction::operator=(const SwitchingFunction &other)
{
    if (this != &other)
    {
        Property::operator=(other);
    
        cutdist = other.cutdist;
        featherdist = other.featherdist;
        
        cut_elec = other.cut_elec;
        feather_elec = other.feather_elec;
        
        cut_vdw = other.cut_vdw;
        feather_vdw = other.feather_vdw;
    }
    
    return *this;
}

/** Comparison operator */
bool SwitchingFunction::operator==(const SwitchingFunction &other) const
{
    return this == &other or
           (cut_elec == other.cut_elec and cut_vdw == other.cut_vdw and
            feather_elec == other.feather_elec and feather_vdw == other.feather_vdw);
}

/** Comparison operator */
bool SwitchingFunction::operator!=(const SwitchingFunction &other) const
{
    return this != &other and
           (cut_elec != other.cut_elec or cut_vdw != other.cut_vdw or
            feather_elec != other.feather_elec or feather_vdw != other.feather_vdw);
}

/////////////
///////////// Implementation of NoCutoff
/////////////

static const RegisterMetaType<NoCutoff> r_nocutoff;

/** Serialise to a binary data stream */
QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, const NoCutoff &nocutoff)
{
    writeHeader(ds, r_nocutoff, 1)
            << static_cast<const SwitchingFunction&>(nocutoff);

    return ds;
}

/** Deserialise from a binary data stream */
QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, NoCutoff &nocutoff)
{
    VersionID v = readHeader(ds, r_nocutoff);

    if (v == 1)
    {
        ds >> static_cast<SwitchingFunction&>(nocutoff);
    }
    else
        throw version_error(v, "1", r_nocutoff, CODELOC);

    return ds;
}

/** Constructor */
NoCutoff::NoCutoff() 
         : ConcreteProperty<NoCutoff,SwitchingFunction>(
                                Length(std::numeric_limits<double>::max()) )
{}

/** Copy constructor */
NoCutoff::NoCutoff(const NoCutoff &other)
         : ConcreteProperty<NoCutoff,SwitchingFunction>(other)
{}

/** Destructor */
NoCutoff::~NoCutoff()
{}

/** Copy assignment operator */
NoCutoff& NoCutoff::operator=(const NoCutoff &other)
{
    SwitchingFunction::operator=(other);
    return *this;
}

/** Comparison operator */
bool NoCutoff::operator==(const NoCutoff&) const
{
    return true;
}

/** Comparison operator */
bool NoCutoff::operator!=(const NoCutoff&) const
{
    return false;
}

/** Return a string representation of this switching function */
QString NoCutoff::toString() const
{
    return QObject::tr("no cutoff");
}

/** Return the scale factor for the electrostatic energies - this
    will always be 1.0, as there are no cutoffs */
double NoCutoff::electrostaticScaleFactor(Length) const
{
    return 1;
}

/** Return the scale factor for the vdw energies - this
    will always be 1.0, as there are no cutoffs */
double NoCutoff::vdwScaleFactor(Length) const
{
    return 1;
}

/** Return the derivative of the electrostatic scale factor - this 
    will always be 0 as there is no cutoff! */
double NoCutoff::dElectrostaticScaleFactor(Length) const
{
    return 0;
}

/** Return the derivative of the VDW scale factor - this will
    always be 0 as there is no cutoff! */
double NoCutoff::dVDWScaleFactor(Length) const
{
    return 0;
}

static SharedPolyPointer<NoCutoff> shared_null;

const NoCutoff& SwitchingFunction::null()
{
    if (shared_null.constData() == 0)
    {
        QMutexLocker lkr( SireBase::globalLock() );
        
        if (shared_null.constData() == 0)
            shared_null = new NoCutoff();
    }
    
    return *(shared_null.constData());
}

const char* NoCutoff::typeName()
{
    return QMetaType::typeName( qMetaTypeId<NoCutoff>() );
}

/////////////
///////////// Implementation of HarmonicSwitchingFunction
/////////////

static const RegisterMetaType<HarmonicSwitchingFunction> r_harm;

/** Internal function used to set the parameters of the harmonic cutoff */
void HarmonicSwitchingFunction::set(double cutelec, double featherelec,
                                    double cutvdw, double feathervdw)
{
    cut_elec = std::abs(cutelec);
    cut_elec2 = SireMaths::pow_2(cutelec);
    feather_elec = qMin( cut_elec, std::abs(featherelec) );

    if (cut_elec != feather_elec)
        norm_elec = 1.0 / (cut_elec2 - SireMaths::pow_2(feather_elec));
    else
        norm_elec = 0;

    cut_vdw = std::abs(cutvdw);
    cut_vdw2 = SireMaths::pow_2(cutvdw);
    feather_vdw = qMin( cut_vdw, std::abs(feathervdw) );

    if (cut_vdw != feather_vdw)
        norm_vdw = 1.0 / (cut_vdw2 - SireMaths::pow_2(feather_vdw));
    else
        norm_vdw = 0;
        
    cutdist = qMax( cut_elec, cut_vdw );
    featherdist = qMax( feather_elec, feather_vdw );
}

/** Serialise to a binary data stream */
QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds,
                                      const HarmonicSwitchingFunction &harm)
{
    writeHeader(ds, r_harm, 1)
          << static_cast<const SwitchingFunction&>(harm);

    return ds;
}

/** Deserialise from a binary data stream */
QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds,
                                      HarmonicSwitchingFunction &harm)
{
    VersionID v = readHeader(ds, r_harm);

    if (v == 1)
    {
        ds >> static_cast<SwitchingFunction&>(harm);

        harm.set( harm.cut_elec, harm.feather_elec,
                  harm.cut_vdw, harm.feather_vdw );
    }
    else
        throw version_error(v, "1", r_harm, CODELOC);

    return ds;
}

/** Construct a null harmonic switching function (no cutoff) */
HarmonicSwitchingFunction::HarmonicSwitchingFunction()
    : ConcreteProperty<HarmonicSwitchingFunction,SwitchingFunction>()
{
    cut_elec2 = cut_elec;
    cut_vdw2 = cut_vdw;
    
    norm_elec = 0;
    norm_vdw = 0;
}

/** Construct an harmonic switching function which represents a hard
    cutoff of both the electrostatic and vdw interactions at a distance
    of 'cutoffdist' */
HarmonicSwitchingFunction::HarmonicSwitchingFunction(Length cutoffdist)
     : ConcreteProperty<HarmonicSwitchingFunction,SwitchingFunction>()
{
    this->set(cutoffdist, cutoffdist, cutoffdist, cutoffdist);
}

/** Construct an harmonic switching function which represents the smoothed
    cutoff, with the cutoff at 'cutoffdist', but smoothed down using an
    harmonic function from 'featherdist'. If featherdist >= cutoffdist
    then this represents a hard cutoff */
HarmonicSwitchingFunction::HarmonicSwitchingFunction(Length cutoffdist,
                                                     Length featherdist)
    : ConcreteProperty<HarmonicSwitchingFunction,SwitchingFunction>()
{
    this->set(cutoffdist, featherdist, cutoffdist, featherdist);
}

/** Construct an harmonic switching function which represents the smoothed
    cutoff, with the cutoff at 'cutoffdist', but with the electrostatic
    interactions smoothed down using a harmonic function from
    'elecfeather' and the vdw interactions smoothed down using
    an harmonic function from 'vdwfeather'. If either feather distance
    is greater than the cutoff, then a hard cutoff for that interaction
    will be used. */
HarmonicSwitchingFunction::HarmonicSwitchingFunction(Length cutoffdist,
                                                     Length elecfeather,
                                                     Length vdwfeather)
     : ConcreteProperty<HarmonicSwitchingFunction,SwitchingFunction>()
{
    this->set(cutdist,elecfeather,cutdist,vdwfeather);
}

/** Construct an harmonic switching function which represents the smoothed
    cutoff, with the electrostatic interactions cutoff at 'eleccutoff', and
    smoothed down using an harmonic function from 'elecfeather', and the
    vdw interactions cutoff at 'vdwcutoff' and smoothed down using an
    harmonic function from 'vdwfeather'. If either feather distance is
    greater than the corresponding cutoff, then a hard cutoff will be
    used for that interaction. */
HarmonicSwitchingFunction::HarmonicSwitchingFunction(Length eleccutoff,
                                                     Length elecfeather,
                                                     Length vdwcutoff,
                                                     Length vdwfeather)
    : ConcreteProperty<HarmonicSwitchingFunction,SwitchingFunction>()
{
    this->set(eleccutoff, elecfeather, vdwcutoff, vdwfeather);
}

/** Copy constructor */
HarmonicSwitchingFunction::HarmonicSwitchingFunction(
                                const HarmonicSwitchingFunction &other)
    : ConcreteProperty<HarmonicSwitchingFunction,SwitchingFunction>(other),
      cut_elec2(other.cut_elec2), norm_elec(other.norm_elec),
      cut_vdw2(other.cut_vdw2), norm_vdw(other.norm_vdw)
{}

/** Destructor */
HarmonicSwitchingFunction::~HarmonicSwitchingFunction()
{}

/** Copy assignment operator */
HarmonicSwitchingFunction& HarmonicSwitchingFunction::operator=(
                                                const HarmonicSwitchingFunction &other)
{
    if (this != &other)
    {
        cut_elec2 = other.cut_elec2;
        norm_elec = other.norm_elec;
        cut_vdw2 = other.cut_vdw2;
        norm_vdw = other.norm_vdw;
        
        SwitchingFunction::operator=(other);
    }

    return *this;
}

/** Comparison operator */
bool HarmonicSwitchingFunction::operator==(const HarmonicSwitchingFunction &other) const
{
    return SwitchingFunction::operator==(other);
}

/** Comparison operator */
bool HarmonicSwitchingFunction::operator!=(const HarmonicSwitchingFunction &other) const
{
    return SwitchingFunction::operator!=(other);
}                                      

/** Return a string representation of this switching function */
QString HarmonicSwitchingFunction::toString() const
{
    return QObject::tr("HarmonicSwitchingFunction( elec = %1 A, vdw = %2 A, "
                       "feather: elec = %3 A, vdw = %4 A )")
                            .arg( electrostaticCutoffDistance().to(angstrom) )
                            .arg( vdwCutoffDistance().to(angstrom) )
                            .arg( electrostaticFeatherDistance().to(angstrom) )
                            .arg( vdwFeatherDistance().to(angstrom) );
}

/** Return the scale factor for the electrostatic interaction for the
    distance 'dist'. This returns;

    dist <= feather_elec           : 1
    feather_elec < dist < cut_elec : (cut_elec^2 - dist^2) / (cut_elec^2 - feather_elec^2)
    dist >= cut_elec               : 0
*/
double HarmonicSwitchingFunction::electrostaticScaleFactor(Length dist) const
{
    if (dist < feather_elec)
        return 1;
    else if (dist > cut_elec)
        return 0;
    else
    {
        return norm_elec * (cut_elec2 - SireMaths::pow_2( double(dist) ));
    }
}

/** Return the scale factor for the vdw interaction for the
    distance 'dist' This returns;

    dist <= feather_vdw           : 1
    feather_vdw < dist < cut_vdw  : (cut_vdw^2 - dist^2) / (cut_vdw^2 - feather_vdw^2)
    dist >= cut_vdw               : 0
*/
double HarmonicSwitchingFunction::vdwScaleFactor(Length dist) const
{
    if (dist < feather_vdw)
        return 1;
    else if (dist > cut_vdw)
        return 0;
    else
    {
        return norm_vdw * (cut_vdw2 - SireMaths::pow_2( double(dist) ));
    }
}

/** This throws an exception as the harmonic switching function has
    a discontinuous first derivative, so is not suitable for force evaluations
    
    \throw SireFF::missing_derivative
*/
double HarmonicSwitchingFunction::dElectrostaticScaleFactor(Length) const
{
    throw SireFF::missing_derivative( QObject::tr(
                "The HarmonicSwitchingFunction has a discontinuous "
                "first derivative so is not suitable for force evaluations. "
                "Try a switching function with a continuous first derivative, "
                "e.g. CharmmSwitchingFunction.") );
}

/** This throws an exception as the harmonic switching function has
    a discontinuous first derivative, so is not suitable for force evaluations
    
    \throw SireFF::missing_derivative
*/
double HarmonicSwitchingFunction::dVDWScaleFactor(Length) const
{
    throw SireFF::missing_derivative( QObject::tr(
                "The HarmonicSwitchingFunction has a discontinuous "
                "first derivative so is not suitable for force evaluations. "
                "Try a switching function with a continuous first derivative, "
                "e.g. CharmmSwitchingFunction.") );
}

const char* HarmonicSwitchingFunction::typeName()
{
    return QMetaType::typeName( qMetaTypeId<HarmonicSwitchingFunction>() );
}

/////////////
///////////// Implementation of CHARMMSwitchingFunction
/////////////

static const RegisterMetaType<CHARMMSwitchingFunction> r_charmm;

/** Internal function used to set the parameters of the cutoff */
void CHARMMSwitchingFunction::set(double cutelec, double featherelec,
                                  double cutvdw, double feathervdw)
{
    cut_elec = std::abs(cutelec);
    cut_elec2 = SireMaths::pow_2(cut_elec);
    feather_elec = qMin( cut_elec, std::abs(featherelec) );
    feather_elec2 = SireMaths::pow_2(feather_elec);

    if (cut_elec != feather_elec)
        //norm_elec = 1.0 / SireMaths::pow_3( cut_elec - feather_elec );
        // JM Feb 11 - above is wrong !
        norm_elec = 1.0 / SireMaths::pow_3( cut_elec2 - feather_elec2 );
    else
        norm_elec = 0;

    cut_vdw = std::abs(cutvdw);
    cut_vdw2 = SireMaths::pow_2(cut_vdw);
    feather_vdw = qMin( cut_vdw, std::abs(feathervdw) );
    feather_vdw2 = SireMaths::pow_2(feather_vdw);

    if (cut_vdw != feather_vdw)
        //norm_vdw = 1.0 / SireMaths::pow_3( cut_vdw - feather_vdw );
        // JM Feb 11
        norm_vdw = 1.0 / SireMaths::pow_3( cut_vdw2 - feather_vdw2 );
    else
        norm_vdw = 0;
        
    cutdist = qMax( cut_elec, cut_vdw );
    featherdist = qMax( feather_elec, feather_vdw );
}

/** Serialise to a binary data stream */
QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds,
                                      const CHARMMSwitchingFunction &charmm)
{
    writeHeader(ds, r_charmm, 1)
          << static_cast<const SwitchingFunction&>(charmm);

    return ds;
}

/** Deserialise from a binary data stream */
QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds,
                                      CHARMMSwitchingFunction &charmm)
{
    VersionID v = readHeader(ds, r_charmm);

    if (v == 1)
    {
        ds >> static_cast<SwitchingFunction&>(charmm);

        charmm.set( charmm.cut_elec, charmm.feather_elec,
                    charmm.cut_vdw, charmm.feather_vdw );
    }
    else
        throw version_error(v, "1", r_charmm, CODELOC);

    return ds;
}

/** Construct a null harmonic switching function (no cutoff) */
CHARMMSwitchingFunction::CHARMMSwitchingFunction()
    : ConcreteProperty<CHARMMSwitchingFunction,SwitchingFunction>()
{
    cut_elec2 = cut_elec;
    feather_elec2 = feather_elec;

    cut_vdw2 = cut_vdw;
    feather_vdw2 = feather_vdw;
    
    norm_elec = 0;
    norm_vdw = 0;
}

/** Construct an harmonic switching function which represents a hard
    cutoff of both the electrostatic and vdw interactions at a distance
    of 'cutoffdist' */
CHARMMSwitchingFunction::CHARMMSwitchingFunction(Length cutoffdist)
     : ConcreteProperty<CHARMMSwitchingFunction,SwitchingFunction>(cutoffdist)
{
    this->set( cutoffdist, cutoffdist, cutoffdist, cutoffdist );
}

/** Construct an harmonic switching function which represents the smoothed
    cutoff, with the cutoff at 'cutoffdist', but smoothed down using an
    harmonic function from 'featherdist'. If featherdist >= cutoffdist
    then this represents a hard cutoff */
CHARMMSwitchingFunction::CHARMMSwitchingFunction(Length cutoffdist,
                                                 Length featherdist)
    : ConcreteProperty<CHARMMSwitchingFunction,SwitchingFunction>(cutoffdist)
{
    this->set(cutoffdist, featherdist, cutoffdist, featherdist);
}

/** Construct an harmonic switching function which represents the smoothed
    cutoff, with the cutoff at 'cutoffdist', but with the electrostatic
    interactions smoothed down using a harmonic function from
    'elecfeather' and the vdw interactions smoothed down using
    an harmonic function from 'vdwfeather'. If either feather distance
    is greater than the cutoff, then a hard cutoff for that interaction
    will be used. */
CHARMMSwitchingFunction::CHARMMSwitchingFunction(Length cutoffdist,
                                                 Length elecfeather,
                                                 Length vdwfeather)
     : ConcreteProperty<CHARMMSwitchingFunction,SwitchingFunction>(cutoffdist)
{
    this->set(cutdist,elecfeather,cutdist,vdwfeather);
}

/** Construct an harmonic switching function which represents the smoothed
    cutoff, with the electrostatic interactions cutoff at 'eleccutoff', and
    smoothed down using an harmonic function from 'elecfeather', and the
    vdw interactions cutoff at 'vdwcutoff' and smoothed down using an
    harmonic function from 'vdwfeather'. If either feather distance is
    greater than the corresponding cutoff, then a hard cutoff will be
    used for that interaction. */
CHARMMSwitchingFunction::CHARMMSwitchingFunction(Length eleccutoff,
                                                 Length elecfeather,
                                                 Length vdwcutoff,
                                                 Length vdwfeather)
    : ConcreteProperty<CHARMMSwitchingFunction,SwitchingFunction>()
{
    this->set(eleccutoff, elecfeather, vdwcutoff, vdwfeather);
}

/** Copy constructor */
CHARMMSwitchingFunction::CHARMMSwitchingFunction(
                                const CHARMMSwitchingFunction &other)
    : ConcreteProperty<CHARMMSwitchingFunction,SwitchingFunction>(other),
      cut_elec2(other.cut_elec2), feather_elec2(other.feather_elec2), 
      norm_elec(other.norm_elec),
      cut_vdw2(other.cut_vdw2), feather_vdw2(other.feather_vdw2),
      norm_vdw(other.norm_vdw)
{}

/** Destructor */
CHARMMSwitchingFunction::~CHARMMSwitchingFunction()
{}

/** Copy assignment operator */
CHARMMSwitchingFunction& CHARMMSwitchingFunction::operator=(
                                                const CHARMMSwitchingFunction &other)
{
    if (this != &other)
    {
        cut_elec2 = other.cut_elec2;
        feather_elec2 = other.feather_elec2;
        norm_elec = other.norm_elec;
        cut_vdw2 = other.cut_vdw2;
        feather_vdw2 = other.feather_vdw2;
        norm_vdw = other.norm_vdw;
        
        SwitchingFunction::operator=(other);
    }

    return *this;
}

/** Comparison operator */
bool CHARMMSwitchingFunction::operator==(const CHARMMSwitchingFunction &other) const
{
    return SwitchingFunction::operator==(other);
}

/** Comparison operator */
bool CHARMMSwitchingFunction::operator!=(const CHARMMSwitchingFunction &other) const
{
    return SwitchingFunction::operator!=(other);
}                                      

/** Return a string representation of this switching function */
QString CHARMMSwitchingFunction::toString() const
{
    return QObject::tr("CHARMMSwitchingFunction( elec = %1 A, vdw = %2 A, "
                       "feather: elec = %3 A, vdw = %4 A )")
                            .arg( electrostaticCutoffDistance().to(angstrom) )
                            .arg( vdwCutoffDistance().to(angstrom) )
                            .arg( electrostaticFeatherDistance().to(angstrom) )
                            .arg( vdwFeatherDistance().to(angstrom) );
}

/** Return the scale factor for the electrostatic interaction for the
    distance 'dist'. This returns;

    dist <= feather_elec           : 1
    feather_elec < dist < cut_elec : (cut_elec^2 - r^2)^3 *
                                       (cut_elec^2 + 2 r ^2 - 3 feather_elec^2 ) /
                                         (cut_elec^2 - feather_elec^2)^3
    JM Feb Correct expresion is  (cut_elec^2 - r^2)^2 *
                                 (cut_elec^2 + 2 r ^2 - 3 feather_elec^2 ) /
                                 (cut_elec^2 - feather_elec^2)^3


    dist >= cut_elec               : 0
*/
double CHARMMSwitchingFunction::electrostaticScaleFactor(Length dist) const
{
    if (dist < feather_elec)
        return 1;
    else if (dist > cut_elec)
        return 0;
    else
    {
        double dist2 = SireMaths::pow_2( double(dist) );
        
	//return SireMaths::pow_3(cut_elec2 - dist2) *
	return SireMaths::pow_2(cut_elec2 - dist2) *
                 (cut_elec2 + 2*dist2 - 3*feather_elec2) * norm_elec;
    }
}

/** Return the scale factor for the vdw interaction for the
    distance 'dist'. This returns;

    dist <= feather_vdw           : 1
    feather_vdw < dist < cut_vdw  : (cut_vdw^2 - r^2)^3 *
                                       (cut_vdw^2 + 2 r ^2 - 3 feather_vdw^2 ) /
                                         (cut_vdw^2 - feather_vdw^2)^3
    dist >= cut_vdw               : 0
*/
double CHARMMSwitchingFunction::vdwScaleFactor(Length dist) const
{
    if (dist < feather_vdw)
        return 1;
    else if (dist > cut_vdw)
        return 0;
    else
    {
        double dist2 = SireMaths::pow_2( double(dist) );
        
        //return SireMaths::pow_3(cut_vdw2 - dist2) *
	return SireMaths::pow_2(cut_vdw2 - dist2) *
                 (cut_vdw2 + 2*dist2 - 3*feather_vdw2) * norm_vdw;
    }
}

/** Return the derivative of the electrostatic scaling factor as
    a function of distance. This returns;
    
    dist <= feather_elec           : 0
    feather_elec < dist < cut_elec : 12 r (cut_elec^2 - r^2)(feather_elec^2 - r^2) /
                                            (cut_elec^2 - feather_elec^2)^3
    dist >= cut_elec               : 0
    
    The derivative is continuous in distance, but the second derivative
    is discontinuous.
*/
double CHARMMSwitchingFunction::dElectrostaticScaleFactor(Length dist) const
{
    if (dist < feather_elec)
        return 0;
    else if (dist > cut_elec)
        return 0;
    else
    {
        double dist2 = SireMaths::pow_2( double(dist) );
        
        return 12 * dist * (cut_elec2 - dist2) * (feather_elec2 - dist2) * norm_elec;
    }
}

/** Return the derivative of the vdw scaling factor as
    a function of distance. This returns;
    
    dist <= feather_vdw           : 0
    feather_vdw < dist < cut_vdw  : 12 r (cut_vdw^2 - r^2)(feather_vdw^2 - r^2) /
                                            (cut_vdw^2 - feather_vdw^2)^3
    dist >= cut_vdw               : 0
    
    The derivative is continuous in distance, but the second derivative
    is discontinuous.
*/
double CHARMMSwitchingFunction::dVDWScaleFactor(Length dist) const
{
    if (dist < feather_vdw)
        return 0;
    else if (dist > cut_vdw)
        return 0;
    else
    {
        double dist2 = SireMaths::pow_2( double(dist) );
        
        return 12 * dist * (cut_vdw2 - dist2) * (feather_vdw2 - dist2) * norm_vdw;
    }
}

const char* CHARMMSwitchingFunction::typeName()
{
    return QMetaType::typeName( qMetaTypeId<CHARMMSwitchingFunction>() );
}
