/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2006  Christopher Woods
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

#ifndef SIREMM_SWITCHINGFUNCTION_H
#define SIREMM_SWITCHINGFUNCTION_H

#include "SireBase/sharedpolypointer.hpp"
#include "SireBase/property.h"

#include "SireUnits/dimensions.h"

SIRE_BEGIN_HEADER

namespace SireMM
{
class SwitchingFunction;

class NoCutoff;
class HarmonicSwitchingFunction;
class CHARMMSwitchingFunction;
}

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::SwitchingFunction&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::SwitchingFunction&);

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::NoCutoff&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::NoCutoff&);

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::HarmonicSwitchingFunction&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::HarmonicSwitchingFunction&);

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::CHARMMSwitchingFunction&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::CHARMMSwitchingFunction&);

namespace SireMM
{

using SireUnits::Dimension::Length;

/** This is the virtual base class of all switching functions. These 
    return scale factors based on the supplied distance

    @author Christopher Woods
*/
class SIREMM_EXPORT SwitchingFunction : public SireBase::Property
{

friend SIREMM_EXPORT QDataStream& ::operator<<(QDataStream&, const SwitchingFunction&);
friend SIREMM_EXPORT QDataStream& ::operator>>(QDataStream&, SwitchingFunction&);

public:
    ~SwitchingFunction();

    static const char* typeName()
    {
        return "SireMM::SwitchingFunction";
    }

    virtual SwitchingFunction* clone() const=0;

    /** Return the electrostatic scale factor for the distance 'dist' */
    virtual double electrostaticScaleFactor(Length dist) const=0;

    /** Return the VDW scale factor for the distance 'dist' */
    virtual double vdwScaleFactor(Length dist) const=0;

    /** Return the derivative (gradient) of the electrostatic
        scale factor at the distance 'dist' */
    virtual double dElectrostaticScaleFactor(Length dist) const=0;

    /** Return the derivative (gradient) of the VDW
        scale factor at the distance 'dist' */
    virtual double dVDWScaleFactor(Length dist) const=0;

    /** Return a string representation of this switching function */
    virtual QString toString() const=0;

    Length cutoffDistance() const;
    Length featherDistance() const;

    Length electrostaticCutoffDistance() const;
    Length electrostaticFeatherDistance() const;
    
    Length vdwCutoffDistance() const;
    Length vdwFeatherDistance() const;

    static const NoCutoff& null();

protected:
    SwitchingFunction();
    SwitchingFunction(Length cutdistance);
    SwitchingFunction(Length cutdistance, Length featherdistance);
    SwitchingFunction(Length eleccut, Length elecfeather,
                      Length vdwcut, Length vdwfeather);

    SwitchingFunction(const SwitchingFunction &other);

    SwitchingFunction& operator=(const SwitchingFunction &other);

    bool operator==(const SwitchingFunction &other) const;
    bool operator!=(const SwitchingFunction &other) const;

    /** The maximum cutoff distance - both the electrostatic
        and vdw energies are scaled to zero beyond this distance */
    double cutdist;

    /** The maximum feather distance - feathering of the electrostatic
        and vdw interaction is *not* performed below this distance */
    double featherdist;

    /** The electrostatic cutoff distance */
    double cut_elec;
    /** The electrostatic feather distance */
    double feather_elec;

    /** The vdw cutoff distance */
    double cut_vdw;
    /** The vdw feather distance */
    double feather_vdw;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Return the cutoff distance beyond which both the electrostatic
    and vdw energies are scaled to zero */
SIRE_ALWAYS_INLINE Length SwitchingFunction::cutoffDistance() const
{
    return Length(cutdist);
}

/** Return the feather distance, below which feathering of the 
    electrostatic and vdw interactions is *not* performed */
SIRE_ALWAYS_INLINE Length SwitchingFunction::featherDistance() const
{
    return Length(featherdist);
}

/** Return the distance beyond which the electrostatic interaction
    is not evaluated */
SIRE_ALWAYS_INLINE Length SwitchingFunction::electrostaticCutoffDistance() const
{
    return Length(cut_elec);
}

/** Return the distance below which the electrostatic interaction
    is *not* feathered */
SIRE_ALWAYS_INLINE Length SwitchingFunction::electrostaticFeatherDistance() const
{
    return Length(feather_elec);
}

/** Return the distance beyond which the VDW interaction is not evaluated */
SIRE_ALWAYS_INLINE Length SwitchingFunction::vdwCutoffDistance() const
{
    return Length(cut_vdw);
}

/** Return the distance below which the VDW interaction is 
    *not* feathered */
SIRE_ALWAYS_INLINE Length SwitchingFunction::vdwFeatherDistance() const
{
    return Length(feather_vdw);
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

/** This class implements no cutoffs (e.g. there is no cutoff, 
    and no switching function!).

    @author Christopher Woods
*/
class SIREMM_EXPORT NoCutoff 
        : public SireBase::ConcreteProperty<NoCutoff,SwitchingFunction>
{

friend SIREMM_EXPORT QDataStream& ::operator<<(QDataStream&, const NoCutoff&);
friend SIREMM_EXPORT QDataStream& ::operator>>(QDataStream&, NoCutoff&);

public:
    NoCutoff();

    NoCutoff(const NoCutoff &other);

    ~NoCutoff();

    static const char* typeName();
    
    NoCutoff& operator=(const NoCutoff &other);
    
    bool operator==(const NoCutoff &other) const;
    bool operator!=(const NoCutoff &other) const;
    
    QString toString() const;
    
    double electrostaticScaleFactor(Length dist) const;
    double vdwScaleFactor(Length dist) const;
    
    double dElectrostaticScaleFactor(Length dist) const;
    double dVDWScaleFactor(Length dist) const;
};

/** This class implements harmonic switching functions - these scale the energy
    harmonically down to zero.

    @author Christopher Woods
*/
class SIREMM_EXPORT HarmonicSwitchingFunction 
         : public SireBase::ConcreteProperty<HarmonicSwitchingFunction,
                                             SwitchingFunction>
{

friend SIREMM_EXPORT QDataStream& ::operator<<(QDataStream&, const HarmonicSwitchingFunction&);
friend SIREMM_EXPORT QDataStream& ::operator>>(QDataStream&, HarmonicSwitchingFunction&);

public:
    HarmonicSwitchingFunction();

    HarmonicSwitchingFunction(Length cutoffdist);
    HarmonicSwitchingFunction(Length cutoffdist, Length featherdist);
    HarmonicSwitchingFunction(Length cutoffdist,
                              Length elecfeather, Length vdwfeather);
    HarmonicSwitchingFunction(Length eleccutoff, Length elecfeather,
                              Length vdwcutoff, Length vdwfeather);

    HarmonicSwitchingFunction(const HarmonicSwitchingFunction &other);

    ~HarmonicSwitchingFunction();

    static const char* typeName();

    QString toString() const;

    HarmonicSwitchingFunction& operator=(const HarmonicSwitchingFunction &other);
    
    bool operator==(const HarmonicSwitchingFunction &other) const;
    bool operator!=(const HarmonicSwitchingFunction &other) const;
    
    double electrostaticScaleFactor(Length dist) const;
    double vdwScaleFactor(Length dist) const;

    double dElectrostaticScaleFactor(Length dist) const;
    double dVDWScaleFactor(Length dist) const;
    
protected:
    void set(double cutelec, double featherelec,
             double cutvdw, double feathervdw);

    /** Square of the cutoff distance */
    double cut_elec2;
    /** Normalisation factor for the electrostatic cutoff */
    double norm_elec;

    /** Square of the cutoff distance */
    double cut_vdw2;
    /** Normalisation factor for the vdw cutoff */
    double norm_vdw;
};

/** This class implements the CHARMMM switching function - these scale the energy
    to zero is such a way that the first derivative of the switching
    function is continuous (and therefore this function can be used
    in a dynamics simulation)

    This is the switching function reported in;
    
    Steinbach and Brooks, "New spherical cutoff methods for long-range
                           forces in macromolecular simulaton"
    
    J. Comp. Chem., 15, 7, pp667-683, 1994

    @author Christopher Woods
*/
class SIREMM_EXPORT CHARMMSwitchingFunction 
         : public SireBase::ConcreteProperty<CHARMMSwitchingFunction,
                                             SwitchingFunction>
{

friend SIREMM_EXPORT QDataStream& ::operator<<(QDataStream&, const CHARMMSwitchingFunction&);
friend SIREMM_EXPORT QDataStream& ::operator>>(QDataStream&, CHARMMSwitchingFunction&);

public:
    CHARMMSwitchingFunction();

    CHARMMSwitchingFunction(Length cutoffdist);
    CHARMMSwitchingFunction(Length cutoffdist, Length featherdist);
    CHARMMSwitchingFunction(Length cutoffdist,
                            Length elecfeather, Length vdwfeather);
    CHARMMSwitchingFunction(Length eleccutoff, Length elecfeather,
                            Length vdwcutoff, Length vdwfeather);

    CHARMMSwitchingFunction(const CHARMMSwitchingFunction &other);

    ~CHARMMSwitchingFunction();

    static const char* typeName();

    CHARMMSwitchingFunction& operator=(const CHARMMSwitchingFunction &other);
    
    bool operator==(const CHARMMSwitchingFunction &other) const;
    bool operator!=(const CHARMMSwitchingFunction &other) const;
    
    QString toString() const;
    
    double electrostaticScaleFactor(Length dist) const;
    double vdwScaleFactor(Length dist) const;

    double dElectrostaticScaleFactor(Length dist) const;
    double dVDWScaleFactor(Length dist) const;
    
protected:
    void set(double cutelec, double featherelec,
             double cutvdw, double feathervdw);

    /** Square of the cutoff distance */
    double cut_elec2;
    
    /** Square of the feather distance */
    double feather_elec2;
    
    /** Normalisation factor for the electrostatic cutoff */
    double norm_elec;

    /** Square of the cutoff distance */
    double cut_vdw2;
    
    /** Square of the feather distance */
    double feather_vdw2;
    
    /** Normalisation factor for the vdw cutoff */
    double norm_vdw;
};

typedef SireBase::PropPtr<SwitchingFunction> SwitchFuncPtr;

}

Q_DECLARE_METATYPE(SireMM::NoCutoff)
Q_DECLARE_METATYPE(SireMM::HarmonicSwitchingFunction)
Q_DECLARE_METATYPE(SireMM::CHARMMSwitchingFunction)

SIRE_EXPOSE_CLASS( SireMM::SwitchingFunction )
SIRE_EXPOSE_CLASS( SireMM::NoCutoff )
SIRE_EXPOSE_CLASS( SireMM::HarmonicSwitchingFunction )
SIRE_EXPOSE_CLASS( SireMM::CHARMMSwitchingFunction )

SIRE_EXPOSE_PROPERTY( SireMM::SwitchFuncPtr, SireMM::SwitchingFunction )

SIRE_END_HEADER

#endif
