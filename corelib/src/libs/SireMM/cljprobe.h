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

#ifndef SIREMM_CLJPROBE_H
#define SIREMM_CLJPROBE_H

#include "ljparameter.h"

#include "SireFF/probe.h"

#include "SireUnits/dimensions.h"

SIRE_BEGIN_HEADER

namespace SireMM
{
class CoulombProbe;
class LJProbe;
class CLJProbe;
}

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::CoulombProbe&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::CoulombProbe&);

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::CLJProbe&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::CLJProbe&);

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::LJProbe&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::LJProbe&);

namespace SireMM
{

/** This is a probe charge that can be used to probe the 
    coulomb electric field or potential in a forcefield */
class SIREMM_EXPORT CoulombProbe 
        : public SireBase::ConcreteProperty<CoulombProbe,SireFF::Probe>
{

friend QDataStream& ::operator<<(QDataStream&, const CoulombProbe&);
friend QDataStream& ::operator>>(QDataStream&, CoulombProbe&);

public:
    CoulombProbe();
    CoulombProbe(SireUnits::Dimension::Charge charge);
    CoulombProbe(const CLJProbe &cljprobe);
    CoulombProbe(const SireFF::Probe &probe);
    
    CoulombProbe(const CoulombProbe &other);
    
    ~CoulombProbe();
    
    CoulombProbe& operator=(const CoulombProbe &other);
    
    bool operator==(const CoulombProbe &other) const;
    bool operator!=(const CoulombProbe &other) const;
    
    static const char* typeName();
    
    SireUnits::Dimension::Charge charge() const;
    double reducedCharge() const;

private:
    /** The actual charge on the probe */
    SireUnits::Dimension::Charge chg;
    
    /** The charge in reduced units */
    double reduced_chg;
};

/** This is a probe charge that can be used to probe the 
    LJ field or potential in a forcefield */
class SIREMM_EXPORT LJProbe 
        : public SireBase::ConcreteProperty<LJProbe,SireFF::Probe>
{

friend QDataStream& ::operator<<(QDataStream&, const LJProbe&);
friend QDataStream& ::operator>>(QDataStream&, LJProbe&);

public:
    LJProbe();
    LJProbe(const LJParameter &ljparam);
    LJProbe(const CLJProbe &cljprobe);
    LJProbe(const SireFF::Probe &probe);
    
    LJProbe(const LJProbe &other);
    
    ~LJProbe();
    
    LJProbe& operator=(const LJProbe &other);
    
    bool operator==(const LJProbe &other) const;
    bool operator!=(const LJProbe &other) const;
    
    static const char* typeName();
    
    const LJParameter& lj() const;

private:
    /** The LJ parameter that represents the probe */
    LJParameter ljparam;
};

/** This is a probe used to probe the coulomb+LJ field
    or potential at points in a forcefield */
class SIREMM_EXPORT CLJProbe 
        : public SireBase::ConcreteProperty<CLJProbe,SireFF::Probe>
{

friend QDataStream& ::operator<<(QDataStream&, const CLJProbe&);
friend QDataStream& ::operator>>(QDataStream&, CLJProbe&);

public:
    CLJProbe();
    CLJProbe(SireUnits::Dimension::Charge charge);
    CLJProbe(const LJParameter &ljparam);
    CLJProbe(SireUnits::Dimension::Charge charge, const LJParameter &ljparam);

    CLJProbe(const CoulombProbe &probe);
    CLJProbe(const LJProbe &probe);
    CLJProbe(const SireFF::Probe &probe);
    
    CLJProbe(const CLJProbe &cljprobe);
    
    ~CLJProbe();
    
    CLJProbe& operator=(const CLJProbe &other);
    
    bool operator==(const CLJProbe &other) const;
    bool operator!=(const CLJProbe &other) const;
    
    static const char* typeName();
    
    SireUnits::Dimension::Charge charge() const;
    double reducedCharge() const;
    
    const LJParameter& lj() const;

private:
    /** The LJ parameter that represents the probe */
    LJParameter ljparam;

    /** The charge on the probe */
    SireUnits::Dimension::Charge chg;
    
    /** The charge in reduced units */
    double reduced_chg;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Return the charge on the probe */
inline SireUnits::Dimension::Charge CoulombProbe::charge() const
{
    return chg;
}

/** Return the reduced charge on the probe */
inline double CoulombProbe::reducedCharge() const
{
    return reduced_chg;
}

/** Return the LJ parameters for this probe */
inline const LJParameter& LJProbe::lj() const
{
    return ljparam;
}

/** Return the charge on the probe */
inline SireUnits::Dimension::Charge CLJProbe::charge() const
{
    return chg;
}

/** Return the reduced charge on the probe */
inline double CLJProbe::reducedCharge() const
{
    return reduced_chg;
}

/** Return the LJ parameters for this probe */
inline const LJParameter& CLJProbe::lj() const
{
    return ljparam;
}


#endif // SIRE_SKIP_INLINE_FUNCTIONS

}

Q_DECLARE_METATYPE( SireMM::CoulombProbe )
Q_DECLARE_METATYPE( SireMM::LJProbe )
Q_DECLARE_METATYPE( SireMM::CLJProbe )

SIRE_EXPOSE_CLASS( SireMM::CoulombProbe )
SIRE_EXPOSE_CLASS( SireMM::LJProbe )
SIRE_EXPOSE_CLASS( SireMM::CLJProbe )

SIRE_END_HEADER

#endif
