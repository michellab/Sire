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

#ifndef SIREMM_CLJCOMPONENT_H
#define SIREMM_CLJCOMPONENT_H

#include "SireFF/ffcomponent.h"
#include "SireFF/ff.h"

#ifdef SIRE_USE_SSE
    #ifdef __SSE__
        #include <emmintrin.h>   // CONDITIONAL_INCLUDE
    #else
        #undef SIRE_USE_SSE
    #endif
#endif

SIRE_BEGIN_HEADER

namespace SireMM
{
class CoulombComponent;
class LJComponent;
class CLJComponent;
}

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::CoulombComponent&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::CoulombComponent&);

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::LJComponent&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::LJComponent&);

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::CLJComponent&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::CLJComponent&);

namespace SireFF
{
class FF;
}

namespace SireMM
{

using SireFF::FF;
using SireFF::FFName;

class CoulombComponent;
class LJComponent;

typedef SireFF::ComponentEnergy<CoulombComponent> CoulombEnergy;
typedef SireFF::ComponentEnergy<LJComponent> LJEnergy;

class CLJEnergy;

/** This class represents a Coulomb component of a forcefield */
class SIREMM_EXPORT CoulombComponent : public SireFF::FFComponent
{
public:
    CoulombComponent(const FFName &ffname = FFName());
    CoulombComponent(const FFName &ffname, const QString &suffix);
    
    CoulombComponent(const SireCAS::Symbol &symbol);
    
    CoulombComponent(const CoulombComponent &other);
    
    ~CoulombComponent();
    
    static const char* typeName();
    
    const char* what() const
    {
        return CoulombComponent::typeName();
    }
    
    CoulombComponent* clone() const
    {
        return new CoulombComponent(*this);
    }
    
    const CoulombComponent& total() const
    {
        return *this;
    }

    void setEnergy(FF &ff, const CoulombEnergy &ljnrg) const;
    void changeEnergy(FF &ff, const CoulombEnergy &ljnrg) const;
    
    SireCAS::Symbols symbols() const
    {
        return *this;
    }
};

/** This class represents a LJ component of a forcefield */
class SIREMM_EXPORT LJComponent : public SireFF::FFComponent
{
public:
    LJComponent(const FFName &ffname = FFName());
    LJComponent(const FFName &ffname, const QString &suffix);
    
    LJComponent(const SireCAS::Symbol &symbol);
    
    LJComponent(const LJComponent &other);
    
    ~LJComponent();
    
    static const char* typeName();
    
    const char* what() const
    {
        return LJComponent::typeName();
    }
    
    LJComponent* clone() const
    {
        return new LJComponent(*this);
    }
    
    const LJComponent& total() const
    {
        return *this;
    }

    void setEnergy(FF &ff, const LJEnergy &ljnrg) const;
    void changeEnergy(FF &ff, const LJEnergy &ljnrg) const;
    
    SireCAS::Symbols symbols() const
    {
        return *this;
    }
};

/** This class represents the sum of the coulomb and LJ components
    of the forcefield */
class SIREMM_EXPORT CLJComponent : public SireFF::FFComponent
{

friend QDataStream& ::operator<<(QDataStream&, const CLJComponent&);
friend QDataStream& ::operator>>(QDataStream&, CLJComponent&);

public:
    CLJComponent(const FFName &name = FFName());
    CLJComponent(const FFName &name, const QString &suffix);
    
    CLJComponent(const SireCAS::Symbol &symbol);
    
    CLJComponent(const CLJComponent &other);
    
    ~CLJComponent();
    
    const CoulombComponent& coulomb() const
    {
        return coul_component;
    }
    
    const LJComponent& lj() const
    {
        return lj_component;
    }
    
    const CLJComponent& total() const
    {
        return *this;
    }
    
    static const char* typeName();
    
    const char* what() const
    {
        return CLJComponent::typeName();
    }
    
    CLJComponent* clone() const
    {
        return new CLJComponent(*this);
    }

    void setEnergy(FF &ff, const CLJEnergy &cljnrg) const;
    void changeEnergy(FF &ff, const CLJEnergy &cljnrg) const;
    
    SireCAS::Symbols symbols() const;

protected:
    /** The coulomb component */
    CoulombComponent coul_component;
    
    /** The LJ component */
    LJComponent lj_component;
};

/** This class holds the coulomb and Lennard-Jones (LJ) components
    of the energy.
*/
class SIREMM_EXPORT CLJEnergy
{
public:
    typedef CLJComponent Components;

    CLJEnergy(double cnrg=0, double ljnrg=0)
    {
        #ifdef SIRE_USE_SSE
        nrgs = _mm_setr_pd(cnrg, ljnrg);
        #else
        icnrg = cnrg;
        iljnrg = ljnrg;
        #endif
    }
    
    CLJEnergy(const CLJEnergy &other)
          #ifdef SIRE_USE_SSE
          : nrgs(other.nrgs)
          #else
          : icnrg(other.icnrg), iljnrg(other.iljnrg)
          #endif
    {}
    
    ~CLJEnergy()
    {}
    
    static const char* typeName()
    {
        return "SireMM::CLJEnergy";
    }
    
    const char* what() const
    {
        return CLJEnergy::typeName();
    }
    
    CLJEnergy& operator+=(const CLJEnergy &other)
    {
        #ifdef SIRE_USE_SSE
        nrgs = _mm_add_pd( nrgs, other.nrgs );
        #else
        icnrg += other.icnrg;
        iljnrg += other.iljnrg;
        #endif
        
        return *this;
    }
    
    CLJEnergy& operator-=(const CLJEnergy &other)
    {
        #ifdef SIRE_USE_SSE
        nrgs = _mm_sub_pd(nrgs, other.nrgs);
        #else
        icnrg -= other.icnrg;
        iljnrg -= other.iljnrg;
        #endif
        
        return *this;
    }
    
    CLJEnergy operator+(const CLJEnergy &other) const
    {
        CLJEnergy ret(*this);
        ret += other;
        return ret;
    }
    
    CLJEnergy operator-(const CLJEnergy &other) const
    {
        CLJEnergy ret(*this);
        ret -= other;
        return ret;
    }
    
    Components components() const
    {
        return Components();
    }
    
    double coulomb() const
    {
        #ifdef SIRE_USE_SSE
        return *((const double*)&nrgs);
        #else
        return icnrg;
        #endif
    }
    
    double lj() const
    {
        #ifdef SIRE_USE_SSE
        return *( ((const double*)&nrgs) + 1 );
        #else
        return iljnrg;
        #endif
    }
    
    double total() const
    {
        return coulomb() + lj();
    }
    
    double component(const CoulombComponent&) const
    {
        return coulomb();
    }
    
    double component(const LJComponent&) const
    {
        return lj();
    }
    
    double component(const CLJComponent&) const
    {
        return total();
    }
    
    operator double() const
    {
        //return the total energy
        return total();
    }
    
    operator SireUnits::Dimension::MolarEnergy() const
    {
        return SireUnits::Dimension::MolarEnergy(total());
    }
    
    operator CoulombEnergy() const
    {
        return CoulombEnergy(coulomb());
    }

    operator LJEnergy() const
    {
        return LJEnergy(lj());
    }

private:
    /** The coulomb and LJ components of the energy */
    #ifdef SIRE_USE_SSE
    __m128d nrgs;
    #else
    double icnrg, iljnrg;
    #endif
};

} // end of namespace SireMM

SIRE_EXPOSE_CLASS( SireMM::CoulombComponent )
SIRE_EXPOSE_CLASS( SireMM::LJComponent )
SIRE_EXPOSE_CLASS( SireMM::CLJComponent )

Q_DECLARE_METATYPE( SireMM::CoulombComponent )
Q_DECLARE_METATYPE( SireMM::LJComponent )
Q_DECLARE_METATYPE( SireMM::CLJComponent )

SIRE_END_HEADER

#endif
