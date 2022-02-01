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

#ifndef SIREMM_SOFTCLJCOMPONENT_H
#define SIREMM_SOFTCLJCOMPONENT_H

#include "cljcomponent.h"

SIRE_BEGIN_HEADER

namespace SireMM
{
class SoftCLJComponent;
class SoftCLJEnergy;
}

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::SoftCLJComponent&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::SoftCLJComponent&);

namespace SireMM
{

namespace detail
{
/** This is the maximum number of different alpha values that a single
    SoftCLJPotential may use */
const int MAX_ALPHA_VALUES = 6;
}

/** This represents the sum of the coulomb and LJ components
    for a soft-core forcefield, in which multiple soft-core
    alpha values are used (and so multiple coulomb and LJ
    components are available). This combined component gives
    access to each of the individual components, and also
    to the sum of them all

    @author Christopher Woods
*/
class SIREMM_EXPORT SoftCLJComponent : public CLJComponent
{

friend SIREMM_EXPORT QDataStream& ::operator<<(QDataStream&, const SoftCLJComponent&);
friend SIREMM_EXPORT QDataStream& ::operator>>(QDataStream&, SoftCLJComponent&);

public:
    SoftCLJComponent();
    SoftCLJComponent(const FFName &name);
    SoftCLJComponent(const SireCAS::Symbol &symbol);

    SoftCLJComponent(const SoftCLJComponent &other);

    ~SoftCLJComponent();

    SoftCLJComponent& operator=(const SoftCLJComponent &other);

    static const char* typeName();

    const char* what() const
    {
        return SoftCLJComponent::typeName();
    }

    SoftCLJComponent* clone() const
    {
        return new SoftCLJComponent(*this);
    }

    static int nAlphaValues()
    {
        return detail::MAX_ALPHA_VALUES;
    }

    const CoulombComponent& coulomb() const;
    const LJComponent& lj() const;
    const CLJComponent& total() const;

    const CoulombComponent& coulomb(int i) const;
    const LJComponent& lj(int i) const;
    const CLJComponent& total(int i) const;

    void setEnergy(FF &ff, const CLJEnergy &value) const;
    void changeEnergy(FF &ff, const CLJEnergy &delta) const;

    void setEnergy(FF &ff, const SoftCLJEnergy &value) const;
    void changeEnergy(FF &ff, const SoftCLJEnergy &delta) const;

    SireCAS::Symbols symbols() const;

protected:
    void rebuildComponents();

    //the coulomb and LJ components on CLJComponent represent the
    //total energy of all of the soft core components

    /** The individual components for each alpha value in the same
        order as the corresponding alpha values */
    QVector<CLJComponent> alpha_components;
};

/** This class holds all of the energy component values of
    a multi-alpha SoftCLJPotential-derived forcefield

    @author Christopher Woods
*/
class SIREMM_EXPORT SoftCLJEnergy
{
public:
    typedef SoftCLJComponent Components;

    SoftCLJEnergy();

    SoftCLJEnergy(const SoftCLJEnergy &other);

    ~SoftCLJEnergy();

    static const char* typeName()
    {
        return "SireMM::SoftCLJEnergy";
    }

    const char* what() const
    {
        return SoftCLJEnergy::typeName();
    }

    SoftCLJEnergy& operator=(const SoftCLJEnergy &other);

    SoftCLJEnergy& operator+=(const SoftCLJEnergy &other);
    SoftCLJEnergy& operator-=(const SoftCLJEnergy &other);

    SoftCLJEnergy operator+(const SoftCLJEnergy &other) const;
    SoftCLJEnergy operator-(const SoftCLJEnergy &other) const;

    Components components() const
    {
        return Components();
    }

    int count() const
    {
        return detail::MAX_ALPHA_VALUES;
    }

    int size() const
    {
        return this->count();
    }

    void setEnergy(int i, double cnrg, double ljnrg);

    double coulomb() const;
    double lj() const;

    double coulomb(int i) const;
    double lj(int i) const;

    double total() const;
    double total(int i) const;

    double component(const CoulombComponent&) const;
    double component(const LJComponent&) const;

    double component(const CoulombComponent&, int i) const;
    double component(const LJComponent&, int i) const;

    operator double() const
    {
        return total();
    }

    operator SireUnits::Dimension::Energy() const
    {
        return SireUnits::Dimension::Energy(total());
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
    double nrgs[detail::MAX_ALPHA_VALUES][2];
};

}

SIRE_EXPOSE_CLASS( SireMM::SoftCLJComponent )

Q_DECLARE_METATYPE( SireMM::SoftCLJComponent )

SIRE_END_HEADER

#endif
