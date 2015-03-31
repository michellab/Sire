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

#ifndef SIREMM_INTERNALCOMPONENT_H
#define SIREMM_INTERNALCOMPONENT_H

#include <QSet>

#include "SireFF/ffcomponent.h"
#include "SireFF/ff.h"

SIRE_BEGIN_HEADER

namespace SireMM
{
class BondComponent;
class AngleComponent;
class DihedralComponent;

class ImproperComponent;
class UreyBradleyComponent;

class StretchStretchComponent;
class StretchBendComponent;
class BendBendComponent;
class StretchBendTorsionComponent;

class InternalComponent;
}

QDataStream& operator<<(QDataStream&, const SireMM::BondComponent&);
QDataStream& operator>>(QDataStream&, SireMM::BondComponent&);

QDataStream& operator<<(QDataStream&, const SireMM::AngleComponent&);
QDataStream& operator>>(QDataStream&, SireMM::AngleComponent&);

QDataStream& operator<<(QDataStream&, const SireMM::DihedralComponent&);
QDataStream& operator>>(QDataStream&, SireMM::DihedralComponent&);

QDataStream& operator<<(QDataStream&, const SireMM::ImproperComponent&);
QDataStream& operator>>(QDataStream&, SireMM::ImproperComponent&);

QDataStream& operator<<(QDataStream&, const SireMM::UreyBradleyComponent&);
QDataStream& operator>>(QDataStream&, SireMM::UreyBradleyComponent&);

QDataStream& operator<<(QDataStream&, const SireMM::StretchStretchComponent&);
QDataStream& operator>>(QDataStream&, SireMM::StretchStretchComponent&);

QDataStream& operator<<(QDataStream&, const SireMM::StretchBendComponent&);
QDataStream& operator>>(QDataStream&, SireMM::StretchBendComponent&);

QDataStream& operator<<(QDataStream&, const SireMM::BendBendComponent&);
QDataStream& operator>>(QDataStream&, SireMM::BendBendComponent&);

QDataStream& operator<<(QDataStream&, const SireMM::StretchBendTorsionComponent&);
QDataStream& operator>>(QDataStream&, SireMM::StretchBendTorsionComponent&);

QDataStream& operator<<(QDataStream&, const SireMM::InternalComponent&);
QDataStream& operator>>(QDataStream&, SireMM::InternalComponent&);

namespace SireFF
{
class FF;
}

namespace SireMM
{

using SireFF::FF;
using SireFF::FFName;

using SireCAS::Symbol;

typedef SireFF::ComponentEnergy<BondComponent> BondEnergy;
typedef SireFF::ComponentEnergy<AngleComponent> AngleEnergy;
typedef SireFF::ComponentEnergy<DihedralComponent> DihedralEnergy;

typedef SireFF::ComponentEnergy<ImproperComponent> ImproperEnergy;
typedef SireFF::ComponentEnergy<UreyBradleyComponent> UreyBradleyEnergy;

typedef SireFF::ComponentEnergy<StretchStretchComponent> StretchStretchEnergy;
typedef SireFF::ComponentEnergy<StretchBendComponent> StretchBendEnergy;
typedef SireFF::ComponentEnergy<BendBendComponent> BendBendEnergy;
typedef SireFF::ComponentEnergy<StretchBendTorsionComponent> StretchBendTorsionEnergy;

class InternalEnergy;

/** This class represents a Bond component of a forcefield */
class SIREMM_EXPORT BondComponent : public SireFF::FFComponent
{
public:
    BondComponent(const FFName &ffname = FFName());
    BondComponent(const SireCAS::Symbol &symbol);
    
    BondComponent(const BondComponent &other);
    
    ~BondComponent();
    
    static const char* typeName();
    
    const char* what() const
    {
        return BondComponent::typeName();
    }
    
    BondComponent* clone() const
    {
        return new BondComponent(*this);
    }
    
    const BondComponent& total() const
    {
        return *this;
    }

    void setEnergy(FF &ff, const BondEnergy &bondnrg) const;
    void changeEnergy(FF &ff, const BondEnergy &bondnrg) const;

    SireCAS::Symbols symbols() const
    {
        return *this;
    }
};

/** This class represents a Angle component of a forcefield */
class SIREMM_EXPORT AngleComponent : public SireFF::FFComponent
{
public:
    AngleComponent(const FFName &ffname = FFName());
    AngleComponent(const SireCAS::Symbol &symbol);
    
    AngleComponent(const AngleComponent &other);
    
    ~AngleComponent();
    
    static const char* typeName();
    
    const char* what() const
    {
        return AngleComponent::typeName();
    }
    
    AngleComponent* clone() const
    {
        return new AngleComponent(*this);
    }
    
    const AngleComponent& total() const
    {
        return *this;
    }

    void setEnergy(FF &ff, const AngleEnergy &angnrg) const;
    void changeEnergy(FF &ff, const AngleEnergy &angnrg) const;

    SireCAS::Symbols symbols() const
    {
        return *this;
    }
};

/** This class represents a Dihedral component of a forcefield */
class SIREMM_EXPORT DihedralComponent : public SireFF::FFComponent
{
public:
    DihedralComponent(const FFName &ffname = FFName());
    DihedralComponent(const SireCAS::Symbol &symbol);
    
    DihedralComponent(const DihedralComponent &other);
    
    ~DihedralComponent();
    
    static const char* typeName();
    
    const char* what() const
    {
        return DihedralComponent::typeName();
    }
    
    DihedralComponent* clone() const
    {
        return new DihedralComponent(*this);
    }
    
    const DihedralComponent& total() const
    {
        return *this;
    }

    void setEnergy(FF &ff, const DihedralEnergy &dihnrg) const;
    void changeEnergy(FF &ff, const DihedralEnergy &dihnrg) const;

    SireCAS::Symbols symbols() const
    {
        return *this;
    }
};

/** This class represents a Improper component of a forcefield */
class SIREMM_EXPORT ImproperComponent : public SireFF::FFComponent
{
public:
    ImproperComponent(const FFName &ffname = FFName());
    ImproperComponent(const SireCAS::Symbol &symbol);
    
    ImproperComponent(const ImproperComponent &other);
    
    ~ImproperComponent();
    
    static const char* typeName();
    
    const char* what() const
    {
        return ImproperComponent::typeName();
    }
    
    ImproperComponent* clone() const
    {
        return new ImproperComponent(*this);
    }
    
    const ImproperComponent& total() const
    {
        return *this;
    }

    void setEnergy(FF &ff, const ImproperEnergy &impnrg) const;
    void changeEnergy(FF &ff, const ImproperEnergy &impnrg) const;

    SireCAS::Symbols symbols() const
    {
        return *this;
    }
};

/** This class represents a UreyBradley component of a forcefield */
class SIREMM_EXPORT UreyBradleyComponent : public SireFF::FFComponent
{
public:
    UreyBradleyComponent(const FFName &ffname = FFName());
    UreyBradleyComponent(const SireCAS::Symbol &symbol);
    
    UreyBradleyComponent(const UreyBradleyComponent &other);
    
    ~UreyBradleyComponent();
    
    static const char* typeName();
    
    const char* what() const
    {
        return UreyBradleyComponent::typeName();
    }
    
    UreyBradleyComponent* clone() const
    {
        return new UreyBradleyComponent(*this);
    }
    
    const UreyBradleyComponent& total() const
    {
        return *this;
    }

    void setEnergy(FF &ff, const UreyBradleyEnergy &ubnrg) const;
    void changeEnergy(FF &ff, const UreyBradleyEnergy &ubnrg) const;

    SireCAS::Symbols symbols() const
    {
        return *this;
    }
};

/** This class represents a StretchStretch component of a forcefield */
class SIREMM_EXPORT StretchStretchComponent : public SireFF::FFComponent
{
public:
    StretchStretchComponent(const FFName &ffname = FFName());
    StretchStretchComponent(const SireCAS::Symbol &symbol);
    
    StretchStretchComponent(const StretchStretchComponent &other);
    
    ~StretchStretchComponent();
    
    static const char* typeName();
    
    const char* what() const
    {
        return StretchStretchComponent::typeName();
    }
    
    StretchStretchComponent* clone() const
    {
        return new StretchStretchComponent(*this);
    }
    
    const StretchStretchComponent& total() const
    {
        return *this;
    }

    void setEnergy(FF &ff, const StretchStretchEnergy &ssnrg) const;
    void changeEnergy(FF &ff, const StretchStretchEnergy &ssnrg) const;

    SireCAS::Symbols symbols() const
    {
        return *this;
    }
};

/** This class represents a StretchBend component of a forcefield */
class SIREMM_EXPORT StretchBendComponent : public SireFF::FFComponent
{
public:
    StretchBendComponent(const FFName &ffname = FFName());
    StretchBendComponent(const SireCAS::Symbol &symbol);
    
    StretchBendComponent(const StretchBendComponent &other);
    
    ~StretchBendComponent();
    
    static const char* typeName();
    
    const char* what() const
    {
        return StretchBendComponent::typeName();
    }
    
    StretchBendComponent* clone() const
    {
        return new StretchBendComponent(*this);
    }
    
    const StretchBendComponent& total() const
    {
        return *this;
    }

    void setEnergy(FF &ff, const StretchBendEnergy &sbnrg) const;
    void changeEnergy(FF &ff, const StretchBendEnergy &sbnrg) const;

    SireCAS::Symbols symbols() const
    {
        return *this;
    }
};

/** This class represents a BendBend component of a forcefield */
class SIREMM_EXPORT BendBendComponent : public SireFF::FFComponent
{
public:
    BendBendComponent(const FFName &ffname = FFName());
    BendBendComponent(const SireCAS::Symbol &symbol);
    
    BendBendComponent(const BendBendComponent &other);
    
    ~BendBendComponent();
    
    static const char* typeName();
    
    const char* what() const
    {
        return BendBendComponent::typeName();
    }
    
    BendBendComponent* clone() const
    {
        return new BendBendComponent(*this);
    }
    
    const BendBendComponent& total() const
    {
        return *this;
    }

    void setEnergy(FF &ff, const BendBendEnergy &bbnrg) const;
    void changeEnergy(FF &ff, const BendBendEnergy &bbnrg) const;

    SireCAS::Symbols symbols() const
    {
        return *this;
    }
};

/** This class represents a StretchBendTorsion component of a forcefield */
class SIREMM_EXPORT StretchBendTorsionComponent : public SireFF::FFComponent
{
public:
    StretchBendTorsionComponent(const FFName &ffname = FFName());
    StretchBendTorsionComponent(const SireCAS::Symbol &symbol);
    
    StretchBendTorsionComponent(const StretchBendTorsionComponent &other);
    
    ~StretchBendTorsionComponent();
    
    static const char* typeName();
    
    const char* what() const
    {
        return StretchBendTorsionComponent::typeName();
    }
    
    StretchBendTorsionComponent* clone() const
    {
        return new StretchBendTorsionComponent(*this);
    }
    
    const StretchBendTorsionComponent& total() const
    {
        return *this;
    }

    void setEnergy(FF &ff, const StretchBendTorsionEnergy &sbtnrg) const;
    void changeEnergy(FF &ff, const StretchBendTorsionEnergy &sbtnrg) const;

    SireCAS::Symbols symbols() const
    {
        return *this;
    }
};

/** This class represents the sum of the internal MM energy
    components (bond, angle, dihedral, Urey-Bradley) */
class SIREMM_EXPORT InternalComponent : public SireFF::FFComponent
{

friend QDataStream& ::operator<<(QDataStream&, const InternalComponent&);
friend QDataStream& ::operator>>(QDataStream&, InternalComponent&);

public:
    InternalComponent(const FFName &name = FFName());
    InternalComponent(const SireCAS::Symbol &symbol);
    
    InternalComponent(const InternalComponent &other);
    
    ~InternalComponent();
    
    const BondComponent& bond() const
    {
        return bond_component;
    }
    
    const AngleComponent& angle() const
    {
        return angle_component;
    }
    
    const DihedralComponent& dihedral() const
    {
        return dihedral_component;
    }
    
    const ImproperComponent& improper() const
    {
        return improper_component;
    }
    
    const UreyBradleyComponent& ureyBradley() const
    {
        return ub_component;
    }
    
    const StretchStretchComponent& stretchStretch() const
    {
        return ss_component;
    }
    
    const StretchBendComponent& stretchBend() const
    {
        return sb_component;
    }
    
    const BendBendComponent& bendBend() const
    {
        return bb_component;
    }
    
    const StretchBendTorsionComponent& stretchBendTorsion() const
    {
        return sbt_component;
    }
    
    const InternalComponent& total() const
    {
        return *this;
    }
    
    static const char* typeName();
    
    const char* what() const
    {
        return InternalComponent::typeName();
    }
    
    InternalComponent* clone() const
    {
        return new InternalComponent(*this);
    }

    void setEnergy(FF &ff, const InternalEnergy &nrg) const;
    void changeEnergy(FF &ff, const InternalEnergy &nrg) const;

    SireCAS::Symbols symbols() const;

protected:
    /** The bond component */
    BondComponent bond_component;
    /** The angle component */
    AngleComponent angle_component;
    /** The dihedral component */
    DihedralComponent dihedral_component;

    /** The improper component */
    ImproperComponent improper_component;
    /** The Urey-Bradley component */
    UreyBradleyComponent ub_component;
    
    /** The stretch-stretch component */
    StretchStretchComponent ss_component;
    /** The stretch-bend component */
    StretchBendComponent sb_component;
    /** The bend-bend component */
    BendBendComponent bb_component;
    /** The stretch-bend-torsion component */
    StretchBendTorsionComponent sbt_component;
};

/** This class holds the complete molecule mechanics internal
    energy (combination of the bond, angle, dihedral and Urey-Bradley
    energies) */
class SIREMM_EXPORT InternalEnergy
{
public:
    typedef InternalComponent Components;

    InternalEnergy(double bondnrg=0, double anglenrg=0, double dihedralnrg=0, 
                   double impropernrg=0, double ubnrg=0,
                   double ssnrg=0, double sbnrg=0, double bbnrg=0, double sbtnrg=0);
    
    InternalEnergy(const InternalEnergy &other);
    
    ~InternalEnergy();
    
    static const char* typeName();
    
    const char* what() const
    {
        return InternalEnergy::typeName();
    }

    InternalEnergy& operator=(const InternalEnergy &other)
    {
        ibndnrg = other.ibndnrg;
        iangnrg = other.iangnrg;
        idihnrg = other.idihnrg;
        
        iimpnrg = other.iimpnrg;
        iubnrg = other.iubnrg;
        
        issnrg = other.issnrg;
        isbnrg = other.isbnrg;
        ibbnrg = other.ibbnrg;
        isbtnrg = other.isbtnrg;
        
        return *this;
    }
    
    InternalEnergy& operator+=(const InternalEnergy &other)
    {
        ibndnrg += other.ibndnrg;
        iangnrg += other.iangnrg;
        idihnrg += other.idihnrg;
        
        iimpnrg += other.iimpnrg;
        iubnrg += other.iubnrg;
        
        issnrg += other.issnrg;
        isbnrg += other.isbnrg;
        ibbnrg += other.ibbnrg;
        isbtnrg += other.isbtnrg;
    
        return *this;
    }
    
    InternalEnergy& operator-=(const InternalEnergy &other)
    {
        ibndnrg -= other.ibndnrg;
        iangnrg -= other.iangnrg;
        idihnrg -= other.idihnrg;
        
        iimpnrg -= other.iimpnrg;
        iubnrg -= other.iubnrg;
        
        issnrg -= other.issnrg;
        isbnrg -= other.isbnrg;
        ibbnrg -= other.ibbnrg;
        isbtnrg -= other.isbtnrg;
    
        return *this;
    }
    
    InternalEnergy operator+(const InternalEnergy &other) const
    {
        return InternalEnergy( ibndnrg + other.ibndnrg,
                               iangnrg + other.iangnrg,
                               idihnrg + other.idihnrg,
                               iimpnrg + other.iimpnrg,
                               iubnrg + other.iubnrg,
                               issnrg + other.issnrg,
                               isbnrg + other.isbnrg,
                               ibbnrg + other.ibbnrg,
                               isbtnrg + other.isbtnrg );
    }
    
    InternalEnergy operator-(const InternalEnergy &other) const
    {
        return InternalEnergy( ibndnrg - other.ibndnrg,
                               iangnrg - other.iangnrg,
                               idihnrg - other.idihnrg,
                               iimpnrg - other.iimpnrg,
                               iubnrg - other.iubnrg,
                               issnrg - other.issnrg,
                               isbnrg - other.isbnrg,
                               ibbnrg - other.ibbnrg,
                               isbtnrg - other.isbtnrg );
    }
    
    InternalEnergy& operator+=(const BondEnergy &bndnrg)
    {
        ibndnrg += bndnrg.total();
        return *this;
    }
    
    InternalEnergy& operator+=(const AngleEnergy &angnrg)
    {
        iangnrg += angnrg.total();
        return *this;
    }
    
    InternalEnergy& operator+=(const DihedralEnergy &dihnrg)
    {
        idihnrg += dihnrg.total();
        return *this;
    }
    
    InternalEnergy& operator+=(const ImproperEnergy &impnrg)
    {
        iimpnrg += impnrg.total();
        return *this;
    }
    
    InternalEnergy& operator+=(const UreyBradleyEnergy &ubnrg)
    {
        iubnrg += ubnrg.total();
        return *this;
    }
    
    InternalEnergy& operator+=(const StretchStretchEnergy &ssnrg)
    {
        issnrg += ssnrg.total();
        return *this;
    }
    
    InternalEnergy& operator+=(const StretchBendEnergy &sbnrg)
    {
        isbnrg += sbnrg.total();
        return *this;
    }
    
    InternalEnergy& operator+=(const BendBendEnergy &bbnrg)
    {
        ibbnrg += bbnrg.total();
        return *this;
    }
    
    InternalEnergy& operator+=(const StretchBendTorsionEnergy sbtnrg)
    {
        isbtnrg += sbtnrg.total();
        return *this;
    }
    
    InternalEnergy& operator-=(const BondEnergy &bndnrg)
    {
        ibndnrg -= bndnrg.total();
        return *this;
    }
    
    InternalEnergy& operator-=(const AngleEnergy &angnrg)
    {
        iangnrg -= angnrg.total();
        return *this;
    }
    
    InternalEnergy& operator-=(const DihedralEnergy &dihnrg)
    {
        idihnrg -= dihnrg.total();
        return *this;
    }
    
    InternalEnergy& operator-=(const ImproperEnergy &impnrg)
    {
        iimpnrg -= impnrg.total();
        return *this;
    }
    
    InternalEnergy& operator-=(const UreyBradleyEnergy &ubnrg)
    {
        iubnrg -= ubnrg.total();
        return *this;
    }
    
    InternalEnergy& operator-=(const StretchStretchEnergy &ssnrg)
    {
        issnrg -= ssnrg.total();
        return *this;
    }
    
    InternalEnergy& operator-=(const StretchBendEnergy &sbnrg)
    {
        isbnrg -= sbnrg.total();
        return *this;
    }
    
    InternalEnergy& operator-=(const BendBendEnergy &bbnrg)
    {
        ibbnrg -= bbnrg.total();
        return *this;
    }
    
    InternalEnergy& operator-=(const StretchBendTorsionEnergy sbtnrg)
    {
        isbtnrg -= sbtnrg.total();
        return *this;
    }
    
    Components components() const
    {
        return Components();
    }
    
    double component(const BondComponent&) const
    {
        return ibndnrg;
    }
    
    double component(const AngleComponent&) const
    {
        return iangnrg;
    }
    
    double component(const DihedralComponent&) const
    {
        return idihnrg;
    }
    
    double component(const ImproperComponent&) const
    {
        return iimpnrg;
    }
    
    double component(const UreyBradleyComponent&) const
    {
        return iubnrg;
    }
    
    double component(const StretchStretchComponent&) const
    {
        return issnrg;
    }
    
    double component(const StretchBendComponent&) const
    {
        return isbnrg;
    }
    
    double component(const BendBendComponent&) const
    {
        return ibbnrg;
    }
    
    double component(const StretchBendTorsionComponent&) const
    {
        return isbtnrg;
    }
    
    double component(const InternalComponent&) const
    {
        return ibndnrg + iangnrg + idihnrg + 
               iimpnrg + iubnrg +
               issnrg + isbnrg + ibbnrg + isbtnrg;
    }
    
    double bond() const
    {
        return ibndnrg;
    }
    
    double angle() const
    {
        return iangnrg;
    }
    
    double dihedral() const
    {
        return idihnrg;
    }
    
    double improper() const
    {
        return iimpnrg;
    }
    
    double ureyBradley() const
    {
        return iubnrg;
    }
    
    double stretchStretch() const
    {
        return issnrg;
    }
    
    double stretchBend() const
    {
        return isbnrg;
    }
    
    double bendBend() const
    {
        return ibbnrg;
    }
    
    double stretchBendTorsion() const
    {
        return isbtnrg;
    }
    
    double total() const
    {
        return ibndnrg + iangnrg + idihnrg + 
               iimpnrg + iubnrg +
               issnrg + isbnrg + ibbnrg + isbtnrg;
    }
    
    operator double() const
    {
        //return the total energy
        return InternalEnergy::total();
    }
    
    operator SireUnits::Dimension::MolarEnergy() const
    {
        return SireUnits::Dimension::MolarEnergy( InternalEnergy::total() );
    }
    
    operator BondEnergy() const
    {
        return BondEnergy(ibndnrg);
    }

    operator AngleEnergy() const
    {
        return AngleEnergy(iangnrg);
    }
    
    operator DihedralEnergy() const
    {
        return DihedralEnergy(idihnrg);
    }
    
    operator ImproperEnergy() const
    {
        return ImproperEnergy(iimpnrg);
    }
    
    operator UreyBradleyEnergy() const
    {
        return UreyBradleyEnergy(iubnrg);
    }

    operator StretchStretchEnergy() const
    {
        return StretchStretchEnergy(issnrg);
    }

    operator StretchBendEnergy() const
    {
        return StretchBendEnergy(isbnrg);
    }
    
    operator BendBendEnergy() const
    {
        return BendBendEnergy(ibbnrg);
    }
    
    operator StretchBendTorsionEnergy() const
    {
        return StretchBendTorsionEnergy(isbtnrg);
    }

private:
    /** All of the component energies */
    double ibndnrg, iangnrg, idihnrg, 
           iimpnrg, iubnrg,
           issnrg, isbnrg, ibbnrg, isbtnrg;
};

} // end of namespace SireMM

SIRE_EXPOSE_CLASS( SireMM::BondComponent )
SIRE_EXPOSE_CLASS( SireMM::AngleComponent )
SIRE_EXPOSE_CLASS( SireMM::DihedralComponent )
SIRE_EXPOSE_CLASS( SireMM::ImproperComponent )
SIRE_EXPOSE_CLASS( SireMM::UreyBradleyComponent )
SIRE_EXPOSE_CLASS( SireMM::StretchStretchComponent )
SIRE_EXPOSE_CLASS( SireMM::StretchBendComponent )
SIRE_EXPOSE_CLASS( SireMM::BendBendComponent )
SIRE_EXPOSE_CLASS( SireMM::StretchBendTorsionComponent )

SIRE_EXPOSE_CLASS( SireMM::InternalComponent )

Q_DECLARE_METATYPE( SireMM::BondComponent )
Q_DECLARE_METATYPE( SireMM::AngleComponent )
Q_DECLARE_METATYPE( SireMM::DihedralComponent )
Q_DECLARE_METATYPE( SireMM::ImproperComponent )
Q_DECLARE_METATYPE( SireMM::UreyBradleyComponent )
Q_DECLARE_METATYPE( SireMM::StretchStretchComponent )
Q_DECLARE_METATYPE( SireMM::StretchBendComponent )
Q_DECLARE_METATYPE( SireMM::BendBendComponent )
Q_DECLARE_METATYPE( SireMM::StretchBendTorsionComponent )
Q_DECLARE_METATYPE( SireMM::InternalComponent )

SIRE_END_HEADER

#endif

