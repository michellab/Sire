/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2014  Christopher Woods
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

#ifndef SIREMM_MULTICLJCOMPONENT_H
#define SIREMM_MULTICLJCOMPONENT_H

#include "SireFF/ffcomponent.h"
#include "SireMM/cljcomponent.h"

#include <boost/tuple/tuple.hpp>

SIRE_BEGIN_HEADER

namespace SireMM
{
class MultiCLJComponent;
class MultiCLJEnergy;
}

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::MultiCLJComponent&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::MultiCLJComponent&);

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::MultiCLJEnergy&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::MultiCLJEnergy&);

namespace SireMM
{

class MultiCLJEnergy;

boost::tuple<QString,QString> getSubscriptedProperty(QString name);

/** This class is used to hold the symbols for CLJ forcefields
    that allow multiple CLJ functions to be indexed by key
    
    @author Christopher Woods
*/
class SIREMM_EXPORT MultiCLJComponent : public SireFF::FFComponent
{

friend QDataStream& ::operator<<(QDataStream&, const MultiCLJComponent&);
friend QDataStream& ::operator>>(QDataStream&, MultiCLJComponent&);

public:
    MultiCLJComponent(const FFName &name = FFName());
    
    MultiCLJComponent(const MultiCLJComponent &other);
    
    ~MultiCLJComponent();
    
    MultiCLJComponent rename(const FFName &name) const;
    
    QString toString() const;
    
    MultiCLJComponent& operator=(const MultiCLJComponent &other);
    
    bool operator==(const MultiCLJComponent &other) const;
    bool operator!=(const MultiCLJComponent &other) const;
    
    const CoulombComponent& coulomb() const;
    const CoulombComponent& coulomb(QString key) const;

    const LJComponent& lj() const;
    const LJComponent& lj(QString key) const;
    
    const CLJComponent& total() const;
    const CLJComponent& total(QString key) const;
    
    static const char* typeName();
    
    const char* what() const;
    
    MultiCLJComponent* clone() const;

    void setEnergy(FF &ff, const MultiCLJEnergy &cljnrg) const;
    void changeEnergy(FF &ff, const MultiCLJEnergy &cljnrg) const;
    
    SireCAS::Symbols symbols() const;

    int add(QString key);
    int remove(QString key);
    void removeAll();

    bool hasKey(QString key) const;

    int indexOf(QString key) const;

    QStringList keys() const;

    int count() const;
    int size() const;
    int nKeys() const;

private:
    void assertValidKey(QString key) const;

    /** Array of the CLJ components in the order they appear in the forcefield */
    QVector<CLJComponent> comps;
    
    /** Hash mapping from key name to index in the array */
    QHash<QString,quint32> key_to_idx;
};

/** This class is used during a CLJ calculation to hold all of the 
    coulomb and LJ energies, and to then update them in the forcefield
    
    @author Christopher Woods
*/
class SIREMM_EXPORT MultiCLJEnergy : public CLJEnergy
{

friend QDataStream& ::operator<<(QDataStream&, const MultiCLJEnergy&);
friend QDataStream& ::operator>>(QDataStream&, MultiCLJEnergy&);

public:
    typedef MultiCLJComponent Components;

    MultiCLJEnergy(double cnrg=0, double ljnrg=0) : CLJEnergy(cnrg, ljnrg)
    {}
    
    MultiCLJEnergy(const CLJEnergy &cljnrg) : CLJEnergy(cljnrg)
    {}

    MultiCLJEnergy(const QVector<double> &coul_nrgs, const QVector<double> &lj_nrgs)
           : CLJEnergy()
    {
        double cnrg = 0;
        double ljnrg = 0;
    
        if (not coul_nrgs.isEmpty())
        {
            cnrg = coul_nrgs.at(0);
        }
        
        if (not lj_nrgs.isEmpty())
        {
            ljnrg = lj_nrgs.at(0);
        }
        
        CLJEnergy::operator=( CLJEnergy(cnrg,ljnrg) );
        
        cnrgs = coul_nrgs;
        ljnrgs = lj_nrgs;
    }
    
    MultiCLJEnergy(const MultiCLJEnergy &other)
           : CLJEnergy(other), cnrgs(other.cnrgs), ljnrgs(other.ljnrgs)
    {}
    
    ~MultiCLJEnergy()
    {}
    
    static const char* typeName()
    {
        return "SireMM::MultiCLJEnergy";
    }
    
    const char* what() const
    {
        return MultiCLJEnergy::typeName();
    }
    
    MultiCLJEnergy& operator=(const MultiCLJEnergy &other)
    {
        CLJEnergy::operator=(other);
        cnrgs = other.cnrgs;
        ljnrgs = other.ljnrgs;
        return *this;
    }
    
    MultiCLJEnergy& operator+=(const MultiCLJEnergy &other)
    {
        CLJEnergy::operator+=(other);
        
        if (not cnrgs.isEmpty())
        {
            for (int i=0; i<qMin(cnrgs.count(), other.cnrgs.count()); ++i)
            {
                cnrgs.data()[i] += other.cnrgs.data()[i];
            }
        }
        
        if (not ljnrgs.isEmpty())
        {
            for (int i=0; i<qMin(ljnrgs.count(), other.ljnrgs.count()); ++i)
            {
                ljnrgs.data()[i] += other.ljnrgs.data()[i];
            }
        }
        
        return *this;
    }
    
    MultiCLJEnergy& operator-=(const MultiCLJEnergy &other)
    {
        CLJEnergy::operator-=(other);
        
        if (not cnrgs.isEmpty())
        {
            for (int i=0; i<qMin(cnrgs.count(), other.cnrgs.count()); ++i)
            {
                cnrgs.data()[i] -= other.cnrgs.data()[i];
            }
        }
        
        if (not ljnrgs.isEmpty())
        {
            for (int i=0; i<qMin(ljnrgs.count(), other.ljnrgs.count()); ++i)
            {
                ljnrgs.data()[i] -= other.ljnrgs.data()[i];
            }
        }
        
        return *this;
    }
    
    MultiCLJEnergy operator+(const MultiCLJEnergy &other) const
    {
        MultiCLJEnergy ret(*this);
        ret += other;
        return ret;
    }
    
    MultiCLJEnergy operator-(const MultiCLJEnergy &other) const
    {
        MultiCLJEnergy ret(*this);
        ret -= other;
        return ret;
    }
    
    Components components() const
    {
        return Components();
    }
    
    double coulomb() const
    {
        return CLJEnergy::coulomb();
    }
    
    double lj() const
    {
        return CLJEnergy::lj();
    }
    
    double total() const
    {
        return CLJEnergy::total();
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

    double coulomb(quint32 i) const
    {
        if (i == 0)
            return CLJEnergy::coulomb();
        else
        {
            if (i >= cnrgs.count())
                assertValidCoulombIndex(i);
            
            return cnrgs.at(i);
        }
    }
    
    double lj(quint32 i) const
    {
        if (i == 0)
            return CLJEnergy::lj();
        else
        {
            if (i >= ljnrgs.count())
                assertValidLJIndex(i);
            
            return ljnrgs.at(i);
        }
    }
    
    double total(quint32 i) const
    {
        return coulomb(i) + lj(i);
    }
    
    double component(const CoulombComponent&, quint32 i) const
    {
        return coulomb(i);
    }
    
    double component(const LJComponent&, quint32 i) const
    {
        return lj(i);
    }
    
    double component(const CLJComponent&, quint32 i) const
    {
        return total(i);
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
    void assertValidCoulombIndex(quint32 i) const;
    void assertValidLJIndex(quint32 i) const;

    /** The array of energies if multiple CLJ functions were used */
    QVector<double> cnrgs;
    QVector<double> ljnrgs;
};

}

Q_DECLARE_METATYPE( SireMM::MultiCLJComponent )

SIRE_EXPOSE_CLASS( SireMM::MultiCLJComponent )

SIRE_END_HEADER

#endif
