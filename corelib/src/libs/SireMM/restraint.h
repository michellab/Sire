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

#ifndef SIREMM_RESTRAINT_H
#define SIREMM_RESTRAINT_H

#include "SireBase/property.h"
#include "SireBase/propertymap.h"

#include "SireVol/space.h"

#include "SireCAS/expression.h"
#include "SireCAS/values.h"

#include "SireUnits/dimensions.h"

SIRE_BEGIN_HEADER

namespace SireMM
{
class Restraint;
class Restraint3D;

class ExpressionRestraint3D;

class NullRestraint;
}

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::Restraint&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::Restraint&);

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::Restraint3D&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::Restraint3D&);

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::ExpressionRestraint3D&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::ExpressionRestraint3D&);

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::NullRestraint&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::NullRestraint&);

namespace SireMol
{
class MoleculeData;
class Molecules;
class MolNum;
class MolID;
}

namespace SireFF
{
class MolForceTable;
class ForceTable;
}

namespace SireCAS
{
class Symbol;
class Symbols;
}

namespace SireMM
{

using SireBase::PropertyMap;

using SireMol::MoleculeData;
using SireMol::Molecules;
using SireMol::MolNum;
using SireMol::MolID;

using SireFF::MolForceTable;
using SireFF::ForceTable;

using SireCAS::Symbol;
using SireCAS::Symbols;
using SireCAS::Values;
using SireCAS::Expression;

using SireVol::Space;

typedef SireBase::PropPtr<Restraint> RestraintPtr;
typedef SireBase::PropPtr<Restraint3D> Restraint3DPtr;

/** This is the base class of all restraints. A restraint is a
    function that calculates the energy or force acting on
    a molecule caused by external potential, e.g. a harmonic
    restraining potential, or a solvent cap potential, or
    a dihedral restraint potential
    
    @author Christopher Woods
*/
class SIREMM_EXPORT Restraint : public SireBase::Property
{

friend QDataStream& ::operator<<(QDataStream&, const Restraint&);
friend QDataStream& ::operator>>(QDataStream&, Restraint&);

public:
    Restraint();
    
    Restraint(const Restraint &other);
    
    virtual ~Restraint();
    
    static const char* typeName()
    {
        return "SireMM::Restraint";
    }
    
    virtual Restraint* clone() const=0;
    
    virtual QString toString() const=0;
    
    virtual SireUnits::Dimension::MolarEnergy energy() const=0;

    virtual void update(const MoleculeData &moldata)=0;

    virtual void update(const Molecules &molecules)=0;

    virtual Molecules molecules() const=0;
    
    virtual bool contains(MolNum molnum) const=0;
    virtual bool contains(const MolID &molid) const=0;

    virtual bool usesMoleculesIn(const Molecules &molecules) const=0;

    virtual void setValue(const Symbol &symbol, double value)=0;
    virtual double getValue(const Symbol &symbol) const=0;

    virtual bool hasValue(const Symbol &symbol) const=0;

    virtual Symbols symbols() const=0;
    virtual Symbols builtinSymbols() const=0;
    virtual Symbols userSymbols() const=0;
    
    virtual Values values() const=0;
    virtual Values builtinValues() const=0;
    virtual Values userValues() const=0;

    virtual RestraintPtr differentiate(const Symbol &symbol) const=0;

    static const NullRestraint& null();

protected:
    Restraint& operator=(const Restraint &other);
    
    bool operator==(const Restraint &other) const;
    bool operator!=(const Restraint &other) const;
};

/** This is the base class of all restraints that operate in 3 dimensions,
    and so can thus return the force on the molecule caused by the restraint
    
    @author Christopher Woods
*/
class SIREMM_EXPORT Restraint3D : public Restraint
{

friend QDataStream& ::operator<<(QDataStream&, const Restraint3D&);
friend QDataStream& ::operator>>(QDataStream&, Restraint3D&);

public:
    Restraint3D();
    
    Restraint3D(const Restraint3D &other);
    
    virtual ~Restraint3D();

    static const char* typeName()
    {
        return "SireMM::Restraint3D";
    }

    virtual Restraint3D* clone() const=0;

    virtual void force(MolForceTable &forcetable, double scale_force=1) const=0;
    virtual void force(ForceTable &forcetable, double scale_force=1) const=0;
    
    virtual bool usesMoleculesIn(const ForceTable &forcetable) const=0;
    virtual bool usesMoleculesIn(const Molecules &molecules) const=0;
    
    const Space& space() const;
    
    virtual void setSpace(const Space &space);
    
protected:
    Restraint3D& operator=(const Restraint3D &other);
    
    bool operator==(const Restraint3D &other) const;
    bool operator!=(const Restraint3D &other) const;

private:
    /** The 3D space in which this restraint operates */
    SireVol::SpacePtr spce;
};

/** This is the base class of 3D restraints that use a user-supplied
    energy function to calculate the restraint energy */
class SIREMM_EXPORT ExpressionRestraint3D : public Restraint3D
{

friend QDataStream& ::operator<<(QDataStream&, const ExpressionRestraint3D&);
friend QDataStream& ::operator>>(QDataStream&, ExpressionRestraint3D&);

public:
    ExpressionRestraint3D();
    ExpressionRestraint3D(const Expression &nrg_expression,
                          const Values &values = Values());
    
    ExpressionRestraint3D(const ExpressionRestraint3D &other);
    
    ~ExpressionRestraint3D();

    static const char* typeName()
    {
        return "SireMM::ExpressionRestraint3D";
    }

    QString toString() const;

    const Expression& restraintFunction() const;

    void setValue(const Symbol &symbol, double value);
    double getValue(const Symbol &symbol) const;

    bool hasValue(const Symbol &symbol) const;

    Symbols symbols() const;
    Symbols userSymbols() const;
    
    Values values() const;
    Values userValues() const;

    SireUnits::Dimension::MolarEnergy energy() const;

protected:
    ExpressionRestraint3D& operator=(const ExpressionRestraint3D &other);
    
    bool operator==(const ExpressionRestraint3D &other) const;
    bool operator!=(const ExpressionRestraint3D &other) const;

    void _pvt_setValue(const Symbol &symbol, double value);

private:
    /** The energy expression */
    Expression nrg_expression;
    
    /** All values that are plugged into this expression */
    Values vals;
};

/** This is a null restraint, that does not affect the energy
    or force on any molecule
    
    @author Christopher Woods
*/
class SIREMM_EXPORT NullRestraint 
            : public SireBase::ConcreteProperty<NullRestraint,Restraint3D>
{

friend QDataStream& ::operator<<(QDataStream&, const NullRestraint &other);
friend QDataStream& ::operator>>(QDataStream&, NullRestraint &other);

public:
    NullRestraint();
    
    NullRestraint(const NullRestraint &other);
    
    ~NullRestraint();
    
    NullRestraint& operator=(const NullRestraint &other);
    
    bool operator==(const NullRestraint &other) const;
    bool operator!=(const NullRestraint &other) const;
    
    static const char* typeName();

    QString toString() const;
    
    SireUnits::Dimension::MolarEnergy energy() const;

    void setValue(const Symbol &symbol, double value);
    double getValue(const Symbol &symbol) const;

    bool hasValue(const Symbol &symbol) const;

    Symbols symbols() const;
    Symbols userSymbols() const;
    Symbols builtinSymbols() const;

    Values values() const;
    Values userValues() const;
    Values builtinValues() const;

    RestraintPtr differentiate(const Symbol &symbol) const;

    void force(MolForceTable &forcetable, double scale_force=1) const;
    void force(ForceTable &forcetable, double scale_force=1) const;
    
    void update(const MoleculeData &moldata);
    void update(const Molecules &molecules);

    Molecules molecules() const;
    
    bool contains(MolNum molnum) const;
    bool contains(const MolID &molid) const;
    
    bool usesMoleculesIn(const ForceTable &forcetable) const;
    bool usesMoleculesIn(const Molecules &molecules) const;
};

}

Q_DECLARE_METATYPE( SireMM::NullRestraint )

SIRE_EXPOSE_CLASS( SireMM::Restraint )
SIRE_EXPOSE_CLASS( SireMM::Restraint3D )
SIRE_EXPOSE_CLASS( SireMM::NullRestraint)

SIRE_EXPOSE_PROPERTY( SireMM::RestraintPtr, SireMM::Restraint )
SIRE_EXPOSE_PROPERTY( SireMM::Restraint3DPtr, SireMM::Restraint3D )

SIRE_END_HEADER

#endif
