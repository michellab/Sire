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

#include <QSet>

#include "forcefields.h"

#include "ffidx.h"
#include "ffname.h"
#include "ff3d.h"
#include "energytable.h"
#include "forcetable.h"
#include "fieldtable.h"
#include "potentialtable.h"
#include "probe.h"

#include "SireMol/mgnum.h"
#include "SireMol/mgidx.h"

#include "SireMol/moleculeview.h"
#include "SireMol/partialmolecule.h"
#include "SireMol/molecule.h"
#include "SireMol/viewsofmol.h"
#include "SireMol/molecules.h"
#include "SireMol/moleculegroup.h"

#include "SireCAS/identities.h"

#include "SireBase/linktoproperty.h"
#include "SireBase/combineproperties.h"

#include "tostring.h"

#include "SireMol/errors.h"
#include "SireFF/errors.h"
#include "SireBase/errors.h"
#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <boost/shared_ptr.hpp>

#include <tbb/tbb.h>

#include <QDebug>
#include <QTime>

using namespace SireFF;
using namespace SireMol;
using namespace SireBase;
using namespace SireCAS;
using namespace SireStream;

using SireUnits::Dimension::MolarEnergy;

namespace SireFF
{
namespace detail
{

/** This is a private hierarchy of classes that is used just by ForceFields
    to relate a symbol to an energy component, forcefield expression or
    constant */
class FFSymbol
{
public:
    FFSymbol();
    FFSymbol(const Symbol &symbol);
    
    FFSymbol(const FFSymbol &other);
    
    virtual ~FFSymbol();
    
    virtual const char* what() const=0;
    
    virtual void load(QDataStream &ds);
    virtual void save(QDataStream &ds) const;
    
    virtual Expression toExpression() const=0;
    
    virtual bool isConstant() const=0;

    bool isEnergy() const
    {
        return not this->isConstant();
    }
    
    virtual double value(const QHash<Symbol,FFSymbolPtr> &ffsymbols) const=0;
    
    virtual MolarEnergy energy(QVector<FFPtr> &forcefields,
                          const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                          double scale_energy=1) const=0;
    
    virtual void energy(EnergyTable &energytable,
                       QVector<FFPtr> &forcefields,
                       const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                       double scale_energy=1) const=0;
                      
    virtual void force(ForceTable &forcetable,
                       QVector<FFPtr> &forcefields,
                       const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                       double scale_force=1) const=0;
                          
    virtual void field(FieldTable &fieldtable,
                       QVector<FFPtr> &forcefields,
                       const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                       double scale_field=1) const=0;
                          
    virtual void field(FieldTable &fieldtable,
                       const Probe &probe,
                       QVector<FFPtr> &forcefields,
                       const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                       double scale_field=1) const=0;
                          
    virtual void potential(PotentialTable &pottable,
                           QVector<FFPtr> &forcefields,
                           const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                           double scale_potential=1) const=0;
                          
    virtual void potential(PotentialTable &pottable,
                           const Probe &probe,
                           QVector<FFPtr> &forcefields,
                           const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                           double scale_potential=1) const=0;

    const Symbol& symbol() const;

    template<class T>
    bool isA() const
    {
        return dynamic_cast<const T*>(this) != 0;
    }
    
    template<class T>
    const T& asA() const
    {
        return dynamic_cast<const T&>(*this);
    }

    template<class T>
    T& asA()
    {
        return dynamic_cast<T&>(*this);
    }

private:
    /** The symbol that this object represents */
    Symbol s;
};

/** This is an FFSymbol that holds just a single value */
class FFConstantValue : public FFSymbol
{
public:
    FFConstantValue();
    FFConstantValue(const Symbol &symbol, double value);
    
    FFConstantValue(const FFConstantValue &other);
    
    ~FFConstantValue();
    
    const char* what() const
    {
        return "SireFF::FFConstantValue";
    }
    
    void load(QDataStream &ds);
    void save(QDataStream &ds) const;
    
    bool isConstant() const { return true; }
    
    Expression toExpression() const;

    double value() const;
    
    double value(const QHash<Symbol,FFSymbolPtr> &ffsymbols) const;
    
    MolarEnergy energy(QVector<FFPtr> &forcefields,
                       const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                       double scale_energy=1) const;
    
    void energy(EnergyTable &energytable, QVector<FFPtr> &forcefields,
		 const QHash<Symbol,FFSymbolPtr> &ffsymbols,
		 double scale_energy=1) const; 

    void force(ForceTable &forcetable,
               QVector<FFPtr> &forcefields,
               const QHash<Symbol,FFSymbolPtr> &ffsymbols,
               double scale_force=1) const;
    
    void field(FieldTable &fieldtable,
               QVector<FFPtr> &forcefields,
               const QHash<Symbol,FFSymbolPtr> &ffsymbols,
               double scale_field=1) const;
    
    void field(FieldTable &fieldtable,
               const Probe &probe,
               QVector<FFPtr> &forcefields,
               const QHash<Symbol,FFSymbolPtr> &ffsymbols,
               double scale_field=1) const;
    
    void potential(PotentialTable &pottable,
                   QVector<FFPtr> &forcefields,
                   const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                   double scale_potential=1) const;
    
    void potential(PotentialTable &pottable,
                   const Probe &probe,
                   QVector<FFPtr> &forcefields,
                   const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                   double scale_potential=1) const;
    
private:
    /** The value of this symbol */
    double v;
};

/** This is an FFSymbol that holds a constant expression */
class FFConstantExpression : public FFSymbol
{
public:
    FFConstantExpression();
    FFConstantExpression(const Symbol &symbol, const Expression &expression);
    
    FFConstantExpression(const FFConstantExpression &other);
    
    ~FFConstantExpression();
    
    const char* what() const
    {
        return "SireFF::FFConstantExpression";
    }
    
    void load(QDataStream &ds);
    void save(QDataStream &ds) const;

    bool isConstant() const { return true; }
    
    Expression toExpression() const;
    
    void assertNotDepends(const QSet<Symbol> &ffsymbols) const;
    
    double value(const QHash<Symbol,FFSymbolPtr> &ffsymbols) const;
    
    MolarEnergy energy(QVector<FFPtr> &forcefields,
                       const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                       double scale_energy=1) const;
    
    void energy(EnergyTable &energytable, QVector<FFPtr> &forcefields,
		 const QHash<Symbol,FFSymbolPtr> &ffsymbols,
		 double scale_energy=1) const; 

    void force(ForceTable &forcetable,
               QVector<FFPtr> &forcefields,
               const QHash<Symbol,FFSymbolPtr> &ffsymbols,
               double scale_force=1) const;
    
    void field(FieldTable &fieldtable,
               QVector<FFPtr> &forcefields,
               const QHash<Symbol,FFSymbolPtr> &ffsymbols,
               double scale_field=1) const;
    
    void field(FieldTable &fieldtable,
               const Probe &probe,
               QVector<FFPtr> &forcefields,
               const QHash<Symbol,FFSymbolPtr> &ffsymbols,
               double scale_field=1) const;
    
    void potential(PotentialTable &pottable,
                   QVector<FFPtr> &forcefields,
                   const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                   double scale_potential=1) const;
    
    void potential(PotentialTable &pottable,
                   const Probe &probe,
                   QVector<FFPtr> &forcefields,
                   const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                   double scale_potential=1) const;

private:
    Expression expression;
    Symbols syms;
};

/** This is an FFSymbol that holds a single forcefield component */
class FFSymbolFF : public FFSymbol
{
public:
    FFSymbolFF();
    FFSymbolFF(FFIdx ffidx, const Symbol &component);
    
    FFSymbolFF(const FFSymbolFF &other);
    
    ~FFSymbolFF();
    
    static const char* typeName()
    {
        return "SireFF::FFSymbolFF";
    }

    const char* what() const
    {
        return FFSymbolFF::typeName();
    }

    bool isConstant() const { return false; }
    
    void load(QDataStream &ds);
    void save(QDataStream &ds) const;

    Expression toExpression() const;
    
    FFIdx ffIdx() const;
    
    double value(const QHash<Symbol,FFSymbolPtr> &ffsymbols) const;
    
    MolarEnergy energy(QVector<FFPtr> &forcefields,
                       const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                       double scale_energy=1) const;

    void energy(EnergyTable &energytable, QVector<FFPtr> &forcefields,
		 const QHash<Symbol,FFSymbolPtr> &ffsymbols,
		 double scale_energy=1) const;        
  
    void force(ForceTable &forcetable, QVector<FFPtr> &forcefields,
               const QHash<Symbol,FFSymbolPtr> &ffsymbols,
               double scale_force=1) const;
    
    void field(FieldTable &fieldtable,
               QVector<FFPtr> &forcefields,
               const QHash<Symbol,FFSymbolPtr> &ffsymbols,
               double scale_field=1) const;
    
    void field(FieldTable &fieldtable,
               const Probe &probe,
               QVector<FFPtr> &forcefields,
               const QHash<Symbol,FFSymbolPtr> &ffsymbols,
               double scale_field=1) const;
    
    void potential(PotentialTable &pottable,
                   QVector<FFPtr> &forcefields,
                   const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                   double scale_potential=1) const;
    
    void potential(PotentialTable &pottable,
                   const Probe &probe,
                   QVector<FFPtr> &forcefields,
                   const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                   double scale_potential=1) const;

private:
    /** The index of the forcefield that contains this component */
    FFIdx ffidx;
};

/** This is an FFSymbol that represents a complete forcefield expression */
class FFSymbolExpression : public FFSymbol
{
public:
    FFSymbolExpression();
    FFSymbolExpression(const Symbol &symbol, const Expression &expression);
    
    FFSymbolExpression(const FFSymbolExpression &other);
    
    ~FFSymbolExpression();
    
    static const char* typeName()
    {
        return "SireFF::FFSymbolExpression";
    }
    
    const char* what() const
    {
        return FFSymbolExpression::typeName();
    }
    
    void load(QDataStream &ds);
    void save(QDataStream &ds) const;

    bool isConstant() const { return false; }

    Expression toExpression() const;
    
    const Expression& expression() const;
    
    void expandInTermsOf(const QSet<Symbol> &ffsymbols,
                         const QHash<Symbol,FFSymbolPtr> &ffmap);
    
    double value(const QHash<Symbol,FFSymbolPtr> &ffsymbols) const;
    
    MolarEnergy energy(QVector<FFPtr> &forcefields,
                       const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                       double scale_energy=1) const;
    
    void energy(EnergyTable &energytable, QVector<FFPtr> &forcefields,
		 const QHash<Symbol,FFSymbolPtr> &ffsymbols,
		 double scale_energy=1) const; 
              
    void force(ForceTable &forcetable, QVector<FFPtr> &forcefields,
               const QHash<Symbol,FFSymbolPtr> &ffsymbols,
               double scale_force=1) const;
    
    void field(FieldTable &fieldtable,
               QVector<FFPtr> &forcefields,
               const QHash<Symbol,FFSymbolPtr> &ffsymbols,
               double scale_field=1) const;
    
    void field(FieldTable &fieldtable,
               const Probe &probe,
               QVector<FFPtr> &forcefields,
               const QHash<Symbol,FFSymbolPtr> &ffsymbols,
               double scale_field=1) const;
    
    void potential(PotentialTable &pottable,
                   QVector<FFPtr> &forcefields,
                   const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                   double scale_potential=1) const;
    
    void potential(PotentialTable &pottable,
                   const Probe &probe,
                   QVector<FFPtr> &forcefields,
                   const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                   double scale_potential=1) const;

private:
    class Component
    {
    public:
        Component();
        Component(const Expression &scale_factor, const Symbol &component,
                  int ffidx=-1);
        
        Component(const Component &other);
        
        ~Component();
        
        int ffIdx() const;
        int nDependents() const;
        const QVector<Symbol>& dependents() const;
        
        double scalingFactor(const Values &values) const;
        const Expression& scalingExpression() const;
        
        const Symbol& symbol() const;
        
    private:
        /** The symbol for this component */
        Symbol s;
        
        /** The symbols used in the scaling factor */
        QVector<Symbol> deps;
        
        /** The expression for the scaling factor */
        Expression sclfac;
        
        /** The index of the forcefield for this component */
        int ffidx;
    };

    /** The forcefield expression */
    Expression ffexpression;
    
    /** All of the components of this expression */
    QVector<Component> components;
    
    /** All of the components sorted by FFIdx */
    QVector< QVector<Component> > sorted_components;
};

/** This is an FFSymbol that represents just a simple total of the energy
    of the forcefields */
class FFTotalExpression : public FFSymbol
{
public:
    FFTotalExpression();
    FFTotalExpression(const Symbol &symbol);
    
    FFTotalExpression(const FFTotalExpression &other);
    
    ~FFTotalExpression();
    
    static const char* typeName()
    {
        return "SireFF::FFTotalExpression";
    }
    
    const char* what() const
    {
        return FFTotalExpression::typeName();
    }
    
    void load(QDataStream &ds);
    void save(QDataStream &ds) const;

    bool isConstant() const { return false; }
    
    Expression toExpression() const;

    double value(const QHash<Symbol,FFSymbolPtr> &ffsymbols) const;
    
    MolarEnergy energy(QVector<FFPtr> &forcefields,
                       const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                       double scale_energy=1) const;
    
    void energy(EnergyTable &energytable, QVector<FFPtr> &forcefields,
		 const QHash<Symbol,FFSymbolPtr> &ffsymbols,
		 double scale_energy=1) const; 
              
    void force(ForceTable &forcetable, QVector<FFPtr> &forcefields,
               const QHash<Symbol,FFSymbolPtr> &ffsymbols,
               double scale_force=1) const;

    
    void field(FieldTable &fieldtable,
               QVector<FFPtr> &forcefields,
               const QHash<Symbol,FFSymbolPtr> &ffsymbols,
               double scale_field=1) const;
    
    void field(FieldTable &fieldtable,
               const Probe &probe,
               QVector<FFPtr> &forcefields,
               const QHash<Symbol,FFSymbolPtr> &ffsymbols,
               double scale_field=1) const;
    
    void potential(PotentialTable &pottable,
                   QVector<FFPtr> &forcefields,
                   const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                   double scale_potential=1) const;
    
    void potential(PotentialTable &pottable,
                   const Probe &probe,
                   QVector<FFPtr> &forcefields,
                   const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                   double scale_potential=1) const;
};

} // end of namespace detail

} // end of namespace SireFF

using namespace SireFF::detail;

/** Serialise an FFSymbol to a binary datastream */
QDataStream& operator<<(QDataStream &ds, const FFSymbolPtr &ffsymbol)
{
    if (ffsymbol.get() == 0)
    {
        ds << QString::null;
    }
    else
    {
        ds << QString( ffsymbol->what() );
        ffsymbol->save(ds);
    }
    
    return ds;
}

/** Extract an FFSymbol from a binary datastream */
QDataStream& operator>>(QDataStream &ds, FFSymbolPtr &ffsymbol)
{
    QString type_name;
    
    ds >> type_name;
    
    if (type_name.isEmpty())
    {
        ffsymbol.reset();
    }
    else 
    {
        if (type_name == QLatin1String("SireFF::FFSymbolFF"))
        {
            ffsymbol.reset( new FFSymbolFF() );
        }
        else if (type_name == QLatin1String("SireFF::FFConstantValue") or
                 type_name == QLatin1String("SireFF::FFSymbolValue"))
        {
            ffsymbol.reset( new FFConstantValue() );
        }
        else if (type_name == QLatin1String("SireFF::FFConstantExpression"))
        {
            ffsymbol.reset( new FFConstantExpression() );
        }
        else if (type_name == QLatin1String("SireFF::FFSymbolExpression"))
        {
            ffsymbol.reset( new FFSymbolExpression() );
        }
        else if (type_name == QLatin1String("SireFF::FFTotalExpression"))
        {
            ffsymbol.reset( new FFTotalExpression() );
        }
        else
        {
            throw SireError::program_bug( QObject::tr(
                "Internal error with ForceFields - it can't recognise the "
                "FFSymbol type %1.").arg(type_name), CODELOC );
        }
        
        ffsymbol->load(ds);
    }
    
    return ds;
}

///////////
/////////// Implementation of FFSymbol
///////////

FFSymbol::FFSymbol()
{}

FFSymbol::FFSymbol(const Symbol &symbol) : s(symbol)
{}

FFSymbol::FFSymbol(const FFSymbol &other) : s(other.s)
{}

FFSymbol::~FFSymbol()
{}

const Symbol& FFSymbol::symbol() const
{
    return s;
}

void FFSymbol::load(QDataStream &ds)
{
    ds >> s;
}

void FFSymbol::save(QDataStream &ds) const
{
    ds << s;
}

///////////
/////////// Implementation of FFConstantValue
///////////

FFConstantValue::FFConstantValue() : FFSymbol(), v(0)
{}

FFConstantValue::FFConstantValue(const Symbol &symbol, double value) 
              : FFSymbol(symbol), v(value)
{}

FFConstantValue::FFConstantValue(const FFConstantValue &other)
               : FFSymbol(other), v(other.v)
{}

FFConstantValue::~FFConstantValue()
{}

void FFConstantValue::load(QDataStream &ds)
{
    ds >> v;
    FFSymbol::load(ds);
}

void FFConstantValue::save(QDataStream &ds) const
{
    ds << v;
    FFSymbol::save(ds);
}

Expression FFConstantValue::toExpression() const
{
    return Expression(v);
}

double FFConstantValue::value() const
{
    return v;
}

double FFConstantValue::value(const QHash<Symbol,FFSymbolPtr>&) const
{
    return v;
}

MolarEnergy FFConstantValue::energy(QVector<FFPtr> &forcefields,
                                  const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                                  double scale_energy) const
{
    throw SireError::program_bug( QObject::tr(
            "Constant values or expressions do not have an energy, and should "
            "not be called as if they have an energy! %1 == %2")
                .arg(this->symbol().toString()).arg(v), CODELOC );
                
    return MolarEnergy();
}

void FFConstantValue::energy(EnergyTable &energytable,
                          QVector<FFPtr> &forcefields,
                          const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                          double scale_energy) const
{
    throw SireError::program_bug( QObject::tr(
            "Constant values or expressions do not have an energy, and should "
            "not be called as if they have an energy! %1 == %2")
                .arg(this->symbol().toString()).arg(v), CODELOC );
}

void FFConstantValue::force(ForceTable &forcetable,
                          QVector<FFPtr> &forcefields,
                          const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                          double scale_force) const
{
    throw SireError::program_bug( QObject::tr(
            "Constant values or expressions do not cause a force, and should "
            "not be called as if they could cause a force! %1 == %2")
                .arg(this->symbol().toString()).arg(v), CODELOC );
}

void FFConstantValue::field(FieldTable &fieldtable,
                            QVector<FFPtr> &forcefields,
                            const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                            double scale_field) const
{
    throw SireError::program_bug( QObject::tr(
            "Constant values or expressions do not cause a field, and should "
            "not be called as if they could cause a field! %1 == %2")
                .arg(this->symbol().toString()).arg(v), CODELOC );
}

void FFConstantValue::field(FieldTable &fieldtable,
                            const Probe &probe,
                            QVector<FFPtr> &forcefields,
                            const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                            double scale_field) const
{
    FFConstantValue::field(fieldtable, forcefields, ffsymbols, scale_field);
}

void FFConstantValue::potential(PotentialTable &pottable,
                                QVector<FFPtr> &forcefields,
                                const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                                double scale_potential) const
{
    throw SireError::program_bug( QObject::tr(
            "Constant values or expressions do not cause a field, and should "
            "not be called as if they could cause a field! %1 == %2")
                .arg(this->symbol().toString()).arg(v), CODELOC );
}

void FFConstantValue::potential(PotentialTable &pottable,
                                const Probe &probe,
                                QVector<FFPtr> &forcefields,
                                const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                                double scale_potential) const
{
    FFConstantValue::potential(pottable, forcefields, ffsymbols, scale_potential);
}

///////////
/////////// Implementation of FFConstantExpression
///////////

FFConstantExpression::FFConstantExpression() : FFSymbol()
{}

FFConstantExpression::FFConstantExpression(const Symbol &symbol,
                                           const Expression &e) 
                     : FFSymbol(symbol), expression(e), syms(e.symbols())
{}

FFConstantExpression::FFConstantExpression(const FFConstantExpression &other)
                     : FFSymbol(other), expression(other.expression),
                       syms(other.syms)
{}

FFConstantExpression::~FFConstantExpression()
{}

void FFConstantExpression::load(QDataStream &ds)
{
    ds >> expression;
    FFSymbol::load(ds);
    
    syms = expression.symbols();
}

void FFConstantExpression::save(QDataStream &ds) const
{
    ds << expression;
    FFSymbol::save(ds);
}

Expression FFConstantExpression::toExpression() const
{
    return expression;
}

void FFConstantExpression::assertNotDepends(const QSet<Symbol> &ffsyms) const
{
    if (ffsyms.count() <= syms.count())
    {
        foreach (const Symbol &ffsym, ffsyms)
        {
            if (syms.contains(ffsym))
                throw SireError::incompatible_error( QObject::tr(
                    "The constant expression %1 == %2 is not really constant as it "
                    "depends on the energy expression with symbol %3.")
                        .arg(this->symbol().toString(), expression.toString())
                        .arg(ffsym.toString()), CODELOC );
        }
    }
    else
    {
        foreach (const Symbol &sym, syms)
        {
            if (ffsyms.contains(sym))
                throw SireError::incompatible_error( QObject::tr(
                    "The constant expression %1 == %2 is not really constant as it "
                    "depends on the energy expression with symbol %3.")
                        .arg(this->symbol().toString(), expression.toString())
                        .arg(sym.toString()), CODELOC );
        }
    }
}

double FFConstantExpression::value(const QHash<Symbol,FFSymbolPtr> &ffsymbols) const
{
    if (not syms.isEmpty())
    {
        Values vals;
        vals.reserve(syms.count());
        
        foreach (const Symbol &symbol, syms)
        {
            if (ffsymbols.contains(symbol))
                vals.set(symbol, ffsymbols[symbol]->value(ffsymbols));
        }
        
        return expression(vals);
    }
    else
        return expression( Values() );
}

MolarEnergy FFConstantExpression::energy(QVector<FFPtr> &forcefields,
                                         const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                                         double scale_energy) const
{
    throw SireError::program_bug( QObject::tr(
            "Constant values or expressions do not have an energy, and should "
            "not be called as if they have an energy! %1 == %2")
                .arg(this->symbol().toString(), expression.toString()), CODELOC );
                
    return MolarEnergy();
}

void FFConstantExpression::energy(EnergyTable &energytable,
                                 QVector<FFPtr> &forcefields,
                                 const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                                 double scale_energy) const
{
    throw SireError::program_bug( QObject::tr(
            "Constant values or expressions do not have an energy, and should "
            "not be called as if they have an energy! %1 == %2")
                .arg(this->symbol().toString(), expression.toString()), CODELOC );
}

void FFConstantExpression::force(ForceTable &forcetable,
                                 QVector<FFPtr> &forcefields,
                                 const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                                 double scale_force) const
{
    throw SireError::program_bug( QObject::tr(
            "Constant values or expressions do not cause a force, and should "
            "not be called as if they could cause a force! %1 == %2")
                .arg(this->symbol().toString(), expression.toString()), CODELOC );
}

void FFConstantExpression::field(FieldTable &fieldtable,
                                 QVector<FFPtr> &forcefields,
                                 const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                                 double scale_field) const
{
    throw SireError::program_bug( QObject::tr(
            "Constant values or expressions do not cause a field, and should "
            "not be called as if they could cause a field! %1 == %2")
                .arg(this->symbol().toString()).arg(expression.toString()), CODELOC );
}

void FFConstantExpression::field(FieldTable &fieldtable,
                                 const Probe &probe,
                                 QVector<FFPtr> &forcefields,
                                 const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                                 double scale_field) const
{
    FFConstantExpression::field(fieldtable, forcefields, ffsymbols, scale_field);
}

void FFConstantExpression::potential(PotentialTable &pottable,
                                     QVector<FFPtr> &forcefields,
                                     const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                                     double scale_potential) const
{
    throw SireError::program_bug( QObject::tr(
            "Constant values or expressions do not cause a field, and should "
            "not be called as if they could cause a field! %1 == %2")
                .arg(this->symbol().toString()).arg(expression.toString()), CODELOC );
}

void FFConstantExpression::potential(PotentialTable &pottable,
                                     const Probe &probe,
                                     QVector<FFPtr> &forcefields,
                                     const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                                     double scale_potential) const
{
    FFConstantExpression::potential(pottable, forcefields, ffsymbols, scale_potential);
}

///////////
/////////// Implementation of FFSymbolFF
///////////

FFSymbolFF::FFSymbolFF() : FFSymbol(), ffidx(0)
{}

FFSymbolFF::FFSymbolFF(FFIdx ffindex, const Symbol &component)
           : FFSymbol(component), ffidx(ffindex)
{}

FFSymbolFF::FFSymbolFF(const FFSymbolFF &other)
           : FFSymbol(other), ffidx(other.ffidx)
{}

FFSymbolFF::~FFSymbolFF()
{}

void FFSymbolFF::load(QDataStream &ds)
{
    ds >> ffidx;
    FFSymbol::load(ds);
}

void FFSymbolFF::save(QDataStream &ds) const
{
    ds << ffidx;
    FFSymbol::save(ds);
}

FFIdx FFSymbolFF::ffIdx() const
{
    return ffidx;
}

Expression FFSymbolFF::toExpression() const
{
    return Expression(this->symbol());
}

double FFSymbolFF::value(const QHash<Symbol,FFSymbolPtr>&) const
{
    throw SireError::program_bug( QObject::tr(
        "There is no constant value associated with a forcefield (%1, %2)")
            .arg(ffidx).arg(this->symbol().toString()), CODELOC );
            
    return 0;
}

MolarEnergy FFSymbolFF::energy(QVector<FFPtr> &forcefields,
                               const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                               double scale_energy) const
{
    if (scale_energy != 0)
        return scale_energy * forcefields[ffidx].edit().energy(this->symbol());
    else
        return MolarEnergy(0);
}

void FFSymbolFF::energy(EnergyTable &energytable, QVector<FFPtr> &forcefields,
                       const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                       double scale_energy) const
{

    FFPtr &ffield = forcefields[ffidx];
    
    if (not ffield->isA<FF3D>())
        throw SireFF::missing_derivative( QObject::tr(
            "The forcefield of type %1 does not inherit from FF3D so does "
            "not provide an energy function.")
                .arg(ffield->what()), CODELOC );

    ffield.edit().asA<FF3D>().energy(energytable, this->symbol(), scale_energy);
}

void FFSymbolFF::force(ForceTable &forcetable, QVector<FFPtr> &forcefields,
                       const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                       double scale_force) const
{
    FFPtr &ffield = forcefields[ffidx];
    
    if (not ffield->isA<FF3D>())
        throw SireFF::missing_derivative( QObject::tr(
            "The forcefield of type %1 does not inherit from FF3D so does "
            "not provide a force function.")
                .arg(ffield->what()), CODELOC );

    ffield.edit().asA<FF3D>().force(forcetable, this->symbol(), scale_force);
}

void FFSymbolFF::field(FieldTable &fieldtable, QVector<FFPtr> &forcefields,
                       const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                       double scale_field) const
{
    FFPtr &ffield = forcefields[ffidx];
    
    if (not ffield->isA<FF3D>())
        throw SireFF::missing_derivative( QObject::tr(
            "The forcefield of type %1 does not inherit from FF3D so does "
            "not provide a force function.")
                .arg(ffield->what()), CODELOC );

    ffield.edit().asA<FF3D>().field(fieldtable, this->symbol(), scale_field);
}

void FFSymbolFF::field(FieldTable &fieldtable, const Probe &probe, 
                       QVector<FFPtr> &forcefields,
                       const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                       double scale_field) const
{
    FFPtr &ffield = forcefields[ffidx];
    
    if (not ffield->isA<FF3D>())
        throw SireFF::missing_derivative( QObject::tr(
            "The forcefield of type %1 does not inherit from FF3D so does "
            "not provide a force function.")
                .arg(ffield->what()), CODELOC );

    ffield.edit().asA<FF3D>().field(fieldtable, this->symbol(), probe, scale_field);
}

void FFSymbolFF::potential(PotentialTable &pottable, QVector<FFPtr> &forcefields,
                           const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                           double scale_potential) const
{
    FFPtr &ffield = forcefields[ffidx];
    
    if (not ffield->isA<FF3D>())
        throw SireFF::missing_derivative( QObject::tr(
            "The forcefield of type %1 does not inherit from FF3D so does "
            "not provide a force function.")
                .arg(ffield->what()), CODELOC );

    ffield.edit().asA<FF3D>().potential(pottable, this->symbol(), scale_potential);
}

void FFSymbolFF::potential(PotentialTable &pottable, const Probe &probe,
                           QVector<FFPtr> &forcefields,
                           const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                           double scale_potential) const
{
    FFPtr &ffield = forcefields[ffidx];
    
    if (not ffield->isA<FF3D>())
        throw SireFF::missing_derivative( QObject::tr(
            "The forcefield of type %1 does not inherit from FF3D so does "
            "not provide a force function.")
                .arg(ffield->what()), CODELOC );

    ffield.edit().asA<FF3D>().potential(pottable, this->symbol(), 
                                        probe, scale_potential);
}

///////////
/////////// Implementation of FFSymbolExpression
///////////

FFSymbolExpression::Component::Component()
{}

FFSymbolExpression::Component::Component(const Expression &scale_factor, 
                                         const Symbol &component, int idx)
                   : s(component), sclfac(scale_factor), ffidx(idx)
{
    Symbols dependents = scale_factor.symbols();
    
    if (not dependents.isEmpty())
    {
        deps = QVector<Symbol>( dependents.count() );
        Symbol *deps_array = deps.data();
        
        int i = 0;
        foreach (const Symbol &dependent, dependents)
        {
            deps_array[i] = dependent;
            ++i;
        }
    }
}

FFSymbolExpression::Component::Component(const Component &other)
                   : s(other.s), deps(other.deps), sclfac(other.sclfac), ffidx(other.ffidx)
{}

FFSymbolExpression::Component::~Component()
{}

int FFSymbolExpression::Component::ffIdx() const
{
    return ffidx;
}

int FFSymbolExpression::Component::nDependents() const
{
    return deps.count();
}

const QVector<Symbol>& FFSymbolExpression::Component::dependents() const
{
    return deps;
}

const Expression& FFSymbolExpression::Component::scalingExpression() const
{
    return sclfac;
}

double FFSymbolExpression::Component::scalingFactor(const Values &values) const
{
    return sclfac.evaluate(values);
}

const Symbol& FFSymbolExpression::Component::symbol() const
{
    return s;
}

FFSymbolExpression::FFSymbolExpression() : FFSymbol()
{}

FFSymbolExpression::FFSymbolExpression(const Symbol &symbol, 
                                       const Expression &expression)
                   : FFSymbol(symbol), ffexpression(expression)
{}

FFSymbolExpression::FFSymbolExpression(const FFSymbolExpression &other)
                   : FFSymbol(other), ffexpression(other.ffexpression),
                     components(other.components), sorted_components(other.sorted_components)
{}

FFSymbolExpression::~FFSymbolExpression()
{}

void FFSymbolExpression::load(QDataStream &ds)
{
    ds >> ffexpression;
    
    components.clear();
    
    FFSymbol::load(ds);
}

void FFSymbolExpression::save(QDataStream &ds) const
{
    ds << ffexpression;
    FFSymbol::save(ds);
}

const Expression& FFSymbolExpression::expression() const
{
    return ffexpression;
}

void FFSymbolExpression::expandInTermsOf(const QSet<Symbol> &ffsymbols,
                                         const QHash<Symbol,FFSymbolPtr> &ffmap)
{
    //get all of the symbols in which to expand this expression
    components.clear();
    sorted_components.clear();
    
    QSet<Symbol> symbols = ffexpression.symbols();
    symbols.intersect(ffsymbols);

    //now, we need to fully expand the forcefield expression. This is so that we
    //can evaluate all of the forcefields individually in one go, rather than
    //walking down a tree (and thus potentially evaluating the same forcefield
    //in two different walks) - we will expand into a copy of the expression
    Expression expanded = ffexpression;

    while (true)
    {
        bool fully_expanded = true;
        
        foreach (const Symbol &symbol, symbols)
        {
            if (ffmap.contains(symbol))
            {
                if (ffmap[symbol]->isA<FFSymbolExpression>())
                {
                    Expression e = ffmap[symbol]->asA<FFSymbolExpression>().expression();
                
                    expanded = expanded.substitute( symbol == e );
                    
                    symbols = expanded.symbols();
                    symbols.intersect(ffsymbols);
                    
                    fully_expanded = false;
                }
            }
        }
        
        if (fully_expanded)
        {
            break;
        }
    }
    
    Expression remainder = expanded;
    
    foreach (const Symbol &symbol, symbols)
    {
        QList<Factor> factors = expanded.expand(symbol);
        
        foreach (const Factor &factor, factors)
        {
            if (factor.power().isZero())
                continue;
        
            if (factor.power() != Expression(1))
                throw SireError::incompatible_error( QObject::tr(
                    "You cannot raise a forcefield energy component (%1) "
                    "to any power (%2) other than 1 as this is not "
                    "dimensionally correct.")
                        .arg(symbol.toString(), factor.power().toString()),
                            CODELOC );

            foreach (const Symbol &fac_symbol, factor.factor().symbols())
            {
                if (ffsymbols.contains(fac_symbol))
                    throw SireError::incompatible_error( QObject::tr(
                        "You cannot multiply or divide one forcefield energy "
                        "component by another (%1 by %2), as this is not "
                        "dimensionally correct.")
                            .arg(symbol.toString(), fac_symbol.toString()), CODELOC );
            }
            
            //now find out which forcefield (if any) this symbol needs to evaluate
            //the energy
            int ffidx = -1;
            
            if (ffmap.contains(symbol))
            {
                const FFSymbolPtr &ffptr = ffmap[symbol];
                
                if (ffptr->isA<FFSymbolFF>())
                {
                    const FFSymbolFF &ffsymff = ffptr->asA<FFSymbolFF>();
                    
                    ffidx = ffsymff.ffIdx();
                }
            }
            
            Component component = Component(factor.factor(), symbol, ffidx);
            components.append(component);
            
            remainder = (remainder - (symbol * factor.factor())).simplify();

            //if this is a forcefield component, then add this to the sorted
            //list, so that we can collect all symbols that belong to the same forcefield
            if (ffidx != -1)
            {
                bool found = false;
            
                for (int i=0; i<sorted_components.count(); ++i)
                {
                    if (sorted_components[i][0].ffIdx() == ffidx)
                    {
                        sorted_components[i].append(component);
                        found = true;
                        break;
                    }
                }
                
                if (not found)
                {
                    sorted_components.append( QVector<Component>() );
                    sorted_components.last().append(component);
                }
            }
        }
    }

    components.squeeze();

    int ncheck = 0;

    for (int i=0; i<sorted_components.count(); ++i)
    {
        sorted_components[i].squeeze();
        ncheck += sorted_components[i].count();
    }
    
    sorted_components.squeeze();

    if (ncheck != components.count())
    {
        qDebug() << "WARNING - DISAGREEMENT OF NUMBER OF COMPONENTS"
                 << ncheck << components.count();

        throw SireError::program_bug( QObject::tr(
                "Disagreement of number of components. %1 vs. %2")
                    .arg(ncheck).arg(components.count()), CODELOC );
    }

    //normally this should be zero, but there are some good cases (e.g. Nautilus)
    //when you want dimensionless constants to be added to the energy
    remainder = remainder.simplify();
}

Expression FFSymbolExpression::toExpression() const
{
    return ffexpression;
}

double FFSymbolExpression::value(const QHash<Symbol,FFSymbolPtr>&) const
{
    throw SireError::incompatible_error( QObject::tr(
        "There is no constant value associated with an energy expression. "
        "You cannot multiply one forcefield component (%1) by another in "
        "the forcefield expression %2.")
            .arg(symbol().toString(), ffexpression.toString()), CODELOC );
            
    return 0;
}

MolarEnergy FFSymbolExpression::energy(QVector<FFPtr> &forcefields,
                                       const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                                       double scale_energy) const
{
    if (scale_energy == 0)
        return MolarEnergy(0);
    
    int ncomponents = components.count();
    const Component *components_array = components.constData();
    
    //The energy expression takes the form;
    // E = Sum{i=1,n} ( {scale expression}_i * {forcefield energy component}_i )
    
    // where {scale_expression}_i is a constant expression that acts to scale (multiply)
    // the energy of the ith forcefield energy component, {forcefield energy component}_i
    
    // First, we have to calculate all of the values of the {scale expression}s,
    // knowing that these are constants and do not contain any forcefield energy
    // components (i.e. they are just values like 'lambda' or 'alpha')
    
    // we do this by putting all of the values into the below 'values' map
    // (this maps symbol to numeric value)
    Values values;

    // loop over all i=1,n components of the energy expression
    for (int i=0; i<ncomponents; ++i)
    {
        //this is the ith component
        const Component &component = components_array[i];
        
        //evaluate all of the dependent symbols for the scaling constant
        //for the energy in this ith component
        int ndeps = component.nDependents();
        const Symbol *deps_array = component.dependents().constData();
        
        for (int j=0; j<ndeps; ++j)
        {
            const Symbol &symbol = deps_array[j];

            if (not values.contains(symbol))
            {
                //all of the FFSymbols should be constants, so will raise
                //an exception here if they require a forcefield evaluation
                values.set( symbol, ffsymbols[symbol]->value(ffsymbols) );
            }
        }
    }

    //now that we have all of the scaling factors, we need to loop over
    //all dependent forcefields of this expression and will evaluate their
    //energies. To do this, we have to find all of the sub-forcefields,
    //as we must loop over forcefields, not components (else we risk
    //a race condition where the same forcefield is queried simultaneously
    //from two different threads)
    const int nforcefields = sorted_components.count();
    const QVector<Component> *sorted_components_array = sorted_components.constData();
    FFPtr *forcefields_array = forcefields.data();

    //loop over each forcefield in the expression and sum the total energy
    double total_nrg = tbb::parallel_reduce( tbb::blocked_range<int>(0,nforcefields), 0.0,
    [=](const tbb::blocked_range<int> &r, double nrg)->double
    {
        for (int i=r.begin(); i != r.end(); ++i)
        {
            const QVector<Component> &ffcomponents = sorted_components_array[i];
            const int ffidx = ffcomponents[0].ffIdx();

            FF &forcefield = forcefields_array[ffidx].edit();
        
            for (int j=0; j<ffcomponents.count(); ++j)
            {
                const Component &component = ffcomponents.constData()[j];

                const double scl = scale_energy * component.scalingFactor(values);

                if (scl != 0)
                {
                    nrg += forcefield.energy(component.symbol()) * scl;
                }
            }
        }
        
        return nrg;
    },
    []( double x, double y )->double
    {
        return x+y;
    });

    return MolarEnergy(total_nrg);
}

void FFSymbolExpression::energy(EnergyTable &energytable,
                               QVector<FFPtr> &forcefields,
                               const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                               double scale_energy) const
{
    int ncomponents = components.count();
    const Component *components_array = components.constData();
    
    Values values;
    
    for (int i=0; i<ncomponents; ++i)
    {
        const Component &component = components_array[i];
        
        //evaluate all of the dependent symbols...
        int ndeps = component.nDependents();
        const Symbol *deps_array = component.dependents().constData();
        
        for (int j=0; j<ndeps; ++j)
        {
            const Symbol &symbol = deps_array[j];

            if (not values.contains(symbol))
                values.set( symbol, ffsymbols[symbol]->value(ffsymbols) );
        }
        
        //now evaluate the scaling factor...
        double scale = scale_energy * component.scalingFactor(values);
        
        ffsymbols[component.symbol()]->energy(energytable, forcefields, 
					      ffsymbols, scale);
    }
}

void FFSymbolExpression::force(ForceTable &forcetable,
                               QVector<FFPtr> &forcefields,
                               const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                               double scale_force) const
{
    int ncomponents = components.count();
    const Component *components_array = components.constData();
    
    Values values;
    
    for (int i=0; i<ncomponents; ++i)
    {
        const Component &component = components_array[i];
        
        //evaluate all of the dependent symbols...
        int ndeps = component.nDependents();
        const Symbol *deps_array = component.dependents().constData();
        
        for (int j=0; j<ndeps; ++j)
        {
            const Symbol &symbol = deps_array[j];

            if (not values.contains(symbol))
                values.set( symbol, ffsymbols[symbol]->value(ffsymbols) );
        }
        
        //now evaluate the scaling factor...
        double scale = scale_force * component.scalingFactor(values);
        
        ffsymbols[component.symbol()]->force(forcetable, forcefields, 
                                             ffsymbols, scale);
    }
}

void FFSymbolExpression::field(FieldTable &fieldtable,
                               QVector<FFPtr> &forcefields,
                               const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                               double scale_field) const
{
    int ncomponents = components.count();
    const Component *components_array = components.constData();
    
    Values values;
    
    for (int i=0; i<ncomponents; ++i)
    {
        const Component &component = components_array[i];
        
        //evaluate all of the dependent symbols...
        int ndeps = component.nDependents();
        const Symbol *deps_array = component.dependents().constData();
        
        for (int j=0; j<ndeps; ++j)
        {
            const Symbol &symbol = deps_array[j];

            if (not values.contains(symbol))
                values.set( symbol, ffsymbols[symbol]->value(ffsymbols) );
        }
        
        //now evaluate the scaling factor...
        double scale = scale_field * component.scalingFactor(values);
        
        ffsymbols[component.symbol()]->field(fieldtable, forcefields, 
                                             ffsymbols, scale);
    }
}

void FFSymbolExpression::field(FieldTable &fieldtable,
                               const Probe &probe,
                               QVector<FFPtr> &forcefields,
                               const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                               double scale_field) const
{
    int ncomponents = components.count();
    const Component *components_array = components.constData();
    
    Values values;
    
    for (int i=0; i<ncomponents; ++i)
    {
        const Component &component = components_array[i];
        
        //evaluate all of the dependent symbols...
        int ndeps = component.nDependents();
        const Symbol *deps_array = component.dependents().constData();
        
        for (int j=0; j<ndeps; ++j)
        {
            const Symbol &symbol = deps_array[j];

            if (not values.contains(symbol))
                values.set( symbol, ffsymbols[symbol]->value(ffsymbols) );
        }
        
        //now evaluate the scaling factor...
        double scale = scale_field * component.scalingFactor(values);
        
        ffsymbols[component.symbol()]->field(fieldtable, probe, forcefields, 
                                             ffsymbols, scale);
    }
}

void FFSymbolExpression::potential(PotentialTable &pottable,
                                   QVector<FFPtr> &forcefields,
                                   const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                                   double scale_potential) const
{
    int ncomponents = components.count();
    const Component *components_array = components.constData();
    
    Values values;
    
    for (int i=0; i<ncomponents; ++i)
    {
        const Component &component = components_array[i];
        
        //evaluate all of the dependent symbols...
        int ndeps = component.nDependents();
        const Symbol *deps_array = component.dependents().constData();
        
        for (int j=0; j<ndeps; ++j)
        {
            const Symbol &symbol = deps_array[j];

            if (not values.contains(symbol))
                values.set( symbol, ffsymbols[symbol]->value(ffsymbols) );
        }
        
        //now evaluate the scaling factor...
        double scale = scale_potential * component.scalingFactor(values);
        
        ffsymbols[component.symbol()]->potential(pottable, forcefields, 
                                                 ffsymbols, scale);
    }
}

void FFSymbolExpression::potential(PotentialTable &pottable,
                                   const Probe &probe,
                                   QVector<FFPtr> &forcefields,
                                   const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                                   double scale_potential) const
{
    int ncomponents = components.count();
    const Component *components_array = components.constData();
    
    Values values;
    
    for (int i=0; i<ncomponents; ++i)
    {
        const Component &component = components_array[i];
        
        //evaluate all of the dependent symbols...
        int ndeps = component.nDependents();
        const Symbol *deps_array = component.dependents().constData();
        
        for (int j=0; j<ndeps; ++j)
        {
            const Symbol &symbol = deps_array[j];

            if (not values.contains(symbol))
                values.set( symbol, ffsymbols[symbol]->value(ffsymbols) );
        }
        
        //now evaluate the scaling factor...
        double scale = scale_potential * component.scalingFactor(values);
        
        ffsymbols[component.symbol()]->potential(pottable, probe, forcefields, 
                                                 ffsymbols, scale);
    }
}

///////////
/////////// Implementation of FFTotalExpression
///////////

FFTotalExpression::FFTotalExpression() : FFSymbol()
{}

FFTotalExpression::FFTotalExpression(const Symbol &symbol)
                  : FFSymbol(symbol)
{}

FFTotalExpression::FFTotalExpression(const FFTotalExpression &other)
                  : FFSymbol(other)
{}
  
FFTotalExpression::~FFTotalExpression()
{}

void FFTotalExpression::load(QDataStream &ds)
{
    FFSymbol::load(ds);
}

void FFTotalExpression::save(QDataStream &ds) const
{
    FFSymbol::save(ds);
}

Expression FFTotalExpression::toExpression() const
{
    return Expression(this->symbol());
}

double FFTotalExpression::value(const QHash<Symbol,FFSymbolPtr>&) const
{
    throw SireError::program_bug( QObject::tr(
        "An FFTotalExpression does not have a constant value and should never "
        "be used in a situation where its constant value must be determined..."), 
            CODELOC );
        
    return 0;
}

MolarEnergy FFTotalExpression::energy(QVector<FFPtr> &forcefields,
                                      const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                                      double scale_energy) const
{
    int nffields = forcefields.count();
    
    if (nffields <= 0)
    {
        return MolarEnergy(0);
    }
    
    FFPtr *ffields_array = forcefields.data();
    
    QVector<double> nrgs( nffields, 0 );
    double *nrgs_data = nrgs.data();
    
    tbb::parallel_for(0, nffields, 1, [=](int i)
    {
        nrgs_data[i] = ffields_array[i].edit().energy();
    });

    double nrg = 0;
    
    for (int i=0; i<nffields; ++i)
    {
        nrg += nrgs_data[i];
    }
    
    return MolarEnergy(nrg * scale_energy);
}

void FFTotalExpression::energy(EnergyTable &energytable,
                              QVector<FFPtr> &forcefields,
                              const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                              double scale_energy) const
{
    int nffields = forcefields.count();
    FFPtr *ffields_array = forcefields.data();
    
    for (int i=0; i<nffields; ++i)
    {
        FFPtr &ffield = ffields_array[i];
        
        if (ffield->isA<FF3D>())
            ffield.edit().asA<FF3D>().energy(energytable, scale_energy);
    }
}

void FFTotalExpression::force(ForceTable &forcetable,
                              QVector<FFPtr> &forcefields,
                              const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                              double scale_force) const
{
    int nffields = forcefields.count();
    FFPtr *ffields_array = forcefields.data();
    
    for (int i=0; i<nffields; ++i)
    {
        FFPtr &ffield = ffields_array[i];
        
        if (ffield->isA<FF3D>())
            ffield.edit().asA<FF3D>().force(forcetable, scale_force);
    }
}

void FFTotalExpression::field(FieldTable &fieldtable,
                              QVector<FFPtr> &forcefields,
                              const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                              double scale_field) const
{
    int nffields = forcefields.count();
    FFPtr *ffields_array = forcefields.data();
    
    for (int i=0; i<nffields; ++i)
    {
        FFPtr &ffield = ffields_array[i];
        
        if (ffield->isA<FF3D>())
            ffield.edit().asA<FF3D>().field(fieldtable, scale_field);
    }
}

void FFTotalExpression::field(FieldTable &fieldtable,
                              const Probe &probe,
                              QVector<FFPtr> &forcefields,
                              const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                              double scale_field) const
{
    int nffields = forcefields.count();
    FFPtr *ffields_array = forcefields.data();
    
    for (int i=0; i<nffields; ++i)
    {
        FFPtr &ffield = ffields_array[i];
        
        if (ffield->isA<FF3D>())
            ffield.edit().asA<FF3D>().field(fieldtable, probe, scale_field);
    }
}

void FFTotalExpression::potential(PotentialTable &pottable,
                                  QVector<FFPtr> &forcefields,
                                  const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                                  double scale_potential) const
{
    int nffields = forcefields.count();
    FFPtr *ffields_array = forcefields.data();
    
    for (int i=0; i<nffields; ++i)
    {
        FFPtr &ffield = ffields_array[i];
        
        if (ffield->isA<FF3D>())
            ffield.edit().asA<FF3D>().potential(pottable, scale_potential);
    }
}

void FFTotalExpression::potential(PotentialTable &pottable,
                                  const Probe &probe,
                                  QVector<FFPtr> &forcefields,
                                  const QHash<Symbol,FFSymbolPtr> &ffsymbols,
                                  double scale_potential) const
{
    int nffields = forcefields.count();
    FFPtr *ffields_array = forcefields.data();
    
    for (int i=0; i<nffields; ++i)
    {
        FFPtr &ffield = ffields_array[i];
        
        if (ffield->isA<FF3D>())
            ffield.edit().asA<FF3D>().potential(pottable, probe, scale_potential);
    }
}

///////////
/////////// Implementation of ForceFields
///////////

static const RegisterMetaType<ForceFields> r_ffields;

/** Serialise to a binary datastream */
QDataStream SIREFF_EXPORT &operator<<(QDataStream &ds, const ForceFields &ffields)
{
    writeHeader(ds, r_ffields, 2);
    
    SharedDataStream sds(ds);
    
    //write out all of the forcefields
    sds << ffields.ffields_by_idx
        << ffields.ffsymbols
        << ffields.additional_properties
        << ffields.property_aliases
        << ffields.combined_properties;
    
    return ds;
}
 
/** Extract from a binary datastream */
QDataStream SIREFF_EXPORT &operator>>(QDataStream &ds, ForceFields &ffields)
{
    VersionID v = readHeader(ds, r_ffields);
    
    if (v == 2)
    {
        SharedDataStream sds(ds);

        //read into a new object, as a lot can go wrong!
        ForceFields new_ffields;

        //read in the forcefields
        sds >> new_ffields.ffields_by_idx
            >> new_ffields.ffsymbols
            >> new_ffields.additional_properties
            >> new_ffields.property_aliases
            >> new_ffields.combined_properties;

        //rebuild the index
        new_ffields.rebuildIndex();
        
        ffields = new_ffields;
    }
    else if (v == 1)
    {
        SharedDataStream sds(ds);

        //read into a new object, as a lot can go wrong!
        ForceFields new_ffields;

        //read in the forcefields
        sds >> new_ffields.ffields_by_idx
            >> new_ffields.ffsymbols;

        //rebuild the index
        new_ffields.rebuildIndex();
        
        ffields = new_ffields;
    }
    else
        throw version_error(v, "1,2", r_ffields, CODELOC);

    return ds;
}

Symbol ForceFields::total_component("E_{total}");

/** Return the symbol representing the total energy component */
const Symbol& ForceFields::totalComponent()
{
    return total_component;
}

/** Constructor */
ForceFields::ForceFields() : ConcreteProperty<ForceFields,MolGroupsBase>()
{}

/** Internal function used to return the ith forcefield
    (this performs no bounds checking!) */
const FF& ForceFields::_pvt_forceField(int i) const
{
    return ffields_by_idx.constData()[i].read();
}

/** Internal function used to return the ith forcefield
    (this performs no bounds checking!) */
FF& ForceFields::_pvt_forceField(int i)
{
    return ffields_by_idx.data()[i].edit();
}

/** Internal function used to return the forcefield with name 'ffname'.
    This does not check to see if this forcefield exists */
const FF& ForceFields::_pvt_forceField(const FFName &ffname) const
{
    return this->_pvt_forceField( *(ffields_by_name.constFind(ffname)) );
}

/** Internal function used to return the forcefield with name 'ffname'.
    This does not check to see if this forcefield exists */
FF& ForceFields::_pvt_forceField(const FFName &ffname)
{
    return this->_pvt_forceField( *(ffields_by_name.constFind(ffname)) );
}

/** Internal function used to return the forcefield that contains the
    molecule group with number 'mgnum'. This does not check to see 
    if such a group exists */
const FF& ForceFields::_pvt_forceField(const MGNum &mgnum) const
{
    return this->_pvt_forceField( *(mgroups_by_num.constFind(mgnum)) );
}

/** Internal function used to return the forcefield that contains the
    molecule group with number 'mgnum'. This does not check to see 
    if such a group exists */
FF& ForceFields::_pvt_forceField(const MGNum &mgnum)
{
    return this->_pvt_forceField( *(mgroups_by_num.constFind(mgnum)) );
}

/** Reindex the moleculegroups and molecules */
void ForceFields::reindex()
{
    MolGroupsBase::clearIndex();

    int nffields = ffields_by_idx.count();
    const FFPtr *ffields_array = ffields_by_idx.constData();
    
    for (int i=0; i<nffields; ++i)
    {
        const FFPtr &ffield = ffields_array[i];
        
        for (int j=0; j<ffield->nGroups(); ++j)
        {
            MolGroupsBase::addToIndex( ffield->at(MGIdx(j)) );
        }
    }
}

/** Rebuild the index of forcefields and groups */
void ForceFields::rebuildIndex()
{
    ffields_by_name.clear();
    mgroups_by_num.clear();
    
    int nffields = ffields_by_idx.count();
    const FFPtr *ffields_array = ffields_by_idx.constData();
    
    for (int i=0; i<nffields; ++i)
    {
        const FFPtr &ffield = ffields_array[i];
        
        if (ffields_by_name.contains(ffield->name()))
        {
            const FFPtr &old_ffield = this->_pvt_forceField(ffield->name());
        
            throw SireFF::duplicate_forcefield( QObject::tr(
                "Cannot have two forcefields in the same set that both "
                "have the same name! (%1, %2 version %3 vs. %4 version %5)")
                    .arg(ffield->name())
                    .arg(ffield->UID().toString()).arg(ffield->version())
                    .arg(old_ffield->UID().toString()).arg(old_ffield->version()),
                            CODELOC );
        }
        
        ffields_by_name.insert( ffield->name(), i );
        
        foreach (MGNum mgnum, ffield->mgNums())
        {
            if (mgroups_by_num.contains(mgnum))
            {
                const FFPtr &old_ffield = this->_pvt_forceField(mgnum);
                
                throw SireMol::duplicate_group( QObject::tr(
                    "Cannot have two different forcefields containing the same "
                    "molecule group - %1 (%2 version %3 %4 vs. "
                    "%5 version %6 %7)")
                        .arg(mgnum)
                        .arg(ffield->name()).arg(ffield->version())
                        .arg(ffield->UID().toString())
                        .arg(old_ffield->name()).arg(old_ffield->version())
                        .arg(old_ffield->UID().toString()), CODELOC );
            }
        
            mgroups_by_num.insert( mgnum, ffield->name() );
        }
    }
    
    ffields_by_name.squeeze();
    mgroups_by_num.squeeze();

    //now rebuild the index of molecules and molecule groups
    MolGroupsBase::clearIndex();
    
    for (int i=0; i<nffields; ++i)
    {
        const FFPtr &ffield = ffields_array[i];
        
        for (int j=0; j<ffield->nGroups(); ++j)
        {
            MolGroupsBase::addToIndex( ffield->at(MGIdx(j)) );
        }
    }
    
    //now rebuild the index of forcefield expressions and symbols
    QHash<Symbol,FFSymbolPtr> new_symbols;

    //first copy in all of the symbols representing all of the forcefield
    //components
    QSet<Symbol> all_ff_symbols;
    
    for (FFIdx i(0); i<nffields; ++i)
    {
        Symbols symbols = ffields_array[i]->components().symbols();

        foreach (const Symbol &symbol, symbols)
        {
            if (new_symbols.contains(symbol))
                throw SireError::program_bug( QObject::tr(
                    "It should not be possible for two forcefields to have "
                    "the same component symbol... (%1)")
                        .arg(symbol.toString()), CODELOC );
        
            new_symbols.insert( symbol, FFSymbolPtr(new FFSymbolFF(i, symbol)) );
            all_ff_symbols.insert(symbol);
        }
    }

    //copy in the non-forcefield symbols from the old array
    for (QHash<Symbol,FFSymbolPtr>::const_iterator it = ffsymbols.constBegin();
         it != ffsymbols.constEnd();
         ++it)
    {
        if ( not (it.value()->isA<FFSymbolFF>() or 
                  it.value()->isA<FFTotalExpression>()) )
        {
            if (new_symbols.contains(it.key()))
                throw SireFF::duplicate_component( QObject::tr(
                    "You cannot use the symbol %1 to represent a forcefield "
                    "expression or parameter as it already in use to represent "
                    "an energy component of a forcefield.")
                        .arg(it.key().toString()), CODELOC );
            
            new_symbols.insert(it.key(), it.value());
        }
    }
    
    //if there isn't a total energy component, then add the default one
    if (not new_symbols.contains( this->totalComponent() ))
    {
        new_symbols.insert(this->totalComponent(),
                           FFSymbolPtr(new FFTotalExpression()));
    }

    for (QHash<Symbol,FFSymbolPtr>::iterator it = new_symbols.begin();
         it != new_symbols.end();
         ++it)
    {
        if (it.value()->isA<FFSymbolExpression>() or
            it.value()->isA<FFTotalExpression>())
        {
            all_ff_symbols.insert(it.key());
        }
    }
    
    //now process each forcefield expression...
    for (QHash<Symbol,FFSymbolPtr>::iterator it = new_symbols.begin();
         it != new_symbols.end();
         ++it)
    {
        if (it.value()->isA<FFSymbolExpression>())
        {
            it.value()->asA<FFSymbolExpression>().expandInTermsOf(all_ff_symbols,
                                                                  new_symbols);
        }
    }
    
    //now loop through the constant expressions and make sure
    //that they really are constant!
    for (QHash<Symbol,FFSymbolPtr>::iterator it = new_symbols.begin();
         it != new_symbols.end();
         ++it)
    {
        if (it.value()->isA<FFConstantExpression>())
        {
            it.value()->asA<FFConstantExpression>().assertNotDepends(all_ff_symbols);
        }
    }
    
    ffsymbols = new_symbols;
}

/** Construct a group that holds just a single forcefield */
ForceFields::ForceFields(const FF& forcefield)
            : ConcreteProperty<ForceFields,MolGroupsBase>()
{
    ffields_by_idx.append(forcefield);
    this->rebuildIndex();
}

/** Construct a group that holds lots of forcefields */
ForceFields::ForceFields(const QList<FFPtr> &forcefields)
            : ConcreteProperty<ForceFields,MolGroupsBase>()
{
    ffields_by_idx = forcefields.toVector();
    
    int nffields = ffields_by_idx.count();
    
    if (nffields > 1)
    {
        Molecules mols;
        
        FFPtr *ffields_array = ffields_by_idx.data();
        
        for (int i=1; i<nffields; ++i)
        {
            mols += ffields_array[i-1]->molecules();
            ffields_array[i].edit().update(mols);
        }
    }
    
    this->rebuildIndex();
}

/** Construct a group that holds lots of forcefields */
ForceFields::ForceFields(const QVector<FFPtr> &forcefields)
            : ConcreteProperty<ForceFields,MolGroupsBase>(),
              ffields_by_idx(forcefields)
{
    this->rebuildIndex();
}

/** Copy constructor */
ForceFields::ForceFields(const ForceFields &other)
            : ConcreteProperty<ForceFields,MolGroupsBase>(other),
              ffields_by_idx(other.ffields_by_idx),
              ffields_by_name(other.ffields_by_name),
              mgroups_by_num(other.mgroups_by_num),
              ffsymbols(other.ffsymbols),
              additional_properties(other.additional_properties),
              property_aliases(other.property_aliases),
              combined_properties(other.combined_properties)
{}

/** Destructor */
ForceFields::~ForceFields()
{}

/** Copy assignment operator */
ForceFields& ForceFields::operator=(const ForceFields &other)
{
    if (this != &other)
    {
        ffields_by_idx = other.ffields_by_idx;
        ffields_by_name = other.ffields_by_name;
        mgroups_by_num = other.mgroups_by_num;
        ffsymbols = other.ffsymbols;
        additional_properties = other.additional_properties;
        property_aliases = other.property_aliases;
        combined_properties = other.combined_properties;
        
        MolGroupsBase::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool ForceFields::operator==(const ForceFields &other) const
{
    return this == &other or
           ( ffields_by_idx == other.ffields_by_idx and
             additional_properties == other.additional_properties and
             property_aliases == other.property_aliases and
             combined_properties == other.combined_properties );
}

/** Comparison operator */
bool ForceFields::operator!=(const ForceFields &other) const
{
    return this != &other and
           ( ffields_by_idx != other.ffields_by_idx or
             additional_properties != other.additional_properties or
             property_aliases != other.property_aliases or
             combined_properties != other.combined_properties );
}

/** Internal function used to return the group with number 'mgnum' */
const MoleculeGroup& ForceFields::getGroup(MGNum mgnum) const
{
    if (not mgroups_by_num.contains(mgnum))
        throw SireMol::missing_group( QObject::tr(
            "None of the forcefields in this set contain a molecule group "
            "with number %1. Available groups are %2.")
                .arg(mgnum).arg(Sire::toString(mgroups_by_num.keys())), 
                    CODELOC );

    return this->_pvt_forceField(mgnum).group(mgnum);
}

/** Internal function used to get the pointers to lots of groups */
void ForceFields::getGroups(const QList<MGNum> &mgnums,
                            QVarLengthArray<const MoleculeGroup*,10> &groups) const
{
    groups.clear();
    
    foreach (MGNum mgnum, mgnums)
    {
        groups.append( &(this->getGroup(mgnum)) );
    }
}

/** Internal function used to get pointers to all of the groups
    in all of the forcefields of this set */
QHash<MGNum,const MoleculeGroup*> ForceFields::getGroups() const
{
    QHash<MGNum,const MoleculeGroup*> groups;
    
    for (QHash<MGNum,FFName>::const_iterator it = mgroups_by_num.constBegin();
         it != mgroups_by_num.constEnd();
         ++it)
    {
        groups.insert( it.key(), &(this->_pvt_forceField(it.value()).group(it.key())) );
    }
    
    return groups;
}

/** Return the forcefield with name 'ffname'

    \throw SireFF::missing_forcefield
*/
const FF& ForceFields::operator[](const FFName &ffname) const
{
    if (not ffields_by_name.contains(ffname))
        throw SireFF::missing_forcefield( QObject::tr(
            "There is no forcefield called \"%1\" in this set. Available "
            "forcefields are called %2.")
                .arg(ffname).arg( Sire::toString(ffields_by_name.keys()) ),
                    CODELOC );
                    
    return this->_pvt_forceField(ffname);
}

/** Return the forcefield at index 'ffidx'

    \throw SireError::invalid_index
*/
const FF& ForceFields::operator[](const FFIdx &ffidx) const
{
    return this->_pvt_forceField( ffidx.map( ffields_by_idx.count() ) );
}

/** Return the forcefield with ID 'ffid'

    \throw SireFF::missing_forcefield
    \throw SireFF::duplicate_forcefield
    \throw SireError::invalid_index
*/
const FF& ForceFields::operator[](const FFID &ffid) const
{
    return this->_pvt_forceField( this->ffIdx(ffid) );
}

/** Return the forcefield with name 'ffname'

    \throw SireFF::missing_forcefield
*/
const FF& ForceFields::at(const FFName &ffname) const
{
    return this->operator[](ffname);
}

/** Return the forcefield at index 'ffidx'

    \throw SireError::invalid_index
*/
const FF& ForceFields::at(const FFIdx &ffidx) const
{
    return this->operator[](ffidx);
}

/** Return the forcefield with ID 'ffid'

    \throw SireFF::missing_forcefield
    \throw SireFF::duplicate_forcefield
    \throw SireError::invalid_index
*/
const FF& ForceFields::at(const FFID &ffid) const
{
    return this->operator[](ffid);
}

/** Return the forcefield with name 'ffname'

    \throw SireFF::missing_forcefield
*/
const FF& ForceFields::forceField(const FFName &ffname) const
{
    return this->operator[](ffname);
}

/** Return the forcefield at index 'ffidx'

    \throw SireError::invalid_index
*/
const FF& ForceFields::forceField(const FFIdx &ffidx) const
{
    return this->operator[](ffidx);
}

/** Return the forcefield with ID 'ffid'

    \throw SireFF::missing_forcefield
    \throw SireFF::duplicate_forcefield
    \throw SireError::invalid_index
*/
const FF& ForceFields::forceField(const FFID &ffid) const
{
    return this->operator[](ffid);
}

/** Return the forcefield that contains the molecule group with
    number 'mgnum'
    
    \throw SireMol::missing_group
*/
const FF& ForceFields::forceField(const MGNum &mgnum) const
{
    if (not mgroups_by_num.contains(mgnum))
        throw SireMol::missing_group( QObject::tr(
            "None of the forcefields in this set contain a molecule group "
            "with number %1. Available groups are %2.")
                .arg(mgnum).arg(Sire::toString(mgroups_by_num.keys())), 
                    CODELOC );

    return this->_pvt_forceField(mgnum);
}

/** Return the number of forcefields in this set */
int ForceFields::nForceFields() const
{
    return ffields_by_idx.count();
}
    
/** Return the index of the forcefield with name 'ffname'

    \throw SireFF::missing_forcefield 
    \throw SireFF::duplicate_forcefield
*/
FFIdx ForceFields::ffIdx(const FFName &ffname) const
{
    QList<FFIdx> ffidxs = ffname.map(*this);
    
    if (ffidxs.count() > 1)
        throw SireFF::duplicate_forcefield( QObject::tr(
            "There is more than one forcefield that matches the name "
            "\"%1\".").arg(ffname), CODELOC );

    return ffidxs.first();
}

/** Simple function that allows a shortcut for ffIdx(FFIdx)

    \throw SireError::invalid_index
*/
FFIdx ForceFields::ffIdx(const FFIdx &ffidx) const
{
    return FFIdx( ffidx.map(this->nForceFields()) );
}

/** Return the FFIdx of the forcefield that matches the ID 'ffid'

    \throw SireFF::missing_forcefield
    \throw SireFF::duplicate_forcefield
    \throw SireError::invalid_index
*/
FFIdx ForceFields::ffIdx(const FFID &ffid) const
{
    QList<FFIdx> ffidxs = ffid.map(*this);
    
    if (ffidxs.count() > 1)
        throw SireFF::duplicate_forcefield( QObject::tr(
            "More than one forcefield in this set matches the ID %1. "
            "Matching forcefield have indicies %2.")
                .arg(ffid.toString(), Sire::toString(ffidxs)),
                    CODELOC );
                    
    return ffidxs.at(0);
}
    
QList<FFIdx> ForceFields::map(const FFID &ffid) const
{
    return ffid.map(*this);
}

QList<FFIdx> ForceFields::map(const FFIdx &ffidx) const
{
    QList<FFIdx> ffidxs;
    
    ffidxs.append( this->ffIdx(ffidx) );
    
    return ffidxs;
}

QList<FFIdx> ForceFields::map(const FFName &ffname) const
{
    QList<FFIdx> ffidxs;

    if (ffname.isCaseSensitive())
    {
        QHash<QString,int>::const_iterator it = ffields_by_name.constFind(ffname);
        
        if (it != ffields_by_name.constEnd())
            ffidxs.append( FFIdx(it.value()) );
    }
    else
    {
        QString lower_name = QString(ffname).toLower();
        
        for (QHash<QString,int>::const_iterator it = ffields_by_name.constBegin();
             it != ffields_by_name.constEnd();
             ++it)
        {
            if (it.key().toLower() == lower_name)
                ffidxs.append( FFIdx(it.value()) );
        }
    }

    if (ffidxs.isEmpty())
        throw SireFF::missing_forcefield( QObject::tr(
            "There are no forcefields in this set called %1. "
            "Available forcefields are called %2.")
                .arg(ffname).arg( Sire::toString(ffields_by_name.keys()) ),
                    CODELOC );

    return ffidxs;
}

/** Simple function that short cuts ffName(FFName) 

    \throw SireFF::missing_forcefield
*/
const FFName& ForceFields::ffName(const FFName &ffname) const
{
    QList<FFIdx> ffidxs = this->map(ffname);

    if (ffidxs.count() > 1)
        throw SireFF::duplicate_forcefield( QObject::tr(
            "More than one forcefield in this set has the name %1. "
            "Matching forcefield have indicies %2.")
                .arg(ffname.toString(), Sire::toString(ffidxs)),
                    CODELOC );

    return ffname;
}

/** Return the name of the forcefield at index 'ffidx'

    \throw SireError::invalid_index
*/
const FFName& ForceFields::ffName(const FFIdx &ffidx) const
{
    return this->_pvt_forceField( ffidx.map(this->nForceFields()) ).name();
}

/** Return the name of the forcefield that matches the ID 'ffid'

    \throw SireFF::missing_forcefield
    \throw SireFF::duplicate_forcefield
    \throw SireError::invalid_index
*/
const FFName& ForceFields::ffName(const FFID &ffid) const
{
    QList<FFIdx> ffidxs = ffid.map(*this);

    if (ffidxs.count() > 1)
        throw SireFF::duplicate_forcefield( QObject::tr(
            "More than one forcefield in this set matches the ID %1. "
            "Matching forcefield have indicies %2.")
                .arg(ffid.toString(), Sire::toString(ffidxs)),
                    CODELOC );

    return this->_pvt_forceField( ffidxs.at(0).map(this->nForceFields()) ).name();
}

/** Return a string representation of this set */
QString ForceFields::toString() const
{
    return QObject::tr("FFPtr( nForceFields() == %1 )").arg(this->nForceFields());
}
    
/** Return the energy associated with the symbol 'component'. This component 
    may either be a component of one of the constituent forcefields,
    or it may represent a function of the forcefield components
    
    \throw SireFF::missing_component
*/
SireUnits::Dimension::MolarEnergy ForceFields::energy(const Symbol &component)
{
    FFSymbolPtr comp = ffsymbols.value(component);

    if (comp.get() == 0)
        throw SireFF::missing_component( QObject::tr(   
            "There is no component of the energy represented by the "
            "symbol %1. Available components are %2.")
                .arg(component.toString(), Sire::toString(energySymbols())),
                    CODELOC );

    if (comp->isConstant())
        throw SireFF::missing_component( QObject::tr(
                "The component %1 is a constant component (it is not an energy "
                "component). Available energy components are %2.")
                    .arg(component.toString(),
                         Sire::toString(energySymbols())), CODELOC );

    return comp->energy(ffields_by_idx, ffsymbols);
}

/** Return the energy or constant component associated with the symbol 'symbol'

    \throw SireFF::missing_component
*/
double ForceFields::componentValue(const Symbol &component)
{
    FFSymbolPtr comp = ffsymbols.value(component);

    if (comp.get() == 0)
        throw SireFF::missing_component( QObject::tr(   
            "There is no available component represented by the "
            "symbol %1. Available components are %2.")
                .arg(component.toString(), Sire::toString(componentSymbols())),
                    CODELOC );

    if (comp->isConstant())
        return comp->value(ffsymbols);
    else
        return comp->energy(ffields_by_idx, ffsymbols).value();
    
}

/** Return the values of the energy or constant components whose
    symbols are in 'symbols'
    
    \throw SireFF::missing_component
*/
Values ForceFields::componentValues(const QSet<Symbol> &components)
{
    Values vals;
    
    foreach (const Symbol &component, components)
    {
        vals.set( component, this->componentValue(component) );
    }
    
    return vals;
}

/** Return the values of all energy and constant components */
Values ForceFields::componentValues()
{
    Values vals;
    vals.reserve(ffsymbols.count());
    
    QHashIterator<Symbol,FFSymbolPtr> it( ffsymbols );
    
    while (it.hasNext())
    {
        it.next();
        
        if (it.value()->isConstant())
            vals.set(it.key(), it.value()->value(ffsymbols));
        else
            vals.set(it.key(), it.value()->energy(ffields_by_idx, ffsymbols).value());
    }
    
    return vals;
}

/** Return the energy of this set of forcefields. This uses the supplied
    total energy function to calculate the energy, if one exists,
    or it just calculates the sum of the total energies of all of the
    contained forcefields */
SireUnits::Dimension::MolarEnergy ForceFields::energy()
{
    return this->energy( this->totalComponent() );
}

/** Return the energies of all of the energy components of all of the forcefields,
    constants and expressions */
Values ForceFields::energies()
{
    Values vals;
    
    for (QHash<Symbol,FFSymbolPtr>::const_iterator it = ffsymbols.constBegin();
         it != ffsymbols.constEnd();
         ++it)
    {
        if ( (*it)->isEnergy() )
            vals.set(it.key(), this->energy(it.key()).value());
    }

    return vals;
}

/** Return the energies of all of the energy components whose symbols are 
    listed in 'components'
    
    \throw SireFF::missing_component
*/
Values ForceFields::energies(const QSet<Symbol> &components)
{
    Values vals;
    vals.reserve(components.count());
    
    foreach (const Symbol &component, components)
    {
        vals.set(component, this->energy(component).value());
    }
    
    return vals;
}

/** Return whether or not the forcefield component 'component'
    is an energy component
    
    \throw SireFF::missing_component
*/
bool ForceFields::isEnergyComponent(const Symbol &component) const
{
    FFSymbolPtr comp = ffsymbols.value(component);

    if (comp.get() == 0)
        throw SireFF::missing_component( QObject::tr(   
            "There is no component of the energy represented by the "
            "symbol %1. Available components are %2.")
                .arg(component.toString(), Sire::toString(energySymbols())),
                    CODELOC );

    return comp->isEnergy();
}

/** Return whether or not the forcefields have an energy component
    with symbol 'component' */
bool ForceFields::hasEnergyComponent(const Symbol &component) const
{
    FFSymbolPtr comp = ffsymbols.value(component);

    if (comp.get() != 0)
        return comp->isEnergy();
    else
        return false;
}

/** Set the component represented by the symbol 'symbol' equal to the expression
    contained in 'expression'. This replaces any existing constant or
    energy component with this expression. 
    
    Note that this expression must only involve terms that are linear in 
    forcefield components, and there may be no products of forcefield
    components (i.e. each term of the expression must have dimensions
    of energy)
    
    Note that an exception will be raised if you try to replace a 
    component that exists in one of the constituent forcefields
    
    \throw SireFF::duplicate_component
*/
void ForceFields::setEnergyComponent(const Symbol &symbol, const Expression &expression)
{
    ForceFields old_state( *this );
    
    try
    {
        ffsymbols.insert( symbol, 
                          FFSymbolPtr(new FFSymbolExpression(symbol, expression)) );
        this->rebuildIndex();
    }
    catch(...)
    {
        this->operator=(old_state);
        throw;
    }
}

/** Return all of the symbols that represent energy components */
QSet<Symbol> ForceFields::energySymbols() const
{
    QSet<Symbol> syms;
    
    for (QHash<Symbol,FFSymbolPtr>::const_iterator it = ffsymbols.constBegin();
         it != ffsymbols.constEnd();
         ++it)
    {
        if ( (*it)->isEnergy() )
            syms.insert( it.key() );
    }

    return syms;
}

/** Synonym for ForceFields::energies() */
Values ForceFields::energyComponents()
{
    return this->energies();
}

/** Return the energy expression for the symbol 'component'

    \throw SireFF::missing_component
*/
Expression ForceFields::energyExpression(const Symbol &component) const
{
    FFSymbolPtr comp = ffsymbols.value(component);

    if (comp.get() == 0)
        throw SireFF::missing_component( QObject::tr(   
            "There is no component of the energy represented by the "
            "symbol %1. Available components are %2.")
                .arg(component.toString(), Sire::toString(energySymbols())),
                    CODELOC );

    if (comp->isConstant())
        throw SireFF::missing_component( QObject::tr(
                "The component %1 is a constant component (it is not an energy "
                "component). Available energy components are %2.")
                    .arg(component.toString(),
                         Sire::toString(energySymbols())), CODELOC );

    return comp->toExpression();
}

/** Return the energy expressions for the components whose
    symbols are in 'symbols'
    
    \throw SireFF::missing_component
*/
QHash<Symbol,Expression> ForceFields::energyExpressions(
                                        const QSet<Symbol> &symbols) const
{
    QHash<Symbol,Expression> exps;
    exps.reserve( symbols.count() );
    
    foreach (const Symbol &symbol, symbols)
    {
        exps.insert( symbol, this->energyExpression(symbol) );
    }
    
    return exps;
}

/** Return all of the energy expressions in this forcefield */
QHash<Symbol,Expression> ForceFields::energyExpressions() const
{
    QHash<Symbol,Expression> exps;
    
    for (QHash<Symbol,FFSymbolPtr>::const_iterator it = ffsymbols.constBegin();
         it != ffsymbols.constEnd();
         ++it)
    {
        if ( (*it)->isEnergy() )
            exps.insert( it.key(), (*it)->toExpression() );
    }

    return exps;
}
                        
/** Return the value associated with the constant component 'component'

    \throw SireFF::missing_component
*/
double ForceFields::constant(const Symbol &component) const
{
    if (not ffsymbols.contains(component))
    {
        throw SireFF::missing_component( QObject::tr(
            "There is no constant represented by the symbol %1. "
            "Available constants are %2.")
                .arg(component.toString(), 
                        Sire::toString(this->constantSymbols())),
                    CODELOC );
    }
    
    FFSymbolPtr val = ffsymbols.value(component);
    
    if (not val->isConstant())
    {
        throw SireFF::missing_component( QObject::tr(
            "There is no constant represented by the symbol %1. There is an "
            "energy component with this symbol, but this is not constant! "
            "Available constants are %2.")
                .arg(component.toString(), 
                     Sire::toString(this->constantSymbols())),
                        CODELOC );
    }
    
    return val->value(ffsymbols);
}

/** Return the values of all constant components in the forcefields */
Values ForceFields::constants() const
{
    Values constant_values;

    for (QHash<Symbol,FFSymbolPtr>::const_iterator it = ffsymbols.constBegin();
         it != ffsymbols.constEnd();
         ++it)
    {
        if ((*it)->isConstant())
        {
            constant_values.set(it.key(), (*it)->value(ffsymbols));
        }
    }
    
    return constant_values;
}

/** Return the constant values associated with the constant components
    in 'components'

    \throw SireFF::missing_component
*/
Values ForceFields::constants(const QSet<Symbol> &components) const
{
    Values constant_values;
    
    for (QSet<Symbol>::const_iterator it = components.constBegin();
         it != components.constEnd();
         ++it)
    {
        constant_values.set( *it, this->constant(*it) );
    }
    
    return constant_values;
}

/** Return whether or not the forcefield component 'component'
    is a constant component
    
    \throw SireFF::missing_component
*/
bool ForceFields::isConstantComponent(const Symbol &component) const
{
    if (not ffsymbols.contains(component))
        throw SireFF::missing_component( QObject::tr(
                "There is no forcefield component %1 (constant or otherwise!). "
                "Available constant components are %2.")
                    .arg(component.toString(),
                         Sire::toString(constantSymbols())), CODELOC);
                         
    return this->hasConstantComponent(component);
}

/** Return whether or not there is a constant component in the 
    forcefield expressions with symbol 'component' */
bool ForceFields::hasConstantComponent(const Symbol &component) const
{
    if (ffsymbols.contains(component))
    {
        return ffsymbols.value(component)->isConstant();
    }
    else
        return false;
}

/** Set the constant component represented by the symbol 'symbol' equal to the 
    value 'value'. This replaces any existing constant or energy component with 
    this value. Note that an exception will be raised if you try to replace a component
    that exists in one of the constituent forcefields.
    
    \throw SireFF::duplicate_component
*/
void ForceFields::setConstantComponent(const Symbol &symbol, double value)
{
    if (ffsymbols.contains(symbol))
    {
        //if we are just changing the value then we can short-cut the process
        const FFSymbolPtr &ffsym = ffsymbols.constFind(symbol).value();
        
        if (ffsym->isA<FFConstantValue>())
        {
            if (ffsym->asA<FFConstantValue>().value() == value)
                //no need to change anything
                return;
            else
            {
                ffsymbols[symbol] = FFSymbolPtr( new FFConstantValue(symbol,value) );
                return;
            }
        }
    }

    ForceFields old_state( *this );

    try
    {
        ffsymbols.insert( symbol, 
                          FFSymbolPtr(new FFConstantValue(symbol, value)) );
        this->rebuildIndex();
    }
    catch(...)
    {
        this->operator=(old_state);
        throw;
    }
}

/** Set the constant component represented by the symbol 'symbol' equal to the 
    expression 'expression'. This replaces any existing constant or energy component with 
    this value. 
    
    Note that this expression must only involve constants, or other constant
    components. A constant expression may not depend on a forcefield energy
    
    Note that an exception will be raised if you try to replace a component
    that exists in one of the constituent forcefields.
    
    \throw SireFF::duplicate_component
*/
void ForceFields::setConstantComponent(const Symbol &symbol,
                                       const Expression &expression)
{
    if (expression.isConstant())
    {
        this->setConstantComponent(symbol, expression.evaluate(Values()));
        return;
    }

    ForceFields old_state( *this );

    try
    {
        ffsymbols.insert( symbol, 
                          FFSymbolPtr(new FFConstantExpression(symbol, expression)) );
        this->rebuildIndex();
    }
    catch(...)
    {
        this->operator=(old_state);
        throw;
    }
}

/** Return the symbols representing the constant forcefield components */
QSet<Symbol> ForceFields::constantSymbols() const
{
    QSet<Symbol> constant_symbols;

    for (QHash<Symbol,FFSymbolPtr>::const_iterator it = ffsymbols.constBegin();
         it != ffsymbols.constEnd();
         ++it)
    {
        if ((*it)->isConstant())
        {
            constant_symbols.insert(it.key());
        }
    }
    
    return constant_symbols;
}

/** Synonym for ForceFields::constants() */
Values ForceFields::constantComponents() const
{
    return this->constants();
}

/** Return the expression for the constant component 'symbol'

    \throw SireFF::missing_component
*/
Expression ForceFields::constantExpression(const Symbol &component) const
{
    if (not ffsymbols.contains(component))
        throw SireFF::missing_component( QObject::tr(
                "There is no forcefield component %1 (constant or otherwise!). "
                "Available constant components are %2.")
                    .arg(component.toString(),
                         Sire::toString(constantSymbols())), CODELOC);

    FFSymbolPtr val = ffsymbols.value(component);
    
    if (not val->isConstant())
        throw SireFF::missing_component( QObject::tr(
                "The component %1 is not a constant component. Available constant "
                "components are %2.")
                    .arg(component.toString(),
                         Sire::toString(constantSymbols())), CODELOC );

    return val->toExpression();
}

/** Return the expressions for the constant components in 'symbols'

    \throw SireFF::missing_component
*/
QHash<Symbol,Expression> ForceFields::constantExpressions(
                                            const QSet<Symbol> &symbols) const
{
    QHash<Symbol,Expression> exps;
    exps.reserve( symbols.count() );
    
    foreach (const Symbol &symbol, symbols )
    {
        exps.insert( symbol, this->constantExpression(symbol) );
    }
    
    return exps;
}

/** Return all of the constant expressions in the forcefields */
QHash<Symbol,Expression> ForceFields::constantExpressions() const
{
    QHash<Symbol,Expression> exps;

    for (QHash<Symbol,FFSymbolPtr>::const_iterator it = ffsymbols.constBegin();
         it != ffsymbols.constEnd();
         ++it)
    {
        if ((*it)->isConstant())
        {
            exps.insert(it.key(), it.value()->toExpression());
        }
    }
    
    return exps;
}

/** Set the component represented by the symbol 'symbol' equal to the 
    value 'value'. This replaces any existing component with this value.
    Note that an exception will be raised if you try to replace a component
    that exists in one of the constituent forcefields.
    
    This is a convenient shorthand for ForceFields::setConstantComponent(symbol, value)
    
    \throw SireFF::duplicate_component
*/
void ForceFields::setComponent(const Symbol &symbol, double value)
{
    this->setConstantComponent(symbol, value);
}

/** Set the component represented by the symbol 'symbol' equal to the expression
    contained in 'expression'. This replaces any existing constant or
    energy component with this expression. 
    
    Note that this expression must only involve terms that are linear in 
    forcefield components, and there may be no products of forcefield
    components (i.e. each term of the expression must have dimensions
    of energy)
    
    Note that an exception will be raised if you try to replace a 
    component that exists in one of the constituent forcefields
    
    This is a convenient short-hand for 
    ForceFields::setEnergyComponent(symbol,expression)
    
    \throw SireFF::duplicate_component
*/
void ForceFields::setComponent(const Symbol &symbol, const Expression &expression)
{
    if (expression.isConstant())
        this->setConstantComponent(symbol, expression);
    else
        this->setEnergyComponent(symbol, expression);
}

/** Return the symbols representing all of the constant and energy components */
QSet<Symbol> ForceFields::componentSymbols() const
{
    return ffsymbols.keys().toSet();
}

/** Return whether or not there is a constant or energy component with symbol 'symbol' */
bool ForceFields::hasComponent(const Symbol &symbol) const
{
    return ffsymbols.contains(symbol);
}

/** Return the expression for the constant or energy component 
    represented by 'symbol'
    
    \throw SireFF::missing_component
*/
Expression ForceFields::componentExpression(const Symbol &symbol) const
{
    FFSymbolPtr sym = ffsymbols.value(symbol);
    
    if (sym.get() == 0)
        throw SireFF::missing_component( QObject::tr(
                "There is no energy or constant component %1. Available "
                "components are %2.")
                    .arg(symbol.toString(),
                         Sire::toString(ffsymbols.keys())), CODELOC );
                         
    return sym->toExpression();
}

/** Return all of the expressions for the constant or energy
    components whose symbols are in 'symbols'
    
    \throw SireFF::missing_component
*/
QHash<Symbol,Expression> ForceFields::componentExpressions(
                                        const QSet<Symbol> &symbols) const
{
    QHash<Symbol,Expression> exps;
    exps.reserve(symbols.count());
    
    foreach (const Symbol &symbol, symbols)
    {
        exps.insert( symbol, this->componentExpression(symbol) );
    }
    
    return exps;
}

/** Return all of the constant and energy expressions attached
    to these forcefields */
QHash<Symbol,Expression> ForceFields::componentExpressions() const
{
    QHash<Symbol,Expression> exps;
    
    for (QHash<Symbol,FFSymbolPtr>::const_iterator it = ffsymbols.constBegin();
         it != ffsymbols.constEnd();
         ++it)
    {
        exps.insert(it.key(), it.value()->toExpression());
    }
    
    return exps;
}

void ForceFields::energy(EnergyTable &energytable, const Symbol &component,
			 double scale_energy)
{
    FFSymbolPtr comp = ffsymbols.value(component);

    if (comp.get() == 0)
        throw SireFF::missing_component( QObject::tr(   
            "There is no component of the energy represented by the "
            "symbol %1. Available components are %2.")
                .arg(component.toString(), Sire::toString(energySymbols())),
                    CODELOC );

    if (comp->isConstant())
        throw SireFF::missing_component( QObject::tr(
                "The component %1 is a constant component (it is not an energy "
                "component). Available energy components are %2.")
                    .arg(component.toString(),
                         Sire::toString(energySymbols())), CODELOC );
    comp->energy(energytable, ffields_by_idx, ffsymbols, scale_energy);
}

/** Add the energies due to the forcefields in this set to the molecules
    in the energy table 'energytable', scaled by 'scale_energy' */
void ForceFields::energy(EnergyTable &energytable, double scale_energy)
{
    this->energy(energytable, this->totalComponent(), scale_energy);
}


/** Add the force due to the component 'component' to the molecules
    in the force table 'forcetable', scaled by 'scale_force'
    
    \throw SireFF::missing_component
*/
void ForceFields::force(ForceTable &forcetable, const Symbol &component,
                        double scale_force)
{
    FFSymbolPtr comp = ffsymbols.value(component);

    if (comp.get() == 0)
        throw SireFF::missing_component( QObject::tr(   
            "There is no component of the energy represented by the "
            "symbol %1. Available components are %2.")
                .arg(component.toString(), Sire::toString(energySymbols())),
                    CODELOC );

    if (comp->isConstant())
        throw SireFF::missing_component( QObject::tr(
                "The component %1 is a constant component (it is not an energy "
                "component). Available energy components are %2.")
                    .arg(component.toString(),
                         Sire::toString(energySymbols())), CODELOC );

    comp->force(forcetable, ffields_by_idx, 
                ffsymbols, scale_force);
}

/** Add the forces due to the forcefields in this set to the molecules
    in the force table 'forcetable', scaled by 'scale_force' */
void ForceFields::force(ForceTable &forcetable, double scale_force)
{
    this->force(forcetable, this->totalComponent(), scale_force);
}

/** Add the field due to the component 'component' to the molecules
    in the field table 'fieldtable', scaled by 'scale_field'
    
    \throw SireFF::missing_component
*/
void ForceFields::field(FieldTable &fieldtable, const Symbol &component,
                        double scale_field)
{
    FFSymbolPtr comp = ffsymbols.value(component);

    if (comp.get() == 0)
        throw SireFF::missing_component( QObject::tr(   
            "There is no component of the energy represented by the "
            "symbol %1. Available components are %2.")
                .arg(component.toString(), Sire::toString(energySymbols())),
                    CODELOC );

    if (comp->isConstant())
        throw SireFF::missing_component( QObject::tr(
                "The component %1 is a constant component (it is not an energy "
                "component). Available energy components are %2.")
                    .arg(component.toString(),
                         Sire::toString(energySymbols())), CODELOC );

    comp->field(fieldtable, ffields_by_idx, 
                ffsymbols, scale_field);
}

/** Add the fields due to the forcefields in this set to the molecules
    in the field table 'fieldtable', scaled by 'scale_field' */
void ForceFields::field(FieldTable &fieldtable, double scale_field)
{
    this->field(fieldtable, this->totalComponent(), scale_field);
}

/** Add the field due to the component 'component' to the molecules
    in the field table 'fieldtable', scaled by 'scale_field'
    
    \throw SireFF::missing_component
*/
void ForceFields::field(FieldTable &fieldtable, const Symbol &component,
                        const Probe &probe, double scale_field)
{
    FFSymbolPtr comp = ffsymbols.value(component);

    if (comp.get() == 0)
        throw SireFF::missing_component( QObject::tr(   
            "There is no component of the energy represented by the "
            "symbol %1. Available components are %2.")
                .arg(component.toString(), Sire::toString(energySymbols())),
                    CODELOC );

    if (comp->isConstant())
        throw SireFF::missing_component( QObject::tr(
                "The component %1 is a constant component (it is not an energy "
                "component). Available energy components are %2.")
                    .arg(component.toString(),
                         Sire::toString(energySymbols())), CODELOC );

    comp->field(fieldtable, probe, ffields_by_idx, 
                ffsymbols, scale_field);
}

/** Add the fields due to the forcefields in this set to the molecules
    in the field table 'fieldtable', scaled by 'scale_field' */
void ForceFields::field(FieldTable &fieldtable, const Probe &probe, double scale_field)
{
    this->field(fieldtable, this->totalComponent(), probe, scale_field);
}

/** Add the potential due to the component 'component' to the molecules
    in the potential table 'pottable', scaled by 'scale_potential'
    
    \throw SireFF::missing_component
*/
void ForceFields::potential(PotentialTable &pottable, const Symbol &component,
                            double scale_potential)
{
    FFSymbolPtr comp = ffsymbols.value(component);

    if (comp.get() == 0)
        throw SireFF::missing_component( QObject::tr(   
            "There is no component of the energy represented by the "
            "symbol %1. Available components are %2.")
                .arg(component.toString(), Sire::toString(energySymbols())),
                    CODELOC );

    if (comp->isConstant())
        throw SireFF::missing_component( QObject::tr(
                "The component %1 is a constant component (it is not an energy "
                "component). Available energy components are %2.")
                    .arg(component.toString(),
                         Sire::toString(energySymbols())), CODELOC );

    comp->potential(pottable, ffields_by_idx, 
                    ffsymbols, scale_potential);
}

/** Add the potential due to the forcefields in this set to the molecules
    in the potential table 'pottable', scaled by 'scale_potential' */
void ForceFields::potential(PotentialTable &pottable, double scale_potential)
{
    this->potential(pottable, this->totalComponent(), scale_potential);
}

/** Add the potential due to the component 'component' to the molecules
    in the potential table 'pottable', scaled by 'scale_potential'
    
    \throw SireFF::missing_component
*/
void ForceFields::potential(PotentialTable &pottable, const Symbol &component,
                            const Probe &probe, double scale_potential)
{
    FFSymbolPtr comp = ffsymbols.value(component);

    if (comp.get() == 0)
        throw SireFF::missing_component( QObject::tr(   
            "There is no component of the energy represented by the "
            "symbol %1. Available components are %2.")
                .arg(component.toString(), Sire::toString(energySymbols())),
                    CODELOC );

    if (comp->isConstant())
        throw SireFF::missing_component( QObject::tr(
                "The component %1 is a constant component (it is not an energy "
                "component). Available energy components are %2.")
                    .arg(component.toString(),
                         Sire::toString(energySymbols())), CODELOC );

    comp->potential(pottable, probe, ffields_by_idx, 
                    ffsymbols, scale_potential);
}

/** Add the potential due to the forcefields in this set to the molecules
    in the potential table 'pottable', scaled by 'scale_potential' */
void ForceFields::potential(PotentialTable &pottable, 
                            const Probe &probe, double scale_potential)
{
    this->potential(pottable, this->totalComponent(), probe, scale_potential);
}

/** Sanitise the user properties. This removes dangling links,
    invalid combined properties and additional properties
    that have the same name as a built-in property. You need
    to call this whenever you change a property or add or
    remove a forcefield. This is because the because the
    dependencies of properties can become quite complicated...! */
void ForceFields::sanitiseUserProperties()
{
    bool something_changed = false;
    
    //need to keep iterating until nothing changes as there
    //can be quite complicated dependencies caused by links 
    //disappering causing combined properties to fall back
    //to a different real value, which may then make the
    //combined property invalid, which will then break
    //another link...!
    do
    {
        something_changed = false;

        if (not additional_properties.isEmpty())
        {
            //we cannot have an additional property that is also 
            //a built-in property
            foreach (QString key, additional_properties.keys())
            {
                for (QVector<FFPtr>::const_iterator it = ffields_by_idx.constBegin();
                     it != ffields_by_idx.constEnd();
                     ++it)
                {
                    if (it->read().containsProperty(key))
                    {
                        //this is a built-in property!
                        additional_properties.remove(key);
                        something_changed = true;
                        break;
                    }
                }
            }
        }

        if (not property_aliases.isEmpty())
        {
            foreach (QString key, property_aliases.keys())
            {
                if (not this->isValidLink(key))
                {
                    property_aliases.remove(key);
                    something_changed = true;
                }
            }
        }

        if (not combined_properties.isEmpty())
        {
            QSet<QString> updated_combinations;
        
            foreach (QString key, combined_properties.keys())
            {
                try
                {
                    if (not updated_combinations.contains(key))
                        this->updateCombinedProperty(key, true, &updated_combinations);
                }
                catch(...)
                {
                    //we can no longer update this combined property or
                    //one of its dependencies, so lets just remove it
                    combined_properties.remove(key);
                    something_changed = true;
                }
            }
        }
    
    } while (something_changed);
}

/** Remove the property with name 'name'. Note that this can only
    remove user-level properties - it cannot remove built-in properties
    of the forcefields. This does nothing if there is no user-level 
    property with this name */
void ForceFields::removeProperty(const QString &name)
{
    combined_properties.remove(name);
    property_aliases.remove(name);
    additional_properties.remove(name);

    this->sanitiseUserProperties();
}

/** Return whether or not the property 'name' exists and is a compound 
    property (either a link or a combined property) */
bool ForceFields::isCompoundProperty(const QString &name) const
{
    return combined_properties.contains(name) or 
           property_aliases.contains(name);
}

/** Return whether or not the property 'name' exists and is a user 
    supplied property (either a compound property or an extra
    ForceFields property) */
bool ForceFields::isUserProperty(const QString &name) const
{
    return this->isCompoundProperty(name) or additional_properties.contains(name);
}

/** Return whether or not the property 'name' exists and is a builtin
    property of one of the forcefields in this set */
bool ForceFields::isBuiltinProperty(const QString &name) const
{
    return (not this->isUserProperty(name)) and this->containsProperty(name);
}

/** Return the raw compound property with name 'name' - this returns
    the property representing the link, or the combined property,
    and raises an exception if a compound property with this name
    does not exist
    
    \throw SireBase::missing_property
*/
const Property& ForceFields::compoundProperty(const QString &name) const
{
    bool is_link = property_aliases.contains(name);
    bool is_combined = combined_properties.contains(name);
    
    if (is_link and is_combined)
        throw SireError::program_bug( QObject::tr(
            "How can the property %1 be both a link and a combined property?")
                .arg(name), CODELOC );
                
    else if (not (is_link or is_combined))
    {
        QStringList compound_props = property_aliases.keys() + 
                                     combined_properties.keys();
    
        if (not this->containsProperty(name))
            throw SireBase::missing_property( QObject::tr(
                    "There is no property %1 in this set of forcefields. "
                    "Available compound properties are [ %2 ].")
                        .arg(name, Sire::toString(compound_props)), CODELOC );

        else
            throw SireBase::missing_property( QObject::tr(
                    "The property %1 is not a compound property. "
                    "Available compound properties are [ %2 ].")
                        .arg(name, Sire::toString(compound_props)), CODELOC );
    }
    else if (is_link)
    {
        return property_aliases.constFind(name)->read();
    }
    else
        return combined_properties.constFind(name)->read();
}

/** Return the user-supplied property at 'name'. This raises an
    exception if there is no user-supplied property with this name
    
    \throw SireBase::missing_property
*/
const Property& ForceFields::userProperty(const QString &name) const
{
    if (not this->isUserProperty(name))
        throw SireBase::missing_property( QObject::tr(
            "There is no user property with name %1. Available user properties "
            "are [ %2 ].")
                .arg(name, Sire::toString( property_aliases.keys() +
                                           combined_properties.keys() + 
                                           additional_properties.keys() ) ), CODELOC );

    return this->property(name);
}

/** Return the built-in property at 'name'. This will by-pass any
    user-supplied property with this name, and will raise an
    exception if there is no built-in property with this name
    
    \throw SireBase::missing_property
*/
const Property& ForceFields::builtinProperty(const QString &name) const
{
    return this->property(FFIdentifier(), name);
}

/** Return all of the ForceFields level properties that are needed to
    get the value of the property with name 'name' - note that this
    does not include the names of ForceField level ("filtered") properties */
void ForceFields::getDependencies(const QString &name, QSet<QString> &deps) const
{
    if (property_aliases.contains(name))
    {
        const LinkToProperty &link = property_aliases.constFind(name)
                                                     ->read().asA<LinkToProperty>();
                                                     
        if (not link.isFiltered() and link.target().hasSource())
        {
            if (not deps.contains(link.target().source()))
            {
                deps.insert(link.target().source());
                this->getDependencies(link.target().source(), deps);
            }
        }
    }
    else if (combined_properties.contains(name))
    {
        const CombineProperties &combined = combined_properties.constFind(name)
                                                    ->read().asA<CombineProperties>();
        
        for (CombineProperties::const_iterator it = combined.constBegin();
             it != combined.constEnd();
             ++it)
        {
            if (it->hasSource())
            {
                if (not deps.contains(it->source()))
                {
                    deps.insert(it->source());
                    this->getDependencies(it->source(), deps);
                }
            }
        }
    }
}

/** Assert that the property 'name' does not contain a circular
    reference */
void ForceFields::assertNonCircularProperty(const QString &name) const
{
    QSet<QString> deps;
    
    this->getDependencies(name, deps);

    if (deps.contains(name))
        throw SireError::invalid_state( QObject::tr(
            "Circular link detected from %1 to %1. The dependencies "
            "of this property are %2.")
                .arg(name).arg( Sire::toString(deps) ), CODELOC );
}

/** Assert that the link at 'name' is valid */
void ForceFields::assertValidLink(const QString &name) const
{
    if (not property_aliases.contains(name))
        throw SireBase::missing_property( QObject::tr(
                "There is no link with name %1.").arg(name), CODELOC );
        
    const LinkToProperty &link = property_aliases.constFind(name)
                                                ->read().asA<LinkToProperty>();
    
    if (link.target().isNull())
        throw SireError::invalid_state( QObject::tr(
                "Null link detected from %1 -> NULL")
                    .arg(name), CODELOC );

    if (link.isFiltered())
    {
        if (not link.filter().isA<FFID>())
            throw SireError::invalid_cast( QObject::tr(
                    "You cannot not add a link that is filtered by a non FFID-derived "
                    "identifier (as you can only filter on forcefield IDs). "
                    "The link from %1 is %2.")
                        .arg(name, link.toString()), CODELOC );
    
        this->map(link.filter().asA<FFID>());
    }

    //see if we can get this property
    if (link.target().hasSource())
    {
        if (link.isFiltered())
        {
            //ensure that we can actually get the link
            this->property(link.filter().asA<FFID>(), link.target().source());
        }
        else
        {
            //ensure that this property is not a circular link                
            this->assertNonCircularProperty(name);

            //ensure that we can actually get the link
            this->property(link.target().source());
        }
    }
}

/** Return whether or not the specified link is valid */
bool ForceFields::isValidLink(const QString &name) const
{
    try
    {
        this->assertValidLink(name);
        return true;
    }
    catch(...)
    {
        return false;
    }
}

/** Update the combined property with name 'name' - we have
    already ensured that there are no circular references. This
    optionally record the set of combined properties that were updated
    when this was updated. This internal function is not atomic,
    so should only be called by another function that explicitly
    ensures atomicity */
void ForceFields::updateCombinedProperty(const QString &name, 
                                         bool update_dependencies,
                                         QSet<QString> *updated_combinations)
{
    if (updated_combinations)
    {
        if (updated_combinations->contains(name))
            return;
    }

    if (not combined_properties.contains(name))
        throw SireBase::missing_property( QObject::tr(
                "There is no property combination with name %1.").arg(name), CODELOC );
        
    CombineProperties &combined = combined_properties.find(name)
                                           ->edit().asA<CombineProperties>();

    //get all of the dependent properties of this combination
    Properties dependencies;
    
    for (CombineProperties::const_iterator it = combined.constBegin();
         it != combined.constEnd();
         ++it)
    {
        if (it->hasSource())
        {
            if (update_dependencies and combined_properties.contains(it->source()))
                this->updateCombinedProperty(it->source());
            
            dependencies.setProperty(it->source(), this->property(it->source()));
        }
    }

    combined.updateFrom(dependencies);
    
    if (updated_combinations)
        updated_combinations->insert(name);
}

/** Update all of the combined properties */
void ForceFields::updateCombinedProperties()
{
    if (combined_properties.isEmpty())
        return;

    else if (combined_properties.count() == 1)
    {
        this->updateCombinedProperty( combined_properties.constBegin().key(), false );
    }
    else
    {
        QHash<QString,PropertyPtr> old_combined_properties = combined_properties;
    
        try
        {
            QSet<QString> updated_combinations;
    
            foreach (QString key, combined_properties.keys())
            {
                this->updateCombinedProperty(key, true, &updated_combinations);
            }
        }
        catch(...)
        {
            combined_properties = old_combined_properties;
            throw;
        }
    }
}

/** Set the property 'name' to have the value 'value' in *all* of the 
    forcefields contained in this set
    
    \throw SireError::incompatible_error
    \throw SireError::invalid_cast
*/
void ForceFields::setProperty(const QString &name, const Property &value)
{
    if (name.isEmpty())
        return;

    //is this property a link? If so, then we need to forward this
    //setProperty to the target of the link
    if (property_aliases.contains(name))
    {
        const LinkToProperty &link = property_aliases.constFind(name)->read()
                                                        .asA<LinkToProperty>();
                                                        
        if (link.target().hasSource())
        {
            if (link.isFiltered())
                this->setProperty(link.filter().asA<FFID>(), 
                                  link.target().source(), value);
            else
                this->setProperty(link.target().source(), value);
        }
        
        return;
    }

    if (value.isA<LinkToProperty>())
    {
        //remove any existing link or combined property with this name
        QHash<QString,PropertyPtr> old_combined_properties = combined_properties;
        QHash<QString,PropertyPtr> old_property_aliases = property_aliases;
        QHash<QString,PropertyPtr> old_additional_properties = additional_properties;
        
        try
        {
            additional_properties.remove(name);
            combined_properties.remove(name);
            property_aliases.insert(name, value);

            this->assertValidLink(name);
            this->sanitiseUserProperties();
        }
        catch(...)
        {
            additional_properties = old_additional_properties;
            property_aliases = old_property_aliases;
            combined_properties = old_combined_properties;
            throw;
        }
    }
    else if (value.isA<CombineProperties>())
    {
        //remove any existing link or combined property with this name
        QHash<QString,PropertyPtr> old_additional_properties = additional_properties;
        QHash<QString,PropertyPtr> old_combined_properties = combined_properties;
        QHash<QString,PropertyPtr> old_property_aliases = property_aliases;
        
        try
        {
            additional_properties.remove(name);
            combined_properties.insert(name, value);
            property_aliases.remove(name);

            this->assertNonCircularProperty(name);
            this->sanitiseUserProperties();
        }
        catch(...)
        {
            additional_properties = old_additional_properties;
            property_aliases = old_property_aliases;
            combined_properties = old_combined_properties;
            throw;
        }
    }
    else
    {
        QVector<FFPtr> old_ffields_by_idx = ffields_by_idx;
        QHash<QString,PropertyPtr> old_additional_properties = additional_properties;
        QHash<QString,PropertyPtr> old_combined_properties = combined_properties;
        QHash<QString,PropertyPtr> old_property_aliases = property_aliases;

        try
        {
            if (additional_properties.contains(name))
            {
                //we already know that the forcefields don't contain
                //this property - update the additional property
                additional_properties.insert(name, value);
                this->sanitiseUserProperties();
            }
            else
            {
                //this is a new property that has replaced any link or combined property
                combined_properties.remove(name);
                property_aliases.remove(name);
                
                int nffields = ffields_by_idx.count();
                FFPtr *ffields_array = ffields_by_idx.data();
        
                bool changed_property = false;
                bool contains_property = false;
        
                for (int i=0; i<nffields; ++i)
                {
                    if (ffields_array[i].read().containsProperty(name))
                    {
                        contains_property = true;
                        
                        if (ffields_array[i].edit().setProperty(name, value))
                            changed_property = true;
                    }
                }
        
                if (changed_property)
                    this->sanitiseUserProperties();

                else if (not contains_property)
                    //there is no property in the forcefields, so add this
                    //as an additional property - no need to sanitise
                    //user properties as this is a totally new property
                    additional_properties.insert(name, value);
            }
        }
        catch(...)
        {
            ffields_by_idx = old_ffields_by_idx;
            additional_properties = old_additional_properties;
            combined_properties = old_combined_properties;
            property_aliases = old_property_aliases;
            throw;
        }
    }
}

/** Set the built-in property 'name' to have the value 'value' in all of the forcefields
    in this set that match the ID 'ffid'
    
    Note that because this operates on the level of individual forcefields,
    it operates only on built-in properties, not on user-supplied properties

    Note also that if this breaks any links or combined properties then
    the broken links and combined properties (including all those that
    depend on them) will be removed
    
    \throw SireFF::missing_forcefield
    \throw SireError::invalid_index
    \throw SireError::incompatible_error
    \throw SireError::invalid_cast
*/
void ForceFields::setProperty(const FFID &ffid, const QString &name, 
                              const Property &value)
{
    QVector<FFPtr> old_state( ffields_by_idx );

    try
    {
        QList<FFIdx> ffidxs = ffid.map(*this);
        
        bool changed_property = false;
        
        foreach (const FFIdx &ffidx, ffidxs)
        {
            if (this->_pvt_forceField(ffidx).containsProperty(name))
            {
                if (this->_pvt_forceField(ffidx).setProperty(name, value))
                    changed_property = true;
            }
        }
        
        if (changed_property)
            this->sanitiseUserProperties();
    }
    catch(...)
    {
        ffields_by_idx = old_state;
        throw;
    }
}

/** Return the value of the property with name 'name'. This returns
    the property if it exists in at least one forcefield, and
    if all occurances of the property have the same value

    \throw SireBase::duplicate_property
    \throw SireBase::missing_property
*/
const Property& ForceFields::property(const QString &name) const
{
    if (additional_properties.contains(name))
    {
        return additional_properties.constFind(name)->read();
    }
    else if (property_aliases.contains(name))
    {
        const LinkToProperty &link = property_aliases.constFind(name)->read()
                                                     .asA<LinkToProperty>();
    
        if (link.target().isNull())
            throw SireError::program_bug( QObject::tr(
                    "How did a null link from %1 get through???")
                        .arg(name), CODELOC );
    
        if (link.target().hasSource())
        {
            if (link.isFiltered())
            {
                return this->property(link.filter().asA<FFID>(),
                                      link.target().source());
            }
            else
            {
                if (link.target().source() == name)
                    throw SireError::program_bug( QObject::tr(
                        "How did the circular reference from %1 to %1 get through?")
                            .arg(name), CODELOC );

                return this->property(link.target().source());
            }
        }
        else
            return link.target().value();
    }
    else if (combined_properties.contains(name))
    {
        return combined_properties.constFind(name)->read().asA<CombineProperties>()
                                                          .combinedProperty();
    }
    else
    {
        //there is no user property with this name - try to find 
        //a built-in property with this name
    
        int nffields = ffields_by_idx.count();
        const FFPtr *ffields_array = ffields_by_idx.constData();
        
        const Property *p = &(Property::null());
        
        for (int i=0; i<nffields; ++i)
        {
            if (ffields_array[i]->containsProperty(name))
            {
                if ( p->equals(Property::null()) )
                {
                    p = &(ffields_array[i]->property(name));
                }
                else if ( not p->equals(ffields_array[i]->property(name)) )
                {
                    throw SireBase::duplicate_property( QObject::tr(
                        "More than one forcefield contains the property (%1), and "
                        "it has a different value in these forcefields. "
                        "You will have to search for this property "
                        "individually in each forcefield (e.g. using the "
                        "ForceFields::forceFieldsWithProperty(...) member function).")
                            .arg(name), CODELOC );
                }
            }
        }
        
        if ( p->equals(Property::null()) )
        {
            throw SireBase::missing_property( QObject::tr(
                "None of the contained forcefields have a property called %1. "
                "Available properties are %2.")
                    .arg(name, Sire::toString(this->propertyKeys()) ), CODELOC );
        }
        
        return *p;
    }
}

/** Return the value of the property 'name' in the forcefields identified
    by the ID 'ffid'
    
    Note that because this operates on the level of individual forcefields,
    it can only return built-in properties, and ignores any 
    user-supplied properties

    \throw SireBase::duplicate_property
    \throw SireFF::missing_forcefield
    \throw SireFF::duplicate_forcefield
    \throw SireBase::missing_property
*/
const Property& ForceFields::property(const FFID &ffid, const QString &name) const
{
    //get the list of forcefields that match this ID
    QList<FFIdx> ffidxs = ffid.map(*this);

    const Property *p = &(Property::null());

    bool has_property = false;

    foreach (FFIdx ffidx, ffidxs)
    {
        const FF &ff = this->_pvt_forceField(ffidx);
        
        if (ff.containsProperty(name))
        {
            if ( p->equals(Property::null()) )
            {
                p = &(ff.property(name));
                has_property = true;
            }
            else if ( not p->equals(ff.property(name)) )
            {
                throw SireBase::duplicate_property( QObject::tr(
                    "More than one of the forcefields that match the ID %1 "
                    "contain the property %2, and it has a different value "
                    "in at least two of these forcefields. You will need "
                    "to search for this property one forcefield at a time "
                    "(e.g. using the ForceFields::forceFieldsWithProperty(...) "
                    "member function).")
                        .arg(ffid.toString(), name), CODELOC );
            }
        }
    }
    
    if (not has_property)
        throw SireBase::missing_property( QObject::tr(
            "None of the forcefields that match the ID '%1' contain "
            "a property called %2.")
                .arg(ffid.toString(), name), CODELOC );

    return *p;
}

/** Return whether or not any of the forcefields contain the property
    with name 'name' */
bool ForceFields::containsProperty(const QString &name) const
{
    if (this->isUserProperty(name))
    {
        return true;
    }
    else
    {
        //look for a built-in property with this name
        int nffields = ffields_by_idx.count();
        const FFPtr *ffields_array = ffields_by_idx.constData();
    
        for (int i=0; i<nffields; ++i)
        {
            if (ffields_array[i]->containsProperty(name))
                return true;
        }
    
        return false;
    }
}

/** Return whether or not any of the forcefields that match the ID 'ffid'
    contain the property with name 'name'
    
    Note that because this operates on the level of individual forcefields,
    it can only return built-in properties, and ignores any 
    user-supplied properties
    
    \throw SireFF::missing_forcefield
    \throw SireError::invalid_index
*/
bool ForceFields::containsProperty(const FFID &ffid, const QString &name) const
{
    QList<FFIdx> ffidxs = ffid.map(*this);
    
    foreach (const FFIdx &ffidx, ffidxs)
    {
        if (this->_pvt_forceField(ffidx).containsProperty(name))
            return true;
    }
    
    return false;
}

/** Return whether or not any of the forcefields contain the property
    with name 'name' */
bool ForceFields::containsProperty(const PropertyName &name) const
{
    if (name.hasValue())
        return true;
    else
        return this->containsProperty(name.source());
}

/** Return whether or not any of the forcefields that match the ID 'ffid'
    contain the property with name 'name'
    
    Note that because this operates on the level of individual forcefields,
    it can only return built-in properties, and ignores any 
    user-supplied properties
    
    \throw SireFF::missing_forcefield
    \throw SireError::invalid_index
*/
bool ForceFields::containsProperty(const FFID &ffid, const PropertyName &name) const
{
    if (name.hasValue())
        return true;
    else
        return this->containsProperty(ffid, name.source());
}

/** Return the names of all of the properties in all of the forcefields */
QStringList ForceFields::propertyKeys() const
{
    QSet<QString> keys;
    
    keys.unite( property_aliases.keys().toSet() );
    keys.unite( combined_properties.keys().toSet() );
    keys.unite( additional_properties.keys().toSet() );

    int nffields = this->nForceFields();

    if (nffields == 0)
    {
        return QStringList( keys.toList() );
    }
    else if (nffields == 1)
    {
        keys.unite( this->_pvt_forceField(0).propertyKeys().toSet() );
        return QStringList( keys.toList() );
    }
    
    for (int i=0; i<nffields; ++i)
    {
        keys.unite( this->_pvt_forceField(i).propertyKeys().toSet() );
    }
    
    return QStringList( keys.toList() );
}

/** Return the names of all of the properties in the forcefields 
    identified by the ID 'ffid'
    
    Note that because this operates on the level of individual forcefields,
    it can only return built-in properties, and ignores any 
    user-supplied properties
    
    \throw SireFF::missing_forcefield
    \throw SireError::invalid_index
*/
QStringList ForceFields::propertyKeys(const FFID &ffid) const
{
    QList<FFIdx> ffidxs = ffid.map(*this);
    
    if (ffidxs.count() == 1)
    {
        return this->_pvt_forceField(ffidxs.at(0)).propertyKeys();
    }
    
    QSet<QString> keys;
    
    foreach (FFIdx idx, ffidxs)
    {
        keys.unite( this->_pvt_forceField(idx).propertyKeys().toSet() );
    }
    
    return QStringList( keys.toList() );
}

/** Return all of the properties in all of the forcefields. This will raise
    an error if there are properties with the same name in different 
    forcefields that have different values.
    
    \throw SireBase::duplicate_property
*/
Properties ForceFields::properties() const
{
    Properties props;
    
    QStringList keys = this->propertyKeys();
    
    foreach (QString key, keys)
    {
        props.setProperty(key, this->property(key));
    }
    
    return props;
}

/** Return all of the properties in all of the forcefields identified by
    the ID 'ffid'. This will raise an error if there are properties with
    the same name in different forcefields that have different values.
    
    Note that because this operates on the level of individual forcefields,
    it can only return built-in properties, and ignores any 
    user-supplied properties
    
    \throw SireBase::duplicate_property
*/
Properties ForceFields::properties(const FFID &ffid) const
{
    QList<FFIdx> ffidxs = ffid.map(*this);
    
    if (ffidxs.count() == 1)
    {
        return this->_pvt_forceField(ffidxs.at(0)).properties();
    }
    
    QStringList keys = this->propertyKeys(ffid);
    
    Properties props;
    
    foreach (QString key, keys)
    {
        props.setProperty(key, this->property(ffid, key));
    }
    
    return props;
}
    
/** Return all of the user-supplied properties in this set of forcefields */
Properties ForceFields::userProperties() const
{
    Properties props;
    
    foreach (QString key, additional_properties.keys())
    {
        props.setProperty(key, this->property(key));
    }
    
    foreach (QString key, property_aliases.keys())
    {
        props.setProperty(key, this->property(key));
    }
    
    foreach (QString key, combined_properties.keys())
    {
        props.setProperty(key, this->property(key));
    }
    
    return props;
}

/** Return all of the built-in properties of the forcefields in this set */
Properties ForceFields::builtinProperties() const
{
    return this->properties(FFIdentifier());
}

/** Return the list of all forcefields that contain a property with name 'name'
    
    Note that because this operates on the level of individual forcefields,
    it can only return built-in properties, and ignores any 
    user-supplied properties

    \throw SireBase::missing_property
*/
QVector<FFPtr> ForceFields::forceFieldsWithProperty(const QString &name) const
{
    QVector<FFPtr> ffs;
    
    int nffields = this->nForceFields();
    
    for (int i=0; i<nffields; ++i)
    {
        if (this->_pvt_forceField(i).containsProperty(name))
            ffs.append( this->_pvt_forceField(i) );
    }
    
    if (ffs.isEmpty())
        throw SireBase::missing_property( QObject::tr(
            "No forcefields contain the property %1. Available properties "
            "are %2.")
                .arg(name, Sire::toString(this->propertyKeys())), CODELOC );
                
    return ffs;
}

/** Return the list of forcefields that match the ID "ffid" and that
    contain the property with name 'name'
    
    Note that because this operates on the level of individual forcefields,
    it can only return built-in properties, and ignores any 
    user-supplied properties
    
    \throw SireFF::missing_forcefield
    \throw SireBase::missing_property
    \throw SireError::invalid_index
*/
QVector<FFPtr> ForceFields::forceFieldsWithProperty(const FFID &ffid, 
                                                    const QString &name) const
{
    QList<FFIdx> ffidxs = ffid.map(*this);
    
    QVector<FFPtr> ffs;
    
    foreach (FFIdx ffidx, ffidxs)
    {
        if (this->_pvt_forceField(ffidx).containsProperty(name))
            ffs.append( this->_pvt_forceField(ffidx) );
    }
    
    if (ffs.isEmpty())
        throw SireBase::missing_property( QObject::tr(
            "None of the forcefields that match the ID %1 contain the "
            "property called %2. Available properties are %3.")
                .arg(ffid.toString(), name, 
                     Sire::toString(this->propertyKeys(ffid))), CODELOC );
                     
    return ffs;
}
   
/** Return the list of all forcefields that match the ID 'ffid'

    \throw SireFF::missing_forcefield
    \throw SireError::invalid_index
*/
QVector<FFPtr> ForceFields::forceFields(const FFID &ffid) const
{
    QList<FFIdx> ffidxs = ffid.map(*this);
    
    QVector<FFPtr> ffs;
    ffs.reserve(ffidxs.count());
    
    foreach (FFIdx ffidx, ffidxs)
    {
        ffs.append( this->_pvt_forceField(ffidx) );
    }
    
    return ffs;
}

/** Return the names of all of the forcefields that match the ID 'ffid'

    \throw SireFF::missing_forcefield
    \throw SireError::invalid_index
*/
QList<FFName> ForceFields::ffNames(const FFID &ffid) const
{
    QList<FFIdx> ffidxs = ffid.map(*this);

    QList<FFName> names;
    
    foreach (FFIdx ffidx, ffidxs)
    {
        names.append( this->_pvt_forceField(ffidx).name() );
    }
    
    return names;
}

/** Return an array containing all of the forcefields in this set, ordered
    in the same order as they appear in this set */
const QVector<FFPtr>& ForceFields::forceFields() const
{
    return ffields_by_idx;
}

/** Return an array containing all of the forcefields in this set, ordered
    in the same order as they appear in this set */
const QVector<FFPtr>& ForceFields::list() const
{
    return this->forceFields();
}

/** Return a list of all of the forcefield names, ordered in the same
    order as the forcefields appear in this set */
QList<FFName> ForceFields::ffNames() const
{
    QList<FFName> ffnames;
    
    for (QHash<QString,int>::const_iterator it = ffields_by_name.constBegin();
         it != ffields_by_name.constEnd();
         ++it)
    {
        ffnames.append( FFName(it.key()) );
    }

    return ffnames;
}

/** Return a list of all of the forcefield names, ordered in the same
    order as the forcefields appear in this set */
QList<FFName> ForceFields::names() const
{
    return this->ffNames();
}

/** Tell all of the forcefields that they must now recalculate their
    energies from scratch. This is a good way to debug the forcefields,
    but may also speed up cases where you know in advance that you will
    be moving most (or all) of the molecules between energy calculations */
void ForceFields::mustNowRecalculateFromScratch()
{
    int nffields = ffields_by_idx.count();
    FFPtr *ffields_array = ffields_by_idx.data();
    
    for (int i=0; i<nffields; ++i)
    {
        ffields_array[i].edit().mustNowRecalculateFromScratch();
    }
}

/** Return whether or not any of the forcefields in this set are dirty
    (the molecules have changed since the last energy calculation) */
bool ForceFields::isDirty() const
{
    int nffields = ffields_by_idx.count();
    const FFPtr *ffields_array = ffields_by_idx.constData();
    
    for (int i=0; i<nffields; ++i)
    {
        if (ffields_array[i]->isDirty())
            return true;
    }
    
    return false;
}

/** Return whether or not all of the forcefields in this set are clean
    (there have been no changes since the last energy evaluation) */
bool ForceFields::isClean() const
{
    return not this->isDirty();
}

/** Add the forcefield 'forcefield' to this set. This will raise
    an exception if this forcefield (or one with the same name)
    is already present in this set. Note that if the added
    forcefield will be updated to contain the versions of 
    any molecules that are already present in any of the 
    other forcefields.

    \throw SireFF::duplicate_forcefield
    \throw SireMol::duplicate_group
*/
void ForceFields::add(const FF &forcefield)
{
    FFPtr ff( forcefield );
    ff.edit().update( this->matchToExistingVersion(forcefield.molecules()) );

    ForceFields old_state( *this );
    
    try
    {
        ffields_by_idx.append(ff);
        this->rebuildIndex();
        
        //adding the forcefield may break some links
        this->sanitiseUserProperties();
    }
    catch(...)
    {
        this->operator=(old_state);
        throw;
    }
}

/** Internal function used to remove the ith forcefield */
void ForceFields::_pvt_remove(int i)
{
    ForceFields old_state( *this );
    
    try
    {
        ffields_by_idx.remove(i);
        this->rebuildIndex();

        //removing the forcefield may break some links
        this->sanitiseUserProperties();
    }
    catch(...)
    {
        this->operator=(old_state);
        throw;
    }
}

/** Remove the forcefield at index 'ffidx' from this set

    \throw SireError::invalid_index
*/
void ForceFields::remove(const FFIdx &ffidx)
{
    this->_pvt_remove( ffidx.map(ffields_by_idx.count()) );
}

/** Remove the forcefield with name 'ffname'.

    \throw SireFF::missing_forcefield
*/
void ForceFields::remove(const FFName &ffname)
{
    this->_pvt_remove( this->ffIdx(ffname) );
}

/** Remove the forcefield(s) that match the ID 'ffid' 

    \throw SireFF::missing_forcefield
    \throw SireError::invalid_index
*/
void ForceFields::remove(const FFID &ffid)
{
    QList<FFIdx> ffidxs = ffid.map(*this);
    
    ForceFields old_state( *this );
    
    try
    {
        foreach (const FFIdx &ffidx, ffidxs)
        {
            this->_pvt_remove(ffidx);
        }
        
        this->rebuildIndex();

        //removing the forcefield may break some links
        this->sanitiseUserProperties();
    }
    catch(...)
    {
        this->operator=(old_state);
        throw;
    }
}

/** Remove all of the forcefields from this set */
void ForceFields::removeAllForceFields()
{
    ForceFields old_state(*this);
    
    try
    {
        ffields_by_idx.clear();
        ffields_by_name.clear();
        mgroups_by_num.clear();

        property_aliases.clear();
        combined_properties.clear();

        this->clearIndex();
        this->rebuildIndex();
    }
    catch(...)
    {
        this->operator=(old_state);
        throw;
    }
}

/** Return the molecule group that has number 'mgnum'

    \throw SireMol::missing_group
*/
const MoleculeGroup& ForceFields::at(MGNum mgnum) const
{
    if (not mgroups_by_num.contains(mgnum))
        throw SireMol::missing_group( QObject::tr(
            "None of the forcefields in this set contain a molecule group "
            "with number %1. Available groups have numbers %2.")
                .arg(mgnum).arg( Sire::toString(mgroups_by_num.keys()) ),
                    CODELOC );
                    
    return this->_pvt_forceField(mgnum).at(mgnum);
}

/** Add the molecule view 'molview' to the molecule groups identified 
    by 'mgid', using the supplied property map to find the properties
    required for the forcefield
    
    \throw SireMol::missing_group
    \throw SireBase::missing_parameter
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void ForceFields::add(const MoleculeView &molview, const MGID &mgid,
                      const PropertyMap &map)
{
    QList<MGNum> mgnums = mgid.map(*this);

    PartialMolecule view(molview);
    view.update( this->matchToExistingVersion(view.data()) );
    
    ForceFields old_state( *this );
    
    try
    {
        foreach (const MGNum &mgnum, mgnums)
        {
            this->_pvt_forceField(mgnum).add(view, mgnum, map);
            MolGroupsBase::addToIndex(mgnum, view.data().number());
        }
    }
    catch(...)
    {
        this->operator=(old_state);
        throw;
    }
}
 
/** Add the views of the molecule in 'molviews' to the molecule groups identified 
    by 'mgid', using the supplied property map to find the properties
    required for the forcefield
    
    \throw SireMol::missing_group
    \throw SireBase::missing_parameter
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void ForceFields::add(const ViewsOfMol &molviews, const MGID &mgid,
                      const PropertyMap &map)
{
    QList<MGNum> mgnums = mgid.map(*this);

    ViewsOfMol views(molviews);
    views.update( this->matchToExistingVersion(views.data()) );
    
    ForceFields old_state( *this );
    
    try
    {
        foreach (const MGNum mgnum, mgnums)
        {
            this->_pvt_forceField(mgnum).add(views, mgnum, map);
            MolGroupsBase::addToIndex(mgnum, views.number());
        }
    }
    catch(...)
    {
        this->operator=(old_state);
        throw;
    }
}

/** Add the molecules in 'molecules' to the molecule groups identified 
    by 'mgid', using the supplied property map to find the properties
    required for the forcefield
    
    \throw SireMol::missing_group
    \throw SireBase::missing_parameter
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void ForceFields::add(const Molecules &molecules, const MGID &mgid,
                      const PropertyMap &map)
{
    QList<MGNum> mgnums = mgid.map(*this);
    
    Molecules mols = this->matchToExistingVersion(molecules);

    ForceFields old_state( *this );
    
    try
    {
        foreach (const MGNum &mgnum, mgnums)
        {
            this->_pvt_forceField(mgnum).add(mols, mgnum, map);

            MolGroupsBase::addToIndex(mgnum, mols.molNums());
        }
    }
    catch(...)
    {
        this->operator=(old_state);
        throw;
    }
}

/** Add the molecules in the molecule group 'molgroup' to the molecule 
    groups identified by 'mgid', using the supplied property map to find the 
    properties required for the forcefield
    
    \throw SireMol::missing_group
    \throw SireBase::missing_parameter
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void ForceFields::add(const MoleculeGroup &molgroup, const MGID &mgid,
                      const PropertyMap &map)
{
    QList<MGNum> mgnums = mgid.map(*this);
    
    //update the group to match the molecule versions
    //already present in this group...
    MolGroupPtr group(molgroup);
    group.edit().update( this->matchToExistingVersion(group.read().molecules()) );

    ForceFields old_state( *this );
    
    try
    {
        foreach (const MGNum &mgnum, mgnums)
        {
            this->_pvt_forceField(mgnum).add(group, mgnum, map);
            
            MolGroupsBase::addToIndex(mgnum, group->molNums().toList().toSet());
        }
    }
    catch(...)
    {
        this->operator=(old_state);
        throw;
    }
}

/** Add the molecule view 'molview' to the molecule groups identified 
    by 'mgid', using the supplied property map to find the properties
    required for the forcefield. This only adds the view to groups 
    that don't already contain this view.
    
    \throw SireMol::missing_group
    \throw SireBase::missing_parameter
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void ForceFields::addIfUnique(const MoleculeView &molview, const MGID &mgid,
                              const PropertyMap &map)
{
    QList<MGNum> mgnums = mgid.map(*this);

    PartialMolecule view(molview);
    view.update( this->matchToExistingVersion(view.data()) );
    
    ForceFields old_state( *this );
    
    try
    {
        foreach (const MGNum &mgnum, mgnums)
        {
            this->_pvt_forceField(mgnum).addIfUnique(view, mgnum, map);
                                                     
            MolGroupsBase::addToIndex(mgnum, view.data().number());
        }
    }
    catch(...)
    {
        this->operator=(old_state);
        throw;
    }
}

/** Add the molecule views in 'molviews' to the molecule groups identified 
    by 'mgid', using the supplied property map to find the properties
    required for the forcefield. This only adds the view to groups 
    that don't already contain this view.
    
    \throw SireMol::missing_group
    \throw SireBase::missing_parameter
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void ForceFields::addIfUnique(const ViewsOfMol &molviews, const MGID &mgid,
                              const PropertyMap &map)
{
    QList<MGNum> mgnums = mgid.map(*this);

    ViewsOfMol views(molviews);
    views.update( this->matchToExistingVersion(views.data()) );
    
    ForceFields old_state( *this );
    
    try
    {
        foreach (const MGNum &mgnum, mgnums)
        {
            this->_pvt_forceField(mgnum).addIfUnique(views, mgnum, map);
            
            MolGroupsBase::addToIndex(mgnum, views.number());
        }
    }
    catch(...)
    {
        this->operator=(old_state);
        throw;
    }
}
                              
/** Add the molecules in 'molecules' to the molecule groups identified 
    by 'mgid', using the supplied property map to find the properties
    required for the forcefield. This only adds the view to groups 
    that don't already contain this view.
    
    \throw SireMol::missing_group
    \throw SireBase::missing_parameter
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void ForceFields::addIfUnique(const Molecules &molecules, const MGID &mgid,
                              const PropertyMap &map)
{
    QList<MGNum> mgnums = mgid.map(*this);
    
    Molecules mols = this->matchToExistingVersion(molecules);
    QSet<MolNum> molnums = mols.molNums();

    ForceFields old_state( *this );
    
    try
    {
        foreach (const MGNum &mgnum, mgnums)
        {
            this->_pvt_forceField(mgnum).addIfUnique(mols, mgnum, map);
            
            MolGroupsBase::addToIndex(mgnum, mols.molNums());
        }
    }
    catch(...)
    {
        this->operator=(old_state);
        throw;
    }
}

/** Add the molecules in the molecule group 'molgroup' to the molecule 
    groups identified by 'mgid', using the supplied property map to 
    find the properties required for the forcefield. This only adds the 
    view to groups that don't already contain this view.
    
    \throw SireMol::missing_group
    \throw SireBase::missing_parameter
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void ForceFields::addIfUnique(const MoleculeGroup &molgroup, const MGID &mgid,
                              const PropertyMap &map)
{
    QList<MGNum> mgnums = mgid.map(*this);

    //update the group...
    MolGroupPtr group(molgroup);
    group.edit().update( this->matchToExistingVersion(molgroup.molecules()) );
    
    ForceFields old_state( *this );
    
    try
    {
        foreach (const MGNum &mgnum, mgnums)
        {
            this->_pvt_forceField(mgnum).addIfUnique(group, mgnum, map);
            
            MolGroupsBase::addToIndex(mgnum, group.read().molNums().toList().toSet());
        }
    }
    catch(...)
    {
        this->operator=(old_state);
        throw;
    }
}

/** Add the molecule view 'molview' to the molecule groups identified 
    by 'mgid', using the default locations to find any properties 
    required by the forcefields. 
        
    \throw SireMol::missing_group
    \throw SireBase::missing_parameter
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void ForceFields::add(const MoleculeView &molview, const MGID &mgid)
{
    this->add(molview, mgid, PropertyMap());
}

/** Add the molecule views in 'molviews' to the molecule groups identified 
    by 'mgid', using the supplied property map to find the properties
    required for the forcefield
    
    \throw SireMol::missing_group
    \throw SireBase::missing_parameter
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void ForceFields::add(const ViewsOfMol &molviews, const MGID &mgid)
{
    this->add(molviews, mgid, PropertyMap());
}

/** Add the molecules in 'molecules' to the molecule groups identified 
    by 'mgid', using the supplied property map to find the properties
    required for the forcefield
    
    \throw SireMol::missing_group
    \throw SireBase::missing_parameter
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void ForceFields::add(const Molecules &molecules, const MGID &mgid)
{
    this->add(molecules, mgid, PropertyMap());
}

/** Add the molecules in the molecule group 'molgroup' to the molecule 
    groups identified by 'mgid', using the supplied property map to 
    find the properties required for the forcefield
    
    \throw SireMol::missing_group
    \throw SireBase::missing_parameter
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void ForceFields::add(const MoleculeGroup &molgroup, const MGID &mgid)
{
    this->add(molgroup, mgid, PropertyMap());
}

/** Add the molecule view 'molview' to the molecule groups identified 
    by 'mgid', using the supplied property map to find the properties
    required for the forcefield. Only views that aren't already in
    the forcefield are added.
    
    \throw SireMol::missing_group
    \throw SireBase::missing_parameter
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void ForceFields::addIfUnique(const MoleculeView &molview, const MGID &mgid)
{
    this->addIfUnique(molview, mgid, PropertyMap());
}

/** Add the molecule views in 'molviews' to the molecule groups identified 
    by 'mgid', using the supplied property map to find the properties
    required for the forcefield. Only views that aren't already in
    the forcefield are added.
    
    \throw SireMol::missing_group
    \throw SireBase::missing_parameter
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void ForceFields::addIfUnique(const ViewsOfMol &molviews, const MGID &mgid)
{
    this->addIfUnique(molviews, mgid, PropertyMap());
}

/** Add the molecules in 'molecules' to the molecule groups identified 
    by 'mgid', using the supplied property map to find the properties
    required for the forcefield. Only views that aren't already in
    the forcefield are added.
    
    \throw SireMol::missing_group
    \throw SireBase::missing_parameter
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void ForceFields::addIfUnique(const Molecules &molecules, const MGID &mgid)
{
    this->addIfUnique(molecules, mgid, PropertyMap());
}

/** Add the molecules in the molecule group 'molgroup' to the molecule 
    groups identified by 'mgid', using the supplied property map to 
    find the properties required for the forcefield. Only views that 
    aren't already in the forcefield are added.
    
    \throw SireMol::missing_group
    \throw SireBase::missing_parameter
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void ForceFields::addIfUnique(const MoleculeGroup &molgroup, const MGID &mgid)
{
    this->addIfUnique(molgroup, mgid, PropertyMap());
}

/** Remove all of the molecule views that are contained in the molecule
    groups identified by the ID 'mgid'
    
    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
bool ForceFields::removeAll(const MGID &mgid)
{
    QList<MGNum> mgnums = mgid.map(*this);
    
    ForceFields old_state( *this );
    
    bool mols_removed = false;
    
    try
    {
        foreach (const MGNum &mgnum, mgnums)
        {
            if (not mols_removed)
            {
                if (not this->at(mgnum).isEmpty())
                    mols_removed = true;
            }
        
            this->_pvt_forceField(mgnum).removeAll(mgnum);
            
            MolGroupsBase::clearIndex(mgnum);
        }
    }
    catch(...)
    {
        this->operator=(old_state);
        throw;
    }
    
    return mols_removed;
}

/** Remove the view 'molview' from the specified groups in this
    forcefield. Note that this only removes the specific view
    (and indeed only the first copy of this view if there 
    are duplicates) - it does not remove the atoms in this
    view from all of the other views
    
    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
bool ForceFields::remove(const MoleculeView &molview, const MGID &mgid)
{
    QList<MGNum> mgnums = mgid.map(*this);
    
    ForceFields old_state( *this );
    
    bool mols_removed = false;
    
    try
    {
        foreach (const MGNum &mgnum, mgnums)
        {
            FF &ff = this->_pvt_forceField(mgnum);
            
            bool mol_removed = ff.remove(molview, mgnum);
            
            mols_removed = mols_removed or mol_removed;
            
            if (not ff.group(mgnum).contains(molview.data().number()))
            {
                MolGroupsBase::removeFromIndex(mgnum, molview.data().number());
            }
        }
    }
    catch(...)
    {
        this->operator=(old_state);
        throw;
    }
    
    return mols_removed;
}

/** Remove the views in 'molviews' from the specified groups in this
    forcefield. Note that this only removes the specific views
    (and indeed only the first copy of this view if there 
    are duplicates) - it does not remove the atoms in this
    view from all of the other views
    
    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
bool ForceFields::remove(const ViewsOfMol &molviews, const MGID &mgid)
{
    QList<MGNum> mgnums = mgid.map(*this);
    
    ForceFields old_state( *this );
    
    bool mols_removed = false;
    
    try
    {
        foreach (const MGNum &mgnum, mgnums)
        {
            FF &ff = this->_pvt_forceField(mgnum);
            
            bool mol_removed = ff.remove(molviews, mgnum);
            
            mols_removed = mols_removed or mol_removed;
            
            if (not ff.group(mgnum).contains(molviews.number()))
            {
                MolGroupsBase::removeFromIndex(mgnum, molviews.number());
            }
        }
    }
    catch(...)
    {
        this->operator=(old_state);
        throw;
    }
    
    return mols_removed;
}

/** Remove them molecules in 'molecules' from the specified groups in this
    forcefield. Note that this only removes the specific views
    (and indeed only the first copy of this view if there 
    are duplicates) - it does not remove the atoms in this
    view from all of the other views
    
    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
bool ForceFields::remove(const Molecules &molecules, const MGID &mgid)
{
    QList<MGNum> mgnums = mgid.map(*this);
    
    ForceFields old_state( *this );
    
    bool mols_removed = false;
    
    try
    {
        foreach (const MGNum &mgnum, mgnums)
        {
            FF &ff = this->_pvt_forceField(mgnum);
            
            bool mol_removed = ff.remove(molecules, mgnum);
            
            mols_removed = mols_removed or mol_removed;
            
            const MoleculeGroup &molgroup = ff.group(mgnum);
            
            for (Molecules::const_iterator it = molecules.constBegin();
                 it != molecules.constEnd();
                 ++it)
            {
                if (not molgroup.contains(it->number()))
                {
                    MolGroupsBase::removeFromIndex(mgnum, it->number());
                }
            }
        }
    }
    catch(...)
    {
        this->operator=(old_state);
        throw;
    }
    
    return mols_removed;
}

/** Remove the views in the molecule group 'molgroup' from the specified 
    groups in this forcefield. Note that this only removes the specific views
    (and indeed only the first copy of this view if there 
    are duplicates) - it does not remove the atoms in this
    view from all of the other views
    
    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
bool ForceFields::remove(const MoleculeGroup &molgroup, const MGID &mgid)
{
    return this->remove(molgroup.molecules(), mgid);
}

/** Remove the all copies of the view in 'molview' from the specified 
    groups in this forcefield. Note that this only removes the specific views
    - it does not remove the atoms in this view from all of the other views
    
    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
bool ForceFields::removeAll(const MoleculeView &molview, const MGID &mgid)
{
    QList<MGNum> mgnums = mgid.map(*this);
    
    ForceFields old_state( *this );
    
    bool mols_removed = false;
    
    try
    {
        foreach (const MGNum &mgnum, mgnums)
        {
            FF &ff = this->_pvt_forceField(mgnum);
            
            bool mol_removed = ff.removeAll(molview, mgnum);
            
            mols_removed = mols_removed or mol_removed;
            
            if (not ff.group(mgnum).contains(molview.data().number()))
            {
                MolGroupsBase::removeFromIndex(mgnum, molview.data().number());
            }
        }
    }
    catch(...)
    {
        this->operator=(old_state);
        throw;
    }
    
    return mols_removed;
}

/** Remove the all copies of the views in 'molviews' from the specified 
    groups in this forcefield. Note that this only removes the specific views
    - it does not remove the atoms in this view from all of the other views
    
    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
bool ForceFields::removeAll(const ViewsOfMol &molviews, const MGID &mgid)
{
    QList<MGNum> mgnums = mgid.map(*this);
    
    ForceFields old_state( *this );
    
    bool mols_removed = false;
    
    try
    {
        foreach (const MGNum &mgnum, mgnums)
        {
            FF &ff = this->_pvt_forceField(mgnum);
            
            bool mol_removed = ff.removeAll(molviews, mgnum);
            
            mols_removed = mols_removed or mol_removed;
            
            if (not ff.group(mgnum).contains(molviews.number()))
            {
                MolGroupsBase::removeFromIndex(mgnum, molviews.number());
            }
        }
    }
    catch(...)
    {
        this->operator=(old_state);
        throw;
    }
    
    return mols_removed;
}

/** Remove the all copies of the molecules in 'molecules' from the specified 
    groups in this forcefield. Note that this only removes the specific views
    - it does not remove the atoms in this view from all of the other views
    
    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
bool ForceFields::removeAll(const Molecules &molecules, const MGID &mgid)
{
    QList<MGNum> mgnums = mgid.map(*this);
    
    ForceFields old_state( *this );
    
    bool mols_removed = false;
    
    try
    {
        foreach (const MGNum &mgnum, mgnums)
        {
            FF &ff = this->_pvt_forceField(mgnum);
            
            bool mol_removed = ff.removeAll(molecules, mgnum);
            
            mols_removed = mols_removed or mol_removed;
            
            const MoleculeGroup &group = ff.group(mgnum);
            
            for (Molecules::const_iterator it = molecules.constBegin();
                 it != molecules.constEnd();
                 ++it)
            {
                if (not group.contains(it->number()))
                {
                    MolGroupsBase::removeFromIndex(mgnum, it->number());
                }
            }
        }
    }
    catch(...)
    {
        this->operator=(old_state);
        throw;
    }
    
    return mols_removed;
}

/** Remove the all copies of the molecules in the molecule group 'molgroup' 
    from the specified groups in this forcefield. Note that this only removes
    the specific views - it does not remove the atoms in this view from all 
    of the other views
    
    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
bool ForceFields::removeAll(const MoleculeGroup &molgroup, const MGID &mgid)
{
    return this->removeAll(molgroup.molecules(), mgid);
}

/** Remove all views of the molecule with number 'molnum' from the molecule
    groups identified by 'mgid'
    
    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
bool ForceFields::remove(MolNum molnum, const MGID &mgid)
{
    QList<MGNum> mgnums = mgid.map(*this);
    
    if (mgnums.count() == 1)
    {
        MGNum mgnum = mgnums.at(0);
        
        if (this->_pvt_forceField(mgnum).remove(molnum, mgnum))
        {
            MolGroupsBase::removeFromIndex(mgnum, molnum);
            return true;
        }
        else
            return false;
    }
    else
    {
        ForceFields old_state( *this );
    
        bool mols_removed = false;
    
        try
        {
            foreach (const MGNum &mgnum, mgnums)
            {
                if (this->_pvt_forceField(mgnum).remove(molnum, mgnum))
                {
                    MolGroupsBase::removeFromIndex(mgnum, molnum);
                    mols_removed = true;
                }
            }
        }
        catch(...)
        {
            this->operator=(old_state);
            throw;
        }
        
        return mols_removed;
    }
}

/** Remove all of the molecules whose numbers are in 'molnums' from
    all of the molecule groups identified by the ID 'mgid'
    
    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
bool ForceFields::remove(const QSet<MolNum> &molnums, const MGID &mgid)
{
    QList<MGNum> mgnums = mgid.map(*this);
    
    ForceFields old_state( *this );

    bool mols_removed = false;
    
    try
    {
        foreach (const MGNum &mgnum, mgnums)
        {
            if (this->_pvt_forceField(mgnum).remove(molnums, mgnum))
            {
                MolGroupsBase::removeFromIndex(mgnum, molnums);
                mols_removed = true;
            }
        }
    }
    catch(...)
    {
        this->operator=(old_state);
        throw;
    }
    
    return mols_removed;
}

/** Update all of the forcefields so that they use the version of the data 
    of the molecule held in 'moldata'
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void ForceFields::update(const MoleculeData &moldata, bool auto_commit)
{
    if (not this->contains(moldata.number()))
        return;

    if (auto_commit and this->needsAccepting())
        this->accept();
    
    const QList<MGNum> &mgnums = this->groupsContaining(moldata.number());

    BOOST_ASSERT(not mgnums.isEmpty());

    if (mgnums.count() == 1)
    {
        this->_pvt_forceField(mgnums.at(0)).update(moldata,auto_commit);
    }
    else
    {
        foreach (const MGNum &mgnum, mgnums)
        {
            this->_pvt_forceField(mgnum).update(moldata,auto_commit);
        }
    }

    if (auto_commit and this->needsAccepting())
        this->accept();
}

/** Update all of the forcefields in this group so that they have the 
    same version of the molecules that are present in 'molecules'
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void ForceFields::update(const Molecules &molecules, bool auto_commit)
{
    if (molecules.isEmpty())
        return;
    
    else if (molecules.count() == 1)
    {
        this->update( molecules.constBegin()->data(), auto_commit );
    }
    else
    {
        if (auto_commit and this->needsAccepting())
            this->accept();

        int nffields = ffields_by_idx.count();
        FFPtr *ffields_array = ffields_by_idx.data();
        
        for (int i=0; i<nffields; ++i)
        {
            ffields_array[i].edit().update(molecules,auto_commit);
        }

        if (auto_commit and this->needsAccepting())
            this->accept();
    }
}

/** Update all of the forcefields in this group so that they have the 
    same version of the molecules that are present in the molecule
    group 'molgroup'
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void ForceFields::update(const MoleculeGroup &molgroup, bool auto_commit)
{
    this->update( molgroup.molecules(), auto_commit );
}

/** Set the contents of the molecule groups identified by the ID 'mgid'
    so that they only contain the view of the molecule in 'molview'.
    The passed property map is used to find any properties that are
    needed by the forcefields
    
    \throw SireMol::missing_group
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void ForceFields::setContents(const MGID &mgid, const MoleculeView &molview,
                              const PropertyMap &map)
{
    QList<MGNum> mgnums = mgid.map(*this);
    
    PartialMolecule view(molview);
    view.update( this->matchToExistingVersion(view.data()) );

    ForceFields old_state( *this );
    
    try
    {
        foreach (const MGNum &mgnum, mgnums)
        {
            this->_pvt_forceField(mgnum).setContents(mgnum, view, map);
        }
        
        this->rebuildIndex();
    }
    catch(...)
    {
        this->operator=(old_state);
        throw;
    }
}

/** Set the contents of the molecule groups identified by the ID 'mgid'
    so that they only contain the views of the molecule in 'molviews'.
    The passed property map is used to find any properties that are
    needed by the forcefields
    
    \throw SireMol::missing_group
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void ForceFields::setContents(const MGID &mgid, const ViewsOfMol &molviews,
                              const PropertyMap &map)
{
    QList<MGNum> mgnums = mgid.map(*this);
    
    ViewsOfMol views(molviews);
    views.update( this->matchToExistingVersion(views.data()) );

    ForceFields old_state( *this );
    
    try
    {
        foreach (const MGNum &mgnum, mgnums)
        {
            this->_pvt_forceField(mgnum).setContents(mgnum, views, map);
        }
        
        this->rebuildIndex();
    }
    catch(...)
    {
        this->operator=(old_state);
        throw;
    }
}

/** Set the contents of the molecule groups identified by the ID 'mgid'
    so that they only contain the molecules in 'molecules'.
    The passed property map is used to find any properties that are
    needed by the forcefields
    
    \throw SireMol::missing_group
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void ForceFields::setContents(const MGID &mgid, const Molecules &molecules,
                              const PropertyMap &map)
{
    QList<MGNum> mgnums = mgid.map(*this);

    Molecules mols = this->matchToExistingVersion(molecules);
    
    ForceFields old_state( *this );
    
    try
    {
        foreach (const MGNum &mgnum, mgnums)
        {
            this->_pvt_forceField(mgnum).setContents(mgnum, mols, map);
        }
        
        this->rebuildIndex();
    }
    catch(...)
    {
        this->operator=(old_state);
        throw;
    }
}

/** Set the contents of the molecule groups identified by the ID 'mgid'
    so that they only contain the molecules in the group 'molgroup'.
    The passed property map is used to find any properties that are
    needed by the forcefields
    
    \throw SireMol::missing_group
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void ForceFields::setContents(const MGID &mgid, const MoleculeGroup &molgroup,
                              const PropertyMap &map)
{
    QList<MGNum> mgnums = mgid.map(*this);
    
    MoleculeGroup group(molgroup);
    group.update( this->matchToExistingVersion(group.molecules()) );

    ForceFields old_state( *this );
    
    try
    {
        foreach (const MGNum &mgnum, mgnums)
        {
            this->_pvt_forceField(mgnum).setContents(mgnum, group, map);
        }
        
        this->rebuildIndex();
    }
    catch(...)
    {
        this->operator=(old_state);
        throw;
    }
}

/** Set the contents of the molecule groups identified by the ID 'mgid'
    so that they only contain the view of the molecule in 'molview'.
    Properties required by the forcefields are searched for in the
    default properties.
    
    \throw SireMol::missing_group
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void ForceFields::setContents(const MGID &mgid, const MoleculeView &molview)
{
    this->setContents(mgid, molview, PropertyMap());
}

/** Set the contents of the molecule groups identified by the ID 'mgid'
    so that they only contain the views of the molecule in 'molviews'.
    Properties required by the forcefields are searched for in the
    default properties.
    
    \throw SireMol::missing_group
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void ForceFields::setContents(const MGID &mgid, const ViewsOfMol &molviews)
{
    this->setContents(mgid, molviews, PropertyMap());
}

/** Set the contents of the molecule groups identified by the ID 'mgid'
    so that they only contain the molecules in 'molecules'.
    Properties required by the forcefields are searched for in the
    default properties.
    
    \throw SireMol::missing_group
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void ForceFields::setContents(const MGID &mgid, const Molecules &molecules)
{
    this->setContents(mgid, molecules, PropertyMap());
}

/** Set the contents of the molecule groups identified by the ID 'mgid'
    so that they only contain the molecules in the molecule group 'molgroup'.
    Properties required by the forcefields are searched for in the
    default properties.
    
    \throw SireMol::missing_group
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void ForceFields::setContents(const MGID &mgid, const MoleculeGroup &molgroup)
{
    this->setContents(mgid, molgroup, PropertyMap());
}

/** Return whether or not these forcefields are using any temporary workspace that needs
    to be accepted */
bool ForceFields::needsAccepting() const
{
    for (int i=0; i<ffields_by_idx.count(); ++i)
    {
        if (ffields_by_idx.constData()[i].read().needsAccepting())
            return true;
    }
    
    return false;
}

/** Ensure that any forcefields that are using temporary workspace have that accepted */
void ForceFields::accept()
{
    FFPtr *ffs = ffields_by_idx.data();

    for (int i=0; i<ffields_by_idx.count(); ++i)
    {
        if (ffs[i].read().needsAccepting())
        {
            ffs[i].edit().accept();
        }
    }
}

const char* ForceFields::typeName()
{
    return QMetaType::typeName( qMetaTypeId<ForceFields>() );
}
