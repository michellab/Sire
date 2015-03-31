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

#include "function.h"
#include "identities.h"
#include "symbols.h"
#include "functions.h"

#include "SireStream/datastream.h"

#include <boost/assert.hpp>

using namespace SireStream;
using namespace SireCAS;

////////////
//////////// Implementation of FunctionPvt
////////////

static const RegisterMetaType<FunctionPvt> r_functionpvt(MAGIC_ONLY, NO_ROOT,
                                                         "SireCAS::FunctionPvt");

/** Serialise to a binary datastream */
QDataStream& operator<<(QDataStream &ds, const FunctionPvt &f)
{
    writeHeader(ds, r_functionpvt, 1)
              << f.sig << f.syms << f.funcs;

    return ds;
}

/** Deserialise from a binary datastream */
QDataStream& operator>>(QDataStream &ds, FunctionPvt &f)
{
    VersionID v = readHeader(ds, r_functionpvt);

    if (v == 1)
    {
        ds >> f.sig >> f.syms >> f.funcs;
    }
    else
        throw version_error(v, "1", r_functionpvt, CODELOC);

    return ds;
}

/** Constructor */
FunctionPvt::FunctionPvt()
{}

/** Create a function called 'name' */
FunctionPvt::FunctionPvt(const QString &name) : sig(name)
{}

/** Copy constructor */
FunctionPvt::FunctionPvt(const FunctionPvt &other)
            : sig(other.sig), syms(other.syms), funcs(other.funcs)
{}

/** Destructor */
FunctionPvt::~FunctionPvt()
{}

/** Comparison operator */
bool FunctionPvt::operator==(const FunctionPvt &other) const
{
    return sig == other.sig and syms == other.syms and funcs == other.funcs;
}

/** Comparison operator */
bool FunctionPvt::operator!=(const FunctionPvt &other) const
{
    return not operator==(other);
}

/** Return a full string representation of the function */
QString FunctionPvt::toString() const
{
    if (sig.name().isNull() or sig.name().isEmpty() or syms.count() == 0)
        //this is a null function
        return QString::null;
    else
    {
        //create a string for the arguments
        QString args;

        int i=0;
        for (QHash<Symbol,int>::const_iterator it = syms.begin();
             it != syms.end();
             ++it)
        {
            if (i == 0)
                args = it.key().toString();
            else
                args = QString("%1,%2").arg(args, it.key().toString());

            if (it.value() > 0)
            {
                for (int j=0; j < it.value(); ++j)
                    args += "'";
            }
            else if (it.value() < 0)
            {
                for (int j=0; j < -(it.value()); ++j)
                    args += "`";
            }

            i++;
        }

        return QString("%1(%2)").arg(sig.name(),args);
    }
}

/** Add a symbol to the list of symbols in the function */
void FunctionPvt::add(const Symbol &symbol)
{
    if (symbol.ID()==0 or sig.contains(symbol.ID()))
        //we already have this symbol/function, or it is null
        return;
    else
    {
        if (symbol.isA<Function>())
        {
            const Function &func = symbol.asA<Function>();

            if (funcs.contains(func))
                return;

            //save that we are a function of this function
            funcs.insert(func);

            //we need to change our name from f( g(x) ) to f_g( x )
            sig.setName( QString("%1_%2").arg(sig.name(),func.name()) );

            //insert all of the symbols of this function
            for (QHash<Symbol,int>::const_iterator it = func.data().symbols().begin();
                 it != func.data().symbols().end();
                 ++it)
            {
                this->add( it.key() );
            }
        }
        else
        {
            //add this symbol to the list, with no differentiation or integration
            syms.insert(symbol, 0);
            sig.add(symbol.ID());
        }
    }
}

/** Increment the calculus count for the symbol 'symbol' */
void FunctionPvt::differentiate(const Symbol &symbol)
{
    if (syms.contains(symbol))
        syms[symbol] ++;
}

/** Decrement the calculus count for the symbol 'symbol' */
void FunctionPvt::integrate(const Symbol &symbol)
{
    if (syms.contains(symbol))
        syms[symbol] --;
}

//////////////
////////////// Implementation of Function
//////////////

static const RegisterMetaType<Function> r_function;

/** Serialise to a binary datastream */
QDataStream SIRECAS_EXPORT &operator<<(QDataStream &ds, const Function &f)
{
    ds << f.d << static_cast<const Symbol&>(f);

    return ds;
}

/** Deserialise from a binary datastream */
QDataStream SIRECAS_EXPORT &operator>>(QDataStream &ds, Function &f)
{
    ds >> f.d >> static_cast<Symbol&>(f);

    return ds;
}

/** Null constructor */
Function::Function() : Symbol()
{}

/** Construct a function called 'f' */
Function::Function(const QString &f) : Symbol(), d(f)
{}

/** Construct a function called 'f' that is a function of the passed
    symbols */
Function::Function(const QString &f, const Symbols &symbols)
         : Symbol(), d(f)
{
    *this = this->operator()(symbols);
}

/** Construct a function called 'f' that is a function of the passed
    symbols */
Function::Function(const QString &f,
                   const Symbol &sym0)
         : Symbol(), d(f)
{
    *this = this->operator()(sym0);
}

/** Construct a function called 'f' that is a function of the passed
    symbols */
Function::Function(const QString &f,
                   const Symbol &sym0, const Symbol &sym1)
         : Symbol(), d(f)
{
    *this = this->operator()(sym0, sym1);
}

/** Construct a function called 'f' that is a function of the passed
    symbols */
Function::Function(const QString &f,
                   const Symbol &sym0, const Symbol &sym1, const Symbol &sym2)
         : Symbol(), d(f)
{
    *this = this->operator()(sym0, sym1, sym2);
}

/** Construct a function called 'f' that is a function of the passed
    symbols */
Function::Function(const QString &f,
                   const Symbol &sym0, const Symbol &sym1, const Symbol &sym2,
                   const Symbol &sym3)
         : Symbol(), d(f)
{
    *this = this->operator()(sym0, sym1, sym2,
                             sym3);
}

/** Construct a function called 'f' that is a function of the passed
    symbols */
Function::Function(const QString &f,
                   const Symbol &sym0, const Symbol &sym1, const Symbol &sym2,
                   const Symbol &sym3, const Symbol &sym4)
         : Symbol(), d(f)
{
    *this = this->operator()(sym0, sym1, sym2,
                             sym3, sym4);
}

/** Construct a function called 'f' that is a function of the passed
    symbols */
Function::Function(const QString &f,
                   const Symbol &sym0, const Symbol &sym1, const Symbol &sym2,
                   const Symbol &sym3, const Symbol &sym4, const Symbol &sym5)
         : Symbol(), d(f)
{
    *this = this->operator()(sym0, sym1, sym2,
                             sym3, sym4, sym5);
}

/** Construct a function called 'f' that is a function of the passed
    symbols */
Function::Function(const QString &f,
                   const Symbol &sym0, const Symbol &sym1, const Symbol &sym2,
                   const Symbol &sym3, const Symbol &sym4, const Symbol &sym5,
                   const Symbol &sym6)
         : Symbol(), d(f)
{
    *this = this->operator()(sym0, sym1, sym2,
                             sym3, sym4, sym5,
                             sym6);
}

/** Construct a function called 'f' that is a function of the passed
    symbols */
Function::Function(const QString &f,
                   const Symbol &sym0, const Symbol &sym1, const Symbol &sym2,
                   const Symbol &sym3, const Symbol &sym4, const Symbol &sym5,
                   const Symbol &sym6, const Symbol &sym7)
         : Symbol(), d(f)
{
    *this = this->operator()(sym0, sym1, sym2,
                             sym3, sym4, sym5,
                             sym6, sym7);
}

/** Construct a function called 'f' that is a function of the passed
    symbols */
Function::Function(const QString &f,
                   const Symbol &sym0, const Symbol &sym1, const Symbol &sym2,
                   const Symbol &sym3, const Symbol &sym4, const Symbol &sym5,
                   const Symbol &sym6, const Symbol &sym7, const Symbol &sym8)
         : Symbol(), d(f)
{
    *this = this->operator()(sym0, sym1, sym2,
                             sym3, sym4, sym5,
                             sym6, sym7, sym8);
}

/** Construct a function called 'f' that is a function of the passed
    symbols */
Function::Function(const QString &f,
                   const Symbol &sym0, const Symbol &sym1, const Symbol &sym2,
                   const Symbol &sym3, const Symbol &sym4, const Symbol &sym5,
                   const Symbol &sym6, const Symbol &sym7, const Symbol &sym8,
                   const Symbol &sym9)
         : Symbol(), d(f)
{
    *this = this->operator()(sym0, sym1, sym2,
                             sym3, sym4, sym5,
                             sym6, sym7, sym8,
                             sym9);
}

/** Construct a new function from the data contained in the FunctionPvt object */
Function::Function(const FunctionPvt &data)
         : Symbol(data.toString()), d(data)
{}

/** Copy constructor */
Function::Function(const Function &other)
         : Symbol(other), d(other.d)
{}

/** Destructor */
Function::~Function()
{}

/** Comparison operator */
bool Function::operator==(const ExBase &other) const
{
    const Symbol *other_func = dynamic_cast<const Symbol*>(&other);

    return other_func != 0 and ID() == other_func->ID();
}

/** Assignment operator */
Function& Function::operator=(const Function &other)
{
    Symbol::operator=(other);
    d = other.d;

    return *this;
}

/** Return whether or not this is a function of 'symbol' */
bool Function::isFunction(const Symbol &symbol) const
{
    if ( Symbol::isFunction(symbol) or d.symbols().contains(symbol) )
        return true;
    else if (symbol.isA<Function>())
    {
        //check that we are not a function of this function
        for (QSet<Function>::const_iterator it = d.functions().begin();
             it != d.functions().end();
             ++it)
        {
            if (it->isFunction(symbol))
                return true;
        }
    }

    return false;
}

/** Return the differential of this function with respect to 'symbol' */
Expression Function::differentiate(const Symbol &symbol) const
{
    BOOST_ASSERT( not symbol.isA<Function>() );

    if (this->isFunction(symbol))
    {
        FunctionPvt newdata(d);
        newdata.differentiate(symbol);
        return Function(newdata);
    }
    else
        return 0;
}

/** Return the integral of this function with respect to 'symbol' */
Expression Function::integrate(const Symbol &symbol) const
{
    BOOST_ASSERT( not symbol.isA<Function>() );

    if (this->isFunction(symbol))
    {
        FunctionPvt newdata(d);
        newdata.integrate(symbol);
        return Function(newdata);
    }
    else
        return *this * symbol;
}

/** Return a hash of this Function */
uint Function::hash() const
{
    return Symbol::hash();
}

/** Return the expression that matches this function */
Expression Function::match(const Function &func, const Expression &ex) const
{
    if (data() == func.data())
    {
        //the function is the same as 'func', so nothing to do :-)
        return ex;
    }
    else
    {
        FunctionPvt fdata = func.data();

        //go through each symbol in turn...
        for (QHash<Symbol,int>::const_iterator it = fdata.symbols().begin();
             it != fdata.symbols().end();
             ++it)
        {
            const Symbol &symbol = it.key();

            //compare the calculus level of the symbol in 'func' to this
            //function...
            BOOST_ASSERT( data().symbols().contains(symbol) );

            int current_level = data().symbols().value(symbol);

            if ( it.value() > current_level )
            {
                //we need to integrate 'func' and 'ex' with
                //respect to symbol to match this function
                fdata.integrate(symbol);
                Expression integ = ex.integrate(symbol);

                //leave the integration constant on the expression - this will
                //cause problems if we need to integrate again...
                return this->match( Function(fdata), integ );
            }
            else if ( it.value() < current_level )
            {
                //we need to differentiate 'func' and 'ex' with
                //respect to symbol to match this function

                fdata.differentiate(symbol);
                Expression diff = ex.differentiate(symbol);

                return this->match( Function(fdata), diff );
            }
        }

        //all of the symbols are at the same level, so surely data() == func.data()?
        return ex;
    }
}

/** Substitute this function with an identity (it one exists) */
Expression Function::substitute(const Identities &identities) const
{
    if (not identities.contains(*this))
    {
        return Expression(*this);
    }

    //get the form of this function that is stored in the identities
    Function stored_func = identities.function(*this);

    //get the expression that this stored function is equal to...
    Expression ex = identities.expression(stored_func);

    //now differentiate/integrate this expression as required so
    //that it equals this function...
    return this->match(stored_func, ex);
}

/** Return all of the symbols in the function */
Symbols Function::symbols() const
{
    //need to build the set of symbols by hand as we have no operator= in ExBase
    Symbols syms;

    for (QHash<Symbol,int>::const_iterator it = d.symbols().begin();
         it != d.symbols().end();
         ++it)
    {
        syms.insert( it.key() );
    }

    return syms;
}

/** Return this function */
Functions Function::functions() const
{
    return Functions(*this);
}

/** Return this Function as a function of the passed symbols */
Function Function::operator()(const Symbols &symbols) const
{
    if (symbols.isEmpty())
        return *this;

    FunctionPvt newdata(d);

    for (Symbols::const_iterator it = symbols.begin();
         it != symbols.end();
         ++it)
    {
        newdata.add(*it);
    }

    return Function(newdata);
}

/** Return this Function as a function of the passed symbols */
Function Function::operator()(const Symbol &sym0) const
{
    FunctionPvt newdata(d);

    newdata.add(sym0);

    return Function(newdata);
}

/** Return this Function as a function of the passed symbols */
Function Function::operator()(const Symbol &sym0, const Symbol &sym1) const
{
    FunctionPvt newdata(d);

    newdata.add(sym0);
    newdata.add(sym1);

    return Function(newdata);
}

/** Return this Function as a function of the passed symbols */
Function Function::operator()(const Symbol &sym0, const Symbol &sym1, const Symbol &sym2) const
{
    FunctionPvt newdata(d);

    newdata.add(sym0);
    newdata.add(sym1);
    newdata.add(sym2);

    return Function(newdata);
}

/** Return this Function as a function of the passed symbols */
Function Function::operator()(const Symbol &sym0, const Symbol &sym1, const Symbol &sym2,
                              const Symbol &sym3) const
{
    FunctionPvt newdata(d);

    newdata.add(sym0);
    newdata.add(sym1);
    newdata.add(sym2);
    newdata.add(sym3);

    return Function(newdata);
}

/** Return this Function as a function of the passed symbols */
Function Function::operator()(const Symbol &sym0, const Symbol &sym1, const Symbol &sym2,
                              const Symbol &sym3, const Symbol &sym4) const
{
    FunctionPvt newdata(d);

    newdata.add(sym0);
    newdata.add(sym1);
    newdata.add(sym2);
    newdata.add(sym3);
    newdata.add(sym4);

    return Function(newdata);
}

/** Return this Function as a function of the passed symbols */
Function Function::operator()(const Symbol &sym0, const Symbol &sym1, const Symbol &sym2,
                              const Symbol &sym3, const Symbol &sym4, const Symbol &sym5) const
{
    FunctionPvt newdata(d);

    newdata.add(sym0);
    newdata.add(sym1);
    newdata.add(sym2);
    newdata.add(sym3);
    newdata.add(sym4);
    newdata.add(sym5);

    return Function(newdata);
}

/** Return this Function as a function of the passed symbols */
Function Function::operator()(const Symbol &sym0, const Symbol &sym1, const Symbol &sym2,
                              const Symbol &sym3, const Symbol &sym4, const Symbol &sym5,
                              const Symbol &sym6) const
{
    FunctionPvt newdata(d);

    newdata.add(sym0);
    newdata.add(sym1);
    newdata.add(sym2);
    newdata.add(sym3);
    newdata.add(sym4);
    newdata.add(sym5);
    newdata.add(sym6);

    return Function(newdata);
}

/** Return this Function as a function of the passed symbols */
Function Function::operator()(const Symbol &sym0, const Symbol &sym1, const Symbol &sym2,
                              const Symbol &sym3, const Symbol &sym4, const Symbol &sym5,
                              const Symbol &sym6, const Symbol &sym7) const
{
    FunctionPvt newdata(d);

    newdata.add(sym0);
    newdata.add(sym1);
    newdata.add(sym2);
    newdata.add(sym3);
    newdata.add(sym4);
    newdata.add(sym5);
    newdata.add(sym6);
    newdata.add(sym7);

    return Function(newdata);
}

/** Return this Function as a function of the passed symbols */
Function Function::operator()(const Symbol &sym0, const Symbol &sym1, const Symbol &sym2,
                              const Symbol &sym3, const Symbol &sym4, const Symbol &sym5,
                              const Symbol &sym6, const Symbol &sym7, const Symbol &sym8) const
{
    FunctionPvt newdata(d);

    newdata.add(sym0);
    newdata.add(sym1);
    newdata.add(sym2);
    newdata.add(sym3);
    newdata.add(sym4);
    newdata.add(sym5);
    newdata.add(sym6);
    newdata.add(sym7);
    newdata.add(sym8);

    return Function(newdata);
}

/** Return this Function as a function of the passed symbols */
Function Function::operator()(const Symbol &sym0, const Symbol &sym1, const Symbol &sym2,
                              const Symbol &sym3, const Symbol &sym4, const Symbol &sym5,
                              const Symbol &sym6, const Symbol &sym7, const Symbol &sym8,
                              const Symbol &sym9) const
{
    FunctionPvt newdata(d);

    newdata.add(sym0);
    newdata.add(sym1);
    newdata.add(sym2);
    newdata.add(sym3);
    newdata.add(sym4);
    newdata.add(sym5);
    newdata.add(sym6);
    newdata.add(sym7);
    newdata.add(sym8);
    newdata.add(sym9);

    return Function(newdata);
}

const char* Function::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Function>() );
}

Function* Function::clone() const
{
    return new Function(*this);
}

