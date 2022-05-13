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

#include "internalperturbation.h"

#include "twoatomfunctions.h"
#include "threeatomfunctions.h"
#include "fouratomfunctions.h"

#include "SireMol/molecule.h"
#include "SireMol/moleditor.h"
#include "SireMol/mover.hpp"
#include "SireMol/core.h"

#include "SireCAS/values.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMM;
using namespace SireMol;
using namespace SireCAS;
using namespace SireBase;
using namespace SireStream;

////////////
//////////// Implementation of InternalPerturbation
////////////

static const RegisterMetaType<InternalPerturbation> r_intpert( MAGIC_ONLY,
                                                    InternalPerturbation::typeName() );

QDataStream &operator<<(QDataStream &ds,
                                      const InternalPerturbation &intpert)
{
    writeHeader(ds, r_intpert, 1);

    SharedDataStream sds(ds);

    sds << intpert.base_expression << intpert.initial_forms
        << intpert.final_forms << static_cast<const Perturbation&>(intpert);

    return ds;
}

QDataStream &operator>>(QDataStream &ds,
                                      InternalPerturbation &intpert)
{
    VersionID v = readHeader(ds, r_intpert);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> intpert.base_expression >> intpert.initial_forms
            >> intpert.final_forms >> static_cast<Perturbation&>(intpert);

        intpert.buildPerturbExpression();
    }
    else
        throw version_error(v, "1", r_intpert, CODELOC);

    return ds;
}

/** Internal function used to rebuild the perturb expression from
    the mapping function, base expression and initial and final forms */
void InternalPerturbation::buildPerturbExpression()
{
    //build all of the perturbed identities
    Identities perturbed_idents;

    const Expression &mapfunc = mappingFunction();

    foreach (Symbol symbol, initial_forms.symbols())
    {
        if (final_forms.contains(symbol))
        {
            Identities ident;
            ident.set( symbols().initial(), initial_forms[symbol] );
            ident.set( symbols().final(), final_forms[symbol] );

            perturbed_idents.set(symbol, mapfunc.substitute(ident));
        }
        else
        {
            Identities ident;
            ident.set( symbols().initial(), initial_forms[symbol] );
            ident.set( symbols().final(), Expression(0) );

            perturbed_idents.set(symbol, mapfunc.substitute(ident));
        }
    }

    foreach (Symbol symbol, final_forms.symbols())
    {
        if (not initial_forms.contains(symbol))
        {
            Identities ident;
            ident.set( symbols().initial(), Expression(0) );
            ident.set( symbols().final(), final_forms[symbol] );

            perturbed_idents.set(symbol, mapfunc.substitute(ident));
        }
    }

    perturb_expression = base_expression.substitute(perturbed_idents);
}

/** Null constructor */
InternalPerturbation::InternalPerturbation() : Perturbation()
{}

/** Construct to map from one function to another (from initial_function to
    final_function) using the passed mapping function */
InternalPerturbation::InternalPerturbation(const Expression &initial_function,
                                           const Expression &final_function,
                                           const Expression &mapping_function,
                                           const PropertyMap &map)
                     : Perturbation(mapping_function, map)
{
    Symbol f("f");

    initial_forms = (f == initial_function);
    final_forms = (f == final_function);
    base_expression = f;

    this->buildPerturbExpression();
}

/** Construct to map from 'base_expression' populated with the identities
    in 'initial_forms' to 'base_expression' populated with the identities
    in 'final_forms', using the passed mapping function to map the identites
    from initial to final */
InternalPerturbation::InternalPerturbation(const Expression &base,
                                           const Identities &initial,
                                           const Identities &final,
                                           const Expression &mapping_function,
                                           const PropertyMap &map)
                     : Perturbation(mapping_function, map),
                       base_expression(base),
                       initial_forms(initial), final_forms(final)
{
    this->buildPerturbExpression();
}

/** Copy constructor */
InternalPerturbation::InternalPerturbation(const InternalPerturbation &other)
                     : Perturbation(other),
                       base_expression(other.base_expression),
                       perturb_expression(other.perturb_expression),
                       initial_forms(other.initial_forms),
                       final_forms(other.final_forms)
{}

/** Destructor */
InternalPerturbation::~InternalPerturbation()
{}

/** Copy assignment operator */
InternalPerturbation& InternalPerturbation::operator=(const InternalPerturbation &other)
{
    if (this != &other)
    {
        Perturbation::operator=(other);
        base_expression = other.base_expression;
        perturb_expression = other.perturb_expression;
        initial_forms = other.initial_forms;
        final_forms = other.final_forms;
    }

    return *this;
}

/** Comparison operator */
bool InternalPerturbation::operator==(const InternalPerturbation &other) const
{
    return base_expression == other.base_expression and
           initial_forms == other.initial_forms and
           final_forms == other.final_forms and
           Perturbation::operator==(other);
}

/** Comparison operator */
bool InternalPerturbation::operator!=(const InternalPerturbation &other) const
{
    return not InternalPerturbation::operator==(other);
}

const char* InternalPerturbation::typeName()
{
    return "SireMM::InternalPerturbation";
}

/** Return the base expression - this is the expression into which
    the mapped identites are substituted */
const Expression& InternalPerturbation::baseExpression() const
{
    return base_expression;
}

/** Return the perturbed expression - this is the expression that
    is used to calculate the energy */
const Expression& InternalPerturbation::perturbExpression() const
{
    return perturb_expression;
}

/** Return the initial forms - these are the identities that
    are substituted into the base expression at the initial state */
const Identities& InternalPerturbation::initialForms() const
{
    return initial_forms;
}

/** Return the final forms - these are the identities that
    are substituted into the base expression at the final state */
const Identities& InternalPerturbation::finalForms() const
{
    return final_forms;
}

PerturbationPtr InternalPerturbation::recreate(const Expression &expression) const
{
    PerturbationPtr ret = Perturbation::recreate(expression);

    ret.edit().asA<InternalPerturbation>().buildPerturbExpression();

    return ret;
}

PerturbationPtr InternalPerturbation::recreate(const Expression &expression,
                                               const PropertyMap &map) const
{
    PerturbationPtr ret = Perturbation::recreate(expression, map);

    ret.edit().asA<InternalPerturbation>().buildPerturbExpression();

    return ret;
}

PerturbationPtr InternalPerturbation::substitute(const Identities &identities) const
{
    PerturbationPtr ret = Perturbation::substitute(identities);

    ret.edit().asA<InternalPerturbation>().buildPerturbExpression();

    return ret;
}

////////////
//////////// Implementation of TwoAtomPerturbation
////////////

static const RegisterMetaType<TwoAtomPerturbation> r_twoatom;

QDataStream &operator<<(QDataStream &ds,
                                      const TwoAtomPerturbation &twoatom)
{
    writeHeader(ds, r_twoatom, 1);

    SharedDataStream sds(ds);

    sds << twoatom.atm0 << twoatom.atm1
        << static_cast<const InternalPerturbation&>(twoatom);

    return ds;
}

QDataStream &operator>>(QDataStream &ds,
                                      TwoAtomPerturbation &twoatom)
{
    VersionID v = readHeader(ds, r_twoatom);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> twoatom.atm0 >> twoatom.atm1
            >> static_cast<InternalPerturbation&>(twoatom);
    }
    else
        throw version_error(v, "1", r_twoatom, CODELOC);

    return ds;
}

/** Null constructor */
TwoAtomPerturbation::TwoAtomPerturbation()
                    : ConcreteProperty<TwoAtomPerturbation,InternalPerturbation>()
{}

/** Construct to perturb the function between the atoms 'atom0' and 'atom1'
    to use 'initial_form' at the initial state and 'final_form' at the
    final state, where the functions are mapped between these two states
    using the default mapping function */
TwoAtomPerturbation::TwoAtomPerturbation(const AtomID &atom0, const AtomID &atom1,
                                         const Expression &initial_form,
                                         const Expression &final_form,
                                         const PropertyMap &map)
     : ConcreteProperty<TwoAtomPerturbation,InternalPerturbation>(initial_form,
                                                                  final_form,
                                                    Perturbation::defaultFunction(),
                                                                  map),
       atm0(atom0), atm1(atom1)
{}

/** Construct to perturb the function between the atoms 'atom0' and 'atom1'
    to use 'initial_form' at the initial state and 'final_form' at the
    final state, where the functions are mapped between these two states
    using 'mapping_function' */
TwoAtomPerturbation::TwoAtomPerturbation(const AtomID &atom0, const AtomID &atom1,
                                         const Expression &initial_form,
                                         const Expression &final_form,
                                         const Expression &mapping_function,
                                         const PropertyMap &map)
     : ConcreteProperty<TwoAtomPerturbation,InternalPerturbation>(initial_form,
                                                                  final_form,
                                                                  mapping_function,
                                                                  map),
       atm0(atom0), atm1(atom1)
{}

/** Construct to perturb the function between the atoms 'atom0' and 'atom1'
    to use 'base_expression' populated with the identities in 'initial_forms'
    at the initial state, and populated with the identities in 'final_forms'
    at the final state, where the identities are mapped between the initial
    and final states using the default mapping function */
TwoAtomPerturbation::TwoAtomPerturbation(const AtomID &atom0, const AtomID &atom1,
                                         const Expression &base_expression,
                                         const Identities &initial_forms,
                                         const Identities &final_forms,
                                         const PropertyMap &map)
     : ConcreteProperty<TwoAtomPerturbation,InternalPerturbation>(base_expression,
                                                                  initial_forms,
                                                                  final_forms,
                                                    Perturbation::defaultFunction(),
                                                                  map),
       atm0(atom0), atm1(atom1)
{}


/** Construct to perturb the function between the atoms 'atom0' and 'atom1'
    to use 'base_expression' populated with the identities in 'initial_forms'
    at the initial state, and populated with the identities in 'final_forms'
    at the final state, where the identities are mapped between the initial
    and final states using 'mapping_function' */
TwoAtomPerturbation::TwoAtomPerturbation(const AtomID &atom0, const AtomID &atom1,
                                         const Expression &base_expression,
                                         const Identities &initial_forms,
                                         const Identities &final_forms,
                                         const Expression &mapping_function,
                                         const PropertyMap &map)
     : ConcreteProperty<TwoAtomPerturbation,InternalPerturbation>(base_expression,
                                                                  initial_forms,
                                                                  final_forms,
                                                                  mapping_function,
                                                                  map),
       atm0(atom0), atm1(atom1)
{}

/** Copy constructor */
TwoAtomPerturbation::TwoAtomPerturbation(const TwoAtomPerturbation &other)
                    : ConcreteProperty<TwoAtomPerturbation,InternalPerturbation>(other),
                      atm0(other.atm0), atm1(other.atm1)
{}

/** Destructor */
TwoAtomPerturbation::~TwoAtomPerturbation()
{}

/** Copy assignment operator */
TwoAtomPerturbation& TwoAtomPerturbation::operator=(const TwoAtomPerturbation &other)
{
    if (this != &other)
    {
        atm0 = other.atm0;
        atm1 = other.atm1;
        InternalPerturbation::operator=(other);
    }

    return *this;
}

/** Comparison operator */
bool TwoAtomPerturbation::operator==(const TwoAtomPerturbation &other) const
{
    return ( (atm0 == other.atm0 and atm1 == other.atm1) or
             (atm0 == other.atm1 and atm1 == other.atm0) ) and
           InternalPerturbation::operator==(other);
}

/** Comparison operator */
bool TwoAtomPerturbation::operator!=(const TwoAtomPerturbation &other) const
{
    return not TwoAtomPerturbation::operator==(other);
}

const char* TwoAtomPerturbation::typeName()
{
    return QMetaType::typeName( qMetaTypeId<TwoAtomPerturbation>() );
}

/** Return the first of the two atoms whose potential is being changed */
const AtomID& TwoAtomPerturbation::atom0() const
{
    return atm0.base();
}

/** Return the second of the two atoms whose potential is being changed */
const AtomID& TwoAtomPerturbation::atom1() const
{
    return atm1.base();
}

/** Return a string representation of this perturbation */
QString TwoAtomPerturbation::toString() const
{
    return QObject::tr("TwoAtomPerturbation( %1-%2 : %3 )")
            .arg(atm0.toString(), atm1.toString())
            .arg(perturbExpression().toString());
}

/** Return the properties required or changed by this perturbation */
QSet<QString> TwoAtomPerturbation::requiredProperties() const
{
    QSet<QString> props;

    const PropertyName &param_property = propertyMap()["parameters"];

    if (param_property.hasSource())
        props.insert(param_property.source());

    return props;
}

/** Return whether or not this perturbation with the passed values would
    change the molecule 'molecule' */
bool TwoAtomPerturbation::wouldChange(const Molecule &molecule,
                                      const Values &values) const
{
    const PropertyName &param_property = propertyMap()["parameters"];

    if (not molecule.hasProperty(param_property))
        return true;

    const TwoAtomFunctions &funcs = molecule.property(param_property)
                                            .asA<TwoAtomFunctions>();

    Identities idents;

    foreach (Symbol symbol, values.symbols())
    {
        idents.set(symbol, Expression(values[symbol]));
    }

    Expression new_function = perturbExpression().substitute(idents);

    return new_function != funcs.potential(atm0, atm1);
}

/** Perturb the two atom function in passed molecule using the reaction
    coordinate(s) in 'values'

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void TwoAtomPerturbation::perturbMolecule(MolEditor &molecule,
                                          const Values &values) const
{
    const PropertyName &param_property = propertyMap()["parameters"];

    TwoAtomFunctions funcs;

    bool new_property = false;

    if (molecule.hasProperty(param_property))
    {
        funcs = molecule.property(param_property).asA<TwoAtomFunctions>();
    }
    else
    {
        funcs = TwoAtomFunctions(molecule.data());
        new_property = true;
    }

    //calculate the perturbed function
    Identities idents;

    foreach (Symbol symbol, values.symbols())
    {
        idents.set(symbol, Expression(values[symbol]));
    }

    Expression new_function = perturbExpression().substitute(idents);

    Expression old_function = funcs.potential(atm0, atm1);

    if (new_property or (new_function != old_function))
    {
      //funcs.set(atm0, atm1, new_function);
        funcs.set( BondID(atm0, atm1), new_function);
        molecule.setProperty(param_property, funcs);
    }
}

////////////
//////////// Implementation of ThreeAtomPerturbation
////////////

static const RegisterMetaType<ThreeAtomPerturbation> r_threeatom;

QDataStream &operator<<(QDataStream &ds,
                                      const ThreeAtomPerturbation &threeatom)
{
    writeHeader(ds, r_threeatom, 1);

    SharedDataStream sds(ds);

    sds << threeatom.atm0 << threeatom.atm1 << threeatom.atm2
        << static_cast<const InternalPerturbation&>(threeatom);

    return ds;
}

QDataStream &operator>>(QDataStream &ds,
                                      ThreeAtomPerturbation &threeatom)
{
    VersionID v = readHeader(ds, r_threeatom);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> threeatom.atm0 >> threeatom.atm1 >> threeatom.atm2
            >> static_cast<InternalPerturbation&>(threeatom);
    }
    else
        throw version_error(v, "1", r_threeatom, CODELOC);

    return ds;
}

/** Null constructor */
ThreeAtomPerturbation::ThreeAtomPerturbation()
                      : ConcreteProperty<ThreeAtomPerturbation,InternalPerturbation>()
{}

/** Construct to perturb the function between the atoms 'atom0', 'atom1' and 'atom2'
    to use 'initial_form' at the initial state and 'final_form' at the
    final state, where the functions are mapped between these two states
    using the default mapping function */
ThreeAtomPerturbation::ThreeAtomPerturbation(const AtomID &atom0, const AtomID &atom1,
                                             const AtomID &atom2,
                                             const Expression &initial_form,
                                             const Expression &final_form,
                                             const PropertyMap &map)
     : ConcreteProperty<ThreeAtomPerturbation,InternalPerturbation>(initial_form,
                                                                    final_form,
                                                      Perturbation::defaultFunction(),
                                                                    map),
       atm0(atom0), atm1(atom1), atm2(atom2)
{}

/** Construct to perturb the function between the atoms 'atom0', 'atom1' and 'atom2'
    to use 'initial_form' at the initial state and 'final_form' at the
    final state, where the functions are mapped between these two states
    using 'mapping_function' */
ThreeAtomPerturbation::ThreeAtomPerturbation(const AtomID &atom0, const AtomID &atom1,
                                             const AtomID &atom2,
                                             const Expression &initial_form,
                                             const Expression &final_form,
                                             const Expression &mapping_function,
                                             const PropertyMap &map)
     : ConcreteProperty<ThreeAtomPerturbation,InternalPerturbation>(initial_form,
                                                                    final_form,
                                                                    mapping_function,
                                                                    map),
       atm0(atom0), atm1(atom1), atm2(atom2)
{}

/** Construct to perturb the function between the atoms 'atom0', 'atom1' and 'atom2'
    to use 'base_expression' populated with the identities in 'initial_forms'
    at the initial state, and populated with the identities in 'final_forms'
    at the final state, where the identities are mapped between the initial
    and final states using the default mapping function */
ThreeAtomPerturbation::ThreeAtomPerturbation(const AtomID &atom0, const AtomID &atom1,
                                             const AtomID &atom2,
                                             const Expression &base_expression,
                                             const Identities &initial_forms,
                                             const Identities &final_forms,
                                             const PropertyMap &map)
     : ConcreteProperty<ThreeAtomPerturbation,InternalPerturbation>(base_expression,
                                                                    initial_forms,
                                                                    final_forms,
                                                      Perturbation::defaultFunction(),
                                                                    map),
       atm0(atom0), atm1(atom1), atm2(atom2)
{}

/** Construct to perturb the function between the atoms 'atom0', 'atom1' and 'atom2'
    to use 'base_expression' populated with the identities in 'initial_forms'
    at the initial state, and populated with the identities in 'final_forms'
    at the final state, where the identities are mapped between the initial
    and final states using 'mapping_function' */
ThreeAtomPerturbation::ThreeAtomPerturbation(const AtomID &atom0, const AtomID &atom1,
                                             const AtomID &atom2,
                                             const Expression &base_expression,
                                             const Identities &initial_forms,
                                             const Identities &final_forms,
                                             const Expression &mapping_function,
                                             const PropertyMap &map)
     : ConcreteProperty<ThreeAtomPerturbation,InternalPerturbation>(base_expression,
                                                                    initial_forms,
                                                                    final_forms,
                                                                    mapping_function,
                                                                    map),
       atm0(atom0), atm1(atom1), atm2(atom2)
{}

/** Copy constructor */
ThreeAtomPerturbation::ThreeAtomPerturbation(const ThreeAtomPerturbation &other)
             : ConcreteProperty<ThreeAtomPerturbation,InternalPerturbation>(other),
               atm0(other.atm0), atm1(other.atm1), atm2(other.atm2)
{}

/** Destructor */
ThreeAtomPerturbation::~ThreeAtomPerturbation()
{}

/** Copy assignment operator */
ThreeAtomPerturbation& ThreeAtomPerturbation::operator=(const ThreeAtomPerturbation &other)
{
    if (this != &other)
    {
        atm0 = other.atm0;
        atm1 = other.atm1;
        atm2 = other.atm2;
        InternalPerturbation::operator=(other);
    }

    return *this;
}

/** Comparison operator */
bool ThreeAtomPerturbation::operator==(const ThreeAtomPerturbation &other) const
{
    return ( (atm0 == other.atm0 and atm2 == other.atm2) or
             (atm0 == other.atm2 and atm2 == other.atm0) ) and
             atm1 == other.atm1 and
           InternalPerturbation::operator==(other);
}

/** Comparison operator */
bool ThreeAtomPerturbation::operator!=(const ThreeAtomPerturbation &other) const
{
    return not ThreeAtomPerturbation::operator==(other);
}

const char* ThreeAtomPerturbation::typeName()
{
    return QMetaType::typeName( qMetaTypeId<ThreeAtomPerturbation>() );
}

/** Return the first of the three atoms whose potential is being changed */
const AtomID& ThreeAtomPerturbation::atom0() const
{
    return atm0.base();
}

/** Return the second of the three atoms whose potential is being changed */
const AtomID& ThreeAtomPerturbation::atom1() const
{
    return atm1.base();
}

/** Return the third of the three atoms whose potential is being changed */
const AtomID& ThreeAtomPerturbation::atom2() const
{
    return atm2.base();
}

/** Return a string representation of this perturbation */
QString ThreeAtomPerturbation::toString() const
{
    return QObject::tr("ThreeAtomPerturbation( %1-%2-%3 : %4 )")
            .arg(atm0.toString(), atm1.toString(), atm2.toString())
            .arg(perturbExpression().toString());
}

/** Return the properties required or changed by this perturbation */
QSet<QString> ThreeAtomPerturbation::requiredProperties() const
{
    QSet<QString> props;

    const PropertyName &param_property = propertyMap()["parameters"];

    if (param_property.hasSource())
        props.insert(param_property.source());

    return props;
}

/** Return whether or not this perturbation with the passed values would
    change the molecule 'molecule' */
bool ThreeAtomPerturbation::wouldChange(const Molecule &molecule,
                                        const Values &values) const
{
    const PropertyName &param_property = propertyMap()["parameters"];

    if (not molecule.hasProperty(param_property))
        return true;

    const ThreeAtomFunctions &funcs = molecule.property(param_property)
                                              .asA<ThreeAtomFunctions>();

    Identities idents;

    foreach (Symbol symbol, values.symbols())
    {
        idents.set(symbol, Expression(values[symbol]));
    }

    Expression new_function = perturbExpression().substitute(idents);

    return new_function != funcs.potential(atm0, atm1, atm2);
}

/** Perturb the three atom function in passed molecule using the reaction
    coordinate(s) in 'values'

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void ThreeAtomPerturbation::perturbMolecule(MolEditor &molecule,
                                            const Values &values) const
{
    const PropertyName &param_property = propertyMap()["parameters"];

    ThreeAtomFunctions funcs;

    bool new_property = false;

    if (molecule.hasProperty(param_property))
    {
        funcs = molecule.property(param_property).asA<ThreeAtomFunctions>();
    }
    else
    {
        funcs = ThreeAtomFunctions(molecule.data());
        new_property = true;
    }

    //calculate the perturbed function
    Identities idents;

    foreach (Symbol symbol, values.symbols())
    {
        idents.set(symbol, Expression(values[symbol]));
    }

    Expression new_function = perturbExpression().substitute(idents);

    Expression old_function = funcs.potential(atm0, atm1, atm2);

    if (new_property or (new_function != old_function))
    {
      //funcs.set(atm0, atm1, atm2, new_function);
        funcs.set(AngleID(atm0, atm1, atm2), new_function);
        molecule.setProperty(param_property, funcs);
    }
}

////////////
//////////// Implementation of FourAtomPerturbation
////////////

static const RegisterMetaType<FourAtomPerturbation> r_fouratom;

QDataStream &operator<<(QDataStream &ds,
                                      const FourAtomPerturbation &fouratom)
{
    writeHeader(ds, r_fouratom, 1);

    SharedDataStream sds(ds);

    sds << fouratom.atm0 << fouratom.atm1 << fouratom.atm2 << fouratom.atm3
        << static_cast<const InternalPerturbation&>(fouratom);

    return ds;
}

QDataStream &operator>>(QDataStream &ds,
                                      FourAtomPerturbation &fouratom)
{
    VersionID v = readHeader(ds, r_fouratom);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> fouratom.atm0 >> fouratom.atm1 >> fouratom.atm2 >> fouratom.atm3
            >> static_cast<InternalPerturbation&>(fouratom);
    }
    else
        throw version_error(v, "1", r_fouratom, CODELOC);

    return ds;
}

/** Null constructor */
FourAtomPerturbation::FourAtomPerturbation()
                     : ConcreteProperty<FourAtomPerturbation,InternalPerturbation>()
{}

/** Construct to perturb the function between the atoms
    to use 'initial_form' at the initial state and 'final_form' at the
    final state, where the functions are mapped between these two states
    using the default mapping function */
FourAtomPerturbation::FourAtomPerturbation(const AtomID &atom0, const AtomID &atom1,
                                           const AtomID &atom2, const AtomID &atom3,
                                           const Expression &initial_form,
                                           const Expression &final_form,
                                           const PropertyMap &map)
     : ConcreteProperty<FourAtomPerturbation,InternalPerturbation>(initial_form,
                                                                   final_form,
                                                     Perturbation::defaultFunction(),
                                                                   map),
       atm0(atom0), atm1(atom1), atm2(atom2), atm3(atom3)
{}

/** Construct to perturb the function between the atoms
    to use 'initial_form' at the initial state and 'final_form' at the
    final state, where the functions are mapped between these two states
    using 'mapping_function' */
FourAtomPerturbation::FourAtomPerturbation(const AtomID &atom0, const AtomID &atom1,
                                           const AtomID &atom2, const AtomID &atom3,
                                           const Expression &initial_form,
                                           const Expression &final_form,
                                           const Expression &mapping_function,
                                           const PropertyMap &map)
     : ConcreteProperty<FourAtomPerturbation,InternalPerturbation>(initial_form,
                                                                   final_form,
                                                                   mapping_function,
                                                                   map),
       atm0(atom0), atm1(atom1), atm2(atom2), atm3(atom3)
{}

/** Construct to perturb the function between the atoms
    to use 'base_expression' populated with the identities in 'initial_forms'
    at the initial state, and populated with the identities in 'final_forms'
    at the final state, where the identities are mapped between the initial
    and final states using the default mapping function */
FourAtomPerturbation::FourAtomPerturbation(const AtomID &atom0, const AtomID &atom1,
                                           const AtomID &atom2, const AtomID &atom3,
                                           const Expression &base_expression,
                                           const Identities &initial_forms,
                                           const Identities &final_forms,
                                           const PropertyMap &map)
     : ConcreteProperty<FourAtomPerturbation,InternalPerturbation>(base_expression,
                                                                   initial_forms,
                                                                   final_forms,
                                                     Perturbation::defaultFunction(),
                                                                   map),
       atm0(atom0), atm1(atom1), atm2(atom2), atm3(atom3)
{}

/** Construct to perturb the function between the atoms
    to use 'base_expression' populated with the identities in 'initial_forms'
    at the initial state, and populated with the identities in 'final_forms'
    at the final state, where the identities are mapped between the initial
    and final states using 'mapping_function' */
FourAtomPerturbation::FourAtomPerturbation(const AtomID &atom0, const AtomID &atom1,
                                           const AtomID &atom2, const AtomID &atom3,
                                           const Expression &base_expression,
                                           const Identities &initial_forms,
                                           const Identities &final_forms,
                                           const Expression &mapping_function,
                                           const PropertyMap &map)
     : ConcreteProperty<FourAtomPerturbation,InternalPerturbation>(base_expression,
                                                                   initial_forms,
                                                                   final_forms,
                                                                   mapping_function,
                                                                   map),
       atm0(atom0), atm1(atom1), atm2(atom2), atm3(atom3)
{}

/** Copy constructor */
FourAtomPerturbation::FourAtomPerturbation(const FourAtomPerturbation &other)
             : ConcreteProperty<FourAtomPerturbation,InternalPerturbation>(other),
               atm0(other.atm0), atm1(other.atm1), atm2(other.atm2), atm3(other.atm3)
{}

/** Destructor */
FourAtomPerturbation::~FourAtomPerturbation()
{}

/** Copy assignment operator */
FourAtomPerturbation& FourAtomPerturbation::operator=(const FourAtomPerturbation &other)
{
    if (this != &other)
    {
        atm0 = other.atm0;
        atm1 = other.atm1;
        atm2 = other.atm2;
        atm3 = other.atm3;
        InternalPerturbation::operator=(other);
    }

    return *this;
}

/** Comparison operator */
bool FourAtomPerturbation::operator==(const FourAtomPerturbation &other) const
{
    return ( (atm0 == other.atm0 and atm1 == other.atm1 and
              atm2 == other.atm2 and atm3 == other.atm3) or
             (atm0 == other.atm3 and atm1 == other.atm2 and
              atm2 == other.atm1 and atm3 == other.atm0) ) and
           InternalPerturbation::operator==(other);
}

/** Comparison operator */
bool FourAtomPerturbation::operator!=(const FourAtomPerturbation &other) const
{
    return not FourAtomPerturbation::operator==(other);
}

const char* FourAtomPerturbation::typeName()
{
    return QMetaType::typeName( qMetaTypeId<FourAtomPerturbation>() );
}

/** Return the first of the four atoms whose potential is being changed */
const AtomID& FourAtomPerturbation::atom0() const
{
    return atm0.base();
}

/** Return the second of the four atoms whose potential is being changed */
const AtomID& FourAtomPerturbation::atom1() const
{
    return atm1.base();
}

/** Return the third of the four atoms whose potential is being changed */
const AtomID& FourAtomPerturbation::atom2() const
{
    return atm2.base();
}

/** Return the fourth of the four atoms whose potential is being changed */
const AtomID& FourAtomPerturbation::atom3() const
{
    return atm3.base();
}

/** Return a string representation of this perturbation */
QString FourAtomPerturbation::toString() const
{
    return QObject::tr("FourAtomPerturbation( %1-%2-%3-%4 : %5 )")
            .arg(atm0.toString(), atm1.toString(), atm2.toString(), atm3.toString())
            .arg(perturbExpression().toString());
}

/** Return the properties required or changed by this perturbation */
QSet<QString> FourAtomPerturbation::requiredProperties() const
{
    QSet<QString> props;

    const PropertyName &param_property = propertyMap()["parameters"];

    if (param_property.hasSource())
        props.insert(param_property.source());

    return props;
}

/** Return whether or not this perturbation with the passed values would
    change the molecule 'molecule' */
bool FourAtomPerturbation::wouldChange(const Molecule &molecule,
                                       const Values &values) const
{
    const PropertyName &param_property = propertyMap()["parameters"];

    if (not molecule.hasProperty(param_property))
        return true;

    const FourAtomFunctions &funcs = molecule.property(param_property)
                                             .asA<FourAtomFunctions>();

    Identities idents;

    foreach (Symbol symbol, values.symbols())
    {
        idents.set(symbol, Expression(values[symbol]));
    }

    Expression new_function = perturbExpression().substitute(idents);

    // Use DihedralID as this looks up both atm0-atm1-atm2-atm3
    //                                  and atm3-atm2-atm1-atm0
    return new_function != funcs.potential( DihedralID(atm0, atm1, atm2, atm3) );
}

/** Perturb the four atom function in passed molecule using the reaction
    coordinate(s) in 'values'

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void FourAtomPerturbation::perturbMolecule(MolEditor &molecule,
                                           const Values &values) const
{
    const PropertyName &param_property = propertyMap()["parameters"];

    FourAtomFunctions funcs;

    bool new_property = false;

    if (molecule.hasProperty(param_property))
    {
        funcs = molecule.property(param_property).asA<FourAtomFunctions>();
    }
    else
    {
        funcs = FourAtomFunctions(molecule.data());
        new_property = true;
    }

    //calculate the perturbed function
    Identities idents;

    foreach (Symbol symbol, values.symbols())
    {
        idents.set(symbol, Expression(values[symbol]));
    }

    Expression new_function = perturbExpression().substitute(idents);

    // Use DihedralID as this looks up both atm0-atm1-atm2-atm3 and
    //                                      atm3-atm2-atm1-atm0
    Expression old_function = funcs.potential( DihedralID(atm0, atm1, atm2, atm3) );

    if (new_property or (new_function != old_function))
    {
        funcs.set( DihedralID(atm0, atm1, atm2, atm3), new_function);
        molecule.setProperty(param_property, funcs);
    }
}
