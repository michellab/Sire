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

#include "perturbation.h"
#include "geometryperturbation.h"
#include "molecule.h"
#include "moleditor.h"
#include "mover.hpp"

#include "core.h"

#include "SireCAS/values.h"
#include "SireCAS/identities.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMol;
using namespace SireBase;
using namespace SireCAS;
using namespace SireStream;

//////////
////////// Implementation of PerturbationSymbols
//////////

PerturbationSymbols::PerturbationSymbols()
                    : lam("lambda"), init("initial"), finl("final")
{}

PerturbationSymbols::~PerturbationSymbols()
{}

/** Return the symbol used to represent the driving (reaction)
    coordinate */
const Symbol& PerturbationSymbols::lambda() const
{
    return lam;
}

/** Return the symbol used to represent the initial state */
const Symbol& PerturbationSymbols::initial() const
{
    return init;
}

/** Return the symbol used to represent the final state */
const Symbol& PerturbationSymbols::final() const
{
    return finl;
}

//////////
////////// Implementation of Perturbation
//////////

static const RegisterMetaType<Perturbation> r_pert( MAGIC_ONLY,
                                                    Perturbation::typeName() );

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                       const Perturbation &pert)
{
    writeHeader(ds, r_pert, 1);

    SharedDataStream sds(ds);

    sds << pert.mapping_eqn << pert.map << static_cast<const Property&>(pert);

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, Perturbation &pert)
{
    VersionID v = readHeader(ds, r_pert);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> pert.mapping_eqn >> pert.map >> static_cast<Property&>(pert);
    }
    else
        throw version_error(v, "1", r_pert, CODELOC);

    return ds;
}

static Expression *default_equation(0);

Q_GLOBAL_STATIC_WITH_ARGS( QMutex, globalMutex, (QMutex::Recursive) )
Q_GLOBAL_STATIC( PerturbationSymbols, perturbationSymbols )

/** Return the symbols object that contains the symbols used
    by the mapping equation */
const PerturbationSymbols& Perturbation::symbols()
{
    return *(perturbationSymbols());
}

/** Return the default mapping equation for the perturbations - this
    linearly maps from the initial values at lambda=0 to the
    final value at lambda=1 */
const Expression& Perturbation::defaultFunction()
{
    if (not default_equation)
    {
        QMutexLocker lkr( globalMutex() );

        if (not default_equation)
        {
            default_equation = new Expression();

            *default_equation = ((1 - Perturbation::symbols().lambda()) *
                                            Perturbation::symbols().initial() ) +
                                ( Perturbation::symbols().lambda() *
                                            Perturbation::symbols().final() );
        }
    }

    return *default_equation;
}

/** Constructor */
Perturbation::Perturbation() : Property(), mapping_eqn( defaultFunction() )
{}

Perturbation::Perturbation(const PropertyMap &m)
             : Property(), mapping_eqn( defaultFunction() ), map(m)
{}

Perturbation::Perturbation(const Expression &equation, const PropertyMap &m)
             : Property(), mapping_eqn(equation), map(m)
{}

/** Copy constructor */
Perturbation::Perturbation(const Perturbation &other)
             : Property(other), mapping_eqn(other.mapping_eqn), map(other.map)
{}

/** Destructor */
Perturbation::~Perturbation()
{}

/** Copy assignment operator */
Perturbation& Perturbation::operator=(const Perturbation &other)
{
    if (this != &other)
    {
        mapping_eqn = other.mapping_eqn;
        map = other.map;

        Property::operator=(other);
    }

    return *this;
}

/** Comparison operator */
bool Perturbation::operator==(const Perturbation &other) const
{
    return mapping_eqn == other.mapping_eqn and map == other.map and
           Property::operator==(other);
}

/** Comparison operator */
bool Perturbation::operator!=(const Perturbation &other) const
{
    return not Perturbation::operator==(other);
}

/** Recreate this perturbation - this has the same effect as .clone() */
PerturbationPtr Perturbation::recreate() const
{
    return *this;
}

/** Recreate this perturbation, replacing its current mapping function
    with 'mapping_function' */
PerturbationPtr Perturbation::recreate(const Expression &mapping_function) const
{
    if (mapping_eqn == mapping_function)
    {
        return this->recreate();
    }
    else
    {
        PerturbationPtr new_pert(*this);

        new_pert.edit().mapping_eqn = mapping_function;

        return new_pert;
    }
}

/** Recreate this perturbation, replacing the current property map with 'map' */
PerturbationPtr Perturbation::recreate(const PropertyMap &new_map) const
{
    if (new_map == map)
    {
        return this->recreate();
    }
    else
    {
        PerturbationPtr new_pert(*this);

        new_pert.edit().map = new_map;

        return new_pert;
    }
}

/** Recreate this perturbation, replacing both the mapping function and
    the property map */
PerturbationPtr Perturbation::recreate(const Expression &mapping_function,
                                       const PropertyMap &new_map) const
{
    if (mapping_function == mapping_eqn and new_map == map)
    {
        return this->recreate();
    }
    else
    {
        PerturbationPtr new_pert(*this);

        new_pert.edit().mapping_eqn = mapping_function;
        new_pert.edit().map = new_map;

        return new_pert;
    }
}

/** Substitute the identities in 'identities' in all of the mapping functions
    used by this perturbation. This is useful if, for example, you want to
    switch from using 'lambda' to control the perturbation to using 'alpha', e.g.

    alpha_perturbations = lambda_perturbations.substitute( lam == Expression(alpha) );
*/
PerturbationPtr Perturbation::substitute(const Identities &identities) const
{
    Expression new_mapping_eqn = mapping_eqn.substitute(identities);

    if (new_mapping_eqn != mapping_eqn)
    {
        PerturbationPtr new_pert(*this);

        new_pert.edit().mapping_eqn = new_mapping_eqn;
        return new_pert;
    }
    else
        return *this;
}

/** Substitute the symbol 'old_symbol' with the symbol 'new_symbol'
    in all of the mapping functions used by this perturbation. This is
    useful if, for example, you want to switch from using
    'lambda' to control the perturbation to using 'alpha', e.g.

    alpha_perturbations = lambda_perturbations.substitute(lam, alpha);
*/
PerturbationPtr Perturbation::substitute(const Symbol &old_symbol,
                                         const Symbol &new_symbol) const
{
    return this->substitute( old_symbol == Expression(new_symbol) );
}

/** Return all of the child perturbations that make up
    this perturbation */
QList<PerturbationPtr> Perturbation::children() const
{
    QList<PerturbationPtr> ret;
    ret.append( *this );
    return ret;
}

/** Return all of the symbols that need to be supplied
    to the mapping function (i.e. ignoring symbols().initial()
    and symbols().final() ) */
QSet<Symbol> Perturbation::requiredSymbols() const
{
    QSet<Symbol> syms = mapping_eqn.symbols();
    syms.remove( symbols().initial() );
    syms.remove( symbols().final() );

    return syms;
}

/** Return the equation used to control the mapping from the
    the initial value (represented using symbols().initial()) to
    the final value (represented using symbols().final()) as a
    function of the reaction coordinate (which is normally
    represented using symbols().lambda()) */
const Expression& Perturbation::mappingFunction() const
{
    return mapping_eqn;
}

/** Return the property map used to find the properties used,
    and affected by this perturbation */
const PropertyMap& Perturbation::propertyMap() const
{
    return map;
}

/** Perturb the passed molecule, returning the result

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError:incompatible_error
*/
Molecule Perturbation::perturb(const Molecule &molecule, const Values &values) const
{
    MolEditor editor = molecule.edit();
    this->perturbMolecule(editor, values);
    return editor.commit();
}

const char* Perturbation::typeName()
{
    return "SireMol::Perturbation";
}

Q_GLOBAL_STATIC( SharedPolyPointer<Perturbation>, perturbationPtr );

const NullPerturbation& Perturbation::null()
{
    SharedPolyPointer<Perturbation> *ptr = perturbationPtr();

    if (ptr->constData() == 0)
    {
        QMutexLocker lkr( globalMutex() );

        if (ptr->constData() == 0)
            *ptr = static_cast<Perturbation*>(new NullPerturbation());
    }

    return ptr->constData()->asA<NullPerturbation>();
}

//////////
////////// Implementation of NullPerturbation
//////////

static const RegisterMetaType<NullPerturbation> r_nullpert;

QDataStream &operator<<(QDataStream &ds, const NullPerturbation &nullpert)
{
    writeHeader(ds, r_nullpert, 1);

    ds << static_cast<const Perturbation&>(nullpert);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, NullPerturbation &nullpert)
{
    VersionID v = readHeader(ds, r_nullpert);

    if (v == 1)
    {
        ds >> static_cast<Perturbation&>(nullpert);
    }
    else
        throw version_error(v, "1", r_nullpert, CODELOC);

    return ds;
}

NullPerturbation::NullPerturbation() : ConcreteProperty<NullPerturbation,Perturbation>()
{}

NullPerturbation::NullPerturbation(const NullPerturbation &other)
                 : ConcreteProperty<NullPerturbation,Perturbation>(other)
{}

NullPerturbation::~NullPerturbation()
{}

const char* NullPerturbation::typeName()
{
    return QMetaType::typeName( qMetaTypeId<NullPerturbation>() );
}

NullPerturbation& NullPerturbation::operator=(const NullPerturbation &other)
{
    Perturbation::operator=(other);
    return *this;
}

bool NullPerturbation::operator==(const NullPerturbation &other) const
{
    return Perturbation::operator==(other);
}

bool NullPerturbation::operator!=(const NullPerturbation &other) const
{
    return Perturbation::operator!=(other);
}

QSet<Symbol> NullPerturbation::requiredSymbols() const
{
    return QSet<Symbol>();
}

QSet<QString> NullPerturbation::requiredProperties() const
{
    return QSet<QString>();
}

bool NullPerturbation::wouldChange(const Molecule&, const Values&) const
{
    return false;
}

void NullPerturbation::perturbMolecule(MolEditor&, const Values&) const
{}

//////////
////////// Implementation of Perturbations
//////////

static const RegisterMetaType<Perturbations> r_perts;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                       const Perturbations &perts)
{
    writeHeader(ds, r_perts, 1);

    SharedDataStream sds(ds);

    sds << perts.perts << static_cast<const Perturbation&>(perts);

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, Perturbations &perts)
{
    VersionID v = readHeader(ds, r_perts);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> perts.perts >> static_cast<Perturbation&>(perts);

        perts.makeSane();
    }
    else
        throw version_error(v, "1", r_perts, CODELOC);

    return ds;
}

/** This removes all null perturbations and collapses GeometryPerturbation
    objects into a single GeometryPerturbations */
void Perturbations::makeSane()
{
    QList<PerturbationPtr> norm_perts;
    QList<GeomPertPtr> geom_perts;

    QList<PerturbationPtr> kids = this->children();

    for (QList<PerturbationPtr>::const_iterator it = kids.constBegin();
         it != kids.constEnd();
         ++it)
    {
        if ((*it)->isA<NullPerturbation>())
            continue;

        else if ((*it)->isA<GeometryPerturbation>() and
                 not (*it)->isA<GeometryPerturbations>())
        {
            geom_perts.append( (*it)->asA<GeometryPerturbation>() );
        }
        else if (not (*it)->isA<Perturbations>())
        {
            norm_perts.append(*it);
        }
    }

    norm_perts.append( GeometryPerturbations(geom_perts) );

    perts = norm_perts;
}

/** Constructor */
Perturbations::Perturbations()
              : ConcreteProperty<Perturbations,Perturbation>(Expression())
{}

/** Construct just to perform the passed perturbation */
Perturbations::Perturbations(const Perturbation &perturbation)
              : ConcreteProperty<Perturbations,Perturbation>(Expression())
{
    perts.append(perturbation);
    this->makeSane();
}

/** Construct to perform the passed perturbations */
Perturbations::Perturbations(const QList<PerturbationPtr> &perturbations)
              : ConcreteProperty<Perturbations,Perturbation>(Expression()),
                perts(perturbations)
{
    this->makeSane();
}

/** Copy constructor */
Perturbations::Perturbations(const Perturbations &other)
              : ConcreteProperty<Perturbations,Perturbation>(other),
                perts(other.perts)
{}

/** Destructor */
Perturbations::~Perturbations()
{}

/** Copy assignment operator */
Perturbations& Perturbations::operator=(const Perturbations &other)
{
    perts = other.perts;
    Perturbation::operator=(other);
    return *this;
}

/** Comparison operator */
bool Perturbations::operator==(const Perturbations &other) const
{
    return perts == other.perts and Perturbation::operator==(other);
}

/** Comparison operator */
bool Perturbations::operator!=(const Perturbations &other) const
{
    return perts != other.perts or Perturbation::operator!=(other);
}

const char* Perturbations::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Perturbations>() );
}

QString Perturbations::toString() const
{
    QStringList lines;

    if (perts.isEmpty())
        return QObject::tr("Perturbations::null");

    lines.append( QObject::tr("Perturbations:") );

    foreach (PerturbationPtr pert, perts)
    {
        lines.append( QString("   %1").arg(pert->toString()) );
    }

    return lines.join("\n");
}

/** Return a list of all perturbations performed by this object */
QList<PerturbationPtr> Perturbations::perturbations() const
{
    return perts;
}

/** Return a re-created version of this set of perturbations where all child
    perturbations are changed to use the passed mapping function */
PerturbationPtr Perturbations::recreate(const Expression &mapping_function) const
{
    QList<PerturbationPtr> new_perts;

    for (QList<PerturbationPtr>::const_iterator it = perts.constBegin();
         it != perts.constEnd();
         ++it)
    {
        new_perts.append( it->read().recreate(mapping_function) );
    }

    Perturbations ret(*this);
    ret.perts = new_perts;

    return ret;
}

/** Return a re-created version of this set of perturbations where all child
    perturbations are changed to use the passed property map */
PerturbationPtr Perturbations::recreate(const PropertyMap &map) const
{
    QList<PerturbationPtr> new_perts;

    for (QList<PerturbationPtr>::const_iterator it = perts.constBegin();
         it != perts.constEnd();
         ++it)
    {
        new_perts.append( it->read().recreate(map) );
    }

    Perturbations ret(*this);
    ret.perts = new_perts;

    return ret;
}

/** Return a re-created version of this set of perturbations where all child
    perturbations are changed to use the passed mapping function and property map */
PerturbationPtr Perturbations::recreate(const Expression &mapping_function,
                                        const PropertyMap &map) const
{
    QList<PerturbationPtr> new_perts;

    for (QList<PerturbationPtr>::const_iterator it = perts.constBegin();
         it != perts.constEnd();
         ++it)
    {
        new_perts.append( it->read().recreate(mapping_function,map) );
    }

    Perturbations ret(*this);
    ret.perts = new_perts;

    return ret;
}

/** Substitute the identities in 'identities' in all of the mapping functions
    used by this perturbation. This is useful if, for example, you want to
    switch from using 'lambda' to control the perturbation to using 'alpha', e.g.

    alpha_perturbations = lambda_perturbations.substitute( lam == Expression(alpha) );
*/
PerturbationPtr Perturbations::substitute(const Identities &identities) const
{
    QList<PerturbationPtr> new_perts;

    for (QList<PerturbationPtr>::const_iterator it = perts.constBegin();
         it != perts.constEnd();
         ++it)
    {
        new_perts.append( it->read().substitute(identities) );
    }

    Perturbations ret(*this);
    ret.perts = new_perts;

    return ret;
}

PerturbationPtr Perturbations::substitute(const SireCAS::Symbol &old_symbol,
                                          const SireCAS::Symbol &new_symbol) const
{
    return Perturbation::substitute(old_symbol, new_symbol);
}

/** Return a list of all of the children of this perturbation
    (and the children of these children etc.) */
QList<PerturbationPtr> Perturbations::children() const
{
    QList<PerturbationPtr> kids;

    for (QList<PerturbationPtr>::const_iterator it = perts.constBegin();
         it != perts.constEnd();
         ++it)
    {
        kids += it->read().children();
    }

    return kids;
}

/** Return all of the symbols that need to be input to these perturbations */
QSet<Symbol> Perturbations::requiredSymbols() const
{
    QSet<Symbol> syms;

    for (QList<PerturbationPtr>::const_iterator it = perts.constBegin();
         it != perts.constEnd();
         ++it)
    {
        syms += it->read().requiredSymbols();
    }

    return syms;
}

/** Return all of the properties that are needed or affected by
    these perturbations */
QSet<QString> Perturbations::requiredProperties() const
{
    QSet<QString> props;

    for (QList<PerturbationPtr>::const_iterator it = perts.constBegin();
         it != perts.constEnd();
         ++it)
    {
        props += it->read().requiredProperties();
    }

    return props;
}

/** Return whether or not these perturbations with the passed values would
    change the molecule 'molecule' */
bool Perturbations::wouldChange(const Molecule &molecule, const Values &values) const
{
    try
    {
        for (QList<PerturbationPtr>::const_iterator it = perts.constBegin();
             it != perts.constEnd();
             ++it)
        {
            if (it->read().wouldChange(molecule,values))
                return true;
        }

        return false;
    }
    catch(...)
    {
        //if an error occured, then the molecule won't be changed
        return false;
    }
}

/** Apply this perturbation to the passed molecule for the
    specified lambda value */
void Perturbations::perturbMolecule(MolEditor &molecule, const Values &values) const
{
    for (QList<PerturbationPtr>::const_iterator it = perts.constBegin();
         it != perts.constEnd();
         ++it)
    {
        it->read().perturbMolecule(molecule, values);
    }
}
