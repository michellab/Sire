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

#include "constraint.h"
#include "system.h"
#include "delta.h"

#include "SireBase/numberproperty.h"
#include "SireBase/propertylist.h"

#include "SireMaths/maths.h"

#include "SireSystem/errors.h"
#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"
#include "SireStream/streamdata.hpp"

#include <QDebug>

using namespace SireSystem;
using namespace SireBase;
using namespace SireStream;

//////////
////////// Implementation of Constraint
//////////

static const RegisterMetaType<Constraint> r_constraint( MAGIC_ONLY,
                                                        "SireSystem::Constraint" );

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const Constraint &constraint)
{
    writeHeader(ds, r_constraint, 1);
    
    ds << static_cast<const Property&>(constraint);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, Constraint &constraint)
{
    VersionID v = readHeader(ds, r_constraint);
    
    if (v == 1)
    {
        ds >> static_cast<Property&>(constraint);
        constraint.clearLastSystem();
    }
    else
        throw version_error( v, "1", r_constraint, CODELOC );
        
    return ds;
}

/** Constructor */
Constraint::Constraint() : Property(), last_subversion(0), last_was_satisfied(false)
{}

/** Copy constructor */
Constraint::Constraint(const Constraint &other) 
           : Property(other),
             last_sysuid(other.last_sysuid), 
             last_sysversion(other.last_sysversion),
             last_subversion(other.last_subversion),
             last_was_satisfied(other.last_was_satisfied)
{}

/** Destructor */
Constraint::~Constraint()
{}

/** Copy assignment operator */
Constraint& Constraint::operator=(const Constraint &other)
{
    if (this != &other)
    {
        Property::operator=(other);
        last_sysuid = other.last_sysuid;
        last_sysversion = other.last_sysversion;
        last_was_satisfied = other.last_was_satisfied;
        last_subversion = other.last_subversion;
    }
    
    return *this;
}

/** Return whether or not this constraint *may* affect the passed delta */
bool Constraint::mayAffect(const Delta &delta) const
{
    const System &system = delta.deltaSystem();

    if ( Constraint::wasLastSatisfied() and Constraint::wasLastSystem(system) )
    {
        if ( Constraint::wasLastSubVersion(system) )
            return false;
        else
            return this->mayChange(delta, last_subversion);
    }
    else
    {
        return true;
    }
}

/** Internal function to apply this constraint to the passed delta,
    returning whether this constraint changes the delta */
bool Constraint::apply(Delta &delta)
{
    try
    {
        const System &system = delta.deltaSystem();

        if ( not Constraint::wasLastSystem(system) )
            this->setSystem(system);

        bool changed = false;

        if ( not Constraint::wasLastSatisfied() )
        {
            if ( not Constraint::wasLastSubVersion(system) )
            {
                this->setSystem(system);
            
                if (not Constraint::wasLastSatisfied())
                    changed = this->fullApply(delta);
            }
            else
                changed = this->fullApply(delta);
        }
        else if ( not Constraint::wasLastSubVersion(system) )
        {
            if ( this->mayChange(delta, last_subversion) )
                changed = this->deltaApply(delta, last_subversion);
        }

        if (not this->isSatisfied(delta.deltaSystem()))
            throw SireError::program_bug( QObject::tr(
                "Constraint %1 is not satisfied despite having just been applied!!!")
                    .arg(this->toString()), CODELOC );
    
        //the above function *MUST* have set the system to satisfy the constraint
        setSatisfied(delta.deltaSystem(), true);
        
        return changed;
    }
    catch(...)
    {
        clearLastSystem();
        throw;
    }
    
    return false;
}

/** Apply this constraint to the passed system, returning
    a new system in which the constraint is satisfied */
System Constraint::apply(const System &system)
{
    this->setSystem(system);
    
    if (wasLastSatisfied())
        //nothing to do
        return system;
        
    Delta delta(system);

    System new_system = system;

    for (int i=0; i<10; ++i)
    {
        if (this->apply(delta))
        {
            new_system = delta.apply();
        }
        else
        {
            //the system satisfies this constraint
            this->committed(new_system);
            return new_system;
        }
    }
    
    throw SireSystem::constraint_error( QObject::tr(
            "The constraint %1 could not be satisfied in concert with the "
            "constraints already present in the system %2 (%3)")
                .arg(this->toString(), system.toString(), 
                     system.constraints().toString()), CODELOC );
                     
    return System();
}

/** Function called by Delta to say that the passed system has just 
    been committed (which means that it has passed this constraint) */
void Constraint::committed(const System &system)
{
    this->setSatisfied(system, true);
}

/** Return whether or not this constraint is satisfied for 
    the passed system */
bool Constraint::isSatisfied(const System &system) const
{
    if ( wasLastSystem(system) and wasLastSubVersion(system) )
        return last_was_satisfied;

    else
    {
        #ifdef BOOST_NO_CXX11_SMART_PTR
          std::auto_ptr<Constraint> copy( this->clone() );
        #else
          std::unique_ptr<Constraint> copy( this->clone() );
        #endif 
       
        copy->setSystem(system);
        
        if (not (copy->wasLastSystem(system) and copy->wasLastSubVersion(system)))
        {
            throw SireError::program_bug( QObject::tr(
                    "Error with setSystem for constraint %1. This should "
                    "set the last system to %2:%3.%4, but is has instead "
                    "set it to %5:%6.%7 (%8).")
                        .arg(copy->toString())
                        .arg(system.UID().toString())
                        .arg(system.version().toString())
                        .arg(system.subVersion())
                        .arg(copy->last_sysuid.toString())
                        .arg(copy->last_sysversion.toString())
                        .arg(copy->last_subversion)
                        .arg(copy->last_was_satisfied), CODELOC );
        }
                 
        return copy->last_was_satisfied;
    }
}

/** Internal function called by the constraint to say 
    whether or not it is satisfied on the passed system */
void Constraint::setSatisfied(const System &system, bool is_satisfied)
{
    last_was_satisfied = is_satisfied;
    last_sysuid = system.UID();
    last_sysversion = system.version();
    last_subversion = system.subVersion();
}

/** Clear the cache for the last system */
void Constraint::clearLastSystem()
{
    last_was_satisfied = false;
    last_sysuid = QUuid();
    last_sysversion = Version();
    last_subversion = 0;
}

/** Return whether or not the passed system was the last one
    on which this constraint was applied - note that this
    ignores the subversion number - use
    wasLastSystem(system) and wasLastSubVersion(system)
    to test if this is exactly the same as the last system */
bool Constraint::wasLastSystem(const System &system) const
{
    return system.UID() == last_sysuid and
           system.version() == last_sysversion;
}

/** Return whether or not the passed system has the same subversion
    as the last system set to this constraint */
bool Constraint::wasLastSubVersion(const System &system) const
{
    return system.subVersion() == last_subversion;
}

/** Return whether or not the last system satisfied this constraint */
bool Constraint::wasLastSatisfied() const
{
    return last_was_satisfied;
}

/** Assert that the constraint is satisfied in the passed system 

    \throw SireSystem::constraint_error
*/
void Constraint::assertSatisfied(const System &system) const
{
    if (not this->isSatisfied(system))
        throw SireSystem::constraint_error( QObject::tr(
            "The constraint %1 is not maintained in the system %2.")
                .arg(this->toString())
                .arg(system.toString()), CODELOC );
}

//////////
////////// Implementation of NullConstraint
//////////

static const RegisterMetaType<NullConstraint> r_nullconstraint;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, 
                                          const NullConstraint &nullconstraint)
{
    writeHeader(ds, r_nullconstraint, 1);
    
    ds << static_cast<const Constraint&>(nullconstraint);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, 
                                          NullConstraint &nullconstraint)
{
    VersionID v = readHeader(ds, r_nullconstraint);
    
    if (v == 1)
    {
        ds >> static_cast<Constraint&>(nullconstraint);
    }
    else
        throw version_error(v, "1", r_nullconstraint, CODELOC);
        
    return ds;
}

/** Null constructor */
NullConstraint::NullConstraint() : ConcreteProperty<NullConstraint,Constraint>()
{}

/** Copy constructor */
NullConstraint::NullConstraint(const NullConstraint &other)
               : ConcreteProperty<NullConstraint,Constraint>(other)
{}

/** Destructor */
NullConstraint::~NullConstraint()
{}

/** Copy assignment operator */
NullConstraint& NullConstraint::operator=(const NullConstraint &other)
{
    Constraint::operator=(other);
    return *this;
}

/** Comparison operator */
bool NullConstraint::operator==(const NullConstraint&) const
{
    return true;
}

/** Comparison operator */
bool NullConstraint::operator!=(const NullConstraint&) const
{
    return false;
}

/** Return a string representation */
QString NullConstraint::toString() const
{
    return QObject::tr("NullConstraint");
}

/** Set the baseline system for the constraint - this is 
    used to pre-calculate everything for the system
    and to check if the constraint is satisfied */
void NullConstraint::setSystem(const System &system)
{
    Constraint::setSatisfied(system, true);
}

/** Return whether or not the changes in the passed
    delta *may* have changed the system since the system
    with subversion 'last_subversion' */
bool NullConstraint::mayChange(const Delta&, quint32) const
{
    return false;
}

/** Fully apply this constraint on the passed delta - this returns
    whether or not this constraint affects the delta */
bool NullConstraint::fullApply(Delta&)
{
    return false;
}

/** Apply this constraint based on the delta, knowing that the 
    last application of this constraint was on this system, 
    at subversion number last_subversion */
bool NullConstraint::deltaApply(Delta&, quint32)
{
    return false;
}

const NullConstraint& Constraint::null()
{
    return *(create_shared_null<NullConstraint>());
}

const char* NullConstraint::typeName()
{
    return QMetaType::typeName( qMetaTypeId<NullConstraint>() );
}

//////////
////////// Implementation of PropertyConstraint
//////////

static const RegisterMetaType<PropertyConstraint> r_propconstraint;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, 
                                          const PropertyConstraint &constraint)
{
    writeHeader(ds, r_propconstraint, 1);
    
    SharedDataStream sds(ds);
    
    sds << constraint.ffid << constraint.propname
        << constraint.eqn
        << static_cast<const Constraint&>(constraint);
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds,
                                          PropertyConstraint &constraint)
{
    VersionID v = readHeader(ds, r_propconstraint);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);

        //completely clear any cached data
        constraint = PropertyConstraint();
        
        sds >> constraint.ffid >> constraint.propname
            >> constraint.eqn
            >> static_cast<Constraint&>(constraint);
            
        constraint.syms = constraint.eqn.symbols();
    }
    else
        throw version_error(v, "1", r_propconstraint, CODELOC);
        
    return ds;
}

/** Null constructor */
PropertyConstraint::PropertyConstraint()
                   : ConcreteProperty<PropertyConstraint,Constraint>(),
                     target_value(0)
{}

/** Construct to constrain the property with name 'name' in all forcefields
    to have the value resulting from the expression 'expression' */
PropertyConstraint::PropertyConstraint(const QString &name, 
                                       const SireCAS::Expression &expression)
                   : ConcreteProperty<PropertyConstraint,Constraint>(),
                     propname(name), eqn(expression), syms(expression.symbols()),
                     target_value(0)
{}

/** Construct to constrain the property with name 'name' in the forcefield(s)
    that match the ID 'ffid' to have the value resulting the expression
    'expression' */
PropertyConstraint::PropertyConstraint(const QString &name, const SireFF::FFID &id,
                                       const SireCAS::Expression &expression)
                   : ConcreteProperty<PropertyConstraint,Constraint>(),
                     ffid(id), propname(name), eqn(expression),
                     syms(expression.symbols()), target_value(0)
{}

/** Copy constructor */
PropertyConstraint::PropertyConstraint(const PropertyConstraint &other)
                   : ConcreteProperty<PropertyConstraint,Constraint>(other),
                     ffid(other.ffid), ffidxs(other.ffidxs),
                     propname(other.propname), eqn(other.eqn),
                     syms(other.syms),
                     component_vals(other.component_vals),
                     constrained_value(other.constrained_value),
                     target_value(other.target_value)
{}

/** Destructor */
PropertyConstraint::~PropertyConstraint()
{}

/** Copy assignment operator */
PropertyConstraint& PropertyConstraint::operator=(const PropertyConstraint &other)
{
    if (this != &other)
    {
        Constraint::operator=(other);
        ffid = other.ffid;
        ffidxs = other.ffidxs;
        propname = other.propname;
        eqn = other.eqn;
        syms = other.syms;
        component_vals = other.component_vals;
        constrained_value = other.constrained_value;
        target_value = other.target_value;
    }
    
    return *this;
}

/** Comparison operator */
bool PropertyConstraint::operator==(const PropertyConstraint &other) const
{
    return ffid == other.ffid and propname == other.propname and 
           eqn == other.eqn;
}

/** Comparison operator */
bool PropertyConstraint::operator!=(const PropertyConstraint &other) const
{
    return not this->operator==(other);
}

/** Return a string representation of the constraint */
QString PropertyConstraint::toString() const
{
    return QObject::tr("PropertyConstraint( FFID=%1 property=%2 expression=%3 )")
                .arg(ffid.toString(), propname, eqn.toString());
}

static PropertyPtr getProperty(const System &system, const QString &propname,
                               const QList<FFIdx> &ffidxs)
{
    bool is_first = true;
    
    PropertyPtr property;
    
    foreach (FFIdx ffidx, ffidxs)
    {
        if (system.containsProperty(ffidx, propname))
        {
            const Property &p = system.property(ffidx, propname);
            
            if (is_first)
            {
                property = p;
                is_first = false;
            }
            else
            {
                if (property.isNull())
                    return PropertyPtr();
                
                else if (not property.read().equals(p))
                    return PropertyPtr();
            }
        }
        else if (not is_first)
        {
            if (not property.isNull())
                return PropertyPtr();
        }
    }
    
    return property;
}

/** Return whether or not the changes in the passed
    delta *may* have changed the system since the last
    subversion 'subversion' */
bool PropertyConstraint::mayChange(const Delta &delta, quint32 last_subversion) const
{
    BOOST_ASSERT( Constraint::wasLastSystem(delta.deltaSystem()) );

    return (not Constraint::wasLastSatisfied()) or
           delta.sinceChanged(propname, last_subversion) or 
           delta.sinceChanged(syms, last_subversion);
}

/** Set the baseline system for the constraint - this is 
    used to pre-calculate everything for the system
    and to check if the constraint is satisfied */
void PropertyConstraint::setSystem(const System &system)
{
    if (Constraint::wasLastSystem(system) and Constraint::wasLastSubVersion(system))
        return;
    
    Constraint::clearLastSystem();
    
    if (not ffid.isNull())
        ffidxs = ffid.map(system.forceFields());
        
    constrained_value = PropertyPtr();
    
    if (not ffidxs.isEmpty())
    {
        constrained_value = ::getProperty(system, propname, ffidxs);
    }
    else 
    {
        if (system.containsProperty(propname))
            constrained_value = system.property(propname);
    }

    double constrained_val = 0;
    bool has_constrained_value = false;

    if (not constrained_value.isNull())
    {
        constrained_val = constrained_value.read().asADouble();
        has_constrained_value = true;
    }
    
    component_vals = system.constants(syms);
    
    target_value = eqn(component_vals);

    Constraint::setSatisfied(system, has_constrained_value and 
                                     SireMaths::areEqual(constrained_val, target_value));
}

/** Fully apply this constraint on the passed delta - this returns
    whether or not this constraint affects the delta */
bool PropertyConstraint::fullApply(Delta &delta)
{
    BOOST_ASSERT( Constraint::wasLastSystem(delta.deltaSystem()) and
                  Constraint::wasLastSubVersion(delta.deltaSystem()) );

    bool changed = false;

    if (ffidxs.isEmpty())
        changed = delta.update(propname, wrap(target_value));
    else
        changed = delta.update(propname, ffidxs, wrap(target_value));
    
    if (changed)
        this->setSystem( delta.deltaSystem() );
    
    return changed;
}

/** Apply this constraint based on the delta, knowing that the 
    last application of this constraint was on this system, 
    at subversion number last_subversion */
bool PropertyConstraint::deltaApply(Delta &delta, quint32 last_subversion)
{
    BOOST_ASSERT( Constraint::wasLastSystem(delta.deltaSystem()) );

    bool changed_prop = delta.sinceChanged(propname, last_subversion);
    bool changed_syms = delta.sinceChanged(syms, last_subversion);
    bool changed_target = false;

    if (changed_prop or changed_syms)
    {
        const System &system = delta.deltaSystem();

        if (changed_syms)
        {
            component_vals = system.constants(syms);
            double new_target = eqn(component_vals);
                
            if (new_target != target_value)
            {
                target_value = new_target;
                changed_target = true;
            }
        }
            
        if (changed_prop)
        {
            if (ffidxs.isEmpty())
            {
                if (system.containsProperty(propname))
                {
                    const Property &new_property = system.property(propname);
                        
                    if (constrained_value.isNull() or 
                          (not constrained_value.read().equals(new_property)))
                    {
                        constrained_value = new_property;
                    }
                    else
                        changed_prop = false;
                }
                else
                    constrained_value = PropertyPtr();
            }
            else
            {
                PropertyPtr new_property = ::getProperty(system, propname, ffidxs);
                                        
                if ( constrained_value.isNull() or new_property.isNull() or
                      (not constrained_value.read().equals(new_property)) )
                {
                    constrained_value = new_property;
                }
                else
                    changed_prop = false;
            }
        }
            
        if (changed_target or changed_prop)
        {
            //get the current value of the property
            if (not constrained_value.isNull())
            {
                if (target_value == constrained_value.read().asADouble())
                {
                    return false;
                }
            }
            
            if (ffidxs.isEmpty())
                return delta.update(propname, wrap(target_value));
            else
                return delta.update(propname, ffidxs, wrap(target_value));
        }
    }
    else if (not this->isSatisfied(delta.deltaSystem()))
    {
        return this->fullApply(delta);
    }

    return false;
}

const char* PropertyConstraint::typeName()
{
    return QMetaType::typeName( qMetaTypeId<PropertyConstraint>() );
}

//////////
////////// Implementation of ComponentConstraint
//////////

static const RegisterMetaType<ComponentConstraint> r_compconstraint;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, 
                                          const ComponentConstraint &constraint)
{
    writeHeader(ds, r_compconstraint, 1);
    
    SharedDataStream sds(ds);
    
    sds << constraint.constrained_component
        << constraint.eqn
        << static_cast<const Constraint&>(constraint);
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds,
                                          ComponentConstraint &constraint)
{
    VersionID v = readHeader(ds, r_compconstraint);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        constraint = ComponentConstraint();
        
        sds >> constraint.constrained_component
            >> constraint.eqn
            >> static_cast<Constraint&>(constraint);
            
        constraint.syms = constraint.eqn.symbols();
    }
    else
        throw version_error(v, "1", r_compconstraint, CODELOC);
        
    return ds;
}

/** Null constructor */
ComponentConstraint::ComponentConstraint()
                    : ConcreteProperty<ComponentConstraint,Constraint>(),
                      constrained_value(0), target_value(0),
                      has_constrained_value(false)
{}

/** Construct to constrain the component with symbol 'component'
    to have the value resulting from the expression 'expression' */
ComponentConstraint::ComponentConstraint(const Symbol &component,
                                         const SireCAS::Expression &expression)
                   : ConcreteProperty<ComponentConstraint,Constraint>(),
                     constrained_component(component), eqn(expression),
                     syms(expression.symbols()),
                     constrained_value(0), target_value(0),
                     has_constrained_value(false)
{}

/** Copy constructor */
ComponentConstraint::ComponentConstraint(const ComponentConstraint &other)
                   : ConcreteProperty<ComponentConstraint,Constraint>(other),
                     constrained_component(other.constrained_component), eqn(other.eqn),
                     syms(other.syms), component_vals(other.component_vals),
                     constrained_value(other.constrained_value),
                     target_value(other.target_value),
                     has_constrained_value(other.has_constrained_value)
{}

/** Destructor */
ComponentConstraint::~ComponentConstraint()
{}

/** Copy assignment operator */
ComponentConstraint& ComponentConstraint::operator=(const ComponentConstraint &other)
{
    if (this != &other)
    {
        Constraint::operator=(other);
        constrained_component = other.constrained_component;
        eqn = other.eqn;
        syms = other.syms;
        component_vals = other.component_vals;
        constrained_value = other.constrained_value;
        target_value = other.target_value;
        has_constrained_value = other.has_constrained_value;
    }
    
    return *this;
}

/** Comparison operator */
bool ComponentConstraint::operator==(const ComponentConstraint &other) const
{
    return constrained_component == other.constrained_component and 
           eqn == other.eqn;
}

/** Comparison operator */
bool ComponentConstraint::operator!=(const ComponentConstraint &other) const
{
    return not this->operator==(other);
}

/** Return a string representation of the constraint */
QString ComponentConstraint::toString() const
{
    return QObject::tr("ComponentConstraint( component=%1 expression=%2 )")
                .arg(constrained_component.toString(), eqn.toString());
}

/** Return the symbol representing the component being constrained */
const Symbol& ComponentConstraint::component() const
{
    return constrained_component;
}

/** Return the expression used to evaluate the constraint */
const Expression& ComponentConstraint::expression() const
{
    return eqn;
}

/** Return whether or not the changes in the passed
    delta *may* have changed the system since the last
    subversion 'subversion' */
bool ComponentConstraint::mayChange(const Delta &delta, quint32 last_subversion) const
{
    BOOST_ASSERT( Constraint::wasLastSystem(delta.deltaSystem()) );

    return true;

    return (not Constraint::wasLastSatisfied()) or
           delta.sinceChanged(constrained_component, last_subversion) or
           delta.sinceChanged(syms, last_subversion);
}

/** Set the baseline system for the constraint - this is 
    used to pre-calculate everything for the system
    and to check if the constraint is satisfied */
void ComponentConstraint::setSystem(const System &system)
{
    if ( Constraint::wasLastSystem(system) and Constraint::wasLastSubVersion(system) )
        return;
    
    Constraint::clearLastSystem();
    constrained_value = 0;
    target_value = 0;
    has_constrained_value = false;

    if (system.hasConstantComponent(constrained_component))
    {
        constrained_value = system.constant(constrained_component);
        has_constrained_value = true;
    }
        
    component_vals = system.constants(syms);
    
    target_value = eqn(component_vals);
    
    Constraint::setSatisfied(system, has_constrained_value and
                                     SireMaths::areEqual(target_value, constrained_value));
}

/** Fully apply this constraint on the passed delta - this returns
    whether or not this constraint affects the delta */
bool ComponentConstraint::fullApply(Delta &delta)
{
    BOOST_ASSERT( Constraint::wasLastSystem(delta.deltaSystem()) and
                  Constraint::wasLastSubVersion(delta.deltaSystem()) );

    bool changed = delta.update(constrained_component, target_value);
    
    if (changed)
        this->setSystem(delta.deltaSystem());
    
    return changed;
}

/** Apply this constraint based on the delta, knowing that the 
    last application of this constraint was on this system, 
    at subversion number last_subversion */
bool ComponentConstraint::deltaApply(Delta &delta, quint32 last_subversion)
{
    BOOST_ASSERT( Constraint::wasLastSystem(delta.deltaSystem()) );

    bool changed_comp = delta.sinceChanged(constrained_component, last_subversion);
    bool changed_syms = delta.sinceChanged(syms, last_subversion);
    bool changed_target = false;
    
    if (changed_comp or changed_syms)
    {
        const System &system = delta.deltaSystem();
        
        if (changed_syms)
        {
            component_vals = system.constants(syms);
            double new_target = eqn(component_vals);
            
            if (new_target != target_value)
            {
                target_value = new_target;
                changed_target = true;
            }
            
            if (system.hasConstantComponent(constrained_component))
            {
                double new_comp = system.constant(constrained_component);
                
                if ( (not has_constrained_value) or constrained_value != new_comp )
                {
                    has_constrained_value = true;
                    constrained_value = new_comp;
                    changed_comp = true;
                }
                else
                    changed_comp = false;
            }
            else
            {
                has_constrained_value = false;
                constrained_value = 0;
                changed_comp = true;
            }
        }
        else if (changed_comp)
        {
            if (system.hasConstantComponent(constrained_component))
            {
                double new_comp = system.constant(constrained_component);
                
                if ( (not has_constrained_value) or constrained_value != new_comp )
                {
                    has_constrained_value = true;
                    constrained_value = new_comp;
                    changed_comp = true;
                }
                else
                    changed_comp = false;
            }
            else
            {
                has_constrained_value = false;
                constrained_value = 0;
                changed_comp = true;
            }
        }

        if ( (changed_target or changed_comp) and not
             (has_constrained_value and (constrained_value == target_value)) )
        {
            return delta.update(constrained_component, target_value);
        }
    }
    else if (not this->isSatisfied(delta.deltaSystem()))
    {
        return this->fullApply(delta);
    }
    
    return false;
}

const char* ComponentConstraint::typeName()
{
    return QMetaType::typeName( qMetaTypeId<ComponentConstraint>() );
}

//////////
////////// Implementation of WindowedComponent
//////////

static const RegisterMetaType<WindowedComponent> r_windowedcomp;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, 
                                          const WindowedComponent &windowedcomp)
{
    writeHeader(ds, r_windowedcomp, 1);
    
    SharedDataStream sds(ds);
    
    sds << windowedcomp.constrained_component
        << windowedcomp.reference_component
        << windowedcomp.window_values
        << windowedcomp.step_size
        << static_cast<const Constraint&>(windowedcomp);
        
    return ds;
}
              
/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, 
                                          WindowedComponent &windowedcomp)
{
    VersionID v = readHeader(ds, r_windowedcomp);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        windowedcomp = WindowedComponent();
        
        sds >> windowedcomp.constrained_component
            >> windowedcomp.reference_component
            >> windowedcomp.window_values
            >> windowedcomp.step_size
            >> static_cast<Constraint&>(windowedcomp);
    }
    else
        throw version_error(v, "1", r_windowedcomp, CODELOC);
        
    return ds;
}

/** Null constructor */
WindowedComponent::WindowedComponent()
                  : ConcreteProperty<WindowedComponent,Constraint>(),
                    step_size(0), component_val(0), constrained_value(0),
                    target_value(0), has_constrained_value(false)
{}

/** Construct a WindowedConstraint that constrains the component represented
    by the symbol 'component' to lie 'step_size' windows above (or below if
    'step_size' is negative) the window in which the component represented
    by the symbol 'reference' resides - where 'values' contains all of
    the values of the windows, in the order that you have arranged
    them. Whilst this will not sort 'window_values', it will remove
    all duplicate values */
WindowedComponent::WindowedComponent(const SireCAS::Symbol &component,
                                     const SireCAS::Symbol &reference,
                                     const QVector<double> &values,
                                     int step)
                  : ConcreteProperty<WindowedComponent,Constraint>(),
                    constrained_component(component),
                    reference_component(reference),
                    window_values(values),
                    step_size(step),
                    component_val(0),
                    constrained_value(0), target_value(0),
                    has_constrained_value(false)
{
    if (component == reference)
        throw SireError::invalid_arg( QObject::tr(
                "You cannot use the component %1 as both the constrained "
                "and reference components in a WindowedComponent.")
                    .arg(component.toString()), CODELOC );

    int i = 0;

    while(true)
    {
        if (i >= window_values.count())
            break;
    
        while (window_values.count( window_values[i] ) > 1)
        {
            window_values.remove( window_values.lastIndexOf(window_values[i]) );
        }
        
        ++i;
    }
} 

/** Copy constructor */
WindowedComponent::WindowedComponent(const WindowedComponent &other)
                  : ConcreteProperty<WindowedComponent,Constraint>(other),
                    constrained_component(other.constrained_component),
                    reference_component(other.reference_component),
                    window_values(other.window_values),
                    step_size(other.step_size),
                    component_val(other.component_val),
                    constrained_value(other.constrained_value),
                    target_value(other.target_value),
                    has_constrained_value(other.has_constrained_value)
{}

/** Destructor */
WindowedComponent::~WindowedComponent()
{}

/** Copy assignment operator */
WindowedComponent& WindowedComponent::operator=(const WindowedComponent &other)
{
    if (this != &other)
    {
        constrained_component = other.constrained_component;
        reference_component = other.reference_component;
        window_values = other.window_values;
        step_size = other.step_size;
        component_val = other.component_val;
        constrained_value = other.constrained_value;
        target_value = other.target_value;
        has_constrained_value = other.has_constrained_value;
            
        Constraint::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool WindowedComponent::operator==(const WindowedComponent &other) const
{
    return constrained_component == other.constrained_component and
           reference_component == other.reference_component and
           window_values == other.window_values and
           step_size == other.step_size;
}

/** Comparison operator */
bool WindowedComponent::operator!=(const WindowedComponent &other) const
{
    return not this->operator==(other);
}

/** Return a string representation of this constraint */
QString WindowedComponent::toString() const
{
    return QObject::tr("WindowedComponent( component=%1 reference=%2 step_size=%3 )\n"
                       "   windowed_values = %4")
                .arg(constrained_component.toString(), reference_component.toString())
                .arg(step_size).arg(Sire::toString(window_values));
}

/** Return the symbol representing the component being constrained */
const SireCAS::Symbol& WindowedComponent::component() const
{
    return constrained_component;
}

/** Return the symbol representing the reference component */
const SireCAS::Symbol& WindowedComponent::referenceComponent() const
{
    return reference_component;
}

/** Return the values of all of the windows */
const QVector<double>& WindowedComponent::windowValues() const
{
    return window_values;
}

/** Return the step size for this windows - this is the number of 
    windows above (or below if step_size is negative) for this 
    window compared to the window containing the reference component */
int WindowedComponent::stepSize() const
{
    return step_size;
}

static double getNextWindow(double reference_value,
                            const QVector<double> &window_values,
                            int step_size)
{
    int idx = window_values.indexOf(reference_value);
    
    if (idx == -1)
    {
        //find the closest value
        double del2 = std::numeric_limits<double>::max();
        
        for (int i=0; i<window_values.count(); ++i)
        {
            double my_del2 = SireMaths::pow_2( reference_value - window_values.at(i) );
            
            if (my_del2 < del2)
            {
                del2 = my_del2;
                idx = i;
            }
        }
    }
    
    BOOST_ASSERT( idx != -1 );
    
    idx += step_size;
    
    if (idx < 0)
        idx = 0;
    else if (idx >= window_values.count())
        idx = window_values.count() - 1;

    return window_values.at(idx);
}

/** Return whether or not the changes in the passed
    delta *may* have changed the system since the last
    subversion 'subversion' */
bool WindowedComponent::mayChange(const Delta &delta, quint32 last_subversion) const
{
    BOOST_ASSERT( Constraint::wasLastSystem(delta.deltaSystem()) );

    return (not Constraint::wasLastSatisfied()) or
           delta.sinceChanged(reference_component, last_subversion) or
           delta.sinceChanged(constrained_component, last_subversion);
}

/** Set the baseline system for the constraint - this is 
    used to pre-calculate everything for the system
    and to check if the constraint is satisfied */
void WindowedComponent::setSystem(const System &system)
{
    if ( Constraint::wasLastSystem(system) and Constraint::wasLastSubVersion(system) )
        return;
    
    Constraint::clearLastSystem();
    
    if (window_values.isEmpty())
    {
        Constraint::setSatisfied(system, true);
        return;
    }
    
    has_constrained_value = false;
    constrained_value = 0;
    target_value = 0;
    
    if (system.hasConstantComponent(constrained_component))
    {
        has_constrained_value = true;
        constrained_value = system.constant(constrained_component);
    }
    
    component_val = system.constant(reference_component);
    
    if (window_values.count() == 1)
        target_value = window_values.at(0);
    else
        target_value = ::getNextWindow(component_val, window_values, step_size);
    
    Constraint::setSatisfied( system, has_constrained_value and
                                      SireMaths::areEqual(constrained_value,target_value) );
}

/** Fully apply this constraint on the passed delta - this returns
    whether or not this constraint affects the delta */
bool WindowedComponent::fullApply(Delta &delta)
{
    if (window_values.isEmpty())
        return false;

    BOOST_ASSERT( Constraint::wasLastSystem(delta.deltaSystem()) and
                  Constraint::wasLastSubVersion(delta.deltaSystem()) );
    
    constrained_value = target_value;
    
    bool changed = delta.update(constrained_component, target_value);
    
    if (changed)
    {
        this->setSystem(delta.deltaSystem());
    }
    
    return changed;
}

/** Apply this constraint based on the delta, knowing that the 
    last application of this constraint was on this system, 
    at subversion number last_subversion */
bool WindowedComponent::deltaApply(Delta &delta, quint32 last_subversion)
{
    if (window_values.isEmpty())
        return false;

    BOOST_ASSERT( Constraint::wasLastSystem(delta.deltaSystem()) );

    bool changed_comp = delta.sinceChanged(constrained_component, last_subversion);
    bool changed_sym = delta.sinceChanged(reference_component, last_subversion);
    bool changed_target = false;
    
    if (changed_comp or changed_sym or (constrained_value != target_value))
    {
        const System &system = delta.deltaSystem();
        
        if (changed_sym)
        {
            component_val = system.constant(reference_component);
            
            double new_target;

            if (window_values.count() == 1)
                new_target = window_values.at(0);
            else
                new_target = ::getNextWindow(component_val, window_values, step_size);            

            if (new_target != target_value)
            {
                target_value = new_target;
                changed_target = true;
            }
        }
        
        if (changed_comp)
        {
            if (system.hasConstantComponent(constrained_component))
            {
                double new_comp = system.constant(constrained_component);
                
                if ( (not has_constrained_value) or constrained_value != new_comp )
                {
                    has_constrained_value = true;
                    constrained_value = new_comp;
                    changed_comp = true;
                }
                else
                    changed_comp = false;
            }
            else
            {
                has_constrained_value = false;
                constrained_value = 0;
                changed_comp = true;
            }
        }

        if ( (changed_target or changed_comp) and not
             (has_constrained_value and (constrained_value == target_value)) )
        {
            bool changed = delta.update(constrained_component, target_value);
            constrained_value = delta.deltaSystem().constant(constrained_component);
            return changed;
        }
    }
    else if (not this->isSatisfied(delta.deltaSystem()))
    {
        return this->fullApply(delta);
    }

    return false;
}

const char* WindowedComponent::typeName()
{
    return QMetaType::typeName( qMetaTypeId<WindowedComponent>() );
}

NullConstraint* NullConstraint::clone() const
{
    return new NullConstraint(*this);
}

