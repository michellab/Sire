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

#include "constraints.h"
#include "moleculeconstraint.h"
#include "system.h"
#include "delta.h"

#include "SireMol/molecules.h"

#include "SireID/index.h"

#include "SireBase/savestate.h"

#include "SireSystem/errors.h" 
#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>

using namespace SireSystem;
using namespace SireMol;
using namespace SireID;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<Constraints> r_constraints;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, 
                                          const Constraints &constraints)
{
    writeHeader(ds, r_constraints, 3);
    
    SharedDataStream sds(ds);
    
    sds << constraints.cons
        << static_cast<const Property&>(constraints);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, Constraints &constraints)
{
    VersionID v = readHeader(ds, r_constraints);
    
    if (v == 3)
    {
        SharedDataStream sds(ds);
        
        sds >> constraints.cons >> static_cast<Property&>(constraints);
    }
    else if (v == 2)
    {
        SharedDataStream sds(ds);
        
        QVector<ConstraintPtr> molcons;
        
        sds >> constraints.cons >> molcons
            >> static_cast<Property&>(constraints);
            
        constraints.cons += molcons;
    }
    else if (v == 1)
    {
        SharedDataStream sds(ds);

        QList<ConstraintPtr> cons;

        sds >> cons;
        
        constraints.cons = cons.toVector();
    }
    else
        throw version_error(v, "1,2,3", r_constraints, CODELOC);
        
    return ds;
}

/** Null constructor */
Constraints::Constraints() : ConcreteProperty<Constraints,Property>()
{}

/** Construct to contain just the constraint 'constraint' */
Constraints::Constraints(const Constraint &constraint)
            : ConcreteProperty<Constraints,Property>()
{
    this->add(constraint);
}

/** Construct from the passed list of constraints */
Constraints::Constraints(const QVector<ConstraintPtr> &constraints)
            : ConcreteProperty<Constraints,Property>()
{
    foreach (const ConstraintPtr constraint, constraints)
    {
        cons.append(constraint);
    }
}

/** Construct from the passed list of constraints */
Constraints::Constraints(const QList<ConstraintPtr> &constraints)
            : ConcreteProperty<Constraints,Property>()
{
    foreach (const ConstraintPtr constraint, constraints)
    {
        cons.append(constraint);
    }
}

/** Copy constructor */
Constraints::Constraints(const Constraints &other)
            : ConcreteProperty<Constraints,Property>(other),
              cons(other.cons)
{}

/** Destructor */
Constraints::~Constraints()
{}

/** Copy assignment operator */
Constraints& Constraints::operator=(const Constraints &other)
{
    if (this != &other)
    {
        cons = other.cons;
        Property::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool Constraints::operator==(const Constraints &other) const
{
    return cons == other.cons;
}

/** Comparison operator */
bool Constraints::operator!=(const Constraints &other) const
{
    return not this->operator==(other);
}

/** Return the ith constraint in this set

    \throw SireError::invalid_index
*/
const Constraint& Constraints::operator[](int i) const
{
    i = Index(i).map( this->nConstraints() );
    
    return cons.at( i ).read();
}

/** Syntactic sugar for Constraints::add */
Constraints& Constraints::operator+=(const Constraint &constraint)
{
    this->add(constraint);
    return *this;
}

/** Syntactic sugar for Constraints::add */
Constraints& Constraints::operator+=(const Constraints &constraints)
{
    this->add(constraints);
    return *this;
}

/** Syntactic sugar for Constraints::remove */
Constraints& Constraints::operator-=(const Constraint &constraint)
{
    this->remove(constraint);
    return *this;
}

/** Syntactic sugar for Constraints::remove */
Constraints& Constraints::operator-=(const Constraints &constraints)
{
    this->remove(constraints);
    return *this;
}

/** Return the number of constraints in this set */
int Constraints::nConstraints() const
{
    return cons.count();
}

/** Return the number of constraints in this set */
int Constraints::count() const
{
    return this->nConstraints();
}

/** Return the number of constraints in this set */
int Constraints::size() const
{
    return this->nConstraints();
}

/** Return whether this is empty (contains no constraints) */
bool Constraints::isEmpty() const
{
    return cons.isEmpty();
}

/** Return the list of all of the constraints in this set */
QVector<ConstraintPtr> Constraints::constraints() const
{
    return cons;
}

/** Add the passed constraint to this set - this is only added
    if it does not exist in this set */
void Constraints::add(const Constraint &constraint)
{
    foreach (const ConstraintPtr con, cons)
    {
        if (con->equals(constraint))
           return;
    }
        
    cons.append(constraint);
    cons.squeeze();
}

/** Add all of the passed constraints to this set. This only
    adds the constraints that are not already part of this set */
void Constraints::add(const Constraints &constraints)
{
    for (QVector<ConstraintPtr>::const_iterator it = constraints.cons.constBegin();
         it != constraints.cons.constEnd();
         ++it)
    {
        this->add(it->read());
    }
}

/** Remove the constraint 'constraint' from this set - this
    does nothing if this constraint is not part of this set */
void Constraints::remove(const Constraint &constraint)
{
    QMutableVectorIterator<ConstraintPtr> it(cons);

    while (it.hasNext())
    {
        it.next();
    
        if (it.value()->equals(constraint))
            it.remove();
    }

    cons.squeeze();
}

/** Remove all of the constraints in 'constraints' from this
    set - this ignores constraints that are not in this set */
void Constraints::remove(const Constraints &constraints)
{
    for (QVector<ConstraintPtr>::const_iterator it = constraints.cons.constBegin();
         it != constraints.cons.constEnd();
         ++it)
    {
        this->remove( it->read() );
    }
}

/** Remove the constraint at index 'i' */
void Constraints::removeAt(int i)
{
    i = Index(i).map( this->nConstraints() );
    
    cons.remove(i);
    cons.squeeze();
}

/** Return whether or not all of the constraints in this set are
    satisfied in the passed system */
bool Constraints::areSatisfied(const System &system) const
{
    for (QVector<ConstraintPtr>::const_iterator it = cons.constBegin();
         it != cons.constEnd();
         ++it)
    {
        if (not it->read().isSatisfied(system))
        {
            return false;
        }
    }
    
    return true;
}

/** Assert that all of the constraints in this set are satisfied
    in the passed system
    
    \throw SireSystem::constraint_error
*/
void Constraints::assertSatisfied(const System &system) const
{
    QStringList broken_constraints;
    
    for (QVector<ConstraintPtr>::const_iterator it = cons.constBegin();
         it != cons.constEnd();
         ++it)
    {
        if (not it->read().isSatisfied(system))
        {
            broken_constraints.append( QObject::tr("%1 : %2")
                        .arg(broken_constraints.count() + 1)
                        .arg(it->read().toString()) );
        }
    }
    
    if (not broken_constraints.isEmpty())
    {
        throw SireSystem::constraint_error( QObject::tr(
            "Some of the constraints are not satisfied in the system %1. "
            "The number of unsatisfied constraints in %2. Here they are;\n%3")
                .arg(system.name())
                .arg(broken_constraints.count())
                .arg(broken_constraints.join("\n")), CODELOC );
    }
}

/** Apply all of the constraints to the passed system. This
    returns a system that satisfies all of the constraints */
System Constraints::apply(const System &system)
{
    if (cons.isEmpty())
        return system;

    Delta delta(system);
    
    System new_system(system);
    
    for (int i=0; i<10; ++i)
    {
        if (this->apply(delta))
            new_system = delta.apply();
        else
            return new_system;
    }
    
    throw SireSystem::constraint_error( QObject::tr(
            "The constraints %1 cannot be satisfied in connection with "
            "the constraints in the system %2, %3.")
                .arg(this->toString(), system.toString(), 
                     system.constraints().toString()), CODELOC );
                     
    return System();
}

/** Internal function used to apply all of the constraints in this system to 
    the passed delta - this returns whether or not any of the constraints
    changed the system */
bool Constraints::apply(Delta &delta)
{
    if (cons.isEmpty())
        return false;
        
    else if (cons.count() == 1)
    {
        if (cons.at(0).read().mayAffect(delta))
            return cons[0].edit().apply(delta);
        else
            return false;
    }
    else
    {
        bool system_changed = false;
        
        for (int i=0; i<10; ++i)
        {
            bool something_changed = false;
        
            for (int j=0; j<cons.count(); ++j)
            {
                if (cons.at(j).read().mayAffect(delta))
                {
                    bool this_changed = cons[j].edit().apply(delta);
                    something_changed = something_changed or this_changed;

                    system_changed = system_changed or this_changed;
                }
            }
            
            if (not something_changed)
            {
                bool all_satisfied = true;
            
                for (int j=0; j<cons.count(); ++j)
                {
                    if (not cons.at(j).read().isSatisfied(delta.deltaSystem()))
                    {
                        qDebug() << "Constraint" << cons.at(j).read().toString()
                                 << "not satisfied. Need to iterate constraints...";
                        qDebug() << Sire::toString(delta.deltaSystem().constants());
                        all_satisfied = false;
                        break;
                    }
                }
            
                if (all_satisfied)
                    return system_changed;
            }
        }

        bool all_satisfied = true;

        for (int j=0; j<cons.count(); ++j)
        {
            if (not cons.at(j).read().isSatisfied(delta.deltaSystem()))
            {
                 all_satisfied = false;
                 break;
            }
        }

        if (all_satisfied)
            return system_changed;
    
        //the constraints couldn't be satisfied - get the list
        //of unsatisfied constraints
        QStringList unsatisfied_constraints;
        
        for (int j=0; j<cons.count(); ++j)
        {
            if (not cons.at(j).read().isSatisfied(delta.deltaSystem()))
            {
                unsatisfied_constraints.append(cons[j].read().toString());
            }
        }
        
        throw SireSystem::constraint_error( QObject::tr(
                "Cannot simultaneously satisfy the following constraints "
                "on the system %1\n%2")
                    .arg(delta.deltaSystem().toString(), 
                         unsatisfied_constraints.join("\n")),
                            CODELOC );
    
        return true;
    }
}

void Constraints::committed(const System &system)
{
    if (cons.isEmpty())
        return;
    
    else
    {
        for (int i=0; i<cons.count(); ++i)
        {
            cons[i].edit().committed(system);
        }
    }
}

const char* Constraints::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Constraints>() );
}
