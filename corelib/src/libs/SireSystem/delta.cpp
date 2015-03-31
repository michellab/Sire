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

#include "delta.h"
#include "constraints.h"

#include "SireFF/point.h"

#include "SireError/errors.h"
#include "SireBase/errors.h"

using namespace SireSystem;
using namespace SireMol;
using namespace SireFF;
using namespace SireBase;
using namespace SireCAS;
using namespace SireStream;

/** Null constructor */
Delta::Delta() : last_change(0), last_mol_change(0),
                 last_comp_change(0), last_prop_change(0), auto_commit(true)
{}

/** Construct to begin applying to the passed system */
Delta::Delta(const System &system, bool autoc)
      : delta_system(system), last_change(0),
        last_mol_change(0), last_comp_change(0),
        last_prop_change(0), auto_commit(autoc)
{
    if (system.subVersion() != 0)
        throw SireError::program_bug( QObject::tr(
                "It is a mistake to construct a delta for a System (%1) "
                "that is already undergoing a staged change!")
                    .arg(system.toString()), CODELOC );
}

/** Copy constructor */
Delta::Delta(const Delta &other)
      : delta_system(other.delta_system), changed_mols(other.changed_mols),
        changed_comps(other.changed_comps), changed_props(other.changed_props),
        last_change(other.last_change), last_mol_change(other.last_mol_change),
        last_comp_change(other.last_comp_change), 
        last_prop_change(other.last_prop_change),
        auto_commit(other.auto_commit)
{}

/** Destructor */
Delta::~Delta()
{}

/** Copy assignment operator */
Delta& Delta::operator=(const Delta &other)
{
    if (this != &other)
    {
        delta_system = other.delta_system;
        changed_mols = other.changed_mols;
        changed_comps = other.changed_comps;
        changed_props = other.changed_props;
        last_change = other.last_change;
        last_mol_change = other.last_mol_change;
        last_comp_change = other.last_comp_change;
        last_prop_change = other.last_prop_change;
        auto_commit = other.auto_commit;
    }
    
    return *this;
}

/** Comparison operator */
bool Delta::operator==(const Delta &other) const
{
    return this == &other or
           (delta_system == other.delta_system and
            changed_mols == other.changed_mols and
            changed_comps == other.changed_comps and
            changed_props == other.changed_props and
            last_change == other.last_change and
            last_mol_change == other.last_mol_change and
            last_comp_change == other.last_comp_change and
            last_prop_change == other.last_prop_change and
            auto_commit == other.auto_commit);
}

/** Comparison operator */
bool Delta::operator!=(const Delta &other) const
{
    return not Delta::operator==(other);
}

/** Return a string representation of this delta */
QString Delta::toString() const
{
    return QObject::tr("Delta( %1, last change %2 [ %3, %4, %5 ] )")
                .arg(delta_system.toString())
                .arg(last_change)
                .arg(last_mol_change)
                .arg(last_comp_change)
                .arg(last_prop_change);
}

/** Return whether or not this delta will auto-commit each update */
bool Delta::willAutoCommit() const
{
    return auto_commit;
}

/** Return whether or not this is an empty delta */
bool Delta::isEmpty() const
{
    return delta_system.UID().isNull();
}

/** Return whether or not this is a null delta */
bool Delta::isNull() const
{
    return isEmpty();
}

/** Return the system in its current in-between state. 
    Note that this system is not likely to be in a sensible
    state! */
const System& Delta::deltaSystem() const
{
    return delta_system;
}

/** Return whether or not this delta represents a change */
bool Delta::hasChange() const
{
    return last_change != 0;
}

/** Return whether or not this delta represents a major
    change (addition or removal of molecules, change
    of system components or properties) */
bool Delta::hasMajorChange() const
{
    return last_comp_change != 0 or last_prop_change != 0;
}

/** Return whether or not this delta represents a minor
    change (possibly together with a major change), e.g.
    a change of molecules only */
bool Delta::hasMinorChange() const
{
    return last_mol_change != 0;
}

/** Return whether or not this delta has a change of molecule */
bool Delta::hasMoleculeChange() const
{
    return last_mol_change != 0;
}

/** Return whether or not this delta has a change in system
    components */
bool Delta::hasComponentChange() const
{
    return last_comp_change != 0;
}

/** Return whether or not this delta has a change in 
    system property */
bool Delta::hasPropertyChange() const
{
    return last_prop_change != 0;
}

/** Return whether or not this delta has a change in 
    the system since subversion 'subversion' */
bool Delta::hasChangeSince(quint32 subversion) const
{
    return last_change > subversion;
}

/** Return whether or not this delta has a major change in 
    the system since subversion 'subversion' */
bool Delta::hasMajorChangeSince(quint32 subversion) const
{
    return last_comp_change > subversion or last_prop_change > subversion;
}

/** Return whether or not this delta has a minor change in 
    the system since subversion 'subversion' */
bool Delta::hasMinorChangeSince(quint32 subversion) const
{
    return last_mol_change > subversion;
}

/** Return whether or not this delta has a change in 
    molecules since subversion 'subversion' */
bool Delta::hasMoleculeChangeSince(quint32 subversion) const
{
    return last_mol_change > subversion;
}

/** Return whether or not this delta has a change in 
    the system components since subversion 'subversion' */
bool Delta::hasComponentChangeSince(quint32 subversion) const
{
    return last_comp_change > subversion;
}

/** Return whether or not this delta has a change in 
    the system properties since subversion 'subversion' */
bool Delta::hasPropertyChangeSince(quint32 subversion) const
{
    return last_prop_change > subversion;
}

/** Return whether or not this delta has a change in 
    the molecule with number 'molnum' */
bool Delta::changed(MolNum molnum) const
{
    return changed_mols.contains(molnum);
}

/** Return whether or not the molecule viewed in 'molview' has
    changed */
bool Delta::changed(const MoleculeView &molview) const
{
    return changed_mols.contains(molview.data().number());
}

/** Return whether or not any of the molecules in 'molecules' 
    have been changed */
bool Delta::changed(const Molecules &molecules) const
{
    if (changed_mols.isEmpty())
        return false;

    if (molecules.count() <= changed_mols.count())
    {
        for (Molecules::const_iterator it = molecules.constBegin();
             it != molecules.constEnd();
             ++it)
        {
            if (changed_mols.contains(it.key()))
                return true;
        }
        
        return false;
    }
    else
    {
        for (QHash<MolNum,quint32>::const_iterator it = changed_mols.constBegin();
             it != changed_mols.constEnd();
             ++it)
        {
            Molecules::const_iterator it2 = molecules.constFind(it.key());
            
            if (it2 != molecules.constEnd())
                return true;
        }
        
        return false;
    }
}

/** Return whether or not the component 'component' has changed */
bool Delta::changed(const Symbol &component) const
{
    return changed_comps.contains(component);
}

/** Return whether or not any of the components whose symbols are
    in 'components' have been changed */
bool Delta::changed(const QSet<Symbol> &components) const
{
    if (changed_comps.isEmpty())
        return false;

    else if (components.count() <= changed_comps.count())
    {
        foreach (const Symbol &component, components)
        {
            if (changed_comps.contains(component))
                return true;
        }
        
        return false;
    }
    else
    {
        for (QHash<Symbol,quint32>::const_iterator it = changed_comps.constBegin();
             it != changed_comps.constEnd();
             ++it)
        {
            if (components.contains(it.key()))
                return true;
        }
        
        return false;
    }
}

/** Return whether or not any of the components whose values
    are in 'values' have changed */
bool Delta::changed(const Values &values) const
{
    if (changed_comps.isEmpty())
        return false;
    
    else if (values.count() <= changed_comps.count())
    {
        for (Values::const_iterator it = values.constBegin();
             it != values.constEnd();
             ++it)
        {
            if (changed_comps.contains( Symbol(it.key()) ))
                return true;
        }
        
        return false;
    }
    else
    {
        for (QHash<Symbol,quint32>::const_iterator it = changed_comps.constBegin();
             it != changed_comps.constEnd();
             ++it)
        {
            if (values.contains(it.key()))
                return true;
        }
        
        return false;
    }
}

/** Return whether or not the system property 'property' has changed */
bool Delta::changed(const QString &property) const
{
    return changed_props.contains(property);
}

/** Return whether or not the system property 'property' has changed */
bool Delta::changed(const PropertyName &property) const
{
    if (property.hasValue())
        return false;
    else 
        return this->changed(property.source());
}

/** Return if any of the system properties in 'properties' have changed */
bool Delta::changed(const QSet<QString> &properties) const
{
    if (changed_props.isEmpty())
        return false;
        
    else if (properties.count() <= changed_props.count())
    {
        foreach (const QString &property, properties)
        {
            if (changed_props.contains(property))
                return true;
        }
        
        return false;
    }
    else
    {
        for (QHash<QString,quint32>::const_iterator it = changed_props.constBegin();
             it != changed_props.constEnd();
             ++it)
        {
            if (properties.contains(it.key()))
                return true;
        }
        
        return false;
    }
}

/** Return if any of the system properties in 'properties' have changed */
bool Delta::changed(const QList<PropertyName> &properties) const
{
    if (changed_props.isEmpty())
        return false;
        
    else
    {
        foreach (const PropertyName &property, properties)
        {
            if ( (not property.hasValue()) and changed_props.contains(property.source()) )
                return true;
        }
        
        return false;
    }
}

/** Return if any of the properties in 'properties' have changed */
bool Delta::changed(const Properties &properties) const
{
    if (changed_props.isEmpty())
        return false;
        
    else if (properties.count() <= changed_props.count())
    {
        for (Properties::const_iterator it = properties.constBegin();
             it != properties.constEnd();
             ++it)
        {
            if (changed_props.contains(it.key()))
                return true;
        }
        
        return false;
    }
    else
    {
        for (QHash<QString,quint32>::const_iterator it = changed_props.constBegin();
             it != changed_props.constEnd();
             ++it)
        {
            if (properties.hasProperty(it.key()))
                return true;
        }
        
        return false;
    }
}

/** Return whether or not the passed point will have been changed by this delta */
bool Delta::changed(const Point &point) const
{
    if (point.nMolecules() == 0)
        return false;
    
    else
    {
        for (QHash<MolNum,quint32>::const_iterator it = changed_mols.constBegin();
             it != changed_mols.constEnd();
             ++it)
        {
            if (point.contains(it.key()))
                return true;
        }
        
        return false;
    }
}

/** Return whether or not the molecule with number 'molnum' has changed
    since subversion number 'subversion' */
bool Delta::sinceChanged(MolNum molnum, quint32 subversion) const
{
    if (subversion == 0)
        return this->changed(molnum);

    else if (changed_mols.isEmpty() or subversion <= last_mol_change)
        return false;
    
    else
        return changed_mols.value(molnum, 0) > subversion;
}

/** Return whether or not the molecule viewed in 'molview' has changed
    since subversion number 'subversion' */
bool Delta::sinceChanged(const MoleculeView &molview, quint32 subversion) const
{
    return this->sinceChanged(molview.data().number(), subversion);
}

/** Return whether or not any of the molecules in 'molecules' have
    changed since subversion number 'subversion' */
bool Delta::sinceChanged(const Molecules &molecules, quint32 subversion) const
{
    if (subversion == 0)
        return this->changed(molecules);
    
    else if (changed_mols.isEmpty() or subversion <= last_mol_change)
        return false;
    
    else if (molecules.nMolecules() <= changed_mols.count())
    {
        for (Molecules::const_iterator it = molecules.constBegin();
             it != molecules.constEnd();
             ++it)
        {
            if (changed_mols.value(it.key(), 0) > subversion)
                return true;
        }
        
        return false;
    }
    else
    {
        for (QHash<MolNum,quint32>::const_iterator it = changed_mols.constBegin();
             it != changed_mols.constEnd();
             ++it)
        {
            if (it.value() > subversion and molecules.contains(it.key()))
                return true;
        }
        
        return false;
    }
}

/** Return whether or not the constant system component 'component' has
    changed since subversion number 'subversion' */
bool Delta::sinceChanged(const Symbol &component, quint32 subversion) const
{
    if (subversion == 0)
        return this->changed(component);

    else if (changed_comps.isEmpty() or last_comp_change <= subversion)
        return false;
    
    else
        return changed_comps.value(component, 0) > subversion;
}

/** Return whether or not any of the system constant components whose
    symbols are in 'symbols' have changed since subversion number 'subversion' */
bool Delta::sinceChanged(const QSet<Symbol> &symbols, quint32 subversion) const
{
    if (subversion == 0)
        return this->changed(symbols);

    else if (changed_comps.isEmpty() or last_comp_change <= subversion)
        return false;
        
    else if (symbols.count() <= changed_comps.count())
    {
        foreach (const Symbol &symbol, symbols)
        {
            if (changed_comps.value(symbol, 0) > subversion)
                return true;
        }
        
        return false;
    }
    else
    {
        for (QHash<Symbol,quint32>::const_iterator it = changed_comps.constBegin(); 
             it != changed_comps.constEnd();
             ++it)
        {
            if (it.value() > subversion and symbols.contains(it.key()))
                return true;
        }
        
        return false;
    }
}

/** Return whether or not the system constant components from 'values' 
    have changed since subversion number 'subversion' */
bool Delta::sinceChanged(const Values &values, quint32 subversion) const
{
    if (subversion == 0)
        return this->changed(values);
    
    else if (changed_comps.isEmpty() or last_comp_change <= subversion)
        return false;
    
    else if (values.count() <= changed_comps.count())
    {
        for (Values::const_iterator it = values.constBegin();
             it != values.constEnd();
             ++it)
        {
            if (changed_comps.value(Symbol(it.key()), 0) > subversion)
                return true;
        }
        
        return false;
    }
    else
    {
        for (QHash<Symbol,quint32>::const_iterator it = changed_comps.constBegin();
             it != changed_comps.constEnd();
             ++it)
        {
            if (it.value() > subversion and values.contains(it.key()))
                return true;
        }
        
        return false;
    }
}

/** Return whether or not the property 'property' has changed since 
    subversion number 'subversion' */
bool Delta::sinceChanged(const QString &property, quint32 subversion) const
{
    if (subversion == 0)
        return this->changed(property);

    else if (changed_props.isEmpty() or last_prop_change <= subversion)
        return false;
    
    else
        return changed_props.value(property, 0) > subversion;
}

/** Return whether or not the property 'property' has changed since 
    subversion number 'subversion' */
bool Delta::sinceChanged(const PropertyName &property, quint32 subversion) const
{
    if (property.hasValue())
        return false;
    else
        return this->sinceChanged(property.source(), subversion);
}

/** Return whether or not any of the properties in 'properties' have changed
    since subversion number 'subversion' */
bool Delta::sinceChanged(const QSet<QString> &properties, quint32 subversion) const
{
    if (subversion == 0)
        return this->changed(properties);

    else if (changed_props.isEmpty() or last_prop_change <= subversion)
        return false;

    else if (properties.count() <= changed_props.count())
    {
        foreach (const QString &property, properties)
        {
            if (changed_props.value(property, 0) > subversion)
                return true;
        }
        
        return false;
    }
    else
    {
        for (QHash<QString,quint32>::const_iterator it = changed_props.constBegin();
             it != changed_props.constEnd();
             ++it)
        {
            if (it.value() > subversion and properties.contains(it.key()))
                return true;
        }
        
        return false;
    }
}

/** Return whether or not any of the properties in 'properties' have changed
    since subversion number 'subversion' */
bool Delta::sinceChanged(const QList<PropertyName> &properties, quint32 subversion) const
{
    if (subversion == 0)
        return this->changed(properties);

    else if (changed_props.isEmpty() or last_prop_change <= subversion)
        return false;

    else
    {
        foreach (const PropertyName &property, properties)
        {
            if ( (not property.hasValue()) and 
                        changed_props.value(property.source(), 0) > subversion )
            {
                return true;
            }
        }
        
        return false;
    }
}

/** Return whether or not any of the properties in 'properties' have changed
    since subversion number 'subversion' */
bool Delta::sinceChanged(const Properties &properties, quint32 subversion) const
{
    if (subversion == 0)
        return this->changed(properties);

    else if (changed_props.isEmpty() or last_prop_change <= subversion)
        return false;
        
    else if (properties.count() <= changed_props.count())
    {
        for (Properties::const_iterator it = properties.constBegin();
             it != properties.constEnd();
             ++it)
        {
            if (changed_props.value(it.key(), 0) > subversion)
                return true;
        }
        
        return false;
    }
    else
    {
        for (QHash<QString,quint32>::const_iterator it = changed_props.constBegin();
             it != changed_props.constEnd();
             ++it)
        {
            if (it.value() > subversion and properties.hasProperty(it.key()))
                return true;
        }
        
        return false;
    }
}

/** Return whether or not the passed point will have been changed by this delta 
    since version 'subversion' */
bool Delta::sinceChanged(const Point &point, quint32 subversion) const
{
    if (point.nMolecules() == 0 or last_mol_change <= subversion)
        return false;
    
    else if (subversion == 0)
        return this->changed(point);
    
    else
    {
        for (QHash<MolNum,quint32>::const_iterator it = changed_mols.constBegin();
             it != changed_mols.constEnd();
             ++it)
        {
            if (it.value() > subversion and point.contains(it.key()))
                return true;
        }
        
        return false;
    }
}
    
/** Return the numbers of the molecules that have changed in this delta */
QList<MolNum> Delta::changedMolecules() const
{
    return changed_mols.keys();
}

/** Return the numbers of the molecules that have changed 
    in this delta since subversion 'subversion' */
QList<MolNum> Delta::changedMoleculesSince(quint32 subversion) const
{
    if (subversion == 0)
        return this->changedMolecules();

    QList<MolNum> molnums;
    
    if (last_mol_change > subversion)
    {
        for (QHash<MolNum,quint32>::const_iterator it = changed_mols.constBegin();
             it != changed_mols.constEnd();
             ++it)
        {
            if (it.value() > subversion)
                molnums.append( it.key() );
        }
    }
    
    return molnums;
}

/** Return the numbers of the molecules from 'molecules' that have 
    been changed by this delta */
QList<MolNum> Delta::changedMolecules(const Molecules &molecules) const
{
    if (changed_mols.isEmpty() or molecules.isEmpty())
        return QList<MolNum>();
    
    else if (changed_mols.count() <= molecules.count())
    {
        QList<MolNum> molnums;
        
        for (QHash<MolNum,quint32>::const_iterator it = changed_mols.constBegin();
             it != changed_mols.constEnd();
             ++it)
        {
            if (molecules.contains(it.key()))
                molnums.append(it.key());
        }
        
        return molnums;
    }
    else
    {
        QList<MolNum> molnums;
        
        for (Molecules::const_iterator it = molecules.constBegin();
             it != molecules.constEnd();
             ++it)
        {
            if (changed_mols.contains(it.key()))
                molnums.append(it.key());
        }
        
        return molnums;
    }
}

/** Return the numbers of the molecules from 'molecules' that have
    changed in this delta since subversion 'subversion' */
QList<MolNum> Delta::changedMoleculesSince(const Molecules &molecules,
                                           quint32 subversion) const
{
    if (subversion == 0)
        return this->changedMolecules(molecules);
    
    QList<MolNum> molnums;
    
    if (last_mol_change > subversion)
    {
        if (changed_mols.count() <= molecules.count())
        {
            for (QHash<MolNum,quint32>::const_iterator it = changed_mols.constBegin();
                 it != changed_mols.constEnd();
                 ++it)
            {
                if (it.value() > subversion and molecules.contains(it.key()))
                    molnums.append(it.key());
            }
        }
        else
        {
            for (Molecules::const_iterator it = molecules.constBegin();
                 it != molecules.constEnd();
                 ++it)
            {
                if (changed_mols.value(it.key(), 0) > subversion)
                    molnums.append(it.key());
            }
        }
    }

    return molnums;
}

/** Return the symbols of all changed system components */
QList<Symbol> Delta::changedComponents() const
{
    return changed_comps.keys();
}

/** Return the symbols of all changed system components in
    this delta since subversion 'subversion' */
QList<Symbol> Delta::changedComponentsSince(quint32 subversion) const
{
    if (subversion == 0)
        return this->changedComponents();
    
    QList<Symbol> comps;
    
    if (last_comp_change > subversion)
    {
        for (QHash<Symbol,quint32>::const_iterator it = changed_comps.constBegin(); 
             it != changed_comps.constEnd();
             ++it)
        {
            if (it.value() > subversion)
                comps.append( it.key() );
        }
    }
    
    return comps;
}

/** Return all of the changed system properties in this delta */
QList<QString> Delta::changedProperties() const
{
    return changed_props.keys();
}

/** Return all of the changed system properties in this delta since
    subversion 'subversion' */
QList<QString> Delta::changedPropertiesSince(quint32 subversion) const
{
    if (subversion == 0)
        return this->changedProperties();
    
    QList<QString> props;
    
    if (last_prop_change > subversion)
    {
        for (QHash<QString,quint32>::const_iterator it = changed_props.constBegin();
             it != changed_props.constEnd();
             ++it)
        {
            if (it.value() > subversion)
                props.append( it.key() );
        }
    }
    
    return props;
}

/** Update the contained system to match the molecule version
    in 'moldata' */
bool Delta::update(const MoleculeData &moldata)
{
    if (delta_system.deltaUpdate(moldata,auto_commit))
    {
        last_change = delta_system.subVersion();
        last_mol_change = last_change;
        changed_mols.insert(moldata.number(), last_mol_change);
        
        return true;
    }
    else
        return false;
}

/** Update the contained system to match the version of the molecule
    contained in 'molview' */
bool Delta::update(const MoleculeView &molview)
{
    if (delta_system.deltaUpdate(molview.data(),auto_commit))
    {
        last_change = delta_system.subVersion();
        last_mol_change = last_change;
        changed_mols.insert(molview.data().number(), last_mol_change);
        
        return true;
    }
    else
        return false;
}

/** Update the contained system to match the versions of the molecules
    contained in 'molecules' */
bool Delta::update(const Molecules &molecules)
{
    QList<MolNum> molnums = delta_system.deltaUpdate(molecules,auto_commit);
    
    if (not molnums.isEmpty())
    {
        last_change = delta_system.subVersion();
        last_mol_change = last_change;
        
        if (changed_mols.count() < molnums.count())
            changed_mols.reserve(molnums.count());
        
        foreach (const MolNum &molnum, molnums)
        {
            changed_mols.insert(molnum, last_mol_change);
        }
        
        return true;
    }
    else
        return false;
}

/** Update the contained system to set the value of the component
    'component' to the value 'value' */
bool Delta::update(const Symbol &component, double value)
{
    if (delta_system.constant(component) != value)
    {
        if (delta_system.deltaUpdate(component, value))
        {
            last_change = delta_system.subVersion();
            last_comp_change = last_change;
            
            changed_comps.insert(component, last_comp_change);
            
            return true;
        }
    }
    
    return false;
}

/** Update the contained system to set the value of the property
    'property' to 'value' */
bool Delta::update(const QString &property, const Property &value)
{
    if (delta_system.deltaUpdate(property,value))
    {
        last_change = delta_system.subVersion();
        last_prop_change = last_change;
        
        changed_props.insert(property, last_prop_change);
        
        return true;
    }
    else
        return false;
}

/** Update the contained system to set the value of the property
    'property' to 'value' */
bool Delta::update(const PropertyName &property, const Property &value)
{
    if (property.hasValue())
        return false;
    else
        return this->update(property.source(), value);
}

/** Update the contained system to set the value of the property
    in the forcefield(s) matching 'ffid' to the value 'value' */
bool Delta::update(const QString &property, const FFID &ffid, const Property &value)
{
    if (delta_system.deltaUpdate(property, ffid, value))
    {
        last_change = delta_system.subVersion();
        last_prop_change = last_change;
        
        changed_props.insert(property, last_prop_change);
        
        return true;
    }
    else
        return false;
}

/** Update the contained system to set the value of the property
    in the forcefields whose indicies are in 'ffidxs' to the value 'value' */
bool Delta::update(const QString &property, const QList<FFIdx> &ffidxs,
                   const Property &value)
{
    bool changed_prop = false;

    foreach (const FFIdx &ffidx, ffidxs)
    {
        if (delta_system.deltaUpdate(property, ffidx, value))
        {
            last_change = delta_system.subVersion();
            last_prop_change = last_change;
        
            changed_props.insert(property, last_prop_change);
            changed_prop = true;
        }
    }
    
    return changed_prop;
}

/** Update the contained system to set the value of the property
    in the forcefield(s) matching 'ffid' to the value 'value' */
bool Delta::update(const PropertyName &property, const FFID &ffid, const Property &value)
{
    if (property.hasValue())
        return false;
    else
        return this->update(property.source(), ffid, value);
}
/** Update the contained system to set the value of the property
    in the forcefields whose indicies are in 'ffidxs' to the value 'value' */
bool Delta::update(const PropertyName &property, const QList<FFIdx> &ffidxs,
                   const Property &value)
{
    if (property.hasValue())
        return false;
    else
        return this->update(property.source(), ffidxs, value);
}

/** Apply the delta - this resolves all of the constraints, cleans
    up all of the changes and returns a valid system that incorporates
    everything in this delta */
System Delta::apply()
{
    Constraints constraints = delta_system.constraints();

    //apply these constraints
    constraints.apply(*this);
    delta_system.commitDelta(constraints, hasMinorChange(), hasMajorChange());

    if (last_change != 0)
    {
        last_change = 0;
        last_mol_change = 0;
        last_comp_change = 0;
        last_prop_change = 0;
    
        changed_mols = QHash<MolNum,quint32>();
        changed_comps = QHash<Symbol,quint32>();
        changed_props = QHash<QString,quint32>();
    }
    
    return delta_system;
}
