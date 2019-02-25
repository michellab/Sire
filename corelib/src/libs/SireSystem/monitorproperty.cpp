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

#include <QFile>
#include <QTextStream>

#include "monitorproperty.h"

#include "SireMol/moleculegroup.h"

#include "SireFF/ff.h"

#include "system.h"

#include "SireError/errors.h"
#include "SireMol/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireSystem;
using namespace SireMol;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<MonitorProperty> r_monprop;

QDataStream &operator<<(QDataStream &ds, 
                                          const MonitorProperty &monprop)
{
    writeHeader(ds, r_monprop, 2);
    
    SharedDataStream sds(ds);
    
    sds << monprop.what_is_monitored;
    
    switch( monprop.what_is_monitored )
    {
        case MonitorProperty::SYSTEM_PROPERTY:
            sds << monprop.prop << monprop.props;
            break;
        case MonitorProperty::FORCEFIELD_PROPERTY:
            sds << monprop.prop << monprop.ffid << monprop.props;
            break;
        case MonitorProperty::MOLECULE_PROPERTY:
            sds << monprop.prop << monprop.mgid << monprop.molprops;
            break;
    }
    
    sds << static_cast<const SystemMonitor&>(monprop);
    
    return ds;
}

QDataStream &operator>>(QDataStream &ds,
                                          MonitorProperty &monprop)
{
    VersionID v = readHeader(ds, r_monprop);
    
    if (v == 2)
    {
        SharedDataStream sds(ds);
        
        MonitorProperty mon;
        sds >> mon.what_is_monitored;
        
        switch(mon.what_is_monitored)
        {
        case MonitorProperty::SYSTEM_PROPERTY:
            sds >> mon.prop >> mon.props;
            break;
        case MonitorProperty::FORCEFIELD_PROPERTY:
            sds >> mon.prop >> mon.ffid >> mon.props;
            break;
        case MonitorProperty::MOLECULE_PROPERTY:
            sds >> mon.prop >> mon.mgid >> mon.molprops;
            break;
        }
        
        sds >> static_cast<SystemMonitor&>(mon);
        
        monprop = mon;
    }
    else if (v == 1)
    {
        SharedDataStream sds(ds);
        
        MonitorProperty mon;
        sds >> mon.what_is_monitored;
        
        QVector<PropertyPtr> props;
        QHash< MolNum,QVector<PropertyPtr> > molprops;
        
        switch(mon.what_is_monitored)
        {
        case MonitorProperty::SYSTEM_PROPERTY:
            sds >> mon.prop >> props;
            mon.props = ChunkedVector<PropertyPtr,2048>::fromVector(props);
            break;
        case MonitorProperty::FORCEFIELD_PROPERTY:
            sds >> mon.prop >> mon.ffid >> props;
            mon.props = ChunkedVector<PropertyPtr,2048>::fromVector(props);
            break;
        case MonitorProperty::MOLECULE_PROPERTY:
            sds >> mon.prop >> mon.mgid >> molprops;
        
            mon.molprops.reserve(molprops.count());
        
            for (QHash< MolNum,QVector<PropertyPtr> >::const_iterator
                                                it = molprops.constBegin();
                 it != molprops.constEnd();
                 ++it)
            {
                mon.molprops.insert(it.key(), 
                                    ChunkedVector<PropertyPtr,2048>::fromVector(*it));
            }

            break;
        }
        
        sds >> static_cast<SystemMonitor&>(mon);
        
        monprop = mon;
    }
    else
        throw version_error(v, "1", r_monprop, CODELOC);
        
    return ds;
}

/** Null constructor */
MonitorProperty::MonitorProperty() 
                : ConcreteProperty<MonitorProperty,SystemMonitor>(),
                  what_is_monitored(IS_NULL)
{}

/** Construct to monitor the system property 'property' */
MonitorProperty::MonitorProperty(const QString &property)
                : ConcreteProperty<MonitorProperty,SystemMonitor>(),
                  prop(property), what_is_monitored(SYSTEM_PROPERTY)
{}

/** Construct to monitor the property 'property' of the molecule(s) contained
    in the molecule group 'molgroup' */
MonitorProperty::MonitorProperty(const QString &property, const MoleculeGroup &molgroup)
                : ConcreteProperty<MonitorProperty,SystemMonitor>(),
                  prop(property), mgid(molgroup.number()), 
                  what_is_monitored(MOLECULE_PROPERTY)
{}

/** Construct to monitor the property 'property' of the molecule(s) contained
    in the molecule group(s) that matches the ID 'mgid' */
MonitorProperty::MonitorProperty(const QString &property, const MGID &id)
                : ConcreteProperty<MonitorProperty,SystemMonitor>(),
                  prop(property), mgid(id), 
                  what_is_monitored(MOLECULE_PROPERTY)
{}

/** Construct to monitor the property 'property' of the passed forcefield */
MonitorProperty::MonitorProperty(const QString &property, const FF &forcefield)
                : ConcreteProperty<MonitorProperty,SystemMonitor>(),
                  prop(property), ffid(forcefield.name()),
                  what_is_monitored(FORCEFIELD_PROPERTY)
{}

/** Construct to monitor the property 'property' of the forcefield(s)
    that match the passed ID - note that the property must be the same
    in each of these forcefields at every step of the simulation */
MonitorProperty::MonitorProperty(const QString &property, const FFID &id)
                : ConcreteProperty<MonitorProperty,SystemMonitor>(),
                  prop(property), ffid(id), 
                  what_is_monitored(FORCEFIELD_PROPERTY)
{}

/** Copy constructor */
MonitorProperty::MonitorProperty(const MonitorProperty &other)
                : ConcreteProperty<MonitorProperty,SystemMonitor>(other),
                  prop(other.prop), mgid(other.mgid), ffid(other.ffid),
                  what_is_monitored(other.what_is_monitored),
                  props(other.props), molprops(other.molprops)
{}

/** Destructor */
MonitorProperty::~MonitorProperty()
{}

/** Copy assignment operator */
MonitorProperty& MonitorProperty::operator=(const MonitorProperty &other)
{
    if (this != &other)
    {
        SystemMonitor::operator=(other);
        prop = other.prop;
        mgid = other.mgid;
        ffid = other.ffid;
        what_is_monitored = other.what_is_monitored;
        props = other.props;
        molprops = other.molprops;
    }
    
    return *this;
}

const char* MonitorProperty::typeName()
{
    return QMetaType::typeName( qMetaTypeId<MonitorProperty>() );
}

/** Comparison operator */
bool MonitorProperty::operator==(const MonitorProperty &other) const
{
    return this == &other or
           (prop == other.prop and what_is_monitored == other.what_is_monitored and
            ffid == other.ffid and mgid == other.mgid and
            props == other.props and molprops == other.molprops);
}

/** Comparison operator */
bool MonitorProperty::operator!=(const MonitorProperty &other) const
{
    return not MonitorProperty::operator==(other);
}

/** Return whether or not this is monitoring a system property */
bool MonitorProperty::monitoringSystemProperty() const
{
    return what_is_monitored == SYSTEM_PROPERTY;
}

/** Return whether or not this is monitoring a forcefield property */
bool MonitorProperty::monitoringForceFieldProperty() const
{
    return what_is_monitored == FORCEFIELD_PROPERTY;
}

/** Return whether or not this is monitoring a molecule property
    of molecules in molecule groups */
bool MonitorProperty::monitoringMoleculeProperty() const
{
    return what_is_monitored == MOLECULE_PROPERTY;
}

/** Return the name of the property being monitored */
const QString& MonitorProperty::property() const
{
    return prop;
}

QString MonitorProperty::toString() const
{
    switch( what_is_monitored )
    {
    case IS_NULL:
        return QObject::tr("MonitorProperty::null");

    case SYSTEM_PROPERTY:
        return QObject::tr("MonitorProperty( system property \"%1\" )")
                                .arg(prop);

    case FORCEFIELD_PROPERTY:
        return QObject::tr("MonitorProperty( forcefield (%1) property \"%2\" )")
                    .arg(ffid.toString(), prop);
                    
    case MOLECULE_PROPERTY:
        return QObject::tr("MonitorProperty( molecule group (%1) property \"%2\" )")
                    .arg(mgid.toString(), prop);
    }
    
    return QObject::tr("MonitorProperty( BUGGY! )");
}

/** Return the ID of the molecule group(s) whose molecules
    are being monitored */
const MGID& MonitorProperty::mgID() const
{
    if (not monitoringMoleculeProperty())
        throw SireError::incompatible_error( QObject::tr(
                "This MonitorProperty (%1) is not monitoring the property "
                "of molecules in molecule group(s).")
                    .arg( this->toString() ), CODELOC );

    return mgid.base();
}

/** Return the ID of the forcefield(s) that are being monitored */
const FFID& MonitorProperty::ffID() const
{
    if (not monitoringForceFieldProperty())
        throw SireError::incompatible_error( QObject::tr(
                "This MonitorProperty (%1) is not monitoring the property "
                "of forcefield(s).")
                    .arg( this->toString() ), CODELOC );

    return ffid.base();
}

/** Clear all statistics */
void MonitorProperty::clearStatistics()
{
    props.clear();
    molprops.clear();
}

/** Return the values of the monitored system or forcefield properties */
QVector<PropertyPtr> MonitorProperty::properties() const
{
    if (not (monitoringSystemProperty() or monitoringForceFieldProperty()))
        throw SireError::incompatible_error( QObject::tr(
                "This MonitorProperty (%1) is not monitoring the property "
                "of the system or forcefield(s).")
                    .arg( this->toString() ), CODELOC );
    
    return props.toVector();
}

/** Return the values of the monitored molecule properties for molecule 'molnum'.

    \throw SireMol::missing_molecule
*/
QVector<PropertyPtr> MonitorProperty::properties(MolNum molnum) const
{
    if (not monitoringMoleculeProperty())
        throw SireError::incompatible_error( QObject::tr(
                "This MonitorProperty (%1) is not monitoring the property "
                "of molecules in molecule group(s).")
                    .arg( this->toString() ), CODELOC );

    QHash< MolNum,ChunkedVector<PropertyPtr,2048> >::const_iterator 
                                                    it = molprops.constFind(molnum);
    
    if (it == molprops.constEnd())
        throw SireMol::missing_molecule( QObject::tr(
                "No properties were monitored for the molecule with number %1. "
                "Use the \"monitoredMolecules()\" function to get the list of "
                "monitored molecule numbers.")
                    .arg(molnum.toString()), CODELOC );

    return it.value().toVector();
}

/** Return the numbers of molecules whose properties have been monitored */
QList<MolNum> MonitorProperty::monitoredMolecules() const
{
    if (not monitoringMoleculeProperty())
        throw SireError::incompatible_error( QObject::tr(
                "This MonitorProperty (%1) is not monitoring the property "
                "of molecules in molecule group(s).")
                    .arg( this->toString() ), CODELOC );

    return molprops.keys();
}

/** Write all of the properties to disk in text format to the 
    file 'filename'. The aim of this function is to provide something
    quick and dirty to follow the properties. To properly save the properties
    you should use the streaming functions */
void MonitorProperty::writeToDisk(const QString &filename)
{
    if (what_is_monitored == IS_NULL or (props.isEmpty() and molprops.isEmpty()))
        return;
        
    QFile f( filename );
    
    if (not f.open(QIODevice::WriteOnly))
        throw SireError::file_error(f, CODELOC);
        
    QTextStream ts(&f);
    
    switch (what_is_monitored)
    {
    case SYSTEM_PROPERTY:
        ts << QObject::tr("# System property \"%1\"\n").arg(prop);
        break;
        
    case FORCEFIELD_PROPERTY:
        ts << QObject::tr("# Forcefield property \"%1\" for %2\n")
                        .arg(prop, ffid.toString());
        break;
        
    case MOLECULE_PROPERTY:
        ts << QObject::tr("# Molecule property \"%1\" for %2\n")
                        .arg(prop, mgid.toString());
        break;
    }
    
    if (not props.isEmpty())
    {
        for (int i=0; i<props.count(); ++i)
        {
            ts << i << " : " << props.at(i).read().toString() << "\n";
        }
    }
    else if (not molprops.isEmpty())
    {
        for (QHash< MolNum,ChunkedVector<PropertyPtr,2048> >::const_iterator
                                    it = molprops.constBegin();
             it != molprops.constEnd();
             ++it)
        {
            ts << QObject::tr("\n ---- Molecule %1 ----\n")
                        .arg(it.key());
                        
            const ChunkedVector<PropertyPtr,2048> &p = *it;
            
            for (int i=0; i<p.count(); ++i)
            {
                ts << i << " : " << p.at(i).read().toString() << "\n";
            }
        }
    }
    
    f.close();
}

void MonitorProperty::monitor(const Molecules &molecules)
{
    int nsteps = 0;
    
    if (not molprops.isEmpty())
        nsteps = molprops.constBegin()->count();
    
    for (Molecules::const_iterator it = molecules.constBegin();
         it != molecules.constEnd();
         ++it)
    {
        if (nsteps > 0)
        {
            if (not molprops.contains(it.key()))
                molprops.insert(it.key(), ChunkedVector<PropertyPtr,2048>(nsteps));
        }

        const MoleculeData &moldata = it->data();
        
        PropertyPtr property;
        
        if (moldata.hasProperty(prop))
            property = moldata.property(prop);
            
        molprops[it.key()].append(property);
    }
}

/** Monitor the system */
void MonitorProperty::monitor(System &system)
{
    switch( what_is_monitored )
    {
    case IS_NULL:
        return;
        
    case SYSTEM_PROPERTY:
        if (system.containsProperty(prop))
            props.append( system.property(prop) );
        else
            props.append( PropertyPtr() );
            
        break;
        
    case FORCEFIELD_PROPERTY:
        if (system.containsProperty(ffid, prop))
            props.append( system.property(ffid, prop) );
        else
            props.append( PropertyPtr() );
            
        break;
        
    case MOLECULE_PROPERTY:
        this->monitor( system.molecules(mgid) );
        break;
    }
}
