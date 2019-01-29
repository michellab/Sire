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

#ifndef SIRESYSTEM_MONITORPROPERTY_H
#define SIRESYSTEM_MONITORPROPERTY_H

#include "systemmonitor.h"

#include "SireFF/ffidentifier.h"

#include "SireMol/mgidentifier.h"
#include "SireMol/molnum.h"

#include "SireBase/chunkedvector.hpp"

SIRE_BEGIN_HEADER

namespace SireSystem
{
class MonitorProperty;
}

SIRESYSTEM_EXPORT QDataStream& operator<<(QDataStream&, const SireSystem::MonitorProperty&);
SIRESYSTEM_EXPORT QDataStream& operator>>(QDataStream&, SireSystem::MonitorProperty&);

namespace SireFF
{
class FF;
}

namespace SireMol
{
class MoleculeGroup;
}

namespace SireSystem 
{

using SireFF::FFID;
using SireFF::FFIdentifier;
using SireFF::FF;

using SireMol::MoleculeGroup;
using SireMol::MGID;
using SireMol::MGIdentifier;
using SireMol::MolNum;
using SireMol::Molecules;

using SireBase::PropertyPtr;

/** This monitor is used to monitor the value of system, forcefield or
    molecule properties during a simulation
    
    @author Christopher Woods
*/
class SIRESYSTEM_EXPORT MonitorProperty 
        : public SireBase::ConcreteProperty<MonitorProperty,SystemMonitor>
{

friend QDataStream& ::operator<<(QDataStream&, const MonitorProperty&);
friend QDataStream& ::operator>>(QDataStream&, MonitorProperty&);

public:
    MonitorProperty();
    MonitorProperty(const QString &property);
    MonitorProperty(const QString &property, const MoleculeGroup &molgroup);
    MonitorProperty(const QString &property, const MGID &mgid);
    MonitorProperty(const QString &property, const FF &forcefield);
    MonitorProperty(const QString &property, const FFID &ffid);
                     
    MonitorProperty(const MonitorProperty &other);
    
    ~MonitorProperty();
    
    MonitorProperty& operator=(const MonitorProperty &other);
    
    static const char* typeName();
    
    bool operator==(const MonitorProperty &other) const;
    bool operator!=(const MonitorProperty &other) const;
    
    QString toString() const;
    
    bool monitoringSystemProperty() const;
    bool monitoringForceFieldProperty() const;
    bool monitoringMoleculeProperty() const;
    
    const QString& property() const;

    const MGID& mgID() const;
    const FFID& ffID() const;

    void clearStatistics();

    QVector<PropertyPtr> properties() const;
    QVector<PropertyPtr> properties(MolNum molnum) const;

    QList<MolNum> monitoredMolecules() const;

    void writeToDisk(const QString &filename);

    void monitor(System &system);

private:
    void monitor(const Molecules &molecules);

    /** The name of the property being monitored */
    QString prop;
    
    /** The ID of the molecule group(s) whose molecules 
        are being monitored */
    MGIdentifier mgid;
    
    /** The ID of the forcefield(s) whose properties are being monitored */
    FFIdentifier ffid;
    
    enum { IS_NULL = 0,
           SYSTEM_PROPERTY = 1,
           MOLECULE_PROPERTY = 2,
           FORCEFIELD_PROPERTY = 3 };
           
    /** What is being monitored */
    quint32 what_is_monitored;
    
    /** The system or forcefield properties */
    SireBase::ChunkedVector<PropertyPtr,2048> props;
    
    /** The molecule properties */
    QHash< MolNum,SireBase::ChunkedVector<PropertyPtr,2048> > molprops;
};

}

Q_DECLARE_METATYPE( SireSystem::MonitorProperty )

SIRE_EXPOSE_CLASS( SireSystem::MonitorProperty )

SIRE_END_HEADER

#endif
