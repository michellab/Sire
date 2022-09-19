/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2006  Christopher Woods
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

#ifndef SIRESYSTEM_SYSTEM_H
#define SIRESYSTEM_SYSTEM_H

#include <QUuid>

#include "sysname.h"
#include "systemmonitors.h"
#include "constraints.h"

#include "SireVol/space.h"

#include "SireMol/moleculegroup.h"
#include "SireMol/moleculegroups.h"
#include "SireMol/mgnum.h"

#include "SireFF/forcefields.h"

#include "SireBase/majorminorversion.h"

SIRE_BEGIN_HEADER

namespace SireSystem
{
class System;
}

SIRESYSTEM_EXPORT QDataStream& operator<<(QDataStream&, const SireSystem::System&);
SIRESYSTEM_EXPORT QDataStream& operator>>(QDataStream&, SireSystem::System&);

namespace SireSystem
{

using SireFF::FFPtr;
using SireFF::ForceFields;
using SireFF::FF;
using SireFF::FFID;
using SireFF::FFIdx;
using SireFF::FFName;
using SireFF::EnergyTable;
using SireFF::ForceTable;
using SireFF::FieldTable;
using SireFF::PotentialTable;
using SireFF::Probe;

using SireMol::MolGroupsBase;
using SireMol::MoleculeGroup;
using SireMol::MolGroupsPtr;
using SireMol::MoleculeGroups;
using SireMol::MGNum;
using SireMol::MGIdx;
using SireMol::MGName;
using SireMol::MGID;
using SireMol::MolID;
using SireMol::MolNum;
using SireMol::MoleculeData;
using SireMol::MoleculeView;
using SireMol::ViewsOfMol;
using SireMol::Molecules;

using SireVol::Space;
using SireVol::Space;

using SireCAS::Symbol;
using SireCAS::Symbols;
using SireCAS::Expression;
using SireCAS::Values;

using SireBase::Property;
using SireBase::PropertyPtr;
using SireBase::Properties;
using SireBase::PropertyMap;
using SireBase::PropertyName;
using SireBase::MajorMinorVersion;
using SireBase::Version;

/** This is a simulation system. If contains molecules, forcefields that
    provide energy functions of those molecules, and monitors that
    can monitor the changing state of the system

    @author Christopher Woods
*/
class SIRESYSTEM_EXPORT System
            : public SireBase::ConcreteProperty<System,MolGroupsBase>
{

friend SIRESYSTEM_EXPORT QDataStream& ::operator<<(QDataStream&, const System&);
friend SIRESYSTEM_EXPORT QDataStream& ::operator>>(QDataStream&, System&);

public:
    System();
    System(const QString &name);

    System(const System &other);

    ~System();

    static const char* typeName();

    System& operator=(const System &other);

    bool operator==(const System &other) const;
    bool operator!=(const System &other) const;

    using SireMol::MolGroupsBase::operator[];

    const FF& operator[](const FFID &ffid) const;
    const SystemMonitor& operator[](const MonitorID &monid) const;
    const MoleculeGroup& operator[](const MGID &mgid) const;

    ViewsOfMol operator[](int i) const;
    ViewsOfMol operator[](const QString &name) const;
    QList<SireMol::MolViewPtr> operator[](const SireBase::Slice &slice) const;

    ViewsOfMol operator[](MolNum molnum) const;
    ViewsOfMol operator[](const MolID &molid) const;

    System& operator+=(const FF &forcefield);
    System& operator+=(const MoleculeGroup &molgroup);
    System& operator+=(const Constraint &constraint);
    System& operator+=(const Constraints &constraints);

    System& operator-=(const FF &forcefield);
    System& operator-=(const MoleculeGroup &molgroup);
    System& operator-=(const FFID &ffid);
    System& operator-=(const MGID &mgid);
    System& operator-=(const MolID &molid);
    System& operator-=(const Constraint &constraint);
    System& operator-=(const Constraints &constraints);

    const QUuid& UID() const;
    const SysName& name() const;
    const Version& version() const;
    quint32 subVersion() const;

    void setName(const SysName &newname);
    void setName(const QString &newname);

    void collectStats();

    void applyConstraints();
    bool constraintsSatisfied() const;

    using SireMol::MolGroupsBase::at;

    const FF& at(const FFID &ffid) const;
    const SystemMonitor& at(const MonitorID &monid) const;

    const FF& forceField(const FFID &ffid) const;
    const FF& forceField(const MGID &mgid) const;

    const SystemMonitor& monitor(const MonitorID &monid) const;

    QList<SysMonPtr> monitors(const MonitorID &monid) const;

    int nForceFields() const;
    int nMonitors() const;
    int nConstraints() const;

    FFIdx ffIdx(const FFID &ffid) const;

    const FFName& ffName(const FFID &ffid) const;

    MonitorName monitorName(const MonitorID &monid) const;

    QString toString() const;

    const Symbol& totalComponent() const;

    SireUnits::Dimension::MolarEnergy energy();
    SireUnits::Dimension::MolarEnergy energy(const Symbol &component);

    Values energies();
    Values energies(const QSet<Symbol> &components);

    bool isEnergyComponent(const Symbol &component) const;
    bool hasEnergyComponent(const Symbol &component) const;

    void setEnergyComponent(const Symbol &symbol,
                            const SireCAS::Expression &expression);

    QSet<Symbol> energySymbols() const;
    Values energyComponents();

    SireCAS::Expression energyExpression(const Symbol &expression) const;
    QHash<Symbol,SireCAS::Expression> energyExpressions(
                                            const QSet<Symbol> &symbols) const;
    QHash<Symbol,SireCAS::Expression> energyExpressions() const;

    double constant(const Symbol &component) const;

    Values constants() const;
    Values constants(const QSet<Symbol> &components) const;

    bool isConstantComponent(const Symbol &component) const;
    bool hasConstantComponent(const Symbol &component) const;

    void setConstant(const Symbol &symbol, double value);
    void setConstant(const Symbol &symbol, const SireCAS::Expression &expression);

    void setConstantComponent(const Symbol &symbol, double value);
    void setConstantComponent(const Symbol &symbol,
                              const SireCAS::Expression &expression);

    QSet<Symbol> constantSymbols() const;
    Values constantComponents() const;

    SireCAS::Expression constantExpression(const Symbol &symbol) const;
    QHash<Symbol,SireCAS::Expression> constantExpressions(
                                            const QSet<Symbol> &symbols) const;
    QHash<Symbol,SireCAS::Expression> constantExpressions() const;

    void setComponent(const Symbol &symbol, double value);
    void setComponent(const Symbol &symbol, const SireCAS::Expression &expression);

    QSet<Symbol> componentSymbols() const;
    Values componentValues();
    Values componentValues(const QSet<Symbol> &symbols);

    bool hasComponent(const Symbol &symbol) const;
    double componentValue(const Symbol &symbol);

    SireCAS::Expression componentExpression(const Symbol &symbol) const;
    QHash<Symbol,SireCAS::Expression> componentExpressions(
                                            const QSet<Symbol> &symbols) const;
    QHash<Symbol,SireCAS::Expression> componentExpressions() const;

    void energy(EnergyTable &energytable, double scale_energy=1);
    void energy(EnergyTable &energytable, const Symbol &component, double scale_energy=1);

    void force(ForceTable &forcetable, double scale_force=1);
    void force(ForceTable &forcetable, const Symbol &component,
               double scale_force=1);

    void field(FieldTable &fieldtable, double scale_field=1);
    void field(FieldTable &fieldtable, const Symbol &component,
               double scale_field=1);

    void field(FieldTable &fieldtable, const Probe &probe, double scale_field=1);
    void field(FieldTable &fieldtable, const Symbol &component,
               const Probe &probe, double scale_field=1);

    void potential(PotentialTable &pottable, const Probe &probe,
                   double scale_potential=1);
    void potential(PotentialTable &pottable, const Symbol &component,
                   const Probe &probe, double scale_potential=1);

    void potential(PotentialTable &pottable, double scale_potential=1);
    void potential(PotentialTable &pottable, const Symbol &component,
                   double scale_potential=1);

    void setProperty(const QString &name, const Property &value);
    void setProperty(const FFID &ffid, const QString &name, const Property &value);

    void removeProperty(const QString &name);

    bool isCompoundProperty(const QString &name) const;
    bool isUserProperty(const QString &name) const;
    bool isBuiltinProperty(const QString &name) const;

    const Property& compoundProperty(const QString &name) const;
    const Property& userProperty(const QString &name) const;
    const Property& builtinProperty(const QString &name) const;

    const Property& property(const PropertyName &name) const;

    const Property& property(const FFID &ffid, const PropertyName &name) const;

    QStringList propertyKeys() const;
    QStringList propertyKeys(const FFID &ffid) const;

    bool containsProperty(const QString &name) const;
    bool containsProperty(const FFID &ffid, const QString &name) const;

    bool containsProperty(const PropertyName &name) const;
    bool containsProperty(const FFID &ffid, const PropertyName &name) const;

    Properties properties() const;
    Properties properties(const FFID &ffid) const;

    Properties userProperties() const;
    Properties builtinProperties() const;

    const SystemMonitors& monitors() const;
    const ForceFields& forceFields() const;
    const MoleculeGroups& extraGroups() const;
    const Constraints& constraints() const;

    void clearStatistics();
    void clearStatistics(const MonitorID &monid);

    void mustNowRecalculateFromScratch();

    void accept();
    bool needsAccepting() const;

    bool isDirty() const;
    bool isClean() const;

    using SireMol::MolGroupsBase::add;
    using SireMol::MolGroupsBase::remove;
    using SireMol::MolGroupsBase::update;

    void add(const QString &name, const SystemMonitor &monitor,
             int frequency = 1);

    void add(const SystemMonitors &monitors);
    void add(const SystemMonitors &monitors, int frequency);

    void setMonitors(const SystemMonitors &monitors);
    void setMonitors(const SystemMonitors &monitors, int frequency);

    void add(const FF &forcefield);
    void add(const MoleculeGroup &molgroup);

    void add(const Constraint &constraint);
    void add(const Constraints &constraints);

    void setConstraints(const Constraints &constraints);

    void remove(const MonitorID &monid);

    void remove(const FFID &ffid);
    void remove(const FF &ff);

    bool remove(const MGID &mgid);

    bool remove(const MoleculeGroup &molgroup);
    bool remove(const MolID &molid);

    void remove(const Constraint &constraint);
    void remove(const Constraints &constraints);

    bool removeAllMolecules();

    void removeAllMoleculeGroups();
    void removeAllForceFields();
    void removeAllMonitors();
    void removeAllConstraints();

    //overloading MolGroupsBase virtual functions
    const MoleculeGroup& at(MGNum mgnum) const;

    void add(const MoleculeView &molview, const MGID &mgid,
             const PropertyMap &map);
    void add(const ViewsOfMol &molviews, const MGID &mgid,
             const PropertyMap &map);
    void add(const Molecules &molecules, const MGID &mgid,
             const PropertyMap &map);
    void add(const MoleculeGroup &molgroup, const MGID &mgid,
             const PropertyMap &map);

    void addIfUnique(const MoleculeView &molview, const MGID &mgid,
                     const PropertyMap &map);
    void addIfUnique(const ViewsOfMol &molviews, const MGID &mgid,
                     const PropertyMap &map);
    void addIfUnique(const Molecules &molecules, const MGID &mgid,
                     const PropertyMap &map);
    void addIfUnique(const MoleculeGroup &molgroup, const MGID &mgid,
                     const PropertyMap &map);

    void add(const MoleculeView &molview, const MGID &mgid);
    void add(const ViewsOfMol &molviews, const MGID &mgid);
    void add(const Molecules &molecules, const MGID &mgid);
    void add(const MoleculeGroup &molgroup, const MGID &mgid);

    void addIfUnique(const MoleculeView &molview, const MGID &mgid);
    void addIfUnique(const ViewsOfMol &molviews, const MGID &mgid);
    void addIfUnique(const Molecules &molecules, const MGID &mgid);
    void addIfUnique(const MoleculeGroup &molgroup, const MGID &mgid);

    using MolGroupsBase::removeAll;

    bool removeAll(const MGID &mgid);
    bool remove(const MoleculeView &molview, const MGID &mgid);
    bool remove(const ViewsOfMol &molviews, const MGID &mgid);
    bool remove(const Molecules &molecules, const MGID &mgid);
    bool remove(const MoleculeGroup &molgroup, const MGID &mgid);

    bool removeAll(const MoleculeView &molview, const MGID &mgid);
    bool removeAll(const ViewsOfMol &molviews, const MGID &mgid);
    bool removeAll(const Molecules &molecules, const MGID &mgid);
    bool removeAll(const MoleculeGroup &molgroup, const MGID &mgid);

    bool remove(MolNum molnum, const MGID &mgid);
    bool remove(const QSet<MolNum> &molnums, const MGID &mgid);

    void update(const MoleculeData &moldata, bool auto_commit=true);
    void update(const Molecules &molecules, bool auto_commit=true);
    void update(const MoleculeGroup &molgroup, bool auto_commit=true);

    void setContents(const MGID &mgid, const MoleculeView &molview,
                     const PropertyMap &map);
    void setContents(const MGID &mgid, const ViewsOfMol &molviews,
                     const PropertyMap &map);
    void setContents(const MGID &mgid, const Molecules &molecules,
                     const PropertyMap &map);
    void setContents(const MGID &mgid, const MoleculeGroup &molgroup,
                     const PropertyMap &map);

    void setContents(const MGID &mgid, const MoleculeView &molview);
    void setContents(const MGID &mgid, const ViewsOfMol &molviews);
    void setContents(const MGID &mgid, const Molecules &molecules);
    void setContents(const MGID &mgid, const MoleculeGroup &molgroup);

    int nFrames() const;
    int nFrames(const SireBase::PropertyMap &map) const;

    void loadFrame(int frame);
    void saveFrame(int frame);
    void saveFrame();
    void deleteFrame(int frame);

    void loadFrame(int frame, const SireBase::PropertyMap &map);
    void saveFrame(int frame, const SireBase::PropertyMap &map);
    void saveFrame(const SireBase::PropertyMap &map);
    void deleteFrame(int frame, const SireBase::PropertyMap &map);

    static const System& null();

protected:
    const MoleculeGroup& getGroup(MGNum mgnum) const;

    void getGroups(const QList<MGNum> &mgnums,
                   QVarLengthArray<const MoleculeGroup*,10> &groups) const;

    QHash<MGNum,const MoleculeGroup*> getGroups() const;

    void reindex();

    friend class Delta; // so can call below functions
    bool deltaUpdate(const MoleculeData &moldata, bool auto_commit);
    QList<MolNum> deltaUpdate(const Molecules &molecules, bool auto_commit);

    void applyAllConstraints();

    bool deltaUpdate(const Symbol &component, double value);
    bool deltaUpdate(const QString &property, const Property &value);
    bool deltaUpdate(const QString &property, const FFID &ffid,
                     const Property &value);
    bool deltaUpdate(const QString &property, const QList<FFIdx> &ffidxs,
                     const Property &value);

    void commitDelta(const Constraints &constraints,
                     bool is_minor_change,
                     bool is_major_change);

private:
    void rebuildIndex();

    ForceFields& _pvt_forceFields();

    const ForceFields& _pvt_constForceFields() const;
    const ForceFields& _pvt_forceFields() const;

    MoleculeGroups& _pvt_moleculeGroups();

    const MoleculeGroups& _pvt_constMoleculeGroups() const;
    const MoleculeGroups& _pvt_moleculeGroups() const;

    MolGroupsBase& _pvt_moleculeGroups(MGNum mgnum);

    const MolGroupsBase& _pvt_moleculeGroups(MGNum mgnum) const;
    const MolGroupsBase& _pvt_constMoleculeGroups(MGNum mgnum) const;

    const MoleculeGroup& _pvt_moleculeGroup(MGNum mgnum) const;

    void _pvt_throwMissingGroup(MGNum mgnum) const;

    void _pvt_applyMoleculeConstraints();
    void _pvt_applyMoleculeConstraints(MolNum molnum);
    void _pvt_applyMoleculeConstraints(const Molecules &molecules);

    /** The unique ID for this system */
    QUuid uid;

    /** The name of this system */
    SysName sysname;

    /** The version number of this system */
    MajorMinorVersion sysversion;

    /** The molecule groups in this system. These are divided
        into two - the first set are actually the forcefields,
        while the second are the non-forcefield groups */
    MolGroupsPtr molgroups[2];

    /** All of the monitors that monitor this system */
    SystemMonitors sysmonitors;

    /** All of the constraints that are applied to this system */
    Constraints cons;

    /** The index of which of the two set of MoleculeGroups each
        individual molecule group in this set is in */
    QHash<MGNum,int> mgroups_by_num;

    /** The subversion of this system - this is incremented when
        delta updates are being applied. A system with non-zero
        subversion is not guaranteed to be in a valid state */
    quint32 subversion;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Return the unique ID of this system */
SIRE_ALWAYS_INLINE const QUuid& System::UID() const
{
    return uid;
}

/** Return the name of this system */
SIRE_ALWAYS_INLINE const SysName& System::name() const
{
    return sysname;
}

/** Return a reference to the version of this system */
SIRE_ALWAYS_INLINE const Version& System::version() const
{
    return sysversion.version();
}

/** Return the subversion number of this system */
SIRE_ALWAYS_INLINE quint32 System::subVersion() const
{
    return subversion;
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

Q_DECLARE_METATYPE( SireSystem::System )

SIRE_EXPOSE_CLASS( SireSystem::System )

SIRE_END_HEADER

#endif
