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

#ifndef SIREFF_FORCEFIELDS_H
#define SIREFF_FORCEFIELDS_H

#include "forcefield.h"

#include <boost/shared_ptr.hpp>

SIRE_BEGIN_HEADER

namespace SireFF
{
class ForceFields;
}

SIREFF_EXPORT QDataStream& operator<<(QDataStream&, const SireFF::ForceFields&);
SIREFF_EXPORT QDataStream& operator>>(QDataStream&, SireFF::ForceFields&);

namespace SireFF
{

class EnergyTable;
class ForceTable;
class FieldTable;
class PotentialTable;
class Probe;

using SireMol::MoleculeGroup;
using SireMol::MGNum;
using SireMol::MGID;
using SireMol::Molecules;
using SireMol::ViewsOfMol;
using SireMol::MoleculeView;

using SireBase::PropertyPtr;
using SireBase::PropertyName;

using SireCAS::Symbols;

namespace detail
{
class FFSymbol;
typedef boost::shared_ptr<FFSymbol> FFSymbolPtr;
}

/** A ForceFields object contains a collection of forcefields,
    together with energy functions that allow energies/forces
    of a combination of forcefield components to be evaluated.
    
    @author Christopher Woods
*/
class SIREFF_EXPORT ForceFields 
        : public SireBase::ConcreteProperty<ForceFields,SireMol::MolGroupsBase>
{

friend QDataStream& ::operator<<(QDataStream&, const ForceFields&);
friend QDataStream& ::operator>>(QDataStream&, ForceFields&);

public:
    ForceFields();
    ForceFields(const FF& forcefield);
    ForceFields(const QList<FFPtr> &forcefields);
    ForceFields(const QVector<FFPtr> &forcefields);
    
    ForceFields(const ForceFields &other);
    
    ~ForceFields();
    
    static const char* typeName();

    ForceFields& operator=(const ForceFields &other);
    
    bool operator==(const ForceFields &other) const;
    bool operator!=(const ForceFields &other) const;
    
    using SireMol::MolGroupsBase::operator[];
    using SireMol::MolGroupsBase::at;
    
    const FF& operator[](const FFName &ffname) const;
    const FF& operator[](const FFIdx &ffidx) const;
    const FF& operator[](const FFID &ffid) const;
    
    const FF& at(const FFName &ffname) const;
    const FF& at(const FFIdx &ffidx) const;
    const FF& at(const FFID &ffid) const;
    
    const FF& forceField(const FFName &ffname) const;
    const FF& forceField(const FFIdx &ffidx) const;
    const FF& forceField(const FFID &ffid) const;
    const FF& forceField(const MGNum &mgnum) const;
    
    int nForceFields() const;
    
    FFIdx ffIdx(const FFName &ffname) const;
    FFIdx ffIdx(const FFIdx &ffidx) const;
    FFIdx ffIdx(const FFID &ffid) const;
    
    const FFName& ffName(const FFName &ffname) const;
    const FFName& ffName(const FFIdx &ffidx) const;
    const FFName& ffName(const FFID &ffid) const;
    
    QList<FFIdx> map(const FFID &ffid) const;
    QList<FFIdx> map(const FFIdx &ffidx) const;
    QList<FFIdx> map(const FFName &ffname) const;
    
    QString toString() const;
    
    static const Symbol& totalComponent();
    
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
    
    void setConstantComponent(const Symbol &symbol, double value);
    void setConstantComponent(const Symbol &symbol,
                              const SireCAS::Expression &expression);

    QSet<Symbol> constantSymbols() const;
    Values constantComponents() const;
    
    double componentValue(const Symbol &symbol);
    Values componentValues(const QSet<Symbol> &symbols);
    Values componentValues();
    
    SireCAS::Expression constantExpression(const Symbol &symbol) const;
    QHash<Symbol,SireCAS::Expression> constantExpressions(
                                            const QSet<Symbol> &symbols) const;
    QHash<Symbol,SireCAS::Expression> constantExpressions() const;
    
    void setComponent(const Symbol &symbol, double value);
    void setComponent(const Symbol &symbol, const SireCAS::Expression &expression);
    
    QSet<Symbol> componentSymbols() const;
    
    bool hasComponent(const Symbol &symbol) const;
    
    SireCAS::Expression componentExpression(const Symbol &symbol) const;
    QHash<Symbol,SireCAS::Expression> componentExpressions(
                                            const QSet<Symbol> &symbols) const;
    QHash<Symbol,SireCAS::Expression> componentExpressions() const;
    
    void energy(EnergyTable &energytable, double scale_energy=1);
    void energy(EnergyTable &energytable, const Symbol &component,
		double scale_energy=1);

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
    
    const Property& property(const QString &name) const;

    const Property& property(const FFID &ffid, const QString &name) const;
 
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
    
    QVector<FFPtr> forceFieldsWithProperty(const QString &name) const;
    QVector<FFPtr> forceFieldsWithProperty(const FFID &ffid, const QString &name) const;
    
    const QVector<FFPtr>& forceFields() const;
    QList<FFName> ffNames() const;
    
    QVector<FFPtr> forceFields(const FFID &ffid) const;
    QList<FFName> ffNames(const FFID &ffid) const;
    
    const QVector<FFPtr>& list() const;
    QList<FFName> names() const;
    
    void mustNowRecalculateFromScratch();
    
    bool isDirty() const;
    bool isClean() const;
    
    using SireMol::MolGroupsBase::add;
    using SireMol::MolGroupsBase::remove;
    using SireMol::MolGroupsBase::update;
    
    void add(const FF &forcefield);
    
    void remove(const FFIdx &ffidx);
    void remove(const FFName &ffname);
    void remove(const FFID &ffid);

    void removeAllForceFields();

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

    bool needsAccepting() const;
    void accept();

protected:
    const MoleculeGroup& getGroup(MGNum mgnum) const;
    
    void getGroups(const QList<MGNum> &mgnums,
                   QVarLengthArray<const MoleculeGroup*,10> &groups) const;

    QHash<MGNum,const MoleculeGroup*> getGroups() const;

    void reindex();

private:
    void rebuildIndex();

    void sanitiseUserProperties();

    void getDependencies(const QString &name, QSet<QString> &deps) const;
    void assertNonCircularProperty(const QString &name) const;
    
    void assertValidLink(const QString &name) const;
    bool isValidLink(const QString &name) const;
    
    void updateCombinedProperty(const QString &name, bool update_dependencies=false,
                                QSet<QString> *updated_combinations=0);
                                
    void updateCombinedProperties();

    const FF& _pvt_forceField(int idx) const;
    FF& _pvt_forceField(int idx);

    const FF& _pvt_forceField(const FFName &ffname) const;
    FF& _pvt_forceField(const FFName &ffname);

    const FF& _pvt_forceField(const MGNum &mgnum) const;
    FF& _pvt_forceField(const MGNum &mgnum);

    void _pvt_remove(int i);

    /** The global symbol used to refer to the total energy of a collection
        of forcefields */
    static Symbol total_component;

    /** All of the forcefields arranged by FFIdx */
    QVector<FFPtr> ffields_by_idx;
    
    /** Map from forcefield name to its index */
    QHash<QString,int> ffields_by_name;
    
    /** Map from molecule group number to the name of the 
        forcefield that contains that group */
    QHash<MGNum, FFName> mgroups_by_num;
    
    /** All of the energy components, expressions and constants
        available in this collection of forcefields */
    QHash<Symbol, detail::FFSymbolPtr> ffsymbols;
    
    /** All of the additional normal properties. These are
        used to store additional information about the forcefields */
    QHash<QString, PropertyPtr> additional_properties;
    
    /** All of the property aliases, indexed by the property name
        they alias */
    QHash<QString, PropertyPtr> property_aliases;
    
    /** All of the combined properties, indexed by their name */
    QHash<QString, PropertyPtr> combined_properties;
};

}

Q_DECLARE_METATYPE( SireFF::ForceFields )

SIRE_EXPOSE_CLASS( SireFF::ForceFields )

SIRE_END_HEADER

#endif
