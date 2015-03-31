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

#ifndef SIRESYSTEM_DELTA_H
#define SIRESYSTEM_DELTA_H

#include "constraint.h"
#include "system.h"

#include "SireMol/molecule.h"
#include "SireMol/molecules.h"

#include "SireCAS/values.h"

SIRE_BEGIN_HEADER

namespace SireFF
{
class Point;
}

namespace SireSystem
{

using SireBase::Property;
using SireBase::Properties;

using SireCAS::Symbol;
using SireCAS::Values;

using SireMol::Molecule;
using SireMol::Molecules;
using SireMol::MolNum;
using SireMol::MoleculeView;

using SireFF::FFIdx;
using SireFF::FFID;
using SireFF::Point;

/** This class records and can apply a change from one system
    state to another. The purpose is to allow the state of 
    the system to be changed in small steps, and then
    committed as a chunk. This allows changes to be aggregated
    before constaints are applied or major/minor version
    numbers are incremented
    
    @author Christopher Woods
*/
class SIRESYSTEM_EXPORT Delta
{
public:
    Delta();
    Delta(const System &system, bool auto_commit=true);
          
    Delta(const Delta &other);
    
    ~Delta();
    
    Delta& operator=(const Delta &other);
    
    bool operator==(const Delta &other) const;
    bool operator!=(const Delta &other) const;

    QString toString() const;

    bool isEmpty() const;
    bool isNull() const;

    const System& deltaSystem() const;
    
    bool hasChange() const;
    
    bool hasMajorChange() const;
    bool hasMinorChange() const;
    
    bool hasMoleculeChange() const;
    bool hasComponentChange() const;
    bool hasPropertyChange() const;
    
    bool hasChangeSince(quint32 subversion) const;
    
    bool hasMajorChangeSince(quint32 subversion) const;
    bool hasMinorChangeSince(quint32 subversion) const;
    
    bool hasMoleculeChangeSince(quint32 subversion) const;
    bool hasComponentChangeSince(quint32 subversion) const;
    bool hasPropertyChangeSince(quint32 subversion) const;
    
    bool changed(MolNum molnum) const;
    bool changed(const MoleculeView &molview) const;
    bool changed(const Molecules &molecules) const;
    
    bool changed(const Symbol &component) const;
    bool changed(const QSet<Symbol> &symbols) const;
    bool changed(const Values &values) const;
    
    bool changed(const QString &property) const;
    bool changed(const PropertyName &property) const;
    bool changed(const QSet<QString> &properties) const;
    bool changed(const QList<PropertyName> &properties) const;
    bool changed(const Properties &properties) const;
    
    bool changed(const Point &point) const;
    
    bool sinceChanged(MolNum molnum, quint32 subversion) const;
    bool sinceChanged(const MoleculeView &molview, quint32 subversion) const;
    bool sinceChanged(const Molecules &molecules, quint32 subversion) const;

    bool sinceChanged(const Symbol &component, quint32 subversion) const;
    bool sinceChanged(const QSet<Symbol> &symbols, quint32 subversion) const;
    bool sinceChanged(const Values &values, quint32 subversion) const;

    bool sinceChanged(const QString &property, quint32 subversion) const;
    bool sinceChanged(const PropertyName &property, quint32 subversion) const;
    bool sinceChanged(const QSet<QString> &properties, quint32 subversion) const;
    bool sinceChanged(const QList<PropertyName> &properties, quint32 subversion) const;
    bool sinceChanged(const Properties &properties, quint32 subversion) const;
    
    bool sinceChanged(const Point &point, quint32 subversion) const;
    
    QList<MolNum> changedMolecules() const;
    QList<MolNum> changedMoleculesSince(quint32 subversion) const;
    
    QList<MolNum> changedMolecules(const Molecules &molecules) const;
    QList<MolNum> changedMoleculesSince(const Molecules &molecules,
                                        quint32 subversion) const;
    
    QList<Symbol> changedComponents() const;
    QList<Symbol> changedComponentsSince(quint32 subversion) const;
    
    QList<QString> changedProperties() const;
    QList<QString> changedPropertiesSince(quint32 subversion) const;
    
    bool update(const MoleculeView &molview);
    bool update(const MoleculeData &moldata);
    bool update(const Molecules &molecules);

    bool update(const Symbol &component, double value);

    bool update(const QString &property, const Property &value);
    bool update(const QString &property, const FFID &ffid, const Property &value);
    bool update(const QString &property, const QList<FFIdx> &ffidxs,
                const Property &value);
    
    bool update(const PropertyName &property, const Property &value);
    bool update(const PropertyName &property, const FFID &ffid,
                const Property &value);
    bool update(const PropertyName &property, const QList<FFIdx> &ffidxs,
                const Property &value);
    
    bool willAutoCommit() const;
    
    System apply();
    
private:
    /** The system being changed */
    System delta_system;
    
    /** Numbers of all of the changed molecules, together with 
        the subversion number when they were changed */
    QHash<MolNum,quint32> changed_mols;
    
    /** Symbols of all of the changed components, together with
        the subversion number when they were changed */
    QHash<Symbol,quint32> changed_comps;
    
    /** Names of all of the changed properties, together with
        the subversion number when they were changed */
    QHash<QString,quint32> changed_props;
    
    /** The subversion of the last change (0 for no change) */
    quint32 last_change;
    
    /** The subversion of the last molecule change (0 for no change) */
    quint32 last_mol_change;
    
    /** The subversion of the last component change (0 for no change) */
    quint32 last_comp_change;
    
    /** The subversion of the last property change (0 for no change) */
    quint32 last_prop_change;
    
    /** Whether or not to auto-commit after each update */
    bool auto_commit;
};

}

SIRE_END_HEADER

#endif
