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

#ifndef SireFF_FF_H
#define SireFF_FF_H

#include <QUuid>

#include <boost/shared_ptr.hpp>

#include "ffname.h"
#include "ffcomponent.h"

#include "SireBase/incremint.h"

#include "SireMol/moleculegroups.h"
#include "SireBase/properties.h"
#include "SireCAS/values.h"

SIRE_BEGIN_HEADER

namespace SireFF
{
class FF;
}

QDataStream& operator<<(QDataStream&, const SireFF::FF&);
QDataStream& operator>>(QDataStream&, SireFF::FF&);

namespace SireMol
{
class Molecule;
}

namespace SireFF
{

class FFMolGroup;

using SireBase::Properties;
using SireBase::Property;
using SireBase::PropertyMap;

using SireCAS::Symbol;
using SireCAS::Values;

using SireMol::MoleculeView;
using SireMol::ViewsOfMol;
using SireMol::Molecules;
using SireMol::MoleculeGroup;
using SireMol::Molecule;
using SireMol::MoleculeData;
using SireMol::MGID;
using SireMol::MGIdx;
using SireMol::MGNum;
using SireMol::MolNum;

class NullFF;

namespace detail
{
void throwForceFieldRestoreBug(const char *this_what, const char *ffield_what);
}

/** This class is the base class of all of the forcefield classes. 
    A forcefield is a collection of molecule groups (SireMol::MolGroups) 
    that had additional code to allow the energies, forces (and
    further derivatives) to be calculated for the molecules in those
    groups. The forcefield can, optionally, break this energy
    down into a set of components, and the energy, force (and/or
    further derivatives) of the components for the group(s) of
    molecules can be queried individually.

    FF derived objects are derived from MolGroupsBase, and hold
    molecule group objects that are derived from MoleculeGroup, e.g.
    
    FFGroupPvt : public MoleculeGroup  (used internally)
    
    FFGroup : public MoleculeGroup (used externally - needed so can
                               hold a copy of the FF) - FFGroupPvt
                               auto-converts to FFGroup when copied

    This class is derived into a series of types, e.g.
    
    G1FF - forcefield that only contains a single group of molecules
    G2FF - forcefield that contains two groups of molecules
    
    MultiFF - forcefield that is an amalgamation of other forcefields
    
    There are also additional pure virtual interfaces, e.g.
    
    FF2D - interface for all 2-dimensional forcefields
           (can be used to calculate forces in 2D space)
    
    FF3D - interface for all 3-dimensional forcefields
           (can be used to calculate forces in 3D space)
    
    InterCLJFF  == Inter2B3DG1FF<InterCLJPotential>
    IntraCLJFF  == Intra2B3DG1FF<IntraCLJPotential>

    @author Christopher Woods
*/
class SIREFF_EXPORT FF : public SireMol::MolGroupsBase
{

friend QDataStream& ::operator<<(QDataStream&, const FF&);
friend QDataStream& ::operator>>(QDataStream&, FF&);

friend class FFMolGroup; // so can call the group_???() functions
friend class FFComponent; // so can call setComponent and changeComponent

public:
    virtual ~FF();

    static const char* typeName()
    {
        return "SireFF::FF";
    }

    bool operator==(const FF &other) const;
    bool operator!=(const FF &other) const;

    /** Return the class name of the forcefield */
    virtual const char* what() const=0;

    /** Return a clone of this forcefield */
    virtual FF* clone() const=0;

    /** Return an object that describes all of the components
        of this forcefield (complete with SireCAS::Symbols to 
        uniquely ID each component) */
    virtual const FFComponent& components() const=0;

    virtual QString toString() const;

    const QUuid& UID() const;
    quint64 version() const;

    const FFName& name() const;
    void setName(const QString &name);

    void setName(const MGID &mgid, const QString &name);

    SireUnits::Dimension::MolarEnergy energy();
    SireUnits::Dimension::MolarEnergy energy(const Symbol &component);

    Values energies(const QSet<Symbol> &components);
    Values energies();

    /** Set the property with name 'name' to have the value 'value'
    
        \throw SireBase::missing_property
        \throw SireError::invalid_cast
        \throw SireError::incompatible_error
    */
    virtual bool setProperty(const QString &name, const Property &value)=0;
    
    /** Return the property that has the name 'name'
    
        \throw SireBase::missing_property
    */
    virtual const Property& property(const QString &name) const=0;

    /** Return whether or not this forcefield contains a property
        called 'name' */
    virtual bool containsProperty(const QString &name) const=0;

    /** Return the values of all of the properties of this forcefield */
    virtual const Properties& properties() const=0;

    virtual QStringList propertyKeys() const;

    /** Tell the forcefield that it has to recalculate everything from
        scratch */
    virtual void mustNowRecalculateFromScratch()=0;

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
    void update(const MoleculeView &molview, bool auto_commit=true);
    
    void update(const Molecules &molecules, bool auto_commit=true);
    void update(const MoleculeGroup &molgroup, bool auto_commit=true);

    void setContents(const MGID &mgid, const MoleculeView &molview);
    void setContents(const MGID &mgid, const ViewsOfMol &molview);
    void setContents(const MGID &mgid, const Molecules &molecules);
    void setContents(const MGID &mgid, const MoleculeGroup &molgroup);

    void setContents(const MGID &mgid, const MoleculeView &molview, 
                     const PropertyMap &map);
    void setContents(const MGID &mgid, const ViewsOfMol &molviews, 
                     const PropertyMap &map);
    void setContents(const MGID &mgid, const Molecules &molecules, 
                     const PropertyMap &map);
    void setContents(const MGID &mgid, const MoleculeGroup &molgroup, 
                     const PropertyMap &map);

    bool isDirty() const;
    bool isClean() const;

    static const NullFF& null();

protected:
    FF();
    FF(const QString &name);

    FF(const FF &other);

    FF& operator=(const FF &other);

    void setComponent(const Symbol &component, double nrg);
    void changeComponent(const Symbol &component, double delta);

    void incrementVersion();

    void setDirty();
    void setClean();

    const Values& currentEnergies() const;

    /** Virtual function used to trigger a recalculation of the total energy
        and of all of the component energies */
    virtual void recalculateEnergy()=0;

    /** Virtual functions used to add, remove or modify molecules in 
        this forcefield */
    virtual void group_setName(quint32 i, const QString &new_name)=0;
        
    virtual void group_add(quint32 i, const MoleculeView &molview,
                           const PropertyMap &map)=0;
    virtual void group_add(quint32 i, const ViewsOfMol &molviews, 
                           const PropertyMap &map)=0;
    virtual void group_add(quint32 i, const Molecules &molecules, 
                           const PropertyMap &map)=0;
    virtual void group_add(quint32 i, const MoleculeGroup &molgroup, 
                           const PropertyMap &map)=0;
    
    virtual bool group_addIfUnique(quint32 i, const MoleculeView &molview, 
                                   const PropertyMap &map)=0;
    virtual ViewsOfMol group_addIfUnique(quint32 i, const ViewsOfMol &molviews, 
                                         const PropertyMap &map)=0;
    virtual QList<ViewsOfMol> group_addIfUnique(quint32 i, const Molecules &molecules, 
                                                const PropertyMap &map)=0;
    virtual QList<ViewsOfMol> group_addIfUnique(quint32 i, const MoleculeGroup &molgroup, 
                                                const PropertyMap &map)=0;

    virtual bool group_remove(quint32 i, const MoleculeView &molview)=0;
    virtual ViewsOfMol group_remove(quint32 i, const ViewsOfMol &molviews)=0;
    virtual QList<ViewsOfMol> group_remove(quint32 i, const Molecules &molecules)=0;
    virtual QList<ViewsOfMol> group_remove(quint32 i, const MoleculeGroup &molgroup)=0;
    
    virtual bool group_removeAll(quint32 i, const MoleculeView &molview)=0;
    virtual ViewsOfMol group_removeAll(quint32 i, const ViewsOfMol &molviews)=0;
    virtual QList<ViewsOfMol> group_removeAll(quint32 i, const Molecules &molecules)=0;
    virtual QList<ViewsOfMol> group_removeAll(quint32 i, const MoleculeGroup &molgroup)=0;

    virtual ViewsOfMol group_remove(quint32 i, MolNum molnum)=0;
    virtual QList<ViewsOfMol> group_remove(quint32 i, const QSet<MolNum> &molnums)=0;

    virtual void group_removeAll(quint32 i)=0;

    virtual bool group_update(quint32 i, const MoleculeData &moldata, bool auto_commit)=0;

    virtual QList<Molecule> group_update(quint32 i, const Molecules &molecules,
                                         bool auto_commit)=0;
    virtual QList<Molecule> group_update(quint32 i, const MoleculeGroup &molgroup,
                                         bool auto_commit)=0;
    
    virtual bool group_setContents(quint32 i, const MoleculeView &molview, 
                                   const PropertyMap &map)=0;
    virtual bool group_setContents(quint32 i, const ViewsOfMol &molviews, 
                                   const PropertyMap &map)=0;
    virtual bool group_setContents(quint32 i, const Molecules &molecules, 
                                   const PropertyMap &map)=0;
    virtual bool group_setContents(quint32 i, const MoleculeGroup &molgroup, 
                                   const PropertyMap &map)=0;

    virtual void _pvt_updateName()=0;

private:

    /** The unique ID for this forcefield - this uniquely identifies
        this forcefield */
    QUuid uid;
    
    /** The version number of this forcefield - this is incremented every
        time this forcefield is changed. The combination of the uid and 
        version number is guaranteed to be unique, so that if two
        forcefields have the same uid and version number, then they
        are guaranteed to be identical (although there will be rollover
        in the version number after about 2^63 moves on 64 bit 
        machines, or 2^31 moves on 32 bit machines) */
    quint64 versn;

    /** Pointer to the incremint used to increment the version for
        this forcefield */
    boost::shared_ptr<SireBase::Incremint> version_ptr;

    /** The name of this forcefield */
    FFName ffname;

    /** All of the cached energy components in this forcefield, indexed
        by their symbol ID number (includes the total energy) */
    Values nrg_components;

    /** Whether or not this forcefield is dirty (requires an update) */
    bool isdirty;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Set the energy value of the component 'comp' */
inline void FF::setComponent(const Symbol &component, double nrg)
{
    nrg_components.set(component,nrg);
}

/** Change the existing value of the component 'comp' by delta */
inline void FF::changeComponent(const Symbol &component, double delta)
{
    nrg_components.set( component, delta + nrg_components.value(component) );
}

/** Return whether or not the forcefield is dirty (the energy
    needs to be recalculated) */
inline bool FF::isDirty() const
{
    return isdirty;
}

/** Return whether or not the forcefield is clean (the energy
    does not need to be recalculated) */
inline bool FF::isClean() const
{
    return not isDirty();
}

/** Record that this forcefield is now clean (has calculated the 
    energy of the current state) */
inline void FF::setClean()
{
    isdirty = false;
}

/** Record that this forcefield is dirty (the energy of the 
    current state is unknown) */
inline void FF::setDirty()
{
    isdirty = true;
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

SIRE_EXPOSE_CLASS( SireFF::FF )

SIRE_END_HEADER

#endif
