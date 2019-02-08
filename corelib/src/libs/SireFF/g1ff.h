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

#ifndef SIREFF_G1FF_H
#define SIREFF_G1FF_H

#include "ff.h"
#include "ffmolgroup.h"

SIRE_BEGIN_HEADER

namespace SireFF
{
class G1FF;
}

SIREFF_EXPORT QDataStream& operator<<(QDataStream&, const SireFF::G1FF&);
SIREFF_EXPORT QDataStream& operator>>(QDataStream&, SireFF::G1FF&);

namespace SireMol
{
class PartialMolecule;
}

namespace SireFF
{

using SireMol::PartialMolecule;

/** This is the base class of all forcefields that hold just
    a single group of molecules
    
    @author Christopher Woods
*/
class SIREFF_EXPORT G1FF : public FF
{

friend SIREFF_EXPORT QDataStream& ::operator<<(QDataStream&, const G1FF&);
friend SIREFF_EXPORT QDataStream& ::operator>>(QDataStream&, G1FF&);

public:
    ~G1FF();

    const MoleculeGroup& at(MGNum mgnum) const;

    void assertContains(MGNum mgnum) const;

    using FF::add;
    using FF::addIfUnique;
    using FF::removeAll;
    using FF::remove;
    using FF::setContents;

    void add(const MoleculeView &molview, const PropertyMap &map);
    void add(const ViewsOfMol &molviews, const PropertyMap &map);
    void add(const Molecules &molecules, const PropertyMap &map);
    void add(const MoleculeGroup &molgroup, const PropertyMap &map);
    
    void addIfUnique(const MoleculeView &molview, const PropertyMap &map);
    void addIfUnique(const ViewsOfMol &molviews, const PropertyMap &map);
    void addIfUnique(const Molecules &molecules, const PropertyMap &map);
    void addIfUnique(const MoleculeGroup &molgroup, const PropertyMap &map);

    void add(const MoleculeView &molview);
    void add(const ViewsOfMol &molviews);
    void add(const Molecules &molecules);
    void add(const MoleculeGroup &molgroup);
    
    void addIfUnique(const MoleculeView &molview);
    void addIfUnique(const ViewsOfMol &molviews);
    void addIfUnique(const Molecules &molecules);
    void addIfUnique(const MoleculeGroup &molgroup);

    bool removeAll();
    
    bool remove(const MoleculeView &molview);
    bool remove(const ViewsOfMol &molviews);
    bool remove(const Molecules &molecules);
    bool remove(const MoleculeGroup &molgroup);
    
    bool removeAll(const MoleculeView &molview);
    bool removeAll(const ViewsOfMol &molviews);
    bool removeAll(const Molecules &molecules);
    bool removeAll(const MoleculeGroup &molgroup);

    bool remove(MolNum molnum);
    bool remove(const QSet<MolNum> &molnums);

    void setContents(const MoleculeView &molview);
    void setContents(const ViewsOfMol &molview);
    void setContents(const Molecules &molecules);
    void setContents(const MoleculeGroup &molgroup);

    void setContents(const MoleculeView &molview, const PropertyMap &map);
    void setContents(const ViewsOfMol &molviews, const PropertyMap &map);
    void setContents(const Molecules &molecules, const PropertyMap &map);
    void setContents(const MoleculeGroup &molgroup, const PropertyMap &map);

    bool needsAccepting() const;
    void accept();

protected:
    G1FF(bool allow_overlap_of_atoms=false);
    G1FF(const G1FF &other);
    
    G1FF& operator=(const G1FF &other);

    ////
    //// Virtual functions that must be implemented
    //// by derived forcefields
    ////
    virtual void _pvt_added(const PartialMolecule &mol,
                            const PropertyMap &map)=0;

    virtual void _pvt_removed(const PartialMolecule &mol)=0;

    virtual void _pvt_changed(const Molecule &molecule, bool auto_commit)=0;
    virtual void _pvt_changed(const QList<Molecule> &molecules, bool auto_commit)=0;
    
    virtual void _pvt_removedAll()=0;
        
    virtual bool _pvt_wouldChangeProperties(MolNum molnum, 
                                            const PropertyMap &map) const=0;
        
    ////
    //// Virtual functions that must be changed if this 
    //// forcefield allows overlapping atoms (most don't!)
    ////
    virtual void _pvt_added(const ViewsOfMol &mol,
                            const PropertyMap &map);
                            
    virtual void _pvt_removed(const ViewsOfMol &mol);

    virtual void _pvt_removedAll(const PartialMolecule &mol);
    virtual void _pvt_removedAll(const ViewsOfMol &mol);

    ////
    //// Implementation of FF virtual functions
    ////

    const MoleculeGroup& getGroup(MGNum mgnum) const;
    
    void getGroups(const QList<MGNum> &mgnums,
                   QVarLengthArray<const MoleculeGroup*,10> &groups) const;

    QHash<MGNum,const MoleculeGroup*> getGroups() const;
    
    void group_setName(quint32 i, const QString &new_name);
        
    void group_add(quint32 i, const MoleculeView &molview,
                           const PropertyMap &map);
    void group_add(quint32 i, const ViewsOfMol &molviews, 
                           const PropertyMap &map);
    void group_add(quint32 i, const Molecules &molecules, 
                           const PropertyMap &map);
    void group_add(quint32 i, const MoleculeGroup &molgroup, 
                           const PropertyMap &map);
    
    bool group_addIfUnique(quint32 i, const MoleculeView &molview, 
                           const PropertyMap &map);
    ViewsOfMol group_addIfUnique(quint32 i, const ViewsOfMol &molviews, 
                                 const PropertyMap &map);
    QList<ViewsOfMol> group_addIfUnique(quint32 i, const Molecules &molecules, 
                                        const PropertyMap &map);
    QList<ViewsOfMol> group_addIfUnique(quint32 i, const MoleculeGroup &molgroup, 
                                        const PropertyMap &map);

    bool group_remove(quint32 i, const MoleculeView &molview);
    ViewsOfMol group_remove(quint32 i, const ViewsOfMol &molviews);
    QList<ViewsOfMol> group_remove(quint32 i, const Molecules &molecules);
    QList<ViewsOfMol> group_remove(quint32 i, const MoleculeGroup &molgroup);
    
    bool group_removeAll(quint32 i, const MoleculeView &molview);
    ViewsOfMol group_removeAll(quint32 i, const ViewsOfMol &molviews);
    QList<ViewsOfMol> group_removeAll(quint32 i, const Molecules &molecules);
    QList<ViewsOfMol> group_removeAll(quint32 i, const MoleculeGroup &molgroup);

    ViewsOfMol group_remove(quint32 i, MolNum molnum);
    QList<ViewsOfMol> group_remove(quint32 i, const QSet<MolNum> &molnums);

    void group_removeAll(quint32 i);

    bool group_update(quint32 i, const MoleculeData &moldata, bool auto_commit);

    QList<Molecule> group_update(quint32 i, const Molecules &molecules, bool auto_commit);
    QList<Molecule> group_update(quint32 i, const MoleculeGroup &molgroup, bool auto_commit);
    
    bool group_setContents(quint32 i, const MoleculeView &molview, 
                           const PropertyMap &map);
    bool group_setContents(quint32 i, const ViewsOfMol &molviews, 
                           const PropertyMap &map);
    bool group_setContents(quint32 i, const Molecules &molecules, 
                           const PropertyMap &map);
    bool group_setContents(quint32 i, const MoleculeGroup &molgroup, 
                           const PropertyMap &map);

    void _pvt_updateName();

    // implementation of MolGroupsBase virtual functions
    void reindex();

private:
    void assertValidGroup(quint32 i) const;

    void assertNoOverlap(const MoleculeView &molview) const;
    void assertNoOverlap(const MoleculeGroup &molgroup,
                         const MoleculeView &molview) const;

    /** The single group of molecules in this forcefield */
    detail::FFMolGroupPvt molgroup;
    
    /** Whether or not this forcefield allows overlap of atoms
        (e.g. allowing an atom to appear several times in the forcefield) */
    bool allow_overlap_of_atoms;
};

}

SIRE_EXPOSE_CLASS( SireFF::G1FF )

SIRE_END_HEADER

#endif
