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

#ifndef SIREFF_G2FF_H
#define SIREFF_G2FF_H

#include "ff.h"
#include "ffmolgroup.h"

SIRE_BEGIN_HEADER

namespace SireFF
{
class G2FF;
}

QDataStream& operator<<(QDataStream&, const SireFF::G2FF&);
QDataStream& operator>>(QDataStream&, SireFF::G2FF&);

namespace SireMol
{
class PartialMolecule;
}

namespace SireFF
{

using SireMol::PartialMolecule;

/** This is the base class of all forcefields that hold
    two groups of molecules, e.g. InterGroupCLJFF, that calculates
    the CLJ energy between two groups of molecules, or
    QMMMFF, that calculates the QM/MM energy of a QM group
    of molecules interacting with an MM group
    
    @author Christopher Woods
*/
class SIREFF_EXPORT G2FF : public FF
{

friend QDataStream& ::operator<<(QDataStream&, const G2FF&);
friend QDataStream& ::operator>>(QDataStream&, G2FF&);

public:
    ~G2FF();
    
    const MoleculeGroup& at(MGNum mgnum) const;
    
    void assertContains(MGNum mgnum) const;
    
    void accept();
    bool needsAccepting() const;
    
protected:
    G2FF(bool allow_overlap_of_atoms=false);
    G2FF(const G2FF &other);

    G2FF& operator=(const G2FF &other);

    ////
    //// Virtual functions that must be implemented
    //// by derived forcefields
    ////
    virtual void _pvt_added(quint32 groupid,
                            const PartialMolecule &mol,
                            const PropertyMap &map)=0;

    virtual void _pvt_removed(quint32 groupid,
                              const PartialMolecule &mol)=0;

    virtual void _pvt_changed(quint32 groupid,
                              const Molecule &molecule,
                              bool auto_commit)=0;
                              
    virtual void _pvt_changed(quint32 groupid, 
                              const QList<Molecule> &molecules,
                              bool auto_commit)=0;
    
    virtual void _pvt_removedAll(quint32 groupid)=0;
        
    virtual bool _pvt_wouldChangeProperties(quint32 groupid,
                                            MolNum molnum, 
                                            const PropertyMap &map) const=0;

    ////
    //// Virtual functions that must be changed if this 
    //// forcefield allows overlapping atoms (most don't!)
    ////
    virtual void _pvt_added(quint32 groupid,
                            const ViewsOfMol &mol,
                            const PropertyMap &map);
                            
    virtual void _pvt_removed(quint32 groupid,
                              const ViewsOfMol &mol);

    virtual void _pvt_removedAll(quint32 groupid,
                                 const PartialMolecule &mol);

    virtual void _pvt_removedAll(quint32 groupid,
                                 const ViewsOfMol &mol);

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

    /// implementation of virtual MolGroupsBase functions
    void reindex();

private:
    void assertValidGroup(quint32 i) const;

    void assertNoOverlap(const MoleculeView &molview) const;
    void assertNoOverlap(const MoleculeGroup &molgroup,
                         const MoleculeView &molview) const;

    void assertNoOverlap(const MoleculeGroup &molgroup,
                         const Molecules &molecules) const;

    /** The two groups of molecules in this forcefield */
    detail::FFMolGroupPvt molgroup[2];
    
    /** Whether or not this forcefield allows overlap of atoms
        (e.g. allowing an atom to appear several times in the forcefield) */
    bool allow_overlap_of_atoms;
};

}

SIRE_EXPOSE_CLASS( SireFF::G2FF )

SIRE_END_HEADER

#endif
