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

#ifndef SIREFF_FFMOLGROUP_H
#define SIREFF_FFMOLGROUP_H

#include "SireMol/moleculegroup.h"
#include "SireMol/mgidx.h"

#include "SireFF/ff.h"
#include "SireFF/forcefield.h"

SIRE_BEGIN_HEADER

namespace SireFF
{
class FFMolGroup;

namespace detail
{
class FFMolGroupPvt;
}

}

SIREFF_EXPORT QDataStream& operator<<(QDataStream&, const SireFF::FFMolGroup&);
SIREFF_EXPORT QDataStream& operator>>(QDataStream&, SireFF::FFMolGroup&);

SIREFF_EXPORT QDataStream& operator<<(QDataStream&, const SireFF::detail::FFMolGroupPvt&);
SIREFF_EXPORT QDataStream& operator>>(QDataStream&, SireFF::detail::FFMolGroupPvt&);

namespace SireFF
{

namespace detail{ class FFMolGroupPvt; }

using SireBase::PropertyMap;

using SireMol::MGIdx;
using SireMol::MolNum;
using SireMol::MoleculeView;
using SireMol::ViewsOfMol;
using SireMol::Molecules;
using SireMol::MoleculeGroup;
using SireMol::MoleculeData;
using SireMol::Molecule;

/** This class holds a molecule group that is part of a forcefield.
    This actually holds a copy of the forcefield that contains
    this molecule group, so that changes to this group change
    the actual forcefield to which this group belongs. This is
    the publicly available class that corresponds to the private
    SireFF::detail::FFMolGroupPvt class.
    
    @author Christopher Woods
*/
class SIREFF_EXPORT FFMolGroup 
                : public SireBase::ConcreteProperty<FFMolGroup,MoleculeGroup>
{

friend SIREFF_EXPORT QDataStream& ::operator<<(QDataStream&, const FFMolGroup&);
friend SIREFF_EXPORT QDataStream& ::operator>>(QDataStream&, FFMolGroup&);

public:
    FFMolGroup();
    
    FFMolGroup(const detail::FFMolGroupPvt &ffmolgroup);

    FFMolGroup(const MoleculeGroup &other);
    
    FFMolGroup(const FFMolGroup &other);
    
    virtual ~FFMolGroup();
    
    static const char* typeName();
    
    const char* what() const
    {
        return FFMolGroup::typeName();
    }
    
    FFMolGroup* clone() const;
    
    FFMolGroup& operator=(const FFMolGroup &other);
    FFMolGroup& operator=(const MoleculeGroup &other);

    const FF& forceField() const;
    
    MGIdx index() const;

    void setName(const QString &name);

    void add(const MoleculeView &molview);
    void add(const ViewsOfMol &molviews);
    void add(const Molecules &molecules);
    void add(const MoleculeGroup &molgroup);
    
    void add(const MoleculeView &molview, const PropertyMap &map);
    void add(const ViewsOfMol &molviews, const PropertyMap &map);
    void add(const Molecules &molecules, const PropertyMap &map);
    void add(const MoleculeGroup &molgroup, const PropertyMap &map);
    
    bool addIfUnique(const MoleculeView &molview);
    ViewsOfMol addIfUnique(const ViewsOfMol &molviews);
    QList<ViewsOfMol> addIfUnique(const Molecules &molecules);
    QList<ViewsOfMol> addIfUnique(const MoleculeGroup &molgroup);

    bool addIfUnique(const MoleculeView &molview, const PropertyMap &map);
    ViewsOfMol addIfUnique(const ViewsOfMol &molviews, const PropertyMap &map);
    QList<ViewsOfMol> addIfUnique(const Molecules &molecules, const PropertyMap &map);
    QList<ViewsOfMol> addIfUnique(const MoleculeGroup &molgroup, const PropertyMap &map);

    bool remove(const MoleculeView &molview);
    ViewsOfMol remove(const ViewsOfMol &molviews);
    QList<ViewsOfMol> remove(const Molecules &molecules);
    QList<ViewsOfMol> remove(const MoleculeGroup &molgroup);
    
    bool removeAll(const MoleculeView &molview);
    ViewsOfMol removeAll(const ViewsOfMol &molviews);
    QList<ViewsOfMol> removeAll(const Molecules &molecules);
    QList<ViewsOfMol> removeAll(const MoleculeGroup &molgroup);

    ViewsOfMol remove(MolNum molnum);
    QList<ViewsOfMol> remove(const QSet<MolNum> &molnums);

    void removeAll();

    bool update(const MoleculeData &moldata, bool auto_commit=true);

    QList<Molecule> update(const Molecules &molecules, bool auto_commit=true);
    QList<Molecule> update(const MoleculeGroup &molgroup, bool auto_commit=true);
    
    bool setContents(const MoleculeView &molview);
    bool setContents(const ViewsOfMol &molviews);
    bool setContents(const Molecules &molecules);
    bool setContents(const MoleculeGroup &molgroup);
    
    bool setContents(const MoleculeView &molview, const PropertyMap &map);
    bool setContents(const ViewsOfMol &molviews, const PropertyMap &map);
    bool setContents(const Molecules &molecules, const PropertyMap &map);
    bool setContents(const MoleculeGroup &molgroup, const PropertyMap &map);
    
    bool needsAccepting() const;
    void accept();
    
private:
    void assertNotNull() const;
    void updateGroup();

    /** The index of this group in the forcefield */
    MGIdx mgidx;

    /** A copy of the forcefield that contains this molecule group */
    ForceField ffield;
};

namespace detail
{

/** Private class used by the G1FF, G2FF etc. to hold the molecule
    group in a forcefield. This is a private, internal class that
    should not be used in other code. The 'clone()' function of this
    class returns an FFMolGroup, as it is not possible to copy
    this class (as it is internal to the forcefield)
    
    @author Christopher Woods
*/
class SIREFF_EXPORT FFMolGroupPvt : public MoleculeGroup
{

friend SIREFF_EXPORT QDataStream& ::operator<<(QDataStream&, const FFMolGroupPvt&);
friend SIREFF_EXPORT QDataStream& ::operator>>(QDataStream&, FFMolGroupPvt&);

public:
    FFMolGroupPvt();
    
    FFMolGroupPvt(const QString &name, quint32 i, FF *ff);
    
    FFMolGroupPvt(const FFMolGroupPvt &other);
    
    virtual ~FFMolGroupPvt();

    static const char* typeName()
    {
        return "SireFF::detail::FFMolGroupPvt";
    }
    
    const char* what() const
    {
        return FFMolGroupPvt::typeName();
    }

    FFMolGroupPvt& operator=(const FFMolGroupPvt &other);

    MoleculeGroup* clone() const;

    const FF& forceField() const
    {
        assertNotNull();
        return *ffield;
    }
    
    FF& forceField()
    {
        assertNotNull();
        return *ffield;
    }

    MGIdx index() const
    {
        assertNotNull();
        return mgidx;
    }

    void setParent(FF *new_parent);
    void setIndex(quint32 new_idx);

private:
    void assertNotNull() const;

    /** The index of this group in the forcefield - this uniquely
        identifies this group in the forcefield */
    MGIdx mgidx;

    /** Pointer to the forcefield that contains this group
          - this allows the forcefield that contains this
            group to be copied out of this group */
    FF *ffield;
};

} // end of namespace detail

} // end of namespace SireFF

Q_DECLARE_METATYPE(SireFF::FFMolGroup);

SIRE_EXPOSE_CLASS( SireFF::FFMolGroup )

SIRE_END_HEADER

#endif
