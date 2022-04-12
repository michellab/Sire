/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2007  Christopher Woods
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

#ifndef SIREMOL_MOLECULEGROUPS_H
#define SIREMOL_MOLECULEGROUPS_H

#include <QVarLengthArray>
#include <QHash>
#include <QList>

#include <limits>

#include "SireBase/property.h"

#include "moleculegroup.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class MolGroupsBase;
class MoleculeGroups;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::MolGroupsBase&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::MolGroupsBase&);

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::MoleculeGroups&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::MoleculeGroups&);

namespace SireBase
{
class Slice;
}

namespace SireMol
{

using SireBase::Property;
using SireBase::ConcreteProperty;

template<class T>
class Selector;

class Molecules;
class ViewsOfMol;
class MoleculeView;
class MoleculeData;
class Segment;
class Chain;
class Residue;
class CutGroup;
class Atom;

class MGNum;
class MGName;
class MGIdx;
class MGID;

class MolNum;
class MolName;
class MolIdx;
class MolID;

class SegID;
class ChainID;
class ResID;
class CGID;
class AtomID;

/** This is the base class of all MoleculeGroups objects.
    These are containers for MoleculeGroup objects, thereby
    allowing lots of MoleculeGroup objects to be collected
    together and indexed (e.g. so that you can find all
    CA atoms in the "proteins" group). This is the virtual
    base class of the hierarchy - see MolGroups for a simple
    concrete instantiation.

    @author Christopher Woods
*/
class SIREMOL_EXPORT MolGroupsBase : public Property
{

friend SIREMOL_EXPORT QDataStream& ::operator<<(QDataStream&, const MolGroupsBase&);
friend SIREMOL_EXPORT QDataStream& ::operator>>(QDataStream&, MolGroupsBase&);

public:
    virtual ~MolGroupsBase();

    const MoleculeGroup& operator[](MGNum mgnum) const;
    const MoleculeGroup& operator[](const MGName &mgname) const;
    const MoleculeGroup& operator[](MGIdx mgidx) const;
    const MoleculeGroup& operator[](const MGID &mgid) const;

    ViewsOfMol operator[](int i) const;
    ViewsOfMol operator[](const QString &name) const;

    ViewsOfMol operator[](MolNum molnum) const;
    ViewsOfMol operator[](const MolID &molid) const;

    Segment operator[](const SegID &segid) const;
    Chain operator[](const ChainID &chainid) const;
    Residue operator[](const ResID &resid) const;
    CutGroup operator[](const CGID &cgid) const;
    Atom operator[](const AtomID &atomid) const;
    QList<MolViewPtr> operator[](const SireBase::Slice &slice) const;

    virtual MolGroupsBase* clone() const=0;

    static const char* typeName()
    {
        return "SireMol::MolGroupsBase";
    }

    MGNum getGroupNumber(MGNum mgnum) const;
    MGNum getGroupNumber(MGIdx mgidx) const;
    MGNum getGroupNumber(const MGName &mgname) const;
    MGNum getGroupNumber(const MGID &mgid) const;

    MGIdx mgIdx(MGNum mgnum) const;

    QList<MGNum> map(MGNum mgnum) const;
    QList<MGNum> map(MGIdx mgidx) const;
    QList<MGNum> map(const MGName &mgname) const;
    QList<MGNum> map(const MGID &mgid) const;

    MolNum getMoleculeNumber(MolNum molnum) const;
    MolNum getMoleculeNumber(MolIdx molidx) const;
    MolNum getMoleculeNumber(const MolName &molname) const;
    MolNum getMoleculeNumber(const MolID &molid) const;

    quint64 getMoleculeVersion(MolNum molnum) const;
    quint64 getMoleculeVersion(const MolID &molid) const;

    QList<MolNum> molNums() const;
    QList<MolNum> getMoleculeNumbers() const;

    QList<MolNum> map(MolNum molnum) const;
    QList<MolNum> map(MolIdx molidx) const;
    QList<MolNum> map(const MolName &molname) const;
    QList<MolNum> map(const MolID &molid) const;

    virtual const MoleculeGroup& at(MGNum mgnum) const=0;

    const MoleculeGroup& at(MGIdx mgidx) const;
    const MoleculeGroup& at(const MGName &mgname) const;
    const MoleculeGroup& at(const MGID &mgid) const;

    ViewsOfMol at(MolNum molnum) const;
    ViewsOfMol at(const MolID &molid) const;

    Segment at(const SegID &segid) const;
    Chain at(const ChainID &chainid) const;
    Residue at(const ResID &resid) const;
    CutGroup at(const CGID &cgid) const;
    Atom at(const AtomID &atomid) const;

    SelectResult search(const QString &search_string) const;

    const MoleculeGroup& select(const MGID &mgid) const;
    ViewsOfMol select(const MolID &molid) const;
    Segment select(const SegID &segid) const;
    Chain select(const ChainID &chainid) const;
    Residue select(const ResID &resid) const;
    CutGroup select(const CGID &cgid) const;
    Atom select(const AtomID &atomid) const;

    QList<MolGroupPtr> selectAll() const;

    QList<MolGroupPtr> selectAll(MGNum mgnum) const;
    QList<MolGroupPtr> selectAll(MGIdx mgidx) const;
    QList<MolGroupPtr> selectAll(const MGName &mgname) const;
    QList<MolGroupPtr> selectAll(const MGID &mgid) const;

    QList<ViewsOfMol> selectAll(const MolID &molid) const;

    QHash< MolNum,Selector<Segment> > selectAll(const SegID &segid) const;
    QHash< MolNum,Selector<Chain> > selectAll(const ChainID &chainid) const;
    QHash< MolNum,Selector<Residue> > selectAll(const ResID &resid) const;
    QHash< MolNum,Selector<CutGroup> > selectAll(const CGID &cgid) const;
    QHash< MolNum,Selector<Atom> > selectAll(const AtomID &atomid) const;

    const MoleculeGroup& group(MGNum mgnum) const;
    const MoleculeGroup& group(const MGName &mgname) const;
    const MoleculeGroup& group(MGIdx mgidx) const;
    const MoleculeGroup& group(const MGID &mgid) const;

    QList<MGNum> groupNumbers() const;
    QList<MGName> groupNames() const;

    QList<MolGroupPtr> groups() const;

    QList<MolGroupPtr> groups(MGNum mgnum) const;
    QList<MolGroupPtr> groups(MGIdx mgidx) const;
    QList<MolGroupPtr> groups(const MGName &mgname) const;
    QList<MolGroupPtr> groups(const MGID &mgid) const;

    ViewsOfMol molecule(MolNum molnum) const;
    ViewsOfMol molecule(const MolID &molid) const;

    QList<ViewsOfMol> molecules(MolNum molnum) const;
    QList<ViewsOfMol> molecules(const MolID &molid) const;

    Segment segment(const SegID &segid) const;
    Chain chain(const ChainID &chainid) const;
    Residue residue(const ResID &resid) const;
    CutGroup cutGroup(const CGID &cgid) const;
    Atom atom(const AtomID &atomid) const;

    QHash< MolNum,Selector<Segment> > segments(const SegID &segid) const;
    QHash< MolNum,Selector<Chain> > chains(const ChainID &chainid) const;
    QHash< MolNum,Selector<Residue> > residues(const ResID &resid) const;
    QHash< MolNum,Selector<CutGroup> > cutGroups(const CGID &cgid) const;
    QHash< MolNum,Selector<Atom> > atoms(const AtomID &atomid) const;

    bool contains(MGNum mgnum) const;

    bool contains(MolNum molnum) const;
    bool contains(const QList<MolNum> &molnums) const;

    bool contains(const MoleculeView &molview) const;
    bool contains(const ViewsOfMol &molviews) const;
    bool contains(const Molecules &molecules) const;

    bool intersects(const MoleculeView &molview) const;
    bool intersects(const Molecules &other) const;

    const QList<MGNum>& groupsContaining(MolNum molnum) const;

    int nGroups() const;
    int count() const;

    int nMolecules() const;

    int nViews() const;
    int nViews(MolNum molnum) const;

    bool isEmpty() const;

    virtual void accept()=0;
    virtual bool needsAccepting() const=0;

    Molecules molecules() const;
    Molecules molecules(const MGID &mgid) const;

    QList<MGNum> mgNums() const;
    QList<MGName> mgNames() const;

    void assertContains(MolNum molnum) const;
    void assertContains(const MolID &molid) const;

    void assertContains(MGNum mgnum) const;
    void assertContains(const MGID &mgid) const;

    virtual void add(const MoleculeView &molview, const MGID &mgid)=0;
    virtual void add(const ViewsOfMol &molviews, const MGID &mgid)=0;
    virtual void add(const Molecules &molecules, const MGID &mgid)=0;
    virtual void add(const MoleculeGroup &molgroup, const MGID &mgid)=0;

    virtual void addIfUnique(const MoleculeView &molview,
                             const MGID &mgid)=0;
    virtual void addIfUnique(const ViewsOfMol &molviews,
                             const MGID &mgid)=0;
    virtual void addIfUnique(const Molecules &molecules,
                             const MGID &mgid)=0;
    virtual void addIfUnique(const MoleculeGroup &molgroup,
                             const MGID &mgid)=0;

    void unite(const MoleculeView &molview, const MGID &mgid);
    void unite(const ViewsOfMol &molviews, const MGID &mgid);
    void unite(const Molecules &molecules, const MGID &mgid);
    void unite(const MoleculeGroup &molgroup, const MGID &mgid);

    virtual bool remove(const MoleculeGroup &molgroup);

    virtual bool remove(const MoleculeView &molview);
    virtual bool remove(const ViewsOfMol &molviews);
    virtual bool remove(const Molecules &molecules);

    virtual bool removeAll(const MoleculeView &molview);
    virtual bool removeAll(const ViewsOfMol &molviews);
    virtual bool removeAll(const Molecules &molecules);
    virtual bool removeAll(const MoleculeGroup &molgroup);

    virtual bool remove(MolNum molnum);
    virtual bool remove(const QSet<MolNum> &molnums);

    virtual bool remove(const MolID &molid);
    virtual bool remove(const MGID &mgid);

    virtual bool removeAll(const MGID &mgid)=0;

    virtual bool removeAll();

    virtual bool remove(const MoleculeView &molview, const MGID &mgid)=0;
    virtual bool remove(const ViewsOfMol &molviews, const MGID &mgid)=0;
    virtual bool remove(const Molecules &molecules, const MGID &mgid)=0;
    virtual bool remove(const MoleculeGroup &molgroup, const MGID &mgid)=0;

    virtual bool removeAll(const MoleculeView &molview, const MGID &mgid)=0;
    virtual bool removeAll(const ViewsOfMol &molviews, const MGID &mgid)=0;
    virtual bool removeAll(const Molecules &molecules, const MGID &mgid)=0;
    virtual bool removeAll(const MoleculeGroup &molgroup, const MGID &mgid)=0;

    virtual bool remove(MolNum molnum, const MGID &mgid)=0;
    virtual bool remove(const QSet<MolNum> &molnums, const MGID &mgid)=0;

    virtual void update(const MoleculeData &moldata, bool auto_commit=true)=0;
    void update(const MoleculeView &molview, bool auto_commit=true);

    virtual void update(const Molecules &molecules, bool auto_commit=true)=0;
    virtual void update(const MoleculeGroup &molgroup, bool auto_commit=true)=0;

    virtual void setContents(const MGID &mgid, const MoleculeView &molview)=0;
    virtual void setContents(const MGID &mgid, const ViewsOfMol &molviews)=0;
    virtual void setContents(const MGID &mgid, const Molecules &molecules)=0;
    virtual void setContents(const MGID &mgid, const MoleculeGroup &molgroup)=0;

    static const MoleculeGroups& null();

protected:
    MolGroupsBase();

    MolGroupsBase(const MolGroupsBase &other);

    MolGroupsBase& operator=(const MolGroupsBase &other);

    virtual const MoleculeGroup& getGroup(MGNum mgnum) const=0;

    virtual void getGroups(const QList<MGNum> &mgnums,
                           QVarLengthArray<const MoleculeGroup*,10> &groups) const=0;

    virtual QHash<MGNum,const MoleculeGroup*> getGroups() const=0;

    bool needToUpdate(const MoleculeData &moldata) const;

    const MoleculeData& matchToExistingVersion(const MoleculeData &moldata) const;

    Molecules matchToExistingVersion(const Molecules &molecules) const;

    virtual void reindex()=0;

    void addToIndex(const MoleculeGroup &molgroup);
    void addToIndex(MGNum mgnum, MolNum molnum);
    void addToIndex(MGNum mgnum, const QSet<MolNum> &molnums);
    void addToIndex(MGNum mgnum, const QList<MolNum> &molnums);

    void removeFromIndex(MGNum mgnum);
    void removeFromIndex(MolNum molnum);

    void removeFromIndex(MGNum mgnum, MolNum molnum);
    void removeFromIndex(MGNum mgnum, const QSet<MolNum> &molnums);

    void changeNameIndex(MGNum mgnum, const MGName &old_name,
                         const MGName &new_name);

    void clearIndex(MGNum mgnum);
    void clearIndex();

private:
    /** This index keeps an order of MoleculeGroup objects */
    QList<MGNum> mgidx_to_num;

    /** This index maps the names of the MoleculeGroup objects */
    QHash< MGName, QList<MGNum> > mgname_to_mgnum;

    /** This is an index of which groups contain which molecules */
    QHash< MolNum, QList<MGNum> > molnum_to_mgnum;
};

/** This class holds a collection of MoleculeGroup objects. This
    allows multiple groups to be themselves grouped together.
    This is a virtual class, which can hold the virtual
    MoleculeGroup class objects. This can be used usefully in
    several situations, e.g.;

    System is derived from MolGroupsBase, and uses the MolGroups
    code to manage the indexing and version management of all
    of the molecules in the system.

    The forcefields are also derived from MolGroupsBase, allowing
    the MolGroups code to do the indexing and version management
    of molecules in a forcefield. Also, this allows easy
    management of multiple groups in a forcefield, e.g.
    QM molecules and MM molecules, or group A and group B.

    While System and the forcefields are derived from MolGroupsBase,
    this class, MolGroups, provides a concrete class that allows
    the user to easily group together different MoleculeGroups.

    @author Christopher Woods
*/
class SIREMOL_EXPORT MoleculeGroups
            : public ConcreteProperty<MoleculeGroups,MolGroupsBase>
{

friend SIREMOL_EXPORT QDataStream& ::operator<<(QDataStream&, const MoleculeGroups&);
friend SIREMOL_EXPORT QDataStream& ::operator>>(QDataStream&, MoleculeGroups&);

public:
    MoleculeGroups();

    MoleculeGroups(const MoleculeGroup &molgroup);

    MoleculeGroups(const QList<MolGroupPtr> &molgroups);

    MoleculeGroups(const MoleculeGroups &other);

    ~MoleculeGroups();

    MoleculeGroups* clone() const;

    static const char* typeName();

    MoleculeGroups& operator=(const MoleculeGroups &other);

    bool operator==(const MoleculeGroups &other) const;
    bool operator!=(const MoleculeGroups &other) const;

    MoleculeGroups& operator+=(const MoleculeGroup &molgroup);
    MoleculeGroups& operator-=(const MoleculeGroup &molgroup);

    MoleculeGroups& operator-=(const MGID &mgid);

    MoleculeGroups& operator-=(const Molecules &molecules);
    MoleculeGroups& operator-=(const MolID &molid);

    void add(const MoleculeGroup &molgroup);

    bool remove(const MGID &mgid);

    ///////////////////////////////////////////////
    /// Pure virtual functions of MoleculeGroupsBase ///
    ///////////////////////////////////////////////

    const MoleculeGroup& at(MGNum mgnum) const;

    void add(const MoleculeView &molview, const MGID &mgid);
    void add(const ViewsOfMol &molviews, const MGID &mgid);
    void add(const Molecules &molecules, const MGID &mgid);
    void add(const MoleculeGroup &molgroup, const MGID &mgid);

    void addIfUnique(const MoleculeView &molview, const MGID &mgid);
    void addIfUnique(const ViewsOfMol &molviews, const MGID &mgid);
    void addIfUnique(const Molecules &molecules, const MGID &mgid);
    void addIfUnique(const MoleculeGroup &molgroup, const MGID &mgid);

    bool remove(const MoleculeGroup &molgroup);

    bool remove(const MolID &molid);

    bool remove(const MoleculeView &molview, const MGID &mgid);
    bool remove(const ViewsOfMol &molviews, const MGID &mgid);
    bool remove(const Molecules &molecules, const MGID &mgid);
    bool remove(const MoleculeGroup &molgroup, const MGID &mgid);

    bool removeAll(const MoleculeView &molview, const MGID &mgid);
    bool removeAll(const ViewsOfMol &molviews, const MGID &mgid);
    bool removeAll(const Molecules &molecules, const MGID &mgid);
    bool removeAll(const MoleculeGroup &molgroup, const MGID &mgid);

    bool removeAll(const MGID &mgid);

    bool remove(MolNum molnum, const MGID &mgid);
    bool remove(const QSet<MolNum> &molnums, const MGID &mgid);

    void update(const MoleculeData &moldata, bool auto_commit=true);

    void update(const Molecules &molecules, bool auto_commit=true);
    void update(const MoleculeGroup &molgroup, bool auto_commit=true);

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

    bool _pvt_remove(MGNum mgnum);

private:
    /** All of the MoleculeGroup objects in this collection,
        indexed by molecule group number */
    QHash<MGNum,MolGroupPtr> mgroups;
};

typedef SireBase::PropPtr<MolGroupsBase> MolGroupsPtr;

}

Q_DECLARE_METATYPE(SireMol::MoleculeGroups);

SIRE_EXPOSE_CLASS( SireMol::MolGroupsBase )
SIRE_EXPOSE_CLASS( SireMol::MoleculeGroups )

SIRE_EXPOSE_PROPERTY( SireMol::MolGroupsPtr, SireMol::MolGroupsBase )

SIRE_END_HEADER

#endif
