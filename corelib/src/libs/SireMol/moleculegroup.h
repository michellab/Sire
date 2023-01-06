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

#ifndef SIREMOL_MOLECULEGROUP_H
#define SIREMOL_MOLECULEGROUP_H

#include <QList>

#include <boost/tuple/tuple.hpp>

#include "SireBase/property.h"
#include "SireBase/shareddatapointer.hpp"

#include "molecules.h"
#include "molgroupworkspace.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class MoleculeGroup;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::MoleculeGroup&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::MoleculeGroup&);

namespace SireID
{
class Index;
}

namespace SireBase
{
class Version;
}

namespace SireMol
{

namespace detail
{
class MolGroupPvt;
}

using SireID::Index;

class MoleculeData;
class Molecules;
class Molecule;
class MoleculeView;
class PartialMolecule;
class ViewsOfMol;
class SelectResult;

class MolNum;
class MolNumViewIdx;

class MGName;
class MGNum;

using SireBase::ConcreteProperty;
using SireBase::Property;
using SireBase::Version;

/** This is the virtual base class of all MoleculeGroup type
    objects. Molecule groups are groups of molecules that also
    provide full indexing, versioning and identification support.

    Molecule group form the foundation of forcefields (which use
    groups to hold the molecules), Systems (again, use groups
    to hold the molecules) and Moves (use groups to select which
    molecules should be moved).

    Molecule groups provide the common interface for indexing,
    searching and managing groups of molecules.

    @author Christopher Woods
*/
class SIREMOL_EXPORT MoleculeGroup : public ConcreteProperty<MoleculeGroup,Property>
{

friend SIREMOL_EXPORT QDataStream& ::operator<<(QDataStream&, const MoleculeGroup&);
friend SIREMOL_EXPORT QDataStream& ::operator>>(QDataStream&, MoleculeGroup&);

public:

    typedef Molecules::const_iterator const_iterator;
    typedef Molecules::iterator iterator;

    MoleculeGroup();

    MoleculeGroup(const Molecules &molecules);

    MoleculeGroup(const QString &name);

    MoleculeGroup(const QString &name, const MoleculeView &molview);
    MoleculeGroup(const QString &name, const Molecules &molecules);
    MoleculeGroup(const QString &name, const MoleculeGroup &other);

    MoleculeGroup(const SelectResult &result);
    MoleculeGroup(const QString &name, const SelectResult &result);

    MoleculeGroup(const MoleculeGroup &other);

    virtual ~MoleculeGroup();

    static const char* typeName();

    virtual const char* what() const
    {
        return MoleculeGroup::typeName();
    }

    virtual MoleculeGroup* clone() const;

    MoleculeGroup& operator=(const MoleculeGroup &other);

    virtual bool operator==(const MoleculeGroup &other) const;
    virtual bool operator!=(const MoleculeGroup &other) const;

    const ViewsOfMol& operator[](MolNum molnum) const;
    const ViewsOfMol& operator[](MolIdx molidx) const;
    const ViewsOfMol& operator[](const MolName &molname) const;
    const ViewsOfMol& operator[](const MolID &molid) const;

    PartialMolecule operator[](const boost::tuple<MolNum,Index> &viewidx) const;
    PartialMolecule operator[](const boost::tuple<MolIdentifier,Index> &viewidx) const;

    MoleculeGroup& operator+=(const Molecules &molecules);
    MoleculeGroup& operator-=(const Molecules &molecules);

    quint64 getMoleculeVersion(MolNum molnum) const;
    quint64 getMoleculeVersion(const MolID &molid) const;

    const ViewsOfMol& at(MolNum molnum) const;
    const ViewsOfMol& at(MolIdx molidx) const;
    const ViewsOfMol& at(const MolName &molname) const;
    const ViewsOfMol& at(const MolID &molid) const;

    PartialMolecule at(const boost::tuple<MolNum,Index> &viewidx) const;
    PartialMolecule at(const boost::tuple<MolIdentifier,Index> &viewidx) const;

    PartialMolecule at(MolNum molnum, int viewidx) const;
    PartialMolecule at(const MolID &molid, int viewidx) const;

    const ViewsOfMol& moleculeAt(int idx) const;
    PartialMolecule viewAt(int idx) const;

    MolNum molNumAt(int idx) const;
    const boost::tuple<MolNum,SireID::Index>& molViewIndexAt(int idx) const;

    int indexOf(const MoleculeView &molview) const;
    int indexOf(MolNum molnum) const;

    const QVector<MolNum>& molNums() const;
    const QVector< boost::tuple<MolNum,SireID::Index> >& molViewIndicies() const;

    const ViewsOfMol& molecule(MolNum molnum) const;
    const ViewsOfMol& molecule(MolIdx molidx) const;
    const ViewsOfMol& molecule(const MolName &molname) const;
    const ViewsOfMol& molecule(const MolID &molid) const;

    Molecules molecules(const MolID &molid) const;

    MolNum getMoleculeNumber(MolNum molnum) const;
    MolNum getMoleculeNumber(MolIdx molidx) const;
    MolNum getMoleculeNumber(const MolName &molname) const;
    MolNum getMoleculeNumber(const MolID &molid) const;

    QList<MolNum> map(MolNum molnum) const;
    QList<MolNum> map(MolIdx molidx) const;
    QList<MolNum> map(const MolName &molname) const;
    QList<MolNum> map(const MolID &molid) const;

    bool contains(MolNum molnum) const;
    bool contains(MolIdx molidx) const;
    bool contains(const MolName &molname) const;
    bool contains(const MolID &molid) const;

    bool contains(const MoleculeView &molview) const;
    bool contains(const ViewsOfMol &molviews) const;
    bool contains(const Molecules &molecules) const;
    bool contains(const MoleculeGroup &MoleculeGroup) const;

    bool intersects(const MoleculeView &molview) const;
    bool intersects(const Molecules &other) const;
    bool intersects(const MoleculeGroup &MoleculeGroup) const;

    int nAtoms() const;
    int nResidues() const;
    int nChains() const;
    int nSegments() const;
    int nMolecules() const;

    int nViews() const;

    int nViews(MolNum molnum) const;
    int nViews(const MolID &molid) const;

    int nViews(Index idx) const;

    bool isEmpty() const;

    const Molecules& molecules() const;

    const ViewsOfMol& first() const;
    const ViewsOfMol& last() const;

    const ViewsOfMol& front() const;
    const ViewsOfMol& back() const;

    const_iterator begin() const;
    const_iterator end() const;

    const_iterator constBegin() const;
    const_iterator constEnd() const;

    const_iterator find(MolNum molnum) const;
    const_iterator constFind(MolNum molnum) const;

    const_iterator find(const MolID &molid) const;
    const_iterator constFind(const MolID &molid) const;

    SelectResult search(const QString &search_term) const;

    QSet<MolName> molNames() const;

    void assertContains(MolNum molnum) const;
    void assertContains(const MolName &molname) const;

    const MGName& name() const;
    MGNum number() const;
    const Version& version() const;

    virtual QString toString() const;

    virtual void setName(const QString &new_name);
    virtual void setNumber(quint32 new_number);
    void setNewNumber();

    quint64 majorVersion() const;
    quint64 minorVersion() const;

    virtual void add(const MoleculeView &molview);
    virtual void add(const ViewsOfMol &molviews);
    virtual void add(const Molecules &molecules);
    virtual void add(const MoleculeGroup &MoleculeGroup);

    virtual bool addIfUnique(const MoleculeView &molview);
    virtual ViewsOfMol addIfUnique(const ViewsOfMol &molviews);
    virtual QList<ViewsOfMol> addIfUnique(const Molecules &molecules);
    virtual QList<ViewsOfMol> addIfUnique(const MoleculeGroup &MoleculeGroup);

    bool unite(const MoleculeView &molview);
    ViewsOfMol unite(const ViewsOfMol &molviews);
    QList<ViewsOfMol> unite(const Molecules &molecules);
    QList<ViewsOfMol> unite(const MoleculeGroup &MoleculeGroup);

    virtual bool remove(const MoleculeView &molview);
    virtual ViewsOfMol remove(const ViewsOfMol &molviews);
    virtual QList<ViewsOfMol> remove(const Molecules &molecules);
    virtual QList<ViewsOfMol> remove(const MoleculeGroup &MoleculeGroup);

    virtual bool removeAll(const MoleculeView &molview);
    virtual ViewsOfMol removeAll(const ViewsOfMol &molviews);
    virtual QList<ViewsOfMol> removeAll(const Molecules &molecules);
    virtual QList<ViewsOfMol> removeAll(const MoleculeGroup &MoleculeGroup);

    virtual ViewsOfMol remove(MolNum molnum);
    virtual QList<ViewsOfMol> remove(const QSet<MolNum> &molnums);

    virtual void removeAll();

    virtual bool update(const MoleculeData &moldata, bool auto_commit=true);
    bool update(const MoleculeView &molview, bool auto_commit=true);

    virtual QList<Molecule> update(const Molecules &molecules, bool auto_commit=true);
    virtual QList<Molecule> update(const MoleculeGroup &MoleculeGroup, bool auto_commit=true);

    virtual bool setContents(const MoleculeView &molview);
    virtual bool setContents(const ViewsOfMol &molviews);
    virtual bool setContents(const Molecules &molecules);
    virtual bool setContents(const MoleculeGroup &MoleculeGroup);

    virtual void accept();
    virtual bool needsAccepting() const;

    virtual int nFrames() const;
    virtual int nFrames(const SireBase::PropertyMap &map) const;

    virtual void loadFrame(int frame);
    virtual void saveFrame(int frame);
    virtual void saveFrame();
    virtual void deleteFrame(int frame);

    virtual void loadFrame(int frame, const SireBase::PropertyMap &map);
    virtual void saveFrame(int frame, const SireBase::PropertyMap &map);
    virtual void saveFrame(const SireBase::PropertyMap &map);
    virtual void deleteFrame(int frame, const SireBase::PropertyMap &map);

    static const MoleculeGroup& null();

private:
    bool _pvt_remove(const MoleculeView &molview);
    ViewsOfMol _pvt_remove(const ViewsOfMol &molviews);

    ViewsOfMol _pvt_remove(MolNum molnum);

    bool _pvt_removeAll(const MoleculeView &molview);
    ViewsOfMol _pvt_removeAll(const ViewsOfMol &molviews);

    void _pvt_setContents(const Molecules &molecules);

    /** Implicitly shared pointer to the contents and index
        of this group */
    SireBase::SharedDataPointer<detail::MolGroupPvt> d;

    /** The workspace used to cache updates, thus preventing
        excessive re-allocation of memory during, e.g. MC moves */
    MolGroupWorkspace workspace;
};

typedef SireBase::PropPtr<MoleculeGroup> MolGroupPtr;

}

Q_DECLARE_METATYPE(SireMol::MoleculeGroup);

SIRE_EXPOSE_CLASS( SireMol::MoleculeGroup )

SIRE_EXPOSE_PROPERTY( SireMol::MolGroupPtr, SireMol::MoleculeGroup )

SIRE_END_HEADER

#endif
