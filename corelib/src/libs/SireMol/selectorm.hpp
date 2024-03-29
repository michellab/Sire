/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2022  Christopher Woods
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

#ifndef SIREMOL_SELECTORM_HPP
#define SIREMOL_SELECTORM_HPP

#include "selectormol.h"

#include "core.h"

#include "SireBase/parallel.h"
#include "SireBase/booleanproperty.h"

#include "SireMol/errors.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
template<class T>
class SelectorM;

class SelectorMol;
}

template<class T>
SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::SelectorM<T>&);

template<class T>
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::SelectorM<T>&);

namespace SireMol
{

/** This is an analogue of the Selector<T> class that is designed
    to hold views from multiple molecules
*/
template<class T>
class SIREMOL_EXPORT SelectorM
    : public SireBase::ConcreteProperty<SelectorM<T>,SireBase::Property>
{

friend SIREMOL_EXPORT QDataStream& ::operator<<<>(QDataStream&, const SelectorM<T>&);
friend SIREMOL_EXPORT QDataStream& ::operator>><>(QDataStream&, SelectorM<T>&);

public:
    typedef typename QList< Selector<T> >::const_iterator iterator;
    typedef typename QList< Selector<T> >::const_iterator const_iterator;

    SelectorM();
    SelectorM(const T &view);
    SelectorM(const Selector<T> &views);
    SelectorM(const Molecules &mols);
    SelectorM(const MoleculeGroup &mols);
    SelectorM(const MolGroupsBase &mols);
    SelectorM(const SelectResult &mols);

    SelectorM(const SelectorMol &mols);
    SelectorM(const SelectorMol &mols, const SireBase::Slice &slice);
    SelectorM(const SelectorMol &mols, const QList<qint64> &idxs);
    SelectorM(const SelectorMol &mols, const QString &name,
              const SireBase::PropertyMap &map = SireBase::PropertyMap());
    SelectorM(const SelectorMol &mols, const typename T::ID &id);

    template<class U>
    SelectorM(const SelectorM<U> &other);
    template<class U>
    SelectorM(const SelectorM<U> &other, const SireBase::Slice &slice);
    template<class U>
    SelectorM(const SelectorM<U> &other, const QList<qint64> &idxs);
    template<class U>
    SelectorM(const SelectorM<U> &other, const QString &name);
    template<class U>
    SelectorM(const SelectorM<U> &other, const typename T::ID &id);

    SelectorM(const SelectorM<T> &other);

    virtual ~SelectorM();

    static const char* typeName();

    virtual SelectorM<T>* clone() const
    {
        return new SelectorM<T>(*this);
    }

    SelectorM<T>& operator=(const SelectorM<T> &other);

    bool operator==(const SelectorM<T> &other) const;
    bool operator!=(const SelectorM<T> &other) const;

    T operator[](int i) const;
    SelectorM<T> operator[](const SireBase::Slice &slice) const;
    SelectorM<T> operator[](const QList<qint64> &idxs) const;
    T operator[](const QString &name) const;
    T operator[](const typename T::ID &id) const;

    T operator()(int i) const;
    T operator()(const QString &name) const;
    T operator()(const typename T::ID &id) const;

    QList<MolViewPtr> toList() const;
    Molecules toMolecules() const;

    int count() const;
    int size() const;

    void update(const Molecules &molecules);

    EvaluatorM evaluate() const;

    MoleculeGroup toMoleculeGroup() const;
    SelectResult toSelectResult() const;

    bool isSelector() const;

    SelectorM<T> add(const SelectorM<T> &other) const;
    SelectorM<T> add(const Selector<T> &views) const;
    SelectorM<T> add(const T &view) const;

    SelectorM<T> subtract(const SelectorM<T> &other) const;
    SelectorM<T> subtract(const Selector<T> &views) const;
    SelectorM<T> subtract(const T &view) const;

    SelectorM<T> intersection(const SelectorM<T> &other) const;
    SelectorM<T> intersection(const Selector<T> &views) const;
    SelectorM<T> intersection(const T &view) const;

    SelectorM<T> invert() const;

    bool intersects(const SelectorM<T> &other) const;
    bool intersects(const Selector<T> &view) const;
    bool intersects(const T &view) const;

    bool contains(const SelectorM<T> &other) const;
    bool contains(const Selector<T> &view) const;
    bool contains(const T &view) const;

    Molecule molecule(int i,
                      const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
    Molecule molecule(const QString &name,
                      const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
    Molecule molecule(const MolID &molid,
                      const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SelectorMol molecules() const;
    SelectorMol molecules(int i,
                          const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
    SelectorMol molecules(const SireBase::Slice &slice,
                          const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
    SelectorMol molecules(const QList<qint64> &idxs,
                          const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
    SelectorMol molecules(const QString &name,
                          const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
    SelectorMol molecules(const MolID &molid,
                          const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;

    Atom atom(int i,
              const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
    Atom atom(const QString &name,
              const SireBase::PropertyMap &map=SireBase::PropertyMap()) const;
    Atom atom(const AtomID &atomid,
              const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;

    Residue residue(int i,
                    const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
    Residue residue(const QString &name,
                    const SireBase::PropertyMap &map=SireBase::PropertyMap()) const;
    Residue residue(const ResID &resid,
                    const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;

    Chain chain(int i,
                const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
    Chain chain(const QString &name,
                const SireBase::PropertyMap &map=SireBase::PropertyMap()) const;
    Chain chain(const ChainID &chainid,
                const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;

    Segment segment(int i,
                    const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
    Segment segment(const QString &name,
                    const SireBase::PropertyMap &map=SireBase::PropertyMap()) const;
    Segment segment(const SegID &segid,
                    const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;

    CutGroup cutGroup(int i,
                      const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
    CutGroup cutGroup(const QString &name,
                      const SireBase::PropertyMap &map=SireBase::PropertyMap()) const;
    CutGroup cutGroup(const CGID &cgid,
                      const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;

    SelectorM<Atom> atoms() const;
    SelectorM<Atom> atoms(int i,
                          const SireBase::PropertyMap &map=SireBase::PropertyMap()) const;
    SelectorM<Atom> atoms(const SireBase::Slice &slice,
                          const SireBase::PropertyMap &map=SireBase::PropertyMap()) const;
    SelectorM<Atom> atoms(const QList<qint64> &idxs,
                          const SireBase::PropertyMap &map=SireBase::PropertyMap()) const;
    SelectorM<Atom> atoms(const QString &name,
                          const SireBase::PropertyMap &map=SireBase::PropertyMap()) const;
    SelectorM<Atom> atoms(const AtomID &atomid,
                          const SireBase::PropertyMap &map=SireBase::PropertyMap()) const;

    SelectorM<Residue> residues() const;
    SelectorM<Residue> residues(int i,
                                const SireBase::PropertyMap &map=SireBase::PropertyMap()) const;
    SelectorM<Residue> residues(const SireBase::Slice &slice,
                                const SireBase::PropertyMap &map=SireBase::PropertyMap()) const;
    SelectorM<Residue> residues(const QList<qint64> &idxs,
                                const SireBase::PropertyMap &map=SireBase::PropertyMap()) const;
    SelectorM<Residue> residues(const QString &name,
                                const SireBase::PropertyMap &map=SireBase::PropertyMap()) const;
    SelectorM<Residue> residues(const ResID &resid,
                                const SireBase::PropertyMap &map=SireBase::PropertyMap()) const;

    SelectorM<Chain> chains() const;
    SelectorM<Chain> chains(int i,
                            const SireBase::PropertyMap &map=SireBase::PropertyMap()) const;
    SelectorM<Chain> chains(const SireBase::Slice &slice,
                            const SireBase::PropertyMap &map=SireBase::PropertyMap()) const;
    SelectorM<Chain> chains(const QList<qint64> &idxs,
                            const SireBase::PropertyMap &map=SireBase::PropertyMap()) const;
    SelectorM<Chain> chains(const QString &name,
                            const SireBase::PropertyMap &map=SireBase::PropertyMap()) const;
    SelectorM<Chain> chains(const ChainID &chainid,
                            const SireBase::PropertyMap &map=SireBase::PropertyMap()) const;

    SelectorM<Segment> segments() const;
    SelectorM<Segment> segments(int i,
                                const SireBase::PropertyMap &map=SireBase::PropertyMap()) const;
    SelectorM<Segment> segments(const SireBase::Slice &slice,
                                const SireBase::PropertyMap &map=SireBase::PropertyMap()) const;
    SelectorM<Segment> segments(const QList<qint64> &idxs,
                                const SireBase::PropertyMap &map=SireBase::PropertyMap()) const;
    SelectorM<Segment> segments(const QString &name,
                                const SireBase::PropertyMap &map=SireBase::PropertyMap()) const;
    SelectorM<Segment> segments(const SegID &segid,
                                const SireBase::PropertyMap &map=SireBase::PropertyMap()) const;

    SelectorM<CutGroup> cutGroups() const;
    SelectorM<CutGroup> cutGroups(int i,
                                  const SireBase::PropertyMap &map=SireBase::PropertyMap()) const;
    SelectorM<CutGroup> cutGroups(const SireBase::Slice &slice,
                                  const SireBase::PropertyMap &map=SireBase::PropertyMap()) const;
    SelectorM<CutGroup> cutGroups(const QList<qint64> &idxs,
                                  const SireBase::PropertyMap &map=SireBase::PropertyMap()) const;
    SelectorM<CutGroup> cutGroups(const QString &name,
                                  const SireBase::PropertyMap &map=SireBase::PropertyMap()) const;
    SelectorM<CutGroup> cutGroups(const CGID &cgid,
                                  const SireBase::PropertyMap &map=SireBase::PropertyMap()) const;

    SelectResult search(const QString &search_string) const;
    SelectResult search(const QString &search_string,
                        const SireBase::PropertyMap &map) const;

    QList<typename T::Index> IDs() const;
    QList<typename T::Index> indexes() const;
    QList<typename T::Number> numbers() const;
    QList<typename T::Name> names() const;

    int nAtoms() const;
    int nResidues() const;
    int nChains() const;
    int nSegments() const;
    int nCutGroups() const;
    int nMolecules() const;

    bool isEmpty() const;

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

    const_iterator begin() const;
    const_iterator end() const;

    const_iterator constBegin() const;
    const_iterator constEnd() const;

    virtual QString toString() const;

protected:
    void _append(const T &view);

    /** The actual views */
    QList< Selector<T> > vws;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

#ifndef GCCXML_PARSE
namespace detail
{

template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool _usesParallel(const QList<T> vws, const SireBase::PropertyMap &map)
{
    if (vws.count() < 16)
        return false;

    else if (map["parallel"].hasValue())
    {
        return map["parallel"].value().asA<SireBase::BooleanProperty>().value();
    }

    return true;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
int _nFrames(const QList<T> &vws, const SireBase::PropertyMap &map)
{
    int nframes = std::numeric_limits<int>::max();

    if (_usesParallel(vws, map))
    {
        nframes = tbb::parallel_reduce( tbb::blocked_range<int>(0, vws.count()),
                                        std::numeric_limits<int>::max(),
                                        [&](tbb::blocked_range<int> r, int my_nframes)
        {
            for (int i=r.begin(); i<r.end(); ++i)
            {
                my_nframes = std::min(my_nframes, vws.at(i).nFrames(map));
            }

            return my_nframes;

        }, [](int a, int b){ return std::min(a, b);});
    }
    else
    {
        for (const auto &view : vws)
        {
            nframes = std::min(nframes, view.nFrames(map));

            if (nframes == 0)
                break;
        }
    }

    return nframes;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
void _loadFrame(QList<T> &vws, int frame, const SireBase::PropertyMap &map)
{
    const int nframes = _nFrames(vws, map);

    if (not (nframes == 0 and frame == 0))
        frame = SireID::Index(frame).map(nframes);

    vws.detach();
    const int n = vws.count();

    if (_usesParallel(vws, map))
    {
        tbb::parallel_for(tbb::blocked_range<int>(0, n),
                          [&](tbb::blocked_range<int> r)
        {
            for (int i=r.begin(); i<r.end(); ++i)
            {
                vws[i].loadFrame(frame, map);
            }
        });
    }
    else
    {
        for (int i=0; i<n; ++i)
        {
            vws[i].loadFrame(frame, map);
        }
    }
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
void _saveFrame(QList<T> &vws, int frame, const SireBase::PropertyMap &map)
{
    frame = SireID::Index(frame).map(_nFrames(vws, map));

    vws.detach();
    const int n = vws.count();

    if (_usesParallel(vws, map))
    {
        tbb::parallel_for(tbb::blocked_range<int>(0, n),
                          [&](tbb::blocked_range<int> r)
        {
            for (int i=r.begin(); i<r.end(); ++i)
            {
                vws[i].saveFrame(frame, map);
            }
        });
    }
    else
    {
        for (int i=0; i<n; ++i)
        {
            vws[i].saveFrame(frame, map);
        }
    }
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
void _saveFrame(QList<T> &vws, const SireBase::PropertyMap &map)
{
    vws.detach();
    const int n = vws.count();

    if (_usesParallel(vws, map))
    {
        tbb::parallel_for(tbb::blocked_range<int>(0, n),
                          [&](tbb::blocked_range<int> r)
        {
            for (int i=r.begin(); i<r.end(); ++i)
            {
                vws[i].saveFrame(map);
            }
        });
    }
    else
    {
        for (int i=0; i<n; ++i)
        {
            vws[i].saveFrame(map);
        }
    }
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
void _deleteFrame(QList<T> &vws, int frame, const SireBase::PropertyMap &map)
{
    frame = SireID::Index(frame).map(_nFrames(vws, map));

    vws.detach();
    const int n = vws.count();

    if (_usesParallel(vws, map))
    {
        tbb::parallel_for(tbb::blocked_range<int>(0, n),
                          [&](tbb::blocked_range<int> r)
        {
            for (int i=r.begin(); i<r.end(); ++i)
            {
                vws[i].deleteFrame(frame, map);
            }
        });
    }
    else
    {
        for (int i=0; i<n; ++i)
        {
            vws[i].deleteFrame(frame, map);
        }
    }
}

} // end of namespace detail
#endif // GCCXML_PARSE


template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<T>::SelectorM()
             : SireBase::ConcreteProperty<SelectorM<T>,SireBase::Property>()
{}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<T>::SelectorM(const T &view)
             : SireBase::ConcreteProperty<SelectorM<T>,SireBase::Property>()
{
    this->vws.append(Selector<T>(view));
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<T>::SelectorM(const Selector<T> &views)
             : SireBase::ConcreteProperty<SelectorM<T>,SireBase::Property>()
{
    if (not views.isEmpty())
        this->vws.append(views);
}

template<class T>
struct _get_view
{
    template<class C>
    static int count(const C &mols)
    {
        return mols.nAtoms();
    }

    template<class C>
    static T at(const C &mols, int idx)
    {
        return mols.atom(idx);
    }

    template<class C>
    static Selector<Atom> get(const C &view)
    {
        if (view.template isA< Selector<Atom> >())
            return view.template asA< Selector<Atom> >();
        else
            return view.atoms();
    }

    static SelectorM<T> getName(const SelectorM<T> &view, const QString &name)
    {
        return view.atoms(name);
    }

    static SelectorM<T> getID(const SelectorM<T> &view,
                              const typename T::ID &id)
    {
        return view.atoms(id);
    }

    template<class C, class ID>
    static Selector<T> get(const C &view, const ID &id)
    {
        return view.atoms(id);
    }

    static void raise_duplicate(const QString &id, int n)
    {
        throw SireMol::duplicate_atom(QObject::tr(
            "Multiple atoms matched %1. Number of matches is %2.")
                .arg(id).arg(n), CODELOC);
    }

    static void raise_missing(const QString &id)
    {
        throw SireMol::missing_atom(QObject::tr(
            "No atom matches %1.").arg(id), CODELOC);
    }
};

template<>
struct _get_view<Residue>
{
    template<class C>
    static int count(const C &mols)
    {
        return mols.nResidues();
    }

    template<class C>
    static Residue at(const C &mols, int idx)
    {
        return mols.residue(idx);
    }

    static SelectorM<Residue> getName(const SelectorM<Residue> &view,
                                      const QString &name)
    {
        return view.residues(name);
    }

    static SelectorM<Residue> getID(const SelectorM<Residue> &view,
                                    const Residue::ID &id)
    {
        return view.residues(id);
    }

    template<class C>
    static Selector<Residue> get(const C &view)
    {
        if (view.template isA< Selector<Residue> >())
            return view.template asA< Selector<Residue> >();
        else
            return view.residues();
    }

    template<class C, class ID>
    static Selector<Residue> get(const C &view, const ID &id)
    {
        return view.residues(id);
    }

    static void raise_duplicate(const QString &id, int n)
    {
        throw SireMol::duplicate_residue(QObject::tr(
            "Multiple residues matched %1. Number of matches is %2.")
                .arg(id).arg(n), CODELOC);
    }

    static void raise_missing(const QString &id)
    {
        throw SireMol::missing_residue(QObject::tr(
            "No residue matches %1.").arg(id), CODELOC);
    }
};

template<>
struct _get_view<Chain>
{
    template<class C>
    static int count(const C &mols)
    {
        return mols.nChains();
    }

    template<class C>
    static Chain at(const C &mols, int idx)
    {
        return mols.chain(idx);
    }

    static SelectorM<Chain> getName(const SelectorM<Chain> &view,
                                    const QString &name)
    {
        return view.chains(name);
    }

    static SelectorM<Chain> getID(const SelectorM<Chain> &view,
                                  const Chain::ID &id)
    {
        return view.chains(id);
    }

    template<class C>
    static Selector<Chain> get(const C &view)
    {
        if (view.template isA< Selector<Chain> >())
            return view.template asA< Selector<Chain> >();
        else
            return view.chains();
    }

    template<class C, class ID>
    static Selector<Chain> get(const C &view, const ID &id)
    {
        return view.chains(id);
    }

    static void raise_duplicate(const QString &id, int n)
    {
        throw SireMol::duplicate_chain(QObject::tr(
            "Multiple chains matched %1. Number of matches is %2.")
                .arg(id).arg(n), CODELOC);
    }

    static void raise_missing(const QString &id)
    {
        throw SireMol::missing_chain(QObject::tr(
            "No chain matches %1.").arg(id), CODELOC);
    }
};

template<>
struct _get_view<CutGroup>
{
    template<class C>
    static int count(const C &mols)
    {
        return mols.nCutGroups();
    }

    template<class C>
    static CutGroup at(const C &mols, int idx)
    {
        return mols.cutGroup(idx);
    }

    static SelectorM<CutGroup> getName(const SelectorM<CutGroup> &view,
                                       const QString &name)
    {
        return view.cutGroup(name);
    }

    static SelectorM<CutGroup> getID(const SelectorM<Residue> &view,
                                     const CutGroup::ID &id)
    {
        return view.cutGroups(id);
    }

    template<class C>
    static Selector<CutGroup> get(const C &view)
    {
        if (view.template isA< Selector<CutGroup> >())
            return view.template asA< Selector<CutGroup> >();
        else
            return view.cutGroups();
    }

    template<class C, class ID>
    static Selector<CutGroup> get(const C &view, const ID &id)
    {
        return view.cutGroups(id);
    }

    static void raise_duplicate(const QString &id, int n)
    {
        throw SireMol::duplicate_cutgroup(QObject::tr(
            "Multiple CutGroups matched %1. Number of matches is %2.")
                .arg(id).arg(n), CODELOC);
    }

    static void raise_missing(const QString &id)
    {
        throw SireMol::missing_cutgroup(QObject::tr(
            "No CutGroup matches %1.").arg(id), CODELOC);
    }
};

template<>
struct _get_view<Segment>
{
    template<class C>
    static int count(const C &mols)
    {
        return mols.nSegments();
    }

    template<class C>
    static Segment at(const C &mols, int idx)
    {
        return mols.segment(idx);
    }

    static SelectorM<Segment> getName(const SelectorM<Segment> &view,
                                      const QString &name)
    {
        return view.segments(name);
    }

    static SelectorM<Segment> getID(const SelectorM<Segment> &view,
                                    const Segment::ID &id)
    {
        return view.segments(id);
    }

    template<class C>
    static Selector<Segment> get(const C &view)
    {
        if (view.template isA< Selector<Segment> >())
            return view.template asA< Selector<Segment> >();
        else
            return view.segments();
    }

    template<class C, class ID>
    static Selector<Segment> get(const C &view, const ID &id)
    {
        return view.segments(id);
    }

    static void raise_duplicate(const QString &id, int n)
    {
        throw SireMol::duplicate_segment(QObject::tr(
            "Multiple segments matched %1. Number of matches is %2.")
                .arg(id).arg(n), CODELOC);
    }

    static void raise_missing(const QString &id)
    {
        throw SireMol::missing_segment(QObject::tr(
            "No segment matches %1.").arg(id), CODELOC);
    }
};

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<T>::SelectorM(const SelectorMol &mols)
             : SireBase::ConcreteProperty<SelectorM<T>,SireBase::Property>()
{
    for (const auto &mol : mols)
    {
        this->vws += _get_view<T>::get(mol);
    }
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
void SelectorM<T>::_append(const T &view)
{
    if (this->vws.isEmpty())
    {
        this->vws.append(Selector<T>(view));
    }
    else if (this->vws.last().data().number() != view.data().number())
    {
        // new molecule
        this->vws.append(Selector<T>(view));
    }
    else
    {
        // a new view in the current molecule
        this->vws.last() = this->vws.last().add(view);
    }
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<T>::SelectorM(const Molecules &molecules)
             : SireBase::ConcreteProperty<SelectorM<T>,SireBase::Property>()
{
    if (not molecules.isEmpty())
    {
        auto toList = [](const QSet<MolNum> &molnums)
        {
            //return QList<MolNum>(molnums.constBegin(), molnums.constEnd());
            return molnums.values();
        };

        auto molnums = toList(molecules.molNums());

        //sort them, as this is also likely the order the molecules
        //were read in from a file, and so more likely to be the
        //order the user would expect
        std::sort(molnums.begin(), molnums.end());

        this->vws.reserve(molnums.count());

        for (const auto &molnum : molnums)
        {
            this->vws.append(_get_view<T>::get(molecules.at(molnum)));
        }
    }
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<T>::SelectorM(const MoleculeGroup &molecules)
             : SireBase::ConcreteProperty<SelectorM<T>,SireBase::Property>()
{
    if (not molecules.isEmpty())
    {
        const auto molnums = molecules.molNums();
        this->vws.reserve(molnums.count());

        for (const auto &molnum : molnums)
        {
            this->vws.append(_get_view<T>::get(molecules.at(molnum)));
        }
    }
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<T>::SelectorM(const MolGroupsBase &molecules)
             : SireBase::ConcreteProperty<SelectorM<T>,SireBase::Property>()
{
    if (not molecules.isEmpty())
    {
        const auto molnums = molecules.molNums();
        this->vws.reserve(molnums.count());

        for (const auto &molnum : molnums)
        {
            this->vws.append(_get_view<T>::get(molecules.at(molnum)));
        }
    }
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<T>::SelectorM(const SelectResult &molecules)
             : SireBase::ConcreteProperty<SelectorM<T>,SireBase::Property>()
{
    if (not molecules.isEmpty())
    {
        this->vws.reserve(molecules.count());

        for (const auto &mol : molecules)
        {
            this->vws.append(_get_view<T>::get(*mol));
        }
    }
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<T>::SelectorM(const SelectorMol &mols, const SireBase::Slice &slice)
             : SireBase::ConcreteProperty<SelectorM<T>,SireBase::Property>()
{
    for (auto it = slice.begin(_get_view<T>::count(mols));
         not it.atEnd(); it.next())
    {
        this->_append(_get_view<T>::at(mols, it.value()));
    }
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<T>::SelectorM(const SelectorMol &mols, const QList<qint64> &idxs)
             : SireBase::ConcreteProperty<SelectorM<T>,SireBase::Property>()
{
    for (const auto &idx : idxs)
    {
        this->_append(_get_view<T>::at(mols, idx));
    }
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<T>::SelectorM(const SelectorMol &mols, const QString &name,
                        const SireBase::PropertyMap &map)
             : SireBase::ConcreteProperty<SelectorM<T>,SireBase::Property>()
{
    for (const auto &mol : mols)
    {
        auto flag = SireError::exception::enableFastExceptions();

        try
        {
            this->vws.append(_get_view<T>::get(mol, typename T::Name(name)));
        }
        catch(...)
        {}
    }

    if (this->vws.isEmpty())
    {
        // try a search
        try
        {
            this->operator=(SelectorM<T>(mols.search(name, map)));
        }
        catch(...)
        {
            if (name.length() < 5)
                //likely a name error
                _get_view<T>::raise_missing(name);
            else
                //likely a syntax error
                throw;
        }
    }
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<T>::SelectorM(const SelectorMol &mols, const typename T::ID &id)
             : SireBase::ConcreteProperty<SelectorM<T>,SireBase::Property>()
{
    for (const auto &mol : mols)
    {
        auto flag = SireError::exception::enableFastExceptions();

        try
        {
            this->vws.append(_get_view<T>::get(mol, id));
        }
        catch(...)
        {}
    }

    if (this->vws.isEmpty())
        _get_view<T>::raise_missing(id.toString());
}

template<class T>
template<class U>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<T>::SelectorM(const SelectorM<U> &other)
             : SireBase::ConcreteProperty<SelectorM<T>,SireBase::Property>()
{
    for (const auto &o : other)
    {
        this->vws.append(_get_view<T>::get(o));
    }
}

template<class T>
template<class U>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<T>::SelectorM(const SelectorM<U> &other, const SireBase::Slice &slice)
             : SireBase::ConcreteProperty<SelectorM<T>,SireBase::Property>()
{
    for (auto it = slice.begin(_get_view<T>::count(other));
         not it.atEnd(); it.next())
    {
        this->_append(_get_view<T>::at(other, it.value()));
    }
}

template<class T>
template<class U>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<T>::SelectorM(const SelectorM<U> &other, const QList<qint64> &idxs)
             : SireBase::ConcreteProperty<SelectorM<T>,SireBase::Property>()
{
    for (const auto &idx : idxs)
    {
        this->_append(_get_view<T>::at(other, idx));
    }
}

template<class T>
template<class U>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<T>::SelectorM(const SelectorM<U> &other, const QString &name)
             : SireBase::ConcreteProperty<SelectorM<T>,SireBase::Property>()
{
    for (const auto &view : other)
    {
        auto flag = SireError::exception::enableFastExceptions();

        try
        {
            this->vws.append(_get_view<T>::get(view, name));
        }
        catch(...)
        {}
    }

    if (this->vws.isEmpty())
    {
        // try a search
        try
        {
            this->operator=(SelectorM<T>(other.search(name)));
        }
        catch(...)
        {
            if (name.length() < 5)
                //likely a name error
                _get_view<T>::raise_missing(name);
            else
                //likely a syntax error
                throw;
        }
    }
}

template<class T>
template<class U>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<T>::SelectorM(const SelectorM<U> &other, const typename T::ID &id)
             : SireBase::ConcreteProperty<SelectorM<T>,SireBase::Property>()
{
    for (const auto &view : other)
    {
        auto flag = SireError::exception::enableFastExceptions();

        try
        {
            this->vws.append(_get_view<T>::get(view, id));
        }
        catch(...)
        {}
    }

    if (this->vws.isEmpty())
        _get_view<T>::raise_missing(id.toString());
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<T>::SelectorM(const SelectorM<T> &other)
             : SireBase::ConcreteProperty<SelectorM<T>,SireBase::Property>(),
               vws(other.vws)
{}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<T>::~SelectorM()
{}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
const char* SelectorM<T>::typeName()
{
    return QMetaType::typeName( qMetaTypeId< SelectorM<T> >() );
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<T>& SelectorM<T>::operator=(const SelectorM<T> &other)
{
    if (this != &other)
    {
        this->vws = other.vws;
        Property::operator=(other);
    }

    return *this;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool SelectorM<T>::operator==(const SelectorM<T> &other) const
{
    return vws == other.vws;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool SelectorM<T>::operator!=(const SelectorM<T> &other) const
{
    return not this->operator==(other);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
T SelectorM<T>::operator[](int i) const
{
    i = SireID::Index(i).map(this->count());

    for (const auto &v : vws)
    {
        if (i < v.count())
        {
            return v(i);
        }
        else
        {
            i -= v.count();
        }
    }

    throw SireError::program_bug(QObject::tr("Should not get here!"), CODELOC);

    return T();
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<T> SelectorM<T>::operator[](const SireBase::Slice &slice) const
{
    return SelectorM<T>(*this, slice);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<T> SelectorM<T>::operator[](const QList<qint64> &idxs) const
{
    return SelectorM<T>(*this,idxs);
}


template<class T>
SIRE_OUTOFLINE_TEMPLATE
T SelectorM<T>::operator[](const QString &name) const
{
    auto all = _get_view<T>::getName(*this, name);

    if (all.count() > 1)
    {
        _get_view<T>::raise_duplicate(name, all.count());
    }

    BOOST_ASSERT( not all.isEmpty() );

    return all(0);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
T SelectorM<T>::operator[](const typename T::ID &id) const
{
    auto all = _get_view<T>::getID(*this, id);

    if (all.count() > 1)
    {
        _get_view<T>::raise_duplicate(id.toString(), all.count());
    }

    BOOST_ASSERT( not all.isEmpty() );

    return all(0);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
T SelectorM<T>::operator()(int i) const
{
    return this->operator[](i);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
T SelectorM<T>::operator()(const QString &name) const
{
    return this->operator[](name);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
T SelectorM<T>::operator()(const typename T::ID &id) const
{
    return this->operator[](id);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
QList<MolViewPtr> SelectorM<T>::toList() const
{
    QList<MolViewPtr> l;
    l.reserve(vws.count());

    for (const auto &v : vws)
    {
        l.append(MolViewPtr(v.clone()));
    }

    return l;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
Molecules SelectorM<T>::toMolecules() const
{
    return Molecules(this->vws);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectResult SelectorM<T>::toSelectResult() const
{
    return SelectResult(this->vws);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool SelectorM<T>::isSelector() const
{
    return true;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<T> SelectorM<T>::add(const SelectorM<T> &other) const
{
    SelectorM<T> ret(*this);

    for (int i=0; i<other.count(); ++i)
    {
        ret._append(other(i));
    }

    return ret;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<T> SelectorM<T>::add(const Selector<T> &other) const
{
    SelectorM<T> ret(*this);

    for (int i=0; i<other.count(); ++i)
    {
        ret._append(other(i));
    }

    return ret;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<T> SelectorM<T>::add(const T &other) const
{
    SelectorM<T> ret(*this);
    ret._append(other);
    return ret;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<T> SelectorM<T>::subtract(const SelectorM<T> &other) const
{
    SelectorM<T> ret;

    for (const auto &view : this->vws)
    {
        auto result = view;

        for (const auto &other_view : other.vws)
        {
            if (result.data().number() == other_view.data().number())
            {
                result = result.subtract(other_view);
            }
        }

        if (not result.isEmpty())
            ret.vws.append(result);
    }

    return ret;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<T> SelectorM<T>::subtract(const Selector<T> &other) const
{
    SelectorM<T> ret;

    for (const auto &view : this->vws)
    {
        if (view.data().number() == other.data().number())
        {
            auto result = view.subtract(other);

            if (not result.isEmpty())
                ret.vws.append(result);
        }
        else
            ret.vws.append(view);
    }

    return ret;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<T> SelectorM<T>::subtract(const T &other) const
{
    SelectorM<T> ret;

    for (const auto &view : this->vws)
    {
        if (view.data().number() == other.data().number())
        {
            auto result = view.subtract(other);

            if (not result.isEmpty())
                ret.vws.append(result);
        }
        else
            ret.vws.append(view);
    }

    return ret;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<T> SelectorM<T>::intersection(const SelectorM<T> &other) const
{
    SelectorM<T> ret;

    for (const auto &view : this->vws)
    {
        for (const auto &other_view : other.vws)
        {
            if (view.data().number() == other_view.data().number())
            {
                auto intersect = view.intersection(other_view);

                if (not intersect.isEmpty())
                {
                    ret.vws.append(intersect);
                }
            }
        }
    }

    return ret;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<T> SelectorM<T>::intersection(const Selector<T> &views) const
{
    SelectorM<T> ret;

    for (const auto &view : this->vws)
    {
        if (view.data().number() == view.data().number())
        {
            auto intersect = view.intersection(views);

            if (not intersect.isEmpty())
            {
                ret.vws.append(intersect);
            }
        }
    }

    return ret;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<T> SelectorM<T>::intersection(const T &view) const
{
    if (this->contains(view))
    {
        return SelectorM<T>(view);
    }
    else
        return SelectorM<T>();
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<T> SelectorM<T>::invert() const
{
    SelectorM<T> ret;

    for (const auto &view : this->vws)
    {
        auto inverted = view.invert();

        if (not inverted.isEmpty())
        {
            ret.vws.append(inverted);
        }
    }

    return ret;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool SelectorM<T>::intersects(const SelectorM<T> &other) const
{
    for (int i=0; i<other.count(); ++i)
    {
        if (this->contains(other(i)))
        {
            return true;
        }
    }

    return false;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool SelectorM<T>::intersects(const Selector<T> &other) const
{
    for (int i=0; i<other.count(); ++i)
    {
        if (this->contains(other(i)))
        {
            return true;
        }
    }

    return false;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool SelectorM<T>::intersects(const T &view) const
{
    return this->contains(view);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool SelectorM<T>::contains(const SelectorM<T> &other) const
{
    for (int i=0; i<other.count(); ++i)
    {
        if (not this->contains(other(i)))
        {
            return false;
        }
    }

    return true;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool SelectorM<T>::contains(const Selector<T> &other) const
{
    for (int i=0; i<other.count(); ++i)
    {
        if (not this->contains(other(i)))
        {
            return false;
        }
    }

    return true;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool SelectorM<T>::contains(const T &other) const
{
    for (const auto &view : this->vws)
    {
        if (view == other)
            return true;
    }

    return false;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
void SelectorM<T>::update(const Molecules &molecules)
{
    // better to create a map from MolNum to index here
    QMultiHash<MolNum, int> molnum_to_idx;
    molnum_to_idx.reserve(this->vws.count());

    int i = 0;

    for (const auto &view : this->vws)
    {
        molnum_to_idx.insert(view.data().number(), i);
        i += 1;
    }

    for (const auto &mol : molecules)
    {
        const auto molnum = mol.data().number();

        auto it = molnum_to_idx.constFind(molnum);

        while (it != molnum_to_idx.constEnd() && it.key() == molnum)
        {
            this->vws[it.value()].update(mol.data());
        }
    }
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectResult SelectorM<T>::search(const QString &search_term) const
{
    Select search(search_term);
    return search(this->toSelectResult());
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectResult SelectorM<T>::search(const QString &search_term,
                                  const SireBase::PropertyMap &map) const
{
    Select search(search_term);
    return search(this->toSelectResult(), map);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
int SelectorM<T>::count() const
{
    int n = 0;

    for (const auto &v : this->vws)
    {
        n += v.count();
    }

    return n;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
int SelectorM<T>::size() const
{
    return this->count();
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
MoleculeGroup SelectorM<T>::toMoleculeGroup() const
{
    MoleculeGroup grp;

    for (const auto &view : this->vws)
    {
        grp.add(view);
    }

    return grp;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
Molecule SelectorM<T>::molecule(int i,
                                const SireBase::PropertyMap&) const
{
    i = SireID::Index(i).map(this->vws.count());

    return this->vws.at(i).molecule();
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
Molecule SelectorM<T>::molecule(const QString &name,
                                const SireBase::PropertyMap &map) const
{
    auto mols = this->molecules(name, map);

    if (mols.count() > 1)
        throw SireMol::duplicate_molecule(QObject::tr(
            "More than one molecule matches '%1'. Number of matches is %2.")
                .arg(name).arg(mols.count()), CODELOC);

    BOOST_ASSERT(not mols.isEmpty());

    return mols[0];
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
Molecule SelectorM<T>::molecule(const MolID &molid,
                                const SireBase::PropertyMap &map)
{
    auto mols = this->molecules(molid, map);

    if (mols.count() > 1)
        throw SireMol::duplicate_molecule(QObject::tr(
            "More than one molecule matches '%1'. Number of matches is %2.")
                .arg(molid.toString()).arg(mols.count()), CODELOC);

    BOOST_ASSERT(not mols.isEmpty());

    return mols[0];
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorMol SelectorM<T>::molecules() const
{
    return SelectorMol(*this);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorMol SelectorM<T>::molecules(int i,
                                    const SireBase::PropertyMap &map) const
{
    return SelectorMol(this->molecule(i, map));
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorMol SelectorM<T>::molecules(const SireBase::Slice &slice,
                                    const SireBase::PropertyMap&) const
{
    return SelectorMol(*this, slice);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorMol SelectorM<T>::molecules(const QList<qint64> &idxs,
                                    const SireBase::PropertyMap&) const
{
    return SelectorMol(*this, idxs);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorMol SelectorM<T>::molecules(const QString &name,
                                    const SireBase::PropertyMap &map) const
{
    return SelectorMol(*this, name, map);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorMol SelectorM<T>::molecules(const MolID &molid,
                                    const SireBase::PropertyMap&) const
{
    return SelectorMol(*this, molid);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
Atom SelectorM<T>::atom(int i, const SireBase::PropertyMap&) const
{
    i = SireID::Index(i).map(this->nAtoms());

    for (const auto &v : this->vws)
    {
        if (i < v.nAtoms())
        {
            return v.atom(i);
        }
        else
        {
            i -= v.nAtoms();
        }
    }

    throw SireError::program_bug(QObject::tr("Should not get here!"), CODELOC);

    return Atom();
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
Atom SelectorM<T>::atom(const QString &name,
                        const SireBase::PropertyMap &map) const
{
    auto all = this->atoms(name, map);

    if (all.count() > 1)
        _get_view<Atom>::raise_duplicate(name, all.count());

    BOOST_ASSERT( not all.isEmpty() );

    return all[0];
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
Atom SelectorM<T>::atom(const AtomID &atomid,
                        const SireBase::PropertyMap &map) const
{
    auto all = this->atoms(atomid, map);

    if (all.count() > 1)
        _get_view<Atom>::raise_duplicate(atomid.toString(), all.count());

    BOOST_ASSERT( not all.isEmpty() );

    return all[0];
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
Residue SelectorM<T>::residue(int i, const SireBase::PropertyMap&) const
{
    i = SireID::Index(i).map(this->nResidues());

    for (const auto &v : this->vws)
    {
        if (i < v.nResidues())
        {
            return v.residue(i);
        }
        else
        {
            i -= v.nResidues();
        }
    }

    throw SireError::program_bug(QObject::tr("Should not get here!"), CODELOC);

    return Residue();
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
Residue SelectorM<T>::residue(const QString &name,
                              const SireBase::PropertyMap &map) const
{
    auto all = this->residues(name, map);

    if (all.count() > 1)
        _get_view<Residue>::raise_duplicate(name, all.count());

    BOOST_ASSERT( not all.isEmpty() );

    return all[0];
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
Residue SelectorM<T>::residue(const ResID &resid,
                              const SireBase::PropertyMap &map) const
{
    auto all = this->residues(resid, map);

    if (all.count() > 1)
        _get_view<Residue>::raise_duplicate(resid.toString(), all.count());

    BOOST_ASSERT( not all.isEmpty() );

    return all[0];
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
Chain SelectorM<T>::chain(int i, const SireBase::PropertyMap&) const
{
    i = SireID::Index(i).map(this->nChains());

    for (const auto &v : this->vws)
    {
        if (i < v.nChains())
        {
            return v.chain(i);
        }
        else
        {
            i -= v.nChains();
        }
    }

    throw SireError::program_bug(QObject::tr("Should not get here!"), CODELOC);

    return Chain();
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
Chain SelectorM<T>::chain(const QString &name,
                          const SireBase::PropertyMap &map) const
{
    auto all = this->chains(name, map);

    if (all.count() > 1)
        _get_view<Chain>::raise_duplicate(name, all.count());

    BOOST_ASSERT( not all.isEmpty() );

    return all[0];
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
Chain SelectorM<T>::chain(const ChainID &chainid,
                          const SireBase::PropertyMap &map) const
{
    auto all = this->chains(chainid);

    if (all.count() > 1)
        _get_view<Chain>::raise_duplicate(chainid.toString(), all.count());

    BOOST_ASSERT( not all.isEmpty() );

    return all[0];
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
Segment SelectorM<T>::segment(int i, const SireBase::PropertyMap&) const
{
    i = SireID::Index(i).map(this->nSegments());

    for (const auto &v : this->vws)
    {
        if (i < v.nSegments())
        {
            return v.segment(i);
        }
        else
        {
            i -= v.nSegments();
        }
    }

    throw SireError::program_bug(QObject::tr("Should not get here!"), CODELOC);

    return Segment();
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
Segment SelectorM<T>::segment(const QString &name,
                              const SireBase::PropertyMap &map) const
{
    auto all = this->segments(name, map);

    if (all.count() > 1)
        _get_view<Segment>::raise_duplicate(name, all.count());

    BOOST_ASSERT( not all.isEmpty() );

    return all[0];
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
Segment SelectorM<T>::segment(const SegID &segid,
                              const SireBase::PropertyMap &map) const
{
    auto all = this->segments(segid, map);

    if (all.count() > 1)
        _get_view<Segment>::raise_duplicate(segid.toString(), all.count());

    BOOST_ASSERT( not all.isEmpty() );

    return all[0];
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
CutGroup SelectorM<T>::cutGroup(int i, const SireBase::PropertyMap&) const
{
    i = SireID::Index(i).map(this->nCutGroups());

    for (const auto &v : this->vws)
    {
        if (i < v.nCutGroups())
        {
            return v.cutGroup(i);
        }
        else
        {
            i -= v.nCutGroups();
        }
    }

    throw SireError::program_bug(QObject::tr("Should not get here!"), CODELOC);

    return CutGroup();
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
CutGroup SelectorM<T>::cutGroup(const QString &name,
                                const SireBase::PropertyMap &map) const
{
    auto all = this->cutGroups(name, map);

    if (all.count() > 1)
        _get_view<CutGroup>::raise_duplicate(name, all.count());

    BOOST_ASSERT( not all.isEmpty() );

    return all[0];
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
CutGroup SelectorM<T>::cutGroup(const CGID &cgid,
                                const SireBase::PropertyMap &map) const
{
    auto all = this->cutGroups(cgid, map);

    if (all.count() > 1)
        _get_view<CutGroup>::raise_duplicate(cgid.toString(), all.count());

    BOOST_ASSERT( not all.isEmpty() );

    return all[0];
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<Atom> SelectorM<T>::atoms() const
{
    return SelectorM<Atom>(*this);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<Atom> SelectorM<T>::atoms(int i, const SireBase::PropertyMap &map) const
{
    return SelectorM<Atom>(this->atom(i,map));
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<Atom> SelectorM<T>::atoms(const SireBase::Slice &slice,
                                    const SireBase::PropertyMap&) const
{
    return SelectorM<Atom>(*this, slice);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<Atom> SelectorM<T>::atoms(const QList<qint64> &idxs,
                                    const SireBase::PropertyMap&) const
{
    return SelectorM<Atom>(*this, idxs);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<Atom> SelectorM<T>::atoms(const QString &name,
                                    const SireBase::PropertyMap &map) const
{
    return SelectorM<Atom>(*this, name, map);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<Atom> SelectorM<T>::atoms(const AtomID &atomid,
                                    const SireBase::PropertyMap&) const
{
    return SelectorM<Atom>(*this, atomid);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<Residue> SelectorM<T>::residues() const
{
    return SelectorM<Residue>(*this);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<Residue> SelectorM<T>::residues(int i,
                                          const SireBase::PropertyMap &map) const
{
    return SelectorM<Residue>(this->residue(i, map));
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<Residue> SelectorM<T>::residues(const SireBase::Slice &slice,
                                          const SireBase::PropertyMap&) const
{
    return SelectorM<Residue>(*this, slice);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<Residue> SelectorM<T>::residues(const QList<qint64> &idxs,
                                          const SireBase::PropertyMap&) const
{
    return SelectorM<Residue>(*this, idxs);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<Residue> SelectorM<T>::residues(const QString &name,
                                          const SireBase::PropertyMap &map) const
{
    return SelectorM<Residue>(*this, name, map);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<Residue> SelectorM<T>::residues(const ResID &resid,
                                          const SireBase::PropertyMap&) const
{
    return SelectorM<Residue>(*this, resid);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<Chain> SelectorM<T>::chains() const
{
    return SelectorM<Chain>(*this);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<Chain> SelectorM<T>::chains(int i,
                                      const SireBase::PropertyMap &map) const
{
    return SelectorM<Chain>(this->chain(i, map));
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<Chain> SelectorM<T>::chains(const SireBase::Slice &slice,
                                      const SireBase::PropertyMap&) const
{
    return SelectorM<Chain>(*this, slice);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<Chain> SelectorM<T>::chains(const QList<qint64> &idxs,
                                      const SireBase::PropertyMap&) const
{
    return SelectorM<Chain>(*this, idxs);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<Chain> SelectorM<T>::chains(const QString &name,
                                      const SireBase::PropertyMap &map) const
{
    return SelectorM<Chain>(*this, name, map);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<Chain> SelectorM<T>::chains(const ChainID &chainid,
                                      const SireBase::PropertyMap&) const
{
    return SelectorM<Chain>(*this, chainid);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<Segment> SelectorM<T>::segments() const
{
    return SelectorM<Segment>(*this);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<Segment> SelectorM<T>::segments(int i,
                                          const SireBase::PropertyMap &map) const
{
    return SelectorM<Segment>(this->segment(i, map));
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<Segment> SelectorM<T>::segments(const SireBase::Slice &slice,
                                          const SireBase::PropertyMap&) const
{
    return SelectorM<Segment>(*this, slice);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<Segment> SelectorM<T>::segments(const QList<qint64> &idxs,
                                          const SireBase::PropertyMap&) const
{
    return SelectorM<Segment>(*this, idxs);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<Segment> SelectorM<T>::segments(const QString &name,
                                          const SireBase::PropertyMap &map) const
{
    return SelectorM<Segment>(*this, name, map);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<Segment> SelectorM<T>::segments(const SegID &segid,
                                          const SireBase::PropertyMap&) const
{
    return SelectorM<Segment>(*this, segid);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<CutGroup> SelectorM<T>::cutGroups() const
{
    return SelectorM<CutGroup>(*this);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<CutGroup> SelectorM<T>::cutGroups(int i,
                                            const SireBase::PropertyMap &map) const
{
    return SelectorM<CutGroup>(this->cutGroup(i, map));
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<CutGroup> SelectorM<T>::cutGroups(const SireBase::Slice &slice,
                                            const SireBase::PropertyMap&) const
{
    return SelectorM<CutGroup>(*this, slice);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<CutGroup> SelectorM<T>::cutGroups(const QList<qint64> &idxs,
                                            const SireBase::PropertyMap&) const
{
    return SelectorM<CutGroup>(*this, idxs);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<CutGroup> SelectorM<T>::cutGroups(const QString &name,
                                            const SireBase::PropertyMap &map) const
{
    return SelectorM<CutGroup>(*this, name, map);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorM<CutGroup> SelectorM<T>::cutGroups(const CGID &cgid,
                                            const SireBase::PropertyMap&) const
{
    return SelectorM<CutGroup>(*this, cgid);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
QList<typename T::Index> SelectorM<T>::IDs() const
{
    QList<typename T::Index> ids;

    for (const auto &v : this->vws)
    {
        ids += v.IDs();
    }

    return ids;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
QList<typename T::Index> SelectorM<T>::indexes() const
{
    QList<typename T::Index> idxs;

    for (const auto &v : this->vws)
    {
        idxs += v.indexes();
    }

    return idxs;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
QList<typename T::Number> SelectorM<T>::numbers() const
{
    QList<typename T::Number> nums;

    for (const auto &v : this->vws)
    {
        nums += v.numbers();
    }

    return nums;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
QList<typename T::Name> SelectorM<T>::names() const
{
    QList<typename T::Name> nmes;

    for (const auto &v : this->vws)
    {
        nmes += v.names();
    }

    return nmes;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
int SelectorM<T>::nAtoms() const
{
    int n = 0;

    for (const auto &v : this->vws)
    {
        n += v.nAtoms();
    }

    return n;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
int SelectorM<T>::nResidues() const
{
    int n = 0;

    for (const auto &v : this->vws)
    {
        n += v.nResidues();
    }

    return n;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
int SelectorM<T>::nChains() const
{
    int n = 0;

    for (const auto &v : this->vws)
    {
        n += v.nChains();
    }

    return n;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
int SelectorM<T>::nSegments() const
{
    int n = 0;

    for (const auto &v : this->vws)
    {
        n += v.nSegments();
    }

    return n;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
int SelectorM<T>::nCutGroups() const
{
    int n = 0;

    for (const auto &v : this->vws)
    {
        n += v.nCutGroups();
    }

    return n;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
int SelectorM<T>::nMolecules() const
{
    return this->vws.count();
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
int SelectorM<T>::nFrames() const
{
    return this->nFrames(PropertyMap());
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
int SelectorM<T>::nFrames(const SireBase::PropertyMap &map) const
{
    return SireMol::detail::_nFrames(this->vws, map);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
void SelectorM<T>::loadFrame(int frame)
{
    this->loadFrame(frame, PropertyMap());
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
void SelectorM<T>::saveFrame(int frame)
{
    this->saveFrame(frame, PropertyMap());
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
void SelectorM<T>::saveFrame()
{
    this->saveFrame(PropertyMap());
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
void SelectorM<T>::deleteFrame(int frame)
{
    this->deleteFrame(frame, PropertyMap());
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
void SelectorM<T>::loadFrame(int frame, const SireBase::PropertyMap &map)
{
    SireMol::detail::_loadFrame(this->vws, frame, map);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
void SelectorM<T>::saveFrame(int frame, const SireBase::PropertyMap &map)
{
    SireMol::detail::_saveFrame(this->vws, frame, map);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
void SelectorM<T>::saveFrame(const SireBase::PropertyMap &map)
{
    SireMol::detail::_saveFrame(this->vws, map);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
void SelectorM<T>::deleteFrame(int frame, const SireBase::PropertyMap &map)
{
    SireMol::detail::_deleteFrame(this->vws, frame, map);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool SelectorM<T>::isEmpty() const
{
    return this->vws.isEmpty();
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
typename SelectorM<T>::const_iterator SelectorM<T>::begin() const
{
    return this->vws.constBegin();
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
typename SelectorM<T>::const_iterator SelectorM<T>::end() const
{
    return this->vws.constEnd();
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
typename SelectorM<T>::const_iterator SelectorM<T>::constBegin() const
{
    return this->vws.constBegin();
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
typename SelectorM<T>::const_iterator SelectorM<T>::constEnd() const
{
    return this->vws.constEnd();
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
QString SelectorM<T>::toString() const
{
    if (this->isEmpty())
    {
        return QObject::tr("%1::empty").arg(this->what());
    }
    else
    {
        QStringList parts;

        const auto n = this->count();

        if (n <= 10)
        {
            for (int i=0; i<n; ++i)
            {
                const auto view = this->operator[](i);

                parts.append(QString("%1: %2 %3")
                    .arg(i).arg(view.data().number().toString())
                    .arg(view.toString()));
            }
        }
        else
        {
            for (int i=0; i<5; ++i)
            {
                const auto view = this->operator[](i);

                parts.append(QString("%1: %2 %3")
                    .arg(i).arg(view.data().number().toString())
                    .arg(view.toString()));
            }

            parts.append("...");

            for (int i=n-5; i<n; ++i)
            {
                const auto view = this->operator[](i);

                parts.append(QString("%1: %2 %3")
                    .arg(i).arg(view.data().number().toString())
                    .arg(view.toString()));
            }
        }

        return QObject::tr("%1( size=%2\n%3\n)")
                    .arg(this->what()).arg(n)
                    .arg(parts.join("\n"));
    }
}

#endif // SIRE_SKIP_INLINE_FUNCTIONS

} // end of namespace SireMol

#include "evaluatorm.h"

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

namespace SireMol
{

template<class T>
SIRE_OUTOFLINE_TEMPLATE
EvaluatorM::EvaluatorM(const SelectorM<T> &views)
           : SireBase::ConcreteProperty<EvaluatorM,SireBase::Property>()
{
    this->vws.reserve(views.count());

    for (const auto &view : views)
    {
        this->vws.append(PartialMolecule(view));
    }
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
EvaluatorM SelectorM<T>::evaluate() const
{
    return EvaluatorM(*this);
}

} // end of namespace SireMol

#endif // SIRE_SKIP_INLINE_FUNCTIONS

SIRE_END_HEADER

#endif
