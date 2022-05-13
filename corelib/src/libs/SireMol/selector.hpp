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

#ifndef SIREMOL_SELECTOR_HPP
#define SIREMOL_SELECTOR_HPP

#include "moleculeview.h"
#include "evaluator.h"

#include "atom.h"
#include "cutgroup.h"
#include "residue.h"
#include "chain.h"
#include "segment.h"

#include "SireBase/slice.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"
#include "SireStream/errors.h"

#include "SireMol/errors.h"

#include "tostring.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
template<class T>
class Selector;
}

template<class T>
SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::Selector<T>&);

template<class T>
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::Selector<T>&);

namespace SireMol
{

class MoleculeData;

namespace detail
{

SIREMOL_EXPORT bool has_property(const Atom*, const MoleculeData &moldata,
                  const PropertyName &key);
SIREMOL_EXPORT bool has_property(const Chain*, const MoleculeData &moldata,
                  const PropertyName &key);
SIREMOL_EXPORT bool has_property(const CutGroup*, const MoleculeData &moldata,
                  const PropertyName &key);
SIREMOL_EXPORT bool has_property(const Residue*, const MoleculeData &moldata,
                  const PropertyName &key);
SIREMOL_EXPORT bool has_property(const Segment*, const MoleculeData &moldata,
                  const PropertyName &key);

SIREMOL_EXPORT bool has_metadata(const Atom*, const MoleculeData &moldata,
                  const PropertyName &metakey);
SIREMOL_EXPORT bool has_metadata(const Chain*, const MoleculeData &moldata,
                  const PropertyName &metakey);
SIREMOL_EXPORT bool has_metadata(const CutGroup*, const MoleculeData &moldata,
                  const PropertyName &metakey);
SIREMOL_EXPORT bool has_metadata(const Residue*, const MoleculeData &moldata,
                  const PropertyName &metakey);
SIREMOL_EXPORT bool has_metadata(const Segment*, const MoleculeData &moldata,
                  const PropertyName &metakey);

SIREMOL_EXPORT bool has_metadata(const Atom*, const MoleculeData &moldata,
                  const PropertyName &key, const PropertyName &metakey);
SIREMOL_EXPORT bool has_metadata(const Chain*, const MoleculeData &moldata,
                  const PropertyName &key, const PropertyName &metakey);
SIREMOL_EXPORT bool has_metadata(const CutGroup*, const MoleculeData &moldata,
                  const PropertyName &key, const PropertyName &metakey);
SIREMOL_EXPORT bool has_metadata(const Residue*, const MoleculeData &moldata,
                  const PropertyName &key, const PropertyName &metakey);
SIREMOL_EXPORT bool has_metadata(const Segment*, const MoleculeData &moldata,
                  const PropertyName &key, const PropertyName &metakey);

} // end of namespace detail

/** This template class provides a way to manipulate the selection
    of parts of a a molecule, e.g. Selector<Atom> provides a
    MoleculeView that has a user-configurable view of an
    arbitrary set of atoms. This class is used by the code
    whenever multiple groups from within the same molecule
    are returned, e.g. when returning all residues that
    have a particular name

    @author Christopher Woods
*/
template<class T>
class SIREMOL_EXPORT Selector : public SireBase::ConcreteProperty<Selector<T>,MoleculeView>
{

friend SIREMOL_EXPORT QDataStream& ::operator<<<>(QDataStream&, const Selector<T>&);
friend SIREMOL_EXPORT QDataStream& ::operator>><>(QDataStream&, Selector<T>&);

public:
    Selector();

    Selector(const T &view);

    Selector(const MoleculeData &moldata);

    Selector(const MoleculeData &moldata,
             const AtomSelection &selected_atoms);

    Selector(const MoleculeData &moldata, const typename T::ID &viewid);

    Selector(const MoleculeData &moldata, const QList<typename T::Index> &idxs);

    Selector(const Selector<T> &other);

    ~Selector();

    Selector<T>& operator=(const Selector<T> &other);
    Selector<T>& operator=(const T &view);

    static const char* typeName();

    Selector<T>* clone() const;

    QList<typename T::Index> IDs() const;
    QList<typename T::Index> indexes() const;
    QList<typename T::Name> names() const;
    QList<typename T::Number> numbers() const;

    bool operator==(const Selector<T> &other) const;
    bool operator!=(const Selector<T> &other) const;

    Selector<T> operator+(const Selector<T> &other) const;
    Selector<T> operator-(const Selector<T> &other) const;

    Selector<T> operator+(const typename T::ID &id) const;
    Selector<T> operator-(const typename T::ID &id) const;

    Selector<T> operator+(const T &view) const;
    Selector<T> operator-(const T &view) const;

    int nViews() const;

    MolViewPtr operator[](int i) const;
    MolViewPtr operator[](const QString &name) const;
    MolViewPtr operator[](const SireBase::Slice &slice) const;
    MolViewPtr operator[](const QList<qint64> &idxs) const;
    MolViewPtr operator[](const AtomID &atomid) const;
    MolViewPtr operator[](const ResID &resid) const;
    MolViewPtr operator[](const CGID &cgid) const;
    MolViewPtr operator[](const ChainID &chainid) const;
    MolViewPtr operator[](const SegID &segid) const;
    MolViewPtr operator[](const SireID::Index &idx) const;

    T operator()(int i) const;
    T operator()(const SireID::Index &idx) const;
    Selector<T> operator()(int i, int j) const;
    Selector<T> operator()(const SireBase::Slice &slice) const;
    Selector<T> operator()(const QList<qint64> &idxs) const;
    Selector<Atom> operator()(const AtomID &atomid) const;
    Selector<Residue> operator()(const ResID &resid) const;
    Selector<CutGroup> operator()(const CGID &cgid) const;
    Selector<Chain> operator()(const ChainID &chainid) const;
    Selector<Segment> operator()(const SegID &segid) const;

    typename T::Index index(int i) const;

    QString toString() const;

    bool isEmpty() const;
    bool selectedAll() const;

    Selector<T> add(const Selector<T> &other) const;
    Selector<T> add(const T &view) const;
    Selector<T> add(const typename T::ID &id) const;

    Selector<T> subtract(const Selector<T> &other) const;
    Selector<T> subtract(const T &view) const;
    Selector<T> subtract(const typename T::ID &id) const;

    Selector<T> intersection(const Selector<T> &other) const;
    Selector<T> intersection(const T &view) const;
    Selector<T> intersection(const typename T::ID &id) const;

    Selector<T> invert() const;

    bool intersects(const Selector<T> &other) const;
    bool intersects(const T &view) const;
    bool intersects(const typename T::ID &id) const;

    bool contains(const Selector<T> &other) const;
    bool contains(const T &view) const;
    bool contains(const typename T::ID &id) const;

    AtomSelection selection() const;
    AtomSelection selection(int i) const;
    AtomSelection selection(int i, int j) const;

    Mover< Selector<T> > move() const;
    Mover< Selector<T> > move(int i) const;
    Mover< Selector<T> > move(int i, int j) const;

    Evaluator evaluate() const;
    Evaluator evaluate(int i) const;
    Evaluator evaluate(int i, int j) const;

    Selector<T> selector() const;
    Selector<T> selector(int i) const;
    Selector<T> selector(int i, int j) const;

    bool hasProperty(const PropertyName &key) const;
    bool hasMetadata(const PropertyName &metakey) const;
    bool hasMetadata(const PropertyName &key,
                     const PropertyName &metakey) const;

    QStringList propertyKeys() const;
    QStringList metadataKeys() const;
    QStringList metadataKeys(const PropertyName &key) const;

    template<class V>
    QList<V> property(const PropertyName &key) const;

    template<class V>
    QList<V> metadata(const PropertyName &metakey) const;

    template<class V>
    QList<V> metadata(const PropertyName &key,
                      const PropertyName &metakey) const;

protected:
    template<class V>
    void setProperty(const QString &key, const QList<V> &values);

    template<class V>
    void setMetadata(const QString &metakey, const QList<V> &values);

    template<class V>
    void setMetadata(const QString &key, const QString &metakey,
                     const QList<V> &values);

    template<class V>
    void setProperty(const QString &key, const V &value);

    template<class V>
    void setMetadata(const QString &metakey, V &value);

    template<class V>
    void setMetadata(const QString &key, const QString &metakey,
                     const V &value);

private:
    /** The list of indicies of the selected parts
        of the molecule. This is empty if the entire
        molecule is selected
    */
    QList<typename T::Index> idxs;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Null constructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Selector<T>::Selector() : SireBase::ConcreteProperty<Selector<T>,MoleculeView>()
{}

/** Construct the set of all groups for the molecule whose
    data is in 'moldata' */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Selector<T>::Selector(const MoleculeData &moldata)
            : SireBase::ConcreteProperty<Selector<T>,MoleculeView>(moldata)
{}

/** Construct the set of all groups for the molecule whose data is in 'moldata'
    and with the specfied groups selected */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Selector<T>::Selector(const MoleculeData &moldata,
                      const QList<typename T::Index> &indexes)
            : SireBase::ConcreteProperty<Selector<T>,MoleculeView>(moldata)
{
    int n = detail::getCount<T>(moldata.info());

    //ensure all of the indexes are valid!
    foreach (typename T::Index idx, indexes)
    {
        idx.map(n);
    }

    idxs = indexes;
}

/** Construct the set of all groups that contain at least
    one of the atoms selected in 'selected_atoms'

    \throw SireError::incompatible_error
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Selector<T>::Selector(const MoleculeData &moldata,
                      const AtomSelection &selected_atoms)
            : SireBase::ConcreteProperty<Selector<T>,MoleculeView>(moldata)
{
    if (selected_atoms.selectedNone())
    {
        this->operator=(Selector<T>());
    }
    else if (not selected_atoms.selectedAll())
    {
        idxs = detail::getAll<T>(moldata.info(), selected_atoms);
    }
}

/** Construct the set that contains only the view 'view' */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Selector<T>::Selector(const T &view)
            : SireBase::ConcreteProperty<Selector<T>,MoleculeView>(view)
{
    idxs.append(view.index());
}

/** Construct the set of parts that match the ID 'id' from the
    molecule whose data is in 'moldata'

    \throw SireMol::missing_atom
    \throw SireMol::missing_residue
    \throw SireMol::missing_chain
    \throw SireMol::missing_segment
    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Selector<T>::Selector(const MoleculeData &moldata,
                      const typename T::ID &viewid)
            : SireBase::ConcreteProperty<Selector<T>,MoleculeView>(moldata)
{
    idxs = viewid.map(moldata.info());

    if (idxs.count() == detail::getCount<T>(moldata.info()))
    {
        idxs.clear();
    }
    else if (idxs.isEmpty())
    {
        this->operator=(Selector<T>());
    }
}

/** Copy constructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Selector<T>::Selector(const Selector<T> &other)
            : SireBase::ConcreteProperty<Selector<T>,MoleculeView>(other),
              idxs(other.idxs)
{}

/** Destructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Selector<T>::~Selector()
{}

/** Copy assignment operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Selector<T>& Selector<T>::operator=(const Selector<T> &other)
{
    MoleculeView::operator=(other);
    idxs = other.idxs;

    return *this;
}

/** Copy so that this contains just 'view' */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Selector<T>& Selector<T>::operator=(const T &view)
{
    MoleculeView::operator=(view);

    idxs.clear();
    idxs.append(view.index());

    return *this;
}

/** Comparison operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool Selector<T>::operator==(const Selector<T> &other) const
{
    return idxs == other.idxs and MoleculeView::operator==(other);
}

/** Comparison operator */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool Selector<T>::operator!=(const Selector<T> &other) const
{
    return idxs != other.idxs or MoleculeView::operator!=(other);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
const char* Selector<T>::typeName()
{
    return QMetaType::typeName( qMetaTypeId< Selector<T> >() );
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
Selector<T>* Selector<T>::clone() const
{
    return new Selector<T>(*this);
}

/** Return the IDs of all of the items in this selector, in the
 *  order they appear in the selector
 */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QList<typename T::Index> Selector<T>::IDs() const
{
    if (this->isEmpty())
        return QList<typename T::Index>();

    else if (this->selectedAll())
    {
        const int n = detail::getCount<T>(this->data().info());

        QList<typename T::Index> ret;
        ret.reserve(n);

        for (int i=0; i<n; ++i)
        {
            ret.append(typename T::Index(i));
        }

        return ret;
    }
    else
        return idxs;
}

/** Return the indexes of all of the items in this selector, in the
 *  order they appear in the selector
 */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QList<typename T::Index> Selector<T>::indexes() const
{
    return this->IDs();
}

/** Return the names of all of the items in this selector, in the
 *  order they appear in the selector
 */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QList<typename T::Name> Selector<T>::names() const
{
    auto indexes = this->IDs();

    QList<typename T::Name> n;
    n.reserve(indexes.count());

    const auto &molinfo = this->data().info();

    for (const auto &index : indexes)
    {
        n.append( molinfo.name(index) );
    }

    return n;
}

/** Return the numbers of all of the items in this selector, in the
 *  order they appear in the selector
 */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QList<typename T::Number> Selector<T>::numbers() const
{
    auto indexes = this->IDs();

    QList<typename T::Number> n;
    n.reserve(indexes.count());

    const auto &molinfo = this->data().info();

    for (const auto &index : indexes)
    {
        n.append( molinfo.number(index) );
    }

    return n;
}

/** Return whether this set is empty */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool Selector<T>::isEmpty() const
{
    return this->data().isEmpty();
}

/** Return whether or not the entire molecule is selected */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool Selector<T>::selectedAll() const
{
    // have no indexes if we have selected all of the molecule
    return (not this->isEmpty()) and idxs.isEmpty();
}

/** Return a string representation of this selector */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QString Selector<T>::toString() const
{
    if (this->isEmpty())
    {
        return QObject::tr("Selector<%1>::empty").arg(T::typeName());
    }
    else
    {
        QStringList parts;

        if (this->size() <= 10)
        {
            for (int i=0; i<this->size(); ++i)
            {
                parts.append(QString("%1:  %2").arg(i).arg(this->operator()(i).toString()));
            }
        }
        else
        {
            for (int i=0; i<5; ++i)
            {
                parts.append(QString("%1:  %2").arg(i).arg(this->operator()(i).toString()));
            }

            parts.append("...");

            for (int i=this->size()-5; i<this->size(); ++i)
            {
                parts.append(QString("%1:  %3").arg(i).arg(this->operator()(i).toString()));
            }
        }

        return QObject::tr( "Selector<%1>( size=%2\n%3\n)")
                        .arg(T::typeName()).arg(this->count())
                        .arg(parts.join("\n"));
    }
}

/** Return the sum of this selection with 'other' - items from
    'other' are appended to this set only if they are not already
    in this set

    \throw SireError::incompatible_error
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Selector<T> Selector<T>::add(const Selector<T> &other) const
{
    if (this->selectedAll())
    {
        if (other.isEmpty())
            return *this;

        MoleculeView::assertSameMolecule(other);

        return *this;
    }
    else if (other.selectedAll())
    {
        if (this->isEmpty())
            return other;

        MoleculeView::assertSameMolecule(other);

        // do this so that we always use the version from the left hand molecule
        Selector<T> ret(*this);
        ret.idxs.clear();
        return ret;
    }

    MoleculeView::assertSameMolecule(other);

    Selector<T> ret(*this);

    QSet<typename T::Index> seen_idxs = _list_to_set(this->idxs);

    const int n = detail::getCount<T>(this->data().info());

    foreach (typename T::Index idx, other.idxs)
    {
        if (not seen_idxs.contains(idx))
        {
            ret.idxs.append(idx);
            if (ret.idxs.count() == n)
            {
                ret.idxs.clear();
                break;
            }
        }
    }

    return ret;
}

/** Return the view 'view' to this set */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Selector<T> Selector<T>::add(const T &view) const
{
    return this->add(Selector<T>(view));
}

/** Return the set that has all views that match the ID 'id'
    added

    \throw SireMol::missing_atom
    \throw SireMol::missing_residue
    \throw SireMol::missing_chain
    \throw SireMol::missing_segment
    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Selector<T> Selector<T>::add(const typename T::ID &id) const
{
    return this->operator+(Selector<T>(this->data(), id));
}

/** Subtract all of the views in 'other' from this set

    \throw SireError::incompatible_error
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Selector<T> Selector<T>::subtract(const Selector<T> &other) const
{
    if (this->isEmpty())
        return Selector<T>();
    else if (other.selectedAll())
    {
        MoleculeView::assertSameMolecule(other);
        return Selector<T>();
    }

    Selector<T> ret(*this);

    const int n = detail::getCount<T>(this->data().info());

    if (ret.idxs.isEmpty())
    {

        for (int i=0; i<n; ++i)
        {
            ret.idxs.append(typename T::Index(i));
        }
    }

    foreach (typename T::Index idx, other.idxs)
    {
        ret.idxs.removeAll(idx);

        if (ret.idxs.isEmpty())
            return Selector<T>();
    }

    if (ret.idxs.count() == n)
    {
        // nothing was removed?
        ret.idxs.clear();
    }

    return ret;
}

/** Remove the view 'view' from this set

    \throw SireError::incompatible_error
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Selector<T> Selector<T>::subtract(const T &view) const
{
    return this->subtract(Selector<T>(view));
}

/** Subtract all of the views that match the ID 'id' from
    this set

    \throw SireMol::missing_atom
    \throw SireMol::missing_residue
    \throw SireMol::missing_chain
    \throw SireMol::missing_segment
    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Selector<T> Selector<T>::subtract(const typename T::ID &id) const
{
    return this->subtract(Selector<T>(this->data(), id));
}

/** Syntactic sugar for add() */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Selector<T> Selector<T>::operator+(const Selector<T> &other) const
{
    return this->add(other);
}

/** Syntactic sugar for subtract() */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Selector<T> Selector<T>::operator-(const Selector<T> &other) const
{
    return this->subtract(other);
}

/** Syntactic sugar for add() */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Selector<T> Selector<T>::operator+(const typename T::ID &id) const
{
    return this->add(id);
}

/** Syntactic sugar for subtract() */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Selector<T> Selector<T>::operator-(const typename T::ID &id) const
{
    return this->subtract(id);
}

/** Syntactic sugar for add() */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Selector<T> Selector<T>::operator+(const T &view) const
{
    return this->add(view);
}

/** Syntactic sugar for subtract() */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Selector<T> Selector<T>::operator-(const T &view) const
{
    return this->subtract(view);
}

/** Return the ith view in this set (this supports negative indexing!)

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
MolViewPtr Selector<T>::operator[](int i) const
{
    return this->operator()(i);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
MolViewPtr Selector<T>::operator[](const QString &name) const
{
    auto match = this->operator()(AtomName(name));

    if (match.nAtoms() == 1)
        return match(0);
    else
        return match;
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
MolViewPtr Selector<T>::operator[](const SireBase::Slice &slice) const
{
    return this->operator()(slice);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
MolViewPtr Selector<T>::operator[](const QList<qint64> &idxs) const
{
    return this->operator()(idxs);
}

/** Exposing MoleculeView::operator[] */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
MolViewPtr Selector<T>::operator[](const AtomID &atom) const
{
    auto match = this->operator()(atom);

    if (match.nAtoms() == 1)
        return match(0);
    else
        return match;
}

/** Exposing MoleculeView::operator[] */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
MolViewPtr Selector<T>::operator[](const ResID &res) const
{
    auto match = this->operator()(res);

    if (match.nResidues() == 1)
        return match(0);
    else
        return match;
}

/** Exposing MoleculeView::operator[] */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
MolViewPtr Selector<T>::operator[](const CGID &cg) const
{
    auto match = this->operator()(cg);

    if (match.nCutGroups() == 1)
        return match(0);
    else
        return match;
}

/** Exposing MoleculeView::operator[] */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
MolViewPtr Selector<T>::operator[](const ChainID &chain) const
{
    auto match = this->operator()(chain);

    if (match.nChains() == 1)
        return match(0);
    else
        return match;
}

/** Exposing MoleculeView::operator[] */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
MolViewPtr Selector<T>::operator[](const SegID &seg) const
{
    auto match = this->operator()(seg);

    if (match.nSegments() == 1)
        return match(0);
    else
        return match;
}

/** Exposing MoleculeView::operator[] */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
MolViewPtr Selector<T>::operator[](const SireID::Index &idx) const
{
    return this->operator()(idx.value());
}

/** Return the index of the ith view in this set (this supports negative indexing)

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
typename T::Index Selector<T>::index(int i) const
{
    if (not this->isEmpty())
    {
        if (idxs.isEmpty())
        {
            const int n = detail::getCount<T>(this->data().info());
            return typename T::Index(SireID::Index(i).map(n));
        }
        else
        {
            return idxs.at( SireID::Index(i).map(idxs.count()));
        }
    }

    return typename T::Index(SireID::Index(i).map(0));
}

/** Return the ith view in this set (this supports negative indexing!)

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
T Selector<T>::operator()(int i) const
{
    return T(this->data(), this->index(i));
}

/** Return the range of views from index i to j in this set. This
    supports negative indexing, and also, if j is less then i, then
    the order of views in the returned set is reversed

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Selector<T> Selector<T>::operator()(int i, int j) const
{
    if (this->isEmpty())
        // raise an index exception
        i = SireID::Index(i).map(0);

    if (this->selectedAll())
    {
        const int n = detail::getCount<T>(this->data().info());

        i = SireID::Index(i).map(n);
        j = SireID::Index(j).map(n);

        Selector<T> ret(*this);
        ret.idxs.clear();

        if (i <= j)
        {
            if (i == 0 and j == n-1)
                // selected all
                return ret;

            for ( ; i<=j; ++i)
            {
                ret.idxs.append(typename T::Index(i));
            }
        }
        else
        {
            if (i == n-1 and j == 0)
                // selected all
                return ret;

            for ( ; i >= j; --i)
            {
                ret.idxs.append(typename T::Index(i));
            }
        }

        if (ret.idxs.count() == n)
            ret.idxs.clear();

        return ret;
    }
    else
    {
        i = SireID::Index(i).map(idxs.count());
        j = SireID::Index(j).map(idxs.count());

        Selector<T> ret(*this);
        ret.idxs.clear();

        if (i <= j)
        {
            for ( ; i<=j; ++i)
            {
                ret.idxs.append(this->idxs.at(i));
            }
        }
        else
        {
            for ( ; i >= j; --i)
            {
                ret.idxs.append(this->idxs.at(i));
            }
        }

        return ret;
    }
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
Selector<T> Selector<T>::operator()(const SireBase::Slice &slice) const
{
    if (this->isEmpty())
        //raise an exception
        slice.begin(0);

    if (this->selectedAll())
    {
        const int n = detail::getCount<T>(this->data().info());

        Selector<T> ret(*this);
        ret.idxs.clear();

        for (auto it = slice.begin(n); not it.atEnd(); it.next())
        {
            ret.idxs.append(typename T::Index(it.value()));
        }

        if (ret.idxs.count() == n)
            ret.idxs.clear();

        return ret;
    }
    else
    {
        Selector<T> ret(*this);
        ret.idxs.clear();

        for (auto it = slice.begin(this->idxs.count());
             not it.atEnd(); it.next())
        {
            ret.idxs.append(this->idxs.at(it.value()));
        }

        return ret;
    }
}

/** Return the range of views from whose indicies are in idxs in this set.

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Selector<T> Selector<T>::operator()(const QList<qint64> &idxs) const
{
    if (idxs.isEmpty())
        return Selector<T>();

    else if (this->isEmpty())
    {
        //raise an exception
        SireID::Index(idxs.at(0)).map(0);
        return Selector<T>();
    }
    else if (this->selectedAll())
    {
        Selector<T> ret(*this);
        ret.idxs.clear();

        const int n = detail::getCount<T>(this->data().info());

        QSet<typename T::Index> seen;

        for (const auto &idx : idxs)
        {
            typename T::Index index(SireID::Index(idx).map(n));

            if (not seen.contains(index))
            {
                ret.idxs.append(index);
                seen.insert(index);
            }
        }

        if (ret.idxs.count() >= n)
            ret.idxs.clear();

        return ret;
    }
    else
    {
        Selector<T> ret(*this);
        ret.idxs.clear();

        QSet<typename T::Index> seen;

        for (const auto &idx : idxs)
        {
            auto index = this->idxs.at(SireID::Index(idx).map(this->idxs.count()));

            if (not seen.contains(index))
            {
                ret.idxs.append(index);
                seen.insert(index);
            }
        }

        return ret;
    }
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
Selector<Atom> Selector<T>::operator()(const AtomID &atomid) const
{
    return this->atoms(atomid);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
Selector<Residue> Selector<T>::operator()(const ResID &resid) const
{
    return this->residues(resid);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
Selector<CutGroup> Selector<T>::operator()(const CGID &cgid) const
{
    return this->cutGroups(cgid);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
Selector<Chain> Selector<T>::operator()(const ChainID &chainid) const
{
    return this->chains(chainid);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
Selector<Segment> Selector<T>::operator()(const SegID &segid) const
{
    return this->segments(segid);
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
T Selector<T>::operator()(const SireID::Index &idx) const
{
    return this->operator()(idx.value());
}

/** Return the number of views in this set */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
int Selector<T>::nViews() const
{
    if (this->isEmpty())
        return 0;
    else if (this->selectedAll())
        return detail::getCount<T>(this->data().info());
    else
        return this->idxs.count();
}

template<class T>
SIRE_INLINE_TEMPLATE
QSet<T> _list_to_set(const QList<T> &vals)
{
    #if QT_VERSION >= QT_VERSION_CHECK(5, 1, 4)
        return QSet<T>(vals.constBegin(), vals.constEnd());
    #else
        return vals.toSet();
    #endif
}

/** Return the intersection of this set with 'other' - the
    views in both sets are returned, in the order that they
    appear in this set

    \throw SireError::incompatible_error
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Selector<T> Selector<T>::intersection(const Selector<T> &other) const
{
    if (this->isEmpty() or other.isEmpty())
        return Selector<T>();

    else if (this->selectedAll())
    {
        MoleculeView::assertSameMolecule(other);
        Selector<T> ret(*this);
        ret.idxs = other.idxs;
        return ret;
    }
    else if (other.selectedAll())
    {
        MoleculeView::assertSameMolecule(other);
        return *this;
    }
    else
    {
        MoleculeView::assertSameMolecule(other);

        Selector<T> ret(*this);

        auto seen = _list_to_set(other.idxs);

        for (const auto &idx : this->idxs)
        {
            if (seen.contains(idx))
                ret.idxs.append(idx);
        }

        return ret;
    }
}

/** Return the intersection of this set with 'view'

    \throw SireError::incompatible_error
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Selector<T> Selector<T>::intersection(const T &view) const
{
    return this->intersection(Selector<T>(view));
}

/** Return the intersection of this set with the views identified
    by the ID 'id'

    \throw SireMol::missing_atom
    \throw SireMol::missing_residue
    \throw SireMol::missing_chain
    \throw SireMol::missing_segment
    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Selector<T> Selector<T>::intersection(const typename T::ID &id) const
{
    return this->intersection(Selector<T>(this->data(), id));
}

/** Return the set that has a completely inverted selection */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Selector<T> Selector<T>::invert() const
{
    if (this->isEmpty())
    {
        throw SireMol::missing_atom(QObject::tr(
            "Cannot invert an empty selection!"), CODELOC);
    }
    else if (this->selectedAll())
    {
        return Selector<T>();
    }
    else
    {
        Selector<T> ret(*this);
        ret.idxs.clear();

        const int n = detail::getCount<T>(this->data().info());

        if (this->idxs.count() == 1)
        {
            auto seen = this->idxs.at(0);

            for (int i=0; i<n; ++i)
            {
                typename T::Index idx(i);

                if (idx != seen)
                    ret.idxs.append(idx);
            }
        }
        else
        {
            auto seen = _list_to_set(this->idxs);

            for (int i=0; i<n; ++i)
            {
                typename T::Index idx(i);

                if (not seen.contains(idx))
                    ret.idxs.append(idx);
            }
        }

        return ret;
    }
}

/** Return whether this set contains all of the views
    in 'other' */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool Selector<T>::contains(const Selector<T> &other) const
{
    if (this->isEmpty() or other.isEmpty())
        return false;

    MoleculeView::assertSameMolecule(other);

    if (this->selectedAll())
        return true;
    else if (other.selectedAll())
        return false;

    else
    {
        auto seen = _list_to_set(this->idxs);

        for (const auto &idx : other.idxs)
        {
            if (not seen.contains(idx))
                return false;
        }

        return true;
    }
}

/** Return whether this set contains some of the view
    in other

    \throw SireError::incompatible_error
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool Selector<T>::intersects(const Selector<T> &other) const
{
    if (this->isEmpty() or other.isEmpty())
        return false;

    MoleculeView::assertSameMolecule(other);

    if (this->selectedAll() or other.selectedAll())
        return true;

    auto seen = _list_to_set(this->idxs);

    for (const auto &idx : other.idxs)
    {
        if (seen.contains(idx))
            return true;
    }

    return false;
}

/** Return whether or not this set contains the view 'view' */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool Selector<T>::contains(const T &view) const
{
    return this->contains(Selector<T>(view));
}

/** Return whether or not this set contains all of the
    view identified by the ID 'id' */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool Selector<T>::contains(const typename T::ID &id) const
{
    return this->contains(Selector<T>(this->data(), id));
}

/** Return whether or not this set contains the view 'view' */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool Selector<T>::intersects(const T &view) const
{
    return this->intersects(Selector<T>(view));
}

/** Return whether this set contains some of the views
    identified by the ID 'id' */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool Selector<T>::intersects(const typename T::ID &id) const
{
    return this->intersects(Selector<T>(this->data(), id));
}

/** Return all of the atoms selected in this set */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomSelection Selector<T>::selection() const
{
    if (this->isEmpty())
        return AtomSelection();

    else if (this->selectedAll())
    {
        return AtomSelection(this->data());
    }

    AtomSelection selected_atoms(this->data());

    selected_atoms.selectNone();

    for (const auto &idx : this->idxs)
    {
        selected_atoms.select(idx);
    }

    return selected_atoms;
}

/** Return the selection of the atoms in the ith view

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomSelection Selector<T>::selection(int i) const
{
    return this->operator()(i).selection();
}

/** Return the selection of the atoms in the ith to jth views

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomSelection Selector<T>::selection(int i, int j) const
{
    return this->operator()(i, j).selection();
}

/** Return an object that can move a copy of all of the views
    in this set */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Mover< Selector<T> > Selector<T>::move() const
{
    return Mover< Selector<T> >(*this);
}

/** Return an object that can move the ith view in this set */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Mover< Selector<T> > Selector<T>::move(int i) const
{
    return Mover< Selector<T> >(Selector<T>(this->operator()(i)));
}

/** Return an object that can move the ith to jth views in this set */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Mover< Selector<T> > Selector<T>::move(int i, int j) const
{
    return Mover< Selector<T> >(this->operator()(i,j));
}

/** Return an evaluator that can evaluate properties over all
    of the views in this set */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Evaluator Selector<T>::evaluate() const
{
    return Evaluator(*this);
}

/** Return an evaluator that can evaluate properties over the
    ith view in this set */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Evaluator Selector<T>::evaluate(int i) const
{
    return Evaluator(this->operator()(i));
}

/** Return an evaluator that can evaluate properties over
    the ith to jth views in this set */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Evaluator Selector<T>::evaluate(int i, int j) const
{
    return Evaluator(this->operator()(i,j));
}

/** Return a selector that can change the selection of this view */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Selector<T> Selector<T>::selector() const
{
    return *this;
}

/** Return a selector for the ith view that can change the
    selection of that view

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Selector<T> Selector<T>::selector(int i) const
{
    return Selector<T>(this->operator()(i));
}

/** Return a selector for the ith to jth views

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Selector<T> Selector<T>::selector(int i, int j) const
{
    return this->operator()(i, j);
}

/** Return a list of the values of the property called 'key' for each
    of the views in this set, in the order that they appear in this set.

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
template<class T>
template<class V>
SIRE_OUTOFLINE_TEMPLATE
QList<V> Selector<T>::property(const PropertyName &key) const
{
    T *ptr = 0;
    return detail::get_property<V>(ptr, this->data(), this->IDs(), key);
}

/** Return a list of all of the metadata called 'metakey' for each
    of the views in this set, in the order they appear in this set.

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
template<class T>
template<class V>
SIRE_OUTOFLINE_TEMPLATE
QList<V> Selector<T>::metadata(const PropertyName &metakey) const
{
    T *ptr = 0;
    return detail::get_metadata<V>(ptr, this->data(),
                                   this->IDs(), metakey);
}

/** Return a list of all of the metadata called 'key'/'metakey' for each
    of the views in this set, in the order they appear in this set.

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
template<class T>
template<class V>
SIRE_OUTOFLINE_TEMPLATE
QList<V> Selector<T>::metadata(const PropertyName &key,
                               const PropertyName &metakey) const
{
    T *ptr = 0;
    return detail::get_metadata<V>(ptr, this->data(), this->IDs(),
                                   key, metakey);
}

/** Set the property at key 'key' for all of the views in this set
    to the values in 'values', with the values given in the same
    order as the views in this set.

    \throw SireError::incompatible_error
    \throw SireError::invalid_cast
*/
template<class T>
template<class V>
SIRE_OUTOFLINE_TEMPLATE
void Selector<T>::setProperty(const QString &key, const QList<V> &values)
{
    T *ptr = 0;
    detail::set_property<V>(ptr, this->data(), this->IDs(), key, values);
}

/** Set the metadata at key 'metakey' for all of the views in this set
    to the values in 'values', with the values given in the same
    order as the views in this set.

    \throw SireError::incompatible_error
    \throw SireError::invalid_cast
*/
template<class T>
template<class V>
SIRE_OUTOFLINE_TEMPLATE
void Selector<T>::setMetadata(const QString &metakey, const QList<V> &values)
{
    T *ptr = 0;
    detail::set_metadata<V>(ptr, this->data(), this->IDs(), metakey, values);
}


/** Set the metadata at key 'key'/'metakey' for all of the views in this set
    to the values in 'values', with the values given in the same
    order as the views in this set.

    \throw SireError::incompatible_error
    \throw SireError::invalid_cast
*/
template<class T>
template<class V>
SIRE_OUTOFLINE_TEMPLATE
void Selector<T>::setMetadata(const QString &key, const QString &metakey,
                              const QList<V> &values)
{
    T *ptr = 0;
    detail::set_metadata<V>(ptr, this->data(), this->IDs(), key, metakey, values);
}

/** Set the property at key 'key' for all of the views in this set
    to the value 'value'

    \throw SireError::invalid_cast
*/
template<class T>
template<class V>
SIRE_OUTOFLINE_TEMPLATE
void Selector<T>::setProperty(const QString &key, const V &value)
{
    T *ptr = 0;
    detail::set_property<V>(ptr, this->data(), this->IDs(), key, value);
}

/** Set the metadata at metakey 'metakey' for all of the views in this set
    to the value 'value'

    \throw SireError::invalid_cast
*/
template<class T>
template<class V>
SIRE_OUTOFLINE_TEMPLATE
void Selector<T>::setMetadata(const QString &metakey, V &value)
{
    T *ptr = 0;
    detail::set_metadata<V>(ptr, this->data(), this->IDs(), metakey, value);
}

/** Set the metadata at key 'key'/'metakey' for all of the views in this set
    to the value 'value'

    \throw SireError::invalid_cast
*/
template<class T>
template<class V>
SIRE_OUTOFLINE_TEMPLATE
void Selector<T>::setMetadata(const QString &key, const QString &metakey,
                              const V &value)
{
    T *ptr = 0;
    detail::set_metadata<V>(ptr, this->data(), this->IDs(), key, metakey, value);
}

/** Return whether or not the views of this selector has a property
    at key 'key' */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool Selector<T>::hasProperty(const PropertyName &key) const
{
    T *ptr = 0;
    return detail::has_property(ptr, this->data(), key);
}

/** Return whether or not the views of this selector have metadata
    at metakey 'metakey' */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool Selector<T>::hasMetadata(const PropertyName &metakey) const
{
    T *ptr = 0;
    return detail::has_metadata(ptr, this->data(), metakey);
}

/** Return whether or not the property at key 'key' for the views
    of this selector has metadata at metakey 'metakey'

    \throw SireBase::missing_property
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool Selector<T>::hasMetadata(const PropertyName &key,
                              const PropertyName &metakey) const
{
    T *ptr = 0;
    return detail::has_metadata(ptr, this->data(), key, metakey);
}

/** Return all of the keys for the properties attached to the
    parts of this selection */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QStringList Selector<T>::propertyKeys() const
{
    if (this->isEmpty())
        return QStringList();
    else
        return this->operator()(0).propertyKeys();
}

/** Return all of the metakeys for the metadata attached to the
    parts of this selection */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QStringList Selector<T>::metadataKeys() const
{
    if (this->isEmpty())
        return QStringList();
    else
        return this->operator()(0).metadataKeys();
}

/** Return all of the metakeys for the metadata attached to the
    parts of this selection for the property at key 'key'

    \throw SireBase::missing_property
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QStringList Selector<T>::metadataKeys(const PropertyName &key) const
{
    if (this->isEmpty())
        return QStringList();
    else
        return this->operator()(0).metadataKeys(key);
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

} // end of namespace SireMol

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Extract from a binary datastream */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QDataStream &operator>>(QDataStream &ds, SireMol::Selector<T> &views)
{
    QString cls;
    SireStream::VersionID version;

    ds >> cls;

    if (cls != QString(SireMol::Selector<T>::typeName()))
    {
        throw SireStream::corrupted_data(QObject::tr(
            "Found the wrong class (%1) when trying to read a %2.")
                .arg(cls).arg(SireMol::Selector<T>::typeName()), CODELOC);
    }

    ds >> version;

    if (version == 2)
    {
        SireStream::SharedDataStream sds(ds);
        sds >> views.idxs >> static_cast<SireMol::MoleculeView&>(views);
    }
    else if (version == 1)
    {
        SireStream::SharedDataStream sds(ds);

        QList<typename T::Index> idxs;
        SireMol::MoleculeData moldata;
        sds >> moldata >> idxs;

        views = SireMol::Selector<T>(moldata, idxs);
    }
    else
    {
        throw SireStream::version_error(version, "1",
                                        SireMol::Selector<T>::typeName(),
                                        CODELOC);
    }

    return ds;
}

/** Serialise to a binary datastream */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QDataStream &operator<<(QDataStream &ds, const SireMol::Selector<T> &views)
{
    ds << QString(SireMol::Selector<T>::typeName()) << SireStream::VersionID(2);

    SireStream::SharedDataStream sds(ds);

    sds << views.idxs << static_cast<const SireMol::MoleculeView&>(views);

    return ds;
}

#endif // SIRE_SKIP_INLINE_FUNCTIONS

SIRE_END_HEADER

#endif
