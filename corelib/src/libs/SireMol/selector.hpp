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

#include "tostring.h"

SIRE_BEGIN_HEADER

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
    MolViewPtr operator[](const AtomID &atomid) const;
    MolViewPtr operator[](const ResID &resid) const;
    MolViewPtr operator[](const CGID &cgid) const;
    MolViewPtr operator[](const ChainID &chainid) const;
    MolViewPtr operator[](const SegID &segid) const;
    MolViewPtr operator[](const SireID::Index &idx) const;

    T operator()(int i) const;
    Selector<T> operator()(int i, int j) const;

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
    void _pvt_add(typename T::Index idx);
    void _pvt_sub(typename T::Index idx);

    /** The list of indicies of the selected parts
        of the molecule */
    QList<typename T::Index> idxs;

    /** The set of indicies of the selected parts - used
        for rapid searching */
    QSet<typename T::Index> idxs_set;
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
{
    idxs = detail::getAll<T>(moldata.info());
    idxs_set = convert_to_qset(idxs);
}

/** Construct the set of all groups for the molecule whose data is in 'moldata'
    and with the specfied groups selected */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Selector<T>::Selector(const MoleculeData &moldata,
                      const QList<typename T::Index> &indexes)
            : SireBase::ConcreteProperty<Selector<T>,MoleculeView>(moldata)
{
    int nats = moldata.info().nAtoms();

    //ensure all of the indexes are valid!
    foreach (typename T::Index idx, indexes)
    {
        idx.map(nats);
    }

    idxs = indexes;
    idxs_set = convert_to_qset(idxs);
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
    idxs = detail::getAll<T>(moldata.info(), selected_atoms);
    idxs_set = convert_to_qset(idxs);
}

/** Construct the set that contains only the view 'view' */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Selector<T>::Selector(const T &view)
            : SireBase::ConcreteProperty<Selector<T>,MoleculeView>(view)
{
    idxs.append(view.index());
    idxs_set.insert(view.index());
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
    idxs_set = convert_to_qset(idxs);
}

/** Copy constructor */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Selector<T>::Selector(const Selector<T> &other)
            : SireBase::ConcreteProperty<Selector<T>,MoleculeView>(other),
              idxs(other.idxs), idxs_set(other.idxs_set)
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
    idxs_set = other.idxs_set;

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

    idxs_set.clear();
    idxs_set.insert(view.index());

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
    return idxs;
}

/** Return whether this set is empty */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool Selector<T>::isEmpty() const
{
    return idxs.isEmpty();
}

/** Return whether or not the entire molecule is selected */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool Selector<T>::selectedAll() const
{
    //there has to be a better way of doing this...
    return this->selection().selectedAll();
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
                parts.append(QString("%1:  %2").arg(i).arg(this->at(i).read().toString()));
            }
        }
        else
        {
            for (int i=0; i<5; ++i)
            {
                parts.append(QString("%1:  %2").arg(i).arg(this->at(i).read().toString()));
            }

            parts.append("...");

            for (int i=this->size()-5; i<this->size(); ++i)
            {
                parts.append(QString("%1:  %3").arg(i).arg(this->at(i).read().toString()));
            }
        }

        return QObject::tr( "Selector<%1>( size=%2\n%3\n)")
                        .arg(T::typeName()).arg(this->count())
                        .arg(parts.join("\n"));
    }
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
void Selector<T>::_pvt_add(typename T::Index idx)
{
    if (not idxs_set.contains(idx))
    {
        idxs.append(idx);
        idxs_set.insert(idx);
    }
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
void Selector<T>::_pvt_sub(typename T::Index idx)
{
    if (idxs_set.contains(idx))
    {
        idxs.removeAll(idx);
        idxs_set.remove(idx);
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
    MoleculeView::assertSameMolecule(other);

    Selector<T> ret(*this);

    foreach (typename T::Index idx, other.idxs)
    {
        ret._pvt_add(idx);
    }

    return ret;
}

/** Return the view 'view' to this set */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Selector<T> Selector<T>::add(const T &view) const
{
    MoleculeView::assertSameMolecule(view);

    Selector<T> ret(*this);
    ret._pvt_add(view.index());
    return ret;
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
    Selector<T> ret(*this);
    QList<typename T::Index> idxs = id.map(this->d->info());

    foreach (typename T::Index idx, idxs)
    {
        ret._pvt_add(idx);
    }

    return ret;
}

/** Subtract all of the views in 'other' from this set

    \throw SireError::incompatible_error
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Selector<T> Selector<T>::subtract(const Selector<T> &other) const
{
    MoleculeView::assertSameMolecule(other);

    Selector<T> ret(*this);

    foreach (typename T::Index idx, other.idxs)
    {
        ret._pvt_sub(idx);
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
    MoleculeView::assertSameMolecule(view);

    Selector<T> ret(*this);
    ret._pvt_sub(view.index());

    return ret;
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
    QList<typename T::Index> idxs = id.map(this->d->info());

    Selector<T> ret(*this);

    foreach (typename T::Index idx, idxs)
    {
        ret._pvt_sub(idx);
    }

    return ret;
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
    return T( this->data(), idxs.at( SireID::Index(i).map(idxs.count()) ) );
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
MolViewPtr Selector<T>::operator[](const QString &name) const
{
    return this->operator[](AtomName(name));
}

template<class T>
SIRE_OUTOFLINE_TEMPLATE
MolViewPtr Selector<T>::operator[](const SireBase::Slice &slice) const
{
    QList<typename T::Index> new_idxs;

    for (auto it = slice.begin(idxs.count()); not it.atEnd(); it.next())
    {
        new_idxs.append(idxs.at(it.value()));
    }

    return Selector<T>( this->data(), new_idxs );
}

/** Exposing MoleculeView::operator[] */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
MolViewPtr Selector<T>::operator[](const AtomID &atom) const
{
    return MoleculeView::operator[](atom);
}

/** Exposing MoleculeView::operator[] */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
MolViewPtr Selector<T>::operator[](const ResID &res) const
{
    return MoleculeView::operator[](res);
}

/** Exposing MoleculeView::operator[] */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
MolViewPtr Selector<T>::operator[](const CGID &cg) const
{
    return MoleculeView::operator[](cg);
}

/** Exposing MoleculeView::operator[] */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
MolViewPtr Selector<T>::operator[](const ChainID &chain) const
{
    return MoleculeView::operator[](chain);
}

/** Exposing MoleculeView::operator[] */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
MolViewPtr Selector<T>::operator[](const SegID &seg) const
{
    return MoleculeView::operator[](seg);
}

/** Exposing MoleculeView::operator[] */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
MolViewPtr Selector<T>::operator[](const SireID::Index &idx) const
{
    return MoleculeView::operator[](idx);
}

/** Return the index of the ith view in this set (this supports negative indexing)

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
typename T::Index Selector<T>::index(int i) const
{
    return idxs.at( SireID::Index(i).map(idxs.count()) );
}

/** Return the ith view in this set (this supports negative indexing!)

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
T Selector<T>::operator()(int i) const
{
    return T( this->data(), idxs.at( SireID::Index(i).map(idxs.count()) ) );
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
    i = SireID::Index(i).map(idxs.count());
    j = SireID::Index(j).map(idxs.count());

    Selector<T> ret(*this);
    ret.idxs.clear();
    ret.idxs_set.clear();

    if (i <= j)
    {
        for ( ; i<=j; ++i)
        {
            ret._pvt_add( typename T::Index(i) );
        }
    }
    else
    {
        for ( ; i >= j; --i)
        {
            ret._pvt_add( typename T::Index(i) );
        }
    }

    return ret;
}

/** Return the number of views in this set */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
int Selector<T>::nViews() const
{
    return idxs.count();
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
    MoleculeView::assertSameMolecule(other);

    Selector<T> ret(*this);

    if (idxs.count() <= other.idxs.count())
    {
        foreach (typename T::Index idx, idxs)
        {
            if (not other.idxs_set.contains(idx))
                ret._pvt_sub(idx);
        }
    }
    else
    {
        foreach (typename T::Index idx, other.idxs)
        {
            if (not idxs_set.contains(idx))
                ret._pvt_sub(idx);
        }
    }

    return ret;
}

/** Return the intersection of this set with 'view'

    \throw SireError::incompatible_error
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Selector<T> Selector<T>::intersection(const T &view) const
{
    MoleculeView::assertSameMolecule(view);

    if (idxs_set.contains(view.index()))
    {
        Selector<T> ret(*this);
        ret.idxs.clear();
        ret.idxs_set.clear();
        ret._pvt_add(view.index());
        return ret;
    }
    else
        return Selector<T>(this->data());
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
    QList<typename T::Index> other_idxs = id.map(this->d->info());

    Selector<T> ret(*this);

    if (idxs.count() <= other_idxs.count())
    {
        auto other_idxs_set = convert_to_qset(other_idxs);

        foreach (typename T::Index idx, idxs)
        {
            if (not other_idxs_set.contains(idx))
                ret._pvt_sub(idx);
        }
    }
    else
    {
        foreach (typename T::Index idx, other_idxs)
        {
            if (not idxs_set.contains(idx))
                ret._pvt_sub(idx);
        }
    }

    return ret;
}

/** Return the set that has a completely inverted selection */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Selector<T> Selector<T>::invert() const
{
    QList<typename T::Index> all_idxs = detail::getAll<T>(this->d->info());

    Selector<T> ret(*this);
    ret.idxs.clear();
    ret.idxs_set.clear();

    foreach (typename T::Index idx, all_idxs)
    {
        if (not idxs_set.contains(idx))
            ret._pvt_add(idx);
    }

    return ret;
}

/** Return whether this set contains all of the views
    in 'other' */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool Selector<T>::contains(const Selector<T> &other) const
{
    if (idxs_set.count() < other.idxs_set.count())
        return false;
    else
    {
        foreach (typename T::Index idx, other.idxs)
        {
            if (not idxs_set.contains(idx))
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
    MoleculeView::assertSameMolecule(other);

    if (idxs_set.count() <= other.idxs_set.count())
    {
        foreach (typename T::Index idx, idxs)
        {
            if (other.idxs_set.contains(idx))
                return true;
        }
    }
    else
    {
        foreach (typename T::Index idx, other.idxs)
        {
            if (idxs_set.contains(idx))
                return true;
        }
    }

    return false;
}

/** Return whether or not this set contains the view 'view' */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool Selector<T>::contains(const T &view) const
{
    MoleculeView::assertSameMolecule(view);
    return idxs_set.contains(view.index());
}

/** Return whether or not this set contains all of the
    view identified by the ID 'id' */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool Selector<T>::contains(const typename T::ID &id) const
{
    QList<typename T::Index> id_idxs = id.map(this->d->info());

    if (idxs.count() < id_idxs.count())
        return false;
    else
    {
        foreach (typename T::Index idx, id_idxs)
        {
            if (not idxs_set.contains(idx))
                return false;
        }

        return true;
    }
}

/** Return whether or not this set contains the view 'view' */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool Selector<T>::intersects(const T &view) const
{
    return this->isSameMolecule(view) and idxs.contains(view.index());
}

/** Return whether this set contains some of the views
    identified by the ID 'id' */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool Selector<T>::intersects(const typename T::ID &id) const
{
    QList<typename T::Index> id_idxs = id.map(this->d->info());

    if (idxs.count() <= id_idxs.count())
    {
        auto id_idxs_set = convert_to_qset(id_idxs);

        foreach (typename T::Index idx, idxs)
        {
            if (id_idxs_set.contains(idx))
                return true;
        }
    }
    else
    {
        foreach (typename T::Index idx, id_idxs)
        {
            if (idxs_set.contains(idx))
                return true;
        }
    }

    return false;
}

/** Return all of the atoms selected in this set */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomSelection Selector<T>::selection() const
{
    AtomSelection selected_atoms(*(this->d));

    selected_atoms.selectNone();

    foreach (typename T::Index idx, idxs)
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
    AtomSelection selected_atoms(*(this->d));

    selected_atoms.selectOnly( idxs.at( SireID::Index(i).map(idxs.count()) ) );

    return selected_atoms;
}

/** Return the selection of the atoms in the ith to jth views

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
AtomSelection Selector<T>::selection(int i, int j) const
{
    AtomSelection selected_atoms(*(this->d));
    selected_atoms.selectNone();

    i = SireID::Index(i).map(idxs.count());
    j = SireID::Index(j).map(idxs.count());

    if (i > j)
        qSwap(i,j);

    for ( ; i<=j; ++i)
    {
        selected_atoms = selected_atoms.select( idxs.at(i) );
    }

    return selected_atoms;
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
    return Mover< Selector<T> >(*this, this->selection(i));
}

/** Return an object that can move the ith to jth views in this set */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Mover< Selector<T> > Selector<T>::move(int i, int j) const
{
    return Mover< Selector<T> >(*this, this->selection(i,j));
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
    return Evaluator(*this, this->selection(i));
}

/** Return an evaluator that can evaluate properties over
    the ith to jth views in this set */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Evaluator Selector<T>::evaluate(int i, int j) const
{
    return Evaluator(*this, this->selection(i,j));
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
    return Selector<T>(*(this->d), this->selection(i));
}

/** Return a selector for the ith to jth views

    \throw SireError::invalid_index
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
Selector<T> Selector<T>::selector(int i, int j) const
{
    return Selector<T>(*(this->d), this->selection(i,j));
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
    return detail::get_property<V>(ptr, this->data(), this->idxs, key);
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
                                   this->idxs, metakey);
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
    return detail::get_metadata<V>(ptr, this->data(), this->idxs,
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
    detail::set_property<V>(ptr, this->data(), this->idxs, key, values);
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
    detail::set_metadata<V>(ptr, this->data(), this->idxs, metakey, values);
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
    detail::set_metadata<V>(ptr, this->data(), this->idxs, key, metakey, values);
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
    detail::set_property<V>(ptr, this->data(), this->idxs, key, value);
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
    detail::set_metadata<V>(ptr, this->data(), this->idxs, metakey, value);
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
    detail::set_metadata<V>(ptr, this->data(), this->idxs, key, metakey, value);
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

}

SIRE_END_HEADER

#endif
