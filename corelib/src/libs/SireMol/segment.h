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

#ifndef SIREMOL_SEGMENT_H
#define SIREMOL_SEGMENT_H

#include "moleculeview.h"
#include "segproperty.hpp"
#include "atomselection.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class Segment;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::Segment&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::Segment&);

namespace SireMol
{

class SegID;
class SegIdx;
class SegName;

class Evaluator;
class SegEditor;

template<class T>
class Mover;

template<class T>
class Selector;

class Atom;
class Molecule;

/** This is a view of a single segment within a molecule

    @author Christopher Woods
*/
class SIREMOL_EXPORT Segment : public SireBase::ConcreteProperty<Segment,MoleculeView>
{

friend SIREMOL_EXPORT QDataStream& ::operator<<(QDataStream&, const Segment&);
friend SIREMOL_EXPORT QDataStream& ::operator>>(QDataStream&, Segment&);

public:
    typedef SegID ID;
    typedef SegIdx Index;
    typedef SegName Name;
    typedef SegIdx Number;

    Segment();

    Segment(const MoleculeData &data, const SegID &segid);

    Segment(const Segment &other);

    ~Segment();

    Segment& operator=(const Segment &other);

    bool operator==(const Segment &other) const;
    bool operator!=(const Segment &other) const;

    static const char* typeName();

    Segment* clone() const;

    QString toString() const;

    MolViewPtr toSelector() const;

    bool isEmpty() const;
    bool selectedAll() const;

    AtomSelection selection() const;

    void update(const MoleculeData &moldata);

    const SegName& name() const;
    SegIdx index() const;
    SegIdx number() const;

    bool hasProperty(const PropertyName &key) const;
    bool hasMetadata(const PropertyName &metakey) const;
    bool hasMetadata(const PropertyName &key,
                     const PropertyName &metakey) const;

    QStringList propertyKeys() const;
    QStringList metadataKeys() const;
    QStringList metadataKeys(const PropertyName &key) const;

    QVariant propertyAsVariant(const PropertyName &key) const;
    SireBase::PropertyPtr propertyAsProperty(const PropertyName &key) const;

    template<class T>
    const T& property(const PropertyName &key) const;

    template<class T>
    const T& metadata(const PropertyName &metakey) const;

    template<class T>
    const T& metadata(const PropertyName &key,
                      const PropertyName &metakey) const;

    Mover<Segment> move() const;
    Evaluator evaluate() const;
    SegEditor edit() const;
    Selector<Segment> selector() const;
    Selector<Segment> invert() const;

    int nAtoms() const;

    const QList<AtomIdx>& atomIdxs() const;

    bool contains(AtomIdx atomidx) const;
    bool contains(const AtomID &atomid) const;
    bool intersects(const AtomID &atomid) const;

    void assertContainsProperty(const PropertyName &key) const;

    void assertContainsMetadata(const PropertyName &metakey) const;
    void assertContainsMetadata(const PropertyName &key,
                                const PropertyName &metakey) const;

protected:
    template<class T>
    void setProperty(const QString &key, const T &value);

    template<class T>
    void setMetadata(const QString &metakey, const T &value);

    template<class T>
    void setMetadata(const QString &key, const QString &metakey,
                     const T &value);

private:
    /** The index of the segment in the molecule */
    SegIdx segidx;

    /** The atoms that are part of this segment */
    AtomSelection selected_atoms;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Return the property (of type T) at key 'key' that is
    specifically assigned to this segment. This will only work
    if the property at this key is a segment property (i.e.
    has one value for every segment) and that it can be
    cast to type T

    \throw SireMol::missing_property
    \throw SireError::invalid_cast
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T& Segment::property(const PropertyName &key) const
{
    const Property &property = d->property(key);
    const SegProperty<T> &seg_props = property.asA< SegProperty<T> >();
    return seg_props.at(this->index());
}

/** Return the metadata at metakey 'metakey' for this residue

    \throw SireMol::missing_property
    \throw SireError::invalid_cast
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T& Segment::metadata(const PropertyName &metakey) const
{
    const Property &property = d->metadata(metakey);
    const SegProperty<T> &seg_props = property.asA< SegProperty<T> >();
    return seg_props.at(this->index());
}

/** Return the metadata at metakey 'metakey' for the property
    at key 'key'

    \throw SireMol::missing_property
    \throw SireError::invalid_cast
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T& Segment::metadata(const PropertyName &key,
                           const PropertyName &metakey) const
{
    const Property &property = d->metadata(key, metakey);
    const SegProperty<T> &seg_props = property.asA< SegProperty<T> >();
    return seg_props.at(this->index());
}

/** Set the property (of type T) at key 'key' for this
    segment to be equal to 'value'. This works by creating
    a SrgProperty<T> for this molecule, and assigning
    the value for this segment to 'value'. If there is already
    a property at key 'key', then it must be of type
    SegProperty<T> for this to work

    \throw SireMol::invalid_cast
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void Segment::setProperty(const QString &key, const T &value)
{
    MoleculeView::setProperty<SegIdx,SegProperty<T>,T>(*d, key, this->index(),
                                                       value);
}

/** Set the metadata at metakey 'metakey' to the value 'value'
    for this residue

    \throw SireError::invalid_cast
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void Segment::setMetadata(const QString &metakey, const T &value)
{
    MoleculeView::setMetadata<SegIdx,SegProperty<T>,T>(*d, metakey, this->index(),
                                                       value);
}

/** Set the metadata at metakey 'metakey' for the property at key
    'key' to the value 'value'

    \throw SireError::invalid_cast
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void Segment::setMetadata(const QString &key, const QString &metakey,
                           const T &value)
{
    MoleculeView::setMetadata<SegIdx,SegProperty<T>,T>(*d, key, metakey, this->index(),
                                                       value);
}

namespace detail
{

void assertSameSize(Segment*, int nres, int nprops);

template<>
SIRE_ALWAYS_INLINE int getCount<Segment>(const MolInfo &molinfo)
{
    return molinfo.nSegments();
}

template<>
SIRE_ALWAYS_INLINE QList<SegIdx> getAll<Segment>(const MolInfo &molinfo)
{
    return molinfo.getSegments();
}

template<>
SIRE_ALWAYS_INLINE QList<SegIdx> getAll<Segment>(const MolInfo &molinfo,
                                     const AtomSelection &selected_atoms)
{
    molinfo.assertCompatibleWith(selected_atoms);
    return selected_atoms.selectedSegments();
}

template<class V>
SIRE_OUTOFLINE_TEMPLATE
QList<V> get_property(Segment*, const MoleculeData &moldata,
                      const QList<Segment::Index> &idxs,
                      const PropertyName &key)
{
    return get_property<SegProperty<V>,Segment::Index,V>(moldata,idxs,key);
}

template<class V>
SIRE_OUTOFLINE_TEMPLATE
QList<V> get_metadata(Segment*, const MoleculeData &moldata,
                      const QList<Segment::Index> &idxs,
                      const PropertyName &metakey)
{
    return get_metadata<SegProperty<V>,Segment::Index,V>(moldata,idxs,metakey);
}

template<class V>
SIRE_OUTOFLINE_TEMPLATE
QList<V> get_metadata(Segment*, const MoleculeData &moldata,
                      const QList<Segment::Index> &idxs,
                      const PropertyName &key,
                      const PropertyName &metakey)
{
    return get_metadata<SegProperty<V>,Segment::Index,V>(moldata,idxs,key,metakey);
}

template<class V>
SIRE_OUTOFLINE_TEMPLATE
void set_property(Segment *ptr, MoleculeData &moldata,
                  const QList<Segment::Index> &idxs,
                  const QString &key,
                  const QList<V> &values)
{
    assertSameSize(ptr, idxs.count(), values.count());
    set_property<SegProperty<V>,Segment::Index,V>(moldata,idxs,key,values);
}

template<class V>
SIRE_OUTOFLINE_TEMPLATE
void set_metadata(Segment *ptr, MoleculeData &moldata,
                  const QList<Segment::Index> &idxs,
                  const QString &metakey,
                  const QList<V> &values)
{
    assertSameSize(ptr, idxs.count(), values.count());
    set_metadata<SegProperty<V>,Segment::Index,V>(moldata,idxs,metakey,values);
}

template<class V>
SIRE_OUTOFLINE_TEMPLATE
void set_metadata(Segment *ptr, MoleculeData &moldata,
                  const QList<Segment::Index> &idxs,
                  const QString &key, const QString &metakey,
                  const QList<V> &values)
{
    assertSameSize(ptr, idxs.count(), values.count());

    set_metadata<SegProperty<V>,Segment::Index,V>(moldata,idxs,key,metakey,values);
}

template<class V>
SIRE_OUTOFLINE_TEMPLATE
void set_property(Segment*, MoleculeData &moldata,
                  const QList<Segment::Index> &idxs,
                  const QString &key,
                  const V &value)
{
    set_property<SegProperty<V>,Segment::Index,V>(moldata,idxs,key,value);
}

template<class V>
SIRE_OUTOFLINE_TEMPLATE
void set_metadata(Segment*, MoleculeData &moldata,
                  const QList<Segment::Index> &idxs,
                  const QString &metakey,
                  const V &value)
{
    set_metadata<SegProperty<V>,Segment::Index,V>(moldata,idxs,metakey,value);
}

template<class V>
SIRE_OUTOFLINE_TEMPLATE
void set_metadata(Segment*, MoleculeData &moldata,
                  const QList<Segment::Index> &idxs,
                  const QString &key, const QString &metakey,
                  const V &value)
{
    set_metadata<SegProperty<V>,Segment::Index,V>(moldata,idxs,key,metakey,value);
}

SIREMOL_EXPORT bool has_property(const Segment*, const MoleculeData &moldata,
                  const PropertyName &key);

SIREMOL_EXPORT bool has_metadata(const Segment*, const MoleculeData &moldata,
                  const PropertyName &metakey);

SIREMOL_EXPORT bool has_metadata(const Segment*, const MoleculeData &moldata,
                  const PropertyName &key, const PropertyName &metakey);

} //end of namespace detail

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

Q_DECLARE_METATYPE(SireMol::Segment);
Q_DECLARE_METATYPE(SireMol::Mover<SireMol::Segment>);
Q_DECLARE_METATYPE(SireMol::Selector<SireMol::Segment>);

//Q_DECLARE_METATYPE(SireMol::Mover< SireMol::Selector<SireMol::Segment> >);

SIRE_EXPOSE_CLASS( SireMol::Segment )

SIRE_EXPOSE_ALIAS( SireMol::Mover<SireMol::Segment>, SireMol::Mover_Segment_ )
SIRE_EXPOSE_ALIAS( SireMol::Selector<SireMol::Segment>, SireMol::Selector_Segment_ )

#ifdef SIRE_INSTANTIATE_TEMPLATES

#include "mover.hpp"

template class SireMol::Mover<SireMol::Segment>;

#endif //SIRE_INSTANTIATE_TEMPLATES

SIRE_END_HEADER

#endif
