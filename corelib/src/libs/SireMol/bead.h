/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2010  Christopher Woods
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

#ifndef SIREMOL_BEAD_H
#define SIREMOL_BEAD_H

#include "moleculeview.h"
#include "beading.h"
#include "beadidx.h"
#include "beadproperty.hpp"
#include "atomselection.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class Bead;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::Bead&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::Bead&);

namespace SireMol
{

class Evaluator;
class BeadEditor;

template<class T>
class Mover;

template<class T>
class Selector;

class Atom;
class Molecule;

using SireBase::Property;

/** A Bead is a group of atoms (defined using a SireMol::Beading function)
    within a molecule. Beads can be used for coarse-graining, or for
    implementing group-based cutoffs

    @author Christopher Woods
*/
class SIREMOL_EXPORT Bead : public SireBase::ConcreteProperty<Bead,MoleculeView>
{

friend SIREMOL_EXPORT QDataStream& ::operator<<(QDataStream&, const Bead&);
friend SIREMOL_EXPORT QDataStream& ::operator>>(QDataStream&, Bead&);

public:
    typedef BeadID ID;
    typedef BeadIdx Index;

    Bead();
    Bead(const MoleculeData &moldata, const BeadIdx &bead,
         const PropertyMap &map = PropertyMap());

    Bead(const Bead &other);

    ~Bead();

    Bead& operator=(const Bead &other);

    bool operator==(const Bead &other) const;
    bool operator!=(const Bead &other) const;

    static const char* typeName();

    Bead* clone() const;

    MolViewPtr operator[](int i) const;
    int nViews() const;

    Atom atom(int i) const;

    QString toString() const;

    bool isEmpty() const;
    bool selectedAll() const;

    MolViewPtr toSelector() const;

    AtomSelection selection() const;

    void update(const MoleculeData &moldata);

    BeadIdx index() const;

    Beads beads() const;

    bool hasProperty(const PropertyName &key) const;
    bool hasMetadata(const PropertyName &metakey) const;
    bool hasMetadata(const PropertyName &key,
                     const PropertyName &metakey) const;

    QStringList propertyKeys() const;
    QStringList metadataKeys() const;
    QStringList metadataKeys(const PropertyName &key) const;

    template<class T>
    const T& property(const PropertyName &key) const;

    template<class T>
    const T& metadata(const PropertyName &metakey) const;

    template<class T>
    const T& metadata(const PropertyName &key,
                      const PropertyName &metakey) const;

    Mover<Bead> move() const;
    Evaluator evaluate() const;
    BeadEditor edit() const;

    int nAtoms() const;

    QList<AtomIdx> atomIdxs() const;

    const Beading& beading() const;

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

    friend class Beads;
    Bead(const MoleculeData &moldata, BeadIdx beadidx,
         const Beading &beading, const PropertyName &beading_property);

private:
    /** The index of the bead in the molecule */
    BeadIdx beadidx;

    /** The beading used to create the beads */
    BeadingPtr bdng;

    /** The location of the beading property */
    PropertyName beading_property;

    /** The atoms that are selected as part of this residue */
    AtomSelection selected_atoms;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Return the property (of type T) at key 'key' that is
    specifically assigned to this bead. This will only work
    if the property at this key is a bead property (i.e.
    has one value for every bead) and that it can be
    cast to type T

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T& Bead::property(const PropertyName &key) const
{
    const Property &property = d->property(key);
    const BeadProperty<T> &bead_props = property.asA< BeadProperty<T> >();
    return bead_props.at(this->index());
}

/** Return the metadata at metakey 'metakey' for this bead

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T& Bead::metadata(const PropertyName &metakey) const
{
    const Property &property = d->metadata(metakey);
    const BeadProperty<T> &bead_props = property.asA< BeadProperty<T> >();
    return bead_props.at(this->index());
}

/** Return the metadata at metakey 'metakey' for the property
    at key 'key'

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T& Bead::metadata(const PropertyName &key,
                        const PropertyName &metakey) const
{
    const Property &property = d->metadata(key, metakey);
    const BeadProperty<T> &bead_props = property.asA< BeadProperty<T> >();
    return bead_props.at(this->index());
}

/** Set the property (of type T) at key 'key' for this
    residue to be equal to 'value'. This works by creating
    a ResProperty<T> for this molecule, and assigning
    the value for this residue to 'value'. If there is already
    a property at key 'key', then it must be of type
    ResProperty<T> for this to work

    \throw SireMol::invalid_cast
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void Bead::setProperty(const QString &key, const T &value)
{
    BeadProperty<T> props;
    MoleculeData &data = *d;

    if (data.hasProperty(key))
    {
        //take the property to prevent unnecessary copying caused
        //by implicit sharing of the property
        PropertyPtr old_property = data.takeProperty(key);
        props = old_property->asA< BeadProperty<T> >();
    }
    else
        props = BeadProperty<T>(data.info(), bdng.read());

    props.set(this->index(), value);

    data.setProperty(key, props);
}

/** Set the metadata at metakey 'metakey' to the value 'value'
    for this residue

    \throw SireError::invalid_cast
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void Bead::setMetadata(const QString &metakey, const T &value)
{
    BeadProperty<T> props;
    MoleculeData &data = *d;

    if (data.hasMetadata(metakey))
    {
        //take the metadata to prevent unnecessary copying caused
        //by implicit sharing of the property
        PropertyPtr old_metadata = data.takeMetadata(metakey);
        props = old_metadata->asA< BeadProperty<T> >();
    }
    else
        props = BeadProperty<T>(data.info(), bdng.read());

    props.set(this->index(), value);

    data.setMetadata(metakey, props);
}

/** Set the metadata at metakey 'metakey' for the property at key
    'key' to the value 'value'

    \throw SireError::invalid_cast
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void Bead::setMetadata(const QString &key, const QString &metakey,
                       const T &value)
{
    BeadProperty<T> props;
    MoleculeData &data = *d;

    if (data.hasMetadata(key, metakey))
    {
        //take the metadata to prevent unnecessary copying caused
        //by implicit sharing of the property
        PropertyPtr old_metadata = data.takeMetadata(key, metakey);
        props = old_metadata->asA< BeadProperty<T> >();
    }
    else
        props = BeadProperty<T>(data.info(),bdng.read());

    props.set(this->index(), value);

    data.setMetadata(key, metakey, props);
}

namespace detail
{

void assertSameSize(Bead*, int nres, int nprops);

template<class V>
SIRE_OUTOFLINE_TEMPLATE
QList<V> get_property(Bead*, const MoleculeData &moldata,
                      const QList<Bead::Index> &idxs,
                      const PropertyName &key)
{
    return get_property<BeadProperty<V>,Bead::Index,V>(moldata,idxs,key);
}

template<class V>
SIRE_OUTOFLINE_TEMPLATE
QList<V> get_metadata(Bead*, const MoleculeData &moldata,
                      const QList<Bead::Index> &idxs,
                      const PropertyName &metakey)
{
    return get_metadata<BeadProperty<V>,Bead::Index,V>(moldata,idxs,metakey);
}

template<class V>
SIRE_OUTOFLINE_TEMPLATE
QList<V> get_metadata(Bead*, const MoleculeData &moldata,
                      const QList<Bead::Index> &idxs,
                      const PropertyName &key,
                      const PropertyName &metakey)
{
    return get_metadata<BeadProperty<V>,Bead::Index,V>(moldata,idxs,key,metakey);
}

template<class V>
SIRE_OUTOFLINE_TEMPLATE
void set_property(Bead *ptr, MoleculeData &moldata,
                  const QList<Bead::Index> &idxs,
                  const QString &key,
                  const QList<V> &values)
{
    assertSameSize(ptr, idxs.count(), values.count());
    set_property<BeadProperty<V>,Bead::Index,V>(moldata,idxs,key,values);
}

template<class V>
SIRE_OUTOFLINE_TEMPLATE
void set_metadata(Bead *ptr, MoleculeData &moldata,
                  const QList<Bead::Index> &idxs,
                  const QString &metakey,
                  const QList<V> &values)
{
    assertSameSize(ptr, idxs.count(), values.count());
    set_metadata<BeadProperty<V>,Bead::Index,V>(moldata,idxs,metakey,values);
}

template<class V>
SIRE_OUTOFLINE_TEMPLATE
void set_metadata(Bead *ptr, MoleculeData &moldata,
                  const QList<Bead::Index> &idxs,
                  const QString &key, const QString &metakey,
                  const QList<V> &values)
{
    assertSameSize(ptr, idxs.count(), values.count());

    set_metadata<BeadProperty<V>,Bead::Index,V>(moldata,idxs,key,metakey,values);
}

template<class V>
SIRE_OUTOFLINE_TEMPLATE
void set_property(Bead*, MoleculeData &moldata,
                  const QList<Bead::Index> &idxs,
                  const QString &key,
                  const V &value)
{
    set_property<BeadProperty<V>,Bead::Index,V>(moldata,idxs,key,value);
}

template<class V>
SIRE_OUTOFLINE_TEMPLATE
void set_metadata(Bead*, MoleculeData &moldata,
                  const QList<Bead::Index> &idxs,
                  const QString &metakey,
                  const V &value)
{
    set_metadata<BeadProperty<V>,Bead::Index,V>(moldata,idxs,metakey,value);
}

template<class V>
SIRE_OUTOFLINE_TEMPLATE
void set_metadata(Bead*, MoleculeData &moldata,
                  const QList<Bead::Index> &idxs,
                  const QString &key, const QString &metakey,
                  const V &value)
{
    set_metadata<BeadProperty<V>,Bead::Index,V>(moldata,idxs,key,metakey,value);
}

SIREMOL_EXPORT bool has_property(const Bead*, const MoleculeData &moldata,
                  const PropertyName &key);

SIREMOL_EXPORT bool has_metadata(const Bead*, const MoleculeData &moldata,
                  const PropertyName &metakey);

SIREMOL_EXPORT bool has_metadata(const Bead*, const MoleculeData &moldata,
                  const PropertyName &key, const PropertyName &metakey);


} //end of namespace detail

#endif //SIRE_SKIP_INLINE_FUNCTIONS

} //end of namespace SireMol

Q_DECLARE_METATYPE(SireMol::Bead);
Q_DECLARE_METATYPE(SireMol::Mover<SireMol::Bead>);

SIRE_EXPOSE_CLASS( SireMol::Bead )

SIRE_EXPOSE_ALIAS( SireMol::Mover<SireMol::Bead>, SireMol::Mover_Bead_ )

#ifdef SIRE_INSTANTIATE_TEMPLATES

#include "mover.hpp"

template class SireMol::Mover<SireMol::Bead>;

#endif //SIRE_INSTANTIATE_TEMPLATES

#endif
