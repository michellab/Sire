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

#ifndef SIREMOL_CHAIN_H
#define SIREMOL_CHAIN_H

#include "moleculeview.h"
#include "chainproperty.hpp"
#include "atomselection.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class Chain;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::Chain&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::Chain&);

namespace SireMol
{

class ChainID;
class ChainIdx;
class ChainName;
class ChainEditor;

class Evaluator;

template<class T>
class Mover;

template<class T>
class Selector;

class Atom;
class Residue;
class Molecule;

/** This class represents a Chain in a Molecule.

    @author Christopher Woods
*/
class SIREMOL_EXPORT Chain : public ConcreteProperty<Chain,MoleculeView>
{

friend SIREMOL_EXPORT QDataStream& ::operator<<(QDataStream&, const Chain&);
friend SIREMOL_EXPORT QDataStream& ::operator>>(QDataStream&, Chain&);

public:
    typedef ChainID ID;
    typedef ChainIdx Index;
    typedef ChainName Name;

    Chain();
    
    Chain(const MoleculeData &moldata, const ChainID &chainid);

    Chain(const Chain &other);

    ~Chain();

    Chain& operator=(const Chain &other);
    
    static const char* typeName();
    
    Chain* clone() const;
    
    bool operator==(const Chain &other) const;
    bool operator!=(const Chain &other) const;

    QString toString() const;
    
    bool isEmpty() const;
    bool selectedAll() const;
    
    AtomSelection selection() const;
    
    void update(const MoleculeData &moldata);
    
    ChainName name() const;
    ChainIdx index() const;
    
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
    
    Mover<Chain> move() const;
    Evaluator evaluate() const;
    ChainEditor edit() const;
    Selector<Chain> selector() const;
    
    int nAtoms() const;
    
    QList<AtomIdx> atomIdxs() const;
    
    bool contains(AtomIdx atomidx) const;
    bool contains(const AtomID &atomid) const;
    bool intersects(const AtomID &atomid) const;
    
    int nResidues() const;
    
    const QList<ResIdx>& resIdxs() const;
    
    bool contains(ResIdx residx) const;
    bool contains(const ResID &resid) const;
    bool intersects(const ResID &resid) const;

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
    /** The index of the Chain in the molecule */
    ChainIdx chainidx;
    
    /** The atoms that are selected as part of this Chain */
    AtomSelection selected_atoms;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Return the property (of type T) at key 'key' that is 
    specifically assigned to this chain. This will only work
    if the property at this key is a chain property (i.e.
    has one value for every chain) and that it can be
    cast to type T
    
    \throw SireMol::missing_property
    \throw SireError::invalid_cast
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T& Chain::property(const PropertyName &key) const
{
    const Property &property = d->property(key);
    const ChainProperty<T> &chain_props = property.asA< ChainProperty<T> >();
    return chain_props.at(this->index());
}

/** Return the metadata at metakey 'metakey' for this residue

    \throw SireMol::missing_property
    \throw SireError::invalid_cast
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T& Chain::metadata(const PropertyName &key) const
{
    const Property &property = d->metadata(key);
    const ChainProperty<T> &chain_props = property.asA< ChainProperty<T> >();
    return chain_props.at(this->index());
}

/** Return the metadata at metakey 'metakey' for the property
    at key 'key'
    
    \throw SireMol::missing_property
    \throw SireError::invalid_cast
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
const T& Chain::metadata(const PropertyName &key, const PropertyName &metakey) const
{
    const Property &property = d->metadata(key, metakey);
    const ChainProperty<T> &chain_props = property.asA< ChainProperty<T> >();
                                                
    return chain_props.at(this->index());
}

/** Set the property (of type T) at key 'key' for this
    chain to be equal to 'value'. This works by creating
    a ChainProperty<T> for this molecule, and assigning
    the value for this chain to 'value'. If there is already
    a property at key 'key', then it must be of type 
    ChainProperty<T> for this to work
    
    \throw SireMol::invalid_cast
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void Chain::setProperty(const QString &key, const T &value)
{
    MoleculeView::setProperty<ChainIdx,ChainProperty<T>,T>(*d, key, this->index(),
                                                           value);
}

/** Set the metadata at metakey 'metakey' to the value 'value' 
    for this residue
    
    \throw SireError::invalid_cast
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void Chain::setMetadata(const QString &metakey, const T &value)
{
    MoleculeView::setMetadata<ChainIdx,ChainProperty<T>,T>(*d, metakey, this->index(),
                                                           value);
}

/** Set the metadata at metakey 'metakey' for the property at key
    'key' to the value 'value'
    
    \throw SireError::invalid_cast
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
void Chain::setMetadata(const QString &key, const QString &metakey,
                        const T &value)
{
    MoleculeView::setMetadata<ChainIdx,ChainProperty<T>,T>(*d, key, metakey, 
                                                           this->index(), value);
}

namespace detail
{

void assertSameSize(Chain*, int nres, int nprops);
    
template<>
SIRE_ALWAYS_INLINE QList<ChainIdx> getAll<Chain>(const MolInfo &molinfo)
{
    return molinfo.getChains();
}

template<>
SIRE_ALWAYS_INLINE QList<ChainIdx> getAll<Chain>(const MolInfo &molinfo,
                                     const AtomSelection &selected_atoms)
{
    molinfo.assertCompatibleWith(selected_atoms);
    return selected_atoms.selectedChains();
}

template<class V>
SIRE_OUTOFLINE_TEMPLATE
QList<V> get_property(Chain*, const MoleculeData &moldata,
                      const QList<Chain::Index> &idxs,
                      const PropertyName &key)
{
    return get_property<ChainProperty<V>,Chain::Index,V>(moldata,idxs,key);
}

template<class V>
SIRE_OUTOFLINE_TEMPLATE
QList<V> get_metadata(Chain*, const MoleculeData &moldata,
                      const QList<Chain::Index> &idxs,
                      const PropertyName &metakey)
{
    return get_metadata<ChainProperty<V>,Chain::Index,V>(moldata,idxs,metakey);
}

template<class V>
SIRE_OUTOFLINE_TEMPLATE
QList<V> get_metadata(Chain*, const MoleculeData &moldata,
                      const QList<Chain::Index> &idxs,
                      const PropertyName &key,
                      const PropertyName &metakey)
{
    return get_metadata<ChainProperty<V>,Chain::Index,V>(moldata,idxs,key,metakey);
}

template<class V>
SIRE_OUTOFLINE_TEMPLATE
void set_property(Chain *ptr, MoleculeData &moldata,
                  const QList<Chain::Index> &idxs,
                  const QString &key,
                  const QList<V> &values)
{
    assertSameSize(ptr, idxs.count(), values.count());
    set_property<ChainProperty<V>,Chain::Index,V>(moldata,idxs,key,values);
}

template<class V>
SIRE_OUTOFLINE_TEMPLATE
void set_metadata(Chain *ptr, MoleculeData &moldata,
                  const QList<Chain::Index> &idxs,
                  const QString &metakey,
                  const QList<V> &values)
{
    assertSameSize(ptr, idxs.count(), values.count());
    set_metadata<ChainProperty<V>,Chain::Index,V>(moldata,idxs,metakey,values);
}

template<class V>
SIRE_OUTOFLINE_TEMPLATE
void set_metadata(Chain *ptr, MoleculeData &moldata,
                  const QList<Chain::Index> &idxs,
                  const QString &key, const QString &metakey,
                  const QList<V> &values)
{
    assertSameSize(ptr, idxs.count(), values.count());
    
    set_metadata<ChainProperty<V>,Chain::Index,V>(moldata,idxs,key,metakey,values);
}

template<class V>
SIRE_OUTOFLINE_TEMPLATE
void set_property(Chain*, MoleculeData &moldata,
                  const QList<Chain::Index> &idxs,
                  const QString &key,
                  const V &value)
{
    set_property<ChainProperty<V>,Chain::Index,V>(moldata,idxs,key,value);
}

template<class V>
SIRE_OUTOFLINE_TEMPLATE
void set_metadata(Chain*, MoleculeData &moldata,
                  const QList<Chain::Index> &idxs,
                  const QString &metakey,
                  const V &value)
{
    set_metadata<ChainProperty<V>,Chain::Index,V>(moldata,idxs,metakey,value);
}

template<class V>
SIRE_OUTOFLINE_TEMPLATE
void set_metadata(Chain*, MoleculeData &moldata,
                  const QList<Chain::Index> &idxs,
                  const QString &key, const QString &metakey,
                  const V &value)
{
    set_metadata<ChainProperty<V>,Chain::Index,V>(moldata,idxs,key,metakey,value);
}

SIREMOL_EXPORT bool has_property(const Chain*, const MoleculeData &moldata,
                  const PropertyName &key);
                  
SIREMOL_EXPORT bool has_metadata(const Chain*, const MoleculeData &moldata,
                  const PropertyName &metakey);
                  
SIREMOL_EXPORT bool has_metadata(const Chain*, const MoleculeData &moldata,
                  const PropertyName &key, const PropertyName &metakey);                 

} //end of namespace detail

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

Q_DECLARE_METATYPE(SireMol::Chain);
Q_DECLARE_METATYPE(SireMol::Mover<SireMol::Chain>);
Q_DECLARE_METATYPE(SireMol::Selector<SireMol::Chain>);

//Q_DECLARE_METATYPE(SireMol::Mover< SireMol::Selector<SireMol::Chain> >);

SIRE_EXPOSE_CLASS( SireMol::Chain )

SIRE_EXPOSE_ALIAS( SireMol::Mover<SireMol::Chain>, SireMol::Mover_Chain_ )
SIRE_EXPOSE_ALIAS( SireMol::Selector<SireMol::Chain>, SireMol::Selector_Chain_ )

#ifdef SIRE_INSTANTIATE_TEMPLATES

#include "mover.hpp"

template class SireMol::Mover<SireMol::Chain>;

#endif //SIRE_INSTANTIATE_TEMPLATES


SIRE_END_HEADER

#endif

