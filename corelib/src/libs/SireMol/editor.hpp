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

#ifndef SIREMOL_EDITOR_HPP
#define SIREMOL_EDITOR_HPP

#include "moleculeview.h"

SIRE_BEGIN_HEADER

namespace SireMol
{

class AtomEditor;
class CGEditor;
class ResEditor;
class ChainEditor;
class SegEditor;
class MolEditor;

/** This is the class that is used to edit the molecule view
    of type 'T'

    @author Christopher Woods
*/
template<class Parent, class T>
class Editor : public T
{
public:
    ~Editor();

    Editor<Parent,T>& operator=(const Editor<Parent,T> &other);
    Editor<Parent,T>& operator=(const T &other);

    MolViewPtr operator[](int i);
    MolViewPtr operator[](const QString &name);
    MolViewPtr operator[](const AtomID &atomid);
    MolViewPtr operator[](const ResID &resid);
    MolViewPtr operator[](const CGID &cgid);
    MolViewPtr operator[](const ChainID &chainid);
    MolViewPtr operator[](const SegID &segid);
    MolViewPtr operator[](const SireID::Index &idx);

    template<class V>
    Parent& setProperty(const SireBase::PropertyName &key, const V &value);

    template<class V>
    Parent& setMetadata(const SireBase::PropertyName &metakey, const V &value);

    template<class V>
    Parent& setMetadata(const SireBase::PropertyName &key,
                        const SireBase::PropertyName &metakey,
                        const V &value);

    Parent& removeProperty(const SireBase::PropertyName &key);
    Parent& removeMetadata(const SireBase::PropertyName &metakey);
    Parent& removeMetadata(const SireBase::PropertyName &key,
                           const SireBase::PropertyName &metakey);

    AtomEditor atom();
    AtomEditor atom(int i,
                    const PropertyMap &map = PropertyMap());
    AtomEditor atom(const QString &name,
                    const PropertyMap &map = PropertyMap());
    AtomEditor atom(const AtomID &atomid,
                    const PropertyMap &map = PropertyMap());

    CGEditor cutGroup();
    CGEditor cutGroup(int i,
                      const PropertyMap &map = PropertyMap());
    CGEditor cutGroup(const QString &name,
                      const PropertyMap &map = PropertyMap());
    CGEditor cutGroup(const CGID &cgid,
                      const PropertyMap &map = PropertyMap());

    ResEditor residue();
    ResEditor residue(int i,
                      const PropertyMap &map = PropertyMap());
    ResEditor residue(const QString &name,
                      const PropertyMap &map = PropertyMap());
    ResEditor residue(const ResID &resid,
                      const PropertyMap &map = PropertyMap());

    ChainEditor chain();
    ChainEditor chain(int i,
                      const PropertyMap &map = PropertyMap());
    ChainEditor chain(const QString &name,
                      const PropertyMap &map = PropertyMap());
    ChainEditor chain(const ChainID &chainid,
                      const PropertyMap &map = PropertyMap());

    SegEditor segment();
    SegEditor segment(int i,
                      const PropertyMap &map = PropertyMap());
    SegEditor segment(const QString &name,
                      const PropertyMap &map = PropertyMap());
    SegEditor segment(const SegID &segid,
                      const PropertyMap &map = PropertyMap());

    MolEditor molecule();

    AtomEditor select(const AtomID &atomid,
                      const PropertyMap &map = PropertyMap());

    CGEditor select(const CGID &cgid,
                    const PropertyMap &map = PropertyMap());

    ResEditor select(const ResID &resid,
                     const PropertyMap &map = PropertyMap());

    ChainEditor select(const ChainID &chainid,
                       const PropertyMap &map = PropertyMap());

    SegEditor select(const SegID &segid,
                     const PropertyMap &map = PropertyMap());

protected:
    Editor();
    Editor(const T &view);

    Editor(const Editor<Parent,T> &other);
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

} // end of namespace SireMol

#include "atomeditor.h"
#include "cgeditor.h"
#include "reseditor.h"
#include "chaineditor.h"
#include "segeditor.h"
#include "moleditor.h"

namespace SireMol
{

/** Null constructor */
template<class Parent, class T>
SIRE_OUTOFLINE_TEMPLATE
Editor<Parent, T>::Editor() : T()
{}

/** Construct an editor of the passed view */
template<class Parent, class T>
SIRE_OUTOFLINE_TEMPLATE
Editor<Parent, T>::Editor(const T &view) : T(view)
{}

/** Copy constructor */
template<class Parent, class T>
SIRE_OUTOFLINE_TEMPLATE
Editor<Parent, T>::Editor(const Editor<Parent,T> &other) : T(other)
{}

/** Destructor */
template<class Parent, class T>
SIRE_OUTOFLINE_TEMPLATE
Editor<Parent, T>::~Editor()
{}

/** Copy assignment operator */
template<class Parent, class T>
SIRE_OUTOFLINE_TEMPLATE
Editor<Parent, T>& Editor<Parent, T>::operator=(const Editor<Parent,T> &other)
{
    T::operator=(other);

    return *this;
}

/** Copy assignment from an object of type T */
template<class Parent, class T>
SIRE_OUTOFLINE_TEMPLATE
Editor<Parent, T>& Editor<Parent, T>::operator=(const T &other)
{
    T::operator=(other);

    return *this;
}

template<class Parent, class T>
SIRE_OUTOFLINE_TEMPLATE
MolViewPtr Editor<Parent, T>::operator[](int i)
{
    return this->atom(i);
}

template<class Parent, class T>
SIRE_OUTOFLINE_TEMPLATE
MolViewPtr Editor<Parent, T>::operator[](const QString &key)
{
    return this->atom(key);
}

template<class Parent, class T>
SIRE_OUTOFLINE_TEMPLATE
MolViewPtr Editor<Parent, T>::operator[](const AtomID &atomid)
{
    return this->atom(atomid);
}

/** Return the residue(s) that match 'resid' in this view of the molecule */
template<class Parent, class T>
SIRE_OUTOFLINE_TEMPLATE
MolViewPtr Editor<Parent, T>::operator[](const ResID &resid)
{
    return this->residue(resid);
}

/** Return the CutGroups(s) that match 'resid' in this view of the molecule */
template<class Parent, class T>
SIRE_OUTOFLINE_TEMPLATE
MolViewPtr Editor<Parent, T>::operator[](const CGID &cgid)
{
    return this->cutGroup(cgid);
}

/** Return the residue(s) that match 'resid' in this view of the molecule */
template<class Parent, class T>
SIRE_OUTOFLINE_TEMPLATE
MolViewPtr Editor<Parent, T>::operator[](const ChainID &chainid)
{
    return this->chain(chainid);
}

/** Return the residue(s) that match 'resid' in this view of the molecule */
template<class Parent, class T>
SIRE_OUTOFLINE_TEMPLATE
MolViewPtr Editor<Parent, T>::operator[](const SegID &segid)
{
    return this->segment(segid);
}

/** This is an overload of operator[](int), allowing a SireID::Index to be used
    as the int */
template<class Parent, class T>
SIRE_OUTOFLINE_TEMPLATE
MolViewPtr Editor<Parent, T>::operator[](const SireID::Index &idx)
{
    return this->operator[](idx.value());
}

/** Expose the protected 'T::setProperty()' function */
template<class Parent, class T>
template<class V>
SIRE_OUTOFLINE_TEMPLATE
Parent& Editor<Parent, T>::setProperty(const SireBase::PropertyName &key,
                                       const V &value)
{
    if (key.hasSource())
        T::setProperty(key.source(), value);

    return static_cast<Parent&>(*this);
}

/** Expose the protected 'T::setMetadata()' function */
template<class Parent, class T>
template<class V>
SIRE_OUTOFLINE_TEMPLATE
Parent& Editor<Parent, T>::setMetadata(const SireBase::PropertyName &metakey,
                                       const V &value)
{
    if (metakey.hasSource())
        T::setMetadata(metakey.source(), value);

    return static_cast<Parent&>(*this);
}

/** Expose the protected 'T::setMetadata()' function

    \throw SireBase::missing_property
*/
template<class Parent, class T>
template<class V>
SIRE_OUTOFLINE_TEMPLATE
Parent& Editor<Parent, T>::setMetadata(const SireBase::PropertyName &key,
                                       const SireBase::PropertyName &metakey,
                                       const V &value)
{
    if (key.hasSource() and metakey.hasSource())
        T::setMetadata(key.source(), metakey.source(), value);

    return static_cast<Parent&>(*this);
}

/** Completely remove the property 'key', if this is valid
    property for this view. Note that this will remove this
    property for *all* views, e.g. if this is a Mover<Atom>,
    then this will remove the property if it is an AtomProp,
    and it will remove the property for *all* atoms.

    \throw SireBase::missing_property
*/
template<class Parent, class T>
SIRE_OUTOFLINE_TEMPLATE
Parent& Editor<Parent, T>::removeProperty(const SireBase::PropertyName &key)
{
    if (key.hasSource())
    {
        T::assertContainsProperty(key.source());
        this->d->removeProperty(key.source());
    }

    return static_cast<Parent&>(*this);
}

/** Completely remove the metadata 'metakey', if this is valid
    property for this view. Note that this will remove this
    property for *all* views, e.g. if this is a Mover<Atom>,
    then this will remove the property if it is an AtomProp,
    and it will remove the property for *all* atoms.

    \throw SireBase::missing_property
*/
template<class Parent, class T>
SIRE_OUTOFLINE_TEMPLATE
Parent& Editor<Parent, T>::removeMetadata(const SireBase::PropertyName &metakey)
{
    if (metakey.hasSource())
    {
        T::assertContainsMetadata(metakey.source());
        this->d->removeMetadata(metakey.source());
    }

    return static_cast<Parent&>(*this);
}

/** Completely remove metadata with metakey 'metakey' from
    the property with 'key', if this is valid
    property for this view. Note that this will remove this
    property for *all* views, e.g. if this is a Mover<Atom>,
    then this will remove the property if it is an AtomProp,
    and it will remove the property for *all* atoms.

    \throw SireBase::missing_property
*/
template<class Parent, class T>
SIRE_OUTOFLINE_TEMPLATE
Parent& Editor<Parent, T>::removeMetadata(const SireBase::PropertyName &key,
                                          const SireBase::PropertyName &metakey)
{
    if (key.hasSource() and metakey.hasSource())
    {
        T::assertContainsMetadata(key.source(), metakey.source());
        this->d->removeMetadata(key.source(), metakey.source());
    }

    return static_cast<Parent&>(*this);
}

/** Return the atom of this view - for this to work, only a single
    atom should be contained in this view

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
*/
template<class Parent, class T>
SIRE_OUTOFLINE_TEMPLATE
AtomEditor Editor<Parent, T>::atom()
{
    return AtomEditor( MoleculeView::atom() );
}

template<class Parent, class T>
SIRE_OUTOFLINE_TEMPLATE
AtomEditor Editor<Parent, T>::atom(int i, const PropertyMap &map)
{
    return AtomEditor( MoleculeView::atom(i, map) );
}

template<class Parent, class T>
SIRE_OUTOFLINE_TEMPLATE
AtomEditor Editor<Parent, T>::atom(const QString &name, const PropertyMap &map)
{
    return AtomEditor( MoleculeView::atom(name, map) );
}

/** Return the atom from this view that matches the ID 'atomid'

    \throw SireMol::missing_atom
    \throw SireError::invalid_index
    \throw SireMol::duplicate_atom
*/
template<class Parent, class T>
SIRE_OUTOFLINE_TEMPLATE
AtomEditor Editor<Parent, T>::atom(const AtomID &atomid, const PropertyMap &map)
{
    return AtomEditor( MoleculeView::atom(atomid, map) );
}

/** Return the CutGroup involved with this view - for this
    to work, only a single CutGroup should be involved in this view

    \throw SireMol::missing_cutgroup
    \throw SireMol::duplicate_cutgroup
*/
template<class Parent, class T>
SIRE_OUTOFLINE_TEMPLATE
CGEditor Editor<Parent, T>::cutGroup()
{
    return CGEditor( MoleculeView::cutGroup() );
}

template<class Parent, class T>
SIRE_OUTOFLINE_TEMPLATE
CGEditor Editor<Parent, T>::cutGroup(int i, const PropertyMap &map)
{
    return CGEditor( MoleculeView::cutGroup(i, map) );
}

template<class Parent, class T>
SIRE_OUTOFLINE_TEMPLATE
CGEditor Editor<Parent, T>::cutGroup(const QString &name, const PropertyMap &map)
{
    return CGEditor( MoleculeView::cutGroup(name, map) );
}

/** Return the CutGroups from this view that match the ID 'cgid'

    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
    \throw SireMol::duplicate_cutgroup
*/
template<class Parent, class T>
SIRE_OUTOFLINE_TEMPLATE
CGEditor Editor<Parent, T>::cutGroup(const CGID &cgid, const PropertyMap &map)
{
    return CGEditor( MoleculeView::cutGroup(cgid, map) );
}

/** Return the residue involved with this view - for this
    to work, only a single residue should be involved in this view

    \throw SireMol::missing_residue
    \throw SireMol::duplicate_residue
*/
template<class Parent, class T>
SIRE_OUTOFLINE_TEMPLATE
ResEditor Editor<Parent, T>::residue()
{
    return ResEditor( MoleculeView::residue() );
}

template<class Parent, class T>
SIRE_OUTOFLINE_TEMPLATE
ResEditor Editor<Parent, T>::residue(int i, const PropertyMap &map)
{
    return ResEditor( MoleculeView::residue(i, map) );
}

template<class Parent, class T>
SIRE_OUTOFLINE_TEMPLATE
ResEditor Editor<Parent, T>::residue(const QString &name, const PropertyMap &map)
{
    return ResEditor( MoleculeView::residue(name, map) );
}

/** Return the residues from this view that match the ID 'resid'

    \throw SireMol::missing_residue
    \throw SireError::invalid_index
    \throw SireMol::duplicate_residue
*/
template<class Parent, class T>
SIRE_OUTOFLINE_TEMPLATE
ResEditor Editor<Parent, T>::residue(const ResID &resid, const PropertyMap &map)
{
    return ResEditor( MoleculeView::residue(resid, map) );
}

/** Return the chain involved with this view - for this
    to work, only a single chain should be involved in this view

    \throw SireMol::missing_chain
    \throw SireMol::duplicate_chain
*/
template<class Parent, class T>
SIRE_OUTOFLINE_TEMPLATE
ChainEditor Editor<Parent, T>::chain()
{
    return ChainEditor( MoleculeView::chain() );
}

template<class Parent, class T>
SIRE_OUTOFLINE_TEMPLATE
ChainEditor Editor<Parent, T>::chain(int i, const PropertyMap &map)
{
    return ChainEditor( MoleculeView::chain(i, map) );
}

template<class Parent, class T>
SIRE_OUTOFLINE_TEMPLATE
ChainEditor Editor<Parent, T>::chain(const QString &name, const PropertyMap &map)
{
    return ChainEditor( MoleculeView::chain(name, map) );
}

/** Return the chains from this view that match the ID 'chainid'

    \throw SireMol::missing_chain
    \throw SireError::invalid_index
    \throw SireMol::duplicate_chain
*/
template<class Parent, class T>
SIRE_OUTOFLINE_TEMPLATE
ChainEditor Editor<Parent, T>::chain(const ChainID &chainid, const PropertyMap &map)
{
    return ChainEditor( MoleculeView::chain(chainid, map) );
}

/** Return the segment involved with this view - for this
    to work, only a single segment should be involved in this view

    \throw SireMol::missing_segment
    \throw SireMol::duplicate_segment
*/
template<class Parent, class T>
SIRE_OUTOFLINE_TEMPLATE
SegEditor Editor<Parent, T>::segment()
{
    return SegEditor( MoleculeView::segment() );
}

template<class Parent, class T>
SIRE_OUTOFLINE_TEMPLATE
SegEditor Editor<Parent, T>::segment(int i, const PropertyMap &map)
{
    return SegEditor( MoleculeView::segment(i, map) );
}

template<class Parent, class T>
SIRE_OUTOFLINE_TEMPLATE
SegEditor Editor<Parent, T>::segment(const QString &name, const PropertyMap &map)
{
    return SegEditor( MoleculeView::segment(name, map) );
}

/** Return the segments from this view that match the ID 'segid'

    \throw SireMol::missing_segment
    \throw SireError::invalid_index
    \throw SireMol::duplicate_segment
*/
template<class Parent, class T>
SIRE_OUTOFLINE_TEMPLATE
SegEditor Editor<Parent, T>::segment(const SegID &segid, const PropertyMap &map)
{
    return SegEditor( MoleculeView::segment(segid, map) );
}

/** Return the editor for the molecule that is viewed */
template<class Parent, class T>
SIRE_OUTOFLINE_TEMPLATE
MolEditor Editor<Parent, T>::molecule()
{
    return MolEditor( MoleculeView::molecule() );
}

/** Return the atom from this view that matches the ID 'atomid'

    \throw SireMol::missing_atom
    \throw SireError::invalid_index
    \throw SireMol::duplicate_atom
*/
template<class Parent, class T>
SIRE_OUTOFLINE_TEMPLATE
AtomEditor Editor<Parent, T>::select(const AtomID &atomid, const PropertyMap &map)
{
    return Editor<Parent,T>::atom(atomid, map);
}

/** Return the CutGroups from this view that match the ID 'cgid'

    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
    \throw SireMol::duplicate_cutgroup
*/
template<class Parent, class T>
SIRE_OUTOFLINE_TEMPLATE
CGEditor Editor<Parent, T>::select(const CGID &cgid, const PropertyMap &map)
{
    return Editor<Parent,T>::cutGroup(cgid, map);
}

/** Return the residues from this view that match the ID 'resid'

    \throw SireMol::missing_residue
    \throw SireError::invalid_index
    \throw SireMol::duplicate_residue
*/
template<class Parent, class T>
SIRE_OUTOFLINE_TEMPLATE
ResEditor Editor<Parent, T>::select(const ResID &resid, const PropertyMap &map)
{
    return Editor<Parent,T>::residue(resid, map);
}

/** Return the chains from this view that match the ID 'chainid'

    \throw SireMol::missing_chain
    \throw SireError::invalid_index
    \throw SireMol::duplicate_chain
*/
template<class Parent, class T>
SIRE_OUTOFLINE_TEMPLATE
ChainEditor Editor<Parent, T>::select(const ChainID &chainid, const PropertyMap &map)
{
    return Editor<Parent,T>::chain(chainid, map);
}

/** Return the segments from this view that match the ID 'segid'

    \throw SireMol::missing_segment
    \throw SireError::invalid_index
    \throw SireMol::duplicate_segment
*/
template<class Parent, class T>
SIRE_OUTOFLINE_TEMPLATE
SegEditor Editor<Parent, T>::select(const SegID &segid, const PropertyMap &map)
{
    return Editor<Parent,T>::segment(segid, map);
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

SIRE_END_HEADER

#endif
