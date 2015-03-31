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

#include "editor.hpp"

#ifndef SIREMOL_ATOMEDITOR_H
#define SIREMOL_ATOMEDITOR_H

#include "structureeditor.h"
#include "atom.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class AtomEditor;
class AtomStructureEditor;
}

QDataStream &operator<<(QDataStream&, const SireMol::AtomEditor&);
QDataStream &operator>>(QDataStream&, SireMol::AtomEditor&);

QDataStream &operator<<(QDataStream&, const SireMol::AtomStructureEditor&);
QDataStream &operator>>(QDataStream&, SireMol::AtomStructureEditor&);

namespace SireMol
{

class MolStructureEditor;
class SegStructureEditor;
class ChainStructureEditor;
class ResStructureEditor;
class CGStructureEditor;

class MolEditor;
class SegEditor;
class ChainEditor;
class ResEditor;
class CGEditor;

class CGIdx;
class CGID;
class ResIdx;
class ResID;

class AtomEditor;
typedef Editor<AtomEditor, Atom> AtomEditorBase;

/** This class is used to edit an atom in a molecule. This
    class is able to edit everything about the molecule
    *except* for its relationship to other parts of the
    molecule. To do that, you need an AtomStructureEditor
    (which is created automatically by the 'reparent()'
    function)

    @author Christopher Woods
*/
class SIREMOL_EXPORT AtomEditor 
            : public SireBase::ConcreteProperty<AtomEditor,AtomEditorBase>
{

friend QDataStream& ::operator<<(QDataStream&, const AtomEditor&);
friend QDataStream& ::operator>>(QDataStream&, AtomEditor&);

public:
    AtomEditor();
    
    AtomEditor(const Atom &atom);
    
    AtomEditor(const AtomEditor &other);
    
    ~AtomEditor();
    
    AtomEditor& operator=(const Atom &other);
    AtomEditor& operator=(const AtomEditor &other);

    static const char* typeName();
    
    const char* what() const
    {
        return AtomEditor::typeName();
    }

    QString toString() const;

    AtomEditor& rename(const AtomName &name);
    AtomEditor& renumber(AtomNum number);
    
    AtomStructureEditor reindex(AtomIdx atomidx) const;

    MolStructureEditor remove() const;
    
    AtomStructureEditor reparent(CGIdx cgidx) const;
    AtomStructureEditor reparent(const CGID &cgid) const;
    
    AtomStructureEditor reparent(ResIdx residx) const;
    AtomStructureEditor reparent(const ResID &resid) const;

    AtomStructureEditor reparent(SegIdx segidx) const;
    AtomStructureEditor reparent(const SegID &segid) const;
};

/** This class is used to edit an atom's relationship to 
    other parts of the molecule (e.g. which CutGroup it
    is in, or which Residue it is in)
    
    @author Christopher Woods
*/
class SIREMOL_EXPORT AtomStructureEditor : public StructureEditor
{

friend QDataStream& ::operator<<(QDataStream&, const AtomStructureEditor&);
friend QDataStream& ::operator>>(QDataStream&, AtomStructureEditor&); 

public:
    AtomStructureEditor();
    AtomStructureEditor(const Atom &atom);
    AtomStructureEditor(const StructureEditor &data, AtomIdx atomidx);

    AtomStructureEditor(const AtomStructureEditor &other);
    
    ~AtomStructureEditor();

    static const char* typeName();
    
    const char* what() const
    {
        return AtomStructureEditor::typeName();
    }
    
    AtomStructureEditor* clone() const;

    AtomStructureEditor& operator=(const Atom &atom);
    AtomStructureEditor& operator=(const AtomStructureEditor &other);

    bool selectedAll() const;

    QString toString() const;

    ResStructureEditor residue();
    CGStructureEditor cutGroup();
    ChainStructureEditor chain();
    SegStructureEditor segment();
    MolStructureEditor molecule();

    const AtomName& name() const;
    AtomNum number() const;
    AtomIdx index() const;

    AtomStructureEditor& rename(const AtomName &name);
    AtomStructureEditor& renumber(AtomNum number);
    AtomStructureEditor& reindex(AtomIdx idx);
    
    MolStructureEditor remove();

    AtomStructureEditor& reparent(CGIdx cgidx);
    AtomStructureEditor& reparent(const CGID &cgid);
    
    AtomStructureEditor& reparent(ResIdx residx);
    AtomStructureEditor& reparent(const ResID &resid);
    
    AtomStructureEditor& reparent(SegIdx segidx);
    AtomStructureEditor& reparent(const SegID &segid);

    template<class T>
    T property(const QString &key) const;

    template<class T>
    T metadata(const QString &metakey) const;

    template<class T>
    T metadata(const QString &key, const QString &metakey) const;

    template<class T>
    AtomStructureEditor& setProperty(const QString &key, const T &value);

    template<class T>
    AtomStructureEditor& setMetadata(const QString &metakey, const T &value);
    
    template<class T>
    AtomStructureEditor& setMetadata(const QString &key, const QString &metakey,
                                     const T &value);

    Atom commit() const;

    operator Atom() const;

private:
    /** The unique (temporary) ID of this atom in the molecule.
        A temporary and private ID is used as this atom can have
        all of its other ID tokens changed, so it would be difficult
        to keep track of if it didn't have a private, non-editable ID */
    quint32 uid;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Return the value for this atom of the property at key 'key'.
    Note that this property *must* be of type AtomProperty<T> for
    this to work!
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
template<class T>
T AtomStructureEditor::property(const QString &key) const
{
    const QVariant &value = this->getAtomProperty(uid, key);
    return this->_pvt_getProperty<T>(key, value);
}

/** Return the value for this atom of the metadata at metakey 'metakey'.
    Note that this property *must* be of type AtomProperty<T> for
    this to work!
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
template<class T>
T AtomStructureEditor::metadata(const QString &metakey) const
{
    const QVariant &value = this->getAtomMetadata(uid, metakey);
    return this->_pvt_getMetadata<T>(metakey, value);
}

/** Return the value for this atom of the metadata at metakey 'metakey'
    for the property at key 'key'.
    
    Note that this property *must* be of type AtomProperty<T> for
    this to work!
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
template<class T>
T AtomStructureEditor::metadata(const QString &key, 
                                const QString &metakey) const
{
    const QVariant &value = this->getAtomMetadata(uid, key, metakey);
    return this->_pvt_getMetadata<T>(key, metakey, value);
}

/** Set the property at key 'key' to have the value 'value' for this atom.
    Note that an exception will be thrown if an existing property for
    this key is not of type AtomProperty<T>

    \throw SireError::invalid_cast
*/
template<class T>
AtomStructureEditor& AtomStructureEditor::setProperty(const QString &key, 
                                                      const T &value)
{
    this->assertValidAtom(uid);

    //create space for this property
    this->_pvt_createSpaceForProperty< AtomProperty<T> >(key);
    
    //now set the value of the property
    this->_pvt_setAtomProperty(uid, key, QVariant::fromValue<T>(value));
    
    return *this;
}

/** Set the metadata at metakey 'metakey' to have the value 'value' for
    this atom. Note that an exception will be thrown if an existing 
    property for this metakey is not of type AtomProperty<T>
    
    \throw SireError::invalid_cast
*/
template<class T>
AtomStructureEditor& AtomStructureEditor::setMetadata(const QString &metakey,
                                                      const T &value)
{
    this->assertValidAtom(uid);
    
    //create space for this metadata
    this->_pvt_createSpaceForMetadata< AtomProperty<T> >(metakey);
    
    //now set the value of this metadata
    this->_pvt_setAtomMetadata(uid, metakey, QVariant::fromValue<T>(value));
    
    return *this;
}

/** Set the metadata at metakey 'metakey' for the property at
    key 'key' to have the value 'value' for
    this atom. Note that an exception will be thrown if an existing 
    property for this metakey is not of type AtomProperty<T>
    
    \throw SireError::invalid_cast
*/
template<class T>
AtomStructureEditor& AtomStructureEditor::setMetadata(const QString &key, 
                                                      const QString &metakey,
                                                      const T &value)
{
    this->assertValidAtom(uid);
    
    //create space for this metadata
    this->_pvt_createSpaceForMetadata< AtomProperty<T> >(key, metakey);
    
    //now set the value of this metadata
    this->_pvt_setAtomMetadata(uid, key, metakey, QVariant::fromValue<T>(value));
    
    return *this;
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

Q_DECLARE_METATYPE( SireMol::AtomEditor );
Q_DECLARE_METATYPE( SireMol::AtomStructureEditor );

SIRE_EXPOSE_ALIAS( (SireMol::Editor<SireMol::AtomEditor, SireMol::Atom>),
                    SireMol::AtomEditorBase )

SIRE_EXPOSE_CLASS( SireMol::AtomEditor )
SIRE_EXPOSE_CLASS( SireMol::AtomStructureEditor )

SIRE_END_HEADER

#endif
