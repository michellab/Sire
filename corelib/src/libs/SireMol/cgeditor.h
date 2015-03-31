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

#ifndef SIREMOL_CGEDITOR_H
#define SIREMOL_CGEDITOR_H

#include "structureeditor.h"
#include "cutgroup.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class CGEditor;
class CGStructureEditor;
}

QDataStream& operator<<(QDataStream&, const SireMol::CGEditor&);
QDataStream& operator>>(QDataStream&, SireMol::CGEditor&);

QDataStream& operator<<(QDataStream&, const SireMol::CGStructureEditor&);
QDataStream& operator>>(QDataStream&, SireMol::CGStructureEditor&);

namespace std
{
class slice;
}

namespace SireMol
{

class MolStructureEditor;
class SegStructureEditor;
class ChainStructureEditor;
class ResStructureEditor;
class CGStructureEditor;
class AtomStructureEditor;

class MolEditor;
class SegEditor;
class ChainEditor;
class ResEditor;
class CGEditor;
class AtomEditor;

class CGIdx;
class CGID;
class ResIdx;
class ResID;

class CGEditor;
typedef Editor<CGEditor, CutGroup> CGEditorBase;

/** This class is used to edit the non-structural parts of a CutGroup

    @author Christopher Woods
*/
class SIREMOL_EXPORT CGEditor 
        : public SireBase::ConcreteProperty< CGEditor, Editor<CGEditor, CutGroup> >
{

friend QDataStream& ::operator<<(QDataStream&, const CGEditor&);
friend QDataStream& ::operator>>(QDataStream&, CGEditor&);

public:
    CGEditor();
    
    CGEditor(const CutGroup &cutgroup);
    
    CGEditor(const CGEditor &other);
    
    ~CGEditor();
    
    CGEditor& operator=(const CutGroup &cutgroup);
    CGEditor& operator=(const CGEditor &other);
    
    static const char* typeName();
    
    QString toString() const;
    
    CGEditor& rename(const CGName &name);
    CGStructureEditor reindex(CGIdx index) const;
    
    MolStructureEditor remove() const;

    AtomStructureEditor add(const AtomName &atomname) const;
    AtomStructureEditor add(AtomNum atomnum) const;
    
    CGStructureEditor remove(const AtomID &atomid) const;

    CGStructureEditor remove(int i) const;
    
    CGStructureEditor transfer(const AtomID &atomid, const CGID &cgid) const;
    CGStructureEditor transfer(int i, const CGID &cgid) const;
    
    CGStructureEditor transferAll(const CGID &cgid) const;
    
    CutGroup commit() const;
};

/** This is the class used to edit a CutGroup's structure 

    @author Christopher Woods
*/
class SIREMOL_EXPORT CGStructureEditor : public StructureEditor
{

friend QDataStream& ::operator<<(QDataStream&, const CGStructureEditor&);
friend QDataStream& ::operator>>(QDataStream&, CGStructureEditor&);

public:
    CGStructureEditor();
    CGStructureEditor(const CutGroup &cutgroup);
    CGStructureEditor(const StructureEditor &data, CGIdx cgidx);
    
    CGStructureEditor(const CGStructureEditor &other);
    
    ~CGStructureEditor();
    
    static const char* typeName();
    
    const char* what() const
    {
        return CGStructureEditor::typeName();
    }
    
    CGStructureEditor* clone() const;
    
    CGStructureEditor& operator=(const CutGroup &cutgroup);
    CGStructureEditor& operator=(const CGStructureEditor &other);
    
    bool selectedAll() const;
    
    QString toString() const;
    
    const CGName& name() const;
    CGIdx index() const;
    
    int nAtoms() const;
    
    MolStructureEditor molecule();
    
    AtomStructureEditor atom(int i);
    AtomStructureEditor atom(const AtomID &atomid);

    AtomStructureEditor select(int i);
    AtomStructureEditor select(const AtomID &atomid);
    
    CGStructureEditor& rename(const CGName &name);
    CGStructureEditor& reindex(CGIdx index);
    
    MolStructureEditor remove();

    AtomStructureEditor add(const AtomName &atomname);
    AtomStructureEditor add(AtomNum atomnum);
    
    CGStructureEditor& remove(const AtomID &atomid);

    CGStructureEditor& remove(int i);
    
    CGStructureEditor& transfer(const AtomID &atomid, const CGID &cgid);
    CGStructureEditor& transfer(int i, const CGID &cgid);
    
    CGStructureEditor& transferAll(const CGID &cgid);

    template<class T>
    T property(const QString &key) const;

    template<class T>
    T metadata(const QString &metakey) const;

    template<class T>
    T metadata(const QString &key, const QString &metakey) const;

    template<class T>
    CGStructureEditor& setProperty(const QString &key, const T &value);

    template<class T>
    CGStructureEditor& setMetadata(const QString &metakey, const T &value);
    
    template<class T>
    CGStructureEditor& setMetadata(const QString &key, const QString &metakey,
                                   const T &value);
    
    CutGroup commit() const;
    
    operator CutGroup() const;

private:
    /** Unique ID number for this CutGroup in the molecule editor */
    quint32 uid;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Return the value for this CutGroup of the property at key 'key'.
    Note that this property *must* be of type CGProperty<T> for
    this to work!
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
template<class T>
T CGStructureEditor::property(const QString &key) const
{
    const QVariant &value = this->getCGProperty(uid, key);
    return this->_pvt_getProperty<T>(key, value);
}

/** Return the value for this CutGroup of the metadata at metakey 'metakey'.
    Note that this property *must* be of type CGProperty<T> for
    this to work!
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
template<class T>
T CGStructureEditor::metadata(const QString &metakey) const
{
    const QVariant &value = this->getCGMetadata(uid, metakey);
    return this->_pvt_getMetadata<T>(metakey, value);
}

/** Return the value for this CutGroup of the metadata at metakey 'metakey'
    for the property at key 'key'.
    
    Note that this property *must* be of type CGProperty<T> for
    this to work!
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
template<class T>
T CGStructureEditor::metadata(const QString &key, 
                              const QString &metakey) const
{
    const QVariant &value = this->getCGMetadata(uid, key, metakey);
    return this->_pvt_getMetadata<T>(key, metakey, value);
}

/** Set the property at key 'key' to have the value 'value' for this CutGroup.
    Note that an exception will be thrown if an existing property for
    this key is not of type CGProperty<T>

    \throw SireError::invalid_cast
*/
template<class T>
CGStructureEditor& CGStructureEditor::setProperty(const QString &key, 
                                                  const T &value)
{
    this->assertValidCutGroup(uid);

    //create space for this property
    this->_pvt_createSpaceForProperty< CGProperty<T> >(key);
    
    //now set the value of the property
    this->_pvt_setCGProperty(uid, key, QVariant::fromValue<T>(value));
    
    return *this;
}

/** Set the metadata at metakey 'metakey' to have the value 'value' for
    this CutGroup. Note that an exception will be thrown if an existing 
    property for this metakey is not of type CGProperty<T>
    
    \throw SireError::invalid_cast
*/
template<class T>
CGStructureEditor& CGStructureEditor::setMetadata(const QString &metakey,
                                                  const T &value)
{
    this->assertValidCutGroup(uid);
    
    //create space for this metadata
    this->_pvt_createSpaceForMetadata< CGProperty<T> >(metakey);
    
    //now set the value of this metadata
    this->_pvt_setAtomMetadata(uid, metakey, QVariant::fromValue<T>(value));
    
    return *this;
}

/** Set the metadata at metakey 'metakey' for the property at
    key 'key' to have the value 'value' for
    this CutGroup. Note that an exception will be thrown if an existing 
    property for this metakey is not of type CGProperty<T>
    
    \throw SireError::invalid_cast
*/
template<class T>
CGStructureEditor& CGStructureEditor::setMetadata(const QString &key, 
                                                  const QString &metakey,
                                                  const T &value)
{
    this->assertValidCutGroup(uid);
    
    //create space for this metadata
    this->_pvt_createSpaceForMetadata< CGProperty<T> >(key, metakey);
    
    //now set the value of this metadata
    this->_pvt_setAtomMetadata(uid, key, metakey, QVariant::fromValue<T>(value));
    
    return *this;
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

Q_DECLARE_METATYPE( SireMol::CGEditor );
Q_DECLARE_METATYPE( SireMol::CGStructureEditor );

SIRE_EXPOSE_ALIAS( (SireMol::Editor<SireMol::CGEditor, SireMol::CutGroup>),
                    SireMol::CGEditorBase )

SIRE_EXPOSE_CLASS( SireMol::CGEditor )
SIRE_EXPOSE_CLASS( SireMol::CGStructureEditor )

SIRE_END_HEADER

#endif
