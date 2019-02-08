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

#ifndef SIREMOL_SEGEDITOR_H
#define SIREMOL_SEGEDITOR_H

#include "structureeditor.h"
#include "segment.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class SegEditor;
class SegStructureEditor;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::SegEditor&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::SegEditor&);

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::SegStructureEditor&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::SegStructureEditor&);

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

class SegEditor;
typedef Editor<SegEditor, Segment> SegEditorBase;

/** This class is used to edit the non-structural parts of a segment

    @author Christopher Woods
*/
class SIREMOL_EXPORT SegEditor 
        : public SireBase::ConcreteProperty< SegEditor,Editor<SegEditor,Segment> >
{

friend SIREMOL_EXPORT QDataStream& ::operator<<(QDataStream&, const SegEditor&);
friend SIREMOL_EXPORT QDataStream& ::operator>>(QDataStream&, SegEditor&);

public:
    SegEditor();
    
    SegEditor(const Segment &residue);
    
    SegEditor(const SegEditor &other);
    
    ~SegEditor();
    
    SegEditor& operator=(const Segment &residue);
    SegEditor& operator=(const SegEditor &other);
    
    static const char* typeName();

    QString toString() const;
    
    SegEditor& rename(const SegName &name);
    SegStructureEditor reindex(SegIdx index) const;
    
    MolStructureEditor remove() const;

    AtomStructureEditor add(const AtomName &atomname) const;
    AtomStructureEditor add(AtomNum atomnum) const;
    
    SegStructureEditor remove(const AtomID &atomid) const;

    SegStructureEditor remove(int i) const;
    
    SegStructureEditor transfer(const AtomID &atomid, const SegID &segid) const;
    SegStructureEditor transfer(int i, const SegID &segid) const;
    
    SegStructureEditor transferAll(const SegID &segid) const;
    
    Segment commit() const;
};

/** This is the class used to edit a segment's structure 

    @author Christopher Woods
*/
class SIREMOL_EXPORT SegStructureEditor : public StructureEditor
{

friend SIREMOL_EXPORT QDataStream& ::operator<<(QDataStream&, const SegStructureEditor&);
friend SIREMOL_EXPORT QDataStream& ::operator>>(QDataStream&, SegStructureEditor&);

public:
    SegStructureEditor();
    SegStructureEditor(const Segment &residue);
    SegStructureEditor(const StructureEditor &data, SegIdx residx);
    
    SegStructureEditor(const SegStructureEditor &other);
    
    ~SegStructureEditor();
    
    static const char* typeName();
    
    const char* what() const
    {
        return SegStructureEditor::typeName();
    }
    
    SegStructureEditor* clone() const;
    
    SegStructureEditor& operator=(const Segment &residue);
    SegStructureEditor& operator=(const SegStructureEditor &other);
    
    QString toString() const;
    
    bool selectedAll() const;
    
    const SegName &name() const;
    SegIdx index() const;
    
    int nAtoms() const;
    
    MolStructureEditor molecule();
    
    AtomStructureEditor atom(int i);
    AtomStructureEditor atom(const AtomID &atomid);

    AtomStructureEditor select(int i);
    AtomStructureEditor select(const AtomID &atomid);
    
    SegStructureEditor& rename(const SegName &name);
    SegStructureEditor& reindex(SegIdx index);
    
    MolStructureEditor remove();

    AtomStructureEditor add(const AtomName &atomname);
    AtomStructureEditor add(AtomNum atomnum);
    
    SegStructureEditor& remove(const AtomID &atomid);

    SegStructureEditor& remove(int i);
    
    SegStructureEditor& transfer(const AtomID &atomid, const SegID &segid);
    SegStructureEditor& transfer(int i, const SegID &segid);
    
    SegStructureEditor& transferAll(const SegID &segid);

    template<class T>
    T property(const QString &key) const;

    template<class T>
    T metadata(const QString &metakey) const;

    template<class T>
    T metadata(const QString &key, const QString &metakey) const;

    template<class T>
    SegStructureEditor& setProperty(const QString &key, const T &value);

    template<class T>
    SegStructureEditor& setMetadata(const QString &metakey, const T &value);
    
    template<class T>
    SegStructureEditor& setMetadata(const QString &key, const QString &metakey,
                                    const T &value);
    
    Segment commit() const;
    
    operator Segment() const;

private:
    /** The unique ID for this segment in the molecule editor */
    quint32 uid;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Return the value for this segment of the property at key 'key'.
    Note that this property *must* be of type SegProperty<T> for
    this to work!
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
template<class T>
T SegStructureEditor::property(const QString &key) const
{
    const QVariant &value = this->getSegProperty(uid, key);
    return this->_pvt_getProperty<T>(key, value);
}

/** Return the value for this segment of the metadata at metakey 'metakey'.
    Note that this property *must* be of type SegProperty<T> for
    this to work!
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
template<class T>
T SegStructureEditor::metadata(const QString &metakey) const
{
    const QVariant &value = this->getSegMetadata(uid, metakey);
    return this->_pvt_getMetadata<T>(metakey, value);
}

/** Return the value for this segment of the metadata at metakey 'metakey'
    for the property at key 'key'.
    
    Note that this property *must* be of type SegProperty<T> for
    this to work!
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
template<class T>
T SegStructureEditor::metadata(const QString &key, 
                               const QString &metakey) const
{
    const QVariant &value = this->getSegMetadata(uid, key, metakey);
    return this->_pvt_getMetadata<T>(key, metakey, value);
}

/** Set the property at key 'key' to have the value 'value' for this segment.
    Note that an exception will be thrown if an existing property for
    this key is not of type SegProperty<T>

    \throw SireError::invalid_cast
*/
template<class T>
SegStructureEditor& SegStructureEditor::setProperty(const QString &key, 
                                                    const T &value)
{
    this->assertValidSegment(uid);

    //create space for this property
    this->_pvt_createSpaceForProperty< SegProperty<T> >(key);
    
    //now set the value of the property
    this->_pvt_setSegProperty(uid, key, QVariant::fromValue<T>(value));
    
    return *this;
}

/** Set the metadata at metakey 'metakey' to have the value 'value' for
    this segment. Note that an exception will be thrown if an existing 
    property for this metakey is not of type SegProperty<T>
    
    \throw SireError::invalid_cast
*/
template<class T>
SegStructureEditor& SegStructureEditor::setMetadata(const QString &metakey,
                                                    const T &value)
{
    this->assertValidSegment(uid);
    
    //create space for this metadata
    this->_pvt_createSpaceForMetadata< SegProperty<T> >(metakey);
    
    //now set the value of this metadata
    this->_pvt_setAtomMetadata(uid, metakey, QVariant::fromValue<T>(value));
    
    return *this;
}

/** Set the metadata at metakey 'metakey' for the property at
    key 'key' to have the value 'value' for
    this segment. Note that an exception will be thrown if an existing 
    property for this metakey is not of type SegProperty<T>
    
    \throw SireError::invalid_cast
*/
template<class T>
SegStructureEditor& SegStructureEditor::setMetadata(const QString &key, 
                                                    const QString &metakey,
                                                    const T &value)
{
    this->assertValidSegment(uid);
    
    //create space for this metadata
    this->_pvt_createSpaceForMetadata< SegProperty<T> >(key, metakey);
    
    //now set the value of this metadata
    this->_pvt_setAtomMetadata(uid, key, metakey, QVariant::fromValue<T>(value));
    
    return *this;
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

Q_DECLARE_METATYPE( SireMol::SegEditor );
Q_DECLARE_METATYPE( SireMol::SegStructureEditor );

SIRE_EXPOSE_CLASS( SireMol::SegEditor )
SIRE_EXPOSE_CLASS( SireMol::SegStructureEditor )

SIRE_EXPOSE_ALIAS( (SireMol::Editor<SireMol::SegEditor, SireMol::Segment>),
                    SireMol::SegEditorBase )

SIRE_END_HEADER

#endif
