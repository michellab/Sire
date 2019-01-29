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

#ifndef SIREMOL_CHAINEDITOR_H
#define SIREMOL_CHAINEDITOR_H

#include "structureeditor.h"
#include "chain.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class ChainEditor;
class ChainStructureEditor;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::ChainEditor&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::ChainEditor&);

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::ChainStructureEditor&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::ChainStructureEditor&);

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

class ChainEditor;
typedef Editor<ChainEditor, Chain> ChainEditorBase;

/** This class is used to edit the non-structural parts of a chain

    @author Christopher Woods
*/
class SIREMOL_EXPORT ChainEditor 
            : public SireBase::ConcreteProperty< ChainEditor,Editor<ChainEditor,Chain> >
{

friend QDataStream& ::operator<<(QDataStream&, const ChainEditor&);
friend QDataStream& ::operator>>(QDataStream&, ChainEditor&);

public:
    ChainEditor();
    
    ChainEditor(const Chain &chain);
    
    ChainEditor(const ChainEditor &other);
    
    ~ChainEditor();
    
    ChainEditor& operator=(const Chain &chain);
    ChainEditor& operator=(const ChainEditor &other);
    
    static const char* typeName();

    QString toString() const;
    
    ChainEditor& rename(const ChainName &name);
    ChainStructureEditor reindex(ChainIdx index) const;
    
    MolStructureEditor remove() const;

    ResStructureEditor add(const ResName &resname) const;
    ResStructureEditor add(ResNum atomnum) const;
    
    ChainStructureEditor remove(const AtomID &atomid) const;
    ChainStructureEditor remove(const ResID &resid) const;

    ChainStructureEditor remove(int i) const;
    
    ChainStructureEditor transfer(const ResID &resid, const ChainID &chainid) const;
    ChainStructureEditor transfer(int i, const ChainID &chainid) const;
    
    ChainStructureEditor transferAll(const ChainID &chainid) const;
    
    Chain commit() const;
};

/** This is the class used to edit a chain's structure 

    @author Christopher Woods
*/
class SIREMOL_EXPORT ChainStructureEditor : public StructureEditor
{

friend QDataStream& ::operator<<(QDataStream&, const ChainStructureEditor&);
friend QDataStream& ::operator>>(QDataStream&, ChainStructureEditor&);

public:
    ChainStructureEditor();
    ChainStructureEditor(const Chain &chain);
    ChainStructureEditor(const StructureEditor &data, ChainIdx chainidx);
    
    ChainStructureEditor(const ChainStructureEditor &other);
    
    ~ChainStructureEditor();
    
    static const char* typeName();
    
    const char* what() const
    {
        return ChainStructureEditor::typeName();
    }
    
    ChainStructureEditor* clone() const;
    
    QString toString() const;
    
    bool selectedAll() const;
    
    const ChainName& name() const;
    ChainIdx index() const;
    
    int nAtoms() const;
    int nResidues() const;
    
    ChainStructureEditor& operator=(const Chain &chain);
    ChainStructureEditor& operator=(const ChainStructureEditor &other);
    
    MolStructureEditor molecule();
    
    AtomStructureEditor atom(const AtomID &atomid);

    ResStructureEditor residue(int i);
    ResStructureEditor residue(const ResID &resid);

    AtomStructureEditor select(const AtomID &atomid);

    ResStructureEditor select(int i);
    ResStructureEditor select(const ResID &resid);
    
    ChainStructureEditor& rename(const ChainName &name);
    
    ChainStructureEditor& reindex(ChainIdx index);
    
    MolStructureEditor remove();

    ResStructureEditor add(const ResName &resname);
    ResStructureEditor add(ResNum resnum);
    
    ChainStructureEditor& remove(const ResID &resid);
    ChainStructureEditor& remove(const AtomID &atomid);

    ChainStructureEditor& remove(int i);
    
    ChainStructureEditor& transfer(const ResID &resid, const ChainID &chainid);
    ChainStructureEditor& transfer(int i, const ChainID &resid);
    
    ChainStructureEditor& transferAll(const ChainID &chainid);

    template<class T>
    T property(const QString &key) const;

    template<class T>
    T metadata(const QString &metakey) const;

    template<class T>
    T metadata(const QString &key, const QString &metakey) const;

    template<class T>
    ChainStructureEditor& setProperty(const QString &key, const T &value);

    template<class T>
    ChainStructureEditor& setMetadata(const QString &metakey, const T &value);
    
    template<class T>
    ChainStructureEditor& setMetadata(const QString &key, const QString &metakey,
                                      const T &value);
    
    Chain commit() const;
    
    operator Chain() const;

private:
    /** Unique ID to identify this chain in this molecule editor */
    quint32 uid;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Return the value for this chain of the property at key 'key'.
    Note that this property *must* be of type ChainProperty<T> for
    this to work!
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
template<class T>
T ChainStructureEditor::property(const QString &key) const
{
    const QVariant &value = this->getChainProperty(uid, key);
    return this->_pvt_getProperty<T>(key, value);
}

/** Return the value for this chain of the metadata at metakey 'metakey'.
    Note that this property *must* be of type ChainProperty<T> for
    this to work!
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
template<class T>
T ChainStructureEditor::metadata(const QString &metakey) const
{
    const QVariant &value = this->getChainMetadata(uid, metakey);
    return this->_pvt_getMetadata<T>(metakey, value);
}

/** Return the value for this chain of the metadata at metakey 'metakey'
    for the property at key 'key'.
    
    Note that this property *must* be of type ChainProperty<T> for
    this to work!
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
*/
template<class T>
T ChainStructureEditor::metadata(const QString &key, 
                                 const QString &metakey) const
{
    const QVariant &value = this->getChainMetadata(uid, key, metakey);
    return this->_pvt_getMetadata<T>(key, metakey, value);
}

/** Set the property at key 'key' to have the value 'value' for this chain.
    Note that an exception will be thrown if an existing property for
    this key is not of type ChainProperty<T>

    \throw SireError::invalid_cast
*/
template<class T>
ChainStructureEditor& ChainStructureEditor::setProperty(const QString &key, 
                                                  const T &value)
{
    this->assertValidChain(uid);

    //create space for this property
    this->_pvt_createSpaceForProperty< ChainProperty<T> >(key);
    
    //now set the value of the property
    this->_pvt_setChainProperty(uid, key, QVariant::fromValue<T>(value));
    
    return *this;
}

/** Set the metadata at metakey 'metakey' to have the value 'value' for
    this chain. Note that an exception will be thrown if an existing 
    property for this metakey is not of type ChainProperty<T>
    
    \throw SireError::invalid_cast
*/
template<class T>
ChainStructureEditor& ChainStructureEditor::setMetadata(const QString &metakey,
                                                        const T &value)
{
    this->assertValidChain(uid);
    
    //create space for this metadata
    this->_pvt_createSpaceForMetadata< ChainProperty<T> >(metakey);
    
    //now set the value of this metadata
    this->_pvt_setAtomMetadata(uid, metakey, QVariant::fromValue<T>(value));
    
    return *this;
}

/** Set the metadata at metakey 'metakey' for the property at
    key 'key' to have the value 'value' for
    this chain. Note that an exception will be thrown if an existing 
    property for this metakey is not of type ChainProperty<T>
    
    \throw SireError::invalid_cast
*/
template<class T>
ChainStructureEditor& ChainStructureEditor::setMetadata(const QString &key, 
                                                        const QString &metakey,
                                                        const T &value)
{
    this->assertValidChain(uid);
    
    //create space for this metadata
    this->_pvt_createSpaceForMetadata< ChainProperty<T> >(key, metakey);
    
    //now set the value of this metadata
    this->_pvt_setAtomMetadata(uid, key, metakey, QVariant::fromValue<T>(value));
    
    return *this;
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

Q_DECLARE_METATYPE( SireMol::ChainEditor );
Q_DECLARE_METATYPE( SireMol::ChainStructureEditor );

SIRE_EXPOSE_ALIAS( (SireMol::Editor<SireMol::ChainEditor, SireMol::Chain>),
                    SireMol::ChainEditorBase )

SIRE_EXPOSE_CLASS( SireMol::ChainEditor )
SIRE_EXPOSE_CLASS( SireMol::ChainStructureEditor )

SIRE_END_HEADER

#endif
