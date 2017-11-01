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

#include "atomeditor.h"
#include "atom.h"

#include "molecule.h"
#include "segment.h"
#include "chain.h"
#include "residue.h"
#include "cutgroup.h"

#include "moleditor.h"
#include "segeditor.h"
#include "chaineditor.h"
#include "reseditor.h"
#include "cgeditor.h"

#include "selector.hpp"
#include "mover.hpp"

#include "residx.h"
#include "cgidx.h"
#include "segidx.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMol;
using namespace SireStream;

namespace SireMol
{
    //fully instantiate the Editor<Atom> and Editor< Selector<Atom> > classes
    template class Editor<AtomEditor, Atom>;
    //template class Editor< Selector<Atom> >;
}

/////////
///////// Implementation of AtomEditor
/////////

static const RegisterMetaType<AtomEditor> r_atomeditor;

/** Serialise to a binary datastream */
QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds,
                                       const AtomEditor &atomeditor)
{
    writeHeader(ds, r_atomeditor, 1);
    
    ds << static_cast<const Editor<AtomEditor,Atom>&>(atomeditor);
    
    return ds;
}

/** Deserialise from a binary datastream */
QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds,
                                       AtomEditor &atomeditor)
{
    VersionID v = readHeader(ds, r_atomeditor);
    
    if (v == 1)
    {
        ds >> static_cast<Editor<AtomEditor,Atom>&>(atomeditor);
    }
    else
        throw version_error(v, "1", r_atomeditor, CODELOC);

    return ds;
}

/** Null constructor */
AtomEditor::AtomEditor() : ConcreteProperty< AtomEditor,Editor<AtomEditor,Atom> >()
{}

/** Construct an editor that edits a copy of 'atom' */
AtomEditor::AtomEditor(const Atom &atom)
           : ConcreteProperty< AtomEditor,Editor<AtomEditor,Atom> >(atom)
{}

/** Copy constructor */
AtomEditor::AtomEditor(const AtomEditor &other)
           : ConcreteProperty< AtomEditor,Editor<AtomEditor,Atom> >(other)
{}

/** Destructor */
AtomEditor::~AtomEditor()
{}

/** Set this editor so that it is editing a copy of 'atom' */
AtomEditor& AtomEditor::operator=(const Atom &atom)
{
    Editor<AtomEditor,Atom>::operator=(atom);
    return *this;
}

/** Copy assignment operator */
AtomEditor& AtomEditor::operator=(const AtomEditor &other)
{
    Editor<AtomEditor,Atom>::operator=(other);
    return *this;
}

/** Return a string representation of this editor */
QString AtomEditor::toString() const
{
    return QObject::tr( "Editor{ %1 }" ).arg( Atom::toString() );
}

/** Rename this atom so that it is called 'newname' */
AtomEditor& AtomEditor::rename(const AtomName &newname)
{
    if (newname == this->name())
        //nothing needs to be done
        return *this;
        
    d->rename( this->index(), newname );
    
    return *this;
}

/** Renumber this atom so that it has number 'newnum' */
AtomEditor& AtomEditor::renumber(AtomNum newnum)
{
    if (newnum == this->number())
        //nothing needs to be done
        return *this;
        
    d->renumber( this->index(), newnum );
    
    return *this;
}

/** Reindex this atom so that it lies at index 'newidx'. Note
    that if 'newidx' is greater than the number of atoms, then
    this will move this atom to be the last in the list */
AtomStructureEditor AtomEditor::reindex(AtomIdx newidx) const
{
    AtomStructureEditor editor(*this);
    editor.reindex(newidx);
    return editor;
}

/** Remove this atom from the molecule, returning an editor
    that can further edit the structure of the molecule */
MolStructureEditor AtomEditor::remove() const
{
    MolStructureEditor moleditor(*this);    
    moleditor.remove(this->index());
    
    return moleditor;
}

/** Reparent this atom so that it will be placed into the CutGroup
    with index 'cgidx' - this returns the updated atom in 
    an AtomStructureEditor, which is optimised for further
    editing of the molecule structure
    
    \throw SireError::invalid_index
*/
AtomStructureEditor AtomEditor::reparent(CGIdx cgidx) const
{
    AtomStructureEditor editor(*this);
    editor.reparent( cgidx );
    return editor;
}

/** Reparent this atom so that it will be placed into the CutGroup
    with ID 'cgid' - this returns the updated atom in 
    an AtomStructureEditor, which is optimised for further
    editing of the molecule structure
    
    \throw SireMol::missing_cutgroup
    \throw SireMol::duplicate_cutgroup
    \throw SireError::invalid_index
*/
AtomStructureEditor AtomEditor::reparent(const CGID &cgid) const
{
    AtomStructureEditor editor(*this);
    editor.reparent(cgid);
    return editor;
}

/** Reparent this atom so that it will be placed into the residue
    with index 'residx' - this returns the updated atom in 
    an AtomStructureEditor, which is optimised for further
    editing of the molecule structure
    
    \throw SireError::invalid_index
*/
AtomStructureEditor AtomEditor::reparent(ResIdx residx) const
{
    AtomStructureEditor editor(*this);
    editor.reparent(residx);
    return editor;
}

/** Reparent this atom so that it will be placed into the residue
    with ID 'resid' - this returns the updated atom in 
    an AtomStructureEditor, which is optimised for further
    editing of the molecule structure
    
    \throw SireMol::missing_residue
    \throw SireMol::duplicate_residue
    \throw SireError::invalid_index
*/
AtomStructureEditor AtomEditor::reparent(const ResID &resid) const
{
    AtomStructureEditor editor(*this);
    editor.reparent(resid);
    return editor;
}

/** Reparent this atom so that it will be placed into the segment
    with index 'segidx' - this returns the updated atom in 
    an AtomStructureEditor, which is optimised for further
    editing of the molecule structure
    
    \throw SireError::invalid_index
*/
AtomStructureEditor AtomEditor::reparent(SegIdx segidx) const
{
    AtomStructureEditor editor(*this);
    editor.reparent(segidx);
    return editor;
}

/** Reparent this atom so that it will be placed into the segment
    with ID 'segid' - this returns the updated atom in 
    an AtomStructureEditor, which is optimised for further
    editing of the molecule structure
    
    \throw SireMol::missing_segment
    \throw SireMol::duplicate_segment
    \throw SireError::invalid_index
*/
AtomStructureEditor AtomEditor::reparent(const SegID &segid) const
{
    AtomStructureEditor editor(*this);
    editor.reparent(segid);
    return editor;
}

const char* AtomEditor::typeName()
{
    return QMetaType::typeName( qMetaTypeId<AtomEditor>() );
}

/////////
///////// Implementation of AtomStructureEditor
/////////

static const RegisterMetaType<AtomStructureEditor> r_atomstructeditor;

/** Serialise to a binary datastream */
QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds,
                                       const AtomStructureEditor &atomeditor)
{
    writeHeader(ds, r_atomstructeditor, 1);
    
    ds << atomeditor.uid
       << static_cast<const StructureEditor&>(atomeditor);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds,
                                       AtomStructureEditor &atomeditor)
{
    VersionID v = readHeader(ds, r_atomstructeditor);
    
    if (v == 1)
    {
        ds >> atomeditor.uid
           >> static_cast<StructureEditor&>(atomeditor);
    }
    else
        throw version_error( v, "1", r_atomstructeditor, CODELOC );

    return ds;
}

/** Null constructor */
AtomStructureEditor::AtomStructureEditor() : StructureEditor()
{
    this->operator=( StructureEditor::addAtom() );
}

/** Construct from an Atom */
AtomStructureEditor::AtomStructureEditor(const Atom &atom)
                    : StructureEditor(atom.data())
{
    uid = this->getUID(atom.index());
}

/** Construct for the atom at index 'idx' in the molecule whose data
    is being edited in 'moldata'
    
    \throw SireError::invalid_index
*/
AtomStructureEditor::AtomStructureEditor(const StructureEditor &moldata,
                                         AtomIdx idx)
                    : StructureEditor(moldata)
{
    uid = this->getUID(idx);
}

/** Copy constructor */
AtomStructureEditor::AtomStructureEditor(const AtomStructureEditor &other)
                    : StructureEditor(other), uid(other.uid)
{}

/** Destructor */
AtomStructureEditor::~AtomStructureEditor()
{}

/** Assign to edit the structure of a copy of 'atom' */
AtomStructureEditor& AtomStructureEditor::operator=(const Atom &atom)
{
    StructureEditor::operator=(atom.data());
    uid = this->getUID(atom.index());
    
    return *this;
}

/** Copy assignment operator */
AtomStructureEditor& AtomStructureEditor::operator=(const AtomStructureEditor &other)
{
    StructureEditor::operator=(other);
    uid = other.uid;
    return *this;
}

/** Return a string representation of this editor */
QString AtomStructureEditor::toString() const
{
    return QObject::tr( "StructureEditor{ Atom( %1 : %2 ) }" )
                    .arg( this->name() )
                    .arg( this->number() );
}

/** Return whether or not this contains the whole molecule */
bool AtomStructureEditor::selectedAll() const
{
    return StructureEditor::nAtomsInMolecule() == 1;
}

/** Return the name of this atom */
const AtomName& AtomStructureEditor::name() const
{
    return this->atomName(uid);
}

/** Return the number of this atom */
AtomNum AtomStructureEditor::number() const
{
    return this->atomNum(uid);
}

/** Return the index number of this atom in the molecule */
AtomIdx AtomStructureEditor::index() const
{
    return this->atomIdx(uid);
}

/** Return the editor for the residue that contains this atom */
ResStructureEditor AtomStructureEditor::residue()
{
    return ResStructureEditor(*this, resIdx(residueParentOfAtom(uid)));
}

/** Return the editor for the CutGrop that contains this atom */
CGStructureEditor AtomStructureEditor::cutGroup()
{
    return CGStructureEditor(*this, cgIdx(cutGroupParentOfAtom(uid)));
}

/** Return the editor for the chain that contains this atom */
ChainStructureEditor AtomStructureEditor::chain()
{
    return ChainStructureEditor(*this, chainIdx(chainParentOfAtom(uid)));
}

/** Return the editor for the segment that contain this atom */
SegStructureEditor AtomStructureEditor::segment()
{
    return SegStructureEditor(*this, segIdx(segmentParentOfAtom(uid)));
}

/** Return the editor for the molecule that contains this atom */
MolStructureEditor AtomStructureEditor::molecule()
{
    return MolStructureEditor(*this);
}

/** Rename this atom to 'newname' */
AtomStructureEditor& AtomStructureEditor::rename(const AtomName &newname)
{
    this->renameAtom(uid, newname);
    return *this;
}

/** Renumber this atom to 'newnum' */
AtomStructureEditor& AtomStructureEditor::renumber(AtomNum newnum)
{
    this->renumberAtom(uid, newnum);
    return *this;
}

/** Reindex this atom to 'newidx' - this will move the atom to 
    the end if 'newidx' is greater than the number of atoms
    in the molecule */
AtomStructureEditor& AtomStructureEditor::reindex(AtomIdx newidx)
{
    this->reindexAtom(uid, newidx);
    return *this;
}

/** Completely remove this atom from the molecule and return
    a MolStructureEditor that can be used to continue editing
    the molecule */
MolStructureEditor AtomStructureEditor::remove()
{
    this->removeAtom(uid);
    return MolStructureEditor(*this);
}

/** Reparent this atom so that it is now in the CutGroup at index 'cgidx'

    \throw SireError::invalid_index
*/
AtomStructureEditor& AtomStructureEditor::reparent(CGIdx cgidx)
{
    this->reparentAtom(uid, cgidx);
    return *this;
}

/** Reparent this atom so that it is now in the CutGroup identified
    by ID 'cgid' 
    
    \throw SireMol::missing_cutgroup
    \throw SireMol::duplicate_cutgroup
    \throw SireMol::invalid_index
*/
AtomStructureEditor& AtomStructureEditor::reparent(const CGID &cgid)
{
    return this->reparent( this->cgIdx(cgid) );
}

/** Reparent this atom so that it is now in the residue at index 'residx'

    \throw SireError::invalid_index
*/
AtomStructureEditor& AtomStructureEditor::reparent(ResIdx residx)
{
    this->reparentAtom(uid, residx);
    return *this;
}

/** Reparent this atom so that it is now in the residue identified
    by ID 'resid' 
    
    \throw SireMol::missing_residue
    \throw SireMol::duplicate_residue
    \throw SireMol::invalid_index
*/
AtomStructureEditor& AtomStructureEditor::reparent(const ResID &resid)
{
    return this->reparent( this->resIdx(resid) );
}

/** Reparent this atom so that it is now in the segment at index 'segidx'

    \throw SireError::invalid_index
*/
AtomStructureEditor& AtomStructureEditor::reparent(SegIdx segidx)
{
    this->reparentAtom(uid, segidx);
    return *this;
}

/** Reparent this atom so that it is now in the segment identified
    by ID 'segid' 
    
    \throw SireMol::missing_segment
    \throw SireMol::duplicate_segment
    \throw SireMol::invalid_index
*/
AtomStructureEditor& AtomStructureEditor::reparent(const SegID &segid)
{
    return this->reparent( this->segIdx(segid) );
}

/** Commit all of the changes, returning the uneditable Atom */
Atom AtomStructureEditor::commit() const
{
    return Atom( this->commitChanges(), this->index() );
}

/** Allow automatic casting to an Atom() */
AtomStructureEditor::operator Atom() const
{
    return this->commit();
}

const char* AtomStructureEditor::typeName()
{
    return QMetaType::typeName( qMetaTypeId<AtomStructureEditor>() );
}

AtomStructureEditor* AtomStructureEditor::clone() const
{
    return new AtomStructureEditor(*this);
}
