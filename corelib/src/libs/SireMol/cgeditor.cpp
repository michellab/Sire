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
#include "cgeditor.h"
#include "reseditor.h"
#include "chaineditor.h"
#include "segeditor.h"
#include "moleditor.h"

#include "atom.h"
#include "cutgroup.h"
#include "residue.h"
#include "chain.h"
#include "segment.h"
#include "molecule.h"
#include "mover.hpp"
#include "selector.hpp"

#include "groupatomids.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"

using namespace SireMol;
using namespace SireStream;

namespace SireMol
{
    // fully instantiate Editor<CutGroup>
    template class Editor<CGEditor,CutGroup>;
    //template class Selector< Editor<CutGroup> >;
}

////////
//////// Implementation of CGEditor
////////

static const RegisterMetaType<CGEditor> r_cgeditor;

/** Serialise to a binary datastream */
QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds,
                                       const CGEditor &cgeditor)
{
    writeHeader(ds, r_cgeditor, 1);
    
    ds << static_cast<const Editor<CGEditor,CutGroup>&>(cgeditor);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds,
                                       CGEditor &cgeditor)
{
    VersionID v = readHeader(ds, r_cgeditor);
    
    if (v == 1)
    {
        ds >> static_cast<Editor<CGEditor,CutGroup>&>(cgeditor);
    }
    else
        throw version_error( v, "1", r_cgeditor, CODELOC );
        
    return ds;
}

/** Null constructor */
CGEditor::CGEditor() : ConcreteProperty< CGEditor, Editor<CGEditor,CutGroup> >()
{}

/** Construct to edit a copy of the CutGroup 'cutgroup' */
CGEditor::CGEditor(const CutGroup &cutgroup) 
         : ConcreteProperty< CGEditor, Editor<CGEditor,CutGroup> >(cutgroup)
{}

/** Copy constructor */
CGEditor::CGEditor(const CGEditor &other) 
         : ConcreteProperty< CGEditor, Editor<CGEditor,CutGroup> >(other)
{}

/** Destructor */
CGEditor::~CGEditor()
{}

/** Assign this editor so that it edits a copy of the CutGroup 'cutgroup' */
CGEditor& CGEditor::operator=(const CutGroup &cutgroup)
{
    Editor<CGEditor,CutGroup>::operator=(cutgroup);
    return *this;
}

/** Copy assignment operator */
CGEditor& CGEditor::operator=(const CGEditor &other)
{
    Editor<CGEditor,CutGroup>::operator=(other);
    return *this;
}

/** Return a string representation of this editor */
QString CGEditor::toString() const
{
    return QObject::tr( "Editor{ %1 }" ).arg(CutGroup::toString());
}

/** Rename this CutGroup to 'newname' */
CGEditor& CGEditor::rename(const CGName &newname)
{
    if (newname == this->name())
        //nothing to do
        return *this;
        
    d->rename( this->index(), newname );
    
    return *this;
}

/** Move this CutGroup to index 'newidx' - this will move it
    to the start or end if this index is out of range */
CGStructureEditor CGEditor::reindex(CGIdx newidx) const
{
    CGStructureEditor editor(*this);
    editor.reindex(newidx);
    
    return editor;
}

/** Complete remove this CutGroup, and return an editor
    for the molecule that contained it */
MolStructureEditor CGEditor::remove() const
{
    CGStructureEditor editor(*this);
    
    return editor.remove();
}

/** Add an atom called 'atomname' to this CutGroup and return
    an editor for that atom */
AtomStructureEditor CGEditor::add(const AtomName &atomname) const
{
    CGStructureEditor editor(*this);
    
    return editor.add(atomname);
}

/** Add an atom with number 'atomnum' to this CutGroup and return
    an editor for that atom */
AtomStructureEditor CGEditor::add(AtomNum atomnum) const
{
    CGStructureEditor editor(*this);
    
    return editor.add(atomnum);
}

/** Completely remove all atoms that match the ID 'atomid' from 
    this CutGroup 
    
    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
CGStructureEditor CGEditor::remove(const AtomID &atomid) const
{
    CGStructureEditor editor(*this);
    editor.remove(atomid);
    
    return editor;
}

/** Remove the ith atom from this CutGroup

    \throw SireError::invalid_index
*/
CGStructureEditor CGEditor::remove(int i) const
{
    CGStructureEditor editor(*this);
    editor.remove(i);
    
    return editor;
}

/** Transfer all atoms that match the ID 'atomid' in this CutGroup 
    to the CutGroup that matches the ID 'cgid'
    
    \throw SireMol::missing_atom
    \throw SireMol::missing_cutgroup
    \throw SireMol::duplicate_cutgroup
    \throw SireError::invalid_index
*/
CGStructureEditor CGEditor::transfer(const AtomID &atomid, 
                                     const CGID &cgid) const
{
    CGStructureEditor editor(*this);
    editor.transfer(atomid, cgid);
    
    return editor;
}
                                     
/** Transfer the ith atom of this CutGroup into the CutGroup identified
    by the ID 'cgid'
    
    \throw SireError::invalid_index
*/
CGStructureEditor CGEditor::transfer(int i, const CGID &cgid) const
{
    CGStructureEditor editor(*this);
    editor.transfer(i, cgid);
    
    return editor;
}

/** Completely transfer all of the atoms in this CutGroup to 
    the CutGroup that matches the ID 'cgid'
    
    \throw SireMol::missing_cutgroup
    \throw SireMol::duplicate_cutgroup
    \throw SireError::invalid_index
*/
CGStructureEditor CGEditor::transferAll(const CGID &cgid) const
{
    CGStructureEditor editor(*this);
    editor.transferAll(cgid);
    
    return editor;
}

/** Commit the changes made by this editor and return the 
    updated CutGroup */
CutGroup CGEditor::commit() const
{
    return *this;
}

const char* CGEditor::typeName()
{
    return QMetaType::typeName( qMetaTypeId<CGEditor>() );
}

////////
//////// Implementation of CGStructureEditor
////////

static const RegisterMetaType<CGStructureEditor> r_cgstructeditor;

/** Serialise to a binary datastream */
QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds,
                                       const CGStructureEditor &cgeditor)
{
    writeHeader(ds, r_cgstructeditor, 1);
    
    ds << cgeditor.uid
       << static_cast<const StructureEditor&>(cgeditor);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds,
                                       CGStructureEditor &cgeditor)
{
    VersionID v = readHeader(ds, r_cgstructeditor);
    
    if (v == 1)
    {
        ds >> cgeditor.uid
           >> static_cast<StructureEditor&>(cgeditor);
    }
    else
        throw version_error( v, "1", r_cgstructeditor, CODELOC );
        
    return ds;
}

/** Null constructor */
CGStructureEditor::CGStructureEditor()
                  : StructureEditor(), uid(0)
{}

/** Construct to edit a copy of the CutGroup 'cutgroup' */
CGStructureEditor::CGStructureEditor(const CutGroup &cutgroup)
                  : StructureEditor(cutgroup.data())
{
    uid = this->getUID(cutgroup.index());
}

/** Construct to edit the CutGroup at index 'cgidx' in the molecule 
    also being edited in 'data'
    
    \throw SireError::invalid_index
*/
CGStructureEditor::CGStructureEditor(const StructureEditor &data, CGIdx cgidx)
                  : StructureEditor(data)
{
    uid = this->getUID(cgidx);
}

/** Copy constructor */
CGStructureEditor::CGStructureEditor(const CGStructureEditor &other)
                  : StructureEditor(other), uid(other.uid)
{}

/** Destructor */
CGStructureEditor::~CGStructureEditor()
{}

/** Assign this editor so that it edits a copy of the CutGroup 'cutgroup' */
CGStructureEditor& CGStructureEditor::operator=(const CutGroup &cutgroup)
{
    StructureEditor::operator=(cutgroup.data());
    uid = this->getUID(cutgroup.index());
    
    return *this;
}

/** Copy assignment operator */
CGStructureEditor& CGStructureEditor::operator=(const CGStructureEditor &other)
{
    StructureEditor::operator=(other);
    uid = other.uid;
    
    return *this;
}

/** Return a string representation of this editor */
QString CGStructureEditor::toString() const
{
    return QObject::tr( "StructureEditor{ CutGroup( %1 ) }" )
                .arg( this->name() );
}

/** Return the name of this CutGroup */
const CGName& CGStructureEditor::name() const
{
    return this->cgName(uid);
}

/** Does this hold the entire molecule */
bool CGStructureEditor::selectedAll() const
{
    return StructureEditor::nCutGroupsInMolecule() == 1;
}

/** Return the index of this CutGroup in the molecule */
CGIdx CGStructureEditor::index() const
{
    return this->cgIdx(uid);
}

/** Return the number of atoms in this CutGroup (could be zero!) */
int CGStructureEditor::nAtoms() const
{
    return this->nAtomsInCutGroup(uid);
}

/** Return an editor for the molecule that contains this CutGroup */
MolStructureEditor CGStructureEditor::molecule()
{
    return MolStructureEditor(*this);
}

/** Return an editor for the ith atom of this CutGroup

    \throw SireError::invalid_index
*/
AtomStructureEditor CGStructureEditor::atom(int i)
{
    return AtomStructureEditor(*this, atomIdx( atomInCutGroup(uid,i) ));
}

/** Return an editor for the atom that matches the ID 'atomid' in
    this CutGroup
    
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
AtomStructureEditor CGStructureEditor::atom(const AtomID &atomid)
{
    return AtomStructureEditor(*this, atomIdx(this->index() + atomid));
}

/** Return an editor for the ith atom of this CutGroup

    \throw SireError::invalid_index
*/
AtomStructureEditor CGStructureEditor::select(int i)
{
    return AtomStructureEditor(*this, atomIdx( atomInCutGroup(uid,i) ));
}

/** Return an editor for the atom that matches the ID 'atomid' in
    this CutGroup
    
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
AtomStructureEditor CGStructureEditor::select(const AtomID &atomid)
{
    return AtomStructureEditor(*this, atomIdx( this->index() + atomid ));
}

/** Rename this CutGroup to 'newname' */
CGStructureEditor& CGStructureEditor::rename(const CGName &newname)
{
    this->renameCutGroup(uid, newname);
    return *this;
}

/** Move this CutGroup to index 'newidx' - this will move it
    to the start or end if this index is out of range */
CGStructureEditor& CGStructureEditor::reindex(CGIdx newidx)
{
    this->reindexCutGroup(uid, newidx);
    return *this;
}

/** Complete remove this CutGroup, and return an editor
    for the molecule that contained it */
MolStructureEditor CGStructureEditor::remove()
{
    this->removeCutGroup(uid);
    return MolStructureEditor(*this);
}

/** Add an atom called 'atomname' to this CutGroup and return
    an editor for that atom */
AtomStructureEditor CGStructureEditor::add(const AtomName &atomname)
{
    this->assertValidCutGroup(uid);

    AtomStructureEditor atom = this->addAtom();
    atom.rename(atomname);
    atom.reparent( this->index() );
    
    return atom;
}

/** Add an atom with number 'atomnum' to this CutGroup and return
    an editor for that atom */
AtomStructureEditor CGStructureEditor::add(AtomNum atomnum)
{
    this->assertValidCutGroup(uid);

    AtomStructureEditor atom = this->addAtom();
    atom.renumber(atomnum);
    atom.reparent( this->index() );
    
    return atom;
}

/** Completely remove all atoms that match the ID 'atomid' from 
    this CutGroup 
    
    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
CGStructureEditor& CGStructureEditor::remove(const AtomID &atomid)
{
    this->removeAtoms( this->index() + atomid );
    return *this;
}

/** Remove the ith atom from this CutGroup

    \throw SireError::invalid_index
*/
CGStructureEditor& CGStructureEditor::remove(int i)
{
    this->removeAtoms( atomIdx(atomInCutGroup(uid,i)) );
    return *this;
}

/** Transfer all atoms that match the ID 'atomid' in this CutGroup 
    to the CutGroup that matches the ID 'cgid'
    
    \throw SireMol::missing_atom
    \throw SireMol::missing_cutgroup
    \throw SireMol::duplicate_cutgroup
    \throw SireError::invalid_index
*/
CGStructureEditor& CGStructureEditor::transfer(const AtomID &atomid, 
                                               const CGID &cgid)
{
    this->reparentAtom( this->getUID(this->index() + atomid),
                        this->cgIdx(cgid) );
                        
    return *this;
}

/** Transfer the ith atom of this CutGroup into the CutGroup identified
    by the ID 'cgid'
    
    \throw SireError::invalid_index
*/
CGStructureEditor& CGStructureEditor::transfer(int i, const CGID &cgid)
{
    this->reparentAtom( atomInCutGroup(uid,i), this->cgIdx(cgid) );
    return *this;
}

/** Completely transfer all of the atoms in this CutGroup to 
    the CutGroup that matches the ID 'cgid'
    
    \throw SireMol::missing_cutgroup
    \throw SireMol::duplicate_cutgroup
    \throw SireError::invalid_index
*/
CGStructureEditor& CGStructureEditor::transferAll(const CGID &cgid)
{
    CGIdx cgidx = this->cgIdx(cgid);
    
    int nats = this->nAtoms();
    
    for (int i=0; i<nats; ++i)
    {
        this->reparentAtom( atomInCutGroup(uid,i), cgidx );
    }
    
    return *this;
}

/** Commit the changes made by this editor and return the 
    updated CutGroup */
CutGroup CGStructureEditor::commit() const
{
    return CutGroup( this->commitChanges(), this->index() );
}

/** Allow automatic casting of this editor to a CutGroup */
CGStructureEditor::operator CutGroup() const
{
    return this->commit();
}

const char* CGStructureEditor::typeName()
{
    return QMetaType::typeName( qMetaTypeId<CGStructureEditor>() );
}

CGStructureEditor* CGStructureEditor::clone() const
{
    return new CGStructureEditor(*this);
}
