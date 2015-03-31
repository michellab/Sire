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
#include "segeditor.h"
#include "reseditor.h"
#include "chaineditor.h"
#include "segeditor.h"
#include "moleditor.h"

#include "atom.h"
#include "segment.h"
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
    // fully instantiate Editor<Segment>
    template class Editor<SegEditor,Segment>;
    //template class Selector< Editor<Segment> >;
}

////////
//////// Implementation of SegEditor
////////

static const RegisterMetaType<SegEditor> r_segeditor;

/** Serialise to a binary datastream */
QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds,
                                       const SegEditor &segeditor)
{
    writeHeader(ds, r_segeditor, 1);
    
    ds << static_cast<const Editor<SegEditor,Segment>&>(segeditor);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds,
                                       SegEditor &segeditor)
{
    VersionID v = readHeader(ds, r_segeditor);
    
    if (v == 1)
    {
        ds >> static_cast<Editor<SegEditor,Segment>&>(segeditor);
    }
    else
        throw version_error( v, "1", r_segeditor, CODELOC );
        
    return ds;
}

/** Null constructor */
SegEditor::SegEditor() : ConcreteProperty< SegEditor,Editor<SegEditor,Segment> >()
{}

/** Construct to edit a copy of the Segment 'segment' */
SegEditor::SegEditor(const Segment &segment) 
          : ConcreteProperty< SegEditor,Editor<SegEditor,Segment> >(segment)
{}

/** Copy constructor */
SegEditor::SegEditor(const SegEditor &other) 
          : ConcreteProperty< SegEditor,Editor<SegEditor,Segment> >(other)
{}

/** Destructor */
SegEditor::~SegEditor()
{}

/** Assign this editor so that it edits a copy of the Segment 'segment' */
SegEditor& SegEditor::operator=(const Segment &segment)
{
    Editor<SegEditor,Segment>::operator=(segment);
    return *this;
}

/** Copy assignment operator */
SegEditor& SegEditor::operator=(const SegEditor &other)
{
    Editor<SegEditor,Segment>::operator=(other);
    return *this;
}

/** Return a string representation of this editor */
QString SegEditor::toString() const
{
    return QObject::tr( "Editor{ %1 }" )
                .arg( Segment::toString() );
}

/** Rename this Segment to 'newname' */
SegEditor& SegEditor::rename(const SegName &newname)
{
    if (newname == this->name())
        //nothing to do
        return *this;
        
    throw SireError::incomplete_code( CODELOC );
    
    return *this;
}

/** Move this Segment to index 'newidx' - this will move it
    to the start or end if this index is out of range */
SegStructureEditor SegEditor::reindex(SegIdx newidx) const
{
    SegStructureEditor editor(*this);
    editor.reindex(newidx);
    
    return editor;
}

/** Complete remove this Segment, and return an editor
    for the molecule that contained it */
MolStructureEditor SegEditor::remove() const
{
    SegStructureEditor editor(*this);
    
    return editor.remove();
}

/** Add an atom called 'atomname' to this Segment and return
    an editor for that atom */
AtomStructureEditor SegEditor::add(const AtomName &atomname) const
{
    SegStructureEditor editor(*this);
    
    return editor.add(atomname);
}

/** Add an atom with number 'atomnum' to this Segment and return
    an editor for that atom */
AtomStructureEditor SegEditor::add(AtomNum atomnum) const
{
    SegStructureEditor editor(*this);
    
    return editor.add(atomnum);
}

/** Completely remove all atoms that match the ID 'atomid' from 
    this Segment 
    
    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
SegStructureEditor SegEditor::remove(const AtomID &atomid) const
{
    SegStructureEditor editor(*this);
    editor.remove(atomid);
    
    return editor;
}

/** Remove the ith atom from this Segment

    \throw SireError::invalid_index
*/
SegStructureEditor SegEditor::remove(int i) const
{
    SegStructureEditor editor(*this);
    editor.remove(i);
    
    return editor;
}

/** Transfer all atoms that match the ID 'atomid' in this Segment 
    to the Segment that matches the ID 'cgid'
    
    \throw SireMol::missing_atom
    \throw SireMol::missing_segment
    \throw SireMol::duplicate_segment
    \throw SireError::invalid_index
*/
SegStructureEditor SegEditor::transfer(const AtomID &atomid, 
                                     const SegID &cgid) const
{
    SegStructureEditor editor(*this);
    editor.transfer(atomid, cgid);
    
    return editor;
}
                                     
/** Transfer the ith atom of this Segment into the Segment identified
    by the ID 'cgid'
    
    \throw SireError::invalid_index
*/
SegStructureEditor SegEditor::transfer(int i, const SegID &cgid) const
{
    SegStructureEditor editor(*this);
    editor.transfer(i, cgid);
    
    return editor;
}

/** Completely transfer all of the atoms in this Segment to 
    the Segment that matches the ID 'cgid'
    
    \throw SireMol::missing_segment
    \throw SireMol::duplicate_segment
    \throw SireError::invalid_index
*/
SegStructureEditor SegEditor::transferAll(const SegID &cgid) const
{
    SegStructureEditor editor(*this);
    editor.transferAll(cgid);
    
    return editor;
}

/** Commit the changes made by this editor and return the 
    updated Segment */
Segment SegEditor::commit() const
{
    return *this;
}

const char* SegEditor::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Segment>() );
}

////////
//////// Implementation of SegStructureEditor
////////

static const RegisterMetaType<SegStructureEditor> r_segstructeditor;

/** Serialise to a binary datastream */
QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds,
                                       const SegStructureEditor &segeditor)
{
    writeHeader(ds, r_segstructeditor, 1);
    
    ds << segeditor.uid
       << static_cast<const StructureEditor&>(segeditor);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds,
                                       SegStructureEditor &segeditor)
{
    VersionID v = readHeader(ds, r_segstructeditor);
    
    if (v == 1)
    {
        ds >> segeditor.uid
           >> static_cast<StructureEditor&>(segeditor);
    }
    else
        throw version_error( v, "1", r_segstructeditor, CODELOC );
        
    return ds;
}

/** Null constructor */
SegStructureEditor::SegStructureEditor()
                  : StructureEditor(), uid(0)
{}

/** Construct to edit a copy of the Segment 'segment' */
SegStructureEditor::SegStructureEditor(const Segment &segment)
                  : StructureEditor(segment.data())
{
    uid = this->getUID(segment.index());
}

/** Construct to edit the Segment at index 'cgidx' in the molecule 
    also being edited in 'data'
    
    \throw SireError::invalid_index
*/
SegStructureEditor::SegStructureEditor(const StructureEditor &data, SegIdx cgidx)
                  : StructureEditor(data)
{
    uid = this->getUID(cgidx);
}

/** Copy constructor */
SegStructureEditor::SegStructureEditor(const SegStructureEditor &other)
                  : StructureEditor(other), uid(other.uid)
{}

/** Destructor */
SegStructureEditor::~SegStructureEditor()
{}

/** Assign this editor so that it edits a copy of the Segment 'segment' */
SegStructureEditor& SegStructureEditor::operator=(const Segment &segment)
{
    StructureEditor::operator=(segment.data());
    uid = this->getUID(segment.index());
    
    return *this;
}

/** Copy assignment operator */
SegStructureEditor& SegStructureEditor::operator=(const SegStructureEditor &other)
{
    StructureEditor::operator=(other);
    uid = other.uid;
    
    return *this;
}

/** Return a string representation of this editor */
QString SegStructureEditor::toString() const
{
    return QObject::tr( "StructureEditor{ Segment( %1 ) }" )
                .arg( this->name() );
}

/** Return whether or not this segment is the whole molecule */
bool SegStructureEditor::selectedAll() const
{
    return StructureEditor::nSegmentsInMolecule() == 1;
}

/** Return the name of this Segment */
const SegName& SegStructureEditor::name() const
{
    return this->segName(uid);
}

/** Return the index of this Segment in the molecule */
SegIdx SegStructureEditor::index() const
{
    return this->segIdx(uid);
}

/** Return the number of atoms in this Segment (could be zero!) */
int SegStructureEditor::nAtoms() const
{
    return this->nAtomsInSegment(uid);
}

/** Return an editor for the molecule that contains this Segment */
MolStructureEditor SegStructureEditor::molecule()
{
    return MolStructureEditor(*this);
}

/** Return an editor for the ith atom of this Segment

    \throw SireError::invalid_index
*/
AtomStructureEditor SegStructureEditor::atom(int i)
{
    return AtomStructureEditor(*this, atomIdx( atomInSegment(uid,i) ));
}

/** Return an editor for the atom that matches the ID 'atomid' in
    this Segment
    
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
AtomStructureEditor SegStructureEditor::atom(const AtomID &atomid)
{
    return AtomStructureEditor(*this, atomIdx(this->index() + atomid));
}

/** Return an editor for the ith atom of this Segment

    \throw SireError::invalid_index
*/
AtomStructureEditor SegStructureEditor::select(int i)
{
    return AtomStructureEditor(*this, atomIdx( atomInSegment(uid,i) ));
}

/** Return an editor for the atom that matches the ID 'atomid' in
    this Segment
    
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
AtomStructureEditor SegStructureEditor::select(const AtomID &atomid)
{
    return AtomStructureEditor(*this, atomIdx( this->index() + atomid ));
}

/** Rename this Segment to 'newname' */
SegStructureEditor& SegStructureEditor::rename(const SegName &newname)
{
    this->renameSegment(uid, newname);
    return *this;
}

/** Move this Segment to index 'newidx' - this will move it
    to the start or end if this index is out of range */
SegStructureEditor& SegStructureEditor::reindex(SegIdx newidx)
{
    this->reindexSegment(uid, newidx);
    return *this;
}

/** Complete remove this Segment, and return an editor
    for the molecule that contained it */
MolStructureEditor SegStructureEditor::remove()
{
    this->removeSegment(uid);
    return MolStructureEditor(*this);
}

/** Add an atom called 'atomname' to this Segment and return
    an editor for that atom */
AtomStructureEditor SegStructureEditor::add(const AtomName &atomname)
{
    this->assertValidSegment(uid);

    AtomStructureEditor atom = this->addAtom();
    atom.rename(atomname);
    atom.reparent( this->index() );
    
    return atom;
}

/** Add an atom with number 'atomnum' to this Segment and return
    an editor for that atom */
AtomStructureEditor SegStructureEditor::add(AtomNum atomnum)
{
    this->assertValidSegment(uid);

    AtomStructureEditor atom = this->addAtom();
    atom.renumber(atomnum);
    atom.reparent( this->index() );
    
    return atom;
}

/** Completely remove all atoms that match the ID 'atomid' from 
    this Segment 
    
    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
SegStructureEditor& SegStructureEditor::remove(const AtomID &atomid)
{
    this->removeAtoms( this->index() + atomid );
    return *this;
}

/** Remove the ith atom from this Segment

    \throw SireError::invalid_index
*/
SegStructureEditor& SegStructureEditor::remove(int i)
{
    this->removeAtoms( atomIdx(atomInSegment(uid,i)) );
    return *this;
}

/** Transfer all atoms that match the ID 'atomid' in this Segment 
    to the Segment that matches the ID 'segid'
    
    \throw SireMol::missing_atom
    \throw SireMol::missing_segment
    \throw SireMol::duplicate_segment
    \throw SireError::invalid_index
*/
SegStructureEditor& SegStructureEditor::transfer(const AtomID &atomid, 
                                                 const SegID &segid)
{
    this->reparentAtom( this->getUID(this->index() + atomid),
                        this->segIdx(segid) );
                        
    return *this;
}

/** Transfer the ith atom of this Segment into the Segment identified
    by the ID 'segid'
    
    \throw SireError::invalid_index
*/
SegStructureEditor& SegStructureEditor::transfer(int i, const SegID &segid)
{
    this->reparentAtom( atomInSegment(uid,i), this->segIdx(segid) );
    return *this;
}

/** Completely transfer all of the atoms in this Segment to 
    the Segment that matches the ID 'segid'
    
    \throw SireMol::missing_segment
    \throw SireMol::duplicate_segment
    \throw SireError::invalid_index
*/
SegStructureEditor& SegStructureEditor::transferAll(const SegID &segid)
{
    SegIdx segidx = this->segIdx(segid);
    
    int nats = this->nAtoms();
    
    for (int i=0; i<nats; ++i)
    {
        this->reparentAtom( atomInSegment(uid,i), segidx );
    }
    
    return *this;
}

/** Commit the changes made by this editor and return the 
    updated Segment */
Segment SegStructureEditor::commit() const
{
    return Segment( this->commitChanges(), this->index() );
}

/** Allow automatic casting of this editor to a Segment */
SegStructureEditor::operator Segment() const
{
    return this->commit();
}

const char* SegStructureEditor::typeName()
{
    return QMetaType::typeName( qMetaTypeId<SegStructureEditor>() );
}

SegStructureEditor* SegStructureEditor::clone() const
{
    return new SegStructureEditor(*this);
}

Segment* Segment::clone() const
{
    return new Segment(*this);
}

