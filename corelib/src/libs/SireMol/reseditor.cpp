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

#include "reseditor.h"

#include "atomeditor.h"
#include "cgeditor.h"
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

#include "SireStream/datastream.h"

using namespace SireMol;
using namespace SireStream;

namespace SireMol
{
    //instantiate Editor<Residue> fully
    template class Editor<ResEditor,Residue>;
    //template class Editor< Selector<Residue> >;
}

//////////
////////// Implementation of ResEditor
//////////

static const RegisterMetaType<ResEditor> r_reseditor;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                       const ResEditor &reseditor)
{
    writeHeader(ds, r_reseditor, 1);

    ds << static_cast<const Editor<ResEditor,Residue>&>(reseditor);

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds,
                                       ResEditor &reseditor)
{
    VersionID v = readHeader(ds, r_reseditor);

    if (v == 1)
    {
        ds >> static_cast<Editor<ResEditor,Residue>&>(reseditor);
    }
    else
        throw version_error( v, "1", r_reseditor, CODELOC );

    return ds;
}

/** Null constructor */
ResEditor::ResEditor() : ConcreteProperty< ResEditor,Editor<ResEditor,Residue> >()
{}

/** Construct an editor that edits a copy of the residue 'residue' */
ResEditor::ResEditor(const Residue &residue)
          : ConcreteProperty< ResEditor,Editor<ResEditor,Residue> >(residue)
{}

/** Copy constructor */
ResEditor::ResEditor(const ResEditor &other)
          : ConcreteProperty< ResEditor,Editor<ResEditor,Residue> >(other)
{}

/** Destructor */
ResEditor::~ResEditor()
{}

/** Copy assignment so that this will edit a copy of 'residue */
ResEditor& ResEditor::operator=(const Residue &residue)
{
    Editor<ResEditor,Residue>::operator=(residue);
    return *this;
}

/** Copy assignment operator */
ResEditor& ResEditor::operator=(const ResEditor &other)
{
    Editor<ResEditor,Residue>::operator=(other);
    return *this;
}

/** Return a string representation of this editor */
QString ResEditor::toString() const
{
    return QObject::tr( "Editor{ %1 }" ).arg(Residue::toString());
}

/** Rename this residue to 'newname' */
ResEditor& ResEditor::rename(const ResName &newname)
{
    if (newname == this->name())
        //nothing to do
        return *this;

    d->rename( this->index(), newname );

    return *this;
}

/** Renumber this residue to 'newnum' */
ResEditor& ResEditor::renumber(ResNum newnum)
{
    if (newnum == this->number())
        //nothing to do
        return *this;

    d->renumber( this->index(), newnum );

    return *this;
}

/** Change the index of this residue to 'newidx'. If this
    is larger than the number of residues in the molecule
    then this residue is moved to the end */
ResStructureEditor ResEditor::reindex(ResIdx newidx) const
{
    ResStructureEditor editor(*this);
    editor.reindex(newidx);

    return editor;
}

/** Completely remove this residue from the molecule - this returns
    a MolStructureEditor that can be used to further edit the molecule */
MolStructureEditor ResEditor::remove() const
{
    ResStructureEditor editor(*this);
    return editor.remove();
}

/** Move this residue into the chain with ID 'chainid'

    \throw SireMol::missing_chain
    \throw SireMol::duplicate_chain
    \throw SireError::invalid_index
*/
ResStructureEditor ResEditor::reparent(const ChainID &chainid) const
{
    ResStructureEditor editor(*this);
    editor.reparent(chainid);
    return editor;
}

/** Add a new atom called 'name' to this residue - this returns
    an editor that can be used to further edit this atom */
AtomStructureEditor ResEditor::add(const AtomName &name) const
{
    ResStructureEditor editor(*this);
    return editor.add(name);
}

/** Add a new atom with the number 'number' to this residue - this
    returns an editor that can be used to further edit this atom */
AtomStructureEditor ResEditor::add(AtomNum number) const
{
    ResStructureEditor editor(*this);
    return editor.add(number);
}

/** Remove all atoms with ID 'atomid' from this residue

    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
ResStructureEditor ResEditor::remove(const AtomID &atomid) const
{
    ResStructureEditor editor(*this);
    editor.remove(atomid);
    return editor;
}

/** Remove the ith atom from this residue

    \throw SireError::invalid_index
*/
ResStructureEditor ResEditor::remove(int i) const
{
    ResStructureEditor editor(*this);
    editor.remove(i);
    return editor;
}

/** Transfer all atoms that match the ID 'atomid' into the residue that
    matches the ID 'resid'

    \throw SireMol::missing_atom
    \throw SireMol::missing_residue
    \throw SireMol::duplicate_residue
    \throw SireError::invalid_index
*/
ResStructureEditor ResEditor::transfer(const AtomID &atomid,
                                       const ResID &resid) const
{
    ResStructureEditor editor(*this);
    editor.transfer(atomid, resid);
    return editor;
}

/** Transfer the ith atom from this residue into the residue that
    matches the ID 'resid'

    \throw SireMol::missing_residue
    \throw SireMol::duplicate_residue
    \throw SireError::invalid_index
*/
ResStructureEditor ResEditor::transfer(int i, const ResID &resid) const
{
    ResStructureEditor editor(*this);
    editor.transfer(i, resid);
    return editor;
}

/** Transfer all atoms from this residue into the residue with ID 'resid'

    \throw SireMol::missing_residue
    \throw SireMol::duplicate_residue
    \throw SireError::invalid_index
*/
ResStructureEditor ResEditor::transferAll(const ResID &resid) const
{
    ResStructureEditor editor(*this);
    editor.transferAll(resid);
    return editor;
}

/** Commit the changes made by this editor and return the updated Residue */
Residue ResEditor::commit() const
{
    return *this;
}

const char* ResEditor::typeName()
{
    return QMetaType::typeName( qMetaTypeId<ResEditor>() );
}

//////////
////////// Implementation of ResStructureEditor
//////////

static const RegisterMetaType<ResStructureEditor> r_resstructeditor;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                       const ResStructureEditor &reseditor)
{
    writeHeader(ds, r_resstructeditor, 1);

    ds << reseditor.uid
       << static_cast<const StructureEditor&>(reseditor);

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds,
                                       ResStructureEditor &reseditor)
{
    VersionID v = readHeader(ds, r_resstructeditor);

    if (v == 1)
    {
        ds >> reseditor.uid
           >> static_cast<StructureEditor&>(reseditor);
    }
    else
        throw version_error( v, "1", r_resstructeditor, CODELOC );

    return ds;
}

/** Null constructor */
ResStructureEditor::ResStructureEditor()
                   : StructureEditor(), uid(0)
{
    this->operator=( StructureEditor::addResidue() );
}

/** Construct an editor to edit the structure of a copy of the
    residue 'residue' */
ResStructureEditor::ResStructureEditor(const Residue &residue)
                   : StructureEditor(residue.data())
{
    uid = this->getUID(residue.index());
}

/** Construct an editor to edit the residue at index 'residx' in the
    editor whose data is in 'data'

    \throw SireError::invalid_index
*/
ResStructureEditor::ResStructureEditor(const StructureEditor &data, ResIdx residx)
                   : StructureEditor(data)
{
    uid = this->getUID(residx);
}

/** Copy constructor */
ResStructureEditor::ResStructureEditor(const ResStructureEditor &other)
                   : StructureEditor(other), uid(other.uid)
{}

/** Destructor */
ResStructureEditor::~ResStructureEditor()
{}

/** Assign so that this edits a copy of 'residue' */
ResStructureEditor& ResStructureEditor::operator=(const Residue &residue)
{
    StructureEditor::operator=(residue.data());
    uid = this->getUID(residue.index());

    return *this;
}

/** Copy assignment operator */
ResStructureEditor& ResStructureEditor::operator=(const ResStructureEditor &other)
{
    StructureEditor::operator=(other);
    uid = other.uid;

    return *this;
}

/** Return a string representation of this editor */
QString ResStructureEditor::toString() const
{
    return QObject::tr( "StructureEditor{ Residue( %1 : %2 ) }" )
                .arg( this->name() )
                .arg( this->number() );
}

/** Is this editor editing the entire molecule? */
bool ResStructureEditor::selectedAll() const
{
    return StructureEditor::nResiduesInMolecule() == 1;
}

/** Return the name of this residue */
const ResName& ResStructureEditor::name() const
{
    return this->resName(uid);
}

/** Return the number of this residue */
ResNum ResStructureEditor::number() const
{
    return this->resNum(uid);
}

/** Return the index of this residue in the molecule */
ResIdx ResStructureEditor::index() const
{
    return this->resIdx(uid);
}

/** Return the number of atoms in this residue - this may be zero! */
int ResStructureEditor::nAtoms() const
{
    return this->nAtomsInResidue(uid);
}

/** Return an editor for the chain that contains this residue */
ChainStructureEditor ResStructureEditor::chain()
{
    return ChainStructureEditor(*this, chainIdx(chainParentOfResidue(uid)));
}

/** Return an editor for the molecule that contains this residue */
MolStructureEditor ResStructureEditor::molecule()
{
    return MolStructureEditor(*this);
}

/** Return an editor for the ith atom in this residue

    \throw SireError::invalid_index
*/
AtomStructureEditor ResStructureEditor::atom(int i)
{
    return AtomStructureEditor(*this, atomIdx(atomInResidue(uid,i)));
}

/** Return an editor for the atom with ID == 'atomid' in
    this residue

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
AtomStructureEditor ResStructureEditor::atom(const AtomID &atomid)
{
    return AtomStructureEditor(*this, atomIdx(atomid));
}

/** Return an editor for the ith atom in this residue

    \throw SireError::invalid_index
*/
AtomStructureEditor ResStructureEditor::select(int i)
{
    return AtomStructureEditor(*this, atomIdx(atomInResidue(uid,i)));
}

/** Return an editor for the atom with ID == 'atomid' in
    this residue

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
AtomStructureEditor ResStructureEditor::select(const AtomID &atomid)
{
    return AtomStructureEditor(*this, atomIdx(atomid));
}

/** Rename this residue to 'newname' */
ResStructureEditor& ResStructureEditor::rename(const ResName &newname)
{
    this->renameResidue(uid, newname);
    return *this;
}

/** Renumber this residue to 'newnum' */
ResStructureEditor& ResStructureEditor::renumber(ResNum newnum)
{
    this->renumberResidue(uid, newnum);
    return *this;
}

/** Change the index of this residue to 'newidx'. If this
    is larger than the number of residues in the molecule
    then this residue is moved to the end */
ResStructureEditor& ResStructureEditor::reindex(ResIdx newidx)
{
    this->reindexResidue(uid, newidx);
    return *this;
}

/** Completely remove this residue from the molecule - this returns
    a MolStructureEditor that can be used to further edit the molecule */
MolStructureEditor ResStructureEditor::remove()
{
    this->removeResidue(uid);
    return MolStructureEditor(*this);
}

/** Move this residue into the chain with ID 'chainid'

    \throw SireMol::missing_chain
    \throw SireMol::duplicate_chain
    \throw SireError::invalid_index
*/
ResStructureEditor& ResStructureEditor::reparent(const ChainID &chainid)
{
    this->reparentResidue(uid, this->chainIdx(chainid));
    return *this;
}

/** Add a new atom called 'name' to this residue - this returns
    an editor that can be used to further edit this atom */
AtomStructureEditor ResStructureEditor::add(const AtomName &atomname)
{
    this->assertValidResidue(uid);

    AtomStructureEditor atom = this->addAtom();
    atom = atom.rename(atomname);
    atom = atom.reparent(this->index());

    return atom;
}

/** Add a new atom with the number 'number' to this residue - this
    returns an editor that can be used to further edit this atom */
AtomStructureEditor ResStructureEditor::add(AtomNum atomnum)
{
    this->assertValidResidue(uid);

    AtomStructureEditor atom = this->addAtom();
    atom = atom.renumber(atomnum);
    atom = atom.reparent(this->index());

    return atom;
}

/** Remove all atoms with ID 'atomid' from this residue

    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
ResStructureEditor& ResStructureEditor::remove(const AtomID &atomid)
{
    this->removeAtoms( this->index() + atomid );
    return *this;
}

/** Remove the ith atom from this residue

    \throw SireError::invalid_index
*/
ResStructureEditor& ResStructureEditor::remove(int i)
{
    this->removeAtoms( atomIdx(atomInResidue(uid, i)) );
    return *this;
}

/** Transfer all atoms that match the ID 'atomid' into the residue that
    matches the ID 'resid'

    \throw SireMol::missing_atom
    \throw SireMol::missing_residue
    \throw SireMol::duplicate_residue
    \throw SireError::invalid_index
*/
ResStructureEditor& ResStructureEditor::transfer(const AtomID &atomid,
                                                 const ResID &resid)
{
    this->reparentAtom( this->getUID(this->index() + atomid),
                        this->resIdx(resid) );

    return *this;
}

/** Transfer the ith atom from this residue into the residue that
    matches the ID 'resid'

    \throw SireMol::missing_residue
    \throw SireMol::duplicate_residue
    \throw SireError::invalid_index
*/
ResStructureEditor& ResStructureEditor::transfer(int i, const ResID &resid)
{
    this->reparentAtom( atomInResidue(uid, i), this->resIdx(resid) );
    return *this;
}

/** Transfer all atoms from this residue into the residue with ID 'resid'

    \throw SireMol::missing_residue
    \throw SireMol::duplicate_residue
    \throw SireError::invalid_index
*/
ResStructureEditor& ResStructureEditor::transferAll(const ResID &resid)
{
    ResIdx residx = this->resIdx(resid);

    int nats = this->nAtoms();

    for (int i=0; i<nats; ++i)
    {
        this->reparentAtom( atomInResidue(uid,i), residx );
    }

    return *this;
}

/** Commit the changes made by this editor and return the updated residue */
Residue ResStructureEditor::commit() const
{
    return Residue( this->commitChanges(), this->index() );
}

/** Allow automatic casting to a Residue() */
ResStructureEditor::operator Residue() const
{
    return this->commit();
}

const char* ResStructureEditor::typeName()
{
    return QMetaType::typeName( qMetaTypeId<ResStructureEditor>() );
}

ResStructureEditor* ResStructureEditor::clone() const
{
    return new ResStructureEditor(*this);
}
