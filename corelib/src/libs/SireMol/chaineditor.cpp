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
#include "chaineditor.h"
#include "reseditor.h"
#include "chaineditor.h"
#include "segeditor.h"
#include "moleditor.h"

#include "atom.h"
#include "chain.h"
#include "residue.h"
#include "chain.h"
#include "segment.h"
#include "molecule.h"
#include "mover.hpp"
#include "selector.hpp"

#include "chainresid.h"
#include "groupatomids.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"

using namespace SireMol;
using namespace SireStream;

namespace SireMol
{
    // fully instantiate Editor<Chain>
    template class Editor<ChainEditor,Chain>;
    //template class Selector< Editor<Chain> >;
}

////////
//////// Implementation of ChainEditor
////////

static const RegisterMetaType<ChainEditor> r_chaineditor;

/** Serialise to a binary datastream */
QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds,
                                       const ChainEditor &chaineditor)
{
    writeHeader(ds, r_chaineditor, 1);
    
    ds << static_cast<const Editor<ChainEditor,Chain>&>(chaineditor);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds,
                                       ChainEditor &chaineditor)
{
    VersionID v = readHeader(ds, r_chaineditor);
    
    if (v == 1)
    {
        ds >> static_cast<Editor<ChainEditor,Chain>&>(chaineditor);
    }
    else
        throw version_error( v, "1", r_chaineditor, CODELOC );
        
    return ds;
}

/** Null constructor */
ChainEditor::ChainEditor() : ConcreteProperty< ChainEditor,Editor<ChainEditor,Chain> >()
{}

/** Construct to edit a copy of the Chain 'chain' */
ChainEditor::ChainEditor(const Chain &chain) 
            : ConcreteProperty< ChainEditor,Editor<ChainEditor,Chain> >(chain)
{}

/** Copy constructor */
ChainEditor::ChainEditor(const ChainEditor &other) 
            : ConcreteProperty< ChainEditor,Editor<ChainEditor,Chain> >(other)
{}

/** Destructor */
ChainEditor::~ChainEditor()
{}

/** Assign this editor so that it edits a copy of the Chain 'chain' */
ChainEditor& ChainEditor::operator=(const Chain &chain)
{
    Editor<ChainEditor,Chain>::operator=(chain);
    return *this;
}

/** Copy assignment operator */
ChainEditor& ChainEditor::operator=(const ChainEditor &other)
{
    Editor<ChainEditor,Chain>::operator=(other);
    return *this;
}

/** Return a string representation of this chain */
QString ChainEditor::toString() const
{
    return QObject::tr( "Editor{ %1 }" ).arg( Chain::toString() );
}

/** Rename this Chain to 'newname' */
ChainEditor& ChainEditor::rename(const ChainName &newname)
{
    if (newname == this->name())
        //nothing to do
        return *this;
        
    d->rename( this->index(), newname );
    
    return *this;
}

/** Move this Chain to index 'newidx' - this will move it
    to the start or end if this index is out of range */
ChainStructureEditor ChainEditor::reindex(ChainIdx newidx) const
{
    ChainStructureEditor editor(*this);
    editor.reindex(newidx);
    
    return editor;
}

/** Complete remove this Chain, and return an editor
    for the molecule that contained it */
MolStructureEditor ChainEditor::remove() const
{
    ChainStructureEditor editor(*this);
    
    return editor.remove();
}

/** Add a residue called 'resname' to this Chain and return
    an editor for that residue */
ResStructureEditor ChainEditor::add(const ResName &resname) const
{
    ChainStructureEditor editor(*this);
    
    return editor.add(resname);
}

/** Add a residue with number 'resnum' to this Chain and return
    an editor for that residue */
ResStructureEditor ChainEditor::add(ResNum resnum) const
{
    ChainStructureEditor editor(*this);
    
    return editor.add(resnum);
}

/** Completely remove all residues that match the ID 'resid' from 
    this Chain 
    
    \throw SireMol::missing_residue
    \throw SireError::invalid_index
*/
ChainStructureEditor ChainEditor::remove(const ResID &resid) const
{
    ChainStructureEditor editor(*this);
    editor.remove(resid);
    
    return editor;
}

/** Remove the ith residue from this Chain

    \throw SireError::invalid_index
*/
ChainStructureEditor ChainEditor::remove(int i) const
{
    ChainStructureEditor editor(*this);
    editor.remove(i);
    
    return editor;
}

/** Remove the atom that matches the ID 'atomid' from this chain

    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
ChainStructureEditor ChainEditor::remove(const AtomID &atomid) const
{
    ChainStructureEditor editor(*this);
    editor.remove(atomid);
    
    return editor;
}

/** Transfer all residues that match the ID 'resid' in this Chain 
    to the Chain that matches the ID 'cgid'
    
    \throw SireMol::missing_residue
    \throw SireMol::missing_chain
    \throw SireMol::duplicate_chain
    \throw SireError::invalid_index
*/
ChainStructureEditor ChainEditor::transfer(const ResID &resid, 
                                           const ChainID &cgid) const
{
    ChainStructureEditor editor(*this);
    editor.transfer(resid, cgid);
    
    return editor;
}
                                     
/** Transfer the ith residue of this Chain into the Chain identified
    by the ID 'cgid'
    
    \throw SireError::invalid_index
*/
ChainStructureEditor ChainEditor::transfer(int i, const ChainID &cgid) const
{
    ChainStructureEditor editor(*this);
    editor.transfer(i, cgid);
    
    return editor;
}

/** Completely transfer all of the residues in this Chain to 
    the Chain that matches the ID 'cgid'
    
    \throw SireMol::missing_chain
    \throw SireMol::duplicate_chain
    \throw SireError::invalid_index
*/
ChainStructureEditor ChainEditor::transferAll(const ChainID &cgid) const
{
    ChainStructureEditor editor(*this);
    editor.transferAll(cgid);
    
    return editor;
}

/** Commit the changes made by this editor and return the 
    updated Chain */
Chain ChainEditor::commit() const
{
    return *this;
}

const char* ChainEditor::typeName()
{
    return QMetaType::typeName( qMetaTypeId<ChainEditor>() );
}

////////
//////// Implementation of ChainStructureEditor
////////

static const RegisterMetaType<ChainStructureEditor> r_cgstructeditor;

/** Serialise to a binary datastream */
QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds,
                                       const ChainStructureEditor &chaineditor)
{
    writeHeader(ds, r_cgstructeditor, 1);
    
    ds << chaineditor.uid
       << static_cast<const StructureEditor&>(chaineditor);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds,
                                       ChainStructureEditor &chaineditor)
{
    VersionID v = readHeader(ds, r_cgstructeditor);
    
    if (v == 1)
    {
        ds >> chaineditor.uid
           >> static_cast<StructureEditor&>(chaineditor);
    }
    else
        throw version_error( v, "1", r_cgstructeditor, CODELOC );
        
    return ds;
}

/** Null constructor */
ChainStructureEditor::ChainStructureEditor()
                  : StructureEditor(), uid(0)
{}

/** Construct to edit a copy of the Chain 'chain' */
ChainStructureEditor::ChainStructureEditor(const Chain &chain)
                  : StructureEditor(chain.data())
{
    uid = this->getUID(chain.index());
}

/** Construct to edit the Chain at index 'cgidx' in the molecule 
    also being edited in 'data'
    
    \throw SireError::invalid_index
*/
ChainStructureEditor::ChainStructureEditor(const StructureEditor &data, ChainIdx cgidx)
                  : StructureEditor(data)
{
    uid = this->getUID(cgidx);
}

/** Copy constructor */
ChainStructureEditor::ChainStructureEditor(const ChainStructureEditor &other)
                  : StructureEditor(other), uid(other.uid)
{}

/** Destructor */
ChainStructureEditor::~ChainStructureEditor()
{}

/** Assign this editor so that it edits a copy of the Chain 'chain' */
ChainStructureEditor& ChainStructureEditor::operator=(const Chain &chain)
{
    StructureEditor::operator=(chain.data());
    uid = this->getUID(chain.index());
    
    return *this;
}

/** Copy assignment operator */
ChainStructureEditor& ChainStructureEditor::operator=(const ChainStructureEditor &other)
{
    StructureEditor::operator=(other);
    uid = other.uid;
    
    return *this;
}

/** Return a string representation of this chain */
QString ChainStructureEditor::toString() const
{
    return QObject::tr( "StructureEditor{ Chain( %1 ) }" )
                .arg( this->name() );
}

/** Return whether or not this chain is the whole molecule */
bool ChainStructureEditor::selectedAll() const
{
    return StructureEditor::nChainsInMolecule() == 1;
}

/** Return the name of this Chain */
const ChainName& ChainStructureEditor::name() const
{
    return this->chainName(uid);
}

/** Return the index of this Chain in the molecule */
ChainIdx ChainStructureEditor::index() const
{
    return this->chainIdx(uid);
}

/** Return the number of atoms in this Chain (could be zero!) */
int ChainStructureEditor::nAtoms() const
{
    return this->nAtomsInChain(uid);
}

/** Return the number of residues in this Chain (could be zero!) */
int ChainStructureEditor::nResidues() const
{
    return this->nResiduesInChain(uid);
}

/** Return an editor for the molecule that contains this Chain */
MolStructureEditor ChainStructureEditor::molecule()
{
    return MolStructureEditor(*this);
}

/** Return an editor for the ith residue of this Chain

    \throw SireError::invalid_index
*/
ResStructureEditor ChainStructureEditor::residue(int i)
{
    return ResStructureEditor(*this, resIdx( residueInChain(uid,i) ));
}

/** Return an editor for the residue that matches the ID 'resid' in
    this chain
    
    \throw SireMol::missing_residue
    \throw SireMol::duplicate_residue
    \throw SireError::invalid_index
*/
ResStructureEditor ChainStructureEditor::residue(const ResID &resid)
{
    return ResStructureEditor(*this, resIdx(this->index() + resid));
}

/** Return an editor for the atom that matches the ID 'atomid' in
    this Chain
    
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
AtomStructureEditor ChainStructureEditor::atom(const AtomID &atomid)
{
    return AtomStructureEditor(*this, atomIdx(this->index() + atomid));
}

/** Return an editor for the ith residue of this Chain

    \throw SireError::invalid_index
*/
ResStructureEditor ChainStructureEditor::select(int i)
{
    return ResStructureEditor(*this, resIdx( residueInChain(uid,i) ));
}

/** Return an editor for the atom that matches the ID 'atomid' in
    this Chain
    
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
AtomStructureEditor ChainStructureEditor::select(const AtomID &atomid)
{
    return AtomStructureEditor(*this, atomIdx( this->index() + atomid ));
}

/** Return an editor for the residue that matches the ID 'resid' in
    this chain
    
    \throw SireMol::missing_residue
    \throw SireMol::duplicate_residue
    \throw SireError::invalid_index
*/
ResStructureEditor ChainStructureEditor::select(const ResID &resid)
{
    return ResStructureEditor(*this, resIdx( this->index() + resid ));
}

/** Rename this Chain to 'newname' */
ChainStructureEditor& ChainStructureEditor::rename(const ChainName &newname)
{
    this->renameChain(uid, newname);
    return *this;
}

/** Move this Chain to index 'newidx' - this will move it
    to the start or end if this index is out of range */
ChainStructureEditor& ChainStructureEditor::reindex(ChainIdx newidx)
{
    this->reindexChain(uid, newidx);
    return *this;
}

/** Complete remove this Chain, and return an editor
    for the molecule that contained it */
MolStructureEditor ChainStructureEditor::remove()
{
    this->removeChain(uid);
    return MolStructureEditor(*this);
}

/** Add a residue called 'resname' to this Chain and return
    an editor for that residue */
ResStructureEditor ChainStructureEditor::add(const ResName &resname)
{
    this->assertValidChain(uid);

    ResStructureEditor residue = this->addResidue();
    residue.rename(resname);
    residue.reparent( this->index() );
    
    return residue;
}

/** Add a residue with number 'resnum' to this Chain and return
    an editor for that residue */
ResStructureEditor ChainStructureEditor::add(ResNum resnum)
{
    this->assertValidChain(uid);

    ResStructureEditor residue = this->addResidue();
    residue.renumber(resnum);
    residue.reparent( this->index() );
    
    return residue;
}

/** Completely remove all residues that match the ID 'resid' from 
    this Chain 
    
    \throw SireMol::missing_residue
    \throw SireError::invalid_index
*/
ChainStructureEditor& ChainStructureEditor::remove(const ResID &resid)
{
    this->removeResidues( this->index() + resid );
    return *this;
}

/** Remove the ith residue from this Chain

    \throw SireError::invalid_index
*/
ChainStructureEditor& ChainStructureEditor::remove(int i)
{
    this->removeResidues( resIdx(residueInChain(uid,i)) );
    return *this;
}

/** Completely remove all atoms that match the ID 'atomid' from
    all of the residues that are part of this chain
    
    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
ChainStructureEditor& ChainStructureEditor::remove(const AtomID &atomid)
{
    this->removeAtoms( this->index() + atomid );
    return *this;
}

/** Transfer all residues that match the ID 'resid' in this Chain 
    to the Chain that matches the ID 'chainid'
    
    \throw SireMol::missing_residue
    \throw SireMol::missing_chain
    \throw SireMol::duplicate_chain
    \throw SireError::invalid_index
*/
ChainStructureEditor& ChainStructureEditor::transfer(const ResID &resid, 
                                                     const ChainID &chainid)
{
    this->reparentResidue( this->getUID(this->index() + resid),
                           this->chainIdx(chainid) );
                        
    return *this;
}

/** Transfer the ith residue of this Chain into the Chain identified
    by the ID 'chainid'
    
    \throw SireError::invalid_index
*/
ChainStructureEditor& ChainStructureEditor::transfer(int i, 
                                                     const ChainID &chainid)
{
    this->reparentResidue( residueInChain(uid,i), this->chainIdx(chainid) );
    return *this;
}

/** Completely transfer all of the residues in this Chain to 
    the Chain that matches the ID 'cgid'
    
    \throw SireMol::missing_chain
    \throw SireMol::duplicate_chain
    \throw SireError::invalid_index
*/
ChainStructureEditor& ChainStructureEditor::transferAll(const ChainID &chainid)
{
    ChainIdx chainidx = this->chainIdx(chainid);
    
    int nres = this->nResidues();
    
    for (int i=0; i<nres; ++i)
    {
        this->reparentResidue( residueInChain(uid,i), chainidx );
    }
    
    return *this;
}

/** Commit the changes made by this editor and return the 
    updated Chain */
Chain ChainStructureEditor::commit() const
{
    return Chain( this->commitChanges(), this->index() );
}

/** Allow automatic casting of this editor to a Chain */
ChainStructureEditor::operator Chain() const
{
    return this->commit();
}

const char* ChainStructureEditor::typeName()
{
    return QMetaType::typeName( qMetaTypeId<ChainStructureEditor>() );
}

ChainStructureEditor* ChainStructureEditor::clone() const
{
    return new ChainStructureEditor(*this);
}
