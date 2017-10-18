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

#include "atomselection.h"
#include "moleculeview.h"
#include "moleculedata.h"

#include "moleculeinfodata.h"

#include "SireError/errors.h"
#include "SireMol/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>

using namespace SireMol;
using namespace SireBase;
using namespace SireStream;

RegisterMetaType<AtomSelection> r_selection;

/** Serialise to a binary datastream */
QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds,
                                       const AtomSelection &selection)
{
    writeHeader(ds, r_selection, 2);
    
    SharedDataStream sds(ds);

    sds << selection.d << selection.selected_atoms 
        << selection.nselected;
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds,
                                       AtomSelection &selection)
{
    VersionID v = readHeader(ds, r_selection);
    
    if (v == 2)
    {
        SharedDataStream sds(ds);
        
        sds >> selection.d >> selection.selected_atoms
            >> selection.nselected;
    }
    else if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> selection.d >> selection.selected_atoms
            >> selection.nselected;
        
        if (not (selection.selectedAll() or selection.selectedNone()))
        {
            //we need to invert the CutGroup selections.
            QHash< CGIdx, QSet<Index> > s;
            
            for (CGIdx i(0); i<selection.d.read().nCutGroups(); ++i)
            {
                if (not selection.selected_atoms.contains(i))
                {
                    //all atoms are selected in this CutGroup
                    //this is represneted now as an empty set
                    s.insert(i, QSet<Index>());
                }
                else
                {
                    const auto cgatoms = selection.selected_atoms.value(i);
                    
                    if (not cgatoms.isEmpty())
                    {
                        //some of the atoms are selected
                        s.insert(i, cgatoms);
                    }
                }
            }
            
            selection.selected_atoms = s;
        }
    }
    else
        throw version_error(v, "1,2", r_selection, CODELOC);
        
    return ds;
}

/** Null constructor */
AtomSelection::AtomSelection() 
              : ConcreteProperty<AtomSelection,MoleculeProperty>(),
                d(MoleculeInfoData::null()),
                nselected(0)
{}

/** Return the info object for the molecule whose atoms
    are being selected */
const MoleculeInfoData& AtomSelection::info() const
{
    return *d;
}

/** Construct a selection of all of the atoms in the view 
    'molecule' */
AtomSelection::AtomSelection(const MoleculeView &molecule)
              : ConcreteProperty<AtomSelection,MoleculeProperty>()
{
    this->operator=(molecule.selection());
}

/** Construct a selection of all of the atoms in the 
    molecule whose data is in 'moldata' */
AtomSelection::AtomSelection(const MoleculeData &moldata)
              : ConcreteProperty<AtomSelection,MoleculeProperty>(),
                d(moldata.info()), nselected(moldata.info().nAtoms())
{}              

/** Construct a selection of all of the atoms in the
    molecule that is described by the info object 'molinfo' */
AtomSelection::AtomSelection(const MoleculeInfoData &molinfo)
              : ConcreteProperty<AtomSelection,MoleculeProperty>(),
                d(molinfo), nselected(molinfo.nAtoms())
{}

/** Copy constructor */
AtomSelection::AtomSelection(const AtomSelection &other)
              : ConcreteProperty<AtomSelection,MoleculeProperty>(other),
                selected_atoms(other.selected_atoms),
                d(other.d), nselected(other.nselected)
{}

/** Destructor */
AtomSelection::~AtomSelection()
{}

/** Copy assignment operator */
AtomSelection& AtomSelection::operator=(const AtomSelection &other)
{
    MoleculeProperty::operator=(other);

    selected_atoms = other.selected_atoms;
    d = other.d;
    nselected = other.nselected;
    
    return *this;
}

/** Comparison operator */
bool AtomSelection::operator==(const AtomSelection &other) const
{
    return this == &other or
           ( nselected == other.nselected and 
             (d == other.d or *d == *(other.d)) and
             selected_atoms == other.selected_atoms );
}

/** Comparison operator */
bool AtomSelection::operator!=(const AtomSelection &other) const
{
    return this != &other and
           ( nselected != other.nselected or
             (d != other.d and *d != *(other.d)) or
             selected_atoms != other.selected_atoms );
}

/** Return wheter no atoms are selected */
bool AtomSelection::isEmpty() const
{
    return nselected == 0;
}

/** Return whether or not this is a null selection */
bool AtomSelection::isNull() const
{
    return nselected == 0 and info() == MoleculeInfoData::null();
}

/** Return the number of selected atoms */
int AtomSelection::nSelected() const
{
    return nselected;
}

bool AtomSelection::_pvt_selected(const CGAtomIdx &cgatomidx) const
{
    auto it = selected_atoms.constFind( cgatomidx.cutGroup() );
    
    if (it != selected_atoms.constEnd())
    {
        return (it->isEmpty() or it->contains(cgatomidx.atom()));
    }
    else
        return this->selectedAll();
}

bool AtomSelection::_pvt_selected(AtomIdx atomidx) const
{
    return this->_pvt_selected( d->cgAtomIdx(atomidx) );
}

/** Return whether or not the atom at index 'cgatomidx' has
    been selected
    
    \throw SireError::invalid_index
*/
bool AtomSelection::selected(const CGAtomIdx &cgatomidx) const
{
    CGAtomIdx sane_idx( CGIdx(cgatomidx.cutGroup().map(info().nCutGroups())),
                        Index(cgatomidx.atom()
                                        .map(info().nAtoms(cgatomidx.cutGroup()))) );
                                              
    return this->_pvt_selected(sane_idx);
}

/** Return whether the atom at index 'atomidx' has been selected 

    \throw SireError::invalid_index
*/
bool AtomSelection::selected(AtomIdx atomidx) const
{
    return this->_pvt_selected( info().cgAtomIdx(atomidx) );
}

/** Return whether any of the atom(s) identified by the ID 'atomid'
    have been selected
    
    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
bool AtomSelection::selected(const AtomID &atomid) const
{
    for (auto atomidx : info().cgAtomIdxs(atomid))
    {
        if (this->_pvt_selected(atomidx))
            return true;
    }
    
    return false;
}

/** Return whether or not any atom in the CutGroup
    at index 'cgidx' has been selected
    
    \throw SireError::invalid_index
*/
bool AtomSelection::selected(CGIdx cgidx) const
{
    cgidx = CGIdx( cgidx.map(info().nCutGroups()) );
    return this->selectedAll() or selected_atoms.contains(cgidx);
}

/** Return whether or not any atoms in the residue
    at index 'residx' has been selected
    
    \throw SireError::invalid_index
*/
bool AtomSelection::selected(ResIdx residx) const
{
    residx = ResIdx( residx.map(info().nResidues()) );
    
    if (selected_atoms.isEmpty())
        return nselected > 0;
    else
    {
        for ( const auto atomidx : info().cgAtomIdxs(residx) )
        {
            if (this->_pvt_selected(atomidx))
                return true;
        }
        
        return false;
    }
}

/** Return whether or not any atoms in the chain
    at index 'chainidx' have been selected
    
    \throw SireError::invalid_index
*/
bool AtomSelection::selected(ChainIdx chainidx) const
{
    chainidx = ChainIdx( chainidx.map(info().nChains()) );
    
    if (selected_atoms.isEmpty())
        return nselected > 0;
    else
    {
        foreach( ResIdx residx, info().getResiduesIn(chainidx) )
        {
            foreach( const CGAtomIdx &atomidx, info().cgAtomIdxs(residx) )
            {
                if (this->_pvt_selected(atomidx))
                    return true;
            }
        }
        
        return false;
    }
}

/** Return whether or not any atoms in the segment at 
    index 'segidx' have been selected
    
    \throw SireError::invalid_index
*/
bool AtomSelection::selected(SegIdx segidx) const
{
    segidx = SegIdx( segidx.map(info().nSegments()) );
    
    if (selected_atoms.isEmpty())
        return nselected > 0;
    else
    {
        foreach( const CGAtomIdx &atomidx, info().cgAtomIdxs(segidx) )
        {
            if (this->_pvt_selected(atomidx))
                return true;
        }
        
        return false;
    }
}

/** Return whether any atoms in the CutGroup(s) identified
    by 'cgid' have been selected 
    
    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
*/
bool AtomSelection::selected(const CGID &cgid) const
{
    foreach (CGIdx cgidx, cgid.map(info()))
    {
        if (this->selected(cgidx))
            return true;
    }
    
    return false;
}

/** Return whether any atoms in the residue(s) identified
    by 'resid' have been selected
    
    \throw SireMol::missing_residue
    \throw SireError::invalid_index
*/
bool AtomSelection::selected(const ResID &resid) const
{
    foreach (ResIdx residx, resid.map(info()))
    {
        if (this->selected(residx))
            return true;
    }
    
    return false;
}

/** Return whether any atoms in the chain(s) identified
    by 'chainid' have been selected 
    
    \throw SireMol::missing_chain
    \throw SireError::invalid_index
*/
bool AtomSelection::selected(const ChainID &chainid) const
{
    foreach (ChainIdx chainidx, chainid.map(info()))
    {
        if (this->selected(chainidx))
            return true;
    }
    
    return false;
}

/** Return whether any atoms in the segment(s) identified
    by 'segid' have been selected 
    
    \throw SireMol::missing_segment
    \throw SireError::invalid_index
*/
bool AtomSelection::selected(const SegID &segid) const
{
    foreach (SegIdx segidx, segid.map(info()))
    {
        if (this->selected(segidx))
            return true;
    }
    
    return false;
}

/** Return whether or not any of the atoms selected in 'selection'
    are also selected in this set 
    
    \throw SireError::incompatible_error
*/
bool AtomSelection::selected(const AtomSelection &selection) const
{
    this->info().assertEqualTo( selection.info() );
    
    //get rid of the easy cases...
    if (this->isEmpty() or selection.isEmpty())
        return false;
    else if (this->selectedAll() or selection.selectedAll())
        return true;
            
    //now do the hard stuff!
    if (selection.selectedAllCutGroups())
    {
        for (CGIdx i(0); i < this->info().nCutGroups(); ++i)
        {
            if (this->selectedAll(i))
                return true;
            else if (selection.selectedAll(i) and not this->selectedNone(i))
                return true;
            else if (not this->selectedNone(i))
            {
                const QSet<Index> &atoms = *(selected_atoms.find(i));
                
                foreach (Index idx, selection.selectedAtoms(i))
                {
                    if (atoms.contains(idx))
                        return true;
                }
            }
        }
    }
    else
    {
        QList<CGIdx> cgidxs = selection.selectedCutGroups();
        
        foreach (CGIdx i, cgidxs)
        {
            if (this->selectedAll(i))
                return true;
            else if (selection.selectedAll(i) and not this->selectedNone(i))
                return true;
            else if (not this->selectedNone(i))
            {
                const QSet<Index> &atoms = *(selected_atoms.find(i));
                
                foreach (Index idx, selection.selectedAtoms(i))
                {
                    if (atoms.contains(idx))
                        return true;
                }
            }
        }
    }
    
    return false;
}

/** Return whether all atoms have been selected */
bool AtomSelection::selectedAllAtoms() const
{
    return this->selectedAll();
}

/** Return whether all CutGroups contain at least
    one selected atom */
bool AtomSelection::selectedAllCutGroups() const
{
    if (selected_atoms.isEmpty())
        return nselected > 0;
    else
        return selected_atoms.count() == info().nCutGroups();
}

/** Return whether all residues contain at least
    one selected atom */
bool AtomSelection::selectedAllResidues() const
{
    if (selected_atoms.isEmpty())
        return nselected > 0;
    else
    {
        for (ResIdx i(0); i<info().nResidues(); ++i)
        {
            if (not this->selected(i))
                return false;
        }
        
        return true;
    }
}

/** Return whether all chains contain at least
    one selected atom */
bool AtomSelection::selectedAllChains() const
{
    if (selected_atoms.isEmpty())
        return nselected > 0;
    else
    {
        for (ChainIdx i(0); i<info().nChains(); ++i)
        {
            if (not this->selected(i))
                return false;
        }
        
        return true;
    }
}

/** Return whether all segments contain at least
    one selected atom */
bool AtomSelection::selectedAllSegments() const
{
    if (selected_atoms.isEmpty())
        return nselected > 0;
    else
    {
        for (SegIdx i(0); i<info().nSegments(); ++i)
        {
            if (not this->selected(i))
                return false;
        }
        
        return true;
    }
}

/** Return whether or not no atoms have been selected */
bool AtomSelection::selectedNone() const
{
    return this->isEmpty();
}

/** Return whether the atom at index 'atomidx' has not been selected

    \throw SireError::invalid_index
*/
bool AtomSelection::selectedNone(AtomIdx atomidx) const
{
    return not this->selected(atomidx);
}

/** Return whether none of the atoms in the CutGroup at
    index 'cgidx' have been selected
    
    \throw SireError::invalid_index
*/
bool AtomSelection::selectedNone(CGIdx cgidx) const
{
    return not this->selected(cgidx);
}

/** Return whether none of the atoms in the residue at
    index 'residx' have been selected
    
    \throw SireError::invalid_index
*/
bool AtomSelection::selectedNone(ResIdx residx) const
{
    return not this->selected(residx);
}

/** Return whether none of the atoms in the chain at
    index 'chainidx' have been selected
    
    \throw SireError::invalid_index
*/
bool AtomSelection::selectedNone(ChainIdx chainidx) const
{
    return not this->selected(chainidx);
}

/** Return whether none of the atoms in the segment at
    index 'segidx' have been selected
    
    \throw SireError::invalid_index
*/
bool AtomSelection::selectedNone(SegIdx segidx) const
{
    return not this->selected(segidx);
}

/** Return whether none of the atoms identified by 'atomid'
    have been selected
    
    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
bool AtomSelection::selectedNone(const AtomID &atomid) const
{
    return not this->selected(atomid);
}

/** Return whether none of the atoms in the CutGroup(s) 
    identified by 'cgid' have been selected 
    
    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
*/
bool AtomSelection::selectedNone(const CGID &cgid) const
{
    return not this->selected(cgid);
}

/** Return whether none of the atoms in the residue(s) 
    identified by 'resid' have been selected 
    
    \throw SireMol::missing_residue
    \throw SireError::invalid_index
*/
bool AtomSelection::selectedNone(const ResID &resid) const
{
    return not this->selected(resid);
}

/** Return whether none of the atoms in the chain(s) 
    identified by 'chainid' have been selected 
    
    \throw SireMol::missing_chain
    \throw SireError::invalid_index
*/
bool AtomSelection::selectedNone(const ChainID &chainid) const
{
    return not this->selected(chainid);
}

/** Return whether none of the atoms in the segment(s) 
    identified by 'segid' have been selected 
    
    \throw SireMol::missing_segment
    \throw SireError::invalid_index
*/
bool AtomSelection::selectedNone(const SegID &segid) const
{
    return not this->selected(segid);
}

/** Return whether none of the atoms selected in 'selection' have
    been selected in this set
    
    \throw SireError::incompatible_error
*/
bool AtomSelection::selectedNone(const AtomSelection &selection) const
{
    return not this->selected(selection);
}

/** Return whether or not all of the atoms are selected */
bool AtomSelection::selectedAll() const
{
    return nselected != 0 and selected_atoms.isEmpty();
}

/** Return whether or not the atom at index 'atomidx' is selected 

    \throw SireError::invalid_index
*/
bool AtomSelection::selectedAll(AtomIdx atomidx) const
{
    return this->selected(atomidx);
}

bool AtomSelection::_pvt_selectedAll(CGIdx cgidx) const
{
    if (nselected > 0)
    {
        //have all atoms been selected?
        if (selected_atoms.isEmpty())
            return true;
    
        const auto it = selected_atoms.constFind(cgidx);
        
        //have all atoms in the cutgroup been selected?
        if (it != selected_atoms.constEnd())
        {
            //all atoms selected if this list is empty (it wouldn't
            //exist if none of the atoms are selected)
            return it->isEmpty();
        }
    }

    return false;
}

/** Return whether or not all of the atoms in the CutGroup
    at index 'cgidx' have been selected
    
    \throw SireError::invalid_index
*/
bool AtomSelection::selectedAll(CGIdx cgidx) const
{
    return this->_pvt_selectedAll( CGIdx(cgidx.map(info().nCutGroups())) );
}

bool AtomSelection::_pvt_selectedAll(const QVector<CGAtomIdx> &atomidxs) const
{
    if (this->selectedAll())
        return true;
    else if (this->selectedNone() or nselected < atomidxs.count())
        return false;

    const CGAtomIdx *atomidxs_array = atomidxs.constData();
    int nats = atomidxs.count();

    for (int i=0; i<nats; ++i)
    {
        if (not this->_pvt_selected(atomidxs_array[i]))
            return false;
    }
    
    return true;
}

/** Return whether or not all of the atoms in the residue
    at index 'residx' have been selected 
    
    \throw SireError::invalid_index
*/
bool AtomSelection::selectedAll(ResIdx residx) const
{
    return this->selectedAll() or this->_pvt_selectedAll( info().cgAtomIdxs(residx) );
}

/** Return whether or not all of the atoms in the chain
    at index 'chainidx' have been selected 
    
    \throw SireError::invalid_index
*/
bool AtomSelection::selectedAll(ChainIdx chainidx) const
{
    return this->selectedAll() or this->_pvt_selectedAll( info().cgAtomIdxs(chainidx) );
}

/** Return whether or not all of the atoms in the segment
    at index 'segidx' have been selected
    
    \throw SireError::invalid_index
*/
bool AtomSelection::selectedAll(SegIdx segidx) const
{
    return this->selectedAll() or this->_pvt_selectedAll( info().cgAtomIdxs(segidx) );
}

/** Return whether or not all of the atoms matching the 
    ID 'atomid' have been selected
    
    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
bool AtomSelection::selectedAll(const AtomID &atomid) const
{
    return this->selectedAll() or this->_pvt_selectedAll( info().cgAtomIdxs(atomid) );
}

/** Return whether or not all of the atoms in the CutGroups matching the 
    ID 'cgid' have been selected
    
    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
*/
bool AtomSelection::selectedAll(const CGID &cgid) const
{
    foreach (CGIdx cgidx, cgid.map(this->info()))
    {
        if (not this->_pvt_selectedAll(cgidx))
            return false;
    }
    
    return true;
}

/** Return whether or not all of the atoms in the residues matching the 
    ID 'resid' have been selected
    
    \throw SireMol::missing_residue
    \throw SireError::invalid_index
*/
bool AtomSelection::selectedAll(const ResID &resid) const
{
    return this->_pvt_selectedAll( info().cgAtomIdxs(resid) );
}

/** Return whether or not all of the atoms in the chains matching the 
    ID 'atomid' have been selected
    
    \throw SireMol::missing_chain
    \throw SireError::invalid_index
*/
bool AtomSelection::selectedAll(const ChainID &chainid) const
{
    return this->_pvt_selectedAll( info().cgAtomIdxs(chainid) );
}

/** Return whether or not all of the atoms in the segments matching the 
    ID 'atomid' have been selected
    
    \throw SireMol::missing_segment
    \throw SireError::invalid_index
*/
bool AtomSelection::selectedAll(const SegID &segid) const
{
    return this->_pvt_selectedAll( info().cgAtomIdxs(segid) );
}

/** Return the set of indicies of the atoms in the CutGroup
    at index 'cgidx' that are selected within this CutGroup
    
    \throw SireError::invalid_index
*/
QSet<Index> AtomSelection::selectedAtoms(CGIdx cgidx) const
{
    cgidx = CGIdx( cgidx.map(info().nCutGroups()) );
    
    if (this->selectedAll())
    {
        //all atoms selected
        QSet<Index> ret;
        const int nats = info().nAtoms(cgidx);
        ret.reserve(nats);
        
        for (Index i(0); i<nats; ++i)
        {
            ret.insert(i);
        }
        
        return ret;
    }
    else
    {
        const auto it = selected_atoms.constFind(cgidx);
    
        if (it != selected_atoms.constEnd())
        {
            if (it->isEmpty())
            {
                //all atoms selected
                QSet<Index> ret;
                const int nats = info().nAtoms(cgidx);
                ret.reserve(nats);
                
                for (Index i(0); i<nats; ++i)
                {
                    ret.insert(i);
                }
                
                return ret;
            }
            else
                return *it;
        }
        else
            //the CutGroup is not selected
            return QSet<Index>();
    }
}

/** Return the list of indicies of CutGroups that contain at least
    one selected atom */
QList<CGIdx> AtomSelection::selectedCutGroups() const
{
    if (this->selectedAll())
    {
        return info().getCutGroups();
    }
    else
    {
        QList<CGIdx> keys = selected_atoms.keys();
        qSort(keys);
        return keys;
    }
}

/** Return the list of residues that contain at least one selected atom */
QList<ResIdx> AtomSelection::selectedResidues() const
{
    if (this->selectedAll())
    {
        return info().getResidues();
    }
    else if (this->selectedNone())
    {
        return QList<ResIdx>();
    }
    else
    {
        //run over all of the residues and see if they have
        //any selected atoms
        int nres = info().nResidues();
        QList<ResIdx> selected_res;
        
        for (ResIdx i(0); i<nres; ++i)
        {
            int nats = info().nAtoms(i);
            
            for (int j(0); j<nats; ++j)
            {
                if (this->selected( info().getAtom(i,j) ))
                {
                    selected_res.append(i);
                    break;
                }
            }
        }
        
        return selected_res;
    }
}

/** Return the list of chains that contain at least one selected atom */
QList<ChainIdx> AtomSelection::selectedChains() const
{
    if (this->selectedAll())
    {
        return info().getChains();
    }
    else if (this->selectedNone())
    {
        return QList<ChainIdx>();
    }
    else
    {
        //run over all of the chains and see if they have
        //any selected atoms
        int nchains = info().nChains();
        QList<ChainIdx> selected_chains;
        
        for (ChainIdx i(0); i<nchains; ++i)
        {
            bool has_atom = false;
            
            for (ResIdx j : info().getResiduesIn(i))
            {
                int nats = info().nAtoms(j);
            
                for (int k(0); k<nats; ++k)
                {
                    if (this->selected( info().getAtom(j,k) ))
                    {
                        has_atom = true;
                        break;
                    }
                }
                
                if (has_atom)
                    break;
            }
            
            if (has_atom)
                selected_chains.append(i);
        }
        
        return selected_chains;
    }
}

/** Return the list of segments that contain at least one selected atom */
QList<SegIdx> AtomSelection::selectedSegments() const
{
    if (this->selectedAll())
    {
        return info().getSegments();
    }
    else if (this->selectedNone())
    {
        return QList<SegIdx>();
    }
    else
    {
        //run over all of the segments and see if they have
        //any selected atoms
        int nseg = info().nSegments();
        QList<SegIdx> selected_seg;
        
        for (SegIdx i(0); i<nseg; ++i)
        {
            int nats = info().nAtoms(i);
            
            for (int j(0); j<nats; ++j)
            {
                if (this->selected( info().getAtom(i,j) ))
                {
                    selected_seg.append(i);
                    break;
                }
            }
        }
        
        return selected_seg;
    }
}

/** Return whether or not all of the atoms selected in 'selection' 
    have also been selected in this object
    
    \throw SireError::incompatible_error
*/
bool AtomSelection::selectedAll(const AtomSelection &selection) const
{
    this->info().assertEqualTo( selection.info() );
    
    //get rid of the easy cases...
    if (this->isEmpty() or selection.isEmpty())
        return false;
    else if (this->selectedAll())
        return true;
    else if (selection.selectedAll())
        return false;
    else if (selection.selectedAllCutGroups() and not this->selectedAllCutGroups())
        return false;
            
    //now do the hard stuff!
    if (selection.selectedAllCutGroups())
    {
        for (CGIdx i(0); i < this->info().nCutGroups(); ++i)
        {
            bool selected_all = this->selectedAll(i);
        
            if (selection.selectedAll(i) and not selected_all)
                return false;
            else if (selection.nSelected(i) > this->nSelected(i))
                return false;
            else if (not selected_all)
            {
                const QSet<Index> &atoms = *(selected_atoms.find(i));
            
                foreach (Index idx, selection.selectedAtoms(i))
                {
                    if (not atoms.contains(idx))
                        return false;
                }
            }
        }
    }
    else
    {
        QList<CGIdx> cgidxs = selection.selectedCutGroups();
        
        foreach (CGIdx i, cgidxs)
        {
            if (this->selectedNone(i))
                return false;
        
            bool selected_all = this->selectedAll(i);
            
            if (selection.selectedAll(i) and not selected_all)
                return false;
            else if (selection.nSelected(i) > this->nSelected(i))
                return false;
            else
            {
                const QSet<Index> &atoms = *(selected_atoms.find(i));
                
                foreach (Index idx, selection.selectedAtoms(i))
                {
                    if (not atoms.contains(idx))
                        return false;
                }
            }
        }
    }
    
    return true;
}

/** Return the number of atoms selected in the CutGroup at 
    index 'cgidx' */
int AtomSelection::nSelected(CGIdx cgidx) const
{
    if (this->selectedAll(cgidx))
        return info().nAtoms(cgidx);
    else
    {
        return selected_atoms.value( CGIdx(cgidx.map(info().nCutGroups())) ).count();
    }
}

/** Return whether the atom at index atomidx has been selected 

    \throw SireError::invalid_index
*/
int AtomSelection::nSelected(AtomIdx atomidx) const
{
    if (this->_pvt_selected( AtomIdx(atomidx.map(info().nAtoms())) ))
        return 1;
    else
        return 0;
}

int AtomSelection::_pvt_nSelected(ResIdx residx) const
{
    if (this->isEmpty())
        return 0;
    else if (this->selectedAll())
        return info().nAtoms(residx);

    int nats = 0;

    foreach (const CGAtomIdx &atomidx, info().cgAtomIdxs(residx))
    {
        if (this->_pvt_selected(atomidx))
            ++nats;
    }
    
    return nats;
}

/** Return the number of atoms from the residue at index 'residx'
    that have been selected 
    
    \throw SireError::invalid_index
*/
int AtomSelection::nSelected(ResIdx residx) const
{
    return this->_pvt_nSelected( ResIdx(residx.map(info().nResidues())) );
}

/** Return the number of atoms from the chain at index 'chainidx'
    that have been selected 
    
    \throw SireError::invalid_index
*/
int AtomSelection::nSelected(ChainIdx chainidx) const
{
    chainidx = ChainIdx( chainidx.map(info().nChains()) );
    
    if (this->isEmpty())
        return 0;
    else if (this->selectedAll())
        return info().nChains();
    
    int nats = 0;
    
    foreach (ResIdx residx, info().getResiduesIn(chainidx))
    {
        nats += this->_pvt_nSelected(residx);
    }
    
    return nats;
}

/** Return the number of atoms from the segment at index 'segidx'
    that have been selected 
    
    \throw SireError::invalid_index
*/
int AtomSelection::nSelected(SegIdx segidx) const
{
    segidx = SegIdx( segidx.map(info().nSegments()) );
    
    if (this->isEmpty())
        return 0;
    else if (this->selectedAll())
        return info().nAtoms(segidx);
     
    int nats = 0;
          
    foreach (const CGAtomIdx &atomidx, info().cgAtomIdxs(segidx))
    {
        if (this->_pvt_selected(atomidx))
            ++nats;
    }
    
    return nats;
}

/** Return the number of atoms from the CutGroups identified
    by 'cgid' that have been selected 
    
    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
*/
int AtomSelection::nSelected(const CGID &cgid) const
{
    int nats = 0;

    foreach (CGIdx cgidx, cgid.map(info()))
    {
        nats += this->nSelected(cgidx);
    }
    
    return nats;
}

/** Return the number of atoms that are identified by
    'atomid' that have been selected
    
    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
int AtomSelection::nSelected(const AtomID &atomid) const
{
    int nats = 0;
    
    foreach (const CGAtomIdx &atomidx, info().cgAtomIdxs(atomid))
    {
        if (this->_pvt_selected(atomidx))
            ++nats;
    }
    
    return nats;
}

/** Return the number of atoms from the residues identified
    by 'resid' that have been selected
    
    \throw SireMol::missing_residue
    \throw SireError::invalid_index
*/
int AtomSelection::nSelected(const ResID &resid) const
{
    int nats = 0;
    
    foreach (ResIdx residx, resid.map(info()))
    {
        nats += this->nSelected(residx);
    }
    
    return nats;
}

/** Return the number of atoms from the chain(s) identified
    by 'chainid' that have been selected
    
    \throw SireMol::missing_chain
    \throw SireError::invalid_index
*/
int AtomSelection::nSelected(const ChainID &chainid) const
{
    int nats = 0;
    
    foreach (ChainIdx chainidx, chainid.map(info()))
    {
        nats += this->nSelected(chainidx);
    }
    
    return nats;
}

/** Return the number of atoms from the segment(s) 
    identified by 'segid' that have been selected
    
    \throw SireMol::missing_segment
    \throw SireError::invalid_index
*/
int AtomSelection::nSelected(const SegID &segid) const
{
    int nats = 0;
    
    foreach (SegIdx segidx, segid.map(info()))
    {
        nats += this->nSelected(segidx);
    }
    
    return nats;
}

/** Return the number of atoms from the passed selection
    that have also been selected in this selection
    
    \throw SireError::incompatible_error
*/
int AtomSelection::nSelected(const AtomSelection &selection) const
{
    info().assertEqualTo(selection.info());
    
    if (selection.isEmpty() or this->isEmpty())
        return 0;
    else if (selection.selectedAll())
        return this->nSelected();
    else if (this->selectedAll())
        return selection.nSelected();
     
    int nats = 0;
              
    if (selection.selectedAllCutGroups())
    {
        for (CGIdx i(0); i < info().nCutGroups(); ++i)
        {
            if (this->selectedAll(i))
                nats += selection.nSelected(i);
            else if (selection.selectedAll(i))
                nats += this->nSelected(i);
            else if (not this->selectedNone(i))
            {
                const QSet<Index> &atoms = *(selected_atoms.find(i));
            
                foreach (Index idx, selection.selectedAtoms(i))
                {
                    if (atoms.contains(idx))
                        ++nats;
                }
            }
        }
    }
    else
    {
        foreach (CGIdx i, selection.selectedCutGroups())
        {
            if (this->selectedAll(i))
                nats += selection.nSelected(i);
            else if (selection.selectedAll(i))
                nats += this->nSelected(i);
            else if (not this->selectedNone(i))
            {
                const QSet<Index> &atoms = *(selected_atoms.find(i));
            
                foreach (Index idx, selection.selectedAtoms(i))
                {
                    if (atoms.contains(idx))
                        ++nats;
                }
            }
        }
    }
    
    return nats;
}

/** Return the total number of selected atoms */
int AtomSelection::nSelectedAtoms() const
{
    return nselected;
}

/** Return the number of CutGroups that have at least
    one selected atom */
int AtomSelection::nSelectedCutGroups() const
{
    if (this->selectedAll())
        return info().nCutGroups();
    else
    {
        return selected_atoms.count();
    }
}

/** Return the number of residues that contain at
    least one selected atom */
int AtomSelection::nSelectedResidues() const
{
    if (this->selectedAll())
        return info().nResidues();
    else if (this->selectedNone())
        return 0;
    else
    {
        int nres = 0;
        
        for (ResIdx i(0); i<info().nResidues(); ++i)
        {
            foreach (const CGAtomIdx &atomidx, info().cgAtomIdxs(i))
            {
                if (this->_pvt_selected(atomidx))
                {
                    ++nres;
                    break;
                }
            }
        }
        
        return nres;
    }
}

/** Return the number of chains that have at least one selected atom */
int AtomSelection::nSelectedChains() const
{
    if (this->selectedAll())
        return info().nChains();
    else if (this->selectedNone())
        return 0;
    else
    {
        int nchains = 0;
        
        for (ChainIdx i(0); i<info().nChains(); ++i)
        {
            if (this->selected(i))
                ++nchains;
        }
        
        return nchains;
    }
}

/** Return the number of segments that contain at 
    least one selected atom */
int AtomSelection::nSelectedSegments() const
{
    if (this->selectedAll())
        return info().nSegments();
    else if (this->selectedNone())
        return 0;
    else
    {
        int nsegs = 0;
        
        for (SegIdx i(0); i<info().nSegments(); ++i)
        {
            if (this->selected(i))
                ++nsegs;
        }
        
        return nsegs;
    }
}

/** Return the total number of atoms in the molecule */
int AtomSelection::nAtoms() const
{
    return info().nAtoms();
}

/** Return the total number of CutGroups in the molecule */
int AtomSelection::nCutGroups() const
{
    return info().nCutGroups();
}

/** Return the total number of residues in the molecule */
int AtomSelection::nResidues() const
{
    return info().nResidues();
}

/** Return the total number of chains in the molecule */
int AtomSelection::nChains() const
{
    return info().nChains();
}

/** Return the total number of segments in this molecule */
int AtomSelection::nSegments() const
{
    return info().nSegments();
}

/** Return a selection that has all of the atoms selected */
AtomSelection& AtomSelection::selectAll()
{
    selected_atoms.clear();
    nselected = info().nAtoms();
    
    return *this;
}

/** Return a selection that has none of the atoms selected */
AtomSelection& AtomSelection::deselectAll()
{
    selected_atoms.clear();
    nselected = 0;
    
    return *this;
}

/** Return a selection that has none of the atoms selected */
AtomSelection& AtomSelection::selectNone()
{
    return this->deselectAll();
}

void AtomSelection::_pvt_select(const CGAtomIdx &cgatomidx)
{
    if (this->_pvt_selected(cgatomidx))
        return;
    
    if (this->nSelected() == (info().nAtoms() - 1))
    {
        //this would result in selecting all atoms
        this->selectAll();
        return;
    }
    
    //if the cutgroup is already in selected_atoms then it has been selected before
    auto it = selected_atoms.find(cgatomidx.cutGroup());
    
    if (it == selected_atoms.end())
    {
        //this is the first atom selected in this CutGroup
        QSet<Index> cgatoms;
        cgatoms.insert(cgatomidx.atom());
        selected_atoms.insert(cgatomidx.cutGroup(), cgatoms);
    }
    else
    {
        if (it->count() == (info().nAtoms(cgatomidx.cutGroup()) - 1))
        {
            //we have now selected all atoms of this CutGroup - set the
            //set equal to an empty set, as this indicates that all atoms are here
            *it = QSet<Index>();
        }
        else
            it->insert(cgatomidx.atom());
    }
    
    ++nselected;
}

void AtomSelection::_pvt_select(const QVector<CGAtomIdx> &cgatomidxs)
{
    if (this->selectedAll())
        return;

    const CGAtomIdx *cgatomidxs_array = cgatomidxs.constData();
    int nats = cgatomidxs.count();
    
    for (int i=0; i<nats; ++i)
    {
        this->_pvt_select(cgatomidxs_array[i]);
    }
}

void AtomSelection::_pvt_select(CGIdx cgidx, const QSet<Index> &atoms)
{
    if (this->selectedAll())
        return;
        
    foreach (Index atom, atoms)
    {
        this->_pvt_select( CGAtomIdx(cgidx,atom) );
    }
}

void AtomSelection::_pvt_select(AtomIdx atomidx)
{
    this->_pvt_select( info().cgAtomIdx(atomidx) );
}

/** Select the atom at index 'atomidx' 

    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::select(AtomIdx atomidx)
{
    this->_pvt_select(atomidx);
    return *this;
}

void AtomSelection::_pvt_deselect(const CGAtomIdx &cgatomidx)
{
    if (not this->selected(cgatomidx))
        return;
    
    else if (nselected == 1)
    {
        //we have just removed the last selected atom
        selected_atoms.clear();
        nselected = 0;
        return;
    }
    else if ( this->selectedAll() )
    {
        //we need to create space for all CutGroups
        selected_atoms.clear();
        selected_atoms.reserve( info().nCutGroups() );
        
        for (CGIdx i(0); i<info().nCutGroups(); ++i)
        {
            if (i != cgatomidx.cutGroup())
            {
                selected_atoms.insert( i, QSet<Index>() );
            }
            else
            {
                int nats = info().nAtoms(i);
            
                QSet<Index> atoms;
                atoms.reserve(nats-1);
                
                for (Index j(0); j<nats; ++j)
                {
                    if (j != cgatomidx.atom())
                        atoms.insert(j);
                }
                
                selected_atoms.insert(i, atoms);
            }
        }
    }
    else if ( this->selectedAll(cgatomidx.cutGroup()) )
    {
        //we need to recreate space for this CutGroup
        int nats = info().nAtoms(cgatomidx.cutGroup());
        
        if (nats > 1)
        {
            QSet<Index> atoms;
            atoms.reserve(nats-1);
            
            for (Index i(0); i<nats; ++i)
            {
                if (i != cgatomidx.atom())
                    atoms.insert(i);
            }
            
            selected_atoms[cgatomidx.cutGroup()] = atoms;
        }
        else
        {
            //removed the one and only atom from the CutGroup
            selected_atoms.remove(cgatomidx.cutGroup());
        }
    }
    else
    {
        selected_atoms[cgatomidx.cutGroup()].remove(cgatomidx.atom());

        if (selected_atoms[cgatomidx.cutGroup()].isEmpty())
        {
            selected_atoms.remove(cgatomidx.cutGroup());
            
            if (selected_atoms.isEmpty() and nselected > 1)
            {
                throw SireError::program_bug( QObject::tr(
                    "Invalid state: nselected not 1 when removing last atom? %1, %2")
                        .arg(nselected).arg(cgatomidx.toString()), CODELOC );
            }
        }
    }
    
    --nselected;

    if (selected_atoms.isEmpty() and nselected > 0)
    {
        throw SireError::program_bug( QObject::tr(
            "Invalid state: nselected not zero when removed last atoms? %1, %2")
                .arg(nselected).arg(cgatomidx.toString()), CODELOC );
    }
    else if (nselected < 0)
    {
        throw SireError::program_bug( QObject::tr(
            "Invalid state: nselected has fallen below zero! %1, %2")
                .arg(nselected).arg(cgatomidx.toString()), CODELOC );
    }
}

void AtomSelection::_pvt_deselect(const QVector<CGAtomIdx> &cgatomidxs)
{
    if (this->selectedNone())
        return;

    const CGAtomIdx *cgatomidxs_array = cgatomidxs.constData();
    int nats = cgatomidxs.count();
    
    for (int i=0; i<nats; ++i)
    {
        this->_pvt_deselect( cgatomidxs_array[i] );
    }
}

void AtomSelection::_pvt_deselect(AtomIdx atomidx)
{
    this->_pvt_deselect( info().cgAtomIdx(atomidx) );
}

/** Deselect the atom at index 'atomidx' 

    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::deselect(AtomIdx atomidx)
{
    this->_pvt_deselect(atomidx);
    return *this;
}

/** Select only the atom at index 'atomidx'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::selectOnly(AtomIdx atomidx)
{
    info().assertContains(atomidx);
    return this->deselectAll().select(atomidx);
}

void AtomSelection::_pvt_select(CGIdx cgidx)
{
    if (this->selectedAll(cgidx))
        return;
    
    else if (this->selectedNone())
    {
        if (info().nCutGroups() == 1)
            nselected = info().nAtoms();
        else
        {
            selected_atoms.insert( cgidx, QSet<Index>() );
            nselected = info().nAtoms(cgidx);
        }
    }
    else
    {
        QSet<Index> atoms = selected_atoms.take(cgidx);

        int nadded = info().nAtoms(cgidx) - atoms.count();
        
        nselected += nadded;
        
        if (nselected == info().nAtoms())
        {
            //we have now selected all atoms
            selected_atoms.clear();
        }
        else
        {
            //we have selected all atoms in this CutGroup, which is indicated
            //with an empty set
            selected_atoms.insert(cgidx, QSet<Index>());
        }
    }
}

void AtomSelection::_pvt_deselect(CGIdx cgidx)
{
    //have we selected this CutGroup?
    if (this->selectedAll())
    {
        selected_atoms.clear();
        nselected -= info().nAtoms(cgidx);
        
        if (info().nCutGroups() > 1)
        {
            //we need to indicate that all other CutGroups are selected by adding
            //in extra empty sets
            selected_atoms.reserve( info().nCutGroups()-1 );
            
            for (CGIdx i(0); i<info().nCutGroups(); ++i)
            {
                if (i != cgidx)
                {
                    selected_atoms.insert( i, QSet<Index>() );
                }
            }
        }
    }
    else
    {
        //have we selected any of the atoms?
        if (selected_atoms.contains(cgidx))
        {
            auto atoms = selected_atoms.take(cgidx);
            
            if (atoms.isEmpty())
            {
                //the whole cutgroup was selected
                nselected -= info().nAtoms(cgidx);
            }
            else
            {
                //only part of the CutGroup was selected
                nselected -= atoms.count();
            }
        }
    }
    
    if (selected_atoms.isEmpty() and nselected > 0)
    {
        throw SireError::program_bug( QObject::tr(
            "Invalid state: nselected not zero when removed last atoms? %1, %2")
                .arg(nselected).arg(cgidx.toString()), CODELOC );
    }
    else if (nselected < 0)
    {
        throw SireError::program_bug( QObject::tr(
            "Invalid state: nselected has fallen below zero! %1, %2")
                .arg(nselected).arg(cgidx.toString()), CODELOC );
    }
}

/** Select the CutGroup at index 'cgidx'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::select(CGIdx cgidx)
{
    this->_pvt_select(cgidx);
    return *this;
}

/** Deselect the CutGroup at index 'cgidx'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::deselect(CGIdx cgidx)
{
    this->_pvt_deselect(cgidx);
    return *this;
}

/** Select only the CutGroup at index 'cgidx'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::selectOnly(CGIdx cgidx)
{
    info().assertContains(cgidx);
    return this->deselectAll().select(cgidx);
}

template<class IDXS>
void AtomSelection::_pvt_selectAtoms(const IDXS &atoms)
{
    foreach (const AtomIdx &atom, atoms)
    {
        this->_pvt_select(atom);
    }
}

template<class IDXS>
void AtomSelection::_pvt_deselectAtoms(const IDXS &atoms)
{
    foreach (const AtomIdx &atom, atoms)
    {
        this->_pvt_deselect(atom);
    }
}

/** Select the atoms in the residue at index 'residx'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::select(ResIdx residx)
{
    if (info().isResidueCutting(residx))
    {
        this->_pvt_select(info().cgIdx(residx));
    }
    else
    {
        this->_pvt_selectAtoms( info().getAtomsIn(residx) );
    }
    
    return *this;
}

/** Deselect the atoms in the residue at index 'residx'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::deselect(ResIdx residx)
{
    if (info().isResidueCutting(residx))
    {
        this->_pvt_deselect(info().cgIdx(residx));
    }
    else
    {
        this->_pvt_deselectAtoms( info().getAtomsIn(residx) );
    }
    
    return *this;
}

/** Select only the atoms in the residue at index 'residx'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::selectOnly(ResIdx residx)
{
    info().assertContains(residx);
    return this->deselectAll().select(residx);
}

/** Select all of the atoms in the chain at index 'chainidx'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::select(ChainIdx chainidx)
{
    foreach (const ResIdx residx, d->getResiduesIn(chainidx) )
    {
        this->select(residx);
    }
    
    return *this;
}

/** Deselect all of the atoms in the chain at index 'chainidx'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::deselect(ChainIdx chainidx)
{
    foreach (const ResIdx residx, d->getResiduesIn(chainidx))
    {
        this->deselect(residx);
    }
    
    return *this;
}

/** Select only the atoms that are in the chain at index 'chainidx'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::selectOnly(ChainIdx chainidx)
{
    info().assertContains(chainidx);
    return this->deselectAll().select(chainidx);
}

/** Select all of the atoms in the segment at index 'segidx'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::select(SegIdx segidx)
{
    this->_pvt_selectAtoms( d->getAtomsIn(segidx) );
    return *this;
}

/** Deselect all of the atoms in the segment at index 'segidx'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::deselect(SegIdx segidx)
{
    this->_pvt_deselectAtoms( d->getAtomsIn(segidx) );
    return *this;
}

/** Select only the atoms that are in the segment at index 'segidx'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::selectOnly(SegIdx segidx)
{
    info().assertContains(segidx);
    return this->deselectAll().select(segidx);
}

/** Deselect all of the atoms whose indicies are in 'atomidxs'

    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::select(const QSet<AtomIdx> &atomidxs)
{
    //check the indicies are sane!
    foreach (AtomIdx atomidx, atomidxs)
    {
        info().assertContains(atomidx);
    }

    this->_pvt_selectAtoms(atomidxs);
    
    return *this;
}

/** Deselect all of the atoms whose indicies are in 'atomidxs'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::deselect(const QSet<AtomIdx> &atomidxs)
{
    //check the indicies are sane
    foreach (AtomIdx atomidx, atomidxs)
    {
        info().assertContains(atomidx);
    }

    this->_pvt_deselectAtoms(atomidxs);
    
    return *this;
}

/** Select only the atoms whose indicies are in 'atomidxs'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::selectOnly(const QSet<AtomIdx> &atomidxs)
{
    //check the indicies are sane
    foreach (AtomIdx atomidx, atomidxs)
    {
        info().assertContains(atomidx);
    }

    this->deselectAll();
    this->_pvt_selectAtoms(atomidxs);
    
    return *this;
}

/** Deselect all of the atoms whose indicies are in 'atomidxs'

    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::select(const QList<AtomIdx> &atomidxs)
{
    //check that the indicies are sane
    foreach (AtomIdx atomidx, atomidxs)
    {
        info().assertContains(atomidx);
    }
    
    this->_pvt_selectAtoms(atomidxs);
    
    return *this;
}

/** Deselect all of the atoms whose indicies are in 'atomidxs'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::deselect(const QList<AtomIdx> &atomidxs)
{
    //check that the indicies are sane
    foreach (AtomIdx atomidx, atomidxs)
    {
        info().assertContains(atomidx);
    }

    this->_pvt_deselectAtoms(atomidxs);
    
    return *this;
}

/** Select only the atoms whose indicies are in 'atomidxs'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::selectOnly(const QList<AtomIdx> &atomidxs)
{
    //check that the indicies are sane
    foreach (AtomIdx atomidx, atomidxs)
    {
        info().assertContains(atomidx);
    }

    this->deselectAll();
    this->_pvt_selectAtoms(atomidxs);
    
    return *this;
}

/** Select all of the CutGroups whose indicies are in 'cgidxs'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::select(const QList<CGIdx> &cgidxs)
{
    //check that the indicies are sane
    foreach (CGIdx cgidx, cgidxs)
    {
        info().assertContains(cgidx);
    }
    
    int ncg = info().nCutGroups();
    
    foreach (const CGIdx &cgidx, cgidxs)
    {
        this->_pvt_select( CGIdx(cgidx.map(ncg)) );
    }
    
    return *this;
}

/** Deselect all of the CutGroups whose indicies are in 'cgidxs'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::deselect(const QList<CGIdx> &cgidxs)
{
    //check that the indicies are sane
    foreach (CGIdx cgidx, cgidxs)
    {
        info().assertContains(cgidx);
    }
    
    int ncg = info().nCutGroups();
    
    foreach (const CGIdx &cgidx, cgidxs)
    {
        this->_pvt_deselect( CGIdx(cgidx.map(ncg)) );
    }
    
    return *this;
}

/** Select only the atoms in the CutGroups whose indicies are in 'cgidxs'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::selectOnly(const QList<CGIdx> &cgidxs)
{
    //check that the indicies are sane
    foreach (CGIdx cgidx, cgidxs)
    {
        info().assertContains(cgidx);
    }
    
    int ncg = info().nCutGroups();
    this->deselectAll();
    
    foreach (const CGIdx &cgidx, cgidxs)
    {
        this->_pvt_select( CGIdx(cgidx.map(ncg)) );
    }
    
    return *this;
}

/** Select the atoms in the residues whose indicies are in 'residxs'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::select(const QList<ResIdx> &residxs)
{
    foreach (ResIdx residx, residxs)
    {
        info().assertContains(residx);
    }

    foreach (const ResIdx &residx, residxs)
    {
        this->_pvt_selectAtoms( d->getAtomsIn(residx) );
    }
    
    return *this;
}

/** Deselect all of the atoms that are in the residues whose
    indicies are in 'residxs'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::deselect(const QList<ResIdx> &residxs)
{
    foreach (ResIdx residx, residxs)
    {
        info().assertContains(residx);
    }

    foreach (const ResIdx &residx, residxs)
    {
        this->_pvt_deselectAtoms( d->getAtomsIn(residx) );
    }
    
    return *this;
}

/** Select only the atoms that in the residues whose indicies are
    in 'residxs'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::selectOnly(const QList<ResIdx> &residxs)
{
    foreach (ResIdx residx, residxs)
    {
        info().assertContains(residx);
    }
    
    this->selectNone();
    
    foreach (const ResIdx &residx, residxs)
    {
        this->_pvt_selectAtoms( d->getAtomsIn(residx) );
    }
    
    return *this;
}

/** Select the atoms that are in the chains whose indicies are in 'chainidxs'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::select(const QList<ChainIdx> &chainidxs)
{
    foreach (ChainIdx chainidx, chainidxs)
    {
        info().assertContains(chainidx);
    }
    
    foreach (const ChainIdx &chainidx, chainidxs)
    {
        foreach (const ResIdx &residx, d->getResiduesIn(chainidx))
        {
            this->_pvt_selectAtoms( d->getAtomsIn(residx) );
        }
    }
    
    return *this;
}

/** Deselect the atoms that are in the chains whose indicies
    are in 'chainidxs'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::deselect(const QList<ChainIdx> &chainidxs)
{
    foreach (ChainIdx chainidx, chainidxs)
    {
        info().assertContains(chainidx);
    }
    
    foreach (const ChainIdx &chainidx, chainidxs)
    {
        foreach (const ResIdx &residx, d->getResiduesIn(chainidx))
        {
            this->_pvt_deselectAtoms( d->getAtomsIn(residx) );
        }
    }
    
    return *this;
}

/** Select only the atoms that are in the chains whose indicies are
    in 'chainidxs'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::selectOnly(const QList<ChainIdx> &chainidxs)
{
    foreach (ChainIdx chainidx, chainidxs)
    {
        info().assertContains(chainidx);
    }

    this->deselectAll();
    
    foreach (const ChainIdx &chainidx, chainidxs)
    {
        foreach (const ResIdx &residx, d->getResiduesIn(chainidx))
        {
            this->_pvt_selectAtoms( d->getAtomsIn(residx) );
        }
    }
    
    return *this;
}

/** Select the atoms that are in the segments whose indicies
    are in 'segidxs'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::select(const QList<SegIdx> &segidxs)
{
    foreach (SegIdx segidx, segidxs)
    {
        info().assertContains(segidx);
    }
    
    foreach (const SegIdx &segidx, segidxs)
    {
        this->_pvt_selectAtoms( d->getAtomsIn(segidx) );
    }
    
    return *this;
}

/** Deselect all of the atoms in the segments whose indicies are
    in 'segidxs'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::deselect(const QList<SegIdx> &segidxs)
{
    foreach (SegIdx segidx, segidxs)
    {
        info().assertContains(segidx);
    }
    
    foreach (const SegIdx &segidx, segidxs)
    {
        this->_pvt_deselectAtoms( d->getAtomsIn(segidx) );
    }
    
    return *this;
}

/** Select only the atoms in the segments whose indicies are
    in 'segidxs'
        
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::selectOnly(const QList<SegIdx> &segidxs)
{
    foreach (SegIdx segidx, segidxs)
    {
        info().assertContains(segidx);
    }

    this->deselectAll();
    
    foreach (const SegIdx &segidx, segidxs)
    {
        this->_pvt_selectAtoms( d->getAtomsIn(segidx) );
    }
    
    return *this;
}

/** Select all of the atoms that match 'atomid'
    
    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::select(const AtomID &atomid)
{
    this->_pvt_select( info().cgAtomIdxs(atomid) );
    return *this;
}

/** Deselect all of the atoms that match 'atomid'
    
    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::deselect(const AtomID &atomid)
{
    this->_pvt_deselect( info().cgAtomIdxs(atomid) );
    
    return *this;
}

/** Select only that atoms that match 'atomid'
    
    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::selectOnly(const AtomID &atomid)
{
    QList<AtomIdx> atomidxs = atomid.map(info());
    
    this->deselectAll();
    
    foreach (AtomIdx atomidx, atomidxs)
    {
        this->_pvt_select(atomidx);
    }
    
    return *this;
}

/** Select the atoms in the CutGroups that match 'cgid'
    
    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::select(const CGID &cgid)
{
    this->_pvt_selectAtoms( d->getAtomsIn(cgid) );

    return *this;
}

/** Deselect the atoms in the CutGroups that match 'cgid'
    
    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::deselect(const CGID &cgid)
{
    this->_pvt_deselectAtoms( d->getAtomsIn(cgid) );
    
    return *this;
}

/** Select only the atoms that are in the CutGroups that
    match 'cgid'
    
    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::selectOnly(const CGID &cgid)
{
    QList<AtomIdx> atomidxs = info().getAtomsIn(cgid);
    
    this->deselectAll();
    this->_pvt_selectAtoms(atomidxs);
    
    return *this;
}

/** Select the atoms in the residues that match 'resid'
    
    \throw SireMol::missing_residue
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::select(const ResID &resid)
{
    this->_pvt_selectAtoms( d->getAtomsIn(resid) );

    return *this;
}

/** Deselect the atoms in the residues that match 'resid'
    
    \throw SireMol::missing_residue
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::deselect(const ResID &resid)
{
    this->_pvt_deselectAtoms( d->getAtomsIn(resid) );
    
    return *this;
}

/** Select only the atoms in the residues that match 'resid'
    
    \throw SireMol::missing_residue
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::selectOnly(const ResID &resid)
{
    QList<AtomIdx> atomidxs = info().getAtomsIn(resid);

    this->deselectAll();
    this->_pvt_selectAtoms(atomidxs);
    
    return *this;
}

/** Select the atoms in the chains that match 'chainid'
    
    \throw SireMol::missing_chain
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::select(const ChainID &chainid)
{
    this->_pvt_selectAtoms( d->getAtomsIn(chainid) );
    
    return *this;
}

/** Deselect the atoms in the chains that match 'chainid'
    
    \throw SireMol::missing_chain
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::deselect(const ChainID &chainid)
{
    this->_pvt_deselectAtoms( d->getAtomsIn(chainid) );
    
    return *this;
}

/** Select only the atoms in the chains that match 'chainid'
    
    \throw SireMol::missing_chain
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::selectOnly(const ChainID &chainid)
{
    QList<AtomIdx> atomidxs = info().getAtomsIn(chainid);
    
    this->deselectAll();
    this->_pvt_selectAtoms( d->getAtomsIn(chainid) );
    
    return *this;
}

/** Select all of the atoms in the segments that match 'segid'
    
    \throw SireMol::missing_segment
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::select(const SegID &segid)
{
    this->_pvt_selectAtoms( d->getAtomsIn(segid) );
    return *this;
}

/** Deselect all of the atoms in the segments that match 'segid'
    
    \throw SireMol::missing_segment
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::deselect(const SegID &segid)
{
    this->_pvt_deselectAtoms( d->getAtomsIn(segid) );
    
    return *this;
}

/** Select only the atoms in the segments that match 'segid'
    
    \throw SireMol::missing_segment
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::selectOnly(const SegID &segid)
{
    QList<AtomIdx> atomidxs = info().getAtomsIn(segid);

    this->deselectAll();
    this->_pvt_selectAtoms( d->getAtomsIn(segid) );
    
    return *this;
}

void AtomSelection::_pvt_select(const AtomSelection &selection)
{
    if (this == &selection)
        return;

    info().assertEqualTo(selection.info());
    
    if (this->selectedAll())
        return;
    else if (selection.selectedAll())
    {
        selected_atoms.clear();
        nselected = info().nAtoms();
        return;
    }
    
    if (selection.selectedAllCutGroups())
    {
        int ncg = info().nCutGroups();
        
        for (CGIdx i(0); i<ncg; ++i)
        {
            if (selection.selectedAll(i))
                this->_pvt_select(i);
            else if (not this->selectedAll(i))
            {
                foreach (Index idx, selection.selectedAtoms(i))
                {
                    this->_pvt_select( CGAtomIdx(i,idx) );
                }
            }
        }
    }
    else
    {
        foreach (CGIdx i, selection.selectedCutGroups())
        {
            if (selection.selectedAll(i))
                this->_pvt_select(i);
            else if (not this->selectedAll(i))
            {
                foreach(Index idx, selection.selectedAtoms(i))
                {
                    this->_pvt_select( CGAtomIdx(i,idx) );
                }
            }
        }
    }
}

/** Select all of the atoms in 'selection'
    
    \throw SireError::incompatible_error
*/
AtomSelection& AtomSelection::select(const AtomSelection &selection)
{
    this->_pvt_select(selection);
      
    return *this;
}

/** Deselect all of the atoms in 'selection'
    
    \throw SireError::incompatible_error
*/
AtomSelection& AtomSelection::deselect(const AtomSelection &selection)
{
    info().assertEqualTo(selection.info());
    
    if (this->selectedNone())
        return *this;
    else if (selection.selectedAll())
    {
        this->deselectAll();
        return *this;
    }
    
    if (selection.selectedAllCutGroups())
    {
        int ncg = info().nCutGroups();
        
        for (CGIdx i(0); i<ncg; ++i)
        {
            if (selection.selectedAll(i))
            {
                this->_pvt_deselect(i);
            }
            else if (not this->selectedNone(i))
            {
                foreach (Index idx, selection.selectedAtoms(i))
                {
                    this->_pvt_deselect( CGAtomIdx(i,idx) );
                }
            }
        }
    }
    else
    {
        foreach (CGIdx i, selection.selectedCutGroups())
        {
            if (selection.selectedAll(i))
            {
                this->_pvt_deselect(i);
            }
            else if (not this->selectedNone(i))
            {
                foreach (Index idx, selection.selectedAtoms(i))
                {
                    this->_pvt_deselect( CGAtomIdx(i,idx) );
                }
            }
        }
    }
    
    return *this;
}

/** Select only the atoms that are selected in 'selection'
    
    \throw SireError::incompatible_error
*/
AtomSelection& AtomSelection::selectOnly(const AtomSelection &selection)
{
    info().assertEqualTo(selection.info());
    
    selected_atoms = selection.selected_atoms;
    nselected = selection.nselected;
    
    return *this;
}

/** Invert this selection */
AtomSelection& AtomSelection::invert()
{
    if (this->selectedAll())
        return this->selectNone();
    else if (this->selectedNone())
        return this->selectAll();
        
    for (CGIdx i(0); i<info().nCutGroups(); ++i)
    {
        if (this->selectedAll(i))
            this->_pvt_deselect(i);
        else if (this->selectedNone(i))
            this->_pvt_select(i);
        else
        {
            QSet<Index> atoms = this->selectedAtoms(i);
            int nats = info().nAtoms(i);

            this->_pvt_deselect(i);
            
            for (Index j(0); j<nats; ++j)
            {
                if (not atoms.contains(j))
                    this->_pvt_select( CGAtomIdx(i,j) );
            }
        }
    }
    
    return *this;
}

/** Return whether or not this contains the atom at index 'atomidx'

    \throw SireError::invalid_index
*/
bool AtomSelection::intersects(AtomIdx atomidx) const
{
    return this->selected(atomidx);
}

/** Return whether or not the CutGroup at index 'cgidx' contains
    some atoms that have been selected 
    
    \throw SireError::invalid_index
*/
bool AtomSelection::intersects(CGIdx cgidx) const
{
    return this->selected(cgidx);
}

/** Return whether or not the residue at index 'residx' contains
    some atoms that have been selected 
    
    \throw SireError::invalid_index
*/
bool AtomSelection::intersects(ResIdx residx) const
{
    return this->selected(residx);
}

/** Return whether or not the chain at index 'chainidx' contains
    some atoms that have been selected 
    
    \throw SireError::invalid_index
*/
bool AtomSelection::intersects(ChainIdx chainidx) const
{
    return this->selected(chainidx);
}

/** Return whether or not the segment at index 'segidx' contains
    some atoms that have been selected 
    
    \throw SireError::invalid_index
*/
bool AtomSelection::intersects(SegIdx segidx) const
{
    return this->selected(segidx);
}

/** Return whether or not any of the atoms identified 
    by 'atomid' have been selected
    
    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
bool AtomSelection::intersects(const AtomID &atomid) const
{
    return this->selected(atomid);
}

/** Return whether or not any of the atoms in the CutGroup(s) identified 
    by 'cgid' have been selected
    
    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
*/
bool AtomSelection::intersects(const CGID &cgid) const
{
    return this->selected(cgid);
}

/** Return whether or not any of the atoms in the residue(s) identified 
    by 'resid' have been selected
    
    \throw SireMol::missing_residue
    \throw SireError::invalid_index
*/
bool AtomSelection::intersects(const ResID &resid) const
{
    return this->selected(resid);
}

/** Return whether or not any of the atoms in the chain(s) identified 
    by 'chainid' have been selected
    
    \throw SireMol::missing_chain
    \throw SireError::invalid_index
*/
bool AtomSelection::intersects(const ChainID &chainid) const
{
    return this->selected(chainid);
}

/** Return whether or not any of the atoms in the segment(s) identified 
    by 'segid' have been selected
    
    \throw SireMol::missing_segment
    \throw SireError::invalid_index
*/
bool AtomSelection::intersects(const SegID &segid) const
{
    return this->selected(segid);
}

/** Return whether any of the atoms selected in 'selection'
    are also selected in this set
    
    \throw SireError::incompatible_error
*/
bool AtomSelection::intersects(const AtomSelection &selection) const
{
    return this->selected(selection);
}

/** Return whether the atom at index 'atomidx' is selected

    \throw SireError::invalid_index
*/
bool AtomSelection::contains(AtomIdx atomidx) const
{
    return this->selected(atomidx);
}

/** Return whether all of the atoms in the CutGroup at 
    index 'cgidx' have been selected
    
    \throw SireError::invalid_index
*/
bool AtomSelection::contains(CGIdx cgidx) const
{
    return this->selectedAll(cgidx);
}

/** Return whether all of the atoms in the residue at 
    index 'residx' have been selected
    
    \throw SireError::invalid_index
*/
bool AtomSelection::contains(ResIdx residx) const
{
    return this->selectedAll(residx);
}

/** Return whether all of the atoms in the chain at 
    index 'chainidx' have been selected
    
    \throw SireError::invalid_index
*/
bool AtomSelection::contains(ChainIdx chainidx) const
{
    return this->selectedAll(chainidx);
}

/** Return whether all of the atoms in the segments at 
    index 'segidx' have been selected
    
    \throw SireError::invalid_index
*/
bool AtomSelection::contains(SegIdx segidx) const
{
    return this->selectedAll(segidx);
}

/** Return whether all of the atoms identified by
    'atomid' have been selected
    
    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
bool AtomSelection::contains(const AtomID &atomid) const
{
    return this->selectedAll(atomid);
}

/** Return whether all of the atoms in the CutGroup(s) identified by
    'cgid' have been selected
    
    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
*/
bool AtomSelection::contains(const CGID &cgid) const
{
    return this->selectedAll(cgid);
}

/** Return whether all of the atoms in the residue(s) identified by
    'resid' have been selected
    
    \throw SireMol::missing_residue
    \throw SireError::invalid_index
*/
bool AtomSelection::contains(const ResID &resid) const
{
    return this->selectedAll(resid);
}

/** Return whether all of the atoms in the chain(s) identified by
    'chainid' have been selected
    
    \throw SireMol::missing_chain
    \throw SireError::invalid_index
*/
bool AtomSelection::contains(const ChainID &chainid) const
{
    return this->selectedAll(chainid);
}

/** Return whether all of the atoms in the segment(s) identified by
    'segid' have been selected
    
    \throw SireMol::missing_segment
    \throw SireError::invalid_index
*/
bool AtomSelection::contains(const SegID &segid) const
{
    return this->selectedAll(segid);
}

/** Return whether all of the atoms selected in 'selection' are
    selected in this set
    
    \throw SireError::incompatible_error
*/
bool AtomSelection::contains(const AtomSelection &selection) const
{
    return this->selectedAll(selection);
}

/** Intersect this selection with 'other'

    \throw SireError::incompatible_error
*/
AtomSelection& AtomSelection::intersect(const AtomSelection &other)
{
    if (this == &other)
        return *this;

    info().assertEqualTo(other.info());

    if (this->selectedNone() or other.selectedNone())
    {
        return this->selectNone();
    }
    else if (this->selectedAll() and other.selectedAll())
    {
        return *this;
    }
    else
    {
        if (this->selectedAllCutGroups() and other.selectedAllCutGroups())
        {
            int ncg = info().nCutGroups();
            
            for (CGIdx i(0); i<ncg; ++i)
            {
                bool this_all = this->selectedAll(i);
                bool other_all = other.selectedAll(i);
            
                if (this_all and other_all)
                {
                    //nothing to do
                    continue;
                }
                else if (this_all)
                {
                    this->_pvt_deselect(i);
                    this->_pvt_select(i, other.selectedAtoms(i));
                }
                else if (other_all)
                {
                    //nothing to do
                    continue;
                }
                else
                {
                    QSet<Index> atoms = this->selectedAtoms(i);
                    atoms.intersect(other.selectedAtoms(i));

                    this->_pvt_deselect(i);
                    this->_pvt_select(i, atoms);
                } 
            }
        }
        else
        {
            QSet<CGIdx> cgidxs = this->selectedCutGroups().toSet();
            cgidxs.intersect( other.selectedCutGroups().toSet() );
            
            int ncg = info().nCutGroups();
            
            for (CGIdx i(0); i<ncg; ++i)
            {
                if (not cgidxs.contains(i))
                {
                    this->_pvt_deselect(i);
                    continue;
                }
            
                bool this_all = this->selectedAll(i);
                bool other_all = other.selectedAll(i);
            
                if (this_all and other_all)
                {
                    //nothing to do
                    continue;
                }
                else if (this_all)
                {
                    this->deselect(i);
                    this->_pvt_select(i, other.selectedAtoms(i));
                }
                else if (other_all)
                {
                    //nothing to do
                    continue;
                }
                else
                {
                    QSet<Index> atoms = this->selectedAtoms(i);
                    atoms.intersect(other.selectedAtoms(i));
                    
                    this->_pvt_deselect(i);
                    this->_pvt_select(i, atoms);
                } 
            }
        }
        
        return *this;
    }
}

/** Intersect this selection with the index 'atomidx'

    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::intersect(AtomIdx atomidx)
{
    if (this->selected(atomidx))
        return this->selectOnly(atomidx);
    else
        return this->deselectAll();
}

/** Intersect this selection with the index 'cgidx'

    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::intersect(CGIdx cgidx)
{
    if (not this->intersects(cgidx))
        return this->deselectAll();

    else
    {
        int ncg = info().nCutGroups();
        
        for (CGIdx i(0); i<ncg; ++i)
        {
            if (i != cgidx)
                this->_pvt_deselect(i);
        }
        
        return *this;
    }
}

/** Intersect this with the atoms whose indicies are in 'atomidxs'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::intersect(const QList<AtomIdx> &atomidxs)
{
    AtomSelection copy(*this);
    copy.selectOnly(atomidxs);

    return this->intersect(copy);
}

/** Intersect this with the atoms in the CutGroups whose indicies
    are in 'cgidxs'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::intersect(const QList<CGIdx> &cgidxs)
{
    if (cgidxs.isEmpty())
    {
        //nothing is selected
        this->deselectAll();
        return *this;
    }
    else if (cgidxs.count() == 1)
    {
        return this->intersect(cgidxs.first());
    }
    else
    {
        QSet<CGIdx> selected;
        selected.reserve(cgidxs.count());

        foreach (CGIdx idx, cgidxs)
        {
            selected.insert(idx);
        }

        int ncg = info().nCutGroups();

        for (CGIdx idx(0); idx < ncg; ++idx)
        {
            if (not selected.contains(idx))
            {
                this->_pvt_deselect(idx);
            }
        }

        return *this;
    }
}

/** Intersect this selection with the index 'residx'

    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::intersect(ResIdx residx)
{
    return this->intersect( info().getAtomsIn(residx) );
}

/** Intersect this selection with the index 'chainidx'

    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::intersect(ChainIdx chainidx)
{
    return this->intersect( info().getAtomsIn(chainidx) );
}

/** Intersect this selection with the index 'segidx' 

    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::intersect(SegIdx segidx)
{
    return this->intersect( info().getAtomsIn(segidx) );
}

/** Intersect this with the atoms in the residues whose indicies are
    in 'residxs'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::intersect(const QList<ResIdx> &residxs)
{
    AtomSelection copy(*this);
    copy.selectOnly(residxs);

    return this->intersect(copy);
}

/** Intersect this with the atoms in the chains whose indicies are
    in 'chainidxs'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::intersect(const QList<ChainIdx> &chainidxs)
{
    AtomSelection copy(*this);
    copy.selectOnly(chainidxs);
    return this->intersect(copy);
}

/** Intersect this with the atoms in the segments whose indicies are
    in 'segidxs'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::intersect(const QList<SegIdx> &segidxs)
{
    AtomSelection copy(*this);
    copy.selectOnly(segidxs);
    return this->intersect(copy);
}

/** Intersect this set with the atoms that match
    the ID 'atomid'
    
    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::intersect(const AtomID &atomid)
{
    return this->intersect( atomid.map(info()) );
}

/** Intersect this set with the atoms in the 
    CutGroup(s) that match the ID 'cgid'
    
    \throw SireMol::missing_cutgroup
*/
AtomSelection& AtomSelection::intersect(const CGID &cgid)
{
    return this->intersect( cgid.map(info()) );
}

/** Intersect this set with the atoms in the 
    residue(s) that match the ID 'resid'
    
    \throw SireMol::missing_residue
*/
AtomSelection& AtomSelection::intersect(const ResID &resid)
{
    return this->intersect( resid.map(info()) );
}

/** Intersect this set with the atoms in the 
    chain(s) that match the ID 'chainid'
    
    \throw SireMol::missing_chain
*/
AtomSelection& AtomSelection::intersect(const ChainID &chainid)
{
    return this->intersect( chainid.map(info()) );
}

/** Intersect this set with the atoms in the 
    segment(s) that match the ID 'segid'
    
    \throw SireMol::missing_segment
*/
AtomSelection& AtomSelection::intersect(const SegID &segid)
{
    return this->intersect( segid.map(info()) );
}

/** Unite this set with the atom at index 'atomidx'

    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::unite(AtomIdx atomidx)
{
    return this->select(atomidx);
}

/** Unite this set with the atoms in the   
    CutGroup at index 'cgidx'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::unite(CGIdx cgidx)
{
    return this->select(cgidx);
}

/** Unite this set with the atoms in the   
    residue at index 'residx'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::unite(ResIdx residx)
{
    return this->select(residx);
}

/** Unite this set with the atoms in the   
    chain at index 'chainidx'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::unite(ChainIdx chainidx)
{
    return this->select(chainidx);
}

/** Unite this set with the atoms in the   
    segment at index 'segidx'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::unite(SegIdx segidx)
{
    return this->select(segidx);
}

/** Unite this set with the atoms whose   
    indicies are in 'atomidxs'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::unite(const QList<AtomIdx> &atomidxs)
{
    return this->select(atomidxs);
}

/** Unite this set with the atoms in the   
    CutGroups whose indicies are in 'cgidxs'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::unite(const QList<CGIdx> &cgidxs)
{
    return this->select(cgidxs);
}

/** Unite this set with the atoms in the   
    residues whose indicies are in 'residxs'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::unite(const QList<ResIdx> &residxs)
{
    return this->select(residxs);
}

/** Unite this set with the atoms in the   
    chains whose indicies are in 'chainidxs'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::unite(const QList<ChainIdx> &chainidxs)
{
    return this->select(chainidxs);
}

/** Unite this set with the atoms in the   
    segments whose indicies are in 'segidxs'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::unite(const QList<SegIdx> &segidxs)
{
    return this->select(segidxs);
}

/** Unite this set with the atoms identified by
    'atomid'
    
    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::unite(const AtomID &atomid)
{
    return this->select(atomid);
}

/** Unite this set with the atoms in the 
    CutGroups identified by 'cgid'
    
    \throw SireMol::missing_cutgroup
*/
AtomSelection& AtomSelection::unite(const CGID &cgid)
{
    return this->select(cgid);
}

/** Unite this set with the atoms in the 
    residues identified by 'resid'
    
    \throw SireMol::missing_residue
*/
AtomSelection& AtomSelection::unite(const ResID &resid)
{
    return this->select(resid);
}

/** Unite this set with the atoms in the 
    chains identified by 'chainid'
    
    \throw SireMol::missing_chain
*/
AtomSelection& AtomSelection::unite(const ChainID &chainid)
{
    return this->select(chainid);
}

/** Unite this set with the atoms in the 
    segments identified by 'segid'
    
    \throw SireMol::missing_segment
*/
AtomSelection& AtomSelection::unite(const SegID &segid)
{
    return this->select(segid);
}

/** Unite this set with the atoms selected 
    in 'selection'
    
    \throw SireError::incompatible_error
*/
AtomSelection& AtomSelection::unite(const AtomSelection &selection)
{
    return this->select(selection);
}

/** Unite this set with all of the passed selections

    \throw SireError::incompatible_error
*/
AtomSelection& AtomSelection::unite(const QList<AtomSelection> &selections)
{
    //check sanity!
    foreach (const AtomSelection &selection, selections)
    {
        info().assertEqualTo( selection.info() );
    }

    foreach (const AtomSelection &selection, selections)
    {
        if (this->selectedAll())
            break;
            
        this->_pvt_select(selection);
    }
    
    return *this;
}

/** Subtract the atom at index 'atomidx' 
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::subtract(AtomIdx atomidx)
{
    return this->deselect(atomidx);
}

/** Subtract the atoms in the CutGroup
    at index 'cgidx'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::subtract(CGIdx cgidx)
{
    return this->deselect(cgidx);
}

/** Subtract the atoms in the residue
    at index 'residx'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::subtract(ResIdx residx)
{
    return this->deselect(residx);
}

/** Subtract the atoms in the chain
    at index 'chainidx'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::subtract(ChainIdx chainidx)
{
    return this->deselect(chainidx);
}

/** Subtract the atoms in the segment
    at index 'segidx' have been subtracted
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::subtract(SegIdx segidx)
{
    return this->deselect(segidx);
}

/** Subtract the atoms whose indicies
    are in 'atomidxs'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::subtract(const QList<AtomIdx> &atomidxs)
{
    return this->deselect(atomidxs);
}

/** Subtract the atoms that
    are in the CutGroups whose indicies are in 'cgidxs'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::subtract(const QList<CGIdx> &cgidxs)
{
    return this->deselect(cgidxs);
}

/** Subtract the atoms that
    are in the residues whose indicies are in 'residxs'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::subtract(const QList<ResIdx> &residxs)
{
    return this->deselect(residxs);
}

/** Subtract all of the atoms that
    are in the chains whose indicies are in 'chainidxs'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::subtract(const QList<ChainIdx> &chainidxs)
{
    return this->deselect(chainidxs);
}

/** Subtract all of the atoms that
    are in the segments whose indicies are in 'segidxs'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::subtract(const QList<SegIdx> &segidxs)
{
    return this->deselect(segidxs);
}

/** Subtract all of the atoms that 
    match the ID 'atomid'
    
    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::subtract(const AtomID &atomid)
{
    return this->deselect(atomid);
}

/** Subtract all of the atoms that
    are in the CutGroup(s) that match the ID 'cgid'
    
    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::subtract(const CGID &cgid)
{
    return this->deselect(cgid);
}

/** Subtract all of the atoms that
    are in the residue(s) that match the ID 'resid'
    
    \throw SireMol::missing_residue
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::subtract(const ResID &resid)
{
    return this->deselect(resid);
}

/** Subtract all of the atoms that
    are in the chain(s) that match the ID 'resid'
    
    \throw SireMol::missing_residue
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::subtract(const ChainID &chainid)
{
    return this->deselect(chainid);
}

/** Subtract all of the atoms that
    are in the segment(s) that match the ID 'segid'
    
    \throw SireMol::missing_segment
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::subtract(const SegID &segid)
{
    return this->deselect(segid);
}

/** Subtract all of the atoms selected
    in 'selection'
    
    \throw SireError::incompatible_error
*/
AtomSelection& AtomSelection::subtract(const AtomSelection &selection)
{
    return this->deselect(selection);
}

/** Mask this selection by the 
    atom at index 'atomidx'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::mask(AtomIdx atomidx)
{
    return this->intersect(atomidx);
}

/** Mask this selection by the 
    CutGroup at index 'cgidx'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::mask(CGIdx cgidx)
{
    return this->intersect(cgidx);
}

/** Mask this selection by the 
    residue at index 'residx' 
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::mask(ResIdx residx)
{
    return this->intersect(residx);
}

/** Mask this selection by the 
    chain at index 'chainidx' 
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::mask(ChainIdx chainidx)
{
    return this->intersect(chainidx);
}

/** Mask this selection by the 
    segment at index 'segidx' 
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::mask(SegIdx segidx)
{
    return this->intersect(segidx);
}

/** Mask this selection by the 
    indicies are in 'atomidxs'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::mask(const QList<AtomIdx> &atomidxs)
{
    return this->intersect(atomidxs);
}

/** Mask this selection by the 
    CutGroups whose indicies are in 'cgidxs'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::mask(const QList<CGIdx> &cgidxs)
{
    return this->intersect(cgidxs);
}

/** Mask this selection by the 
    residues whose indicies are in 'residxs'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::mask(const QList<ResIdx> &residxs)
{
    return this->intersect(residxs);
}

/** Mask this selection by the 
    chains whose indicies are in 'chainidxs'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::mask(const QList<ChainIdx> &chainidxs)
{
    return this->intersect(chainidxs);
}

/** Mask this selection by the 
    segments whose indicies are in 'segidxs'
    
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::mask(const QList<SegIdx> &segidxs)
{
    return this->intersect(segidxs);
}

/** Mask this selection by the 
    atoms that match the ID 'atomid'
    
    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::mask(const AtomID &atomid)
{
    return this->intersect(atomid);
}

/** Mask this selection by the 
    CutGroups that match the ID 'cgid'
    
    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::mask(const CGID &cgid)
{
    return this->intersect(cgid);
}

/** Mask this selection by the 
    residues that match the ID 'atomid'
    
    \throw SireMol::missing_residue
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::mask(const ResID &resid)
{
    return this->intersect(resid);
}

/** Mask this selection by the 
    chains that match the ID 'chainid'
    
    \throw SireMol::missing_chain
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::mask(const ChainID &chainid)
{
    return this->intersect(chainid);
}

/** Mask this selection by the 
    segments that match the ID 'atomid'
    
    \throw SireMol::missing_segments
    \throw SireError::invalid_index
*/
AtomSelection& AtomSelection::mask(const SegID &segid)
{
    return this->intersect(segid);
}

/** Mask this selection by 'other'

    \throw SireError::incompatible_error
*/
AtomSelection& AtomSelection::mask(const AtomSelection &other)
{
    return this->intersect(other);
}

/** Return the list of indicies of all of the atoms that
    have been selected */
QVector<AtomIdx> AtomSelection::selectedAtoms() const
{
    if (nselected == 0)
        return QVector<AtomIdx>();

    QVector<AtomIdx> ret(nselected);
    
    AtomIdx *ret_array = ret.data();
    
    if (this->selectedAll())
    {
        if (nselected != info().nAtoms())
        {
            throw SireError::program_bug( QObject::tr(
                "SERIOUS BUG: Disagreement of the number of atoms when selected all... "
                "%1 versus %2").arg(nselected).arg(info().nAtoms()), CODELOC );
        }
    
        for (AtomIdx i(0); i<info().nAtoms(); ++i)
        {
            ret_array[i] = i;
        }
    }
    else if (this->selectedAllCutGroups())
    {
        int count = 0;
        
        for (CGIdx i(0); i<info().nCutGroups(); ++i)
        {
            if (this->selectedAll(i))
            {
                for (Index j(0); j<info().nAtoms(i); ++i)
                {
                    ret_array[count] = info().atomIdx(CGAtomIdx(i,j));
                    ++count;
                }
            }
            else
            {
                foreach (Index j, this->selectedAtoms(i))
                {
                    ret_array[count] = info().atomIdx(CGAtomIdx(i,j));
                    ++count;
                }
            }
        }
    }
    else
    {
        int count = 0;
    
        foreach (CGIdx i, this->selectedCutGroups())
        {
            if (this->selectedAll(i))
            {
                for (Index j(0); j<info().nAtoms(i); ++j)
                {
                    ret_array[count] = info().atomIdx(CGAtomIdx(i,j));
                    ++count;
                }
            }
            else
            {
                foreach (Index j, this->selectedAtoms(i))
                {
                    ret_array[count] = info().atomIdx(CGAtomIdx(i,j));
                    ++count;
                }
            }
        }
    }

    qSort(ret);
    
    return ret;
}

/** Assert that this selection contains the atom at index 'atomidx'

    \throw SireError::invalid_index
*/
void AtomSelection::assertSelected(AtomIdx atomidx) const
{
    if (not this->contains(atomidx))
        throw SireError::invalid_index( QObject::tr(
            "This selection does not contain the atom at index %1.")
                .arg(atomidx), CODELOC );
}

/** Assert that this selection contains all of the atoms identified 
    by the ID 'atomid'
    
    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
void AtomSelection::assertSelected(const AtomID &atomid) const
{
    if (not this->contains(atomid))
        throw SireMol::missing_atom( QObject::tr(
            "This selection does not contain all of the atoms "
            "identified by the ID %1.")
                .arg(atomid.toString()), CODELOC );
}

/** Assert that this selection is compatible with the molecule info
    in 'molinfo'
    
    \throw SireError::incompatible_error
*/
void AtomSelection::assertCompatibleWith(const MoleculeInfoData &molinfo) const
{
    if (*d != molinfo)
        throw SireError::incompatible_error( QObject::tr(
            "The layout \"%1\" is incompatible with this selection, "
            "which is for the layout \"%2\".")
                .arg(molinfo.UID().toString()).arg(d->UID().toString()), CODELOC );
}

/** Assert that this selection is compatible with the molecule whose
    data is in 'moldata'
    
    \throw SireError::incompatible_error
*/
void AtomSelection::assertCompatibleWith(const MoleculeData &moldata) const
{
    if (*d != moldata.info())
        throw SireError::incompatible_error( QObject::tr(
            "The molecule \"%1\" (layout %2) is incompatible with this selection, "
            "which is for the layout \"%3\".")
                .arg(moldata.name()).arg(moldata.info().UID().toString())
                .arg(d->UID().toString()), CODELOC );
}

/** Return whether or not this selection is compatible with the molecule info
    in 'molinfo' */
bool AtomSelection::isCompatibleWith(const MoleculeInfoData &molinfo) const
{
    return molinfo == *d;
}

/** Assert that this selection is compatible with the molecule viewed
    in 'molview'
    
    \throw SireError::incompatible_error
*/
void AtomSelection::assertCompatibleWith(const MoleculeView &molview) const
{
    this->assertCompatibleWith(molview.data());
}

/** Assert that this selection is compatible with 'other'

    \throw SireError::incompatible_error
*/
void AtomSelection::assertCompatibleWith(const AtomSelection &other) const
{
    if (*d != other.info())
        throw SireError::incompatible_error( QObject::tr(
            "The selection for layout \"%1\" is incompatible with this selection, "
            "which is for the layout \"%2\".")
                .arg(other.info().UID().toString())
                .arg(d->UID().toString()), CODELOC );
}

AtomSelection& AtomSelection::select(const QSet<SireMol::ResIdx> &ids)
{
    return this->select(ids.toList());
}

AtomSelection& AtomSelection::selectOnly(const QSet<SireMol::ResIdx> &ids)
{
    return this->selectOnly(ids.toList());
}

AtomSelection& AtomSelection::selectOnly(const QSet<SireMol::CGIdx> &ids)
{
    return this->selectOnly(ids.toList());
}

AtomSelection& AtomSelection::intersect(const QSet<SireMol::CGIdx> &ids)
{
    return this->intersect(ids.toList());
}

AtomSelection& AtomSelection::subtract(const QSet<SireMol::ResIdx> &ids)
{
    return this->subtract(ids.toList());
}

AtomSelection& AtomSelection::subtract(const QSet<SireMol::AtomIdx> &ids)
{
    return this->subtract(ids.toList());
}

AtomSelection& AtomSelection::unite(const QSet<SireMol::SegIdx> &ids)
{
    return this->unite(ids.toList());
}

AtomSelection& AtomSelection::mask(const QSet<SireMol::AtomIdx> &ids)
{
    return this->select(ids.toList());
}

AtomSelection& AtomSelection::selectOnly(const QSet<SireMol::ChainIdx> &ids)
{
    return this->selectOnly(ids.toList());
}

AtomSelection& AtomSelection::select(const QSet<SireMol::ChainIdx> &ids)
{
    return this->select(ids.toList());
}

AtomSelection& AtomSelection::mask(const QSet<SireMol::SegIdx> &ids)
{
    return this->mask(ids.toList());
}

AtomSelection& AtomSelection::intersect(const QSet<SireMol::SegIdx> &ids)
{
    return this->intersect(ids.toList());
}

AtomSelection& AtomSelection::deselect(const QSet<SireMol::SegIdx> &ids)
{
    return this->deselect(ids.toList());
}

AtomSelection& AtomSelection::unite(const QSet<SireMol::ResIdx> &ids)
{
    return this->unite(ids.toList());
}

AtomSelection& AtomSelection::select(const QSet<SireMol::CGIdx> &ids)
{
    return this->select(ids.toList());
}

AtomSelection& AtomSelection::mask(const QSet<SireMol::ResIdx> &ids)
{
    return this->mask(ids.toList());
}

AtomSelection& AtomSelection::intersect(const QSet<SireMol::ResIdx> &ids)
{
    return this->intersect(ids.toList());
}

AtomSelection& AtomSelection::deselect(const QSet<SireMol::ResIdx> &ids)
{
    return this->deselect(ids.toList());
}

AtomSelection& AtomSelection::intersect(const QSet<SireMol::AtomIdx> &ids)
{
    return this->intersect(ids.toList());
}

AtomSelection& AtomSelection::intersect(const QSet<SireMol::ChainIdx> &ids)
{
    return this->intersect(ids.toList());
}

AtomSelection& AtomSelection::deselect(const QSet<SireMol::CGIdx> &ids)
{
    return this->deselect(ids.toList());
}

AtomSelection& AtomSelection::unite(const QSet<SireMol::ChainIdx> &ids)
{
    return this->unite(ids.toList());
}

AtomSelection& AtomSelection::select(const QSet<SireMol::SegIdx> &ids)
{
    return this->select(ids.toList());
}

AtomSelection& AtomSelection::selectOnly(const QSet<SireMol::SegIdx> &ids)
{
    return this->selectOnly(ids.toList());
}

AtomSelection& AtomSelection::subtract(const QSet<SireMol::CGIdx> &ids)
{
    return this->subtract(ids.toList());
}

AtomSelection& AtomSelection::deselect(const QSet<SireMol::ChainIdx> &ids)
{
    return this->deselect(ids.toList());
}

AtomSelection& AtomSelection::subtract(const QSet<SireMol::SegIdx> &ids)
{
    return this->subtract(ids.toList());
}

AtomSelection& AtomSelection::unite(const QSet<SireMol::CGIdx> &ids)
{
    return this->unite(ids.toList());
}

AtomSelection& AtomSelection::mask(const QSet<SireMol::ChainIdx> &ids)
{
    return this->mask(ids.toList());
}

AtomSelection& AtomSelection::mask(const QSet<SireMol::CGIdx> &ids)
{
    return this->mask(ids.toList());
}

AtomSelection& AtomSelection::subtract(const QSet<SireMol::ChainIdx> &ids)
{
    return this->subtract(ids.toList());
}

AtomSelection& AtomSelection::unite(const QSet<SireMol::AtomIdx> &ids)
{
    return this->unite(ids.toList());
}

const char* AtomSelection::typeName()
{
    return QMetaType::typeName( qMetaTypeId<AtomSelection>() );
}
