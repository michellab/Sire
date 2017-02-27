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

#include <QMutex>

#include "moleculeinfodata.h"
#include "structureeditor.h"
#include "atomselection.h"

#include "tostring.h"

#include "SireError/errors.h"
#include "SireMol/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <boost/tuple/tuple.hpp>

#include <QDebug>

using namespace SireMol;
using namespace SireID;
using namespace SireBase;
using namespace SireStream;

using boost::tuple;

//////////
////////// Implementation of MolInfo
//////////

MolInfo::MolInfo()
{}

MolInfo::~MolInfo()
{}

//////////
////////// Implementation of MolInfoRegistry
//////////

namespace SireMol
{
namespace detail
{

class ResInfo
{
public:
    ResInfo();
    ~ResInfo();
    
    bool operator==(const ResInfo &other) const;
    bool operator!=(const ResInfo &other) const;
    
    /** The name of this residue */
    ResName name;
    
    /** The number of this residue */
    ResNum number;

    /** The index of the chain this residue is in */
    ChainIdx chainidx;

    /** The list of the indicies of all atoms that are in this residue
        in the order they were added to this residue */
    QList<AtomIdx> atom_indicies;
    
    /** Hash mapping the name of each atom in this residue
        to the indicies of the atoms in the molecule
        (there may be more than one atom with the same name) */
    QMultiHash<QString,AtomIdx> atoms_by_name;
};

class ChainInfo
{
public:
    ChainInfo();
    ~ChainInfo();
    
    bool operator==(const ChainInfo &other) const;
    bool operator!=(const ChainInfo &other) const;
    
    /** The name of this chain */
    ChainName name;
    
    /** The list of all indicies of the residues that 
        are contained in this chain in the order they
        were added to the chain */
    QList<ResIdx> res_indicies;
};

class SegInfo
{
public:
    SegInfo();
    ~SegInfo();

    bool operator==(const SegInfo &other) const;
    bool operator!=(const SegInfo &other) const;

    /** The name of this segment */
    SegName name;
    
    /** The list of all indicies of the atoms that are in this segment
        in the order in which they were added to this segment */
    QList<AtomIdx> atom_indicies;
};

class CGInfo
{
public:
    CGInfo();
    ~CGInfo();
    
    bool operator==(const CGInfo &other) const;
    bool operator!=(const CGInfo &other) const;
    
    /** The name of this CutGroup */
    CGName name;
    
    /** The list of all indicies of all of the atoms 
        that are in this CutGroup in the order that they
        were added to this CutGroup */
    QList<AtomIdx> atom_indicies;
};

class AtomInfo
{
public:
    AtomInfo();
    ~AtomInfo();
    
    bool operator==(const AtomInfo &other) const;
    bool operator!=(const AtomInfo &other) const;
    
    /** The name of this atom */
    AtomName name;

    /** The number of this atom */
    AtomNum number;

    /** Index of the residue this atom is in
        (this also tells us the index of the chain
         this atom is in) */
    ResIdx residx;
    
    /** Index of the segment this atom is in */
    SegIdx segidx;
    
    /** The CGAtomIdx index of the atom in the molecule
         - this says which CutGroup the atom is in,
           which atom in that CutGroup this is */
    CGAtomIdx cgatomidx;
};

} // end of namespace detail
} // end of namespace SireMol

using namespace SireMol::detail;

////////
//////// Implementation of MoleculeInfoData::AtomInfo
////////

QDataStream& operator<<(QDataStream &ds, 
                        const AtomInfo &atominfo)
{
    ds << atominfo.name << atominfo.number
       << atominfo.residx << atominfo.segidx
       << atominfo.cgatomidx;
       
    return ds;
}

QDataStream& operator>>(QDataStream &ds,
                        AtomInfo &atominfo)
{
    ds >> atominfo.name >> atominfo.number
       >> atominfo.residx >> atominfo.segidx
       >> atominfo.cgatomidx;
       
    return ds;
}

AtomInfo::AtomInfo()
         : name( QString::null ), number( AtomNum::null() ),
           residx( ResIdx::null() ), segidx( SegIdx::null() )
{}

AtomInfo::~AtomInfo()
{}

bool AtomInfo::operator==(const AtomInfo &other) const
{
    return this == &other or
           (name == other.name and number == other.number and 
            residx == other.residx and segidx == other.segidx and
            cgatomidx == other.cgatomidx);
}

bool AtomInfo::operator!=(const AtomInfo &other) const
{
    return this != &other and
           (name != other.name or number != other.number or
            residx != other.residx or segidx != other.segidx or
            cgatomidx != other.cgatomidx);
}

////////
//////// Implementation of MoleculeInfoData::CGInfo
////////

QDataStream& operator<<(QDataStream &ds, 
                        const CGInfo &cginfo)
{
    ds << cginfo.name << cginfo.atom_indicies;
       
    return ds;
}

QDataStream& operator>>(QDataStream &ds,
                        CGInfo &cginfo)
{
    ds >> cginfo.name >> cginfo.atom_indicies;
       
    return ds;
}

CGInfo::CGInfo() : name( QString::null )
{}

CGInfo::~CGInfo()
{}

bool CGInfo::operator==(const CGInfo &other) const
{
    return name == other.name and atom_indicies == other.atom_indicies;
}

bool CGInfo::operator!=(const CGInfo &other) const
{
    return name != other.name or atom_indicies != other.atom_indicies;
}

////////
//////// Implementation of MoleculeInfoData::ResInfo;
////////

QDataStream& operator<<(QDataStream &ds, 
                        const ResInfo &resinfo)
{
    ds << resinfo.name << resinfo.number << resinfo.chainidx
       << resinfo.atom_indicies;
       
    return ds;
}

QDataStream& operator>>(QDataStream &ds,
                        ResInfo &resinfo)
{
    ds >> resinfo.name >> resinfo.number >> resinfo.chainidx
       >> resinfo.atom_indicies;
       
    return ds;
}

ResInfo::ResInfo()
        : name( QString::null ), number( ResNum::null() ),
          chainidx( ChainIdx::null() )
{}

ResInfo::~ResInfo()
{}

bool ResInfo::operator==(const ResInfo &other) const
{
    return this == &other or
           (name == other.name and number == other.number and
            chainidx == other.chainidx and atom_indicies == other.atom_indicies);
}

bool ResInfo::operator!=(const ResInfo &other) const
{
    return not this->operator==(other);
}

////////
//////// Implementation of MoleculeInfoData::ChainInfo
////////

QDataStream& operator<<(QDataStream &ds, 
                        const ChainInfo &chaininfo)
{
    ds << chaininfo.name << chaininfo.res_indicies;
       
    return ds;
}

QDataStream& operator>>(QDataStream &ds,
                        ChainInfo &chaininfo)
{
    ds >> chaininfo.name >> chaininfo.res_indicies;
       
    return ds;
}

ChainInfo::ChainInfo() : name( QString::null )
{}

ChainInfo::~ChainInfo()
{}

bool ChainInfo::operator==(const ChainInfo &other) const
{
    return name == other.name and res_indicies == other.res_indicies;
}

bool ChainInfo::operator!=(const ChainInfo &other) const
{
    return name != other.name or res_indicies != other.res_indicies;
}

////////
//////// Implementation of MoleculeInfoData::SegInfo
////////

QDataStream& operator<<(QDataStream &ds, 
                        const SegInfo &seginfo)
{
    ds << seginfo.name << seginfo.atom_indicies;
       
    return ds;
}

QDataStream& operator>>(QDataStream &ds,
                        SegInfo &seginfo)
{
    ds >> seginfo.name >> seginfo.atom_indicies;
       
    return ds;
}

SegInfo::SegInfo() : name( QString::null )
{}

SegInfo::~SegInfo()
{}

bool SegInfo::operator==(const SegInfo &other) const
{
    return name == other.name and atom_indicies == other.atom_indicies;
}

bool SegInfo::operator!=(const SegInfo &other) const
{
    return name != other.name or atom_indicies != other.atom_indicies;
}

////////
//////// Implementation of MoleculeInfoData
////////

static const RegisterMetaType<MoleculeInfoData> r_molinfo(NO_ROOT);

/** Serialise to a binary datastream */
QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds,
                                       const MoleculeInfoData &molinfo)
{
    writeHeader(ds, r_molinfo, 1);
    
    SharedDataStream sds(ds);
    
    sds << molinfo.uid
        << molinfo.atoms_by_index
        << molinfo.res_by_index << molinfo.chains_by_index
        << molinfo.seg_by_index << molinfo.cg_by_index;
        
    return ds;
}

void MoleculeInfoData::rebuildNameAndNumberIndexes()
{
    //now rebuild the other indexes
    atoms_by_name.clear();
    atoms_by_num.clear();
    
    int nats = atoms_by_index.count();
    
    if (nats > 0)
    {
        const AtomInfo *atoms_by_index_array = atoms_by_index.constData();
    
        for (AtomIdx i(0); i < nats; ++i)
        {
            const AtomInfo &atominfo = atoms_by_index_array[i];
        
            if (not atominfo.name.isNull())
                atoms_by_name.insert(atominfo.name, i);
            
            if (not atominfo.number.isNull())
                atoms_by_num.insert(atominfo.number, i);
        }
    }
    
    res_by_name.clear();
    res_by_num.clear();
    
    int nres = res_by_index.count();
    
    if (nres > 0)
    {
        ResInfo *res_by_index_array = res_by_index.data();
        const AtomInfo *atoms_by_index_array = atoms_by_index.constData();
    
        for (ResIdx i(0); i < nres; ++i)
        {
            ResInfo &resinfo = res_by_index_array[i];
    
            if (not resinfo.name.isNull())
                res_by_name.insert(resinfo.name, i);
        
            if (not resinfo.number.isNull())
                res_by_num.insert(resinfo.number, i);
            
            //index the names of all of the atoms in this residue
            foreach (AtomIdx atomidx, resinfo.atom_indicies)
            {
                BOOST_ASSERT( atomidx >= 0 and atomidx < nats );
            
                const AtomInfo &atominfo = atoms_by_index_array[atomidx];
            
                if (not atominfo.name.isNull())
                    resinfo.atoms_by_name.insert(atominfo.name, atomidx);
            }
            
            if (not resinfo.atoms_by_name.isEmpty())
                resinfo.atoms_by_name.squeeze();
        }
    }
 
    chains_by_name.clear();

    int nchains = chains_by_index.count();
    const ChainInfo *chains_by_index_array = chains_by_index.constData();
    
    if (nchains > 0)
    {
        for (ChainIdx i(0); i < nchains; ++i)
        {
            const ChainInfo &chaininfo = chains_by_index_array[i];
        
            if (not chaininfo.name.isNull())
                chains_by_name.insert(chaininfo.name, i);
        }
    }
    
    seg_by_name.clear();
    
    int nseg = seg_by_index.count();
    
    if (nseg > 0)
    {
        const SegInfo *seg_by_index_array = seg_by_index.constData();
    
        for (SegIdx i(0); i < nseg; ++i)
        {
            const SegInfo &seginfo = seg_by_index_array[i];
        
            if (not seginfo.name.isNull())
                seg_by_name.insert(seginfo.name, i);
        }
    }
    
    
    cg_by_name.clear();
    int ncg = cg_by_index.count();
    
    if (ncg > 0)
    {
        const CGInfo *cg_by_index_array = cg_by_index.constData();
        
        for (CGIdx i(0); i < ncg; ++i)
        {
            const CGInfo &cginfo = cg_by_index_array[i];
        
            if (not cginfo.name.isNull())
                cg_by_name.insert(cginfo.name, i);
        }
    }
    
    if (not atoms_by_name.isEmpty())
        atoms_by_name.squeeze();
    
    if (not atoms_by_num.isEmpty())
        atoms_by_num.squeeze();
    
    if (not res_by_name.isEmpty())
        res_by_name.squeeze();
    
    if (not res_by_num.isEmpty())
        res_by_num.squeeze();
    
    if (not chains_by_name.isEmpty())
        chains_by_name.squeeze();
    
    if (not seg_by_name.isEmpty())
        seg_by_name.squeeze();
    
    if (not cg_by_name.isEmpty())
        cg_by_name.squeeze();
}

/** Extract from a binary datastream */
QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds,
                                       MoleculeInfoData &molinfo)
{
    VersionID v = readHeader(ds, r_molinfo);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> molinfo.uid
            >> molinfo.atoms_by_index
            >> molinfo.res_by_index >> molinfo.chains_by_index
            >> molinfo.seg_by_index >> molinfo.cg_by_index;

        //reconstruct the name and number indexes
        molinfo.rebuildNameAndNumberIndexes();
    }
    else
        throw version_error( v, "1", r_molinfo, CODELOC );
        
    return ds;
}

/** Null constructor */
MoleculeInfoData::MoleculeInfoData() : RefCountData()
{}
    
const MoleculeInfoData& MoleculeInfoData::null()
{
    return *(create_shared_null<MoleculeInfoData>());
}
    
/** Construct from the passed StructureEditor */
MoleculeInfoData::MoleculeInfoData(const StructureEditor &editor)
                 : RefCountData()
{
    if (not editor.needsInfoRebuild())
    {
        this->operator=( editor.info() );
        return;
    }

    //first lets allocate memory for all of the parts of the molecule

    //first the atoms...
    int nats = editor.nAtomsInMolecule();
    
    if (nats == 0)
    {
        //this is an empty molecule!
        uid = 0;
        return;
    }
    
    atoms_by_index.resize(nats);
    atoms_by_index.squeeze();
    
    //now the CutGroups
    int ncg = editor.nCutGroupsInMolecule();
    
    if (ncg == 0)
    {
        //we cannot have a molecule with no CutGroups!
        throw SireMol::missing_cutgroup( QObject::tr(
            "You have not created any CutGroups in which to place the atoms. "
            "ALL atoms must be placed into a CutGroup!"), CODELOC );
    }
    
    cg_by_index.resize(ncg);
    cg_by_index.squeeze();
    
    //now the residues
    int nres = editor.nResiduesInMolecule();
    
    if (nres > 0)
    {
        res_by_index.resize(nres);
        res_by_index.squeeze();
    }
    
    //now the chains
    int nchains = editor.nChainsInMolecule();
    
    if (nchains > 0)
    {
        chains_by_index.resize(nchains);
        chains_by_index.squeeze();
    }
    
    //finally, the segments
    int nsegs = editor.nSegmentsInMolecule();
    
    if (nsegs > 0)
    {
        seg_by_index.resize(nsegs);
        seg_by_index.squeeze();
    }
    
    ///
    /// Now we must transfer the information from the editor
    /// to this object. Exceptions will be thrown if either
    /// there are atoms that are not part of any CutGroup,
    /// or if there are residues, cutgroups, chains or segments
    /// that are empty (don't contain any atoms)
    ///
    
    /// First, transfer the atom data
    ///
    AtomInfo* atoms_by_index_array = atoms_by_index.data();
    QList<AtomIdx> atoms_missing_cutgroups;
    
    for (AtomIdx i(0); i<nats; ++i)
    {
        //create the data for atom 'i'
        AtomInfo &atom = atoms_by_index_array[i];
        
        tuple<AtomName,AtomNum,CGAtomIdx,ResIdx,SegIdx> 
                                atomdata = editor.getAtomData(i);
                                
        atom.name = atomdata.get<0>();
        atom.number = atomdata.get<1>();
        atom.cgatomidx = atomdata.get<2>();
        atom.residx = atomdata.get<3>();
        atom.segidx = atomdata.get<4>();
        
        if (atom.cgatomidx.cutGroup().isNull())
            atoms_missing_cutgroups.append(i);
    }
    
    //ALL atoms must be placed into a CutGroup
    if (not atoms_missing_cutgroups.isEmpty())
        throw SireMol::missing_cutgroup( QObject::tr(
            "All atoms in the molecule MUST be placed into a CutGroup. "
            "Atoms missing a CutGroup are %1.")
                .arg( Sire::toString(atoms_missing_cutgroups) ), CODELOC );

    /// Now transfer the CutGroup data
    ///
    QList<CGIdx> empty_cutgroups;
    CGInfo *cg_by_index_array = cg_by_index.data();

    for (CGIdx i(0); i<ncg; ++i)
    {
        CGInfo &cutgroup = cg_by_index_array[i];
        
        tuple< CGName,QList<AtomIdx> > cgdata = editor.getCGData(i);
        
        cutgroup.name = cgdata.get<0>();
        cutgroup.atom_indicies = cgdata.get<1>();
        
        if (cutgroup.atom_indicies.isEmpty())
            empty_cutgroups.append(i);
    }
    
    //ALL CutGroups must contain at least one atom
    if (not empty_cutgroups.isEmpty())
        throw SireMol::missing_atom( QObject::tr(
            "ALL CutGroups must contain at least one atom. "
            "CutGroups missing atoms are %1.")
                .arg( Sire::toString(empty_cutgroups) ), CODELOC );

    /// Now transfer the residue data
    ///
    QList<ResIdx> empty_residues;
    ResInfo *res_by_index_array = res_by_index.data();
    
    for (ResIdx i(0); i<nres; ++i)
    {
        ResInfo &residue = res_by_index_array[i];
        
        tuple< ResName,ResNum,ChainIdx,QList<AtomIdx> >
                                resdata = editor.getResData(i);
                                
        residue.name = resdata.get<0>();
        residue.number = resdata.get<1>();
        residue.chainidx = resdata.get<2>();
        residue.atom_indicies = resdata.get<3>();
        
        if (residue.atom_indicies.isEmpty())
            empty_residues.append(i);
    }
    
    //ALL residues must contain at least one atom
    if (not empty_residues.isEmpty())
        throw SireMol::missing_atom( QObject::tr(
            "ALL residues must contain at least one atom. "
            "Residues missing atoms are %1.")
                .arg( Sire::toString(empty_residues) ), CODELOC );

    /// Now transfer the chain data
    ///
    QList<ChainIdx> empty_chains;
    ChainInfo *chains_by_index_array = chains_by_index.data();
    
    for (ChainIdx i(0); i<nchains; ++i)
    {
        ChainInfo &chain = chains_by_index_array[i];
        
        tuple< ChainName,QList<ResIdx> > chaindata = editor.getChainData(i);
        
        chain.name = chaindata.get<0>();
        chain.res_indicies = chaindata.get<1>();
        
        if (chain.res_indicies.isEmpty())
            empty_chains.append(i);
    }
    
    //ALL chains must contain at least one residue
    if (not empty_chains.isEmpty())
        throw SireMol::missing_residue( QObject::tr(
            "ALL chains must contain at least one residue. "
            "Chains missing residues are %1.")
                .arg( Sire::toString(empty_chains) ), CODELOC );
    
    /// Finally convert the segments
    ///
    QList<SegIdx> empty_segments;
    SegInfo *seg_by_index_array = seg_by_index.data();
    
    for (SegIdx i(0); i<nsegs; ++i)
    {
        SegInfo &segment = seg_by_index_array[i];
        
        tuple< SegName,QList<AtomIdx> > segdata = editor.getSegData(i);
        
        segment.name = segdata.get<0>();
        segment.atom_indicies = segdata.get<1>();
        
        if (segment.atom_indicies.isEmpty())
            empty_segments.append(i);
    }
    
    //ALL segments must contain at least one atom
    if (not empty_segments.isEmpty())
        throw SireMol::missing_atom( QObject::tr(
            "ALL segments must contain at least one atom. "
            "Segments missing atoms are %1.")
                .arg( Sire::toString(empty_segments) ), CODELOC );

    //now the conversion is complete, we finish by creating
    //the additional indexes that speed up searches
    this->rebuildNameAndNumberIndexes();
    
    //... and of course get a UID for this layout
    uid = QUuid::createUuid();
}
    
/** Copy constructor */
MoleculeInfoData::MoleculeInfoData(const MoleculeInfoData &other)
                 : RefCountData(),
                   uid(other.uid),
                   atoms_by_index(other.atoms_by_index),
                   atoms_by_name(other.atoms_by_name),
                   atoms_by_num(other.atoms_by_num),
                   res_by_index(other.res_by_index),
                   res_by_name(other.res_by_name),
                   res_by_num(other.res_by_num),
                   chains_by_index(other.chains_by_index),
                   chains_by_name(other.chains_by_name),
                   seg_by_index(other.seg_by_index),
                   seg_by_name(other.seg_by_name),
                   cg_by_index(other.cg_by_index),
                   cg_by_name(other.cg_by_name)
{}
  
/** Destructor */  
MoleculeInfoData::~MoleculeInfoData()
{}

/** Copy assignment operator */
MoleculeInfoData& MoleculeInfoData::operator=(const MoleculeInfoData &other)
{
    if (&other != this)
    {
        uid = other.uid;
        atoms_by_index = other.atoms_by_index;
        atoms_by_name = other.atoms_by_name;
        atoms_by_num = other.atoms_by_num;
        res_by_index = other.res_by_index;
        res_by_name = other.res_by_name;
        res_by_num = other.res_by_num;
        chains_by_index = other.chains_by_index;
        chains_by_name = other.chains_by_name;
        seg_by_index = other.seg_by_index;
        seg_by_name = other.seg_by_name;
        cg_by_index = other.cg_by_index;
        cg_by_name = other.cg_by_name;
    }
    
    return *this;
}

/** Assert that this molinfo to equal to 'other'

    \throw SireError::incompatible_error
*/
void MoleculeInfoData::assertEqualTo(const MoleculeInfoData &other) const
{
    if (this->operator!=(other))
        throw SireError::incompatible_error( QObject::tr(
            "The two passed molecule layouts are incompatible "
            "as they have different UIDs! (%1 vs. %2)")
                .arg(uid.toString()).arg(other.uid.toString()), CODELOC );
}

/** Use this function to minimise memory usage - this function
    compares the shared data in this info with 'other', and where
    they are equal it copies the data from 'other', thereby reducing
    wastage caused by duplicated storage
*/
void MoleculeInfoData::squeeze(const MoleculeInfoData &other) const
{
    if (this == &other)
        return;

    MoleculeInfoData *squeezed_this = const_cast<MoleculeInfoData*>(this);
    
    if (seg_by_index == other.seg_by_index)
    {
        squeezed_this->seg_by_index = other.seg_by_index;
        squeezed_this->seg_by_name = other.seg_by_name;
    }
        
    if (chains_by_index == other.chains_by_index)
    {
        squeezed_this->chains_by_index = other.chains_by_index;
        squeezed_this->chains_by_name = other.chains_by_name;
    }
        
    if (cg_by_index == other.cg_by_index)
    {
        squeezed_this->cg_by_index = other.cg_by_index;
        squeezed_this->cg_by_name = other.cg_by_name;
    }

    //ok, now only need to compare the atom and residue names
    int nres = res_by_index.count();

    const ResInfo *this_res_array = res_by_index.constData();
    const ResInfo *other_res_array = other.res_by_index.constData();

    if (this_res_array == other_res_array)
    {
        //both arrays are already identical!
        squeezed_this->res_by_index = other.res_by_index;
        squeezed_this->res_by_name = other.res_by_name;
        squeezed_this->res_by_num = other.res_by_num;
    }   
    else if (nres == other.res_by_index.count())
    {
        bool same_res_numbers = true;
        bool same_res_names = true;
        bool identical_residues = true;
        
        //make sure the residues have the same names and contain the
        //same atoms
        for (int i=0; i<nres; ++i)
        {
            const ResInfo &this_res = this_res_array[i];
            const ResInfo &other_res = other_res_array[i];
        
            same_res_names = same_res_names and
                             this_res.name == other_res.name;
                             
            same_res_numbers = same_res_numbers and
                               this_res.number == other_res.number;

            identical_residues = identical_residues and
                                 same_res_names and same_res_numbers and
                                 this_res.chainidx == other_res.chainidx and
                                 this_res.atom_indicies == other_res.atom_indicies and
                                 this_res.atoms_by_name == other_res.atoms_by_name;
            
            if ( not (same_res_names or same_res_numbers or identical_residues) )
                break;
        }
        
        if (identical_residues)
        {
            //have the same residue layout and names and numbers :-)
            squeezed_this->res_by_index = other.res_by_index;
            squeezed_this->res_by_name = other.res_by_name;
            squeezed_this->res_by_num = other.res_by_num;
        }
        else
        {
            if (same_res_names)
                //residues have the same names
                squeezed_this->res_by_name = other.res_by_name;
            
            if (same_res_numbers)
                //residues have the same numbers;
                squeezed_this->res_by_num = other.res_by_num;
        }
    }
        
    // now process the atoms
    int nats = atoms_by_index.count();

    const AtomInfo *this_atom_array = atoms_by_index.constData();
    const AtomInfo *other_atom_array = atoms_by_index.constData();

    if (this_atom_array == other_atom_array)
    {
        //both arrays are identical
        squeezed_this->atoms_by_index = other.atoms_by_index;
        squeezed_this->atoms_by_name = other.atoms_by_name;
        squeezed_this->atoms_by_num = other.atoms_by_num;
    }
    else if (nats == other.atoms_by_index.count())
    {
        bool same_atom_nums = true;
        bool same_atom_names = true;
        bool identical_atoms = true;
        
        for (int i=0; i<nats; ++i)
        {
            const AtomInfo &this_atom = this_atom_array[i];
            const AtomInfo &other_atom = other_atom_array[i];
            
            same_atom_nums = same_atom_nums and
                             this_atom.number == other_atom.number;
                             
            same_atom_names = same_atom_names and
                              this_atom.name == other_atom.name;
                              
            identical_atoms = identical_atoms and
                              same_atom_nums and same_atom_names and
                              this_atom.residx == other_atom.residx and
                              this_atom.segidx == other_atom.segidx and
                              this_atom.cgatomidx == other_atom.cgatomidx;
                              
            if ( not (same_atom_nums or same_atom_names or identical_atoms) )
                break;
        }
        
        if (identical_atoms)
        {
            squeezed_this->atoms_by_index = other.atoms_by_index;
            squeezed_this->atoms_by_name = other.atoms_by_name;
            squeezed_this->atoms_by_num = other.atoms_by_num;
        }
        else
        {
            if (same_atom_names)
                squeezed_this->atoms_by_name = other.atoms_by_name;
                
            if (same_atom_nums)
                squeezed_this->atoms_by_num = other.atoms_by_num;
        }
    }
}

/** Comparison operator - two molinfos are equal if they have the same UID */
bool MoleculeInfoData::operator==(const MoleculeInfoData &other) const
{
    return uid == other.uid;
}

/** Comparison operator - two molinfos are equal if they have the same UID */
bool MoleculeInfoData::operator!=(const MoleculeInfoData &other) const
{
    return uid != other.uid;
}

/** Return the unique ID number for this molecule layout */
const QUuid& MoleculeInfoData::UID() const
{
    return uid;
}

/** Return the name of the chain that matches the ID 'chainid' */
const ChainName& MoleculeInfoData::name(const ChainID &chainid) const
{
    return chains_by_index[ this->chainIdx(chainid) ].name;
}

/** Return the name of the chain that has the index 'chainidx' */
const ChainName& MoleculeInfoData::name(ChainIdx chainidx) const
{
    chainidx = chainidx.map(chains_by_index.count());
    return chains_by_index[chainidx].name;
}

/** Return the name of the segment with ID 'segid' */
const SegName& MoleculeInfoData::name(const SegID &segid) const
{
    return seg_by_index[ this->segIdx(segid) ].name;
}

/** Return the name of the segment with index 'segidx' */
const SegName& MoleculeInfoData::name(SegIdx segidx) const
{
    segidx = segidx.map(seg_by_index.count());
    return seg_by_index[segidx].name;
}

/** Return the name of the residue with ID 'resid' */
const ResName& MoleculeInfoData::name(const ResID &resid) const
{
    return res_by_index[ this->resIdx(resid) ].name;
}

/** Return the name of the residue with index 'residx' */
const ResName& MoleculeInfoData::name(ResIdx residx) const
{
    residx = residx.map(res_by_index.count());
    return res_by_index[residx].name;
}

/** Return the name of the CutGroup with ID 'cgid' */
const CGName& MoleculeInfoData::name(const CGID &cgid) const
{
    return cg_by_index[ this->cgIdx(cgid) ].name;
}

/** Return the name of the CutGroup with index 'cgidx' */
const CGName& MoleculeInfoData::name(CGIdx cgidx) const
{
    cgidx = cgidx.map(cg_by_index.count());
    return cg_by_index[cgidx].name;
}

/** Return the name of the atom with ID 'atomid' */
const AtomName& MoleculeInfoData::name(const AtomID &atomid) const
{
    return atoms_by_index[ this->atomIdx(atomid) ].name;
}

/** Return the name of the atom at index 'atomidx' */
const AtomName& MoleculeInfoData::name(AtomIdx atomidx) const
{
    atomidx = atomidx.map(atoms_by_index.count());
    return atoms_by_index[atomidx].name;
}

/** Return the number of the residue at ID 'resid' */
ResNum MoleculeInfoData::number(const ResID &resid) const
{
    return res_by_index[ this->resIdx(resid) ].number;
}

/** Return the number of the residue at index 'residx' */
ResNum MoleculeInfoData::number(ResIdx residx) const
{
    residx = residx.map(res_by_index.count());
    return res_by_index[residx].number;
}

/** Return the number of the atom with ID 'atomid' */
AtomNum MoleculeInfoData::number(const AtomID &atomid) const
{
    return atoms_by_index[ this->atomIdx(atomid) ].number;
}

/** Return the number of the atom at index 'atomidx' */
AtomNum MoleculeInfoData::number(AtomIdx atomidx) const
{
    atomidx = atomidx.map(atoms_by_index.count());
    return atoms_by_index[atomidx].number;
}

#if QT_VERSION < 0x040300
    /** This function provides the same functionality as the Qt 4.3 QMultiHash::remove(Key, Value) function */
    template<class Key, class T>
    void remove_from_hash(QMultiHash<Key,T> &hash, const Key &key, const T &value) 
    {
        //get all items with this key
        QList<T> values = hash.values(key);

        //remove all item with this value
        if (values.removeAll(value) > 0)
        {
             //some items were removed
             hash.remove(key);

             //put them back into the hash in the right order
             for (int i=values.count()-1; i>=0; ++i){ hash.insertMulti(key, values.at(i)); }
        }
    }
#endif

/** Rename the atom at index 'atomidx' to have the new name 'newname'. 
    
    \throw SireError::invalid_index
*/
MoleculeInfoData MoleculeInfoData::rename(AtomIdx atomidx, 
                                          const AtomName &newname) const
{
    atomidx = atomidx.map( this->nAtoms() );

    if (atoms_by_index.at(atomidx).name == newname)
        //the name hasn't changed
        return *this;

    //make the change in a copy of this object
    MoleculeInfoData newinfo(*this);
    
    AtomInfo &atom = newinfo.atoms_by_index[atomidx];

    #if QT_VERSION >= 0x040300
        newinfo.atoms_by_name.remove(atom.name, atomidx);
    #else
        ::remove_from_hash(newinfo.atoms_by_name, QString(atom.name), atomidx);
    #endif   
 
    if (not atom.residx.isNull())
    {
        ResInfo &residue = newinfo.res_by_index[atom.residx];

        #if QT_VERSION >= 0x040300
            residue.atoms_by_name.remove(atom.name, atomidx);
        #else
            ::remove_from_hash(residue.atoms_by_name, QString(atom.name), atomidx);
        #endif

        if (not newname.isNull())
            residue.atoms_by_name.insert(newname, atomidx);
    }

    if (not newname.isNull())
        newinfo.atoms_by_name.insert(newname, atomidx);

    atom.name = newname;
    
    //ok, we must now get a new uid for this info object
    newinfo.uid = QUuid::createUuid();
    
    return newinfo;
}

/** Renumber the atom at index 'atomidx' to 'newnum'. 
    
    \throw SireError::invalid_index
*/
MoleculeInfoData MoleculeInfoData::renumber(AtomIdx atomidx, 
                                            const AtomNum &newnum) const
{
    atomidx = atomidx.map( this->nAtoms() );
    
    if (atoms_by_index.at(atomidx).number == newnum)
        return *this;
        
    //make the change in a copy of this object
    MoleculeInfoData newinfo(*this);
    
    AtomInfo &atom = newinfo.atoms_by_index[atomidx];
    
    #if QT_VERSION >= 0x040300
        newinfo.atoms_by_num.remove(atom.number, atomidx);
    #else
        ::remove_from_hash(newinfo.atoms_by_num, atom.number, atomidx);
    #endif
    
    if (not newnum.isNull())
        newinfo.atoms_by_num.insert(newnum, atomidx);
    
    atom.number = newnum;
    
    newinfo.uid = QUuid::createUuid();
    
    return newinfo;
}

/** Rename the residue at index 'residx' with the name 'newname'. 

    \throw SireError::invalid_index
*/
MoleculeInfoData MoleculeInfoData::rename(ResIdx residx, 
                                          const ResName &newname) const
{
    residx = residx.map( this->nResidues() );
    
    if (res_by_index.at(residx).name == newname)
        return *this;
        
    //make the change in a copy of this object
    MoleculeInfoData newinfo(*this);
    
    ResInfo &residue = newinfo.res_by_index[residx];
    
    #if QT_VERSION >= 0x040300
        newinfo.res_by_name.remove(residue.name, residx);
    #else
        ::remove_from_hash(newinfo.res_by_name, QString(residue.name), residx);
    #endif

    if (not newname.isNull())
        newinfo.res_by_name.insert(newname, residx);
        
    residue.name = newname;
    
    //update the fingerprint
    newinfo.uid = QUuid::createUuid();
    
    return newinfo;
}

/** Renumber the residue at index 'residx' to 'newnum'.
    
    \throw SireError::invalid_index
*/
MoleculeInfoData MoleculeInfoData::renumber(ResIdx residx, 
                                            const ResNum &newnum) const
{
    residx = residx.map( this->nResidues() );
    
    if (res_by_index.at(residx).number == newnum)
        return *this;
        
    //make the change in a copy of this object
    MoleculeInfoData newinfo(*this);
    
    ResInfo &residue = newinfo.res_by_index[residx];
    
    #if QT_VERSION >= 0x040300
        newinfo.res_by_num.remove(residue.number, residx);
    #else
        ::remove_from_hash(newinfo.res_by_num, residue.number, residx);
    #endif
    
    if (not newnum.isNull())
        newinfo.res_by_num.insert(newnum, residx);
        
    residue.number = newnum;
    
    newinfo.uid = QUuid::createUuid();
    
    return newinfo;
}

/** Rename the CutGroup at index 'cgidx' to 'newname'.
    
    \throw SireError::invalid_index
*/
MoleculeInfoData MoleculeInfoData::rename(CGIdx cgidx, 
                                          const CGName &newname) const
{
    cgidx = cgidx.map( this->nCutGroups() );
    
    if (cg_by_index.at(cgidx).name == newname)
        return *this;
    
    //make the change in a copy
    MoleculeInfoData newinfo(*this);
    
    CGInfo &cgroup = newinfo.cg_by_index[cgidx];
    
    #if QT_VERSION >= 0x040300
        newinfo.cg_by_name.remove(cgroup.name, cgidx);
    #else
        ::remove_from_hash(newinfo.cg_by_name, QString(cgroup.name), cgidx);
    #endif
    
    if (not newname.isNull())
        newinfo.cg_by_name.insert(newname, cgidx);
        
    cgroup.name = newname;
    
    newinfo.uid = QUuid::createUuid();
    
    return newinfo;
}
                                          
/** Rename the chain at index 'chainidx' to 'newname'.
    
    \throw SireError::invalid_index
*/
MoleculeInfoData MoleculeInfoData::rename(ChainIdx chainidx, 
                                          const ChainName &newname) const
{
    chainidx = chainidx.map( this->nChains() );
    
    if (chains_by_index.at(chainidx).name == newname)
        return *this;
    
    //make the change in a copy
    MoleculeInfoData newinfo(*this);
    
    ChainInfo &chain = newinfo.chains_by_index[chainidx];
    
    #if QT_VERSION >= 0x040300
        newinfo.chains_by_name.remove(chain.name, chainidx);
    #else
        ::remove_from_hash(newinfo.chains_by_name, QString(chain.name), chainidx);
    #endif
    
    if (not newname.isNull())
        newinfo.chains_by_name.insert(newname, chainidx);
        
    chain.name = newname;
    
    newinfo.uid = QUuid::createUuid();
    
    return newinfo;
}

/** Rename the segment at index 'segidx' to 'newname'.
    
    \throw SireError::invalid_index
*/
MoleculeInfoData MoleculeInfoData::rename(SegIdx segidx, 
                                          const SegName &newname) const
{
    segidx = segidx.map( this->nSegments() );
    
    if (seg_by_index.at(segidx).name == newname)
        return *this;
    
    //make the change in a copy
    MoleculeInfoData newinfo(*this);
    
    SegInfo &segment = newinfo.seg_by_index[segidx];
    
    #if QT_VERSION >= 0x040300
        newinfo.seg_by_name.remove(segment.name, segidx);
    #else
        ::remove_from_hash(newinfo.seg_by_name, QString(segment.name), segidx);
    #endif
    
    if (not newname.isNull())
        newinfo.seg_by_name.insert(newname, segidx);
        
    segment.name = newname;
    
    newinfo.uid = QUuid::createUuid();
    
    return newinfo;
}

/** Return the CGAtomIdx of the atom at index 'atomidx' */
const CGAtomIdx& MoleculeInfoData::cgAtomIdx(AtomIdx atomidx) const
{
    atomidx = atomidx.map(atoms_by_index.count());
    return atoms_by_index.constData()[atomidx].cgatomidx;
}

/** Return the CGAtomIdx of the atom with ID 'atomid' */
const CGAtomIdx& MoleculeInfoData::cgAtomIdx(const AtomID &atomid) const
{
    return atoms_by_index[ this->atomIdx(atomid) ].cgatomidx;
}

/** Assert that there is the index for just one atom 
   
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
*/
void MolInfo::assertSingleAtom(const QList<AtomIdx> &atomidxs) const
{
    if (atomidxs.isEmpty())
        throw SireMol::missing_atom( QObject::tr(
            "Could not find a matching atom."), CODELOC );
                
    else if (atomidxs.count() > 1)
        throw SireMol::duplicate_atom( QObject::tr(
            "Too many atoms have matched (%1).")
                .arg( Sire::toString(atomidxs) ), CODELOC );
}

/** Assert that there is the index for just one residue
   
    \throw SireMol::missing_residue
    \throw SireMol::duplicate_residue
*/
void MolInfo::assertSingleResidue(const QList<ResIdx> &residxs) const
{
    if (residxs.isEmpty())
        throw SireMol::missing_residue( QObject::tr(
            "Could not find a matching residue."), CODELOC );
                
    else if (residxs.count() > 1)
        throw SireMol::duplicate_residue( QObject::tr(
            "Too many residues have matched (%1).")
                .arg( Sire::toString(residxs) ), CODELOC );
}

/** Assert that there is the index for just one chain 
   
    \throw SireMol::missing_chain
    \throw SireMol::duplicate_chain
*/
void MolInfo::assertSingleChain(const QList<ChainIdx> &chainidxs) const
{
    if (chainidxs.isEmpty())
        throw SireMol::missing_chain( QObject::tr(
            "Could not find a matching chain."), CODELOC );
                
    else if (chainidxs.count() > 1)
        throw SireMol::duplicate_chain( QObject::tr(
            "Too many chains have matched (%1).")
                .arg( Sire::toString(chainidxs) ), CODELOC );
}

/** Assert that there is the index for just one CutGroup 
   
    \throw SireMol::missing_cutgroup
    \throw SireMol::duplicate_cutgroup
*/
void MolInfo::assertSingleCutGroup(const QList<CGIdx> &cgidxs) const
{
    if (cgidxs.isEmpty())
        throw SireMol::missing_cutgroup( QObject::tr(
            "Could not find a matching CutGroup."), CODELOC );
                
    else if (cgidxs.count() > 1)
        throw SireMol::duplicate_cutgroup( QObject::tr(
            "Too many CutGroups have matched (%1).")
                .arg( Sire::toString(cgidxs) ), CODELOC );
}

/** Assert that there is the index for just one segment 
   
    \throw SireMol::missing_segment
    \throw SireMol::duplicate_segment
*/
void MolInfo::assertSingleSegment(const QList<SegIdx> &segidxs) const
{
    if (segidxs.isEmpty())
        throw SireMol::missing_segment( QObject::tr(
            "Could not find a matching segment."), CODELOC );
                
    else if (segidxs.count() > 1)
        throw SireMol::duplicate_segment( QObject::tr(
            "Too many segments have matched (%1).")
                .arg( Sire::toString(segidxs) ), CODELOC );
}

/** Return the index of the atom with ID 'atomid' */
AtomIdx MoleculeInfoData::atomIdx(const AtomID &atomid) const
{
    QList<AtomIdx> atomidxs = atomid.map(*this);
    
    this->assertSingleAtom(atomidxs);
    
    return atomidxs.first();
}

/** Return the index of the atom with CGAtomIdx 'cgatomidx' */
AtomIdx MoleculeInfoData::atomIdx(const CGAtomIdx &cgatomidx) const
{
    const QList<AtomIdx> &atomidxs 
            = cg_by_index[cgatomidx.cutGroup().map(cg_by_index.count())].atom_indicies;
    
    return atomidxs[cgatomidx.atom().map(atomidxs.count())];
}

/** Return the index of the residue that matches the ID 'resid' */
ResIdx MoleculeInfoData::resIdx(const ResID &resid) const
{
    QList<ResIdx> residxs = resid.map(*this);
    
    this->assertSingleResidue(residxs);
    
    return residxs.first();
}

/** Return the index of the chain that matches the ID 'chainid' */
ChainIdx MoleculeInfoData::chainIdx(const ChainID &chainid) const
{
    QList<ChainIdx> chainidxs = chainid.map(*this);
    
    this->assertSingleChain(chainidxs);
    
    return chainidxs.first();
}

/** Return the index of the segment that matches the ID 'segid' */
SegIdx MoleculeInfoData::segIdx(const SegID &segid) const
{
    QList<SegIdx> segidxs = segid.map(*this);
    
    this->assertSingleSegment(segidxs);
    
    return segidxs.first();
}

/** Return the index of the CutGroup that matches the ID 'cgid' */
CGIdx MoleculeInfoData::cgIdx(const CGID &cgid) const
{
    QList<CGIdx> cgidxs = cgid.map(*this);
    
    this->assertSingleCutGroup(cgidxs);
    
    return cgidxs.first();
}

/** Return the indicies of all of the segments in this molecule */
QList<SegIdx> MoleculeInfoData::getSegments() const
{
    QList<SegIdx> segidxs;
    
    for (int i=0; i<seg_by_index.count(); ++i)
    {
        segidxs.append(SegIdx(i));
    }
    
    return segidxs;
}

/** Return the indicies of all of the CutGroups in this molecule */
QList<CGIdx> MoleculeInfoData::getCutGroups() const
{
    QList<CGIdx> cgidxs;
    
    for (int i=0; i<cg_by_index.count(); ++i)
    {
        cgidxs.append(CGIdx(i));
    }
    
    return cgidxs;
}

/** Return the indicies of all of the chains in this molecule */
QList<ChainIdx> MoleculeInfoData::getChains() const
{
    QList<ChainIdx> chainidxs;
    
    for (int i=0; i<chains_by_index.count(); ++i)
    {
        chainidxs.append(ChainIdx(i));
    }
    
    return chainidxs;
}

/** Return the indicies of all of the residues in this molecule */
QList<ResIdx> MoleculeInfoData::getResidues() const
{
    QList<ResIdx> residxs;
    
    for (int i=0; i<res_by_index.count(); ++i)
    {
        residxs.append(ResIdx(i));
    }
    
    return residxs;
}

/** Return the indicies of all of the residues in the chain at index 'chainidx' */
const QList<ResIdx>& MoleculeInfoData::getResiduesIn(ChainIdx chainidx) const
{
    return chains_by_index[chainidx.map(chains_by_index.count())]
                    .res_indicies;
}

/** Return the indicies of all of the residues in the chains that match
    the ID 'chainid' */
QList<ResIdx> MoleculeInfoData::getResiduesIn(const ChainID &chainid) const
{
    QList<ChainIdx> chainidxs = chainid.map(*this);
    
    QList<ResIdx> residxs;
    
    foreach (ChainIdx chainidx, chainidxs)
    {
        residxs += chains_by_index[chainidx].res_indicies;
    }
    
    if (chainidxs.count() > 1)
        qSort(residxs);
        
    return residxs;
}

/** Return the index of the ith atom in the CutGroup at index 'cgidx'

    \throw SireError::invalid_index
*/
AtomIdx MoleculeInfoData::getAtom(CGIdx cgidx, int i) const
{
    const QList<AtomIdx> &atomidxs = this->getAtomsIn(cgidx);
    
    return atomidxs.at( Index(i).map(atomidxs.count()) );
}

/** Return the index of the ith atom in the residue at index 'residx'

    \throw SireError::invalid_index
*/
AtomIdx MoleculeInfoData::getAtom(ResIdx residx, int i) const
{
    const QList<AtomIdx> &atomidxs = this->getAtomsIn(residx);
    
    return atomidxs.at( Index(i).map(atomidxs.count()) );
}

/** Return the index of the ith atom in the chain at index 'chainidx'

    \throw SireError::invalid_index
*/
AtomIdx MoleculeInfoData::getAtom(ChainIdx chainidx, int i) const
{
    QList<AtomIdx> atomidxs = this->getAtomsIn(chainidx);
    
    return atomidxs.at( Index(i).map(atomidxs.count()) );
}

/** Return the index of the ith atom in the segment at index 'cgidx'

    \throw SireError::invalid_index
*/
AtomIdx MoleculeInfoData::getAtom(SegIdx segidx, int i) const
{
    const QList<AtomIdx> &atomidxs = this->getAtomsIn(segidx);
    
    return atomidxs.at( Index(i).map(atomidxs.count()) );
}

/** Return the index of the ith residue in the chain at index 'chainidx'

    \throw SireError::invalid_index
*/
ResIdx MoleculeInfoData::getResidue(ChainIdx chainidx, int i) const
{
    const QList<ResIdx> &residxs = this->getResiduesIn(chainidx);
    
    return residxs.at( Index(i).map(residxs.count()) );
}

/** Return the indicies of all of the atoms in this molecule */
QList<AtomIdx> MoleculeInfoData::getAtoms() const
{
    QList<AtomIdx> atomidxs;
    
    for (int i=0; i<atoms_by_index.count(); ++i)
    {
        atomidxs.append( AtomIdx(i) );
    }
    
    return atomidxs;
}

/** Return the indicies of all of the atoms in the residue at index 'residx' */
const QList<AtomIdx>& MoleculeInfoData::getAtomsIn(ResIdx residx) const
{
    return res_by_index[residx.map(res_by_index.count())].atom_indicies;
}

QList<AtomIdx> MoleculeInfoData::_pvt_getAtomsIn(const QList<ResIdx> &residxs) const
{
    QList<AtomIdx> atomidxs;
    
    foreach (ResIdx residx, residxs)
    {
        atomidxs += res_by_index[residx].atom_indicies;
    }
    
    if (residxs.count() > 1)
        qSort(atomidxs);
       
    return atomidxs;
}

/** Return the indicies of all of the atoms in the residues that match
    the ID 'resid' */
QList<AtomIdx> MoleculeInfoData::getAtomsIn(const ResID &resid) const
{
    QList<ResIdx> residxs = resid.map(*this);
    QList<AtomIdx> atomidxs = this->_pvt_getAtomsIn(residxs);
    
    BOOST_ASSERT( not atomidxs.isEmpty() );
    
    return atomidxs;
}

/** Return the indicies of the atoms in the residue at index 'residx'
    that are called 'name' */
QList<AtomIdx> MoleculeInfoData::getAtomsIn(ResIdx residx, 
                                            const AtomName &name) const
{
    QList<AtomIdx> atomidxs 
            = res_by_index[residx.map(res_by_index.count())].atoms_by_name.values(name);
    
    if (atomidxs.isEmpty())
        throw SireMol::missing_atom( QObject::tr(
            "There is no atom called \"%1\" in the residue at index %2 "
            "in the layout \"%3\".")
                .arg(name).arg(residx).arg(uid.toString()), CODELOC );
                
    if (atomidxs.count() > 1)
        qSort(atomidxs);
        
    return atomidxs;
}

QList<AtomIdx> MoleculeInfoData::_pvt_getAtomsIn(const QList<ResIdx> &residxs,
                                                 const AtomName &name) const
{
    QList<AtomIdx> atomidxs;
    
    foreach (ResIdx residx, residxs)
    {
        atomidxs += res_by_index[residx].atoms_by_name.values(name);
    }
                
    if (atomidxs.count() > 1)
        qSort(atomidxs);
    
    return atomidxs;
}

/** Return the indicies of the atoms in the residues that match the ID 'resid' 
    that are called 'name' */
QList<AtomIdx> MoleculeInfoData::getAtomsIn(const ResID &resid,
                                            const AtomName &name) const
{
    QList<ResIdx> residxs = resid.map(*this);

    QList<AtomIdx> atomidxs = this->_pvt_getAtomsIn(residxs,name);
    
    if (atomidxs.isEmpty())
        throw SireMol::missing_atom( QObject::tr(
            "There are no atoms called \"%1\" in the residues matching "
            "%2 in the layout \"%3\".")
                .arg(name, resid.toString()).arg(uid.toString()), CODELOC );

    return atomidxs;
}                                            

/** Return the list of atoms in the chain at index 'chainidx' */
QList<AtomIdx> MoleculeInfoData::getAtomsIn(ChainIdx chainidx) const
{
    QList<ResIdx> residxs = this->getResiduesIn(chainidx);
    QList<AtomIdx> atomidxs = this->_pvt_getAtomsIn(residxs);

    BOOST_ASSERT( not atomidxs.isEmpty() );
    
    return atomidxs;
}

/** Return the list of atoms in the chains identified by 'chainid' */
QList<AtomIdx> MoleculeInfoData::getAtomsIn(const ChainID &chainid) const
{
    QList<ResIdx> residxs = this->getResiduesIn(chainid);
    QList<AtomIdx> atomidxs = this->_pvt_getAtomsIn(residxs);

    BOOST_ASSERT( not atomidxs.isEmpty() );
    
    return atomidxs;
}

/** Return the list of atoms in the chain at index 'chainidx' that
    are also called 'name' */
QList<AtomIdx> MoleculeInfoData::getAtomsIn(ChainIdx chainidx, 
                                            const AtomName &name) const
{
    QList<ResIdx> residxs = this->getResiduesIn(chainidx);
    QList<AtomIdx> atomidxs = this->_pvt_getAtomsIn(residxs, name);
    
    if (atomidxs.isEmpty())
        throw SireMol::missing_atom( QObject::tr(
            "There are no atoms called \"%1\" in the chain at index "
            "%2 in the layout \"%3\".")
                .arg(name).arg(chainidx).arg(uid.toString()), CODELOC );
                
    return atomidxs;
}
                                            
/** Return the list of atoms in the chains identified by the ID 'chainid'
    that are also called 'name' */
QList<AtomIdx> MoleculeInfoData::getAtomsIn(const ChainID &chainid,
                                            const AtomName &name) const
{
    QList<ResIdx> residxs = this->getResiduesIn(chainid);
    QList<AtomIdx> atomidxs = this->_pvt_getAtomsIn(residxs, name);
    
    if (atomidxs.isEmpty())
        throw SireMol::missing_atom( QObject::tr(
            "There are no atoms called \"%1\" in the chains matching "
            "%2 in the layout \"%3\".")
                .arg(name, chainid.toString()).arg(uid.toString()), CODELOC );
    
    return atomidxs;
}
                          
/** Return the indicies of atoms in the CutGroup at index 'cgidx' */
const QList<AtomIdx>& MoleculeInfoData::getAtomsIn(CGIdx cgidx) const
{
    return cg_by_index[cgidx.map(cg_by_index.count())].atom_indicies;
}

/** Return the indicies of atoms in the CutGroups identified by ID 'cgid' */
QList<AtomIdx> MoleculeInfoData::getAtomsIn(const CGID &cgid) const
{
    QList<CGIdx> cgidxs = cgid.map(*this);
    
    QList<AtomIdx> atomidxs;
    
    foreach (CGIdx cgidx, cgidxs)
    {
        atomidxs += cg_by_index[cgidx].atom_indicies;
    }
    
    BOOST_ASSERT( not atomidxs.isEmpty() );
    
    if (cgidxs.count() > 1)
        qSort(atomidxs);
        
    return atomidxs;
}

/** Return the indicies of the atoms in the segment at index 'segidx' */
const QList<AtomIdx>& MoleculeInfoData::getAtomsIn(SegIdx segidx) const
{
    return seg_by_index[segidx.map(seg_by_index.count())].atom_indicies;
}

/** Return the indicies of atoms in the segments that match the ID 'segid' */
QList<AtomIdx> MoleculeInfoData::getAtomsIn(const SegID &segid) const
{
    QList<SegIdx> segidxs = segid.map(*this);
    
    QList<AtomIdx> atomidxs;
    
    foreach (SegIdx segidx, segidxs)
    {
        atomidxs += seg_by_index[segidx].atom_indicies;
    }
    
    BOOST_ASSERT( not atomidxs.isEmpty() );
    
    if (segidxs.count() > 1)
        qSort(atomidxs);
        
    return atomidxs;
}

QVector<CGAtomIdx> 
MoleculeInfoData::_pvt_cgAtomIdxs(const QList<AtomIdx> &atomidxs) const
{
    int nats = atomidxs.count();
    
    QVector<CGAtomIdx> cgatomidxs(nats);
    
    CGAtomIdx *cgatomidxs_array = cgatomidxs.data();
    const AtomInfo *atoms_by_index_array = atoms_by_index.constData();
    
    for (int i=0; i<nats; ++i)
    {
        cgatomidxs_array[i] = atoms_by_index_array[atomidxs[i]].cgatomidx;
    }
    
    return cgatomidxs;
}

/** Simple function that returns the CGAtomIdx of the atom at index 'atomidx'

    \throw SireError::invalid_index
*/
QVector<CGAtomIdx> MoleculeInfoData::cgAtomIdxs(AtomIdx atomidx) const
{
    return QVector<CGAtomIdx>(1, this->cgAtomIdx(atomidx));
}

/** Return the CGAtomIdxs of all of the atoms that match the ID 'atomid' */
QVector<CGAtomIdx> MoleculeInfoData::cgAtomIdxs(const AtomID &atomid) const
{
    return this->_pvt_cgAtomIdxs( atomid.map(*this) );
}

/** Return the CGAtomIdxs of all of the atoms in the CutGroup
    at index 'cgidx'
    
    \throw SireError::invalid_index
*/
QVector<CGAtomIdx> MoleculeInfoData::cgAtomIdxs(CGIdx cgidx) const
{
    return this->_pvt_cgAtomIdxs( this->getAtomsIn(cgidx) );
}

/** Return the CGAtomIdxs of all of the atoms in the residue
    at index 'residx'
    
    \throw SireError::invalid_index
*/
QVector<CGAtomIdx> MoleculeInfoData::cgAtomIdxs(ResIdx residx) const
{
    return this->_pvt_cgAtomIdxs( this->getAtomsIn(residx) );
}

/** Return the CGAtomIdxs of all of the atoms in the chain 
    at index 'chainidx'
    
    \throw SireError::invalid_index
*/
QVector<CGAtomIdx> MoleculeInfoData::cgAtomIdxs(ChainIdx chainidx) const
{
    return this->_pvt_cgAtomIdxs( this->getAtomsIn(chainidx) );
}

/** Return the CGAtomIdxs of all of the atoms in the segment
    at index 'segidx'
    
    \throw SireError::invalid_index
*/
QVector<CGAtomIdx> MoleculeInfoData::cgAtomIdxs(SegIdx segidx) const
{
    return this->_pvt_cgAtomIdxs( this->getAtomsIn(segidx) );
}

/** Return the CGAtomIdxs of all of the atoms in the CutGroup(s)
    identified by 'cgid'
    
    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
*/
QVector<CGAtomIdx> MoleculeInfoData::cgAtomIdxs(const CGID &cgid) const
{
    return this->_pvt_cgAtomIdxs( this->getAtomsIn(cgid) );
}

/** Return the CGAtomIdxs of all of the atoms in the residues(s)
    identified by 'resid'
    
    \throw SireMol::missing_residue
    \throw SireError::invalid_index
*/
QVector<CGAtomIdx> MoleculeInfoData::cgAtomIdxs(const ResID &resid) const
{
    return this->_pvt_cgAtomIdxs( this->getAtomsIn(resid) );
}

/** Return the CGAtomIdxs of all of the atoms in the chain(s)
    identified by 'chainid'
    
    \throw SireMol::missing_chain
    \throw SireError::invalid_index
*/
QVector<CGAtomIdx> MoleculeInfoData::cgAtomIdxs(const ChainID &chainid) const
{
    return this->_pvt_cgAtomIdxs( this->getAtomsIn(chainid) );
}

/** Return the CGAtomIdxs of all of the atoms in the segment(s)
    identified by 'segid'
    
    \throw SireMol::missing_segment
    \throw SireError::invalid_index
*/
QVector<CGAtomIdx> MoleculeInfoData::cgAtomIdxs(const SegID &segid) const
{
    return this->_pvt_cgAtomIdxs( this->getAtomsIn(segid) );
}

/** Return whether or not the atom with the specified index is
    part of a residue
    
    \throw SireError::invalid_index
*/
bool MoleculeInfoData::isWithinResidue(AtomIdx atomidx) const
{
    ResIdx residx = atoms_by_index[atomidx.map(atoms_by_index.count())].residx;
    return not residx.isNull();
}

/** Return whether or not the atom identified by the ID 'atomid' is
    part of a residue
    
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
bool MoleculeInfoData::isWithinResidue(const AtomID &atomid) const
{
    ResIdx residx = atoms_by_index[ this->atomIdx(atomid) ].residx;
    return not residx.isNull();
}

/** Return whether or not the residue at index 'residx' is 
    part of a chain
    
    \throw SireError::invalid_index
*/
bool MoleculeInfoData::isWithinChain(ResIdx residx) const
{
    ChainIdx chainidx = res_by_index[residx.map(res_by_index.count())].chainidx;
    return not chainidx.isNull();
}

/** Return whether or not the residue identified by the ID 'resid'
    is part of a chain
    
    \throw SireMol::missing_residue
    \throw SireError::duplicate_residue
    \throw SireError::invalid_index
*/
bool MoleculeInfoData::isWithinChain(const ResID &resid) const
{
    ChainIdx chainidx = res_by_index[ this->resIdx(resid) ].chainidx;
    return not chainidx.isNull();
}

/** Return whether or not the atom with the specified index is
    part of a chain
    
    \throw SireError::invalid_index
*/
bool MoleculeInfoData::isWithinChain(AtomIdx atomidx) const
{
    ResIdx residx = atoms_by_index[atomidx.map(atoms_by_index.count())].residx;

    if (not residx.isNull())
        return this->isWithinChain(residx);
        
    else
        return false;
}

/** Return whether or not the atom identified by the ID 'atomid' is
    part of a chain
    
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
bool MoleculeInfoData::isWithinChain(const AtomID &atomid) const
{
    ResIdx residx = atoms_by_index[ this->atomIdx(atomid) ].residx;
    
    if (not residx.isNull())
        return this->isWithinChain(residx);
        
    else    
        return false;
}

/** Return whether or not the atom with the specified index is
    part of a segment
    
    \throw SireError::invalid_index
*/
bool MoleculeInfoData::isWithinSegment(AtomIdx atomidx) const
{
    SegIdx segidx = atoms_by_index[atomidx.map(atoms_by_index.count())].segidx;
    return not segidx.isNull();
}

/** Return whether or not the atom identified by the ID 'atomid' is
    part of a segment
    
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
bool MoleculeInfoData::isWithinSegment(const AtomID &atomid) const
{
    SegIdx segidx = atoms_by_index[ this->atomIdx(atomid) ].segidx;
    return not segidx.isNull();
}

/** Return the index of the chain that contains the residue with index 'residx'

    \throw SireMol::missing_chain
    \throw SireError::invalid_index
*/
ChainIdx MoleculeInfoData::parentChain(ResIdx residx) const
{
    ChainIdx chainidx = res_by_index[residx.map(res_by_index.count())].chainidx;
    
    if (chainidx.isNull())
        throw SireMol::missing_chain( QObject::tr(
            "The residue at index 'residx' has not been added to any chain!")
                .arg(residx), CODELOC );
                
    return chainidx;
}

/** Return the index of the chain that contains the residue
    identified by the ID 'resid' 
    
    \throw SireMol::missing_residue
    \throw SireMol::missing_chain
    \throw SireError::invalid_index
*/
ChainIdx MoleculeInfoData::parentChain(const ResID &resid) const
{
    ChainIdx chainidx = res_by_index[ this->resIdx(resid) ].chainidx;
    
    if (chainidx.isNull())
        throw SireMol::missing_chain( QObject::tr(
            "The residue with ID \"%1\" has not been added to any chain!")
                .arg(resid.toString()), CODELOC );
                
    return chainidx;
}

/** Return the parent residue of the atom at index 'atomidx'

    \throw SireMol::missing_residue
    \throw SireError::invalid_index
*/
ResIdx MoleculeInfoData::parentResidue(AtomIdx atomidx) const
{
    ResIdx residx = atoms_by_index[atomidx.map(atoms_by_index.count())].residx;
    
    if (residx.isNull())
        throw SireMol::missing_residue( QObject::tr(
            "The atom at index %1 has not been added to any residue!")
                .arg(atomidx), CODELOC );
                
    return residx;
}

/** Return the parent residue of the atom identified by 'atomid'

    \throw SireMol::missing_residue
    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
ResIdx MoleculeInfoData::parentResidue(const AtomID &atomid) const
{
    ResIdx residx = atoms_by_index[ this->atomIdx(atomid) ].residx;
    
    if (residx.isNull())
        throw SireMol::missing_atom( QObject::tr(
            "The atom with ID \"%1\" has not been added to any residue!")
                .arg(atomid.toString()), CODELOC );
                
    return residx;
}

/** Return the parent chain of the atom at index 'atomidx'

    \throw SireError::invalid_index
*/
ChainIdx MoleculeInfoData::parentChain(AtomIdx atomidx) const
{
    return this->parentChain( this->parentResidue(atomidx) );
}

/** Return the parent chain of the atom identified by the ID 'atomid'

    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
ChainIdx MoleculeInfoData::parentChain(const AtomID &atomid) const
{
    return this->parentChain( this->parentResidue(atomid) );
}

/** Return the parent segment of the atom at index 'atomidx' 

    \throw SireMol::missing_segment
    \throw SireError::invalid_index 
*/
SegIdx MoleculeInfoData::parentSegment(AtomIdx atomidx) const
{
    SegIdx segidx = atoms_by_index[atomidx.map(atoms_by_index.count())].segidx;
    
    if (segidx.isNull())
        throw SireMol::missing_segment( QObject::tr(
            "The atom at index %1 has not been added to any segment!")
                .arg(atomidx), CODELOC );

    return segidx;
}

/** Return the parent segment of the atom identified by the ID 'atomid'

    \throw SireMol::missing_segment
    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
SegIdx MoleculeInfoData::parentSegment(const AtomID &atomid) const
{
    SegIdx segidx = atoms_by_index[ this->atomIdx(atomid) ].segidx;
    
    if (segidx.isNull())
        throw SireMol::missing_segment( QObject::tr(
            "The atom at ID \"%1\" has not been added to any segment!")
                .arg(atomid.toString()), CODELOC );
                
    return segidx;
}

/** Return the parent CutGroup of the atom at index 'atomid'. All
    atoms *have* to be placed into a CutGroup, so this will always
    return a valid result.

    \throw SireError::invalid_index
*/
CGIdx MoleculeInfoData::parentCutGroup(AtomIdx atomidx) const
{
    return atoms_by_index[atomidx.map(atoms_by_index.count())].cgatomidx.cutGroup();
}

/** Return the parent CutGroup of the atom identified by 'atomid'.
    All atoms *have* to be placed into a CutGroup, so this will always
    return a valid result.

    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
CGIdx MoleculeInfoData::parentCutGroup(const AtomID &atomid) const
{
    return atoms_by_index[this->atomIdx(atomid)].cgatomidx.cutGroup();
}

/** Return whether the residue at index 'residx' contains the atom 
    at index 'atomidx'
    
    \throw SireError::invalid_index
*/
bool MoleculeInfoData::contains(ResIdx residx, AtomIdx atomidx) const
{
    return res_by_index[residx.map(res_by_index.count())]
               .atom_indicies.contains( AtomIdx(atomidx.map(atoms_by_index.count())) );
}

/** Return whether the chain at index 'chainidx' contains the atom
    at index 'atomidx'
    
    \throw SireError::invalid_index
*/
bool MoleculeInfoData::contains(ChainIdx chainidx, AtomIdx atomidx) const
{
    foreach (ResIdx residx, chains_by_index[chainidx.map(chains_by_index.count())]
                                    .res_indicies)
    {
        if (this->contains(residx, atomidx))
            return true;
    }
    
    return false;
}

/** Return whether the segment at index 'segidx' contains the atom
    at index 'atomidx' 
    
    \throw SireError::invalid_index
*/
bool MoleculeInfoData::contains(SegIdx segidx, AtomIdx atomidx) const
{
    return seg_by_index[segidx.map(seg_by_index.count())]
              .atom_indicies.contains( AtomIdx(atomidx.map(atoms_by_index.count())) );
}

/** Return whether the CutGroup at index 'cgidx' contains the atom at 
    index 'atomidx'
    
    \throw SireError::invalid_index
*/
bool MoleculeInfoData::contains(CGIdx cgidx, AtomIdx atomidx) const
{
    return cg_by_index[cgidx.map(cg_by_index.count())]
            .atom_indicies.contains( AtomIdx(atomidx.map(atoms_by_index.count())) );
}

/** Return whether the chain at index 'chainidx' contains the residue
    at index 'residx'
    
    \throw SireError::invalid_index
*/
bool MoleculeInfoData::contains(ChainIdx chainidx, ResIdx residx) const
{
    return chains_by_index[chainidx.map(chains_by_index.count())]
            .res_indicies.contains( ResIdx(residx.map(res_by_index.count())) );
}

/** Return whether the residue at index 'residx' contains all of the 
    atoms that match the ID 'atomid'
    
    \throw SireError::invalid_index
*/
bool MoleculeInfoData::contains(ResIdx residx, const AtomID &atomid) const
{
    residx = ResIdx( residx.map(res_by_index.count()) );

    try
    {
        foreach(AtomIdx atomidx, atomid.map(*this))
        {
            if (not this->contains(residx, atomidx))
                return false;
        }
        
        return true;
    }
    catch(...)
    {
        return false;
    }
}

/** Return whether the chain at index 'chainidx' contains all of the 
    atoms identified by 'atomid'
    
    \throw SireError::invalid_index
*/
bool MoleculeInfoData::contains(ChainIdx chainidx, const AtomID &atomid) const
{
    chainidx = ChainIdx( chainidx.map(chains_by_index.count()) );
    
    try
    {
        foreach (AtomIdx atomidx, atomid.map(*this))
        {
            if (not this->contains(chainidx,atomidx))
                return false;
        }
        
        return true;
    }
    catch(...)
    {
        return false;
    }
}

/** Return whether the segment at index 'segid' contains all of the 
    atoms identified by 'atomid'
    
    \throw SireError::invalid_index
*/
bool MoleculeInfoData::contains(SegIdx segidx, const AtomID &atomid) const
{
    segidx = SegIdx( segidx.map(seg_by_index.count()) );
    
    try
    {
        foreach( AtomIdx atomidx, atomid.map(*this) )
        {
            if (not this->contains(segidx, atomidx))
                return false;
        }
        
        return true;
    }
    catch(...)
    {
        return false;
    }
}

/** Return whether the CutGroup at index 'cgidx' contains all of the
    atom identified by 'atomid' 
    
    \throw SireError::invalid_index
*/
bool MoleculeInfoData::contains(CGIdx cgidx, const AtomID &atomid) const
{
    cgidx = CGIdx( cgidx.map(cg_by_index.count()) );
    
    try
    {
        foreach(AtomIdx atomidx, atomid.map(*this))
        {
            if (not this->contains(cgidx, atomidx))
                return false;
        }
        
        return true;
    }
    catch(...)
    {
        return false;
    }
}

/** Return whether the chain at index 'chainidx' contains all of 
    the residues identified by 'resid'
    
    \throw SireError::invalid_index
*/
bool MoleculeInfoData::contains(ChainIdx chainidx, const ResID &resid) const
{
    chainidx = ChainIdx( chainidx.map(chains_by_index.count()) );
    
    try
    {
        foreach (ResIdx residx, resid.map(*this))
        {
            if (not this->contains(chainidx, residx))
                return false;
        }
        
        return true;
    }
    catch(...)
    {
        return false;
    }
}

/** Return whether or not the residue at index 'residx' contains any of
    the atoms identified by the ID 'atomid' */
bool MoleculeInfoData::intersects(ResIdx residx, const AtomID &atomid) const
{
    try
    {
        foreach (AtomIdx atomidx, atomid.map(*this))
        {
            if (this->contains(residx, atomidx))
                return true;
        }
        
        return false;
    }
    catch(...)
    {
        return false;
    }
}

/** Return whether or not the chain at index 'chainidx' contains any
    of the atoms identified by 'atomid' */
bool MoleculeInfoData::intersects(ChainIdx chainidx, const AtomID &atomid) const
{
    try
    {
        foreach (AtomIdx atomidx, atomid.map(*this))
        {
            if (this->contains(chainidx, atomidx))
                return true;
        }
        
        return false;
    }
    catch(...)
    {
        return false;
    }
}

/** Return whether or not the segment at index segidx contains any
    of the atoms identified by 'atomid' */
bool MoleculeInfoData::intersects(SegIdx segidx, const AtomID &atomid) const
{
    try
    {
        foreach (AtomIdx atomidx, atomid.map(*this))
        {
            if (this->contains(segidx, atomidx))
                return true;
        }
        
        return false;
    }
    catch(...)
    {
        return false;
    }
}

/** Return whether the CutGroup at index cgidx contains any of the 
    atoms identified by 'atomid' */
bool MoleculeInfoData::intersects(CGIdx cgidx, const AtomID &atomid) const
{
    try
    {
        foreach (AtomIdx atomidx, atomid.map(*this))
        {
            if (this->contains(cgidx, atomidx))
                return true;
        }
        
        return false;
    }
    catch(...)
    {
        return false;
    }
}

/** Return whether the chain at index chainidx contains any of residues
    identified by 'resid' */
bool MoleculeInfoData::intersects(ChainIdx chainidx, const ResID &resid) const
{
    try
    {
        foreach (ResIdx residx, resid.map(*this))
        {
            if (this->contains(chainidx, residx))
                return true;
        }
        
        return false;
    }
    catch(...)
    {
        return false;
    }
}

/** Return the number of atoms in the molecule */
int MoleculeInfoData::nAtoms() const
{
    return atoms_by_index.count();
}

int MoleculeInfoData::_pvt_nAtoms(const QVector<ResIdx> &residxs) const
{
    int nats = 0;
    
    foreach (ResIdx residx, residxs)
    {
        nats += res_by_index[residx].atom_indicies.count();
    }
    
    return nats;
}

int MoleculeInfoData::_pvt_nAtoms(const QList<ResIdx> &residxs) const
{
    int nats = 0;
    
    foreach (ResIdx residx, residxs)
    {
        nats += res_by_index[residx].atom_indicies.count();
    }
    
    return nats;
}

int MoleculeInfoData::_pvt_nAtoms(ChainIdx chainidx) const
{
    return this->_pvt_nAtoms( chains_by_index[chainidx].res_indicies );
}

/** Return the number of atoms in the chains identified by 'chainid'

    \throw SireMol::missing_chain
    \throw SireError::invalid_index
*/
int MoleculeInfoData::nAtoms(const ChainID &chainid) const
{
    QList<ChainIdx> chainidxs = chainid.map(*this);
    
    int nats = 0;
    
    foreach (ChainIdx chainidx, chainidxs)
    {
        nats += this->_pvt_nAtoms(chainidx);
    }
    
    return nats;
}

/** Return the number of atoms in the chain at index 'chainidx' 

    \throw SireError::invalid_index
*/
int MoleculeInfoData::nAtoms(ChainIdx chainidx) const
{
    if (chainidx.isNull())
        return this->nAtoms();

    chainidx = ChainIdx(chainidx.map(chains_by_index.count()));
    
    return this->_pvt_nAtoms(chainidx);
}

/** Return the number of atoms in the residue(s) identified by 'resid' 

    \throw SireMol::missing_residue
    \throw SireError::invalid_index
*/
int MoleculeInfoData::nAtoms(const ResID &resid) const
{
    return this->_pvt_nAtoms( resid.map(*this) );
}

/** Return the number of atoms in the residue at index 'residx'

    \throw SireError::invalid_index
*/
int MoleculeInfoData::nAtoms(ResIdx residx) const
{
    if (residx.isNull())
        return this->nAtoms();

    return res_by_index[residx.map(res_by_index.count())].atom_indicies.count();
}

/** Return the number of residues in the segment identified by 'segid'

    \throw SireMol::missing_segment
    \throw SireError::invalid_index
*/
int MoleculeInfoData::nAtoms(const SegID &segid) const
{
    QList<SegIdx> segidxs = segid.map(*this);
    
    int nats = 0;
    
    foreach (SegIdx segidx, segidxs)
    {
        nats += seg_by_index[segidx].atom_indicies.count();
    }
    
    return nats;
}

/** Return the number of atoms in the segment at index 'segidx' 

    \throw SireError::invalid_index
*/
int MoleculeInfoData::nAtoms(SegIdx segidx) const
{
    if (segidx.isNull())
        return this->nAtoms();
        
    return seg_by_index[segidx.map(seg_by_index.count())].atom_indicies.count();
}

/** Return the number of atoms in the CutGroup(s) that are 
    identified by 'cgid'
    
    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
*/
int MoleculeInfoData::nAtoms(const CGID &cgid) const
{
    QList<CGIdx> cgidxs = cgid.map(*this);
    
    int nats = 0;
    
    foreach (CGIdx cgidx, cgidxs)
    {
        nats += cg_by_index[cgidx].atom_indicies.count();
    }
    
    return nats;
}

/** Return the number of atoms in the CutGroup at index 'cgidx'

    \throw SireError::invalid_index
*/
int MoleculeInfoData::nAtoms(CGIdx cgidx) const
{
    if (cgidx.isNull())
        return this->nAtoms();
        
    return cg_by_index[cgidx.map(cg_by_index.count())].atom_indicies.count();
}

/** Return the number of residues in this molecule */
int MoleculeInfoData::nResidues() const
{
    return res_by_index.count();
}

/** Return the number of residues in the chain identified by 'chainid'

    \throw SireMol::missing_chain
    \throw SireError::invalid_index
*/
int MoleculeInfoData::nResidues(const ChainID &chainid) const
{
    QList<ChainIdx> chainidxs = chainid.map(*this);
    
    int nres = 0;
    
    foreach (ChainIdx chainidx, chainidxs)
    {
        nres += chains_by_index[chainidx].res_indicies.count();
    }
    
    return nres;
}

/** Return the number of residues in the chain at index 'chainidx'

    \throw SireError::invalid_index
*/
int MoleculeInfoData::nResidues(ChainIdx chainidx) const
{
    if (chainidx.isNull())
        return this->nResidues();
    
    return chains_by_index[chainidx.map(chains_by_index.count())].res_indicies.count();
}

/** Return the number of chains in this molecule */
int MoleculeInfoData::nChains() const
{
    return chains_by_index.count();
}

/** Return the number of CutGroups in this molecule */
int MoleculeInfoData::nCutGroups() const
{
    return cg_by_index.count();
}

/** Return the number of segments in this molecule */
int MoleculeInfoData::nSegments() const
{
    return seg_by_index.count();
}

/** Return the indicies of all of the residues in this molecule that
    are called 'name' - this returns the ResIdx indicies of the 
    residues sorted in the order that they appear in the molecule.
    This raises an exception if there are no residues with this
    name in this molecule.
    
    \throw SireMol::missing_residue
*/
QList<ResIdx> MoleculeInfoData::map(const ResName &name) const
{
    if (name.isNull())
        return this->getResidues();
    
    QList<ResIdx> residxs;
    
    if (name.isCaseSensitive())
    {
        residxs = res_by_name.values(name);
    }
    else
    {
        //search manually
        QString lower_name = QString(name).toLower();
        
        for (QMultiHash<QString,ResIdx>::const_iterator it = res_by_name.constBegin();
             it != res_by_name.constEnd();
             ++it)
        {
            if (it.key().toLower() == lower_name)
                residxs.append(it.value());
        }
    }
    
    if (residxs.isEmpty())
        throw SireMol::missing_residue( QObject::tr(
            "There is no residue called \"%1\" in the layout \"%2\". "
            "Available residues are %3.")
                .arg(name).arg(uid.toString())
                .arg(Sire::toString(res_by_name.keys())), CODELOC );
                
    qSort(residxs);
    return residxs;
}

/** Return the indicies of all of the residues in this molecule that
    have the number 'num' - this returns the ResIdx indicies of the 
    residues sorted in the order that they appear in the molecule.
    This raises an exception if there are no residues with this
    number in this molecule.
    
    \throw SireMol::missing_residue
*/
QList<ResIdx> MoleculeInfoData::map(ResNum num) const
{
    if (num.isNull())
        return this->getResidues();

    QList<ResIdx> residxs = res_by_num.values(num);
    
    if (residxs.isEmpty())
        throw SireMol::missing_residue( QObject::tr(
            "There is no residue with the number \"%1\" in the layout \"%2\".")
                .arg(num).arg(uid.toString()), CODELOC );
                
    qSort(residxs);
    return residxs;
}

/** Obvious overload function that allows the trivial mapping of  
    a ResIdx to a list of ResIdxs (that just contains the passed ResIdx.
    This does map the ResIdx into this molecule, and raises an exception
    if this in an invalid ResIdx
    
    \throw SireError::invalid_index
*/
QList<ResIdx> MoleculeInfoData::map(ResIdx idx) const
{
    if (idx.isNull())
        return this->getResidues();

    QList<ResIdx> residxs;
    residxs.append( ResIdx(idx.map(res_by_index.count())) );
    
    return residxs;
}

/** Return the indicies of residues that match the generic residue ID 'resid'.
    This returns the indicies sorted in the order the residues appear
    in the molecule, and raises an exception if there is no residue
    in the molecule that matches this ID.
    
    \throw SireMol::missing_residue
*/
QList<ResIdx> MoleculeInfoData::map(const ResID &resid) const
{
    if (resid.isNull())
        return this->getResidues();
        
    return resid.map(*this);
}
    
/** Return the indicies of chains that match the name 'name'. This
    returns the indicies sorted in the order the chains appear
    in the molecule, and raises an exception if there is no chain
    with this name.
    
    \throw SireMol::missing_chain
*/
QList<ChainIdx> MoleculeInfoData::map(const ChainName &name) const
{
    if (name.isNull())
        return this->getChains();

    QList<ChainIdx> chainidxs;

    if (name.isCaseSensitive())
    {
        chainidxs = chains_by_name.values(name);
    }
    else
    {
        //search manually
        QString lower_name = QString(name).toLower();
        
        for (QMultiHash<QString,ChainIdx>::const_iterator 
                                                it = chains_by_name.constBegin();
             it != chains_by_name.constEnd();
             ++it)
        {
            if (it.key().toLower() == lower_name)
                chainidxs.append(it.value());
        }
    }
    
    if (chainidxs.isEmpty())
        throw SireMol::missing_chain( QObject::tr(
            "There is no chain called \"%1\" in the layout \"%2\".")
                .arg(name).arg(uid.toString()), CODELOC );
                
    qSort(chainidxs);
    return chainidxs;
}

/** Obvious function that maps the ChainIdx to a list of ChainIdx
    objects
    
    \throw SireError::invalid_index
*/
QList<ChainIdx> MoleculeInfoData::map(ChainIdx idx) const
{
    if (idx.isNull())
        return this->getChains();

    QList<ChainIdx> chainidxs;
    chainidxs.append( ChainIdx(idx.map(chains_by_index.count())) );
    
    return chainidxs;
}

/** Return the indicies of chains that match the ChainID 'chainid'

    \throw SireMol::missing_chain
    \throw SireError::invalid_index
*/
QList<ChainIdx> MoleculeInfoData::map(const ChainID &chainid) const
{
    if (chainid.isNull())
        return this->getChains();

    return chainid.map(*this);
}

/** Return the indicies of the segments that have the name 'name'

    \throw SireMol::missing_segment
*/
QList<SegIdx> MoleculeInfoData::map(const SegName &name) const
{
    if (name.isNull())
        return this->getSegments();

    QList<SegIdx> segidxs;
    
    if (name.isCaseSensitive())
    {
        segidxs = seg_by_name.values(name);
    }
    else
    {
        //search manually
        QString lower_name = QString(name).toLower();
        
        for (QMultiHash<QString,SegIdx>::const_iterator it = seg_by_name.constBegin();
             it != seg_by_name.constEnd();
             ++it)
        {
            if (it.key().toLower() == lower_name)
                segidxs.append(it.value());
        }
    }
    
    if (segidxs.isEmpty())
        throw SireMol::missing_segment( QObject::tr(
            "There is no segment called \"%1\" in the layout \"%2\".")
                .arg(name).arg(uid.toString()), CODELOC );
                
    qSort(segidxs);
    return segidxs;
}

/** Obvious function that maps the SegIdx to a list of SegIdx
    objects
    
    \throw SireError::invalid_index
*/
QList<SegIdx> MoleculeInfoData::map(SegIdx idx) const
{
    if (idx.isNull())
        return this->getSegments();

    QList<SegIdx> segidxs;
    segidxs.append( SegIdx(idx.map(seg_by_index.count())) );
    
    return segidxs;
}

/** Return the indicies of segments that match the ID 'segid'

    \throw SireMol::missing_segment
    \throw SireError::invalid_index
*/
QList<SegIdx> MoleculeInfoData::map(const SegID &segid) const
{
    if (segid.isNull())
        return this->getSegments();

    return segid.map(*this);
}

/** Return the indicies of CutGroups that have the name 'name'

    \throw SireMol::missing_cutgroup
*/
QList<CGIdx> MoleculeInfoData::map(const CGName &name) const
{
    if (name.isNull())
        return this->getCutGroups();
    
    QList<CGIdx> cgidxs;
    
    if (name.isCaseSensitive())
    {
        cgidxs = cg_by_name.values(name);
    }
    else
    {
        //search manually
        QString lower_name = QString(name).toLower();
        
        for (QMultiHash<QString,CGIdx>::const_iterator it = cg_by_name.constBegin();
             it != cg_by_name.constEnd();
             ++it)
        {
            if (it.key().toLower() == lower_name)
                cgidxs.append(it.value());
        }
    }
    
    if (cgidxs.isEmpty())
        throw SireMol::missing_cutgroup( QObject::tr(
            "There is no CutGroup called \"%1\" in the layout \"%2\".")
                .arg(name).arg(uid.toString()), CODELOC );
                
    qSort(cgidxs);
    return cgidxs;
}

/** Obvious function that maps the CGIdx to a list of CGIdx
    objects
    
    \throw SireError::invalid_index
*/
QList<CGIdx> MoleculeInfoData::map(CGIdx idx) const
{
    if (idx.isNull())
        return this->getCutGroups();

    QList<CGIdx> cgidxs;
    cgidxs.append( CGIdx(idx.map(cg_by_index.count())) );
    
    return cgidxs;
}

/** Return the indicies of CutGroups that match the ID 'cgid'

    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
*/
QList<CGIdx> MoleculeInfoData::map(const CGID &cgid) const
{
    if (cgid.isNull())
        return this->getCutGroups();

    return cgid.map(*this);
}
    
/** Return the indicies of atoms that are called 'name'
 
    \throw SireMol::missing_atom
*/
QList<AtomIdx> MoleculeInfoData::map(const AtomName &name) const
{
    if (name.isNull())
        return this->getAtoms();

    QList<AtomIdx> atomidxs;

    if (name.isCaseSensitive())
    {
        atomidxs = atoms_by_name.values(name);
    }
    else
    {
        //search manually...
        QString lower_name = QString(name).toLower();
                
        for (QMultiHash<QString,AtomIdx>::const_iterator it = atoms_by_name.constBegin();
             it != atoms_by_name.constEnd();
             ++it)
        {
            if (it.key().toLower() == lower_name)
                atomidxs.append( it.value() );
        }
    }
    
    if (atomidxs.isEmpty())
        throw SireMol::missing_atom( QObject::tr(
            "There is no atom called \"%1\" in the layout \"%2\".")
                .arg(name).arg(uid.toString()), CODELOC );
                
    qSort(atomidxs);
    return atomidxs;
}

/** Return the indicies of atoms that have the number 'num'.

    \throw SireMol::missing_atom
*/
QList<AtomIdx> MoleculeInfoData::map(AtomNum num) const
{
    if (num.isNull())
        return this->getAtoms();

    QList<AtomIdx> atomidxs = atoms_by_num.values(num);
    
    if (atomidxs.isEmpty())
        throw SireMol::missing_atom( QObject::tr(
            "There is no atom with the number \"%1\" in the layout \"%2\".")
                .arg(num).arg(uid.toString()), CODELOC );
                
    qSort(atomidxs);
    return atomidxs;
}

/** Obvious function that maps the AtomIdx to a list of AtomIdx
    objects
    
    \throw SireError::invalid_index
*/
QList<AtomIdx> MoleculeInfoData::map(AtomIdx idx) const
{
    if (idx.isNull())
        return this->getAtoms();

    QList<AtomIdx> atomidxs;
    atomidxs.append( AtomIdx(idx.map(atoms_by_index.count())) );
    
    return atomidxs;
}

/** Return the indicies of atoms that match the ID 'atomid'

    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
QList<AtomIdx> MoleculeInfoData::map(const AtomID &atomid) const
{
    if (atomid.isNull())
        return this->getAtoms();

    return atomid.map(*this);
}

/** Assert that this molecule contains the atom at index 'atomidx' 

    \throw SireError::invalid_index
*/
void MoleculeInfoData::assertContains(AtomIdx atomidx) const
{
    try
    {
        atomidx.map(this->nAtoms());
    }
    catch(const SireError::invalid_index&)
    {
        throw SireError::invalid_index( QObject::tr(
            "There is no atom at index %1. nAtoms() == %2 for the "
            "molecule layout with UID %3.")
                .arg(atomidx).arg(nAtoms()).arg(UID().toString()), CODELOC );
    }
}

/** Assert that this molecule contains the CutGroup at index 'cgidx' 

    \throw SireError::invalid_index
*/
void MoleculeInfoData::assertContains(CGIdx cgidx) const
{
    try
    {
        cgidx.map(this->nCutGroups());
    }
    catch(const SireError::invalid_index&)
    {
        throw SireError::invalid_index( QObject::tr(
            "There is no CutGroup at index %1. nCutGroups() == %2 for the "
            "molecule layout with UID %3.")
                .arg(cgidx).arg(nCutGroups()).arg(UID().toString()), CODELOC );
    }
}

/** Assert that this molecule contains the residue at index 'residx' 

    \throw SireError::invalid_index
*/
void MoleculeInfoData::assertContains(ResIdx residx) const
{
    try
    {
        residx.map(this->nResidues());
    }
    catch(const SireError::invalid_index&)
    {
        throw SireError::invalid_index( QObject::tr(
            "There is no residue at index %1. nResidues() == %2 for the "
            "molecule layout with UID %3.")
                .arg(residx).arg(nResidues()).arg(UID().toString()), CODELOC );
    }
}

/** Assert that this molecule contains the chain at index 'chainidx' 

    \throw SireError::invalid_index
*/
void MoleculeInfoData::assertContains(ChainIdx chainidx) const
{
    try
    {
        chainidx.map(this->nChains());
    }
    catch(const SireError::invalid_index&)
    {
        throw SireError::invalid_index( QObject::tr(
            "There is no chain at index %1. nChains() == %2 for the "
            "molecule layout with UID %3.")
                .arg(chainidx).arg(nChains()).arg(UID().toString()), CODELOC );
    }
}

/** Assert that this molecule contains the atom at index 'segidx' 

    \throw SireError::invalid_index
*/
void MoleculeInfoData::assertContains(SegIdx segidx) const
{
    try
    {
        segidx.map(this->nSegments());
    }
    catch(const SireError::invalid_index&)
    {
        throw SireError::invalid_index( QObject::tr(
            "There is no atom at index %1. nSegments() == %2 for the "
            "molecule layout with UID %3.")
                .arg(segidx).arg(nSegments()).arg(UID().toString()), CODELOC );
    }
}

void MoleculeInfoData::assertCompatibleWith(
                             const AtomSelection &selected_atoms) const
{
    return selected_atoms.assertCompatibleWith(*this);
}

const char* MoleculeInfoData::typeName()
{
    return QMetaType::typeName( qMetaTypeId<MoleculeInfoData>() );
}
