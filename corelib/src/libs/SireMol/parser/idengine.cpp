/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2018  Christopher Woods
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

#include "idengine.h"

#include "SireBase/parallel.h"
#include "SireBase/booleanproperty.h"

#include "tostring.h"

using namespace SireMol;
using namespace SireBase;
using namespace parser_idengine;

////////
//////// Implementation of the IDNameEngine
////////

IDNameEngine::IDNameEngine() : SelectEngine()
{}

SelectEnginePtr IDNameEngine::construct(IDObject o, NameValues vals)
{
    IDNameEngine *ptr = new IDNameEngine();
    ptr->obj = o;

    try
    {
        for (const auto val : vals)
        {
            if (val.value.which() == 0)
            {
                RegExpValue v = boost::get<RegExpValue>(val.value);
                QString r = QString::fromStdString(v.value);
                
                QRegularExpression regexp;
                
                if (v.is_case_sensitive)
                    regexp = QRegularExpression(r);
                else
                    regexp = QRegularExpression(r, QRegularExpression::CaseInsensitiveOption);
                
                if (not regexp.isValid())
                {
                    throw SireMol::parse_error( QObject::tr("Failed to interpret the "
                      "regular expression '%1' (escaped version is '%2'). Error is '%3'")
                        .arg( QString::fromStdString(v.value) )
                        .arg(r)
                        .arg(regexp.errorString()), CODELOC );
                }
                
                //optimise (JIT-compile) the regular expression now as it will
                //be used a lot
                #if QT_VERSION >= 0x504000
                    regexp.optimize();  // introduced in Qt 5.4
                #endif
                
                ptr->regexps.append(regexp);
            }
            else if (val.value.which() == 1)
                ptr->names.append( QString::fromStdString(boost::get<std::string>(val.value)) );
        }
    }
    catch(...)
    {
        delete ptr;
        throw;
    }
    
    return makePtr(ptr);
}

IDNameEngine::~IDNameEngine()
{}

/** Function used to find the matching parts in mols using 'searchFunction' */
template<class T>
SelectResult searchForMatch(const SelectResult &mols, const T &searchFunction, bool uses_parallel)
{
    //now loop through all of the molecules and find the matching atoms
    SelectResult::Container result;
    
    if (uses_parallel)
    {
        QVector<ViewsOfMol> tmpresult(mols.count());
        const auto molviews = mols.views();

        tbb::parallel_for( tbb::blocked_range<int>(0,mols.count()),
                           [&](const tbb::blocked_range<int> &r)
        {
            for (int i=r.begin(); i<r.end(); ++i)
            {
                tmpresult[i] = searchFunction(molviews[i]);
            }
        });
        
        for (const auto &r : tmpresult)
        {
            if (not r.isEmpty())
                result.append(r);
        }
    }
    else
    {
        for (const auto mol : mols)
        {
            auto match = searchFunction(mol);
            
            if (not match.isEmpty())
                result.append(match);
        }
    }

    return result;
}

/** Function used to find all of the atoms that match by name */
SelectResult IDNameEngine::selectAtoms(const SelectResult &mols, bool uses_parallel) const
{
    // function that finds all of the atoms that have been selected from the molecule
    auto selectFromMol = [&](const ViewsOfMol &mol)
    {
        const auto molinfo = mol.data().info();

        // function that tests whether or not the passed atom has been selected
        auto selectFromAtom = [&](const AtomIdx &idx)
        {
            const auto atomname = molinfo.name(idx).value();

            //try all of the fixed names
            for (const auto name : names)
            {
                if (name == atomname)
                {
                    //name matches exactly
                    return true;
                }
            }
            
            //now try all of the regexps
            for (const auto regexp : regexps)
            {
                auto match = regexp.match(atomname);
                
                if (match.hasMatch())
                {
                    //we have a regexp match :-)
                    return true;
                }
            }
            
            return false;
        };

        QList<AtomIdx> selected_atoms;

        if (mol.selectedAll())
        {
            const int natoms = molinfo.nAtoms();
        
            if (uses_parallel)
            {
                QMutex mutex;
            
                tbb::parallel_for( tbb::blocked_range<int>(0,natoms),
                                   [&](const tbb::blocked_range<int> &r)
                {
                    QList<AtomIdx> atoms;
                    
                    for (int i=r.begin(); i<r.end(); ++i)
                    {
                        const AtomIdx idx(i);
                        
                        if (selectFromAtom(idx))
                            atoms.append(idx);
                    }
                    
                    QMutexLocker lkr(&mutex);
                    selected_atoms += atoms;
                });
            }
            else
            {
                for (int i=0; i<natoms; ++i)
                {
                    const AtomIdx idx(i);
                
                    if (selectFromAtom(idx))
                        selected_atoms.append(idx);
                }
            }
        }
        else
        {
            const auto view_atoms = mol.selection().selectedAtoms();
            
            if (uses_parallel)
            {
                QMutex mutex;
            
                tbb::parallel_for( tbb::blocked_range<int>(0,view_atoms.count()),
                                   [&](const tbb::blocked_range<int> &r)
                {
                    QList<AtomIdx> atoms;
                    
                    for (int i=r.begin(); i<r.end(); ++i)
                    {
                        const auto atom = view_atoms.constData()[i];
                    
                        if (selectFromAtom(atom))
                            atoms.append(atom);
                    }
                    
                    QMutexLocker lkr(&mutex);
                    selected_atoms += atoms;
                });
            }
            else
            {
                for (const auto atom : view_atoms)
                {
                    if (selectFromAtom(atom))
                        selected_atoms.append(atom);
                }
            }
        }
        
        if (selected_atoms.isEmpty())
            //no atoms matched
            return ViewsOfMol();
        else if (selected_atoms.count() == 1)
            //only a single atom matched
            return ViewsOfMol( Atom(mol.data(),selected_atoms[0]) );
        else if (selected_atoms.count() == molinfo.nAtoms())
            //the entire molecule matched
            return ViewsOfMol( mol.molecule() );
        else
        {
            //a subset of the molecule matches
            return ViewsOfMol( mol.data(),
                               Selector<Atom>(mol.data(),selected_atoms).selection() );
        }
    };
    
    return searchForMatch(mols, selectFromMol, uses_parallel);
}

SelectResult IDNameEngine::selectCutGroups(const SelectResult &mols, bool uses_parallel) const
{
    // function that finds all of the cutgroups that have been selected from the molecule
    auto selectFromMol = [&](const ViewsOfMol &mol)
    {
        const auto molinfo = mol.data().info();

        // function that tests whether or not the passed cutgroup has been selected
        auto selectFromCutGroup = [&](const CGIdx &idx)
        {
            const auto cgname = molinfo.name(idx).value();

            //try all of the fixed names
            for (const auto name : names)
            {
                if (name == cgname)
                {
                    //name matches exactly
                    return true;
                }
            }
            
            //now try all of the regexps
            for (const auto regexp : regexps)
            {
                auto match = regexp.match(cgname);
                
                if (match.hasMatch())
                {
                    //we have a regexp match :-)
                    return true;
                }
            }
            
            return false;
        };

        QList<CGIdx> selected_cutgroups;

        if (mol.selectedAll())
        {
            const int ncgs = molinfo.nCutGroups();
        
            if (uses_parallel)
            {
                QMutex mutex;
            
                tbb::parallel_for( tbb::blocked_range<int>(0,ncgs),
                                   [&](const tbb::blocked_range<int> &r)
                {
                    QList<CGIdx> cutgroups;
                    
                    for (int i=r.begin(); i<r.end(); ++i)
                    {
                        const CGIdx idx(i);
                        
                        if (selectFromCutGroup(idx))
                            cutgroups.append(idx);
                    }
                    
                    QMutexLocker lkr(&mutex);
                    selected_cutgroups += cutgroups;
                });
            }
            else
            {
                for (int i=0; i<ncgs; ++i)
                {
                    const CGIdx idx(i);
                
                    if (selectFromCutGroup(idx))
                        selected_cutgroups.append(idx);
                }
            }
        }
        else
        {
            const auto view_cutgroups = mol.selection().selectedCutGroups();
            
            if (uses_parallel)
            {
                QMutex mutex;
            
                tbb::parallel_for( tbb::blocked_range<int>(0,view_cutgroups.count()),
                                   [&](const tbb::blocked_range<int> &r)
                {
                    QList<CGIdx> cutgroups;
                    
                    for (int i=r.begin(); i<r.end(); ++i)
                    {
                        const auto cutgroup = view_cutgroups[i];
                    
                        if (selectFromCutGroup(cutgroup))
                            cutgroups.append(cutgroup);
                    }
                    
                    QMutexLocker lkr(&mutex);
                    selected_cutgroups += cutgroups;
                });
            }
            else
            {
                for (const auto cutgroup : view_cutgroups)
                {
                    if (selectFromCutGroup(cutgroup))
                        selected_cutgroups.append(cutgroup);
                }
            }
        }
        
        if (selected_cutgroups.isEmpty())
            //no cutgroups matched
            return ViewsOfMol();
        else if (selected_cutgroups.count() == 1)
            //only a single cutgroup matched
            return ViewsOfMol( CutGroup(mol.data(),selected_cutgroups[0]) );
        else if (selected_cutgroups.count() == molinfo.nCutGroups())
            //the entire molecule matched
            return ViewsOfMol( mol.molecule() );
        else
        {
            //a subset of the molecule matches
            return ViewsOfMol( mol.data(),
                               Selector<CutGroup>(mol.data(),selected_cutgroups).selection() );
        }
    };
    
    return searchForMatch(mols, selectFromMol, uses_parallel);
}

SelectResult IDNameEngine::selectResidues(const SelectResult &mols, bool uses_parallel) const
{
    // function that finds all of the atoms that have been selected from the molecule
    auto selectFromMol = [&](const ViewsOfMol &mol)
    {
        const auto molinfo = mol.data().info();

        // function that tests whether or not the passed residue has been selected
        auto selectFromResidue = [&](const ResIdx &idx)
        {
            const auto resname = molinfo.name(idx).value();

            //try all of the fixed names
            for (const auto name : names)
            {
                if (name == resname)
                {
                    //name matches exactly
                    return true;
                }
            }
            
            //now try all of the regexps
            for (const auto regexp : regexps)
            {
                auto match = regexp.match(resname);
                
                if (match.hasMatch())
                {
                    //we have a regexp match :-)
                    return true;
                }
            }
            
            return false;
        };

        QList<ResIdx> selected_residues;

        if (mol.selectedAll())
        {
            const int nres = molinfo.nResidues();
        
            if (uses_parallel)
            {
                QMutex mutex;
            
                tbb::parallel_for( tbb::blocked_range<int>(0,nres),
                                   [&](const tbb::blocked_range<int> &r)
                {
                    QList<ResIdx> residues;
                    
                    for (int i=r.begin(); i<r.end(); ++i)
                    {
                        const ResIdx idx(i);
                        
                        if (selectFromResidue(idx))
                            residues.append(idx);
                    }
                    
                    QMutexLocker lkr(&mutex);
                    selected_residues += residues;
                });
            }
            else
            {
                for (int i=0; i<nres; ++i)
                {
                    const ResIdx idx(i);
                
                    if (selectFromResidue(idx))
                        selected_residues.append(idx);
                }
            }
        }
        else
        {
            const auto view_residues = mol.selection().selectedResidues();
            
            if (uses_parallel)
            {
                QMutex mutex;
            
                tbb::parallel_for( tbb::blocked_range<int>(0,view_residues.count()),
                                   [&](const tbb::blocked_range<int> &r)
                {
                    QList<ResIdx> residues;
                    
                    for (int i=r.begin(); i<r.end(); ++i)
                    {
                        const auto residue = view_residues[i];
                    
                        if (selectFromResidue(residue))
                            residues.append(residue);
                    }
                    
                    QMutexLocker lkr(&mutex);
                    selected_residues += residues;
                });
            }
            else
            {
                for (const auto residue : view_residues)
                {
                    if (selectFromResidue(residue))
                        selected_residues.append(residue);
                }
            }
        }
        
        if (selected_residues.isEmpty())
            //no residues matched
            return ViewsOfMol();
        else if (selected_residues.count() == 1)
            //only a single residue matched
            return ViewsOfMol( Residue(mol.data(),selected_residues[0]) );
        else if (selected_residues.count() == molinfo.nResidues())
            //the entire molecule matched
            return ViewsOfMol( mol.molecule() );
        else
        {
            //a subset of the molecule matches
            return ViewsOfMol( mol.data(),
                               Selector<Residue>(mol.data(),selected_residues).selection() );
        }
    };
    
    return searchForMatch(mols, selectFromMol, uses_parallel);
}

SelectResult IDNameEngine::selectChains(const SelectResult &mols, bool uses_parallel) const
{
    // function that finds all of the atoms that have been selected from the molecule
    auto selectFromMol = [&](const ViewsOfMol &mol)
    {
        const auto molinfo = mol.data().info();

        // function that tests whether or not the passed chain has been selected
        auto selectFromChain = [&](const ChainIdx &idx)
        {
            const auto chainname = molinfo.name(idx).value();

            //try all of the fixed names
            for (const auto name : names)
            {
                if (name == chainname)
                {
                    //name matches exactly
                    return true;
                }
            }
            
            //now try all of the regexps
            for (const auto regexp : regexps)
            {
                auto match = regexp.match(chainname);
                
                if (match.hasMatch())
                {
                    //we have a regexp match :-)
                    return true;
                }
            }
            
            return false;
        };

        QList<ChainIdx> selected_chains;

        if (mol.selectedAll())
        {
            const int nchains = molinfo.nChains();
        
            if (uses_parallel)
            {
                QMutex mutex;
            
                tbb::parallel_for( tbb::blocked_range<int>(0,nchains),
                                   [&](const tbb::blocked_range<int> &r)
                {
                    QList<ChainIdx> chains;
                    
                    for (int i=r.begin(); i<r.end(); ++i)
                    {
                        const ChainIdx idx(i);
                        
                        if (selectFromChain(idx))
                            chains.append(idx);
                    }
                    
                    QMutexLocker lkr(&mutex);
                    selected_chains += chains;
                });
            }
            else
            {
                for (int i=0; i<nchains; ++i)
                {
                    const ChainIdx idx(i);
                
                    if (selectFromChain(idx))
                        selected_chains.append(idx);
                }
            }
        }
        else
        {
            const auto view_chains = mol.selection().selectedChains();
            
            if (uses_parallel)
            {
                QMutex mutex;
            
                tbb::parallel_for( tbb::blocked_range<int>(0,view_chains.count()),
                                   [&](const tbb::blocked_range<int> &r)
                {
                    QList<ChainIdx> chains;
                    
                    for (int i=r.begin(); i<r.end(); ++i)
                    {
                        const auto chain = view_chains[i];
                    
                        if (selectFromChain(chain))
                            chains.append(chain);
                    }
                    
                    QMutexLocker lkr(&mutex);
                    selected_chains += chains;
                });
            }
            else
            {
                for (const auto chain : view_chains)
                {
                    if (selectFromChain(chain))
                        selected_chains.append(chain);
                }
            }
        }
        
        if (selected_chains.isEmpty())
            //no chains matched
            return ViewsOfMol();
        else if (selected_chains.count() == 1)
            //only a single chain matched
            return ViewsOfMol( Chain(mol.data(),selected_chains[0]) );
        else if (selected_chains.count() == molinfo.nChains())
            //the entire molecule matched
            return ViewsOfMol( mol.molecule() );
        else
        {
            //a subset of the molecule matches
            return ViewsOfMol( mol.data(),
                               Selector<Chain>(mol.data(),selected_chains).selection() );
        }
    };

    return searchForMatch(mols, selectFromMol, uses_parallel);
}

SelectResult IDNameEngine::selectSegments(const SelectResult &mols, bool uses_parallel) const
{
    // function that finds all of the segments that have been selected from the molecule
    auto selectFromMol = [&](const ViewsOfMol &mol)
    {
        const auto molinfo = mol.data().info();

        // function that tests whether or not the passed segment has been selected
        auto selectFromSegment = [&](const SegIdx &idx)
        {
            const auto segname = molinfo.name(idx).value();

            //try all of the fixed names
            for (const auto name : names)
            {
                if (name == segname)
                {
                    //name matches exactly
                    return true;
                }
            }
            
            //now try all of the regexps
            for (const auto regexp : regexps)
            {
                auto match = regexp.match(segname);
                
                if (match.hasMatch())
                {
                    //we have a regexp match :-)
                    return true;
                }
            }
            
            return false;
        };

        QList<SegIdx> selected_segments;

        if (mol.selectedAll())
        {
            const int nsegs = molinfo.nSegments();
        
            if (uses_parallel)
            {
                QMutex mutex;
            
                tbb::parallel_for( tbb::blocked_range<int>(0,nsegs),
                                   [&](const tbb::blocked_range<int> &r)
                {
                    QList<SegIdx> segments;
                    
                    for (int i=r.begin(); i<r.end(); ++i)
                    {
                        const SegIdx idx(i);
                        
                        if (selectFromSegment(idx))
                            segments.append(idx);
                    }
                    
                    QMutexLocker lkr(&mutex);
                    selected_segments += segments;
                });
            }
            else
            {
                for (int i=0; i<nsegs; ++i)
                {
                    const SegIdx idx(i);
                
                    if (selectFromSegment(idx))
                        selected_segments.append(idx);
                }
            }
        }
        else
        {
            const auto view_segments = mol.selection().selectedSegments();
            
            if (uses_parallel)
            {
                QMutex mutex;
            
                tbb::parallel_for( tbb::blocked_range<int>(0,view_segments.count()),
                                   [&](const tbb::blocked_range<int> &r)
                {
                    QList<SegIdx> segments;
                    
                    for (int i=r.begin(); i<r.end(); ++i)
                    {
                        const auto segment = view_segments[i];
                    
                        if (selectFromSegment(segment))
                            segments.append(segment);
                    }
                    
                    QMutexLocker lkr(&mutex);
                    selected_segments += segments;
                });
            }
            else
            {
                for (const auto segment : view_segments)
                {
                    if (selectFromSegment(segment))
                        selected_segments.append(segment);
                }
            }
        }
        
        if (selected_segments.isEmpty())
            //no segments matched
            return ViewsOfMol();
        else if (selected_segments.count() == 1)
            //only a single segments matched
            return ViewsOfMol( Segment(mol.data(),selected_segments[0]) );
        else if (selected_segments.count() == molinfo.nSegments())
            //the entire molecule matched
            return ViewsOfMol( mol.molecule() );
        else
        {
            //a subset of the molecule matches
            return ViewsOfMol( mol.data(),
                               Selector<Segment>(mol.data(),selected_segments).selection() );
        }
    };

    return searchForMatch(mols, selectFromMol, uses_parallel);
}

SelectResult IDNameEngine::selectMolecules(const SelectResult &mols, bool uses_parallel) const
{
    // function that finds all of the molecules that have been selected
    auto selectFromMol = [&](const ViewsOfMol &mol)
    {
        const auto molname = mol.data().name().value();
        
        //try all of the fixed names
        for (const auto name : names)
        {
            if (name == molname)
            {
                //name matches exactly
                return ViewsOfMol(mol.molecule());
            }
        }
        
        //now try all of the regexps
        for (const auto regexp : regexps)
        {
            auto match = regexp.match(molname);
            
            if (match.hasMatch())
            {
                //we have a regexp match :-)
                return ViewsOfMol(mol.molecule());
            }
        }
        
        //no match
        return ViewsOfMol();
    };

    return searchForMatch(mols, selectFromMol, uses_parallel);
}

SelectResult IDNameEngine::select(const SelectResult &mols, const PropertyMap &map) const
{
    SelectResult::Container result;

    bool uses_parallel = true;
    
    if (map["parallel"].hasValue())
    {
        uses_parallel = map["parallel"].value().asA<BooleanProperty>().value();
    }
    
    switch(obj)
    {
    case AST::ATOM:
        return selectAtoms(mols, uses_parallel);
    case AST::CUTGROUP:
        return selectCutGroups(mols, uses_parallel);
    case AST::RESIDUE:
        return selectResidues(mols, uses_parallel);
    case AST::CHAIN:
        return selectChains(mols, uses_parallel);
    case AST::SEGMENT:
        return selectSegments(mols, uses_parallel);
    case AST::MOLECULE:
        return selectMolecules(mols, uses_parallel);
    default:
        return SelectResult();
    }
}

SelectEngine::ObjType IDNameEngine::objectType() const
{
    switch(obj)
    {
    case AST::ATOM:
        return SelectEngine::ATOM;
    case AST::CUTGROUP:
        return SelectEngine::CUTGROUP;
    case AST::RESIDUE:
        return SelectEngine::RESIDUE;
    case AST::CHAIN:
        return SelectEngine::CHAIN;
    case AST::SEGMENT:
        return SelectEngine::SEGMENT;
    case AST::MOLECULE:
        return SelectEngine::MOLECULE;
    default:
        return SelectEngine::COMPLEX;
    }
}

////////
//////// Implementation of the IDNumberEngine
////////

/** Internal function used to see if the passed integer matches any of the 
    range values */
static bool match(const int idx, const RangeValues &vals)
{
    for (const auto val : vals)
    {
        if (val.which() == 0)
        {
            auto v = boost::get<RangeValue>(val);
            
            if (v.start <= v.end)
            {
                for (int i=v.start; i<=v.end; i+=v.step)
                {
                    if (i == idx)
                        return true;
                }
            }
            else
            {
                for (int i=v.start; i>=v.end; i-=v.step)
                {
                    if (i == idx)
                        return true;
                }
            }
        }
        else
        {
            auto v = boost::get<CompareValue>(val);
            
            switch(v.compare)
            {
            case ID_CMP_LT:
                if (idx < v.value)
                    return true;
                break;
            case ID_CMP_LE:
                if (idx <= v.value)
                    return true;
                break;
            case ID_CMP_EQ:
                if (idx == v.value)
                    return true;
                break;
            case ID_CMP_NE:
                if (idx != v.value)
                    return true;
                break;
            case ID_CMP_GE:
                if (idx >= v.value)
                    return true;
                break;
            case ID_CMP_GT:
                if (idx > v.value)
                    return true;
                break;
            default:
                return false;
            }
        }
    }

    return false;
}

IDNumberEngine::IDNumberEngine()
{}

SelectEnginePtr IDNumberEngine::construct( IDObject o, RangeValues v )
{
    IDNumberEngine *ptr = new IDNumberEngine();
    ptr->obj = o;
    ptr->vals = v;

    return makePtr(ptr);
}

IDNumberEngine::~IDNumberEngine()
{}

template<class IDX, class OBJ, class T, class U, class V>
ViewsOfMol genericSelect(const ViewsOfMol &mol, const T &selectFunc,
                         const U &countFunc, const V &getSelectedFunc,
                         bool uses_parallel)
{
    const auto molinfo = mol.data().info();

    QList<IDX> selected;

    if (mol.selectedAll())
    {
        const int count = countFunc(molinfo);
    
        if (uses_parallel)
        {
            QMutex mutex;
        
            tbb::parallel_for( tbb::blocked_range<int>(0,count),
                               [&](const tbb::blocked_range<int> &r)
            {
                QList<IDX> thread_selected;
                
                for (int i=r.begin(); i<r.end(); ++i)
                {
                    const IDX idx(i);
                    
                    if (selectFunc(molinfo, idx))
                        thread_selected.append(idx);
                }
                
                QMutexLocker lkr(&mutex);
                selected += thread_selected;
            });
        }
        else
        {
            for (int i=0; i<count; ++i)
            {
                const IDX idx(i);
            
                if (selectFunc(molinfo, idx))
                    selected.append(idx);
            }
        }
    }
    else
    {
        const auto viewed = getSelectedFunc(mol);
        
        if (uses_parallel)
        {
            QMutex mutex;
        
            tbb::parallel_for( tbb::blocked_range<int>(0,viewed.count()),
                               [&](const tbb::blocked_range<int> &r)
            {
                QList<IDX> thread_selected;
                
                for (int i=r.begin(); i<r.end(); ++i)
                {
                    const auto idx = viewed[i];
                
                    if (selectFunc(molinfo, idx))
                        thread_selected.append(idx);
                }
                
                QMutexLocker lkr(&mutex);
                selected += thread_selected;
            });
        }
        else
        {
            for (const auto idx : viewed)
            {
                if (selectFunc(molinfo, idx))
                    selected.append(idx);
            }
        }
    }
    
    if (selected.isEmpty())
        //nothing matched
        return ViewsOfMol();
    else if (selected.count() == 1)
        //only a single object matched
        return ViewsOfMol( OBJ(mol.data(),selected[0]) );
    else if (selected.count() == countFunc(molinfo))
        //the entire molecule matched
        return ViewsOfMol( mol.molecule() );
    else
    {
        //a subset of the molecule matches
        return ViewsOfMol( mol.data(),
                           Selector<OBJ>(mol.data(),selected).selection() );
    }
}

static auto countAtoms = [](const MoleculeInfoData &molinfo)
{
    return molinfo.nAtoms();
};

static auto getSelectedAtoms = [](const ViewsOfMol &mol)
{
    return mol.selection().selectedAtoms();
};

static auto countCutGroups = [](const MoleculeInfoData &molinfo)
{
    return molinfo.nCutGroups();
};

static auto getSelectedCutGroups = [](const ViewsOfMol &mol)
{
    return mol.selection().selectedCutGroups();
};

static auto countResidues = [](const MoleculeInfoData &molinfo)
{
    return molinfo.nResidues();
};

static auto getSelectedResidues = [](const ViewsOfMol &mol)
{
    return mol.selection().selectedResidues();
};

static auto countChains = [](const MoleculeInfoData &molinfo)
{
    return molinfo.nChains();
};

static auto getSelectedChains = [](const ViewsOfMol &mol)
{
    return mol.selection().selectedChains();
};

static auto countSegments = [](const MoleculeInfoData &molinfo)
{
    return molinfo.nSegments();
};

static auto getSelectedSegments = [](const ViewsOfMol &mol)
{
    return mol.selection().selectedSegments();
};

/** Function used to find all of the atoms that match by number */
SelectResult IDNumberEngine::selectAtoms(const SelectResult &mols, bool uses_parallel) const
{
    auto matchAtom = [&](const MoleculeInfoData &molinfo, const AtomIdx &idx)
    {
        const auto number = molinfo.number(idx).value();
        return match(number, vals);
    };

    auto selectFromMol = [&](const ViewsOfMol &mol)
    {
        return ::genericSelect<AtomIdx,Atom>(mol, matchAtom, countAtoms,
                                             getSelectedAtoms, uses_parallel);
    };
    
    return searchForMatch(mols, selectFromMol, uses_parallel);
}

/** Function used to find all of the residues that match by number */
SelectResult IDNumberEngine::selectResidues(const SelectResult &mols, bool uses_parallel) const
{
    auto matchResidue = [&](const MoleculeInfoData &molinfo, const ResIdx &idx)
    {
        const auto number = molinfo.number(idx).value();
        qDebug() << "COMPARE" << number;
        return match(number, vals);
    };

    auto selectFromMol = [&](const ViewsOfMol &mol)
    {
        return ::genericSelect<ResIdx,Residue>(mol, matchResidue, countResidues,
                                               getSelectedResidues, uses_parallel);
    };
    
    return searchForMatch(mols, selectFromMol, uses_parallel);
}

/** Function used to find all of the molecules that match by number */
SelectResult IDNumberEngine::selectMolecules(const SelectResult &mols, bool uses_parallel) const
{
    auto matchMolecule = [&](const ViewsOfMol &mol)
    {
        const auto number = mol.data().number().value();
        return match(number, vals);
    };

    SelectResult::Container result;
    const auto m = mols.views();

    if (uses_parallel)
    {
        QVector<ViewsOfMol> views(m.count());
        
        tbb::parallel_for( tbb::blocked_range<int>(0,m.count()),
                           [&](const tbb::blocked_range<int> &r)
        {
            for (int i=r.begin(); i<r.end(); ++i)
            {
                const auto mol = m[i];
            
                if (matchMolecule(mol))
                    views[i] = ViewsOfMol(mol.molecule());
            }
        });
        
        for (const auto view : views)
        {
            if (not view.isEmpty())
                result.append(view);
        }
    }
    else
    {
        for (const auto mol : m)
        {
            if (matchMolecule(mol))
                result.append( ViewsOfMol(mol.molecule()) );
        }
    }
    
    return result;
}

SelectResult IDNumberEngine::select(const SelectResult &mols, const PropertyMap &map) const
{
    SelectResult::Container result;

    bool uses_parallel = true;
    
    if (map["parallel"].hasValue())
    {
        uses_parallel = map["parallel"].value().asA<BooleanProperty>().value();
    }
    
    switch(obj)
    {
    case AST::ATOM:
        return selectAtoms(mols, uses_parallel);
    case AST::RESIDUE:
        return selectResidues(mols, uses_parallel);
    case AST::MOLECULE:
        return selectMolecules(mols, uses_parallel);
    default:
        return SelectResult();
    }
}

SelectEngine::ObjType IDNumberEngine::objectType() const
{
    switch(obj)
    {
    case AST::ATOM:
        return SelectEngine::ATOM;
    case AST::CUTGROUP:
        return SelectEngine::CUTGROUP;
    case AST::RESIDUE:
        return SelectEngine::RESIDUE;
    case AST::CHAIN:
        return SelectEngine::CHAIN;
    case AST::SEGMENT:
        return SelectEngine::SEGMENT;
    case AST::MOLECULE:
        return SelectEngine::MOLECULE;
    default:
        return SelectEngine::COMPLEX;
    }
}

////////
//////// Implementation of the IDIndexEngine
////////

IDIndexEngine::IDIndexEngine()
{}

SelectEnginePtr IDIndexEngine::construct( IDObject o, RangeValues v )
{
    IDIndexEngine *ptr = new IDIndexEngine();
    ptr->obj = o;
    ptr->vals = v;
    
    return makePtr(ptr);
}

IDIndexEngine::~IDIndexEngine()
{}

SelectResult IDIndexEngine::selectAtoms(const SelectResult &mols, bool use_parallel) const
{
    return mols;
}

SelectResult IDIndexEngine::selectCutGroups(const SelectResult &mols, bool use_parallel) const
{
    return mols;
}

SelectResult IDIndexEngine::selectResidues(const SelectResult &mols, bool use_parallel) const
{
    return mols;
}

SelectResult IDIndexEngine::selectChains(const SelectResult &mols, bool use_parallel) const
{
    return mols;
}

SelectResult IDIndexEngine::selectSegments(const SelectResult &mols, bool use_parallel) const
{
    return mols;
}

SelectResult IDIndexEngine::selectMolecules(const SelectResult &mols, bool use_parallel) const
{
    return mols;
}

SelectResult IDIndexEngine::select(const SelectResult &mols, const PropertyMap &map) const
{
    SelectResult::Container result;

    bool uses_parallel = true;
    
    if (map["parallel"].hasValue())
    {
        uses_parallel = map["parallel"].value().asA<BooleanProperty>().value();
    }
    
    switch(obj)
    {
    case AST::ATOM:
        return selectAtoms(mols, uses_parallel);
    case AST::CUTGROUP:
        return selectCutGroups(mols, uses_parallel);
    case AST::RESIDUE:
        return selectResidues(mols, uses_parallel);
    case AST::CHAIN:
        return selectChains(mols, uses_parallel);
    case AST::SEGMENT:
        return selectSegments(mols, uses_parallel);
    case AST::MOLECULE:
        return selectMolecules(mols, uses_parallel);
    default:
        return SelectResult();
    }
}

SelectEngine::ObjType IDIndexEngine::objectType() const
{
    switch(obj)
    {
    case AST::ATOM:
        return SelectEngine::ATOM;
    case AST::CUTGROUP:
        return SelectEngine::CUTGROUP;
    case AST::RESIDUE:
        return SelectEngine::RESIDUE;
    case AST::CHAIN:
        return SelectEngine::CHAIN;
    case AST::SEGMENT:
        return SelectEngine::SEGMENT;
    case AST::MOLECULE:
        return SelectEngine::MOLECULE;
    default:
        return SelectEngine::COMPLEX;
    }
}

////////
//////// Implementation of the IDAndEngine
////////

IDAndEngine::IDAndEngine()
{}

SelectEnginePtr IDAndEngine::construct(SelectEnginePtr p0, SelectEnginePtr p1)
{
    IDAndEngine *ptr = new IDAndEngine();
    
    auto p = makePtr(ptr);
    
    if (p0)
        p0->setParent(p);

    if (p1)
        p1->setParent(p);

    ptr->part0 = p0;
    ptr->part1 = p1;
    
    return p;
}

IDAndEngine::~IDAndEngine()
{}

SelectResult IDAndEngine::select(const SelectResult &mols, const PropertyMap &map) const
{
    SelectResult::Container result;
    
    return result;
}

SelectEnginePtr IDAndEngine::simplify()
{
    if (part0.get())
        part0 = part0->simplify();
    
    if (part1.get())
        part1 = part1->simplify();
    
    if (part0.get() == 0)
        return part1;
    else if (part1.get() == 0)
        return part0;
    else
        return selfptr.lock();
}

SelectEngine::ObjType IDAndEngine::objectType() const
{
    if (part0.get() and part1.get())
    {
        auto o0 = part0->objectType();
        auto o1 = part1->objectType();
        
        if (o0 == o1)
            return o0;
        else
        {
            //the object type is always the smallest, e.g. atom and residue == atom
            return qMin(o0,o1);
        }
    }
    else if (part0.get())
    {
        return part0->objectType();
    }
    else if (part1.get())
    {
        return part1->objectType();
    }
    else
        return SelectEngine::COMPLEX;
}

////////
//////// Implementation of the IDOrEngine
////////

IDOrEngine::IDOrEngine()
{}

SelectEnginePtr IDOrEngine::construct(SelectEnginePtr part0, SelectEnginePtr part1)
{
    IDOrEngine *ptr = new IDOrEngine();
    
    auto p = makePtr(ptr);
    
    if (part0)
        part0->setParent(p);

    if (part1)
        part1->setParent(p);
    
    ptr->parts.append(part0);
    ptr->parts.append(part1);

    return p;
}

SelectEnginePtr IDOrEngine::construct(QList<SelectEnginePtr> parts)
{
    IDOrEngine *ptr = new IDOrEngine();
    
    auto p = makePtr(ptr);
    
    for (auto &part : parts)
    {
        if (part)
            part->setParent(p);
    }
    
    ptr->parts = parts;
    
    return p;
}

IDOrEngine::~IDOrEngine()
{}

SelectResult IDOrEngine::select(const SelectResult &mols, const PropertyMap &map) const
{
    SelectResult::Container result;
    
    return result;
}

SelectEnginePtr IDOrEngine::simplify()
{
    for (auto &part : parts)
    {
        part = part->simplify();
    }
    
    return selfptr.lock();
}

SelectEngine::ObjType IDOrEngine::objectType() const
{
    bool set = false;
    auto o = SelectEngine::COMPLEX;

    for (auto part : parts)
    {
        if (part.get())
        {
            if (not set)
            {
                o = part->objectType();
                set = true;
            }
            else if (o != part->objectType())
            {
                return SelectEngine::COMPLEX;
            }
        }
    }
    
    return o;
}

////////
//////// Implementation of the IDNotEngine
////////

IDNotEngine::IDNotEngine()
{}

SelectEnginePtr IDNotEngine::construct(SelectEnginePtr part)
{
    IDNotEngine *ptr = new IDNotEngine();
    auto p = makePtr(ptr);
    
    if (part)
        part->setParent(p);
    
    ptr->part = part;
    
    return p;
}

IDNotEngine::~IDNotEngine()
{}

SelectResult IDNotEngine::select(const SelectResult &mols, const PropertyMap &map) const
{
    SelectResult::Container result;
    
    return result;
}

SelectEnginePtr IDNotEngine::simplify()
{
    if (part.get())
        part = part->simplify();
    
    return selfptr.lock();
}

SelectEngine::ObjType IDNotEngine::objectType() const
{
    if (part)
        return part->objectType();
    else
        return SelectEngine::COMPLEX;
}

////////
//////// Implementation of the IDJoinEngine
////////

IDJoinEngine::IDJoinEngine()
{}

SelectEnginePtr IDJoinEngine::construct(SelectEnginePtr part)
{
    IDJoinEngine *ptr = new IDJoinEngine();
    auto p = makePtr(ptr);
    
    if (part)
        part->setParent(p);
    
    ptr->part = part;
    
    return p;
}

IDJoinEngine::~IDJoinEngine()
{}

SelectResult IDJoinEngine::select(const SelectResult &mols, const PropertyMap &map) const
{
    SelectResult::Container result;
    
    return result;
}

SelectEnginePtr IDJoinEngine::simplify()
{
    if (part.get())
        part = part->simplify();
    
    return selfptr.lock();
}

SelectEngine::ObjType IDJoinEngine::objectType() const
{
    return SelectEngine::COMPLEX;
}

////////
//////// Implementation of the IDSubScriptEngine
////////

IDSubScriptEngine::IDSubScriptEngine()
{}

SelectEnginePtr IDSubScriptEngine::construct(SelectEnginePtr part, const RangeValue &val)
{
    IDSubScriptEngine *ptr = new IDSubScriptEngine();
    auto p = makePtr(ptr);
    
    if (part)
        part->setParent(p);
    
    ptr->part = part;
    ptr->val = val;
    
    return p;
}

IDSubScriptEngine::~IDSubScriptEngine()
{}

SelectResult IDSubScriptEngine::select(const SelectResult &mols, const PropertyMap &map) const
{
    SelectResult::Container result;
    
    return result;
}

SelectEnginePtr IDSubScriptEngine::simplify()
{
    if (part.get())
        part = part->simplify();
    
    return selfptr.lock();
}

SelectEngine::ObjType IDSubScriptEngine::objectType() const
{
    if (part)
        return part->objectType();
    else
        return SelectEngine::COMPLEX;
}

////////
//////// Implementation of the IDWithEngine
////////

IDWithEngine::IDWithEngine()
{}

SelectEnginePtr IDWithEngine::construct( IDObject obj, IDToken token, SelectEnginePtr part)
{
    IDWithEngine *ptr = new IDWithEngine();
    auto p = makePtr(ptr);
    
    if (part)
        part->setParent(p);
    
    ptr->part = part;
    ptr->obj = obj;
    ptr->token = token;
    
    return p;
}

IDWithEngine::~IDWithEngine()
{}

SelectResult IDWithEngine::select(const SelectResult &mols, const PropertyMap &map) const
{
    SelectResult::Container result;
    
    return result;
}

SelectEnginePtr IDWithEngine::simplify()
{
    if (part.get())
        part = part->simplify();
    
    return selfptr.lock();
}

SelectEngine::ObjType IDWithEngine::objectType() const
{
    switch(obj)
    {
    case AST::ATOM:
        return SelectEngine::ATOM;
    case AST::CUTGROUP:
        return SelectEngine::CUTGROUP;
    case AST::RESIDUE:
        return SelectEngine::RESIDUE;
    case AST::CHAIN:
        return SelectEngine::CHAIN;
    case AST::SEGMENT:
        return SelectEngine::SEGMENT;
    case AST::MOLECULE:
        return SelectEngine::MOLECULE;
    default:
        return SelectEngine::COMPLEX;
    }
}
