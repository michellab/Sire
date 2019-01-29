/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2014  Christopher Woods
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

#include "cljworkspace.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QVarLengthArray>
#include <QThreadStorage>

using namespace SireMM;
using namespace SireStream;
using namespace SireMol;

namespace SireMM
{
namespace detail
{
    /** This class holds all of the data of a CLJWorkspace. This is held
        as an implicitly shared pointer, with the caveat that copying the 
        pointer from an empty workspace will reset the original to point
        to an empty copy, while taking this copy to the detached (new) 
        object
        
        @author Christopher Woods
    */
    class CLJWorkspaceData
    {
    public:
        CLJWorkspaceData()
        {}
        
        CLJWorkspaceData(const CLJWorkspaceData &other)
                : deltas(other.deltas), same_ids(other.same_ids)
        {}
        
        ~CLJWorkspaceData()
        {}
        
        CLJWorkspaceData& operator=(const CLJWorkspaceData &other)
        {
            if (this != &other)
            {
                deltas = other.deltas;
                same_ids = other.same_ids;
            }
            
            return *this;
        }
        
        QVarLengthArray<CLJDelta,4> deltas;
        
        QVarLengthArray<QVector<CLJBoxIndex>, 4> same_ids;
        
        void clear()
        {
            deltas.resize(0);
            same_ids.resize(0);
        }
        
        bool isEmpty() const
        {
            return deltas.isEmpty() and same_ids.isEmpty();
        }
        
        /** Return whether or not all of the atoms in the deltas
            have a single ID. If they have a single ID, then we don't
            need to remove the old atoms before calculating the energy
            delta (as atoms with the same ID don't interact with one another) */
        bool isSingleID() const
        {
            if (isEmpty())
                return true;
            
            qint32 id = CLJAtoms::idOfDummy()[0];
            const qint32 dummy_id = id;
            
            for (int i=0; i<deltas.count(); ++i)
            {
                //loop over the old and new atoms
                {
                    const QVector<MultiInt> &ids = deltas.constData()[i].newAtoms().ID();
                
                    for (int j=0; j<ids.count(); ++j)
                    {
                        for (int k=0; k<MultiInt::count(); ++k)
                        {
                            if (ids[j][k] != dummy_id)
                            {
                                if (id == dummy_id)
                                    id = ids[j][k];

                                else if (id != ids[j][k])
                                {
                                    //we found an atom with a different ID
                                    return false;
                                }
                            }
                        }
                    }
                }

                //loop over the old and new atoms
                {
                    const QVector<MultiInt> &ids = deltas.constData()[i].oldAtoms().ID();
                
                    for (int j=0; j<ids.count(); ++j)
                    {
                        for (int k=0; k<MultiInt::count(); ++k)
                        {
                            if (ids[j][k] != dummy_id)
                            {
                                if (id == dummy_id)
                                    id = ids[j][k];

                                else if (id != ids[j][k])
                                {
                                    //we found an atom with a different ID
                                    return false;
                                }
                            }
                        }
                    }
                }
            }
            
            return true;
        }

        void removeSameIDAtoms(CLJBoxes &boxes)
        {
            for (int i=0; i<same_ids.count(); ++i)
            {
                boxes.remove(same_ids.at(i));
            }
            
            same_ids.clear();
        }

        void mustRecalculateFromScratch(CLJBoxes &boxes)
        {
            removeSameIDAtoms(boxes);
            deltas.clear();
        }

        void accept(CLJBoxes &boxes)
        {
            removeSameIDAtoms(boxes);
            
            for (int i=0; i<deltas.count(); ++i)
            {
                if (not deltas.at(i).isNull())
                {
                    throw SireError::program_bug( QObject::tr(
                            "Accepting a CLJWorkspace when some of the deltas have not "
                            "been applied: %1 %2. This is likely to cause an energy leak!")
                                .arg(i).arg(deltas.at(i).toString()), CODELOC );
                }
            }
            
            deltas.clear();
        }

        CLJDelta push(CLJBoxes &boxes, const QVector<CLJBoxIndex> &old_atoms,
                      const CLJAtoms &new_atoms, const CLJDelta &old_delta)
        {
            if (deltas.isEmpty())
            {
                if (not old_delta.isNull())
                    throw SireError::program_bug( QObject::tr(
                            "How can we have an old delta (%1) when the current "
                            "CLJWorkspace is empty? %2")
                                .arg(old_delta.toString()).arg(deltas.count()), CODELOC );
            
                //this is the first change being made to the boxes
                CLJAtoms old_cljatoms = boxes.get(old_atoms);
                
                deltas.append( CLJDelta(0, old_cljatoms, new_atoms) );
                
                if (this->isSingleID())
                {
                    //we don't need to remove the old atoms from the
                    //CLJBoxes as they will not interfere with any calculation
                    same_ids.append( old_atoms );
                }
                else
                {
                    //we do need to remove the old atoms from the
                    //CLJBoxes as they will interfere with the calculation
                    boxes.remove(old_atoms);
                }
                
                return deltas.last();
            }
            else if (not same_ids.isEmpty())
            {
                //this is another change to the boxes, but all changes
                //so far have involved atoms with the same ID. Is this
                //going to continue?
                CLJAtoms old_cljatoms = boxes.get(old_atoms);
                
                CLJDelta new_delta;
                
                if (old_delta.isNull())
                {
                    new_delta = CLJDelta(deltas.count(), old_cljatoms, new_atoms);
                    deltas.append(new_delta);
                }
                else
                {
                    if (old_delta.ID() < 0 or old_delta.ID() >= deltas.count())
                        throw SireError::program_bug( QObject::tr(
                                "How can we have the CLJDelta %1 with an ID of %2 "
                                "when this CLJWorkspace only has %3 deltas?")
                                    .arg(old_delta.toString()).arg(old_delta.ID())
                                        .arg(deltas.count()), CODELOC );
                
                    if (old_delta.oldAtoms() != old_cljatoms)
                    {
                        QStringList differences;
                        
                        for (int i=0; i<qMin(old_delta.oldAtoms().count(), old_cljatoms.count());
                                ++i)
                        {
                            CLJAtom old0 = old_delta.oldAtoms().at(i);
                            CLJAtom old1 = old_cljatoms.at(i);
                            
                            if (old0 != old1)
                            {
                                differences.append( QObject::tr(
                                    "Disagreement of atom %1: %2 vs. %3")
                                        .arg(i).arg(old0.toString(),old1.toString()) );
                            }
                        }
                        
                        if (old_delta.oldAtoms().nAtoms() != old_cljatoms.nAtoms())
                        {
                            differences.append( QObject::tr(
                                    "Disagreement of number of atoms: %1 vs. %2")
                                        .arg(old_delta.oldAtoms().nAtoms())
                                        .arg(old_cljatoms.nAtoms()) );
                        }
                        
                        if (not differences.isEmpty())
                            throw SireError::program_bug( QObject::tr(
                                    "Disagreement in the old atoms when changing "
                                    "the CLJAtoms for a second time:\n%1\n")
                                        .arg(differences.join("\n")), CODELOC );
                    }
                    
                    new_delta = CLJDelta(old_delta.ID(), old_cljatoms, new_atoms);
                    deltas[old_delta.ID()] = new_delta;
                }

                if (this->isSingleID())
                {
                    //yes, this is another change to the same ID
                    same_ids.append( old_atoms );
                }
                else
                {
                    //no, we are now changing ID. Need to remove all of the
                    //old atoms
                    boxes.remove(old_atoms);
                    removeSameIDAtoms(boxes);
                }
                
                return new_delta;
            }
            else
            {
                //we have to remove the atoms as changes now involve
                //more than one ID number
                if (old_delta.isNull())
                {
                    //we have not moved this group of CLJAtoms before
                    CLJAtoms old_cljatoms = boxes.take(old_atoms);
                    deltas.append( CLJDelta(deltas.count(),old_cljatoms,new_atoms) );
                    return deltas.last();
                }
                else
                {
                    //we have changed this group of CLJAtoms before.
                    //Because the old atoms have been removed from the CLJBoxes
                    //we need to trust that the ones stored in the delta are correct...
                    if (old_delta.ID() < 0 or old_delta.ID() >= deltas.count())
                        throw SireError::program_bug( QObject::tr(
                                "How can we have the CLJDelta %1 with an ID of %2 "
                                "when this CLJWorkspace only has %3 deltas?")
                                    .arg(old_delta.toString()).arg(old_delta.ID())
                                        .arg(deltas.count()), CODELOC );

                    CLJDelta new_delta(old_delta.ID(), old_delta.oldAtoms(), new_atoms);
                    deltas[old_delta.ID()] = new_delta;
                    return new_delta;
                }
            }
        }
        
        tuple<CLJAtoms,CLJAtoms,CLJAtoms> merge() const
        {
            if (deltas.isEmpty())
                return tuple<CLJAtoms,CLJAtoms,CLJAtoms>();
            else
            {
                return CLJDelta::merge(deltas.constData(), deltas.count());
            }
        }
        
        CLJAtoms changedAtoms() const
        {
            if (deltas.isEmpty())
                return CLJAtoms();
            else if (deltas.count() == 1)
                return deltas.constData()[0].changedAtoms();
            else
            {
                return CLJDelta::mergeChanged(deltas.constData(), deltas.count());
            }
        }

        CLJAtoms newAtoms() const
        {
            if (deltas.isEmpty())
                return CLJAtoms();
            else
            {
                return CLJDelta::mergeNew(deltas.constData(), deltas.count());
            }
        }

        CLJAtoms oldAtoms() const
        {
            if (deltas.isEmpty())
                return CLJAtoms();
            else
            {
                return CLJDelta::mergeOld(deltas.constData(), deltas.count());
            }
        }
        
        bool needsAccepting() const
        {
            if (not same_ids.isEmpty())
                return true;
            else
            {
                for (int i=0; i<deltas.count(); ++i)
                {
                    if (not deltas.at(i).isNull())
                        return true;
                }
                
                return false;
            }
        }
        
        QVector<CLJBoxIndex> commit(CLJBoxes &boxes, const CLJDelta &delta)
        {
            removeSameIDAtoms(boxes);
            
            if (delta.ID() < 0 or delta.ID() >= deltas.count())
                //invalid delta, so return nothing
                return QVector<CLJBoxIndex>();
            
            else if (deltas.at(delta.ID()).isNull())
                //we have already committed this delta
                return QVector<CLJBoxIndex>();
            
            else
            {
                //make sure that we have agreement regarding the delta...
                deltas.at(delta.ID()).assertIdenticalTo(delta);
                deltas[delta.ID()] = CLJDelta();
                return boxes.add(delta.newAtoms());
            }
        }
        
        QVector<CLJBoxIndex> revert(CLJBoxes &boxes, const CLJDelta &delta)
        {
            removeSameIDAtoms(boxes);
            
            if (delta.ID() < 0 or delta.ID() >= deltas.count())
                //invalid delta, so return nothing
                return QVector<CLJBoxIndex>();
            
            else if (deltas.at(delta.ID()).isNull())
                //we have already committed this delta
                return QVector<CLJBoxIndex>();
            
            else
            {
                //make sure that we have agreement regarding the delta...
                deltas.at(delta.ID()).assertIdenticalTo(delta);
                deltas[delta.ID()] = CLJDelta();
                return boxes.add(delta.oldAtoms());
            }
        }
    };
} // end of namespace detail
}

static const RegisterMetaType<CLJWorkspace> r_workspace(NO_ROOT);

QDataStream  &operator<<(QDataStream &ds, const SireMM::detail::CLJWorkspaceData &ws)
{
    ds << quint32( ws.deltas.count() );
    
    for (int i=0; i<ws.deltas.count(); ++i)
    {
        ds << ws.deltas.at(i);
    }
    
    return ds;
}

QDataStream  &operator>>(QDataStream &ds, SireMM::detail::CLJWorkspaceData &ws)
{
    quint32 n;
    ds >> n;
    
    ws.deltas.resize(n);
    
    for (quint32 i=0; i<n; ++i)
    {
        ds >> ws.deltas[i];
    }
    
    return ds;
}

QDataStream &operator<<(QDataStream &ds, const CLJWorkspace &ws)
{
    writeHeader(ds, r_workspace, 2);
    
    ds << ws.recalc_from_scratch << ws.isEmpty();
    
    if (not ws.isEmpty())
    {
        SharedDataStream sds(ds);
        sds << ws.d;
    }
    
    return ds;
}

QDataStream &operator>>(QDataStream &ds, CLJWorkspace &ws)
{
    VersionID v = readHeader(ds, r_workspace);
    
    if (v <= 2)
    {
        bool recalc_from_scratch = false;
        bool is_empty = true;
        
        if (v == 2)
            ds >> recalc_from_scratch;
        
        ds >> is_empty;
        
        if (is_empty)
        {
            ws = CLJWorkspace();
            ws.recalc_from_scratch = recalc_from_scratch;
        }
        else
        {
            SharedDataStream sds(ds);
            sds >> ws.d;
            ws.recalc_from_scratch = recalc_from_scratch;
        }
    }
    else
        throw version_error(v, "1,2", r_workspace, CODELOC);
    
    return ds;
}

typedef QVarLengthArray<boost::shared_ptr<SireMM::detail::CLJWorkspaceData>,32> CLJWorkspaceCache;
static QThreadStorage<CLJWorkspaceCache*> cache;

/** Call this to return the memory allocated in this object back to the memory pool */
void CLJWorkspace::returnToMemoryPool()
{
    if (d.get() == 0)
        return;
    
    if (not d.unique())
    {
        d.reset();
        return;
    }
    
    if (not cache.hasLocalData())
    {
        cache.setLocalData( new CLJWorkspaceCache() );
    }
    
    if (cache.localData()->count() < 32)
    {
        d->clear();
        cache.localData()->push_back(d);
    }
    
    d.reset();
}

/** Construct this by using memory from the pool */
void CLJWorkspace::createFromMemoryPool()
{
    if (d.get() != 0)
    {
        if (d.unique())
        {
            d->clear();
            return;
        }
        else
            d.reset();
    }
    
    if (cache.hasLocalData())
    {
        if (not cache.localData()->isEmpty())
        {
            d = cache.localData()->back();
            cache.localData()->pop_back();
            //qDebug() << "TAKEN" << quintptr(d.get()) << "FROM THE POOL" << d->isEmpty();
            return;
        }
    }
    
    //no available value in the pool
    //qDebug() << "CREATING AS NOT AVAILABLE IN THE POOL";
    d.reset( new SireMM::detail::CLJWorkspaceData() );
}

/** Constructor */
CLJWorkspace::CLJWorkspace() : recalc_from_scratch(false)
{}

/** Copy constructor */
CLJWorkspace::CLJWorkspace(const CLJWorkspace &other)
             : recalc_from_scratch(other.recalc_from_scratch)
{
    if (not other.isEmpty())
    {
        d = other.d;
    }
}

/** Destructor */
CLJWorkspace::~CLJWorkspace()
{
    returnToMemoryPool();
}

/** Copy assignment operator. The new copy will get the memory of an empty workspace */
CLJWorkspace& CLJWorkspace::operator=(const CLJWorkspace &other)
{
    if (d.get() != other.d.get())
    {
        returnToMemoryPool();
        d = other.d;
    }

    recalc_from_scratch = other.recalc_from_scratch;
    
    return *this;
}

/** Comparison operator */
bool CLJWorkspace::operator==(const CLJWorkspace &other) const
{
    return (recalc_from_scratch == other.recalc_from_scratch) and
             ((isEmpty() and other.isEmpty()) or (d.get() == other.d.get()));
}

/** Comparison operator */
bool CLJWorkspace::operator!=(const CLJWorkspace &other) const
{
    return not operator==(other);
}

const char* CLJWorkspace::typeName()
{
    return QMetaType::typeName( qMetaTypeId<CLJWorkspace>() );
}

const char* CLJWorkspace::what() const
{
    return CLJWorkspace::typeName();
}

QString CLJWorkspace::toString() const
{
    if (this->isEmpty())
        return QObject::tr("CLJWorkspace::empty");
    else
        return QObject::tr("CLJWorkspace( nDeltas() == %1 )")
                    .arg(nDeltas());
}

/** Return the ith delta */
const CLJDelta& CLJWorkspace::operator[](int i) const
{
    if (isEmpty())
        throw SireError::invalid_index( QObject::tr(
                "Cannot return the delta at index %1 at this is an empty workspace.")
                        .arg(i), CODELOC );
    
    if (i < 0 or i >= d->deltas.count())
        throw SireError::invalid_index( QObject::tr(
                "Cannot return the delta at index %1 as there number of deltas is only %2.")
                    .arg(i).arg(d->deltas.count()), CODELOC );
    
    return d->deltas.at(i);
}

/** Return the ith delta */
const CLJDelta& CLJWorkspace::at(int i) const
{
    return operator[](i);
}

/** Return the ith delta */
CLJDelta CLJWorkspace::getitem(int i) const
{
    return operator[](i);
}

/** Return a raw pointer to the array of deltas */
const CLJDelta* CLJWorkspace::data() const
{
    if (d.get() == 0)
        return 0;
    else
        return d->deltas.constData();
}

/** Return a raw pointer to the array of deltas */
const CLJDelta* CLJWorkspace::constData() const
{
    return data();
}

/** Return the number of deltas */
int CLJWorkspace::count() const
{
    if (d.get() == 0)
        return 0;
    else
        return d->deltas.count();
}

/** Return the number of deltas */
int CLJWorkspace::size() const
{
    return count();
}

void CLJWorkspace::detach()
{
    if (d.get() != 0)
    {
        if (not d.unique())
        {
            boost::shared_ptr<detail::CLJWorkspaceData> d2 = d;
            d.reset();
            createFromMemoryPool();
            
            d->operator=(*d2);
        }
    }
}

/** Internal function used to fully remove atoms that have not been
    removed because they all have the same ID. This is used when the 
    optimisation of not removing same ID atoms would break the energy
    calculation, e.g. if we have multiple CLJGroups and have multiple
    changed IDs across these groups */
void CLJWorkspace::removeSameIDAtoms(CLJBoxes &boxes)
{
    if (d.get() != 0)
    {
        detach();
        d->removeSameIDAtoms(boxes);
    }
}

/** Push the passed change onto the workspace. This changes the atoms in the 
    passed CLJBoxes from their values in 'old_atoms' to their values in the
    passed 'new_atoms'. The last CLJDelta used for these atoms is supplied as
    'old_delta', and this returns the new CLJDelta for these atoms. */
CLJDelta CLJWorkspace::push(CLJBoxes &boxes, const QVector<CLJBoxIndex> &old_atoms,
                            const CLJAtoms &new_atoms, const CLJDelta &old_delta)
{
    if (not recalc_from_scratch)
    {
        detach();
        
        if (d.get() == 0)
        {
            createFromMemoryPool();
        }
        
        return d->push(boxes, old_atoms, new_atoms, old_delta);
    }
    else
    {
        if (old_delta.isNull())
        {
            CLJAtoms old_cljatoms = boxes.take(old_atoms);
            return CLJDelta(0,old_cljatoms,new_atoms);
        }
        else
        {
            return CLJDelta(old_delta.ID(), old_delta.oldAtoms(), new_atoms);
        }
    }
}

/** Commit the changes in the passed delta into the passed CLJBoxes */
QVector<CLJBoxIndex> CLJWorkspace::commit(CLJBoxes &boxes, const CLJDelta &delta)
{
    if (d.get() == 0)
    {
        return boxes.add(delta.newAtoms());
    }
    else
    {
        detach();
        return d->commit(boxes, delta);
    }
}

/** Revert the changes supplied in the passed delta */
QVector<CLJBoxIndex> CLJWorkspace::revert(CLJBoxes &boxes, const CLJDelta &delta)
{
    if (d.get() == 0)
    {
        return boxes.add(delta.oldAtoms());
    }
    else
    {
        detach();
        return d->revert(boxes, delta);
    }
}

/** Tell the workspace that we are now recalculating everything from scratch */
void CLJWorkspace::mustRecalculateFromScratch(CLJBoxes &boxes)
{
    recalc_from_scratch = true;

    if (d.get() != 0)
    {
        detach();
        d->mustRecalculateFromScratch(boxes);
        returnToMemoryPool();
    }
}

/** Return whether or not we are recalculating everything from scratch */
bool CLJWorkspace::recalculatingFromScratch() const
{
    return recalc_from_scratch;
}

/** Return whether or not this workspace needs accepting */
bool CLJWorkspace::needsAccepting() const
{
    if (d.get() != 0)
    {
        return d->needsAccepting();
    }
    else
        return recalc_from_scratch;
}

/** Accept this workspace - this clears the 'recalc_from_scratch' flag as it
    signals that we have put the CLJBoxes into a sane state */
void CLJWorkspace::accept(CLJBoxes &boxes)
{
    recalc_from_scratch = false;

    if (d.get() != 0)
    {
        detach();
        d->accept(boxes);
        returnToMemoryPool();
    }
}

/** Return the number of deltas in this workspace */
int CLJWorkspace::nDeltas() const
{
    return count();
}

/** Return whether or not this workspace contains deltas with a single
    ID (a single CLJAtoms ID) */
bool CLJWorkspace::isSingleID() const
{
    if (d.get() == 0)
        return false;
    else
        return d->isSingleID();
}

/** Return whether or not this workspace is empty */
bool CLJWorkspace::isEmpty() const
{
    if (d.get() == 0)
        return true;
    else
        return d->isEmpty();
}

/** Merge all of the deltas together into a single set of changed CLJAtoms that
    can be used for the change in energy calculation */
CLJAtoms CLJWorkspace::changedAtoms() const
{
    if (isEmpty())
        return CLJAtoms();
    else
        return d->changedAtoms();
}

/** Merge all of the old atoms from the deltas together into a single 
    set of old CLJAtoms that can be used for the change in energy calculation */
CLJAtoms CLJWorkspace::oldAtoms() const
{
    if (isEmpty())
        return CLJAtoms();
    else
        return d->oldAtoms();
}

/** Merge all of the new atoms from the deltas together into a single
    set of new CLJAtoms that can be used for the change in energy calculation */
CLJAtoms CLJWorkspace::newAtoms() const
{
    if (isEmpty())
        return CLJAtoms();
    else
        return d->newAtoms();
}

/** Merge all of the deltas together to return the tuple of the 
    changed, old and new atoms. This is equivalent to calling
    changedAtoms(), oldAtoms() and newAtoms() and placing them
    into a tuple. This is more efficient than three separate calls */
tuple<CLJAtoms,CLJAtoms,CLJAtoms> CLJWorkspace::merge() const
{
    if (isEmpty())
        return tuple<CLJAtoms,CLJAtoms,CLJAtoms>();
    else
        return d->merge();
}

/** Clear this workspace */
void CLJWorkspace::clear()
{
    returnToMemoryPool();
    recalc_from_scratch = true;
}
