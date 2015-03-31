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

#include "molgroupworkspace.h"

#include "SireBase/majorminorversion.h"

#include "SireMol/moleculegroup.h"
#include "SireMol/molecule.h"
#include "SireMol/viewsofmol.h"
#include "SireMol/partialmolecule.h"

#include "SireError/errors.h"

#include <QDebug>
#include <QVarLengthArray>
#include <QThreadStorage>

using namespace SireMol;

namespace SireMol
{
    namespace detail
    {
        class MolGroupWorkspaceData
        {
        public:
            MolGroupWorkspaceData()
            {
                //qDebug() << "MolGroupWorkspaceData::MolGroupWorkspaceData()" << quintptr(this);
            }
            
            MolGroupWorkspaceData(const MolGroupWorkspaceData &other)
                    : ver(other.ver), mols(other.mols), views(other.views)
            {
                //qDebug() << "MolGroupWorkspaceData::MolGroupWorkspaceData(copy)"
                //         << quintptr(this) << quintptr(&other);
            }
            
            ~MolGroupWorkspaceData()
            {
                //qDebug() << "MolGroupWorkspaceData::~MolGroupWorkspaceData()" << quintptr(this);
            }
            
            MolGroupWorkspaceData& operator=(const MolGroupWorkspaceData &other)
            {
                ver = other.ver;
                mols = other.mols;
                views = other.views;
                return *this;
            }
            
            void clear()
            {
                mols.resize(0);
                views.clear();
                ver = MajorMinorVersion();
            }
            
            bool contains(const MolNum &molnum) const
            {
                for (int i=0; i<mols.count(); ++i)
                {
                    if (mols.constData()[i].read().number() == molnum)
                        return true;
                }
                
                return false;
            }
            
            bool isEmpty() const
            {
                return mols.isEmpty() and ver.majorVersion() == 0 and ver.minorVersion() == 0;
            }
            
            const SharedDataPointer<MoleculeData>* data() const
            {
                return mols.data();
            }
            
            int count() const
            {
                return mols.count();
            }
            
            void push(const Molecule &molecule)
            {
                mols.push_back(molecule.data());
            }
            
            const ViewsOfMol& getUpdated(const ViewsOfMol &oldmol)
            {
                for (int i=0; i<mols.count(); ++i)
                {
                    const MoleculeData &newmol = mols.constData()[i].read();
                    
                    if (newmol.number() == oldmol.number())
                    {
                        ViewsOfMol mol(oldmol);
                        mol.update(newmol);
                        views.push_back(mol);
                        //qDebug() << "RETURNING UPDATED";
                        return views.last();
                    }
                }
                
                //qDebug() << "RETURNING ORIGINAL?";
                return oldmol;
            }
            
            void setVersion(const MajorMinorVersion &version)
            {
                ver = version;
            }
            
            const MajorMinorVersion& version() const
            {
                return ver;
            }
            
            void incrementMajor()
            {
                if (not isEmpty())
                    ver.incrementMajor();
            }
            
            void incrementMinor()
            {
                if (not isEmpty())
                    ver.incrementMinor();
            }
            
        private:
            MajorMinorVersion ver;
            QVarLengthArray<SireBase::SharedDataPointer<MoleculeData>,10> mols;
            QVarLengthArray<ViewsOfMol,10> views;
        };
    }
}

typedef QVarLengthArray<boost::shared_ptr<SireMol::detail::MolGroupWorkspaceData>,32>
                                                                        MolGroupWorkspaceCache;

static QThreadStorage<MolGroupWorkspaceCache*> cache;

/** Call this to return the memory allocated in this object back to the memory pool */
void MolGroupWorkspace::returnToMemoryPool()
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
        cache.setLocalData( new MolGroupWorkspaceCache() );
    }
    
    if (cache.localData()->count() < 32)
    {
        d->clear();
        //qDebug() << "PUSHED" << quintptr(d.get()) << "INTO THE POOL";
        cache.localData()->push_back(d);
    }
    else
        qDebug() << "DELETING AS THE CACHE IS FULL";
    
    d.reset();
}

/** Construct this by using memory from the pool */
void MolGroupWorkspace::createFromMemoryPool()
{
    if (d.get() != 0)
    {
        d->clear();
    }
    else
    {
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
        d.reset( new SireMol::detail::MolGroupWorkspaceData() );
    }
}

void MolGroupWorkspace::detach()
{
    if (d.get() != 0)
    {
        if (not d.unique())
        {
            boost::shared_ptr<detail::MolGroupWorkspaceData> d2 = d;
            
            d.reset();
            createFromMemoryPool();
            d->operator=(*d2);
        }
    }
}

/** Constructor */
MolGroupWorkspace::MolGroupWorkspace()
{}

/** Copy constructor */
MolGroupWorkspace::MolGroupWorkspace(const MolGroupWorkspace &other)
{
    if (not other.isEmpty())
    {
        d = other.d;
    }
}

/** Destructor */
MolGroupWorkspace::~MolGroupWorkspace()
{
    returnToMemoryPool();
}

/** Copy assignment operator */
MolGroupWorkspace& MolGroupWorkspace::operator=(const MolGroupWorkspace &other)
{
    if (d.get() != other.d.get())
    {
        returnToMemoryPool();
        d = other.d;
    }
    
    return *this;
}

/** Comparison operator */
bool MolGroupWorkspace::operator==(const MolGroupWorkspace &other) const
{
    return (this->isEmpty() and other.isEmpty()) or
           d.get() == other.d.get();
}

/** Comparison operator */
bool MolGroupWorkspace::operator!=(const MolGroupWorkspace &other) const
{
    return not operator==(other);
}

/** Return whether or not this workspace is empty */
bool MolGroupWorkspace::isEmpty() const
{
    if (d.get())
        return d->isEmpty();
    else
        return true;
}

/** Return the number of molecules stored in this workspace */
int MolGroupWorkspace::count() const
{
    if (d.get())
        return d->count();
    else
        return 0;
}

/** Return the number of molecules stored in this workspace */
int MolGroupWorkspace::size() const
{
    return count();
}

/** Return the ith molecule stored in this workspace

    \throw SireError::invalid_index
*/
Molecule MolGroupWorkspace::operator[](int i) const
{
    if (i < 0 or i >= this->count())
        throw SireError::invalid_index( QObject::tr(
                "Cannot return the Molecule at index %1 as the number of molecules is %2.")
                    .arg(i).arg(count()), CODELOC );

    return Molecule(d->data()[i].read());
}

/** Return the ith molecule stored in this workspace

    \throw SireError::invalid_index
*/
Molecule MolGroupWorkspace::at(int i) const
{
    return operator[](i);
}

/** Return the ith molecule stored in this workspace

    \throw SireError::invalid_index
*/
Molecule MolGroupWorkspace::getitem(int i) const
{
    return operator[](i);
}

/** Return a raw pointer to the array for the molecules in the workspace */
const SireBase::SharedDataPointer<MoleculeData>* MolGroupWorkspace::data() const
{
    if (d.get())
        return d->data();
    else
        return 0;
}

/** Return a raw pointer to the array for the molecules in the workspace */
const SireBase::SharedDataPointer<MoleculeData>* MolGroupWorkspace::constData() const
{
    return data();
}

/** Push a molecule into this workspace */
void MolGroupWorkspace::push(const MoleculeData &molecule)
{
    detach();

    if (not d.get())
    {
        this->createFromMemoryPool();
    }
    
    d->push(molecule);
}

/** Return the updated view of the passed molecule */
const ViewsOfMol& MolGroupWorkspace::getUpdated(const ViewsOfMol &oldmol) const
{
    //qDebug() << "MolGroupWorkspace::getUpdated(" << oldmol.number().value() << ")";

    if (d.get())
    {
        if (d->contains(oldmol.number()))
        {
            //qDebug() << "DOING SOME WORK";
            MolGroupWorkspace *nonconst_this = const_cast<MolGroupWorkspace*>(this);
            nonconst_this->detach();
            return nonconst_this->d->getUpdated(oldmol);
        }
    }

    return oldmol;
}

/** Return the updated view of the passed molecule */
PartialMolecule MolGroupWorkspace::getUpdated(const PartialMolecule &oldmol) const
{
    if (d.get())
    {
        if (d->contains(oldmol.number()))
        {
            for (int i=0; i<d->count(); ++i)
            {
                const MoleculeData &mol = data()[i].read();
                
                if (mol.number() == oldmol.number())
                {
                    PartialMolecule newmol(oldmol);
                    newmol.update(mol);
                    return newmol;
                }
            }
        }
    }
    
    return oldmol;
}

/** Return the number of molecules stored in this workspace */
int MolGroupWorkspace::nMolecules() const
{
    return count();
}

/** Clear this workspace */
void MolGroupWorkspace::clear()
{
    returnToMemoryPool();
}

/** Set the version object used to track the version of a MoleculeGroup */
void MolGroupWorkspace::setVersion(const MajorMinorVersion &version)
{
    if (d.get())
    {
        detach();
        d->setVersion(version);
    }
}

/** Return the version of the current workspace */
MajorMinorVersion MolGroupWorkspace::version() const
{
    if (d.get())
        return d->version();
    else
        return MajorMinorVersion();
}

/** Increment the major version */
void MolGroupWorkspace::incrementMajor()
{
    if (d.get())
    {
        detach();
        d->incrementMajor();
    }
}

/** Increment the minor version */
void MolGroupWorkspace::incrementMinor()
{
    if (d.get())
    {
        detach();
        d->incrementMinor();
    }
}
