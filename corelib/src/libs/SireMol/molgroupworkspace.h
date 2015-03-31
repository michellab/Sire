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

#ifndef SIREMOL_MOLGROUPWORKSPACE_H
#define SIREMOL_MOLGROUPWORKSPACE_H

#include "SireBase/shareddatapointer.hpp"

#include <boost/shared_ptr.hpp>

SIRE_BEGIN_HEADER

namespace SireBase
{
class MajorMinorVersion;
}

namespace SireMol
{

using SireBase::MajorMinorVersion;

class Molecule;
class PartialMolecule;
class ViewsOfMol;
class MoleculeData;

namespace detail{ class MolGroupWorkspaceData; }

/** This is a utility class that is used by MoleculeGroup as a temporary 
    workspace used when updating molecules held in the group. This class
    is not part of the public API and is used to allow memory 
    allocation / deallocation to be minimised during a
    Monte Carlo simulation
    
    @author Christopher Woods
*/
class MolGroupWorkspace
{
public:
    MolGroupWorkspace();
    MolGroupWorkspace(const MolGroupWorkspace &other);
    
    ~MolGroupWorkspace();
    
    MolGroupWorkspace& operator=(const MolGroupWorkspace &other);
    
    bool operator==(const MolGroupWorkspace &other) const;
    bool operator!=(const MolGroupWorkspace &other) const;
    
    bool isEmpty() const;
    
    int count() const;
    int size() const;
    
    Molecule operator[](int i) const;
    
    Molecule at(int i) const;
    
    Molecule getitem(int i) const;
    
    const SireBase::SharedDataPointer<MoleculeData>* data() const;
    const SireBase::SharedDataPointer<MoleculeData>* constData() const;
    
    void push(const MoleculeData &molecule);
    
    const ViewsOfMol& getUpdated(const ViewsOfMol &old_molecule) const;
    PartialMolecule getUpdated(const PartialMolecule &old_molecule) const;
    
    int nMolecules() const;
    
    void clear();
    
    MajorMinorVersion version() const;
    
    void incrementMinor();
    void incrementMajor();
    
    void setVersion(const MajorMinorVersion &version);
    
private:
    void returnToMemoryPool();
    void createFromMemoryPool();
    
    void detach();

    boost::shared_ptr<detail::MolGroupWorkspaceData> d;
};

}

SIRE_END_HEADER

#endif
