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

#ifndef SIREMM_CLJGROUP_H
#define SIREMM_CLJGROUP_H

#include "cljextractor.h"

#include "SireMol/moleculegroup.h"
#include "SireBase/chunkedhash.hpp"

SIRE_BEGIN_HEADER

namespace SireMM
{
class CLJGroup;
}

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::CLJGroup&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::CLJGroup&);

namespace SireMM
{

using SireMol::MolNum;
using SireMol::MoleculeGroup;

/** This class holds and manages a group of molecules that have been 
    added to a CLJBoxes object. 
    
    @author Christopher Woods
*/
class SIREMM_EXPORT CLJGroup
{

friend SIREMM_EXPORT QDataStream& ::operator<<(QDataStream&, const SireMM::CLJGroup&);
friend SIREMM_EXPORT QDataStream& ::operator>>(QDataStream&, SireMM::CLJGroup&);

public:
    CLJGroup();
    CLJGroup(CLJAtoms::ID_SOURCE id_source);
    CLJGroup(CLJExtractor::EXTRACT_SOURCE extract_source);
    CLJGroup(CLJAtoms::ID_SOURCE id_source, CLJExtractor::EXTRACT_SOURCE extract_source);

    CLJGroup(const CLJGroup &other);
    
    ~CLJGroup();
    
    CLJGroup& operator=(const CLJGroup &other);
    
    bool operator==(const CLJGroup &other) const;
    bool operator!=(const CLJGroup &other) const;
    
    static const char* typeName();
    
    const char* what() const;
    
    QString toString() const;
    
    Length boxLength() const;
    void setBoxLength(Length box_length);
    
    bool isEmpty() const;
    
    Molecules molecules() const;
    
    PropertyMap mapForMolecule(MolNum molnum) const;
    
    void add(const MoleculeView &molview, const PropertyMap &map = PropertyMap());
    void add(const Molecules &molecules, const PropertyMap &map = PropertyMap());
    void add(const MoleculeGroup &molgroup, const PropertyMap &map = PropertyMap());
    
    void update(const MoleculeView &molview);
    void update(const Molecules &molecules);
    void update(const MoleculeGroup &molecules);
    
    void updatedConnectedGroup();
    
    void remove(const MoleculeView &molview);
    void remove(const Molecules &molecules);
    void remove(const MoleculeGroup &molecules);
    
    void remove(MolNum molnum);
    
    void removeAll();
    
    bool needsAccepting() const;
    void accept();
    
    const CLJBoxes& cljBoxes() const;
    
    CLJAtoms changedAtoms() const;
    CLJAtoms newAtoms() const;
    CLJAtoms oldAtoms() const;
    
    bool isSingleIDChange() const;
    tuple<CLJAtoms,CLJAtoms,CLJAtoms> mergeChanges() const;
    
    int nChangedMolecules() const;
    
    Molecules changedMolecules() const;
    
    void mustRecalculateFromScratch();
    
    void mustReallyRecalculateFromScratch();
    
    bool recalculatingFromScratch() const;
    
private:
    /** All of the extractors that manage extracting the charge and LJ 
        properties from all of the molecules */
    SireBase::ChunkedHash<MolNum,CLJExtractor> cljexts;
    
    /** The boxes to which the CLJAtoms from each molecule in this group
        have been added */
    CLJBoxes cljboxes;
    
    /** The workspace used to cache information about moves */
    CLJWorkspace cljworkspace;
    
    /** All of the changed molecules */
    QHash<MolNum,CLJExtractor> changed_mols;
    
    /** The set of all of the property maps used for each molecule */
    QHash<MolNum,PropertyMap> props;
    
    /** The source for the ID numbers for the atoms */
    CLJAtoms::ID_SOURCE id_source;
    
    /** How to extract atoms in the CLJExtractor */
    CLJExtractor::EXTRACT_SOURCE extract_source;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Return the CLJBoxes that contains all of the atoms in this group */
inline const CLJBoxes& CLJGroup::cljBoxes() const
{
    return cljboxes;
}

/** Return whether or not we are recalculating the energy from scratch */
inline bool CLJGroup::recalculatingFromScratch() const
{
    return cljworkspace.recalculatingFromScratch();
}

/** Return whether or not we don't have any molecules */
inline bool CLJGroup::isEmpty() const
{
    return cljexts.isEmpty() and changed_mols.isEmpty();
}

#endif // SIRE_SKIP_INLINE_FUNCTIONS

}

Q_DECLARE_METATYPE( SireMM::CLJGroup )

SIRE_EXPOSE_CLASS( SireMM::CLJGroup )

SIRE_END_HEADER

#endif
