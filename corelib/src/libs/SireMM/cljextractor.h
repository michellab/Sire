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

#ifndef SIREMM_CLJEXTRACTOR_H
#define SIREMM_CLJEXTRACTOR_H

#include "cljatoms.h"
#include "cljdelta.h"
#include "cljworkspace.h"

#include "SireMol/molecule.h"
#include "SireMol/partialmolecule.h"
#include "SireMol/atomselection.h"

#include "SireBase/propertymap.h"

SIRE_BEGIN_HEADER

namespace SireMM
{
class CLJExtractor;
}

QDataStream& operator<<(QDataStream&, const SireMM::CLJExtractor&);
QDataStream& operator>>(QDataStream&, SireMM::CLJExtractor&);

namespace SireMM
{

using SireMol::MoleculeView;
using SireMol::Molecule;
using SireMol::AtomSelection;
using SireMol::PartialMolecule;

using SireBase::PropertyMap;
using SireBase::PropertyName;

/** This class is used to convert from a MoleculeView into a set of CLJAtoms
    objects. This class manages the extraction of data, recording of molecular
    properties and holding of all metadata needed to minimise the work of
    extracting data, e.g. as the molecule is updated or changed during
    a Monte Carlo simulation 
    
    @author Christopher Woods
*/
class SIREMM_EXPORT CLJExtractor
{

friend QDataStream& ::operator<<(QDataStream&, const CLJExtractor&);
friend QDataStream& ::operator>>(QDataStream&, CLJExtractor&);

public:
    enum EXTRACT_SOURCE
    {
        EXTRACT_BY_CUTGROUP = 0,
        EXTRACT_BY_RESIDUE = 1,
        EXTRACT_BY_MOLECULE = 2
    };

    CLJExtractor();
    CLJExtractor(const MoleculeView &mol, const PropertyMap &map = PropertyMap());
    CLJExtractor(const MoleculeView &mol, EXTRACT_SOURCE extract_source,
                 const PropertyMap &map = PropertyMap());

    CLJExtractor(const MoleculeView &mol, CLJAtoms::ID_SOURCE id_source,
                 const PropertyMap &map = PropertyMap());
    CLJExtractor(const MoleculeView &mol, CLJAtoms::ID_SOURCE id_source,
                 EXTRACT_SOURCE extract_source,
                 const PropertyMap &map = PropertyMap());
    
    CLJExtractor(const CLJExtractor &other);
    
    ~CLJExtractor();
    
    CLJExtractor& operator=(const CLJExtractor &other);
    
    bool operator==(const CLJExtractor &other) const;
    bool operator!=(const CLJExtractor &other) const;
    
    static const char* typeName();
    
    const char* what() const;
    
    QString toString() const;
    
    bool changed() const;
    bool needsCommitting() const;
    
    bool isEmpty() const;
    bool isNull() const;
    
    PartialMolecule oldMolecule() const;
    PartialMolecule newMolecule() const;
    
    PropertyMap propertyMap() const;
    
    PropertyName coordinatesProperty() const;
    PropertyName chargeProperty() const;
    PropertyName ljProperty() const;
    
    bool extractingByCutGroup() const;
    bool extractingByResidue() const;
    bool extractingByMolecule() const;
    
    CLJAtoms::ID_SOURCE idSource() const;
    
    void add(const MoleculeView &new_molecule, CLJBoxes &boxes, CLJWorkspace &workspace);
    void add(const AtomSelection &new_selection, CLJBoxes &boxes, CLJWorkspace &workspace);
    
    void updateSelection(const AtomSelection &selection,
                         CLJBoxes &boxes, CLJWorkspace &workspace);
    
    void update(const MoleculeView &new_molecule,
                CLJBoxes &boxes, CLJWorkspace &workspace);

    void remove(const AtomSelection &new_selection,
                CLJBoxes &boxes, CLJWorkspace &workspace);
    void remove(const MoleculeView &new_molecule,
                CLJBoxes &boxes, CLJWorkspace &workspace);
                
    void removeAll(CLJBoxes &boxes, CLJWorkspace &workspace);

    void commit(CLJBoxes &boxes, CLJWorkspace &workspace);
    void revert(CLJBoxes &boxes, CLJWorkspace &workspace);

private:
    void initialise(CLJBoxes &boxes, CLJWorkspace &workspace);

    /** Copy of the molecule itself */
    Molecule mol;
    
    /** Copy of the current atom selection (empty if we have
        selected the entire molecule) */
    AtomSelection selected_atoms;
    
    /** Copy of the new molecule  */
    Molecule newmol;
    
    /** Copy of the new selection */
    AtomSelection new_selected_atoms;
    
    /** The property map used to extract data (empty if we are
        using default properties) */
    PropertyMap props;
    
    /** The indicies of all of the CLJAtoms in the CLJBoxes */
    QVector< QVector<CLJBoxIndex> > cljidxs;

    /** A copy of all of the deltas from updating or adding to the molecule */
    QVector<CLJDelta> cljdeltas;

    /** The source of the ID property in CLJAtoms */
    CLJAtoms::ID_SOURCE id_source;

    /** How are we extracting atoms? */
    EXTRACT_SOURCE extract_source;
};

}

Q_DECLARE_METATYPE( SireMM::CLJExtractor )

SIRE_EXPOSE_CLASS( SireMM::CLJExtractor )

SIRE_END_HEADER

#endif
