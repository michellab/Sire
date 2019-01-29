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

#ifndef SIREMM_CLJ14GROUP_H
#define SIREMM_CLJ14GROUP_H

#include "cljfunction.h"

#include "SireMol/partialmolecule.h"

#include <boost/tuple/tuple.hpp>

SIRE_BEGIN_HEADER

namespace SireMM
{
class CLJ14Group;
}

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::CLJ14Group&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::CLJ14Group&);

namespace SireMM
{

namespace detail
{
class CLJ14PairData;
}

using SireMol::PartialMolecule;
using SireMol::MoleculeView;
using SireMol::CGIdx;
using SireMol::AtomSelection;

using SireBase::PropertyMap;

/** This class holds all of the information needed to calculate
    the 14-nonbonded energy for a molecule
    
    @author Christopher Woods
*/
class SIREMM_EXPORT CLJ14Group
{

friend QDataStream& ::operator<<(QDataStream&, const CLJ14Group&);
friend QDataStream& ::operator>>(QDataStream&, CLJ14Group&);

public:
    CLJ14Group();
    CLJ14Group(const MoleculeView &molecule, const PropertyMap &map = PropertyMap());
    CLJ14Group(const MoleculeView &molecule, CLJFunction::COMBINING_RULES combining_rules,
               bool is_strict, const PropertyMap &map = PropertyMap());
    CLJ14Group(const CLJ14Group &other);
    
    ~CLJ14Group();
    
    CLJ14Group& operator=(const CLJ14Group &other);
    
    bool operator==(const CLJ14Group &other) const;
    bool operator!=(const CLJ14Group &other) const;
    
    static const char* typeName();
    
    const char* what() const;
    
    QString toString() const;
    
    bool isNull() const;

    bool setStrict(bool isstrict);
    bool isStrict() const;
    
    const MoleculeView& molecule() const;
    PropertyMap propertyMap() const;

    void add(const MoleculeView &new_molecule);
    void add(const AtomSelection &new_selection);
    
    void updateSelection(const AtomSelection &selection);
    
    void update(const MoleculeView &new_molecule);

    void remove(const AtomSelection &new_selection);
    void remove(const MoleculeView &new_molecule);
    
    boost::tuple<double,double> energy();

    bool recalculatingFromScratch() const;
    void mustNowRecalculateFromScratch();
    void mustReallyRecalculateFromScratch();
    
    void setArithmeticCombiningRules(bool on);
    void setGeometricCombiningRules(bool on);
    
    CLJFunction::COMBINING_RULES combiningRules() const;
    void setCombiningRules(CLJFunction::COMBINING_RULES rules);
    
    bool usingArithmeticCombiningRules() const;
    bool usingGeometricCombiningRules() const;
    
    bool wouldChangeProperties(const PropertyMap &map) const;
    
private:
    void addCGData(CGIdx cg0, CGIdx cg1, const QVector<detail::CLJ14PairData> &cgdata);
    void reextract();

    /** The molecule whose 14 energy is being calculated */
    PartialMolecule mol;
    
    /** The new molecule if we have been updated */
    PartialMolecule newmol;
    
    /** The property map to use to extract parameters */
    PropertyMap propmap;
    
    /** The 14 pair data for all CutGroups, and between all 
        bonded pairs of CutGroups */
    QVector< QVector<detail::CLJ14PairData> > data_for_pair;

    /** The indicies of all elements in 'cgroups' that contains 14 data
        about the CutGroup at CGIdx 'cgidx' */
    QHash< quint32, QSet<quint32> > cgidx_to_idx;
    
    /** Combining rules to use for the LJ calculation */
    CLJFunction::COMBINING_RULES combining_rules;
    
    /** The current coulomb and LJ energies */
    double total_cnrg, total_ljnrg;
    
    /** Whether or not the energy needs to be recalculated */
    bool needs_energy;
    
    /** Whether or not this group is strict
        (1-4 values are only calculated when both atoms in the 
         pair are contained in the group, as opposed to when only
         one of the atoms are in the group) */
    bool is_strict;
};

}

Q_DECLARE_METATYPE( SireMM::CLJ14Group )

SIRE_EXPOSE_CLASS( SireMM::CLJ14Group )

SIRE_END_HEADER

#endif
