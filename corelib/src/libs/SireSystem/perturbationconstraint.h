/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2009  Christopher Woods
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

#ifndef SIRESYSTEM_PERTURBATIONCONSTRAINT_H
#define SIRESYSTEM_PERTURBATIONCONSTRAINT_H

#include "moleculeconstraint.h"

#include "SireBase/shareddatapointer.hpp"
#include "SireBase/refcountdata.h"

#include "SireMol/perturbation.h"
#include "SireMol/moleculegroup.h"

SIRE_BEGIN_HEADER

namespace SireSystem
{
class PerturbationConstraint;
}

SIRESYSTEM_EXPORT QDataStream& operator<<(QDataStream&, const SireSystem::PerturbationConstraint&);
SIRESYSTEM_EXPORT QDataStream& operator>>(QDataStream&, SireSystem::PerturbationConstraint&);

namespace SireSystem
{

using SireMol::Molecule;
using SireMol::MoleculeGroup;
using SireMol::MolGroupPtr;

using SireMol::Perturbation;
using SireMol::PerturbationPtr;

using SireCAS::Values;

using SireBase::PropertyMap;

namespace detail
{

class PerturbationData : public SireBase::RefCountData
{
public:
    PerturbationData();
    PerturbationData(const PerturbationPtr &perturbation);
    
    PerturbationData(const PerturbationData &other);
    
    ~PerturbationData();
    
    bool wouldChange(const Molecule &molecule, const Values &values) const;

    Molecule perturb(const Molecule &molecule, const Values &values);
    
private:
    /** The actual perturbation */
    PerturbationPtr pert;
    
    /** The properties required, and the version of 
        the property in the molecule the last time it
        was updated */
    QHash<QString,quint64> props;
};

}

/** This constraint is used to constrain part or parts of a z-matrix
    to particular values - this is useful during single-topology
    free energy simulations
    
    @author Christopher Woods
*/
class SIRESYSTEM_EXPORT PerturbationConstraint
         : public SireBase::ConcreteProperty<PerturbationConstraint,MoleculeConstraint>
{

friend SIRESYSTEM_EXPORT QDataStream& ::operator<<(QDataStream&, const PerturbationConstraint&);
friend SIRESYSTEM_EXPORT QDataStream& ::operator>>(QDataStream&, PerturbationConstraint&);

public:
    PerturbationConstraint();
    PerturbationConstraint(const MoleculeGroup &molgroup, 
                           const PropertyMap &map = PropertyMap());
                   
    PerturbationConstraint(const PerturbationConstraint &other);
    
    ~PerturbationConstraint();
    
    PerturbationConstraint& operator=(const PerturbationConstraint &other);
    
    bool operator==(const PerturbationConstraint &other) const;
    bool operator!=(const PerturbationConstraint &other) const;
    
    static const char* typeName();
    
    QString toString() const;
    
    const MoleculeGroup& moleculeGroup() const;

    SireBase::PropertyName perturbationProperty() const;

protected:
    void setSystem(const System &system);
    bool mayChange(const Delta &delta, quint32 last_subversion) const;

    bool fullApply(Delta &delta);
    bool deltaApply(Delta &delta, quint32 last_subversion);

private:
    bool pvt_update(Molecule &molecule, const SireCAS::Values &values);

    /** The molecule group containing the molecules that are affected
        by these perturbations */
    MolGroupPtr molgroup;

    /** The name of the property containing the perturbations
        to be applied to each molecule in the group */
    SireBase::PropertyName perts_property;

    typedef QList< SireBase::SharedDataPointer<detail::PerturbationData> > PertDataList;
    typedef QHash<MolNum,PertDataList> PertDataHash;

    /** Information about all of the perturbations about each molecule */
    PertDataHash pertdata;

    /** All of the symbols of components used by all of the perturbations */
    QSet<SireCAS::Symbol> all_pert_syms;

    /** The symbols used by the perturbations for each molecule, 
        and their current values */
    QHash< MolNum,QSet<SireCAS::Symbol> > pert_syms;

    /** The molecules that must change to maintain the constraint */
    Molecules changed_mols;
};

}

Q_DECLARE_METATYPE( SireSystem::PerturbationConstraint )

SIRE_EXPOSE_CLASS( SireSystem::PerturbationConstraint )

SIRE_END_HEADER

#endif
