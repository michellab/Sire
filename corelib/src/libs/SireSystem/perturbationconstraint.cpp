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

#include "perturbationconstraint.h"

#include "SireCAS/symbol.h"
#include "SireCAS/values.h"

#include "SireMol/moleculegroup.h"
#include "SireMol/molecule.h"
#include "SireMol/perturbation.h"
#include "SireMol/mgname.h"
#include "SireMol/mgnum.h"

#include "SireSystem/system.h"
#include "SireSystem/delta.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireSystem;
using namespace SireSystem::detail;
using namespace SireMol;
using namespace SireCAS;
using namespace SireBase;
using namespace SireStream;

/////////
///////// Implementation of PerturbationData
/////////

PerturbationData::PerturbationData()
{}

PerturbationData::PerturbationData(const PerturbationPtr &perturbation)
                 : pert(perturbation)
{}
    
PerturbationData::PerturbationData(const PerturbationData &other)
                 : pert(other.pert), props(other.props)
{}
    
PerturbationData::~PerturbationData()
{}
    
bool PerturbationData::wouldChange(const Molecule &molecule,
                                   const Values &values) const
{
    if (props.isEmpty())
        return true;
    else
    {
        for (QHash<QString,quint64>::const_iterator it = props.constBegin();
             it != props.constEnd();
             ++it)
        {
            if (molecule.version(it.key()) != it.value())
            {
                //one of the properties needed by this perturbation has
                //changed - see if this will change the molecule
                return pert.read().wouldChange(molecule, values);
            }
        }
        
        return false;
    }
}

Molecule PerturbationData::perturb(const Molecule &molecule, const Values &values)
{
    if (props.isEmpty())
    {
        if (not pert.read().wouldChange(molecule, values))
        {
            //now save the versions of the properties used by this perturbation
            foreach (QString property, pert.read().requiredProperties())
            {
                props.insert(property, molecule.version(property));
            }   
            
            return molecule;
        }
    }

    //apply the perturbation
    Molecule perturbed_mol = pert.read().perturb(molecule, values);
    
    //now save the versions of the properties used by this perturbation
    foreach (QString property, pert.read().requiredProperties())
    {
        props.insert(property, perturbed_mol.version(property));
    }
    
    return perturbed_mol;
}

/////////
///////// Implementation of PerturbationConstraint
/////////

static const RegisterMetaType<PerturbationConstraint> r_pertcons;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                          const PerturbationConstraint &pertcons)
{
    writeHeader(ds, r_pertcons, 1);
    
    SharedDataStream sds(ds);
    
    sds << pertcons.molgroup << pertcons.perts_property
        << static_cast<const MoleculeConstraint&>(pertcons);
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds,
                                          PerturbationConstraint &pertcons)
{
    VersionID v = readHeader(ds, r_pertcons);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);

        pertcons = PerturbationConstraint();
        
        sds >> pertcons.molgroup >> pertcons.perts_property
            >> static_cast<MoleculeConstraint&>(pertcons);
    }
    else
        throw version_error(v, "1", r_pertcons, CODELOC);
        
    return ds;
}

/** Constructor */
PerturbationConstraint::PerturbationConstraint()
                       : ConcreteProperty<PerturbationConstraint,MoleculeConstraint>(),
                         perts_property("perturbations")
{}

/** Construct specifying the molecule group and perturbation property */
PerturbationConstraint::PerturbationConstraint(const MoleculeGroup &mgroup, 
                                               const PropertyMap &map)
                       : ConcreteProperty<PerturbationConstraint,MoleculeConstraint>(),
                         molgroup(mgroup), 
                         perts_property( map["perturbations"] )
{}
            
/** Copy constructor */
PerturbationConstraint::PerturbationConstraint(const PerturbationConstraint &other)
           : ConcreteProperty<PerturbationConstraint,MoleculeConstraint>(other),
             molgroup(other.molgroup), perts_property(other.perts_property),
             pertdata(other.pertdata), all_pert_syms(other.all_pert_syms), 
             pert_syms(other.pert_syms), changed_mols(other.changed_mols)
{}

/** Destructor */
PerturbationConstraint::~PerturbationConstraint()
{}

/** Copy assignment operator */
PerturbationConstraint& 
PerturbationConstraint::operator=(const PerturbationConstraint &other)
{
    if (this != &other)
    {
        MoleculeConstraint::operator=(other);
        molgroup = other.molgroup;
        perts_property = other.perts_property;
        pertdata = other.pertdata;
        all_pert_syms = other.all_pert_syms;
        pert_syms = other.pert_syms;
        changed_mols = other.changed_mols;
    }
    
    return *this;
}

/** Comparison operator */
bool PerturbationConstraint::operator==(const PerturbationConstraint &other) const
{
    return molgroup == other.molgroup and perts_property == other.perts_property;
}

/** Comparison operator */
bool PerturbationConstraint::operator!=(const PerturbationConstraint &other) const
{
    return not PerturbationConstraint::operator==(other);
}

const char* PerturbationConstraint::typeName()
{
    return QMetaType::typeName( qMetaTypeId<PerturbationConstraint>() );
}

/** Return a string representation of this constraint */
QString PerturbationConstraint::toString() const
{
    if (molgroup.constData() == 0)
        return QObject::tr("PerturbationConstraint::null");
    else
        return QObject::tr("PerturbationConstraint( [%1:%2], %3 )")
                .arg(molgroup.read().name().value())
                .arg(molgroup.read().number().value())
                .arg(perts_property.toString());
}

/** Return the molecule group that is acted on by this constraint */
const MoleculeGroup& PerturbationConstraint::moleculeGroup() const
{
    return *molgroup;
}

/** Return the property used to find the perturbations to apply
    to the molecules in this constraint */
PropertyName PerturbationConstraint::perturbationProperty() const
{
    return perts_property;
}

/** This is called once for each new molecule in the group to update
    it and put the perturbed molecule (if it needs perturbing) into
    the perturbed_mols hash */
bool PerturbationConstraint::pvt_update(Molecule &molecule, const Values &values)
{
    //the molecule has changed, but does it still obey the constraints?
    {
        PertDataList perts = pertdata.value(molecule.number());
                
        bool would_change = false;
                
        for (PertDataList::const_iterator it = perts.constBegin();
             it != perts.constEnd();
             ++it)
        {
            if ( (*it)->wouldChange(molecule, values) )
            {
                would_change = true;
                break;
            }
        }
                
        if (not would_change)
            return false;
    }
        
    //ok - we now know that the molecule will need to be updated
    // - apply the perturbation
    {
        PertDataList &perts = pertdata[molecule.number()];

        quint64 old_version = molecule.version();
            
        for (PertDataList::iterator it = perts.begin();
             it != perts.end();
             ++it)
        {
            if ((*it)->wouldChange(molecule, values))
                molecule = (*it)->perturb(molecule,values);
        }
            
        return molecule.version() != old_version;
    }
}

void PerturbationConstraint::setSystem(const System &system)
{
    if (Constraint::wasLastSystem(system) and Constraint::wasLastSubVersion(system))
        return;
        
    changed_mols = Molecules();
        
    if (molgroup.isNull())
    {
        Constraint::setSatisfied(system, true);
        return;
    }
    
    if ( system.contains(molgroup.read().number()) )
    {
        molgroup = system[molgroup.read().number()];
    }
    else
    {
        molgroup.edit().update(system.molecules());
    }
        
    const Molecules &molecules = molgroup.read().molecules();
        
    pertdata.clear();
    pert_syms.clear();
    all_pert_syms.clear();

    pertdata.reserve(molecules.nMolecules());
    pert_syms.reserve(molecules.nMolecules());
        
    //now update all of the contained molecules
    for (Molecules::const_iterator it = molecules.constBegin();
         it != molecules.constEnd();
         ++it)
    {
        Molecule molecule = it.value().molecule();
        MolNum molnum = it.key();

        const Perturbation &perturbation = molecule.property(perts_property)
                                                   .asA<Perturbation>();

        QSet<Symbol> perturbation_symbols = perturbation.requiredSymbols();

        Values values = system.constants(perturbation_symbols);
        pert_syms.insert(molnum, perturbation_symbols);
        all_pert_syms += perturbation_symbols;
        
        //perturb the molecule
        Molecule perturbed_mol = molecule;
        
        if (perturbation.wouldChange(molecule, values))
            perturbed_mol = perturbation.perturb(molecule, values);
        
        //now save information about all of the perturbations
        //(so that they can be applied individually in future)
        QList<PerturbationPtr> perts = perturbation.children();
                                               
        PertDataList &pertlist = pertdata[molnum];
        
        for (QList<PerturbationPtr>::const_iterator it2 = perts.constBegin();
             it2 != perts.constEnd();
             ++it2)
        {
            SharedDataPointer<PerturbationData> d( new PerturbationData(*it2) );
            pertlist.append(d);
        }
        
        if (molecule.version() != perturbed_mol.version())
            changed_mols.add(perturbed_mol);
    }
    
    Constraint::setSatisfied(system, changed_mols.isEmpty());
}

bool PerturbationConstraint::mayChange(const Delta &delta, quint32 last_subversion) const
{
    if (molgroup.isNull())
        return false;
    else
        return delta.sinceChanged(all_pert_syms, last_subversion) or
               delta.sinceChanged(molgroup.read().molecules(), last_subversion);
}

bool PerturbationConstraint::fullApply(Delta &delta)
{
    this->setSystem(delta.deltaSystem());
    
    if (not changed_mols.isEmpty())
    {
        bool changed = delta.update(changed_mols);
        changed_mols = Molecules();
        return changed;
    }
    else
        return false;
}

bool PerturbationConstraint::deltaApply(Delta &delta, quint32 last_subversion)
{
    if (molgroup.isNull())
        return false;

    else if (not changed_mols.isEmpty())
        return this->fullApply(delta);

    const System &system = delta.deltaSystem();

    //get the list of molecules that have changed
    QList<MolNum> changed_molnums = delta.changedMoleculesSince(
                                                molgroup.read().molecules(),
                                                last_subversion);

    if (not changed_molnums.isEmpty())
    {
        //we need to update our copy of the molecule group
        if (system.contains(molgroup.read().number()))
            molgroup = system[molgroup.read().number()];
        else
            molgroup.edit().update(system.molecules());
    }

    //now add to this molecules that must change because their
    //dependent components have changed
    if (delta.hasComponentChangeSince(last_subversion))
    {
        for (QHash< MolNum,QSet<Symbol> >::const_iterator it = pert_syms.constBegin();
             it != pert_syms.constEnd();
             ++it)
        {
            if ( (not changed_molnums.contains(it.key())) and 
                 delta.sinceChanged(it.value(), last_subversion) )
            {
                changed_molnums.append(it.key());
            }
        }
    }

    changed_mols = Molecules();
    const Molecules &molecules = molgroup.read().molecules();

    //change all of the molecules that need changing
    foreach (MolNum changed_molnum, changed_molnums)
    {
        Molecule molecule = molecules[changed_molnum].molecule();
        
        if (this->pvt_update( molecule, 
                              system.constants(pert_syms.value(changed_molnum))) )
        {
            changed_mols.add(molecule);
        }
    }
    
    if (not changed_mols.isEmpty())
    {
        bool changed = delta.update(changed_mols);
        
        if (delta.deltaSystem().contains(molgroup.read().number()))
            molgroup = delta.deltaSystem()[molgroup.read().number()];
        else
            molgroup.edit().update(changed_mols);
        
        changed_mols = Molecules();
        return changed;
    }
    else
        return false;
}
