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

#include "qmchargeconstraint.h"

#include "SireMol/molecules.h"
#include "SireMol/partialmolecule.h"
#include "SireMol/molecule.h"
#include "SireMol/atomcharges.h"
#include "SireMol/mgname.h"
#include "SireMol/moleculegroup.h"
#include "SireMol/moleditor.h"

#include "SireSystem/system.h"
#include "SireSystem/delta.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace Squire;
using namespace SireSystem;
using namespace SireMol;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<QMChargeConstraint> r_qmchgconstraint;

/** Serialise to a binary datastream */
QDataStream SQUIRE_EXPORT &operator<<(QDataStream &ds,
                                      const QMChargeConstraint &qmchgconstraint)
{
    writeHeader(ds, r_qmchgconstraint, 1);
    
    SharedDataStream sds(ds);
    
    sds << qmchgconstraint.charge_calculator
        << static_cast<const ChargeConstraint&>(qmchgconstraint);
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream SQUIRE_EXPORT &operator>>(QDataStream &ds,
                                      QMChargeConstraint &qmchgconstraint)
{
    VersionID v = readHeader(ds, r_qmchgconstraint);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> qmchgconstraint.charge_calculator
            >> static_cast<ChargeConstraint&>(qmchgconstraint);
            
        qmchgconstraint.must_recalc_from_scratch = true;
    }
    else
        throw version_error(v, "1", r_qmchgconstraint, CODELOC);
        
    return ds;
}

/** Constructor */
QMChargeConstraint::QMChargeConstraint()
                   : ConcreteProperty<QMChargeConstraint,ChargeConstraint>(),
                     must_recalc_from_scratch(true)
{}

/** Construct to constrain the charges for the molecules in the 
    molecule group 'molgroup' using the optionally supplied property
    map to find the necessary properteis */
QMChargeConstraint::QMChargeConstraint(const MoleculeGroup &molgroup,
                                       const PropertyMap &map)
                   : ConcreteProperty<QMChargeConstraint,ChargeConstraint>(molgroup,
                                                                           map),
                     must_recalc_from_scratch(true)
{}

/** Construct to constrain the charges for the molecules in the     
    molecule group 'molgroup' to those calculated using the 
    QM charge calculator 'chargecalculator', using the optionally
    supplied property map to find the necessary properties */
QMChargeConstraint::QMChargeConstraint(const MoleculeGroup &molgroup,
                                       const QMChargeCalculator &chargecalculator,
                                       const PropertyMap &map)
                   : ConcreteProperty<QMChargeConstraint,ChargeConstraint>(molgroup,
                                                                           map),
                     charge_calculator(chargecalculator),
                     must_recalc_from_scratch(true)
{}
                                                                           
/** Copy constructor */
QMChargeConstraint::QMChargeConstraint(const QMChargeConstraint &other)
                   : ConcreteProperty<QMChargeConstraint,ChargeConstraint>(other),
                     charge_calculator(other.charge_calculator),
                     must_recalc_from_scratch(other.must_recalc_from_scratch)
{}

/** Destructor */
QMChargeConstraint::~QMChargeConstraint()
{}

const char* QMChargeConstraint::typeName()
{
    return QMetaType::typeName( qMetaTypeId<QMChargeConstraint>() );
}

/** Copy assignment operator */
QMChargeConstraint& QMChargeConstraint::operator=(const QMChargeConstraint &other)
{
    if (this != &other)
    {
        charge_calculator = other.charge_calculator;
        must_recalc_from_scratch = other.must_recalc_from_scratch;
        
        ChargeConstraint::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool QMChargeConstraint::operator==(const QMChargeConstraint &other) const
{
    return charge_calculator == other.charge_calculator and
           must_recalc_from_scratch == other.must_recalc_from_scratch and
           ChargeConstraint::operator==(other);
}

/** Comparison operator */
bool QMChargeConstraint::operator!=(const QMChargeConstraint &other) const
{
    return not QMChargeConstraint::operator==(other);
}

/** Return a string representation of this constraint */
QString QMChargeConstraint::toString() const
{
    return QObject::tr("QMChargeConstraint( moleculeGroup() : %1, "
                       "chargeCalculator() : %2 )")
                            .arg(this->moleculeGroup().name())
                            .arg(this->chargeCalculator().toString());
}

/** Return the calculator used to calculate atomic partial charges
    from an underlying QM calculation */
const QMChargeCalculator& QMChargeConstraint::chargeCalculator() const
{
    return charge_calculator.read();
}

/** Set the charge calculator used to calculate atomic partial charges
    from an underlying QM calculation */
void QMChargeConstraint::setChargeCalculator(const QMChargeCalculator &chargecalculator)
{
    charge_calculator = chargecalculator;
    must_recalc_from_scratch = true;
}

/** Internal function used to calculate the charges of the passed
    molecule. This returns the molecule with the new charges, or
    an empty molecule if the charges haven't changed */
Molecule QMChargeConstraint::_pvt_calculateCharges(const PartialMolecule &molecule) const
{
    AtomCharges new_chgs = charge_calculator.read().calculate(molecule,
                                                              this->propertyMap());
    
    const PropertyName &charge_property = this->propertyMap()["charge"];

    if (molecule.hasProperty(charge_property))
    {
        const Property &old_chgs = molecule.property(charge_property);
        
        if (old_chgs.isA<AtomCharges>())
        {  
            if (old_chgs.asA<AtomCharges>() == new_chgs)
                //the charges haven't changed
                return Molecule();
        }
    }
    
    return molecule.molecule().edit()
                              .setProperty(charge_property, new_chgs)
                              .commit();
}

void QMChargeConstraint::setSystem(const System &system)
{
    if ( (not must_recalc_from_scratch) and
         Constraint::wasLastSystem(system) and Constraint::wasLastSubVersion(system) )
    {
        return;
    }

    this->updateGroup(system);

    mols_to_change = Molecules();

    const Molecules &mols = this->moleculeGroup().molecules();
        
    for (Molecules::const_iterator it = mols.constBegin();
         it != mols.constEnd();
         ++it)
    {
        Molecule new_mol = this->_pvt_calculateCharges(*it);
           
        if (not new_mol.isEmpty())
            mols_to_change.add(new_mol);
    }
        
    must_recalc_from_scratch = false;

    Constraint::setSatisfied(system, mols_to_change.isEmpty());
}

bool QMChargeConstraint::mayChange(const Delta &delta, quint32 last_subversion) const
{
    return delta.sinceChanged( moleculeGroup().molecules(), last_subversion );
}

bool QMChargeConstraint::fullApply(Delta &delta)
{
    if (must_recalc_from_scratch)
        this->setSystem(delta.deltaSystem());
        
    if (mols_to_change.isEmpty())
        return false;
    
    else
    {
        bool changed = delta.update(mols_to_change);
        mols_to_change = Molecules();
        
        if (changed)
            this->updateGroup(delta.deltaSystem());
        
        return changed;
    }
}

bool QMChargeConstraint::deltaApply(Delta &delta, quint32 last_subversion)
{
    if ( must_recalc_from_scratch or (not mols_to_change.isEmpty()) )
        return this->fullApply(delta);
        
    QList<MolNum> changed_mols = delta.changedMoleculesSince(last_subversion);
    
    if (changed_mols.isEmpty())
        return false;
        
    this->updateGroup(delta.deltaSystem());
    
    const Molecules &mols = moleculeGroup().molecules();
    
    foreach (MolNum changed_mol, changed_mols)
    {
        mols_to_change.remove(changed_mol);

        Molecules::const_iterator it = mols.constFind(changed_mol);
        
        if (it != mols.constEnd())
        {
            Molecule new_mol = this->_pvt_calculateCharges(*it);

            if (not new_mol.isEmpty())
                mols_to_change.add(new_mol);
        }
    }
    
    if (mols_to_change.isEmpty())
        return false;

    else
    {
        bool changed = delta.update(mols_to_change);
        mols_to_change = Molecules();
        
        if (changed)
            this->updateGroup(delta.deltaSystem());
            
        return changed;
    }
}
