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

#ifndef SQUIRE_QMCHARGECONSTRAINT_H
#define SQUIRE_QMCHARGECONSTRAINT_H

#include "SireSystem/chargeconstraint.h"

#include "qmchargecalculator.h"

SIRE_BEGIN_HEADER

namespace Squire
{
class QMChargeConstraint;
}

SQUIRE_EXPORT QDataStream& operator<<(QDataStream&, const Squire::QMChargeConstraint&);
SQUIRE_EXPORT QDataStream& operator>>(QDataStream&, Squire::QMChargeConstraint&);

namespace Squire
{

using SireSystem::ChargeConstraint;
using SireSystem::System;
using SireSystem::Delta;

using SireMol::MoleculeGroup;
using SireMol::Molecules;
using SireMol::Molecule;

using SireBase::PropertyMap;

/** This is a charge constraint that constrains the charges
    of molecules to equal those calculated from QM calculations
    (e.g. the charges can be constrained to equal those
    from AM1-BCC calculations)
    
    @author Christopher Woods
*/
class SQUIRE_EXPORT QMChargeConstraint
          : public SireBase::ConcreteProperty<QMChargeConstraint,ChargeConstraint>
{

friend QDataStream& ::operator<<(QDataStream&, const QMChargeConstraint&);
friend QDataStream& ::operator>>(QDataStream&, QMChargeConstraint&);

public:
    QMChargeConstraint();
    
    QMChargeConstraint(const MoleculeGroup &molgroup,
                       const PropertyMap &map = PropertyMap());
    
    QMChargeConstraint(const MoleculeGroup &molgroup,
                       const QMChargeCalculator &chargecalculator,
                       const PropertyMap &map = PropertyMap());
                       
    QMChargeConstraint(const QMChargeConstraint &other);
    
    ~QMChargeConstraint();
    
    static const char* typeName();
    
    QMChargeConstraint& operator=(const QMChargeConstraint &other);
    
    bool operator==(const QMChargeConstraint &other) const;
    bool operator!=(const QMChargeConstraint &other) const;
    
    QString toString() const;
    
    const QMChargeCalculator& chargeCalculator() const;
    
    void setChargeCalculator(const QMChargeCalculator &chargecalculator);

protected:
    void setSystem(const System &system);
    bool mayChange(const Delta &delta, quint32 last_subversion) const;

    bool fullApply(Delta &delta);
    bool deltaApply(Delta &delta, quint32 last_subversion);

private:
    Molecule _pvt_calculateCharges(const PartialMolecule &molecule) const;

    /** The charge calculator used to calculate the charges */
    QMChargeCalculatorPtr charge_calculator;

    /** The set of molecules that need to change to maintain this constraint */
    Molecules mols_to_change;
    
    /** If this flag is set, then all of the charges need
        to be recalculated */
    bool must_recalc_from_scratch;
};

}

Q_DECLARE_METATYPE( Squire::QMChargeConstraint )

SIRE_EXPOSE_CLASS( Squire::QMChargeConstraint )

SIRE_END_HEADER

#endif

