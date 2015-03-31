/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2010  Christopher Woods
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

#ifndef SIRESYSTEM_POLARISECHARGES_H
#define SIRESYSTEM_POLARISECHARGES_H

#include "chargeconstraint.h"

#include "SireMM/cljprobe.h"

#include "SireFF/g1ff.h"

#include "SireCAS/symbol.h"

SIRE_BEGIN_HEADER

namespace SireSystem
{
class PolariseCharges;
class PolariseChargesFF;
}

QDataStream& operator<<(QDataStream&, const SireSystem::PolariseCharges&);
QDataStream& operator>>(QDataStream&, SireSystem::PolariseCharges&);

QDataStream& operator<<(QDataStream&, const SireSystem::PolariseChargesFF&);
QDataStream& operator>>(QDataStream&, SireSystem::PolariseChargesFF&);

namespace SireSystem
{

namespace detail
{
class PolariseChargesData;
}

using SireFF::SingleComponent;

using SireMol::ViewsOfMol;
using SireMol::PartialMolecule;

using SireBase::Properties;

/** This charge constraint adjusts the partial charges of contained
    molecules to give the impression that the molecule contains
    polarisable dipoles. This is based on the method developed
    by Reynolds et al. (see ...)
    
    @author Christopher Woods
*/
class SIRESYSTEM_EXPORT PolariseCharges
         : public SireBase::ConcreteProperty<PolariseCharges,ChargeConstraint>
{

friend QDataStream& ::operator<<(QDataStream&, const PolariseCharges&);
friend QDataStream& ::operator>>(QDataStream&, PolariseCharges&);

public:
    PolariseCharges();
    PolariseCharges(const MoleculeGroup &molgroup,
                    const PropertyMap &map = PropertyMap());
    PolariseCharges(const MoleculeGroup &molgroup,
                    const SireFF::Probe &probe,
                    const PropertyMap &map = PropertyMap());

    PolariseCharges(const MoleculeGroup &molgroup,
                    const SireCAS::Symbol &field_component,
                    const PropertyMap &map = PropertyMap());
    PolariseCharges(const MoleculeGroup &molgroup,
                    const SireCAS::Symbol &field_component,
                    const SireFF::Probe &probe,
                    const PropertyMap &map = PropertyMap());
    
    PolariseCharges(const PolariseCharges &other);
    
    ~PolariseCharges();
    
    PolariseCharges& operator=(const PolariseCharges &other);
    
    bool operator==(const PolariseCharges &other) const;
    bool operator!=(const PolariseCharges &other) const;
    
    static const char* typeName();
    
    PolariseCharges* clone() const;
    
    QString toString() const;

    void setConvergenceLimit(double limit);
    
    double convergenceLimit() const;

    const SireCAS::Symbol& fieldComponent() const;

    const SireMM::CoulombProbe& probe() const;

    PolariseChargesFF selfEnergyFF() const;

protected:
    void setSystem(const System &system);

    bool mayChange(const Delta &delta, quint32 last_subversion) const;
    
    bool fullApply(Delta &delta);
    bool deltaApply(Delta &delta, quint32 last_subversion);

private:
    void setProbe(const SireFF::Probe &probe);

    /** The forcefield component that is used to calculate 
        the potential on the atoms to be polarised */
    SireCAS::Symbol field_component;
    
    /** The probe used to calculate the potential on the atoms */
    SireMM::CoulombProbe field_probe;
    
    /** The information about each molecule that is needed
        to calculate the polarisability */
    QHash<SireMol::MolNum,
            QSharedDataPointer<detail::PolariseChargesData> > moldata;

    /** The collection of molecules that have been changed by this constraint */
    Molecules changed_mols;
    
    /** The convergence limit - charges are only updated if 
        they change by more than this limit */
    double convergence_limit;
};

/** This class implements the forcefield that is used to calculate
    the self-energy of polarising the charges. This is a companion
    forcefield to the PolariseCharges constraint and is not
    designed to be used on its own
    
    @author Christopher Woods
*/
class SIRESYSTEM_EXPORT PolariseChargesFF
            : public SireBase::ConcreteProperty<PolariseChargesFF,SireFF::G1FF>
{

friend QDataStream& ::operator<<(QDataStream&, const PolariseChargesFF&);
friend QDataStream& ::operator>>(QDataStream&, PolariseChargesFF&);

public:
    PolariseChargesFF();
    PolariseChargesFF(const PolariseCharges &constraint);
    PolariseChargesFF(const QString &name, const PolariseCharges &constraint);
    
    PolariseChargesFF(const PolariseChargesFF &other);
    
    ~PolariseChargesFF();
    
    static const char* typeName();
    
    PolariseChargesFF& operator=(const PolariseChargesFF &other);
    
    bool operator==(const PolariseChargesFF &other) const;
    bool operator!=(const PolariseChargesFF &other) const;
    
    PolariseChargesFF* clone() const;
    
    const SingleComponent& components() const;
    
    bool setProperty(const QString &name, const Property &property);
    const Property& property(const QString &name) const;
    bool containsProperty(const QString &name) const;
    const Properties& properties() const;

    using G1FF::add;
    using G1FF::remove;
    using G1FF::contains;

    void mustNowRecalculateFromScratch();    

protected:

    ////
    //// Virtual functions from SireFF::FF
    ////

    const SingleComponent& _pvt_components() const;
    
    void recalculateEnergy();
    
    void _pvt_updateName();
    
    ////
    //// Virtual functions from SireFF::G1FF
    ////

    void _pvt_added(const SireMol::PartialMolecule &mol, 
                    const SireBase::PropertyMap&);
                    
    void _pvt_removed(const SireMol::PartialMolecule &mol);
    
    void _pvt_changed(const SireMol::Molecule &mol, bool auto_update);
    
    void _pvt_changed(const QList<SireMol::Molecule> &mols, bool auto_update);
    
    void _pvt_removedAll();
    
    bool _pvt_wouldChangeProperties(SireMol::MolNum molnum, 
                                    const SireBase::PropertyMap &map) const;

    void _pvt_added(const ViewsOfMol &mol, const PropertyMap &map);
                            
    void _pvt_removed(const ViewsOfMol &mol);

    void _pvt_removedAll(const PartialMolecule &mol);
    void _pvt_removedAll(const ViewsOfMol &mol);

private:
    /** The symbol for the single component of the energy */
    SingleComponent ffcomponent;
    
    /** The location of the self polarisation energy property - the
        constraint adds this energy as a property of the polarised molecules
        so that this forcefield can then extract the property to calculate
        the total self-polarisation energy */
    PropertyName energy_property;
    
    /** The cache of the energy for each molecule in the forcefield */
    QHash<MolNum,SireUnits::Dimension::MolarEnergy> molnrg;
};

}

Q_DECLARE_METATYPE( SireSystem::PolariseCharges )
Q_DECLARE_METATYPE( SireSystem::PolariseChargesFF )

SIRE_EXPOSE_CLASS( SireSystem::PolariseCharges )
SIRE_EXPOSE_CLASS( SireSystem::PolariseChargesFF )

SIRE_END_HEADER

#endif
