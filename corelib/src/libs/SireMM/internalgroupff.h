/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2022  Christopher Woods
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

#ifndef SIREMM_INTERNALGROUPFF_H
#define SIREMM_INTERNALGROUPFF_H

#include "SireBase/propertymap.h"

#include "SireFF/g2ff.h"
#include "SireFF/ff3d.h"
#include "SireMM/internalff.h"

namespace SireMM
{
class InternalGroupFF;
}

SIREMM_EXPORT QDataStream &operator<<(QDataStream&, const SireMM::InternalGroupFF&);
SIREMM_EXPORT QDataStream &operator>>(QDataStream&, SireMM::InternalGroupFF&);

namespace SireMM
{

using SireFF::G2FF;

/** This is a forcefield that calculates the energies and forces
    caused by molecular mechanics style internal intramolecular
    potentials, e.g. bond, angle, dihedral, urey bradley terms
    between two groups (i.e. to calculate bond energy between
    two groups within a molecule)

    @author Christopher Woods
*/
class SIREMM_EXPORT InternalGroupFF
            : public SireBase::ConcreteProperty<InternalGroupFF,G2FF>,
              public FF3D,
              protected InternalPotential
{

friend SIREMM_EXPORT QDataStream& ::operator<<(QDataStream&, const InternalGroupFF&);
friend SIREMM_EXPORT QDataStream& ::operator>>(QDataStream&, InternalGroupFF&);

public:
    typedef InternalPotential::Components Components;
    typedef InternalPotential::ParameterNames ParameterNames;

    InternalGroupFF();
    InternalGroupFF(const QString &name);

    InternalGroupFF(const InternalGroupFF &other);

    ~InternalGroupFF();

    static const char* typeName();

    InternalGroupFF* clone() const;

    InternalGroupFF& operator=(const InternalGroupFF &other);

    bool operator==(const InternalGroupFF &other) const;
    bool operator!=(const InternalGroupFF &other) const;

    const InternalSymbols& symbols() const;
    const ParameterNames& parameters() const;
    const Components& components() const;

    void setArithmeticCombiningRules(bool on);
    void setGeometricCombiningRules(bool on);

    CLJFunction::COMBINING_RULES combiningRules() const;
    bool setCombiningRules(CLJFunction::COMBINING_RULES rules);

    bool usingArithmeticCombiningRules() const;
    bool usingGeometricCombiningRules() const;

    void enable14Calculation();
    void disable14Calculation();
    bool setUse14Calculation(bool on);
    bool uses14Calculation() const;

    bool setProperty(const QString &name, const Property &property);
    const Property& property(const QString &name) const;
    bool containsProperty(const QString &name) const;
    const Properties& properties() const;

    SireUnits::Dimension::MolarEnergy energy();
    SireUnits::Dimension::MolarEnergy energy(const Symbol &component);

    void energy(EnergyTable &energytable, double scale_energy=1);

    void energy(EnergyTable &energytable, const Symbol &symbol,
		double scale_energy=1);

    void force(ForceTable &forcetable, double scale_force=1);

    void force(ForceTable &forcetable, const Symbol &symbol,
               double scale_force=1);

    void field(FieldTable &fieldtable, double scale_field=1);

    void field(FieldTable &fieldtable, const Symbol &component,
               double scale_field=1);

    void potential(PotentialTable &potentialtable, double scale_potential=1);

    void potential(PotentialTable &potentialtable, const Symbol &component,
                   double scale_potential=1);

    void field(FieldTable &fieldtable, const Probe &probe, double scale_field=1);

    void field(FieldTable &fieldtable, const Symbol &component,
               const Probe &probe, double scale_field=1);

    void potential(PotentialTable &potentialtable, const Probe &probe,
                   double scale_potential=1);

    void potential(PotentialTable &potentialtable, const Symbol &component,
                   const Probe &probe, double scale_potential=1);

    void mustNowRecalculateFromScratch();

protected:
    const FFComponent& _pvt_components() const;

    void recalculateEnergy();

    void _pvt_added(quint32 group_id,
                    const PartialMolecule &mol, const PropertyMap &map);

    void _pvt_removed(quint32 group_id, const PartialMolecule &mol);

    void _pvt_changed(quint32 group_id,
                      const SireMol::Molecule &molecule, bool auto_update);
    void _pvt_changed(quint32 group_id,
                      const QList<SireMol::Molecule> &molecules, bool auto_update);

    void _pvt_removedAll(quint32 group_id);

    bool _pvt_wouldChangeProperties(quint32 group_id, MolNum molnum,
                                    const PropertyMap &map) const;

    void _pvt_extractMoleculeData(MolNum molnum,
                                  const PropertyMap &map);

    void _pvt_updateName();

private:
    typedef InternalPotential::Molecule Molecule;
    typedef InternalPotential::Molecules Molecules;

    /** All of the molecules currently in this forcefield */
    Molecules mols;

    /** All of the molecules currently in group 0 */
    Molecules mols0;

    /** All of the molecules currently in group 1 */
    Molecules mols1;

    /** The CLJ14Group that is used to calculate the 1-4 nonbonded energy
        of all contained molecules */
    QHash<MolNum,CLJ14Group> cljgroups;

    /** And then the same for group 0 */
    QHash<MolNum,CLJ14Group> cljgroups0;

    /** And then the same for group 1 */
    QHash<MolNum,CLJ14Group> cljgroups1;

    /** All of the (non-default) property maps for the molecules */
    QHash<MolNum,PropertyMap> propmaps;

    /** The energy components available for this forcefield */
    Components ffcomponents;

    /** The properties of this forcefield */
    Properties props;

    /** Whether or not to calculate 1-4 nonbonded energies */
    bool calc_14_nrgs;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

////////
//////// Inline functions of InternalFF
////////

SIRE_ALWAYS_INLINE SireUnits::Dimension::MolarEnergy InternalGroupFF::energy()
{
    return G2FF::energy();
}

SIRE_ALWAYS_INLINE SireUnits::Dimension::MolarEnergy InternalGroupFF::energy(const Symbol &component)
{
    return G2FF::energy(component);
}

/** Return all of the symbols used in the internal energy functions */
SIRE_ALWAYS_INLINE const InternalSymbols& InternalGroupFF::symbols() const
{
    return InternalPotential::symbols();
}

/** Return the names of all of the properties used to store the
    parameters for this potential */
SIRE_ALWAYS_INLINE const InternalFF::ParameterNames& InternalGroupFF::parameters() const
{
    return InternalPotential::parameters();
}

/** Return all of the symbols representing the components
    of this forcefield */
SIRE_ALWAYS_INLINE const InternalFF::Components& InternalGroupFF::components() const
{
    return ffcomponents;
}

#endif // SIRE_SKIP_INLINE_FUNCTIONS

}

Q_DECLARE_METATYPE( SireMM::InternalGroupFF );

SIRE_EXPOSE_CLASS( SireMM::InternalGroupFF )

#endif
