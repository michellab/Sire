/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2008  Christopher Woods
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

#ifndef SQUIRE_QMFF_H
#define SQUIRE_QMFF_H

#include "qmpotential.h"

#include "SireFF/g1ff.h"
#include "SireFF/ff3d.h"

SIRE_BEGIN_HEADER

namespace Squire
{
class QMFF;
}

QDataStream& operator<<(QDataStream&, const Squire::QMFF&);
QDataStream& operator>>(QDataStream&, Squire::QMFF&);

namespace Squire
{

using SireBase::Property;
using SireBase::Properties;

using SireFF::EnergyTable;
using SireFF::ForceTable;
using SireFF::FieldTable;
using SireFF::PotentialTable;
using SireFF::Probe;

using SireCAS::Symbol;

/** This is a forcefield that uses an external Quantum Chemical program
    to calculate the quantum mechanics energy and / or force on the
    contained molecules.
    
    @author Christopher Woods
*/
class SQUIRE_EXPORT QMFF : public SireBase::ConcreteProperty<QMFF,SireFF::G1FF>,
                           public SireFF::FF3D,
                           protected QMPotential
{

friend QDataStream& ::operator<<(QDataStream&, const QMFF&);
friend QDataStream& ::operator>>(QDataStream&, QMFF&);

public:
    typedef QMPotential::Parameters Parameters;

    QMFF();
    QMFF(const QString &name);
    
    QMFF(const QMFF &other);
    
    ~QMFF();
    
    static const char* typeName();
    
    const char* what() const
    {
        return QMFF::typeName();
    }
    
    QMFF& operator=(const QMFF &other);
    
    bool operator==(const QMFF &other) const;
    bool operator!=(const QMFF &other) const;
    
    const Components& components() const;

    Parameters parameters() const
    {
        return Parameters();
    }

    const Space& space() const;
    const QMProgram& quantumProgram() const;
    SireUnits::Dimension::MolarEnergy zeroEnergy() const;
    
    bool setSpace(const Space &space);
    bool setQuantumProgram(const QMProgram &qmprog);
    bool setZeroEnergy(SireUnits::Dimension::MolarEnergy zero_energy);

    bool setProperty(const QString &name, const Property &property);
    const Property& property(const QString &name) const;
    bool containsProperty(const QString &name) const;
    const Properties& properties() const;

    void mustNowRecalculateFromScratch();    

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

    void field(FieldTable &fieldtable, const SireFF::Probe &probe, double scale_field=1);
    
    void field(FieldTable &fieldtable, const Symbol &component,
               const SireFF::Probe &probe, double scale_field=1);
               
    void potential(PotentialTable &potentialtable, const SireFF::Probe &probe,
                   double scale_potential=1);
    
    void potential(PotentialTable &potentialtable, const Symbol &component,
                   const SireFF::Probe &probe, double scale_potential=1);

    QString energyCommandFile() const;
    QString forceCommandFile(const ForceTable &forcetable) const;

    QString fieldCommandFile(const FieldTable &fieldtable) const;
    QString fieldCommandFile(const FieldTable &fieldtable, 
                             const SireFF::Probe &probe) const;

    QString potentialCommandFile(const PotentialTable &pottable) const;
    QString potentialCommandFile(const PotentialTable &pottable,
                                 const SireFF::Probe &probe) const;

protected:

    ////
    //// Virtual functions from SireFF::FF
    ////

    const Components& _pvt_components() const;
    
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

    ////
    //// Virtual functions of QMFF
    ////
    
    void changedPotential();

private:
    /** The components of the energy */
    Components ffcomponents;
    
    /** All of the molecules in this forcefield */
    QMPotential::Molecules qmmols;
};

} // end of namespace Squire

Q_DECLARE_METATYPE( Squire::QMFF )

SIRE_EXPOSE_CLASS( Squire::QMFF )

SIRE_END_HEADER

#endif
