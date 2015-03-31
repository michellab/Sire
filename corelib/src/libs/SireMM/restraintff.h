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

#ifndef SIREMM_RESTRAINTFF_H
#define SIREMM_RESTRAINTFF_H

#include "restraint.h"
#include "restraintcomponent.h"

#include "SireFF/g1ff.h"
#include "SireFF/ff3d.h"

#include "SireVol/space.h"

#include "SireBase/properties.h"

SIRE_BEGIN_HEADER

namespace SireMM
{
class RestraintFF;
}

QDataStream& operator<<(QDataStream&, const SireMM::RestraintFF&);
QDataStream& operator>>(QDataStream&, SireMM::RestraintFF&);

namespace SireFF
{
class EnergyTable;
class ForceTable;
class FieldTable;
class PotentialTable;
class Probe;
}

namespace SireMM
{

using SireBase::Properties;

using SireMol::Molecule;
using SireMol::MoleculeData;
using SireMol::PartialMolecule;
using SireMol::ViewsOfMol;
using SireMol::Molecules;

using SireFF::EnergyTable;
using SireFF::ForceTable;
using SireFF::FieldTable;
using SireFF::PotentialTable;
using SireFF::Probe;

using SireVol::Space;

using SireCAS::Symbol;
using SireCAS::Symbols;
using SireCAS::Values;

/** This is a forcefield that holds and evaluates a collection of 
    Restraint3D restraints.
    
    @author Christopher Woods
*/
class SIREMM_EXPORT RestraintFF 
            : public SireBase::ConcreteProperty<RestraintFF,SireFF::G1FF>,
              public SireFF::FF3D
{

friend QDataStream& ::operator<<(QDataStream&, const RestraintFF&);
friend QDataStream& ::operator>>(QDataStream&, RestraintFF&);

public:
    RestraintFF();
    RestraintFF(const QString &name);
    
    RestraintFF(const RestraintFF &other);
    
    ~RestraintFF();
    
    static const char* typeName();
    
    RestraintFF& operator=(const RestraintFF &other);
    
    bool operator==(const RestraintFF &other) const;
    bool operator!=(const RestraintFF &other) const;
    
    RestraintFF* clone() const;
    
    const RestraintComponent& components() const;
    
    const Space& space() const;
    
    bool setSpace(const Space &space);
    
    bool setValue(const Symbol &symbol, double value);
    double getValue(const Symbol &symbol) const;
    bool hasValue(const Symbol &symbol) const;
    
    bool setProperty(const QString &name, const Property &property);
    const Property& property(const QString &name) const;
    bool containsProperty(const QString &name) const;
    const Properties& properties() const;
    
    Symbols symbols() const;
    Symbols userSymbols() const;
    Symbols builtinSymbols() const;
    
    Values userValues() const;
    
    RestraintFF differentiate(const Symbol &symbol) const;

    using G1FF::add;
    using G1FF::remove;
    using G1FF::contains;
    
    bool add(const Restraint3D &restraint);
    
    bool contains(const Restraint3D &restraint) const;
    
    bool remove(const Restraint3D &restraint);

    void removeRestraintAt(int i);
    
    QVector<Restraint3DPtr> restraints() const;
    
    const Restraint3D& restraintAt(int i) const;
    
    int nRestraints() const;

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

    void field(FieldTable &fieldtable, const Probe &probe, double scale_field=1);
    
    void field(FieldTable &fieldtable, const Symbol &component,
               const Probe &probe, double scale_field=1);
               
    void potential(PotentialTable &potentialtable, const Probe &probe,
                   double scale_potential=1);
    
    void potential(PotentialTable &potentialtable, const Symbol &component,
                   const Probe &probe, double scale_potential=1);

protected:

    ////
    //// Virtual functions from SireFF::FF
    ////

    const RestraintComponent& _pvt_components() const;
    
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
    void reindexRestraints();
    void rebuildProperties();

    void updateRestraints(const MoleculeData &moldata);
    void updateRestraints(const Molecules &molecules);

    /** The components of the energy */
    RestraintComponent ffcomponents;
    
    /** All of the restraints in this forcefield */
    QVector<Restraint3DPtr> restraints_by_idx;
    
    /** All of the old restraints - these are used when recalculating
        the energy */
    QHash<quint32,Restraint3DPtr> old_restraints_by_idx;
    
    /** Index of which restraints involve which molecules */
    QHash< MolNum, QList<quint32> > restraints_by_molnum;

    /** The space in which the restraints are evaluated */
    SireVol::SpacePtr spce;

    /** All of the values of the user symbols in the restraints */
    Values user_values;
    
    /** All of the built-in symbols of the restraints */
    Symbols builtin_symbols;
    
    /** All of the properties of this forcefield */
    Properties props;
    
    /** Whether or not to recalculate everything from scratch */
    bool recalc_from_scratch;
};


}

Q_DECLARE_METATYPE( SireMM::RestraintFF )

SIRE_EXPOSE_CLASS( SireMM::RestraintFF )

SIRE_END_HEADER

#endif
