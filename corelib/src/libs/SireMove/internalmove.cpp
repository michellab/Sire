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

#include "internalmove.h"
#include "flexibility.h"
#include "ensemble.h"

#include "SireSystem/system.h"

#include "SireMol/partialmolecule.h"
#include "SireMol/molecule.h"
#include "SireMol/moleditor.h"
#include "SireMol/mover.hpp"
#include "SireMol/atomidx.h"
#include "SireMol/connectivity.h"
#include "SireMol/bondid.h"
#include "SireMol/angleid.h"
#include "SireMol/dihedralid.h"

#include "SireMaths/vectorproperty.h"

#include "SireUnits/dimensions.h"
#include "SireUnits/units.h"
#include "SireUnits/temperature.h"
#include "SireUnits/units.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>
#include <QTime>

using namespace SireMove;
using namespace SireMol;
using namespace SireSystem;
using namespace SireMaths;
using namespace SireUnits;
using namespace SireUnits::Dimension;
using namespace SireStream;

static const RegisterMetaType<InternalMove> r_internalmove;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const InternalMove &internalmove)
{
    writeHeader(ds, r_internalmove, 2);
    
    SharedDataStream sds(ds);
    
    sds << internalmove.smplr
        << internalmove.center_function
        << internalmove.flexibility_property
        << static_cast<const MonteCarlo&>(internalmove);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, InternalMove &internalmove)
{
    VersionID v = readHeader(ds, r_internalmove);

    if (v == 2)
    {
        SharedDataStream sds(ds);
    
        sds >> internalmove.smplr
            >> internalmove.center_function
            >> internalmove.flexibility_property
            >> static_cast<MonteCarlo&>(internalmove);
    }
    else if (v == 1)
    {
        SharedDataStream sds(ds);
    
        internalmove.center_function = GetCOMPoint();
    
        sds >> internalmove.smplr 
            >> internalmove.flexibility_property
            >> static_cast<MonteCarlo&>(internalmove);
    }
    else
        throw version_error(v, "1", r_internalmove, CODELOC);

    return ds;
}

/** Null constructor */
InternalMove::InternalMove(const PropertyMap &map) 
             : ConcreteProperty<InternalMove,MonteCarlo>(map)
{
    flexibility_property = map["flexibility"];
    MonteCarlo::setEnsemble( Ensemble::NVT(25*celsius) );
    center_function = GetCOMPoint();
}

/** Construct the mover move for the passed group of molecules */
InternalMove::InternalMove(const MoleculeGroup &molgroup, const PropertyMap &map)
             : ConcreteProperty<InternalMove,MonteCarlo>(),
               smplr( UniformSampler(molgroup) )
{
    flexibility_property = map["flexibility"];
    MonteCarlo::setEnsemble( Ensemble::NVT(25*celsius) );
    smplr.edit().setGenerator( this->generator() );
    center_function = GetCOMPoint();
}

/** Construct the mover move that samples molecules from the
    passed sampler */
InternalMove::InternalMove(const Sampler &sampler, const PropertyMap &map)
             : ConcreteProperty<InternalMove,MonteCarlo>(),
               smplr(sampler)
{
    flexibility_property = map["flexibility"];
    MonteCarlo::setEnsemble( Ensemble::NVT(25*celsius) );
    smplr.edit().setGenerator( this->generator() );
    center_function = GetCOMPoint();
}

/** Copy constructor */
InternalMove::InternalMove(const InternalMove &other)
             : ConcreteProperty<InternalMove,MonteCarlo>(other),
               smplr(other.smplr), center_function(other.center_function),
               flexibility_property(other.flexibility_property)
{}

/** Destructor */
InternalMove::~InternalMove()
{}

/** Copy assignment operator */
InternalMove& InternalMove::operator=(const InternalMove &other)
{
    if (this != &other)
    {
        MonteCarlo::operator=(other);
        smplr = other.smplr;
        center_function = other.center_function;
        flexibility_property = other.flexibility_property;
    }
    
    return *this;
}

/** Comparison operator */
bool InternalMove::operator==(const InternalMove &other) const
{
    return MonteCarlo::operator==(other) and smplr == other.smplr and
           center_function == other.center_function and
           flexibility_property == other.flexibility_property;
}

/** Comparison operator */
bool InternalMove::operator!=(const InternalMove &other) const
{
    return not InternalMove::operator==(other);
}

/** Return a string representation of this move */
QString InternalMove::toString() const
{
    return QObject::tr("InternalMove( nAccepted() = %1 nRejected() == %2 )")
                            .arg( this->nAccepted() )
                            .arg( this->nRejected() );
}

/** Set the sampler used to sample molecules for this move */
void InternalMove::setSampler(const Sampler &sampler)
{
    smplr = sampler;
    smplr.edit().setGenerator( this->generator() );
}

/** Set the sampler so that it draws molecules uniformly from the
    molecule group 'molgroup' */
void InternalMove::setSampler(const MoleculeGroup &molgroup)
{
    this->setSampler( UniformSampler(molgroup) );
}

/** Return the sampler used to sample molecules to move */
const Sampler& InternalMove::sampler() const
{
    return smplr;
}

/** Return the molecule group that is sampled for this move */
const MoleculeGroup& InternalMove::moleculeGroup() const
{
    return smplr->group();
}

/** Return the property used to find the flexibility of each molecule*/
const PropertyName& InternalMove::flexibilityProperty() const
{
    return flexibility_property;
}

/** Set the name of the property used to find the flexibility of each molecule */
void InternalMove::setFlexibilityProperty(const PropertyName &property)
{
    flexibility_property = property;
    Move::setProperty("flexibility", flexibility_property);
}

/** Set the function used to find the center of the molecule */
void InternalMove::setCenterOfMolecule(const GetPoint &cf)
{
    center_function = cf;
}

/** Return the function used to find the center of the molecule */
const GetPoint& InternalMove::centerOfMolecule() const
{
    return center_function.read();
}

/** Set the random number generator used to generate the random
    number used for this move */
void InternalMove::setGenerator(const RanGenerator &rangenerator)
{
    MonteCarlo::setGenerator(rangenerator);
    smplr.edit().setGenerator(this->generator());
}

/** Internal function used to set the ensemble based on the
    passed temperature */
void InternalMove::_pvt_setTemperature(const Temperature &temperature)
{
    MonteCarlo::setEnsemble( Ensemble::NVT(temperature) );
}

double clipDouble(double value)
{
    static const double mindbl = std::numeric_limits<double>::min() * 0.1;
    static const double maxdbl = std::numeric_limits<double>::max() * 0.1;

    return (std::min)(maxdbl, (std::max)(mindbl, value));
}

/** Actually perform 'nmoves' moves of the molecules in the 
    system 'system', optionally recording simulation statistics
    if 'record_stats' is true */
void InternalMove::move(System &system, int nmoves, bool record_stats)
{
    if (nmoves <= 0)
        return;

    InternalMove old_state(*this);
    System old_system_state(system);

    try
    {
        const PropertyMap &map = Move::propertyMap();
        const PropertyName &center_property = map["center"];

        for (int i=0; i<nmoves; ++i)
        {
            double old_nrg = system.energy( this->energyComponent() );
            System old_system(system);
            SamplerPtr old_sampler(smplr);
      
            double old_bias = 1;
            double new_bias = 1;
      
            //move one molecule
            //update the sampler with the latest version of the molecules
            smplr.edit().updateFrom(system);

            //this will randomly select one molecule
            tuple<PartialMolecule,double> mol_and_bias = smplr.read().sample();
	  
            const PartialMolecule &oldmol = mol_and_bias.get<0>();
            
            if (smplr.read().isBiased())
                old_bias = mol_and_bias.get<1>();
	
            Flexibility flex = oldmol.property(flexibility_property).asA<Flexibility>();

            // Select the dofs to move. 
            QList<BondID> moved_bonds;
            QList<AngleID> moved_angles;
            QList<DihedralID> moved_dihedrals;
    
            QList<BondID> flex_bonds = flex.flexibleBonds();
            QList<AngleID> flex_angs = flex.flexibleAngles();
            QList<DihedralID> flex_dihs = flex.flexibleDihedrals();

            int maxbondvar = flex.maximumBondVar();
            int maxanglevar = flex.maximumAngleVar();
            int maxdihedralvar = flex.maximumDihedralVar();

            int nbonds = flex_bonds.count();
            int nangles = flex_angs.count();
            int ndihedrals = flex_dihs.count();

            // Select bonds
            if ( nbonds == 0  || maxbondvar < 0  || maxbondvar >= nbonds )
            {
                moved_bonds = flex_bonds;
            }
            else
            {
                int movecount = 0;
                while (movecount < maxbondvar)
                {
                    int rand = this->generator().randInt(0, nbonds - 1);
                    const BondID &bond = flex_bonds.at( rand );
                    if ( not moved_bonds.contains(bond) )
                    {
                        moved_bonds.append(bond);
                        ++movecount;
                    }
                }
            }
            
            // Select angles
            if ( nangles == 0  || maxanglevar < 0  || maxanglevar >= nangles )
            {
                moved_angles = flex_angs;
            }
            else
            {
                int movecount = 0;
                while (movecount < maxanglevar)
                {
                    int rand = this->generator().randInt(0, nangles - 1);
                    const AngleID &angle = flex_angs.at( rand );
                    if ( not moved_angles.contains(angle) )
                    {
                        moved_angles.append(angle);
                        ++movecount;
                    }
                }
            }
            
            // Select dihedrals
            if ( ndihedrals == 0  || maxdihedralvar < 0  || maxdihedralvar >= ndihedrals )
            {
                moved_dihedrals = flex_dihs;
            }
            else
            {
                int movecount = 0;
                while (movecount < maxdihedralvar)
                {
                    int rand = this->generator().randInt(0, ndihedrals - 1);
                    const DihedralID &dihedral = flex_dihs.at( rand );
                    if ( not moved_dihedrals.contains(dihedral) )
                    {
                        moved_dihedrals.append(dihedral);
                        ++movecount;
                    }
                }
            }

            // Now save the current centre of the molecule, so that we can
            // ensure that this doesn't change as a result of this move
            // (we can't have an internal move cause rigid body translation
            //  of the molecule)
            const bool has_center_property = (oldmol.selectedAll() and
                                              oldmol.hasProperty(center_property));
    
            Vector old_center;
            
            if (has_center_property)
            {
                old_center = oldmol.property(center_property).asA<VectorProperty>();
            }
            else
            {
                old_center = center_function.read()(oldmol,map);
            }

            // Now actually move the selected dofs
            Mover<Molecule> mol_mover = oldmol.molecule().move();

            // move the bonds of this molecule
            Length bond_delta;
            
            foreach (const BondID &bond, moved_bonds)
            {
                //const Length bond_delta_value = flex.bond_deltas[bond];
                double bond_delta_value = clipDouble(flex.delta(bond));
                bond_delta = Length( this->generator().rand(-bond_delta_value,
                                                             bond_delta_value) );

                mol_mover.change(bond, bond_delta, map);
            }

            // and the angles
            Angle angle_delta;
            
            foreach (const AngleID &angle, moved_angles)
            {
                //const Angle angle_delta_value = flex.angle_deltas[angle];
                double angle_delta_value = clipDouble(flex.delta(angle));
                angle_delta = Angle( this->generator().rand(-angle_delta_value,
                                                             angle_delta_value) );

                mol_mover.change(angle, angle_delta, map);
            }
	  
            // and the torsions
            Angle dihedral_delta;
            
            foreach (const DihedralID &dihedral, moved_dihedrals)
            {
                double angle_delta_value = clipDouble(flex.delta(dihedral));
                dihedral_delta =  Angle( this->generator().rand(-angle_delta_value,
                                                                 angle_delta_value) );

                // 50% chance to either rotate around central bond or to just change that dihedral
                if (this->generator().randBool())
                {                                 
                    //move only this specific dihedral
                    mol_mover.change(dihedral, dihedral_delta, map);
                }
                else
                {
                    BondID centralbond;
                    centralbond = BondID(dihedral.atom1(), dihedral.atom2());
                    mol_mover.change(centralbond, dihedral_delta, map);
                }

            }

            Molecule newmol = mol_mover.commit();

            //evaluate the new center of geometry, and translate the molecule
            //to correct for any motion
            Vector new_center;
            
            if (has_center_property)
            {
                new_center = newmol.property(center_property).asA<VectorProperty>();
            }
            else
            {
                new_center = center_function.read()(newmol,map);
            }
            
            if (old_center != new_center)
            {
                newmol = newmol.move().translate( old_center-new_center, map ).commit();
            }

            //update the system with the new coordinates
            system.update(newmol);

            //get the new bias on this molecule
            smplr.edit().updateFrom(system);
	  
            if (smplr.read().isBiased())
                new_bias = smplr.read().probabilityOf( PartialMolecule(newmol,
                                                                       oldmol.selection()) );

            //calculate the energy of the system
            double new_nrg = system.energy( this->energyComponent() );
    
            //accept or reject the move based on the change of energy
            //and the biasing factors
            if (not this->test(new_nrg, old_nrg, new_bias, old_bias))
            {
                //the move has been rejected - reset the state
                smplr = old_sampler;
                system = old_system;
            }

            if (record_stats)
            {
                system.collectStats();
            }
        }
    }
    catch(...)
    {
        system = old_system_state;
        this->operator=(old_state);
        throw;
    }
}

const char* InternalMove::typeName()
{
    return QMetaType::typeName( qMetaTypeId<InternalMove>() );
}
