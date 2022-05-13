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

#include "internalmovesingle.h"
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
#include "SireMol/core.h"

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
using namespace SireUnits;
using namespace SireUnits::Dimension;
using namespace SireStream;

static const RegisterMetaType<InternalMoveSingle> r_internalmovesingle;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const InternalMoveSingle &internalmove)
{
    writeHeader(ds, r_internalmovesingle, 1);

    SharedDataStream sds(ds);

    sds << internalmove.smplr
        << internalmove.flexibility_property
	<< internalmove.synched_molgroup
        << static_cast<const MonteCarlo&>(internalmove);

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, InternalMoveSingle &internalmove)
{
    VersionID v = readHeader(ds, r_internalmovesingle);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> internalmove.smplr
            >> internalmove.flexibility_property
	    >> internalmove.synched_molgroup
            >> static_cast<MonteCarlo&>(internalmove);
    }
    else
        throw version_error(v, "1", r_internalmovesingle, CODELOC);

    return ds;
}

/** Null constructor */
InternalMoveSingle::InternalMoveSingle(const PropertyMap &map)
             : ConcreteProperty<InternalMoveSingle,MonteCarlo>(map)
{
    flexibility_property = map["flexibility"];
    MonteCarlo::setEnsemble( Ensemble::NVT(25*celsius) );
}

/** Construct the internal move for the passed group of molecules */
InternalMoveSingle::InternalMoveSingle(const MoleculeGroup &molgroup, const PropertyMap &map)
             : ConcreteProperty<InternalMoveSingle,MonteCarlo>(),
               smplr( UniformSampler(molgroup) )
{
    flexibility_property = map["flexibility"];
    MonteCarlo::setEnsemble( Ensemble::NVT(25*celsius) );
    smplr.edit().setGenerator( this->generator() );
}

/** Construct the mover move that samples molecules from the
    passed sampler */
InternalMoveSingle::InternalMoveSingle(const Sampler &sampler, const PropertyMap &map)
             : ConcreteProperty<InternalMoveSingle,MonteCarlo>(),
               smplr(sampler)
{
    flexibility_property = map["flexibility"];
    MonteCarlo::setEnsemble( Ensemble::NVT(25*celsius) );
    smplr.edit().setGenerator( this->generator() );
}

/** Copy constructor */
InternalMoveSingle::InternalMoveSingle(const InternalMoveSingle &other)
             : ConcreteProperty<InternalMoveSingle,MonteCarlo>(other),
               smplr(other.smplr),
               flexibility_property(other.flexibility_property),
	       synched_molgroup(other.synched_molgroup)
{}

/** Destructor */
InternalMoveSingle::~InternalMoveSingle()
{}

/** Copy assignment operator */
InternalMoveSingle& InternalMoveSingle::operator=(const InternalMoveSingle &other)
{
    if (this != &other)
    {
        MonteCarlo::operator=(other);
        smplr = other.smplr;
        flexibility_property = other.flexibility_property;
	synched_molgroup = other.synched_molgroup;
    }

    return *this;
}

/** Comparison operator */
bool InternalMoveSingle::operator==(const InternalMoveSingle &other) const
{
    return MonteCarlo::operator==(other) and smplr == other.smplr and
           flexibility_property == other.flexibility_property and
           synched_molgroup == other.synched_molgroup ;
}

/** Comparison operator */
bool InternalMoveSingle::operator!=(const InternalMoveSingle &other) const
{
    return not InternalMoveSingle::operator==(other);
}

/** Return a string representation of this move */
QString InternalMoveSingle::toString() const
{
    return QObject::tr("InternalMoveSingle( nAccepted() = %1 nRejected() == %2 )")
                            .arg( this->nAccepted() )
                            .arg( this->nRejected() );
}

/** Set the sampler used to sample molecules for this move */
void InternalMoveSingle::setSampler(const Sampler &sampler)
{
    smplr = sampler;
    smplr.edit().setGenerator( this->generator() );
}

/** Set the sampler so that it draws molecules uniformly from the
    molecule group 'molgroup' */
void InternalMoveSingle::setSampler(const MoleculeGroup &molgroup)
{
    this->setSampler( UniformSampler(molgroup) );
}

/** Return the sampler used to sample molecules to move */
const Sampler& InternalMoveSingle::sampler() const
{
    return smplr;
}

/** Return the molecule group that is sampled for this move */
const MoleculeGroup& InternalMoveSingle::moleculeGroup() const
{
    return smplr->group();
}

/** Return the property used to find the flexibility of each molecule*/
const PropertyName& InternalMoveSingle::flexibilityProperty() const
{
    return flexibility_property;
}

/** Set the name of the property used to find the flexibility of each molecule */
void InternalMoveSingle::setFlexibilityProperty(const PropertyName &property)
{
    flexibility_property = property;
    Move::setProperty("flexibility", flexibility_property);
}

/** Set the random number generator used to generate the random
    number used for this move */
void InternalMoveSingle::setGenerator(const RanGenerator &rangenerator)
{
    MonteCarlo::setGenerator(rangenerator);
    smplr.edit().setGenerator(this->generator());
}

void InternalMoveSingle::setSynchronisedCoordinates(const MoleculeGroup &molgroup)
{
  synched_molgroup = molgroup;
}

const MoleculeGroup& InternalMoveSingle::synchronisedMols() const
{
  return synched_molgroup;
}

/** Internal function used to set the ensemble based on the
    passed temperature */
void InternalMoveSingle::_pvt_setTemperature(const Temperature &temperature)
{
    MonteCarlo::setEnsemble( Ensemble::NVT(temperature) );
}

/** Actually perform 'nmoves' moves of the molecules in the
    system 'system', optionally recording simulation statistics
    if 'record_stats' is true */
void InternalMoveSingle::move(System &system, int nmoves, bool record_stats)
{
    if (nmoves <= 0)
        return;

    InternalMoveSingle old_state(*this);
    System old_system_state(system);

    try
    {
        PropertyMap map;
        map.set("coordinates", this->coordinatesProperty());

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
//             int ndofs = nbonds + nangles + ndihedrals ;

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

//             if ( ndofs == 0 || maxvar < 0 || maxvar >= ndofs  )
//             {
//                 // We move everything in these cases..
//                 moved_bonds = flex_bonds;
//                 moved_angles = flex_angs;
//                 moved_dihedrals = flex_dihs;
//             }
//             else
//             {
//                 // draw a random number [0, bonds.size()+angles.size()+dihedrals.size())
//                 // to find matching dof.
//                 // Add it to moved_bonds or moved_angles or moved_dihedrals
//                 // if not already present
//                 int movecount = 0;

//                 while (movecount < maxvar)
//                 {
//                     int rand = this->generator().randInt(0, ndofs - 1);
//                     //qDebug() << " rand is " << rand;

//                     if ( rand < nbonds )
//                     {
//                         // it is a bond...
//                         const BondID &bond = flex_bonds.at( rand );
//                         if ( not moved_bonds.contains(bond) )
//                         {
//                             //qDebug() << " adding bond " << rand;
//                             moved_bonds.append(bond);
//                             ++movecount;
//                         }
//                     }
//                     else if ( rand < ( nbonds + nangles ) )
//                     {
//                         // it is an angle...
//                         const AngleID &angle = flex_angs.at( rand - nbonds );

//                         if ( not moved_angles.contains(angle) )
//                         {
//                             //qDebug() << " adding angle " << rand - nbonds;
//                             moved_angles.append(angle);
//                             ++movecount;
//                         }
//                     }
//                     else
//                     {
//                         // it is a dihedral...
//                         const DihedralID &dihedral = flex_dihs.at(
//                                                             rand - nbonds - nangles );

//                         if ( not moved_dihedrals.contains(dihedral) )
//                         {
//                             //qDebug() << " adding dihedral " << rand - nbonds - nangles;
//                             moved_dihedrals.append(dihedral);
//                             ++movecount;
//                         }
//                     }
//                 }
//             }

            // Now actually move the selected dofs
            Mover<Molecule> mol_mover = oldmol.molecule().move();

            // move the bonds of this molecule
            Length bond_delta;

            foreach (const BondID &bond, moved_bonds)
            {
                //const Length bond_delta_value = flex.bond_deltas[bond];
                double bond_delta_value = flex.delta(bond);
                bond_delta = Length( this->generator().rand(-bond_delta_value,
                                                             bond_delta_value) );

                mol_mover.change(bond, bond_delta);
            }

            // and the angles
            Angle angle_delta;

            foreach (const AngleID &angle, moved_angles)
            {
                //const Angle angle_delta_value = flex.angle_deltas[angle];
                double angle_delta_value = flex.delta(angle);
                angle_delta = Angle( this->generator().rand(-angle_delta_value,
                                                             angle_delta_value) );

                mol_mover.change(angle, angle_delta);
            }

            // and the torsions
            Angle dihedral_delta;

            foreach (const DihedralID &dihedral, moved_dihedrals)
            {
                // We rotate by picking the central bond of the dihedral to
                // handle concerted motions
                //BondID centralbond;
                //centralbond = BondID(dihedral.atom1(), dihedral.atom2());

                //const Angle angle_delta_value = flex.angle_deltas[dihedral];
                double angle_delta_value = flex.delta(dihedral);
                dihedral_delta =  Angle( this->generator().rand(-angle_delta_value,
                                                                 angle_delta_value) );
                //mol_mover.change(centralbond, dihedral_delta);
                if (this->generator().randBool())
                {
                    //move only this specific dihedral
                    mol_mover.change(dihedral, dihedral_delta);
                }
                else
                {
                    BondID centralbond;
                    centralbond = BondID(dihedral.atom1(), dihedral.atom2());
                    mol_mover.change(centralbond, dihedral_delta);
                }
            }

            //update the system with the new coordinates
            Molecule newmol = mol_mover.commit();

            //system.update(newmol);

	    //JM 03/11 test if molecule coordinates synched with others
	    // if so, for each synched mol, set coordinates to those of newmol
	    // This only works if all molecules have identical topology
	    // --> This should be a constraint ?

	    //qDebug() << "HELLO " << synched_molgroup.molecules().toString();

	    // Get the synched_molgroup from system using the name of the stored synched_molgroup
	    // We have to look up by names at numbers may have changed if for instance
	    // two replicas have been swapped

	    const MGName &sync_name = synched_molgroup.name();

	    synched_molgroup = system[ sync_name ];

	    if (not synched_molgroup.molecules().isEmpty())
	      {
		const PropertyName &coords_property = map["coordinates"];
		const AtomCoords &coords = newmol.data().property(coords_property)
                                         		  .asA<AtomCoords>();

		Molecules new_molecules = synched_molgroup.molecules();

		for (Molecules::const_iterator it = synched_molgroup.molecules().constBegin();
		     it != synched_molgroup.molecules().constEnd();
		     ++it)
		  {
		    //qDebug() << it;
		    new_molecules.update( it->molecule().edit()
					  .setProperty( map["coordinates"],
							coords )
					  .commit() ) ;
		    //qDebug() << molecule.toString();
		    //MolEditor editmol = molecule.edit()
		  }

		new_molecules.add(newmol);

		system.update(new_molecules);
	      }
	    else
	      {
		//update the system with the new coordinates
		system.update(newmol);
	      }

            //get the new bias on this molecule
            smplr.edit().updateFrom(system);

            new_bias = smplr.read().probabilityOf( PartialMolecule(newmol,
			                                                       oldmol.selection()) );

	    // JM DEBUG THERE IS A BUG AND THE INTRAGROUP FORCEFIELDS ENERGIES OF THE MOVED MOLECULE
	    // ARE NOT UPDATED ( BUT THE PERTURBED ENERGIES ARE CORRECT ! )
	    //system.mustNowRecalculateFromScratch();
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

const char* InternalMoveSingle::typeName()
{
    return QMetaType::typeName( qMetaTypeId<InternalMoveSingle>() );
}
