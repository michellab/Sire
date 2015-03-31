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

#include "zmatmove.h"

#include "zmatrix.h"
#include "ensemble.h"

#include "SireSystem/system.h"

#include "SireMol/partialmolecule.h"
#include "SireMol/molecule.h"
#include "SireMol/moleditor.h"
#include "SireMol/mover.hpp"
#include "SireMol/atomidx.h"

#include "SireUnits/dimensions.h"
#include "SireUnits/temperature.h"

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

static const RegisterMetaType<ZMatMove> r_zmatmove;

/** Serialise to a binary datastream */
QDataStream SIREMOVE_EXPORT &operator<<(QDataStream &ds, const ZMatMove &zmatmove)
{
    writeHeader(ds, r_zmatmove, 2);
    
    SharedDataStream sds(ds);
    
    sds << zmatmove.smplr << zmatmove.zmatrix_property
        << zmatmove.sync_bonds << zmatmove.sync_angles
        << zmatmove.sync_dihedrals
        << static_cast<const MonteCarlo&>(zmatmove);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOVE_EXPORT &operator>>(QDataStream &ds, ZMatMove &zmatmove)
{
    VersionID v = readHeader(ds, r_zmatmove);
    
    if (v == 2)
    {
        SharedDataStream sds(ds);
        
        sds >> zmatmove.smplr >> zmatmove.zmatrix_property
            >> zmatmove.sync_bonds >> zmatmove.sync_angles
            >> zmatmove.sync_dihedrals
            >> static_cast<MonteCarlo&>(zmatmove);
    }
    else if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> zmatmove.smplr >> zmatmove.zmatrix_property
            >> static_cast<MonteCarlo&>(zmatmove);
            
        zmatmove.sync_bonds = false;
        zmatmove.sync_angles = false;
        zmatmove.sync_dihedrals = false;
    }
    else
        throw version_error(v, "1", r_zmatmove, CODELOC);
        
    return ds;
}

/** Null constructor */
ZMatMove::ZMatMove(const PropertyMap &map) 
         : ConcreteProperty<ZMatMove,MonteCarlo>(map),
           sync_bonds(false), sync_angles(false),
           sync_dihedrals(false)
{
    zmatrix_property = map["z-matrix"];
    MonteCarlo::setEnsemble( Ensemble::NVT(25*celsius) );
}

/** Construct the z-matrix move for the passed group of molecules */
ZMatMove::ZMatMove(const MoleculeGroup &molgroup, const PropertyMap &map)
         : ConcreteProperty<ZMatMove,MonteCarlo>(),
           smplr( UniformSampler(molgroup) ),
           sync_bonds(false), sync_angles(false),
           sync_dihedrals(false)
{
    zmatrix_property = map["z-matrix"];
    MonteCarlo::setEnsemble( Ensemble::NVT(25*celsius) );
    smplr.edit().setGenerator( this->generator() );
}

/** Construct the z-matrix move that samples molecules from the
    passed sampler */
ZMatMove::ZMatMove(const Sampler &sampler, const PropertyMap &map)
         : ConcreteProperty<ZMatMove,MonteCarlo>(),
           smplr(sampler),
           sync_bonds(false), sync_angles(false),
           sync_dihedrals(false)
{
    zmatrix_property = map["z-matrix"];
    MonteCarlo::setEnsemble( Ensemble::NVT(25*celsius) );
    smplr.edit().setGenerator( this->generator() );
}

/** Copy constructor */
ZMatMove::ZMatMove(const ZMatMove &other)
         : ConcreteProperty<ZMatMove,MonteCarlo>(other),
           smplr(other.smplr),
           zmatrix_property(other.zmatrix_property),
           sync_bonds(other.sync_bonds), sync_angles(other.sync_angles),
           sync_dihedrals(other.sync_dihedrals)
{}

/** Destructor */
ZMatMove::~ZMatMove()
{}

/** Copy assignment operator */
ZMatMove& ZMatMove::operator=(const ZMatMove &other)
{
    if (this != &other)
    {
        MonteCarlo::operator=(other);
        smplr = other.smplr;
        zmatrix_property = other.zmatrix_property;
        sync_bonds = other.sync_bonds;
        sync_angles = other.sync_angles;
        sync_dihedrals = other.sync_dihedrals;
    }
    
    return *this;
}

/** Comparison operator */
bool ZMatMove::operator==(const ZMatMove &other) const
{
    return MonteCarlo::operator==(other) and smplr == other.smplr and
           sync_bonds == other.sync_bonds and sync_angles == other.sync_angles and
           sync_dihedrals == other.sync_dihedrals and
           zmatrix_property == other.zmatrix_property;
}

/** Comparison operator */
bool ZMatMove::operator!=(const ZMatMove &other) const
{
    return not ZMatMove::operator==(other);
}

/** Return a string representation of this move */
QString ZMatMove::toString() const
{
    return QObject::tr("ZMatMove( nAccepted() = %1 nRejected() == %2 )")
                .arg( this->nAccepted() )
                .arg( this->nRejected() );
}

/** Set the sampler used to sample molecules for this move */
void ZMatMove::setSampler(const Sampler &sampler)
{
    smplr = sampler;
    smplr.edit().setGenerator( this->generator() );
}

/** Set the sampler so that it draws molecules uniformly from the
    molecule group 'molgroup' */
void ZMatMove::setSampler(const MoleculeGroup &molgroup)
{
    this->setSampler( UniformSampler(molgroup) );
}

/** Return the sampler used to sample molecules to move */
const Sampler& ZMatMove::sampler() const
{
    return smplr;
}

/** Return the molecule group that is sampled for this move */
const MoleculeGroup& ZMatMove::moleculeGroup() const
{
    return smplr->group();
}

/** Return the property used to find the z-matrix of each molecule */
const PropertyName& ZMatMove::zmatrixProperty() const
{
    return zmatrix_property;
}
    
/** Set the name of the property used to find the z-matrix of each molecule */
void ZMatMove::setZMatrixProperty(const PropertyName &property)
{
    zmatrix_property = property;
    Move::setProperty("z-matrix", zmatrix_property);
}

/** Set whether or not to synchronise all motion for all molecules 
    in the group */
void ZMatMove::setSynchronisedMotion(bool on)
{
    sync_bonds = on;
    sync_angles = on;
    sync_dihedrals = on;
}

/** Set whether or not to synchronise all bond moves for all molecules */
void ZMatMove::setSynchronisedBonds(bool on)
{
    sync_bonds = on;
}

/** Set whether or not to synchronise all angle moves for all molecules */
void ZMatMove::setSynchronisedAngles(bool on)
{
    sync_angles = on;
}

/** Set whether or not to synchronise all dihedral moves for all molecules */
void ZMatMove::setSynchronisedDihedrals(bool on)
{
    sync_dihedrals = on;
}

/** Return whether or not all moves for all molecules are synchronised */
bool ZMatMove::synchronisedMotion() const
{
    return sync_bonds and sync_angles and sync_dihedrals;
}

/** Return whether or not all bond moves for all molecules
    are synchronised */
bool ZMatMove::synchronisedBonds() const
{
    return sync_bonds;
}

/** Return whether or not all angle moves for all molecules
    are synchronised */
bool ZMatMove::synchronisedAngles() const
{
    return sync_angles;
}

/** Return whether or not all dihedral moves for all molecules
    are synchronised */
bool ZMatMove::synchronisedDihedrals() const
{
    return sync_dihedrals;
}

/** Set the random number generator used to generate the random
    number used for this move */
void ZMatMove::setGenerator(const RanGenerator &rangenerator)
{
    MonteCarlo::setGenerator(rangenerator);
    smplr.edit().setGenerator(this->generator());
}

/** Internal function used to set the ensemble based on the
    passed temperature */
void ZMatMove::_pvt_setTemperature(const Temperature &temperature)
{
    MonteCarlo::setEnsemble( Ensemble::NVT(temperature) );
}

/** Internal function used to move the bond, angle and dihedral
    that is used to build the atom 'atom' in the z-matrix 'zmatrix' */
void ZMatMove::move(AtomIdx atom, ZMatrixCoords &zmatrix,
                    QHash< AtomIdx, tuple<Length,Angle,Angle> > &saved_deltas)
{
    // first generate the amounts by which to change the 
    // bond, angle and dihedral values
    Length bonddelta;
    Angle angledelta, dihedraldelta;

    if ( sync_bonds or sync_angles or sync_dihedrals )
    {
        //we are synchronising bonds, angles or dihedrals, so
        //we may need to look up previous values for previous molecules
    
        if (saved_deltas.contains(atom))
        {
            bonddelta = saved_deltas.value(atom).get<0>();
            angledelta = saved_deltas.value(atom).get<1>();
            dihedraldelta = saved_deltas.value(atom).get<2>();
            
            if ((not sync_bonds) and bonddelta.value() != 0)
            {
                bonddelta = Length( this->generator().rand(-bonddelta.value(),
                                                            bonddelta.value() ) );
            }
            
            if ((not sync_angles) and angledelta.value() != 0)
            {
                angledelta = Angle( this->generator().rand(-angledelta.value(),
                                                            angledelta.value() ) );
            }
            
            if ((not sync_dihedrals) and dihedraldelta.value() != 0)
            {
                dihedraldelta = Angle( this->generator().rand(-dihedraldelta.value(),
                                                                dihedraldelta.value() ) );
            }
        }
        else
        {
            bonddelta = zmatrix.bondDelta(atom);
            angledelta = zmatrix.angleDelta(atom);
            dihedraldelta = zmatrix.dihedralDelta(atom);
            
            if (sync_bonds and bonddelta.value() != 0)
            {
                bonddelta = Length( this->generator().rand(-bonddelta.value(),
                                                            bonddelta.value() ) );
            }
            
            if (sync_angles and angledelta.value() != 0)
            {
                angledelta = Angle( this->generator().rand(-angledelta.value(),
                                                            angledelta.value() ) );
            }
            
            if (sync_dihedrals and dihedraldelta.value() != 0)
            {
                dihedraldelta = Angle( this->generator().rand(-dihedraldelta.value(),
                                                                dihedraldelta.value() ) );
            }
            
            saved_deltas[atom] = tuple<Length,Angle,Angle>(bonddelta, angledelta,
                                                           dihedraldelta);
        }
    }
    else
    {
        bonddelta = zmatrix.bondDelta(atom);
        angledelta = zmatrix.angleDelta(atom);
        dihedraldelta = zmatrix.dihedralDelta(atom);
        
        if (bonddelta.value() != 0)
        {
            bonddelta = Length( this->generator().rand(-bonddelta.value(),
                                                        bonddelta.value() ) );
        }
        
        if (angledelta.value() != 0)
        {
            angledelta = Angle( this->generator().rand(-angledelta.value(),
                                                        angledelta.value() ) );
        }
        
        if (dihedraldelta.value() != 0)
        {
            dihedraldelta = Angle( this->generator().rand(-dihedraldelta.value(),
                                                              dihedraldelta.value() ) );
        }
    }

    // now that the bond, angle and dihedral delta have been evaluated,
    // change the zmatrix
    if (bonddelta.value() != 0)
    {
        zmatrix.moveBond(atom, bonddelta);
    }
    
    if (angledelta.value() != 0)
    {
        zmatrix.moveAngle(atom, angledelta);
    }

    if (dihedraldelta.value() != 0)
    {
        zmatrix.moveDihedral(atom, dihedraldelta);
    }
}

/** Actually perform 'nmoves' moves of the molecules in the 
    system 'system', optionally recording simulation statistics
    if 'record_stats' is true */
void ZMatMove::move(System &system, int nmoves, bool record_stats)
{
    if (nmoves <= 0)
        return;
      
    //save our, and the system's, current state
    ZMatMove old_state(*this);
    
    System old_system_state(system);
    
    try
    {
        const PropertyMap &map = Move::propertyMap();
        
        for (int i=0; i<nmoves; ++i)
        {
            //get the old energy of the system
            double old_nrg = system.energy( this->energyComponent() );
                
            //save the old system and sampler
            System old_system(system);
            SamplerPtr old_sampler(smplr);

            QHash< AtomIdx,tuple<Length,Angle,Angle> > saved_deltas;

            double old_bias = 1;
            double new_bias = 1;

            if (sync_bonds and sync_angles and sync_dihedrals)
            {
                //move all of everything!
                const Molecules &molecules = smplr.read().group().molecules();
                
                Molecules new_molecules = molecules;
                
                for (Molecules::const_iterator it = molecules.constBegin();
                     it != molecules.constEnd();
                     ++it)
                {
                    ZMatrixCoords zmatrix( *it, map );

                    //move the internal coordinates of selected atoms in the 
                    //z-matrix
                    AtomSelection selected_atoms = it->selection();
            
                    if (selected_atoms.selectedAll())
                    {
                        //move everything
                        for (QHash<AtomIdx,int>::const_iterator 
                                                    it2 = zmatrix.index().constBegin();
                             it2 != zmatrix.index().constEnd();
                             ++it2)
                        {
                            this->move(it2.key(), zmatrix, saved_deltas);
                        }
                    }
                    else
                    {
                        //move only the selected atoms
                        for (QHash<AtomIdx,int>::const_iterator 
                                                    it2 = zmatrix.index().constBegin();
                             it2 != zmatrix.index().constEnd();
                             ++it2)
                        {
                            if (selected_atoms.selected(it2.key()))
                                this->move(it2.key(), zmatrix, saved_deltas);
                        }
                    }

                    new_molecules.update( it->molecule().edit()
                                            .setProperty( map["coordinates"], 
                                                          zmatrix.toCartesian() )
                                            .commit() );
                }
                
                system.update(new_molecules);
            }
            else if (sync_bonds or sync_angles or sync_dihedrals)
            {
                //move some of everything, and all of just one molecule
                throw SireError::incomplete_code( QObject::tr(
                        "This code needs to be written!!!"), CODELOC );
            }
            else
            {
                //move all of just one molecule
    
                //update the sampler with the latest version of the molecules
                smplr.edit().updateFrom(system);

                //randomly select a molecule to move
                tuple<PartialMolecule,double> mol_and_bias = smplr.read().sample();

                const PartialMolecule &oldmol = mol_and_bias.get<0>();
                old_bias = mol_and_bias.get<1>();

                ZMatrixCoords zmatrix( oldmol.molecule(), map );

                //move the internal coordinates of selected atoms in the 
                //z-matrix
                AtomSelection selected_atoms = oldmol.selection();
                AtomCoords coords;
            
                if (selected_atoms.selectedAll())
                {
                    //move everything
                    for (QHash<AtomIdx,int>::const_iterator 
                                                    it = zmatrix.index().constBegin();
                         it != zmatrix.index().constEnd();
                         ++it)
                    {
                        this->move(it.key(), zmatrix, saved_deltas);
                    }
                    
                    coords = zmatrix.toCartesian();
                }
                else
                {
                    //move only the selected atoms
                    for (QHash<AtomIdx,int>::const_iterator 
                                                it = zmatrix.index().constBegin();
                         it != zmatrix.index().constEnd();
                         ++it)
                    {
                        if (selected_atoms.selected(it.key()))
                        {
                            this->move(it.key(), zmatrix, saved_deltas);
                        }
                    }
                    
                    AtomCoords oldcoords = oldmol.property(map["coordinates"]).asA<AtomCoords>();
                    AtomCoords newcoords = zmatrix.toCartesian();
                    
                    foreach (CGIdx cgidx, selected_atoms.selectedCutGroups())
                    {
                        if (selected_atoms.selectedAll(cgidx))
                        {
                            oldcoords = oldcoords.set(cgidx, newcoords[cgidx]);
                        }
                        else
                        {
                            foreach (Index i, selected_atoms.selectedAtoms(cgidx))
                            {
                                CGAtomIdx cgatomidx(cgidx,i);
                                oldcoords = oldcoords.set(cgatomidx, newcoords[cgatomidx]);
                            }
                        }
                    }
                    
                    coords = oldcoords;
                }

                Molecule newmol = oldmol.molecule().edit()
                                        .setProperty( map["coordinates"], coords )
                                        .commit();
                        
                //update the system with the new coordinates
                system.update(newmol);

                //get the new bias on this molecule
                smplr.edit().updateFrom(system);
        
                new_bias = smplr.read().probabilityOf( PartialMolecule(newmol,
                                                       oldmol.selection()) );
            }

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

const char* ZMatMove::typeName()
{
    return QMetaType::typeName( qMetaTypeId<ZMatMove>() );
}
