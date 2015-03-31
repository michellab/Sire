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

#include "qmmmelecembedpotential.h"
#include "qmprogram.h"

#include "latticecharges.h"

#include "SireUnits/units.h"

#include "SireMol/atomelements.h"

#include "SireBase/numberproperty.h"

#include "SireFF/errors.h"
#include "SireBase/errors.h"

#include "SireStream/datastream.h"

#include <QDebug>

using boost::tuples::tuple;

using namespace Squire;
using namespace SireMM;
using namespace SireMol;
using namespace SireVol;
using namespace SireBase;
using namespace SireUnits;
using namespace SireUnits::Dimension;
using namespace SireMM::detail;
using namespace SireStream;

namespace Squire
{
    template class QMMMPotential<QMPotential,InterCoulombPotential>;
}

static const RegisterMetaType<QMMMElecEmbedPotential> r_qmmm( MAGIC_ONLY, NO_ROOT,
                                                "Squire::QMMMElecEmbedPotential" );
                                                
/** Serialise to a binary datastream */
QDataStream SQUIRE_EXPORT &operator<<(QDataStream &ds,
                                      const QMMMElecEmbedPotential &qmmm)
{
    writeHeader(ds, r_qmmm, 2);
    
    ds << static_cast<const QMMMPotential<QMPotential,InterCoulombPotential>&>(qmmm)
       << qmmm.chg_sclfac;
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SQUIRE_EXPORT &operator>>(QDataStream &ds,
                                      QMMMElecEmbedPotential &qmmm)
{
    VersionID v = readHeader(ds, r_qmmm);
    
    if (v <= 2)
    {
        SharedDataStream sds(ds);
    
        ds >> static_cast<QMMMPotential<QMPotential,InterCoulombPotential>&>(qmmm);

        if (v == 2)
            ds >> qmmm.chg_sclfac;
        else
            qmmm.chg_sclfac = 1;

        qmmm.mergeProperties();
    }
    else
        throw version_error(v, "1", r_qmmm, CODELOC);
    
    return ds;
}

/** Merge the properties of the QM and MM forcefields */
void QMMMElecEmbedPotential::mergeProperties()
{
    //we are only interested in the space, qm program and switching function
    //properties
    props = Properties();
    
    props.setProperty("space", this->space());
    props.setProperty("switchingFunction", this->switchingFunction());
    props.setProperty("quantum program", this->quantumProgram());
    props.setProperty("zero energy", QMPotential::properties().property("zero energy"));
    props.setProperty("chargeScalingFactor", NumberProperty(chg_sclfac));
}

/** Constructor */
QMMMElecEmbedPotential::QMMMElecEmbedPotential()
                       : QMMMPotential<QMPotential,InterCoulombPotential>(),
                         chg_sclfac(1.0)
{
    this->mergeProperties();
}

/** Copy constructor */
QMMMElecEmbedPotential::QMMMElecEmbedPotential(const QMMMElecEmbedPotential &other)
                       : QMMMPotential<QMPotential,InterCoulombPotential>(other),
                         props(other.props), chg_sclfac(other.chg_sclfac)
{}

/** Destructor */
QMMMElecEmbedPotential::~QMMMElecEmbedPotential()
{}

/** Copy assignment operator */
QMMMElecEmbedPotential& 
QMMMElecEmbedPotential::operator=(const QMMMElecEmbedPotential &other)
{
    QMMMPotential<QMPotential,InterCoulombPotential>::operator=(other);
    props = other.props;
    chg_sclfac = other.chg_sclfac;
    
    return *this;
}

/** Return the space within which the molecules in this potential exist */
const Space& QMMMElecEmbedPotential::space() const
{
    return QMPotential::space();
}

/** Return the switching function that is used to implement the non-bonded
    cutoff */
const SwitchingFunction& QMMMElecEmbedPotential::switchingFunction() const
{
    return MMPotential::switchingFunction();
}

/** Return the handle to the quantum chemical program that is used 
    by this potential to calculate the QM energies and forces */
const QMProgram& QMMMElecEmbedPotential::quantumProgram() const
{
    return QMPotential::quantumProgram();
}

/** Return the absolute value of the energy which is considered
    as zero (on the relative energy scale used by this potential).
    A relative scale is used so that the QM energy can be shifted
    so that it is comparable to an MM energy */
MolarEnergy QMMMElecEmbedPotential::zeroEnergy() const
{
    return QMPotential::zeroEnergy();
}

/** Return the amount by which the MM charges are scaled in the QM/MM interaction */
double QMMMElecEmbedPotential::chargeScalingFactor() const
{
    return chg_sclfac;
}

/** Set the space within which all of the molecules in this potential
    will exist. This returns whether or not this changes the
    potential. */
bool QMMMElecEmbedPotential::setSpace(const Space &space)
{
    if (QMPotential::setSpace(space))
    {
        this->mergeProperties();
        return true;
    }
    else
        return false;
}

/** Set the switching function that will be used to implement the 
    non-bonded cutoff in the QM/MM interface */
bool QMMMElecEmbedPotential::setSwitchingFunction(const SwitchingFunction &switchfunc)
{
    if (MMPotential::setSwitchingFunction(switchfunc))
    {
        this->mergeProperties();
        return true;
    }
    else
        return false;
}

/** Set the handle to the quantum chemical program that will be
    used by this potential to calculate the QM energies and forces.
    This returns whether or not this changes this potential */
bool QMMMElecEmbedPotential::setQuantumProgram(const QMProgram &program)
{
    if (QMPotential::setQuantumProgram(program))
    {
        this->mergeProperties();
        return true;
    }
    else
        return false;
}

/** Set the absolute value of the energy which is considered
    as zero (on the relative energy scale used by this potential).
    A relative scale is used so that the QM energy can be shifted
    so that it is comparable to an MM energy */
bool QMMMElecEmbedPotential::setZeroEnergy(MolarEnergy zero_energy)
{
    if (QMPotential::setZeroEnergy(zero_energy))
    {
        this->mergeProperties();
        return true;
    }
    else
        return false;
}

/** Set the scaling factor for the MM charges in the QM/MM intermolecular interaction */
bool QMMMElecEmbedPotential::setChargeScalingFactor(double scale_factor)
{
    if (scale_factor != chg_sclfac)
    {
        chg_sclfac = scale_factor;
        this->mergeProperties();
        return true;
    }
    else
        return false;
}

/** Set the property called 'name' to the value 'value'

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
bool QMMMElecEmbedPotential::setProperty(const QString &name, const Property &value)
{
    if (not this->containsProperty(name))
        throw SireBase::missing_property( QObject::tr(
            "There is no property called \"%1\" in the potential %2. "
            "Available properties are [ %3 ].")
                .arg(name).arg(this->what())
                .arg(Sire::toString(props.propertyKeys())),
                    CODELOC );

    if (name == "chargeScalingFactor")
    {
        return this->setChargeScalingFactor( value.asA<NumberProperty>().value() );
    }
    else if (QMPotential::containsProperty(name))
    {
        if (QMPotential::setProperty(name, value))
        {
            this->mergeProperties();
            return true;
        }
        else
            return false;
    }
    else
    {
        if (MMPotential::setProperty(name, value))
        {
            this->mergeProperties();
            return true;
        }
        else
            return false;
    }
}

/** Return the property called 'name'

    \throw SireBase::missing_property
*/
const Property& QMMMElecEmbedPotential::property(const QString &name) const
{
    return props.property(name);
}

/** Return whether or not this potential contains the property at name 'name' */
bool QMMMElecEmbedPotential::containsProperty(const QString &name) const
{
    return props.hasProperty(name);
}

/** Return all of the properties of this potential */
const Properties& QMMMElecEmbedPotential::properties() const
{
    return props;
}

/** This converts the MM molecules in 'mmmols' into a set of lattice charges
    that surround the QM molecules in 'qmmols' */
LatticeCharges QMMMElecEmbedPotential::getLatticeCharges(const QMMolecules &qmmols,
                                                         const MMMolecules &mmmols,
                              QHash<MolNum,AtomIntProperty> *lattice_indicies) const
{
    if (qmmols.isEmpty() or mmmols.isEmpty())
    {
        return LatticeCharges();
    }

    //the QM molecules are already in this space
    
    //merge all of the atoms from the QM molecules into a single CoordGroup
    CoordGroup qmgroup;
    
    if (qmmols.count() == 1)
    {
        const QMMolecule &qmmol = qmmols.moleculesByIndex()[0];
        qmgroup = qmmol.coordinates().merge();
    }
    else
    {
        int nqmmols = qmmols.count();
        QVector<CoordGroup> qmgroups(nqmmols);
        
        const ChunkedVector<QMMolecule> &qmmols_array = qmmols.moleculesByIndex();
        CoordGroup *qmgroups_array = qmgroups.data();
        
        for (int i=0; i<nqmmols; ++i)
        {
            qmgroups_array[i] = qmmols_array[i].coordinates().merge();
        }
        
        qmgroup = CoordGroupArray(qmgroups).merge();
    }

    //now map all of the MM molecules into the same space as the QM molecules
    const Space &spce = this->space();
    const SwitchingFunction &switchfunc = this->switchingFunction();
    
    double cutoff = switchfunc.electrostaticCutoffDistance();
    
    int nmols = mmmols.count();
    const ChunkedVector<MMMolecule> &mmmols_array = mmmols.moleculesByIndex();
    
    //try to reserve enough space
    int nats = 0;
    for (int i=0; i<nmols; ++i)
    {
        nats += mmmols_array[i].coordinates().nCoords();
    }
    
    LatticeCharges lattice_charges;
    lattice_charges.reserve(nats);

    //now place the molecules' charges onto the lattice, recording
    //the lattice index of the atoms of the closest CutGroup copy
    if (lattice_indicies != 0)
    {
        *lattice_indicies = QHash<MolNum,AtomIntProperty>();
        lattice_indicies->reserve(nmols);
    }
    
    for (int i=0; i<nmols; ++i)
    {
        const MMMolecule &mmmol = mmmols_array[i];
        
        int ngroups = mmmol.coordinates().nCoordGroups();
        const CoordGroup *cgroup_array = mmmol.coordinates().constData();
        const MMParameters::Array *charge_array = mmmol.parameters()
                                                       .atomicParameters().constData();

        // nasty code - I need the atom elements and am going to have to assume they
        // are correct. I need to update the QMMM potential to properly get the
        // charge and element property from each atom so that this nasty hack is not needed.
        const AtomElements &elems = mmmol.molecule().molecule().property("element")
                                         .asA<AtomElements>();
        
        BOOST_ASSERT( ngroups == mmmol.parameters().atomicParameters().nArrays() );
        
        AtomIntProperty lattice_idxs;
        
        if (lattice_indicies != 0)
            lattice_idxs = AtomIntProperty(mmmol.molecule().data().info(), -1);
        
        for (int j=0; j<ngroups; ++j)
        {
            //get all copies of this molecule within the cutoff distance
            //of any QM atom
            QList< tuple<double,CoordGroup> > mapped_groups = 
                                   spce.getCopiesWithin(cgroup_array[j], qmgroup,
                                                        cutoff);

            const MMParameters::Array &group_charges = charge_array[j];
            const ChargeParameter *group_charges_array = group_charges.constData();

            double mindist = std::numeric_limits<double>::max();

            const CGIdx cgidx = mmmol.cgIdx(j);

            for (QList< tuple<double,CoordGroup> >::const_iterator
                                                        it = mapped_groups.constBegin();
                 it != mapped_groups.constEnd();
                 ++it)
            {
                const double sqrt_4pieps0 = std::sqrt(SireUnits::four_pi_eps0);

                //get any scaling feather factor for this group (and to convert
                //the charge from reduced units to mod_electrons)
                double scl = switchfunc.electrostaticScaleFactor( Length(it->get<0>()) )
                                   * sqrt_4pieps0 * chg_sclfac;
                                   
                if (scl == 0)
                    continue;
                
                //add the coordinates and charges
                const CoordGroup &mapped_group = it->get<1>();
                
                BOOST_ASSERT(mapped_group.count() == group_charges.count());
                
                const Vector *mapped_group_array = mapped_group.constData();

                bool index_this_group = false;
                if (it->get<0>() < mindist)
                {
                    mindist = it->get<0>();
                    index_this_group = true;
                }
                
                for (int k=0; k<mapped_group.count(); ++k)
                {
                    double chg = scl * group_charges_array[k].reduced_charge;
                    
                    if (chg != 0)
                    {
                        const CGAtomIdx cgatomidx(cgidx, Index(k));
                    
                        if (index_this_group and (lattice_indicies != 0))
                            lattice_idxs.set( cgatomidx, lattice_charges.count() );
                        
                        //lattice charges are electron charges, with coordinates
                        //in angstroms
                        lattice_charges.add( 
                                LatticeCharge(mapped_group_array[k],
                                              chg, elems[cgatomidx]) );
                    }
                }
            }
        }
        
        if (lattice_indicies != 0)
            lattice_indicies->insert(mmmol.molecule().number(), lattice_idxs);
    }
    
    //get the limit of the number of MM atoms for the used number of QM atoms
    const int num_mm_limit = quantumProgram().numberOfMMAtomsLimit(qmgroup.count());

    if (num_mm_limit > 0 and num_mm_limit < lattice_charges.count())
    {
        //create a QMap indexed by distance, as this will sort by distance
        QMultiMap<float,int> distances;
        
        Cartesian space;
        
        for (int i=0; i<lattice_charges.count(); ++i)
        {
            const Vector coords( lattice_charges[i].x(), lattice_charges[i].y(),
                                 lattice_charges[i].z() );
        
            float dist = space.minimumDistance( qmgroup.aaBox(), coords );

            distances.insert( dist, i );
        }
        
        int n_to_remove = lattice_charges.count() - num_mm_limit;
        
        QMultiMap<float,int>::const_iterator it = distances.constEnd();
        
        float max_distance = 0;
        
        for (int i=0; i<n_to_remove; ++i)
        {
            --it;
            lattice_charges.setCharge(it.value(), 0.0);
            max_distance = it.key();
        }

        //there are too many MM atoms. We have to remove MM atoms, starting
        //from the furthest ones out, until we are under the limit
        qDebug() << "WARNING: Number of MM atoms is too high." << lattice_charges.count()
                 << num_mm_limit;
        qDebug() << "All MM atoms beyond a distance of" << max_distance
                 << "A have had their charges set to 0.";
    }
    
    return lattice_charges;
}

/** Calculate the QM forces on the molecules in 'molecules' and add them 
    onto the forces in the passed force table (scaled by the optional
    scaling constant 'scale_force') */
void QMMMElecEmbedPotential::calculateForce(const QMMolecules &qmmols, 
                                            const MMMolecules &mmmols,
                                            ForceTable &forcetable, 
                                            double scale_force) const
{
    if (scale_force == 0)
        return;

    //map all of the molecules so that they are in this space
    QMMolecules mapped_qmmols = QMPotential::mapIntoSpace(qmmols);
    
    QHash<MolNum,AtomIntProperty> lattice_indicies;

    LatticeCharges charges = this->getLatticeCharges(mapped_qmmols, mmmols,
                                                     &lattice_indicies);

    QVector<Vector> lattice_forces = quantumProgram().calculateForce(
                                                    mapped_qmmols, charges,
                                                    forcetable, scale_force);

    //loop over all MMMolecules and see if they are in the forcetable

    //map the lattice forces back to the potentials on the molecules
    qDebug() << "WARNING - NEED TO MAP LATTICE FORCES BACK TO MM ATOMS";
    qDebug() << "YOUR SIMULATION IS BROKEN!!!";
}
                    
/** Calculate the QM forces on the molecules in 'molecules' and add them 
    onto the forces in the passed force table (scaled by the optional
    scaling constant 'scale_force') */
void QMMMElecEmbedPotential::calculateForce(const QMMolecules &qmmols,
                                            const MMMolecules &mmmols, 
                                            ForceTable &forcetable,
                                            const Symbol &symbol, 
                                            const Components &components,
                                            double scale_force) const
{
    if (symbol == components.total())
        this->calculateForce(qmmols, mmmols, forcetable, scale_force);
        
    else
        throw SireFF::missing_component( QObject::tr(
            "There is no force component in potential %1 - available "
            "components are %2.")
                .arg(this->what())
                .arg(components.total().toString()), CODELOC );
}

/** Calculate the QM energy of the molecules in 'qmmols' in the electrostatic
    field of the molecules in 'mmmols' and add it on to the energy 'nrg', optionally
    scaled by the scaling constant 'scale_energy' */
void QMMMElecEmbedPotential::calculateEnergy(const QMMolecules &qmmols, 
                                             const MMMolecules &mmmols,
                                             Energy &nrg, double scale_energy) const
{
    if (scale_energy == 0)
        return;

    //map all of the molecules so that they are in this space
    QMMolecules mapped_qmmols = QMPotential::mapIntoSpace(qmmols);

    LatticeCharges charges = this->getLatticeCharges(mapped_qmmols, mmmols);

    double qmnrg = this->quantumProgram().calculateEnergy(mapped_qmmols, charges);
    
    nrg += Energy(scale_energy * (qmnrg - QMPotential::zeroEnergy()));
}

/** Return the contents of the QM program command file that will be used
    to calculate the QM energy of the molecules in 'qmmols' surrounded
    in the field of the molecules in 'mmols' */
QString QMMMElecEmbedPotential::energyCommandFile(
                          const QMMMElecEmbedPotential::QMMolecules &qmmols,
                          const QMMMElecEmbedPotential::MMMolecules &mmmols) const
{
    //map all of the molecules so that they are in this space
    QMMolecules mapped_qmmols = QMPotential::mapIntoSpace(qmmols);

    LatticeCharges charges = this->getLatticeCharges(mapped_qmmols, mmmols);

    return this->quantumProgram().energyCommandFile(mapped_qmmols, charges);
}

/** Return the contents of the QM program command file that will be used
    to calculate the QM forces on the molecules 'qmmols', and on the 
    molecules 'mmols', where 'qmmols' are the QM molecules that are
    in the field of the MM molecules 'mmmols' */
QString QMMMElecEmbedPotential::forceCommandFile(
                          const QMMMElecEmbedPotential::QMMolecules &qmmols,
                          const QMMMElecEmbedPotential::MMMolecules &mmmols,
                          const ForceTable &forcetable) const
{
    //map all of the molecules so that they are in this space
    QMMolecules mapped_qmmols = QMPotential::mapIntoSpace(qmmols);

    LatticeCharges charges = this->getLatticeCharges(mapped_qmmols, mmmols);

    return this->quantumProgram().forceCommandFile(mapped_qmmols, charges, forcetable);
}

/** Return the contents of the QM program command file that will be used
    to calculate the QM fields around the molecules 'qmmols', and on the 
    molecules 'mmols', where 'qmmols' are the QM molecules that are
    in the field of the MM molecules 'mmmols' */
QString QMMMElecEmbedPotential::fieldCommandFile(
                          const QMMMElecEmbedPotential::QMMolecules &qmmols,
                          const QMMMElecEmbedPotential::MMMolecules &mmmols,
                          const FieldTable &fieldtable,
                          const SireFF::Probe &probe) const
{
    //map all of the molecules so that they are in this space
    QMMolecules mapped_qmmols = QMPotential::mapIntoSpace(qmmols);

    LatticeCharges charges = this->getLatticeCharges(mapped_qmmols, mmmols);

    return this->quantumProgram().fieldCommandFile(mapped_qmmols, charges, 
                                                   fieldtable, probe);
}

/** Return the contents of the QM program command file that will be used
    to calculate the QM potentials around the molecules 'qmmols', and on the 
    molecules 'mmols', where 'qmmols' are the QM molecules that are
    in the field of the MM molecules 'mmmols' */
QString QMMMElecEmbedPotential::potentialCommandFile(
                          const QMMMElecEmbedPotential::QMMolecules &qmmols,
                          const QMMMElecEmbedPotential::MMMolecules &mmmols,
                          const PotentialTable &pottable,
                          const SireFF::Probe &probe) const
{
    //map all of the molecules so that they are in this space
    QMMolecules mapped_qmmols = QMPotential::mapIntoSpace(qmmols);

    LatticeCharges charges = this->getLatticeCharges(mapped_qmmols, mmmols);

    return this->quantumProgram().potentialCommandFile(mapped_qmmols, charges, 
                                                       pottable, probe);
}

/** Calculate the QM field at the points in the passed potential table */
void QMMMElecEmbedPotential::calculateField(const QMMolecules &qmmols, 
                                            const MMMolecules &mmmols,
                                            FieldTable &fieldtable,
                                            const SireFF::Probe &probe,
                                            double scale_field) const
{
    if (scale_field == 0)
        return;

    //map all of the molecules so that they are in this space
    QMMolecules mapped_qmmols = QMPotential::mapIntoSpace(qmmols);

    LatticeCharges charges = this->getLatticeCharges(mapped_qmmols, mmmols);
    
    QVector<Vector> lattice_fields = quantumProgram().calculateField(
                                                    mapped_qmmols, charges, fieldtable, 
                                                    probe, scale_field);

    //map the lattice potentials back to the potentials on the molecules
    qDebug() << "WARNING - NEED TO MAP LATTICE FIELDS BACK TO MM ATOMS";
    qDebug() << "YOUR SIMULATION IS BROKEN!!!";
}
    
/** Calculate the QM field at the points in the passed potential table */
void QMMMElecEmbedPotential::calculateField(const QMMolecules &qmmols,
                                            const MMMolecules &mmmols,
                                            FieldTable &fieldtable,
                                            const SireFF::Probe &probe,
                                            const Symbol &symbol,
                                            const Components &components,
                                            double scale_field) const
{
    if (symbol == components.total())
        this->calculateField(qmmols, mmmols, fieldtable, probe, scale_field);
        
    else
        throw SireFF::missing_component( QObject::tr(
            "There is no field component in potential %1 - available "
            "components are %2.")
                .arg(this->what())
                .arg(components.total().toString()), CODELOC );
}

/** Calculate the QM potential at the points in the passed potential table */
void QMMMElecEmbedPotential::calculatePotential(const QMMolecules &qmmols, 
                                                const MMMolecules &mmmols,
                                                PotentialTable &pottable,
                                                const SireFF::Probe &probe,
                                                double scale_potential) const
{
    if (scale_potential == 0)
        return;

    //map all of the molecules so that they are in this space
    QMMolecules mapped_qmmols = QMPotential::mapIntoSpace(qmmols);

    QHash<MolNum,AtomIntProperty> lattice_indicies;

    LatticeCharges charges = this->getLatticeCharges(mapped_qmmols, mmmols,
                                                     &lattice_indicies);
    
    QVector<MolarEnergy> lattice_potentials = quantumProgram().calculatePotential(
                                                    mapped_qmmols, charges, pottable, 
                                                    probe, scale_potential);

    //loop over the MMMolecules and see if they are in the potential table
    int nmols = mmmols.count();
    const ChunkedVector<MMMolecule> &mmmols_array = mmmols.moleculesByIndex();
                      
    const MolarEnergy *potentials_array = lattice_potentials.constData();
                                        
    for (int i=0; i<nmols; ++i)
    {
        const MMMolecule &mmmol = mmmols_array[i];
        MolNum molnum = mmmol.number();
    
        if (lattice_indicies.contains(molnum) and pottable.contains(molnum))
        {
            const AtomIntProperty &indicies = *(lattice_indicies.constFind(molnum));
            
            MolPotentialTable &moltable = pottable.getTable(molnum);
            
            for (CGIdx cgidx(0); cgidx<indicies.nCutGroups(); ++cgidx)
            {
                for (Index atomidx(0); atomidx<indicies.nAtoms(cgidx); ++atomidx)
                {
                    CGAtomIdx cgatomidx(cgidx, atomidx);
                    
                    int idx = indicies.at(cgatomidx);
                    
                    if (idx >= 0)
                        moltable.add(cgatomidx, potentials_array[idx]);
                }
            }
        }
    }
}
    
/** Calculate the QM potential at the points in the passed potential table */
void QMMMElecEmbedPotential::calculatePotential(const QMMolecules &qmmols,
                                                const MMMolecules &mmmols,
                                                PotentialTable &pottable,
                                                const SireFF::Probe &probe,
                                                const Symbol &symbol,
                                                const Components &components,
                                                double scale_potential) const
{
    if (symbol == components.total())
        this->calculatePotential(qmmols, mmmols, pottable, probe, scale_potential);
        
    else
        throw SireFF::missing_component( QObject::tr(
            "There is no potential component in potential %1 - available "
            "components are %2.")
                .arg(this->what())
                .arg(components.total().toString()), CODELOC );
}
