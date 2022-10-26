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

#include "polarisecharges.h"
#include "delta.h"

#include "SireMaths/nmatrix.h"
#include "SireMaths/nvector.h"

#include "SireMol/connectivity.h"
#include "SireMol/atomcoords.h"
#include "SireMol/atomcharges.h"
#include "SireMol/atomenergies.h"
#include "SireMol/atompolarisabilities.h"
#include "SireMol/molecule.h"
#include "SireMol/moleditor.h"
#include "SireMol/moleculedata.h"
#include "SireMol/atomselection.h"
#include "SireMol/mgname.h"
#include "SireMol/core.h"

#include "SireFF/probe.h"
#include "SireFF/potentialtable.h"

#include "SireMM/cljprobe.h"

#include "SireSystem/system.h"

#include "SireUnits/units.h"
#include "SireUnits/convert.h"

#include "SireBase/refcountdata.h"

#include "SireBase/errors.h"
#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireSystem;
using namespace SireMol;
using namespace SireFF;
using namespace SireMM;
using namespace SireUnits;
using namespace SireUnits::Dimension;
using namespace SireMaths;
using namespace SireBase;
using namespace SireStream;

/////////////
///////////// Implementation of PolariseChargesData
/////////////

namespace SireSystem
{
    namespace detail
    {
        class PolariseChargesData : public RefCountData
        {
        public:
            PolariseChargesData() : RefCountData()
            {}

            PolariseChargesData(const MoleculeView &molview,
                                const PropertyName &coords_property,
                                const PropertyName &connectivity_property,
                                const PropertyName &polarise_property);

            PolariseChargesData(const PolariseChargesData &other)
                  : RefCountData(),
                    xx_matricies(other.xx_matricies),
                    inv_xx_matricies(other.inv_xx_matricies),
                    connectivity(other.connectivity),
                    selected_atoms(other.selected_atoms),
                    molversion(other.molversion)
            {}

            ~PolariseChargesData()
            {}

            PolariseChargesData& operator=(const PolariseChargesData &other)
            {
                if (this != &other)
                {
                    xx_matricies = other.xx_matricies;
                    inv_xx_matricies = other.inv_xx_matricies;
                    connectivity = other.connectivity;
                    selected_atoms = other.selected_atoms;
                    molversion = other.molversion;
                }

                return *this;
            }

            /** The matrix holding (1/alpha) * (X X^T) */
            QVector<NMatrix> xx_matricies;

            /** The matrix holding alpha * (X X^T)**-1 */
            QVector<NMatrix> inv_xx_matricies;

            /** The connectivity of the molecule */
            Connectivity connectivity;

            /** The selection of atoms to be polarised - this
                is empty if all of the atoms are selected */
            AtomSelection selected_atoms;

            /** The version number of the molecule for which this
                data has been calculated */
            quint64 molversion;
        };
    }
}

using namespace SireSystem::detail;

static void calculateAtomMatrix(AtomIdx atomidx, const AtomCoords &coords,
                                const Volume &polarisability,
                                const Connectivity &connectivity,
                                const MoleculeInfoData &molinfo,
                                QVarLengthArray<Vector,10> &X,
                                NMatrix &xx, NMatrix &inv_xx)
{
    //get the atoms bonded to this atom
    const QSet<AtomIdx> &bonded_atoms = connectivity.connectionsTo(atomidx);

    const int nbonded = bonded_atoms.count();

    //qDebug() << atomidx << "Bonded" << bonded_atoms;

    const double alpha = four_pi_eps0 * polarisability.value();

    if (nbonded == 0 or alpha < 1e-6)
        return;

    //qDebug() << "coords" << coords[ molinfo.cgAtomIdx(atomidx) ].toString();
    //qDebug() << "alpha" << polarisability.toString() << four_pi_eps0 << alpha;

    //construct the X matrix (matrix of vectors from connected
    //atoms to the this atom
    {
        X.resize(nbonded);

        const Vector &atom_coords = coords[ molinfo.cgAtomIdx(atomidx) ];

        int i = 0;
        foreach (AtomIdx bonded_atom, bonded_atoms)
        {
            X[i] = coords[ molinfo.cgAtomIdx(bonded_atom) ] - atom_coords;

            ++i;
        }
    }

    //qDebug() << "X";
    //for (int i=0; i<X.count(); ++i)
    //{
    //    qDebug() << i << X[i].toString();
    //}

    //now construct (1/alpha)( X^T X )
    xx = NMatrix(nbonded, nbonded);

    const double one_over_alpha = 1 / alpha;

    for (int i=0; i<nbonded; ++i)
    {
        xx(i,i) = one_over_alpha * Vector::dot(X[i], X[i]);

        for (int j=i+1; j<nbonded; ++j)
        {
            const double i_dot_j = one_over_alpha * Vector::dot(X[i], X[j]);
            xx(i,j) = i_dot_j;
            xx(j,i) = i_dot_j;
        }
    }

    //qDebug() << "(1/alpha) * X^T X\n" << xx.toString();
    //qDebug() << "(X^T X)-1\n" << xx.inverse().toString();

    //now construct alpha ( X^T X )**-1
    inv_xx = xx.inverse();

    //qDebug() << "1/alpha (X^T X)\n" << xx.toString();
    //qDebug() << "alpha (X^T X)-1\n" << inv_xx.toString();
}

PolariseChargesData::PolariseChargesData(const MoleculeView &molview,
                                         const PropertyName &coords_property,
                                         const PropertyName &connectivity_property,
                                         const PropertyName &polarise_property)
                    : RefCountData()
{
    const AtomCoords &coords = molview.data().property(coords_property)
                                             .asA<AtomCoords>();

    const MoleculeInfoData &molinfo = molview.data().info();

    connectivity = molview.data().property(connectivity_property)
                                 .asA<Connectivity>();

    const AtomPolarisabilities &polarise = molview.data().property(polarise_property)
                                                         .asA<AtomPolarisabilities>();

    xx_matricies = QVector<NMatrix>(coords.nAtoms());
    inv_xx_matricies = QVector<NMatrix>(coords.nAtoms());
    xx_matricies.squeeze();
    inv_xx_matricies.squeeze();

    //this is the matrix of r vectors
    QVarLengthArray<Vector,10> X;

    if (molview.selectedAll())
    {
        //loop over all of the atoms
        int nats = coords.nAtoms();

        for (AtomIdx i(0); i<nats; ++i)
        {
            calculateAtomMatrix(i, coords, polarise[molinfo.cgAtomIdx(i)],
                                connectivity, molinfo, X,
                                xx_matricies[i], inv_xx_matricies[i]);
        }
    }
    else
    {
        selected_atoms = molview.selection();

        int ncgroups = selected_atoms.nCutGroups();

        const MoleculeInfoData &molinfo = molview.data().info();

        if (selected_atoms.selectedAllCutGroups())
        {
            for (CGIdx i(0); i<ncgroups; ++i)
            {
                if (selected_atoms.selectedAll(i))
                {
                    for (Index j(0); j<molinfo.nAtoms(i); ++j)
                    {
                        CGAtomIdx cgatomidx(i,j);

                        AtomIdx atomidx = molinfo.atomIdx(cgatomidx);

                        calculateAtomMatrix(atomidx, coords, polarise[cgatomidx],
                                            connectivity, molinfo, X,
                                            xx_matricies[atomidx],
                                            inv_xx_matricies[atomidx]);
                    }
                }
                else
                {
                    foreach (Index j, selected_atoms.selectedAtoms(i))
                    {
                        CGAtomIdx cgatomidx(i,j);

                        AtomIdx atomidx = molinfo.atomIdx(cgatomidx);

                        calculateAtomMatrix(atomidx, coords, polarise[cgatomidx],
                                            connectivity, molinfo, X,
                                            xx_matricies[atomidx],
                                            inv_xx_matricies[atomidx]);
                    }
                }
            }
        }
        else
        {
            foreach (CGIdx i, selected_atoms.selectedCutGroups())
            {
                if (selected_atoms.selectedAll(i))
                {
                    for (Index j(0); j<molinfo.nAtoms(i); ++j)
                    {
                        CGAtomIdx cgatomidx(i,j);

                        AtomIdx atomidx = molinfo.atomIdx(cgatomidx);

                        calculateAtomMatrix(atomidx, coords, polarise[cgatomidx],
                                            connectivity, molinfo, X,
                                            xx_matricies[atomidx],
                                            inv_xx_matricies[atomidx]);
                    }
                }
                else
                {
                    foreach (Index j, selected_atoms.selectedAtoms(i))
                    {
                        CGAtomIdx cgatomidx(i,j);

                        AtomIdx atomidx = molinfo.atomIdx(cgatomidx);

                        calculateAtomMatrix(atomidx, coords, polarise[cgatomidx],
                                            connectivity, molinfo, X,
                                            xx_matricies[atomidx],
                                            inv_xx_matricies[atomidx]);
                    }
                }
            }
        }
    }
}

/////////////
///////////// Implementation of PolariseCharges
/////////////

static const RegisterMetaType<PolariseCharges> r_polarise_charges;

QDataStream &operator<<(QDataStream &ds,
                                          const PolariseCharges &polchgs)
{
    writeHeader(ds, r_polarise_charges, 2);

    SharedDataStream sds(ds);

    sds << polchgs.field_component << polchgs.field_probe
        << polchgs.convergence_limit
        << static_cast<const ChargeConstraint&>(polchgs);

    return ds;
}

QDataStream &operator>>(QDataStream &ds,
                                          PolariseCharges &polchgs)
{
    VersionID v = readHeader(ds, r_polarise_charges);

    if (v == 2)
    {
        SharedDataStream sds(ds);

        polchgs = PolariseCharges();

        sds >> polchgs.field_component >> polchgs.field_probe
            >> polchgs.convergence_limit
            >> static_cast<ChargeConstraint&>(polchgs);
    }
    else if (v == 1)
    {
        SharedDataStream sds(ds);

        polchgs = PolariseCharges();

        sds >> polchgs.field_component >> polchgs.field_probe
            >> static_cast<ChargeConstraint&>(polchgs);

        polchgs.convergence_limit = 1e-3;
    }
    else
        throw version_error( v, "1", r_polarise_charges, CODELOC );

    return ds;
}

void PolariseCharges::setProbe(const Probe &probe)
{
    if (probe.isA<CoulombProbe>())
    {
        field_probe = probe.asA<CoulombProbe>();
    }
    else if (probe.isA<CLJProbe>())
    {
        field_probe = CoulombProbe( probe.asA<CLJProbe>() );
    }
    else
        throw SireError::incompatible_error( QObject::tr(
                "You can only use a CoulombProbe or CLJProbe with "
                "the PolariseCharges constraint - you cannot use %1.")
                    .arg(probe.toString()), CODELOC );
}

/** Null constructor */
PolariseCharges::PolariseCharges()
                : ConcreteProperty<PolariseCharges,ChargeConstraint>(),
                  convergence_limit(1e-3)
{}

/** Construct a constraint that uses the total energy field and a
    single unit charge probe to polarise the molecules in 'molgroup' */
PolariseCharges::PolariseCharges(const MoleculeGroup &molgroup,
                                 const PropertyMap &map)
                : ConcreteProperty<PolariseCharges,ChargeConstraint>(molgroup, map),
                  field_component(ForceFields::totalComponent()),
                  convergence_limit(1e-3)
{
    this->setProbe( CoulombProbe( 1*mod_electron ) );
}

/** Construct a constraint that uses the total energy field and
    the passed probe to polarise the molecules in 'molgroup' */
PolariseCharges::PolariseCharges(const MoleculeGroup &molgroup,
                                 const Probe &probe, const PropertyMap &map)
                : ConcreteProperty<PolariseCharges,ChargeConstraint>(molgroup, map),
                  field_component(ForceFields::totalComponent()),
                  convergence_limit(1e-3)
{
    this->setProbe(probe);
}

/** Construct a constraint that uses the field represented by 'field_component'
    and a single unit charge to polarise the molecules in 'molgroup' */
PolariseCharges::PolariseCharges(const MoleculeGroup &molgroup,
                                 const Symbol &fieldcomp, const PropertyMap &map)
                : ConcreteProperty<PolariseCharges,ChargeConstraint>(molgroup, map),
                  field_component(fieldcomp),
                  convergence_limit(1e-3)
{
    this->setProbe( CoulombProbe( 1*mod_electron ) );
}

/** Construct a constraint that uses the field represented by 'field_component'
    and the passed probe to polarise the molecules in 'molgroup' */
PolariseCharges::PolariseCharges(const MoleculeGroup &molgroup,
                                 const Symbol &fieldcomp, const Probe &probe,
                                 const PropertyMap &map)
                : ConcreteProperty<PolariseCharges,ChargeConstraint>(molgroup, map),
                  field_component(fieldcomp), convergence_limit(1e-3)
{
    this->setProbe(probe);
}

/** Copy constructor */
PolariseCharges::PolariseCharges(const PolariseCharges &other)
                : ConcreteProperty<PolariseCharges,ChargeConstraint>(other),
                  field_component(other.field_component),
                  field_probe(other.field_probe),
                  moldata(other.moldata), changed_mols(other.changed_mols),
                  convergence_limit(other.convergence_limit)
{}

/** Destructor */
PolariseCharges::~PolariseCharges()
{}

/** Copy assignment operator */
PolariseCharges& PolariseCharges::operator=(const PolariseCharges &other)
{
    if (this != &other)
    {
        field_component = other.field_component;
        field_probe = other.field_probe;
        moldata = other.moldata;
        changed_mols = other.changed_mols;
        convergence_limit = other.convergence_limit;

        ChargeConstraint::operator=(other);
    }

    return *this;
}

/** Comparison operator */
bool PolariseCharges::operator==(const PolariseCharges &other) const
{
    return (this == &other) or
           (field_component == other.field_component and
            field_probe == other.field_probe and
            convergence_limit == other.convergence_limit and
            ChargeConstraint::operator==(other));
}

/** Comparison operator */
bool PolariseCharges::operator!=(const PolariseCharges &other) const
{
    return not PolariseCharges::operator==(other);
}

const char* PolariseCharges::typeName()
{
    return QMetaType::typeName( qMetaTypeId<PolariseCharges>() );
}

PolariseCharges* PolariseCharges::clone() const
{
    return new PolariseCharges(*this);
}

/** Return a string representation of this constraint */
QString PolariseCharges::toString() const
{
    return "SireSystem::PolariseCharges";
}

/** Return the component of the forcefield that is used to
    calculate the electrostatic field on the atoms to be
    polarised */
const Symbol& PolariseCharges::fieldComponent() const
{
    return field_component;
}

/** Return the probe that is used to calculate the electrostatic
    field on the atoms to be polarised */
const CoulombProbe& PolariseCharges::probe() const
{
    return field_probe;
}

/** Set the convergence limit of the calculation */
void PolariseCharges::setConvergenceLimit(double limit)
{
    convergence_limit = limit;
}

/** Return the convergence limit of the calculation */
double PolariseCharges::convergenceLimit() const
{
    return convergence_limit;
}

static void calculateCharges(AtomIdx atomidx,
                             const MolPotentialTable &moltable,
                             const NMatrix &inv_alpha_XX,
                             const NMatrix &alpha_inv_XX,
                             const Connectivity &connectivity,
                             const MoleculeInfoData &molinfo,
                             AtomCharges &induced_charges,
                             AtomEnergies &selfpol_nrgs)
{
    if (alpha_inv_XX.nRows() == 0)
        //this atom has a zero polarisability, or no bonded atoms
        return;

    CGAtomIdx cgatomidx = molinfo.cgAtomIdx(atomidx);

    const QSet<AtomIdx> &bonded_atoms = connectivity.connectionsTo(atomidx);

    int nbonded = bonded_atoms.count();

    BOOST_ASSERT(nbonded == alpha_inv_XX.nRows());

    double phi_a = moltable[cgatomidx.cutGroup()][cgatomidx.atom()].value();

    //qDebug() << "atom" << atomidx << "phi" << phi_a;

    NVector delta_phi(nbonded);
    {
        int i = 0;
        foreach (AtomIdx bonded_atom, bonded_atoms)
        {
            CGAtomIdx bonded_cgatom = molinfo.cgAtomIdx(bonded_atom);

            delta_phi[i] = phi_a - moltable[bonded_cgatom.cutGroup()]
                                           [bonded_cgatom.atom()].value();

            ++i;
        }
    }

    //qDebug() << "delta_phi\n" << delta_phi.toString();

    NVector delta_q = alpha_inv_XX * delta_phi;

    //qDebug() << "delta_q\n" << delta_q.toString();

    //calculate the self energy - this is (1 / 2 alpha) p^T X X^T p
    MolarEnergy self_nrg( 0.5 * delta_q.dot( inv_alpha_XX * delta_q ) );

    //qDebug() << "self_nrg" << self_nrg.toString();

    induced_charges.set(cgatomidx,
                        induced_charges[cgatomidx] + Charge( -(delta_q.sum()) ) );

    //qDebug() << "SUM charges" << delta_q.sum() << induced_charges[cgatomidx];

    selfpol_nrgs.set(cgatomidx, self_nrg);

    {
        int i = 0;

        foreach (AtomIdx bonded_atom, bonded_atoms)
        {
            CGAtomIdx bonded_cgatom = molinfo.cgAtomIdx(bonded_atom);
            induced_charges.set( bonded_cgatom,
                                  induced_charges[bonded_cgatom] +
                                                Charge( delta_q[i] ) );
            ++i;
        }
    }
}

static QPair<AtomCharges,AtomEnergies>
calculateCharges(const MoleculeView &molview, const PolariseChargesData &poldata,
                 const MolPotentialTable &moltable)
{
    const AtomSelection &selected_atoms = poldata.selected_atoms;
    const NMatrix *xx_mat_array = poldata.xx_matricies.constData();
    const NMatrix *inv_xx_mat_array = poldata.inv_xx_matricies.constData();
    const MoleculeInfoData &molinfo = molview.data().info();
    const Connectivity &connectivity = poldata.connectivity;

    AtomCharges induced_charges(molinfo);
    AtomEnergies selfpol_nrgs(molinfo);

    if (selected_atoms.isEmpty())
    {
        //all of the atoms have been selected
        int nats = moltable.nValues();
        BOOST_ASSERT( nats == poldata.inv_xx_matricies.count() );

        for (AtomIdx i(0); i<nats; ++i)
        {
            calculateCharges(i, moltable, xx_mat_array[i],
                             inv_xx_mat_array[i], connectivity,
                             molinfo, induced_charges, selfpol_nrgs);


        }
    }
    else if (selected_atoms.selectedAllCutGroups())
    {
        for (CGIdx i(0); i<molinfo.nCutGroups(); ++i)
        {
            if (selected_atoms.selectedAll(i))
            {
                for (Index j(0); j<molinfo.nAtoms(i); ++j)
                {
                    CGAtomIdx cgatomidx(i,j);
                    AtomIdx atomidx = molinfo.atomIdx(cgatomidx);

                    calculateCharges(atomidx, moltable,
                                     xx_mat_array[atomidx],
                                     inv_xx_mat_array[atomidx],
                                     connectivity, molinfo, induced_charges,
                                     selfpol_nrgs);
                }
            }
            else
            {
                foreach (Index j, selected_atoms.selectedAtoms(i))
                {
                    CGAtomIdx cgatomidx(i,j);
                    AtomIdx atomidx = molinfo.atomIdx(cgatomidx);

                    calculateCharges(atomidx, moltable,
                                     xx_mat_array[atomidx],
                                     inv_xx_mat_array[atomidx],
                                     connectivity, molinfo, induced_charges,
                                     selfpol_nrgs);
                }
            }
        }
    }
    else
    {
        foreach (CGIdx i, selected_atoms.selectedCutGroups())
        {
            if (selected_atoms.selectedAll(i))
            {
                for (Index j(0); j<molinfo.nAtoms(i); ++j)
                {
                    CGAtomIdx cgatomidx(i,j);
                    AtomIdx atomidx = molinfo.atomIdx(cgatomidx);

                    calculateCharges(atomidx, moltable,
                                     xx_mat_array[atomidx],
                                     inv_xx_mat_array[atomidx],
                                     connectivity, molinfo, induced_charges,
                                     selfpol_nrgs);
                }
            }
            else
            {
                foreach (Index j, selected_atoms.selectedAtoms(i))
                {
                    CGAtomIdx cgatomidx(i,j);
                    AtomIdx atomidx = molinfo.atomIdx(cgatomidx);

                    calculateCharges(atomidx, moltable,
                                     xx_mat_array[atomidx],
                                     inv_xx_mat_array[atomidx],
                                     connectivity, molinfo, induced_charges,
                                     selfpol_nrgs);
                }
            }
        }
    }

    return QPair<AtomCharges,AtomEnergies>(induced_charges, selfpol_nrgs);
}

static bool haveConverged(const AtomCharges &new_charges,
                          const AtomCharges &old_charges,
                          double convergence_limit)
{
    //calculate the root mean square change in charge
    double msd = 0;

    const int nchgs = new_charges.nAtoms();
    BOOST_ASSERT( old_charges.nAtoms() == nchgs );

    const Charge *new_array = new_charges.array().constValueData();
    const Charge *old_array = old_charges.array().constValueData();

    for (int i=0; i<nchgs; ++i)
    {
        msd += SireMaths::pow_2( new_array[i].value() - old_array[i].value() );
    }

    //qDebug() << "Charge RMSD ==" << std::sqrt( msd / nchgs );

    return (std::sqrt( msd / nchgs ) < convergence_limit);
}

/** Set the baseline system for the constraint - this is
    used to pre-calculate everything for the system
    and to check if the constraint is satisfied */
void PolariseCharges::setSystem(const System &system)
{
    if (Constraint::wasLastSystem(system) and Constraint::wasLastSubVersion(system))
        return;

    Constraint::clearLastSystem();
    changed_mols = Molecules();

    //update the molecule group - clear the current list of molecules
    //if molecules have been added or removed from the group]
    {
        Version old_version = this->moleculeGroup().version();

        this->updateGroup(system);

        if (old_version.majorVersion() != this->moleculeGroup().version().majorVersion())
            moldata.clear();
    }

    const Molecules &molecules = this->moleculeGroup().molecules();

    const PropertyMap &map = this->propertyMap();
    const PropertyName coords_property = map["coordinates"];
    const PropertyName connectivity_property = map["connectivity"];
    const PropertyName polarise_property = map["polarisability"];
    const PropertyName induced_charges_property = map["induced_charge"];
    const PropertyName fixed_charges_property = map["fixed_charge"];
    const PropertyName charges_property = map["charge"];
    const PropertyName energies_property = map["self_polnrg"];

    //now calculate the potential on each molecule
    PotentialTable potentials(this->moleculeGroup());
    System new_system(system);
    new_system.potential(potentials, field_component, field_probe);

    //now calculate the induced charges
    for (Molecules::const_iterator it = molecules.constBegin();
         it != molecules.constEnd();
         ++it)
    {
        QHash< MolNum,SharedDataPointer<PolariseChargesData> >::const_iterator
                                        it2 = moldata.constFind(it.key());

        if (it2 == moldata.constEnd())
        {
            moldata.insert( it.key(), SharedDataPointer<PolariseChargesData>(
                                            new PolariseChargesData(
                                                    it.value(),
                                                    coords_property,
                                                    connectivity_property,
                                                    polarise_property) ) );
        }
        else
        {
            if (it2.value()->molversion != it.value().version())
            {
                //this molecule has changed and needs updating
                *(moldata[it.key()]) = PolariseChargesData(it.value(), coords_property,
                                                           connectivity_property,
                                                           polarise_property);
            }
        }

        it2 = moldata.constFind(it.key());

        const MolPotentialTable &moltable = potentials.getTable(it.key());

        QPair<AtomCharges,AtomEnergies> result
                                = calculateCharges(it.value(), *(it2.value().constData()),
                                                   moltable);

        const AtomCharges &new_charges = result.first;
        const AtomEnergies &selfpol_nrgs = result.second;

        if (not new_charges.isEmpty())
        {
            Molecule new_mol(it.value());

            if (new_mol.hasProperty(induced_charges_property))
            {
                const Property &p = new_mol.property(induced_charges_property);

                if (p.isA<AtomCharges>())
                {
                    //have the charges converged?
                    if ( haveConverged(new_charges, p.asA<AtomCharges>(),
                                       convergence_limit) )
                        continue;
                }
            }

            if (new_mol.hasProperty(fixed_charges_property))
            {
                PackedArray2D<Charge> charges = new_mol.property(fixed_charges_property)
                                                       .asA<AtomCharges>()
                                                       .array();

                Charge *charges_array = charges.valueData();
                const Charge *new_charges_array = new_charges.array().constValueData();

                BOOST_ASSERT( charges.nValues() == new_charges.array().nValues() );

                for (int i=0; i<charges.nValues(); ++i)
                {
                    charges_array[i] += new_charges_array[i];
                }

                new_mol = new_mol.edit()
                                 .setProperty(induced_charges_property, new_charges)
                                 .setProperty(charges_property, AtomCharges(charges))
                                 .setProperty(energies_property, selfpol_nrgs)
                                 .commit();
            }
            else
            {
                new_mol = new_mol.edit()
                                 .setProperty(induced_charges_property, new_charges)
                                 .setProperty(energies_property, selfpol_nrgs)
                                 .commit();
            }

            changed_mols.add(new_mol);
        }
    }

    Constraint::setSatisfied(system, changed_mols.isEmpty());
}

/** Return whether or not the changes in the passed
    delta *may* have changed the system since the last
    subversion 'subversion' */
bool PolariseCharges::mayChange(const Delta &delta, quint32 last_subversion) const
{
    //this constraint will need to be applied if we change any molecule
    //in the system, or if we change properties that are used to calculate
    //the energy, or if we change components used in energy expressions.
    // This is too complex to work out, so we will just say that this constraint
    // has to be applied after *any* change.
    return true;
}

/** Fully apply this constraint on the passed delta - this returns
    whether or not this constraint affects the delta */
bool PolariseCharges::fullApply(Delta &delta)
{
    this->setSystem(delta.deltaSystem());

    if (not changed_mols.isEmpty())
    {
        //qDebug() << "Iteration 1";
        bool changed = delta.update(changed_mols);

        if (changed)
        {
            if (this->moleculeGroup().nMolecules() > 1)
            {
                //we have to iterate to ensure self-consistent polarisation
                int n_iterations = 1;

                while (changed and n_iterations < 10)
                {
                    changed = false;
                    ++n_iterations;

                    //qDebug() << "Iteration" << n_iterations;
                    this->setSystem(delta.deltaSystem());

                    if (not changed_mols.isEmpty())
                        changed = delta.update(changed_mols);
                }
            }

            //qDebug() << "Complete!";

            return true;
        }
    }

    return false;
}

/** Apply this constraint based on the delta, knowing that the
    last application of this constraint was on this system,
    at subversion number last_subversion */
bool PolariseCharges::deltaApply(Delta &delta, quint32 last_subversion)
{
    return this->fullApply(delta);
}

/** Return the forcefield that is used to calculate the self-energy of
    polarising the charges. This must be added to any system to which
    this constraint is applied, as maintaining the constraint
    (by polarising the charges) costs energy, which must be part
    of the system Hamiltonian */
PolariseChargesFF PolariseCharges::selfEnergyFF() const
{
    return PolariseChargesFF(QObject::tr("self_polarise_%1")
                               .arg(this->moleculeGroup().name().value()), *this);
}

/////////////
///////////// Implementation of PolariseChargesFF
/////////////

static const RegisterMetaType<PolariseChargesFF> r_polchgff;

QDataStream &operator<<(QDataStream &ds,
                                          const PolariseChargesFF &polchgff)
{
    writeHeader(ds, r_polchgff, 1);

    SharedDataStream sds(ds);

    sds << polchgff.energy_property << static_cast<const G1FF&>(polchgff);

    return ds;
}

QDataStream &operator>>(QDataStream &ds,
                                          PolariseChargesFF &polchgff)
{
    VersionID v = readHeader(ds, r_polchgff);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> polchgff.energy_property >> static_cast<G1FF&>(polchgff);

        polchgff.molnrg = QHash<MolNum,MolarEnergy>();

        polchgff._pvt_updateName();
    }
    else
        throw version_error(v, "1", r_polchgff, CODELOC);

    return ds;
}

/** Null constructor */
PolariseChargesFF::PolariseChargesFF() : ConcreteProperty<PolariseChargesFF,G1FF>(true)
{
    this->_pvt_updateName();
}

/** Construct to calculate the self energy of the molecules affected
    by the passed constraint - note that this forcefield won't notice
    if molecules are added or removed from the constraint, so you must
    make sure that you add or remove molecules from this forcefield whenever
    you add or remove molecules from this constraint */
PolariseChargesFF::PolariseChargesFF(const PolariseCharges &constraint)
                  : ConcreteProperty<PolariseChargesFF,G1FF>(true)
{
    this->_pvt_updateName();
    energy_property = constraint.propertyMap()["self_polnrg"];
    this->add( constraint.moleculeGroup() );
}

/** Construct to calculate the self energy of the molecules affected
    by the passed constraint - note that this forcefield won't notice
    if molecules are added or removed from the constraint, so you must
    make sure that you add or remove molecules from this forcefield whenever
    you add or remove molecules from this constraint */
PolariseChargesFF::PolariseChargesFF(const QString &name,
                                     const PolariseCharges &constraint)
                  : ConcreteProperty<PolariseChargesFF,G1FF>(true)
{
    FF::setName(name);
    energy_property = constraint.propertyMap()["self_polnrg"];
    this->add( constraint.moleculeGroup() );
}

/** Copy constructor */
PolariseChargesFF::PolariseChargesFF(const PolariseChargesFF &other)
                  : ConcreteProperty<PolariseChargesFF,G1FF>(other),
                    ffcomponent(other.ffcomponent),
                    energy_property(other.energy_property),
                    molnrg(other.molnrg)
{}

/** Destructor */
PolariseChargesFF::~PolariseChargesFF()
{}

const char* PolariseChargesFF::typeName()
{
    return QMetaType::typeName( qMetaTypeId<PolariseChargesFF>() );
}

/** Copy assignment operator */
PolariseChargesFF& PolariseChargesFF::operator=(const PolariseChargesFF &other)
{
    if (this != &other)
    {
        G1FF::operator=(other);
        ffcomponent = other.ffcomponent;
        energy_property = other.energy_property;
        molnrg = other.molnrg;
    }

    return *this;
}

/** Comparison operator */
bool PolariseChargesFF::operator==(const PolariseChargesFF &other) const
{
    return this == &other or
           (energy_property == other.energy_property and
            G1FF::operator==(other));
}

/** Comparison operator */
bool PolariseChargesFF::operator!=(const PolariseChargesFF &other) const
{
    return not PolariseChargesFF::operator==(other);
}

PolariseChargesFF* PolariseChargesFF::clone() const
{
    return new PolariseChargesFF(*this);
}

/** Return the components of this forcefield */
const SingleComponent& PolariseChargesFF::components() const
{
    return ffcomponent;
}

/** You cannot set any properties of this forcefield

    \throw SireError::incompatible_error
*/
bool PolariseChargesFF::setProperty(const QString &name, const Property &property)
{
    throw SireError::incompatible_error( QObject::tr(
                "The PolariseChargesFF forcefield cannot have any properties "
                "set, so you cannot set the property \"%1\" to the "
                "value %2.")
                    .arg(name, property.toString()), CODELOC );

    return false;
}

/** This forcefield doesn't contain any properties */
bool PolariseChargesFF::containsProperty(const QString &name) const
{
    return false;
}

Q_GLOBAL_STATIC( Properties, nullProperties );

/** This forcefield doesn't contain any properties */
const Properties& PolariseChargesFF::properties() const
{
    return *(nullProperties());
}

/** This forcefield doesn't contain any properties

    \throw SireBase::missing_property
*/
const Property& PolariseChargesFF::property(const QString &name) const
{
    throw SireBase::missing_property( QObject::tr(
            "The PolariseChargesFF does not contain any properties, so "
            "it definitely doesn't contain the property \"%1\".")
                .arg(name), CODELOC );

    return Property::null();
}

/** Tell the forcefield that the energy must now be recalculated
    from scratch */
void PolariseChargesFF::mustNowRecalculateFromScratch()
{
    molnrg.clear();
}

/** Return the components of this forcefield */
const SingleComponent& PolariseChargesFF::_pvt_components() const
{
    return ffcomponent;
}

/** Calculate the energy of this forcefield */
void PolariseChargesFF::recalculateEnergy()
{
    MolarEnergy nrg(0);

    //loop over all of the molecules in this forcefield
    Molecules mols = this->molecules();

    for (Molecules::const_iterator it = mols.constBegin();
         it != mols.constEnd();
         ++it)
    {
        QHash<MolNum,MolarEnergy>::const_iterator it2 = molnrg.constFind(it.key());

        if (it2 != molnrg.constEnd())
        {
            nrg += it2.value();
        }
        else
        {
            if (not it.value().data().hasProperty(energy_property))
            {
                molnrg.insert(it.key(), MolarEnergy(0));
                continue;
            }

            //need to calculate this energy
            const AtomEnergies &selfpol = it.value().data()
                                            .property(energy_property)
                                            .asA<AtomEnergies>();

            if (it->selectedAll())
            {
                const MolarEnergy *nrgs_array = selfpol.array().constValueData();
                int nats = selfpol.nAtoms();

                MolarEnergy this_molnrg(0);

                for (int i=0; i<nats; ++i)
                {
                    this_molnrg += nrgs_array[i];
                }

                nrg += this_molnrg;
                molnrg.insert(it.key(), this_molnrg);
            }
            else
            {
                QVector<MolarEnergy> nrgs = selfpol.toVector(it->selection());

                const MolarEnergy *nrgs_array = nrgs.constData();
                int nats = nrgs.count();

                MolarEnergy this_molnrg(0);

                for (int i=0; i<nats; ++i)
                {
                    this_molnrg += nrgs_array[i];
                }

                nrg += this_molnrg;
                molnrg.insert(it.key(), this_molnrg);
            }
        }
    }

    this->components().setEnergy(*this, SingleEnergy(nrg));
}

/** Update from the changed name of this forcefield */
void PolariseChargesFF::_pvt_updateName()
{
    ffcomponent = SingleComponent(this->name());
    G1FF::_pvt_updateName();
}

/** Called when a molecule is added to this forcefield */
void PolariseChargesFF::_pvt_added(const PartialMolecule&, const PropertyMap&)
{
    this->mustNowRecalculateFromScratch();
}

/** Called when a molecule is remove from this forcefield */
void PolariseChargesFF::_pvt_removed(const PartialMolecule&)
{
    this->mustNowRecalculateFromScratch();
}

/** Called when a molecule is changed in this forcefield */
void PolariseChargesFF::_pvt_changed(const Molecule &mol, bool auto_update)
{
    molnrg.remove(mol.number());
}

/** Called when a molecule is changed in this forcefield */
void PolariseChargesFF::_pvt_changed(const QList<Molecule> &mols, bool auto_update)
{
    if (2*mols.count() > molnrg.count())
    {
        molnrg.clear();
    }
    else
    {
        for (QList<Molecule>::const_iterator it = mols.constBegin();
             it != mols.constEnd();
             ++it)
        {
            molnrg.remove(it->number());
        }
    }
}

/** Called when all molecules are removed from this forcefield */
void PolariseChargesFF::_pvt_removedAll()
{
    molnrg.clear();
}

/** Called to test if the passed property map would change the
    properties used for the passed molecule - in this case
    it wouldn't, as the properties are fixed */
bool PolariseChargesFF::_pvt_wouldChangeProperties(SireMol::MolNum,
                                                   const PropertyMap &map) const
{
    return false;
}

/** Called when a molecule is added to this forcefield */
void PolariseChargesFF::_pvt_added(const ViewsOfMol&, const PropertyMap&)
{
    this->mustNowRecalculateFromScratch();
}

/** Called when a molecule is removed from this forcefield */
void PolariseChargesFF::_pvt_removed(const ViewsOfMol&)
{
    this->mustNowRecalculateFromScratch();
}

/** Called when a molecule is removed from this forcefield */
void PolariseChargesFF::_pvt_removedAll(const PartialMolecule&)
{
    this->mustNowRecalculateFromScratch();
}

/** Called when a molecule is removed from this forcefield */
void PolariseChargesFF::_pvt_removedAll(const ViewsOfMol&)
{
    this->mustNowRecalculateFromScratch();
}
