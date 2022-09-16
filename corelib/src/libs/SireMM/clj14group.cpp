/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2014  Christopher Woods
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

#include "clj14group.h"
#include "cljnbpairs.h"

#include "SireMol/atomcoords.h"
#include "SireMol/atomcharges.h"
#include "SireMM/atomljs.h"
#include "SireMol/mover.hpp"
#include "SireMol/editor.hpp"
#include "SireMol/core.h"

#include "SireUnits/units.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMM;
using namespace SireMol;
using namespace SireBase;
using namespace SireStream;

namespace SireMM
{
    namespace detail
    {
        class CLJ14PairData
        {
        public:
            /** index of the two atoms, in sorted CGAtomIdx order */
            CGAtomIdx atom0, atom1;

            /** The reduced charges of both atoms */
            double chg0, chg1;

            /** The reduced sigma parameters of both atoms */
            double sig0, sig1;

            /** The reduced epsilon parameters for both atoms */
            double eps0, eps1;

            /** The 14 coulomb and LJ scale factors between these atoms */
            double coul14scl, lj14scl;

            CLJ14PairData()
                : atom0( CGAtomIdx::null() ), atom1( CGAtomIdx::null() ),
                  chg0(0), chg1(0), sig0(0), sig1(0),
                  eps0(0), eps1(0), coul14scl(0), lj14scl(0)
            {}

            CLJ14PairData(const CGAtomIdx &idx0, const CGAtomIdx &idx1,
                          const Charge &charge0, const Charge &charge1,
                          const LJParameter &lj0, const LJParameter &lj1,
                          const CLJScaleFactor &cljscl)
                : atom0(idx0), atom1(idx1),
                  coul14scl(cljscl.coulomb()), lj14scl(cljscl.lj())
            {
                //store the reduced charge and reduced sigma and epsilon parameters
                chg0 = charge0.value() * std::sqrt(SireUnits::one_over_four_pi_eps0);
                chg1 = charge1.value() * std::sqrt(SireUnits::one_over_four_pi_eps0);

                sig0 = std::sqrt(lj0.sigma());
                sig1 = std::sqrt(lj1.sigma());

                eps0 = std::sqrt(4.0 * lj0.epsilon());
                eps1 = std::sqrt(4.0 * lj1.epsilon());
            }

            CLJ14PairData(const CLJ14PairData &other)
                : atom0(other.atom0), atom1(other.atom1),
                  chg0(other.chg0), chg1(other.chg1),
                  sig0(other.sig0), sig1(other.sig1),
                  eps0(other.eps0), eps1(other.eps1),
                  coul14scl(other.coul14scl), lj14scl(other.lj14scl)
            {}

            ~CLJ14PairData()
            {}

            CLJ14PairData& operator=(const CLJ14PairData &other)
            {
                if (this != &other)
                {
                    atom0 = other.atom0;
                    atom1 = other.atom1;
                    chg0 = other.chg0;
                    chg1 = other.chg1;
                    sig0 = other.sig0;
                    sig1 = other.sig1;
                    eps0 = other.eps0;
                    eps1 = other.eps1;
                    coul14scl = other.coul14scl;
                    lj14scl = other.lj14scl;
                }

                return *this;
            }

            bool operator==(const CLJ14PairData &other) const
            {
                return atom0 == other.atom0 and
                       atom1 == other.atom1 and
                       chg0 == other.chg0 and
                       chg1 == other.chg1 and
                       sig0 == other.sig0 and
                       sig1 == other.sig1 and
                       eps0 == other.eps0 and
                       eps1 == other.eps1 and
                       coul14scl == other.coul14scl and
                       lj14scl == other.lj14scl;
            }
        };
    }
}

QDataStream& operator<<(QDataStream &ds, const SireMM::detail::CLJ14PairData &atom)
{
    quint32 version = 1;

    ds << version
       << atom.atom0 << atom.atom1
       << atom.chg0 << atom.chg1
       << atom.sig0 << atom.sig1
       << atom.eps0 << atom.eps1
       << atom.coul14scl << atom.lj14scl;

    return ds;
}

QDataStream& operator>>(QDataStream &ds, SireMM::detail::CLJ14PairData &atom)
{
    quint32 version;

    ds >> version;

    if (version == 1)
    {
        ds >> atom.atom0 >> atom.atom1
           >> atom.chg0 >> atom.chg1
           >> atom.sig0 >> atom.sig1
           >> atom.eps0 >> atom.eps1
           >> atom.coul14scl >> atom.lj14scl;
    }
    else
        throw version_error( QObject::tr(
                "You are trying to load an unsupported version of CLJ14PairData. "
                "You are loading version %1, but this code only supports version 1.")
                    .arg(version), CODELOC );

    return ds;
}

static const RegisterMetaType<CLJ14Group> r_group( NO_ROOT );

QDataStream &operator<<(QDataStream &ds, const CLJ14Group &group)
{
    writeHeader(ds, r_group, 3);

    SharedDataStream sds(ds);

    sds << group.mol << group.newmol << group.propmap
        << group.data_for_pair
        << qint32(group.combining_rules)
        << group.total_cnrg << group.total_ljnrg
        << group.needs_energy << group.is_strict;

    return ds;
}

QDataStream &operator>>(QDataStream &ds, CLJ14Group &group)
{
    VersionID v = readHeader(ds, r_group);

    if (v == 3)
    {
        SharedDataStream sds(ds);

        qint32 combining_rules;

        sds >> group.mol >> group.newmol >> group.propmap
            >> group.data_for_pair
            >> combining_rules
            >> group.total_cnrg >> group.total_ljnrg
            >> group.needs_energy >> group.is_strict;

        group.is_strict = true;

        group.combining_rules = CLJFunction::COMBINING_RULES(combining_rules);
    }
    else
        throw version_error(v, "3", r_group, CODELOC );

    return ds;
}

/** Constructor */
CLJ14Group::CLJ14Group()
           : combining_rules(CLJFunction::ARITHMETIC),
             total_cnrg(0), total_ljnrg(0),
             needs_energy(false), is_strict(false)
{}

/** Construct to calculate the 14-energy of the passed molecule */
CLJ14Group::CLJ14Group(const MoleculeView &molecule, const PropertyMap &map)
           : propmap(map),
             combining_rules(CLJFunction::ARITHMETIC),
             total_cnrg(0), total_ljnrg(0),
             needs_energy(false), is_strict(false)
{
    this->add(molecule);
}

/** Construct to calculate the 14-energy of the passed molecule using the
    supplied combining rules and strict mode */
CLJ14Group::CLJ14Group(const MoleculeView &molecule, CLJFunction::COMBINING_RULES rules,
                       bool strict, const PropertyMap &map)
           : propmap(map),
             combining_rules(rules),
             total_cnrg(0), total_ljnrg(0),
             needs_energy(false), is_strict(strict)
{
    this->add(molecule);
}

/** Copy constructor */
CLJ14Group::CLJ14Group(const CLJ14Group &other)
           : mol(other.mol), newmol(other.newmol), propmap(other.propmap),
             data_for_pair(other.data_for_pair),
             combining_rules(other.combining_rules),
             total_cnrg(other.total_cnrg), total_ljnrg(other.total_ljnrg),
             needs_energy(other.needs_energy), is_strict(other.is_strict)
{}

/** Destructor */
CLJ14Group::~CLJ14Group()
{}

/** Copy assignment operator */
CLJ14Group& CLJ14Group::operator=(const CLJ14Group &other)
{
    if (this != &other)
    {
        mol = other.mol;
        newmol = other.newmol;
        propmap = other.propmap;
        data_for_pair = other.data_for_pair;
        combining_rules = other.combining_rules;
        total_cnrg = other.total_cnrg;
        total_ljnrg = other.total_ljnrg;
        needs_energy = other.needs_energy;
        is_strict = other.is_strict;
    }

    return *this;
}

/** Comparison operator */
bool CLJ14Group::operator==(const CLJ14Group &other) const
{
    return mol == other.mol and
           newmol == other.newmol and
           combining_rules == other.combining_rules and
           propmap == other.propmap and
           is_strict == other.is_strict;
}

/** Comparison operator */
bool CLJ14Group::operator!=(const CLJ14Group &other) const
{
    return not operator==(other);
}

const char* CLJ14Group::typeName()
{
    return QMetaType::typeName( qMetaTypeId<CLJ14Group>() );
}

const char* CLJ14Group::what() const
{
    return CLJ14Group::typeName();
}

bool CLJ14Group::isNull() const
{
    return newmol.isNull();
}

QString CLJ14Group::toString() const
{
    if (isNull())
        return QObject::tr("CLJ14Group::null");
    else
        return QObject::tr("CLJ14Group( molecule = %1 )").arg(newmol.toString());
}

/** Return the molecule that is in this group */
const MoleculeView& CLJ14Group::molecule() const
{
    return newmol;
}

/** Return the property map used to find the properties needed to
    calculate the 14 energy ("coordinates", "charge", "LJ"
    and "intrascale") */
PropertyMap CLJ14Group::propertyMap() const
{
    return propmap;
}

/** Return whether or not we are recalculating things from scratch */
bool CLJ14Group::recalculatingFromScratch() const
{
    return (mol.version() == newmol.version()) and needs_energy;
}

/** Set the flag to say that we must recalculate everything
    from scratch */
void CLJ14Group::mustNowRecalculateFromScratch()
{
    mol = newmol;
    total_cnrg = 0;
    total_ljnrg = 0;
    needs_energy = true;
}

/** Set the flag to ensure that the energy is really completely
    recalculated from scratch */
void CLJ14Group::mustReallyRecalculateFromScratch()
{
    mol = newmol;
    total_cnrg = 0;
    total_ljnrg = 0;
    data_for_pair.clear();
    needs_energy = true;
}

/** Set whether or not 'strict' mode is on. If 'strict' mode is on,
    then this means that the 1-4 energy is calculated only if both of the
    atoms are selected. If 'strict' mode is off, then the 1-4 energy
    is calculated when at least one of the atoms is selected. */
bool CLJ14Group::setStrict(bool isstrict)
{
    if (not isstrict == is_strict)
    {
        is_strict = isstrict;
        this->mustReallyRecalculateFromScratch();
        return true;
    }
    else
        return false;
}

/** Return whether or not 'strict' mode is on */
bool CLJ14Group::isStrict() const
{
    return is_strict;
}

/** Switch on or off the use of arithmetic combining rules */
void CLJ14Group::setArithmeticCombiningRules(bool on)
{
    if (on)
        setCombiningRules(CLJFunction::ARITHMETIC);
    else
        setCombiningRules(CLJFunction::GEOMETRIC);
}

/** Switch on or off the use og geometric combining rules */
void CLJ14Group::setGeometricCombiningRules(bool on)
{
    if (on)
        setCombiningRules(CLJFunction::GEOMETRIC);
    else
        setCombiningRules(CLJFunction::ARITHMETIC);
}

/** Return the type of combining rules in place */
CLJFunction::COMBINING_RULES CLJ14Group::combiningRules() const
{
    return combining_rules;
}

/** Set the combining rules to 'rules' */
void CLJ14Group::setCombiningRules(CLJFunction::COMBINING_RULES rules)
{
    if (combining_rules != rules)
    {
        combining_rules = rules;
        this->mustNowRecalculateFromScratch();
    }
}

/** Return whether or not arithmetic combining rules are used */
bool CLJ14Group::usingArithmeticCombiningRules() const
{
    return combining_rules == CLJFunction::ARITHMETIC;
}

/** Return whether or not geometric combining rules are used */
bool CLJ14Group::usingGeometricCombiningRules() const
{
    return combining_rules == CLJFunction::GEOMETRIC;
}

/** Add the passed molecule to this group */
void CLJ14Group::add(const MoleculeView &new_molecule)
{
    if (isNull())
    {
        mustReallyRecalculateFromScratch();
        newmol = new_molecule;
    }
    else
    {
        mustReallyRecalculateFromScratch();
        AtomSelection selected_atoms = newmol.selection();
        selected_atoms = selected_atoms.select(new_molecule.selection());
        this->updateSelection(selected_atoms);
        this->update(new_molecule);
    }
}

/** Add the passed selection onto this group */
void CLJ14Group::add(const AtomSelection &new_selection)
{
    AtomSelection selected_atoms = newmol.selection();
    selected_atoms = selected_atoms.select(new_selection);
    this->updateSelection(selected_atoms);
}

/** Update the selection to the new passed value */
void CLJ14Group::updateSelection(const AtomSelection &selection)
{
    if (isNull())
        return;

    else if (selection.selectedNone())
    {
        PropertyMap map = propmap;
        this->operator=( CLJ14Group() );
        propmap = map;
        return;
    }

    if (newmol.selection() != selection)
    {
        newmol = PartialMolecule(newmol.data(), selection);
        mustReallyRecalculateFromScratch();
    }
}

/** Update the contained molecule to the newest version */
void CLJ14Group::update(const MoleculeView &new_molecule)
{
    if (isNull())
        return;

    if (new_molecule.data().number() != newmol.number())
        throw SireError::incompatible_error( QObject::tr(
                    "You cannot update the new molecule %1 as it is not the same "
                    "molecule as the old molecule %2.")
                        .arg(new_molecule.toString())
                        .arg(mol.toString()), CODELOC );

    if (recalculatingFromScratch())
    {
        newmol = PartialMolecule(new_molecule.data(), newmol.selection());
        mol = new_molecule;
        needs_energy = true;
    }
    else
    {
        newmol = PartialMolecule(new_molecule.data(), newmol.selection());
        needs_energy = (newmol.version() != mol.version());
    }
}

/** Remove the passed selection from the group */
void CLJ14Group::remove(const AtomSelection &new_selection)
{
    if (not isNull())
    {
        AtomSelection selected_atoms = newmol.selection();
        selected_atoms = selected_atoms.deselect(new_selection);
        this->updateSelection(selected_atoms);
    }
}

/** Remove the passed molecule from this group */
void CLJ14Group::remove(const MoleculeView &new_molecule)
{
    if (not isNull())
    {
        this->remove(new_molecule.selection());
        this->update(new_molecule);
    }
}

/** Return whether or not the passed property map would change properties that
    are used by this calculation */
bool CLJ14Group::wouldChangeProperties(const PropertyMap &map) const
{
    if (propmap != map)
    {
        return (map["coordinates"] != propmap["coordinates"]) or
               (map["charge"] != propmap["charge"]) or
               (map["LJ"] != propmap["LJ"]) or
               (map["intrascale"] != propmap["intrascale"]);
    }
    else
        return false;
}

/** Internal function used to extract all of the 1-4 pairs */
void CLJ14Group::reextract()
{
    data_for_pair.clear();
    total_cnrg = 0;
    total_ljnrg = 0;
    needs_energy = true;

    const MoleculeInfoData &molinfo = newmol.data().info();
    const AtomSelection selected_atoms = newmol.selection();

    if (selected_atoms.selectedNone())
    {
        //don't have enough atoms for 1-4 pairs!
        needs_energy = false;
        return;
    }

    const CLJNBPairs &pairs = newmol.data().property( propmap["intrascale"] ).asA<CLJNBPairs>();

    const AtomCharges &chgs = newmol.data().property( propmap["charge"] ).asA<AtomCharges>();
    const AtomLJs &ljs = newmol.data().property( propmap["LJ"] ).asA<AtomLJs>();

    const auto &connectivity = newmol.data().property( propmap["connectivity"] ).asA<Connectivity>();

    const auto bond_matrix = connectivity.getBondMatrix(4,4);

    const int nats = bond_matrix.count();

    for (int i=0; i<nats-1; ++i)
    {
        const auto row = bond_matrix.constData()[i];

        const AtomIdx atomidx0(i);
        const auto cgatomidx0 = molinfo.cgAtomIdx(atomidx0);

        if (selected_atoms.selectedAll() or selected_atoms.selected(cgatomidx0))
        {
            for (int j=i+1; j<nats; ++j)
            {
                if (row.constData()[j])
                {
                    const AtomIdx atomidx1(j);
                    const auto cgatomidx1 = molinfo.cgAtomIdx(atomidx1);

                    //these are 1-4 bonded :-)
                    if ((not is_strict) or (selected_atoms.selectedAll()) or (selected_atoms.selected(cgatomidx1)))
                    {
                        const auto scl = pairs.get(cgatomidx0, cgatomidx1);

                        if (scl != CLJScaleFactor(0))
                        {
                            data_for_pair.append( SireMM::detail::CLJ14PairData(
                                cgatomidx0, cgatomidx1,
                                chgs[cgatomidx0], chgs[cgatomidx1],
                                ljs[cgatomidx0], ljs[cgatomidx1],
                                scl) );
                        }
                    }
                }
            }
        }
    }

    data_for_pair.squeeze();
}

/** Internal function used to calculate the energy of the 1-4 interactions for the
    atom pairs in 'pairdata' and to return it in the passed 'cnrg' and 'ljnrg' variables
    (using arithmetic combining rules) */
static void calculateAriEnergy(const AtomCoords &coords,
                               const QVector<SireMM::detail::CLJ14PairData> &pairdata,
                               double &cnrg, double &ljnrg)
{
    cnrg = 0;
    ljnrg = 0;

    for (int i=0; i<pairdata.count(); ++i)
    {
        const SireMM::detail::CLJ14PairData &pair = pairdata.constData()[i];

        const Vector &atom0 = coords[pair.atom0];
        const Vector &atom1 = coords[pair.atom1];

        const double one_over_r = Vector::invDistance(atom0, atom1);

        cnrg += (pair.coul14scl * pair.chg0 * pair.chg1 * one_over_r);

        //arithmetic combining rules
        const double sig = 0.5 * ( (pair.sig0*pair.sig0) + (pair.sig1*pair.sig1) );

        const double sig2_over_r2 = SireMaths::pow_2( sig * one_over_r );
        const double sig6_over_r6 = SireMaths::pow_3( sig2_over_r2 );
        const double sig12_over_r12 = SireMaths::pow_2( sig6_over_r6 );

        ljnrg += (pair.lj14scl * pair.eps0 * pair.eps1) * (sig12_over_r12 - sig6_over_r6);
    }
}

/** Internal function used to calculate the energy of the 1-4 interactions for the
    atom pairs in 'pairdata' and to return it in the passed 'cnrg' and 'ljnrg' variables
    (using geometric combining rules) */
static void calculateGeoEnergy(const AtomCoords &coords,
                               const QVector<SireMM::detail::CLJ14PairData> &pairdata,
                               double &cnrg, double &ljnrg)
{
    cnrg = 0;
    ljnrg = 0;

    for (int i=0; i<pairdata.count(); ++i)
    {
        const SireMM::detail::CLJ14PairData &pair = pairdata.constData()[i];

        const Vector &atom0 = coords[pair.atom0];
        const Vector &atom1 = coords[pair.atom1];

        const float one_over_r = Vector::invDistance(atom0, atom1);

        cnrg += (pair.coul14scl * pair.chg0 * pair.chg1 * one_over_r);

        //geometric combining rules
        const float sig = pair.sig0*pair.sig1;

        const float sig2_over_r2 = SireMaths::pow_2( sig * one_over_r );
        const float sig6_over_r6 = SireMaths::pow_3( sig2_over_r2 );
        const float sig12_over_r12 = SireMaths::pow_2( sig6_over_r6 );

        ljnrg += (pair.lj14scl * pair.eps0 * pair.eps1) * (sig12_over_r12 - sig6_over_r6);
    }

    ljnrg *= 4.0;
}

/** Internal function used to calculate the energy of the 1-4 interactions for the
    atom pairs in 'pairdata' and to return it in the passed 'cnrg' and 'ljnrg' variables
    (using the specified combining rules) */
static void calculateEnergy(CLJFunction::COMBINING_RULES combining_rules,
                            const AtomCoords &coords,
                            const QVector<SireMM::detail::CLJ14PairData> &pairdata,
                            double &cnrg, double &ljnrg)
{
    switch (combining_rules)
    {
    case CLJFunction::ARITHMETIC:
        calculateAriEnergy(coords, pairdata, cnrg, ljnrg);
        return;
    case CLJFunction::GEOMETRIC:
        calculateGeoEnergy(coords, pairdata, cnrg, ljnrg);
        return;
    default:
        throw SireError::unsupported( QObject::tr(
                "Unknown combining rules requested! %1")
                    .arg(combining_rules), CODELOC );
    }
}

/** Calculate and return the coulomb and LJ 14 energy */
boost::tuple<double,double> CLJ14Group::energy()
{
    if (not needs_energy)
    {
        return boost::tuple<double,double>(total_cnrg,total_ljnrg);
    }

    if (data_for_pair.isEmpty())
    {
        this->reextract();
    }

    total_cnrg = 0;
    total_ljnrg = 0;

    ::calculateEnergy(combining_rules,
                      newmol.data().property(propmap["coordinates"])
                                   .asA<AtomCoords>(),
                      data_for_pair, total_cnrg, total_ljnrg);

    mol = newmol;
    needs_energy = false;

    return boost::tuple<double,double>(total_cnrg,total_ljnrg);
}
