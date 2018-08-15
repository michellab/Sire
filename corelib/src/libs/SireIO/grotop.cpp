/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2017  Christopher Woods
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

#include "SireIO/grotop.h"

#include "SireSystem/system.h"

#include "SireError/errors.h"
#include "SireIO/errors.h"

#include "SireMol/molecule.h"
#include "SireMol/moleditor.h"
#include "SireMol/moleculegroup.h"
#include "SireMol/errors.h"
#include "SireMol/atomcharges.h"
#include "SireMol/atommasses.h"
#include "SireMol/atomelements.h"
#include "SireMol/connectivity.h"

#include "SireMM/internalff.h"
#include "SireMM/atomljs.h"
#include "SireMM/twoatomfunctions.h"
#include "SireMM/threeatomfunctions.h"
#include "SireMM/fouratomfunctions.h"
#include "SireMM/cljnbpairs.h"

#include "SireBase/parallel.h"
#include "SireBase/stringproperty.h"
#include "SireBase/booleanproperty.h"

#include "SireUnits/units.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QRegularExpression>
#include <QFileInfo>
#include <QDateTime>
#include <QDir>
#include <QElapsedTimer>

using namespace SireIO;
using namespace SireUnits;
using namespace SireMol;
using namespace SireMM;
using namespace SireFF;
using namespace SireBase;
using namespace SireSystem;
using namespace SireStream;

////////////////
//////////////// Implementation of GroAtom
////////////////

static const RegisterMetaType<GroAtom> r_groatom(NO_ROOT);

QDataStream SIREIO_EXPORT &operator<<(QDataStream &ds, const GroAtom &atom)
{
    writeHeader(ds, r_groatom, 2);

    SharedDataStream sds(ds);

    sds << atom.atmname << atom.resname << atom.atmtyp << atom.bndtyp
        << atom.atmnum << atom.resnum << atom.chggrp
        << atom.chg.to(mod_electron) << atom.mss.to(g_per_mol);

    return ds;
}

QDataStream SIREIO_EXPORT &operator>>(QDataStream &ds, GroAtom &atom)
{
    VersionID v = readHeader(ds, r_groatom);

    if (v == 2)
    {
        SharedDataStream sds(ds);

        double chg, mass;

        sds >> atom.atmname >> atom.resname >> atom.atmtyp >> atom.bndtyp
            >> atom.atmnum >> atom.resnum >> atom.chggrp
            >> chg >> mass;

        atom.chg = chg*mod_electron;
        atom.mss = mass*g_per_mol;
    }
    else if (v == 1)
    {
        SharedDataStream sds(ds);

        double chg, mass;

        sds >> atom.atmname >> atom.resname >> atom.atmtyp
            >> atom.atmnum >> atom.resnum >> atom.chggrp
            >> chg >> mass;

        atom.bndtyp = atom.atmtyp;
        atom.chg = chg*mod_electron;
        atom.mss = mass*g_per_mol;
    }
    else
        throw version_error(v, "1,2", r_groatom, CODELOC);

    return ds;
}

/** Constructor */
GroAtom::GroAtom() : atmnum(-1), resnum(-1), chggrp(-1), chg(0), mss(0)
{}

/** Copy constructor */
GroAtom::GroAtom(const GroAtom &other)
        : atmname(other.atmname), resname(other.resname),
          atmtyp(other.atmtyp), bndtyp(other.bndtyp), atmnum(other.atmnum),
          resnum(other.resnum), chggrp(other.chggrp),
          chg(other.chg), mss(other.mss)
{}

/** Destructor */
GroAtom::~GroAtom()
{}

/** Copy assignment operator */
GroAtom& GroAtom::operator=(const GroAtom &other)
{
    if (this != &other)
    {
        atmname = other.atmname;
        resname = other.resname;
        atmtyp = other.atmtyp;
        bndtyp = other.bndtyp;
        atmnum = other.atmnum;
        resnum = other.resnum;
        chggrp = other.chggrp;
        chg = other.chg;
        mss = other.mss;
    }

    return *this;
}

/** Comparison operator */
bool GroAtom::operator==(const GroAtom &other) const
{
    return atmname == other.atmname and
           resname == other.resname and
           atmtyp == other.atmtyp and
           bndtyp == other.bndtyp and
           atmnum == other.atmnum and
           resnum == other.resnum and
           chggrp == other.chggrp and
           chg == other.chg and
           mss == other.mss;
}

/** Comparison operator */
bool GroAtom::operator!=(const GroAtom &other) const
{
    return not operator==(other);
}

const char* GroAtom::typeName()
{
    return QMetaType::typeName( qMetaTypeId<GroAtom>() );
}

const char* GroAtom::what() const
{
    return GroAtom::typeName();
}

QString GroAtom::toString() const
{
    if (isNull())
        return QObject::tr( "GroAtom::null");
    else
        return QObject::tr("GroAtom( name() = %1, number() = %2 )")
                    .arg(atmname).arg(atmnum);
}

/** Return whether or not this atom is null */
bool GroAtom::isNull() const
{
    return operator==( GroAtom() );
}

/** Return the name of the atom */
AtomName GroAtom::name() const
{
    return AtomName(atmname);
}

/** Return the number of the atom */
AtomNum GroAtom::number() const
{
    return AtomNum(atmnum);
}

/** Return the name of the residue that contains this atom */
ResName GroAtom::residueName() const
{
    return ResName(resname);
}

/** Return the number of the residue that contains this atom */
ResNum GroAtom::residueNumber() const
{
    return ResNum(resnum);
}

/** Return the charge group of this atom */
qint64 GroAtom::chargeGroup() const
{
    return chggrp;
}

/** Return the atom type of this atom */
QString GroAtom::atomType() const
{
    return atmtyp;
}

/** Return the bond type of this atom. This is normally the same as the atom type */
QString GroAtom::bondType() const
{
    return bndtyp;
}

/** Return the charge on this atom */
SireUnits::Dimension::Charge GroAtom::charge() const
{
    return chg;
}

/** Return the mass of this atom */
SireUnits::Dimension::MolarMass GroAtom::mass() const
{
    return mss;
}

/** Set the name of this atom */
void GroAtom::setName(const QString &name)
{
    atmname = name;
}

/** Set the number of this atom */
void GroAtom::setNumber(qint64 number)
{
    if (number >= 0)
        atmnum = number;
}

/** Set the name of the residue containing this atom */
void GroAtom::setResidueName(const QString &name)
{
    resname = name;
}

/** Set the number of the residue containing this atom */
void GroAtom::setResidueNumber(qint64 number)
{
    resnum = number;
}

/** Set the charge group of this atom */
void GroAtom::setChargeGroup(qint64 grp)
{
    if (grp >= 0)
        chggrp = grp;
}

/** Set the atom type and bond type of this atom. To set
    the bond type separately, you need to set it after calling
    this function */
void GroAtom::setAtomType(const QString &atomtype)
{
    atmtyp = atomtype;
    bndtyp = atomtype;
}

/** Set the bond type of this atom */
void GroAtom::setBondType(const QString &bondtype)
{
    bndtyp = bondtype;
}

/** Set the charge on this atom */
void GroAtom::setCharge(SireUnits::Dimension::Charge charge)
{
    chg = charge;
}

/** Set the mass of this atom */
void GroAtom::setMass(SireUnits::Dimension::MolarMass mass)
{
    if (mass.value() > 0)
        mss = mass;
}

////////////////
//////////////// Implementation of GroMolType
////////////////

static const RegisterMetaType<GroMolType> r_gromoltyp(NO_ROOT);

QDataStream SIREIO_EXPORT &operator<<(QDataStream &ds, const GroMolType &moltyp)
{
    writeHeader(ds, r_gromoltyp, 2);

    SharedDataStream sds(ds);

    sds << moltyp.nme << moltyp.warns << moltyp.atms0 << moltyp.atms1
        << moltyp.bnds0 << moltyp.bnds1 << moltyp.angs0 << moltyp.angs1
        << moltyp.dihs0 << moltyp.dihs1 << moltyp.first_atoms0 << moltyp.first_atoms1
        << moltyp.ffield0 << moltyp.ffield1 << moltyp.nexcl0 << moltyp.nexcl1
        << moltyp.is_perturbable;

    return ds;
}

QDataStream SIREIO_EXPORT &operator>>(QDataStream &ds, GroMolType &moltyp)
{
    VersionID v = readHeader(ds, r_gromoltyp);

    if (v == 1 or v == 2)
    {
        SharedDataStream sds(ds);

        sds >> moltyp.nme >> moltyp.warns >> moltyp.atms0 >> moltyp.atms1
            >> moltyp.bnds0 >> moltyp.bnds1 >> moltyp.angs0 >> moltyp.angs1
            >> moltyp.dihs0 >> moltyp.dihs1 >> moltyp.first_atoms0 >> moltyp.first_atoms1;

        if (v == 2)
            sds >> moltyp.ffield0 >> moltyp.ffield1;
        else
        {
            moltyp.ffield0 = MMDetail();
            moltyp.ffield1 = MMDetail();
        }

        sds >> moltyp.nexcl0 >> moltyp.nexcl1 >> moltyp.is_perturbable;
    }
    else
        throw version_error(v, "1,2", r_gromoltyp, CODELOC);

    return ds;
}

/** Constructor */
GroMolType::GroMolType() : nexcl0(3), nexcl1(3), // default to 3 as this is normal for most molecules
                           is_perturbable(false) // default to a non-perturbable molecule
{}

/** Construct from the passed molecule */
GroMolType::GroMolType(const SireMol::Molecule &mol, const PropertyMap &map)
           : nexcl0(3), nexcl1(3),  // default to '3' as this is normal for most molecules
             is_perturbable(false)  // default to a non-perturbable molecule
{
    if (mol.nAtoms() == 0)
        return;

    // Try to see if this molecule is perturbable.
    try
    {
        is_perturbable = mol.property("is_perturbable").asABoolean();
    }
    catch (...)
    {}

    // Perturbable molecule.
    if (is_perturbable)
    {
        // For perturbable molecules we don't user the user PropertyMap to extract
        // properties since the naming must be consistent for the properties at
        // lambda = 0 and lambda = 1, e.g. "charge0" and "charge1".

        //get the name either from the molecule name or the name of the first
        //residue
        nme = mol.name();

        if (nme.isEmpty())
        {
            nme = mol.residue(ResIdx(0)).name();
        }

        //get the forcefields for this molecule
        try
        {
            ffield0 = mol.property("forcefield0").asA<MMDetail>();
        }
        catch(...)
        {
            warns.append( QObject::tr("Cannot find a valid MM forcefield for this molecule at lambda = 0!") );
        }
        try
        {
            ffield1 = mol.property("forcefield1").asA<MMDetail>();
        }
        catch(...)
        {
            warns.append( QObject::tr("Cannot find a valid MM forcefield for this molecule at lambda = 1!") );
        }

        const auto molinfo = mol.info();

        bool uses_parallel = true;
        if (map["parallel"].hasValue())
        {
            uses_parallel = map["parallel"].value().asA<BooleanProperty>().value();
        }

        //get information about all atoms in this molecule
        auto extract_atoms = [&](bool is_lambda1)
        {
            if (is_lambda1)
                atms1 = QVector<GroAtom>(molinfo.nAtoms());
            else
                atms0 = QVector<GroAtom>(molinfo.nAtoms());

            AtomMasses masses;
            AtomElements elements;
            AtomCharges charges;
            AtomIntProperty groups;
            AtomStringProperty atomtypes;
            AtomStringProperty bondtypes;

            bool has_mass(false), has_elem(false), has_chg(false), has_group(false), has_type(false);
            bool has_bondtype(false);

            try
            {
                if (is_lambda1)
                    masses = mol.property("mass1").asA<AtomMasses>();
                else
                    masses = mol.property("mass0").asA<AtomMasses>();
                has_mass = true;
            }
            catch(...)
            {}

            if (not has_mass)
            {
                try
                {
                    if (is_lambda1)
                        elements = mol.property("element1").asA<AtomElements>();
                    else
                        elements = mol.property("element0").asA<AtomElements>();
                    has_elem = true;
                }
                catch(...)
                {}
            }

            try
            {
                if (is_lambda1)
                    charges = mol.property("charge1").asA<AtomCharges>();
                else
                    charges = mol.property("charge0").asA<AtomCharges>();
                has_chg = true;
            }
            catch(...)
            {}

            try
            {
                if (is_lambda1)
                    groups = mol.property("charge_group1").asA<AtomIntProperty>();
                else
                    groups = mol.property("charge_group0").asA<AtomIntProperty>();
                has_group = true;
            }
            catch(...)
            {}

            try
            {
                if (is_lambda1)
                    atomtypes = mol.property("atomtype1").asA<AtomStringProperty>();
                else
                    atomtypes = mol.property("atomtype0").asA<AtomStringProperty>();
                has_type = true;
            }
            catch(...)
            {}

            try
            {
                if (is_lambda1)
                    bondtypes = mol.property("bondtype1").asA<AtomStringProperty>();
                else
                    bondtypes = mol.property("bondtype0").asA<AtomStringProperty>();
                has_bondtype = true;
            }
            catch(...)
            {}

            if (not (has_chg and has_type and (has_elem or has_mass)))
            {
                warns.append( QObject::tr("Cannot find valid charge, atomtype and (element or mass) "
                "properties for the molecule. These are needed! "
                "has_charge=%1, has_atomtype=%2, has_mass=%3, has_element=%4")
                    .arg(has_chg).arg(has_type).arg(has_mass).arg(has_elem) );
                return;
            }

            //run through the atoms in AtomIdx order
            auto extract_atom = [&](int iatm, bool is_lambda1)
            {
                AtomIdx i(iatm);

                const auto cgatomidx = molinfo.cgAtomIdx(i);
                const auto residx = molinfo.parentResidue(i);

                //atom numbers have to count up sequentially from 1
                int atomnum = i+1;
                QString atomnam = molinfo.name(i);

                //assuming that residues are in the same order as the atoms
                int resnum = residx + 1;
                QString resnam = molinfo.name(residx);

                //people like to preserve the residue numbers of ligands and
                //proteins. This is very challenging for the gromacs topology,
                //as it would force a different topology for every solvent molecule,
                //so deciding on the difference between protein/ligand and solvent
                //is tough. Will preserve the residue number if the number of
                //residues is greater than 1 and the number of atoms is greater
                //than 32 (so octanol is a solvent)
                if (molinfo.nResidues() > 1 or molinfo.nAtoms() > 32)
                {
                    resnum = molinfo.number(residx).value();
                }

                int group = atomnum;

                if (has_group)
                {
                    group = groups[cgatomidx];
                }

                auto charge = charges[cgatomidx];

                SireUnits::Dimension::MolarMass mass;

                if (has_mass)
                {
                    mass = masses[cgatomidx];
                }
                else
                {
                    mass = elements[cgatomidx].mass();
                }

                if (mass <= 0)
                {
                    //not allowed to have a zero or negative mass
                    mass = 1.0 * g_per_mol;
                }

                QString atomtype = atomtypes[cgatomidx];

                if (is_lambda1)
                {
                    auto &atom = atms1[i];

                    atom.setName(atomnam);
                    atom.setNumber(atomnum);
                    atom.setResidueName(resnam);
                    atom.setResidueNumber(resnum);
                    atom.setChargeGroup(group);
                    atom.setCharge(charge);
                    atom.setMass(mass);
                    atom.setAtomType(atomtype);

                    if (has_bondtype)
                    {
                        atom.setBondType(bondtypes[cgatomidx]);
                    }
                    else
                    {
                        atom.setBondType(atomtype);
                    }
                }
                else
                {
                    auto &atom = atms0[i];

                    atom.setName(atomnam);
                    atom.setNumber(atomnum);
                    atom.setResidueName(resnam);
                    atom.setResidueNumber(resnum);
                    atom.setChargeGroup(group);
                    atom.setCharge(charge);
                    atom.setMass(mass);
                    atom.setAtomType(atomtype);

                    if (has_bondtype)
                    {
                        atom.setBondType(bondtypes[cgatomidx]);
                    }
                    else
                    {
                        atom.setBondType(atomtype);
                    }
                }
            };

            if (uses_parallel)
            {
                tbb::parallel_for( tbb::blocked_range<int>(0,molinfo.nAtoms()),
                                [&](const tbb::blocked_range<int> &r)
                {
                    for (int i=r.begin(); i<r.end(); ++i)
                    {
                        extract_atom(i, is_lambda1);
                    }
                });
            }
            else
            {
                for (int i=0; i<molinfo.nAtoms(); ++i)
                {
                    extract_atom(i, is_lambda1);
                }
            }
        };

        //get all of the bonds in this molecule
        auto extract_bonds = [&](bool is_lambda1)
        {
            bool has_conn(false), has_funcs(false);

            TwoAtomFunctions funcs;
            Connectivity conn;

            const auto R = InternalPotential::symbols().bond().r();

            try
            {
                if (is_lambda1)
                    funcs = mol.property("bond1").asA<TwoAtomFunctions>();
                else
                    funcs = mol.property("bond0").asA<TwoAtomFunctions>();
                has_funcs = true;
            }
            catch(...)
            {}

            try
            {
                conn = mol.property("connectivity").asA<Connectivity>();
                has_conn = true;
            }
            catch(...)
            {}

            //get the bond potentials first
            if (has_funcs)
            {
                for (const auto bond : funcs.potentials())
                {
                    AtomIdx atom0 = molinfo.atomIdx(bond.atom0());
                    AtomIdx atom1 = molinfo.atomIdx(bond.atom1());

                    if (atom0 > atom1)
                        qSwap(atom0,atom1);

                    if (is_lambda1)
                        bnds1.insert( BondID(atom0,atom1), GromacsBond(bond.function(), R) );
                    else
                        bnds0.insert( BondID(atom0,atom1), GromacsBond(bond.function(), R) );
                }
            }

            //now fill in any missing bonded atoms with null bonds
            if (has_conn)
            {
                for (const auto bond : conn.getBonds())
                {
                    AtomIdx atom0 = molinfo.atomIdx(bond.atom0());
                    AtomIdx atom1 = molinfo.atomIdx(bond.atom1());

                    if (atom0 > atom1)
                        qSwap(atom0,atom1);

                    BondID b(atom0,atom1);

                    if (is_lambda1)
                    {
                        if (not bnds1.contains(b))
                            bnds1.insert(b, GromacsBond(5));  // function 5 is a simple connection
                    }
                    else
                    {
                        if (not bnds0.contains(b))
                            bnds0.insert(b, GromacsBond(5));  // function 5 is a simple connection
                    }
                }
            }
        };

        //get all of the angles in this molecule
        auto extract_angles = [&](bool is_lambda1)
        {
            bool has_funcs(false);

            ThreeAtomFunctions funcs;

            const auto theta = InternalPotential::symbols().angle().theta();

            try
            {
                if (is_lambda1)
                    funcs = mol.property("angle1").asA<ThreeAtomFunctions>();
                else
                    funcs = mol.property("angle0").asA<ThreeAtomFunctions>();
                has_funcs = true;
            }
            catch(...)
            {}

            if (has_funcs)
            {
                for (const auto angle : funcs.potentials())
                {
                    AtomIdx atom0 = molinfo.atomIdx(angle.atom0());
                    AtomIdx atom1 = molinfo.atomIdx(angle.atom1());
                    AtomIdx atom2 = molinfo.atomIdx(angle.atom2());

                    if (atom0 > atom2)
                        qSwap(atom0,atom2);

                    if (is_lambda1)
                    {
                        angs1.insert( AngleID(atom0,atom1,atom2),
                                    GromacsAngle(angle.function(), theta) );
                    }
                    else
                    {
                        angs0.insert( AngleID(atom0,atom1,atom2),
                                    GromacsAngle(angle.function(), theta) );
                    }
                }
            }
        };

        //get all of the dihedrals in this molecule
        auto extract_dihedrals = [&](bool is_lambda1)
        {
            bool has_funcs(false);

            FourAtomFunctions funcs;

            const auto phi = InternalPotential::symbols().dihedral().phi();

            try
            {
                if (is_lambda1)
                    funcs = mol.property("dihedral1").asA<FourAtomFunctions>();
                else
                    funcs = mol.property("dihedral0").asA<FourAtomFunctions>();
                has_funcs = true;
            }
            catch(...)
            {}

            if (has_funcs)
            {
                for (const auto dihedral : funcs.potentials())
                {
                    AtomIdx atom0 = molinfo.atomIdx(dihedral.atom0());
                    AtomIdx atom1 = molinfo.atomIdx(dihedral.atom1());
                    AtomIdx atom2 = molinfo.atomIdx(dihedral.atom2());
                    AtomIdx atom3 = molinfo.atomIdx(dihedral.atom3());

                    if (atom0 > atom3)
                    {
                        qSwap(atom0,atom3);
                        qSwap(atom1,atom2);
                    }

                    //get all of the dihedral terms (could be a lot)
                    auto parts = GromacsDihedral::construct(dihedral.function(), phi);

                    DihedralID dihid(atom0,atom1,atom2,atom3);

                    for (const auto part : parts)
                    {
                        if (is_lambda1)
                            dihs1.insertMulti(dihid, part);
                        else
                            dihs0.insertMulti(dihid, part);
                    }
                }
            }

            bool has_imps(false);

            FourAtomFunctions imps;

            try
            {
                if (is_lambda1)
                    imps = mol.property("improper1").asA<FourAtomFunctions>();
                else
                    imps = mol.property("improper0").asA<FourAtomFunctions>();
                has_imps = true;
            }
            catch(...)
            {}

            if (has_imps)
            {
                for (const auto improper : imps.potentials())
                {
                    AtomIdx atom0 = molinfo.atomIdx(improper.atom0());
                    AtomIdx atom1 = molinfo.atomIdx(improper.atom1());
                    AtomIdx atom2 = molinfo.atomIdx(improper.atom2());
                    AtomIdx atom3 = molinfo.atomIdx(improper.atom3());

                    //get all of the dihedral terms (could be a lot)
                    auto parts = GromacsDihedral::constructImproper(improper.function(), phi);

                    DihedralID impid(atom0,atom1,atom2,atom3);

                    for (const auto part : parts)
                    {
                        if (is_lambda1)
                            dihs1.insertMulti(impid, part);
                        else
                            dihs0.insertMulti(impid, part);
                    }
                }
            }
        };

        const QVector< std::function< void(bool) > > functions = { extract_atoms, extract_bonds,
                                                            extract_angles, extract_dihedrals };

        if (uses_parallel)
        {
            tbb::parallel_for( tbb::blocked_range<int>(0,functions.count(),1),
                            [&](const tbb::blocked_range<int> &r)
            {
                for (int i=r.begin(); i<r.end(); ++i)
                {
                    functions[i](false);
                    functions[i](true);
                }
            });
        }
        else
        {
            for (int i=0; i<functions.count(); ++i)
            {
                functions[i](false);
                functions[i](true);
            }
        }

        //sanitise this object
        this->_pvt_sanitise();
        this->_pvt_sanitise(true);
    }

    // Regular molecule.
    else
    {
        //get the name either from the molecule name or the name of the first
        //residue
        nme = mol.name();

        if (nme.isEmpty())
        {
            nme = mol.residue(ResIdx(0)).name();
        }

        //get the forcefield for this molecule
        try
        {
            ffield0 = mol.property(map["forcefield"]).asA<MMDetail>();
        }
        catch(...)
        {
            warns.append( QObject::tr("Cannot find a valid MM forcefield for this molecule!") );
        }

        const auto molinfo = mol.info();

        bool uses_parallel = true;
        if (map["parallel"].hasValue())
        {
            uses_parallel = map["parallel"].value().asA<BooleanProperty>().value();
        }

        //get information about all atoms in this molecule
        auto extract_atoms = [&]()
        {
            atms0 = QVector<GroAtom>(molinfo.nAtoms());

            AtomMasses masses;
            AtomElements elements;
            AtomCharges charges;
            AtomIntProperty groups;
            AtomStringProperty atomtypes;
            AtomStringProperty bondtypes;

            bool has_mass(false), has_elem(false), has_chg(false), has_group(false), has_type(false);
            bool has_bondtype(false);

            try
            {
                masses = mol.property(map["mass"]).asA<AtomMasses>();
                has_mass = true;
            }
            catch(...)
            {}

            if (not has_mass)
            {
                try
                {
                    elements = mol.property(map["element"]).asA<AtomElements>();
                    has_elem = true;
                }
                catch(...)
                {}
            }

            try
            {
                charges = mol.property(map["charge"]).asA<AtomCharges>();
                has_chg = true;
            }
            catch(...)
            {}

            try
            {
                groups = mol.property(map["charge_group"]).asA<AtomIntProperty>();
                has_group = true;
            }
            catch(...)
            {}

            try
            {
                atomtypes = mol.property(map["atomtype"]).asA<AtomStringProperty>();
                has_type = true;
            }
            catch(...)
            {}

            try
            {
                bondtypes = mol.property(map["bondtype"]).asA<AtomStringProperty>();
                has_bondtype = true;
            }
            catch(...)
            {}

            if (not (has_chg and has_type and (has_elem or has_mass)))
            {
                warns.append( QObject::tr("Cannot find valid charge, atomtype and (element or mass) "
                "properties for the molecule. These are needed! "
                "has_charge=%1, has_atomtype=%2, has_mass=%3, has_element=%4")
                    .arg(has_chg).arg(has_type).arg(has_mass).arg(has_elem) );
                return;
            }

            //run through the atoms in AtomIdx order
            auto extract_atom = [&](int iatm)
            {
                AtomIdx i(iatm);

                const auto cgatomidx = molinfo.cgAtomIdx(i);
                const auto residx = molinfo.parentResidue(i);

                //atom numbers have to count up sequentially from 1
                int atomnum = i+1;
                QString atomnam = molinfo.name(i);

                //assuming that residues are in the same order as the atoms
                int resnum = residx + 1;
                QString resnam = molinfo.name(residx);

                //people like to preserve the residue numbers of ligands and
                //proteins. This is very challenging for the gromacs topology,
                //as it would force a different topology for every solvent molecule,
                //so deciding on the difference between protein/ligand and solvent
                //is tough. Will preserve the residue number if the number of
                //residues is greater than 1 and the number of atoms is greater
                //than 32 (so octanol is a solvent)
                if (molinfo.nResidues() > 1 or molinfo.nAtoms() > 32)
                {
                    resnum = molinfo.number(residx).value();
                }

                int group = atomnum;

                if (has_group)
                {
                    group = groups[cgatomidx];
                }

                auto charge = charges[cgatomidx];

                SireUnits::Dimension::MolarMass mass;

                if (has_mass)
                {
                    mass = masses[cgatomidx];
                }
                else
                {
                    mass = elements[cgatomidx].mass();
                }

                if (mass <= 0)
                {
                    //not allowed to have a zero or negative mass
                    mass = 1.0 * g_per_mol;
                }

                QString atomtype = atomtypes[cgatomidx];

                auto &atom = atms0[i];
                atom.setName(atomnam);
                atom.setNumber(atomnum);
                atom.setResidueName(resnam);
                atom.setResidueNumber(resnum);
                atom.setChargeGroup(group);
                atom.setCharge(charge);
                atom.setMass(mass);
                atom.setAtomType(atomtype);

                if (has_bondtype)
                {
                    atom.setBondType(bondtypes[cgatomidx]);
                }
                else
                {
                    atom.setBondType(atomtype);
                }
            };

            if (uses_parallel)
            {
                tbb::parallel_for( tbb::blocked_range<int>(0,molinfo.nAtoms()),
                                [&](const tbb::blocked_range<int> &r)
                {
                    for (int i=r.begin(); i<r.end(); ++i)
                    {
                        extract_atom(i);
                    }
                });
            }
            else
            {
                for (int i=0; i<molinfo.nAtoms(); ++i)
                {
                    extract_atom(i);
                }
            }
        };

        //get all of the bonds in this molecule
        auto extract_bonds = [&]()
        {
            bool has_conn(false), has_funcs(false);

            TwoAtomFunctions funcs;
            Connectivity conn;

            const auto R = InternalPotential::symbols().bond().r();

            try
            {
                funcs = mol.property(map["bond"]).asA<TwoAtomFunctions>();
                has_funcs = true;
            }
            catch(...)
            {}

            try
            {
                conn = mol.property(map["connectivity"]).asA<Connectivity>();
                has_conn = true;
            }
            catch(...)
            {}

            //get the bond potentials first
            if (has_funcs)
            {
                for (const auto bond : funcs.potentials())
                {
                    AtomIdx atom0 = molinfo.atomIdx(bond.atom0());
                    AtomIdx atom1 = molinfo.atomIdx(bond.atom1());

                    if (atom0 > atom1)
                        qSwap(atom0,atom1);

                    bnds0.insert( BondID(atom0,atom1), GromacsBond(bond.function(), R) );
                }
            }

            //now fill in any missing bonded atoms with null bonds
            if (has_conn)
            {
                for (const auto bond : conn.getBonds())
                {
                    AtomIdx atom0 = molinfo.atomIdx(bond.atom0());
                    AtomIdx atom1 = molinfo.atomIdx(bond.atom1());

                    if (atom0 > atom1)
                        qSwap(atom0,atom1);

                    BondID b(atom0,atom1);

                    if (not bnds0.contains(b))
                    {
                        bnds0.insert(b, GromacsBond(5));  // function 5 is a simple connection
                    }
                }
            }
        };

        //get all of the angles in this molecule
        auto extract_angles = [&]()
        {
            bool has_funcs(false);

            ThreeAtomFunctions funcs;

            const auto theta = InternalPotential::symbols().angle().theta();

            try
            {
                funcs = mol.property(map["angle"]).asA<ThreeAtomFunctions>();
                has_funcs = true;
            }
            catch(...)
            {}

            if (has_funcs)
            {
                for (const auto angle : funcs.potentials())
                {
                    AtomIdx atom0 = molinfo.atomIdx(angle.atom0());
                    AtomIdx atom1 = molinfo.atomIdx(angle.atom1());
                    AtomIdx atom2 = molinfo.atomIdx(angle.atom2());

                    if (atom0 > atom2)
                        qSwap(atom0,atom2);

                    angs0.insert( AngleID(atom0,atom1,atom2),
                                GromacsAngle(angle.function(), theta) );
                }
            }
        };

        //get all of the dihedrals in this molecule
        auto extract_dihedrals = [&]()
        {
            bool has_funcs(false);

            FourAtomFunctions funcs;

            const auto phi = InternalPotential::symbols().dihedral().phi();

            try
            {
                funcs = mol.property(map["dihedral"]).asA<FourAtomFunctions>();
                has_funcs = true;
            }
            catch(...)
            {}

            if (has_funcs)
            {
                for (const auto dihedral : funcs.potentials())
                {
                    AtomIdx atom0 = molinfo.atomIdx(dihedral.atom0());
                    AtomIdx atom1 = molinfo.atomIdx(dihedral.atom1());
                    AtomIdx atom2 = molinfo.atomIdx(dihedral.atom2());
                    AtomIdx atom3 = molinfo.atomIdx(dihedral.atom3());

                    if (atom0 > atom3)
                    {
                        qSwap(atom0,atom3);
                        qSwap(atom1,atom2);
                    }

                    //get all of the dihedral terms (could be a lot)
                    auto parts = GromacsDihedral::construct(dihedral.function(), phi);

                    DihedralID dihid(atom0,atom1,atom2,atom3);

                    for (const auto part : parts)
                    {
                        dihs0.insertMulti(dihid, part);
                    }
                }
            }

            bool has_imps(false);

            FourAtomFunctions imps;

            try
            {
                imps = mol.property(map["improper"]).asA<FourAtomFunctions>();
                has_imps = true;
            }
            catch(...)
            {}

            if (has_imps)
            {
                for (const auto improper : imps.potentials())
                {
                    AtomIdx atom0 = molinfo.atomIdx(improper.atom0());
                    AtomIdx atom1 = molinfo.atomIdx(improper.atom1());
                    AtomIdx atom2 = molinfo.atomIdx(improper.atom2());
                    AtomIdx atom3 = molinfo.atomIdx(improper.atom3());

                    //get all of the dihedral terms (could be a lot)
                    auto parts = GromacsDihedral::constructImproper(improper.function(), phi);

                    DihedralID impid(atom0,atom1,atom2,atom3);

                    for (const auto part : parts)
                    {
                        dihs0.insertMulti(impid, part);
                    }
                }
            }
        };

        const QVector< std::function< void() > > functions = { extract_atoms, extract_bonds,
                                                            extract_angles, extract_dihedrals };

        if (uses_parallel)
        {
            tbb::parallel_for( tbb::blocked_range<int>(0,functions.count(),1),
                            [&](const tbb::blocked_range<int> &r)
            {
                for (int i=r.begin(); i<r.end(); ++i)
                {
                    functions[i]();
                }
            });
        }
        else
        {
            for (int i=0; i<functions.count(); ++i)
            {
                functions[i]();
            }
        }

        //sanitise this object
        this->_pvt_sanitise();
    }
}

/** Copy constructor */
GroMolType::GroMolType(const GroMolType &other)
           : nme(other.nme), warns(other.warns),
             atms0(other.atms0), atms1(other.atms1), first_atoms0(other.first_atoms0),
             first_atoms1(other.first_atoms1), bnds0(other.bnds0), bnds1(other.bnds1),
             angs0(other.angs0), angs1(other.angs1), dihs0(other.dihs0), dihs1(other.dihs1),
             ffield0(other.ffield0), ffield1(other.ffield1), nexcl0(other.nexcl0),
             nexcl1(other.nexcl1), is_perturbable(other.is_perturbable)
{}

/** Destructor */
GroMolType::~GroMolType()
{}

/** Copy assignment operator */
GroMolType& GroMolType::operator=(const GroMolType &other)
{
    if (this != &other)
    {
        nme = other.nme;
        warns = other.warns;
        atms0 = other.atms0;
        atms1 = other.atms1;
        first_atoms0 = other.first_atoms0;
        first_atoms1 = other.first_atoms1;
        bnds0 = other.bnds0;
        bnds1 = other.bnds1;
        angs0 = other.angs0;
        angs1 = other.angs1;
        dihs0 = other.dihs0;
        dihs1 = other.dihs1;
        ffield0 = other.ffield0;
        ffield1 = other.ffield1;
        nexcl0 = other.nexcl0;
        nexcl0 = other.nexcl1;
        is_perturbable = other.is_perturbable;
    }

    return *this;
}

/** Comparison operator */
bool GroMolType::operator==(const GroMolType &other) const
{
    return nme == other.nme and
           warns == other.warns and
           atms0 == other.atms0 and
           atms1 == other.atms1 and
           first_atoms0 == other.first_atoms0 and
           first_atoms1 == other.first_atoms1 and
           bnds0 == other.bnds0 and
           bnds1 == other.bnds1 and
           angs0 == other.angs0 and
           angs1 == other.angs1 and
           dihs0 == other.dihs0 and
           dihs1 == other.dihs1 and
           nexcl0 == other.nexcl0 and
           nexcl1 == other.nexcl1 and
           is_perturbable == other.is_perturbable;
}

/** Comparison operator */
bool GroMolType::operator!=(const GroMolType &other) const
{
    return not operator==(other);
}

const char* GroMolType::typeName()
{
    return QMetaType::typeName( qMetaTypeId<GroMolType>() );
}

const char* GroMolType::what() const
{
    return GroMolType::typeName();
}

/** Return whether or not this object is null */
bool GroMolType::isNull() const
{
    return this->operator==( GroMolType() );
}

/** Return whether or not this molecule is perturbable. */
bool GroMolType::isPerturbable() const
{
    return this->is_perturbable;
}

/** Return a string form for this object */
QString GroMolType::toString() const
{
    if (this->isNull())
        return QObject::tr("GroMolType::null");

    return QObject::tr("GroMolType( name() = '%1', nExcludedAtoms() = %2 )")
                    .arg(name()).arg(nExcludedAtoms());
}

/** Set the name of this moleculetype */
void GroMolType::setName(const QString &name)
{
    nme = name;
}

/** Return the name of this moleculetype */
QString GroMolType::name() const
{
    return nme;
}

/** Return the guessed forcefield for this molecule type */
MMDetail GroMolType::forcefield(bool is_lambda1) const
{
    // The molecule is not perturbable!
    if (is_lambda1 and not this->is_perturbable)
        throw SireError::incompatible_error(QObject::tr("The molecule isn't perturbable!"));

    if (is_lambda1)
        return ffield1;
    else
        return ffield0;
}

/** Set the number of excluded atoms */
void GroMolType::setNExcludedAtoms(qint64 n, bool is_lambda1)
{
    // The molecule is not perturbable!
    if (is_lambda1 and not this->is_perturbable)
        throw SireError::incompatible_error(QObject::tr("The molecule isn't perturbable!"));

    if (n >= 0)
        if (is_lambda1)
            nexcl1 = n;
        else
            nexcl0 = n;
    else
    {
        if (is_lambda1)
            nexcl1 = 0;
        else
            nexcl0 = 0;
    }
}

/** Return the number of excluded atoms */
qint64 GroMolType::nExcludedAtoms(bool is_lambda1) const
{
    // The molecule is not perturbable!
    if (is_lambda1 and not this->is_perturbable)
        throw SireError::incompatible_error(QObject::tr("The molecule isn't perturbable!"));

    if (is_lambda1)
        return nexcl1;
    else
        return nexcl0;
}

/** Add an atom to this moleculetype, with specified atom type, residue number,
    residue name, atom name, charge group, charge and mass */
void GroMolType::addAtom(const GroAtom &atom, bool is_lambda1)
{
    // The molecule is not perturbable!
    if (is_lambda1 and not this->is_perturbable)
        throw SireError::incompatible_error(QObject::tr("The molecule isn't perturbable!"));

    if (not atom.isNull())
    {
        if (is_lambda1)
            atms1.append( atom );
        else
            atms0.append( atom );
    }
}

/** Return whether or not this molecule needs sanitising */
bool GroMolType::needsSanitising(bool is_lambda1) const
{
    // The molecule is not perturbable!
    if (is_lambda1 and not this->is_perturbable)
        throw SireError::incompatible_error(QObject::tr("The molecule isn't perturbable!"));

    if (is_lambda1)
    {
        if (atms1.isEmpty())
            return false;
        else
            return ffield1.isNull() or first_atoms1.isEmpty();
    }
    else
    {
        if (atms0.isEmpty())
            return false;
        else
            return ffield0.isNull() or first_atoms0.isEmpty();
    }
}

/** Return the number of atoms in this molecule */
int GroMolType::nAtoms(bool is_lambda1) const
{
    // The molecule is not perturbable!
    if (is_lambda1 and not this->is_perturbable)
        throw SireError::incompatible_error(QObject::tr("The molecule isn't perturbable!"));

    if (needsSanitising(is_lambda1))
    {
        GroMolType other(*this);
        other._pvt_sanitise(is_lambda1);
        return other.nAtoms(is_lambda1);
    }
    else
    {
        if (is_lambda1)
            return atms1.count();
        else
            return atms0.count();
    }
}

/** Return the number of residues in this molecule */
int GroMolType::nResidues(bool is_lambda1) const
{
    // The molecule is not perturbable!
    if (is_lambda1 and not this->is_perturbable)
        throw SireError::incompatible_error(QObject::tr("The molecule isn't perturbable!"));

    if (needsSanitising(is_lambda1))
    {
        GroMolType other(*this);
        other._pvt_sanitise(is_lambda1);
        return other.nResidues(is_lambda1);
    }
    else
    {
        if (is_lambda1)
            return first_atoms1.count();
        else
            return first_atoms0.count();
    }
}

/** Return the atom at index 'atomidx' */
GroAtom GroMolType::atom(const AtomIdx &atomidx, bool is_lambda1) const
{
    // The molecule is not perturbable!
    if (is_lambda1 and not this->is_perturbable)
        throw SireError::incompatible_error(QObject::tr("The molecule isn't perturbable!"));

    if (needsSanitising(is_lambda1))
    {
        GroMolType other(*this);
        other._pvt_sanitise(is_lambda1);
        return other.atom(atomidx, is_lambda1);
    }
    else
    {
        if (is_lambda1)
        {
            int i = atomidx.map(atms1.count());
            return atms1.constData()[i];
        }
        else
        {
            int i = atomidx.map(atms0.count());
            return atms0.constData()[i];
        }
    }
}

/** Return the atom with number 'atomnum' */
GroAtom GroMolType::atom(const AtomNum &atomnum, bool is_lambda1) const
{
    // The molecule is not perturbable!
    if (is_lambda1 and not this->is_perturbable)
        throw SireError::incompatible_error(QObject::tr("The molecule isn't perturbable!"));

    if (needsSanitising(is_lambda1))
    {
        GroMolType other(*this);
        other._pvt_sanitise(is_lambda1);
        return other.atom(atomnum, is_lambda1);
    }

    if (is_lambda1)
    {
        for (int i=0; i<atms1.count(); ++i)
        {
            if (atms1.constData()[i].number() == atomnum)
            {
                return atms1.constData()[i];
            }
        }
    }
    else
    {
        for (int i=0; i<atms0.count(); ++i)
        {
            if (atms0.constData()[i].number() == atomnum)
            {
                return atms0.constData()[i];
            }
        }
    }

    throw SireMol::missing_atom( QObject::tr(
        "There is no atom with number '%1' in this molecule '%2'")
            .arg(atomnum.toString()).arg(this->toString()), CODELOC );
}

/** Return the first atom with name 'atomnam'. If you want all atoms
    with this name then call 'atoms(AtomName atomname)' */
GroAtom GroMolType::atom(const AtomName &atomnam, bool is_lambda1) const
{
    // The molecule is not perturbable!
    if (is_lambda1 and not this->is_perturbable)
        throw SireError::incompatible_error(QObject::tr("The molecule isn't perturbable!"));

    if (needsSanitising(is_lambda1))
    {
        GroMolType other(*this);
        other._pvt_sanitise(is_lambda1);
        return other.atom(atomnam, is_lambda1);
    }

    if (is_lambda1)
    {
        for (int i=0; i<atms1.count(); ++i)
        {
            if (atms1.constData()[i].name() == atomnam)
            {
                return atms1.constData()[i];
            }
        }
    }
    else
    {
        for (int i=0; i<atms0.count(); ++i)
        {
            if (atms0.constData()[i].name() == atomnam)
            {
                return atms0.constData()[i];
            }
        }
    }

    throw SireMol::missing_atom( QObject::tr(
        "There is no atom with name '%1' in this molecule '%2'")
            .arg(atomnam.toString()).arg(this->toString()), CODELOC );
}

/** Return all atoms that have the passed name. Returns an empty
    list if there are no atoms with this name */
QVector<GroAtom> GroMolType::atoms(const AtomName &atomnam, bool is_lambda1) const
{
    // The molecule is not perturbable!
    if (is_lambda1 and not this->is_perturbable)
        throw SireError::incompatible_error(QObject::tr("The molecule isn't perturbable!"));

    if (needsSanitising(is_lambda1))
    {
        GroMolType other(*this);
        other._pvt_sanitise(is_lambda1);
        return other.atoms(atomnam, is_lambda1);
    }

    QVector<GroAtom> ret;

    if (is_lambda1)
    {
        for (int i=0; i<atms1.count(); ++i)
        {
            if (atms1.constData()[i].name() == atomnam)
            {
                ret.append( atms1.constData()[i] );
            }
        }
    }
    else
    {
        for (int i=0; i<atms0.count(); ++i)
        {
            if (atms0.constData()[i].name() == atomnam)
            {
                ret.append( atms0.constData()[i] );
            }
        }
    }

    return ret;
}

/** Return all of the atoms in this molecule */
QVector<GroAtom> GroMolType::atoms(bool is_lambda1) const
{
    // The molecule is not perturbable!
    if (is_lambda1 and not this->is_perturbable)
        throw SireError::incompatible_error(QObject::tr("The molecule isn't perturbable!"));

    if (needsSanitising(is_lambda1))
    {
        GroMolType other(*this);
        other._pvt_sanitise(is_lambda1);
        return other.atoms(is_lambda1);
    }

    if (is_lambda1)
        return atms1;
    else
        return atms0;
}

/** Return all of the atoms in the specified residue */
QVector<GroAtom> GroMolType::atoms(const ResIdx &residx, bool is_lambda1) const
{
    // The molecule is not perturbable!
    if (is_lambda1 and not this->is_perturbable)
        throw SireError::incompatible_error(QObject::tr("The molecule isn't perturbable!"));

    if (needsSanitising(is_lambda1))
    {
        GroMolType other(*this);
        other._pvt_sanitise(is_lambda1);
        return other.atoms(residx, is_lambda1);
    }

    if (is_lambda1)
    {
        int ires = residx.map( first_atoms1.count() );

        int start = first_atoms1.constData()[ires];
        int end = atms1.count();

        if (ires+1 < first_atoms1.count())
        {
            end = first_atoms1.constData()[ires+1];
        }

        return atms1.mid(start, end);
    }
    else
    {
        int ires = residx.map( first_atoms0.count() );

        int start = first_atoms0.constData()[ires];
        int end = atms0.count();

        if (ires+1 < first_atoms0.count())
        {
            end = first_atoms0.constData()[ires+1];
        }

        return atms0.mid(start, end);
    }
}

/** Return all of the atoms in the specified residue(s) */
QVector<GroAtom> GroMolType::atoms(const ResNum &resnum, bool is_lambda1) const
{
    // The molecule is not perturbable!
    if (is_lambda1 and not this->is_perturbable)
        throw SireError::incompatible_error(QObject::tr("The molecule isn't perturbable!"));

    if (needsSanitising(is_lambda1))
    {
        GroMolType other(*this);
        other._pvt_sanitise(is_lambda1);
        return other.atoms(resnum, is_lambda1);
    }

    //find the indicies all all matching residues
    QList<ResIdx> idxs;

    if (is_lambda1)
    {
        for (int idx=0; idx<first_atoms1.count(); ++idx)
        {
            if (atms1[first_atoms1.at(idx)].residueNumber() == resnum)
            {
                idxs.append( ResIdx(idx) );
            }
        }
    }
    else
    {
        for (int idx=0; idx<first_atoms0.count(); ++idx)
        {
            if (atms0[first_atoms0.at(idx)].residueNumber() == resnum)
            {
                idxs.append( ResIdx(idx) );
            }
        }
    }

    QVector<GroAtom> ret;

    for (const auto idx : idxs)
    {
        ret += this->atoms(idx, is_lambda1);
    }

    return ret;
}

/** Return all of the atoms in the specified residue(s) */
QVector<GroAtom> GroMolType::atoms(const ResName &resnam, bool is_lambda1) const
{
    // The molecule is not perturbable!
    if (is_lambda1 and not this->is_perturbable)
        throw SireError::incompatible_error(QObject::tr("The molecule isn't perturbable!"));

    if (needsSanitising(is_lambda1))
    {
        GroMolType other(*this);
        other._pvt_sanitise(is_lambda1);
        return other.atoms(resnam, is_lambda1);
    }

    //find the indicies all all matching residues
    QList<ResIdx> idxs;

    if (is_lambda1)
    {
        for (int idx=0; idx<first_atoms1.count(); ++idx)
        {
            if (atms1[first_atoms1.at(idx)].residueName() == resnam)
            {
                idxs.append( ResIdx(idx) );
            }
        }
    }
    else
    {
        for (int idx=0; idx<first_atoms0.count(); ++idx)
        {
            if (atms0[first_atoms0.at(idx)].residueName() == resnam)
            {
                idxs.append( ResIdx(idx) );
            }
        }
    }

    QVector<GroAtom> ret;

    for (const auto idx : idxs)
    {
        ret += this->atoms(idx, is_lambda1);
    }

    return ret;
}

/** Internal function to do the non-forcefield parts of sanitising */
void GroMolType::_pvt_sanitise(bool is_lambda1)
{
    // The molecule is not perturbable!
    if (is_lambda1 and not this->is_perturbable)
        throw SireError::incompatible_error(QObject::tr("The molecule isn't perturbable!"));

    //sort the atoms so that they are in residue number / atom number order, and
    //we check and remove duplicate atom numbers

    if (is_lambda1)
        first_atoms1.append(0);
    else
        first_atoms0.append(0);
}

/** Sanitise this moleculetype. This assumes that the moleculetype has
    been fully specified, so it collects everything together and checks that the
    molecule makes sense. Any warnings generated can be retrieved using the
    'warnings' function. It also uses the passed defaults from the top file,
    together with the information in the molecule to guess the forcefield for
    the molecule */
void GroMolType::sanitise(QString elecstyle, QString vdwstyle, QString combrule,
                          double elec14, double vdw14, bool is_lambda1)
{
    // The molecule is not perturbable!
    if (is_lambda1 and not this->is_perturbable)
        throw SireError::incompatible_error(QObject::tr("The molecule isn't perturbable!"));

    if (not needsSanitising(is_lambda1))
        return;

    this->_pvt_sanitise(is_lambda1);

    //also check that the bonds/angles/dihedrals all refer to actual atoms...

    //work out the bond, angle and dihedral function styles. We will
    //do this assuming that anything other than simple harmonic/cosine is
    //an "interesting" gromacs-style forcefield
    QString bondstyle = "harmonic";

    if (is_lambda1)
    {
        for (auto it = bnds1.constBegin(); it != bnds1.constEnd(); ++it)
        {
            if (not (it.value().isSimple() and it.value().isHarmonic()))
            {
                bondstyle = "gromacs";
                break;
            }
        }

        QString anglestyle = "harmonic";

        for (auto it = angs1.constBegin(); it != angs1.constEnd(); ++it)
        {
            if (not (it.value().isSimple() and it.value().isHarmonic()))
            {
                anglestyle = "gromacs";
                break;
            }
        }

        QString dihedralstyle = "cosine";

        for (auto it = dihs1.constBegin(); it != dihs1.constEnd(); ++it)
        {
            if (not (it.value().isSimple() and it.value().isCosine()))
            {
                dihedralstyle = "gromacs";
                break;
            }
        }

        //finally generate a forcefield description for this molecule based on the
        //passed defaults and the functional forms of the internals
        ffield1 = MMDetail::guessFrom(combrule, elecstyle, vdwstyle, elec14, vdw14,
                                      bondstyle, anglestyle, dihedralstyle);
    }
    else
    {
        for (auto it = bnds0.constBegin(); it != bnds0.constEnd(); ++it)
        {
            if (not (it.value().isSimple() and it.value().isHarmonic()))
            {
                bondstyle = "gromacs";
                break;
            }
        }

        QString anglestyle = "harmonic";

        for (auto it = angs0.constBegin(); it != angs0.constEnd(); ++it)
        {
            if (not (it.value().isSimple() and it.value().isHarmonic()))
            {
                anglestyle = "gromacs";
                break;
            }
        }

        QString dihedralstyle = "cosine";

        for (auto it = dihs0.constBegin(); it != dihs0.constEnd(); ++it)
        {
            if (not (it.value().isSimple() and it.value().isCosine()))
            {
                dihedralstyle = "gromacs";
                break;
            }
        }

        //finally generate a forcefield description for this molecule based on the
        //passed defaults and the functional forms of the internals
        ffield0 = MMDetail::guessFrom(combrule, elecstyle, vdwstyle, elec14, vdw14,
                                      bondstyle, anglestyle, dihedralstyle);
    }
}

/** Add a warning that has been generated while parsing or creatig this object */
void GroMolType::addWarning(const QString &warning)
{
    warns.append(warning);
}

/** Return any warnings associated with this moleculetype */
QStringList GroMolType::warnings() const
{
    return warns;
}

/** Add the passed bond to the molecule */
void GroMolType::addBond(const BondID &bond, const GromacsBond &param, bool is_lambda1)
{
    // The molecule is not perturbable!
    if (is_lambda1 and not this->is_perturbable)
        throw SireError::incompatible_error(QObject::tr("The molecule isn't perturbable!"));

    if (is_lambda1)
        bnds1.insertMulti(bond, param);
    else
        bnds0.insertMulti(bond, param);
}

/** Add the passed angle to the molecule */
void GroMolType::addAngle(const AngleID &angle, const GromacsAngle &param,
                          bool is_lambda1)
{
    // The molecule is not perturbable!
    if (is_lambda1 and not this->is_perturbable)
        throw SireError::incompatible_error(QObject::tr("The molecule isn't perturbable!"));

    if (is_lambda1)
        angs1.insertMulti(angle, param);
    else
        angs0.insertMulti(angle, param);
}

/** Add the passed dihedral to the molecule */
void GroMolType::addDihedral(const DihedralID &dihedral, const GromacsDihedral &param,
                             bool is_lambda1)
{
    // The molecule is not perturbable!
    if (is_lambda1 and not this->is_perturbable)
        throw SireError::incompatible_error(QObject::tr("The molecule isn't perturbable!"));

    if (is_lambda1)
        dihs1.insertMulti(dihedral, param);
    else
        dihs0.insertMulti(dihedral, param);
}

/** Add the passed bonds to the molecule */
void GroMolType::addBonds(const QMultiHash<BondID,GromacsBond> &bonds,
                          bool is_lambda1)
{
    // The molecule is not perturbable!
    if (is_lambda1 and not this->is_perturbable)
        throw SireError::incompatible_error(QObject::tr("The molecule isn't perturbable!"));

    if (is_lambda1)
        bnds1 += bonds;
    else
        bnds0 += bonds;
}

/** Add the passed angles to the molecule */
void GroMolType::addAngles(const QMultiHash<AngleID,GromacsAngle> &angles,
                           bool is_lambda1)
{
    // The molecule is not perturbable!
    if (is_lambda1 and not this->is_perturbable)
        throw SireError::incompatible_error(QObject::tr("The molecule isn't perturbable!"));

    if (is_lambda1)
        angs1 += angles;
    else
        angs0 += angles;
}

/** Add the passed dihedrals to the molecule */
void GroMolType::addDihedrals(const QMultiHash<DihedralID,GromacsDihedral> &dihedrals,
                              bool is_lambda1)
{
    // The molecule is not perturbable!
    if (is_lambda1 and not this->is_perturbable)
        throw SireError::incompatible_error(QObject::tr("The molecule isn't perturbable!"));

    if (is_lambda1)
        dihs1 += dihedrals;
    else
        dihs0 += dihedrals;
}

/** Return all of the bonds */
QMultiHash<BondID,GromacsBond> GroMolType::bonds(bool is_lambda1) const
{
    // The molecule is not perturbable!
    if (is_lambda1 and not this->is_perturbable)
        throw SireError::incompatible_error(QObject::tr("The molecule isn't perturbable!"));

    if (is_lambda1)
        return bnds1;
    else
        return bnds0;
}

/** Return all of the angles */
QMultiHash<AngleID,GromacsAngle> GroMolType::angles(bool is_lambda1) const
{
    // The molecule is not perturbable!
    if (is_lambda1 and not this->is_perturbable)
        throw SireError::incompatible_error(QObject::tr("The molecule isn't perturbable!"));

    if (is_lambda1)
        return angs1;
    else
        return angs0;
}

/** Return all of the dihedrals */
QMultiHash<DihedralID,GromacsDihedral> GroMolType::dihedrals(bool is_lambda1) const
{
    // The molecule is not perturbable!
    if (is_lambda1 and not this->is_perturbable)
        throw SireError::incompatible_error(QObject::tr("The molecule isn't perturbable!"));

    if (is_lambda1)
        return dihs1;
    else
        return dihs0;
}

/** Return whether or not this is a topology for water. This should
    return true for all water models (including TIP4P and TIP5P) */
bool GroMolType::isWater(bool is_lambda1) const
{
    // The molecule is not perturbable!
    if (is_lambda1 and not this->is_perturbable)
        throw SireError::incompatible_error(QObject::tr("The molecule isn't perturbable!"));

    if (nResidues(is_lambda1) == 1)
    {
        if (nAtoms(is_lambda1) >= 3 and nAtoms(is_lambda1) <= 5) // catch SPC/TIP3P to TIP5P
        {
            //the total mass of the molecule should be 18 (rounded)
            //and the number of oxygens should be 1 and hydrogens should be 2
            int noxy = 0;
            int nhyd = 0;
            int total_mass = 0;

            if (is_lambda1)
            {
                for (const auto &atm : atms1)
                {
                    //round down the mass to the nearest integer unit, so to
                    //exclude isotopes
                    const int mass = int( std::floor(atm.mass().value()) );

                    total_mass += mass;

                    if (total_mass > 18)
                        return false;

                    if (mass == 16)
                    {
                        noxy += 1;
                        if (noxy > 1)
                            return false;
                    }
                    else if (mass == 1)
                    {
                        nhyd += 1;
                        if (nhyd > 2)
                            return false;
                    }
                    else
                        //not an oxygen or hydrogen
                        return false;
                }

                //this is a water :-)
                return true;
            }
            else
            {
                for (const auto &atm : atms0)
                {
                    //round down the mass to the nearest integer unit, so to
                    //exclude isotopes
                    const int mass = int( std::floor(atm.mass().value()) );

                    total_mass += mass;

                    if (total_mass > 18)
                        return false;

                    if (mass == 16)
                    {
                        noxy += 1;
                        if (noxy > 1)
                            return false;
                    }
                    else if (mass == 1)
                    {
                        nhyd += 1;
                        if (nhyd > 2)
                            return false;
                    }
                    else
                        //not an oxygen or hydrogen
                        return false;
                }

                //this is a water :-)
                return true;
            }
        }
    }

    return false;
}

/** Return the settles lines for this molecule. This currently only returns
    settles lines for water molecules. These lines are used to constrain the
    bonds/angles of the water molecule */
QStringList GroMolType::settlesLines(bool is_lambda1) const
{
    // The molecule is not perturbable!
    if (is_lambda1 and not this->is_perturbable)
        throw SireError::incompatible_error(QObject::tr("The molecule isn't perturbable!"));

    if (not this->isWater(is_lambda1))
        return QStringList();

    QStringList lines;

    lines.append("[ settles ]");
    lines.append("; OW    funct   doh dhh");

    //find the OH and HH bonds to get the equilibrium OH and HH bond length
    //for this water molecule - if we don't have it, then use these as default (TIP3P)
    double hh_length = 0.15136000;
    double oh_length = 0.09572000;

    //there should only be two bonds - OH and HH. The longer one is HH
    if (is_lambda1)
    {
        if (bnds1.count() == 2)
        {
            auto it = bnds1.begin();

            double hh_length = it.value().equilibriumLength().to(nanometer);
            ++it;
            double oh_length = it.value().equilibriumLength().to(nanometer);

            if (oh_length > hh_length)
                qSwap(oh_length, hh_length);
        }
    }
    else
    {
        if (bnds0.count() == 2)
        {
            auto it = bnds0.begin();

            double hh_length = it.value().equilibriumLength().to(nanometer);
            ++it;
            double oh_length = it.value().equilibriumLength().to(nanometer);

            if (oh_length > hh_length)
                qSwap(oh_length, hh_length);
        }
    }

    lines.append( QString("1       1       %1 %2").arg(oh_length, 7, 'f', 5)
                                                  .arg(hh_length, 7, 'f', 5) );

    lines.append("");
    lines.append("[ exclusions ]");
    lines.append("1   2   3");
    lines.append("2   1   3");
    lines.append("3   1   2");

    return lines;
}

////////////////
//////////////// Implementation of GroSystem
////////////////

static const RegisterMetaType<GroSystem> r_grosys(NO_ROOT);

QDataStream SIREIO_EXPORT &operator<<(QDataStream &ds, const GroSystem &grosys)
{
    writeHeader(ds, r_grosys, 1);

    SharedDataStream sds(ds);

    sds << grosys.nme << grosys.moltypes << grosys.nmols;

    return ds;
}

QDataStream SIREIO_EXPORT &operator>>(QDataStream &ds, GroSystem &grosys)
{
    VersionID v = readHeader(ds, r_grosys);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> grosys.nme >> grosys.moltypes >> grosys.nmols;

        grosys.total_nmols = 0;

        for (auto it = grosys.nmols.constBegin();
             it != grosys.nmols.constEnd(); ++it)
        {
            grosys.total_nmols += *it;
        }
    }
    else
        throw version_error(v, "1", r_grosys, CODELOC);

    return ds;
}

/** Construct a null GroSystem */
GroSystem::GroSystem() : total_nmols(0)
{}

/** Construct a GroSystem with the passed name */
GroSystem::GroSystem(const QString &name) : nme(name), total_nmols(0)
{}

/** Copy constructor */
GroSystem::GroSystem(const GroSystem &other)
          : nme(other.nme), moltypes(other.moltypes), nmols(other.nmols),
            total_nmols(other.total_nmols)
{}

/** Destructor */
GroSystem::~GroSystem()
{}

/** Copy assignment operator */
GroSystem& GroSystem::operator=(const GroSystem &other)
{
    nme = other.nme;
    moltypes = other.moltypes;
    nmols = other.nmols;
    total_nmols = other.total_nmols;
    return *this;
}

/** Comparison operator */
bool GroSystem::operator==(const GroSystem &other) const
{
    return nme == other.nme and
           total_nmols == other.total_nmols and
           moltypes == other.moltypes and
           nmols == other.nmols;
}

/** Comparison operator */
bool GroSystem::operator!=(const GroSystem &other) const
{
    return not operator==(other);
}

/** Return the molecule type of the ith molecule */
QString GroSystem::operator[](int i) const
{
    i = Index(i).map(total_nmols);

    auto it2 = moltypes.constBegin();
    for (auto it = nmols.constBegin(); it != nmols.constEnd(); ++it)
    {
        if (i < *it)
        {
            return *it2;
        }
        else
        {
            i -= *it;
            ++it2;
        }
    }

    //we should never get here...
    throw SireError::program_bug( QObject::tr(
        "How did we get here? %1 : %2 : %3")
            .arg(i).arg( Sire::toString(moltypes) ).arg( Sire::toString(nmols) ), CODELOC );

    return QString();
}

/** Return the molecule type of the ith molecule */
QString GroSystem::at(int i) const
{
    return operator[](i);
}

/** Return the number of molecules in the system */
int GroSystem::size() const
{
    return total_nmols;
}

/** Return the number of molecules in the system */
int GroSystem::count() const
{
    return size();
}

/** Return the number of molecules in the system */
int GroSystem::nMolecules() const
{
    return size();
}

const char* GroSystem::typeName()
{
    return QMetaType::typeName( qMetaTypeId<GroSystem>() );
}

const char* GroSystem::what() const
{
    return GroSystem::typeName();
}

/** Return the name of the system */
QString GroSystem::name() const
{
    return nme;
}

/** Set the name of the system */
void GroSystem::setName( QString name )
{
    nme = name;
}

/** Return a string representation of this system */
QString GroSystem::toString() const
{
    if (this->isNull())
    {
        return QObject::tr( "GroSystem::null" );
    }
    else if (this->isEmpty())
    {
        return QObject::tr( "GroSystem( %1 : empty )" ).arg(this->name());
    }
    else
    {
        return QObject::tr( "GroSystem( %1 : nMolecules()=%2 )")
                    .arg(this->name()).arg(this->nMolecules());
    }
}

/** Return whether or not this is a null GroSystem */
bool GroSystem::isNull() const
{
    return nme.isNull() and total_nmols == 0;
}

/** Return whether or not this is an empty system (no molecules) */
bool GroSystem::isEmpty() const
{
    return total_nmols == 0;
}

/** Return the list of unique molecule types held in the system */
QStringList GroSystem::uniqueTypes() const
{
    QStringList typs;

    for (const auto moltype : moltypes)
    {
        if (not typs.contains(moltype))
        {
            typs.append(moltype);
        }
    }

    return typs;
}

/** Add (optionally ncopies) copies of the molecule with type 'moltype'
    to the system */
void GroSystem::add(QString moltype, int ncopies)
{
    if (ncopies <= 0)
        return;

    if (total_nmols > 0)
    {
        if (moltypes.back() == moltype)
        {
            nmols.back() += ncopies;
            total_nmols += ncopies;
            return;
        }
    }

    moltypes.append(moltype);
    nmols.append(ncopies);
    total_nmols += ncopies;
}

////////////////
//////////////// Implementation of GroTop
////////////////

const RegisterParser<GroTop> register_grotop;
static const RegisterMetaType<GroTop> r_grotop;

QDataStream SIREIO_EXPORT &operator<<(QDataStream &ds, const GroTop &grotop)
{
    writeHeader(ds, r_grotop, 1);

    SharedDataStream sds(ds);

    sds << grotop.include_path << grotop.included_files
        << grotop.expanded_lines
        << grotop.atom_types
        << grotop.bond_potentials
        << grotop.ang_potentials
        << grotop.dih_potentials
        << grotop.moltypes
        << grotop.grosys
        << grotop.nb_func_type
        << grotop.combining_rule << grotop.fudge_lj
        << grotop.fudge_qq << grotop.parse_warnings << grotop.generate_pairs
        << static_cast<const MoleculeParser&>(grotop);

    return ds;
}

QDataStream SIREIO_EXPORT &operator>>(QDataStream &ds, GroTop &grotop)
{
    VersionID v = readHeader(ds, r_grotop);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> grotop.include_path >> grotop.included_files
            >> grotop.expanded_lines
            >> grotop.atom_types
            >> grotop.bond_potentials
            >> grotop.ang_potentials
            >> grotop.dih_potentials
            >> grotop.moltypes
            >> grotop.grosys
            >> grotop.nb_func_type
            >> grotop.combining_rule >> grotop.fudge_lj
            >> grotop.fudge_qq >> grotop.parse_warnings >> grotop.generate_pairs
            >> static_cast<MoleculeParser&>(grotop);
    }
    else
        throw version_error(v, "1", r_grotop, CODELOC);

    return ds;
}

//first thing is to parse in the gromacs files. These use #include, #define, #if etc.
//so we need to pull all of them together into a single set of lines

/** Internal function to return a LJParameter from the passed W and V values
    for the passed Gromacs combining rule */
static LJParameter toLJParameter(double v, double w, int rule)
{
    if (rule == 2 or rule == 3)
    {
        // v = sigma in nm, and w = epsilon in kJ mol-1
        return LJParameter::fromSigmaAndEpsilon( v * nanometer, w * kJ_per_mol );
    }
    else
    {
        // v = 4 epsilon sigma^6 in kJ mol-1 nm^6, w = 4 epsilon sigma^12 in kJ mol-1 nm^12
        // so sigma = (w/v)^1/6 and epsilon = v^2 / 4w
        return LJParameter::fromSigmaAndEpsilon( std::pow(w/v, 1.0/6.0) * nanometer,
                                                 (v*v / (4.0*w)) * kJ_per_mol );
    }
}

/** Internal function to convert a LJParameter to V and W based on the passed gromacs
    combining rule */
static std::tuple<double,double> fromLJParameter(const LJParameter &lj, int rule)
{
    const double sigma = lj.sigma().to(nanometer);
    const double epsilon = lj.epsilon().to(kJ_per_mol);

    if (rule == 2 or rule == 3)
    {
        return std::make_tuple(sigma, epsilon);
    }
    else
    {
        double sig6 = SireMaths::pow(sigma, 6);
        double v = 4.0 * epsilon * sig6;
        double w = v * sig6;

        return std::make_tuple(v, w);
    }
}

/** Internal function to create a string version of the LJ function type */
static QString _getVDWStyle(int type)
{
    if (type == 1)
        return "lj";
    else if (type == 2)
        return "buckingham";
    else
        throw SireError::invalid_arg( QObject::tr(
            "Cannot find the VDW function type from value '%1'. Should be 1 or 2.")
                .arg(type), CODELOC );

    return QString();
}

/** Internal function to convert a MMDetail description of the LJ function type back
    to the gromacs integer */
static int _getVDWStyleFromFF(const MMDetail &ffield)
{
    if (ffield.usesLJTerm())
        return 1;
    else if (ffield.usesBuckinghamTerm())
        return 2;
    else
        throw SireError::invalid_arg( QObject::tr(
            "Cannot find the VDW function type for forcefield\n%1\n. "
            "This writer only support LJ or Buckingham VDW terms.").arg(ffield.toString()),
                CODELOC );

    return 0;
}

/** Internal function to create the string version of the combining rules */
static QString _getCombiningRules(int type)
{
    if (type == 1 or type == 3)
        return "geometric";
    else if (type == 2)
        return "arithmetic";
    else
        throw SireError::invalid_arg( QObject::tr(
            "Cannot find the combining rules type from value '%1'. Should be 1, 2 or 3.")
                .arg(type), CODELOC );

    return QString();
}

/** Internal function to create the combining rules from the passed forcefield */
int _getCombiningRulesFromFF(const MMDetail &ffield)
{
    if (ffield.usesGeometricCombiningRules())
        return 1;   //I don't know what 3 is...
    else if (ffield.usesArithmeticCombiningRules())
        return 2;
    else
        throw SireError::invalid_arg( QObject::tr(
            "Cannot find the combining rules to match the forcefield\n%1\n"
            "Valid options are arithmetic or geometric.").arg(ffield.toString()), CODELOC );

    return 0;
}

/** Constructor */
GroTop::GroTop()
       : ConcreteProperty<GroTop,MoleculeParser>(),
         nb_func_type(0), combining_rule(0), fudge_lj(0), fudge_qq(0),
         generate_pairs(false)
{}

/** This function gets the gromacs include path from the passed property map,
    as well as the current system environment */
void GroTop::getIncludePath(const PropertyMap &map)
{
    QStringList path;

    //now, see if the path is given in "GROMACS_PATH" in map
    try
    {
        const auto p = map["GROMACS_PATH"];

        if (p.hasValue())
        {
            path += p.value().asA<StringProperty>().toString().split(":", QString::SkipEmptyParts);
        }
        else if (p.source() != "GROMACS_PATH")
        {
            path += p.source().split(":", QString::SkipEmptyParts);
        }
    }
    catch(...)
    {}

    //now, see if the path is given in the "GROMACS_PATH" environment variable
    QString val = QString::fromLocal8Bit( qgetenv("GROMACS_PATH") );

    if (not val.isEmpty())
    {
        path += val.split(":", QString::SkipEmptyParts);
    }

    //now go through each path and convert it into an absolute path based on the
    //current directory
    for (const auto p : path)
    {
        include_path.append( QFileInfo(p).canonicalFilePath() );
    }
}

/** Construct to read in the data from the file called 'filename'. The
    passed property map can be used to pass extra parameters to control
    the parsing */
GroTop::GroTop(const QString &filename, const PropertyMap &map)
       : ConcreteProperty<GroTop,MoleculeParser>(filename,map),
         nb_func_type(0), combining_rule(0), fudge_lj(0), fudge_qq(0),
         generate_pairs(false)
{
    this->getIncludePath(map);

    //parse the data in the parse function, passing in the absolute path
    //to the directory that contains this file
    this->parseLines( QFileInfo(filename).absolutePath(), map);

    //now make sure that everything is correct with this object
    this->assertSane();
}

/** Construct to read in the data from the passed text lines. The
    passed property map can be used to pass extra parameters to control
    the parsing */
GroTop::GroTop(const QStringList &lines, const PropertyMap &map)
       : ConcreteProperty<GroTop,MoleculeParser>(lines,map),
         nb_func_type(0), combining_rule(0), fudge_lj(0), fudge_qq(0),
         generate_pairs(false)
{
    this->getIncludePath(map);

    //parse the data in the parse function, assuming the file has
    //come from the current directory
    this->parseLines(QDir::current().absolutePath(), map);

    //now make sure that everything is correct with this object
    this->assertSane();
}

/** Internal function used to generate the lines for the defaults section */
static QStringList writeDefaults(const MMDetail &ffield)
{
    QStringList lines;
    lines.append("; Gromacs Topology File written by Sire");
    lines.append(QString("; File written %1")
        .arg( QDateTime::currentDateTime().toString("MM/dd/yy  hh:mm:ss") ) );

    lines.append("[ defaults ]");
    lines.append("; nbfunc      comb-rule       gen-pairs      fudgeLJ     fudgeQQ");

    //all forcefields we support have gen-pairs = true (gromacs only understands 'yes' and 'no')
    const QString gen_pairs = "yes";

    lines.append( QString("  %1           %2               %3            %4         %5")
                        .arg( _getVDWStyleFromFF(ffield) )
                        .arg( _getCombiningRulesFromFF(ffield) )
                        .arg(gen_pairs)
                        .arg(ffield.vdw14ScaleFactor())
                        .arg(ffield.electrostatic14ScaleFactor()) );

    lines.append("");

    return lines;
}

/** Internal function used to write all of the atom types. This function requires
    that the atom types for all molecules are all consistent, i.e. atom type
    X has the same mass, vdw parameters, element type etc. for all molecules */
static QStringList writeAtomTypes(const QHash<QString,GroMolType> &moltyps,
                                  const QHash<QString,Molecule> &molecules,
                                  const MMDetail &ffield,
                                  const PropertyMap &map)
{
    //first, build up a dictionary of all of the unique atom types
    QHash<QString,QString> atomtypes;

    auto elemprop = map["element"];
    auto massprop = map["mass"];
    auto ljprop = map["LJ"];

    //get the combining rules - these determine the format of the LJ parameter in the file
    const int combining_rules = _getCombiningRulesFromFF(ffield);

    for (auto it = moltyps.constBegin(); it != moltyps.constEnd(); ++it)
    {
        const auto moltyp = it.value();

        // Store whether the molecule is perturbable.
        const auto is_perturbable = moltyp.isPerturbable();

        // Rename property keys.
        if (is_perturbable)
        {
            elemprop = "element0";
            massprop = "mass0";
            ljprop = "LJ0";
        }
        else
        {
            elemprop = map["element"];
            massprop = map["mass"];
            ljprop = map["LJ"];
        }

        const auto &atoms = moltyp.atoms();

        for (int i=0; i<atoms.count(); ++i)
        {
            const auto &atom = atoms[i];
            const auto atomtype = atom.atomType();

            if (atomtypes.contains(atomtype))
                continue;

            //we haven't seen this atom type before. Get the corresponding atom
            //in the molecule
            const auto mol = molecules[it.key()];
            const auto cgatomidx = mol.info().cgAtomIdx( AtomIdx(i) );

            //now get the corresponding Element and LJ properties for this atom
            Element elem;

            try
            {
                elem = mol.property(elemprop).asA<AtomElements>()[cgatomidx];
            }
            catch(...)
            {
                elem = Element::elementWithMass(
                            mol.property(massprop).asA<AtomMasses>()[cgatomidx] );
            }

            double chg = 0;  // always use a zero charge as this will be supplied with the atom

            auto lj = mol.property(ljprop).asA<AtomLJs>()[cgatomidx];
            auto ljparams = ::fromLJParameter(lj, combining_rules);

            QString particle_type = "A";  // A is for Atom

            if (elem.nProtons() == 0 and lj.isDummy())
            {
                particle_type = "D"; //this atomtype is a Dummy
            }

            atomtypes.insert( atomtype, QString("  %1        %2  %3  %4  %5  %6  %7")
                    .arg(atomtype, 4)
                    .arg(elem.nProtons(), 4)
                    .arg(elem.mass().to(g_per_mol), 10, 'f', 6)
                    .arg(chg, 10, 'f', 6)
                    .arg(particle_type, 3)
                    .arg(std::get<0>(ljparams), 10, 'f', 6)
                    .arg(std::get<1>(ljparams), 10, 'f', 6) );

        }

        // Add additional atom types from lambda = 1.
        if (is_perturbable)
        {
            const auto &atoms = moltyp.atoms(true);

            for (int i=0; i<atoms.count(); ++i)
            {
                const auto &atom = atoms[i];
                const auto atomtype = atom.atomType();

                if (atomtypes.contains(atomtype))
                    continue;

                //we haven't seen this atom type before. Get the corresponding atom
                //in the molecule
                const auto mol = molecules[it.key()];
                const auto cgatomidx = mol.info().cgAtomIdx( AtomIdx(i) );

                //now get the corresponding Element and LJ properties for this atom
                Element elem;

                try
                {
                    elem = mol.property("element1").asA<AtomElements>()[cgatomidx];
                }
                catch(...)
                {
                    elem = Element::elementWithMass(
                                mol.property("mass1").asA<AtomMasses>()[cgatomidx] );
                }

                double chg = 0;  // always use a zero charge as this will be supplied with the atom

                auto lj = mol.property("LJ1").asA<AtomLJs>()[cgatomidx];
                auto ljparams = ::fromLJParameter(lj, combining_rules);

                QString particle_type = "A";  // A is for Atom

                if (elem.nProtons() == 0 and lj.isDummy())
                {
                    particle_type = "D"; //this atomtype is a Dummy
                }

                atomtypes.insert( atomtype, QString("  %1        %2  %3  %4  %5  %6  %7")
                        .arg(atomtype, 4)
                        .arg(elem.nProtons(), 4)
                        .arg(elem.mass().to(g_per_mol), 10, 'f', 6)
                        .arg(chg, 10, 'f', 6)
                        .arg(particle_type, 3)
                        .arg(std::get<0>(ljparams), 10, 'f', 6)
                        .arg(std::get<1>(ljparams), 10, 'f', 6) );
            }
        }
    }

    //now sort and write all of the atomtypes
    QStringList lines;
    auto keys = atomtypes.keys();
    qSort(keys);

    lines.append( "[ atomtypes ]" );
    lines.append( "; name      at.num   mass         charge     ptype      sigma      epsilon" );

    for (const auto key : keys )
    {
        lines.append( atomtypes[key] );
    }

    lines.append("");

    return lines;
}

/** Internal function used to convert a Gromacs Moltyp to a set of lines */
static QStringList writeMolType(const QString &name, const GroMolType &moltype,
                                const Molecule &mol, bool uses_parallel)
{
    QStringList lines;

    lines.append( "[ moleculetype ]" );
    lines.append( "; name  nrexcl" );
    lines.append( QString("%1  %2").arg(name).arg(moltype.nExcludedAtoms()) );
    lines.append( "" );

    QStringList atomlines, bondlines, anglines, dihlines, scllines;

    // Store whether the molecule is perturbable.
    const auto is_perturbable = moltype.isPerturbable();

    //write all of the atoms
    auto write_atoms = [&]()
    {
        if (is_perturbable)
        {
            // Get the atoms from the molecule.
            const auto &atoms0 = moltype.atoms();
            const auto &atoms1 = moltype.atoms(true);

            // Loop over all of the atoms.
            for (int i=0; i<atoms0.count(); ++i)
            {
                const auto &atom0 = atoms0[i];
                const auto &atom1 = atoms1[i];

                atomlines.append( QString("%1   %2 %3    %4  %5   %6 %7   %8   %9 %10   %11")
                         .arg(atom0.number().value(), 6)
                         .arg(atom0.atomType(), 4)
                         .arg(atom0.residueNumber().value(), 6)
                         .arg(atom0.residueName().value(), 4)
                         .arg(atom0.name().value(), 4)
                         .arg(atom0.chargeGroup(), 4)
                         .arg(atom0.charge().to(mod_electron), 10, 'f', 6)
                         .arg(atom0.mass().to(g_per_mol), 10, 'f', 6)
                         .arg(atom1.atomType(), 4)
                         .arg(atom1.charge().to(mod_electron), 10, 'f', 6)
                         .arg(atom1.mass().to(g_per_mol), 10, 'f', 6) );
            }
        }
        else
        {
            // Get the atoms from the molecule.
            const auto &atoms = moltype.atoms();

            // Loop over all of the atoms.
            for (int i=0; i<atoms.count(); ++i)
            {
                const auto &atom = atoms[i];

                atomlines.append( QString("%1   %2 %3    %4  %5   %6 %7   %8")
                         .arg(atom.number().value(), 6)
                         .arg(atom.atomType(), 4)
                         .arg(atom.residueNumber().value(), 6)
                         .arg(atom.residueName().value(), 4)
                         .arg(atom.name().value(), 4)
                         .arg(atom.chargeGroup(), 4)
                         .arg(atom.charge().to(mod_electron), 10, 'f', 6)
                         .arg(atom.mass().to(g_per_mol), 10, 'f', 6) );
            }
        }

        atomlines.append( "" );
    };

    //write all of the bonds
    auto write_bonds = [&]()
    {
        if (is_perturbable)
        {
            // Get the bonds from the molecule.
            const auto &bonds0 = moltype.bonds();
            const auto &bonds1 = moltype.bonds(true);

            for (auto it = bonds0.constBegin(); it != bonds0.constEnd(); ++it)
            {
                const auto &bond = it.key();
                const auto &param = it.value();

                //AtomID is AtomIdx. Add 1, as gromacs is 1-indexed
                int atom0 = bond.atom0().asA<AtomIdx>().value() + 1;
                int atom1 = bond.atom1().asA<AtomIdx>().value() + 1;

                QStringList params0;
                for (const auto p : param.parameters())
                {
                    params0.append( QString::number(p) );
                }

                // Find the corresponding bond parameters at lambda = 1.
                const auto &param1 = bonds1.find(bond).value();

                QStringList params1;
                for (const auto p : param1.parameters())
                {
                    params1.append( QString::number(p) );
                }

                bondlines.append( QString("%1 %2 %3  %4  %5")
                         .arg(atom0,6).arg(atom1,6).arg(param.functionType(),6)
                         .arg(params0.join("  "))
                         .arg(params1.join("  ")) );
            }
        }
        else
        {
            // Get the bonds from the molecule.
            const auto &bonds = moltype.bonds();

            for (auto it = bonds.constBegin(); it != bonds.constEnd(); ++it)
            {
                const auto &bond = it.key();
                const auto &param = it.value();

                //AtomID is AtomIdx. Add 1, as gromacs is 1-indexed
                int atom0 = bond.atom0().asA<AtomIdx>().value() + 1;
                int atom1 = bond.atom1().asA<AtomIdx>().value() + 1;

                QStringList params;
                for (const auto p : param.parameters())
                    params.append( QString::number(p) );

                bondlines.append( QString("%1 %2 %3  %4")
                         .arg(atom0,6).arg(atom1,6).arg(param.functionType(),6)
                         .arg(params.join("  ")) );
            }
        }

        qSort(bondlines);
    };

    //write all of the angles
    auto write_angs = [&]()
    {
        if (is_perturbable)
        {
            // Get the angles from the molecule.
            const auto &angles0 = moltype.angles();
            const auto &angles1 = moltype.angles(true);

            // Sets to contain the AngleIDs at lambda = 0 and lambda = 1.
            QSet<AngleID> angles0_idx;
            QSet<AngleID> angles1_idx;

            // Loop over all angles at lambda = 0.
            for (const auto &idx : angles0.uniqueKeys())
                angles0_idx.insert(idx);

            // Loop over all angles at lambda = 1.
            for (const auto &idx : angles1.uniqueKeys())
            {
                if (angles0_idx.contains(idx.mirror()))
                    angles1_idx.insert(idx.mirror());
                else
                    angles1_idx.insert(idx);
            }

            // Now work out the AngleIDs that are unique at lambda = 0 and lambda = 1,
            // as well as those that are shared.
            QSet<AngleID> angles0_uniq_idx;
            QSet<AngleID> angles1_uniq_idx;
            QSet<AngleID> angles_shared_idx;

            // lambda = 0
            for (const auto &idx : angles0_idx)
            {
                if (not angles1_idx.contains(idx))
                    angles0_uniq_idx.insert(idx);
                else
                    angles_shared_idx.insert(idx);
            }

            // lambda = 1
            for (const auto &idx : angles1_idx)
            {
                if (not angles0_idx.contains(idx))
                    angles1_uniq_idx.insert(idx);
                else
                    angles_shared_idx.insert(idx);
            }

            // First create parameter records for the angles unique to lambda = 0/1.

            // lambda = 0
            for (const auto &idx : angles0_uniq_idx)
            {
                //AtomID is AtomIdx. Add 1, as gromacs is 1-indexed
                int atom0 = idx.atom0().asA<AtomIdx>().value() + 1;
                int atom1 = idx.atom1().asA<AtomIdx>().value() + 1;
                int atom2 = idx.atom2().asA<AtomIdx>().value() + 1;

                // Get all of the parameters for this AngleID.
                const auto &params = angles0.values(idx);

                // Loop over all of the parameters.
                for (const auto &param : params)
                {
                    QStringList param_string;
                    for (const auto p : param.parameters())
                        param_string.append( QString::number(p) );

                    anglines.append( QString("%1 %2 %3 %4   %5  0  0")
                            .arg(atom0,6).arg(atom1,6).arg(atom2,6).arg(param.functionType(),7)
                            .arg(param_string.join("  ")) );
                }
            }

            // lambda = 1
            for (const auto &idx : angles1_uniq_idx)
            {
                //AtomID is AtomIdx. Add 1, as gromacs is 1-indexed
                int atom0 = idx.atom0().asA<AtomIdx>().value() + 1;
                int atom1 = idx.atom1().asA<AtomIdx>().value() + 1;
                int atom2 = idx.atom2().asA<AtomIdx>().value() + 1;

                // Get all of the parameters for this AngleID.
                const auto &params = angles1.values(idx);

                // Loop over all of the parameters.
                for (const auto &param : params)
                {
                    QStringList param_string;
                    for (const auto p : param.parameters())
                        param_string.append( QString::number(p) );

                    anglines.append( QString("%1 %2 %3 %4   0  0  %5")
                            .arg(atom0,6).arg(atom1,6).arg(atom2,6).arg(param.functionType(),7)
                            .arg(param_string.join("  ")) );
                }
            }

            // Next add the shared angle parameters.

            for (auto idx : angles_shared_idx)
            {
                //AtomID is AtomIdx. Add 1, as gromacs is 1-indexed
                int atom0 = idx.atom0().asA<AtomIdx>().value() + 1;
                int atom1 = idx.atom1().asA<AtomIdx>().value() + 1;
                int atom2 = idx.atom2().asA<AtomIdx>().value() + 1;

                // Get a list of the parameters at lambda = 0.
                const auto &params0 = angles0.values(idx);

                // Invert the index.
                if (not angles1.contains(idx))
                    idx = idx.mirror();

                // Get a list of the parameters at lambda = 1.
                const auto &params1 = angles1.values(idx);

                // More or same number of records at lambda = 0.
                if (params0.count() >= params1.count())
                {
                    for (int i=0; i<params1.count(); ++i)
                    {
                        QStringList param_string0;
                        for (const auto p : params0[i].parameters())
                            param_string0.append( QString::number(p) );

                        QStringList param_string1;
                        for (const auto p : params1[i].parameters())
                            param_string1.append( QString::number(p) );

                        anglines.append( QString("%1 %2 %3 %4   %5  %6")
                                .arg(atom0,6).arg(atom1,6).arg(atom2,6).arg(params0[i].functionType(),7)
                                .arg(param_string0.join("  "))
                                .arg(param_string1.join("  ")) );
                    }

                    // Now add parameters for which there is no matching record
                    // at lambda = 1.
                    for (int i=params1.count(); i<params0.count(); ++i)
                    {
                        QStringList param_string;
                        for (const auto p : params0[i].parameters())
                            param_string.append( QString::number(p) );

                        anglines.append( QString("%1 %2 %3 %4   %5  0  0")
                                .arg(atom0,6).arg(atom1,6).arg(atom2,6).arg(params0[i].functionType(),7)
                                .arg(param_string.join("  ")) );
                    }
                }

                // More records at lambda = 1.
                else
                {
                    for (int i=0; i<params0.count(); ++i)
                    {
                        QStringList param_string0;
                        for (const auto p : params0[i].parameters())
                            param_string0.append( QString::number(p) );

                        QStringList param_string1;
                        for (const auto p : params1[i].parameters())
                            param_string1.append( QString::number(p) );

                        anglines.append( QString("%1 %2 %3 %4   %5  %6")
                                .arg(atom0,6).arg(atom1,6).arg(atom2,6).arg(params1[i].functionType(),7)
                                .arg(param_string0.join("  "))
                                .arg(param_string1.join("  ")) );
                    }

                    // Now add parameters for which there is no matching record
                    // at lambda = 0.
                    for (int i=params0.count(); i<params1.count(); ++i)
                    {
                        QStringList param_string;
                        for (const auto p : params1[i].parameters())
                            param_string.append( QString::number(p) );

                        anglines.append( QString("%1 %2 %3 %4   0  0  %5")
                                .arg(atom0,6).arg(atom1,6).arg(atom2,6).arg(params1[i].functionType(),7)
                                .arg(param_string.join("  ")) );
                    }
                }
            }
        }
        else
        {
            // Get the angles from the molecule.
            const auto &angles = moltype.angles();

            for (auto it = angles.constBegin(); it != angles.constEnd(); ++it)
            {
                const auto &angle = it.key();
                const auto &param = it.value();

                //AtomID is AtomIdx. Add 1, as gromacs is 1-indexed
                int atom0 = angle.atom0().asA<AtomIdx>().value() + 1;
                int atom1 = angle.atom1().asA<AtomIdx>().value() + 1;
                int atom2 = angle.atom2().asA<AtomIdx>().value() + 1;

                QStringList params;
                for (const auto p : param.parameters())
                    params.append( QString::number(p) );

                anglines.append( QString("%1 %2 %3 %4   %5")
                        .arg(atom0,6).arg(atom1,6).arg(atom2,6).arg(param.functionType(),7)
                        .arg(params.join("  ")) );
            }
        }

        qSort(anglines);
    };

    //write all of the dihedrals/impropers (they are merged)
    auto write_dihs = [&]()
    {
        if (is_perturbable)
        {
            // Get the dihedrals from the molecule.
            const auto &dihedrals0 = moltype.dihedrals();
            const auto &dihedrals1 = moltype.dihedrals(true);

            // Sets to contain the DihedralID at lambda = 0 and lambda = 1.
            QSet<DihedralID> dihedrals0_idx;
            QSet<DihedralID> dihedrals1_idx;

            // Loop over all dihedrals at lambda = 0.
            for (const auto &idx : dihedrals0.uniqueKeys())
                dihedrals0_idx.insert(idx);

            // Loop over all dihedrals at lambda = 1.
            for (const auto &idx : dihedrals1.uniqueKeys())
            {
                if (dihedrals0_idx.contains(idx.mirror()))
                    dihedrals1_idx.insert(idx.mirror());
                else
                    dihedrals1_idx.insert(idx);
            }

            // Now work out the DihedralIDs that are unique at lambda = 0 and lambda = 1,
            // as well as those that are shared.
            QSet<DihedralID> dihedrals0_uniq_idx;
            QSet<DihedralID> dihedrals1_uniq_idx;
            QSet<DihedralID> dihedrals_shared_idx;

            // lambda = 0
            for (const auto &idx : dihedrals0_idx)
            {
                if (not dihedrals1_idx.contains(idx))
                    dihedrals0_uniq_idx.insert(idx);
                else
                    dihedrals_shared_idx.insert(idx);
            }

            // lambda = 1
            for (const auto &idx : dihedrals1_idx)
            {
                if (not dihedrals0_idx.contains(idx))
                    dihedrals1_uniq_idx.insert(idx);
                else
                    dihedrals_shared_idx.insert(idx);
            }

            // First create parameter records for the dihedrals unique to lambda = 0/1.

            // lambda = 0
            for (const auto &idx : dihedrals0_uniq_idx)
            {
                //AtomID is AtomIdx. Add 1, as gromacs is 1-indexed
                int atom0 = idx.atom0().asA<AtomIdx>().value() + 1;
                int atom1 = idx.atom1().asA<AtomIdx>().value() + 1;
                int atom2 = idx.atom2().asA<AtomIdx>().value() + 1;
                int atom3 = idx.atom3().asA<AtomIdx>().value() + 1;

                // Get all of the parameters for this DihedralID.
                const auto &params = dihedrals0.values(idx);

                // Loop over all of the parameters.
                for (const auto &param : params)
                {
                    QStringList param_string;
                    for (const auto p : param.parameters())
                        param_string.append( QString::number(p) );

                    dihlines.append( QString("%1 %2 %3 %4 %5  %6  0  0  0")
                            .arg(atom0,6).arg(atom1,6)
                            .arg(atom2,6).arg(atom3,6).arg(param.functionType(),6)
                            .arg(param_string.join("  ")) );
                }
            }

            // lambda = 1
            for (const auto &idx : dihedrals1_uniq_idx)
            {
                //AtomID is AtomIdx. Add 1, as gromacs is 1-indexed
                int atom0 = idx.atom0().asA<AtomIdx>().value() + 1;
                int atom1 = idx.atom1().asA<AtomIdx>().value() + 1;
                int atom2 = idx.atom2().asA<AtomIdx>().value() + 1;
                int atom3 = idx.atom3().asA<AtomIdx>().value() + 1;

                // Get all of the parameters for this AngleID.
                const auto &params = dihedrals1.values(idx);

                // Loop over all of the parameters.
                for (const auto &param : params)
                {
                    QStringList param_string;
                    for (const auto p : param.parameters())
                        param_string.append( QString::number(p) );

                    dihlines.append( QString("%1 %2 %3 %4 %5  0  0  0  %6")
                            .arg(atom0,6).arg(atom1,6)
                            .arg(atom2,6).arg(atom3,6).arg(param.functionType(),6)
                            .arg(param_string.join("  ")) );
                }
            }

            // Next add the shared angle parameters.

            for (auto idx : dihedrals_shared_idx)
            {
                //AtomID is AtomIdx. Add 1, as gromacs is 1-indexed
                int atom0 = idx.atom0().asA<AtomIdx>().value() + 1;
                int atom1 = idx.atom1().asA<AtomIdx>().value() + 1;
                int atom2 = idx.atom2().asA<AtomIdx>().value() + 1;
                int atom3 = idx.atom3().asA<AtomIdx>().value() + 1;

                // Get a list of the parameters at lambda = 0.
                const auto &params0 = dihedrals0.values(idx);

                // Invert the index.
                if (not dihedrals1.contains(idx))
                    idx = idx.mirror();

                // Get a list of the parameters at lambda = 1.
                const auto &params1 = dihedrals1.values(idx);

                // More or same number of records at lambda = 0.
                if (params0.count() >= params1.count())
                {
                    for (int i=0; i<params1.count(); ++i)
                    {
                        QStringList param_string0;
                        for (const auto p : params0[i].parameters())
                            param_string0.append( QString::number(p) );

                        QStringList param_string1;
                        for (const auto p : params1[i].parameters())
                            param_string1.append( QString::number(p) );

                        dihlines.append( QString("%1 %2 %3 %4 %5  %6  %7")
                                .arg(atom0,6).arg(atom1,6)
                                .arg(atom2,6).arg(atom3,6).arg(params0[i].functionType(),6)
                                .arg(param_string0.join("  "))
                                .arg(param_string1.join("  ")) );
                    }

                    // Now add parameters for which there is no matching record
                    // at lambda = 1.
                    for (int i=params1.count(); i<params0.count(); ++i)
                    {
                        QStringList param_string;
                        for (const auto p : params0[i].parameters())
                            param_string.append( QString::number(p) );

                        dihlines.append( QString("%1 %2 %3 %4 %5  %6  0  0  0")
                                .arg(atom0,6).arg(atom1,6)
                                .arg(atom2,6).arg(atom3,6).arg(params0[i].functionType(),6)
                                .arg(param_string.join("  ")) );
                    }
                }

                // More records at lambda = 1.
                else
                {
                    for (int i=0; i<params0.count(); ++i)
                    {
                        QStringList param_string0;
                        for (const auto p : params0[i].parameters())
                            param_string0.append( QString::number(p) );

                        QStringList param_string1;
                        for (const auto p : params1[i].parameters())
                            param_string1.append( QString::number(p) );

                        dihlines.append( QString("%1 %2 %3 %4 %5  %6  %7")
                                .arg(atom0,6).arg(atom1,6)
                                .arg(atom2,6).arg(atom3,6).arg(params1[i].functionType(),6)
                                .arg(param_string0.join("  "))
                                .arg(param_string1.join("  ")) );
                    }

                    // Now add parameters for which there is no matching record
                    // at lambda = 0.
                    for (int i=params0.count(); i<params1.count(); ++i)
                    {
                        QStringList param_string;
                        for (const auto p : params1[i].parameters())
                            param_string.append( QString::number(p) );

                        dihlines.append( QString("%1 %2 %3 %4 %5  0  0  0  %6")
                                .arg(atom0,6).arg(atom1,6)
                                .arg(atom2,6).arg(atom3,6).arg(params1[i].functionType(),6)
                                .arg(param_string.join("  ")) );
                    }
                }
            }
        }
        else
        {
            // Get the dihedrals from the molecule.
            const auto &dihedrals = moltype.dihedrals();

            for (auto it = dihedrals.constBegin(); it != dihedrals.constEnd();
                ++it)
            {
                const auto &dihedral = it.key();
                const auto &param = it.value();

                //AtomID is AtomIdx. Add 1, as gromacs is 1-indexed
                int atom0 = dihedral.atom0().asA<AtomIdx>().value() + 1;
                int atom1 = dihedral.atom1().asA<AtomIdx>().value() + 1;
                int atom2 = dihedral.atom2().asA<AtomIdx>().value() + 1;
                int atom3 = dihedral.atom3().asA<AtomIdx>().value() + 1;

                QStringList params;
                for (const auto p : param.parameters())
                    params.append( QString::number(p) );

                dihlines.append( QString("%1 %2 %3 %4 %5  %6")
                        .arg(atom0,6).arg(atom1,6)
                        .arg(atom2,6).arg(atom3,6).arg(param.functionType(),6)
                        .arg(params.join("  ")) );
            }
        }

        qSort(dihlines);
    };

    //write all of the pairs (1-4 scaling factors). This is needed even though
    //we have set autogenerate pairs to "yes"
    auto write_pairs = [&]()
    {
        if (is_perturbable)
        {
            CLJNBPairs scl0;
            CLJNBPairs scl1;

            try
            {
                scl0 = mol.property("intrascale0").asA<CLJNBPairs>();
            }
            catch(...)
            {
                return;
            }

            try
            {
                scl1 = mol.property("intrascale1").asA<CLJNBPairs>();
            }
            catch(...)
            {
                return;
            }

            // must record every pair that has a non-default scaling factor
            for (int i=0; i<scl0.nAtoms()-1; ++i)
            {
                for (int j=i+1; j<scl0.nAtoms(); ++j)
                {
                    const auto s0 = scl0.get( AtomIdx(i), AtomIdx(j) );
                    const auto s1 = scl1.get( AtomIdx(i), AtomIdx(j) );

                    if ( not ( (s0.coulomb() == 0 and s0.lj() == 0 and s1.coulomb() == 0 and s1.lj() == 0) or
                               (s0.coulomb() == 1 and s0.lj() == 1 and s1.coulomb() == 1 and s1.lj() == 1) ) )
                    {
                        //this is a non-default pair
                        scllines.append( QString("%1 %2     1")
                                            .arg(i+1,6).arg(j+1,6) );
                    }
                }
            }
        }
        else
        {
            CLJNBPairs scl;

            try
            {
                scl = mol.property("intrascale").asA<CLJNBPairs>();
            }
            catch(...)
            {
                return;
            }

            // must record every pair that has a non-default scaling factor
            for (int i=0; i<scl.nAtoms()-1; ++i)
            {
                for (int j=i+1; j<scl.nAtoms(); ++j)
                {
                    const auto s = scl.get( AtomIdx(i), AtomIdx(j) );

                    if ( not ( (s.coulomb() == 0 and s.lj() == 0) or
                            (s.coulomb() == 1 and s.lj() == 1) ) )
                    {
                        //this is a non-default pair
                        scllines.append( QString("%1 %2     1")
                                            .arg(i+1,6).arg(j+1,6) );
                    }
                }
            }
        }
    };

    const QVector< std::function<void()> > funcs =
                 { write_atoms, write_bonds, write_angs, write_dihs, write_pairs };

    if (uses_parallel)
    {
        tbb::parallel_for( tbb::blocked_range<int>(0, funcs.count(), 1),
                           [&](const tbb::blocked_range<int> &r)
        {
            for (int i=r.begin(); i<r.end(); ++i)
            {
                funcs[i]();
            }
        });
    }
    else
    {
        for (int i=0; i<funcs.count(); ++i)
        {
            funcs[i]();
        }
    }

    lines.append( "[ atoms ]" );
    if (is_perturbable)
        lines.append(";   nr  type0  resnr residue  atom   cgnr    charge0        mass0  type1    charge1        mass1");
    else
        lines.append(";   nr   type  resnr residue  atom   cgnr     charge         mass");
    lines.append(atomlines);

    // we need to detect whether this is a water molecule. If so, then we
    // need to add in "settles" lines to constrain bonds / angles of the
    // water molecule
    const bool is_water = moltype.isWater();

    if (is_water)
    {
        lines.append( "#ifdef FLEXIBLE" );
    }

    if (not bondlines.isEmpty())
    {
        lines.append( "[ bonds ]" );
        lines.append( ";   ai     aj  funct  parameters" );
        lines += bondlines;
        lines.append("");
    }

    if (not scllines.isEmpty())
    {
        lines.append( "[ pairs ]" );
        lines.append( ";  ai    aj funct " );
        lines += scllines;
        lines.append("");
    }

    if (not anglines.isEmpty())
    {
        lines.append( "[ angles ]" );
        lines.append( ";   ai     aj     ak   funct   parameters" );
        lines += anglines;
        lines.append("");
    }

    if (not dihlines.isEmpty())
    {
        lines.append( "[ dihedrals ]" );
        lines.append( ";   ai     aj     ak     al  funct  parameters" );
        lines += dihlines;
        lines.append("");
    }

    if (is_water)
    {
        lines.append("#else");
        lines.append("");
        lines += moltype.settlesLines();
        lines.append("");
        lines.append("#endif");
    }

    return lines;
}

/** Internal function used to convert an array of Gromacs Moltyps into
    lines of a Gromacs topology file */
static QStringList writeMolTypes(const QHash<QString,GroMolType> &moltyps,
                                 const QHash<QString,Molecule> &examples,
                                 bool uses_parallel)
{
    QHash<QString,QStringList> typs;

    if (uses_parallel)
    {
        const QVector<QString> keys = moltyps.keys().toVector();
        QMutex mutex;

        tbb::parallel_for( tbb::blocked_range<int>(0,keys.count(),1),
                           [&](const tbb::blocked_range<int> &r)
        {
            for (int i=r.begin(); i<r.end(); ++i)
            {
                QStringList typlines = ::writeMolType(keys[i], moltyps[keys[i]],
                                                      examples[keys[i]], uses_parallel);

                QMutexLocker lkr(&mutex);
                typs.insert(keys[i], typlines);
            }
        });
    }
    else
    {
        for (auto it = moltyps.constBegin(); it != moltyps.constEnd(); ++it)
        {
            typs.insert( it.key(), ::writeMolType(it.key(), it.value(),
                                                  examples[it.key()], uses_parallel) );
        }
    }

    QStringList keys = moltyps.keys();
    keys.sort();

    QStringList lines;

    for (const auto key : keys)
    {
        lines += typs[key];
        lines += "";
    }

    return lines;
}

/** Internal function used to write the system part of the gromacs file */
static QStringList writeSystem(QString name, const QStringList &mol_to_moltype)
{
    QStringList lines;
    lines.append( "[ system ]" );
    lines.append( name );
    lines.append( "" );
    lines.append( "[ molecules ]" );
    lines.append( ";molecule name    nr." );

    QString lastmol;

    int count = 0;

    for (auto it = mol_to_moltype.constBegin(); it != mol_to_moltype.constEnd(); ++it)
    {
        if (*it != lastmol)
        {
            if (lastmol.isNull())
            {
                lastmol = *it;
                count = 1;
            }
            else
            {
                lines.append( QString("%1 %2").arg(lastmol,14).arg(count,6) );
                lastmol = *it;
                count = 1;
            }
        }
        else
            count += 1;
    }

    lines.append( QString("%1 %2").arg(lastmol,14).arg(count,6) );

    lines.append("");

    return lines;
}

/** Construct this parser by extracting all necessary information from the
    passed SireSystem::System, looking for the properties that are specified
    in the passed property map */
GroTop::GroTop(const SireSystem::System &system, const PropertyMap &map)
       : ConcreteProperty<GroTop,MoleculeParser>(map),
         nb_func_type(0), combining_rule(0), fudge_lj(0), fudge_qq(0),
         generate_pairs(false)
{
    //get the MolNums of each molecule in the System - this returns the
    //numbers in MolIdx order
    const QVector<MolNum> molnums = system.getMoleculeNumbers().toVector();

    if (molnums.isEmpty())
    {
        //no molecules in the system
        this->operator=(GroTop());
        return;
    }

    //generate the GroMolType object for each molecule in the system. It is likely
    //the multiple molecules will have the same GroMolType. We will be a bit slow
    //now and generate this for all molecules independently, and will then consolidate
    //them later
    QVector<GroMolType> mtyps(molnums.count());

    if (usesParallel())
    {
        tbb::parallel_for( tbb::blocked_range<int>(0,molnums.count()),
                           [&](const tbb::blocked_range<int> r)
        {
            for (int i=r.begin(); i<r.end(); ++i)
            {
                mtyps[i] = GroMolType(system[molnums[i]].molecule(),map);
            }
        });
    }
    else
    {
        for (int i=0; i<molnums.count(); ++i)
        {
            mtyps[i] = GroMolType(system[molnums[i]].molecule(),map);
        }
    }

    //now go through each type, remove duplicates, and record the name of
    //each molecule so that we can write this in the [system] section
    QStringList mol_to_moltype;
    QHash<QString,GroMolType> name_to_mtyp;
    QHash<QString,Molecule> name_to_example;

    for (int i=0; i<mtyps.count(); ++i)
    {
        const auto moltype = mtyps[i];
        QString name = moltype.name();

        if (name_to_mtyp.contains(name))
        {
            if (moltype != name_to_mtyp[name])
            {
                //this has the same name but different details. Give this a new
                //name
                int j = 0;

                while(true)
                {
                    j++;
                    name = QString("%1_%2").arg(moltype.name()).arg(j);

                    if (name_to_mtyp.contains(name))
                    {
                        if (moltype == name_to_mtyp[name])
                            //match :-)
                            break;
                    }
                    else
                    {
                        //new moltype
                        name_to_mtyp.insert(name, moltype);

                        //save an example of this molecule so that we can
                        //extract any other details necessary
                        name_to_example.insert(name, system[molnums[i]].molecule());

                        break;
                    }

                    //we have got here, meaning that we need to try a different name
                }
            }
        }
        else
        {
            name_to_mtyp.insert(name, moltype);
            name_to_example.insert(name, system[molnums[i]].molecule());
        }

        mol_to_moltype.append(name);
    }

    QStringList errors;

    //first, we need to extract the common forcefield from the molecules
    MMDetail ffield = name_to_mtyp.constBegin()->forcefield();

    for (auto it = name_to_mtyp.constBegin(); it != name_to_mtyp.constEnd(); ++it)
    {
        if (not ffield.isCompatibleWith(it.value().forcefield()))
        {
            errors.append( QObject::tr( "The forcefield for molecule '%1' is not "
              "compatible with that for other molecules.\n%1 versus\n%2")
                .arg(it.key()).arg(it.value().forcefield().toString()).arg(ffield.toString()) );
        }
    }

    if (not errors.isEmpty())
    {
        throw SireError::incompatible_error( QObject::tr(
            "Cannot write this system to a Gromacs Top file as the forcefields of the "
            "molecules are incompatible with one another.\n%1")
                .arg(errors.join("\n\n")), CODELOC );
    }

    //next, we need to write the defaults section of the file
    QStringList lines = ::writeDefaults(ffield);

    //next, we need to extract and write all of the atom types from all of
    //the molecules
    lines += ::writeAtomTypes(name_to_mtyp, name_to_example, ffield, map);

    //now convert these into text lines that can be written as the file
    lines += ::writeMolTypes(name_to_mtyp, name_to_example, usesParallel());

    //now write the system part
    lines += ::writeSystem(system.name(), mol_to_moltype);

    if (not errors.isEmpty())
    {
        throw SireIO::parse_error( QObject::tr(
            "Errors converting the system to a Gromacs Top format...\n%1")
                .arg(lines.join("\n")), CODELOC );
    }

    //we don't need params any more, so free the memory
    mtyps.clear();
    name_to_mtyp.clear();
    mol_to_moltype.clear();

    //now we have the lines, reparse them to make sure that they are correct
    //and we have a fully-constructed and sane GroTop object
    GroTop parsed( lines, map );

    this->operator=(parsed);
}

/** Copy constructor */
GroTop::GroTop(const GroTop &other)
       : ConcreteProperty<GroTop,MoleculeParser>(other),
         include_path(other.include_path), included_files(other.included_files),
         expanded_lines(other.expanded_lines),
         atom_types(other.atom_types),
         bond_potentials(other.bond_potentials),
         ang_potentials(other.ang_potentials),
         dih_potentials(other.dih_potentials),
         moltypes(other.moltypes), grosys(other.grosys),
         nb_func_type(other.nb_func_type), combining_rule(other.combining_rule),
         fudge_lj(other.fudge_lj), fudge_qq(other.fudge_qq),
         parse_warnings(other.parse_warnings),
         generate_pairs(other.generate_pairs)
{}

/** Destructor */
GroTop::~GroTop()
{}

/** Copy assignment operator */
GroTop& GroTop::operator=(const GroTop &other)
{
    if (this != &other)
    {
        include_path = other.include_path;
        included_files = other.included_files;
        expanded_lines = other.expanded_lines;
        atom_types = other.atom_types;
        bond_potentials = other.bond_potentials;
        ang_potentials = other.ang_potentials;
        dih_potentials = other.dih_potentials;
        moltypes = other.moltypes;
        grosys = other.grosys;
        nb_func_type = other.nb_func_type;
        combining_rule = other.combining_rule;
        fudge_lj = other.fudge_lj;
        fudge_qq = other.fudge_qq;
        parse_warnings = other.parse_warnings;
        generate_pairs = other.generate_pairs;
        MoleculeParser::operator=(other);
    }

    return *this;
}

/** Comparison operator */
bool GroTop::operator==(const GroTop &other) const
{
    return include_path == other.include_path and
           included_files == other.included_files and
           expanded_lines == other.expanded_lines and
           MoleculeParser::operator==(other);
}

/** Comparison operator */
bool GroTop::operator!=(const GroTop &other) const
{
    return not operator==(other);
}

/** Return the C++ name for this class */
const char* GroTop::typeName()
{
    return QMetaType::typeName( qMetaTypeId<GroTop>() );
}

/** Return the C++ name for this class */
const char* GroTop::what() const
{
    return GroTop::typeName();
}

/** Return whether or not this is a lead parser. The lead parser is responsible
    for starting the process of turning the parsed file into the System. There
    must be one and one-only lead parser in a set of parsers creating a System */
bool GroTop::isLead() const
{
    return true;
}

/** Return whether or not this parser can follow another lead parser, and add
    data to an existing molecular system. The GroTop parser cannot follow. */
bool GroTop::canFollow() const
{
    return false;
}

/** Return the list of names of directories in which to search for
    include files. The directories are either absolute, or relative
    to the current directory. If "absolute_paths" is true then
    the full absolute paths for directories that exist on this
    machine will be returned */
QStringList GroTop::includePath(bool absolute_paths) const
{
    if (absolute_paths)
    {
        QStringList abspaths;

        for (const auto path : include_path)
        {
            QFileInfo file(path);

            if (file.exists())
                abspaths.append( file.absoluteFilePath() );
        }

        return abspaths;
    }
    else
        return include_path;
}

/** Return the list of names of files that were included when reading or
    writing this file. The files are relative. If "absolute_paths"
    is true then the full absolute paths for the files will be
    used */
QStringList GroTop::includedFiles(bool absolute_paths) const
{
    //first, go through the list of included files
    QStringList files;

    for (auto it = included_files.constBegin(); it != included_files.constEnd(); ++it)
    {
        files += it.value();
    }

    if (absolute_paths)
    {
        //these are already absolute filenames
        return files;
    }
    else
    {
        //subtract any paths that relate to the current directory or GROMACS_PATH
        QString curpath = QDir::current().absolutePath();

        for (auto it = files.begin(); it != files.end(); ++it)
        {
            if ( it->startsWith(curpath) )
            {
                *it = it->mid(curpath.length()+1);
            }
            else
            {
                for (const auto path : include_path)
                {
                    if (it->startsWith(path))
                    {
                        *it = it->mid(path.length()+1);
                    }
                }
            }
        }

        return files;
    }
}

/** Return the parser that has been constructed by reading in the passed
    file using the passed properties */
MoleculeParserPtr GroTop::construct(const QString &filename,
                                  const PropertyMap &map) const
{
    return GroTop(filename,map);
}

/** Return the parser that has been constructed by reading in the passed
    text lines using the passed properties */
MoleculeParserPtr GroTop::construct(const QStringList &lines,
                                  const PropertyMap &map) const
{
    return GroTop(lines,map);
}

/** Return the parser that has been constructed by extract all necessary
    data from the passed SireSystem::System using the specified properties */
MoleculeParserPtr GroTop::construct(const SireSystem::System &system,
                                    const PropertyMap &map) const
{
    return GroTop(system,map);
}

/** Return a string representation of this parser */
QString GroTop::toString() const
{
    return QObject::tr("GroTop( includePath() = [%1], includedFiles() = [%2] )")
                .arg(includePath().join(", "))
                .arg(includedFiles().join(", "));
}

/** Return the format name that is used to identify this file format within Sire */
QString GroTop::formatName() const
{
    return "GroTop";
}

/** Return a description of the file format */
QString GroTop::formatDescription() const
{
    return QObject::tr("Gromacs Topology format files.");
}

/** Return the suffixes that these files are normally associated with */
QStringList GroTop::formatSuffix() const
{
    static const QStringList suffixes = { "top" };
    return suffixes;
}

/** Function that is called to assert that this object is sane. This
    should raise an exception if the parser is in an invalid state */
void GroTop::assertSane() const
{
    //check state, raise SireError::program_bug if we are in an invalid state
}

/** Return the atom type data for the passed atom type. This returns
    null data if it is not present */
GromacsAtomType GroTop::atomType(const QString &atm) const
{
    return atom_types.value(atm, GromacsAtomType());
}

/** Return the ID string for the bond atom types 'atm0' 'atm1'. This
    creates the string 'atm0;atm1' or 'atm1;atm0' depending on which
    of the atoms is lower. The ';' character is used as a separator
    as it cannot be in the atom names, as it is used as a comment
    character in the Gromacs Top file */
static QString get_bond_id(const QString &atm0, const QString &atm1, int func_type)
{
    if (func_type == 0) // default type
        func_type = 1;

    if (atm0 < atm1)
    {
        return QString("%1;%2;%3").arg(atm0,atm1).arg(func_type);
    }
    else
    {
        return QString("%1;%2;%3").arg(atm1,atm0).arg(func_type);
    }
}

/** Return the ID string for the angle atom types 'atm0' 'atm1' 'atm2'. This
    creates the string 'atm0;atm1;atm2' or 'atm2;atm1;atm0' depending on which
    of the atoms is lower. The ';' character is used as a separator
    as it cannot be in the atom names, as it is used as a comment
    character in the Gromacs Top file */
static QString get_angle_id(const QString &atm0, const QString &atm1, const QString &atm2,
                            int func_type)
{
    if (func_type == 0)
        func_type = 1;  //default type

    if (atm0 < atm2)
    {
        return QString("%1;%2;%3;%4").arg(atm0,atm1,atm2).arg(func_type);
    }
    else
    {
        return QString("%1;%2;%3;%4").arg(atm2,atm1,atm0).arg(func_type);
    }
}

/** Return the ID string for the dihedral atom types 'atm0' 'atm1' 'atm2' 'atm3'. This
    creates the string 'atm0;atm1;atm2;atm3' or 'atm3;atm2;atm1;atm0' depending on which
    of the atoms is lower. The ';' character is used as a separator
    as it cannot be in the atom names, as it is used as a comment
    character in the Gromacs Top file */
static QString get_dihedral_id(const QString &atm0, const QString &atm1,
                               const QString &atm2, const QString &atm3,
                               int func_type)
{
    if (atm0 < atm3)
    {
        return QString("%1;%2;%3;%4;%5").arg(atm0,atm1,atm2,atm3).arg(func_type);
    }
    else
    {
        return QString("%1;%2;%3;%4;%5").arg(atm3,atm2,atm1,atm0).arg(func_type);
    }
}

/** Return the Gromacs System that describes the list of molecules that should
    be contained */
GroSystem GroTop::groSystem() const
{
    return grosys;
}

/** Return the bond potential data for the passed pair of atoms. This only returns
    the most recently inserted parameter for this pair. Use 'bonds' if you want
    to allow for multiple return values */
GromacsBond GroTop::bond(const QString &atm0, const QString &atm1,
                         int func_type) const
{
    return bond_potentials.value( get_bond_id(atm0,atm1,func_type), GromacsBond() );
}

/** Return the bond potential data for the passed pair of atoms. This returns
    a list of all associated parameters */
QList<GromacsBond> GroTop::bonds(const QString &atm0, const QString &atm1, int func_type) const
{
    return bond_potentials.values( get_bond_id(atm0,atm1,func_type) );
}

/** Return the angle potential data for the passed triple of atoms. This only returns
    the most recently inserted parameter for these atoms. Use 'angles' if you want
    to allow for multiple return values */
GromacsAngle GroTop::angle(const QString &atm0, const QString &atm1, const QString &atm2,
                           int func_type) const
{
    return ang_potentials.value( get_angle_id(atm0,atm1,atm2,func_type), GromacsAngle() );
}

/** Return the angle potential data for the passed triple of atoms. This returns
    a list of all associated parameters */
QList<GromacsAngle> GroTop::angles(const QString &atm0, const QString &atm1,
                                   const QString &atm2, int func_type) const
{
    return ang_potentials.values( get_angle_id(atm0,atm1,atm2,func_type) );
}

/** Search for a dihedral type parameter that matches the atom types
    atom0-atom1-atom2-atom3. This will try to find an exact match. If that fails,
    it will then use one of the wildcard matches. Returns a null string if there
    is no match. This will return the key into the dih_potentials dictionary */
QString GroTop::searchForDihType(const QString &atm0, const QString &atm1,
                                 const QString &atm2, const QString &atm3,
                                 int func_type) const
{
    QString key = get_dihedral_id(atm0,atm1,atm2,atm3,func_type);

    //qDebug() << "SEARCHING FOR" << key;

    if (dih_potentials.contains(key))
    {
        //qDebug() << "FOUND" << key;
        return key;
    }

    static const QString wild = "X";

    //look for *-atm1-atm2-atm3
    key = get_dihedral_id(wild, atm1, atm2, atm3, func_type);

    if (dih_potentials.contains(key))
    {
        //qDebug() << "FOUND" << key;
        return key;
    }

    //look for *-atm2-atm1-atm0
    key = get_dihedral_id(wild, atm2, atm1, atm0, func_type);

    if (dih_potentials.contains(key))
    {
        //qDebug() << "FOUND" << key;
        return key;
    }

    //this failed. Look for *-atm1-atm2-* or *-atm2-atm1-*
    key = get_dihedral_id(wild, atm1, atm2, wild, func_type);

    if (dih_potentials.contains(key))
    {
        //qDebug() << "FOUND" << key;
        return key;
    }

    key = get_dihedral_id(wild, atm2, atm1, wild, func_type);

    if (dih_potentials.contains(key))
    {
        //qDebug() << "FOUND" << key;
        return key;
    }

    //look for *-*-atm2-atm3
    key = get_dihedral_id(wild, wild, atm2, atm3, func_type);

    if (dih_potentials.contains(key))
    {
        //qDebug() << "FOUND" << key;
        return key;
    }

    //look for *-*-atm1-atm0
    key = get_dihedral_id(wild, wild, atm1, atm0, func_type);

    if (dih_potentials.contains(key))
    {
        //qDebug() << "FOUND" << key;
        return key;
    }

    //finally look for *-*-*-*
    key = get_dihedral_id(wild, wild, wild, wild, func_type);

    if (dih_potentials.contains(key))
    {
        //qDebug() << "FOUND" << key;
        return key;
    }

    return QString();
}

/** Return the dihedral potential data for the passed quad of atoms. This only returns
    the most recently inserted parameter for these atoms. Use 'dihedrals' if you want
    to allow for multiple return values */
GromacsDihedral GroTop::dihedral(const QString &atm0, const QString &atm1,
                                 const QString &atm2, const QString &atm3,
                                 int func_type) const
{
    return dih_potentials.value( searchForDihType(atm0,atm1,atm2,atm3,func_type),
                                 GromacsDihedral() );
}

/** Return the dihedral potential data for the passed quad of atoms. This returns
    a list of all associated parameters */
QList<GromacsDihedral> GroTop::dihedrals(const QString &atm0, const QString &atm1,
                                         const QString &atm2, const QString &atm3,
                                         int func_type) const
{
    return dih_potentials.values( searchForDihType(atm0,atm1,atm2,atm3,func_type) );
}

/** Return the atom types loaded from this file */
QHash<QString,GromacsAtomType> GroTop::atomTypes() const
{
    return atom_types;
}

/** Return the bond potentials loaded from this file */
QMultiHash<QString,SireMM::GromacsBond> GroTop::bondPotentials() const
{
    return bond_potentials;
}

/** Return the angle potentials loaded from this file */
QMultiHash<QString,SireMM::GromacsAngle> GroTop::anglePotentials() const
{
    return ang_potentials;
}

/** Return the dihedral potentials loaded from this file */
QMultiHash<QString,SireMM::GromacsDihedral> GroTop::dihedralPotentials() const
{
    return dih_potentials;
}

/** Return the moleculetype with name 'name'. This returns an invalid (empty)
    GroMolType if one with this name does not exist */
GroMolType GroTop::moleculeType(const QString &name) const
{
    for (const auto moltype : moltypes)
    {
        if (moltype.name() == name)
            return moltype;
    }

    return GroMolType();
}

/** Return all of the moleculetypes that have been loaded from this file */
QVector<GroMolType> GroTop::moleculeTypes() const
{
    return moltypes;
}

/** Return whether or not the gromacs preprocessor would change these lines */
static bool gromacs_preprocess_would_change(const QVector<QString> &lines,
                                            bool use_parallel,
                                            const QHash<QString,QString> &defines)
{
    //create the regexps that are needed to find all of the
    //data that may be #define'd
    QVector<QRegularExpression> regexps;

    if (not defines.isEmpty())
    {
        regexps.reserve(defines.count());

        for (const auto key : defines.keys())
        {
            regexps.append( QRegularExpression( QString("\\s+%1\\s*").arg(key) ) );
        }
    }

    //function that says whether or not an individual line would change
    auto lineWillChange = [&](const QString &line)
    {
        if (line.indexOf(QLatin1String(";")) != -1 or
            line.indexOf(QLatin1String("#include")) != -1 or
            line.indexOf(QLatin1String("#ifdef")) != -1 or
            line.indexOf(QLatin1String("#ifndef")) != -1 or
            line.indexOf(QLatin1String("#else")) != -1 or
            line.indexOf(QLatin1String("#endif")) != -1 or
            line.indexOf(QLatin1String("#define")) != -1 or
            line.indexOf(QLatin1String("#error")) != -1 )
        {
            return true;
        }
        else
        {
            for (int i=0; i<regexps.count(); ++i)
            {
                if (line.contains( regexps.constData()[i] ))
                    return true;
            }

            if (line.trimmed().endsWith("\\"))
            {
                //this is a continuation line
                return true;
            }

            return false;
        }
    };

    const auto lines_data = lines.constData();

    if (use_parallel)
    {
        QMutex mutex;

        bool must_change = false;

        tbb::parallel_for( tbb::blocked_range<int>(0, lines.count()),
                           [&](const tbb::blocked_range<int> &r)
        {
            if (not must_change)
            {
                for (int i=r.begin(); i<r.end(); ++i)
                {
                    if (lineWillChange(lines_data[i]))
                    {
                        QMutexLocker lkr(&mutex);
                        must_change = true;
                        break;
                    }
                }
            }
        });

        return must_change;
    }
    else
    {
        for (int i=0; i<lines.count(); ++i)
        {
            if (lineWillChange(lines_data[i]))
                return true;
        }
    }

   return false;
}

/** Return the full path to the file 'filename' searching through the
    Gromacs file path. This throws an exception if the file is not found */
QString GroTop::findIncludeFile(QString filename, QString current_dir)
{
    //new file, so first see if this filename is absolute
    QFileInfo file(filename);

    //is the filename absolute?
    if (file.isAbsolute())
    {
        if (not (file.exists() and file.isReadable()))
        {
            throw SireError::io_error( QObject::tr(
                "Cannot find the file '%1'. Please make sure that this file exists "
                "and is readable").arg(filename), CODELOC );
        }

        return filename;
    }

    //does this exist from the current directory?
    file = QFileInfo( QString("%1/%2").arg(current_dir).arg(filename) );

    if (file.exists() and file.isReadable())
        return file.absoluteFilePath();

    //otherwise search the GROMACS_PATH
    for (const auto path : include_path)
    {
        file = QFileInfo( QString("%1/%2").arg(path).arg(filename) );

        if (file.exists() and file.isReadable())
        {
            return file.absoluteFilePath();
        }
    }

    //nothing was found!
    throw SireError::io_error( QObject::tr(
            "Cannot find the file '%1' using GROMACS_PATH = [ %2 ], current directory '%3'. "
            "Please make "
            "sure the file exists and is readable within your GROMACS_PATH from the "
            "current directory '%3' (e.g. "
            "set the GROMACS_PATH environment variable to include the directory "
            "that contains '%1', or copy this file into one of the existing "
            "directories [ %2 ])")
                .arg(filename).arg(include_path.join(", ")).arg(current_dir), CODELOC );

    return QString();
}

/** This function will use the Gromacs search path to find and load the
    passed include file. This will load the file and return the
    un-preprocessed text. The file, together with its QFileInfo, will
    be saved in the 'included_files' hash */
QVector<QString> GroTop::loadInclude(QString filename, QString current_dir)
{
    //try to find the file
    QString absfile = findIncludeFile(filename, current_dir);

    //now load the file
    return MoleculeParser::readTextFile(absfile);
}

/** This function scans through a set of gromacs file lines and expands all
    macros, removes all comments and includes all #included files */
QVector<QString> GroTop::preprocess(const QVector<QString> &lines,
                                    QHash<QString,QString> &defines,
                                    const QString &current_directory,
                                    const QString &parent_file)
{
    //first, scan through to see if anything needs changing
    if (not gromacs_preprocess_would_change(lines, usesParallel(), defines))
    {
        //nothing to do
        return lines;
    }

    //Ok, we have to change the lines...
    QVector<QString> new_lines;
    new_lines.reserve(lines.count());

    //regexps used to parse the files...
    QRegularExpression include_regexp("\\#include\\s*(<([^\"<>|\\b]+)>|\"([^\"<>|\\b]+)\")");

    //loop through all of the lines...
    QVectorIterator<QString> lines_it(lines);

    QList<bool> ifparse;

    while (lines_it.hasNext())
    {
        QString line = lines_it.next();

        //remove any comments
        if (line.indexOf(QLatin1String(";")) != -1)
        {
            line = line.mid(0, line.indexOf(QLatin1String(";"))).simplified();

            //this is just an empty line, so ignore it
            if (line.isEmpty())
            {
                continue;
            }
        }
        else if (line.startsWith("*"))
        {
            //the whole line is a comment
            continue;
        }
        else
        {
            //simplify the line to remove weirdness
            line = line.simplified();
        }

        //now look to see if the line should be joined to the next line
        while (line.endsWith("\\"))
        {
            if (not lines_it.hasNext())
            {
                throw SireIO::parse_error( QObject::tr(
                    "Continuation line on the last line of the Gromacs file! '%1'")
                        .arg(line), CODELOC );
            }

            line += lines_it.next();
            line = line.simplified();
        }

        //first, look to see if the line starts with #error, as this should
        //terminate processing
        if (line.startsWith("#error"))
        {
            //stop processing, and pass the error to the user
            line = line.mid(6).simplified();
            throw SireIO::parse_error( QObject::tr(
                "Error in Gromacs file! '%1'")
                    .arg(line), CODELOC );
        }

        //now look to see if there is an #ifdef
        if (line.startsWith("#ifdef"))
        {
            //we have an ifdef - has it been defined?
            auto symbol = line.split(" ", QString::SkipEmptyParts).last();

            //push the current parse state (whether we parse if or else)
            ifparse.append( defines.value(symbol,"0") != "0" );
            continue;
        }

        //now look to see if there is an #ifndef
        if (line.startsWith("#ifndef"))
        {
            //we have an ifndef - has it been defined?
            auto symbol = line.split(" ", QString::SkipEmptyParts).last();

            //push the current parse state (whether we parse if or else)
            ifparse.append( defines.value(symbol,"0") == "0" );
            continue;
        }

        if (line == "#else")
        {
            //switch the last ifdef state
            if (ifparse.isEmpty())
                throw SireIO::parse_error( QObject::tr(
                    "Unmatched '#else' in the GROMACS file!"), CODELOC );

            ifparse.last() = not ifparse.last();
            continue;
        }

        if (line == "#endif")
        {
            //pop off the last 'ifdef' state
            if (ifparse.isEmpty())
                throw SireIO::parse_error( QObject::tr(
                    "Unmatched '#endif' in the GROMACS file!"), CODELOC );

            ifparse.removeLast();
            continue;
        }

        if (not ifparse.isEmpty())
        {
            //are we allowed to read this?
            if (not ifparse.last())
            {
                //no, this is blocked out
                continue;
            }
        }

        //now look for any #define lines
        if (line.startsWith("#define"))
        {
            auto words = line.split(" ", QString::SkipEmptyParts);

            if (words.count() == 1)
                throw SireIO::parse_error( QObject::tr(
                    "Malformed #define line in Gromacs file? %1").arg(line), CODELOC );

            if (words.count() == 2)
            {
                defines.insert(words[1], "1");
            }
            else
            {
                auto key = words[1];
                words.takeFirst();
                words.takeFirst();
                defines.insert(key, words.join(" "));
            }

            continue;
        }

        //now try to substitute any 'defines' in the line with their defined values
        for (auto it = defines.constBegin(); it != defines.constEnd(); ++it)
        {
            if (line.indexOf(it.key()) != -1)
            {
                auto words = line.split(" ", QString::SkipEmptyParts);

                for (int i=0; i<words.count(); ++i)
                {
                    if (words[i] == it.key())
                    {
                        words[i] = it.value();
                    }
                }

                line = words.join(" ");
            }
        }

        //now look for #include lines
        if (line.startsWith("#include"))
        {
            //now insert the contents of any included files
            auto m = include_regexp.match(line);

            if (not m.hasMatch())
            {
                throw SireIO::parse_error( QObject::tr(
                    "Malformed #include line in Gromacs file? %1").arg(line), CODELOC );
            }

            //we have to include a file
            auto filename = m.captured( m.lastCapturedIndex() );

            //now find the absolute path to the file...
            auto absfile = findIncludeFile(filename, current_directory);

            //now load the file
            auto included_lines = MoleculeParser::readTextFile(absfile);

            //now get the absolute path to the included file
            auto parts = absfile.split("/");
            parts.removeLast();

            //fully preprocess these lines using the current set of defines
            included_lines = preprocess(included_lines, defines,
                                        parts.join("/"), absfile);

            //add these included lines to the set
            new_lines.reserve( new_lines.count() + included_lines.count() );
            new_lines += included_lines;

            //finally, record that this file depends on the included file
            included_files[parent_file].append(absfile);

            continue;
        }

        //finally, make sure that we have not missed any '#' directives...
        if (line.startsWith("#"))
        {
            throw SireIO::parse_error( QObject::tr(
                "Unrecognised directive on Gromacs file line '%1'").arg(line), CODELOC );
        }

        //skip empty lines
        if (not line.isEmpty())
        {
            //otherwise this is a normal line, so append this to the set of new_lines
            new_lines.append(line);
        }
    }

    if (not ifparse.isEmpty())
    {
        throw SireIO::parse_error( QObject::tr(
            "Unmatched #ifdef or #ifndef in Gromacs file!"), CODELOC );
    }

    return new_lines;
}

/** Return the non-bonded function type for the molecules in this file */
int GroTop::nonBondedFunctionType() const
{
    return nb_func_type;
}

/** Return the combining rules to use for the molecules in this file */
int GroTop::combiningRules() const
{
    return combining_rule;
}

/** Return the Lennard Jones fudge factor for the molecules in this file */
double GroTop::fudgeLJ() const
{
    return fudge_lj;
}

/** Return the electrostatic fudge factor for the molecules in this file */
double GroTop::fudgeQQ() const
{
    return fudge_qq;
}

/** Return whether or not the non-bonded pairs should be automatically generated
    for the molecules in this file */
bool GroTop::generateNonBondedPairs() const
{
    return generate_pairs;
}

/** Return the expanded set of lines (after preprocessing) */
const QVector<QString>& GroTop::expandedLines() const
{
    return expanded_lines;
}

/** Public function used to return the list of post-processed lines */
QStringList GroTop::postprocessedLines() const
{
    return expanded_lines.toList();
}

/** Internal function, called by ::interpret() that processes all of the data
    from all of the directives, returning a set of warnings */
QStringList GroTop::processDirectives(const QMap<int,QString> &taglocs,
                                      const QHash<QString,int> &ntags)
{
    //internal function that returns the lines associated with the
    //specified directive
    auto getLines = [&](const QString &directive, int n) -> QStringList
    {
        if (n >= ntags.value(directive,0))
        {
            return QStringList();
        }

        bool found = false;
        int start = 0;
        int end = expandedLines().count();

        //find the tag
        for (auto it = taglocs.constBegin(); it != taglocs.constEnd(); ++it)
        {
            if (it.value() == directive)
            {
                if (n == 0)
                {
                    found = true;
                    start = it.key()+1;

                    ++it;

                    if (it != taglocs.constEnd())
                    {
                        end = it.key();
                    }

                    break;
                }
                else
                    n -= 1;
            }
        }

        if (not found)
            throw SireError::program_bug( QObject::tr(
                "Cannot find tag '%1' at index '%2'. This should not happen!")
                    .arg(directive).arg(n), CODELOC );

        QStringList lines;

        for (int i=start; i<end; ++i)
        {
            lines.append( expandedLines().constData()[i] );
        }

        return lines;
    };

    //return the lines associated with the directive at line 'linenum'
    auto getDirectiveLines = [&](int linenum) -> QStringList
    {
        auto it = taglocs.constFind(linenum);

        if (it == taglocs.constEnd())
            throw SireError::program_bug( QObject::tr(
                "Cannot find a tag associated with line '%1'. This should not happen!")
                    .arg(linenum), CODELOC );

        int start = it.key() + 1;
        int end = expandedLines().count();

        ++it;

        if (it != taglocs.constEnd())
            end = it.key();

        QStringList lines;

        for (int i=start; i<end; ++i)
        {
            lines.append( expandedLines().constData()[i] );
        }

        return lines;
    };

    //return all of the lines associated with all copies of the passed directive
    auto getAllLines = [&](const QString &directive) -> QStringList
    {
        QStringList lines;

        for (int i=0; i<ntags.value(directive,0); ++i)
        {
            lines += getLines(directive, i);
        }

        return lines;
    };

    //interpret a bool from the passed string
    auto gromacs_toBool = [&](const QString &word, bool *ok)
    {
        QString w = word.toLower();

        if (ok) *ok = true;

        if (w == "yes" or w == "y" or w == "true" or w == "1")
        {
            return true;
        }
        else if (w == "no" or w == "n" or w == "false" or w == "0")
        {
            return false;
        }
        else
        {
            if (ok) *ok = false;
            return false;
        }
    };

    //internal function to process the defaults lines
    auto processDefaults = [&]()
    {
        QStringList warnings;

        //there should only be one defaults line
        const auto lines = getLines("defaults", 0);

        if (lines.isEmpty())
            throw SireIO::parse_error( QObject::tr(
                "The required data for the '[defaults]' directive in Gromacs is "
                "not supplied. This is not a valid Gromacs topology file!"), CODELOC );

        auto words = lines[0].split(" ", QString::SkipEmptyParts);

        //there should be five words; non-bonded function type, combinination rule,
        //                            generate pairs, fudge LJ and fudge QQ
        if (words.count() < 5)
        {
            throw SireIO::parse_error( QObject::tr(
                "There is insufficient data for the '[defaults]' line '%1'. This is "
                "not a valid Gromacs topology file!").arg(lines[0]), CODELOC );
        }

        bool ok;
        int nbtyp = words[0].toInt(&ok);

        if (not ok)
            throw SireIO::parse_error( QObject::tr(
                "The first value for the '[defaults]' line '%1' is not an integer. "
                "This is not a valid Gromacs topology file!").arg(lines[0]), CODELOC );

        int combrule = words[1].toInt(&ok);

        if (not ok)
            throw SireIO::parse_error( QObject::tr(
                "The second value for the '[defaults]' line '%1' is not an integer. "
                "This is not a valid Gromacs topology file!").arg(lines[0]), CODELOC );

        bool gen_pairs = gromacs_toBool(words[2], &ok);

        if (not ok)
            throw SireIO::parse_error( QObject::tr(
                "The third value for the '[defaults]' line '%1' is not a yes/no. "
                "This is not a valid Gromacs topology file!").arg(lines[0]), CODELOC );

        double lj = words[3].toDouble(&ok);

        if (not ok)
            throw SireIO::parse_error( QObject::tr(
                "The fourth value for the '[defaults]' line '%1' is not a double. "
                "This is not a valid Gromacs topology file!").arg(lines[0]), CODELOC );

        double qq = words[4].toDouble(&ok);

        if (not ok)
            throw SireIO::parse_error( QObject::tr(
                "The fifth value for the '[defaults]' line '%1' is not a double. "
                "This is not a valid Gromacs topology file!").arg(lines[0]), CODELOC );

        //validate and then save these values
        if (nbtyp <= 0 or nbtyp > 2)
        {
            warnings.append( QObject::tr("A non-supported non-bonded function type (%1) "
              "is requested.").arg(nbtyp) );
        }

        if (combrule <= 0 or combrule > 3)
        {
            warnings.append( QObject::tr("A non-supported combinig rule (%1) is requested!")
                .arg(combrule) );
        }

        if (lj < 0 or lj > 1)
        {
            warnings.append( QObject::tr("An invalid value of fudge_lj (%1) is requested!")
                .arg(lj) );

            if (lj < 0) lj = 0;
            else if (lj > 1) lj = 1;
        }

        if (qq < 0 or qq > 1)
        {
            warnings.append( QObject::tr("An invalid value of fudge_qq (%1) is requested!")
                .arg(qq) );

            if (qq < 0) qq = 0;
            else if (qq > 1) qq = 1;
        }

        // Gromacs uses a non-exact value of the Amber fudge_qq (1.0/1.2). Correct this
        // in case we convert to another file format
        if (std::abs( qq - (1.0/1.2) ) < 0.01)
        {
            qq = 1.0 / 1.2;
        }

        nb_func_type = nbtyp;
        combining_rule = combrule;
        fudge_lj = lj;
        fudge_qq = qq;
        generate_pairs = gen_pairs;

        return warnings;
    };

    //wildcard atomtype (this is 'X' in gromacs files)
    static const QString wildcard_atomtype = "X";

    //internal function to process the atomtypes lines
    auto processAtomTypes = [&]()
    {
        QStringList warnings;

        //get all 'atomtypes' lines
        const auto lines = getAllLines("atomtypes");

        //the database of all atom types
        QHash<QString,GromacsAtomType> typs;

        //now parse each atom
        for (const auto line : lines)
        {
            const auto words = line.split(" ", QString::SkipEmptyParts);

            //should either have 2 words (atom type, mass) or
            //have 6 words; atom type, mass, charge, type, V, W or
            //have 7 words; atom type, atom number, mass, charge, type, V, W
            //have 8 words; atom type, bond type, atom number, mass, charge, type, V, W
            if (words.count() < 2)
            {
                warnings.append( QObject::tr( "There is not enough data for the "
                  "atomtype data '%1'. Skipping this line." ).arg(line) );
                continue;
            }

            GromacsAtomType typ;

            if (words.count() < 6)
            {
                //only getting the atom type and mass
                bool ok_mass;
                double mass = words[1].toDouble( &ok_mass );

                if (not ok_mass)
                {
                    warnings.append( QObject::tr( "Could not interpret the atom type data "
                       "from line '%1'. Skipping this line.").arg(line) );
                    continue;
                }

                typ = GromacsAtomType(words[0], mass*g_per_mol);
            }
            else if (words.count() < 7)
            {
                bool ok_mass, ok_charge, ok_ptyp, ok_v, ok_w;
                double mass = words[1].toDouble(&ok_mass);
                double chg = words[2].toDouble(&ok_charge);
                auto ptyp = GromacsAtomType::toParticleType(words[3], &ok_ptyp);
                double v = words[4].toDouble(&ok_v);
                double w = words[5].toDouble(&ok_w);

                if (not (ok_mass and ok_charge and ok_ptyp and ok_v and ok_w))
                {
                    warnings.append( QObject::tr( "Could not interpret the atom type data "
                      "from line '%1'. Skipping this line.").arg(line) );
                    continue;
                }

                typ = GromacsAtomType(words[0], mass*g_per_mol, chg*mod_electron,
                                      ptyp, ::toLJParameter(v,w,combining_rule));
            }
            else if (words.count() < 8)
            {
                bool ok_mass, ok_elem, ok_charge, ok_ptyp, ok_v, ok_w;
                int nprotons = words[1].toInt(&ok_elem);
                double mass = words[2].toDouble(&ok_mass);
                double chg = words[3].toDouble(&ok_charge);
                auto ptyp = GromacsAtomType::toParticleType(words[4], &ok_ptyp);
                double v = words[5].toDouble(&ok_v);
                double w = words[6].toDouble(&ok_w);

                if (not (ok_mass and ok_charge and ok_ptyp and ok_v and ok_w))
                {
                    warnings.append( QObject::tr( "Could not interpret the atom type data "
                      "from line '%1'. Skipping this line.").arg(line) );
                    continue;
                }

                if (ok_elem)
                {
                    typ = GromacsAtomType(words[0], mass*g_per_mol, chg*mod_electron,
                                          ptyp, ::toLJParameter(v,w,combining_rule),
                                          Element(nprotons));
                }
                else
                {
                    //some gromacs files don't use 'nprotons', but instead give
                    //a "bond_type" ambertype
                    typ = GromacsAtomType(words[0], words[1], mass*g_per_mol, chg*mod_electron,
                                          ptyp, ::toLJParameter(v,w,combining_rule));
                }
            }
            else if (words.count() < 9)
            {
                bool ok_mass, ok_elem, ok_charge, ok_ptyp, ok_v, ok_w;
                int nprotons = words[2].toInt(&ok_elem);
                double mass = words[3].toDouble(&ok_mass);
                double chg = words[4].toDouble(&ok_charge);
                auto ptyp = GromacsAtomType::toParticleType(words[5], &ok_ptyp);
                double v = words[6].toDouble(&ok_v);
                double w = words[7].toDouble(&ok_w);

                if (not (ok_mass and ok_charge and ok_ptyp and ok_v and ok_w))
                {
                    warnings.append( QObject::tr( "Could not interpret the atom type data "
                      "from line '%1'. Skipping this line.").arg(line) );
                    continue;
                }

                typ = GromacsAtomType(words[0], words[1], mass*g_per_mol, chg*mod_electron,
                                      ptyp, ::toLJParameter(v,w,combining_rule),
                                      Element(nprotons));
            }
            else
            {
                warnings.append( QObject::tr( "The atomtype line '%1' contains more data "
                  "than is expected!").arg(line) );
            }

            if (typs.contains(typ.atomType()))
            {
                //only replace if the new type is fully specified
                if (typ.hasMassOnly())
                    continue;

                warnings.append( QObject::tr( "The data for atom type '%1' exists already! "
                 "This will now be replaced with new data from line '%2'")
                    .arg(typ.atomType()).arg(line) );
            }

            typs.insert( typ.atomType(), typ );
        }

        //save the database of types
        atom_types = typs;

        return warnings;
    };

    //internal function to process the bondtypes lines
    auto processBondTypes = [&]()
    {
        QStringList warnings;

        //get all 'bondtypes' lines
        const auto lines = getAllLines("bondtypes");

        //save into a database of bonds
        QMultiHash<QString,GromacsBond> bnds;

        for (const auto line : lines)
        {
            //each line should contain the atom types of the two atoms, then
            //the function type, then the parameters for the function
            const auto words = line.split(" ", QString::SkipEmptyParts);

            if (words.count() < 3)
            {
                warnings.append( QObject::tr("There is not enough data on the "
                  "line '%1' to extract a Gromacs bond parameter. Skipping line.")
                    .arg(line) );
                continue;
            }

            const auto atm0 = words[0];
            const auto atm1 = words[1];

            bool ok;
            int func_type = words[2].toInt(&ok);

            if (not ok)
            {
                warnings.append( QObject::tr("Unable to determine the function type "
                  "for the bond on line '%1'. Skipping line.")
                    .arg(line) );
                continue;
            }

            //now read in all of the remaining values as numbers...
            QList<double> params;

            for (int i=3; i<words.count(); ++i)
            {
                double param = words[i].toDouble(&ok);

                if (ok) params.append(param);
            }

            GromacsBond bond;

            try
            {
                bond = GromacsBond(func_type, params);
            }
            catch(const SireError::exception &e)
            {
                warnings.append( QObject::tr("Unable to extract the correct information "
                  "to form a bond from line '%1'. Error is '%2'")
                    .arg(line).arg(e.error()) );
                continue;
            }

            QString key = get_bond_id(atm0,atm1,bond.functionType());
            bnds.insertMulti(key, bond);
        }

        bond_potentials = bnds;

        return warnings;
    };

    //internal function to process the pairtypes lines
    auto processPairTypes = [&]()
    {
        QStringList warnings;

        //get all 'bondtypes' lines
        const auto lines = getAllLines("pairtypes");

        for (const auto line : lines)
        {
            warnings.append( QString("Ignoring 'pairtypes' %1").arg(line) );
        }

        return warnings;
    };

    //internal function to process the angletypes lines
    auto processAngleTypes = [&]()
    {
        QStringList warnings;

        //get all 'bondtypes' lines
        const auto lines = getAllLines("angletypes");

        //save into a database of angles
        QMultiHash<QString,GromacsAngle> angs;

        for (const auto line : lines)
        {
            //each line should contain the atom types of the three atoms, then
            //the function type, then the parameters for the function
            const auto words = line.split(" ", QString::SkipEmptyParts);

            if (words.count() < 4)
            {
                warnings.append( QObject::tr("There is not enough data on the "
                  "line '%1' to extract a Gromacs angle parameter. Skipping line.")
                    .arg(line) );
                continue;
            }

            const auto atm0 = words[0];
            const auto atm1 = words[1];
            const auto atm2 = words[2];

            bool ok;
            int func_type = words[3].toInt(&ok);

            if (not ok)
            {
                warnings.append( QObject::tr("Unable to determine the function type "
                  "for the angle on line '%1'. Skipping line.")
                    .arg(line) );
                continue;
            }

            //now read in all of the remaining values as numbers...
            QList<double> params;

            for (int i=4; i<words.count(); ++i)
            {
                double param = words[i].toDouble(&ok);

                if (ok) params.append(param);
            }

            GromacsAngle angle;

            try
            {
                angle = GromacsAngle(func_type, params);
            }
            catch(const SireError::exception &e)
            {
                warnings.append( QObject::tr("Unable to extract the correct information "
                  "to form an angle from line '%1'. Error is '%2'")
                    .arg(line).arg(e.error()) );
                continue;
            }

            QString key = get_angle_id(atm0,atm1,atm2,angle.functionType());
            angs.insertMulti(key, angle);
        }

        ang_potentials = angs;

        return warnings;
    };

    //internal function to process the dihedraltypes lines
    auto processDihedralTypes = [&]()
    {
        QStringList warnings;

        //get all 'bondtypes' lines
        const auto lines = getAllLines("dihedraltypes");

        //save into a database of dihedrals
        QMultiHash<QString,GromacsDihedral> dihs;

        for (const auto line : lines)
        {
            //each line should contain the atom types of the four atoms, then
            //the function type, then the parameters for the function.
            //(however, some files have the atom types of just two atoms, which
            // I assume are the two middle atoms of the dihedral...)
            const auto words = line.split(" ", QString::SkipEmptyParts);

            if (words.count() < 3)
            {
                warnings.append( QObject::tr("There is not enough data on the "
                  "line '%1' to extract a Gromacs dihedral parameter. Skipping line.")
                    .arg(line) );
                continue;
            }

            //first, let's try to parse this assuming that it is a 2-atom dihedral line...
            //(most cases this should fail)
            auto atm0 = wildcard_atomtype;
            auto atm1 = words[0];
            auto atm2 = words[1];
            auto atm3 = wildcard_atomtype;

            GromacsDihedral dihedral;

            bool ok;
            int func_type = words[2].toInt(&ok);

            if (ok)
            {
                //this may be a two-atom dihedral - read in the rest of the parameters
                QList<double> params;

                for (int i=3; i<words.count(); ++i)
                {
                    double param = words[i].toDouble(&ok);
                    if (ok) params.append(param);
                }

                try
                {
                    dihedral = GromacsDihedral(func_type, params);
                    ok = true;
                }
                catch(...)
                {
                    ok = false;
                }
            }

            if (not ok)
            {
                //we couldn't parse as a two-atom dihedral, so parse as a four-atom dihedral

                if (words.count() < 5)
                {
                    warnings.append( QObject::tr("There is not enough data on the "
                      "line '%1' to extract a Gromacs dihedral parameter. Skipping line.")
                        .arg(line) );
                    continue;
                }

                atm0 = words[0];
                atm1 = words[1];
                atm2 = words[2];
                atm3 = words[3];

                bool ok;
                int func_type = words[4].toInt(&ok);

                if (not ok)
                {
                    warnings.append( QObject::tr("Unable to determine the function type "
                      "for the dihedral on line '%1'. Skipping line.")
                        .arg(line) );
                    continue;
                }

                //now read in all of the remaining values as numbers...
                QList<double> params;

                for (int i=5; i<words.count(); ++i)
                {
                    double param = words[i].toDouble(&ok);

                    if (ok) params.append(param);
                }

                try
                {
                    dihedral = GromacsDihedral(func_type, params);
                }
                catch(const SireError::exception &e)
                {
                    warnings.append( QObject::tr("Unable to extract the correct information "
                      "to form a dihedral from line '%1'. Error is '%2'")
                        .arg(line).arg(e.error()) );
                    continue;
                }
            }

            QString key = get_dihedral_id(atm0,atm1,atm2,atm3,dihedral.functionType());
            dihs.insertMulti(key, dihedral);
        }

        dih_potentials = dihs;

        return warnings;
    };

    //internal function to process the constrainttypes lines
    auto processConstraintTypes = [&]()
    {
        QStringList warnings;

        //get all 'bondtypes' lines
        const auto lines = getAllLines("constrainttypes");

        for (const auto line : lines)
        {
            warnings.append( QString("Ignoring 'constrainttypes' %1").arg(line) );
        }

        return warnings;
    };

    //internal function to process the nonbond_params lines
    auto processNonBondParams = [&]()
    {
        QStringList warnings;

        //get all 'bondtypes' lines
        const auto lines = getAllLines("nonbond_params");

        for (const auto line : lines)
        {
            warnings.append( QString("Ignoring 'nonbond_params' %1").arg(line) );
        }

        return warnings;
    };

    //internal function to process moleculetype lines
    auto processMoleculeTypes = [&]()
    {
        QStringList warnings;

        //how many moleculetypes are there? Divide them up and get
        //the child tags for each moleculetype
        QList< QMultiHash<QString,int> > moltags;
        {
            //list of tags that are valid within a moleculetype
            const QStringList valid_tags = { "atoms", "bonds", "pairs", "pairs_nb",
                                             "angles", "dihedrals", "exclusions",
                                             "contraints", "settles", "virtual_sites2",
                                             "virtual_sitesn", "position_restraints",
                                             "distance_restraints", "orientation_restraints",
                                             "angle_restraints", "angle_restraints_z" };

            auto it = taglocs.constBegin();

            while (it != taglocs.constEnd())
            {
                if (it.value() == "moleculetype")
                {
                    //we have found another molecule - save the location
                    //of all of its child tags
                    QMultiHash<QString,int> tags;
                    tags.insert(it.value(), it.key());
                    ++it;

                    while (it != taglocs.constEnd())
                    {
                        //save all child tags until we reach the end
                        //of definition of this moleculetype
                        if (valid_tags.contains(it.value()))
                        {
                            //this is a valid child tag - save its location
                            //(note that a tag can exist multiple times!)
                            tags.insertMulti(it.value(), it.key());
                            ++it;
                        }
                        else if (it.value() == "moleculetype")
                        {
                            //this is the next molecule
                            break;
                        }
                        else
                        {
                            //this is the end of the 'moleculetype'
                            ++it;
                            break;
                        }
                    }

                    moltags.append(tags);
                }
                else
                {
                    ++it;
                }
            }
        }

        //now we define a set of functions that are needed to parse the
        //various child tags

        //function that extract the metadata about the moleculetype
        //and returns it as a 'GroMolType' object
        auto getMolType = [&]( int linenum )
        {
            GroMolType moltype;

            //get the directives for this molecule - there should be
            //one line that contains the name and number of excluded atoms
            const auto lines = getDirectiveLines(linenum);

            if (lines.count() != 1)
            {
                moltype.addWarning( QObject::tr("Expecting only one line that "
                  "provides the name and number of excluded atoms for this moleculetype. "
                  "Instead, the number of lines is %1 : [\n%2\n]")
                    .arg(lines.count()).arg(lines.join("\n")) );
            }

            if (lines.count() > 0)
            {
                //try to read the infromation from the first line only
                const auto words = lines[0].split(" ");

                if (words.count() != 2)
                {
                    moltype.addWarning( QObject::tr("Expecting two words for the "
                      "moleculetype line, containing the name and number of excluded "
                      "atoms. Instead we get '%1'").arg(lines[0]) );
                }

                if (words.count() > 0)
                {
                    moltype.setName( words[0] );
                }

                if (words.count() > 1)
                {
                    bool ok;
                    qint64 nexcl = words[1].toInt(&ok);

                    if (not ok)
                    {
                        moltype.addWarning( QObject::tr("Expecting the second word in "
                          "the moleculetype line '%1' to be the number of excluded "
                          "atoms. It isn't!").arg(lines[0]) );
                    }
                    else
                    {
                        moltype.setNExcludedAtoms(nexcl);
                    }
                }
            }

            return moltype;
        };

        //function that extracts all of the information from the 'atoms' lines
        //and adds it to the passed GroMolType
        auto addAtomsTo = [&](GroMolType &moltype, int linenum)
        {
            QStringList lines = getDirectiveLines(linenum);

            for (const auto line : lines)
            {
                //each line should contain index number, atom type,
                //residue number, residue name, atom name, charge group number,
                //charge (mod_electron) and mass (atomic mass)
                const auto words = line.split(" ");

                if (words.count() < 6)
                {
                    moltype.addWarning( QObject::tr("Cannot extract atom information "
                      "from the line '%1' as it should contain at least six words "
                      "(pieces of information)").arg(line) );

                    continue;
                }
                else if (words.count() > 8)
                {
                    moltype.addWarning( QObject::tr("The line containing atom information "
                       "'%1' contains more information than can be parsed. It should only "
                       "contain six-eight words (pieces of information)").arg(line) );
                }

                bool ok_idx, ok_resnum, ok_chggrp, ok_chg, ok_mass;

                const qint64 atomnum = words[0].toInt(&ok_idx);
                const auto atomtyp = words[1];
                const auto resnum = words[2].toInt(&ok_resnum);
                const auto resnam = words[3];
                const auto atmnam = words[4];
                const qint64 chggrp = words[5].toInt(&ok_chggrp);

                double chg = 0; ok_chg = true;
                if (words.count() > 6)
                    chg = words[6].toDouble(&ok_chg);

                double mass = 0; ok_mass = true; bool found_mass = false;
                if (words.count() > 7)
                {
                    mass = words[7].toDouble(&ok_mass);
                    found_mass = true;
                }

                if (not (ok_idx and ok_resnum and ok_chggrp and ok_chg and ok_mass))
                {
                    moltype.addWarning( QObject::tr("Could not interpret the necessary "
                      "atom information from the line '%1' | %2 %3 %4 %5 %6")
                        .arg(line).arg(ok_idx).arg(ok_resnum).arg(ok_chggrp)
                        .arg(ok_chg).arg(ok_mass) );
                    continue;
                }

                GroAtom atom;
                atom.setNumber(atomnum);
                atom.setAtomType(atomtyp);
                atom.setResidueNumber(resnum);
                atom.setResidueName(resnam);
                atom.setName(atmnam);
                atom.setChargeGroup(chggrp);
                atom.setCharge(chg*mod_electron);
                atom.setMass(mass*g_per_mol);

                //we now need to look up the atom type of this atom to see if there
                //is a separate bond_type
                const auto atom_type = atom_types.value(atomtyp);

                if ((not atom_type.isNull()) and atom_type.bondType() != atomtyp)
                {
                    atom.setBondType(atom_type.bondType());
                }

                //now do the same to assign the mass if it has not been given explicitly
                if ( (not found_mass) and (not atom_type.isNull()) )
                {
                    atom.setMass( atom_type.mass() );
                }

                moltype.addAtom(atom);
            }
        };

        //function that extracts all of the information from the 'bonds' lines
        auto addBondsTo = [&](GroMolType &moltype, int linenum)
        {
            QStringList lines = getDirectiveLines(linenum);

            QMultiHash<BondID,GromacsBond> bonds;
            bonds.reserve(lines.count());

            for (const auto line : lines)
            {
                const auto words = line.split(" ");

                if (words.count() < 2)
                {
                    moltype.addWarning( QObject::tr("Cannot extract bond information "
                      "from the line '%1' as it should contain at least two words "
                      "(pieces of information)").arg(line) );
                    continue;
                }

                bool ok0, ok1;

                int atm0 = words[0].toInt(&ok0);
                int atm1 = words[1].toInt(&ok1);

                if (not (ok0 and ok1))
                {
                    moltype.addWarning( QObject::tr("Cannot extract bond information "
                      "from the line '%1' as the first two words need to be integers. ")
                            .arg(line) );
                    continue;
                }

                //now see if any information about the bond is provided...
                GromacsBond bond;

                if (words.count() > 2)
                {
                    bool ok;
                    int func_type = words[2].toInt(&ok);

                    if (not ok)
                    {
                        moltype.addWarning( QObject::tr("Unable to extract the correct "
                           "information to form a bond from line '%1' as the third word "
                           "is not an integer.").arg(line) );
                        continue;
                    }

                    if (words.count() > 3)
                    {
                        //now read in all of the remaining values as numbers...
                        QList<double> params;

                        for (int i=3; i<words.count(); ++i)
                        {
                            double param = words[i].toDouble(&ok);

                            if (ok) params.append(param);
                        }

                        try
                        {
                            bond = GromacsBond(func_type, params);
                        }
                        catch(const SireError::exception &e)
                        {
                            moltype.addWarning( QObject::tr("Unable to extract the correct "
                              "information to form a bond from line '%1'. Error is '%2'")
                                .arg(line).arg(e.error()) );
                            continue;
                        }
                    }
                    else
                    {
                        bond = GromacsBond(func_type);
                    }
                }

                bonds.insertMulti( BondID(AtomNum(atm0),AtomNum(atm1)),
                                   bond );
            }

            //save the bonds in the molecule
            moltype.addBonds(bonds);
        };

        //function that extracts all of the information from the 'angles' lines
        auto addAnglesTo = [&](GroMolType &moltype, int linenum)
        {
            QStringList lines = getDirectiveLines(linenum);

            QMultiHash<AngleID,GromacsAngle> angs;
            angs.reserve(lines.count());

            for (const auto line : lines)
            {
                const auto words = line.split(" ");

                if (words.count() < 3)
                {
                    moltype.addWarning( QObject::tr("Cannot extract angle information "
                      "from the line '%1' as it should contain at least three words "
                      "(pieces of information)").arg(line) );
                    continue;
                }

                bool ok0, ok1, ok2;

                int atm0 = words[0].toInt(&ok0);
                int atm1 = words[1].toInt(&ok1);
                int atm2 = words[2].toInt(&ok2);

                if (not (ok0 and ok1 and ok2))
                {
                    moltype.addWarning( QObject::tr("Cannot extract angle information "
                      "from the line '%1' as the first three words need to be integers. ")
                            .arg(line) );
                    continue;
                }

                //now see if any information about the angle is provided...
                GromacsAngle angle;

                if (words.count() > 3)
                {
                    bool ok;
                    int func_type = words[3].toInt(&ok);

                    if (not ok)
                    {
                        moltype.addWarning( QObject::tr("Unable to extract the correct "
                           "information to form an angle from line '%1' as the fourth word "
                           "is not an integer.").arg(line) );
                        continue;
                    }

                    if (words.count() > 4)
                    {
                        //now read in all of the remaining values as numbers...
                        QList<double> params;

                        for (int i=4; i<words.count(); ++i)
                        {
                            double param = words[i].toDouble(&ok);

                            if (ok) params.append(param);
                        }

                        try
                        {
                            angle = GromacsAngle(func_type, params);
                        }
                        catch(const SireError::exception &e)
                        {
                            moltype.addWarning( QObject::tr("Unable to extract the correct "
                              "information to form an angle from line '%1'. Error is '%2'")
                                .arg(line).arg(e.error()) );
                            continue;
                        }
                    }
                    else
                    {
                        angle = GromacsAngle(func_type);
                    }
                }

                angs.insertMulti( AngleID(AtomNum(atm0),AtomNum(atm1),AtomNum(atm2)),
                                  angle );
            }

            //save the angles in the molecule
            moltype.addAngles(angs);
        };

        //function that extracts all of the information from the 'dihedrals' lines
        auto addDihedralsTo = [&](GroMolType &moltype, int linenum)
        {
            QStringList lines = getDirectiveLines(linenum);

            QMultiHash<DihedralID,GromacsDihedral> dihs;
            dihs.reserve(lines.count());

            for (const auto line : lines)
            {
                const auto words = line.split(" ");

                if (words.count() < 4)
                {
                    moltype.addWarning( QObject::tr("Cannot extract dihedral information "
                      "from the line '%1' as it should contain at least four words "
                      "(pieces of information)").arg(line) );
                    continue;
                }

                bool ok0, ok1, ok2, ok3;

                int atm0 = words[0].toInt(&ok0);
                int atm1 = words[1].toInt(&ok1);
                int atm2 = words[2].toInt(&ok2);
                int atm3 = words[3].toInt(&ok3);

                if (not (ok0 and ok1 and ok2 and ok3))
                {
                    moltype.addWarning( QObject::tr("Cannot extract dihedral information "
                      "from the line '%1' as the first four words need to be integers. ")
                            .arg(line) );
                    continue;
                }

                //now see if any information about the dihedral is provided...
                GromacsDihedral dihedral;

                if (words.count() > 4)
                {
                    bool ok;
                    int func_type = words[4].toInt(&ok);

                    if (not ok)
                    {
                        moltype.addWarning( QObject::tr("Unable to extract the correct "
                           "information to form a dihedral from line '%1' as the fifth word "
                           "is not an integer.").arg(line) );
                        continue;
                    }

                    if (words.count() > 5)
                    {
                        //now read in all of the remaining values as numbers...
                        QList<double> params;

                        for (int i=5; i<words.count(); ++i)
                        {
                            double param = words[i].toDouble(&ok);

                            if (ok) params.append(param);
                        }

                        try
                        {
                            dihedral = GromacsDihedral(func_type, params);
                        }
                        catch(const SireError::exception &e)
                        {
                            moltype.addWarning( QObject::tr("Unable to extract the correct "
                              "information to form a dihedral from line '%1'. Error is '%2'")
                                .arg(line).arg(e.error()) );
                            continue;
                        }
                    }
                    else
                    {
                        dihedral = GromacsDihedral(func_type);
                    }
                }

                dihs.insertMulti( DihedralID(AtomNum(atm0),AtomNum(atm1),
                                             AtomNum(atm2),AtomNum(atm3)),
                                  dihedral );
            }

            //save the dihedrals in the molecule
            moltype.addDihedrals(dihs);
        };

        //interpret the defaults so that the forcefield for each moltype can
        //be determined
        const QString elecstyle = "coulomb";
        const QString vdwstyle = _getVDWStyle(nb_func_type);
        const QString combrules = _getCombiningRules(combining_rule);

        //ok, now we know the location of all child tags of each moleculetype
        auto processMolType = [&](const QHash<QString,int> &moltag)
        {
            auto moltype = getMolType( moltag.value("moleculetype", -1) );

            for (auto linenum : moltag.values("atoms"))
            {
                addAtomsTo( moltype, linenum );
            }

            for (auto linenum : moltag.values("bonds"))
            {
                addBondsTo( moltype, linenum );
            }

            for (auto linenum : moltag.values("angles"))
            {
                addAnglesTo( moltype, linenum );
            }

            for (auto linenum : moltag.values("dihedrals"))
            {
                addDihedralsTo( moltype, linenum );
            }

            //now print out warnings for any lines that are missed...
            const QStringList missed_tags = { "pairs", "pairs_nb",
                                              "exclusions",
                                              "contraints", "settles", "virtual_sites2",
                                              "virtual_sitesn", "position_restraints",
                                              "distance_restraints", "orientation_restraints",
                                              "angle_restraints", "angle_restraints_z" };

            for (const auto tag : missed_tags)
            {
                //not parsed this tag type
                for (auto linenum : moltag.values(tag))
                {
                    moltype.addWarning( QObject::tr("\nSkipping these lines which "
                      "are part of the '%1' tag:\n%2")
                        .arg(tag).arg( getDirectiveLines(linenum).join("\n") ) );
                }
            }

            //should be finished, run some checks that this looks sane
            moltype.sanitise(elecstyle, vdwstyle, combrules, fudge_qq, fudge_lj);

            return moltype;
        };

        //set the size of the array of moltypes
        moltypes = QVector<GroMolType>(moltags.count());
        auto moltypes_array = moltypes.data();

        //load all of the molecule types (in parallel if possible)
        if (usesParallel())
        {
            tbb::parallel_for( tbb::blocked_range<int>(0,moltags.count()),
                               [&](const tbb::blocked_range<int> &r)
            {
                for (int i=r.begin(); i<r.end(); ++i)
                {
                    auto moltype = processMolType( moltags.at(i) );
                    moltypes_array[i] = moltype;
                }
            });
        }
        else
        {
            for (int i=0; i<moltags.count(); ++i)
            {
                auto moltype = processMolType( moltags.at(i) );
                moltypes_array[i] = moltype;
            }
        }

        //now collect any warnings from the types
        for (const auto moltype : moltypes)
        {
            if (not moltype.warnings().isEmpty())
            {
                warnings.append(moltype.warnings());
            }
        }

        return warnings;
    };

    //function used to parse the [system] part of the file
    auto processSystem = [&]
    {
        QStringList warnings;

        //look for the locations of the child tags of [system]
        QList< QMultiHash<QString,int> > systags;
        {
            //list of tags that are valid within a [system]
            const QStringList valid_tags = { "molecules" };

            auto it = taglocs.constBegin();

            while (it != taglocs.constEnd())
            {
                if (it.value() == "system")
                {
                    //we have found another 'system' - save the location
                    //of all of its child tags
                    QMultiHash<QString,int> tags;
                    tags.insert(it.value(), it.key());
                    ++it;

                    while (it != taglocs.constEnd())
                    {
                        //save all child tags until we reach the end
                        //of definition of this system
                        if (valid_tags.contains(it.value()))
                        {
                            //this is a valid child tag - save its location
                            //(note that a tag can exist multiple times!)
                            tags.insertMulti(it.value(), it.key());
                            ++it;
                        }
                        else
                        {
                            //this is the end of the 'system'
                            ++it;
                            break;
                        }
                    }

                    systags.append(tags);
                }
                else
                {
                    ++it;
                }
            }
        }

        //in theory, there should be one, and only one [system]
        if (systags.count() != 1)
        {
            warnings.append( QObject::tr( "There should be one, and only one "
              "[system] section in a Gromacs topology file. The number of "
              "[system] sections equals %1.").arg(systags.count()) );
            return warnings;
        }

        //now parse the two parts of [system]
        const auto tags = systags.at(0);

        if (not (tags.contains("system") and tags.contains("molecules")))
        {
            warnings.append( QObject::tr("The [system] section should contain "
              "both [system] and [molecules]. It contains '%1'")
                .arg( Sire::toString(tags) ) );
            return warnings;
        }

        //process [system] first..
        //each of these lines is part of the title of the system
        GroSystem mysys( getDirectiveLines(tags.value("system")).join(" ") );

        //now process the [molecules]
        for (auto linenum : tags.values("molecules"))
        {
            const auto lines = getDirectiveLines(linenum);

            for (const auto line : lines)
            {
                //each line should be the molecule type name, followed by the number
                const auto words = line.split(" ");

                if (words.count() < 2)
                {
                    warnings.append( QObject::tr("Cannot understand the [molecules] line "
                      "'%1' as it should have two words!").arg(line) );
                    continue;
                }

                if (words.count() > 2)
                {
                    warnings.append( QObject::tr("Ignoring the extraneous information at "
                      "the end of the [molecules] line '%1'").arg(line) );
                }

                bool ok;
                int nmols = words[1].toInt(&ok);

                if (not ok)
                {
                    warnings.append( QObject::tr("Cannot interpret the number of molecules "
                      "from the [molecules] line '%1'. The second word should be an integer "
                      "that gives the number of molecules...").arg(line) );
                    continue;
                }

                mysys.add(words[0], nmols);
            }
        }

        //save the system object to this GroTop
        grosys = mysys;

        return warnings;
    };

    //process the defaults data first, as this affects the rest of the parsing
    auto warnings = processDefaults();

    //next, read in the atom types as these have to be present before
    //reading anything else...
    warnings += processAtomTypes();

    //now we can process the other tags
    const QVector< std::function<QStringList()> > funcs =
                 { processBondTypes, processPairTypes,
                   processAngleTypes, processDihedralTypes, processConstraintTypes,
                   processNonBondParams, processMoleculeTypes, processSystem
                 };

    if (usesParallel())
    {
        QMutex mutex;

        tbb::parallel_for( tbb::blocked_range<int>(0, funcs.count()),
                           [&](const tbb::blocked_range<int> &r)
        {
            QStringList local_warnings;

            for (int i=r.begin(); i<r.end(); ++i)
            {
                local_warnings += funcs[i]();
            }

            if (not local_warnings.isEmpty())
            {
                QMutexLocker lkr(&mutex);
                warnings += local_warnings;
            }
        });
    }
    else
    {
        for (int i=0; i<funcs.count(); ++i)
        {
            warnings += funcs[i]();
        }
    }

    return warnings;
}

/** Interpret the fully expanded set of lines to extract all of the necessary data */
void GroTop::interpret()
{
    //first, go through and find the line numbers of all tags
    const QRegularExpression re("\\[\\s*([\\w\\d]+)\\s*\\]");

    //map giving the type and line number of each directive tag
    QMap<int,QString> taglocs;

    const int nlines = expandedLines().count();
    const auto lines = expandedLines().constData();

    //run through this file to find all of the directives
    if (usesParallel())
    {
        QMutex mutex;

        tbb::parallel_for( tbb::blocked_range<int>(0,nlines),
                           [&](const tbb::blocked_range<int> &r)
        {
            QMap<int,QString> mylocs;

            for (int i=r.begin(); i<r.end(); ++i)
            {
                auto m = re.match(lines[i]);

                if (m.hasMatch())
                {
                    auto tag = m.captured(1);
                    mylocs.insert(i,tag);
                }
            }

            if (not mylocs.isEmpty())
            {
                QMutexLocker lkr(&mutex);

                for (auto it=mylocs.constBegin(); it!=mylocs.constEnd(); ++it)
                {
                    taglocs.insert(it.key(), it.value());
                }
            }
        });
    }
    else
    {
        for (int i=0; i<nlines; ++i)
        {
            auto m = re.match(lines[i]);

            if (m.hasMatch())
            {
                auto tag = m.captured(1);
                taglocs.insert(i,tag);
            }
        }
    }

    //now, validate that this looks like a gromacs top file. Rules are taken
    //from page 138 of the Gromacs 5.1 PDF reference manual

    //first, count up the number of each tag
    QHash<QString,int> ntags;

    for (auto it = taglocs.constBegin(); it != taglocs.constEnd(); ++it)
    {
        if (not ntags.contains(it.value()))
        {
            ntags.insert(it.value(), 1);
        }
        else
        {
            ntags[it.value()] += 1;
        }
    }

    //there should be only one 'defaults' tag
    if (ntags.value("defaults", 0) != 1)
    {
        throw SireIO::parse_error( QObject::tr(
            "This is not a valid GROMACS topology file. Such files contain one, and one "
            "only 'defaults' directive. The number of such directives in this file is %1.")
                .arg(ntags.value("defaults",0)), CODELOC );
    }

    //now process all of the directives
    auto warnings = this->processDirectives(taglocs, ntags);

    if (not warnings.isEmpty())
    {
        parse_warnings = warnings;
    }

    this->setScore(100);
}

/** Return all of the warnings that were raised when parsing the file */
QStringList GroTop::warnings() const
{
    QStringList w = parse_warnings;

    for (const auto moltype : moltypes)
    {
        auto molwarns = moltype.warnings();

        if (not molwarns.isEmpty())
        {
            w.append( QObject::tr("Warnings for molecule type %1").arg(moltype.toString()) );
            w += molwarns;
        }
    }

    return w;
}

/** Internal function that is used to actually parse the data contained
    in the lines of the file */
void GroTop::parseLines(const QString &path, const PropertyMap &map)
{
    //first, see if there are any GROMACS defines in the passed map
    //and then preprocess the lines to create the fully expanded file to parse
    {
        QHash<QString,QString> defines;

        try
        {
            const auto p = map["GROMACS_DEFINE"];

            QStringList d;

            if (p.hasValue())
            {
                d = p.value().asA<StringProperty>().toString().split(":", QString::SkipEmptyParts);
            }
            else if (p.source() != "GROMACS_DEFINE")
            {
                d = p.source().split(":", QString::SkipEmptyParts);
            }

            for (const auto define : d)
            {
                auto words = define.split("=");

                if (words.count() == 1)
                {
                    defines.insert( words[0].simplified(), "1" );
                }
                else
                {
                    defines.insert( words[0].simplified(), words[1].simplified() );
                }
            }
        }
        catch(...)
        {}

        //now go through an expand any macros and include the contents of any
        //included files
        expanded_lines = preprocess(lines(), defines, path, ".");
    }

    //now we know that there are no macros to expand, no other files to
    //include, and everything should be ok... ;-)
    this->interpret();
}

/** This function is used to create a molecule. Any errors should be written
    to the 'errors' QStringList passed as an argument */
Molecule GroTop::createMolecule(const GroMolType &moltype, QStringList &errors,
                                const PropertyMap&) const
{
    MolStructureEditor mol;

    //first go through and create the Molecule layout
    //(the atoms are already sorted into Residues)
    int cgidx = 1;
    ResStructureEditor res;
    CGStructureEditor cgroup;

    for (const auto &atom : moltype.atoms())
    {
        if (cgroup.nAtoms() == 0)
        {
            //this is the first atom in the molecule
            cgroup = mol.add( CGName( QString::number(cgidx) ) );
            cgidx += 1;

            res = mol.add( ResNum(atom.residueNumber()) );
            res = res.rename( atom.residueName() );
        }
        else if (atom.residueNumber() != res.number() or
                 atom.residueName() != res.name())
        {
            //this atom is in a different residue
            cgroup = mol.add( CGName( QString::number(cgidx) ) );
            cgidx += 1;

            res = mol.add( ResNum(atom.residueNumber()) );
            res = res.rename(atom.residueName());
        }

        //add the atom to the residue
        auto a = res.add( AtomName(atom.name()) );
        a = a.renumber(atom.number());
        a = a.reparent(cgroup.name());
    }

    return mol.commit();
}

/** This function is used to return atom properties for the passed molecule */
GroTop::PropsAndErrors GroTop::getAtomProperties(const MoleculeInfo &molinfo,
                                                 const GroMolType &moltype) const
{
    //create space for all of the properties
    AtomStringProperty atom_type(molinfo);
    AtomStringProperty bond_type(molinfo);
    AtomIntProperty charge_group(molinfo);
    AtomCharges atom_chgs(molinfo);
    AtomMasses atom_masses(molinfo);

    AtomLJs atom_ljs(molinfo);
    AtomElements atom_elements(molinfo);
    AtomStringProperty particle_type(molinfo);

    const auto atoms = moltype.atoms();

    QStringList errors;

    bool uses_bondtypes = false;

    //loop over each atom and look up parameters
    for (int i=0; i<atoms.count(); ++i)
    {
        auto cgatomidx = molinfo.cgAtomIdx( AtomIdx(i) );

        //get the template for this atom (templates in same order as atoms)
        const auto atom = atoms.constData()[i];

        //information from the template
        atom_type.set(cgatomidx, atom.atomType());
        bond_type.set(cgatomidx, atom.bondType());

        if (atom.atomType() != atom.bondType())
            uses_bondtypes = true;

        charge_group.set(cgatomidx, atom.chargeGroup());
        atom_chgs.set(cgatomidx, atom.charge());
        atom_masses.set(cgatomidx, atom.mass());

        //information from the atom type
        const auto atmtyp = atom_types.value(atom.atomType());

        if (atmtyp.isNull())
        {
            errors.append( QObject::tr("There are no parameters for the atom "
               "type '%1', needed by atom '%2'")
                    .arg(atom.atomType()).arg(atom.toString()) );
            continue;
        }
        else if (atmtyp.hasMassOnly())
        {
            errors.append( QObject::tr("The parameters for the atom type '%1' needed "
               "by atom '%2' has mass parameters only! '%3'")
                    .arg(atom.atomType()).arg(atom.toString()).arg(atmtyp.toString()) );
            continue;
        }

        atom_ljs.set(cgatomidx, atmtyp.ljParameter());
        atom_elements.set(cgatomidx, atmtyp.element());
        particle_type.set(cgatomidx, atmtyp.particleTypeString());
    }

    Properties props;

    props.setProperty("atomtype", atom_type);

    if (uses_bondtypes)
    {
        //this forcefield uses a different ato
        props.setProperty("bondtype", bond_type);
    }

    props.setProperty("charge_group", charge_group);
    props.setProperty("charge", atom_chgs);
    props.setProperty("mass", atom_masses);
    props.setProperty("LJ", atom_ljs);
    props.setProperty("element", atom_elements);
    props.setProperty("particle_type", particle_type);

    return std::make_tuple(props, errors);
}

/** This internal function is used to return all of the bond properties
    for the passed molecule */
GroTop::PropsAndErrors GroTop::getBondProperties(const MoleculeInfo &molinfo,
                                                 const GroMolType &moltype) const
{
    const auto R = InternalPotential::symbols().bond().r();

    QStringList errors;

    //add in all of the bond functions, together with the connectivity of the
    //molecule
    auto connectivity = Connectivity(molinfo).edit();
    connectivity = connectivity.disconnectAll();

    TwoAtomFunctions bondfuncs(molinfo);

    const auto bonds = moltype.bonds();

    for (auto it = bonds.constBegin(); it != bonds.constEnd(); ++it)
    {
        const auto &bond = it.key();
        auto potential = it.value();

        AtomIdx idx0 = molinfo.atomIdx( bond.atom0() );
        AtomIdx idx1 = molinfo.atomIdx( bond.atom1() );

        if (idx1 < idx0)
        {
            qSwap(idx0, idx1);
        }

        //do we need to resolve this bond parameter (look up the parameters)?
        if (not potential.isResolved())
        {
            //look up the atoms in the molecule template
            const auto atom0 = moltype.atom(idx0);
            const auto atom1 = moltype.atom(idx1);

            //get the bond parameter for these bond types
            auto new_potential = this->bond(atom0.bondType(), atom1.bondType(),
                                            potential.functionType());

            if (not new_potential.isResolved())
            {
                errors.append( QObject::tr("Cannot find the bond parameters for "
                  "the bond between atoms %1-%2 (atom types %3-%4, function type %5).")
                    .arg(atom0.toString()).arg(atom1.toString())
                    .arg(atom0.bondType()).arg(atom1.bondType())
                    .arg(potential.functionType()) );
                continue;
            }

            potential = new_potential;
        }

        //add the connection
        if (potential.atomsAreBonded())
        {
            connectivity.connect(idx0, idx1);
        }

        //now create the bond expression
        auto exp = potential.toExpression(R);

        if (not exp.isZero())
        {
            //add this expression onto any existing expression
            auto oldfunc = bondfuncs.potential(idx0, idx1);

            if (not oldfunc.isZero())
            {
                bondfuncs.set(idx0, idx1, exp+oldfunc);
            }
            else
            {
                bondfuncs.set(idx0, idx1, exp);
            }
        }
    }

    auto conn = connectivity.commit();

    Properties props;
    props.setProperty("connectivity", conn);
    props.setProperty("bond", bondfuncs);

    //if 'generate_pairs' is true, then we need to automatically generate
    //the excluded atom pairs, using fudge_qq and fudge_lj for the 1-4 interactions
    if (generate_pairs)
    {
        if (bonds.isEmpty())
        {
            //there are no bonds, so there cannot be any intramolecular nonbonded
            //energy (don't know how atoms are connected). This likely means
            //that this is a solvent molecule, so set the intrascales to 0
            CLJNBPairs nbpairs(molinfo, CLJScaleFactor(0));
            props.setProperty("intrascale", nbpairs);
        }
        else
        {
            CLJNBPairs nbpairs(conn, CLJScaleFactor(fudge_qq,fudge_lj));
            props.setProperty("intrascale", nbpairs);
        }
    }

    return std::make_tuple(props, errors);
}

/** This internal function is used to return all of the angle properties
    for the passed molecule */
GroTop::PropsAndErrors GroTop::getAngleProperties(const MoleculeInfo &molinfo,
                                                  const GroMolType &moltype) const
{
    const auto THETA = InternalPotential::symbols().angle().theta();

    QStringList errors;

    //add in all of the angle functions
    ThreeAtomFunctions angfuncs(molinfo);

    const auto angles = moltype.angles();

    for (auto it = angles.constBegin(); it != angles.constEnd(); ++it)
    {
        const auto &angle = it.key();
        auto potential = it.value();

        AtomIdx idx0 = molinfo.atomIdx( angle.atom0() );
        AtomIdx idx1 = molinfo.atomIdx( angle.atom1() );
        AtomIdx idx2 = molinfo.atomIdx( angle.atom2() );

        if (idx2 < idx0)
        {
            qSwap(idx0, idx2);
        }

        //do we need to resolve this angle parameter (look up the parameters)?
        if (not potential.isResolved())
        {
            //look up the atoms in the molecule template
            const auto atom0 = moltype.atom(idx0);
            const auto atom1 = moltype.atom(idx1);
            const auto atom2 = moltype.atom(idx2);

            //get the angle parameter for these atom types
            auto new_potential = this->angle(atom0.bondType(), atom1.bondType(),
                                             atom2.bondType(), potential.functionType());

            if (not new_potential.isResolved())
            {
                errors.append( QObject::tr("Cannot find the angle parameters for "
                  "the angle between atoms %1-%2-%3 (atom types %4-%5-%6, "
                  "function type %7).")
                    .arg(atom0.toString()).arg(atom1.toString()).arg(atom2.toString())
                    .arg(atom0.bondType()).arg(atom1.bondType()).arg(atom2.bondType())
                    .arg(potential.functionType()) );
                continue;
            }

            potential = new_potential;
        }

        //now create the angle expression
        auto exp = potential.toExpression(THETA);

        if (not exp.isZero())
        {
            //add this expression onto any existing expression
            auto oldfunc = angfuncs.potential(idx0, idx1, idx2);

            if (not oldfunc.isZero())
            {
                angfuncs.set(idx0, idx1, idx2, exp+oldfunc);
            }
            else
            {
                angfuncs.set(idx0, idx1, idx2, exp);
            }
        }
    }

    Properties props;
    props.setProperty("angle", angfuncs);

    return std::make_tuple(props, errors);
}

/** This internal function is used to return all of the dihedral properties
    for the passed molecule */
GroTop::PropsAndErrors GroTop::getDihedralProperties(const MoleculeInfo &molinfo,
                                                     const GroMolType &moltype) const
{
    const auto PHI = InternalPotential::symbols().dihedral().phi();

    QStringList errors;

    //add in all of the dihedral and improper functions
    FourAtomFunctions dihfuncs(molinfo);

    const auto dihedrals = moltype.dihedrals();

    for (auto it = dihedrals.constBegin(); it != dihedrals.constEnd(); ++it)
    {
        const auto &dihedral = it.key();
        auto potential = it.value();

        AtomIdx idx0 = molinfo.atomIdx( dihedral.atom0() );
        AtomIdx idx1 = molinfo.atomIdx( dihedral.atom1() );
        AtomIdx idx2 = molinfo.atomIdx( dihedral.atom2() );
        AtomIdx idx3 = molinfo.atomIdx( dihedral.atom3() );

        if (idx3 < idx0)
        {
            qSwap(idx0, idx3);
            qSwap(idx2, idx1);
        }

        Expression exp;

        //do we need to resolve this dihedral parameter (look up the parameters)?
        if (not potential.isResolved())
        {
            //look up the atoms in the molecule template
            const auto atom0 = moltype.atom(idx0);
            const auto atom1 = moltype.atom(idx1);
            const auto atom2 = moltype.atom(idx2);
            const auto atom3 = moltype.atom(idx3);

            //get the dihedral parameter for these atom types - could be
            //many, as they will be added together
            auto resolved = this->dihedrals(atom0.bondType(), atom1.bondType(),
                                            atom2.bondType(), atom3.bondType(),
                                            potential.functionType());

            if (resolved.isEmpty())
            {
                errors.append( QObject::tr("Cannot find the dihedral parameters for "
                  "the dihedral between atoms %1-%2-%3-%4 (atom types %5-%6-%7-%8, "
                  "function type %9).")
                    .arg(atom0.toString()).arg(atom1.toString())
                    .arg(atom2.toString()).arg(atom3.toString())
                    .arg(atom0.bondType()).arg(atom1.bondType())
                    .arg(atom2.bondType()).arg(atom3.bondType())
                    .arg(potential.functionType()) );
                continue;
            }

            //sum all of the parts together
            for (const auto r : resolved)
            {
                if (r.isResolved())
                {
                    exp += r.toExpression(PHI);
                }
            }
        }
        else
        {
            //we have a fully-resolved dihedral potential
            exp = potential.toExpression(PHI);
        }

        if (not exp.isZero())
        {
            //add this expression onto any existing expression
            auto oldfunc = dihfuncs.potential(idx0, idx1, idx2, idx3);

            if (not oldfunc.isZero())
            {
                dihfuncs.set(idx0, idx1, idx2, idx3, exp+oldfunc);
            }
            else
            {
                dihfuncs.set(idx0, idx1, idx2, idx3, exp);
            }
        }
    }

    Properties props;
    props.setProperty("dihedral", dihfuncs);

    return std::make_tuple(props, errors);
}

/** This function is used to create a molecule. Any errors should be written
    to the 'errors' QStringList passed as an argument */
Molecule GroTop::createMolecule(QString moltype_name, QStringList &errors,
                                const PropertyMap &map) const
{
    //find the molecular template for this molecule
    int idx = -1;

    for (int i=0; i<moltypes.count(); ++i)
    {
        if (moltypes.constData()[i].name() == moltype_name)
        {
            idx = i;
            break;
        }
    }

    if (idx == -1)
    {
        QStringList typs;

        for (const auto &moltype : moltypes)
        {
            typs.append( moltype.name() );
        }

        errors.append( QObject::tr("There is no molecular template called '%1' "
          "in this Gromacs file. Available templates are [ %2 ]")
            .arg(moltype_name).arg(Sire::toString(typs)) );
        return Molecule();
    }

    const auto moltype = moltypes.constData()[idx];

    //create the underlying molecule
    auto mol = this->createMolecule(moltype, errors, map).edit();
    mol.rename(moltype_name);
    const auto molinfo = mol.info();

    //now get all of the molecule properties
    const QVector< std::function<PropsAndErrors()> > funcs =
                { [&](){ return getAtomProperties(molinfo, moltype); },
                  [&](){ return getBondProperties(molinfo, moltype); },
                  [&](){ return getAngleProperties(molinfo, moltype); },
                  [&](){ return getDihedralProperties(molinfo, moltype); } };

    QVector<PropsAndErrors> props(funcs.count());
    auto props_data = props.data();

    if (usesParallel())
    {
        tbb::parallel_for( tbb::blocked_range<int>(0,funcs.count()),
                           [&](const tbb::blocked_range<int> &r)
        {
            for (int i=r.begin(); i<r.end(); ++i)
            {
                props_data[i] = funcs.at(i)();
            }
        });
    }
    else
    {
        for (int i=0; i<funcs.count(); ++i)
        {
            props_data[i] = funcs.at(i)();
        }
    }

    //assemble all of the properties together
    for (int i=0; i<props.count(); ++i)
    {
        const auto &p = std::get<0>(props.at(i));
        const auto &pe = std::get<1>(props.at(i));

        if (not pe.isEmpty())
        {
            errors += pe;
        }

        for (const auto key : p.propertyKeys())
        {
            const auto mapped = map[key];

            if (mapped.hasValue())
            {
                mol.setProperty(key, mapped.value());
            }
            else
            {
                mol.setProperty(mapped, p.property(key));
            }
        }
    }

    //finally set the forcefield property
    const auto mapped = map["forcefield"];

    if (mapped.hasValue())
    {
        mol.setProperty("forcefield", mapped.value());
    }
    else
    {
        mol.setProperty(mapped, moltype.forcefield());
    }

    return mol.commit();
}

/** Use the data contained in this parser to create a new System of molecules,
    assigning properties based on the mapping in 'map' */
System GroTop::startSystem(const PropertyMap &map) const
{
    if (grosys.isEmpty())
    {
        //there are no molecules to process
        return System();
    }

    //first, create template molecules for each of the unique molecule types
    const auto unique_typs = grosys.uniqueTypes();

    QHash<QString,Molecule> mol_templates;
    QHash<QString,QStringList> template_errors;
    mol_templates.reserve(unique_typs.count());

    //loop over each unique type, creating the associated molecule and storing
    //in mol_templates. If there are any errors, then store them in template_errors
    if (usesParallel())
    {
        QMutex mutex;

        tbb::parallel_for( tbb::blocked_range<int>(0,unique_typs.count()),
                           [&](const tbb::blocked_range<int> &r)
        {
            for (int i=r.begin(); i<r.end(); ++i)
            {
                auto typ = unique_typs[i];

                QStringList errors;
                Molecule mol = this->createMolecule(typ, errors, map);

                QMutexLocker lkr(&mutex);

                if (not errors.isEmpty())
                {
                    template_errors.insert(typ, errors);
                }

                mol_templates.insert(typ, mol);
            }
        });
    }
    else
    {
        for (auto typ : unique_typs)
        {
            QStringList errors;
            Molecule mol = this->createMolecule(typ, errors, map);

            if (not errors.isEmpty())
            {
                template_errors.insert(typ, errors);
            }

            mol_templates.insert(typ,mol);
        }
    }

    //next, we see if there were any errors. If there are, then raise an exception
    if (not template_errors.isEmpty())
    {
        QStringList errors;

        for (auto it = template_errors.constBegin();
             it != template_errors.constEnd(); ++it)
        {
            errors.append( QObject::tr("Error constructing the molecule associated with "
               "template '%1' : %2").arg(it.key()).arg(it.value().join("\n")) );
        }

        throw SireIO::parse_error( QObject::tr(
                "Could not construct a molecule system from the information stored "
                "in this Gromacs topology file. Errors include:\n%1")
                    .arg(errors.join("\n\n")), CODELOC );
    }

    //next, make sure that none of the molecules are empty...
    {
        QStringList errors;

        for (auto it = mol_templates.constBegin();
             it != mol_templates.constEnd(); ++it)
        {
            if (it.value().isNull())
            {
                errors.append( QObject::tr("Error constructing the molecule associated with "
                  "template '%1' : The molecule is empty!")
                    .arg(it.key()) );
            }
        }

        if (not errors.isEmpty())
            throw SireIO::parse_error( QObject::tr(
                "Could not construct a molecule system from the information stored "
                "in this Gromacs topology file. Errors include:\n%1")
                    .arg(errors.join("\n\n")), CODELOC );
    }

    //now that we have the molecules, we just need to duplicate them
    //the correct number of times to create the full system
    MoleculeGroup molgroup("all");

    for (int i=0; i<grosys.nMolecules(); ++i)
    {
        molgroup.add( mol_templates.value(grosys[i]).edit().renumber() );
    }

    System system(grosys.name());
    system.add(molgroup);
    system.setProperty(map["fileformat"].source(), StringProperty(this->formatName()));

    return system;
}
