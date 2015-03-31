
from Sire.System import *
from Sire.IO import *
from Sire.FF import *
from Sire.Maths import *
from Sire.Mol import *
from Sire.Vol import *
from Sire.Units import *
from Sire.MM import *
from Sire.Move import *

print("Loading the molecules...")
solute = PDB().readMolecule("test/io/methanol.pdb")

solute_params = { 3010 : ( 0.117, 3.3997, 0.1094),
                  3011 : ( 0.028, 2.4714, 0.0157),
                  3012 : (-0.598, 3.0665, 0.2104),
                  3013 : ( 0.397, 0.0000, 0.0000) }

solute_template = { 3010 : "C01",
                    3011 : "H02:H03:H04",
                    3012 : "O05",
                    3013 : "H06" }

solvent_params = { 3000 : ( 0.000,  3.15365, 0.1550),
                   3001 : ( 0.520,  0.00000, 0.0000),
                   3002 : (-1.040,  0.00000, 0.0000) }

solvent_template = { 3000 : "O00",
                     3001 : "H01:H02",
                     3002 : "M03" }

identity_atoms = ("C01", "H02", "H03", "H04", "O05", "H06")

for t in solute_template:
    atoms = solute_template[t].split(":")

    for atom in atoms:
        charge = solute_params[t][0] * mod_electron
        lj = LJParameter(solute_params[t][1] * angstrom, solute_params[t][2] * kcal_per_mol)

        solute = solute.atom( AtomName(atom)) \
                                  .edit().setProperty("charge", charge) \
                                         .setProperty("LJ", lj) \
                                         .molecule().commit()

solute = MoleculeGroup("solute", solute)

mols = PDB().read("test/io/water.pdb")

mol = mols.moleculeAt(0).molecule()

for t in solvent_template:
    atoms = solvent_template[t].split(":")

    for atom in atoms:
        charge = solvent_params[t][0] * mod_electron
        lj = LJParameter(solvent_params[t][1] * angstrom, solvent_params[t][2] * kcal_per_mol)

        mol = mol.atom( AtomName(atom)) \
                             .edit().setProperty("charge", charge) \
                                    .setProperty("LJ", lj) \
                                    .molecule().commit()

solvent_charges = mol.property("charge")
solvent_ljs = mol.property("LJ")

solvent = MoleculeGroup("solvent", mols)

for molnum in mols.molNums():
    mol = mols[molnum].molecule().edit().setProperty("charge", solvent_charges) \
                                        .setProperty("LJ", solvent_ljs).commit()

    solvent.update(mol)

system = System()
system.add(solute)
system.add(solvent)

solvent_cljff = InterCLJFF("solvent_cljff")
solvent_cljff.add(solvent)
solute_solvent_cljff = InterGroupCLJFF("solute_solvent_cljff")
solute_solvent_cljff.add( solute, MGIdx(0) )
solute_solvent_cljff.add( solvent, MGIdx(1) )

system.add(solvent_cljff)
system.add(solute_solvent_cljff)

system.setProperty("space", PeriodicBox( Vector(-18.3854,-18.66855,-18.4445),
                                         Vector( 18.3854, 18.66855, 18.4445) ) )

print("Equilibrating the system...")

solvent_move = RigidBodyMC( PrefSampler(solute.moleculeAt(0), solvent, 200*angstrom2) )
solvent_move.setMaximumTranslation( 0.2 * angstrom )
solvent_move.setMaximumRotation( 5 * degrees )

solute_move = RigidBodyMC( solute )
solute_move.setMaximumTranslation( 0.2 * angstrom )
solute_move.setMaximumRotation( 5 * degrees )

moves = WeightedMoves()
moves.add( solvent_move, 100 )
moves.add( solute_move, 1 )

for i in range(1,11):
    system = moves.move(system, 5000, False )
    #system = moves.move(system, 50, False)
    print("Step %d of 10: Energy = %s" % (i, system.energy()))

    PDB().write(system.molecules(), "test_equil_%0004d.pdb" % i)

print("Equilibration complete")

print("Adding energy monitors...")

identity_points = []

for atom in identity_atoms:
    identity_points.append( AtomPoint( solute.moleculeAt(0).atom(AtomName(atom)) ) )

idassigner = IDAssigner(identity_points, solvent)

nrgmon0 = EnergyMonitor(solute, solvent)
nrgmon1 = EnergyMonitor(solute, idassigner)

system.add("solute_solvent", nrgmon0, 100)
system.add("solute_shell", nrgmon1, 100)


def getBeta(nrg):

    nrg = 50.0 + (20.0 * nrg)

    if nrg < 1:
        nrg = 1.0
    elif nrg > 99:
        nrg = 99.0

    return nrg

def view(nrgmon, id):

    views0 = nrgmon.views0()

    coul_group = MoleculeGroup("coulomb")
    lj_group = MoleculeGroup("lj")
    total_group = MoleculeGroup("total")

    for view in views0:
        coul_group.add(view)
        lj_group.add(view)
        total_group.add(view)

    views1 = nrgmon.views1()
    coul_nrgs = nrgmon.coulombEnergies()
    lj_nrgs = nrgmon.ljEnergies()

    for view in views1:
        coul_group.add(view)
        lj_group.add(view)
        total_group.add(view)

    cnrg = 0
    ljnrg = 0

    for i in range(0, len(views1)):
        this_cnrg = coul_nrgs(0,i).average()
        this_ljnrg = lj_nrgs(0,i).average()

        print(views1[i], this_cnrg, this_ljnrg, (this_cnrg+this_ljnrg))

        beta_cnrg = getBeta(this_cnrg)
        beta_ljnrg = getBeta(this_ljnrg)
        beta_total = getBeta(this_cnrg + this_ljnrg)

        print("(",beta_cnrg,beta_ljnrg,beta_total,")")

        cnrg += coul_nrgs(0,i).average()
        ljnrg += lj_nrgs(0,i).average()

        coul_mol = views1[i].molecule()
        lj_mol = views1[i].molecule()
        total_mol = views1[i].molecule()

        atoms = views1[i].atoms()

        for j in range(0,atoms.count()):
            coul_mol = coul_mol.atom(atoms[j].index()).edit() \
                               .setProperty("b-factor",beta_cnrg).molecule().commit()

            lj_mol = lj_mol.atom(atoms[j].index()).edit() \
                           .setProperty("b-factor",beta_ljnrg).molecule().commit()

            total_mol = total_mol.atom(atoms[j].index()).edit() \
                                 .setProperty("b-factor",beta_total).molecule().commit()

        coul_group.update(coul_mol)
        lj_group.update(lj_mol)
        total_group.update(total_mol)

    print("Total energies: %f  %f  %f" % (cnrg, ljnrg, cnrg+ljnrg))

    PDB().write(coul_group, "test_coul_%0004d.pdb" % id)
    PDB().write(lj_group, "test_lj_%0004d.pdb" % id)
    PDB().write(total_group, "test_total_%0004d.pdb" % id)

print("Running a simulation...")

for i in range(1,11):
    system = moves.move(system, 2000, True)

    nrgmon = system.monitor( MonitorName("solute_shell") )

    PDB().write(system.molecules(), "test_sim_%0004d.pdb" % i)

    system.clearStatistics()

    print("Step %d of 10: Energy = %s" % (i, system.energy()))

    view(nrgmon, i)
