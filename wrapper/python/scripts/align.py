
from Sire.IO import *
from Sire.Mol import *
from Sire.Maths import *
description = """
align is a simple app that is used to align one molecule against another. Alignment is based on a graph theory analysis of the bonding of the two molecules, with a maximum common substructure algorithm used to match atoms and to then find the optimum overlap.
Assuming you have two PDB files, “sys0.pdb”, which contains a molecule with residue name “MOL0”, and “sys1.pdb” that contains a molecule with residue name “MOL1”, then;
sire.app/bin/align -p0 sys0.pdb -p1 sys1.pdb -l0 MOL0 -l1 MOL1 -o output.pdb
will find align molecule MOL1 against molecule MOL0 and will write the result to “output.pdb”. 
There are several options that you can use to control alignment. For example, by default, light atoms are excluded from the match (and alignment). To include them, use the “--match-light” option.
The maximum common substructure algorithm has been optimised, but it can be slow (it is a polynomial time algorithm). You can specify the maximum timeout for the match (in seconds) using the “--timeout” option. You can also help the match by manually specifying the names of equivalent atoms, e.g. 
sire.app/bin/align -p0 sys0.pdb -p1 sys1.pdb -l0 MOL0 -l1 MOL1 -m CA:CB CD:CE CG:CH -o output.pdb
tells the app that the atom called “CB” in MOL1 maps to atom “CA” in MOL0, “CE” maps to “CD” and “CH” maps to “CG”.If you need more help understanding or using align then please feel free to get in touch via the Sire users mailing list.
"""

from Sire.Base import *
from Sire.Units import *

import Sire.Config

import argparse
import os
import sys

parser = argparse.ArgumentParser(description="Align ligands/small molecules/fragments against each other.",
                                 epilog="align is built using Sire and is distributed "
                                        "under the GPL. For more information please visit "
                                        "http://siremol.org/align",
                                 prog="align")

parser.add_argument('--description', action="store_true",
                    help="Print a complete description of this program.")

parser.add_argument('--author', action="store_true",
                    help="Get information about the authors of this script.")

parser.add_argument('--version', action="store_true",
                    help="Get version information about this script.")

parser.add_argument('-l0', '--ligand0', nargs=1,
                    help="Supply the name of the ligand in the first PDB file against "
                         "which the ligand in the second PDB file will be aligned.")

parser.add_argument('-l1', '--ligand1', nargs=1,
                    help="Supply the name of the ligand in the second PDB file that will "
                         "be aligned against the ligand in the first PDB file.")

parser.add_argument('-p0', '--pdb0', nargs=1,
                    help="Supply the PDB file containing the first ligand against which "
                         "the second ligand will be aligned.")

parser.add_argument('-p1', '--pdb1', nargs=1,
                    help="Supply the PDB file containing the second ligand, which will be "
                         "aligned against the first ligand.")

parser.add_argument('-t', '--timeout', nargs=1,
                    help="Set the timeout (in seconds) used when searching for the maximum common overlap of "
                         "the molecules. The longer the timeout, the more likely you will find the "
                         "best match.")

parser.add_argument('-m', '--match', nargs="+",
                    help="Manually specify matching atoms in the molecules. You can do this using "
                         "atom name, e.g. -m CA:CB means that the atom called \"CA\" in the first ligand "
                         "is called \"CB\" in the second ligand. You can specify multiple matches using spaces, "
                         "e.g. -m CA:CB CD:CE CG:CH etc.")

parser.add_argument('--match-light', action="store_true",
                    help="Flag that alignment should include matching hydrogen and other light atoms.")

parser.add_argument('-o', '--output', nargs=1,
                    help="Name of the PDB file in which to output the aligned copy of the "
                         "second ligand.")

sys.stdout.write("\n")
args = parser.parse_args()

must_exit = False

if args.description:
    print("%s\n" % description)
    must_exit = True

if args.author:
    print("\nalign was written by Christopher Woods (C) 2014")
    print("It is based on the Molecule::move().align() function distributed in Sire.")
    must_exit = True

if args.version:
    print("\nalign version 0.2")
    print(Sire.Config.versionString())
    must_exit = True

if must_exit:
    sys.exit(0)

if args.timeout:
    timeout = float(args.timeout[0]) * second
else:
    timeout = 1 * second

if args.match_light:
    match_light = True
else:
    match_light = False

if args.match:
    matcher = {}

    for m in args.match:
        parts = m.split(":")
        if len(parts) == 2:
            matcher[parts[0]] = parts[1]

    if len(matcher) > 0:
        matcher = AtomIDMatcher(matcher)
    else:
        matcher = None
else:
    matcher = None

pdb0 = args.pdb0
pdb1 = args.pdb1
ligname0 = args.ligand0
ligname1 = args.ligand1
outfile = args.output

if pdb0 is None or pdb1 is None or ligname0 is None or ligname1 is None or outfile is None:
    parser.print_help()
    sys.exit(0)

pdb0 = args.pdb0[0]
pdb1 = args.pdb1[0]
ligname0 = args.ligand0[0]
ligname1 = args.ligand1[0]
outfile = args.output[0]

print("\nAligning ligand %s from PDB file %s against ligand %s in PDB file %s." % \
            (ligname1,pdb1,ligname0,pdb0))
print("Aligned coordinates of ligand %s will be written to file %s." % (ligname1,outfile))
print("Match will use a timeout of %s s" % timeout.to(second))

if match_light:
    print("Match will include light (hydrogen) atoms.")

if matcher:
    print("Match will enforce that atoms are assigned using: %s" % matcher)

mols0 = PDB().read(pdb0)
mols1 = PDB().read(pdb1)

def getResidueNames(molecule):
    nres = molecule.nResidues()

    resnams = []

    for i in range(0, nres):
        resnams.append( str( molecule.residue(ResIdx(i)).name().value()).upper() )

    return resnams

def findMolecule(molecules, molname):
    molname = molname.upper()

    for molnum in molecules.molNums():
        molecule = molecules[molnum].molecule()

        if str(molecule.name().value()).upper() == molname:
            return molecule

        resnams = getResidueNames(molecule)

        for resnam in resnams:
            if resnam == molname:
                return molecule

    return None

lig0 = findMolecule(mols0, ligname0)
lig1 = findMolecule(mols1, ligname1)

can_align = True

if lig0 is None:
    print("\nWARNING: Cannot find the first ligand (%s) in PDB file %s." % (ligname0,pdb0))
    can_align = False

if lig1 is None:
    print("\nWARNING: Cannot find the second ligand (%s) in PDB file %s." % (ligname1,pdb1))
    can_align = False

if not can_align:
    print("Cannot align the molecule! Exiting...")
    sys.exit(-1)

print("Looking for the maximum common substructure of the two ligands...")

if matcher:
    atommap = AtomMCSMatcher(AtomMatchInverter(matcher),timeout,match_light) \
                    .match(lig1, lig0)
else:
    atommap = AtomMCSMatcher(timeout,match_light).match(lig1, lig0)

keys = atommap.keys()

lines = []

for key in keys:
    lines.append("%s <=> %s" % (lig1.atom(key).name(),lig0.atom(atommap[key]).name()))

lines.sort()

print("\nMapping from %s to %s" % (ligname1,ligname0))
for line in lines:
    print(line)

print("\nNumber of matched atoms == %s" % len(atommap))

print("\nRMSD before alignment == %s" % lig1.evaluate().rmsd(lig0, AtomResultMatcher(atommap)))

print("\nAligning the molecule...")
lig2 = lig1.move().align(lig0, AtomResultMatcher(atommap)).commit()

print("\nRMSD after alignment == %s" % lig2.evaluate().rmsd(lig0, AtomResultMatcher(atommap)))

PDB().write(lig2, outfile)
