description = """
transform is a simple app that is used to perform geometric transformations on a molecule taken from a PDB file.
Assuming you have a PDB file, “sys.pdb” that contains a molecule with residue name “MOL” then;
sire.app/bin/transform -p sys.pdb -l MOL --translate 0*nanometer 3*nanometer -0.5*nanometer -o output.pdb
will translate MOL by 0*nanometers along the x axis, 3 nanometers on the y axis and -0.5 nanometers on the z axis and will write the result to “output.pdb”. transform recognises all of the length units recognised by Sire (although translating by meters is perhaps a bit over the top!).
You can also rotate molecules, using the “--rotate” option. By default this will use the center of mass (via the “--com” option) or center of geometry (via the “--cog” option), or you can manually specify the center using the “--rotcent” option. Rotations can be specified in degrees or radians.If you supply both a translation and rotation, then the rotation is performed first, and the translation is performed second.
If you need more help understanding or using align then please feel free to get in touch via the Sire users mailing list.
"""

from Sire.IO import *
from Sire.Mol import *
from Sire.Maths import *
from Sire.Base import *
from Sire.Units import *

import Sire.Config

import argparse
import os
import sys

parser = argparse.ArgumentParser(description="Transform ligands/small molecules/fragments "
                                             "in a PDB file (e.g. translate/rotate)",
                                 epilog="This program is built using Sire, which is released "
                                        "under the GPL. For more information please visit "
                                        "http://siremol.org/transform",
                                 prog="transform")

parser.add_argument('--description', action="store_true",
                    help="Print a complete description of this program.")

parser.add_argument('--author', action="store_true",
                    help="Get information about the authors of this script.")

parser.add_argument('--version', action="store_true",
                    help="Get version information about this script.")

parser.add_argument('-l', '--ligand', nargs=1,
                    help="Supply the name of the ligand you want to translate/rotate.")

parser.add_argument('-p', '--pdb', nargs=1,
                    help="Supply the PDB file containing the ligand.")

parser.add_argument('-t', '--translate', nargs=3,
                    help="How much to translate the ligand along x, then y, then z, "
                         "e.g. 0*nanometer 3*nanometer -0.5*nanometer. If no units "
                         "are supplied, then they are assumed to be angstroms.")

parser.add_argument('-r', '--rotate', nargs=1,
                    help="How must to rotate the ligand around the rotation axis, e.g. "
                         "30*degrees. If no units are supplied, then they are assumed "
                         "to be degrees.")

parser.add_argument('-ra', '--rotaxis', type=float, nargs=3,
                    help="The rotation axis about which to rotate the molecule, e.g. "
                         "1.0 0.0 0.0. If this is not supplied, then the rotation will "
                         "be about the major principal axis of the molecule.")

parser.add_argument('-rc', '--rotcent', nargs=3,
                    help="The center of rotation (before translation) for the rotation, e.g. "
                         "5*angstrom 5.2*angstrom 4.8*angstrom. If no units are supplied, then "
                         "angstroms are assumed. If this is not specified, then rotation will "
                         "be about the center of mass or center of geometry of the molecule. "
                         "(see --com or --cog).")

parser.add_argument('--cog', action="store_true",
                    help="Whether or not to rotate about the center of geometry of the molecule.")

parser.add_argument('--com', action="store_true",
                    help="Whether or not to rotate about the center of mass of the molecule.")

parser.add_argument('-o', '--output', nargs=1,
                    help="Name of the PDB file in which to output the transformed copy of the "
                         "ligand.")

sys.stdout.write("\n")
args = parser.parse_args()

must_exit = False

if args.description:
    print("%s\n" % description)
    must_exit = True

if args.author:
    print("\ntransform was written by Christopher Woods (C) 2014")
    print("It is based on the Molecule::move().transform() function distributed in Sire.")
    must_exit = True

if args.version:
    print("\ntransform version 0.2")
    print(Sire.Config.versionString())
    must_exit = True

if must_exit:
    sys.exit(0)

pdb = args.pdb
ligname = args.ligand
outfile = args.output

if pdb is None or ligname is None or outfile is None:
    parser.print_help()
    sys.exit(0)

pdb = args.pdb[0]
ligname = args.ligand[0]
outfile = args.output[0]

delta = Vector(0,0,0)
rotvec = None
rotcent = None
rotang = 0
axstring = "major principal axis"

use_cog = True
use_com = False
rotstring = "center of geometry"

def stringToValue(val):
    # see if we have to turn the value from a string into a python object
    try:
        return eval(val).value()
    except:
        try:
            words = val.split("*")
            return float(words[0])
        except:
            return float(val)

def stringToAngle(val):
    # see if we have to turn the value from a string into a python object
    try:
        return eval(val).to(degrees)
    except:
        words = val.split("*")
        return float(words[0])

if args.translate:
    dx = args.translate[0]
    dy = args.translate[1]
    dz = args.translate[2]

    delta = Vector(stringToValue(dx), stringToValue(dy), stringToValue(dz))

if args.rotate:
    rotang = stringToAngle(args.rotate[0])

if args.rotaxis:
    x = stringToValue(args.rotaxis[0])
    y = stringToValue(args.rotaxis[1])
    z = stringToValue(args.rotaxis[2])

    rotvec = Vector(x,y,z)
    axstring = "axis %s" % rotvec
else:
    axstring = "major principal axis"

if args.rotcent:
    x = stringToValue(args.rotcent[0])
    y = stringToValue(args.rotcent[1])
    z = stringToValue(args.rotcent[2])

    rotcent = Vector(x,y,z)                   
    rotstring = "axis %s A" % rotcent

elif args.com:
    use_com = True
    use_cog = False
    rotstring = "center of mass"

elif args.cog:
    use_cog = True
    use_com = False
    rotstring = "center of geometry"

print("\nRotating ligand %s from PDB file %s by %s degrees around %s, "
      "about %s, and then translating by %s A." % (ligname,pdb,rotang,axstring,
                 rotstring, delta))

print("\nTransformed coordinates of ligand %s will be written to file %s." % (ligname,outfile))

mols = PDB().read(pdb)

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

lig = findMolecule(mols, ligname)

if lig is None:
    print("\nWARNING: Cannot find the ligand (%s) in PDB file %s." % (ligname,pdb))
    print("Cannot transform. Exiting...")
    sys.exit(-1)

if rotcent is None:
    if use_com:
        rotcent = lig.evaluate().centerOfMass()
        print("Center of mass is %s" % rotcent)
    else:
        rotcent = lig.evaluate().centerOfGeometry()
        print("Center of geometry is %s" % rotcent)

if rotvec is None:
    axes = lig.evaluate().principalAxes()
    rotvec = axes.matrix().column0()
    print("Major principal axis is %s" % rotvec)

t = Transform(delta, Quaternion(rotang*degrees,rotvec), rotcent)

print("Transform equals %s" % t.toString())

lig = lig.move().transform(t).commit()

PDB().write(lig, outfile)
