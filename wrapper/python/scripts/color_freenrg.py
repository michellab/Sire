description = """
color_freenrg is used to extract the free energy component outputs from a waterswap, ligandswap or proteinswap simulation
and to use them to color-code residues in a protein structure file
"""

import Sire.IO
import Sire.Mol
import Sire.Maths
import Sire.Config

import argparse
import os
import sys
import glob

pd = Sire.try_import("pandas")
from pandas import Series, DataFrame

parser = argparse.ArgumentParser(description="Extract residue-component free energies from Xswap simulations "
                                             "and use them to color-code residues in a protein structure",
                                 epilog="color_freenrg is built using Sire and is distributed "
                                        "under the GPL. For more information please "
                                        "type "
                                        "'color_freenrg --description'",
                                 prog="color_freenrg")

parser.add_argument('--description', action="store_true",
                    help="Print a complete description of this program.")

parser.add_argument('--author', action="store_true",
                    help="Get information about the authors of this script.")

parser.add_argument('--version', action="store_true",
                    help="Get version information about this script.")

parser.add_argument('-d', '--directory', nargs="?",
                    help='Full path to the directory containing the simulation output')

parser.add_argument('-r', '--range', nargs=2,
                    help="Supply the range of iterations over which to average. "
                         "By default, this will be over the last 60 percent of iterations.")

sys.stdout.write("\n")
args = parser.parse_args()

must_exit = False

if args.description:
    print("%s\n" % description)
    must_exit = True

if args.author:
    print("\ncolor_freenrg was written by Christopher Woods (C) 2018")
    must_exit = True

if args.version:
    print("color_freenrg -- from Sire release version <%s>" % Sire.__version__)
    print("This particular release can be downloaded here: "
          "https://github.com/michellab/Sire/releases/tag/v%s" % Sire.__version__)
    must_exit = True


outdir = args.directory

if outdir is None:
    must_exit = True

if must_exit:
    sys.exit(0)
    
if args.range:
    range_start = int(args.range[0])
    range_end = int(args.range[1])
    if range_end < 1:
        range_end = 1
    if range_start < 1:
        range_start = 1
    if range_start > range_end:
        tmp = range_start
        range_start = range_end
        range_end = tmp
else:
    range_start = None
    range_end = None

percent = 60.0

print("Calculating free energy components from the simulation data in %s" % (outdir))

# find all files and work out the highest iteration
filenames = glob.glob("%s/results_????.log" % outdir)

maxit = 0
for filename in filenames:
    it = int(filename[-8:-4])
    if it > maxit:
        maxit = it

if range_start:
    if range_end > maxit:
        range_end = maxit
else:
    range_end = maxit
    range_start = int(0.01*(100-percent)*maxit)

print("Creating average from iterations %s to %s (inclusive)" % (range_start,range_end))

# create the list of filenames to process
filenames = ["%s/results_%04d.log" % (outdir,i) for i in range(range_start,range_end+1)]

def getComponents(filenames):
    """Read all of the residue-based free energy components from the log files produced
       by a waterswap or ligandswap simulation (passed as a list of filenames). Return
       the average components as a pandas DataFrame"""
    avgs = {}
    resids = {}
    
    # Loop over all of the files...
    for filename in filenames:
        has_started=False
        for line in open(filename).readlines():
            # Read from the line "RESIDUE FREE ENERGY COMPONENTS" onwards...
            if line.find("RESIDUE FREE ENERGY COMPONENTS") != -1:
                has_started = True
            
            elif has_started:
                words = line.split()
                if len(words) == 8:
                    resname = words[1]
                    resnum = int(words[3])
                    total = float(words[-3])
                    coul = float(words[-2])
                    lj = float(words[-1])
                    key = "%s:%s" % (resname,resnum)
                    
                    if not key in avgs:
                        avgs[key] = [Sire.Maths.Average(), Sire.Maths.Average(), Sire.Maths.Average()]
                        if not resnum in resids:
                            resids[resnum] = [resname]
                        else:
                            resids[resnum].append(resname)
                    
                    # accumulate the average total, coulomb and LJ free energies
                    avgs[key][0].accumulate(total)
                    avgs[key][1].accumulate(coul)
                    avgs[key][2].accumulate(lj)
                    
                elif line.find("COMPONENTS") != -1:
                    break
    
    # Now sort the data into a pandas DataFrame
    resnums = list(resids.keys())
    resnums.sort()
    resnams = []
    total = []
    coul = []
    lj = []
    
    for resnum in resnums:
        for resname in resids[resnum]:
            key = "%s:%s" % (resname,resnum)
            avg = avgs[key]
            resnams.append(resname)
            total.append(avg[0].average())
            coul.append(avg[1].average())
            lj.append(avg[2].average())
    
    # The data is in lists which can be put into pandas columns. We will index the 
    # DataFrame using the residue number (assuming that they are all unique)
    return DataFrame( index = resnums,
                      data = {"name" : resnams, "total" : total, "coulomb" : coul, "LJ" : lj},
                      columns=["name", "total", "coulomb", "LJ"] )

components = getComponents(filenames)

print("\nAll of the average components are listed below")
with pd.option_context('display.max_rows', None):
    print(components)
print("\n", end="")

# now color-code all of the bound_mobile*.pdb files that we find, writing the output
# to "bound_mobile_*_component.pdb"

def colorProtein(protein, data, column):
    """Color-code the passed protein using the data contained in the passed dataframe, using the
       specified column"""
    
    # first find the maximum absolute value - we will scale linearly from there
    vals = data[column]
    maxval = vals.abs().max()
    
    # now create an AtomFloatProperty that will contain a number for each atom
    # in each residue. This will be from 0-100, with 0 representing -maxval, 
    # 50 representing 0 and 100 representing maxval
    betas = Sire.Mol.AtomFloatProperty(protein, 50.0)
    
    for x in data.index:
        resnum = Sire.Mol.ResNum(int(x))
        resnam = Sire.Mol.ResName(data.name[x])
        value = vals[x]
        
        scaled = 50.0 + 50.0*(value/maxval)
        
        # issues with beta mean it must lie between 0 and 99.99
        if scaled < 0:
            scaled = 0.0
        elif scaled > 99.99:
            scaled = 99.99

        residue = protein[ resnam + resnum ]
        
        for atom in residue.atoms():
            betas.set(atom.cgAtomIdx(), scaled)
            
    # Set the 'beta-factor' property as this is the name used for the 'beta_factor'
    # value by the PDB writer
    protein = protein.edit().setProperty("beta_factor", betas).commit()
    return protein

for filename in glob.glob("%s/bound_mobile*.pdb" % outdir):
    if filename.endswith("_total.pdb") or filename.endswith("_coulomb.pdb") or filename.endswith("_LJ.pdb"):
        continue

    print("Color-coding file %s..." % filename, end="")

    system = Sire.IO.MoleculeParser.read(filename)
    protein = system[Sire.Mol.MolWithResID(Sire.Mol.ResName("ALA"))[0]]

    for c in ["total", "coulomb", "LJ"]:
        protein = colorProtein(protein, components, c)
        system.update(protein)

        newname = filename.replace(".pdb", "_%s.pdb" % c)
        Sire.IO.PDB2(system).writeToFile(newname)
        print("%s written to %s..." % (c,newname), end="")
    
    print("\n")

