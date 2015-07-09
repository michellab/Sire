from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout

from chemistry.amber import *
from chemistry import *

import os,re

#Per tutti i amber_inps file che ho 
#Creo un sistem 
# le force  le creo prima di tutto
# quindi per la riga che leggo 
# aggiungo i restraint e minimizzo
# salvo
# ricomincio il ciclo


###Natural sorting functions####
def try_int(s):
    "Convert to integer if possible."
    try: return int(s)
    except: return s

def natsort_key(s):
    "Used internally to get a tuple by which s is sorted."
    import re
    return map(try_int, re.findall(r'(\d+|\D+)', s))

def natcmp(a, b):
    "Natural string comparison, case sensitive."
    return cmp(natsort_key(a), natsort_key(b))

def natcasecmp(a, b):
    "Natural string comparison, ignores case."
    return natcmp(a.lower(), b.lower())

def natsort(seq, cmp=natcmp):
    "In-place natural string sort."
    seq.sort(cmp)
    
def natsorted(seq, cmp=natcmp):
    "Returns a copy of seq, sorted by natural string sort."
    import copy
    temp = copy.copy(seq)
    natsort(temp, cmp)
    return temp

###########Coordinate generation function ###########

def generateCoordFiles(outputfile, coordinates,outputfolder):
     

     
     string_rev   = re.sub("[^0-9.e-]", " ", coordinates)
     string_split = string_rev.split()
     

     file_dest = outputfolder + "/" + outputfile[0] + "." + outputfile[1]
     singlecoord = open(file_dest, "w")
     singlecoord.write("LIG\n")
     singlecoord.write("    32\n")

     small_file_index=0

     for f in string_split:
         if float(f)<0:
             singlecoord.write("  %.7f" % float(f))   #absolutely important to have 7 decimals otherwise openmm minimizer does not work!!
         elif float(f)>=10:
             singlecoord.write("  %.7f" % float(f))
         else:
             singlecoord.write("   %.7f" % float(f))
         small_file_index+=1
         if not small_file_index%6:
             singlecoord.write("\n")

#####################################

####Main code. Here all the input files are read####

output_pdb_folder = "amber_out_pdb"
output_crd_folder = "amber_out_crd"
if not os.path.exists(output_pdb_folder):
    os.makedirs(output_pdb_folder)
if not os.path.exists(output_crd_folder):
    os.makedirs(output_crd_folder)

input_path = os.getcwd()    
input_files = []


for root, dirs, files in os.walk(input_path):
    if "amber_inps" in root:
        for name in files:
            if "alanine" in name:
                #path = root + "/" + name
                
                input_files.append(name)

sort_list=natsorted(input_files) 


atom_list= open("atom_list.dat", "r")
readfile = atom_list.readlines()
index_file = 0 

platform = Platform.getPlatformByName("Reference")
 
for reading in readfile:

    for files in range(0,(len(sort_list)/len(readfile))): # 3 is equal to len(sort_list)/len(readfile)

        inputfile = "amber_inps/" + sort_list[index_file]
        base = AmberParm("a.prmtop", inputfile)  #then  for cycle and open everything in amber_inps
        index_file+=1
        print("Loading...")
        system = base.createSystem(nonbondedMethod=CutoffNonPeriodic, nonbondedCutoff=1*nanometer,constraints=HBonds) #HBONDS?  no barostat since we are in vacuum
        ##custom force must be declare before creation of simulation object####


        
        force = CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        force.addGlobalParameter("k", 100*kilocalories_per_mole/angstrom**2)
        force.addPerParticleParameter("x0")
        force.addPerParticleParameter("y0")
        force.addPerParticleParameter("z0")


        single_atom = reading.split()
      
        for i, atom_coord in enumerate(base.positions):
            
            if(base.atoms[i].name==single_atom[0] and base.atoms[i].idx == int(single_atom[1])) or \
            (base.atoms[i].name==single_atom[2] and base.atoms[i].idx == int(single_atom[3])) or \
            (base.atoms[i].name==single_atom[4] and base.atoms[i].idx == int(single_atom[5])) or \
            (base.atoms[i].name==single_atom[6] and base.atoms[i].idx == int(single_atom[7])):
                  print("Constraining atoms %s" % base.atoms[i])
                  force.addParticle(i, atom_coord.value_in_unit(nanometers)) #atom_coord.valuevalue_in_unit(u.nanometers)
        system.addForce(force)

        integrator = VerletIntegrator(2.0*femtoseconds)# try leapfrog
        simulation = Simulation(base.topology, system, integrator,platform) # creation of a simulation object
        simulation.context.setPositions(base.positions) #set the particle positions



        #positions = simulation.context.getState(getPositions=True).getPositions()
        print("Minimizing")
        simulation.minimizeEnergy()
        print("Finished minimization.")

        print("Saving PDB")
        positions = simulation.context.getState(getPositions=True).getPositions()
        name_out_file =sort_list[index_file].split(".")
        out_file = output_pdb_folder + "/" + name_out_file[0] + ".pdb"
        PDBFile.writeFile(simulation.topology,positions,open(out_file, "w"))
    
        print("Saving mdcrd")
        coord_string = str(positions.value_in_unit(angstrom))
        generateCoordFiles(name_out_file,coord_string,output_crd_folder)
        print("Done")
        

