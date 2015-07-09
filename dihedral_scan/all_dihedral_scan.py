import os,re, sys, shutil
import math
import numpy


from Sire.IO import *
from Sire.Mol import *
from Sire.CAS import *
from Sire.System import *
from Sire.Move import *
from Sire.MM import *
from Sire.FF import *
from Sire.Units import *
from Sire.Vol import *
from Sire.Maths import *
from Sire.Base import *
from Sire.Qt import *
from Sire.ID import *
from Sire.Config import *
from Sire.Analysis import *

from Sire.Tools.DCDFile import *

from Sire.Tools import Parameter, resolveParameters

import Sire.Stream



def createSystem(molecules):
    
    print("Creating the system...")

    moleculeNumbers = molecules.molNums()
    moleculeList = []

    for moleculeNumber in moleculeNumbers:
        molecule = molecules.molecule(moleculeNumber).molecule()
        moleculeList.append(molecule)

    molecules = MoleculeGroup("molecules")
    ions = MoleculeGroup("ions")
    
    for molecule in moleculeList:
        natoms = molecule.nAtoms()
        if natoms == 1:
            ions.add(molecule)
        else:
            molecules.add(molecule)

    all = MoleculeGroup("all")
    all.add(molecules)
    all.add(ions)

    system = System()

    system.add(all)
    system.add(molecules)
    system.add(ions)
       
    return system

def atomList(dihedral_atom_list,atom_dictionary):

  

    atom_list_file = open("atom_list.dat","w")
    
    index=0
    
    for dihedral in dihedral_atom_list:
        vox = atom_dictionary[str(dihedral)]
        
        for v in vox:
            atom_list_file.write("%s   " %str(v))
            index+=1
            if index==8:
                atom_list_file.write("\n")
                index=0

        

def generateDihedrals(old_solute,bigfile,atom_dictionary): 

    connectivity = old_solute.property("connectivity")
    all_dihedrals = connectivity.getDihedrals()
    
    dihedral_atom_list = []
    
    
    dihedral_len = len(all_dihedrals)
    index_order = 0
    for dihedral in all_dihedrals:
        atm0= dihedral.atom0()
        atm1= dihedral.atom1()
        atm2= dihedral.atom2()
        atm3= dihedral.atom3()
        dihedral_atom_list.append(atm0)
        dihedral_atom_list.append(atm1)
        dihedral_atom_list.append(atm2)
        dihedral_atom_list.append(atm3)
        dihbond = BondID(dihedral.atom1(), dihedral.atom2())
        for deg in range(0,360,10):     #scan every 10 degreees
            new_solute = old_solute.move().change(dihbond,deg*degrees).commit()
           
            generateCoordFiles(bigfile,new_solute.property("coordinates"),deg,index_order)
        index_order+=1
                
    atomList(dihedral_atom_list, atom_dictionary)      
          
def generateCoordFiles(bigfile,coordinates,degrees,idx):
     
     
     coord_toVect = coordinates.toVector()
     string_coord = str(coord_toVect)
     string_rev   = re.sub("[^0-9.e-]", " ", string_coord)
     string_split = string_rev.split()
     
     big_file_index = 0
     for f in string_split:
         if float(f)<0:
             bigfile.write("  %.3f" % float(f))
         elif float(f)>=10:
             bigfile.write("  %.3f" % float(f))
         else:
             bigfile.write("   %.3f" % float(f))
         big_file_index+=1
         if not big_file_index%10:
             bigfile.write("\n")
     bigfile.write("\n")

     small_file_index=0  
     
     file_name = "amber_inps/alanine" + "_" + str(idx) + "_" + str(degrees) + ".mdcrd"
     singlecoord = open(file_name, "w")
     singlecoord.write("LIG\n")
     singlecoord.write("    32\n")
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


 
########################################MAIN#################################################################

if __name__ == "__main__":
 

     try:
         top_file = sys.argv[1]
         crd_file = sys.argv[2]
        
        
     except IndexError:
        
         top_file = "SYSTEM.top"
         crd_file = "SYSTEM.crd"

     amber_inps= "amber_inps"
     if not os.path.exists(amber_inps):
         os.makedirs(amber_inps)

     amber = Amber()
     
     molecules, space = amber.readCrdTop(crd_file, top_file)
     system = createSystem(molecules)
     old_solute = system[MGName("all")].moleculeAt(0).molecule()

   
     bigfile = open("ALANINE.mdcrd", "a")
     bigfile.write("Cpptraj generator\n")
   
     atom_dictionary = {}
     atoms = old_solute.atoms()
     for i in range(0, atoms.count()):
         at = atoms[i]
         at_idx = at.index()
         at_name = at.name().value()
         at_number = int(at.number().value()) -1  #-1 due to mismatcing between chemistry.amber and sire
         atom_dictionary[str(at_idx)] = [ at_name, at_number]

     generateDihedrals(old_solute,bigfile,atom_dictionary) 

   
