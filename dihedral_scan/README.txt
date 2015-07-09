Documentation         9th July 2015
-----------------------------------------

GENERAL

-  all_dihedral_scan.py  is a Sire python script which allows to scan all the dihedrals by 10 degrees interval from 0 to 350 degrees (hard coded for the moment). 
   Output files  are: 
   amber_inps/alanine_index_degrees.mdcrd   all coordinates files which need minimization
   ALANINE.mdcrd 			    one file with all coordinates written in a "vmd" comfortable format
   atom_list.dat			    list of all the dihedral modified, with atom present and atom number. Necessary for dealing with parmed.py
  
   To check your final structure in vmd type in the terminal:   vmd a.prmtop ALANINE.mdcrd 

-  OpenMM_minimizer.py  python script which use parmed.py (import chemistry) and OpenMM python modules. First harmonic restrains are imposed on atoms enlisted in atom_list.dat, secondly 
   a minimization is performed, finally pdb and crd files are saved.
   Output files are:
   amber_out_pdb/			    folder with all pdb files in order to visualize in vmd (in the future this folder can be skipped)
   amber_out_crd/		            folder with all coordinate files. All these files will be used for paramfit.

-  alanine_test_files/ contains a.prmtop and a.inpcrd, topology and coordinate files of a tri-alanine peptide. It can be use for a tutorial/test
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


USAGE


-  in ~/sire.app/bin/python run all_dihedral_scan.py:      ~/sire.app/bin/python 
                                                           %run all_dihedral_scan.py  a.prmtop a.inpcrd
   

-  After all_dihedral_scan.py finished to work simply using normal python:     python OpenMM_minimizer.py
   it will take a while, because it has to perform around 2000 minimizations. 


----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



Paramfit note:


At the moment the idea is to calculate the molecular mechanic energies using paramfit of all the structures present in amber_out_crd and take the lowest energy conformations. 
After that works are still in progress to understand how many dihedrals can be fitted per time, what are the best choices and so on.
      
