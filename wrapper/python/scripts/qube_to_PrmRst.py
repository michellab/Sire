
import argparse
import os
import sys
from Sire.Tools.OpenMMMD import *

# Reads qubekit output files (pdb/xml) and converts them to amber format (prm7/rst7)
# usage: ~/sire.app/bin/qube_to_PrmRst -p PDB_FILE -x XML_FILE

def vsiteListToProperty(list):
    prop = Properties()
    i = 0
    for entry in list:
        for key, value in entry.items():
            prop.setProperty("%s(%d)" % (key,i), VariantProperty(value))
        i += 1
    prop.setProperty("nvirtualsites",VariantProperty(i))
    return prop



def readXmlParameters(pdbfile, xmlfile):
# 1) Read a pdb file describing the system to simulate

    p = PDB2(pdbfile)
    s = p.toSystem()
    molecules = s.molecules()
    #print (molecules)
    with open (pdbfile, "r") as f:
        for line in f:
            if line.split()[0] == "CRYST1" :
                print (line)
                pbc_x = float(line.split()[1])
                pbc_y = float(line.split()[2])
                pbc_z = float(line.split()[3])
                space = PeriodicBox(Vector(pbc_x, pbc_y, pbc_z))
                break
            else:
                space = Cartesian()
    #print("space:", space)

    system = System()

     # 2) Now we read the xml file, and store parameters for each molecule


    import xml.dom.minidom as minidom
    xmldoc = minidom.parse(xmlfile)

    itemlist_type = xmldoc.getElementsByTagName('Type')
    dicts_type = []
    for items in itemlist_type:
        d = {}
        for a in items.attributes.values():
            d[a.name] = a.value
        dicts_type.append(d)
    dicts_tp =  str(dicts_type).split()

    itemlist_atom = xmldoc.getElementsByTagName('Atom')
    dicts_atom = []
    for items in itemlist_atom:
        d = {}
        for a in items.attributes.values():
            d[a.name] = a.value
        dicts_atom.append(d)
    dicts_at =  str(dicts_atom).split()

    itemlist_bond = xmldoc.getElementsByTagName('Bond')
    dicts_bond = []
    for items in itemlist_bond:
        d = {}
        for a in items.attributes.values():
            d[a.name] = a.value
        dicts_bond.append(d)
    dicts_b =  str(dicts_bond).split()
    #print (dicts_b)

    nbond = itemlist_bond.length

    itemlist_angle = xmldoc.getElementsByTagName('Angle')
    dicts_angle = []
    for items in itemlist_angle:
        d = {}
        for a in items.attributes.values():
            d[a.name] = a.value
        dicts_angle.append(d)
    dicts_ang =  str(dicts_angle).split()
    #print (dicts_angle)

    nAngles= itemlist_angle.length

    itemlist_proper = xmldoc.getElementsByTagName('Proper')
    dicts_proper = []
    for items in itemlist_proper:
        d = {}
        for a in items.attributes.values():
            d[a.name] = a.value
        dicts_proper.append(d)
    dicts_pr =  str(dicts_proper).split()
    nProper = itemlist_proper.length


    itemlist_improper = xmldoc.getElementsByTagName('Improper')
    dicts_improper = []
    for items in itemlist_improper:
        d = {}
        for a in items.attributes.values():
            d[a.name] = a.value
        dicts_improper.append(d)
    dicts_impr =  str(dicts_improper).split()
    #print (dicts_impr)
    nImproper = itemlist_improper.length

    itemlist_VirtualSite = xmldoc.getElementsByTagName('VirtualSite')
    dicts_virtualsite = []
    for items in itemlist_VirtualSite:
        d = {}
        for a in items.attributes.values():
            d[a.name] = a.value
        dicts_virtualsite.append(d)

    nVirtualSites = itemlist_VirtualSite.length 

    v_site_CLJ = []
    for i in range(0, int(len(dicts_atom))):
        if dicts_atom[i]['type'][0] == 'v':
            v_site_CLJ = dicts_atom[i]
            dicts_virtualsite.append(v_site_CLJ)

    for i in range(0, len(itemlist_VirtualSite)):
        dicts_virtualsite[i].update(dicts_virtualsite[i+len(itemlist_VirtualSite)])
        dicts_virtualsite[i].update(dicts_virtualsite[i+2*len(itemlist_VirtualSite)]) 

    dict_vs = []
    for i in range(0, len(itemlist_VirtualSite)):
        dicts_virtualsite[i]
        dict_vs.append(dicts_virtualsite[i])

    itemlist_residue = xmldoc.getElementsByTagName('Residue')
    dicts_residue = []
    for items in itemlist_residue:
        d = {}
        for a in items.attributes.values():
            d[a.name] = a.value
        dicts_residue.append(d)
    dicts_res =  str(dicts_residue).split()

    itemlist_nonbond = xmldoc.getElementsByTagName('NonbondedForce')
    dicts_nonb = []
    for items in itemlist_nonbond:
        d = {}
        for a in items.attributes.values():
            d[a.name] = a.value
        dicts_nonb.append(d)
    dicts_nb =  str(dicts_nonb).split()
    #print (dicts_nb)
    nNonBonded = itemlist_nonbond.length

    # 3) Now we create an Amberparameters object for each molecule
    molnums = molecules.molNums()

    newmolecules = Molecules()
    for molnum in molnums:
        mol = molecules.at(molnum)

        editmol = mol.edit()
        mol_params = AmberParameters(editmol) #SireMol::AmberParameters()
        atoms = editmol.atoms()
        # We update atom parameters see setAtomParameters in SireIO/amber.cpp l2122 
        natoms = editmol.nAtoms()
        
        opls=[]
        for i in range (0, int(len(dicts_atom)/2)): 
            opl={} 
            opl = dicts_atom[i]['type'] 
            opls.append(opl) 

        name=[]
        for i in range (0, int(len(dicts_atom)/2)): 
            nm={} 
            nm = dicts_atom[i]['name'] 
            name.append(nm) 

        two=[] 
        #print(len(name)) 
        for i in range(0, len(name)): 
            t=(opls[i],name[i]) 
            two.append(t) 

        import numpy as np 

        atom_sorted = []
        for j in range(0, len(two)): 
            for i in range(int(len(dicts_atom)/2), len(dicts_atom)):   
                if dicts_atom[i]['type'] == two[j][0]: 
                    dic_a = {}
                    dic_a = dicts_atom[i]
                    atom_sorted.append(dic_a)
      
        type_sorted = []
        for j in range(0, len(two)): 
            for i in range(0, int(len(dicts_type))):   
                if dicts_type[i]['name'] == two[j][0]: 
                    dic_t = {}
                    dic_t = dicts_type[i]
                    type_sorted.append(dic_t)
        print(" ")
        print("There are ",natoms," atoms in this molecule. ")
        print("*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*")

        for atom in atoms: 
            editatom = editmol.atom(atom.index())

            i = int(str(atom.number()).split('(')[1].replace(")" , " ")) 

            editatom.setProperty("charge", float(atom_sorted[i-1]['charge']) * mod_electron)
            editatom.setProperty("mass", float(type_sorted[i-1]['mass']) * g_per_mol) 
            editatom.setProperty("LJ", LJParameter( float(atom_sorted[i-1]['sigma'])*10 * angstrom , float(atom_sorted[i-1]['epsilon'])/4.184 * kcal_per_mol))
            editatom.setProperty("ambertype", dicts_atom[i-1]['type'])
           
            editmol = editatom.molecule()

        # Now we create a connectivity see setConnectivity in SireIO/amber.cpp l2144
        # XML data tells us how atoms are bonded in the molecule (Bond 'from' and 'to')
        
        if natoms > 1:
            print("Set up connectivity")

            con = []
            for i in range(0,int(nbond/2)):
                if natoms > 1:   
                    connect_prop= {}
                    connect_prop = dicts_bond[i]['from'], dicts_bond[i]['to']
                con.append(connect_prop)

            conn = Connectivity(editmol.info()).edit()

            for j in range(0,len(con)):
                conn.connect(atoms[int(con[j][0]) ].index(), atoms[int(con[j][1]) ].index()) 
            
            editmol.setProperty("connectivity", conn.commit()).commit()
            mol = editmol.setProperty("connectivity", conn.commit()).commit()
            system.update(mol)

             # Now we add bond parameters to the Sire molecule. We also update amberparameters see SireIO/amber.cpp l2154

            internalff = InternalFF()
                    
            bondfuncs = TwoAtomFunctions(mol)
            r = internalff.symbols().bond().r()

            for j in range(0,len(con)):
                bondfuncs.set(atoms[int(con[j][0]) ].index(), atoms[int(con[j][1]) ].index(), float(dicts_bond[j+len(con)]['k'])/(2*100*4.184)* (float(dicts_bond[j+len(con)]['length'])*10 - r) **2  )
                bond_id = BondID(atoms[int(con[j][0])].index(), atoms[int(con[j][1])].index())
                mol_params.add(bond_id, float(dicts_bond[j+len(con)]['k'])/(2*100*4.184), float(dicts_bond[j+len(con)]['length'])*10 )   
                editmol.setProperty("bonds", bondfuncs).commit()
                molecule = editmol.commit()
             
            mol_params.getAllBonds() 

            editmol.setProperty("amberparameters", mol_params).commit() # Weird, should work - investigate ? 
            molecule = editmol.commit()

        # Now we add angle parameters to the Sire molecule. We also update amberparameters see SireIO/amber.cpp L2172
        if natoms > 2:
            print("Set up angles")

            anglefuncs = ThreeAtomFunctions(mol)
            at1 = []
            for i in range(0, nAngles):
                a1 = {}
                to_str1 = str(re.findall(r"\d+",str(dicts_angle[i]['class1'])))
                if dicts_atom[i]['type'][0] == 'o': #if opls_
                    a1 = int(to_str1.replace("[","").replace("]","").replace("'","") )-800
                else:#if QUBE_
                    a1 = int(to_str1.replace("[","").replace("]","").replace("'","") )
                at1.append(a1)

            at2 = []
            for i in range(0, nAngles):
                a2 = {}
                to_str2 = str(re.findall(r"\d+",str(dicts_angle[i]['class2'])))
                if dicts_atom[i]['type'][0] == 'o': #if opls_
                    a2 = int(to_str2.replace("[","").replace("]","").replace("'","") )-800
                else: #if QUBE_
                    a2 = int(to_str2.replace("[","").replace("]","").replace("'","") )
                at2.append(a2)

            at3 = []
            for i in range(0, nAngles):
                a3 = {}
                to_str3 = str(re.findall(r"\d+",str(dicts_angle[i]['class3'])))
                if dicts_atom[i]['type'][0] == 'o': #if opls_
                    a3 = int(to_str3.replace("[","").replace("]","").replace("'","") )-800
                else: #if QUBE_
                    a3 = int(to_str3.replace("[","").replace("]","").replace("'","") )
                at3.append(a3)

            theta = internalff.symbols().angle().theta()
            for j in range(0,nAngles):
                anglefuncs.set( atoms[at1[j]].index(), atoms[at2[j]].index(), atoms[at3[j]].index(), float(dicts_angle[j]['k'])/(2*4.184) * ( (float(dicts_angle[j]['angle']) - theta )**2 ))
                angle_id = AngleID( atoms[int(at1[j])].index(), atoms[int(at2[j])].index(), atoms[int(at3[j])].index())
                mol_params.add(angle_id, float(dicts_angle[j]['k'])/(2*4.184), float(dicts_angle[j]['angle']) ) 
            
        # Now we add dihedral parameters to the Sire molecule. We also update amberparameters see SireIO/amber.cpp L2190

        if natoms > 3:
            print("Set up dihedrals")
            di1 = []

            for i in range(0, nProper):
                d1 = {}
                to_str1 = str(re.findall(r"\d+",str(dicts_proper[i]['class1'])))
                if dicts_atom[0]['type'][0] == 'o':#if opls_
                    d1 = int(to_str1.replace("[","").replace("]","").replace("'","") )-800
                else:  #if QUBE_
                    d1 = int(to_str1.replace("[","").replace("]","").replace("'","") )
                di1.append(d1)


            di2 = []
            for i in range(0, nProper):
                d2 = {}
                to_str2 = str(re.findall(r"\d+",str(dicts_proper[i]['class2'])))
                if dicts_atom[0]['type'][0] == 'o':#if opls_
                    d2 = int(to_str2.replace("[","").replace("]","").replace("'","") )-800
                else: #if QUBE_
                    d2 = int(to_str2.replace("[","").replace("]","").replace("'","") )
                di2.append(d2)

            di3 = []
            for i in range(0, nProper):
                d3 = {}
                to_str3 = str(re.findall(r"\d+",str(dicts_proper[i]['class3'])))
                if dicts_atom[0]['type'][0] == 'o':#if opls_
                    d3 = int(to_str3.replace("[","").replace("]","").replace("'","") )-800
                else: #if QUBE_
                    d3 = int(to_str3.replace("[","").replace("]","").replace("'","") )
                di3.append(d3)
          

            di4 = []
            for i in range(0, nProper):
                d4 = {}
                to_str4 = str(re.findall(r"\d+",str(dicts_proper[i]['class4'])))
                if dicts_atom[0]['type'][0] == 'o':#if opls_
                    d4 = int(to_str4.replace("[","").replace("]","").replace("'","") )-800
                else: #if QUBE_
                    d4 = int(to_str4.replace("[","").replace("]","").replace("'","") )
                di4.append(d4)

            dihedralfuncs = FourAtomFunctions(mol)
    
            phi = internalff.symbols().dihedral().phi()
            for i in range(0,nProper):  
                if atoms[int(di1[i])].index() != atoms[int(di4[i])].index():
                    dihedral_id = DihedralID( atoms[int(di1[i])].index(), atoms[int(di2[i])].index(), atoms[int(di3[i])].index(), atoms[int(di4[i])].index()) 
                    dih1= float(dicts_proper[i]['k1'])/4.184*(1+Cos(int(dicts_proper[i]['periodicity1'])* phi- float(dicts_proper[i]['phase1'])))
                    dih2= float(dicts_proper[i]['k2'])/4.184*(1+Cos(int(dicts_proper[i]['periodicity2'])* phi- float(dicts_proper[i]['phase2'])))
                    dih3= float(dicts_proper[i]['k3'])/4.184*(1+Cos(int(dicts_proper[i]['periodicity3'])* phi- float(dicts_proper[i]['phase3'])))
                    dih4= float(dicts_proper[i]['k4'])/4.184*(1+Cos(int(dicts_proper[i]['periodicity4'])* phi- float(dicts_proper[i]['phase4'])))
                    dih_fun = dih1 + dih2 +dih3 +dih4
                    dihedralfuncs.set(dihedral_id, dih_fun)

                    for t in range(1,5):
                        mol_params.add(dihedral_id, float(dicts_proper[i]['k%s'%t])/4.184, int(dicts_proper[i]['periodicity%s'%t]), float(dicts_proper[i]['phase%s'%t]) ) 

            print("Set up impropers")

            di_im1 = []
            for i in range(0, nImproper):
                d1 = {}
                to_str1 = str(re.findall(r"\d+",str(dicts_improper[i]['class1'])))
                if dicts_atom[0]['type'][0] == 'o':#if opls_
                    d1 = int(to_str1.replace("[","").replace("]","").replace("'","") )-800
                else:
                    d1 = int(to_str1.replace("[","").replace("]","").replace("'","") )
                di_im1.append(d1)

            di_im2 = []
            for i in range(0, nImproper):
                d2 = {}
                to_str2 = str(re.findall(r"\d+",str(dicts_improper[i]['class2'])))
                if dicts_atom[0]['type'][0] == 'o':#if opls_
                    d2 = int(to_str2.replace("[","").replace("]","").replace("'","") )-800
                else:
                    d2 = int(to_str2.replace("[","").replace("]","").replace("'","") )
                di_im2.append(d2)

            di_im3 = []
            for i in range(0, nImproper):
                d3 = {}
                to_str3 = str(re.findall(r"\d+",str(dicts_improper[i]['class3'])))
                if dicts_atom[0]['type'][0] == 'o':#if opls_
                    d3 = int(to_str3.replace("[","").replace("]","").replace("'","") )-800
                else:
                    d3 = int(to_str3.replace("[","").replace("]","").replace("'","") )
                di_im3.append(d3)

            di_im4 = []
            for i in range(0, nImproper):
                d4 = {}
                to_str4 = str(re.findall(r"\d+",str(dicts_improper[i]['class4'])))
                if dicts_atom[0]['type'][0] == 'o':#if opls_
                    d4 = int(to_str4.replace("[","").replace("]","").replace("'","") )-800
                else:
                    d4 = int(to_str4.replace("[","").replace("]","").replace("'","") )
                di_im4.append(d4)
            
            improperfuncs = FourAtomFunctions(mol)
            phi_im = internalff.symbols().improper().phi()

            for i in range(0,nImproper):  
                improper_id = ImproperID( atoms[int(di_im2[i])].index(), atoms[int(di_im3[i])].index(), atoms[int(di_im1[i])].index(), atoms[int(di_im4[i])].index()) 
                imp1= float(dicts_improper[i]['k1'])*(1/4.184)*(1+Cos(int(dicts_improper[i]['periodicity1'])* phi_im - float(dicts_improper[i]['phase1'])))
                imp2= float(dicts_improper[i]['k2'])*(1/4.184)*(1+Cos(int(dicts_improper[i]['periodicity2'])* phi_im - float(dicts_improper[i]['phase2'])))
                imp3= float(dicts_improper[i]['k3'])*(1/4.184)*(1+Cos(int(dicts_improper[i]['periodicity3'])* phi_im - float(dicts_improper[i]['phase3'])))
                imp4= float(dicts_improper[i]['k4'])*(1/4.184)*(1+Cos(int(dicts_improper[i]['periodicity4'])* phi_im - float(dicts_improper[i]['phase4'])))
                imp_fun = imp1 + imp2 +imp3 +imp4
                improperfuncs.set(improper_id, imp_fun)
                #print(improperfuncs.potentials())

                for t in range(1,5):
                    mol_params.add(improper_id, float(dicts_improper[i]['k%s'%t])*(1/4.184), int(dicts_improper[i]['periodicity%s'%t]), float(dicts_improper[i]['phase%s'%t]) ) 


            mol = editmol.setProperty("bond", bondfuncs).commit()
            mol = editmol.setProperty("angle" , anglefuncs).commit()
            mol = editmol.setProperty("dihedral" , dihedralfuncs).commit()
            mol = editmol.setProperty("improper" , improperfuncs).commit()
            system.update(mol)

        # Now we work out non bonded pairs see SireIO/amber.cpp L2213


            print("Set up non-bonded pairs")
            print("*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*")
            ## Define the bonded pairs in a list that is called are12
            #print("Now calculating 1-2 intercactions")
            are12 = []
            for i in range(0, natoms): 
                for j in range (0, natoms): 
                    if conn.areBonded(atoms[i].index(), atoms[j].index()) == True:
                        #ij = {}
                        ij= (i,j)
                        are12.append(ij)
            are12_bckup = are12[:]


            #print("Now calculating 1-3 intercactions")
            are13 = []
            for i in range(0, natoms): 
                for j in range (0, natoms): 
                    if conn.areAngled(atoms[i].index(), atoms[j].index()) == True:
                        ij = {}
                        ij= (i,j)
                        are13.append(ij)
            are13_bckup = are13[:]

           # print("Now calculating 1-4 intercactions")
            are14 = []
            for i in range(0, natoms): 
                for j in range (0, natoms):
                   
                    if conn.areDihedraled(atoms[i].index(), atoms[j].index()) == True and conn.areAngled(atoms[i].index(), atoms[j].index()) == False:
                        ij = {}
                        ij= (i,j)
                        are14.append(ij)
            are14_bckup = are14[:]

           # print("Now calculating the non-bonded intercactions")
            bonded_pairs_list = are12_bckup + are13_bckup + are14_bckup    
            nb_pair_list =[]

            for i in range(0, natoms): 
                #print("i=",i)
                for j in range (0, natoms):
                    if i != j and (i,j) not in bonded_pairs_list:
                        nb_pair_list.append((i,j))
            are_nb_bckup = nb_pair_list[:]

            nbpairs = CLJNBPairs(editmol.info(), CLJScaleFactor(0,0))

            for i in range(0, len(are12)):
                scale_factor1 = 0
                scale_factor2 = 0
                nbpairs.set(atoms.index( int(are12[i][0])), atoms.index(int(are12[i][1])), CLJScaleFactor(scale_factor1,scale_factor2))

            for i in range(0, len(are13)):
                scale_factor1 = 0
                scale_factor2 = 0
                nbpairs.set(atoms.index( int(are13[i][0])), atoms.index(int(are13[i][1])), CLJScaleFactor(scale_factor1,scale_factor2))
            
            for i in range(0, len(are14)): 
                scale_factor1 = 1/2
                scale_factor2 = 1/2
                nbpairs.set(atoms.index( int(are14[i][0])), atoms.index(int(are14[i][1])), CLJScaleFactor(scale_factor1,scale_factor2))
                mol_params.add14Pair(BondID(atoms.index( int(are14[i][0])), atoms.index( int(are14[i][1]))),scale_factor1 , scale_factor2)
           
            for i in range(0, len(nb_pair_list)):  
                scale_factor1 = 1
                scale_factor2 = 1
                nbpairs.set(atoms.index( int(nb_pair_list[i][0])), atoms.index(int(nb_pair_list[i][1])), CLJScaleFactor(scale_factor1,scale_factor2))
                
            mol = editmol.setProperty("intrascale" , nbpairs).commit()
            system.update(mol)

        molecule = editmol.commit()
        newmolecules.add(molecule)

    return (newmolecules, space)

parser = argparse.ArgumentParser(description="Performs file conversion from pdb and xml to rst7, prm7", prog="qube_to_PrmRst")
parser.add_argument('-p', '--pdbfile', nargs="?", help="The pdb file with the coordinates of the molecule.")
parser.add_argument('-x', '--xmlfile', nargs="?", help="The xml file with the parameters of the molecule.")
args = parser.parse_args()
(lig1, space) = readXmlParameters(args.pdbfile, args.xmlfile)
system = createSystem(lig1) # Sire.System._System.System
system = setupForcefields(system, space)
rst = Sire.IO.AmberRst7(system)
prm = AmberPrm(system)
print("Writing the prm7 and rst7 files")
pdb_name = str(args.pdbfile).split('.')[0]
prm.writeToFile("%s.prm7"%pdb_name)
rst.writeToFile("%s.rst7"%pdb_name)
print("Process completed!")
print()
print("Generated files: %s.prm7 and %s.rst7"%(pdb_name, pdb_name))
print()
