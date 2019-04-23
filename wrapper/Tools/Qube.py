#
#
from Sire.IO import *
from Sire.Mol import *
from Sire.MM import *
from Sire.System import *
from Sire.Units import *
from Sire.CAS import * 
from Sire.Maths import * 
from Sire.Base import *


def readXmlParameters(pdbfile, xmlfile):
    # 1) Read a pdb file describing the system to simulate

    xmlfile = 'pyridine/MOL_extra.xml'
    pdbfile = 'pyridine/MOL.pdb'
    p = PDB2(pdbfile)
    s = p.toSystem()
    molecules = s.molecules()
    print (molecules)
    system = System()

     # 2) Now we read the xml file, and store parameters for each molecule


    import xml.dom.minidom as minidom
    xmldoc = minidom.parse(xmlfile)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~ TAG NAME: TYPE ~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    itemlist_type = xmldoc.getElementsByTagName('Type')
    dicts_type = []
    for items in itemlist_type:
        d = {}
        for a in items.attributes.values():
            d[a.name] = a.value
        dicts_type.append(d)
    dicts_tp =  str(dicts_type).split()
    print (dicts_tp)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~ TAG NAME: ATOM ~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    itemlist_atom = xmldoc.getElementsByTagName('Atom')
    dicts_atom = []
    for items in itemlist_atom:
        d = {}
        for a in items.attributes.values():
            d[a.name] = a.value
        dicts_atom.append(d)
    dicts_at =  str(dicts_atom).split()
    print (dicts_at)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~ TAG NAME: BOND ~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    itemlist_bond = xmldoc.getElementsByTagName('Bond')
    dicts_bond = []
    for items in itemlist_bond:
        d = {}
        for a in items.attributes.values():
            d[a.name] = a.value
        dicts_bond.append(d)
    dicts_b =  str(dicts_bond).split()
    print (dicts_b)

    nbond = itemlist_bond.length

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~ TAG NAME: ANGLE ~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    itemlist_angle = xmldoc.getElementsByTagName('Angle')
    dicts_angle = []
    for items in itemlist_angle:
        d = {}
        for a in items.attributes.values():
            d[a.name] = a.value
        dicts_angle.append(d)
    dicts_ang =  str(dicts_angle).split()
    print (dicts_ang)

    nAngles= itemlist_angle.length
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~ TAG NAME: PROPER ~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    itemlist_proper = xmldoc.getElementsByTagName('Proper')
    dicts_proper = []
    for items in itemlist_proper:
        d = {}
        for a in items.attributes.values():
            d[a.name] = a.value
        dicts_proper.append(d)
    dicts_pr =  str(dicts_proper).split()
    print (dicts_pr)

    nProper = itemlist_proper.length

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~ TAG NAME: IMPROPER ~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    itemlist_improper = xmldoc.getElementsByTagName('Improper')
    dicts_improper = []
    for items in itemlist_improper:
        d = {}
        for a in items.attributes.values():
            d[a.name] = a.value
        dicts_improper.append(d)
    dicts_impr =  str(dicts_improper).split()
    print (dicts_impr)
    nImproper = itemlist_improper.length

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~ TAG NAME: VIRTUAL SITES ~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    itemlist_VirtualSite = xmldoc.getElementsByTagName('VirtualSite')
    dicts_virtualsite = []
    for items in itemlist_VirtualSite:
        d = {}
        for a in items.attributes.values():
            d[a.name] = a.value
        dicts_virtualsite.append(d)
    dicts_vs =  str(dicts_virtualsite).split()
    #print (dicts_vs)
    nVirtualSites = itemlist_VirtualSite.length 
        

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~ TAG NAME: RESIDUE ~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    itemlist_residue = xmldoc.getElementsByTagName('Residue')
    dicts_residue = []
    for items in itemlist_residue:
        d = {}
        for a in items.attributes.values():
            d[a.name] = a.value
        dicts_residue.append(d)
    dicts_res =  str(dicts_residue).split()
    print (dicts_res)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~ TAG NAME: NON BONDED FORCE ~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    itemlist_nonbond = xmldoc.getElementsByTagName('NonbondedForce')
    dicts_nonb = []
    for items in itemlist_nonbond:
        d = {}
        for a in items.attributes.values():
            d[a.name] = a.value
        dicts_nonb.append(d)
    dicts_nb =  str(dicts_nonb).split()
    print (dicts_nb)
    nNonBonded = itemlist_nonbond.length

    # 3) Now we create an Amberparameters object for each molecule
    molnums = molecules.molNums()

    newmolecules = Molecules()
    for molnum in molnums:
        mol = molecules.at(molnum)
        print (mol) 

        # Add potential virtual site parameters
        if len(dicts_virtualsite) > 0:
            mol = mol.edit().setProperty("virtual-sites", vsiteListToProperty(dicts_virtualsite)).commit()

        # We populate the Amberparameters object with a list of bond, angle, dihedrals
        # We look up parameters from the contents of the xml file
        # We also have to set the atomic parameters (q, sigma, epsilon)

        editmol = mol.edit()
        mol_params = AmberParameters(editmol) #SireMol::AmberParameters()
        atoms = editmol.atoms()
        # We update atom parameters see setAtomParameters in SireIO/amber.cpp l2122 
        natoms = editmol.nAtoms()
        print("number of atoms is %s" %natoms)
        
    #natoms don't include the virtual sites! 

    # Loop over each molecule in the molecules object
    
        for atom in atoms: 
            editatom = editmol.atom(atom.index())
            i = int(str(atom.number()).split('(')[1].replace(")" , " "))
            editatom.setProperty("charge", float(dicts_atom[i+natoms+nVirtualSites-1]['charge']) * mod_electron)
            editatom.setProperty("mass", float(dicts_type[i-1]['mass']) * g_per_mol) #
            editatom.setProperty("LJ", LJParameter( float(dicts_atom[i+natoms+nVirtualSites-1]['sigma']) * angstrom, float(dicts_atom[i+natoms+nVirtualSites-1]['epsilon']) * kcal_per_mol))
            editatom.setProperty("ambertype", dicts_atom[i+natoms+nVirtualSites-1]['type'])
           
            editmol = editatom.molecule()
            print(editmol)
            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
         

        # Now we create a connectivity see setConnectivity in SireIO/amber.cpp l2144
        # XML data tells us how atoms are bonded in the molecule (Bond 'from' and 'to')
        
        if natoms > 1:
            print("Set up connectivity")

            c = []
            for i in range(0,int(nbond/2)):
               # i = int(str(atom.number()).split('(')[1].replace(")" , " "))
                if natoms > 1:   
                    connect_prop= {}
                    connect_prop = dicts_bond[i]['from'], dicts_bond[i]['to']
                c.append(connect_prop)
            print(c)

            conn = Connectivity(editmol.info()).edit()

            for j in range(0,len(c)):
                conn.connect(atoms[int(c[j][0]) ].index(), atoms[int(c[j][1]) ].index()) 
                   
            #for atom in atoms:
            editmol.setProperty("connectivity", conn.commit()).commit()
            #       print(editmol.setProperty("connectivity", conn.commit()))
            

             # Now we add bond parameters to the Sire molecule. We also update amberparameters see SireIO/amber.cpp l2154

            internalff = InternalFF()
                    
            bondfuncs = TwoAtomFunctions(mol)
            r = internalff.symbols().bond().r()

            for j in range(0,len(c)):
                bf = bondfuncs.set(atoms[int(c[j][0]) ].index(), atoms[int(c[j][1]) ].index(), float(dicts_bond[j+len(c)]['k'])* (float(dicts_bond[j+len(c)]['length']) - r) **2  )
                bond_id = BondID(atoms[int(c[j][0])].index(), atoms[int(c[j][1])].index())
                print(bond_id)
                for i in range(int(nbond/2), nbond):
                    mol_params.add(bond_id, float(dicts_bond[i]['k']), float(dicts_bond[i]['length']) ) 
            
            editmol.setProperty("bonds", bondfuncs).commit()

            
            mol_params.getAllBonds() 
            print(mol_params.getAllBonds() )


            editmol.setProperty("amberparameters", mol_params).commit() # Weird, should work - investigate ? 
            molecule = editmol.commit()
            #newmolecules.add(molecule)
        
    # return newmolecules

        # Now we add angle parameters to the Sire molecule. We also update amberparameters see SireIO/amber.cpp L2172
        if natoms > 2:
            print("Set up angles")

            anglefuncs = ThreeAtomFunctions(mol)

            at1 = []
            for i in range(0, nAngles):
                for j in range(0,natoms):
                    if dicts_angle[i]['class1']  == dicts_type[j]['class']:
                        a1 = {}
                        a1 = j
                        at1.append(a1)
            print (at1)

            at2 = []
            for i in range(0, nAngles):
                for j in range(0,natoms):
                    if dicts_angle[i]['class2']  == dicts_type[j]['class']:
                        a2 = {}
                        a2 = j
                        at2.append(a2)
            print (at2)

            at3 = []
            for i in range(0, nAngles):
                for j in range(0,natoms):
                    if dicts_angle[i]['class3']  == dicts_type[j]['class']:
                        a3 = {}
                        a3 = j
                        at3.append(a3)
            print (at3)

            theta = internalff.symbols().angle().theta()
            for j in range(0,nAngles):
                anglefuncs.set( atoms[at1[j]].index(), atoms[at2[j]].index(), atoms[at3[j]].index(), float(dicts_angle[j]['k']) * ( (float(dicts_angle[j]['angle']) *degrees.value() - theta )**2 ))
                angle_id = AngleID( atoms[int(at1[j])].index(), atoms[int(at2[j])].index(), atoms[int(at3[j])].index())
                print(angle_id)

                for i in range(0, nAngles):
                        mol_params.add(angle_id, float(dicts_angle[i]['k']), float(dicts_angle[i]['angle']) ) 
            print(mol_params.getAllAngles() )

        # Now we add dihedral parameters to the Sire molecule. We also update amberparameters see SireIO/amber.cpp L2190

        if natoms > 3:
            print("Set up dihedrals")

            di1 = []
            for i in range(0, nProper):
                for j in range(0,natoms):
                    if dicts_proper[i]['class1']  == dicts_type[j]['class']:
                        d1 = {}
                        d1 = j
                        di1.append(d1)
            print (di1)


            di2 = []
            for i in range(0, nProper):
                for j in range(0,natoms):
                    if dicts_proper[i]['class2']  == dicts_type[j]['class']:
                        d2 = {}
                        d2 = j
                        di2.append(d2)
            print (di2)


            di3 = []
            for i in range(0, nProper):
                for j in range(0,natoms):
                    if dicts_proper[i]['class3']  == dicts_type[j]['class']:
                        d3 = {}
                        d3 = j
                        di3.append(d3)
            print (di3)

            di4 = []
            for i in range(0, nProper):
                for j in range(0,natoms):
                    if dicts_proper[i]['class4']  == dicts_type[j]['class']:
                        d4 = {}
                        d4 = j
                        di4.append(d4)
            print (di4)

            dihedralfuncs = FourAtomFunctions(mol)
    
            phi = internalff.symbols().dihedral().phi()

            for i in range(1,5):    
                for j in range(0,nProper):
                    dihedralfuncs.set(atoms[di1[j]].index(), atoms[di2[j]].index(),atoms[di3[j]].index(),atoms[di4[j]].index(),\
                     float(dicts_proper[j]['k%s'%i])*(1+ Cos(int(dicts_proper[j]['periodicity%s'%i]) * phi - float(dicts_proper[j]['phase%s'%i]))))
                    dihedral_id = DihedralID( atoms[int(di1[j])].index(), atoms[int(di2[j])].index(), atoms[int(di3[j])].index(), atoms[int(di4[j])].index())
                    print(dihedral_id) 
                    for l in range(0, nProper ):
                        for t in range(1,5):
                            mol_params.add(dihedral_id, float(dicts_proper[l]['k%s'%t]), int(dicts_proper[l]['periodicity%s'%t]), float(dicts_proper[l]['phase%s'%t]) ) 
            print(mol_params.getAllDihedrals() )

            

            print("Set up impropers")

            di_im1 = []
            for i in range(0, nImproper):
                for j in range(0,natoms):
                    if dicts_improper[i]['class1']  == dicts_type[j]['class']:
                        d1 = {}
                        d1 = j
                        di_im1.append(d1)
            print (di_im1)


            di_im2 = []
            for i in range(0, nImproper):
                for j in range(0,natoms):
                    if dicts_improper[i]['class2']  == dicts_type[j]['class']:
                        d2 = {}
                        d2 = j
                        di_im2.append(d2)
            print (di_im2)


            di_im3 = []
            for i in range(0, nImproper):
                for j in range(0,natoms):
                    if dicts_improper[i]['class3']  == dicts_type[j]['class']:
                        d3 = {}
                        d3 = j
                        di_im3.append(d3)
            print (di_im3)


            di_im4 = []
            for i in range(0, nImproper):
                for j in range(0,natoms):
                    if dicts_improper[i]['class4']  == dicts_type[j]['class']:
                        d4 = {}
                        d4 = j
                        di_im4.append(d4)
            print (di_im4)


            improperfuncs = FourAtomFunctions(mol)

            phi_im = internalff.symbols().improper().phi()
            for i in range(1,5):    
                for j in range(0,nImproper):
                    improperfuncs.set(atoms[di_im1[j]].index(), atoms[di_im2[j]].index(),atoms[di_im3[j]].index(),atoms[di_im4[j]].index(),\
                     float(dicts_improper[j]['k%s'%i])*(1+ Cos(int(dicts_improper[j]['periodicity%s'%i]) * phi_im - float(dicts_improper[j]['phase%s'%i]))))
                    
                    improper_id = ImproperID( atoms[int(di_im1[j])].index(), atoms[int(di_im2[j])].index(), atoms[int(di_im3[j])].index(), atoms[int(di_im4[j])].index())
                    print(improper_id) 
                    for l in range(0, nImproper ):
                        for t in range(1,5):
                            mol_params.add(improper_id, float(dicts_improper[l]['k%s'%t]), int(dicts_improper[l]['periodicity%s'%t]), float(dicts_improper[l]['phase%s'%t]) ) 
            print(mol_params.getAllDihedrals() )

            mol = editmol.setProperty("bond", bondfuncs).commit()
            mol = editmol.setProperty("angle" , anglefuncs).commit()
            mol = editmol.setProperty("dihedral" , dihedralfuncs).commit()
            mol = editmol.setProperty("improper" , improperfuncs).commit()
            system.update(mol)

        # Now we work out non bonded pairs see SireIO/amber.cpp L2213


        print("Set up nbpairs")

        if natoms <= 3:
            import re
            scale_factor1 = 0
            scale_factor2 = 0
            nbpairs = CLJNBPairs(editmol.info(), CLJScaleFactor(scale_factor1,scale_factor2))

            for i in range(0, len(Connectivity(atoms).getBonds())):
                getB = str(Connectivity(atoms).getBonds()[i])
                nbpairs.set(atoms.index( int(re.findall('\d+',getB)[0])), atoms.index(int(re.findall('\d+',getB)[1])), CLJScaleFactor(scale_factor1,scale_factor2))
                mol = editmol.setProperty("intrascale" , nbpairs).commit()
            system.update(mol)

            for i in range(0, len(Connectivity(atoms).getAngles())):
                getA = str(Connectivity(atoms).getAngles()[i])
                nbpairs.set(atoms.index( int(re.findall('\d+',getA)[0])), atoms.index(int(re.findall('\d+',getA)[1])), CLJScaleFactor(scale_factor1,scale_factor2))
                mol = editmol.setProperty("intrascale" , nbpairs).commit()
            system.update(mol)


        elif natoms <= 4:
            import re
            scale_factor1 = float(dicts_nonb[i]['coulomb14scale'])
            scale_factor2 = float(dicts_nonb[i]['lj14scale'])
            nbpairs = CLJNBPairs(editmol.info(), CLJScaleFactor(scale_factor1,scale_factor2))

            for i in range(0, len(Connectivity(atoms).getDihedrals())):
                getD = str(Connectivity(atoms).getDihedrals()[i])
                nbpairs.set(atoms.index( int(re.findall('\d+',getD)[0])), atoms.index(int(re.findall('\d+',getD)[1])), CLJScaleFactor(scale_factor1,scale_factor2))
                mol = editmol.setProperty("intrascale" , nbpairs).commit()
            system.update(mol)


        else: 
            scale_factor1 = 1
            scale_factor2 = 1
            nbpairs = CLJNBPairs(editmol.info(), CLJScaleFactor(scale_factor1,scale_factor2))
            mol = editmol.setProperty("intrascale" , nbpairs).commit()
            system.update(mol)

        molecule = editmol.commit()
        newmolecules.add(molecule)

    return newmolecules
   
    # 4) Print the moleculeand atom properties that were loaded from the xml file
    # in a new file named xml_parameters.txt

    # new_file = open("xml_parameters.txt", "w")
    # new_file.write("\n")
    # new_file.write("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    # new_file.write("\n")
    # new_file.write("~~~~~~~~~~~~~~~~   Testing the xml parameters   ~~~~~~~~~~~~~~~")
    # new_file.write("\n")
    # new_file.write("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    # new_file.write("\n")
    # new_file.write("\n")
    # new_file.write("Connectivity of the molecule:\n")
    # new_file.write("{0}\n".format(conn))
    # new_file.write("\n")
    # new_file.write("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    # new_file.write("\n")
    # new_file.write("The bonds connect the following atoms:\n")
    # new_file.write("{0}\n".format(str(mol_params.getAllBonds()).replace(") ),", ") ),\n")))
    # new_file.write("\n")
    # new_file.write("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    # new_file.write("\n")
    # new_file.write("The angles are between the following atoms:\n")
    # new_file.write("{0}\n".format(str(mol_params.getAllAngles()).replace(") ),", ") ),\n")))
    # new_file.write("\n")
    # new_file.write("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    # new_file.write("\n")
    # new_file.write("The dihedrals are defined from the following atoms:\n")
    # new_file.write("{0}\n".format(str(mol_params.getAllDihedrals()).replace(") ),", ") ),\n")))
    # new_file.write("\n")
    # new_file.write("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    # new_file.write("\n")
    # new_file.write("The impropers are defined from the following atoms:\n")
    # new_file.write("{0}\n".format(str(mol_params.getAllImpropers()).replace(") ),", ") ),\n")))
    # new_file.write("\n")
    # new_file.write("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    # new_file.write("\n")
    # new_file.write("The overal properties of the new molecule are the following: ")
    # new_file.write("{0}\n".format(str(editmol.properties())))
    # new_file.close()

        # By the end of this loop we have a new set of mol that looks
        # exactly like a molecules object returned by Amber().readCrdTop(...) 

if __name__ == '__main__':
    molecules = readXmlParameters('pyridine/MOL.pdb', 'pyridine/MOL_extra.xml')

    m0 = molecules.first().molecule()
    print (m0.properties())
    editmol = m0.edit()