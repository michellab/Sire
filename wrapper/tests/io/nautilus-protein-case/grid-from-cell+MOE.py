#!/usr/bin/python

import os,sys, struct

from math import exp, log

# The input folder where the cell files are stored
cell_dir = "./cell_oldcode"
# The ouput folder where the grid files will be saved
grid_dir = "./grid-oldcode-0.5"
# The grid will be averaged using all the cell files in the time interval below. Set to -1 and 1e10 if you want to use all data. Note the interval is in frame numbers
cell_interval = (0, 22)
# The grid.in file defines the grid location
grid_infile = "./grid.in"
# The grid step defines the grid density in x/y/z
grid_step = 0.5 # Angstroms
# Grid points with less than count will be discarded
grid_count_cutoff = 0
# temperature
temperature = 298 # kelvin
# averageupdate frequency
frequencyupdate = 10
# Whether to overwrite existing output
Overwrite = True
# The water model used for the simulations
waterModel = "TIP4PEW-SireOpenMM"
#########################################
T = temperature
SMALL = 1. #0.00000001
MININF = -99999999.0
NAVOGADRO = 6.022141e23
MASS_WATER = 18.01528 # g/mol
KCAL_TO_J = 4184
KB = 0.001987206500956023# * 1000. * 4.184
## DICTIONNARY OF BULK PROPERTIES PER WATER MODEL ###############
# For consistency, all should be derived using same cutoff/code.# 
#################################################################
# TIP3P-SireOpenMM: 10 ns NPT, 10 Angstrom cutoff, reaction field, 10k snapshots
# RF seems to lower density and make 
# TIP4P-RH From J. Chem. Phys. 126, 064504 2007 
# Forces are in N x 10e-10
# Torques are in N m-1 x 10e-20
# Energies are in kcal/mol
# Densities are in kg/m3

waterProps = {"TIP4P-RH":            { "BULKORI": 2.94,
                                       "BULKFORCES": [ 1.627, 1.518, 1.284],
                                       "BULKTORQUES": [ 1.157, 0.994, 1.371],
                                       "BULKNRG": -9.92, 
                                       "BULKDENSITY": "X.XX"
                                     },
               "TIP4PEW-SireOpenMM": { "BULKORI": 3.305,
                                       "BULKFORCES": [ 1.587, 1.735, 1.334 ],
                                       "BULKTORQUES": [ 1.061, 1.194, 1.453 ],
                                       "BULKNRG": -11.025,
                                       "BULKDENSITY":  995.4
                                     }
               }
###############################################

BULKORI = waterProps[waterModel]['BULKORI']
BULKFORCES = waterProps[waterModel]['BULKFORCES']
BULKTORQUES = waterProps[waterModel]['BULKTORQUES']
BULKNRG = waterProps[waterModel]['BULKNRG']
BULKDENSITY = waterProps[waterModel]['BULKDENSITY']
BULKWATMOLPERVOXEL = (BULKDENSITY * (grid_step * 1e-9) ** 3 / MASS_WATER) * NAVOGADRO


def initGrid( grid_infile, grid_step):
    
    stream = open(grid_infile,'r')
    lines = stream.readlines()
    stream.close()

    #print lines
    for x in lines:
        #print x 
        if x.startswith("#"):
            continue
        elems = x.split()
#center x/y/z then +DX -DX +DY -DY +DZ -DZ

        try:
            centerx = float(elems[0])
            centery = float(elems[1])
            centerz = float(elems[2])
            maxx = float(elems[3])
            minx = float(elems[4])
            maxy = float(elems[5])
            miny = float(elems[6])
            maxz = float(elems[7])
            minz = float(elems[8])
        except:
            print "Problem with format of %s. Please check." % grid_infile
            sys.exit(-1)

        minx = centerx - minx
        miny = centery - miny 
        minz = centerz - minz 
    
        vec_min = [ minx, miny, minz ]

        maxx = centerx + maxx 
        maxy = centery + maxy 
        maxz = centerz + maxz 

        vec_max = [ maxx, maxy, maxz ]

    print "# The grid coordinates range from (%8.3f, %8.3f, %8.3f) to (%8.3f, %8.3f, %8.3f) " % ( minx, miny, minz, maxx, maxy, maxz)

    nx = int ( round ( ( vec_max[0] - vec_min[0] ) / grid_step ) ) 
    ny = int ( round ( ( vec_max[1] - vec_min[1] ) / grid_step ) ) 
    nz = int ( round ( ( vec_max[2] - vec_min[2] ) / grid_step ) ) 

    origin_x = minx + grid_step / 2.0
    origin_y = miny + grid_step / 2.0
    origin_z = minz + grid_step / 2.0

    grid = {}
    grid['origin'] = [ origin_x, origin_y, origin_z ] 
    grid['dimensions'] = ( nx, ny, nz )
    grid['step'] = grid_step
    grid['points'] = {}
    grid['npoints'] = 0
    grid['nsnapshots'] = 0

    grid_counter = 0

    for x in xrange(0,nx):
        grid['points'][x] = {}
        for y in xrange(0,ny):
            grid['points'][x][y] = {}
            for z in xrange(0,nz):
                grid['points'][x][y][z] = {}
                
                grid['points'][x][y][z]['avgnrg'] = [ 0.0, 0.0, 0.0 ] 
                grid['points'][x][y][z]['avgsolnrg'] = [ 0.0, 0.0, 0.0 ]
                grid['points'][x][y][z]['avgforces'] = [ 0.0, 0.0, 0.0 ] 
                grid['points'][x][y][z]['avgtorques'] = [ 0.0, 0.0, 0.0 ] 
                grid['points'][x][y][z]['avgori'] = 0
                grid['points'][x][y][z]['count'] = 0

                
                grid_counter += 1

    print "# There are %s grid points ( dimensions %s %s %s step %s ) " % ( grid_counter, nx, ny, nz, grid_step)        
    
    #print grid        

    #sys.exit(-1)

    grid['npoints'] = grid_counter

    return grid

def getCellData(cell_dir, cell_interval):
    
    cell_files = []

    files = os.listdir(cell_dir)
    
    for f in files:
        if f.startswith("cell") and f.endswith(".bin"):
            elems = f.split("_")
            num = elems[2].strip(".bin")
            num = float(num)
            if num >= cell_interval[0] and num <= cell_interval[1]:
                cell_files.append(os.path.join(cell_dir,f))
    
    cell_files.sort()

    return cell_files

def loadCell( cell_file ):

    waters_energies = [ 0.0, 0.0, 0.0 ]
    waters_solenergies = [0.0, 0.0, 0.0]
    waters_aforces = [0.0, 0.0, 0.0]
    waters_atorques = [0.0, 0.0, 0.0 ]
    waters_aori = 0

    nwaters = 0

    cell = {}
    cell['water'] = {}

    stream = open(cell_file, 'r')

    s = stream.read(8)
    framenumber, volume = struct.unpack("if", s)

    cell['frame'] = framenumber
    cell['volume'] = volume


    totnrg = 0.0
    s = stream.read(56)
    while ( len(s) == 56):
        resnum, watflag, cljnrg, cljsolnrg, fx, fy, fz, tx, ty, tz, ox, oy, oz, ori = struct.unpack('<2i12f', s)
        #print "<<<<", resnum, watflag, cljnrg, cljsolnrg, fx, fy, fz, tx, ty, tz, ox, oy, oz, ori, ">>>>"
        #sys.exit(-1)
        #if ori < 1 and watflag:
        #    #print "# low ori ( %s ) for water %s in cell file %s " % ( ori, resnum, cell_file)
        #    #ori = 1
        #    pass

        if watflag:
            cell['water'][resnum] = {}
            cell['water'][resnum]['energies'] = ( cljnrg, 0, 0 )
            cell['water'][resnum]['solenergies'] = ( cljsolnrg, 0, 0)
            cell['water'][resnum]['forces'] = ( fx, fy, fz ) 
            cell['water'][resnum]['torques'] = ( tx, ty, tz )
            cell['water'][resnum]['o_coords'] = ( ox, oy, oz )
            cell['water'][resnum]['ori'] = ori          
            nwaters += 1
            
            waters_energies[0] += cljnrg
            waters_solenergies[0] += cljsolnrg

            waters_aforces[0] += fx
            waters_aforces[1] += fy
            waters_aforces[2] += fz

            waters_atorques[0] += tx
            waters_atorques[1] += ty
            waters_atorques[2] += tz

            waters_aori += ori        

        s = stream.read(56)

    s = stream.read(56)
    if len(s) != 0:
        print "Cell file %s seems corrupted. Aborting." % cell_file
        sys.exit(-1)

    if nwaters != 0:
       
        for i in xrange(0,3):
            waters_energies[i] /= nwaters
            waters_solenergies[i] /= nwaters
            waters_aforces[i] /= nwaters
            waters_atorques[i] /= nwaters

        waters_aori /= float( nwaters )

    #print (waters_energies, waters_solenergies, waters_aforces)
    #sys.exit(-1)


    cell['waters_energies'] = waters_energies
    cell['waters_solenergies'] = waters_solenergies
    cell['waters_aforces'] = waters_aforces
    cell['waters_atorques'] = waters_atorques
    cell['waters_aori'] = waters_aori
    cell['nwaters'] = nwaters

    return cell


def updateGridPropsCumulative( grid, cell):

    step = grid['step']
    nx = grid['dimensions'][0]
    ny = grid['dimensions'][1]
    nz = grid['dimensions'][2]
    originx = grid['origin'][0]
    originy = grid['origin'][1]
    originz = grid['origin'][2]

    waters = cell['water'].keys()
    
    for water in waters:
        ox = cell['water'][water]['o_coords'][0]
        oy = cell['water'][water]['o_coords'][1]
        oz = cell['water'][water]['o_coords'][2]        

        # Convert cartesian into grid coordinates 
        gridx = int ( round ( ( ox - originx ) / step ) )
        gridy = int ( round ( ( oy - originy ) / step ) )
        gridz = int ( round ( ( oz - originz ) / step ) )

        # Ignore water molecules outside of grid
        if ( gridx >= nx or gridy >= ny or gridz >= nz  or 
             gridx < 0 or gridy < 0 or gridz < 0):
            #print "WARNING water %s outside of grid, it will be ignored " % (water)
            continue
        
        # We accumulate arithmetic mean of forces/torques/coordination numbers, and arithmetic mean of energies
        cumweight = grid['points'][gridx][gridy][gridz]['count'] + 1.0
        grid['points'][gridx][gridy][gridz]['count'] = cumweight
	newweight = 1.0

	bigratio = ( cumweight - newweight ) / cumweight
   	smallratio = 1.0 - bigratio
        
        # Energies
        oldavgnrg =  grid['points'][gridx][gridy][gridz]['avgnrg'][0]
        newnrg =  cell['water'][water]['energies'][0]
        grid['points'][gridx][gridy][gridz]['avgnrg'][0] = oldavgnrg*bigratio + newnrg*smallratio

        oldavgsolnrg =  grid['points'][gridx][gridy][gridz]['avgsolnrg'][0]
        newsolnrg =  cell['water'][water]['solenergies'][0]
        grid['points'][gridx][gridy][gridz]['avgsolnrg'][0] = oldavgsolnrg*bigratio + newsolnrg*smallratio

        # Forces
        oldavgforces = grid['points'][gridx][gridy][gridz]['avgforces']
        newforces = [ 0.0, 0.0, 0.0 ]
        newavgforces = [ 0.0, 0.0, 0.0]
        for m in range(0,3):
            newforces[m] = cell['water'][water]['forces'][m] 
            newavgforces[m] = oldavgforces[m] * bigratio + newforces[m] * smallratio
        grid['points'][gridx][gridy][gridz]['avgforces'] = newavgforces

        # Torques
        oldavgtorques = grid['points'][gridx][gridy][gridz]['avgtorques']
        newtorques = [ 0.0, 0.0, 0.0 ]
        newavgtorques = [ 0.0, 0.0, 0.0]
        for m in range(0,3):
            newtorques[m] = cell['water'][water]['torques'][m]
            newavgtorques[m] = oldavgtorques[m] * bigratio + newtorques[m] * smallratio
        grid['points'][gridx][gridy][gridz]['avgtorques'] = newavgtorques

        # Orientational numbers
        oldavgori = grid['points'][gridx][gridy][gridz]['avgori']
        newori = cell['water'][water]['ori']
        newavgori = oldavgori*bigratio + newori*smallratio
        grid['points'][gridx][gridy][gridz]['avgori'] = newavgori

    grid['nsnapshots'] += 1

def updatePropsCumulative( water_props, cell):

    #
    # We accumulate the arithmetic average of the forces, torques 
    #

    if len(water_props) == 0:
        water_props['density'] = 0.0
        water_props['forces'] = [ 0.0, 0.0, 0.0 ]
        water_props['torques'] =  [ 0.0, 0.0, 0.0 ]
        water_props['ori'] = 0.0

        water_props['nrg'] = 0.0
        water_props['solnrg'] = 0.0
        water_props['cum_weight'] = 0.0
        water_props['avg_w'] = 0.0

    newweight = cell['nwaters']

    water_props['cum_weight'] += newweight

    cumweight = water_props['cum_weight']

    bigratio = ( cumweight - newweight ) / cumweight
    smallratio = 1.0 - bigratio

    oldavgw = water_props['avg_w']
    neww = cell['nwaters']
    newavgw = oldavgw*bigratio + neww*smallratio
    water_props['avg_w'] = newavgw

    oldavgnrg = water_props['nrg']
    newnrg =  cell['waters_energies'][0]
    newavgnrg = oldavgnrg*bigratio + newnrg*smallratio
    water_props['nrg'] = newavgnrg

    oldavgsolnrg = water_props['solnrg']
    newsolnrg =  cell['waters_solenergies'][0]
    newavgsolnrg = oldavgsolnrg*bigratio + newsolnrg*smallratio
    water_props['solnrg'] = newavgsolnrg


    # Arithmetic mean forces/torques/ori/density
    oldavgforces =  water_props['forces']
    newavgforces = [0.0,0.0,0.0]

    # average of forces and torques
    newforces = [ 0.0, 0.0, 0.0 ]
    # print water_num
    for m in range(0,3):
        newforces[m] = cell['waters_aforces'][m]

    oldavgtorques =  water_props['torques']
    newavgtorques = [0.0,0.0,0.0]
    newtorques = [ 0.0, 0.0, 0.0]
    for m in range(0,3):
        newtorques[m] = cell['waters_atorques'][m]

    for m in range(0,3):
        newavgforces[m] = oldavgforces[m] * bigratio + newforces[m] * smallratio
        newavgtorques[m] = oldavgtorques[m] * bigratio + newtorques[m] * smallratio

    water_props['forces'] = newavgforces
    water_props['torques'] = newavgtorques

    oldavgori =  water_props['ori']
    newori = cell['waters_aori']
    newavgori = oldavgori*bigratio + newori*smallratio
    water_props['ori'] = newavgori
    
    # density..
    oldavgdensity = water_props['density']
    newdensity = calcDensity( cell['volume'], cell['nwaters'] ) 
    newavgdensity = oldavgdensity*bigratio + newdensity*smallratio
    water_props['density'] = newavgdensity


def calcAverageProps( waters_props):

    vibratio = 1.
    libratio = 1.
    oriratio = 1.

    DHWW = 0.0
 
    NW = water_props['avg_w']

    # average over all waters
    meanF = [ 0.0, 0.0, 0.0 ]
    meanT = [ 0.0, 0.0, 0.0 ]
    meanori = 0.0

    for m in range(0,3):
       meanF[m] = waters_props['forces'][m] 
       meanT[m] = waters_props['torques'][m] 

    # !!!! The forces from the cell files are in kcal/mol/A 
    # However, for consistency with the literature, the BULK forces are 
    # in N/m x 1e-10 so we need to convert
    conv_factor = ( KCAL_TO_J / NAVOGADRO ) * 1e20 
    for m in range(0,3):
       meanF[m] = meanF[m] * conv_factor
       meanT[m] = meanT[m] * conv_factor

    for m in range(0,3):
        if meanF[m] == 0:
            vibratio = 1.0
            break
        if meanT[m] == 0:
            libratio = 1.0
            break
        vibratio *= ( BULKFORCES[m] / meanF[m] ) 
        libratio *= ( BULKTORQUES[m] / meanT[m] )

    ori = waters_props['ori']
    if ori < SMALL:
        ori = SMALL
    oriratio *=  ( ori / BULKORI )

    #print "vibratio is %s " % vibratio
    #print "libratio is %s " % libratio
    #print "oriratio is %s " % oriratio

    DHSOL = NW * ( waters_props['solnrg'] )
 
    DHSLV = NW * ( waters_props['nrg'] - BULKNRG )

    meanwwnrg = waters_props['nrg']

    meanori = ori

    DS_vib = -T * KB * NW * log(vibratio)

    DS_lib = -T * KB * NW * log(libratio)

    DS_ori = -T * KB * NW * log(oriratio)

    DSW = DS_vib + DS_lib + DS_ori 

    DGW = DHSOL + DHSLV + DSW

    meandensity = water_props['density']
    #print meandensity
    #sys.exit(-1)

    return DS_vib, DS_lib, DS_ori, DSW, DHSOL, DHSLV, DGW, meanF, meanT, meanwwnrg, meanori, meandensity


#def calcPointProps( point, volume, nsnapshots ):
def calcPointProps( point ):

    vibratio = 1.
    libratio = 1.
    oriratio = 1.

    NW = 1.0

    meanF = [ 0.0, 0.0, 0.0 ]
    meanT = [ 0.0, 0.0, 0.0 ]

    for m in range(0,3):
       meanF[m] = point['avgforces'][m]
       meanT[m] = point['avgtorques'][m]

    # !!!! The forces from the cell files are in kcal/mol/A 
    # However, for consistency with the literature, the BULK forces are 
    # in N/m x 1e-10 so we need to convert
    conv_factor = ( KCAL_TO_J / NAVOGADRO ) * 1e20
    for m in range(0,3):
       meanF[m] = meanF[m] * conv_factor
       meanT[m] = meanT[m] * conv_factor

    for m in range(0,3):
        if meanF[m] == 0:
            vibratio = 1
            break
        if meanT[m] == 0:
            libratio = 1
            break
        vibratio *= ( BULKFORCES[m] / meanF[m] ) 
        libratio *= ( BULKTORQUES[m] / meanT[m] )
    
    ori = point['avgori']
    if (ori < SMALL):
        ori = SMALL 
    oriratio *=  ( ori / BULKORI )

    DHSOL = ( point['avgsolnrg'][0] ) * NW

    DHSLV = ( point['avgnrg'][0] - BULKNRG ) * NW
 
    DS_vib = -T * KB * log(vibratio) * NW
    DS_lib = -T * KB * log(libratio) * NW
    DS_ori = -T * KB * log(oriratio) * NW
    
    #density = calcDensity( volume,  point['count'] / nsnapshots )
    #
    #relative_density = density / BULKDENSITY

    #return relative_density, DS_vib, DS_lib, DS_ori, DH
    return DS_vib, DS_lib, DS_ori, DHSOL, DHSLV

def calcDensity(volume, nwats):
    if nwats != 0:
        water_per_cubicangstrom = nwats / volume
        cubicdm_per_mol_water = (1 / water_per_cubicangstrom) * NAVOGADRO * 1e-27
        water_density =  (1/cubicdm_per_mol_water)*MASS_WATER # g/L    
    else:
        water_density = 0
    #print water_density
    #sys.exit(-1)
    return water_density

def outputGrid( grid, grid_dir ):
    #
    # Format a PDB style file with the cartesian coordinate of the grid points. Use columns like 
    # occupancy or B-factors to write the entropies enthalpies etc...
    #

    titlepdbDG='%s/grid_DG.pdb' % (grid_dir)
    filinpdbDG=open( titlepdbDG, 'w')

    titlepdbDH='%s/grid_DH.pdb' % (grid_dir)
    filinpdbDH=open( titlepdbDH, 'w')

    titlepdbDHsol='%s/grid_DHsolute.pdb' % (grid_dir)
    filinpdbDHsol=open( titlepdbDHsol, 'w')

    titlepdbDHslv='%s/grid_DHsolvent.pdb' % (grid_dir)
    filinpdbDHslv=open( titlepdbDHslv, 'w')

    titlepdbDS='%s/grid_minTDS.pdb' % (grid_dir)
    filinpdbDS=open( titlepdbDS, 'w')

    titlepdbDSvib='%s/grid_minTDSvib.pdb' % (grid_dir)
    filinpdbDSvib=open( titlepdbDSvib, 'w')

    titlepdbDSlib='%s/grid_minTDSlib.pdb' % (grid_dir)
    filinpdbDSlib=open( titlepdbDSlib, 'w')

    titlepdbDSori='%s/grid_minTDSori.pdb' % (grid_dir)
    filinpdbDSori=open( titlepdbDSori, 'w')

    titlepdboccupancy='%s/grid_density.pdb' % (grid_dir)
    filinpdboccupancy=open( titlepdboccupancy, 'w')

    nx, ny, nz = grid['dimensions']
    originx, originy, originz = grid['origin']
    step = grid['step']
    point_volume = step * step * step

    nsnapshots = grid['nsnapshots']

    index=0

    for x in range(0,nx):
        for y in range(0,ny):
            for z in range(0,nz):
                count = grid['points'][x][y][z]['count'] 
                if count > grid_count_cutoff:
                    abs_density = calcDensity( point_volume,  count / nsnapshots )
                    density = abs_density / BULKDENSITY
                    minTDSvib, minTDSlib, minTDSori, DHSOL, DHSLV = calcPointProps( grid['points'][x][y][z] )
                else:
                    density = 0.0
                    DH = 0.0
                    minTDSvib = 0.0
                    minTDSlib = 0.0
                    minTDSori = 0.0
                    DHSOL = 0.0
                    DHSLV = 0.0
                
                index += 1
                
                minTDS = minTDSvib + minTDSlib + minTDSori 
                DH = DHSOL + DHSLV
                DG = DH + minTDS

                grid['points'][x][y][z]['delta_G'] = DG
                grid['points'][x][y][z]['delta_H'] = DH
                grid['points'][x][y][z]['delta_H_sol'] = DHSOL
                grid['points'][x][y][z]['delta_H_slv'] = DHSLV
                grid['points'][x][y][z]['minT_delta_S'] = minTDS
                grid['points'][x][y][z]['minT_delta_S_vib'] = minTDSvib
                grid['points'][x][y][z]['minT_delta_S_lib'] = minTDSlib
                grid['points'][x][y][z]['minT_delta_S_ori'] = minTDSori
                grid['points'][x][y][z]['density'] = density
 
                cartx = originx + step * x
                carty = originy + step * y
                cartz = originz + step * z     

                textDG='ATOM%7s%2s%24.2f%8.2f%8.2f%5s %8.4f %10.4f%12s' % (index, 'H', cartx, carty, cartz, '1.00', DG, density, 'H')
                textDH='ATOM%7s%2s%24.2f%8.2f%8.2f%5s %8.4f %10.4f%12s' % (index, 'H', cartx, carty, cartz, '1.00', DH, density, 'H')
                textDHsol='ATOM%7s%2s%24.2f%8.2f%8.2f%5s %8.4f %10.4f%12s' % (index, 'H', cartx, carty, cartz, '1.00', DHSOL, density, 'H')
                textDHslv='ATOM%7s%2s%24.2f%8.2f%8.2f%5s %8.4f %10.4f%12s' % (index, 'H', cartx, carty, cartz, '1.00', DHSLV, density, 'H')
                textDS='ATOM%7s%2s%24.2f%8.2f%8.2f%5s %8.4f %10.4f%12s' % (index, 'H', cartx, carty, cartz, '1.00', minTDS, density, 'H')
                textDSvib='ATOM%7s%2s%24.2f%8.2f%8.2f%5s %8.4f %10.4f%12s' % (index, 'H', cartx, carty, cartz, '1.00', minTDSvib, density, 'H')
                textDSlib='ATOM%7s%2s%24.2f%8.2f%8.2f%5s %8.4f %10.4f%12s' % (index, 'H', cartx, carty, cartz, '1.00', minTDSlib, density, 'H')                    
                textDSori='ATOM%7s%2s%24.2f%8.2f%8.2f%5s %8.4f %10.4f%12s' % (index, 'H', cartx, carty, cartz, '1.00', minTDSori, density, 'H')

                textocc='ATOM%7s%2s%24.2f%8.2f%8.2f%5s %8.4f %10.4f%12s' % (index, 'H', cartx, carty, cartz, '1.00', density, count, 'H') 

                filinpdbDG.write( textDG + "\n" )
                filinpdbDH.write( textDH + "\n" )
                filinpdbDHsol.write( textDHsol + "\n" )
                filinpdbDHslv.write( textDHslv + "\n" )
                filinpdbDS.write( textDS + "\n" )  
                filinpdbDSvib.write( textDSvib + "\n" )  
                filinpdbDSlib.write( textDSlib + "\n" )  
                filinpdbDSori.write( textDSori + "\n" )                      
                  
                filinpdboccupancy.write( textocc + "\n" )
                    
    #sys.exit(-1)

    filinpdbDG.close()
    filinpdbDH.close()
    filinpdbDHsol.close()
    filinpdbDHslv.close()
    filinpdbDS.close()
    filinpdbDSvib.close()
    filinpdbDSlib.close()
    filinpdbDSori.close()
    filinpdboccupancy.close()

    # Format a DX style file with the cartesian coordinate of the grid points. Use columns like 
    # occupancy or B-factors to write the entropies enthalpies etc...
    #
    
    ### keep list of each data file for moe
    dg_moe = ""
    dh_moe = ""
    dhsolute_moe = ""
    dhsolvent_moe = ""
    mintds_moe = ""
    mintdsvib_moe = ""
    mintdslib_moe = ""
    mintdsori_moe = ""
    occupancy_moe = ""


    titledxDG='%s/grid_DG.dx' % (grid_dir)
    filindxDG=open( titledxDG, 'w')

    titledxDH='%s/grid_DH.dx' % (grid_dir)
    filindxDH=open( titledxDH, 'w')

    titledxDHsol='%s/grid_DHsolute.dx' % (grid_dir)
    filindxDHsol=open( titledxDHsol, 'w')

    titledxDHslv='%s/grid_DHsolvent.dx' % (grid_dir)
    filindxDHslv=open( titledxDHslv, 'w')

    titledxDS='%s/grid_minTDS.dx' % (grid_dir)
    filindxDS=open( titledxDS, 'w')

    titledxDSvib='%s/grid_minTDSvib.dx' % (grid_dir)
    filindxDSvib=open( titledxDSvib, 'w')

    titledxDSlib='%s/grid_minTDSlib.dx' % (grid_dir)
    filindxDSlib=open( titledxDSlib, 'w')

    titledxDSori='%s/grid_minTDSori.dx' % (grid_dir)
    filindxDSori=open( titledxDSori, 'w')

    titledxoccupancy='%s/grid_density.dx' % (grid_dir)
    filindxoccupancy=open( titledxoccupancy, 'w')

    nx, ny, nz = grid['dimensions']
    originx, originy, originz = grid['origin']
    step = grid['step']

    index=0

    ######FORMAT an appropriate DX file for both DH and DS ###############
    #####header of DX file ##############
    #object 1 class gridpositions counts 16 17 15
    line1='object 1 class gridpositions counts %i %i %i\n' % (grid['dimensions'][0], grid['dimensions'][1], grid['dimensions'][2])
    #origin -1.0 -1.0  0.0
    line2='origin %4.2f %4.2f %4.2f\n' % (grid['origin'][0], grid['origin'][1], grid['origin'][2])
    #delta 1 0 0
    line3='delta %4.2f 0 0\n' % (grid['step'])
    #delta 0 1 0
    line4='delta 0 %4.2f 0\n' % (grid['step'])
    #delta 0 0 1
    line5='delta 0 0 %4.2f\n' % (grid['step'])
    #object 2 class gridconnections counts 16 17 15 
    line6='object 2 class gridconnections counts %4.2f %4.2f %4.2f\n' % (grid['dimensions'][0], grid['dimensions'][1], grid['dimensions'][2]) 
    #object 3 class array type double rank 0 items 4080 data follows   
    line7='object 3 class array type double rank 0 items %i data follows\n' % (grid['dimensions'][0]*grid['dimensions'][1]*grid['dimensions'][2])  

    filindxDS.write(line1+line2+line3+line4+line5+line6+line7)
    filindxDH.write(line1+line2+line3+line4+line5+line6+line7)
    filindxDHsol.write(line1+line2+line3+line4+line5+line6+line7)
    filindxDHslv.write(line1+line2+line3+line4+line5+line6+line7)
    filindxDG.write(line1+line2+line3+line4+line5+line6+line7)
    filindxDSvib.write(line1+line2+line3+line4+line5+line6+line7)
    filindxDSlib.write(line1+line2+line3+line4+line5+line6+line7)
    filindxDSori.write(line1+line2+line3+line4+line5+line6+line7)
    filindxoccupancy.write(line1+line2+line3+line4+line5+line6+line7)

    strdh=""
    strdhsol=""
    strdhslv=""
    strds=""
    strdg=""
    strocc=""
    strvib=""
    strlib=""
    strori=""

    entries=0
    cond = 0

    for x in range(0,nx):
        for y in range(0,ny):
            for z in range(0,nz):
                index += 1
                count = grid['points'][x][y][z]['count']
                if count > grid_count_cutoff:
                    density = grid['points'][x][y][z]['density']
                    DG = grid['points'][x][y][z]['delta_G'] 
                    DH = grid['points'][x][y][z]['delta_H'] 
                    DHSOL = grid['points'][x][y][z]['delta_H_sol'] 
                    DHSLV = grid['points'][x][y][z]['delta_H_slv'] 
                    minTDS = grid['points'][x][y][z]['minT_delta_S']
                    minTDSvib = grid['points'][x][y][z]['minT_delta_S_vib'] 
                    minTDSlib = grid['points'][x][y][z]['minT_delta_S_lib'] 
                    minTDSori = grid['points'][x][y][z]['minT_delta_S_ori'] 
                else:
                    density = 0.0
                    DH = 0.0
                    DHSOL = 0.0
                    DHSLV = 0.0
                    minTDSvib = 0.0
                    minTDSlib = 0.0
                    minTDSori = 0.0
                    DH = 0.0

                minTDS = minTDSvib + minTDSlib + minTDSori
                DG = DH + minTDS

                strdh += " %8.4f" % DH
                strdhsol += " %8.4f" % DHSOL
                strdhslv += " %8.4f" % DHSLV
                strds += " %8.4f" % minTDS   
                strdg += " %8.4f" % DG
                strocc += " %8.4f" % density   
                strvib += " %8.4f" % minTDSvib   
                strlib += " %8.4f" % minTDSlib
                strori += " %8.4f" % minTDSori  

                if index == 1:
                    dg_moe += "[%8.4f," % DG
                    dh_moe += "[%8.4f," % DH
                    dhsolute_moe += "[%8.4f," % DHSOL
                    dhsolvent_moe += "[%8.4f," % DHSLV
                    mintds_moe += "[%8.4f," % minTDS
                    mintdsvib_moe += "[%8.4f," % minTDSvib
                    mintdslib_moe += "[%8.4f," % minTDSlib
                    mintdsori_moe += "[%8.4f," % minTDSori
                    occupancy_moe += "[%8.4f," % density
                if index > 1:
                    dg_moe += "%8.4f," % DG
                    dh_moe += "%8.4f," % DH
                    dhsolute_moe += "%8.4f," % DHSOL
                    dhsolvent_moe += "%8.4f," % DHSLV
                    mintds_moe += "%8.4f," % minTDS
                    mintdsvib_moe += "%8.4f," % minTDSvib
                    mintdslib_moe += "%8.4f," % minTDSlib
                    mintdsori_moe += "%8.4f," % minTDSori
                    occupancy_moe += "%8.4f," % density
                entries += 1
                if (entries % 3) == 0:
                    strdh += "\n"
                    strdhsol += "\n"
                    strdhslv += "\n"
                    strds += "\n"
                    strdg += "\n"
                    strocc += "\n" 
                    strvib += "\n"
                    strlib += "\n"
                    strori += "\n"

    filindxDG.write( strdg + "\n" )
    filindxDH.write( strdh + "\n" )
    filindxDHsol.write( strdhsol + "\n" )
    filindxDHslv.write( strdhslv + "\n" )
    filindxDS.write( strds + "\n" )  
    filindxDSvib.write( strvib + "\n" )  
    filindxDSlib.write( strlib + "\n" )  
    filindxDSori.write( strori + "\n" )                      
                 
    filindxoccupancy.write( strocc + "\n" )

    filindxDG.close()
    filindxDH.close()
    filindxDHsol.close()
    filindxDHslv.close()
    filindxDS.close()
    filindxDSvib.close()
    filindxDSlib.close()
    filindxDSori.close()
    filindxoccupancy.close()

    #
    # Region file for entire grid
    #
    region = open( os.path.join( grid_dir, "all.region"), 'w' )
    regstr = ""
    for x in range(1,entries+1):
        regstr += " %d " % x
        if (x % 10) == 0:
            regstr += "\n"
    regstr += "\n"
    region.write(regstr)
    region.close()


    ##### ALSO output MOE grid files.  Which has two files one shape file for both weighted and unweighted MOE grid files 
    ##### this contains the [x1,x2,....xn],[y1,y2,.....yn], [z1,z2....zn] of all the grid points
    ##### All data points will be labeled each with simply a file with [d_x1_y1_z1, d_x1_y1_z2....d_xn_yn_zn]
    moe_shape = open( os.path.join( grid_dir, "shape"), 'w' )
    regstr = "[%s,%s,%s]" % (nx,ny,nz)
    moe_shape.write(regstr)
    moe_shape.close()

    moe_outputs = {"dg_moe" : dg_moe, "dh_moe" : dh_moe, "dhsolute_moe" : dhsolute_moe, "dhsolvent_moe" : dhsolvent_moe, "mintds_moe" : mintds_moe, "mintdsvib_moe" : mintdsvib_moe, "mintdslib_moe" : mintdslib_moe, "mintdsori_moe" : mintdsori_moe, "occupancy_moe" : occupancy_moe}

    for i in moe_outputs.keys():
        output = open( os.path.join( grid_dir, i), 'w')
        moestr = moe_outputs[i] + "]\n"
        output.write(moestr)
        output.close()




def outputGridParams( grid, grid_dir  ):

    nx, ny, nz = grid['dimensions']
    originx, originy, originz = grid['origin']
    step = grid['step']
    point_volume = step * step * step
    nsnapshots = grid['nsnapshots']
    
    kcalmolA_to_Nm = ( KCAL_TO_J / NAVOGADRO ) * 1e20

    index=0

    outfile = os.path.join(grid_dir, "grid.forces") 
    gridforces = open( outfile, 'w')

    # grid details
    gridforces.write("#Grid Dimensions Origin Step \n")
    line = "%8d\t%8d\t%8d\t%10.5f\t%10.5f\t%10.5f\t%10.5f\n" % (nx, ny, nz, originx, originy, originz, step)
    gridforces.write(line) 
    # Could add standard deviation of each parameter
    gridforces.write("#Point\tX\tY\tZ\tFx\tFy\tFz\tTx\tTy\tTz\tOri\tDensity\tVolume\tEnergy\tSolEnergy\n")

    for x in range(0,nx):
        for y in range(0,ny):
            for z in range(0,nz):
                index += 1
                count = grid['points'][x][y][z]['count']
                if count > 0:
                    density = calcDensity( point_volume, count / nsnapshots )
                else:
                    density = 0
                meanF = grid['points'][x][y][z]['avgforces']
                meanT = grid['points'][x][y][z]['avgtorques']
                volume = point_volume
                # !!!! The forces from the cell files are in kcal/mol/A 
                # However, for consistency with the literature, the BULK forces are 
                # in N/m x  1e-10 so we need to convert
                for m in range(0,3):
                    meanF[m] = meanF[m] * kcalmolA_to_Nm 
                    meanT[m] = meanT[m] * kcalmolA_to_Nm 
                meanori = grid['points'][x][y][z]['avgori']
                meannrg = grid['points'][x][y][z]['avgnrg'][0]
                meansolnrg = grid['points'][x][y][z]['avgsolnrg'][0]
                #cartx = originx + step * x
                #carty = originy + step * y
                #cartz = originz + step * z
 
                line = "%9d %8d %8d %8d %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %12.3f %10.5f %10.5f %10.5f\n" % (index,x,y,z,meanF[0],meanF[1],meanF[2],meanT[0],meanT[1],meanT[2],meanori,density,volume,meannrg, meansolnrg)
                gridforces.write(line)                
    gridforces.close()

def outputWeightedGrid( grid, grid_dir ):
    #
    # Format a PDB style file with the cartesian coordinate of the grid points. Use columns like 
    # occupancy or B-factors to write the entropies enthalpies etc...
    #

    titlepdbDG='%s/grid_DG_weighted.pdb' % (grid_dir)
    filinpdbDG=open( titlepdbDG, 'w')

    titlepdbDH='%s/grid_DH_weighted.pdb' % (grid_dir)
    filinpdbDH=open( titlepdbDH, 'w')

    titlepdbDHSOL='%s/grid_DHsolute_weighted.pdb' % (grid_dir)
    filinpdbDHSOL=open( titlepdbDHSOL, 'w')

    titlepdbDHSLV='%s/grid_DHsolvent_weighted.pdb' % (grid_dir)
    filinpdbDHSLV=open( titlepdbDHSLV, 'w')

    titlepdbDS='%s/grid_minTDS_weighted.pdb' % (grid_dir)
    filinpdbDS=open( titlepdbDS, 'w')

    titlepdbDSvib='%s/grid_minTDSvib_weighted.pdb' % (grid_dir)
    filinpdbDSvib=open( titlepdbDSvib, 'w')

    titlepdbDSlib='%s/grid_minTDSlib_weighted.pdb' % (grid_dir)
    filinpdbDSlib=open( titlepdbDSlib, 'w')

    titlepdbDSori='%s/grid_minTDSori_weighted.pdb' % (grid_dir)
    filinpdbDSori=open( titlepdbDSori, 'w')

    #titlepdboccupancy='%s/grid_density.pdb' % (grid_dir)
    #filinpdboccupancy=open( titlepdboccupancy, 'w')

    nx, ny, nz = grid['dimensions']
    originx, originy, originz = grid['origin']
    step = grid['step']
    point_volume = step * step * step

    nsnapshots = grid['nsnapshots']

    index=0

    for x in range(0,nx):
        for y in range(0,ny):
            for z in range(0,nz):
                count = grid['points'][x][y][z]['count'] 
                if count > 0:
                    abs_density = calcDensity( point_volume,  count / nsnapshots )
                    density = abs_density / BULKDENSITY
                    molpervoxel = density * BULKWATMOLPERVOXEL
                    minTDSvib, minTDSlib, minTDSori, DHSOL, DHSLV = calcPointProps( grid['points'][x][y][z] )
                else:
                    density = 0.0
                    molpervoxel = 0.0
                    DHSOL = 0.0
                    DHSLV = 0.0
                    minTDSvib = 0.0
                    minTDSlib = 0.0
                    minTDSori = 0.0

                minTDS = minTDSvib + minTDSlib + minTDSori
                DH = DHSOL + DHSLV
                DG = DH + minTDS

                index += 1
                  
                grid['points'][x][y][z]['delta_Gw'] = DG*molpervoxel
                grid['points'][x][y][z]['delta_Hw'] = DH*molpervoxel
                grid['points'][x][y][z]['delta_H_solw'] = DHSOL*molpervoxel
                grid['points'][x][y][z]['delta_H_slvw'] = DHSLV*molpervoxel
                grid['points'][x][y][z]['minT_delta_Sw'] = minTDS*molpervoxel
                grid['points'][x][y][z]['minT_delta_S_vibw'] = minTDSvib*molpervoxel
                grid['points'][x][y][z]['minT_delta_S_libw'] = minTDSlib*molpervoxel
                grid['points'][x][y][z]['minT_delta_S_oriw'] = minTDSori*molpervoxel
                grid['points'][x][y][z]['densityw'] = density       
 
                cartx = originx + step * x
                carty = originy + step * y
                cartz = originz + step * z     

                textDG='ATOM%7s%2s%24.2f%8.2f%8.2f%5s %10.4f %10.4f%12s' % (index, 'H', cartx, carty, cartz, '1.00', DG*molpervoxel, density, 'H')
                textDH='ATOM%7s%2s%24.2f%8.2f%8.2f%5s %8.4f %10.4f%12s' % (index, 'H', cartx, carty, cartz, '1.00', DH*molpervoxel, density, 'H')
                textDHSOL='ATOM%7s%2s%24.2f%8.2f%8.2f%5s %8.4f %10.4f%12s' % (index, 'H', cartx, carty, cartz, '1.00', DHSOL*molpervoxel, density, 'H')
                textDHSLV='ATOM%7s%2s%24.2f%8.2f%8.2f%5s %8.4f %10.4f%12s' % (index, 'H', cartx, carty, cartz, '1.00', DHSLV*molpervoxel, density, 'H')
                textDS='ATOM%7s%2s%24.2f%8.2f%8.2f%5s %8.4f %10.4f%12s' % (index, 'H', cartx, carty, cartz, '1.00', minTDS*molpervoxel, density, 'H')
                textDSvib='ATOM%7s%2s%24.2f%8.2f%8.2f%5s %8.4f %10.4f%12s' % (index, 'H', cartx, carty, cartz, '1.00', minTDSvib*molpervoxel, density, 'H')
                textDSlib='ATOM%7s%2s%24.2f%8.2f%8.2f%5s %8.4f %10.4f%12s' % (index, 'H', cartx, carty, cartz, '1.00', minTDSlib*molpervoxel, density, 'H')
                textDSori='ATOM%7s%2s%24.2f%8.2f%8.2f%5s %8.4f %10.4f%12s' % (index, 'H', cartx, carty, cartz, '1.00', minTDSori*molpervoxel, density, 'H')

                filinpdbDG.write( textDG + "\n" )
                filinpdbDH.write( textDH + "\n" )
                filinpdbDHSOL.write( textDHSOL + "\n" )
                filinpdbDHSLV.write( textDHSLV + "\n" )
                filinpdbDS.write( textDS + "\n" )  
                filinpdbDSvib.write( textDSvib + "\n" )  
                filinpdbDSlib.write( textDSlib + "\n" )  
                filinpdbDSori.write( textDSori + "\n" )                      
                
    #sys.exit(-1)

    filinpdbDG.close()
    filinpdbDH.close()
    filinpdbDHSOL.close()
    filinpdbDHSLV.close()
    filinpdbDS.close()
    filinpdbDSvib.close()
    filinpdbDSlib.close()
    filinpdbDSori.close()

    # Format a DX style file with the cartesian coordinate of the grid points. Use columns like 
    # occupancy or B-factors to write the entropies enthalpies etc...
    #
    titledxDG='%s/grid_DG_weighted.dx' % (grid_dir)
    filindxDG=open( titledxDG, 'w')

    titledxDH='%s/grid_DH_weighted.dx' % (grid_dir)
    filindxDH=open( titledxDH, 'w')

    titledxDHSOL='%s/grid_DHsolute_weighted.dx' % (grid_dir)
    filindxDHSOL=open( titledxDHSOL, 'w')

    titledxDHSLV='%s/grid_DHsolvent_weighted.dx' % (grid_dir)
    filindxDHSLV=open( titledxDHSLV, 'w')

    titledxDS='%s/grid_minTDS_weighted.dx' % (grid_dir)
    filindxDS=open( titledxDS, 'w')

    titledxDSvib='%s/grid_minTDSvib_weighted.dx' % (grid_dir)
    filindxDSvib=open( titledxDSvib, 'w')

    titledxDSlib='%s/grid_minTDSlib_weighted.dx' % (grid_dir)
    filindxDSlib=open( titledxDSlib, 'w')

    titledxDSori='%s/grid_minTDSori_weighted.dx' % (grid_dir)
    filindxDSori=open( titledxDSori, 'w')

    nx, ny, nz = grid['dimensions']
    originx, originy, originz = grid['origin']
    step = grid['step']

    index=0

    ######FORMAT an appropriate DX file for both DH and DS ###############
    #####header of DX file ##############
    #object 1 class gridpositions counts 16 17 15
    line1='object 1 class gridpositions counts %i %i %i\n' % (grid['dimensions'][0], grid['dimensions'][1], grid['dimensions'][2])
    #origin -1.0 -1.0  0.0
    line2='origin %4.2f %4.2f %4.2f\n' % (grid['origin'][0], grid['origin'][1], grid['origin'][2])
    #delta 1 0 0
    line3='delta %4.2f 0 0\n' % (grid['step'])
    #delta 0 1 0
    line4='delta 0 %4.2f 0\n' % (grid['step'])
    #delta 0 0 1
    line5='delta 0 0 %4.2f\n' % (grid['step'])
    #object 2 class gridconnections counts 16 17 15 
    line6='object 2 class gridconnections counts %i %i %i\n' % (grid['dimensions'][0], grid['dimensions'][1], grid['dimensions'][2]) 
    #object 3 class array type double rank 0 items 4080 data follows   
    line7='object 3 class array type double rank 0 items %i data follows\n' % (grid['dimensions'][0]*grid['dimensions'][1]*grid['dimensions'][2])  

    filindxDS.write(line1+line2+line3+line4+line5+line6+line7)
    filindxDH.write(line1+line2+line3+line4+line5+line6+line7)
    filindxDHSOL.write(line1+line2+line3+line4+line5+line6+line7)
    filindxDHSLV.write(line1+line2+line3+line4+line5+line6+line7)
    filindxDG.write(line1+line2+line3+line4+line5+line6+line7)
    filindxDSvib.write(line1+line2+line3+line4+line5+line6+line7)
    filindxDSlib.write(line1+line2+line3+line4+line5+line6+line7)
    filindxDSori.write(line1+line2+line3+line4+line5+line6+line7)

    ### keep list of each data file for moe
    dg_moe = ""
    dh_moe = ""
    dhsolute_moe = ""
    dhsolvent_moe = ""
    mintds_moe = ""
    mintdsvib_moe = ""
    mintdslib_moe = ""
    mintdsori_moe = ""
    occupancy_moe = ""

    strdh=""
    strdhsol=""
    strdhslv=""
    strds=""
    strdg=""
    strvib=""
    strlib=""
    strori=""

    entries=0
    cond = 0

    for x in range(0,nx):
        for y in range(0,ny):
            for z in range(0,nz):
                index += 1
                DG = grid['points'][x][y][z]['delta_Gw'] 
                DH = grid['points'][x][y][z]['delta_Hw'] 
                DHSOL = grid['points'][x][y][z]['delta_H_solw']
                DHSLV = grid['points'][x][y][z]['delta_H_slvw']

                minTDS = grid['points'][x][y][z]['minT_delta_Sw']
                minTDSvib = grid['points'][x][y][z]['minT_delta_S_vibw'] 
                minTDSlib = grid['points'][x][y][z]['minT_delta_S_libw'] 
                minTDSori = grid['points'][x][y][z]['minT_delta_S_oriw'] 
                density = grid['points'][x][y][z]['densityw'] 

                minTDS = minTDSvib + minTDSlib + minTDSori

                strdh += " %8.4f" % (DH)
                strdhsol += " %8.4f" % (DHSOL)
                strdhslv += " %8.4f" % (DHSLV)
                strds += " %8.4f" % (minTDS)
                strdg += " %8.4f" % (DG)
                strvib += " %8.4f" % (minTDSvib)
                strlib += " %8.4f" % (minTDSlib)
                strori += " %8.4f" % (minTDSori)

                if index == 1:
                    dg_moe += "[%8.4f, " % DG
                    dh_moe += "[%8.4f, " % DH
                    dhsolute_moe += "[%8.4f, " % DHSOL
                    dhsolvent_moe += "[%8.4f, " % DHSLV
                    mintds_moe += "[%8.4f, " % minTDS
                    mintdsvib_moe += "[%8.4f, " % minTDSvib
                    mintdslib_moe += "[%8.4f, " % minTDSlib
                    mintdsori_moe += "[%8.4f, " % minTDSori
                    occupancy_moe += "[%8.4f, " % density
                if index > 1:
                    dg_moe += "%8.4f, " % DG
                    dh_moe += "%8.4f, " % DH
                    dhsolute_moe += "%8.4f, " % DHSOL
                    dhsolvent_moe += "%8.4f, " % DHSLV
                    mintds_moe += "%8.4f, " % minTDS
                    mintdsvib_moe += "%8.4f, " % minTDSvib
                    mintdslib_moe += "%8.4f, " % minTDSlib
                    mintdsori_moe += "%8.4f, " % minTDSori
                    occupancy_moe += "%8.4f, " % density

                entries += 1
                if (entries % 3) == 0:
                    strdh += "\n"
                    strdhsol += "\n"
                    strdhslv += "\n"
                    strds += "\n"
                    strdg += "\n"
                    strvib += "\n"
                    strlib += "\n"
                    strori += "\n"

    filindxDG.write( strdg + "\n" )
    filindxDH.write( strdh + "\n" )
    filindxDHSOL.write( strdhsol + "\n" )
    filindxDHSLV.write( strdhslv + "\n" )
    filindxDS.write( strds + "\n" )  
    filindxDSvib.write( strvib + "\n" )  
    filindxDSlib.write( strlib + "\n" )  
    filindxDSori.write( strori + "\n" )                      
                 
    filindxDG.close()
    filindxDH.close()
    filindxDHSOL.close()
    filindxDHSLV.close()
    filindxDS.close()
    filindxDSvib.close()
    filindxDSlib.close()
    filindxDSori.close()

    moe_outputs = {"dg_moe_weighted" : dg_moe, "dh_moe_weighted" : dh_moe, "dhsolute_moe_weighted" : dhsolute_moe, "dhsolvent_moe_weighted" : dhsolvent_moe, "mintds_moe_weighted" : mintds_moe, "mintdsvib_moe_weighted" : mintdsvib_moe, "mintdslib_moe_weighted" : mintdslib_moe, "mintdsori_moe_weighted" : mintdsori_moe}
    
    for i in moe_outputs.keys():
        output = open( os.path.join( grid_dir, i), 'w')
        moestr = moe_outputs[i] + "]\n"
        output.write(moestr)
        output.close()


###################

if __name__ == "__main__":

    if os.path.exists( grid_dir ):
        if Overwrite:
            cmd = "rm -rf %s " % grid_dir
            os.system(cmd)
        else:
            print "Output folder %s exists, but overwrite mode has not been set. Aborting."
            sys.exit(-1)

    if  not os.path.exists( grid_dir ):
        cmd = "mkdir %s " % grid_dir
        #print cmd
        os.system(cmd)

    # Initialise grid
    grid = initGrid( grid_infile, grid_step)
    
    # Get list of cell files
    cell_files = getCellData(cell_dir, cell_interval)

    # Keep track of mean forces/torques/energies and orientational numbers for each water
    water_props = {}
    
    # Update grid with cell data

    print "#Frame       -TDS_vib -TDS_lib -TDS_ori  -TDSW      DHW    DHX    DG (kcal/mol) Forces  (N x m x 1e-10)  Torques (N x m x 1e-20)   MeanWatNrg    OriN   Density (g/L)"

    framecount = 0

    for cell_file in cell_files:
        cell = loadCell( cell_file )
        # sys.exit(-1)
        updatePropsCumulative( water_props, cell)
        updateGridPropsCumulative( grid, cell)
        #sys.exit(-1)
        framecount += 1

        if (framecount % frequencyupdate) == 0:
             framenum = framecount + cell_interval[0]
             DS_vib, DS_lib, DS_ori, DSW, DHSOL, DHSLV, DGW, meanF, meanT, meannrg, meanori, density = calcAverageProps( water_props )
             outstr = "%12d %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f"  % ( framenum, DS_vib, DS_lib, DS_ori, DSW, DHSLV, DHSOL, DGW)
             outstr += "%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f " % ( meanF[0], meanF[1], meanF[2], meanT[0], meanT[1], meanT[2] )
             outstr += "%8.5f %8.5f %8.5f" % (meannrg, meanori, density)
             print outstr
             #sys.exit(-1)
    #sys.exit(-1)
    outputGrid( grid, grid_dir)
    outputWeightedGrid( grid, grid_dir)
    outputGridParams( grid, grid_dir )

