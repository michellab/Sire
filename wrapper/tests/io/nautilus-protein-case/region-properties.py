#!/usr/bin/python
#
# Evaluate thermodynamic properties over a region of the grid
#
import os,sys
from math import log, exp
#
#
# Usage 
# ./script.py my.region grid.forces 
#
# Output
# 
# my_properties.dat
# my_region.pdb

temperature = 298 # Kelvin
waterModel = "TIP4PEW-SireOpenMM"

###################################
# Constants #
waterProps = { "TIP4PEW-SireOpenMM": { "BULKORI": 3.305,
                                       "BULKFORCES": [ 1.587, 1.735, 1.334 ],
                                       "BULKTORQUES": [ 1.061, 1.194, 1.453 ],
                                       "BULKNRG": -11.025,
                                       "BULKDENSITY":  995.4
                                     }
               }

BULKORI = waterProps[waterModel]['BULKORI']
BULKFORCES = waterProps[waterModel]['BULKFORCES']
BULKTORQUES = waterProps[waterModel]['BULKTORQUES']
BULKNRG = waterProps[waterModel]['BULKNRG']
BULKDENSITY = waterProps[waterModel]['BULKDENSITY']

T = temperature
SMALL = 0.000001
MININF = -99999999.0
NAVOGADRO = 6.022141e23
MASS_WATER = 18.01528 # g/mol
KCAL_TO_J = 4184
KB = 0.001987206500956023

#########################################

def loadGrid( gridforces ):
   grid = {}
   
   stream = open(gridforces,'r')
   buffer = stream.readlines()
   stream.close()

   # the grid dimensions/origin/step
   nx, ny, nz, ox, oy, oz, step = buffer[1].split()

   grid['params'] = [ int(nx), int(ny), int(nz), float(ox), float(oy), float(oz), float(step) ]
   grid['points'] = {}
  
   for line in buffer[2:]:
       if line.startswith("#"):
           continue
       elems = line.split()
       index = int( elems[0] )
       x = int( elems[1] )
       y = int( elems[2] )
       z = int( elems[3] )
       fx = float( elems[4] )
       fy = float( elems[5] )
       fz = float( elems[6] )
       tx = float( elems[7] )
       ty = float( elems[8] )
       tz = float( elems[9] )
       ori = float( elems[10] )
       density = float( elems[11] )
       volume = float( elems[12] )
       nrg = float( elems[13] )
       solnrg = float( elems[14] )

       grid['points'][index] = {}
       grid['points'][index]['coords'] = (x, y, z)
       grid['points'][index]['forces'] = (fx, fy, fz)
       grid['points'][index]['torques'] = (tx, ty, tz) 
       grid['points'][index]['ori'] = ori
       grid['points'][index]['density'] = density
       grid['points'][index]['volume'] = volume 
       grid['points'][index]['nrg'] = nrg
       grid['points'][index]['solnrg'] = solnrg

   return grid

def loadRegion( regionfile ):

   region = []
 
   stream = open(regionfile,'r')
   buffer = stream.readlines()
   stream.close()
    
   for line in buffer:
      if line.startswith("#"):
          continue
      elems = line.split()
      for value in elems:
          region.append( int(value) )       

   return region
 
def getProperties( region, grid ):
    #
    # Evaluate thermodynamic properties over the grid regon
    #
    #

    region_DHW = 0.0
    region_DHX = 0.0
    density = 0.0
    region_f = [ 0.0, 0.0, 0.0]
    region_t = [ 0.0, 0.0, 0.0]
    region_ori = 0.0
    region_volume = 0.0
    region_nw = 0.0

    for p in region:
        if not grid.has_key(p):  
            print "Problem ! Point %s in your region file is not present in the grid file !. Abort. " 
            sys.exit(-1)
        density = grid[p]['density']
        volume = grid[p]['volume']
        region_volume += volume
        if density < SMALL:
            # No data on this point, so no contribution
            continue 
        density_prefactor = ( volume * NAVOGADRO * 1e-27) / ( MASS_WATER )
        # Each point contains a variable number of water molecules, so to recover properties over whole volume from grid, we 
        # must weight the contribution of each point by local density (water/A**3)
        density_conv = density * density_prefactor
        #print density, density_conv

        region_DHW += (grid[p]['nrg']-BULKNRG) * density_conv
        region_DHX += (grid[p]['solnrg']) * density_conv
        region_ori += grid[p]['ori'] * density_conv
        for m in range(0,3):
            region_f[m] += grid[p]['forces'][m] * density_conv
            region_t[m] += grid[p]['torques'][m] * density_conv

        #region_volume += volume
        region_nw += density_conv
    #
        
    if region_nw > 0:
       region_ori /= region_nw
       for m in range(0,3):
          region_f[m] /= region_nw
          region_t[m] /= region_nw

    oriratio = 1.0
    vibratio = 1.0
    libratio = 1.0

    if region_nw > 0:
       if region_ori < 1.0:
          region_ori = 1.0
       oriratio *= ( region_ori / BULKORI )
       for m in range(0,3):
          vibratio *= ( BULKFORCES[m] / (region_f[m]) )
          libratio *= ( BULKTORQUES[m] / (region_t[m]) )
    #print vibratio
    region_minTDSvib = -KB*T*region_nw*log( vibratio )
    #print libratio
    region_minTDSlib = -KB*T*region_nw*log( libratio )
    #print oriratio
    region_minTDSori = -KB*T*region_nw*log( oriratio )

    #sys.exit(-1)

    region_DH = region_DHW + region_DHX
    region_minTDS = region_minTDSvib + region_minTDSlib + region_minTDSori
    region_DG = region_DH + region_minTDS
    #print region_nw
    #print region_volume
    if region_nw > 0:
       #print region_nw
       #print region_volume
       #print density_prefactor
       region_density = (region_nw/region_volume)/( ( 1 * NAVOGADRO * 1e-27) / ( MASS_WATER ) )
       #print region_density
    else:
       region_density = 0.0
    #print region_density
    region_density /= BULKDENSITY
    #print region_density

    return (region_DG, region_DH, region_minTDS, region_DHW, region_DHX, region_minTDSvib, region_minTDSlib, region_minTDSori, region_density, region_nw, region_volume)

def writeRegionProperties( regionprops, root ):
      
    region_stats_file = root+"_properties.dat"

    stream = open( region_stats_file, 'w')

    line = "#Excess Thermodynamic properties of water within region\n"
    stream.write(line)
    line = """ 
Delta_G = %8.4f kcal/mol 
Delta_H = %8.4f kcal/mol
-T*Delta_S = %8.4f kcal/mol
Delta_Hw = %8.4f kcal/mol
Delta_Hx = %8.4f kcal/mol
-T*Delta_Svib = %8.4f kcal/mol
-T*Delta_Slib = %8.4f kcal/mol
-T*Delta_Sori = %8.4f kcal/mol
Relative_Density = %8.4f 
Average_number_of_water_molecules = %8.4f 
Volume = %8.4f A^3
""" % ( regionprops[0], regionprops[1], regionprops[2], regionprops[3], regionprops[4], regionprops[5], regionprops[6], regionprops[7], regionprops[8], regionprops[9], regionprops[10] )

    stream.write(line)

    stream.close()

def writeRegionPDB( region, grid,root ):

    region_pdb_file = root+".pdb"

    stream = open( region_pdb_file, 'w')

    originx = grid['params'][3]
    originy = grid['params'][4]
    originz = grid['params'][5]
    step = grid['params'][6]    

    for p in region:
        if not grid['points'].has_key(p):
            print "Problem ! Point %s in your region file is not present in the grid file !. Abort. "
            sys.exit(-1)
        coords = grid['points'][p]['coords']
        cartx = originx + coords[0] * step
        carty = originy + coords[1] * step
        cartz = originz + coords[2] * step
        density = grid['points'][p]['density']
        #volume = grid[p]['volume']
        #density_prefactor = ( volume * NAVOGADRO * 1e-27) / ( MASS_WATER )
        #density_conv = density * density_prefactor

        line='ATOM%7s%3s   SIT %5d    %8.3f%8.3f%8.3f%6.2f%6.2f%12s\n' % (p, 'C', 1, cartx, carty, cartz, 1.00, density/BULKDENSITY, 'C')
        #line='ATOM%7s%2s%24.2f%8.2f%8.2f%5s %8.4f %10.4f%12s\n' % (p, 'C', cartx, carty, cartz, '1.00', 0.00, 0.000, 'H')
        stream.write(line)
    stream.close()

if __name__ == "__main__":
    try:
        regionfile = sys.argv[1]
        gridforces = sys.argv[2]
    except IndexError:
        print "Usage is ./script.py my.region grid.forces"
        sys.exit(-1)
    name = os.path.split(regionfile)[-1]
    root = name.replace(".region","")
    #print name, root
    #sys.exit(-1)
    # Load the grid 
    grid = loadGrid( gridforces )
    # Load the region
    region = loadRegion( regionfile )
    # Get region properties
    region_props = getProperties( region, grid['points'] )
    # Write region properties  
    writeRegionProperties( region_props, root )
    # Write a pdb to visualize the region
    writeRegionPDB( region, grid, root )
