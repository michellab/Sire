
import struct, time, array, os

from Sire.Mol import *
from Sire.IO import *

#
# Adapted from Peter Eastman's code in OpenMM python API to write a DCD file
#

class DCDFile(object):
    """DCDFile provides methods for creating DCD files.
    
    DCD is a file format for storing simulation trajectories.  It is supported by many programs, such
    as CHARMM, NAMD, and X-PLOR.  Note, however, that different programs produce subtly different
    versions of the format.  This class generates the CHARMM version.  Also note that there is no
    standard byte ordering (big-endian or little-endian) for this format.  This class always generates
    files with little-endian ordering.
    
    To use this class, create a DCDFile object, then call writeModel() once for each model in the file."""
    
    def __init__(self, strfile, group, space, dt, firstStep=0, interval=1):
        """Create a DCD file and write out the header.
        
        Parameters:
         - file (file) A file to write to
         - topology (Topology) The Topology defining the molecular system being written
         - dt (time) The time step used in the trajectory
         - firstStep (int=0) The index of the first step in the trajectory
         - interval (int=1) The frequency (measured in time steps) at which states are written to the trajectory
        """

        file = open(strfile,'wb')

        #PDB().write(group, "%s.pdb" % strfile)

        self._file = file
        self._group = group
        self._space = space
        self._firstStep = firstStep
        self._interval = interval
        self._modelCount = 0
        #if is_quantity(dt):
        #    dt = dt.value_in_unit(picoseconds)
        #dt /= 0.04888821
        dt = dt.value()

        natoms = 0
        molecules = group.molecules()
        molnums = molecules.molNums()
        for molnum in molnums:
            mol = molecules.molecule(molnum)[0]
            nat = mol.nAtoms()
            natoms += nat

        print("There are %s atoms in the group " % natoms)

        #sys.exit(-1)
        boxFlag = 0
        if space.isPeriodic():
            boxFlag = 1
        header = struct.pack(b'<i4c9if', 84, b'C', b'O', b'R', b'D', 0, firstStep, interval, 0, 0, 0, 0, 0, 0, dt)
        header += struct.pack(b'<13i', boxFlag, 0, 0, 0, 0, 0, 0, 0, 0, 24, 84, 164, 2)
        header += struct.pack(b'<80s', b'Created by OpenMM')
        header += struct.pack(b'<80s', bytes('Created '+time.asctime(time.localtime(time.time())),"utf-8"))
        header += struct.pack(b'<4i', 164, 4, natoms, 4)
        file.write( header )
    
    def writeModel(self, group, space):
        """Write out a model to the DCD file.
                
        Parameters:
         - positions (list) The list of atomic positions to write
        """
        #if len(list(self._topology.atoms())) != len(positions):
        #    raise ValueError('The number of positions must match the number of atoms') 
        #if is_quantity(positions):
        #    positions = positions.value_in_unit(nanometers)
        file = self._file
        
        # Update the header.
        
        self._modelCount += 1
        file.seek(8, os.SEEK_SET)
        file.write(struct.pack('<i', self._modelCount))
        file.seek(20, os.SEEK_SET)
        file.write(struct.pack('<i', self._firstStep+self._modelCount*self._interval))
        
        # Write the data.
        
        file.seek(0, os.SEEK_END)
        
        if space.isPeriodic():
            boxSize = space.dimensions()
            file.write(struct.pack('<i6di', 48, boxSize[0], 0, boxSize[1], 0, 0, boxSize[2], 48))

        natoms = 0

        for i in range(0,group.nMolecules()):
            mol = group[MolIdx(i)]
            nat = mol.nAtoms()
            natoms += nat

        length = struct.pack('<i', 4*natoms)

        # To get the positions...
        # Loop over that group
        nmols = group.nMolecules()
        
        coords = []

        #spacedims = space.dimensions()

        #wrapmolcoordinates = False

        #wrapatomcoordinates = False

        # JM 10/14 bugfix change of behavior of QSet in QT5
        molnums = group.molNums()
        molnums.sort()

        for i in range(0,group.nMolecules()):
            #mol = group[MolIdx(i)].molecule()
            mol = group[molnums[i]][0]
            #print (mol)
            molcoords = mol.property("coordinates")

            #if wrapmolcoordinates:
            #    molcog = CenterOfGeometry(mol).point()
            #
            #    wrapdelta = Vector( int( math.floor( molcog.x() / spacedims.x() ) ) ,\
            #                        int( math.floor( molcog.y() / spacedims.y() ) ) ,\
            #                        int( math.floor( molcog.z() / spacedims.z() ) ) )
            #
            #    if ( wrapdelta[0] != 0 or wrapdelta[1] != 0 or wrapdelta[2] != 0):
            #        print("Mol %s wrapdelta %s %s %s " % (molnum.toString(), wrapdelta[0], wrapdelta[1], wrapdelta[2]))
            #        print(spacedims)
            #        print(molcoords.toVector())
            #        wrap = Vector( - wrapdelta[0] * spacedims.x() , - wrapdelta[1] * spacedims.y(), -wrapdelta[2] * spacedims.z() )               
            #        molcoords.translate(wrap)
            #        print(molcoords.toVector())
                    
            
            #molcoords.translate(wrapdelta)
            #coords += molcoords
            coords += molcoords.toVector()
           
            #if wrapatomcoordinates:
            #    molvec = molcoords.toVector()
            #    for atvec in molvec:
            #        wrapdelta = Vector( int( math.floor( atvec.x() / spacedims.x() ) ) ,\
            #                            int( math.floor( atvec.y() / spacedims.y() ) ) ,\
            #                            int( math.floor( atvec.z() / spacedims.z() ) ) )                
            #        if ( wrapdelta[0] != 0 or wrapdelta[1] != 0 or wrapdelta[2] != 0):
            #            wrap = Vector( - wrapdelta[0] * spacedims.x() , - wrapdelta[1] * spacedims.y(), -wrapdelta[2] * spacedims.z() )
            #            atvec = atvec + wrap
            #        coords += atvec

        #print coords
        #print len(coords)

        # Have to study that bit...
        for i in range(3):
            file.write(length)
            data = array.array('f', (x[i] for x in coords))
            data.tofile(file)
            file.write(length)

    def writeBufferedModels(self, group, dimensions):
        """Write out a collection of snapshots to the DCD file.
                
        Parameters:
         - positions (list) The list of atomic positions to write
        """
        #if len(list(self._topology.atoms())) != len(positions):
        #    raise ValueError('The number of positions must match the number of atoms') 
        #if is_quantity(positions):
        #    positions = positions.value_in_unit(nanometers)
        file = self._file

        # Find the number of buffered frames we have by inspecting the first molecule in the group
        # assuming all molecules have same number of buffered coordinates...
        mol = group.first()[0]
        molprops = mol.propertyKeys()
        nbuf = 0
        for molprop in molprops:
            if molprop.startswith("buffered_coord"):
                nbuf += 1
        if nbuf <= 0:
            print("Could not find any buffered coordinates in the passed group ! ")
            return
        
        #
        # Should be more efficient to loop over all mols once
        #
        for x in range(0,nbuf):
            # Update the header
            self._modelCount += 1
            file.seek(8, os.SEEK_SET)
            file.write(struct.pack('<i', self._modelCount))
            file.seek(20, os.SEEK_SET)
            file.write(struct.pack('<i', self._firstStep+self._modelCount*self._interval))
          
            # Write the data.
        
            file.seek(0, os.SEEK_END)
            # Get buffered space...
            boxSize = None
            if ("buffered_space_%s" % x) in dimensions:
                boxSize = dimensions["buffered_space_%s" % x].dimensions()
            #print "buffered_space_%s" % x, boxSize
            if boxSize is not None:
                file.write(struct.pack('<i6di', 48, boxSize[0], 0, boxSize[1], 0, 0, boxSize[2], 48))

            natoms = 0
            

            for i in range(0,group.nMolecules()):
                mol = group[MolIdx(i)][0]
                nat = mol.nAtoms()
                natoms += nat

            length = struct.pack('<i', 4*natoms)

            # To get the positions...
            # Loop over that group
            nmols = group.nMolecules()
        
            coords = []

            # JM 10/14 bugfix change of behavior of QSet in QT5
            molnums = group.molNums()
            molnums.sort()
            for i in range(0,group.nMolecules()):
                #mol = group[MolIdx(i)].molecule()
                mol = group[molnums[i]][0]
                molcoords = mol.property("buffered_coord_%s" % x)

                coords += molcoords.toVector()
           
            # Have to study that bit...
            for i in range(3):
                file.write(length)
                data = array.array('f', (x[i] for x in coords))
                data.tofile(file)
                file.write(length)
            #rewind
            file.seek(0, os.SEEK_SET)

