#!/usr/bin/python

import os,sys, math, time

from MDAnalysis import Universe

top = "new.top"
crd = "min00001.rst7"
dcd = "traj000000001.dcd"
grid = "grid.in"

SLEEP = 30

def usage():
    print "USage is script nworkers local/queue maxframes"
    sys.exit(-1)

if __name__ == "__main__":

    try:
        nworkers = sys.argv[1]
        flag = sys.argv[2]
    except IndexError:
        usage()

    try:
        nworkers = int(nworkers)
    except:
        usage()

    if flag != "local" and flag != "serial14.q" and flag != "serial57.q"  :
        usage()

    maxframes = -1
    try:
        maxframes = sys.argv[3]
    except:
        usage()

    maxframes = int(maxframes)

    # This seems to overload the cluster. trying without...
    #dcd_traj = Universe(top, dcd)
    #dcdframes = len(dcd_traj.trajectory)

    #nframes = dcdframes

    #if maxframes > 0:
    #    if maxframes < nframes:
    #        nframes = maxframes
   
    nframes = maxframes    
    dcdframes = maxframes
 
    print "Will process first %s frames out of %s. We have %s workers" % ( nframes, dcdframes, nworkers)

    chunks = []

    ratio = nframes/float(nworkers)
    
    ratio_integer = int( math.floor(ratio) )
    ratio_remainder = ratio - ratio_integer
    
    # Add one more frame every extra steps
    if ratio_remainder > 0:
        extra = int( round( 1/ratio_remainder) )
    else:
        extra = 0

    start = 0
    end = -1
    for x in range(0, nworkers):
        end += ( ratio_integer  )
        if ( extra > 0 ):
            if ( (x % extra) == 0 ):
                end += 1
        chunks.append( ( start, end ) )
        start = end + 1

    #print chunks

    sleeprange = len(chunks)*SLEEP

    count = 1
    for (start_frame, end_frame) in chunks:
        sub_script = """
#!/bin/bash
#
# Script to run a gpu job
#

# SGE submission options
#$ -q %s             # Select the queue
#$ -l h_rt=48:00:00        # Set 48 hour of wall clock time
#$ -cwd                    # Change to current working directory
#$ -V                      # Export environment variables into script

rand=$[ 1 + $[ RANDOM %% %s ]]

echo sleep $rand
sleep $rand

export OMP_NUM_THREADS=1
date
python compute-cell-properties-water.py %s %s %s %s %s %s
date

""" % (flag, sleeprange, top, crd, dcd, grid, start_frame, end_frame)
        stream = open("cell_chunk_%05d.sh" % count,"w")
        stream.writelines(sub_script)
        stream.close()
        count += 1

    # Time for the OS to write all files...
    time.sleep(2)
    
    count = 1
    for (start_frame, end_frame) in chunks:    
        if flag == "local":
            cmd = "bash cell_chunk_%05d.sh 1> /dev/null 2> /dev/null  &" % count
        else:
            cmd = "qsub cell_chunk_%05d.sh" % count
        print cmd
        os.system(cmd)
        count += 1
        time.sleep(1)
        #sys.exit(-1)
