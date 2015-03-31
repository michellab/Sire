# Parameter files
parfile1 /Users/chris/Work/ProtoMS/parameter/amber99.ff
parfile2 /Users/chris/Work/ProtoMS/parameter/solvents.ff
parfile3 /Users/chris/Work/ProtoMS/parameter/gaff.ff
parfile4 ethane2methanol.par
# PDB Files
SIRE_COMPATIBILITY_MODE on
solute1 ethane2methanol.pdb
#solvent1 boxT4P.pdb
#set the output files
streamwarning stdout
streamfatal stdout
streamfatal fatal
streamheader stdout
streaminfo stdout
streamdetail off
streamenergy stdout
streamaccept off
cutoff 10.00
feather 0.5
boundary solvent
nptsim on
pressure 1.0
prefsampling 1 200.0
lambda 0.5 0.0 1.0

chunk1 singlepoint
chunk2 soluteenergy 1
#chunk3 pdb all file=protoms_test.pdb
#chunk4 printparameters

