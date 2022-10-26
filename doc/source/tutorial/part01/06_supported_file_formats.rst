======================
Supported file formats
======================

``sire`` supports reading and writing to many common molecular file formats.
You can print the list of supported formats using

.. code-block:: python

   >>> print(sr.supported_formats())

   ## Parser Gro87 ##
   Supports files: gro
   Gromacs Gro87 structure format files.
   ##################

   ## Parser GroTop ##
   Supports files: top
   Gromacs Topology format files.
   ###################

   ## Parser MOL2 ##
   Supports files: mol2
   Sybyl Mol2 format files.
   #################

   ## Parser PDB ##
   Supports files: pdb
   Protein Data Bank (PDB) format files.
   ################

   ## Parser PRM7 ##
   Supports files: prm7, prm, top7, top, prmtop7, prmtop
   Amber topology/parameter format files supported from Amber 7 upwards.
   #################

   ## Parser PSF ##
   Supports files: psf
   Charmm PSF format files.
   ################

   ## Parser RST ##
   Supports files: rst, crd, trj, traj
   Amber coordinate/velocity binary (netcdf) restart/trajectory files supported since Amber 9, now default since Amber 16.
   ################

   ## Parser RST7 ##
   Supports files: rst7, rst, crd7, crd
   Amber coordinate/velocity text (ascii) restart files supported from Amber 7 upwards.
   #################

   ## Parser SDF ##
   Supports files: sdf, mol
   Structure Data File (SDF) format files.
   ################

   ## Parser SUPPLEMENTARY ##
   Supports files: *
   Files that are supplementary to a lead parser.
   ##########################

Symmetric Input / Output
------------------------

One of our design principles is that molecular file input and output
is symmetrical. This means that ``sire`` can read in and write out the same
amount of information from a file (i.e. it can always read what it writes).

Another design principle is that information should not be lost. As much
as possible, ``sire`` will load and preserve all information it can
read from a molecular file.
