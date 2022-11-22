===============================================
Part 4 - Measurement, Trajectories and Movement
===============================================

The whole purpose of loading molecules into sire is that you can then
do things with them. Typical things you may want to do could include;

* Measuring distances, angles or torsion angles between atoms, residue or
  molecules.
* Making measurements across multiple frames of a trajectory.
* Measuring energies of molecules or subsets of molecules, and
  across multiple frames of a trajectory.
* Moving atoms, residues or molecules, e.g. by direct translation and
  rotation of parts, or even by stretching bonds or rotating angles,
  and saving these to a trajectory.
* Performing energy scans by calculating the energy for different
  configurations of the molecule (e.g. a dihedral scan).

This chapter will teach you to do all of the above. Since measurement,
trajectories and movement are easier if you can see molecules and
create graphs, you will also learn how to visualise sire objects and data.

.. toctree::
   :maxdepth: 1

   part04/01_measure
   part04/02_trajectory
   part04/03_energies
   part04/04_energy_trajectories
   part04/05_movement
   part04/06_move_internals
