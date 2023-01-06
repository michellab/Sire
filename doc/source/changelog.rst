=========
Changelog
=========

There have been four main phases of development of :mod:`sire`
since the project started in 2005.

GitHub (OpenBioSim): January 2023 - current
-------------------------------------------

Development was migrated into the
`OpenBioSim <https://github.com/openbiosim>`__
organisation on `GitHub <https://github.com/openbiosim/sire>`__.

`2023.1.0 <https://github.com/openbiosim/sire/releases/tag/2023.1.0>`__ - January 2023
--------------------------------------------------------------------------------------

* Initial release of the OpenBioSim version of sire. The code has been completely
  refurbished using a tutorial-driven development process.

GitHub (michellab): June 22nd 2015 - January 2023
-------------------------------------------------

Thanks to `@ppxasjsm <https://github.com/ppxasjsm>`__ and
`@jmichel80 <https://github.com/jmichel80>`__ development
was migrated into the `michellab <https://github.com/michellab>`__
organisation on `GitHub <https://github.com/michellab/sire>`__.

This comprised 2495 commits, from developers
`@lohedges <https://github.com/lohedges>`__,
`@chryswoods <https://github.com/chryswoods>`__,
`@ppxasjsm <https://github.com/ppxasjsm>`__,
`@halx <https://github.com/halx>`__,
`@jmichel80 <https://github.com/jmichel80>`__,
`@ptosco <https://github.com/ptosco>`__,
`@SofiaBariami <https://github.com/SofiaBariami>`__,
`@fjclark <https://github.com/fjclark>`__,
`@Steboss <https://github.com/Steboss>`__,
`@nigel-palmer <https://github.com/nigel-palmer>`__,
`@msuruzon <https://github.com/msuruzhon>`__ and
`@kexul <https://github.com/kexul>`__.

Here is the changelog for this stage of development.

..

    [2023.0.3] January 2023: Added the beginnings of a new sphinx-based website
               (in the `doc` folder), which includes the sire API documentation.
               Cleaned up the new sire API via use of `__all__` in all of the
               new modules. The public API is very limited at the moment, but
               will grow as we port in more classes.  However, the aim is that
               users will mostly not create classes directly, but will instead
               implicitly create them as they load molecular systems and call
               functions on those systems. Added a tutorial to this website
               that will be used to demonstrate and teach the new sire API.
               Fully updated the search functionality, making it more robust,
               more consistent and more powerful. Added a detailed guide on the
               search grammar to the new website. Added a set of Cursor classes
               for editing, and made these work consistently with most of the
               property types. Getting and setting properties should now be
               easier, with auto-wrapping and expanding of properties. Made
               the AtomProperty classes behave more like standard python
               containers.  This makes them easier to work with, and is the
               first step to hiding them completely (they will eventually be
               auto-converted to/from standard Python containers or NumPy
               arrays. Added `apply` and `apply_reduce` functions that let you
               map functions across all objects in a molecular container.
               Cleaned up the handling of units - now everything maps into
               GeneralUnit and GeneralUnitProperty, which are auto-converted
               when exposed to Python. Added Python wrapping and
               monkey-patching to sire.maths.Vector so that it has length units.
               Improved the printing of units to the screen (using the correct
               unicode). Added functions that empower the userto choose their
               own default units, e.g. changing angstroms to picometers, or
               switching to full SI units. This only impacts the Python layer
               when rendering the unit, or auto-converting numbers to units,
               so does not break or change the C++ layer. Any view can now be
               assigned a GeneralUnit property. Added Bond, Angle, Dihedral,
               Improper and their related molecule view container classes (e.g.
               SelectorBond, SelectorMBond etc). This allows you to have
               molecule views that represent bonds, angles and dihedrals (or
               collections of these). Added measurement functions so that you
               can easily get their lengths or sizes. Added `.energy()` to let
               you calculate energies of views (or views with views). This uses
               the parameters / forcefield loaded with the molecule(s). You can
               get energies of any views of sub-views. Also created an proper
               return type for energies that embeds the energy components.
               Now `view.energy().components()` works as you would expect.
               Added `.energies()` to molecule containers so that you can get
               the energies of each view in the container. Added support for
               progress bars using Rich so that the user has an indication of
               progress. Added initial support for trajectories. Reworked the
               molecular parser so that multiple "frame" types files will load
               multiple frames of a trajectory (e.g. so that a trajectory can
               be loaded from multiple PDB files, or from multiple DCD or traj
               files). Added a TrajectoryIterator class that lets you easily
               iterate over and query trajectories. Fully documented this in
               the tutorial. You can now do cool things like measure bonds over
               trajectories, or evaluate energies. Added a `.view()` function
               based on NGLView that lets you easily see any molecule view (or
               collection of molecule views). Added a `save_to_string` function
               that writes a text-based molecule file to an in-memory string
               rather than a file (so that you don't have to use temporary
               files with NGLView). Added support for viewing trajectories, so
               that trajectories that are loaded in sire are also playable in
               NGLView. Added movement functions to the Cursor classes so that
               you can more easily move molecules (or molecule views).
               Documented this in the tutorial. Re-worked the way PropertyMap is
               passed via Python. Now have a sire.base.create_map function that
               can create a PropertyMap from anything that is passed. This has
               some examples in its documentation that show how is can be used.
               Made sure that all of the new functionality can use PropertyMap
               and uses `create_map` to support function calls like
               `cursor.translate( (1,2,3), map={"coordinates":"coords2"} )`.
               Speaking of which, also updated `sire.maths.Vector` adding in
               functions that allow auto-conversion of list-like python objects
               to `sire.maths.Vector`. It should almost be the case that a user
               will not have to use this class directly themselves, as things
               should just auto-convert. Added support for creating Vectors
               from plain numbers or length units, using the default length
               unit if plain numbers are used. Fixed lots of bugs and expanded
               the unit test suite to test the above functionality. Removed
               lots of unnecessary files. Moved some files into the website
               docs so that there is a single source of truth. Began the process
               of updating paths and links to point to the new locations in
               OpenBioSim. Fixed CI build issues on Windows by building in the
               right directory. Updated the pythonizing framework so that we
               only pythonize the C++ layer, and avoid the circular dependencies
               that were causing random import errors (particularly on Windows).

    [2023.0.2] December 2023: Fix multiple distance restraint bug in SOMD
               (@fjclark). Add support for PME FEP with SOMD and fix
               associated bugs (@halx, @jmichel80). Fix CI issues so that
               PRs use the correct URL when triggered by external forks.
               Exclude dummy atoms when repartitioning hydrogen masses.
               Deprecate py37.

    [2023.0.1] November 2023: Improve handling of HETATM and TER records in
               PDB files. Fix SOMD selection issues following update to the
               2023 API. Fix writing of steps to SOMD simfile.dat (@fjclark).
               Throw exception when CHAMBER format AMBER topology files are
               detected. Expose toVector() method for the velocity property.
               Match against inverted dihedral records of for A-B-C-A when
               building GROMACS topologies. Fixed calling of static Py++
               functions. Build against conda-forge AmberTools and GROMACS
               packages as host requirements, allowing users to create
               BioSimSpace environments with or without these dependencies
               installed. Added the ability to search on whether or not a
               property exists.  Make sure searches are returned in MolIdx
               order. Ensure Sire is built against packages with the "dev"
               label.

    [2023.0.0] July 2023 - Updated Sire's API to a more pythonic style.
               Module names are in lower case, e.g. `import Sire` becomes
               `import sire`, or `import sire as sr`. Functions are in
               underscore_case. This change is not backwards compatible. To
               support old code, a `sire.use_old_api()` function has been added.
               New functions have been added that make it easier to load
               and save molecules. These can load from URLs. Tests have been
               updated to pytest and now load input data from the sire website.
               The search system has been overhauled, optimised and updated.
               This is described in the new tutorials that are in the process
               of being written in the `doc` directory. This also contains
               the new sphinx website. The `CMakeLists.txt` files and build
               system have been completely reworked. These now use more
               pythonic `setup.py` scripts. These have been updated to fully
               support MacOS M1 and Windows. The conda recipe has been
               updated to use these scripts. Conda packages are now built
               and supported across Linux, MacOS and Windows.

    [2022.3.0] June 2022 - Added support for parsing SDF files (@chryswoods).
               Move conda build process to Miniforge and mambdabuild (boa) to
               avoid timeouts and memory issues. Update GroTop parser to ensure
               new atom types are created when names match but parameters
               differ. Added additional BioSimSpace wrapper to update
               coordinates and velocities in a system, without first requiring
               that it is modified to have unique atom and residue numbers.
               Use -Oz compiler flag rather than -Os for compiling Python
               wrappers to avoid "illegal hardware instruction" error with
               Clang 14 on macOS x86_64. Fixed issue reconstructing triclinic
               box objects from a binary data stream. Added missing streaming
               operators to Sire.Unit.GeneralUnit.

    [2022.2.0] March 2022 - Fixed formatting of SOLVENT_POINTERS flag in
               AmberPrm7 parser. Removed duplicate definition of sigma_av
               in OpenMMFreEnergySt.cpp. Fixed SOMD issues related to
               assumption that perturbable molecule always has MolIdx(1)
               (@fjclark). Fixed wrappers and added significant performance
               enhancements to the SireIO::updateCoordinatesAndVelocities
               function. This significantly (200x) speeds up the remapping
               of coordinates/velocities from SOMD trajectory frames, which
               was a bottleneck for large protein-ligand simulations within
               BioSimSpace. Disabled GSL error handling to avoid a potential
               segmentation fault within a singular value decomposition
               routine called by SireMaths::align.

    [2022.1.0] Jan 2022 - Fixed counting of protons to account for dummy atoms
               when swapping water topology and ensure that original molecular
               properties are preserved. Added a fallback to the BGFS solver
               to improve robustness of FEP analysis (@kexul). Fixed a bug
               that caused distance restraints to be skipped if the ligand
               wasn't the first molecule in the FEP topology (@jmichel80,
               @fjclark). Improved atomic element inference in AMBER parsers.
               Update Sire build to latest versions of dependencies on macOS
               and Linux. This required substantial mini-changes across the
               entire codebase due to changes in APIs and deprecations. This
               includes moving away from qAlgorithm, using the new Qt
               container constructors, moving to OneAPI, switching to the
               conda-forge OpenMM and switching to the new C++ ABI (@chryswoods).
               Simplified Sire wrapper generation using a minimal Docker
               container with the latest Py++ (@chryswoods). Add support for
               native Python pickling of Sire objects (@chryswoods.) Switch
               to GitHub actions for CI. This uses a conda-forge compliant
               conda build, with packages then uploaded to the Anaconda cloud.

    [2021.1.0] Aug 2021 - Added support for multiple combining rules in SOMD
               (@SofiaBariami). Added support for triclinic simulation boxes.
               Convert Ryckaert-Bellememans form dihedral functions from
               GROMACS to Fourier series to allow conversion to AMBER format.
               Updated search functionality to enable searching for objects
               within an arbitrary distance of a point. Fixed PDB2 parser bug
               to ensure that residue names are fixed width. Ensure that
               NUMEXTRA pointer is written so that AMBER topology files can
               be read by tools such as ParmEd. Write NATYP pointer and
               correct number of SOLTY flags. Even though these aren't used,
               incorrect values break external tools, e.g. ParmEd. Added
               support for AMBER TIP5P water topology conversion. Correctly
               flag OPLS style force fields when creating MMDetail object so
               that users can reconstruct OPLS systems written to AMBER format.
               Made build Python 3.8+ compliant. (Python libraries are now
               ABI compatible.) Switched to using std::atomic since
               tbb/atomic.h is now deprecated. Switched to using HTTPS for
               sending analytics. Updated build to be able to link against
               conda version of libcpuid. Added support for generating PDB
               CONECT records from a Sire.Mol.Connectivity object. Fixed
               issue with PMEMD skipping torsions with zero periodicity.
               Fixed random number seeding bug in somd-freenrg, which
               resulted in the OpenMM generator being seeded with the same
               seed for each cycle of the simulation.

    [2020.1.0] July 2020 - Fixed bug in WaterView program to ensure that a
               molecule is extracted from the returned list. Stable sorting
               of dihedrals and other potential terms to allow reproducible
               writing of input files for SOMD (@ptosco). Updated the
               FreeEnergyAnalysis script to support different versions of the
               pymbar API. Significant performance improvement to the GroTop
               parser by looping over cut-groups during non-bonded matrix
               evaluation. Updated Miniconda and conda dependencies to latest
               cross-compatible versions. Fixed minor copiler and runtime
               issues (@nigel-cresset).

    [2019.3.0] November 2019 - Added functionality to restrict the search space
               when finding paths between atoms or searching for rings. Fixed
               performance issue in GroTop parser caused by an N^2 loop over
               atoms when searching the intrascale matrix. We now loop over
               cut-groups, which is far more efficient. Fixed issues with
               Python wrapper generation caused by issues with missing define
               symbols and a bug in the scanheaders.py script.

    [2019.2.1] October 2019 - Updated the Conda recipe to pin the dependencies
               of dependencies that are used at run time since Conda doesn't
               automatically do this for you. Added instructions detailing the
               Azure Pipeline build process and how to create a new release.

    [2019.2.0] September 2019 - Updated the Gromacs topology writer to support
               perturbable molecules containing a variable number of bonds.
               Created a Docker container for building wrappers and updated
               to using CastXML. Added support for running background
               processes on Windows (@ptosco). Updated SOMD Python wrapper
               to write restart files every cycle to simplify system monitoring
               in BioSimSpace. Fixed macOS build issue by not linking against
               libpython. Made sure that Conda dependencies are pinned
               correctly to avoid compatibility issues. Fixed bug that
               prevented upload statistics being sent and added support for
               tracking BioSimSpace usage.

    [2019.1.0] May 2019 - Updates to the Gromacs topology writer to support
               free energy perturbation simulations. The MCS matching
               functionality has been extended to allow matches between heavy
               and light atoms, and the ability to return all current matches,
               rather than just the most recent. Temporarily disabled
               parallelisation in the AmberPrm parser to avoid a threading
               issue. Switched to using Azure Pipelines for continuous
               integration to enable a fully automated build, testing, and
               deployment pipeline. In addition, we finally have created
               a Sire Conda package to simplify the installation and update
               process.

    [2018.2.0] July 2018 - Improvements to the Gromacs topology reader/writer,
               addition of code to improve matching of atoms in proteins,
               fixing compile issues on modern Ubuntu, bugfixes for crashes
               in the AmberPrm reader, added in text-based searching for
               atoms, residues etc. from Systems, MoleculeGroups, Molecules,
               etc. based on boost::spirit, updated boost to latest version,
               bugfixes for quantomm infinite rotation bug for ions,
               general bugfixes.

    [2018.1.1] May 2018 - Small bug fixes to allow single-atom solutes
               and also to fix small issues with some parsers for BioSimSpace

    [2018.1.0] March 2018 - Signficantly improved the Gromacs and Charmm
               parser and  fixed bugs. Can now write with both :-). Fixed
               compilation on Windows 7 and above. Small changes to
               the API of AtomProperties to make them easier to work
               with from Python. Added a script to automatically color
               code swap-based free energies. Fixed localisation
               problems for the PDB writer. Improved mbar analysis
               code. Added code to track forcefield data of a molecue.
               Added code to better manage processes, including redirection
               of standard output and error.

    [2017.3.0] December 2017 - Added new PDB (PDB2), Mol2, Gro87, and CharmmPSF
               parsers, as well as a GroTop parser, all part of the new
               MoleculeParser framework. Updated all of the swap-based
               methods to use this.

               Removed the ViewsOfMol Python wrapper and now have the
               code automatically return the correct python object for
               the molecule (or part of molecule) that is returned from
               the system. This makes simple scripts easier to write.

               General bugfixes and optimisations, including fixing
               bugs with the way that PropertyMap worked, cleaning
               up to/from converters from python objects to automatic
               Property wrappers, and fixing Process so that it can
               redirect to stdout and that isRunning works without the
               user having to call "wait" first!

    [2017.2.0] September 2017 - The MoleculeParser framework has
               been created to support reading and writing of molecules
               in lots of different formats. The first set of formats
               that have been completed are Amber Prmtop, Amber Rst7
               and Amber Rst/Trj. The parsers work in parallel, with
               file formats automatically detected by the parser,
               e.g. system = MoleculeParser.read( "file.prm", "file.rst" )
               will automatically do the right thing.

               Improved automatic compilation on Arch linux.

               Fixed temperature checking and general bugfixing for mbar code.

    [2017.1.0] April 20 2017 - First 2017 release. Included new
            parallel MoleculeParser code for reading molecules,
            and moved fully over to C++-14 style coding for new code.
            Included AVX-512 vectorisation for Intel KNL and can
            now successfully compile and run using GCC 5 and GCC 6,
            as well as Intel 2017 compilers and Clang.

    [2016.3.1] January 9 2017 - Minor patch release that fixes bugs:
        (1) Correctly sets MACOSX_DEPLOYMENT_TARGET to 10.8 so executable works
            on OS X 10.8 (Mountain Lion) and above
        (2) Fixed a parsing bug in Parameter that prevented integer or float
            parameters from being passed to ligandswap, waterswap etc.
        (3) Fixed a small bug in MultiDouble that meant it lost precision when
            swapping individual values
        (4) Fixed a parsing bug in Parameter that meant that windows path names
            were not interpreted correctly
        (5) Fixed the build scripts so that they placed bundled libraries into
            bundled/lib rather than bundled/lib64 (affected SUSE-based distributions)

    [2016.3.0] December 22 2016 - Public release containing full LigandSwap. Uses
     new optimised forcefields for energy calculations, built on top of Intel Threaded
     Building Blocks for parallelisation. New code is significantly faster with better
     scaling.

    [2016.2.0] June 3 2016 - Semi-private release for the CCP-BioSim workshops. Included
     the first version of LigandSwap and general bug fixes

    [2016.1.0] April 29 2016 - Merge of Bristol and Edinburgh codes, moved to miniconda
     and clean packaging system, including OpenMM fully, added in nautilus, somd etc.,
     added optimised forcefields, added a proper unit testing suite.

    [OLD] Updated gradient compuation in openmmfreenergst to finite differece gradients
    Allow the computation of reduced perturbed energies in openmmfreenergst of all computed lambda values
    Separated minimization and equilibration from production run.
    Implemented mass repartitioning for hydrogens atoms to allow for larger integration timesteps
    Added nautilus scripts

Google Code: August 7th 2006 - April 1st 2015
---------------------------------------------

Sire was developed against the subversion repository provided
by Google Code. Here is an
`archive of the repository <https://code.google.com/p/sire>`__.

This comprised 2775 commits, from developers
`@chryswoods <https://github.com/chryswoods>`__,
`@jmichel80 <https://github.com/jmichel80>`__ and
`@nividic73 <mailto:nividic73@googlemail.com>`__.

Local Subversion: February 5th 2005 - July 25th 2006
----------------------------------------------------

Sire was developed against a local subversion repository.
Here is a
`svndump of the original repository <https://siremol.org/f/orig_sire_repository.dump.bz2>`__,
and all of the `commit history <https://siremol.org/f/original_repository_comments.txt>`__.

This comprised 831 commits from developer `@chryswoods <https://github.com/chryswoods>`__.

Sire started as ``ProtoMS 3``, a complete C++ rewrite of
`ProtoMS 2 <https://code.google.com/archive/p/protoms/source/default/commits>`__,
developed originally as a Fortran program
by `@chryswoods <https://github.com/chryswoods>`__ and
`@jmichel80 <https://github.com/jmichel80>`__. ProtoMS has since continued
to be developed by the
`Essex Group <https://www.essexgroup.soton.ac.uk>`__ and is
itself now available as `ProtoMS 3.4 <https://protoms.org>`__.

More detail about the history and parallel development of Sire and
ProtoMS can be `found here <https://www.essexgroup.soton.ac.uk/ProtoMS/FAQ/index.html>`__.
