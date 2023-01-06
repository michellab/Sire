# Sire MoleculeParser System

The Sire MoleculeParser system has been built to simplify the
user experience of loading and saving molecules. At a top level,
the user can load a SireSystem::System of molecules, and then
save them back to the same file format by simply typing;

```python
from Sire.IO import MoleculeParser

system = MoleculeParser.read( filenames... )

new_names = MoleculeParser.save( system, new_filename )
```

where `filenames...` is the individual file to be parsed,
or the list of files that are needed to build the system,
and `new_filename` is the root of the filename to be used
for the output data (with this root extended if multiple files
need to be written). The names of the file(s) that are written
are returned in `new_names`.

For example, to read and write a PDB file, the user would use

```python
from Sire.IO import MoleculeParser

system = MoleculeParser.read( "file.pdb" )

new_pdb = MoleculeParser.save( "new_file" )

print("New PDB filename is %s" % new_pdb)

```

while to read an Amber topology/coordinate file pair, the user
would type

```python
from Sire.IO import MoleculeParser

system = MoleculeParser.read( "file.top", "file.crd" )

# or system = MoleculeParser.read( ["file.top","file.crd"] )

new_files = MoleculeParser.save( "new_file" )

print("New top/crd files are called %s" % str(new_files))
```

## How does this work?

Sire has a number of `MoleculeParser` objects that are all derived
from the `MoleculeParser` base class. The call to `MoleculeParser.read`
will try all of the registered parsers to see if any are able to
read the file(s). The parsers that are successful report a score showing
how well they read the file(s). The parser that reported the highest
score for each file parsed is taken forward as the correct parser.

Once the files have been parsed, they must next be converted into a
`SireSystem::System` object. How this works depends on whether the
parser is a "lead" parser, or a follower.

If `MoleculeParser.isTopology()` is true, then the parser is a lead parser,
and so has the responsibility of leading the whole conversion operation.
The parser creates the `SireSystem::System` and populates it with the
molecules that have been parsed from the file. It then delegates control
to all of the non-lead parsers, who supplement the information in each
molecule and in the System with information that they have parsed from
their respective files. Note, currently Sire only supports creating
systems from a set of parsers in which only a single member of the set
is a lead parser. This requirement may change in future versions.

When a `SireSystem::System` has been loaded, it is tagged with the
IDs of the parsers that were used to create that system (using the
`fileformat` property). The call to `MoleculeParser.save()` will
look for this property and will choose the same parsers to write
out the file. Where more than one parser is used, it will write
multiple files, using the root filename given and the default extension
of the parser to create the output filename. The names of the written
file(s) are returned to the user. If the user wishes, they can
specify the exact parser they want either by specifying the extension
of the file, or by passing the fileformat using a property map,
e.g.

```python
filenames = MoleculeParser.save( system, "output.pdb" )
```

or

```python
filenames = MoleculeParser.save( system, {"fileformat":AmberPrm} )
```

## Three step process for reading

The full procedure for reading information from the file and converting it
into a `SireSystem::System` is a three step process:

1. The file is parsed. This involves reading the data from the file and
   converting it into an easy to access form, e.g. arrays of numbers,
   data indexed by keys etc. The purpose is to parse the data into an
   easy to access datastructure that can be used in the next stage.

2. The parsed data is queried to extract all relevant information.
   The information is extracted in terms of the molecular/system data
   that is provided by the file. This means converting the raw string
   or numeric data from the file into information that has chemical
   meaning e.g. arrays of charges, names of atoms, parameters for bonds
   etc.

3. The relevant chemical information is then used to create `SireMol::Molecule`
   objects, and to populate them with properties that correspond to
   the chemical information extracted from the file. These molecules
   are then assembled into a `SireSystem::System`, which itself
   is populated with properties that represent system information,
   and is then returned to the user.

## Three step process for writing

The same process is used in reverse to write data back to a file.

1. The `SireSystem::System` is queried to extract all of the chemical
   information that is needed for a particular file format, with
   that information extracted into, e.g. arrays of charges, bond parameters
   etc.

2. The chemical information is then converted back into the raw arrays
   of text and numeric data that is needed to go back into the file

3. The raw arrays of data are then back-parsed to write a file of
   the specified name. For text files, you need to write an array
   of `QString` objects, which are then written directly by
   `MoleculeParser` to a file. For binary files you need to handle
   writing yourself, e.g. by using interfaces such as the
   `SireIO::NetCDFFile` interface to write a binary NetCDF file.

## Example three step process for AmberPrm

The `SireIO::AmberPrm` class is used to parse Amber prmtop topology
files. The three steps of the parse process are handled by the
following classes:

1. `SireIO::AmberPrm` is responsible for converting the raw text
   of the file into the arrays of text and numeric information
   present from the file. This provides a nice API that lets
   you query this data, e.g. `AmberPrm.floatData(flag)` will
   return the array of floating point data associated with a particular
   Amber prmtop flag, while `AmberPrm.nAtoms()` will return the number
   of atoms present in the file, and `AmberPrm.nBonds()` will
   return the number of bonds in the file.

2. `SireMM::AmberParams` is responsible for holding all of the
   chemical information that is contained within an amber topology
   file for an individual molecule. This provides a nice API for
   easy access to this raw chemical data, e.g. `AmberParams::charges()`
   returns the raw charge information, `AmberParams::bonds()`
   returns the raw bond information. The `SireMM::AmberParams` class
   comes with associated classes, e.g. `SireMM::AmberBond`. These
   associated classes are used to provide a nice API to access
   aspects of the information, e.g. `SireMM::AmberBond` holds
   the `k` and `theta0` data for an individual amber bond, while
   `SireMM::AmberDihedral` holds the series of `SireMM::AmberDihTerm`
   dihedral terms, which contain the phase, force constant etc.
   for each cosine term in an Amber dihedral.

3. `SireMol::Molecule` and its associated classes represent
    the molecule within Sire. The molecular information is
    held within the `SireMol::Molecule` using the property
    system, e.g. the charges are stored in the `charge` property.

## Detail for reading and writing a file using AmberPrm

The process for writing a file involves extracting the information
from each `SireMol::Molecule` to obtain the set of `SireMM::AmberParams`
objects. These are then sent to `SireIO::AmberPrm` to turn back
into raw arrays of numbers and text, which are then converted
to text and written using `MoleculeParser::writeToFile()`.

The process for reading a file involves `MoleculeParser` reading
in an array of textlines, which are then converted into arrays of numbers
and textual information in `AmberPrm`. This then consructs the `System`
by creating each `Molecule`, via extracting the `AmberParams` from the
raw text/number arrays for each molecule, and then combining it
with a `SireMol::MoleculeInfo` to create the full `SireMol::Molecule`.

Relevant functions to look at are `AmberPrm::getMolecule`, which
creates the `ith` molecule from the file via obtaining the `ith`
`MoleculeInfo` object using `AmberPrm::getMolStructure` and the
`ith` set of `AmberParams` using the `AmberPrm::getAmberParams`
function.

The `AmberParams` conversion functions, e.g. `AmberParams.bondFunctions`
do the work of converting from the Amber-format data in `AmberParams`
to the Sire format data in `SireMol::Molecule`. These are called from
the `AmberPrm::getMolecule()` function, e.g.

```c++
    mol.setProperty(map["ambertype"], amber_params.amberTypes());
    mol.setProperty(map["connectivity"], amber_params.connectivity());
    mol.setProperty(map["bond"],
                    amber_params.bondFunctions(InternalPotential::symbols().bond().r()));
    mol.setProperty(map["angle"],
                    amber_params.angleFunctions(InternalPotential::symbols().angle().theta()));
```

are setting the `ambertype`, `connectivity`, `bond` and `angle` properties
of the molecule by asing `amber_params` to convert them from the data held
within itself using the `AmberParams.amberTypes()`, `AmberParams.bondFunctions()`
functions etc.

These in turn call any sub-class conversion functions, e.g. `AmberBond.toExpression()`
will convert the Amber `k` and `theta0` values into a `SireCAS::Expression`
that represents the equation for the bond potential.

Note that, like all properties, we use `map` to allow the user to provide
an alternative name for that property, e.g. `map["bond"]` will return
the user's chosen name for the bond property, or `bond` by default. Always
make sure to use `map` so that a user's property naming system is respected.

The classes are used in reverse when writing a file. So an `AmberBond`
can be constructed from a `SireCAS::Expression`, which will try to interpret
that expression to extract the `k` and `theta0` values for Amber.
Similarly, an `AmberDihedral` will interpret a `SireCAS::Expression` to
extract the series of `AmberDihTerms` that make up the cosine terms
of the Amber dihedral. These are packaged with charges, LJ parameters
etc. to create an `AmberParams` object for each molecule, which is passed
to `AmberPrm` to assemble into arrays of numbers and text that
represent the parsed Amber prm file. This is then converted back to
text and written. Important functions to look at are the `AmberBond`,
`AmberAngle` and `AngleDihedral` constructors, that show how to
extract information from a `SireCAS::Expression`, and the static
functions in `amberprm.cpp` called `getBondData`, `getAngleData` and
`getDihedralData` that show how the parsed bond, angle and dihedral
data is extracted from arrays of `AmberParams` objects from within
the `toLines` function to create the whole set of lines used to
write the file.

## Curiosities

When you look at many of the parsers you may have a couple of questions;

1. Why do we use three levels of parsing? Why not just parse straight
   from `SireMol::Molecule` to the file?

2. Why does the parser when writing a file immediately re-read that file
   in its constructor?

The answer to (1) is because this is a clean design that separates the
three distinct parts of parsing: (i) extracting structured data from a file,
(ii) formatting that structured data into chemical information that is present
within the file, and (iii) convert that chemical information into a form that is
usable within the program.

By using this design, `AmberParams` is separated from `AmberPrm`, meaning
that `AmberParams` can be re-used with any file format parser that uses
the same chemical information that it contains.

The answer to (2) is because, when we write a file, we should immediately
prove that we can re-read it ourselves. By re-reading the file we ensure
that the data written to disk is consistent with the information that
is contained within the `AmberPrm` object. Any errors when re-reading
are immediately caught, ensuring that bugs are found and the 'writer'
part of the parser is always consistent with its 'reader' part.
