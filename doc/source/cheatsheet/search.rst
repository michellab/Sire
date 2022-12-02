==================
Search Cheat Sheet
==================

Views
-----

Searches are based on the concept of "views":

* ``atom`` - an individual atom
* ``residue`` - an individual residue. Atoms (optionally) belong to residues.
  Also sometimes shortened to ``res``.
* ``chain`` - an individual chain. Residues (optionally) belong to chains.
* ``segment`` - an individual collection of atoms.
  Also sometimes shortened to ``seg``.
* ``molecule`` - a complete molecule
  Also sometimes shortened to ``mol``.

Identifiers
-----------

Searches then use the different types of identifiers for views:

* ``name`` - the name of the view, e.g. ``atomname CA`` would match atoms
  that are called ``CA``, while ``resname ALA`` matches residues called ``ALA``.
* ``num`` - the number of the view, e.g. ``atomnum 35`` matches atoms with
  number ``35``, while ``molnum 1`` matches molecule number ``1``. Note that
  only atoms, residues and molecules have numbers.
* ``idx`` - the index of the view within its container, e.g.
  ``chainidx 0`` would match the first chain in the molecule, while
  ``segidx -1`` would match the last segment.

The names of these identifiers for each view are shown in this table.

+----------+-----------+---------+----------+
| View     | name      | number  | index    |
+==========+===========+=========+==========+
| atom     | atomname  | atomnum | atomidx  |
+----------+-----------+---------+----------+
| residue  | rename    | resnum  | residx   |
+----------+-----------+---------+----------+
| chain    | chainname |  N/A    | chainidx |
+----------+-----------+---------+----------+
| segment  | segname   |  N/A    | segidx   |
+----------+-----------+---------+----------+
| molecule | molname   | molnum  | molidx   |
+----------+-----------+---------+----------+

Numbers
-------

Numbers used in searches can be:

* Individual, e.g. ``resnum 35`` searches for the residue with number ``35``.
* A range, e.g. ``atomidx 0:10`` searches for atoms with indicies in the
  range from 0 to 9, i.e. the half-open range ``[0-10)``.
* A stepped range, e.g. ``molidx 0:2:10`` would search for molecules with
  indicies 0, 2, 4, 6, and 8, i.e. the half-open range from 0 to 10 in
  steps of 2.
* A reverse range, e.g. ``residx 3:0:-1`` would search for residues with
  indicies 3, 2, or 1, i.e. the half-open range from 3 to 0 in steps of -1
* An unbounded range, e.g. ``molidx 10:`` would search for molecules with
  indicies 10 upwards. Similarly, ``residx 3::-1`` search from 3 downwards,
  in steps of -1, so would match 3, 2, 1, or 0.
* A list of numbers, e.g. ``resnum 12,14,20`` would search for residues
  with numbers 12, 14 or 20.
* Any combination of the above, e.g. ``molidx 1,3,5,10:15,100:`` would search
  for molecules with indicies 1, 3 or 5, or in the half-open range from
  ``[10,15)``, or from 100 upwards.

Strings
-------

Strings used in searches can be:

* Individual, e.g. ``resname ALA`` searches only for residues with name ``ALA``.
* A regular expression, e.g. ``resname /AL*/`` matches any residues whose
  name starts with ``AL``. Similarly, ``resname /HI[ESDP]/`` matches any
  residues whose names are ``HI`` followed by a ``E``, ``S``, ``D`` or ``P``, i.e.
  ``HIE``, ``HIS``, ``HID`` or ``HIP``.
* A case-insensitive regular expression, e.g. ``resname /ala/i`` matches any
  residue whose name (in any case) matches ``ala``.
* A list of names, e.g. ``resname ALA,ARG,ASP`` would match any residue whose
  name was ``ALA`` or ``ARG`` or ``ASP``.
* Any combination of the above, e.g. ``atomname N,/C*/i`` would match atoms
  whose names were ``N`` or names that started with ``C`` or ``c``.

Comparisons
-----------

Numbers (and other values as discussed below) can also be used in comparisons
using the standard comparison operators:

* ``==`` - compare equal, e.g. ``resnum == 3`` matches residues with number 3.
  Note you can compare strings too, e.g. ``atomname == CA`` matches atoms
  called ``CA``.
* ``!=`` - compare not equal, e.g. ``atomnum != 5`` matches all atoms which
  don't have a number 5. You can compare strings too, e.g.
  ``resnam != ALA`` matches all residues with names that are not equal
  to ``ALA``.
* ``>`` - compare greater than, e.g. ``molidx > 5`` matches all molecule
  indicies greater than 5.
* ``>=`` - compare greater or equal to, e.g ``atomidx >= 10`` matches all
  atom indicies greater or equal to 10.
* ``<`` - compare less than, e.g. ``segidx < 2`` matches all segment
  indicies that are less than 2.
* ``<=`` - compare less than or equal to, e.g. ``molnum <= 100`` matches
  all molecule numbers less than or equal to 100.
* ``=~`` compare approximately equal to, e.g. ``atom mass =~ 16`` matches
  all atoms whose mass is approximately 16.

Logical operators
-----------------

Searches can be combined using the logical operators ``and``, ``or`` and
``not``. You can combine as many searches as you want. It is best
to use round brackets to specify order when combining lots of searchers.

For example,

* ``not X`` - return views that don't match ``X``, e.g. ``not atomname CA``
  matches all atoms whose name is not ``CA``.
* ``X and Y`` - return views that match search ``X`` and search ``Y``.

