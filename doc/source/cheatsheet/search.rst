==================
Search Cheat Sheet
==================

Views
-----

Searches are based on the concept of "views":

* ``atom`` - an individual atom
* ``bond`` - an individual bond. Pairs of atoms can (optionally) be connected
  together via bonds.
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
| bond     | N/A       | N/A     | N/A      |
+----------+-----------+---------+----------+
| residue  | rename    | resnum  | residx   |
+----------+-----------+---------+----------+
| chain    | chainname |  N/A    | chainidx |
+----------+-----------+---------+----------+
| segment  | segname   |  N/A    | segidx   |
+----------+-----------+---------+----------+
| molecule | molname   | molnum  | molidx   |
+----------+-----------+---------+----------+

Note that ``bond`` views cannot be identified directly by name,
number or index.

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
* A regular expression glob pattern, e.g. ``resname /AL*/`` matches any residues whose
  name starts with ``AL``. Similarly, ``resname /HI[ESDP]/`` matches any
  residues whose names are ``HI`` followed by a ``E``, ``S``, ``D`` or ``P``, i.e.
  ``HIE``, ``HIS``, ``HID`` or ``HIP``. This is implemented using Qt's
  `wildCardToRegularExpression <https://doc.qt.io/qt-5/qregularexpression.html#wildcardToRegularExpression>`__
  function, the syntax of which is `further described here <https://en.wikipedia.org/wiki/Glob_(programming)>`__.
* A case-insensitive regular expression globa pattern, e.g. ``resname /ala/i`` matches any
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
* ``X and Y`` - return views that match search ``X`` and search ``Y``, e.g.
  ``atomname CA and resname ALA`` or ``atomname N and (not resname HIS)``.
* ``X or Y`` - return views that match search ``X`` or search ``Y``, e.g.
  ``atomnum 5 or atomname C`` or ``atomname N and (resname HIS or resnum 5)``.
  Note that passing a list of names or numbers is equivalent to passing
  the names or numbers via ``or``, e.g. ``resname ALA or resname ASP`` is
  the same as ``resname ALA,ASP``.

Returned View Types
-------------------

Every search has a default view return type. This is the type of view
that is returned by that search. Atom-based searches return
:class:`~sire.mol.Atom` views, residue-based searches return
:class:`~sire.mol.Residue` views etc.

When two searches with different return view types are combined, the smallest
view type is returned. So a residue-based search combined with
an atom-based search will return :class:`~sire.mol.Atom` views.
A molecule-based search combined with a chain-based search will return
:class:`~sire.mol.Chain` views.

Segment-based searches introduce complications, because they sit outside
of the Atom < Residue < Chain < Molecule hierarchy. As such, a segment-based
search combined with any smaller view-type search will always return
:class:`~sire.mol.Atom` views.

This table shows the return view type for any combination of two
searches.

+----------+------+-------+---------+---------+----------+----------+
|          | Atom | Bond  | Residue | Chain   | Segment  | Molecule |
+----------+------+-------+---------+---------+----------+----------+
| Atom     | Atom | Atom  | Atom    | Atom    | Atom     | Atom     |
+----------+------+-------+---------+---------+----------+----------+
| Bond     | Atom | Bond  | Bond    | Bond    | Bond     | Bond     |
+----------+------+-------+---------+---------+----------+----------+
| Residue  | Atom | Bond  | Residue | Residue | Atom     | Residue  |
+----------+------+-------+---------+---------+----------+----------+
| Chain    | Atom | Bond  | Residue | Chain   | Atom     | Chain    |
+----------+------+-------+---------+---------+----------+----------+
| Segment  | Atom | Bond  | Atom    | Atom    | Segment  | Segment  |
+----------+------+-------+---------+---------+----------+----------+
| Molecule | Atom | Bond  | Residue | Chain   | Segment  | Molecule |
+----------+------+-------+---------+---------+----------+----------+

Combinations of more than two searches are always decomposed into
pairs of searches, which are combined using the precedence rules
of the grammar, and with a return view type as described in the
above table. To avoid confusion, it is better to indicate
your preferred precedence using round brackets, e.g.
use brackets to choose between ``(atomname CA or atomname C) and resname ALA``
or ``atomname CA or (atomname C and resname ALA)``.

Expansive (`with`) and Contractive (`in`) Searches
--------------------------------------------------

Expansive searches are those which return results that are
at least the same size or larger than the original view being searched. For example, searching
for the molecule from an atom view will give a molecule result. This
is at least the same size (for single-atom molecules), but likely larger
(for multi-atom molecules) than the searched atom.

Contractive searches are those which return results that are
at most the same size or smaller than the original view being searched.
For example, searching for an atom from a molecule view will give an atom
result. This is at least the same size (for single-atom molecules), but
likely smaller (for multi-atom molecules) than the searched molecule.

Whether a search is expansive or contractive depends on the relative
size of the view being searched, and the view that is being returned
as the result. This can be summarised as a table.

TABLE OF EXPANSIVE AND CONTRACTIVE

Another way to think of this is that contractive searches are looking
inside a view, while expansive searches are looking outside a view
(looking for results with the original searched view contained within).

The ``in`` and ``with`` keywords enable you to be explicit about the
contractive or expansive nature of a search.

* ``view_type in X`` - perform a contractive search returning views of type ``view_type``
  within the results of the search ``X``. Examples include; ``atoms in resname ALA``
  (return all of the atom views in residues that have name ``ALA``);
  ``residues in chainidx 0`` (return all of the residue views in the first
  chain); and ``bonds in *`` (return all of the bonds in the current view).

* ``view_type with X`` - perform an expansive search returning views of type ``view_type``
  that contain the results of the search ``X``. Examples include;
  ``residues with atomname CA`` (return all of the residues that contain
  atoms called ``CA``); ``molecules with resname ALA``` (return all of the
  molecules with residues called ``ALA``); and ``bonds with atomname H``

  (return all of the bonds that contain at least one atom called ``H``)

Match All Atoms, Residues, Chains etc.
--------------------------------------

You can match everything, specifying the return view type using the
"match all" keywords:

* ``atoms`` - return all of the atoms in the searched object.
* ``bonds`` - return all of the bonds in the searched object.
* ``residues`` - return all of the residues in the searched object.
* ``chains`` - return all of the chains in the searched object.
* ``segments`` - return all of the segments in the searched object.
* ``molecules`` - return all of the molecules in the searched object.
* ``all`` or ``*`` - everything in the current view being searched.

Whether these are expansive or contractive depends on the view that
is searched, based on the same rules as for normal searches
(e.g. searching for ``residues`` on an atoms view will be expansive,
 as it returns all residues that contain those atoms, while searching
 for ``residues`` on a molecule view will be contractive, as it searches
 for all residues in that molecule)

Advanced ``in`` and ``with`` Searches
-------------------------------------

``with`` and ``in`` can be used for more than just ``view_type in X`` searches.
The general syntax is ``X in Y`` or ``X with Y``.

* ``X in Y`` - search for views that match ``X`` in the results of the views
  that match ``Y``. An example would be ``atomname CA in resname ALA`` (find
  atoms called ``CA`` in residues called ``ALA``)
* ``X with Y`` - search for views that match ``Y`` in the results of the views
  that match ``X``. An example would be ``resname ALA with atomname CA``
  (find residues called ``ALA`` in all of the residues that contain atoms
  called ``CA``)

Understanding this, you can see that ``residues with atomname CA`` is a search
that finds all residues that are with (contain) the results of searching for atoms
with name ``CA``. Similarly ``atoms in resname ALA`` is searching for all atoms
that are in the results of searching for residues with name ``ALA``.
Similarly ``bonds in *`` means search for all bonds that are in the
current view (as ``*`` or ``all`` is returns the view being searched).

Searching for Bonds using ``in``, ``with``, ``from`` and ``to``
---------------------------------------------------------------

Searches involving bonds are more complex as they involve bridging
between (potentially) two views. At the most basic, they always
involve connecting two atoms. But this bond could be entirely within
a residue, or between pairs of residues (and similarly for chains
and segments).

* ``bonds with X`` - an expansive search that finds all bonds that
  contains the result of ``X``. As the only view smaller than a bond
  is an atom, ``X`` can only be a search that returns atom views.
  For example, ``bonds with atomname CA`` would return all bonds
  where at least one of the atoms in the bond was called ``CA``.

* ``bonds in X`` - a contractive search that finds all bonds that
  are contained wholly within the result of searching for ``X``,
  e.g. ``bonds in residx 0`` returns all bonds that are wholly
  within (both atoms within) the first residue.

* ``bonds with atoms in X`` - an expansive search that finds all bonds
   involve any atoms that is returned by the search ``X``. For example,
   ``bonds with atoms in resname ALA`` returns all bonds that have
   any atoms in residues called ``ALA``.

The above two searches find either all the bonds inside ``X``, or
all the bonds that involve atoms in ``X``. We need other keywords
to find the bonds that specifically bridge between two views.
These are ``to`` and ``from .. to``.

* ``bonds to X`` - a search that returns all bonds that have one
  atom contained in ``X`` and one atom that is not contained in ``X``
  (i.e. all bonds that connect ``X`` to another view). For example,
  ``bonds to resnum 1`` returns all bonds that connect residue number
  1 to any other residue (or view).

* ``bonds from X to Y`` - a search that returns all bonds that have
  one atom in ``X`` and one atom in ``Y``, i.e. that connect ``X``
  and ``Y``. For example, ``bonds from resnum 1 to resnum 2`` returns
  all the bonds between residue numbers 1 and 2.

Searching by Chemical Element
-----------------------------

You can search for views that match or contain atoms with specified
chemical elements. You do this using ``element``. This search
returns atom views. For example, ``element C`` would return all atom views
that have the chemical element carbon.

The chemical element can be specified in a number of different ways:

* ``element C`` - specify the chemical element by its IUPAC symbol, e.g.
  ``element H``, ``element Na`` etc.

* ``element c`` - specify the chemical element using its lowercase
  symbol, e.g. ``element h``, ``element na`` etc.

* ``element carbon`` - specify the chemical element using its lowercase
  name, e.g. ``element hydrogen``, ``element sodium`` etc.

* ``element C, H, Na`` - specify a list of elements. The atoms returned
   are those that match any in the list. The list can use any of the
   ways of specifying elements as above, e.g. ``element carbon, hydrogen``
   would be valid, as would ``element carbon, H, na``.

* ``element biological`` - specify the elements that are considered to
  be "biological", i.e. those in period 3 or less, which are not
  halogens or noble gases (the same definition used by the
  :func:`sire.mol.Element.biological` function).
