=========
Searching
=========

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
* A stepped range, e.g. ``molidx 0:10:2`` would search for molecules with
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
* A case-insensitive regular expression glob pattern, e.g. ``resname /ala/i`` matches any
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

Whether a search is expansive  (``with``) or contractive (``in``) depends on the relative
size of the view being searched, and the view that is being returned
as the result. This can be summarised as a table.

+----------+---------------------+---------------------+------------------------+----------------------+------------------------+----------------------+
|          | Atom                | Bond                | Residue                | Chain                | Segment                | Molecule             |
+----------+---------------------+---------------------+------------------------+----------------------+------------------------+----------------------+
| Atom     | N/A                 | atoms in bond       | atoms in residue       | atoms in chain       | atoms in segment       | atoms in molecule    |
+----------+---------------------+---------------------+------------------------+----------------------+------------------------+----------------------+
| Bond     | bonds with atom     | N/A                 | bonds in residue       | bonds in chain       | bonds in segment       | bonds in molecule    |
+----------+---------------------+---------------------+------------------------+----------------------+------------------------+----------------------+
| Residue  | residues with atom  | residues with bond  | N/A                    | residues in chain    | residues in segment    | residues in molecule |
+----------+---------------------+---------------------+------------------------+----------------------+------------------------+----------------------+
| Chain    | chains with atom    | chains with bond    | chains with residue    | N/A                  | chains in segment      | chains in molecule   |
+----------+---------------------+---------------------+------------------------+----------------------+------------------------+----------------------+
| Segment  | segments with atomv | segments with bond  | segments with residue  | segments with chain  | N/A                    | segments in molecule |
+----------+---------------------+---------------------+------------------------+----------------------+------------------------+----------------------+
| Molecule | molecules with atom | molecules with bond | molecules with residue | molecules with chain | molecules with segment | N/A                  |
+----------+---------------------+---------------------+------------------------+----------------------+------------------------+----------------------+

.. note::

   ``X in Y`` is a contractive search that looks for the smaller views of
   type ``X`` within the larger view ``Y``.

   ``Y with X`` is an expansive search that returns larger views of type
   ``X`` that contain the smaller views of type ``Y``.

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

* ``atoms`` - return all of the atoms in the searched object. Also abbreviated
  to ``atom``
* ``bonds`` - return all of the bonds in the searched object. Also abbreviated
  to ``bond``.
* ``residues`` - return all of the residues in the searched object. Also
  abbreviated to ``residue`` or ``res``.
* ``chains`` - return all of the chains in the searched object. Also
  abbreviated to ``chain``.
* ``segments`` - return all of the segments in the searched object.
  Also abbreviated to ``segment``, ``segs`` and ``seg``.
* ``molecules`` - return all of the molecules in the searched object.
  Also abbreviated to ``molecule``, ``mols`` and ``mol``.
* ``all`` or ``*`` - everything in the current view being searched.

Whether these are expansive or contractive depends on the view that
is searched, based on the same rules as for normal searches
(e.g. searching for ``residues`` on an atoms view will be expansive,
as it returns all residues that contain those atoms, while searching
for ``residues`` on a molecule view will be contractive, as it searches
for all residues in that molecule)

Advanced ``with`` and ``in`` Searches
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
   involving any atoms that is returned by the search ``X``. For example,
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
  :func:`sire.mol.Element.biological` function). Note you can use
  the shorthand ``element bio`` to also match biological atoms.

Searching by Count (i.e. Number of Atoms)
-----------------------------------------

You can search by counts, e.g. finding all molecules with more than
one residue, using ``count``.

You can do simple searches of the form ``view_type with/in count(view_type) compare number``, e.g.

* ``molecules with count(residues) > 1`` - match all molecules that
  contain more than on residue
* ``atoms in count(residues) == 1`` - match all atoms that are in molecules
  that contain just one residue
* ``residues in molecules with count(atoms) < 20`` - match all residues that
  are in molecules with less than twenty atoms

You can also use more advanced ``with/in`` searches, e.g.

* ``element C in residues with count(atoms) > 5`` - match all carbon atoms
  in residues that have more than five atoms
* ``element O in molecules with count(atoms) == 3`` - match all oxygen atoms
  in molecules that have three atoms
* ``resname /ALA/i in molecules with count(residues) > 20`` - match all
  residues call ``ALA`` (in any case) in molecules that have more
  than 20 residues.

In the general case, the argument to ``count`` is actually a search too!
So the full syntax is ``X with/in count(Y) compare number``, which
returns views that match ``X`` where the count of views that match ``Y``
in the view being searched compares true with the specified number.
For example;

* ``residues with count(element C) > 2`` - match all residues that contain
  more than two carbon atoms
* ``atoms in residues with count(atomname CA) == 1`` - match all atoms
  in residues that contain a single atom called ``CA``
* ``(molecules with count(element O) == 1) and (molecules with count(element H) == 2) and (molecules with count(atoms) == 3)`` -
  match all molecules that contain there atoms, and that have one oxygen atom
  and have two hydrogen atoms (i.e. are water molecules).

Searching by Charge or Mass
---------------------------

You can search for views by their charge or mass using the ``charge`` or
``mass`` keywords. The grammar is ``X with charge compare number``
or ``X with mass Y compare number``, where ``X`` is the view (or views)
you want to calculate the mass or charge for, ``compare`` is the
comparison operator (e.g. ``==``, ``<=`` etc.), and ``number`` is the
value of the charge or mass you want to compare to. For example,

* ``atoms with charge < 0`` - return all atoms that have a change that
  is less than zero.
* ``atoms with mass > 4`` - return all atoms that have a mass that is
  greater than ``4 g mol-1``.
* ``residues with charge >= 0.8`` - return all residues that have a total
  charge of greater than or equal to ``0.8 |e|`` (modulo electron charge units).
* ``molecules with mass < 50`` - return all molecules whose total mass is
  less than ``50 g mol-1``.

In the general case, ``X`` can be any search. So you could use,

* ``element C with charge > 0.1`` - return all carbon atoms with a charge
  greater than ``0.1 |e|``.
* ``resname LIG with mass > 100`` - return all residues called ``LIG``
  that have a mass greater than ``100 g mol-1``.

You could also combine this with other searches, e.g.

* ``mols with count(residues with charge > 0.5) > 5`` - return all molecules
  with more than five residues that have a charge of more than ``0.5 |e|``.

Numerical imprecision can cause issues with ``charge`` searches.
This is because rounding errors can mean that the sum of charges across,
e.g. a residue or molecule, are non-integer. For example, for
the ``aladip`` example used in the tutorials, the charges of the first
two residues are ``5.48778e-10`` and ``-1.09756e-09``. Searches for
residues with zero charge, or with negative charge would fail for
these two residues.

Searching by Charge using Approximate Comparisons
-------------------------------------------------

The approximate comparison operators can help solve the numerical
imprecision issues experienced when searching for views by charge.

The approximate comparison operators are;

* ``=~`` - approximate equal to. :func:`sire.search.approx_equal`
* ``!~`` - not approximates equal to. :func:`sire.search_approx_not_equal`
* ``>~`` - greater than (but not approximately equal to). :func:`sire.search.approx_greater`
* ``<~`` - less than (but not approximately equal to). :func:`sire.search.approx_less`
* ``>=~`` - greater than or approximately equal to. :func:`sire.search.approx_greater_equal`
* ``<=~`` - less than or approximatley equal to. :func:`sire.search.approx_less_equal`

These operators are logically consistent. They are built from
the approximately equal to operator (``=~``). This returns ``True``
if the two values compared are equal, or are equal to within an
epsilon value. This is set via :func:`sire.search.set_approx_epsilon` and
can be retrieved using :func:`sire.search.get_approx_epsilon`.

These operators can be used to, e.g.

* ``mols with charge =~ 0`` - find all neutral molecules
* ``residues with charge >~ 0`` - find all positively charged residues
* ``(molecules with count(atoms) == 1) with charge =~ 1`` - all ``+1`` charge ions

.. note::

   Brackets were used for the last search to ensure that we look for
   single-atom molecules first, and then, from these, find those that
   have a charge of +1.

Searching by Property
---------------------

You can search for views based on any of the properties of those views.
There are several routes to do this;

* ``X with property name`` - return all views that have a property
  called ``name``, and if so, the value of this property is not
  ``False`` or ``0``, ``0.0``, or the string ``false`` in any case.
  For example, ``molecules with property is_perturbable`` would
  return all molecules that have the ``is_perturbable`` property
  and that it wasn't false. ``atoms with property atomtype``
  would return all atoms that have an ``atomtype`` property, t
  and if so, that it wasn't false.

* ``X with property name compare value`` - return all views that have
  a property called ``name``, and where the value of the property
  compares truthfully using the comparison operator ``compare``
  against the value ``value``. For example,
  ``residues with property is_default == False`` would return
  all residues that have an ``is_default`` property, and this
  property is equal to ``False``. Searching for
  ``atoms with property radius > 0.5`` would return all atoms
  that have a ``residue`` property, and the value of this is
  greater than ``0.5``, while ``atoms with property atomtype == HA``
  would return all atoms with an ``atomtype`` property which
  is equal to ``HA``.

.. note::

   The values in property searches should be in the default units
   for the property being searched (e.g. ``radius > 0.5`` is in
   units of Å as this is the default length unit). Remember this
   if you change the default length unit, i.e. if the default
   length unit is picometers, then the above search would be
   ``radius > 50``.

Finding the Nth View that Matches
---------------------------------

You can use subscripting to pick out the nth view that matches a particular
search. The grammar is ``X[i]`` where ``X`` is the search, and ``i``
is the index of the result you want to match.

* ``element C[0]`` - return the first carbon atom
* ``(resname ALA)[-1]`` - return the last residue called ``ALA``
* ``(bonds with element H)[0:5]`` - return the first five bonds that
  contain hydrogen.

Searching by Distance
---------------------

You can search for views based on their distances to either each other,
or to fixed points in space. There are a few variants of this search.

The first is ``X within D of Y``. This returns views that match ``X`` that are
within the specified distance ``D`` (in default length units) of views that match
``Y``. Examples include ``waters within 5 of resname ALA``, which finds
all water molecules within 5 Å of residues called ``ALA``. Or
``(element H in water) within 3 of (element O, Cl, F in protein)``
which finds all hydrogen atoms in water molecules that are within
3 Å of all oxygen, chlorine or fluorine atoms in protein molecules.

The distance criterion used is whether any atoms in the view ``X`` is
within the specified distance of any atom in the view ``Y``.

You can specify the distance criterion using searches of form
``X where CRITERION within D of Y``. The ``CRITERION`` has many possible
values;

* ``center`` - the geometric center of the atoms in the two views will
  be compared, e.g. ``waters where center within 5 of resname ALA`` finds
  all waters whose geometric centers are within 5 Å of the geometric
  centers of all ``ALA`` residues combined.

* ``max`` - the maximum coordinates of the two views will be compared,
  e.g. ``atoms where max within 5 of protein`` finds all atoms that
  are within 5 Å of the maximum coordinates of all proteins combined.

* ``min`` - the minimum coordinates.

* ``A.x``, ``A.y``, ``A.z`` - compare specifically the x, y or z
  dimensions only, e.g. ``center.x`` compares only the x dimension
  of the centers of the views, while ``max.y`` would compare the
  y dimension of the maximum, and ``min.z`` would be the z dimension
  of the minimum.

You can compare against fixed points in space by passing in a vector
in place of ``Y``. For example, ``atoms within 10 of (0,0,0)``
would find all atoms within 10 Å of the origin, while
``waters within 3 of (5.3, 9.8, -2.1)`` would find all water molecules
within 3 Å of the point ``(5.3, 9.8, -2.1)``.

.. note::

   All distance searches use the first ``space`` property that can be
   found in the views to calculate distances. This should ensure that
   periodic boundaries are accounted for in the search. You can specify
   the space to use by passing this in via a property map.

Searching for Water or Protein Molecules
----------------------------------------

You can search for protein or water molecules using the ``protein``
or ``water`` keywords.

* ``water`` - search for all water molecules. These are defined as
  molecules that contain one oxygen, two hydrogen and any number
  of null-element (dummy) atoms.

* ``protein`` - search for all protein molecules. These are defined
  as molecules that contain a minimum number of residues whose
  names are found in the list of possible protein residue names.

You can get and set the minimum number of protein residues to match
using the :func:`sire.search.get_min_protein_residues` and
:func:`sire.search.set_min_protein_residues` functions. The default
minimum is 5.

You can get and set the list of protein residue names using the
:func:`sire.search.get_protein_residue_names` and
:func:`sire.search.set_protein_residue_names` functions. The names
are searched via a case-insensitive search. The default list
are the standard names of the biological amino acids, including
names commonly used for protonated or deprotonated residues.

Creating Custom Search Tokens
-----------------------------

You can create and use your own search tokens! You set new search tokens
using :func:`sire.search.set_token`, e.g.

>>> sire.search.set_token("CH-bonds", "bonds from element C to element H")

would create a new search token called ``CH-bonds`` that finds all bonds
between carbon and hydrogen atoms. You could then use this token just
like any other token, e.g.

>>> mols["CH-bonds in protein"]

would return all the carbon-hydrogen bonds in protein molecules.

Token names cannot contain spaces, and must start with a letter (i.e.
that can't start with numbers or symbols). Token names cannot be valid
search terms themselves - meaning you can't use custom tokens to redefine
the above tokens or grammar. They also can't start with any existing search
term or token. This means that, once set, you cannot
set the token again, or, indeed, set any token that starts with this token.
To change a token you first have to delete it
using the :func:`sire.search.delete_token` function, e.g.

>>> sire.search.delete_token("CH-bonds")

.. note::

   This function will silently do nothing if the token doesn't exist,
   or if you attempt to delete an built-in token. You cannot delete
   built-in tokens, e.g. like ``protein`` or ``water``.

You can use :func:`sire.search.delete_all_tokens` to delete all
custom search tokens.

The function :func:`sire.search.has_token` will query whether or not
a particular token exists.

The function :func:`sire.search.get_token` will return the search term
that corresponds to the passed token, e.g.

>>> sire.search.set_token("CH-bonds", "bonds from element C to element H")
>>> print(sire.search.get_token("CH-bonds"))
bonds from element C to element H

Tokens are expanded when they are used as part of other tokens, e.g.

>>> sire.search.set_token("p-CH-bonds", "CH-bonds in protein")
>>> print(sire.search.get_token("p-CH-bonds"))
({ CH-bonds => bonds from element C to element H }) in (protein)

This shows that ``CH-bonds`` is ``bonds from element C to element H``.

Deleting or changing the token later does not affect any pre-existing
tokens that used it, e.g.

>>> sire.search.delete_token("CH-bonds")
>>> print(sire.search.get_token("p-CH-bonds"))
({ CH-bonds => bonds from element C to element H }) in (protein)
>>> sire.search.set_token("CH-bonds", "bonds to element C")
>>> print(sire.search.get_token("p-CH-bonds"))
({ CH-bonds => bonds from element C to element H }) in (protein)
