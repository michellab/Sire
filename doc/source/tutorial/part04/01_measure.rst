==============================
Measuring Distances and Angles
==============================

We have already seen that we can measure the lengths of bonds or
sizes of angles using the functions of the
:class:`~sire.mm.Bond`, :code:`~sire.mm.Angle`, :code:`~sire.mm.Dihedral`
and :class:`~sire.mm.Improper` classes.

For example, load up the ``aladip`` system...

>>> import sire as sr
>>> mols = sr.load(sr.expand(sr.tutorial_url, ["ala.top", "ala.crd"]))
>>> mol = mols[0]

and we could then find the lengths of all of the carbon-oxygen bonds using
the :func:`~sire.mm.SelectorBond.lengths` function;

>>> print(mols.bonds("element carbon", "element oxygen").lengths())
[1.20803 angstrom, 1.24385 angstrom]

Similarly, you could get the size of the first five hydrogen-oxygen-hydrogen
angles using the :func:`~sire.mm.SelectorAngle.sizes` function;

>>> print(mols.angles("element H", "element O", "element H")[0:5].sizes())
[104.491 degree, 104.491 degree, 104.491 degree, 104.491 degree, 104.491 degree]

Sire uses the synonym ``measure`` for both lengths and sizes. This lets you
use the same function name for measuring bond lengths, angle sizes or
dihedral torsion sizes. For example, you could use
:func:`~sire.mm.SelectorBond.measures` in place of
:func:`~sire.mm.SelectorBond.lengths` above, e.g.

>>> print(mols.bonds("element carbon", "element oxygen").measures())
[1.20803 angstrom, 1.24385 angstrom]

or :func:`~sire.mm.SelectorAngle.measures` in place of
:func:`~sire.mm.SelectorAngle.sizes`,

>>> print(mols.angles("element H", "element O", "element H")[0:5].measures())
[104.491 degree, 104.491 degree, 104.491 degree, 104.491 degree, 104.491 degree]

.. note::

    You can use ``measures`` instead of ``sizes`` to measure dihedral
    angles and improper angles too.

    For individual bond, angle, dihedral, or improper objects, you could
    use ``measure`` instead of ``length`` or ``size``.

Making measurements between atoms
---------------------------------

So far, you have used ``length``, ``size`` or ``measure`` for measuring
actual class:`~sire.mm.Bond`, :code:`~sire.mm.Angle`, :code:`~sire.mm.Dihedral`
or :class:`~sire.mm.Improper` objects.

There are many cases where you want to measure the distance or angles
between arbitrary atoms.

To do this, you use the :func:`sire.measure` function. For example,
we could measure the distance between the oxygen atoms of the first
two water molecules using;

>>> oxygens = mols["water and element O"]
>>> print(sr.measure(oxygens[0], oxygens[1]))
18.5067 angstrom

The measurement returned depends on the number of items passed to the
:func:`~sire.measure` function. Passing two items, as above, will measure
and return the distance. Passing three items will measure and return
the angle, so

>>> print(sr.measure(oxygens[0], oxygens[1], oxygens[2]))
53.3414 degree

has returned the angle between the oxygens of the first three water
molecules.

Passing in four items will measure the dihedral (torsion) angle, i.e.

>>> print(sr.measure(oxygens[0], oxygens[1], oxygens[2], oxygens[3]))
60.0107 degree

measures the torsion angle between the oxygens of the first four
water molecules.

Improper angles are also measured between four items. Set
``improper_angle`` to ``True`` to get the improper angle instead;

>>> print(sr.measure(oxygens[0], oxygens[1], oxygens[2], oxygens[3],
...                  improper_angle=True))
-44.0118 degree

Passing in only a single item will call the ``.measure()`` function
on that item. This means that this will only work for individual
:class:`~sire.mm.Bond`, :code:`~sire.mm.Angle`, :code:`~sire.mm.Dihedral`
or :class:`~sire.mm.Improper` objects;

>>> bond = mols[0].bonds()[0]
>>> print(bond, bond.measure())
Bond( HH31:1 => CH3:2 ) 1.09 angstrom
>>> print(sr.measure(bond))
1.09 angstrom

Making measurements between arbitray views
------------------------------------------

The :func:`~sire.measure` function calls the ``.coordinates()`` function
on the items that are passed. This means that you can pass in any
object that has a ``.coordinates()`` function. For example, you can
calculate the distance between the first two water molecules using

>>> waters = mols["water"]
>>> print(sr.measure(waters[0], waters[1]))
18.4583 angstrom

This is not the same as the distance between the oxygens of
these water molecules. This is because the ``.coordinates()``
function on a molecule returns the molecule's center of mass.

If you wanted to return the distance between the molecules' centers
of geometry you would use

>>> print(sr.measure(waters[0].evaluate().center_of_geometry(),
...                  waters[1].evaluate().center_of_geometry()))
18.0674 angstrom

... sub views, introduce sire.maths.Vector object.


