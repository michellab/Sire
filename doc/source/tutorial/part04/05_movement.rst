==========================
Moving Atoms and Molecules
==========================

Moving atoms and molecules is very simple. At the most basic, you
can move atoms by simply updating their coordinates property
(just as you would update any other property). For example,
let's load up the ``aladip`` system and translate the first
molecule by 1 Å along the x axis.

First, lets get the molecule and print the coordinates as they
are from the input file...

>>> mols = sr.load(sr.expand(sr.tutorial_url, ["ala.top", "ala.traj"]))
>>> mol = mols[0]
>>> print(mol.property("coordinates"))
AtomCoords( size=22
0: ( 18.9652 Å, 3.62266 Å, 12.5364 Å )
1: ( 19.4344 Å, 3.49272 Å, 13.5116 Å )
2: ( 20.4907 Å, 3.75359 Å, 13.4455 Å )
3: ( 19.3447 Å, 2.46862 Å, 13.8739 Å )
4: ( 18.6496 Å, 4.35419 Å, 14.4994 Å )
...
17: ( 14.8153 Å, 5.52081 Å, 17.2499 Å )
18: ( 13.948 Å, 3.8589 Å, 18.2705 Å )
19: ( 14.5959 Å, 3.37224 Å, 18.9996 Å )
20: ( 13.2209 Å, 4.57526 Å, 18.6529 Å )
21: ( 13.4889 Å, 2.96239 Å, 17.8539 Å )
)

Now let's edit the coordinates. The coordinates are just another property
of the molecule, so can be edited using a cursor.

>>> c = mol.cursor()
>>> for atom in c.atoms():
...     atom["coordinates"] = atom["coordinates"] + (1,0,0)
>>> mol = c.commit()
>>> print(mol.property("coordinates"))
AtomCoords( size=22
0: ( 19.9652 Å, 3.62266 Å, 12.5364 Å )
1: ( 20.4344 Å, 3.49272 Å, 13.5116 Å )
2: ( 21.4907 Å, 3.75359 Å, 13.4455 Å )
3: ( 20.3447 Å, 2.46862 Å, 13.8739 Å )
4: ( 19.6496 Å, 4.35419 Å, 14.4994 Å )
...
17: ( 15.8153 Å, 5.52081 Å, 17.2499 Å )
18: ( 14.948 Å, 3.8589 Å, 18.2705 Å )
19: ( 15.5959 Å, 3.37224 Å, 18.9996 Å )
20: ( 14.2209 Å, 4.57526 Å, 18.6529 Å )
21: ( 14.4889 Å, 2.96239 Å, 17.8539 Å )
)

.. note::

   Note that the tuple, ``(1,0,0)`` is automatically converted to a
   ``Vector`` when added to the coordinates. The units will be the current
   default length units.

In this case, we have translated every atom by 1 Å along the x-axis.
You could set the coordinates directly. For example, here we set the
coordinates of the hydrogen atoms to ``(0,0,0)``.

>>> for atom in c.atoms("element H"):
...     atom["coordinates"] = (0, 0, 0)
>>> mol = c.commit()
>>> print(mol.property("coordinates"))
AtomCoords( size=22
0: ( 0 , 0 , 0  )
1: ( 20.4344 Å, 3.49272 Å, 13.5116 Å )
2: ( 0 , 0 , 0  )
3: ( 0 , 0 , 0  )
4: ( 19.6496 Å, 4.35419 Å, 14.4994 Å )
...
17: ( 0 , 0 , 0  )
18: ( 14.948 Å, 3.8589 Å, 18.2705 Å )
19: ( 0 , 0 , 0  )
20: ( 0 , 0 , 0  )
21: ( 0 , 0 , 0  )
)

Moving molecules using move functions
-------------------------------------

The :class:`~sire.mol.Cursor` has "move" functions that simplify
the process of translating, rotating and moving atoms and molecules.

For example, this is how you can use a cursor more easily translate
the first molecule by 1 Å along the x axis.

>>> mol = mol.cursor().translate( (1,0,0) ).commit()
>>> print(mol.property("coordinates"))

The :func:`~sire.mol.Cursor.translate` function translates all of the
atoms selected by the cursor by the passed vector (or passed x, y and
z components).

For example, you could translate all of the hydrogen atoms by
the vector ``(1, 2, 3)`` using

>>> cursor = mol.cursor()
>>> cursor["element H"].translate(1, 2, 3)
>>> mol = cursor.commit()
>>> print(mol.property("coordinates"))

