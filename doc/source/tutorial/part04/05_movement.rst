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

>>> mol = mols[0]
>>> mol = mol.cursor().translate( (1,0,0) ).commit()
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

The :func:`~sire.mol.Cursor.translate` function translates all of the
atoms selected by the cursor by the passed vector (or passed x, y and
z components).

For example, you could translate all of the hydrogen atoms by
the vector ``(1, 2, 3)`` using

>>> cursor = mol.cursor()
>>> cursor["element H"].translate(1, 2, 3)
>>> mol = cursor.commit()
>>> print(mol.property("coordinates"))
AtomCoords( size=22
0: ( 20.9652 Å, 5.62266 Å, 15.5364 Å )
1: ( 20.4344 Å, 3.49272 Å, 13.5116 Å )
2: ( 22.4907 Å, 5.75359 Å, 16.4455 Å )
3: ( 21.3447 Å, 4.46862 Å, 16.8739 Å )
4: ( 19.6496 Å, 4.35419 Å, 14.4994 Å )
...
17: ( 16.8153 Å, 7.52081 Å, 20.2499 Å )
18: ( 14.948 Å, 3.8589 Å, 18.2705 Å )
19: ( 16.5959 Å, 5.37224 Å, 21.9996 Å )
20: ( 15.2209 Å, 6.57526 Å, 21.6529 Å )
21: ( 15.4889 Å, 4.96239 Å, 20.8539 Å )
)

.. note::

   You can pass in the vector to translate either as arguments, e.g.
   ``translate(1, 2, 3)``, or as a ``Vector``, e.g.
   ``translate(sr.maths.Vector(1,2,3))`` or ``translate((1,2,3))``.
   As for the rest of Sire, the default units are Å, which can be
   changed using, e.g. ``sr.units.set_length_unit``. You can also
   specify the units yourself, e.g. ``translate(1*sr.units.angstrom, 0, 0)``
   or ``translate(sr.maths.Vector(1*sr.units.angstrom,0, 0))``.

You can even translate all of the molecules that have been loaded,
using the cursor for the whole system.

>>> cursor = mols.cursor()
>>> cursor.translate(3,4,5)
>>> mols = cursor.commit()
>>> print(mols[1].property("coordinates"))
AtomCoords( size=3
0: ( 27.9816 Å, 12.1269 Å, 27.3254 Å )
1: ( 28.4376 Å, 12.3231 Å, 26.5069 Å )
2: ( 28.6385 Å, 11.689 Å, 27.8666 Å )
)

You can rotate molecules using a cursor's :func:`~sire.mol.Cursor.rotate`
function.

>>> mol = mols[0]
>>> cursor = mol.cursor()
>>> cursor.rotate(5)
>>> mol = cursor.commit()
>>> print(mol.property("coordinates"))
AtomCoords( size=22
0: ( 22.0096 Å, 7.82943 Å, 17.5364 Å )
1: ( 22.4883 Å, 7.74088 Å, 18.5116 Å )
2: ( 23.5178 Å, 8.09281 Å, 18.4455 Å )
3: ( 22.4883 Å, 6.71286 Å, 18.8739 Å )
4: ( 21.6315 Å, 8.53067 Å, 19.4994 Å )
...
17: ( 17.71 Å, 9.35867 Å, 22.2499 Å )
18: ( 16.9909 Å, 7.62749 Å, 23.2705 Å )
19: ( 17.6787 Å, 7.19915 Å, 23.9996 Å )
20: ( 16.2041 Å, 8.27775 Å, 23.6529 Å )
21: ( 16.6117 Å, 6.69438 Å, 22.8539 Å )
)

In this case, we rotated the molecule by 5° about the z-axis of the molecule,
around its center of mass.

You can specify the units yourself, e.g. ``5 * sr.units.degrees``, and can
also specify the axis and centers of rotation as additional arguments, e.g.

>>> cursor.rotate(0.1*sr.units.radians, (1,0,0))
>>> print(cursor["coordinates"])
AtomCoords( size=22
0: ( 22.0096 Å, 8.13227 Å, 17.511 Å )
1: ( 22.4883 Å, 7.94679 Å, 18.4725 Å )
2: ( 23.5178 Å, 8.30357 Å, 18.4418 Å )
3: ( 22.4883 Å, 6.88774 Å, 18.7304 Å )
4: ( 21.6315 Å, 8.63402 Å, 19.5342 Å )
...
17: ( 17.71 Å, 9.18329 Å, 22.3536 Å )
18: ( 16.9909 Å, 7.35887 Å, 23.1964 Å )
19: ( 17.6787 Å, 6.85988 Å, 23.879 Å )
20: ( 16.2041 Å, 7.9677 Å, 23.6418 Å )
21: ( 16.6117 Å, 6.47201 Å, 22.6886 Å )
)

rotates by 0.1 radians about the x-axis (``(1,0,0)``) around the
molecule's center of mass, while

>>> cursor.rotate(10*sr.units.degrees, (0,1,0), (0,0,0))
>>> print(cursor["coordinates"])
AtomCoords( size=22
0: ( 24.716 Å, 8.13227 Å, 13.423 Å )
1: ( 25.3544 Å, 7.94679 Å, 14.2868 Å )
2: ( 26.3629 Å, 8.30357 Å, 14.0778 Å )
3: ( 25.3991 Å, 6.88774 Å, 14.5408 Å )
4: ( 24.6949 Å, 8.63402 Å, 15.4812 Å )
...
17: ( 21.3227 Å, 9.18329 Å, 18.9387 Å )
18: ( 20.7608 Å, 7.35887 Å, 19.8935 Å )
19: ( 21.5567 Å, 6.85988 Å, 20.4463 Å )
20: ( 20.0633 Å, 7.9677 Å, 20.4688 Å )
21: ( 20.2992 Å, 6.47201 Å, 19.4593 Å )
)

rotates by 10° about the y-axis with the rotation centered on the origin
(``(0,0,0)``).

You can also specify the rotations directly via rotation matrices
(:class:`sire.maths.Matrix`) or quaternions (:class:`sire.maths.Quaternion`).

To do this, pass in the matrix or quaternion that represents the rotation, e.g.

>>> cursor.rotate(sr.maths.Quaternion(5*sr.units.degrees,
...                                   sr.maths.Vector(1,0,0)))
>>> print(cursor["coordinates"])
AtomCoords( size=22
0: ( 24.716 Å, 8.42963 Å, 13.4271 Å )
1: ( 25.3544 Å, 8.16958 Å, 14.2714 Å )
2: ( 26.3629 Å, 8.54321 Å, 14.0943 Å )
3: ( 25.3991 Å, 7.09242 Å, 14.4321 Å )
4: ( 24.6949 Å, 8.75009 Å, 15.5212 Å )
...
17: ( 21.3227 Å, 8.99593 Å, 19.0134 Å )
18: ( 20.7608 Å, 7.09524 Å, 19.8055 Å )
19: ( 21.5567 Å, 6.54996 Å, 20.3128 Å )
20: ( 20.0633 Å, 7.65161 Å, 20.4317 Å )
21: ( 20.2992 Å, 6.24959 Å, 19.2957 Å )
)

or

>>> rotmat = sr.maths.Matrix(1,0,0,
...                          0,0.984808,-0.173648,
...                          0,0.173648,0.984808)
>>> cursor.rotate(rotmat)
>>> print(cursor["coordinates"])

.. note::

   The above rotation matrix rotates by 10° about the x-axis.
   If was generated using the ``to_matrix()`` function of the
   :class:`~sire.maths.Quaternion` that represented this
   rotation.

As before, the center of rotation defaults to the center of mass
of the molecule. You can specify the center of rotation via the
``center`` keyword argument. For example,

>>> cursor.rotate(rotmat, center=(0,0,0))
>>> print(cursor["coordinates"])
AtomCoords( size=22
0: ( 24.716 Å, 6.5342 Å, 14.8733 Å )
1: ( 25.3544 Å, 6.00105 Å, 15.5778 Å )
2: ( 26.3629 Å, 6.41271 Å, 15.5391 Å )
3: ( 25.3991 Å, 4.93388 Å, 15.3604 Å )
4: ( 24.6949 Å, 6.11911 Å, 16.9507 Å )
...
17: ( 21.3227 Å, 5.15571 Å, 20.3164 Å )
18: ( 20.7608 Å, 3.09872 Å, 20.4107 Å )
19: ( 21.5567 Å, 2.41285 Å, 20.7008 Å )
20: ( 20.0633 Å, 3.40739 Å, 21.1894 Å )
21: ( 20.2992 Å, 2.47844 Å, 19.6424 Å )
)

rotates using the passed rotation matrix, centered on the origin.
