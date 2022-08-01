===============
Atom Properties
===============

Many of the properties you will use will be atomic properties. These
are properties that have one value per atom in the molecule.

For example, it is very common that molecules will have atomic coordinates.
By default, these are placed into a property called `coordinates`.

>>> print(mol.property("coordinates"))
AtomCoords( size=22
0: ( 18.4532, 3.49423, 12.4365 )
1: ( 18.9818, 3.44823, 13.3886 )
2: ( 20.0513, 3.63293, 13.2874 )
3: ( 18.798, 2.43076, 13.7337 )
4: ( 18.4805, 4.54971, 14.3514 )
...
17: ( 15.3407, 5.44815, 17.9626 )
18: ( 13.8341, 3.93668, 18.3509 )
19: ( 14.3525, 3.40994, 19.1521 )
20: ( 13.1933, 4.59022, 18.9428 )
21: ( 13.2149, 3.33301, 17.6874 )
)

The coordinates are held in a :class:`sire.mol.AtomCoords` object,
which is an example of a `AtomProperty` object. This behaves like
a python list, e.g.

>>> coords = mol.property("coordinates")
>>> print(coords[0])
( 18.4532, 3.49423, 12.4365 )
>>> print(coords[0:5])
[( 18.4532, 3.49423, 12.4365 ), ( 18.9818, 3.44823, 13.3886 ),
 ( 20.0513, 3.63293, 13.2874 ), ( 18.798, 2.43076, 13.7337 ),
 ( 18.4805, 4.54971, 14.3514 )]
>>> for coord in coords[0:3]:
...     print(coord)
( 18.4532, 3.49423, 12.4365 )
( 18.9818, 3.44823, 13.3886 )
( 20.0513, 3.63293, 13.2874 )

The coordinates themselves are :class:`sire.maths.Vector` objects. You
can get the x, y, and z components using the corresponding functions;

>>> coord = coords[0]
>>> print(coord.x())
18.4532476
>>> print(coord.y(), coord.z())
3.4942278 12.4364968

Accessing an atom property via the molecule will return the complete
:class:`~sire.mol.AtomCoords` object, containing the coordinates for
all of the atoms in the molecule.

You can get the `coordinates` properties for an individual atom using
the `property` function on the atom. For example, to get the coordinates
on the first atom in the molecule you could use;

>>> atom = mol[0]
>>> print(atom.property("coordinates"))
( 18.4532, 3.49423, 12.4365 )

To get the charge on the `CH3` atom in residue number `1` you could use

>>> print(mol["resnum 1"]["CH3"].property("charge"))
-0.3662 |e|

Convenience functions
---------------------

To reduce the amount of typing, there are shorthand, convenience functions
that can be used via the `Atom` view to access commonly-used properties.
The `coordinates` property can be accessed via the `coordinates` or
`coords` functions, e.g.

>>> print(atom.coordinates())
( 18.4532, 3.49423, 12.4365 )
>>> print(atom.coords())
( 18.4532, 3.49423, 12.4365 )

You can also get the `x`, `y`, and `z` components
of the coordinates directly, e.g.

>>> print(atom.x(), atom.y(), atom.z())
18.4532476 3.4942278 12.4364968

Properties that can be accessed this way are `charge`, `coordinates`,
`element`, `lj` (Lennard Jones parameters) and `mass`.

>>> print(atom.charge(), atom.element(), atom.lj(), atom.mass())
0.1123 |e| Hydrogen (H, 1)
LJ( sigma = 2.64953 A, epsilon = 0.0157 kcal mol-1 ) 1.008 g mol-1

These convenience functions can also be used for larger views. However,
in these cases they evaluate a single value that represents that
view from all of the values of the atoms in that view. Calling the `mass()`
on a molecule will return the total mass of all atoms in that molecule;

>>> print(mol.mass())
144.176 g mol-1

Similarly, calling `charge()` on a residue will return the total
charge on that residue;

>>> print(mol["resnum 1"].charge())
5.48778e-10 |e|

while calling `coordinates()` or `coords()` on a view will return
the center of mass of that view

>>> print(mol["resnum 1"].coords())
( 18.9264, 4.47803, 14.1498 )

This works for any view into a molecule, e.g. the total mass of the
first five atoms could be calculated via

>>> print(mol[0:5].mass())
27.044 g mol-1

Accessing atom properties via views
-----------------------------------

While the above convenience functions are useful, there are times when
you will want to get the individual atom properties for all atoms in
a view. You can do this by calling the `property()` function on that
view.

For example, to get the elements of all of the atoms in the first residue
you would use

>>> residue = mol["resnum 1"]
>>> print(residue.property("element"))
[Hydrogen (H, 1), Carbon (C, 6), Hydrogen (H, 1),
 Hydrogen (H, 1), Carbon (C, 6), Oxygen (O, 8)]

or, more directly

>>> print(mol["resnum 1"].property("element"))
[Hydrogen (H, 1), Carbon (C, 6), Hydrogen (H, 1),
 Hydrogen (H, 1), Carbon (C, 6), Oxygen (O, 8)]

This works for collections of views, e.g. to get all of the coordinates
on the first five atoms of the molecule, you would use

>>> print(mol[0:5].property("coordinates"))
[( 18.4532, 3.49423, 12.4365 ), ( 18.9818, 3.44823, 13.3886 ),
 ( 20.0513, 3.63293, 13.2874 ), ( 18.798, 2.43076, 13.7337 ),
 ( 18.4805, 4.54971, 14.3514 )]

or you could get the Lennard Jones parameters of all of the carbon
atoms using

>>> print(mol["element C"].property("LJ"))
[LJ( sigma = 3.39967 A, epsilon = 0.1094 kcal mol-1 ),
 LJ( sigma = 3.39967 A, epsilon = 0.086 kcal mol-1 ),
 LJ( sigma = 3.39967 A, epsilon = 0.1094 kcal mol-1 ),
 LJ( sigma = 3.39967 A, epsilon = 0.1094 kcal mol-1 ),
 LJ( sigma = 3.39967 A, epsilon = 0.086 kcal mol-1 ),
 LJ( sigma = 3.39967 A, epsilon = 0.1094 kcal mol-1 )]

Using apply to get the properties of views in a container
---------------------------------------------------------

Another route to getting the properties is to use the `apply` function.
The `apply` function will call the passed function on all views
within a molecular container. For example, calling the `charge` function
on `mol.atoms()` will return the total charge on the molecule,

>>> print(mol.atoms().charge())
-5.48778e-10 |e|

To call the `charge` function on each atom in the `mol.atoms()` container,
we would use `apply`, e.g.

>>> print(mol.atoms().apply("charge"))
[ 0.1123 |e|, -0.3662 |e|, 0.1123 |e|, 0.1123 |e|, 0.5972 |e|, -0.5679 |e|,
 -0.4157 |e|, 0.2719 |e|, 0.0337 |e|, 0.0823 |e|, -0.1825 |e|, 0.0603 |e|,
  0.0603 |e|, 0.0603 |e|, 0.5973 |e|, -0.5679 |e|, -0.4157 |e|, 0.2719 |e|,
 -0.149 |e|, 0.0976 |e|, 0.0976 |e|, 0.0976 |e|]

Apply calls the specified function on each view in a container, returning
the result as a list. You can either pass in the name of the function
you want to apply, or you can pass in a function yourself. In this case,
we will use `apply` with a lambda expression to get the x coordinates
of all of the atoms in the first residue;

>>> print(mol["resnum 1"].apply(lambda atom: atom.x()))
[18.4532476, 18.9818473, 20.0513009, 18.7980053, 18.4805331, 19.1865897]

You can pass in positional and named arguments to the applied function
as arguments to `apply`. For example, here we will ask for the `mass`
property on each atom by calling the `property` function via an `apply`;

>>> print(mol.atoms().apply("property", "mass"))
[1.008 g mol-1, 12.01 g mol-1, 1.008 g mol-1, 1.008 g mol-1, 12.01 g mol-1,
 16 g mol-1, 14.01 g mol-1, 1.008 g mol-1, 12.01 g mol-1, 1.008 g mol-1,
 12.01 g mol-1, 1.008 g mol-1, 1.008 g mol-1, 1.008 g mol-1, 12.01 g mol-1,
 16 g mol-1, 14.01 g mol-1, 1.008 g mol-1, 12.01 g mol-1, 1.008 g mol-1,
 1.008 g mol-1, 1.008 g mol-1]
