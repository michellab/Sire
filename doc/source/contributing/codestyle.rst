============
Coding Style
============

We really appreciate your help developing :mod:`sire` and welcome
pull requests to the ``devel`` branch. To help us more quickly
review your pull request, and to keep a consistent coding style
throughout, we ask that you please follow the below coding styles
(and that you aren't offended if we modify your submission so
that it meets these styles).

Python
======

Python code should be written to be `PEP8-compliant <https://pep8.org>`__,
ideally autoformatted using a tool such as
`black <https://black.readthedocs.io/en/stable/>`__.

C++
===

To keep the code consistent, there a strict coding style is used in
:mod:`sire`.

Please follow these rules to maintain this style. Also, if you see code
that doesn't follow these rules (we are all not perfect) then please
feel free to either correct that code, or to notify one of the lead
developers. This coding style has evolved over many years and many
developers, and is needed to ensure that the C++ source code looks
like a single homogenous set, which anyone can edit at any point.

All indentation should be multiples of 4 spaces. All tabs should be replaced by 4 spaces
----------------------------------------------------------------------------------------

For example

.. code-block:: c++

   void foo()
   {
       if (true)
       {
           for (int i=0; i<10; ++i)
           {
               //do something
           }
       }
   }

not

.. code-block:: c++

   void foo()
   {
    if (true)
    {
	   for (int i=0; i<10; ++i){ /* do something */ }
    }
   }

Curly brackets should be used for all blocks, with ``{`` on a new line
----------------------------------------------------------------------

For example

.. code-block:: c++

   for (int i=0; i<10; ++i)
   {
       if (condition)
       {
           //do something
       }
   }

not

.. code-block:: c++

   for (int i=0; i<10; ++i){
       if (condition) /* do something */;
   }

Classes should be named using capital letters, using only the letters A-Za-z and numbers 0-9
--------------------------------------------------------------------------------------------

Please do not use underscores.
For example ``BigMolecule`` is acceptable, but ``bigMolecule``, ``Bigmolecule``
``Big_Molecule`` or ``bigmolecule`` are not.

Functions (methods) should be named in the same way as classes, except that the first letter should not be capitalised
----------------------------------------------------------------------------------------------------------------------

For example
``getRadius()`` is acceptable, but ``GetRadius()``, ``getradius()`` or
``get_radius()`` is not.

Variables (member data) should be named using all small case letters or numbers
-------------------------------------------------------------------------------

Underscores should be used to separate
words, and obvious abbreviations are recommended (e.g. ``mol`` for ``molecule``).
For example, ``added_mol`` is acceptable, but ``added_molecule`` should be avoided,
and ``Added_Mol``, ``addedMol``, ``Addedmol`` are all not acceptable

Exceptions are named in the same way as variables
-------------------------------------------------

except abbreviations should not be used,
e.g. ``missing_molecule`` is acceptable, but ``missing_mol``
or ``Missing_Molecule`` or ``MissingMolecule`` is not.

No line should be over 90 characters long
-----------------------------------------

Long lines should be split,
with the extra part indented so that it lines up with the above line, e.g.

.. code-block:: c++

   AtomCoords coords = mol.atom( AtomName("O00") )
                          .property("coordinates")
                          .asA<AtomCoords>();

Always code using a fixed-width font
------------------------------------

The code
uses whitespace and indentation to make things clear, and this is lost
if you use a variable width font

Use whitespace to make the code clean
-------------------------------------

For example, always have a blank
line before a code block (e.g. function, if statement, for loop),
except if it comes directly after an open brace ``{``. For example

.. code-block:: c++

   void foo()
   {
       int a;

       if (a == 5)
       {
           for (int i=0; i<10; ++i)
           {
               if (b == 10)
               {
                   a = 5 * b;

                   for (int j=0; j<11; ++j)
                   {}
               }
           }
       }
   }

Speaking of braces, please use the above style
----------------------------------------------

e.g. braces are on their own
line and line up. This makes it much easier to read.

Sire uses doxygen-style comments
--------------------------------

Sire uses ``Py++`` to extract doxygen-style comments from the
C++ source, which are then added to the Python wrappers,
then extracted by sphinx to create the website.

This means comments
should be written using these rules which are followed at all
times, as the comments are seen in the Python wrappers.

1. All class and function comments should start ``/**`` and end with ``*/``
2. Please do not add ``@author`` information to classes or files,
   as this discourages others from changing your code. Instead, please
   add your name (or GitHub handle) to our :doc:`../contributors` file.
   (note that we are removing existing ``@author`` info, so bear with
   us if you still see any in the code).
3. Use ``//`` for all other comments (even multiline). This is so that it is
   possible to quickly comment out blocks of text using ``/*`` and ``*/``

Finally, keep an eye here as more rules will be written
-------------------------------------------------------

This document will continue to evolve. If you would like to debate
an existing rule or propose a new rule then please
`raise an issue on GitHub <https://github.com/OpenBioSim/sire/issues>`__.
