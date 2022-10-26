========
Sire.CAS
========

This module implements a complete Computer Algebra System. This is
used to create algebraic expressions that can represent energy
expressions (or anything else you want). This is used in Sire to
hold generic expressions for internal energies (e.g. bond, angle or
dihedral potentials). It is also used to build energy expressions
as a combination of components of forcefields. This lets you build
arbitrary energy expressions that involve user-defined parameters
(such as lambda, for free energy simulations).

The module is built from some core classes.

:class:`~sire.cas.Symbol`
    This represents an algebraic symbol.

:class:`~sire.cas.Expression`
    This represents a complete algebraic expression.

:func:`~sire.cas.create_symbols`
    Create symbols from the passed strings.

There are also lots of functions defined, e.g. `~sire.cas.Cos`,
`~sire.cas.Sin`, `~sire.cas.Coth` etc.

Example
-------

.. code-block:: python

   >>> x, theta = create_symbols("x", "theta")
   >>> f = x**2 + 3 * x - 5
   >>> print(f.differentiate(x))
   2 x + 3

   >>> f = (Sin(theta) + Cos(theta))**3
   >>> print(f.evaluate({theta: 0.5}))
   2.4988910432100395

   >>> print(f.differentiate(theta))
   3 [[cos(theta) - sin(theta)] [sin(theta) + cos(theta)]^2]

.. toctree::
   :maxdepth: 3

   index_api_Sire_CAS
