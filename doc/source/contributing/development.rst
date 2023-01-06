=================
Developer's guide
=================

The source code for :mod:`sire`` is available on
`GitHub <https://github.com/openbiosim/sire>`__.

Setting up your computer
=========================

You first need to create an anaconda, miniconda or
mambaforce environment as described in
the :doc:`installation page <../install>`.

We recommend using mambaforge, as this sets the right
priority for the conda-forge channel, and it bundles
mamba, which we find to be much faster than conda.

Virtual environments
--------------------

We recommend that you develop :mod:`sire` in its own
conda environment. For example, you could call this
environment ``openbiosim``.

It is worth activating this environment during development and testing,
e.g. via

.. code-block:: bash

   $ source /path/to/environment/bin/activate

(where ``/path/to/environment`` is the file path to your environment)

or by running

.. code-block:: bash

   $ conda activate environment_name

(where ``environment_name`` should be replaced by the name
of your environment - e.g. ``openbiosim``).

This will update your shell so that all python commands (such as
``python``, ``mamba`` etc.) will use the virtual environment. You can
deactivate the environment and return to the "standard" Python using;

.. code-block:: bash

   $ conda deactivate

Python Coding Style
===================

Sire is written predominantly in C++. This was for speed and memory
efficiency. The Sire C++ objects are wrapped into Python using Py++.

The legacy API was very C++, and thus not very pythonic in nature.

We have engaged in a modernisation program, and now (nearly?) all
Python-exposed or Python-native code in the public API is written
in a Pythonic style. We aim to be fully
`PEP8-compliant <https://pep8.org>`__ and ask that all new
Python code contributed to :mod:`sire` is written to be
`PEP8-compliant <https://pep8.org>`__. We plan to move to
enforcing this by asking contributers to use a PEP8
autoformatter
(e.g. `black <https://black.readthedocs.io/en/stable/>`__)
and will add a linting test to our CI/CD pipeline.

C++ Coding Style
================

C++ code style is used for names, with the code written to be strictly
C++ 2014 conformant (although we welcome requests to move to a newer
C++ standard, if this is justifiable). The code is very portable and
should remain so. We ourselves are running production code on X86-64
and ARM64 processors, on Linux, MacOS and Windows.
We know of people who compile and use :mod:`sire` on PowerPC.

We have a strict C++ coding style, which is
:doc:`described here <codestyle>`.

For ease of installation and support, we require that dependencies are
available in `conda-forge <https://conda-forge.org>`__.
As a last resort, we will vendor dependencies,
but this does introduce a significant extra support burden.

Guidelines
==========

With this in mind, we use the following conventions:

* Packages: Lowercase, singleword
* Classes: CamelCase
* Methods: snake_case for pure Python, lowerCamelCase for C++
* Functions: snake_case for pure Python, lowerCamelCase for C++
* Variables: snake_case for pure Python, should not be used in C++ (variables should be private).
  But public Python variables are discouraged. Private variables should
  be preferred, named using a leading underscore.
* Source Files: snake_case with a leading underscore for pure Python, lowerclassname.cpp / lowerclassname.h for C++
* ``__all__`` should be used in Python to expose the public API of
  a file or module. This is used to control what is seen using
  tab completion in ipython / notebooks, and what is extracted
  by sphinx to form the API documentation on the website.
* Documentation - use doxygen style comments for C++ and
  numpy-style documentation for Python. All functions / classes
  in the public API should be documented.

Functions or variables in Python that are private should be named with a leading
underscore. This prevents them from being prominantly visible in Python's
help and tab completion. Any C++ code should only use private variables,
and should use private or protected as much as possible to reduce the API
of C++ classes.

Workflow
========

Feature branches
----------------

First make sure that you are on the development branch of :mod:`sire`

.. code-block:: bash

   git checkout devel

Now create and switch to a feature branch. This should be prefixed with
``feat``, e.g.

.. code-block:: bash

   git checkout -b feat-process

Testing
=======

When working on your feature it is important to write tests to ensure that it
does what is expected and doesn't break any existing functionality. Tests
should be placed inside the ``tests`` directory, and should be designed
to be run using ``pytest``. Note that you should not place any input
files or structure files in the ``tests`` directory. Instead, they
should be placed on the web, and downloaded using :func:`sire.load`
via their URL. When we accept your pull request we will move your
input files onto the main website and will update your test to
download the files from there.

The test suite is intended to be run using
`pytest <https://docs.pytest.org/en/latest/contents.html>`__.
When run, ``pytest`` searches for tests in all directories and files
below the current directory, collects the tests together, then runs
them. Pytest uses name matching to locate the tests. Valid names start
or end with *test*\ , e.g.:

::

   # Files:
   test_file.py       file_test.py

.. code-block:: python

   # Functions:
   def test_func():
      # code to perform tests...
      return

   def func_test():
      # code to perform tests...
      return

We use the convention of ``test_*`` when naming files and functions.

Running tests
-------------

To run the full test suite, simply run ``pytest`` pointing to
the ``tests`` directory, e.g.

.. code-block:: bash

   pytest tests

Tests for each module are in a directory named after
that module.

To run tests for a specific sub-module, e.g. :mod:`sire.mol`
type

.. code-block:: bash

   pytest tests/mol

To only run the unit tests in a particular file,
e.g. ``tests/mol/test_atomprops.py``, you can type

.. code-block:: bash

   pytest tests/mol/test_atomprops.py

To get more detailed information about each test, run pytests using the
*verbose* flag, e.g.:

.. code-block:: bash

   pytest -v tests

More details regarding how to invoke ``pytest`` can be
found `here <https://docs.pytest.org/en/latest/usage.html>`__.

Writing tests
=============

Basics
------

Try to keep individual unit tests clear and fast. The aim is that they
should test a single part of the code, and should complete in seconds
(if not quicker). Use fixtures to re-use files that have been
downloaded and parsed as much as possible. These are all defined
in the file `tests/conftests.py <https://github.com/OpenBioSim/sire/blob/devel/tests/conftest.py>`__.

Floating point comparisons
--------------------------

Make use of the
`approx <https://docs.pytest.org/en/latest/builtin.html#comparing-floating-point-numbers>`__
function from the ``pytest`` package for performing floating
point comparisons, e.g:

.. code-block:: python

   from pytest import approx

   assert 0.1 + 0.2 == approx(0.3)

By default, the ``approx`` function compares the result using a
relative tolerance of 1e-6. This can be changed by passing a keyword
argument to the function, e.g:

.. code-block:: python

   assert 2 + 3 == approx(7, rel=2)

Skipping tests
--------------

If you are using
`test-driven development <https://en.wikipedia.org/wiki/Test-driven_development>`__
it might be desirable to write your tests before implementing the functionality,
i.e. you are asserting what the *output* of a function should be, not how it should
be *implemented*. In this case, you can make use of
the ``pytest`` *skip* decorator
to flag that a unit test should be skipped, e.g.:

.. code-block:: python

   @pytest.mark.skip(reason="Not yet implemented.")
   def test_new_feature():
       # A unit test for an, as yet, unimplemented feature.
       ...

Parametrizing tests
-------------------

Often it is desirable to run a test for a range of different input parameters.
This can be achieved using the ``parametrize`` decorator, e.g.:

.. code-block:: python

   import pytest
   from operator import mul

   @pytest.mark.parametrize("x", [1, 2])
   @pytest.mark.parametrize("y", [3, 4])
   def test_mul(x, y):
       """ Test the mul function. """
       assert mul(x, y) == mul(y, x)

Here the function test_mul is parametrized with two parameters, ``x`` and ``y``.
By marking the test in this manner it will be executed using all possible
parameter pairs ``(x, y)``\ , i.e. ``(1, 3), (1, 4), (2, 3), (2, 4)``.

Alternatively:

.. code-block:: python

   import pytest
   from operator import sub
   @pytest.mark.parametrize("x, y, expected",
                           [(1, 2, -1),
                            (7, 3,  4),
                            (21, 58, -37)])
   def test_sub(x, y, expected):
       """ Test the sub function. """
       assert sub(x, y) == -sub(y, x) == expected

Here we are passing a list containing different parameter sets, with the names
of the parameters matched against the arguments of the test function.

Testing exceptions
------------------

Pytest provides a way of testing your code for known exceptions. For example,
suppose we had a function that raises an ``IndexError``\ :

.. code-block:: python

   def indexError():
       """ A function that raises an IndexError. """
       a = []
       a[3]

We could then write a test to validate that the error is thrown as expected:

.. code-block:: python

   def test_indexError():
       with pytest.raises(IndexError):
           indexError()

Custom attributes
-----------------

It's possible to mark test functions with any attribute you like. For example:

.. code-block:: python

   @pytest.mark.slow
   def test_slow_function():
       """ A unit test that takes a really long time. """
       ...

Here we have marked the test function with the attribute ``slow`` in order to
indicate that it takes a while to run. From the command line it is possible
to run or skip tests with a particular mark.

.. code-block:: bash

   pytest mypkg -m "slow"        # only run the slow tests
   pytest mypkg -m "not slow"    # skip the slow tests

The custom attribute can just be a label, as in this case, or could be your
own function decorator.

Please do use ``slow`` to mark tests that take more than 3 seconds
to run, and use ``veryslow`` for tests that take more than 10 seconds
to run.

Continuous integration and delivery
-----------------------------------

We use GitHub Actions to run a full continuous integration (CI)
on all pull requests to ``devel`` and
``main``, and all pushes to ``devel`` and ``main``. We will not merge a pull
request until all tests pass. We only accept pull requests to ``devel``.
Only the release managers and accept pull requests to ``devel``.

Only the release managers can make and accept pull requests
from ``devel`` to ``main``, and only as part of creating a new
release of :mod:`sire`. In addition to CI,
we also perform a build of the website on pushes to devel and tags
to ``main``. Finally, we have set up
continuous delivery (CD) on pushes to ``main`` and ``devel``, which
build and upload the conda packages.

Documentation
=============

Sire is fully documented using a combination of hand-written files
(in the ``doc`` folder) and auto-generated api documentation created from
`NumPy <https://numpy.org>`__ style docstrings.
See `here <https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard>`__
for details. The documentation is automatically built using
`Sphinx <http://sphinx-doc.org>`__.

To build the documentation locally you will first need to install some
additional packages as described in the
`requirements.txt <https://github.com/OpenBioSim/sire/blob/devel/doc/requirements.txt>`__

.. code-block:: bash

   mamba install sphinx sphinxcontrib-programoutput sphinx_issues furo

Then move to the ``doc`` directory and run:

.. code-block:: bash

   make

When finished, point your browser to ``build/html/index.html``.

Committing
==========

If you create new tests, please make sure that they pass locally before
commiting. Please also check that all your Python code is formatted
to be `PEP8-compliant <https://pep8.org>`__. This will be easier
if you use an autoformatter such as `black <https://black.readthedocs.io/en/stable/>`__.

When happy, commit your changes, e.g.

.. code-block:: bash

   git commit -a -m "Implementation and test for new feature."

Remember that it is better to make small changes and commit frequently.

Next, make sure that you have no conflicts with the ``devel``
branch. Pull this branch via;

.. code-block:: bash

   git pull origin devel

and resolve any conflicts that appear (ideally by modifying
your code). Please feel free to get in touch if there are many
conflicts or you need to modify lots of other code in :mod:`sire`.

Remember to then recompile your code and check that all
of the unit tests (including your new tests) pass.

If your edits don't change the :mod:`sire` source code, or documentation,
e.g. fixing typos, then please add ``ci skip`` to your commit message, e.g.

.. code-block:: bash

   git commit -a -m "Updating docs [ci skip]"

This will avoid unnecessarily running the
`GitHub Actions <https://github.com/OpenBioSim/sire/actions>`__, e.g.
building a new :mod:`sire`` package, updating the website, etc.
(the GitHub actions are configured in the file
``.github/workflows/main.yaml``).

Next, push your changes to the remote server, e.g.

.. code-block:: bash

   git push

When the feature is complete, create a *pull request* on GitHub so that the
changes can be merged back into the development branch.
For information, see the documentation
`here <https://help.github.com/articles/about-pull-requests>`__.

Thanks
======

First, thanks to you for your interest in :mod:`sire`` and for reading this
far. We hope you enjoy having a play with the code and having a go
at adding new functionality, fixing bugs, writing docs etc.

We would also like to thank Lester Hedges and the
`BioSimSpace <https://biosimspace.org>`__ team who provided great advice
to set up the above, and from whose
`GitHub repo <https://github.com/michellab/biosimspace>`__
most of the procedures, scripts and documentation above is derived.
