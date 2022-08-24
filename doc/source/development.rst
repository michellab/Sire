=================
Developer's guide
=================

The source code for Sire is available on
`GitHub <https://github.com/michellab/Sire>`__.

Setting up your computer
=========================

Sire will install everything it needs when you run the
`compile_sire.sh` script in the top-level Sire directory.

This will automatically download miniconda, and then use that
to install all of Sire's dependencies (including all of
BioSimSpace's dependencies). These include compilers that are
packaged with miniconda. The script will then call cmake
to first compile the core library (corelib), and then all
of the python wrappers. Compilation can take a very long time!

Virtual environments
--------------------

It is recommended that you develop Sire in the miniconda environment
that the `compile_sire.sh` script creates automatically. This will,
by default, be in `$HOME/sire.app`.

It is worth activating this environment during development and testing,
e.g. via

.. code-block:: bash

   $ source $HOME/sire.app/bin/activate

This will update your shell so that all python commands (such as
``python``, ``pip`` etc.) will use the virtual environment. You can
deactivate the environment and return to the "standard" Python using;

.. code-block:: bash

   $ conda deactivate

If you no longer want the environment then you can remove it using

.. code-block:: bash

  rm -rf $HOME/sire.app

Coding Style
============

Sire is written predominantly in C++. This was for speed and memory
efficiency. The Sire C++ objects are wrapped into Python using Py++.

C++ code style is used for names, with the code written to be strictly
C++ 2014 conformant (although we welcome requests to move to a newer
C++ standard, if this is justifiable). The code is very portable and
should remain so. We ourselves are running production code on X86-64
and ARM64 processors, on Linux and MacOS. We know of people who have
used Sire on PowerPC and Windows. We are keen to add first-party
support for a Windows Sire build (when we have time!).

We are in the process of (slowly) improving the ease of use of Sire
by adding a direct Python layer. We aim as much as possible
in this Python layer to follow a
`PEP8 <https://www.python.org/dev/peps/pep-0008/>`__ python coding style and
recommend that developers install and use
a linter such as `flake8 <https://flake8.pycqa.org/en/latest/>`__.

For ease of installation and support, we require that dependencies are
available in stock Anaconda Python, or on conda-forge. As a last resort,
we will bundle dependencies, although this does make our conda build
more complex.

With this in mind, we use the following coding conventions:

Naming
------

We follow a Python style naming convention.

* Packages: Uppercase, singleword (note we are exploring moving to a lowercase style)
* Classes: CamelCase
* Methods: snake_case for pure Python, lowerCamelCase for C++
* Functions: snake_case for pure Python, lowerCamelCase for C++
* Variables: snake_case for pure Python, should not be used in C++ (variables should be private)
* Source Files: snake_case with a leading underscore for pure Python, lowerclassname.cpp / lowerclassname.h for C++

Functions or variables in Python that are private should be named with a leading
underscore. This prevents them from being prominantly visible in Python's
help and tab completion. Any C++ code should only use private variables,
and should use private or protected as much as possible to reduce the API
of C++ classes.

Workflow
========

Feature branches
----------------

First make sure that you are on the development branch of Sire:

.. code-block:: bash

   git checkout devel

Now create and switch to a feature branch. This should be prefixed with
*feat*, e.g.

.. code-block:: bash

   git checkout -b feat-process

Testing
=======

When working on your feature it is important to write tests to ensure that it
does what is expected and doesn't break any existing functionality. Tests
should be placed inside the separate `SireUnitTests <https://github.com/michellab/SireUnitTests>`__
repository, creating an appropriately
named sub-directory for any new modules. Add the tests together with
a guard so that they can detect if your new code is available, and will
be skipped if run on an older version of Sire.

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

To run the full test suite, go to the `SireUnitTests/unittest`
directory and type:

.. code-block:: bash

   pytest .

To run tests for a specific sub-module, e.g. SireMol:

.. code-block:: bash

   pytest SireMol

To only run the unit tests in a particular file, e.g.:

.. code-block:: bash

   pytest SireMol/test_atomselection.py

To get more detailed information about each test, run pytests using the
*verbose* flag, e.g.:

.. code-block:: bash

   pytest -v

More details regarding how to invoke ``pytest`` can be
found `here <https://docs.pytest.org/en/latest/usage.html>`__.

Writing tests
^^^^^^^^^^^^^

Basics
""""""

Try to keep individual unit tests short and clear. Aim to test one thing, and
test it well. Where possible, try to minimise the use of ``assert`` statements
within a unit test. Since the test will return on the first failed assertion,
additional contextual information may be lost.

Floating point comparisons
""""""""""""""""""""""""""

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
""""""""""""""

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
"""""""""""""""""""

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
""""""""""""""""""

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
"""""""""""""""""

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

Continuous integration and delivery
-----------------------------------

We use GitHub Actions to run a full continuous integration (CI)
on all pull requests to devel and
main, and all pushes to devel and main. We will not merge a pull
request until all tests pass. We only accept pull requests to devel.
We only allow pull requests from devel to main. In addition to CI,
we also perform a build of the website on pushes to devel and tags
to main. Finally, we have set up
continuous delivery (CD) on pushes to main and devel, which
build the conda packages and website.

Documentation
=============

Sire is fully documented using a combination of hand-written files
(in the ``doc`` folder) and auto-generated api documentation created from
`NumPy <https://numpy.org>`__ style docstrings.
See `here <https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard>`__
for details. The documentation is automatically built using
`Sphinx <http://sphinx-doc.org>`__ whenever a commit is pushed to devel
that contains the tag "SOMETHING".

To build the documentation locally you will first need to install some
additional packages.

.. code-block:: bash

   pip install sphinx sphinx_issues sphinx_rtd_theme

Then move to the ``doc`` directory and run:

.. code-block:: bash

   make html

When finished, point your browser to ``build/html/index.html``.

Committing
==========

If you create new tests, please make sure that they pass locally before
commiting. When happy, commit your changes, e.g.

.. code-block:: bash

   git commit -a -m "Implementation and test for new feature."

Remember that it is better to make small changes and commit frequently.

If your edits don't change the Sire source code, or documentation,
e.g. fixing typos, then please add ``ci skip`` to your commit message, e.g.

.. code-block:: bash

   git commit -a -m "Updating docs [ci skip]"

This will avoid unnecessarily running the
`GitHub Actions <https://github.com/metawards/MetaWards/actions>`__, e.g.
building a new Sire package, updating the website, etc.
(the GitHub actions are configured in the file
``.github/workflows/main.yaml``). To this end, we
have provided a git hook that will append ``[ci skip]`` if the commit only
modifies files in a blacklist that is specified in the file ``.ciignore``
(analagous to the ``.gitignore`` used to ignore untracked files). To enable
the hook, simply copy it into the ``.git/hooks`` directory:

.. code-block:: bash

    cp git_hooks/commit-msg .git/hooks

Any additional files or paths that shouldn't trigger a re-build can be added
to the ``.ciignore`` file.

Next, push your changes to the remote server, e.g.

.. code-block:: bash

   # Push to the feature branch on the main MetaWards repo, if you have access.
   git push origin feature

   # Push to the feature branch your own fork.
   git push fork feature

When the feature is complete, create a *pull request* on GitHub so that the
changes can be merged back into the development branch.
For information, see the documentation
`here <https://help.github.com/articles/about-pull-requests>`__.

Thanks
======

First, thanks to you for your interest in Sire and for reading this
far. We hope you enjoy having a play with the code and having a go
at adding new functionality, fixing bugs, writing docs etc.

We would also like to thank Lester Hedges and the
`BioSimSpace <https://biosimspace.org>`__ team who provided great advice
to set up the above, and from whose
`GitHub repo <https://github.com/michellab/biosimspace>`__
most of the procedures, scripts and documentation above is derived.
