Sire build instructions
***********************

You only need to build Sire if there isn't a package available for your
computer and operating system on conda, or if you want to make changes
to Sire.

If you want to compile Sire from source, then first download it from
Git, e.g. `git clone https://github.com/michellab/Sire` will download
the (stable) development branch.

Change into this directory (`cd Sire`) to continue with the build.

Dependencies
------------

Sire is built using
`Miniconda <https://docs.conda.io/en/latest/miniconda.html>`__,
`Miniforge <https://github.com/conda-forge/miniforge>`__, or
`Mambaforge <https://github.com/conda-forge/miniforge>`__ creating
an isolated Conda environment into which Sire with be compiled and deployed.
This avoids the need for external dependencies, since everything that is
required is already available within the Conda ecosystem.

The main script used to compile Sire is
`setup.py <https://github.com/michellab/Sire/blob/devel/setup.py>`__.

To compile, activate your conda environment and then install the build
dependencies. These are;

* Mamba - this is an alternative to `conda` that significantly speeds
  up installation of dependencies. Install using `conda install mamba`.

* The conda C and C++ compilers (`mamba install gcc gxx` on Linux,
  `mamba install clang clangxx` on MacOS, and `mamba install conda-build` on Windows).

* Cmake - `mamba install cmake`

* pip-requirements-parser - `mamba install pip-requirements-parser`

Compiling Sire
--------------

With the dependencies installed, you can now compile and install Sire
using the command

::
    python setup.py install


This will download the rest of Sire's dependencies, will compile Sire
and will install it into your current environment.

Additional Options
------------------

There are some additional options that are useful to control compilation,
or recompilation if you are editing Sire.

* `--help` : Show instructions for all additional options
* `--skip-deps` : Skip the testing and installation of conda dependencies.
* `--skip-build` : Skip the build phase - just install the Python wrappers.
* `--install-bss-deps` : Additionally install BioSimSpace's dependencies too.

The available `setup.py` options are;

* `install` : Compile and install Sire in full
* `install_module` : Install only the pure Python module. Useful during
   development if you have only made changes to the pure Python module.
* `build` : Only compile Sire. This will build and install the core C++ library,
  and will build (but not install) the C++ Python wrappers.
* `install_requires` : Install all of Sire's dependencies using `conda` or `mamba`.

Typically, as a developer, if you make any change to `corelib` you should
re-run the compile using `python setup.py --skip-deps install`.

If you only make changes to the wrappers, you can speed things up by
re-running the installation using `python setup.py --skip-deps --skip-build install`.

If you only make changes to the pure Python module, then you can get away
with just running `python setup.py install_module`.

Building the conda package
--------------------------

Alternatively, you can compile Sire using the conda recipe. To do this, first
install the conda build dependencies, including `boa` (supports mamba builds)

::

  mamba install -y -c conda-forge boa anaconda-client pip-requirements-parser


Next, update the recipe by typing

::

    python actions/update_recipe.py

Finally, build the recipe via

::

    conda mambabuild -c conda-forge -c michellab recipes/sire

This set of commands is used by Sire's CI/CD
`GitHub Actions workflow <https://github.com/michellab/Sire/blob/devel/.github/workflows/main.yaml>`__.

We have tested this on Windows, MacOS and Linux, building packages
(currently) for Python 3.7, 3.8 and 3.9.
