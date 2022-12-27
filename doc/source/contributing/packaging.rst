==================
Packaging releases
==================

Sire is now fully tested and deployed using GitHub actions.
The development process should be;

* New features are developed on feature branches, called ``feature-{feature}``,
  either in the `main Sire repository <https://github.com/michellab/Sire>`__
  for authorised developers, or in personal forks for
  new developers.
* Bug fixes or issue fixes are developed on fix branches, called
  ``fix-issue-{number}`` (again in either the main repository or forks).
* Pull requests are issued from these branches to ``devel``. All merge conflicts
  must be fixed in the branch and all tests must pass before the pull
  request can be merged into ``devel``. NOTE THAT ONLY AUTHORISED
  DEVELOPERS CAN ACCEPT THE PULL REQUEST. Authorised developers will
  review the pull request as quickly as they can. They will be greatly
  helped if the feature is accompanied with tests, examples and/or tutorial
  instructions.

The result of this is that "devel" should contain the fully-working and
tested, and most up-to-date version of ``Sire``. However, this
version should not be used for production runs.

.. note::

  The group of developers authorised to have access to the
  `main Sire repository <https://github.com/michellab/Sire>`__
  and to accept pull requests is not fixed,
  and will evolve over time. If you wish to join this group then
  please complete the tutorial and then demostrate your commitment
  by submitting good issues and pull requests from
  a personal fork of the repository. Please get in touch if you find
  this difficult, or follow
  `this workshop <https://chryswoods.com/beginning_git>`__ if you need
  to learn how to use Git, GitHub, feature branching, merging, pull
  requests etc.

Defining a release
------------------

We will release ``Sire`` regularly. Releases aim to be backwards
compatible and capable of being used for production runs, at least for
the functionality that is fully described in the tutorial.

.. note::

  It is the job of the release manager (currently
  `lohedges <https://github.com/lohedges>`__) to decide when it is time
  to create a new release. If you are interested in helping join the release
  management group then please feel free to get in touch.

Creating a release
------------------

To create a release accept a pull request from devel into main. This
will automatically trigger the GitHub Action to build the conda packages
for the platforms and python versions that we currently support.
These will be pushed automatically to conda-forge.

Please remember to update the :doc:`changelog` with details about this release. This
should include the link at the top of the release that shows the commit
differences between versions. This can be easily copied from a previous
release and updated.
