Sire release instructions
*************************

The following instructions describe how to create a new Sire release:

Step 1
======

Update the `CHANGELOG <https://github.com/michellab/Sire/blob/devel/CHANGELOG>`_
file with a summary of the changes for this relase. Feel free to link to
`GitHub issues <https://github.com/michellab/Sire/issues>`_ where relevant
and give credit for specific contributions.

Step 2
======

When you're happy, tag the commit that you want to be associated with the
release. The following with tag the latest commit:

.. code-block:: bash

    git tag -a 2019.1.0 -m "Tagging the 2019.1.0 release of Sire."

Step 3
======

Push the commit and tag to the ``devel`` branch on the remote:

.. code-block:: bash

    git push origin devel --follow-tags

This will trigger a new Azure Pipelines build which will create binaries
and Conda packages for the release. If you make a mistake and want to move
the tag to a later commit, simply delete the tag from the remote:

.. code-block:: bash

    git push origin :refs/tags/2019.1.0

Next, delete the Conda release package from the `Anaconda Cloud <https://anaconda.org/michellab/sire/files>`_.
You can then move the tag to the latest commit:

.. code-block:: bash

    git tag -fa 2019.1.0

Finally, push the new commit and updated tag:

.. code-block:: bash

    git push origin devel --follow-tags

Step 4
======

Once the build has finished you can log into the `Oracle Cloud <https://cloud.oracle.com/home>`__
and create pre-authenticated URLs for the new release binaries. (These files
are currently located in the ``software_releases`` compartment of the ``Object Storage``
menu. When generating a download URL make sure to choose a sensible expiry
date. Copy the URL to your clipboard and then paste it into a new file called
``redirect_sire_YYYY_MAJOR_MINOR_OS.run`` in the ``largefiles/sire_releases``
directory of the `siremol.org <https://github.com/chryswoods/siremol.org/tree/master/largefiles/sire_releases>`_
GitHub repository. Following this you'll also need to update the
`pages/binaries.md <https://github.com/chryswoods/siremol.org/blob/master/pages/binaries.md>`_
file with download links for the new release. Make sure to add a MD5 hashes for
the files and to move the download links for the previous release to the
``Older Versions`` section at the bottom of the page.

Step 5
======

Create a `GitHub release <https://github.com/michellab/Sire/releases>`_. When
drafting the release, simpy choose the tag that you have created.

Step 6
======

Create a `pull request <https://github.com/michellab/Sire/pulls>`_ to merge
``devel`` into the ``master`` branch.

That's it!
