
Building a ``conda``-installable distribution of ``kaipy``
==========================================================

Introduction
------------

This page will walk you through the process of creating a packaged version of
the ``kaipy`` software that can be installed with the ``conda`` Python package
management tool. The ``conda`` tool can install packages from local package
files, as well as from web-based repositories such as the standard Anaconda
repositories.

Setup
-----

The standard python distribution requires additional modules to prepare
``conda``-installable packages. You can install these modules into your
current python environment using the ``conda`` package manager:

.. code-block:: bash

   conda install conda-build
   conda install anaconda-client

Building the package
--------------------

The process of building the ``conda``-installable packages is very simple -
but it assumes you have already performed a ``pip``\ -based package build and
upload to PyPI. This approach allows maintenance of both package flavors with
a minimum of effort. If needed, this procedure can be updated in the future to
remove the ``pip``\ -package prerequisite.

Once you have cloned the ``kaiju`` repository, move into the repository
directory and run the following commands to build the package:

.. code-block:: bash

   cd kaiju
   conda build kaipy

This command reads the file ``kaiju/meta.yaml``\ , and performs all of the
steps required to build a ``conda``\ -installable package. The package files,
by default, are created in a subdirectory of your ``conda``\ -based ``python``
installation. For example, if your ``python`` is installed under:

.. code-block:: bash

   $HOME/miniconda3

you would find the package files in a path similar to:

.. code-block:: bash

   $HOME/miniconda3/envs/PYTHON_ENVIRONMENT_NAME/conda-bld/osx-64/kaipy-0.0.2-py38_0.tar.bz2

where ``PYTHON_ENVIRONMENT_NAME`` is the name of the ``conda``-based
``python`` environment that you used to run the ``conda build`` command. Note
that in this directory you may see more than one file, since previous builds
are not removed. Also, the version numbers you see in the filenames will
increase with time.

Testing the package locally
---------------------------

Once the package has been created locally, perform a test installation into
your local ``python`` environment. Assuming you have run the build described
in the previous step, you can perform the local installation with the command:

.. code-block:: bash

   conda install --use-local $HOME/miniconda3/envs/kaiju-3.8/conda-bld/osx-64/kaipy-0.0.2-py38_0.tar.bz2

where ``kaiju-3.8`` should be replaced with the name of your ``conda``-based
``python`` environment.

Now launch your ``python`` and load the ``kaipy`` package:

.. code-block:: bash

   python -c "import kaipy"

If this command produces no output, it worked - you have installed ``kaipy``
from your local ``pip`` package.

Uploading the package to the Anaconda repository
------------------------------------------------

Once you have built and tested the package locally, you can upload it to the
standard Anaconda repository so that others can easily install ``kaipy``
using ``conda``. Uploading to the Anaconda repository requires the use of a
Anaconda account and password.

.. code-block:: bash

   anaconda login
   # Enter account email and password when prompted.
   anaconda upload $HOME/miniconda3/envs/kaiju-3.8/conda-bld/osx-64/kaipy-0.0.2-py38_0.tar.bz2
   anaconda logout

where ``kaiju-3.8`` should be replaced with the name of your ``conda``-based
``python`` environment.

*IMPORTANT NOTE*\ : Currently, the ``kaipy`` package on Anaconda is managed by
the ``cgsapl`` account, which is attached to the email address
``eric.winter@jhuapl.edu``. If you need to upload a ``conda``-installable
package, please see `Eric Winter <mailto:eric.winter@jhuapl.edu>`_ for the CGS
Anaconda account credentials.

Testing the package from the repository
---------------------------------------

Once the package has been uploaded to the Anaconda repository, perform a test
installation into your local ``conda``\ -based ``python`` environment.
Assuming you have run the steps described above, you can perform the
installation from the Anaconda repository with the command:

.. code-block:: bash

   # Uninstall any previous version for now.
   conda uninstall kaipy
   conda install -c cgsapl kaipy

*NOTE*\ : The ``-c cgsapl`` indicates that the package will be installed from
our "channel" (\ ``cgsapl``\ ) on the Anaconda repository.

Now launch your ``python`` and load the ``kaipy`` package:

.. code-block:: bash

   python -c "import kaipy"

If this command produces no output, it worked - you have installed ``kaipy``
from the Anaconda repository.
