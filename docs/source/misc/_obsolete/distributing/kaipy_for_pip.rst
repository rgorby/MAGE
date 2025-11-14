
Building a ``pip``-installable distribution of ``kaipy``
========================================================

Introduction
------------

This page will walk you through the process of creating a packaged version of
the ``kaipy`` software that can be installed with the ``pip`` Python package
management tool. The ``pip`` tool can install packages from local package
files, as well as from web-based repositories such as
`PyPI (the Python Package Index) <https://pypi.org>`_.

Setup
-----

The standard python distribution already contains the tools needed to build
``pip``-installable packages, but an additional module
(`twine <https://pypi.org/project/twine/>`_) is needed to permit the package
to be uploaded to PyPI so it can be installed by any user.

You can install ``twine`` into your current python environment using the
``pip`` package manager:

.. code-block:: bash

   pip install twine

If you are using a ``conda``-based python environment, you can install
``twine`` using the ``conda`` package manager:

.. code-block:: bash

   conda install twine

Building the package
--------------------

The process of building the ``pip``\ -installable packages is very simple.
Once you have cloned the ``kaiju`` repository, move into the repository
directory and run the ``setup.py`` script like this:

.. code-block:: bash

   cd kaiju
   python setup.py sdist

This will create a subdirectory called ``dist``, containing an installable
tarball for the ``kaipy`` package:

.. code-block:: bash

   (kaiju-3.8) C02FC1QCMD6T-ML:kaiju winteel1$ ls dist
   kaipy-0.0.1.tar.gz      kaipy-0.0.2.tar.gz

Note that you may see more than one file, since previous builds are not
removed. Also, the version numbers you see in the filenames will increase with
time.

Testing the package locally
---------------------------

Once the package has been created locally, perform a test installation into
your local python environment. Assuming you have run the build described in
the previous step, you can perform the local installation with the command:

.. code-block:: bash

   pip install dist/kaipy-0.0.2.tar.gz

Now launch your ``python`` and load the ``kaipy`` package:

.. code-block:: bash

   python -c "import kaipy"

If this command produces no output, it worked - you have installed ``kaipy``
from your local ``pip`` package.

Uploading the package
---------------------

Once you have built and tested the package locally, you can upload it to PyPI
so that others can easily install ``kaipy`` using ``pip``. Uploading to PyPI
requires the use of a PyPI account and password.

.. code-block:: bash

   python -m twine upload dist/*
   # Enter account email and password when prompted.

*IMPORTANT NOTE*\ : Currently, the ``kaipy`` package on PyPI is managed by
the "Center for Geospace Storms" account, which is attached to the email
address ``eric.winter@jhuapl.edu``. If you need to upload a
``pip``-installable package, please see
`Eric Winter <mailto:eric.winter@jhuapl.edu>`_ for the CGS PyPI account
credentials.

Testing the package from the repository
---------------------------------------

Once the package has been uploaded to PyPI, perform a test installation into
your local python environment. Assuming you have run the steps described
above, you can perform the installation from PyPI with the command:

.. code-block:: bash

   # Uninstall any previous version for now.
   pip uninstall kaipy
   pip install kaipy

Now launch your ``python`` and load the ``kaipy`` package:

.. code-block:: bash

   python -c "import kaipy"

If this command produces no output, it worked - you have installed ``kaipy``
from PyPI.
