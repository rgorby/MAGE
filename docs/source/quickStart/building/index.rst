
Building the ``kaiju`` software
===============================

Introduction
------------

This section describes how to build the ``kaiju`` software on two different
supercomputers - ``derecho`` and ``pleiades``. If you are trying to build the
``kaiju`` software on a different system, use these instructions as a starting
point.

In general, building the ``kaiju`` software requires several steps:

    1. **Preparing your software environment.**

    This includes making sure you have the required compilers and external
    libraries, and making sure all environment variables are set properly. On
    HPC systems such as ``derecho`` and ``pleiades``, this also means running
    the appropriate ``module`` commands.

    2. **Build missing prerequisite libraries.**

    This step is required because not all of the prerequisite libraries are
    provided on all platforms. Typically, this problem is limited to the
    `NetCDF library <https://www.unidata.ucar.edu/software/netcdf/>`_ (which
    is needed by the SpacePy module used in ``kaipy`` code which accesses
    satellite data from `CDAWeb <https://cdaweb.gsfc.nasa.gov/>`_). If you do
    not plan to use this capability, then the NetCDF build steps may be
    ignored.

    The `GEOS <https://libgeos.org/>`_ library (which is used by the CartoPy
    module in some ``kaipy`` code which requires mapping functions) is
    available on ``derecho`` as a ``module``, but not on ``pleiades``. We do
    not provide instructions for building the GEOS library. Therefore,
    postprocessing code which uses CartoPy will not run on ``pleiades``.

    3. **Create a python environment.**

    The ``kaiju`` software, and its accompanying ``kaipy`` python package, are
    best used in a separate ``python`` virtual environment. These instructions
    assume the use of the ``conda`` tool for creating and managing virtual
    environments. If you use a different tool, such as the standard ``venv``
    ``python`` module, please adjust these instructions accordingly. The only
    significant differences will be in the commands used to create and
    activate the virtual environments. These instructions should work with
    personal installations of ``python`` (created with the
    `Miniconda <https://docs.anaconda.com/miniconda/>`_ or
    `Anaconda <https://anaconda.org/>`_ ``python`` distributions), or with a
    system-wide installation (such as a ``python`` installation made available
    via the ``module`` command).

    .. important::

        The ``kaiju`` and ``kaipy`` code were developed using ``python`` 3.8,
        which is no longer supported. If you try to run the ``kaiju`` and
        ``kaipy`` code using a different version of ``python``, you will
        probably encounter compatibility issues. We are working to update the
        ``kaiju`` and ``kaipy`` code to support at least ``python`` 3.10.

    4. **Build the kaiju software.**

    This is the step where you actually compile the ``kaiju`` code, which is
    written primarily in `Fortran <https://fortran-lang.org/>`_. The ``kaiju``
    code was developed primarily with the
    `Intel Fortran compiler <https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html>`_,
    but it can also be built using the
    `GNU compiler suite <https://gcc.gnu.org/>`_.

    The source code for the ``kaiju`` software is available on BitBucket. You
    must clone the ``kaiju`` repository before building the code.

    .. important::

        The ``kaiju`` repository on BitBucket uses ``git-lfs`` to support the
        use of large binary files. You *must* make sure ``git-lfs`` is
        available in your ``git`` installation to ensure a complete clone
        of the ``kaiju`` repository.

    The ``kaiju`` software is built using a 2-step process. In the first step,
    the ``cmake`` tool is used to generate the ``Makefile`` for the ``kaiju``
    software. In the second step, the ``make`` tool is used to perform the
    actual build process. The build is *not* performed in the source code
    directory - a separate build directory is always created for these steps.

Table of Contents
-----------------

.. toctree::
    :maxdepth: 1

    buildDerecho
    buildPleiades
