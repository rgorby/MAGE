
Install HDF5 3rd Party Compression Filters
==========================================

In recent builds of HDF5 (> 1.13.3), HDF Group has included a set of 3rd party plugins with HDF5 distributions. You can also install these `plugins <https://github.com/HDFGroup/hdf5_plugins>`_ separately.

If you are using python tools to read compressed hdf5 data (kaipy, h5py). You can use those libraries without modification if your files were compressed with SZIP or ZLIB. For any 3rd party compression plugins, you can simpley import `hdf5plugin <https://hdf5plugin.readthedocs.io/en/stable/>`_ after your h5py import.
