# Project Kaiju
Reorganization of Kaiju models


# Example build

mkdir build
cd build
cmake ..

To force GNU, prior to cmake do
export FC="gfortran"

Can also do "ccmake .." for semi-graphical interface.
"cmake -L" to list various options

# Cheyenne configuration
Cheyenne is somewhat sensitive to the order of modules.  This is the configuration I've used (3/20/19)

Currently Loaded Modules:
  1) intel/18.0.1      4) python/2.7.13   7) h5py/2.7.0    10) ncarenv/1.2
  2) hdf5/1.8.18       5) numpy/1.13.3    8) git/2.10.2    11) ncarcompilers/0.4.1
  3) impi/2018.1.163   6) scipy/0.19.1    9) cmake/3.12.1