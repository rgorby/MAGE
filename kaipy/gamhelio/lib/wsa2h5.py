#!/usr/bin/env python

# vgm, Sep 2019: Script to convert WSA fits files into simple h5 files

from astropy.io import fits
import h5py
import glob

wsaFiles = glob.glob('./*.fits')

for wsaFile in wsaFiles:
    h5file = wsaFile[:-4]+'h5'
    
    with fits.open(wsaFile) as hdul:
        with h5py.File(h5file,'w') as hf:
            # take care of the comments
            # removing the last ones since we're incorporating them into dataset attributes directly
            for key in hdul[0].header[:-4]:
                # stores header data as attributes of the file
                hf.attrs[key] = hdul[0].header[key]
                hf.attrs[key+" COMMENT"] = str(hdul[0].header.comments[key])

            # take care of the datasets
            br = hf.create_dataset("br",data=hdul[0].data[0])
            br.attrs['Units'] = 'nT'
            br.attrs['Comment'] = 'Coronal field at RADOUT'

            vr = hf.create_dataset("vr",data=hdul[0].data[1])
            vr.attrs['Units'] = 'km/s'
            vr.attrs['Comment'] = 'Radial velocity at RADOUT'

