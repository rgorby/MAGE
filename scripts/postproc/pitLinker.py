# Creates a single file that links to a series of slice / chop h5 files generated using Parallel in Time feature

import h5py as h5
import argparse
import os


if __name__=="__main__":

    fTag = 'ebSlc.eb'
    outfname = ""
    numEb = 8

    MainS = """Takes a series of Parallel in Time - generated eb files and makes a new h5 file linking to them
    """

    parser = argparse.ArgumentParser(description="Generates XDMF file from Gamera HDF5 output")
    parser.add_argument('-outname',type=str,default=outfname,help="Name of generated h5 file")
    parser.add_argument('-id',type=str,default=fTag,help=" Input Run ID (default: %(default)s)")
    parser.add_argument('-numEb',type=int,default=numEb,help=" Number of EB files to join (default: %(default)s)")
    args = parser.parse_args()

    fTag = args.id
    outfname = args.outname
    numEb = args.numEb


    """
    myfile = h5py.File('foo.hdf5','w')
    myfile['ext link'] = h5py.ExternalLink("otherfile.hdf5", "/path/to/resource")
    """

