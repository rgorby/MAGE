#!/usr/bin/env python
#Simple script to spit out the step information of a Kaiju H5 file
import argparse
from argparse import RawTextHelpFormatter
import os
import kaiH5 as kh5
import numpy as np

if __name__ == "__main__":
    #Defaults

    MainS = """Identifies the domain (in steps and time) of a Kaiju HDF-5 file"""
    
    parser = argparse.ArgumentParser(description=MainS, formatter_class=RawTextHelpFormatter)
    parser.add_argument('h5F',nargs='+',metavar='Gamera.h5',help="Filename of Gamera HDF5 Output")

    #Finished getting arguments, parse and move on
    args = parser.parse_args()
    #h5F = args.h5F
    for idx, h5F in enumerate(args.h5F):
        print("Reading %s"%(h5F))
        nSteps,sIds = kh5.cntSteps(h5F)
        s0 = sIds.min()
        sE = sIds.max()
        print("\tFound %d steps"%(nSteps))
        print("\tSteps = [%d,%d]"%(s0,sE))
        tMin = kh5.getTs(h5F,np.array([s0]))
        tMax = kh5.getTs(h5F,np.array([sE]))
        print("\tTime = [%f,%f]"%(tMin,tMax))
    #---------------------
