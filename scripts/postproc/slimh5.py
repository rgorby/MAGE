#!/usr/bin/env python
#Takes Gamera/Chimp/xxx file and slims it down based on start:stop:stride

import argparse
import os
import h5py
import numpy as np

def cntSteps(fname):
    with h5py.File(fname,'r') as hf:
        grps = hf.values()
        grpNames = [str(grp.name) for grp in grps]
        #Steps = [stp if "/Step#" in stp for stp in grpNames]
        Steps = [stp for stp in grpNames if "/Step#" in stp]
        nSteps = len(Steps)

        sIds = np.array([str.split(s,"#")[-1] for s in Steps],dtype=np.int)
    return nSteps,sIds

def createfile(iH5,fOut):
    print('Creating new output file:',fOut)
    oH5 = h5py.File(fOut,'w')
#Start by scraping all variables from root
    #Copy root attributes
    for k in iH5.attrs.keys():
        aStr = str(k)
        oH5.attrs.create(k,iH5.attrs[aStr])
    #Copy root groups
    for Q in iH5.keys():
        sQ = str(Q)
        #Don't include stuff that starts with "Step"
        if "Step" not in sQ:
            oH5.create_dataset(sQ,data=iH5[sQ])
    return oH5
if __name__ == "__main__":
    #Set defaults
    ns = 0
    ne = -1 #Proxy for last entry
    parser = argparse.ArgumentParser(description="Slims down an HDF output file")

    parser.add_argument('inH5',metavar='Fat.h5',help="Filename of input fat HDF5 file")
    parser.add_argument('outH5',metavar='Slim.h5',help="Filename of slimmed HDF5 file")
    parser.add_argument('-s',type=int,metavar="start",default=ns,help="Starting slice (default: %(default)s)")
    parser.add_argument('-e',type=int,metavar="end",default=-1,help="Ending slice (default: N)")
    parser.add_argument('-sk',type=int,metavar="nsk",default=1,help="Stride (default: %(default)s)")
    parser.add_argument('-sf',type=int,metavar="nsf",default=250,help="File write stride (default: %(default)s)")
    parser.add_argument('--p',choices=('True','False'))
    
    #Finalize parsing
    args = parser.parse_args()

    Ns = args.s
    Ne = args.e
    Nsk = args.sk
    Nsf = args.sf
    fIn = args.inH5
    fOut = str(Ns)+'-'+str(Nsf)+args.outH5
    p = args.p == 'True'
    
    print(p)

    N,sIds = cntSteps(fIn)
    if (Ne == -1):
        Ne = N

    #Open both files, get to work
    iH5 = h5py.File(fIn,'r')
    # oH5 = h5py.File(fOut,'w')
    oH5=createfile(iH5,fOut)


#Now loop through steps and do same thing
    nOut = 0
    for n in range(Ns,Ne,Nsk):
        if(p):nOut = n
        gIn = "Step#%d"%(n)
        gOut = "Step#%d"%(nOut)
        if(not p): nOut = nOut+1 # use the same group numbers as originally in file - frt

        print("Copying %s to %s"%(gIn,gOut))

        oH5.create_group(gOut)
        #Root atts
        for k in iH5[gIn].attrs.keys():
            aStr = str(k)
            oH5[gOut].attrs.create(k,iH5[gIn].attrs[aStr])
        #Root vars
        for Q in iH5[gIn].keys():
            sQ = str(Q)
            #print("\tCopying %s"%(sQ))
            oH5[gOut].create_dataset(sQ,data=iH5[gIn][sQ])
    # make a new file every Nsf steps
        if(n%Nsf==0):
            oH5.close()
            fOut = str(n)+'-'+str(Nsf+n)+args.outH5
            oH5=createfile(iH5,fOut)

    #Close up
    iH5.close()
    oH5.close()


