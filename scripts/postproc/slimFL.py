#!/usr/bin/env python
#Takes H5 field line file and slims it down

import argparse
import os
import kaipy.kaiH5 as kh5
import h5py
import numpy as np

#Create new file w/ same root vars/attributes as old
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
    parser = argparse.ArgumentParser(description="Slims down a FL file")

    parser.add_argument('inH5',metavar='Fat.h5',help="Filename of input fat HDF5 file")
    parser.add_argument('outH5',metavar='Slim.h5',help="Filename of slimmed HDF5 file")
    parser.add_argument( '-pskip',metavar='pskip',default=1,help="Stride over points on field line")
    parser.add_argument('-flskip',metavar='flskip',default=1,help="Stride over field lines")

    #Finalize parsing
    args = parser.parse_args()

    fIn = args.inH5
    fOut = args.outH5
    pSk = np.int(args.pskip)
    lSk = np.int(args.flskip)

    kh5.CheckOrDie(fIn)

    N,sIDs = kh5.cntSteps(fIn)
    
    nS = sIDs.min()
    nE = sIDs.max()

    #For now assuming same number of lines per step
    gID = "Step#%d"%(nS)

    NumL,FLidS = kh5.cntX(fIn,gID,"Line#")
    
    sFL = FLidS.min()
    eFL = FLidS.max()
    LID = "Line#%d"%(sFL)
    print(sFL,eFL)

    #Now open both files and get to work
    #Open both files, get to work
    iH5 = h5py.File(fIn,'r')
    oH5=createfile(iH5,fOut)

    #Get variables from first line
    vIDs = [str(k) for k in iH5[gID][LID].keys()]
    #Remove weird variables to do manually
    vIDs.remove("xyz")
    vIDs.remove("LCon")
    

    for nStp in range(nS,nE+1):
        print(nStp)
        gStr = "Step#%d"%(nStp)
        #Start by copying attributes from old to new
        oH5.create_group(gStr)
        #Root atts
        for k in iH5[gStr].attrs.keys():
            aStr = str(k)
            oH5[gStr].attrs.create(k,iH5[gStr].attrs[aStr])

        #Now loop over field lines
        nOut = 0
        for nFL in range(sFL,eFL+1,lSk):
            #print(nFL)
            iLine = "Line#%d"%(nFL)
            oLine = "Line#%d"%(nOut)

            #Copy iLine => oLine w/ new stride
            oH5[gStr].create_group(oLine)
            for vID in vIDs:
                Q = iH5[gStr][iLine][vID]
                oH5[gStr][oLine].create_dataset(vID,data=Q[::pSk])
            #Get points on line
            xyz0 = iH5[gStr][iLine]["xyz"]
            xyzN = xyz0[::pSk,:]
            oH5[gStr][oLine].create_dataset("xyz",data=xyzN)
            #Now create connectivity
            NumP = np.int32(xyzN.shape[0])

            LCon = np.zeros((NumP-1,2),dtype=np.int32)
            LCon[:,0] = np.arange(0,NumP-1)
            LCon[:,1] = np.arange(1,NumP)
            oH5[gStr][oLine].create_dataset("LCon",data=LCon)
            #Finish up w/ attributes
            aIDs = [str(k) for k in iH5[gStr][iLine].attrs.keys()]
            aIDs.remove("Np")
            aIDs.remove("n0")
            for aID in aIDs:
                oH5[gStr][oLine].attrs.create(aID,iH5[gStr][iLine].attrs[aID])
            oH5[gStr][oLine].attrs.create("Np",NumP)

            nOut = nOut+1
    #Close up
    iH5.close()
    oH5.close()

    	