#!/usr/bin/env python
#Takes Gamera/Chimp/xxx file and slims it down based on start:stop:stride
#Removes all that janky tecplot stuff that got into the regular slim script

import argparse
import os
import h5py
import numpy as np

def genMPIStr(di,dj,dk,i,j,k,n_pad=4):
    inpList = [di, dj, dk, i, j, k]
    sList = ["{:0>{n}d}".format(s, n=n_pad) for s in inpList]
    mpiStr = '_'.join(sList)
    return mpiStr

def cntSteps(fname):
    with h5py.File(fname,'r') as hf:
        """
        grps = hf.values()
        grpNames = [str(grp.name) for grp in grps]
        #Steps = [stp if "/Step#" in stp for stp in grpNames]
        Steps = [stp for stp in grpNames if "/Step#" in stp]
        nSteps = len(Steps)
        """
        #sIds = np.array([str.split(s,"#")[-1] for s in Steps],dtype=np.int)
        sIds = np.array([str.split(s,"#")[-1] for s in hf.keys() if "Step#" in s],dtype=np.int)
        nSteps = len(sIds)
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
    ns = -1
    ne = -1 #Proxy for last entry

    doMPI = False
    parser = argparse.ArgumentParser(description="Slims down an HDF output file")

    parser.add_argument('-intag' ,metavar='intag' ,default="msphere",help="Run ID of input fat HDF file (default: %(default)s)")
    parser.add_argument('-outtag',metavar='outtag',default="slim",help="Run ID of output slimmed HDF file (default: %(default)s)")
    parser.add_argument('-type'  ,metavar='type',default="gam",help="Model type (default: %(default)s)")
    parser.add_argument('-s',type=int,metavar="start",default=ns,help="Starting slice")
    parser.add_argument('-e',type=int,metavar="end",default=-1,help="Ending slice (default: N)")
    parser.add_argument('-sk',type=int,metavar="nsk",default=1,help="Stride (default: %(default)s)")
    parser.add_argument('-mpi',type=str,metavar="ijk",default="", help="Comma-separated mpi dimensions (example: '4,4,1', default: noMPI)")
    
    #Finalize parsing
    args = parser.parse_args()

    Ns = args.s
    Ne = args.e
    Nsk = args.sk
    inTag = args.intag
    outTag = args.outtag
    mType = args.type
    mpiIn = args.mpi

    if ( (Ns == -1) or (Ne == -1) ):
        fIn = "%s.volt.h5"%(inTag)    
        N,sIds = cntSteps(fIn)

        if (Ns == -1):
            Ns = np.sort(sIds)[0]
        if (Ne == -1):
            Ne = N+np.sort(sIds)[0]

    #Designed for 3-dim gamera mpi decomp
    if mpiIn != "":
        doMPI = True
        spl = [int(x) for x in mpiIn.split(',')]
        if len(spl) != 3:
            print("Need 3 dimensions for MPI decomp, try again")
            quit()
        mi, mj, mk = spl

        inFiles = []
        outFiles = []
        for i in range(mi):
            for j in range(mj):
                for k in range(mk):
                    mpiStr = genMPIStr(mi,mj,mk,i,j,k)
                    #fName = runTag+"_"+mpiStr+'.'+endTag
                    fIn  = inTag  + "_" + mpiStr + ".%s.h5"%(mType)
                    fOut = outTag + "_" + mpiStr + ".%s.h5"%(mType)
                    print("%s to %s"%(fIn,fOut))
                    if os.path.exists(fIn):
                        inFiles.append(fIn)
                        outFiles.append(fOut)


    else:
        inFiles  = [inTag  + '.' + str(mType) + '.h5']
        outFiles = [outTag + '.' + str(mType) + '.h5']


    for i in range(len(inFiles)):
        fOut = outFiles[i]
        print(fOut)

        #Open both files, get to work
        iH5 = h5py.File(inFiles[i],'r')
        oH5=createfile(iH5,outFiles[i])

        #Now loop through steps and do same thing
        nOut = 0
        for n in range(Ns,Ne,Nsk):
            
            gIn  = "Step#%d"%(n)
            gOut = "Step#%d"%(nOut)
            nOut = nOut + 1
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

        #Close up
        iH5.close()
        oH5.close()


