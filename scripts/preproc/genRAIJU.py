#!/usr/bin/env python	
# Generates RAIJU config data
import os
import h5py as h5
import argparse
from argparse import RawTextHelpFormatter

import kaipy.kaiTools as kT
import kaipy.kaijson as kj
import kaipy.kaiH5 as kh5

import kaipy.raiju.lambdautils.AlamData as aD
import kaipy.raiju.lambdautils.AlamParams as aP
import kaipy.raiju.lambdautils.DistTypes as dT
import kaipy.raiju.lambdautils.fileIO as fileIO

from kaipy.raiju.waveModel.wmData import wmParams
import kaipy.raiju.waveModel.genWM as genWM

#import kaipy.rcm.lambdautils.plotter as plotter


EFLAV = 1
PFLAV = 2

EFUDGE = 1./3.
PFUDGE = 0.0

if __name__ == "__main__":

    #Arg parsing
    fOut = "raijuconfig.h5"
    num_e  = 39
    num_p  = 120
    eminp = 1  # [eV]
    emine = 1  # [eV]
    emaxp = 100 # [keV]
    emaxe = 25  # [keV] , 1/4 of 100 keV
    L_kt   = 10
    wolfP1 = 3
    wolfP2 = 1
    plotChoices = ['none', 'spec', 'vs']


    MainS = """Generates RAIJU configuration data
    """
    parser = argparse.ArgumentParser(description=MainS, formatter_class=RawTextHelpFormatter)
    parser.add_argument('-o',type=str,default=fOut,metavar="fOut",help="Output file name (default: %(default)s)")
    parser.add_argument('-ne', type=int,default=num_e, help="Number of electron channels (default: %(default)s)")
    parser.add_argument('-np', type=int,default=num_p, help="Number of proton channels (default: %(default)s)")
    parser.add_argument('-mine', type=float,default=emine, help="Min. energy [eV] for electrons at L_kt (default: %(default)s)")
    parser.add_argument('-minp', type=float,default=eminp, help="Min. energy [eV] for protons at L_kt (default: %(default)s)")
    parser.add_argument('-maxe', type=float,default=emaxe, help="Max. energy [keV] for electrons at L_kt (default: %(default)s)")
    parser.add_argument('-maxp', type=float,default=emaxp, help="Max. energy [keV] for protons at L_kt (default: %(default)s)")
    parser.add_argument('-L', type=float,default=L_kt, help="L shell [R_e] at which kt should be resolved (default: %(default)s [R_e])")
    parser.add_argument('-p1', type=float,default=wolfP1, help="Wolf low-energy  p* (default: %(default)s)")
    parser.add_argument('-p2', type=float,default=wolfP2, help="Wolf high-energy p* (default: %(default)s)")
    parser.add_argument('-plotType', choices=plotChoices,default=plotChoices[0], help="Plot mode (default: %(default)s)")
    parser.add_argument('--nop',action='store_true',default=False,help="Do not add zero energy first channel (default: %(default)s)")
    parser.add_argument('--noWaveModel',action='store_true',default=False, help="Don't use wave models in the electron/ion loss (default: %(default)s)")
    parser.add_argument('--append',action='store_true',default=False, help="Append existing file rather than make a new one (default: %(default)s)")

    # Finalize parsing
    args = parser.parse_args()
    fOut = args.o
    num_e = args.ne
    num_p = args.np
    emine = args.mine
    eminp = args.minp
    emaxe = args.maxe*1e3  # [keV -> eV]
    emaxp = args.maxp*1e3  # [keV -> eV]
    L_kt = float(args.L)
    wolfP1 = args.p1
    wolfP2 = args.p2
    plotType = args.plotType
    noPsph = args.nop
    noWaveModel = args.noWaveModel
    doAppend = args.append


    if doAppend:
        if not os.path.exists(fOut): 
            # Will only append to a file that exists
            # i.e. we will not make a new file if the --append arg was used
            print("Couldn't find {} to append to, quitting.".format(fOut))
            quit()
        # If append == True and file is there, we'll just add any data the file doesn't already have
    else:
        # Othewise, we will make a new file
        # Since all writing functions open in 'a' mode, we need to ensure we have a fresh file to work with
        print("Making new {}, destroying pre-existing file if there".format(fOut))
        with h5.File(fOut, 'w') as tmpF:
            pass


    # Tag the output file with latest info
    with h5.File(fOut, 'a') as f5:
        print("Stamping file with git hash and branch, and script args")
        kh5.StampHash(fOut)
        kh5.StampBranch(fOut)
        f5.attrs['ScriptArgs'] = kj.dumps(vars(args),noIndent=True)


    # Add wavemodel if desired
    if not noWaveModel:
        tauParams = wmParams(dim = 4, nKp = 7, nMLT = 97, nL = 41, nEk = 155)
        genWM.genh5(fOut,tauParams)


    # Now add species
    # Determine lambda limits using desired energy limits at given L
    bVol = kT.L_to_bVol(L_kt)
    vm = bVol**(-2/3)  # [(Rp/nT)^(-2/3)]
    alamMin_p = eminp/vm  # [eV * (Rp/nT)^(2/3)]
    alamMax_p = emaxp/vm
    alamMin_e = -1*emine/vm
    alamMax_e = -1*emaxe/vm

    dtManual = dT.DT_Single()
    dtWolf   = dT.DT_Wolf(p1=wolfP1,p2=wolfP2)  # Lambda channels will have a (slightly modified) Wolf distribution type

    # For reference: SpecParams(n, amin, amax, distType, 
    #           flav, numNuc_p, numNuc_n, q,
    #           fudge, name)

    if not noPsph: sPsphere = aP.SpecParams(1, 0, 0, dtManual, 
                                            0, numNuc_p=1, numNuc_n=0, q=1,
                                            fudge=0, name="0_Plasmasphere") # Zero-energy plasmasphere channel
    sPe = aP.SpecParams(num_e, alamMin_e, alamMax_e, dtWolf, 
                        EFLAV, numNuc_p=0, numNuc_n=0, q=-1, 
                        fudge=EFUDGE, name='Hot Electrons')  # Parameters to create electron channels
    sPp = aP.SpecParams(num_p, alamMin_p, alamMax_p, dtWolf, 
                        PFLAV, numNuc_p=1, numNuc_n=0, q=1, 
                        fudge=PFUDGE, name='Hot Protons'  )  # Parameters to create proton channels

    # Stuff into an AlamData object, which will automatically generate the fully realized species distribution for us
    # NOTE: ORDER MATTERS. They will show up in the 1D lamc in this order
    if noPsph:
        alamData = aD.AlamData([sPe, sPp])
    else:
        alamData = aD.AlamData([sPsphere, sPe, sPp])

    
    # Save lambda data to file (will also save params as a root attribute)
    fileIO.saveAlamData(fOut, alamData)