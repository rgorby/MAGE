#!/usr/bin/env python	
#Generates SIF config data
import numpy as np
import h5py as h5
import argparse
from argparse import RawTextHelpFormatter

import kaipy.kdefs as kdefs
import kaipy.kaiTools as kT
import kaipy.kaijson as kj
import kaipy.kaiH5 as kh5

import kaipy.sif.lambdautils.AlamData as aD
import kaipy.sif.lambdautils.AlamParams as aP
import kaipy.sif.lambdautils.DistTypes as dT
import kaipy.sif.lambdautils.fileIO as fileIO

from kaipy.rcm.wmutils.wmData import wmParams
import kaipy.rcm.wmutils.genWM as genWM

#import kaipy.rcm.lambdautils.plotter as plotter


EFLAV = 1
PFLAV = 2

EFUDGE = 1./3.
PFUDGE = 0.0

if __name__ == "__main__":

    #Arg parsing
    fOut = "sifconfig.h5"
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


    MainS = """Generates SIF configuration data
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
    parser.add_argument('--addWM', action='store_true',default=False, help="Add wave models to an existing rcmconfig file, input file needed to be presented (default: %(default)s)")
    parser.add_argument('-i', type=str,default=fOut,metavar="fIn", help="Input file name when addWM is true (default: %(default)s)")


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
    noPsph = args.nop
    addWM = args.addWM
    noWaveModel = args.noWaveModel
    fIn = args.i	
    plotType = args.plotType


    # Check if we are just here to add the wavemodel to an existing config file
    if addWM:
        tauParams = wmParams(dim = 4, nKp = 7, nMLT = 25, nL = 41, nEk = 155)
        genWM.genh5(fIn,fOut,tauParams,useWMh5 = True)
        quit()

    # Otherwise, procede with generating a file from scratch

    # Determine proton channel limits based on resolving a certain (proton) temperature at given L
    bVol = kT.L_to_bVol(L_kt)
    vm = bVol**(-2/3)
    alamMin_p = eminp/vm
    alamMax_p = emaxp/vm
    alamMin_e = -1*emine/vm
    alamMax_e = -1*emaxe/vm

    dtManual = dT.DT_Single()
    dtWolf   = dT.DT_Wolf(p1=wolfP1,p2=wolfP2)  # Lambda channels will have a (slightly modified) Wolf distribution type

    # Calculate masses in amu
    Me_amu = kdefs.Me_cgs*1e-3/kdefs.dalton
    Mp_amu = kdefs.Mp_cgs*1e-3/kdefs.dalton

    # For reference: SpecParams(n, amin, amax, distType, flav, amu, fudge, name)

    if not noPsph: sPsphere = aP.SpecParams(1, 0, 0, dtManual, 0, Mp_amu, 0, name="0_Plasmasphere") # Zero-energy plasmasphere channel
    sPe = aP.SpecParams(num_e, alamMin_e, alamMax_e, dtWolf, EFLAV, Me_amu, EFUDGE, name='Hot Electrons')  # Parameters to create electron channels
    sPp = aP.SpecParams(num_p, alamMin_p, alamMax_p, dtWolf, PFLAV, Mp_amu, PFUDGE, name='Hot Protons'  )  # Parameters to create proton channels

    # Stuff into an AlamData object, which will automatically generate the fully realized species distribution for us
    # NOTE: ORDER MATTERS. They will show up in the 1D lamc in this order
    alamData = aD.AlamData([sPe, sPp]) if noPsph else aD.AlamData([sPsphere, sPe, sPp])


    # Tag a brand new output file
    with h5.File(fOut, 'w') as f5:
        kh5.StampHash(fOut)
        kh5.StampBranch(fOut)
        f5.attrs['ScriptArgs'] = kj.dumps(vars(args),noIndent=True)
    
    # Save parameter info to file so we can recover later if needed
    fileIO.saveAlamParams(fOut, alamData)
    # Save main data to file
    fileIO.saveAlamData(fOut, alamData)

    # TODO:!! Check that we have the right number of alamis, flavs and fudges