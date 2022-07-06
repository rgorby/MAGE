#!/usr/bin/env python	
#Generates RCM config data
import numpy as np
import argparse
from argparse import RawTextHelpFormatter

import kaipy.rcm.lambdautils.AlamData as aD
import kaipy.rcm.lambdautils.AlamParams as aP
import kaipy.rcm.lambdautils.DistTypes as dT

import kaipy.rcm.lambdautils.genAlam as genAlam
from kaipy.rcm.wmutils.wmData import wmParams
import kaipy.rcm.wmutils.genWM as genWM
import kaipy.rcm.lambdautils.fileIO as fileIO
import kaipy.rcm.lambdautils.plotter as plotter


EFLAV = 1
PFLAV = 2

EFUDGE = 1./3.
PFUDGE = 0.0


def L_to_bVol(L):  # L shell [Re] to V [Re/nT]
    bsurf_nT = 3.11E4
    colat = np.arcsin(np.sqrt(1.0/L))

    cSum = 35*np.cos(colat) - 7*np.cos(3*colat) +(7./5.)*np.cos(5*colat) - (1./7.)*np.cos(7*colat)
    cSum /= 64.
    s8 = np.sin(colat)**8
    V = 2*cSum/s8/bsurf_nT
    return V


if __name__ == "__main__":

	#Arg parsing
	fOut = "rcmconfig.h5"
	num_e  = 39
	num_p  = 120
	aminp  = 10
	amine  = -1
	ktMax  = 15  # [keV]
	L_kt   = 10
	tiote  = 4
	wolfP1 = 3
	wolfP2 = 1
	plotChoices = ['none', 'spec', 'vs']


	MainS = """Generates RCM configuration data
	"""
	parser = argparse.ArgumentParser(description=MainS, formatter_class=RawTextHelpFormatter)
	parser.add_argument('-o',type=str,default=fOut,metavar="fOut",help="Output file name (default: %(default)s)")
	parser.add_argument('-nop',action='store_true',default=False,help="Do not add zero loss first channel (default: %(default)s)")
	parser.add_argument('-ne', type=int,default=num_e, help="Number of electron channels (default: %(default)s)")
	parser.add_argument('-np', type=int,default=num_p, help="Number of proton channels (default: %(default)s)")
	parser.add_argument('-amine', type=float,default=amine, help="Min. lambda for electrons (default: %(default)s)")
	parser.add_argument('-aminp', type=float,default=aminp, help="Min. lambda for protons (default: %(default)s)")
	parser.add_argument('-kt', type=float,default=ktMax, help="Highest thermal energy [keV] RCM should resolve at L_kt (default: %(default)s [keV])")
	parser.add_argument('-L', type=float,default=L_kt, help="L shell [R_e] at which kt should be resolved (default: %(default)s [R_e])")
	parser.add_argument('-tiote', type=float,default=tiote, help="Ratio between temperatures of ions and electrons (default: %(default)s)")
	parser.add_argument('-addWM', type=bool,default=False, help="Add wave models to existing rcmconfig file (default: %(default)s)")
	parser.add_argument('-i', type=str,default=fOut,metavar="fIn", help="Input file name when addWM is true (default: %(default)s)")
	parser.add_argument('-waveModel', type=bool,default=False, help="Use wave models in the electron/ion loss (default: %(default)s)")
	parser.add_argument('-p1', type=float,default=wolfP1, help="Wolf low-energy  p* (default: %(default)s)")
	parser.add_argument('-p2', type=float,default=wolfP2, help="Wolf high-energy p* (default: %(default)s)")
	parser.add_argument('-plotType', choices=plotChoices,default=plotChoices[0], help="Plot mode (default: %(default)s)")


	#Finalize parsing
	args = parser.parse_args()
	fOut = args.o
	num_e = args.ne
	num_p = args.np
	aminp = args.aminp
	amine = args.amine
	ktMax = args.kt*1E3  # [keV to eV]
	L_kt = float(args.L)
	tiote = args.tiote
	addWM = args.addWM
	waveModel = args.waveModel
	fIn = args.i	
	plotType = args.plotType

	if addWM:
		tauParams = wmParams(dim = 4, nKp = 7, nMLT = 25, nL = 41, nEk = 155)
		genWM.genh5(fIn,fOut,tauParams,useWMh5 = True)
	else:
		#Determine proton channel limits based on resolving a certain (proton) temperature at given L
		bVol = L_to_bVol(L_kt)
		vm = bVol**(-2/3)
		alamMin_p = aminp
		alamMax_p = 10*(ktMax/vm)
		#Electron max based on proton channel max and ti/te
		alamMin_e = -1*np.abs(amine)  # Make sure its negative
		alamMax_e = -1*alamMax_p/tiote

		dtWolf = dT.DT_Wolf(p1=3,p2=1)  # Lambda channels will have a (slightly modified) Wolf distribution type

		sPe = aP.SpecParams(num_e, alamMin_e, alamMax_e, dtWolf, EFLAV, EFUDGE, name='Electrons')  # Parameters to create electron channels
		sPp = aP.SpecParams(num_p, alamMin_p, alamMax_p, dtWolf, PFLAV, PFUDGE, name='Protons')  # Parameters to create proton channels
		alamParams = aP.AlamParams(True,[sPe, sPp])  # (doUsePsphere, List[SpecParams])
		alamParams.tiote = tiote
		alamParams.ktMax = ktMax
		alamParams.L_kt = L_kt
		alamData = genAlam.genAlamDataFromParams(alamParams)  # Use AlamParams to generate all of the lambda distributions

		if plotType == 'spec':
			plotter.plotLambdasBySpec(alamData.specs,yscale='log',L=L_kt)
		elif plotType == 'vs':
			plotter.plotLambdas_Val_Spac(alamData.specs,yscale='log',L=L_kt)

		fileIO.saveRCMConfig(alamData,params=alamParams,fname=fOut)
	
		if waveModel == True:
			tauParams = wmParams(dim = 4, nKp = 7, nMLT = 25, nL = 41, nEk = 155)	
			genWM.genh5(fOut,fOut,tauParams,useWMh5 = True)
	
	print("Wrote RCM configuration to %s"%(fOut))

