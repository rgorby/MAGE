#!/usr/bin/env python	
#Generates RCM config data
import numpy as np
import argparse
from argparse import RawTextHelpFormatter

import kaipy.kaiTools as kT

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


if __name__ == "__main__":

	#Arg parsing
	fOut = "rcmconfig.h5"
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


	MainS = """Generates RCM configuration data
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
	parser.add_argument('--nop',action='store_true',default=False,help="Do not add zero loss first channel (default: %(default)s)")
	parser.add_argument('--noWaveModel',action='store_true',default=False, help="Don't use wave models in the electron/ion loss (default: %(default)s)")
	parser.add_argument('--addWM', action='store_true',default=False, help="Add wave models to an existing rcmconfig file, input file needed to be presented (default: %(default)s)")
	parser.add_argument('-i', type=str,default=fOut,metavar="fIn", help="Input file name when addWM is true (default: %(default)s)")


	#Finalize parsing
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
	addWM = args.addWM
	noWaveModel = args.noWaveModel
	fIn = args.i	
	plotType = args.plotType

	if addWM:
		tauParams = wmParams(dim = 4, nKp = 7, nMLT = 25, nL = 41, nEk = 155)
		genWM.genh5(fIn,fOut,tauParams,useWM = True)
	else:
		# Determine proton channel limits based on resolving a certain (proton) temperature at given L
		bVol = kT.L_to_bVol(L_kt)
		vm = bVol**(-2/3)
		alamMin_p = eminp/vm
		alamMax_p = emaxp/vm
		alamMin_e = -1*emine/vm
		alamMax_e = -1*emaxe/vm

		dtWolf = dT.DT_Wolf(p1=wolfP1,p2=wolfP2)  # Lambda channels will have a (slightly modified) Wolf distribution type

		sPe = aP.SpecParams(num_e, alamMin_e, alamMax_e, dtWolf, EFLAV, EFUDGE, name='Electrons')  # Parameters to create electron channels
		sPp = aP.SpecParams(num_p, alamMin_p, alamMax_p, dtWolf, PFLAV, PFUDGE, name='Protons'  )  # Parameters to create proton channels
		alamParams = aP.AlamParams(True,[sPe, sPp])  # (doUsePsphere, List[SpecParams])
		alamParams.emine = emine
		alamParams.eminp = eminp
		alamParams.emaxe = emaxe
		alamParams.emaxp = emaxp
		alamParams.L_kt = L_kt
		alamData = genAlam.genAlamDataFromParams(alamParams)  # Use AlamParams to generate all of the lambda distributions


		# Save
		fileIO.saveRCMConfig(alamData,params=alamParams,fname=fOut)
		# Add data needed for wavemodel
		if not noWaveModel:
			tauParams = wmParams(dim = 4, nKp = 7, nMLT = 25, nL = 41, nEk = 155, dimTDS = 1, nEkTDS = 109)	
			genWM.genh5(fOut,fOut,tauParams,useWM = True)

		print("Wrote RCM configuration to %s"%(fOut))


		# Plotting

		if plotType == 'spec':  # 1 figure per species
			plotter.plotLambdasBySpec(alamData.specs,yscale='log',L=L_kt)
		elif plotType == 'vs':  # 2 figures (value and spacing), group all species
			plotter.plotLambdas_Val_Spac(alamData.specs,yscale='log',L=L_kt)

		if plotType != 'none':  # Show energy range covered (assuming dipole field)
			plotter.plotEnergyRange(alamData.specs, rInner=1.5, rOuter=15, rRes=100)

	


