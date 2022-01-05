#!/usr/bin/env python
#Generates RCM config data
import numpy as np
import argparse
from argparse import RawTextHelpFormatter
#import kaipy.rcm.lambdautils.genAlam as genAlam
from kaipy.rcm.lambdautils.AlamData import AlamParams

import kaipy.rcm.lambdautils2.AlamData as aD
import kaipy.rcm.lambdautils2.AlamParams as aP
import kaipy.rcm.lambdautils2.DistTypes as dT

import kaipy.rcm.lambdautils2.genAlam as genAlam
import kaipy.rcm.lambdautils2.fileIO as fileIO
import kaipy.rcm.lambdautils2.plotter as plotter


EFLAV = 1
PFLAV = 2

EFUDGE = 1./3.
PFUDGE = 0.0


def L_to_bVol(L):  # L shell [Re] to V [Re/nT]
    bsurf_nT = 3.11E4
    print(np.sqrt(1.0/L))
    colat = np.arcsin(np.sqrt(1.0/L))

    cSum = 35*np.cos(colat) - 7*np.cos(3*colat) +(7./5.)*np.cos(5*colat) - (1./7.)*np.cos(7*colat)
    cSum /= 64.
    s8 = np.sin(colat)**8
    V = 2*cSum/s8/bsurf_nT
    return V


if __name__ == "__main__":

#Arg parsing
	fOut = "rcmconfig.h5"
	apDef = AlamParams()  # Default rcm lambda parameters

	MainS = """Generates RCM configuration data
	"""
	parser = argparse.ArgumentParser(description=MainS, formatter_class=RawTextHelpFormatter)
	parser.add_argument('-o',type=str,default=fOut,metavar="fOut",help="Output file name (default: %(default)s)")
	parser.add_argument('-nop',action='store_true',default=False,help="Do not add zero loss first channel (default: %(default)s)")
	parser.add_argument('-ne', type=int,default=apDef.num_e, help="Number of electron channels (default: %(default)s)")
	parser.add_argument('-np', type=int,default=apDef.num_p, help="Number of proton channels (default: %(default)s)")
	parser.add_argument('-amine', type=float,default=apDef.aMin_e, help="Min. lambda for electrons (default: %(default)s)")
	parser.add_argument('-aminp', type=float,default=apDef.aMin_p, help="Min. lambda for protons (default: %(default)s)")
	parser.add_argument('-kt', type=float,default=apDef.ktMax/1E3, help="Highest thermal energy [keV] RCM should resolve at L_kt (default: %(default)s [keV])")
	parser.add_argument('-L', type=float,default=apDef.L_kt, help="L shell [R_e] at which kt should be resolved (default: %(default)s [R_e])")
	parser.add_argument('-tiote', type=float,default=apDef.tiote, help="Ratio between temperatures of ions and electrons (default: %(default)s)")
	parser.add_argument('-p1', type=float,default=apDef.p1, help="(default: %(default)s)")
	parser.add_argument('-p2', type=float,default=apDef.p2, help="(default: %(default)s)")


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

	#Determine proton channel limits based on resolving a certain (proton) temperature at given L
	bVol = L_to_bVol(L_kt)
	vm = bVol**(-2/3)
	alamMin_p = aminp
	alamMax_p = 10*(ktMax/vm)
	#Electron max based on proton channel max and ti/te
	alamMin_e = -1*np.abs(amine)  # Make sure its negative
	alamMax_e = -1*alamMax_p/tiote

	dtWolf = dT.DT_Wolf(p1=2,p2=2)  # Lambda channels will have a (slightly modified) Wolf distribution type

	n1 = 5; join1 = 100  # eV
	#valueSpecs = [dT.ValueSpec(alamMin_p, join1, n=n1, c=10, scaleType='log'),
	#	dT.ValueSpec(join1, alamMax_p, n=num_p-n1, c=10, scaleType='log')]
					
	valueSpecs = [dT.ValueSpec(alamMin_p, alamMax_p, n=num_p, scaleType='spacing_lin')]
	dtValSpec = dT.DT_ValueSpec(specList=valueSpecs)

	sPe = aP.SpecParams(num_p, alamMin_p, alamMax_p, dtWolf, EFLAV, EFUDGE, name='Protons_Wolf')  # Parameters to create electron channels
	sPp = aP.SpecParams(num_p, alamMin_p, alamMax_p, dtValSpec, PFLAV, PFUDGE, name='Protons_ValSpec_scaling-lin')  # Parameters to create proton channels
	alamParams = aP.AlamParams(True,[sPe, sPp])  # (doUsePsphere, List[SpecParams])
	alamData = genAlam.genAlamDataFromParams(alamParams)  # Use AlamParams to generate all of the lambda distributions

	plotter.plotLambdasBySpec(alamData.specs,yscale='log',vm=vm)

	print("Writing RCM configuration to %s"%(fOut))
	fileIO.saveRCMConfig(alamData,params=alamParams,fname=fOut)


