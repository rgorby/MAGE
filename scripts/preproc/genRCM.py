#!/usr/bin/env python	
#Generates RCM config data

import argparse
from argparse import RawTextHelpFormatter
import kaipy.rcm.lambdautils.genAlam as genAlam
from kaipy.rcm.lambdautils.AlamData import AlamParams
from kaipy.rcm.wmutils.wmData import wmParams
import kaipy.rcm.wmutils.genWM as genWM

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

	print("Writing RCM configuration to %s"%(fOut))

	lamParams = AlamParams(num_e=args.ne, num_p=args.np,
						alamMin_e=args.amine, alamMin_p=args.aminp,
						ktMax=args.kt*1E3, L_kt=args.L, tiote=args.tiote,
						p1=args.p1, p2=args.p2,
						addPsphere = (not args.nop))
	genAlam.genh5(fOut, lamParams, doTests=False)
	tauParams = wmParams(dim = 4, nKp = 7, nMLT = 25, nL = 41, nEk = 155)	
	genWM.genh5(fOut,tauParams,useWMh5 = True)





