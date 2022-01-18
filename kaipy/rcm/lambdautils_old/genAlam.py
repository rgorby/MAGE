import numpy as np
import h5py as h5
import argparse
from argparse import RawTextHelpFormatter

from kaipy.rcm.lambdautils.dpetadp import L_to_bVol
import kaipy.rcm.lambdautils.alamTester as alamTester
import kaipy.rcm.lambdautils.plotter as plotter
from kaipy.rcm.lambdautils.AlamData import AlamParams

try:
	import kaipy.rcm.rcminit as rcminit
	hasKaiRCM = True
except:
	print("Couldn't load kaipy.rcm.init, can't add dktable to rcmconfig file.")
	hasKaiRCM = False

EFLAV = 1
PFLAV = 2

EFUDGE = 1./3.
PFUDGE = 0.0

distTypes = {'lin', 'log'}

# Takes full alams, flavs, fudges, and dktable arrays
# Writes arrays to file in rcmconfig.h5 format
def genh5(fname, inputParams, doTests=True, doShowPlot=False):

	alams, flavs, fudges = genAlams(inputParams, doTests=doTests, doShowPlot=doShowPlot)
	attrs = inputParams.getAttrs()

	f5 = h5.File(fname, 'w')
	f5.create_dataset('alamc', data=alams)
	f5.create_dataset('ikflavc', data=flavs)
	f5.create_dataset('fudgec', data=fudges)
	if hasKaiRCM:
		dktab = rcminit.LoadDKTab()
		f5.create_dataset('dktable', data=dktab)
	else:
		print("Couldn't load kaipy.rcm.init, can't add dktable to rcmconfig file.")

	for key in attrs.keys():
		f5.attrs[key] = attrs[key]
	f5.close()

# Adds completely new channel to alamc, ikflavc, and fudges
def addPsphere(alams, flavs, fudges):
	alam_wp = np.concatenate(([0], alams))
	flav_wp = np.concatenate(([1], flavs))
	fudge_wp = np.concatenate(([0], fudges))
	return alam_wp, flav_wp, fudge_wp

# Takes alams as an array, and single values for ikflav and fudge
# Returns arrays ready to write in rcmconfig.h5 format
def genspecData(alams, flav, fudge):
	#Very dumb
	flavs = np.array([flav for i in range(len(alams))])
	fudges = np.array([fudge for i in range(len(alams))])
	return alams, flavs, fudges

def wolfAlam(k, kmin, kmax, lammin, lammax, p1, p2):
	kfrac = (k-kmin)/(kmax-kmin)
	pstar = (1-kfrac)*p1 + kfrac*p2

	lammax2 = lammax - lammin
	#return lammax*((k - kmin + 0.5)/(kmax-kmin + 0.5))**p
	return lammax2*((k - kmin + 0.5)/(kmax-kmin + 0.5))**pstar + lammin

# Main function to get things going
def genAlams(alamParams, doTests=True, doShowPlot=False):

	#Unpack the basics
	bVol = L_to_bVol(alamParams.L_kt)
	vm = bVol**(-2/3)
	alamMin_p = alamParams.aMin_p
	alamMax_p = 10*(alamParams.ktMax/vm)
	numAlam_p = alamParams.num_p

	alamMin_e = alamParams.aMin_e
	alamMax_e = -1*alamMax_p/alamParams.tiote
	numAlam_e = alamParams.num_e

	# Gen desired alam distribution
	lbase = 10
	if alamParams.distType == 'lin':
		alams_e = np.linspace(alamMin_e, alamMax_e, numAlam_e, endpoint=True)
		alams_p = np.linspace(alamMin_p, alamMax_p, numAlam_p, endpoint=True)
	elif alamParams.distType == 'log':
		aminlog_e = np.log(np.abs(alamMin_e))/np.log(lbase)
		amaxlog_e = np.log(np.abs(alamMax_e))/np.log(lbase)
		aminlog_p = np.log(np.abs(alamMin_p))/np.log(lbase)
		amaxlog_p = np.log(np.abs(alamMax_p))/np.log(lbase)
		alams_e = -1*np.logspace(aminlog_e, amaxlog_e, numAlam_e, base=lbase, endpoint=True)
		alams_p = np.logspace(aminlog_p, amaxlog_p, numAlam_p, base=lbase, endpoint=True)
	elif alamParams.distType == 'wolf':
		alams_e = np.array([])
		alams_p = np.array([])
		p1 = alamParams.p1
		p2 = alamParams.p2
		for k in range(numAlam_e):
			a = wolfAlam(k, 0, numAlam_e, alamMin_e, alamMax_e, p1, p2)
			alams_e = np.append(alams_e, a)
		for k in range(numAlam_p):
			a = wolfAlam(k, 0, numAlam_p, alamMin_p, alamMax_p, p1, p2)
			alams_p = np.append(alams_p, a)

	# Questionable implementation
	alams_e, flavs_e, fudges_e = genspecData(alams_e, EFLAV, EFUDGE)
	alams_p, flavs_p, fudges_p = genspecData(alams_p, PFLAV, PFUDGE)

	# Combine all species
	alams  = np.concatenate((alams_e , alams_p ))
	flavs  = np.concatenate((flavs_e , flavs_p ))
	fudges = np.concatenate((fudges_e, fudges_p))

	# Add on plasmasphere channel if desired
	if alamParams.doAddPsphere:
		alams, flavs, fudges = addPsphere(alams, flavs, fudges)

	# Run lambda channel tests
	if doTests:
		alamTester.testAlam_arr(alams)
	
	# TODO: Write superplot, where fig is returned, and we will save it here
	plotter.plotChannels(alams, doShow=doShowPlot)

	return alams, flavs, fudges


"""
	TODO: Maybe need to make a config file in the future to reduce command line arguments
			Especially if we want multi-ion species
"""
if __name__ == "__main__":

	MainS = """Generate rcmconfig.h5 file with custom alam distribution"""

	ofname_in = "tmprcmconfig.h5"
	noPSphere = False
	doShowPlot = False

	dP = AlamParams()  # Default params
	
	parser = argparse.ArgumentParser(description=MainS, formatter_class=RawTextHelpFormatter)
	parser.add_argument('--outfile', type=str, metavar="<filename>",default=ofname_in,
			help="Name of output file (default: %(default)s)")
	parser.add_argument('--L', type=float, default=dP.L_kt,
			help="L shell [R_e] at which tmax should be resolved (default: %(default)s R_E)")
	parser.add_argument('--tmax', type=float,default=dP.ktMax/1E3,
        	help="Max desired kT to resolve at L [keV] (default: %(default)s keV)")
	parser.add_argument('--tiote', type=float,default=dP.tiote,
        	help="T_ion/T_electron (default: %(default)s)")
	parser.add_argument('--nrgmin_e', type=float,default=dP.aMin_e,
        	help="Mimumim electron energy [eV] (excluding plasmasphere) (default: %(default)s eV)")
	parser.add_argument('--nrgmin_p', type=float,default=dP.aMin_p,
        	help="Mimumim proton energy [eV] (default: %(default)s eV)")
	parser.add_argument('--num_e', type=float,default=dP.num_e,
        	help="Number of electron channels  (excluding plasmasphere) (default: %(default)s)")
	parser.add_argument('--num_p', type=float,default=dP.num_p,
        	help="Number of proton channels (default: %(default)s)")
	parser.add_argument('--p1', type=float,default=dP.p1,
        	help="Value for p1 parameter in Wolf distribution (default: %(default)s)")
	parser.add_argument('--p2', type=float,default=dP.p2,
        	help="Value for p2 parameter in Wolf distribution (default: %(default)s)")
	parser.add_argument("--noPsphere", default=False, action='store_true',
			help="Don't add plasmasphere channel (default: %(default)s)")
	parser.add_argument("--noTests", default=False, action='store_true',
			help="Don't run any tests (default: %(default)s)")
	parser.add_argument("--show", default=False, action='store_true',
			help="Show superplot after saving (default: %(default)s)")
	
	args = parser.parse_args()
	ofname = args.outfile
	doTests = not args.noTests
	doShowPlot = args.show

	inputParams = AlamParams(
		distType='wolf',
		num_e = args.num_e,
		num_p = args.num_p,
		alamMin_e = args.nrgmin_e,
		alamMin_p = args.nrgmin_p,
		ktMax = args.tmax*1E3,
		L_kt = args.L,
		tiote = args.tiote,
		p1 = args.p1,
		p2 = args.p2,
		addPsphere = not args.noPsphere  
		)

	alams, flavs, fudges = genAlams(inputParams, doTests=doTests, doShowPlot=doShowPlot)

	genh5(ofname, alams, flavs, fudges)
